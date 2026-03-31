// snuggle_nester.hpp — Genetic Arranger for 3D Print Bed Nesting
//
// Takes voxelized parts and evolves arrangements using a genetic algorithm.
// All parts are bed-locked (Z=0, no tilt, no flip). Only XY + Z-rotation.
//
// POLITENESS GUARANTEES:
//   - Yields CPU every generation
//   - Configurable population size and generation count
//   - Hard timeout with graceful abort
//   - Memory bounded by population size * parts count * 3 floats
//   - Progress callback for UI integration / abort
//   - Single-threaded by default
//
// COLLISION: Hard constraint via feasibility-first selection (Deb 2000).
//   Infeasible individuals ALWAYS lose to feasible ones.
//   Among infeasible: rank by violation count (fewer = better).
//   Among feasible: rank by fitness.

#pragma once

#include "polite_voxelizer.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <functional>
#include <numeric>

namespace snuggle {

// Forward declaration — full type needed by callers that invoke
// evaluator methods. SnuggleArrange.cpp includes gpu_collision.hpp
// before this header to satisfy the dependency.
class CollisionEvaluator;

// ── Configuration ─────────────────────────────────────────
struct NesterConfig {
    // Population & evolution
    size_t population_size   = 512;    // Candidates per generation (keep modest for CPU)
    size_t max_generations   = 100;    // Max evolution cycles
    float  elitism_ratio     = 0.02f;  // Top 2% survive unchanged
    size_t tournament_size   = 4;      // Tournament selection bracket

    // Mutation rates
    float  mutation_nudge_prob    = 0.70f;  // Small XY + rotation tweak
    float  mutation_shuffle_prob  = 0.20f;  // Swap two parts' positions
    float  mutation_wild_prob     = 0.10f;  // Random new position for one part
    float  nudge_xy_mm           = 10.0f;   // Max nudge distance
    float  nudge_rot_deg         = 30.0f;   // Max rotation nudge

    // Bed definition
    float  bed_width_mm      = 256.0f;
    float  bed_height_mm     = 256.0f;
    float  min_gap_mm        = 5.0f;   // Minimum clearance between parts
    float  bed_margin_mm     = 4.0f;   // Safety margin from bed edge (voxel padding + gap)

    // Rotation control
    bool   lock_rotation     = false;  // true = XY only, preserve user's Z rotation
    float  rotation_step_rad = 0.0f;   // 0 = continuous, else snap to this increment
                                       // initial_zrot is always treated as the "home" position

    // Fitness weights
    float  w_compactness     = 1.0f;
    float  w_clustering      = 1.5f;   // penalize spread-out layouts
    float  w_height_center   = 0.3f;

    // Compaction
    bool   compact           = true;   // post-GA jiggle toward center
    double compact_timeout_s = 2.0;    // max seconds for compaction phase
    int    compact_max_sweeps = 20;    // max full sweeps

    // Politeness
    double timeout_seconds   = 30.0;   // Hard timeout
    size_t yield_every_gens  = 1;      // Yield CPU this often
};

// ── Per-part placement (3DOF: x, y, z_rotation) ──────────
struct Placement {
    float x    = 0.0f;  // mm, grid offset (internal to nester)
    float y    = 0.0f;  // mm, grid offset (internal to nester)
    float zrot = 0.0f;  // radians

    // Where the instance origin (mesh-local 0,0) lands on the bed,
    // in bed-relative mm. Computed after placement by the nester.
    float origin_bed_x = 0.0f;
    float origin_bed_y = 0.0f;
};

// ── Individual: one complete arrangement ──────────────────
struct Individual {
    std::vector<Placement> placements;   // One per part

    // Feasibility-first scoring (Deb 2000)
    size_t collision_count   = 0;   // Total overlapping voxels (0 = feasible)
    size_t oob_count         = 0;   // Voxels out of bed bounds
    float  fitness           = 0.0f; // Only meaningful when feasible

    bool is_feasible() const { return collision_count == 0 && oob_count == 0; }

    // Feasibility-first comparison: feasible always beats infeasible
    bool is_better_than(const Individual &other) const {
        if (is_feasible() && !other.is_feasible()) return true;
        if (!is_feasible() && other.is_feasible()) return false;
        if (!is_feasible() && !other.is_feasible()) {
            // Both infeasible: fewer violations wins
            return (collision_count + oob_count) < (other.collision_count + other.oob_count);
        }
        // Both feasible: higher fitness wins
        return fitness > other.fitness;
    }
};

// ── Part descriptor (input to nester) ─────────────────────
struct PartInfo {
    VoxelGrid   grid;           // Voxelized mesh
    float       max_height_mm;  // Tallest point (for height centering)
    float       hull_area_mm2;  // 2D footprint area (for compactness)
    float       initial_zrot;   // User's original Z rotation (radians), used when lock_rotation=true
    std::string name;           // For debug output

    PartInfo() : max_height_mm(0), hull_area_mm2(0), initial_zrot(0) {}
};

// ── Result ────────────────────────────────────────────────
struct NesterResult {
    std::vector<Placement> placements;
    float  fitness          = 0.0f;
    size_t collisions       = 0;
    size_t oob              = 0;
    size_t generations_run  = 0;
    double time_ms          = 0.0;
    bool   timed_out        = false;
    bool   feasible         = false;
};

// ── Progress callback ─────────────────────────────────────
using NesterProgressFn = std::function<bool(size_t gen, size_t max_gen,
                                            const Individual& best)>;

// ── The Nester ────────────────────────────────────────────
class SnuggleNester {
public:
    SnuggleNester(const NesterConfig &cfg, CollisionEvaluator* evaluator = nullptr, uint64_t seed = 42)
        : cfg_(cfg), evaluator_(evaluator), rng_(seed) {}

    NesterResult run(
        const std::vector<PartInfo> &parts,
        NesterProgressFn progress = nullptr)
    {
        auto t_start = std::chrono::steady_clock::now();
        NesterResult result;
        size_t n_parts = parts.size();

        if (n_parts == 0 || cfg_.population_size == 0) {
            result.time_ms = 0;
            return result;
        }

        // ── Build rotation cache (one-time cost) ──────────
        // Always build, even for lock_rotation — the evaluator needs it.
        // When locked, every bin gets the same grid (at the locked angle).
        build_rotation_cache(parts);
        if (evaluator_) {
            evaluator_->upload_grids(parts, rot_cache_);
        }

        // ── Initialize population ──────────────────────────
        std::vector<Individual> pop(cfg_.population_size);
        for (auto &ind : pop) {
            ind.placements.resize(n_parts);
            randomize_placement(ind, parts);
        }

        // Seed first 4 with heuristic strategies
        if (pop.size() >= 4) {
            seed_greedy_center(pop[0], parts);
            seed_greedy_line(pop[1], parts);
            seed_greedy_grid(pop[2], parts);
            seed_bottom_left(pop[3], parts);
        }

        // ── Evolution loop ─────────────────────────────────
        Individual best;
        best.collision_count = SIZE_MAX;

        // Stagnation tracking
        float best_fitness_seen = -1e18f;
        size_t gens_without_improvement = 0;
        float current_mutation_scale = 1.0f;
        size_t gen = 0;

        for (; gen < cfg_.max_generations; gen++) {
            // Timeout check
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - t_start).count();
            if (elapsed > cfg_.timeout_seconds) {
                result.timed_out = true;
                break;
            }

            // Evaluate all individuals
            if (evaluator_) {
                evaluator_->evaluate_batch(pop, parts,
                    cfg_.bed_width_mm, cfg_.bed_height_mm, cfg_.bed_margin_mm);
                for (auto &ind : pop) {
                    if (ind.is_feasible())
                        compute_fitness(ind, parts);
                }
            } else {
                for (auto &ind : pop) {
                    evaluate(ind, parts);
                }
            }

            // Find best
            for (const auto &ind : pop) {
                if (ind.is_better_than(best)) {
                    best = ind;
                }
            }

            // Track stagnation
            if (best.is_feasible() && best.fitness > best_fitness_seen + 0.001f) {
                best_fitness_seen = best.fitness;
                gens_without_improvement = 0;
                current_mutation_scale = 1.0f;
            } else {
                gens_without_improvement++;
            }

            // Adaptive mutation: crank up when stagnant
            if (gens_without_improvement >= 15) {
                current_mutation_scale = std::min(3.0f, current_mutation_scale + 0.3f);
            }

            // Early exit: converged
            if (best.is_feasible() && gens_without_improvement > 45) {
                result.generations_run = gen + 1;
                break;
            }

            // Progress callback
            if (progress && !progress(gen, cfg_.max_generations, best)) {
                break; // User abort
            }

            // ── Selection + Crossover + Mutation ───────────
            std::vector<Individual> next_pop;
            next_pop.reserve(cfg_.population_size);

            // Elitism: keep top N unchanged
            std::vector<size_t> ranking(pop.size());
            std::iota(ranking.begin(), ranking.end(), 0);
            std::sort(ranking.begin(), ranking.end(), [&](size_t a, size_t b) {
                return pop[a].is_better_than(pop[b]);
            });

            size_t n_elite = std::max((size_t)1,
                (size_t)(cfg_.elitism_ratio * cfg_.population_size));
            for (size_t i = 0; i < n_elite && i < ranking.size(); i++) {
                next_pop.push_back(pop[ranking[i]]);
            }

            // Fill rest with offspring
            while (next_pop.size() < cfg_.population_size) {
                // Tournament selection: pick 2 parents
                const Individual &parentA = tournament_select(pop);
                const Individual &parentB = tournament_select(pop);

                // Crossover
                Individual child;
                child.placements.resize(n_parts);
                crossover(parentA, parentB, child, n_parts);

                // Mutation
                mutate(child, parts, current_mutation_scale);

                next_pop.push_back(std::move(child));
            }

            pop = std::move(next_pop);

            // Polite yield
            if ((gen + 1) % cfg_.yield_every_gens == 0) {
                polite_yield();
            }
        }

        // Final evaluation of best
        evaluate(best, parts);

        auto t_end = std::chrono::steady_clock::now();

        result.placements     = best.placements;
        result.fitness        = best.fitness;
        result.collisions     = best.collision_count;
        result.oob            = best.oob_count;
        result.feasible       = best.is_feasible();
        result.generations_run = gen; // actual generations completed
        result.time_ms        = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        // Post-GA compaction: jiggle parts toward cluster center with
        // micro-rotations until gaps close to min_gap. Runs for up to
        // 2 seconds or 20 sweeps, whichever comes first.
        if (cfg_.compact) {
            compact_toward_center(result, parts);
            // Re-snap rotations after compaction (micro-rotations may violate step constraint)
            if (cfg_.rotation_step_rad > 0.001f) {
                for (size_t i = 0; i < result.placements.size() && i < parts.size(); i++)
                    result.placements[i].zrot = snap_rotation(result.placements[i].zrot, parts[i].initial_zrot);
            }
        }

        // Resolve where each part's instance origin ends up on the bed.
        compute_origin_positions(result, parts);

        return result;
    }

private:
    // ── Post-GA compaction ──────────────────────────────────
    // Iteratively jiggle each part toward the cluster centroid.
    // Each sweep: sort by distance from center (farthest first),
    // binary-search the max inward step, optionally micro-rotate.
    // Stops after time limit or convergence.
    void compact_toward_center(
        NesterResult &result,
        const std::vector<PartInfo> &parts)
    {
        size_t n = result.placements.size();
        if (n < 2) return;

        auto t_start = std::chrono::steady_clock::now();
        bool allow_rot = !cfg_.lock_rotation;
        float vs = parts[0].grid.voxel_size;

        for (int sweep = 0; sweep < cfg_.compact_max_sweeps; sweep++) {
            // Timeout check
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration<double>(now - t_start).count() > cfg_.compact_timeout_s)
                break;

            // Compute cluster centroid
            float cx = 0, cy = 0;
            for (size_t i = 0; i < n; i++) {
                cx += result.placements[i].x;
                cy += result.placements[i].y;
            }
            cx /= n; cy /= n;

            // Sort indices by distance from centroid (farthest first)
            std::vector<size_t> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                float da = (result.placements[a].x - cx) * (result.placements[a].x - cx)
                         + (result.placements[a].y - cy) * (result.placements[a].y - cy);
                float db = (result.placements[b].x - cx) * (result.placements[b].x - cx)
                         + (result.placements[b].y - cy) * (result.placements[b].y - cy);
                return da > db;
            });

            int total_moved = 0;

            for (size_t idx : order) {
                auto &pl = result.placements[idx];
                float dx = cx - pl.x;
                float dy = cy - pl.y;
                float dist = std::sqrt(dx*dx + dy*dy);
                if (dist < vs) continue;  // already at center

                float dir_x = dx / dist;
                float dir_y = dy / dist;

                // Build rotated grid for this part
                VoxelGrid rot_i = parts[idx].grid.rotated_copy(pl.zrot);

                // Pre-compute rotated copies for collision partners (avoid redundant copies)
                std::vector<const VoxelGrid*> rot_others(n, nullptr);
                std::vector<VoxelGrid> rot_others_storage;
                rot_others_storage.reserve(n);
                for (size_t j = 0; j < n; j++) {
                    if (j == idx) continue;
                    rot_others_storage.push_back(parts[j].grid.rotated_copy(result.placements[j].zrot));
                    rot_others[j] = &rot_others_storage.back();
                }

                // Binary search: max step toward center without collision
                float lo = 0, hi = dist;
                float best_step = 0;
                float best_rot_delta = 0;

                for (int bs = 0; bs < 12; bs++) {
                    float mid = (lo + hi) * 0.5f;
                    float test_x = pl.x + dir_x * mid;
                    float test_y = pl.y + dir_y * mid;

                    // Check collision with all other parts at test position
                    bool collides = false;
                    Vec3f off_test = {test_x, test_y, 0.0f};
                    for (size_t j = 0; j < n; j++) {
                        if (j == idx || !rot_others[j]) continue;
                        Vec3f off_j = {result.placements[j].x, result.placements[j].y, 0.0f};
                        if (VoxelGrid::collision_count(rot_i, off_test, *rot_others[j], off_j) > 0) {
                            collides = true;
                            break;
                        }
                    }

                    // Also check bed bounds
                    float margin = cfg_.bed_margin_mm;
                    float pmin_x = test_x + rot_i.origin.x;
                    float pmin_y = test_y + rot_i.origin.y;
                    float pmax_x = pmin_x + rot_i.nx * rot_i.voxel_size;
                    float pmax_y = pmin_y + rot_i.ny * rot_i.voxel_size;
                    if (pmin_x < margin || pmin_y < margin ||
                        pmax_x > cfg_.bed_width_mm - margin ||
                        pmax_y > cfg_.bed_height_mm - margin)
                        collides = true;

                    if (!collides) {
                        best_step = mid;
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }

                // Try micro-rotations at the best step position
                if (allow_rot && best_step > vs) {
                    float base_x = pl.x + dir_x * best_step;
                    float base_y = pl.y + dir_y * best_step;
                    float test_angles[] = {-0.05f, 0.05f, -0.1f, 0.1f};

                    for (float da : test_angles) {
                        float test_rot = pl.zrot + da;
                        VoxelGrid rot_test = parts[idx].grid.rotated_copy(test_rot);
                        Vec3f off_test = {base_x, base_y, 0.0f};

                        bool ok = true;
                        for (size_t j = 0; j < n && ok; j++) {
                            if (j == idx) continue;
                            VoxelGrid rot_j = parts[j].grid.rotated_copy(result.placements[j].zrot);
                            Vec3f off_j = {result.placements[j].x, result.placements[j].y, 0.0f};
                            if (VoxelGrid::collision_count(rot_test, off_test, rot_j, off_j) > 0)
                                ok = false;
                        }

                        if (ok) {
                            // Try stepping even further with this rotation
                            float extra_lo = best_step, extra_hi = dist;
                            float extra_best = best_step;
                            for (int bs2 = 0; bs2 < 8; bs2++) {
                                float mid2 = (extra_lo + extra_hi) * 0.5f;
                                float tx = pl.x + dir_x * mid2;
                                float ty = pl.y + dir_y * mid2;
                                VoxelGrid rot2 = parts[idx].grid.rotated_copy(test_rot);
                                Vec3f off2 = {tx, ty, 0.0f};
                                bool col2 = false;
                                for (size_t j = 0; j < n && !col2; j++) {
                                    if (j == idx) continue;
                                    VoxelGrid rj2 = parts[j].grid.rotated_copy(result.placements[j].zrot);
                                    Vec3f oj2 = {result.placements[j].x, result.placements[j].y, 0.0f};
                                    if (VoxelGrid::collision_count(rot2, off2, rj2, oj2) > 0)
                                        col2 = true;
                                }
                                if (!col2) { extra_best = mid2; extra_lo = mid2; }
                                else { extra_hi = mid2; }
                            }
                            if (extra_best > best_step) {
                                best_step = extra_best;
                                best_rot_delta = da;
                            }
                            break;  // take first improving rotation
                        }
                    }
                }

                // Apply the move
                if (best_step > vs) {
                    pl.x += dir_x * best_step;
                    pl.y += dir_y * best_step;
                    pl.zrot += best_rot_delta;
                    total_moved++;
                }
            }

            if (total_moved == 0) break;  // converged
            polite_yield();
        }
    }

    // Compute where each part's instance origin (mesh-local 0,0) lands
    // on the bed, given the grid placement and rotation.
    void compute_origin_positions(
        NesterResult &result,
        const std::vector<PartInfo> &parts)
    {
        for (size_t i = 0; i < result.placements.size() && i < parts.size(); i++) {
            auto &pl = result.placements[i];
            float ox = parts[i].grid.origin.x;
            float oy = parts[i].grid.origin.y;
            float ca = std::cos(pl.zrot);
            float sa = std::sin(pl.zrot);

            // Instance origin (0,0 in mesh space) is at (-ox, -oy) relative
            // to the rotation pivot (grid min corner). Rotate, then add
            // the pivot position and placement offset.
            pl.origin_bed_x = pl.x + ox * (1.0f - ca) + oy * sa;
            pl.origin_bed_y = pl.y + oy * (1.0f - ca) - ox * sa;
        }
    }

    NesterConfig cfg_;
    CollisionEvaluator* evaluator_ = nullptr;  // non-owning, caller manages lifetime
    std::mt19937 rng_;

    // ── Rotation cache ────────────────────────────────────
    // Quantize angles to 1-degree bins. Pre-build on first access.
    // Avoids rebuilding rotated grids every evaluation (~100x speedup).
    static constexpr int ROT_CACHE_BINS = 360;
    std::vector<std::vector<VoxelGrid>> rot_cache_; // [part][angle_bin]
    bool rot_cache_built_ = false;

    void build_rotation_cache(const std::vector<PartInfo> &parts) {
        if (rot_cache_built_) return;
        size_t n = parts.size();
        rot_cache_.resize(n);

        if (cfg_.lock_rotation) {
            // Build one rotated copy per part and replicate to all bins.
            // This avoids 360x the compute while ensuring any bin lookup
            // returns a valid (non-empty) grid — critical for GPU evaluator.
            for (size_t i = 0; i < n; i++) {
                VoxelGrid rotated = parts[i].grid.rotated_copy(parts[i].initial_zrot);
                rot_cache_[i].resize(ROT_CACHE_BINS, rotated); // fill all bins with same copy
                polite_yield();
            }
        } else {
            // Build all 360 bins
            for (size_t i = 0; i < n; i++) {
                rot_cache_[i].resize(ROT_CACHE_BINS);
                for (int bin = 0; bin < ROT_CACHE_BINS; bin++) {
                    float angle = (float)bin * (2.0f * 3.14159265f / ROT_CACHE_BINS);
                    rot_cache_[i][bin] = parts[i].grid.rotated_copy(angle);
                }
                polite_yield();
            }
        }
        rot_cache_built_ = true;
    }

    // Snap rotation to step increment relative to part's initial_zrot.
    // If step is 0 or lock_rotation, returns initial_zrot unchanged.
    // Otherwise quantizes (angle - initial) to nearest step, adds back initial.
    float snap_rotation(float angle_rad, float initial_zrot) const {
        if (cfg_.lock_rotation) return initial_zrot;
        if (cfg_.rotation_step_rad <= 0.001f) return angle_rad; // continuous
        float delta = angle_rad - initial_zrot;
        // Normalize delta to [0, 2pi)
        constexpr float TWO_PI = 2.0f * 3.14159265f;
        while (delta < 0) delta += TWO_PI;
        while (delta >= TWO_PI) delta -= TWO_PI;
        // Snap to nearest step
        float steps = std::round(delta / cfg_.rotation_step_rad);
        return initial_zrot + steps * cfg_.rotation_step_rad;
    }

    const VoxelGrid& cached_rotated(size_t part_idx, float angle_rad) const {
        int bin = (int)std::floor(angle_rad * ROT_CACHE_BINS / (2.0f * 3.14159265f));
        bin = ((bin % ROT_CACHE_BINS) + ROT_CACHE_BINS) % ROT_CACHE_BINS;
        return rot_cache_[part_idx][bin];
    }

    // ── Random float in range ─────────────────────────────
    float randf(float lo, float hi) {
        return std::uniform_real_distribution<float>(lo, hi)(rng_);
    }
    size_t randi(size_t lo, size_t hi) {
        return std::uniform_int_distribution<size_t>(lo, hi)(rng_);
    }

    // ── Randomize an individual's placements ──────────────
    void randomize_placement(Individual &ind, const std::vector<PartInfo> &parts) {
        for (size_t i = 0; i < ind.placements.size(); i++) {
            auto &p = ind.placements[i];
            float part_margin = parts[i].grid.nx * parts[i].grid.voxel_size * 0.5f;
            float lo = std::max(part_margin, cfg_.bed_margin_mm);
            // Ensure lo < hi to avoid UB in uniform_real_distribution
            // Right bound must also account for part margin
            float hi_x = std::max(lo + 0.1f, cfg_.bed_width_mm - std::max(part_margin, cfg_.bed_margin_mm));
            float hi_y = std::max(lo + 0.1f, cfg_.bed_height_mm - std::max(part_margin, cfg_.bed_margin_mm));
            p.x = randf(lo, hi_x);
            p.y = randf(lo, hi_y);
            p.zrot = cfg_.lock_rotation
                ? parts[i].initial_zrot
                : randf(0.0f, 2.0f * 3.14159265f);
        }
    }

    // ── Greedy seed: center placement ─────────────────────
    void seed_greedy_center(Individual &ind, const std::vector<PartInfo> &parts) {
        float cx = cfg_.bed_width_mm * 0.5f;
        float cy = cfg_.bed_height_mm * 0.5f;

        // Sort parts by footprint area (largest first)
        std::vector<size_t> order(parts.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return parts[a].hull_area_mm2 > parts[b].hull_area_mm2;
        });

        // Place largest at center, spiral outward
        float radius = 0.0f;
        float angle = 0.0f;
        for (size_t idx : order) {
            auto &p = ind.placements[idx];
            p.x = cx + radius * std::cos(angle);
            p.y = cy + radius * std::sin(angle);
            p.zrot = cfg_.lock_rotation ? parts[idx].initial_zrot : 0.0f;

            // Clamp to bed (respect bed margin)
            p.x = std::clamp(p.x, cfg_.bed_margin_mm, cfg_.bed_width_mm - cfg_.bed_margin_mm);
            p.y = std::clamp(p.y, cfg_.bed_margin_mm, cfg_.bed_height_mm - cfg_.bed_margin_mm);

            float part_size = std::sqrt(parts[idx].hull_area_mm2) * 0.5f;
            radius += part_size + cfg_.min_gap_mm;
            angle += 2.4f; // Golden angle-ish for spiral
        }
    }

    // ── Greedy seed: line placement ───────────────────────
    void seed_greedy_line(Individual &ind, const std::vector<PartInfo> &parts) {
        float cursor_x = cfg_.min_gap_mm;
        float cursor_y = cfg_.bed_height_mm * 0.5f;

        for (size_t i = 0; i < parts.size(); i++) {
            float part_width = parts[i].grid.nx * parts[i].grid.voxel_size;
            auto &p = ind.placements[i];
            p.x = cursor_x + part_width * 0.5f;
            p.y = cursor_y;
            p.zrot = cfg_.lock_rotation ? parts[i].initial_zrot : 0.0f;

            // Clamp to bed (respect bed margin)
            p.x = std::clamp(p.x, cfg_.bed_margin_mm, cfg_.bed_width_mm - cfg_.bed_margin_mm);

            cursor_x += part_width + cfg_.min_gap_mm;
            // Wrap to next row if needed
            if (cursor_x > cfg_.bed_width_mm - 10.0f) {
                cursor_x = cfg_.min_gap_mm;
                cursor_y += 40.0f; // rough row height
            }
        }
    }

    // ── Greedy seed: grid placement ──────────────────────
    void seed_greedy_grid(Individual &ind, const std::vector<PartInfo> &parts) {
        if (parts.empty()) return;
        std::vector<size_t> order(parts.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return parts[a].hull_area_mm2 > parts[b].hull_area_mm2;
        });
        float avg_size = 0;
        for (const auto &p : parts) avg_size += std::sqrt(p.hull_area_mm2);
        avg_size = avg_size / parts.size() + cfg_.min_gap_mm;
        int cols = std::max(1, (int)(cfg_.bed_width_mm / avg_size));
        float cell_w = cfg_.bed_width_mm / cols;
        float cell_h = avg_size;
        for (size_t idx = 0; idx < order.size(); idx++) {
            size_t pi = order[idx];
            auto &p = ind.placements[pi];
            p.x = ((int)(idx % cols) + 0.5f) * cell_w;
            p.y = ((int)(idx / cols) + 0.5f) * cell_h + cfg_.min_gap_mm;
            p.zrot = cfg_.lock_rotation ? parts[pi].initial_zrot : 0.0f;
            p.x = std::clamp(p.x, cfg_.bed_margin_mm, cfg_.bed_width_mm - cfg_.bed_margin_mm);
            p.y = std::clamp(p.y, cfg_.bed_margin_mm, cfg_.bed_height_mm - cfg_.bed_margin_mm);
        }
    }

    // ── Greedy seed: bottom-left fill ────────────────────
    void seed_bottom_left(Individual &ind, const std::vector<PartInfo> &parts) {
        std::vector<size_t> order(parts.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return parts[a].hull_area_mm2 > parts[b].hull_area_mm2;
        });
        struct PlacedRect { float x, y, w, h; };
        std::vector<PlacedRect> placed;
        for (size_t idx : order) {
            float pw = parts[idx].grid.nx * parts[idx].grid.voxel_size;
            float ph = parts[idx].grid.ny * parts[idx].grid.voxel_size;
            auto &p = ind.placements[idx];
            p.zrot = cfg_.lock_rotation ? parts[idx].initial_zrot : 0.0f;
            bool found = false;
            for (float try_y = cfg_.bed_margin_mm; try_y < cfg_.bed_height_mm - ph; try_y += cfg_.min_gap_mm) {
                for (float try_x = cfg_.bed_margin_mm; try_x < cfg_.bed_width_mm - pw; try_x += cfg_.min_gap_mm) {
                    bool collides = false;
                    for (const auto &pr : placed) {
                        if (try_x < pr.x+pr.w+cfg_.min_gap_mm && try_x+pw > pr.x-cfg_.min_gap_mm &&
                            try_y < pr.y+pr.h+cfg_.min_gap_mm && try_y+ph > pr.y-cfg_.min_gap_mm) {
                            collides = true; try_x = pr.x+pr.w+cfg_.min_gap_mm-cfg_.min_gap_mm; break;
                        }
                    }
                    if (!collides) { p.x=try_x+pw*0.5f; p.y=try_y+ph*0.5f; placed.push_back({try_x,try_y,pw,ph}); found=true; break; }
                }
                if (found) break;
            }
            if (!found) {
                float hx=std::max(pw+0.1f,cfg_.bed_width_mm-pw), hy=std::max(ph+0.1f,cfg_.bed_height_mm-ph);
                p.x=randf(std::min(pw,hx),hx); p.y=randf(std::min(ph,hy),hy);
            }
        }
    }

    // ── Tournament selection ──────────────────────────────
    const Individual& tournament_select(const std::vector<Individual> &pop) {
        size_t best_idx = randi(0, pop.size() - 1);
        for (size_t t = 1; t < cfg_.tournament_size; t++) {
            size_t idx = randi(0, pop.size() - 1);
            if (pop[idx].is_better_than(pop[best_idx])) {
                best_idx = idx;
            }
        }
        return pop[best_idx];
    }

    // ── Crossover: inherit from fitter, mix in some from other ─
    void crossover(const Individual &a, const Individual &b,
                   Individual &child, size_t n_parts)
    {
        const Individual &fitter = a.is_better_than(b) ? a : b;
        const Individual &other  = a.is_better_than(b) ? b : a;
        child.placements = fitter.placements;
        for (size_t i = 0; i < n_parts; i++) {
            if (randf(0, 1) < 0.3f)
                child.placements[i] = other.placements[i];
        }
    }

    // ── Mutation (with adaptive scaling) ───────────────────
    void mutate(Individual &ind, const std::vector<PartInfo> &parts,
                float scale = 1.0f)
    {
        float nudge_range = cfg_.nudge_xy_mm * scale;
        float rot_range = cfg_.nudge_rot_deg * scale;

        for (size_t i = 0; i < ind.placements.size(); i++) {
            float roll = randf(0, 1);
            auto &p = ind.placements[i];

            if (roll < 0.55f) {
                // Nudge: small XY perturbation scaled by stagnation pressure
                p.x += randf(-nudge_range, nudge_range);
                p.y += randf(-nudge_range, nudge_range);
                if (!cfg_.lock_rotation)
                    p.zrot += randf(-rot_range, rot_range) * (3.14159265f / 180.0f);
            } else if (roll < 0.75f) {
                // Shuffle: swap positions with another part
                size_t j = randi(0, ind.placements.size() - 1);
                if (cfg_.lock_rotation) {
                    std::swap(p.x, ind.placements[j].x);
                    std::swap(p.y, ind.placements[j].y);
                } else {
                    std::swap(ind.placements[i], ind.placements[j]);
                }
            } else if (roll < 0.88f) {
                // Push-apart: find nearest neighbor and move away
                float best_dist = 1e18f;
                size_t nearest = i;
                for (size_t j = 0; j < ind.placements.size(); j++) {
                    if (j == i) continue;
                    float dx = ind.placements[j].x - p.x;
                    float dy = ind.placements[j].y - p.y;
                    float d = dx*dx + dy*dy;
                    if (d < best_dist) { best_dist = d; nearest = j; }
                }
                if (nearest != i) {
                    float dx = p.x - ind.placements[nearest].x;
                    float dy = p.y - ind.placements[nearest].y;
                    float len = std::sqrt(dx*dx + dy*dy);
                    if (len > 0.01f) {
                        float push = cfg_.min_gap_mm * scale;
                        p.x += dx/len * push;
                        p.y += dy/len * push;
                    }
                }
            } else {
                // Wildcard: completely random new position
                float wild_margin = std::max(parts[i].grid.nx * parts[i].grid.voxel_size * 0.5f,
                                             cfg_.bed_margin_mm);
                float wx_hi = std::max(wild_margin + 0.1f, cfg_.bed_width_mm - wild_margin);
                float wy_hi = std::max(wild_margin + 0.1f, cfg_.bed_height_mm - wild_margin);
                p.x = randf(wild_margin, wx_hi);
                p.y = randf(wild_margin, wy_hi);
                if (!cfg_.lock_rotation)
                    p.zrot = randf(0, 2.0f * 3.14159265f);
            }

            // Clamp to bed (respect bed margin)
            p.x = std::clamp(p.x, cfg_.bed_margin_mm, cfg_.bed_width_mm - cfg_.bed_margin_mm);
            p.y = std::clamp(p.y, cfg_.bed_margin_mm, cfg_.bed_height_mm - cfg_.bed_margin_mm);

            // Snap rotation to step increment relative to initial position
            p.zrot = snap_rotation(p.zrot, parts[i].initial_zrot);
        }
    }

    // ── Fitness scoring (CPU only, no collision data needed) ─
    //    Computes compactness + height-centering score.
    //    Called separately when using an external evaluator,
    //    or as part of evaluate() on the CPU fallback path.
    void compute_fitness(Individual &ind, const std::vector<PartInfo> &parts) {
        size_t n = parts.size();
        ind.fitness = 0.0f;

        if (!ind.is_feasible()) return;

        // Need rotated grid bounds for bounding-box calculation
        std::vector<const VoxelGrid*> rotated(n);
        std::vector<VoxelGrid> rotated_locked;
        if (cfg_.lock_rotation) {
            rotated_locked.resize(n);
            for (size_t i = 0; i < n; i++) {
                rotated_locked[i] = parts[i].grid.rotated_copy(ind.placements[i].zrot);
                rotated[i] = &rotated_locked[i];
            }
        } else {
            for (size_t i = 0; i < n; i++) {
                rotated[i] = &cached_rotated(i, ind.placements[i].zrot);
            }
        }

        // Compactness: minimize bounding rectangle of all parts
        float min_x = 1e18f, max_x = -1e18f;
        float min_y = 1e18f, max_y = -1e18f;

        for (size_t i = 0; i < n; i++) {
            const auto &pi = ind.placements[i];
            const auto &g = *rotated[i];
            float px_min = pi.x + g.origin.x;
            float py_min = pi.y + g.origin.y;
            float px_max = pi.x + g.origin.x + g.nx * g.voxel_size;
            float py_max = pi.y + g.origin.y + g.ny * g.voxel_size;

            min_x = std::min(min_x, px_min);
            min_y = std::min(min_y, py_min);
            max_x = std::max(max_x, px_max);
            max_y = std::max(max_y, py_max);
        }

        float bbox_area = (max_x - min_x) * (max_y - min_y);
        float bed_area = cfg_.bed_width_mm * cfg_.bed_height_mm;
        float compactness = 1.0f - (bbox_area / bed_area);
        compactness = std::clamp(compactness, 0.0f, 1.0f);

        // Clustering: penalize large average distance between part centers.
        // Prevents the GA from spreading parts in a line across the bed.
        float clustering = 0.0f;
        if (n >= 2) {
            // Compute centroid of all part centers
            float cx_sum = 0, cy_sum = 0;
            for (size_t i = 0; i < n; i++) {
                const auto &g = *rotated[i];
                cx_sum += ind.placements[i].x + g.origin.x + g.nx * g.voxel_size * 0.5f;
                cy_sum += ind.placements[i].y + g.origin.y + g.ny * g.voxel_size * 0.5f;
            }
            float centroid_x = cx_sum / n;
            float centroid_y = cy_sum / n;

            // Average distance from centroid (normalized by bed diagonal)
            float bed_diag = std::sqrt(cfg_.bed_width_mm * cfg_.bed_width_mm +
                                       cfg_.bed_height_mm * cfg_.bed_height_mm);
            float avg_dist = 0;
            for (size_t i = 0; i < n; i++) {
                const auto &g = *rotated[i];
                float pcx = ind.placements[i].x + g.origin.x + g.nx * g.voxel_size * 0.5f;
                float pcy = ind.placements[i].y + g.origin.y + g.ny * g.voxel_size * 0.5f;
                float dx = pcx - centroid_x;
                float dy = pcy - centroid_y;
                avg_dist += std::sqrt(dx*dx + dy*dy);
            }
            avg_dist /= (n * bed_diag);
            clustering = 1.0f - std::clamp(avg_dist, 0.0f, 1.0f);
        }

        // Height centering: tall parts near center score better
        float height_score = 0.0f;
        float bed_cx = cfg_.bed_width_mm * 0.5f;
        float bed_cy = cfg_.bed_height_mm * 0.5f;
        float max_dist = std::sqrt(bed_cx * bed_cx + bed_cy * bed_cy);

        for (size_t i = 0; i < n; i++) {
            float dx = ind.placements[i].x - bed_cx;
            float dy = ind.placements[i].y - bed_cy;
            float dist = std::sqrt(dx*dx + dy*dy);
            float centrality = 1.0f - (dist / max_dist);
            height_score += parts[i].max_height_mm * centrality;
        }
        float max_possible_height = 0;
        for (const auto &p : parts) max_possible_height += p.max_height_mm;
        if (max_possible_height > 0)
            height_score /= max_possible_height;

        // Aspect ratio: reward square-ish bounding boxes over long lines
        float bbox_w = max_x - min_x;
        float bbox_h = max_y - min_y;
        float aspect = (bbox_w > 0.01f && bbox_h > 0.01f)
            ? std::min(bbox_w, bbox_h) / std::max(bbox_w, bbox_h) : 0.0f;

        // Proximity: reward parts that are close but not colliding
        float proximity_score = 0.0f;
        float n_pair_count = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                float dx = ind.placements[i].x - ind.placements[j].x;
                float dy = ind.placements[i].y - ind.placements[j].y;
                float dist = std::sqrt(dx*dx + dy*dy);
                float ri = std::sqrt(parts[i].hull_area_mm2) * 0.5f;
                float rj = std::sqrt(parts[j].hull_area_mm2) * 0.5f;
                float ideal = ri + rj + cfg_.min_gap_mm;
                if (dist > 0.01f) proximity_score += std::min(ideal / dist, 1.0f);
                n_pair_count++;
            }
        }
        if (n_pair_count > 0) proximity_score /= n_pair_count;

        ind.fitness = cfg_.w_compactness * compactness
                    + cfg_.w_clustering * clustering
                    + cfg_.w_compactness * 0.4f * aspect
                    + cfg_.w_compactness * 0.3f * proximity_score
                    + cfg_.w_height_center * height_score;
    }

    // ── Evaluate an individual (CPU fallback) ─────────────
    //    Full inline path: collision + bounds + fitness.
    //    Used when no CollisionEvaluator is provided.
    void evaluate(Individual &ind, const std::vector<PartInfo> &parts) {
        size_t n = parts.size();
        ind.collision_count = 0;
        ind.oob_count = 0;
        ind.fitness = 0.0f;

        // ── Get rotated grids for each part ───────────────
        // When rotation is locked, use the original grid (no rotation applied).
        // When unlocked, use the pre-built rotation cache (quantized to 1° bins).
        // This avoids O(N * grid_volume) rotated_copy calls per evaluation.
        std::vector<const VoxelGrid*> rotated(n);
        std::vector<VoxelGrid> rotated_locked; // storage for lock_rotation case
        if (cfg_.lock_rotation) {
            rotated_locked.resize(n);
            for (size_t i = 0; i < n; i++) {
                rotated_locked[i] = parts[i].grid.rotated_copy(ind.placements[i].zrot);
                rotated[i] = &rotated_locked[i];
            }
        } else {
            for (size_t i = 0; i < n; i++) {
                rotated[i] = &cached_rotated(i, ind.placements[i].zrot);
            }
        }

        // ── Collision: pairwise voxel overlap ──────────────
        for (size_t i = 0; i < n; i++) {
            const auto &pi = ind.placements[i];
            Vec3f off_i = {pi.x, pi.y, 0.0f};

            for (size_t j = i + 1; j < n; j++) {
                const auto &pj = ind.placements[j];
                Vec3f off_j = {pj.x, pj.y, 0.0f};

                size_t c = VoxelGrid::collision_count(
                    *rotated[i], off_i,
                    *rotated[j], off_j);
                ind.collision_count += c;
            }

            // ── Bed bounds check ───────────────────────────
            const auto &g = *rotated[i];
            float part_min_x = pi.x + g.origin.x;
            float part_min_y = pi.y + g.origin.y;
            float part_max_x = pi.x + g.origin.x + g.nx * g.voxel_size;
            float part_max_y = pi.y + g.origin.y + g.ny * g.voxel_size;

            float margin = cfg_.bed_margin_mm;
            if (part_min_x < margin) ind.oob_count += (size_t)((margin - part_min_x) / g.voxel_size);
            if (part_min_y < margin) ind.oob_count += (size_t)((margin - part_min_y) / g.voxel_size);
            if (part_max_x > cfg_.bed_width_mm - margin) ind.oob_count += (size_t)((part_max_x - cfg_.bed_width_mm + margin) / g.voxel_size);
            if (part_max_y > cfg_.bed_height_mm - margin) ind.oob_count += (size_t)((part_max_y - cfg_.bed_height_mm + margin) / g.voxel_size);
        }

        // ── Fitness (only meaningful if feasible) ──────────
        if (ind.is_feasible()) {
            compute_fitness(ind, parts);
        }
    }
};

} // namespace snuggle
