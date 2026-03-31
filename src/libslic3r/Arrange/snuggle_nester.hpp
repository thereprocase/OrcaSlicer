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

    // Rotation control
    bool   lock_rotation     = false;  // true = XY only, preserve user's Z rotation

    // Fitness weights
    float  w_compactness     = 1.0f;
    float  w_height_center   = 0.3f;

    // Politeness
    double timeout_seconds   = 30.0;   // Hard timeout
    size_t yield_every_gens  = 1;      // Yield CPU this often
};

// ── Per-part placement (3DOF: x, y, z_rotation) ──────────
struct Placement {
    float x    = 0.0f;  // mm, bed coordinates
    float y    = 0.0f;  // mm, bed coordinates
    float zrot = 0.0f;  // radians
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
    SnuggleNester(const NesterConfig &cfg, uint64_t seed = 42)
        : cfg_(cfg), rng_(seed) {}

    NesterResult run(
        const std::vector<PartInfo> &parts,
        NesterProgressFn progress = nullptr)
    {
        auto t_start = std::chrono::steady_clock::now();
        NesterResult result;
        size_t n_parts = parts.size();

        if (n_parts == 0) {
            result.time_ms = 0;
            return result;
        }

        // ── Initialize population ──────────────────────────
        std::vector<Individual> pop(cfg_.population_size);
        for (auto &ind : pop) {
            ind.placements.resize(n_parts);
            randomize_placement(ind, parts);
        }

        // Seed a few individuals with a greedy center placement
        if (pop.size() >= 4) {
            seed_greedy_center(pop[0], parts);
            seed_greedy_line(pop[1], parts);
        }

        // ── Evolution loop ─────────────────────────────────
        Individual best;
        best.collision_count = SIZE_MAX;

        for (size_t gen = 0; gen < cfg_.max_generations; gen++) {
            // Timeout check
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - t_start).count();
            if (elapsed > cfg_.timeout_seconds) {
                result.timed_out = true;
                break;
            }

            // Evaluate all individuals
            for (auto &ind : pop) {
                evaluate(ind, parts);
            }

            // Find best
            for (const auto &ind : pop) {
                if (ind.is_better_than(best)) {
                    best = ind;
                }
            }

            // Progress callback
            if (progress && !progress(gen, cfg_.max_generations, best)) {
                break; // User abort
            }

            // Early exit: found a good feasible solution and fitness is stable
            if (gen > 20 && best.is_feasible()) {
                // Check if we've improved in the last 20 generations
                // (tracked via simple heuristic: if fitness hasn't changed much)
                // For now, let it run — the timeout is the real exit
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
                mutate(child, parts);

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
        result.generations_run = cfg_.max_generations;
        result.time_ms        = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        return result;
    }

private:
    NesterConfig cfg_;
    std::mt19937 rng_;

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
            float margin = parts[i].grid.nx * parts[i].grid.voxel_size * 0.5f;
            p.x = randf(margin, cfg_.bed_width_mm - margin);
            p.y = randf(margin, cfg_.bed_height_mm - margin);
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

            // Clamp to bed
            float margin = 10.0f;
            p.x = std::clamp(p.x, margin, cfg_.bed_width_mm - margin);
            p.y = std::clamp(p.y, margin, cfg_.bed_height_mm - margin);

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

            // Clamp to bed
            p.x = std::clamp(p.x, 0.0f, cfg_.bed_width_mm);

            cursor_x += part_width + cfg_.min_gap_mm;
            // Wrap to next row if needed
            if (cursor_x > cfg_.bed_width_mm - 10.0f) {
                cursor_x = cfg_.min_gap_mm;
                cursor_y += 40.0f; // rough row height
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

    // ── Crossover: uniform per-part ───────────────────────
    void crossover(const Individual &a, const Individual &b,
                   Individual &child, size_t n_parts)
    {
        for (size_t i = 0; i < n_parts; i++) {
            if (randf(0, 1) < 0.5f) {
                child.placements[i] = a.placements[i];
            } else {
                child.placements[i] = b.placements[i];
            }
        }
    }

    // ── Mutation ──────────────────────────────────────────
    void mutate(Individual &ind, const std::vector<PartInfo> &parts) {
        for (size_t i = 0; i < ind.placements.size(); i++) {
            float roll = randf(0, 1);
            auto &p = ind.placements[i];

            if (roll < cfg_.mutation_nudge_prob) {
                // Nudge: small XY perturbation (+ rotation if unlocked)
                p.x += randf(-cfg_.nudge_xy_mm, cfg_.nudge_xy_mm);
                p.y += randf(-cfg_.nudge_xy_mm, cfg_.nudge_xy_mm);
                if (!cfg_.lock_rotation)
                    p.zrot += randf(-cfg_.nudge_rot_deg, cfg_.nudge_rot_deg) * (3.14159265f / 180.0f);
            } else if (roll < cfg_.mutation_nudge_prob + cfg_.mutation_shuffle_prob) {
                // Shuffle: swap XY with another random part (rotation stays with original part if locked)
                size_t j = randi(0, ind.placements.size() - 1);
                if (cfg_.lock_rotation) {
                    // Only swap XY, keep each part's locked rotation
                    std::swap(p.x, ind.placements[j].x);
                    std::swap(p.y, ind.placements[j].y);
                } else {
                    std::swap(ind.placements[i], ind.placements[j]);
                }
            } else {
                // Wildcard: completely random new XY (+ rotation if unlocked)
                float margin = parts[i].grid.nx * parts[i].grid.voxel_size * 0.5f;
                p.x = randf(margin, cfg_.bed_width_mm - margin);
                p.y = randf(margin, cfg_.bed_height_mm - margin);
                if (!cfg_.lock_rotation)
                    p.zrot = randf(0, 2.0f * 3.14159265f);
            }

            // Clamp to bed bounds (always)
            p.x = std::clamp(p.x, 0.0f, cfg_.bed_width_mm);
            p.y = std::clamp(p.y, 0.0f, cfg_.bed_height_mm);

            // Normalize rotation
            while (p.zrot < 0) p.zrot += 2.0f * 3.14159265f;
            while (p.zrot > 2.0f * 3.14159265f) p.zrot -= 2.0f * 3.14159265f;
        }
    }

    // ── Evaluate an individual ────────────────────────────
    //    1. Collision check (pairwise voxel AND)
    //    2. Bed bounds check
    //    3. Fitness scoring (only matters if feasible)
    void evaluate(Individual &ind, const std::vector<PartInfo> &parts) {
        size_t n = parts.size();
        ind.collision_count = 0;
        ind.oob_count = 0;
        ind.fitness = 0.0f;

        // ── Collision: pairwise voxel overlap ──────────────
        for (size_t i = 0; i < n; i++) {
            const auto &pi = ind.placements[i];
            Vec3f off_i = {pi.x, pi.y, 0.0f};

            for (size_t j = i + 1; j < n; j++) {
                const auto &pj = ind.placements[j];
                Vec3f off_j = {pj.x, pj.y, 0.0f};

                // Note: Z-rotation not yet applied to voxel grid
                // (we check at grid origin + XY offset for now;
                //  rotation support comes with rotated grid lookup)
                size_t c = VoxelGrid::collision_count(
                    parts[i].grid, off_i,
                    parts[j].grid, off_j);
                ind.collision_count += c;
            }

            // ── Bed bounds check ───────────────────────────
            // Check if any part of the voxel grid extends beyond bed
            const auto &g = parts[i].grid;
            float part_min_x = pi.x + g.origin.x;
            float part_min_y = pi.y + g.origin.y;
            float part_max_x = pi.x + g.origin.x + g.nx * g.voxel_size;
            float part_max_y = pi.y + g.origin.y + g.ny * g.voxel_size;

            if (part_min_x < 0) ind.oob_count += (size_t)(-part_min_x / g.voxel_size);
            if (part_min_y < 0) ind.oob_count += (size_t)(-part_min_y / g.voxel_size);
            if (part_max_x > cfg_.bed_width_mm) ind.oob_count += (size_t)((part_max_x - cfg_.bed_width_mm) / g.voxel_size);
            if (part_max_y > cfg_.bed_height_mm) ind.oob_count += (size_t)((part_max_y - cfg_.bed_height_mm) / g.voxel_size);
        }

        // ── Fitness (only meaningful if feasible) ──────────
        if (ind.is_feasible()) {
            // Compactness: minimize bounding rectangle of all part centers
            float min_x = 1e18f, max_x = -1e18f;
            float min_y = 1e18f, max_y = -1e18f;

            for (size_t i = 0; i < n; i++) {
                const auto &pi = ind.placements[i];
                const auto &g = parts[i].grid;
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
            // Normalize
            float max_possible_height = 0;
            for (const auto &p : parts) max_possible_height += p.max_height_mm;
            if (max_possible_height > 0)
                height_score /= max_possible_height;

            ind.fitness = cfg_.w_compactness * compactness
                        + cfg_.w_height_center * height_score;
        }
    }
};

} // namespace snuggle
