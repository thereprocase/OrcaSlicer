// snuggle_radial.hpp — Deterministic radial expansion nester
//
// Places the tallest part at bed center, then slides each subsequent
// part outward along radial directions until it fits. Tries all
// rotation angles at each distance step. Picks the tightest fit.
//
// No genetic algorithm. No population. No stochastic search.
// Deterministic: same input always produces same output.
// Fast: O(N * directions * max_steps * rotations) collision checks.

#pragma once

#include "polite_voxelizer.hpp"
#include "snuggle_constants.hpp"
#include "snuggle_nester.hpp"
#include "gpu_collision.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <chrono>
#include <boost/log/trivial.hpp>

namespace snuggle {

struct RadialPlacement {
    float x = 0, y = 0;
    float zrot = 0;
    float origin_bed_x = 0, origin_bed_y = 0;
    bool placed = false;
};

struct RadialConfig {
    float bed_width_mm  = 256.0f;
    float bed_height_mm = 256.0f;
    float min_gap_mm    = 1.0f;
    float bed_margin_mm = 3.0f;

    int   n_directions  = 24;    // radial angles to try (360/24 = 15° each)
    int   n_rotations   = 24;    // rotation angles to try (360/24 = 15° each)
    float step_mm       = 2.0f;  // slide distance per step (= voxel size)

    bool  lock_rotation = false;
    double timeout_s    = 15.0;
};

struct RadialResult {
    std::vector<RadialPlacement> placements;
    int    placed_count = 0;
    int    total_count  = 0;
    double time_ms      = 0;
    bool   timed_out    = false;
};

using RadialProgressFn = std::function<bool(int placed, int total, const char* status)>;

inline RadialResult radial_arrange(
    const std::vector<PartInfo>& parts,
    const std::vector<std::vector<VoxelGrid>>& rot_cache,
    const RadialConfig& cfg,
    CollisionEvaluator* evaluator = nullptr,
    RadialProgressFn progress = nullptr)
{
    auto t_start = std::chrono::steady_clock::now();
    RadialResult result;
    size_t n = parts.size();
    result.total_count = (int)n;
    result.placements.resize(n);

    if (n == 0) return result;

    float bed_cx = cfg.bed_width_mm * 0.5f;
    float bed_cy = cfg.bed_height_mm * 0.5f;
    float margin = cfg.bed_margin_mm;
    float max_slide = std::sqrt(bed_cx * bed_cx + bed_cy * bed_cy); // bed diagonal

    // Sort parts tallest-first (by max_height_mm, then by footprint area)
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (std::abs(parts[a].max_height_mm - parts[b].max_height_mm) > 1.0f)
            return parts[a].max_height_mm > parts[b].max_height_mm;
        return parts[a].hull_area_mm2 > parts[b].hull_area_mm2;
    });

    // Helper: get rotated grid from cache
    // Lookup rotated grid from cache. Cache may have any number of bins
    // (24 for radial, 1 for lock_rotation). Index by angle → bin.
    auto get_rot = [&](size_t part_idx, float angle) -> const VoxelGrid& {
        static const VoxelGrid empty;
        if (part_idx >= rot_cache.size() || rot_cache[part_idx].empty())
            return empty;
        int n_bins = (int)rot_cache[part_idx].size();
        int bin = (int)std::floor(angle * n_bins / TWO_PI_F);
        bin = ((bin % n_bins) + n_bins) % n_bins;
        return rot_cache[part_idx][bin];
    };

    // Helper: check if a placement is valid (no collision, within bed)
    auto is_valid = [&](size_t part_idx, float px, float py, float rot,
                        const std::vector<size_t>& placed_indices,
                        const std::vector<RadialPlacement>& placements) -> bool {
        const VoxelGrid& grid = get_rot(part_idx, rot);
        if (grid.nx == 0) return false;

        // Bed bounds check (includes gap as margin from bed edges)
        float gap = cfg.min_gap_mm;
        float pmin_x = px + grid.origin.x;
        float pmin_y = py + grid.origin.y;
        float pmax_x = pmin_x + grid.nx * grid.voxel_size;
        float pmax_y = pmin_y + grid.ny * grid.voxel_size;

        if (pmin_x < margin || pmin_y < margin ||
            pmax_x > cfg.bed_width_mm - margin ||
            pmax_y > cfg.bed_height_mm - margin)
            return false;

        // Collision check with gap enforcement.
        // When gap >= half voxel size, use 9-probe pattern (center + 8 offsets).
        // When gap < half voxel size, sub-voxel shifts don't change the result —
        // use center probe only (Legolas optimization: avoids 8x wasted work).
        bool use_probes = gap >= grid.voxel_size * 0.5f;
        Vec3f off_i = {px, py, 0.0f};
        for (size_t j : placed_indices) {
            const auto& pj = placements[j];
            const VoxelGrid& grid_j = get_rot(j, pj.zrot);
            Vec3f off_j = {pj.x, pj.y, 0.0f};

            if (use_probes) {
                for (float dx : {0.0f, gap, -gap}) {
                    for (float dy : {0.0f, gap, -gap}) {
                        Vec3f off_shifted = {px + dx, py + dy, 0.0f};
                        if (VoxelGrid::collision_count(grid, off_shifted, grid_j, off_j) > 0)
                            return false;
                    }
                }
            } else {
                if (VoxelGrid::collision_count(grid, off_i, grid_j, off_j) > 0)
                    return false;
            }
        }
        return true;
    };

    // Build rotation angle list
    std::vector<float> rot_angles;
    if (cfg.lock_rotation) {
        // Only use each part's initial rotation
        // (handled per-part below)
    } else {
        for (int r = 0; r < cfg.n_rotations; r++)
            rot_angles.push_back((float)r * TWO_PI_F / cfg.n_rotations);
    }

    // Build radial direction list
    int n_dirs = std::max(1, cfg.n_directions);
    std::vector<std::pair<float, float>> directions; // (dx, dy) unit vectors
    for (int d = 0; d < n_dirs; d++) {
        float angle = (float)d * TWO_PI_F / cfg.n_directions;
        directions.push_back({std::cos(angle), std::sin(angle)});
    }

    std::vector<size_t> placed_indices;

    // Active search parameters — can be reduced mid-run if over budget
    int active_directions = n_dirs;
    int active_rotations  = (int)rot_angles.size();
    bool degraded = false;

    for (size_t rank = 0; rank < n; rank++) {
        size_t idx = order[rank];

        // Timeout check
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - t_start).count();
        if (elapsed > cfg.timeout_s) {
            result.timed_out = true;
            break;
        }

        // Graceful degradation: if running over budget, reduce search breadth
        // rather than timing out with unplaced parts.
        float remaining_frac = (float)(n - placed_indices.size()) / (float)n;
        float budget_used = (float)(elapsed / cfg.timeout_s);

        if (budget_used > 0.8f && remaining_frac > 0.4f && active_directions > 6) {
            active_directions = std::max(6, active_directions / 2);
            if (!degraded) {
                BOOST_LOG_TRIVIAL(warning) << "Snuggle: degrading search — "
                    << active_directions << " dirs (budget " << (int)(budget_used * 100) << "%)";
                degraded = true;
            }
        }
        if (budget_used > 0.9f && remaining_frac > 0.2f && active_rotations > 4) {
            active_rotations = std::max(4, active_rotations / 2);
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: degrading rotations — "
                << active_rotations << " rots (budget " << (int)(budget_used * 100) << "%)";
        }

        if (progress && !progress((int)placed_indices.size(), (int)n, parts[idx].name.c_str()))
            break;

        // Per-part rotation list (may be truncated by degradation)
        std::vector<float> part_rots;
        if (cfg.lock_rotation) {
            part_rots = {parts[idx].initial_zrot};
        } else {
            int n_rots = std::min(active_rotations, (int)rot_angles.size());
            part_rots.assign(rot_angles.begin(), rot_angles.begin() + n_rots);
        }

        // First part: try to place at center with each rotation
        if (placed_indices.empty()) {
            for (float rot : part_rots) {
                if (is_valid(idx, bed_cx, bed_cy, rot, placed_indices, result.placements)) {
                    result.placements[idx] = {bed_cx, bed_cy, rot, 0, 0, true};
                    placed_indices.push_back(idx);
                    break;
                }
            }
            if (!result.placements[idx].placed) {
                // Can't even place at center — try origin
                result.placements[idx] = {bed_cx, bed_cy, part_rots[0], 0, 0, false};
            }
            continue;
        }

        // Expand outward from center in concentric rings.
        // At each distance, try ALL directions and ALL rotations.
        // Collect all valid placements at the first distance that works,
        // then pick the tightest fit (minimizes distance to nearest
        // placed part — parts cluster instead of scattering).
        //
        // When an evaluator is provided, batch all candidates per ring
        // for GPU dispatch. Otherwise, use inline collision checks (CPU).
        struct Candidate {
            float x, y, rot, dist;
        };

        std::vector<Candidate> ring_candidates;

        // Build placed part info for evaluator batching
        std::vector<size_t> eval_placed_parts;
        std::vector<RadialCandidate> eval_placed_pos;
        if (evaluator) {
            for (size_t pi : placed_indices) {
                eval_placed_parts.push_back(pi);
                eval_placed_pos.push_back({
                    result.placements[pi].x,
                    result.placements[pi].y,
                    result.placements[pi].zrot
                });
            }
        }

        for (float dist = 0; dist <= max_slide; dist += cfg.step_mm) {
            ring_candidates.clear();

            int dirs_to_try = std::min(active_directions, (int)directions.size());

            if (evaluator) {
                // Batch path: collect candidates at this distance, evaluate together
                std::vector<RadialCandidate> batch;
                for (int d = 0; d < dirs_to_try; d++) {
                    float px = bed_cx + directions[d].first * dist;
                    float py = bed_cy + directions[d].second * dist;
                    for (float rot : part_rots) {
                        batch.push_back({px, py, rot});
                    }
                }

                std::vector<RadialCollisionResult> batch_results;
                evaluator->evaluate_radial(eval_placed_parts, eval_placed_pos,
                                           idx, batch, cfg.min_gap_mm,
                                           cfg.bed_width_mm, cfg.bed_height_mm,
                                           cfg.bed_margin_mm, batch_results);

                // Scan results: for each direction, take the first valid rotation
                size_t bi = 0;
                for (int d = 0; d < dirs_to_try; d++) {
                    for (size_t ri = 0; ri < part_rots.size(); ri++, bi++) {
                        if (bi < batch_results.size() && !batch_results[bi].collides) {
                            ring_candidates.push_back({batch[bi].x, batch[bi].y,
                                                       batch[bi].zrot, dist});
                            bi += part_rots.size() - ri - 1;
                            break;
                        }
                    }
                }
            } else {
                // Inline CPU path (no evaluator)
                for (int d = 0; d < dirs_to_try; d++) {
                    float px = bed_cx + directions[d].first * dist;
                    float py = bed_cy + directions[d].second * dist;

                    for (float rot : part_rots) {
                        if (is_valid(idx, px, py, rot, placed_indices, result.placements)) {
                            ring_candidates.push_back({px, py, rot, dist});
                            break;
                        }
                    }
                }
            }

            if (!ring_candidates.empty()) break;
        }

        if (!ring_candidates.empty()) {
            // Pick the candidate closest to existing parts (tight packing).
            // For each candidate, compute minimum distance to any placed part.
            // Pick the candidate with the smallest such distance.
            size_t best_idx = 0;
            float best_min_dist = max_slide + 1.0f;

            for (size_t ci = 0; ci < ring_candidates.size(); ci++) {
                float min_d = max_slide;
                for (size_t pi : placed_indices) {
                    float ddx = ring_candidates[ci].x - result.placements[pi].x;
                    float ddy = ring_candidates[ci].y - result.placements[pi].y;
                    float d = std::sqrt(ddx * ddx + ddy * ddy);
                    min_d = std::min(min_d, d);
                }
                if (min_d < best_min_dist) {
                    best_min_dist = min_d;
                    best_idx = ci;
                }
            }

            auto& c = ring_candidates[best_idx];
            result.placements[idx] = {c.x, c.y, c.rot, 0, 0, true};
            placed_indices.push_back(idx);
        }
        // else: part doesn't fit, stays unplaced
    }

    // Compute origin_bed: where the instance origin lands on the bed.
    // Compensates for rotation pivot offset (Gandalf Option B).
    for (size_t i = 0; i < n; i++) {
        auto& pl = result.placements[i];
        if (!pl.placed) continue;

        float ox = parts[i].grid.origin.x;
        float oy = parts[i].grid.origin.y;
        float ca = std::cos(pl.zrot);
        float sa = std::sin(pl.zrot);

        pl.origin_bed_x = pl.x + ox * (1.0f - ca) + oy * sa;
        pl.origin_bed_y = pl.y + oy * (1.0f - ca) - ox * sa;
    }

    result.placed_count = (int)placed_indices.size();
    auto t_end = std::chrono::steady_clock::now();
    result.time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    return result;
}

} // namespace snuggle
