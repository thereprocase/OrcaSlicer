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
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <chrono>

namespace snuggle {

// Forward from snuggle_nester.hpp — we reuse PartInfo and Placement
struct PartInfo;

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

    int   n_directions  = 36;    // radial angles to try (360/36 = 10° each)
    int   n_rotations   = 24;    // rotation angles to try (360/24 = 15° each)
    float step_mm       = 2.0f;  // slide distance per step (= voxel size)

    bool  lock_rotation = false;
    double timeout_s    = 30.0;
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
    auto get_rot = [&](size_t part_idx, float angle) -> const VoxelGrid& {
        static const VoxelGrid empty;
        if (part_idx >= rot_cache.size() || rot_cache[part_idx].empty())
            return empty;
        int bin = (int)std::floor(angle * ROT_CACHE_BINS / TWO_PI_F);
        bin = ((bin % ROT_CACHE_BINS) + ROT_CACHE_BINS) % ROT_CACHE_BINS;
        if (bin >= (int)rot_cache[part_idx].size()) bin = 0;
        return rot_cache[part_idx][bin];
    };

    // Helper: check if a placement is valid (no collision, within bed)
    auto is_valid = [&](size_t part_idx, float px, float py, float rot,
                        const std::vector<size_t>& placed_indices,
                        const std::vector<RadialPlacement>& placements) -> bool {
        const VoxelGrid& grid = get_rot(part_idx, rot);
        if (grid.nx == 0) return false;

        // Bed bounds check
        float pmin_x = px + grid.origin.x;
        float pmin_y = py + grid.origin.y;
        float pmax_x = pmin_x + grid.nx * grid.voxel_size;
        float pmax_y = pmin_y + grid.ny * grid.voxel_size;

        if (pmin_x < margin || pmin_y < margin ||
            pmax_x > cfg.bed_width_mm - margin ||
            pmax_y > cfg.bed_height_mm - margin)
            return false;

        // Collision check against all placed parts
        Vec3f off_i = {px, py, 0.0f};
        for (size_t j : placed_indices) {
            const auto& pj = placements[j];
            const VoxelGrid& grid_j = get_rot(j, pj.zrot);
            Vec3f off_j = {pj.x, pj.y, 0.0f};
            if (VoxelGrid::collision_count(grid, off_i, grid_j, off_j) > 0)
                return false;
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
    std::vector<std::pair<float, float>> directions; // (dx, dy) unit vectors
    for (int d = 0; d < cfg.n_directions; d++) {
        float angle = (float)d * TWO_PI_F / cfg.n_directions;
        directions.push_back({std::cos(angle), std::sin(angle)});
    }

    std::vector<size_t> placed_indices;

    for (size_t rank = 0; rank < n; rank++) {
        size_t idx = order[rank];

        // Timeout check
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration<double>(now - t_start).count() > cfg.timeout_s) {
            result.timed_out = true;
            break;
        }

        if (progress && !progress((int)placed_indices.size(), (int)n, parts[idx].name.c_str()))
            break;

        // Per-part rotation list
        std::vector<float> part_rots;
        if (cfg.lock_rotation) {
            part_rots = {parts[idx].initial_zrot};
        } else {
            part_rots = rot_angles;
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
        // First valid distance wins — guarantees tightest packing.
        struct Candidate {
            float x, y, rot, dist;
        };
        Candidate best = {0, 0, 0, max_slide + 1};
        bool found = false;

        for (float dist = 0; dist <= max_slide && !found; dist += cfg.step_mm) {
            for (const auto& [dx, dy] : directions) {
                float px = bed_cx + dx * dist;
                float py = bed_cy + dy * dist;

                for (float rot : part_rots) {
                    if (is_valid(idx, px, py, rot, placed_indices, result.placements)) {
                        best = {px, py, rot, dist};
                        found = true;
                        goto ring_done; // first valid at this distance is good enough
                    }
                }
            }
        }
        ring_done:;

        if (found) {
            result.placements[idx] = {best.x, best.y, best.rot, 0, 0, true};
            placed_indices.push_back(idx);
        }
        // else: part doesn't fit, stays unplaced
    }

    // Compute origin_bed positions (same formula as GA nester)
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
