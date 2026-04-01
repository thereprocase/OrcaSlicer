// auto_snuggle.hpp — Automatic parameter selection for the radial nester
//
// Computes optimal voxel size, direction count, rotation count, and timeout
// from part count, part sizes, bed dimensions, and hardware capabilities.
// Targets a 2-10 second arrange with the best quality that fits the budget.

#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

namespace snuggle {

struct AutoSnuggleConfig {
    // Derived parameters
    float  voxel_mm       = 2.0f;
    int    n_directions    = 24;
    int    n_rotations     = 24;
    float  timeout_s       = 10.0f;
    float  step_mm         = 2.0f;  // = voxel_mm

    // Passed through from user (sacred — never overridden)
    float  min_gap_mm      = 1.5f;
    float  bed_margin_mm   = 3.0f;
    bool   lock_rotation   = false;
    int    rotation_step   = 15;

    // Diagnostic
    enum class Tier { FINE, NORMAL, COARSE, EMERGENCY };
    Tier   quality_tier    = Tier::NORMAL;
    float  estimated_time_s = 0.0f;
    float  cells_across    = 50.0f;  // median part size / voxel

    const char* tier_name() const {
        switch (quality_tier) {
            case Tier::FINE:      return "Fine";
            case Tier::NORMAL:    return "Normal";
            case Tier::COARSE:    return "Coarse";
            case Tier::EMERGENCY: return "Emergency";
        }
        return "Unknown";
    }
};

// Compute optimal parameters from plate state.
// part_max_dims: max(width, depth) of each part's XY bounding box in mm.
// Sacred settings are read from padding_mm, lock_rot, rot_step and passed through.
inline AutoSnuggleConfig compute_auto_config(
    size_t                    n_parts,
    const std::vector<float>& part_max_dims,
    float                     bed_w_mm,
    float                     bed_h_mm,
    float                     padding_mm,
    bool                      lock_rot,
    int                       rot_step,
    bool                      has_gpu)
{
    AutoSnuggleConfig cfg;

    // Sacred settings: copy verbatim
    cfg.min_gap_mm    = std::max(0.5f, padding_mm);
    cfg.bed_margin_mm = cfg.min_gap_mm;
    cfg.lock_rotation = lock_rot || rot_step <= 0;
    cfg.rotation_step = rot_step;

    if (n_parts == 0 || part_max_dims.empty()) return cfg;

    // Compute median and max part dimensions
    std::vector<float> dims = part_max_dims;
    std::sort(dims.begin(), dims.end());
    float median_dim = dims[dims.size() / 2];
    float max_dim    = dims.back();

    // Base voxel: target ~60 cells across the median part
    float base_voxel = median_dim / 60.0f;

    // N-scaling: above 10 parts, coarsen to stay in time budget
    float n_scale = 1.0f;
    if (n_parts > 10)
        n_scale = 1.0f + 0.5f * std::log2((float)n_parts / 10.0f);

    float voxel = base_voxel * n_scale;

    // Vase trap: cap grid at ~120 cells across the largest part
    float max_voxel_for_largest = max_dim / 120.0f;
    voxel = std::max(voxel, max_voxel_for_largest);

    // Clamp to sane range
    voxel = std::clamp(voxel, 0.5f, 5.0f);

    // GPU bonus: can afford finer resolution (only worth it above 5 parts
    // where dispatch overhead is amortized across enough candidates)
    if (has_gpu && n_parts > 5) {
        voxel *= 0.75f;
        voxel = std::clamp(voxel, 0.5f, 5.0f);
    }

    cfg.voxel_mm = voxel;
    cfg.step_mm  = voxel;

    // Direction count: scale inversely with N
    if (n_parts <= 5)       cfg.n_directions = 36;
    else if (n_parts <= 15) cfg.n_directions = 24;
    else if (n_parts <= 30) cfg.n_directions = 16;
    else                    cfg.n_directions = 12;

    // Rotation count: from user's step preference, scaled for N
    if (cfg.lock_rotation) {
        cfg.n_rotations = 1;
    } else {
        int step = std::max(1, cfg.rotation_step);
        int base_rots = 360 / step;
        if (n_parts > 20) base_rots = std::min(base_rots, 12);
        if (n_parts > 35) base_rots = std::min(base_rots, 8);
        cfg.n_rotations = std::max(1, base_rots);
    }

    // Timeout: scale with N, guard at 30s
    cfg.timeout_s = std::clamp(2.0f + 0.3f * (float)n_parts, 2.0f, 30.0f);

    // Quality tier from resolution
    cfg.cells_across = median_dim / voxel;
    if (cfg.cells_across >= 60.0f)      cfg.quality_tier = AutoSnuggleConfig::Tier::FINE;
    else if (cfg.cells_across >= 35.0f) cfg.quality_tier = AutoSnuggleConfig::Tier::NORMAL;
    else if (cfg.cells_across >= 20.0f) cfg.quality_tier = AutoSnuggleConfig::Tier::COARSE;
    else                                cfg.quality_tier = AutoSnuggleConfig::Tier::EMERGENCY;

    // Time estimate (rough cost model)
    float vox_est = (float)n_parts * 0.05f;
    float cache_est = (float)n_parts * (float)cfg.n_rotations * 0.002f;
    float place_est = 0.0f;
    for (size_t i = 0; i < n_parts; i++)
        place_est += (float)(i + 1) * (float)cfg.n_directions * (float)cfg.n_rotations * 0.00005f;
    cfg.estimated_time_s = vox_est + cache_est + place_est;

    return cfg;
}

} // namespace snuggle
