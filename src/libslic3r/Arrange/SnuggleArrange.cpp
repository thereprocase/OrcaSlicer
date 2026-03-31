// SnuggleArrange.cpp — Orchestrator for Snuggle 3D nesting in OrcaSlicer

#include "SnuggleArrange.hpp"
#include "SnuggleTransform.hpp"
#include "SnuggleVoxelize.hpp"
#include "snuggle_nester.hpp"
#include "libslic3r/libslic3r.h"    // SCALING_FACTOR
#include <boost/log/trivial.hpp>
#include <chrono>
#include <cassert>

namespace Slic3r { namespace arrangement {

// Maximum number of parts we'll attempt to nest with Snuggle.
// Beyond this, the genetic algorithm's O(N^2) collision checks become slow.
static constexpr size_t MAX_SNUGGLE_PARTS = 128;

// Default voxel size in mm. Trades resolution for speed.
static constexpr float DEFAULT_VOXEL_SIZE_MM = 1.0f;

void snuggle_arrange(
    ArrangePolygons& items,
    const ArrangePolygons& excludes,
    const Points& bed,
    const ArrangeParams& params,
    const Model& model)
{
    auto t_start = std::chrono::steady_clock::now();

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] === Starting Snuggle 3D nesting ===";
    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] items=" << items.size()
        << " excludes=" << excludes.size()
        << " bed_points=" << bed.size();

    // ── Edge case: nothing to arrange ────────────────────────────────
    if (items.empty()) {
        BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] No items to arrange, returning.";
        return;
    }

    // ── Edge case: too many parts ────────────────────────────────────
    if (items.size() > MAX_SNUGGLE_PARTS) {
        BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Too many parts ("
            << items.size() << " > " << MAX_SNUGGLE_PARTS
            << "). Snuggle nesting skipped — caller should fall back to default arranger.";
        // Leave all items UNARRANGED so the caller knows to fall back
        return;
    }

    // ── Step 1: Convert bed to mm dimensions ─────────────────────────
    float bed_width_mm = 0, bed_height_mm = 0;
    float bed_origin_x_mm = 0, bed_origin_y_mm = 0;
    snuggle_xform::bed_dimensions_mm(bed, bed_width_mm, bed_height_mm,
                                     bed_origin_x_mm, bed_origin_y_mm);

    if (bed_width_mm <= 0 || bed_height_mm <= 0) {
        BOOST_LOG_TRIVIAL(error) << "[SnuggleArrange] Invalid bed dimensions: "
            << bed_width_mm << " x " << bed_height_mm << " mm. Aborting.";
        return;
    }

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Bed: "
        << bed_width_mm << " x " << bed_height_mm << " mm"
        << " origin=(" << bed_origin_x_mm << "," << bed_origin_y_mm << ") mm";

    // ── Step 2: Voxelize parts ───────────────────────────────────────
    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Voxelizing " << items.size() << " parts...";

    auto t_vox_start = std::chrono::steady_clock::now();

    std::vector<SnugglePartData> part_data = voxelize_for_snuggle(
        model, items, DEFAULT_VOXEL_SIZE_MM);

    auto t_vox_end = std::chrono::steady_clock::now();
    double vox_ms = std::chrono::duration<double, std::milli>(t_vox_end - t_vox_start).count();

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Voxelization took " << vox_ms << " ms";

    if (part_data.size() != items.size()) {
        BOOST_LOG_TRIVIAL(error) << "[SnuggleArrange] Voxelization returned "
            << part_data.size() << " parts but expected " << items.size()
            << ". Aborting.";
        return;
    }

    // ── Build PartInfo vector for nester ─────────────────────────────
    std::vector<snuggle::PartInfo> nester_parts;
    nester_parts.reserve(part_data.size());
    for (auto& pd : part_data) {
        nester_parts.push_back(std::move(pd.info));
    }

    // ── Step 3: Configure and run nester ─────────────────────────────
    snuggle::NesterConfig cfg;
    cfg.bed_width_mm  = bed_width_mm;
    cfg.bed_height_mm = bed_height_mm;

    // Convert min_obj_distance from scaled coords to mm via SnuggleTransform
    cfg.min_gap_mm = snuggle_xform::scaled_to_mm(params.min_obj_distance);
    if (cfg.min_gap_mm < 1.0f)
        cfg.min_gap_mm = 1.0f; // Minimum 1mm gap for safety

    // Rotation control: default to LOCKED (safe — XY only, no rotation).
    // Snuggle's voxel collision check doesn't rotate grids, so rotated
    // placements are not collision-verified. Lock rotation unless the user
    // explicitly enables it AND snuggle_lock_rotation is false.
    cfg.lock_rotation = true; // Safe default: XY placement only
    if (params.allow_rotations && !params.snuggle_lock_rotation) {
        cfg.lock_rotation = false; // User explicitly wants rotation
        BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Rotation UNLOCKED — "
            "rotated placements are NOT collision-verified in this version.";
    }

    // Time budget: use 80% of any configured timeout, or 30s default
    cfg.timeout_seconds = 30.0;

    // Scale population/generations by part count for better results
    if (items.size() <= 4) {
        cfg.population_size = 256;
        cfg.max_generations = 80;
    } else if (items.size() <= 16) {
        cfg.population_size = 512;
        cfg.max_generations = 100;
    } else {
        cfg.population_size = 256;  // Fewer candidates, more parts = slower eval
        cfg.max_generations = 60;
    }

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Nester config:"
        << " bed=" << cfg.bed_width_mm << "x" << cfg.bed_height_mm << " mm"
        << " gap=" << cfg.min_gap_mm << " mm"
        << " lock_rot=" << cfg.lock_rotation
        << " pop=" << cfg.population_size
        << " gens=" << cfg.max_generations
        << " timeout=" << cfg.timeout_seconds << "s";

    // Log excludes (not yet integrated as obstacles — future work)
    if (!excludes.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] " << excludes.size()
            << " exclude regions present but NOT yet used as Snuggle obstacles."
            << " Parts may overlap with excluded regions.";
    }

    auto t_nest_start = std::chrono::steady_clock::now();

    snuggle::SnuggleNester nester(cfg);
    snuggle::NesterResult nest_result = nester.run(nester_parts,
        [](size_t gen, size_t max_gen, const snuggle::Individual& best) -> bool {
            if (gen % 10 == 0 || gen == max_gen - 1) {
                BOOST_LOG_TRIVIAL(debug) << "[SnuggleArrange] Gen " << gen << "/" << max_gen
                    << " collisions=" << best.collision_count
                    << " oob=" << best.oob_count
                    << " fitness=" << best.fitness
                    << (best.is_feasible() ? " [FEASIBLE]" : " [infeasible]");
            }
            return true; // continue
        });

    auto t_nest_end = std::chrono::steady_clock::now();
    double nest_ms = std::chrono::duration<double, std::milli>(t_nest_end - t_nest_start).count();

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Nester finished in " << nest_ms << " ms"
        << " (" << nest_result.generations_run << " generations)"
        << " feasible=" << nest_result.feasible
        << " collisions=" << nest_result.collisions
        << " oob=" << nest_result.oob
        << " fitness=" << nest_result.fitness
        << (nest_result.timed_out ? " [TIMED OUT]" : "");

    // ── Step 4: Convert results back to ArrangePolygon fields ────────
    if (nest_result.placements.size() != items.size()) {
        BOOST_LOG_TRIVIAL(error) << "[SnuggleArrange] Nester returned "
            << nest_result.placements.size() << " placements but expected "
            << items.size() << ". Leaving items unarranged.";
        return;
    }

    if (nest_result.feasible) {
        // All parts fit on one bed
        for (size_t i = 0; i < items.size(); ++i) {
            const snuggle::Placement& pl = nest_result.placements[i];

            Vec2crd out_translation;
            double out_rotation;
            snuggle_xform::placement_to_arrange(
                pl.x, pl.y, pl.zrot,
                bed_origin_x_mm, bed_origin_y_mm,
                out_translation, out_rotation);

            items[i].translation = out_translation;
            items[i].rotation    = out_rotation;
            items[i].bed_idx     = 0; // Placed on the physical bed

            BOOST_LOG_TRIVIAL(debug) << "[SnuggleArrange] Part " << i
                << " '" << part_data[i].name << "'"
                << " placed at snuggle=(" << pl.x << "," << pl.y << ") mm"
                << " rot=" << pl.zrot << " rad"
                << " -> scaled=(" << out_translation.x() << "," << out_translation.y() << ")"
                << " bed_idx=0";
        }
    } else {
        // Not all parts fit. Try to identify which ones are OOB and mark them.
        BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Nester did NOT find a feasible solution."
            << " collisions=" << nest_result.collisions
            << " oob=" << nest_result.oob;

        // Apply placements anyway — parts that are within bed bounds get bed_idx=0,
        // parts that are out of bounds get UNARRANGED.
        for (size_t i = 0; i < items.size(); ++i) {
            const snuggle::Placement& pl = nest_result.placements[i];

            // Check if this part's placement is within bed bounds.
            // Use rotation-aware bounding box: compute the rotated AABB diagonal
            // so parts at 45 degrees are correctly bounded.
            const snuggle::VoxelGrid& grid = nester_parts[i].grid;
            float hw = grid.nx * grid.voxel_size * 0.5f; // half-width
            float hh = grid.ny * grid.voxel_size * 0.5f; // half-height
            float cos_r = std::abs(std::cos(pl.zrot));
            float sin_r = std::abs(std::sin(pl.zrot));
            float rot_hw = hw * cos_r + hh * sin_r; // rotated half-width
            float rot_hh = hw * sin_r + hh * cos_r; // rotated half-height
            float cx = pl.x + grid.origin.x + hw;   // center X
            float cy = pl.y + grid.origin.y + hh;   // center Y

            bool within_bed = ((cx - rot_hw) >= -1.0f && (cy - rot_hh) >= -1.0f &&
                               (cx + rot_hw) <= bed_width_mm + 1.0f &&
                               (cy + rot_hh) <= bed_height_mm + 1.0f);

            Vec2crd out_translation;
            double out_rotation;
            snuggle_xform::placement_to_arrange(
                pl.x, pl.y, pl.zrot,
                bed_origin_x_mm, bed_origin_y_mm,
                out_translation, out_rotation);

            items[i].translation = out_translation;
            items[i].rotation    = out_rotation;

            if (within_bed) {
                items[i].bed_idx = 0;
                BOOST_LOG_TRIVIAL(debug) << "[SnuggleArrange] Part " << i
                    << " '" << part_data[i].name << "' placed on bed 0 (within bounds)";
            } else {
                items[i].bed_idx = UNARRANGED;
                BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Part " << i
                    << " '" << part_data[i].name << "' UNARRANGED (out of bed bounds)"
                    << " pos=(" << pl.x << "," << pl.y << ") mm"
                    << " center=(" << cx << "," << cy << ") rot_half=(" << rot_hw << "," << rot_hh << ")"
                    << " bed=[0,0]->[" << bed_width_mm << "," << bed_height_mm << "]";
            }
        }
    }

    // ── Timing summary ───────────────────────────────────────────────
    auto t_end = std::chrono::steady_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] === Snuggle nesting complete ==="
        << " total=" << total_ms << " ms"
        << " (voxelize=" << vox_ms << " ms, nest=" << nest_ms << " ms)";

    size_t placed = 0, unarranged = 0;
    for (const auto& item : items) {
        if (item.bed_idx >= 0) ++placed;
        else ++unarranged;
    }
    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Result: "
        << placed << " placed, " << unarranged << " unarranged";
}

}} // namespace Slic3r::arrangement
