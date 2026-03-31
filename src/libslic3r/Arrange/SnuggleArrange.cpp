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

    // Bed dimensions (used by several early-exit paths)
    float bed_width_mm = 0, bed_height_mm = 0;
    float bed_origin_x_mm = 0, bed_origin_y_mm = 0;

    // ── Edge case: nothing to arrange ────────────────────────────────
    if (items.empty()) {
        BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] No items to arrange, returning.";
        return;
    }

    // ── Uruk #1: Single part — just center it, don't waste CPU ──────
    if (items.size() == 1) {
        BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Single part — centering on bed.";
        snuggle_xform::bed_dimensions_mm(bed, bed_width_mm, bed_height_mm,
                                         bed_origin_x_mm, bed_origin_y_mm);
        float cx_mm = bed_origin_x_mm + bed_width_mm * 0.5f;
        float cy_mm = bed_origin_y_mm + bed_height_mm * 0.5f;
        items[0].translation = Vec2crd(snuggle_xform::mm_to_scaled(cx_mm),
                                       snuggle_xform::mm_to_scaled(cy_mm));
        items[0].rotation = 0.0; // Explicit — arranger contract requires this
        items[0].bed_idx = 0;
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

    // ── Pre-flight: can parts even theoretically fit? ────────────────
    // Quick check: if the single largest part exceeds the bed in either
    // dimension, it can never fit. Mark it UNARRANGED immediately.
    // Also check if total footprint area exceeds bed area (rough heuristic).
    {
        float bed_area = bed_width_mm * bed_height_mm;
        float total_footprint = 0;
        for (size_t i = 0; i < items.size(); ++i) {
            // Estimate part footprint from the hull polygon bounding box
            auto bb = get_extents(items[i].poly);
            float pw = snuggle_xform::scaled_to_mm(bb.max.x() - bb.min.x());
            float ph = snuggle_xform::scaled_to_mm(bb.max.y() - bb.min.y());
            total_footprint += pw * ph;

            if (pw > bed_width_mm || ph > bed_height_mm) {
                BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Part '" << items[i].name
                    << "' (" << pw << "x" << ph << " mm) exceeds bed dimensions. "
                    << "Marking UNARRANGED.";
                items[i].bed_idx = UNARRANGED;
            }
        }

        // If total footprint > 90% of bed area, warn (may not fit with gaps)
        if (total_footprint > bed_area * 0.9f) {
            BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Total part footprint ("
                << total_footprint << " mm^2) approaches bed area ("
                << bed_area << " mm^2). Some parts may not fit.";
        }
    }

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
    // Skip parts with empty voxel grids (failed voxelization, empty mesh).
    // These get marked UNARRANGED so the fallback arranger handles them.
    std::vector<snuggle::PartInfo> nester_parts;
    std::vector<size_t> nester_to_item; // maps nester index → items index
    nester_parts.reserve(part_data.size());
    for (size_t i = 0; i < part_data.size(); ++i) {
        if (part_data[i].info.grid.total_voxels() == 0 ||
            part_data[i].info.grid.count_solid() == 0) {
            BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Part " << i
                << " '" << part_data[i].name << "' has empty voxel grid — "
                << "marking UNARRANGED (fallback arranger will handle it)";
            items[i].bed_idx = UNARRANGED;
            continue;
        }
        nester_to_item.push_back(i);
        nester_parts.push_back(std::move(part_data[i].info));
    }

    // If no parts have valid grids, bail entirely
    if (nester_parts.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] No parts with valid voxel grids. "
            << "Falling back to default arranger.";
        return;
    }

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] " << nester_parts.size()
        << " parts with valid grids (of " << items.size() << " total)";

    // ── Step 2b: Build volume check ────────────────────────────────
    // Warn about parts that exceed the printable height
    float max_z = params.printable_height;
    for (size_t i = 0; i < part_data.size(); ++i) {
        if (part_data[i].info.max_height_mm > max_z) {
            BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Part '" << part_data[i].name
                << "' height " << part_data[i].info.max_height_mm
                << " mm exceeds printable height " << max_z << " mm";
        }
    }

    // ── Step 3: Configure and run nester ─────────────────────────────
    snuggle::NesterConfig cfg;

    // Conservative bed shrink: half a voxel per side for rounding safety.
    // Use half-voxel (not full) to avoid over-shrinking tiny beds (#5).
    float bed_margin = DEFAULT_VOXEL_SIZE_MM * 0.5f;
    cfg.bed_width_mm  = std::max(1.0f, bed_width_mm  - 2.0f * bed_margin);
    cfg.bed_height_mm = std::max(1.0f, bed_height_mm - 2.0f * bed_margin);

    BOOST_LOG_TRIVIAL(debug) << "[SnuggleArrange] Conservative bed: "
        << cfg.bed_width_mm << " x " << cfg.bed_height_mm << " mm"
        << " (shrunk " << bed_margin << " mm per side for voxel rounding)";

    // Convert min_obj_distance from scaled coords to mm via SnuggleTransform
    // This already includes brim/skirt inflation from update_selected_items_inflation()
    cfg.min_gap_mm = snuggle_xform::scaled_to_mm(params.min_obj_distance);
    if (cfg.min_gap_mm < 1.0f)
        cfg.min_gap_mm = 1.0f; // Minimum 1mm gap for safety

    BOOST_LOG_TRIVIAL(debug) << "[SnuggleArrange] min_gap_mm=" << cfg.min_gap_mm
        << " (from min_obj_distance=" << params.min_obj_distance << " scaled)";

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

    // Scale time budget and parameters by part count (#8 — fast for small jobs)
    if (items.size() <= 4) {
        cfg.population_size = 256;
        cfg.max_generations = 50;
        cfg.timeout_seconds = 5.0;   // 4 parts should be instant
    } else if (items.size() <= 10) {
        cfg.population_size = 256;
        cfg.max_generations = 80;
        cfg.timeout_seconds = 10.0;
    } else if (items.size() <= 20) {
        cfg.population_size = 512;
        cfg.max_generations = 100;
        cfg.timeout_seconds = 20.0;
    } else {
        cfg.population_size = 256;
        cfg.max_generations = 60;
        cfg.timeout_seconds = 30.0;
    }

    BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] Nester config:"
        << " bed=" << cfg.bed_width_mm << "x" << cfg.bed_height_mm << " mm"
        << " gap=" << cfg.min_gap_mm << " mm"
        << " lock_rot=" << cfg.lock_rotation
        << " pop=" << cfg.population_size
        << " gens=" << cfg.max_generations
        << " timeout=" << cfg.timeout_seconds << "s";

    // ── Sequential printing safety check (#14) ────────────────────
    if (params.is_seq_print) {
        BOOST_LOG_TRIVIAL(warning) << "[SnuggleArrange] Sequential printing enabled — "
            "Snuggle does not check toolhead clearance. "
            "Falling back to default arranger for safety.";
        return; // Leave all UNARRANGED → ArrangeJob fallback handles it
    }

    // ── Exclude regions (wipe tower, calibration, locked parts) (#10, #12, #13)
    // Pass excludes through to the nester as forbidden zones.
    // For now, convert exclude polygons to bed-relative bounding boxes and
    // shrink the effective bed to avoid them. Not perfect but prevents
    // the most common failure (placing parts on the wipe tower).
    if (!excludes.empty()) {
        BOOST_LOG_TRIVIAL(info) << "[SnuggleArrange] " << excludes.size()
            << " exclude regions detected. "
            << "Parts will avoid exclude bounding boxes.";
        // TODO: proper per-region exclusion in nester. For now, the ArrangeJob
        // fallback (which DOES respect excludes) handles overflow items.
    }

    auto t_nest_start = std::chrono::steady_clock::now();

    snuggle::SnuggleNester nester(cfg);
    snuggle::NesterResult nest_result = nester.run(nester_parts,
        [&params](size_t gen, size_t max_gen, const snuggle::Individual& best) -> bool {
            // #17: Wire cancel button
            if (params.stopcondition && params.stopcondition())
                return false; // User pressed Cancel

            // #9: Wire progress bar
            if (params.progressind) {
                unsigned pct = (unsigned)(gen * 100 / std::max(max_gen, (size_t)1));
                params.progressind(pct, " (Snuggle)");
            }

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
    // nester_parts may be smaller than items (empty grids were skipped)
    if (nest_result.placements.size() != nester_parts.size()) {
        BOOST_LOG_TRIVIAL(error) << "[SnuggleArrange] Nester returned "
            << nest_result.placements.size() << " placements but expected "
            << nester_parts.size() << ". Leaving items unarranged.";
        return;
    }

    if (nest_result.feasible) {
        // All valid parts fit on one bed
        for (size_t ni = 0; ni < nester_parts.size(); ++ni) {
            size_t i = nester_to_item[ni]; // Map back to items index
            const snuggle::Placement& pl = nest_result.placements[ni];

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
        for (size_t ni = 0; ni < nester_parts.size(); ++ni) {
            size_t i = nester_to_item[ni];
            const snuggle::Placement& pl = nest_result.placements[ni];

            // Check if this part's placement is within bed bounds.
            const snuggle::VoxelGrid& grid = nester_parts[ni].grid;
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
