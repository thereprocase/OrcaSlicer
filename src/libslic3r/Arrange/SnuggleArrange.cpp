// SnuggleArrange.cpp — Snuggle 3D-aware arrangement for OrcaSlicer
//
// Loads 3D meshes from the Model, voxelizes them, runs the radial
// expansion nester, writes results back to ArrangePolygons.

#include "SnuggleArrange.hpp"
#include "polite_voxelizer.hpp"
#include "gpu_collision.hpp"
#include "snuggle_nester.hpp"
#include "snuggle_radial.hpp"

#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/QuadricEdgeCollapse.hpp"

#include <boost/log/trivial.hpp>

namespace Slic3r { namespace arrangement {

void snuggle_arrange(
    ArrangePolygons&          items,
    const ArrangePolygons&    excludes,
    const Points&             bed,
    const ArrangeParams&      params,
    const Model&              model)
{
    BOOST_LOG_TRIVIAL(info) << "Snuggle: starting 3D-aware arrangement for " << items.size() << " items";

    if (items.empty()) return;

    // ── Determine bed dimensions ──────────────────────────
    BoundingBox bed_bb(bed);
    float bed_w = unscale_(bed_bb.max.x() - bed_bb.min.x());
    float bed_h = unscale_(bed_bb.max.y() - bed_bb.min.y());
    float bed_origin_x = unscale_(bed_bb.min.x());
    float bed_origin_y = unscale_(bed_bb.min.y());

    // ── Single-part: center on bed accounting for polygon offset ──
    if (items.size() == 1) {
        BOOST_LOG_TRIVIAL(info) << "Snuggle: single part — centering on bed";
        // The polygon is centered at instance origin. Compute the polygon's
        // own centroid so the VISIBLE geometry lands at bed center.
        auto poly_bb = get_extents(items[0].poly);
        float poly_cx = unscale_(poly_bb.min.x() + poly_bb.max.x()) * 0.5f;
        float poly_cy = unscale_(poly_bb.min.y() + poly_bb.max.y()) * 0.5f;
        float bed_cx = bed_origin_x + bed_w * 0.5f;
        float bed_cy = bed_origin_y + bed_h * 0.5f;
        // Place instance origin so that geometry centroid lands at bed center
        items[0].translation = Vec2crd(scaled(bed_cx - poly_cx), scaled(bed_cy - poly_cy));
        // Preserve user's Z rotation
        items[0].bed_idx = 0;
        return;
    }

    // ── Too many parts: fall back to standard arranger ────
    size_t max_parts = std::max(2, params.snuggle_max_parts);
    if (items.size() > max_parts) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: " << items.size() << " parts exceeds limit ("
            << max_parts << "). Falling back to standard arranger.";
        // Mark all UNARRANGED so ArrangeJob overflow fallback triggers
        for (auto& item : items) item.bed_idx = -1;
        return;
    }

    // ── Sequential printing: hard block ───────────────────
    if (params.is_seq_print) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: sequential printing active — "
            "cannot verify toolhead clearance. Falling back to standard arranger.";
        for (auto& item : items) item.bed_idx = -1;
        return;
    }

    // ── Pre-flight: flag oversized parts ──────────────────
    for (size_t i = 0; i < items.size(); ++i) {
        auto bb = get_extents(items[i].poly);
        float pw = unscale_(bb.max.x() - bb.min.x());
        float ph = unscale_(bb.max.y() - bb.min.y());
        if (pw > bed_w || ph > bed_h) {
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: part '" << items[i].name
                << "' (" << pw << "x" << ph << " mm) exceeds bed. Marking off-plate.";
            items[i].bed_idx = -1;
        }
    }

    if (!excludes.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: " << excludes.size()
            << " excluded items (locked parts, wipe tower) — passed as obstacles to overflow arranger";
    }

    BOOST_LOG_TRIVIAL(warning) << "Snuggle: bed " << bed_w << " x " << bed_h
                               << " mm, origin (" << bed_origin_x << ", " << bed_origin_y << ")";

    // ── Voxelize each part ─────────────────────────────────
    // At 2mm voxels, a 100mm part is a 50x50 grid — triangle count
    // doesn't matter much since we iterate the grid, not the mesh.
    // But high-poly meshes (>10K tris) make the SAT test expensive.
    // Decimate to ~5K tris max — sufficient for 2mm voxel accuracy.
    float voxel_size = std::clamp(params.snuggle_voxel_mm, 0.5f, 5.0f);
    constexpr size_t MAX_TRIS_FOR_VOXEL = 5000;

    // Walk the model to match items to instances. Build a flat list of
    // (object, instance) pairs in the same order prepare_all() enumerates,
    // then match by name. This gives us the ModelInstance pointer needed
    // to apply the correct instance transform.
    struct ObjInst { const ModelObject* obj; const ModelInstance* inst; };
    std::vector<ObjInst> all_instances;
    for (const ModelObject* obj : model.objects)
        for (const ModelInstance* inst : obj->instances)
            all_instances.push_back({obj, inst});

    std::vector<snuggle::PartInfo> parts;
    auto t_vox_start = std::chrono::steady_clock::now();

    for (size_t i = 0; i < items.size(); i++) {
        snuggle::PartInfo pi;
        pi.name = items[i].name;
        pi.initial_zrot = (float)items[i].rotation;

        // Find matching instance by name (first unused match)
        const ModelObject* obj = nullptr;
        const ModelInstance* inst = nullptr;
        for (auto& oi : all_instances) {
            if (oi.obj && oi.obj->name == items[i].name) {
                obj = oi.obj;
                inst = oi.inst;
                oi.obj = nullptr;  // mark used
                break;
            }
        }

        if (!obj || !inst) {
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: no model instance for item '"
                << items[i].name << "', marking off-plate";
            items[i].bed_idx = -1; // Let overflow fallback handle it
            pi.max_height_mm = 0;
            pi.hull_area_mm2 = 0;
            parts.push_back(std::move(pi));
            continue;
        }

        // Apply the SAME transform as get_arrange_polygon():
        // Instance scale + mirror + X/Y rotation + Z offset
        // but NOT XY offset (nester handles placement) or Z rotation (nester handles rotation).
        // This matches the ArrangePolygon's polygon coordinate frame exactly.
        Geometry::Transformation t(inst->get_transformation());
        t.set_offset(Vec3d(0, 0, t.get_offset().z()));  // keep Z offset, zero XY
        Vec3d rot = t.get_rotation();
        rot.z() = 0;
        t.set_rotation(rot);  // keep X/Y rotation, zero Z

        TriangleMesh mesh = obj->raw_mesh();
        mesh.transform(t.get_matrix());

        // Bed clipping: handled in voxelizer — voxels below Z=0 are
        // naturally empty because no triangles exist there after the
        // grid origin is clamped to Z >= 0. We just need to ensure
        // the grid doesn't extend below the bed surface.
        // (Vertex clamping would squish geometry, making parts fatter
        // than their actual cross-section at the bed plane.)

        size_t orig_tris = mesh.its.indices.size();
        if (orig_tris > MAX_TRIS_FOR_VOXEL) {
            its_quadric_edge_collapse(mesh.its, (uint32_t)MAX_TRIS_FOR_VOXEL);
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: decimated " << pi.name
                << " from " << orig_tris << " to " << mesh.its.indices.size() << " tris";
        }

        const auto& its = mesh.its;

        std::vector<float> verts(its.vertices.size() * 3);
        std::vector<uint32_t> indices(its.indices.size() * 3);

        for (size_t vi = 0; vi < its.vertices.size(); vi++) {
            verts[vi * 3 + 0] = its.vertices[vi].x();
            verts[vi * 3 + 1] = its.vertices[vi].y();
            verts[vi * 3 + 2] = its.vertices[vi].z();
        }
        for (size_t ti = 0; ti < its.indices.size(); ti++) {
            indices[ti * 3 + 0] = its.indices[ti](0);
            indices[ti * 3 + 1] = its.indices[ti](1);
            indices[ti * 3 + 2] = its.indices[ti](2);
        }

        auto err = snuggle::voxelize_indexed_mesh(
            verts.data(), its.vertices.size(),
            indices.data(), its.indices.size(),
            voxel_size, pi.grid);

        if (err != snuggle::VoxError::OK || pi.grid.count_solid() == 0) {
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: voxelization failed for '"
                << pi.name << "' — marking off-plate, fallback will handle it";
            items[i].bed_idx = -1; // Let overflow fallback handle this part
            pi.max_height_mm = 0;
            pi.hull_area_mm2 = 0;
            parts.push_back(std::move(pi));
            continue;
        } {
            pi.max_height_mm = pi.grid.nz * pi.grid.voxel_size;
            pi.hull_area_mm2 = pi.grid.nx * pi.grid.ny * voxel_size * voxel_size;
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: [" << i << "] " << pi.name
                << " (" << its.indices.size() << " tris) -> "
                << pi.grid.nx << "x" << pi.grid.ny << "x" << pi.grid.nz
                << " grid_origin=(" << pi.grid.origin.x << "," << pi.grid.origin.y << ")";
        }

        parts.push_back(std::move(pi));
    }

    auto t_vox_end = std::chrono::steady_clock::now();
    double vox_ms = std::chrono::duration<double, std::milli>(t_vox_end - t_vox_start).count();
    BOOST_LOG_TRIVIAL(warning) << "Snuggle: voxelized " << parts.size()
                               << " parts in " << (int)vox_ms << "ms";

    // ── Configure radial expansion nester ───────────────────
    snuggle::RadialConfig rcfg;
    rcfg.bed_width_mm  = bed_w;
    rcfg.bed_height_mm = bed_h;
    rcfg.min_gap_mm    = std::max(1.5f, params.snuggle_padding_mm);
    rcfg.bed_margin_mm = rcfg.min_gap_mm;
    rcfg.step_mm       = voxel_size;
    rcfg.timeout_s     = std::max(2.0, (double)params.snuggle_timeout_s);

    if (params.snuggle_rotation_step <= 0 || params.snuggle_lock_rotation) {
        rcfg.lock_rotation = true;
        rcfg.n_rotations = 1;
    } else {
        rcfg.lock_rotation = false;
        rcfg.n_rotations = std::max(1, 360 / params.snuggle_rotation_step);
    }

    BOOST_LOG_TRIVIAL(warning) << "Snuggle radial: " << parts.size() << " parts"
        << ", " << rcfg.n_directions << " dirs x " << rcfg.n_rotations << " rots"
        << ", step=" << rcfg.step_mm << "mm, margin=" << rcfg.bed_margin_mm << "mm"
        << ", lock_rot=" << rcfg.lock_rotation;

    // Build rotation cache for radial nester
    std::vector<std::vector<snuggle::VoxelGrid>> rot_cache(parts.size());
    for (size_t i = 0; i < parts.size(); i++) {
        if (rcfg.lock_rotation) {
            snuggle::VoxelGrid single = parts[i].grid.rotated_copy(parts[i].initial_zrot);
            rot_cache[i].resize(snuggle::ROT_CACHE_BINS, single);
        } else {
            rot_cache[i].resize(snuggle::ROT_CACHE_BINS);
            for (int bin = 0; bin < snuggle::ROT_CACHE_BINS; bin++) {
                float angle = (float)bin * snuggle::TWO_PI_F / snuggle::ROT_CACHE_BINS;
                rot_cache[i][bin] = parts[i].grid.rotated_copy(angle);
            }
        }
        if (params.stopcondition && params.stopcondition()) break;
    }

    // Create collision evaluator — GPU if enabled and available, else CPU
    std::unique_ptr<snuggle::CollisionEvaluator> evaluator;
    if (params.snuggle_use_gpu) {
        evaluator = snuggle::create_collision_evaluator(); // tries GPU, falls back to CPU
    } else {
        evaluator = std::make_unique<snuggle::CpuCollisionEvaluator>();
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: GPU disabled by user, using CPU evaluator";
    }
    evaluator->upload_grids(parts, rot_cache);

    auto result = snuggle::radial_arrange(parts, rot_cache, rcfg,
        evaluator.get(),
        [&](int placed, int total, const char* name) -> bool {
            if (params.progressind) {
                unsigned pct = total > 0 ? (placed * 100 / total) : 0;
                std::string status = " (Placing " + std::to_string(placed)
                    + "/" + std::to_string(total) + " " + name + ")";
                params.progressind(pct, status);
            }
            return !(params.stopcondition && params.stopcondition());
        });

    BOOST_LOG_TRIVIAL(warning) << "Snuggle radial: "
        << result.placed_count << "/" << result.total_count << " placed in "
        << (int)result.time_ms << "ms"
        << (result.timed_out ? " (TIMED OUT)" : "");

    // ── Write results back to ArrangePolygons ──────────────
    // Radial nester already validated each placement. Parts marked
    // placed=true are collision-free and within bed bounds.
    int placed = 0, rejected = 0;

    for (size_t i = 0; i < items.size() && i < result.placements.size(); i++) {
        if (items[i].bed_idx == -1) { rejected++; continue; }

        const auto& pl = result.placements[i];
        if (pl.placed) {
            items[i].translation = Vec2crd(
                scaled(pl.origin_bed_x + bed_origin_x),
                scaled(pl.origin_bed_y + bed_origin_y)
            );
            items[i].rotation = (double)pl.zrot;
            items[i].bed_idx = 0;
            placed++;

            BOOST_LOG_TRIVIAL(warning) << "Snuggle: [" << i << "] " << items[i].name
                << " origin_bed=(" << pl.origin_bed_x << "," << pl.origin_bed_y << ")"
                << " rot=" << (int)(pl.zrot * 180.0f / snuggle::PI_F) << "deg";
        } else {
            items[i].bed_idx = -1;
            rejected++;
        }
    }

    if (rejected > 0)
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: " << placed << " placed, "
                                   << rejected << " off-plate (didn't fit)";
}

}} // namespace Slic3r::arrangement
