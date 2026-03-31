// SnuggleArrange.cpp — Snuggle 3D-aware arrangement for OrcaSlicer
//
// Loads 3D meshes from the Model, voxelizes them, runs the genetic nester,
// and writes results back to ArrangePolygons for the standard finalize path.

#include "SnuggleArrange.hpp"
#include "polite_voxelizer.hpp"
#include "gpu_collision.hpp"
#include "snuggle_nester.hpp"

#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/QuadricEdgeCollapse.hpp"

#include <boost/log/trivial.hpp>
#include <map>

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
    float base_voxel_size = std::clamp(params.snuggle_voxel_mm, 0.5f, 5.0f);
    constexpr size_t MAX_TRIS_FOR_VOXEL = 5000;
    // Adaptive voxel sizing: if a part's grid would exceed safe dimensions
    // at the user's voxel size, coarsen that part until it fits.
    // MAX_GRID_DIM is 512 per axis. Rotation expands by ~sqrt(2), so cap
    // at 350 to leave room for rotated_copy.
    constexpr size_t SAFE_GRID_DIM = 350;
    // Total rotation cache budget: 512MB. Skip cache for parts that don't fit.
    constexpr size_t MAX_CACHE_BYTES = 512ULL * 1024 * 1024;

    // Walk the model to match items to instances. Build a flat list of
    // (object, instance) pairs in the same order prepare_all() enumerates,
    // then match by name. This gives us the ModelInstance pointer needed
    // to apply the correct instance transform.
    struct ObjInst { const ModelObject* obj; const ModelInstance* inst; };
    std::multimap<std::string, ObjInst> instance_lookup;
    for (const ModelObject* obj : model.objects)
        for (const ModelInstance* inst : obj->instances)
            instance_lookup.emplace(obj->name, ObjInst{obj, inst});

    std::vector<snuggle::PartInfo> parts;
    auto t_vox_start = std::chrono::steady_clock::now();

    for (size_t i = 0; i < items.size(); i++) {
        snuggle::PartInfo pi;
        pi.name = items[i].name;
        pi.initial_zrot = (float)items[i].rotation;

        // Find matching instance by name (first unused match, O(log N) lookup)
        const ModelObject* obj = nullptr;
        const ModelInstance* inst = nullptr;
        {
            auto it = instance_lookup.find(items[i].name);
            if (it != instance_lookup.end()) {
                obj = it->second.obj;
                inst = it->second.inst;
                instance_lookup.erase(it);  // consume this instance
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

        // Adaptive voxel size: compute mesh AABB, coarsen if grid too large
        float part_voxel = base_voxel_size;
        {
            float mesh_max_dim = 0;
            for (size_t vi = 0; vi < its.vertices.size(); vi++) {
                mesh_max_dim = std::max(mesh_max_dim, std::abs(its.vertices[vi].x()));
                mesh_max_dim = std::max(mesh_max_dim, std::abs(its.vertices[vi].y()));
                mesh_max_dim = std::max(mesh_max_dim, std::abs(its.vertices[vi].z()));
            }
            float extent = mesh_max_dim * 2.0f + part_voxel * 4.0f; // diameter + padding
            while (extent / part_voxel > SAFE_GRID_DIM && part_voxel < 10.0f) {
                part_voxel *= 1.5f;
            }
            if (part_voxel > base_voxel_size) {
                BOOST_LOG_TRIVIAL(warning) << "Snuggle: coarsened " << pi.name
                    << " voxel from " << base_voxel_size << "mm to " << part_voxel
                    << "mm (mesh too large for " << SAFE_GRID_DIM << " grid)";
            }
        }

        auto err = snuggle::voxelize_indexed_mesh(
            verts.data(), its.vertices.size(),
            indices.data(), its.indices.size(),
            part_voxel, pi.grid);

        if (err != snuggle::VoxError::OK || pi.grid.count_solid() == 0) {
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: voxelization failed for '"
                << pi.name << "' — marking off-plate, fallback will handle it";
            items[i].bed_idx = -1; // Let overflow fallback handle this part
            pi.max_height_mm = 0;
            pi.hull_area_mm2 = 0;
        } else {
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

    // ── Configure and run the genetic nester ───────────────
    snuggle::NesterConfig cfg;
    cfg.bed_width_mm    = bed_w;
    cfg.bed_height_mm   = bed_h;
    cfg.min_gap_mm      = std::max(1.0f, params.snuggle_padding_mm);

    // Direct numeric controls — user sets population, generations, timeout
    cfg.population_size = std::clamp(params.snuggle_population, 16, 1024);
    cfg.max_generations = std::clamp(params.snuggle_generations, 10, 500);
    cfg.timeout_seconds = std::max(2.0, (double)params.snuggle_timeout_s);

    // Rotation step: 0=locked, else degrees per snap increment
    // The initial_zrot on each part is the user's pre-rotation (their starting position)
    if (params.snuggle_rotation_step <= 0) {
        cfg.lock_rotation = true;
    } else {
        cfg.lock_rotation = false;
        // Quantize allowed rotations to the step size
        // The nester's rotation_step_rad limits what angles the GA explores
        cfg.rotation_step_rad = (float)params.snuggle_rotation_step * (snuggle::PI_F / 180.0f);
    }
    // Lock rotation override from the checkbox takes priority
    if (params.snuggle_lock_rotation)
        cfg.lock_rotation = true;

    cfg.compact         = params.snuggle_compact;

    BOOST_LOG_TRIVIAL(info) << "Snuggle config:"
        << " pop=" << cfg.population_size << " gens=" << cfg.max_generations
        << " timeout=" << cfg.timeout_seconds << "s"
        << " voxel=" << voxel_size << "mm"
        << " rot_step=" << params.snuggle_rotation_step << "deg"
        << " lock_rot=" << cfg.lock_rotation;

    // Margin = clearance from grid boundary to bed edge.
    // The rotated grid already includes 1-voxel padding beyond geometry,
    // so the margin only needs the gap clearance. Using max_voxel + gap
    // double-counted the padding and wasted bed space (8mm total with 4mm voxels),
    // making tight layouts infeasible with adaptive per-part voxel sizes.
    cfg.bed_margin_mm   = cfg.min_gap_mm;

    // Wire progress/stop to Orca's callbacks
    if (params.stopcondition) {
        // Check periodically but don't make it the inner loop
        cfg.yield_every_gens = 1;
    }

    auto evaluator = snuggle::create_collision_evaluator();
    std::string backend_name = "CPU";

#ifdef SLIC3R_GUI
    if (auto* gpu = dynamic_cast<snuggle::GpuCollisionEvaluator*>(evaluator.get())) {
        if (gpu->is_available())
            backend_name = "GPU";
    }
#endif

    BOOST_LOG_TRIVIAL(warning) << "Snuggle: using " << backend_name << " collision backend";

    snuggle::SnuggleNester nester(cfg, evaluator.get());

    auto result = nester.run(parts, [&](size_t gen, size_t max_gen,
                                        const snuggle::Individual& best) -> bool {
        // Progress string includes backend, gen count, and collision status
        if (params.progressind) {
            unsigned progress = (unsigned)(gen * 100 / std::max(max_gen, (size_t)1));
            std::string status = best.is_feasible()
                ? " Optimizing... " + std::to_string(progress) + "% (all parts fit)"
                : " Optimizing... " + std::to_string(progress) + "% ("
                  + std::to_string(best.collision_count) + " parts overlapping)";
            params.progressind(progress, status);
        }
        // Check stop condition
        if (params.stopcondition && params.stopcondition())
            return false;
        return true;
    });

    BOOST_LOG_TRIVIAL(warning) << "Snuggle [" << backend_name << "]: "
        << parts.size() << " parts, "
        << result.time_ms << "ms, "
        << result.generations_run << " gens, "
        << "feasible=" << result.feasible
        << ", collisions=" << result.collisions
        << (result.timed_out ? ", TIMED OUT" : "");

    // ── Write results back to ArrangePolygons ──────────────
    // The nester resolved each part's instance origin position on the bed.
    // If the result is feasible, place everything. If not, greedily accept
    // parts that fit and send the rest off-plate (bed_idx = -1).
    int placed = 0, rejected = 0;

    if (result.feasible) {
        // All fit — place everything
        for (size_t i = 0; i < items.size() && i < result.placements.size(); i++) {
            const auto& pl = result.placements[i];
            items[i].translation = Vec2crd(
                scaled(pl.origin_bed_x + bed_origin_x),
                scaled(pl.origin_bed_y + bed_origin_y)
            );
            items[i].rotation = (double)pl.zrot;
            items[i].bed_idx = 0;
            placed++;
        }
    } else {
        // Infeasible — greedy accept: place each part if it doesn't collide
        // with already-accepted parts and stays within bed bounds.
        float margin = cfg.bed_margin_mm;
        std::vector<size_t> accepted;

        for (size_t i = 0; i < items.size() && i < result.placements.size(); i++) {
            // Skip parts with empty grids (failed voxelization, already marked off-plate)
            if (items[i].bed_idx == -1) { rejected++; continue; }

            const auto& pl = result.placements[i];
            const auto& grid_i = parts[i].grid;
            if (grid_i.total_voxels() == 0) { items[i].bed_idx = -1; rejected++; continue; }
            auto rot_i = grid_i.rotated_copy(pl.zrot);

            // Bounds check
            float pmin_x = pl.x + rot_i.origin.x;
            float pmin_y = pl.y + rot_i.origin.y;
            float pmax_x = pmin_x + rot_i.nx * rot_i.voxel_size;
            float pmax_y = pmin_y + rot_i.ny * rot_i.voxel_size;
            bool in_bounds = (pmin_x >= margin && pmin_y >= margin &&
                              pmax_x <= bed_w - margin && pmax_y <= bed_h - margin);

            // Collision check against accepted parts
            bool collides = false;
            if (in_bounds) {
                snuggle::Vec3f off_i = {pl.x, pl.y, 0.0f};
                for (size_t j : accepted) {
                    const auto& pl_j = result.placements[j];
                    auto rot_j = parts[j].grid.rotated_copy(pl_j.zrot);
                    snuggle::Vec3f off_j = {pl_j.x, pl_j.y, 0.0f};
                    if (snuggle::VoxelGrid::collision_count(rot_i, off_i, rot_j, off_j) > 0) {
                        collides = true;
                        break;
                    }
                }
            }

            if (in_bounds && !collides) {
                items[i].translation = Vec2crd(
                    scaled(pl.origin_bed_x + bed_origin_x),
                    scaled(pl.origin_bed_y + bed_origin_y)
                );
                items[i].rotation = (double)pl.zrot;
                items[i].bed_idx = 0;
                accepted.push_back(i);
                placed++;
            } else {
                items[i].bed_idx = -1;  // off-plate
                rejected++;
            }
        }
    }

    for (size_t i = 0; i < items.size() && i < result.placements.size(); i++) {
        if (items[i].bed_idx < 0) continue;
        const auto& pl = result.placements[i];
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: [" << i << "] " << items[i].name
            << " origin_bed=(" << pl.origin_bed_x << "," << pl.origin_bed_y << ")"
            << " rot=" << (int)(pl.zrot * 180.0 / snuggle::PI_F) << "deg";
    }

    if (rejected > 0)
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: " << placed << " placed, "
                                   << rejected << " off-plate (didn't fit)";
}

}} // namespace Slic3r::arrangement
