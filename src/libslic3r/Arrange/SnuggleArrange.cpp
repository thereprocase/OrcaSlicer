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

namespace Slic3r { namespace arrangement {

void snuggle_arrange(
    ArrangePolygons&          items,
    const ArrangePolygons&    excludes,
    const Points&             bed,
    const ArrangeParams&      params,
    const Model&              model)
{
    BOOST_LOG_TRIVIAL(warning) << "Snuggle: starting 3D-aware arrangement for " << items.size() << " items";

    if (items.empty()) return;

    if (!excludes.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: " << excludes.size()
            << " excluded items (locked parts, wipe tower) will be IGNORED — not yet implemented";
    }

    // ── Determine bed dimensions from bed points ───────────
    BoundingBox bed_bb(bed);
    float bed_w = unscale_(bed_bb.max.x() - bed_bb.min.x());
    float bed_h = unscale_(bed_bb.max.y() - bed_bb.min.y());
    float bed_origin_x = unscale_(bed_bb.min.x());
    float bed_origin_y = unscale_(bed_bb.min.y());

    BOOST_LOG_TRIVIAL(warning) << "Snuggle: bed " << bed_w << " x " << bed_h
                               << " mm, origin (" << bed_origin_x << ", " << bed_origin_y << ")";

    // ── Voxelize each part ─────────────────────────────────
    // At 2mm voxels, a 100mm part is a 50x50 grid — triangle count
    // doesn't matter much since we iterate the grid, not the mesh.
    // But high-poly meshes (>10K tris) make the SAT test expensive.
    // Decimate to ~5K tris max — sufficient for 2mm voxel accuracy.
    float voxel_size = 2.0f;
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
                << items[i].name << "', skipping voxelization";
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

        // Clip at Z=0: if user sank the part into the bed, only voxelize
        // the portion above the bed surface. Clamp vertices below Z=0.
        for (auto& v : mesh.its.vertices)
            if (v.z() < 0) v.z() = 0;

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

        if (err != snuggle::VoxError::OK) {
            BOOST_LOG_TRIVIAL(warning) << "Snuggle: voxelization failed for "
                << pi.name << ": " << snuggle::vox_error_str(err);
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

    // Quality 1-10 maps to population and generations
    int quality = std::clamp(params.snuggle_quality, 1, 10);
    cfg.population_size = 64 + quality * 64;    // 128 to 704
    cfg.max_generations = 20 + quality * 10;    // 30 to 120
    cfg.timeout_seconds = 5.0 + quality * 5.0;  // 10s to 55s

    cfg.lock_rotation   = params.snuggle_lock_rotation;

    // Shrink effective bed by voxel padding (1 voxel per side) + min gap
    cfg.bed_margin_mm   = voxel_size + cfg.min_gap_mm;

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
            unsigned progress = (unsigned)(gen * items.size() / max_gen);
            std::string status = " (Snuggle " + backend_name
                + " gen " + std::to_string(gen) + "/" + std::to_string(max_gen)
                + (best.is_feasible() ? " OK" : " " + std::to_string(best.collision_count) + " collisions")
                + ")";
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
            const auto& pl = result.placements[i];
            const auto& grid_i = parts[i].grid;
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
            << " rot=" << (int)(pl.zrot * 180.0 / 3.14159265) << "deg";
    }

    if (rejected > 0)
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: " << placed << " placed, "
                                   << rejected << " off-plate (didn't fit)";
}

}} // namespace Slic3r::arrangement
