// SnuggleArrange.cpp — Snuggle 3D-aware arrangement for OrcaSlicer
//
// Loads 3D meshes from the Model, voxelizes them, runs the genetic nester,
// and writes results back to ArrangePolygons for the standard finalize path.

#include "SnuggleArrange.hpp"
#include "polite_voxelizer.hpp"
#include "snuggle_nester.hpp"
#include "gpu_collision.hpp"

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

    std::vector<snuggle::PartInfo> parts;
    size_t instance_idx = 0;
    auto t_vox_start = std::chrono::steady_clock::now();

    for (const ModelObject* obj : model.objects) {
        for (const ModelInstance* inst : obj->instances) {
            if (instance_idx >= items.size()) break;

            snuggle::PartInfo pi;
            pi.name = obj->name;
            pi.initial_zrot = (float)inst->get_rotation().z();

            TriangleMesh mesh = obj->raw_mesh();

            // Decimate high-poly meshes for faster voxelization.
            // At 2mm voxels, sub-2mm mesh detail is invisible. 5K tris is plenty.
            size_t orig_tris = mesh.its.indices.size();
            if (orig_tris > MAX_TRIS_FOR_VOXEL) {
                its_quadric_edge_collapse(mesh.its, (uint32_t)MAX_TRIS_FOR_VOXEL);
                BOOST_LOG_TRIVIAL(warning) << "Snuggle: decimated " << obj->name
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
                    << obj->name << ": " << snuggle::vox_error_str(err);
                pi.max_height_mm = 0;
                pi.hull_area_mm2 = 0;
            } else {
                pi.max_height_mm = pi.grid.nz * pi.grid.voxel_size;
                pi.hull_area_mm2 = pi.grid.nx * pi.grid.ny * voxel_size * voxel_size;
                BOOST_LOG_TRIVIAL(warning) << "Snuggle: voxelized " << obj->name
                    << " (" << its.indices.size() << " tris) -> "
                    << pi.grid.nx << "x" << pi.grid.ny << "x" << pi.grid.nz
                    << " (" << pi.grid.count_solid() << " solid voxels)";
            }

            parts.push_back(std::move(pi));
            instance_idx++;
        }
    }

    auto t_vox_end = std::chrono::steady_clock::now();
    double vox_ms = std::chrono::duration<double, std::milli>(t_vox_end - t_vox_start).count();
    BOOST_LOG_TRIVIAL(warning) << "Snuggle: voxelized " << parts.size()
                               << " parts in " << (int)vox_ms << "ms";

    // ── Configure and run the genetic nester ───────────────
    snuggle::NesterConfig cfg;
    cfg.bed_width_mm    = bed_w;
    cfg.bed_height_mm   = bed_h;
    cfg.population_size = 256;
    cfg.max_generations = 50;
    cfg.timeout_seconds = 20.0;
    cfg.min_gap_mm      = std::max(1.0f, (float)unscale_(params.min_obj_distance));
    cfg.lock_rotation   = params.snuggle_lock_rotation;

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
    for (size_t i = 0; i < items.size() && i < result.placements.size(); i++) {
        const auto& pl = result.placements[i];

        // Convert from bed-relative mm to Orca's scaled coordinates
        // Snuggle places parts relative to bed (0,0), Orca expects absolute coordinates
        items[i].translation = Vec2crd(
            scaled(pl.x + bed_origin_x),
            scaled(pl.y + bed_origin_y)
        );
        items[i].rotation = (double)pl.zrot;
        items[i].bed_idx = 0; // All on first bed for now

        BOOST_LOG_TRIVIAL(debug) << "Snuggle: " << items[i].name
            << " -> (" << pl.x << ", " << pl.y << ") rot=" << (pl.zrot * 180.0 / 3.14159265) << "deg";
    }
}

}} // namespace Slic3r::arrangement
