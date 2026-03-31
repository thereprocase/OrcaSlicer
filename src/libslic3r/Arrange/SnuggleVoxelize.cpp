// SnuggleVoxelize.cpp — Voxelize OrcaSlicer model parts for Snuggle nesting

#include "SnuggleVoxelize.hpp"
#include "polite_voxelizer.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Geometry.hpp"
#include <boost/log/trivial.hpp>
#include <cmath>
#include <cassert>

namespace Slic3r { namespace arrangement {

// Helper: compute the same transform matrix that get_arrange_polygon uses.
// From Model.cpp:
//   rotation = instance rotation with Z zeroed
//   offset   = Z-only (XY zeroed)
//   Then: trafo = Transformation(offset, rotation, scale, ...).get_matrix()
//
// This ensures our voxel grid is in the SAME coordinate frame as the
// convex hull polygon that ArrangePolygon uses.
static Transform3d build_hull_transform(const ModelInstance& instance)
{
    // Replicate what get_arrange_polygon does:
    //   Vec3d rotation = get_rotation();
    //   rotation.z() = 0.;
    //   Geometry::Transformation t(m_transformation);
    //   t.set_offset(get_offset().z() * Vec3d::UnitZ());
    //   t.set_rotation(rotation);
    //   ... then uses t.get_matrix()

    Geometry::Transformation t(instance.get_transformation());

    // Zero XY offset — the nester will place XY
    t.set_offset(instance.get_offset().z() * Vec3d::UnitZ());

    // Keep XY rotation, zero Z rotation — stored separately
    Vec3d rot = instance.get_rotation();
    rot.z() = 0.0;
    t.set_rotation(rot);

    return t.get_matrix();
}

// Helper: collect all model-part volume meshes for an object, transformed
// into the hull coordinate frame, and flatten into a snuggle Triangle array.
static std::vector<snuggle::Triangle> collect_triangles(
    const ModelObject& object,
    const Transform3d& hull_trafo)
{
    std::vector<snuggle::Triangle> triangles;

    for (const ModelVolume* vol : object.volumes) {
        if (!vol->is_model_part())
            continue;

        const indexed_triangle_set& its = vol->mesh().its;
        if (its.indices.empty())
            continue;

        // Full transform: hull_trafo * volume_trafo
        // This matches convex_hull_2d which does:
        //   (trafo_instance * v->get_matrix()).cast<float>()
        Transform3d full_trafo = hull_trafo * vol->get_matrix();
        Eigen::Matrix3f rot_scale = full_trafo.matrix().block<3,3>(0,0).cast<float>();
        Eigen::Vector3f translation = full_trafo.matrix().block<3,1>(0,3).cast<float>();

        size_t base = triangles.size();
        triangles.resize(base + its.indices.size());

        for (size_t i = 0; i < its.indices.size(); ++i) {
            const Vec3i32& face = its.indices[i];

            // Transform each vertex: result is in mm (hull frame)
            auto xform_vert = [&](const Vec3f& v) -> snuggle::Vec3f {
                Eigen::Vector3f tv = rot_scale * v + translation;
                return { tv.x(), tv.y(), tv.z() };
            };

            triangles[base + i].v0 = xform_vert(its.vertices[face(0)]);
            triangles[base + i].v1 = xform_vert(its.vertices[face(1)]);
            triangles[base + i].v2 = xform_vert(its.vertices[face(2)]);
        }
    }

    return triangles;
}

std::vector<SnugglePartData> voxelize_for_snuggle(
    const Model& model,
    const ArrangePolygons& items,
    float voxel_size_mm,
    std::function<bool()> stop_check)
{
    BOOST_LOG_TRIVIAL(info) << "[SnuggleVoxelize] Starting voxelization of "
        << items.size() << " parts, voxel_size=" << voxel_size_mm << " mm";

    std::vector<SnugglePartData> result;
    result.reserve(items.size());

    // Build a flat list of (object, instance) pairs matching the ArrangePolygons order.
    // ArrangePolygons are built by iterating objects then instances, same order as Model.
    size_t item_idx = 0;
    for (const ModelObject* obj : model.objects) {
        for (const ModelInstance* inst : obj->instances) {
            if (item_idx >= items.size())
                break;

            if (stop_check && stop_check()) {
                BOOST_LOG_TRIVIAL(warning) << "[SnuggleVoxelize] Aborted by stop_check at part " << item_idx;
                return result;
            }

            SnugglePartData part_data;
            part_data.name = obj->name;
            if (obj->instances.size() > 1)
                part_data.name += " #" + std::to_string(item_idx);
            part_data.initial_z_rotation = static_cast<float>(inst->get_rotation().z());

            BOOST_LOG_TRIVIAL(debug) << "[SnuggleVoxelize] Part " << item_idx
                << " '" << part_data.name << "'"
                << " z_rot=" << part_data.initial_z_rotation << " rad"
                << " offset=(" << inst->get_offset().x() << "," << inst->get_offset().y() << ") mm"
                << " (offset NOT applied — nester places XY)";

            // Build the same transform that get_arrange_polygon uses
            Transform3d hull_trafo = build_hull_transform(*inst);

            BOOST_LOG_TRIVIAL(debug) << "[SnuggleVoxelize] Hull transform for '"
                << part_data.name << "':"
                << " rot_xy=(" << inst->get_rotation().x() << "," << inst->get_rotation().y() << ")"
                << " rot_z=0 (zeroed)"
                << " offset_z=" << inst->get_offset().z();

            // Collect and transform triangles
            std::vector<snuggle::Triangle> triangles = collect_triangles(*obj, hull_trafo);

            if (triangles.empty()) {
                BOOST_LOG_TRIVIAL(warning) << "[SnuggleVoxelize] Part '" << part_data.name
                    << "' has no triangles, skipping voxelization";
                // Still push an empty entry to maintain index correspondence
                part_data.info.name = part_data.name;
                part_data.info.initial_zrot = part_data.initial_z_rotation;
                result.push_back(std::move(part_data));
                ++item_idx;
                continue;
            }

            BOOST_LOG_TRIVIAL(debug) << "[SnuggleVoxelize] Part '" << part_data.name
                << "': " << triangles.size() << " triangles collected";

            // Compute mesh height and approximate hull area for nester metadata
            snuggle::AABB bounds;
            for (const auto& tri : triangles) {
                bounds.expand(tri.v0);
                bounds.expand(tri.v1);
                bounds.expand(tri.v2);
            }

            part_data.info.max_height_mm = bounds.max.z - bounds.min.z;
            // Approximate 2D footprint area from XY bounding box
            float dx = bounds.max.x - bounds.min.x;
            float dy = bounds.max.y - bounds.min.y;
            part_data.info.hull_area_mm2 = dx * dy;
            part_data.info.name = part_data.name;
            part_data.info.initial_zrot = part_data.initial_z_rotation;

            BOOST_LOG_TRIVIAL(debug) << "[SnuggleVoxelize] Part '" << part_data.name
                << "' bounds: [" << bounds.min.x << "," << bounds.min.y << "," << bounds.min.z
                << "] -> [" << bounds.max.x << "," << bounds.max.y << "," << bounds.max.z << "]"
                << " height=" << part_data.info.max_height_mm << " mm"
                << " area~=" << part_data.info.hull_area_mm2 << " mm^2";

            // Voxelize
            snuggle::VoxelGrid grid;
            snuggle::VoxError err = snuggle::voxelize_mesh(
                triangles.data(), triangles.size(),
                voxel_size_mm, grid,
                [&](float progress, const char* phase) -> bool {
                    // Check for abort
                    if (stop_check && stop_check())
                        return false;
                    return true;
                });

            if (err != snuggle::VoxError::OK) {
                BOOST_LOG_TRIVIAL(warning) << "[SnuggleVoxelize] Voxelization failed for '"
                    << part_data.name << "': " << snuggle::vox_error_str(err);
                // Push empty entry to maintain ordering
                result.push_back(std::move(part_data));
                ++item_idx;
                continue;
            }

            part_data.info.grid = std::move(grid);

            BOOST_LOG_TRIVIAL(debug) << "[SnuggleVoxelize] Part '" << part_data.name
                << "' voxelized: grid=" << part_data.info.grid.nx
                << "x" << part_data.info.grid.ny << "x" << part_data.info.grid.nz
                << " voxels, " << part_data.info.grid.count_solid() << " solid"
                << " (" << part_data.info.grid.memory_bytes() << " bytes)";

            result.push_back(std::move(part_data));
            ++item_idx;
        }
    }

    BOOST_LOG_TRIVIAL(info) << "[SnuggleVoxelize] Voxelization complete: "
        << result.size() << " parts processed";

    return result;
}

}} // namespace Slic3r::arrangement
