// SnuggleArrange.hpp — Bridge between OrcaSlicer's arrange API and Snuggle
#pragma once

#include "libslic3r/Arrange.hpp"
#include "libslic3r/Model.hpp"

namespace Slic3r { namespace arrangement {

// Run Snuggle 3D-aware arrangement on the given items.
// model_objects: the ModelObjects to get 3D meshes from (indexed by item order)
// items: ArrangePolygons to populate with results (translation, rotation, bed_idx)
// bed: bed boundary points
// params: arrangement parameters (uses use_snuggle, snuggle_lock_rotation, min_obj_distance)
void snuggle_arrange(
    ArrangePolygons&          items,
    const ArrangePolygons&    excludes,
    const Points&             bed,
    const ArrangeParams&      params,
    const Model&              model,
    std::string*              out_gpu_renderer = nullptr,
    float*                    out_gpu_ratio = nullptr);

}} // namespace Slic3r::arrangement
