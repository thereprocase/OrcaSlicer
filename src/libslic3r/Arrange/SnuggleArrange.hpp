// SnuggleArrange.hpp — Main entry point for Snuggle 3D nesting in OrcaSlicer
//
// Replaces the default 2D arranger with Snuggle's voxel-based genetic nester.
// Handles the full pipeline: voxelize -> nest -> convert results back.

#pragma once

#include "libslic3r/Arrange.hpp"    // ArrangePolygons, ArrangeParams, Points
#include "libslic3r/Model.hpp"

namespace Slic3r { namespace arrangement {

/// Run Snuggle 3D nesting on the given items.
///
/// @param items    ArrangePolygons to arrange. On return, each item's
///                 translation and rotation fields are updated with the
///                 nesting result (in Orca scaled coords + radians).
///                 bed_idx is set to 0 for successfully placed parts,
///                 or UNARRANGED for parts that could not be placed.
/// @param excludes Fixed items on the bed (e.g. wipe tower). Their positions
///                 are respected but not changed. (Currently logged but not
///                 used as obstacles in Snuggle — future work.)
/// @param bed      Bed polygon vertices in Orca scaled coordinates.
/// @param params   Arrange parameters (min_obj_distance, allow_rotations, etc.).
/// @param model    The full Model, needed to extract meshes for voxelization.
///
/// COORDINATE FLOW:
///   1. bed (scaled) -> bed_dimensions_mm via SnuggleTransform
///   2. Model meshes -> voxel grids via SnuggleVoxelize (mm coords)
///   3. Snuggle nester runs entirely in mm
///   4. Results (mm) -> ArrangePolygon fields (scaled) via SnuggleTransform
void snuggle_arrange(
    ArrangePolygons& items,
    const ArrangePolygons& excludes,
    const Points& bed,
    const ArrangeParams& params,
    const Model& model);

}} // namespace Slic3r::arrangement
