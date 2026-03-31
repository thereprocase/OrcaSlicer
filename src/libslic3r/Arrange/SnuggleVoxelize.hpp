// SnuggleVoxelize.hpp — Voxelize OrcaSlicer model parts for Snuggle nesting
//
// Extracts meshes from the Model, applies the SAME transform that
// get_arrange_polygon uses (XY rotation but NOT XY offset, NOT Z offset),
// then feeds to polite_voxelizer. This ensures the voxel grid is in the
// same coordinate frame as the convex hull used by ArrangePolygon.
//
// The XY offset is NOT applied because that's the nester's job — it places
// parts on the bed. The Z rotation is zeroed because get_arrange_polygon
// zeros it and stores it separately in ArrangePolygon::rotation.

#pragma once

#include "snuggle_nester.hpp"       // snuggle::PartInfo
#include "libslic3r/Model.hpp"
#include "libslic3r/Arrange.hpp"    // ArrangePolygons
#include <functional>
#include <string>
#include <vector>

namespace Slic3r { namespace arrangement {

/// Per-part data ready for the Snuggle nester.
/// Matches 1:1 with the ArrangePolygons vector ordering.
struct SnugglePartData {
    snuggle::PartInfo info;         ///< Voxel grid + metadata for the nester
    std::string       name;         ///< Human-readable part name for logging
    float             initial_z_rotation; ///< radians, from ModelInstance (passed to nester for lock_rotation mode)
};

/// Voxelize all arrangeable parts for Snuggle.
///
/// @param model        The OrcaSlicer Model containing objects/instances.
/// @param items        ArrangePolygons — used to match ordering and to know
///                     which instances to process (size determines part count).
/// @param voxel_size_mm Voxel resolution in mm (e.g. 1.0 for coarse, 0.5 for fine).
/// @param stop_check   Optional callback; returns true to abort voxelization.
/// @return             Vector of SnugglePartData, one per item, in same order.
///
/// COORDINATE FRAME: The mesh vertices are transformed by the same matrix
/// that get_arrange_polygon uses:
///   - Translation: Z offset only (XY zeroed — nester places XY).
///   - Rotation:    XY rotation applied, Z rotation zeroed (stored separately).
///   - Scale:       Full instance scale applied.
/// The resulting vertices are in MM (not scaled coords), suitable for Snuggle.
std::vector<SnugglePartData> voxelize_for_snuggle(
    const Model& model,
    const ArrangePolygons& items,
    float voxel_size_mm,
    std::function<bool()> stop_check = nullptr);

}} // namespace Slic3r::arrangement
