// SnuggleTransform.hpp — Coordinate transform helpers for Snuggle <-> OrcaSlicer
//
// OrcaSlicer internal coordinate spaces:
//   - "MM coords":     floating point millimeters (ModelInstance offset, user-facing)
//   - "Scaled coords": integer coord_t = mm / SCALING_FACTOR (ArrangePolygon, internal geometry)
//
// Snuggle operates entirely in MM (float).
//
// This header provides explicit, logged, testable conversions so that
// no inline scaling arithmetic appears anywhere in the integration code.

#pragma once

#include "libslic3r/Point.hpp"      // coord_t, Vec2crd, Points, BoundingBox
#include "libslic3r/BoundingBox.hpp"

namespace Slic3r { namespace arrangement { namespace snuggle_xform {

// ── Orca scaled coords (coord_t) -> Snuggle mm (float) ──────────────
// Input:  val in Orca scaled integer coordinates (coord_t).
// Output: equivalent value in millimeters (float).
// Formula: mm = val * SCALING_FACTOR
float scaled_to_mm(coord_t val);

// ── Snuggle mm (float) -> Orca scaled coords (coord_t) ──────────────
// Input:  val_mm in millimeters (float), from Snuggle placement result.
// Output: equivalent value in Orca scaled integer coordinates (coord_t).
// Formula: scaled = mm / SCALING_FACTOR
coord_t mm_to_scaled(double val_mm);

// ── Bed points (scaled) -> bed dimensions in mm ─────────────────────
// Input:  bed — polygon vertices in Orca scaled coordinates (coord_t).
//         These come from the arrange infrastructure (already scaled).
// Output: width_mm, height_mm — bed size in millimeters.
//         origin_x_mm, origin_y_mm — min corner of bed bounding box in mm.
//
// Computes bounding box of the bed polygon, converts to mm.
void bed_dimensions_mm(const Points& bed,
                       float& width_mm, float& height_mm,
                       float& origin_x_mm, float& origin_y_mm);

// ── Snuggle placement result -> ArrangePolygon fields ───────────────
// Input:  snuggle_x_mm, snuggle_y_mm — placement center in Snuggle mm coords
//           (origin at bed min corner, i.e. relative to bed_origin).
//         snuggle_rot_rad — Z rotation in radians from Snuggle.
//         bed_origin_x_mm, bed_origin_y_mm — bed min corner in mm
//           (from bed_dimensions_mm; needed because Snuggle's coordinate
//            system has origin at the bed's min corner, while Orca's
//            ArrangePolygon translation is in absolute scaled coords).
// Output: out_translation — Vec2crd in Orca scaled coordinates (absolute).
//         out_rotation — Z rotation in radians (pass-through, no scaling).
//
// The translation is: scaled(snuggle_pos_mm + bed_origin_mm)
// The rotation passes through unchanged (radians are radians).
void placement_to_arrange(float snuggle_x_mm, float snuggle_y_mm, float snuggle_rot_rad,
                          float bed_origin_x_mm, float bed_origin_y_mm,
                          Vec2crd& out_translation, double& out_rotation);

}}} // namespace Slic3r::arrangement::snuggle_xform
