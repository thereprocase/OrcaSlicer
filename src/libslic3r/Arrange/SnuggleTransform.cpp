// SnuggleTransform.cpp — Coordinate transform helpers for Snuggle <-> OrcaSlicer

#include "SnuggleTransform.hpp"
#include "libslic3r/libslic3r.h"   // SCALING_FACTOR, scale_(), unscale_()
#include <boost/log/trivial.hpp>
#include <cassert>
#include <cmath>
#include <limits>

namespace Slic3r { namespace arrangement { namespace snuggle_xform {

// ── Orca scaled coords -> Snuggle mm ────────────────────────────────
// Input space:  Orca scaled integer (coord_t). 1 unit = SCALING_FACTOR mm.
// Output space: Snuggle mm (float).
float scaled_to_mm(coord_t val)
{
    // unscale: val * SCALING_FACTOR -> mm
    double mm = unscale<double>(val);
    float result = static_cast<float>(mm);

    BOOST_LOG_TRIVIAL(trace) << "[SnuggleXform] scaled_to_mm: "
        << val << " scaled -> " << result << " mm"
        << " (SCALING_FACTOR=" << SCALING_FACTOR << ")";

    // Sanity: result should not be astronomical
    assert(std::isfinite(result));
    assert(std::abs(result) < 100000.0f); // 100 meters is absurd for a print bed

    return result;
}

// ── Snuggle mm -> Orca scaled coords ────────────────────────────────
// Input space:  Snuggle mm (double for precision on large beds).
// Output space: Orca scaled integer (coord_t). 1 unit = SCALING_FACTOR mm.
coord_t mm_to_scaled(double val_mm)
{
    // scale: val_mm / SCALING_FACTOR -> scaled integer
    coord_t result = scaled<coord_t>(val_mm);

    BOOST_LOG_TRIVIAL(trace) << "[SnuggleXform] mm_to_scaled: "
        << val_mm << " mm -> " << result << " scaled"
        << " (SCALING_FACTOR=" << SCALING_FACTOR << ")";

    return result;
}

// ── Bed polygon -> dimensions in mm ─────────────────────────────────
// Input space:  bed points in Orca scaled coords (coord_t).
// Output space: width/height/origin in mm (float).
void bed_dimensions_mm(const Points& bed,
                       float& width_mm, float& height_mm,
                       float& origin_x_mm, float& origin_y_mm)
{
    if (bed.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "[SnuggleXform] bed_dimensions_mm: empty bed polygon!";
        width_mm = height_mm = origin_x_mm = origin_y_mm = 0.0f;
        return;
    }

    // Compute bounding box in scaled coords
    BoundingBox bb(bed);

    // Convert each corner to mm
    origin_x_mm = scaled_to_mm(bb.min.x());
    origin_y_mm = scaled_to_mm(bb.min.y());
    float max_x_mm = scaled_to_mm(bb.max.x());
    float max_y_mm = scaled_to_mm(bb.max.y());

    width_mm  = max_x_mm - origin_x_mm;
    height_mm = max_y_mm - origin_y_mm;

    BOOST_LOG_TRIVIAL(debug) << "[SnuggleXform] bed_dimensions_mm: "
        << "scaled bb=[" << bb.min.x() << "," << bb.min.y()
        << "]->[" << bb.max.x() << "," << bb.max.y() << "]"
        << " -> mm origin=(" << origin_x_mm << "," << origin_y_mm << ")"
        << " size=(" << width_mm << " x " << height_mm << ") mm";

    // Sanity
    assert(width_mm >= 0.0f);
    assert(height_mm >= 0.0f);
}

// ── Snuggle placement -> ArrangePolygon fields ──────────────────────
// Input space:  Snuggle mm coords (placement relative to bed origin).
// Output space: Orca scaled coords (absolute) + radians.
void placement_to_arrange(float snuggle_x_mm, float snuggle_y_mm, float snuggle_rot_rad,
                          float bed_origin_x_mm, float bed_origin_y_mm,
                          Vec2crd& out_translation, double& out_rotation)
{
    // Convert Snuggle's bed-relative mm position to absolute mm,
    // then to Orca scaled coords.
    // Use double for intermediate sum to preserve precision on large beds
    double abs_x_mm = (double)snuggle_x_mm + (double)bed_origin_x_mm;
    double abs_y_mm = (double)snuggle_y_mm + (double)bed_origin_y_mm;

    out_translation.x() = mm_to_scaled(abs_x_mm);
    out_translation.y() = mm_to_scaled(abs_y_mm);

    // Rotation is already in radians — no conversion needed.
    // Snuggle radians == Orca radians.
    out_rotation = static_cast<double>(snuggle_rot_rad);

    BOOST_LOG_TRIVIAL(trace) << "[SnuggleXform] placement_to_arrange: "
        << "snuggle=(" << snuggle_x_mm << "," << snuggle_y_mm << ") mm"
        << " + bed_origin=(" << bed_origin_x_mm << "," << bed_origin_y_mm << ") mm"
        << " -> abs_mm=(" << abs_x_mm << "," << abs_y_mm << ")"
        << " -> scaled=(" << out_translation.x() << "," << out_translation.y() << ")"
        << " rot=" << out_rotation << " rad";
}

}}} // namespace Slic3r::arrangement::snuggle_xform
