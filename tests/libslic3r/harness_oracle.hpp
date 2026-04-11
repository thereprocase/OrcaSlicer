#pragma once

// harness_oracle.hpp
//
// Provably-correct geometry primitives used as ground-truth oracles
// for validating faster methods (bitmap AND, Boost.Geometry
// intersection, FFT correlation feasibility). See
// docs/PR_C2_TEST_HARNESS_VALIDATION.md for the full strategy.
//
// Rules:
//   - These helpers must never appear in production code paths. They
//     are test-only fixtures for validating the fast methods.
//   - Never optimize. Correctness over speed.
//   - Every fast method we trust must be validated against this
//     oracle on analytical test cases before we believe any quality
//     measurement that depends on it.
//
// Algorithm: winding number (Sunday's algorithm). For each edge of
// the polygon, count upward crossings of a horizontal ray cast from
// the query point. The winding number is the net count. Nonzero =
// inside; zero = outside. Works correctly on simple polygons of any
// winding direction and on self-intersecting polygons (which give
// fractional coverage that odd-even rules mis-handle).

#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"

#include <cmath>
#include <cstddef>

namespace Slic3r { namespace arrangement { namespace test_utils {

// ---------------------------------------------------------------------------
// Boundary policy
// ---------------------------------------------------------------------------
//
// A query point that lies exactly on a polygon edge is ambiguous: it
// belongs to neither the interior nor the exterior strictly. Our
// policy, chosen once and enforced consistently:
//
//   - Points strictly inside (winding != 0) are interior.
//   - Points on an edge or vertex are NOT interior (wn == 0 by the
//     upward-crossing discriminator).
//   - This matches the behavior of the downstream bitmap collision
//     check, which treats touching pixels as non-overlapping, and
//     matches Orca's convention that touching parts are legal.
//
// If the caller needs the opposite policy (touching == inside), they
// can explicitly test for on-boundary via a separate point-to-edge
// distance check. The oracle does not provide that because the
// bitmap pipeline never needs it.

// Returns true if (x_mm, y_mm) lies exactly on the segment (a, b).
// "Exactly" here means the signed area of the triangle (a, b, p) is
// within a tolerance of zero AND the point's projection onto the
// segment lies within [0, 1]. Used by the on-boundary short-circuit
// in point_in_polygon_oracle to implement the boundary policy
// ("on edge/vertex = exterior") — Sunday's winding algorithm alone
// does not handle boundary points consistently.
inline bool point_on_segment(double x_mm, double y_mm,
                             const Point& a, const Point& b,
                             double eps_mm = 1e-9)
{
    double ax = unscaled<double>(a.x());
    double ay = unscaled<double>(a.y());
    double bx = unscaled<double>(b.x());
    double by = unscaled<double>(b.y());
    double dx = bx - ax;
    double dy = by - ay;
    double len2 = dx * dx + dy * dy;
    if (len2 <= eps_mm * eps_mm) {
        // Degenerate segment: treat as a point.
        return std::abs(x_mm - ax) <= eps_mm && std::abs(y_mm - ay) <= eps_mm;
    }
    // Signed area of triangle (a, b, p) scaled by 2.
    double cross = dx * (y_mm - ay) - dy * (x_mm - ax);
    if (std::abs(cross) > eps_mm * std::sqrt(len2))
        return false;  // not on the infinite line through a-b
    // Parameter t of the projection of p onto a-b. t in [0, 1] = on segment.
    double t = ((x_mm - ax) * dx + (y_mm - ay) * dy) / len2;
    return t >= -eps_mm && t <= 1.0 + eps_mm;
}

// Returns true if (x_mm, y_mm) lies on any edge of the polygon.
// Linear scan — slow but correct.
inline bool point_on_polygon_boundary(double x_mm, double y_mm,
                                      const Points& pts,
                                      double eps_mm = 1e-9)
{
    const std::size_t n = pts.size();
    if (n < 2) return false;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = pts[i];
        const Point& b = pts[(i + 1) % n];
        if (point_on_segment(x_mm, y_mm, a, b, eps_mm))
            return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// winding_number
// ---------------------------------------------------------------------------
//
// Sunday's algorithm adapted for Slic3r's scaled integer coordinate
// system. Query point is in mm (unscaled double precision). Polygon
// vertices are in scaled integer (coord_t). Returns the net winding
// number — positive for CCW containment, negative for CW containment,
// zero for exterior.
//
// NOTE: this function does NOT handle boundary points specially —
// Sunday's algorithm alone produces inconsistent results for points
// that fall exactly on edges or vertices. Callers that need the
// stated boundary policy ("on edge/vertex = exterior") should use
// point_in_polygon_oracle(), which short-circuits on the boundary
// check before running the winding-number algorithm.
//
// This is the slow oracle. Do not call from production code.
inline int winding_number(double x_mm, double y_mm, const Points& pts)
{
    const std::size_t n = pts.size();
    if (n < 3) return 0;
    int wn = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = pts[i];
        const Point& b = pts[(i + 1) % n];
        double ax = unscaled<double>(a.x());
        double ay = unscaled<double>(a.y());
        double bx = unscaled<double>(b.x());
        double by = unscaled<double>(b.y());

        // Upward-crossing rule. Edge goes from a to b.
        //   if ay <= y_mm AND by > y_mm: upward crossing — count if
        //   the point is strictly to the LEFT of the edge.
        //   if ay > y_mm AND by <= y_mm: downward crossing — count if
        //   the point is strictly to the RIGHT of the edge.
        //
        // "Strictly to the left" is computed as the sign of the cross
        // product (b - a) x (p - a). Positive = left; negative = right;
        // zero = on the line extending the edge. Points on the edge
        // itself produce a zero and are NOT counted — this implements
        // the boundary policy above.
        if (ay <= y_mm) {
            if (by > y_mm) {
                double cross = (bx - ax) * (y_mm - ay) - (x_mm - ax) * (by - ay);
                if (cross > 0) ++wn;
            }
        } else {
            if (by <= y_mm) {
                double cross = (bx - ax) * (y_mm - ay) - (x_mm - ax) * (by - ay);
                if (cross < 0) --wn;
            }
        }
    }
    return wn;
}

// Convenience wrapper: implements the stated boundary policy
// (on edge/vertex = exterior). Short-circuits on the boundary check
// before delegating to the winding-number algorithm. Any point
// strictly inside the polygon returns true; any point on the boundary
// or outside returns false.
inline bool point_in_polygon_oracle(double x_mm, double y_mm,
                                    const Polygon& poly)
{
    if (point_on_polygon_boundary(x_mm, y_mm, poly.points)) return false;
    return winding_number(x_mm, y_mm, poly.points) != 0;
}

// ExPolygon variant. A point is interior to an ExPolygon if it is
// strictly interior to the outer contour AND NOT strictly interior
// to any hole AND not on any boundary (outer or hole).
inline bool point_in_expolygon_oracle(double x_mm, double y_mm,
                                      const ExPolygon& ep)
{
    // On any boundary → exterior per boundary policy.
    if (point_on_polygon_boundary(x_mm, y_mm, ep.contour.points)) return false;
    for (const Polygon& hole : ep.holes)
        if (point_on_polygon_boundary(x_mm, y_mm, hole.points)) return false;

    if (winding_number(x_mm, y_mm, ep.contour.points) == 0)
        return false;
    for (const Polygon& hole : ep.holes) {
        if (winding_number(x_mm, y_mm, hole.points) != 0)
            return false;  // inside a hole = not interior
    }
    return true;
}

// ---------------------------------------------------------------------------
// expolygons_overlap_oracle
// ---------------------------------------------------------------------------
//
// Grid-sampled overlap test. Walks a regular grid over the bounding
// box of expolygon A and tests each sample point against B via the
// winding-number oracle. Returns true if ANY sample of A's interior
// falls inside B's interior.
//
// Inherently one-directional: this asks "is any point of A inside B?"
// For a complete symmetric check use expolygons_overlap_oracle_both
// which runs the sweep in both directions.
//
// Grid step controls the sampling density. Default 0.1 mm is fine
// for shapes larger than ~10 mm; drop to 0.05 or 0.02 for smaller.
// Runtime is O(bbox_area_mm² / grid_step²) per direction — use
// sparingly, this is a correctness oracle not a fast method.
inline bool expolygons_overlap_oracle(const ExPolygon& a, const ExPolygon& b,
                                      double grid_step_mm = 0.1)
{
    if (a.contour.points.empty() || b.contour.points.empty()) return false;
    BoundingBox bba = get_extents(a);
    double xmin = unscaled<double>(bba.min.x());
    double xmax = unscaled<double>(bba.max.x());
    double ymin = unscaled<double>(bba.min.y());
    double ymax = unscaled<double>(bba.max.y());
    if (grid_step_mm <= 0.0) grid_step_mm = 0.1;
    for (double y = ymin; y <= ymax; y += grid_step_mm) {
        for (double x = xmin; x <= xmax; x += grid_step_mm) {
            if (point_in_expolygon_oracle(x, y, a) &&
                point_in_expolygon_oracle(x, y, b))
                return true;
        }
    }
    return false;
}

// Symmetric variant: runs both A-into-B and B-into-A sweeps. For
// most realistic shapes A-into-B is sufficient (if A overlaps B then
// some sample of A is inside B), but for pathological cases where
// A's sample grid misses a sliver intersection, the reverse sweep
// catches it. Slower but more robust.
inline bool expolygons_overlap_oracle_symmetric(const ExPolygon& a,
                                                const ExPolygon& b,
                                                double grid_step_mm = 0.1)
{
    return expolygons_overlap_oracle(a, b, grid_step_mm) ||
           expolygons_overlap_oracle(b, a, grid_step_mm);
}

// ---------------------------------------------------------------------------
// Relationship to the C1 no_overlap helper
// ---------------------------------------------------------------------------
//
// The existing `test_utils::no_overlap` helper in bitmap_test_utils.hpp
// uses polygon intersection from Clipper/Boost.Geometry with a tiny
// area tolerance. It's fast but has known edge cases (sub-pixel
// overlap fragility, discussed in the 2026-04-11 session diary).
//
// The oracle in this file is the ground truth. When no_overlap and
// this oracle disagree on a test case, the oracle wins and no_overlap
// is broken. Conversely, when they agree across many cases, both are
// trustworthy for that class of input.
//
// Every algorithmic change that touches the C2 scoring path should
// re-run a cross-check suite (see test_harness_oracle.cpp) that
// compares no_overlap vs this oracle on a battery of analytical
// cases. Any divergence is a harness bug, not a solver bug.

}}} // namespace Slic3r::arrangement::test_utils
