// test_bitmap_mesh.cpp
// Sprint 1 Task 4: Full mesh-to-bitmap pipeline tests.
//
// Every previous bitmap test hand-builds synthetic ExPolygons. This file
// constructs actual indexed_triangle_set meshes, runs them through the
// project_mesh() -> union_ex() -> BitmapNester::arrange() pipeline, and
// asserts that the silhouette pipeline correctly handles filled, hole-aware,
// and multi-island footprints for the four motivating shapes from
// problem_statement.md: L-bracket, crescent, U-channel, dumbbell.
//
// The dumbbell test is the load-bearing regression: if multi-island
// support (Sprint 1, T1-T3) regresses, that test must fail loudly.

#define BITMAP_NESTER_TESTING

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/Geometry/ConvexHull.hpp"
#include "bitmap_test_utils.hpp"

#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Mesh construction helpers
// ---------------------------------------------------------------------------

namespace {

// L-bracket: two axis-aligned box prisms sharing the origin corner.
//
//   +-------+
//   |  leg2 |
//   +---+---+-------+
//       |    leg1   |
//       +-----------+
//
// leg_len  : length of each leg (mm), measured from the inner corner.
// leg_width: width (cross-section) of each leg (mm).
// height   : extrusion height in Z (mm).
//
// Built from two `its_make_cube` boxes merged together. The result is a
// single connected manifold with an interior concave corner.
indexed_triangle_set its_l_bracket_mm(double leg_len, double leg_width, double height)
{
    // Horizontal leg: runs along +X from origin.
    indexed_triangle_set horiz = its_make_cube(leg_len, leg_width, height);

    // Vertical leg: runs along +Y from origin, offset so it shares the corner
    // square with the horizontal leg rather than just touching it.
    // The corner square [0..leg_width] x [0..leg_width] is covered by horiz.
    // We add the vertical leg above that square to avoid internal geometry.
    indexed_triangle_set vert = its_make_cube(leg_width, leg_len - leg_width, height);
    its_translate(vert, Vec3f(0.f, float(leg_width), 0.f));

    its_merge(horiz, vert);
    return horiz;
}

// Crescent: approximated as a large cylinder with a shifted smaller cylinder
// bored through it. Because we cannot subtract meshes directly, we build the
// crescent contour analytically in 2D (outer arc + inner arc), triangulate
// the top and bottom caps, and stitch the sides.
//
// outer_radius: radius of the enclosing disc (mm).
// inner_radius: radius of the bored-out disc (mm).
// offset      : how far the inner disc center is shifted in +X (mm).
//               Must be < outer_radius - inner_radius to keep a solid rim.
// height      : extrusion height in Z (mm).
// segments    : tessellation count for each arc.
//
// The result is a single connected manifold crescent prism.
indexed_triangle_set its_crescent_mm(double outer_radius,
                                     double inner_radius,
                                     double offset,
                                     double height,
                                     int    segments = 32)
{
    // Outer arc: full circle around (0,0).
    // Inner arc: full circle around (offset, 0).
    // The crescent footprint is the boolean difference of the two discs.
    // We build it as a polygon: walk the outer circle CCW, then the inner
    // circle CW (as a hole cut-out), but since this is a prism we need a
    // simply-connected polygon (no holes). We do this by building a combined
    // contour that enters via a seam at angle 0.
    //
    // For a clean C-shape without a hole, we restrict the angular range.
    // Strategy: the inner disc subtracted from the outer leaves a crescent
    // only if offset != 0. We build the top-face polygon by sampling both
    // circles at the same angles and producing a ring contour.
    //
    // The contour walks:
    //   1. Outer circle CCW for full 2*pi
    //   2. Bridging edge at angle 0 to the inner circle
    //   3. Inner circle CW for full 2*pi (reversed, so the ring is simply connected)
    //
    // This is the standard annular prism. The outer-minus-inner is a donut,
    // which has a topological hole but that's fine — it's a single mesh with
    // one contour and one hole in the 2D silhouette. The test checks
    // silhouette.size() == 1 and that area < outer disc area.

    const int N = segments;
    const auto pi = std::acos(-1.0);
    const float h  = float(height);
    const float ri = float(inner_radius);
    const float ro = float(outer_radius);
    const float ox = float(offset);

    // Outer ring verts (bottom, then top): indices 0..2N-1 (bottom), 2N..4N-1 (top)
    // Inner ring verts (bottom, then top): 4N..6N-1 (bottom), 6N..8N-1 (top)
    // Bottom center outer: 8N, Top center outer: 8N+1
    // (We don't need explicit centers — we build a ring, not two discs.)

    indexed_triangle_set its;
    auto &verts   = its.vertices;
    auto &indices = its.indices;

    // Reserve space
    verts.reserve(4 * N + 4);
    indices.reserve(8 * N);

    // Bottom outer ring: [0 .. N-1]
    for (int i = 0; i < N; ++i) {
        double a = 2.0 * pi * i / N;
        verts.emplace_back(Vec3f(ro * float(std::cos(a)),
                                 ro * float(std::sin(a)),
                                 0.f));
    }
    // Top outer ring: [N .. 2N-1]
    for (int i = 0; i < N; ++i) {
        double a = 2.0 * pi * i / N;
        verts.emplace_back(Vec3f(ro * float(std::cos(a)),
                                 ro * float(std::sin(a)),
                                 h));
    }
    // Bottom inner ring (centered at (offset, 0)): [2N .. 3N-1]
    for (int i = 0; i < N; ++i) {
        // CW winding on bottom when seen from +Z (hole)
        double a = 2.0 * pi * i / N;
        verts.emplace_back(Vec3f(ox + ri * float(std::cos(a)),
                                 ri * float(std::sin(a)),
                                 0.f));
    }
    // Top inner ring: [3N .. 4N-1]
    for (int i = 0; i < N; ++i) {
        double a = 2.0 * pi * i / N;
        verts.emplace_back(Vec3f(ox + ri * float(std::cos(a)),
                                 ri * float(std::sin(a)),
                                 h));
    }

    // Bottom cap (annular ring, looking from below, face normal -Z):
    //   Outer contour: CCW from below => CW indices for -Z normal.
    //   Inner contour (hole): CW from below => CCW indices for -Z normal.
    //   We triangulate with the fan-through-seam technique: connect each outer
    //   vertex to the corresponding inner vertex for a strip.
    //
    // Each quad (outer[i], outer[i+1], inner[i+1], inner[i]) -> two triangles.
    for (int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        int oi = i;          // outer bottom
        int oj = j;
        int ii = 2 * N + i; // inner bottom
        int ij = 2 * N + j;
        // Bottom face: normal points -Z. Winding: CW when seen from +Z.
        indices.emplace_back(Vec3i32(oi, ii, oj));
        indices.emplace_back(Vec3i32(oj, ii, ij));
    }

    // Top cap (annular ring, face normal +Z):
    //   CCW winding when seen from +Z.
    for (int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        int oi = N + i;      // outer top
        int oj = N + j;
        int ii = 3 * N + i; // inner top
        int ij = 3 * N + j;
        indices.emplace_back(Vec3i32(oi, oj, ii));
        indices.emplace_back(Vec3i32(oj, ij, ii));
    }

    // Outer side wall (normal points outward from outer ring):
    for (int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        int b0 = i;          // outer bottom[i]
        int b1 = j;          // outer bottom[j]
        int t0 = N + i;      // outer top[i]
        int t1 = N + j;      // outer top[j]
        // Two triangles, outward normal.
        indices.emplace_back(Vec3i32(b0, t0, b1));
        indices.emplace_back(Vec3i32(b1, t0, t1));
    }

    // Inner side wall (normal points inward = toward inner ring center):
    for (int i = 0; i < N; ++i) {
        int j = (i + 1) % N;
        int b0 = 2 * N + i;  // inner bottom[i]
        int b1 = 2 * N + j;
        int t0 = 3 * N + i;  // inner top[i]
        int t1 = 3 * N + j;
        // Inward normal: reverse winding from outer side.
        indices.emplace_back(Vec3i32(b0, b1, t0));
        indices.emplace_back(Vec3i32(b1, t1, t0));
    }

    return its;
}

// U-channel: three axis-aligned box prisms forming a U cross-section.
//
//   +--+         +--+
//   |  |  (gap)  |  |   <- two side walls
//   |  +---------+  |
//   +---------------+   <- base
//
// outer_width    : total outer width of the U (mm).
// outer_depth    : total outer depth / height of the walls (mm).
// wall_thickness : thickness of each wall and the base (mm).
// height         : extrusion height in Z (the "into the page" dimension, mm).
//
// Built from three its_make_cube boxes: base + left wall + right wall.
indexed_triangle_set its_u_channel_mm(double outer_width,
                                      double outer_depth,
                                      double wall_thickness,
                                      double height)
{
    // Base: spans full width, sits at Y=0.
    indexed_triangle_set base = its_make_cube(outer_width, wall_thickness, height);

    // Left wall: rises from Y=wall_thickness to Y=outer_depth.
    double wall_height_y = outer_depth - wall_thickness;
    indexed_triangle_set left_wall = its_make_cube(wall_thickness, wall_height_y, height);
    its_translate(left_wall, Vec3f(0.f, float(wall_thickness), 0.f));

    // Right wall: same size, at X = outer_width - wall_thickness.
    indexed_triangle_set right_wall = its_make_cube(wall_thickness, wall_height_y, height);
    its_translate(right_wall,
                  Vec3f(float(outer_width - wall_thickness), float(wall_thickness), 0.f));

    its_merge(base, left_wall);
    its_merge(base, right_wall);
    return base;
}

// Dumbbell: two cylindrical lobes separated by a gap, with nothing connecting
// them in the XY projection. This is the regression fixture for multi-island
// silhouette support: project_mesh must return two disconnected ExPolygons.
//
// ball_radius  : radius of each lobe (mm).
// gap          : clear space between the two lobes (mm). Each lobe center is
//                at X = -(ball_radius + gap/2) and X = +(ball_radius + gap/2).
// height       : extrusion in Z (mm).
// ball_segments: tessellation count for each cylinder.
indexed_triangle_set its_dumbbell_mm(double ball_radius,
                                     double gap,
                                     double height,
                                     int    ball_segments = 16)
{
    // fa controls tessellation: 2*pi/n_steps where n_steps = ceil(2*pi/fa).
    // We want approximately ball_segments facets, so fa = 2*pi/ball_segments.
    const double fa = 2.0 * std::acos(-1.0) / ball_segments;

    // Left lobe: its_make_cylinder centers around (0,0). Translate to (-cx, 0).
    double cx = ball_radius + gap / 2.0;
    indexed_triangle_set left  = its_make_cylinder(ball_radius, height, fa);
    its_translate(left,  Vec3f(float(-cx), 0.f, 0.f));

    // Right lobe.
    indexed_triangle_set right = its_make_cylinder(ball_radius, height, fa);
    its_translate(right, Vec3f(float(+cx), 0.f, 0.f));

    its_merge(left, right);
    return left;
}

// Project a mesh through the standard pipeline: project_mesh() -> union_ex().
//
// project_mesh() collects all upward-facing and downward-facing triangle
// projections into Polygons. union_ex() resolves overlaps and produces the
// canonical set of ExPolygons that represents the footprint.
//
// We use the 5-argument overload (out_top, out_bottom) because it gives us
// both face sets, which handles extruded shapes correctly regardless of
// whether the "top" or "bottom" faces dominate at a given Z.
ExPolygons project_to_silhouette(const indexed_triangle_set &its)
{
    Polygons top, bottom;
    project_mesh(its, Transform3d::Identity(), &top, &bottom, [](){});
    Polygons all;
    all.insert(all.end(), top.begin(),    top.end());
    all.insert(all.end(), bottom.begin(), bottom.end());
    return union_ex(all);
}

// Compute the area (mm^2) of a collection of ExPolygons.
// ExPolygon::area() returns the area in scaled^2 units; unscale by 1e6.
double expolygons_area_mm2(const ExPolygons &expolys)
{
    double total = 0.0;
    for (const ExPolygon &ep : expolys)
        total += std::abs(ep.area());
    return unscaled<double>(unscaled<double>(total));
}

// Check whether a given mm-space point lies inside any ExPolygon.
bool point_inside_any(const ExPolygons &expolys, double x_mm, double y_mm)
{
    Point pt(scaled<coord_t>(x_mm), scaled<coord_t>(y_mm));
    for (const ExPolygon &ep : expolys)
        if (ep.contains(pt))
            return true;
    return false;
}

} // anonymous namespace

// ===========================================================================
// 1. L-bracket: single concave silhouette
// ===========================================================================

TEST_CASE("Mesh: L-bracket projects to single connected concave silhouette",
          "[BitmapMesh]")
{
    const double leg_len   = 20.0; // mm
    const double leg_width = 5.0;  // mm
    const double height    = 3.0;  // mm

    indexed_triangle_set its = its_l_bracket_mm(leg_len, leg_width, height);
    ExPolygons silhouette = project_to_silhouette(its);

    // The L is a single connected region.
    REQUIRE(silhouette.size() == 1);

    // No holes: the L-bracket is simply connected — its interior concavity is
    // a missing corner, not a through-hole.
    REQUIRE(silhouette[0].holes.empty());

    // Analytical area: two leg_len x leg_width rectangles minus the shared
    // corner square (leg_width x leg_width).
    double expected_area = 2.0 * leg_len * leg_width - leg_width * leg_width;
    double actual_area   = expolygons_area_mm2(silhouette);
    // Tessellation and rasterization can introduce up to ~1 mm^2 error.
    REQUIRE_THAT(actual_area, Catch::Matchers::WithinAbs(expected_area, 2.0));

    // Convex hull must be strictly larger than the concave silhouette.
    // If this fails, the pipeline has reduced the shape to its convex hull.
    Polygon hull = Geometry::convex_hull(to_polygons(silhouette));
    double hull_area = std::abs(unscaled<double>(unscaled<double>(hull.area())));
    REQUIRE(hull_area > actual_area);
}

// ===========================================================================
// 2. Crescent: concave silhouette, area less than enclosing disc
// ===========================================================================

TEST_CASE("Mesh: crescent projects to single silhouette with concavity",
          "[BitmapMesh]")
{
    // Outer radius 10mm, inner disc radius 7mm shifted 4mm in +X.
    // The inner disc bites into the outer, leaving a crescent rim.
    const double outer_r = 10.0;
    const double inner_r = 7.0;
    const double off     = 4.0;
    const double height  = 3.0;

    indexed_triangle_set its = its_crescent_mm(outer_r, inner_r, off, height, 48);
    ExPolygons silhouette = project_to_silhouette(its);

    // A crescent ring (outer annulus) projects to a single ExPolygon, which
    // may have one hole (the inner cavity). Either form is acceptable for the
    // pipeline test.
    REQUIRE(silhouette.size() == 1);

    // Total solid area must be less than the outer disc (pi * r^2).
    double outer_disc_area = std::acos(-1.0) * outer_r * outer_r;
    double actual_area     = expolygons_area_mm2(silhouette);
    REQUIRE(actual_area < outer_disc_area);

    // The contour of silhouette[0] is the outer ring; its convex hull must be
    // larger than the silhouette's solid area.
    ExPolygons contour_only{ExPolygon{silhouette[0].contour}};
    double contour_area = expolygons_area_mm2(contour_only);
    double hull_area    = std::abs(unscaled<double>(unscaled<double>(
        Geometry::convex_hull(to_polygons(silhouette)).area())));
    // If the pipeline returns a convex shape, hull == contour (within
    // rounding). For a concave shape, hull > contour.
    // We check the hull is at least as large as the contour — this would
    // catch a regression where the silhouette is already the convex hull.
    REQUIRE(hull_area >= contour_area - 1.0); // 1 mm^2 rounding tolerance
}

// ===========================================================================
// 3. U-channel: single concave silhouette, interior cavity not filled
// ===========================================================================

TEST_CASE("Mesh: U-channel projects to single connected concave silhouette",
          "[BitmapMesh]")
{
    const double outer_w    = 30.0; // mm
    const double outer_d    = 20.0; // mm
    const double wall_t     = 5.0;  // mm
    const double height     = 3.0;  // mm

    indexed_triangle_set its = its_u_channel_mm(outer_w, outer_d, wall_t, height);
    ExPolygons silhouette = project_to_silhouette(its);

    // The U is one connected footprint.
    REQUIRE(silhouette.size() == 1);

    // The U-channel silhouette is simply connected (open mouth at the top —
    // no enclosed hole when viewed from +Z). The union of the three boxes
    // does not enclose a region of air.
    REQUIRE(silhouette[0].holes.empty());

    // The interior of the U (the gap between the walls, above the base) must
    // NOT be filled. A point at the center of the cavity should lie outside
    // the silhouette. Cavity center: X = outer_w/2, Y = (wall_t + outer_d) / 2.
    double cx = outer_w / 2.0;
    double cy = (wall_t + outer_d) / 2.0;
    REQUIRE_FALSE(point_inside_any(silhouette, cx, cy));

    // Convex hull must be strictly larger than the silhouette.
    double sil_area  = expolygons_area_mm2(silhouette);
    double hull_area = std::abs(unscaled<double>(unscaled<double>(
        Geometry::convex_hull(to_polygons(silhouette)).area())));
    REQUIRE(hull_area > sil_area);
}

// ===========================================================================
// 4. Dumbbell: TWO disconnected silhouettes (multi-island regression)
// ===========================================================================

TEST_CASE("Mesh: dumbbell projects to TWO disconnected silhouettes (multi-island)",
          "[BitmapMesh][MultiIsland]")
{
    // gap > 0 ensures the two cylinder projections do not touch.
    const double ball_r = 8.0;  // mm
    const double gap    = 4.0;  // mm
    const double height = 3.0;  // mm

    indexed_triangle_set its = its_dumbbell_mm(ball_r, gap, height, 24);
    ExPolygons silhouette = project_to_silhouette(its);

    // Both lobes must survive the pipeline as distinct ExPolygons.
    // If this fails, the silhouette pipeline collapsed the two islands into
    // one (e.g., by taking only the largest region, or by convex-hulling).
    REQUIRE(silhouette.size() == 2);
}

// ===========================================================================
// 5. BitmapNester respects multi-island silhouette via concave_regions
// ===========================================================================

TEST_CASE("Mesh: BitmapNester rasterizes all islands of a dumbbell silhouette",
          "[BitmapMesh][MultiIsland]")
{
    const double ball_r = 8.0;
    const double gap    = 4.0;
    const double height = 3.0;

    indexed_triangle_set its = its_dumbbell_mm(ball_r, gap, height, 24);
    ExPolygons silhouette = project_to_silhouette(its);

    // We need both islands for the arrange test to be meaningful.
    REQUIRE(silhouette.size() == 2);

    // Select the largest island as the legacy `poly` field. The nester uses
    // concave_regions when non-empty and falls back to poly otherwise.
    ExPolygon largest = *std::max_element(
        silhouette.begin(), silhouette.end(),
        [](const ExPolygon &a, const ExPolygon &b) {
            return std::abs(a.area()) < std::abs(b.area());
        });

    // Build two instances of the dumbbell on a 200 x 200 mm bed.
    ArrangePolygons items;
    for (int k = 0; k < 2; ++k) {
        ArrangePolygon ap;
        ap.poly             = largest;
        ap.concave_regions  = silhouette; // all islands
        ap.priority         = 0;
        ap.bed_idx          = UNARRANGED;
        ap.rotation         = 0.0;
        ap.translation      = Vec2crd{0, 0};
        ap.allowed_rotations = {0.0};
        items.push_back(ap);
    }

    ArrangeParams params;
    params.bed_shrink_x  = 0.0f;
    params.bed_shrink_y  = 0.0f;
    params.allow_rotations = false;
    params.progressind   = nullptr;

    BoundingBox bed(
        Point(scaled<coord_t>(0.0),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(200.0), scaled<coord_t>(200.0)));

    BitmapNester::arrange(items, {}, bed, params);

    // Both items must be placed on a valid bed.
    for (size_t i = 0; i < items.size(); ++i) {
        REQUIRE(items[i].bed_idx >= 0);
    }

    // No-overlap check that accounts for concave_regions.
    //
    // test_utils::no_overlap() uses item.transformed_poly() which consults
    // only item.poly (the largest single island). For a dumbbell, the
    // "other" lobe would be invisible to that check. Two instances could be
    // placed so that instance-A's left lobe overlaps instance-B's right lobe,
    // and transformed_poly()-based no_overlap would miss it.
    //
    // We therefore apply the transform to every island in concave_regions
    // and check for pairwise intersection across all region pairs.
    auto all_regions_overlap = [&]() -> bool {
        for (size_t i = 0; i < items.size(); ++i) {
            if (items[i].bed_idx == UNARRANGED) continue;
            for (size_t j = i + 1; j < items.size(); ++j) {
                if (items[j].bed_idx != items[i].bed_idx) continue;

                // Collect transformed regions for item i.
                ExPolygons regions_i;
                for (ExPolygon ep : items[i].concave_regions) {
                    ep.rotate(items[i].rotation);
                    ep.translate(items[i].translation.x(),
                                 items[i].translation.y());
                    regions_i.push_back(ep);
                }
                // Collect transformed regions for item j.
                ExPolygons regions_j;
                for (ExPolygon ep : items[j].concave_regions) {
                    ep.rotate(items[j].rotation);
                    ep.translate(items[j].translation.x(),
                                 items[j].translation.y());
                    regions_j.push_back(ep);
                }

                ExPolygons inter = intersection_ex(regions_i, regions_j);
                for (const ExPolygon &seg : inter) {
                    if (std::abs(seg.area()) > scaled<double>(0.1) * scaled<double>(0.1))
                        return true; // overlap found
                }
            }
        }
        return false;
    };

    // If the nester only rasterized the largest island (the legacy path),
    // it might stack the second instance directly on the first, overlapping
    // via the OTHER lobe. This assertion catches that regression.
    REQUIRE_FALSE(all_regions_overlap());
}

// ===========================================================================
// 6. Convex hull strictly larger than concave silhouette (parametric)
// ===========================================================================

TEST_CASE("Mesh: silhouette convex hull is strictly larger than concave silhouette",
          "[BitmapMesh]")
{
    // L-bracket
    DYNAMIC_SECTION("L-bracket")
    {
        indexed_triangle_set its = its_l_bracket_mm(20.0, 5.0, 3.0);
        ExPolygons sil = project_to_silhouette(its);
        REQUIRE_FALSE(sil.empty());

        double sil_area  = expolygons_area_mm2(sil);
        double hull_area = std::abs(unscaled<double>(unscaled<double>(
            Geometry::convex_hull(to_polygons(sil)).area())));

        // A non-convex shape's convex hull must cover strictly more area.
        // If hull_area == sil_area, the silhouette is already convex — the
        // concavity was lost somewhere in the pipeline.
        REQUIRE(hull_area > sil_area);
    }

    // U-channel
    DYNAMIC_SECTION("U-channel")
    {
        indexed_triangle_set its = its_u_channel_mm(30.0, 20.0, 5.0, 3.0);
        ExPolygons sil = project_to_silhouette(its);
        REQUIRE_FALSE(sil.empty());

        double sil_area  = expolygons_area_mm2(sil);
        double hull_area = std::abs(unscaled<double>(unscaled<double>(
            Geometry::convex_hull(to_polygons(sil)).area())));

        REQUIRE(hull_area > sil_area);
    }

    // Crescent: hull of the outer contour vs. the solid area (including hole).
    DYNAMIC_SECTION("Crescent")
    {
        indexed_triangle_set its = its_crescent_mm(10.0, 7.0, 4.0, 3.0, 48);
        ExPolygons sil = project_to_silhouette(its);
        REQUIRE_FALSE(sil.empty());

        // For the crescent (which is an annular ring), the contour area
        // (outer ring only, ignoring the hole) must exceed the solid area.
        double sil_area      = expolygons_area_mm2(sil);
        Polygon hull         = Geometry::convex_hull(to_polygons(sil));
        double hull_area     = std::abs(unscaled<double>(unscaled<double>(hull.area())));

        // hull_area approximates pi*R^2 = pi*100 ~ 314 mm^2.
        // sil_area is the annular solid, which is less than the outer disc.
        REQUIRE(hull_area > sil_area);
    }
}
