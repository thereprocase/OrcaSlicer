// test_harness_oracle.cpp
//
// Phase 1 of the C2 harness validation plan
// (docs/PR_C2_TEST_HARNESS_VALIDATION.md). Exhaustively tests the
// winding-number point-in-polygon oracle against analytically known
// cases, then cross-checks the oracle vs the existing no_overlap
// helper on known-overlapping and known-non-overlapping pairs.
//
// Policy: points strictly on polygon edges and vertices are NOT
// considered interior. See the boundary policy note in
// harness_oracle.hpp.

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"

#include "harness_oracle.hpp"
#include "bitmap_test_utils.hpp"

#include <algorithm>
#include <cmath>

using namespace Slic3r;
using namespace Slic3r::arrangement;
using namespace Slic3r::arrangement::test_utils;

namespace {

// Axis-aligned unit-square at origin.
Polygon ho_square_mm(double x, double y, double w, double h)
{
    return Polygon(Points{
        Point(scaled<coord_t>(x),     scaled<coord_t>(y)),
        Point(scaled<coord_t>(x + w), scaled<coord_t>(y)),
        Point(scaled<coord_t>(x + w), scaled<coord_t>(y + h)),
        Point(scaled<coord_t>(x),     scaled<coord_t>(y + h))
    });
}

// L-shape with outer 40x40 bbox and a 20x20 notch in the top-right.
// The notch corner is at (20, 20)-(40, 40). The L occupies the
// bottom row (0..40 x 0..20) and the left column (0..20 x 0..40).
// A point at (30, 30) is INSIDE the bbox but OUTSIDE the L.
Polygon ho_l_shape_mm()
{
    return Polygon(Points{
        Point(scaled<coord_t>( 0.0), scaled<coord_t>( 0.0)),
        Point(scaled<coord_t>(40.0), scaled<coord_t>( 0.0)),
        Point(scaled<coord_t>(40.0), scaled<coord_t>(20.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>(20.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>(40.0)),
        Point(scaled<coord_t>( 0.0), scaled<coord_t>(40.0))
    });
}

// Annulus (ring) with outer radius 20, inner radius 10, approximated
// as regular polygons with 64 vertices each.
ExPolygon ho_annulus_mm()
{
    constexpr int N = 64;
    Points outer, inner;
    for (int i = 0; i < N; ++i) {
        double t = 2.0 * M_PI * (double)i / (double)N;
        outer.emplace_back(scaled<coord_t>(20.0 * std::cos(t)),
                           scaled<coord_t>(20.0 * std::sin(t)));
        inner.emplace_back(scaled<coord_t>(10.0 * std::cos(t)),
                           scaled<coord_t>(10.0 * std::sin(t)));
    }
    // Hole winding must be opposite of contour.
    std::reverse(inner.begin(), inner.end());
    ExPolygon ep;
    ep.contour = Polygon(outer);
    ep.holes   = Polygons{Polygon(inner)};
    return ep;
}

} // namespace

// ===========================================================================
// Phase 1.1 — point_in_polygon_oracle on analytical cases
// ===========================================================================

TEST_CASE("Oracle: point inside unit square is interior",
          "[HarnessOracle][c2-m2-oracle]")
{
    Polygon sq = ho_square_mm(0, 0, 10, 10);
    REQUIRE(point_in_polygon_oracle(5.0, 5.0, sq));     // center
    REQUIRE(point_in_polygon_oracle(0.1, 0.1, sq));     // near BL corner
    REQUIRE(point_in_polygon_oracle(9.9, 9.9, sq));     // near TR corner
    REQUIRE(point_in_polygon_oracle(9.9, 0.1, sq));     // near BR corner
    REQUIRE(point_in_polygon_oracle(0.1, 9.9, sq));     // near TL corner
}

TEST_CASE("Oracle: point outside unit square is not interior",
          "[HarnessOracle][c2-m2-oracle]")
{
    Polygon sq = ho_square_mm(0, 0, 10, 10);
    REQUIRE_FALSE(point_in_polygon_oracle(-1.0,  5.0, sq));  // left
    REQUIRE_FALSE(point_in_polygon_oracle(11.0,  5.0, sq));  // right
    REQUIRE_FALSE(point_in_polygon_oracle( 5.0, -1.0, sq));  // below
    REQUIRE_FALSE(point_in_polygon_oracle( 5.0, 11.0, sq));  // above
    REQUIRE_FALSE(point_in_polygon_oracle(-5.0, -5.0, sq));  // diagonal
    REQUIRE_FALSE(point_in_polygon_oracle(50.0, 50.0, sq));  // far
}

TEST_CASE("Oracle: point on unit-square vertex follows boundary policy (exterior)",
          "[HarnessOracle][c2-m2-oracle]")
{
    Polygon sq = ho_square_mm(0, 0, 10, 10);
    // Per policy, on-boundary points are not interior.
    REQUIRE_FALSE(point_in_polygon_oracle( 0.0,  0.0, sq));  // BL vertex
    REQUIRE_FALSE(point_in_polygon_oracle(10.0,  0.0, sq));  // BR vertex
    REQUIRE_FALSE(point_in_polygon_oracle(10.0, 10.0, sq));  // TR vertex
    REQUIRE_FALSE(point_in_polygon_oracle( 0.0, 10.0, sq));  // TL vertex
}

TEST_CASE("Oracle: point on unit-square edge follows boundary policy (exterior)",
          "[HarnessOracle][c2-m2-oracle]")
{
    Polygon sq = ho_square_mm(0, 0, 10, 10);
    // Mid-edge samples (exact rasterization lands here in practice).
    REQUIRE_FALSE(point_in_polygon_oracle( 5.0,  0.0, sq));  // bottom edge
    REQUIRE_FALSE(point_in_polygon_oracle( 5.0, 10.0, sq));  // top edge
    REQUIRE_FALSE(point_in_polygon_oracle( 0.0,  5.0, sq));  // left edge
    REQUIRE_FALSE(point_in_polygon_oracle(10.0,  5.0, sq));  // right edge
}

TEST_CASE("Oracle: point in L-shape notch is NOT interior despite being in bbox",
          "[HarnessOracle][c2-m2-oracle]")
{
    Polygon L = ho_l_shape_mm();
    // Inside the 40x40 bounding box but in the removed top-right quadrant.
    REQUIRE_FALSE(point_in_polygon_oracle(30.0, 30.0, L));
    REQUIRE_FALSE(point_in_polygon_oracle(25.0, 25.0, L));
    REQUIRE_FALSE(point_in_polygon_oracle(39.0, 39.0, L));
}

TEST_CASE("Oracle: point in L-shape solid body is interior",
          "[HarnessOracle][c2-m2-oracle]")
{
    Polygon L = ho_l_shape_mm();
    // Bottom row — solid for 0..20 y.
    REQUIRE(point_in_polygon_oracle(10.0, 10.0, L));
    REQUIRE(point_in_polygon_oracle(30.0, 10.0, L));
    // Left column — solid for 0..20 x.
    REQUIRE(point_in_polygon_oracle(10.0, 30.0, L));
    REQUIRE(point_in_polygon_oracle( 5.0, 35.0, L));
}

// ===========================================================================
// Phase 1.2 — annulus (hole detection)
// ===========================================================================

TEST_CASE("Oracle: point in annulus solid ring is interior",
          "[HarnessOracle][c2-m2-oracle]")
{
    ExPolygon ring = ho_annulus_mm();
    // Between inner r=10 and outer r=20 — should be interior.
    REQUIRE(point_in_expolygon_oracle(15.0,  0.0, ring));
    REQUIRE(point_in_expolygon_oracle( 0.0, 15.0, ring));
    REQUIRE(point_in_expolygon_oracle(-15.0, 0.0, ring));
    REQUIRE(point_in_expolygon_oracle( 0.0, -15.0, ring));
    // Near outer boundary (just inside).
    REQUIRE(point_in_expolygon_oracle(19.0,  0.0, ring));
}

TEST_CASE("Oracle: point in annulus center hole is NOT interior",
          "[HarnessOracle][c2-m2-oracle]")
{
    ExPolygon ring = ho_annulus_mm();
    // Inside the r=10 hole — should be exterior.
    REQUIRE_FALSE(point_in_expolygon_oracle( 0.0,  0.0, ring));
    REQUIRE_FALSE(point_in_expolygon_oracle( 5.0,  0.0, ring));
    REQUIRE_FALSE(point_in_expolygon_oracle( 0.0,  5.0, ring));
    REQUIRE_FALSE(point_in_expolygon_oracle(-5.0,  5.0, ring));
    // Just inside the inner boundary (the hole side).
    REQUIRE_FALSE(point_in_expolygon_oracle( 9.0,  0.0, ring));
}

TEST_CASE("Oracle: point outside annulus outer boundary is not interior",
          "[HarnessOracle][c2-m2-oracle]")
{
    ExPolygon ring = ho_annulus_mm();
    REQUIRE_FALSE(point_in_expolygon_oracle(25.0,  0.0, ring));
    REQUIRE_FALSE(point_in_expolygon_oracle( 0.0, 25.0, ring));
    REQUIRE_FALSE(point_in_expolygon_oracle(-25.0, 0.0, ring));
    REQUIRE_FALSE(point_in_expolygon_oracle(30.0, 30.0, ring));
}

// ===========================================================================
// Phase 1.3 — grid-sampled overlap oracle
// ===========================================================================

TEST_CASE("Oracle: non-overlapping squares report no overlap",
          "[HarnessOracle][c2-m2-oracle]")
{
    ExPolygon a;
    a.contour = ho_square_mm(0, 0, 10, 10);
    ExPolygon b;
    b.contour = ho_square_mm(20, 0, 10, 10);  // 10mm gap
    REQUIRE_FALSE(expolygons_overlap_oracle(a, b, 0.5));
    REQUIRE_FALSE(expolygons_overlap_oracle_symmetric(a, b, 0.5));
}

TEST_CASE("Oracle: touching squares do not overlap under boundary policy",
          "[HarnessOracle][c2-m2-oracle]")
{
    // Two squares sharing the edge x=10. Per policy, points on the
    // shared edge are exterior to both, so there is no interior
    // sample that lies in both.
    ExPolygon a;
    a.contour = ho_square_mm(0, 0, 10, 10);
    ExPolygon b;
    b.contour = ho_square_mm(10, 0, 10, 10);
    REQUIRE_FALSE(expolygons_overlap_oracle(a, b, 0.5));
}

TEST_CASE("Oracle: overlapping squares are detected",
          "[HarnessOracle][c2-m2-oracle]")
{
    // A at (0,0)-(10,10); B at (5,5)-(15,15). Intersection is 5x5.
    ExPolygon a;
    a.contour = ho_square_mm(0, 0, 10, 10);
    ExPolygon b;
    b.contour = ho_square_mm(5, 5, 10, 10);
    REQUIRE(expolygons_overlap_oracle(a, b, 0.5));
    REQUIRE(expolygons_overlap_oracle_symmetric(a, b, 0.5));
}

TEST_CASE("Oracle: containment counts as overlap",
          "[HarnessOracle][c2-m2-oracle]")
{
    // B is fully inside A.
    ExPolygon a;
    a.contour = ho_square_mm(0, 0, 20, 20);
    ExPolygon b;
    b.contour = ho_square_mm(5, 5, 10, 10);
    REQUIRE(expolygons_overlap_oracle(a, b, 0.5));
    REQUIRE(expolygons_overlap_oracle_symmetric(a, b, 0.5));
}

TEST_CASE("Oracle: L-shape notch allows a small square to fit non-overlapping",
          "[HarnessOracle][c2-m2-oracle]")
{
    // L-shape has a 20x20 notch at (20, 20)-(40, 40). A 10x10 square
    // placed inside the notch (say, at (25, 25)-(35, 35)) doesn't
    // overlap the L.
    ExPolygon L;
    L.contour = ho_l_shape_mm();
    ExPolygon inner;
    inner.contour = ho_square_mm(25, 25, 10, 10);
    REQUIRE_FALSE(expolygons_overlap_oracle(L, inner, 0.25));
    REQUIRE_FALSE(expolygons_overlap_oracle_symmetric(L, inner, 0.25));
}

TEST_CASE("Oracle: L-shape concavity detection — square in body overlaps",
          "[HarnessOracle][c2-m2-oracle]")
{
    // A square placed in the SOLID part of the L (bottom row) should
    // overlap the L.
    ExPolygon L;
    L.contour = ho_l_shape_mm();
    ExPolygon inside;
    inside.contour = ho_square_mm(5, 5, 10, 10);
    REQUIRE(expolygons_overlap_oracle(L, inside, 0.25));
}

// ===========================================================================
// Phase 1.4 — cross-check: oracle vs no_overlap helper
// ===========================================================================
//
// These tests validate that test_utils::no_overlap (the C1 polygon-
// intersection-based helper) agrees with the oracle across a battery
// of analytical cases. If they diverge, the C1 helper has a bug that
// every downstream test silently trusts.

TEST_CASE("Cross-check: no_overlap agrees with oracle on non-overlap cases",
          "[HarnessOracle][c2-m2-oracle][cross-check]")
{
    auto make_ap = [](const ExPolygon& poly, int bed_idx, double tx_mm, double ty_mm) {
        ArrangePolygon ap;
        ap.poly = poly;
        ap.bed_idx = bed_idx;
        ap.rotation = 0.0;
        ap.translation = Vec2crd{scaled<coord_t>(tx_mm), scaled<coord_t>(ty_mm)};
        return ap;
    };

    // Five non-overlapping pairs.
    struct Case { double tx_b; double ty_b; };
    Case cases[] = {
        {20.0, 0.0},    // gap on x
        {0.0, 20.0},    // gap on y
        {15.0, 15.0},   // diagonal gap
        {100.0, 0.0},   // far
        {-20.0, 0.0},   // negative gap
    };
    ExPolygon sq;
    sq.contour = ho_square_mm(0, 0, 10, 10);
    for (const auto& c : cases) {
        ArrangePolygons items;
        items.push_back(make_ap(sq, 0, 0.0, 0.0));
        items.push_back(make_ap(sq, 0, c.tx_b, c.ty_b));
        // Both helpers should report no overlap.
        REQUIRE(test_utils::no_overlap(items));

        // Oracle check on the same pair (transformed polygons).
        ExPolygon a_t = items[0].transformed_poly();
        ExPolygon b_t = items[1].transformed_poly();
        REQUIRE_FALSE(expolygons_overlap_oracle(a_t, b_t, 0.5));
    }
}

TEST_CASE("Cross-check: no_overlap agrees with oracle on overlap cases",
          "[HarnessOracle][c2-m2-oracle][cross-check]")
{
    auto make_ap = [](const ExPolygon& poly, int bed_idx, double tx_mm, double ty_mm) {
        ArrangePolygon ap;
        ap.poly = poly;
        ap.bed_idx = bed_idx;
        ap.rotation = 0.0;
        ap.translation = Vec2crd{scaled<coord_t>(tx_mm), scaled<coord_t>(ty_mm)};
        return ap;
    };

    // Five overlapping pairs.
    struct Case { double tx_b; double ty_b; };
    Case cases[] = {
        {5.0, 0.0},      // halfway overlap on x
        {0.0, 5.0},      // halfway overlap on y
        {2.0, 2.0},      // heavy overlap
        {0.0, 0.0},      // exact overlap
        {8.0, 8.0},      // corner nibble
    };
    ExPolygon sq;
    sq.contour = ho_square_mm(0, 0, 10, 10);
    for (const auto& c : cases) {
        ArrangePolygons items;
        items.push_back(make_ap(sq, 0, 0.0, 0.0));
        items.push_back(make_ap(sq, 0, c.tx_b, c.ty_b));
        // Both helpers should report overlap.
        REQUIRE_FALSE(test_utils::no_overlap(items));

        ExPolygon a_t = items[0].transformed_poly();
        ExPolygon b_t = items[1].transformed_poly();
        REQUIRE(expolygons_overlap_oracle(a_t, b_t, 0.5));
    }
}

TEST_CASE("Cross-check: L-shape concavity is correctly exploited by oracle",
          "[HarnessOracle][c2-m2-oracle][cross-check]")
{
    // Two L-shapes rotated 180° and placed so their notches mesh
    // into each other's bodies. Classical interlock. The oracle
    // should report NO overlap (they fit into each other's
    // concavities) and so should no_overlap.
    ExPolygon L1;
    L1.contour = ho_l_shape_mm();

    // Second L is the same shape rotated 180° and offset.
    // Rotated L's notch is at the bottom-left; placing it offset
    // so the two notches face each other produces a rectangle.
    // L1 occupies (0,0)-(40,40) minus (20,20)-(40,40) notch.
    // Rotated L at (40,40) rotated 180° would put its notch back at
    // (20,20)-(40,40) relative to its original frame — but after the
    // rotation+offset, the second L occupies the TOP-RIGHT slot of
    // a 40x40 square.
    //
    // Easier construction: place L2 by hand so it fills L1's notch
    // directly. L2 is a 20x20 solid square at (20, 20)-(40, 40) —
    // not an L at all, but the complement of L1's notch. The
    // oracle should report no overlap of L1 with this filler, and
    // L1 ∪ filler = a full 40x40 square.
    ExPolygon filler;
    filler.contour = ho_square_mm(20, 20, 20, 20);

    // Policy: shared edge is not interior, so the oracle reports
    // no overlap of L1's boundary with filler's boundary.
    REQUIRE_FALSE(expolygons_overlap_oracle(L1, filler, 0.25));

    // And no_overlap (polygon intersection) agrees.
    auto make_ap = [](const ExPolygon& poly) {
        ArrangePolygon ap;
        ap.poly = poly;
        ap.bed_idx = 0;
        ap.rotation = 0.0;
        ap.translation = Vec2crd{0, 0};
        return ap;
    };
    ArrangePolygons items{make_ap(L1), make_ap(filler)};
    REQUIRE(test_utils::no_overlap(items));
}
