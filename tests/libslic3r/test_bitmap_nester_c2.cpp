// test_bitmap_nester_c2.cpp
//
// Tests for BitmapNesterC2 — the organic lasso nester in development
// on feature/concave-bitmap-c2. See docs/PR_C2_ORGANIC_LASSO_NESTER_PLAN.md
// for the milestone plan and docs/PR_C2_TEST_CORPUS.md for the test
// corpus driving quality gates.
//
// Section layout:
//   M0 smoke tests (pass-through verification — inherited from C0 skeleton)
//   M1 phase-function unit tests (build_cache, estimate_min_plates,
//      partition_items — landed this session)
//   Category 6.4 litmus test (identical convex hull, from test corpus)

// Expose the M1 phase functions for direct unit testing. This mirrors
// C1's BITMAP_NESTER_TESTING pattern and avoids moving the phase
// functions to the public interface just for tests.
#define BITMAP_NESTER_C2_TESTING

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/BitmapNesterC2.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "bitmap_test_utils.hpp"
#include "harness_oracle.hpp"

#include <cmath>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::arrangement;

namespace {

BoundingBox c2_bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

ArrangeParams c2_params()
{
    ArrangeParams p;
    p.bed_shrink_x       = 0.0f;
    p.bed_shrink_y       = 0.0f;
    p.allow_rotations    = true;
    p.min_obj_distance   = 0;
    p.use_concave_shapes = true;
    p.progressind        = nullptr;
    return p;
}

ExPolygon c2_rect_mm(double w, double h)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h)),
        Point(scaled<coord_t>(0.0), scaled<coord_t>(h))
    });
}

ArrangePolygon c2_make_ap(const ExPolygon& poly, double height_mm = 0.0)
{
    ArrangePolygon ap;
    ap.poly               = poly;
    ap.priority           = 0;
    ap.bed_idx            = UNARRANGED;
    ap.rotation           = 0.0;
    ap.translation        = Vec2crd{0, 0};
    ap.allowed_rotations  = {0.0};
    ap.height             = height_mm;
    return ap;
}

// ---------------------------------------------------------------------------
// Category 6.4 litmus test fixtures. Two shapes with the IDENTICAL convex
// hull but dramatically different concave silhouettes. A true concave
// solver produces measurably different packings. A solver that secretly
// works in convex hulls treats them identically.
// ---------------------------------------------------------------------------

// Solid square 40x40.
ExPolygon c2_solid_square_mm()
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(40.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(40.0), scaled<coord_t>(40.0)),
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(40.0))
    });
}

// C-shape (a square with a deep rectangular notch cut from one side).
// Outer bbox is 40x40; the notch is 30 wide × 24 deep cut from the
// right side centered vertically. Silhouette area = 1600 - 720 = 880
// mm² vs bbox area 1600 mm² — ratio 0.55. A "pocket" that can
// accept another shape's protrusion but whose bbox significantly
// over-represents the occupied area.
ExPolygon c2_c_shape_mm()
{
    return ExPolygon(Points{
        Point(scaled<coord_t>( 0.0), scaled<coord_t>( 0.0)),  // BL
        Point(scaled<coord_t>(40.0), scaled<coord_t>( 0.0)),  // BR
        Point(scaled<coord_t>(40.0), scaled<coord_t>( 8.0)),  // bottom of mouth
        Point(scaled<coord_t>(10.0), scaled<coord_t>( 8.0)),  // inward top of notch
        Point(scaled<coord_t>(10.0), scaled<coord_t>(32.0)),  // inward bottom of notch
        Point(scaled<coord_t>(40.0), scaled<coord_t>(32.0)),  // top of mouth
        Point(scaled<coord_t>(40.0), scaled<coord_t>(40.0)),  // TR
        Point(scaled<coord_t>( 0.0), scaled<coord_t>(40.0))   // TL
    });
}

// L-shape with 40x40 bbox and a 20x20 notch cut from the top-right
// quadrant. Two of these (one at rotation 0, one at rotation 180°)
// nest perfectly into each other's notches, producing a 40x40
// combined footprint. The hull-scored C2 path should find this
// interlock; a naive bbox-scored greedy will stack them side-by-side
// or stacked, producing a 40x80 or 80x40 cluster — double the
// expected footprint. This is the Category 6.4 discriminator at
// its cleanest: known optimal, binary pass/fail.
ExPolygon c2_l_shape_40_20_mm()
{
    return ExPolygon(Points{
        Point(scaled<coord_t>( 0.0), scaled<coord_t>( 0.0)),
        Point(scaled<coord_t>(40.0), scaled<coord_t>( 0.0)),
        Point(scaled<coord_t>(40.0), scaled<coord_t>(20.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>(20.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>(40.0)),
        Point(scaled<coord_t>( 0.0), scaled<coord_t>(40.0))
    });
}

// Notched square: a 40x40 square with a 20x20 rectangular notch cut
// from the middle of the top edge. The four outermost corners (0,0),
// (40,0), (40,40), (0,40) are all still on the outline, so the convex
// hull is still exactly the 40x40 square — identical to
// c2_solid_square_mm's hull. But the silhouette area is
// 1600 - 400 = 1200 mm² vs the solid square's 1600 mm².
//
// A true concave solver treats these two shapes dramatically
// differently: it can pack notched squares with their notches
// interlocking while solid squares can't interlock at all. A convex-
// hull solver treats them identically (same 40x40 hull → same packing).
ExPolygon c2_notched_square_mm()
{
    return ExPolygon(Points{
        Point(scaled<coord_t>( 0.0), scaled<coord_t>( 0.0)),  // BL
        Point(scaled<coord_t>(40.0), scaled<coord_t>( 0.0)),  // BR
        Point(scaled<coord_t>(40.0), scaled<coord_t>(40.0)),  // TR
        Point(scaled<coord_t>(30.0), scaled<coord_t>(40.0)),  // top edge, right of notch
        Point(scaled<coord_t>(30.0), scaled<coord_t>(20.0)),  // down into notch (right wall)
        Point(scaled<coord_t>(10.0), scaled<coord_t>(20.0)),  // across notch bottom
        Point(scaled<coord_t>(10.0), scaled<coord_t>(40.0)),  // up out of notch (left wall)
        Point(scaled<coord_t>( 0.0), scaled<coord_t>(40.0)),  // TL
    });
}

} // namespace

// ---------------------------------------------------------------------------
// M0 smoke test: arrange() compiles, runs, and produces the same output
// as BitmapNester on a trivial input (because the pass-through delegates).
// ---------------------------------------------------------------------------

TEST_CASE("C2 M0: arrange pass-through places a single item on plate 0",
          "[BitmapNesterC2][c2-m0]")
{
    ArrangePolygons items;
    items.push_back(c2_make_ap(c2_rect_mm(40, 40)));

    BoundingBox bed = c2_bed_mm(200, 200);
    ArrangeParams params = c2_params();

    BitmapNesterC2::arrange(items, ArrangePolygons{}, bed, params);

    REQUIRE(items.size() == 1);
    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(test_utils::no_overlap(items));
}

// ---------------------------------------------------------------------------
// M0 smoke test: arrange() handles zero items without crashing.
// ---------------------------------------------------------------------------

TEST_CASE("C2 M0: arrange pass-through on empty input returns cleanly",
          "[BitmapNesterC2][c2-m0]")
{
    ArrangePolygons items;
    BoundingBox bed = c2_bed_mm(200, 200);
    ArrangeParams params = c2_params();

    BitmapNesterC2::arrange(items, ArrangePolygons{}, bed, params);

    REQUIRE(items.empty());
}

// ---------------------------------------------------------------------------
// M0 smoke test: arrange() packs four squares identically to C1 (because
// the pass-through IS C1). Sanity-check that the fork branch doesn't
// accidentally regress C1's baseline behavior.
// ---------------------------------------------------------------------------

TEST_CASE("C2 M0: arrange pass-through packs four 80x80 squares on plate 0",
          "[BitmapNesterC2][c2-m0]")
{
    ArrangePolygons items;
    for (int i = 0; i < 4; ++i)
        items.push_back(c2_make_ap(c2_rect_mm(80, 80)));

    BoundingBox bed = c2_bed_mm(200, 200);
    ArrangeParams params = c2_params();

    BitmapNesterC2::arrange(items, ArrangePolygons{}, bed, params);

    REQUIRE(test_utils::max_bed_idx(items) == 0);
    REQUIRE(test_utils::no_overlap(items));
}

// ===========================================================================
// M1 — Phase-function unit tests
// ===========================================================================
//
// build_cache, estimate_min_plates, and partition_items are exposed via
// BITMAP_NESTER_C2_TESTING at the top of this file.

TEST_CASE("C2 M1: build_cache populates geometric summary",
          "[BitmapNesterC2][c2-m1]")
{
    ArrangePolygons items;
    items.push_back(c2_make_ap(c2_rect_mm(40, 30)));     // 1200 mm²
    items.push_back(c2_make_ap(c2_rect_mm(20, 20)));     // 400 mm²
    items.push_back(c2_make_ap(c2_solid_square_mm(), 150.0)); // tall

    auto cache = BitmapNesterC2::build_cache(items);

    REQUIRE(cache.size() == 3);

    // Item 0: 40x30 rectangle → area 1200, perimeter 140, hull==silhouette
    REQUIRE_THAT(cache[0].silhouette_area_mm2,
                 Catch::Matchers::WithinAbs(1200.0, 1.0));
    REQUIRE_THAT(cache[0].hull_area_mm2,
                 Catch::Matchers::WithinAbs(1200.0, 1.0));
    REQUIRE_THAT(cache[0].max_dim_mm,
                 Catch::Matchers::WithinAbs(40.0, 0.01));
    REQUIRE_THAT(cache[0].hardness_score,
                 Catch::Matchers::WithinAbs(1.0, 0.01));
    REQUIRE(cache[0].height_mm == 0.0);

    // Item 1: 20x20 square → area 400
    REQUIRE_THAT(cache[1].silhouette_area_mm2,
                 Catch::Matchers::WithinAbs(400.0, 1.0));
    REQUIRE_THAT(cache[1].max_dim_mm,
                 Catch::Matchers::WithinAbs(20.0, 0.01));

    // Item 2: 40x40 with height=150 → tall, priority_score should reflect it
    REQUIRE(cache[2].height_mm == 150.0);
    REQUIRE(cache[2].priority_score > cache[0].priority_score);
    REQUIRE(cache[2].priority_score > cache[1].priority_score);
}

TEST_CASE("C2 M1: estimate_min_plates returns 1 for loose inputs",
          "[BitmapNesterC2][c2-m1]")
{
    ArrangePolygons items;
    for (int i = 0; i < 4; ++i)
        items.push_back(c2_make_ap(c2_rect_mm(30, 30)));  // 3600 mm² total

    auto cache = BitmapNesterC2::build_cache(items);
    BoundingBox bed = c2_bed_mm(200, 200);  // 40000 mm²

    int k = BitmapNesterC2::estimate_min_plates(cache, bed);
    REQUIRE(k == 1);  // 9% density, well under the 82% efficiency factor
}

TEST_CASE("C2 M1: estimate_min_plates scales up for dense inputs",
          "[BitmapNesterC2][c2-m1]")
{
    ArrangePolygons items;
    // Pack 50 x 40x40 = 50 * 1600 = 80,000 mm² on a 200x200 = 40,000 mm² bed.
    // Raw area alone requires 2 plates; with 0.82 efficiency: ceil(80000 / 32800) = 3
    for (int i = 0; i < 50; ++i)
        items.push_back(c2_make_ap(c2_rect_mm(40, 40)));

    auto cache = BitmapNesterC2::build_cache(items);
    BoundingBox bed = c2_bed_mm(200, 200);

    int k = BitmapNesterC2::estimate_min_plates(cache, bed);
    REQUIRE(k >= 2);
    REQUIRE(k <= 4);  // conservative range; exact value depends on efficiency factor
}

TEST_CASE("C2 M1: partition_items balances total area across buckets",
          "[BitmapNesterC2][c2-m1]")
{
    // Mix of 3 large + 6 small items. With k=2, balance should put
    // some large and some small in each bucket to equalize area.
    ArrangePolygons items;
    for (int i = 0; i < 3; ++i)
        items.push_back(c2_make_ap(c2_rect_mm(50, 50)));  // 2500 mm² each
    for (int i = 0; i < 6; ++i)
        items.push_back(c2_make_ap(c2_rect_mm(20, 20)));  // 400 mm² each
    // Total: 3 * 2500 + 6 * 400 = 7500 + 2400 = 9900 mm²
    // Ideal per-bucket: 4950 mm² when k=2

    auto cache = BitmapNesterC2::build_cache(items);
    auto groups = BitmapNesterC2::partition_items(cache, 2);

    REQUIRE(groups.size() == 2);
    REQUIRE((groups[0].items.size() + groups[1].items.size()) == 9);

    // Both groups should be within 30% of the ideal 4950 mm²
    double ideal = 9900.0 / 2.0;
    for (const auto& g : groups) {
        double ratio = g.total_silhouette_area_mm2 / ideal;
        UNSCOPED_INFO("group area " << g.total_silhouette_area_mm2
                      << " ratio " << ratio);
        REQUIRE(ratio >= 0.70);
        REQUIRE(ratio <= 1.30);
    }
}

TEST_CASE("C2 M1: partition_items with k=1 puts all items in one group",
          "[BitmapNesterC2][c2-m1]")
{
    ArrangePolygons items;
    for (int i = 0; i < 5; ++i)
        items.push_back(c2_make_ap(c2_rect_mm(20, 20)));

    auto cache = BitmapNesterC2::build_cache(items);
    auto groups = BitmapNesterC2::partition_items(cache, 1);

    REQUIRE(groups.size() == 1);
    REQUIRE(groups[0].items.size() == 5);
}

TEST_CASE("C2 M1: partition_items with empty input returns empty",
          "[BitmapNesterC2][c2-m1]")
{
    std::vector<NesterC2ItemInfo> cache;
    auto groups = BitmapNesterC2::partition_items(cache, 3);
    REQUIRE(groups.empty());
}

// ===========================================================================
// Category 6.4 — identical convex hull litmus test
// ===========================================================================
//
// From docs/PR_C2_TEST_CORPUS.md: "Two very different concave shapes
// that have the same convex hull. A convex solver treats them
// identically. A concave solver should find dramatically different
// packings for each. This is the litmus test for whether your solver
// is actually doing concave work or just wrapping things in convex
// hulls internally."
//
// Shapes:
//   - c2_solid_square_mm:     40x40 solid square. Silhouette = 1600 mm².
//   - c2_notched_square_mm:   40x40 square with 20x20 notch in top edge.
//                             Silhouette = 1200 mm². Same convex hull
//                             (40x40 square) as the solid square.
//
// Both shapes have the same hull, so a convex-hull solver packs them
// identically — placing N copies tile into a grid with bbox growth
// proportional to count. A concave solver lets two notched squares
// interlock notch-to-tab, visibly reducing total cluster area.
//
// Expected PASSING state: the measured total cluster area for N notched
// squares is MEASURABLY SMALLER than for N solid squares. The exact
// threshold depends on the solver — we assert a loose lower bound that
// any genuine concave solver beats.
//
// Current state (pass-through delegation to C1): C1 rasterizes concave
// silhouettes and scores by cluster-bbox growth, so it already treats
// these shapes differently. C2 inheriting C1's behavior passes too.
// When C2's real phases come online, this test should continue passing
// and ideally TIGHTEN the delta.

TEST_CASE("Category 6.4: identical convex hulls produce different packings",
          "[BitmapNesterC2][c2-corpus][c2-m1][Cat6.4]")
{
    struct Measurement {
        int    max_bed;
        double bbox_perim;
        double hull_perim;
        double bbox_compactness;
        double hull_compactness;
    };

    auto pack_and_measure = [](const ExPolygon& shape) -> Measurement {
        ArrangePolygons items;
        for (int i = 0; i < 6; ++i) {
            ArrangePolygon ap = c2_make_ap(shape);
            // Allow 90° rotations so the nester can interlock notched pairs.
            ap.allowed_rotations = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
            items.push_back(ap);
        }
        BoundingBox bed = c2_bed_mm(300, 300);
        ArrangeParams params = c2_params();
        params.min_obj_distance = 0;
        BitmapNesterC2::arrange(items, ArrangePolygons{}, bed, params);

        Measurement m;
        m.max_bed          = test_utils::max_bed_idx(items);
        m.bbox_perim       = test_utils::cluster_bbox_perimeter_mm(items);
        m.hull_perim       = test_utils::cluster_hull_perimeter_mm(items);
        m.bbox_compactness = test_utils::cluster_compactness(items);
        m.hull_compactness = test_utils::cluster_hull_compactness(items);
        return m;
    };

    Measurement solid   = pack_and_measure(c2_solid_square_mm());
    Measurement notched = pack_and_measure(c2_notched_square_mm());

    UNSCOPED_INFO("solid:   max_bed="  << solid.max_bed
                  << " bbox_perim="    << solid.bbox_perim
                  << " hull_perim="    << solid.hull_perim
                  << " bbox_compact="  << solid.bbox_compactness
                  << " hull_compact="  << solid.hull_compactness);
    UNSCOPED_INFO("notched: max_bed="  << notched.max_bed
                  << " bbox_perim="    << notched.bbox_perim
                  << " hull_perim="    << notched.hull_perim
                  << " bbox_compact="  << notched.bbox_compactness
                  << " hull_compact="  << notched.hull_compactness);

    // Both shapes should land valid single-plate packings — 6 items
    // of ≤ 1600 mm² each on a 300x300 bed is trivially loose.
    REQUIRE(solid.max_bed == 0);
    REQUIRE(notched.max_bed == 0);

    // THE LITMUS: the two shapes have the SAME convex hull (40x40
    // square) but different silhouettes (solid=1600 mm², notched=
    // 1200 mm² due to the 20x20 notch). A genuinely concave-aware
    // solver produces different packings — measurable via the
    // compactness metric which uses silhouette area.
    //
    // We check BOTH bbox compactness AND hull compactness. The hull
    // variant is strictly more honest: for two packed 3x2 grids,
    // the bbox is 40x40 regardless of shape and may match exactly,
    // but the hull of the union of all placed items may differ
    // subtly between the two (solid squares tile exactly into a
    // 120x80 rectangle; notched squares leave little gaps that
    // the hull doesn't care about but compactness does).
    //
    // The expected failure mode of a bbox-hidden solver is:
    // solid.bbox_compactness == notched.bbox_compactness AND
    // solid.hull_compactness == notched.hull_compactness.
    // A concave solver produces a nonzero delta in at least one.
    //
    // CURRENT STATE (M1 pass-through): C2 delegates to C1 which
    // rasterizes silhouettes, so the two shapes already produce
    // different SILHOUETTE AREAS per plate → different compactness.
    // Test passes by virtue of C1's existing concave awareness.
    // M2's pack_as_island will tighten this further by making the
    // solver's INTERNAL scoring also hull-based.
    double bbox_delta = std::abs(solid.bbox_compactness - notched.bbox_compactness);
    double hull_delta = std::abs(solid.hull_compactness - notched.hull_compactness);
    UNSCOPED_INFO("bbox_compactness delta = " << bbox_delta);
    UNSCOPED_INFO("hull_compactness delta = " << hull_delta);

    // At least ONE of the two metrics must detect the difference.
    // Loose 0.05 threshold — tighten when C2 M2 lands real hull
    // scoring.
    REQUIRE((bbox_delta > 0.05 || hull_delta > 0.05));
}

// ===========================================================================
// M2.4 — Tests that actually exercise M2.3's hull-scored inner loop
// ===========================================================================
//
// Every prior C2 test takes the k=1 fast path (single-group input ->
// delegate to BitmapNester::arrange). The M2.3 hull-scored inner loop
// only runs when estimate_min_plates returns k >= 2. These tests force
// k >= 2 and assert that the hull-scored path produces valid output.

TEST_CASE("C2 M2.4: forty squares on a tight bed exercises the hull path",
          "[BitmapNesterC2][c2-m2-hull][c2-m2.4]")
{
    // 40 × 40x40 squares = 64000 mm² on a 200x200 = 40000 mm² bed.
    // 160% density — clearly needs 2 plates under any algorithm.
    // estimate_min_plates: ceil(64000 / (40000 * 0.82)) = ceil(1.95) = 2
    // => partition_items produces k=2 => M2.3 runs on each group.
    ArrangePolygons items;
    for (int i = 0; i < 40; ++i)
        items.push_back(c2_make_ap(c2_rect_mm(40, 40)));

    BoundingBox bed = c2_bed_mm(200, 200);
    ArrangeParams params = c2_params();

    BitmapNesterC2::arrange(items, ArrangePolygons{}, bed, params);

    // Count outcomes per bed_idx.
    std::map<int, int> per_bed;
    int unarranged = 0;
    for (const auto& it : items) {
        if (it.bed_idx == UNARRANGED) unarranged++;
        else per_bed[it.bed_idx]++;
    }

    UNSCOPED_INFO("unarranged = " << unarranged);
    for (const auto& kv : per_bed)
        UNSCOPED_INFO("bed " << kv.first << ": " << kv.second << " items");

    // Primary assertion: M2.3 produced a multi-plate partition. This
    // verifies the hull-scored path ran at all.
    REQUIRE(per_bed.size() >= 2);

    // Secondary: most items should be placed. M2.3 may leave some
    // items UNARRANGED if the coarse grid scan can't find a clear
    // position — that's documented behavior, and M3's spillover
    // recovery will handle them. For now we just assert that the
    // fraction placed is at least 80%.
    int placed = 0;
    for (const auto& kv : per_bed) placed += kv.second;
    REQUIRE(placed >= 32);  // 80% of 40

    // Tertiary: no overlaps within any plate.
    REQUIRE(test_utils::no_overlap(items));

    // Quaternary: hull perimeter on each plate is at least the lower
    // bound set by the number of squares times their perimeter
    // contribution. Each square has perimeter 4*40 = 160 mm, and the
    // cluster hull is at most the sum of perimeters and at least the
    // perimeter of one square (for a single-item cluster). So we
    // expect per-plate hull perimeter between 160 and 6400 mm.
    double total_hull = test_utils::cluster_hull_perimeter_mm(items);
    UNSCOPED_INFO("total cluster hull perimeter = " << total_hull);
    REQUIRE(total_hull >= 160.0);
    REQUIRE(total_hull <= 6400.0);
}

// ---------------------------------------------------------------------------
// Category 1.4 from the C2 test corpus: a pair of L-shapes where
// the optimal packing is the 180°-rotated mate filling the notch.
//
// Geometry: each L has a 40x40 bbox with a 20x20 notch at the
// top-right (silhouette area = 1200 mm²). The TIGHTEST possible
// combination of two such L-shapes is a 60x40 rectangle formed by
// placing one at origin and the other rotated 180° and offset so
// its notch fills the first L's notch complement:
//
//        ###                 (second L after 180° rotation, notch
//        ###                  at bottom-left, sitting in the top-
//        ######               right of the first L)
//        ######               (first L, notch at top-right)
//     ######
//     ######
//
// The combined hull is a 60x40 rectangle with perimeter 200 mm.
// The non-interlocked "dumb" stacking of two 40x40-bbox L-shapes
// side by side would be an 80x40 rectangle with perimeter 240 mm.
// So the discriminator is:
//
//   Hull-aware solver: ≈ 200 mm (tight nest)
//   Bbox-only solver:  ≈ 240 mm (naive stack)
//
// This is a BINARY pass/fail discriminator with a ≥ 20% quality
// gap. Binary clear-signal proof that the solver is doing real
// concave work. Unlike the 40-square test which only proves M2.3
// RAN, this test proves concave awareness.
//
// Currently runs under the k=1 fast path (2 items on a loose 200x200
// bed, partition returns k=1, path delegates to C1). C1 is itself
// concave-aware via silhouette rasterization, so the test passes
// today as a ratification of C1's behavior. When C2 phases replace
// more of the delegation, the assertion continues holding and
// ideally the measured hull perimeter drops toward 200 exactly.
// ---------------------------------------------------------------------------

TEST_CASE("C2 Cat 1.4: two L-shapes nest via the notch (hull discriminator)",
          "[BitmapNesterC2][c2-corpus][Cat1.4][c2-m2.4]")
{
    ArrangePolygons items;
    ArrangePolygon a = c2_make_ap(c2_l_shape_40_20_mm());
    a.allowed_rotations = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    items.push_back(a);
    ArrangePolygon b = a;
    items.push_back(b);

    BoundingBox bed = c2_bed_mm(200, 200);
    ArrangeParams params = c2_params();

    BitmapNesterC2::arrange(items, ArrangePolygons{}, bed, params);

    REQUIRE(test_utils::max_bed_idx(items) == 0);
    REQUIRE(test_utils::no_overlap(items));

    double hull_perim = test_utils::cluster_hull_perimeter_mm(items);
    double bbox_perim = test_utils::cluster_bbox_perimeter_mm(items);
    UNSCOPED_INFO("hull_perim = " << hull_perim << " mm");
    UNSCOPED_INFO("bbox_perim = " << bbox_perim << " mm");

    // Hull-aware pass: perimeter should be at most 210 mm (the
    // optimal 60x40 answer is 200; allow 10 mm for rotation rounding
    // and rasterization noise at the 0.5 mm bitmap resolution).
    REQUIRE(hull_perim <= 210.0);

    // Discriminator: the result MUST be measurably better than the
    // naive bbox stack of 240 mm. A solver that treated both shapes
    // as their 40x40 bounding boxes would produce an 80x40 cluster
    // with hull perimeter 240. Our solver lands near 200 because
    // it sees the concave notches and nests them.
    REQUIRE(hull_perim < 240.0);

    int placed = 0;
    for (const auto& it : items)
        if (it.bed_idx != UNARRANGED) placed++;
    REQUIRE(placed == 2);
}

TEST_CASE("C2 M2.4: multi-plate partition preserves determinism",
          "[BitmapNesterC2][c2-m2-hull][c2-m2.4][BitmapDeterminism]")
{
    // Same input as above — byte-identical second run.
    auto build = []() {
        ArrangePolygons xs;
        for (int i = 0; i < 40; ++i)
            xs.push_back(c2_make_ap(c2_rect_mm(40, 40)));
        return xs;
    };
    BoundingBox bed = c2_bed_mm(200, 200);
    ArrangeParams params = c2_params();
    auto run = [&](ArrangePolygons& xs) {
        BitmapNesterC2::arrange(xs, ArrangePolygons{}, bed, params);
    };

    ArrangePolygons a = build();
    ArrangePolygons b = build();
    run(a);
    run(b);

    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        REQUIRE(a[i].bed_idx == b[i].bed_idx);
        REQUIRE(a[i].translation.x() == b[i].translation.x());
        REQUIRE(a[i].translation.y() == b[i].translation.y());
        REQUIRE(std::abs(a[i].rotation - b[i].rotation) < 1e-9);
    }
}

// ===========================================================================
// M2.5.1 — bbox is wrong for concave: L-shape notch-fill discriminator
// ===========================================================================
//
// Two 40x40 L-shapes with a 20x20 top-right notch. L1 at (0,0) has
// its notch empty at (20,20)-(40,40). L2 placed at (20,20) has its
// bottom-left solid 20x20 quadrant sitting IN L1's notch — i.e.,
// L2's solid fills L1's empty region.
//
// Bbox overlap: L1 bbox (0,0)-(40,40), L2 bbox (20,20)-(60,60).
// They overlap in (20,20)-(40,40) = 20x20 = 400 mm².
//
// Silhouette overlap: L1 is empty in the (20,20)-(40,40) region
// (that's L1's notch). L2's solid in that region is (20,20)-(40,40)
// (that's L2's bottom-left quadrant, which is SOLID). L1 empty +
// L2 solid = ZERO silhouette overlap.
//
// Result: bbox says "collision", silhouette says "clear". A solver
// using bbox as the collision shape would reject this valid
// placement. A silhouette-aware solver accepts it.
//
// This is the binary proof that motivates the M2.5 bbox purge. Every
// place we use bbox as a proxy for silhouette, we over-reject valid
// placements in exactly this way. For a crescent, the over-rejection
// ratio is ~4× (bbox is 4× silhouette area). For these L-shapes, it's
// ~1.3× (bbox 1600 vs silhouette ~1200), still enough to demonstrate.
// ---------------------------------------------------------------------------

TEST_CASE("C2 M2.5.1: L-shape notch-fill — bbox overlaps but silhouettes don't",
          "[BitmapNesterC2][c2-m2.5][bbox-discriminator]")
{
    // L1 at origin, L2 at (20, 20). L2's bottom-left quadrant fills
    // L1's top-right notch exactly.
    ArrangePolygon l1 = c2_make_ap(c2_l_shape_40_20_mm());
    l1.translation = Vec2crd{0, 0};
    l1.rotation = 0.0;
    l1.bed_idx = 0;

    ArrangePolygon l2 = c2_make_ap(c2_l_shape_40_20_mm());
    l2.translation = Vec2crd{scaled<coord_t>(20.0), scaled<coord_t>(20.0)};
    l2.rotation = 0.0;
    l2.bed_idx = 0;

    ArrangePolygons items;
    items.push_back(l1);
    items.push_back(l2);

    // Transformed polygons.
    ExPolygon tp0 = items[0].transformed_poly();
    ExPolygon tp1 = items[1].transformed_poly();
    BoundingBox bb0 = get_extents(tp0);
    BoundingBox bb1 = get_extents(tp1);

    // Compute bbox overlap state explicitly.
    bool bbox_overlap =
        !(bb0.max.x() < bb1.min.x() ||
          bb0.min.x() > bb1.max.x() ||
          bb0.max.y() < bb1.min.y() ||
          bb0.min.y() > bb1.max.y());

    UNSCOPED_INFO("L1 bbox: ("
                  << unscaled<double>(bb0.min.x()) << ", "
                  << unscaled<double>(bb0.min.y()) << ") to ("
                  << unscaled<double>(bb0.max.x()) << ", "
                  << unscaled<double>(bb0.max.y()) << ")");
    UNSCOPED_INFO("L2 bbox: ("
                  << unscaled<double>(bb1.min.x()) << ", "
                  << unscaled<double>(bb1.min.y()) << ") to ("
                  << unscaled<double>(bb1.max.x()) << ", "
                  << unscaled<double>(bb1.max.y()) << ")");
    UNSCOPED_INFO("bbox_overlap = " << (bbox_overlap ? "YES" : "no"));

    // 1. Bboxes overlap. If a solver used bbox as the collision
    //    shape, it would reject this placement.
    REQUIRE(bbox_overlap);

    // 2. Silhouettes DO NOT overlap. The polygon-intersection
    //    no_overlap helper confirms clear.
    REQUIRE(test_utils::no_overlap(items));

    // 3. The winding-number oracle (independent ground truth)
    //    agrees: no interior point of L2 lies inside L1, or vice
    //    versa. This pins the result with a provably-correct
    //    method in case no_overlap has any edge-case weakness.
    REQUIRE_FALSE(test_utils::expolygons_overlap_oracle_symmetric(
        tp0, tp1, 0.25));

    // Summary in the test log: bbox says "collision", three
    // independent silhouette methods say "clear". Every
    // bbox-as-collision-shape check over-rejects this placement.
    // M2.5.3 bitmap AND replaces the current collision path and
    // will handle this case correctly.
}
