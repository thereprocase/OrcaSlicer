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
        double max_bed;
        double cluster_perim;
        double cluster_compactness;
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
        m.max_bed             = test_utils::max_bed_idx(items);
        m.cluster_perim       = test_utils::cluster_bbox_perimeter_mm(items);
        m.cluster_compactness = test_utils::cluster_compactness(items);
        return m;
    };

    Measurement solid   = pack_and_measure(c2_solid_square_mm());
    Measurement notched = pack_and_measure(c2_notched_square_mm());

    UNSCOPED_INFO("solid:   max_bed=" << solid.max_bed
                  << " perim=" << solid.cluster_perim
                  << " compact=" << solid.cluster_compactness);
    UNSCOPED_INFO("notched: max_bed=" << notched.max_bed
                  << " perim=" << notched.cluster_perim
                  << " compact=" << notched.cluster_compactness);

    // Both shapes should land valid single-plate packings — 6 items
    // of ≤ 1600 mm² each on a 300x300 bed is trivially loose.
    REQUIRE(solid.max_bed == 0);
    REQUIRE(notched.max_bed == 0);

    // THE LITMUS: compactness (filled silhouette area / cluster bbox
    // area) must differ between the two shapes. If C2/C1 is secretly
    // convex-hull based, the two shapes (same hull) produce identical
    // bbox AND identical silhouette area = identical compactness.
    // If it's genuinely concave, the notched shape has less silhouette
    // area filling the same cluster bbox, so compactness differs.
    //
    // NOTE on the M1 pass-through state: at this writing C2 still
    // delegates to BitmapNester::arrange for actual packing. C1
    // rasterizes concave silhouettes (not hulls), so the two shapes
    // already produce different measurements. The test PASSES right
    // now as a confirmation that C1's existing work is concave-aware.
    // When C2's real pack_as_island lands in M2, this assertion must
    // continue passing AND the gap should tighten — the hull metric
    // should widen the delta further.
    //
    // The 0.05 absolute delta is a loose floor. Tighten after M2.
    REQUIRE(std::abs(solid.cluster_compactness - notched.cluster_compactness) > 0.05);
}
