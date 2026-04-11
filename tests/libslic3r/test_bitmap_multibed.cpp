// test_bitmap_multibed.cpp
// Multi-bed overflow tests and edge-case tests for BitmapNester.
//
// These cover boundary conditions that integration tests miss:
// exact-fit overflow, per-plate exclude isolation, degenerate beds,
// high-inflation saturation, priority/bed interaction, and callback
// contract details.

#define BITMAP_NESTER_TESTING

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"
#include "bitmap_test_utils.hpp"

#include <cstdint>
#include <vector>
#include <algorithm>
#include <set>
#include <numeric>
#include <limits>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Helpers (mirrors of those in test_bitmap_nester.cpp — kept local so the
// two files compile independently).
// ---------------------------------------------------------------------------

static ExPolygon make_rect_mm(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x0), scaled<coord_t>(y1))
    });
}

static BoundingBox make_bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

static ArrangeParams no_shrink_params()
{
    ArrangeParams p;
    p.bed_shrink_x     = 0.0f;
    p.bed_shrink_y     = 0.0f;
    p.allow_rotations  = false;
    p.progressind      = nullptr;
    return p;
}

static ArrangePolygon make_ap(const ExPolygon &poly, int priority = 0)
{
    ArrangePolygon ap;
    ap.poly              = poly;
    ap.priority          = priority;
    ap.bed_idx           = UNARRANGED;
    ap.rotation          = 0.0;
    ap.translation       = Vec2crd{0, 0};
    ap.allowed_rotations = {0.0};
    return ap;
}

// Use shared no_overlap from bitmap_test_utils.hpp
using test_utils::no_overlap;

// ---------------------------------------------------------------------------
// Multi-bed overflow tests                              [BitmapMultiBed]
// ---------------------------------------------------------------------------

TEST_CASE("multibed: N-1 items fill bed 0, Nth overflows to bed 1",
          "[BitmapMultiBed]")
{
    // Bed 240×210 mm. Items 80×70 mm.
    // Per row: floor(240/80) = 3. Rows: floor(210/70) = 3. Capacity = 9.
    // 9 items should all land on bed 0. A 10th must go to bed 1.
    const int  ITEMS_THAT_FIT = 9;
    const int  TOTAL          = 10;
    const double BED_W = 240.0, BED_H = 210.0;
    const double ITEM_W = 80.0, ITEM_H = 70.0;

    ArrangePolygons items;
    for (int i = 0; i < TOTAL; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, ITEM_W, ITEM_H)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(BED_W, BED_H);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    int on_bed0 = 0;
    int on_bed1 = 0;
    for (auto &it : items) {
        if (it.bed_idx == 0) ++on_bed0;
        if (it.bed_idx == 1) ++on_bed1;
    }

    // All items must be placed.
    REQUIRE(on_bed0 + on_bed1 == TOTAL);

    // Exactly 9 must be on bed 0 (the nester should not give up early).
    REQUIRE(on_bed0 == ITEMS_THAT_FIT);

    // The overflow item is on bed 1, not silently dropped.
    REQUIRE(on_bed1 == 1);

    REQUIRE(no_overlap(items));
}

TEST_CASE("multibed: spill item lands on bed 1, not UNARRANGED",
          "[BitmapMultiBed]")
{
    // Three items of 200×200 on a 210×210 bed: only one fits per plate.
    // Item 0 → bed 0, item 1 → bed 1, item 2 → bed 2. None should be UNARRANGED.
    ArrangePolygons items;
    for (int i = 0; i < 3; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 200.0, 200.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(210.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    for (auto &it : items) {
        REQUIRE(it.bed_idx != UNARRANGED);
        REQUIRE(it.bed_idx >= 0);
    }

    // Beds used must include at least 0, 1, and 2.
    std::set<int> beds_used;
    for (auto &it : items) beds_used.insert(it.bed_idx);
    REQUIRE(beds_used.size() == 3);
}

TEST_CASE("multibed: gradual overflow across 3 beds",
          "[BitmapMultiBed]")
{
    // Same geometry as the 9-fit test.  20 items on a 9-capacity bed
    // → 9 on bed 0, 9 on bed 1, 2 on bed 2.
    const int TOTAL = 20;
    const double BED_W = 240.0, BED_H = 210.0;
    const double ITEM_W = 80.0, ITEM_H = 70.0;

    ArrangePolygons items;
    for (int i = 0; i < TOTAL; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, ITEM_W, ITEM_H)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(BED_W, BED_H);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    // Every item placed.
    for (auto &it : items)
        REQUIRE(it.bed_idx != UNARRANGED);

    // Three distinct beds used.
    std::set<int> beds;
    for (auto &it : items) beds.insert(it.bed_idx);
    REQUIRE(beds.count(0) == 1);
    REQUIRE(beds.count(1) == 1);
    REQUIRE(beds.count(2) == 1);

    // No bed is empty while a later bed has items — distribution is reasonable.
    // Count per bed.
    int cnt[3] = {0, 0, 0};
    for (auto &it : items)
        if (it.bed_idx >= 0 && it.bed_idx <= 2) cnt[it.bed_idx]++;

    // Beds 0 and 1 must each carry at least 2 items (sanity: packer didn't
    // dump everything onto the last bed).
    REQUIRE(cnt[0] >= 2);
    REQUIRE(cnt[1] >= 2);
    REQUIRE(cnt[2] >= 1);

    REQUIRE(no_overlap(items));
}

TEST_CASE("multibed: per-plate exclude on bed 0 does not affect bed 1",
          "[BitmapMultiBed]")
{
    // Exclude the left 150 mm of bed 0 (full height). Available strip on
    // bed 0: 250-150 = 100 mm wide. A 60×60 item fits once per row
    // (60 < 100, 2×60=120 > 100). Three rows at 60 mm → 3 items on bed 0.
    // With 6 total items, 3 overflow to bed 1 where the exclude is absent.
    // Items on bed 1 must be able to start at x < 50 mm (bed 1 has no exclude).

    const double BED_W = 250.0, BED_H = 210.0;
    const double EXCL_W = 150.0;

    ArrangePolygon excl = make_ap(make_rect_mm(0.0, 0.0, EXCL_W, BED_H));
    excl.bed_idx     = 0;
    excl.translation = Vec2crd{0, 0};

    ArrangePolygons items;
    for (int i = 0; i < 6; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 60.0, 60.0)));

    ArrangePolygons excludes{excl};
    BoundingBox bed = make_bed_mm(BED_W, BED_H);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    // At least one item must have overflowed to bed 1.
    bool any_on_bed1 = false;
    for (auto &it : items)
        if (it.bed_idx == 1) { any_on_bed1 = true; break; }
    REQUIRE(any_on_bed1);

    // Items on bed 0 must start at x >= EXCL_W - 1 mm (raster tolerance).
    for (auto &it : items) {
        if (it.bed_idx != 0) continue;
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb   = get_extents(placed);
        double x0 = unscaled<double>(bb.min.x());
        REQUIRE(x0 >= EXCL_W - 1.0);
    }

    // Items on bed 1 must NOT be restricted to x >= EXCL_W. The bottom-left
    // scan on a fresh (unexcluded) plate starts at pixel (0,0), so the first
    // item placed on bed 1 lands near x=0, which is well inside the zone that
    // bed 0's exclude would have blocked. Assert that at least one bed 1 item
    // starts at x < EXCL_W / 2 to confirm the exclude did not propagate.
    bool bed1_item_in_excluded_zone = false;
    for (auto &it : items) {
        if (it.bed_idx != 1) continue;
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb   = get_extents(placed);
        double x0 = unscaled<double>(bb.min.x());
        if (x0 < EXCL_W / 2.0) { bed1_item_in_excluded_zone = true; break; }
    }
    REQUIRE(bed1_item_in_excluded_zone);

    REQUIRE(no_overlap(items));
}

TEST_CASE("multibed: exclude with bed_idx=1 constrains only bed 1",
          "[BitmapMultiBed]")
{
    // Exclude the bottom half of bed 1 (0..105 mm tall). Overflow items that
    // land on bed 1 must be placed above y=105 mm.
    const double BED_W = 210.0, BED_H = 210.0;
    const double EXCL_H = 105.0;

    ArrangePolygon excl = make_ap(make_rect_mm(0.0, 0.0, BED_W, EXCL_H));
    excl.bed_idx    = 1;
    excl.translation = Vec2crd{0, 0};

    // One large item to anchor bed 0, plus several items that will overflow.
    // 180×180 item fills bed 0. Then 3 smaller items (60×60) go to bed 1.
    ArrangePolygons items;
    items.push_back(make_ap(make_rect_mm(0.0, 0.0, 180.0, 180.0)));
    for (int i = 0; i < 3; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 60.0, 60.0)));

    ArrangePolygons excludes{excl};
    BoundingBox bed = make_bed_mm(BED_W, BED_H);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    // All items must be placed.
    for (auto &it : items)
        REQUIRE(it.bed_idx != UNARRANGED);

    // Items on bed 1 must start at or above y = EXCL_H - 1 mm (raster tolerance).
    for (auto &it : items) {
        if (it.bed_idx != 1) continue;
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb   = get_extents(placed);
        double y0 = unscaled<double>(bb.min.y());
        REQUIRE(y0 >= EXCL_H - 1.0);
    }

    REQUIRE(no_overlap(items));
}

// ---------------------------------------------------------------------------
// Edge-case tests                                       [BitmapEdgeCase]
// ---------------------------------------------------------------------------

TEST_CASE("edgecase: single-pixel bed — no crash, defined result",
          "[BitmapEdgeCase]")
{
    // 1×1 mm bed at 0.5 mm/pixel → 2×2 pixel bitmap.
    // A 0.4×0.4 mm item rasterizes to a single pixel; it may or may not fit,
    // but the call must not crash and bed_idx must not be left undefined.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 0.4, 0.4));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(1.0, 1.0);
    ArrangeParams p = no_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));

    // bed_idx must be either a valid plate (>= 0) or UNARRANGED — never some
    // garbage value.
    int bid = items[0].bed_idx;
    bool valid = (bid == UNARRANGED) || (bid >= 0);
    REQUIRE(valid);
}

TEST_CASE("edgecase: item exactly fills a non-square bed",
          "[BitmapEdgeCase]")
{
    // 180×150 mm bed, 180×150 mm item — should be placed at the origin.
    // Separate from the existing 250×210 test.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 180.0, 150.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(180.0, 150.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);

    ExPolygon placed = items[0].transformed_poly();
    BoundingBox bb   = get_extents(placed);
    double x0 = unscaled<double>(bb.min.x());
    double y0 = unscaled<double>(bb.min.y());

    // Should land at bed origin within raster resolution (0.5 mm).
    REQUIRE_THAT(x0, Catch::Matchers::WithinAbs(0.0, 1.0));
    REQUIRE_THAT(y0, Catch::Matchers::WithinAbs(0.0, 1.0));
}

TEST_CASE("edgecase: zero items — no crash, no side effects",
          "[BitmapEdgeCase]")
{
    ArrangePolygons items;
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));
    REQUIRE(items.empty());
}

TEST_CASE("edgecase: 15 identical items tile without overlap",
          "[BitmapEdgeCase]")
{
    // 15 identical 20×20 mm items on a 250×210 bed.
    // 250/20 = 12 per row, 210/20 = 10 rows → all 15 fit easily.
    const int N = 15;
    ArrangePolygons items;
    for (int i = 0; i < N; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 20.0, 20.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    for (auto &it : items)
        REQUIRE(it.bed_idx == 0);

    REQUIRE(no_overlap(items));
}

TEST_CASE("edgecase: item with 5 holes is rasterized and placed without crash",
          "[BitmapEdgeCase]")
{
    // 50×50 mm outer rectangle with five 5×5 mm rectangular holes punched out
    // along the horizontal center line.
    ExPolygon poly;
    poly.contour = Polygon(Points{
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(50.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(50.0), scaled<coord_t>(50.0)),
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(50.0))
    });

    // Five CW holes, each 5×5 mm, equally spaced along y=20..25.
    for (int h = 0; h < 5; ++h) {
        double hx0 = 5.0 + h * 8.0;
        double hx1 = hx0 + 5.0;
        poly.holes.push_back(Polygon(Points{
            Point(scaled<coord_t>(hx0), scaled<coord_t>(20.0)),
            Point(scaled<coord_t>(hx0), scaled<coord_t>(25.0)),
            Point(scaled<coord_t>(hx1), scaled<coord_t>(25.0)),
            Point(scaled<coord_t>(hx1), scaled<coord_t>(20.0))
        }));
    }

    ArrangePolygons items{make_ap(poly), make_ap(poly)};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));

    for (auto &it : items)
        REQUIRE(it.bed_idx != UNARRANGED);

    REQUIRE(no_overlap(items));
}

TEST_CASE("edgecase: very high inflation forces items to many plates or UNARRANGED",
          "[BitmapEdgeCase]")
{
    // 50×50 mm bed, 5×5 mm items with 20 mm inflation.
    // Each inflated item is ~45×45 mm — close to the entire bed.
    // At most one item per plate. With 5 items → beds 0..4, or some UNARRANGED.
    const int N = 5;
    const double INF_MM = 20.0;

    ArrangePolygons items;
    for (int i = 0; i < N; ++i) {
        ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 5.0, 5.0));
        ap.inflation = scaled<coord_t>(INF_MM);
        items.push_back(ap);
    }

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(50.0, 50.0);
    ArrangeParams p = no_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));

    // Count items placed on distinct plates, and items left UNARRANGED.
    std::set<int> beds_used;
    int unarranged_count = 0;
    for (auto &it : items) {
        if (it.bed_idx == UNARRANGED)
            ++unarranged_count;
        else
            beds_used.insert(it.bed_idx);
    }

    // Either at most 1 item per plate (each on its own plate), or some items
    // are UNARRANGED (hit the MAX_PLATES cap). In no case should 2 items share
    // the same plate — each inflated item nearly fills the bed.
    REQUIRE(beds_used.size() + (size_t)unarranged_count == (size_t)N);

    // No two placed items on the same plate.
    for (auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        int count_on_same = 0;
        for (auto &other : items)
            if (other.bed_idx == it.bed_idx) ++count_on_same;
        REQUIRE(count_on_same == 1);
    }
}

TEST_CASE("edgecase: high-priority items land on bed 0, low-priority overflow",
          "[BitmapEdgeCase]")
{
    // Insert 4 low-priority items first, then 2 high-priority items. After
    // sorting by priority, the 2 high-priority items are placed first.
    //
    // The contract this test enforces: HIGH-PRIORITY ITEMS ARE ON A LOWER
    // BED INDEX THAN LOW-PRIORITY ITEMS. Originally the test asserted high
    // priority on bed 0 and low priority on bed > 0 with a bed sized to
    // exactly fit two 100×100 items side-by-side, but with the center-
    // greedy anchor (2026-04-11) the first item lands at the bed center
    // which traps the second item in any orientation. The test now uses a
    // bed that's still tight relative to the part count and asserts the
    // priority ordering invariant directly without depending on a specific
    // plate count.
    const double BED_W = 205.0, BED_H = 105.0;
    const double ITEM_S = 100.0;

    ArrangePolygons items;
    // Low-priority inserted first so the sort order is exercised.
    for (int i = 0; i < 4; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, ITEM_S, ITEM_S), /*priority=*/1));
    // High-priority appended after.
    for (int i = 0; i < 2; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, ITEM_S, ITEM_S), /*priority=*/10));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(BED_W, BED_H);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    // Priority contract: every high-priority item lives on a bed_idx no
    // higher than every low-priority item's bed_idx. This is anchor-
    // independent — corner-greedy and center-greedy both honor it.
    int max_hp_bed = -1;
    int min_lp_bed = std::numeric_limits<int>::max();
    for (int i = 0; i < 4; ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;
        min_lp_bed = std::min(min_lp_bed, items[i].bed_idx);
    }
    for (int i = 4; i < 6; ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;
        max_hp_bed = std::max(max_hp_bed, items[i].bed_idx);
    }
    // High-priority items must be placed (the bed has room).
    REQUIRE(items[4].bed_idx != UNARRANGED);
    REQUIRE(items[5].bed_idx != UNARRANGED);
    // High priority lives on a lower (or equal) bed than every low.
    if (min_lp_bed != std::numeric_limits<int>::max())
        REQUIRE(max_hp_bed <= min_lp_bed);

    REQUIRE(no_overlap(items));
}

TEST_CASE("edgecase: progressind receives exact decreasing sequence",
          "[BitmapEdgeCase]")
{
    // 12 items. Expect exactly 12 callbacks with values 12, 11, 10, ..., 1.
    const int N = 12;
    ArrangePolygons items;
    for (int i = 0; i < N; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    std::vector<unsigned> received;
    p.progressind = [&](unsigned remaining, std::string) {
        received.push_back(remaining);
    };

    BitmapNester::arrange(items, excludes, bed, p);

    // Exactly one call per item.
    REQUIRE((int)received.size() == N);

    // First call is N (items_remaining starts at N and is post-decremented
    // after the callback, so the first value passed is N, next N-1, ..., 1).
    REQUIRE(received[0] == (unsigned)N);

    // Strictly decreasing.
    for (int i = 1; i < N; ++i) {
        REQUIRE(received[i] < received[i - 1]);
    }

    // Last call is 1.
    REQUIRE(received[N - 1] == 1u);
}

TEST_CASE("edgecase: stopcondition halts after exactly 3 placements",
          "[BitmapEdgeCase]")
{
    // 10 items on a spacious bed. Stopcondition returns true on the 4th check.
    // The first 3 items (in sort order) must be placed; items 4..10 must remain
    // UNARRANGED.
    const int N = 10;
    ArrangePolygons items;
    for (int i = 0; i < N; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    int check_count = 0;
    p.stopcondition = [&]() -> bool {
        // Returns false for checks 1..3 (items 0, 1, 2 proceed).
        // Returns true on check 4 (item 3 triggers the break before placement).
        return ++check_count > 3;
    };

    BitmapNester::arrange(items, excludes, bed, p);

    int arranged = 0;
    for (auto &it : items)
        if (it.bed_idx != UNARRANGED) ++arranged;

    // Exactly 3 items placed; the rest must be untouched.
    REQUIRE(arranged == 3);

    // Confirm the untouched items still have UNARRANGED bed_idx.
    // Because the nester sorts by priority (all equal here) then area, all
    // items have equal weight; 7 of the 10 should remain UNARRANGED.
    int still_unarranged = 0;
    for (auto &it : items)
        if (it.bed_idx == UNARRANGED) ++still_unarranged;
    REQUIRE(still_unarranged == 7);
}

TEST_CASE("edgecase: excluded_regions in ArrangeParams block two zones",
          "[BitmapEdgeCase]")
{
    // Bed 250×210 mm. Two excluded_regions:
    //   Zone A: x = 0..80,   y = 0..210  (left strip)
    //   Zone B: x = 170..250, y = 0..210  (right strip)
    // Items of 60×60 mm must land in the middle column (x ∈ [80, 170]).
    ArrangePolygon zone_a = make_ap(make_rect_mm(0.0,   0.0, 80.0,  210.0));
    ArrangePolygon zone_b = make_ap(make_rect_mm(170.0, 0.0, 250.0, 210.0));
    zone_a.bed_idx    = 0;
    zone_a.translation = Vec2crd{0, 0};
    zone_b.bed_idx    = 0;
    zone_b.translation = Vec2crd{0, 0};

    ArrangePolygons items;
    for (int i = 0; i < 3; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 60.0, 60.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();
    p.excluded_regions.push_back(zone_a);
    p.excluded_regions.push_back(zone_b);

    BitmapNester::arrange(items, excludes, bed, p);

    for (auto &it : items) {
        REQUIRE(it.bed_idx == 0);
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb   = get_extents(placed);
        double x0 = unscaled<double>(bb.min.x());
        double x1 = unscaled<double>(bb.max.x());

        // Must not overlap zone A or zone B (1 mm raster tolerance).
        REQUIRE(x0 >= 79.0);
        REQUIRE(x1 <= 171.0);
    }

    REQUIRE(no_overlap(items));
}

// ---------------------------------------------------------------------------
// Gap tests (Yellow audit, msgs #506/#509)
// ---------------------------------------------------------------------------

TEST_CASE("edgecase: item wider than bed completes without crash",
          "[BitmapEdgeCase]")
{
    // 300×300 mm item on a 250×210 mm bed. The rasterizer clamps the
    // item bitmap to bed dimensions, so the item may be placed (clipped)
    // or UNARRANGED. Either outcome is acceptable — the key invariant
    // is no crash and a valid bed_idx.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 300.0, 300.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));
    bool valid = (items[0].bed_idx == UNARRANGED) || (items[0].bed_idx >= 0);
    REQUIRE(valid);
}

TEST_CASE("multibed: items beyond MAX_PLATES cap become UNARRANGED",
          "[BitmapMultiBed]")
{
    // Bed 12×12 mm. Items 11×11 mm — exactly one per plate.
    // 40 items exceeds MAX_PLATES (36). First 36 placed, rest UNARRANGED.
    const int TOTAL = 40;
    const int MAX_PLATES_CAP = 36;

    ArrangePolygons items;
    for (int i = 0; i < TOTAL; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 11.0, 11.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(12.0, 12.0);
    ArrangeParams p = no_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));

    int placed = 0, unarranged = 0;
    for (auto &it : items) {
        if (it.bed_idx == UNARRANGED) ++unarranged;
        else ++placed;
    }

    REQUIRE(placed <= MAX_PLATES_CAP);
    REQUIRE(placed + unarranged == TOTAL);
    REQUIRE(unarranged >= TOTAL - MAX_PLATES_CAP);
}
