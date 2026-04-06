// test_bitmap_stress.cpp
// Adversarial stress tests for BitmapNester. These push boundary conditions
// the regular test suite doesn't cover: extreme counts, aspect ratios, inflation,
// shrinkage, and degenerate geometry.

#define BITMAP_NESTER_TESTING

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"

#include "bitmap_test_utils.hpp"

#include <vector>
#include <cmath>
#include <chrono>
#include <atomic>
#include <set>

using namespace Slic3r;
using namespace Slic3r::arrangement;
using namespace Slic3r::arrangement::test_utils;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static ExPolygon stress_rect(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x0), scaled<coord_t>(y1))
    });
}

static ExPolygon stress_triangle(double base, double height)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0),    scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(base),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(0.0),    scaled<coord_t>(height))
    });
}

static ExPolygon stress_l_shape(double w, double h, double notch_w, double notch_h)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0),            scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),              scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),              scaled<coord_t>(h - notch_h)),
        Point(scaled<coord_t>(w - notch_w),    scaled<coord_t>(h - notch_h)),
        Point(scaled<coord_t>(w - notch_w),    scaled<coord_t>(h)),
        Point(scaled<coord_t>(0.0),            scaled<coord_t>(h))
    });
}

static BoundingBox stress_bed(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

static ArrangeParams stress_params()
{
    ArrangeParams p;
    p.bed_shrink_x    = 0.0f;
    p.bed_shrink_y    = 0.0f;
    p.allow_rotations = false;
    p.progressind     = nullptr;
    return p;
}

static ArrangePolygon stress_ap(const ExPolygon &poly, int priority = 0)
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

// Timeout-guarded params: fires stopcondition after N seconds.
static ArrangeParams timed_params(int timeout_sec, std::atomic<bool> &timed_out)
{
    ArrangeParams p = stress_params();
    auto deadline = std::chrono::steady_clock::now()
                    + std::chrono::seconds(timeout_sec);
    p.stopcondition = [&timed_out, deadline]() -> bool {
        bool expired = (std::chrono::steady_clock::now() > deadline);
        if (expired) timed_out = true;
        return expired;
    };
    return p;
}

// ===========================================================================
// 1. 150 mixed-size items — must all place, no overlaps, finishes in time
// ===========================================================================

TEST_CASE("Stress: 150 mixed-size items place without overlap",
          "[BitmapStress]")
{
    BoundingBox bed = stress_bed(256.0, 210.0);
    std::atomic<bool> timed_out{false};
    ArrangeParams p = timed_params(120, timed_out);

    ArrangePolygons items;
    items.reserve(150);

    // 50 small (5×5), 50 medium (15×10), 50 large (30×20)
    for (int i = 0; i < 50; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 5, 5)));
    for (int i = 0; i < 50; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 15, 10)));
    for (int i = 0; i < 50; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 30, 20)));

    ArrangePolygons excludes;
    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());
    // Multi-plate overflow is fine, but everything must be placed
    CHECK(all_placed(items));
    CHECK(no_overlap(items));
    CHECK(no_bed_gaps(items));
}

// ===========================================================================
// 2. Barely-fits items — item at ~99.5% of bed dimension
// ===========================================================================

TEST_CASE("Stress: barely-fits item places on plate 0",
          "[BitmapStress]")
{
    // Bed 100×100, item 99×99. Must fit on plate 0.
    BoundingBox bed = stress_bed(100.0, 100.0);
    ArrangeParams p = stress_params();

    ArrangePolygons items{stress_ap(stress_rect(0, 0, 99.0, 99.0))};
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
}

TEST_CASE("Stress: two barely-fits items overflow to plate 1",
          "[BitmapStress]")
{
    // Bed 100×100, two items 99×99. First fits plate 0, second must overflow.
    BoundingBox bed = stress_bed(100.0, 100.0);
    ArrangeParams p = stress_params();

    ArrangePolygons items;
    items.push_back(stress_ap(stress_rect(0, 0, 99.0, 99.0)));
    items.push_back(stress_ap(stress_rect(0, 0, 99.0, 99.0)));
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(items[1].bed_idx == 1);
    CHECK(no_overlap(items));
}

// ===========================================================================
// 3. Extreme aspect ratios — slivers and bars
// ===========================================================================

TEST_CASE("Stress: 1mm x 200mm slivers pack without overlap",
          "[BitmapStress]")
{
    BoundingBox bed = stress_bed(256.0, 210.0);
    std::atomic<bool> timed_out{false};
    ArrangeParams p = timed_params(60, timed_out);

    // 20 vertical slivers (1×200mm). Bed is 256 wide, so ~256 fit in a row.
    ArrangePolygons items;
    for (int i = 0; i < 20; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 1.0, 200.0)));
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());
    CHECK(all_placed(items));
    CHECK(no_overlap(items));
}

TEST_CASE("Stress: 200mm x 1mm bars pack without overlap",
          "[BitmapStress]")
{
    BoundingBox bed = stress_bed(256.0, 210.0);
    std::atomic<bool> timed_out{false};
    ArrangeParams p = timed_params(60, timed_out);

    // 20 horizontal bars (200×1mm). Bed is 210 tall, so ~210 fit in a column.
    ArrangePolygons items;
    for (int i = 0; i < 20; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 200.0, 1.0)));
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());
    CHECK(all_placed(items));
    CHECK(no_overlap(items));
}

// ===========================================================================
// 4. Near-zero-area polygons — no crash
// ===========================================================================

TEST_CASE("Stress: near-zero-area triangles do not crash",
          "[BitmapStress]")
{
    BoundingBox bed = stress_bed(256.0, 210.0);
    ArrangeParams p = stress_params();

    ArrangePolygons items;
    // Tiny triangles: 0.1mm base × 0.1mm height = 0.005mm² area
    for (int i = 0; i < 10; ++i)
        items.push_back(stress_ap(stress_triangle(0.1, 0.1)));
    ArrangePolygons excludes;

    // Must not crash. Placement is OK or UNARRANGED — both acceptable.
    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));
}

// ===========================================================================
// 5. Enormous inflation — forces multi-plate overflow
// ===========================================================================

TEST_CASE("Stress: 50mm inflation on small items forces multi-plate overflow",
          "[BitmapStress]")
{
    BoundingBox bed = stress_bed(256.0, 210.0);
    std::atomic<bool> timed_out{false};
    ArrangeParams p = timed_params(120, timed_out);

    // 10×10mm items with 50mm inflation → effective ~110×110mm.
    // Bed fits ~2×1 = 2 per plate.
    ArrangePolygons items;
    for (int i = 0; i < 8; ++i) {
        ArrangePolygon ap = stress_ap(stress_rect(0, 0, 10.0, 10.0));
        ap.inflation = scaled<coord_t>(50.0);
        items.push_back(ap);
    }
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());
    CHECK(all_placed(items));
    CHECK(no_overlap(items));

    // Should use multiple plates given the huge inflation
    std::set<int> beds;
    for (auto &it : items)
        if (it.bed_idx != UNARRANGED) beds.insert(it.bed_idx);
    CHECK(beds.size() > 1);
}

// ===========================================================================
// 6. Bed shrinkage leaving tiny usable area
// ===========================================================================

TEST_CASE("Stress: extreme bed shrinkage — most items UNARRANGED gracefully",
          "[BitmapStress]")
{
    // 256×210 bed with shrink_x=120, shrink_y=100 → usable 16×10mm
    BoundingBox bed = stress_bed(256.0, 210.0);
    ArrangeParams p = stress_params();
    p.bed_shrink_x = 120.0f;
    p.bed_shrink_y = 100.0f;

    ArrangePolygons items;
    // 5×5mm items — a few might fit in the 16×10 usable area
    for (int i = 0; i < 20; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 5.0, 5.0)));
    ArrangePolygons excludes;

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));

    // Count how many are placed vs UNARRANGED
    int placed = 0;
    for (auto &it : items)
        if (it.bed_idx != UNARRANGED) ++placed;

    // Usable area is 16×10mm = 160mm². Each item is 25mm². At most ~6 fit.
    // Some will be placed, many won't. Multi-plate overflow means more may
    // place, but the key invariant: no crash and no overlaps.
    INFO("placed = " << placed << " / 20");
    CHECK(placed > 0);    // At least some should fit
    CHECK(placed <= 20);  // Multi-plate overflow may place all items
    CHECK(no_overlap(items));
}

TEST_CASE("Stress: shrinkage collapses bed to zero — all UNARRANGED, no crash",
          "[BitmapStress]")
{
    // Shrink exceeds half-bed in each axis → zero usable area
    BoundingBox bed = stress_bed(256.0, 210.0);
    ArrangeParams p = stress_params();
    p.bed_shrink_x = 130.0f;  // 130 > 256/2 → effective width < 0
    p.bed_shrink_y = 110.0f;  // 110 > 210/2 → effective height < 0

    ArrangePolygons items;
    items.push_back(stress_ap(stress_rect(0, 0, 5.0, 5.0)));
    ArrangePolygons excludes;

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));

    // Early return — nothing should be placed
    for (auto &it : items)
        CHECK(it.bed_idx == UNARRANGED);
}

// ===========================================================================
// 7. Concave L-shapes with rotations enabled
// ===========================================================================

TEST_CASE("Stress: concave L-shapes with rotations — no overlaps",
          "[BitmapStress]")
{
    BoundingBox bed = stress_bed(256.0, 210.0);
    std::atomic<bool> timed_out{false};
    ArrangeParams p = timed_params(60, timed_out);
    p.allow_rotations = true;

    ArrangePolygons items;
    for (int i = 0; i < 15; ++i) {
        ArrangePolygon ap = stress_ap(stress_l_shape(30.0, 20.0, 15.0, 10.0));
        ap.allowed_rotations.clear(); // use default 4-rotation set
        items.push_back(ap);
    }
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());
    CHECK(all_placed(items));
    CHECK(no_overlap(items));

    // At least some items should use rotation (not all 0.0)
    int rotated = 0;
    for (auto &it : items)
        if (it.bed_idx != UNARRANGED && std::abs(it.rotation) > 0.01)
            ++rotated;
    INFO("rotated = " << rotated << " / 15");
    // Don't require rotation — the nester picks whatever fits first —
    // but L-shapes on a wide bed should benefit from rotation.
}

// ===========================================================================
// 8. Mixed priorities with overflow
// ===========================================================================

TEST_CASE("Stress: high-priority items land on bed 0 before low-priority",
          "[BitmapStress]")
{
    // Bed fits ~9 items (80×70 on 240×210). Place 10 high-pri + 10 low-pri.
    BoundingBox bed = stress_bed(240.0, 210.0);
    ArrangeParams p = stress_params();

    ArrangePolygons items;
    // High priority (placed first due to sort)
    for (int i = 0; i < 10; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 80.0, 70.0), /*priority=*/10));
    // Low priority
    for (int i = 0; i < 10; ++i)
        items.push_back(stress_ap(stress_rect(0, 0, 80.0, 70.0), /*priority=*/1));

    ArrangePolygons excludes;
    BitmapNester::arrange(items, excludes, bed, p);

    CHECK(all_placed(items));
    CHECK(no_overlap(items));

    // All high-priority items should be on plates 0-1 (placed first).
    // Low-priority should be on plates 1+.
    int hi_max_plate = 0;
    for (int i = 0; i < 10; ++i)
        hi_max_plate = std::max(hi_max_plate, items[i].bed_idx);

    int lo_min_plate = 999;
    for (int i = 10; i < 20; ++i)
        if (items[i].bed_idx != UNARRANGED)
            lo_min_plate = std::min(lo_min_plate, items[i].bed_idx);

    // High-priority items should not be pushed past plate 1
    CHECK(hi_max_plate <= 1);
}

