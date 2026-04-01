// test_bitmap_determinism.cpp
// Determinism gates and known-good regression fixtures for BitmapNester.
//
// These tests catch uninitialized memory, iterator order assumptions, and
// platform differences. They protect against heuristic changes silently
// degrading arrangement quality.

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
#include <cmath>

using namespace Slic3r;
using namespace Slic3r::arrangement;
using namespace Slic3r::arrangement::test_utils;

// ---------------------------------------------------------------------------
// Local helpers (not duplicated from test_bitmap_nester.cpp)
// ---------------------------------------------------------------------------

static ExPolygon rect_mm(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x0), scaled<coord_t>(y1))
    });
}

static ExPolygon triangle_mm(double base, double height)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0),    scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(base),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(0.0),    scaled<coord_t>(height))
    });
}

// Build an L-shape: outer rectangle with a rectangular notch cut from the
// top-right corner. w/h are the overall dimensions; notch_w/notch_h are the
// size of the removed corner.
static ExPolygon l_shape_mm(double w, double h, double notch_w, double notch_h)
{
    // CCW contour tracing the L:
    //   (0,0) -> (w,0) -> (w, h-notch_h) -> (w-notch_w, h-notch_h) -> (w-notch_w, h) -> (0,h)
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0),            scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),              scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),              scaled<coord_t>(h - notch_h)),
        Point(scaled<coord_t>(w - notch_w),    scaled<coord_t>(h - notch_h)),
        Point(scaled<coord_t>(w - notch_w),    scaled<coord_t>(h)),
        Point(scaled<coord_t>(0.0),            scaled<coord_t>(h))
    });
}

static BoundingBox bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

static ArrangeParams zero_shrink_params()
{
    ArrangeParams p;
    p.bed_shrink_x    = 0.0f;
    p.bed_shrink_y    = 0.0f;
    p.allow_rotations = false;
    p.progressind     = nullptr;
    return p;
}

static ArrangePolygon make_item(const ExPolygon &poly, int priority = 0)
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

// Run the nester and return a snapshot of (translation, rotation, bed_idx)
// for every item, so we can compare two runs.
struct PlacementRecord {
    Vec2crd translation;
    double  rotation;
    int     bed_idx;
};

static std::vector<PlacementRecord> run_and_snapshot(
    ArrangePolygons items, // intentional copy
    const ArrangePolygons &excludes,
    const BoundingBox &bed,
    const ArrangeParams &params)
{
    BitmapNester::arrange(items, excludes, bed, params);
    std::vector<PlacementRecord> out;
    out.reserve(items.size());
    for (auto &it : items)
        out.push_back({it.translation, it.rotation, it.bed_idx});
    return out;
}

// Assert two snapshots are identical.
static void require_same(const std::vector<PlacementRecord> &a,
                          const std::vector<PlacementRecord> &b)
{
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        INFO("Item index: " << i);
        REQUIRE(a[i].bed_idx == b[i].bed_idx);
        REQUIRE(a[i].translation.x() == b[i].translation.x());
        REQUIRE(a[i].translation.y() == b[i].translation.y());
        REQUIRE(a[i].rotation == b[i].rotation);
    }
}

// ---------------------------------------------------------------------------
// Section 1: Determinism gates
// ---------------------------------------------------------------------------

TEST_CASE("Determinism: 5 identical rectangles", "[BitmapDeterminism]")
{
    ArrangePolygons items;
    for (int i = 0; i < 5; ++i)
        items.push_back(make_item(rect_mm(0, 0, 30, 20)));

    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    auto snap1 = run_and_snapshot(items, excludes, bed, p);
    auto snap2 = run_and_snapshot(items, excludes, bed, p);

    require_same(snap1, snap2);
}

TEST_CASE("Determinism: 10 mixed-size rectangles", "[BitmapDeterminism]")
{
    ArrangePolygons items;
    // 10 rectangles with varied dimensions
    double sizes[][2] = {
        {10, 15}, {20, 30}, {50, 25}, {15, 40}, {35, 10},
        {8, 8},   {60, 20}, {25, 25}, {12, 50}, {40, 35}
    };
    for (int i = 0; i < 10; ++i)
        items.push_back(make_item(rect_mm(0, 0, sizes[i][0], sizes[i][1])));

    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    auto snap1 = run_and_snapshot(items, excludes, bed, p);
    auto snap2 = run_and_snapshot(items, excludes, bed, p);

    require_same(snap1, snap2);
}

TEST_CASE("Determinism: 3 concave L-shapes with different sizes", "[BitmapDeterminism]")
{
    ArrangePolygons items;
    items.push_back(make_item(l_shape_mm(40, 30, 15, 10)));
    items.push_back(make_item(l_shape_mm(60, 50, 25, 20)));
    items.push_back(make_item(l_shape_mm(30, 25, 10, 8)));

    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    auto snap1 = run_and_snapshot(items, excludes, bed, p);
    auto snap2 = run_and_snapshot(items, excludes, bed, p);

    require_same(snap1, snap2);
}

// ---------------------------------------------------------------------------
// Section 2: Known-good regression fixtures
// ---------------------------------------------------------------------------

TEST_CASE("Regression: single item on bed", "[BitmapRegression]")
{
    ArrangePolygons items{make_item(rect_mm(0, 0, 20, 15))};
    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(all_within_bounds(items, bed));
}

TEST_CASE("Regression: tightly packed homogeneous set", "[BitmapRegression]")
{
    // 20 identical 30x30mm squares on a 256x210mm bed
    ArrangePolygons items;
    for (int i = 0; i < 20; ++i)
        items.push_back(make_item(rect_mm(0, 0, 30, 30)));

    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(all_placed(items));
    REQUIRE(no_overlap(items));
    REQUIRE(all_within_bounds(items, bed));
    REQUIRE(no_bed_gaps(items));

    // 20 squares at 30x30 = 18000 mm^2; bed = 53760 mm^2.
    // Should fit on 1 plate (8 across x 7 down). At worst 2 plates.
    int max_bed = 0;
    for (auto &it : items)
        max_bed = std::max(max_bed, it.bed_idx);
    REQUIRE(max_bed <= 1);
}

TEST_CASE("Regression: mixed concave/convex shapes", "[BitmapRegression]")
{
    ArrangePolygons items;
    // 3 rectangles
    items.push_back(make_item(rect_mm(0, 0, 40, 25)));
    items.push_back(make_item(rect_mm(0, 0, 20, 50)));
    items.push_back(make_item(rect_mm(0, 0, 35, 15)));
    // 2 L-shapes
    items.push_back(make_item(l_shape_mm(30, 25, 12, 10)));
    items.push_back(make_item(l_shape_mm(45, 35, 18, 15)));
    // 1 triangle
    items.push_back(make_item(triangle_mm(40, 30)));

    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(all_placed(items));
    REQUIRE(no_overlap(items));
    REQUIRE(all_within_bounds(items, bed));
    REQUIRE(no_bed_gaps(items));
}

TEST_CASE("Regression: one item barely fits bed", "[BitmapRegression]")
{
    // Item is 99% of bed width, should still fit on plate 0
    double bed_w = 256.0;
    double item_w = bed_w * 0.99;
    ArrangePolygons items{make_item(rect_mm(0, 0, item_w, 30))};
    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(bed_w, 210);
    ArrangeParams p = zero_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(all_within_bounds(items, bed));
}

TEST_CASE("Regression: item cannot fit — zero-size bed after shrinkage", "[BitmapRegression]")
{
    ArrangePolygons items{make_item(rect_mm(0, 0, 10, 10))};
    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(20, 20);
    ArrangeParams p = zero_shrink_params();
    // Shrink 15mm on each side: effective bed goes negative
    p.bed_shrink_x = 15.0f;
    p.bed_shrink_y = 15.0f;

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == UNARRANGED);
}

TEST_CASE("Regression: empty bed — zero items", "[BitmapRegression]")
{
    ArrangePolygons items;
    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));
    REQUIRE(items.empty());
}

TEST_CASE("Regression: single large + many small items", "[BitmapRegression]")
{
    ArrangePolygons items;

    // 1 large item filling ~60% of bed width
    items.push_back(make_item(rect_mm(0, 0, 150, 120)));

    // 20 tiny items (10x10 each)
    for (int i = 0; i < 20; ++i)
        items.push_back(make_item(rect_mm(0, 0, 10, 10)));

    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(all_placed(items));
    REQUIRE(no_overlap(items));
    REQUIRE(all_within_bounds(items, bed));
    REQUIRE(no_bed_gaps(items));
}

TEST_CASE("Regression: rotation lock — all items keep rotation 0", "[BitmapRegression]")
{
    ArrangePolygons items;
    for (int i = 0; i < 10; ++i) {
        auto ap = make_item(rect_mm(0, 0, 20 + i * 3, 15 + i * 2));
        ap.allowed_rotations = {0.0};
        items.push_back(ap);
    }

    ArrangePolygons excludes;
    BoundingBox bed = bed_mm(256, 210);
    ArrangeParams p = zero_shrink_params();
    // Even if global rotations allowed, per-item lock should win
    p.allow_rotations = true;

    BitmapNester::arrange(items, excludes, bed, p);

    for (size_t i = 0; i < items.size(); ++i) {
        INFO("Item index: " << i);
        if (items[i].bed_idx != UNARRANGED) {
            REQUIRE(items[i].rotation == 0.0);
        }
    }

    REQUIRE(all_placed(items));
    REQUIRE(no_overlap(items));
}
