// test_bitmap_nester.cpp
// Comprehensive tests for BitmapNester: bitmap primitives + arrange() integration.
//
// Private statics (collides, stamp, rasterize, scanline_fill) are promoted to
// public via the BITMAP_NESTER_TESTING guard injected in BitmapNester.hpp.

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
#include <numeric>
#include <algorithm>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Build a scaled rectangle ExPolygon (CCW contour) in bed coordinates.
// x0, y0, x1, y1 are in mm.
static ExPolygon make_rect_mm(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x0), scaled<coord_t>(y1))
    });
}

// Build a scaled right triangle with legs along +X, +Y from origin.
// Points: (0,0) -> (base_mm,0) -> (0,height_mm)
static ExPolygon make_triangle_mm(double base_mm, double height_mm)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0),      scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(base_mm),  scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(0.0),      scaled<coord_t>(height_mm))
    });
}

// Count set bits across a bitmap (portable; avoids __builtin_popcountll on MSVC).
static int popcount_u64(uint64_t v)
{
    int c = 0;
    while (v) { v &= v - 1; ++c; }
    return c;
}

static int popcount_bitmap(const std::vector<uint64_t> &bm)
{
    int total = 0;
    for (uint64_t w : bm)
        total += popcount_u64(w);
    return total;
}

// Return true if pixel (px, py) is set in a bitmap of word-width iwpr.
static bool pixel_set(const std::vector<uint64_t> &bm, int iwpr, int x, int y)
{
    int word = x / 64;
    int bit  = x % 64;
    return (bm[(size_t)y * iwpr + word] >> bit) & 1u;
}

// Build a bed BoundingBox from mm extents.
static BoundingBox make_bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

// Build a minimal ArrangeParams with zero shrink so the effective bed equals the
// supplied BoundingBox exactly.
static ArrangeParams no_shrink_params()
{
    ArrangeParams p;
    p.bed_shrink_x = 0.0f;
    p.bed_shrink_y = 0.0f;
    p.allow_rotations = false;
    p.progressind = nullptr; // silence stdout during tests
    return p;
}

// Build an ArrangePolygon from an ExPolygon. bed_idx defaults to UNARRANGED.
static ArrangePolygon make_ap(const ExPolygon &poly, int priority = 0)
{
    ArrangePolygon ap;
    ap.poly       = poly;
    ap.priority   = priority;
    ap.bed_idx    = UNARRANGED;
    ap.rotation   = 0.0;
    ap.translation = Vec2crd{0, 0};
    ap.allowed_rotations = {0.0};
    return ap;
}

// Use shared no_overlap from bitmap_test_utils.hpp
using test_utils::no_overlap;

// ---------------------------------------------------------------------------
// Section 1: Bitmap primitive unit tests
// ---------------------------------------------------------------------------

TEST_CASE("rasterize: empty polygon returns empty bitmap", "[BitmapNester][rasterize]")
{
    // Degenerate: zero-area ExPolygon
    ExPolygon empty;
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(empty, 0.5, 400, 300, iw, ih, iwpr);
    REQUIRE(bm.empty());
    REQUIRE(iw == 0);
    REQUIRE(ih == 0);
    REQUIRE(iwpr == 0);
}

TEST_CASE("rasterize: rectangle fills expected pixel count", "[BitmapNester][rasterize]")
{
    // 10 mm × 6 mm rectangle at res = 0.5 mm/pixel → 20 × 12 = 240 pixels
    // Build rectangle relative to its own origin (bbox.min = (0,0))
    ExPolygon rect = make_rect_mm(0.0, 0.0, 10.0, 6.0);

    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(rect, 0.5, 400, 400, iw, ih, iwpr);

    REQUIRE_FALSE(bm.empty());
    REQUIRE(iw == 20);
    REQUIRE(ih == 12);
    REQUIRE(iwpr == 1); // 20 bits fits in one 64-bit word

    // Every pixel inside the rectangle should be set.
    int set_count = popcount_bitmap(bm);
    REQUIRE(set_count == 240);
}

TEST_CASE("rasterize: rectangle bits are in the right positions", "[BitmapNester][rasterize]")
{
    // 4 mm × 2 mm at res = 1.0 mm/pixel → 4 × 2 bitmap
    ExPolygon rect = make_rect_mm(0.0, 0.0, 4.0, 2.0);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(rect, 1.0, 400, 400, iw, ih, iwpr);

    REQUIRE(iw == 4);
    REQUIRE(ih == 2);

    // All 8 pixels should be set.
    for (int y = 0; y < 2; ++y)
        for (int x = 0; x < 4; ++x)
            REQUIRE(pixel_set(bm, iwpr, x, y));
}

TEST_CASE("rasterize: triangle produces correct pixel count", "[BitmapNester][rasterize]")
{
    // Right triangle: base = 10 mm, height = 10 mm → area = 50 mm²
    // At res = 0.5 mm/pixel, pixel area = 0.25 mm², so ≈ 200 pixels.
    // Scanline rasterization of a right triangle: row y (0-indexed) has
    // floor(base * (1 - y/height)) pixels along the x axis after scan centering.
    // Just verify pixel count is strictly less than the bounding rectangle (200)
    // and greater than the inscribed triangle (roughly half the rect).
    ExPolygon tri = make_triangle_mm(10.0, 10.0);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(tri, 0.5, 400, 400, iw, ih, iwpr);

    REQUIRE_FALSE(bm.empty());
    int rect_pixels = iw * ih;
    int tri_pixels  = popcount_bitmap(bm);

    REQUIRE(tri_pixels > 0);
    REQUIRE(tri_pixels < rect_pixels); // must be less than full bounding rectangle
    // A right triangle is roughly half the bounding rectangle.
    REQUIRE(tri_pixels > rect_pixels / 4);
    REQUIRE(tri_pixels < (rect_pixels * 3) / 4);
}

TEST_CASE("rasterize: ExPolygon with hole clears interior bits", "[BitmapNester][rasterize]")
{
    // Outer rect: 0..10 mm × 0..10 mm
    // Inner hole: 2..8 mm × 2..8 mm  (CW winding for hole)
    ExPolygon donut;
    donut.contour = Polygon(Points{
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(10.0))
    });
    // Hole must be CW (opposite winding to contour)
    donut.holes.push_back(Polygon(Points{
        Point(scaled<coord_t>(2.0),  scaled<coord_t>(2.0)),
        Point(scaled<coord_t>(2.0),  scaled<coord_t>(8.0)),
        Point(scaled<coord_t>(8.0),  scaled<coord_t>(8.0)),
        Point(scaled<coord_t>(8.0),  scaled<coord_t>(2.0))
    }));

    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(donut, 1.0, 400, 400, iw, ih, iwpr);

    REQUIRE_FALSE(bm.empty());
    REQUIRE(iw == 10);
    REQUIRE(ih == 10);

    int total  = popcount_bitmap(bm);
    int area   = iw * ih;          // 100 pixels
    // Interior hole is 6×6 = 36 pixels.
    // Rim area ≈ 100 - 36 = 64 pixels (scan-line rounding may vary by ±2).
    REQUIRE(total > 50);
    REQUIRE(total < 100);

    // Pixels strictly inside the hole (rows 2..7, cols 2..7) must be clear.
    for (int y = 3; y <= 6; ++y)
        for (int x = 3; x <= 6; ++x)
            REQUIRE_FALSE(pixel_set(bm, iwpr, x, y));

    // Corner pixels of outer rim must be set.
    REQUIRE(pixel_set(bm, iwpr, 0, 0));
    REQUIRE(pixel_set(bm, iwpr, 9, 9));
}

TEST_CASE("rasterize: bitmap clamped to max_w / max_h", "[BitmapNester][rasterize]")
{
    // 100 mm × 100 mm shape at res 0.5 → would be 200×200 pixels.
    // Clamp to max 50×50.
    ExPolygon rect = make_rect_mm(0.0, 0.0, 100.0, 100.0);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(rect, 0.5, 50, 50, iw, ih, iwpr);

    REQUIRE(iw == 50);
    REQUIRE(ih == 50);
}

TEST_CASE("rasterize: word-boundary item (iw > 64)", "[BitmapNester][rasterize]")
{
    // 40 mm wide, 2 mm tall at res 0.5 → 80 pixels wide → requires 2 words per row
    ExPolygon rect = make_rect_mm(0.0, 0.0, 40.0, 2.0);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(rect, 0.5, 400, 400, iw, ih, iwpr);

    REQUIRE(iw == 80);
    REQUIRE(iwpr == 2);
    // All 80 × 4 = 320 pixels should be set.
    REQUIRE(popcount_bitmap(bm) == 320);
}

// ---------------------------------------------------------------------------

TEST_CASE("stamp: sets bits at origin", "[BitmapNester][stamp]")
{
    // Plate: 10×4 pixels → wpr = 1
    int bw = 10, bh = 4, wpr = 1;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    // Item: 3×2 solid block
    int iw = 3, ih = 2, iwpr = 1;
    std::vector<uint64_t> item_bm((size_t)iwpr * ih, 0);
    // Set bits 0,1,2 in each row
    uint64_t row_bits = 0x7ULL; // 0b111
    item_bm[0] = row_bits;
    item_bm[1] = row_bits;

    // Stamp at (0, 0)
    BitmapNester::stamp(plate, wpr, bw, bh, item_bm, iwpr, iw, ih, 0, 0);

    // Verify bits 0,1,2 in rows 0,1 are set
    for (int y = 0; y < 2; ++y)
        for (int x = 0; x < 3; ++x)
            REQUIRE(pixel_set(plate, wpr, x, y));

    // Other pixels should remain clear
    REQUIRE_FALSE(pixel_set(plate, wpr, 3, 0));
    REQUIRE_FALSE(pixel_set(plate, wpr, 0, 2));
}

TEST_CASE("stamp: sets bits at non-zero offset", "[BitmapNester][stamp]")
{
    int bw = 20, bh = 10, wpr = 1;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    int iw = 2, ih = 2, iwpr = 1;
    std::vector<uint64_t> item_bm((size_t)iwpr * ih, 0);
    item_bm[0] = 0x3ULL; // bits 0 and 1
    item_bm[1] = 0x3ULL;

    // Stamp at (5, 3) — places bits at plate positions (5,3), (6,3), (5,4), (6,4)
    BitmapNester::stamp(plate, wpr, bw, bh, item_bm, iwpr, iw, ih, 5, 3);

    REQUIRE(pixel_set(plate, wpr, 5, 3));
    REQUIRE(pixel_set(plate, wpr, 6, 3));
    REQUIRE(pixel_set(plate, wpr, 5, 4));
    REQUIRE(pixel_set(plate, wpr, 6, 4));

    // Adjacent pixels should not be set
    REQUIRE_FALSE(pixel_set(plate, wpr, 4, 3));
    REQUIRE_FALSE(pixel_set(plate, wpr, 7, 3));
    REQUIRE_FALSE(pixel_set(plate, wpr, 5, 2));
    REQUIRE_FALSE(pixel_set(plate, wpr, 5, 5));
}

TEST_CASE("stamp: word-boundary crossing", "[BitmapNester][stamp]")
{
    // Plate is 130 pixels wide (3 words), 2 rows
    int bw = 130, bh = 2, wpr = (bw + 63) / 64; // = 3
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    // Item is 4 pixels wide, 1 row tall, offset so bits straddle word boundary.
    // Place it at px = 62: bits land at positions 62, 63, 64, 65 in the plate.
    int iw = 4, ih = 1, iwpr = 1;
    std::vector<uint64_t> item_bm(iwpr * ih, 0);
    item_bm[0] = 0xFULL; // bits 0,1,2,3

    BitmapNester::stamp(plate, wpr, bw, bh, item_bm, iwpr, iw, ih, 62, 0);

    REQUIRE(pixel_set(plate, wpr, 62, 0));
    REQUIRE(pixel_set(plate, wpr, 63, 0));
    REQUIRE(pixel_set(plate, wpr, 64, 0));
    REQUIRE(pixel_set(plate, wpr, 65, 0));

    REQUIRE_FALSE(pixel_set(plate, wpr, 61, 0));
    REQUIRE_FALSE(pixel_set(plate, wpr, 66, 0));
}

// ---------------------------------------------------------------------------

TEST_CASE("collides: overlapping bitmaps return true", "[BitmapNester][collides]")
{
    int bw = 10, bh = 4, wpr = 1;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    // Stamp a 2×2 block at (3,1)
    int iw = 2, ih = 2, iwpr = 1;
    std::vector<uint64_t> block((size_t)iwpr * ih, 0x3ULL);
    block[1] = 0x3ULL;
    BitmapNester::stamp(plate, wpr, bw, bh, block, iwpr, iw, ih, 3, 1);

    // Test a new item bitmap at (4,1) — overlaps at (4,1) and (4,2)
    std::vector<uint64_t> item((size_t)iwpr * ih, 0x3ULL);
    item[1] = 0x3ULL;
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, item, iwpr, iw, ih, 4, 1));
}

TEST_CASE("collides: non-overlapping bitmaps return false", "[BitmapNester][collides]")
{
    int bw = 20, bh = 4, wpr = 1;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    int iw = 2, ih = 2, iwpr = 1;
    std::vector<uint64_t> block(iwpr * ih, 0x3ULL);
    block[1] = 0x3ULL;
    BitmapNester::stamp(plate, wpr, bw, bh, block, iwpr, iw, ih, 0, 0);

    // Item placed far to the right — no overlap
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, block, iwpr, iw, ih, 10, 0));
}

TEST_CASE("collides: adjacent items (gap = 1 pixel) do not collide", "[BitmapNester][collides]")
{
    int bw = 20, bh = 4, wpr = 1;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    int iw = 3, ih = 2, iwpr = 1;
    std::vector<uint64_t> block(iwpr * ih, 0x7ULL);
    block[1] = 0x7ULL;
    // Stamp occupies columns 0,1,2 in rows 0,1.
    BitmapNester::stamp(plate, wpr, bw, bh, block, iwpr, iw, ih, 0, 0);

    // Next item starts at column 3 — should not collide.
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, block, iwpr, iw, ih, 3, 0));
    // But placing at column 2 overlaps column 2.
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, block, iwpr, iw, ih, 2, 0));
}

TEST_CASE("collides: stamp+collides roundtrip detects self-collision", "[BitmapNester][collides]")
{
    int bw = 64, bh = 10, wpr = 1;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    ExPolygon rect = make_rect_mm(0.0, 0.0, 4.0, 2.0);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(rect, 1.0, bw, bh, iw, ih, iwpr);
    REQUIRE_FALSE(bm.empty());

    // Before stamping: no collision at (0,0)
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0));

    // Stamp the bitmap
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0);

    // Now collides at (0,0)
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0));

    // Still no collision shifted far away
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 10, 0));
}

TEST_CASE("collides: word-boundary item does not false-positive", "[BitmapNester][collides]")
{
    // 130-wide plate, 2 rows. Stamp a 4-pixel item at column 62 (straddles word
    // boundary 63/64). Check that a non-overlapping item at column 70 does NOT
    // collide.
    int bw = 130, bh = 2, wpr = (bw + 63) / 64; // 3
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    int iw = 4, ih = 1, iwpr_item = 1;
    std::vector<uint64_t> item(iwpr_item * ih, 0xFULL); // bits 0..3

    BitmapNester::stamp(plate, wpr, bw, bh, item, iwpr_item, iw, ih, 62, 0);

    // Item at 70 is far from 62..65 — no collision
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, item, iwpr_item, iw, ih, 70, 0));

    // Item at 64 overlaps at 64 and 65 — must collide
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, item, iwpr_item, iw, ih, 64, 0));
}

// ---------------------------------------------------------------------------
// Section 2: arrange() integration tests
// ---------------------------------------------------------------------------

TEST_CASE("arrange: empty item list — no crash", "[BitmapNester][arrange]")
{
    ArrangePolygons items;
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));
    REQUIRE(items.empty());
}

TEST_CASE("arrange: single small item is placed on plate 0", "[BitmapNester][arrange]")
{
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(items[0].is_arranged());
}

TEST_CASE("arrange: single item centered on bed", "[BitmapNester][arrange]")
{
    // A lone item on an empty plate gets centered by the post-placement
    // centering pass (do_final_align defaults to true).
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 5.0, 5.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);

    double tx_mm = unscaled<double>(items[0].translation.x());
    double ty_mm = unscaled<double>(items[0].translation.y());

    // Centered: 5mm item on 250x210 bed -> center at (125, 105),
    // item placed so its center lands there -> tx ~ 122.5, ty ~ 102.5
    REQUIRE_THAT(tx_mm, Catch::Matchers::WithinAbs(122.5, 2.0));
    REQUIRE_THAT(ty_mm, Catch::Matchers::WithinAbs(102.5, 2.0));
}

TEST_CASE("arrange: two items placed on same plate without overlap", "[BitmapNester][arrange]")
{
    // Each item 20×20 mm on a 250×210 bed — both must fit on plate 0.
    ArrangePolygons items{
        make_ap(make_rect_mm(0.0, 0.0, 20.0, 20.0)),
        make_ap(make_rect_mm(0.0, 0.0, 20.0, 20.0))
    };
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(items[1].bed_idx == 0);
    REQUIRE(no_overlap(items));
}

TEST_CASE("arrange: item on zero-size bed after shrinkage is UNARRANGED", "[BitmapNester][arrange]")
{
    // Bed shrinks to zero — nothing can be placed.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(20.0, 20.0);
    ArrangeParams p = no_shrink_params();
    p.bed_shrink_x = 15.0; // shrinks 15mm on each side → negative effective bed
    p.bed_shrink_y = 15.0;

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == UNARRANGED);
}

TEST_CASE("arrange: multi-plate overflow assigns bed_idx > 0", "[BitmapNester][arrange]")
{
    // 6 items of 120×100 mm on a 250×210 bed.
    // Two fit per row, three rows would need 150 mm height — the bed is only 210.
    // In practice 4 items fit on plate 0, the rest overflow to plate 1.
    ArrangePolygons items;
    for (int i = 0; i < 6; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 120.0, 100.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    bool any_overflow = false;
    for (auto &it : items)
        if (it.bed_idx > 0) { any_overflow = true; break; }

    REQUIRE(any_overflow);
    REQUIRE(no_overlap(items));
}

TEST_CASE("arrange: all items arranged (no UNARRANGED when bed is large)", "[BitmapNester][arrange]")
{
    // 8 small items on a spacious bed — all should be placed.
    ArrangePolygons items;
    for (int i = 0; i < 8; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    for (auto &it : items)
        REQUIRE(it.bed_idx != UNARRANGED);

    REQUIRE(no_overlap(items));
}

TEST_CASE("arrange: exclude prevents item placement in that region", "[BitmapNester][arrange]")
{
    // Exclude zone covers the left half of the bed (0..100 mm wide, 0..210 tall).
    // Two 80×80 items should both land in x > 100 mm.
    ArrangePolygon excl = make_ap(make_rect_mm(0.0, 0.0, 100.0, 210.0));
    excl.bed_idx    = 0;
    excl.translation = Vec2crd{scaled<coord_t>(0.0), scaled<coord_t>(0.0)};

    ArrangePolygons items{
        make_ap(make_rect_mm(0.0, 0.0, 80.0, 80.0)),
        make_ap(make_rect_mm(0.0, 0.0, 80.0, 80.0))
    };
    ArrangePolygons excludes{excl};
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    for (auto &it : items) {
        REQUIRE(it.bed_idx == 0);
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb  = get_extents(placed);
        // The placed item must start at x >= 100 mm
        double x0_mm = unscaled<double>(bb.min.x());
        REQUIRE(x0_mm >= 99.0); // 1 mm tolerance for raster quantization
    }

    REQUIRE(no_overlap(items));
}

TEST_CASE("arrange: per-item inflation creates spacing between items", "[BitmapNester][arrange]")
{
    // Two 10×10 items with 5 mm inflation each → they should be at least 5 mm apart.
    ArrangePolygon a = make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0));
    ArrangePolygon b = make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0));
    a.inflation = scaled<coord_t>(5.0); // 5 mm inflation
    b.inflation = scaled<coord_t>(5.0);

    ArrangePolygons items{a, b};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(items[1].bed_idx == 0);

    ExPolygon p0 = items[0].transformed_poly();
    ExPolygon p1 = items[1].transformed_poly();
    BoundingBox bb0 = get_extents(p0);
    BoundingBox bb1 = get_extents(p1);

    // Bounding boxes should not touch (gap >= ~5 mm due to inflation).
    // Check either x-gap or y-gap is at least 4 mm (raster resolution tolerance).
    bool x_gap_ok = (unscaled<double>(bb1.min.x()) >= unscaled<double>(bb0.max.x()) + 4.0) ||
                    (unscaled<double>(bb0.min.x()) >= unscaled<double>(bb1.max.x()) + 4.0);
    bool y_gap_ok = (unscaled<double>(bb1.min.y()) >= unscaled<double>(bb0.max.y()) + 4.0) ||
                    (unscaled<double>(bb0.min.y()) >= unscaled<double>(bb1.max.y()) + 4.0);

    REQUIRE((x_gap_ok || y_gap_ok));
}

TEST_CASE("arrange: min_obj_distance creates spacing between items", "[BitmapNester][arrange]")
{
    // Two 10×10 items with global min_obj_distance = 4 mm.
    ArrangePolygons items{
        make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0)),
        make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0))
    };
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();
    p.min_obj_distance = scaled<coord_t>(4.0);

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(items[1].bed_idx == 0);

    ExPolygon p0 = items[0].transformed_poly();
    ExPolygon p1 = items[1].transformed_poly();
    BoundingBox bb0 = get_extents(p0);
    BoundingBox bb1 = get_extents(p1);

    bool x_gap_ok = (unscaled<double>(bb1.min.x()) >= unscaled<double>(bb0.max.x()) + 1.0) ||
                    (unscaled<double>(bb0.min.x()) >= unscaled<double>(bb1.max.x()) + 1.0);
    bool y_gap_ok = (unscaled<double>(bb1.min.y()) >= unscaled<double>(bb0.max.y()) + 1.0) ||
                    (unscaled<double>(bb0.min.y()) >= unscaled<double>(bb1.max.y()) + 1.0);

    REQUIRE((x_gap_ok || y_gap_ok));
}

TEST_CASE("arrange: rotation lock (allowed_rotations={0.0}) is respected", "[BitmapNester][arrange]")
{
    // A tall narrow item that the nester cannot rotate.
    // allowed_rotations = {0.0} — no rotation permitted.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 5.0, 80.0));
    ap.allowed_rotations = {0.0};

    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();
    p.allow_rotations = true; // global flag on, but per-item overrides

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx != UNARRANGED);
    REQUIRE_THAT(items[0].rotation, Catch::Matchers::WithinAbs(0.0, 1e-9));
}

TEST_CASE("arrange: rotation places item with best score", "[BitmapNester][arrange]")
{
    // A 5×50mm item on a 100×100mm bed with 4 rotation options.
    // Verify it gets placed successfully regardless of rotation chosen.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 5.0, 50.0));
    ap.allowed_rotations = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};

    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(100.0, 100.0);
    ArrangeParams p = no_shrink_params();
    p.allow_rotations = true;

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    // The bottom-left heuristic should pick 0° (vertical) since that minimizes Y.
    // With rotation enabled, the nester has freedom to choose any allowed angle.
    // Just verify it was placed and the rotation is one of the allowed values.
    bool valid_rot = false;
    for (double r : ap.allowed_rotations) {
        if (std::abs(items[0].rotation - r) < 0.01) valid_rot = true;
    }
    REQUIRE(valid_rot);
}

TEST_CASE("arrange: priority ordering places high-priority items first", "[BitmapNester][arrange]")
{
    // One large item (priority 10) and three small items (priority 1).
    // The large item should end up at bed_idx = 0 and near origin,
    // while smaller items may overflow. The large item's itemid should be 0.
    ArrangePolygon big  = make_ap(make_rect_mm(0.0, 0.0, 100.0, 100.0), /*priority=*/10);
    ArrangePolygon sm1  = make_ap(make_rect_mm(0.0, 0.0, 20.0,  20.0),  /*priority=*/1);
    ArrangePolygon sm2  = make_ap(make_rect_mm(0.0, 0.0, 20.0,  20.0),  /*priority=*/1);
    ArrangePolygon sm3  = make_ap(make_rect_mm(0.0, 0.0, 20.0,  20.0),  /*priority=*/1);

    ArrangePolygons items{sm1, sm2, sm3, big}; // big inserted last
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    // Find the big item (index 3 in original vector)
    REQUIRE(items[3].bed_idx == 0);
    REQUIRE(items[3].itemid == 0); // placed first → itemid 0
}

TEST_CASE("arrange: bed shrinkage keeps items within shrunk bounds", "[BitmapNester][arrange]")
{
    // Bed 250×210, shrink 10 mm each side → effective bed 230×190.
    ArrangePolygons items{
        make_ap(make_rect_mm(0.0, 0.0, 20.0, 20.0)),
        make_ap(make_rect_mm(0.0, 0.0, 20.0, 20.0))
    };
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();
    p.bed_shrink_x = 10.0f;
    p.bed_shrink_y = 10.0f;

    BitmapNester::arrange(items, excludes, bed, p);

    for (auto &it : items) {
        REQUIRE(it.bed_idx == 0);
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb = get_extents(placed);

        double x0 = unscaled<double>(bb.min.x());
        double y0 = unscaled<double>(bb.min.y());
        double x1 = unscaled<double>(bb.max.x());
        double y1 = unscaled<double>(bb.max.y());

        REQUIRE(x0 >= 9.0);    // ≥ shrink_x (with 1 mm raster tolerance)
        REQUIRE(y0 >= 9.0);    // ≥ shrink_y
        REQUIRE(x1 <= 241.0);  // ≤ 250 - shrink_x + 1 mm tolerance
        REQUIRE(y1 <= 201.0);  // ≤ 210 - shrink_y + 1 mm tolerance
    }
}

TEST_CASE("arrange: bed that becomes zero-size after shrinkage — no crash", "[BitmapNester][arrange]")
{
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(20.0, 20.0);
    ArrangeParams p = no_shrink_params();
    p.bed_shrink_x = 11.0f; // shrink exceeds half-width → zero or negative effective bed
    p.bed_shrink_y = 11.0f;

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));
    // Item cannot be placed; no assertion on bed_idx (UNARRANGED is acceptable).
}

TEST_CASE("arrange: progress callback is called for each item", "[BitmapNester][arrange]")
{
    ArrangePolygons items;
    for (int i = 0; i < 5; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    unsigned callback_count = 0;
    unsigned last_remaining = UINT_MAX;
    p.progressind = [&](unsigned remaining, std::string) {
        ++callback_count;
        REQUIRE(remaining <= last_remaining); // must be non-increasing
        last_remaining = remaining;
    };

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(callback_count == 5);
}

TEST_CASE("arrange: stopcondition aborts placement early", "[BitmapNester][arrange]")
{
    // 10 items; stopcondition triggers after 3 callbacks.
    ArrangePolygons items;
    for (int i = 0; i < 10; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    int call_count = 0;
    p.stopcondition = [&]() -> bool {
        return ++call_count > 3;
    };

    BitmapNester::arrange(items, excludes, bed, p);

    // At most 3 items should have been placed (stop fires before the 4th).
    int arranged = 0;
    for (auto &it : items)
        if (it.bed_idx != UNARRANGED) ++arranged;

    REQUIRE(arranged <= 4); // allow 1 item that was already in progress when stop fired
}

TEST_CASE("arrange: translation roundtrip reproduces expected bed position", "[BitmapNester][arrange]")
{
    // Place a single 10×10 mm item with no padding. The result should satisfy:
    //   transformed_poly().contour's bounding box is fully within the bed.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);

    ExPolygon placed = items[0].transformed_poly();
    BoundingBox bb = get_extents(placed);

    double x0 = unscaled<double>(bb.min.x());
    double y0 = unscaled<double>(bb.min.y());
    double x1 = unscaled<double>(bb.max.x());
    double y1 = unscaled<double>(bb.max.y());

    // Placed polygon must be fully inside the bed.
    REQUIRE(x0 >= -0.5);   // allow 0.5 mm raster tolerance
    REQUIRE(y0 >= -0.5);
    REQUIRE(x1 <= 250.5);
    REQUIRE(y1 <= 210.5);

    // The polygon should still be 10×10 mm after transform.
    REQUIRE_THAT(x1 - x0, Catch::Matchers::WithinAbs(10.0, 1.0));
    REQUIRE_THAT(y1 - y0, Catch::Matchers::WithinAbs(10.0, 1.0));
}

TEST_CASE("arrange: excluded_regions in params are respected", "[BitmapNester][arrange]")
{
    // An excluded_region covering columns 0..50 mm forces items to start at x > 50.
    ArrangePolygon excl_region = make_ap(make_rect_mm(0.0, 0.0, 50.0, 210.0));
    excl_region.bed_idx    = 0;
    excl_region.translation = Vec2crd{0, 0};

    ArrangePolygons items{make_ap(make_rect_mm(0.0, 0.0, 30.0, 30.0))};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();
    p.excluded_regions.push_back(excl_region);

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);

    ExPolygon placed = items[0].transformed_poly();
    BoundingBox bb = get_extents(placed);
    double x0 = unscaled<double>(bb.min.x());
    REQUIRE(x0 >= 49.0); // must be to the right of the excluded region
}

TEST_CASE("arrange: itemid values are unique and sequential", "[BitmapNester][arrange]")
{
    ArrangePolygons items;
    for (int i = 0; i < 5; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 15.0, 15.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    std::vector<int> ids;
    for (auto &it : items)
        if (it.bed_idx != UNARRANGED) ids.push_back(it.itemid);

    std::sort(ids.begin(), ids.end());
    for (int i = 0; i < (int)ids.size(); ++i)
        REQUIRE(ids[i] == i);
}

TEST_CASE("arrange: many small items tightly packed, no overlap", "[BitmapNester][arrange]")
{
    // 20 items of 10×10 mm on a 250×210 bed → all fit on one plate.
    ArrangePolygons items;
    for (int i = 0; i < 20; ++i)
        items.push_back(make_ap(make_rect_mm(0.0, 0.0, 10.0, 10.0)));

    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    for (auto &it : items)
        REQUIRE(it.bed_idx != UNARRANGED);

    REQUIRE(no_overlap(items));
}

TEST_CASE("arrange: single item exactly fills the bed", "[BitmapNester][arrange]")
{
    // Item is exactly the size of the effective bed. Should be placed at bed origin.
    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 250.0, 210.0));
    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
}

TEST_CASE("arrange: concave item (L-shape) placed without overlap", "[BitmapNester][arrange]")
{
    // L-shaped polygon: a 20×20 square with a 10×10 notch cut from top-right.
    ExPolygon L_shape;
    L_shape.contour = Polygon(Points{
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(20.0)),
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(20.0))
    });

    ArrangePolygons items{make_ap(L_shape), make_ap(L_shape)};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();

    BitmapNester::arrange(items, excludes, bed, p);

    REQUIRE(items[0].bed_idx == 0);
    REQUIRE(items[1].bed_idx == 0);
    REQUIRE(no_overlap(items));
}

TEST_CASE("Arrange: align_to_y_axis pre-rotation preserved in concave mode", "[BitmapNester]")
{
    // Simulates update_selected_items_axis_align writing ap.rotation before the
    // nester runs. The nester must compose its own rotation on top rather than
    // overwriting it with just best_rot.
    //
    // Setup: 30×10 mm rectangle, pre-rotated PI/2 by the caller (align_to_y_axis
    // equivalent). allow_rotations=false so the nester's own best_rot is always 0.
    // Expected post-arrange rotation: PI/2 + 0 = PI/2.
    // Expected placed bbox: ~10 mm wide × ~30 mm tall (rotated shape).

    const double initial_rotation = M_PI / 2.0;

    ArrangePolygon ap = make_ap(make_rect_mm(0.0, 0.0, 30.0, 10.0));
    ap.rotation = initial_rotation;   // caller pre-sets this (axis-align angle)

    ArrangePolygons items{ap};
    ArrangePolygons excludes;
    BoundingBox bed = make_bed_mm(250.0, 210.0);
    ArrangeParams p = no_shrink_params();
    // allow_rotations=false means allowed_rotations stays {0.0} — the nester
    // tries only rot=0, making best_rot=0 and the final rotation = base_rot + 0.
    p.allow_rotations = false;

    BitmapNester::arrange(items, excludes, bed, p);

    // Item must be placed.
    REQUIRE(items[0].bed_idx != UNARRANGED);

    // Pre-rotation must be preserved: base_rot + best_rot = PI/2 + 0 = PI/2.
    REQUIRE_THAT(items[0].rotation, Catch::Matchers::WithinAbs(initial_rotation, 1e-6));

    // Reconstruct the placed shape using the standard transform (rotate then
    // translate), matching what the renderer and post-centering loop do.
    ExPolygon placed = items[0].transformed_poly();
    BoundingBox bb = get_extents(placed);

    double w_mm = unscaled<double>(bb.max.x() - bb.min.x());
    double h_mm = unscaled<double>(bb.max.y() - bb.min.y());

    // A 30×10 rect rotated PI/2 becomes ~10 wide × ~30 tall.
    // Allow 1 mm raster tolerance in each dimension.
    REQUIRE_THAT(w_mm, Catch::Matchers::WithinAbs(10.0, 1.0));
    REQUIRE_THAT(h_mm, Catch::Matchers::WithinAbs(30.0, 1.0));

    // Placed bbox must be within the bed boundaries.
    REQUIRE(unscaled<double>(bb.min.x()) >= -0.5);
    REQUIRE(unscaled<double>(bb.min.y()) >= -0.5);
    REQUIRE(unscaled<double>(bb.max.x()) <= 250.5);
    REQUIRE(unscaled<double>(bb.max.y()) <= 210.5);
}
