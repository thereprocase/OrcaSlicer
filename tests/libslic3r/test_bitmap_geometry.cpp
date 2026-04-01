// test_bitmap_geometry.cpp
// Geometry primitive tests for BitmapNester rasterization and collision layers.
// Every test has a known analytical answer. If these fail, integration tests
// built on top become unfalsifiable.

#define BITMAP_NESTER_TESTING

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"

#include <cstdint>
#include <vector>
#include <cmath>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Helpers (scoped to this TU to avoid ODR collisions with test_bitmap_nester)
// ---------------------------------------------------------------------------

namespace {

ExPolygon rect_mm(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x0), scaled<coord_t>(y1))
    });
}

ExPolygon triangle_mm(double x0, double y0,
                      double x1, double y1,
                      double x2, double y2)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x2), scaled<coord_t>(y2))
    });
}

int popcnt64(uint64_t v)
{
    int c = 0;
    while (v) { v &= v - 1; ++c; }
    return c;
}

int popcnt_bm(const std::vector<uint64_t> &bm)
{
    int total = 0;
    for (uint64_t w : bm) total += popcnt64(w);
    return total;
}

bool px_set(const std::vector<uint64_t> &bm, int iwpr, int x, int y)
{
    return (bm[(size_t)y * iwpr + x / 64] >> (x % 64)) & 1u;
}

} // anonymous namespace

// ===========================================================================
// 1. Scanline rasterizer correctness
// ===========================================================================

TEST_CASE("scanline: unit square at 0.5mm res yields exactly 4 pixels",
          "[BitmapGeometry][scanline]")
{
    // 1mm x 1mm at 0.5 mm/pixel -> 2x2 = 4 pixels.
    // Scanline centers at y=0.25 and y=0.75 (both inside [0,1]).
    // Each scanline spans x pixels at centers 0.25 and 0.75 (both inside [0,1]).
    ExPolygon sq = rect_mm(0.0, 0.0, 1.0, 1.0);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 0.5, 400, 400, iw, ih, iwpr);

    REQUIRE(iw == 2);
    REQUIRE(ih == 2);
    REQUIRE(popcnt_bm(bm) == 4);

    // Verify exact bit positions: all 4 pixels set.
    REQUIRE(px_set(bm, iwpr, 0, 0));
    REQUIRE(px_set(bm, iwpr, 1, 0));
    REQUIRE(px_set(bm, iwpr, 0, 1));
    REQUIRE(px_set(bm, iwpr, 1, 1));
}

TEST_CASE("scanline: right triangle pixel count matches area within tolerance",
          "[BitmapGeometry][scanline]")
{
    // Vertices (0,0), (10,0), (0,10). Area = 50 mm^2.
    // At 1 mm/pixel, pixel area = 1 mm^2 -> ~50 pixels.
    // Boundary quantization gives +-5 tolerance.
    ExPolygon tri = triangle_mm(0, 0, 10, 0, 0, 10);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(tri, 1.0, 400, 400, iw, ih, iwpr);

    REQUIRE_FALSE(bm.empty());
    int count = popcnt_bm(bm);
    REQUIRE(count >= 45);
    REQUIRE(count <= 55);
}

TEST_CASE("scanline: degenerate collinear polygon yields zero pixels",
          "[BitmapGeometry][scanline]")
{
    // Three collinear points form a line, not an area.
    // BoundingBox height is zero -> rasterize returns empty.
    ExPolygon line(Points{
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(5.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(0.0))
    });
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(line, 1.0, 400, 400, iw, ih, iwpr);

    // Zero-height bbox -> empty bitmap.
    REQUIRE(bm.empty());
}

TEST_CASE("scanline: sub-pixel polygon produces 0 or 1 pixels without crash",
          "[BitmapGeometry][scanline]")
{
    // Tiny triangle: 0.1mm x 0.1mm at 1mm resolution.
    // The bbox is 0.1mm on each side. ceil(0.1/1.0) = 1 pixel per axis.
    // Whether the scanline center (at y=0.5*res from bbox.min) lands inside
    // depends on geometry. Either 0 or 1 pixel, but no crash.
    ExPolygon tiny = triangle_mm(0, 0, 0.1, 0, 0, 0.1);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(tiny, 1.0, 400, 400, iw, ih, iwpr);

    if (!bm.empty()) {
        int count = popcnt_bm(bm);
        REQUIRE(count <= 1);
    }
    // If empty, that's also acceptable.
    SUCCEED();
}

TEST_CASE("scanline: thin spike with acute angle rasterizes without crash",
          "[BitmapGeometry][scanline]")
{
    // Very thin isoceles spike: base 0.5mm, height 20mm.
    // This challenges scanline intersection pairing because the two edges
    // are nearly parallel for most scanlines.
    ExPolygon spike(Points{
        Point(scaled<coord_t>(0.0),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(0.5),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(0.25),  scaled<coord_t>(20.0))
    });
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(spike, 0.5, 400, 400, iw, ih, iwpr);

    REQUIRE_FALSE(bm.empty());
    int count = popcnt_bm(bm);
    // Area = 0.5 * 0.5 * 20 = 5 mm^2 -> at 0.5mm/px -> 5/0.25 = 20 pixels ideally.
    // Scanline fill overestimates thin features. Accept 10 to 50.
    REQUIRE(count >= 10);
    REQUIRE(count <= 50);
}

// ===========================================================================
// 2. Hole subtraction
// ===========================================================================

TEST_CASE("hole: square with centered square hole",
          "[BitmapGeometry][hole]")
{
    // Outer: 10mm x 10mm -> at 0.5mm/px -> 20x20 = 400 pixels.
    // Hole: 4mm x 4mm centered at (3..7, 3..7) -> 8x8 = 64 pixels.
    // Net ~ 400 - 64 = 336, with boundary tolerance +-10.
    ExPolygon donut;
    donut.contour = Polygon(Points{
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(10.0))
    });
    // Hole: CW winding. Centered 3..7 on both axes.
    donut.holes.push_back(Polygon(Points{
        Point(scaled<coord_t>(3.0), scaled<coord_t>(3.0)),
        Point(scaled<coord_t>(3.0), scaled<coord_t>(7.0)),
        Point(scaled<coord_t>(7.0), scaled<coord_t>(7.0)),
        Point(scaled<coord_t>(7.0), scaled<coord_t>(3.0))
    }));

    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(donut, 0.5, 400, 400, iw, ih, iwpr);

    REQUIRE(iw == 20);
    REQUIRE(ih == 20);

    int total = popcnt_bm(bm);
    REQUIRE(total >= 326);
    REQUIRE(total <= 346);

    // Pixels deep inside the hole (pixel coords 7..12 on each axis, i.e. mm 3.5..6.0)
    // must be clear.
    for (int y = 7; y <= 12; ++y)
        for (int x = 7; x <= 12; ++x)
            REQUIRE_FALSE(px_set(bm, iwpr, x, y));

    // Corner pixels of outer must be set.
    REQUIRE(px_set(bm, iwpr, 0, 0));
    REQUIRE(px_set(bm, iwpr, 19, 19));
}

TEST_CASE("hole: hole touching outer boundary",
          "[BitmapGeometry][hole]")
{
    // Outer: 10mm x 10mm. Hole: 0..5, 0..10 (left half, touches boundary).
    // This leaves only the right half: ~200 pixels at 0.5mm/px.
    ExPolygon half;
    half.contour = Polygon(Points{
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(0.0),  scaled<coord_t>(10.0))
    });
    // Hole covers [0,5] x [0,10], CW winding.
    half.holes.push_back(Polygon(Points{
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(0.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(5.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(5.0), scaled<coord_t>(0.0))
    }));

    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(half, 0.5, 400, 400, iw, ih, iwpr);

    REQUIRE_FALSE(bm.empty());
    int total = popcnt_bm(bm);
    // Right half: 10x20 = 200 pixels. Boundary effects +-25.
    REQUIRE(total >= 175);
    REQUIRE(total <= 225);
    // No negative pixel counts, no crash. The left side should be mostly clear.
    for (int y = 1; y < ih - 1; ++y)
        for (int x = 0; x < 9; ++x) // well inside hole region (px 0..8 = mm 0..4.5)
            REQUIRE_FALSE(px_set(bm, iwpr, x, y));
}

// ===========================================================================
// 3. Collision detection with known geometry
// ===========================================================================

TEST_CASE("collision: two identical squares at same position collide",
          "[BitmapGeometry][collision]")
{
    // Rasterize a 5x5mm square at 1mm/px -> 5x5 bitmap.
    ExPolygon sq = rect_mm(0, 0, 5, 5);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 1.0, 200, 200, iw, ih, iwpr);
    REQUIRE(iw == 5);
    REQUIRE(ih == 5);

    // Plate: 200x200, stamp square at (0,0).
    int bw = 200, bh = 200, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0);

    // Same shape at same position: must collide.
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0));
}

TEST_CASE("collision: two 5x5 squares separated by 1 pixel do not collide",
          "[BitmapGeometry][collision]")
{
    ExPolygon sq = rect_mm(0, 0, 5, 5);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 1.0, 200, 200, iw, ih, iwpr);

    int bw = 200, bh = 200, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0);

    // Place second square at px=6 (1 pixel gap from px=5 which is first open column).
    // First square occupies cols 0..4. px=5 is clear. px=6 starts second square.
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 6, 0));
}

TEST_CASE("collision: two 5x5 squares overlapping by 1 pixel",
          "[BitmapGeometry][collision]")
{
    ExPolygon sq = rect_mm(0, 0, 5, 5);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 1.0, 200, 200, iw, ih, iwpr);

    int bw = 200, bh = 200, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0);

    // First square: cols 0..4. Place second at px=4 -> second occupies 4..8.
    // Column 4 overlaps.
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 4, 0));

    // Place at px=5 -> second occupies 5..9. No overlap.
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 5, 0));
}

TEST_CASE("collision: across 64-bit word boundary",
          "[BitmapGeometry][collision]")
{
    // Item is 5 pixels wide. Place at px=62 -> occupies 62,63,64,65,66.
    // This spans the word boundary at bit 63/64.
    ExPolygon sq = rect_mm(0, 0, 5, 2);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 1.0, 200, 200, iw, ih, iwpr);
    REQUIRE(iw == 5);
    REQUIRE(ih == 2);

    int bw = 200, bh = 200, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 62, 0);

    // Verify stamped bits span the word boundary.
    REQUIRE(px_set(plate, wpr, 62, 0));
    REQUIRE(px_set(plate, wpr, 63, 0));
    REQUIRE(px_set(plate, wpr, 64, 0));
    REQUIRE(px_set(plate, wpr, 65, 0));
    REQUIRE(px_set(plate, wpr, 66, 0));

    // Collision at the exact same position must detect it.
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 62, 0));

    // Item at px=64 overlaps at 64,65,66. Must collide.
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 64, 0));

    // Item at px=67 starts after 66. Must not collide.
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 67, 0));
}

TEST_CASE("collision: item at bed corners does not access out of bounds",
          "[BitmapGeometry][collision]")
{
    ExPolygon sq = rect_mm(0, 0, 3, 3);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 1.0, 100, 100, iw, ih, iwpr);
    REQUIRE(iw == 3);
    REQUIRE(ih == 3);

    int bw = 100, bh = 100, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    // Top-left corner: px=0, py=0.
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0));
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0);
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0));

    // Bottom-right corner: px = bw-iw, py = bh-ih.
    int corner_px = bw - iw;
    int corner_py = bh - ih;
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, corner_px, corner_py));
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, corner_px, corner_py);
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, corner_px, corner_py));
}

// ===========================================================================
// 4. Stamp + collide roundtrip
// ===========================================================================

TEST_CASE("roundtrip: stamp then collide same shape at same position",
          "[BitmapGeometry][roundtrip]")
{
    ExPolygon sq = rect_mm(0, 0, 8, 6);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 0.5, 400, 400, iw, ih, iwpr);

    int bw = 400, bh = 400, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 10, 10);
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 10, 10));
}

TEST_CASE("roundtrip: stamp shape, different shape at non-overlapping position does not collide",
          "[BitmapGeometry][roundtrip]")
{
    ExPolygon sq = rect_mm(0, 0, 5, 5);
    int iw1, ih1, iwpr1;
    auto bm1 = BitmapNester::rasterize(sq, 1.0, 200, 200, iw1, ih1, iwpr1);

    ExPolygon sq2 = rect_mm(0, 0, 3, 3);
    int iw2, ih2, iwpr2;
    auto bm2 = BitmapNester::rasterize(sq2, 1.0, 200, 200, iw2, ih2, iwpr2);

    int bw = 200, bh = 200, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    // Stamp first shape at origin.
    BitmapNester::stamp(plate, wpr, bw, bh, bm1, iwpr1, iw1, ih1, 0, 0);

    // Second shape at (50,50) — far away.
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm2, iwpr2, iw2, ih2, 50, 50));
}

TEST_CASE("roundtrip: stamp two shapes, third collides with either",
          "[BitmapGeometry][roundtrip]")
{
    ExPolygon sq = rect_mm(0, 0, 5, 5);
    int iw, ih, iwpr;
    auto bm = BitmapNester::rasterize(sq, 1.0, 200, 200, iw, ih, iwpr);

    int bw = 200, bh = 200, wpr = (bw + 63) / 64;
    std::vector<uint64_t> plate((size_t)wpr * bh, 0);

    // Stamp shape A at (0,0) and shape B at (20,20).
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 0, 0);
    BitmapNester::stamp(plate, wpr, bw, bh, bm, iwpr, iw, ih, 20, 20);

    // Third shape overlapping A at (3,3) -> cols 3..7, rows 3..7. Overlaps A at (3,3) and (4,4).
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 3, 3));

    // Third shape overlapping B at (22,22).
    REQUIRE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 22, 22));

    // Third shape at (10,10) — between A and B, no overlap.
    REQUIRE_FALSE(BitmapNester::collides(plate, wpr, bw, bh, bm, iwpr, iw, ih, 10, 10));
}

// ===========================================================================
// 5. Bitmap rotation consistency
// ===========================================================================

TEST_CASE("rotation: 10x5 rect at 0 deg vs 90 deg swaps width and height",
          "[BitmapGeometry][rotation]")
{
    // 10mm x 5mm at 1mm/px -> 10 x 5 at 0 degrees.
    ExPolygon rect_0 = rect_mm(0, 0, 10, 5);
    int iw0, ih0, iwpr0;
    auto bm0 = BitmapNester::rasterize(rect_0, 1.0, 200, 200, iw0, ih0, iwpr0);
    REQUIRE(iw0 == 10);
    REQUIRE(ih0 == 5);

    // Rotate 90 degrees. After rotation the bbox should be ~5 x 10.
    ExPolygon rect_90 = rect_mm(0, 0, 10, 5);
    rect_90.rotate(M_PI / 2.0);
    // Translate so bbox.min is at origin (rasterize uses bbox internally, this
    // just ensures the polygon isn't in negative coords which doesn't matter
    // for rasterize — it uses its own bbox. But let's be explicit.)
    int iw90, ih90, iwpr90;
    auto bm90 = BitmapNester::rasterize(rect_90, 1.0, 200, 200, iw90, ih90, iwpr90);

    REQUIRE(iw90 == 5);
    REQUIRE(ih90 == 10);

    // Pixel counts should be equal (same area, just rotated).
    REQUIRE(popcnt_bm(bm0) == popcnt_bm(bm90));
}

// ===========================================================================
// 6. Inflation correctness
// ===========================================================================

TEST_CASE("inflation: 10mm square inflated by 2mm yields ~14x14mm bitmap",
          "[BitmapGeometry][inflation]")
{
    // Start with 10mm x 10mm square.
    ExPolygon sq = rect_mm(0, 0, 10, 10);

    // Inflate by 2mm using offset_ex (the same path arrange() uses).
    ExPolygons inflated = offset_ex(sq, scaled<coord_t>(2.0));
    REQUIRE_FALSE(inflated.empty());
    ExPolygon isq = inflated.front();

    // Rasterize both at 0.5mm/px.
    int iw_orig, ih_orig, iwpr_orig;
    auto bm_orig = BitmapNester::rasterize(sq, 0.5, 400, 400, iw_orig, ih_orig, iwpr_orig);

    int iw_infl, ih_infl, iwpr_infl;
    auto bm_infl = BitmapNester::rasterize(isq, 0.5, 400, 400, iw_infl, ih_infl, iwpr_infl);

    // Original: 20x20 pixels.
    REQUIRE(iw_orig == 20);
    REQUIRE(ih_orig == 20);

    // Inflated: should be approximately 28x28 (14mm / 0.5mm/px).
    // Clipper rounds corners, so width/height might be 27-29.
    REQUIRE(iw_infl >= 27);
    REQUIRE(iw_infl <= 29);
    REQUIRE(ih_infl >= 27);
    REQUIRE(ih_infl <= 29);

    // Inflated pixel count must be strictly larger than original.
    int orig_count = popcnt_bm(bm_orig);
    int infl_count = popcnt_bm(bm_infl);
    REQUIRE(infl_count > orig_count);

    // Original is 400 pixels. Inflated square area ~196 mm^2 -> ~784 px.
    // Clipper rounds corners, so slightly less. Accept 700-800.
    REQUIRE(infl_count >= 700);
    REQUIRE(infl_count <= 800);
}
