// test_bitmap_overlap.cpp
// Overlap detection infrastructure and exhaustive overlap tests for BitmapNester.
//
// This file owns the overlap utility functions. Gandalf is writing the shared
// header for other cross-file utilities; these helpers stay local here per the
// task spec — do not factor them out until that header exists.

#define BITMAP_NESTER_TESTING

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"

#include <vector>
#include <cmath>
#include <string>
#include <sstream>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ============================================================================
// Overlap detection utilities
// ============================================================================

/// Apply translation + rotation stored in an ArrangePolygon and return the
/// resulting ExPolygon in absolute bed coordinates.
/// This mirrors ArrangePolygon::transformed_poly() but is spelled out here
/// so the logic is explicit and reviewable.
static ExPolygon apply_transform(const ArrangePolygon &ap)
{
    ExPolygon result = ap.poly;
    if (ap.rotation != 0.0)
        result.rotate(ap.rotation);
    result.translate(ap.translation.x(), ap.translation.y());
    return result;
}

/// Returns true when two placed ArrangePolygons overlap by more than
/// tolerance_mm² of area. Handles concave polygons by using Clipper's
/// intersection_ex — not just bounding-box overlap.
///
/// @param a            First item (must have valid translation/rotation set).
/// @param b            Second item (must have valid translation/rotation set).
/// @param tolerance_mm Minimum overlap area (mm²) to consider a collision.
///                     Default 0.01 mm² absorbs Clipper numerical slivers.
static bool items_overlap(const ArrangePolygon &a,
                          const ArrangePolygon &b,
                          double tolerance_mm = 0.01)
{
    ExPolygon pa = apply_transform(a);
    ExPolygon pb = apply_transform(b);

    ExPolygons inter = intersection_ex(ExPolygons{pa}, ExPolygons{pb});

    // Scale factor: 1 mm = scaled<double>(1.0) internal units.
    // Area in internal units is scaled² relative to mm².
    double tol_scaled_sq = scaled<double>(1.0) * scaled<double>(1.0) * tolerance_mm;

    for (const ExPolygon &seg : inter) {
        if (std::abs(seg.area()) > tol_scaled_sq)
            return true;
    }
    return false;
}

/// Check all N*(N-1)/2 pairs of arranged items on the same plate.
/// Skips UNARRANGED items and items on different plates.
/// Calls FAIL with a descriptive message on the first overlap found.
static void assert_no_overlaps(const ArrangePolygons &items,
                                double tolerance_mm = 0.01)
{
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;
        for (size_t j = i + 1; j < items.size(); ++j) {
            if (items[j].bed_idx != items[i].bed_idx) continue;
            if (items_overlap(items[i], items[j], tolerance_mm)) {
                std::ostringstream msg;
                msg << "Items " << i << " and " << j
                    << " overlap on plate " << items[i].bed_idx;
                FAIL(msg.str());
            }
        }
    }
}

/// Verify that every placed item's transformed bounding box is contained
/// within the bed bounding box, with up to tolerance_mm tolerance to absorb
/// raster quantization (the nester works at 0.5 mm/pixel).
static void assert_all_in_bounds(const ArrangePolygons &items,
                                  const BoundingBox &bed,
                                  double tolerance_mm = 1.0)
{
    coord_t tol = scaled<coord_t>(tolerance_mm);

    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;

        ExPolygon placed = apply_transform(items[i]);
        BoundingBox bb   = get_extents(placed);

        INFO("Item " << i << " bbox: ["
             << unscaled<double>(bb.min.x()) << ", "
             << unscaled<double>(bb.min.y()) << "] .. ["
             << unscaled<double>(bb.max.x()) << ", "
             << unscaled<double>(bb.max.y()) << "]");

        REQUIRE(bb.min.x() >= bed.min.x() - tol);
        REQUIRE(bb.min.y() >= bed.min.y() - tol);
        REQUIRE(bb.max.x() <= bed.max.x() + tol);
        REQUIRE(bb.max.y() <= bed.max.y() + tol);
    }
}

// ============================================================================
// Local shape builders
// ============================================================================

static ExPolygon ov_rect_mm(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x0), scaled<coord_t>(y1))
    });
}

/// Build an L-shaped ExPolygon (concave) in mm.
/// The L occupies a 20×20 bounding box with a 10×10 notch cut from the
/// top-right corner.
///
///   (0,20) +-------+ (10,20)
///          |       |
///   (0,10) +---+   + (10,10) — notch starts here
///              |   |
///              |   |
///   (10,0) +---+---+ (20,0)   <-- wait, let me keep it simpler
///
/// Actual shape: outer 20×20 minus top-right 10×10 notch.
///   Vertices (CCW): (0,0) (20,0) (20,10) (10,10) (10,20) (0,20)
static ExPolygon ov_L_shape_mm()
{
    return ExPolygon(Points{
        Point(scaled<coord_t>( 0.0), scaled<coord_t>( 0.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>( 0.0)),
        Point(scaled<coord_t>(20.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(10.0)),
        Point(scaled<coord_t>(10.0), scaled<coord_t>(20.0)),
        Point(scaled<coord_t>( 0.0), scaled<coord_t>(20.0))
    });
}

static BoundingBox ov_bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

static ArrangeParams ov_silent_params()
{
    ArrangeParams p;
    p.bed_shrink_x    = 0.0f;
    p.bed_shrink_y    = 0.0f;
    p.allow_rotations = false;
    p.progressind     = nullptr;
    return p;
}

static ArrangePolygon ov_ap(const ExPolygon &poly, int priority = 0)
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

// Helper: build an ArrangePolygon with explicit translation (mm) and no
// arrangement needed — used to set up detector self-tests directly.
static ArrangePolygon placed_ap(const ExPolygon &poly, double tx_mm, double ty_mm)
{
    ArrangePolygon ap = ov_ap(poly);
    ap.bed_idx    = 0;
    ap.translation = Vec2crd{scaled<coord_t>(tx_mm), scaled<coord_t>(ty_mm)};
    return ap;
}

// ============================================================================
// Section 1: Overlap detector self-tests
// ============================================================================

TEST_CASE("BitmapOverlap: detector finds overlap between two overlapping squares",
          "[BitmapOverlap]")
{
    // Two 10×10 squares, second placed at (5, 5) → 25 mm² overlap.
    ArrangePolygon a = placed_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0),  0.0, 0.0);
    ArrangePolygon b = placed_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0),  5.0, 5.0);

    REQUIRE(items_overlap(a, b));
}

TEST_CASE("BitmapOverlap: detector finds no overlap between separated squares",
          "[BitmapOverlap]")
{
    // Two 10×10 squares, 5 mm apart → no overlap.
    ArrangePolygon a = placed_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0),  0.0, 0.0);
    ArrangePolygon b = placed_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0), 15.0, 0.0);

    REQUIRE_FALSE(items_overlap(a, b));
}

TEST_CASE("BitmapOverlap: detector accepts edge-touching squares as non-overlapping",
          "[BitmapOverlap]")
{
    // Second square starts exactly where the first ends — shared edge, zero area.
    ArrangePolygon a = placed_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0),  0.0, 0.0);
    ArrangePolygon b = placed_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0), 10.0, 0.0);

    // Intersection area is zero (or a negligible sliver) — not an overlap.
    REQUIRE_FALSE(items_overlap(a, b));
}

TEST_CASE("BitmapOverlap: detector handles concave L-shapes that interlock without overlapping",
          "[BitmapOverlap]")
{
    // Two L-shapes that nest into each other's notch.
    //
    // L1 occupies: (0,0)-(20,10) bottom bar + (0,10)-(10,20) left bar.
    // L2 is the mirror: rotated 180° and translated so its notch aligns with L1's
    // corner piece. After rotation by π, the L sits differently so we place it
    // such that its "spike" fills L1's notch without overlapping.
    //
    // Concretely:
    //   L1 at (0,0): occupies cols 0-20 / rows 0-10, plus cols 0-10 / rows 10-20.
    //   L2 at (10,10): occupies cols 10-30 / rows 10-20, plus cols 20-30 / rows 0-10.
    //
    // The only shared region would be at x=10-20, y=10-20 — L1 has nothing there
    // (it's the notch) and L2 starts at (10,10) with its own notch in the
    // bottom-left. Let's verify they don't overlap.

    ExPolygon L = ov_L_shape_mm();

    // L1: untransformed at origin → occupies (0-20, 0-10) ∪ (0-10, 10-20).
    ArrangePolygon L1 = placed_ap(L, 0.0, 0.0);

    // L2: rotated 180°, then translated so it fits the complementary region.
    // After 180° rotation, the L's origin moves. The rotated shape spans
    // (-20, -20) .. (0, 0) in local coords, so translate by (30, 20) to land
    // it in (10, 0) .. (30, 20). The two L-shapes share no area.
    ExPolygon L_rot = L;
    L_rot.rotate(M_PI);  // rotate around (0,0)
    // Rotated bbox: (-20,-20)..(0,0). Translate to sit at x=10, y=0:
    //   translate by (10+20, 0+20) = (30, 20)
    ArrangePolygon L2;
    L2.poly              = L_rot;
    L2.bed_idx           = 0;
    L2.rotation          = 0.0; // rotation already baked into poly
    L2.translation       = Vec2crd{scaled<coord_t>(30.0), scaled<coord_t>(20.0)};
    L2.allowed_rotations = {0.0};

    REQUIRE_FALSE(items_overlap(L1, L2));
}

TEST_CASE("BitmapOverlap: detector catches overlapping concave L-shapes",
          "[BitmapOverlap]")
{
    // Two L-shapes at the same position — clearly overlapping.
    ExPolygon L = ov_L_shape_mm();

    ArrangePolygon L1 = placed_ap(L,  0.0, 0.0);
    ArrangePolygon L2 = placed_ap(L,  2.0, 2.0); // shifted only 2 mm — substantial overlap

    REQUIRE(items_overlap(L1, L2));
}

// ============================================================================
// Section 2: Multi-bed overflow stress test
// ============================================================================

TEST_CASE("BitmapOverlap: 30 mixed-size items — no intra-plate overlaps",
          "[BitmapOverlap]")
{
    // Items of three different sizes to exercise plate overflow.
    BoundingBox bed = ov_bed_mm(256.0, 210.0);
    ArrangeParams p = ov_silent_params();

    ArrangePolygons items;
    items.reserve(30);
    for (int i = 0; i < 10; ++i)
        items.push_back(ov_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0)));
    for (int i = 0; i < 10; ++i)
        items.push_back(ov_ap(ov_rect_mm(0.0, 0.0, 30.0, 20.0)));
    for (int i = 0; i < 10; ++i)
        items.push_back(ov_ap(ov_rect_mm(0.0, 0.0, 15.0, 15.0)));

    ArrangePolygons excludes;
    BitmapNester::arrange(items, excludes, bed, p);

    // Every item must be either arranged (bed_idx >= 0) or UNARRANGED — never
    // some garbage value.
    for (size_t i = 0; i < items.size(); ++i) {
        INFO("Item " << i << " bed_idx = " << items[i].bed_idx);
        REQUIRE((items[i].bed_idx >= 0 || items[i].bed_idx == UNARRANGED));
    }

    // No overlap within any single plate.
    assert_no_overlaps(items);

    // All items placed (bed + overflow plates are large enough).
    for (size_t i = 0; i < items.size(); ++i) {
        INFO("Item " << i << " unarranged");
        REQUIRE(items[i].bed_idx != UNARRANGED);
    }
}

// ============================================================================
// Section 3: Exclude zone overlap
// ============================================================================

TEST_CASE("BitmapOverlap: items do not overlap exclude zones",
          "[BitmapOverlap]")
{
    // Five exclude zones in the center of a 256×210 mm bed.
    // 20 items arranged around them must not land on any exclude.

    BoundingBox bed = ov_bed_mm(256.0, 210.0);
    ArrangeParams p = ov_silent_params();

    ArrangePolygons excludes;
    excludes.reserve(5);
    // Excludes: five 20×20 zones in a row across the middle.
    for (int i = 0; i < 5; ++i) {
        ArrangePolygon ex = ov_ap(ov_rect_mm(0.0, 0.0, 20.0, 20.0));
        ex.bed_idx    = 0;
        double ex_x = 30.0 + i * 40.0;   // at 30, 70, 110, 150, 190
        double ex_y = 90.0;               // vertical center of bed
        ex.translation = Vec2crd{scaled<coord_t>(ex_x), scaled<coord_t>(ex_y)};
        excludes.push_back(ex);
    }

    ArrangePolygons items;
    items.reserve(20);
    for (int i = 0; i < 20; ++i)
        items.push_back(ov_ap(ov_rect_mm(0.0, 0.0, 10.0, 10.0)));

    BitmapNester::arrange(items, excludes, bed, p);

    // Each placed item must not overlap any exclude zone.
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].bed_idx != 0) continue; // excludes only on plate 0

        for (size_t e = 0; e < excludes.size(); ++e) {
            // excludes[e] is already placed (has its translation set).
            bool overlap = items_overlap(items[i], excludes[e]);
            if (overlap) {
                std::ostringstream msg;
                msg << "Item " << i << " overlaps exclude zone " << e;
                FAIL(msg.str());
            }
        }
    }

    // Also verify items don't overlap each other.
    assert_no_overlaps(items);

    // Count placed items so the test registers an assertion.
    int placed = 0;
    for (auto &it : items)
        if (it.bed_idx >= 0) ++placed;
    REQUIRE(placed > 0);
}

// ============================================================================
// Section 4: assert_all_in_bounds spot-check
// ============================================================================

TEST_CASE("BitmapOverlap: arranged items stay within bed bounds",
          "[BitmapOverlap]")
{
    BoundingBox bed = ov_bed_mm(256.0, 210.0);
    ArrangeParams p = ov_silent_params();

    ArrangePolygons items;
    for (int i = 0; i < 15; ++i)
        items.push_back(ov_ap(ov_rect_mm(0.0, 0.0, 20.0, 20.0)));

    ArrangePolygons excludes;
    BitmapNester::arrange(items, excludes, bed, p);

    assert_all_in_bounds(items, bed, /*tolerance_mm=*/1.0);
}
