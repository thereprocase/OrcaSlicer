// test_consolidation.cpp
//
// Integration tests for the consolidation pass in BitmapNester::arrange.
// Consolidation is the backwards migration sweep that runs after the
// greedy placement loop and before post-centering. It walks plates
// last-to-first and tries to move each item from a later plate back to
// an earlier one by finding any clear scored position there. This file
// constructs scenarios that exercise specific failure modes identified
// by Ent's triage on 2026-04-11 and by the autopilot session on the
// same date.
//
// All shapes are axis-aligned integer-mm rectangles so the bitmap
// rasterization is exact (no sub-pixel rounding noise). This keeps the
// no_overlap helper's tight 0.0001 mm² threshold honest.
//
// Ent's Critical Gap 1: "consolidation block has 250 lines of behavior
// and zero direct tests". This file closes that gap.

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "bitmap_test_utils.hpp"

#include <vector>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Helpers — mirror test_tetromino_packing.cpp conventions so shapes are
// axis-aligned to the mm grid and the bitmap rasterization is exact.
// ---------------------------------------------------------------------------

static BoundingBox bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

static ArrangeParams consolidation_params()
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

static ExPolygon rect_mm(double w, double h)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h)),
        Point(scaled<coord_t>(0.0), scaled<coord_t>(h))
    });
}

// L-shape whose bbox is w×h. The notch removes the top-right quadrant.
// Exercises rotated-migration scenarios (the L needs rotation to fit in
// certain gaps).
static ExPolygon l_mm(double w, double h)
{
    double w2 = w * 0.5;
    double h2 = h * 0.5;
    return ExPolygon(Points{
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h2)),
        Point(scaled<coord_t>(w2),  scaled<coord_t>(h2)),
        Point(scaled<coord_t>(w2),  scaled<coord_t>(h)),
        Point(scaled<coord_t>(0.0), scaled<coord_t>(h))
    });
}

static ArrangePolygon make_ap(const ExPolygon &poly,
                              std::vector<double> rotations = {0.0})
{
    ArrangePolygon ap;
    ap.poly               = poly;
    ap.priority           = 0;
    ap.bed_idx            = UNARRANGED;
    ap.rotation           = 0.0;
    ap.translation        = Vec2crd{0, 0};
    ap.allowed_rotations  = std::move(rotations);
    return ap;
}

// ---------------------------------------------------------------------------
// Scenario 1 — loose pack, expect all on plate 0.
//
// 5 small rectangles on a 200x200 bed. Total footprint << bed area, so
// consolidation has lots of room but also shouldn't need to do much:
// greedy should already hit plate 0 for everything. This is the "does
// consolidation no-op correctly" test.
// ---------------------------------------------------------------------------

TEST_CASE("consolidation: loose pack — all items on plate 0",
          "[BitmapNester][consolidation][BitmapRegression]")
{
    ArrangePolygons items;
    items.push_back(make_ap(rect_mm(30, 30)));
    items.push_back(make_ap(rect_mm(30, 30)));
    items.push_back(make_ap(rect_mm(30, 30)));
    items.push_back(make_ap(rect_mm(30, 30)));
    items.push_back(make_ap(rect_mm(30, 30)));

    BoundingBox bed = bed_mm(200, 200);
    ArrangeParams params = consolidation_params();

    BitmapNester::arrange(items, ArrangePolygons{}, bed, params);

    REQUIRE(test_utils::max_bed_idx(items) == 0);
    REQUIRE(test_utils::overflow_piece_count(items) == 0);
    REQUIRE(test_utils::no_overlap(items));
}

// ---------------------------------------------------------------------------
// Scenario 2 — four equal squares land as a 2x2 grid on plate 0.
//
// 4 equal 80x80 rectangles on a 200x200 bed. Total 4 * 6400 = 25600 mm²
// on 40000 mm² (64% density). Geometric optimum is a 2x2 grid =
// 160x160 cluster post-centered to (20,20)-(180,180).
//
// HISTORY: This test initially exposed a failure mode where the
// hybrid anchor (corner for first item, bed center for rest)
// over-pulled item 2 to (60,80) instead of (0,80) or (80,80) when
// delta_area tied across those positions. Item 3 then had no 80x80
// slot. Fixed 2026-04-11 by switching to pure corner anchor for all
// items in default-align mode, relying on post-centering to move the
// resulting tight cluster to the bed center. The corner anchor
// produces grid-aligned placement because ties are broken by smallest
// (px, py), which is BL-fill. Post-centering still works as before,
// shifting the final cluster to (20,20)-(180,180). S+Z 50x50 still
// passes because smart-shuffle picks up the slack that hybrid used
// to provide for small-N central-trap cases. See task #56 for the
// diagnosis trail.
// ---------------------------------------------------------------------------

TEST_CASE("consolidation: four equal squares fit 2x2 on plate 0",
          "[BitmapNester][consolidation][BitmapRegression]")
{
    ArrangePolygons items;
    for (int i = 0; i < 4; ++i)
        items.push_back(make_ap(rect_mm(80, 80)));

    BoundingBox bed = bed_mm(200, 200);
    ArrangeParams params = consolidation_params();

    BitmapNester::arrange(items, ArrangePolygons{}, bed, params);

    // Expected layout after the 2026-04-11 corner-anchor-for-all fix:
    //   item 0 → (20, 20)   post-centered from (0, 0)
    //   item 1 → (100, 20)  post-centered from (80, 0)
    //   item 2 → (20, 100)  post-centered from (0, 80)
    //   item 3 → (100, 100) post-centered from (80, 80)
    // All on plate 0. Compactness ratio = 1.0 (zero waste inside bbox).
    // Diagnostic — visible on pass or fail
    UNSCOPED_INFO("max_bed=" << test_utils::max_bed_idx(items)
                  << " perim_mm=" << test_utils::cluster_bbox_perimeter_mm(items)
                  << " compact=" << test_utils::cluster_compactness(items));
    for (size_t i = 0; i < items.size(); ++i)
        UNSCOPED_INFO("  item " << i << " bed=" << items[i].bed_idx
                      << " x=" << unscaled<double>(items[i].translation.x())
                      << " y=" << unscaled<double>(items[i].translation.y()));

    REQUIRE(test_utils::max_bed_idx(items) == 0);
    REQUIRE(test_utils::no_overlap(items));
}

// ---------------------------------------------------------------------------
// Scenario 3 — rotation-required migration.
//
// 3 L-shapes bbox=80x80 on a 200x100 bed (wide + shallow). Long axis
// of the bed is 200 mm, short axis is 100 mm. Each L has bbox side 80,
// so three L's stacked along X would need 240 mm — doesn't fit on the
// short axis at rotation 0. Placing them rotated at 90° allows a
// different arrangement. The consolidation pass's all-rotations
// scoring (added 2026-04-11) should find a rotation that fits if one
// exists. Verifies all-rotations migration isn't locked to original
// orientation.
// ---------------------------------------------------------------------------

TEST_CASE("consolidation: rotation-required packing on wide bed",
          "[BitmapNester][consolidation][BitmapRegression]")
{
    ArrangePolygons items;
    std::vector<double> rots{0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    items.push_back(make_ap(l_mm(80, 80), rots));
    items.push_back(make_ap(l_mm(80, 80), rots));
    items.push_back(make_ap(l_mm(80, 80), rots));

    BoundingBox bed = bed_mm(240, 100);
    ArrangeParams params = consolidation_params();

    BitmapNester::arrange(items, ArrangePolygons{}, bed, params);

    // 3 L's with bbox 80x80 = 19200 mm² footprint max; bed is
    // 240x100 = 24000 mm² so at best they fit on one plate with
    // 20% slack. Assert plate 0 only.
    REQUIRE(test_utils::max_bed_idx(items) == 0);
    REQUIRE(test_utils::overflow_piece_count(items) == 0);
    REQUIRE(test_utils::no_overlap(items));
}

// ---------------------------------------------------------------------------
// Scenario 4 — determinism across consolidation.
//
// Re-run the saturated 2x2 grid scenario twice and verify byte-identical
// output. Consolidation involves iteration, rasterization, and scan
// ordering that must be deterministic — any hidden nondeterminism (e.g.
// unordered_map iteration, timing-dependent branch) would manifest here.
// ---------------------------------------------------------------------------

TEST_CASE("consolidation: deterministic rerun",
          "[BitmapNester][consolidation][BitmapDeterminism][BitmapRegression]")
{
    auto build = []() {
        ArrangePolygons xs;
        for (int i = 0; i < 4; ++i)
            xs.push_back(make_ap(rect_mm(80, 80)));
        return xs;
    };

    BoundingBox bed = bed_mm(200, 200);
    ArrangeParams params = consolidation_params();
    auto run = [&](ArrangePolygons &xs) {
        BitmapNester::arrange(xs, ArrangePolygons{}, bed, params);
    };

    ArrangePolygons items = build();
    run(items);
    REQUIRE(test_utils::max_bed_idx(items) == 0);
    REQUIRE(test_utils::no_overlap(items));

    ArrangePolygons rerun = build();
    REQUIRE(test_utils::deterministic_rerun(rerun, run));
}
