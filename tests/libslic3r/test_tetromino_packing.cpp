// test_tetromino_packing.cpp
//
// Packing regression tests for the concave BitmapNester using the seven
// tetrominoes as probe shapes. See docs/PACKING_METRICS.md for the
// metric philosophy and the Three Seers' vote that adopted plate-count
// regression + determinism + overflow piece count.
//
// Each scenario:
//   1. Builds the piece multiset.
//   2. Scatters starting translations with a fixed-seed std::mt19937.
//   3. Sizes the bed loosely enough that the parts MUST fit on plate 0
//      under any reasonable nester. The bed size is the implicit
//      threshold — if the algorithm regresses, it overflows.
//   4. Asserts:
//      - max_bed_idx == 0           (#12, plate-count regression)
//      - overflow_piece_count == 0  (#4, no spillover)
//      - no_overlap                 (sanity)
//      - deterministic_rerun        (#13, byte-identical second run)
//
// We deliberately do NOT measure cluster bbox density. That metric was
// rejected by all three seers because it penalizes round clusters even
// when they fit identically — see docs/PACKING_METRICS.md for the full
// deliberation.

#include <catch2/catch_all.hpp>

// ClipperUtils.hpp must precede BitmapNester.hpp — the nester uses
// offset_ex inline in its placement loop. Match the include order in
// test_bitmap_nester.cpp.
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "bitmap_test_utils.hpp"

#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Local helpers (kept local rather than shared with test_bitmap_nester.cpp
// to avoid link-order issues on `static` helpers across translation units).
// ---------------------------------------------------------------------------

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
    p.bed_shrink_x = 0.0f;
    p.bed_shrink_y = 0.0f;
    p.allow_rotations = false;
    p.progressind = nullptr;
    return p;
}

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

// ---------------------------------------------------------------------------
// Tetromino factories. Unit = 10 mm. Canonical orientations.
// Each piece is exactly 4 unit squares = 400 mm^2.
// ---------------------------------------------------------------------------

static constexpr double U = 10.0; // mm per unit square
static constexpr double TETRO_AREA_MM2 = 4.0 * U * U;

static ExPolygon poly_from_mm(std::initializer_list<std::pair<double, double>> pts)
{
    Points ps;
    ps.reserve(pts.size());
    for (auto &p : pts)
        ps.emplace_back(scaled<coord_t>(p.first), scaled<coord_t>(p.second));
    return ExPolygon(ps);
}

static ExPolygon tetro_I() {
    return poly_from_mm({{0,0}, {4*U,0}, {4*U,U}, {0,U}});
}
static ExPolygon tetro_O() {
    return poly_from_mm({{0,0}, {2*U,0}, {2*U,2*U}, {0,2*U}});
}
static ExPolygon tetro_T() {
    return poly_from_mm({
        {0,U}, {U,U}, {U,0}, {2*U,0}, {2*U,U}, {3*U,U}, {3*U,2*U}, {0,2*U}
    });
}
static ExPolygon tetro_S() {
    return poly_from_mm({
        {0,0}, {2*U,0}, {2*U,U}, {3*U,U}, {3*U,2*U}, {U,2*U}, {U,U}, {0,U}
    });
}
static ExPolygon tetro_Z() {
    return poly_from_mm({
        {0,U}, {U,U}, {U,0}, {3*U,0}, {3*U,U}, {2*U,U}, {2*U,2*U}, {0,2*U}
    });
}
static ExPolygon tetro_J() {
    return poly_from_mm({
        {0,0}, {U,0}, {U,U}, {3*U,U}, {3*U,2*U}, {0,2*U}
    });
}
static ExPolygon tetro_L() {
    return poly_from_mm({
        {0,U}, {2*U,U}, {2*U,0}, {3*U,0}, {3*U,2*U}, {0,2*U}
    });
}

// ---------------------------------------------------------------------------
// Scatter + driver.
// ---------------------------------------------------------------------------

// Randomize starting translations with a fixed seed so results are
// reproducible across runs.
static void scatter(ArrangePolygons &items, uint32_t seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, 300.0);
    for (auto &it : items) {
        double tx = dist(rng);
        double ty = dist(rng);
        it.translation = Vec2crd{scaled<coord_t>(tx), scaled<coord_t>(ty)};
    }
}

static void enable_quarter_rotations(ArrangePolygons &items)
{
    std::vector<double> rots{0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    for (auto &it : items)
        it.allowed_rotations = rots;
}

static ArrangeParams benchmark_params()
{
    ArrangeParams p = no_shrink_params();
    p.min_obj_distance  = 0;     // measure pure algorithmic packing, no brim
    p.allow_rotations   = true;
    p.use_concave_shapes = true;
    // Leave do_final_align at its default (true) so determinism covers
    // the post-centering pass too.
    return p;
}

// Common driver: build the piece set, scatter, arrange on the given bed,
// then assert the Three Seers' metrics in order.
template <class Factory>
static void verify_packs_on_one_plate(const std::string &name,
                                      Factory build_pieces,
                                      double bed_mm_w,
                                      double bed_mm_h,
                                      uint32_t seed)
{
    auto build_input = [&]() {
        ArrangePolygons xs;
        for (auto &poly : build_pieces())
            xs.push_back(make_ap(poly));
        enable_quarter_rotations(xs);
        scatter(xs, seed);
        return xs;
    };

    BoundingBox bed = make_bed_mm(bed_mm_w, bed_mm_h);
    ArrangeParams params = benchmark_params();

    auto run = [&](ArrangePolygons &xs) {
        BitmapNester::arrange(xs, ArrangePolygons{}, bed, params);
    };

    // Primary run.
    ArrangePolygons items = build_input();
    run(items);

    UNSCOPED_INFO("Scenario:    " << name);
    UNSCOPED_INFO("  bed       = " << bed_mm_w << " x " << bed_mm_h << " mm");
    UNSCOPED_INFO("  pieces    = " << items.size());
    UNSCOPED_INFO("  area mm^2 = " << items.size() * TETRO_AREA_MM2);
    UNSCOPED_INFO("  max bed   = " << test_utils::max_bed_idx(items));
    UNSCOPED_INFO("  overflow  = " << test_utils::overflow_piece_count(items));

    // Metric #12: every piece on plate 0.
    REQUIRE(test_utils::max_bed_idx(items) == 0);
    // Metric #4: explicit overflow tripwire.
    REQUIRE(test_utils::overflow_piece_count(items) == 0);
    // Sanity: no two items occupy the same pixel.
    REQUIRE(test_utils::no_overlap(items));

    // Metric #13: determinism. Run a second time on a fresh copy of the
    // same input and assert byte-identical (bed_idx, rotation, translation,
    // itemid) tuples.
    ArrangePolygons rerun_input = build_input();
    REQUIRE(test_utils::deterministic_rerun(rerun_input, run));
}

// ---------------------------------------------------------------------------
// Scenario 1 — "Mix 14": two of each tetromino on a 140x140 bed.
// 14 * 400 = 5600 mm^2 of material on 19600 mm^2 of bed = 28% density.
// Loose enough that center-greedy + consolidation should fit everything
// regardless of starting scatter.
// ---------------------------------------------------------------------------

TEST_CASE("tetromino: 14 mixed pieces fit on plate 0",
          "[BitmapNester][tetromino][packing]")
{
    auto build = []() {
        std::vector<ExPolygon> v;
        for (int i = 0; i < 2; ++i) {
            v.push_back(tetro_I());
            v.push_back(tetro_O());
            v.push_back(tetro_T());
            v.push_back(tetro_S());
            v.push_back(tetro_Z());
            v.push_back(tetro_J());
            v.push_back(tetro_L());
        }
        return v;
    };
    verify_packs_on_one_plate("Mix 14", build, 140.0, 140.0, 0xC0FFEEu);
}

// ---------------------------------------------------------------------------
// Scenario 2 — "Four I-pieces": four 1x4 bars on a 50x50 bed.
// 1600 / 2500 = 64% density. The classical 4xI tile of a 4x4 square.
// ---------------------------------------------------------------------------

TEST_CASE("tetromino: four I-pieces fit on plate 0",
          "[BitmapNester][tetromino][packing]")
{
    auto build = []() {
        return std::vector<ExPolygon>{tetro_I(), tetro_I(), tetro_I(), tetro_I()};
    };
    verify_packs_on_one_plate("Four I", build, 50.0, 50.0, 0xBADA55u);
}

// ---------------------------------------------------------------------------
// Scenario 3 — "Two L-pieces": same chirality, can't tile a rectangle.
// Best achievable cluster bbox is 60x20 = 1200 mm^2 for 800 mm^2 of
// material. We just need them both on plate 0; the bed is 80x40 (3200
// mm^2 = 25% density) so any arrangement fits.
// ---------------------------------------------------------------------------

TEST_CASE("tetromino: two L-pieces fit on plate 0",
          "[BitmapNester][tetromino][packing]")
{
    auto build = []() {
        return std::vector<ExPolygon>{tetro_L(), tetro_L()};
    };
    verify_packs_on_one_plate("Two L", build, 80.0, 40.0, 0xDEADBEEFu);
}

// ---------------------------------------------------------------------------
// Scenario 4 — "S and Z": two S + two Z on a 50x50 bed (64% density).
// This is the smart-shuffle stress test. Center-greedy placement fragments
// the surrounding space when the first piece lands dead-center, and a
// single-pass greedy can't recover. The smart-shuffle restart driver
// re-runs placement with shuffled within-priority orderings until one
// fits — N=4 means 8 attempts max, all in milliseconds.
//
// If this test fails, smart-shuffle is broken or insufficient for tight
// 4-piece interlocks.
// ---------------------------------------------------------------------------

TEST_CASE("tetromino: two S + two Z fit on plate 0",
          "[BitmapNester][tetromino][packing]")
{
    auto build = []() {
        return std::vector<ExPolygon>{tetro_S(), tetro_S(), tetro_Z(), tetro_Z()};
    };
    verify_packs_on_one_plate("S and Z", build, 50.0, 50.0, 0xF00DCAFEu);
}

// ---------------------------------------------------------------------------
// Scenario 5 — "Mix 40": 40 mixed tetrominoes on a 220x220 bed.
// 16000 / 48400 = 33% density. The stress test.
// ---------------------------------------------------------------------------

TEST_CASE("tetromino: 40 mixed pieces fit on plate 0",
          "[BitmapNester][tetromino][packing]")
{
    auto build = []() {
        std::vector<ExPolygon> v;
        for (int i = 0; i < 6; ++i) v.push_back(tetro_I());
        for (int i = 0; i < 6; ++i) v.push_back(tetro_O());
        for (int i = 0; i < 6; ++i) v.push_back(tetro_T());
        for (int i = 0; i < 5; ++i) v.push_back(tetro_S());
        for (int i = 0; i < 5; ++i) v.push_back(tetro_Z());
        for (int i = 0; i < 6; ++i) v.push_back(tetro_J());
        for (int i = 0; i < 6; ++i) v.push_back(tetro_L());
        return v;
    };
    verify_packs_on_one_plate("Mix 40", build, 220.0, 220.0, 0x5EED1234u);
}

// ---------------------------------------------------------------------------
// Scenario 6 — "Tetris 125 balanced": 125 tetrominoes on a 240x240 bed.
// Mirrors tests/data/tetris_plate_240x240_balanced.stl (manual visual
// fixture). 125 * 400 = 50,000 mm^2 of material on 57,600 mm^2 of bed =
// 86.8% density. K_min = 1. This is THE regression gate for the
// pack-as-best work — the 2026-04-11 investigation identified heavy
// spillover (~25-30 pieces on plate 02) on the current binary and
// hypothesized four root causes (smart-shuffle disabled at N > 16,
// consolidation cstride too coarse, no refine pass, missing boundary
// positions). This test exists to prove those fixes land.
//
// Distribution: 18,18,18,18,18,18,17 across I,O,T,S,Z,J,L = 125 total.
// ---------------------------------------------------------------------------

// Floor assertion for the 125-piece regression gate. True goal is
// `max_bed_idx == 0 && overflow == 0` (K_min = 1 at 86.8% density) but
// the current algorithm tops out near 80.6% density (~9 piece overflow).
// The floor is set just above the current measured value so any
// regression trips it, and it tightens as we improve. See the Sprint-2
// log: each improvement should ratchet these thresholds downward.
TEST_CASE("tetromino: 125 balanced pieces fit on 240x240 plate",
          "[BitmapNester][tetromino][packing][BitmapRegression]")
{
    static constexpr int TETRIS125_MAX_OVERFLOW = 10;
    static constexpr int TETRIS125_MAX_BED_IDX = 1;

    auto build_pieces = []() {
        std::vector<ExPolygon> v;
        v.reserve(125);
        for (int i = 0; i < 18; ++i) v.push_back(tetro_I());
        for (int i = 0; i < 18; ++i) v.push_back(tetro_O());
        for (int i = 0; i < 18; ++i) v.push_back(tetro_T());
        for (int i = 0; i < 18; ++i) v.push_back(tetro_S());
        for (int i = 0; i < 18; ++i) v.push_back(tetro_Z());
        for (int i = 0; i < 18; ++i) v.push_back(tetro_J());
        for (int i = 0; i < 17; ++i) v.push_back(tetro_L());
        return v;
    };

    auto build_input = [&]() {
        ArrangePolygons xs;
        for (auto &poly : build_pieces())
            xs.push_back(make_ap(poly));
        enable_quarter_rotations(xs);
        scatter(xs, 0x7E771500u);
        return xs;
    };

    BoundingBox bed = make_bed_mm(240.0, 240.0);
    ArrangeParams params = benchmark_params();
    auto run = [&](ArrangePolygons &xs) {
        BitmapNester::arrange(xs, ArrangePolygons{}, bed, params);
    };

    ArrangePolygons items = build_input();
    run(items);

    int max_bed = test_utils::max_bed_idx(items);
    int overflow = test_utils::overflow_piece_count(items);

    UNSCOPED_INFO("Scenario:    Tetris 125");
    UNSCOPED_INFO("  bed       = 240 x 240 mm");
    UNSCOPED_INFO("  pieces    = 125");
    UNSCOPED_INFO("  area mm^2 = 50000");
    UNSCOPED_INFO("  max bed   = " << max_bed);
    UNSCOPED_INFO("  overflow  = " << overflow);
    UNSCOPED_INFO("  floor     = max_bed<=" << TETRIS125_MAX_BED_IDX
                  << " overflow<=" << TETRIS125_MAX_OVERFLOW);

    // Regression floor. Tighten these constants as the algorithm improves.
    REQUIRE(max_bed <= TETRIS125_MAX_BED_IDX);
    REQUIRE(overflow <= TETRIS125_MAX_OVERFLOW);
    REQUIRE(test_utils::no_overlap(items));

    // Determinism always required.
    ArrangePolygons rerun_input = build_input();
    REQUIRE(test_utils::deterministic_rerun(rerun_input, run));
}

// ---------------------------------------------------------------------------
// Sanity: every tetromino factory produces a 4-square-area polygon.
// ---------------------------------------------------------------------------

TEST_CASE("tetromino factories have exact area = 4 unit squares",
          "[BitmapNester][tetromino][packing]")
{
    struct Entry { const char *name; ExPolygon (*fn)(); };
    Entry entries[] = {
        {"I", tetro_I}, {"O", tetro_O}, {"T", tetro_T},
        {"S", tetro_S}, {"Z", tetro_Z}, {"J", tetro_J}, {"L", tetro_L},
    };
    for (auto &e : entries) {
        DYNAMIC_SECTION("piece " << e.name) {
            ExPolygon p = e.fn();
            double area_mm2 = unscaled<double>(unscaled<double>(p.area()));
            REQUIRE_THAT(area_mm2,
                         Catch::Matchers::WithinAbs(TETRO_AREA_MM2, 1e-6));
        }
    }
}
