// test_real_mix_packing.cpp
//
// Regression gate using real-ish concave OBJ fixtures from tests/data/.
// The tetris125 fixture is a synthetic extremum (equal area, clumped
// input, interlock-best-when-mixed) that punishes FFD-style heuristics
// that work well on realistic mixed workloads. This test exists so
// algorithmic changes can be validated against a mixed-size mixed-
// concavity distribution instead of only tetris125 — see
// memory/feedback_tetrominoes_quant_not_direction.md.
//
// For each fixture we load the OBJ, project its vertices to XY, take
// the 2D convex hull (cheap and deterministic), and feed that as the
// ArrangePolygon::poly. Not identical to the production path (which
// populates concave_regions from per-layer slicing), but it exercises
// the same BitmapNester code path and gives us a handful of honest
// mixed-size shapes to pack.

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Geometry/ConvexHull.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "../test_utils.hpp"
#include "bitmap_test_utils.hpp"

#include <random>
#include <string>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static BoundingBox make_bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

static ArrangeParams benchmark_params()
{
    ArrangeParams p;
    p.bed_shrink_x = 0.0f;
    p.bed_shrink_y = 0.0f;
    p.allow_rotations = true;
    // Small brim so the nester pads items by ~0.5 mm. Real concave
    // polygon edges fall off the 0.5 mm pixel grid, introducing
    // sub-pixel rasterization rounding. A zero-brim arrangement can
    // produce items that are non-overlapping in the bitmap (which is
    // what the nester checks) but overlap each other in the exact
    // polygon representation by a few mm^2 — above the no_overlap
    // helper's 0.0001 mm^2 threshold. One pixel of padding gives the
    // nester guaranteed bitmap separation that dominates the noise.
    p.min_obj_distance = scaled<coord_t>(1.0);
    p.use_concave_shapes = true;
    p.progressind = nullptr;
    return p;
}

// Load an OBJ fixture and reduce its vertices to a 2D ExPolygon (convex
// hull of the XY projection). Returns an empty ExPolygon if the load
// fails. The test SKIPs in that case — we don't want CI flakes just
// because a fixture path changed, but any successful load should
// produce a non-trivial polygon.
static ExPolygon load_fixture_footprint(const std::string &obj_name)
{
    TriangleMesh mesh = load_model(obj_name);
    if (mesh.its.vertices.empty()) return ExPolygon{};

    Points pts;
    pts.reserve(mesh.its.vertices.size());
    for (const auto &v : mesh.its.vertices) {
        pts.emplace_back(scaled<coord_t>((double)v.x()),
                         scaled<coord_t>((double)v.y()));
    }

    Polygon hull = Geometry::convex_hull(pts);
    if (hull.points.size() < 3) return ExPolygon{};

    // Shift so min corner is at (0, 0), matching the tetromino test's
    // convention. Makes the arrange input positionally neutral.
    BoundingBox bb = hull.bounding_box();
    hull.translate(-bb.min.x(), -bb.min.y());

    ExPolygon ep;
    ep.contour = hull;
    return ep;
}

static ArrangePolygon make_ap(const ExPolygon &poly)
{
    ArrangePolygon ap;
    ap.poly       = poly;
    ap.priority   = 0;
    ap.bed_idx    = UNARRANGED;
    ap.rotation   = 0.0;
    ap.translation = Vec2crd{0, 0};
    ap.allowed_rotations = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    return ap;
}

static void scatter(ArrangePolygons &items, uint32_t seed, double spread_mm)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, spread_mm);
    for (auto &it : items) {
        it.translation = Vec2crd{scaled<coord_t>(dist(rng)),
                                 scaled<coord_t>(dist(rng))};
    }
}

// ---------------------------------------------------------------------------
// Scenario: mixed concave fixtures on a 256 x 256 bed.
//
// Counts are tuned so total footprint area is well under 50% bed area —
// the test isn't measuring maximum density, it's measuring that a mixed
// real-world distribution lands on the primary plate without
// catastrophic spillover. An algorithmic change that pushes these parts
// to plate 1 needs a very good reason.
// ---------------------------------------------------------------------------

TEST_CASE("real mix: concave OBJ fixtures fit on 256x256 plate",
          "[BitmapNester][real-mix][packing][BitmapRegression]")
{
    // Fixture list. Each entry: (filename, copies). Files live in
    // tests/data/ and load_model() resolves TEST_DATA_DIR for us.
    struct Entry { const char *name; int copies; };
    const Entry entries[] = {
        {"ipadstand.obj",            4},
        {"extruder_idler.obj",       4},
        {"frog_legs.obj",            3},
        {"cube_with_concave_hole.obj", 6},
        {"small_dorito.obj",         6},
    };

    // Load footprints first so we can skip the test cleanly if any
    // fixture is missing (CI portability).
    std::vector<ExPolygon> shapes;
    std::vector<int>       copies;
    for (const auto &e : entries) {
        ExPolygon fp = load_fixture_footprint(e.name);
        if (fp.contour.points.empty()) {
            WARN("real mix: fixture not loadable: " << e.name);
            return; // skip — no regression catch, but no false failure
        }
        shapes.push_back(std::move(fp));
        copies.push_back(e.copies);
    }

    auto build_input = [&]() {
        ArrangePolygons xs;
        for (size_t i = 0; i < shapes.size(); ++i)
            for (int c = 0; c < copies[i]; ++c)
                xs.push_back(make_ap(shapes[i]));
        scatter(xs, 0xD1CE0A11u, 400.0);
        return xs;
    };

    BoundingBox bed = make_bed_mm(256.0, 256.0);
    ArrangeParams params = benchmark_params();
    auto run = [&](ArrangePolygons &xs) {
        BitmapNester::arrange(xs, ArrangePolygons{}, bed, params);
    };

    ArrangePolygons items = build_input();
    run(items);

    int max_bed  = test_utils::max_bed_idx(items);
    int overflow = test_utils::overflow_piece_count(items);
    int total    = (int)items.size();

    UNSCOPED_INFO("Scenario:  real mix");
    UNSCOPED_INFO("  bed    = 256 x 256 mm");
    UNSCOPED_INFO("  pieces = " << total);
    UNSCOPED_INFO("  max bed  = " << max_bed);
    UNSCOPED_INFO("  overflow = " << overflow);

    // Primary gate: no spillover at all.  This is a realistic workload
    // with bed area in reserve — any spillover is a regression.
    REQUIRE(max_bed <= 0);
    REQUIRE(overflow == 0);
    REQUIRE(test_utils::no_overlap(items));

    // Determinism: byte-identical second run.
    ArrangePolygons rerun = build_input();
    REQUIRE(test_utils::deterministic_rerun(rerun, run));
}
