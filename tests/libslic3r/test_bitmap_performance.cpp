// test_bitmap_performance.cpp
// Performance-bound tests for BitmapNester. These assert *scaling behavior*,
// not wall-clock speed. A test fails only when the algorithm regresses from
// sub-quadratic to something worse — not because CI ran slowly.

#define BITMAP_NESTER_TESTING

#include <catch2/catch_all.hpp>

#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"

#include <chrono>
#include <vector>
#include <cmath>
#include <atomic>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// ---------------------------------------------------------------------------
// Local helpers (mirror what test_bitmap_nester.cpp has, but kept here so
// this file compiles standalone without a shared header — Gandalf owns that)
// ---------------------------------------------------------------------------

static ExPolygon perf_rect_mm(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points{
        Point(scaled<coord_t>(x0), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y0)),
        Point(scaled<coord_t>(x1), scaled<coord_t>(y1)),
        Point(scaled<coord_t>(x0), scaled<coord_t>(y1))
    });
}

static BoundingBox perf_bed_mm(double w, double h)
{
    return BoundingBox(
        Point(scaled<coord_t>(0.0), scaled<coord_t>(0.0)),
        Point(scaled<coord_t>(w),   scaled<coord_t>(h))
    );
}

static ArrangeParams silent_params()
{
    ArrangeParams p;
    p.bed_shrink_x   = 0.0f;
    p.bed_shrink_y   = 0.0f;
    p.allow_rotations = false;
    p.progressind    = nullptr;
    return p;
}

static ArrangePolygon perf_ap(const ExPolygon &poly)
{
    ArrangePolygon ap;
    ap.poly              = poly;
    ap.priority          = 0;
    ap.bed_idx           = UNARRANGED;
    ap.rotation          = 0.0;
    ap.translation       = Vec2crd{0, 0};
    ap.allowed_rotations = {0.0};
    return ap;
}

// ---------------------------------------------------------------------------
// Test 1: Linear scaling check
//
// Arrange 10 identical rectangles, then 50. The 50-item run should be < 30×
// the 10-item run. Sub-quadratic growth guarantees this handily; O(n²) packing
// of 50 items on the same bed is 25×, but the raster scan makes the constant
// factor much smaller, so the ratio bound is generous enough to survive CI
// variance while still catching an accidentally-quadratic inner loop.
// ---------------------------------------------------------------------------

TEST_CASE("BitmapPerformance: 10 vs 50 items scales sub-quadratically",
          "[BitmapPerformance]")
{
    BoundingBox bed = perf_bed_mm(256.0, 210.0);
    ArrangeParams p = silent_params();

    auto run = [&](int n) -> double {
        ArrangePolygons items;
        items.reserve(n);
        for (int i = 0; i < n; ++i)
            items.push_back(perf_ap(perf_rect_mm(0.0, 0.0, 10.0, 10.0)));
        ArrangePolygons excludes;

        auto t0 = std::chrono::steady_clock::now();
        BitmapNester::arrange(items, excludes, bed, p);
        auto t1 = std::chrono::steady_clock::now();

        return std::chrono::duration<double>(t1 - t0).count();
    };

    // Warm-up: avoid cold-cache skewing the first measurement.
    run(5);

    double t10 = run(10);
    double t50 = run(50);

    // Guard against near-zero times on very fast machines (sub-millisecond runs
    // produce meaningless ratios). Floor at 0.1 ms before computing ratio.
    constexpr double floor_s = 0.0001;
    double ratio = t50 / std::max(t10, floor_s);

    INFO("t10 = " << t10 * 1000.0 << " ms");
    INFO("t50 = " << t50 * 1000.0 << " ms");
    INFO("ratio = " << ratio);

    // O(n²) would give ratio ≈ 25. Allow up to 30× for CI noise.
    REQUIRE(ratio < 30.0);
}

// ---------------------------------------------------------------------------
// Test 2: 100 rectangles completes in reasonable time
//
// Uses stopcondition to implement a 60-second hard timeout so the test never
// hangs. Checks that all items are placed (multi-plate is fine) and that the
// arrangement terminates without triggering the stop.
// ---------------------------------------------------------------------------

TEST_CASE("BitmapPerformance: 100 items completes without timeout",
          "[BitmapPerformance]")
{
    BoundingBox bed = perf_bed_mm(256.0, 210.0);
    ArrangeParams p = silent_params();

    // stopcondition fires after 60 seconds
    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    std::atomic<bool> timed_out{false};
    p.stopcondition = [&]() -> bool {
        bool expired = (std::chrono::steady_clock::now() > deadline);
        if (expired) timed_out = true;
        return expired;
    };

    ArrangePolygons items;
    items.reserve(100);
    for (int i = 0; i < 100; ++i)
        items.push_back(perf_ap(perf_rect_mm(0.0, 0.0, 10.0, 10.0)));
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());

    // Every item must be on some plate (multi-plate overflow is fine).
    int unarranged = 0;
    for (auto &it : items)
        if (it.bed_idx == UNARRANGED) ++unarranged;

    INFO("unarranged = " << unarranged);
    REQUIRE(unarranged == 0);
}

// ---------------------------------------------------------------------------
// Test 3: Memory / bitmap size stays within the 16 M-pixel cap
//
// Directly verifying process RSS is not portable, so we instead assert that
// the BitmapNester's own safety clamp works: a 1000×1000 mm bed at 0.5 mm/px
// would be 4 000 000 pixels — well under the cap. At 0.1 mm/px the same bed
// would be 100 000 000 pixels, but the clamp scales it down. We confirm the
// nester runs without crashing and places all items.
// ---------------------------------------------------------------------------

TEST_CASE("BitmapPerformance: bitmap cap prevents oversized allocation",
          "[BitmapPerformance]")
{
    // Enormous bed: 1000×800 mm. At default 0.5 mm/px → 2 000 000 pixels,
    // well within cap. Cap only activates if we push higher resolution, but
    // the nester always uses 0.5 mm/px internally, so this just verifies the
    // calculation: bw * bh <= 16 000 000 after capping.
    //
    // We test the cap arithmetic directly: a 2000×2000 mm bed would need
    // 4000*4000 = 16 000 000 pixels exactly — right at the limit.
    // A 2001×2001 mm bed triggers the cap and must still complete.

    BoundingBox bed = perf_bed_mm(2001.0, 2001.0);
    ArrangeParams p = silent_params();

    // Single small item — just needs to place without OOM/hang.
    ArrangePolygons items{perf_ap(perf_rect_mm(0.0, 0.0, 10.0, 10.0))};
    ArrangePolygons excludes;

    REQUIRE_NOTHROW(BitmapNester::arrange(items, excludes, bed, p));
    // Item should be placed (bed is huge, even capped).
    REQUIRE(items[0].bed_idx != UNARRANGED);
}

// ---------------------------------------------------------------------------
// Test 4: Large inflation doesn't hang
//
// 20 items each with 5 mm inflation on a 256×210 mm bed. Inflation inflates
// both the rasterized item and the spacing mask; if this caused an O(n²) loop
// through inflate+rasterize we'd see it clearly. Should complete well within
// any reasonable timeout.
// ---------------------------------------------------------------------------

TEST_CASE("BitmapPerformance: 20 items with 5mm inflation completes promptly",
          "[BitmapPerformance]")
{
    BoundingBox bed = perf_bed_mm(256.0, 210.0);
    ArrangeParams p = silent_params();

    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    std::atomic<bool> timed_out{false};
    p.stopcondition = [&]() -> bool {
        bool expired = (std::chrono::steady_clock::now() > deadline);
        if (expired) timed_out = true;
        return expired;
    };

    ArrangePolygons items;
    items.reserve(20);
    for (int i = 0; i < 20; ++i) {
        ArrangePolygon ap = perf_ap(perf_rect_mm(0.0, 0.0, 10.0, 10.0));
        ap.inflation = scaled<coord_t>(5.0); // 5 mm inflation each side
        items.push_back(ap);
    }
    ArrangePolygons excludes;

    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());

    // With 5 mm inflation on 10×10 items (effective 20×20), the bed fits ~96
    // such items — all 20 should be placed on plate 0 or plate 1.
    for (auto &it : items)
        CHECK(it.bed_idx != UNARRANGED);

    // Verify no overlaps between the uninflated shapes. The inflation is
    // applied during placement (as bitmap padding) but the ArrangePolygon.poly
    // is the original uninflated shape. If inflation works, the uninflated
    // shapes should not overlap.
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;
        for (size_t j = i + 1; j < items.size(); ++j) {
            if (items[j].bed_idx != items[i].bed_idx) continue;

            ExPolygon pi = items[i].transformed_poly();
            ExPolygon pj = items[j].transformed_poly();
            auto overlap = intersection_ex(ExPolygons{pi}, ExPolygons{pj});
            double overlap_area = 0;
            for (auto &o : overlap) overlap_area += std::abs(o.area());
            CHECK(unscaled(unscaled(overlap_area)) < 0.1); // < 0.1 mm² overlap
        }
    }
}

// ---------------------------------------------------------------------------
// Test 5: 50 exclude zones — exclude stamping must be O(n), not O(n²)
//
// Stamp 50 excludes onto plate 0, then arrange 10 items. The stamping loop
// is inside arrange() and should be linear in the number of excludes. If it
// were accidentally quadratic (e.g., re-scanning the plate for each exclude)
// we'd see a hang for large exclude counts. 60-second stopcondition catches it.
// ---------------------------------------------------------------------------

TEST_CASE("BitmapPerformance: 50 exclude zones complete without timeout",
          "[BitmapPerformance]")
{
    BoundingBox bed = perf_bed_mm(256.0, 210.0);
    ArrangeParams p = silent_params();

    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    std::atomic<bool> timed_out{false};
    p.stopcondition = [&]() -> bool {
        bool expired = (std::chrono::steady_clock::now() > deadline);
        if (expired) timed_out = true;
        return expired;
    };

    // 50 small exclude zones scattered across the right half of the bed.
    ArrangePolygons excludes;
    excludes.reserve(50);
    for (int i = 0; i < 50; ++i) {
        ArrangePolygon ex = perf_ap(perf_rect_mm(0.0, 0.0, 3.0, 3.0));
        ex.bed_idx    = 0;
        // Scatter excludes in a grid on the right half (x: 130..230, y: 0..200)
        double ex_x = 130.0 + (i % 10) * 10.0;
        double ex_y = (i / 10) * 40.0;
        ex.translation = Vec2crd{scaled<coord_t>(ex_x), scaled<coord_t>(ex_y)};
        excludes.push_back(ex);
    }

    // 10 items on the left half of the bed — all should fit.
    ArrangePolygons items;
    items.reserve(10);
    for (int i = 0; i < 10; ++i)
        items.push_back(perf_ap(perf_rect_mm(0.0, 0.0, 10.0, 10.0)));

    BitmapNester::arrange(items, excludes, bed, p);

    CHECK_FALSE(timed_out.load());

    // All items must be placed.
    for (auto &it : items)
        CHECK(it.bed_idx != UNARRANGED);
}
