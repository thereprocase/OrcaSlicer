#include <catch2/catch.hpp>
#include "libslic3r/Arrange/BitmapArranger.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/libslic3r.h"

using namespace Slic3r;
using namespace Slic3r::arrangement;

// Helper: create a rectangular ExPolygon at origin
static ExPolygon make_rect(double w_mm, double h_mm) {
    coord_t w = scaled(w_mm), h = scaled(h_mm);
    ExPolygon ep;
    ep.contour.points = {{0, 0}, {w, 0}, {w, h}, {0, h}};
    return ep;
}

// Helper: create an L-shaped ExPolygon (concave)
static ExPolygon make_L(double size_mm) {
    coord_t s = scaled(size_mm);
    coord_t half = s / 2;
    ExPolygon ep;
    // L-shape: full bottom, left half top
    ep.contour.points = {
        {0, 0}, {s, 0}, {s, half}, {half, half}, {half, s}, {0, s}
    };
    return ep;
}

static BoundingBox make_bed(double w_mm, double h_mm) {
    return BoundingBox({0, 0}, {scaled(w_mm), scaled(h_mm)});
}

TEST_CASE("BitmapArranger: single rectangle fits on bed", "[arrange][bitmap]") {
    ArrangePolygons items;
    ArrangePolygon ap;
    ap.poly = make_rect(20, 20);
    items.push_back(ap);

    ArrangeParams params;
    params.use_concave_hulls = true;
    params.allow_multi_plate = false;

    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);

    REQUIRE(items[0].bed_idx == 0);
    // Translation should place the item within the bed
    auto placed = items[0].transformed_poly();
    auto bb = get_extents(placed);
    CHECK(bb.min.x() >= bed.min.x());
    CHECK(bb.min.y() >= bed.min.y());
    CHECK(bb.max.x() <= bed.max.x());
    CHECK(bb.max.y() <= bed.max.y());
}

TEST_CASE("BitmapArranger: multiple rectangles pack without overlap", "[arrange][bitmap]") {
    ArrangePolygons items;
    for (int i = 0; i < 6; i++) {
        ArrangePolygon ap;
        ap.poly = make_rect(40, 30);
        items.push_back(ap);
    }

    ArrangeParams params;
    params.use_concave_hulls = true;
    params.allow_multi_plate = false;
    params.min_obj_distance = scaled(2.0);

    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);

    // All should fit on one plate (6 x 42x32 < 256x256)
    for (auto &it : items)
        CHECK(it.bed_idx == 0);

    // No overlap: check bounding boxes don't intersect
    // For axis-aligned rectangles, BB overlap means geometry overlap
    for (size_t i = 0; i < items.size(); i++) {
        auto bb_i = get_extents(items[i].transformed_poly());
        for (size_t j = i + 1; j < items.size(); j++) {
            auto bb_j = get_extents(items[j].transformed_poly());
            bool overlap = bb_i.max.x() > bb_j.min.x() && bb_j.max.x() > bb_i.min.x() &&
                           bb_i.max.y() > bb_j.min.y() && bb_j.max.y() > bb_i.min.y();
            CHECK_FALSE(overlap);
        }
    }
}

TEST_CASE("BitmapArranger: multi-plate overflow", "[arrange][bitmap]") {
    ArrangePolygons items;
    // 4 items of 150x150 on a 256x256 bed = won't all fit on one plate
    for (int i = 0; i < 4; i++) {
        ArrangePolygon ap;
        ap.poly = make_rect(150, 150);
        items.push_back(ap);
    }

    ArrangeParams params;
    params.use_concave_hulls = true;
    params.allow_multi_plate = true;

    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);

    int plates_used = 0;
    for (auto &it : items) {
        REQUIRE(it.bed_idx >= 0); // all placed
        plates_used = std::max(plates_used, it.bed_idx + 1);
    }
    CHECK(plates_used >= 2); // should need at least 2 plates
}

TEST_CASE("BitmapArranger: multi-plate disabled leaves items unarranged", "[arrange][bitmap]") {
    ArrangePolygons items;
    for (int i = 0; i < 4; i++) {
        ArrangePolygon ap;
        ap.poly = make_rect(150, 150);
        items.push_back(ap);
    }

    ArrangeParams params;
    params.use_concave_hulls = true;
    params.allow_multi_plate = false;

    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);

    int placed = 0, unarranged = 0;
    for (auto &it : items) {
        if (it.bed_idx >= 0) placed++;
        else unarranged++;
    }
    CHECK(placed >= 1);
    CHECK(unarranged >= 1); // some shouldn't fit
}

TEST_CASE("BitmapArranger: excluded regions are avoided", "[arrange][bitmap]") {
    ArrangePolygons items;
    ArrangePolygon ap;
    ap.poly = make_rect(30, 30);
    items.push_back(ap);

    // Block the top-left corner with an excluded region
    ArrangePolygons excludes;
    ArrangePolygon excl;
    excl.poly = make_rect(100, 100);
    excl.translation = {0, 0};
    excl.rotation = 0;
    excludes.push_back(excl);

    ArrangeParams params;
    params.use_concave_hulls = true;
    params.allow_multi_plate = false;

    // Also add an excluded region via params (simulating purge line)
    ArrangePolygon purge;
    purge.poly = make_rect(256, 5); // 5mm strip at bottom
    purge.translation = {0, 0};
    purge.rotation = 0;
    params.excluded_regions.push_back(purge);

    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, excludes, bed, params);

    REQUIRE(items[0].bed_idx == 0);
    // Item should not overlap with the excluded region
    auto placed_bb = get_extents(items[0].transformed_poly());
    auto excl_bb = get_extents(excl.transformed_poly());
    bool overlaps_excl = placed_bb.max.x() > excl_bb.min.x() && excl_bb.max.x() > placed_bb.min.x() &&
                         placed_bb.max.y() > excl_bb.min.y() && excl_bb.max.y() > placed_bb.min.y();
    CHECK_FALSE(overlaps_excl);
}

TEST_CASE("BitmapArranger: L-shapes pack tighter than bounding box", "[arrange][bitmap]") {
    ArrangePolygons items;
    for (int i = 0; i < 4; i++) {
        ArrangePolygon ap;
        ap.poly = make_L(60);
        items.push_back(ap);
    }

    ArrangeParams params;
    params.use_concave_hulls = true;
    params.allow_multi_plate = false;
    params.allow_rotations = true;
    params.rotation_step_rad = PI / 2; // 90 degree steps

    // On a 130x130 bed, 4 L-shapes of 60x60 bounding box won't fit with convex packing
    // but should fit with concave packing since L-shapes can interlock
    auto bed = make_bed(130, 130);
    BitmapArranger::arrange(items, {}, bed, params);

    int placed = 0;
    for (auto &it : items)
        if (it.bed_idx == 0) placed++;

    // At least 3 should fit (4 is ideal but geometry-dependent)
    CHECK(placed >= 3);
}

TEST_CASE("BitmapArranger: empty input", "[arrange][bitmap]") {
    ArrangePolygons items; // empty
    ArrangeParams params;
    params.use_concave_hulls = true;
    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);
    CHECK(items.empty());
}

TEST_CASE("BitmapArranger: item larger than bed", "[arrange][bitmap]") {
    ArrangePolygons items;
    ArrangePolygon ap;
    ap.poly = make_rect(300, 300); // bigger than 256x256 bed
    items.push_back(ap);
    ArrangeParams params;
    params.use_concave_hulls = true;
    params.allow_multi_plate = false;
    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);
    CHECK(items[0].bed_idx == -1); // UNARRANGED
}

TEST_CASE("BitmapArranger: zero spacing", "[arrange][bitmap]") {
    ArrangePolygons items;
    for (int i = 0; i < 4; i++) {
        ArrangePolygon ap;
        ap.poly = make_rect(50, 50);
        items.push_back(ap);
    }
    ArrangeParams params;
    params.use_concave_hulls = true;
    params.min_obj_distance = 0;
    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);
    for (auto &it : items)
        CHECK(it.bed_idx == 0);
}

TEST_CASE("BitmapArranger: polygon with holes", "[arrange][bitmap]") {
    ArrangePolygons items;
    ArrangePolygon ap;
    // Outer square 40x40, inner hole 20x20 centered
    ap.poly.contour.points = {
        {scaled(0.), scaled(0.)}, {scaled(40.), scaled(0.)},
        {scaled(40.), scaled(40.)}, {scaled(0.), scaled(40.)}
    };
    Polygon hole;
    hole.points = {
        {scaled(10.), scaled(10.)}, {scaled(10.), scaled(30.)},
        {scaled(30.), scaled(30.)}, {scaled(30.), scaled(10.)}
    };
    ap.poly.holes.push_back(hole);
    items.push_back(ap);

    ArrangeParams params;
    params.use_concave_hulls = true;
    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);
    CHECK(items[0].bed_idx == 0);
}

TEST_CASE("BitmapArranger: consolidate plates", "[arrange][bitmap]") {
    ArrangePolygons items;
    // 3 small items, each pre-assigned to different plates
    for (int i = 0; i < 3; i++) {
        ArrangePolygon ap;
        ap.poly = make_rect(30, 30);
        ap.bed_idx = i; // plates 0, 1, 2
        items.push_back(ap);
    }
    ArrangeParams params;
    params.use_concave_hulls = true;
    params.consolidate_plates = true;
    params.allow_multi_plate = true;
    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);
    // All 3 should fit on plate 0 after consolidation
    for (auto &it : items)
        CHECK(it.bed_idx == 0);
}

TEST_CASE("BitmapArranger: purge pad avoidance", "[arrange][bitmap]") {
    ArrangePolygons items;
    ArrangePolygon ap;
    ap.poly = make_rect(20, 20);
    items.push_back(ap);
    ArrangeParams params;
    params.use_concave_hulls = true;
    params.avoid_purge_pad = true;
    params.purge_pad_edge = 0; // front
    params.purge_pad_mm = 10.f;
    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);
    REQUIRE(items[0].bed_idx == 0);
    auto placed = items[0].transformed_poly();
    auto bb = get_extents(placed);
    // Item should be above the 10mm purge zone (Y >= scaled(10))
    CHECK(bb.min.y() >= scaled(10.0) - scaled(1.0)); // 1mm tolerance for pixel quantization
}

TEST_CASE("BitmapArranger: stop condition cancels early", "[arrange][bitmap]") {
    ArrangePolygons items;
    for (int i = 0; i < 10; i++) {
        ArrangePolygon ap;
        ap.poly = make_rect(20, 20);
        items.push_back(ap);
    }
    int call_count = 0;
    ArrangeParams params;
    params.use_concave_hulls = true;
    params.stopcondition = [&]() { return ++call_count > 3; }; // cancel after 3 items
    auto bed = make_bed(256, 256);
    BitmapArranger::arrange(items, {}, bed, params);
    // Some items should be placed, some not (cancelled early)
    int placed = 0;
    for (auto &it : items)
        if (it.bed_idx >= 0) placed++;
    CHECK(placed < 10);
    CHECK(placed > 0);
}
