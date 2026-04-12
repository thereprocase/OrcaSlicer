#pragma once

// Shared invariant-checking utilities for BitmapNester tests.
// Include after ClipperUtils.hpp and BitmapNester.hpp.

#include "libslic3r/Arrange.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Geometry/ConvexHull.hpp"

#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>

namespace Slic3r { namespace arrangement { namespace test_utils {

// ─── Quality metrics adopted by the Three Seers vote 2026-04-11 ────
// See docs/PACKING_METRICS.md for the full deliberation. The seers
// rejected bbox density, convex hull compactness, adjacency, and
// libnest2d head-to-head. They unanimously adopted plate-count
// regression and determinism, plus Frodo's overflow piece count.
//
// Helpers below are intentionally minimal — they read only the
// observable output (items[].bed_idx) and never touch internal
// nester state.

// Highest plate index used by any placed item. Returns -1 if no
// item is placed (all UNARRANGED). For a single-plate fixture this
// returns 0; for a fixture that overflowed by one plate, returns 1.
inline int max_bed_idx(const ArrangePolygons &items)
{
    int m = -1;
    for (const auto &it : items)
        if (it.bed_idx != UNARRANGED && it.bed_idx > m) m = it.bed_idx;
    return m;
}

// Number of items placed on a plate other than plate 0. The "regret"
// metric: how many pieces wish they were on plate 0. Frodo's pick.
inline int overflow_piece_count(const ArrangePolygons &items)
{
    int n = 0;
    for (const auto &it : items)
        if (it.bed_idx != UNARRANGED && it.bed_idx > 0) ++n;
    return n;
}

// Sum of cluster bounding-box PERIMETER across all plates, in mm.
// A direct proxy for worst-case print-head XY travel: a compact
// cluster has shorter perimeter, so the head's coverage envelope is
// smaller per plate. Two layouts with the same plate count and zero
// overflow can have very different perimeters (a 100x100 cluster has
// perimeter 400 mm; a 400x25 cluster of the same area has perimeter
// 850 mm). The cluster-bbox-growth primary scorer already optimizes
// for perimeter implicitly; this helper lets tests assert bounds
// directly so an algorithm change that trades off tighter-but-longer
// layouts can be detected.
//
// UNARRANGED items are skipped. Perimeter is computed per-plate from
// the axis-aligned bounding box of all placed items' transformed_poly.
// Returns 0.0 for an empty arrangement.
inline double cluster_bbox_perimeter_mm(const ArrangePolygons &items)
{
    // Gather items per plate, compute per-plate union bbox, sum perimeters.
    std::map<int, BoundingBox> per_plate;
    for (const auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        ExPolygon tp = it.transformed_poly();
        BoundingBox bb = get_extents(tp);
        auto it_map = per_plate.find(it.bed_idx);
        if (it_map == per_plate.end()) {
            per_plate.emplace(it.bed_idx, bb);
        } else {
            it_map->second.merge(bb);
        }
    }
    double total = 0.0;
    for (const auto &kv : per_plate) {
        const BoundingBox &bb = kv.second;
        double w = unscaled<double>(bb.size().x());
        double h = unscaled<double>(bb.size().y());
        total += 2.0 * (w + h);
    }
    return total;
}

// Compactness ratio: actual cluster area divided by the minimum
// enclosing axis-aligned rectangle area. 1.0 = cluster fills its
// bbox exactly; 0.5 = cluster wastes half its bbox. Higher is better
// for tight packing. Computed per-plate and averaged by placed-item
// count so a plate with a few tight items and a plate with many
// tight items both score near 1.
//
// "Actual area" here is the sum of item transformed_poly areas (sum
// of per-item footprints), NOT the unioned area. For non-overlapping
// items these are the same; for overlapping-by-rounding items it's
// a slight overcount but still comparable between layouts.
inline double cluster_compactness(const ArrangePolygons &items)
{
    struct PlateAgg {
        double area_mm2 = 0.0;
        BoundingBox bb;
        bool bb_init = false;
        int count = 0;
    };
    std::map<int, PlateAgg> per_plate;
    for (const auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        ExPolygon tp = it.transformed_poly();
        BoundingBox bb = get_extents(tp);
        double a = unscaled<double>(unscaled<double>(std::abs(tp.area())));
        auto &agg = per_plate[it.bed_idx];
        agg.area_mm2 += a;
        if (!agg.bb_init) { agg.bb = bb; agg.bb_init = true; }
        else               agg.bb.merge(bb);
        agg.count++;
    }
    if (per_plate.empty()) return 0.0;
    double weighted_sum = 0.0;
    int total_count = 0;
    for (const auto &kv : per_plate) {
        const PlateAgg &agg = kv.second;
        double w = unscaled<double>(agg.bb.size().x());
        double h = unscaled<double>(agg.bb.size().y());
        double bbox_area = w * h;
        if (bbox_area <= 0.0) continue;
        double ratio = agg.area_mm2 / bbox_area;
        if (ratio > 1.0) ratio = 1.0;  // clamp rounding overcounts
        weighted_sum += ratio * agg.count;
        total_count  += agg.count;
    }
    if (total_count == 0) return 0.0;
    return weighted_sum / (double)total_count;
}

// ─── Hull-based cluster metrics (C2 M2.1) ─────────────────────────────
//
// The bbox-based helpers above treat every cluster as its axis-
// aligned bounding rectangle. That's convenient but loses information
// whenever the actual pack is more compact than its bbox — a diagonal
// strip, a round cluster, or any interlocked concave layout. These
// helpers compute the convex HULL of all placed item vertices per
// plate, then measure the hull instead of the bbox. Hull perimeter
// is a much better proxy for print-head travel than bbox perimeter
// because it hugs the actual cluster outline.
//
// Policy: hull is computed across ALL vertices of ALL placed items
// on a given plate. Holes are ignored (the hull of a shape with a
// hole is the hull of its outer contour). For items rotated at an
// angle, the transformed_poly() is used so the hull reflects the
// actual placed geometry.
//
// Returns 0.0 for empty arrangements.

// Convex hull perimeter across all placed items, summed per-plate.
inline double cluster_hull_perimeter_mm(const ArrangePolygons &items)
{
    std::map<int, Points> per_plate_pts;
    for (const auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        ExPolygon tp = it.transformed_poly();
        auto &bucket = per_plate_pts[it.bed_idx];
        for (const Point &p : tp.contour.points) bucket.push_back(p);
    }
    double total = 0.0;
    for (auto &kv : per_plate_pts) {
        Points &pts = kv.second;
        if (pts.size() < 3) continue;
        Polygon hull = Geometry::convex_hull(pts);
        if (hull.points.size() < 3) continue;
        total += unscaled<double>(hull.length());
    }
    return total;
}

// Convex hull area across all placed items, summed per-plate. Useful
// as a denominator in "packing efficiency vs hull" ratios — the
// inverse of "how much air is inside the cluster's tightest
// wrapping". Differs from cluster bbox area whenever the cluster
// shape is not a rectangle (almost always).
inline double cluster_hull_area_mm2(const ArrangePolygons &items)
{
    std::map<int, Points> per_plate_pts;
    for (const auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        ExPolygon tp = it.transformed_poly();
        auto &bucket = per_plate_pts[it.bed_idx];
        for (const Point &p : tp.contour.points) bucket.push_back(p);
    }
    double total = 0.0;
    for (auto &kv : per_plate_pts) {
        Points &pts = kv.second;
        if (pts.size() < 3) continue;
        Polygon hull = Geometry::convex_hull(pts);
        if (hull.points.size() < 3) continue;
        total += unscaled<double>(unscaled<double>(std::abs(hull.area())));
    }
    return total;
}

// Hull-based compactness ratio: sum of placed silhouette areas
// divided by hull area, per plate, weighted by item count.
// Equals 1.0 when the cluster's outline IS the union of its items
// (impossible in practice — hull always >= union area). For
// concave shapes the hull is typically 5-20% larger than the
// silhouette union, so this ratio lands around 0.8-0.95 for tight
// packs. Lower = air trapped inside the hull. Higher = tight pack.
inline double cluster_hull_compactness(const ArrangePolygons &items)
{
    struct PlateAgg {
        double silhouette_area_mm2 = 0.0;
        Points hull_pts;
        int count = 0;
    };
    std::map<int, PlateAgg> per_plate;
    for (const auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        ExPolygon tp = it.transformed_poly();
        double a = unscaled<double>(unscaled<double>(std::abs(tp.area())));
        auto &agg = per_plate[it.bed_idx];
        agg.silhouette_area_mm2 += a;
        for (const Point &p : tp.contour.points) agg.hull_pts.push_back(p);
        agg.count++;
    }
    if (per_plate.empty()) return 0.0;
    double weighted_sum = 0.0;
    int total_count = 0;
    for (auto &kv : per_plate) {
        PlateAgg &agg = kv.second;
        if (agg.hull_pts.size() < 3) continue;
        Polygon hull = Geometry::convex_hull(agg.hull_pts);
        if (hull.points.size() < 3) continue;
        double hull_area = unscaled<double>(unscaled<double>(std::abs(hull.area())));
        if (hull_area <= 0.0) continue;
        double ratio = agg.silhouette_area_mm2 / hull_area;
        if (ratio > 1.0) ratio = 1.0;  // clamp rounding overcount
        weighted_sum += ratio * agg.count;
        total_count  += agg.count;
    }
    if (total_count == 0) return 0.0;
    return weighted_sum / (double)total_count;
}

// Run an arrange call twice on a fresh copy of the input items each
// time, then verify the two outputs are byte-identical on the
// observable fields (bed_idx, rotation, translation, itemid).
// Returns true if deterministic, false if any field diverges.
//
// Determinism is the cheapest possible defense against state-machine
// bugs in the consolidation pass — if a read-before-write or
// uninitialized cell sneaks in, the second run diverges. Cost: one
// extra arrange call per fixture.
template <class ArrangeFn>
inline bool deterministic_rerun(const ArrangePolygons &input,
                                ArrangeFn &&run)
{
    ArrangePolygons a = input;
    ArrangePolygons b = input;
    run(a);
    run(b);
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].bed_idx != b[i].bed_idx) return false;
        if (a[i].itemid != b[i].itemid) return false;
        if (std::abs(a[i].rotation - b[i].rotation) > 1e-9) return false;
        if (a[i].translation.x() != b[i].translation.x()) return false;
        if (a[i].translation.y() != b[i].translation.y()) return false;
    }
    return true;
}

// Check every placed item's transformed polygon lies inside the bed bounds.
// Returns false if any placed item protrudes outside the bed.
inline bool all_within_bounds(const ArrangePolygons &items, const BoundingBox &bed)
{
    for (auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb = get_extents(placed);
        // Allow 1mm tolerance for raster quantization
        coord_t tol = scaled<coord_t>(1.0);
        if (bb.min.x() < bed.min.x() - tol) return false;
        if (bb.min.y() < bed.min.y() - tol) return false;
        if (bb.max.x() > bed.max.x() + tol) return false;
        if (bb.max.y() > bed.max.y() + tol) return false;
    }
    return true;
}

// Check bed_idx values are contiguous: 0,1,2,... with no gaps.
// Ignores UNARRANGED items. Returns false on gaps like {0, 2}.
inline bool no_bed_gaps(const ArrangePolygons &items)
{
    std::set<int> beds;
    for (auto &it : items) {
        if (it.bed_idx != UNARRANGED)
            beds.insert(it.bed_idx);
    }
    if (beds.empty()) return true;
    int expected = 0;
    for (int idx : beds) {
        if (idx != expected) return false;
        ++expected;
    }
    return true;
}

// Check no item is UNARRANGED. Returns false if any item was not placed.
inline bool all_placed(const ArrangePolygons &items)
{
    for (auto &it : items) {
        if (it.bed_idx == UNARRANGED) return false;
    }
    return true;
}

// Check no two placed items on the same bed overlap.
// Allows a tiny sliver (< 0.01 mm^2) for numerical noise.
inline bool no_overlap(const ArrangePolygons &items)
{
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;
        for (size_t j = i + 1; j < items.size(); ++j) {
            if (items[j].bed_idx != items[i].bed_idx) continue;

            ExPolygon pi = items[i].transformed_poly();
            ExPolygon pj = items[j].transformed_poly();

            ExPolygons inter = intersection_ex(ExPolygons{pi}, ExPolygons{pj});
            for (auto &seg : inter) {
                if (std::abs(seg.area()) > scaled<double>(0.01) * scaled<double>(0.01))
                    return false;
            }
        }
    }
    return true;
}

// ─── Visual render: dump placement to ASCII art ──────────────────
//
// Renders placed items onto a bed-sized grid as a compact text file.
// Each cell (5mm × 5mm at default) shows a character: '.' for empty,
// 'A'-'Z' for items (cycling), '#' for overlap, 'X' for outside bed.
// The bed boundary is drawn with '+', '-', '|'.
//
// Output is ~40×40 characters for a 200×200 bed at 5mm/cell — fits
// in a single screen and costs ~200 tokens for Claude to review.
//
// Usage:
//   test_utils::dump_placement_ascii(items, bed, "test_name.txt");

inline void dump_placement_ascii(const ArrangePolygons& items,
                                  const BoundingBox& bed,
                                  const std::string& filename,
                                  double cell_mm = 5.0)
{
    double bed_w = unscaled<double>(bed.size().x());
    double bed_h = unscaled<double>(bed.size().y());
    int cols = (int)(bed_w / cell_mm) + 2;  // +2 for border
    int rows = (int)(bed_h / cell_mm) + 2;
    if (cols > 200 || rows > 200) return;

    // Grid: 0 = empty, 1-26 = item A-Z, -1 = overlap
    std::vector<int> grid(cols * rows, 0);

    double bed_min_x = unscaled<double>(bed.min.x());
    double bed_min_y = unscaled<double>(bed.min.y());

    for (std::size_t i = 0; i < items.size(); ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;

        ExPolygon poly = items[i].poly;
        if (items[i].rotation != 0.0) poly.rotate(items[i].rotation);
        poly.translate(items[i].translation.x(), items[i].translation.y());

        int label = (int)(i % 26) + 1;

        for (int r = 0; r < rows - 2; ++r) {
            for (int c = 0; c < cols - 2; ++c) {
                double wx = bed_min_x + (c + 0.5) * cell_mm;
                double wy = bed_min_y + (r + 0.5) * cell_mm;
                Point pt(scaled<coord_t>(wx), scaled<coord_t>(wy));
                if (poly.contains(pt)) {
                    int gi = (r + 1) * cols + (c + 1);
                    if (grid[gi] != 0 && grid[gi] != label)
                        grid[gi] = -1;  // overlap
                    else
                        grid[gi] = label;
                }
            }
        }
    }

    // Write text. Top row = Y max (flip for natural orientation).
    FILE* f = fopen(filename.c_str(), "w");
    if (!f) return;

    // Header.
    fprintf(f, "# C2 placement render: %dx%d mm bed, %.0f mm/cell\n",
            (int)bed_w, (int)bed_h, cell_mm);
    fprintf(f, "# Items: %d placed\n", (int)items.size());

    // Top border.
    fprintf(f, "+");
    for (int c = 0; c < cols - 2; ++c) fprintf(f, "-");
    fprintf(f, "+\n");

    // Grid rows (top-down = high Y first).
    for (int r = rows - 3; r >= 0; --r) {
        fprintf(f, "|");
        for (int c = 0; c < cols - 2; ++c) {
            int v = grid[(r + 1) * cols + (c + 1)];
            if (v == 0)       fprintf(f, ".");
            else if (v == -1) fprintf(f, "#");
            else              fprintf(f, "%c", 'A' + (v - 1));
        }
        fprintf(f, "|\n");
    }

    // Bottom border.
    fprintf(f, "+");
    for (int c = 0; c < cols - 2; ++c) fprintf(f, "-");
    fprintf(f, "+\n");

    fclose(f);
}

}}} // namespace Slic3r::arrangement::test_utils
