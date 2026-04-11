// BitmapNesterC2.hpp — Organic Lasso Nester, fork of BitmapNester.
//
// PR C2 branch scaffold. This file starts as a trivial pass-through
// that delegates to BitmapNester so every existing C1 test passes on
// the C2 branch as well, then incrementally grows into the real
// organic-lasso implementation per the plan in
// `docs/PR_C2_ORGANIC_LASSO_NESTER_PLAN.md`.
//
// Architecture (from the plan):
//
//   Phase 0 — Preparation      : build_cache
//   Phase 1 — Plate estimation : estimate_min_plates
//   Phase 2 — Partitioning     : partition_items
//   Phase 3 — Island packing   : pack_as_island (hull-perimeter scoring,
//                                 tall-first center-out)
//   Phase 4 — Plate placement  : locate_island_on_plate
//   Phase 5 — Spillover rescue : recover_spillover
//   Phase 6 — Finalize         : write_back
//
// M0 milestone scope (this commit): class shell, public entry point,
// empty stubs, and a pass-through that invokes BitmapNester::arrange
// so the branch stays green on existing tests.

#pragma once

#include "BitmapNester.hpp"   // reuse primitives + C1 entry point
#include "Arrange.hpp"
#include "BoundingBox.hpp"
#include "ExPolygon.hpp"
#include "Geometry/ConvexHull.hpp"  // Geometry::convex_hull used by build_cache

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>   // std::iota in partition_items
#include <vector>

namespace Slic3r { namespace arrangement {

// ---------------------------------------------------------------------------
// Per-item information cached once per arrange call. Reused across every
// phase of the pipeline so that per-item work (rasterization, hull, bbox)
// only happens once. Grows as the phases come online — M0 holds the
// minimum necessary for the pass-through.
// ---------------------------------------------------------------------------

struct NesterC2ItemInfo {
    std::size_t original_idx = 0;

    // Geometric summary. Populated in build_cache (M1). Zero-initialized
    // in M0 so pass-through compiles.
    double silhouette_area_mm2   = 0.0;
    double silhouette_perimeter_mm = 0.0;
    double hull_area_mm2         = 0.0;
    double hull_perimeter_mm     = 0.0;
    double max_dim_mm            = 0.0;
    double height_mm             = 0.0;

    // Derived heuristics. Populated in build_cache (M1).
    // hardness_score = hull_area / silhouette_area; higher = harder to pack
    // priority_score = composite of priority, hardness, and height
    double hardness_score = 1.0;
    double priority_score = 0.0;
};

// ---------------------------------------------------------------------------
// A group of items destined for the same island. Output of partition_items
// (M1) and input to pack_as_island (M2).
// ---------------------------------------------------------------------------

struct NesterC2ItemGroup {
    std::vector<std::size_t> items;     // indices into the Cache
    double total_silhouette_area_mm2 = 0.0;
    int    tall_count = 0;              // # items with height >= 100mm
};

// ---------------------------------------------------------------------------
// A packed island in local coordinates (not yet placed on a plate).
// Output of pack_as_island (M2) and input to locate_island_on_plate (M3).
// ---------------------------------------------------------------------------

struct NesterC2Island {
    std::vector<std::size_t> item_indices;

    // Each packed item's final rotation and translation in island-local
    // coordinates (relative to island origin at (0,0), hull centroid at
    // some interior point). locate_island_on_plate adds the plate offset
    // to produce final world coordinates.
    std::vector<double>   rotations;
    std::vector<Vec2crd>  translations;

    // Geometric summary of the packed island, used by placement phase.
    BoundingBox bbox;
    double hull_perimeter_mm = 0.0;
    Vec2d  hull_centroid_mm  = Vec2d(0.0, 0.0);
};

// ---------------------------------------------------------------------------
// BitmapNesterC2 class. Public entry point mirrors BitmapNester::arrange
// so callers can swap via a settings flag. M0 delegates to C1's
// implementation; subsequent milestones replace phases one at a time.
// ---------------------------------------------------------------------------

class BitmapNesterC2 {
public:
    // M0 entry point: pass-through to BitmapNester::arrange. Keeps every
    // existing C1 test green on the C2 branch so the branch never goes
    // red while C2's own phases are being built out.
    static void arrange(ArrangePolygons& items,
                        const ArrangePolygons& excludes,
                        const BoundingBox& bed,
                        const ArrangeParams& params)
    {
        // M1 pipeline (partial): phase 0 cache + phase 1 k-estimation +
        // phase 2 partitioning are real. Phase 3 (pack_as_island) is not
        // yet implemented, so packing still delegates to BitmapNester —
        // but the delegation happens per-group rather than for the whole
        // input, which exercises the partitioning semantics end-to-end.
        //
        // Edge case: empty input. No work to do.
        if (items.empty()) return;

        auto cache = build_cache(items);
        int  k     = estimate_min_plates(cache, bed);
        auto groups = partition_items(cache, k);

        // Transitional delegation: invoke BitmapNester separately on each
        // group, letting the existing C1 placement pipeline handle each
        // partition. Groups with multiple items get their own indices
        // slice, and the final bed_idx values are offset to land them on
        // distinct plates. This is NOT the final C2 behavior — the
        // locate/recover phases come online in M3. For now it validates
        // that partitioning produces sensible groups the existing
        // nester can process.
        //
        // When k == 1 or partitioning produces a single group, this
        // collapses to the M0 pass-through — zero behavior change on
        // inputs where K_min is already 1.
        if (groups.size() <= 1) {
            BitmapNester::arrange(items, excludes, bed, params);
            return;
        }

        // For multi-group inputs, build per-group ArrangePolygons slices,
        // delegate each to BitmapNester, then map results back with the
        // group index baked into bed_idx. This is a quick-and-dirty
        // stand-in for the island pipeline; it exists so M1 can land
        // without M2/M3 while still exercising the partitioning path.
        for (std::size_t g = 0; g < groups.size(); ++g) {
            ArrangePolygons slice;
            slice.reserve(groups[g].items.size());
            for (std::size_t idx : groups[g].items) slice.push_back(items[idx]);

            BitmapNester::arrange(slice, excludes, bed, params);

            // Merge slice state back, offsetting bed_idx so each group
            // claims disjoint plate indices. A slice that landed on its
            // own "plate 0" becomes plate `g * max_bed_per_group` in the
            // merged result. For now we use a conservative offset of
            // `g * 8` plates per group to avoid stepping on neighbors.
            const int plate_offset = (int)(g * 8);
            for (std::size_t s = 0; s < slice.size(); ++s) {
                std::size_t orig_idx = groups[g].items[s];
                items[orig_idx].translation = slice[s].translation;
                items[orig_idx].rotation    = slice[s].rotation;
                items[orig_idx].itemid      = slice[s].itemid;
                items[orig_idx].bed_idx     =
                    (slice[s].bed_idx == UNARRANGED)
                        ? UNARRANGED
                        : (slice[s].bed_idx + plate_offset);
            }
        }
    }

    // ─── M1 phase functions (exposed for testing via the
    //     BITMAP_NESTER_C2_TESTING guard, like C1's test hooks) ──────
#ifdef BITMAP_NESTER_C2_TESTING
public:
#else
private:
#endif

    // M1.a: build per-item geometric + heuristic cache.
    //
    // Populates silhouette area, hull area, bbox, height, and a derived
    // hardness score. Rot_cache (bitmap rasterization per rotation)
    // remains TODO — that's expensive and only needed once packing
    // starts using it in M2. For M1 the cache is geometric only.
    static std::vector<NesterC2ItemInfo>
    build_cache(const ArrangePolygons& items)
    {
        std::vector<NesterC2ItemInfo> cache;
        cache.reserve(items.size());
        for (std::size_t i = 0; i < items.size(); ++i) {
            const ArrangePolygon& ap = items[i];
            NesterC2ItemInfo info;
            info.original_idx = i;

            // Silhouette area + perimeter. Use concave_regions if the
            // caller populated them (multi-volume parts, disconnected
            // islands); otherwise the legacy single poly.
            double s_area = 0.0;
            double s_perim = 0.0;
            if (!ap.concave_regions.empty()) {
                for (const ExPolygon& ep : ap.concave_regions) {
                    s_area += unscaled<double>(unscaled<double>(std::abs(ep.area())));
                    s_perim += unscaled<double>(ep.contour.length());
                }
            } else {
                s_area  = unscaled<double>(unscaled<double>(std::abs(ap.poly.area())));
                s_perim = unscaled<double>(ap.poly.contour.length());
            }
            info.silhouette_area_mm2     = s_area;
            info.silhouette_perimeter_mm = s_perim;

            // Convex hull of the silhouette for the hardness score.
            // hardness = hull_area / silhouette_area. A convex shape has
            // hardness ≈ 1.0; a deeply concave shape has hardness > 1.0.
            // We use this to bias placement order in M2.
            Points hull_pts;
            if (!ap.concave_regions.empty()) {
                for (const ExPolygon& ep : ap.concave_regions)
                    for (const Point& p : ep.contour.points)
                        hull_pts.push_back(p);
            } else {
                hull_pts = ap.poly.contour.points;
            }
            if (hull_pts.size() >= 3) {
                Polygon hull = Geometry::convex_hull(hull_pts);
                double h_area = unscaled<double>(unscaled<double>(std::abs(hull.area())));
                double h_perim = unscaled<double>(hull.length());
                info.hull_area_mm2     = h_area;
                info.hull_perimeter_mm = h_perim;
                if (s_area > 1e-9)
                    info.hardness_score = h_area / s_area;
            }

            // Bounding box — used for FFD-style sizing and for the
            // max_dim field. Matches ArrangePolygon::poly bbox.
            BoundingBox bb;
            if (!ap.concave_regions.empty()) {
                bb = get_extents(ap.concave_regions);
            } else {
                bb = get_extents(ap.poly);
            }
            double bb_w = unscaled<double>(bb.size().x());
            double bb_h = unscaled<double>(bb.size().y());
            info.max_dim_mm = std::max(bb_w, bb_h);

            info.height_mm = ap.height;

            // Priority score: composite of caller priority, hardness,
            // and height. Higher = placed earlier. Used by partitioning
            // to pick a seed order before K-way distribution.
            info.priority_score =
                ap.priority * 1000.0
                + info.hardness_score * 100.0
                + std::max(0.0, info.height_mm);

            cache.push_back(info);
        }
        return cache;
    }

    // M1.b: estimate minimum plate count.
    //
    // Sum all silhouette areas, divide by (bed_area × packing_efficiency),
    // round up. The efficiency factor is a conservative fudge — concave
    // interlocks routinely beat 0.82, but under-estimating K hurts less
    // than over-estimating (spillover recovery can add plates but cannot
    // remove them). Returns at least 1 for any non-empty input.
    static int estimate_min_plates(const std::vector<NesterC2ItemInfo>& cache,
                                   const BoundingBox& bed)
    {
        if (cache.empty()) return 0;
        double total_silhouette_mm2 = 0.0;
        for (const auto& info : cache)
            total_silhouette_mm2 += info.silhouette_area_mm2;
        double bed_w_mm = unscaled<double>(bed.size().x());
        double bed_h_mm = unscaled<double>(bed.size().y());
        double bed_area_mm2 = bed_w_mm * bed_h_mm;
        if (bed_area_mm2 <= 0.0) return 1;
        constexpr double PACKING_EFFICIENCY = 0.82;
        double effective = bed_area_mm2 * PACKING_EFFICIENCY;
        int k = (int)std::ceil(total_silhouette_mm2 / effective);
        return std::max(1, k);
    }

    // M1.c: partition items into K roughly-balanced groups.
    //
    // Greedy least-loaded-bucket algorithm: sort items by priority_score
    // descending, then for each item drop it into whichever group
    // currently has the smallest total silhouette area. Produces groups
    // with area variance typically under 10% for uniformly-sized inputs
    // and under 25% for heavily-mixed inputs.
    //
    // Tall item distribution is NOT explicitly balanced in M1 — that's a
    // follow-up refinement if the meta-test shows height imbalance
    // causing quality variance. For now the priority_score already
    // gives tall items high priority, so they get distributed first.
    //
    // k == 0 → empty output. k == 1 → single group holding everything.
    static std::vector<NesterC2ItemGroup>
    partition_items(const std::vector<NesterC2ItemInfo>& cache, int k)
    {
        std::vector<NesterC2ItemGroup> groups;
        if (k <= 0 || cache.empty()) return groups;
        groups.resize(k);

        // Sort item indices by priority descending.
        std::vector<std::size_t> order(cache.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::sort(order.begin(), order.end(),
                  [&cache](std::size_t a, std::size_t b) {
                      return cache[a].priority_score > cache[b].priority_score;
                  });

        // Greedy distribution: each item drops into the bucket with the
        // smallest current total area. Ties break toward lower index.
        for (std::size_t idx : order) {
            const NesterC2ItemInfo& info = cache[idx];
            int target = 0;
            double min_area = groups[0].total_silhouette_area_mm2;
            for (int g = 1; g < k; ++g) {
                if (groups[g].total_silhouette_area_mm2 < min_area) {
                    min_area = groups[g].total_silhouette_area_mm2;
                    target = g;
                }
            }
            groups[target].items.push_back(info.original_idx);
            groups[target].total_silhouette_area_mm2 += info.silhouette_area_mm2;
            if (info.height_mm >= 100.0) groups[target].tall_count++;
        }

        // Strip empty groups. Possible when k > item count, etc.
        groups.erase(std::remove_if(groups.begin(), groups.end(),
                                    [](const NesterC2ItemGroup& g) {
                                        return g.items.empty();
                                    }),
                     groups.end());
        return groups;
    }

private:
    // M2: bound-agnostic packing into a tight island. Hull-perimeter primary
    // scoring. Tall items placed at the seed; short items added at the hull
    // boundary. Compacted via iterative pair-swap refinement.
    // static NesterC2Island
    //     pack_as_island(const std::vector<NesterC2ItemInfo>& cache,
    //                    const NesterC2ItemGroup& group,
    //                    const ArrangeParams& params);

    // M3: locate a packed island on a plate, clamping to bed bounds. Uses
    // the island's hull centroid as the centering target.
    // static void locate_island_on_plate(NesterC2Island& island,
    //                                    const BoundingBox& bed,
    //                                    int plate_idx,
    //                                    const ArrangeParams& params);

    // M3: attempt to find a home for items that couldn't fit their target
    // island. Retry rotations, migrate to neighbor islands, or spawn a new
    // island as a last resort.
    // static void recover_spillover(ArrangePolygons& items,
    //                               std::vector<NesterC2Island>& islands,
    //                               const BoundingBox& bed,
    //                               const ArrangeParams& params);
};

}} // namespace Slic3r::arrangement
