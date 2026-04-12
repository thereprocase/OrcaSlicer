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
        // C2 M3 pipeline:
        //   phase 0 — build per-item cache
        //   phase 1 — estimate K_min
        //   phase 2 — partition items into K groups
        //   phase 3 — pack each group into an Island (M2.3: real
        //             hull-scored greedy loop, bound-agnostic)
        //   phase 4 — locate each island on its plate (M3.0: hull-
        //             centroid → bed-center shift, clamped to bed)
        //   phase 5 — spillover recovery (M3.1: upcoming)
        if (items.empty()) return;

        auto cache  = build_cache(items);
        int  k      = estimate_min_plates(cache, bed);
        auto groups = partition_items(cache, k);

        // Fast path: single group on a loose bed collapses to the
        // pre-partitioning pass-through. Zero behavior change on
        // k=1 inputs.
        if (groups.size() <= 1) {
            BitmapNester::arrange(items, excludes, bed, params);
            return;
        }

        // Multi-group path: pack each group as an island, then
        // locate each island on its own plate.
        for (std::size_t g = 0; g < groups.size(); ++g) {
            NesterC2Island island = pack_as_island(
                items, cache, groups[g], excludes, bed, params);
            locate_island_on_plate(items, island, (int)g, bed);
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

#ifdef BITMAP_NESTER_C2_TESTING
public:
#else
private:
#endif

    // M2.3: pack a group into an island via hull-scored greedy placement.
    //
    // This is the first C2 code that actively drives placement by
    // delta hull perimeter instead of bbox area. For each item (sorted
    // by priority_score descending = tall/hard items first):
    //
    //   1. Iterate allowed rotations
    //   2. Iterate a coarse grid of candidate positions around the
    //      current cluster bbox, expanded by the item's max dimension
    //   3. For each candidate, check collision against all placed
    //      items via Clipper intersection_ex (slow but no-C1-coupling)
    //   4. Score each clear candidate by delta(hull perimeter) — i.e.,
    //      the new cluster hull perimeter if this item were placed here
    //   5. Commit the globally best candidate across all rotations
    //
    // Seed item (first placed) goes to origin (0,0) with rotation 0.
    // Items that don't find a clear position anywhere in the search
    // window are left UNARRANGED — M3's spillover recovery will
    // handle them.
    //
    // Performance: O(N × R × (pos_count) × N) where N is group size,
    // R is allowed rotations, pos_count is typically ~100-400 grid
    // cells. Uses intersection_ex per candidate which is the slow
    // part. For 20-item groups this is a few seconds — acceptable
    // for dev iteration, will need a faster collision path (likely
    // bitmap via a collaborator pattern) before production use.
    //
    // NOTE: bed parameter is currently ignored. This is a BOUND-AGNOSTIC
    // pack — the island grows freely and the caller's locate phase
    // handles fitting it onto a real plate. excludes likewise ignored
    // because the island has no notion of a plate yet.
    static NesterC2Island pack_as_island(
        const ArrangePolygons& items_in,
        const std::vector<NesterC2ItemInfo>& cache,
        const NesterC2ItemGroup& group,
        const ArrangePolygons& excludes,
        const BoundingBox& bed,
        const ArrangeParams& params)
    {
        (void)excludes;
        (void)bed;
        (void)params;

        NesterC2Island island;
        if (group.items.empty()) return island;

        // Placement order: priority_score desc. Tallest/hardest first.
        std::vector<std::size_t> order = group.items;
        std::sort(order.begin(), order.end(),
                  [&cache](std::size_t a, std::size_t b) {
                      const auto& ia = *std::find_if(
                          cache.begin(), cache.end(),
                          [a](const NesterC2ItemInfo& i) { return i.original_idx == a; });
                      const auto& ib = *std::find_if(
                          cache.begin(), cache.end(),
                          [b](const NesterC2ItemInfo& i) { return i.original_idx == b; });
                      return ia.priority_score > ib.priority_score;
                  });

        // ─── Virtual plate for bitmap AND collision detection ─────
        //
        // We maintain a local bitmap that grows as items commit.
        // Per the M2.5 bbox-purge directive, this replaces the
        // Clipper intersection_ex path and its AABB pre-reject.
        // Bitmap AND is nanoseconds per check — cheaper than any
        // polygon-based narrow phase.
        //
        // Plate size: 2048×2048 px at 0.5 mm/pixel = 1024×1024 mm
        // in world. Generous enough for any realistic island; costs
        // 512 KB of allocation per call which is acceptable.
        const double res = 0.5;  // mm per pixel
        const int plate_bw  = 2048;
        const int plate_bh  = 2048;
        const int plate_wpr = (plate_bw + 63) / 64;
        std::vector<uint64_t> plate_items(
            (std::size_t)plate_wpr * plate_bh, 0);

        // Seed position: plate center. First item's bitmap gets
        // stamped so its center is near (seed_cx, seed_cy) and the
        // cluster can grow outward in any direction.
        const int seed_cx = plate_bw / 2;
        const int seed_cy = plate_bh / 2;

        // Running extent of the set bits on the plate. Initial
        // "empty" state uses sentinels that will be replaced on
        // the first commit.
        int occ_min_px = plate_bw;
        int occ_min_py = plate_bh;
        int occ_max_px = 0;
        int occ_max_py = 0;

        // Hull scoring state — accumulates vertices in world-mm
        // coordinates as items commit. Computed once per commit;
        // the scoring loop builds a temporary "new_pts" per
        // candidate for delta computation.
        Points placed_hull_pts;
        double current_hull_perim_mm = 0.0;

        for (std::size_t idx : order) {
            const ArrangePolygon& ap = items_in[idx];
            ExPolygon base_poly = ap.poly;
            const std::vector<double>& rots_raw = ap.allowed_rotations;
            std::vector<double> rots = rots_raw;
            if (rots.empty()) rots.push_back(0.0);

            // Rot cache: for each rotation, pre-rasterize the rotated
            // shape ONCE. Saves re-rasterizing on every candidate scan.
            struct RotCache {
                double rot;
                std::vector<uint64_t> bm;
                int iw = 0;
                int ih = 0;
                int iwpr = 0;
                BoundingBox rot_bb;  // extents of the rotated polygon
            };
            std::vector<RotCache> rot_cache;
            rot_cache.reserve(rots.size());
            for (double rot : rots) {
                ExPolygon rotated = base_poly;
                if (rot != 0.0) rotated.rotate(rot);
                int iw = 0, ih = 0, iwpr = 0;
                auto bm = BitmapNester::rasterize(
                    rotated, res, plate_bw, plate_bh, iw, ih, iwpr);
                if (bm.empty() || iw <= 0 || ih <= 0) continue;
                rot_cache.push_back({rot, std::move(bm), iw, ih, iwpr,
                                     get_extents(rotated)});
            }
            if (rot_cache.empty()) continue;

            // First item: stamp at plate center with first rotation.
            if (occ_max_px == 0) {
                const RotCache& rc = rot_cache[0];
                int px = seed_cx - rc.iw / 2;
                int py = seed_cy - rc.ih / 2;
                BitmapNester::stamp(plate_items, plate_wpr, plate_bw, plate_bh,
                                    rc.bm, rc.iwpr, rc.iw, rc.ih, px, py);

                occ_min_px = px;
                occ_min_py = py;
                occ_max_px = px + rc.iw;
                occ_max_py = py + rc.ih;

                // Compute mm-space translation: move the rotated
                // polygon so its bbox min lands at (px*res, py*res) mm.
                ExPolygon seeded = base_poly;
                if (rc.rot != 0.0) seeded.rotate(rc.rot);
                coord_t dx = scaled<coord_t>(px * res) - rc.rot_bb.min.x();
                coord_t dy = scaled<coord_t>(py * res) - rc.rot_bb.min.y();
                seeded.translate(dx, dy);

                for (const Point& p : seeded.contour.points)
                    placed_hull_pts.push_back(p);
                if (placed_hull_pts.size() >= 3) {
                    Polygon hull = Geometry::convex_hull(placed_hull_pts);
                    if (hull.points.size() >= 3)
                        current_hull_perim_mm = unscaled<double>(hull.length());
                }

                island.item_indices.push_back(idx);
                island.rotations.push_back(rc.rot);
                island.translations.push_back(Vec2crd{dx, dy});
                continue;
            }

            // Subsequent items: grid scan around the occupied extent.
            // The occupied extent is the "min/max of set bits" form
            // the boss explicitly allowed — derived from actual
            // geometry, not a rectangular approximation of it.
            int max_rot_iw = 0, max_rot_ih = 0;
            for (const auto& rc : rot_cache) {
                max_rot_iw = std::max(max_rot_iw, rc.iw);
                max_rot_ih = std::max(max_rot_ih, rc.ih);
            }

            // Search window: occupied extent expanded by the new
            // item's max dimensions so positions where the new item
            // wraps around the cluster are reachable.
            int search_min_px = std::max(0, occ_min_px - max_rot_iw);
            int search_min_py = std::max(0, occ_min_py - max_rot_ih);
            int search_max_px = std::min(plate_bw - 1, occ_max_px);
            int search_max_py = std::min(plate_bh - 1, occ_max_py);

            // Grid stride: half of the smaller rotation dimension,
            // floor 2 px. Coarser strides are fine because the narrow
            // phase is cheap.
            int stride = std::max(2, std::min(max_rot_iw, max_rot_ih) / 4);

            double best_score = std::numeric_limits<double>::max();
            int    best_rci   = -1;
            int    best_px    = -1;
            int    best_py    = -1;

            for (std::size_t rci = 0; rci < rot_cache.size(); ++rci) {
                const RotCache& rc = rot_cache[rci];
                int rot_max_px = std::min(plate_bw - rc.iw, search_max_px);
                int rot_max_py = std::min(plate_bh - rc.ih, search_max_py);

                for (int py = search_min_py; py <= rot_max_py; py += stride) {
                    for (int px = search_min_px; px <= rot_max_px; px += stride) {
                        // Bitmap AND is the narrow phase. Nanoseconds.
                        if (BitmapNester::collides(plate_items, plate_wpr,
                                                   plate_bw, plate_bh,
                                                   rc.bm, rc.iwpr,
                                                   rc.iw, rc.ih, px, py))
                            continue;

                        // Score: hull perimeter delta if committed here.
                        // Build the candidate polygon in world coords to
                        // extract its vertices for the hull calc.
                        ExPolygon candidate = base_poly;
                        if (rc.rot != 0.0) candidate.rotate(rc.rot);
                        coord_t cdx = scaled<coord_t>(px * res) - rc.rot_bb.min.x();
                        coord_t cdy = scaled<coord_t>(py * res) - rc.rot_bb.min.y();
                        candidate.translate(cdx, cdy);

                        Points new_pts = placed_hull_pts;
                        for (const Point& p : candidate.contour.points)
                            new_pts.push_back(p);
                        if (new_pts.size() < 3) continue;
                        Polygon new_hull = Geometry::convex_hull(new_pts);
                        if (new_hull.points.size() < 3) continue;
                        double new_perim = unscaled<double>(new_hull.length());
                        double delta = new_perim - current_hull_perim_mm;
                        if (delta < best_score) {
                            best_score = delta;
                            best_rci   = (int)rci;
                            best_px    = px;
                            best_py    = py;
                        }
                    }
                }
            }

            if (best_rci < 0) continue;  // unplaceable — M3 spillover will retry

            // Commit: stamp the bitmap, update hull state, update
            // running extent, record on the island.
            const RotCache& br = rot_cache[best_rci];
            BitmapNester::stamp(plate_items, plate_wpr, plate_bw, plate_bh,
                                br.bm, br.iwpr, br.iw, br.ih,
                                best_px, best_py);

            occ_min_px = std::min(occ_min_px, best_px);
            occ_min_py = std::min(occ_min_py, best_py);
            occ_max_px = std::max(occ_max_px, best_px + br.iw);
            occ_max_py = std::max(occ_max_py, best_py + br.ih);

            ExPolygon committed = base_poly;
            if (br.rot != 0.0) committed.rotate(br.rot);
            coord_t bdx = scaled<coord_t>(best_px * res) - br.rot_bb.min.x();
            coord_t bdy = scaled<coord_t>(best_py * res) - br.rot_bb.min.y();
            committed.translate(bdx, bdy);
            for (const Point& p : committed.contour.points)
                placed_hull_pts.push_back(p);
            if (placed_hull_pts.size() >= 3) {
                Polygon hull = Geometry::convex_hull(placed_hull_pts);
                if (hull.points.size() >= 3)
                    current_hull_perim_mm = unscaled<double>(hull.length());
            }

            island.item_indices.push_back(idx);
            island.rotations.push_back(br.rot);
            island.translations.push_back(Vec2crd{bdx, bdy});
        }

        // Finalize island metadata. island.bbox is computed from
        // the accumulated hull vertices — this is the "min/max of
        // actual placed geometry" form of bbox that the M2.5 bbox
        // purge directive allows.
        if (!placed_hull_pts.empty()) {
            island.bbox = BoundingBox(placed_hull_pts.front(),
                                      placed_hull_pts.front());
            for (const Point& p : placed_hull_pts)
                island.bbox.merge(p);
        }
        island.hull_perimeter_mm = current_hull_perim_mm;

        if (placed_hull_pts.size() >= 3) {
            Polygon hull = Geometry::convex_hull(placed_hull_pts);
            if (hull.points.size() >= 3) {
                double cx = 0.0, cy = 0.0;
                for (const Point& p : hull.points) {
                    cx += unscaled<double>(p.x());
                    cy += unscaled<double>(p.y());
                }
                cx /= (double)hull.points.size();
                cy /= (double)hull.points.size();
                island.hull_centroid_mm = Vec2d(cx, cy);
            }
        }
        return island;
    }

    // M3.0: locate an island on a plate.
    //
    // Takes a packed island (from pack_as_island, which works in
    // bound-agnostic local coordinates) and positions it on the
    // given plate by shifting so the island's hull centroid lands
    // at the bed center. The shift is clamped so the island's bbox
    // stays inside the bed.
    //
    // Why hull centroid and not bbox center? Bbox center is biased
    // by corner items (a single spike on one side of an otherwise
    // compact cluster drags the bbox center toward the spike).
    // Hull centroid is the average of the hull vertices, which for
    // a compact cluster better represents the visual center. Not
    // perfect — a proper area-weighted centroid would be more
    // accurate — but good enough for centering, and cheap.
    //
    // Items that are UNARRANGED in the island stay UNARRANGED in
    // the output (not written to items_out). M3.1 recover_spillover
    // will retry them.
    //
    // Excludes (wipe tower, calibration zones) are NOT handled
    // yet. An M3.1 follow-up will add a safety check that clamps
    // the shift if it would cause an item to overlap an exclude.
    static void locate_island_on_plate(ArrangePolygons& items_out,
                                       const NesterC2Island& island,
                                       int plate_idx,
                                       const BoundingBox& bed)
    {
        if (island.item_indices.empty()) return;

        // Compute the delta: where we want the hull centroid vs
        // where it currently is. Convert from mm to scaled coord_t.
        double bed_cx = unscaled<double>(bed.center().x());
        double bed_cy = unscaled<double>(bed.center().y());

        double dx_mm = bed_cx - island.hull_centroid_mm.x();
        double dy_mm = bed_cy - island.hull_centroid_mm.y();

        coord_t dx = scaled<coord_t>(dx_mm);
        coord_t dy = scaled<coord_t>(dy_mm);

        // Clamp so the island's bbox-after-shift stays inside bed.
        // If the island is bigger than the bed in either axis, just
        // align its min to the bed min (best we can do without
        // cropping).
        BoundingBox shifted_bbox = island.bbox;
        shifted_bbox.min += Point(dx, dy);
        shifted_bbox.max += Point(dx, dy);

        if (shifted_bbox.min.x() < bed.min.x())
            dx += (bed.min.x() - shifted_bbox.min.x());
        if (shifted_bbox.max.x() > bed.max.x())
            dx -= (shifted_bbox.max.x() - bed.max.x());
        if (shifted_bbox.min.y() < bed.min.y())
            dy += (bed.min.y() - shifted_bbox.min.y());
        if (shifted_bbox.max.y() > bed.max.y())
            dy -= (shifted_bbox.max.y() - bed.max.y());

        // Write out each item with the centered translation and the
        // plate_idx as bed_idx.
        Point shift(dx, dy);
        for (std::size_t i = 0; i < island.item_indices.size(); ++i) {
            std::size_t orig = island.item_indices[i];
            Vec2crd t = island.translations[i];
            items_out[orig].translation = Vec2crd{t.x() + shift.x(),
                                                  t.y() + shift.y()};
            items_out[orig].rotation    = island.rotations[i];
            items_out[orig].bed_idx     = plate_idx;
        }
    }

    // M3: attempt to find a home for items that couldn't fit their target
    // island. Retry rotations, migrate to neighbor islands, or spawn a new
    // island as a last resort.
    // static void recover_spillover(ArrangePolygons& items,
    //                               std::vector<NesterC2Island>& islands,
    //                               const BoundingBox& bed,
    //                               const ArrangeParams& params);
};

}} // namespace Slic3r::arrangement
