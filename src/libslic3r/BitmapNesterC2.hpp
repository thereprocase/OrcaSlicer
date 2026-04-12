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
        // C2 M2 pipeline:
        //   phase 0 — build per-item cache
        //   phase 1 — estimate K_min
        //   phase 2 — partition items into K groups
        //   phase 3 — pack each group into an Island (M2.2:
        //             delegation-based stub that still calls
        //             BitmapNester inside; M2.3 replaces with real
        //             hull-scored inner loop)
        //   phase 4 — locate each island on its plate (M2.2 stub:
        //             direct write-out using the island's absolute
        //             coordinates; M3 replaces with hull-centroid
        //             centering and spillover recovery)
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

        // Multi-group path: pack each group as an island and locate.
        for (std::size_t g = 0; g < groups.size(); ++g) {
            NesterC2Island island = pack_as_island(
                items, cache, groups[g], excludes, bed, params);
            locate_island_on_plate(items, island, (int)g);
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

        // Placed polygons in local (island) coordinates.
        ExPolygons placed_polys;
        Points     placed_hull_pts;  // accumulated vertices for hull calcs
        double     current_hull_perim_mm = 0.0;

        for (std::size_t idx : order) {
            const ArrangePolygon& ap = items_in[idx];
            ExPolygon base_poly = ap.poly;
            const std::vector<double>& rots_raw = ap.allowed_rotations;
            std::vector<double> rots = rots_raw;
            if (rots.empty()) rots.push_back(0.0);

            // First item: seed at origin, rotation 0.
            if (placed_polys.empty()) {
                ExPolygon seeded = base_poly;
                // Translate so bbox min lands at (0, 0) — predictable
                // seed position regardless of input polygon offsets.
                BoundingBox sbb = get_extents(seeded);
                seeded.translate(-sbb.min.x(), -sbb.min.y());
                placed_polys.push_back(seeded);

                for (const Point& p : seeded.contour.points)
                    placed_hull_pts.push_back(p);
                if (placed_hull_pts.size() >= 3) {
                    Polygon hull = Geometry::convex_hull(placed_hull_pts);
                    if (hull.points.size() >= 3)
                        current_hull_perim_mm = unscaled<double>(hull.length());
                }
                island.item_indices.push_back(idx);
                island.rotations.push_back(0.0);
                island.translations.push_back(Vec2crd{-sbb.min.x(), -sbb.min.y()});
                continue;
            }

            // Compute current cluster bbox for grid search window.
            BoundingBox cluster_bb = get_extents(placed_polys);

            // Grid stride: max_dim / 4 mm, floor 2 mm. Coarse enough to
            // keep candidate counts manageable.
            double max_dim_mm = 20.0;
            for (const auto& info : cache) {
                if (info.original_idx == idx) {
                    max_dim_mm = info.max_dim_mm;
                    break;
                }
            }
            double stride_mm = std::max(2.0, max_dim_mm / 4.0);
            coord_t stride = scaled<coord_t>(stride_mm);

            // Search window: cluster bbox expanded by max_dim on all
            // sides. Items can wrap around the cluster.
            coord_t margin = scaled<coord_t>(max_dim_mm * 1.25);
            BoundingBox search_bb = cluster_bb;
            search_bb.min -= Point(margin, margin);
            search_bb.max += Point(margin, margin);

            // Track the globally best (rot, px, py) across all rotations.
            double best_score  = std::numeric_limits<double>::max();
            double best_rot    = 0.0;
            Vec2crd best_trans = Vec2crd{0, 0};
            ExPolygon best_committed;
            bool best_found = false;

            for (double rot : rots) {
                ExPolygon rotated = base_poly;
                if (rot != 0.0) rotated.rotate(rot);
                BoundingBox rot_bb = get_extents(rotated);
                coord_t rot_w = rot_bb.size().x();
                coord_t rot_h = rot_bb.size().y();
                (void)rot_w; (void)rot_h;  // reserved for early bbox-vs-gap pruning

                for (coord_t y = search_bb.min.y();
                     y <= search_bb.max.y(); y += stride) {
                    for (coord_t x = search_bb.min.x();
                         x <= search_bb.max.x(); x += stride) {
                        // Candidate: rotated shape translated so its
                        // bbox min is at (x, y).
                        ExPolygon candidate = rotated;
                        candidate.translate(x - rot_bb.min.x(),
                                            y - rot_bb.min.y());

                        // Collision check vs all placed.
                        ExPolygons cand_vec{candidate};
                        bool collides = false;
                        for (const ExPolygon& placed : placed_polys) {
                            ExPolygons inter = intersection_ex(
                                cand_vec, ExPolygons{placed});
                            if (!inter.empty()) { collides = true; break; }
                        }
                        if (collides) continue;

                        // Score: new hull perimeter if this commits.
                        Points new_pts = placed_hull_pts;
                        for (const Point& p : candidate.contour.points)
                            new_pts.push_back(p);
                        if (new_pts.size() < 3) continue;
                        Polygon new_hull = Geometry::convex_hull(new_pts);
                        if (new_hull.points.size() < 3) continue;
                        double new_perim = unscaled<double>(new_hull.length());
                        double delta = new_perim - current_hull_perim_mm;
                        if (delta < best_score) {
                            best_score     = delta;
                            best_rot       = rot;
                            best_trans     = Vec2crd{
                                x - rot_bb.min.x(), y - rot_bb.min.y()};
                            best_committed = candidate;
                            best_found     = true;
                        }
                    }
                }
            }

            if (!best_found) {
                // Couldn't place this item in the search window. Mark
                // it unplaced; M3 recover_spillover will try harder.
                // (In M2.3's current state this just means the item's
                // bed_idx stays UNARRANGED in items_out.)
                continue;
            }

            // Commit: update placed_polys + hull state, record on island.
            placed_polys.push_back(best_committed);
            for (const Point& p : best_committed.contour.points)
                placed_hull_pts.push_back(p);
            if (placed_hull_pts.size() >= 3) {
                Polygon hull = Geometry::convex_hull(placed_hull_pts);
                if (hull.points.size() >= 3)
                    current_hull_perim_mm = unscaled<double>(hull.length());
            }
            island.item_indices.push_back(idx);
            island.rotations.push_back(best_rot);
            island.translations.push_back(best_trans);
        }

        // Finalize island metadata.
        if (!placed_polys.empty()) {
            island.bbox = get_extents(placed_polys);
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

    // M2.2 stub: locate an island on a plate.
    //
    // CURRENT STATE: writes the island's items directly back to the
    // caller's items[] vector with their bed_idx set to plate_idx.
    // No coordinate shift — the M2.2 delegation path uses world
    // coordinates (not relative-to-seed), so items are already where
    // BitmapNester::arrange placed them.
    //
    // M3 replaces this with a real two-step: shift so hull centroid
    // lands at bed center (clamped to bed bounds), then write out.
    // M3 also handles excludes (wipe tower collision avoidance).
    static void locate_island_on_plate(ArrangePolygons& items_out,
                                       const NesterC2Island& island,
                                       int plate_idx)
    {
        for (std::size_t i = 0; i < island.item_indices.size(); ++i) {
            std::size_t orig = island.item_indices[i];
            items_out[orig].translation = island.translations[i];
            items_out[orig].rotation    = island.rotations[i];
            // plate_idx here is the group index; for single-plate
            // fits we write 0, for multi-plate cases the groups
            // get distinct plate numbers.
            //
            // Items that couldn't place inside pack_as_island keep
            // their UNARRANGED state. They're candidates for M3's
            // spillover recovery.
            items_out[orig].bed_idx = (plate_idx);
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
