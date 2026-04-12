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

// Portable popcount for uint64_t — avoids __builtin_popcountll (GCC)
// and __popcnt64 (MSVC) portability issues.
inline int popcount64(uint64_t x) {
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (int)((x * 0x0101010101010101ULL) >> 56);
}

// ---------------------------------------------------------------------------
// Per-item information cached once per arrange call. Reused across every
// phase of the pipeline so that per-item work (rasterization, hull, bbox)
// only happens once. Grows as the phases come online — M0 holds the
// minimum necessary for the pass-through.
// ---------------------------------------------------------------------------

struct NesterC2ItemInfo {
    std::size_t original_idx = 0;

    // Geometric summary. Populated in build_cache (M1+).
    // silhouette_area_mm2: actual polygon area (pixel-count proxy
    //   at the 0.5 mm bitmap resolution) — use this instead of bbox
    //   area for sort/partition per the M2.5 bbox purge directive.
    // hull_area_mm2, hull_perimeter_mm: convex hull metrics.
    // height_mm: Z extent from ArrangePolygon::height (populated by
    //   Orca's ModelArrange.cpp from the 3D mesh).
    //
    // NOTE: max_dim_mm WAS here in M1 but was removed in M2.5.4 as
    // part of the bbox purge — it was a bbox-derived value and no
    // longer has any consumer in the C2 pipeline.
    double silhouette_area_mm2     = 0.0;
    double silhouette_perimeter_mm = 0.0;
    double hull_area_mm2           = 0.0;
    double hull_perimeter_mm       = 0.0;
    double height_mm               = 0.0;

    // Derived heuristics. Populated in build_cache.
    // hardness_score: hull_area / silhouette_area; higher = more
    //   concavity, harder to pack around neighbors, should go first.
    // priority_score: composite sort key. See build_cache for the
    //   weighting formula — silhouette area dominates, hardness and
    //   height are secondary. Larger parts placed first gives the
    //   greedy more room to find interlocks for the small parts.
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

// ---------------------------------------------------------------------------
// Per-item compaction state. Holds the winning rotation's bitmap (moved
// from RotCache at greedy commit time — zero alloc, zero copy). Used by
// compact_on_plate (Phase 4.5) to XOR-unstamp/restamp parts during the
// trash compactor iteration loop.
// ---------------------------------------------------------------------------

struct CompactItem {
    std::size_t           original_idx = 0;
    std::vector<uint64_t> bm;          // item bitmap (moved from RotCache)
    int                   iw = 0;
    int                   ih = 0;
    int                   iwpr = 0;
    int                   px = 0;      // current pixel position on plate
    int                   py = 0;
    double                rot = 0.0;   // committed rotation
};

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

    // Per-item compaction data. Populated during pack_as_island by
    // std::move from the winning RotCache entry at commit time. Used
    // by compact_on_plate (Phase 4.5) for XOR unstamp/restamp.
    std::vector<CompactItem> compact_items;

    // The plate bitmap resolution and dimensions used during packing.
    // Carried forward so compact_on_plate can rebuild the composite.
    double bitmap_res = 0.5;
    int    plate_bw   = 2048;
    int    plate_bh   = 2048;
    int    plate_wpr  = 32;
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

        // Full C2 pipeline for all group counts including k=1.
        // The k=1 C1 fallback was removed — C2 now handles all
        // cases through its own pack → locate → compact → spill
        // pipeline. This is required for the trash compactor and
        // height-desc sort to actually run on single-plate jobs
        // (which is most real-world usage).
        std::vector<NesterC2Island> islands;
        islands.reserve(groups.size());
        for (std::size_t g = 0; g < groups.size(); ++g) {
            NesterC2Island island = pack_as_island(
                items, cache, groups[g], excludes, bed, params);
            locate_island_on_plate(items, island, (int)g, bed);
            islands.push_back(std::move(island));
        }

        // M4.5 trash compactor: push parts inward from bed walls.
        // Runs after locate, before spillover. Primary goal: get
        // parts ON the plate.
        for (std::size_t g = 0; g < islands.size(); ++g) {
            compact_on_plate(items, islands[g], (int)g, bed);
        }

        // M3.1 spillover recovery: items that pack_as_island couldn't
        // fit (best_rci stayed -1) still have bed_idx == UNARRANGED.
        // Each gets its own plate centered on the bed. This is the
        // minimum viable spillover — guarantees every item ends up
        // placed somewhere, no items silently vanish. Fitting
        // spillovers onto existing plates with remaining space is an
        // optimization for a future milestone.
        recover_spillover(items, cache, (int)groups.size(), bed);
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

            info.height_mm = ap.height;

            // Priority score: composite sort key for placement order.
            // Silhouette area (pixel-count proxy) dominates — bigger
            // parts placed first gives the greedy more freedom to fit
            // small parts around them. Hardness and height are
            // secondary tiebreakers.
            //
            // M2.5.4 note: this replaces the old formula that weighted
            // hardness most. Silhouette area is a better primary key
            // per the boss's directive — a crescent and a square with
            // the same bbox have very different pixel counts, and the
            // crescent should go first because it's harder to fit.
            // Using silhouette area (which IS the pixel count up to
            // rasterization error) as the primary sort key achieves
            // this without adding a separate bitmap_pixel_count field.
            //
            // Weighting:
            //   ap.priority × 1,000,000  — caller priority dominates
            //   silhouette_area × 100    — pixel count is primary
            //   hardness × 10            — concavity ratio is secondary
            //   height × 1               — Z extent is tertiary
            //
            // For typical sizes (100-10000 mm²), silhouette area × 100
            // lands in 10k-1M range — much larger than hardness (1-2)
            // × 10 and typical heights (10-200). Area dominates.
            info.priority_score =
                ap.priority * 1000000.0
                + info.silhouette_area_mm2 * 100.0
                + info.hardness_score * 10.0
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

        // Placement order: height_mm desc primary, priority_score desc
        // tiebreak. M2.6 Council of Elrond / Sauron arbitration: tallest
        // items place first so the greedy hull-perimeter scorer naturally
        // clusters them at the local island origin; shorter items accrete
        // around the growing perimeter. Produces a passive height gradient
        // ("tall in middle, short at edges") for free — no new phase,
        // no new cache fields, no post-process compaction.
        //
        // Stability: std::stable_sort + continuous height key + secondary
        // priority_score tiebreak means a 0.01mm perturbation on a single
        // item moves at most that item in the ordering; equal-height
        // neighbors do not swap. Layout flinch is bounded by the input
        // perturbation.
        //
        // Invariants preserved:
        //   - Plate count minimum: sort changes trial order only, not
        //     the set of items accepted.
        //   - Bound-agnostic: sort key is intrinsic to each item.
        //   - No new allocations: sort is over an index vector.
        //   - Complexity O(N log N), dominated by the placement loop.
        //
        // Perf fix (also M2.6): previous comparator used std::find_if
        // over `cache` per comparison — O(N·M log N) for the full sort.
        // Replaced with a direct lookup table keyed on original_idx.
        std::vector<std::size_t> order = group.items;
        std::vector<std::size_t> cache_idx_by_orig(items_in.size(),
                                                   (std::size_t)-1);
        for (std::size_t ci = 0; ci < cache.size(); ++ci)
            cache_idx_by_orig[cache[ci].original_idx] = ci;
        std::stable_sort(
            order.begin(), order.end(),
            [&cache, &cache_idx_by_orig](std::size_t a, std::size_t b) {
                const auto& ia = cache[cache_idx_by_orig[a]];
                const auto& ib = cache[cache_idx_by_orig[b]];
                if (ia.height_mm != ib.height_mm)
                    return ia.height_mm > ib.height_mm;
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
        // Bitmap resolution = spacing / 2 (Nyquist: need ≥ 2 samples
        // across the minimum gap to distinguish gap from contact).
        // Clamped to [0.1, 0.5] mm/px.
        double spacing_mm = unscaled<double>(params.min_obj_distance);
        if (spacing_mm <= 0.0) spacing_mm = 1.0;  // sane default
        const double res = std::max(0.1, std::min(0.5, spacing_mm / 2.0));
        // Virtual plate = 2× bed at the Nyquist resolution.
        // 2× bed so the island can grow freely before locate clips.
        // Width rounded up to 64-bit word boundary.
        //
        // plate_px = 2 * bed_mm / res = 2 * bed_mm * 2 / spacing
        //          = 4 * bed_mm / spacing
        //
        // At 256mm bed, 0.5mm spacing: 4 * 256 / 0.5 = 2048 px.
        // The old hardcoded 2048 was this formula all along.
        double bed_w_mm = unscaled<double>(bed.size().x());
        double bed_h_mm = unscaled<double>(bed.size().y());
        const int plate_bw  = ((int)(bed_w_mm * 2.0 / res) + 63) & ~63;
        const int plate_bh  = (int)(bed_h_mm * 2.0 / res) + 1;
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
            // If the item is larger than the virtual plate, skip it —
            // it can't fit in imagination space, so it definitely can't
            // fit on the bed. Spillover will catch it.
            if (occ_max_px == 0) {
                const RotCache& rc = rot_cache[0];
                int px = seed_cx - rc.iw / 2;
                int py = seed_cy - rc.ih / 2;
                if (px < 0 || py < 0 ||
                    px + rc.iw > plate_bw || py + rc.ih > plate_bh)
                    continue;
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

                // M4.5: save winning bitmap for trash compactor.
                island.compact_items.push_back({idx,
                    std::move(rot_cache[0].bm),
                    rc.iw, rc.ih, rc.iwpr, px, py, rc.rot});
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

            // Refine pass: 1px stride in a window around the coarse
            // winner. The coarse stride can miss the exact interlock
            // position by up to stride pixels; the refine pass finds
            // the pixel-perfect optimum. Window: ±stride in each axis.
            {
                int ref_min_px = std::max(0, best_px - stride);
                int ref_min_py = std::max(0, best_py - stride);
                int ref_max_px = std::min(plate_bw - 1, best_px + stride);
                int ref_max_py = std::min(plate_bh - 1, best_py + stride);
                const RotCache& rc = rot_cache[best_rci];
                ref_max_px = std::min(ref_max_px, plate_bw - rc.iw);
                ref_max_py = std::min(ref_max_py, plate_bh - rc.ih);
                for (int ry = ref_min_py; ry <= ref_max_py; ++ry) {
                    for (int rx = ref_min_px; rx <= ref_max_px; ++rx) {
                        if (BitmapNester::collides(plate_items, plate_wpr,
                                                    plate_bw, plate_bh,
                                                    rc.bm, rc.iwpr,
                                                    rc.iw, rc.ih, rx, ry))
                            continue;
                        ExPolygon candidate = base_poly;
                        if (rc.rot != 0.0) candidate.rotate(rc.rot);
                        coord_t cdx = scaled<coord_t>(rx * res) - rc.rot_bb.min.x();
                        coord_t cdy = scaled<coord_t>(ry * res) - rc.rot_bb.min.y();
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
                            best_px    = rx;
                            best_py    = ry;
                        }
                    }
                }
            }

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

            // M4.5: save winning bitmap for trash compactor.
            island.compact_items.push_back({idx,
                std::move(rot_cache[best_rci].bm),
                br.iw, br.ih, br.iwpr, best_px, best_py, br.rot});
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
        island.bitmap_res = res;
        island.plate_bw   = plate_bw;
        island.plate_bh   = plate_bh;
        island.plate_wpr  = plate_wpr;
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

    // M3.1: spillover recovery. Items that pack_as_island couldn't
    // fit (no clear position in the scan window) stay UNARRANGED.
    // This gives each one its own plate, centered on the bed.
    //
    // Future optimization: try fitting spillovers onto existing
    // plates with remaining space before spawning new plates.
    static void recover_spillover(ArrangePolygons& items,
                                   const std::vector<NesterC2ItemInfo>& cache,
                                   int next_plate_idx,
                                   const BoundingBox& bed)
    {
        coord_t cx = bed.center().x();
        coord_t cy = bed.center().y();

        for (std::size_t i = 0; i < items.size(); ++i) {
            if (items[i].bed_idx != UNARRANGED) continue;

            // Guard against degenerate polygons (zero-area items that
            // rasterize to nothing). get_extents on an empty contour
            // produces an undefined bbox; calling center() on that is
            // UB. Just assign a plate without centering.
            // (Aragorn audit 2026-04-12)
            if (items[i].poly.contour.points.empty()) {
                items[i].translation = Vec2crd{0, 0};
                items[i].rotation    = 0.0;
                items[i].bed_idx     = next_plate_idx++;
                continue;
            }

            // Center the item on a new plate. Translation moves the
            // item's poly center to the bed center.
            BoundingBox item_bb = get_extents(items[i].poly);
            coord_t dx = cx - item_bb.center().x();
            coord_t dy = cy - item_bb.center().y();

            items[i].translation = Vec2crd{dx, dy};
            items[i].rotation    = 0.0;
            items[i].bed_idx     = next_plate_idx++;
        }
    }
    // ─── M4.5: Trash compactor primitives ────────────────────────────

    // XOR-remove an item bitmap from the composite. Inverse of stamp.
    // Precondition: the item's bits are set in the composite (i.e.,
    // it was stamped there). XOR clears exactly those bits. Safe
    // because the greedy placement loop checks collides() before
    // stamp(), guaranteeing zero overlap between items.
    static void xor_remove(std::vector<uint64_t>& plate,
                           int wpr, int bw, int bh,
                           const std::vector<uint64_t>& item_bm,
                           int iwpr, int iw, int ih,
                           int px, int py)
    {
        for (int iy = 0; iy < ih; ++iy) {
            int by = py + iy;
            if (by < 0 || by >= bh) continue;
            for (int wx = 0; wx < iwpr; ++wx) {
                uint64_t word = item_bm[(std::size_t)iy * iwpr + wx];
                if (word == 0) continue;
                int bed_bit = px + wx * 64;
                int bed_word = bed_bit / 64;
                int shift = bed_bit % 64;
                if (bed_word >= 0 && bed_word < wpr)
                    plate[(std::size_t)by * wpr + bed_word] ^= (word << shift);
                if (shift > 0 && bed_word + 1 >= 0 && bed_word + 1 < wpr)
                    plate[(std::size_t)by * wpr + bed_word + 1] ^= (word >> (64 - shift));
            }
        }
    }

    // M4.5: Trash compactor — bitmap compaction on the BED.
    //
    // Runs after locate_island_on_plate, before recover_spillover.
    // Primary goal: get parts ON the plate by pushing them inward
    // from the bed walls. Secondary: fill corners, rectangular shape.
    //
    // Algorithm:
    //   1. Rebuild composite bitmap from per-part CompactItems
    //      (adjusted for the locate shift)
    //   2. For each iteration (50 max, early-exit on convergence):
    //      a. For each part (height-desc, deterministic):
    //         - XOR-remove from composite
    //         - Compute wall pressure: inward from nearest bed wall
    //         - Try 1px step in wall direction
    //         - On collision: probe void direction (4 rays toward
    //           wall, 64px cap), try tangent step
    //         - XOR-stamp back at best position
    //   3. Write final positions back to items[].translation
    //
    // Quality gate: pixel overflow count. If overflow_after >
    // overflow_before, reject (should never happen with inward-only).
    static void compact_on_plate(ArrangePolygons& items,
                                  NesterC2Island& island,
                                  int plate_idx,
                                  const BoundingBox& bed)
    {
        if (island.compact_items.size() < 2) return;

        const double res = island.bitmap_res;
        const int bw  = island.plate_bw;
        const int bh  = island.plate_bh;
        const int wpr = island.plate_wpr;

        // Bed boundaries in pixel space.
        int bed_min_px = (int)(unscaled<double>(bed.min.x()) / res);
        int bed_min_py = (int)(unscaled<double>(bed.min.y()) / res);
        int bed_max_px = (int)(unscaled<double>(bed.max.x()) / res);
        int bed_max_py = (int)(unscaled<double>(bed.max.y()) / res);

        // Adjust CompactItem pixel positions for the locate shift.
        // locate_island_on_plate added a shift to each item's
        // translation. We need to convert current mm-space translations
        // back to pixel positions for bitmap operations.
        for (auto& ci : island.compact_items) {
            std::size_t orig = ci.original_idx;
            // The item's current translation in scaled coords → mm → px
            double tx_mm = unscaled<double>(items[orig].translation.x());
            double ty_mm = unscaled<double>(items[orig].translation.y());
            // Find the rotation's bbox offset (the same math used at
            // commit time to derive translation from pixel position).
            ExPolygon rotated = items[orig].poly;
            if (ci.rot != 0.0) rotated.rotate(ci.rot);
            BoundingBox rot_bb = get_extents(rotated);
            // px = (translation_mm + rot_bb.min_mm) / res
            // Use lround, not (int) cast — truncation toward zero
            // can drift by 1px due to IEEE 754 roundtrip through
            // scaled/unscaled, leaving ghost bits in the composite.
            // (Sauron audit 2026-04-12)
            ci.px = (int)std::lround((tx_mm + unscaled<double>(rot_bb.min.x())) / res);
            ci.py = (int)std::lround((ty_mm + unscaled<double>(rot_bb.min.y())) / res);
        }

        // Build fresh composite from per-part bitmaps.
        std::vector<uint64_t> composite((std::size_t)wpr * bh, 0);
        for (const auto& ci : island.compact_items) {
            if (ci.bm.empty()) continue;
            BitmapNester::stamp(composite, wpr, bw, bh,
                                ci.bm, ci.iwpr, ci.iw, ci.ih,
                                ci.px, ci.py);
        }

        // Count initial overflow (set bits outside bed).
        auto count_overflow = [&](const std::vector<uint64_t>& plate) {
            int overflow = 0;
            for (int y = 0; y < bh; ++y) {
                if (y < bed_min_py || y >= bed_max_py) {
                    for (int w = 0; w < wpr; ++w)
                        overflow += popcount64(
                            plate[(std::size_t)y * wpr + w]);
                    continue;
                }
                for (int w = 0; w < wpr; ++w) {
                    int bit_start = w * 64;
                    int bit_end   = bit_start + 64;
                    uint64_t word = plate[(std::size_t)y * wpr + w];
                    if (word == 0) continue;
                    if (bit_start >= bed_min_px && bit_end <= bed_max_px)
                        continue;  // fully inside bed
                    // Partially outside — count out-of-bed bits
                    for (int b = 0; b < 64; ++b) {
                        if (!(word & (uint64_t(1) << b))) continue;
                        int px = bit_start + b;
                        if (px < bed_min_px || px >= bed_max_px)
                            ++overflow;
                    }
                }
            }
            return overflow;
        };

        int overflow_before = count_overflow(composite);

        // Save positions for quality gate rollback.
        struct SavedPos { int px, py; };
        std::vector<SavedPos> saved;
        saved.reserve(island.compact_items.size());
        for (const auto& ci : island.compact_items)
            saved.push_back({ci.px, ci.py});

        // ─── Main compaction loop ────────────────────────────────
        constexpr int MAX_ITER = 50;
        for (int iter = 0; iter < MAX_ITER; ++iter) {
            bool any_moved = false;

            for (auto& ci : island.compact_items) {
                if (ci.bm.empty()) continue;

                // 1. Wall pressure: inward from nearest bed wall.
                int cx = ci.px + ci.iw / 2;
                int cy = ci.py + ci.ih / 2;
                int dl = cx - bed_min_px;
                int dr = bed_max_px - cx;
                int dt = cy - bed_min_py;
                int db = bed_max_py - cy;
                int dmin = std::min({dl, dr, dt, db});

                int step_x = 0, step_y = 0;
                if (dmin <= 0 || dmin == dl)      step_x =  1;
                else if (dmin == dr)              step_x = -1;
                else if (dmin == dt)              step_y =  1;
                else if (dmin == db)              step_y = -1;

                // Guard: part at exact center → no wall pressure.
                if (step_x == 0 && step_y == 0) continue;

                // 2. XOR-remove from composite.
                xor_remove(composite, wpr, bw, bh,
                           ci.bm, ci.iwpr, ci.iw, ci.ih,
                           ci.px, ci.py);

                int new_px = ci.px + step_x;
                int new_py = ci.py + step_y;
                bool moved = false;

                // 3. Try wall-pressure step.
                if (new_px >= 0 && new_py >= 0 &&
                    new_px + ci.iw <= bw && new_py + ci.ih <= bh &&
                    !BitmapNester::collides(composite, wpr, bw, bh,
                                            ci.bm, ci.iwpr, ci.iw, ci.ih,
                                            new_px, new_py)) {
                    ci.px = new_px;
                    ci.py = new_py;
                    moved = true;
                } else {
                    // 4. Void-pull: probe 4 directions toward wall
                    //    for tangent sliding.
                    int best_void_px = ci.px;
                    int best_void_py = ci.py;

                    // Try perpendicular steps (tangent to wall).
                    int tangents[][2] = {{0, 1}, {0, -1},
                                         {1, 0}, {-1, 0}};
                    int best_void_len = 0;
                    for (auto& t : tangents) {
                        // Skip the direction we already tried.
                        if (t[0] == step_x && t[1] == step_y) continue;
                        int tp = ci.px + t[0];
                        int tq = ci.py + t[1];
                        if (tp < 0 || tq < 0 ||
                            tp + ci.iw > bw || tq + ci.ih > bh)
                            continue;
                        if (!BitmapNester::collides(
                                composite, wpr, bw, bh,
                                ci.bm, ci.iwpr, ci.iw, ci.ih,
                                tp, tq)) {
                            // Probe void depth in this direction.
                            int void_len = 0;
                            int probe_x = tp + t[0];
                            int probe_y = tq + t[1];
                            while (void_len < 64 &&
                                   probe_x >= 0 && probe_y >= 0 &&
                                   probe_x + ci.iw <= bw &&
                                   probe_y + ci.ih <= bh &&
                                   !BitmapNester::collides(
                                       composite, wpr, bw, bh,
                                       ci.bm, ci.iwpr, ci.iw, ci.ih,
                                       probe_x, probe_y)) {
                                ++void_len;
                                probe_x += t[0];
                                probe_y += t[1];
                            }
                            if (void_len > best_void_len) {
                                best_void_len = void_len;
                                best_void_px = tp;
                                best_void_py = tq;
                            }
                        }
                    }
                    if (best_void_px != ci.px || best_void_py != ci.py) {
                        ci.px = best_void_px;
                        ci.py = best_void_py;
                        moved = true;
                    }
                }

                // 5. Re-stamp at (possibly new) position.
                BitmapNester::stamp(composite, wpr, bw, bh,
                                    ci.bm, ci.iwpr, ci.iw, ci.ih,
                                    ci.px, ci.py);
                if (moved) any_moved = true;
            }

            if (!any_moved) break;
        }

        // ─── Quality gate: pixel overflow count ──────────────────
        int overflow_after = count_overflow(composite);
        if (overflow_after > overflow_before) {
            // Reject compaction — restore saved positions.
            for (std::size_t i = 0; i < island.compact_items.size(); ++i) {
                island.compact_items[i].px = saved[i].px;
                island.compact_items[i].py = saved[i].py;
            }
            return;
        }

        // ─── Grid-snap post-pass ─────────────────────────────────
        int grid_px = std::max(1, (int)(2.0 / res));
        // Rebuild composite for snap checks.
        std::fill(composite.begin(), composite.end(), uint64_t(0));
        for (const auto& ci : island.compact_items) {
            if (ci.bm.empty()) continue;
            BitmapNester::stamp(composite, wpr, bw, bh,
                                ci.bm, ci.iwpr, ci.iw, ci.ih,
                                ci.px, ci.py);
        }
        for (auto& ci : island.compact_items) {
            if (ci.bm.empty()) continue;
            int snap_x = (ci.px / grid_px) * grid_px;
            int snap_y = (ci.py / grid_px) * grid_px;
            if (snap_x == ci.px && snap_y == ci.py) continue;
            xor_remove(composite, wpr, bw, bh,
                       ci.bm, ci.iwpr, ci.iw, ci.ih,
                       ci.px, ci.py);
            if (snap_x >= 0 && snap_y >= 0 &&
                snap_x + ci.iw <= bw && snap_y + ci.ih <= bh &&
                !BitmapNester::collides(composite, wpr, bw, bh,
                                        ci.bm, ci.iwpr, ci.iw, ci.ih,
                                        snap_x, snap_y)) {
                ci.px = snap_x;
                ci.py = snap_y;
            }
            BitmapNester::stamp(composite, wpr, bw, bh,
                                ci.bm, ci.iwpr, ci.iw, ci.ih,
                                ci.px, ci.py);
        }

        // ─── Write back to items[].translation ───────────────────
        for (const auto& ci : island.compact_items) {
            std::size_t orig = ci.original_idx;
            ExPolygon rotated = items[orig].poly;
            if (ci.rot != 0.0) rotated.rotate(ci.rot);
            BoundingBox rot_bb = get_extents(rotated);
            coord_t tx = scaled<coord_t>(ci.px * res) - rot_bb.min.x();
            coord_t ty = scaled<coord_t>(ci.py * res) - rot_bb.min.y();
            items[orig].translation = Vec2crd{tx, ty};
        }
    }
};

}} // namespace Slic3r::arrangement
