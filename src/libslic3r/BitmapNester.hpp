#pragma once

// BitmapNester: 2D bitmap-based concave shape nesting for OrcaSlicer
//
// Drop-in replacement for the libnest2d arrange path when "Use actual part
// shape" is enabled. Handles every scenario the default arranger handles:
//
//   - Per-item inflation (brim, tree support, user spacing)
//   - Per-plate excludes (wipe tower, calibration zones, fixed items)
//   - Multi-plate overflow with exclude propagation
//   - Priority ordering (higher priority placed first)
//   - Per-item allowed rotations
//   - Item ID assignment for sequential print ordering
//   - Bed shrinkage (skirt, brim, sequential clearance)
//
// When this nester is off, the existing libnest2d convex hull path is used
// with zero code changes.

#include "ExPolygon.hpp"
#include "Arrange.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"   // intersection_ex used by post-centering safety check
#include "Geometry/ConvexHull.hpp"
#include <vector>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <numeric>
#include <cmath>
#include <cassert>
#include <random>     // smart-shuffle restart driver uses std::mt19937
#include <utility>

namespace Slic3r { namespace arrangement {

class BitmapNester {
public:
    static void arrange(ArrangePolygons &items,
                        const ArrangePolygons &excludes,
                        const BoundingBox &bed,
                        const ArrangeParams &params)
    {
        if (items.empty()) return;

        double res = 0.5; // mm per pixel (nominal, recomputed after cap)

        // Effective bed after shrinkage
        BoundingBox ebed = bed;
        coord_t shrink_x = scaled(params.bed_shrink_x);
        coord_t shrink_y = scaled(params.bed_shrink_y);
        ebed.min.x() += shrink_x;
        ebed.min.y() += shrink_y;
        ebed.max.x() -= shrink_x;
        ebed.max.y() -= shrink_y;

        if (ebed.min.x() >= ebed.max.x() || ebed.min.y() >= ebed.max.y())
            return;

        coord_t bed_w = ebed.max.x() - ebed.min.x();
        coord_t bed_h = ebed.max.y() - ebed.min.y();
        double bed_w_mm = unscaled<double>(bed_w);
        double bed_h_mm = unscaled<double>(bed_h);
        int bw = std::max(1, (int)std::ceil(bed_w_mm / res));
        int bh = std::max(1, (int)std::ceil(bed_h_mm / res));

        // Cap bitmap size to prevent runaway allocation (L1: recompute res after cap)
        if ((int64_t)bw * bh > 16'000'000) {
            double sc = std::sqrt(16'000'000.0 / ((double)bw * bh));
            bw = std::max(1, (int)(bw * sc));
            bh = std::max(1, (int)(bh * sc));
        }
        // Recompute effective resolution from actual pixel dimensions
        res = std::max(bed_w_mm / bw, bed_h_mm / bh);

        int wpr = (bw + 63) / 64;

        // Sort items by priority (descending), then longest-dimension
        // (descending), then area (descending) as a tiebreaker.
        //
        // FFD (First-Fit-Decreasing) heuristic: place the items whose
        // bounding box has the largest max dimension first. For mixed
        // real-world plates this puts tall/long parts down while the
        // bed is empty, anchoring the cluster outline. Small parts then
        // tuck around them. This is also what the user noticed matters
        // for print-head travel: the big parts dominate head motion, so
        // committing their positions first (where they fit tightest)
        // matters more than committing a small part's position first.
        //
        // Max dim is invariant under 90/180/270 rotation (which is all
        // we support), so the ordering is stable against allowed_rotations.
        auto max_dim_of = [](const ArrangePolygon &it) -> double {
            BoundingBox bb = get_extents(it.poly);
            return std::max(unscaled<double>(bb.size().x()),
                            unscaled<double>(bb.size().y()));
        };
        std::vector<size_t> order(items.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (items[a].priority != items[b].priority)
                return items[a].priority > items[b].priority;
            double da = max_dim_of(items[a]);
            double db = max_dim_of(items[b]);
            if (da != db) return da > db;
            return std::abs(items[a].poly.area()) > std::abs(items[b].poly.area());
        });

        // Snapshot each item's caller-supplied rotation BEFORE the placement
        // loop touches it. The greedy placement loop overwrites
        // `item.rotation` with `base_rot + best_rot` at commit time, which
        // means by the time consolidation runs, `item.rotation` is no
        // longer the caller's pre-rotation — it's that plus whatever the
        // greedy chose. The consolidation pass needs the original value to
        // compose new rotations against (otherwise it would discard the
        // caller's `align_to_y_axis` pre-rotation when migrating). Save
        // them up front in a parallel vector indexed by items[].
        std::vector<double> base_rotations(items.size(), 0.0);
        for (size_t i = 0; i < items.size(); ++i)
            base_rotations[i] = items[i].rotation;

        // Build rotation list for items that allow any rotation
        std::vector<double> default_rotations;
        if (params.allow_rotations) {
            for (int i = 0; i < 4; ++i)
                default_rotations.push_back(i * M_PI / 2.0);
        } else {
            default_rotations.push_back(0.0);
        }

        constexpr int MAX_PLATES = 36;

        // Per-plate forbidden bitmaps. Two parallel vectors, one plate entry per
        // index:
        //   plate_items[p]    — stamps of placed item footprints (CONCAVE). Items
        //                       test their concave bitmap against this so they
        //                       can legitimately interlock into each other's
        //                       concave cavities (the whole point of the branch).
        //   plate_excludes[p] — stamps of exclude regions (wipe tower, cali
        //                       zones, fixed unselected items). Items test their
        //                       CONVEX HULL bitmap against this so placements
        //                       agree with PartPlate::check_outside, which uses
        //                       instance->convex_hull_2d() for exclude
        //                       intersection. See war-council-report.md C15.
        std::vector<std::vector<uint64_t>> plate_items;
        std::vector<std::vector<uint64_t>> plate_excludes;
        // Per-plate running cluster bounding box in pixel coordinates.
        // Empty plates hold an inverted-range sentinel (minx > maxx) so the
        // first placement computes a zero-area starting cluster. Used by the
        // bottom-left-fill scoring function: a candidate placement's cost is
        // the increase in cluster-bbox area it would cause, with anchor
        // distance as a tiebreaker. This is what lets concave shapes
        // interlock — a rotated L slotted into another L's notch has zero
        // delta-area and beats any non-interlocking position.
        struct ClusterBB {
            int minx, miny, maxx, maxy; // inclusive, pixel space; empty if minx > maxx
            bool empty() const { return minx > maxx || miny > maxy; }
            int64_t area() const {
                if (empty()) return 0;
                return (int64_t)(maxx - minx + 1) * (int64_t)(maxy - miny + 1);
            }
        };
        std::vector<ClusterBB> cluster_bb;
        auto ensure_plate = [&](int idx) {
            while ((int)plate_items.size() <= idx) {
                plate_items.emplace_back((size_t)wpr * bh, uint64_t(0));
                plate_excludes.emplace_back((size_t)wpr * bh, uint64_t(0));
                cluster_bb.push_back(ClusterBB{std::numeric_limits<int>::max(), std::numeric_limits<int>::max(),
                                       std::numeric_limits<int>::min(), std::numeric_limits<int>::min()});
            }
        };
        ensure_plate(0);

        // Stamp an ExPolygon (in scaled bed-relative coords) onto a plate's
        // exclude bitmap. Only used during setup to register static obstacles
        // (wipe tower, unselected items on a specific plate, etc.). Item
        // stamps at placement time go into plate_items directly.
        auto stamp_exclude_poly = [&](int plate_idx, const ExPolygon &poly_bed_rel) {
            ensure_plate(plate_idx);
            int iw, ih, iwpr;
            auto bm = rasterize(poly_bed_rel, res, bw, bh, iw, ih, iwpr);
            if (bm.empty()) return;
            BoundingBox pbb = get_extents(poly_bed_rel);
            int px = std::max(0, (int)(unscaled<double>(pbb.min.x()) / res));
            int py = std::max(0, (int)(unscaled<double>(pbb.min.y()) / res));
            stamp(plate_excludes[plate_idx], wpr, bw, bh, bm, iwpr, iw, ih, px, py);
        };

        // Stamp excludes onto their respective plates (L7: skip UNARRANGED excludes)
        for (auto &ex : excludes) {
            if (ex.bed_idx < 0) continue;
            int plate = ex.bed_idx;
            ExPolygon epoly = ex.poly;
            if (ex.rotation != 0.0) epoly.rotate(ex.rotation);
            epoly.translate(ex.translation.x() - ebed.min.x(),
                           ex.translation.y() - ebed.min.y());
            if (ex.inflation > 0) {
                ExPolygons infl = offset_ex(epoly, ex.inflation);
                if (!infl.empty()) epoly = infl.front();
            }
            stamp_exclude_poly(plate, epoly);
        }

        // Pre-rasterize bed exclusion zones (calibration areas, etc.) so they
        // can be stamped onto every plate, including overflow plates created later.
        struct BedZone { std::vector<uint64_t> bm; int iw, ih, iwpr, px, py; };
        std::vector<BedZone> bed_zones;
        for (auto &ex : params.excluded_regions) {
            ExPolygon epoly = ex.poly;
            if (ex.rotation != 0.0) epoly.rotate(ex.rotation);
            epoly.translate(ex.translation.x() - ebed.min.x(),
                           ex.translation.y() - ebed.min.y());
            int ziw, zih, ziwpr;
            auto zbm = rasterize(epoly, res, bw, bh, ziw, zih, ziwpr);
            if (zbm.empty()) continue;
            BoundingBox zbb = get_extents(epoly);
            int zpx = std::max(0, (int)(unscaled<double>(zbb.min.x()) / res));
            int zpy = std::max(0, (int)(unscaled<double>(zbb.min.y()) / res));
            stamp(plate_excludes[0], wpr, bw, bh, zbm, ziwpr, ziw, zih, zpx, zpy);
            bed_zones.push_back({std::move(zbm), ziw, zih, ziwpr, zpx, zpy});
        }

        // Patch ensure_plate to stamp bed exclusion zones onto new plates
        auto ensure_plate_with_zones = [&](int idx) {
            int old_count = (int)plate_items.size();
            ensure_plate(idx);
            for (int p = old_count; p <= idx; ++p) {
                for (auto &z : bed_zones)
                    stamp(plate_excludes[p], wpr, bw, bh, z.bm, z.iwpr, z.iw, z.ih, z.px, z.py);
            }
        };

        int current_plate = 0;
        for (auto &ex : excludes) {
            if (ex.bed_idx > 0)
                current_plate = std::max(current_plate, ex.bed_idx);
        }

        // ── Smart-shuffle restart driver ─────────────────────────────
        // For low part counts the placement order has a large effect
        // on whether everything fits. Center-greedy in particular can
        // trap the second item when the first lands at the bed center
        // and fragments the surrounding space. Random-permutation
        // restarts is the standard fix (see docs/FFT_NESTER_IDEA.md
        // Section 3 for the parallelism design that scales this up).
        // Here we run them sequentially and only when N is small
        // enough that 4-8 retries are cheap. The first attempt always
        // uses the priority+area-desc baseline; subsequent attempts
        // shuffle within priority buckets using deterministic seeds.
        struct InputState { Vec2crd t; double r; int b; int id; };
        std::vector<InputState> input_state(items.size());
        for (size_t i = 0; i < items.size(); ++i)
            input_state[i] = {items[i].translation, base_rotations[i],
                              items[i].bed_idx, items[i].itemid};

        // Snapshot the post-bed-setup plate state so each restart can
        // rebuild from this point without re-running the setup phase.
        const int initial_plate_count    = (int)plate_items.size();
        const int initial_current_plate  = current_plate;
        const auto initial_plate_excludes = plate_excludes;

        // Best result tracking. Score = (placed_count descending,
        // max_bed_idx ascending, total_cluster_area ascending). The
        // placed-count primary is critical for stopcondition-aware
        // correctness: an attempt where the user's stopcondition
        // halted placement early must NEVER beat an earlier attempt
        // that placed more items, even if its max_bed_idx looks
        // numerically better.
        struct BestState {
            bool    valid      = false;
            int     placed     = -1;
            int     max_bed    = std::numeric_limits<int>::max();
            int64_t total_area = std::numeric_limits<int64_t>::max();
            std::vector<Vec2crd> trans;
            std::vector<double>  rots;
            std::vector<int>     beds;
            std::vector<int>     iids;
        };
        BestState best;

        auto reset_for_restart = [&]() {
            for (size_t i = 0; i < items.size(); ++i) {
                items[i].translation = input_state[i].t;
                items[i].rotation    = input_state[i].r;
                items[i].bed_idx     = input_state[i].b;
                items[i].itemid      = input_state[i].id;
            }
            plate_items.assign(initial_plate_count,
                               std::vector<uint64_t>((size_t)wpr * bh, 0));
            plate_excludes = initial_plate_excludes;
            cluster_bb.assign(initial_plate_count, ClusterBB{
                std::numeric_limits<int>::max(),
                std::numeric_limits<int>::max(),
                std::numeric_limits<int>::min(),
                std::numeric_limits<int>::min()});
            current_plate = initial_current_plate;
        };

        // Returns (placed_count, max_bed_idx, total_cluster_area).
        // Lex order with placed_count flipped to descending (more is
        // better) for the comparator below.
        struct AttemptScore { int placed; int max_bed; int64_t area; };
        auto score_current = [&]() -> AttemptScore {
            int placed = 0;
            int max_bed = -1;
            for (const auto &it : items) {
                if (it.bed_idx != UNARRANGED) {
                    ++placed;
                    if (it.bed_idx > max_bed) max_bed = it.bed_idx;
                }
            }
            int64_t area = 0;
            for (const auto &c : cluster_bb)
                if (!c.empty()) area += c.area();
            return {placed, max_bed, area};
        };

        auto save_best = [&](const AttemptScore &s) {
            best.valid      = true;
            best.placed     = s.placed;
            best.max_bed    = s.max_bed;
            best.total_area = s.area;
            best.trans.resize(items.size());
            best.rots.resize(items.size());
            best.beds.resize(items.size());
            best.iids.resize(items.size());
            for (size_t i = 0; i < items.size(); ++i) {
                best.trans[i] = items[i].translation;
                best.rots[i]  = items[i].rotation;
                best.beds[i]  = items[i].bed_idx;
                best.iids[i]  = items[i].itemid;
            }
        };

        // Number of attempts based on N. Restarts are sequential and
        // each runs the full placement+consolidation+centering pipeline,
        // so the runtime is num_attempts × baseline_time. We cap at 8
        // for very small N (still milliseconds total) and 4 for medium
        // N. For large N we run 3 attempts — the order space is too big
        // to sample exhaustively, but 3 shuffles is cheap insurance
        // against a pathological ordering and typically runs in a few
        // seconds. Short-circuits as soon as a single-plate solution
        // with all items placed is found (see early-exit below).
        const int num_attempts =
            (items.size() <= 6)  ? 8 :
            (items.size() <= 16) ? 4 : 3;

        // Snapshot the original `order` for perturbation seeds.
        const std::vector<size_t> base_order = order;

        // Lambda wrapping the placement+consolidation+post-centering
        // pipeline. Captures everything by reference. Each invocation
        // mutates items[], plate_items[], cluster_bb[], current_plate.
        // Reset between calls via reset_for_restart().
        auto do_one_pass = [&]() {

        int item_sequence = 0;

        // Place each item
        unsigned items_remaining = (unsigned)items.size();
        for (size_t oi = 0; oi < order.size(); ++oi) {
            size_t idx = order[oi];
            auto &item = items[idx];

            if (params.stopcondition && params.stopcondition()) break;

            if (params.progressind)
                params.progressind(items_remaining--, "(placing parts)");

            // Per-item inflation
            int pad_px = 0;
            if (item.inflation > 0)
                pad_px = std::max(1, (int)std::ceil(unscaled<double>(item.inflation) / res));
            else if (params.min_obj_distance > 0)
                pad_px = std::max(1, (int)std::ceil(unscaled<double>(params.min_obj_distance / 2) / res));

            const auto &rots = !item.allowed_rotations.empty()
                               ? item.allowed_rotations : default_rotations;

            // Capture any pre-arrange rotation the caller placed on item.rotation
            // (e.g. from update_selected_items_axis_align). The nester's own
            // allowed_rotations are ADDITIONAL rotations tried on top of this.
            // We apply base_rot to the shapes once before the inner loop so the
            // rot_cache entries each represent base_rot + rot, not rot alone.
            const double base_rot = item.rotation;

            // L14: pre-rasterize all rotations for this item so we don't
            // re-rasterize on every plate retry. For 75 items × 4 rotations × 6
            // plates, this cuts rasterize calls from ~1000 to ~300.
            //
            // bm          — concave silhouette bitmap, tested against plate_items.
            // convex_bm   — convex hull of the same silhouette, tested against
            //               plate_excludes so that placements agree with Orca's
            //               PartPlate::check_outside validation (which uses
            //               instance->convex_hull_2d() for exclude intersection).
            //               Both bitmaps share iw/ih/iwpr because the convex
            //               hull of a polygon shares its bounding box.
            struct RotCache {
                double rot;
                std::vector<uint64_t> bm;
                std::vector<uint64_t> convex_bm;
                BoundingBox inflated_bb;
                int iw, ih, iwpr, coarse;
            };
            std::vector<RotCache> rot_cache;
            rot_cache.reserve(rots.size());

            // Build the silhouette input for this item. If the caller populated
            // item.concave_regions (ArrangeJob does this when use_concave_shapes
            // is on), use every island — multi-volume parts, dumbbells, and
            // shapes whose top/bottom projections diverge have legitimately
            // disconnected footprints that must all be collision-checked.
            // Otherwise fall back to the single-region item.poly (libnest2d-
            // compatible convex hull or simple silhouette).
            ExPolygons base_shapes;
            if (!item.concave_regions.empty()) {
                base_shapes = item.concave_regions;
            } else {
                base_shapes.push_back(item.poly);
            }

            // Apply the caller's pre-rotation once so that every subsequent
            // per-rot iteration composes on top of it rather than ignoring it.
            if (base_rot != 0.0) {
                for (ExPolygon &s : base_shapes) s.rotate(base_rot);
            }

            for (double rot : rots) {
                ExPolygons rshapes = base_shapes;
                if (rot != 0.0) {
                    for (ExPolygon &s : rshapes) s.rotate(rot);
                }

                // Rasterize the un-inflated polygon(s), then dilate the bitmap.
                // This avoids Clipper's offset_ex entirely in the hot path.
                int raw_iw, raw_ih, raw_iwpr;
                auto raw_bm = rasterize(rshapes, res, bw, bh, raw_iw, raw_ih, raw_iwpr);
                if (raw_bm.empty()) continue;

                // Build the convex hull of the same rotated silhouette and
                // rasterize it at the same bbox so the two bitmaps share
                // iw/ih/iwpr. Geometry::convex_hull on ExPolygons returns a
                // Polygon; we wrap it as an ExPolygon (no holes) to feed the
                // rasterizer. The hull's bbox is identical to rshapes' bbox
                // (extreme points are shared), so dimensions match.
                ExPolygon hull_expoly;
                hull_expoly.contour = Geometry::convex_hull(rshapes);
                if (hull_expoly.contour.empty()) continue;
                int hull_raw_iw, hull_raw_ih, hull_raw_iwpr;
                auto hull_raw_bm = rasterize(hull_expoly, res, bw, bh,
                                             hull_raw_iw, hull_raw_ih, hull_raw_iwpr);
                if (hull_raw_bm.empty()) continue;

                // The hull share bbox with rshapes so raw dimensions should
                // agree. If they don't (one-pixel rounding at the edge), take
                // the MAX so the hull mask is guaranteed to cover the concave
                // one — being slightly more conservative for the exclude check
                // is safe; being slightly less conservative would be a bug.
                if (hull_raw_iw != raw_iw || hull_raw_ih != raw_ih) {
                    continue;  // mismatched bbox: skip this rotation rather
                               // than feed inconsistent bitmaps to the scan.
                }

                // Dilate both bitmaps by pad_px in all directions. Same
                // dilation geometry so the convex mask stays aligned with
                // the concave mask at the same (px, py) placement position.
                int iw = raw_iw, ih = raw_ih, iwpr_item = raw_iwpr;
                std::vector<uint64_t> ibm;
                std::vector<uint64_t> hull_ibm;
                BoundingBox inflated_bb = get_extents(rshapes);

                if (pad_px > 0) {
                    iw = raw_iw + 2 * pad_px;
                    ih = raw_ih + 2 * pad_px;
                    if (iw > bw || ih > bh) continue;
                    iwpr_item = (iw + 63) / 64;
                    ibm = dilate_bitmap(raw_bm, raw_iwpr, raw_iw, raw_ih,
                                        pad_px, iwpr_item, iw, ih);
                    hull_ibm = dilate_bitmap(hull_raw_bm, raw_iwpr, raw_iw, raw_ih,
                                             pad_px, iwpr_item, iw, ih);
                    // Adjust inflated bbox to account for dilation
                    coord_t pad_sc = scaled(pad_px * res);
                    inflated_bb.min -= Vec2crd(pad_sc, pad_sc);
                    inflated_bb.max += Vec2crd(pad_sc, pad_sc);
                } else {
                    ibm = std::move(raw_bm);
                    hull_ibm = std::move(hull_raw_bm);
                    if (iw > bw || ih > bh) continue;
                }

                int coarse = std::clamp(std::min(iw, ih) / 8, 8, 128);
                rot_cache.push_back({rot, std::move(ibm), std::move(hull_ibm),
                                     inflated_bb, iw, ih, iwpr_item, coarse});
            }

            // Dual-bitmap collision check: a candidate position is valid only
            // when BOTH the concave bitmap clears plate_items AND the convex
            // hull bitmap clears plate_excludes at the same (px, py).
            auto position_clear = [&](int plate_idx, const RotCache &rc,
                                      int px, int py) -> bool {
                if (collides(plate_items[plate_idx], wpr, bw, bh,
                             rc.bm, rc.iwpr, rc.iw, rc.ih, px, py))
                    return false;
                if (collides(plate_excludes[plate_idx], wpr, bw, bh,
                             rc.convex_bm, rc.iwpr, rc.iw, rc.ih, px, py))
                    return false;
                return true;
            };

            bool placed = false;
            int best_px = 0, best_py = 0;
            double best_rot = 0.0;
            int best_plate = -1;
            // Index into rot_cache rather than copying the bitmap — the cache
            // lives for the lifetime of the placement scan and this avoids a
            // full bitmap copy on every successful placement (which becomes a
            // hot-loop allocation once Sprint 1 Task 8 converts this to a
            // scored scan with best-fit selection).
            int best_rc_idx = -1;

            // Anchor resolution per plate. Default arranger clusters items
            // around the wipe tower (libnest2d::on_preload shifts the starting
            // point to itm.boundingBox().center()). We reproduce that effect
            // by anchoring the scoring function on the wipe tower center when
            // one exists on the plate — items score best when their bbox
            // center is near the anchor, and collisions keep them from
            // actually overlapping the tower.
            //
            // When no wipe tower is present AND align_center is the default
            // (0.5, 0.5), the anchor falls back to (0, 0) in pixel space.
            // A squared-distance score against a (0, 0) anchor is minimized
            // by the smallest-px+py placement, which is the classical
            // first-fit top-left scan. Post-placement centering (still run
            // unchanged) shifts the packed cluster to the bed center when
            // permitted. This preserves the old first-fit pack density for
            // the common case while giving wipe-tower-bearing plates the
            // clustering behavior the War Council audit identified as
            // defining parity failure #1 (C1).
            //
            // When the user has explicitly set align_center != (0.5, 0.5)
            // the anchor goes straight to their requested target — post-
            // centering cannot accurately re-target off-center alignment
            // because it moves the whole cluster as a rigid block and is
            // disabled on exclude-bearing plates; the scored scan is the
            // only code path that actually honors off-center alignment in
            // that case.
            auto resolve_anchor_px = [&](int plate_idx) -> std::pair<double, double> {
                for (const auto &ex : excludes) {
                    if (!ex.is_wipe_tower || ex.bed_idx != plate_idx) continue;
                    ExPolygon wt = ex.poly;
                    if (ex.rotation != 0.0) wt.rotate(ex.rotation);
                    wt.translate(ex.translation.x() - ebed.min.x(),
                                 ex.translation.y() - ebed.min.y());
                    BoundingBox wbb = get_extents(wt);
                    double cx_mm = (unscaled<double>(wbb.min.x())
                                    + unscaled<double>(wbb.max.x())) / 2.0;
                    double cy_mm = (unscaled<double>(wbb.min.y())
                                    + unscaled<double>(wbb.max.y())) / 2.0;
                    return {cx_mm / res, cy_mm / res};
                }
                // Hybrid anchor (default-align case only): corner seed
                // for the first item on a plate, bed center for items
                // 2..N. The corner seed avoids the central-trap
                // failure mode where center anchoring puts the first
                // item dead-center and fragments the surrounding
                // space into strips too thin for any subsequent item.
                // Subsequent items grow toward the bed center, which
                // gives a visually balanced cluster without sacrificing
                // tight packing.
                //
                // For non-default align_center, the user has explicitly
                // chosen a target position and the hybrid does not
                // apply — every item anchors at the chosen target so
                // the user's intent is honored. The scored scan is the
                // only code path that respects off-center alignment
                // (post-centering can't accurately re-target a rigid
                // cluster shift).
                double bw_mm = unscaled<double>(ebed.max.x() - ebed.min.x());
                double bh_mm = unscaled<double>(ebed.max.y() - ebed.min.y());
                constexpr double ALIGN_DEFAULT_EPS = 1e-9;
                bool center_default =
                    std::abs(params.align_center.x() - 0.5) < ALIGN_DEFAULT_EPS &&
                    std::abs(params.align_center.y() - 0.5) < ALIGN_DEFAULT_EPS;
                if (center_default && cluster_bb[plate_idx].empty()) {
                    return {0.0, 0.0};
                }
                return {bw_mm * params.align_center.x() / res,
                        bh_mm * params.align_center.y() / res};
            };

            // Try each existing plate, then overflow to a new one. For each
            // plate we run a SCORED scan: enumerate every coarse-stride
            // position across every rotation that passes position_clear,
            // score each by squared distance from the item bbox center to the
            // resolved anchor, keep the winner. Then refine within ±coarse of
            // the winner (also scored) to pick the pixel-precise position.
            // Items naturally cluster around the anchor, giving the default
            // arranger's wipe-tower-centric behavior (Sprint 1 Task 8 / the
            // defining wipe-tower-invisibility fix from the War Council).
            int max_plate = std::min(current_plate + 1, MAX_PLATES - 1);
            for (int plate_idx = 0; plate_idx <= max_plate && !placed; ++plate_idx) {
                ensure_plate_with_zones(plate_idx);

                auto anchor = resolve_anchor_px(plate_idx);
                const double ax = anchor.first;
                const double ay = anchor.second;

                const ClusterBB &cbb = cluster_bb[plate_idx];
                const int64_t cbb_area = cbb.area();
                const bool cbb_empty = cbb.empty();

                // Score a candidate placement as:
                //   primary   = growth in cluster bounding box area (int64)
                //   secondary = squared pixel distance of candidate bbox
                //               center from resolved anchor
                //
                // A rotated L-shape slotted into an existing L's concave notch
                // has its bbox already inside the cluster bbox, so delta_area
                // is zero and it wins over any non-interlocking position.
                // When the cluster is empty, all positions tie on delta_area
                // and the anchor tiebreaker decides (first item clusters
                // around the wipe tower; in the default case the anchor is
                // (0, 0) so first item lands top-left, preserving the old
                // first-fit baseline pack density).
                //
                // We return a pair so lexicographic comparison does the right
                // thing without numerical scaling games.
                auto score_at = [&](const RotCache &rc, int px, int py)
                    -> std::pair<int64_t, double>
                {
                    int nminx = px;
                    int nminy = py;
                    int nmaxx = px + rc.iw - 1;
                    int nmaxy = py + rc.ih - 1;
                    int64_t new_area;
                    if (cbb_empty) {
                        new_area = (int64_t)rc.iw * (int64_t)rc.ih;
                    } else {
                        int mnx = std::min(cbb.minx, nminx);
                        int mny = std::min(cbb.miny, nminy);
                        int mxx = std::max(cbb.maxx, nmaxx);
                        int mxy = std::max(cbb.maxy, nmaxy);
                        new_area = (int64_t)(mxx - mnx + 1) * (int64_t)(mxy - mny + 1);
                    }
                    int64_t delta = new_area - cbb_area;
                    double cx = (double)px + rc.iw * 0.5 - ax;
                    double cy = (double)py + rc.ih * 0.5 - ay;
                    double dist2 = cx * cx + cy * cy;
                    return {delta, dist2};
                };

                // Per-rotation coarse winner — NOT a single global winner.
                // The old code kept only the global coarse best, then refined
                // that rotation alone. That structurally prevented the refine
                // pass from discovering interlocking positions in rotations
                // that lost the coarse-grid lottery by a hair to a
                // non-interlocking rotation. For concave shapes where bbox is
                // insensitive to rotation (any L), the coarse winner is
                // effectively random across rotations and the single-rotation
                // refine misses obvious interlocks. Fix: remember the best
                // coarse (px, py) per rotation, then run refine for EVERY
                // rotation and keep the global best across all refine results.
                struct CoarseBest {
                    std::pair<int64_t, double> score{
                        std::numeric_limits<int64_t>::max(),
                        std::numeric_limits<double>::infinity()};
                    int px = -1;
                    int py = -1;
                };
                std::vector<CoarseBest> per_rot_coarse(rot_cache.size());

                // Coarse scan wrapped in a lambda so we can invoke it twice:
                // first with normal stride, and on failure retry once with
                // halved stride ("desperation pass") before falling through
                // to the next plate. The halved stride catches items whose
                // only valid slot falls between normal-stride samples —
                // typical for narrow gaps in a nearly-full plate where the
                // consolidation pass lives. Only runs on failure, so the
                // fast path pays nothing.
                auto run_coarse_scan = [&](int stride_div) -> bool {
                    per_rot_coarse.assign(rot_cache.size(), CoarseBest{});

                    for (size_t rci = 0; rci < rot_cache.size(); ++rci) {
                        const auto &rc = rot_cache[rci];
                        int stride = std::max(2, rc.coarse / stride_div);
                        int max_py = bh - rc.ih;
                        int max_px = bw - rc.iw;
                        if (max_py < 0 || max_px < 0) continue;

                        std::vector<int> pys;
                        for (int py = 0; py <= max_py; py += stride) pys.push_back(py);
                        if (pys.empty() || pys.back() != max_py) pys.push_back(max_py);

                        std::vector<int> pxs;
                        for (int px = 0; px <= max_px; px += stride) pxs.push_back(px);
                        if (pxs.empty() || pxs.back() != max_px) pxs.push_back(max_px);

                        CoarseBest &cb = per_rot_coarse[rci];
                        for (int py : pys) {
                            for (int px : pxs) {
                                if (!position_clear(plate_idx, rc, px, py)) continue;
                                auto s = score_at(rc, px, py);
                                if (s < cb.score) {
                                    cb.score = s;
                                    cb.px    = px;
                                    cb.py    = py;
                                }
                            }
                        }
                    }

                    for (const auto &cb : per_rot_coarse)
                        if (cb.px >= 0) return true;
                    return false;
                };

                // Normal stride first, then desperation retry on failure.
                if (!run_coarse_scan(1)) {
                    if (!run_coarse_scan(2)) continue;  // nothing fits on this plate
                }

                // Refine: for each rotation that had a coarse winner, scan
                // every pixel in a ±coarse window around that rotation's best
                // coarse position. Score every clear position with the same
                // score_at and track the GLOBAL best across all rotations.
                // This is the fix for bug #3 — an interlocking rotation that
                // lost the coarse lottery still gets a chance at refinement.
                std::pair<int64_t, double> refine_score{
                    std::numeric_limits<int64_t>::max(),
                    std::numeric_limits<double>::infinity()};
                int refine_px = -1, refine_py = -1;
                int refine_rc_idx = -1;

                for (size_t rci = 0; rci < rot_cache.size(); ++rci) {
                    const CoarseBest &cb = per_rot_coarse[rci];
                    if (cb.px < 0) continue;
                    const auto &rc = rot_cache[rci];
                    int ry0 = std::max(0, cb.py - rc.coarse);
                    int ry1 = std::min(bh - rc.ih, cb.py + rc.coarse);
                    int rx0 = std::max(0, cb.px - rc.coarse);
                    int rx1 = std::min(bw - rc.iw, cb.px + rc.coarse);

                    for (int py = ry0; py <= ry1; ++py) {
                        for (int px = rx0; px <= rx1; ++px) {
                            if (!position_clear(plate_idx, rc, px, py)) continue;
                            auto s = score_at(rc, px, py);
                            if (s < refine_score) {
                                refine_score = s;
                                refine_px    = px;
                                refine_py    = py;
                                refine_rc_idx = (int)rci;
                            }
                        }
                    }
                }

                if (refine_rc_idx < 0) continue;  // shouldn't happen but guard
                best_px     = refine_px;
                best_py     = refine_py;
                best_rot    = rot_cache[refine_rc_idx].rot;
                best_plate  = plate_idx;
                best_rc_idx = refine_rc_idx;
                placed      = true;
            }

            if (placed) {
                if (best_plate > current_plate)
                    current_plate = best_plate;

                const RotCache &best_rc = rot_cache[best_rc_idx];
                const BoundingBox &best_inflated_bb = best_rc.inflated_bb;

                // Stamp the concave bitmap into plate_items only. plate_excludes
                // is immutable after setup — items contribute to the items
                // obstacle map, not to the exclude map.
                stamp(plate_items[best_plate], wpr, bw, bh,
                      best_rc.bm, best_rc.iwpr, best_rc.iw, best_rc.ih,
                      best_px, best_py);

                // Extend the running cluster bbox with this placement's
                // pixel footprint so the next item's score_at sees the
                // updated obstacle extent. Uses the rotation cache bbox
                // dimensions, not the raw polygon bbox, so dilation is
                // included.
                {
                    ClusterBB &c = cluster_bb[best_plate];
                    int nminx = best_px;
                    int nminy = best_py;
                    int nmaxx = best_px + best_rc.iw - 1;
                    int nmaxy = best_py + best_rc.ih - 1;
                    if (c.empty()) {
                        c.minx = nminx;
                        c.miny = nminy;
                        c.maxx = nmaxx;
                        c.maxy = nmaxy;
                    } else {
                        c.minx = std::min(c.minx, nminx);
                        c.miny = std::min(c.miny, nminy);
                        c.maxx = std::max(c.maxx, nmaxx);
                        c.maxy = std::max(c.maxy, nmaxy);
                    }
                }

                // L12: use inflated polygon's actual bbox for translation.
                // The rasterizer places inflated_bb.min at pixel (best_px, best_py).
                // The polygon origin (0,0) is at -inflated_bb.min from that corner.
                double origin_x = best_px * res
                                  - unscaled<double>(best_inflated_bb.min.x())
                                  + unscaled<double>(ebed.min.x());
                double origin_y = best_py * res
                                  - unscaled<double>(best_inflated_bb.min.y())
                                  + unscaled<double>(ebed.min.y());

                item.translation = Vec2crd{scaled(origin_x), scaled(origin_y)};
                // Compose the nester's chosen rotation on top of the caller's
                // pre-rotation so the axis-align angle is not discarded.
                item.rotation = base_rot + best_rot;
                item.bed_idx = best_plate;
                item.itemid = item_sequence++;
                if (params.on_packed)
                    params.on_packed(item);
            } else {
                item.bed_idx = UNARRANGED;
            }
        }

        // ── Plate consolidation (backwards pass) ─────────────────────
        // After greedy placement, walk plates last-to-first and try to
        // migrate each item from a later plate back to an earlier plate
        // by finding any clear position there. Most plates after a
        // reasonable greedy pass have small holes that can swallow one
        // or two pieces from the next plate, eliminating the overflow.
        // The single-pass greedy can't see these holes because it only
        // visits each plate once and never revisits after later items
        // change the cluster shape.
        //
        // Cost: each successful migration triggers a re-stamp of the
        // source plate (rebuild plate_items[from_plate] by re-rasterizing
        // every surviving item). Most attempts fail at the cheap bbox
        // check or first position scan; the expensive case (a real
        // migration) is exactly when a plate is about to be eliminated.
        if (current_plate > 0 && !(params.stopcondition && params.stopcondition())) {
            int safety_iter = 0;
            bool any_migration = true;
            while (any_migration && safety_iter++ < 4) {
                any_migration = false;
                for (int from_plate = current_plate; from_plate > 0; --from_plate) {
                    if (params.stopcondition && params.stopcondition()) break;

                    std::vector<size_t> from_indices;
                    for (size_t i = 0; i < items.size(); ++i)
                        if (items[i].bed_idx == from_plate)
                            from_indices.push_back(i);
                    if (from_indices.empty()) continue;

                    bool plate_dirty = false;

                    for (size_t item_idx : from_indices) {
                        if (params.stopcondition && params.stopcondition()) break;
                        auto &it = items[item_idx];

                        const auto &rots = !it.allowed_rotations.empty()
                                           ? it.allowed_rotations : default_rotations;

                        // Caller's original pre-rotation, captured before the
                        // greedy loop ran. Mirror the greedy code path: pre-
                        // rotate base_shapes by base_rot, then compose
                        // rotations from `rots` ON TOP, and at commit write
                        // `item.rotation = base_rot + chosen_rot` so the
                        // caller's axis-align angle survives migration.
                        const double base_rot = base_rotations[item_idx];

                        ExPolygons base_shapes;
                        if (!it.concave_regions.empty())
                            base_shapes = it.concave_regions;
                        else
                            base_shapes.push_back(it.poly);
                        if (base_rot != 0.0)
                            for (ExPolygon &s : base_shapes) s.rotate(base_rot);

                        int pad_px = 0;
                        if (it.inflation > 0)
                            pad_px = std::max(1, (int)std::ceil(unscaled<double>(it.inflation) / res));
                        else if (params.min_obj_distance > 0)
                            pad_px = std::max(1, (int)std::ceil(unscaled<double>(params.min_obj_distance / 2) / res));

                        // Build rot_cache once per item (not per plate × rot
                        // as the previous code did). Each entry has the
                        // rasterized + optionally dilated bitmap pair, the
                        // dimensions, and the inflated bbox. The per-to_plate
                        // loop below scores all rotations against that plate
                        // and picks the best — mirroring the per-rotation
                        // refine pass in the main placement loop. The old
                        // code picked the first rotation that found any valid
                        // position, which is the same asymmetry BUG3 fixed in
                        // the placement loop earlier this session.
                        struct MigRot {
                            double rot;
                            int iw, ih, iwpr_item;
                            std::vector<uint64_t> ibm, hbm;
                            BoundingBox inflated_bb;
                        };
                        std::vector<MigRot> mig_cache;
                        mig_cache.reserve(rots.size());
                        for (double rot : rots) {
                            ExPolygons rshapes = base_shapes;
                            if (rot != 0.0)
                                for (ExPolygon &s : rshapes) s.rotate(rot);

                            int raw_iw, raw_ih, raw_iwpr;
                            auto raw_bm = rasterize(rshapes, res, bw, bh, raw_iw, raw_ih, raw_iwpr);
                            if (raw_bm.empty()) continue;

                            ExPolygon hull_expoly;
                            hull_expoly.contour = Geometry::convex_hull(rshapes);
                            if (hull_expoly.contour.empty()) continue;
                            int h_iw, h_ih, h_iwpr;
                            auto h_raw = rasterize(hull_expoly, res, bw, bh, h_iw, h_ih, h_iwpr);
                            if (h_raw.empty() || h_iw != raw_iw || h_ih != raw_ih) continue;

                            MigRot mr;
                            mr.rot = rot;
                            mr.iw = raw_iw;
                            mr.ih = raw_ih;
                            mr.iwpr_item = raw_iwpr;
                            mr.inflated_bb = get_extents(rshapes);
                            if (pad_px > 0) {
                                mr.iw = raw_iw + 2 * pad_px;
                                mr.ih = raw_ih + 2 * pad_px;
                                if (mr.iw > bw || mr.ih > bh) continue;
                                mr.iwpr_item = (mr.iw + 63) / 64;
                                mr.ibm = dilate_bitmap(raw_bm, raw_iwpr, raw_iw, raw_ih,
                                                       pad_px, mr.iwpr_item, mr.iw, mr.ih);
                                mr.hbm = dilate_bitmap(h_raw, raw_iwpr, raw_iw, raw_ih,
                                                       pad_px, mr.iwpr_item, mr.iw, mr.ih);
                                coord_t pad_sc = scaled(pad_px * res);
                                mr.inflated_bb.min -= Vec2crd(pad_sc, pad_sc);
                                mr.inflated_bb.max += Vec2crd(pad_sc, pad_sc);
                            } else {
                                if (mr.iw > bw || mr.ih > bh) continue;
                                mr.ibm = std::move(raw_bm);
                                mr.hbm = std::move(h_raw);
                            }
                            mig_cache.push_back(std::move(mr));
                        }

                        if (mig_cache.empty()) continue;

                        bool migrated = false;
                        for (int to_plate = 0; to_plate < from_plate && !migrated; ++to_plate) {
                            // Per-to_plate: score every rotation's best
                            // position (coarse + refine), keep the global
                            // best-scoring (rotation, px, py). Commit exactly
                            // that winner.

                            // Destination cluster bbox snapshot for scoring —
                            // same for every rotation tried on this plate.
                            const ClusterBB &dst_cbb = cluster_bb[to_plate];
                            bool dst_empty = dst_cbb.empty();
                            int64_t dst_area = dst_empty ? 0 :
                                (int64_t)(dst_cbb.maxx - dst_cbb.minx + 1) *
                                (int64_t)(dst_cbb.maxy - dst_cbb.miny + 1);

                            int     best_rci   = -1;
                            int64_t best_score = std::numeric_limits<int64_t>::max();
                            int     best_px    = -1;
                            int     best_py    = -1;

                            for (size_t rci = 0; rci < mig_cache.size(); ++rci) {
                                const MigRot &mr = mig_cache[rci];
                                int max_py = bh - mr.ih;
                                int max_px = bw - mr.iw;
                                if (max_py < 0 || max_px < 0) continue;

                                int cstride = std::clamp(std::min(mr.iw, mr.ih) / 8, 4, 64);

                                auto position_clear_dst = [&](int px, int py) -> bool {
                                    if (collides(plate_items[to_plate], wpr, bw, bh,
                                                 mr.ibm, mr.iwpr_item, mr.iw, mr.ih, px, py)) return false;
                                    if (collides(plate_excludes[to_plate], wpr, bw, bh,
                                                 mr.hbm, mr.iwpr_item, mr.iw, mr.ih, px, py)) return false;
                                    return true;
                                };

                                auto score_dst = [&](int px, int py) -> int64_t {
                                    int nminx = px;
                                    int nminy = py;
                                    int nmaxx = px + mr.iw - 1;
                                    int nmaxy = py + mr.ih - 1;
                                    int64_t new_area;
                                    if (dst_empty) {
                                        new_area = (int64_t)mr.iw * (int64_t)mr.ih;
                                    } else {
                                        int mnx = std::min(dst_cbb.minx, nminx);
                                        int mny = std::min(dst_cbb.miny, nminy);
                                        int mxx = std::max(dst_cbb.maxx, nmaxx);
                                        int mxy = std::max(dst_cbb.maxy, nmaxy);
                                        new_area = (int64_t)(mxx - mnx + 1) *
                                                   (int64_t)(mxy - mny + 1);
                                    }
                                    return new_area - dst_area;
                                };

                                // Coarse scan with explicit boundary positions.
                                std::vector<int> pys;
                                for (int py = 0; py <= max_py; py += cstride) pys.push_back(py);
                                if (pys.empty() || pys.back() != max_py) pys.push_back(max_py);
                                std::vector<int> pxs;
                                for (int px = 0; px <= max_px; px += cstride) pxs.push_back(px);
                                if (pxs.empty() || pxs.back() != max_px) pxs.push_back(max_px);

                                int64_t coarse_score = std::numeric_limits<int64_t>::max();
                                int coarse_px = -1, coarse_py = -1;
                                for (int py : pys) {
                                    for (int px : pxs) {
                                        if (!position_clear_dst(px, py)) continue;
                                        int64_t s = score_dst(px, py);
                                        if (s < coarse_score) {
                                            coarse_score = s;
                                            coarse_px    = px;
                                            coarse_py    = py;
                                        }
                                    }
                                }
                                if (coarse_px < 0) continue;

                                // Refine ±cstride around the coarse winner.
                                int ry0 = std::max(0, coarse_py - cstride);
                                int ry1 = std::min(max_py, coarse_py + cstride);
                                int rx0 = std::max(0, coarse_px - cstride);
                                int rx1 = std::min(max_px, coarse_px + cstride);
                                int64_t rot_score = coarse_score;
                                int rot_px = coarse_px;
                                int rot_py = coarse_py;
                                for (int py = ry0; py <= ry1; ++py) {
                                    for (int px = rx0; px <= rx1; ++px) {
                                        if (!position_clear_dst(px, py)) continue;
                                        int64_t s = score_dst(px, py);
                                        if (s < rot_score) {
                                            rot_score = s;
                                            rot_px = px;
                                            rot_py = py;
                                        }
                                    }
                                }

                                if (rot_score < best_score) {
                                    best_score = rot_score;
                                    best_rci   = (int)rci;
                                    best_px    = rot_px;
                                    best_py    = rot_py;
                                }
                            }

                            if (best_rci < 0) continue;  // no rotation had a valid position on this plate

                            const MigRot &best = mig_cache[best_rci];

                            // Commit the migration with the best-scoring
                            // rotation's bitmap and position.
                            stamp(plate_items[to_plate], wpr, bw, bh,
                                  best.ibm, best.iwpr_item, best.iw, best.ih,
                                  best_px, best_py);

                            // Update destination cluster bbox.
                            {
                                ClusterBB &c = cluster_bb[to_plate];
                                int nminx = best_px, nminy = best_py;
                                int nmaxx = best_px + best.iw - 1;
                                int nmaxy = best_py + best.ih - 1;
                                if (c.empty()) {
                                    c.minx = nminx; c.miny = nminy;
                                    c.maxx = nmaxx; c.maxy = nmaxy;
                                } else {
                                    c.minx = std::min(c.minx, nminx);
                                    c.miny = std::min(c.miny, nminy);
                                    c.maxx = std::max(c.maxx, nmaxx);
                                    c.maxy = std::max(c.maxy, nmaxy);
                                }
                            }

                            // Update item state.
                            double origin_x = best_px * res
                                              - unscaled<double>(best.inflated_bb.min.x())
                                              + unscaled<double>(ebed.min.x());
                            double origin_y = best_py * res
                                              - unscaled<double>(best.inflated_bb.min.y())
                                              + unscaled<double>(ebed.min.y());
                            it.translation = Vec2crd{scaled(origin_x), scaled(origin_y)};
                            // Compose: base_rot is the caller's original
                            // pre-rotation; best.rot is the consolidation
                            // pass's chosen alternative orientation. Mirrors
                            // the greedy commit at line ~638.
                            it.rotation    = base_rot + best.rot;
                            it.bed_idx     = to_plate;

                            migrated = true;
                            any_migration = true;
                            plate_dirty = true;
                        }
                    }

                    // Rebuild the source plate's bitmap by re-stamping
                    // every surviving item. We can't selectively un-stamp
                    // because the composite is OR'd (no per-item ownership).
                    if (plate_dirty) {
                        std::fill(plate_items[from_plate].begin(),
                                  plate_items[from_plate].end(), 0);
                        cluster_bb[from_plate] = ClusterBB{
                            std::numeric_limits<int>::max(),
                            std::numeric_limits<int>::max(),
                            std::numeric_limits<int>::min(),
                            std::numeric_limits<int>::min()};
                        for (auto &it2 : items) {
                            if (it2.bed_idx != from_plate) continue;
                            ExPolygons base = !it2.concave_regions.empty()
                                               ? it2.concave_regions
                                               : ExPolygons{it2.poly};
                            if (it2.rotation != 0.0)
                                for (ExPolygon &s : base) s.rotate(it2.rotation);

                            int pad2 = 0;
                            if (it2.inflation > 0)
                                pad2 = std::max(1, (int)std::ceil(unscaled<double>(it2.inflation) / res));
                            else if (params.min_obj_distance > 0)
                                pad2 = std::max(1, (int)std::ceil(unscaled<double>(params.min_obj_distance / 2) / res));

                            int r_iw, r_ih, r_iwpr;
                            auto r_bm = rasterize(base, res, bw, bh, r_iw, r_ih, r_iwpr);
                            if (r_bm.empty()) continue;

                            int s_iw = r_iw, s_ih = r_ih, s_iwpr = r_iwpr;
                            std::vector<uint64_t> s_bm;
                            BoundingBox sbb = get_extents(base);
                            if (pad2 > 0) {
                                s_iw = r_iw + 2 * pad2;
                                s_ih = r_ih + 2 * pad2;
                                s_iwpr = (s_iw + 63) / 64;
                                s_bm = dilate_bitmap(r_bm, r_iwpr, r_iw, r_ih,
                                                     pad2, s_iwpr, s_iw, s_ih);
                                coord_t pad_sc2 = scaled(pad2 * res);
                                sbb.min -= Vec2crd(pad_sc2, pad_sc2);
                                sbb.max += Vec2crd(pad_sc2, pad_sc2);
                            } else {
                                s_bm = std::move(r_bm);
                            }

                            // Recover the pixel position from the world
                            // translation: invert the placement commit math.
                            double ox = unscaled<double>(it2.translation.x())
                                        + unscaled<double>(sbb.min.x())
                                        - unscaled<double>(ebed.min.x());
                            double oy = unscaled<double>(it2.translation.y())
                                        + unscaled<double>(sbb.min.y())
                                        - unscaled<double>(ebed.min.y());
                            int spx = (int)std::round(ox / res);
                            int spy = (int)std::round(oy / res);
                            if (spx < 0 || spy < 0) continue;
                            if (spx + s_iw > bw || spy + s_ih > bh) continue;

                            stamp(plate_items[from_plate], wpr, bw, bh,
                                  s_bm, s_iwpr, s_iw, s_ih, spx, spy);

                            ClusterBB &c = cluster_bb[from_plate];
                            int nminx = spx, nminy = spy;
                            int nmaxx = spx + s_iw - 1, nmaxy = spy + s_ih - 1;
                            if (c.empty()) {
                                c.minx = nminx; c.miny = nminy;
                                c.maxx = nmaxx; c.maxy = nmaxy;
                            } else {
                                c.minx = std::min(c.minx, nminx);
                                c.miny = std::min(c.miny, nminy);
                                c.maxx = std::max(c.maxx, nmaxx);
                                c.maxy = std::max(c.maxy, nmaxy);
                            }
                        }
                    }

                    // Shrink current_plate if from_plate emptied out.
                    bool plate_empty = true;
                    for (auto &it3 : items) {
                        if (it3.bed_idx == from_plate) { plate_empty = false; break; }
                    }
                    if (plate_empty && from_plate == current_plate)
                        current_plate--;
                }
            }
        }

        // Post-placement centering: shift placed items so the cluster
        // center lands on align_center (default: bed center).
        // Now runs on every plate including those with excludes — the
        // safety check inside the loop validates the shift before
        // committing it (refuses the move if any item would collide
        // with an exclude after the shift).
        //
        // Honor stopcondition here as well: when the user cancels mid-arrange,
        // the placement loop breaks out of its own scan, but the centering
        // pass would otherwise keep moving items that were already written
        // with their pre-center coordinates. Bail out cleanly so cancellation
        // leaves placed items at their raw scan positions instead of a
        // half-shifted cluster.
        if (params.do_final_align && !(params.stopcondition && params.stopcondition())) {
            // Build set of plates that have any exclude
            std::vector<bool> plate_has_exclude(current_plate + 1, false);
            for (auto &ex : excludes) {
                if (ex.bed_idx >= 0 && ex.bed_idx <= current_plate)
                    plate_has_exclude[ex.bed_idx] = true;
            }
            if (!params.excluded_regions.empty()) {
                for (int p = 0; p <= current_plate; ++p)
                    plate_has_exclude[p] = true;
            }

            for (int plate = 0; plate <= current_plate; ++plate) {
                if (params.stopcondition && params.stopcondition()) break;

                BoundingBox cluster_bb;
                bool has_item = false;
                for (auto &item : items) {
                    if (item.bed_idx != plate) continue;
                    ExPolygon placed = item.poly;
                    if (item.rotation != 0.0) placed.rotate(item.rotation);
                    placed.translate(item.translation.x(), item.translation.y());
                    BoundingBox ibb = get_extents(placed);
                    if (!has_item) { cluster_bb = ibb; has_item = true; }
                    else cluster_bb.merge(ibb);
                }
                if (!has_item) continue;

                double tx = bed.min.x() + (bed.max.x() - bed.min.x()) * params.align_center.x();
                double ty = bed.min.y() + (bed.max.y() - bed.min.y()) * params.align_center.y();
                double cx = (cluster_bb.min.x() + cluster_bb.max.x()) / 2.0;
                double cy = (cluster_bb.min.y() + cluster_bb.max.y()) / 2.0;
                coord_t dx = (coord_t)(tx - cx);
                coord_t dy = (coord_t)(ty - cy);

                if (cluster_bb.min.x() + dx < bed.min.x())
                    dx = bed.min.x() - cluster_bb.min.x();
                if (cluster_bb.max.x() + dx > bed.max.x())
                    dx = bed.max.x() - cluster_bb.max.x();
                if (cluster_bb.min.y() + dy < bed.min.y())
                    dy = bed.min.y() - cluster_bb.min.y();
                if (cluster_bb.max.y() + dy > bed.max.y())
                    dy = bed.max.y() - cluster_bb.max.y();

                // Safety check on exclude-bearing plates: would the shift
                // cause any item to collide with an exclude? If so, refuse
                // to move (dx=dy=0). Previously this whole branch was
                // skipped on exclude plates; now we attempt the centering
                // and only abort if the result would actually overlap.
                if (plate_has_exclude[plate] && (dx != 0 || dy != 0)) {
                    bool collision_after_shift = false;
                    for (const auto &ex : excludes) {
                        if (ex.bed_idx != plate) continue;
                        ExPolygon ex_poly = ex.poly;
                        if (ex.rotation != 0.0) ex_poly.rotate(ex.rotation);
                        ex_poly.translate(ex.translation.x(), ex.translation.y());
                        for (const auto &item : items) {
                            if (item.bed_idx != plate) continue;
                            ExPolygon shifted = item.poly;
                            if (item.rotation != 0.0) shifted.rotate(item.rotation);
                            shifted.translate(item.translation.x() + dx,
                                              item.translation.y() + dy);
                            if (!intersection_ex(ExPolygons{ex_poly},
                                                 ExPolygons{shifted}).empty()) {
                                collision_after_shift = true;
                                break;
                            }
                        }
                        if (collision_after_shift) break;
                    }
                    if (collision_after_shift) { dx = 0; dy = 0; }
                }

                for (auto &item : items) {
                    if (item.bed_idx != plate) continue;
                    item.translation += Vec2crd(dx, dy);
                }
            }
        }

        }; // end of do_one_pass lambda

        // ── Run attempts ─────────────────────────────────────────────
        for (int attempt = 0; attempt < num_attempts; ++attempt) {
            // Stopcondition check is gated on `attempt > 0`. Calling
            // params.stopcondition() may have side effects (the test
            // suite uses a counter-based predicate, and the user's
            // GUI cancel button increments a token). The first attempt
            // must always run untouched so its placement loop's own
            // internal stopcondition checks have the budget the
            // caller intended.
            if (attempt > 0 && params.stopcondition && params.stopcondition()) break;

            if (attempt > 0) {
                reset_for_restart();
                // Perturb: shuffle within priority buckets so the
                // priority contract is preserved (high priority items
                // still placed first) but the within-bucket order
                // varies between attempts.
                order = base_order;
                std::mt19937 rng(0x5EED1234u + (uint32_t)attempt);
                for (size_t i = 0; i < order.size();) {
                    size_t j = i + 1;
                    while (j < order.size() &&
                           items[order[j]].priority == items[order[i]].priority)
                        ++j;
                    std::shuffle(order.begin() + i, order.begin() + j, rng);
                    i = j;
                }
            }

            do_one_pass();

            AttemptScore sc = score_current();
            // Lex order: more placed > fewer placed; then fewer plates
            // > more plates; then smaller cluster area > larger.
            bool is_better = !best.valid ||
                             sc.placed > best.placed ||
                             (sc.placed == best.placed && sc.max_bed < best.max_bed) ||
                             (sc.placed == best.placed && sc.max_bed == best.max_bed
                              && sc.area < best.total_area);
            if (is_better) save_best(sc);

            // Early termination: a fully-placed result on a single
            // plate is the optimum we'd nominally try to improve.
            // Save the rest of the budget.
            if (best.placed == (int)items.size() && best.max_bed == 0) break;
        }

        // Restore the best result into items[].
        if (best.valid) {
            for (size_t i = 0; i < items.size(); ++i) {
                items[i].translation = best.trans[i];
                items[i].rotation    = best.rots[i];
                items[i].bed_idx     = best.beds[i];
                items[i].itemid      = best.iids[i];
            }
        }
    }

#ifdef BITMAP_NESTER_TESTING
public:
#else
private:
#endif
    // Dilate a bitmap by pad pixels in all directions.
    // The output bitmap is (iw + 2*pad) × (ih + 2*pad), pre-allocated by caller.
    // Algorithm: horizontal spread via word-level shift-and-OR, then vertical
    // propagation. O(oh × owpr × pad) — no per-pixel loops, no Clipper.
    static std::vector<uint64_t> dilate_bitmap(
        const std::vector<uint64_t> &src, int swpr, int sw, int sh,
        int pad, int dwpr, int dw, int dh)
    {
        std::vector<uint64_t> dst((size_t)dwpr * dh, 0);

        // Step 1: copy source into center of destination (offset by pad,pad)
        // Word-level copy with bit-shift for the horizontal offset.
        {
            int word_off = pad / 64;
            int bit_off  = pad % 64;
            for (int y = 0; y < sh; ++y) {
                size_t srow = (size_t)y * swpr;
                size_t drow = (size_t)(y + pad) * dwpr + word_off;
                if (bit_off == 0) {
                    for (int w = 0; w < swpr; ++w)
                        dst[drow + w] |= src[srow + w];
                } else {
                    uint64_t carry = 0;
                    for (int w = 0; w < swpr; ++w) {
                        uint64_t val = src[srow + w];
                        dst[drow + w] |= (val << bit_off) | carry;
                        carry = val >> (64 - bit_off);
                    }
                    if (carry && (word_off + swpr) < dwpr)
                        dst[drow + swpr] |= carry;
                }
            }
        }

        // Step 2: horizontal dilation — spread each row left and right by pad.
        // Doubling strategy: spread by 1, then OR-shift by 2, 4, 8...
        // O(dh × dwpr × log2(pad)) instead of O(dh × dwpr × pad).
        for (int y = 0; y < dh; ++y) {
            size_t row_off = (size_t)y * dwpr;
            // Spread right by pad
            int remaining = pad;
            for (int step = 1; remaining > 0; step *= 2) {
                int shift = std::min(step, remaining);
                int word_shift = shift / 64;
                int bit_shift  = shift % 64;
                // OR row with itself shifted right by 'shift' bits
                for (int w = dwpr - 1; w >= 0; --w) {
                    int sw = w - word_shift;
                    uint64_t lo = (sw >= 0) ? dst[row_off + sw] : 0;
                    uint64_t hi = (sw - 1 >= 0) ? dst[row_off + sw - 1] : 0;
                    uint64_t shifted = (bit_shift == 0) ? lo
                        : (lo << bit_shift) | (hi >> (64 - bit_shift));
                    dst[row_off + w] |= shifted;
                }
                remaining -= shift;
            }
            // Spread left by pad
            remaining = pad;
            for (int step = 1; remaining > 0; step *= 2) {
                int shift = std::min(step, remaining);
                int word_shift = shift / 64;
                int bit_shift  = shift % 64;
                for (int w = 0; w < dwpr; ++w) {
                    int sw = w + word_shift;
                    uint64_t lo = (sw < dwpr) ? dst[row_off + sw] : 0;
                    uint64_t hi = (sw + 1 < dwpr) ? dst[row_off + sw + 1] : 0;
                    uint64_t shifted = (bit_shift == 0) ? lo
                        : (lo >> bit_shift) | (hi << (64 - bit_shift));
                    dst[row_off + w] |= shifted;
                }
                remaining -= shift;
            }
        }

        // Step 3: vertical dilation — spread columns up and down by pad.
        // Same doubling strategy: O(dh × dwpr × log2(pad)).
        {
            int remaining = pad;
            for (int step = 1; remaining > 0; step *= 2) {
                int shift = std::min(step, remaining);
                // Spread down by 'shift' rows
                for (int y = dh - 1 - shift; y >= 0; --y) {
                    size_t src_off = (size_t)y * dwpr;
                    size_t dst_off = (size_t)(y + shift) * dwpr;
                    for (int w = 0; w < dwpr; ++w)
                        dst[dst_off + w] |= dst[src_off + w];
                }
                // Spread up by 'shift' rows
                for (int y = shift; y < dh; ++y) {
                    size_t src_off = (size_t)y * dwpr;
                    size_t dst_off = (size_t)(y - shift) * dwpr;
                    for (int w = 0; w < dwpr; ++w)
                        dst[dst_off + w] |= dst[src_off + w];
                }
                remaining -= shift;
            }
        }

        return dst;
    }

    // Collision check: does item bitmap at (px, py) overlap the plate bitmap?
    // L9: px and py must be non-negative (enforced by the placement scan)
    static bool collides(const std::vector<uint64_t> &plate,
                         int wpr, int bw, int bh,
                         const std::vector<uint64_t> &item_bm,
                         int iwpr, int iw, int ih,
                         int px, int py)
    {
        assert(px >= 0 && py >= 0);
        for (int y = 0; y < ih; ++y) {
            int by = py + y;
            if (by >= bh) break;
            for (int wx = 0; wx < iwpr; ++wx) {
                uint64_t word = item_bm[(size_t)y * iwpr + wx];
                if (word == 0) continue;

                int bit_offset = px + wx * 64;
                if (bit_offset >= bw) continue;
                int bed_word = bit_offset / 64;
                int shift = bit_offset % 64;

                if (bed_word < wpr) {
                    if (plate[(size_t)by * wpr + bed_word] & (word << shift))
                        return true;
                }
                if (shift > 0 && bed_word + 1 < wpr) {
                    uint64_t shifted = word >> (64 - shift);
                    if (shifted && (plate[(size_t)by * wpr + bed_word + 1] & shifted))
                        return true;
                }
            }
        }
        return false;
    }

    // Stamp item bitmap onto plate bitmap at position (px, py)
    static void stamp(std::vector<uint64_t> &plate,
                      int wpr, int bw, int bh,
                      const std::vector<uint64_t> &item_bm,
                      int iwpr, int iw, int ih,
                      int px, int py)
    {
        assert(px >= 0 && py >= 0);
        for (int y = 0; y < ih; ++y) {
            int by = py + y;
            if (by >= bh) break;
            for (int wx = 0; wx < iwpr; ++wx) {
                uint64_t word = item_bm[(size_t)y * iwpr + wx];
                if (word == 0) continue;

                int bit_offset = px + wx * 64;
                if (bit_offset >= bw) continue;
                int bed_word = bit_offset / 64;
                int shift = bit_offset % 64;

                if (bed_word < wpr) {
                    plate[(size_t)by * wpr + bed_word] |= (word << shift);
                }
                if (shift > 0 && bed_word + 1 < wpr) {
                    uint64_t hi = word >> (64 - shift);
                    if (hi)
                        plate[(size_t)by * wpr + bed_word + 1] |= hi;
                }
            }
        }
    }

    // Rasterize an ExPolygon into a bitmap covering its bounding box.
    static std::vector<uint64_t> rasterize(const ExPolygon &poly,
                                           double res_mm,
                                           int max_w, int max_h,
                                           int &iw, int &ih, int &iwpr)
    {
        BoundingBox bb = get_extents(poly);
        if (bb.min.x() >= bb.max.x() || bb.min.y() >= bb.max.y()) {
            iw = ih = iwpr = 0;
            return {};
        }

        iw = std::min(max_w, std::max(1, (int)std::ceil(unscaled<double>(bb.size().x()) / res_mm)));
        ih = std::min(max_h, std::max(1, (int)std::ceil(unscaled<double>(bb.size().y()) / res_mm)));
        iwpr = (iw + 63) / 64;

        std::vector<uint64_t> bm((size_t)iwpr * ih, 0);

        scanline_fill(poly.contour, bb, res_mm, iw, ih, iwpr, bm, true);

        for (auto &hole : poly.holes)
            scanline_fill(hole, bb, res_mm, iw, ih, iwpr, bm, false);

        return bm;
    }

    // Rasterize multiple ExPolygons into a single bitmap covering their union
    // bbox. Each region's contour is filled; each region's holes are cleared.
    // Disconnected islands pack into one bitmap so downstream collision and
    // dilation code can treat the whole silhouette as a single stamp.
    static std::vector<uint64_t> rasterize(const ExPolygons &polys,
                                           double res_mm,
                                           int max_w, int max_h,
                                           int &iw, int &ih, int &iwpr)
    {
        if (polys.empty()) {
            iw = ih = iwpr = 0;
            return {};
        }

        BoundingBox bb = get_extents(polys);
        if (bb.min.x() >= bb.max.x() || bb.min.y() >= bb.max.y()) {
            iw = ih = iwpr = 0;
            return {};
        }

        iw = std::min(max_w, std::max(1, (int)std::ceil(unscaled<double>(bb.size().x()) / res_mm)));
        ih = std::min(max_h, std::max(1, (int)std::ceil(unscaled<double>(bb.size().y()) / res_mm)));
        iwpr = (iw + 63) / 64;

        std::vector<uint64_t> bm((size_t)iwpr * ih, 0);

        for (const ExPolygon &poly : polys) {
            scanline_fill(poly.contour, bb, res_mm, iw, ih, iwpr, bm, true);
            for (auto &hole : poly.holes)
                scanline_fill(hole, bb, res_mm, iw, ih, iwpr, bm, false);
        }

        return bm;
    }

    // Scanline rasterizer for a single polygon ring.
    // L13: xs vector hoisted outside y loop to avoid per-scanline heap allocation.
    // L15: word-at-a-time span fill for runs wider than a single word.
    static void scanline_fill(const Polygon &poly,
                              const BoundingBox &bb,
                              double res_mm,
                              int iw, int ih, int iwpr,
                              std::vector<uint64_t> &bm,
                              bool set)
    {
        if (poly.points.size() < 3) return;

        double ox = unscaled<double>(bb.min.x());
        double oy = unscaled<double>(bb.min.y());

        std::vector<double> xs;
        for (int y = 0; y < ih; ++y) {
            double scan_y = oy + (y + 0.5) * res_mm;

            xs.clear();
            for (size_t i = 0; i < poly.points.size(); ++i) {
                size_t j = (i + 1) % poly.points.size();
                double y0 = unscaled<double>(poly.points[i].y());
                double y1 = unscaled<double>(poly.points[j].y());

                if ((y0 <= scan_y && y1 > scan_y) || (y1 <= scan_y && y0 > scan_y)) {
                    double x0 = unscaled<double>(poly.points[i].x());
                    double x1 = unscaled<double>(poly.points[j].x());
                    double t = (scan_y - y0) / (y1 - y0);
                    xs.push_back(x0 + t * (x1 - x0));
                }
            }

            // L3: defend against odd intersection count from degenerate geometry
            if (xs.size() % 2 != 0) xs.pop_back();

            std::sort(xs.begin(), xs.end());

            for (size_t k = 0; k + 1 < xs.size(); k += 2) {
                // Conservative floor/floor fill rule: pixel i is filled when
                // i <= floor((xs[k+1]-ox)/res_mm). When xs[k+1] lands exactly
                // on a pixel boundary this includes one extra pixel at that
                // boundary — a deliberate off-by-one toward MORE material,
                // consistent with the rasterizer being fuzzy at sub-pixel
                // scales anyway. Fine for concave nesting because the cost
                // is a ~res_mm fringe on one edge, not a closed notch.
                int x_start = std::max(0,      (int)std::floor((xs[k]     - ox) / res_mm));
                int x_end   = std::min(iw - 1, (int)std::floor((xs[k + 1] - ox) / res_mm));

                if (x_start > x_end) continue;

                // L15: word-at-a-time span fill
                int w0 = x_start / 64;
                int b0 = x_start % 64;
                int w1 = x_end / 64;
                int b1 = x_end % 64;

                if (w0 == w1) {
                    // Span fits in a single word
                    uint64_t mask = ((b1 - b0 + 1) == 64)
                        ? ~uint64_t(0)
                        : ((uint64_t(1) << (b1 - b0 + 1)) - 1) << b0;
                    size_t idx = (size_t)y * iwpr + w0;
                    if (set) bm[idx] |= mask;
                    else     bm[idx] &= ~mask;
                } else {
                    // First partial word
                    uint64_t head = ~uint64_t(0) << b0;
                    size_t idx = (size_t)y * iwpr + w0;
                    if (set) bm[idx] |= head;
                    else     bm[idx] &= ~head;

                    // Full middle words
                    for (int w = w0 + 1; w < w1; ++w) {
                        idx = (size_t)y * iwpr + w;
                        if (set) bm[idx] = ~uint64_t(0);
                        else     bm[idx] = 0;
                    }

                    // Last partial word
                    uint64_t tail = (b1 == 63)
                        ? ~uint64_t(0)
                        : (uint64_t(1) << (b1 + 1)) - 1;
                    idx = (size_t)y * iwpr + w1;
                    if (set) bm[idx] |= tail;
                    else     bm[idx] &= ~tail;
                }
            }
        }
    }
};

}} // namespace Slic3r::arrangement
