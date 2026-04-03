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
#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

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

        // Sort items by priority (descending), then by area (descending)
        std::vector<size_t> order(items.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (items[a].priority != items[b].priority)
                return items[a].priority > items[b].priority;
            return std::abs(items[a].poly.area()) > std::abs(items[b].poly.area());
        });

        // Build rotation list for items that allow any rotation
        std::vector<double> default_rotations;
        if (params.allow_rotations) {
            for (int i = 0; i < 4; ++i)
                default_rotations.push_back(i * M_PI / 2.0);
        } else {
            default_rotations.push_back(0.0);
        }

        constexpr int MAX_PLATES = 36;

        // Per-plate forbidden bitmaps
        std::vector<std::vector<uint64_t>> plates;
        auto ensure_plate = [&](int idx) {
            while ((int)plates.size() <= idx) {
                plates.emplace_back((size_t)wpr * bh, uint64_t(0));
            }
        };
        ensure_plate(0);

        // Stamp an ExPolygon (in scaled bed-relative coords) onto a plate bitmap
        auto stamp_poly = [&](int plate_idx, const ExPolygon &poly_bed_rel) {
            ensure_plate(plate_idx);
            int iw, ih, iwpr;
            auto bm = rasterize(poly_bed_rel, res, bw, bh, iw, ih, iwpr);
            if (bm.empty()) return;
            BoundingBox pbb = get_extents(poly_bed_rel);
            int px = std::max(0, (int)(unscaled<double>(pbb.min.x()) / res));
            int py = std::max(0, (int)(unscaled<double>(pbb.min.y()) / res));
            stamp(plates[plate_idx], wpr, bw, bh, bm, iwpr, iw, ih, px, py);
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
            stamp_poly(plate, epoly);
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
            stamp(plates[0], wpr, bw, bh, zbm, ziwpr, ziw, zih, zpx, zpy);
            bed_zones.push_back({std::move(zbm), ziw, zih, ziwpr, zpx, zpy});
        }

        // Patch ensure_plate to stamp bed exclusion zones onto new plates
        auto ensure_plate_with_zones = [&](int idx) {
            int old_count = (int)plates.size();
            ensure_plate(idx);
            for (int p = old_count; p <= idx; ++p) {
                for (auto &z : bed_zones)
                    stamp(plates[p], wpr, bw, bh, z.bm, z.iwpr, z.iw, z.ih, z.px, z.py);
            }
        };

        int current_plate = 0;
        for (auto &ex : excludes) {
            if (ex.bed_idx > 0)
                current_plate = std::max(current_plate, ex.bed_idx);
        }

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

            // L14: pre-rasterize all rotations for this item so we don't
            // re-rasterize on every plate retry. For 75 items × 4 rotations × 6
            // plates, this cuts rasterize calls from ~1000 to ~300.
            struct RotCache {
                double rot;
                std::vector<uint64_t> bm;
                BoundingBox inflated_bb;
                int iw, ih, iwpr, coarse;
            };
            std::vector<RotCache> rot_cache;
            rot_cache.reserve(rots.size());

            for (double rot : rots) {
                ExPolygon rpoly = item.poly;
                if (rot != 0.0) rpoly.rotate(rot);

                // Rasterize the un-inflated polygon, then dilate the bitmap.
                // This avoids Clipper's offset_ex entirely in the hot path.
                int raw_iw, raw_ih, raw_iwpr;
                auto raw_bm = rasterize(rpoly, res, bw, bh, raw_iw, raw_ih, raw_iwpr);
                if (raw_bm.empty()) continue;

                // Dilate bitmap by pad_px in all directions
                int iw = raw_iw, ih = raw_ih, iwpr_item = raw_iwpr;
                std::vector<uint64_t> ibm;
                BoundingBox inflated_bb = get_extents(rpoly);

                if (pad_px > 0) {
                    iw = raw_iw + 2 * pad_px;
                    ih = raw_ih + 2 * pad_px;
                    if (iw > bw || ih > bh) continue;
                    iwpr_item = (iw + 63) / 64;
                    ibm = dilate_bitmap(raw_bm, raw_iwpr, raw_iw, raw_ih,
                                        pad_px, iwpr_item, iw, ih);
                    // Adjust inflated bbox to account for dilation
                    coord_t pad_sc = scaled(pad_px * res);
                    inflated_bb.min -= Vec2crd(pad_sc, pad_sc);
                    inflated_bb.max += Vec2crd(pad_sc, pad_sc);
                } else {
                    ibm = std::move(raw_bm);
                    if (iw > bw || ih > bh) continue;
                }

                int coarse = std::clamp(std::min(iw, ih) / 8, 8, 128);
                rot_cache.push_back({rot, std::move(ibm), inflated_bb,
                                     iw, ih, iwpr_item, coarse});
            }

            bool placed = false;
            int best_px = 0, best_py = 0;
            double best_rot = 0.0;
            int best_plate = -1;
            int best_iw = 0, best_ih = 0, best_iwpr = 0;
            std::vector<uint64_t> best_bm;
            BoundingBox best_inflated_bb;

            // Try each existing plate, then overflow to a new one
            int max_plate = std::min(current_plate + 1, MAX_PLATES - 1);
            for (int plate_idx = 0; plate_idx <= max_plate && !placed; ++plate_idx) {
                ensure_plate_with_zones(plate_idx);

                for (auto &rc : rot_cache) {
                    int found_px = -1, found_py = -1;

                    // Phase 1: coarse scan (always includes boundary positions)
                    {
                        int max_py = bh - rc.ih;
                        int max_px = bw - rc.iw;
                        for (int py = 0; py <= max_py; py += rc.coarse) {
                            int test_py = py;
                            for (int pass_y = 0; pass_y < 2; ++pass_y) {
                                if (pass_y == 1) {
                                    // On second pass, test the boundary row if stride skipped it
                                    if (py + rc.coarse > max_py && py < max_py)
                                        test_py = max_py;
                                    else break;
                                }
                                for (int px = 0; px <= max_px; px += rc.coarse) {
                                    if (!collides(plates[plate_idx], wpr, bw, bh,
                                                 rc.bm, rc.iwpr, rc.iw, rc.ih, px, test_py)) {
                                        found_px = px;
                                        found_py = test_py;
                                        goto coarse_hit;
                                    }
                                    // Also test boundary column if stride skips it
                                    if (px + rc.coarse > max_px && px < max_px) {
                                        if (!collides(plates[plate_idx], wpr, bw, bh,
                                                     rc.bm, rc.iwpr, rc.iw, rc.ih, max_px, test_py)) {
                                            found_px = max_px;
                                            found_py = test_py;
                                            goto coarse_hit;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    goto no_fit_this_rotation;
                    coarse_hit:

                    // Phase 2: refine within ±coarse of the coarse hit
                    {
                        int ry0 = std::max(0, found_py - rc.coarse);
                        int ry1 = std::min(bh - rc.ih, found_py + rc.coarse);
                        int rx0 = std::max(0, found_px - rc.coarse);
                        int rx1 = std::min(bw - rc.iw, found_px + rc.coarse);
                        for (int py = ry0; py <= ry1; ++py) {
                            for (int px = rx0; px <= rx1; ++px) {
                                if (!collides(plates[plate_idx], wpr, bw, bh,
                                             rc.bm, rc.iwpr, rc.iw, rc.ih, px, py)) {
                                    best_px = px;
                                    best_py = py;
                                    best_rot = rc.rot;
                                    best_plate = plate_idx;
                                    best_iw = rc.iw;
                                    best_ih = rc.ih;
                                    best_iwpr = rc.iwpr;
                                    best_bm = rc.bm;
                                    best_inflated_bb = rc.inflated_bb;
                                    placed = true;
                                    goto done_searching;
                                }
                            }
                        }
                    }
                    no_fit_this_rotation:;
                }
            }
            done_searching:

            if (placed) {
                if (best_plate > current_plate)
                    current_plate = best_plate;

                // Stamp the inflated shape to block future items
                stamp(plates[best_plate], wpr, bw, bh,
                      best_bm, best_iwpr, best_iw, best_ih,
                      best_px, best_py);

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
                item.rotation = best_rot;
                item.bed_idx = best_plate;
                item.itemid = item_sequence++;
                if (params.on_packed)
                    params.on_packed(item);
            } else {
                item.bed_idx = UNARRANGED;
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
                int x_start = std::max(0, (int)std::floor((xs[k] - ox) / res_mm));
                int x_end = std::min(iw - 1, (int)std::floor((xs[k + 1] - ox) / res_mm));

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
