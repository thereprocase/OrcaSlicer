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

namespace Slic3r { namespace arrangement {

class BitmapNester {
public:
    static void arrange(ArrangePolygons &items,
                        const ArrangePolygons &excludes,
                        const BoundingBox &bed,
                        const ArrangeParams &params)
    {
        if (items.empty()) return;

        const double res = 0.5; // mm per pixel

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
        int bw = std::max(1, (int)std::ceil(unscaled<double>(bed_w) / res));
        int bh = std::max(1, (int)std::ceil(unscaled<double>(bed_h) / res));

        // Cap bitmap size to prevent runaway allocation
        if ((size_t)bw * bh > 16'000'000) {
            double sc = std::sqrt(16'000'000.0 / ((size_t)bw * bh));
            bw = std::max(1, (int)(bw * sc));
            bh = std::max(1, (int)(bh * sc));
        }

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
            int px = (int)(unscaled<double>(pbb.min.x()) / res);
            int py = (int)(unscaled<double>(pbb.min.y()) / res);
            stamp(plates[plate_idx], wpr, bw, bh, bm, iwpr, iw, ih, px, py);
        };

        // Stamp excludes onto their respective plates
        for (auto &ex : excludes) {
            int plate = std::max(0, ex.bed_idx);
            ExPolygon epoly = ex.poly;
            // Apply the exclude's own transform
            if (ex.rotation != 0.0) epoly.rotate(ex.rotation);
            epoly.translate(ex.translation.x() - ebed.min.x(),
                           ex.translation.y() - ebed.min.y());
            // Inflate excludes by their per-item inflation
            if (ex.inflation > 0) {
                ExPolygons infl = offset_ex(epoly, ex.inflation);
                if (!infl.empty()) epoly = infl.front();
            }
            stamp_poly(plate, epoly);
        }

        // Also stamp params.excluded_regions (bed exclusion zones) onto plate 0
        for (auto &ex : params.excluded_regions) {
            ExPolygon epoly = ex.poly;
            if (ex.rotation != 0.0) epoly.rotate(ex.rotation);
            epoly.translate(ex.translation.x() - ebed.min.x(),
                           ex.translation.y() - ebed.min.y());
            stamp_poly(0, epoly);
        }

        int current_plate = 0;
        // Find highest plate index from excludes
        for (auto &ex : excludes)
            current_plate = std::max(current_plate, std::max(0, ex.bed_idx));

        int item_sequence = 0; // for sequential print ordering

        // Place each item
        unsigned items_remaining = (unsigned)items.size();
        for (size_t oi = 0; oi < order.size(); ++oi) {
            size_t idx = order[oi];
            auto &item = items[idx];

            if (params.stopcondition && params.stopcondition()) break;

            if (params.progressind)
                params.progressind(items_remaining--, "(placing parts)");

            // Per-item inflation (from brim, tree support, or user spacing)
            int pad_px = 0;
            if (item.inflation > 0)
                pad_px = std::max(1, (int)std::ceil(unscaled<double>(item.inflation) / res));
            else if (params.min_obj_distance > 0)
                pad_px = std::max(1, (int)std::ceil(unscaled<double>(params.min_obj_distance / 2) / res));

            // Use item's allowed rotations if specified, otherwise default list.
            // {0.0} means "rotation locked to 0" — respect it, don't override.
            const auto &rots = !item.allowed_rotations.empty()
                               ? item.allowed_rotations : default_rotations;

            bool placed = false;
            int best_px = 0, best_py = 0;
            double best_rot = 0.0;
            int best_plate = -1;
            int best_iw = 0, best_ih = 0, best_iwpr = 0;
            std::vector<uint64_t> best_bm;

            // Try each existing plate, then overflow to a new one
            int max_plate = std::min(current_plate + 1, MAX_PLATES - 1);
            for (int plate_idx = 0; plate_idx <= max_plate; ++plate_idx) {
                ensure_plate(plate_idx);

                for (double rot : rots) {
                    ExPolygon rpoly = item.poly;
                    if (rot != 0.0) rpoly.rotate(rot);

                    // Inflate for spacing
                    if (pad_px > 0) {
                        ExPolygons inflated = offset_ex(rpoly, scaled(pad_px * res));
                        if (!inflated.empty()) rpoly = inflated.front();
                    }

                    int iw, ih, iwpr_item;
                    auto ibm = rasterize(rpoly, res, bw, bh, iw, ih, iwpr_item);
                    if (ibm.empty()) continue;

                    // Bottom-left scan: first valid position wins
                    for (int py = 0; py <= bh - ih; ++py) {
                        for (int px = 0; px <= bw - iw; ++px) {
                            if (!collides(plates[plate_idx], wpr, bw, bh,
                                         ibm, iwpr_item, iw, ih, px, py)) {
                                best_px = px;
                                best_py = py;
                                best_rot = rot;
                                best_plate = plate_idx;
                                best_iw = iw;
                                best_ih = ih;
                                best_iwpr = iwpr_item;
                                best_bm = ibm;
                                placed = true;
                                goto done_searching;
                            }
                        }
                    }
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

                // Compute the polygon's origin (0,0) in absolute bed coords.
                // The rasterizer places the inflated polygon's bbox.min at pixel
                // (best_px, best_py). The actual polygon's bbox.min is pad_px
                // pixels inward. The origin is at -bbox.min from the bbox corner.
                ExPolygon rot_poly = item.poly;
                if (best_rot != 0.0) rot_poly.rotate(best_rot);
                BoundingBox rot_bb = get_extents(rot_poly);

                double pad_mm = pad_px * res;
                double origin_x = (best_px * res + pad_mm)
                                  - unscaled<double>(rot_bb.min.x())
                                  + unscaled<double>(ebed.min.x());
                double origin_y = (best_py * res + pad_mm)
                                  - unscaled<double>(rot_bb.min.y())
                                  + unscaled<double>(ebed.min.y());

                item.translation = Vec2crd{scaled(origin_x), scaled(origin_y)};
                item.rotation = best_rot;
                item.bed_idx = best_plate;
                item.itemid = item_sequence++;
            } else {
                item.bed_idx = UNARRANGED;
            }
        }
    }

private:
    // Collision check: does item bitmap at (px, py) overlap the plate bitmap?
    static bool collides(const std::vector<uint64_t> &plate,
                         int wpr, int bw, int bh,
                         const std::vector<uint64_t> &item_bm,
                         int iwpr, int iw, int ih,
                         int px, int py)
    {
        for (int y = 0; y < ih; ++y) {
            int by = py + y;
            if (by < 0 || by >= bh) continue;
            for (int wx = 0; wx < iwpr; ++wx) {
                uint64_t word = item_bm[(size_t)y * iwpr + wx];
                if (word == 0) continue;

                int bit_offset = px + wx * 64;
                if (bit_offset + 63 < 0 || bit_offset >= bw) continue;
                int bed_word = (bit_offset >= 0) ? bit_offset / 64 : (bit_offset - 63) / 64;
                int shift = bit_offset - bed_word * 64;
                if (shift < 0 || shift >= 64) continue;

                if (bed_word >= 0 && bed_word < wpr) {
                    uint64_t shifted = word << shift;
                    if (plate[(size_t)by * wpr + bed_word] & shifted)
                        return true;
                }
                if (shift > 0 && bed_word + 1 >= 0 && bed_word + 1 < wpr) {
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
        for (int y = 0; y < ih; ++y) {
            int by = py + y;
            if (by < 0 || by >= bh) continue;
            for (int wx = 0; wx < iwpr; ++wx) {
                uint64_t word = item_bm[(size_t)y * iwpr + wx];
                if (word == 0) continue;

                int bit_offset = px + wx * 64;
                if (bit_offset + 63 < 0 || bit_offset >= bw) continue;
                int bed_word = (bit_offset >= 0) ? bit_offset / 64 : (bit_offset - 63) / 64;
                int shift = bit_offset - bed_word * 64;
                if (shift < 0 || shift >= 64) continue;

                if (bed_word >= 0 && bed_word < wpr) {
                    plate[(size_t)by * wpr + bed_word] |= (word << shift);
                }
                if (shift > 0 && bed_word + 1 >= 0 && bed_word + 1 < wpr) {
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

        for (int y = 0; y < ih; ++y) {
            double scan_y = oy + (y + 0.5) * res_mm;

            std::vector<double> xs;
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

            std::sort(xs.begin(), xs.end());

            for (size_t k = 0; k + 1 < xs.size(); k += 2) {
                int x_start = std::max(0, (int)std::floor((xs[k] - ox) / res_mm));
                int x_end = std::min(iw - 1, (int)std::floor((xs[k + 1] - ox) / res_mm));

                for (int x = x_start; x <= x_end; ++x) {
                    int word = x / 64;
                    int bit = x % 64;
                    if (set)
                        bm[(size_t)y * iwpr + word] |= (uint64_t(1) << bit);
                    else
                        bm[(size_t)y * iwpr + word] &= ~(uint64_t(1) << bit);
                }
            }
        }
    }
};

}} // namespace Slic3r::arrangement
