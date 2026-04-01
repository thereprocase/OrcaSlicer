#pragma once

// BitmapNester: 2D bitmap-based concave shape nesting for OrcaSlicer
//
// Rasterizes ExPolygon silhouettes to bitmaps, uses word-parallel AND
// for collision detection, and skyline-based placement for tight packing.
// Supports multi-plate overflow and rotation.
//
// Zero external dependencies beyond what libslic3r already provides.

#include "ExPolygon.hpp"
#include "Arrange.hpp"
#include "BoundingBox.hpp"
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Slic3r { namespace arrangement {

class BitmapNester {
public:
    struct Config {
        double   resolution_mm  = 0.5;    // mm per pixel
        int      rotation_steps = 1;      // 1 = no rotation, 4 = 90° steps, etc.
        coord_t  spacing        = 0;      // min distance between items (scaled)
        double   bed_shrink_x   = 0;      // bed margin x (mm)
        double   bed_shrink_y   = 0;      // bed margin y (mm)
        int      max_plates     = 36;
        bool     allow_multi_materials_on_same_plate = true;
        std::function<void(unsigned, std::string)> progress;
        std::function<bool()> stopcondition;
    };

    static void arrange(ArrangePolygons &items,
                        const ArrangePolygons &excludes,
                        const BoundingBox &bed,
                        const Config &cfg)
    {
        if (items.empty()) return;

        // Effective bed after shrinkage
        BoundingBox ebed = bed;
        coord_t shrink_x = scaled(cfg.bed_shrink_x);
        coord_t shrink_y = scaled(cfg.bed_shrink_y);
        ebed.min.x() += shrink_x;
        ebed.min.y() += shrink_y;
        ebed.max.x() -= shrink_x;
        ebed.max.y() -= shrink_y;

        if (ebed.min.x() >= ebed.max.x() || ebed.min.y() >= ebed.max.y())
            return;

        coord_t bed_w = ebed.max.x() - ebed.min.x();
        coord_t bed_h = ebed.max.y() - ebed.min.y();
        double res = cfg.resolution_mm;
        int bw = std::max(1, (int)std::ceil(unscaled<double>(bed_w) / res));
        int bh = std::max(1, (int)std::ceil(unscaled<double>(bed_h) / res));

        // Cap bitmap size to prevent runaway allocation
        if ((size_t)bw * bh > 16'000'000) {
            double scale = std::sqrt(16'000'000.0 / ((size_t)bw * bh));
            bw = std::max(1, (int)(bw * scale));
            bh = std::max(1, (int)(bh * scale));
        }

        int words_per_row = (bw + 63) / 64;

        // Sort items largest-first by polygon area
        std::vector<size_t> order(items.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return std::abs(items[a].poly.area()) > std::abs(items[b].poly.area());
        });

        // Build rotation list
        std::vector<double> rotations;
        if (cfg.rotation_steps <= 1) {
            rotations.push_back(0.0);
        } else {
            for (int i = 0; i < cfg.rotation_steps; ++i)
                rotations.push_back(i * 2.0 * M_PI / cfg.rotation_steps);
        }

        int current_plate = 0;

        // Per-plate state: a forbidden bitmap tracking occupied pixels
        struct PlateState {
            std::vector<uint64_t> forbidden;
            int bw, bh, words_per_row;

            void init(int w, int h, int wpr) {
                bw = w; bh = h; words_per_row = wpr;
                forbidden.assign((size_t)wpr * h, 0);
            }

            bool collides(const std::vector<uint64_t> &item_bm,
                          int iw, int ih, int iwpr,
                          int px, int py) const
            {
                for (int y = 0; y < ih; ++y) {
                    int by = py + y;
                    if (by < 0 || by >= bh) continue;
                    for (int wx = 0; wx < iwpr; ++wx) {
                        uint64_t word = item_bm[(size_t)y * iwpr + wx];
                        if (word == 0) continue;

                        // Bit-shift to align item bitmap to bed position
                        int bit_offset = px + wx * 64;
                        if (bit_offset < 0 || bit_offset >= bw) continue;
                        int bed_word = bit_offset / 64;
                        int shift = bit_offset % 64;

                        if (bed_word >= 0 && bed_word < words_per_row) {
                            uint64_t shifted = word << shift;
                            if (forbidden[(size_t)by * words_per_row + bed_word] & shifted)
                                return true;
                        }
                        if (shift > 0 && bed_word + 1 >= 0 && bed_word + 1 < words_per_row) {
                            uint64_t shifted = word >> (64 - shift);
                            if (shifted && (forbidden[(size_t)by * words_per_row + bed_word + 1] & shifted))
                                return true;
                        }
                    }
                }
                return false;
            }

            void stamp(const std::vector<uint64_t> &item_bm,
                       int iw, int ih, int iwpr,
                       int px, int py)
            {
                for (int y = 0; y < ih; ++y) {
                    int by = py + y;
                    if (by < 0 || by >= bh) continue;
                    for (int wx = 0; wx < iwpr; ++wx) {
                        uint64_t word = item_bm[(size_t)y * iwpr + wx];
                        if (word == 0) continue;

                        int bit_offset = px + wx * 64;
                        if (bit_offset < 0 || bit_offset >= bw) continue;
                        int bed_word = bit_offset / 64;
                        int shift = bit_offset % 64;

                        if (bed_word >= 0 && bed_word < words_per_row) {
                            forbidden[(size_t)by * words_per_row + bed_word] |= (word << shift);
                        }
                        if (shift > 0 && bed_word + 1 >= 0 && bed_word + 1 < words_per_row) {
                            uint64_t hi = word >> (64 - shift);
                            if (hi)
                                forbidden[(size_t)by * words_per_row + bed_word + 1] |= hi;
                        }
                    }
                }
            }
        };

        // Stamp excludes onto plate 0
        std::vector<PlateState> plates(1);
        plates[0].init(bw, bh, words_per_row);

        for (auto &ex : excludes) {
            ExPolygon epoly = ex.poly;
            epoly.translate(ex.translation.x() - ebed.min.x(),
                           ex.translation.y() - ebed.min.y());
            if (ex.rotation != 0.0) epoly.rotate(ex.rotation);

            int ew, eh, ewpr;
            auto ebm = rasterize(epoly, res, bw, bh, ew, eh, ewpr);
            BoundingBox ebb = get_extents(epoly);
            int epx = (int)(unscaled<double>(ebb.min.x()) / res);
            int epy = (int)(unscaled<double>(ebb.min.y()) / res);
            plates[0].stamp(ebm, ew, eh, ewpr, epx, epy);
        }

        // Inflate spacing into a padding bitmap ring
        int pad_px = (cfg.spacing > 0) ? std::max(1, (int)std::ceil(unscaled<double>(cfg.spacing) / res)) : 0;

        // Place each item
        unsigned items_remaining = (unsigned)items.size();
        for (size_t oi = 0; oi < order.size(); ++oi) {
            size_t idx = order[oi];
            auto &item = items[idx];

            if (cfg.stopcondition && cfg.stopcondition()) break;

            if (cfg.progress)
                cfg.progress(items_remaining--, "(placing parts)");

            bool placed = false;
            double best_score = std::numeric_limits<double>::max();
            int best_px = 0, best_py = 0;
            double best_rot = 0.0;
            int best_plate = -1;
            int best_iw = 0, best_ih = 0, best_iwpr = 0;
            std::vector<uint64_t> best_bm;

            // Check allowed rotations from the item, or use our rotation list
            const auto &rots = (item.allowed_rotations.size() > 1 ||
                               (item.allowed_rotations.size() == 1 && item.allowed_rotations[0] != 0.0))
                               ? item.allowed_rotations : rotations;

            for (int plate_idx = 0; plate_idx <= current_plate && plate_idx < cfg.max_plates; ++plate_idx) {
                if ((int)plates.size() <= plate_idx) {
                    plates.emplace_back();
                    plates.back().init(bw, bh, words_per_row);
                }

                for (double rot : rots) {
                    // Rasterize the polygon at this rotation
                    ExPolygon rpoly = item.poly;
                    if (rot != 0.0) rpoly.rotate(rot);

                    // Inflate for spacing
                    ExPolygons inflated;
                    if (pad_px > 0) {
                        inflated = offset_ex(rpoly, scaled(pad_px * res));
                        if (!inflated.empty()) rpoly = inflated.front();
                    }

                    int iw, ih, iwpr;
                    auto ibm = rasterize(rpoly, res, bw, bh, iw, ih, iwpr);
                    if (ibm.empty()) continue;

                    // Skyline scan: bottom-left placement
                    // Scan Y from bottom, X from left
                    for (int py = 0; py <= bh - ih; ++py) {
                        for (int px = 0; px <= bw - iw; ++px) {
                            if (!plates[plate_idx].collides(ibm, iw, ih, iwpr, px, py)) {
                                // Score: prefer bottom-left (lower Y, then lower X)
                                double score = (double)py * 2.0 + (double)px * 0.001;
                                if (score < best_score) {
                                    best_score = score;
                                    best_px = px;
                                    best_py = py;
                                    best_rot = rot;
                                    best_plate = plate_idx;
                                    best_iw = iw;
                                    best_ih = ih;
                                    best_iwpr = iwpr;
                                    best_bm = ibm;
                                    placed = true;
                                    goto found_on_plate;
                                }
                            }
                        }
                        if (placed) break;
                    }
                found_on_plate:;
                }
                if (placed) break;
            }

            if (!placed && current_plate + 1 < cfg.max_plates) {
                // Try a new plate
                current_plate++;
                plates.emplace_back();
                plates.back().init(bw, bh, words_per_row);

                for (double rot : rots) {
                    ExPolygon rpoly = item.poly;
                    if (rot != 0.0) rpoly.rotate(rot);

                    ExPolygons inflated;
                    if (pad_px > 0) {
                        inflated = offset_ex(rpoly, scaled(pad_px * res));
                        if (!inflated.empty()) rpoly = inflated.front();
                    }

                    int iw, ih, iwpr;
                    auto ibm = rasterize(rpoly, res, bw, bh, iw, ih, iwpr);
                    if (ibm.empty()) continue;

                    for (int py = 0; py <= bh - ih; ++py) {
                        for (int px = 0; px <= bw - iw; ++px) {
                            if (!plates[current_plate].collides(ibm, iw, ih, iwpr, px, py)) {
                                best_px = px;
                                best_py = py;
                                best_rot = rot;
                                best_plate = current_plate;
                                best_iw = iw;
                                best_ih = ih;
                                best_iwpr = iwpr;
                                best_bm = ibm;
                                placed = true;
                                goto found_new_plate;
                            }
                        }
                        if (placed) break;
                    }
                found_new_plate:;
                    if (placed) break;
                }
            }

            if (placed) {
                // Stamp placed item (without inflation — stamp the actual shape)
                ExPolygon actual_poly = item.poly;
                if (best_rot != 0.0) actual_poly.rotate(best_rot);

                int aw, ah, awpr;
                auto abm = rasterize(actual_poly, res, bw, bh, aw, ah, awpr);

                // Compute actual stamp position (re-derive from placed position)
                BoundingBox rpbb = get_extents(actual_poly);
                BoundingBox inflated_poly = item.poly;
                if (best_rot != 0.0) {
                    ExPolygon tmp = item.poly;
                    tmp.rotate(best_rot);
                    ExPolygons infl = offset_ex(tmp, scaled(pad_px * res));
                    if (!infl.empty())
                        inflated_poly = get_extents(infl.front());
                    else
                        inflated_poly = get_extents(tmp);
                }

                // The rasterize function uses the polygon's own bounding box.
                // best_px/best_py are the top-left of the inflated raster.
                // The actual shape's offset within the inflated raster:
                int dx = 0, dy = 0; // actual vs inflated offset in pixels
                if (pad_px > 0) {
                    dx = pad_px; // inflation adds pad_px pixels on each side
                    dy = pad_px;
                }

                // Stamp with padding to maintain spacing
                plates[best_plate].stamp(best_bm, best_iw, best_ih, best_iwpr,
                                        best_px, best_py);

                // Convert pixel position back to the polygon's origin in absolute
                // bed coordinates. apply_arrange_result() interprets translation as
                // the absolute position of the polygon's local (0,0) point.
                //
                // The rasterizer places the inflated polygon's bbox.min at pixel
                // (best_px, best_py). The actual (uninflated, rotated) polygon's
                // bbox.min is offset inward from the inflated bbox by pad_px pixels.
                // The polygon's (0,0) is at -bbox.min relative to the bbox corner.
                //
                // So: origin_in_bed_mm = pixel_to_mm(best_px) + pad_offset_mm
                //                        - unscaled(rotated_bbox.min) + unscaled(ebed.min)

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
            } else {
                item.bed_idx = UNARRANGED;
            }
        }
    }

private:
    // Rasterize an ExPolygon into a bitmap.
    // Returns a vector of uint64 words, row-major, with bit 0 = leftmost pixel.
    // The bitmap covers the polygon's bounding box.
    // iw, ih, iwpr are output: width/height in pixels, words per row.
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

        // Scanline fill the contour
        scanline_fill(poly.contour, bb, res_mm, iw, ih, iwpr, bm, true);

        // Subtract holes
        for (auto &hole : poly.holes)
            scanline_fill(hole, bb, res_mm, iw, ih, iwpr, bm, false);

        return bm;
    }

    // Scanline rasterizer for a single polygon (contour or hole).
    // If 'set' is true, sets bits; if false, clears them.
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

            // Find X intersections with all edges
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

            // Fill between pairs
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
