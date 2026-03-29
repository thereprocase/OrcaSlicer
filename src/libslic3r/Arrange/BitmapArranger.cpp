#include "BitmapArranger.hpp"
#include "libslic3r/ClipperUtils.hpp"

#include <boost/log/trivial.hpp>

// Coordinate systems used in this module:
//   Scaled coords (coord_t): Slic3r internal units. 1mm = scaled(1.0).
//   Pixel coords (int): px = (scaled_coord - bed_origin) / resolution.
//     Pixel (0,0) = effective bed min corner (after purge pad shrink).
//   Placement: ArrangePolygon.translation = bed_min + (px * res) - item_bb_min
//     Applied AFTER ArrangePolygon.rotation.

// 64M words * 8 bytes = 512MB max bitmap allocation per plate
constexpr size_t MAX_BITMAP_WORDS = 64'000'000u;

#ifdef _MSC_VER
#include <intrin.h>
// Count trailing zeros. Precondition: x != 0 (callers must guard).
static inline int ctz64(uint64_t x) {
    unsigned long idx;
    if (!_BitScanForward64(&idx, x)) return 0;
    return (int)idx;
}
#else
// Count trailing zeros. Precondition: x != 0 (callers must guard).
static inline int ctz64(uint64_t x) { return x ? __builtin_ctzll(x) : 0; }
#endif

namespace Slic3r { namespace arrangement {

BitmapArranger::BitmapItem BitmapArranger::rasterize(
    const ExPolygon &poly, coord_t inflation, coord_t res)
{
    BitmapItem item;

    // Inflate the polygon for spacing
    ExPolygons inflated;
    if (inflation > 0)
        inflated = offset_ex(poly, inflation);
    else
        inflated = {poly};

    if (inflated.empty())
        return item;

    const ExPolygon &ep = inflated.front();
    BoundingBox bb = get_extents(ep);

    item.offset_x = bb.min.x();
    item.offset_y = bb.min.y();
    // +1 prevents the rightmost/topmost edge of the polygon from being clipped
    item.width_px  = (int)std::ceil((double)(bb.max.x() - bb.min.x()) / res) + 1;
    item.height_px = (int)std::ceil((double)(bb.max.y() - bb.min.y()) / res) + 1;

    // Degenerate polygon after inflation can produce negative dimensions
    if (item.width_px <= 0 || item.height_px <= 0) {
        item.width_px = item.height_px = item.width_words = 0;
        return item;
    }

    item.width_words = (item.width_px + 63) / 64;

    item.bits.assign(item.width_words * item.height_px, 0);

    // Scanline rasterization
    std::vector<coord_t> x_intersections;
    for (int py = 0; py < item.height_px; py++) {
        coord_t y = bb.min.y() + (coord_t)(py * res + res / 2);

        // Collect x-intersections with all edges (contour + holes)
        x_intersections.clear();

        auto scan_ring = [&](const Polygon &ring) {
            size_t n = ring.points.size();
            for (size_t i = 0; i < n; i++) {
                const Point &p0 = ring.points[i];
                const Point &p1 = ring.points[(i + 1) % n];
                if ((p0.y() <= y && p1.y() > y) || (p1.y() <= y && p0.y() > y)) {
                    double t = (double)(y - p0.y()) / (double)(p1.y() - p0.y());
                    coord_t x = p0.x() + (coord_t)(t * (p1.x() - p0.x()));
                    x_intersections.push_back(x);
                }
            }
        };

        scan_ring(ep.contour);
        for (const Polygon &hole : ep.holes)
            scan_ring(hole);

        std::sort(x_intersections.begin(), x_intersections.end());

        // Fill between pairs of intersections (even-odd rule)
        for (size_t i = 0; i + 1 < x_intersections.size(); i += 2) {
            int px_start = std::max(0, (int)((x_intersections[i] - bb.min.x()) / res));
            int px_end   = std::min(item.width_px - 1,
                                    (int)((x_intersections[i + 1] - bb.min.x()) / res));
            int sw = px_start / 64, ew = px_end / 64;
            int sb = px_start % 64, eb = px_end % 64;
            if (sw == ew) {
                uint64_t mask = ((uint64_t(2) << eb) - 1) & ~((uint64_t(1) << sb) - 1);
                item.bits[py * item.width_words + sw] |= mask;
            } else {
                item.bits[py * item.width_words + sw] |= ~((uint64_t(1) << sb) - 1);
                for (int w = sw + 1; w < ew; w++)
                    item.bits[py * item.width_words + w] = ~uint64_t(0);
                item.bits[py * item.width_words + ew] |= (uint64_t(2) << eb) - 1;
            }
        }
    }

    return item;
}

BitmapArranger::BitmapItem BitmapArranger::rasterize_triangles(
    const Points &tri_verts, coord_t inflation, coord_t res)
{
    BitmapItem item;
    if (tri_verts.size() < 3) return item;

    // Bounding box of all triangle vertices, expanded by inflation
    BoundingBox bb(tri_verts);
    bb.min -= Vec2crd(inflation, inflation);
    bb.max += Vec2crd(inflation, inflation);

    item.offset_x = bb.min.x();
    item.offset_y = bb.min.y();
    item.width_px  = (int)std::ceil((double)(bb.max.x() - bb.min.x()) / res) + 1;
    item.height_px = (int)std::ceil((double)(bb.max.y() - bb.min.y()) / res) + 1;

    if (item.width_px <= 0 || item.height_px <= 0) {
        item.width_px = item.height_px = item.width_words = 0;
        return item;
    }

    item.width_words = (item.width_px + 63) / 64;
    item.bits.assign(item.width_words * item.height_px, 0);

    auto set_pixel = [&](int px, int py) {
        if (px < 0 || px >= item.width_px || py < 0 || py >= item.height_px) return;
        item.bits[(size_t)py * item.width_words + px / 64] |= (uint64_t(1) << (px % 64));
    };

    double inv_res = 1.0 / res;

    // Scanline-fill each triangle directly into the bitmap.
    // Same algorithm GPUs use — for each scanline row that intersects the triangle,
    // find the left and right edges and fill between them.
    for (size_t ti = 0; ti + 2 < tri_verts.size(); ti += 3) {
        // Convert triangle vertices to pixel coordinates
        double x0 = (tri_verts[ti].x()     - bb.min.x()) * inv_res;
        double y0 = (tri_verts[ti].y()     - bb.min.y()) * inv_res;
        double x1 = (tri_verts[ti + 1].x() - bb.min.x()) * inv_res;
        double y1 = (tri_verts[ti + 1].y() - bb.min.y()) * inv_res;
        double x2 = (tri_verts[ti + 2].x() - bb.min.x()) * inv_res;
        double y2 = (tri_verts[ti + 2].y() - bb.min.y()) * inv_res;

        // Sort vertices by Y
        if (y0 > y1) { std::swap(x0, x1); std::swap(y0, y1); }
        if (y0 > y2) { std::swap(x0, x2); std::swap(y0, y2); }
        if (y1 > y2) { std::swap(x1, x2); std::swap(y1, y2); }

        int py_min = std::max(0, (int)y0);
        int py_max = std::min(item.height_px - 1, (int)y2);

        for (int py = py_min; py <= py_max; py++) {
            double row = py + 0.5;
            double left = 1e18, right = -1e18;

            // Intersect scanline with each of the 3 edges
            auto edge_intersect = [&](double ax, double ay, double bx, double by) {
                if ((ay <= row && by > row) || (by <= row && ay > row)) {
                    double t = (row - ay) / (by - ay);
                    double ix = ax + t * (bx - ax);
                    left = std::min(left, ix);
                    right = std::max(right, ix);
                }
            };

            edge_intersect(x0, y0, x1, y1);
            edge_intersect(x1, y1, x2, y2);
            edge_intersect(x0, y0, x2, y2);

            if (left > right) continue;

            int px_left  = std::max(0, (int)left);
            int px_right = std::min(item.width_px - 1, (int)right);

            // Fill the span using word-level operations for speed
            if (px_left <= px_right) {
                int sw = px_left / 64, ew = px_right / 64;
                int sb = px_left % 64, eb = px_right % 64;
                if (sw == ew) {
                    uint64_t mask = ((uint64_t(2) << eb) - 1) & ~((uint64_t(1) << sb) - 1);
                    item.bits[(size_t)py * item.width_words + sw] |= mask;
                } else {
                    item.bits[(size_t)py * item.width_words + sw] |= ~((uint64_t(1) << sb) - 1);
                    for (int w = sw + 1; w < ew; w++)
                        item.bits[(size_t)py * item.width_words + w] = ~uint64_t(0);
                    item.bits[(size_t)py * item.width_words + ew] |= (uint64_t(2) << eb) - 1;
                }
            }
        }
    }

    // If inflation > 0, dilate the bitmap by inflate_px pixels.
    // Simple dilation: for each set bit, set its neighbors.
    if (inflation > 0) {
        int inflate_px = std::max(1, (int)(inflation * inv_res));
        std::vector<uint64_t> dilated = item.bits; // copy
        for (int py = 0; py < item.height_px; py++) {
            for (int wx = 0; wx < item.width_words; wx++) {
                uint64_t word = item.bits[(size_t)py * item.width_words + wx];
                if (word == 0) continue;
                // Expand each set bit in all directions
                while (word) {
                    int bit = ctz64(word);
                    int px = wx * 64 + bit;
                    for (int dy = -inflate_px; dy <= inflate_px; dy++) {
                        int ny = py + dy;
                        if (ny < 0 || ny >= item.height_px) continue;
                        for (int dx = -inflate_px; dx <= inflate_px; dx++) {
                            int nx = px + dx;
                            if (nx < 0 || nx >= item.width_px) continue;
                            dilated[(size_t)ny * item.width_words + nx / 64] |= (uint64_t(1) << (nx % 64));
                        }
                    }
                    word &= word - 1;
                }
            }
        }
        item.bits = std::move(dilated);
    }

    return item;
}

BitmapArranger::BitmapItem BitmapArranger::rotate_bitmap(
    const BitmapItem &src, double angle_rad, coord_t res)
{
    if (src.width_px <= 0 || src.height_px <= 0)
        return src;

    double cos_a = std::cos(angle_rad), sin_a = std::sin(angle_rad);

    // Compute the bounding box of the rotated source bitmap.
    // The four corners of the source, rotated, give the new extents.
    double hw = src.width_px * 0.5, hh = src.height_px * 0.5;
    double corners_x[4] = { -hw, hw, hw, -hw };
    double corners_y[4] = { -hh, -hh, hh, hh };
    double min_x = 1e18, max_x = -1e18, min_y = 1e18, max_y = -1e18;
    for (int i = 0; i < 4; i++) {
        double rx = corners_x[i] * cos_a - corners_y[i] * sin_a;
        double ry = corners_x[i] * sin_a + corners_y[i] * cos_a;
        min_x = std::min(min_x, rx); max_x = std::max(max_x, rx);
        min_y = std::min(min_y, ry); max_y = std::max(max_y, ry);
    }

    BitmapItem dst;
    dst.width_px = (int)std::ceil(max_x - min_x) + 1;
    dst.height_px = (int)std::ceil(max_y - min_y) + 1;
    dst.width_words = (dst.width_px + 63) / 64;
    dst.bits.assign(dst.width_words * dst.height_px, 0);

    // New offset: rotate the source center and compute the new BB min in scaled coords.
    double src_cx = src.offset_x + hw * res;
    double src_cy = src.offset_y + hh * res;
    // The center stays the same after rotation (we rotate around the object's own center)
    dst.offset_x = (coord_t)(src_cx - (dst.width_px * 0.5) * res);
    dst.offset_y = (coord_t)(src_cy - (dst.height_px * 0.5) * res);

    // Reverse mapping: for each destination pixel, find the source pixel.
    // Rotate the destination coordinate backwards by -angle to find the source.
    double dst_cx = dst.width_px * 0.5;
    double dst_cy = dst.height_px * 0.5;

    for (int dy = 0; dy < dst.height_px; dy++) {
        double fy = dy - dst_cy;
        for (int dx = 0; dx < dst.width_px; dx++) {
            double fx = dx - dst_cx;
            // Inverse rotation
            int sx = (int)(fx * cos_a + fy * sin_a + hw);
            int sy = (int)(-fx * sin_a + fy * cos_a + hh);
            if (sx >= 0 && sx < src.width_px && sy >= 0 && sy < src.height_px) {
                if (src.bits[(size_t)sy * src.width_words + sx / 64] & (uint64_t(1) << (sx % 64))) {
                    dst.bits[(size_t)dy * dst.width_words + dx / 64] |= (uint64_t(1) << (dx % 64));
                }
            }
        }
    }

    return dst;
}

// ============================================================
// 3D-aware nesting: Z-slice collision
// ============================================================

BitmapArranger::SliceStack BitmapArranger::rasterize_slices(
    const Points &tri_verts, const std::vector<float> &tri_z,
    coord_t inflation, coord_t res,
    float slice_height_mm, float z_clearance_mm)
{
    SliceStack stack;
    if (tri_verts.size() < 3 || tri_z.size() != tri_verts.size()) return stack;

    // Find max Z to determine number of slices
    float max_z = 0;
    for (float z : tri_z) max_z = std::max(max_z, z);
    stack.n_slices = std::max(1, (int)std::ceil(max_z / slice_height_mm));

    // Bin triangles into slices. Each triangle goes to every slice its Z range
    // intersects, padded by z_clearance.
    std::vector<Points> slice_tris(stack.n_slices);
    for (size_t ti = 0; ti + 2 < tri_verts.size(); ti += 3) {
        float z_min = std::min({tri_z[ti], tri_z[ti + 1], tri_z[ti + 2]});
        float z_max = std::max({tri_z[ti], tri_z[ti + 1], tri_z[ti + 2]});
        int s_lo = std::max(0, (int)((z_min - z_clearance_mm) / slice_height_mm));
        int s_hi = std::min(stack.n_slices - 1, (int)((z_max + z_clearance_mm) / slice_height_mm));
        for (int s = s_lo; s <= s_hi; s++) {
            slice_tris[s].push_back(tri_verts[ti]);
            slice_tris[s].push_back(tri_verts[ti + 1]);
            slice_tris[s].push_back(tri_verts[ti + 2]);
        }
    }

    // Rasterize each slice using the existing triangle rasterizer
    stack.slices.resize(stack.n_slices);
    for (int s = 0; s < stack.n_slices; s++) {
        if (!slice_tris[s].empty())
            stack.slices[s] = rasterize_triangles(slice_tris[s], inflation, res);
    }

    return stack;
}

bool BitmapArranger::collides_3d(
    const SliceStack &bed_stack, int bed_w, int bed_w_px, int bed_h,
    const SliceStack &item_stack, int ox, int oy)
{
    for (int z = 0; z < item_stack.n_slices; z++) {
        if (item_stack.slices[z].width_px <= 0) continue;
        if (z < bed_stack.n_slices && bed_stack.slices[z].width_px > 0) {
            // Check against occupied pixels on this bed slice
            if (collides(bed_stack.slices[z].bits, bed_w, bed_w_px, bed_h,
                         item_stack.slices[z], ox, oy))
                return true;
        }
        // Upper slices beyond bed stack: bed edges are still enforced by
        // collides() bounds checking (out-of-bitmap = collision), so items
        // can't extend past the build volume. No explicit check needed here
        // because find_placement constrains ox/oy to valid bitmap positions.
    }
    return false;
}

void BitmapArranger::stamp_3d(
    SliceStack &bed_stack, int bed_w, int bed_w_px, int bed_h,
    const SliceStack &item_stack, int ox, int oy,
    const std::function<void(std::vector<uint64_t>&)> &init_slice)
{
    // Extend bed stack if item is taller
    while (bed_stack.n_slices < item_stack.n_slices) {
        BitmapItem new_slice;
        new_slice.bits.assign(bed_w * bed_h, 0);
        new_slice.width_words = bed_w;
        new_slice.width_px = bed_w_px;
        new_slice.height_px = bed_h;
        // Stamp physical obstacles (cut notch, fixed items) on new slices.
        // Bed edges are enforced by bitmap dimensions in collides().
        if (init_slice) init_slice(new_slice.bits);
        bed_stack.slices.push_back(std::move(new_slice));
        bed_stack.n_slices++;
    }
    for (int z = 0; z < item_stack.n_slices; z++) {
        if (item_stack.slices[z].width_px <= 0) continue;
        if (bed_stack.slices[z].bits.empty()) {
            bed_stack.slices[z].bits.assign(bed_w * bed_h, 0);
            bed_stack.slices[z].width_words = bed_w;
            bed_stack.slices[z].width_px = bed_w_px;
            bed_stack.slices[z].height_px = bed_h;
            if (init_slice) init_slice(bed_stack.slices[z].bits);
        }
        stamp(bed_stack.slices[z].bits, bed_stack.slices[z].width_words,
              bed_w_px, bed_h, item_stack.slices[z], ox, oy);
    }
}

std::optional<std::pair<int,int>> BitmapArranger::find_placement_3d(
    const SliceStack &bed_stack, int bed_w, int bed_w_px, int bed_h,
    const SliceStack &item_stack, int step)
{
    // Center-out scan, same as 2D but using 3D collision
    // Use the largest slice's dimensions for bounds
    int item_w = 0, item_h = 0;
    for (const auto &s : item_stack.slices) {
        item_w = std::max(item_w, s.width_px);
        item_h = std::max(item_h, s.height_px);
    }
    int max_x = bed_w_px - item_w;
    int max_y = bed_h - item_h;
    if (max_x < 0 || max_y < 0) return std::nullopt;

    int cx = max_x / 2, cy = max_y / 2;
    int max_radius = std::max(max_x, max_y);

    for (int r = 0; r <= max_radius; r += step) {
        int y_lo = std::max(0, cy - r), y_hi = std::min(max_y, cy + r);
        int x_lo = std::max(0, cx - r), x_hi = std::min(max_x, cx + r);

        for (int y = y_lo; y <= y_hi; y += step) {
            for (int x = x_lo; x <= x_hi; x += step) {
                if (r > step && x > x_lo + step && x < x_hi - step &&
                    y > y_lo + step && y < y_hi - step)
                    continue;

                if (!collides_3d(bed_stack, bed_w, bed_w_px, bed_h, item_stack, x, y)) {
                    // Refine for tightest center-ward placement
                    int best_x = x, best_y = y;
                    int best_dist = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                    for (int ry = std::max(0, y - step + 1); ry <= std::min(max_y, y + step - 1); ry++) {
                        for (int rx = std::max(0, x - step + 1); rx <= std::min(max_x, x + step - 1); rx++) {
                            int d = (rx - cx) * (rx - cx) + (ry - cy) * (ry - cy);
                            if (d < best_dist && !collides_3d(bed_stack, bed_w, bed_w_px, bed_h, item_stack, rx, ry)) {
                                best_x = rx; best_y = ry; best_dist = d;
                            }
                        }
                    }
                    return std::make_pair(best_x, best_y);
                }
            }
        }
    }
    return std::nullopt;
}

BitmapArranger::SliceStack BitmapArranger::rotate_stack(
    const SliceStack &src, double angle_rad, coord_t res)
{
    SliceStack dst;
    dst.n_slices = src.n_slices;
    dst.slices.reserve(src.n_slices);
    for (const auto &slice : src.slices) {
        if (slice.width_px > 0)
            dst.slices.push_back(rotate_bitmap(slice, angle_rad, res));
        else
            dst.slices.push_back({});
    }
    return dst;
}

void BitmapArranger::stamp(
    std::vector<uint64_t> &bed_bits, int bed_w, int bed_w_px, int bed_h,
    const BitmapItem &item, int ox, int oy)
{
    for (int py = 0; py < item.height_px; py++) {
        int by = oy + py;
        if (by < 0 || by >= bed_h) continue;
        for (int wx = 0; wx < item.width_words; wx++) {
            uint64_t iword = item.bits[py * item.width_words + wx];
            if (iword == 0) continue;
            int base_px = ox + wx * 64;
            if (base_px + 63 < 0) continue;
            if (base_px >= bed_w_px) break;

            if (base_px < 0) {
                // Bit-by-bit for left-edge partial overlap
                while (iword) {
                    int bit = ctz64(iword);
                    int bx = base_px + bit;
                    if (bx >= 0 && bx < bed_w_px) {
                        bed_bits[(size_t)by * bed_w + bx / 64] |= (uint64_t(1) << (bx % 64));
                    }
                    iword &= iword - 1;
                }
                continue;
            }

            int bed_word_idx = base_px / 64;
            int shift = base_px % 64;
            if (shift == 0) {
                if (bed_word_idx < bed_w)
                    bed_bits[(size_t)by * bed_w + bed_word_idx] |= iword;
            } else {
                // Item pixels span two bed words when not 64-bit aligned.
                // lo = bits landing in bed_word[idx], hi = overflow into bed_word[idx+1].
                uint64_t lo = iword << shift;
                uint64_t hi = (shift < 64) ? (iword >> (64 - shift)) : 0;
                if (bed_word_idx < bed_w)
                    bed_bits[(size_t)by * bed_w + bed_word_idx] |= lo;
                if (bed_word_idx + 1 < bed_w)
                    bed_bits[(size_t)by * bed_w + bed_word_idx + 1] |= hi;
            }
        }
    }
}

bool BitmapArranger::collides(
    const std::vector<uint64_t> &bed_bits, int bed_w, int bed_w_px, int bed_h,
    const BitmapItem &item, int ox, int oy)
{
    for (int py = 0; py < item.height_px; py++) {
        int by = oy + py;
        if (by < 0 || by >= bed_h) return true;

        for (int wx = 0; wx < item.width_words; wx++) {
            uint64_t iword = item.bits[py * item.width_words + wx];
            if (iword == 0) continue;

            int base_px = ox + wx * 64;

            if (base_px + 63 < 0) continue;
            if (base_px >= bed_w_px) return true;

            // Negative base_px: bit-by-bit fallback (shift by negative is UB)
            if (base_px < 0) {
                for (int b = 0; b < 64; b++) {
                    if (!(iword & (uint64_t(1) << b))) continue;
                    int bx = base_px + b;
                    if (bx < 0) continue;
                    if (bx >= bed_w_px) return true;
                    int bword = bx / 64;
                    int bbit  = bx % 64;
                    if (bed_bits[(size_t)by * bed_w + bword] & (uint64_t(1) << bbit))
                        return true;
                }
                continue;
            }

            int bed_word_idx = base_px / 64;
            int shift = base_px % 64;

            if (shift == 0) {
                if (bed_word_idx < bed_w) {
                    if (bed_bits[(size_t)by * bed_w + bed_word_idx] & iword)
                        return true;
                }
            } else {
                // Item pixels span two bed words when not 64-bit aligned.
                // lo = bits landing in bed_word[idx], hi = overflow into bed_word[idx+1].
                uint64_t lo = iword << shift;
                uint64_t hi = (shift < 64) ? (iword >> (64 - shift)) : 0;

                if (bed_word_idx < bed_w) {
                    if (bed_bits[(size_t)by * bed_w + bed_word_idx] & lo)
                        return true;
                }
                if (bed_word_idx + 1 < bed_w) {
                    if (hi && (bed_bits[(size_t)by * bed_w + bed_word_idx + 1] & hi))
                        return true;
                }
            }
        }
    }
    return false;
}

// ============================================================
// Skyline placement — O(bed_width) per item
// ============================================================

BitmapArranger::ItemProfile BitmapArranger::compute_profile(
    const BitmapItem &item, int bed_h)
{
    ItemProfile prof;
    prof.bw = item.width_px;
    prof.bh = item.height_px;
    prof.max_y = bed_h - item.height_px;
    if (prof.max_y < 0) return prof;

    // Bottom profile: lowest set pixel per column. Sentinel (bed_h+1) for inactive columns.
    prof.bottom.assign(prof.bw, bed_h + 1);
    for (int x = 0; x < prof.bw; x++) {
        for (int y = 0; y < prof.bh; y++) {
            int word = x / 64;
            int bit = x % 64;
            if (item.bits[(size_t)y * item.width_words + word] & (uint64_t(1) << bit)) {
                prof.bottom[x] = y;
                break;
            }
        }
    }

    // Top profile: highest set pixel + 1 per active column (for skyline update)
    for (int x = 0; x < prof.bw; x++) {
        if (prof.bottom[x] > bed_h) continue; // inactive column
        for (int y = prof.bh - 1; y >= 0; y--) {
            int word = x / 64;
            int bit = x % 64;
            if (item.bits[(size_t)y * item.width_words + word] & (uint64_t(1) << bit)) {
                prof.top_pairs.push_back({x, y + 1});
                break;
            }
        }
    }

    return prof;
}

std::optional<std::pair<int,int>> BitmapArranger::find_placement_skyline(
    const std::vector<int> &skyline, int bed_w_px, int bed_h,
    const ItemProfile &profile)
{
    if (profile.max_y < 0 || profile.bw > bed_w_px)
        return std::nullopt;

    int n_pos = bed_w_px - profile.bw + 1;
    if (n_pos <= 0) return std::nullopt;

    constexpr int STRIDE = 8;
    constexpr int REFINE = 10;

    // Coarse pass: evaluate every STRIDE-th position
    int best_x = -1, best_y = INT_MAX;

    for (int x = 0; x < n_pos; x += STRIDE) {
        // Compute placement Y: max(skyline[x+col] - bottom[col]) for all columns
        int y = 0;
        for (int c = 0; c < profile.bw; c++) {
            int diff = skyline[x + c] - profile.bottom[c];
            if (diff > y) y = diff;
        }
        if (y < 0) y = 0;
        if (y < best_y) {
            best_y = y;
            best_x = x;
            if (y == 0) break; // floor placement — can't do better
        }
    }

    if (best_x < 0 || best_y > profile.max_y)
        return std::nullopt;

    // Refine around coarse winner
    if (best_y > 0) {
        int ref_lo = std::max(0, best_x - REFINE);
        int ref_hi = std::min(n_pos, best_x + REFINE + 1);
        for (int x = ref_lo; x < ref_hi; x++) {
            int y = 0;
            for (int c = 0; c < profile.bw; c++) {
                int diff = skyline[x + c] - profile.bottom[c];
                if (diff > y) y = diff;
            }
            if (y < 0) y = 0;
            if (y < best_y) {
                best_y = y;
                best_x = x;
            }
        }
    }

    if (best_y > profile.max_y)
        return std::nullopt;

    return std::make_pair(best_x, best_y);
}

std::optional<std::pair<int,int>> BitmapArranger::find_placement(
    const std::vector<uint64_t> &bed_bits, int bed_w, int bed_w_px, int bed_h,
    const BitmapItem &item, int step)
{
    // Multi-tier search: coarse → medium → fine, center-out at each tier.
    // Like Google Maps: continent → city → street.
    int max_x = bed_w_px - item.width_px;
    int max_y = bed_h - item.height_px;
    if (max_x < 0 || max_y < 0) return std::nullopt;

    int cx = max_x / 2;
    int cy = max_y / 2;

    // Tier 1: Very coarse scan (4x step) to find candidate region
    int coarse_step = step * 4;
    int hit_x = -1, hit_y = -1;
    int best_dist = INT_MAX;

    // Center-out spiral at coarse resolution
    int max_radius = std::max(max_x, max_y);
    for (int r = 0; r <= max_radius; r += coarse_step) {
        int y_lo = std::max(0, cy - r);
        int y_hi = std::min(max_y, cy + r);
        int x_lo = std::max(0, cx - r);
        int x_hi = std::min(max_x, cx + r);

        for (int y = y_lo; y <= y_hi; y += coarse_step) {
            for (int x = x_lo; x <= x_hi; x += coarse_step) {
                if (r > coarse_step && x > x_lo + coarse_step && x < x_hi - coarse_step &&
                    y > y_lo + coarse_step && y < y_hi - coarse_step)
                    continue;

                if (!collides(bed_bits, bed_w, bed_w_px, bed_h, item, x, y)) {
                    int d = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                    if (d < best_dist) {
                        hit_x = x; hit_y = y; best_dist = d;
                    }
                    goto tier2; // found a region, refine it
                }
            }
        }
    }
    // Coarse scan found nothing — item doesn't fit
    return std::nullopt;

tier2:
    // Tier 2: Medium scan around the coarse hit (within coarse_step radius)
    {
        int search_r = coarse_step;
        int y_lo = std::max(0, hit_y - search_r);
        int y_hi = std::min(max_y, hit_y + search_r);
        int x_lo = std::max(0, hit_x - search_r);
        int x_hi = std::min(max_x, hit_x + search_r);

        best_dist = INT_MAX;
        for (int y = y_lo; y <= y_hi; y += step) {
            for (int x = x_lo; x <= x_hi; x += step) {
                if (!collides(bed_bits, bed_w, bed_w_px, bed_h, item, x, y)) {
                    int d = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                    if (d < best_dist) {
                        hit_x = x; hit_y = y; best_dist = d;
                    }
                }
            }
        }
    }

    // Tier 3: Pixel-level refinement around the medium hit (within one step)
    {
        int best_x = hit_x, best_y = hit_y;
        best_dist = (hit_x - cx) * (hit_x - cx) + (hit_y - cy) * (hit_y - cy);
        for (int ry = std::max(0, hit_y - step); ry <= std::min(max_y, hit_y + step); ry++) {
            for (int rx = std::max(0, hit_x - step); rx <= std::min(max_x, hit_x + step); rx++) {
                int d = (rx - cx) * (rx - cx) + (ry - cy) * (ry - cy);
                if (d < best_dist && !collides(bed_bits, bed_w, bed_w_px, bed_h, item, rx, ry)) {
                    best_x = rx; best_y = ry; best_dist = d;
                }
            }
        }
        return std::make_pair(best_x, best_y);
    }
}

void BitmapArranger::arrange(
    ArrangePolygons &arrangables,
    const ArrangePolygons &excludes,
    const BoundingBox &bed,
    const ArrangeParams &params)
{
    coord_t res = scaled<coord_t>(std::clamp(params.bitmap_resolution_mm, 0.1f, 2.0f));
    coord_t inflation = params.min_obj_distance / 2;

    // Apply purge pad: shrink the effective bed along one edge
    BoundingBox effective_bed = bed;
    if (params.avoid_purge_pad && params.purge_pad_mm > 0.f) {
        coord_t pad = scaled(params.purge_pad_mm);
        switch (params.purge_pad_edge) {
        case 0: effective_bed.min.y() += pad; break;  // front (Y min)
        case 1: effective_bed.max.y() -= pad; break;  // back (Y max)
        case 2: effective_bed.min.x() += pad; break;  // left (X min)
        case 3: effective_bed.max.x() -= pad; break;  // right (X max)
        default: break;
        }
    }

    // Bed dimensions in pixels. +1 prevents edge pixels from clipping:
    // a polygon whose extent equals the bed size must still fit entirely.
    int bed_w_px = (int)std::ceil((double)(effective_bed.max.x() - effective_bed.min.x()) / res) + 1;
    int bed_h_px = (int)std::ceil((double)(effective_bed.max.y() - effective_bed.min.y()) / res) + 1;
    int bed_w_words = (bed_w_px + 63) / 64;

    // Cap bitmap size to prevent unbounded allocation from pathological beds
    if ((size_t)bed_w_words * bed_h_px > MAX_BITMAP_WORDS) {
        BOOST_LOG_TRIVIAL(error) << "BitmapArranger: bed too large for bitmap arrangement. Use convex hull mode.";
        for (auto &ap : arrangables)
            ap.bed_idx = UNARRANGED;
        return;
    }

    BOOST_LOG_TRIVIAL(info) << "BitmapArranger: bed " << bed_w_px << "x" << bed_h_px
                            << " px, resolution " << params.bitmap_resolution_mm << " mm/px"
                            << ", " << arrangables.size() << " items"
                            << (params.consolidate_plates ? " (consolidate)" : "");

    // Consolidate plates: forget existing plate assignments so every item
    // gets re-packed from scratch into the minimum number of plates.
    if (params.consolidate_plates) {
        for (auto &ap : arrangables)
            ap.bed_idx = UNARRANGED;
    }

    struct ItemEntry {
        size_t orig_idx;
        ExPolygon poly;            // convex hull (for fallback / display)
        Points concave_triangles;  // projected triangle verts (3 per tri)
        std::vector<float> concave_z; // Z coords parallel to concave_triangles
        double rotation;
        coord_t inflation;
        int filament_temp_type = -1;
    };

    std::vector<ItemEntry> entries;
    entries.reserve(arrangables.size());
    for (size_t i = 0; i < arrangables.size(); i++) {
        ItemEntry e;
        e.orig_idx = i;
        e.poly = arrangables[i].poly;
        e.concave_triangles = std::move(arrangables[i].concave_triangles);
        e.concave_z = std::move(arrangables[i].concave_z);
        e.rotation = arrangables[i].rotation;
        e.inflation = std::max(inflation, arrangables[i].inflation);
        e.filament_temp_type = arrangables[i].filament_temp_type;
        entries.push_back(std::move(e));
    }

    // Sort largest-first by area: standard BLF companion heuristic.
    // Large items are hardest to place, so they go first while the bed is empty.
    std::sort(entries.begin(), entries.end(), [](const ItemEntry &a, const ItemEntry &b) {
        return std::abs(a.poly.area()) > std::abs(b.poly.area());
    });

    // Build rotation candidates — always 5° steps when rotations enabled.
    // Bitmap rotation is cheap (pixel shuffle from 0° source), so fine granularity is free.
    std::vector<double> rotations = {0.};
    if (params.allow_rotations) {
        rotations.clear();
        double rot_step = PI / 36.; // 5 degrees
        int n_rot = (int)std::round(2.0 * PI / rot_step);
        for (int i = 0; i < n_rot; i++)
            rotations.push_back(i * rot_step);
    }

    // Stamp fixed items and excluded regions (purge line, calibration area) onto bed
    auto stamp_excludes = [&](std::vector<uint64_t> &bits) {
        auto stamp_polys = [&](const ArrangePolygons &polys) {
            for (const auto &excl : polys) {
                ExPolygon ep = excl.poly;
                ep.rotate(excl.rotation);
                ep.translate(excl.translation.x() - effective_bed.min.x(),
                             excl.translation.y() - effective_bed.min.y());
                auto bmp = rasterize(ep, 0, res);
                int ox = (int)std::floor((double)bmp.offset_x / res);
                int oy = (int)std::floor((double)bmp.offset_y / res);
                stamp(bits, bed_w_words, bed_w_px, bed_h_px, bmp, ox, oy);
            }
        };
        stamp_polys(excludes);
        stamp_polys(params.excluded_regions);
    };

    // Plate bitmap vector — one entry per plate, enabling first-fit across
    // all plates instead of forward-only assignment.
    struct PlateState {
        std::vector<uint64_t> bits;
        std::vector<int> skyline;  // 1D height per column for O(bed_w) placement
        int material_group = -1;
        int free_px = 0;
    };

    struct PlateState3D {
        SliceStack stack;
        int material_group = -1;
        int free_px = 0; // based on slice 0
    };

    int total_bed_px = bed_w_px * bed_h_px;

    auto count_free_px = [&](const std::vector<uint64_t> &bits) -> int {
        int used = 0;
        for (auto w : bits) {
            // popcount: count set bits
            #ifdef _MSC_VER
            used += (int)__popcnt64(w);
            #else
            used += __builtin_popcountll(w);
            #endif
        }
        return total_bed_px - used;
    };

    std::vector<PlateState> plates;
    plates.push_back({std::vector<uint64_t>(bed_w_words * bed_h_px, 0), std::vector<int>(bed_w_px, 0), -1, 0});
    stamp_excludes(plates[0].bits);
    plates[0].free_px = count_free_px(plates[0].bits);
    // Initialize skyline from excludes — set skyline height where excludes are stamped
    for (int x = 0; x < bed_w_px; x++) {
        for (int y = bed_h_px - 1; y >= 0; y--) {
            int word = x / 64, bit = x % 64;
            if (plates[0].bits[(size_t)y * bed_w_words + word] & (uint64_t(1) << bit)) {
                plates[0].skyline[x] = y + 1;
                break;
            }
        }
    }

    std::vector<PlateState3D> plates_3d;
    if (params.nesting_3d) {
        PlateState3D p3d;
        p3d.stack.n_slices = 1;
        p3d.stack.slices.push_back(BitmapItem());
        auto &s0 = p3d.stack.slices[0];
        s0.bits.assign(bed_w_words * bed_h_px, 0);
        s0.width_words = bed_w_words;
        s0.width_px = bed_w_px;
        s0.height_px = bed_h_px;
        stamp_excludes(s0.bits);
        plates_3d.push_back(std::move(p3d));
    }

    // Coarse step (~1mm): scan the bed in large strides, then refine within
    // one step of the first collision-free spot. Balances speed vs. packing quality.
    int scan_step = std::max(1, (int)(scaled<coord_t>(1.0) / res));

    // Check whether an item's filament_temp_type is compatible with a plate's
    // assigned material_group. HighTemp (0) and LowTemp (1) cannot coexist.
    // HighLowCompatible (2) and unassigned (-1) are compatible with anything.
    auto is_material_compatible = [&](int plate_group, int item_type) -> bool {
        if (params.allow_multi_materials_on_same_plate)
            return true;
        if (plate_group < 0 || item_type < 0)
            return true; // unassigned is compatible with everything
        // HighLowCompatible (2) and Undefine (3) are compatible with any group
        if (plate_group >= 2 || item_type >= 2)
            return true;
        // Both HighTemp (0) or both LowTemp (1): compatible
        // One HighTemp and one LowTemp: incompatible
        return plate_group == item_type;
    };

    // Try placing an item on a single plate using SKYLINE placement.
    // O(bed_width) per rotation instead of O(bed_width × bed_height).
    auto try_place_on_plate = [&](const ItemEntry &entry,
                                  const std::vector<std::pair<BitmapItem, double>> &rot_bmps,
                                  int plate_idx) -> bool {
        auto &plate = plates[plate_idx];
        for (const auto &[bmp, rot] : rot_bmps) {
            if (bmp.width_px <= 0 || bmp.height_px <= 0) continue;

            // Compute bottom/top profiles for this rotated bitmap
            auto profile = compute_profile(bmp, bed_h_px);
            if (profile.max_y < 0) continue;

            // Skyline placement — O(bed_width)
            auto result = find_placement_skyline(plate.skyline, bed_w_px, bed_h_px, profile);
            if (result) {
                auto [px, py] = *result;

                // Verify with bitmap collision (catches concave edge cases skyline misses)
                if (collides(plate.bits, bed_w_words, bed_w_px, bed_h_px, bmp, px, py)) {
                    // Skyline said OK but bitmap says overlap — skip this rotation.
                    // This can happen with concave shapes where the skyline overestimates free space.
                    continue;
                }

                arrangables[entry.orig_idx].translation = {
                    effective_bed.min.x() + (coord_t)(px * res) - bmp.offset_x,
                    effective_bed.min.y() + (coord_t)(py * res) - bmp.offset_y
                };
                arrangables[entry.orig_idx].rotation = rot;
                arrangables[entry.orig_idx].bed_idx = plate_idx;

                // Stamp bitmap
                stamp(plate.bits, bed_w_words, bed_w_px, bed_h_px, bmp, px, py);

                // Update skyline
                for (const auto &[col, top] : profile.top_pairs) {
                    int c = px + col;
                    int v = py + top;
                    if (c < bed_w_px && v > plate.skyline[c])
                        plate.skyline[c] = v;
                }

                if (plate.material_group < 0)
                    plate.material_group = entry.filament_temp_type;
                return true;
            }
        }
        return false;
    };

    // 3D placement: try all rotated stacks on a 3D plate
    auto try_place_on_plate_3d = [&](const ItemEntry &entry,
                                      const std::vector<std::pair<SliceStack, double>> &rot_stacks,
                                      int plate_idx) -> bool {
        auto &bed_stack = plates_3d[plate_idx].stack;
        for (const auto &[stack, rot] : rot_stacks) {
            auto result = find_placement_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, scan_step);
            if (result) {
                auto [px, py] = *result;
                const auto &s0 = stack.slices[0];
                arrangables[entry.orig_idx].translation = {
                    effective_bed.min.x() + (coord_t)(px * res) - s0.offset_x,
                    effective_bed.min.y() + (coord_t)(py * res) - s0.offset_y
                };
                arrangables[entry.orig_idx].rotation = rot;
                arrangables[entry.orig_idx].bed_idx = plate_idx;
                stamp_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, px, py, stamp_excludes);
                if (plates_3d[plate_idx].material_group < 0)
                    plates_3d[plate_idx].material_group = entry.filament_temp_type;
                return true;
            }
        }
        return false;
    };

    constexpr int MAX_PLATES = 100;
    int n = (int)entries.size();

    // ============================================================
    // PHASE 1: Compute shapes — rasterize each item's triangles at 0°
    // ============================================================
    std::vector<BitmapItem> base_bitmaps(n);
    std::vector<SliceStack> base_stacks(n);
    for (int i = 0; i < n; i++) {
        auto &entry = entries[i];
        if (params.nesting_3d && !entry.concave_z.empty()) {
            // 3D path: rasterize into a multi-slice stack
            base_stacks[i] = rasterize_slices(entry.concave_triangles, entry.concave_z,
                                              entry.inflation, res,
                                              params.slice_height_mm, params.z_clearance_mm);
            // Also produce a 2D bitmap for pixel-count estimates
            base_bitmaps[i] = rasterize_triangles(entry.concave_triangles, entry.inflation, res);
        } else if (!entry.concave_triangles.empty()) {
            base_bitmaps[i] = rasterize_triangles(entry.concave_triangles, entry.inflation, res);
        } else {
            base_bitmaps[i] = rasterize(entry.poly, entry.inflation, res);
        }
        if (params.progressind)
            params.progressind(i + 1, " (computing shapes)");
        if (params.stopcondition && params.stopcondition())
            return;
    }

    // ============================================================
    // PHASE 2: Prepare rotations — rotate each 0° bitmap to all candidate angles
    // ============================================================
    // rot_bmps_all[i] = vector of (BitmapItem, angle) for item i
    std::vector<std::vector<std::pair<BitmapItem, double>>> rot_bmps_all(n);
    // rot_stacks_all[i] = vector of (SliceStack, angle) for 3D items
    std::vector<std::vector<std::pair<SliceStack, double>>> rot_stacks_all(n);
    for (int i = 0; i < n; i++) {
        if (!base_stacks[i].slices.empty()) {
            // 3D item: rotate the stack instead of the bitmap
            rot_stacks_all[i].push_back({base_stacks[i], 0.});
            for (size_t ri = 1; ri < rotations.size(); ri++) {
                auto rstack = rotate_stack(base_stacks[i], rotations[ri], res);
                rot_stacks_all[i].push_back({std::move(rstack), rotations[ri]});
            }
            // Still build 2D rotations for pixel-count estimates
            auto &base = base_bitmaps[i];
            if (base.width_px > 0 && base.height_px > 0) {
                rot_bmps_all[i].push_back({base, 0.});
            }
        } else {
            auto &base = base_bitmaps[i];
            if (base.width_px <= 0 || base.height_px <= 0) continue;

            rot_bmps_all[i].push_back({base, 0.});
            for (size_t ri = 1; ri < rotations.size(); ri++) {
                auto rbmp = rotate_bitmap(base, rotations[ri], res);
                if (rbmp.width_px > 0 && rbmp.height_px > 0)
                    rot_bmps_all[i].push_back({std::move(rbmp), rotations[ri]});
            }
        }

        if (params.progressind)
            params.progressind(i + 1, " (preparing rotations)");
        if (params.stopcondition && params.stopcondition())
            return;
    }

    // Free the triangle and Z data — no longer needed after rasterization
    for (auto &entry : entries) {
        entry.concave_triangles.clear();
        entry.concave_z.clear();
    }

    // Precompute pixel count estimates for quick plate skip
    std::vector<int> item_est_px(n, 0);
    for (int i = 0; i < n; i++) {
        if (rot_bmps_all[i].empty()) continue;
        for (auto w : rot_bmps_all[i][0].first.bits) {
            #ifdef _MSC_VER
            item_est_px[i] += (int)__popcnt64(w);
            #else
            item_est_px[i] += __builtin_popcountll(w);
            #endif
        }
    }

    // ============================================================
    // PHASE 3: Place parts — PLATE-CENTRIC filling
    // Fill one plate completely (largest to smallest), then move on.
    // Once the smallest remaining item can't fit, the plate is full.
    // No backtracking, no re-scanning closed plates.
    // ============================================================
    int placed_count = 0;
    int failed_count = 0;

    // Track which items still need placement
    std::vector<bool> item_placed(n, false);
    for (int i = 0; i < n; i++) {
        if (rot_bmps_all[i].empty() && rot_stacks_all[i].empty()) {
            arrangables[entries[i].orig_idx].bed_idx = -1;
            item_placed[i] = true;
            failed_count++;
        }
    }

    // Helper to create a new 2D plate
    auto new_2d_plate = [&]() -> int {
        int idx = (int)plates.size();
        plates.push_back({std::vector<uint64_t>(bed_w_words * bed_h_px, 0), std::vector<int>(bed_w_px, 0), -1, 0});
        stamp_excludes(plates[idx].bits);
        plates[idx].free_px = count_free_px(plates[idx].bits);
        // Initialize skyline from excludes
        for (int x = 0; x < bed_w_px; x++) {
            for (int y = bed_h_px - 1; y >= 0; y--) {
                int word = x / 64, bit = x % 64;
                if (plates[idx].bits[(size_t)y * bed_w_words + word] & (uint64_t(1) << bit)) {
                    plates[idx].skyline[x] = y + 1;
                    break;
                }
            }
        }
        return idx;
    };

    // Helper to create a new 3D plate
    auto new_3d_plate = [&]() -> int {
        PlateState3D p3d;
        p3d.stack.n_slices = 1;
        BitmapItem s0;
        s0.bits.assign(bed_w_words * bed_h_px, 0);
        s0.width_words = bed_w_words;
        s0.width_px = bed_w_px;
        s0.height_px = bed_h_px;
        stamp_excludes(s0.bits);
        p3d.stack.slices.push_back(std::move(s0));
        int idx = (int)plates_3d.size();
        plates_3d.push_back(std::move(p3d));
        return idx;
    };

    bool use_3d = params.nesting_3d;
    int current_plate = use_3d ? 0 : 0; // first plate already exists

    for (;;) {
        bool any_placed_this_plate = false;
        bool smallest_failed = false;

        // Sweep through all remaining items (already sorted largest-first)
        for (int i = 0; i < n; i++) {
            if (item_placed[i]) continue;
            if (params.stopcondition && params.stopcondition()) goto done;

            auto &entry = entries[i];
            bool placed = false;

            if (use_3d && !rot_stacks_all[i].empty()) {
                if (!is_material_compatible(plates_3d[current_plate].material_group, entry.filament_temp_type))
                    continue;
                placed = try_place_on_plate_3d(entry, rot_stacks_all[i], current_plate);
            } else {
                if (!is_material_compatible(plates[current_plate].material_group, entry.filament_temp_type))
                    continue;
                placed = try_place_on_plate(entry, rot_bmps_all[i], current_plate);
            }

            if (placed) {
                item_placed[i] = true;
                placed_count++;
                any_placed_this_plate = true;
            }
            // If not placed, skip — a smaller item might still fit

            if (params.progressind)
                params.progressind(placed_count + failed_count, " (placing parts)");
        }

        // Check: did we place anything this sweep? If not, plate is full.
        // Also check: are there remaining items?
        int remaining = 0;
        for (int i = 0; i < n; i++)
            if (!item_placed[i]) remaining++;

        if (remaining == 0) break; // all placed

        if (!any_placed_this_plate) {
            // Nothing fit on this plate. If multi-plate, create a new one.
            if (!params.allow_multi_plate || (use_3d ? (int)plates_3d.size() : (int)plates.size()) >= MAX_PLATES) {
                // Can't create more plates — mark remaining as unarranged
                for (int i = 0; i < n; i++) {
                    if (!item_placed[i]) {
                        arrangables[entries[i].orig_idx].bed_idx = -1;
                        item_placed[i] = true;
                        failed_count++;
                    }
                }
                break;
            }
            current_plate = use_3d ? new_3d_plate() : new_2d_plate();
        }
        // If we placed things but some remain, sweep the same plate again —
        // smaller items might now fit in gaps left by the items we just placed.
        // But if nothing new was placed, we already move to a new plate above.
    }
    done:

    // ============================================================
    // PHASE 4: Compaction — configurable strategies to reduce plate count
    // ============================================================

    // Helper: try to place an item on a plate using adaptive-stride bitmap scan.
    // Stride scales with item bounding box: large items use fine stride (few positions),
    // small items use coarse stride (many positions but gaps are plentiful).
    // Based on tournament-winning quadtree/adaptive approach.
    auto try_compact_bitmap = [&](int idx, const BitmapItem &bmp, double rot, int plate_idx) -> bool {
        int max_x = bed_w_px - bmp.width_px;
        int max_y = bed_h_px - bmp.height_px;
        if (max_x < 0 || max_y < 0) return false;

        // Adaptive stride: scale with position space
        int pos_space = (max_x + 1) * (max_y + 1);
        int step;
        if (pos_space < 40000)       step = std::max(2, (int)(scaled<coord_t>(2.0) / res));
        else if (pos_space < 120000) step = std::max(2, (int)(scaled<coord_t>(4.0) / res));
        else                         step = std::max(2, (int)(scaled<coord_t>(6.0) / res));

        // Precompute row-level free pixel counts for band rejection
        int item_px = 0;
        for (int py = 0; py < bmp.height_px; py++)
            for (int wx = 0; wx < bmp.width_words; wx++) {
                uint64_t w = bmp.bits[(size_t)py * bmp.width_words + wx];
                while (w) { item_px++; w &= w - 1; } // popcount
            }

        for (int y = 0; y <= max_y; y += step) {
            // Band free-pixel check: skip if rows can't possibly hold the item
            int band_free = 0;
            for (int dy = 0; dy < bmp.height_px && band_free < item_px; dy++) {
                for (int wx = 0; wx < bed_w_words; wx++) {
                    uint64_t w = ~plates[plate_idx].bits[(size_t)(y + dy) * bed_w_words + wx];
                    while (w) { band_free++; w &= w - 1; }
                    if (band_free >= item_px) break;
                }
            }
            if (band_free < item_px) continue;

            for (int x = 0; x <= max_x; x += step) {
                if (!collides(plates[plate_idx].bits, bed_w_words, bed_w_px, bed_h_px, bmp, x, y)) {
                    arrangables[entries[idx].orig_idx].translation = {
                        effective_bed.min.x() + (coord_t)(x * res) - bmp.offset_x,
                        effective_bed.min.y() + (coord_t)(y * res) - bmp.offset_y
                    };
                    arrangables[entries[idx].orig_idx].rotation = rot;
                    arrangables[entries[idx].orig_idx].bed_idx = plate_idx;
                    stamp(plates[plate_idx].bits, bed_w_words, bed_w_px, bed_h_px, bmp, x, y);
                    auto prof = compute_profile(bmp, bed_h_px);
                    for (const auto &[col, top] : prof.top_pairs) {
                        int c = x + col;
                        int v = y + top;
                        if (c < bed_w_px && v > plates[plate_idx].skyline[c])
                            plates[plate_idx].skyline[c] = v;
                    }
                    return true;
                }
            }
        }
        return false;
    };

    // Helper: try to place using reverse skyline (scan from bottom, find lowest free row per column)
    auto try_compact_reverse_skyline = [&](int idx, const BitmapItem &bmp, double rot, int plate_idx) -> bool {
        if (bmp.width_px <= 0 || bmp.height_px <= 0) return false;
        auto profile = compute_profile(bmp, bed_h_px);
        if (profile.max_y < 0) return false;

        // Compute reverse skyline: lowest free pixel per column
        std::vector<int> rev_sky(bed_w_px, 0);
        for (int x = 0; x < bed_w_px; x++) {
            rev_sky[x] = plates[plate_idx].skyline[x]; // start from forward skyline
            // Scan downward for any gap
            for (int y = 0; y < plates[plate_idx].skyline[x]; y++) {
                int word = x / 64, bit = x % 64;
                if (!(plates[plate_idx].bits[(size_t)y * bed_w_words + word] & (uint64_t(1) << bit))) {
                    rev_sky[x] = y;
                    break;
                }
            }
        }

        // Slide item top profile across reverse skyline
        int n_pos = bed_w_px - profile.bw + 1;
        if (n_pos <= 0) return false;

        int best_x = -1, best_y = INT_MAX;
        for (int x = 0; x < n_pos; x += 8) { // stride 8
            int y = 0;
            for (int c = 0; c < profile.bw; c++) {
                int diff = rev_sky[x + c] - profile.bottom[c];
                if (diff > y) y = diff;
            }
            if (y < 0) y = 0;
            if (y <= profile.max_y && y < best_y) {
                best_y = y; best_x = x;
            }
        }

        if (best_x < 0 || best_y > profile.max_y) return false;

        // Verify with bitmap collision
        if (collides(plates[plate_idx].bits, bed_w_words, bed_w_px, bed_h_px, bmp, best_x, best_y))
            return false;

        arrangables[entries[idx].orig_idx].translation = {
            effective_bed.min.x() + (coord_t)(best_x * res) - bmp.offset_x,
            effective_bed.min.y() + (coord_t)(best_y * res) - bmp.offset_y
        };
        arrangables[entries[idx].orig_idx].rotation = rot;
        arrangables[entries[idx].orig_idx].bed_idx = plate_idx;
        stamp(plates[plate_idx].bits, bed_w_words, bed_w_px, bed_h_px, bmp, best_x, best_y);
        for (const auto &[col, top] : profile.top_pairs) {
            int c = best_x + col;
            int v = best_y + top;
            if (c < bed_w_px && v > plates[plate_idx].skyline[c])
                plates[plate_idx].skyline[c] = v;
        }
        return true;
    };

    // Run compaction if enabled (works for both 2D and 3D modes)
    // For 3D mode, compaction uses 2D bitmap collision on the bottom slice —
    // if the bottom slice doesn't collide, the item fits (conservative but fast).
    if (plates.size() > 1 && params.compaction_mode > 0) {
        bool compacted = true;
        while (compacted && plates.size() > 1) {
            compacted = false;
            int last_plate = (int)plates.size() - 1;

            std::vector<int> last_plate_items;
            for (int i = n - 1; i >= 0; i--) {
                if (!item_placed[i]) continue;
                if (arrangables[entries[i].orig_idx].bed_idx == last_plate)
                    last_plate_items.push_back(i);
            }
            if (last_plate_items.empty()) break;

            int moved = 0;
            for (int idx : last_plate_items) {
                bool relocated = false;

                // Best-fit: try all plates, pick lowest Y. First-fit: take first match.
                int best_plate = -1, best_y_score = INT_MAX;
                int best_bmp_idx = -1;
                double best_rot = 0;

                for (int pi = 0; pi < last_plate; pi++) {
                    if (!is_material_compatible(plates[pi].material_group, entries[idx].filament_temp_type))
                        continue;

                    for (size_t ri = 0; ri < rot_bmps_all[idx].size(); ri++) {
                        const auto &[bmp, rot] = rot_bmps_all[idx][ri];
                        if (bmp.width_px <= 0 || bmp.height_px <= 0) continue;

                        bool found = false;

                        // Strategy 2: reverse skyline
                        if (params.compaction_mode == 2 || params.compaction_mode == 3) {
                            // For best-fit we'd need to not commit yet — but reverse skyline
                            // commits on success. For now, use first-fit with reverse skyline.
                            if (!params.best_fit_compact) {
                                if (try_compact_reverse_skyline(idx, bmp, rot, pi)) {
                                    relocated = true;
                                    moved++;
                                    found = true;
                                }
                            }
                        }

                        // Strategy 1: coarse bitmap scan
                        if (!found && (params.compaction_mode == 1 || params.compaction_mode == 3)) {
                            if (!params.best_fit_compact) {
                                if (try_compact_bitmap(idx, bmp, rot, pi)) {
                                    relocated = true;
                                    moved++;
                                    found = true;
                                }
                            }
                        }

                        // For best-fit: just check skyline placement Y without committing
                        if (params.best_fit_compact && !relocated) {
                            auto profile = compute_profile(bmp, bed_h_px);
                            if (profile.max_y >= 0) {
                                auto result = find_placement_skyline(plates[pi].skyline, bed_w_px, bed_h_px, profile);
                                if (result) {
                                    auto [px, py] = *result;
                                    if (py < best_y_score) {
                                        best_y_score = py;
                                        best_plate = pi;
                                        best_bmp_idx = (int)ri;
                                        best_rot = rot;
                                    }
                                }
                            }
                        }

                        if (relocated) break;
                    }
                    if (relocated) break;
                }

                // Best-fit: commit the best placement found
                if (!relocated && params.best_fit_compact && best_plate >= 0) {
                    const auto &[bmp, rot] = rot_bmps_all[idx][best_bmp_idx];
                    // Try reverse skyline first, then bitmap scan on the best plate
                    if (params.compaction_mode >= 2)
                        relocated = try_compact_reverse_skyline(idx, bmp, rot, best_plate);
                    if (!relocated && params.compaction_mode != 2)
                        relocated = try_compact_bitmap(idx, bmp, rot, best_plate);
                    if (relocated) moved++;
                }
            }

            if (moved > 0 && moved == (int)last_plate_items.size()) {
                plates.pop_back();
                compacted = true;
                BOOST_LOG_TRIVIAL(info) << "BitmapArranger: compacted — removed plate " << last_plate;
            } else {
                compacted = false;
            }
        }
    }

    // ============================================================
    // PHASE 5: Gravity — slide parts toward Y=0 to tighten layout
    // ============================================================
    if (!use_3d && params.gravity_compact) {
        for (int i = 0; i < n; i++) {
            if (!item_placed[i]) continue;
            int pi = arrangables[entries[i].orig_idx].bed_idx;
            if (pi < 0 || pi >= (int)plates.size()) continue;
            if (rot_bmps_all[i].empty()) continue;

            // Find the bitmap that matches the current rotation
            double cur_rot = arrangables[entries[i].orig_idx].rotation;
            const BitmapItem *cur_bmp = nullptr;
            for (const auto &[bmp, rot] : rot_bmps_all[i]) {
                if (std::abs(rot - cur_rot) < 0.01) { cur_bmp = &bmp; break; }
            }
            if (!cur_bmp || cur_bmp->width_px <= 0) continue;

            // Current position in pixels
            coord_t tx = arrangables[entries[i].orig_idx].translation.x();
            coord_t ty = arrangables[entries[i].orig_idx].translation.y();
            int cur_px = (int)((tx + cur_bmp->offset_x - effective_bed.min.x()) / res);
            int cur_py = (int)((ty + cur_bmp->offset_y - effective_bed.min.y()) / res);

            // Erase from bed bitmap
            // (simple approach: rebuild bed from scratch for this plate — expensive but correct)
            // Instead: just try lower Y positions and re-stamp if we find one
            for (int try_y = 0; try_y < cur_py; try_y++) {
                if (!collides(plates[pi].bits, bed_w_words, bed_w_px, bed_h_px, *cur_bmp, cur_px, try_y)) {
                    // Found a lower valid position — but we can't move without un-stamping first
                    // Skip gravity for now if it requires un-stamping (complex)
                    // TODO: implement proper un-stamp + re-stamp
                    break;
                }
            }
        }
    }

    // ============================================================
    // PHASE 6: Overlap safety check — paranoid failsafe
    // ============================================================
    // Re-rasterize all placed items per plate and verify no pixel overlaps.
    // If any overlap is found, move the offending item to a new plate.
    // This should never trigger — if it does, there's a bug upstream.
    // But overlapping parts on a 14-hour print is a catastrophic failure,
    // so we check anyway. Cost: O(items × item_pixels), negligible.
    {
        int n_plates_check = use_3d ? (int)plates_3d.size() : (int)plates.size();
        // Rebuild bed bitmaps from scratch per plate
        for (int pi = 0; pi < n_plates_check; pi++) {
            std::vector<uint64_t> verify_bits(bed_w_words * bed_h_px, 0);
            stamp_excludes(verify_bits);

            for (int i = 0; i < n; i++) {
                if (!item_placed[i]) continue;
                int item_plate = arrangables[entries[i].orig_idx].bed_idx;
                if (item_plate != pi) continue;
                if (rot_bmps_all[i].empty()) continue;

                // Find the bitmap matching the committed rotation
                double cur_rot = arrangables[entries[i].orig_idx].rotation;
                const BitmapItem *cur_bmp = nullptr;
                for (const auto &[bmp, rot] : rot_bmps_all[i]) {
                    if (std::abs(rot - cur_rot) < 0.01) { cur_bmp = &bmp; break; }
                }
                if (!cur_bmp || cur_bmp->width_px <= 0) continue;

                // Compute pixel position from translation
                coord_t tx = arrangables[entries[i].orig_idx].translation.x();
                coord_t ty = arrangables[entries[i].orig_idx].translation.y();
                int px = (int)((tx + cur_bmp->offset_x - effective_bed.min.x()) / res);
                int py = (int)((ty + cur_bmp->offset_y - effective_bed.min.y()) / res);

                if (collides(verify_bits, bed_w_words, bed_w_px, bed_h_px, *cur_bmp, px, py)) {
                    // Overlap detected — move item to unarranged
                    BOOST_LOG_TRIVIAL(error) << "BitmapArranger: OVERLAP detected for item "
                                             << entries[i].orig_idx << " on plate " << pi
                                             << " — moving to new plate as failsafe";
                    arrangables[entries[i].orig_idx].bed_idx = n_plates_check;
                } else {
                    stamp(verify_bits, bed_w_words, bed_w_px, bed_h_px, *cur_bmp, px, py);
                }
            }
        }
    }

    int total_plates = params.nesting_3d ? (int)plates_3d.size() : (int)plates.size();
    BOOST_LOG_TRIVIAL(info) << "BitmapArranger: placed " << placed_count
                            << "/" << arrangables.size() << " on "
                            << total_plates << " plate(s)"
                            << (params.nesting_3d ? " (3D nesting)" : "");
}

}} // namespace Slic3r::arrangement
