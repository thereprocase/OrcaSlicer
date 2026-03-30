#include "BitmapArranger.hpp"
#include "libslic3r/ClipperUtils.hpp"

#include <boost/log/trivial.hpp>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <atomic>
#include <mutex>
#include <set>
#include <chrono>
#define ARRANGE_LOG(msg) BOOST_LOG_TRIVIAL(info) << msg

// Coordinate systems used in this module:
//   Scaled coords (coord_t): Slic3r internal units. 1mm = scaled(1.0).
//   Pixel coords (int): px = (scaled_coord - bed_origin) / resolution.
//     Pixel (0,0) = effective bed min corner (after purge pad shrink).
//   Placement: ArrangePolygon.translation = bed_min + (px * res) - item_bb_min
//     Applied AFTER ArrangePolygon.rotation.

// 64M × 8 bytes = 512 MB max per plate. Larger beds must use convex-hull mode.
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

// Fill bits [px_left..px_right] in row py of a bitmap with width_words.
static inline void fill_span(std::vector<uint64_t> &bits, int width_words, int py, int px_left, int px_right) {
    int sw = px_left / 64, ew = px_right / 64;
    int sb = px_left % 64, eb = px_right % 64;
    size_t row = (size_t)py * width_words;
    if (sw == ew) {
        bits[row + sw] |= ((uint64_t(2) << eb) - 1) & ~((uint64_t(1) << sb) - 1);
    } else {
        bits[row + sw] |= ~((uint64_t(1) << sb) - 1);
        for (int w = sw + 1; w < ew; w++) bits[row + w] = ~uint64_t(0);
        bits[row + ew] |= (uint64_t(2) << eb) - 1;
    }
}

// Dilate all set pixels by inflate_px in all directions.
static void dilate_bitmap(std::vector<uint64_t> &bits, int width_words, int width_px, int height_px, int inflate_px) {
    std::vector<uint64_t> dilated = bits;
    for (int py = 0; py < height_px; py++) {
        for (int wx = 0; wx < width_words; wx++) {
            uint64_t word = bits[(size_t)py * width_words + wx];
            if (word == 0) continue;
            while (word) {
                int bit = ctz64(word);
                int px = wx * 64 + bit;
                for (int dy = -inflate_px; dy <= inflate_px; dy++) {
                    int ny = py + dy;
                    if (ny < 0 || ny >= height_px) continue;
                    for (int dx = -inflate_px; dx <= inflate_px; dx++) {
                        int nx = px + dx;
                        if (nx < 0 || nx >= width_px) continue;
                        dilated[(size_t)ny * width_words + nx / 64] |= (uint64_t(1) << (nx % 64));
                    }
                }
                word &= word - 1;
            }
        }
    }
    bits = std::move(dilated);
}

namespace Slic3r { namespace arrangement {

BitmapArranger::BitmapItem BitmapArranger::rasterize(
    const ExPolygon &poly, coord_t inflation, coord_t res)
{
    BitmapItem item;

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

    if (item.width_px <= 0 || item.height_px <= 0) {
        item.width_px = item.height_px = item.width_words = 0;
        return item;
    }

    item.width_words = (item.width_px + 63) / 64;

    item.bits.assign(item.width_words * item.height_px, 0);

    std::vector<coord_t> x_intersections;
    for (int py = 0; py < item.height_px; py++) {
        coord_t y = bb.min.y() + (coord_t)(py * res + res / 2);

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

        for (size_t i = 0; i + 1 < x_intersections.size(); i += 2) {
            int px_start = std::max(0, (int)((x_intersections[i] - bb.min.x()) / res));
            int px_end   = std::min(item.width_px - 1,
                                    (int)((x_intersections[i + 1] - bb.min.x()) / res));
            fill_span(item.bits, item.width_words, py, px_start, px_end);
        }
    }

    return item;
}

BitmapArranger::BitmapItem BitmapArranger::rasterize_triangles(
    const Points &tri_verts, coord_t inflation, coord_t res)
{
    BitmapItem item;
    if (tri_verts.size() < 3) return item;

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

    for (size_t ti = 0; ti + 2 < tri_verts.size(); ti += 3) {
        double x0 = (tri_verts[ti].x()     - bb.min.x()) * inv_res;
        double y0 = (tri_verts[ti].y()     - bb.min.y()) * inv_res;
        double x1 = (tri_verts[ti + 1].x() - bb.min.x()) * inv_res;
        double y1 = (tri_verts[ti + 1].y() - bb.min.y()) * inv_res;
        double x2 = (tri_verts[ti + 2].x() - bb.min.x()) * inv_res;
        double y2 = (tri_verts[ti + 2].y() - bb.min.y()) * inv_res;

        if (y0 > y1) { std::swap(x0, x1); std::swap(y0, y1); }
        if (y0 > y2) { std::swap(x0, x2); std::swap(y0, y2); }
        if (y1 > y2) { std::swap(x1, x2); std::swap(y1, y2); }

        int py_min = std::max(0, (int)y0);
        int py_max = std::min(item.height_px - 1, (int)y2);

        for (int py = py_min; py <= py_max; py++) {
            double row = py + 0.5;
            double left = 1e18, right = -1e18;

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

            if (px_left <= px_right)
                fill_span(item.bits, item.width_words, py, px_left, px_right);
        }
    }

    if (inflation > 0) {
        int inflate_px = std::max(1, (int)(inflation * inv_res));
        dilate_bitmap(item.bits, item.width_words, item.width_px, item.height_px, inflate_px);
    }

    return item;
}

BitmapArranger::BitmapItem BitmapArranger::rotate_bitmap(
    const BitmapItem &src, double angle_rad, coord_t res)
{
    if (src.width_px <= 0 || src.height_px <= 0)
        return src;

    double cos_a = std::cos(angle_rad), sin_a = std::sin(angle_rad);

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

    // Rotation is around the object's own center, so the center stays the same
    double src_cx = src.offset_x + hw * res;
    double src_cy = src.offset_y + hh * res;
    dst.offset_x = (coord_t)(src_cx - (dst.width_px * 0.5) * res);
    dst.offset_y = (coord_t)(src_cy - (dst.height_px * 0.5) * res);

    double dst_cx = dst.width_px * 0.5;
    double dst_cy = dst.height_px * 0.5;

    for (int dy = 0; dy < dst.height_px; dy++) {
        double fy = dy - dst_cy;
        for (int dx = 0; dx < dst.width_px; dx++) {
            double fx = dx - dst_cx;
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

// --- 3D-aware nesting: Z-slice collision ---

BitmapArranger::SliceStack BitmapArranger::rasterize_slices(
    const Points &tri_verts, const std::vector<float> &tri_z,
    coord_t inflation, coord_t res,
    float slice_height_mm, float z_clearance_mm)
{
    SliceStack stack;
    if (tri_verts.size() < 3 || tri_z.size() != tri_verts.size()) return stack;

    float max_z = 0;
    for (float z : tri_z) max_z = std::max(max_z, z);
    stack.n_slices = std::max(1, (int)std::ceil(max_z / slice_height_mm));

    // All slices must share this same coordinate frame so that collides_3d
    // can pass the same (ox, oy) to collides() for every slice. Without this,
    // each slice has a different offset/size and the collision check is misaligned.
    BoundingBox global_bb(tri_verts);
    global_bb.min -= Vec2crd(inflation, inflation);
    global_bb.max += Vec2crd(inflation, inflation);

    int global_w = (int)std::ceil((double)(global_bb.max.x() - global_bb.min.x()) / res) + 1;
    int global_h = (int)std::ceil((double)(global_bb.max.y() - global_bb.min.y()) / res) + 1;
    int global_width_words = (global_w + 63) / 64;

    if (global_w <= 0 || global_h <= 0) return stack;

    // Bin triangles into Z slices, expanded by z_clearance
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

    double inv_res = 1.0 / res;

    auto rasterize_into_global = [&](const Points &tris) -> BitmapItem {
        BitmapItem item;
        item.offset_x = global_bb.min.x();
        item.offset_y = global_bb.min.y();
        item.width_px = global_w;
        item.height_px = global_h;
        item.width_words = global_width_words;
        item.bits.assign((size_t)global_width_words * global_h, 0);

        // Why not call rasterize_triangles()? It computes its own bounding box per call,
        // giving each slice a different coordinate frame. collides_3d needs all slices in the same frame.
        for (size_t ti = 0; ti + 2 < tris.size(); ti += 3) {
            double x0 = (tris[ti].x()     - global_bb.min.x()) * inv_res;
            double y0 = (tris[ti].y()     - global_bb.min.y()) * inv_res;
            double x1 = (tris[ti + 1].x() - global_bb.min.x()) * inv_res;
            double y1 = (tris[ti + 1].y() - global_bb.min.y()) * inv_res;
            double x2 = (tris[ti + 2].x() - global_bb.min.x()) * inv_res;
            double y2 = (tris[ti + 2].y() - global_bb.min.y()) * inv_res;

            if (y0 > y1) { std::swap(x0, x1); std::swap(y0, y1); }
            if (y0 > y2) { std::swap(x0, x2); std::swap(y0, y2); }
            if (y1 > y2) { std::swap(x1, x2); std::swap(y1, y2); }

            int py_min = std::max(0, (int)y0);
            int py_max = std::min(global_h - 1, (int)y2);

            for (int py = py_min; py <= py_max; py++) {
                double row = py + 0.5;
                double left = 1e18, right = -1e18;

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
                int px_right = std::min(global_w - 1, (int)right);

                if (px_left <= px_right)
                    fill_span(item.bits, global_width_words, py, px_left, px_right);
            }
        }

        if (inflation > 0) {
            int inflate_px = std::max(1, (int)(inflation * inv_res));
            dilate_bitmap(item.bits, global_width_words, global_w, global_h, inflate_px);
        }

        return item;
    };

    stack.slices.resize(stack.n_slices);
    if (stack.n_slices >= 4) {
        tbb::parallel_for(tbb::blocked_range<int>(0, stack.n_slices),
            [&](const tbb::blocked_range<int> &range) {
                for (int s = range.begin(); s < range.end(); s++) {
                    if (!slice_tris[s].empty())
                        stack.slices[s] = rasterize_into_global(slice_tris[s]);
                    else {
                        // Empty slice still needs correct dimensions for collides_3d
                        stack.slices[s].offset_x = global_bb.min.x();
                        stack.slices[s].offset_y = global_bb.min.y();
                        stack.slices[s].width_px = global_w;
                        stack.slices[s].height_px = global_h;
                        stack.slices[s].width_words = global_width_words;
                    }
                }
            }
        );
    } else {
        for (int s = 0; s < stack.n_slices; s++) {
            if (!slice_tris[s].empty())
                stack.slices[s] = rasterize_into_global(slice_tris[s]);
            else {
                stack.slices[s].offset_x = global_bb.min.x();
                stack.slices[s].offset_y = global_bb.min.y();
                stack.slices[s].width_px = global_w;
                stack.slices[s].height_px = global_h;
                stack.slices[s].width_words = global_width_words;
            }
        }
    }

    // Collapse identical slices to 1 — prismatic parts (cubes, cylinders) don't need
    // per-slice collision. Compare bits, not triangle counts (top/bottom faces differ).
    if (stack.n_slices > 1 && !stack.slices[0].bits.empty()) {
        bool all_same = true;
        const auto &ref_bits = stack.slices[0].bits;
        int ref_w = stack.slices[0].width_px;
        int ref_h = stack.slices[0].height_px;
        for (int s = 1; s < stack.n_slices && all_same; s++) {
            if (stack.slices[s].width_px != ref_w ||
                stack.slices[s].height_px != ref_h ||
                stack.slices[s].bits.size() != ref_bits.size()) {
                all_same = false;
                break;
            }
            if (stack.slices[s].bits != ref_bits)
                all_same = false;
        }
        if (all_same) {
            BitmapItem s0 = std::move(stack.slices[0]);
            stack.slices.clear();
            stack.slices.push_back(std::move(s0));
            BOOST_LOG_TRIVIAL(info) << "BitmapArranger: collapsed " << stack.n_slices
                                    << " identical slices to 1";
        }
    }

    return stack;
}

bool BitmapArranger::collides_3d(
    const SliceStack &bed_stack, int bed_w, int bed_w_px, int bed_h,
    const SliceStack &item_stack, int ox, int oy)
{
    bool item_collapsed = (int)item_stack.slices.size() < item_stack.n_slices;
    bool bed_collapsed = (int)bed_stack.slices.size() < bed_stack.n_slices;

    // Fast path: both collapsed → one 2D check covers all Z levels
    if (item_collapsed && bed_collapsed) {
        if (!item_stack.slices.empty() && item_stack.slices[0].width_px > 0 &&
            !bed_stack.slices.empty() && bed_stack.slices[0].width_px > 0) {
            return collides(bed_stack.slices[0].bits, bed_w, bed_w_px, bed_h,
                            item_stack.slices[0], ox, oy);
        }
        return false;
    }

    int max_z = std::min(item_stack.n_slices,
                         std::max((int)item_stack.slices.size(), (int)bed_stack.slices.size()));
    for (int z = 0; z < max_z; z++) {
        int item_z = item_collapsed ? 0 : z;
        int bed_z = bed_collapsed ? 0 : z;
        if (item_z >= (int)item_stack.slices.size()) break;
        if (item_stack.slices[item_z].width_px <= 0) continue;
        if (bed_z < (int)bed_stack.slices.size() && bed_stack.slices[bed_z].width_px > 0) {
            if (collides(bed_stack.slices[bed_z].bits, bed_w, bed_w_px, bed_h,
                         item_stack.slices[item_z], ox, oy))
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
    while (bed_stack.n_slices < item_stack.n_slices) {
        BitmapItem new_slice;
        new_slice.bits.assign(bed_w * bed_h, 0);
        new_slice.width_words = bed_w;
        new_slice.width_px = bed_w_px;
        new_slice.height_px = bed_h;
        // Physical obstacles (cutter, fixed items) apply at all Z levels.
        // Bed boundary is enforced implicitly by bitmap bounds in collides().
        if (init_slice) init_slice(new_slice.bits);
        bed_stack.slices.push_back(std::move(new_slice));
        bed_stack.n_slices++;
    }

    bool item_collapsed = (int)item_stack.slices.size() < item_stack.n_slices;
    for (int z = 0; z < item_stack.n_slices; z++) {
        int item_z = item_collapsed ? 0 : z;
        if (item_z >= (int)item_stack.slices.size()) break;
        if (item_stack.slices[item_z].width_px <= 0) continue;
        if (bed_stack.slices[z].bits.empty()) {
            bed_stack.slices[z].bits.assign(bed_w * bed_h, 0);
            bed_stack.slices[z].width_words = bed_w;
            bed_stack.slices[z].width_px = bed_w_px;
            bed_stack.slices[z].height_px = bed_h;
            if (init_slice) init_slice(bed_stack.slices[z].bits);
        }
        stamp(bed_stack.slices[z].bits, bed_stack.slices[z].width_words,
              bed_w_px, bed_h, item_stack.slices[item_z], ox, oy);
    }
}

std::optional<std::pair<int,int>> BitmapArranger::find_placement_3d(
    const SliceStack &bed_stack, int bed_w, int bed_w_px, int bed_h,
    const SliceStack &item_stack, int step)
{
    // Only call expensive collides_3d() at positions that already pass
    // the cheap 2D collision check on slice 0.
    if (item_stack.slices.empty() || bed_stack.slices.empty())
        return std::nullopt;

    int item_w = 0, item_h = 0;
    for (const auto &s : item_stack.slices) {
        item_w = std::max(item_w, s.width_px);
        item_h = std::max(item_h, s.height_px);
    }
    int max_x = bed_w_px - item_w;
    int max_y = bed_h - item_h;
    if (max_x < 0 || max_y < 0) return std::nullopt;

    const auto &bed_s0 = bed_stack.slices[0];
    const auto &item_s0 = item_stack.slices[0];

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

                // 2D pre-filter on slice 0 avoids expensive full-3D check
                if (collides(bed_s0.bits, bed_w, bed_w_px, bed_h, item_s0, x, y))
                    continue;

                if (!collides_3d(bed_stack, bed_w, bed_w_px, bed_h, item_stack, x, y)) {
                    int best_x = x, best_y = y;
                    int best_dist = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                    for (int ry = std::max(0, y - step + 1); ry <= std::min(max_y, y + step - 1); ry++) {
                        for (int rx = std::max(0, x - step + 1); rx <= std::min(max_x, x + step - 1); rx++) {
                            int d = (rx - cx) * (rx - cx) + (ry - cy) * (ry - cy);
                            if (d < best_dist &&
                                !collides(bed_s0.bits, bed_w, bed_w_px, bed_h, item_s0, rx, ry) &&
                                !collides_3d(bed_stack, bed_w, bed_w_px, bed_h, item_stack, rx, ry)) {
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

    // Collapsed stacks (prismatic parts) have only one physical slice — rotate it once.
    int actual_slices = (int)src.slices.size();
    dst.slices.resize(actual_slices);

    if (actual_slices >= 4) {
        tbb::parallel_for(tbb::blocked_range<int>(0, actual_slices),
            [&](const tbb::blocked_range<int> &range) {
                for (int s = range.begin(); s < range.end(); s++) {
                    if (src.slices[s].width_px > 0)
                        dst.slices[s] = rotate_bitmap(src.slices[s], angle_rad, res);
                }
            }
        );
    } else {
        for (int s = 0; s < actual_slices; s++) {
            if (src.slices[s].width_px > 0)
                dst.slices[s] = rotate_bitmap(src.slices[s], angle_rad, res);
        }
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
                // Left-edge: shifting by a negative amount is UB, fall back to per-bit
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

            // Left-edge: shifting by a negative amount is UB, fall back to per-bit
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

// --- Skyline placement — O(bed_width) per item ---

BitmapArranger::ItemProfile BitmapArranger::compute_profile(
    const BitmapItem &item, int bed_h)
{
    ItemProfile prof;
    prof.bw = item.width_px;
    prof.bh = item.height_px;
    prof.max_y = bed_h - item.height_px;
    if (prof.max_y < 0) return prof;

    prof.bottom.assign(prof.bw, bed_h + 1); // bed_h+1 = sentinel for empty columns
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

    // Exclusive top: one past the highest set pixel (matches skyline convention)
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
    const ItemProfile &profile, int placement_bias)
{
    if (profile.max_y < 0 || profile.bw > bed_w_px)
        return std::nullopt;

    int n_pos = bed_w_px - profile.bw + 1;
    if (n_pos <= 0) return std::nullopt;

    constexpr int STRIDE = 8;
    constexpr int REFINE = 10;

    auto eval_y = [&](int x) -> int {
        int y = 0;
        for (int c = 0; c < profile.bw; c++) {
            int diff = skyline[x + c] - profile.bottom[c];
            if (diff > y) y = diff;
        }
        return std::max(y, 0);
    };

    if (placement_bias == 1) {
        int best_x = -1, best_y = INT_MAX;

        for (int x = 0; x < n_pos; x += STRIDE) {
            int y = eval_y(x);
            if (y < best_y) {
                best_y = y;
                best_x = x;
                if (y == 0) break; // floor placement, can't improve
            }
        }

        if (best_x < 0 || best_y > profile.max_y)
            return std::nullopt;

        if (best_y > 0) {
            int ref_lo = std::max(0, best_x - REFINE);
            int ref_hi = std::min(n_pos, best_x + REFINE + 1);
            for (int x = ref_lo; x < ref_hi; x++) {
                int y = eval_y(x);
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

    int bed_cx = bed_w_px / 2;

    struct PosY { int x, y; };
    std::vector<PosY> coarse;
    coarse.reserve(n_pos / STRIDE + 1);
    int global_min_y = INT_MAX;

    for (int x = 0; x < n_pos; x += STRIDE) {
        int y = eval_y(x);
        if (y <= profile.max_y) {
            coarse.push_back({x, y});
            if (y < global_min_y) global_min_y = y;
        }
    }

    if (coarse.empty())
        return std::nullopt;

    // Allow candidates within a small Y tolerance of the minimum.
    // This lets center-biased selection avoid a 1-pixel-better corner
    // position when a near-equal center position exists.
    int y_tolerance = std::max(2, global_min_y / 10);
    int y_threshold = global_min_y + y_tolerance;

    int best_x = -1, best_y = INT_MAX;
    int best_dist = INT_MAX;

    for (const auto &p : coarse) {
        if (p.y > y_threshold) continue;
        int item_cx = p.x + profile.bw / 2; // item center X at this position
        int dist = std::abs(item_cx - bed_cx);
        if (dist < best_dist || (dist == best_dist && p.y < best_y)) {
            best_dist = dist;
            best_y = p.y;
            best_x = p.x;
        }
    }

    if (best_x < 0)
        return std::nullopt;

    {
        int ref_lo = std::max(0, best_x - REFINE);
        int ref_hi = std::min(n_pos, best_x + REFINE + 1);
        for (int x = ref_lo; x < ref_hi; x++) {
            int y = eval_y(x);
            if (y > profile.max_y) continue;
            int item_cx = x + profile.bw / 2;
            int dist = std::abs(item_cx - bed_cx);
            if (dist < best_dist || (dist == best_dist && y < best_y)) {
                best_dist = dist;
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
    int max_x = bed_w_px - item.width_px;
    int max_y = bed_h - item.height_px;
    if (max_x < 0 || max_y < 0) return std::nullopt;

    int cx = max_x / 2;
    int cy = max_y / 2;

    int coarse_step = step * 4;
    int hit_x = -1, hit_y = -1;
    int best_dist = INT_MAX;

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
    return std::nullopt;

tier2:
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
    const Polygon &bed_shape,
    const ArrangeParams &params)
{
    BoundingBox bb = get_extents(bed_shape);

    Polygon bb_rect;
    bb_rect.points = {
        {bb.min.x(), bb.min.y()},
        {bb.max.x(), bb.min.y()},
        {bb.max.x(), bb.max.y()},
        {bb.min.x(), bb.max.y()}
    };

    Polygon bed_hole = bed_shape;
    bed_hole.make_clockwise(); // holes must be clockwise

    ExPolygon outside_bed(bb_rect);
    outside_bed.holes.push_back(bed_hole);

    ArrangePolygons augmented_excludes = excludes;
    ArrangePolygon bed_mask_exclude;
    bed_mask_exclude.poly = outside_bed;
    bed_mask_exclude.translation = {0, 0};
    bed_mask_exclude.rotation = 0;
    bed_mask_exclude.bed_idx = 0;
    augmented_excludes.push_back(bed_mask_exclude);

    arrange(arrangables, augmented_excludes, bb, params);
}

void BitmapArranger::arrange(
    ArrangePolygons &arrangables,
    const ArrangePolygons &excludes,
    const BoundingBox &bed,
    const ArrangeParams &params)
{
    coord_t res = scaled<coord_t>(std::clamp(params.bitmap_resolution_mm, 0.1f, 2.0f));
    if (res < scaled<coord_t>(0.3))
        res = scaled<coord_t>(0.3);
    coord_t inflation = params.min_obj_distance / 2;

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

    // Bed dimensions in pixels. No +1 here — the bed bitmap must represent
    // exactly the physical bed area. Item bitmaps use +1 to capture their
    // rightmost edge, but the bed must not extend past its boundary or items
    // placed at the rightmost pixel will overshoot the physical bed edge.
    int bed_w_px = (int)std::ceil((double)(effective_bed.max.x() - effective_bed.min.x()) / res);
    int bed_h_px = (int)std::ceil((double)(effective_bed.max.y() - effective_bed.min.y()) / res);
    int bed_w_words = (bed_w_px + 63) / 64;

    // Guard against pathological bed sizes (e.g. unit mismatch sending mm as scaled coords)
    if ((size_t)bed_w_words * bed_h_px > MAX_BITMAP_WORDS) {
        BOOST_LOG_TRIVIAL(error) << "BitmapArranger: bed too large for bitmap arrangement. Use convex hull mode.";
        for (auto &ap : arrangables)
            ap.bed_idx = UNARRANGED;
        return;
    }

    ARRANGE_LOG("bed " << bed_w_px << "x" << bed_h_px
                << " px, res " << params.bitmap_resolution_mm << " mm/px"
                << ", " << arrangables.size() << " items"
                << (params.nesting_3d ? ", 3D" : ", 2D")
                << (params.allow_rotations ? (", rot " + std::to_string((int)std::round(params.rotation_step_rad * 180.0 / PI)) + " deg") : ", no rot")
                << ", compact=" << params.compaction_mode
                << (params.allow_multi_materials_on_same_plate ? ", multi-mat" : ", split-mat")
                << (params.consolidate_plates ? ", consolidate" : ""));

    auto t_total_start = std::chrono::high_resolution_clock::now();
    auto t_phase_start = t_total_start;

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
        std::set<int> extrude_ids; // extruder indices used by this object
        int preferred_plate = -1; // keep-plates mode: only place on this plate
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
        e.extrude_ids = std::set<int>(arrangables[i].extrude_ids.begin(),
                                       arrangables[i].extrude_ids.end());
        e.preferred_plate = arrangables[i].bed_idx;
        entries.push_back(std::move(e));
    }

    // Largest-first: standard BLF heuristic — large items are hardest to place.
    std::sort(entries.begin(), entries.end(), [](const ItemEntry &a, const ItemEntry &b) {
        return std::abs(a.poly.area()) > std::abs(b.poly.area());
    });

    // Minimum 5 degrees to prevent absurd rotation counts
    std::vector<double> rotations = {0.};
    if (params.allow_rotations) {
        rotations.clear();
        double rot_step = std::max(params.rotation_step_rad, PI / 36.); // minimum 5°
        int n_rot = std::max(1, (int)std::round(2.0 * PI / rot_step));
        for (int i = 0; i < n_rot; i++)
            rotations.push_back(i * rot_step);
    }

    auto stamp_excludes = [&](std::vector<uint64_t> &bits) {
        auto stamp_polys = [&](const ArrangePolygons &polys) {
            for (const auto &excl : polys) {
                ExPolygon ep = excl.poly;
                ep.rotate(excl.rotation);
                ep.translate(excl.translation.x() - effective_bed.min.x(),
                             excl.translation.y() - effective_bed.min.y());
                auto bmp = rasterize(ep, scaled<coord_t>(1.0), res); // 1mm physical padding around exclusion zones (resolution-dependent)
                int ox = (int)std::floor((double)bmp.offset_x / res);
                int oy = (int)std::floor((double)bmp.offset_y / res);
                stamp(bits, bed_w_words, bed_w_px, bed_h_px, bmp, ox, oy);
            }
        };
        stamp_polys(excludes);
        stamp_polys(params.excluded_regions);
    };

    struct PlateState {
        std::vector<uint64_t> bits;
        std::vector<int> skyline;  // 1D height per column for O(bed_w) placement
        int material_group = -1;
        int free_px = 0;
    };

    struct PlateState3D {
        SliceStack stack;
        std::vector<int> skyline;  // cached from bottom slice; updated incrementally after stamp_3d
        int material_group = -1;
        int free_px = 0; // based on slice 0
    };

    int total_bed_px = bed_w_px * bed_h_px;

    auto count_free_px = [&](const std::vector<uint64_t> &bits) -> int {
        int used = 0;
        for (auto w : bits) {
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
        p3d.skyline.assign(bed_w_px, 0);
        for (int x = 0; x < bed_w_px; x++) {
            for (int y = bed_h_px - 1; y >= 0; y--) {
                int word = x / 64, bit = x % 64;
                if (s0.bits[(size_t)y * bed_w_words + word] & (uint64_t(1) << bit)) {
                    p3d.skyline[x] = y + 1;
                    break;
                }
            }
        }
        plates_3d.push_back(std::move(p3d));
    }

    // Ensure a 3D plate's bed stack has at least n_slices entries, each
    // stamped with excludes. Exclusion zones (filament cutter, LIDAR, purge
    // line) are physical structures at ALL heights — the bed must have
    // exclude data at every Z level an item might occupy. Without this,
    // collides_3d skips height levels where the bed has no data, allowing
    // items to be placed over exclusion zones at height.
    auto ensure_bed_height = [&](SliceStack &bed_stack, int required_slices) {
        while ((int)bed_stack.slices.size() < required_slices) {
            BitmapItem new_slice;
            new_slice.bits.assign(bed_w_words * bed_h_px, 0);
            new_slice.width_words = bed_w_words;
            new_slice.width_px = bed_w_px;
            new_slice.height_px = bed_h_px;
            stamp_excludes(new_slice.bits);
            bed_stack.slices.push_back(std::move(new_slice));
        }
        bed_stack.n_slices = std::max(bed_stack.n_slices, required_slices);
    };

    int scan_step = std::max(1, (int)(scaled<coord_t>(1.0) / res));
    int scan_step_3d = std::max(2, (int)(scaled<coord_t>(2.0) / res));

    struct PlateExtruders {
        std::set<int> extruder_ids;
    };
    std::vector<PlateExtruders> plate_extruders(1); // one per plate, grows with plates

    // Must match libnest2d's behavior: check both temperature group and extruder identity
    auto is_material_compatible = [&](int plate_group, int item_type,
                                      int plate_idx, const std::set<int> &item_extruders) -> bool {
        if (params.allow_multi_materials_on_same_plate)
            return true;

        if (plate_group >= 0 && item_type >= 0 &&
            plate_group < 2 && item_type < 2 &&
            plate_group != item_type)
            return false;

        // Extruder identity: item's extruders must be a subset of (or superset of)
        // the plate's extruders, matching the libnest2d behavior.
        if (plate_idx >= 0 && plate_idx < (int)plate_extruders.size()) {
            const auto &pe = plate_extruders[plate_idx].extruder_ids;
            if (!pe.empty() && !item_extruders.empty()) {
                bool item_subset = std::includes(pe.begin(), pe.end(),
                                                  item_extruders.begin(), item_extruders.end());
                bool plate_subset = std::includes(item_extruders.begin(), item_extruders.end(),
                                                   pe.begin(), pe.end());
                if (!item_subset && !plate_subset)
                    return false;
            }
        }
        return true;
    };

    struct MinMaxBmps {
        BitmapItem min_bmp, max_bmp;
        ItemProfile min_profile, max_profile;
        bool valid = false;
    };
    std::vector<MinMaxBmps> item_minmax; // sized to n after Phase 2

    auto find_skyline_topk = [&](const std::vector<int> &skyline,
                                  const ItemProfile &profile, int k) {
        std::vector<std::pair<int,int>> candidates;
        if (profile.max_y < 0 || profile.bw > bed_w_px) return candidates;
        int n_pos = bed_w_px - profile.bw + 1;
        if (n_pos <= 0) return candidates;
        constexpr int STRIDE = 8;
        struct PosY { int x, y; };
        std::vector<PosY> valid;
        valid.reserve(n_pos / STRIDE + 1);
        for (int x = 0; x < n_pos; x += STRIDE) {
            int y = 0;
            for (int c = 0; c < profile.bw; c++) {
                int diff = skyline[x + c] - profile.bottom[c];
                if (diff > y) y = diff;
            }
            if (y < 0) y = 0;
            if (y <= profile.max_y)
                valid.push_back({x, y});
        }

        if (params.placement_bias == 0 && !valid.empty()) {
            // Near-minimum Y candidates first, then prefer closest to center
            int min_y = INT_MAX;
            for (const auto &p : valid)
                if (p.y < min_y) min_y = p.y;
            int y_tol = std::max(2, min_y / 10);
            int bed_cx = bed_w_px / 2;
            std::sort(valid.begin(), valid.end(), [&](const PosY &a, const PosY &b) {
                bool a_near = (a.y <= min_y + y_tol);
                bool b_near = (b.y <= min_y + y_tol);
                if (a_near != b_near) return a_near; // near-minimum candidates first
                if (a_near && b_near) {
                    int da = std::abs(a.x + profile.bw / 2 - bed_cx);
                    int db = std::abs(b.x + profile.bw / 2 - bed_cx);
                    if (da != db) return da < db;
                    return a.y < b.y; // tie-break by Y
                }
                return a.y < b.y; // both far from minimum: sort by Y
            });
        } else {
            std::sort(valid.begin(), valid.end(), [](const PosY &a, const PosY &b) {
                return a.y < b.y;
            });
        }

        int take = std::min(k, (int)valid.size());
        for (int i = 0; i < take; i++) {
            int best_x = valid[i].x, best_y = valid[i].y;
            int ref_lo = std::max(0, best_x - 10);
            int ref_hi = std::min(n_pos, best_x + 11);
            for (int x = ref_lo; x < ref_hi; x++) {
                int y = 0;
                for (int c = 0; c < profile.bw; c++) {
                    int diff = skyline[x + c] - profile.bottom[c];
                    if (diff > y) y = diff;
                }
                if (y < 0) y = 0;
                if (y < best_y && y <= profile.max_y) {
                    best_y = y; best_x = x;
                }
            }
            candidates.push_back({best_x, best_y});
        }
        return candidates;
    };

    auto register_extruders = [&](int plate_idx, const std::set<int> &item_extruders) {
        while (plate_idx >= (int)plate_extruders.size())
            plate_extruders.push_back({});
        plate_extruders[plate_idx].extruder_ids.insert(item_extruders.begin(), item_extruders.end());
    };

    // Selects lowest-Y rotation (base tight to cluster keeps overhang inboard)
    auto try_place_on_plate = [&](const ItemEntry &entry, int entry_idx,
                                  const std::vector<std::pair<BitmapItem, double>> &rot_bmps,
                                  int plate_idx) -> bool {
        auto &plate = plates[plate_idx];
        const auto &mm = item_minmax[entry_idx];

        if (mm.valid && mm.max_profile.max_y >= 0) {
            auto candidates = find_skyline_topk(plate.skyline, mm.max_profile, 5);
            for (auto [cx, cy] : candidates) {
                if (collides(plate.bits, bed_w_words, bed_w_px, bed_h_px, mm.max_bmp, cx, cy))
                    continue;

                int best_y = INT_MAX, best_ri = -1;
                bool min_clear = mm.min_bmp.width_px > 0 &&
                    !collides(plate.bits, bed_w_words, bed_w_px, bed_h_px, mm.min_bmp, cx, cy);

                for (size_t ri = 0; ri < rot_bmps.size(); ri++) {
                    const auto &[bmp, rot] = rot_bmps[ri];
                    if (bmp.width_px <= 0) continue;

                    auto profile = compute_profile(bmp, bed_h_px);
                    if (profile.max_y < 0) continue;

                    int ry = 0;
                    for (int c = 0; c < profile.bw && cx + c < bed_w_px; c++) {
                        int diff = plate.skyline[cx + c] - profile.bottom[c];
                        if (diff > ry) ry = diff;
                    }
                    if (ry < 0) ry = 0;
                    if (ry > profile.max_y) continue;

                    if (ry < best_y) {
                        if (min_clear || !collides(plate.bits, bed_w_words, bed_w_px, bed_h_px, bmp, cx, ry)) {
                            best_y = ry;
                            best_ri = (int)ri;
                        }
                    }
                }

                if (best_ri >= 0) {
                    const auto &[bmp, rot] = rot_bmps[best_ri];
                    auto profile = compute_profile(bmp, bed_h_px);

                    arrangables[entry.orig_idx].translation = {
                        effective_bed.min.x() + (coord_t)(cx * res) - bmp.offset_x,
                        effective_bed.min.y() + (coord_t)(best_y * res) - bmp.offset_y
                    };
                    arrangables[entry.orig_idx].rotation = rot;
                    arrangables[entry.orig_idx].bed_idx = plate_idx;
                    stamp(plate.bits, bed_w_words, bed_w_px, bed_h_px, bmp, cx, best_y);
                    for (const auto &[col, top] : profile.top_pairs) {
                        int c = cx + col;
                        int v = best_y + top;
                        if (c < bed_w_px && v > plate.skyline[c])
                            plate.skyline[c] = v;
                    }
                    if (plate.material_group < 0)
                        plate.material_group = entry.filament_temp_type;
                    return true;
                }
            }
        }

        for (const auto &[bmp, rot] : rot_bmps) {
            if (bmp.width_px <= 0 || bmp.height_px <= 0) continue;
            auto profile = compute_profile(bmp, bed_h_px);
            if (profile.max_y < 0) continue;

            auto result = find_placement_skyline(plate.skyline, bed_w_px, bed_h_px, profile, params.placement_bias);
            if (result) {
                auto [px, py] = *result;
                if (collides(plate.bits, bed_w_words, bed_w_px, bed_h_px, bmp, px, py))
                    continue;

                arrangables[entry.orig_idx].translation = {
                    effective_bed.min.x() + (coord_t)(px * res) - bmp.offset_x,
                    effective_bed.min.y() + (coord_t)(py * res) - bmp.offset_y
                };
                arrangables[entry.orig_idx].rotation = rot;
                arrangables[entry.orig_idx].bed_idx = plate_idx;
                stamp(plate.bits, bed_w_words, bed_w_px, bed_h_px, bmp, px, py);
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

    // Rejects go to batch Pass 2 rather than expensive per-rotation center-out fallback
    auto try_place_on_plate_3d = [&](const ItemEntry &entry,
                                      const std::vector<std::pair<SliceStack, double>> &rot_stacks,
                                      int plate_idx) -> bool {
        auto &bed_stack = plates_3d[plate_idx].stack;
        if (bed_stack.slices.empty()) return false;

        int max_item_slices = 0;
        for (const auto &[stack, rot] : rot_stacks)
            max_item_slices = std::max(max_item_slices, stack.n_slices);
        if (max_item_slices > 0)
            ensure_bed_height(bed_stack, max_item_slices);

        const auto &sky = plates_3d[plate_idx].skyline;

        int best_cx = -1, best_cy = INT_MAX, best_ri = -1;

        for (size_t ri = 0; ri < rot_stacks.size(); ri++) {
            const auto &[stack, rot] = rot_stacks[ri];
            if (stack.slices.empty()) continue;
            const auto &item_s0 = stack.slices[0];
            if (item_s0.width_px <= 0 || item_s0.height_px <= 0) continue;

            auto profile = compute_profile(item_s0, bed_h_px);
            if (profile.max_y < 0) continue;

            auto result = find_placement_skyline(sky, bed_w_px, bed_h_px, profile, params.placement_bias);
            if (result) {
                auto [px, py] = *result;
                if (px + item_s0.width_px > bed_w_px) continue;
                if (py + item_s0.height_px > bed_h_px) continue;
                if (py < best_cy && !collides_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, px, py)) {
                    best_cx = px;
                    best_cy = py;
                    best_ri = (int)ri;
                }
            }
        }

        if (best_ri >= 0) {
            const auto &[stack, rot] = rot_stacks[best_ri];
            const auto &item_s0 = stack.slices[0];
            arrangables[entry.orig_idx].translation = {
                effective_bed.min.x() + (coord_t)(best_cx * res) - item_s0.offset_x,
                effective_bed.min.y() + (coord_t)(best_cy * res) - item_s0.offset_y
            };
            arrangables[entry.orig_idx].rotation = rot;
            arrangables[entry.orig_idx].bed_idx = plate_idx;
            stamp_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, best_cx, best_cy, stamp_excludes);

            auto profile = compute_profile(item_s0, bed_h_px);
            for (const auto &[col, top] : profile.top_pairs) {
                int c = best_cx + col;
                int v = best_cy + top;
                if (c < bed_w_px && v > plates_3d[plate_idx].skyline[c])
                    plates_3d[plate_idx].skyline[c] = v;
            }
            if (plates_3d[plate_idx].material_group < 0)
                plates_3d[plate_idx].material_group = entry.filament_temp_type;
            return true;
        }

        for (const auto &[stack, rot] : rot_stacks) {
            if (stack.slices.empty()) continue;
            auto result = find_placement_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, scan_step_3d);
            if (result) {
                auto [px, py] = *result;
                const auto &item_s0 = stack.slices[0];
                arrangables[entry.orig_idx].translation = {
                    effective_bed.min.x() + (coord_t)(px * res) - item_s0.offset_x,
                    effective_bed.min.y() + (coord_t)(py * res) - item_s0.offset_y
                };
                arrangables[entry.orig_idx].rotation = rot;
                arrangables[entry.orig_idx].bed_idx = plate_idx;
                stamp_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, px, py, stamp_excludes);
    
                auto prof = compute_profile(item_s0, bed_h_px);
                for (const auto &[col, top] : prof.top_pairs) {
                    int c = px + col;
                    int v = py + top;
                    if (c < bed_w_px && v > plates_3d[plate_idx].skyline[c])
                        plates_3d[plate_idx].skyline[c] = v;
                }
                if (plates_3d[plate_idx].material_group < 0)
                    plates_3d[plate_idx].material_group = entry.filament_temp_type;
                return true;
            }
        }

        return false;
    };

    // Must match MAX_PLATE_COUNT in src/slic3r/GUI/PartPlate.hpp.
    // Can't include that header from libslic3r (dependency direction).
    constexpr int MAX_PLATES = 36;
    int n = (int)entries.size();

    // --- Phase 1: Compute shapes --- convert triangle meshes to collision bitmaps
    std::vector<BitmapItem> base_bitmaps(n);
    std::vector<SliceStack> base_stacks(n);
    {
        std::atomic<int> rast_progress{0};
        std::atomic<bool> rast_cancelled{false};

        tbb::parallel_for(tbb::blocked_range<int>(0, n),
            [&](const tbb::blocked_range<int> &range) {
                for (int i = range.begin(); i < range.end(); i++) {
                    if (rast_cancelled.load(std::memory_order_relaxed)) break;

                    auto &entry = entries[i];
                    if (params.nesting_3d && !entry.concave_z.empty()) {
                        base_stacks[i] = rasterize_slices(entry.concave_triangles, entry.concave_z,
                                                          entry.inflation, res,
                                                          params.slice_height_mm, params.z_clearance_mm);
                        base_bitmaps[i] = rasterize_triangles(entry.concave_triangles, entry.inflation, res);
                    } else if (!entry.concave_triangles.empty()) {
                        base_bitmaps[i] = rasterize_triangles(entry.concave_triangles, entry.inflation, res);
                    } else {
                        base_bitmaps[i] = rasterize(entry.poly, entry.inflation, res);
                    }

                    rast_progress.fetch_add(1, std::memory_order_relaxed);
                    if (params.stopcondition && params.stopcondition())
                        rast_cancelled.store(true, std::memory_order_relaxed);
                }
            }
        );

        if (rast_cancelled.load()) return;

        // wxWidgets is not thread-safe — report progress from main thread only
        if (params.progressind)
            params.progressind(n, " (computing shapes)");
    }

    {
        auto t_now = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_phase_start).count();
        int n_3d = 0;
        for (int i = 0; i < n; i++)
            if (!base_stacks[i].slices.empty()) n_3d++;
        ARRANGE_LOG("Phase 1 rasterize: " << us << " us, "
                    << n << " items (" << n_3d << " 3D)");
        t_phase_start = t_now;
    }

    // --- Phase 2: Prepare rotations --- pre-rotate bitmaps so placement only tests positions
    std::vector<std::vector<std::pair<BitmapItem, double>>> rot_bmps_all(n);
    std::vector<std::vector<std::pair<SliceStack, double>>> rot_stacks_all(n);

    std::atomic<int> rot_progress{0};
    std::atomic<bool> rot_cancelled{false};

    tbb::parallel_for(tbb::blocked_range<int>(0, n),
        [&](const tbb::blocked_range<int> &range) {
            for (int i = range.begin(); i < range.end(); i++) {
                if (rot_cancelled.load(std::memory_order_relaxed)) break;

                if (!base_stacks[i].slices.empty()) {
                    rot_stacks_all[i].push_back({base_stacks[i], 0.});
                    for (size_t ri = 1; ri < rotations.size(); ri++) {
                        auto rstack = rotate_stack(base_stacks[i], rotations[ri], res);
                        rot_stacks_all[i].push_back({std::move(rstack), rotations[ri]});
                    }
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

                rot_progress.fetch_add(1, std::memory_order_relaxed);
                if (params.stopcondition && params.stopcondition())
                    rot_cancelled.store(true, std::memory_order_relaxed);
            }
        }
    );

    if (rot_cancelled.load()) return;

    if (params.progressind)
        params.progressind(n, " (preparing rotations)");

    {
        auto t_now = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_phase_start).count();
        int total_rots = 0;
        for (int i = 0; i < n; i++)
            total_rots += (int)std::max(rot_bmps_all[i].size(), rot_stacks_all[i].size());
        ARRANGE_LOG("Phase 2 rotate: " << us << " us, "
                    << rotations.size() << " angles x " << n << " items = "
                    << total_rots << " bitmaps");
        t_phase_start = t_now;
    }

    for (auto &entry : entries) {
        entry.concave_triangles.clear();
        entry.concave_z.clear();
    }

    // --- Phase 2.5: Precompute min/max envelopes --- eliminates ~80% of per-rotation collision checks
    // min_bmp = AND of all rotations (core), max_bmp = OR of all rotations (envelope)
    item_minmax.resize(n);

    for (int i = 0; i < n; i++) {
        if (rot_bmps_all[i].size() < 2) continue;

        int env_w = 0, env_h = 0;
        for (const auto &[bmp, rot] : rot_bmps_all[i]) {
            env_w = std::max(env_w, bmp.width_px);
            env_h = std::max(env_h, bmp.height_px);
        }
        if (env_w <= 0 || env_h <= 0) continue;
        int env_ww = (env_w + 63) / 64;

        std::vector<uint64_t> max_bits((size_t)env_ww * env_h, 0);
        std::vector<uint64_t> min_bits((size_t)env_ww * env_h, ~uint64_t(0));

        for (const auto &[bmp, rot] : rot_bmps_all[i]) {
            int dx = (env_w - bmp.width_px) / 2;
            int dy = (env_h - bmp.height_px) / 2;

            std::vector<uint64_t> centered((size_t)env_ww * env_h, 0);
            for (int py = 0; py < bmp.height_px; py++) {
                int ey = py + dy;
                if (ey < 0 || ey >= env_h) continue;
                for (int px = 0; px < bmp.width_px; px++) {
                    int ex = px + dx;
                    if (ex < 0 || ex >= env_w) continue;
                    if (bmp.bits[(size_t)py * bmp.width_words + px / 64] & (uint64_t(1) << (px % 64))) {
                        centered[(size_t)ey * env_ww + ex / 64] |= (uint64_t(1) << (ex % 64));
                    }
                }
            }

            for (size_t w = 0; w < centered.size(); w++) {
                max_bits[w] |= centered[w];
                min_bits[w] &= centered[w];
            }
        }

        auto &ref_bmp = rot_bmps_all[i][0].first;
        coord_t center_x = ref_bmp.offset_x + (coord_t)(ref_bmp.width_px * res / 2);
        coord_t center_y = ref_bmp.offset_y + (coord_t)(ref_bmp.height_px * res / 2);

        BitmapItem max_bmp;
        max_bmp.width_px = env_w;
        max_bmp.height_px = env_h;
        max_bmp.width_words = env_ww;
        max_bmp.bits = std::move(max_bits);
        max_bmp.offset_x = center_x - (coord_t)(env_w * res / 2);
        max_bmp.offset_y = center_y - (coord_t)(env_h * res / 2);

        BitmapItem min_bmp;
        min_bmp.width_px = env_w;
        min_bmp.height_px = env_h;
        min_bmp.width_words = env_ww;
        min_bmp.bits = std::move(min_bits);
        min_bmp.offset_x = max_bmp.offset_x;
        min_bmp.offset_y = max_bmp.offset_y;

        item_minmax[i].max_bmp = std::move(max_bmp);
        item_minmax[i].min_bmp = std::move(min_bmp);
        item_minmax[i].max_profile = compute_profile(item_minmax[i].max_bmp, bed_h_px);
        item_minmax[i].min_profile = compute_profile(item_minmax[i].min_bmp, bed_h_px);
        item_minmax[i].valid = true;
    }

    {
        auto t_now = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_phase_start).count();
        int n_valid = 0;
        for (int i = 0; i < n; i++) if (item_minmax[i].valid) n_valid++;
        ARRANGE_LOG("Phase 2.5 min/max: " << us << " us, " << n_valid << "/" << n << " items");
        t_phase_start = t_now;
    }

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

    // --- Phase 3: Place parts --- fill one plate completely before moving to the next
    int placed_count = 0;
    int failed_count = 0;

    std::vector<bool> item_placed(n, false);
    for (int i = 0; i < n; i++) {
        if (rot_bmps_all[i].empty() && rot_stacks_all[i].empty()) {
            arrangables[entries[i].orig_idx].bed_idx = -1;
            item_placed[i] = true;
            failed_count++;
        }
    }

    auto new_2d_plate = [&]() -> int {
        int idx = (int)plates.size();
        plates.push_back({std::vector<uint64_t>(bed_w_words * bed_h_px, 0), std::vector<int>(bed_w_px, 0), -1, 0});
        stamp_excludes(plates[idx].bits);
        plates[idx].free_px = count_free_px(plates[idx].bits);
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

    auto new_3d_plate = [&]() -> int {
        PlateState3D p3d;
        p3d.stack.n_slices = 1;
        BitmapItem s0;
        s0.bits.assign(bed_w_words * bed_h_px, 0);
        s0.width_words = bed_w_words;
        s0.width_px = bed_w_px;
        s0.height_px = bed_h_px;
        stamp_excludes(s0.bits);
        p3d.skyline.assign(bed_w_px, 0);
        for (int x = 0; x < bed_w_px; x++) {
            for (int y = bed_h_px - 1; y >= 0; y--) {
                int word = x / 64, bit = x % 64;
                if (s0.bits[(size_t)y * bed_w_words + word] & (uint64_t(1) << bit)) {
                    p3d.skyline[x] = y + 1;
                    break;
                }
            }
        }
        p3d.stack.slices.push_back(std::move(s0));
        int idx = (int)plates_3d.size();
        plates_3d.push_back(std::move(p3d));
        return idx;
    };

    bool use_3d = params.nesting_3d;

    int current_plate = 0;

    for (;;) {
        bool any_placed_this_plate = false;

        for (int i = 0; i < n; i++) {
            if (item_placed[i]) continue;
            if (params.stopcondition && params.stopcondition()) goto done;

            auto &entry = entries[i];
            bool placed = false;

            // keep-plates: skip plates that don't match the preferred plate
            if (entry.preferred_plate >= 0 && current_plate != entry.preferred_plate)
                continue;

            if (use_3d && !rot_stacks_all[i].empty()) {
                if (is_material_compatible(
                        plates_3d[current_plate].material_group,
                        entry.filament_temp_type, current_plate, entry.extrude_ids))
                    placed = try_place_on_plate_3d(entry, rot_stacks_all[i], current_plate);
            } else {
                if (is_material_compatible(
                        plates[current_plate].material_group,
                        entry.filament_temp_type, current_plate, entry.extrude_ids))
                    placed = try_place_on_plate(entry, i, rot_bmps_all[i], current_plate);
            }

            if (placed) {
                item_placed[i] = true;
                placed_count++;
                any_placed_this_plate = true;
                register_extruders(current_plate, entry.extrude_ids);
            }

            if (params.progressind)
                params.progressind(placed_count + failed_count, " (placing parts)");
        }

        int remaining = 0;
        for (int i = 0; i < n; i++)
            if (!item_placed[i]) remaining++;

        if (remaining == 0) break;

        if (!any_placed_this_plate) {
            if (!params.allow_multi_plate ||
                (use_3d ? (int)plates_3d.size() : (int)plates.size()) >= MAX_PLATES) {
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
    }
    done:

    // Keep-plates overflow: items that had a preferred plate but couldn't be placed
    for (int i = 0; i < n; i++) {
        if (!item_placed[i] && entries[i].preferred_plate >= 0) {
            auto &entry = entries[i];
            arrangables[entry.orig_idx].bed_idx = entry.preferred_plate;
            arrangables[entry.orig_idx].translation = {effective_bed.min.x(), effective_bed.min.y()};
            item_placed[i] = true;
            failed_count++;
        }
    }

    {
        auto t_now = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_phase_start).count();
        int n_plates_placed = use_3d ? (int)plates_3d.size() : (int)plates.size();
        int n_rejects = 0;
        for (int i = 0; i < n; i++) if (!item_placed[i]) n_rejects++;
        ARRANGE_LOG("Phase 3a skyline: " << us << " us, "
                    << placed_count << " placed, " << n_rejects << " rejects, "
                    << n_plates_placed << " plates");
        t_phase_start = t_now;
    }

    // --- Phase 3b: Batch reject placement --- center-out fallback for skyline failures
    {
        std::vector<int> rejects;
        for (int i = 0; i < n; i++)
            if (!item_placed[i]) rejects.push_back(i);

        if (!rejects.empty()) {
            ARRANGE_LOG("Phase 3b: " << rejects.size() << " rejects, batch placement");

            for (int idx : rejects) {
                if (params.stopcondition && params.stopcondition()) break;
                auto &entry = entries[idx];
                bool placed = false;

                int n_plates_try = use_3d ? (int)plates_3d.size() : (int)plates.size();
                for (int pi = 0; pi < n_plates_try && !placed; pi++) {
                    if (use_3d && !rot_stacks_all[idx].empty()) {
                        if (!is_material_compatible(
                                plates_3d[pi].material_group,
                                entry.filament_temp_type, pi, entry.extrude_ids))
                            continue;

                        auto &bed_stack = plates_3d[pi].stack;
                        int max_slices = 0;
                        for (const auto &[stack, rot] : rot_stacks_all[idx])
                            max_slices = std::max(max_slices, stack.n_slices);
                        if (max_slices > 0) ensure_bed_height(bed_stack, max_slices);

                        for (const auto &[stack, rot] : rot_stacks_all[idx]) {
                            if (stack.slices.empty()) continue;
                            auto result = find_placement_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, scan_step_3d);
                            if (result) {
                                auto [px, py] = *result;
                                const auto &item_s0 = stack.slices[0];
                                arrangables[entry.orig_idx].translation = {
                                    effective_bed.min.x() + (coord_t)(px * res) - item_s0.offset_x,
                                    effective_bed.min.y() + (coord_t)(py * res) - item_s0.offset_y
                                };
                                arrangables[entry.orig_idx].rotation = rot;
                                arrangables[entry.orig_idx].bed_idx = pi;
                                stamp_3d(bed_stack, bed_w_words, bed_w_px, bed_h_px, stack, px, py, stamp_excludes);
                
                                auto prof = compute_profile(item_s0, bed_h_px);
                                for (const auto &[col, top] : prof.top_pairs) {
                                    int c = px + col;
                                    int v = py + top;
                                    if (c < bed_w_px && v > plates_3d[pi].skyline[c])
                                        plates_3d[pi].skyline[c] = v;
                                }
                                if (plates_3d[pi].material_group < 0)
                                    plates_3d[pi].material_group = entry.filament_temp_type;
                                register_extruders(pi, entry.extrude_ids);
                                placed = true;
                                break;
                            }
                        }
                    } else if (!rot_bmps_all[idx].empty()) {
                        if (!is_material_compatible(
                                plates[pi].material_group,
                                entry.filament_temp_type, pi, entry.extrude_ids))
                            continue;

                        for (const auto &[bmp, rot] : rot_bmps_all[idx]) {
                            if (bmp.width_px <= 0) continue;
                            auto result = find_placement(plates[pi].bits, bed_w_words, bed_w_px, bed_h_px, bmp, scan_step);
                            if (result) {
                                auto [px, py] = *result;
                                arrangables[entry.orig_idx].translation = {
                                    effective_bed.min.x() + (coord_t)(px * res) - bmp.offset_x,
                                    effective_bed.min.y() + (coord_t)(py * res) - bmp.offset_y
                                };
                                arrangables[entry.orig_idx].rotation = rot;
                                arrangables[entry.orig_idx].bed_idx = pi;
                                stamp(plates[pi].bits, bed_w_words, bed_w_px, bed_h_px, bmp, px, py);
                                auto profile = compute_profile(bmp, bed_h_px);
                                for (const auto &[col, top] : profile.top_pairs) {
                                    int c = px + col;
                                    int v = py + top;
                                    if (c < bed_w_px && v > plates[pi].skyline[c])
                                        plates[pi].skyline[c] = v;
                                }
                                if (plates[pi].material_group < 0)
                                    plates[pi].material_group = entry.filament_temp_type;
                                register_extruders(pi, entry.extrude_ids);
                                placed = true;
                                break;
                            }
                        }
                    }
                }

                if (placed) {
                    item_placed[idx] = true;
                    placed_count++;
                } else if (params.allow_multi_plate) {
                    int n_cur = use_3d ? (int)plates_3d.size() : (int)plates.size();
                    if (n_cur >= MAX_PLATES) continue; // all plates full
                    int new_pi = use_3d ? new_3d_plate() : new_2d_plate();
                    if (use_3d && !rot_stacks_all[idx].empty()) {
                        if (try_place_on_plate_3d(entry, rot_stacks_all[idx], new_pi)) {
                            item_placed[idx] = true;
                            placed_count++;
                            register_extruders(new_pi, entry.extrude_ids);
                        }
                    } else if (!rot_bmps_all[idx].empty()) {
                        if (try_place_on_plate(entry, idx, rot_bmps_all[idx], new_pi)) {
                            item_placed[idx] = true;
                            placed_count++;
                            register_extruders(new_pi, entry.extrude_ids);
                        }
                    }
                }

                if (params.progressind)
                    params.progressind(placed_count + failed_count, " (placing rejects)");
            }
        }

        auto t_now = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_phase_start).count();
        ARRANGE_LOG("Phase 3b rejects: " << us << " us");
        t_phase_start = t_now;
    }

    // --- Phase 4: Compaction --- relocate items from last plate to earlier plates to free it
    // Thread safety: find phase is read-only (parallel), commit phase is serial.

    auto find_compact_bitmap = [&](const BitmapItem &bmp, int plate_idx)
        -> std::optional<std::pair<int,int>>
    {
        int max_x = bed_w_px - bmp.width_px;
        int max_y = bed_h_px - bmp.height_px;
        if (max_x < 0 || max_y < 0) return std::nullopt;

        int pos_space = (max_x + 1) * (max_y + 1);
        int step;
        if (pos_space < 40000)       step = std::max(2, (int)(scaled<coord_t>(2.0) / res));
        else if (pos_space < 120000) step = std::max(2, (int)(scaled<coord_t>(4.0) / res));
        else                         step = std::max(2, (int)(scaled<coord_t>(6.0) / res));

        int item_px = 0;
        for (int py = 0; py < bmp.height_px; py++)
            for (int wx = 0; wx < bmp.width_words; wx++) {
                uint64_t w = bmp.bits[(size_t)py * bmp.width_words + wx];
                while (w) { item_px++; w &= w - 1; }
            }

        for (int y = 0; y <= max_y; y += step) {
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
                if (!collides(plates[plate_idx].bits, bed_w_words, bed_w_px, bed_h_px, bmp, x, y))
                    return std::make_pair(x, y);
            }
        }
        return std::nullopt;
    };

    auto find_compact_reverse_skyline = [&](const BitmapItem &bmp, int plate_idx)
        -> std::optional<std::pair<int,int>>
    {
        if (bmp.width_px <= 0 || bmp.height_px <= 0) return std::nullopt;
        auto profile = compute_profile(bmp, bed_h_px);
        if (profile.max_y < 0) return std::nullopt;

        std::vector<int> rev_sky(bed_w_px, 0);
        for (int x = 0; x < bed_w_px; x++) {
            rev_sky[x] = plates[plate_idx].skyline[x];
            for (int y = 0; y < plates[plate_idx].skyline[x]; y++) {
                int word = x / 64, bit = x % 64;
                if (!(plates[plate_idx].bits[(size_t)y * bed_w_words + word] & (uint64_t(1) << bit))) {
                    rev_sky[x] = y;
                    break;
                }
            }
        }

        int n_pos = bed_w_px - profile.bw + 1;
        if (n_pos <= 0) return std::nullopt;

        int best_x = -1, best_y = INT_MAX;
        for (int x = 0; x < n_pos; x += 8) {
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

        if (best_x < 0 || best_y > profile.max_y) return std::nullopt;
        if (collides(plates[plate_idx].bits, bed_w_words, bed_w_px, bed_h_px, bmp, best_x, best_y))
            return std::nullopt;

        return std::make_pair(best_x, best_y);
    };

    auto commit_compact = [&](int idx, const BitmapItem &bmp, double rot, int plate_idx, int px, int py) {
        register_extruders(plate_idx, entries[idx].extrude_ids);
        arrangables[entries[idx].orig_idx].translation = {
            effective_bed.min.x() + (coord_t)(px * res) - bmp.offset_x,
            effective_bed.min.y() + (coord_t)(py * res) - bmp.offset_y
        };
        arrangables[entries[idx].orig_idx].rotation = rot;
        arrangables[entries[idx].orig_idx].bed_idx = plate_idx;
        stamp(plates[plate_idx].bits, bed_w_words, bed_w_px, bed_h_px, bmp, px, py);
        auto prof = compute_profile(bmp, bed_h_px);
        for (const auto &[col, top] : prof.top_pairs) {
            int c = px + col;
            int v = py + top;
            if (c < bed_w_px && v > plates[plate_idx].skyline[c])
                plates[plate_idx].skyline[c] = v;
        }
    };

    if (use_3d && plates_3d.size() > 1 && params.compaction_mode > 0) {
        bool compacted = true;
        while (compacted && plates_3d.size() > 1) {
            compacted = false;
            int last_plate = (int)plates_3d.size() - 1;

            std::vector<int> last_plate_items;
            for (int i = n - 1; i >= 0; i--) {
                if (!item_placed[i]) continue;
                if (arrangables[entries[i].orig_idx].bed_idx == last_plate)
                    last_plate_items.push_back(i);
            }
            if (last_plate_items.empty()) break;

            int moved = 0;
            for (int idx : last_plate_items) {
                if (params.stopcondition && params.stopcondition()) break;
                if (rot_stacks_all[idx].empty()) continue;

                bool relocated = false;

                int max_item_slices = 0;
                for (const auto &[stack, rot] : rot_stacks_all[idx])
                    max_item_slices = std::max(max_item_slices, stack.n_slices);

                for (int pi = 0; pi < last_plate && !relocated; pi++) {
                    if (!is_material_compatible(plates_3d[pi].material_group, entries[idx].filament_temp_type, pi, entries[idx].extrude_ids))
                        continue;

                    ensure_bed_height(plates_3d[pi].stack, max_item_slices);

                    for (const auto &[stack, rot] : rot_stacks_all[idx]) {
                        if (stack.slices.empty()) continue;
                        const auto &item_s0 = stack.slices[0];
                        if (item_s0.width_px <= 0 || item_s0.height_px <= 0) continue;

                        auto profile = compute_profile(item_s0, bed_h_px);
                        if (profile.max_y < 0) continue;

                        auto result = find_placement_skyline(plates_3d[pi].skyline, bed_w_px, bed_h_px, profile, params.placement_bias);
                        if (result) {
                            auto [px, py] = *result;
                            if (!collides_3d(plates_3d[pi].stack, bed_w_words, bed_w_px, bed_h_px, stack, px, py)) {
                                arrangables[entries[idx].orig_idx].translation = {
                                    effective_bed.min.x() + (coord_t)(px * res) - item_s0.offset_x,
                                    effective_bed.min.y() + (coord_t)(py * res) - item_s0.offset_y
                                };
                                arrangables[entries[idx].orig_idx].rotation = rot;
                                arrangables[entries[idx].orig_idx].bed_idx = pi;
                                stamp_3d(plates_3d[pi].stack, bed_w_words, bed_w_px, bed_h_px, stack, px, py, stamp_excludes);
                
                                for (const auto &[col, top] : profile.top_pairs) {
                                    int c = px + col;
                                    int v = py + top;
                                    if (c < bed_w_px && v > plates_3d[pi].skyline[c])
                                        plates_3d[pi].skyline[c] = v;
                                }
                                register_extruders(pi, entries[idx].extrude_ids);
                                relocated = true;
                                moved++;
                                break;
                            }
                        }

                        auto fallback = find_placement_3d(plates_3d[pi].stack, bed_w_words, bed_w_px, bed_h_px, stack, scan_step_3d);
                        if (fallback) {
                            auto [px, py] = *fallback;
                            arrangables[entries[idx].orig_idx].translation = {
                                effective_bed.min.x() + (coord_t)(px * res) - item_s0.offset_x,
                                effective_bed.min.y() + (coord_t)(py * res) - item_s0.offset_y
                            };
                            arrangables[entries[idx].orig_idx].rotation = rot;
                            arrangables[entries[idx].orig_idx].bed_idx = pi;
                            stamp_3d(plates_3d[pi].stack, bed_w_words, bed_w_px, bed_h_px, stack, px, py, stamp_excludes);
            
                            for (const auto &[col, top] : profile.top_pairs) {
                                int c = px + col;
                                int v = py + top;
                                if (c < bed_w_px && v > plates_3d[pi].skyline[c])
                                    plates_3d[pi].skyline[c] = v;
                            }
                            register_extruders(pi, entries[idx].extrude_ids);
                            relocated = true;
                            moved++;
                            break;
                        }
                    }
                }

                if (params.progressind)
                    params.progressind(placed_count + failed_count, " (compacting 3D)");
            }

            if (moved > 0 && moved == (int)last_plate_items.size()) {
                plates_3d.pop_back();
                compacted = true;
                BOOST_LOG_TRIVIAL(info) << "BitmapArranger: 3D compacted — removed plate " << last_plate;
            } else {
                compacted = false;
            }
        }
    }

    if (!use_3d && plates.size() > 1 && params.compaction_mode > 0) {
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
                if (params.stopcondition && params.stopcondition()) break;

                struct Candidate { const BitmapItem *bmp; double rot; size_t ri; };
                std::vector<Candidate> candidates;
                for (size_t ri = 0; ri < rot_bmps_all[idx].size(); ri++) {
                    const auto &[bmp, rot] = rot_bmps_all[idx][ri];
                    if (bmp.width_px > 0 && bmp.height_px > 0)
                        candidates.push_back({&bmp, rot, ri});
                }
                if (candidates.empty()) continue;

                std::vector<int> compat_plates;
                for (int pi = 0; pi < last_plate; pi++) {
                    if (is_material_compatible(plates[pi].material_group, entries[idx].filament_temp_type, pi, entries[idx].extrude_ids))
                        compat_plates.push_back(pi);
                }
                if (compat_plates.empty()) continue;

                struct SearchResult {
                    int plate_idx = -1;
                    int px = 0, py = 0;
                    size_t ri = 0;
                    double rot = 0;
                };

                bool relocated = false;

                if (params.best_fit_compact) {
                    std::vector<SearchResult> results(compat_plates.size());
                    std::atomic<bool> any_found{false};

                    tbb::parallel_for(tbb::blocked_range<size_t>(0, compat_plates.size()),
                        [&](const tbb::blocked_range<size_t> &range) {
                            for (size_t pi_idx = range.begin(); pi_idx < range.end(); pi_idx++) {
                                int pi = compat_plates[pi_idx];
                                for (const auto &cand : candidates) {
                                    std::optional<std::pair<int,int>> pos;
                                    if (params.compaction_mode >= 2)
                                        pos = find_compact_reverse_skyline(*cand.bmp, pi);
                                    if (!pos && params.compaction_mode != 2)
                                        pos = find_compact_bitmap(*cand.bmp, pi);
                                    if (pos) {
                                        results[pi_idx] = {pi, pos->first, pos->second, cand.ri, cand.rot};
                                        any_found.store(true, std::memory_order_relaxed);
                                        break; // found on this plate, try next plate for best-fit
                                    }
                                }
                            }
                        }
                    );

                    if (any_found.load()) {
                        int best_idx = -1, best_y = INT_MAX;
                        for (size_t i = 0; i < results.size(); i++) {
                            if (results[i].plate_idx >= 0 && results[i].py < best_y) {
                                best_y = results[i].py;
                                best_idx = (int)i;
                            }
                        }
                        if (best_idx >= 0) {
                            auto &r = results[best_idx];
    
                            for (const auto &cand : candidates) {
                                if (cand.ri == r.ri) {
                                    if (!collides(plates[r.plate_idx].bits, bed_w_words, bed_w_px, bed_h_px,
                                                  *cand.bmp, r.px, r.py)) {
                                        commit_compact(idx, *cand.bmp, r.rot, r.plate_idx, r.px, r.py);
                                        relocated = true;
                                        moved++;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                } else {
                    std::atomic<bool> found_flag{false};
                    SearchResult winner;
                    std::mutex winner_mutex;

                    tbb::parallel_for(tbb::blocked_range<size_t>(0, compat_plates.size()),
                        [&](const tbb::blocked_range<size_t> &range) {
                            for (size_t pi_idx = range.begin(); pi_idx < range.end(); pi_idx++) {
                                if (found_flag.load(std::memory_order_relaxed)) return;
                                int pi = compat_plates[pi_idx];
                                for (const auto &cand : candidates) {
                                    if (found_flag.load(std::memory_order_relaxed)) return;
                                    std::optional<std::pair<int,int>> pos;
                                    if (params.compaction_mode >= 2)
                                        pos = find_compact_reverse_skyline(*cand.bmp, pi);
                                    if (!pos && params.compaction_mode != 2)
                                        pos = find_compact_bitmap(*cand.bmp, pi);
                                    if (pos) {
                                        std::lock_guard<std::mutex> lock(winner_mutex);
                                        if (!found_flag.load()) {
                                            winner = {pi, pos->first, pos->second, cand.ri, cand.rot};
                                            found_flag.store(true, std::memory_order_release);
                                        }
                                        return;
                                    }
                                }
                            }
                        }
                    );

                    if (found_flag.load()) {
                        for (const auto &cand : candidates) {
                            if (cand.ri == winner.ri) {

                                if (!collides(plates[winner.plate_idx].bits, bed_w_words, bed_w_px, bed_h_px,
                                              *cand.bmp, winner.px, winner.py)) {
                                    commit_compact(idx, *cand.bmp, winner.rot, winner.plate_idx, winner.px, winner.py);
                                    relocated = true;
                                    moved++;
                                }
                                break;
                            }
                        }
                    }
                }

                if (params.progressind)
                    params.progressind(placed_count + failed_count, " (compacting)");
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

    {
        auto t_now = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_phase_start).count();
        int n_plates_after = use_3d ? (int)plates_3d.size() : (int)plates.size();
        ARRANGE_LOG("Phase 4 compact: " << us << " us, "
                    << n_plates_after << " plates after compaction");
        t_phase_start = t_now;
    }

    // PHASE 5: Gravity compaction — not yet implemented.
    // Requires un-stamp + re-stamp to move items without leaving ghost pixels.
    // Placeholder: the UI option exists but does nothing until this is built.

    // --- Phase 6: Overlap safety check --- re-verify all placements to catch compaction bugs
    // Re-rasterize all placed items per plate and verify no pixel overlaps.
    // If any overlap is found, mark the item as unarranged (Phase 7 rescues it).
    {
        int n_plates_check = use_3d ? (int)plates_3d.size() : (int)plates.size();
        for (int pi = 0; pi < n_plates_check; pi++) {
            if (use_3d) {
                SliceStack verify_stack;
                verify_stack.n_slices = 1;
                BitmapItem vs0;
                vs0.bits.assign(bed_w_words * bed_h_px, 0);
                vs0.width_words = bed_w_words;
                vs0.width_px = bed_w_px;
                vs0.height_px = bed_h_px;
                stamp_excludes(vs0.bits);
                verify_stack.slices.push_back(std::move(vs0));

                for (int i = 0; i < n; i++) {
                    if (!item_placed[i]) continue;
                    int item_plate = arrangables[entries[i].orig_idx].bed_idx;
                    if (item_plate != pi) continue;
                    if (rot_stacks_all[i].empty()) continue;

                    double cur_rot = arrangables[entries[i].orig_idx].rotation;
                    const SliceStack *cur_stack = nullptr;
                    for (const auto &[stack, rot] : rot_stacks_all[i]) {
                        if (std::abs(rot - cur_rot) < 0.01) { cur_stack = &stack; break; }
                    }
                    if (!cur_stack || cur_stack->slices.empty()) continue;
                    const auto &s0 = cur_stack->slices[0];
                    if (s0.width_px <= 0) continue;

                    coord_t tx = arrangables[entries[i].orig_idx].translation.x();
                    coord_t ty = arrangables[entries[i].orig_idx].translation.y();
                    int px = (int)((tx + s0.offset_x - effective_bed.min.x()) / res);
                    int py = (int)((ty + s0.offset_y - effective_bed.min.y()) / res);

                    if (collides_3d(verify_stack, bed_w_words, bed_w_px, bed_h_px, *cur_stack, px, py)) {
                        BOOST_LOG_TRIVIAL(error) << "BitmapArranger: 3D OVERLAP detected for item "
                                                 << entries[i].orig_idx << " on plate " << pi;
                        arrangables[entries[i].orig_idx].bed_idx = -1;
                    } else {
                        stamp_3d(verify_stack, bed_w_words, bed_w_px, bed_h_px, *cur_stack, px, py, stamp_excludes);
                    }
                }
            } else {
                std::vector<uint64_t> verify_bits(bed_w_words * bed_h_px, 0);
                stamp_excludes(verify_bits);

                for (int i = 0; i < n; i++) {
                    if (!item_placed[i]) continue;
                    int item_plate = arrangables[entries[i].orig_idx].bed_idx;
                    if (item_plate != pi) continue;
                    if (rot_bmps_all[i].empty()) continue;

                    double cur_rot = arrangables[entries[i].orig_idx].rotation;
                    const BitmapItem *cur_bmp = nullptr;
                    for (const auto &[bmp, rot] : rot_bmps_all[i]) {
                        if (std::abs(rot - cur_rot) < 0.01) { cur_bmp = &bmp; break; }
                    }
                    if (!cur_bmp || cur_bmp->width_px <= 0) continue;

                    coord_t tx = arrangables[entries[i].orig_idx].translation.x();
                    coord_t ty = arrangables[entries[i].orig_idx].translation.y();
                    int px = (int)((tx + cur_bmp->offset_x - effective_bed.min.x()) / res);
                    int py = (int)((ty + cur_bmp->offset_y - effective_bed.min.y()) / res);

                    if (collides(verify_bits, bed_w_words, bed_w_px, bed_h_px, *cur_bmp, px, py)) {
                        BOOST_LOG_TRIVIAL(error) << "BitmapArranger: OVERLAP detected for item "
                                                 << entries[i].orig_idx << " on plate " << pi;
                        arrangables[entries[i].orig_idx].bed_idx = -1;
                    } else {
                        stamp(verify_bits, bed_w_words, bed_w_px, bed_h_px, *cur_bmp, px, py);
                    }
                }
            }
        }
    }

    // --- Phase 7: Nothing left behind --- rescue orphaned items onto overflow plates
    // If multi-plate is enabled, every item MUST end up on a plate.
    // Any item with bed_idx = -1 gets placed on a new overflow plate
    // via the standard placement path.
    if (params.allow_multi_plate) {
        std::vector<int> orphans;
        for (size_t i = 0; i < arrangables.size(); i++) {
            if (arrangables[i].bed_idx < 0)
                orphans.push_back((int)i);
        }
        if (!orphans.empty()) {
            BOOST_LOG_TRIVIAL(warning) << "BitmapArranger: " << orphans.size()
                                       << " orphaned items — placing on overflow plates";
            for (int orig_idx : orphans) {
                int entry_idx = -1;
                for (int j = 0; j < n; j++) {
                    if ((int)entries[j].orig_idx == orig_idx) { entry_idx = j; break; }
                }
                if (entry_idx < 0) continue;

                bool placed = false;
                int n_existing = use_3d ? (int)plates_3d.size() : (int)plates.size();
                for (int pi = 0; pi < n_existing && !placed; pi++) {
                    if (use_3d && !rot_stacks_all[entry_idx].empty()) {
                        placed = try_place_on_plate_3d(entries[entry_idx], rot_stacks_all[entry_idx], pi);
                    } else if (!rot_bmps_all[entry_idx].empty()) {
                        placed = try_place_on_plate(entries[entry_idx], entry_idx, rot_bmps_all[entry_idx], pi);
                    }
                }

                if (!placed) {
                    int n_cur = use_3d ? (int)plates_3d.size() : (int)plates.size();
                    if (n_cur < MAX_PLATES) {
                        int new_pi = use_3d ? new_3d_plate() : new_2d_plate();
                        if (use_3d && !rot_stacks_all[entry_idx].empty()) {
                            placed = try_place_on_plate_3d(entries[entry_idx], rot_stacks_all[entry_idx], new_pi);
                        } else if (!rot_bmps_all[entry_idx].empty()) {
                            placed = try_place_on_plate(entries[entry_idx], entry_idx, rot_bmps_all[entry_idx], new_pi);
                        }
                    }
                }

                if (!placed) {
                    // MAX_PLATES reached or item physically too large — leave as unarranged.
                    BOOST_LOG_TRIVIAL(warning) << "BitmapArranger: item " << orig_idx
                                               << " — all plates full, leaving unarranged";
                }
            }
        }
    }

    // --- Post-placement: Cluster repositioning --- shift packed cluster to match user preference (center/corner)
    {
        int n_plates_shift = use_3d ? (int)plates_3d.size() : (int)plates.size();
        for (int pi = 0; pi < n_plates_shift; pi++) {
            coord_t min_left = effective_bed.max.x();
            coord_t min_bottom = effective_bed.max.y();
            coord_t max_right = effective_bed.min.x();
            coord_t max_top = effective_bed.min.y();
            bool has_items = false;

            for (size_t i = 0; i < arrangables.size(); i++) {
                if (arrangables[i].bed_idx != pi) continue;
                has_items = true;
                ExPolygon rotated = arrangables[i].poly;
                rotated.rotate(arrangables[i].rotation);
                BoundingBox bb = get_extents(rotated);
                coord_t left   = arrangables[i].translation.x() + bb.min.x();
                coord_t bottom = arrangables[i].translation.y() + bb.min.y();
                coord_t right  = arrangables[i].translation.x() + bb.max.x();
                coord_t top    = arrangables[i].translation.y() + bb.max.y();
                min_left   = std::min(min_left,   left);
                min_bottom = std::min(min_bottom, bottom);
                max_right  = std::max(max_right,  right);
                max_top    = std::max(max_top,    top);
            }
            if (!has_items) continue;

            coord_t shift_x = 0, shift_y = 0;
            if (params.placement_bias == 1) {
                shift_x = effective_bed.max.x() - max_right;
                shift_y = effective_bed.max.y() - max_top;
            } else {
                coord_t cluster_cx = (min_left + max_right) / 2;
                coord_t cluster_cy = (min_bottom + max_top) / 2;
                coord_t bed_cx = (effective_bed.min.x() + effective_bed.max.x()) / 2;
                coord_t bed_cy = (effective_bed.min.y() + effective_bed.max.y()) / 2;
                shift_x = bed_cx - cluster_cx;
                shift_y = bed_cy - cluster_cy;
                shift_x = std::max(shift_x, effective_bed.min.x() - min_left);
                shift_x = std::min(shift_x, effective_bed.max.x() - max_right);
                shift_y = std::max(shift_y, effective_bed.min.y() - min_bottom);
                shift_y = std::min(shift_y, effective_bed.max.y() - max_top);
            }

            for (auto &arr : arrangables) {
                if (arr.bed_idx == pi) {
                    arr.translation.x() += shift_x;
                    arr.translation.y() += shift_y;
                }
            }
        }
        ARRANGE_LOG("Placement shift: " << n_plates_shift << " plate(s), bias="
                    << params.placement_bias);
    }

    int total_plates = params.nesting_3d ? (int)plates_3d.size() : (int)plates.size();
    {
        auto t_now = std::chrono::high_resolution_clock::now();
        auto total_us = std::chrono::duration_cast<std::chrono::microseconds>(t_now - t_total_start).count();
        ARRANGE_LOG("DONE: " << placed_count << "/" << arrangables.size()
                    << " placed on " << total_plates << " plate(s)"
                    << ", total " << total_us << " us ("
                    << (total_us / 1000) << " ms)");
    }
}

}} // namespace Slic3r::arrangement
