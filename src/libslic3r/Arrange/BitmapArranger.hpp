#ifndef BITMAP_ARRANGER_HPP
#define BITMAP_ARRANGER_HPP

#include "libslic3r/Arrange.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ClipperUtils.hpp"

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <functional>
#include <optional>

namespace Slic3r { namespace arrangement {

// Bitmap-based arranger for concave polygon packing.
// Rasterizes part silhouettes into bitmaps and uses bottom-left-fill
// placement with bitwise collision detection.
class BitmapArranger {
public:
    // Default resolution. Overridden by ArrangeParams::bitmap_resolution_mm.
    static constexpr double DEFAULT_RESOLUTION_MM = 0.5;

    // Bitmap encoding:
    //   bits[py * width_words + (px / 64)], bit (px % 64) = 1 means occupied.
    //   Bit 0 (LSB) = leftmost pixel in each 64-pixel word.
    //   Row 0 = bottom edge (Y-min). Row-major, bottom-to-top.
    //   offset_x/y: scaled coordinate of pixel (0,0) = bounding box min of the rasterized polygon.
    struct BitmapItem {
        std::vector<uint64_t> bits;   // row-major, 64 pixels per word
        int width_words = 0;          // width in 64-bit words
        int width_px    = 0;          // width in pixels
        int height_px   = 0;          // height in pixels
        coord_t offset_x = 0;        // origin offset from bounding box min
        coord_t offset_y = 0;
    };

    // Run the bitmap arrangement.
    // Modifies arrangables in-place (translation, rotation, bed_idx).
    static void arrange(
        ArrangePolygons &arrangables,
        const ArrangePolygons &excludes,
        const BoundingBox &bed,
        const ArrangeParams &params);

private:
    // Rasterize an ExPolygon into a bitmap (for excludes and non-concave items)
    static BitmapItem rasterize(const ExPolygon &poly, coord_t inflation, coord_t res);

    // Rasterize projected 2D triangles directly into a bitmap (concave mode).
    // Every 3 consecutive Points form one triangle. Each triangle is scanline-filled
    // into the bitmap — same as GPU rasterization. No Clipper, O(triangles × pixels_per_tri).
    static BitmapItem rasterize_triangles(const Points &tri_verts, coord_t inflation, coord_t res);

    // Rotate a bitmap by an arbitrary angle. Always rotates from the source bitmap
    // (not iteratively) to avoid cumulative quality loss. Uses reverse-mapping with
    // nearest-neighbor sampling. The rotated bitmap may be larger than the source.
    static BitmapItem rotate_bitmap(const BitmapItem &src, double angle_rad, coord_t res);

    // Rasterize onto the bed bitmap (stamp an item at position).
    // bed_w = width in 64-bit words, bed_w_px = width in pixels (for bounds).
    static void stamp(std::vector<uint64_t> &bed_bits,
                      int bed_w, int bed_w_px, int bed_h,
                      const BitmapItem &item, int ox, int oy);

    // Check if placing item at (ox, oy) on bed collides.
    static bool collides(const std::vector<uint64_t> &bed_bits,
                         int bed_w, int bed_w_px, int bed_h,
                         const BitmapItem &item, int ox, int oy);

    // Find best placement for an item on the bed using bottom-left-fill.
    static std::optional<std::pair<int,int>> find_placement(
        const std::vector<uint64_t> &bed_bits,
        int bed_w, int bed_w_px, int bed_h,
        const BitmapItem &item,
        int step);
};

}} // namespace Slic3r::arrangement

#endif // BITMAP_ARRANGER_HPP
