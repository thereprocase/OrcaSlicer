#pragma once

// Shared invariant-checking utilities for BitmapNester tests.
// Include after ClipperUtils.hpp and BitmapNester.hpp.

#include "libslic3r/Arrange.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"

#include <vector>
#include <set>
#include <cmath>
#include <algorithm>

namespace Slic3r { namespace arrangement { namespace test_utils {

// Check every placed item's transformed polygon lies inside the bed bounds.
// Returns false if any placed item protrudes outside the bed.
inline bool all_within_bounds(const ArrangePolygons &items, const BoundingBox &bed)
{
    for (auto &it : items) {
        if (it.bed_idx == UNARRANGED) continue;
        ExPolygon placed = it.transformed_poly();
        BoundingBox bb = get_extents(placed);
        // Allow 1mm tolerance for raster quantization
        coord_t tol = scaled<coord_t>(1.0);
        if (bb.min.x() < bed.min.x() - tol) return false;
        if (bb.min.y() < bed.min.y() - tol) return false;
        if (bb.max.x() > bed.max.x() + tol) return false;
        if (bb.max.y() > bed.max.y() + tol) return false;
    }
    return true;
}

// Check bed_idx values are contiguous: 0,1,2,... with no gaps.
// Ignores UNARRANGED items. Returns false on gaps like {0, 2}.
inline bool no_bed_gaps(const ArrangePolygons &items)
{
    std::set<int> beds;
    for (auto &it : items) {
        if (it.bed_idx != UNARRANGED)
            beds.insert(it.bed_idx);
    }
    if (beds.empty()) return true;
    int expected = 0;
    for (int idx : beds) {
        if (idx != expected) return false;
        ++expected;
    }
    return true;
}

// Check no item is UNARRANGED. Returns false if any item was not placed.
inline bool all_placed(const ArrangePolygons &items)
{
    for (auto &it : items) {
        if (it.bed_idx == UNARRANGED) return false;
    }
    return true;
}

// Check no two placed items on the same bed overlap.
// Allows a tiny sliver (< 0.01 mm^2) for numerical noise.
inline bool no_overlap(const ArrangePolygons &items)
{
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].bed_idx == UNARRANGED) continue;
        for (size_t j = i + 1; j < items.size(); ++j) {
            if (items[j].bed_idx != items[i].bed_idx) continue;

            ExPolygon pi = items[i].transformed_poly();
            ExPolygon pj = items[j].transformed_poly();

            ExPolygons inter = intersection_ex(ExPolygons{pi}, ExPolygons{pj});
            for (auto &seg : inter) {
                if (std::abs(seg.area()) > scaled<double>(0.01) * scaled<double>(0.01))
                    return false;
            }
        }
    }
    return true;
}

}}} // namespace Slic3r::arrangement::test_utils
