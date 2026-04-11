#pragma once

// Shared invariant-checking utilities for BitmapNester tests.
// Include after ClipperUtils.hpp and BitmapNester.hpp.

#include "libslic3r/Arrange.hpp"
#include "libslic3r/BitmapNester.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"

#include <vector>
#include <set>
#include <cmath>
#include <algorithm>

namespace Slic3r { namespace arrangement { namespace test_utils {

// ─── Quality metrics adopted by the Three Seers vote 2026-04-11 ────
// See docs/PACKING_METRICS.md for the full deliberation. The seers
// rejected bbox density, convex hull compactness, adjacency, and
// libnest2d head-to-head. They unanimously adopted plate-count
// regression and determinism, plus Frodo's overflow piece count.
//
// Helpers below are intentionally minimal — they read only the
// observable output (items[].bed_idx) and never touch internal
// nester state.

// Highest plate index used by any placed item. Returns -1 if no
// item is placed (all UNARRANGED). For a single-plate fixture this
// returns 0; for a fixture that overflowed by one plate, returns 1.
inline int max_bed_idx(const ArrangePolygons &items)
{
    int m = -1;
    for (const auto &it : items)
        if (it.bed_idx != UNARRANGED && it.bed_idx > m) m = it.bed_idx;
    return m;
}

// Number of items placed on a plate other than plate 0. The "regret"
// metric: how many pieces wish they were on plate 0. Frodo's pick.
inline int overflow_piece_count(const ArrangePolygons &items)
{
    int n = 0;
    for (const auto &it : items)
        if (it.bed_idx != UNARRANGED && it.bed_idx > 0) ++n;
    return n;
}

// Run an arrange call twice on a fresh copy of the input items each
// time, then verify the two outputs are byte-identical on the
// observable fields (bed_idx, rotation, translation, itemid).
// Returns true if deterministic, false if any field diverges.
//
// Determinism is the cheapest possible defense against state-machine
// bugs in the consolidation pass — if a read-before-write or
// uninitialized cell sneaks in, the second run diverges. Cost: one
// extra arrange call per fixture.
template <class ArrangeFn>
inline bool deterministic_rerun(const ArrangePolygons &input,
                                ArrangeFn &&run)
{
    ArrangePolygons a = input;
    ArrangePolygons b = input;
    run(a);
    run(b);
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].bed_idx != b[i].bed_idx) return false;
        if (a[i].itemid != b[i].itemid) return false;
        if (std::abs(a[i].rotation - b[i].rotation) > 1e-9) return false;
        if (a[i].translation.x() != b[i].translation.x()) return false;
        if (a[i].translation.y() != b[i].translation.y()) return false;
    }
    return true;
}

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
