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

#include <vector>
#include <cstddef>

namespace Slic3r { namespace arrangement {

// ---------------------------------------------------------------------------
// Per-item information cached once per arrange call. Reused across every
// phase of the pipeline so that per-item work (rasterization, hull, bbox)
// only happens once. Grows as the phases come online — M0 holds the
// minimum necessary for the pass-through.
// ---------------------------------------------------------------------------

struct NesterC2ItemInfo {
    std::size_t original_idx = 0;

    // Geometric summary. Populated in build_cache (M1). Zero-initialized
    // in M0 so pass-through compiles.
    double silhouette_area_mm2   = 0.0;
    double silhouette_perimeter_mm = 0.0;
    double hull_area_mm2         = 0.0;
    double hull_perimeter_mm     = 0.0;
    double max_dim_mm            = 0.0;
    double height_mm             = 0.0;

    // Derived heuristics. Populated in build_cache (M1).
    // hardness_score = hull_area / silhouette_area; higher = harder to pack
    // priority_score = composite of priority, hardness, and height
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
        // M0: delegate. M1+ will replace this with the phase pipeline:
        //   Cache cache = build_cache(items, params);
        //   int k = estimate_min_plates(cache, bed);
        //   auto groups = partition_items(cache, k);
        //   std::vector<Island> islands;
        //   for (auto& g : groups) islands.push_back(pack_as_island(cache, g));
        //   for (size_t i = 0; i < islands.size(); ++i)
        //       locate_island_on_plate(islands[i], bed, i, params);
        //   recover_spillover(items, islands, bed, params);
        //   write_back(items, cache, islands);
        BitmapNester::arrange(items, excludes, bed, params);
    }

private:
    // ----- Phase stubs. Signatures will stabilize in M1; bodies come
    // online per the milestone schedule in PR_C2_ORGANIC_LASSO_NESTER_PLAN.md.

    // M1: build per-item geometric + heuristic cache.
    // static std::vector<NesterC2ItemInfo>
    //     build_cache(const arrangement::ArrangePolygons& items,
    //                 const arrangement::ArrangeParams& params);

    // M1: estimate minimum plate count from total silhouette area + bed area
    // and a packing-efficiency constant (tuned empirically).
    // static int estimate_min_plates(const std::vector<NesterC2ItemInfo>& cache,
    //                                const BoundingBox& bed);

    // M1: partition items into K groups with balanced area + tall distribution.
    // static std::vector<NesterC2ItemGroup>
    //     partition_items(const std::vector<NesterC2ItemInfo>& cache, int k);

    // M2: bound-agnostic packing into a tight island. Hull-perimeter primary
    // scoring. Tall items placed at the seed; short items added at the hull
    // boundary. Compacted via iterative pair-swap refinement.
    // static NesterC2Island
    //     pack_as_island(const std::vector<NesterC2ItemInfo>& cache,
    //                    const NesterC2ItemGroup& group,
    //                    const arrangement::ArrangeParams& params);

    // M3: locate a packed island on a plate, clamping to bed bounds. Uses
    // the island's hull centroid as the centering target.
    // static void locate_island_on_plate(NesterC2Island& island,
    //                                    const BoundingBox& bed,
    //                                    int plate_idx,
    //                                    const arrangement::ArrangeParams& params);

    // M3: attempt to find a home for items that couldn't fit their target
    // island. Retry rotations, migrate to neighbor islands, or spawn a new
    // island as a last resort.
    // static void recover_spillover(arrangement::ArrangePolygons& items,
    //                               std::vector<NesterC2Island>& islands,
    //                               const BoundingBox& bed,
    //                               const arrangement::ArrangeParams& params);
};

}} // namespace Slic3r::arrangement
