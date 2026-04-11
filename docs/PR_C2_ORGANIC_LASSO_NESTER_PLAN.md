# PR C2 — Organic Lasso Nester (implementation plan)

**Status:** Plan only. No code in this commit. Branch `feature/concave-bitmap-c2`
will be created from `feature/concave-bitmap` HEAD as the starting point.
**Date:** 2026-04-11
**Prompted by:** user directive "write a detailed plan for our new advanced
alternative, PR C2, new fork, does everything I wish this one did."

## Executive summary

PR C1 (`feature/concave-bitmap`) delivered a concave-aware bitmap nester
that beats upstream libnest2d by one plate on concave-heavy fixtures. Its
scoring primary is axis-aligned bbox area growth, its tall-parts-centered
behavior is a weak tiebreaker, and its consolidation + placement are
tangled into a single 900-line lambda.

PR C2 is the clean-slate follow-up that does everything C1 compromised on:

1. **Organic lasso scoring.** Primary metric is the convex hull perimeter
   of the cluster, not the bbox. "Boxes lie" — two layouts with the same
   bbox have very different hull perimeters and very different travel
   costs. The hull is the truth.
2. **Tall parts in the middle as a primary placement driver.** Not a
   tiebreaker. Tall items are PLACED at target positions near the cluster
   geometric center from the start, not retroactively shifted there.
3. **Two-phase architecture.** Phase 1: nest parts into a tight island
   without bed constraint. Phase 2: locate the island on the plate(s).
   Decouples packing shape from bed shape, which is conceptually cleaner
   and opens the door to better splitting.
4. **Proactive multi-plate splitting.** If K_min > 1 for the input, split
   items into roughly equal sub-clusters BEFORE packing each. Don't pack
   plate 1 to maximum density then struggle with plate 2's leftovers.
5. **Spillover recovery.** If after phase 1 an item couldn't fit in any
   island, try harder — rotate, reshuffle, split further.
6. **Reuse C1's primitives.** Rasterization, dilation, collision,
   silhouette extraction, test infrastructure. The low-level bitmap
   machinery is solid and portable; C2 just builds a new algorithm on top.

PR C2 is **not a replacement for PR C1**. C1 should ship first because it's
a real improvement over upstream today, and the fork gives us a clean
implementation surface for experimenting with the more invasive ideas.

---

## Background — why C1's compromises forced C2

The 2026-04-11 autopilot session landed C1's final form. Along the way it
tried and measured several ideas that hit limits:

### Why bbox scoring hits a wall

- **4 equal 80x80 squares on 200x200 bed, old hybrid anchor:** item 2
  landed at (60,80) instead of (0,80) or (80,80) because the bed-center
  tiebreaker preferred x=60. Fixed by reverting to pure corner anchor,
  but the underlying issue — bbox-area ties are common and the anchor
  tiebreaker is arbitrary — remains.
- **Diag² as primary:** tried as a direct worst-case-travel proxy,
  shifted positions enough to break a pre-existing test and got reverted.
  Would need a tuple refactor + test loosening that aren't worth the
  risk in C1.
- **Perimeter as lex secondary:** landed successfully in C1, but it's
  still bbox perimeter, not hull perimeter. It only helps when bbox
  areas tie, which misses the more common case where areas differ
  slightly but hulls differ more.

The right metric is **convex hull perimeter of the cluster after
adding the candidate item**. Expensive to compute per candidate — each
placement attempt recomputes the cluster hull from scratch — but it
captures travel cost directly and doesn't have the "boxes lie" failure
mode.

### Why tall-parts-centered as a tiebreaker doesn't work

- The scoring bias is a tertiary tiebreaker that kicks in only when
  area + perimeter both tie. In a 6-item test (2 tall, 4 short, equal
  footprints) the tall items BL-filled to cluster corners because FFD
  placed them first and their first-placement choice is deterministic
  by scan order, not by bias.
- Height-weighted post-centering was attempted in the autopilot session
  but bounds-clamping defeated it for clusters that fill the plate
  width. Cluster can't shift beyond bed, so shifting to put tall COM
  at bed center doesn't happen.
- **The real fix is placement-time, not centering-time.** Tall items
  need to be placed AT cluster center, not at cluster edges. That
  requires inverting the current FFD-first-then-BL-fill ordering into
  something that reserves central space for tall items.

### Why the 900-line do_one_pass lambda is a structural problem

- C1's placement + consolidation + post-centering are all inside a
  single lambda captured by reference. The smart-shuffle restart driver
  needs to reset state and replay the whole pipeline, which forced the
  monolithic structure.
- Small-batch lookahead (tried and abandoned) would require extracting
  "try placing item X against state Y" into a reusable helper. That's
  impossible while the placement body is 400 lines deep inside a
  capture-by-ref lambda.
- C2 starts with the extraction: `class Nester` with methods for each
  phase, explicit state objects, and the restart driver as an outer
  controller.

---

## C2 goals (restated for the plan)

From the user's directive:

1. **Tall parts in the middle.**
2. **Minimum circumference scoring using organic lasso (convex hull
   perimeter) of nested parts.** No bbox math. Boxes lie.
3. **Nest parts into the smallest island.** Packing phase has no bed
   constraint — optimize for cluster tightness, not for "does it fit
   plate 0".
4. **Tall in the middle, then fit it on the plate.** Locate the island
   on the plate as a post-process.
5. **If anything falls off, find it a home.** Spillover recovery.
6. **If two plates obviously needed, split into equal islands first.**
   Proactive K_min-aware splitting before packing.

---

## Architecture

### Phase structure

```
NesterC2::arrange(items, excludes, bed, params) {
    // 0. Preparation
    Cache cache = build_cache(items);
        // per-item: silhouette, rotations, bbox, area, perimeter,
        // convex hull, Z height, "hardness" (convex hull / area ratio)

    // 1. Target plate count estimation
    int k_min = estimate_min_plates(cache, bed);
        // ceil(sum_silhouette_area / bed_area × 0.85 packing_factor)

    // 2. Partition items into target-plate groups (K-means-style)
    std::vector<ItemGroup> groups = partition_items(cache, k_min);
        // each group targets roughly equal total area
        // groups try to balance item-type mix, rotation diversity
        // and tall-item distribution

    // 3. For each group, pack it as a tight organic island
    std::vector<Island> islands;
    for (auto& group : groups) {
        Island island = pack_as_island(cache, group);
            // bound-agnostic placement optimizing hull perimeter
            // tall items placed at island center FIRST
            // short items placed on the island perimeter last
            // returns an unlocated island (items have relative positions)
        islands.push_back(island);
    }

    // 4. Locate each island on a plate
    for (size_t i = 0; i < islands.size(); ++i) {
        Plate& plate = allocate_plate(i);
        locate_island_on_plate(islands[i], plate);
            // centers the island on the plate
            // bounds-checks against excludes
            // falls back to center-of-bed if island is oversized
    }

    // 5. Spillover recovery
    recover_spillover(items, islands, bed);
        // for any item that couldn't fit its island, try:
        //   (a) a different rotation
        //   (b) migration to a neighbor island (if fits)
        //   (c) a new island (new plate)

    // 6. Finalize
    write_back(items, islands);
}
```

### Phase 1 — Cache building

Builds a per-item bundle reusable across all phases:

```cpp
struct ItemInfo {
    size_t original_idx;       // index in input items[]
    ExPolygons silhouette;     // concave silhouette(s)
    std::vector<RasterCache> rot_cache;  // one per allowed rotation

    // Geometric summary
    double silhouette_area_mm2;
    double silhouette_perimeter_mm;
    Polygon convex_hull;
    double hull_area_mm2;
    double hull_perimeter_mm;
    double max_dim_mm;         // longest bbox side across rotations

    // 3D info
    double height_mm;          // Z extent from caller

    // Derived heuristics
    double hardness_score;     // convex_hull_area / silhouette_area
                               // higher = harder to pack around
    double priority_score;     // composite: priority + hardness + height
};

struct RasterCache {
    double rotation_rad;
    std::vector<uint64_t> silhouette_bitmap;
    std::vector<uint64_t> hull_bitmap;
    int w_px, h_px, iwpr;
    BoundingBox bb;
};
```

### Phase 2 — Minimum plate estimation

Simple heuristic with a packing efficiency factor:

```cpp
int estimate_min_plates(const Cache& cache, BoundingBox bed) {
    double total_silhouette = 0;
    for (auto& info : cache) total_silhouette += info.silhouette_area_mm2;
    double bed_area = bed.area_mm2();
    double packing_efficiency = 0.82;  // conservative
    return (int)std::ceil(total_silhouette / (bed_area * packing_efficiency));
}
```

The `0.82` factor is a fudge based on empirical packing densities of
concave shapes. Tuneable per-input type. For tetrominoes the achievable
density is ~87% so this underestimates; for real STL mixes ~70% so it
slightly overestimates but that's safer.

### Phase 3 — Item partitioning

Split N items into K groups where each group's total silhouette area
is roughly `total / K`. Also balance:

- Item-type distribution (a rough similarity hash on silhouette shape)
- Rotation diversity (equal numbers of each rotation)
- Tall-item distribution (each group should have ~equal count of tall items)

Algorithm sketch:

```cpp
std::vector<ItemGroup> partition_items(Cache& cache, int k) {
    // Sort items by a priority score combining area, hardness, height
    sort_by_priority(cache);

    // Round-robin deal into K buckets, accounting for each bucket's
    // running area to keep them balanced
    std::vector<ItemGroup> groups(k);
    std::vector<double> bucket_area(k, 0.0);
    for (auto& info : cache) {
        // Find the bucket with the smallest current area
        int target = std::min_element(bucket_area.begin(),
                                       bucket_area.end()) - bucket_area.begin();
        groups[target].items.push_back(info.original_idx);
        bucket_area[target] += info.silhouette_area_mm2;
    }

    // Second pass: fix tall-item distribution if unbalanced
    rebalance_tall_items(groups, cache);

    return groups;
}
```

This is a **global** partitioning step, run once at the start. It's
fundamentally different from C1's per-item greedy placement: C2
decides ahead of time which items go together.

### Phase 4 — Island packing

The heart of C2. Takes a group of items and produces a tight cluster
with no bed constraint. This is where "organic lasso" scoring lives.

```cpp
Island pack_as_island(Cache& cache, const ItemGroup& group) {
    Island island;

    // Step 1: Place the tallest item at the origin as the seed
    const ItemInfo& seed = find_tallest(group, cache);
    island.add(seed, Vec2d(0, 0), seed.preferred_rotation());

    // Step 2: Place remaining tall items around the seed, biased toward
    // the current cluster center (forming the island's tall core)
    for (const auto& info : tall_items_in_group(group, cache)) {
        if (info == seed) continue;
        Placement p = find_best_placement(info, island,
                                          PlaceStrategy::CenterBias);
        island.add(info, p.pos, p.rot);
    }

    // Step 3: Place short items around the tall core
    for (const auto& info : short_items_in_group(group, cache)) {
        Placement p = find_best_placement(info, island,
                                          PlaceStrategy::HullBoundary);
        island.add(info, p.pos, p.rot);
    }

    // Step 4: Compact the island via iterative pair-swap
    // "Can I swap any two items to reduce hull perimeter?"
    compact_by_swap(island);

    return island;
}
```

The key method: `find_best_placement`.

```cpp
Placement find_best_placement(const ItemInfo& item, const Island& island,
                              PlaceStrategy strategy) {
    // Candidates: all positions where the item could sit without
    // colliding with any existing island item. Computed via NFP-like
    // bitmap scan against the island's composite bitmap.
    std::vector<Placement> candidates = enumerate_clear_positions(item, island);

    // Score each candidate by:
    //   primary:   delta(convex hull perimeter of island + item)
    //   secondary: for CenterBias: distance from candidate center
    //              to current cluster centroid. For HullBoundary:
    //              inverse — prefer positions NEAR the hull edge.
    //   tertiary:  rotation preference (item.preferred_rotation)
    for (auto& cand : candidates) {
        Polygon new_hull = compute_hull_with(island.hull, item, cand);
        cand.delta_hull_perimeter = new_hull.perimeter() - island.hull_perimeter;
        cand.secondary = compute_secondary(cand, island, strategy);
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const Placement& a, const Placement& b) {
                  if (a.delta_hull_perimeter != b.delta_hull_perimeter)
                      return a.delta_hull_perimeter < b.delta_hull_perimeter;
                  return a.secondary < b.secondary;
              });

    return candidates.front();
}
```

### Phase 5 — Island placement on bed

```cpp
void locate_island_on_plate(const Island& island, Plate& plate) {
    // Center the island by its hull centroid at the bed center
    Vec2d centroid = island.hull.centroid();
    Vec2d bed_center = plate.bed.center();
    Vec2d shift = bed_center - centroid;

    // Clamp to keep island inside bed
    shift = clamp_to_fit(island.bbox() + shift, plate.bed);

    // Apply shift to every item in the island
    for (auto& item : island.items)
        item.translation += shift;

    plate.items.append(island.items);
}
```

If the island is too big for the plate, this is where the "will obviously
be 2 plates" check belongs — before we get here, phase 2 should have
estimated enough plates. If we still overflow, phase 6 (spillover
recovery) handles it.

### Phase 6 — Spillover recovery

For any item that couldn't fit in its target island, try:

1. **Rotation retry.** Try every allowed rotation on the same island.
2. **Migration to a neighbor island.** Maybe a less-full island can
   absorb it.
3. **New island allocation.** Spin up a new plate's worth.

```cpp
void recover_spillover(ArrangePolygons& items,
                       std::vector<Island>& islands, BoundingBox bed) {
    std::vector<size_t> stragglers = find_unfit_items(items);
    for (size_t idx : stragglers) {
        if (try_any_rotation(items[idx], islands[target_island(idx)])) continue;
        if (try_neighbor_island(items[idx], islands)) continue;
        islands.push_back(Island::from_single_item(items[idx]));
    }
}
```

---

## What we reuse from C1

- `rasterize`, `dilate_bitmap`, `collides`, `stamp` — the bitmap
  primitives. Fast, correct, well-tested.
- Silhouette extraction via `item.concave_regions` or `item.poly` fallback
- Exclude zone handling (dual-bitmap collision check)
- `ArrangePolygon` data model (translation, rotation, bed_idx, height)
- Test infrastructure: `max_bed_idx`, `overflow_piece_count`, `no_overlap`,
  `deterministic_rerun`, `cluster_bbox_perimeter_mm`, `cluster_compactness`
- `BITMAP_NESTER_TESTING` ifdef guard pattern for exposing internals

## What C2 replaces

- **Scoring primary:** bbox area → convex hull perimeter
- **Placement order:** FFD corner-seeded → tall-seeded center-out
- **Structure:** 900-line monolithic lambda → phase-based class with
  explicit state objects
- **Consolidation:** backwards migration → proactive spillover recovery
- **Post-centering:** bbox-center-to-bed-center → island-hull-centroid-
  to-bed-center
- **Plate count:** reactive overflow → proactive estimation + partitioning

## What C2 explicitly does NOT do

- **No NFP math.** C2 stays bitmap-based for collision detection. The
  hull-perimeter scoring is a per-candidate addition, not a replacement
  of the scan primitive.
- **No FFT.** That's a separate future fork (`feature/concave-bitmap-fft`
  already has the design doc). C2 can integrate FFT later if needed
  but isn't gated on it.
- **No slicer ordering awareness.** C2 still outputs a layout, not a
  print-time plan. Travel-cost scoring is a proxy, not the real
  measurement.
- **No small-batch permutation search.** C2's proactive partitioning
  is conceptually adjacent to small-batch but simpler — split, pack,
  recover. No per-batch permutation loops.

---

## Implementation phases with milestones

### M0 — branch and skeleton (1 day)

- Create `feature/concave-bitmap-c2` from current HEAD
- Copy `BitmapNester.hpp` to `BitmapNesterC2.hpp` (keep C1 unaffected)
- Add a new entry point `BitmapNesterC2::arrange` with the phase
  stubs (empty functions that return the input untouched)
- Add a new test target that exercises C2 separately from C1
- Milestone: C2 compiles, tests for it pass trivially, existing C1
  tests unaffected

### M1 — Cache building and partitioning (2-3 days)

- Implement `build_cache` — reuse C1's rot_cache building
- Implement `estimate_min_plates`, `partition_items`
- Add tests: 4 items on a small bed with known optimal partition,
  N items with balanced tall distribution, edge cases (zero items,
  one oversized item)
- Milestone: phases 0-2 work correctly; items get divided into groups

### M2 — Island packing with hull scoring (4-6 days)

- Implement `pack_as_island` with the tall-seeded center-out
  strategy
- Implement `find_best_placement` with hull-perimeter scoring
- Implement hull recomputation incrementally (don't recompute from
  scratch every candidate — update the hull as items are added)
- Implement `compact_by_swap` iterative refinement
- Add tests: 4 equal squares pack as 2x2, tall-parts-in-middle
  assertion with positional verification, realistic OBJ mix
- Milestone: island packing produces tight clusters; hull-perimeter
  scoring visibly better than bbox on same inputs

### M3 — Island placement + spillover (1-2 days)

- Implement `locate_island_on_plate`
- Implement `recover_spillover`
- Add tests: island-too-big-for-plate, stragglers that fit in
  neighbors, stragglers that need new plates
- Milestone: full pipeline end-to-end, beats C1 on matched fixtures

### M4 — Performance and polish (2-3 days)

- Profile hull recomputation cost; optimize if needed
- Profile partitioning for large N; optimize if needed
- Split `BitmapNesterC2.hpp` into `.hpp` + `.cpp`
- Add docs: design notes, algorithm reference, example outputs
- Milestone: performance acceptable (within 2x of C1), code
  structured for upstream review

### M5 — Comparison benchmarks and PR prep (2 days)

- Benchmark suite: C1 vs C2 on 10+ fixtures (synthetic + real STLs)
- Visual comparison screenshots: upstream, C1, C2 on same inputs
- Update `docs/PR_C2_ORGANIC_LASSO_NESTER_PLAN.md` with measured
  results
- Write PR C2 research report (parallel to PR_C1)
- Milestone: PR C2 ready for user review and upstream consideration

Total estimate: **12-17 working days** for a full implementation.
Could shrink to 6-8 days for an MVP that cuts scope (no M4 perf
optimization, no M5 benchmark parity).

---

## Risks and mitigations

### Risk: Hull recomputation dominates runtime

Computing a convex hull of N points is O(N log N). Doing it per
candidate placement is O(N × candidates × log N). For N=125 and
1000 candidates that's ~2 million operations per item, 250 million
for the whole pack. May be too slow.

**Mitigation 1:** incremental hull update. When adding a point to a
hull, only recompute the affected edges — O(log N) amortized. This
is a known technique (Graham scan or Andrew monotone chain with
local updates).

**Mitigation 2:** approximate hull with a bounding polygon of fixed
vertex count (e.g., 16-sided). Trades accuracy for speed. Good
first cut; refine to true hull if quality is insufficient.

**Mitigation 3:** FFT-based bitmap correlation for the scan (separate
future fork idea). Not required for C2 but available if perf blocks
us.

### Risk: Partitioning produces unbalanced groups that can't pack

Partitioning is done WITHOUT placement feedback. A group might have
items that can't fit together even though its total area suggests
they should.

**Mitigation:** the spillover recovery phase handles this. If a
group's items don't all fit, the stragglers migrate. In the worst
case, we end up with more plates than K_min estimated. Graceful
degradation.

### Risk: Tall-parts-in-middle is harder than it looks

The naive "place tallest at origin, other tall around it, short at
periphery" doesn't guarantee the tall items end up at the ISLAND
CENTER. The island center is emergent from all item positions.

**Mitigation 1:** after placing all items, compute the island centroid
and shift all positions so tall items' weighted centroid matches
the island centroid. This is a pure translation — doesn't change
the pack quality.

**Mitigation 2:** explicit tall-item pre-allocation. Reserve a central
region sized for the tall items first, then pack short items around
it. More complex but guaranteed tall-centered.

### Risk: C2 underperforms C1 on cases C1 handles well

C1 is production-tested on the user's real fixtures and delivers
visible wins. If C2 regresses on any of those, we've made the fork
worse without justification.

**Mitigation:** keep C1 in the tree alongside C2. Let users toggle
via a settings flag. Ship both, measure both, merge the winner into
mainline. C2 has to beat C1 on matched benchmarks before it
replaces it.

### Risk: Organic lasso scoring has pathologies we haven't thought of

Hull perimeter is a direct travel proxy but edge cases exist: a
single degenerate item could have a huge hull relative to its area,
a shape with many small disconnected islands has a weird hull,
rotation can drastically change hull perimeter for elongated parts.

**Mitigation:** adversarial test suite that tries edge cases early.
If hull scoring hits a pathology, fall back to a lexicographic
(hull perim, bbox perim, bbox area) tuple so the weaker metrics
handle degenerate cases.

### Risk: 12-17 days is a long scope for a single fork

The user might want something visible sooner.

**Mitigation:** the phase plan above has natural stopping points.
M0-M2 alone is 6-9 days and delivers a working (if rough) C2. M3-M5
are polish and benchmarking. Ship M2 as a preview and iterate.

---

## Test strategy

### Reuse from C1

All C1 tests stay active on both branches. The test suite is the
regression gate — anything that worked on C1 must work on C2.

### New C2-specific tests

- **Hull perimeter regression assertions.** For each scenario, after
  arranging, compute `cluster_hull_perimeter_mm` (new helper) and
  assert it's ≤ some floor. Floors tighten as the algorithm improves.
- **Tall-parts-in-middle positional test.** 2 tall + 4 short on a
  large bed (enough slack so bounds clamping doesn't interfere),
  assert tall items' centroid is within 10mm of bed center.
- **Partitioning balance test.** Given N items, check that partition
  produces groups within X% of each other's total area.
- **Spillover recovery test.** Force a group that can't fit, verify
  stragglers find homes.
- **Proactive multi-plate test.** Input that obviously needs 2 plates
  at K_min, verify phase 2 estimates 2 and partitions before packing
  rather than reactive overflow.
- **C1 vs C2 benchmark test.** Runs both arrangers on a fixed
  fixture, compares hull perimeters and plate counts, documents
  the delta.

### Benchmark corpus

Fixed set of STL files used for A/B comparison:
- `tetris_plate_240x240_balanced.stl` — 125 tetrominoes
- User's real L-bracket sets (25 and 80)
- A mix of 20-30 real user parts across size classes
- An adversarial case: 10 identical elongated bars that C1 handles
  poorly
- A curved-part case: 5 gears or cogs with minor concave features
  (tests 45° rotation if extended)

Each benchmark produces: plate count, total hull perimeter, total
travel-estimate (sum of inter-part centroids in a nearest-neighbor
tour), runtime, screenshot.

---

## Success metrics for PR C2

Hard targets (must hit to land):

1. Matches or beats C1 on every matched-fixture benchmark
2. Produces visibly tall-centered layouts on any real STL with mixed
   heights
3. Tests pass with the same 130+ scenarios from C1 (no regressions)
4. Runtime within 2x of C1 for N ≤ 125

Soft targets (would be nice):

5. Beats C1 by ≥1 plate on the 125-tetromino fixture
6. Visibly smaller hull perimeter on realistic 20-30-part plates
7. Code organized for an upstream reviewer without the 1700-line
   header complaint
8. Usable as a drop-in `NesterC2::arrange` via a settings flag, so
   A/B testing is possible in production

---

## Dependencies from PR C1

- `f02a0ecdb8` (C1 tip) is the base for C2's branch
- `bitmap_test_utils.hpp` (shared test helpers)
- `BitmapNester.hpp` low-level primitives (rasterize, dilate, etc)
- `docs/PACKING_METRICS.md` — Three Seers vote, metric philosophy
- `docs/PR_C1_CONCAVE_BITMAP_NESTER.md` — the full C1 narrative, which
  C2 needs to reference for "what we inherit"
- `docs/FFT_NESTER_IDEA.md` — parallel future-fork speculation

## Related memories and notes

- `memory/feedback_arrange_priorities.md` — travel first, plate count
  second, tall centered third (this plan realizes the priority order)
- `memory/feedback_organic_lassos_not_bboxes.md` — the exact user
  framing that motivated this plan
- `memory/project_bitmap_nester_value_prop.md` — one-plate-better
  framing for upstream narrative (C2 should aim for a larger delta
  on realistic cases)
- `dev/sessions/2026-04-11/state-of-the-union.md` — the
  "in a perfect world" section is the template for C2's architecture
- `dev/sessions/2026-04-11/design-gravitational-stack.md` — the
  two-phase pack-then-locate idea formalized in this plan

---

## Next action (when the user greenlights)

1. Create branch: `git branch feature/concave-bitmap-c2 HEAD`
2. Commit this plan doc to both branches
3. M0 skeleton: stub `BitmapNesterC2.hpp` with phase entry points
4. First test: "C2 arrange() on zero items returns no plates" — the
   trivial green bar that proves the pipeline is alive
5. Iterate through M1-M5 as time allows

This plan is the scaffolding. Implementation is deferred until the
user reviews it and says go.
