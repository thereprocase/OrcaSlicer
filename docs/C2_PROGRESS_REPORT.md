# C2 Organic Lasso Nester — Progress Report

**Branch:** `feature/concave-bitmap-c2`
**Date:** 2026-04-12
**HEAD:** `d986d65bec`
**Tests:** 22 cases / 255 assertions — all green
**DLL Snapshots:** 2 (m4.5-trash-compactor + full-pipeline-adaptive-res building)

---

## Where We've Been

The C1 bitmap nester -- `BitmapNester` -- proved the core thesis: raster-scan placement on a pixel grid can beat OrcaSlicer's default libnest2d arranger by exploiting concave geometry instead of discarding it behind convex hulls. C1 shipped with 71 green tests and a headline result of 2 plates vs 3 on the tetris125 + 80-L-bracket fixture. It works. The question was whether it could work *well enough* to ship as a real alternative.

The answer was no, not without architectural surgery. C1's scoring function evolved through a series of increasingly uncomfortable patches -- squared-diagonal primary scores, perimeter-delta secondaries, desperation stride retries, tall-parts-centered bias hacks. Each one fixed a specific fixture but introduced coupling between the scoring logic and the placement search. The equal-squares gap (`fe50b2e8`) demonstrated the core problem: bbox-derived metrics cannot distinguish a crescent from a square of the same extent, so the scorer makes degenerate choices on shapes with high concavity. The revert at `ca61b17` -- backing out squared-diagonal scoring within hours of landing it -- was the clearest signal that C1's scoring surface had become unnavigable.

C2 forked from C1 at `ff658fe` with a detailed plan document (`docs/PR_C2_ORGANIC_LASSO_NESTER_PLAN.md`) and a strict constraint: hull-based metrics everywhere, bbox nowhere. The M0 skeleton was a pass-through that delegated every call to C1, keeping all 71 tests green while the new pipeline grew behind it. This was deliberate -- the branch never went red during construction.

The milestone sequence tells the story. M1 built the per-item cache and partitioner. M2.0-M2.1 introduced a winding-number point-in-polygon oracle and analytical test battery -- the harness validation work that the project rules explicitly prioritize over solver changes. M2.3 landed the real hull-perimeter-scored greedy loop. M2.5 was a four-commit bbox purge: bitmap AND collision replaced Clipper `intersection_ex`, silhouette area replaced max_dim_mm as the primary sort key, and the friend-class pattern gave C2 access to C1's rasterization and stamping primitives without duplicating them. M2.6 added height-descending sort for passive tall-parts-centering via Council of Elrond + Sauron arbitration. M3.0-M3.1 handled plate placement and spillover. M4.5 introduced the trash compactor -- a post-placement bitmap compaction pass that pushes parts inward from bed walls using XOR unstamp/restamp, also designed through Council + Sauron arbitration with user corrections on phase placement and force model.

The most consequential design decision was the bbox purge directive. Nine points, systematically executed across M2.5.1 through M2.5.4. Every place C1 used a bounding box -- scan region, sort key, pre-filter, collision detection -- C2 replaced with either hull-derived or silhouette-derived metrics. The one bbox that remains (`island.bbox`) is explicitly documented as "min/max of actual placed geometry," computed from hull vertices, not from the polygon's axis-aligned extent. This is the architectural spine of C2.

---

## Current Architecture

`BitmapNesterC2` lives entirely in a single header (`BitmapNesterC2.hpp`, ~1250 lines). The public API is one static method -- `arrange()` -- matching `BitmapNester`'s signature so callers can swap via a settings flag.

The pipeline has six phases, executed sequentially in `arrange()`:

**Phase 0 -- `build_cache`:** Iterates all `ArrangePolygon` items and populates a `NesterC2ItemInfo` per item: silhouette area and perimeter, convex hull area and perimeter, height from the 3D mesh, hardness score (hull_area / silhouette_area), and a composite priority score weighted as `caller_priority * 1M + silhouette_area * 100 + hardness * 10 + height * 1`.

**Pre-filter:** Items whose silhouette cannot be inscribed in the bed at any allowed rotation are removed from the cache. They stay UNARRANGED and go straight to spillover.

**Phase 1 -- `estimate_min_plates`:** Sum of silhouette areas divided by bed area times 0.82 packing efficiency, rounded up.

**Phase 2 -- `partition_items`:** Greedy least-loaded-bucket: items sorted by priority descending, each dropped into the group with smallest total silhouette area.

**Phase 3 -- `pack_as_island`:** The core placement loop. Items sorted height-descending (`std::stable_sort` with priority tiebreak). For each item: pre-rasterize all allowed rotations, coarse grid scan with hull-perimeter scoring, pixel-level refine pass around the coarse winner. Collision is bitmap AND via `BitmapNester::collides` (friend-class access). Bitmap resolution derives from the spacing parameter: `res = clamp(spacing / 2, 0.1, 0.5)` mm/px. Virtual plate = 2x bed at this resolution.

**Phase 4 -- `locate_island_on_plate`:** Shifts the island so its hull centroid lands at bed center, clamped so the island stays inside the bed.

**Phase 4.5 -- `compact_on_plate`:** The trash compactor. Rebuilds a composite bitmap from per-item `CompactItem` bitmaps (moved from `RotCache` at commit time -- zero-copy). Up to 50 iterations: XOR-remove each part, push inward from nearest bed wall, probe tangent directions with void-depth sensing on collision. Quality gate: pixel overflow count. Grid-snap post-pass (2mm).

**Phase 5 -- `recover_spillover`:** Items still UNARRANGED each get their own plate, centered on the bed.

---

## Where We're At

The C2 nester has a working pipeline from cache-build through placement, multi-plate partition, spillover recovery, and compaction. If you loaded parts and clicked Arrange today, you would get a valid arrangement. Parts would not overlap. They would land on the bed. The tall ones would cluster together. It would not crash.

### Test Suite

22 cases covering the full milestone ladder: M0 smoke tests, M1 phase-function unit tests, M2.4 hull-scored placement, M2.5 bbox-vs-silhouette proof, M2.6 tall-parts ordering and stability, M3.1 spillover recovery, M4.5 trash compactor (xor_remove roundtrip, neighbor preservation), and corpus discriminators (Category 6.4 identical-hull litmus, Category 1.4 L-shape interlock). The Cat 6.4 test proves solid squares and notched squares with identical convex hulls produce measurably different packings. The Cat 1.4 test proves two L-shapes nest into a 200mm hull perimeter instead of the naive 240mm bbox stack.

### Visual Renders

Grayscale PNG renders with multiplicative-darken compositing at 4x Nyquist oversampling confirm: zero overlaps in greedy placement, zero overlaps after compaction. The L-shapes render shows two L-shapes close but not perfectly interlocked -- a known limitation of the coarse hull-perimeter scoring at current stride. The mixed-squares render shows visible tightening from the compactor.

### What Works

- Hull-perimeter-scored greedy placement with bitmap AND collision
- Height-descending tall-parts grouping (Council-arbitrated)
- Trash compactor with wall pressure + void-pull tangent sliding
- Adaptive bitmap resolution (0.1-0.5mm/px from spacing setting)
- Clean overflow for oversized items
- Pre-filter for unfittable parts
- Spillover recovery

### What's NOT Wired Up

- **Exclude zones** -- no test covers placement around bed obstacles (purge tower, wipe station)
- **Multi-plate visual testing** -- renders only show single-plate results
- **Real STL fixtures** -- every test uses synthetic geometry. The tetris125 + 80 L-bracket benchmark that proved C1's value has not been run through C2
- **GUI integration** -- C2 is not selectable in the OrcaSlicer UI

### DLL Snapshots

Two exist in `builds/c2/`:
1. `20260412-092535-607a04c633-dirty-m4.5-trash-compactor` (from M4.5 commit)
2. `full-pipeline-adaptive-res` (building from `d986d65bec` — current HEAD)

---

## Performance Profile

**Phase 3 (placement) dominates runtime.** O(N x R x positions x hull_call) where positions ≈ 300 per item and hull calls touch 400-600 vertices at N=20 committed items.

| Spacing | res (mm/px) | Plate (256mm bed) | Memory | Est. time (40 parts) |
|---------|-------------|-------------------|--------|---------------------|
| 1.0 mm  | 0.5         | 1024x1024         | 128 KB | ~4 ms placement, ~9 ms total |
| 0.5 mm  | 0.25        | 2048x2048         | 512 KB | ~61 ms total |
| 0.2 mm  | 0.1         | 5120x5120         | 3.2 MB | ~1.6 s total |

**Compaction (Phase 4.5):** Converges in 20-40 iterations for N=50, ~22ms total. Per-part bitmap storage: ~80KB for 50 items. Zero allocations inside the loop.

**The biggest performance lever is scan region replacement** (replacing the rectangular search window with a dilated bitmap scan). The rectangular window wastes most candidate evaluations on interior positions that will collide. Dilation visits only perimeter-adjacent positions — a 3-5x reduction in position count.

---

## Where We're Going

### Path to PR (in order)

1. **M5 — A/B benchmark on real fixtures.** C2 vs C1 vs upstream libnest2d on tetris125 + 80 L-brackets and at least one realistic mixed-part set. The headline from C1 was "2 vs 3 plates." C2 needs to match or beat that.

2. **Scan region dilation.** Replace the rectangular search window with a dilated occupied bitmap. Eliminates the last bbox-derived decision in the placement loop. Also the biggest perf optimization (3-5x fewer candidate evaluations).

3. **SUM-bitmap integrity check.** After pack_as_island, assert no pixel in the composite exceeds 1. Catches accumulated-overlap bugs that no pairwise test finds.

4. **Compaction test battery.** Seven tests from the design memo: zero-overlap invariant, monotonic convergence, spacing preserved, determinism, iteration cap, quality gate, void-pull effect.

5. **GUI wiring.** C2 behind the `use_concave_shapes` flag so the arrange button invokes it. Half-day task. Produces the screenshot that closes the "show me it works" gap.

6. **PR draft.** With benchmark numbers, visual proof, and LOTC review report.

### Known Deferred Work

- FFT correlation for complete feasibility map (replaces grid scan entirely)
- Exclude zone collision (wipe tower, calibration zones)
- Physics compaction with contact sliding (the "full trash compactor" with tangent probing is landed; the advanced version with multi-body simultaneous resolution is post-PR)
- Split BitmapNester.hpp into .hpp/.cpp (Gandalf-crit from Sprint 2)
- Squash commit history into 6-8 logical commits for upstream review

---

*Report generated by LOTC formation: Gandalf (architecture) + Frodo (UX) + Legolas (performance). Collated 2026-04-12.*
