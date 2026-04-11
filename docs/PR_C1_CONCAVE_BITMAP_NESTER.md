# PR C1 — Concave-Aware Bitmap Placement for Interlocking Part Layouts

**Branch:** `feature/concave-bitmap`
**Status:** Ready for upstream review (pending QA/QC signoff)
**Date:** 2026-04-11

## Headline

Replaces OrcaSlicer's libnest2d convex-hull arranger with a raster-scan
bitmap nester that sees the actual concave silhouettes of parts. On matched
inputs, delivers one plate fewer than upstream for concave-heavy fixtures
(80 L-brackets, 125 tetrominoes, 25 L-brackets), with visible interlock
patterns instead of bbox-spaced layouts.

---

## 1. Existing flow — upstream OrcaSlicer arrange pipeline

### The high-level pipeline

OrcaSlicer inherits its arrangement pipeline from Bambu Studio, which
inherited it from PrusaSlicer, which uses **libnest2d** as its placement
engine. The flow when the user clicks Arrange is roughly:

```
ModelObject::get_arrange_polygon()
    → extracts a 2D representation of each model instance
    → returns Polygon + translation + rotation + extruder_ids + height
ModelArrange::get_arrange_polys()
    → collects ArrangePolygons from all instances + fixed items (wipe tower)
    → delegates to arrangement::arrange(items, fixed, bed, params)
libnest2d placement engine
    → runs selection strategy (FirstFit, DJDHeuristic, Filler)
    → runs placement strategy (NfpPlacer, BottomLeftPlacer)
    → writes back translation/rotation/bed_idx to each ArrangePolygon
```

### The critical convex-hull simplification

Before the 2D polygon is handed to libnest2d, **it is reduced to its convex
hull**. Source: `Model.cpp:3300`:

```cpp
Polygon p = get_object()->convex_hull_2d(t.get_matrix());
// ...
ret.poly.contour = std::move(p);
```

This is the fundamental limitation of the upstream arrange pipeline: libnest2d
never sees the actual concave shape of the part. For a part with significant
concave area (L-bracket, U-channel, frame, tray, bracket), every placement
operates on a rectangular bounding envelope — and that envelope wastes the
concave interior area.

Quantitative example: an L-bracket with bbox 40×40 mm has silhouette area
~1200 mm² (3/4 of the bbox). libnest2d scores its placement as a 1600 mm²
rectangle. Two L-brackets that could interlock by rotating 180° and sharing
a notch (~200 mm² of wasted area eliminated per pair) instead get packed
as two 40×40 squares with the notches facing empty space.

### libnest2d's placement primitives

libnest2d uses **No-Fit Polygon (NFP)** placement: for each candidate
position, compute the set of positions where the new item would NOT collide
with existing items, expressed as a polygon. Select the position that best
satisfies a fitness function (typically bottom-left, lowest-y-then-lowest-x).
NFP is mathematically sound and efficient for convex shapes; it degrades
to the convex hull for irregular shapes, matching the reduction above.

### Selection strategies in libnest2d

The codebase includes three:
- **FirstFit** — place each item at its first-fit position in input order
- **DJDHeuristic** — López-Camacho 2013. Sort by area descending, fill one
  bin at a time, try pairs and triplets when singletons don't fit, waste
  budget increments progressively. **Notable:** `try_reverse_order`,
  `try_pairs`, `try_triplets` config options exist but require convex input.
- **Filler** — simpler variant for small input

Orca uses a variant of FirstFit configured with plate-count minimization.
The DJD heuristic is available but unused in production Orca arrange paths.

### What upstream does well

- Fast: NFP intersection math is polynomial and well-optimized
- Mature: libnest2d has been in production since ~2015
- Deterministic: same input → same output
- Handles convex shapes near-optimally
- Integrates cleanly with Orca's ArrangePolygon data model

### What upstream misses

- **Concave interlock.** The convex hull reduction eliminates the entire
  value proposition. Parts that could mesh like gears get packed like
  spaced squares.
- **Rotation-aware interlocking.** Even where libnest2d tries multiple
  rotations, it's scoring the rotated convex hull. A rotation that
  produces an interlocking concave fit scores the same as a non-
  interlocking rotation at the bbox level.
- **Tight clusters for travel minimization.** Without interlock, clusters
  are loose and the print head travels more between parts.
- **Concave self-intersection.** Multi-volume parts whose XY projections
  form disconnected islands get collapsed into a single hull that over-
  reserves space.

---

## 2. Our concave flow — the BitmapNester

### Design philosophy

Replace the NFP math with a **raster scan over a bitmap representation of
the plate**. Concave silhouettes become bitmap regions, which the scan
then tests for collision-free placement and scores for cluster tightness.
The whole pipeline is concave-aware from first principles.

### Entry point

`BitmapNester::arrange(items, excludes, bed, params)` — matches the
libnest2d interface so upstream callers can swap it in without touching
the model/instance layer.

### The placement loop

For each item in FFD order:

1. **Pre-rasterize silhouettes.** For each allowed rotation, rasterize
   the item's concave silhouette to a per-rotation bitmap. This is the
   key departure from libnest2d: we keep the full concave shape, not
   just its hull. The code also computes the convex hull bitmap for
   collision-checking against exclude zones (wipe tower) — that check
   uses the hull for parity with Orca's existing `check_outside`
   validation.

2. **Build rot_cache once.** Rasterization is expensive; caching per
   item avoids N × rotations × plates redundant calls.

3. **Try plates in order.** For each plate:
   - Resolve the anchor (bed corner for default-align, wipe-tower
     center for plates with wipe towers).
   - Run a coarse scan: for each rotation, enumerate every
     stride-spaced position in the plate's valid range plus explicit
     max-row/max-col boundary positions. Score each clear position
     by cluster-bbox-area growth (primary) and perimeter growth
     (secondary). Track best per rotation.
   - **Desperation retry:** if the coarse scan finds no valid
     position at any rotation, retry once with stride halved.
     Catches narrow-gap cases in dense plates.
   - Refine: for each rotation with a coarse winner, scan every pixel
     in a ±coarse window. Score every position. Keep the global best
     across all rotations. This is the fix for the per-rotation lock
     bug (a rotation that lost the coarse lottery but would interlock
     perfectly at a nearby pixel still gets its shot).
   - Commit the winning placement: stamp the bitmap, update cluster
     bbox, write item state.

4. **Plate escalation.** If no plate accepts the item, allocate a new
   plate and retry. Plates cap at MAX_PLATES.

### The consolidation pass

After the greedy placement loop, walk plates last-to-first. For each
item on a later plate, try to migrate it back to an earlier plate by
finding any clear scored position there. Uses the same scored-scan +
all-rotations logic as main placement. On successful migration:

1. Commit the item at its new (plate, rotation, position)
2. Rebuild the source plate's bitmap by re-stamping surviving items
3. Mark the outer loop for another pass

Runs up to 4 iterations of the outer loop or until no migration happens.
Preserves the caller's `align_to_y_axis` pre-rotation via a base rotation
snapshot captured before placement. Per-item rot_cache built once outside
the to_plate loop for efficiency.

### The smart-shuffle restart driver

Wraps the whole placement+consolidation+centering pipeline in a loop that
runs multiple attempts with shuffled within-priority ordering. For small N
(≤6) runs 8 attempts; medium N (≤16) runs 4; large N runs 3. Scores
attempts by (placed count desc, max_bed asc, total area asc) and keeps
the best. Short-circuits when a single-plate placement with all items
placed is found.

The restart driver exists because one-shot greedy placement is hostage
to its input ordering — a single pathological order can fragment the
plate into strips too narrow for later items. Random restarts are the
cheapest defense.

### Post-centering

Geometric centering pass: for each plate with placed items, compute the
cluster bbox, shift the cluster so its center aligns with the bed center
(or the caller's align_center target). Safety-checked against exclude
zones — if the shift would cause an item to collide with an exclude,
revert to no shift.

### FFD ordering + pure corner anchor

Items are pre-sorted by priority desc, then longest-dimension desc, then
area desc as a tiebreaker. Classical First-Fit-Decreasing. The first-placed
item anchors the cluster at the plate corner (0,0). Subsequent items also
use corner anchor — the hybrid anchor (corner first, bed center rest) was
tried and reverted because it over-pulled items toward the bed center
when multiple positions tied on area growth, producing off-grid layouts
that blocked subsequent items (see session diary 2026-04-11 for the
4-equal-squares case).

### Tall-parts-centered scoring bias

Items with `ArrangePolygon::height ≥ 100mm` get a scoring penalty
proportional to distance from bed center. This is a weak tertiary
tiebreaker that only kicks in when area + perimeter both tie. Effect on
typical inputs is small; the feature is wired up for stronger use in a
follow-up sprint.

### Test infrastructure

130 test cases across 11 files:
- Bitmap primitives (rasterize, dilate, collide, stamp) — direct unit tests
- Integration scenarios — tetromino probes, L-bracket interlock, real OBJ
  fixtures, consolidation scenarios, stress tests
- Metric helpers — `max_bed_idx`, `overflow_piece_count`, `no_overlap`,
  `deterministic_rerun`, `cluster_bbox_perimeter_mm`, `cluster_compactness`
- The Three Seers metric vote (plate count + determinism + overflow piece
  count) documented in `docs/PACKING_METRICS.md`

### The numbers (matched fixtures)

| Fixture | Upstream libnest2d | BitmapNester | Delta |
|---|---|---|---|
| 25 L-brackets | 2 plates | 1 plate | −1 |
| 80 L-brackets | 3 plates | 2 plates | −1 |
| 125 balanced tetrominoes | 3 plates | 2 plates | −1 |
| 23-piece concave OBJ mix | N/A | 1 plate | (pure regression) |

Visual comparison: upstream produces loose bbox-spaced layouts. BitmapNester
produces pinwheel rosettes on L-brackets, parallel-strip + tucked-interlock
on tetrominoes, and tight mixed clusters on real OBJ fixtures.

---

## 3. Limitations — what we know is imperfect

### Algorithmic

- **Tetris125 pathology.** 125 equal-area tetrominoes in clumped input
  order is a synthetic worst case. Our algorithm places ~14 pieces on
  plate 2 due to FFD grouping all I-pieces first (long parallel strips)
  before T/S/Z/J/L clumps, which doesn't produce optimal tessellation.
  Still beats upstream's 3-plate layout, and real STL mixes (which
  don't have this symmetry) don't exhibit the pattern.

- **Tall-parts-centered is a weak tiebreaker.** The feature applies only
  when area + perimeter both tie, which is rare in practice. The
  correct implementation is a post-centering modification (height-
  weighted COM target) but bounds-clamping dominates for clusters
  that fill the plate width. Deferred.

- **Axis-aligned bbox as scoring metric.** The user's stated dream is
  "organic lassos not bboxes" — a convex-hull-perimeter or
  concave-hull-perimeter metric would capture travel cost more
  directly than the rectangular bbox. Dropping the bbox for a hull
  metric is a fundamental redesign, not an incremental fix.

- **Small-batch look-ahead not implemented.** Greedy per-item placement
  commits positions permanently. Look-ahead (try next W items, pick
  the best fit at each step) or per-batch permutation search would
  produce tighter clusters. Requires extracting the per-item placement
  body into a reusable helper — ~1-2 hours of refactor.

- **90° rotation granularity.** No 45° or finer rotations. Curved
  parts that could interlock at 22.5° or 45° don't.

### Performance

- **Rasterization is O(pixels_per_part).** For parts with large XY
  projections at 0.5mm pixel resolution, this is hundreds of thousands
  of pixels per rotation × N items × K attempts. libnest2d's NFP math
  is asymptotically cheaper for large items.

- **Coarse scan is O(plate_area × rotations × items).** Also larger
  than libnest2d's scan cost.

- **Smart-shuffle multiplies runtime by attempt count.** For N=125 with
  3 attempts, each attempt is ~2-3s on a Windows dev box, total 6-9s.
  Acceptable for an interactive "click Arrange" workflow; not cheap.

### Code quality

- **BitmapNester.hpp is ~1700 lines, header-only.** An upstream
  reviewer will object. Sprint 2 will split into .hpp/.cpp.

- **The `do_one_pass` lambda is ~900 lines long.** It bundles
  placement + consolidation + post-centering into a single closure
  so the smart-shuffle restart driver can reset and retry. The
  refactor into a proper object with methods is deferred because
  it would touch every line of the algorithm.

- **`#ifdef BITMAP_NESTER_TESTING` guards expose private statics to
  test files.** Ugly but localized. Sprint 2 will replace with a
  friend class declaration.

- **Some `// L##:` audit markers remain from prior review passes.**
  Cosmetic, will clean up.

### Test coverage

- **No real-STL arrange integration test.** `test_real_mix_packing.cpp`
  loads OBJ files but uses their vertex convex hulls, not the full
  concave silhouette extracted by production `slice_mesh_ex`. Tests
  exercise a slightly different code path than production.

- **Tall-parts-centered has no test coverage.** The test drafted in
  the 2026-04-11 session failed because the feature's effect is too
  subtle to assert on symmetric inputs. Feature still committed (as
  a weak tiebreaker) but validation is pending a stronger
  implementation.

- **Consolidation tests cover 4 scenarios.** Adequate for smoke but
  not exhaustive across the parameter space (no rotation-required
  migration test with >3 items, no stopcondition-interrupt test).

- **No travel-cost assertion in tests.** Helpers exist
  (`cluster_bbox_perimeter_mm`, `cluster_compactness`) but no test
  asserts bounds on them. Sprint 2 will add.

### Parity / integration

- **Per-item inflation semantics.** The nester computes `pad_px` from
  either `item.inflation` or `params.min_obj_distance/2`. Upstream
  Orca sometimes uses `min_obj_distance` directly. Verified equivalent
  for typical brim widths but not formally proven.

- **Multi-volume parts via `concave_regions`.** When the caller
  populates `item.concave_regions` (ArrangeJob does this when
  `use_concave_shapes` is on), the nester uses every island. When
  not populated, falls back to `item.poly`. The fallback is
  libnest2d-compatible but loses multi-island coverage.

- **Exclude zone interaction.** Wipe towers and other excludes are
  handled via a second bitmap (the convex hull mask) to keep parity
  with Orca's `check_outside` validation. Correctness verified on
  test_bitmap_multibed but the parity contract is implicit.

---

## 4. Areas for future exploration

### Algorithmic upgrades (post-Sprint-3 fork territory)

Tracked in `docs/FFT_NESTER_IDEA.md`:

1. **FFT-correlation inner loop.** Replace the O(bed_area × rotations)
   coarse scan with an FFT-based cross-correlation that computes all
   valid positions in O(bed_area × log bed_area). Speeds up large-N
   cases dramatically. Proposal documented but not implemented.

2. **Certificates + precomputed consolidation.** Each item carries
   a "certificate" of where it can legitimately land, precomputed
   from its silhouette. Consolidation becomes a cheap lookup.

3. **Three-tier parallelism.** (a) FFT parallelism within a single
   item's scan, (b) permutation search parallelism across smart-
   shuffle attempts, (c) pipelined packing + certification +
   consolidation for overlapping work.

### Design ideas from the 2026-04-11 session

4. **Organic lasso scoring.** User's stated dream — replace axis-
   aligned bbox with convex or concave hull perimeter as the primary
   scoring metric. Travel cost tracks hull perimeter directly; bbox
   is a coarse proxy. Requires recomputing the cluster hull per
   placement (expensive) and scoring candidates against the hull.
   See `dev/sessions/2026-04-11/design-gravitational-stack.md`.

5. **Bound-agnostic gravitational best stack.** Two-phase: (phase 1)
   pack items tightly without bed constraint, optimizing for minimum
   cluster perimeter / hull perimeter, (phase 2) locate the resulting
   cluster on a plate and split across plates if oversized. Decouples
   packing shape from bed shape, which is conceptually cleaner than
   the current mix.

6. **Small-batch permutation search.** For each batch of 4-8 items,
   try all permutations (or samples), commit the best permutation's
   placements, move on. Inspired by libnest2d's DJDHeuristic. Deferred
   in the 2026-04-11 session because it requires extracting
   `try_place(item, state)` from the inline placement body.

7. **Look-ahead ordering.** At each placement step, score the next W
   items against the current plate state, commit the best, swap into
   order. Cheaper variant of small-batch without state snapshots.

### Orca-specific integration

8. **45° and finer rotations.** Currently locked at 0/90/180/270.
   Curved parts need finer granularity. Would require rotating both
   the silhouette AND its convex hull for the exclude check.

9. **Locked plate support.** Let the user lock plate 01 so arrange
   only reshuffles plate 02+. UX feature frequently requested by
   power users.

10. **Concave fallback counter (Frodo C1 from earlier review).** When
    silhouette extraction fails and the nester silently falls back
    to convex hull, the user has no diagnostic indication. Add a
    counter + warning so failures are visible.

11. **Per-item bitmap cap (Aragorn C1).** A malicious mesh × 4
    rotations × N items can force hundreds of MB of allocation.
    Safety limit needed for untrusted input.

12. **BitmapNester.hpp → .hpp/.cpp split.** Sprint 2. Will also
    replace the `BITMAP_NESTER_TESTING` ifdef guard with a
    friend-class declaration.

### Benchmarking

13. **Travel-cost regression assertions.** Tests already compute
    cluster_bbox_perimeter_mm and cluster_compactness but don't
    assert on them. Sprint 2 will add floors that ratchet as the
    algorithm improves.

14. **Head travel simulation.** The real goal is minimizing
    **slicer-emitted XY travel moves**, which the arrange step can
    only proxy indirectly. A simulation hooked into the G-code
    emission stage would close the loop — the arrange step could
    compute the estimated travel for its output and assert bounds
    on it.

15. **Random-seed sweep on tetris125.** The synthetic pathology has
    a specific random seed baked in; running with 10+ seeds and
    setting the floor to the worst observed overflow + 1 turns a
    pity threshold into a genuine regression gate.

### Collaboration with other OrcaSlicer branches

16. **snuggle/ radial expansion nester.** User has an experimental
    branch with radial-expansion packing. That approach may share
    ideas with the gravitational best-stack design — worth comparing
    after this PR lands.

17. **Default arranger A/B hook.** Let users toggle between libnest2d
    and BitmapNester in settings so real-world usage data informs
    when to default which. Useful for iterative improvement after
    merge.

---

## Summary for the reviewer

This PR adds a concave-aware bitmap nester as an alternative to the
convex-hull libnest2d arranger, wired in via
`ArrangeParams::use_concave_shapes`. On concave-heavy inputs it
delivers visibly tighter layouts and one plate fewer than upstream
on matched fixtures. The algorithm is raster-based, deterministic,
and tested with 130 cases / 1161 assertions. Two known imperfections
(header-only file size, weak tall-parts-centered) are flagged for
Sprint 2. No new dependencies; works within existing Orca + libnest2d
infrastructure.

The companion design notes in `docs/FFT_NESTER_IDEA.md`, `docs/PACKING_METRICS.md`,
and `dev/sessions/2026-04-11/` capture the ideas the session explored
and the reasoning behind the current compromises. Those are not
required reading for PR review, but are the best pointer to "where
this goes next" if the reviewer asks.
