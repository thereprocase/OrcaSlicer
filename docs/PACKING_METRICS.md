# Packing Quality Metrics — Brainstorm

Captured 2026-04-11. The previous test harness (`test_tetromino_packing.cpp`)
asserted on `material_area / cluster_bbox_area`. That metric is wrong: it
penalizes clusters with rounder shapes (e.g. center-greedy growth) even
when the actual bed coverage and plate count are identical to a tighter
rectangular cluster. The user said "bbox is the least interesting possible
metric" and asked for alternatives. This document is the brainstorm.

The Three Seers (Sauron, Gandalf, Frodo) were asked to vote on which
combination to adopt for the regression suite. Their recommendation
governs the test rewrite.

## What we're trying to measure

The end-user is a 3D-printer operator. **Most arrange operations are
single-plate** — the user adds parts that fit on the bed and clicks
arrange. Plate count is binary in this case (1) and uninteresting as a
quality metric. The thing the user actually feels every time is **print
head travel** between parts during the print, and that is directly
proportional to **cluster compactness**: the tighter the parts pack
together, the less the head moves between them.

**Primary goal: minimize print head travel within a plate.** This is
proxied by geometric compactness measures — convex hull tightness,
adjacency, perimeter, convex deficiency. Plate count is a secondary
sanity-check for the overflow case (when parts genuinely don't fit).

Secondary goals:
- Visually balanced layout (parts centered, not corner-hugging)
- Robust under input perturbation (same parts in different order →
  same answer)
- Honest comparability against the default libnest2d nester
- Catch the overflow regression when it happens (plate count metric)

## Plate-level metrics (the goal-aligned ones)

### 1. Plate count
The integer number of plates used. The single number that actually
matters for a printer operator. Test: `arrange_count(items) == K_min`.

**Pros:** matches user goal exactly. Boolean per scenario.
**Cons:** insensitive to "almost made it" cases — a one-piece overflow
gets the same score as a half-bed overflow.

### 2. K_min ratio
`K_min / actual_plates_used`. 1.0 = optimal. 0.5 = doubled the plate
count.

**Pros:** unitless, no shape dependence, comparable across scenarios.
**Cons:** discrete at small N — for K_min=1, you get either 1.0 or 0.5,
no granularity.

### 3. Plate-0 utilization
`sum(part_area on plate 0) / bed_area`. The "did we fill the first
plate" metric.

**Pros:** continuous, sensitive to packing tightness, no shape penalty.
**Cons:** depends on bed size relative to part sum — needs tight bed
sizing per scenario.

### 4. Overflow piece count
`count(items where bed_idx > 0)`. The "regret" metric: pieces you wish
were on plate 0.

**Pros:** integer, intuitive, directly tied to user complaint
("the spillover").
**Cons:** treats a 1-overflow and a 100-overflow as different by 100,
not by importance.

### 5. Last-plate utilization
% full of the LAST plate (highest bed_idx). Diagnostic: if it's 5%,
you almost made it on N-1 plates and one piece spoiled it. If it's 95%,
the last plate was a real one anyway.

**Pros:** distinguishes "near miss" from "genuine overflow."
**Cons:** confusing for K_min=1 case where there is no last plate.

## Geometric metrics (cluster shape, not goal-aligned but informative)

### 6. Convex-hull compactness
`sum(part_area) / convex_hull(cluster).area`. The convex hull instead
of the bounding box; round and rectangular tight clusters score the
same.

**Pros:** rotation-invariant, shape-fair.
**Cons:** computing the convex hull is O(N log N) per assertion;
still penalizes legitimately scattered placements.

### 7. Inscribed-empty-rectangle vs largest unplaced part
Find the biggest rectangle that fits inside the empty regions of plate 0.
If it's larger than any unplaced part, **the consolidation pass failed**.
Direct test of consolidation effectiveness.

**Pros:** explicitly tests the consolidation property we care about.
**Cons:** computing the largest empty rectangle is non-trivial (O(N²)
naive); requires bitmap analysis.

### 8. Wasted-space largest blob
Area of the biggest contiguous empty region inside the cluster outline.
Diagnoses "I have a notch-shaped hole exactly the size of the missing
piece."

**Pros:** intuitive; visual.
**Cons:** same complexity issue as #7; requires connected-component
analysis on the bitmap.

### 9. Adjacency / edge-contact length
Total mm where one part's outline touches another's. The literal
"interlocking" measure. Bitmap-friendly: count border bits that touch
other parts' bits.

**Pros:** measures the actual property the concave nester optimizes
for. Cheap on bitmaps.
**Cons:** depends on part shape; tetrominoes naturally have lots of
contact regardless of nester quality.

### 10. Convex deficiency
`convex_hull_area - filled_area`. Measures wasted space *inside* the
cluster outline. Goes to zero for a perfect tessellation in any rotation.

**Pros:** rotation-invariant, sensitive to interior gaps.
**Cons:** doesn't distinguish "one big interior hole" from "many small
ones."

## Comparative metrics (baseline vs us)

### 11. Beat-libnest2d ratio
Run libnest2d on identical input, assert `our_plate_count <= libnest2d_plate_count`.

**Pros:** honest comparison, no thresholds to calibrate, survives
algorithm changes. The number that goes in the PR description.
**Cons:** requires libnest2d baseline plumbing (currently broken in
test_tetromino_packing.cpp — `baseline_placed == 0` failure mode);
makes test runtime ~2× longer.

### 12. Per-fixture plate-count regression
For each curated test fixture (e.g. `tetris_plate_240x240_balanced.stl`),
record the achieved plate count once and assert `<=` going forward.

**Pros:** zero baseline plumbing required, catches all regressions
including subtle ones.
**Cons:** baseline number is locked in once, never improves automatically.

## Robustness metrics (orthogonal to single-run quality)

### 13. Determinism
Same input → byte-identical output. Run twice, assert positions match.

**Pros:** crucial for reproducibility, debugging, CI.
**Cons:** doesn't measure quality, only consistency.

### 14. Order-invariance
Same items in randomized input order → same plate count (not
necessarily same positions). Distinguishes "we got lucky on input
order" from "we found a real packing."

**Pros:** detects fragility.
**Cons:** the algorithm is single-pass greedy and IS order-sensitive;
this metric would currently fail.

### 15. Rotation-invariance
Rotate the entire input set by 30°, re-arrange. Plate count should be
the same.

**Pros:** detects directional biases (e.g. "cluster always grows
rightward").
**Cons:** requires rotated input fixtures; doesn't apply to allowed-
rotations sets that include the chosen rotations.

## Diagnostic (post-fork, FFT branch)

### 16. Gap-fillability count
For each overflow piece, run the bbox-cofit certificate against every
earlier plate's free region. Count how many "should have fit" (cofit
says yes) but didn't. A non-zero number means the consolidation pass
missed something.

**Pros:** quantifies consolidation effectiveness directly.
**Cons:** requires the bbox-cofit primitive and free-region calculation;
not yet implemented.

### 17. Wall-clock time per arrange
Bounded budget. Useful for the FFT fork's parallelism work.

**Pros:** measures cost.
**Cons:** machine-dependent; needs CI normalization.

## The Seers' vote (2026-04-11)

| Metric | Sauron | Gandalf | Frodo | Tally |
|---|:---:|:---:|:---:|:---:|
| #12 Per-fixture plate-count regression | ✓ | ✓ (primary) | ✓ | **3** |
| #13 Determinism | ✓ | ✓ | ✓ | **3** |
| #6 Convex-hull compactness | ✓ | rejected | rejected ("the new bbox") | 1 |
| #10 Convex deficiency | ✓ | rejected ("jargon") | rejected | 1 |
| #11 Beat-libnest2d | rejected ("baseline broken") | ✓ (headline) | rejected ("doubles runtime") | 1 |
| #4 Overflow piece count | (subsumed by #12) | (subsumed by #12) | ✓ (user-pain) | 1 |
| #9 Adjacency / edge-contact | rejected ("rasterizer not nester") | rejected ("FFT-fragile") | ✓ (loose-grid catch) | 1 |

**Unanimous winners:** #12 and #13. Every other pick was rejected by at
least one seer. The geometric compactness metrics (#6, #10) that I
initially favored were dismissed by Frodo as "the new bbox" — same
failure mode as the metric we're replacing. The libnest2d baseline
(#11) is gated on plumbing that's currently broken. Adjacency (#9) is
considered rasterizer-coupled by Sauron and FFT-fragile by Gandalf.

## Decision

Adopt three metrics for the regression suite:

1. **#12 Per-fixture plate-count regression** (primary, unanimous)
   - Lock the achieved plate count per curated fixture; assert `<=`
     going forward. One integer per fixture in a baseline table.
   - Catches overflow regressions (the original loose-grid bug) and
     any future regression that increases plate usage.
   - Zero new infrastructure. Survives algorithm rewrites trivially.
   - Headline fixture: `tests/data/tetris_plate_240x240_balanced.stl`
     locked at the value the current build achieves on first run.

2. **#13 Determinism** (invariant, unanimous)
   - Run arrange twice on identical input. Assert byte-identical
     `(bed_idx, rotation, translation)` tuples for every item.
   - Catches state-machine bugs in consolidation, RNG leaks,
     unordered_map iteration, thread races. Cost: one extra arrange
     call per test fixture, one tuple comparison.
   - Required for any future parallelization work — without
     determinism the FFT-fork's parallel restarts can't be tested.

3. **#4 Overflow piece count** (user-aligned, Frodo's championed pick)
   - `count(items where bed_idx > 0)`. Integer. Zero or not zero.
   - More granular than #12: catches "we used the same plate count
     but more pieces overflowed" (only possible if a fixture also
     overflows on the baseline). For tetromino fixtures sized to
     fit on plate 0, asserts `== 0`.
   - Frodo's argument: this is the literal "spillover" complaint
     a user files a bug about. Both other seers consider it
     subsumed by #12 but neither rejects it.

## What this means for the tetromino tests

Drop the `bbox_density >= X` assertions entirely. Replace with:
- `REQUIRE(arrange_max_bed_idx(items) == 0)` (#12 + #4 combined)
- `REQUIRE(deterministic_rerun(items, params, bed))` (#13)

For each scenario, choose a bed size that fits the parts. The bed
size becomes the implicit threshold — if the algorithm regresses,
it overflows.

## What this means for the tetris STL fixture

`tests/data/tetris_plate_240x240_balanced.stl` becomes a `[BitmapMesh]`
test. Load mesh → arrange → assert `current_plate == 0` (per #12).
Lock that as the baseline.

## What we're explicitly NOT doing

- No bbox density of any kind, ever again
- No convex hull or compactness math (Sauron lost this vote 1-2)
- No libnest2d head-to-head until the baseline plumbing is fixed
  (Sprint 2 follow-up)
- No adjacency / edge-contact (rasterizer-coupled, FFT-fragile)
- No order-invariance or rotation-invariance tests yet (the algorithm
  is admittedly order-sensitive; these would currently fail)
- No wall-clock time assertions (machine-dependent flake)

