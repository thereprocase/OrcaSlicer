# PR C2 — Test Harness Validation Strategy

**Status:** Canonical harness validation strategy for C2 development.
Written from user direction 2026-04-11 (second session). Complements
`PR_C2_TEST_CORPUS.md`: the corpus is WHAT to test (7 categories of
packing quality); this document is HOW TO VALIDATE THE TESTER ITSELF
before trusting any quality measurement.

---

## The premise

**The test harness is more dangerous than the solver.**

If the overlap detector has a false positive, it rejects valid
placements and you get worse packing but never know why. If it has a
false negative, parts fuse on the print bed. Both failure modes are
silent. The solver looks like it's working — it's just working against
broken ground truth.

C1's test suite doesn't validate its own primitives. It asserts "no
overlap" by running the same `no_overlap` helper that's an approximation
with known tolerances, then asserts on a solver output that was
produced by the same rasterization code path the helper is checking.
Circular validation. This document is the plan to break the circularity
in C2 before trusting C2's own quality measurements.

---

## Phase 1 — Overlap detector in isolation

Before testing any packing, establish that **"these two shapes at these
positions do/don't overlap"** is a correct function. Build an exhaustive
test suite for collision detection separate from all packing tests.

### Analytically solvable cases

Test the overlap detector on shapes whose answer is computable by hand:

- **Two unit squares, one at origin, other at (0.5, 0):** overlap
- **One at origin, other at (1.0, 0):** touching — **decide a policy
  (overlap or not) and assert that policy consistently**
- **One at origin, other at (1.001, 0):** no overlap. At 0.5 mm
  resolution this is the bite point: 1.001 and 1.0 rasterize to the
  same pixel. Test at the resolution boundary deliberately.

The decision on "touching counts as overlap" is load-bearing —
document it once and enforce it in every test.

### Winding-number point-in-polygon as the oracle

The canonical validation primitive is **winding-number point-in-polygon**.
It's slow but provably correct for any simple polygon. Implement it
once, never optimize, never use in production. It's the oracle.

For any candidate placement of two shapes, sample a grid of points from
shape A and test each against shape B using winding number. If any
point is interior, overlap is confirmed. Every other overlap method —
bitmap AND, Separating Axis Theorem, FFT feasibility — gets validated
against this oracle on every test case in the suite.

This is the ground-truth generator. Everything downstream trusts it.

---

## Phase 2 — Rasterizer in isolation

The bitmap rasterizer needs its own tests independent of the overlap
detector, because every downstream calculation (area accounting, FFT
correlations, cluster perimeter) inherits whatever error the rasterizer
introduces.

### Known-area tests

- **10×10 mm square at 0.5 mm resolution** → count set bits. Must be
  exactly 400. Off-by-one means the scanline fill has a boundary bug.
- **Circle of known radius R** → count set bits, compare to πR². The
  error should be proportional to circumference divided by resolution,
  O(r/res). If it's worse than that, the rasterizer is broken. If it's
  exactly that, you have a quantified error bound you can propagate
  through every downstream calculation.
- **Shape + one-pixel-offset copy → AND** → overlap area should equal
  original area minus one row or column of pixels. Validates that the
  bitmap coordinate system is consistent between rasterization and
  collision — a surprisingly common bug is rasterizing in one
  coordinate convention and testing overlap in another.

### False-positive boundary

Generate two shapes with known analytical clearance — two circles of
radius R with centers `2R + gap` apart. Sweep the gap from 10 mm down
through the resolution limit. At each step the oracle says "no
overlap." The bitmap detector should agree down to some gap threshold
near the resolution.

**Record where it starts disagreeing.** That's the false-positive
boundary. It should be **exactly 1 pixel width** (the resolution). If
it's 2 pixels, there's a dilation bug — probably inflating the bitmap
by one pixel during rasterization, or using `≥` instead of `>`
somewhere.

### False-negative boundary

Two shapes that overlap by exactly 1 pixel. The oracle says overlap.
Does the bitmap agree?

- **0.1 mm overlap at 0.5 mm resolution:** the bitmap might miss it.
  Acceptable IF DOCUMENTED — but you need to know. The allowed
  false-negative envelope must be part of the harness spec.
- **0.6 mm overlap (just above resolution):** this MUST be caught. If
  it isn't, the rasterizer is wrong, not just imprecise.

---

## Phase 3 — Roundtrip consistency

These catch bugs that no individual unit test finds, because they check
the WHOLE pipeline (rasterize → place → score → commit → measure)
against a closed invariant.

### The bitmap SUM test

Take a known perfect tiling (the tetromino plate, where coverage is
provably 100%). Rasterize every piece. Instead of OR'ing the bitmaps,
**SUM the bitmaps as integers**. Every pixel must be exactly 1:

- Any pixel with value 0 is a **gap** (the tiling has a hole)
- Any pixel with value 2+ is an **overlap** (two parts share a pixel)

OR-based composition makes overlap invisible — you get the correct
union but lose the information that pieces shared pixels. SUM composition
preserves the multiplicity, which is the only way to catch "solver
placed two parts on top of each other and the bitmap AND still said
collision-free because of a one-pixel rounding."

**This is the single most important validation test for the whole
system.** Run it on every packing result during development — it's
cheap and it catches silent integration bugs.

### The complement test

Take the composite bitmap of a packing, invert it, count the empty
pixels, multiply by `res²` to get empty area. Compare to
`plate_area - sum_of_part_areas`. The numbers should agree to within
the rasterization error bound established in Phase 2.

If they disagree by more, either:
- A part was rasterized at the wrong size
- A part was placed at the wrong position
- The area accounting is wrong

Catches coordinate-system bugs that the overlap detector never sees
because the overlap detector is operating in relative coordinates.

### Topological containment test

Every placed part's bitmap must be **fully contained** within the
plate boundary bitmap. Computed as: `part_bitmap AND NOT plate_bitmap`
must be all-zero.

Every placed part's bitmap must have zero set bits in the exclusion
zones (wipe tower, bed clips, bed thermistor, etc.). Computed as:
`part_bitmap AND exclusion_bitmap` must be all-zero.

Simple AND/comparison checks. Run on every packing result.

---

## Phase 4 — FFT/NFP self-consistency gauntlet

Relevant once C2 introduces faster feasibility-map computation (FFT
correlation or NFP sweeps).

**If the FFT correlation says position (x, y) is feasible
(correlation value = 0), then placing the part there and running the
direct bitmap AND overlap check must confirm no overlap.**

Test exhaustively: for a given part and plate state, compute the full
feasibility map via FFT, then spot-check 100 random "feasible" positions
and 100 random "infeasible" positions with the direct bitmap AND
method. Any disagreement means the FFT pipeline has a bug — usually:
- A sign convention error in the conjugation (complex conjugate of
  one of the transforms)
- An off-by-one in the coordinate mapping between FFT output indices
  and plate positions (FFT output uses wrapped coordinates; plate
  uses natural coordinates)

Run this test every time the FFT path is touched. Never trust the FFT
output without a cross-check.

---

## Phase 5 — Rotation fidelity

### Area invariance

Rotate a square by 45°. The resulting bitmap should be a diamond. The
area in pixels should be within the established rasterization error of
the analytical area (which equals the original square's area — rotation
preserves area).

If the rotated area differs by more than the error bound, the
rotation→rasterize pipeline is introducing area errors. Common cause:
transforming vertices as doubles, then rounding, accumulates drift.

### Rotation roundtrip

Rotate a shape by θ, then rotate by -θ. XOR the result with the
original. Any set bits in the XOR are rotation-roundtrip errors.
These tend to accumulate if you're doing multiple rotations, so test
at the exact angles the solver will actually use (0, 45, 90, 135,
180, 225, 270, 315 if the solver supports eighth-turn rotations).

---

## Phase 6 — Adversarial geometry for the rasterizer

Shapes that break naive scanline fill:

- **Bowtie (self-intersecting polygon).** The rasterizer should either
  reject it or handle it consistently. If it silently produces garbage,
  everything downstream is suspect.
- **Polygon with a horizontal edge exactly on a scanline.** Even-odd
  rule ambiguity — the classic bug.
- **Polygon with a vertex exactly on a pixel center.** Is the boundary
  inclusive or exclusive? Document the decision and assert it.
- **A very thin sliver — one pixel wide at resolution but analytically
  nonzero area.** Does the rasterizer produce at least one set pixel,
  or does the shape vanish entirely? **Vanishing parts is a
  catastrophic false negative** — the solver thinks the plate is empty
  where a part actually exists.
- **Coincident vertices (degenerate triangle).** Must not crash.
- **CW vs CCW winding.** Must produce the same bitmap. If the
  rasterizer is winding-sensitive and the mesh projector doesn't
  guarantee winding, you get random failures.

---

## Phase 7 — Coordinate space mismatches

The pipeline has at least three coordinate spaces:

1. **Mesh coordinates** (the STL / 3MF model)
2. **Plate coordinates** (the ArrangePolygon world, in scaled coord_t)
3. **Bitmap coordinates** (pixel indices)

Every transition between spaces is a potential sign flip, offset
error, or scale error. The definitive test:

**Place a single L-shaped part at a known plate position. Read back
the bitmap. Manually verify that the set pixels correspond to the
right physical location. Do this for ALL FOUR CORNERS OF THE PLATE.**

Corner cases — literally — catch axis transpositions that
center-of-plate tests miss. (E.g. an off-by-one in the x coordinate
produces a 1-pixel shift at every position; placing a part at the
top-right corner makes the shift visible because it now falls outside
the plate boundary.)

---

## Continuous integration strategy

The test harness runs at different frequencies depending on its cost
and what class of bug it catches:

### Every commit
**Analytically solvable tests** — fast, catch regressions immediately.
- Overlap detector on axis-aligned shapes with known answers
- Rasterizer area tests on squares and circles
- Rotation roundtrip on the solver's exact rotation angles

### Every packing result (including production)
**Roundtrip consistency tests** — cheap, catch integration bugs that
unit tests miss.
- Bitmap SUM test (no gaps, no overlaps — every pixel exactly 1)
- Complement area test (empty area accounts correctly)
- Topological containment test (every part inside plate, outside
  excludes)

These run inside `arrange()` itself, gated by a debug flag, so every
real arrange call in dev builds validates its own output. Production
builds strip them out.

### Nightly
**Adversarial geometry suite + random polygon meta-test** — slower,
catches edge cases that only matter when real-world meshes produce
degenerate projections (which they eventually will, because users
upload arbitrary STLs).

---

## Relationship to the C2 test corpus

`PR_C2_TEST_CORPUS.md` defines the WHAT: 7 categories of packing
quality benchmarks plus the random-polygon meta-test.

`PR_C2_TEST_HARNESS_VALIDATION.md` (this file) defines the PRE-WHAT:
the harness primitives that every corpus test implicitly trusts. A
category-6.4 litmus test that says "the concave shape packs
differently than its convex hull" is meaningless if the overlap
detector says two touching squares overlap. You'd just be measuring
harness noise.

**Implementation order:**
1. Phase 1 (overlap oracle) before any category-1 test
2. Phase 2 (rasterizer area) before any efficiency assertion in any
   category
3. Phase 3 (SUM + complement) before trusting any multi-item result
4. Phase 5 (rotation) before any category that uses rotations
5. Phase 4 (FFT self-consistency) only when FFT lands (post-M5)
6. Phase 6 (adversarial) as a nightly suite, can lag
7. Phase 7 (coordinate spaces) as a one-time sweep when the pipeline
   structure stabilizes

---

## Current state vs. this plan

C1 has essentially none of Phase 1-4. Its `no_overlap` helper is a
polygon-intersection check with a 0.0001 mm² tolerance, which is its
own approximation layer. Its rasterizer tests check individual
operations but not the SUM invariant or the complement invariant. Its
coordinate-space sanity is implicit rather than explicit.

C2 inherits C1's rasterizer and primitives (by design — don't rewrite
what works), but C2 must ADD the harness validation layer before its
own quality assertions can be trusted. The implementation order above
is the M1/M2 prerequisite list — every C2 milestone has harness
validation as a gate.

---

## Memory marker

**Test harness bugs are silent and downstream. Always validate the
harness before validating the solver.** See
`memory/feedback_harness_over_solver.md`.
