# PR C2 Test Corpus — the driver

**Status:** Canonical test strategy for PR C2 (organic lasso nester).
Written from user direction 2026-04-11: "The old test suite is peanuts.
This is the driver."

The existing 130-test suite validated correctness (no crashes, no
overlaps, deterministic reruns) but it didn't measure **how good** the
packer is. That was fine for PR C1, which aimed to beat libnest2d by
one plate on matched inputs. C2 aims higher — approach known-optimal
packings on shapes where the optimum is provable, and measure closeness
to optimal on shapes where it isn't.

Every test case in this corpus has a **known-optimal or provably-bounded
answer**, so efficiency is measurable as a single number per fixture.
That turns the corpus into a leaderboard and turns development into a
numerical optimization problem instead of a vibes check.

---

## Scoring methodology

For every test case, the harness records three numbers:

1. **Packing efficiency** = `filled_area / bounding_rectangle_of_arrangement`.
   The fraction of the arrangement's bbox that's actually covered by
   parts. Higher is tighter. This is the PRIMARY quality metric.

2. **Time to solve**. Wall-clock ms from `arrange()` entry to exit.
   Secondary metric — quality wins, but a 10-minute solver is a
   non-starter for interactive UX.

3. **Ratio of achieved efficiency to known optimum**. For categories
   with provable optima (circles, crosses, L-pairs), this is the
   ground truth. For others, it's left as NaN and the leaderboard
   tracks absolute efficiency across runs.

**Weighting is category-dependent.** If you're tuning rotation search,
weight categories 1 and 6 heavily. If you're tuning the placement
scan or FFT inner loop, category 4 (the island shapes) stresses
bitmap resolution and cavity-filling opportunities. Every category
tells the solver something different about where it's weak.

---

## Category 1: Perfect-fit pairs

**Purpose:** Binary pass/fail cases. The solver either finds the perfect
fit or it doesn't — no ambiguity about optimality, no partial credit.

**Progression matters:** if the solver can't find a sine-wave jigsaw
fit, nothing else in the corpus matters yet. If it nails jigsaw pairs
but fails yin-yang, the rotation search is weak. Run this category
first in every benchmark — it gates the rest.

### 1.1 Jigsaw pairs

Take a rectangle, cut it with an arbitrary curve, you get two shapes
that reconstitute the rectangle. Progression:

- **1.1.a Simple sine cut.** A 100×60 rectangle split by a single
  sine period of amplitude 10 mm. The two halves fit back together
  with zero waste if the solver rotates one 180° and aligns.
- **1.1.b Multi-period sine cut.** Same rectangle, 3 periods,
  amplitude 15. More constraint on the rotation.
- **1.1.c Square-wave cut.** Same rectangle, square wave profile
  with 5 teeth. Interlocks at one specific position only.
- **1.1.d Fractal cut.** Rectangle split by a fractal boundary
  (Koch snowflake edge, 3 iterations). The hard version.

Pass criterion: packing efficiency ≥ 0.99 (zero waste within
rasterization noise).

### 1.2 Yin-yang

Two identical S-curves that nest into a circle. Tests rotation search
— only one relative angle works.

- **1.2.a Classical yin-yang.** Two comma shapes that form a disk.
  Only a 180° rotation produces the fit.
- **1.2.b Offset yin-yang.** Same shape but with the comma dots
  removed. Simpler geometry, same rotation requirement.

Pass criterion: packing efficiency ≥ 0.95.

### 1.3 Dovetail

A trapezoid with a triangular notch and its exact complement. Classic
joinery shape — the concavity is the whole point. Convex-hull solvers
are blind to the mating fit.

- **1.3.a Symmetric dovetail.** Trapezoid base 60, top 40, height 30,
  notch 20×15. Complement is the inverse.
- **1.3.b Asymmetric dovetail.** Off-center notch, non-square
  complement.

Pass criterion: packing efficiency ≥ 0.98.

### 1.4 Tetromino complements

An L and its complement within a 3×2 rectangle. These already exist
in `test_tetromino_packing.cpp` as an interlock test but are extended
here with their complement forms and pass/fail efficiency
assertions.

Pass criterion: packing efficiency = 1.00 exactly.

---

## Category 2: Known-optimal tilings

**Purpose:** Shapes with mathematically proven optimal packing
densities. The solver's ratio-to-optimum is a hard number.

### 2.1 Circles

Optimal hexagonal packing density = π/(2√3) ≈ 0.9069.

- **2.1.a N=16 identical circles** on a bed sized to their hexagonal
  optimum. Grid packing gets ~0.785; hexagonal gets ~0.907. The
  solver must do better than grid to pass.
- **2.1.b N=50 identical circles.** Larger N, more rotation-
  insensitive. Verifies scalability.
- **2.1.c Mixed radii.** 10 large + 20 small. No closed-form
  optimum but the solver should clearly beat grid packing.

Pass criterion: efficiency ≥ 0.85 (beats grid packing, approaches
hexagonal).

### 2.2 Regular pentagons

Pentagons don't tile the plane. Optimal packing density ≈ 0.921,
involves two orientations. Tests whether the solver discovers the
two-orientation trick.

Pass criterion: efficiency ≥ 0.85. Stretch: ≥ 0.90.

### 2.3 Crosses (plus signs)

Plus-sign polyominoes tile at 100% in a brick-lay pattern. Any gap
is solver failure.

- **2.3.a Standard cross.** Five-square plus sign, 10 copies.
- **2.3.b Fat cross.** 3-square arms, 10 copies.

Pass criterion: efficiency ≥ 0.98 (essentially 1.00 with boundary
tolerance).

### 2.4 Polyhedron nets

Any polyhedron net has a complex outline with multiple concavities
and known tilings with its mirror image.

- **2.4.a Cube net** (6 squares in a cross pattern). Tiles with
  its complement.
- **2.4.b Herschel enneahedron net.** The most concave of the
  classical polyhedra nets.
- **2.4.c Dodecahedron net.** High concavity, curved outline
  approximated as polygon.

Pass criterion: efficiency ≥ 0.80 (well above any convex-hull
baseline).

---

## Category 3: Concavity depth spectrum

**Purpose:** Same base shape with progressively deeper concavities.
Tells the solver where its concavity exploitation breaks down.

### 3.1 Square-with-notch spectrum

Base shape: 40×40 square. Cut a notch in one side. Notch depth and
width both vary across a 7×7 matrix:

- `notch_depth_frac ∈ {0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.80}`
  times the side length
- `notch_width_frac ∈ {0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.80}`
  times the side length

49 test shapes total. For each, run the solver on N=10 copies and
record the packing efficiency.

**Expected curve:** efficiency roughly flat at shallow depths (notch
barely matters), dropping sharply in a transition region (medium
notch — solver should discover notch-mating nesting), recovering at
extreme depths (deep U-shapes become qualitatively different and
require cavity filling, category 4 territory).

Plot efficiency vs notch depth. The **dropoff location** tells you
where the solver's concavity exploitation stops working.

Pass criterion: efficiency ≥ 0.70 uniformly across the matrix (no
catastrophic dropoff in any cell). Stretch: ≥ 0.80 in the medium-
notch region where notch-mating should apply.

### 3.2 Rectangular slot spectrum

Base: 60×20 rectangle. Cut a rectangular slot through from one
long edge. Same 7×7 depth×width matrix. Tests the same dropoff but
for asymmetric bases.

---

## Category 4: Islands and archipelagos

**Purpose:** Shapes where the bounding box is much larger than the
actual silhouette area. Ratio `silhouette_area / bbox_area` measures
how much a convex-hull solver leaves on the table. A concave solver
must exploit the gap; a convex solver can't.

### 4.1 Crescents

Two crescents nest beautifully opening-to-opening if rotated.

- **4.1.a Classic crescent.** Circle minus offset circle. Area/bbox
  ratio ≈ 0.35.
- **4.1.b N=10 identical crescents.** Tight packing tests the
  interlock pattern at scale.

Pass criterion: efficiency ≥ 0.80 (well above 0.35 bbox baseline).

### 4.2 Aggressive C-shape

The open channel is wider than the walls. Optimal pack puts another
part INSIDE the C.

- **4.2.a Fat-walled C.** Walls 10mm, channel 20mm wide. Another
  C fits neatly.
- **4.2.b Thin-walled C.** Walls 5mm, channel 30mm wide. Even more
  cavity to exploit.

Pass criterion: efficiency ≥ 0.65.

### 4.3 Ring / annulus

The ultimate island shape — the center is literally empty. Can the
solver put a small part INSIDE the hole? This is qualitatively
different nesting behavior — not edge-to-edge concavity mating,
it's cavity filling.

- **4.3.a Large ring + N small circles.** Ring with inner diameter
  30, outer 60; small circles of diameter 20 each. Pack 1 ring + 5
  circles: the optimum puts one circle inside the ring.
- **4.3.b Stacked rings.** 3 concentric rings of decreasing radius.
  Russian-doll arrangement.

Pass criterion (4.3.a): the solver places at least one circle in
the ring's interior. This is a BINARY test — it's qualitatively
different behavior than edge nesting.

### 4.4 Starfish / asterisk

Multiple concavities radiating outward. The fingers of one starfish
should interleave with another's.

- **4.4.a 5-point starfish.** Two copies should interlock via
  finger-in-gap.
- **4.4.b 8-point asterisk.** More interlock opportunities.

Pass criterion: efficiency ≥ 0.75 on N=2 pair, ≥ 0.70 on N=10.

### 4.5 Spiral

Flat Archimedean spiral with wall thickness. Two spirals can
theoretically interlock. **Probably the hardest possible case for
any solver** — concavity is topologically complex and the valid
nesting requires precise rotation.

Pass criterion: exploratory. Record absolute efficiency without a
threshold. If any solver version hits efficiency ≥ 0.40, that's a
notable milestone.

---

## Category 5: Mixed-scale packing

**Purpose:** Real workloads have parts of wildly different sizes.
The solver needs to fill gaps left by big parts with small parts.

### 5.1 Large C + cavity-filling squares

One large C-shape plus N small squares that exactly fill the C's
cavity. Does the solver discover that the squares go inside the C?

- **5.1.a C with 4 cavity squares.** C channel dimensions match 4
  squares exactly.
- **5.1.b C with 1 perfect square.** One square that fills the
  cavity perfectly.

Pass criterion: solver places the small parts inside the C's
cavity. BINARY test — checks for the "put small thing inside big
thing" behavior.

### 5.2 Circles in circles

One large circle plus smaller circles that achieve the proven
optimal packing of N circles within a circle. Published solutions
exist for N=2 through ~30.

- **5.2.a N=2 small circles inside 1 large.** Optimum is two
  circles with diameter ≈ 0.414 × R.
- **5.2.b N=7 small circles inside 1 large.** Optimum is a central
  circle + 6 surrounding circles (hexagonal pattern inside the
  boundary).
- **5.2.c N=19 small circles inside 1 large.** Complex packing
  with multiple shells.

Pass criterion: efficiency ≥ 0.85 of the known optimum.

### 5.3 Graduated crescents

Same crescent at 100%, 75%, 50%, 25% scale. Smaller ones nest
inside larger ones' concavities, creating a Russian-doll
arrangement.

Pass criterion: BINARY test — each smaller crescent placed inside
the next-larger one's concavity.

---

## Category 6: Adversarial cases designed to break specific heuristics

**Purpose:** Each case targets a specific heuristic weakness. If the
solver passes everything else but fails a category-6 test, you know
exactly which component needs work.

### 6.1 Anti-bottom-left

Two L-shapes where the optimal nesting requires placing one at the
top-right and rotating 180°. A pure bottom-left heuristic will
never find this. Proves the solver isn't just a fancy gravity drop.

Pass criterion: solver places the L at the top-right rotated
position, achieving efficiency > 0.80. A BL-only solver scores
~0.50.

### 6.2 Rotation-sensitive

Asymmetric triangle at ~73° packing. The optimal packing exists
only at an irrational angle (not a multiple of any sensible step
size). Tests whether rotation resolution is sufficient or whether
the solver can do continuous-angle refinement.

- **6.2.a Fixed-angle test.** Triangle whose optimum is at exactly
  73.5°. Solver rotating in 90° steps fails trivially.
- **6.2.b Graduated-angle test.** Same triangle, record efficiency
  vs rotation resolution (90°, 45°, 15°, 5°, 1°, continuous).

Pass criterion (6.2.a): efficiency within 5% of the best possible
at the solver's available rotation resolution.

### 6.3 Order-sensitive

Three shapes A, B, C where packing in order A→B→C wastes 40% but
B→A→C wastes 5%. Tests whether the permutation search is actually
working.

Pass criterion: solver discovers the B→A→C order (or equivalent)
via its smart-shuffle / small-batch / other permutation mechanism.

### 6.4 Identical convex hull

**The litmus test.** Two very different concave shapes that have
the same convex hull. A convex solver treats them identically. A
concave solver should find dramatically different packings for
each.

- **6.4.a Hollow box vs solid box.** Same convex hull. Hollow box
  has an inner cavity that changes everything.
- **6.4.b Plus-sign vs square.** Same convex hull (a square). Plus
  sign tiles perfectly; square grid-packs with waste in any
  orientation.

Pass criterion: efficiency difference between the two shapes >
20%. If C1 and C2 score identically on both, the solver is
secretly wrapping things in convex hulls. **If this passes, the
solver is genuinely doing concave work.** This test alone gates
"are we solving the problem we think we're solving."

---

## Category 7: Real-world regression shapes

**Purpose:** Realistic complexity, no mathematical optimum known,
but a leaderboard tracks the best the solver has ever achieved.

### 7.1 Organic figurine base

Trace the silhouette of a miniature figure's footprint. Typically
an irregular blob with one or two concavities (between legs,
between arm and body).

Pass criterion: leaderboard tracking only. Initial floor set at
the efficiency the current solver produces; future runs must not
regress.

### 7.2 Drone propeller

Long, thin, curved. The bounding box is huge relative to area.
Tests the same failure mode as crescents but with an asymmetric
shape.

Pass criterion: efficiency ≥ 0.50 (any solver that handles
crescents should handle this).

### 7.3 Phone case

Rounded rectangle with camera cutout. Nearly convex but with one
small concavity. Tests whether the solver wastes time on a
concavity too small to exploit.

Pass criterion: efficiency ≥ 0.85. Time within 2x of a
pure-rectangle baseline (no wasted search on the cutout).

### 7.4 Gear

Many teeth, many small concavities. **The optimal packing of two
gears is to mesh them.** If the solver discovers gear-meshing,
it's genuinely doing concavity-aware nesting.

- **7.4.a Spur gear with 20 teeth.** Pair packing.
- **7.4.b Mixed gears.** Different tooth counts, different
  interlocks possible.

Pass criterion (7.4.a): solver meshes the gear pair. Efficiency >
0.85. BINARY test — either meshed or not.

---

## The meta-test

**Purpose:** Quantify solver robustness across random inputs. More
informative for development velocity than any individual test case.

Generate 50 random concave polygons by:

1. Start with a circle of random radius.
2. Apply N random radial perturbations (amplitude ≤ 0.3 × R).
3. Randomly rotate and scale each polygon.

Run the solver on batches of 5, 10, 20 polygons chosen at random.
No known optimum — but run the solver 1000 times with different
random seeds and track:

- **Mean efficiency** — average packing quality
- **Variance in efficiency** — spread across seeds
- **Min efficiency** — worst-case seed
- **Mean time** — runtime distribution
- **Variance in time** — consistency of runtime

Interpretation:
- **Low efficiency variance + high mean** — solver is robust and
  good. Ship it.
- **High efficiency variance + high mean** — solver is lucky
  sometimes. Smart-shuffle is finding good seeds but the core
  algorithm is weak. Needs a better primary strategy.
- **Low variance + low mean** — solver is consistently mediocre.
  Fundamental rethink needed.
- **High time variance** — solver has pathological cases where it
  falls into expensive fallbacks. Profile those specifically.

**This is the development-velocity test.** Run it on every
algorithm change, plot mean + variance vs time, and you have an
objective signal that tells you whether the change was an
improvement or noise.

---

## How this corpus drives C2 development

Each milestone in PR_C2_ORGANIC_LASSO_NESTER_PLAN.md has a test
quality gate derived from this corpus:

- **M0 skeleton** — category 1.1.a (sine jigsaw) must return a
  valid arrangement. No efficiency bar yet.
- **M1 cache + partition** — meta-test on N=10 batches must show
  low variance in partitioning output. Time < 1s per batch.
- **M2 island packing** — category 1 (all subcategories) must
  hit their pass criteria. Category 2.1 (circles) must beat
  grid packing.
- **M3 island placement + spillover** — category 5 (mixed-scale)
  must demonstrate cavity filling on 5.1.a and 5.2.a.
- **M4 performance polish** — meta-test time variance must be
  within 20% of mean (no pathological cases).
- **M5 PR prep** — category 6.4 (identical convex hull) must pass,
  proving C2 is actually doing concave work. Category 7 (real-
  world) efficiency must exceed C1's on the same inputs.

**Category 6.4 is the most important single test in the entire
corpus.** If C2 fails it, the whole organic lasso approach is
secretly still bbox-based and the fork has no reason to exist.
Write this test FIRST when M1 starts.

---

## Implementation notes

- **Fixture storage.** Test fixtures are deterministic generators,
  not STL files. Keeps the test suite self-contained and avoids
  asset-path hassles. Fractal cuts, random polygons, and spirals
  are all parametric.
- **Efficiency measurement.** The `cluster_bbox_perimeter_mm` and
  `cluster_compactness` helpers in `bitmap_test_utils.hpp` are the
  starting point; we'll add a `packing_efficiency(items, params)`
  helper that computes filled area divided by arrangement bbox
  area.
- **Time measurement.** `std::chrono::steady_clock` around the
  `arrange()` call. Record mean + stdev across 5 runs per fixture
  to smooth GC noise.
- **Leaderboard storage.** One JSON file per test case with
  historical best results, timestamps, and commit hashes. Gets
  updated whenever a new best score lands. `docs/c2-leaderboard/`.
- **Visual diffs.** Every test case produces an SVG showing the
  final layout. Side-by-side with the leaderboard's previous best
  makes regressions obvious.

---

## Relationship to the existing test suite

**The 130 C1 tests stay.** They're the correctness gate (no
crashes, no overlaps, deterministic). They're peanuts on quality
but they're the regression floor on behavior.

**This corpus is the C2 development driver.** Every M2 through M5
milestone has a quality gate from a specific category. Early C2
runs will hit < 50% on most categories; that's expected and
useful signal. The goal is to watch the efficiency curves climb
as the algorithm matures.

**C1 never runs this corpus.** C1 is locked at its current
quality. The corpus gates C2 exclusively, so failures on the
corpus are noise on the C1 branch but real regressions on C2.

**The meta-test runs on both.** It's the universal benchmark that
lets C2 be directly compared to C1 across thousands of inputs
without picking favorites.
