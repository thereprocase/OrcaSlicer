# Tetromino Packing Benchmark

A reproducible efficiency harness for the concave BitmapNester, using the seven
tetrominoes (I, O, T, S, Z, J, L) as probe shapes.

## Why tetrominoes

Every tetromino is exactly 4 unit squares, so piece area is always known and
integer. Entire rectangles can be cut into tetrominoes and reassembled with
density 1.0, which gives every scenario a rigorous upper bound. Five of the
seven pieces (T, S, Z, J, L) are concave, so their convex hulls waste area —
that waste is the number the bitmap nester has to beat.

## Conventions

- **Unit cell**: 10 mm per square. A unit tetromino is 4 squares of 10 mm x 10 mm = 400 mm^2.
- **Coordinate origin**: each piece is built with its local bbox anchored at (0, 0).
  Rotation around the local origin is handled by the nester and `allowed_rotations`.
- **Efficiency metric**: `density = total_piece_area / cluster_bbox_area`, where the
  cluster bbox is the axis-aligned bounding box of every placed, transformed piece
  on bed 0. Values above 1.0 are impossible; values near 1.0 mean tight nesting.
- **Baseline**: same items fed through `arrangement::arrange(...)` with
  `use_concave_shapes=false`, which routes to libnest2d and packs convex hulls. The
  bitmap nester must beat this density by a scenario-specific margin.
- **Reproducibility**: every random scatter uses a fixed `std::mt19937` seed so the
  benchmark is deterministic across runs.

## Piece definitions (mm, unit=10)

All seven pieces in their canonical rotation. Each vertex list is CCW. Vertices
can be plugged directly into `ExPolygon(Points{Point(scaled<coord_t>(x), ...)})`.

### I — 1 x 4 bar
```
(0,0) (40,0) (40,10) (0,10)
```
Convex. Hull area = piece area = 400 mm^2. No interlocking gain.

### O — 2 x 2 square
```
(0,0) (20,0) (20,20) (0,20)
```
Convex. Hull area = 400 mm^2. Baseline nester handles this fine.

### T — T-piece
```
(0,10) (10,10) (10,0) (20,0) (20,10) (30,10) (30,20) (0,20)
```
Concave. Hull is a 30 x 20 trapezoid-ish shape. Hull area 600 mm^2; piece 400.
Hull waste: 33%.

### S — S-piece
```
(0,0) (20,0) (20,10) (30,10) (30,20) (10,20) (10,10) (0,10)
```
Concave. Hull area = 500 mm^2 (30 x 20 minus two corner triangles).
Hull waste: ~20%.

### Z — Z-piece (mirror of S)
```
(0,10) (10,10) (10,0) (30,0) (30,10) (20,10) (20,20) (0,20)
```
Concave. Same 500 mm^2 hull, same ~20% hull waste.

### J — J-piece
```
(0,0) (10,0) (10,10) (30,10) (30,20) (0,20)
```
Concave. Hull area = 500 mm^2 (a trapezoid). Hull waste: ~20%.

### L — L-piece (mirror of J)
```
(0,10) (20,10) (20,0) (30,0) (30,20) (0,20)
```
Concave. Hull area = 500 mm^2. Hull waste: ~20%.

## Scenarios

Each scenario lists:
- **Pieces**: the multiset of tetrominoes used.
- **Ideal bbox**: the rectangle those pieces can tile with zero waste.
- **Ideal density**: always 1.0 by construction.
- **Convex baseline**: rough density the libnest2d (convex-hull) path can achieve,
  estimated from summed hull areas divided by the smallest convex-hull packing
  rectangle. This is a pessimistic ceiling used to set the "must beat" threshold.
- **Bitmap threshold**: density the concave nester must exceed. Set well above
  the convex baseline but well below the ideal — the raster grid and scored scan
  prevent anyone from actually hitting 1.0.

### Scenario 1 — "Tetris 4x4" (one of each, rotated)

Not all 7 tetrominoes fit into a single rectangle (I, T, S, Z, J, L, O have a
combined area of 28 squares, but a 4x7 rectangle with one of each is impossible
due to coloring-parity: S+Z+T together cover an odd number of black squares on
a checkerboard — there is no tiling). Instead this scenario uses **two copies
of each piece** (14 pieces, 56 squares) and an **8 x 7 = 56 square** target
(80 mm x 70 mm). This target is known to be tileable by duplicated tetrominoes.

- **Pieces**: 2x each of I, O, T, S, Z, J, L (14 pieces, 5600 mm^2)
- **Ideal bbox**: 80 mm x 70 mm = 5600 mm^2, density 1.0
- **Convex baseline**: ~0.70 (hulls of T/S/Z/J/L waste 20-33% and libnest2d
  won't tile them perfectly either — empirically we expect ~0.60-0.75)
- **Bitmap threshold**: density >= 0.85
- **Bitmap must beat baseline by**: >= 10 percentage points

### Scenario 2 — "Rectangular 4x4" (I only)

Four I-tetrominoes stacked tile a 4 x 4 square exactly (40 mm x 40 mm).

- **Pieces**: 4x I (1600 mm^2)
- **Ideal bbox**: 40 mm x 40 mm = 1600 mm^2, density 1.0
- **Convex baseline**: ~1.0 — I-pieces are convex, libnest2d will also tile this.
- **Bitmap threshold**: density >= 0.95
- **Bitmap must beat baseline by**: not applicable — this is a sanity check.
  Both nesters should hit > 0.95.

### Scenario 3 — "L-only 2x4"

**CORRECTION (2026-04-11):** the original spec for this scenario claimed
two L-tetrominoes of the same chirality tile a 4x2 rectangle. They do not.
The classical 4x2 tiling requires an L + J pair (mirror images of each
other). With pure rotations only, two same-chirality L's CANNOT tile any
rectangle — the best they can achieve is two bboxes side-by-side for
density 800 / 1200 = **0.6667 exactly**. The test threshold is set just
below that floor; anything less indicates the nester regressed to a
non-touching grid (which is the original loose-grid bug this benchmark
exists to catch).

For a real "interlock or fail" test, scenario 4 (S+Z) is the better
discriminator — its 4x4 tiling is achievable with pure rotations.

- **Pieces**: 2x L (800 mm^2)
- **Achievable bbox** (same chirality): 60 mm x 20 mm = 1200 mm^2, density 0.6667
- **Convex baseline**: same — libnest2d also produces side-by-side hulls
- **Bitmap threshold**: density >= 0.65 (regression guard, not optimality test)

### Scenario 4 — "S-and-Z"

A single S and a single Z cannot tile a rectangle, but 2x S + 2x Z can tile a
4 x 4 square (40 mm x 40 mm). Both pieces have the worst convex-hull overhead
that can be mitigated by interlocking.

- **Pieces**: 2x S, 2x Z (1600 mm^2)
- **Ideal bbox**: 40 mm x 40 mm = 1600 mm^2, density 1.0
- **Convex baseline**: ~0.65 — hulls of S/Z don't interlock at all in libnest2d.
- **Bitmap threshold**: density >= 0.85
- **Bitmap must beat baseline by**: >= 15 percentage points

### Scenario 5 — "Mixed 40-piece set"

A 10 x 16 square target (100 mm x 160 mm) cut into 40 tetrominoes. Weighted
toward concave pieces so the bitmap nester has room to shine.

- **Pieces**: 6x I, 6x O, 6x T, 5x S, 5x Z, 6x J, 6x L (40 pieces, 16000 mm^2)
- **Ideal bbox**: 100 mm x 160 mm = 16000 mm^2, density 1.0
- **Convex baseline**: ~0.70
- **Bitmap threshold**: density >= 0.80
- **Bitmap must beat baseline by**: >= 8 percentage points

## Efficiency thresholds (summary)

| Scenario | Pieces | Ideal area (mm^2) | Bitmap min density | Convex baseline min beat |
|---|---|---:|---:|---:|
| 1 Tetris 4x4 | 14 | 5600 | 0.85 | +0.10 |
| 2 Rectangular 4x4 | 4 | 1600 | 0.95 | n/a (sanity) |
| 3 L-only 2x4 | 2 | 800 | 0.90 | +0.15 |
| 4 S-and-Z | 4 | 1600 | 0.85 | +0.15 |
| 5 Mixed 40-piece | 40 | 16000 | 0.80 | +0.08 |

The thresholds are starting points. Real measured numbers should replace them
after the first green run — the test file's `REQUIRE(density >= X)` lines are
the place to lock in what the nester actually achieves.

## How the test measures

1. Build the piece multiset for the scenario.
2. For each piece, construct an `ArrangePolygon`:
   - `ap.poly = tetromino_expolygon`
   - `ap.allowed_rotations = {0, PI/2, PI, 3*PI/2}` (so the nester can rotate)
   - `ap.translation = {random x, random y}` inside a 300 mm x 300 mm scatter box,
     seeded with `std::mt19937(0xC0FFEE)`.
3. Call `BitmapNester::arrange(items, excludes, bed, params)` directly. Bed is
   generous (500 mm x 500 mm) so the nester is not bed-limited.
4. Params: `no_shrink_params()` with `p.min_obj_distance = 0` and
   `p.allow_rotations = true`.
5. After arrange, iterate placed items, call `it.transformed_poly()`, union the
   bounding boxes. Compute density = sum(piece_area) / cluster_bbox_area.
6. For the baseline, build a fresh item list the same way, then call
   `arrangement::arrange(items, excludes, bed, params)` with
   `params.use_concave_shapes = false`. Measure its density the same way.
7. Assert: `density_bitmap >= scenario_threshold` AND
   `density_bitmap >= density_baseline + required_margin`.

## Notes on raster tolerance

At the default raster resolution (~0.5 mm/pixel inside the nester) a single
tetromino spans 20x20 to 20x60 pixels. One pixel of error per edge gives up to
~5% area error on the smallest pieces. Thresholds above are set with 5-10%
headroom so raster quantization alone won't fail a test.
