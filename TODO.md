# Bitmap Nester — TODO

Branch: `feature/concave-bitmap`
Last session: 2026-04-11 (Hotfix session + War Council round 2)

## Hotfix 2026-04-11 — L-shapes weren't interlocking

Visual regression: user's Sprint 1 binary arranged L-brackets in a loose grid.
War Council round 2 (8 agents) found multiple reasons.

- [x] **BUG A — anchor-distance scoring structurally ties rotations.** Fixed by
  replacing with cluster-bbox-growth scoring (BLF). Rotations now compete on
  actual pack density instead of always losing to rotation 0.
- [x] **BUG B — refine scan locked to the coarse-winner's rotation.** Fixed
  by tracking per-rotation coarse best and refining every rotation, picking
  the global winner.
- [x] **Build filter bug.** `[BitmapRegression]` tag was missing from the
  `build.sh` test allowlist — 8 regression tests silently skipped on every
  run. This is the process root cause for how a broken binary shipped with
  "green tests." Fixed by replacing the allowlist with a `[Bitmap*]` glob.
- [ ] **Scanline OBOE at pixel-aligned edges** (Uruk-hai). Deliberately LEFT
  AS-IS per user: conservative off-by-one (widens shape) is fine for a fuzzy
  rasterizer. Documented in BitmapNester.hpp comment. Not a bug anymore.
- [ ] **Pre-inflation welds notches shut** (Sauron, deferred). `ap.inflation =
  min_obj_distance/2` with pad-per-side sums to the correct `min_obj_distance`
  gap, but at default 6 mm spacing a 10 mm notch can't admit a 10 mm arm.
  Semantic, not a math bug. Revisit if user's retest shows interlocking still
  broken at low spacing — the fix is to halve pad_px in concave mode.

## Sprint 2 — PR Polish (current target)

**New items from War Council round 2** (in addition to the items below):

- [ ] **Per-item bitmap cap** (Aragorn C1). Bed cap is 16M px but each item's
  rot_cache entry is uncapped, enabling a malicious mesh to force 4 rotations
  × N items × up to 128 MB per entry. Add `MAX_ITEM_PIXELS = 1'000'000` and
  reject oversized items before rasterize.
- [ ] **`dilate_bitmap` dst size cap** (Aragorn C2). Check `dwpr * dh` at call
  site before invoking; assert-and-continue rather than allocate blindly.
- [ ] **NaN/Inf guards** on `base_rot`, polygon vertices, `bed_shrink_*`
  (Aragorn W3,W4,N8). Reject non-finite inputs with a log line.
- [ ] **Cap `concave_regions` island count** at ~64 per item (Aragorn W7) and
  total polygon count out of `project_mesh` (Aragorn W6).
- [ ] **Surface concave fallback count** in `ArrangeJob.cpp:747-749`
  (Frodo C1). `"Arranging done (8/10 concave, 2 convex fallback)"` rather
  than silent log-only warning.
- [ ] **Checkbox guardrail** — when `use_concave_shapes` is true, show a yellow
  hint "Spacing > 2 mm disables interlocking" in the ArrangeSettings panel
  (Frodo C2).
- [ ] **INester interface** (Gandalf W3) — pull `BitmapNester` and libnest2d
  dispatch into a strategy selected once at settings-resolution time, rather
  than a branch inside templated `arrange<BedT>`. Enables purge-line cherry
  pick to slot in without touching `Arrange.cpp`.
- [ ] **Kill inflation fallback branch** (Gandalf C4) at `BitmapNester.hpp:218`.
  Make `ap.inflation` authoritative; both nesters honor it identically.
- [ ] **`ObstacleSource` abstraction** (Gandalf W6) — unify excludes, bed
  zones, and the future purge-line regions under one registration path.
- [ ] **`make_concave_ap()` test helper** (Gandalf W5) that populates
  `concave_regions` and 4-rotation by default, deprecating the bare
  `make_ap({0.0})` default.
- [ ] **Smoke test: arrange + install + bbox check** (Gimli C2). Catches ship
  regressions that unit tests miss.
- [ ] **MD5 check BEFORE snapshot copy** (Gimli C3) against
  `$INSTALL_DIR/OrcaSlicer.dll` — detect stale install before write, not
  after.
- [ ] **Strengthen `no_bed_gaps`** — current helper only checks plate index
  contiguity, not actual density. Rename existing helper to
  `plate_indices_contiguous` and add a real `density_above(items, threshold)`
  helper for the regression tests.
- [ ] **Strip 9 stale audit markers** (Gollum 1-7): L1, L3, L7, L9, L12, L13,
  L14, L15. All implemented, all obsolete.
- [ ] **Fix malformed comment** at `BitmapNester.hpp:238` — starts with `\`
  instead of `//` (Gollum 8).
- [ ] **Doc strings on public API**: `rasterize()`, `dilate_bitmap`, `stamp`,
  `collides`, `arrange` (Gollum 16-19).
- [ ] **Tetromino packing benchmark** — test harness + `docs/PACKING_BENCHMARK.md`
  with tier-1 L-tromino 2x3 (100%), tier-2 four-L-tetromino 4x4 (100%),
  tier-3 12-pentomino 6x10 (100%). Tests must measure actual packing density
  and compare against libnest2d convex-hull baseline. Agent-generated in
  this session; verify before shipping.

### Sprint 2 priority tasks (original)

The branch is correct after Sprint 1 but has structural warts that would
bounce on first upstream review. Sprint 2's goal is to land the branch as
the upstream PR target against SoftFever/OrcaSlicer issue #10909.

### Sprint 2 priority tasks

- [ ] **Rebase onto upstream `c948d87`** (1 commit behind — flush_multiplier type fix, shouldn't conflict)
- [ ] **Split `BitmapNester.hpp` into `.hpp` + `.cpp`** — 700-line header-only is the first thing an upstream reviewer will object to (C14)
- [ ] **Remove `#ifdef BITMAP_NESTER_TESTING`** — replace with `friend class BitmapNesterTest;` (C11)
- [ ] **Restructure `goto` topology** — `coarse_hit` / `no_fit_this_rotation` / `done_searching` → small `std::optional`-returning helpers (though T8 already removed `done_searching`)
- [ ] **Strip `L#` audit markers** from comments (L1, L3, L7, L9, L12, L13, L14, L15 — all 2026-04-03 audit shorthand)
- [ ] **Name magic numbers**: `BITMAP_PIXEL_LIMIT`, `DEFAULT_RES_MM`, `COARSE_DIVISOR`, `COARSE_MIN`, `COARSE_MAX`
- [ ] **Fix `sw` shadow** in `dilate_bitmap` (GCC/Clang `-Wshadow`, MSVC C4456)
- [ ] **Add explicit `#include "ClipperUtils.hpp"`** — currently works via the force-included PCH only
- [ ] **Honor `locked_plate`** (fixed obstacle treated as placed on its plate)
- [ ] **Honor `is_extrusion_cali_object`** (zero-inflation branch)
- [ ] **Add `bed_temp` as secondary sort key**
- [ ] **Add 45° rotations** to default rotation set
- [ ] **Map `accuracy` param** to coarse stride
- [ ] **`offset_ex(...)` loop** over all result polygons for exclude inflation (don't drop disjoint pieces)
- [ ] **Apply `-2*EPSILON` shrinkage** parity with libnest2d exclude handling
- [ ] **Bitmap cap fix:** recompute `bw/bh` after `res` update so dimensions stay consistent
- [ ] **Post-centering memory fix:** use stored `inflated_bb` instead of `ExPolygon` copy per item (Legolas's 57MB warning)
- [ ] **NaN/Inf rotation validation** before `rotate()` calls
- [ ] **`BOOST_LOG_TRIVIAL(warning)` on silent UNARRANGED drops** + MAX_PLATES overflow
- [ ] **Improve "can't fit" error** to distinguish too-big / MAX_PLATES / inflation-dropped / keepout-rejected (C10)
- [ ] **GUI tooltip** on "Use actual part shape" checkbox (Frodo)
- [ ] **GUI Reset button desync fix** in `GLCanvas3D.cpp` (Frodo)
- [ ] **GUI log line** distinguishing concave vs default arrange mode (Frodo)
- [ ] **Purge line avoidance cherry-pick:** new `ConfigOptionPoints purge_line_region` key mirroring upstream `bd066e7f` wrapping-detection pattern. Zero `BitmapNester.hpp` changes.
- [ ] **Build harness:** patch `build.sh` to run `cmake --install` before snapshot step so install-dir staleness doesn't keep failing MD5 check. Or document the manual workflow in `docs/BUILD.md`.

### Sprint 2 exit criteria

- `BitmapNester.hpp` is a thin header; `.cpp` holds impl
- No `#ifdef BITMAP_NESTER_TESTING` anywhere
- `-Wshadow` + `/W4` clean
- Purge-line key visible in Tab UI; both nesters honor it
- ~100 tests total green (Sprint 1's 109 + ~15 new Sprint 2 tests)
- Upstream PR title: `feat(arrange): bitmap-based concave nester + purge line region`
- Upstream target: SoftFever/OrcaSlicer issue #10909

## Sprint 3 — Workflow features (after upstream merge)

Deferred until Sprint 2 lands. Split into three follow-up PRs:

- Material grouping pre-partition + nonprefered soft-penalty bitmap + BBL bailout escape hatch (~310 LOC)
- v2 arrange modes port: Keep Plates + This Plate + Stragglers + three-way repack fallback (~430 LOC, archaeology-validated)
- Consolidation post-pass from v1 simplified (~150 LOC, fixes known bug #1)

## Known Bugs (from 2026-04-03 audit, still open)

- [ ] **Plate consolidation overlap** — parts on plate 3 with empty plates 1-2 get condensed onto plate 1 but overlap each other. The multi-plate consolidation path isn't re-running placement after moving items between plates. Sprint 3 consolidation post-pass will fix this.
- [ ] **Keepout zone padding** — excludes may need additional padding around wipe tower / cali zones. Sprint 1 C15 dual-bitmap fix already addresses the convex-vs-concave mismatch that caused phantom "part outside plate" warnings, but the raw pad distance itself may still be tight.
- [ ] **Plate ID global/local mismatch** — likely to surface when the plate grid layout changes (e.g. 4→5 plates). Plate IDs may be indexed against the old layout.

## Post-Sprint-3 fork: FFT-correlation nester

Reference design: `docs/FFT_NESTER_IDEA.md`. Pin for after the upstream PR
lands. Two qualitative wins over the current scored scan: (a) FFT-based
placement search is independent of part vertex complexity and asymptotically
better at high resolution + many rotations, (b) random-permutation outer
loop with parallel restarts has well-known density wins over single-pass
greedy. Estimated effort: 5 phases starting with a 64x64 hand-FFT
correctness test, ending with head-to-head benchmarks against the scored
scan on the tetromino harness. Branch name: `feature/concave-bitmap-fft`.

## Deferred indefinitely (cut from 3-sprint plan per Gandalf)

- Full scored placement extensions (same-size alignment bonus, big/small orbit, height clustering) — anchor bias gets ~80% of visible win for 10% of the LOC. Revisit only if a user files a specific complaint anchor bias can't explain.
- `parallel` arrange param — weak effect, threading complexity
- Extruder-change minimization — transitively handled by Sprint 3 material grouping
- Silhouette pipeline hardening (worker-thread move, Clipper hang guard, morphological-close fix, `throw_on_cancel`) — separate follow-up session focused on silhouette pipeline as a unit

## Sprint 1 — Done (2026-04-10)

- [x] Doc reorg: split in-tree `CLAUDE.md` into focused guides (`be321a18`)
- [x] Fix stress assertion contradicting overflow behavior (`2e2df815`)
- [x] Silhouette multi-island, cache key, silent-failure fixes (C5, C6, C7) (`efb3997e`)
- [x] Multi-island rasterize, `align_to_y_axis` composition, stopcondition in centering (C2, C13) (`a53282d3`)
- [x] End-to-end mesh fixture tests (6 new tests, `[BitmapMesh]` tag) (`f76440de`)
- [x] C15 dual-bitmap: convex-hull vs plate_excludes agrees with `PartPlate::check_outside` (C12, C15) (`56607f91`)
- [x] Anchor-biased scored placement for wipe tower clustering (C1) (`167eab35a7`)
- [x] 109/109 tests green, 1082 assertions
- [x] End-of-sprint snapshot: `orca/builds/clean/20260410-210803-167eab35a7-sprint1-complete/`
