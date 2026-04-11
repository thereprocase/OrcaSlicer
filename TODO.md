# Bitmap Nester — TODO

Branch: `feature/concave-bitmap` @ `167eab35a7`
Last session: 2026-04-10 (Sprint 1 complete)

## Sprint 2 — PR Polish (current target)

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
