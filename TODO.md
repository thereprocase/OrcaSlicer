# Bitmap Nester MVP — TODO

Branch: `feature/concave-bitmap` @ `8bc11dd`
Last session: 2026-04-03

## Must-do

- [ ] Verify stress tests pass — file committed at `beec87f` but full run unconfirmed
- [ ] Rebase onto upstream `c948d87` (flush_multiplier type fix, shouldn't conflict)
- [ ] Run full test suite after rebase — all bitmap tags including `[BitmapStress]`
- [ ] Snapshot after rebase
- [ ] Clean up `orca/rainbow-test/` — reports migrated to `orca/dev/sessions/2026-04-03/`, originals can be deleted

## Should-do

- [ ] Lower coarse stride minimum from 8→4 for tiny items (convergent finding from Pink + Green — escalates to must-do if stress tests fail on thin items -- can we have a minimum pad that fixes this? or some sort of thickness aware padding to buff up thin sections?)
- [ ] Verify `on_packed` callback has downstream consumers (added but nobody checked if anything uses it)
- [ ] Note `goto` topology and `#ifdef BITMAP_NESTER_TESTING` in eventual PR description (preempt reviewer questions)
- [ ] Note exclude inflation `.front()` simplification in PR description (Pink finding #2)

## Build harness (local tooling, not committed — upstream doesn't ship one)

- [x] `orca/build.sh` — updated with `test` target, `--preflight`, per-worktree snapshots, full-install from `build/OrcaSlicer/`
- [x] `orca/check-build.sh` — pre-flight validation with worktree shortcuts, smart junction checks
- These live at `orca/` level, shared across worktrees. Not part of any branch.

## Known Bugs (investigate)

- [ ] **Plate consolidation overlap** — parts on plate 3 with empty plates 1-2 get condensed onto plate 1 but overlap each other. The multi-plate consolidation path isn't re-running placement after moving items between plates.
- [ ] **Keepout zone padding** — exclude zones (wipe tower, calibration areas) need additional padding. Items placed right at the edge of a keepout boundary.
- [ ] **Plate ID global/local mismatch** — likely to surface when the plate grid layout changes (e.g. 4→5 plates). Plate IDs may be indexed against the old layout, causing items to land on wrong plates or overlap after a layout reflow.

## Deferred (post-MVP)

- [ ] Material grouping (`allow_multi_materials_on_same_plate`) — complex, needs architectural changes
- [ ] Inflation double-counting — conservative (extra spacing), matches convex nester behavior
- [ ] Non-preferred regions (`nonprefered_regions`) — cosmetic only
- [ ] Y-axis alignment (`align_to_y_axis`) — aesthetic only
- [ ] Parallel placement — single-threaded is fast enough for now

## Done (this session)

- [x] seq_print bailout → `80094dd`
- [x] Exclude-safe centering pass → `7130a02`
- [x] on_packed callback → `8b3d618`
- [x] Test expectation update for centering → `a7ea1a6`
- [x] Gap tests (oversized items, MAX_PLATES cap) → `35998bd`
- [x] Stress test suite (11 adversarial cases) → `beec87f`
- [x] Build guide → `8bc11dd` (`docs/BUILD.md`)
- [x] 6-agent audit — zero bugs found across 4 independent reviews
- [x] Full-install snapshot → `orca/builds/clean/20260403-093123-beec87f137-mvp-fixes/`
- [x] Build harness upgrades (build.sh + check-build.sh)
- [x] Name scrub across all report files
- [x] Dev notes repo → `orca/dev/` pushed to `thereprocase/orca-dev` (private)
- [x] Session reports, worktree map, CURRENT.md, CLAUDE.md for AI ramp-up
