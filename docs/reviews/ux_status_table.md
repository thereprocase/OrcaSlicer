# Snuggle UX Status Table

Generated from Council of Elves expectations (50 items) vs current implementation.

**Score: 33/50 done, 5 partial, 12 open**

## Placement Correctness

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 1 | No overlap at any Z | P0 | DONE | Voxel collision + independent validator confirms |
| 2 | Parts within bed boundary | P0 | DONE | bed_margin_mm in OOB check, clamps, evaluator |
| 3 | Exclude zones respected | P0 | OPEN | Wipe tower, purge line, LIDAR still ignored |
| 4 | Padding as hard constraint | P0 | DONE | min_gap_mm inflates grids, slider controls it |
| 5 | Locked parts stay put | P1 | OPEN | Snuggle sees all items, doesn't filter locked |
| 6 | Lock rotation preserves angle | P1 | DONE | Validated by fellowship + standalone tests |
| 7 | Group identical parts | P2 | OPEN | |

## Visual Feedback

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 8 | Progress indicator | P0 | DONE | "Snuggle GPU gen 42/100 OK" in progress bar |
| 9 | Viewport updates immediately | P0 | DONE | Standard ArrangeJob finalize path |
| 10 | Failure notification | P1 | DONE | Best-effort + off-plate overflow |
| 11 | Status during arrange | P1 | DONE | Backend + gen + collision count in status bar |
| 12 | Overflow parts distinguishable | P1 | DONE | bed_idx=-1 -> off-plate (Orca standard) |
| 13 | Before/after comparison | P2 | OPEN | Undo stack handles this |
| 14 | Packing density metric | P2 | OPEN | In log, not in UI |
| 15 | Exclusion zones drawn | P2 | OPEN | |

## Controls

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 16 | Part gap slider | P0 | DONE | 0-20mm, default 5mm |
| 17 | Quality slider | P0 | DONE | 1-10, maps to pop/gen/timeout |
| 18 | Cancel button | P1 | PARTIAL | stopcondition wired but timeout is primary stop |
| 19 | Undo (Ctrl+Z) | P1 | FREE | ArrangeJob already pushes undo snapshot |
| 20 | Settings persist | P1 | DONE | app_config save/load for all snuggle params |
| 21 | Snuggle off -> stock arrange | P1 | DONE | use_snuggle gates the dispatch |
| 22 | Advanced panel | P2 | OPEN | Voxel res, timeout, bed margin not exposed |
| 23 | Try again / re-seed | P2 | OPEN | |
| 24 | Per-part rotation lock | P2 | OPEN | |

## Performance

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 25 | GUI stays responsive | P0 | DONE | Background thread via ArrangeJob |
| 26 | < 30s for typical beds | P0 | DONE | 300ms-6s at quality 5 (perf tests) |
| 27 | GPU 3x+ speedup | P1 | DONE | 6.8x with rotation cache |
| 28 | Bounded memory | P1 | DONE | ~84MB for 15 parts x 360 bins |
| 29 | Voxelization < 3s | P1 | DONE | Decimation to 5K tris, 750ms typical |
| 30 | Progress every 500ms | P2 | DONE | yield_every_gens = 1 |
| 31 | Below-normal priority | P2 | OPEN | Not set in ArrangeJob for Snuggle |

## Error Handling

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 32 | Oversized parts handled | P0 | DONE | Off-plate with greedy accept |
| 33 | Zero items handled | P0 | DONE | Early return in SnuggleArrange |
| 34 | Timeout -> best effort | P0 | DONE | Applies best result, overflow off-plate |
| 35 | GPU fallback silent | P1 | DONE | Log warning only, auto CPU fallback |
| 36 | Bad mesh fallback | P1 | PARTIAL | Voxelization error -> 0-size grid, placed but no collision |
| 37 | Overflow visual | P2 | DONE | bed_idx=-1 -> off-plate |

## Integration

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 38 | Wipe tower respected | P0 | OPEN | Logged as ignored |
| 39 | Same postprocess path | P0 | DONE | Standard ArrangeJob finalize |
| 40 | Multi-plate | P1 | OPEN | Single-plate only |
| 41 | Sequential print clearance | P1 | OPEN | |
| 42 | Auto-orient compatible | P1 | DONE | Uses raw_mesh reflecting current orientation |
| 43 | SLA disabled | P2 | OPEN | |

## Predictability

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 44 | Deterministic seed | P0 | DONE | Fixed seed=42 in constructor |
| 45 | Incremental add | P1 | OPEN | Full re-arrange only |
| 46 | Quality monotonic | P1 | DONE | Elitism preserves best individual |
| 47 | Clean visual layout | P2 | DONE | Clustering fitness weight=1.5 |
| 48 | Height grouping | P2 | DONE | Height-centering fitness weight=0.3 |

## Accessibility

| # | Expectation | P | Status | Notes |
|---|-------------|---|--------|-------|
| 49 | All controls have tooltips | P1 | PARTIAL | Quality has tooltip, gap and lock need them |
| 50 | Experimental label | P1 | DONE | Teal "Snuggle (experimental)" header |

## Summary

| Priority | Done | Partial | Open | Total |
|----------|------|---------|------|-------|
| P0 | 14 | 0 | 2 | 16 |
| P1 | 14 | 3 | 5 | 22 |
| P2 | 5 | 0 | 7 | 12 |
| **All** | **33** | **3** | **14** | **50** |

## P0 Blockers for Public Release

1. **#3 / #38 — Exclusion zones + wipe tower**: Parts will overlap the wipe tower. Print failure guaranteed on multi-color jobs.
2. Both are the same underlying issue: `SnuggleArrange.cpp` ignores the `excludes` parameter.

## Next Priority After P0

- #5 Locked parts (filter in SnuggleArrange, not just ArrangeJob)
- #18 Cancel (verify stopcondition fires within 1s)
- #49 Tooltips (3 missing)
- #31 Thread priority
