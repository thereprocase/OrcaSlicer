# Commit Build/Run Status — feature/concave-arrange

Last updated: 2026-03-29

## Status Key
- GOOD = built AND launched AND tested by Repro
- BUILDS = compiled successfully, not confirmed to launch
- UNTESTED = not built from this specific commit
- CRASHED = launches then crashes, or doesn't launch at all

## Commits (newest first)

| Hash | Status | Notes |
|------|--------|-------|
| `d4887e407a` | CRASHED | Comment-only changes. Crash was stale .obj from b810ec (GLCanvas3D.obj had dead version() call baked in). Clean rebuild may fix. |
| `a223aca341` | UNTESTED | Just adds TODO.md. Should be identical binary to 3f401c. |
| `3f401c638b` | **GOOD** | Last confirmed working build. Good 3D density, no boundary errors. Snapshot: `3f401c638b_density_fix/` |
| `cbef94e14d` | GOOD | Fixed 3D boundary overshoot. Per-rotation bounds check. |
| `92189a7b8e` | GOOD | Megapatch. Non-rect beds, fill-bed, 36 cap. 3D density regressed (fixed in 3f401c). |
| `00f06e4451` | GOOD | Removed +1 from bed dimensions. Fixed right-edge overshoot. |
| `c29871f8d3` | GOOD | Min/max bitmap centering fix. 2-pixel bed margin (later removed). |
| `11fdbda856` | GOOD | Min/max pre-filter, best-rotation, batch rejects. |
| `982984dcbf` | **GOOD** | Major milestone. Multi-material, rotation step, UI cleanup, debug log. Repro soak-tested this extensively. |
| `3da63b37e0` | GOOD | Exclusion zones at all Z heights. |
| `7f6e93ad47` | GOOD | 3D compaction, 3D overlap check. "FUCK YEAH BUDDY" |
| `2afc4b3ccc` | GOOD | 3D slice coordinate alignment. Fixed yeeted parts. |
| `cbe7c414f9` | GOOD | Code review fixes (6 issues). |
| `f2342fa346` | GOOD | 3D bitmap collapse + cached skyline. |
| `81ebffe5bf` | GOOD | Reverted pipeline placement (vector invalidation race). |
| `a53852bb0b` | BUILDS | 3D slice collapse. Not independently tested. |
| `282ecb9cc6` | CRASHED | Pipeline placement had vector invalidation race. Parts fell off bed. |
| `59df5b687d` | GOOD | Parallel TBB compaction. |
| `4571b6cf9c` | GOOD | Adaptive-stride compaction, overlap safety check, UI fix. |
| `b86a60acae` | GOOD | Compaction UI controls. |
| `d777b97350` | GOOD | Configurable compaction strategies. |
| `20ec0bcf28` | GOOD | Coarse bitmap scan compaction. |
| `f3a015961b` | GOOD | First compaction via bitmap scanning. |
| `d2a6a5f15d` | GOOD | Compaction sweep. |
| `77c43d67f7` | GOOD | Skyline placement. Major speed improvement. |
| `c593306da6` | GOOD | 3-tier placement search. |
| `f3345a59c1` | GOOD | Plate-centric filling. |
| `438c7eac9b` | GOOD | 3D boundary fixes. |
| `def6f52390` | GOOD | 3D-aware nesting via Z-slices. |
| `ab7043ed69` | GOOD | 3-phase batch with progress. |
| `636c745935` | GOOD | Bitmap rotation, center-out placement. |
| `d224a57cd2` | GOOD | Triangle rasterization. First working concave build. |

## End-of-session note (2026-03-29 evening)
ALL builds stopped launching — including previously confirmed working ones.
System showed signs of instability (explorer hangs, .NET runtime failures,
Excel crash). Likely needs a reboot. If old builds launch after reboot,
the code is fine and we just had stale .obj contamination + system issues.
