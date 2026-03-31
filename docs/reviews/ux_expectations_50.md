# Snuggle 3D Nester: 50 Core UI/UX Expectations

Organized by category. Priority levels:
- **[P0]** Must have -- broken without it
- **[P1]** Should have -- users will complain
- **[P2]** Nice to have -- polish

Current Snuggle state: genetic algorithm with voxel collision (CPU + GPU), single-plate, no exclusion zone handling, has padding/quality/lock-rotation controls, no undo, no preview, no multi-plate overflow.

---

## 1. Placement Correctness

**1.** [P0] After arrange completes, no two parts overlap at any Z height. Specifically: for every pair of parts, the intersection of their voxelized volumes (at the configured voxel resolution) is empty after applying the minimum gap.

**2.** [P0] Every part's bounding box (including padding) lies entirely within the printable bed region. No vertex of any part's convex hull at Z=0 may exceed the bed boundary minus `bed_margin_mm`.

**3.** [P0] Exclusion zones (wipe tower footprint, purge line strip, LIDAR calibration region) are treated as immovable occupied rectangles. No part may overlap these zones after arrange. *Currently violated -- Snuggle logs "excluded items will be IGNORED".*

**4.** [P0] The configured `snuggle_padding_mm` gap is enforced as a hard constraint, not a soft fitness penalty. A 5mm gap setting means no two parts' surfaces are closer than 5mm in XY at any Z height.

**5.** [P1] Parts placed by the user before arrange and marked as "locked" (via the plate lock mechanism) remain at their exact position and rotation. Other parts arrange around them.

**6.** [P1] When `snuggle_lock_rotation` is enabled, every part's Z-rotation after arrange equals its Z-rotation before arrange, within 0.01 degrees.

**7.** [P2] Parts with identical geometry (multiple instances of the same STL) may optionally be grouped to maintain consistent orientation, so prints look uniform.

---

## 2. Visual Feedback

**8.** [P0] While arrange is running, a progress indicator is visible. At minimum: a spinner or progress bar. Ideal: percentage + generation count + best fitness, updated every 0.5s or fewer.

**9.** [P0] When arrange finishes, the 3D viewport updates to show the new part positions immediately. No stale ghost positions from the previous arrangement.

**10.** [P1] If arrange fails entirely (all parts too large, timeout with no feasible solution), a notification banner appears stating the failure reason in user-facing language ("Parts do not fit on the bed" not "infeasible individual count > 0").

**11.** [P1] During arrange, the status bar or notification area shows elapsed time and a brief status string ("Snuggle: generation 47/100, best gap 5.0mm").

**12.** [P1] After arrange completes, parts that could not be placed (too large, no room) are visually distinguishable -- red tint, warning icon, or moved to an "overflow" area -- rather than silently dropped or left at origin.

**13.** [P2] A "before/after" comparison: the user can toggle or split-view the arrangement vs. the previous layout to evaluate whether the result is better.

**14.** [P2] The notification after arrange reports a packing density metric ("72% bed utilization") so the user can judge quality without eyeballing it.

**15.** [P2] Exclusion zones (wipe tower, purge strip, LIDAR region) are drawn as shaded rectangles on the bed during arrange preview, so the user sees the actual available space.

---

## 3. Controls

**16.** [P0] "Part gap" slider range (currently 0-20mm) covers the useful range. 0mm minimum must still enforce voxel-resolution clearance (no literal overlap). Default 5mm is reasonable for most FDM printers.

**17.** [P0] "Quality" slider (currently 1-10, integer) maps to generation count and population size in a way the user doesn't need to understand. Higher number = visibly better packing at the cost of proportionally longer runtime. The mapping should be documented internally (e.g., quality 1 = 20 generations, quality 10 = 200 generations).

**18.** [P1] A "Cancel" button is accessible during arrange execution. Pressing it stops the genetic algorithm within 1 second and restores the pre-arrange layout (no partial result applied).

**19.** [P1] Ctrl+Z (undo) after arrange restores the previous part positions. This requires the arrange operation to push a snapshot onto OrcaSlicer's undo stack before modifying positions.

**20.** [P1] The "Use Snuggle 3D nester" checkbox persists across sessions (currently implemented via app_config -- verified in code).

**21.** [P1] When the user toggles Snuggle off, the stock OrcaSlicer arrange (libnest2d NFP) runs instead, with no residual Snuggle state affecting it. The "Defaults" button correctly resets all four Snuggle parameters.

**22.** [P2] An "advanced" expandable section exposing: voxel resolution (default 2mm), timeout (default 30s), bed margin (default 4mm). Power users who print small parts at 0.1mm nozzle need finer voxels; large-format users need longer timeouts.

**23.** [P2] A "seed" or "try again" button that re-runs arrange with a different random seed, producing a different (possibly better) layout. The genetic algorithm is stochastic; users should be able to explore alternatives.

**24.** [P2] Per-part rotation lock: right-click a part to lock its rotation individually, separate from the global `snuggle_lock_rotation` toggle. Useful for parts with orientation-dependent print quality.

---

## 4. Performance

**25.** [P0] The GUI remains responsive (repaints, mouse events, cancel button) during arrange. Arrange must run on a background thread. Currently ArrangeJob inherits from Job, which runs on a worker thread -- verify Snuggle doesn't call OpenGL or wx from the worker.

**26.** [P0] For a typical bed (256x256mm, 5-15 parts), arrange completes within the configured timeout (default 30s). Quality=5 should finish in under 15 seconds on a mid-range CPU (Ryzen 5600 / i5-12400).

**27.** [P1] GPU-accelerated collision detection (when available) provides at least a 3x speedup over CPU-only. If GPU init fails, fallback to CPU is silent and automatic -- no error dialog, just a log line.

**28.** [P1] Memory usage stays bounded. Population of 512 individuals x 50 parts x 3 floats = ~300KB. Voxel grids at 2mm resolution for a 256mm bed = 128x128 per layer. Total voxel memory should not exceed 100MB for typical part counts.

**29.** [P1] Voxelization (mesh decimation + grid fill) completes in under 3 seconds for up to 20 parts. The current decimation target (5000 tris) is reasonable; verify it doesn't dominate runtime for already-low-poly meshes.

**30.** [P2] Progress callbacks fire at least every 500ms so the UI can update generation count and fitness. The current `yield_every_gens = 1` in NesterConfig enables this, but the callback must not block on UI thread synchronization.

**31.** [P2] Build priority: the arrange thread runs at below-normal CPU priority so it doesn't starve the UI or other background tasks (slicing, preview generation).

---

## 5. Error Handling

**32.** [P0] If a single part is larger than the printable bed area (bounding box exceeds bed minus margins), arrange reports which part(s) don't fit by name, and either skips them or places them centered with a warning. It must not crash or hang.

**33.** [P0] If zero arrangeable items are selected (all locked or all unprintable), arrange shows a notification and returns immediately. *Currently implemented in ArrangeJob for the stock path -- verify Snuggle also handles this.*

**34.** [P0] If the genetic algorithm hits the timeout (default 30s) without finding a feasible solution, it returns the best infeasible solution found and warns: "Could not fit all parts -- N parts overlap. Try reducing part count or increasing quality." It must not return an empty or zero-position result.

**35.** [P1] If GPU collision init fails (missing OpenCL/CUDA, driver too old, out of VRAM), Snuggle falls back to CPU collision without user intervention. A single log-level warning is emitted, not a modal dialog.

**36.** [P1] If the model contains non-manifold or degenerate meshes that can't be voxelized, Snuggle falls back to 2D convex hull for that specific part and proceeds. The user sees a warning naming the problem part.

**37.** [P2] If only some parts fit and others don't, the result clearly separates placed parts (on bed, green) from overflow parts (off bed or flagged red). Ideally, overflow parts move to a second plate. *Multi-plate is currently missing.*

---

## 6. Integration

**38.** [P0] Snuggle respects the wipe tower position. When multi-color printing is active, the wipe tower's footprint (from `get_wipe_tower_info`) is an immovable exclusion rectangle. *Currently not implemented -- the code explicitly logs that excludes are ignored.*

**39.** [P0] Snuggle's output is consumed by the same `postprocess_arrange_polygon` / `apply()` path as stock arrange. Translations and rotations are in the same coordinate system (scaled integers, bed-relative). Verify by checking that slicing after Snuggle arrange produces correct G-code positions.

**40.** [P1] Snuggle works correctly with OrcaSlicer's multi-plate system. At minimum: arrange on the current plate only, leaving other plates' parts untouched. Ideal: overflow parts auto-assigned to plate 2+ via `bed_idx`.

**41.** [P1] Sequential print mode: when "Print sequence = By object" is active, Snuggle must enforce the toolhead clearance envelope between parts (extruder + gantry collision box), not just the part geometry. *The stock arrange has `is_seq_print` handling -- Snuggle needs equivalent logic.*

**42.** [P1] Snuggle result is compatible with the "auto-orient" feature. If the user runs auto-orient then arrange, Snuggle uses the oriented mesh. If arrange then orient, the orient step doesn't invalidate the arrangement. (In practice, users should orient first.)

**43.** [P2] Support for SLA (resin) printers: the arrange settings panel disables Snuggle for SLA profiles since 3D voxel nesting is irrelevant for platforms where parts are oriented arbitrarily in XYZ. The `m_arrange_settings_sla` struct should not carry Snuggle fields.

---

## 7. Predictability

**44.** [P0] Running arrange twice with identical inputs (same parts, same settings, same seed) produces the same result. The genetic algorithm must use a deterministic seed derived from input hash unless the user explicitly requests randomization.

**45.** [P1] Adding one small part to a previously arranged bed does not rearrange all existing parts. Ideal behavior: existing parts stay in place, new part is placed in the best available gap. *This requires an incremental mode, which the current full-population GA does not support.*

**46.** [P1] The quality slider produces monotonically better results: quality=7 should never produce a worse packing than quality=5 for the same input. (Statistically, higher generations should converge to equal or better fitness. Verify by ensuring elitism preserves the best individual.)

**47.** [P2] Arrange result is visually "clean" -- parts align to a grid or cluster neatly rather than scattering randomly across the bed. The `w_clustering` fitness weight (currently 1.5) should produce perceptually organized layouts, not technically-optimal-but-messy ones.

**48.** [P2] Parts of similar height are placed near each other when possible, reducing Z-hop travel. The `w_height_center` weight (currently 0.3) contributes to this, but it should be tunable.

---

## 8. Accessibility

**49.** [P1] Every Snuggle control has a tooltip. Current state: "Quality" slider has a tooltip ("Higher = better packing, slower") but "Part gap" and "Lock rotation" do not. Add tooltips:
  - **Part gap**: "Minimum clearance between parts in millimeters. Increase for easier part removal."
  - **Lock rotation**: "Keep each part's current Z-rotation. Disable to let the nester rotate parts for tighter packing."
  - **Use Snuggle 3D nester**: "Experimental 3D-aware arranger using genetic optimization. Considers part height for better nesting of tall and short parts."

**50.** [P1] The "Snuggle (experimental)" label should include a visual indicator (icon or color -- currently teal text, which is good) distinguishing it from the stable arrange controls above. If Snuggle becomes non-experimental, the label and color should update to match the standard UI style.

---

## Gap Analysis vs. Current Implementation

| Area | Status | Blocking Issues |
|------|--------|-----------------|
| Placement correctness | Partial | Exclusion zones ignored (#3, #38) |
| Visual feedback | Minimal | No progress callback wired to UI (#8, #11) |
| Controls | Basic | No cancel (#18), no undo (#19), no advanced panel (#22) |
| Performance | Functional | GPU fallback untested (#35), priority not set (#31) |
| Error handling | Weak | Timeout behavior unclear (#34), no per-part fallback (#36) |
| Integration | Partial | Wipe tower ignored (#38), no multi-plate (#40), no seq-print clearance (#41) |
| Predictability | Unknown | Deterministic seeding not verified (#44), no incremental mode (#45) |
| Accessibility | Partial | Missing tooltips (#49) |

Top 5 items to fix before any public release:
1. Exclusion zone handling (#3, #38) -- wipe tower collision will destroy prints
2. Progress/cancel UI (#8, #18) -- 30-second hang with no feedback is unacceptable
3. Undo support (#19) -- users expect Ctrl+Z to work after any operation
4. Timeout error reporting (#34) -- silent failure when parts don't fit causes waste
5. Deterministic seeding (#44) -- non-reproducible results erode trust
