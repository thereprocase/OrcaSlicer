# WAR COUNCIL LEDGER
## Snuggle Nester — Reconciled Findings
**Date:** 2026-03-31
**Reviewers:** Sauron, Gandalf, Frodo, Aragorn, Legolas, Gimli (absent), Gollum, Uruk-Hai ×4 (Voxelizer, Nester, GPU, Bridge+UI)
**Source files:** `SnuggleArrange.cpp`, `snuggle_nester.hpp`, `polite_voxelizer.hpp`, `gpu_collision.cpp/.hpp`, `GLCanvas3D.cpp`

---

## PART I — CONFIRMED BUGS

### CRITICAL

---

**W-001**
**Severity:** CRITICAL
**File:line:** `SnuggleArrange.cpp:196,205`
**Description:** Dangling else-brace causes voxelization-failure path to push an empty `PartInfo` twice — once inside the error branch and once unconditionally at line 205 — corrupting `parts.size()` vs `items.size()` and cascading all placement writeback into off-by-one index errors.
**Found by:** Uruk-Hai (Bridge+UI)
**Fix:** Replace the dangling brace with a proper `if/else` block and a single unconditional `parts.push_back(std::move(pi))` after the branch.

---

**W-002**
**Severity:** CRITICAL
**File:line:** `polite_voxelizer.hpp:552–558`
**Description:** Mesh index values `i0/i1/i2` from the 3MF `indices` array are used to index `vertices[]` with no bounds check against `num_verts`; a crafted mesh with `indices[k] = 0xFFFFFFFF` produces an out-of-bounds heap read.
**Found by:** Aragorn
**Fix:** Add `if (i0 >= num_verts || i1 >= num_verts || i2 >= num_verts) continue;` before the vertex dereference.

---

**W-003**
**Severity:** CRITICAL
**File:line:** `gpu_collision.cpp:249–332` (init_context)
**Description:** On failure paths after `GetDC()` succeeds — specifically when `wglCreateContext`, `glewInit`, or the OpenGL version check fail — `ReleaseDC()` is never called, leaking the device context.
**Found by:** Uruk-Hai (GPU)
**Fix:** Track `hwnd`/`hdc` as locals and call `ReleaseDC(hwnd, hdc)` (and `DestroyWindow`) on every error return path, or wrap in RAII.

---

### HIGH

---

**W-004**
**Severity:** HIGH
**File:line:** `gpu_collision.cpp:477` (evaluate_batch) + shader lines ~174–177
**Description:** `bed_margin` is silently discarded with `(void)bed_margin`; the GPU shader checks out-of-bounds against raw `u_bed_w/u_bed_h` with a zero margin, so the GPU and CPU evaluators produce different feasibility verdicts for identical populations — GPU accepts placements the CPU would reject.
**Found by:** Gandalf, Aragorn (M5), Legolas (issue #8), Gollum
**Fix:** Pass `bed_margin` as a `uniform float u_bed_margin` to the shader and apply it in the OOB condition (`pmin_x < u_bed_margin`, etc.).

---

**W-005**
**Severity:** HIGH
**File:line:** `gpu_collision.cpp:96–101,417` (GridMeta::data_offset)
**Description:** `data_offset` in `GridMeta` is `uint32_t`, but the accumulation variable `offset` is `size_t`; with sufficient parts the cast silently wraps at 4 GB, making the shader read voxel data from the wrong offset and reporting zero collisions for genuinely overlapping parts.
**Found by:** Aragorn (H2)
**Fix:** Change `GridMeta::data_offset` to `uint64_t` and update the shader layout binding to match, or add an explicit overflow assertion before the cast.

---

**W-006**
**Severity:** HIGH
**File:line:** `polite_voxelizer.hpp:161–165` (world_to_grid)
**Description:** Division by `voxel_size` with no guard; if a `VoxelGrid` with `voxel_size == 0` reaches this function the result is undefined behavior (NaN/Inf cast to int at three call sites: `collision_count` lines 313/319, `voxelize_mesh` lines 471/475/479).
**Found by:** Aragorn (M1 — identifies the loop consequence), Uruk-Hai (Voxelizer, Bug #1/#2)
**Fix:** Add `if (voxel_size <= 0.0f) return;` (or return sentinel) at the top of `world_to_grid`, and add `if (vs <= 0.0f) return 0;` before the triple loop in `collision_count`.

---

**W-007**
**Severity:** HIGH
**File:line:** `SnuggleArrange.cpp:122–128`
**Description:** Part-to-instance matching is done by object name; two instances of the same object (or imports with identical names) cause the second instance to never match, silently dropping it from the arrangement — the same fragility also creates O(N²) string comparisons exploitable for DoS with many same-named objects.
**Found by:** Gandalf, Aragorn (M4), Uruk-Hai (Bridge+UI, Bug #2)
**Fix:** Build a `std::multimap<std::string, ObjInst*>` once before the loop for O(log N) lookups; or, preferred long-term, store a stable instance ID in `ArrangePolygon` at creation time and match by ID.

---

**W-008**
**Severity:** HIGH
**File:line:** `gpu_collision.cpp:37–41` (CpuCollisionEvaluator::get_rotated)
**Description:** `rot_cache_` is a raw pointer set by `upload_grids`; `get_rotated` has no null guard and no bounds check on `part_idx` against `rot_cache_->size()` — a call before `upload_grids` or with a mismatched part count crashes with a null or OOB dereference.
**Found by:** Aragorn (H1)
**Fix:** Assert or return a default grid if `rot_cache_ == nullptr || part_idx >= rot_cache_->size()`.

---

**W-009**
**Severity:** HIGH
**File:line:** `GLCanvas3D.cpp:5945–6092` (disabled_begin/end scope)
**Description:** The `disabled_begin(true)` / `disabled_end()` pair is guarded by the *current* value of `use_snuggle`; if the user enables Snuggle while the Advanced TreeNode is open, `disabled_end()` is never called on the subsequent frame, permanently corrupting ImGui's disabled-stack for all downstream controls.
**Found by:** Uruk-Hai (Bridge+UI, Bug #5)
**Fix:** Capture `bool snuggle_disabled = !settings_out.use_snuggle` before the block and use that captured value for both `disabled_begin` and `disabled_end`.

---

### MEDIUM

---

**W-010**
**Severity:** MEDIUM
**File:line:** `snuggle_nester.hpp:567–570` (snap_rotation)
**Description:** Rotation normalization uses `while (delta < 0) delta += TWO_PI; while (delta >= TWO_PI) delta -= TWO_PI;`; a crafted extreme float (e.g., `1e18f`) from a 3MF's rotation field causes ~2×10¹⁷ loop iterations, stalling the UI thread.
**Found by:** Aragorn (M2), Legolas (issue #16)
**Fix:** Replace both while-loops with `delta = fmodf(delta, TWO_PI); if (delta < 0) delta += TWO_PI;`.

---

**W-011**
**Severity:** MEDIUM
**File:line:** `snuggle_nester.hpp:889,894` (compute_fitness)
**Description:** `bed_diag` is computed as `sqrt(width² + height²)`; if both bed dimensions are zero, `bed_diag = 0` and the subsequent `avg_dist /= (n * bed_diag)` is a division by zero producing NaN that poisons all fitness values.
**Found by:** Uruk-Hai (Nester, Bug #3)
**Fix:** Guard: `if (bed_diag > 0.0001f) { avg_dist /= (n * bed_diag); ... } else { clustering = 0.0f; }`.

---

**W-012**
**Severity:** MEDIUM
**File:line:** `SnuggleArrange.cpp:323–341` (greedy accept)
**Description:** `parts[j].grid.rotated_copy()` is called on previously-placed parts after all grids were moved into the `parts` vector via `push_back`; if `VoxelGrid` holds heap-allocated data and its move constructor leaves the source in a valid-but-empty state, calling `rotated_copy` on a moved-from grid is undefined behavior.
**Found by:** Uruk-Hai (Bridge+UI, Bug #3)
**Fix:** Confirm `VoxelGrid` move leaves the source in a safely-callable state (document this contract), or store grids by value in a non-moving container and access by index.

---

**W-013**
**Severity:** MEDIUM
**File:line:** `GLCanvas3D.cpp:1177–1179` (config load)
**Description:** On config load, `snuggle_rotation_step` is snapped to the nearest valid dropdown value, but `snuggle_lock_rotation` is not re-derived from the snapped value — leaving the two fields out of sync until the user touches the dropdown, which causes the GA to run with the wrong rotation mode.
**Found by:** Uruk-Hai (Bridge+UI, Bug #6)
**Fix:** After the snap on load, add `m_arrange_settings_fff.snuggle_lock_rotation = (m_arrange_settings_fff.snuggle_rotation_step == 0);`.

---

**W-014**
**Severity:** MEDIUM
**File:line:** `snuggle_nester.hpp:376–490` (compact_toward_center)
**Description:** Inside the micro-rotation sub-loop, `rotated_copy(test_rot)` for the current part and `rotated_copy()` for every other part are called inside the `bs2` binary-search loop (8 iterations per angle per part per sweep); the angle does not change inside `bs2`, so all 8 copies are identical — this is the single most expensive unnecessary allocation in the codebase.
**Found by:** Legolas (issues #1/#5)
**Fix:** Hoist `rotated_copy(test_rot)` and `rot_others` construction to before the `bs2` loop; the existing `rot_others_storage` vector is the right mechanism but is rebuilt instead of reused.

---

**W-015**
**Severity:** MEDIUM
**File:line:** `snuggle_nester.hpp:336` (compact_toward_center)
**Description:** `parts[0].grid.voxel_size` is accessed unconditionally; if `parts` is empty while `result.placements` is somehow non-empty, this is an out-of-bounds access.
**Found by:** Uruk-Hai (Nester, Bug #2)
**Fix:** Add `if (parts.empty()) return;` before `parts[0]` is accessed.

---

**W-016**
**Severity:** MEDIUM
**File:line:** `polite_voxelizer.hpp:215–222` (rotated_copy fallback)
**Description:** When allocation for the rotated grid fails, `rotated_copy` silently returns an unrotated copy; the caller receives wrong collision geometry with no error indication, potentially accepting placements that actually overlap.
**Found by:** Sauron, Gollum
**Fix:** Add `BOOST_LOG_TRIVIAL(warning)` on the fallback path and document the silent-degradation contract in the header comment.

---

**W-017**
**Severity:** MEDIUM
**File:line:** `snuggle_nester.hpp:702–719` (seed_bottom_left)
**Description:** Triple-nested loop (bed_area/gap² × n_parts) with no yield and no cap; with small `min_gap_mm` and many parts, this runs O(n³) synchronously before the GA timeout begins counting, stalling the UI thread for seconds.
**Found by:** Aragorn (M3), Legolas (issue #3)
**Fix:** Add a `polite_yield()` every N outer-loop iterations and enforce a `max_seed_tries` cap.

---

### LOW

---

**W-018**
**Severity:** LOW
**File:line:** `gpu_collision.cpp:254–262` (init_context — window class name)
**Description:** The window class name `"SnuggleGPUCollision"` is a static string; if a different module registers the same name first with its own `WndProc`, `ERROR_CLASS_ALREADY_EXISTS` is swallowed and `CreateWindowA` silently uses the foreign proc — a low-exploitation window-squatting vector in shared-process or plugin scenarios.
**Found by:** Aragorn (H3)
**Fix:** Append the process ID to the class name: `"SnuggleGPUCollision_" + std::to_string(GetCurrentProcessId())`.

---

**W-019**
**Severity:** LOW
**File:line:** `gpu_collision.cpp:549–580` (cleanup)
**Description:** `cleanup()` calls `DestroyWindow` and `wglDeleteContext` but never calls `UnregisterClassA`; in a process that creates and destroys `GpuCollisionEvaluator` multiple times, the class registration leaks and the `ERROR_CLASS_ALREADY_EXISTS` guard in `init_context` prevents re-registration forever after.
**Found by:** Aragorn (L1)
**Fix:** Add `UnregisterClassA("SnuggleGPUCollision...", GetModuleHandle(nullptr))` in `cleanup()`.

---

**W-020**
**Severity:** LOW
**File:line:** `gpu_collision.cpp:549–562` (cleanup — null guard)
**Description:** `cleanup()` accesses `gl_dc_` unconditionally when `gl_context_` is set; `gl_dc_` could theoretically be null if `init_context` failed before setting it, making `wglMakeCurrent(nullptr, ctx)` have platform-defined behavior.
**Found by:** Uruk-Hai (GPU, Bug #4)
**Fix:** Change the guard to `if (gl_context_ && gl_dc_)`.

---

**W-021**
**Severity:** LOW
**File:line:** `gpu_collision.cpp:319` (GL_RENDERER logging)
**Description:** `glGetString(GL_RENDERER)` is an untrusted driver-supplied string logged directly into the Boost log stream; a crafted or spoofed driver could inject log-newline or ANSI escape characters.
**Found by:** Aragorn (H4)
**Fix:** Sanitize or wrap the string: strip non-printable characters before logging.

---

**W-022**
**Severity:** LOW
**File:line:** `snuggle_nester.hpp:526` (build_rotation_cache)
**Description:** `rot_cache_built_` is set but never reset; reusing a `SnuggleNester` instance with different parts or config silently uses the stale cache.
**Found by:** Aragorn (L2)
**Fix:** Document that `SnuggleNester` is single-use (current practice), or add cache invalidation tied to the `parts` vector identity.

---

**W-023**
**Severity:** LOW
**File:line:** `SnuggleArrange.cpp:73–82`
**Description:** Parts that fail the oversized check get `bed_idx = -1` but still proceed through full voxelization and rotation-cache building before being discarded — wasting CPU time proportional to part size.
**Found by:** Sauron
**Fix:** `continue` after marking `bed_idx = -1` to skip voxelization entirely for oversized parts.

---

**W-024**
**Severity:** LOW
**File:line:** `gpu_collision.cpp:349–352` (compile_shader)
**Description:** Shader info log uses a fixed 2048-byte stack buffer; verbose Mesa drivers can produce error messages exceeding 2 KB, truncating silently.
**Found by:** Aragorn (L3)
**Fix:** Use `glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len)` to query actual length before allocating.

---

**W-025**
**Severity:** LOW
**File:line:** `GLCanvas3D.cpp:5956,6011–6013,6028,6042,6056,6067` (sliders/inputs)
**Description:** Several `BBLDragFloat`/`InputInt` widgets treat their min/max as suggested limits, not hard limits; a user can drag outside the range before the post-widget `std::clamp` fires, and if the value is committed before sync the out-of-range value can persist in the saved config.
**Found by:** Uruk-Hai (Bridge+UI, Bug #4)
**Fix:** Apply the clamp before reading from the widget, or enforce hard limits via widget flags.

---

**W-026**
**Severity:** LOW
**File:line:** `polite_voxelizer.hpp:126–132` (VoxelGrid::allocate)
**Description:** `nx_ * ny_ * nz_` is an unchecked multiplication that relies implicitly on the `MAX_GRID_DIM = 512` per-dimension cap to stay within `size_t`; if `MAX_GRID_DIM` is ever raised without reanalysis, the product can overflow and bypass the size guard.
**Found by:** Aragorn (C2), Uruk-Hai (Voxelizer, Bug #3)
**Fix:** Add an explicit pre-multiplication overflow check (e.g., `if (nx_ > MAX_TOTAL_VOXELS / ny_ / nz_) return GRID_TOO_LARGE;`) so the guard is self-contained.

---

---

## PART II — ARCHITECTURAL / STYLE CONCERNS

### HIGH (architecture — fix before upstream PR)

---

**W-027**
**Severity:** HIGH (architecture)
**File:line:** `gpu_collision.hpp` (`ROT_BINS`) + `snuggle_nester.hpp` (`ROT_CACHE_BINS`) + `gpu_collision.cpp` (`ROT_CACHE_BINS`)
**Description:** The rotation-bin count of 360 is defined three times under two different names (`ROT_CACHE_BINS` in CPU evaluator and nester; `ROT_BINS` in GPU evaluator); changing one without the others silently breaks the angle-to-bin mapping.
**Found by:** Gandalf, Gollum
**Fix:** Rename `GpuCollisionEvaluator::ROT_BINS` to `ROT_CACHE_BINS` and promote both to a single `constexpr` in a shared header (e.g., `snuggle_types.hpp` or the `CollisionEvaluator` base class).

---

**W-028**
**Severity:** HIGH (architecture)
**File:line:** `snuggle_nester.hpp` (multiple lines: 544, 567, 576, 603, 763, 803) + `SnuggleArrange.cpp:232,371` + `gpu_collision.cpp:38`
**Description:** The value `3.14159265f` (π) appears as a bare literal in nine separate locations across three files with no named constant, and the shader uses `6.2831853` (2π, different precision) for the same purpose.
**Found by:** Gandalf, Gollum
**Fix:** Define `constexpr float PI_F = 3.14159265358979f;` and `constexpr float TWO_PI_F = 2.0f * PI_F;` once in a shared snuggle header; use `TWO_PI_F` in the shader source string at full precision.

---

**W-029**
**Severity:** MEDIUM (architecture)
**File:line:** `snuggle_nester.hpp` (entire file, ~1000 LOC)
**Description:** The nester is entirely header-only; while harmless with one consumer today, adding a second consumer (tests, benchmark, alternate bridge) multiplies compile times and leaks implementation details into the interface.
**Found by:** Gandalf
**Fix:** Split into `snuggle_nester.hpp` (class declaration + public types only) and `snuggle_nester.cpp` (all method bodies) before upstream PR submission.

---

**W-030**
**Severity:** MEDIUM (architecture)
**File:line:** `polite_voxelizer.hpp` (type definitions section)
**Description:** `VoxelGrid`, `Vec3f`, and `VoxError` — core types used by both the voxelizer and the nester — are defined in the voxelizer header, creating an unnecessary dependency between the nester and the voxelizer implementation detail.
**Found by:** Gandalf
**Fix:** Extract these types into a minimal `snuggle_types.hpp` that both headers include.

---

**W-031**
**Severity:** MEDIUM (architecture)
**File:line:** `gpu_collision.cpp:55–57` (evaluate_batch CPU — inner loop)
**Description:** A fresh `std::vector<const VoxelGrid*>` of size n is heap-allocated for every individual in every generation (~51,200 allocations over a full run); the vector is the same size for every individual.
**Found by:** Legolas (issue #4)
**Fix:** Declare the vector once outside the individual loop and reuse it with `.clear()`.

---

### MEDIUM (architecture)

---

**W-032**
**Severity:** MEDIUM (architecture)
**File:line:** `gpu_collision.cpp:496–529` (GpuCollisionEvaluator::evaluate_batch)
**Description:** Every generation does: upload placement SSBO (orphaning the old buffer), dispatch, `glMemoryBarrier`, then `glGetBufferSubData` — a synchronous CPU-GPU sync that sits the CPU idle for ~0.5–2 ms per generation; additionally, `glGetUniformLocation` is queried every generation instead of being cached at link time.
**Found by:** Legolas (issue #7)
**Fix:** Cache all uniform locations after `glLinkProgram`; double-buffer the placement SSBO with fence sync to eliminate the re-orphan stall; replace `glGetBufferSubData` with `glMapBufferRange` + `glClientWaitSync` for driver scheduling freedom.

---

**W-033**
**Severity:** MEDIUM (architecture)
**File:line:** `snuggle_nester.hpp:920–932` (compute_fitness proximity loop)
**Description:** Part centroid positions are recomputed three times within a single `compute_fitness` call (centroid, avg_dist, proximity loop); with a feasible population of 512 individuals this is ~97 K redundant sqrt() calls per generation.
**Found by:** Legolas (issue #2)
**Fix:** Compute part center coordinates once into a local array and reuse across all three passes.

---

**W-034**
**Severity:** MEDIUM (architecture)
**File:line:** `gpu_collision.cpp` shader, line ~116
**Description:** GPU shader uses `local_size_x = 1` — one invocation evaluates all pairs serially, wasting 31/32 of each GPU warp's capacity; for upstream PR this is flagged as a known limitation in a comment but will require redesign for acceptable GPU performance.
**Found by:** Gandalf, Legolas (issue #12)
**Fix:** Redesign the work-group to assign pairs to lanes and use a parallel reduction for the final collision count — defer to a dedicated GPU pass, but document current limitation clearly.

---

**W-035**
**Severity:** MEDIUM (architecture)
**File:line:** `polite_voxelizer.hpp:308–329` (collision_count triple loop)
**Description:** The overlap region is iterated in floating-point world-space coordinates with `world_to_grid()` (three `floor` + three subtracts) called on every voxel; the loop also calls `polite_yield()` every 10K iterations, adding unnecessary syscall overhead inside the innermost collision loop.
**Found by:** Legolas (issue #9)
**Fix:** Convert loop bounds to integer grid coordinates once before the loops; remove `polite_yield` from inside `collision_count` (the outer generation loop already provides responsiveness).

---

### LOW / STYLE

---

**W-036**
**Severity:** STYLE
**File:line:** `SnuggleArrange.cpp:197`
**Description:** Variable `hull_area_mm2` is actually the grid bounding-rectangle area (`nx * ny * voxel_size²`), not the convex hull area; a future maintainer could build incorrect logic on this assumption.
**Found by:** Sauron
**Fix:** Rename to `grid_bbox_area_mm2` or add a comment explaining the approximation.

---

**W-037**
**Severity:** STYLE
**File:line:** `polite_voxelizer.hpp:268–274` (count_solid)
**Description:** Manual bit-serial popcount (`while (b) { count += b & 1; b >>= 1; }`) is replaceable with `__builtin_popcount` / `std::popcount` (C++20) for a single-instruction implementation.
**Found by:** Gandalf, Legolas (issue #14)
**Fix:** Replace with `std::popcount` (C++20) or `__builtin_popcount` guarded by a `__has_builtin` check.

---

**W-038**
**Severity:** STYLE
**File:line:** `snuggle_nester.hpp:708`
**Description:** `try_x = pr.x + pr.w + cfg_.min_gap_mm - cfg_.min_gap_mm` simplifies to `pr.x + pr.w` — the subtraction cancels the addition, likely a copy-paste artifact; the correct skip target to maintain the gap should be `pr.x + pr.w + min_gap_mm`.
**Found by:** Sauron
**Fix:** Remove the `- cfg_.min_gap_mm` term.

---

**W-039**
**Severity:** STYLE
**File:line:** `GLCanvas3D.cpp` (Snuggle panel, lines ~5950–6074)
**Description:** Four controls have no tooltip: Part gap slider (units unknown), Timeout input (behavior-on-fire unknown), Max parts input (silent fallback not explained), and Multi-plate overflow checkbox (behavior ambiguous); two controls (voxel resolution, part gap) use format strings without "mm" unit suffix.
**Found by:** Frodo (P1, P5, P6, P7, P8)
**Fix:** Add `ImGui::SetTooltip(...)` blocks for all four controls; add " mm" to format strings for gap and resolution.

---

**W-040**
**Severity:** STYLE
**File:line:** `GLCanvas3D.cpp:5964` and `GLCanvas3D.cpp:6080`
**Description:** "Compact after arrange" checkbox appears twice (main section and Advanced section) controlling the same setting with conflicting tooltip guidance; the advanced tooltip discourages use while the main one presents it neutrally.
**Found by:** Frodo (P2)
**Fix:** Remove the main-section duplicate; keep only the Advanced-section instance with the more informative tooltip.

---

**W-041**
**Severity:** STYLE
**File:line:** `SnuggleArrange.cpp` (multiple), `snuggle_nester.hpp` (multiple)
**Description:** Non-ASCII Unicode em-dash separators (`──`) appear in 13+ comment dividers; these may cause issues in diff tools, syntax highlighters, or build logs with ASCII-only encoding.
**Found by:** Gollum
**Fix:** Replace with ASCII separators (`//----` or `//====`).

---

**W-042**
**Severity:** STYLE
**File:line:** `GLCanvas3D.cpp:~5922`
**Description:** The panel section header "Snuggle 3D Arrangement" carries no "(experimental)" or "(beta)" label; new users have no signal this feature is non-production and no brief subtitle explaining what it does.
**Found by:** Frodo (P3)
**Fix:** Change header to "Snuggle (experimental) — 3D-aware part nesting" or equivalent.

---

**W-043**
**Severity:** STYLE
**File:line:** `SnuggleArrange.cpp` (progress callback, status bar messages)
**Description:** Status bar progress messages show developer-facing detail: "Snuggle CPU gen 12/30 3 collisions" — "gen", "CPU", and collision count are implementation noise to end users; "3 collisions" in particular reads as alarming.
**Found by:** Frodo (P4)
**Fix:** Replace with user-facing copy: "Arranging parts... (NN%)"; move generation/collision detail to the debug log only.

---

**W-044**
**Severity:** STYLE
**File:line:** `SnuggleArrange.cpp:~57` (max-parts guard), `SnuggleArrange.cpp:~188–191` (voxelization failure)
**Description:** Both the max-parts limit and per-part voxelization failure are silent to the user — only a `BOOST_LOG_TRIVIAL(warning)` is emitted; users see either a normal-looking arrangement (max-parts hit) or a missing part (voxelization fail) with no explanation.
**Found by:** Frodo (P6, P9)
**Fix:** Surface both as visible status-bar or notification messages matching OrcaSlicer's existing warning pattern.

---

---

## SUMMARY COUNTS

| Severity | Count |
|---|---|
| CRITICAL (confirmed bugs) | 3 |
| HIGH (confirmed bugs) | 6 |
| MEDIUM (confirmed bugs) | 8 |
| LOW (confirmed bugs) | 9 |
| HIGH (architecture) | 2 |
| MEDIUM (architecture) | 6 |
| STYLE | 9 |
| **Total entries** | **43** |

---

## FIX BEFORE NEXT BUILD

These items either produce wrong output, corrupt state, or crash. Fix all of them before the next soak-test cycle.

| ID | Why it cannot wait |
|---|---|
| W-001 | Double-push corrupts parts[] — every placement writeback is wrong |
| W-002 | OOB heap read from any 3MF with a bad index — reachable today |
| W-003 | DC leak on GPU init failure — resource exhaustion in long sessions |
| W-004 | GPU and CPU report different feasibility — GA result is unreliable on GPU path |
| W-005 | uint32 overflow silently passes infeasible arrangements as feasible |
| W-006 | Division by zero in world_to_grid — crash or UB at multiple call sites |
| W-007 | Silent part drop on duplicate names — user loses parts with no warning |
| W-008 | Null deref if get_rotated called before upload_grids |
| W-009 | ImGui disabled-stack corruption — entire panel freezes |
| W-010 | While-loop DoS from crafted rotation value in 3MF |
| W-011 | NaN fitness poisons entire population if bed dims are zero |
| W-013 | lock_rotation/rotation_step desync on config load — wrong GA behavior on startup |
| W-027 | ROT_CACHE_BINS triplicated — silent mismatch breaks all rotation lookups if any copy changes |
| W-028 | Pi literals scattered — correctness and maintainability; fix alongside W-027 |

---

## DEFER TO POLISH PASS

Safe to merge without fixing; address before upstream PR submission.

**Architecture:** W-012, W-014, W-015, W-016, W-017, W-026, W-029, W-030, W-031, W-032, W-033, W-034, W-035

**Low/safety:** W-018, W-019, W-020, W-021, W-022, W-023, W-024, W-025

**Style/UX:** W-036, W-037, W-038, W-039, W-040, W-041, W-042, W-043, W-044

---

## ADDENDUM — Late Arrivals (Legolas, Gimli)

**W-045**
**Severity:** HIGH
**File:line:** `snuggle_nester.hpp:340-469` (compact_toward_center)
**Description:** Micro-rotation section rebuilds `rotated_copy()` for every other part, every candidate angle, every binary search step, every part, every sweep. Hundreds of redundant 3D grid copies. Worst hot path in the codebase.
**Found by:** Legolas
**Fix:** Pre-build rotated grids for all parts at sweep start, reuse across binary search and micro-rotation attempts.

**W-046**
**Severity:** MEDIUM
**File:line:** `gpu_collision.cpp:507-530` (evaluate_batch)
**Description:** `glGetBufferSubData` is a hard CPU-GPU sync — CPU blocks until GPU drains. Combined with SSBO re-orphan via `glBufferData` each generation, creates unnecessary stalls.
**Found by:** Legolas
**Fix:** Cache uniform locations at link time, double-buffer placement SSBO, use fence + `glMapBufferRange` instead of `glGetBufferSubData`.

**W-047**
**Severity:** MEDIUM
**File:line:** `gpu_collision.cpp:55-59` (evaluate_batch inner loop)
**Description:** `vector<const VoxelGrid*>` allocated inside per-individual loop — 512 heap allocs per generation. Should be one allocation outside the loop.
**Found by:** Legolas
**Fix:** Move vector declaration outside the individual loop, `.resize()` once.

**W-048**
**Severity:** LOW
**File:line:** `gpu_collision.cpp:260-310` (init_context)
**Description:** Uses legacy `wglCreateContext()` instead of `wglCreateContextAttribsARB` with explicit GL 4.3 core profile request. Works on most drivers but relies on driver goodwill.
**Found by:** Gimli
**Fix:** Use `wglCreateContextAttribsARB` with `WGL_CONTEXT_MAJOR_VERSION_ARB=4, WGL_CONTEXT_MINOR_VERSION_ARB=3`.

**W-049**
**Severity:** LOW
**File:line:** `gpu_collision.hpp:19-21`
**Description:** Forward-declares `PartInfo` and `Individual` but never includes `snuggle_nester.hpp`. Any standalone includer (test file) gets unresolved types.
**Found by:** Gimli
**Fix:** Document the include-order contract in a comment, or move the struct definitions to a shared types header.

---

## UPDATED COUNTS

| Severity | Bugs | Architecture | Style | Total |
|----------|------|--------------|-------|-------|
| CRITICAL | 3 | — | — | 3 |
| HIGH | 7 | 2 | — | 9 |
| MEDIUM | 10 | 6 | — | 16 |
| LOW | 11 | — | 9 | 20 |
| **Total** | **31** | **8** | **9** | **49** |
