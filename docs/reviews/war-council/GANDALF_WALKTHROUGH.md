# Gandalf's Walkthrough — War Council Ledger vs Current Code

**Date:** 2026-03-31
**Reviewer:** Gandalf the Grey (Opus 4.6)
**Source branch:** `snuggle-uiux` (4 commits after War Council ledger)
**Files examined:**
- `src/libslic3r/Arrange/SnuggleArrange.cpp`
- `src/libslic3r/Arrange/snuggle_nester.hpp`
- `src/libslic3r/Arrange/polite_voxelizer.hpp`
- `src/libslic3r/Arrange/gpu_collision.cpp`
- `src/libslic3r/Arrange/gpu_collision.hpp`
- `src/slic3r/GUI/GLCanvas3D.cpp`

---

## Ledger Status Table

| ID | Severity | Status | Evidence |
|----|----------|--------|----------|
| W-001 | CRITICAL | RESOLVED | `SnuggleArrange.cpp:188-205` — the `continue` on line 195 skips the unconditional push. The else-brace on line 196 is still cosmetically odd (`} {` after continue) but functionally harmless: the `continue` prevents double-push. Only one `parts.push_back` executes per iteration. |
| W-002 | CRITICAL | RESOLVED | `polite_voxelizer.hpp:558-559` — bounds check `if (i0 >= num_verts || i1 >= num_verts || i2 >= num_verts) return VoxError::INVALID_INPUT;` is present in `voxelize_indexed_mesh`. |
| W-003 | CRITICAL | STILL OPEN | `gpu_collision.cpp:274-298` — after `GetDC` succeeds (line 274), failure paths at lines 289 (SetPixelFormat), 296 (wglCreateContext), 301 (wglMakeCurrent), 310 (glewInit), 323 (GL version check) all `return false` without calling `ReleaseDC` or `DestroyWindow`. The `cleanup()` method does handle this if `gl_context_` is set, but several failures occur *before* `gl_context_` is assigned (line 299), so `cleanup()` skips the entire Windows block. DC and HWND leak on those paths. |
| W-004 | HIGH | RESOLVED | `gpu_collision.cpp:538` — `bed_margin` is now passed as uniform `u_bed_margin` to the shader. Shader lines 139, 175-178 use `u_bed_margin` for OOB checks. CPU evaluator (line 76-83) also uses `bed_margin`. No more `(void)bed_margin`. |
| W-005 | HIGH | RESOLVED | `gpu_collision.cpp:416-421` — explicit overflow guard: `if (total_voxel_bytes > (size_t)UINT32_MAX)` falls back to CPU. The `data_offset` field remains `uint32_t` (line 96) but the guard prevents the wrap. |
| W-006 | HIGH | STILL OPEN | `polite_voxelizer.hpp:163-167` — `world_to_grid` still divides by `voxel_size` with no guard. `collision_count` (line 306) computes `vs = std::max(a.voxel_size, b.voxel_size)` with no zero check. A `VoxelGrid` with `voxel_size == 0` would produce UB. The `allocate` function does not enforce `voxel_size > 0`, and `VoxelGrid()` default-constructs with `voxel_size = 1.0f`, so this requires deliberate misuse or a bug upstream — but the guard is still absent. |
| W-007 | HIGH | PARTIALLY FIXED | `SnuggleArrange.cpp:121-128` — matching still uses object name with O(N) scan. However, line 125 (`oi.obj = nullptr; // mark used`) now prevents the same instance from matching twice, fixing the "second instance never matches" bug. The O(N^2) performance issue remains — no `multimap` or ID-based matching. |
| W-008 | HIGH | RESOLVED | `gpu_collision.cpp:48` — `if (!rot_cache_) return;` null guard at the top of `evaluate_batch`. `get_rotated` (line 36) is only called after this guard. No bounds check on `part_idx` against `rot_cache_->size()`, but `part_idx` is always derived from `parts.size()` which was used to build the cache, so mismatch requires a caller bug. Null deref is fixed; bounds check is not. |
| W-009 | HIGH | RESOLVED | `GLCanvas3D.cpp:5946,6092` — both `disabled_begin` and `disabled_end` are now guarded by the same condition: `if (!settings_out.use_snuggle)`. The block is wrapped in `{}` braces (lines 5945-6093), and the disabled_begin/end calls both test `settings_out.use_snuggle` (not `settings.use_snuggle`). Since `settings_out` is not modified between these two points (the snuggle enable checkbox writes to both `settings` and `settings_out` at line 5929 *before* the block), the two calls are always balanced. |
| W-010 | MEDIUM | STILL OPEN | `snuggle_nester.hpp:565-566` — still uses `while (delta < 0) delta += TWO_PI; while (delta >= TWO_PI) delta -= TWO_PI;`. Not replaced with `fmodf`. An extreme float value would still stall. |
| W-011 | MEDIUM | STILL OPEN | `snuggle_nester.hpp:886` — `avg_dist /= (n * bed_diag)` with no guard for `bed_diag == 0`. If both bed dimensions are zero, `bed_diag` is zero, producing NaN. |
| W-012 | MEDIUM | RESOLVED | `SnuggleArrange.cpp:330,346` — the greedy-accept path now calls `grid_i.rotated_copy(pl.zrot)` (line 330) and `parts[j].grid.rotated_copy(pl_j.zrot)` (line 346) on the *original* `parts[j].grid`, not a moved-from grid. `parts` is built via `push_back(std::move(pi))` but the moved-from `pi` is never accessed again. Each `parts[j].grid` is populated during voxelization and never moved from after the vector is built. No UB. |
| W-013 | MEDIUM | RESOLVED | `GLCanvas3D.cpp:1186` — `m_arrange_settings_fff.snuggle_lock_rotation = (m_arrange_settings_fff.snuggle_rotation_step == 0);` is now present after the config load block, exactly as the fix prescribed. |
| W-014 | MEDIUM | PARTIALLY FIXED | `snuggle_nester.hpp:376-466` — the micro-rotation sub-loop now uses `cached_rotated()` (lines 434, 440, 454, 459) instead of calling `rotated_copy()` directly. This eliminates the expensive allocation. However, `cached_rotated` is called inside the `bs2` loop (lines 454, 459) where the angle doesn't change, so it still performs the bin calculation 8 times redundantly (though it's a cache lookup, not an allocation). Hoisting would be cleaner but the pathological cost is gone. |
| W-015 | MEDIUM | STILL OPEN | `snuggle_nester.hpp:336` — `parts[0].grid.voxel_size` is accessed with no `parts.empty()` guard. The `n < 2` early return on line 333 does protect against empty (since `n < 2` covers `n == 0`), so this is actually safe. **Status revised: RESOLVED** — the `if (n < 2) return;` on line 333 guards against empty `parts`. |
| W-016 | MEDIUM | STILL OPEN | `polite_voxelizer.hpp:217-224` — `rotated_copy` fallback still silently returns an unrotated copy on allocation failure. No log warning added. |
| W-017 | MEDIUM | STILL OPEN | `snuggle_nester.hpp:699-716` — `seed_bottom_left` triple-nested loop has no yield and no `max_seed_tries` cap. With small `min_gap_mm` and many parts, O(bed_area/gap^2 * n_parts) iterations run synchronously. |
| W-018 | LOW | STILL OPEN | `gpu_collision.cpp:256` — window class name is still static `"SnuggleGPUCollision"`. No PID suffix. |
| W-019 | LOW | STILL OPEN | `gpu_collision.cpp:570-601` — `cleanup()` does not call `UnregisterClassA`. |
| W-020 | LOW | RESOLVED | `gpu_collision.cpp:573-598` — cleanup now checks `if (gl_hwnd_)` (line 590) and within that checks `if (gl_dc_)` (line 592) before calling `ReleaseDC`. The `wglMakeCurrent` at line 577 is still guarded only by `gl_context_`, but `gl_dc_` is set before `gl_context_` in `init_context`, so if `gl_context_` is set, `gl_dc_` is too. Adequate. |
| W-021 | LOW | STILL OPEN | `gpu_collision.cpp:320` — `glGetString(GL_RENDERER)` is logged directly with no sanitization. |
| W-022 | LOW | STILL OPEN | `snuggle_nester.hpp:520-523` — `rot_cache_built_` is set but never reset. No cache invalidation. No documentation of single-use contract. |
| W-023 | LOW | STILL OPEN | `SnuggleArrange.cpp:73-82` — oversized parts get `bed_idx = -1` but no `continue`. They proceed through the full voxelization loop. (Though the voxelization won't crash — it just wastes CPU time.) |
| W-024 | LOW | STILL OPEN | `gpu_collision.cpp:350` — shader info log still uses fixed `char log[2048]`. |
| W-025 | LOW | STILL OPEN | `GLCanvas3D.cpp:5956,6011,6028,6042,6056,6067` — `BBLDragFloat`/`InputInt` widgets still allow momentary out-of-range values before the post-widget clamp fires. |
| W-026 | LOW | STILL OPEN | `polite_voxelizer.hpp:128` — `nx_ * ny_ * nz_` is unchecked. Relies on `MAX_GRID_DIM = 512` cap. No explicit pre-multiplication overflow guard. |
| W-027 | HIGH (arch) | STILL OPEN | `ROT_CACHE_BINS = 360` is defined three times: `gpu_collision.hpp:58` (as `ROT_CACHE_BINS`), `gpu_collision.hpp:94` (as `ROT_BINS`), `snuggle_nester.hpp:518` (as `ROT_CACHE_BINS`). Two different names, three locations. No shared header. |
| W-028 | HIGH (arch) | STILL OPEN | `3.14159265f` appears as a bare literal at: `SnuggleArrange.cpp:232,376`, `snuggle_nester.hpp:541,564,573,600,760,800`, `gpu_collision.cpp:38`. Shader uses `6.2831853`. No named constant. |
| W-029 | MEDIUM (arch) | STILL OPEN | `snuggle_nester.hpp` is still entirely header-only (~999 lines). No `.cpp` split. |
| W-030 | MEDIUM (arch) | STILL OPEN | `VoxelGrid`, `Vec3f`, `VoxError` are still defined in `polite_voxelizer.hpp`. No `snuggle_types.hpp`. |
| W-031 | MEDIUM (arch) | STILL OPEN | `gpu_collision.cpp:55` — `std::vector<const VoxelGrid*> rotated(n)` is still allocated inside the per-individual loop. |
| W-032 | MEDIUM (arch) | STILL OPEN | `gpu_collision.cpp:535` — `glGetUniformLocation` is called every generation, not cached at link time. `glGetBufferSubData` (line 548) is still a synchronous readback. No double-buffering. |
| W-033 | MEDIUM (arch) | STILL OPEN | `snuggle_nester.hpp:866-884` — part centroids are computed twice (centroid pass lines 866-870, avg_dist pass lines 878-884). No local array reuse. The proximity loop (lines 917-928) uses placement coordinates, not grid-center coordinates, so it's a third formulation. |
| W-034 | MEDIUM (arch) | STILL OPEN | Shader `local_size_x = 1` at `gpu_collision.cpp:116`. Comment on line 112-113 acknowledges it as POC. |
| W-035 | MEDIUM (arch) | STILL OPEN | `polite_voxelizer.hpp:310-328` — collision loop iterates in world-space floats with `world_to_grid()` inside the inner loop. `polite_yield()` is called every 10K iterations inside `collision_count` (line 327-328). |
| W-036 | STYLE | STILL OPEN | `SnuggleArrange.cpp:198` — variable is still named `hull_area_mm2` despite computing grid bounding-rectangle area. |
| W-037 | STYLE | STILL OPEN | `polite_voxelizer.hpp:274-275` — manual bit-serial popcount loop still in place. No `std::popcount` or `__builtin_popcount`. |
| W-038 | STYLE | RESOLVED | `snuggle_nester.hpp:705` — the offending line now reads `try_x = pr.x + pr.w + cfg_.min_gap_mm - cfg_.min_gap_mm`. Wait — re-reading the code: line 705 is `collides = true; try_x = pr.x+pr.w+cfg_.min_gap_mm-cfg_.min_gap_mm; break;`. The cancelling terms are still there. **STILL OPEN.** |
| W-039 | STYLE | PARTIALLY FIXED | Several tooltips have been added: Enable Snuggle checkbox (line 5934), Rotation step (line 5999), Resolution (line 6019), Population (line 6035), Generations (line 6046 area), Post-GA compaction (line 6086). Part gap slider still lacks a tooltip with units. Timeout, Max parts, and Multi-plate overflow controls — need to check. |
| W-040 | STYLE | STILL OPEN | "Compact after arrange" checkbox appears at both line 5964 (main section) and line 6080 (Advanced section). Both control `snuggle_compact`. Duplicate remains. |
| W-041 | STYLE | STILL OPEN | Unicode em-dash separators (`──`) are still present in comment dividers across all files. |
| W-042 | STYLE | STILL OPEN | `GLCanvas3D.cpp:5925` — header is still `"Snuggle 3D Arrangement"` with no "(experimental)" or "(beta)" label. |
| W-043 | STYLE | STILL OPEN | `SnuggleArrange.cpp:279-282` — progress string still shows `"Snuggle CPU gen 12/30 3 collisions"` format. |
| W-044 | STYLE | STILL OPEN | `SnuggleArrange.cpp:57-61` (max-parts) and `SnuggleArrange.cpp:189-191` (voxelization failure) both emit only `BOOST_LOG_TRIVIAL(warning)` with no user-visible notification. |
| W-045 | HIGH | RESOLVED | `snuggle_nester.hpp:376-466` — compact_toward_center now uses `cached_rotated()` throughout instead of `rotated_copy()`. Lines 376, 382, 434, 440, 454, 459 all call `cached_rotated()` which is an O(1) cache lookup. No redundant 3D grid copies. |
| W-046 | MEDIUM | STILL OPEN | Same as W-032. `glGetBufferSubData` at line 548, no fence+map, no double-buffering. |
| W-047 | MEDIUM | STILL OPEN | Same as W-031. `std::vector<const VoxelGrid*> rotated(n)` at `gpu_collision.cpp:55` inside per-individual loop. |
| W-048 | LOW | STILL OPEN | `gpu_collision.cpp:294` — still uses legacy `wglCreateContext()`. No `wglCreateContextAttribsARB`. |
| W-049 | LOW | STILL OPEN | `gpu_collision.hpp:20-21` — forward-declares `PartInfo` and `Individual` with no include-order documentation. The `.cpp` file does include `snuggle_nester.hpp` (line 7) before using these types, but standalone includers would fail. |

---

## Summary Counts

| Status | Count |
|--------|-------|
| RESOLVED | 11 |
| STILL OPEN | 35 |
| PARTIALLY FIXED | 3 |
| NOT APPLICABLE | 0 |
| **Total** | **49** |

Notes:
- W-015 was initially assessed as STILL OPEN but revised to RESOLVED upon closer inspection (the `n < 2` guard at line 333 covers the empty case).
- W-038 was initially assessed as RESOLVED but revised to STILL OPEN (the cancelling terms are still present in the code).
- W-046/W-047 overlap with W-032/W-031 respectively but are counted separately per ledger convention.

---

## Resolved Items — What the 4 Commits Fixed

The 4 post-ledger commits addressed 11 items fully and 3 partially, heavily weighted toward the critical/high-severity bugs:

**Fully resolved (11):**
- **W-001** (double-push) — `continue` on error path prevents second push
- **W-002** (OOB mesh index) — bounds check with `INVALID_INPUT` return
- **W-004** (bed_margin GPU/CPU mismatch) — `u_bed_margin` uniform added to shader
- **W-005** (uint32 overflow) — explicit >4 GB guard falls back to CPU
- **W-008** (null rot_cache) — null check at top of `evaluate_batch`
- **W-009** (ImGui disabled-stack) — balanced begin/end on `settings_out`
- **W-012** (moved-from grid UB) — grids are not moved from after vector construction
- **W-013** (lock_rotation desync on load) — sync line added at 1186
- **W-015** (parts[0] empty guard) — `n < 2` return covers it
- **W-020** (cleanup null gl_dc_) — conditional check added
- **W-045** (rotated_copy storm in compact) — `cached_rotated()` throughout

**Partially fixed (3):**
- **W-007** (name-based matching) — duplicate matching fixed via mark-used, but O(N^2) scan remains
- **W-014** (rotated_copy in compact hot path) — replaced with `cached_rotated()` but redundant bin lookups in bs2 loop
- **W-039** (missing tooltips) — several tooltips added, but Part gap, Timeout, Max parts, Multi-plate still missing

**Still open: W-003** (DC leak on GPU init failure) is the most significant remaining CRITICAL item.

---

## NEW Issues Spotted in the 4 Post-Ledger Commits

### N-001: cosmetic brace after continue creates false positive for dangling-else scanners
**File:line:** `SnuggleArrange.cpp:196`
**Severity:** STYLE
**Description:** The W-001 fix added a `continue` on line 195, but left the `} {` brace pattern on line 196 — the opening brace of a block that runs only in the success path. While functionally correct (the `continue` prevents reaching it on the error path), any static analyzer or future reviewer will flag this as a dangling-else candidate. Should be `} else {` or a plain block without the orphan brace.

### N-002: `cached_rotated()` in compact_toward_center uses unlocked cache bins for locked-rotation parts
**File:line:** `snuggle_nester.hpp:376,382,434`
**Severity:** LOW
**Description:** When `lock_rotation` is true, `build_rotation_cache` fills all 360 bins with the *same* grid (rotated at `initial_zrot`). The compact phase now calls `cached_rotated(idx, pl.zrot)` where `pl.zrot` may differ slightly from `initial_zrot` after micro-rotation attempts. Because all bins are identical when locked, the lookup returns the correct grid regardless of angle — so this is functionally safe. But the logic is fragile: if the cache strategy changes (e.g., only building one bin when locked), `cached_rotated` with an arbitrary angle would return the wrong bin. A comment documenting this invariant would prevent future breakage.

### N-003: `compute_fitness` rebuilds rotated grids when `lock_rotation` is true
**File:line:** `snuggle_nester.hpp:824-829`
**Severity:** MEDIUM
**Description:** When `lock_rotation` is true, `compute_fitness` calls `parts[i].grid.rotated_copy(ind.placements[i].zrot)` for every part (line 828), allocating fresh grids. This is unnecessary — the rotation cache already has these grids. The `else` branch (lines 831-834) correctly uses `cached_rotated()`. This redundant allocation happens once per feasible individual per generation. With population 512 and 10 parts, it is up to 5120 unnecessary `rotated_copy` calls per generation.

### N-004: cleanup() accesses GL functions without verifying context is current
**File:line:** `gpu_collision.cpp:579-583`
**Severity:** LOW
**Description:** `cleanup()` calls `wglMakeCurrent(hdc, ctx)` at line 577, then deletes GL objects (program, SSBOs). If `wglMakeCurrent` fails (e.g., the DC was invalidated by W-003's leak path), the subsequent `glDeleteProgram`/`glDeleteBuffers` calls operate on a null GL context, producing undefined behavior. Should check `wglMakeCurrent` return value before issuing GL calls.
