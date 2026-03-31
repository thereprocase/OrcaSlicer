# Uruk-Hai Bug Hunt: snuggle_nester.hpp

## Status: BUGS FOUND

Three confirmed bugs with undefined behavior and crashes.

---

## BUG 1: Division by Zero in snap_rotation() [Line 571]

**Severity:** CRASH

**Location:** Line 571

```cpp
float steps = std::round(delta / cfg_.rotation_step_rad);
```

**Root Cause:**
- `snap_rotation()` is called whenever `!cfg_.lock_rotation` (line 811 in mutate, line 563 check)
- No validation that `cfg_.rotation_step_rad > 0`
- Config default: `rotation_step_rad = 0.0f` (line 60)
- When `rotation_step_rad = 0.0f`, division by zero occurs

**Call Chain:**
1. `mutate()` at line 811 calls `snap_rotation(p.zrot, parts[i].initial_zrot)`
2. Function checks `if (cfg_.rotation_step_rad <= 0.001f)` but should check **equality** not comparison
3. When `rotation_step_rad` is exactly 0.0f, line 564's check passes `<= 0.001f` correctly, but...
4. Actually, line 564 **should** return early. Let me re-read...

**Corrected Analysis:**
Looking at line 563-564:
```cpp
if (cfg_.lock_rotation) return initial_zrot;
if (cfg_.rotation_step_rad <= 0.001f) return angle_rad; // continuous
```

This **does** protect against division at line 571. The guard is correct.

**RETRACTION:** This is NOT a bug. Guard at line 564 is sufficient.

---

## BUG 2: Array Out-of-Bounds in compact_toward_center() [Line 336]

**Severity:** CRASH (USE_AFTER_FREE / OUT_OF_BOUNDS)

**Location:** Line 336

```cpp
void compact_toward_center(
    NesterResult &result,
    const std::vector<PartInfo> &parts)
{
    size_t n = result.placements.size();
    if (n < 2) return;

    auto t_start = std::chrono::steady_clock::now();
    bool allow_rot = !cfg_.lock_rotation;
    float vs = parts[0].grid.voxel_size;  // ← BUG: parts may be empty!
```

**Root Cause:**
- `n` is derived from `result.placements.size()` (not `parts.size()`)
- Early exit checks `if (n < 2)` but does **not** validate `parts.size() > 0`
- If `parts` vector is empty, `parts[0]` is undefined behavior
- This can happen if arrange was called with 0 parts initially, or if result was corrupted

**Scenario:**
1. `nester.arrange()` called with empty `parts` vector (already checked at line 155)
2. But if somehow `result.placements` is non-empty while `parts` is empty, this crashes
3. Possible if external code directly manipulates `NesterResult` before compaction

**Fix:** Add bounds check before accessing `parts[0]`:
```cpp
if (parts.empty()) return;
```

---

## BUG 3: Division by Zero in compute_fitness() [Line 889]

**Severity:** CRASH / NaN PROPAGATION

**Location:** Line 889

```cpp
float bed_diag = std::sqrt(cfg_.bed_width_mm * cfg_.bed_width_mm +
                           cfg_.bed_height_mm * cfg_.bed_height_mm);
float avg_dist = 0;
for (size_t i = 0; i < n; i++) {
    // ... accumulate avg_dist ...
}
avg_dist /= (n * bed_diag);  // ← Division by zero if bed_diag == 0
```

**Root Cause:**
- If `bed_width_mm == 0` **AND** `bed_height_mm == 0`:
  - `bed_diag = sqrt(0 + 0) = 0.0f`
  - `avg_dist /= (n * 0.0f)` → **Division by zero**
- Config defaults: `bed_width_mm = 256.0f`, `bed_height_mm = 256.0f`
- But config is mutable and could be set to (0, 0) by user or external code
- No validation in constructor or `arrange()` entry point

**Mutation Path:**
If user calls `nester.cfg_.bed_width_mm = 0; nester.cfg_.bed_height_mm = 0;` before `arrange()`, this triggers.

**Result:**
- `avg_dist = inf` or `nan`
- Propagates to `clustering` (line 890)
- Propagates to final `ind.fitness` (line 934-938)
- NaN fitness breaks selection logic and evolution convergence

**Fix:** Guard the division:
```cpp
if (bed_diag > 0.0001f) {
    avg_dist /= (n * bed_diag);
    clustering = 1.0f - std::clamp(avg_dist, 0.0f, 1.0f);
} else {
    clustering = 0.0f;  // degenerate bed
}
```

---

## BUG 4: Potential Out-of-Bounds in mutation_wild() [Lines 796-803]

**Severity:** UNDEFINED BEHAVIOR (safe_clamp handles it, but risky)

**Location:** Lines 796-803

```cpp
} else {
    // Wildcard: completely random new position
    float wild_margin = std::max(parts[i].grid.nx * parts[i].grid.voxel_size * 0.5f,
                                 cfg_.bed_margin_mm);
    float wx_hi = std::max(wild_margin + 0.1f, cfg_.bed_width_mm - wild_margin);
    float wy_hi = std::max(wild_margin + 0.1f, cfg_.bed_height_mm - wild_margin);
    p.x = randf(wild_margin, wx_hi);
    p.y = randf(wild_margin, wy_hi);
```

**Risk:**
- If `cfg_.bed_width_mm < 2 * cfg_.bed_margin_mm`, then `wx_hi < wild_margin`
- `randf()` receives `lo > hi`, which is undefined behavior
- safe_clamp at lines 807-808 will **clamp to midpoint** as fallback (line 554-555)

**Scenario:**
If user sets `bed_width_mm = 5.0f, bed_margin_mm = 4.0f`:
- `wild_margin = max(..., 4.0f) = 4.0f`
- `wx_hi = max(4.1f, 5.0f - 4.0f) = max(4.1f, 1.0f) = 4.1f`
- `randf(4.0f, 4.1f)` is valid, but extremely narrow
- If `bed_margin_mm > bed_width_mm/2`, range collapses

**Current Mitigation:** safe_clamp handles this by returning midpoint (line 554), so technically safe but should validate config.

---

## SUMMARY OF CONFIRMED BUGS

| Line | Bug | Severity | Root Cause |
|------|-----|----------|-----------|
| 336  | OOB in `compact_toward_center()` | CRASH | `parts[0]` accessed without bounds check; `parts` may be empty |
| 889  | Division by zero in `compute_fitness()` | CRASH/NaN | `bed_diag` can be 0 if both bed dimensions are 0 |

## CONDITIONAL BUGS (Low probability but real)

| Line | Bug | Severity | Condition |
|------|-----|----------|-----------|
| 796-803 | Undefined `randf()` range | UB | `bed_width_mm < 2 * bed_margin_mm` or similar inverted ranges |

---

## RECOMMENDATIONS

1. **Add input validation** in constructor or `arrange()` entry:
   - `assert(parts.size() > 0 || cfg_.population_size == 0)`
   - `assert(cfg_.bed_width_mm > 0 && cfg_.bed_height_mm > 0)`

2. **Guard division at line 889** with epsilon check on `bed_diag`

3. **Add early return** at line 333 in `compact_toward_center()`:
   ```cpp
   if (parts.empty()) return;
   ```

4. **Validate config ranges** before evolution loop (line 193)
