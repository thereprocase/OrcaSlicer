# War Council — Uruk-Hai Bug Review
## File: polite_voxelizer.hpp
## Reviewer: Uruk-Hai (Adversarial Bug Hunter)

---

## CONFIRMED BUGS

### BUG #1: Division by Zero in `world_to_grid()`
**Location:** Lines 161-165
**Severity:** HIGH — Undefined behavior → crash or garbage results

```cpp
void world_to_grid(const Vec3f &world, int &gx, int &gy, int &gz) const {
    gx = (int)std::floor((world.x - origin.x) / voxel_size);  // ← Division by zero if voxel_size == 0
    gy = (int)std::floor((world.y - origin.y) / voxel_size);
    gz = (int)std::floor((world.z - origin.z) / voxel_size);
}
```

**Callsites:**
- Line 313: `a.world_to_grid()` in `collision_count()`
- Line 319: `b.world_to_grid()` in `collision_count()`
- Lines 471, 475, 479: in `voxelize_mesh()`

**Root cause:** No guard against `voxel_size <= 0`. If grid initialization fails silently or voxel_size is corrupted, division by zero produces NaN/Inf, then cast to int is undefined behavior.

**Fix:** Add at entry:
```cpp
if (voxel_size <= 0.0f) return 0;  // return sentinel
```

---

### BUG #2: Infinite Loop in `collision_count()` with Zero Voxel Size
**Location:** Lines 304, 308-310
**Severity:** HIGH — Infinite loop if grids have voxel_size == 0

```cpp
float vs = std::max(a.voxel_size, b.voxel_size);  // If both are 0, vs = 0
// ...
for (float wz = oz_min + vs * 0.5f; wz < oz_max; wz += vs) {
    for (float wy = oy_min + vs * 0.5f; wy < oy_max; wy += vs) {
        for (float wx = ox_min + vs * 0.5f; wx < ox_max; wx += vs) {
            // ...
            a.world_to_grid({wx - a_offset.x, ...}, ax, ay, az);  // Divides by zero here too
```

**Condition:** Both grids have `voxel_size == 0.0f` and overlap region is non-empty (`ox_min < ox_max`).

**Fix:** Guard voxel_size in `world_to_grid()` (see Bug #1), which also prevents this path. Alternatively, add guard before loop:
```cpp
if (vs <= 0.0f) return 0;
```

---

### BUG #3: Integer Overflow Risk in `allocate()` — Fragile Guard
**Location:** Lines 121-128
**Severity:** MEDIUM — Current guards are sufficient, but pattern is fragile

```cpp
VoxError allocate(size_t nx_, size_t ny_, size_t nz_) {
    // Bounds check — per-dimension only
    if (nx_ > MAX_GRID_DIM || ny_ > MAX_GRID_DIM || nz_ > MAX_GRID_DIM)
        return VoxError::GRID_TOO_LARGE;

    size_t total = nx_ * ny_ * nz_;  // ← Can overflow if guard at line 123 is removed
    if (total > MAX_TOTAL_VOXELS)
        return VoxError::GRID_TOO_LARGE;
```

**Current status:** Safe because `nx_, ny_, nz_ <= 512` implies product `<= 512^3 = 134M < MAX_TOTAL_VOXELS (64M)` check will catch oversized products.

**Fragility:** If future code removes the per-dimension checks (lines 123-124) to relax constraints, the multiplication becomes vulnerable to wrap-around. Example: `nx_=1024, ny_=1024, nz_1024` overflows `size_t` and wraps to small value, bypassing line 127 check.

**Fix:** Add pre-multiplication guard:
```cpp
if (nx_ > 0 && ny_ > 0 && nz_ > 0 &&
    nx_ > MAX_TOTAL_VOXELS / ny_ ||
    nx_ * ny_ > MAX_TOTAL_VOXELS / nz_) {
    return VoxError::GRID_TOO_LARGE;
}
size_t total = nx_ * ny_ * nz_;
```

---

## CLEAN FINDINGS

### `get()` and `set()` Bounds Checking — Lines 147-157
**Status:** CLEAN

Both functions have explicit bounds checks that return safely on out-of-bounds access. The bitwise operations (`idx >> 3`, `idx & 7`) are safe.

### `rotated_copy()` Edge Cases — Lines 179-261

**Zero rotation (angle < 1e-6):**
Status: CLEAN. Creates exact copy of all fields.

**360-degree rotation (angle ≈ 2π):**
Status: CLEAN. Math is symmetric; AABB computation takes min/max over all corners.

**Negative angles:**
Status: CLEAN. Trigonometric operations are symmetric; AABB min/max guards preserve correctness.

**Note:** All three cases inherit the `voxel_size` division risk from lines 210-211 if `voxel_size == 0`, but this depends on Bug #1 remaining unpatched.

### `collision_count()` with Empty Grids
**Status:** CLEAN

Empty grid detection works correctly: if `nx = 0` or `ny = 0`, bounds collapse and overlap check at line 300 returns 0 safely. Assumes `voxel_size > 0` (Bug #1).

---

## SEVERITY RANKING

1. **BUG #1 (Division by Zero in `world_to_grid()`)**: HIGH — Affects multiple call sites, crashes likely
2. **BUG #2 (Infinite Loop in `collision_count()`)**: HIGH — Depends on Bug #1 + zero voxel_size condition
3. **BUG #3 (Integer Overflow Fragility)**: MEDIUM — Currently mitigated by guards, but unsafe pattern

---

## RECOMMENDATIONS

- **Immediate:** Add `voxel_size` validation guard to `world_to_grid()` (fixes Bugs #1 and #2)
- **Short-term:** Add pre-multiplication overflow check to `allocate()` (fixes Bug #3)
- **Testing:** Add fuzz tests with zero/negative voxel_size, zero grid dimensions, and extreme overlap cases

---

**Reviewed by:** Uruk-Hai
**Date:** 2026-03-31
**Code Quality:** 7/10 — Politeness guarantees well-documented, but assumes caller honesty (trusts voxel_size > 0)
