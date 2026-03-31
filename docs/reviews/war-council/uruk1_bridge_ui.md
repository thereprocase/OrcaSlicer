# Uruk-Hai Bug Hunt — SnuggleArrange Bridge + GLCanvas3D UI

## SnuggleArrange.cpp

### BUG #1: Silent grid discard on voxelization failure (Line 196)

**Location:** Line 196
**Severity:** CRITICAL — silent data loss, downstream crash risk

```cpp
if (err != snuggle::VoxError::OK || pi.grid.count_solid() == 0) {
    // ... mark item off-plate, push stub PartInfo with empty grid
    parts.push_back(std::move(pi));
    continue;
} {  // ← DANGLING BRACE
    pi.max_height_mm = pi.grid.nz * pi.grid.voxel_size;
    pi.hull_area_mm2 = pi.grid.nx * pi.grid.ny * voxel_size * voxel_size;
    // ... log successful grid
}

parts.push_back(std::move(pi));  // ← PUSHED AGAIN (same broken pi)
```

**Problem:** The dangling brace on line 196 creates unreachable code. When voxelization fails, pi is pushed with `max_height_mm=0` and `hull_area_mm2=0`. The success branch never executes. But critically, the second `parts.push_back(std::move(pi))` at line 205 **pushes the same empty pi again** because `parts[i]` now refers to an already-moved object.

After the loop:
- `items` has the correct count (some entries marked `bed_idx=-1`)
- `parts` has **double-counted failed parts** (empty grids pushed twice)
- `result.placements.size()` matches `parts.size()`, not `items.size()`

**Downstream crash:** Lines 302, 318, 366 iterate `i < items.size() && i < result.placements.size()`. When items[k] failed voxelization, `parts.size() > items.size()`, so placement writes start misaligning. Part 0 placement applies to item 0, part 1 (which is a duplicate empty) tries to apply to item 1, etc. **Off-by-one index error cascades through result application.**

**Fix:** Delete the dangling brace. Restructure as:
```cpp
if (err != snuggle::VoxError::OK || pi.grid.count_solid() == 0) {
    BOOST_LOG_TRIVIAL(warning) << "Snuggle: voxelization failed for '" << pi.name << "'";
    items[i].bed_idx = -1;
    pi.max_height_mm = 0;
    pi.hull_area_mm2 = 0;
} else {
    pi.max_height_mm = pi.grid.nz * pi.grid.voxel_size;
    pi.hull_area_mm2 = pi.grid.nx * pi.grid.ny * voxel_size * voxel_size;
    BOOST_LOG_TRIVIAL(warning) << "Snuggle: [" << i << "] " << pi.name << "...";
}
parts.push_back(std::move(pi));
```

---

### BUG #2: Name-based mesh lookup with duplicate names (Lines 122-128)

**Location:** Lines 122–128
**Severity:** HIGH — silent data loss on duplicate names

```cpp
for (auto& oi : all_instances) {
    if (oi.obj && oi.obj->name == items[i].name) {
        obj = oi.obj;
        inst = oi.inst;
        oi.obj = nullptr;  // mark used
        break;
    }
}
```

**Problem:** Lookup by `oi.obj->name` (object name) against `items[i].name` (which is the *part* name computed by the arrange pipeline). If two instances of the *same* object are in the model, or if object names collide due to import/duplication, the first match is always taken. The `oi.obj = nullptr` sentinel marks it used, but:

1. If the same object is imported twice with identical names, the second instance never finds a match (first one is null'd out).
2. The logic assumes a stable correspondence between the order of `all_instances` (constructed by walking model.objects) and `items` (computed by `prepare_all()`). If part enumeration happens in a different order, names match but instances are wrong.

**Consequence:** Part with no matching instance is marked off-plate and skipped. If multiple parts should use different instances of the same object, only the first one gets processed correctly.

**Fix:** Validate that `items[i].name` truly corresponds to an *object* name (not a part identifier). If multipart objects exist, include instance index in the match (e.g., `name_instance_idx`). Or: store a persistent instance ID at ArrangePolygon creation time and carry it through the arrange pipeline.

---

### BUG #3: Grid index reuse after move (Line 323–324)

**Location:** Lines 323–341
**Severity:** MEDIUM — undefined behavior when greedy accept runs

```cpp
const auto& grid_i = parts[i].grid;  // ← parts[i] after move
if (grid_i.total_voxels() == 0) { items[i].bed_idx = -1; rejected++; continue; }
auto rot_i = grid_i.rotated_copy(pl.zrot);
// ...
for (size_t j : accepted) {
    const auto& pl_j = result.placements[j];
    auto rot_j = parts[j].grid.rotated_copy(pl_j.zrot);  // ← parts[j] UB if j > 0
```

**Problem:** After the voxelization loop (line 206), all `PartInfo` objects are moved into `parts` via `push_back`. In the greedy accept phase, `parts[i]` is read safely for the *current* part, but `parts[j]` for previously-accepted parts is accessed after their grids have been moved. If `grid` contains pointers or heap allocations, this is use-after-move.

If `VoxelGrid` contains a `std::unique_ptr` or similar, calling `rotated_copy(pl_j.zrot)` on a moved grid is undefined behavior (likely crash or garbled collision data).

**Consequence:** Collision detection fails intermittently. Parts reported as "no collision" actually collide. Or: crash on grid dereference.

**Fix:** Store grids immutably after voxelization. Reference them by index without move semantics. Or ensure `VoxelGrid::rotated_copy` is safe on moved grids.

---

## GLCanvas3D.cpp (snuggle section only)

### BUG #4: Slider value not clamped after direct drag input (Line 6011–6012)

**Location:** Lines 6011–6013
**Severity:** LOW — config may persist out-of-range values

```cpp
bool b_voxel_input = ImGui::BBLDragFloat("##voxel_input", &settings.snuggle_voxel_mm, 0.1f, 0.5f, 5.0f, "%.1f");
if (b_voxel || b_voxel_input) {
    settings.snuggle_voxel_mm = std::clamp(settings.snuggle_voxel_mm, 0.5f, 5.0f);
```

**Problem:** The drag operation uses min/max (0.5f, 5.0f) as *suggested* limits, not hard limits. User can drag outside [0.5, 5.0] before the clamp. The clamp is applied *after* UI feedback, but if the user clicks away before it syncs, `settings.snuggle_voxel_mm` may briefly hold an out-of-range value. When saved to config, it persists.

**Same issue at:**
- Line 5956 (`snuggle_padding_mm`, range 0.0–20.0)
- Line 6028 (`snuggle_population`, range 16–1024, uses InputInt)
- Line 6042 (`snuggle_generations`, range 10–500)
- Line 6056 (`snuggle_timeout_s`, range 2.0–300.0)
- Line 6067 (`snuggle_max_parts`, range 2–500)

**Fix:** Apply clamp *before* reading from the widget, or ensure BBLDragFloat / InputInt enforce hard limits internally.

---

### BUG #5: No disabled_end match for rotation step dropdown (Lines 5945–6092)

**Location:** Lines 5945–6092
**Severity:** MEDIUM — ImGui state corruption

```cpp
{
    if (!settings_out.use_snuggle) { imgui->disabled_begin(true); }  // Line 5946

    // ... many controls ...

    if (ImGui::TreeNode(...)) {
        // ... Advanced section ...
        ImGui::TreePop();  // Line 6089
    }

    if (!settings_out.use_snuggle) { imgui->disabled_end(); }  // Line 6092
}
```

**Problem:** If the user *opens* the TreeNode and `use_snuggle` is false:
1. `disabled_begin(true)` is called (line 5946)
2. TreeNode is skipped (disabled)
3. User checks/enables `use_snuggle` (line 5928)
4. Next frame: `disabled_end()` is not called because the condition on line 6092 is now false
5. All subsequent ImGui state is corrupted — remaining controls in the panel are disabled indefinitely

**Consequence:** If user enables Snuggle after opening the "Effort / Advanced" tree while it was disabled, the UI becomes unresponsive. Buttons, checkboxes downstream stop working.

**Fix:** Unconditionally match the disabled_begin/end:
```cpp
bool snuggle_disabled = !settings_out.use_snuggle;
if (snuggle_disabled) { imgui->disabled_begin(true); }

// ... all controls ...

if (snuggle_disabled) { imgui->disabled_end(); }
```

Or refactor the disabled scope to not span TreeNode interactions.

---

### BUG #6: Config load mismatch on snuggle_rotation_step (Line 1179)

**Location:** Lines 1177–1179 (load), 5992–5995 (save)
**Severity:** MEDIUM — config persistence broken

```cpp
// Load (GLCanvas3D.cpp:1179)
try { m_arrange_settings_fff.snuggle_rotation_step = std::clamp(std::stoi(snuggle_rot_step_str), 0, 90); } catch (...) {}

// Save (GLCanvas3D.cpp:5992)
appcfg->set("arrange", "snuggle_rotation_step", std::to_string(settings_out.snuggle_rotation_step));

// Dropdown sync (GLCanvas3D.cpp:5989–5997)
if (ImGui::Combo("##SnuggleRotStep", &rot_idx, rot_labels, 6)) {
    settings.snuggle_rotation_step = rot_values[rot_idx];
    settings_out.snuggle_rotation_step = settings.snuggle_rotation_step;
    appcfg->set("arrange", "snuggle_rotation_step", std::to_string(settings_out.snuggle_rotation_step));
    settings_out.snuggle_lock_rotation = (settings_out.snuggle_rotation_step == 0);
    appcfg->set("arrange", "snuggle_lock_rotation", settings_out.snuggle_lock_rotation ? "1" : "0");
}
```

**Problem:** Load validates clamp(val, 0, 90). The UI dropdown restricts to [0, 1, 5, 15, 45, 90]. But the config can be edited externally to (say) 37 degrees. On load:
1. snuggle_rot_step_str = "37"
2. Clamped to 37 (within [0, 90]) — passed
3. Passed to dropdown snap logic (lines 5979–5987)
4. Snapped to nearest valid [0, 1, 5, 15, 45, 90]

This is intentional, but there's **no verification that snuggle_lock_rotation persists correctly when rot_step is snapped**. If config has `snuggle_rotation_step=37` (invalid) and `snuggle_lock_rotation=0`:
- On load, rot_step snaps to (say) 45
- `snuggle_lock_rotation` is **not** re-derived; it stays 0 from config
- Line 5994 only sets lock_rotation when Combo changes, not on load
- After load, lock_rotation is **out of sync with rotation_step** (should be false if step != 0)

**Consequence:** User loads config → rotation_step is valid but lock_rotation doesn't match → GA uses wrong rotation mode.

**Fix:** After snapping rotation_step on load (line 1179), re-derive lock_rotation:
```cpp
m_arrange_settings_fff.snuggle_rotation_step = std::clamp(std::stoi(snuggle_rot_step_str), 0, 90);
m_arrange_settings_fff.snuggle_lock_rotation = (m_arrange_settings_fff.snuggle_rotation_step == 0);
```

---

## Summary

| Bug | File | Line(s) | Severity | Type |
|-----|------|---------|----------|------|
| #1: Dangling brace double-push | SnuggleArrange.cpp | 196, 205 | CRITICAL | Data corruption |
| #2: Name-based lookup ambiguity | SnuggleArrange.cpp | 122–128 | HIGH | Silent data loss |
| #3: Use-after-move on grid | SnuggleArrange.cpp | 323–341 | MEDIUM | UB, crash risk |
| #4: Slider clamp timing | GLCanvas3D.cpp | 5956, 6011–6013, 6028, 6042, 6056, 6067 | LOW | Config pollution |
| #5: disabled_begin/end mismatch | GLCanvas3D.cpp | 5945–6092 | MEDIUM | UI corruption |
| #6: lock_rotation not synced on load | GLCanvas3D.cpp | 1179, 5994 | MEDIUM | Config inconsistency |

