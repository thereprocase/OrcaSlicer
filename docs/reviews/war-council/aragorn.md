# War Council — Security Audit: Aragorn the Ranger

**Scope:** Snuggle subsystem, read-only analysis.
**Files reviewed:**
- `src/libslic3r/Arrange/polite_voxelizer.hpp`
- `src/libslic3r/Arrange/snuggle_nester.hpp`
- `src/libslic3r/Arrange/gpu_collision.cpp`
- `src/libslic3r/Arrange/SnuggleArrange.cpp`

**Threat model:** Crafted 3MF/mesh input that causes memory corruption, unbounded
allocation, resource exhaustion, or incorrect placement silently accepted as feasible.

---

## Findings by Severity

### CRITICAL

#### C1 — Integer overflow in `voxelize_indexed_mesh`: unchecked index into vertex array

**File:** `polite_voxelizer.hpp`, lines 552–558

```cpp
uint32_t i0 = indices[i * 3 + 0];
// ...
tris[i].v0 = {vertices[i0*3], vertices[i0*3+1], vertices[i0*3+2]};
```

The index values `i0`, `i1`, `i2` come directly from the mesh `indices` array. There is
**no bounds check** against `num_verts`. A crafted mesh can set `i0 >= num_verts`, causing
an out-of-bounds read of the `vertices` array. On a 32-bit index, `i0 * 3` can also
overflow `size_t` on a 32-bit address space (not the common case on 64-bit, but the
unchecked read is the primary issue regardless of pointer width).

**Attack:** Load a 3MF with `indices[k] = 0xFFFFFFFF`. This reads 12 bytes far past the
`verts` vector's allocation. Depending on heap layout, this can leak floats from adjacent
allocations or, if vertices is near a page boundary, segfault.

**Fix:** Add `if (i0 >= num_verts || i1 >= num_verts || i2 >= num_verts) { skip or error; }`
before accessing `vertices[i0*3]`.

---

#### C2 — Integer overflow in `VoxelGrid::allocate` memory size check

**File:** `polite_voxelizer.hpp`, lines 126–132

```cpp
size_t total = nx_ * ny_ * nz_;
if (total > MAX_TOTAL_VOXELS) return VoxError::GRID_TOO_LARGE;

size_t bytes = (total + 7) / 8;
if (bytes > MAX_VOXEL_MEMORY_MB * 1024 * 1024)
    return VoxError::MEMORY_CAP_EXCEEDED;
```

`MAX_VOXEL_MEMORY_MB * 1024 * 1024` is computed as integer arithmetic on `size_t`
constants. On a 64-bit target this is fine. But `nx_ * ny_ * nz_` is **also a
`size_t` multiplication with no overflow guard**. Each dimension is already checked
against `MAX_GRID_DIM = 512` before reaching this line, so `512 * 512 * 512 = 134M`,
which fits in 64-bit cleanly. This is safe today.

However, `MAX_GRID_DIM` is a `size_t` constant with no enforcement relationship to
`MAX_TOTAL_VOXELS`. If `MAX_GRID_DIM` is ever raised (e.g., to 1024 for higher-res
nesting) without updating the product overflow analysis, `1024^3 = 1G` still fits in
64-bit but `(1024*1024*1024 + 7) / 8 = 128MB` — which passes the bytes check at
`MAX_VOXEL_MEMORY_MB = 256`. That is fine, but the **order of checks creates a silent
dependency**: the MAX_GRID_DIM cap is what actually prevents overflow, not any explicit
product overflow guard. This is a latent maintenance hazard rather than an active
exploit, but worth naming.

---

### HIGH

#### H1 — No bounds check on `rot_cache_` access in `CpuCollisionEvaluator::get_rotated`

**File:** `gpu_collision.cpp`, lines 37–41

```cpp
int bin = (int)std::floor(angle * ROT_CACHE_BINS / (2.0f * 3.14159265f));
bin = ((bin % ROT_CACHE_BINS) + ROT_CACHE_BINS) % ROT_CACHE_BINS;
return (*rot_cache_)[part_idx][bin];
```

`rot_cache_` is set by `upload_grids` via a raw pointer. If `upload_grids` was never
called (or if `evaluator_` is used before the cache is built), `rot_cache_` is `nullptr`
and the dereference crashes. The early return in `evaluate_batch` (`if (!rot_cache_) return;`)
guards `evaluate_batch` itself, but `get_rotated` has no such guard. Any code path
that calls `get_rotated` directly (possible in subclass or test code) will segfault.

More importantly: `part_idx` is never checked against `rot_cache_->size()`. A mismatch
between the number of parts passed to `upload_grids` and the number passed to
`evaluate_batch` causes an out-of-bounds vector access. The caller in `snuggle_nester.hpp`
always passes the same `parts` vector, but nothing enforces this structurally.

---

#### H2 — `data_offset` in `GridMeta` is stored as `uint32_t` but total voxel bytes can exceed 4 GB

**File:** `gpu_collision.cpp`, lines 96–101 and 417

```cpp
struct alignas(4) GridMeta {
    uint32_t data_offset;  // byte offset into voxel SSBO
    // ...
};
// ...
meta_data[meta_idx].data_offset = (uint32_t)offset;
```

`offset` accumulates as `size_t`. With 100 parts × 360 rotation bins × 512^3 voxels
at 1-bit packing, total data can exceed `UINT32_MAX = 4GB`. The cast to `uint32_t`
silently wraps, causing the shader to read from the wrong offset — producing incorrect
collision results **without any error signal**. The nester would then report zero
collisions for genuinely overlapping parts, accepting an infeasible arrangement as
feasible.

The `MAX_VOXEL_MEMORY_MB = 256` cap per grid limits a single part to 256MB, but
`n_parts * ROT_BINS * 256MB` can easily exceed 4GB with modest part counts (>= 16 parts).
In practice the GPU SSBO size limit would OOM first, but the truncation happens
before the OOM check.

**Fix:** Use `uint64_t data_offset` in `GridMeta` and match the shader layout, or add
an explicit overflow check before the cast with a hard error if exceeded.

---

#### H3 — GL context window class name is a static string: registration collision

**File:** `gpu_collision.cpp`, lines 254–262

```cpp
wc.lpszClassName = "SnuggleGPUCollision";
if (!RegisterClassA(&wc)) {
    DWORD err = GetLastError();
    if (err != ERROR_CLASS_ALREADY_EXISTS) {
        return false;
    }
}
```

`ERROR_CLASS_ALREADY_EXISTS` is treated as non-fatal, which is correct for the common
multi-instance case. However: when the class already exists from a **different module**
with a different `lpfnWndProc`, `CreateWindowA` with the old class will silently use
the other proc. An adversarial DLL that registers "SnuggleGPUCollision" first with its
own `WndProc` would receive all messages sent to the GL helper window. This is a
window-handle squatting attack (low exploitability in a desktop slicer context, but
real in a shared process or plugin scenario).

**Fix:** Append the process ID or a GUID to the class name, or use `RegisterClassExA`
with `CS_CLASSDC` and verify the registered proc address before proceeding.

---

#### H4 — `GL_RENDERER` string logged without sanitization

**File:** `gpu_collision.cpp`, line 319

```cpp
BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: OpenGL " << major << "." << minor
                        << " (" << glGetString(GL_RENDERER) << ")";
```

`glGetString(GL_RENDERER)` returns a driver-supplied string. On a compromised or
spoofed driver (e.g., via MESA override, virtual machine, or a crafted OpenGL wrapper),
this string is untrusted and can contain log injection characters (`\n`, `\r`, ANSI
escape codes). Depending on the log consumer, this could corrupt log entries or
interfere with log-parsing tools.

Low severity on its own, but it is an untrusted data path into a formatted log stream.

---

### MEDIUM

#### M1 — Floating-point loop termination in `collision_count` can become infinite

**File:** `polite_voxelizer.hpp`, lines 308–329

```cpp
float vs = std::max(a.voxel_size, b.voxel_size);
for (float wz = oz_min + vs * 0.5f; wz < oz_max; wz += vs) {
    for (float wy = oy_min + vs * 0.5f; wy < oy_max; wy += vs) {
        for (float wx = ox_min + vs * 0.5f; wx < ox_max; wx += vs) {
```

If `vs` is zero or subnormal (from a zero-size voxel grid somehow reaching this path),
the loops never terminate. The same pattern is replicated in the GLSL shader (lines
206–208), where GPU divergence would hang the workgroup. The voxelizer rejects
`voxel_size_mm <= 0.01f` at entry, so in practice `vs` is safe after normal
voxelization — but `rotated_copy` can return a grid with `voxel_size` inherited from
a corrupted source, and `VoxelGrid`'s default constructor leaves `voxel_size = 1.0f`
(safe). The risk is low but worth noting because the failure mode is a hang, not an error.

**Fix:** Assert or clamp `vs > 0` before the loops, return 0 early if not.

---

#### M2 — `snap_rotation` uses `while` loops that can be infinite for extreme inputs

**File:** `snuggle_nester.hpp`, lines 567–570

```cpp
while (delta < 0) delta += TWO_PI;
while (delta >= TWO_PI) delta -= TWO_PI;
```

For a `delta` that is NaN, both conditions are false and the loops terminate
immediately — safe. For a finite but extremely large value (e.g., `1e18f`), the
loop runs ~2.4×10^17 iterations. This is DoS via crafted `initial_zrot` input.
`initial_zrot` comes from `items[i].rotation` in `SnuggleArrange.cpp` line 116, which
comes from the ArrangePolygon's rotation field. A 3MF file can set this to any float.

**Fix:** Replace the `while` loops with `fmod`: `delta = fmodf(delta, TWO_PI); if (delta < 0) delta += TWO_PI;`.

---

#### M3 — `seed_bottom_left` inner loop lacks a yield and can stall on large part counts

**File:** `snuggle_nester.hpp`, lines 702–712

```cpp
for (float try_y = cfg_.bed_margin_mm; try_y < cfg_.bed_height_mm - ph; try_y += cfg_.min_gap_mm) {
    for (float try_x = cfg_.bed_margin_mm; try_x < cfg_.bed_width_mm - pw; try_x += cfg_.min_gap_mm) {
        for (const auto &pr : placed) { ... }
```

Triple nested loop: O(bed_area / gap^2 × n_parts). With `min_gap_mm = 1.0f` and a
256mm bed, the outer two loops alone reach ~65K iterations per part. With 50 parts, the
O(n²) `placed` check makes this ~160M iterations for the last part before it gives up.
There is no timeout check and no yield inside this loop. This is a seeding function, not
the main GA, but it runs synchronously on the UI thread before the GA starts.

A crafted scenario: 50 tiny parts, `min_gap_mm = 0.1mm`, no yield, main thread stalls
for several seconds before the GA timeout even starts counting.

**Fix:** Add a yield every N outer iterations, and cap `max_tries` to a bounded count.

---

#### M4 — Name-based model lookup is O(n²) and exploitable for DoS

**File:** `SnuggleArrange.cpp`, lines 120–128

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

For N items, this is O(N²) string comparisons. A 3MF with 1000 identically-named
objects triggers worst-case matching: every item scans the entire `all_instances`
list before finding or failing to find a match. String comparison on long names
amplifies the cost. Not a crash, but a noticeable stall before the GA even starts,
which could be used to defeat the timeout enforcement.

**Fix:** Build a `std::multimap<std::string, ObjInst*>` keyed by name once, then do
O(log N) lookups. Alternatively cap the search after the first miss.

---

#### M5 — GPU `evaluate_batch` ignores `bed_margin` entirely

**File:** `gpu_collision.cpp`, line 477

```cpp
// TODO: pass bed_margin to shader as uniform for GPU OOB check
(void)bed_margin;
```

The GPU path computes OOB against raw `u_bed_w / u_bed_h` with zero margin. The CPU
path applies `bed_margin` correctly. When the GPU backend is active, parts can be
placed within the margin zone and the evaluator reports them as feasible. The GA
evolves toward arrangements that the CPU would call infeasible. Result: parts placed
closer to the bed edge than `min_gap_mm` allows.

This is correctness-affecting, not memory-safety, but it means the GPU path's
"feasible" signal is unreliable. The mismatch is the bug.

---

### LOW

#### L1 — GL context cleanup does not unregister the window class

**File:** `gpu_collision.cpp`, `cleanup()`, lines 549–580

`DestroyWindow` and `wglDeleteContext` are called, but `UnregisterClassA` is never
called. In a long-running process that creates and destroys `GpuCollisionEvaluator`
multiple times (e.g., per-arrange-job), the class registration leaks. The
`ERROR_CLASS_ALREADY_EXISTS` guard in `init_context` masks re-registration failures
forever after the first instance. This is a resource leak, not a safety issue.

---

#### L2 — Rotation cache built only once, not invalidated on config change

**File:** `snuggle_nester.hpp`, line 526

```cpp
void build_rotation_cache(const std::vector<PartInfo> &parts) {
    if (rot_cache_built_) return;
```

`rot_cache_built_` is never reset. If a `SnuggleNester` instance is reused with
different parts (or a different `lock_rotation` setting), the stale cache is used
silently. The current call site in `SnuggleArrange.cpp` creates a fresh nester per
arrange invocation, so this is not exploitable today — but it is a correctness time
bomb if the nester is ever pooled or reused.

---

#### L3 — Shader info log uses a fixed 2048-byte stack buffer; driver output may be truncated

**File:** `gpu_collision.cpp`, lines 349–352

```cpp
char log[2048] = {};
glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
```

Driver error messages from verbose compilers (e.g., Mesa) can exceed 2KB. The truncation
is silent — the logged warning may be misleading about the actual compile error.
Not a security issue; a debuggability issue.

---

## Summary Table

| ID | Severity | Where | Impact |
|----|----------|-------|--------|
| C1 | Critical | `voxelize_indexed_mesh` | OOB read from crafted mesh index |
| C2 | Critical (latent) | `VoxelGrid::allocate` | Overflow dependency on MAX_GRID_DIM cap |
| H1 | High | `CpuCollisionEvaluator::get_rotated` | Null deref / OOB if cache not initialized |
| H2 | High | `GridMeta::data_offset` | Silent uint32 truncation → wrong collision results |
| H3 | High | `init_context` | Window class squatting in shared-process context |
| H4 | High | `init_context` | Log injection from driver-supplied GL_RENDERER string |
| M1 | Medium | `collision_count` + shader | Infinite loop if voxel_size is zero |
| M2 | Medium | `snap_rotation` | Infinite loop from crafted extreme rotation float |
| M3 | Medium | `seed_bottom_left` | UI thread stall, O(n² × bed_area) without yield |
| M4 | Medium | `SnuggleArrange.cpp` | O(n²) name lookup, DoS with many same-named objects |
| M5 | Medium | GPU `evaluate_batch` | Margin ignored on GPU path; incorrect feasibility |
| L1 | Low | `cleanup()` | Window class not unregistered — resource leak |
| L2 | Low | `build_rotation_cache` | Stale cache if nester reused with different parts |
| L3 | Low | `compile_shader` | Truncated error log from verbose GPU drivers |

---

## Priority Recommendation

Fix **C1** first — it is an unconditional OOB read reachable from any 3MF with a
bad index, requires no special conditions, and is straightforward to fix with two
lines of bounds checking.

Fix **H2** next — the silent `uint32_t` truncation corrupts collision detection
without any error signal. The GPU path reports infeasible arrangements as feasible,
which undermines the entire safety guarantee of the nester.

Fix **M2** and **M5** in the same pass — both are one-liners and both affect
correctness of the output, not just memory safety.
