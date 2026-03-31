# Snuggle GPU Collision: Design Document

GPU-accelerated voxel collision evaluation for the Snuggle 3D genetic nester. Replaces the CPU inner loop in `snuggle_nester.hpp::evaluate()` with OpenGL 4.3 compute shaders.

## Status Quo

The CPU path in `evaluate()` does this for each of `pop_size` individuals:

1. Look up N pre-rotated voxel grids from the rotation cache (360 bins per part)
2. For each of `N*(N-1)/2` pairs: compute world-space overlap box, iterate voxels in that box, AND both grids, sum collisions
3. For each of N parts: check bounding box against bed edges, accumulate OOB penalty
4. If feasible: compute bounding-box fitness

With `pop_size=256` and `N=15`, that is 256 x 105 = 26,880 pairwise collision checks per generation, each iterating ~1000 voxels in the overlap region. The CPU spends nearly all its time here.

---

## Architecture

### Data Layout: SSBOs, Not 3D Textures

3D textures would require unpacking bit-packed data into RGBA voxels (8x memory inflation) and can't exploit the existing uint32 bitwise AND we need. SSBOs let us upload the raw `uint8_t` bit-packed grids and operate on them with integer arithmetic in the shader.

**Voxel Data SSBO** — one large buffer holding every rotation variant of every part, packed end-to-end:

```
[part_0_rot_0][part_0_rot_1]...[part_0_rot_359][part_1_rot_0]...
```

Each grid entry is `ceil(nx*ny*nz / 8)` bytes of bit-packed voxel data, byte-aligned.

**Grid Metadata SSBO** — one struct per (part, rotation_bin):

```glsl
struct GridMeta {
    uint   data_offset;   // byte offset into voxel data SSBO
    uint   nx, ny, nz;    // grid dimensions in voxels
    float  voxel_size;    // mm per voxel
    float  origin_x;      // world-space origin of this rotated grid
    float  origin_y;
    float  origin_z;
    uint   _pad;           // align to 48 bytes
};
```

Total: `n_parts * 360` entries. At 48 bytes each, 15 parts = 259,200 bytes (~253 KB).

**Placement SSBO** — the GA's current population:

```glsl
struct Placement {
    float x, y, zrot;     // 12 bytes per part
    uint  _pad;            // align to 16 bytes
};
```

Layout: `pop_size * n_parts` entries, written by CPU each generation.

**Results SSBO** — one entry per individual:

```glsl
struct Result {
    uint collision_count;
    uint oob_count;
};
```

Layout: `pop_size` entries, read back by CPU after dispatch.

### Why Not Texture Arrays or Atlases

The rotation cache grids vary in dimension (rotation changes the bounding box). A texture array requires uniform dimensions per layer. An atlas requires packing logic and wastes memory on padding. SSBOs with per-grid metadata offsets handle variable sizes naturally and match the bit-packed representation we already have.

---

## Memory Budget

| Resource | Size | Notes |
|---|---|---|
| Voxel data | ~84 MB | 15 parts x 360 rotations x ~16 KB avg per grid |
| Grid metadata | ~253 KB | 15 x 360 x 48 bytes |
| Placements | ~61 KB | 256 x 15 x 16 bytes |
| Results | 2 KB | 256 x 8 bytes |
| **Total** | **~85 MB** | Fits any GPU with >= 256 MB VRAM |

For larger part counts, memory scales linearly. 30 parts would be ~170 MB. The 360-bin rotation cache is the dominant cost. If memory becomes a concern, we can reduce to 90 bins (4-degree resolution) for a 4x reduction, though the CPU cache already uses 360 bins so this would mean divergent behavior.

---

## Compute Shader Design

### Work Group Layout

One work group per individual. Each work group handles all `N*(N-1)/2` pairwise checks plus N bounds checks for its individual.

```
Dispatch: glDispatchCompute(pop_size, 1, 1)
Work group size: (N_PAIRS, 1, 1) where N_PAIRS = N*(N-1)/2
```

Problem: N_PAIRS can be up to 105 (for 15 parts), which is fine for a work group. But each invocation within the group still needs to iterate ~1000 voxels, which is serial per-invocation work.

Better approach — **two-level parallelism**:

```
Work group size: (128, 1, 1)  // or 64/256 depending on occupancy tuning
gl_WorkGroupID.x  = individual index [0..pop_size)
gl_LocalInvocationID.x = work item within this individual
```

Each work group processes one individual. The 128 invocations split the N*(N-1)/2 pair checks among themselves. Each invocation accumulates a local collision count, then we use `atomicAdd` to a shared variable, and the first invocation writes the final result.

For bounds checking: invocations 0..N-1 each handle one part's bounds check.

### Shader Pseudocode

```glsl
#version 430

layout(local_size_x = 128) in;

// ── Bindings ──────────────────────────────────────────────
layout(std430, binding = 0) readonly buffer VoxelData {
    uint voxels[];  // bit-packed, accessed as uint32 words
};

layout(std430, binding = 1) readonly buffer GridMetadata {
    // Per (part, rotation_bin): 12 uints = 48 bytes
    // [data_offset, nx, ny, nz, voxel_size_bits, ox, oy, oz, pad...]
    uint grid_meta[];
};

layout(std430, binding = 2) readonly buffer Placements {
    // Per (individual, part): x, y, zrot, pad
    vec4 placements[];
};

layout(std430, binding = 3) writeonly buffer Results {
    uvec2 results[];  // .x = collision_count, .y = oob_count
};

uniform uint u_n_parts;
uniform uint u_n_pairs;       // n_parts * (n_parts - 1) / 2
uniform uint u_pop_size;
uniform float u_bed_width;
uniform float u_bed_height;

// ── Shared memory for reduction ───────────────────────────
shared uint s_collision_count;
shared uint s_oob_count;

// ── Helper: decode grid metadata ──────────────────────────
struct GMeta {
    uint  data_offset;
    uint  nx, ny, nz;
    float voxel_size;
    float ox, oy, oz;
};

GMeta load_meta(uint part_idx, uint rot_bin) {
    uint base = (part_idx * 360u + rot_bin) * 12u;
    GMeta m;
    m.data_offset = grid_meta[base + 0];
    m.nx          = grid_meta[base + 1];
    m.ny          = grid_meta[base + 2];
    m.nz          = grid_meta[base + 3];
    m.voxel_size  = uintBitsToFloat(grid_meta[base + 4]);
    m.ox          = uintBitsToFloat(grid_meta[base + 5]);
    m.oy          = uintBitsToFloat(grid_meta[base + 6]);
    m.oz          = uintBitsToFloat(grid_meta[base + 7]);
    return m;
}

// ── Helper: read one bit from voxel data ──────────────────
bool voxel_get(GMeta m, int gx, int gy, int gz) {
    if (gx < 0 || gy < 0 || gz < 0) return false;
    if (uint(gx) >= m.nx || uint(gy) >= m.ny || uint(gz) >= m.nz) return false;
    uint idx = uint(gx) + uint(gy) * m.nx + uint(gz) * m.nx * m.ny;
    uint word_offset = m.data_offset / 4u + idx / 32u;
    uint bit = idx % 32u;
    return (voxels[word_offset] & (1u << bit)) != 0u;
}

// ── Helper: angle to rotation bin ─────────────────────────
uint angle_to_bin(float zrot) {
    int bin = int(floor(zrot * 360.0 / 6.28318530718));
    return uint(((bin % 360) + 360) % 360);
}

// ── Helper: pair index (i,j) from flat pair ID ────────────
//    Pairs enumerated as: (0,1),(0,2),...,(0,N-1),(1,2),...
void pair_from_flat(uint flat, uint n, out uint i, out uint j) {
    // Inverse of the triangular number: i = floor((2N-1 - sqrt((2N-1)^2 - 8*flat)) / 2)
    // Simpler: just iterate. N <= 30 so this is trivial.
    i = 0; j = 0;
    uint acc = 0;
    for (i = 0; i < n - 1u; i++) {
        uint row_len = n - 1u - i;
        if (flat < acc + row_len) {
            j = i + 1u + (flat - acc);
            return;
        }
        acc += row_len;
    }
}

void main() {
    uint ind_idx = gl_WorkGroupID.x;
    uint local_id = gl_LocalInvocationID.x;
    uint local_size = gl_WorkGroupSize.x;

    // Initialize shared counters
    if (local_id == 0u) {
        s_collision_count = 0u;
        s_oob_count = 0u;
    }
    barrier();

    // ── Phase 1: Pairwise collision checks ────────────────
    // Each invocation handles ceil(n_pairs / local_size) pairs
    for (uint pair_id = local_id; pair_id < u_n_pairs; pair_id += local_size) {
        uint pi, pj;
        pair_from_flat(pair_id, u_n_parts, pi, pj);

        // Load placements for this individual
        uint base_i = ind_idx * u_n_parts + pi;
        uint base_j = ind_idx * u_n_parts + pj;
        vec4 pl_i = placements[base_i];
        vec4 pl_j = placements[base_j];

        // Get rotation bins
        uint rot_i = angle_to_bin(pl_i.z);
        uint rot_j = angle_to_bin(pl_j.z);

        // Load grid metadata
        GMeta mi = load_meta(pi, rot_i);
        GMeta mj = load_meta(pj, rot_j);

        // World-space bounds for each part
        float ai_min_x = mi.ox + pl_i.x;
        float ai_min_y = mi.oy + pl_i.y;
        float ai_min_z = mi.oz;
        float ai_max_x = ai_min_x + float(mi.nx) * mi.voxel_size;
        float ai_max_y = ai_min_y + float(mi.ny) * mi.voxel_size;
        float ai_max_z = ai_min_z + float(mi.nz) * mi.voxel_size;

        float aj_min_x = mj.ox + pl_j.x;
        float aj_min_y = mj.oy + pl_j.y;
        float aj_min_z = mj.oz;
        float aj_max_x = aj_min_x + float(mj.nx) * mj.voxel_size;
        float aj_max_y = aj_min_y + float(mj.ny) * mj.voxel_size;
        float aj_max_z = aj_min_z + float(mj.nz) * mj.voxel_size;

        // Overlap box
        float ox_min = max(ai_min_x, aj_min_x);
        float oy_min = max(ai_min_y, aj_min_y);
        float oz_min = max(ai_min_z, aj_min_z);
        float ox_max = min(ai_max_x, aj_max_x);
        float oy_max = min(ai_max_y, aj_max_y);
        float oz_max = min(ai_max_z, aj_max_z);

        if (ox_min < ox_max && oy_min < oy_max && oz_min < oz_max) {
            float vs = max(mi.voxel_size, mj.voxel_size);
            uint local_collisions = 0u;

            for (float wz = oz_min + vs * 0.5; wz < oz_max; wz += vs) {
                for (float wy = oy_min + vs * 0.5; wy < oy_max; wy += vs) {
                    for (float wx = ox_min + vs * 0.5; wx < ox_max; wx += vs) {
                        int ax = int(floor((wx - pl_i.x - mi.ox) / mi.voxel_size));
                        int ay = int(floor((wy - pl_i.y - mi.oy) / mi.voxel_size));
                        int az = int(floor((wz - mi.oz) / mi.voxel_size));

                        if (voxel_get(mi, ax, ay, az)) {
                            int bx = int(floor((wx - pl_j.x - mj.ox) / mj.voxel_size));
                            int by = int(floor((wy - pl_j.y - mj.oy) / mj.voxel_size));
                            int bz = int(floor((wz - mj.oz) / mj.voxel_size));

                            if (voxel_get(mj, bx, by, bz)) {
                                local_collisions++;
                            }
                        }
                    }
                }
            }

            if (local_collisions > 0u) {
                atomicAdd(s_collision_count, local_collisions);
            }
        }
    }

    // ── Phase 2: Bounds checks ────────────────────────────
    // First N invocations each check one part
    if (local_id < u_n_parts) {
        uint base = ind_idx * u_n_parts + local_id;
        vec4 pl = placements[base];
        uint rot = angle_to_bin(pl.z);
        GMeta m = load_meta(local_id, rot);

        float part_min_x = pl.x + m.ox;
        float part_min_y = pl.y + m.oy;
        float part_max_x = part_min_x + float(m.nx) * m.voxel_size;
        float part_max_y = part_min_y + float(m.ny) * m.voxel_size;

        uint oob = 0u;
        if (part_min_x < 0.0) oob += uint(-part_min_x / m.voxel_size);
        if (part_min_y < 0.0) oob += uint(-part_min_y / m.voxel_size);
        if (part_max_x > u_bed_width) oob += uint((part_max_x - u_bed_width) / m.voxel_size);
        if (part_max_y > u_bed_height) oob += uint((part_max_y - u_bed_height) / m.voxel_size);

        if (oob > 0u) {
            atomicAdd(s_oob_count, oob);
        }
    }

    barrier();

    // ── Write result ──────────────────────────────────────
    if (local_id == 0u) {
        results[ind_idx] = uvec2(s_collision_count, s_oob_count);
    }
}
```

### Optimization: Word-Level AND

The per-bit `voxel_get` in the inner loop is the obvious first target. When two grids share the same voxel size and their origins align to voxel boundaries, we can AND entire uint32 words instead of individual bits — a 32x throughput improvement on the inner loop.

This requires pre-computing per-pair alignment at upload time, or a runtime check in the shader. Worth doing as a second pass once the naive version works.

---

## Dispatch Strategy

**One dispatch per generation.** The compute shader evaluates all `pop_size` individuals in parallel:

```
glDispatchCompute(pop_size, 1, 1);  // 256 work groups
glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
```

After dispatch, map the results SSBO and read back `pop_size * 2` uints. The CPU then does fitness scoring (bounding box calculation) only for feasible individuals — this is cheap and doesn't benefit from GPU acceleration.

**Sync model:**

```
CPU: write placements -> glBufferSubData
GPU: dispatch collision shader
CPU: glGetBufferSubData results (implicit sync — blocks until dispatch completes)
CPU: score fitness for feasible individuals
CPU: selection, crossover, mutation
CPU: write next generation's placements
repeat
```

The `glGetBufferSubData` call provides the synchronization point. No explicit fence or glFinish needed.

---

## Performance Estimate

### Theoretical

- Modern GPU: ~5000 shader cores at ~1.5 GHz
- Each collision check: ~1000 voxel AND operations
- Total work: 256 individuals x 105 pairs x 1000 ops = 26.9M operations
- At 5000 cores: ~5,380 ops per core = ~3.6 microseconds
- Per-generation GPU time: **~4 microseconds** (compute-bound estimate)

### Realistic

Memory latency dominates. Each voxel AND requires two random-access reads from the voxel SSBO. L2 cache on a modern GPU is 4-6 MB, which fits ~4 parts' rotation grids. Cache misses will dominate for the other 11 parts.

- SSBO random access latency: ~300-500 cycles per miss
- Expected cache hit rate for 15-part problem: ~30% (working set ~84 MB >> L2 size)
- Effective ops per pair: ~1000 x (0.3 x 4 + 0.7 x 400) = ~281,200 cycles
- Per individual: 105 pairs x 281K cycles / 5000 cores = ~5,900 cycles = ~4 microseconds
- But: work groups don't perfectly distribute. Realistic occupancy ~60%.
- **Expected per-generation: ~100-500 microseconds**

Even at the pessimistic end, this is 200-1000x faster than the CPU path, which runs ~100ms per generation for 15 parts at pop_size=256.

### Upload/Readback Overhead

- Placement upload: 61 KB via `glBufferSubData` — ~10 microseconds
- Results readback: 2 KB via `glGetBufferSubData` — ~5 microseconds
- Voxel data upload: 84 MB, **once** at nester startup — ~20ms over PCIe 3.0

Upload/readback is negligible per generation. The one-time voxel upload is amortized across 50-100 generations.

---

## The GL Context Thread Problem

This is the hard part. OpenGL contexts are thread-affine — a context can only be current on one thread at a time. OrcaSlicer's rendering context lives on the main (UI) thread. The arrange job runs on a Boost worker thread (`BoostThreadWorker`).

### Option A: Dedicated Compute Context on Worker Thread (Recommended)

Create a second OpenGL context on the worker thread, sharing resources with the main render context via `wglShareLists` (Windows) or equivalent.

**How it works:**

1. At startup (main thread): create a hidden `wxGLCanvas` + second `wxGLContext` with `shareWith` pointing to the existing render context
2. Before the arrange job launches: pass the second context and hidden canvas to the nester
3. On the worker thread (first call): `wglMakeCurrent(hidden_hdc, compute_context)`, then load compute shader, create SSBOs
4. All GPU dispatch happens on the worker thread using the compute context
5. On job completion: `wglMakeCurrent(NULL, NULL)` to release

**Shared resource implications:** The voxel SSBO is created on the compute context. Since contexts share texture/buffer namespaces via `wglShareLists`, both contexts can see the same GL objects. But we don't actually need sharing — the compute context's resources are private to the arrange job.

**Why this is cleanest:**
- Zero contention with the render thread — the UI never blocks on collision checks
- No marshaling, no message queues, no synchronization beyond GL's internal fence
- The worker thread owns the context for its entire lifetime, which is the GL threading model working as designed

**Risks:**
- Some drivers handle shared contexts poorly (Intel HD 4000-era, mostly dead by now)
- Need to handle the case where compute shader compilation fails (fall back to CPU)
- wxWidgets' GL context sharing on Windows requires the hidden canvas to remain alive

**Implementation sketch (C++ side):**

```cpp
// In OpenGLManager or a new SnuggleGPU helper:
class SnuggleGPUContext {
    wxGLCanvas*  m_hidden_canvas = nullptr;
    wxGLContext* m_compute_ctx   = nullptr;
    bool         m_initialized   = false;

public:
    // Called on main thread during app init
    bool create_shared(wxWindow* parent, wxGLContext* share_with);

    // Called on worker thread
    bool make_current();
    void release();
    bool has_compute_support() const;  // GL 4.3+ check
};
```

### Option B: Main-Thread Dispatch with Synchronization

Keep all GL calls on the main thread. The worker thread posts a "please evaluate these placements" request, blocks until the main thread dispatches and reads back results.

**Why this is worse:**
- The main thread is the UI thread. Blocking it for GPU work (even ~0.5ms) causes jank during arrange
- Requires a producer-consumer queue between worker and main threads
- The existing `BoostThreadWorker` has a `MainThreadCallData` mechanism for exactly this — but it's designed for rare, fast operations like progress updates, not 50-100 round trips per second

Only consider this if shared GL contexts prove unreliable on target hardware.

### Option C: Vulkan Compute (Future)

Vulkan has no thread-affinity constraint — command buffers can be recorded on any thread and submitted to any queue. OrcaSlicer doesn't currently use Vulkan, so this would mean adding a Vulkan dependency for one feature. Not worth it now, but worth noting as the long-term answer if GPU collision becomes a core feature.

### Recommendation

**Option A.** Create a shared compute context. Fall back to CPU if:
- GL version < 4.3 (no compute shaders)
- Context creation fails
- Shader compilation fails

The CPU path remains the default. GPU is opt-in, detected at runtime.

---

## Integration Plan

### Where to Hook

The GPU evaluator replaces the inner loop of `SnuggleNester::evaluate()`. The nester itself doesn't change — we add a `GpuCollisionEvaluator` class that the nester calls instead of inline collision code.

```
snuggle_nester.hpp  — SnuggleNester::evaluate() calls evaluator interface
gpu_collision.hpp   — GpuCollisionEvaluator (new file, GL compute path)
gpu_collision.cpp   — Shader source, SSBO management, dispatch
snuggle_nester.hpp  — CPU fallback remains inline (current code)
```

The evaluator interface:

```cpp
// In snuggle_nester.hpp or a new header
class CollisionEvaluator {
public:
    virtual ~CollisionEvaluator() = default;

    // Evaluate all individuals in the population.
    // Writes collision_count and oob_count into each Individual.
    virtual void evaluate_batch(
        std::vector<Individual>& population,
        const std::vector<PartInfo>& parts,
        const NesterConfig& cfg) = 0;
};

class CpuCollisionEvaluator : public CollisionEvaluator {
    // Current inline code, extracted into a class
    // Uses the rotation cache from SnuggleNester
};

class GpuCollisionEvaluator : public CollisionEvaluator {
    // GL compute shader path
    // Uploads rotation cache at construction time
    // Each evaluate_batch: upload placements, dispatch, readback
};
```

`SnuggleArrange.cpp` creates the evaluator based on capability detection:

```cpp
std::unique_ptr<CollisionEvaluator> evaluator;
if (gpu_ctx && gpu_ctx->has_compute_support()) {
    evaluator = std::make_unique<GpuCollisionEvaluator>(gpu_ctx, parts, cfg);
} else {
    evaluator = std::make_unique<CpuCollisionEvaluator>(parts, cfg);
}
snuggle::SnuggleNester nester(cfg, evaluator.get());
```

### Initialization Sequence

1. **App startup (main thread):** `SnuggleGPUContext::create_shared()` creates hidden canvas + shared GL context. Stash in `GUI_App` or a singleton.
2. **Arrange job starts (worker thread):** `SnuggleArrange.cpp` checks if GPU context exists, calls `make_current()`, creates `GpuCollisionEvaluator` which compiles shader + uploads voxel data.
3. **Per generation (worker thread):** `evaluate_batch()` uploads placements, dispatches, reads back results.
4. **Arrange job ends (worker thread):** Evaluator destructor deletes GL objects. `release()` unbinds context.

### Fallback Logic

```cpp
bool try_gpu = snuggle_gpu_ctx != nullptr
            && snuggle_gpu_ctx->make_current()
            && snuggle_gpu_ctx->has_compute_support();

if (try_gpu) {
    try {
        evaluator = std::make_unique<GpuCollisionEvaluator>(...);
    } catch (...) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: GPU init failed, falling back to CPU";
        snuggle_gpu_ctx->release();
        try_gpu = false;
    }
}
if (!try_gpu) {
    evaluator = std::make_unique<CpuCollisionEvaluator>(...);
}
```

---

## Fitness Scoring

Fitness scoring (bounding box compactness, height centering) stays on the CPU. It only runs for feasible individuals, involves simple floating-point arithmetic over N placements, and touches no voxel data. Moving it to GPU would save microseconds while adding complexity. Not worth it.

---

## What Could Go Wrong

**Driver bugs with shared contexts.** Some older AMD drivers mishandle `wglShareLists`. Mitigation: if context creation succeeds but first dispatch produces garbage, detect via a known-answer test (upload a simple 2-voxel collision case, verify result = 1) and fall back to CPU.

**SSBO size limits.** GL 4.3 guarantees `GL_MAX_SHADER_STORAGE_BLOCK_SIZE >= 128 MB`. Our 84 MB voxel buffer fits. For 30+ parts, we might exceed this. Mitigation: split into multiple SSBOs (e.g., one per part) or reduce rotation bins.

**Shader compilation time.** First-time compilation of the compute shader can take 100-500ms on some drivers. This is amortized over the entire arrange job (20+ seconds) so it doesn't matter, but we should compile asynchronously or during context setup rather than blocking the first `evaluate_batch`.

**Integer overflow in collision count.** Each pair can produce at most ~16K collisions (overlap region). With 105 pairs, max is ~1.7M. Fits in uint32 with room to spare.

**Voxel data alignment.** The bit-packed format stores 8 voxels per byte. When reading as uint32 words in the shader, byte offsets must be 4-byte aligned. The upload code should pad each grid's data to a 4-byte boundary.

---

## Open Questions

1. **Population size scaling.** At pop_size=512 (current default), we dispatch 512 work groups. Each group needs 128 invocations. That's 65,536 total invocations — well within GPU capacity but may underutilize very large GPUs. Bumping pop_size to 1024 or 2048 is free on the GPU side and may improve GA convergence.

2. **Async readback.** The current design uses synchronous readback (`glGetBufferSubData`), which stalls the worker thread until the GPU finishes. An async path using persistent mapped buffers (`GL_MAP_PERSISTENT_BIT`) with fences would let the CPU start mutation/crossover while the GPU finishes the previous generation's evaluation. This is a second-pass optimization.

3. **2D-only mode.** Most prints are bed-locked with nz=1. In this case, the voxel grids are 2D bitmaps and the inner loop degenerates to 2D overlap. A specialized 2D shader could use `imageLoad` on 2D textures and bitwise AND on full texels. Worth investigating if 2D nesting remains the common case.

4. **Min-gap enforcement.** The current CPU code inflates grids by `min_gap_mm` during voxelization. If we want to change the gap without re-voxelizing, the GPU shader could dilate the overlap check region by gap/voxel_size voxels in each direction. Not needed for v1 — just inflate during voxelization as before.
