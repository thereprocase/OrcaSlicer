# Legolas Scout Report — Snuggle Performance Audit

**Scope:** `snuggle_nester.hpp`, `polite_voxelizer.hpp`, `gpu_collision.cpp`, `SnuggleArrange.cpp`
**Date:** 2026-03-31
**Reviewer:** Legolas (Scout)

---

## Summary

The architecture is solid and the rotation cache is the right call. The big problems are in
`compact_toward_center`, which is structurally an O(n² * V) disaster per part per sweep, and
in the GPU path, which uploads the placement buffer and immediately stalls on `glGetBufferSubData`
with no pipeline overlap. Several smaller issues compound on the hot paths.

---

## CRITICAL — O(n²) or worse

### 1. `compact_toward_center`: O(sweeps × n² × V) with repeated `rotated_copy` inside the inner loop

**File:** `snuggle_nester.hpp`, lines ~376–490

This is the worst offender in the codebase. The structure is:

```
for sweep in 0..compact_max_sweeps:         // up to 20
  for idx in sorted_parts:                   // n parts
    rot_i = parts[idx].grid.rotated_copy(…)  // allocates + stamps full 3D grid
    for j in 0..n:                           // n-1 others
      rot_others_storage.push_back(parts[j].grid.rotated_copy(…))  // allocates again

    for bs in 0..12 (binary search):
      for j in 0..n:
        VoxelGrid::collision_count(…)        // walks overlap volume V

    for da in test_angles (4 angles):
      for j in 0..n:
        rot_j = parts[j].grid.rotated_copy(…)   // ANOTHER full copy per angle per part
        VoxelGrid::collision_count(…)
      if ok:
        for bs2 in 0..8:
          rot2 = parts[idx].grid.rotated_copy(test_rot)  // ANOTHER copy inside inner loop
          for j in 0..n:
            rj2 = parts[j].grid.rotated_copy(…)          // AND ANOTHER
            VoxelGrid::collision_count(…)
```

The outermost structure is O(sweeps × n), the collision checks are O(n × V), and the rotation
copies are O(grid_volume). On each part's turn, it:

- Builds `rot_others` by allocating and copying the full voxel grid of every other part. That
  vector is rebuilt from scratch for every part in every sweep. With 20 parts and 20 sweeps,
  that's 400 full grid copies for the `rot_others` block alone.

- Then, inside the micro-rotation loop, builds *more* rotated copies of every other part again
  (lines ~441–464). The `rot_j` and `rj2` copies are redundant — the rotation for part `j` does
  not change during idx's move, so these copies are identical to the ones already computed two
  dozen lines above.

- Inside the `bs2` binary-search loop (8 iterations), calls
  `parts[idx].grid.rotated_copy(test_rot)` at line ~457. Same angle, same grid, 8 allocations
  in a tight loop. The result is identical every iteration.

**Impact:** With 20 parts at 50×50×50 voxels (125K voxels per grid, 16 KB bit-packed),
rotated_copy is not free — it iterates nx×ny×nz=125K voxels per call. The inner bs2 loop
alone does 8 × n full copies per angle per part per sweep.

**Fix:** Pre-compute and cache `rot_others` once per sweep. Hoist `rotated_copy(test_rot)` out
of the `bs2` loop (the angle does not change inside bs2). The existing `rot_others_storage`
vector is a step in the right direction but is rebuilt instead of reused.

---

### 2. `compute_fitness`: O(n²) proximity loop on every feasible individual every generation

**File:** `snuggle_nester.hpp`, lines ~920–932

```cpp
for (size_t i = 0; i < n; i++)
    for (size_t j = i + 1; j < n; j++)
        proximity_score += …;
```

This runs once per *feasible* individual in the population, every generation. With pop=512 and
n=20 parts, a fully feasible population costs 512 × 190 pair distance computations = ~97K
sqrt() calls per generation. The loop itself is fine at small n, but `compute_fitness` is also
called redundantly during `evaluate()` on the CPU fallback path and again on the external
evaluator path, which both call it separately. There is no caching of part center positions
across calls — the center x/y for each part is recomputed twice (once for the centroid, once
for the avg_dist loop, and once more in the proximity loop).

**Fix:** Minor: compute part centers once, store in a local array, reuse for all three passes.

---

### 3. `seed_bottom_left`: O(n² × bed²/gap²) worst case

**File:** `snuggle_nester.hpp`, lines ~702–719

The outer loop is over parts. The inner double loop is over the bed grid at `min_gap_mm`
resolution. For each candidate position the code scans all previously placed rects. With a
5mm gap on a 256mm bed, the grid is ~50×50=2500 candidate positions per part. With 20 parts,
worst case is 20 × 2500 × 19 = ~950K rect tests. This is only called once at startup for the
`seed_bottom_left` individual, so it is not a per-generation cost — but it is an O(n³) seeder
when placements fail to find a spot early.

---

## HIGH — Per-frame allocations in hot loops

### 4. `evaluate_batch` (CPU): `std::vector<const VoxelGrid*> rotated(n)` allocated per individual

**File:** `gpu_collision.cpp`, lines ~55–57

Inside `CpuCollisionEvaluator::evaluate_batch`, a fresh `std::vector<const VoxelGrid*>` of
size n is allocated on the heap for every individual in the population, every generation. With
pop=512, that is 512 heap allocations per generation × up to 100 generations = 51,200
small-vector allocations. Each is freed immediately after the individual is processed. These
will hammer the allocator.

**Fix:** Allocate `rotated` once outside the individual loop and reuse it. It is the same size
for every individual.

### 5. `compact_toward_center`: `rot_others_storage` reserve + push_back per part per sweep

**File:** `snuggle_nester.hpp`, lines ~380–385

`rot_others_storage.reserve(n)` and `rot_others_storage.push_back(…)` are inside the outer
part loop. The vector is declared inside the loop body, so it is constructed and destroyed every
iteration. With n=20 parts × 20 sweeps = 400 construction/destruction cycles, each holding up
to 19 full VoxelGrid objects with their heap-allocated bit vectors.

**Fix:** Declare `rot_others_storage` outside the `for (size_t idx : order)` loop and call
`.clear()` + `.resize(n)` at the top of each iteration. Or better: cache rotated grids for
the entire sweep before the per-part loop begins, since part rotations do not change during
a sweep.

### 6. `SnuggleArrange.cpp`: intermediate `std::vector<float> verts` and `std::vector<uint32_t> indices` per part

**File:** `SnuggleArrange.cpp`, lines ~169–181

For every part, two transient vectors are allocated to reformat mesh data before passing it to
`voxelize_indexed_mesh`. The vertices and indices already exist in `its.vertices` and
`its.indices` in Eigen format. The reformat is a straight copy with coordinate extraction.
This is a one-time cost per part, not per-generation, so it is not a hot-path issue — but it
doubles peak memory usage for the vertex buffer.

**Fix:** Teach `voxelize_indexed_mesh` to accept Eigen `Vector3f` vertices directly, or write
a thin adapter that reads from `its` without copying.

---

## HIGH — GPU sync stall

### 7. `GpuCollisionEvaluator::evaluate_batch`: placement upload blocks, then immediate readback stalls

**File:** `gpu_collision.cpp`, lines ~496–529

The sequence each generation is:

```
glBufferData(placement_ssbo_)     // Upload placements — implicit pipeline stall if GPU is busy
glBufferData(results_ssbo_)       // Allocate results buffer — another potential stall
glDispatchCompute(pop_size, 1, 1) // Dispatch
glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)
glGetBufferSubData(results_ssbo_) // FULL CPU-GPU SYNC — blocks until GPU finishes
```

`glGetBufferSubData` is a synchronous readback. It inserts a pipeline bubble: the CPU sits idle
waiting for the GPU to drain. On a modern discrete GPU the round-trip latency is ~0.5–2ms per
generation, which dwarfs the actual shader execution time for small populations.

Two additional issues:

- `glBufferData` on the placement SSBO uses `GL_DYNAMIC_DRAW` and re-orphans the buffer every
  call. This triggers implicit synchronization if the previous frame's dispatch is still in
  flight. The driver cannot return until it knows the old buffer is no longer in use.
  Double-buffering the placement SSBO would eliminate this stall.

- `glUniform1ui(glGetUniformLocation(program_, "u_n_parts"), …)` at line ~515 performs a
  uniform location query every generation. `glGetUniformLocation` is a hash-map lookup inside
  the driver. The locations are constant — they should be cached at `compile_shader()` time.

**Fix:**
  - Cache uniform locations after linking.
  - Double-buffer the placement SSBO (ping-pong between two buffers, use fence sync to check
    the previous one is safe to overwrite).
  - Replace `glGetBufferSubData` with `glMapBufferRange(GL_MAP_READ_BIT)` after a
    `glClientWaitSync` on a fence inserted after dispatch. This gives the driver more
    scheduling freedom.

---

### 8. `GpuCollisionEvaluator::evaluate_batch`: `bed_margin` parameter silently ignored

**File:** `gpu_collision.cpp`, line ~477

```cpp
// TODO: pass bed_margin to shader as uniform for GPU OOB check
(void)bed_margin;
```

The GPU OOB check uses a hardcoded comparison against `0.0` and `u_bed_w`/`u_bed_h` (shader
lines ~174–177) rather than the configured margin. This means the GPU evaluator's OOB count
diverges from the CPU evaluator's — the GPU reports fewer OOB violations. If the GPU backend
is selected, the feasibility-first selection will accept individuals that the CPU path would
have rejected, producing layouts that clip the bed margin.

This is a correctness bug that masquerades as a performance win.

---

## MEDIUM — Cache-unfriendly access patterns

### 9. `collision_count`: floating-point triple-loop over world space with per-iteration integer converts

**File:** `polite_voxelizer.hpp`, lines ~308–330

The overlap iteration uses floating-point world-space coordinates (`wx`, `wy`, `wz`) as loop
variables and converts them to grid indices via `world_to_grid` on every iteration. Each
`world_to_grid` call does three `std::floor` and three subtractions. For a 50×50×50 voxel
overlap, that is 125K float-to-int conversions per pair.

More important: the voxel access pattern on grid A is sequential in X (matching the `x + y*nx
+ z*nx*ny` linear layout), but the access pattern on grid B is shifted by the offset between
the two grids' origins. If the two parts are misaligned by an odd offset, B's access pattern
is a non-unit stride in the bit array, producing scattered byte reads with no vectorization.

The `polite_yield()` call every 10,000 iterations (line ~325) inserts a `SwitchToThread()`
inside the collision inner loop. This is a non-trivial syscall. For a 125K-iteration overlap
that's 12 yields per pair per evaluation — unnecessary overhead when the collision loop is
already completing in < 1ms.

**Fix:** Convert the loop bounds to integer grid coordinates once before the loop, iterate in
integer grid space, and access both grids with integer index arithmetic. Eliminate the yield
from the collision inner loop (the outer evaluate loop already yields per-generation).

### 10. `rotated_copy`: forward scatter into destination, not gather from source

**File:** `polite_voxelizer.hpp`, lines ~235–258

`rotated_copy` iterates source voxels in XYZ order (good for source access) and scatters each
set bit into a computed destination position. Destination writes are non-sequential — the
rotated index `(dx, dy, z)` does not follow the destination grid's linear order. This means
each `rot.set()` is a random write into the destination bit array, breaking cache locality on
the destination side.

For the rotation cache build (360 bins × n parts), this is a one-time cost per nesting run,
so the impact is bounded. For the compaction loop where `rotated_copy` is called repeatedly
without the cache (see issue #1 and #4), this degrades every compaction step.

### 11. `voxelize_mesh`: inner triangle loop recomputes `grid_to_world_center` every voxel

**File:** `polite_voxelizer.hpp`, lines ~492–505

```cpp
for (int iz …) for (int iy …) for (int ix …) {
    Vec3f center = out_grid.grid_to_world_center(ix, iy, iz);
    Vec3f v0 = tri.v0 - center;
    …
    triangle_aabb_overlap(v0, v1, v2, half);
}
```

`grid_to_world_center` computes `origin + (idx + 0.5) * voxel_size` for all three axes every
iteration. The Z component only changes in the outermost loop, and the Y component only changes
in the middle loop. Both are being redundantly recomputed in the innermost loop.

**Fix:** Hoist the Z and Y components: compute `center_z` in the iz-loop, `center_y` in the
iy-loop, and only compute `center_x` in the ix-loop.

---

## MEDIUM — Shader architecture

### 12. GPU shader: one invocation = one individual, no work-group parallelism

**File:** `gpu_collision.cpp`, shader at line ~116

```glsl
layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;
```

Each compute invocation evaluates one entire individual serially — all n parts, all n(n-1)/2
pairs. A GPU thread handles this with no intra-individual parallelism. With local_size=1 and
pop=512, the GPU dispatches 512 single-threaded invocations. Most GPUs execute warps of 32
threads; single-thread work groups waste 31/32 of each warp's capacity.

The shader comment at line ~117 acknowledges this: "A future version would use work-group
parallelism across pairs." This is the right direction. Pairs can be assigned to work-group
lanes; a final parallel reduction collapses the collision counts.

This is a correctness-neutral performance issue — it means the GPU path is roughly 32× slower
than it could be on the pairwise evaluation itself.

### 13. GPU shader: redundant `metas` lookup for part `i` on every inner-loop iteration

**File:** `gpu_collision.cpp`, shader lines ~164–168 and ~180–185

`meta_idx_i`, `mi`, `pmin_x`, `pmin_y`, `pmax_x`, `pmax_y`, and the OOB check for part `i`
are all computed correctly *before* the inner `j` loop. However, `ai_min_x` and `ai_min_y`
are recomputed inside the j-loop from `pi.x + mi.origin_x` (lines ~188–189). These are
identical to `pmin_x` and `pmin_y` computed five lines earlier. Minor redundancy.

---

## LOW — Miscellaneous

### 14. `count_solid`: manual popcount instead of `__builtin_popcount`

**File:** `polite_voxelizer.hpp`, lines ~268–274

```cpp
while (b) { count += b & 1; b >>= 1; }
```

This is a bit-serial loop over each byte. `__builtin_popcount` (or `std::popcount` in C++20)
compiles to a single instruction on x86. `count_solid` is not on the hot path but is called
after voxelization to check for empty parts (SnuggleArrange.cpp line ~188), so the impact is
minor.

### 15. `build_rotation_cache` in `lock_rotation` mode: 360 copies of the same grid

**File:** `snuggle_nester.hpp`, lines ~534–537

```cpp
rot_cache_[i].resize(ROT_CACHE_BINS, rotated); // fill all bins with same copy
```

`std::vector::resize(N, value)` with a `VoxelGrid` value copy-constructs N=360 copies. Each
copy duplicates the bit vector heap allocation. For a 50×50×50 grid, that is 360 × 16KB =
5.6MB per part, all holding identical data. This is intentional (GPU evaluator requires all
bins to be populated), but 360× the memory and 360 heap allocations per locked part is
expensive for something that could be a pointer-to-one-copy with an indirection.

### 16. `snap_rotation`: while-loop normalization instead of fmod

**File:** `snuggle_nester.hpp`, lines ~567–570

```cpp
while (delta < 0) delta += TWO_PI;
while (delta >= TWO_PI) delta -= TWO_PI;
```

`std::fmod(delta, TWO_PI)` with a clamp for the negative case is O(1). The while loops are
also O(1) in practice (angles are rarely far out of range), but this is an unnecessary
pattern that would misbehave if a mutation accumulated many full rotations.

---

## Complexity Summary Table

| Location | Operation | Complexity | Hot? |
|---|---|---|---|
| `compact_toward_center` inner micro-rotation loop | rotated_copy × n per angle per bs2 iteration | O(sweeps × n² × V × 12 × 4 × 8) | Yes |
| `compact_toward_center` rot_others build | rotated_copy × n-1 per part per sweep | O(sweeps × n² × V) | Yes |
| `evaluate_batch` CPU | pairwise collision | O(pop × n² × V) | Yes (per gen) |
| `compute_fitness` proximity | pairwise distance | O(pop × n²) | Yes (per gen, feasible) |
| `voxelize_mesh` rasterization | triangles × voxels-per-tri | O(T × (extent/vs)³) | Once |
| `build_rotation_cache` | 360 × n rotated_copy | O(360 × n × V) | Once |
| `seed_bottom_left` | O(n² × bed²/gap²) | O(n³) worst | Once |

V = voxel grid volume (nx × ny × nz)

---

## Priority Hits

1. **`compact_toward_center` micro-rotation loop** (issues #1, #5): most expensive single block.
   Cache `rot_others` once per sweep before the per-part loop. Hoist `rotated_copy(test_rot)`
   out of `bs2`. Together these eliminate the majority of the allocation pressure.

2. **GPU sync stall + uniform location query** (issue #7): double-buffer placement SSBO, cache
   uniform locations. Should halve per-generation GPU overhead.

3. **GPU `bed_margin` not passed to shader** (issue #8): correctness bug. Fix before trusting
   GPU results.

4. **`evaluate_batch` CPU rotated vector allocation** (issue #4): one-line fix, 512
   allocations/gen eliminated.

5. **`collision_count` float loop + yield** (issue #9): convert to integer loop, remove inner
   yield. Measurable improvement on CPU path.
