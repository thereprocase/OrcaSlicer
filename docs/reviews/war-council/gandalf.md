# War Council: Gandalf the Grey -- Architecture Review

*"I have watched codebases grow from careful seedlings into tangled forests. This one shows the marks of a builder who thinks before reaching for the axe."*

---

## Summary Verdict

The Snuggle nester is architecturally sound for a proof-of-concept that has matured into a working feature. The layer separation is clean: voxelizer, nester, collision evaluator, and bridge each own a distinct responsibility. Several structural decisions are correct and non-obvious (feasibility-first selection, rotation cache, conservative voxelization). The code is not yet PR-ready for upstream OrcaSlicer, but the remaining work is refinement, not redesign.

---

## 1. Header-Only Nester: Should It Be Split?

**Verdict: Yes, split it. But not urgently.**

`snuggle_nester.hpp` is ~1000 lines of implementation in a header. Every translation unit that includes it gets the full GA, all four seed strategies, the compaction pass, fitness scoring, and evaluation. Today only `SnuggleArrange.cpp` includes it, so the cost is zero. But:

- If a second consumer appears (tests, a standalone benchmark, a different bridge), compile times multiply.
- The header already has a forward declaration of `CollisionEvaluator` and depends on `polite_voxelizer.hpp` -- these are implementation details leaking into the interface.
- The public surface is tiny: `SnuggleNester(cfg, evaluator, seed)` and `run(parts, progress)`. Everything else is private. A clean split would put the class declaration + `NesterConfig` + `PartInfo` + `NesterResult` + `Individual` + `Placement` in the `.hpp`, and all method bodies in the `.cpp`.

**Risk of not splitting**: Low today, grows with adoption. Worth doing before upstream PR.

---

## 2. CollisionEvaluator Interface

**Verdict: Clean interface, solid CPU/GPU split. Two issues.**

The abstract interface is well-designed:

```cpp
virtual void upload_grids(parts, rot_cache) = 0;
virtual void evaluate_batch(pop, parts, bed_w, bed_h, margin) = 0;
```

Two methods, clear contract: "I own collision/bounds, you own fitness." The factory pattern (`create_collision_evaluator()`) with graceful GPU-to-CPU fallback is correct.

### Issue 1: GPU evaluator ignores `bed_margin`

```cpp
// TODO: pass bed_margin to shader as uniform for GPU OOB check
(void)bed_margin;
```

The GPU shader checks `pmin_x < 0.0` instead of `pmin_x < bed_margin`. This means GPU and CPU evaluators produce different results for the same population. The GA could converge on a GPU-feasible solution that the CPU post-validation (greedy accept in `SnuggleArrange.cpp`) rejects. This is a correctness bug, not just a TODO.

### Issue 2: Rotation cache coupling

Both `CpuCollisionEvaluator` and `GpuCollisionEvaluator` hardcode `ROT_CACHE_BINS = 360` independently, and the nester also has `ROT_CACHE_BINS = 360`. If anyone changes one without the others, the angle-to-bin mapping silently breaks. This constant should live in exactly one place -- either in the `CollisionEvaluator` base class or as a shared constant in a common header.

### Minor: `CpuCollisionEvaluator` stores a raw pointer to `rot_cache_`

```cpp
const std::vector<std::vector<VoxelGrid>>* rot_cache_ = nullptr;
```

This is a non-owning pointer to data owned by `SnuggleNester`. The lifetime contract is correct (nester outlives evaluator during `run()`), but fragile. A comment documenting the lifetime requirement would help future maintainers. Alternatively, `upload_grids` could copy what it needs.

---

## 3. Bridge Layer (SnuggleArrange.cpp)

**Verdict: Doing the right amount of work. Voxelization belongs here.**

The bridge handles:
1. Bed geometry extraction
2. Single-part fast path
3. Part-count and sequential-print guards
4. Mesh extraction + transform + decimation
5. Voxelization
6. Nester configuration
7. Result writeback with greedy fallback

This is a lot of code (~380 lines), but each step is necessary glue between OrcaSlicer's data model and the nester's clean input format. Moving voxelization into a separate translation unit would reduce the bridge's line count, but the voxelizer already lives in its own header. The bridge calls it -- that's the right relationship.

### Structural concern: instance matching by name

```cpp
if (oi.obj && oi.obj->name == items[i].name) {
    obj = oi.obj;
    inst = oi.inst;
    oi.obj = nullptr;  // mark used
    break;
}
```

This walks a flat list of `(object, instance)` pairs and matches by name, consuming matches. If two objects share a name (user renames parts, or imports duplicates), this silently mismatches meshes to silhouettes. The comment acknowledges this is matching `prepare_all()` enumeration order -- but the match-by-name strategy is fragile. A positional match (same index) or a hash of the mesh would be more robust.

### The greedy-accept fallback is O(n^2) in collision checks

When the GA result is infeasible, the bridge does pairwise `VoxelGrid::collision_count()` against all accepted parts for each candidate. For the 200-part limit this is 20K collision checks, each potentially iterating a 3D overlap volume. Acceptable for the current part-count ceiling, but the ceiling should not be raised without addressing this.

---

## 4. Abstraction Leaks

### 4a. `snuggle::Vec3f` vs Slic3r's `Vec3f`

The voxelizer defines its own `Vec3f` in the `snuggle` namespace. The bridge manually converts between Slic3r's `Vec3d`/`Vec3crd` types and `snuggle::Vec3f`. This is deliberate isolation (the snuggle namespace has no Slic3r dependencies), which is architecturally correct for a self-contained module. But it means the bridge does a lot of format conversion in the voxelization loop. Not a leak per se, but a friction point.

### 4b. `polite_voxelizer.hpp` defines types used by `snuggle_nester.hpp`

`VoxelGrid`, `Vec3f`, `VoxError` -- these are defined in the voxelizer header but are core types for the entire snuggle system. The voxelizer is really two things in one file: a type library and a voxelization algorithm. Splitting the types into `snuggle_types.hpp` would make the dependency graph cleaner.

### 4c. `#ifdef SLIC3R_GUI` in the evaluator

The GPU evaluator is `#ifdef`'d out for headless builds. The factory handles this gracefully. But the nester's `run()` method has this:

```cpp
#ifdef SLIC3R_GUI
    if (auto* gpu = dynamic_cast<snuggle::GpuCollisionEvaluator*>(evaluator.get())) {
```

This is in `SnuggleArrange.cpp`, which is the right place for it. The nester itself has no `#ifdef` guards -- also correct. The evaluator interface fully hides the GPU/CPU distinction from the algorithm. Well done.

---

## 5. PR Readiness for Upstream OrcaSlicer

**Not yet. Here's what needs attention:**

### Must-fix

1. **GPU `bed_margin` bug** -- CPU and GPU evaluators must produce identical results for the same inputs. The TODO in `gpu_collision.cpp:477` is a correctness issue.

2. **`ROT_CACHE_BINS` triplicated** -- Single source of truth needed.

3. **`pi` constant scattered as `3.14159265f`** -- Use `M_PI` or define `constexpr float TWO_PI` once. Seven separate literals in the nester alone.

### Should-fix

4. **Split `snuggle_nester.hpp` into `.hpp`/`.cpp`** -- Standard practice for upstream contributions. Header-only is a red flag for reviewers.

5. **Instance matching by name** -- Document the fragility or switch to positional matching.

6. **`polite_yield()` in collision inner loop** -- The `collision_count()` method yields every 10K iterations. During GA evaluation with hundreds of individuals, each checking N*(N-1)/2 pairs, this adds up to thousands of yield calls per generation. Profile whether this actually helps or just adds syscall overhead. The polite_yield in the *generation* loop is sufficient for responsiveness.

7. **No unit tests** -- The nester, voxelizer, and collision evaluator are all testable in isolation. A PR without tests will be rejected by any serious project.

### Nice-to-have

8. **`VoxelGrid::count_solid()` uses manual popcount** -- `__builtin_popcount` / `_mm_popcnt_u32` / `std::popcount` (C++20) would be faster and shorter.

9. **The compaction phase rebuilds `rotated_copy()` for every part in every sweep** -- The rotation doesn't change during compaction (only position does), so the rotated grids could be built once and reused across sweeps. The micro-rotation path does need fresh copies, but only for the part being tested.

10. **GPU shader uses `local_size_x = 1`** -- One invocation per individual, serial pair iteration within each. This underutilizes GPU parallelism. The comment acknowledges this ("POC version"). For upstream, either parallelize across pairs within each work group, or document why the current approach is sufficient.

---

## 6. What's Done Right

I want to be clear about what this architecture gets right, because these are the decisions a maintainer will thank you for in two years:

- **Clean namespace isolation.** The entire snuggle system lives in `namespace snuggle` with no Slic3r includes. The bridge is the only file that touches both worlds. This means the nester could be extracted into an independent library with zero refactoring.

- **Feasibility-first selection (Deb 2000).** Correct implementation. Infeasible individuals always lose to feasible ones. This is the right approach for constrained optimization and avoids the common pitfall of penalty-weighted fitness where the GA learns to minimize penalty rather than find feasible solutions.

- **Conservative voxelization.** Surface-only, outward-rounding. Parts are "fatter" in voxel space than reality, which means collision detection is conservative (no false negatives). The comment explicitly explains this tradeoff. Good.

- **Graceful degradation chain.** Single part -> center. Too many parts -> standard arranger. Sequential print -> standard arranger. Voxelization failure -> off-plate. GA infeasible -> greedy accept + overflow. Every failure mode has a defined behavior.

- **The "politeness" contract.** Memory caps, yield points, bounded loops, progress callbacks, timeouts. This is the work of someone who has shipped code that runs on other people's machines. The contract is documented in the header comments and enforced in the implementation.

---

## Architecture Diagram (as I understand it)

```
ArrangeJob.cpp (GUI thread)
    |
    v
SnuggleArrange.cpp (bridge)
    |-- reads Model objects, extracts meshes
    |-- decimates high-poly meshes
    |-- calls polite_voxelizer to build VoxelGrids
    |-- configures NesterConfig from ArrangeParams
    |-- creates CollisionEvaluator via factory
    |-- instantiates SnuggleNester, calls run()
    |-- writes Placement results back to ArrangePolygons
    |-- handles overflow fallback
    |
    v
snuggle_nester.hpp (GA core)
    |-- builds rotation cache (360 bins per part)
    |-- uploads to CollisionEvaluator
    |-- evolves population: select, crossover, mutate, evaluate
    |-- post-GA compaction
    |-- returns NesterResult
    |
    +-- CollisionEvaluator (interface)
        |-- CpuCollisionEvaluator (always available)
        |-- GpuCollisionEvaluator (OpenGL 4.3, Windows only, SLIC3R_GUI)
```

The layers are clean. Data flows downward. The evaluator interface is the only horizontal seam, and it's well-defined.

---

*"A wizard is never late, nor is he early. He delivers precisely when the architecture is sound. This one is close."*
