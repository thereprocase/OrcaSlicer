# Build System Audit — Snuggle Feature Branch
**Reviewer:** Gimli the Craftsman
**Date:** 2026-03-31
**Scope:** CMakeLists.txt, gpu_collision.cpp/hpp, polite_voxelizer.hpp, snuggle_nester.hpp, SnuggleArrange.cpp
**Verdict:** Stone is mostly sound. Five cracks found — two are load-bearing.

---

## 1. CMakeLists.txt

### File registration
All five new files are listed in `lisbslic3r_sources` (note: this variable name has a typo — `lisbslic3r` not `libslic3r` — but this is a pre-existing issue, not new):

```
Arrange/SnuggleArrange.cpp
Arrange/SnuggleArrange.hpp
Arrange/gpu_collision.cpp
Arrange/gpu_collision.hpp
Arrange/polite_voxelizer.hpp
Arrange/snuggle_nester.hpp
```

Header-only files (`polite_voxelizer.hpp`, `snuggle_nester.hpp`) are listed alongside `.cpp` files. This is cosmetic — CMake ignores headers in a `STATIC` target's source list for compilation purposes, but it is conventional and harmless. Source group indexing benefits from the listing.

### Conditional GPU linking

```cmake
if(SLIC3R_GUI)
    target_link_libraries(libslic3r PRIVATE GLEW::GLEW OpenGL::GL)
endif()
```

This is correct. `GLEW::GLEW` and `OpenGL::GL` are only linked when `SLIC3R_GUI` is set, which mirrors the `#ifdef SLIC3R_GUI` guard in `gpu_collision.cpp`. The root `CMakeLists.txt` confirms `find_package(OpenGL REQUIRED)` and `find_package(GLEW REQUIRED)` are called unconditionally at project scope (lines 700, 714), so the targets exist whenever needed. The `PRIVATE` linkage is correct — downstream consumers of `libslic3r` do not need GL headers.

### GLEW find_package scope
`find_package(GLEW)` and `find_package(OpenGL)` live in the root `CMakeLists.txt`, not in `src/libslic3r/CMakeLists.txt`. This is the existing project convention. The `GLEW::GLEW` and `OpenGL::GL` targets are therefore available when the libslic3r CMakeLists runs. No issue.

### gpu_collision.cpp is NOT conditionally compiled
**FINDING — MINOR.** `gpu_collision.cpp` is always compiled into `libslic3r` regardless of `SLIC3R_GUI`. In headless builds, `CpuCollisionEvaluator` and `create_collision_evaluator()` compile fine because they carry no GL dependency. The GPU class body is inside `#ifdef SLIC3R_GUI` in both `.hpp` and `.cpp`. So the `.cpp` compiles cleanly in both build modes. This is intentional and correct, but a future engineer reading the CMakeLists might expect to see something like:

```cmake
if(SLIC3R_GUI)
    list(APPEND lisbslic3r_sources Arrange/gpu_collision.cpp)
endif()
```

The current approach works but is less explicit about intent. Not a defect; a documentation gap.

---

## 2. gpu_collision.cpp — `#ifdef SLIC3R_GUI` guards

**Structure:**
- Lines 13–21: `NOMINMAX`, `windows.h`, `GL/glew.h` — all inside `#ifdef SLIC3R_GUI`. Clean.
- Lines 92–582: `GpuCollisionEvaluator` class body — all inside `#ifdef SLIC3R_GUI`. Clean.
- Lines 588–601: `create_collision_evaluator()` factory — outside the guard, but references `GpuCollisionEvaluator` only inside its own `#ifdef SLIC3R_GUI` block. Clean.

**FINDING — DEFECT (LOAD-BEARING).** The `init_context()` method only implements the `_WIN32` path. The non-Windows branch returns `false`, which causes graceful CPU fallback. However, the Win32 context creation uses `wglCreateContext` (legacy OpenGL 1.x context), **not** `wglCreateContextAttribsARB`. This means the created context may be a compatibility profile context, not a core profile. OpenGL 4.3 compute shaders are available in both profiles, so this likely works in practice — but `wglCreateContextAttribsARB` with `WGL_CONTEXT_MAJOR_VERSION_ARB = 4`, `WGL_CONTEXT_MINOR_VERSION_ARB = 3` would be the correct approach. The GL version check at lines 321–325 guards against the worst case (runtime graceful fallback), but a driver that provides `wglCreateContext` for 4.3 in legacy-compat mode is an implementation detail, not guaranteed behavior.

**FINDING — MINOR.** The window created by `init_context()` (`"SnuggleGPUCollision"`) is a hidden 1x1 pixel window. `DestroyWindow` is called in `cleanup()`, but `UnregisterClassA` is never called. In long-running processes that create/destroy `GpuCollisionEvaluator` instances multiple times, the `RegisterClassA` call will return `ERROR_CLASS_ALREADY_EXISTS` on the second call. The code handles this gracefully (the `if (err != ERROR_CLASS_ALREADY_EXISTS)` check), so it does not fail — but the window class leaks. Not a crash risk, cosmetic.

**Headless build check:** The `CpuCollisionEvaluator` and factory function compile without any GL or Windows headers. Headless builds will work. Verified.

---

## 3. polite_voxelizer.hpp — Platform guards

```cpp
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sched.h>
#endif
```

This is correctly guarded. `windows.h` only on Win32, `sched.h` only on non-Win32.

`NOMINMAX` guard is present and correct. Order is right: `NOMINMAX` defined before `windows.h` inclusion.

`polite_yield()` correctly uses `SwitchToThread()` on Win32 and `sched_yield()` on POSIX. Both are appropriate — they yield to any ready thread without sleeping.

No issues found in this file.

---

## 4. snuggle_nester.hpp — Include chain and forward declarations

**Include chain:**
```
snuggle_nester.hpp
  -> polite_voxelizer.hpp (direct include, line 21)
```

`gpu_collision.hpp` is **not** included from `snuggle_nester.hpp`. Instead, the header uses a forward declaration for `CollisionEvaluator`:

```cpp
// Forward declaration — full type needed by callers that invoke
// evaluator methods. SnuggleArrange.cpp includes gpu_collision.hpp
// before this header to satisfy the dependency.
class CollisionEvaluator;
```

This is architecturally clean — `snuggle_nester.hpp` stores a `CollisionEvaluator*` (non-owning pointer), so the forward declaration is sufficient. The full definition is not needed in the header.

**FINDING — DEFECT (LOAD-BEARING).** The comment says "SnuggleArrange.cpp includes gpu_collision.hpp before this header." The actual include order in `SnuggleArrange.cpp` is:

```cpp
#include "SnuggleArrange.hpp"      // line 6
#include "polite_voxelizer.hpp"    // line 7
#include "gpu_collision.hpp"       // line 8
#include "snuggle_nester.hpp"      // line 9
```

`gpu_collision.hpp` is included at line 8, before `snuggle_nester.hpp` at line 9. **This works correctly.** However `gpu_collision.hpp` includes `polite_voxelizer.hpp`, and `snuggle_nester.hpp` also includes `polite_voxelizer.hpp`. Both use `#pragma once`, so double-inclusion is safe. No circular dependency exists.

But: `gpu_collision.hpp` forward-declares `PartInfo` and `Individual` (lines 19–21), while `snuggle_nester.hpp` defines them. When `gpu_collision.hpp` is included before `snuggle_nester.hpp` (as it is in `SnuggleArrange.cpp`), the forward declarations resolve when `snuggle_nester.hpp` is processed. This is fine in `SnuggleArrange.cpp`. However, any other translation unit that includes only `gpu_collision.hpp` (without `snuggle_nester.hpp`) will have unresolved forward declarations for `PartInfo` and `Individual`. Since `gpu_collision.hpp`'s virtual methods take these types by reference, **any file that calls `upload_grids()` or `evaluate_batch()` directly must also include `snuggle_nester.hpp`**. This is an implicit contract not enforced by the headers themselves. A future test file or integration point will likely break here.

**ROT_CACHE_BINS vs ROT_BINS:** Two different constant names for the same value (360):

- `snuggle_nester.hpp` (line 521): `static constexpr int ROT_CACHE_BINS = 360;`
- `gpu_collision.hpp` (line 57): `static constexpr int ROT_CACHE_BINS = 360;` (in `CpuCollisionEvaluator`)
- `gpu_collision.hpp` (line 93): `static constexpr int ROT_BINS = 360;` (in `GpuCollisionEvaluator`)

All three are 360. They are independent class-scoped constants, so they cannot fall out of sync via the preprocessor, but they can fall out of sync via human edit. If someone changes `ROT_CACHE_BINS` in the nester and forgets `ROT_BINS` in the GPU evaluator, the rotation bins will mismatch and lookups will be wrong. This should be a single shared constant.

---

## 5. SnuggleArrange.cpp — Include order and transitive dependencies

**Include order:**
```cpp
#include "SnuggleArrange.hpp"        // -> Arrange.hpp -> ExPolygon, PrintConfig, Print
                                     // -> Model.hpp -> Geometry.hpp (Transformation)
#include "polite_voxelizer.hpp"
#include "gpu_collision.hpp"
#include "snuggle_nester.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/QuadricEdgeCollapse.hpp"
#include <boost/log/trivial.hpp>
```

`SnuggleArrange.hpp` includes `libslic3r/Model.hpp`, and `Model.hpp` includes `Geometry.hpp` (confirmed at Model.hpp line 6). `Geometry.hpp` brings in `Geometry::Transformation`, which is used directly in `SnuggleArrange.cpp` at line 144. This transitive path is intact. No missing include for `Geometry.hpp`.

`unscale_()` is a macro defined in `libslic3r.h` (line 94). `libslic3r.h` is transitively pulled in through the `Model.hpp` → `libslic3r.h` chain. Verified present.

`scaled()` is a function defined in `Point.hpp`. Also transitively available via the `Model.hpp` include chain.

`get_extents()` on `Polygon` is from `BoundingBox.hpp`, which is explicitly included. Clean.

`its_quadric_edge_collapse()` is from `QuadricEdgeCollapse.hpp`, explicitly included. Clean.

No missing includes. No circular dependencies in the include graph.

---

## Summary Table

| Item | File | Severity | Description |
|------|------|----------|-------------|
| Legacy GL context | gpu_collision.cpp | **DEFECT** | `wglCreateContext` creates legacy-compat context, not explicit 4.3 core. Works on most drivers but not guaranteed. |
| `PartInfo`/`Individual` forward-decl contract | gpu_collision.hpp | **DEFECT** | Files including `gpu_collision.hpp` without `snuggle_nester.hpp` will see unresolved forward decls. No enforced include. |
| Dual ROT_BINS constants | gpu_collision.hpp, snuggle_nester.hpp | **WARNING** | `ROT_CACHE_BINS` (nester + CPU evaluator) and `ROT_BINS` (GPU evaluator) are independent 360-valued constants. A one-sided edit silently produces bin-mismatch bugs. Should be a single shared constant in a common header. |
| Window class not unregistered | gpu_collision.cpp | **MINOR** | `UnregisterClassA` never called in `cleanup()`. Class leak if evaluator is constructed multiple times. |
| gpu_collision.cpp always compiled | CMakeLists.txt | **MINOR** | Not conditionally compiled, but works because non-GUI code is outside the guard. Intent is unclear to readers. |

---

## Verdict

The foundation is sound enough to build on. The two load-bearing cracks are the legacy WGL context creation (GPU path may silently get a wrong context type on strict drivers) and the forward-declaration contract in `gpu_collision.hpp` (will ambush the first person who writes a standalone unit test for the GPU evaluator). The ROT_BINS duplication is a maintenance trap waiting to be sprung. Fix those three before this goes to a wider review.
