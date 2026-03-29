# Breaking Changes Investigation — b810ec3459

## What happened
Build `b810ec3459_reviewed_final` crashes on launch: "Unhandled unknown exception; terminating the application." during "Loading configuration..."

## What was in the commit
The commit (`b810ec3459`) combined:

1. **Comment cleanup** (~50 edits) — comment-only changes, no code
2. **`version()` static method** added to BitmapArranger.hpp (public)
3. **`#include "BitmapArranger.hpp"`** added to GLCanvas3D.cpp
4. **`BitmapArranger::version()` call** in GLCanvas3D.cpp UI rendering
5. **Defaults button** in GLCanvas3D.cpp (ImGui, no cross-module calls)
6. **Resolution floor** (0.3mm clamp in UI slider + backend)
7. **Exclusion padding** changed from `res` to `scaled(0.5)`
8. **3D skyline cache** (PlateState3D gains `skyline` vector, init/update)

## What likely broke it
Items 2-4 are the prime suspects. Adding a public method to the hpp caused a
full rebuild of everything that includes Arrange headers (transitively, most of
libslic3r_gui). The crash happens before any arrange code runs (during config
load), suggesting a static initialization or DLL linking issue.

Possible causes:
- Cross-module static access: `BITMAP_ARRANGE_VERSION` is file-scoped static
  in BitmapArranger.cpp, `version()` returns a pointer to it, GLCanvas3D.cpp
  (in libslic3r_gui.lib) calls it. If the DLL initialization order doesn't
  guarantee BitmapArranger.cpp's static init before GLCanvas3D's first use,
  the pointer could be garbage.
- Header inclusion cascade: BitmapArranger.hpp pulls in Arrange.hpp,
  BoundingBox.hpp, ExPolygon.hpp, ClipperUtils.hpp, etc. into GLCanvas3D's
  translation unit. Possible ODR violation or conflicting definitions.
- Simply a build artifact: the DLL was copied correctly (same size, 1 min
  apart) but maybe the PDB or some other artifact is stale.

## Items 5-8 were working
Build `3f401c638b_skyline_cache` (the test build before the commit) had items
5-8 working. It was built from uncommitted changes on top of `3f401c638b`. The
crash was introduced by adding items 2-4 in the commit.

## What to do next
1. Apply ONLY the comment changes (item 1) — already in progress
2. Apply items 5-8 one at a time, building and testing each
3. For the version tag: use a simple hardcoded string in GLCanvas3D.cpp
   instead of a cross-module call. No hpp change, no include.
4. If that works, the root cause was the cross-module static access.
   Document and close.
