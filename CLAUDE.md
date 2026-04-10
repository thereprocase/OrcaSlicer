# OrcaSlicer — In-Tree Dev Guide

C++ slicer forked from Bambu Studio. wxWidgets GUI, CMake build, Catch2 tests.

## Build

```bash
# Windows
cmake --build . --config Release --target ALL_BUILD -- -m

# macOS
cmake --build build/arm64 --config RelWithDebInfo --target all

# Linux
cmake --build build --config RelWithDebInfo --target all
```

- CMake 3.13–3.31.x, primary build dir: `build/`, deps in `deps/build/`
- Windows: Visual Studio generators. macOS: Xcode (or Ninja with -x). Linux: Ninja.

## Test

```bash
cd build && ctest                          # all tests
cd build && ctest --output-on-failure      # verbose
ctest --test-dir ./tests/libslic3r        # single suite
```

Test dirs: `tests/libslic3r/` (21 files), `tests/fff_print/` (12), `tests/sla_print/` (4), `tests/libnest2d/`, `tests/slic3rutils/`.

## Code Style
- C++17 with selective C++20. PascalCase classes, snake_case functions/variables.
- `#pragma once`. Smart pointers, RAII. TBB for parallelization.

## Key Entry Points
- `src/OrcaSlicer.cpp` — application startup
- `libslic3r/Print.cpp` — slicing pipeline orchestrator
- `PrintConfig.cpp` — all print/printer/material settings

## Architecture
For full architecture details (libraries, algorithms, file formats, dependencies, development workflow): read `ORCASLICER_REFERENCE.md`.
