# Building the Bitmap Nester (feature/concave-bitmap)

## Quick Reference

```bash
# Pre-flight (run after any session break or workspace change)
bash orca/check-build.sh clean

# Build and run tests
BUILD_CONSENT=yes bash orca/build.sh clean test

# Build full app + runnable snapshot
BUILD_CONSENT=yes bash orca/build.sh clean dll --snapshot --note=description
```

All commands run from the workspace root (`F:/Claude/`). Snapshots land in `orca/builds/clean/YYYYMMDD-HHMMSS-commithash-description/`.

---

## Worktree Layout

| Shortname | Path | Branch |
|-----------|------|--------|
| `clean` | `orca/clean/` | feature/concave-bitmap |
| `v1` | `orca/v1/` | feature/concave-arrange |
| `v2` | `orca/v2/` | feature/concave-arrange-v2 |
| `snuggle` | `orca/snuggle/` | feature/snuggle-radial |

Each worktree shares a single build harness (`orca/build.sh`) and a single snapshot root (`orca/builds/`). Snapshots are organized by worktree:

```
orca/builds/
  clean/YYYYMMDD-HHMMSS-hash-note/    ← full runnable installs
  v1/...
  snuggle/...
  (legacy snapshots at root level)     ← pre-consolidation, hash_note naming
```

## Dependencies

The clean worktree borrows pre-built dependencies from the snuggle worktree:

```
orca/snuggle/deps/build/OrcaSlicer_dep/usr/local
```

These deps were originally built at an older path. Filesystem junctions map old paths to new locations so baked-in cmake references still resolve. Run `check-build.sh --fix-junctions` if needed.

## Pre-flight: `check-build.sh`

Run before building, especially after a session break or workspace change:

```bash
bash orca/check-build.sh clean              # basic checks
bash orca/check-build.sh clean test         # also verify BUILD_TESTS=ON
bash orca/check-build.sh clean --fix-junctions  # auto-create missing junctions
```

Checks:
1. Build directory and .sln exist (cmake was configured)
2. `CMAKE_PREFIX_PATH` resolves to a real directory
3. `BUILD_TESTS=ON` (when target is test)
4. Junctions exist for any old-path deps referenced in the cmake cache
5. No competing msbuild/cl processes
6. Full install directory exists (exe + resources)

## Build Targets

| Target | What | Time (incremental) |
|--------|------|-------------------|
| `test` | Build + run libslic3r_tests | ~2 min |
| `lib` | libslic3r only | ~30 sec |
| `dll` | OrcaSlicer.dll | ~2 min |
| `gui` | Full app (orca-slicer.exe) | ~5 min |
| `all` | Entire solution | ~10 min |

A from-scratch build after cache nuke takes ~10 minutes with `-m:4`.

## Snapshots

`--snapshot` copies the full install directory — exe, DLLs, VC runtimes, and `resources/`. The snapshot is a complete runnable OrcaSlicer install.

```bash
# Snapshot to default location (orca/builds/clean/)
BUILD_CONSENT=yes bash orca/build.sh clean dll --snapshot --note=baseline

# Snapshot to custom location
BUILD_CONSENT=yes bash orca/build.sh clean dll --snapshot --snapshot-dir=some/other/path --note=baseline
```

The harness:
- Auto-upgrades `lib`/`dll` targets to `gui` (snapshot needs orca-slicer.exe)
- Prefers the full install dir (`build/OrcaSlicer/`) over the raw Release dir
- Copies VC runtime DLLs for standalone execution
- Removes build artifacts (.exp, .lib, .pdb)
- MD5-verifies `OrcaSlicer.dll` against source
- Fails loud if exe, DLL, or resources/ are missing

## Fixing a Broken Build

### Stale cmake paths (most common after workspace moves)

Symptom: LNK1181 errors referencing old paths like `orcaSlicer-snuggle`.

**Fast fix — junctions:**
```bash
bash orca/check-build.sh clean --fix-junctions
```

**Clean fix — nuke and reconfigure:**
```bash
cd orca/clean/build
rm -rf CMakeCache.txt CMakeFiles/
cmake .. -DCMAKE_PREFIX_PATH=../../snuggle-radial/deps/build/OrcaSlicer_dep/usr/local -DBUILD_TESTS=ON
```

Editing `CMakeCache.txt` alone is NOT sufficient — the generated `.vcxproj` files have their own baked-in paths from configure time.

### BUILD_TESTS is OFF

Symptom: test project not in solution, `libslic3r_tests` target not found.

```bash
cd orca/clean/build
cmake .. -DBUILD_TESTS=ON
```

The variable is `BUILD_TESTS`, not `BUILD_TESTING`.

## Rules

- **One builder at a time.** Concurrent builds on the same directory corrupt PDB locks and outputs.
- **Commit before building.** If the build fails, you need to distinguish your changes from environment problems.
- **Snapshot after every milestone.** Full install, not just test binaries.
- **`BUILD_CONSENT=yes` is required.** Safety gate against accidental builds.
- **Use `-m:4` (default).** Keeps the machine usable. `--afk` for all cores when unattended.
