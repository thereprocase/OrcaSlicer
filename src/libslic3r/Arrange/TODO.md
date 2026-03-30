# BitmapArranger — Open Items

## Performance

### Min/max pre-filter for 3D items
Phase 2.5 builds min/max bitmaps (AND/OR across rotations) to skip
per-rotation collision checks. Works for 2D items (75/75 valid). Fails
for 3D items (0/111 valid) because rotated SliceStacks have different
bitmap dimensions — the centering logic only handles BitmapItem, not
SliceStack. Fix: center each rotated stack's bottom slice in the
envelope frame, same as the 2D path does for rotated bitmaps.
Impact: ~5-10x speedup on 3D Phase 3a for large item counts.

### 3D center-out fallback dominates on full plates
When skyline placement succeeds on the bottom slice but collides_3d
rejects it (higher slices conflict), we fall through to find_placement_3d
which does O(bed^2 x n_slices) center-out scan. This fires frequently
on crowded plates and dominates Phase 3a time. Options:
- Try multiple skyline candidates before falling back (top-K per rotation)
- Use 2D pre-filter in the fallback to skip positions that fail on slice 0
  (already implemented but the fallback itself is still O(bed^2))
- Accept that 3D nesting of non-prismatic parts on full plates is inherently
  slower than 2D — the Z-slice collision is the cost of the feature

### Phase 4 compaction for 3D is slow
3D compaction rebuilds skylines per-plate per-item. The maintained skyline
helps placement but compaction still does per-rotation collision with
collides_3d. Parallel plate search (TBB) helps but doesn't eliminate the
per-slice cost.

## Features

### Fill bed with copies
BitmapArranger is well suited: take one part, repeat-stamp via skyline
until find_placement_skyline returns nullopt. Add fill_mode flag to
ArrangeParams. Support "fill with N copies" (cap) and "fill all plates"
(overflow). FillBedJob already calls our arrange pipeline — just needs
the estimation and candidate creation changes (partially done).

### Non-rectangular bed testing
Polygon bed overload is implemented (rasterize outside-bed as exclude).
CircleBed approximated as 64-gon. Not yet tested with actual delta
printer profiles. Need a test case with a circular bed to verify the
mask doesn't clip items at the curved boundary.
