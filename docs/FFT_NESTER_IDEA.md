# FFT-Correlation Bitmap Nester — Reference Idea

**Status:** Pinned for post-Sprint-3 exploration. Not on the current PR-target
roadmap. Captured 2026-04-11 from a design discussion. The current
`BitmapNester` shares the foundational ideas (rasterize everything, AND
collision, plate composite bitmap) but uses a brute-force scored scan for
placement search and a single greedy outer loop. The FFT approach below is
qualitatively different in two places — placement search and outer-loop
search strategy — and is worth a dedicated fork to evaluate.

## The proposal verbatim

> Rasterize everything and stay rasterized. If you commit to bitmaps, you
> never need polygon booleans at all. The silhouette projection becomes: for
> each mesh triangle, project to XY, rasterize the 2D triangle into a bitmap.
> Standard scanline fill per triangle — linear in triangle count, no boolean
> unions, no Clipper. You end up with a 1-bit bitmap per part. Done.
>
> Collision detection is bitwise AND. Two parts overlap if their bitmaps,
> offset to candidate positions, have any overlapping set bits. That's a
> single pass of `&` across scanlines. At 0.5mm resolution on a 240mm plate
> that's 480x480 = ~29K uint64s per part. Collision check is ~450 AND
> operations. That's nanoseconds. You can brute-force thousands of candidate
> placements per part without blinking.
>
> NFP equivalent in bitmap space: the summed-area trick. Instead of computing
> geometric NFPs, do this: for all already-placed parts, maintain a single
> composite "occupied" bitmap. For the new part's bitmap, compute the
> feasible-placement mask by sliding the part bitmap across the plate
> bitmap. This is a 2D correlation — and here's where it gets good: 2D
> correlation via FFT is O(N^2 log N) where N is the plate resolution,
> independent of part complexity. Part has 10 vertices or 10,000 vertices?
> Same cost. The FFT doesn't care about geometric complexity at all.
>
> Concretely:
>   1. Plate bitmap P (480x480, 1-bit but stored as float for FFT)
>   2. Part bitmap Q (small, padded to plate size)
>   3. Correlation C = IFFT(FFT(P) * conj(FFT(Q)))
>   4. Any pixel in C where the value equals zero is a valid placement position
>
> That gives you the complete feasible placement map for one part in two
> FFTs. No polygon intersection, no NFP, no Clipper. And it's the same cost
> regardless of how many parts are already placed — they're all baked into
> one bitmap.
>
> Placement strategy stays simple. Scan the feasible-placement mask bottom-
> left (or whatever heuristic) for the first valid pixel. That's your
> placement. Update the composite bitmap with an OR. Next part.
>
> Rotation handling without blowup. Precompute each part's bitmap at your N
> candidate rotations (say 12). That's 12 rasterizations per part, done once
> upfront, linear in part count. During placement search, you run the FFT
> correlation for each rotation candidate and pick the one that gives the
> tightest bottom-left placement. 12 FFTs per part. Still fast — a 480x480
> FFT is ~1ms on any modern CPU.
>
> The whole pipeline is near-linear in part count:
>   1. Rasterize all parts: O(P) where P = total triangle count across all parts
>   2. For each of N parts, do R rotation candidates x 1 FFT correlation:
>      O(N * R * G^2 log G) where G = grid resolution
>   3. Total: O(P + N * R * G^2 log G)
>
> No polygon intersections. No convex decomposition. No Minkowski sums. No
> combinatorial search. Nothing exponential. The only quadratic term is the
> plate resolution squared, which is a fixed constant you choose once.
>
> Where this falls down and what to do about it. FFT correlation gives you
> the valid placement mask but not the optimal placement. "Optimal" in the
> packing sense means minimizing wasted space, which is still a hard
> combinatorial problem over placement order. But the practical answer is
> the same one that works: try a few random permutations of placement order
> (say 50-200), run the greedy bottom-left-via-FFT for each, keep the best.
> Each permutation is O(N * R * G^2 log G). Totally parallelizable. Fixed
> time budget, predictable wall-clock, quality scales with time you're
> willing to spend.
>
> One implementation note on the bitmap resolution question. You don't need
> uniform resolution. Rasterize at coarse (1mm) for the placement search,
> then do a single fine-resolution (0.1mm) validation pass on the final
> arrangement to check for sub-pixel overlaps. The validation is one AND per
> part pair — linear.
>
> This whole approach is basically what chip placement in EDA has converged
> on, just adapted to irregular 2D parts instead of rectangular cells.

## How close is the current `BitmapNester` to this?

| Layer | Current `clean/src/libslic3r/BitmapNester.hpp` | FFT proposal | Status |
|---|---|---|---|
| Rasterize all parts up front | Yes, per rotation, into `RotCache` | Same | Match |
| Collision = bitmap AND | Yes, `collides()` is word-level uint64 AND | Same | Match |
| Single composite plate bitmap | Yes, `plate_items[p]` is the OR-accumulator | Same | Match |
| Skip polygon booleans entirely | **No** — `ArrangeJob::project_mesh` -> `union_ex` (Clipper) -> `concave_regions` -> rasterize. The Clipper union is the only polygon-boolean step. | Yes — direct per-triangle rasterize, no union | **Different** |
| Placement search | Brute-force coarse-stride scored scan with refine. Score = cluster-bbox area growth + anchor distance tiebreaker. ~16k position tests per part per rotation at default stride. | FFT correlation produces the entire feasible-placement mask in two FFTs per rotation. Pick the bottom-leftmost (or anchor-nearest) zero in the mask. | **Different** |
| Cost scaling | O(R * stride_positions * part_words). Linear in part complexity (part_words). | O(R * G^2 log G). **Independent of part complexity.** | Different scaling regime |
| Outer loop | **Single greedy pass**, sort by descending area. | **Random-permutation restarts** (50-200), keep best of all. Embarrassingly parallel. | **Different** |
| Multi-resolution | Coarse stride + ±coarse refine in one pass. Same resolution throughout. | Coarse-resolution placement search (1mm) + fine-resolution validation (0.1mm) at the end. | Partial |
| Per-rotation pre-rasterization | Yes (4 rotations, would extend to 8 or 12 cheaply) | Yes (12 rotations) | Match |

**Verdict:** the foundation is the same. The two qualitative differences are
the inner loop (brute scan vs FFT correlation) and the outer loop (single
greedy vs random restarts). Both differences are isolated enough that they
could be developed in a separate fork without disrupting the upstream PR.

## Why this is worth a fork (and why it isn't blocking now)

**Worth doing:**
- Asymptotically wins at high resolution and many rotations. At 1024x1024
  grid + 12 rotations + 100 parts, our scored scan does ~256 stride
  positions x 12 rotations x ~16K word ANDs per position = ~50M ops per
  part. FFT does 12 FFTs of 1024x1024 ≈ 12 * 10M = 120M ops per part. Same
  order today. At 2048x2048 + 24 rotations the scored scan blows up
  quadratically while FFT only grows by log factor.
- Random-restart outer loop has well-known quality lifts over single greedy
  for irregular bin packing — published BLF + random-restart hybrids beat
  pure greedy by 5-15% on standard ESICUP benchmarks.
- The composite-plate bitmap means cost is independent of items already
  placed, which matters for big plates.
- Predictable wall-clock budget. Our scored scan's runtime depends on how
  early it finds a clear position; FFT runtime is constant per attempt.

**Why not now:**
- Adds an FFT dependency. The current branch's contract with upstream is
  "no new external libraries" — we'd need to use a stdlib FFT or write a
  small one. Real implementations want FFTW or PocketFFT. Either way it's a
  new dep that has to be argued for.
- Random-restart parallelism wants TBB integration through the arrange job,
  which is already not great in the current code (cancel-handling, progress
  reporting, status updates). Best done as a clean rewrite, not a graft.
- The current scored scan + cluster-bbox-growth scoring is already producing
  reasonable density on the tetromino benchmarks (scenarios 1, 2, 4, 5
  above 0.65-0.95). The user-visible "loose grid" complaint is fixed.
  Sprint 2's job is to land the PR, not to chase the next performance
  ceiling.
- The Clipper `union_ex` step that the FFT approach skips is already cached
  per-instance in `ArrangeJob::silhouette_cache`, so its amortized cost
  across an arrange call is small. The "skip polygon booleans entirely"
  benefit is mostly aesthetic for our usage pattern.

## Suggested fork structure (when we get there)

Branch name: `feature/fft-nester` or `feature/concave-bitmap-fft`.

Phase 1: prove the core math.
- Add a `tests/libslic3r/test_fft_nester.cpp` that hard-codes a 64x64 grid
  with two simple shapes, runs FFT correlation, asserts the feasible mask
  matches the brute-force AND result for every position.
- Use a small handwritten radix-2 FFT (no external dep) at 64x64 to
  validate the algorithm. ~200 lines.

Phase 2: scale up.
- Add PocketFFT or KISS FFT (both header-only) as a tracked third-party
  source under `deps/`. Argue it with upstream.
- Rewrite the placement loop in a new file `FftNester.cpp` keeping the same
  `ArrangePolygons` API as `BitmapNester::arrange`. Same dispatch in
  `Arrange.cpp` adding a third strategy.

Phase 3: outer loop.
- Random-permutation restarts with TBB parallel-for, fixed time budget, keep
  best by total cluster-bbox area. Default 50 restarts.

Phase 4: head-to-head benchmarks.
- Run both nesters against the tetromino harness (`docs/PACKING_BENCHMARK.md`)
  on identical inputs. Publish density and wall-clock numbers in
  `docs/FFT_VS_SCAN.md`.

Phase 5: decide.
- If FFT nester wins on density AND is competitive on wall-clock, deprecate
  the scored scan. Otherwise keep both as user-selectable strategies.

## Notes from current code that the fork should preserve

- `ArrangePolygon::concave_regions` (multi-island) — required for parts
  whose top/bottom projections diverge.
- The `convex_bm` companion bitmap that agrees with
  `PartPlate::check_outside` for exclude validation. The FFT path can
  compute this once per rotation the same way.
- `align_to_y_axis` pre-rotation composition (`base_rot + best_rot`).
- `stopcondition` cancellation respect.
- `MAX_PLATES` overflow handling.
- The log-line distinction between "extraction failed" and "no matching
  ModelInstance" so silent fallbacks remain visible in the GUI log.

These are functional contracts the upstream arranger expects, not nester
internals. Any new strategy must honor them.
