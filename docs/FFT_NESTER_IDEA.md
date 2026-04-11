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

## Companion idea: cheap impossibility certificates + plate consolidation

Captured 2026-04-11 from a follow-up design discussion. This is the
"when not to think" half of the FFT fork — the FFT proposal above
covers *how* to search the placement space; this section covers *when*
to call the search and when to declare "obviously done" without
running it. Both belong in the same fork.

### The core insight

The human isn't optimizing — they're recognizing when optimization is
pointless. Computers do the opposite: they optimize indiscriminately
because the framework gives them no cheap "is this worth thinking
about" test. Add that test, layered cheapest-first, and most
refinement work is skipped before it starts.

### Start with a bound, not a solution

Compute the theoretical minimum plate count `K_min` before placing
anything:

```
K_min = ceil(sum(part_silhouette_area) / plate_area)
```

That's a hard lower bound — you cannot do better. After the greedy
pack, if the producer hits exactly `K_min` plates, you are provably
done. No rearrangement pass, no second guessing, no FFT calls. This
single check catches both the easy case (lots of small parts) and
the trivially-hard case (one huge part per plate) for free.

### Pair-level "obviously won't fit" certificate

The "two diagonal rectangles" intuition is a constant-time bounding-
box test. For each pair of parts on separate plates, check whether
*any* combination of 90-degree rotations places their bounding boxes
on a single plate. The check is four inequalities per rotation
combo (16 total for 4 rotations each):

```
bbox_cofit(A, B, plate) =
  any rot_a in {0,90,180,270}, any rot_b in {0,90,180,270} where:
    (rot(A).w + rot(B).w <= plate.w  AND  max(rot(A).h, rot(B).h) <= plate.h)
    OR
    (rot(A).h + rot(B).h <= plate.h  AND  max(rot(A).w, rot(B).w) <= plate.w)
```

If no rotation combo passes, the exact silhouettes cannot co-fit
either. This is the "I can see these don't fit together" moment a
human has instantly. It is constant-time per pair and prunes the
expensive FFT work before it starts.

### Layered impossibility certificate

Cheapest first. Each layer says either "definitely impossible, skip"
or "maybe, ask the next layer."

**Layer 0 — area saturation.** If the current plate's utilization is
above 85-90 percent, no rearrangement gains a meaningful win. Skip
the plate. Tunable threshold. This is the "I can see this plate is
basically full" check.

**Layer 1 — bounding-box bin pack.** Ignore silhouette shapes
entirely. Can the axis-aligned bounding boxes of the parts on the
*next* plate fit into the free bounding-box of the *current* plate?
This is 2D bin packing of rectangles, NP-hard in general but for
2-5 candidate parts it brute-forces in microseconds. If the bboxes
don't fit, the real shapes can't either.

**Layer 2 — FFT correlation.** Only run this if layers 0 and 1
return "maybe." Compute the feasible-placement mask for each
candidate part on the target plate at every rotation. If every mask
is empty (every pixel is non-zero), the placement is certified
impossible. Move on.

The point is that layer 2 — the only expensive layer — runs on
maybe one in twenty pair attempts. The other nineteen are killed by
layers 0 and 1 in nanoseconds.

### Suspicion-weighted refinement budget

After the greedy pass, compute a "suspicion score" per plate:

```
suspicion(plate) =
    plate.utilization is in [50%, 80%]
  AND
    next_plate has at least one part whose bbox bbox_cofit_passes(plate)
```

Plates above 95 percent utilization score zero — they're not worth
revisiting. Plates between 50 and 80 percent with migratable
candidates score highest. Spend the refinement time budget
proportional to suspicion, not uniformly. Plates that already look
"nearly full" get no thinking time at all.

### Plate consolidation loop

Work backwards from the last plate. Try to empty each plate by
migrating its parts to earlier plates. The moment a plate empties,
plate count drops by one — that is the biggest visible win a user
notices.

```
for plate_i from K-1 down to 1:
    if bbox_cofit_test(plate_i.parts, any earlier plate) == IMPOSSIBLE_FOR_ALL:
        continue        # Layer 1 prune
    for each part P on plate_i, in descending area order:
        for each earlier plate_j in [0, i):
            if layer0(plate_j) returns SATURATED: continue
            if bbox_cofit_test(P, plate_j) == IMPOSSIBLE: continue
            mask = fft_correlate(P, plate_j, all_rotations)
            if mask has any zero pixel:
                move P to plate_j
                if plate_i is now empty:
                    eliminate plate_i; K -= 1
                    break
    if time_budget_exceeded: break
```

The two `continue` lines are the cheap pruning — most iterations
skip the FFT call entirely. Within-plate density improvements are
secondary; eliminating an entire plate is the headline metric.

### Termination

Three stopping conditions, whichever hits first:

1. `K == K_min`. Provably can't improve. Area-theoretic optimum.
2. Every remaining plate pair fails the bbox co-fit test. Provably
   can't consolidate further. This is the "human glances at the
   layout and shrugs" condition.
3. Wall-clock budget exhausted. Default 2-5 seconds for refinement.

Condition 2 is the most important one to mechanize — it is what a
human does without thinking, and it costs almost nothing
computationally.

### Near-miss detection (the one place local search pays off)

The FFT correlation map gives you more than a binary "fits / doesn't
fit." The numerical values are the count of overlapping pixels at
each candidate position. Near-zero values mean "almost fits, just a
few pixels of overlap somewhere."

When the best mask value for a candidate placement is small (say
below a few percent of part area), do *not* give up immediately.
Instead:

1. Identify the blocking already-placed part (the one whose pixels
   the candidate would overlap).
2. Try nudging the blocking part by a small amount in 4-8 directions,
   each time re-running the candidate's correlation.
3. Bail after at most 3 nudges per candidate.

This is the "if I just shifted that other one by 5 mm..." reflex.
Bounded local search, costs almost nothing, occasionally produces
the move that consolidates a plate.

### Combined complexity profile

- Greedy initial pack: O(N * R * G^2 log G) — linear in parts
- K_min bound: O(N) — single pass
- Bbox co-fit certificates: O(K^2 * N) constant-time checks — tiny
- FFT refinement (only on layer-0+1 survivors): O(M * R * G^2 log G)
  where M is migration attempts, bounded by time budget
- Near-miss nudges: at most 3 * R per candidate, capped

**Total wall clock: greedy time + fixed refinement budget (2-5 s).**

Nothing exponential. Nothing combinatorial. The certificate layers
prune the expensive FFT work so aggressively that *most refinement
calls return without ever running an FFT*. The cost of asking "is
this worth thinking about" is negligible compared to the cost of
one unnecessary correlation.

### Why this matters more than the FFT itself

The FFT proposal in the previous section is about the inner loop
running faster. This section is about not running the inner loop at
all when there's nothing to gain. In practice the second win is
larger — most plates after a reasonable greedy pass are either
already optimal (catch with K_min and layer 0) or provably
unfixable (catch with bbox co-fit). The FFT only earns its keep on
the small middle slice where rearrangement might actually work, and
even there the certificate stack has already filtered down to the
moves worth attempting.

The fork's pitch lines up nicely:

- **What is faster:** FFT correlation vs scored brute-force scan.
- **What is smarter:** layered certificates + plate consolidation
  vs single greedy pass.
- **What is the user-visible result:** plate count drops by one
  more often than the current nester achieves, because rearrangement
  actually runs on the cases where it can win, and refinement time
  is spent where it pays off.

## Companion idea (3 of 3): three-tier parallelism strategy

Captured 2026-04-11 from a follow-up discussion. The first two
companion sections cover the algorithm (FFT correlation) and the
control flow (certificate-pruned consolidation). This third section
covers how to schedule all of it onto multiple cores so the
"refinement budget" is hidden inside the packing wall-clock instead
of added to it.

### The fundamental tension

Greedy packing is sequential — each placement depends on the
previous one. You cannot parallelize within a single packing pass
without changing the algorithm. But you can parallelize *across*
strategies and *within* the expensive subroutines.

### Tier 1 — FFT parallelism (free, automatic)

Each FFT correlation for a candidate rotation is independent. With
12 rotations per part, that's 12 independent FFTs. FFTW or PocketFFT
will thread within a single transform, but the bigger win is
dispatching all 12 as independent tasks. On an 8-core machine, all
12 rotations complete in the wall-clock time of 2.

This requires zero algorithm changes — the greedy packer just gets
its answers faster. Every approach below benefits automatically.

### Tier 2 — parallel permutation search (the big quality win)

The greedy packer's quality depends on placement order. Different
orderings produce different results. They are completely
independent — perfect parallel work units, no shared state during
execution.

**Primary path.** One thread runs the deterministic "smart" ordering
(area-descending, the best deterministic heuristic). It starts
immediately and is the guaranteed baseline.

**Speculative paths.** The remaining N-1 threads each run a seeded
*perturbation* of the smart ordering — swap 2-3 parts, shuffle
within size-class buckets, randomize ties. These are exploring the
ordering neighborhood. They are not fully random — the smart
ordering is the center of mass and the perturbations are nearby.

**No communication during execution.** No locks, no sync, no shared
state. Each thread has its own bitmap, its own composite plate
state, its own answer. They only converge at the end.

**Convergence.** When all paths complete, compare results by total
plate count (primary metric) and within-plate density (tiebreaker).
Adopt the best. The wall-clock cost is exactly one packing pass —
because they all ran simultaneously — but you've effectively
searched N orderings.

### Tier 3 — pipelined packing + certification + consolidation

This is the subtle one. The "is this worth improving" certificate
work and the actual consolidation work happen *while* the packer is
still running on later plates. There is no separate refinement
phase that the user waits for.

```
Thread A (packer):       plate 1 -> plate 2 -> plate 3 -> plate 4 -> ...
Thread B (certifier):    idle    -> cert 1   -> cert 2   -> cert 3   -> ...
Thread C (consolidator): idle    -> idle     -> try 2->1 -> try 3->? -> ...
```

**Thread A** is the greedy packer, churning through parts plate by
plate.

**Thread B** trails by one plate. When plate K is "done" (the
packer has moved on), Thread B computes plate K's utilization, the
free-region bounding box, and pre-computes the feasible-placement
mask for the free space at all rotations. This is the certificate
work from the section above. Each finished plate is fully
characterized within milliseconds of the packer leaving it.

**Thread C** trails by two plates. It takes Thread B's
characterization and attempts consolidation — can any part from
plate K+1 move backwards into an earlier plate? Because Thread B
has already cached the feasibility data, the migration check is
just a lookup. If a migration succeeds, the plate count drops.

**The pipeline means consolidation results are ready by the time
the packer finishes.** The user-perceived wall-clock for "pack +
refine" equals the wall-clock for "pack" alone. The refinement is
free in clock time, bounded in cores.

### Handling shared state without locks

The pipeline has one data dependency: if Thread C migrates a part
from plate K+1 to plate K, that changes plate K's bitmap, which
invalidates Thread B's characterization of plate K.

Use **versioning, not locks**. Each plate carries a version
counter. Thread B's characterization is tagged with the version it
was computed against. If Thread C bumps the version (by mutating
the plate), Thread B's stale tag tells it to re-characterize on the
next pass.

In practice, successful migrations are rare (most attempts fail at
the bounding-box certificate, layer 1 from the previous section),
so re-characterization runs almost never. The common case is a
single linear pipeline pass with no rework.

### Launch strategy by core count

Detect at startup. Let `C` be the available core count.

```
C == 1:    sequential greedy, no parallelism, skip refinement
C in 2..3: 1 primary + 1-2 speculative permutations
           refinement runs serially after packing (no pipeline)
C in 4..7: 1 primary + 1 certifier + remaining cores on speculative
           consolidator runs sequentially after pipeline drains
C >= 8:    1 primary + 1 certifier + 1 consolidator pipeline
           remaining cores on speculative permutations
```

The progression maps cores to the highest-leverage parallel work
that's still worth the bookkeeping cost.

### Work unit balancing

The packer for `N` parts takes roughly `N * R * FFT_cost`. A
speculative permutation takes the same. The certifier per plate is
much cheaper — one utilization computation plus a few bounding-box
checks. The consolidator per plate is variable — might be zero work
(bbox says impossible) or one FFT (try the migration).

**Risk:** the certifier and consolidator starve because they're
waiting on the packer to finish plates.

**Fix:** the packer emits "plate done" signals into a lock-free
queue. The certifier pulls from the queue and always has work as
long as the packer is running. If the certifier runs out of work
(packer is stalled on a hard plate), it steals a speculative
permutation task. Work-stealing keeps all cores busy.

### Cross-thread early termination

If the primary packer achieves `K_min` plates (the area-theoretic
optimum from the previous section), broadcast a cancel to all
speculative threads. Their work is now waste.

Symmetrically, if a speculative thread finds a `K_min` solution
before the primary finishes, cancel the primary.

The cancel is a single atomic flag, polled at plate boundaries.
No complex synchronization, no thread state machine.

### What the user sees

Packing result appears in the same wall-clock time as a single
greedy pass — because the speculative paths and the pipeline both
overlap with the primary packer. But the result is the *best of N*
orderings with consolidation already applied. The "few extra
seconds for refinement" the user wanted cost zero extra seconds
because they ran in parallel with the packing they were already
waiting for.

### Memory footprint caveat

Each 480x480 complex-float bitmap is ~1.8 MB. With `C` threads
each holding `R` rotation bitmaps and a plate bitmap, total memory
is roughly `C * R * 2MB`. At 8 cores * 12 rotations = ~200 MB.
Fine on a desktop, worth watching on constrained environments.

If memory is tight: share the *immutable* part rotation bitmaps
read-only across speculative threads (each speculative thread
still needs its own *occupied* plate bitmap because they're
packing independently), and drop speculative-thread rotation
candidates from 12 to 4. Speculative threads are already
exploring a different ordering, so rotation diversity inside each
thread matters less.

### Combined wall-clock budget

Sequential greedy was: pack_time + refine_time
This proposal is: max(pack_time, refine_time + pipeline_lag) ≈ pack_time

The refinement is functionally free in clock time. The cost is
linear in cores (and capped by memory). On any modern desktop the
nester finishes in the same time it does today but with strictly
better answers.

### Why this matters relative to the FFT itself

These three companion sections form the actual fork pitch:

| Section | Lever | User-visible win |
|---|---|---|
| 1 — FFT inner loop | Algorithm | Independent of part complexity, scales to high resolution |
| 2 — Certificates + consolidation | Control flow | Plate count drops by one more often, refinement runs only where it pays |
| 3 — Three-tier parallelism | Scheduling | The whole improvement happens within the existing wall-clock budget |

Section 1 alone does not justify the fork — at our current
resolution and rotation count it's roughly cost-neutral against the
scored scan. Section 2 alone helps but is limited by single-threaded
greedy time. Section 3 is what makes the combination feel
qualitatively different to the user — the same wait time produces
visibly tighter plates with provable lower bounds backing them up.

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
