# Frodo Review: Snuggle UI and User-Facing Behavior

Reviewer: Frodo (regular user, zero patience for confusion)
Files reviewed:
- `src/slic3r/GUI/GLCanvas3D.cpp` (lines 5922-6106, Snuggle UI panel)
- `src/libslic3r/Arrange/SnuggleArrange.cpp` (progress, errors, fallback)
- `src/slic3r/GUI/Jobs/ArrangeJob.cpp` (integration, cancellation, finish messages)

---

## Good

1. **Enable checkbox tooltip is excellent.** "Pack parts using their full 3D shapes instead of flat outlines. Produces tighter blob-shaped clusters. Best with 2-30 parts. Takes a few seconds." -- tells me exactly what it does, when to use it, and how long to wait. Perfect.

2. **Sequential printing warning is visible and clear.** Orange text, warning icon, tells me it will fall back. I won't be confused about why my arrangement looks normal.

3. **Rotation dropdown is well-designed.** Discrete options (Locked, 90, 45, 15, 5, 1 degree) instead of a free-entry field. Tooltip explains "relative to starting orientation." No way to enter nonsense.

4. **Resolution slider tooltip is helpful.** "2.0 mm is a good default. 0.5 mm for tight packing." Gives me a reference point.

5. **Overflow fallback is seamless.** If Snuggle can't place everything, leftover parts silently go to the standard arranger. The user sees "Snuggle placed 7 of 10 parts. Some overflow to next plate." -- clear, not alarming.

6. **Completion message includes bed utilization percentage.** Nice touch. Gives me a sense of how well it packed.

7. **Advanced settings hidden in a collapsible tree.** Population, Generations, Timeout, Max parts -- all tucked away. Casual users never see them.

8. **Cancellation works.** Stop condition is wired through, checked every generation. Cancel shows "Arranging canceled." -- standard OrcaSlicer behavior.

---

## Problems

### P1: "Part gap" slider has no tooltip

The slider at line 5950-5962 has no `IsItemHovered` / `SetTooltip` block. Every other control in the section has one. A user dragging this slider has no idea what the units are (mm? pixels? percentage?). The slider format string is `"%4.1f"` -- no unit suffix. Compare with the Resolution slider which shows `"%3.1f mm"`.

**What I'd expect:** Tooltip saying something like "Minimum gap between parts in mm. 0 = touching. 5 mm default." And the slider format should show "mm" like Resolution does.

### P2: "Compact after arrange" checkbox appears TWICE

It appears at line 5964 (main Snuggle section) AND at line 6080 (inside the "Effort / Advanced" tree node). Both control the same `snuggle_compact` setting. The tooltips differ:

- Main section: "Jiggle parts toward center after arrangement to close gaps."
- Advanced section: "Jiggle parts toward center after GA. Usually not needed -- the GA's fitness function already optimizes for tight clusters."

This is confusing. If I check the top one, does the bottom one update? (Probably yes, since they share the same setting variable.) But the advanced tooltip actively discourages its use ("usually not needed") while the main one presents it neutrally. Pick one location. The advanced section's tooltip is more informative -- put it there and remove the main-section duplicate.

### P3: "Snuggle" is not a self-explanatory name

The section header says "Snuggle 3D Arrangement" in teal. No "(experimental)" or "(beta)" label. No version indicator. A user who's never heard of this feature has to hover over the checkbox to learn what it does. The header alone tells me nothing -- "3D Arrangement" could mean anything.

Suggestion: "Snuggle (experimental) -- 3D-aware part nesting" or at minimum add a brief subtitle.

### P4: Progress bar messages are developer-speak

During arrangement, the status bar shows messages like:

```
Arranging (Snuggle CPU gen 12/30 OK)
Arranging (Snuggle CPU gen 12/30 3 collisions)
```

"gen 12/30" means nothing to a normal user. "3 collisions" is alarming -- am I doing something wrong? "CPU" is an implementation detail (there's apparently a GPU backend too).

**What I'd expect:** Something like "Arranging parts... 40%" or "Arranging parts... optimizing layout (40%)". The generation count and collision status belong in the log, not the status bar.

### P5: No tooltip on "Part gap" means users don't know 0 = touching

Setting the gap to 0.0 is valid (the clamp range is 0.0-20.0). But with no tooltip, a user might not realize 0 means literal contact between parts. For FDM printing, 0 gap is almost never what you want -- parts will fuse. A tooltip warning "0 mm = parts may touch. For FDM, 2-5 mm recommended." would prevent support tickets.

### P6: "Max parts" limit has no tooltip and no user-facing feedback when hit

In the Advanced section, "Max parts" defaults to 200 (clamp 2-500). When you exceed this limit, Snuggle silently marks all parts as unarranged and falls back to the standard arranger. The only evidence is a BOOST_LOG_TRIVIAL(warning) in SnuggleArrange.cpp line 57. The user sees... nothing different. Their arrangement just looks like regular arrange, not Snuggle.

This should show a visible warning: "Too many parts for Snuggle (limit: 200). Using standard arrangement."

### P7: "Multi-plate overflow" checkbox has no tooltip

Line 6074. What does this do? Does it let Snuggle spill to a second plate? Does it prevent it? The name is ambiguous. A tooltip is mandatory.

### P8: Timeout has no tooltip

Line 6056. The Timeout input (2-300 seconds) has no tooltip explaining what happens when it fires. Does it use whatever the GA has found so far? Does it fall back? Does it cancel? The user is left guessing.

### P9: Voxelization failure is invisible

When a part fails to voxelize (SnuggleArrange.cpp line 188-191), it gets silently marked `bed_idx = -1` and falls to the overflow arranger. The user never learns that their part couldn't be 3D-nested. For weird geometry, this could happen reliably and the user would think Snuggle is broken when it's just one troublesome part.

### P10: The teal header color may be invisible in dark themes

`ImVec4(0.00f, 0.59f, 0.53f, 1.00f)` is a dark teal. On dark backgrounds it reads fine. On some custom themes with dark-teal or cyan backgrounds, it could vanish. Minor, but worth noting -- the sequential-printing warning uses orange (0.92, 0.60, 0.0) which is more universally visible.

---

## Summary

The UI is well-structured overall. The collapsible advanced section, rotation dropdown, and overflow fallback are all good design. The main gaps are:

1. **Missing tooltips** on Part gap, Max parts, Multi-plate overflow, and Timeout (4 controls with no hover help)
2. **Duplicate checkbox** (Compact after arrange appears twice with conflicting guidance)
3. **Developer-facing progress messages** that will confuse normal users
4. **Silent failures** when voxelization fails or max-parts limit is hit -- no user-visible feedback

None of these are blockers, but items P1, P2, and P4 would meaningfully improve the experience with minimal code changes.
