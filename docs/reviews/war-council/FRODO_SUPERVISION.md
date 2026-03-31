# Frodo Supervision Report -- Snuggle UI Review

Reviewer: Frodo (regular user perspective)
Date: 2026-03-31
Files reviewed:
- `D:\ClauDe\orcaSlicer-snuggle\src\slic3r\GUI\GLCanvas3D.cpp` (lines 5922-6158)
- `D:\ClauDe\orcaSlicer-snuggle\src\libslic3r\Arrange\SnuggleArrange.cpp` (all)
- `D:\ClauDe\orcaSlicer-snuggle\src\slic3r\GUI\Jobs\ArrangeJob.cpp` (snuggle sections)

---

## 1. Controls and Tooltips

### Present and good
- **Enable Snuggle** checkbox (line 5928): Has tooltip. Clear, mentions "2-30 parts" and "few seconds." Good.
- **Compact after arrange** checkbox (line 5964): Has tooltip. "Jiggle parts toward center" is plain language.
- **Rotation** dropdown (line 5989): Has tooltip explaining Locked and 90-degree options. Good.
- **Resolution** slider (line 6008): Has tooltip. "Lower = more precise but slower" is helpful.
- **Population** input (line 6028): Has tooltip. Fine for an advanced section.
- **Generations** input (line 6042): Has tooltip. Fine.

### MISSING tooltips
- **Part gap** slider (line 5950-5962): NO TOOLTIP. User has no idea what units this is in (mm) or what 0 vs 20 means. This is the most important snuggle setting after the enable checkbox.
- **Timeout (s)** input (line 6053-6061): NO TOOLTIP. User doesn't know what happens when timeout hits -- does it fail? Use best result so far?
- **Max parts** input (line 6064-6072): NO TOOLTIP. User doesn't know what happens above the limit -- fallback? Error?
- **Multi-plate overflow** checkbox (line 6074-6078): NO TOOLTIP. What does this do? If unchecked, do extra parts just vanish?

---

## 2. Progress Messages -- Developer-Speak

The progress callback in `SnuggleArrange.cpp` lines 279-282 produces strings like:

```
(Snuggle GPU gen 42/100 OK)
(Snuggle CPU gen 12/100 3 collisions)
```

Problems:
- **"gen"** means nothing to a normal user. Should say "step" or "pass."
- **"collisions"** is developer jargon. A user would understand "overlaps" or "parts still overlapping."
- **"OK"** is ambiguous. "All parts fit" or "no overlaps" would be clearer.
- The parentheses around the entire string look odd in a status bar.

---

## 3. Failure Notifications

**Verdict: adequate but inconsistent in tone.**

- When arrange is cancelled: "Arranging canceled." (line 635) -- fine.
- When some parts overflow: "Snuggle placed 5 of 8 parts. Some overflow to next plate." (line 655) -- good, clear.
- When all parts placed: "Snuggle complete -- 5 parts arranged. Bed usage: 42%" (line 666) -- good.
- When snuggle fully falls back (seq print, too many parts): "Arranging done." (line 672) -- **misleading**. The user enabled Snuggle, Snuggle did nothing, standard arrange ran instead, and the message doesn't mention this. Should say something like "Snuggle skipped -- standard arrange used."

The sequential printing warning in the UI panel (line 6099) is good -- amber text, clear message. But this is only visible if the panel is open. If the user closes the panel and hits Arrange, they get a silent fallback with no notification.

---

## 4. Slider Ranges and Defaults

| Control | Range | Default | Verdict |
|---------|-------|---------|---------|
| Part gap | 0.0 - 20.0 mm | 1.0 mm | **BUG**: UI allows 0.0 but backend clamps to 1.0 (SnuggleArrange.cpp line 217: `std::max(1.0f, params.snuggle_padding_mm)`). User sets 0, expects touching parts, gets 1mm gap. Slider minimum should be 1.0 or backend should honor 0. |
| Resolution | 0.5 - 5.0 mm | 2.0 mm | Fine |
| Population | 16 - 1024 | 64 | Fine |
| Generations | 10 - 500 | 30 | Fine |
| Timeout | 2.0 - 300.0 s | 5.0 s | Fine, though 300s = 5 minutes is a long wait with no progress bar |
| Max parts | 2 - 500 | 200 | Fine |
| Rotation | Locked/90/45/15/5/1 | 15 deg | Fine |

---

## 5. Compact Checkbox

Yes, there is a "Compact after arrange" checkbox at line 5964 with a tooltip. However:

**BUG: Duplicate compact checkbox.** There is a SECOND compact checkbox labeled "Post-GA compaction" at line 6080 inside the "Effort / Advanced" collapsible section. Both read and write `settings.snuggle_compact`. A user can:
1. Check "Compact after arrange" (visible by default)
2. Open "Effort / Advanced" and see "Post-GA compaction" already checked
3. Uncheck "Post-GA compaction"
4. Close the tree node
5. "Compact after arrange" is now unchecked too, but the user doesn't see the change

This is confusing. One of them should be removed. The top-level "Compact after arrange" is in the right place with the right label. The "Post-GA compaction" in the advanced section should be deleted -- "Post-GA" is developer jargon anyway.

---

## 6. Confusing Elements for Normal Users

1. **"Snuggle controls rotation -- see Lock Rotation below"** (line 5941): There is no "Lock Rotation" label below. The dropdown is labeled "Rotation" and one of its options is "Locked." The hint text references a control name that doesn't exist in the UI.

2. **"Genetic 3D nester -- voxel collision detection"** (line 6103): This description line at the bottom means nothing to a user. "Genetic," "nester," "voxel" -- pure developer vocabulary. Should say something like "3D shape-aware packing" or just be removed.

3. **"Effort / Advanced"** tree node label (line 6023): "Effort" is an unusual word choice. "Advanced settings" alone would be conventional.

4. **The button is labeled "Reset" (line 6117), not "Defaults."** The question asks about a "Defaults button" -- there is no such button. "Reset" is fine, but it resets ALL arrange settings, not just snuggle settings. A user who tweaked the standard spacing slider, then enabled snuggle, then clicks Reset loses their standard spacing value too. This may be intentional but is worth noting.

---

## 7. disabled_begin/end Mismatch Check

- Line 5946: `if (!settings_out.use_snuggle) { imgui->disabled_begin(true); }`
- Line 6092: `if (!settings_out.use_snuggle) { imgui->disabled_end(); }`

These are structurally matched -- both guarded by the same condition, within the same enclosing block (lines 5945-6093). No mismatch. The disable/enable will execute in pairs. **No lockup risk here.**

---

## 8. Defaults (Reset) Button Coverage

Lines 6130-6155 reset these snuggle settings:

| Setting | Reset value | Persisted? |
|---------|------------|-----------|
| use_snuggle | false | Yes |
| snuggle_lock_rotation | false | Yes |
| snuggle_padding_mm | 1.0 | Yes |
| snuggle_population | 64 | Yes |
| snuggle_generations | 30 | Yes |
| snuggle_voxel_mm | 2.0 | Yes |
| snuggle_compact | false | Yes |
| snuggle_max_parts | 200 | Yes |
| snuggle_timeout_s | 5.0 | Yes |
| snuggle_rotation_step | 15 | Yes |
| snuggle_multi_plate | true | Yes |

**All snuggle settings are covered by Reset.** Both in-memory values and persisted config keys are written. No gaps.

---

## Summary of Issues by Severity

### Must fix (user-visible bugs)
1. **Part gap slider allows 0mm but backend enforces 1mm minimum** -- silent lie to user (GLCanvas3D.cpp:5953 vs SnuggleArrange.cpp:217)
2. **Duplicate compact checkbox** -- two controls for one setting, confusing state sync (GLCanvas3D.cpp:5964 and 6080)

### Should fix (confusing UX)
3. **Missing tooltip on Part gap slider** -- most important slider has no explanation (GLCanvas3D.cpp:5950-5962)
4. **Missing tooltips on Timeout, Max parts, Multi-plate overflow** (GLCanvas3D.cpp:6053-6078)
5. **Progress messages use developer jargon** -- "gen", "collisions" (SnuggleArrange.cpp:279-282)
6. **"See Lock Rotation below"** references nonexistent UI label (GLCanvas3D.cpp:5941)
7. **Silent fallback message** -- "Arranging done." when Snuggle was enabled but skipped (ArrangeJob.cpp:672)

### Nice to fix (polish)
8. **"Genetic 3D nester -- voxel collision detection"** is developer-speak (GLCanvas3D.cpp:6103)
9. **"Post-GA compaction"** label uses jargon (GLCanvas3D.cpp:6080) -- moot if duplicate is removed
10. **"Effort / Advanced"** label is unconventional (GLCanvas3D.cpp:6023)
