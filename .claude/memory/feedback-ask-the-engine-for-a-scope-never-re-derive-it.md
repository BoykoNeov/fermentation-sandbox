---
name: feedback-ask-the-engine-for-a-scope-never-re-derive-it
description: "A scope the engine already resolves must be ASKED FOR, not re-derived from declared reads — the walk misses schedule reads and seed reads (12 of 97), and the raw scope then over-counts pinned draws"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 86025c03-b1ae-45b8-b84a-18332d969846
  modified: 2026-09-01T12:15:29.945Z
---

**A count you print as a promise has to come from the code that will keep it.** Re-deriving a
sampling scope by walking the active mechanisms and unioning their declared `reads` looks
complete and is not: the ensemble also samples what the *schedule* reads and what the compile
seam read to build `y0` (D-241's seed reads), and neither is visible from a `reads` walk by
construction.

D-261. The console prints "ranking which number drives the uncertainty needs more re-runs than
there are varying numbers, so at least N here" *before* the user spends a minute of compute.
The `reads`-walk projection came back **85 where the engine sampled 97** on a wine and **81 where
it sampled 89** on a beer — always low, never high, so anyone following the printed figure got
the underdetermined error anyway, after the wait. Replacing it with a call to the engine's own
`_resolve_sample_names` (one compile, no integration) fixed the direction — and immediately
exposed the *opposite* error: the raw scope includes parameters pinned to a single value, which
the attribution drops for having zero variance and which therefore do not count against the
budget. The honest number is **the engine's scope minus what cannot move**: 89 for the wine,
which is exactly the figure the engine's own underdetermined error names.

**Why:** both errors are silent and neither is visible on screen. Under-count and the promise
fails after the cost is paid; over-count and the user waits longer than needed, having just been
"corrected". A re-derivation also rots on a schedule nobody controls — every future channel the
engine folds into its scope is one the walk will miss, with no signal.

**How to apply:** when a UI or a report states a quantity the engine computes, call the engine's
resolver even if it is private and even if it costs a compile — an import of a private helper in
the same repo is cheaper than a copy that drifts (same reason `app/library.py` reads the allowed
initial keys and the verb tables from `compile.py` rather than listing them). Then pin it: run a
tiny 4-member ensemble in a test and assert the projection equals what was actually drawn, and
walk the boundary end-to-end at a deliberately narrowed scope — N re-runs refused, N+1 accepted.
Related: [[feedback-an-empty-combine-returns-the-top-mark]],
[[feedback-verify-latest-state-not-breadcrumbs]].
