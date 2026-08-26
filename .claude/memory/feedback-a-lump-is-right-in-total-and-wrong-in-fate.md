---
name: a-lump-is-right-in-total-and-wrong-in-fate
description: "A constant back-solved to reproduce a total absorbs every unmodelled contributor to that total - including ones whose FATE differs from the lump's"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 80611b07-f4ed-492d-8b4b-03c3866a0e21
  modified: 2026-08-26T18:41:00.312Z
---

A constant back-solved so the model reproduces a **published total** silently absorbs every
contributor to that total the model does not carry. The fit is then right at the calibration
point and wrong everywhere the absorbed contributors have a **different fate** from the lump.

`peptide_buffer_capacity_beer` is rooted so beer's wort reproduces Peyer's BC = 1.18. The
back-solve ran on eight organic acids and the lump and nothing else, so a real wort's free
amino-acid buffering had nowhere to go but into the lump. The lump is permanent; the amino acids
are eaten inside two days. Beer's wort therefore buffered correctly at t=0 and went on buffering
like a wort for the rest of the ferment (D-239). Nothing about the shipped 1.5481 g/L was
"wrong" — it conflated two pools with two fates.

**Why:** a calibration constrains the model at ONE state. Fitting it hides a structural error
whenever the absorbed share is dynamic and the absorber is not, and the error grows with time
rather than showing up as a bad fit.

**How to apply:** for every back-solved lump, ask *what else is really in the quantity I fitted
to, and does it leave?* Where a contributor is sourced and has a different fate, the fix is a
**re-partition at constant total** — split it out, hold the calibration exactly, and price the
trajectory difference — not a re-fit. Expect the split to be invisible at the anchor and to grow
afterwards; and check the sign at BOTH ends, because while the drained pool is still partly
present the split can move the model the OTHER way
[[feedback-a-calibrated-level-decays-when-anything-upstream-moves]],
[[feedback-a-derived-yield-encodes-its-rate-law]].
