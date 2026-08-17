---
name: feedback-two-measurements-one-figure-can-disagree
description: Two series read off the same paper can be mutually inconsistent in your frame; check the target against the OTHER observable before building toward it
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 319bf67b-4e83-4c11-9ffb-068bbcb68911
  modified: 2026-08-17T15:39:46.151Z
---

A beat was picked because a measured target sat on disk: beer's acetic course, transcribed from
Tyrell Fig. 13, against which the model was short **64.87 mg/L at day 1**. The obvious build was
"close that gap".

Solving instead for what the **pH** admits — the *other* series, Fig. 4 of the same paper, same
four ferments — gave a **window**: acetic 94.24-141.22 / 90.20-136.82 / 86.11-132.38 mg/L at the
three arms of the nitrogen-charge band. **The measured 145.0 sits ABOVE that window at every
arm.** The model needed only **+6 to +14 mg/L**, a fifth of the apparent gap. Closing the gap
fully would have driven day-1 pH *below* the day-1 pH measured from the same figure — a build
that scored as a success on one observable while breaking the other.

**Why:** a target series constrains your parameter only *through your model's frame*. Two series
from one paper are consistent in the author's world, not automatically in yours — every term you
have not built sits between them. So "the model is 64.87 short of the measurement" is a statement
about the *pair*, and the size of the real defect is set by whichever observable is tighter, not
by the one you happened to be looking at. This is [[feedback-a-hit-can-be-two-errors-cancelling]]
caught *before* the build instead of after, and the same shape as
[[feedback-fit-the-observable-not-the-consequence]].

**How to apply:** before building toward a transcribed target, **solve for the range of that
quantity which every OTHER measured observable admits**, at all band arms, and check the target
is inside it. If it is outside, say so first — the gap you were going to close is not the defect.
Two cheap habits that made this work: the check is a *window* (solve for both band edges), not a
point; and confirm the discrepancy is specific — here every other charged acid was inside its own
measured band, which is what ruled out "the day-1 profile is generally low" and made the pair
claim assertable. Cf. [[feedback-a-summary-statistic-is-not-the-curve]],
[[feedback-pin-the-band-not-the-nominal]].
