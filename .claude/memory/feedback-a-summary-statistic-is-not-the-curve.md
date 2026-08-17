---
name: feedback-a-summary-statistic-is-not-the-curve
description: "An endpoint/aggregate metric can PASS while the trajectory behind it is badly wrong; and an OUT call needs a measured tolerance, not a chosen one"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 09a91935-982e-429a-ba2b-c094a06612d5
  modified: 2026-08-17T08:35:06.040Z
---

A validation that scores one **summary number** is not a validation of the **curve**, and the
gap can be enormous while the test stays green. D-207: beer's pH agreement was a single
fraction-of-the-total-drop (87.1 %, inside its own 77.6-97.0 % band) for four beats. Reading
the daily course from a figure **already rendered to disk** showed the model completes
**85.2 %** of its pH fall by day 1 where the measurement has done **42 %** — 0.195 pH out,
**8.1×** the read tolerance — and that it is BELOW the measured band early and ABOVE it late.
It overshoots, then stalls. The endpoint cannot see any of that, by construction.

**Why:** an aggregate is a projection; a wrong trajectory and a right one map to the same
scalar. So "the test passes" is evidence only about the projection, and the acceptance width
compounds it — a band wide enough to pass is a band wide enough to hide a shape. The tell is a
test whose *docstring* says "trajectory" while its *assert* reads two endpoints.

**How to apply:**
- Before trusting a validation, ask **what shape of quantity it asserts** — a point, a
  difference, or a curve. If the source has a curve and the test has a scalar, the curve is
  unread work, not settled work. Check the source's other panels: D-180 transcribed one figure
  day-by-day and took two numbers off the panel beside it.
- **An IN/OUT call needs a tolerance you measured.** My first draft claimed "outside on five of
  seven days" using an unjustified ±0.03 slack. Measured, the tolerance is **0.024 pH** — set
  not by extraction precision (0.0028) but by how far *two independent reads of the same
  figure* disagree. That re-scored four of the five days to ~2×, leaving **one** day carrying
  the finding. Quantify before weighting, or the record overclaims and a later beat walks it back.
- **A frame difference is a fork, not a caveat.** Whether the measurement was taken under the
  same conditions the model reports (here: degassed sample or not) can INVERT the diagnosis —
  the two frames gave "falls too fast" and "falls too slow" with different culprits. Resolve it
  or report both; never pick the one with the cleaner headline.
- Pair with [[feedback-re-read-the-source-you-already-mined]] (this is its second instance) and
  [[feedback-conceded-caveats-are-not-coverage]]. When the mutation comes back green, weigh
  [[feedback-mutate-the-premise-before-building-the-guard]] against the "fights the fix" test: a
  pin on the defect's own shape is one a correct fix must delete, so ship the measurement and
  the DATA instead, and say in the record that no guard was built and why.
