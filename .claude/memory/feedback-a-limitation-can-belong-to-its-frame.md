---
name: feedback-a-limitation-can-belong-to-its-frame
description: A limitation blamed on the mechanism can be a property of the scenario it was measured in — re-measure it on the frame that has data before building anything
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 870e1cf7-7f86-43eb-b8df-19e2af8180be
  modified: 2026-08-26T07:39:20.560Z
---

Three surfaces (two decision records and the memory ledger) named one inherited limitation as
the next thing to build: "growth stops dead at day 0.92 with 81 % of the sugar unfermented, so
all the aroma is made on day one." It was **true, and it was a property of the test scenario's
inoculum** — that frame pitched a flat 1.0 g/L, the same back-computed residual an earlier beat
had already retired and removed from the *other* scenario in the same medium. On the one trial
with measured cell counts, the model's growth curve was already inside the measured spread at
every informative day. Correcting the frame widened the window 1.53× and moved **zero** of the
eight calibrated levels (6e-7 relative), so nothing in the suite could see it either way.

**Why:** a limitation gets attributed to the mechanism because that is where it was noticed, and
the attribution then rides forward in every copy-forward "still open" list. A scenario input is
not part of the mechanism, is rarely re-derived, and can carry a number that was retired
elsewhere in the same repo — so the same defect survives in the frame nobody re-read.

**How to apply:** before building a mechanism to fix a named limitation, re-measure the
limitation on the frame that actually has data, and separately check the complaining frame's own
inputs are sourced rather than conventional. If the two frames disagree, the limitation is the
frame's. Two corollaries: sample the model the way the trial sampled it (a continuous `argmin`/
`argmax` against a daily-sampled measurement invents a gap — here 2.32 d vs day 3, which agreed
once both were read daily, cf. [[feedback-read-a-fast-curve-on-a-fixed-grid]] on the MEASURED
side); and a correction every existing guard is blind to is exactly the one that owes a new guard
[[feedback-prefer-the-variant-your-guards-can-see]]. Related: [[feedback-agreement-can-be-a-frame-difference]],
[[feedback-check-the-blocker-is-still-blocking]].
