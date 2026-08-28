---
name: feedback-read-a-channels-timing-against-its-own-run
description: "A sub-process that looks \"1.5x fast\" against a paper's landmark may be slower than the whole run containing it — compare on the same trajectory and divide the clock out before proposing a fix"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6025169e-c573-4f1a-af77-db80fef3cfd0
  modified: 2026-08-28T07:44:44.733Z
---

A record left the nitrogen exhaustion time "open and unattributed" at **~1.5× fast** against the
paper's own landmark, and the owner picked it as a beat. Read against the *fermentation the
nitrogen sits inside* — same trajectory, only the observable differing — the whole run was
**1.92×** fast where the nitrogen channel was 1.59×. The channel was **slower than the run
containing it**, so it could not be what made the run fast. Dividing the clock out entirely
(each landmark as a fraction of that run's own duration) the model exhausted its nitrogen at
0.225 of its run against the paper's 0.187 — *later*, not earlier. The whole-run gap turned out
to be an already-characterised, already-guarded cross-study difference from 190 records back.

**Why:** a landmark comparison has a single number on each side, so the sub-process gets the
whole discrepancy attributed to it by default. Nothing in "1.5× fast against her 28 h" names a
denominator, and the run containing the channel is usually not measured in the same breath — so
a global effect gets booked against whichever local channel was under the microscope that week.

**How to apply:** before pricing a repair to any *timing*, measure the enclosing process on the
same trajectory and compare the two gaps. Then compute the time-free version — the landmark as a
fraction of the run's own duration — because that survives multiplying either time axis by any
constant and cannot be argued away by a rate calibration. Check whether the enclosing gap is
already characterised somewhere in the repo before treating it as new; and prefer a
pitch/frame-invariant statement, since a claim that holds at two frames does not depend on a
constant nobody sourced [[feedback-a-limitation-can-belong-to-its-frame]]. Corollary on the
other side: a residue that is genuinely local will be **time-free** (here, biomass-against-
nitrogen), and that is the one worth building on. Related:
[[feedback-compute-the-clean-fix-before-adopting-it]],
[[feedback-an-expectation-transfers-with-the-channel-not-the-subject]].
