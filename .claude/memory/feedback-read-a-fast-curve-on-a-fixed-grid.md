---
name: read-a-fast-curve-on-a-fixed-grid
description: "argmin over solve_ivp's adaptive output samples wherever the solver stepped, not the time you asked for - it put a headline number 26% off and reproduced it exactly for two records"
metadata:
  node_type: memory
  type: feedback
---

**`argmin |t − t_target|` over an adaptive trajectory reads the curve wherever the SOLVER stepped.**
D-211 §9 reported beer's day-1 pH sitting **0.034** above its ceiling. D-214 re-ran it: the
adaptive read landed on **t = 23.6382 h** — 22 minutes early, on the steepest part of the pH
course — and the correct fixed-grid value is **0.0274**, a 26 % difference in a headline number
that two records then quoted. The day-7 figure moved only 0.0004, because that part of the curve
is flat: **the error scales with the local slope**, so it hides on every quantity except the one
you are usually measuring.

**Why:** this is not solver noise and no tolerance protects against it
([[feedback-pin-tolerance-vs-solver-tolerance]] is about a different failure). Passing `t_eval`, or
interpolating with `np.interp`, samples the time you asked for; omitting it samples the time BDF
found convenient, which shifts whenever anything perturbs stepping. Worse, it is **stable per
configuration** — it reproduced 0.034 and 0.0082 to the digit — so it looks like a settled
measurement rather than a read artefact, and a later beat that fixes the grid appears to be
contradicting the model instead of the probe.

**How to apply:** always give a trajectory a fixed grid (`t_eval=`) or interpolate onto the exact
time before reading a value out of it, in probes as well as tests — the shipped suite here already
does the right thing via `np.interp`, and the damage was done entirely in throwaway probe scripts
that nobody diffed. When a new measurement disagrees with a recorded one by a few thousandths,
**suspect the read before the model**: run the same scenario both ways and print the timestamp the
adaptive read actually landed on. That check is two lines and it converts "my numbers disagree with
the archive" into a named correction. A control arm comparing one model to itself along two paths
is what surfaces this at all — D-214 only noticed because a null arm that should have printed 0.00
printed 1.04e-2.
