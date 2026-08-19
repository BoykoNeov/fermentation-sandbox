---
name: feedback-a-locked-pair-repairs-or-drifts-together
description: "Two quantities sharing a rate law are locked exactly, so one's history transfers with no archaeology - and a repair that skips the other leaves a guaranteed defect"
metadata:
  type: feedback
---

Beer's three esters share one rate law, one activation energy and one stripping constant, so
their ratios are fixed by their `k` alone — measured invariant **to 7 significant figures**
across both speed knobs' band edges and a 15 °C run that moved every pool 2.7×. Two things
follow, and I used only the first at D-225 before noticing the second.

**Forwards, it is free evidence.** A `git log -S` pickaxe over the archive timed out twice; the
lock made it unnecessary. The repaired pool's already-measured history transferred by one ratio
and dated the drift exactly — 0.2326 mg/L at D-99, 0.1731 after D-223 — with **no old tree
checked out**. When a quantity you need the history of is locked to one whose history is already
in the record, the archaeology is arithmetic.

**Backwards, it is a liability.** D-224 diagnosed the drift mechanism correctly and repaired
**seven of the eight** constants that mechanism condemns, enumerating them from a hand-written
list. Because the eighth was locked to one of the seven, its defect was not a risk — it was
**guaranteed**, at exactly the same factor, and it will recur on the next speed change. A locked
set has no partial repair: you fix all of it or the ones you skipped keep the whole error.

**Why:** a shared rate law is a shared fate. Reasoning about the members one at a time hides
both the cheap inference and the certain defect, and a guard built from a list of members cannot
tell you which ones you forgot — see [[feedback-grep-finds-claims-not-guards]].

**How to apply:** when two quantities look like they might share a rate law, **test the lock**
before doing anything else — vary the upstream knobs and a temperature and check the ratio to
several significant figures. If it holds: transfer the history instead of digging for it, and
treat the set as **atomic** for any repair, enumerating it from the code registry rather than
from a list ([[feedback-a-magnitude-defence-can-be-scoped-to-the-smallest-member]] is the same
error with a different denominator). If the ratio moves, you have a coupling worth its own line.
Related: [[feedback-a-calibrated-level-decays-when-anything-upstream-moves]],
[[feedback-gate-both-halves-of-a-pair]].
