---
name: feedback-a-guards-teeth-are-in-its-literals
description: A circularity objection names an ACT, not a tool — check which literals a guard holds against each other before believing a fix would hollow it out
metadata:
  type: feedback
---

D-233 §8 declined to make beer's peptide capacity/pKa pair coherent per member, on the ground
that moving the BC back-solve into `src` "would make the round-trip test compare the root-finder
against itself". That reason was carried forward unexamined by D-234, D-235 and D-236, and the
repair sat open for four records.

It was aimed at one act — **deriving the shipped constant**. The repair that was actually needed
re-roots the per-member *seed* and leaves the YAML literal alone. The guard's teeth were never in
owning a copy of the arithmetic; they were in holding **two literals** against each other (the
shipped constant, and Peyer's published 1.18), and neither is produced by the solver. Both failure
modes it exists for still fire after the change.

**Why:** "that fix would make the test circular" sounds decisive and is cheap to repeat. It is a
claim about *which value the test derives*, and it is checkable in one minute: name the two things
the assertion compares and ask whether the proposed change produces either one.

**How to apply:** Before accepting an archived objection as a blocker, re-derive it. Write down the
guard's anchors. A change that touches neither anchor cannot hollow it out, whatever it moves
elsewhere. And when you record such an objection yourself, name the **act** it forbids, not the
tool — "deriving the shipped constant in `src`", never "putting a root-find in `src`".

Related: [[feedback-verify-latest-state-not-breadcrumbs]],
[[feedback-a-defect-pin-can-outlive-its-defect-by-driving-another-path]].
