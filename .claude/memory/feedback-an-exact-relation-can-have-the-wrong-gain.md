---
name: an-exact-relation-can-have-the-wrong-gain
description: "Before reporting an estimate from an exact mass balance, compute how precisely its inputs must be known to discriminate - an exact lever with the wrong gain is a null result"
metadata:
  node_type: memory
  type: feedback
---

`X0 = X_final - YAN/f_N` is an exact identity and it cannot settle the inoculum. The
inoculum is a few per cent of final biomass at either candidate, so the entire candidate range moves
the implied yield by **8.0 %**, and separating the two candidates needs `f_N` known to **6.7 %**
where its own declared band spans a factor of **1.75**. Off by more than an order of magnitude. The
point estimate was frame-dependent too, and the two frames landed on opposite sides of the value
being tested.

**Why:** exactness is about the algebra, not the discrimination. A relation can be true, computable
and completely unable to tell two hypotheses apart - and it will still print a confident number.

**How to apply:** before quoting the estimate, compute the sensitivity: how far apart are the
candidates in the OUTPUT, and how well is each input known? If the required precision is inside the
inputs' bands, the estimate discriminates; if not, ship the null with the required gain as the
finding, and pin the gain rather than the conclusion so the null reopens if a band ever narrows.
Related: [[feedback-a-band-is-per-parameter-a-claim-is-joint]], [[feedback-an-always-agreeing-landmark-is-not-evidence]].
