---
name: a-positive-schema-check-is-not-a-licence
description: "A slot existing proves the term is expressible, not that it is buildable - where the slot is READ decides its sign, and a green schema check reads as a go-ahead"
metadata:
  node_type: memory
  type: feedback
---

**A positive schema check is not a licence to build.** Opening a beat by confirming the quantity
is representable answers *can the model say this*, never *what happens when it does*. At D-214
`peptide_buffer` was a real state slot, so trub settling looked expressible and the check was
recorded as a premise. It was the trap: that pool is the **t=0 cation back-solve's counter-anion**,
so a Process draining it during fermentation does not model less buffering — it unbalances the
charge equation. A 20 % cut at 6 h took beer to **pH 7.08**; a full one to **11.66**.

**Why:** the archive's usual failure is the *negative* schema claim — "the model can't represent X"
inferred from one Process not reading X ([[feedback-check-the-schema-not-the-caller]]). The
positive direction has the opposite asymmetry and no guard: a slot that exists says nothing about
**where in the pipeline it is written**, and a quantity fixed by an anchor behaves nothing like the
same quantity evolving freely. D-205's *expressible ≠ identifiable* is the same lesson one level
down; this is *expressible ≠ buildable*. Both times the cheap check came back green and the
expensive one refused.

**How to apply:** after confirming a slot exists, ask **who writes it and when** before designing
anything on it. Specifically: is it seeded by a compile-time solve (an anchor, a back-solve, a
calibration)? If so, a runtime term touching it is fighting that solve, not extending it, and the
sign is likely to invert. Say which side of the anchor the mechanism belongs on **in the
pre-registration**, and if the answer is "before", the honest target is usually the *parameter*,
not a new Process. Cheapest discriminator: cut the pool mid-run and look at the sign — one probe,
and it either matches the mechanism you meant or it exposes that you are modelling something else.
Relatedly, when a mechanism turns out to sit before t=0, check whether the calibration already
absorbs it before calling it an omission — D-214's trub was inside Peyer's control-wort anchor all
along.
