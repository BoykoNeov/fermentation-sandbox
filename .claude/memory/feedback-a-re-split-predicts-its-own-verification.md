---
name: feedback-a-re-split-predicts-its-own-verification
description: "A repair that only re-splits a conserved quantity across two slots predicts bit-identity of everything else; pre-register that and the diff finds your unenumerated consumer"
metadata:
  node_type: memory
  type: feedback
---

D-250 moved nitrogen from one slot to two. The two share a source and every sink, so the **sum**
had to follow exactly the trajectory the single slot followed before. That gives a prediction
strong enough to be the verification: **every observable except the one the repair targets must be
bit-identical.**

Measured, one process per tree: **beer 0.00e+00 on everything**, wine ≤**1.9e-7** on biomass,
sugar, ethanol, all eight amino-acid pools, three fusels, H2S and the carbon park. Only pH moved.

**It is diagnostic, not just reassuring.** Anything else moving means exactly one of two things:
a consumer of the old slot was not enumerated and is now reading a collapsed pool, or a paired
flux (here the D-32 swap's refund of growth's draw) is not split on the same rule as its partner.
Both are real hazards — the unsplit refund sends net `dN/dt` positive and re-creates the very
artefact the repair removes, with the repaired Process innocent.

**Write the residual off to the schema only by falsifying it.** The wine's 1e-7 is BDF step
selection reacting to 98 → 99 slots, and the proof is the **undosed** arm, where the Process is
disabled at compile and cannot contribute — it moves 1.17e-7 too.

**How to apply:** whenever a change conserves a quantity and only redistributes it, pre-register
bit-identity of everything downstream **before** editing, then diff across trees. It costs one
probe, it is the cheapest diagnostic available, and it is what turns "I think I enumerated the
readers" into a check. Related: [[feedback-a-new-state-slot-moves-every-tight-pin]],
[[feedback-the-setting-where-a-change-is-exact-is-the-control]], [[feedback-pre-register-the-cheap-prediction]].
