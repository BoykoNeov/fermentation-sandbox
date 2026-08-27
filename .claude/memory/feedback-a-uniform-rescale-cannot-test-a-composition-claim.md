---
name: feedback-a-uniform-rescale-cannot-test-a-composition-claim
description: "A probe that scales one shared constant changes nothing RELATIVE — it cannot be evidence for a defect about composition, however well its number lands"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 739662e9-dbab-46bf-ab23-e6df88fbf32e
  modified: 2026-08-27T16:39:13.513Z
---

**Before reading a probe's gain as evidence for a mechanism, ask what it varied *relative to
what*.** A knob that multiplies a shared constant moves every term together. If the claim under
test is that some terms are mis-scaled *against each other*, that probe cannot speak to it —
whatever number comes out.

D-246 §6 (corrected at D-247). The availability gate reads `aa_i / (K·f_i + aa_i)`, `f_i` the
must-spectrum share. A per-species override makes the pools stop matching that spectrum, which is
a real commensurability defect. The record probed it by scaling `K_amino_acids` by 0.7155 — the
ratio of two pool masses — got propanol from 0.7963 to 0.8028, across a sourced floor, and wrote
that the defect "spans the whole of what is left of propanol's miss".

**It could not have.** The gate only ever reads the *product* `K·f_i`, so scaling `K` by `a` is
the same run as scaling **all eight** shares by `a` — measured identical to 12 digits. A uniform
scaling moves no share against any other, so nothing about composition was tested. The correction
the record *described* — each share re-referenced to the composition the run really holds, with
the spectrum sum preserved so the level cannot move — is worth **−0.000258**: 6.9 % of the gap, in
the wrong direction, stable across four decades of solver tolerance.

**Why:** the probe's ratio was built from a real pool mass over a *declared dose* the fixture only
carried to keep an isolability gate open. That folded the run's nitrogen level into a number
presented as being about composition, and the level is where all the movement was. A ratio of two
quantities is not automatically a test of the thing the ratio is named after.

**How to apply:** write the correction as an explicit per-term multiplier set and check its
weighted mean — a composition-only correction averages to 1 by construction, and a probe whose
multipliers are all the same value is a level change wearing a composition argument. Then pin the
identity (`scale the constant` ≡ `scale every share`) as an assertion, not a paragraph: it is what
stops the reading coming back. Compare tolerances at `rel=1e-12`, never bit-for-bit — the two runs
differ only in float associativity ([[feedback-assert-the-property-not-the-bit-pattern]]).
Related: [[feedback-a-normalisation-is-a-free-parameter]] (the reference of a ratio-form term is a
parameter nobody counts — and when nothing sources *either* candidate reference, no measurement
adjudicates and the answer is a sourcing ask), [[feedback-a-ratio-guard-cannot-see-a-common-factor]],
[[feedback-two-measured-quantities-do-not-locate-a-model-defect]].
