---
name: feedback-a-shape-change-is-a-change-to-the-pair
description: "Changing one half of a source/sink pair changes the invariant the PAIR maintains, not just the term you edited (D-189)"
metadata:
  node_type: memory
  type: feedback
---

When a quantity is maintained by a **pair** of opposing terms, the unit of a shape change is the
**property the pair conserves**, never the term being edited. At D-189 the item had been framed
for 82 records entirely in terms of the source — "growth-coupled *excretion*", deferred because
"nothing reads the peak", diagnosed as "the *excretion* SHAPE". Not one framing mentioned the
sink. But the thing the pool exists to deliver — a persistent residual — is not produced by the
source at all: it is what the pair leaves behind, and it survives only because **both halves die
at the same moment** (dryness). Moving one driver made the residual stop being a ratio and become
an exponential of the window that opened between the two deaths: `exp(-k · ∫flux)` ≈ 2.4e-5, so
2.000 mg/L → 0.000 and the downstream aroma with it.

**Why:** the edit looks local and is not. Every framing that named only the source made the change
sound scoped, and each one was written by someone who understood the model. The measured price
appears only when you ask what the pair *guarantees* — here a residual pegged to end-of-fermentation
and independent of the ferment — and check whether the edit still delivers it. A second signature:
the rescue that restores the number may be **output-identical** and still wrong, because it buys
the number by asserting something unsourced about the other half.

**How to apply:** before editing one term of a balance, write down the invariant the pair holds and
what pins it (which two events coincide, which ratio survives). Then measure that invariant, not
the edited term. Check whether the parameter file already anticipates the failure — at D-189 the
provenance note said the freeze existed *"rather than draining it to ~0 over the long tail"*, which
is exactly what happened, and nobody had connected that sentence to the proposal. And price the
invariant's robustness across a condition sweep, not at one operating point: the shipped pair held
the residual to 5.65e-7 across 15–28 °C, both rescues to 2.03× and 16.33×.

Sibling: [[feedback-mutate-the-premise-before-building-the-guard]] — at D-189 mutating first also
showed **five existing asserts** already forbade the change, so no new guard was owed.
