---
name: feedback-a-containment-claim-is-scoped-to-where-the-ranges-nest
description: "\"Band A contains band B, so B's uncertainty is already reported\" is measured at ONE evaluation point; if A is fixed literals and B slides with the scenario, they nest only over an interval — and the verdict gets restated as a property of the pair"
metadata:
  node_type: memory
  type: feedback
---

A "this band subsumes that one" argument is the strongest reason to *decline* to sample a
parameter, so it deserves the scrutiny of a decision rather than of a check. Its weak spot is that
it is almost always evaluated at one scenario, then written down as a fact about the two bands.

**Why:** D-241 declined to draw two regression coefficients because the parameter they derive is
itself sampled over a band that contained their implied range and was 2.11x wider. That was
measured honestly, at nine corners — **at the battery wine's YAN of 250 mg/L**. The verdict then
entered the census registry with no qualifier at all. But the containing band is two hard-coded
literals while the implied range slides with the evaluation point, so the nesting holds only over
`YAN ∈ [66.0, 324.8] mg N/L`, and the width ratio runs **7.16x at YAN=50 and 1.22x at 350** — the
"2.11x" is 250's value, not the pair's. The repo already runs a scenario at 50. The escapes were
small and the verdict survived in substance, so this **sharpened** rather than overturned it — but
only because someone swept the axis.

**How to apply:** whenever you write "A contains B", finish the sentence with *"for X in [..]"*. If
one side is constant and the other is a function of anything scenario-specific, the claim is an
interval, so solve for the crossings rather than sampling a point — `brentq` on
`implied_edge(x) − bracket_edge` costs one line. Guard **both** directions: assert containment
inside the interval *and* its failure just outside, or the guard passes equally well on an interval
claimed to be the whole real line [[feedback-a-non-vacuity-check-can-itself-be-vacuous]]. And when
a red arrives, do not widen the containing band to restore it — that band is what the ensemble
actually draws, so widening it changes the reported spread of the very quantity under discussion.
