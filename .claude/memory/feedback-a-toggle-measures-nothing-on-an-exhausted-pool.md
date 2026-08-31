---
name: feedback-a-toggle-measures-nothing-on-an-exhausted-pool
description: "Switching a route off and diffing the state returns ~0 when its precursor is supply-limited - which reads as 'the route does nothing' rather than 'this method cannot see it'"
metadata:
  node_type: memory
  type: feedback
---

The exact-state-difference counterfactual — turn a route off, diff the end state, call the
difference that route's contribution — **silently returns zero whenever the pool it draws from is
fully consumed either way.** The answer looks like a finding ("the route contributes nothing")
when it is actually the method failing to have any signal.

**Why:** at D-254 I tried to measure the valine→isoamyl branch by zeroing its share parameter and
diffing. First attempt diffed **isoamyl** and got ~1e-14, because the route sets carbon
*provenance* and the rate law sets the *amount* — wrong slot entirely. Second attempt diffed
**valine**, the right slot, and got ~1e-14 again — because valine exhausts in both arms, so its
total consumption is fixed by supply and only the *fate split* moved. That split lives in no state
slot at all. Both zeros were convincing and both meant "ask a different way". D-112 had already
recorded the underlying fact ("where a precursor is fully consumed the total drawn is fixed by
SUPPLY"), which is exactly the sentence that should have stopped attempt two before it ran.

**How to apply:** before trusting a toggle diff, ask two questions. *Which slot would actually
move?* — a provenance change and an amount change land in different places, and for a re-routing
Process it is almost never the product. *Does the source pool exhaust?* — if it ends at ~0 in both
arms, the diff is structurally zero and the method is blind, so reach for integrating the applied
draw along the trajectory instead, with a grid-refinement check and a cross-check of the
quadrature against some exact state difference you *can* form. And treat a suspiciously exact zero
as a signal about the instrument, not a result: two of them in a row is the instrument.

Related: [[feedback-a-ratio-guard-can-pass-on-an-overproduction]],
[[feedback-a-clamp-sits-between-designed-and-realised]].
