---
name: feedback-a-control-belongs-where-the-number-is-made
description: "A validity control asserted in one consumer does not protect a shared fixture's other consumers — put it where the number is computed, or a mutation passes most of your guards"
metadata:
  node_type: memory
  type: feedback
---

When a fixture computes a number whose validity depends on a control, assert that control **inside
the fixture**, not in one of the tests that reads it. A control living in a consumer protects that
consumer and nothing else.

**Why:** D-259's fixture returned a realised split computed as a ratio of two integrals recomputed
along the trajectory. Both integrands are *formulae*, so the fixture will happily produce a split
even when the Process that actually ran was not the counterfactual one — and that split means
nothing. The quadrature/closure check that ties the ratio to the run was written as an assert in
the one test about the bracket. Mutation arm A (swap the counterfactual back for the shipped sink)
therefore **passed three of the four guards**: they read the same poisoned fixture and never
checked it. Moved into the fixture, the same arm errors all four at setup.

The tell is structural, not stylistic: **the number and its validity condition were computed in the
same function and asserted in different ones.** Anything that can produce a plausible value from
the wrong run needs its check co-located with the computation.

**How to apply:** when a fixture derives a quantity that only means something under a condition
(a closure, a conservation residual, a convergence check, "the arm I think ran actually ran"),
assert the condition where the quantity is built. A consumer may *also* assert it for visibility —
that is fine and costs nothing. Then run the mutation that invalidates the condition and confirm
**every** consumer goes red, not just the one you wrote the check in; a mutation that reds only
some of a shared fixture's guards is telling you the control is in the wrong place.

Related: [[feedback-a-shared-fixture-has-two-consumers]],
[[feedback-a-non-vacuity-check-can-itself-be-vacuous]],
[[feedback-mutate-the-premise-before-building-the-guard]],
[[feedback-check-the-published-test-can-fail]].
