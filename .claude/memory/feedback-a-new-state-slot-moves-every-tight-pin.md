---
name: feedback-a-new-state-slot-moves-every-tight-pin
description: "Adding a state slot shifts BDF step selection with no model change; falsify a pin move by disabling the Process while keeping the slot"
metadata:
  node_type: memory
  type: feedback
---

D-248 added one carbon-park slot to the wine schema (97 → 98). Aging pins asserted at 1e-4
relative went red, and the tests' own failure text says "the burst leaked into the default build
— a finding for the decision record, not a tolerance to widen".

**Two causes had to be separated, and only one was the model.** Re-running with the new
Process's capacity set to **0** — inert, but the slot still in the schema — every entry sat
within **7.5e-5**, inside tolerance. The Process owned the 1.0–1.6e-4 that actually broke them.

**Why:** `solve_ivp`'s BDF error norm is RMS-weighted over the state vector, so `n → n+1`
changes the norm and therefore step selection, even for a component that is identically zero
with a zero derivative. The trajectory is mathematically unchanged and numerically is not. The
project already recorded this for the `quinone` slot at `test_oxidative_cascade_guards._PIN_RTOL`;
D-248 is the second instance, so it is a pattern rather than a one-off.

**How to apply:** when a beat that adds a slot breaks a tight pin, **falsify the model reading
first** by neutering the new Process and keeping the slot — that separates "my schema moved the
solver" from "my physics moved the answer", and it is what licenses re-pinning in the D-182
old → new style instead of widening a tolerance. Only re-pin the entries that actually broke.
Related: [[feedback-a-screen-is-not-idempotent-under-its-own-repair]].
