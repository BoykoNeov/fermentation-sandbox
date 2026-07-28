---
name: feedback-a-null-result-needs-a-positive-control
description: "\"Nothing moved\" is not evidence of unreachability — run a control through the identical harness, or the null confirms itself"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ca267d52-61b6-45bb-bbf9-09b1270dbd86
  modified: 2026-07-28T19:32:18.857Z
---

When a test's claim is a **null** — "perturbing this changes nothing", "no
call sites exist", "the sweep found none" — that null is worthless until a
**positive control runs through the identical harness at identical settings**
and comes back non-null. Design the control before trusting the null.

**Why:** D-159 planned to pin structural undrawability as "force the parameter
into the sample with `only=`, assert the trajectories are identical". The
exemplars all came back `max|dy| == 0.0`. So did **five genuinely
Process-read parameters** — `k_copper_multiplier` (the `copper` slot is never
written), `ethanol_tolerance` (its term is clamped off below tolerance), and
`k_so2_oxidation` / `k_browning_base` / `k_ethanol_oxidation` (supply-limited:
no O₂ in that scenario). "Frozen" conflates *unreachable* with *reachable but
zero in this scenario*, so the test would have passed for the wrong reason on
every supply-limited parameter and shipped as a decoration. Only the control
exposed it — the first two controls picked were themselves frozen, which is
how the whole problem surfaced.

**How to apply:** Pick the control so it *cannot* be null for an incidental
reason — prefer a quantity active from `t0` over one gated behind a supply,
a dose, or a slot something else must write. If the control is also null, the
harness is broken; do not interpret the exemplars. Then strengthen the null
into a **pair**: *consequential* (the value demonstrably reaches the model, at
a named site) **and** *inert under the lever you are testing*. Prefer the
consequential half be **mechanism-shaped** — an identity against the slot or
the write it lands in — since a mechanism assertion cannot be satisfied by a
term that happens to be zero, and an integration-shaped one can.

Same family as [[feedback-count-and-print-your-skips]] (a silent denominator
reporting "5 of 5 clean"), [[feedback-verify-the-restore-between-mutation-arms]]
(always design one arm to be GREEN) and
[[feedback-mutate-the-premise-before-building-the-guard]]. All four are one
rule: **an experiment that cannot produce the other answer has not measured
anything.**

Its sibling trap is **vacuity by shared derivation**: if both sides of an
assertion come from the same declaration, the test is set arithmetic, not a
fact about the code. D-159 refused a population test for exactly this —
`structural = banded − declared_reads` versus `undrawn = banded −
_schedule_reads(...)` read the same `reads` tuples. See
[[feedback-check-the-schema-not-the-caller]].
