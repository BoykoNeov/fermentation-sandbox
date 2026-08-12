---
name: feedback-a-threshold-cannot-separate-same-sign-regimes
description: A one-sided magnitude assert passes for the wrong reason when the broken premise moves the quantity the SAME way as the healthy one — and a green mutation is what OWES a guard
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 53cf9b2e-d59e-43b1-8918-019200d03f50
  modified: 2026-08-12T10:05:09.222Z
---

`assert residual_sugar > 50.0` was meant to witness *"the wine is sweet"*. When D-192's brake
silently turned the scenario into a must that never finished fermenting, the assert **still
passed** — because an unfermented must has *more* sugar than a stuck one. All 98 tests in the file
stayed green across a total change of regime.

**Why:** a one-sided threshold on an unqualified pool only sees magnitude. If the failure mode
pushes the quantity in the **same direction** as health does, the assert cannot separate them, and
its name goes on telling you it did. Assert counts don't help either: D-195's pyruvate had 27
asserts naming it and D-189's α-ketobutyrate had none, and *both* mutations went red without a new
guard, while D-194's — with eleven tests supposedly about the thing — came back green.

**How to apply:** when an assert is the witness for a *premise* (this is a sweet wine; this ferment
finished; this run reached steady state), ask what the broken premise does to the number's **sign**,
not its size. If broken and healthy agree in sign, add a predicate the two regimes disagree on —
here, "arrived" (`E` at the ethanol ceiling at the breakpoint) and "stuck" (relative drift across
the tail), which read 156.05/0.0028 healthy against 3.93/0.1238 broken.

And treat the mutation's result as the **licensing condition**, not a formality:
[[feedback-mutate-the-premise-before-building-the-guard]] says break it first — the addition here
is that **GREEN is the answer that owes a guard**, and it is not predictable in advance. Pair it
with [[feedback-verify-the-restore-between-mutation-arms]]: keep the broken value as a permanent
positive control inside the new test, so a threshold discriminating nothing can never read like one
that discriminates. (D-194)
