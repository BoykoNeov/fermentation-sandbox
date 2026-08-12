---
name: feedback-a-ratio-guard-cannot-see-a-common-factor
description: A guard asserting a RATIO is blind to any term that scales both sides; and a red from a crashed fixture names nothing
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c5f901db-8b1f-4b2e-a63e-b0739efff44d
  modified: 2026-08-12T13:56:37.328Z
---

D-196 shipped `test_the_fenton_limb_draws_one_o2_per_acetaldehyde`, asserting `n_o2 ==
n_acetaldehyde` at machine tolerance — exact stoichiometry, the strongest kind of claim. It is
**blind to the term D-196's own §5 named**: a hydroperoxyl-recycling limb multiplies that limb's
whole flux by `1/(1−s)`, scaling numerator and denominator alike, so the ratio stays 1. Measured,
not assumed — a capped route-1 arm ran **GREEN** on it (D-198).

**Why:** a ratio is invariant under a common factor, which is exactly what a rate-law change
looks like. So the strongest-looking guard in the file — exact, machine-tolerance,
mechanistically stated — had a blind spot shaped precisely like the open item sitting next to it.
Choosing the assertion's *quantity* matters more than its tightness: a loose absolute check would
have caught what a 1e-12 ratio check could not.

**How to apply:** ask what transformation leaves your assertion's value unchanged, then check
that no candidate mechanism has that shape. If a term can multiply both sides, assert one side
against its **own inputs** instead — D-198's replacement pins the acetaldehyde flux to
`share × activation`, which only holds while the limb is a pure consumer of the node.

**And check what a RED names, not that it is red.** The uncapped route-1 arm was red on both
existing guards, which looked like coverage. It was the **fixture crashing** on the series' pole
— trajectory truncated, assertion never executed. That is why the capped arm exists: to make the
assertion actually run so the GREEN could be observed. A mutation arm that dies before the assert
tells you nothing, and reads exactly like a guard working.
Related: [[feedback-grep-finds-claims-not-guards]],
[[feedback-pair-the-red-with-an-ordering-preserving-baseline]],
[[feedback-prefer-the-variant-your-guards-can-see]].
