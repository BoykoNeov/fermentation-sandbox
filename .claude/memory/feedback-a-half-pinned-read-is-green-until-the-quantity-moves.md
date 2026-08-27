---
name: feedback-a-half-pinned-read-is-green-until-the-quantity-moves
description: A guard that reads a shared baseline instead of the per-item one stays green for as long as the two agree; adding the first thing that moves them apart surfaces the old bug as if it were your new breakage
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ea3f3770-ac9d-4299-8476-1248aca8b862
  modified: 2026-08-27T08:45:59.973Z
---

A test that compares a per-item result against a **shared** baseline is only correct while the two
are equal. It can be wrong for records at a time and green the whole way, and it goes red on the
beat that first moves the quantity it reads — which makes an old defect arrive looking like a
regression in the change you just made.

**Why:** `test_scheduled_ensemble_conserves_across_jumps_per_member` checks
`final == initial + Σ flows` for every ensemble member, carefully using that member's own sampled
accounting fractions — and read `initial` off the **compiled** `y0` rather than the member's. That
was already wrong when D-233 began re-deriving part of `y0` per member, and wrong again after
D-236, and green both times, because the slots those records re-seed (`cation_charge`, `copper`)
carry neither carbon nor nitrogen. D-241 drew `must_fermentable_fraction`, which moves `S[0]` —
the carbon denominator itself — and the guard reported 99.043 against 97.834 g C/L: the whole
sugar band reading as a conservation breach. The model was never wrong. The guard was comparing
two different members.

**How to apply:** when you extend a per-item mechanism to cover a new quantity, grep every
consumer of that mechanism for reads of the *shared* baseline and ask whether each one should be
per-item — do this before running the suite, so a red is a confirmation rather than a diagnosis.
Then close the hole for good: after fixing such a read, add an assertion that the per-item values
actually **differ** from the shared one, because otherwise the fixed line is textually
indistinguishable from the line it replaced and will pass again the next time the two happen to
agree [[feedback-a-control-needs-mechanical-reach]]. Also treat a red of this shape as a
*question*, not a verdict — "did I break conservation, or was this guard measuring the wrong
thing" — since the reflex fix is to relax the tolerance, which would have buried a correct model
under a broken test.
