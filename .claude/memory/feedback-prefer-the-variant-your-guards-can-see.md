---
name: feedback-prefer-the-variant-your-guards-can-see
description: "When two versions of a fix differ only in size, the one that leaves the existing assertion green is usually still carrying part of the bug"
metadata:
  node_type: memory
  type: feedback
---

When choosing between two versions of a fix, **prefer the one the existing guards can SEE**. A
variant that leaves the old assertion green is usually still carrying part of the defect — and
"the test will need strengthening anyway" is a verdict on the design, not a chore to work around.

**Why:** at D-193 the pre-registered design scaled routed sulfur by an existing 0.95 retention
fraction ("one particle, one fate" — principled-sounding). Its own prediction table recorded that
under it, `assert c_of(flow.delta) < 0.0` **stays green**, and filed that as a caveat. It was the
tell: the variant kept 5 % of the very mechanism the beat existed to discharge (an invented loss to
filtering/racking, operations the event does not model). The 1.0 variant made the same assert go
red and inverted it into a machine-precision pin at zero — **for 0.05 µg/L of output difference**,
so the choice was never about size, only about which claim the code makes.

**How to apply:** when a candidate is invisible to the tests that exist *for exactly this
behaviour*, ask what it still shares with the thing being fixed and name that share in one clause
before adopting it. If the clause is indefensible on its own terms, the variant is wrong. Not
"always take the bigger change" — guard-detectability is evidence about the *variant*, not merely
about the suite. Corollary, cheaper and sharper: **an assertion that passes both before and after
a fix is not a guard for either behaviour** — D-193's `< 0.0` check on a flow the fix turns into
float noise gave +6.7e-23 in one scenario and −1.2e-23 in another, so its verdict was rounding, not
physics. Pin at zero with a tolerance plus a positive control, never flip a sign check.

Related: [[feedback-compute-the-clean-fix-before-adopting-it]],
[[feedback-mutate-the-premise-before-building-the-guard]],
[[feedback-grep-finds-claims-not-guards]], [[feedback-pin-the-band-not-the-nominal]].
