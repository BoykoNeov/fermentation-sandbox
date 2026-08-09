---
name: feedback-compute-the-clean-fix-before-adopting-it
description: The structurally-cleanest candidate fix can be arithmetically wrong — compute its consequence on every affected case before adopting it
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 941cd08a-82a1-4c4f-a54f-48f6a4d187e1
  modified: 2026-08-09T11:30:33.612Z
---

When choosing between candidate fixes, **structural cleanliness is not evidence
of correctness**. A candidate that satisfies every stated constraint, removes the
defect by construction, introduces no new constants, and makes an existing
docstring literally true can still be **worse than what ships**. Compute its
actual consequence, on every affected case, before adopting it.

**Why:** D-164. The class-(d) defect was an epistemic band gating a scenario
override. The obvious fix — have the override mint a band recentred
multiplicatively on the stipulated value, preserving the base's *relative* width
— passed every check I could name in advance: the value stays inside its band
(so the sampler's `triangular(low, value, high)` precondition holds), the band
stays non-degenerate (so the parameter is still sampled), no number enters the
code (so prime directive 2 holds without a schema change). It was recommended,
and it was wrong. `k_autolysis`'s band is log-symmetric (×0.1 … ×10); recentring
that relative width on a high-edge override yields `[1e-3, 1e-1]`, whose
triangular mean is **3.71×** what the scenario asked for, against the shipped
**0.67×**. Only the arithmetic said so — every structural argument pointed the
other way.

**How to apply:** Implement the candidate as a scratch calculation (not in the
repo) and put its number next to the shipped number, per affected case, in one
table. Two cases minimum: a candidate that helps case A and hurts case B is the
normal outcome, and one case cannot show you that. Watch specifically for
**asymmetric or log-scaled bands** — a rule that reads as symmetric ("preserve
the relative width") behaves wildly differently on `[0.8×, 2×]` than on
`[0.1×, 10×]`.

Second half, equally load-bearing: **run the control before attributing a number
to your mechanism.** I nearly recorded "the ensemble mean is 3.07 against a
requested 2.2" as an override defect. Measured without any override, the same
statistic is *worse* (the cap is non-biting in 76.8 % of members vs 68.6 % with
the override) — it is a property of wide speculative bands, not of the override.
The consequential half of a finding needs its own baseline, or you ship a true
number attached to the wrong cause. See
[[feedback-measure-which-side-before-building]] (measure the sign before
building) and [[feedback-conceded-caveats-are-not-coverage]].
