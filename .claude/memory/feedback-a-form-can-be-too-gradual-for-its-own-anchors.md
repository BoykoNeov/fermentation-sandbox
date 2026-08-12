---
name: feedback-a-form-can-be-too-gradual-for-its-own-anchors
description: "Before fitting a cited functional form, check its SHAPE can span the source's own anchors — arithmetic on the form, not a sweep"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2619e25d-f5ec-4f22-8c72-48d1de2dc5a3
  modified: 2026-08-12T08:23:18.984Z
---

A source can state two things that its own cited functional form cannot hold at once. At D-192
the Handbook said inhibition is negligible at 200 g/L *and* that a 300 g/L must can yield less
alcohol than a 200 g/L one. The literature's own substrate-inhibition form (Ghose & Tyagi's
Haldane) has inhibition group `S²/(K_S + S)`, which grows **1.51×** across that span. Holding
both statements needs **~19×** — a power-law steepness of about n = 7.3. No constant exists that
satisfies both, and a 2.8-order sweep of that constant could only ever show one end or the other.

**Why:** I spent a whole sweep looking for a value when the problem was the *shape*. The sweep
looked productive — it produced a clean monotone table — but every row was answering a question
that had no solution. One line of arithmetic on the form's own group, evaluated at the source's
two anchors, would have ruled the family out before any integration ran. The advisor had steered
me into that family and the arithmetic is what retracted it.

**How to apply:** When a source supplies both a form and two or more anchors, **before fitting
anything**, evaluate the form's characteristic group at each anchor and take the ratio. Compare
it to the ratio the anchors demand. If the form cannot span it, no parameter search will help —
say so, and go looking for a different family (usually: a *thresholded* one-sided form, which
buys arbitrary steepness by giving up smoothness at one point, and which also makes inertness
below the threshold structural rather than statistical).

The corollary that made D-192 shippable: a threshold placed at the edge of a keystone's
validated envelope turns "does this disturb the calibration?" from a tolerance argument into a
bitwise one — `max|dS| = 0.000e+00`, and the calibration RMSE unchanged to nine decimals. It
also kills a tuning hazard: a mild smooth brake *improved* agreement with the reference, and a
threshold means that improvement can never be quietly banked.

Related: [[feedback-a-band-is-per-parameter-a-claim-is-joint]],
[[feedback-compute-the-clean-fix-before-adopting-it]],
[[prohibitions-osmotic-high-sugar]].
