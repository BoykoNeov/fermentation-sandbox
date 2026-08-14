---
name: feedback-a-pair-constrains-a-response
description: A published PAIR of values constrains a rate-free RESPONSE where a single value constrains only a level - but an accepted deviation on the level can eat the range the response needs
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 760f4220-cd88-4a64-91d5-a5002cca53b9
  modified: 2026-08-14T06:34:13.908Z
---

Two halves of one idea, both measured at D-202 on the quinone branching.

**A pair of published values constrains a response, and a response needs no absolute rate.** The
branching had sat "blocked on structure" for ~4 beats behind the sentence *"nothing in hand
adjudicates which"* — because every candidate was being scored against a **single** number, which
only ever constrains a level. A printed table cell giving the same statistic **with and without**
one added species turned out to constrain the *ratio between the two*, and that ratio is
independent of the absolute rate the model was blocked on. It ruled out one of the two candidate
closures outright, across that constant's whole declared band.

**But an accepted magnitude deviation can consume the dynamic range the response lives in.** The
published effect was a **halving** (2:1 → 1:1). The model's un-dosed baseline already sat at 1.108
— at the bottom of that gap before the new species existed — so the largest response it could ever
show was ~10 %, and pushing the new rate constant **100× above its printed value** still only
reached 0.898. The prediction that a ~30× constant would close it was wrong for a further reason:
the arithmetic behind it assumed a **fixed** concentration, and the real pool depletes when the
dose (0.34 mmol/L) is comparable to the total flux it must intercept (~0.42 mmol/L).

**Why:** "we accept this magnitude deviation" is normally read as a bounded, local concession. It
is not: it silently forecloses every *later* phenomenon whose signal lives inside the range that
deviation ate. And the reverse is the opportunity — when a level is unreachable, a ratio of two
published levels may still be reachable and may still discriminate.

**How to apply:** When a subject is "blocked on an absolute rate", go looking for a **pair** of
published values that differ by one intervention, and score the model on the *response* instead of
the level. Before predicting that response, check the baseline has room to show it — compute the
model's current position inside the published gap first, because a model already sitting at the
"after" value can never demonstrate the "before → after" move regardless of its rate constants.
And when a share argument is used to size a rate, check whether the pool **depletes** over the
window: a fixed-concentration share is an upper bound that stops transferring the moment the dose
is comparable to the flux. When you accept a magnitude deviation, record **what range it eats**,
not just its size. Related: [[feedback-a-hit-can-be-two-errors-cancelling]],
[[feedback-check-the-blocker-is-still-blocking]], [[feedback-a-named-pull-may-not-answer-the-question]],
[[ascorbate-quinone-route]].
