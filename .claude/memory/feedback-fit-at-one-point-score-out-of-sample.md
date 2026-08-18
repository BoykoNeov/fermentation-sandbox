---
name: feedback-fit-at-one-point-score-out-of-sample
description: "Fitting one constant at ONE condition and scoring it at the conditions you withheld turns a fit into a test of the model's response; at D-223 two withheld temperatures landed within 6%"
metadata:
  type: feedback
---

When a source gives the same measurement at several conditions, fit the constant at **one** and
score it at the others. At D-223 beer's uptake rate was fitted only to Foster's 15 °C course,
and the model then landed **0.945 / 0.973×** the measured mean at 12 °C and 22 °C — conditions
the fit never saw. A joint fit across all of them would have produced better-looking agreement
and proved nothing.

**Why:** a constant fitted to every condition absorbs whatever the model's *response* gets
wrong, so the residual stops being diagnostic — the same failure mode as fitting a downstream
observable ([[feedback-fit-the-observable-not-the-consequence]]). Held-out conditions are the
only thing that separates "the level was wrong" from "the response is wrong", and that
distinction is usually what decides whether a *different* parameter has to be re-sourced. Here
it is what let D-217's refusal to re-source the activation energy stand.

**How to apply:** name in advance which condition is the fit and which are the score, and write
the prediction down first. Score the withheld ones with their own assert, separate from the
fitted one, so the test can fail for a reason other than arithmetic. State plainly in the record
that the fitted condition is not corroboration. And read the condition where agreement gets
WORSE as a finding, not a nuisance — D-223's 30 °C column inverted to 0.66× and that inversion
is the honest price of the fix, not a detail.
