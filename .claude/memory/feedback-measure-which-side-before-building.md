---
name: feedback-measure-which-side-before-building
description: "Before building a corrective, measure which side of it the model is already on — a one-directional instrument only helps in one direction"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 08c245e5-68a2-45fd-8996-93fa173519d6
  modified: 2026-07-20T13:38:30.674Z
---

Before building any **corrective** — a cap, a floor, a clamp, a gate — measure
which side of it the model is **already on**. A one-directional instrument only
helps in one direction, so "the literature says X is mostly Y" never by itself
justifies building the thing that pushes toward Y.

**Why:** D-120 refused a scoped, requested, twice-carried-forward build
(`f_de_novo_isoamyl_alcohol`) in one measurement. The de-novo cap is
`gate *= (1 − f_de_novo)`, which can only *reduce* amino-acid sourcing — and the
model already attributed *less* to amino acids than the source measured, for all
three alcohols. The premise ("all three are de-novo dominated in-study, so it's a
class of error") was true and still did not imply the build. Two beats carried it
forward because nobody measured the sign.

**How to apply:** Before writing the parameter, run the model and put its current
value next to the target in one table with an explicit direction column. Then ask
two more questions the sign alone won't catch: (1) **reach** — can the mechanism
touch every branch the sourced number is defined over? (isoamyl's cap could not
reach its *larger* valine branch); (2) **instrument** — is it a rate knob on a
supply-limited quantity? Where a pool exhausts, total draw is fixed by supply and
no rate multiplier moves it (D-112's `(1−f)` ceiling; it generalises).

Corollary worth its own check: a parameter can be **load-bearing for one quantity
and inert for another**. The shipped 2-PE cap moves the realised share by ~0 yet
is essential to the instantaneous carbon-refund guard. When prose justifies a
parameter by effect A, verify it was measured on A — see
[[project-fermentation-sandbox]] and [[feedback-rejected-values-must-be-unreachable]].
