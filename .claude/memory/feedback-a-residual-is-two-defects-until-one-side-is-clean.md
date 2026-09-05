---
name: feedback-a-residual-is-two-defects-until-one-side-is-clean
description: "Model-minus-source is the SUM of both sides' defects, so its per-item sign can be a property of neither — score a rate law in the frame of the quantity it is paid on, and treat one identical modelled number across several items as the tell"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 18f4cecb-8304-41a0-86db-da312a32a929
  modified: 2026-09-05T22:20:59.873Z
---

**A residual between model and source is not a finding about the thing you are studying until
the other side of it is clean.** `model − source` is the sum of every defect on both sides, and
when the second defect is shared while the first varies per item, the residual's *sign* can be a
property of neither.

D-215 §3 scored beer's three produced acids — lactic, malic, succinic — by comparing each
measured course against the modelled one, day by day, and read the sign off the difference:
succinic late by 25 points, malic early by 25, lactic nearly right. It concluded the timing
errors **oppose** and therefore that no single correction could help all three, and that
conclusion sat in a record, two test docstrings and a data comment for 58 beats.

It was the frame. The modelled column carried the engine's own ~3× day-2 speed deficit — a
separate, known, deliberately parked defect. Scored against the *source's own* measured sugar
curve instead (both halves transcribed from the same four ferments of the same wort, and both
already in the repo), **all three acids trail the sugar, in one direction**, ordered succinic
< lactic < malic: +13.5 / +44.8 / +63.5 points at day 2. The shared speed deficit was +40.7,
which lands *inside* that span — so subtracting it turned the smallest lag negative and left the
largest positive, and manufactured an opposition that exists nowhere in the data (D-273).

**Why:** the arithmetic is a tautology — `a − c = (a − b) + (b − c)` for any three numbers — so
nothing about the residual warns you it is composite. It runs, it has a plausible size, and it
has a *sign*, which is the most convincing-looking thing a measurement can hand you. Worse, the
composite residual is stable and reproducible, so re-running it only confirms it. D-215 caught
this exact error one section later (§7, where a measured acid course was paired with a measured
flux when the argument needed the model's) and still shipped it in §3 pointing the other way.

**How to apply:**

* **Score a rate law in the frame of the quantity it is paid on, not on the calendar.** A
  producer written as `Y · flux` makes a claim about *sugar*, not about *time*. Plot it against
  the fermented fraction and the model's speed drops out of the comparison entirely — which
  matters most when that speed is a defect you are not allowed to fix.
* **Ask what the source can be scored against itself.** If a paper published two curves for the
  same ferments, comparing them contains no model quantity at all, and no defect of yours can
  reach it. That is the strongest evidence available and it costs one arithmetic pass.
* **The tell: one identical modelled number across several items.** D-215's modelled column read
  20.5 % for *all three* acids, because they share one rate law and are literally one curve
  (1.1e−13 points apart). A column that carries no per-item information cannot support a per-item
  conclusion — every bit of the variation in the residual came from the source, minus a constant.
  If your model column is one repeated number, stop before reading signs off the difference.
* **When you do explain a residual by decomposing it, guard the part that could have gone the
  other way.** "The deficit lands inside the span of the lags" is falsifiable and was asserted;
  the decomposition itself is not, and must not be offered as evidence.

Related: [[feedback-pair-the-arm-with-its-baseline]],
[[feedback-a-gap-needs-both-sides-from-the-same-run]],
[[feedback-a-limitation-can-belong-to-its-frame]],
[[feedback-a-summary-statistic-is-not-the-curve]].
