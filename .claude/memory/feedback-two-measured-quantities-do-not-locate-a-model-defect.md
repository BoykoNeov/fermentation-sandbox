---
name: feedback-two-measured-quantities-do-not-locate-a-model-defect
description: "D-215 - comparing two MEASURED curves diagnosed a model defect that was not there; the model tracks its OWN version of the baseline, so compute that first"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c639a01c-94f5-44c3-a91b-491833d1c5c9
  modified: 2026-08-17T18:40:25.478Z
---

**A defect inferred from two source numbers is a claim about the model that never touched the
model.** At D-215 I read the measured acid courses against the measured extract curve, saw the
acids arriving later than the sugar was consumed, and concluded the model *front-loads* its
flux-linked acids. I told the owner so, with a plain-language table. Then I measured it: the
model books **20.5 %** of each acid's rise by day 2 where the measurement says 14.6 % (lactic),
45.9 % (succinic), −4.1 % (malic). Two of my three predictions missed at 90 % confidence, and the
one load-bearing claim built on them — *"the model is being flattered, so the real early deficit
is bigger"* — had the **opposite sign** and sat 3.8× inside the read-noise floor.

**Why:** the inference "measured-acid lags measured-flux ⇒ the model is early" silently assumes
**the model reproduces the measured flux**. This engine ferments that wort ~2.8× too slowly at
day 2, so its acids track a much slower baseline and the whole comparison was between the
source's clock and the source's clock. This is [[feedback-pair-the-arm-with-its-baseline]] in the
form where the wrong baseline is a *published curve* rather than another random draw — harder to
see, because both halves are impeccably sourced and the arithmetic is right.

**How to apply:** before saying "the model gets X's timing wrong", compute **the model's own X and
the model's own driver**, and put all four numbers in one table. If a diagnosis can be written
without running the model, it is a statement about the literature, not about the code. The good
side-effect: doing this is what surfaced the real, unscored defect (the extract-course gap) —
[[feedback-a-summary-statistic-is-not-the-curve]] one panel over. And when the surprise lands,
[[feedback-pre-register-the-cheap-prediction]] is what makes it a diagnosis instead of a
retrofit: the miss was legible only because the prediction was on disk first.
