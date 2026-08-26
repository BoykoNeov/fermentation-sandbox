---
name: feedback-check-whose-fit-the-other-side-is
description: "Before calling a comparison 'two independent sources', check whether one side is your own repo's fit of the other — grep the source text for the number"
metadata:
  node_type: memory
  type: feedback
---

At D-232 I found that the engine's wine biomass yield and the beer trial it is scored against
disagree ~2x on cells per gram of assimilated nitrogen, and wrote it up as "two trials disagree".
A reviewer asked whose number the wine side actually was. It is **this repo's** regression, fitted
at D-13/D-14 from the paper's table and figure and carrying a **published-typo correction** to one
exponent. Checked against the paper's full text, already on disk: the value **"10.06" appears
nowhere in it**, and the author's own prose says only that the yield *"can be estimated from the
relationships shown in Fig. 4"*. D-219's table had described the same number as "his own Fig. 4
regression", which reads as the author's when it is ours.

**Why:** "two independent sources agree/disagree" is a much stronger claim than "our fit of one
source disagrees with another source", and the strength is what a later beat will act on. A
derived quantity re-derived through the same fit is one measurement counted twice — the same trap
as a cross-check that shares an input. Here the finding survived (the disagreement is real and its
sign holds across every band), but the *class* of the claim had to change, and the guard's wording
with it.

**How to apply:** whenever a comparison has a literature number on one side, **grep the source
text for that exact number before pairing it**. If it is not printed there, the side is a
construction: say so in the record, in the guard docstring, and in the memory row. Name what makes
it a construction (fitted, unit-converted, typo-corrected, read off a figure). Prefer the narrower
sentence — "our X, read as a count, disagrees with their Y" — over the flattering one. Related:
[[feedback-a-transcription-answers-more-than-its-purpose]], [[feedback-paywalled-is-one-host]].
