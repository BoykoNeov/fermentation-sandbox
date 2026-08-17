---
name: feedback-agreement-can-be-a-frame-difference
description: "A published number is measured in a frame; a term that improves agreement may only have moved the comparison, not the model"
metadata:
  node_type: memory
  type: feedback
---

**A measurement carries a FRAME, and agreement bought without moving the model is a frame
difference.** D-182 added dissolved CO₂ to beer's charge balance and the agreement with Tyrell's
pH rose from 42.7-62.2 % to 77.6-97.0 %, pre-registered and celebrated. D-208 read the method
Tyrell cite: MEBAK II 2.14 is *"pH (EBC)"* and Analytica-EBC 9.35's scope line is *"the
determination of pH at 20 °C of **decarbonated** beer"*. The published number has no dissolved
CO₂ in it. Re-scoring in that frame gives **43.2-62.9 %**, within ~0.7 pp of the pre-D-182 value,
and **0 of 59049 joint corners reach** where 7332 did. The shipped `> 0.70` floor needed the
sample to have retained **65 %** of saturation.

**Why:** the term was right and the comparison was wrong, which is the hardest version of this
error to see — nothing in the model is broken, the chemistry has a source, and the number moves the
way you predicted. Four beats quoted the improved headline. It also inverted a *diagnosis*: D-207
concluded the model "overshoots the fall on day 1"; in the settled frame day 1 is out by the same
amount in the **opposite** direction, and 85 % of a drop attributed to one term was attributed
off-frame.

**How to apply:**
- Before scoring a model output against a published one, **read the cited method's own scope line**,
  not the paper's prose. Standards publish their scope: MEBAK's site had the decarbonation sentence
  while the 24-text corpus had zero hits for "MEBAK" — the corpus was the wrong place, not a block
  ([[feedback-paywalled-is-one-host]], [[feedback-a-refusal-reads-like-a-gap]]).
- When a newly built term improves agreement, ask **which side moved**: the modelled quantity, or
  the frame it is being compared in. A term that shifts the model *and* the comparison in the same
  direction will look like a validation.
- Fix it by **naming the second frame**, not by deleting the term: rates must read the in-vessel
  quantity. Pin BOTH frames' numbers so they cannot be re-conflated ([[feedback-pin-the-band-not-the-nominal]]).
- A frame is usually a **BOUND, not a point** (a real decarbonation leaves a residue). Walk the
  one-parameter family; here days 4-7 were unreachable at every retained fraction, which is what
  made the shortfall frame-robust instead of arguable. Do **not** re-anchor a guard to the bound's
  edge — that pins your own choice ([[feedback-prefer-the-variant-your-guards-can-see]]); an
  `xfail(strict=True)` on the level states the gap without inventing a target.
