---
name: feedback-a-counterfactual-must-carry-its-anchors-modifiers
description: "A draw anchored to a rate is not anchored to the rate that RAN: rate modifiers are attached per Process name, so a counterfactual gets none of them and silently measures the pre-modifier frame"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: fed512ed-8463-4a86-8d85-68ff40d1b22f
  modified: 2026-09-01T08:54:44.872Z
---

**If you anchor a new draw to another Process's rate, you have anchored it to that rate's
*base* value, not to the rate the solver actually applied.** Multiplicative modifiers
(Arrhenius, carrying capacity, osmotic) are attached by Process **name**, so a Process written
this week carries none of them, while the quantity it claims to feed grows at the modified rate.
The result is a measurement in a frame nobody declared, and it is invisible: every balance still
closes, because a scalar on a conserving vector conserves.

D-260. D-259's growth-anchored counterfactual draws `w_i · base_dx · gate`. `for_growth(*also_scales)`
scales the growth Process plus the D-32 amino-acid swap; the carrying-capacity modifier names
the same pair; **neither names the precursor sink**. The frame mismatch was inside D-259's own
record — its protein-demand figure (580-1088 µM) is in the realised frame, its counterfactual
draw in the pre-modifier one. Attaching the modifiers moved the leucine bracket 13.1-22.0 % →
**21.3-33.7 %**, and D-104's 20.9 % went from sitting at the TOP edge to sitting below the
BOTTOM one. The refusal and the order correction survived; every number in them moved.

**Why:** a counterfactual exists to be compared against the shipped form, and a comparison across
two rate frames is not one. D-32 exists in this repo precisely to stop a growth-anchored refund
outrunning growth's realised draw — the same hazard, one Process over — so the discipline was
already written down and was simply not applied to a test-only Process.

**How to apply:** when a Process anchors to another's rate, enumerate the modifiers that name
that other Process and attach yours to them (`m.modifies` is editable per compile, and the
factories build fresh objects so the edit cannot leak). Then reconstruct the draw **with** the
factor product when you integrate it off the trajectory, and iterate `active_modifiers`, not the
whole dict — a disabled modifier is still in the set and the solver never applied it. A closure
control against the pool the solver actually emptied is what catches both mistakes; see
[[feedback-a-mutation-harness-must-snapshot-to-disk]] and
[[feedback-a-control-belongs-where-the-number-is-made]]. Related:
[[feedback-a-freedom-can-be-an-artefact-of-the-frame]], [[feedback-validate-calibrations-in-the-frame-that-binds]].
