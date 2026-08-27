---
name: feedback-two-errors-in-the-same-frame-cancel
description: "A declaration defect and a release defect in the SAME frame cancel exactly; repairing one half alone is a regression"
metadata:
  node_type: memory
  type: feedback
---

`yan_mgl` was asked to be both total nitrogen and assimilable nitrogen (D-246 §2), and the
obvious repair was the compile seam: carve the amino-acid pools out of the declaration at their
*assimilable* nitrogen so the field means what a lab means. Self-contained, one function.

**It makes the outcome worse.** The seam carves out at TOTAL nitrogen and every in-run
deamination releases that same TOTAL — the same frame on both sides, so what the run makes
available equals what was declared, **exactly, for any dose**. Repairing only the declaration
leaves more ammonium behind while the untouched release frame still delivers all of it: measured
+15.28 mg N/L at a 0.5 g/L dose, 6.1 % of the declaration (D-248 §8).

**Why:** two errors of equal size and opposite sign are one error apart from being invisible.
Fixing the visible half converts a wrong-intermediate/right-answer into a right-intermediate/
wrong-answer, and the second is worse because nothing downstream is guarding it.

**How to apply:** before repairing a bookkeeping frame, find every place the same frame is
applied and check whether they cancel — and test the identity against the **user-facing
declaration**, not against a re-derivation of the internal channels, because that is the quantity
the cancellation is a promise about. If they cancel, the repair is all-or-nothing: price the
complete version and refuse the half. Related: [[feedback-price-the-refusal-not-just-the-error]],
[[feedback-closer-to-reality-decides]].
