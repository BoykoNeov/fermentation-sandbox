---
name: feedback-negligible-sensitivity-is-not-no-freedom
description: "A term whose sensitivity is 2.1% still had a 0.022 pH free normalisation, because the constant's admissible RANGE multiplies a large concentration"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 41ccdfd6-1486-4600-8cb4-43288814b4cf
  modified: 2026-08-17T12:43:50.206Z
---

"This term barely responds, so a constant will do" is **two different claims**, and measuring the
first does not settle the second. Sensitivity is *how much the term moves when conditions move*;
freedom is *how far the answer moves across the admissible values of the constant you would have to
pick*. A flat term with a large magnitude has tiny sensitivity and a large freedom.

**Why:** at D-210 the dosed phosphate's buffering was measured at **0.97 mEq/L/pH against wine's own
~47 — 2.1 %**, which is a real argument for replacing a state slot plus two pKa parameters with one
negative increment to an existing slot: one line, no ledger churn, and an existing structural
guarantee left intact. But the cheap form needs the anion charge per phosphate **chosen**, and its
admissible range is 0.876 (pH 3.0) to 1.000 (fully dissociated). That 12 % of a *large* charge
(8.3 mmol/L against a wine holding ~42) is **0.022 pH** of end-state spread — resolvable, and a free
normalisation of the kind D-205 refused a whole term for. The sourced pKas have no such freedom, so
the measurement licensed the more expensive form rather than my preference doing it.

**How to apply:** when weighing a cheap constant against a sourced functional form, run **two**
numbers: the sensitivity (which tells you whether the physics matters) and the **span of the answer
across the constant's full admissible range** (which tells you whether the *choice* matters). Report
both. If the span is above the noise floor, the form is what removes the freedom and the extra slot
is paid for; if it is below, the constant is not free in any consequential sense and the cheap form
is the honest one. Never let "it barely buffers" stand in for "any constant will do".
[[feedback-a-normalisation-is-a-free-parameter]] is the same failure one level up — a ratio hiding
its reference. [[feedback-compute-the-clean-fix-before-adopting-it]] is why you price both.
