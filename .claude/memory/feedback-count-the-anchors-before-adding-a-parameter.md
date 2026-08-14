---
name: feedback-count-the-anchors-before-adding-a-parameter
description: "A second source into a node with ONE observable and one fitted constant is unidentifiable -- and 'it breaks the calibration' is a different, weaker objection that must be checked against the calibration frame's actual contents"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 76396680-20bf-4e86-af29-a50e3b00ae04
  modified: 2026-08-14T07:12:31.722Z
---

Before adding a source to an existing node, **count the node's observables against its free
parameters.** D-203 was going to refuse the ascorbate→sotolon route because
`k_sotolon_aldol` is *"CALIBRATED … against the ~2 mg/L alpha-ketobutyrate residual"* — i.e.
"a second source breaks the fit". **Measured, that was false**: the calibration wine carries
`ascorbate` = 0.0 at t=0 *and* at the end (the pool is dose-only, default-0), so the new route
contributes **exactly zero** in the frame the constant was fitted in. A rate constant fitted at
one concentration and evaluated at another is a rate constant working, not a violation.

**The objection that survives is IDENTIFIABILITY, and it is stronger:** the node has **one**
observable (a 5-20 µg/L anchor, model at 7.4461) and building the route creates **two** free
parameters (the fitted constant, itself only an author estimate, and the new yield). Every pair
reproducing the anchor is admissible ⇒ the yield is not identified by anything the model can see,
and the yield is exactly what sets the answer.

**Why:** the two objections point at different fixes. "Breaks the calibration" says *re-fit*;
"unidentifiable" says *get a second anchor*, and names which one. Shipping the wrong one invites a
correction and sends the next beat re-fitting a constant that was never disturbed.

**How to apply:** when a `notes:`/`uncertainty:` field says a constant was calibrated *against*
some pool, do **not** read that as a prohibition on touching the pool — it is usually about band
composition ([[feedback-a-notes-field-is-unchecked-storage]]). Run the calibration scenario and
print the pool the new source would feed; if it is zero there, the calibration objection is dead.
Then count anchors vs parameters. Pair the refusal with the **crossing value** — the yield/size at
which the unknown would start to matter — so a future source settles it in one step
([[feedback-relocate-the-unsourced-factor]], D-200's pricing idiom). Related:
[[feedback-validate-calibrations-in-the-frame-that-binds]],
[[feedback-a-derived-yield-encodes-its-rate-law]],
[[feedback-check-the-blocker-is-still-blocking]].
