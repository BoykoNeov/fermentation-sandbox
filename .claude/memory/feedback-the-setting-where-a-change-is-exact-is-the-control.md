---
name: feedback-the-setting-where-a-change-is-exact-is-the-control
description: "When a rate-law change is a pure reparameterisation in one setting and a real mechanism change in another, solve BOTH — the exact one proves the claim the fitted one can only assert"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ae98d8d8-bbb4-4e60-be15-48e20c8776e7
  modified: 2026-08-25T19:40:33.769Z
---

When a change to a rate law is an exact reparameterisation in one setting and a genuine mechanism
change in another, **solve both**. The exact setting is the control: it proves the claim that the
other setting can only assert.

**Why:** D-227 swapped the ester sink's driver from a Monod stand-in to the real evolved-CO₂ rate,
which forced its constant to be re-solved. Beer's factor came out 3.34686 by bisecting the engine —
a fitted number, defensible but unverifiable, and precisely the kind of thing that later reads as
tuning. Wine has one sugar slot and no catabolite repression, so there the two drivers differ by a
closed form: `1/(q_sugar_max · co2_yield · scale)` = 2.5252809. Wine's bisection returned 2.52530.
Agreeing to **six significant figures** with a value computed on paper is what turned "the constant
folded a speed knob into itself" from a claim in a comment into a measurement. And the same pair
sized the other half honestly: beer's solve sits 17.7 % off *its* analytic value, and that gap is
the mechanism — three sugars, repression, different CO₂ yields per gram. Without wine, the 17.7 %
would have been unattributable; without beer, the change would have looked like bookkeeping.

**How to apply:** before solving a constant numerically, ask whether some configuration of the same
model makes the change analytically exact — one species instead of many, a term switched off, a
degenerate limit. Solve that one too, compare against the closed form, and report the pair. The
agreement licenses the derivation; the *disagreement* in the other setting is the size of the
physics you added, and it is only quotable because the first half pinned the arithmetic. Related:
[[feedback-a-constant-lives-in-its-rate-laws-coordinates]] (a constant's meaning is its rate law's)
and [[feedback-compute-the-clean-fix-before-adopting-it]].
