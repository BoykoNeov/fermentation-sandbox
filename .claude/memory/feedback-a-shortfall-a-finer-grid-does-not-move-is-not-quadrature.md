---
name: feedback-a-shortfall-a-finer-grid-does-not-move-is-not-quadrature
description: "Refine the grid FIRST when a closure control fails: an error that survives a 5x finer mesh is a missing term in the reconstruction, not a trapezoid error, and the two have opposite fixes"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: fed512ed-8463-4a86-8d85-68ff40d1b22f
  modified: 2026-09-01T08:55:19.978Z
---

**A failing closure control has two very different causes and one cheap discriminator.** If the
integrated draws account for less than the depletion the solver realised, either the mesh is too
coarse for a sharp transient (D-103's trapezoid trap) or the reconstruction is missing part of a
draw. Refine the grid once: if the number barely moves, stop looking at quadrature — you are not
reconstructing what ran. Reaching for a looser tolerance at that point buries a real fault.

D-260. A λ=5 over-draw arm closed at **0.671**; the fix was a missing rate-modifier factor, which
took it to **0.9747**. A 5× finer grid then moved it to 0.9738 — i.e. nothing. That flat response
is what said the residual 1.8 % was a *second* missing term, and it was: a modifier that sat in
the set but was **disabled**, folded into the factor product though the solver never applied it.
Both faults were in the beat's own arithmetic, not in the model, and both would have shipped as
measurements. The pool was checked for a negative excursion first (0.000 µM, so not that) — one
cheap elimination before the expensive one.

**Why:** the two causes look identical in the number, and the wrong fix for each is available and
tempting — widen the tolerance, or refine forever. A closure control only earns its keep if a
failure is *diagnosed* rather than accommodated; the value of asserting it in the fixture is that
it fires before any consumer reads the number.

**How to apply:** on a closure failure, in order — (i) check the pool for a negative excursion
(that inflates the denominator and reads as a shortfall); (ii) re-run at 5× the mesh and compare;
(iii) if flat, diff the reconstruction against what the run applies, term by term — modifier
factors, per-Process enable flags, any parameter the Process overrides internally. Only tolerate
the residual *after* those, and say in the assert message what the tolerance is buying. Related:
[[feedback-a-control-belongs-where-the-number-is-made]],
[[feedback-a-counterfactual-must-carry-its-anchors-modifiers]],
[[feedback-pin-tolerance-vs-solver-tolerance]], [[feedback-a-non-vacuity-check-can-itself-be-vacuous]].
