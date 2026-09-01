---
name: feedback-two-sourced-targets-can-be-one-knob
description: "Before calling one side of a ratio 'untested', check whether the other sourced observable is the SAME number read from the other end — a pinned total plus a pinned denominator makes two targets one knob"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: fed512ed-8463-4a86-8d85-68ff40d1b22f
  modified: 2026-09-01T08:54:28.315Z
---

**A repair that moves a split may be moving a second sourced observable in lockstep, in the
opposite direction.** When the pool a split divides is fully consumed whatever draws it, and the
denominator the other observable is measured against does not respond to the repair, the two
observables are one degree of freedom: `B = (1 − A) · consumed / denominator`. Then "fix A
without touching B" is not a hard repair — it is arithmetically impossible, and no mechanism
that only re-anchors the same draw can escape it. Check the two invariances **before** designing
the repair or declaring the other side unexplored.

D-260. D-259 §5 refused un-inverting the protein split *by cutting the Ehrlich draw* (it breaks
Rollero's ¹³C-leucine tracer) and recorded the numerator side — raising growth's own draw — as
"untested", expecting it to move the split "while leaving the Ehrlich draw, and therefore the
tracer, untouched". Measured: the leucine pool ends at 0.00 % in **every** arm and the isoamyl
total moves < 0.01 % under every numerator arm, so the split and the tracer trade one-for-one.
The shipped sink sits at the Crépin end (split 81.50 %, tracer 1.507 % against a measured
3.4-8.2 %); the growth-anchored form sits at the Rollero end (27.58 %, 5.900 %). An over-draw
of λ=5 reaches Crépin's band at 73.29 % and pays 2.176 % on the tracer, exactly on the line.
**Neither form is the other's repair.** What would relieve it is the DENOMINATOR — joint
satisfaction needs ≤ 1170 µM of isoamyl where the fixture makes 2123 and the source's own
ferments print 793-1365 — i.e. a commensurability problem, not a split problem.

**Why:** a fence and a degeneracy read the same in a summary ("you can't have both") and imply
opposite next beats. A fence says *find another route to the same target*; a degeneracy says
*the target pair is unreachable until something else moves*, and names what. Recording the first
where the second is true sends the next beat looking for a mechanism that cannot exist.

**How to apply:** write the identity relating the two observables, then measure its two inputs
across every arm — is the shared pool fully consumed, and is the denominator inert? Pin both in
the guard, not the identity itself (the identity alone is a tautology if you read the split off
the same pool you read the total from — build the split from the two draws' own integrals and
let a closure control tie it back). Report the exchange rate and the denominator value that
would dissolve the collision. Related: [[feedback-a-tautology-can-smuggle-an-attribution]],
[[feedback-a-summary-statistic-is-not-the-curve]], [[feedback-two-measured-quantities-do-not-locate-a-model-defect]],
[[feedback-a-bracket-on-one-axis-is-not-reachability]].
