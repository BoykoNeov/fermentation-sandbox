---
name: feedback-a-mole-table-is-two-mass-frames
description: "A mol% composition converted to mass has two readings; the consumer's own molar masses decide which, not the literature"
metadata:
  node_type: memory
  type: feedback
---

A published composition in **mol %** does not become a mass fraction until you
choose what each mole weighs, and for a polymer there are two defensible
answers: the **residue** it contributes to the chain, or the **free monomer**
that has to leave the pool to contribute it. They differ by the condensation
water — ~16 % for amino acids — and the choice is not a judgement about the
paper. **It is decided by whatever consumes the constant.**

D-267: Lange & Heijnen 2001's Table IV is mol %. The constant it answers,
`_D259_RESIDUE_SHARE_OF_PROTEIN`, is declared in residue mass, but the Process
that reads it subtracts its draw from the **free** amino-acid pool and the
engine's `MOLAR_MASS` is the free acid (leucine 131.175, not 113.160). In the
residue frame the source contradicts two of the five bracket edges; in the frame
the draw actually uses it contradicts **all five**, by 22–43 % on the mid
composition three records quote as their headline. Same table, opposite verdict
on three species.

**Why:** the frame question is invisible while the table and the constant are
both "composition", and it does not show up as a unit error — both readings are
grams per 100 g of protein and both are internally consistent. The only thing
that separates them is what the number is multiplied into downstream, which
lives in the code, not the source.

**How to apply:** when converting a mole composition to mass, find the consumer
before you convert, and derive the frame from **its** constants — then check
that derivation against the code in a test rather than asserting it in a
comment (D-267 checks `MOLAR_MASS` equals the free acid, and that free minus
residue is exactly one water). Ship one reading and pin the loser beside it:
this is [[feedback-a-units-fork-is-not-a-band]], and crossing the two into a
"band" would give a spread that is entirely the fork. Report the frame the old
constant was in as a finding — it says which earlier numbers were low and by
how much. Related: [[feedback-derive-the-papers-nitrogen-frame]],
[[feedback-agreement-can-be-a-frame-difference]].
