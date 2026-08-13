---
name: feedback-relocate-the-unsourced-factor
description: "When one coefficient in a rate law is unsourced, choose the algebraic form that puts it on the side that cannot reach your headline — don't assume it away"
metadata:
  node_type: memory
  type: feedback
---

**If a rate law has one coefficient nothing sources, write the law so that coefficient multiplies
the quantity you are NOT reporting.** An unsourced factor is not automatically a blocker and it is
not something to assert past — it is a constraint on the *algebraic form*, and the form is usually
free.

**The case (D-201).** The quinone–H₂S sink is 1:1 on the sulfide side (sourced: *Understanding Wine
Chemistry* puts thiols in the "condensation … to yield an adduct" category, not the "reduction by
bisulfite/ascorbate" one, so the two-electron redox route that would make it 2:1 is a different
reaction). But H₂S is **divalent**, so the adduct's residual –SH could take a *second* quinone, and
**nothing in the corpus addresses that.** Written the obvious way — quinone-first, `d(h2s)` derived
from `d(quinone)`, which is how the first probe wrote it and how the shipped sibling
`QuinoneSulfonation` reads — the 2× ambiguity lands directly on the reported sulfide removal and
doubles the headline. Written **sulfide-first**, the same factor multiplies only the quinone draw,
where it moves a 0.003 % share of a node nobody is reporting, and `d(h2s)/dt` is *exactly*
independent of it. Same chemistry, same numbers, one of them shippable.

**Why:** the alternative moves are both bad. Asserting the factor is 1 launders a guess into
`source:` ([[feedback-a-derived-yield-encodes-its-rate-law]]); declaring the beat blocked forfeits
a result the ambiguity does not actually touch. Relocation is the third option and it is usually
available, because a bilinear law can be written from either reactant's side.

**How to apply:** before writing the law, ask *which single number here is unsourced* and *which
number am I going to report*. If they are on the same side, flip the form. Then **guard the
property, not the value**: mutate the unsourced coefficient (1 → 2), verify the mutation LANDED on
the side it should move, and assert the reported quantity's rate law is bitwise unchanged. Beware
asserting the *trajectory* is bitwise unchanged — it usually is not, because the two sides couple
back through the shared pool; the honest bound is the Process's own share of that pool (D-201
measured 8.2e-7 against a 3e-5 share). And say in the docstring which side carries the ambiguity
and why, or the next person "simplifies" the law back to the obvious form and silently re-arms it
([[feedback-prefer-the-variant-your-guards-can-see]]).
