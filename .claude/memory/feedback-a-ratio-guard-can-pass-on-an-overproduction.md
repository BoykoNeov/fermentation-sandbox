---
name: feedback-a-ratio-guard-can-pass-on-an-overproduction
description: D-245 - a sourced share went red after a fidelity fix because the fix removed the 2x over-production that had been inflating its denominator; check which leg moved before calling it a regression
metadata:
  node_type: memory
  type: feedback
---

**When a fidelity fix turns a sourced ratio red, ask whether the fix removed the thing that was
passing it.** D-244 corrected where a nitrogen-dependent yield fit is evaluated, which roughly
halved biomass, and six fusel guards went red — read at the time as "the corrected yield costs the
sourced de-novo floor". D-245 measured it on one knob. The amino-acid draw (the numerator) did not
move **at all**: those pools exhaust whatever the biomass, so it is supply-limited and pinned to
five decimal places across the whole sweep. What moved was the **denominator** — and the same
correction took every higher alcohol in that must from ~2x its own meta-analytic anchor to ~1.1x.
The surplus was de-novo, sugar-sourced carbon sitting under the ratio. **The floor had been cleared
by an over-production, not by the supply structure the guard names.**

**Why:** a share is two quantities, and a guard on it is green whenever their ratio is right —
including when both are wrong in the same direction. That makes an over-production a *hiding place*
for a sourcing defect, and it makes the day you fix the over-production look like the day you broke
the model. The absolute levels were sitting in the suite the whole time, anchored per molecule; the
ratio guard simply never looked at them. Same family as
[[feedback-a-ratio-guard-cannot-see-a-common-factor]], one level up: there a common factor cancels,
here a common inflation is invisible because only one leg is asserted.

**How to apply:** when a ratio guard reddens after a change, **decompose it before diagnosing** —
print numerator and denominator separately across the change, and say which one moved. If the
denominator moved toward its own independent anchor, the RED is the guard finally becoming honest
and the correct action is to record the miss, not to restore the old behaviour. And where a ratio
is the asserted quantity, keep a **level** guard beside it on at least one leg; an undosed control
fixture that the change leaves bit-for-bit is the cleanest form of it, because the pair states what
the change actually did.
