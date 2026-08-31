---
name: feedback-a-bracket-on-one-axis-is-not-reachability
description: "Before slowing a rate to fix its timing, find out what ELSE that rate is load-bearing for — the repair can land its target number exactly and still be disqualified by what it strands"
metadata:
  node_type: memory
  type: feedback
---

**A counterfactual that brackets the measurement tells you the target is reachable *on that
axis*. It says nothing about whether the mechanism you are about to slow down is carrying
something else.** Find that out before building, by asking what consumes or depends on the rate
you are about to change — not after, from a red suite.

D-257. D-256 bracketed a timing defect (model 100 %, ungated 16-24 %, measured 42-54 %) and called
a mixed form reachable. It is: a blended nitrogen gate `f·N/(K_n+N) + (1−f)` at `f = 0.79` lands
**48.1 %** on the source's own middle must and holds the level anchor at 172.3 mg/L. It was still
reverted. Holding the level meant scaling the five rate constants ×0.4033, and that rate turned out
to drive the model's *only* consumer of the must's phenylalanine — so the amino acid stopped being
consumed (20.3 % → 65.8 % left) and a downstream aging aldehyde rose ~15×. **31 tests red.** The
escape route was measured too and is closed: holding the constants keeps the amino acid but caps
the repair at 90.7 % against the measured 44-52 %.

**Why:** a bracket is a statement about one output's range. Rates in a coupled model are almost
never load-bearing for one output only, and the couplings that matter are the *implausible* ones —
here, protein synthesis anchored to higher-alcohol production, which no one would have designed
and which only a sweep reveals. A repair validated on its target axis alone is validated against
the half of the model you were already thinking about.

**How to apply:** before building, sweep the constant you intend to move — over the factor the
repair will actually apply — and diff the WHOLE end state, not the target pool: precursor pools,
conservation margins, anything with a regression pin. Then check whether the escape route (hold
the constant, move only the new parameter) reaches the target inside the constant's own sourced
band; if it does not, the repair is blocked and that is the deliverable. Report the built-and-
reverted repair with its numbers — it converts "the owner should authorise this" into "this is
unreachable until X", which is a specific next beat. Keep the patch. Related:
[[feedback-a-refused-counterfactual-still-buys-a-bracket]] (which this corrects — the bracket was
necessary, not sufficient), [[feedback-closer-to-reality-decides]],
[[feedback-compute-the-clean-fix-before-adopting-it]].
