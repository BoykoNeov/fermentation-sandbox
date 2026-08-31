---
name: feedback-a-refused-counterfactual-still-buys-a-bracket
description: "Run the naive repair even when you will not ship it: the shipped form and the crude one BRACKET the measurement, which decides reachability before anyone designs a replacement"
metadata:
  node_type: memory
  type: feedback
---

**When a defect's repair is out of scope, still run the crudest version of it — not to propose
it, but to find out whether the measurement is reachable at all.** The shipped form and the crude
form are two points; if the source's number lies *between* them, some intermediate form can hit
it, and that is a result the owner can act on without a design.

D-256. The model makes 100 % of its isoamyl alcohol before the nitrogen gate shuts; Rollero
measures 42-54 % made by then. Deleting the gate (`nitrogen_gate = 1.0`) is not a candidate — it
overshoots the total 11.5× because `k_isoamyl_alcohol` is anchored against the gated form — but
run anyway it puts the fraction at **16-24 %**. The measured 42-54 % sits between the two, so a
mixed gated/de-novo form can reach it. The same run also showed the nitrogen *response* going
2.27× → 0.99× against a measured 0.76×, which attributes the level miss to the same gate as the
timing miss and cost one command to learn.

**Why:** "this needs a rate-law change, which is the owner's call" is a correct place to stop and
a bad place to stop *empty-handed*. Without the bracket the owner is being asked to authorise a
build whose target may not be reachable by any form of it — and a beat that refuses a repair
without pricing it cannot tell the difference between "expensive" and "impossible".

**How to apply:** define the crude counterfactual as the one-line removal of the suspected cause,
run it, and report three columns — shipped, counterfactual, measured — for every quantity the
finding names. State the overshoot that makes it unshippable in the same table, so nobody reads
the bracket as a recommendation. Keep the counterfactual out of the tree: apply it, measure it,
restore, and confirm `git diff` on `src/` is empty before committing. Related:
[[feedback-measure-which-side-before-building]], [[feedback-compute-the-clean-fix-before-adopting-it]],
[[feedback-closer-to-reality-decides]].
