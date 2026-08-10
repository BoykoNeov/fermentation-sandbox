---
name: feedback-pair-the-red-with-an-ordering-preserving-baseline
description: "A RED from a mutation arm is only evidence about the thing you flipped if the same-sized move WITHOUT the flip stays green — and attribute at ASSERT granularity, not node"
metadata:
  node_type: memory
  type: feedback
---

**When a mutation arm goes RED, it is not evidence about the property you broke until you re-run the
same-sized perturbation with that property INTACT.** Keep the perturbation, drop the flip.

**Why (D-171).** Round 2 of D-171 crossed five sourced orderings by moving parameters onto opposing
band edges, and three arms went red. Every one of those REDs was **magnitude, not the inversion**:

* `y_acetaldehyde_per_tannin` 0.06 → 0.12 (crossing 0.09) killed two trajectory pins. Halving it to
  **0.03**, which leaves the ordering intact, killed the **identical two**.
* Crossing `ethanol_tolerance_mlf` above `ethanol_tolerance` killed three MLF tests. Moving the MLF
  side *alone* — ordering preserved — killed the same three with **identical assertion values**
  (`2.3545520666619044`), and the other side alone was green.

Had I stopped at the red, the record would have claimed nominal-scoped coverage that does not exist,
which is [[feedback-name-guards-for-what-they-forbid]] in reverse: a test that never mentions the
sibling, credited with forbidding the inversion.

**How to apply.**

* For every RED arm, run one baseline of comparable size with the property preserved. **Baseline RED ⇒
  the arm is INCONCLUSIVE about the property — say so; do not upgrade it to coverage.** Baseline GREEN
  ⇒ the flip is what bit. Cost is zero when the arms are green, so make it conditional.
* When the arm had to move **two** parameters (because a single-parameter crossing is schema-illegal —
  see [[feedback-verify-the-restore-between-mutation-arms]]), the sharp decomposition is **each half
  alone**: each half preserves the property, only together do they cross. A failure is attributable to
  the property only if it fires under the pair and under **neither** half.
* **Attribute at ASSERT granularity, not node.** Set-differencing failing *node* names is only sound
  if one node has one cause. `test_integrated_wine_aroma_temperature_directions` carries four asserts
  with three causes; node-level subtraction would have deleted it as "magnitude" and hidden a real
  finding. Parse the `E  assert …` lines and the `file:line` of each failure, per arm.
* Run the **clean tree** over the same node set first. A node set that is not green unmutated cannot
  attribute anything ([[feedback-a-null-result-needs-a-positive-control]]).
* GREEN arms need no baseline — a null has no attribution problem.

**The general shape, which recurred three times in one beat at three different levels:** a RED was
credited to the mechanism under test when it belonged to something else — first a schema rejection,
then a magnitude, then a co-located assert. Each was caught by the same move: **check what the RED
names, not that it is red.** See [[feedback-mutate-the-premise-before-building-the-guard]] and
[[feedback-count-and-print-your-skips]].
