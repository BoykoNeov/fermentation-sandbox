---
name: feedback-a-control-needs-mechanical-reach
description: "A control with no mechanical reach on the quantity proves nothing even when it 'passes' — check it COULD have reproduced the effect"
metadata:
  node_type: memory
  type: feedback
---

Before reading a control as evidence, ask whether it was **capable** of reproducing the effect.
A control that cannot mechanically move the quantity is not a passed control — it is a tautology
wearing one.

**Why:** D-196 pre-registered D-195's control — scale the rate constant by the same factor, no
new term — and it moved the benchmark ratio **0.06 %** against the mechanism's 2.23 %. That looked
like "the shape does work". It was worthless: the beat had *already measured* that a hermetic O2
pool is **supply-limited**, so a *rate* change cannot move a ratio set by *stoichiometry*. The
control had zero reach by construction and would have "passed" against any mechanism whatsoever.

The real control had to move the same quantity by the same average amount through a different
structure: charge the identical extra O2 as a **flat fraction** instead of a state-dependent one.
That one agreed to **0.21 %** on the headline number — and **broke both asymptotes** the shipped
form preserves. Only then was "the shape is load-bearing" a measurement.

**How to apply:** name the channel the control acts through and the channel the mechanism acts
through. If they differ in *kind* (rate vs stoichiometry, magnitude vs timing, level vs slope),
the control is not ordering-preserving and a null from it means nothing. Build the control that
travels the **same channel** at the **same size**. Sibling to
[[feedback-pair-the-red-with-an-ordering-preserving-baseline]] and
[[feedback-a-tautology-can-smuggle-an-attribution]]: those catch a *stated* caveat being waved
through, this catches a control that was never in the running.

Corollary from the same beat: a mutation harness that **captures** a baseline (`_ORIG = Cls.method`
at import) is valid only until that baseline moves. Once the source shipped, the capture *became*
the new code, the "before" arm silently became "after", and nothing went red — it surfaced only
because an identity read 0.503 where it had to read 1.0. Define arms as **suppressions of the live
attribute**, not as captures. [[feedback-verify-the-restore-between-mutation-arms]]
