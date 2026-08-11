---
name: feedback-grep-finds-claims-not-guards
description: "Grepping for tests that pin the restriction you're removing returns prose and assertions indistinguishably; four 'expected red' hits were comments, and the beat's only real defect lived in that gap"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 99326074-907e-4a86-84c2-8cc80588bab0
  modified: 2026-08-11T12:46:28.356Z
---

When a change **removes a restriction**, grep for tests asserting that restriction before
running — that is [[feedback-mutate-the-premise-before-building-the-guard]]'s companion and it
is right. But the grep returns **claims**, not **guards**. Check each hit is an `assert` and not
a comment before counting it as coverage.

At D-179 I grepped for "beer has no pH system" and found four hits, listing all four as expected
reds. **None was an assertion**: three were comments, and the fourth asserted something else
that was still true. So the red list read as coverage where there was none.

**Why it matters:** the beat's only live defect lived exactly in that gap. Giving beer a
`cation_charge` slot opened a gate keyed on slot presence for *every* beer, including ones that
never asked for a pH — whose empty charge balance solves to pure water, pH 7.0, making a rate
factor `10^(3.3-7) ≈ 2e-4`. A **5000×** change from a scenario that supplied no pH. The suite
stayed green through all of it.

**How to apply:** after grepping, open the hits. A line that *describes* an invariant proves
nothing; only an executed assertion does. If a behaviour you are changing has no assertion,
that is the finding — write the test before the change, not after. I found this one by reading
the *passing* tests to understand **why** they passed, which is worth doing whenever a predicted
red comes back green: a green you predicted red is either a wrong prediction or a missing test,
and both are worth knowing which.

**Corollary, same beat: a guard is only as broad as the registry it names.**
`test_no_shipped_acid_reaches_the_general_branch` was written one beat earlier to fail the day a
3-pKa acid shipped. It shipped, and the test **passed** — because the test read `ACID_STATE`,
which the same commit had quietly re-scoped from "every acid" to "wine's acids". Re-scoping a
registry silently re-scopes every guard that names it. A false green is worse than a red: it
reports that the tripwire held.

Related: [[feedback-full-suite-before-green]], [[feedback-a-null-result-needs-a-positive-control]],
[[feedback-count-and-print-your-skips]], [[feedback-name-the-field-your-predicate-read]].
