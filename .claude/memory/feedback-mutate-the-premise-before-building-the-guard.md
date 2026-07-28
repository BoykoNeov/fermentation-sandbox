---
name: feedback-mutate-the-premise-before-building-the-guard
description: "Before building a guard, mutate the thing it would protect and see if the suite already goes red — the matrix decides whether to build or refuse"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: b06cf888-62c7-48a2-a1df-778fa08411fe
  modified: 2026-07-28T14:56:01.582Z
---

When a record says "nothing guards X" and names building that guard as work, **do not build
it — measure it first.** Break X in the working tree (mutate the value, the scope, the
default), run the full suite, and let the result decide. This is the owner's own instruction:
*"given how this beat went, I'd run a mutation matrix on its premise before building it too."*

**Why:** "I read the tests and didn't see a guard" and "no guard fires" are different claims,
and the archive has now been wrong in **both** directions. D-153 asserted an unguarded margin
from a reading of the source; D-155's matrix caught every arm, so the proposed guard would have
been a decoration and the beat closed as a **REFUSAL** costing one line instead of a day. D-156
ran the same apparatus on the sampling-surface split and got four **GREENs**, so that guard was
genuinely owed and shipped. One method, opposite verdicts — which is exactly why it has to be
run rather than predicted.

**How to apply:**

- **Get the polarity right before reading the matrix.** For a *value* mutation, RED = the
  property is already held ⇒ refuse. For a *scoping* mutation (a file joins a merge list, a
  default flips), GREEN = nothing catches it ⇒ build. Write the polarity down before running,
  or the table reads backwards.
- **A RED only refuses the guard if it comes from a test that ASSERTS the property**, not one
  that merely consumes the content. A scenario that raises because a name went missing is
  content-consumption catching a scoping change by accident — one-directional, and it does not
  discharge the other direction.
- **One instance cannot license a claim about a family.** If a mutation matrix drops one of N
  files and goes red, sweep all N before writing "this is already caught."
- **Ship only what flips in its OWN named test.** Re-run each arm against just the new file and
  require the failure to land in the test written for it; record any arm caught only by a
  sibling rather than glossing it. This is what separates a guard from a decoration.
- **Mutate by editing the tree and reverting with `git checkout -- .`**, asserting the exact
  line content before replacing it (a drifted line number must abort, not silently mutate
  elsewhere), and verify the revert with `git diff --exit-code`.

Related: [[feedback-rejected-values-must-be-unreachable]],
[[feedback-validate-calibrations-in-the-frame-that-binds]],
[[feedback-conceded-caveats-are-not-coverage]], [[project-fermentation-sandbox]].
