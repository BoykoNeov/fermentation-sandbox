---
name: feedback-a-guard-about-a-literal-must-read-the-literal
description: "A guard whose subject IS a hard-coded constant, but which keeps its own copy of that constant, compares two numbers neither of which the engine uses — it stays green through every mutation of the real one"
metadata:
  node_type: memory
  type: feedback
---

When the thing a test is *about* is a literal in the source, copying that literal into the test
makes the test structurally unable to notice it moving. Both numbers change together only if a
human remembers to change both, and nothing ever fails to remind them.

**Why:** D-243 guarded a claim about the bracket `[0.03, 0.15]` that `_apply_nitrogen_dependent_yield`
hard-codes as the band an ensemble draws `biomass_N_fraction` over. The guard opened with
`_OVERRIDE_BRACKET = (0.03, 0.15)` as a module constant. Two of its three mutation arms moved the
seam's own `low=`/`high=` — the exact edit the guard exists to catch — and it stayed **GREEN both
times**, because it was comparing a test-side copy against a range computed from parameters, with
the engine's value participating in neither side. Only the arm that widened a *parameter band* went
red, which made the guard look 1-for-3 rather than blind. Reading the override back off a compiled
scenario turned it 3-for-3 with no change to a single assertion.

This is [[feedback-a-half-pinned-read-is-green-until-the-quantity-moves]] one layer up: there the
guard read a shared baseline instead of a per-item one; here it reads a *copy* instead of the
original. Same signature — a green that means "did not look" rather than "looked and agreed".

**How to apply:** before writing a guard, name the quantity it forbids moving, then ask by what
call the test **obtains** that quantity. If the answer is "it is written in the test too", the
guard has no teeth on it — reach for it through the real path (compile the scenario and read the
attribute, load the parameter and read `.uncertainty`, import the constant) even when that costs a
compile. A literal in a test is fine as an *expected value* on the other side of an assertion; it
is never acceptable as the *subject*. And design at least one mutation arm that moves the literal
itself rather than something upstream of it — an arm that only perturbs inputs cannot distinguish a
guard that reads the constant from one that copies it.
