---
name: feedback-assert-the-property-not-the-bit-pattern
description: "Exact float equality on a solver's output asserts your toolchain, not your model: state the structural guarantee separately and bound the numerical one in ULPs"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ea3f3770-ac9d-4299-8476-1248aca8b862
  modified: 2026-08-27T09:26:31.484Z
---

`assert solved == shipped` on the output of a root-find, an optimiser or any iterative solver is
a claim that **two toolchains take the same number of steps and round identically**. Nobody gets
that. A different libm, BLAS or scipy build enters the same bracket and stops a step away, so the
test pins the machine the literal was produced on and fails everywhere else.

**Why:** D-238 shipped a back-solved constant as a YAML literal and guarded it with exact equality
against the runtime re-derivation, with the message *"do not paper over this with a tolerance"*.
The instruction was protecting something real — a loose relative tolerance would pass a solver
that found a *different* root, or a literal that had drifted. But the assertion was not that
thing: it returned 0 ULP on Windows and 1–3 ULP on both CI Pythons, so it had **never passed on
CI**, from the record that introduced it.

**How to apply:** split the claim in two, because they have different guarantees.

* **The structural half** — what actually makes the result exact. Here it was an early `return`
  that skips the re-solve entirely at the nominal draw, so the value is kept by control flow and
  is bit-identical on every platform. That is the one to assert first, and D-238's own code
  comment had already said so; the *test* just had not caught up with the *rule*.
* **The numerical half** — bound it in **ULPs**, not in a relative tolerance. ULPs are the unit
  "same root, different toolchain" lives in, and the budget is defensible only if you *measure*
  the distance to the nearest real defect. Measured here: 4 ULP budget, and the smallest plausible
  defect (the target nudged one part in a billion) landed **4.4 million ULP** away. With a
  million-fold margin the choice is not delicate and the bound is a restatement, not a loosening.

Two traps that follow. **Assert that the structural guarantee is what delivered the result**, not
merely that the result is right — on the platform where the solver *does* reproduce the literal,
deleting the skip leaves the first assertion green. Making the solver raise and requiring the call
to succeed anyway settles it with no platform argument
[[feedback-a-control-needs-mechanical-reach]]. And when you restate a guard whose message forbids
loosening it, **say in the record which half of that instruction you kept** — the objection was
correct about what it protected, and a reader who sees only "the tolerance was widened" is right
to be suspicious.
