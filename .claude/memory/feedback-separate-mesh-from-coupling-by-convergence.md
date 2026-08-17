---
name: feedback-separate-mesh-from-coupling-by-convergence
description: A threshold cannot tell a small real coupling from an adaptive-solver mesh artifact - tighten the tolerance and check the difference CONVERGES to zero
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 319bf67b-4e83-4c11-9ffb-068bbcb68911
  modified: 2026-08-17T16:17:39.671Z
---

Adding an inert state slot (beer's wort oxygen) shifted **every other column** of the
trajectory. Worst was `h2s` at 1.1e-9 g/L — which sounds like nothing until you scale it: that
pool peaks at 3.4 µg/L, so the shift was **3.3e-4 relative, 300× the solver's `rtol`**. A
relative-tolerance assertion failed, and the obvious repair — loosen the bound until it passes —
would have been worthless, because **a genuine small coupling and a mesh artifact both fit under
any threshold chosen generously enough to pass.**

What separates them is not size but **behaviour under refinement**. Re-running at tighter solver
tolerance:

| `rtol` | worst difference |
|---|---|
| 1e-6 (shipped) | 1.105e-9 g/L |
| 1e-9 | 3.395e-12 |
| 1e-11 | 1.055e-13 |

It tracks `rtol` essentially linearly ⇒ converges to zero ⇒ **adaptive mesh, not a pathway.** A
real coupling would have flattened out at a non-zero floor.

**Why:** an isolability claim is the project's prime directive 3, and it is normally written as
"byte-for-byte". That phrase is true of the **derivative** (the Process writes one slot, checked
exactly) and **false of the integrated trajectory** — adding an equation changes error-controlled
step selection, so every column moves a little. Writing "byte-for-byte" of a trajectory is a claim
that will fail the moment someone measures it, and the failure looks like a physics bug.

**How to apply:** state the RHS claim and the trajectory claim **separately**. For the trajectory,
assert **convergence across two decades of `rtol`**, not a magnitude — that is the assertion that
actually discriminates, and it stays meaningful if the model later grows stiffer. Keep a loose
absolute ceiling alongside it, chosen where a *real* coupling would land (here 1e-4 g/L, two
orders below the 6-14 mg/L a real acid pathway would have moved), so the test still fails loudly
for the right reason. Watch for **relative** metrics on small pools: a nanogram on a microgram
pool reads as a catastrophe. Related: [[feedback-pin-tolerance-vs-solver-tolerance]],
[[feedback-a-gate-is-a-discontinuity-the-solver-probes]],
[[feedback-a-null-result-needs-a-positive-control]].
