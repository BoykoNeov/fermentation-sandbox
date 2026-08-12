---
name: feedback-a-returned-intermediate-is-a-rate-law-change
description: "In a quasi-steady-state pool, a limb that returns its own input feeds its production flux — a geometric series, not a stoichiometric tweak; \"upper bound\" inverted to lower"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c5f901db-8b1f-4b2e-a63e-b0739efff44d
  modified: 2026-08-12T13:56:58.670Z
---

D-196 wrote that `_O2_PER_ACETALDEHYDE = 1.0` was "an **UPPER BOUND** on the NET draw", because
the limb's leftover radical might hand an oxidising equivalent back. D-198 measured it: for the
fate three wine texts describe, the recycled H₂O₂ re-enters **the same quasi-steady-state node
the limb drew from**, so the balance is `F = A + s·F`, i.e. `F = A/(1−s)`. That draws
`s²/(1−s)` **more** O₂ than the shipped constant at every `s > 0`. The bound was a *lower* bound.

**Why:** the phrase "hands an equivalent back" invites arithmetic on one turn — subtract
something from a constant, get a smaller constant. But if the returned species is an
*intermediate whose pool is quasi-steady-state*, it does not reduce the per-turn cost; it
increases the **number of turns**. That is a change to which step is rate-determining, and every
calibration anchored on the old rate-determining step (here `k_o2_depletion_total`) goes with it.
The wrong framing also invites the wrong edit: dialling the constant down, which models nothing.

**How to apply:** when a candidate term returns any species that is *already produced upstream in
the same loop*, write the balance equation before estimating a magnitude. If the species has no
state slot — quasi-steady-state, partitioned by branch fractions — the return closes as a
geometric series, so look for `1/(1 − c·s)` and then **find its pole**: `s = 1` here is
unsulfited wine, a real operating point, where the arm stops being integrable at all. Check the
closed form against a brute-force turn-by-turn sum; they use different methods, so agreement is
a genuine control where "it looks right" is not.

Corollary on reporting: a number produced under a **clamp you put in yourself** is yours, not the
model's. The unsulfited draw printed 0.0 mg/L; varying that clamp over six orders of magnitude
moved it across 21.4 mg/L and the solver had failed outright. Report **undefined**.
Related: [[feedback-a-derived-yield-encodes-its-rate-law]],
[[feedback-a-shape-change-is-a-change-to-the-pair]], [[feedback-measure-which-side-before-building]].
