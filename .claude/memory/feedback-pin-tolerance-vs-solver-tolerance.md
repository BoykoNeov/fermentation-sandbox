---
name: feedback-pin-tolerance-vs-solver-tolerance
description: "A regression pin's tolerance is meaningless until measured against the integrator's own rtol; and pytest.approx passes on EITHER bound, so a small atol never binds"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6d4266a9-7185-4ae1-a691-3d3e542c4883
  modified: 2026-07-27T12:19:30.345Z
---

Before pinning numeric output of an ODE run, **read the integrator's tolerance and measure
the noise floor**. `simulate_scheduled` runs **BDF at rtol=1e-6 / atol=1e-9**, so the
"obviously safe, four orders of margin" pin of rtol=1e-6 I shipped in D-140 sat at *exactly*
the solver's error budget — zero margin. Measure by re-integrating four orders tighter and
differencing; that difference IS the shipped run's error.

Two traps that follow:

- **`pytest.approx(x, rel=R, abs=A)` passes on EITHER bound.** An `abs` smaller than
  `rel * x` can never bind. My `atol=1e-15`, written to protect a 2.3e-9 slot, was dead code
  under a 2.3e-15 relative band.
- **Near-exhausted pools break relative tolerances.** A slot that has faded to ~1e-6 g/L
  (or ~1e-9) carries relative solver error of 1e-3, not 1e-6, while its *absolute* error is
  tiny. Those slots need absolute floors (~10× measured absolute noise), not a looser global
  rtol that would blunt every well-behaved slot.

**Why:** an over-tight pin is not a strict guard, it is a flake — and a flake gets loosened
by whoever it blocks, which destroys the guard exactly when it was supposed to bite. Loosening
100× cost nothing here: rtol=1e-4 still fires 8 of 18 pins on a **0.1%** parameter change.

**How to apply:** state the measured noise floor in the file next to the constant, and say
what the tolerance must catch. Also expect state-vector *length* changes to perturb results
with no model change — BDF's error norm is RMS-weighted over the vector, so adding one slot
shifts step selection. Related: [[feedback-rejected-values-must-be-unreachable]],
[[feedback-conceded-caveats-are-not-coverage]].
