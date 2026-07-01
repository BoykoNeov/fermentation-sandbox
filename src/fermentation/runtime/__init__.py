"""Runtime: drives the deterministic core through time.

Between interventions the state evolves continuously (ODE integration); interventions
(add SO2, pitch MLF, step temperature, dose oxygen) are timed events that mutate
state or change the active Process set. This package will grow an event queue,
phase switching, and a stochastic ensemble wrapper (handoff sections 1.4, 1.6).

Today it provides two pieces the architecture needs to be runnable end-to-end: a
deterministic :func:`simulate` over a fixed Process set, and a :func:`simulate_ensemble`
stochastic wrapper (handoff §1.6, decision D-24) that samples parameters within their
provenance uncertainty bands and reports median + spread. Both layer on top of the pure
core without changing it; the event loop is still to come.
"""

from fermentation.runtime.ensemble import Band, Ensemble, sample_parameters, simulate_ensemble
from fermentation.runtime.integrate import Trajectory, simulate

__all__ = [
    "Band",
    "Ensemble",
    "Trajectory",
    "sample_parameters",
    "simulate",
    "simulate_ensemble",
]
