"""Runtime: drives the deterministic core through time.

Between interventions the state evolves continuously (ODE integration); interventions
(add SO2, pitch MLF, step temperature, dose oxygen) are timed events that mutate
state or change the active Process set (handoff §1.4).

Three pieces layer on top of the pure core without changing it:

* :func:`simulate` — deterministic integration over a *fixed* Process set and parameter
  map (the thin ``solve_ivp`` wrapper everything else builds on);
* :func:`simulate_scheduled` — the event-driven loop (decision D-35): segments the run at
  timed :class:`ScheduledEvent` breakpoints, restarting :func:`simulate` after each state
  jump / Process-set reconfiguration / parameter change, and tracking the external-flow
  conservation ledger across the discontinuities. This is how a piecewise temperature
  schedule and discrete winemaking interventions (DAP, SO₂, racking, pitching) are driven;
* :func:`simulate_ensemble` — the stochastic wrapper (handoff §1.6, decision D-24) that
  samples parameters within their provenance uncertainty bands and reports median + spread.
"""

from fermentation.runtime.ensemble import Band, Ensemble, sample_parameters, simulate_ensemble
from fermentation.runtime.integrate import Trajectory, simulate
from fermentation.runtime.schedule import (
    ExternalFlow,
    ScheduledEvent,
    ScheduledTrajectory,
    simulate_scheduled,
)

__all__ = [
    "Band",
    "Ensemble",
    "ExternalFlow",
    "ScheduledEvent",
    "ScheduledTrajectory",
    "Trajectory",
    "sample_parameters",
    "simulate",
    "simulate_ensemble",
    "simulate_scheduled",
]
