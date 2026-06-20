"""Runtime: drives the deterministic core through time.

Between interventions the state evolves continuously (ODE integration); interventions
(add SO2, pitch MLF, step temperature, dose oxygen) are timed events that mutate
state or change the active Process set. This package will grow an event queue,
phase switching, and a stochastic ensemble wrapper (handoff sections 1.4, 1.6).

Today it provides the minimal piece the architecture needs to be runnable
end-to-end: a deterministic :func:`simulate` over a fixed Process set. The event
loop and stochastic wrapper layer on top of this without changing the core.
"""

from fermentation.runtime.integrate import Trajectory, simulate

__all__ = ["Trajectory", "simulate"]
