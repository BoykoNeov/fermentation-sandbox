"""Kinetics-agnostic conservation and sanity checks over a trajectory.

The harness knows nothing about specific chemistry. A model supplies a
``quantity_fn(state) -> float`` that should be conserved (e.g. total carbon
across sugar, ethanol, CO2, and biomass), and these helpers assert it stays
constant along the trajectory. Encoding conservation as runtime/test assertions
is how we catch a model that quietly creates or destroys mass.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from fermentation.core.state import FloatArray
from fermentation.runtime.integrate import Trajectory

QuantityFn = Callable[[FloatArray], float]


def _evaluate(traj: Trajectory, quantity_fn: QuantityFn) -> FloatArray:
    return np.array([quantity_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])


def max_drift(traj: Trajectory, quantity_fn: QuantityFn) -> float:
    """Maximum absolute deviation of ``quantity_fn`` from its initial value."""
    q = _evaluate(traj, quantity_fn)
    return float(np.max(np.abs(q - q[0]))) if q.size else 0.0


def assert_conserved(
    traj: Trajectory,
    quantity_fn: QuantityFn,
    *,
    rtol: float = 1e-6,
    atol: float = 1e-9,
    label: str = "quantity",
) -> None:
    """Assert ``quantity_fn`` stays within tolerance of its initial value.

    Tolerance is relative to the initial magnitude plus an absolute floor, so it
    behaves sensibly when the conserved quantity is near zero. Raises
    ``AssertionError`` (so it reads naturally inside tests and runtime checks).
    """
    q = _evaluate(traj, quantity_fn)
    if q.size == 0:
        return
    tol = atol + rtol * abs(q[0])
    drift = np.max(np.abs(q - q[0]))
    if drift > tol:
        worst = int(np.argmax(np.abs(q - q[0])))
        raise AssertionError(
            f"{label} not conserved: drift {drift:.3e} > tol {tol:.3e} "
            f"(initial {q[0]:.6g}, at t={traj.t[worst]:.3g}h -> {q[worst]:.6g})"
        )


def assert_nonnegative(traj: Trajectory, variables: tuple[str, ...], *, atol: float = 1e-9) -> None:
    """Assert the named variables never go meaningfully negative.

    Concentrations and biomass are physical and must stay >= 0; a small negative
    excursion within ``atol`` is tolerated as solver noise.
    """
    for name in variables:
        series = traj.series(name)
        worst = float(np.min(series))
        if worst < -atol:
            idx = int(np.argmin(series))
            raise AssertionError(f"{name} went negative: {worst:.3e} at t={traj.t[idx]:.3g}h")
