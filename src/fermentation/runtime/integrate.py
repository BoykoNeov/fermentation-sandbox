"""Deterministic continuous integration of a Process set.

Stiffness is guaranteed in fermentation (fast CO2/acetaldehyde transients beside
months of aging), so the default solver is implicit and adaptive (BDF). This
wrapper is intentionally thin; the event-driven loop and stochastic ensembles
build on top of it.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


@dataclass(frozen=True)
class Trajectory:
    """Result of a :func:`simulate` run.

    ``t`` has shape ``(n_times,)`` in internal hours; ``y`` has shape
    ``(n_vars, n_times)`` matching ``schema``. ``tier_map`` records the derived
    confidence tier of each variable for honest downstream reporting — capped by
    the tiers of the parameters each Process read when ``simulate`` was given
    ``param_tiers`` (parameter-tier propagation, D-1), otherwise structural only.
    """

    schema: StateSchema
    t: FloatArray
    y: FloatArray
    success: bool
    message: str
    tier_map: Mapping[str, Tier]

    def series(self, name: str) -> FloatArray:
        """Time series of one (scalar) variable."""
        sl = self.schema.slice(name)
        block = self.y[sl, :]
        return block[0] if block.shape[0] == 1 else block

    def final(self) -> dict[str, float | FloatArray]:
        """The last state as a name -> value(s) mapping."""
        return self.schema.unpack(self.y[:, -1])

    def overall_tier(self) -> Tier:
        return min(self.tier_map.values(), default=Tier.VALIDATED)


def simulate(
    process_set: ProcessSet,
    params: Mapping[str, float],
    y0: FloatArray,
    t_span: tuple[float, float],
    *,
    param_tiers: Mapping[str, Tier] | None = None,
    t_eval: FloatArray | None = None,
    method: str = "BDF",
    rtol: float = 1e-6,
    atol: float = 1e-9,
    max_step: float = np.inf,
) -> Trajectory:
    """Integrate ``d(state)/dt = process_set.total_derivatives(t, y, params)``.

    Parameters mirror ``scipy.integrate.solve_ivp``; the implicit ``BDF`` default
    is deliberate (see module docstring). Returns a :class:`Trajectory` carrying
    the per-variable tier map derived from the active Processes.

    ``param_tiers`` (a ``{name: Tier}`` map, e.g. ``ParameterSet.tier_map()``) caps
    each variable's reported tier by the tiers of the parameters its Processes read
    (parameter-tier propagation, D-1). Omit it and the tier map is structural only
    — the Processes' own tiers, which over-reports confidence when they run on
    speculative parameters. Pass it on any path whose tier a user will trust.
    """
    schema = process_set.schema
    y0 = np.asarray(y0, dtype=np.float64)
    if y0.shape != (schema.size,):
        raise ValueError(f"y0 has shape {y0.shape}, expected ({schema.size},)")

    def rhs(t: float, y: FloatArray) -> FloatArray:
        return process_set.total_derivatives(t, y, params)

    sol = solve_ivp(
        rhs,
        t_span,
        y0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        dense_output=False,
    )
    return Trajectory(
        schema=schema,
        t=sol.t,
        y=sol.y,
        success=bool(sol.success),
        message=str(sol.message),
        tier_map=process_set.tier_map(param_tiers),
    )
