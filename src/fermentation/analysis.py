"""Trajectory-level derived observables — the pH/TA/SO₂ readout layer (decisions D-18, D-22).

The scalar pH/TA/molecular-SO₂ functions are *pure* and live in the core
(:mod:`fermentation.core.acidbase`); these *series* helpers map them over a
:class:`~fermentation.runtime.integrate.Trajectory`, which is a runtime type, so they
sit one layer up — a thin observable module mirroring how :mod:`fermentation.units`
provides scalar conversions and how the benchmarks map ABV over ``traj.series("E")``.
The one-directional dependency rule holds: this top layer imports the core solver and
the runtime ``Trajectory``; nothing in ``core``/``runtime`` imports back.

pH carries no ``dpH/dt`` — it is reconstructed from each stored state column on demand,
so these helpers add no integration cost and stay honest about pH being a derived
algebraic function of the acid state, not an integrated variable. Report the tier with
:func:`fermentation.core.acidbase.ph_tier` (pass ``ParameterSet.tier_map()``).

In D-18 the acid slots are constant, but the pH series is **not flat**: ``Byp`` is core
state (the D-16 realised-yield diversion grows it 0 → ~3 g/L over a wine ferment) and is
read by the charge balance (include-by-reading), so with the cation frozen at its pitch
value pH **drifts mildly down** (~0.05–0.1 units) as fermentation proceeds — a realistic,
emergent demonstration that the solver responds to acid dynamics with no scripting.

CAVEAT on :func:`titratable_acidity_series` — the *initial* (must) TA is the fidelity-grade
value; the series **rises** ~3–4 g/L as ``Byp`` accumulates (the whole pool read as
titratable diprotic succinic), which runs *backwards* to real wine (TA flat-to-declining
during ferment). That is an upstream pool sizing/booking artifact (D-16/D-19), not the
solver; treat the end-of-ferment TA as an over-estimate / directional only. See
:func:`fermentation.core.acidbase.titratable_acidity`.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from fermentation.core import acidbase
from fermentation.core.state import FloatArray
from fermentation.runtime.integrate import Trajectory


def ph_series(traj: Trajectory, params: Mapping[str, float]) -> FloatArray:
    """pH at each stored time, by solving the charge balance per state column.

    ``params`` is the resolved ``{name: float}`` map (e.g. ``CompiledScenario.param_values``
    or ``ParameterSet.resolve()``); it must include the pKa parameters from
    ``acidbase.yaml``. Returns an array matching ``traj.t``.
    """
    return np.array(
        [acidbase.ph_of_state(traj.y[:, i], traj.schema, params) for i in range(traj.y.shape[1])]
    )


def titratable_acidity_series(traj: Trajectory, params: Mapping[str, float]) -> FloatArray:
    """Titratable acidity [g/L tartaric-equivalent] at each stored time.

    The TA counterpart to :func:`ph_series`, over the same acid state (no new state).
    """
    return np.array(
        [
            acidbase.titratable_acidity(traj.y[:, i], traj.schema, params)
            for i in range(traj.y.shape[1])
        ]
    )


def molecular_so2_series(traj: Trajectory, params: Mapping[str, float]) -> FloatArray:
    """Molecular (antimicrobial) SO₂ [g/L] at each stored time (decision D-22).

    Maps :func:`fermentation.core.acidbase.molecular_so2` over the trajectory: at every
    column it solves pH from the organic acids and partitions the dosed free SO₂, so the
    molecular fraction tracks the (mildly drifting) pH with no scripting — the SO₂
    counterpart to :func:`ph_series`. Returns g/L; convert to the conventional mg/L with
    :func:`fermentation.units.convert.gpl_to_mgl`. Report the tier with
    :func:`fermentation.core.acidbase.molecular_so2_tier`. With no SO₂ dosed (``so2_free``
    ≡ 0) this is identically zero — and, because SO₂ is carbon-free and outside the
    charge balance, dosing it leaves :func:`ph_series` and ``total_carbon`` unchanged.
    """
    return np.array(
        [acidbase.molecular_so2(traj.y[:, i], traj.schema, params) for i in range(traj.y.shape[1])]
    )
