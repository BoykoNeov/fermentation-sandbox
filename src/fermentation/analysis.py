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
from dataclasses import dataclass

import numpy as np
from scipy.stats import rankdata

from fermentation.core import acidbase
from fermentation.core.state import FloatArray
from fermentation.core.tiers import Tier
from fermentation.runtime.ensemble import Ensemble
from fermentation.runtime.integrate import Trajectory

#: Sensitivity methods for :func:`attribute_spread`.
ATTRIBUTION_METHODS = ("src", "srrc")


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
    """Molecular (antimicrobial) SO₂ [g/L] at each stored time (decisions D-22, D-28).

    Maps :func:`fermentation.core.acidbase.molecular_so2` over the trajectory: at every
    column it solves pH from the organic acids, splits the dosed *total* SO₂ into
    acetaldehyde-bound vs free (D-28), and partitions the free share into its molecular
    fraction — so the readout tracks both the (mildly drifting) pH and the transient
    acetaldehyde-binding dip with no scripting. Returns g/L; convert to the conventional
    mg/L with :func:`fermentation.units.convert.gpl_to_mgl`. Report the tier with
    :func:`fermentation.core.acidbase.molecular_so2_tier`. With no SO₂ dosed (``so2_total``
    ≡ 0) this is identically zero — and, because SO₂ is carbon-free and outside the
    charge balance, dosing it leaves :func:`ph_series` and ``total_carbon`` unchanged.
    """
    return np.array(
        [acidbase.molecular_so2(traj.y[:, i], traj.schema, params) for i in range(traj.y.shape[1])]
    )


def free_so2_series(traj: Trajectory, params: Mapping[str, float]) -> FloatArray:
    """Free SO₂ [g/L] at each stored time — total minus acetaldehyde-bound (decision D-28).

    The analytically-measured free SO₂: maps :func:`fermentation.core.acidbase.free_so2`
    over the trajectory. During active fermentation the acetaldehyde peak sequesters SO₂ so
    free **dips** below the dosed total, then recovers as acetaldehyde is reduced — the
    emergent D-28 coupling. Equals the (constant) total whenever acetaldehyde is 0. g/L;
    convert with :func:`fermentation.units.convert.gpl_to_mgl`.
    """
    return np.array(
        [acidbase.free_so2(traj.y[:, i], traj.schema, params) for i in range(traj.y.shape[1])]
    )


def bound_so2_series(traj: Trajectory, params: Mapping[str, float]) -> FloatArray:
    """Acetaldehyde-bound SO₂ [g/L] at each stored time (decision D-28).

    The complement of :func:`free_so2_series`: maps
    :func:`fermentation.core.acidbase.bound_so2` over the trajectory. Rises with the early
    acetaldehyde peak and relaxes as it is reduced; ``free + bound`` equals the conserved
    dosed total at every column. g/L. CAVEAT (D-28): acetaldehyde-only binder, so this
    **under**-estimates real bound SO₂ (pyruvate / α-ketoglutarate / sugars unmodelled).
    """
    return np.array(
        [acidbase.bound_so2(traj.y[:, i], traj.schema, params) for i in range(traj.y.shape[1])]
    )


# -- ensemble spread attribution (sensitivity) --------------------------------
#
# Which sampled parameters *drive* an ensemble's output spread, and how does that
# spread partition across confidence tiers? This is a first-order variance
# decomposition computed post-hoc from a single :class:`Ensemble` — no extra
# integrations. Because the ensemble samples each parameter *independently* (the
# D-24 assumption, asserted there), the standardized-regression coefficients (SRC)
# are near-orthogonal, so their squares approximately sum to the regression R² and
# give a genuine variance split. What R² leaves on the table (``1 − R²``) is the
# nonlinear/interaction remainder — reported explicitly so the budget never reads
# as "everything explained" when the model (Monod, logistic gates, Arrhenius) is
# nonlinear. It belongs one layer up (top-level observable over a runtime
# ``Ensemble``), like the pH/SO₂ series above.


@dataclass(frozen=True)
class SpreadAttribution:
    """First-order attribution of an ensemble output's spread to its sampled inputs.

    ``variable``/``slot``/``time_index`` identify the scalar output analysed (a state
    variable, a chosen ``S`` slot, at one time column — final by default). ``method`` is
    ``"src"`` (standardized regression on the raw draws) or ``"srrc"`` (rank-transformed
    first — the robust fallback when the response is monotone-but-curved and plain SRC's
    R² is poor).

    The variance budget: ``per_param`` maps each *varying* parameter to its SRC² (a
    non-negative share of output variance); ``per_tier`` rolls those up by the parameter's
    confidence :class:`~fermentation.core.tiers.Tier`; ``unexplained`` is ``max(0, 1 −
    r_squared)`` — the interaction/nonlinearity remainder. ``per_param`` values plus
    ``unexplained`` sum to ≈ 1 (exactly under a linear, independent-input model), so the
    budget is readable as fractions of total spread. ``per_param_signed`` keeps the signed
    SRC so a caller can see *direction* (does raising this parameter raise or lower the
    output). Pinned parameters (zero-variance draws) are excluded — they explain nothing.

    **Degenerate (empty-budget) case.** When nothing varied — no parameter was sampled, or
    every sampled band was pinned — there is no variance to attribute: ``per_param`` and
    ``per_tier`` are empty and both ``r_squared`` and ``unexplained`` are ``0.0``. This is a
    deliberate exception to the sum-to-1 invariant above (an empty budget, not a budget that
    is entirely unexplained). ``unexplained`` is not keyed off exact-float output spread,
    because byte-identical member inputs need not produce byte-identical solver output across
    platforms (threaded-BLAS reduction order), and a hair of numerical noise is not real spread.
    """

    variable: str
    slot: int
    time_index: int
    method: str
    r_squared: float
    n_members: int
    per_param: Mapping[str, float]
    per_param_signed: Mapping[str, float]
    per_tier: Mapping[Tier, float]
    unexplained: float

    def ranked(self) -> list[tuple[str, float]]:
        """``(name, SRC²)`` pairs, largest contributor first."""
        return sorted(self.per_param.items(), key=lambda kv: kv[1], reverse=True)


def attribute_spread(
    ensemble: Ensemble,
    variable: str,
    param_tiers: Mapping[str, Tier],
    *,
    slot: int = 0,
    time_index: int = -1,
    method: str = "src",
) -> SpreadAttribution:
    """Attribute an ensemble output's spread to its sampled parameters (and their tiers).

    Reads the already-computed ``ensemble.members`` and ``ensemble.member_params`` — no
    new integrations. For each *varying* sampled parameter it standardizes the drawn
    values and the scalar output ``variable[slot]`` at ``time_index`` across members, fits
    an ordinary least-squares regression, and reports each coefficient's square (SRC²) as
    that parameter's share of output variance. ``param_tiers`` (e.g.
    ``ParameterSet.tier_map()``) maps each parameter to its confidence tier for the
    per-tier rollup — the Ensemble's own ``tier_map`` is per *state variable*, not per
    parameter, so the parameter tiers are passed in.

    This is a **first-order screen**: it captures the monotone, roughly-linear part of the
    response and lumps the rest into ``unexplained`` (``1 − R²``). It needs a decently
    sized ensemble to be stable — ``n_members`` ≳ 50–100; with only a handful of members
    the SRCs are noise, and with fewer members than varying parameters the fit is
    underdetermined (this raises). If ``R²`` is low on a strongly curved response, retry
    with ``method="srrc"`` (rank-transformed), which recovers monotone nonlinearity.
    """
    if method not in ATTRIBUTION_METHODS:
        raise ValueError(f"unknown method {method!r}; expected one of {ATTRIBUTION_METHODS}")
    if variable not in ensemble.schema:
        raise KeyError(f"no variable {variable!r} in schema")

    sl = ensemble.schema.slice(variable)
    width = sl.stop - sl.start
    if not 0 <= slot < width:
        raise ValueError(f"slot {slot} out of range for {variable!r} (has {width} slot(s))")
    row = sl.start + slot

    n = ensemble.n_succeeded
    # Output vector: the chosen scalar for every surviving member.
    y = ensemble.members[:, row, time_index].astype(float)

    # Design matrix: one column per parameter that actually varied across members. A pinned
    # parameter (in ``sampled_names`` but drawn constant, e.g. a zero-width band) has zero
    # variance and explains nothing, so it is dropped rather than dividing by zero.
    names: list[str] = []
    cols: list[FloatArray] = []
    for name in ensemble.sampled_names:
        x = np.array([mp[name] for mp in ensemble.member_params], dtype=float)
        if np.std(x) > 0.0:
            names.append(name)
            cols.append(x)

    empty_tiers: dict[Tier, float] = {}
    if not names or np.std(y) == 0.0:
        # Nothing varied, or the output has no spread — no variance to attribute, so the
        # budget is empty and unexplained is 0.0 (an empty budget, not an all-unexplained
        # one). Do NOT key unexplained off np.std(y) > 0: with no sampled parameter varying,
        # any residual y-spread is solver noise across byte-identical inputs (threaded-BLAS
        # reduction order differs by platform), not attributable variance.
        return SpreadAttribution(
            variable=variable,
            slot=slot,
            time_index=time_index,
            method=method,
            r_squared=0.0,
            n_members=n,
            per_param={},
            per_param_signed={},
            per_tier=empty_tiers,
            unexplained=0.0,
        )
    if n <= len(names):
        raise ValueError(
            f"underdetermined attribution: {n} members <= {len(names)} varying parameters; "
            "use more members (n_members >= 50-100 is recommended for a stable fit)"
        )

    x_mat = np.column_stack(cols)
    if method == "srrc":  # rank-transform inputs and output first (robust to curvature)
        x_mat = np.column_stack([rankdata(c) for c in cols])
        y = rankdata(y)

    # Standardize both sides (mean 0, unit sample std) so the OLS coefficients ARE the
    # SRCs; centring removes the need for an intercept term.
    x_std = (x_mat - x_mat.mean(axis=0)) / x_mat.std(axis=0, ddof=1)
    y_std = (y - y.mean()) / y.std(ddof=1)
    beta, *_ = np.linalg.lstsq(x_std, y_std, rcond=None)

    resid = y_std - x_std @ beta
    ss_tot = float(y_std @ y_std)  # == n - 1 after standardization
    r_squared = 1.0 - float(resid @ resid) / ss_tot

    per_param_signed = {name: float(b) for name, b in zip(names, beta, strict=True)}
    per_param = {name: float(b * b) for name, b in zip(names, beta, strict=True)}
    per_tier: dict[Tier, float] = {}
    for name, share in per_param.items():
        tier = param_tiers.get(name, Tier.SPECULATIVE)
        per_tier[tier] = per_tier.get(tier, 0.0) + share

    return SpreadAttribution(
        variable=variable,
        slot=slot,
        time_index=time_index,
        method=method,
        r_squared=r_squared,
        n_members=n,
        per_param=per_param,
        per_param_signed=per_param_signed,
        per_tier=per_tier,
        unexplained=max(0.0, 1.0 - r_squared),
    )
