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


def ibu_series(traj: Trajectory) -> FloatArray:
    """International Bitterness Units (IBU) at each stored time (decision D-64).

    IBU is, by definition, ~1 mg/L of iso-alpha-acids, so this is simply the ``iso_alpha``
    state (g/L) times 1000. The trajectory starts at the boil-derived iso-alpha wired in at the
    compile seam and **declines** as :class:`~fermentation.core.kinetics.hops.IsoAlphaAcidLoss`
    adsorbs it onto viable yeast during fermentation — so the finished-beer IBU is below the
    end-of-boil value, the ~5-20% wort-to-beer bitterness drop. Requires a beer trajectory (the
    ``iso_alpha`` slot is beer-only); an unhopped beer is identically zero.

    TIER (decision D-64, derived not asserted — pass ``ParameterSet.tier_map()`` to
    ``ProcessSet.tier_of('iso_alpha', ...)`` for the reported tier): the boil isomerization
    kinetics are sourced (Malowicki 2005, plausible), but the finished value also depends on the
    speculative ``hop_utilization_efficiency`` and ``IsoAlphaAcidLoss``, so parameter-tier
    propagation (D-1) caps the finished IBU at **speculative**. Unlike carbon, iso-alpha is off
    the conservation ledger, so this readout adds no invariant and dosing hops leaves
    ``total_carbon`` byte-for-byte unchanged.
    """
    return np.asarray(traj.series("iso_alpha"), dtype=np.float64) * 1000.0


def astringency_series(traj: Trajectory) -> FloatArray:
    """Free-tannin astringency at each stored time, as mg/L (decisions D-78/D-79).

    Astringency is a **taste** (a tactile/mouthfeel percept), *not* an aroma — so, exactly like
    ``iso_alpha``/IBU (:func:`ibu_series`, D-64), it is deliberately **excluded from the D-67 OAV
    odor lens** and read out here instead. And exactly like ``ibu_series``, this reads **no
    threshold**: it reports the **free (harsh) tannin** directly as mg/L, astringency being monotone
    in it. A calibrated perception-intensity mapping (an "astringency index" against a taste
    threshold) is the deferred refinement — reporting the concentration is the honest Tier-3 v1 (the
    aging axis is directional, not magnitudes).

    **Two harsh tannins, summed (D-79 extends the D-78 oak-only readout).** Free tannin is the sum
    of
    the grape **condensed** tannin (``tannin``, the dominant red-wine astringency) and the oak
    **hydrolysable** tannin (``ellagitannin``, D-78) — both taste astringent, so ``astringency =
    (tannin + ellagitannin) × 1000`` (g/L → mg/L, gallic/ellagic/catechin-equivalent). The
    **polymeric pigment** the D-79 condensation forms is *soft* and deliberately **excluded** — that
    exclusion is exactly what makes softening emerge (see below). A red wine with no oak reads grape
    ``tannin`` alone; an oaked white reads ``ellagitannin`` alone; both 0 ⇒ identically zero.

    The trajectory **softens** three ways. (1) The DOMINANT mechanism (D-79):
    :class:`~fermentation.core.kinetics.aging.TanninAnthocyaninCondensation` condenses grape
    ``tannin`` (with ``anthocyanin``) into soft polymeric pigment, drawing the free-tannin pool —
    and
    the astringency — down. (2) The acetaldehyde-bridged route (D-80):
    :class:`~fermentation.core.kinetics.aging.AcetaldehydeBridgedCondensation` *also* draws
    ``tannin``
    down, bridging it to ``anthocyanin`` via a dissolved-O₂ acetaldehyde linker — so
    **micro-oxygenation
    softens astringency** (the emergent O₂ → colour/mouthfeel link). (3) The oak contributor (D-78):
    :class:`~fermentation.core.kinetics.aging.EllagitanninOxidation` consumes ``ellagitannin`` to
    scavenge O₂ (oak protection). Oak ``ellagitannin`` also **rises** first as
    :class:`~fermentation.core.kinetics.aging.OakExtraction` diffuses it in toward its
    ``add_oak``-set
    ceiling. Even with all three mechanisms this remains **one-directional-per-pool** honest: grape
    tannin
    self-polymerization and tannin–ethyl–tannin bridging (the *other* softeners) are
    further-deferred
    beats, so anthocyanin is the limiting reagent and A–T condensation softens only modestly (the
    D-78/D-79/D-80 scope).

    TIER (derived not asserted — pass ``ParameterSet.tier_map()`` to ``ProcessSet.tier_of(...)`` for
    each pool's reported tier): the whole aging axis is speculative, so parameter-tier propagation
    (D-1) caps this readout at **speculative**. Like ``iso_alpha``, both ``tannin`` and
    ``ellagitannin`` are off every conservation ledger (grape-/wood-derived), so this readout adds
    no
    invariant and dosing them leaves ``total_carbon`` byte-for-byte unchanged.
    """
    tannin = np.asarray(traj.series("tannin"), dtype=np.float64)
    ellagitannin = np.asarray(traj.series("ellagitannin"), dtype=np.float64)
    return (tannin + ellagitannin) * 1000.0


def polymeric_pigment_series(traj: Trajectory) -> FloatArray:
    """Stable polymeric pigment formed by tannin–anthocyanin condensation, mg/L (decision D-79).

    As a red wine ages, :class:`~fermentation.core.kinetics.aging.TanninAnthocyaninCondensation`
    condenses free monomeric ``anthocyanin`` (with grape ``tannin``) into **polymeric pigment** — a
    softer-tasting, SO₂/pH-**stable** colour form. This readout reports that stable pigment as
    **mg/L
    of the anthocyanin condensed into it**: since in v1 condensation is the *sole* fate of
    anthocyanin, the pigment is exactly the anthocyanin the pool has lost, ``(anthocyanin₀ −
    anthocyanin(t)) × 1000``. It **rises** monotonically from 0 as the monomeric → polymeric shift
    proceeds — the young-purple → aged-brick-red evolution's stable fraction.

    **A post-hoc readout, not a state slot (the A420 discriminator, D-74).** Unlike ``A420`` — whose
    O₂ driver has *competing* sinks so its browning share cannot be reconstructed and had to be
    integrated as a slot — anthocyanin here has a **single** fate (→ pigment), so the pigment is
    reconstructible from the pool's own drawdown (the ``iso_alpha``/IBU readout pattern).
    ``anthocyanin₀`` is taken as the trajectory's first stored anthocyanin value (the compiled
    initial condition; the only processes that move ``anthocyanin`` are the two condensation routes,
    and it is a t0 grape must input). This stays valid through the D-80 acetaldehyde-bridge beat:
    that
    route adds a second *formation* pathway (:class:`~fermentation.core.kinetics.aging.\
    AcetaldehydeBridgedCondensation`, tannin–ethyl–anthocyanin) but anthocyanin's sole fate is still
    → pigment, so ``anthocyanin₀ − anthocyanin`` still equals the total pigment (direct + bridged).
    Only a future **bleaching** beat (a second anthocyanin fate to a *colourless* form) would break
    the identity and promote the pigment to a slot.

    Requires a wine trajectory carrying the ``anthocyanin`` slot (wine-only, D-79); a white / no-red
    wine (anthocyanin ≡ 0) reads identically zero. TIER: **speculative** (the condensation params
    are
    speculative; parameter-tier propagation, D-1). Off every ledger (grape-derived), so this readout
    adds no conservation invariant.
    """
    anthocyanin = np.asarray(traj.series("anthocyanin"), dtype=np.float64)
    anthocyanin_0 = anthocyanin[0] if anthocyanin.size else 0.0
    return (anthocyanin_0 - anthocyanin) * 1000.0


def color_series(traj: Trajectory) -> FloatArray:
    """Total anthocyanin-derived red colour at each stored time, mg/L (decision D-79).

    The point of the D-79 condensation beat is colour **stabilization**, so this index counts
    **both** colour-bearing forms: free monomeric ``anthocyanin`` (bright but bleachable) **and**
    the
    :func:`polymeric_pigment_series` it condenses into (SO₂/pH-**stable**). As the wine ages the
    free
    anthocyanin declines but the stable pigment rises, so **total red colour is retained as its form
    shifts labile → stable** — the correct young-purple → aged-brick-red picture (reporting only
    free
    anthocyanin would show colour wrongly *vanishing*). Reported as ``(anthocyanin +
    polymeric_pigment) × 1000`` mg/L, counting each condensed anthocyanin unit as contributing the
    same colour it did when free (just now stable — the equal-absorptivity v1 simplification).

    **CAVEAT — in v1 this is an algebraic identity, ``≡ anthocyanin₀ × 1000`` (a flat line).**
    Because
    the pigment is a *reconstructed* readout (``anthocyanin₀ − anthocyanin``, not an independently
    integrated slot) and condensation is anthocyanin's sole fate, the sum collapses to the constant
    initial anthocyanin — so this series does **not** independently verify the condensation Process
    (that is what :func:`~fermentation.core.kinetics.aging.TanninAnthocyaninCondensation`'s
    closed-form derivative tests do) and it plots as a flat line, not a curve. It is shipped for two
    honest reasons: it states the v1 *physics* claim (direct condensation loses no colour, only
    stabilizes it — the observable *dynamic* is the monomeric → polymeric **shift**,
    :func:`polymeric_pigment_series`), and it is **future-ready** — a later **SO₂/pH bleaching**
    beat
    (a second anthocyanin fate → colourless) makes the total genuinely *decline*, and promoting the
    pigment to an integrated state slot then makes this a non-tautological sum of two independent
    quantities (with a testable ``anthocyanin + polymeric ≡ anthocyanin₀`` conservation invariant).

    Requires a wine trajectory carrying ``anthocyanin`` (wine-only, D-79); a white wine reads
    identically zero. TIER **speculative**; off every ledger (grape-derived), adds no invariant.
    """
    anthocyanin = np.asarray(traj.series("anthocyanin"), dtype=np.float64)
    return anthocyanin * 1000.0 + polymeric_pigment_series(traj)


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
