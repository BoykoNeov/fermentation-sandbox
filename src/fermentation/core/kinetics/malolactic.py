"""Malolactic fermentation (MLF) — malate conversion (D-23) + citrate → diacetyl (D-31).

This module holds the *Oenococcus oeni* Processes. :class:`MalolacticConversion` (D-23) is the
malate → lactate + CO2 deacidification. :class:`MalolacticCitrateMetabolism` and
:class:`OenococcusDiacetylReduction` (D-31) add **MLF-derived diacetyl**: the bacterium
co-metabolises citric acid into α-acetolactate (feeding the shared VDK reservoir, so diacetyl
emerges from the existing D-26 machinery) and also reduces diacetyl on the lees. All three share
the *O. oeni* environmental gate (:func:`malolactic_environmental_gate`) and the compile-seam
isolability: they are disabled unless *O. oeni* is pitched, so an un-pitched wine run is
byte-for-byte the validated core and the ``malic``/``lactic``/``citrate`` slots keep their
VALIDATED tier. See :class:`MalolacticCitrateMetabolism` for the citrate-pool rationale (its
carbon must survive past sugar-dryness) and the lumped carbon-closing stoichiometry.

*Oenococcus oeni* converts L-malic acid (C4, diprotic) to L-lactic acid (C3,
monoprotic) plus CO2, mole-for-mole. That single reaction deacidifies the wine
(pH rises ~0.1–0.3, the D-18 headline coupling) and softens its perceived acidity.
This module is the first **RHS consumer** of the pH charge-balance keystone (D-18)
and of the molecular-SO₂ readout (D-22): the conversion rate is gated by the *solved*
pH, by molecular (antimicrobial) SO₂, by ethanol, and by a temperature optimum — so
the deacidification feedback (pH ↑ ⇒ rate ↑, self-limited as malate depletes) and the
SO₂/ethanol arrest of MLF *emerge* from the model rather than being scripted.

**Conversion was built first (v1, decision D-23); MLF-growth is now landed as a clean
extension.** *O. oeni* builds biomass mostly from amino acids/peptides, but the lumped
``N`` (YAN) is carbon-free in :func:`total_carbon` (D-19) and is driven to ~0 within
~1.3 d of the AF pitch *regardless of dose* (the empirical finding that settles D-23), so
there is no free nitrogen at the MLF pitch point to fund bacterial growth. Modelling
MLF-growth honestly therefore needed a separate amino-acid ledger (D-32) *and* an
autolytic-peptide refill source (D-34) — both since landed. So :class:`Malolactic
Conversion` treats the bacterium as a **catalyst that scales the conversion rate**, and
:class:`MalolacticGrowth` (the deferred growth beat, now built) makes ``X_mlf`` *dynamic*:
it consumes the ``amino_acids`` pool to build bacterial biomass, which — since conversion
is linear in ``X_mlf`` — *accelerates* deacidification autocatalytically. Growth was a
pure add-a-Process extension (no refactor of the conversion), exactly as v1 promised.

**X_mlf is promoted from an inert carbon-/nitrogen-free catalyst to real biomass.** With
:class:`MalolacticGrowth` active, ``X_mlf`` gains carbon and nitrogen, so it is now
weighted in ``total_carbon``/``total_nitrogen`` at the biomass fractions (bacterial
elemental composition approximated by the yeast's, a documented v1 simplification — the
same fractions the growth stoichiometry draws against, so conservation closes exactly).
A co-inoculation dose or a ``pitch_mlf`` intervention therefore now adds bacterial-biomass
carbon/nitrogen to the run (booked as an external flow for the intervention), superseding
the v1 "dosing X_mlf leaves total_carbon byte-for-byte" claim (decision D-38).

**Carbon (and mass) close on the existing ledger — no new conservation code.** With
``r`` [mol/L/h] the malate turnover,

    d(malic)/dt  = −r · M_malic
    d(lactic)/dt = +r · M_lactic
    d(CO2)/dt    = +r · M_CO2

and since malic (4 C) = lactic (3 C) + CO2 (1 C), carbon closes mole-for-mole; mass
closes too (134.087 = 90.078 + 44.009 g/mol — a clean decarboxylation, no water term).
:func:`total_carbon` has weighted ``malic``/``lactic``/``CO2`` since D-18 precisely so
this conversion is carbon-closing, so the Process touches only those three slots and
adds nothing to the conservation harness. (``total_mass`` does *not* count ``malic``, so
it is not a valid check on a dosed-MLF run — use ``total_carbon``, the same scoping as
the glycerol-on caveat, D-16.)

**Rate form — substrate-limited, catalyst-scaled, multiplicatively gated.**

    r = k_mlf · X_mlf · [malate]/(K_mlf + [malate]) · g_pH · g_EtOH · g_SO₂ · γ(T)

* ``k_mlf · X_mlf`` — specific malolactic activity times the (constant) bacterial
  concentration. Linear in ``X_mlf`` so the dose sets the timescale and the growth beat
  drops in by making ``X_mlf`` dynamic.
* **Michaelis–Menten in malate** (``K_mlf`` a low half-saturation, ~mM): the malolactic
  enzyme has high malate affinity, so the term stays ≈1 until malate is nearly exhausted,
  giving a clean near-complete conversion rather than a long tail.
* **pH gate** ``g_pH = 1/(1 + 10^(pH_half − pH))`` — a smooth logistic rising with pH
  (≈0.5 at ``pH_half_mlf``). MLF is strongly inhibited below ~pH 3.0. As malate→lactate
  raises pH this gate *rises*: realistic self-reinforcing feedback, bounded above by 1 and
  self-limited by malate depletion, and C∞ in pH (good for the BDF solver).
* **ethanol gate** ``g_EtOH = max(0, 1 − E/E_max)^n`` — the Luong wall reused from
  :class:`~fermentation.core.kinetics.inhibition.EthanolInhibition` (``E_max =
  ethanol_tolerance_mlf`` ≈ *O. oeni*'s ~13–14 % ABV tolerance, below the yeast's). v1 is
  **co-inoculation** (bacteria from t=0), so ethanol is low while malate converts and
  this gate is ≈1 early — the regime where co-inoculated MLF actually runs.
* **molecular-SO₂ gate** ``g_SO₂ = exp(−[SO₂]_molecular / s)`` — molecular SO₂ is the
  antimicrobial species (D-22), partitioned at the *same* solved pH (computed once, see
  below). Exponential so it is smooth, in (0, 1], and ≈1 when undosed.
* **temperature** ``γ(T)`` — a **cardinal-temperature optimum** (:func:`cardinal_\
  temperature_factor`, Rosso et al. 1993), peaking at ``T_opt_mlf`` and falling to 0 at
  ``T_min_mlf``/``T_max_mlf``. Unlike the monotone-Arrhenius byproduct Processes, MLF has
  a genuine optimum and *declines* in the warm — an Arrhenius factor would be
  qualitatively wrong above ~25 °C (decision D-23 explicitly scoped "a temperature
  optimum"). All terms are temperature *differences*, so the form is identical in K or °C
  (the cardinals are stored in K for canonical consistency, sourced in °C).

**Performance / isolability — guard before the pH solve.** Reading the solved pH costs a
``brentq`` per RHS evaluation, so the Process returns a zero contribution *before* solving
pH whenever ``X_mlf ≤ 0`` (undosed) or there is no malate left. That keeps every undosed
wine run (and the §2.2 wine benchmark) from paying the solve for nothing, and makes the
undosed contribution byte-for-byte zero. Value-isolability is thus structural; **tier**
isolability is handled at the compile seam, which *disables* this Process when MLF is not
pitched so the inert ``malic``/``lactic`` slots keep their VALIDATED tier (nothing active
touches them) rather than being dragged to speculative by an enabled-but-zero Process
(decision D-23; ``ProcessSet.tier_of`` counts enabled processes, not nonzero ones).

Tier: **speculative**. The conversion stoichiometry and the *directions* of the gates are
sound, but the rate magnitude ``k_mlf`` and every gate/cardinal constant are
order-of-magnitude estimates (no kinetic model of this exact per-catalyst flux form is
sourced), and v1 knowingly omits bacterial growth/lag. Parameter-tier propagation (D-1)
caps the ``malic``/``lactic``/``CO2`` outputs at speculative regardless.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

from fermentation.core.acidbase import (
    SO2_STATE_KEY,
    molecular_so2_at_ph,
    ph_of_state,
)
from fermentation.core.chemistry import (
    M_ACETOLACTATE,
    M_BUTANEDIOL,
    M_CITRIC,
    M_CO2,
    M_DIACETYL,
    M_LACTIC,
    M_MALIC,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import draw_carbon_from_sugar
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


def cardinal_temperature_factor(temp: float, t_min: float, t_opt: float, t_max: float) -> float:
    """Cardinal-temperature growth factor — Rosso et al. (1993) CTMI, peak 1 at ``t_opt``.

    The standard predictive-microbiology temperature response: a smooth unimodal curve
    that is 0 at and outside the cardinal bounds ``[t_min, t_max]`` and exactly 1 at the
    optimum ``t_opt`` (so the rate constant is read unscaled at the optimum, like the
    Arrhenius factor at ``T_ref``). The functional form is

        γ(T) = (T − T_max)(T − T_min)² /
               { (T_opt − T_min)·[ (T_opt − T_min)(T − T_opt)
                                   − (T_opt − T_max)(T_opt + T_min − 2T) ] }

    Every term is a temperature *difference*, so the value is identical whether ``temp``
    and the cardinals are all in K or all in °C — we pass Kelvin (canonical units, D-3).
    Returns 0 outside ``(t_min, t_max)`` (no extrapolation past the cardinals), giving the
    biologically-correct decline in the warm that a monotone Arrhenius factor cannot.

    Reference: Rosso, Lobry & Flandrois (1993), J. Theor. Biol. 162:447–463 — the
    cardinal-temperature model with inflection (CTMI), the canonical bacterial form.
    """
    if temp <= t_min or temp >= t_max:
        return 0.0
    numerator = (temp - t_max) * (temp - t_min) ** 2
    denominator = (t_opt - t_min) * (
        (t_opt - t_min) * (temp - t_opt) - (t_opt - t_max) * (t_opt + t_min - 2.0 * temp)
    )
    return float(numerator / denominator)


#: The *O. oeni* environmental-gate parameters, shared by every MLF Process (the malate
#: conversion and the citrate → diacetyl branch, D-31). Declared once so the two Processes'
#: ``reads`` tuples and :func:`malolactic_environmental_gate` cannot drift apart.
_MLF_GATE_READS: tuple[str, ...] = (
    "pH_half_mlf",
    "ethanol_tolerance_mlf",
    "mlf_ethanol_exponent",
    "molecular_so2_inhib_mlf",
    "T_min_mlf",
    "T_opt_mlf",
    "T_max_mlf",
)


def malolactic_environmental_gate(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], ph: float
) -> float:
    """The shared *O. oeni* environmental gate ``g_pH · g_EtOH · g_SO₂ · γ(T)`` ∈ [0, 1].

    Every *O. oeni* activity — malate conversion (D-23) and citrate → diacetyl (D-31) — is
    throttled by the *same* environment, so both Processes multiply their rate by this one
    factor (a shared helper, decision D-31): the pH logistic (rises with pH), the Luong
    ethanol wall, the molecular-SO₂ exponential, and the Rosso cardinal-temperature optimum.
    The caller passes the *already-solved* ``ph`` (each Process solves it once via
    :func:`ph_of_state`, after its cheap ``X_mlf``/substrate guards) so this helper never
    triggers a second ``brentq``; molecular SO₂ is partitioned at that same pH. See the
    module docstring for each term's form and sourcing.
    """
    gate_ph = 1.0 / (1.0 + 10.0 ** (params["pH_half_mlf"] - ph))

    total_so2 = float(y[schema.slice(SO2_STATE_KEY)][0]) if SO2_STATE_KEY in schema else 0.0
    if total_so2 > 0.0:
        molecular_so2 = molecular_so2_at_ph(y, schema, params, ph)
        gate_so2 = math.exp(-molecular_so2 / params["molecular_so2_inhib_mlf"])
    else:
        gate_so2 = 1.0

    e = max(float(y[schema.slice("E")][0]), 0.0)
    remaining = 1.0 - e / params["ethanol_tolerance_mlf"]
    gate_eth = remaining ** params["mlf_ethanol_exponent"] if remaining > 0.0 else 0.0

    temp = float(y[schema.slice("T")][0])
    gamma_t = cardinal_temperature_factor(
        temp, params["T_min_mlf"], params["T_opt_mlf"], params["T_max_mlf"]
    )
    return float(gate_ph * gate_eth * gate_so2 * gamma_t)


class MalolacticConversion(Process):
    """Malolactic fermentation v1 — malate → lactate + CO2, catalyst-scaled and gated.

    ``d(malic)/dt = −r·M_malic``, ``d(lactic)/dt = +r·M_lactic``, ``d(CO2)/dt = +r·M_CO2``
    with the molar turnover ``r = k_mlf·X_mlf·[malate]/(K_mlf+[malate])·g_pH·g_EtOH·g_SO₂·
    γ(T)`` (see the module docstring). Carbon and mass close on the existing ledger
    (4 C = 3 C + 1 C). Touches only ``malic``/``lactic``/``CO2``; reads the (constant)
    catalyst ``X_mlf``, the solved pH, molecular SO₂, ethanol ``E`` and ``T`` from state.

    Returns a zero contribution before solving pH when undosed (``X_mlf ≤ 0``) or when
    malate is exhausted — structural value-isolability and no wasted ``brentq`` (the
    compile seam additionally disables the Process when MLF is not pitched, for tier
    isolability; decision D-23). Tier **speculative** (rate/gate magnitudes are estimates).
    """

    name = "malolactic_conversion"
    tier = Tier.SPECULATIVE
    touches = ("malic", "lactic", "CO2")
    #: Kinetic + gate parameters (all wine, *O. oeni*; see wine_generic.yaml). Their tiers
    #: cap the malic/lactic/CO2 output tier via parameter-tier propagation (D-1). The pH /
    #: SO₂ pKa parameters the gates read through ``acidbase`` are NOT listed (they are all
    #: plausible, and this Process is speculative, so the combined floor is unchanged —
    #: matching the convention that pKa dependence is reported via ``acidbase.ph_tier``).
    reads: tuple[str, ...] = ("k_mlf", "K_mlf", *_MLF_GATE_READS)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE the pH brentq solve (D-23): no catalyst or no substrate ⇒ no
        # conversion, and undosed runs must not pay a per-RHS pH solve for a zero result.
        x_mlf = max(float(y[schema.slice("X_mlf")][0]), 0.0) if "X_mlf" in schema else 0.0
        if x_mlf <= 0.0:
            return d
        malic_gpl = max(float(y[schema.slice("malic")][0]), 0.0)
        if malic_gpl <= 0.0:
            return d

        # pH is solved once and reused by both the pH gate and the SO₂ partition (D-22/D-28):
        # the shared gate takes the solved pH so acidbase.speciate_so2 needs no second brentq.
        # The gate weakens with the early acetaldehyde peak (which sequesters free SO₂, D-28).
        ph = ph_of_state(y, schema, params)
        gate = malolactic_environmental_gate(y, schema, params, ph)

        malate_molar = malic_gpl / M_MALIC
        monod = malate_molar / (params["K_mlf"] + malate_molar)
        r = params["k_mlf"] * x_mlf * monod * gate

        d[schema.slice("malic")] = -r * M_MALIC
        d[schema.slice("lactic")] = r * M_LACTIC
        d[schema.slice("CO2")] = r * M_CO2
        return d


class MalolacticCitrateMetabolism(Process):
    """MLF-derived diacetyl (decision D-31) — *O. oeni* citrate → α-acetolactate + CO2.

    The real coupling MLF unlocks: alongside malate, *Oenococcus oeni* co-metabolises
    **citric acid**, overflowing α-acetolactate that (non-enzymatically) decarboxylates to
    **diacetyl** — the buttery note post-MLF wines carry. This Process is the citrate source;
    the diacetyl then *emerges* from the existing VDK machinery (decision D-26), which needs no
    change: the always-on :class:`~fermentation.core.kinetics.vicinal_diketones.\
    AcetolactateDecarboxylation` converts the reservoir → diacetyl, and both the yeast
    :class:`~fermentation.core.kinetics.vicinal_diketones.DiacetylReduction` (viable ``X``) and
    the bacterial :class:`OenococcusDiacetylReduction` (``X_mlf``) clear it. Feeding the shared
    α-acetolactate reservoir is the *genuine* topology, not merely DRY: the reservoir → diacetyl
    step is a spontaneous property of the identical molecule regardless of whether the
    α-acetolactate came from yeast valine overflow or bacterial citrate metabolism.

    **Why a citrate pool at all (the load-bearing scope decision, D-31).** MLF-diacetyl is a
    late-MLF, often *post-dryness* phenomenon, so its carbon cannot come from sugar: the yeast
    VDK stand-in draws α-acetolactate carbon out of ``S`` (:func:`~fermentation.core.kinetics.\
    carbon_routing.draw_carbon_from_sugar`), which correctly no-ops once ``S`` is exhausted — it
    would either strand carbon (breaking ``total_carbon`` closure) or stop diacetyl production
    exactly when this beat needs it. Citrate is present independent of sugar, so it is the
    honest source (a dosed must input, ``citrate`` slot).

    **Stoichiometry — a lumped, carbon-closing stand-in (own the fiction).** With ``r_c``
    [mol/L/h] the citrate turnover, ``d(citrate)/dt = −r_c·M_citric``,
    ``d(acetolactate)/dt = +r_c·M_acetolactate``, ``d(CO2)/dt = +r_c·M_CO2``. Citric acid (6 C)
    → α-acetolactate (5 C) + CO2 (1 C), so **carbon closes mole-for-mole (6 = 5 + 1)** on the
    existing ledger, exactly like malic → lactic + CO2 (D-23). *Mass* carries a small gap
    (192.124 ≠ 132.116 + 44.009; the real pathway's untracked acetate/H₂O/redox), so carbon is
    the invariant, as for beer's hydrolysis water (D-8) and the VDK decarb. CAVEAT: real citrate
    metabolism is ``citrate → acetate + oxaloacetate → pyruvate + CO2``, ~2 citrate per
    α-acetolactate, with **acetate** (a volatile-acidity contributor) the *dominant* co-product.
    The single-reaction stand-in drops the acetate/lactate branches; the rate ``k_citrate`` is
    held low so citrate stays **mostly unconsumed** — the *trace diacetyl branch only* — which
    keeps the fiction honest (we do not claim to resolve citrate's full fate).

    **Rate — citrate's own Monod × the shared MLF environmental gate (NOT malate's ``r``).**
    ``r_c = k_citrate · X_mlf · [citrate]/(K_citrate + [citrate]) · g_pH·g_EtOH·g_SO₂·γ(T)``.
    Coupling to citrate (not the malate turnover) is deliberate: malate's rate → 0 at malate
    depletion, which would kill exactly the *post-malate* diacetyl this pool exists to capture.
    Citrate is metabolised while bacteria are active and depletes on its own timescale, so
    diacetyl accumulates through and past malate conversion, then falls as reduction clears it —
    the realistic late peak. The environmental gate is the *same* one MLF conversion uses
    (:func:`malolactic_environmental_gate`), so SO₂/ethanol/low-pH arrest citrate metabolism
    just as they arrest MLF. Guards mirror MLF conversion: a zero contribution *before* the pH
    ``brentq`` when undosed (``X_mlf ≤ 0``) or citrate-exhausted (value + perf isolability; the
    compile seam additionally disables the Process when MLF is not pitched, for tier
    isolability, so the dosed ``citrate`` slot keeps its VALIDATED tier undosed — like
    ``malic``/``lactic``). Tier **speculative** (rate/gate magnitudes are estimates).
    """

    name = "malolactic_citrate_metabolism"
    tier = Tier.SPECULATIVE
    touches = ("citrate", "acetolactate", "CO2")
    #: The citrate Monod pair plus the shared *O. oeni* environmental-gate parameters. Their
    #: tiers cap the citrate/acetolactate/CO2 output tier via parameter-tier propagation (D-1);
    #: ``acetolactate``/``CO2`` are already speculative (the always-on yeast VDK pathway), so
    #: this adds no new tier headline. pKa/SO₂ params read via ``acidbase`` are omitted for the
    #: same reason as MalolacticConversion (all plausible; the Process is already speculative).
    reads: tuple[str, ...] = ("k_citrate", "K_citrate", *_MLF_GATE_READS)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE the pH brentq solve (mirroring MalolacticConversion): no catalyst
        # or no citrate ⇒ no diacetyl branch, and undosed runs must not pay a per-RHS pH solve.
        x_mlf = max(float(y[schema.slice("X_mlf")][0]), 0.0) if "X_mlf" in schema else 0.0
        if x_mlf <= 0.0:
            return d
        citrate_gpl = max(float(y[schema.slice("citrate")][0]), 0.0) if "citrate" in schema else 0.0
        if citrate_gpl <= 0.0:
            return d

        ph = ph_of_state(y, schema, params)
        gate = malolactic_environmental_gate(y, schema, params, ph)

        citrate_molar = citrate_gpl / M_CITRIC
        monod = citrate_molar / (params["K_citrate"] + citrate_molar)
        r_c = params["k_citrate"] * x_mlf * monod * gate  # citrate turnover, mol/L/h

        d[schema.slice("citrate")] = -r_c * M_CITRIC
        d[schema.slice("acetolactate")] = r_c * M_ACETOLACTATE  # feeds the shared VDK reservoir
        d[schema.slice("CO2")] = (
            r_c * M_CO2
        )  # citrate C6 → acetolactate C5 + CO2 C1 (carbon-closing)
        return d


class OenococcusDiacetylReduction(Process):
    """Bacterial diacetyl → 2,3-butanediol by *O. oeni* — lees-contact clearing (D-31).

    ``d(diacetyl)/dt = −L``, ``d(butanediol)/dt = +L·M_butanediol/M_diacetyl`` with the mass
    loss ``L = k_mlf_diacetyl_reduction · X_mlf · f(T) · [diacetyl]`` and ``f(T) =
    arrhenius_factor(T, E_a_reduction, T_ref)``. A mole-for-mole C4 → C4 transfer (both weighted
    at their own carbon fraction), so it is carbon-neutral like the yeast reduction (D-26); no
    sugar draw. Acetoin is lumped into ``butanediol``.

    **Why add it (the owner's D-31 call).** *O. oeni* reduces diacetyl too, so leaving wine on
    the lees with active bacteria lowers diacetyl — the real reason a completed MLF left on lees
    cleans up, and the reason an early post-MLF SO₂ addition (which would kill the bacteria)
    *locks in* diacetyl. It complements the yeast :class:`~fermentation.core.kinetics.\
    vicinal_diketones.DiacetylReduction`: in co-inoculation the yeast is still viable and clears
    diacetyl fast; this bacterial reducer keeps clearing it after the yeast is
    ethanol-inactivated, as long as *O. oeni* is present.

    **Gated on ``X_mlf`` and temperature only (v1 simplification).** Like the MLF conversion
    catalyst, ``X_mlf`` is a constant, inert dose in v1 (no Process grows or kills it), so this
    reduction does not carry the environmental (ethanol/SO₂/pH) arrest gates — bacterial *death*
    is deferred with the MLF-growth beat. Consequence to flag: with *O. oeni* dosed, MLF-diacetyl
    is **not** permanently stranded (the realistic lees-contact clearing); the "package/rack
    early ⇒ diacetyl locked in" case needs a racking event to remove ``X_mlf``, deferred to the
    event loop (decision D-23/D-31). ``diacetyl`` and ``X_mlf`` are clamped ≥ 0 against solver
    undershoot. Tier **speculative** (rate magnitude estimate).
    """

    name = "oenococcus_diacetyl_reduction"
    tier = Tier.SPECULATIVE
    touches = ("diacetyl", "butanediol")
    #: ``k_mlf_diacetyl_reduction`` sets the bacterial reductase magnitude; ``E_a_reduction``
    #: (shared with the yeast reduction — a generic oxidoreductase activation energy) and
    #: ``T_ref`` set the temperature shape. Reads the constant catalyst ``X_mlf`` from state.
    reads: tuple[str, ...] = ("k_mlf_diacetyl_reduction", "E_a_reduction", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        diacetyl = max(float(y[schema.slice("diacetyl")][0]), 0.0)
        if diacetyl <= 0.0:  # nothing to reduce
            return d
        x_mlf = max(float(y[schema.slice("X_mlf")][0]), 0.0) if "X_mlf" in schema else 0.0
        if x_mlf <= 0.0:  # no bacteria ⇒ no bacterial reduction
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_reduction"], params["T_ref"])
        loss = params["k_mlf_diacetyl_reduction"] * x_mlf * f_t * diacetyl
        d[schema.slice("diacetyl")] = -loss
        d[schema.slice("butanediol")] = loss * M_BUTANEDIOL / M_DIACETYL  # mole-for-mole C4→C4
        return d


class MalolacticGrowth(Process):
    """*Oenococcus oeni* biomass growth — the deferred MLF-growth beat (decision D-23).

    Makes the ``X_mlf`` catalyst *dynamic*: it grows on the assimilable ``amino_acids`` pool
    (D-32, refilled by autolysis D-34), and because :class:`MalolacticConversion` is linear in
    ``X_mlf`` the deacidification then **accelerates autocatalytically** as bacteria multiply —
    the realistic MLF the constant-``X_mlf`` v1 could not produce. This is the pure
    add-a-Process extension v1 promised: no change to the conversion kinetics.

    **Rate — the bacterial growth law, environmentally gated.**

        dX_mlf/dt = μ_max_mlf · X_mlf · aa/(K_aa_mlf + aa) · S/(K_s + S) · g_pH·g_EtOH·g_SO₂·γ(T)

    Michaelis–Menten in the amino-acid pool (the nitrogen/carbon fuel) *and* in total sugar
    (the energy source O. oeni ferments during co-fermentation, so growth self-limits to the
    sugar-present window and never draws its carbon shortfall from an empty ``S`` — see below),
    scaled by the *same* :func:`malolactic_environmental_gate` the conversion uses. The ethanol
    wall in that gate is what makes co-inoculation the *dominant* MLF-growth mode **emergently**:
    a post-AF pitch into a high-ABV must lands past the O. oeni ethanol tolerance, so γ·g_EtOH ≈ 0
    and bacteria cannot build up — while a normal-ABV sequential MLF, where g_EtOH is small but
    nonzero, still grows. This is left to the gate, not hard-coded (the compile seam gates only on
    amino-acid fuel, mirroring how conversion trusts its ethanol gate rather than a pitch rule;
    decision D-38).

    **Conservation — nitrogen-anchored, carbon shortfall from sugar (decision D-38).**
    New bacterial biomass needs ``f_N·dX_mlf`` nitrogen and ``f_C·dX_mlf`` carbon (``f_N``/``f_C``
    the biomass fractions ``X_mlf`` is weighted at, so this closes exactly). All the nitrogen is
    taken from amino acids, consuming ``ρ = f_N·dX_mlf/y_N`` of the pool (``y_N``/``y_C`` =
    arginine's nitrogen/carbon mass fractions). That arginine carries only ``ρ·y_C`` carbon —
    less than biomass needs, because arginine (mass C:N ≈ 1.29) is far more nitrogen-rich than
    biomass (C:N ≈ 4–11) — so the **shortfall** ``f_C·dX_mlf − ρ·y_C = dX_mlf·(f_C − f_N·y_C/y_N)``
    is drawn from sugar (:func:`~fermentation.core.kinetics.carbon_routing.draw_carbon_from_\
    sugar`). This is the mirror image of yeast growth (nitrogen from a nitrogen pool, carbon from
    sugar; :class:`~fermentation.core.kinetics.growth.GrowthNitrogenLimited`) and of D-34 autolysis
    (which routes the *excess* carbon of the reverse split to debris). The shortfall coefficient
    ``f_C − f_N·y_C/y_N`` is **structurally positive** across Coleman's whole ``f_N`` range (biomass
    C:N always exceeds arginine's), so no clamp and no C⁰ kink for the BDF solver. Carbon closes
    (``X_mlf`` gains ``f_C·dX_mlf`` = amino-acid carbon ``ρ·y_C`` + sugar shortfall); nitrogen
    closes (``X_mlf`` gains ``f_N·dX_mlf`` = amino-acid nitrogen ``ρ·y_N``). Touches
    ``(X_mlf, amino_acids, S)`` — notably **not** ``N``: unlike the C-anchored alternative it
    releases no artificial ammonium (the anchoring fork, decision D-38).

    Guards mirror the conversion: a zero contribution *before* the pH ``brentq`` when there is no
    catalyst (``X_mlf ≤ 0``), no amino-acid fuel, or no sugar (so the shortfall never targets an
    empty ``S``). Tier **speculative** (``μ_max_mlf``/``K_aa_mlf`` are author estimates and the
    bacterial-≈-yeast composition is a simplification).
    """

    name = "malolactic_growth"
    tier = Tier.SPECULATIVE
    #: Builds bacterial biomass ``X_mlf`` from the ``amino_acids`` pool, drawing the carbon
    #: shortfall from ``S``. Does NOT touch ``N`` (nitrogen is sourced from amino acids, not the
    #: ammonium pool — the D-38 anchoring choice).
    touches = ("X_mlf", "amino_acids", "S")
    #: ``mu_max_mlf``/``K_aa_mlf`` are the bacterial growth rate + amino-acid half-saturation;
    #: ``K_s`` (reused from yeast growth) sets the sugar Monod; ``biomass_N_fraction``/
    #: ``biomass_C_fraction`` are the composition the biomass is built at (and weighted at in the
    #: conservation checks — the same fractions, so carbon/nitrogen close). The shared MLF
    #: environmental-gate parameters throttle growth by pH/ethanol/SO₂/temperature. Their tiers
    #: cap the ``X_mlf``/``amino_acids``/``S`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "mu_max_mlf",
        "K_aa_mlf",
        "K_s",
        "biomass_N_fraction",
        "biomass_C_fraction",
        *_MLF_GATE_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE the pH brentq solve (mirroring MalolacticConversion): no catalyst,
        # no amino-acid fuel, or no sugar (the carbon shortfall must never target an empty S) ⇒
        # no growth, and no wasted per-RHS pH solve on an undosed/unpitched run.
        x_mlf = max(float(y[schema.slice("X_mlf")][0]), 0.0) if "X_mlf" in schema else 0.0
        if x_mlf <= 0.0:
            return d
        aa = max(float(y[schema.slice("amino_acids")][0]), 0.0) if "amino_acids" in schema else 0.0
        if aa <= 0.0:
            return d
        s_total = max(float(y[schema.slice("S")].sum()), 0.0)
        if s_total <= 0.0:
            return d

        ph = ph_of_state(y, schema, params)
        gate = malolactic_environmental_gate(y, schema, params, ph)

        aa_monod = aa / (params["K_aa_mlf"] + aa)
        s_monod = s_total / (params["K_s"] + s_total)
        dx_mlf = params["mu_max_mlf"] * x_mlf * aa_monod * s_monod * gate  # [g X_mlf/L/h]
        if dx_mlf <= 0.0:
            return d

        f_n = params["biomass_N_fraction"]
        f_c = params["biomass_C_fraction"]
        y_n = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
        y_c = carbon_mass_fraction(AMINO_ACID_SPECIES)
        # Nitrogen-anchored: consume the arginine that carries the new biomass's nitrogen. Its
        # carbon (rho*y_C) falls short of the biomass carbon demand (f_C*dx_mlf) because arginine
        # is N-rich, so the positive shortfall is drawn from sugar. Both ledgers close exactly:
        # X_mlf gains f_N*dx_mlf N (= rho*y_N, the amino-acid N) and f_C*dx_mlf C (= rho*y_C amino
        # acid + shortfall sugar). The shortfall coefficient (f_C - f_N*y_C/y_N) > 0 structurally
        # (biomass C:N >> arginine's), so no clamp is needed (decision D-38).
        rho = f_n * dx_mlf / y_n  # [g arginine/L/h] consumed to supply the biomass nitrogen
        d[schema.slice("X_mlf")] = dx_mlf
        d[schema.slice("amino_acids")] = -rho
        shortfall = f_c * dx_mlf - rho * y_c  # [g C/L/h] biomass carbon not covered by arginine
        draw_carbon_from_sugar(d, y, schema, shortfall)
        return d
