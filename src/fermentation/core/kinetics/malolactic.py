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


#: The *toxicity* half of the *O. oeni* environmental gate: the pH logistic, the Luong ethanol
#: wall, and the molecular-SO₂ exponential — every term that goes to 0 under an *adverse* chemical
#: environment (low pH, high ethanol, high SO₂). Split out from the temperature cardinals (D-39)
#: because :class:`MalolacticDeath` is driven by ``1 − toxicity`` (stress kills bacteria) but must
#: NOT reuse the cardinal γ(T): γ(T)→0 in the *cold*, which would make cold *kill* bacteria when it
#: in fact preserves them — death carries its own Arrhenius factor instead (see that class).
_MLF_TOXICITY_READS: tuple[str, ...] = (
    "pH_half_mlf",
    "ethanol_tolerance_mlf",
    "mlf_ethanol_exponent",
    "molecular_so2_inhib_mlf",
)
#: The cardinal-temperature parameters of the growth/conversion gate (Rosso optimum γ(T)).
_MLF_TEMP_READS: tuple[str, ...] = ("T_min_mlf", "T_opt_mlf", "T_max_mlf")
#: The full *O. oeni* environmental-gate parameters, shared by the growth/conversion Processes (the
#: malate conversion and the citrate → diacetyl branch, D-31). Declared once so those Processes'
#: ``reads`` tuples and :func:`malolactic_environmental_gate` cannot drift apart. Order preserved
#: from before the D-39 toxicity/temperature split so every consumer's ``reads`` is byte-identical.
_MLF_GATE_READS: tuple[str, ...] = (*_MLF_TOXICITY_READS, *_MLF_TEMP_READS)


def malolactic_toxicity_gate(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], ph: float
) -> float:
    """The *chemical-toxicity* factor ``g_pH · g_EtOH · g_SO₂`` ∈ (0, 1] (decision D-39).

    The temperature-independent half of the *O. oeni* environmental gate: the pH logistic (rises
    with pH), the Luong ethanol wall, and the molecular-SO₂ exponential. It is ≈1 in a benign
    environment (near-neutral pH, low ethanol, no SO₂) and →0 as any stressor bites. Growth and
    conversion multiply it by the cardinal γ(T) (:func:`malolactic_environmental_gate`); bacterial
    *death* keys off ``1 − toxicity`` — high ethanol / low pH / SO₂ raise the death rate — with its
    own Arrhenius temperature factor rather than γ(T) (which would spuriously kill in the cold).
    The caller passes the *already-solved* ``ph`` so molecular SO₂ is partitioned without a second
    ``brentq``. See the module docstring for each term's form and sourcing.
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

    return float(gate_ph * gate_eth * gate_so2)


def malolactic_environmental_gate(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], ph: float
) -> float:
    """The shared *O. oeni* environmental gate ``g_pH · g_EtOH · g_SO₂ · γ(T)`` ∈ [0, 1].

    Every *O. oeni* growth/conversion activity — malate conversion (D-23), citrate → diacetyl
    (D-31) and bacterial growth (D-38) — is throttled by the *same* environment, so they multiply
    their rate by this one factor (a shared helper, decision D-31): the :func:`malolactic_\
    toxicity_gate` chemistry (pH · ethanol · SO₂) times the Rosso cardinal-temperature optimum
    γ(T). The caller passes the *already-solved* ``ph`` (each Process solves it once via
    :func:`ph_of_state`, after its cheap ``X_mlf``/substrate guards) so this helper never
    triggers a second ``brentq``; molecular SO₂ is partitioned at that same pH. See the
    module docstring for each term's form and sourcing.

    Since the D-39 split this is ``toxicity · γ(T)`` with the multiplication grouped exactly as
    before (``((g_pH·g_EtOH)·g_SO₂)·γ(T)``), so the three growth/conversion consumers' rates are
    byte-for-byte unchanged.
    """
    toxicity = malolactic_toxicity_gate(y, schema, params, ph)
    temp = float(y[schema.slice("T")][0])
    gamma_t = cardinal_temperature_factor(
        temp, params["T_min_mlf"], params["T_opt_mlf"], params["T_max_mlf"]
    )
    return float(toxicity * gamma_t)


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


class MalolacticDeath(Process):
    """*Oenococcus oeni* death — the **SO₂-driven** bacterial kill (decision D-39).

    The counterpart to :class:`MalolacticGrowth` that completes the MLF arc (D-23 → D-31 → D-38):
    it moves viable ``X_mlf`` into the non-viable ``X_mlf_dead`` pool, so bacterial biomass
    *declines* when the wine is sulfited and the *O. oeni* activities that scale with ``X_mlf`` —
    malate conversion (D-23), citrate → diacetyl and, above all,
    :class:`OenococcusDiacetylReduction` (D-31) — wind down as the bacteria die. That is the
    mechanism behind the winemaking headline the D-31 reducer flagged as deferred: an early
    post-MLF **SO₂ addition kills the bacteria**, so diacetyl stops being cleared on the lees and is
    **locked in** (a rack that draws the bacteria off the lees does the same physically — decision
    D-39, the ``rack`` extension).

    **Rate — molecular-SO₂-driven, Arrhenius temperature.**

        r_death = k_death_mlf · X_mlf · (1 − g_SO₂) · arrhenius(T, E_a_death_mlf, T_ref)
        g_SO₂   = exp(−[SO₂]_molecular / molecular_so2_inhib_mlf)

    The driver is ``1 − g_SO₂`` — the *molecular-SO₂* term alone, the **same** ``g_SO₂`` the
    conversion/growth gate uses (D-22 antimicrobial readout, partitioned at the solved pH). Death
    is **exactly 0 without SO₂** (``g_SO₂ = 1``) and rises toward its Arrhenius ceiling as molecular
    SO₂ accumulates. Temperature enters through its own Arrhenius factor (enzymatic/chemical
    mortality, faster when warm — the same ``E_a``/``T_ref`` shape autolysis uses), **not** the
    cardinal γ(T): γ(T) →0 in the cold, which would make cold *kill* bacteria, whereas cold in fact
    preserves them. So warm accelerates the SO₂ kill and cold slows it to dormancy (decision D-39).

    **Why SO₂ alone, not the full toxicity gate (the D-39 crux, empirically settled).** An earlier
    draft drove death by ``1 − toxicity`` (pH·ethanol·SO₂). A probe killed it: the Luong ethanol
    wall already drives ``1 − toxicity`` to ~0.92 at ordinary post-AF ethanol (~75 g/L), so death
    was near-maximal *from ethanol alone* — bacteria died in ~1 week, when in reality *O. oeni*
    persists for weeks-to-months in dry wine and is cleared deliberately by **SO₂ (or racking)**.
    Ethanol's wall is a "can't grow" signal, not a "dying" signal; coupling death to it was the bug.
    No power transform of 0.92 can be both ~0 (slow baseline) and clearly below the SO₂-elevated
    value, so the fix is a driver that contains **no ethanol term** — molecular SO₂ only (decision
    D-39). The co-inoculation-vs-post-AF-pitch dominance now rests entirely on the *growth* gate's
    ``g_EtOH`` (high-ABV post-AF ⇒ γ·g_EtOH ≈ 0 ⇒ no growth), which is where it belongs — a
    high-ABV post-AF pitch simply sits **inert** (no growth, no conversion, no death until SO₂).

    **Scope: SO₂-driven acute kill only — the slow baseline decline is now :class:`Malolactic
    Senescence` (v2, D-41).** This Process is the *deliberate-action* lever: death is **exactly 0
    without SO₂**, so it stays unconfounded by ethanol (the D-39 crux) and keeps ``k_death_mlf`` at
    its true SO₂-kill magnitude. The slow ethanol/age decline of *O. oeni* over weeks-to-months —
    which v1 deferred, letting unsulfited bacteria persist forever — is now supplied by the separate
    :class:`MalolacticSenescence` baseline mortality (built as its own isolable Process so this SO₂
    kill remains byte-for-byte as D-39 built it). Total *O. oeni* mortality is therefore
    ``r_sen + r_death`` (benign baseline + SO₂-induced); a stabilizing SO₂ dose still dominates,
    crashing the population in ~1–3 d on top of the ~2-month benign half-life.

    **Conservation — a carbon/nitrogen-neutral transfer (the D-13 pattern).** Since D-38 both
    ``X_mlf`` and ``X_mlf_dead`` are weighted in ``total_carbon``/``total_nitrogen`` at the *same*
    biomass fractions, so moving a gram from one to the other (``d[X_mlf] = −r``,
    ``d[X_mlf_dead] = +r``) is carbon- and nitrogen-neutral **by construction** — exactly like the
    yeast :class:`~fermentation.core.kinetics.inactivation.EthanolInactivation` transfer
    ``X → X_dead``. No new conservation code, no sugar draw; touches only ``(X_mlf, X_mlf_dead)``.

    **Isolability.** ``X_mlf ≤ 0`` (undosed / un-pitched) *or* no total SO₂ returns a zero
    contribution *before* the pH ``brentq`` — the ``so2_total ≤ 0`` guard is exact (death is
    identically 0 without SO₂, mirroring the ``total_so2 > 0`` shortcut inside the toxicity gate),
    so a pitched-but-unsulfited run pays no per-RHS pH solve and its contribution is byte-for-byte
    zero. The compile seam enables this Process with the other pitch-gated MLF Processes
    (:data:`~fermentation.core.media._MLF_PROCESSES`) — it needs no amino acids (bacteria die
    whether or not they were growing), so it is pitch-gated, not amino-acid-gated like growth.
    Consequence: on any *pitched* run ``X_mlf``/``X_mlf_dead`` report **speculative** (this enabled
    Process touches them) — honest, since a pitched population that can be sulfited has a
    speculative trajectory. Tier **speculative** (``k_death_mlf``/``E_a_death_mlf`` are estimates,
    and the SO₂-driven form with arrest-scale reused as kill-scale is a modelling choice).
    """

    name = "malolactic_death"
    tier = Tier.SPECULATIVE
    #: Viable bacteria leave ``X_mlf``; the same mass enters the non-viable ``X_mlf_dead`` pool.
    #: Declaring both keeps the carbon/nitrogen-neutral transfer inside the ``touches`` contract.
    touches = ("X_mlf", "X_mlf_dead")
    #: ``k_death_mlf`` sets the death magnitude at full SO₂ kill; ``E_a_death_mlf``/``T_ref`` its
    #: Arrhenius temperature shape; ``molecular_so2_inhib_mlf`` is the SO₂ decay scale (reused from
    #: the conversion gate — arrest-scale = kill-scale, a documented simplification). NOT the
    #: pH/ethanol toxicity params (death drops them, D-39). Their tiers cap the ``X_mlf``/
    #: ``X_mlf_dead`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_death_mlf", "E_a_death_mlf", "T_ref", "molecular_so2_inhib_mlf")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE the pH brentq (the SO₂ partition reads the solved pH): no catalyst ⇒
        # no bacteria to kill; no total SO₂ ⇒ death is identically 0 (g_SO₂ = 1, so 1 − g_SO₂ = 0).
        # The so2_total ≤ 0 guard is EXACT, not an approximation, so an unsulfited pitched run pays
        # no per-RHS pH solve and its contribution is byte-for-byte zero (mirrors the tox gate).
        x_mlf = max(float(y[schema.slice("X_mlf")][0]), 0.0) if "X_mlf" in schema else 0.0
        if x_mlf <= 0.0:
            return d
        total_so2 = float(y[schema.slice(SO2_STATE_KEY)][0]) if SO2_STATE_KEY in schema else 0.0
        if total_so2 <= 0.0:
            return d

        ph = ph_of_state(y, schema, params)
        molecular_so2 = molecular_so2_at_ph(y, schema, params, ph)
        g_so2 = math.exp(-molecular_so2 / params["molecular_so2_inhib_mlf"])  # (0, 1]; 1 = no SO₂
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_death_mlf"], params["T_ref"])
        # SO₂ drives death: 1 − g_SO₂ ∈ [0, 1) rises with molecular SO₂. Arrhenius (not the cardinal
        # γ(T)) carries temperature, so cold preserves rather than kills (D-39).
        r_death = params["k_death_mlf"] * x_mlf * (1.0 - g_so2) * f_t  # [g X_mlf/L/h]
        d[schema.slice("X_mlf")] = -r_death
        d[schema.slice("X_mlf_dead")] = r_death  # carbon/nitrogen-neutral: same biomass fractions
        return d


class MalolacticSenescence(Process):
    """*Oenococcus oeni* benign senescence — slow baseline mortality (MLF **v2**, D-41/D-52/D-53).

    Lifts the owned v1 tradeoff of :class:`MalolacticDeath` (D-39): *"without SO₂, bacteria never
    die."* This Process is a small, always-on (when pitched) baseline mortality that moves viable
    ``X_mlf`` into the same non-viable ``X_mlf_dead`` pool the SO₂ kill uses, so a pitched wine left
    alone can in principle lose its bacteria (and the ``X_mlf``-scaled activities — conversion,
    citrate → diacetyl, lees-contact diacetyl reduction — fade with them) instead of holding a
    viable culture indefinitely. **D-53 correction (2026-07-07):** a real-wine literature check
    found no support for the D-41 "weeks-to-months spontaneous decline" premise — Windholtz et al.
    2025 (OENO One) and the Millet 2001 thesis both show O. oeni populations *stable* for 3–5
    months in real SO₂-free wine, and the steep decline the original citations pointed at turns out
    to be SO₂-driven (:class:`MalolacticDeath`'s territory), not spontaneous. The magnitude below is
    corrected accordingly — see the ``k_senescence_mlf`` provenance for full sourcing.

    **Rate — a baseline rate scaled by a bounded ethanol/starvation stress multiplier, Arrhenius
    temperature (D-52).**

        r_sen = k_senescence_mlf · X_mlf · arrhenius(T, E_a_death_mlf, T_ref) · stress
        stress = 1 + k_senescence_ethanol_scale·[E/(E+ethanol_tolerance_mlf)]
                   + k_senescence_starvation_scale·[K_aa_mlf/(K_aa_mlf+amino_acids)]

    * **Constant baseline ``k_senescence_mlf``** (t½ ≈ 7.9 years at ``T_ref``, ``stress`` = 1;
      D-53-corrected from D-41's original ~8-week value — see below) — deliberately tiny, well
      below the full-SO₂-kill ``k_death_mlf``.
    * **The D-52 stress multiplier — lifts the D-41 "environment-free" deferral, without
      reintroducing the D-39 wipeout bug.** D-39's crux was a *large* rate (``k_death_mlf``,
      calibrated to represent a full SO₂ kill) multiplied by ``1 − toxicity`` ≈ 0.92 at ordinary
      post-AF ethanol — near-maximal death *from ethanol alone*. Here the multiplier scales the
      **tiny** senescence baseline instead, and both stress terms are smooth Monod-type factors
      bounded in **[0, 1)** by construction (no clamp needed, C¹ for the BDF solver) rather than the
      Luong wall's near-binary "1 at zero stress, 0 at the tolerance wall" shape. ``stress`` is
      therefore hard-capped at ``1 + k_senescence_ethanol_scale + k_senescence_starvation_scale``
      regardless of how far ``E`` or nutrient depletion runs — at the shipped values (1.0/0.5) that
      ceiling is 2.5×. **Post-D-53** even that worst case gives a multi-year half-life (~3.16 y):
      the stress mechanism is unchanged from D-52, but it now scales a baseline small enough that
      senescence is honestly negligible on every timescale this model simulates — the D-39
      wipeout regime (~1 week) is nowhere close on any axis. Reuses ``ethanol_tolerance_mlf`` and
      ``K_aa_mlf`` as the two terms'
      half-saturation points — the same "arrest-scale reused as a death-adjacent scale"
      simplification :class:`MalolacticDeath` already makes with ``molecular_so2_inhib_mlf`` — so no
      new concentration-scale parameters are introduced, only the two dimensionless ceilings.
    * **Starvation term reuses the growth fuel pool** (``amino_acids``, the same pool
      :class:`MalolacticGrowth` draws on): it is ≈1 (near-max stress) once amino acids are
      exhausted — which the D-23 finding places at ~1.3 d post-pitch regardless of dose — and falls
      back toward 0 whenever autolysis (D-34) refills the pool, so the term is not a flat add-on but
      tracks the real nutrient-refill dynamic already in the model.
    * **Arrhenius ``arrhenius(T, E_a_death_mlf, T_ref)`` — NOT the cardinal γ(T)** (the reused D-39
      choice). Warm accelerates senescence, cold slows it to dormancy — the physically correct
      direction. The cardinal γ(T) peaks at ``T_opt_mlf`` (23 °C) and vanishes past ``T_max_mlf``,
      which would make senescence *maximal at the growth optimum* and *switch off* in the warm —
      exactly backwards for a decline. Reuses ``E_a_death_mlf``/``T_ref`` (no new temperature
      params); the factor is 1 at the 20 °C benchmark, like every other Arrhenius rate.
    * The SO₂-driven acute kill remains :class:`MalolacticDeath`'s job — total *O. oeni* mortality
      is therefore ``r_sen + r_death`` (stress-modulated baseline + SO₂-induced), the two built as
      **separate isolable Processes** (prime directive #3) so the SO₂ lever stays byte-for-byte as
      D-39 built it and this baseline toggles off independently.

    **Conservation — the carbon/nitrogen-neutral transfer, no new code (the D-13/D-39 pattern).**
    Both ``X_mlf`` and ``X_mlf_dead`` are weighted in ``total_carbon``/``total_nitrogen`` at the
    *same* biomass fractions (since D-38/D-39), so ``d[X_mlf] = −r_sen``, ``d[X_mlf_dead] = +r_sen``
    is C- and N-neutral by construction — identical to the SO₂ kill and the yeast ``X → X_dead``
    inactivation. ``X_mlf_dead`` is a **terminal sink** here: :class:`~fermentation.core.kinetics.\
autolysis.YeastAutolysis` reads only the yeast ``X_dead`` pool, so senescing bacteria do **not**
    refuel the ``amino_acids`` pool (no self-cancelling recycling loop). Touches ``(X_mlf,
    X_mlf_dead)`` only — reading ``E``/``amino_acids`` for the D-52 stress terms adds no new touched
    state.

    **Isolability + performance.** ``X_mlf ≤ 0`` (undosed / un-pitched) returns a zero contribution.
    This Process still reads **no SO₂ and no pH**, so it never triggers a ``brentq`` — it remains
    strictly cheaper than the SO₂ kill even with the D-52 stress terms (``E``/``amino_acids`` are
    read directly off state, no equilibrium solve). Pitch-gated at the compile seam (enabled with
    the other ``_MLF_PROCESSES`` when ``mlf_pitch_gpl > 0``), NOT amino-acid-gated: bacteria age
    whether or not they were growing. **Supersedes the v1 "no-SO₂ pitched run is byte-for-byte
    inert" property in structure, not in observable magnitude**: a pitched, unsulfited run shows a
    slow monotone ``X_mlf`` decline in principle, but post-D-53 that decline is honestly negligible
    at the timescales this model simulates — matching the real-wine finding of no detectable
    spontaneous die-off within 3–5 months. Tier **speculative** (``k_senescence_mlf`` and the D-52
    stress-ceiling parameters are author estimates; direction — ethanol/starvation stress
    accelerates decline — is sourced, magnitude is upper-bound-derived, not fitted).
    """

    name = "malolactic_senescence"
    tier = Tier.SPECULATIVE
    #: Viable bacteria leave ``X_mlf`` for the same non-viable ``X_mlf_dead`` pool the SO₂ kill
    #: fills. Declaring both keeps the carbon/nitrogen-neutral transfer in the ``touches`` contract.
    touches = ("X_mlf", "X_mlf_dead")
    #: ``k_senescence_mlf`` sets the baseline mortality magnitude; ``E_a_death_mlf``/``T_ref`` its
    #: Arrhenius temperature shape (shared with :class:`MalolacticDeath`). The D-52 stress terms
    #: read ``ethanol_tolerance_mlf``/``K_aa_mlf`` (half-saturation scales reused from the
    #: conversion/growth gates) and their own dimensionless ceilings
    #: ``k_senescence_ethanol_scale``/``k_senescence_starvation_scale``. Still NO SO₂/pH params.
    #: Their tiers cap the ``X_mlf``/``X_mlf_dead`` output tiers via parameter-tier propagation
    #: (D-1).
    reads: tuple[str, ...] = (
        "k_senescence_mlf",
        "E_a_death_mlf",
        "T_ref",
        "k_senescence_ethanol_scale",
        "ethanol_tolerance_mlf",
        "k_senescence_starvation_scale",
        "K_aa_mlf",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Single guard: no catalyst ⇒ no bacteria to age. No SO₂/pH read ⇒ NO brentq (unlike the SO₂
        # kill), so an undosed run is byte-for-byte zero at zero solve cost.
        x_mlf = max(float(y[schema.slice("X_mlf")][0]), 0.0) if "X_mlf" in schema else 0.0
        if x_mlf <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_death_mlf"], params["T_ref"])

        # D-52 bounded stress multiplier: two smooth Monod-type terms in [0, 1), each capped by its
        # own dimensionless ceiling, so `stress` cannot exceed 1 + ethanol_scale + starvation_scale
        # regardless of how far E or nutrient depletion runs (no clamp needed — the wipeout guard).
        e = max(float(y[schema.slice("E")][0]), 0.0)
        ethanol_stress = e / (e + params["ethanol_tolerance_mlf"])
        aa = max(float(y[schema.slice("amino_acids")][0]), 0.0) if "amino_acids" in schema else 0.0
        starvation_stress = params["K_aa_mlf"] / (params["K_aa_mlf"] + aa)
        stress = (
            1.0
            + params["k_senescence_ethanol_scale"] * ethanol_stress
            + params["k_senescence_starvation_scale"] * starvation_stress
        )

        # Stress-modulated baseline mortality: warm-accelerated by Arrhenius (NOT γ(T), which would
        # spuriously peak at the growth optimum), scaled up (bounded) by ethanol/starvation stress.
        r_sen = params["k_senescence_mlf"] * x_mlf * f_t * stress  # [g X_mlf/L/h]
        d[schema.slice("X_mlf")] = -r_sen
        d[schema.slice("X_mlf_dead")] = r_sen  # carbon/nitrogen-neutral: same biomass fractions
        return d
