"""Malolactic fermentation (MLF) v1 — conversion-only (decision D-23).

*Oenococcus oeni* converts L-malic acid (C4, diprotic) to L-lactic acid (C3,
monoprotic) plus CO2, mole-for-mole. That single reaction deacidifies the wine
(pH rises ~0.1–0.3, the D-18 headline coupling) and softens its perceived acidity.
This module is the first **RHS consumer** of the pH charge-balance keystone (D-18)
and of the molecular-SO₂ readout (D-22): the conversion rate is gated by the *solved*
pH, by molecular (antimicrobial) SO₂, by ethanol, and by a temperature optimum — so
the deacidification feedback (pH ↑ ⇒ rate ↑, self-limited as malate depletes) and the
SO₂/ethanol arrest of MLF *emerge* from the model rather than being scripted.

**v1 is conversion-only (decision D-23).** *O. oeni* builds biomass mostly from amino
acids/peptides, but the lumped ``N`` (YAN) is carbon-free in :func:`total_carbon`
(D-19) and is driven to ~0 within ~1.3 d of the AF pitch *regardless of dose* (the
empirical finding that settles D-23), so there is no nitrogen at the MLF pitch point to
fund bacterial growth. Modelling MLF-growth honestly therefore needs a separate
amino-acid ledger *and* an autolytic-peptide refill source, both deferred. So in v1 the
bacterium is a **dosed-but-inert catalyst**: ``X_mlf`` is a constant concentration that
*scales* the conversion rate, and **no Process grows or kills it**. The later growth beat
is then a clean extension (add a Process touching ``X_mlf``), not a refactor.

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
from fermentation.core.chemistry import M_CO2, M_LACTIC, M_MALIC
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
    reads: tuple[str, ...] = (
        "k_mlf",
        "K_mlf",
        "pH_half_mlf",
        "ethanol_tolerance_mlf",
        "mlf_ethanol_exponent",
        "molecular_so2_inhib_mlf",
        "T_min_mlf",
        "T_opt_mlf",
        "T_max_mlf",
    )

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
        # molecular_so2_at_ph(…, ph) avoids a second brentq solve inside acidbase.speciate_so2.
        ph = ph_of_state(y, schema, params)
        gate_ph = 1.0 / (1.0 + 10.0 ** (params["pH_half_mlf"] - ph))

        # Antimicrobial suppression is by MOLECULAR SO₂, the undissociated share of FREE SO₂.
        # Under D-28 the dosed slot is *total* SO₂; free = total − acetaldehyde-bound, so as
        # acetaldehyde peaks it sequesters SO₂ and the suppression correctly weakens (bound SO₂
        # is not antimicrobial). At acetaldehyde = 0 this equals the D-22 free × fraction(pH).
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

        malate_molar = malic_gpl / M_MALIC
        monod = malate_molar / (params["K_mlf"] + malate_molar)
        r = params["k_mlf"] * x_mlf * monod * gate_ph * gate_eth * gate_so2 * gamma_t

        d[schema.slice("malic")] = -r * M_MALIC
        d[schema.slice("lactic")] = r * M_LACTIC
        d[schema.slice("CO2")] = r * M_CO2
        return d
