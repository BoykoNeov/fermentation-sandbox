"""*Brettanomyces / Dekkera* — the volatile-phenol spoilage pathway (decision D-40).

This module holds the *Brettanomyces bruxellensis* Processes, the mixed-culture beat that
closes Milestone 2. Brett is the canonical wine spoilage yeast: it decarboxylates grape-must
**hydroxycinnamic acids** (p-coumaric, ferulic) to **vinylphenols**, then reduces those to the
**ethylphenols** (4-ethylphenol "horse-sweat/barnyard", 4-ethylguaiacol "clove/smoky") that
define Brett character. Both precursor pair and product pair are **lumped** (fork-2, D-40):
``hydroxycinnamics`` (booked as p-coumaric), ``vinylphenols``, ``ethylphenols``.

**Brett carries BOTH enzymes — that is why it spoils normal wine on its own.** Cinnamate
decarboxylase (:class:`BrettDecarboxylation`) *and* vinylphenol reductase
(:class:`BrettVinylphenolReduction`) both live in Brett, so a dosed Brett culture takes must
hydroxycinnamics all the way to 4-ethylphenol with no help. *S. cerevisiae* has at most the
decarboxylase (and only if POF+, the uncommon case), never the reductase — so a POF+ yeast can
fill the shared ``vinylphenols`` reservoir during AF but cannot clear it; if Brett is absent the
vinylphenol *strands* (the emergent yeast/Brett coupling, the α-acetolactate-reservoir parallel of
D-26/D-31). The POF+ yeast decarboxylase is a separate opt-in strain Process (decision D-40 pt4);
this module is Brett, which spoils POF-negative wine unaided.

**Carbon closes on the existing ledger — no new conservation code.** With ``r`` [mol/L/h] the
decarboxylase turnover and ``L`` [g/L/h] the reductase mass flux,

    d(hydroxycinnamics)/dt = −r·M_p_coumaric   d(vinylphenols)/dt = +r·M_vinylphenol − L
    d(CO2)/dt              = +r·M_CO2           d(ethylphenols)/dt = +L·M_ethylphenol/M_vinylphenol

p-coumaric (9 C) → vinylphenol (8 C) + CO2 (1 C) closes mole-for-mole (9 = 8 + 1, the malic →
lactic + CO2 idiom, D-23); vinylphenol (8 C) → ethylphenol (8 C) is a mole-for-mole C8 → C8
transfer between two weighted pools (the diacetyl → butanediol idiom, D-26). :func:`~fermentation.
validation.conservation.total_carbon` weights all three phenol pools, so the Processes touch only
``hydroxycinnamics``/``vinylphenols``/``ethylphenols``/``CO2`` and add nothing to the harness.

**Environmental gate — SO₂ and temperature only (decision D-40).** Unlike *O. oeni*, Brett is
markedly **acid-tolerant** (it spoils low-pH wine) and **ethanol-tolerant** (a full-strength-wine
aging spoiler), so — heeding the explicit design warning not to copy the MLF gate — the Brett gate
carries **no pH logistic and no Luong ethanol wall**: those would spuriously arrest Brett exactly
where it thrives. What controls Brett in the cellar is **molecular SO₂** (the antimicrobial
readout, D-22) and **temperature** (a cardinal optimum, warmer than *O. oeni*'s). So

    gate = g_SO₂ · γ(T),   g_SO₂ = exp(−[SO₂]_molecular / molecular_so2_inhib_brett)

with γ(T) the Rosso cardinal factor (:func:`~fermentation.core.kinetics.malolactic.cardinal_\
temperature_factor`, reused). SO₂ suppresses phenol production metabolically here; it also *kills*
Brett biomass (:class:`BrettDeath`, D-40 pt3) — the two mechanisms the winemaker's SO₂ addition
combines. Molecular SO₂ needs the solved pH, so a dosed-SO₂ run pays one ``brentq`` per RHS; an
unsulfited run skips it (the ``so2_total ≤ 0`` shortcut) and gate = γ(T).

**Isolability.** ``X_brett`` is a dosed catalyst in v1 pt1 (constant; :class:`BrettGrowth`, pt2,
makes it dynamic). Each Process returns a zero contribution before any SO₂/pH work when there is no
catalyst (``X_brett ≤ 0``) or no substrate, so an unpitched wine run is byte-for-byte the validated
core; the compile seam additionally *disables* the Processes when Brett is not pitched, so the
inert ``hydroxycinnamics``/``vinylphenols``/``ethylphenols`` slots keep their VALIDATED tier
(``tier_of`` counts enabled, not nonzero, Processes — the D-23 MLF pattern). Wine-only.

Tier: **speculative**. The decarboxylase → reductase topology and the gate *directions* are sound,
but every rate/gate/cardinal magnitude is an author estimate (no per-catalyst kinetic model of this
flux form is sourced), and v1 lumps the two phenol branches. Parameter-tier propagation (D-1) caps
the phenol-pool outputs at speculative regardless.
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
    M_CO2,
    M_ETHYLPHENOL,
    M_P_COUMARIC,
    M_VINYLPHENOL,
)
from fermentation.core.kinetics.malolactic import cardinal_temperature_factor
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The *Brettanomyces* environmental-gate parameters, shared by the decarboxylase and reductase
#: Processes (and, from D-40 pt3, the SO₂ scale by :class:`BrettDeath`). Declared once so the
#: Processes' ``reads`` tuples and :func:`brett_environmental_gate` cannot drift apart. Only SO₂ and
#: the temperature cardinals — no pH/ethanol terms (Brett is acid- and ethanol-tolerant, D-40).
_BRETT_GATE_READS: tuple[str, ...] = (
    "molecular_so2_inhib_brett",
    "T_min_brett",
    "T_opt_brett",
    "T_max_brett",
)


def brett_environmental_gate(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], ph: float
) -> float:
    """The shared *Brettanomyces* environmental gate ``g_SO₂ · γ(T)`` ∈ [0, 1] (decision D-40).

    Every Brett metabolic activity — hydroxycinnamate decarboxylation, vinylphenol reduction and
    (pt3) growth — is throttled by the *same* environment, so they multiply their rate by this one
    factor: the molecular-SO₂ exponential (the antimicrobial readout, D-22, partitioned at the
    *already-solved* ``ph`` the caller passes so there is no second ``brentq``) times the Rosso
    cardinal-temperature optimum γ(T). **No pH or ethanol term** — Brett is acid- and
    ethanol-tolerant across the wine range, so those MLF gate terms are deliberately absent
    (copying them would arrest Brett where it in fact thrives; decision D-40). ≈1 in a warm,
    unsulfited cellar and →0 as SO₂ accumulates or temperature leaves the cardinal window.
    """
    total_so2 = float(y[schema.slice(SO2_STATE_KEY)][0]) if SO2_STATE_KEY in schema else 0.0
    if total_so2 > 0.0:
        molecular_so2 = molecular_so2_at_ph(y, schema, params, ph)
        gate_so2 = math.exp(-molecular_so2 / params["molecular_so2_inhib_brett"])
    else:
        gate_so2 = 1.0

    temp = float(y[schema.slice("T")][0])
    gamma_t = cardinal_temperature_factor(
        temp, params["T_min_brett"], params["T_opt_brett"], params["T_max_brett"]
    )
    return float(gate_so2 * gamma_t)


def _needs_ph_solve(y: FloatArray, schema: StateSchema) -> bool:
    """Whether the SO₂ partition (and thus a ``brentq``) is needed — true iff SO₂ is dosed.

    The gate's only pH-dependent term is molecular SO₂; with no total SO₂ the gate is just γ(T),
    so an unsulfited run must not pay a per-RHS pH solve (mirrors :class:`~fermentation.core.\
    kinetics.malolactic.MalolacticDeath`'s exact ``so2_total ≤ 0`` shortcut).
    """
    total_so2 = float(y[schema.slice(SO2_STATE_KEY)][0]) if SO2_STATE_KEY in schema else 0.0
    return total_so2 > 0.0


class BrettDecarboxylation(Process):
    """*Brettanomyces* cinnamate decarboxylase — hydroxycinnamics → vinylphenols + CO2 (D-40).

    ``d(hydroxycinnamics)/dt = −r·M_p_coumaric``, ``d(vinylphenols)/dt = +r·M_vinylphenol``,
    ``d(CO2)/dt = +r·M_CO2`` with the molar turnover ``r = k_brett_decarb · X_brett · [hc]/(K_
    hydroxycinnamic + [hc]) · g_SO₂ · γ(T)`` (Michaelis–Menten in the precursor, catalyst-scaled,
    gated). p-coumaric (9 C) → vinylphenol (8 C) + CO2 (1 C) closes carbon mole-for-mole on the
    existing ledger. The vinylphenol it makes feeds the shared reservoir
    :class:`BrettVinylphenolReduction` (and a POF+ yeast, D-40 pt4) drains. Touches only
    ``hydroxycinnamics``/``vinylphenols``/``CO2``; reads the (constant, in pt1) catalyst
    ``X_brett`` plus the gate state (T, SO₂, solved pH).

    Returns a zero contribution before any pH work when undosed (``X_brett ≤ 0``) or when the
    precursor is exhausted — structural value-isolability and no wasted ``brentq`` (the compile
    seam additionally disables the Process when Brett is not pitched, for tier isolability).
    Tier **speculative** (rate/gate magnitudes are estimates).
    """

    name = "brett_decarboxylation"
    tier = Tier.SPECULATIVE
    touches = ("hydroxycinnamics", "vinylphenols", "CO2")
    #: The decarboxylase Monod pair plus the shared Brett environmental-gate parameters. Their
    #: tiers cap the hydroxycinnamics/vinylphenols/CO2 output tier via parameter-tier propagation
    #: (D-1). CO2 is already speculative (the always-on VDK decarboxylation), so this adds no new
    #: tier headline. pKa/SO₂ params read via ``acidbase`` are omitted (all plausible; the Process
    #: is already speculative — the MalolacticConversion convention).
    reads: tuple[str, ...] = ("k_brett_decarb", "K_hydroxycinnamic", *_BRETT_GATE_READS)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE any pH solve: no catalyst or no precursor ⇒ no conversion, and an
        # unsulfited/undosed run must not pay a per-RHS pH solve for a zero (or γ(T)-only) result.
        x_brett = max(float(y[schema.slice("X_brett")][0]), 0.0) if "X_brett" in schema else 0.0
        if x_brett <= 0.0:
            return d
        hc_gpl = max(float(y[schema.slice("hydroxycinnamics")][0]), 0.0)
        if hc_gpl <= 0.0:
            return d

        ph = ph_of_state(y, schema, params) if _needs_ph_solve(y, schema) else 0.0
        gate = brett_environmental_gate(y, schema, params, ph)

        hc_molar = hc_gpl / M_P_COUMARIC
        monod = hc_molar / (params["K_hydroxycinnamic"] + hc_molar)
        r = params["k_brett_decarb"] * x_brett * monod * gate  # decarboxylase turnover, mol/L/h

        d[schema.slice("hydroxycinnamics")] = -r * M_P_COUMARIC
        d[schema.slice("vinylphenols")] = r * M_VINYLPHENOL  # feeds the shared reductase reservoir
        d[schema.slice("CO2")] = (
            r * M_CO2
        )  # p-coumaric C9 → vinylphenol C8 + CO2 C1 (carbon-closing)
        return d


class BrettVinylphenolReduction(Process):
    """*Brettanomyces* vinylphenol reductase — vinylphenols → ethylphenols (decision D-40).

    ``d(vinylphenols)/dt = −L``, ``d(ethylphenols)/dt = +L·M_ethylphenol/M_vinylphenol`` with the
    mass flux ``L = k_brett_reduction · X_brett · g_SO₂ · γ(T) · [vinylphenols]`` (first-order in
    the intermediate, catalyst-scaled, gated by the *same* Brett environment as the decarboxylase).
    A mole-for-mole C8 → C8 transfer (both pools weighted at their own carbon fraction), so it is
    carbon-neutral like the diacetyl → butanediol reduction (D-26); no sugar draw. This is the
    step *S. cerevisiae* lacks (even POF+), so it is what makes the shared ``vinylphenols``
    reservoir clear only when Brett is present — the emergent coupling (D-40 pt4).

    Guarded like the decarboxylase: a zero contribution before any pH work when there is no
    catalyst (``X_brett ≤ 0``) or no vinylphenol to reduce (value + perf isolability; compile-seam
    disable gives tier isolability). Tier **speculative** (rate magnitude estimate).
    """

    name = "brett_vinylphenol_reduction"
    tier = Tier.SPECULATIVE
    touches = ("vinylphenols", "ethylphenols")
    #: ``k_brett_reduction`` sets the reductase magnitude; the shared Brett environmental-gate
    #: parameters throttle it by SO₂/temperature. Their tiers cap the vinylphenols/ethylphenols
    #: output tier via parameter-tier propagation (D-1). Reads the constant catalyst ``X_brett``.
    reads: tuple[str, ...] = ("k_brett_reduction", *_BRETT_GATE_READS)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        vinylphenols = max(float(y[schema.slice("vinylphenols")][0]), 0.0)
        if vinylphenols <= 0.0:  # nothing to reduce
            return d
        x_brett = max(float(y[schema.slice("X_brett")][0]), 0.0) if "X_brett" in schema else 0.0
        if x_brett <= 0.0:  # no Brett ⇒ no reduction (the reservoir strands)
            return d

        ph = ph_of_state(y, schema, params) if _needs_ph_solve(y, schema) else 0.0
        gate = brett_environmental_gate(y, schema, params, ph)

        loss = params["k_brett_reduction"] * x_brett * gate * vinylphenols  # [g vinylphenol/L/h]
        d[schema.slice("vinylphenols")] = -loss
        d[schema.slice("ethylphenols")] = (
            loss * M_ETHYLPHENOL / M_VINYLPHENOL
        )  # mole-for-mole C8→C8
        return d
