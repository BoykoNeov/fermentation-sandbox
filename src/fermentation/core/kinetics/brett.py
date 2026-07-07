"""*Brettanomyces / Dekkera* — the volatile-phenol spoilage pathway (decision D-40, split D-55).

This module holds the *Brettanomyces bruxellensis* Processes, the mixed-culture beat that
closes Milestone 2. Brett is the canonical wine spoilage yeast: it decarboxylates grape-must
**hydroxycinnamic acids** (p-coumaric, ferulic) to **vinylphenols** (4-vinylphenol,
4-vinylguaiacol), then reduces those to the **ethylphenols** (4-ethylphenol "horse-sweat/
barnyard", 4-ethylguaiacol "clove/smoky") that define Brett character. D-40 originally **lumped**
both the precursor pair and the product pair into one p-coumaric-booked chain
(``hydroxycinnamics``/``vinylphenols``/``ethylphenols``); decision D-55 split it into two
genuinely distinct parallel chains — p-coumaric (``hydroxycinnamics``/``vinylphenols``/
``ethylphenols``) and ferulic (``ferulic_acid``/``vinylguaiacols``/``ethylguaiacols``) — because
ferulic acid is a different-carbon-count molecule (10 C vs p-coumaric's 9 C) that cannot be
represented as a fixed-ratio split of the p-coumaric flow without breaking carbon closure. Both
chains share the same enzymes, catalysts and gates; only the substrate-specific rate constants
differ (ratio-derived from paired literature kinetics, not independently re-estimated).

**Brett carries BOTH enzymes — that is why it spoils normal wine on its own.** Cinnamate
decarboxylase (:class:`BrettDecarboxylation`) *and* vinylphenol reductase
(:class:`BrettVinylphenolReduction`) both live in Brett, so a dosed Brett culture takes must
hydroxycinnamics all the way to 4-ethylphenol with no help. *S. cerevisiae* has at most the
decarboxylase (and only if POF+, the uncommon case), never the reductase — so a POF+ yeast can
fill the shared ``vinylphenols`` reservoir during AF but cannot clear it; if Brett is absent the
vinylphenol *strands* (the emergent yeast/Brett coupling, the α-acetolactate-reservoir parallel of
D-26/D-31). That POF+ yeast decarboxylase is :class:`YeastPOFDecarboxylation` (decision D-40 pt4) —
a **separate opt-in strain** Process living in this module: it shares the phenol species and the
carbon routing with :class:`BrettDecarboxylation`, but its catalyst is *viable yeast* (coupled to
the fermentative flux) rather than ``X_brett``, and it is enabled by a POF+ strain opt-in, wholly
independent of the Brett pitch. The rest of this module is Brett, which spoils POF-negative wine
unaided.

**Carbon closes on the existing ledger — no new conservation code.** With ``r`` [mol/L/h] the
decarboxylase turnover and ``L`` [g/L/h] the reductase mass flux,

    d(hydroxycinnamics)/dt = −r·M_p_coumaric   d(vinylphenols)/dt = +r·M_vinylphenol − L
    d(CO2)/dt              = +r·M_CO2           d(ethylphenols)/dt = +L·M_ethylphenol/M_vinylphenol

p-coumaric (9 C) → vinylphenol (8 C) + CO2 (1 C) closes mole-for-mole (9 = 8 + 1, the malic →
lactic + CO2 idiom, D-23); vinylphenol (8 C) → ethylphenol (8 C) is a mole-for-mole C8 → C8
transfer between two weighted pools (the diacetyl → butanediol idiom, D-26). The ferulic branch
(decision D-55) is the identical shape on its own substrate: ferulic (10 C) → vinylguaiacol (9 C)
+ CO2 (1 C) closes 10 = 9 + 1, and vinylguaiacol (9 C) → ethylguaiacol (9 C) is a mole-for-mole
C9 → C9 transfer. :func:`~fermentation.validation.conservation.total_carbon` weights all six
phenol pools, so the Processes touch only ``hydroxycinnamics``/``vinylphenols``/``ferulic_acid``/
``vinylguaiacols``/``ethylphenols``/``ethylguaiacols``/``CO2`` and add nothing to the harness.

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
catalyst (``X_brett ≤ 0``) or no substrate in EITHER branch, so an unpitched wine run is
byte-for-byte the validated core; the compile seam additionally *disables* the Processes when
Brett is not pitched, so the inert phenol slots (both branches) keep their VALIDATED tier
(``tier_of`` counts enabled, not nonzero, Processes — the D-23 MLF pattern). Wine-only.

Tier: **speculative**. The decarboxylase → reductase topology and the gate *directions* are sound,
but every rate/gate/cardinal magnitude is an author estimate (no per-catalyst kinetic model of this
flux form is sourced). The D-40 lumping is now (D-55) split into two genuinely distinct parallel
chains with paired-literature-sourced *relative* rates, though each chain's absolute magnitude
remains an author estimate. Parameter-tier propagation (D-1) caps the phenol-pool outputs at
speculative regardless.
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
    M_ETHYLGUAIACOL,
    M_ETHYLPHENOL,
    M_FERULIC,
    M_P_COUMARIC,
    M_VINYLGUAIACOL,
    M_VINYLPHENOL,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import fermentative_flux_shape
from fermentation.core.kinetics.malolactic import cardinal_temperature_factor
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The *Brettanomyces* environmental-gate parameters, shared by the decarboxylase, reductase and
#: growth Processes. Declared once so the Processes' ``reads`` tuples and
#: :func:`brett_environmental_gate` cannot drift apart. Only SO₂ and the temperature cardinals — no
#: pH/ethanol terms (Brett is acid- and ethanol-tolerant, D-40). :class:`BrettDeath` (pt3) does NOT
#: splat this tuple: it uses an *Arrhenius* temperature factor, not the cardinal γ(T), so it reads
#: only ``molecular_so2_inhib_brett`` (the SO₂ scale) explicitly, not the ``T_*_brett`` cardinals.
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


#: The ethanol-toxicity survival-factor parameters, shared by :class:`BrettGrowth` (as an upper
#: growth wall) and :class:`BrettEthanolToxicity` (as the death driver, decision D-58). Declared
#: once so the two Processes' ``reads`` tuples and :func:`brett_ethanol_survival_factor` cannot
#: drift apart — the ``_BRETT_GATE_READS`` pattern.
_BRETT_ETHANOL_TOXICITY_READS: tuple[str, ...] = (
    "brett_ethanol_toxicity_onset",
    "brett_ethanol_toxicity_ceiling",
    "brett_ethanol_toxicity_exponent",
)


def brett_ethanol_survival_factor(e: float, params: Mapping[str, float]) -> float:
    """The *Brettanomyces* ethanol-toxicity survival factor ∈ [0, 1] (decision D-58).

    Sourced from Barata et al. 2008 (Int. J. Food Microbiol. 121(2):201–207): Brett grows normally
    up to ``brett_ethanol_toxicity_onset`` (~14% v/v, 110 g/L) and is fully arrested by
    ``brett_ethanol_toxicity_ceiling`` (~15% v/v, 118 g/L) — a THRESHOLD effect, deliberately
    unlike the O. oeni Luong wall (:func:`~fermentation.core.kinetics.malolactic.\
    malolactic_toxicity_gate`), which decays continuously from ``E = 0``. A whole-range Luong wall
    would be the wrong functional form here: Brett is documented (D-40) as markedly MORE
    ethanol-tolerant than O. oeni, so a wall with a similar shape would spuriously suppress Brett
    across the entire normal wine ethanol range (~90–105 g/L) where it in fact thrives — exactly
    the mistake the Brett environmental gate already avoids by carrying no ethanol term at all
    (D-40). This factor is 1 (no effect) for ``E ≤ onset`` — so it is a genuine zero contribution,
    not just a small one, across ordinary wine strength — and eases smoothly (C1, no BDF kink, the
    D-40 pt2 shadow idiom) to 0 by the ceiling.

    ``BrettGrowth`` multiplies this into its rate as an upper wall (reconciling ethanol's dual
    role — carbon SOURCE at low concentration via the existing ``E/(K_E_brett+E)`` Monod, toxin at
    high concentration via this factor — on the SAME state variable). ``BrettEthanolToxicity``
    drives death off ``1 − this factor`` (the ``BrettDeath`` ``1 − g_SO₂`` idiom).
    """
    onset = params["brett_ethanol_toxicity_onset"]
    span = params["brett_ethanol_toxicity_ceiling"] - onset
    remaining = 1.0 - max(0.0, e - onset) / span
    return remaining ** params["brett_ethanol_toxicity_exponent"] if remaining > 0.0 else 0.0


def _needs_ph_solve(y: FloatArray, schema: StateSchema) -> bool:
    """Whether the SO₂ partition (and thus a ``brentq``) is needed — true iff SO₂ is dosed.

    The gate's only pH-dependent term is molecular SO₂; with no total SO₂ the gate is just γ(T),
    so an unsulfited run must not pay a per-RHS pH solve (mirrors :class:`~fermentation.core.\
    kinetics.malolactic.MalolacticDeath`'s exact ``so2_total ≤ 0`` shortcut).
    """
    total_so2 = float(y[schema.slice(SO2_STATE_KEY)][0]) if SO2_STATE_KEY in schema else 0.0
    return total_so2 > 0.0


def _decarboxylation_branch(
    precursor_gpl: float,
    precursor_molar_mass: float,
    product_molar_mass: float,
    k: float,
    k_half_saturation: float,
    activity: float,
) -> tuple[float, float, float]:
    """One hydroxycinnamate-decarboxylase branch: precursor → product + CO2 (decision D-55).

    Shared by both :class:`BrettDecarboxylation` and :class:`YeastPOFDecarboxylation`, and by both
    the p-coumaric branch (``hydroxycinnamics`` → ``vinylphenols``) and the ferulic branch
    (``ferulic_acid`` → ``vinylguaiacols``) within each — the two branches are the *same* reaction
    on a different substrate, differing only in molar masses and the rate/half-saturation
    parameters (:class:`BrettDecarboxylation`'s docstring covers the D-55 fork for why this is a
    genuine second precursor pool, not a fixed-ratio split of the existing one).

    ``activity`` folds in everything upstream of the precursor Monod that both branches share
    equally (the catalyst/flux magnitude, and any gate/Arrhenius factor) — ``X_brett · gate`` for
    Brett, ``flux · arrhenius(T, E_a_pof)`` for POF — so this helper only computes the Monod
    turnover and the resulting carbon-closing mass flux: ``r = k · activity · [precursor]/
    (k_half_saturation + [precursor])`` (mol/L/h), returned as
    ``(d_precursor, d_product, d_CO2)`` in g/L/h.

    Returns exactly ``(0.0, 0.0, 0.0)`` if the precursor is exhausted or ``activity ≤ 0`` — the
    same value+perf isolability guard each caller already applies before this helper runs.
    """
    if precursor_gpl <= 0.0 or activity <= 0.0:
        return 0.0, 0.0, 0.0
    precursor_molar = precursor_gpl / precursor_molar_mass
    monod = precursor_molar / (k_half_saturation + precursor_molar)
    r = k * activity * monod  # mol/L/h
    return (
        -r * precursor_molar_mass,
        r * product_molar_mass,
        r * M_CO2,
    )


class BrettDecarboxylation(Process):
    """*Brettanomyces* cinnamate decarboxylase — two branches to vinylphenols/CO2 (D-40, D-55).

    ``d(hydroxycinnamics)/dt = −r·M_p_coumaric``, ``d(vinylphenols)/dt = +r·M_vinylphenol``,
    ``d(CO2)/dt = +r·M_CO2`` with the molar turnover ``r = k_brett_decarb · X_brett · [hc]/(K_
    hydroxycinnamic + [hc]) · g_SO₂ · γ(T)`` (Michaelis–Menten in the precursor, catalyst-scaled,
    gated). p-coumaric (9 C) → vinylphenol (8 C) + CO2 (1 C) closes carbon mole-for-mole on the
    existing ledger. The vinylphenol it makes feeds the shared reservoir
    :class:`BrettVinylphenolReduction` (and a POF+ yeast, D-40 pt4) drains.

    **Ferulic-acid branch (decision D-55) — a genuine second precursor, not a split of the
    first.** ``hydroxycinnamics`` is booked as p-coumaric acid specifically (9 carbons); ferulic
    acid is a distinct 10-carbon molecule whose decarboxylation (10 = 9 + 1, to 4-vinylguaiacol)
    cannot be represented as a fixed fraction of the p-coumaric flow without breaking carbon
    closure (a 9-carbon precursor cannot yield a 9-carbon product + CO2). So the ferulic branch
    reads its own state (``ferulic_acid``), its own Monod pair (``k_brett_decarb_ferulic``/
    ``K_hydroxycinnamic_ferulic`` — ratio-derived from Edlin et al. 1998's paired kinetics on the
    same enzyme, not independently re-estimated), and writes its own product pool
    (``vinylguaiacols``), which :class:`BrettVinylphenolReduction` also drains (Tchobanov et al.
    2008 confirm the same reductase acts on both vinylguaiacol and vinylphenol). Both branches
    share the *same* catalyst/gate (``X_brett · gate`` — the enzyme and its environmental
    sensitivity do not depend on which substrate it happens to be processing), computed once and
    passed to :func:`_decarboxylation_branch` for each substrate.

    Touches ``hydroxycinnamics``/``vinylphenols``/``ferulic_acid``/``vinylguaiacols``/``CO2``;
    reads the (constant, in pt1) catalyst ``X_brett`` plus the gate state (T, SO₂, solved pH).

    Returns a zero contribution before any pH work when undosed (``X_brett ≤ 0``) or when BOTH
    precursors are exhausted — structural value-isolability and no wasted ``brentq`` (the compile
    seam additionally disables the Process when Brett is not pitched, for tier isolability).
    Tier **speculative** (rate/gate magnitudes are estimates).
    """

    name = "brett_decarboxylation"
    tier = Tier.SPECULATIVE
    touches = ("hydroxycinnamics", "vinylphenols", "ferulic_acid", "vinylguaiacols", "CO2")
    #: The decarboxylase Monod pairs (p-coumaric and, D-55, ferulic) plus the shared Brett
    #: environmental-gate parameters. Their tiers cap the touched-pool output tiers via
    #: parameter-tier propagation (D-1). CO2 is already speculative (the always-on VDK
    #: decarboxylation), so this adds no new tier headline. pKa/SO₂ params read via ``acidbase``
    #: are omitted (all plausible; the Process is already speculative — the MalolacticConversion
    #: convention).
    reads: tuple[str, ...] = (
        "k_brett_decarb",
        "K_hydroxycinnamic",
        "k_brett_decarb_ferulic",
        "K_hydroxycinnamic_ferulic",
        *_BRETT_GATE_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE any pH solve: no catalyst, or BOTH precursors exhausted, ⇒ no
        # conversion, and an unsulfited/undosed run must not pay a per-RHS pH solve for a zero (or
        # γ(T)-only) result.
        x_brett = max(float(y[schema.slice("X_brett")][0]), 0.0) if "X_brett" in schema else 0.0
        if x_brett <= 0.0:
            return d
        hc_gpl = max(float(y[schema.slice("hydroxycinnamics")][0]), 0.0)
        fa_gpl = max(float(y[schema.slice("ferulic_acid")][0]), 0.0)
        if hc_gpl <= 0.0 and fa_gpl <= 0.0:
            return d

        ph = ph_of_state(y, schema, params) if _needs_ph_solve(y, schema) else 0.0
        gate = brett_environmental_gate(y, schema, params, ph)
        activity = x_brett * gate  # shared by both branches — same catalyst, same environment

        d_hc, d_vp, d_co2_pc = _decarboxylation_branch(
            hc_gpl,
            M_P_COUMARIC,
            M_VINYLPHENOL,
            params["k_brett_decarb"],
            params["K_hydroxycinnamic"],
            activity,
        )
        d_fa, d_vg, d_co2_fer = _decarboxylation_branch(
            fa_gpl,
            M_FERULIC,
            M_VINYLGUAIACOL,
            params["k_brett_decarb_ferulic"],
            params["K_hydroxycinnamic_ferulic"],
            activity,
        )

        d[schema.slice("hydroxycinnamics")] = d_hc
        d[schema.slice("vinylphenols")] = d_vp  # feeds the shared reservoir
        d[schema.slice("ferulic_acid")] = d_fa
        d[schema.slice("vinylguaiacols")] = d_vg  # feeds the ferulic-branch shared reservoir
        # p-coumaric C9 → vinylphenol C8 + CO2 C1, and ferulic C10 → vinylguaiacol C9 + CO2 C1
        # (both carbon-closing); CO2 sums both branches' contributions.
        d[schema.slice("CO2")] = d_co2_pc + d_co2_fer
        return d


def _reduction_branch(
    intermediate_gpl: float,
    intermediate_molar_mass: float,
    product_molar_mass: float,
    k: float,
    activity: float,
) -> tuple[float, float]:
    """One vinylphenol-reductase branch: intermediate → product, first-order (decision D-55).

    Shared by both the vinylphenol → ethylphenol branch and the vinylguaiacol → ethylguaiacol
    branch within :class:`BrettVinylphenolReduction` — the same reductase (Tchobanov et al. 2008
    confirm it acts on both substrates) reducing a different intermediate, differing only in molar
    masses. ``activity`` folds in the catalyst/gate term both branches share equally
    (``X_brett · gate``). Returns ``(d_intermediate, d_product)`` in g/L/h, exactly
    ``(0.0, 0.0)`` if there is nothing to reduce or ``activity ≤ 0``.
    """
    if intermediate_gpl <= 0.0 or activity <= 0.0:
        return 0.0, 0.0
    loss = k * activity * intermediate_gpl  # [g intermediate/L/h]
    return -loss, loss * product_molar_mass / intermediate_molar_mass


class BrettVinylphenolReduction(Process):
    """*Brettanomyces* vinylphenol reductase — two branches to ethylphenols (D-40, D-55).

    ``d(vinylphenols)/dt = −L``, ``d(ethylphenols)/dt = +L·M_ethylphenol/M_vinylphenol`` with the
    mass flux ``L = k_brett_reduction · X_brett · g_SO₂ · γ(T) · [vinylphenols]`` (first-order in
    the intermediate, catalyst-scaled, gated by the *same* Brett environment as the decarboxylase).
    A mole-for-mole C8 → C8 transfer (both pools weighted at their own carbon fraction), so it is
    carbon-neutral like the diacetyl → butanediol reduction (D-26); no sugar draw. This is the
    step *S. cerevisiae* lacks (even POF+), so it is what makes the shared ``vinylphenols``
    reservoir clear only when Brett is present — the emergent coupling (D-40 pt4).

    **Ferulic-branch reduction (decision D-55) — sourced enzyme identity, unsourced rate ratio.**
    Reduces ``vinylguaiacols → ethylguaiacols`` (a mole-for-mole C9 → C9 transfer, same form as
    vinylphenols → ethylphenols) alongside the p-coumaric branch above, via the shared
    :func:`_reduction_branch` helper. Tchobanov et al. 2008 directly confirm Brett's vinylphenol
    reductase acts on *both* 4-vinylguaiacol and 4-vinylphenol — so, unlike the decarboxylase
    branches (which have a sourced *relative* rate from Edlin et al. 1998), this reuses the *same*
    ``k_brett_reduction`` for both substrates: the enzyme-identity claim is sourced, but no paired
    vinylguaiacol-vs-vinylphenol rate comparison was found, so the reuse is a documented
    simplification (not yet a sourced ratio like the decarboxylase branches). Both branches share
    the same catalyst/gate ``activity = X_brett · gate``, computed once.

    Guarded like the decarboxylase: a zero contribution before any pH work when there is no
    catalyst (``X_brett ≤ 0``) or nothing to reduce in EITHER branch (value + perf isolability;
    compile-seam disable gives tier isolability). Tier **speculative** (rate magnitude estimate).
    """

    name = "brett_vinylphenol_reduction"
    tier = Tier.SPECULATIVE
    touches = ("vinylphenols", "ethylphenols", "vinylguaiacols", "ethylguaiacols")
    #: ``k_brett_reduction`` sets the reductase magnitude (reused for both branches, decision
    #: D-55 — see the class docstring for why); the shared Brett environmental-gate parameters
    #: throttle it by SO₂/temperature. Their tiers cap the touched-pool output tiers via
    #: parameter-tier propagation (D-1). Reads the constant catalyst ``X_brett``.
    reads: tuple[str, ...] = ("k_brett_reduction", *_BRETT_GATE_READS)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        vinylphenols = max(float(y[schema.slice("vinylphenols")][0]), 0.0)
        vinylguaiacols = max(float(y[schema.slice("vinylguaiacols")][0]), 0.0)
        if vinylphenols <= 0.0 and vinylguaiacols <= 0.0:  # nothing to reduce in either branch
            return d
        x_brett = max(float(y[schema.slice("X_brett")][0]), 0.0) if "X_brett" in schema else 0.0
        if x_brett <= 0.0:  # no Brett ⇒ no reduction (both reservoirs strand)
            return d

        ph = ph_of_state(y, schema, params) if _needs_ph_solve(y, schema) else 0.0
        gate = brett_environmental_gate(y, schema, params, ph)
        activity = x_brett * gate  # shared by both branches — same catalyst, same environment

        d_vp, d_ep = _reduction_branch(
            vinylphenols, M_VINYLPHENOL, M_ETHYLPHENOL, params["k_brett_reduction"], activity
        )
        d_vg, d_eg = _reduction_branch(
            vinylguaiacols,
            M_VINYLGUAIACOL,
            M_ETHYLGUAIACOL,
            params["k_brett_reduction"],
            activity,
        )

        d[schema.slice("vinylphenols")] = d_vp  # mole-for-mole C8→C8
        d[schema.slice("ethylphenols")] = d_ep
        d[schema.slice("vinylguaiacols")] = d_vg  # mole-for-mole C9→C9
        d[schema.slice("ethylguaiacols")] = d_eg
        return d


class BrettGrowth(Process):
    """*Brettanomyces* biomass growth — makes ``X_brett`` dynamic (decision D-40 pt2).

    The Brett twin of :class:`~fermentation.core.kinetics.malolactic.MalolacticGrowth`, with one
    load-bearing difference: **Brett grows on ETHANOL, not sugar**, so it can build up in a *dry,
    finished* wine — its actual niche as a post-AF/barrel spoiler. Because the decarboxylase and
    reductase are linear in ``X_brett``, a growing population makes the volatile-phenol spoilage
    **accelerate autocatalytically** over the months a barrel sits — the "it gets worse the longer
    you leave it" dynamic a constant catalyst (pt1) cannot produce.

    **Rate — the growth law, environmentally gated, carrying-capacity-braked.**

        dX_brett/dt = μ_max_brett · X_brett · aa/(K_aa_brett + aa) · E/(K_E_brett + E)
                                    · g_SO₂ · γ(T) · (1 − X_brett/K)

    Michaelis–Menten in the ``amino_acids`` pool (the nitrogen fuel, refilled post-AF by autolysis,
    D-34 — the honest way a dry wine feeds Brett) **and** in ethanol ``E`` (its carbon/energy
    source, below), scaled by the shared Brett gate. **No sugar Monod** (unlike MLF): Brett does not
    need sugar present, so growth is *not* self-limited to the sugar window — the ethanol Monod
    ``E/(K_E_brett + E)`` is ≈1 across the whole normal-wine ethanol range (``K_E_brett`` ~ a few
    g/L), so a finished dry wine feeds Brett fully while an unfermented must (E ≈ 0) does not. That
    ethanol availability, plus the amino-acid fuel and SO₂/temperature environment, is what makes
    ``pitch_brett`` into a high-ethanol post-AF wine grow rather than sit inert. (The ethanol Monod
    is also the *smooth shadow* of the ``E ≤ 0`` guard — it drives the rate to zero continuously as
    ``E → 0`` so the BDF Jacobian never straddles an on/off step during primary AF, mirroring how
    the amino-acid Monod shadows the ``aa`` guard and the ``(1 − X/K)`` brake shadows the cap.)

    **Ethanol-toxicity upper wall (decision D-58) — ethanol's dual role on one state variable.**
    Real Brett tolerance is bounded: Barata et al. 2008 measured growth at ~8% v/v ethanol and
    full arrest by ~14.5–15% (:func:`brett_ethanol_survival_factor`). The rate multiplies by this
    shared factor — 1 (no effect) up to ``brett_ethanol_toxicity_onset`` (above normal wine
    strength, so ordinary finished wine is unaffected), easing to 0 by
    ``brett_ethanol_toxicity_ceiling``. Combined with ``e_monod`` above, ethanol's net effect on
    growth is a *hump*: rises as a carbon source at low concentration, flat across the normal wine
    range, eases toward 0 near Barata's ceiling — the reconciliation of "ethanol as fuel" and
    "ethanol as toxin" on the same state variable, rather than two independent, potentially
    double-counting mechanisms. :class:`BrettEthanolToxicity` reuses the same factor (as
    ``1 −`` it) to drive a matching death term.

    **Carrying-capacity brake — the load-bearing difference from MLF (decision D-40 pt2).** MLF
    growth is *self-arresting*: its sugar Monod ``S/(K_s + S)`` vanishes as sugar is consumed and
    its gate carries an ethanol wall, so ``O. oeni`` cannot run away. Brett deliberately has
    **neither** (that is its dry-wine, ethanol-tolerant niche), so amino-acid Monod alone is *not*
    a brake — a barrel with a refilled amino-acid pool would grow ``X_brett`` exponentially with no
    ceiling, driving the amino-acid pool negative on a solver overshoot. So Brett carries an
    intrinsic **logistic carrying capacity** ``(1 − X_brett/K)`` (``K = brett_carrying_capacity``),
    the same lumped form as the D-30 yeast :class:`~fermentation.core.kinetics.carrying_capacity.\
    BiomassCarryingCapacity` — real Brett saturates at a finite cell density (nutrient/quorum
    limits). Linear ``1 − X/K`` (not a smoothed power) is deliberate: ``X_brett`` self-limits, so
    growth → 0 as ``X_brett → K`` and the state never gets driven past the wall — no derivative
    kink for the BDF solver, and no runaway. The factor is clamped ``≥ 0`` so a solver excursion
    ``X_brett > K`` cannot flip growth into a biomass/nitrogen *source*. Because ``K`` bounds
    ``X_brett`` small (realistic Brett biomass), the amino-acid draw rate ``ρ ∝ dX_brett`` stays
    small, so the pool depletes *smoothly* to a positive residual rather than overshooting negative
    — the brake is what makes the nitrogen ledger physically honest (it bounds X, MLF's environment
    does the same job there). This is intrinsic and always-on, **not** the opt-in isolable modifier
    D-30 is (which exists to depart from the Coleman anchor); it is a plain factor in the rate.

    **Conservation — nitrogen-anchored, carbon shortfall from ETHANOL (decision D-40, owner fork).**
    New biomass needs ``f_N·dX_brett`` nitrogen and ``f_C·dX_brett`` carbon. All the nitrogen comes
    from amino acids, consuming ``ρ = f_N·dX_brett/y_N`` of the pool (``y_N``/``y_C`` = arginine's
    nitrogen/carbon mass fractions). That arginine carries only ``ρ·y_C`` carbon — less than biomass
    needs (arginine is N-rich) — so the **shortfall** ``f_C·dX_brett − ρ·y_C`` is drawn from
    **ethanol** ``E`` (``d[E] −= shortfall / c_ethanol``, where ``c_ethanol`` is ethanol's carbon
    fraction, so exactly ``shortfall`` grams of carbon leave ``E``). This mirrors MLF's growth
    stoichiometry but sources the shortfall from ethanol instead of sugar — the mechanistic reason
    Brett thrives where the wine is dry. Carbon closes (``X_brett`` gains ``f_C·dX_brett`` =
    amino-acid carbon ``ρ·y_C`` + ethanol shortfall); nitrogen closes (``X_brett`` gains
    ``f_N·dX_brett`` = amino-acid nitrogen ``ρ·y_N``). The shortfall coefficient
    ``f_C − f_N·y_C/y_N`` is **structurally positive** (biomass C:N ≫ arginine's), so no clamp and
    no C⁰ kink. Touches ``(X_brett, amino_acids, E)`` — **not** ``S`` (Brett skips sugar) and
    **not** ``N`` (nitrogen from amino acids, no ammonium release; the D-38 anchoring choice).

    **v1 simplification (owned).** Real Brett *oxidizes* most of the ethanol it consumes to
    acetaldehyde/**acetic acid** (its volatile-acidity/"vinegar" side), assimilating only a fraction
    as biomass — this Process models the biomass-assimilation branch only (carbon-closing), and the
    acetic-acid overflow is a deferred pool. So the ethanol *drawdown* here is a lower bound on
    Brett's true ethanol consumption, and no volatile acidity is produced yet (decision D-40).

    Guards mirror MLF growth: a zero contribution before any pH work when there is no catalyst
    (``X_brett ≤ 0``), no amino-acid fuel, or no ethanol (the shortfall must never target an empty
    ``E``). Tier **speculative** (``μ_max_brett``/``K_aa_brett`` are author estimates and the
    bacterial-≈-yeast composition is a simplification).
    """

    name = "brett_growth"
    tier = Tier.SPECULATIVE
    #: Builds ``X_brett`` from the ``amino_acids`` pool, drawing the carbon shortfall from ``E``.
    #: Does NOT touch ``S`` (Brett grows in dry wine) or ``N`` (nitrogen from amino acids).
    touches = ("X_brett", "amino_acids", "E")
    #: ``mu_max_brett``/``K_aa_brett`` are the Brett growth rate + amino-acid half-saturation;
    #: ``biomass_N_fraction``/``biomass_C_fraction`` are the composition biomass is built (and
    #: weighted) at, so carbon/nitrogen close. The shared Brett gate params throttle growth by
    #: SO₂/temperature. Their tiers cap the ``X_brett``/``amino_acids``/``E`` output tiers (D-1).
    reads: tuple[str, ...] = (
        "mu_max_brett",
        "K_aa_brett",
        "K_E_brett",
        "brett_carrying_capacity",
        "biomass_N_fraction",
        "biomass_C_fraction",
        *_BRETT_GATE_READS,
        *_BRETT_ETHANOL_TOXICITY_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE any pH solve: no catalyst, no amino-acid fuel, or no ethanol (the
        # carbon shortfall must never target an empty E) ⇒ no growth, and no wasted pH solve.
        x_brett = max(float(y[schema.slice("X_brett")][0]), 0.0) if "X_brett" in schema else 0.0
        if x_brett <= 0.0:
            return d
        aa = max(float(y[schema.slice("amino_acids")][0]), 0.0) if "amino_acids" in schema else 0.0
        if aa <= 0.0:
            return d
        e = max(float(y[schema.slice("E")][0]), 0.0)
        if e <= 0.0:
            return d

        ph = ph_of_state(y, schema, params) if _needs_ph_solve(y, schema) else 0.0
        gate = brett_environmental_gate(y, schema, params, ph)

        # Logistic carrying-capacity brake (decision D-40 pt2): unlike MLF, Brett has no sugar
        # Monod / ethanol wall to self-arrest growth, so this is the ONLY ceiling. Clamp >= 0 so a
        # solver excursion X_brett > K cannot flip growth into a biomass/nitrogen source; growth
        # eases to 0 as X_brett -> K (no kink, no runaway, no negative amino-acid overshoot).
        brake = 1.0 - x_brett / params["brett_carrying_capacity"]
        if brake <= 0.0:  # at/above the carrying capacity, growth is fully shut down
            return d

        aa_monod = aa / (params["K_aa_brett"] + aa)
        # Ethanol-availability Monod (decision D-40 pt2). Brett grows ON ethanol, so growth scales
        # with how much ethanol is present: ~0 during early AF (E small) and ~1 in a finished wine
        # (E >> K_E_brett). This is ALSO the smooth "shadow" the hard `if e <= 0` guard needs: like
        # the amino-acid Monod and the (1 - X/K) brake, it drives the rate to zero SMOOTHLY as E→0,
        # so BDF's finite-difference Jacobian never straddles an on/off step at E=0 (which otherwise
        # corrupts the implicit solve into an autocatalytic X_brett blow-up during primary AF, when
        # E rises through zero — RK45/LSODA build no Jacobian and never saw it). K_E_brett is small
        # (~a few g/L), so this is ≈1 across the whole normal-wine ethanol range and only smooths
        # the near-zero crossing — Brett's dry-finished-wine niche is preserved.
        e_monod = e / (params["K_E_brett"] + e)
        # Ethanol-toxicity upper wall (decision D-58): reconciles ethanol's dual role on this SAME
        # state variable — carbon source at low concentration (e_monod, above) and toxin at high
        # concentration (this factor). Exactly 1 (no effect) across ordinary wine strength (Barata's
        # onset ~14% v/v is above typical finished-wine ethanol), easing to 0 by ~15% v/v — the
        # combined shape is the classic "hump": rises via e_monod, flat across the normal range,
        # eases toward 0 near Barata's measured growth ceiling. See brett_ethanol_survival_factor.
        e_toxicity_wall = brett_ethanol_survival_factor(e, params)
        dx_brett = (
            params["mu_max_brett"] * x_brett * aa_monod * e_monod * e_toxicity_wall * gate * brake
        )  # [g X_brett/L/h]
        if dx_brett <= 0.0:
            return d

        f_n = params["biomass_N_fraction"]
        f_c = params["biomass_C_fraction"]
        y_n = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
        y_c = carbon_mass_fraction(AMINO_ACID_SPECIES)
        # Nitrogen-anchored: consume the arginine that carries the new biomass's nitrogen. Its
        # carbon (rho*y_C) falls short of the biomass carbon demand (f_C*dx_brett) because arginine
        # is N-rich, so the positive shortfall is drawn from ETHANOL (not sugar — Brett's niche is
        # the dry wine). Both ledgers close exactly: X_brett gains f_N*dx_brett N (= rho*y_N) and
        # f_C*dx_brett C (= rho*y_C amino acid + ethanol shortfall). Shortfall coeff > 0
        # structurally (biomass C:N >> arginine's), so no clamp needed (decision D-40).
        rho = f_n * dx_brett / y_n  # [g arginine/L/h] consumed to supply the biomass nitrogen
        shortfall = f_c * dx_brett - rho * y_c  # [g C/L/h] biomass carbon not covered by arginine
        d[schema.slice("X_brett")] = dx_brett
        d[schema.slice("amino_acids")] = -rho
        # Remove `shortfall` grams of carbon from ethanol: g ethanol = g C / (g C/g ethanol).
        d[schema.slice("E")] = -shortfall / carbon_mass_fraction("ethanol")
        return d


class BrettDeath(Process):
    """*Brettanomyces* death — the **SO₂-driven** spoilage-yeast kill (decision D-40 pt3).

    The counterpart to :class:`BrettGrowth` that completes the Brett arc (pt1 pathway → pt2 growth →
    pt3 death): it moves viable ``X_brett`` into the non-viable ``X_brett_dead`` pool, so the
    spoilage population *declines* when the wine is sulfited and the volatile-phenol activities that
    scale with ``X_brett`` — :class:`BrettDecarboxylation` and :class:`BrettVinylphenolReduction` —
    wind down as Brett dies. That is the mechanism behind the winemaker's headline lever: a
    **molecular-SO₂ addition kills Brett**, so 4-EP production *stops rising* (the produced-only
    ``ethylphenols`` readout has no sink, so killing Brett simply halts the accrual — a clean,
    unconfounded lever, unlike the MLF diacetyl "lock-in" which kills both a source and a sink). A
    rack that draws Brett off the lees does the same physically (:data:`~fermentation.scenario.\
    compile._LEES_SLOTS` already carries ``X_brett``/``X_brett_dead`` since pt1) — the two ways a
    winemaker removes an established Brett contamination.

    **Rate — molecular-SO₂-driven, Arrhenius temperature (the :class:`~fermentation.core.kinetics.\
    malolactic.MalolacticDeath` form).**

        r_death = k_death_brett · X_brett · (1 − g_SO₂) · arrhenius(T, E_a_death_brett, T_ref)
        g_SO₂   = exp(−[SO₂]_molecular / molecular_so2_inhib_brett)

    The driver is ``1 − g_SO₂`` — the **same** molecular-SO₂ term the decarboxylase/reductase gate
    uses (D-22 antimicrobial readout, partitioned at the solved pH). Death is **exactly 0 without
    SO₂** (``g_SO₂ = 1``) and rises toward its Arrhenius ceiling as molecular SO₂ accumulates.
    Temperature enters through its own **Arrhenius** factor (enzymatic/chemical mortality, faster
    when warm — the shape autolysis and MLF death share), **not** the cardinal γ(T): γ(T) → 0 in the
    *cold*, which would make cold *kill* Brett, whereas cold in fact **preserves** it (part of why
    Brett is so hard to eradicate from a cool cellar). So warm accelerates the SO₂ kill and cold
    slows it toward dormancy.

    **Why SO₂ alone is the natural driver for Brett (contrast :class:`MalolacticDeath`).** MLF death
    had to *drop* an ethanol/pH toxicity driver because *O. oeni*'s Luong ethanol wall spuriously
    made bacteria "die" from ordinary post-AF ethanol (the D-39 crux). Brett has **no such wall** —
    the Brett gate (:func:`brett_environmental_gate`) carries no ethanol or pH term at all — because
    Brett is ethanol- and acid-tolerant. So "molecular SO₂ alone kills Brett" is not a
    correction of a confounder but the *directly correct* physics: the winemaker's ~0.5–0.8 mg/L
    molecular-SO₂ Brett-control target is the real-world expression of this term. Without SO₂ (or a
    rack) Brett persists indefinitely in v1 — the honest reflection of how tenacious a barrel Brett
    infection is; a slow benign-environment senescence is a deferred v2 refinement (see
    ``k_death_brett`` provenance).

    **Conservation — a carbon/nitrogen-neutral transfer (the D-13 pattern).** Since pt2 both
    ``X_brett`` and ``X_brett_dead`` are weighted in ``total_carbon``/``total_nitrogen`` at the
    *same* biomass fractions, so moving a gram from one to the other (``d[X_brett] = −r``,
    ``d[X_brett_dead] = +r``) is carbon- and nitrogen-neutral **by construction** — exactly like the
    yeast :class:`~fermentation.core.kinetics.inactivation.EthanolInactivation` transfer
    ``X → X_dead`` and the bacterial ``X_mlf → X_mlf_dead`` kill. No new conservation code, and no
    ``S``/``E`` draw; touches only ``(X_brett, X_brett_dead)``.

    **Isolability.** ``X_brett ≤ 0`` (undosed / unpitched) *or* no total SO₂ returns a zero
    contribution *before* the pH ``brentq`` — the ``so2_total ≤ 0`` guard is exact (death is
    identically 0 without SO₂, since ``1 − g_SO₂ = 0`` there), so a pitched-but-unsulfited run pays
    no per-RHS pH solve and its contribution is byte-for-byte zero. The compile seam enables this
    Process with the other **pitch-gated** Brett Processes (:data:`~fermentation.scenario.compile.\
    _BRETT_GATED_PROCESSES`) — it needs no amino acids (Brett dies whether or not it grew), so it
    is pitch-gated, not amino-acid-gated like :class:`BrettGrowth`. Consequence: on any *pitched*
    run ``X_brett``/``X_brett_dead`` report **speculative** (this enabled Process touches them) —
    honest, since a pitched population that can be sulfited has a speculative trajectory. Tier
    **speculative** (``k_death_brett``/``E_a_death_brett`` are estimates, and the SO₂-driven form —
    with the arrest-scale reused as kill-scale — is a modelling choice).
    """

    name = "brett_death"
    tier = Tier.SPECULATIVE
    #: Viable Brett leaves ``X_brett``; the same mass enters the non-viable ``X_brett_dead`` pool.
    #: Declaring both keeps the carbon/nitrogen-neutral transfer inside the ``touches`` contract.
    touches = ("X_brett", "X_brett_dead")
    #: ``k_death_brett`` sets the magnitude at full kill; ``E_a_death_brett``/``T_ref`` set the
    #: Arrhenius temperature shape; ``molecular_so2_inhib_brett`` is the SO₂ decay scale (reused
    #: from the metabolic gate — arrest-scale = kill-scale, a documented simplification). NOT the
    #: ``_BRETT_GATE_READS`` cardinals: death uses Arrhenius, not γ(T), so it must not read
    #: ``T_min/T_opt/T_max_brett`` (mirrors :class:`MalolacticDeath.reads`). Their tiers cap the
    #: ``X_brett``/``X_brett_dead`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_death_brett",
        "E_a_death_brett",
        "T_ref",
        "molecular_so2_inhib_brett",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards BEFORE any pH solve (the SO₂ partition reads the solved pH): no catalyst ⇒ no
        # Brett to kill; no total SO₂ ⇒ death is identically 0 (g_SO₂ = 1, so 1 − g_SO₂ = 0). The
        # so2_total ≤ 0 guard is EXACT, not an approximation, so an unsulfited pitched run pays no
        # per-RHS pH solve and its contribution is byte-for-byte zero (mirrors MalolacticDeath).
        x_brett = max(float(y[schema.slice("X_brett")][0]), 0.0) if "X_brett" in schema else 0.0
        if x_brett <= 0.0:
            return d
        total_so2 = float(y[schema.slice(SO2_STATE_KEY)][0]) if SO2_STATE_KEY in schema else 0.0
        if total_so2 <= 0.0:
            return d

        ph = ph_of_state(y, schema, params)
        molecular_so2 = molecular_so2_at_ph(y, schema, params, ph)
        g_so2 = math.exp(-molecular_so2 / params["molecular_so2_inhib_brett"])  # (0, 1]; 1 = no SO₂
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_death_brett"], params["T_ref"])
        # SO₂ drives death: 1 − g_SO₂ ∈ [0, 1) rises with molecular SO₂. Arrhenius (not the cardinal
        # γ(T)) carries temperature, so cold preserves rather than kills Brett (decision D-40 pt3).
        r_death = params["k_death_brett"] * x_brett * (1.0 - g_so2) * f_t  # [g X_brett/L/h]
        d[schema.slice("X_brett")] = -r_death
        d[schema.slice("X_brett_dead")] = r_death  # carbon/nitrogen-neutral: same biomass fractions
        return d


class BrettEthanolToxicity(Process):
    """*Brettanomyces* ethanol-toxicity death — the sourced alternative to a declined senescence
    (decision D-58).

    D-52 declined a generic ``BrettSenescence`` twin of :class:`~fermentation.core.kinetics.\
    malolactic.MalolacticSenescence`: no literature source shows Brett declining from elapsed time
    alone. D-58's follow-up research (two independent literature agents) re-confirmed that, but
    surfaced a real, DIFFERENT, sourced mechanism this Process closes: Barata et al. 2008 measured
    Brett growing normally up to ~14% v/v ethanol and fully arrested by ~14.5–15%, in closed-system
    model wine WITHOUT SO₂ — a threshold effect on *ethanol concentration*, not a function of
    elapsed time, so it is not the senescence D-52 declined.

    **Rate — the same ``1 − survival`` idiom as :class:`BrettDeath`'s ``1 − g_SO₂``.**

        r_death = k_death_brett · X_brett · (1 − survival(E)) · arrhenius(T, E_a_death_brett, T_ref)
        survival(E) = brett_ethanol_survival_factor(E, params)  — see that function's docstring

    Exactly 0 for ``E ≤ brett_ethanol_toxicity_onset`` (~110 g/L, ~14% v/v) — so ordinary
    finished wine (typically ~90–105 g/L) sees no contribution at all, not merely a small one —
    rising to the full ``k_death_brett`` rate by ``brett_ethanol_toxicity_ceiling`` (~118 g/L,
    ~15% v/v). **Reuses** ``k_death_brett``/``E_a_death_brett``/``T_ref`` rather than introducing
    new death-magnitude/temperature params — Barata's data was measured at one fixed 25 °C, so no
    independent activation energy is sourced, and reusing the existing SO₂-kill scale as the
    ethanol-kill scale mirrors :class:`BrettDeath`'s own documented "arrest-scale = kill-scale"
    simplification. Needs **no** SO₂ — unlike :class:`BrettDeath`, this fires on an unsulfited
    high-ethanol wine, which is exactly its point (the mechanism Barata measured requires no
    sulfite).

    **Scope limitation (documented, not silently assumed).** Barata's most-cited headline
    number — a 12% v/v ethanol, no-SO₂, 50-day population crash — is explicitly described in the
    source as bloom-on-trace-carbon *then* starvation-plus-ethanol-stress: a confounded result
    mixing substrate exhaustion (an unrelated, unmodelled mechanism here) with ethanol toxicity.
    12% v/v (~95 g/L) is *below* ``brett_ethanol_toxicity_onset``, so this Process alone predicts
    no decline at that concentration — it models only the distinct, unconfounded per-concentration
    boundary data (growth at ~8%, death onset ~14%, ceiling ~14.5–15%), not the
    starvation-confounded 12% result. A starvation-driven decline mechanism, if ever wanted, is a
    separate, not-yet-scoped addition.

    **Conservation — carbon/nitrogen-neutral transfer (the** :class:`BrettDeath` **pattern).** Same
    ``(X_brett, X_brett_dead)`` transfer at the same biomass fractions, so no new conservation code.

    **Isolability.** ``X_brett ≤ 0`` *or* ``E ≤ brett_ethanol_toxicity_onset`` returns a zero
    contribution — the ethanol guard is EXACT (survival ≡ 1 at or below onset), so an ordinary-
    strength wine (pitched or not) pays no cost beyond the guard check and its contribution is
    byte-for-byte zero. No pH solve is ever needed (no SO₂ term). Pitch-gated (added to
    :data:`~fermentation.scenario.compile._BRETT_GATED_PROCESSES` alongside :class:`BrettDeath` —
    Brett dies whether or not it grew), NOT amino-acid-gated. Tier **speculative** (the reused
    death-rate/temperature params are estimates; the survival-factor boundaries are sourced but the
    functional form/exponent are modelling choices).
    """

    name = "brett_ethanol_toxicity"
    tier = Tier.SPECULATIVE
    #: Same carbon/nitrogen-neutral transfer as :class:`BrettDeath` — both pools already weighted
    #: at the biomass fractions (D-40 pt3), so declaring both keeps the transfer inside ``touches``.
    touches = ("X_brett", "X_brett_dead")
    #: ``k_death_brett``/``E_a_death_brett``/``T_ref`` are REUSED from :class:`BrettDeath`
    #: (documented simplification — no independent ethanol-death rate/temperature law is sourced);
    #: the ethanol survival-factor params are new (D-58). Their tiers cap the
    #: ``X_brett``/``X_brett_dead`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_death_brett",
        "E_a_death_brett",
        "T_ref",
        *_BRETT_ETHANOL_TOXICITY_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Early guards: no catalyst ⇒ no Brett to kill; E at/below the toxicity onset ⇒ survival ≡ 1
        # so death is identically 0 (EXACT guard, not an approximation — ordinary-strength wine is
        # unaffected). No pH solve is ever needed (no SO₂ term in this Process).
        x_brett = max(float(y[schema.slice("X_brett")][0]), 0.0) if "X_brett" in schema else 0.0
        if x_brett <= 0.0:
            return d
        e = max(float(y[schema.slice("E")][0]), 0.0)
        if e <= params["brett_ethanol_toxicity_onset"]:
            return d

        survival = brett_ethanol_survival_factor(e, params)
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_death_brett"], params["T_ref"])
        # Ethanol toxicity drives death: 1 - survival rises from 0 at the onset to 1 at the ceiling.
        # Same Arrhenius temperature factor as BrettDeath (reused, decision D-58).
        r_death = params["k_death_brett"] * x_brett * (1.0 - survival) * f_t  # [g X_brett/L/h]
        d[schema.slice("X_brett")] = -r_death
        d[schema.slice("X_brett_dead")] = r_death  # carbon/nitrogen-neutral: same biomass fractions
        return d


class YeastPOFDecarboxylation(Process):
    """POF+ *S. cerevisiae* cinnamate decarboxylase — hydroxycinnamics → vinylphenols + CO2 (pt4).

    The *yeast* half of the volatile-phenol story (decision D-40 pt4). A **POF+**
    (phenolic-off-flavour-positive) primary-fermentation strain carries the cinnamate
    decarboxylase — the *same* enzyme :class:`BrettDecarboxylation` models — but **not** the
    vinylphenol reductase, so during AF it takes must ``hydroxycinnamics`` → ``vinylphenols`` + CO2
    and there it **stops**: it fills the shared reductase reservoir it cannot drain. If Brett is
    absent the ``vinylphenols`` *strand* (nothing reduces them to ethylphenols — the emergent
    yeast/Brett coupling this 3-pool design was chosen for, the α-acetolactate-reservoir parallel
    of D-26/D-31, DECISIONS D-40); if Brett arrives later, :class:`BrettVinylphenolReduction` drains
    that pre-filled reservoir, so a POF+ AF gives a subsequent Brett contamination a *head start*.

    **Same reaction as Brett's decarboxylase, different catalyst.** The chemistry — p-coumaric (9 C)
    → vinylphenol (8 C) + CO2 (1 C), carbon-closing mole-for-mole (9 = 8 + 1) — and the carbon
    routing are **identical** to :class:`BrettDecarboxylation` (it reuses ``M_P_COUMARIC``/
    ``M_VINYLPHENOL``/``M_CO2`` and :func:`_decarboxylation_branch`), so it touches only
    ``hydroxycinnamics``/``vinylphenols``/``ferulic_acid``/``vinylguaiacols``/``CO2`` and closes
    on the existing ledger with no new conservation code; when both this and Brett are active they
    draw the *same* ``hydroxycinnamics``/``ferulic_acid`` pools (both close 9 = 8 + 1 and
    10 = 9 + 1). The difference is the catalyst and its rate law:

        r = k_pof_decarb · flux(T) · S_total/(K_sugar_uptake + S_total)
              · [hc]/(K_hydroxycinnamic + [hc])

    — **flux-coupled** to active fermentation via the shared :func:`~fermentation.core.kinetics.\
    carbon_routing.fermentative_flux_shape` (the ester/α-acetolactate idiom, D-19/D-26), not scaled
    by ``X_brett``. POF decarboxylation is a *primary-fermentation* phenomenon: the flux term makes
    production track the yeast's fermentative activity and **stop at dryness** (``S → 0`` ⇒ rate 0),
    leaving whatever hydroxycinnamics/ferulic acid remain for a later Brett. The precursor Monod
    (shared ``K_hydroxycinnamic``/``K_hydroxycinnamic_ferulic``, decision D-55) rolls production off
    as each pool is consumed. **No Brett SO₂/temperature gate** — this is yeast metabolism during
    AF, before any Brett or sulfite lever applies.

    **Ferulic-acid branch (decision D-55).** Same fork as :class:`BrettDecarboxylation`'s: a genuine
    second precursor pool (``ferulic_acid``, 10 carbons), not a fixed-ratio split of
    ``hydroxycinnamics`` (9 carbons) — the two decarboxylate to different-carbon-count products, so
    only a real second pool stays carbon-exact. Uses ``k_pof_decarb_ferulic``/
    ``K_hydroxycinnamic_ferulic`` (the same Edlin et al. 1998 ratio applied to
    :class:`BrettDecarboxylation`'s ferulic branch — same enzyme family, Pad1/Fdc1), writing to
    ``vinylguaiacols`` (the ferulic-branch counterpart to ``vinylphenols``). Both branches share the
    *same* flux/Arrhenius activity (``flux · arrhenius(T, E_a_pof)``), computed once and passed to
    :func:`_decarboxylation_branch` for each substrate.

    **Temperature-dependent (v2, decision D-54) — net conversion FALLS with warmer fermentation.**
    The rate now carries its own Arrhenius factor ``arrhenius(T, E_a_pof, T_ref)`` — a real enzyme
    genuinely speeds up with warmth (``E_a_pof > 0``; Edlin et al. 1998 puts a homologous
    hydroxycinnamate decarboxylase's own thermal optimum at 40 °C, well above any wine/beer ferment
    temperature). But because this Process is **flux-coupled** (``r ∝ fermentative_flux_shape``,
    which itself rides the sugar-uptake Arrhenius ``E_a_uptake``), a warmer ferment also finishes
    *faster* — shrinking the time window POF has to act. The D-19 flux-coupled-byproduct ordering
    constraint (the same one governing ``E_a_esters``/``E_a_fusels``) says the **net** finished-wine
    conversion (the time-integral of ``r`` to dryness) scales as
    ``exp(-((E_a_pof − E_a_uptake)/R)·(1/T − 1/T_ref))`` — so whether total vinylphenol production
    rises or falls with temperature depends on ``E_a_pof`` **relative to** ``E_a_uptake``, not on
    ``E_a_pof`` alone. ``E_a_pof`` is set BELOW ``E_a_uptake`` (55,100 J/mol) — the *opposite*
    ordering from esters/fusels (set above ``E_a_uptake`` so their totals rise with T) — because the
    sourced real-world direction here is opposite theirs: brewing practice on this *exact* enzyme
    (Pad1/Fdc1 hydroxycinnamate decarboxylase, the same POF+ trait) is unambiguous — cooler
    wheat-beer fermentation retains more clove/4-vinylguaiacol character, warmer fermentation
    favours esters over phenolics. This supersedes D-40 pt4's "temperature-flat" choice: that was a
    reasoned v1 simplification (no pt4 behaviour needed POF's intrinsic direction, and the implicit
    ``E_a_pof = 0 < E_a_uptake`` already fell with warmer T through the flux-window effect alone, so
    the v1 *direction* was accidentally already right) — v2 replaces the implicit T-invariant
    placeholder with a genuine (if still speculative) intrinsic enzyme temperature term, sized to
    preserve and reinforce that same direction rather than risk reversing it.

    **Isolability — a separate opt-in strain, wholly independent of the Brett pitch.** POF+ is a
    binary *strain* trait, so it is enabled by its own compile-seam opt-in (``pof_positive``), NOT
    by ``brett_pitch_gpl`` — a POF+ ferment need not have Brett, and a POF-negative wine (the
    default) must make **no** vinylphenol. The Process is wired into the wine medium's own
    ``_POF_PROCESSES`` tuple and **disabled** unless the strain is opted in, so a default (POF−) run
    is byte-for-byte the validated core and the phenol slots keep their VALIDATED tier (the
    Brett-unpitched pattern; ``tier_of`` counts enabled, not nonzero, Processes). Returns a zero
    contribution before any work when there is no precursor or no fermentative flux (post-AF). On a
    POF+ run ``vinylphenols`` reports **speculative** (this enabled Process touches it) — honest —
    while ``ethylphenols`` stays VALIDATED at 0 unless Brett (the only reductase) is also present.

    Tier **speculative**: ``k_pof_decarb`` is an author estimate (no per-strain POF decarboxylase
    rate of this flux-coupled form is sourced), and the reaction is lumped as for Brett. Wine-only.
    """

    name = "yeast_pof_decarboxylation"
    tier = Tier.SPECULATIVE
    touches = ("hydroxycinnamics", "vinylphenols", "ferulic_acid", "vinylguaiacols", "CO2")
    #: ``k_pof_decarb``/``k_pof_decarb_ferulic`` (D-55) set the POF decarboxylase magnitude on
    #: each branch; ``K_sugar_uptake`` (shared with the fermentative-uptake flux this tracks) and
    #: ``K_hydroxycinnamic``/``K_hydroxycinnamic_ferulic`` (shared with Brett's decarboxylase —
    #: same whole-cell precursor affinities) shape them. ``E_a_pof``/``T_ref`` (D-54) give the
    #: decarboxylase its own Arrhenius temperature term (shared by both branches — see the class
    #: docstring for why this is set BELOW ``E_a_uptake`` rather than read from it here). Their
    #: tiers cap the touched-pool output tiers via parameter-tier propagation (D-1). CO2 is already
    #: speculative (the always-on VDK decarboxylation), so this adds no new tier headline.
    reads: tuple[str, ...] = (
        "k_pof_decarb",
        "K_hydroxycinnamic",
        "k_pof_decarb_ferulic",
        "K_hydroxycinnamic_ferulic",
        "K_sugar_uptake",
        "E_a_pof",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # No precursor in EITHER branch ⇒ nothing to decarboxylate; no fermentative flux (post-AF,
        # S→0, or crashed yeast) ⇒ POF production has stopped and the reservoirs simply wait for
        # Brett.
        hc_gpl = max(float(y[schema.slice("hydroxycinnamics")][0]), 0.0)
        fa_gpl = max(float(y[schema.slice("ferulic_acid")][0]), 0.0)
        if hc_gpl <= 0.0 and fa_gpl <= 0.0:
            return d
        flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d

        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_pof"], params["T_ref"])
        activity = flux * f_t  # shared by both branches — same catalyst, same temperature term

        d_hc, d_vp, d_co2_pc = _decarboxylation_branch(
            hc_gpl,
            M_P_COUMARIC,
            M_VINYLPHENOL,
            params["k_pof_decarb"],
            params["K_hydroxycinnamic"],
            activity,
        )
        d_fa, d_vg, d_co2_fer = _decarboxylation_branch(
            fa_gpl,
            M_FERULIC,
            M_VINYLGUAIACOL,
            params["k_pof_decarb_ferulic"],
            params["K_hydroxycinnamic_ferulic"],
            activity,
        )

        d[schema.slice("hydroxycinnamics")] = d_hc
        d[schema.slice("vinylphenols")] = d_vp  # fills the shared reductase reservoir
        d[schema.slice("ferulic_acid")] = d_fa
        d[schema.slice("vinylguaiacols")] = d_vg  # fills the ferulic-branch shared reservoir
        # p-coumaric C9 → vinylphenol C8 + CO2 C1, and ferulic C10 → vinylguaiacol C9 + CO2 C1
        # (both carbon-closing, same as Brett's decarboxylase); CO2 sums both branches.
        d[schema.slice("CO2")] = d_co2_pc + d_co2_fer
        return d
