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
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.kinetics.arrhenius import arrhenius_factor
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
        dx_brett = (
            params["mu_max_brett"] * x_brett * aa_monod * e_monod * gate * brake
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
