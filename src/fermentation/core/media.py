"""Medium definitions — the named state layouts the validated core models.

A *medium* (wine, beer, …) fixes two things the rest of the engine builds on:

  * its :class:`~fermentation.core.state.StateSchema` — how many sugar slots, in
    what order, alongside biomass / ethanol / nitrogen / temperature / CO2; and
  * the Processes that act on that state (the kinetics).

Both are data here, not physics. This module declares *what a wine or beer state
looks like* and *which Processes apply*; the Processes themselves are ordinary
:class:`~fermentation.core.process.Process` subclasses elsewhere in the core, and
the industry-unit conversion boundary lives in ``fermentation.scenario.compile``.
Keeping the layout in the core gives the Processes (which reference variable names
like ``"S"`` and ``"N"``) and the scenario→core compile seam a single source of
truth to agree on.

The shared variables (decisions D-B / D-4):

    X      viable biomass        g/L (dry cell weight)
    S      sugar                 g/L — a *vector*: 1 slot for wine, 3 for beer
    E      ethanol               g/L
    N      yeast-assimilable N   g/L
    T      temperature           K
    CO2    evolved CO2           g/L
    X_dead ethanol-inactivated   g/L (non-viable biomass; carbon/nitrogen still
                                 counted, but no longer catalytic — decision D-13)
    Gly    glycerol              g/L (realised-yield byproduct sink — decision D-16)
    Byp    minor byproducts      g/L (lumped organic acids / higher alcohols,
                                 carbon-accounted as succinic acid — decision D-16)
    esters esters                g/L (aroma byproducts; lumped produced-only pool)
    fusels fusel/higher alcohols g/L (Ehrlich pathway; lumped produced-only pool)
    acetolactate α-acetolactate  g/L (vicinal-diketone precursor reservoir — decision D-26)
    diacetyl diacetyl (VDK)      g/L (buttery off-note; produced then reabsorbed — D-26)
    butanediol 2,3-butanediol    g/L (flavour-inactive diacetyl-reduction product — D-26)

Sugar is always a vector so beer's sequential glucose → maltose → maltotriose
uptake needs no structural change to also support wine's single lumped sugar.
``X_dead``, ``Gly``, ``Byp``, ``esters``, ``fusels`` and the VDK pools
(``acetolactate``/``diacetyl``/``butanediol``) start at zero at pitch and are only
accumulated by the kinetics, so they declare a default initial of 0
(`VarSpec.default`) and need not be named at every initial-condition call site. The
``esters``/``fusels`` pools are filled by the Tier-2 byproduct Processes wired below;
the three VDK pools by the diacetyl-pathway Processes (decision D-26).
Under **decision D-19 (option a1)** those Processes route the aroma carbon *out of
``S``* and ``total_carbon`` weights the pools (as ethyl acetate / isoamyl alcohol), so
``esters``/``fusels`` are real carbon-accounted state alongside ``Gly``/``Byp`` — not
diagnostic re-expressions. The former ``Byp`` double-count (it once lumped higher
alcohols) is resolved by carving them out of ``Y_byproduct_sugar``; the draw touches
only ``S`` (never ``E``/``CO2``), so turning the byproducts on perturbs the core only
by the trace sugar they consume. See ``docs/plans/milestone-2-tasks.md`` and the
``kinetics.byproducts`` module docstring.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from fermentation.core.kinetics import (
    ArrheniusTemperature,
    EsterSynthesis,
    EsterVolatilization,
    EthanolInactivation,
    FuselAlcoholsEhrlich,
    GrowthNitrogenLimited,
    MalolacticConversion,
    SugarUptakeToEthanolCO2,
)
from fermentation.core.process import Process, ProcessSet, RateModifier
from fermentation.core.state import StateSchema, VarSpec


def _common_specs(sugar: VarSpec) -> list[VarSpec]:
    """The six state variables every medium shares, around its ``sugar`` spec.

    Order fixes the flat-array layout, so it is part of the contract: biomass
    first, then the (scalar or vector) sugar, then ethanol / nitrogen /
    temperature / evolved CO2.
    """
    return [
        VarSpec("X", "g/L", description="viable biomass (dry cell weight)"),
        sugar,
        VarSpec("E", "g/L", description="ethanol"),
        VarSpec("N", "g/L", description="yeast-assimilable nitrogen"),
        VarSpec("T", "K", description="temperature"),
        VarSpec("CO2", "g/L", description="evolved CO2"),
        VarSpec(
            "X_dead", "g/L", default=0.0, description="ethanol-inactivated (non-viable) biomass"
        ),
        VarSpec("Gly", "g/L", default=0.0, description="glycerol (realised-yield byproduct)"),
        VarSpec(
            "Byp",
            "g/L",
            default=0.0,
            description="minor byproducts (organic acids/higher alcohols; succinic-equivalent)",
        ),
        VarSpec(
            "esters",
            "g/L",
            default=0.0,
            description="esters (fermentation aroma; lumped produced-only pool)",
        ),
        VarSpec(
            "fusels",
            "g/L",
            default=0.0,
            description="fusel / higher alcohols (Ehrlich pathway; lumped produced-only pool)",
        ),
        VarSpec(
            "esters_gas",
            "g/L",
            default=0.0,
            description="esters lost to the headspace by CO2 stripping (volatilized; "
            "carbon-bookkeeping pool, decision D-20)",
        ),
        VarSpec(
            "acetolactate",
            "g/L",
            default=0.0,
            description="alpha-acetolactate — vicinal-diketone precursor reservoir "
            "(spontaneously decarboxylates to diacetyl; decision D-26)",
        ),
        VarSpec(
            "diacetyl",
            "g/L",
            default=0.0,
            description="diacetyl (2,3-butanedione) — buttery vicinal diketone; "
            "produced then yeast-reabsorbed (the diacetyl rest, decision D-26)",
        ),
        VarSpec(
            "butanediol",
            "g/L",
            default=0.0,
            description="2,3-butanediol — flavour-inactive terminal product of "
            "diacetyl reduction by viable yeast (decision D-26)",
        ),
    ]


def wine_schema() -> StateSchema:
    """Wine state layout: a single lumped fermentable sugar slot, plus the wine-only
    charge-active acid + strong-cation slots the pH charge-balance solver reads
    (decision D-18), the free-SO₂ pool the molecular-SO₂ readout reads (decision D-22),
    and the ``X_mlf`` malolactic-catalyst slot (decision D-23).

    These six slots are appended to ``wine_schema`` only (not ``_common_specs``), so
    ``beer_schema`` is untouched — beer's pH is a phosphate-buffered different acid
    system with no sourced data yet, explicitly deferred. ``default=0.0`` is
    load-bearing: existing wine scenarios/tests that name no acids still compile (all
    six → 0), and with acids, cation, SO₂ and ``X_mlf`` at 0 the slots are inert — they
    contribute 0 to every conservation sum, so the validated core and its tests are
    untouched (prime directive #3). The acid/cation/SO₂ slots have no Process touching
    them in D-18/D-22; under D-23 :class:`~fermentation.core.kinetics.malolactic.\
    MalolacticConversion` depletes ``malic`` / grows ``lactic`` / evolves ``CO2`` *only
    when ``X_mlf`` is dosed* (and is disabled at the compile seam otherwise), so undosed
    wine runs keep a constant acid trajectory. ``X_mlf`` itself is inert in v1 (no Process
    grows or kills it) and carbon-free in ``total_carbon`` (constant ⇒ 0 drift); it enters
    the carbon ledger only when the later MLF-growth beat lands. pH is
    simply not meaningful for a no-acid scenario and is only *computed* when requested
    (``fermentation.analysis``). ``cation_charge`` is a charge density (mol⁺/L), not a
    mass concentration — state is already heterogeneous (``T`` in K) — back-solved from
    the scenario's measured ``initial_ph`` at compile and held constant (D-18).
    ``so2_free`` (g/L of SO₂-equivalent) is a dosed input read by ``acidbase.molecular_so2``
    to partition the antimicrobial molecular fraction at the solved pH; it is **not** in
    the charge balance (readout-only, D-22) and is carbon-free, so it leaves both pH and
    ``total_carbon`` unchanged.
    """
    specs = _common_specs(VarSpec("S", "g/L", description="fermentable sugar"))
    specs += [
        VarSpec("tartaric", "g/L", default=0.0, description="tartaric acid (must input; diprotic)"),
        VarSpec(
            "malic",
            "g/L",
            default=0.0,
            description="L-malic acid (must input; diprotic; MLF substrate)",
        ),
        VarSpec(
            "lactic", "g/L", default=0.0, description="L-lactic acid (produced-only; MLF product)"
        ),
        VarSpec(
            "cation_charge",
            "mol/L",
            default=0.0,
            description="net strong-cation charge (K+-dominant), constant; "
            "back-solved from initial_ph (D-18)",
        ),
        VarSpec(
            "so2_free",
            "g/L",
            default=0.0,
            description="free SO2 (molecular+bisulfite+sulfite, as SO2); dosed input, "
            "inert; pH-driven molecular-fraction readout (D-22)",
        ),
        VarSpec(
            "X_mlf",
            "g/L",
            default=0.0,
            description="Oenococcus oeni biomass — dosed-but-inert MLF catalyst "
            "(scales the malolactic rate; no growth/death in v1, decision D-23)",
        ),
    ]
    return StateSchema(specs)


def beer_schema() -> StateSchema:
    """Beer state layout: three sugars consumed sequentially.

    Glucose is taken up first, then maltose, then maltotriose — the order the
    ``components`` tuple records and the sugar-uptake Process will honour.
    """
    return StateSchema(
        _common_specs(
            VarSpec(
                "S",
                "g/L",
                size=3,
                description="fermentable sugars (sequential uptake)",
                components=("glucose", "maltose", "maltotriose"),
            )
        )
    )


@dataclass(frozen=True)
class Medium:
    """A named beverage family: its state schema plus the kinetics that act on it.

    ``process_factories`` are zero-argument callables that each build one additive
    :class:`Process`; ``modifier_factories`` likewise build the multiplicative
    :class:`RateModifier` objects (ethanol inhibition, Arrhenius temperature
    dependence) that scale those Processes. Both are *factories* rather than shared
    instances so every ``build_process_set`` call gets fresh objects — two media (or
    two runs) never share a mutable Process/modifier. Kinetics read their parameters
    at ``derivatives``/``factor`` time, not construction time, so the factories need
    no arguments.

    An empty pair of tuples integrates to a constant trajectory — the honest
    "no kinetics" baseline a bare :class:`Medium` still provides.
    """

    name: str
    schema: StateSchema
    process_factories: tuple[Callable[[], Process], ...] = ()
    modifier_factories: tuple[Callable[[], RateModifier], ...] = field(default=())

    def build_process_set(self, *, strict: bool = False) -> ProcessSet:
        """Assemble this medium's Processes and modifiers into a :class:`ProcessSet`."""
        return ProcessSet(
            self.schema,
            [factory() for factory in self.process_factories],
            modifiers=[factory() for factory in self.modifier_factories],
            strict=strict,
        )


#: The validated-core primary-fermentation kinetics, as zero-argument factories.
#: Wine and beer share the *same* mechanism set — biomass growth, fermentative
#: sugar uptake, and ethanol-driven cell inactivation (the cumulative viability
#: brake that sets the fermentation timescale, Coleman 2007), with per-rate
#: Arrhenius temperature dependence scaling growth and uptake. The only structural
#: difference between the two media is the sugar vector (1 slot vs 3): beer's
#: sequential glucose→maltose→maltotriose uptake is handled *inside*
#: :class:`~fermentation.core.kinetics.uptake.SugarUptakeToEthanolCO2` via catabolite
#: repression, so it needs no extra Process here.
#:
#: The instantaneous Luong ethanol wall (``EthanolInhibition``) is **not** wired in:
#: the cumulative inactivation Process is the mechanistically-correct ethanol brake,
#: and stacking an instantaneous wall on top would double-count ethanol toxicity
#: (decision D-13). The class is retained for optional/strain use.
_PRIMARY_FERMENTATION_PROCESSES: tuple[Callable[[], Process], ...] = (
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
    EthanolInactivation,
)
_PRIMARY_FERMENTATION_MODIFIERS: tuple[Callable[[], RateModifier], ...] = (
    ArrheniusTemperature.for_growth,
    ArrheniusTemperature.for_uptake,
)

#: Tier-2 temperature-/metabolism-driven aroma byproducts (Milestone 2, decision
#: D-18/D-19): ester synthesis and Ehrlich-pathway fusel alcohols. Kept as a
#: *separate* tuple from the validated-core primary set so the speculative beat stays
#: **isolable** (prime directive #3): building a ProcessSet without this tuple is the
#: pure validated core. Under D-19 (option a1) they route aroma carbon out of ``S``
#: and ``total_carbon`` weights the ``esters``/``fusels`` pools, so they no longer
#: leave the core byte-for-byte when enabled — turning them on draws a *trace* of
#: sugar (~0.2 % of ``S0``), perturbing only ``dS`` (never ``dE``/``dCO2``). Carbon
#: still closes to machine precision with them on, and the §2.2 trio stays in band.
#: See D-19 / milestone-2-tasks.md.
#:
#: :class:`EsterVolatilization` (decision D-20) is the gas-stripping sink that moves
#: liquid ``esters`` into the bookkeeping ``esters_gas`` headspace pool as CO2 sparges
#: the must — the physics behind wine's "warmer ⇒ *less* liquid ester" (Rollero 2014):
#: with ``E_a_ester_volatil`` set *per medium* it is held **above** ``E_a_esters`` for
#: wine (stripping outruns synthesis, liquid esters fall with T) and **below** it for
#: beer (synthesis dominates, esters rise with T — de Andrés-Toro). The transfer is
#: carbon-neutral (``esters`` → ``esters_gas``, both booked as ethyl acetate), so it is
#: in this isolable tuple too and ``total_carbon`` still closes to machine precision.
_BYPRODUCT_PROCESSES: tuple[Callable[[], Process], ...] = (
    EsterSynthesis,
    FuselAlcoholsEhrlich,
    EsterVolatilization,
)

#: Malolactic fermentation (wine-only, decision D-23): the *Oenococcus oeni* malate →
#: lactate + CO2 conversion, the first RHS consumer of the D-18 pH solver and the D-22
#: molecular-SO₂ readout. Kept as its own tuple so it stays **isolable** (prime directive
#: #3): the conversion contributes zero before the pH solve whenever ``X_mlf`` is undosed
#: (structural *value* isolability), and the compile seam *disables* it when MLF is not
#: pitched so the inert ``malic``/``lactic`` slots keep their VALIDATED tier rather than
#: being dragged to speculative by an enabled-but-zero Process (*tier* isolability —
#: ``ProcessSet.tier_of`` counts enabled, not nonzero, Processes). Wine-only: beer has no
#: ``malic``/``lactic`` slots, so it is never wired there.
_MLF_PROCESSES: tuple[Callable[[], Process], ...] = (MalolacticConversion,)


#: The registry of known media. Adding a beverage family = adding an entry here
#: (and, at the I/O boundary, an initial-composition vocabulary in
#: ``fermentation.scenario.compile``).
MEDIA: dict[str, Medium] = {
    "wine": Medium(
        name="wine",
        schema=wine_schema(),
        process_factories=_PRIMARY_FERMENTATION_PROCESSES + _BYPRODUCT_PROCESSES + _MLF_PROCESSES,
        modifier_factories=_PRIMARY_FERMENTATION_MODIFIERS,
    ),
    "beer": Medium(
        name="beer",
        schema=beer_schema(),
        process_factories=_PRIMARY_FERMENTATION_PROCESSES + _BYPRODUCT_PROCESSES,
        modifier_factories=_PRIMARY_FERMENTATION_MODIFIERS,
    ),
}


def get_medium(name: str) -> Medium:
    """Look up a registered :class:`Medium` by name."""
    try:
        return MEDIA[name]
    except KeyError:
        raise KeyError(f"Unknown medium {name!r}; known media: {sorted(MEDIA)}") from None
