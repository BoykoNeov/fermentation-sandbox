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

Sugar is always a vector so beer's sequential glucose → maltose → maltotriose
uptake needs no structural change to also support wine's single lumped sugar.
``X_dead``, ``Gly``, ``Byp``, ``esters`` and ``fusels`` start at zero at pitch and are
only accumulated by the kinetics, so they declare a default initial of 0
(`VarSpec.default`) and need not be named at every initial-condition call site. The
``esters``/``fusels`` pools are filled by the Tier-2 byproduct Processes wired below.
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
    ]


def wine_schema() -> StateSchema:
    """Wine state layout: a single lumped fermentable sugar slot, plus the wine-only
    charge-active acid + strong-cation slots the pH charge-balance solver reads
    (decision D-18).

    These four slots are appended to ``wine_schema`` only (not ``_common_specs``), so
    ``beer_schema`` is untouched — beer's pH is a phosphate-buffered different acid
    system with no sourced data yet, explicitly deferred. ``default=0.0`` is
    load-bearing: existing wine scenarios/tests that name no acids still compile (all
    four → 0), and with acids and cation at 0 the slots are inert — they contribute 0
    to every conservation sum, so the validated core and its tests are untouched (prime
    directive #3). No Process touches them in D-18, so their derivatives are 0 and the
    trajectory is constant; they exist so a future MLF Process can deplete ``malic`` /
    grow ``lactic`` and so the charge balance and ``total_carbon`` can read them. pH is
    simply not meaningful for a no-acid scenario and is only *computed* when requested
    (``fermentation.analysis``). ``cation_charge`` is a charge density (mol⁺/L), not a
    mass concentration — state is already heterogeneous (``T`` in K) — back-solved from
    the scenario's measured ``initial_ph`` at compile and held constant (D-18).
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


#: The registry of known media. Adding a beverage family = adding an entry here
#: (and, at the I/O boundary, an initial-composition vocabulary in
#: ``fermentation.scenario.compile``).
MEDIA: dict[str, Medium] = {
    "wine": Medium(
        name="wine",
        schema=wine_schema(),
        process_factories=_PRIMARY_FERMENTATION_PROCESSES + _BYPRODUCT_PROCESSES,
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
