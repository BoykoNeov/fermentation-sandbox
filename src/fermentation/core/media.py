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
``X_dead``, ``Gly``, ``Byp``, ``esters`` and ``fusels`` are *produced-only* pools —
always zero at pitch and only accumulated by the kinetics — so they declare a default
initial of 0 (`VarSpec.default`) and need not be named at every initial-condition call
site. The ``esters``/``fusels`` pools are added in the Milestone-2 byproducts beat as
empty slots; the Processes that fill them are wired below. Their (trace) carbon is left
out of ``total_carbon`` under **interim accounting (b)** to avoid double-counting the
``Byp`` succinic sink (which already books the higher alcohols); routing it from sugar
and weighting the pools (the agreed **option (a)**) is planned for a future session —
see ``docs/plans/milestone-2-tasks.md`` and the ``kinetics.byproducts`` module
docstring. Settled as decision D-19 once (a) lands.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from fermentation.core.kinetics import (
    ArrheniusTemperature,
    EsterSynthesis,
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
    ]


def wine_schema() -> StateSchema:
    """Wine state layout: a single lumped fermentable sugar slot."""
    return StateSchema(_common_specs(VarSpec("S", "g/L", description="fermentable sugar")))


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
#: *separate* tuple from the validated-core primary set so the speculative beat
#: stays isolable — disabling these by name leaves the core byte-for-byte (prime
#: directive #3). They are additive, produced-only Processes that touch only the
#: ``esters``/``fusels`` pools, so wiring them in does not perturb the §2.2 trio or
#: carbon conservation — under **interim accounting (b)** their carbon is left out of
#: ``total_carbon`` (booked against ``Byp``); the agreed **option (a)** (route from
#: sugar, weight the pools) is planned next session. See D-19 / milestone-2-tasks.md.
_BYPRODUCT_PROCESSES: tuple[Callable[[], Process], ...] = (
    EsterSynthesis,
    FuselAlcoholsEhrlich,
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
