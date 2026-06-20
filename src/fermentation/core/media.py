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

    X    biomass               g/L (dry cell weight)
    S    sugar                 g/L — a *vector*: 1 slot for wine, 3 for beer
    E    ethanol               g/L
    N    yeast-assimilable N   g/L
    T    temperature           K
    CO2  evolved CO2           g/L

Sugar is always a vector so beer's sequential glucose → maltose → maltotriose
uptake needs no structural change to also support wine's single lumped sugar.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import StateSchema, VarSpec


def _common_specs(sugar: VarSpec) -> list[VarSpec]:
    """The six state variables every medium shares, around its ``sugar`` spec.

    Order fixes the flat-array layout, so it is part of the contract: biomass
    first, then the (scalar or vector) sugar, then ethanol / nitrogen /
    temperature / evolved CO2.
    """
    return [
        VarSpec("X", "g/L", description="biomass (dry cell weight)"),
        sugar,
        VarSpec("E", "g/L", description="ethanol"),
        VarSpec("N", "g/L", description="yeast-assimilable nitrogen"),
        VarSpec("T", "K", description="temperature"),
        VarSpec("CO2", "g/L", description="evolved CO2"),
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
    """A named beverage family: its state schema plus the Processes that act on it.

    ``process_factories`` are zero-argument callables that each build one
    :class:`Process` (kinetics read their parameters at ``derivatives`` time, not
    construction time, so no arguments are needed here). The tuple is empty until
    the validated-core Processes land in Milestone 1; an empty set integrates to a
    constant trajectory, which is the honest "no kinetics yet" baseline.
    """

    name: str
    schema: StateSchema
    process_factories: tuple[Callable[[], Process], ...] = ()

    def build_process_set(self, *, strict: bool = False) -> ProcessSet:
        """Assemble this medium's Processes into a :class:`ProcessSet`."""
        return ProcessSet(
            self.schema, [factory() for factory in self.process_factories], strict=strict
        )


#: The registry of known media. Adding a beverage family = adding an entry here
#: (and, at the I/O boundary, an initial-composition vocabulary in
#: ``fermentation.scenario.compile``).
MEDIA: dict[str, Medium] = {
    "wine": Medium(name="wine", schema=wine_schema()),
    "beer": Medium(name="beer", schema=beer_schema()),
}


def get_medium(name: str) -> Medium:
    """Look up a registered :class:`Medium` by name."""
    try:
        return MEDIA[name]
    except KeyError:
        raise KeyError(f"Unknown medium {name!r}; known media: {sorted(MEDIA)}") from None
