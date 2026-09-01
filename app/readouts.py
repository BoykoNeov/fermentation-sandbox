"""What the console is allowed to draw, and what it must say while drawing it.

A **readout** is one line on a chart: where the numbers come from, what unit they are in,
how far they can be trusted, and — the part that matters — any warning the project has
recorded about it. A warning here is not a footnote. It is printed beside the chart every
single time, because a line the project already knows is misleading is a line someone will
quote the moment it appears without one.

Two of the warnings come straight from the engine's own notes. Total acidity climbs over a
ferment here where a real one holds steady or falls, and bound sulfur dioxide is counted
short because only one of the compounds that binds it is modelled. Neither is a fault of
this screen, and neither may be quietly smoothed over by it.

Nothing in this file draws anything. It describes the readouts; the app and the written
report both read the same description, so they cannot end up disagreeing.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

import fermentation.analysis as analysis
from app.runner import RunResult
from fermentation.core.acidbase import ALL_ACIDS
from fermentation.core.state import FloatArray
from fermentation.core.tiers import Tier, combine
from fermentation.units import abv_from_ethanol, kelvin_to_celsius, sugar_gpl_to_brix


@dataclass(frozen=True)
class Readout:
    """One plottable quantity, with its unit, its tier sources and its caveat."""

    key: str
    label: str
    unit: str
    #: State variables this readout is computed from. Its reported tier is the lowest of
    #: theirs — the same min-combine rule the engine uses everywhere else.
    sources: tuple[str, ...]
    #: Turns a finished run into the series. Given the run so it can reach the resolved
    #: parameters for the readouts (pH, SO2) that need them.
    compute: Callable[[RunResult], FloatArray]
    #: Printed on every chart this readout appears on. Never suppressed.
    caveat: str = ""
    #: Media this readout exists for. A wine-only readout is simply absent from a beer run.
    media: tuple[str, ...] = ("wine", "beer")
    #: Slot index for a vector variable (sugar). ``None`` ⇒ scalar or already reduced.
    slot: int | None = None

    def available(self, result: RunResult) -> bool:
        """True when this run's schema carries at least one of the readout's sources.

        "At least one" rather than "all", because the derived readouts list every slot that
        *could* feed them across both media — pH is solved from whichever acids the medium
        actually has — and a wine-only acid missing from a beer schema is not a reason to
        hide the beer's pH.
        """
        return any(name in result.schema for name in self.sources)

    def inert(self, result: RunResult) -> bool:
        """True when no mechanism in the run writes to any of this readout's sources.

        An inert readout is a flat line at whatever the compile seam seeded — an unhopped
        beer's bitterness, browning in a run that never ages. It has no confidence tier
        because no mechanism made a claim about it, and it must never be badged
        ``validated``, which is what the engine's own empty-list tier combine would return.
        """
        touched = result.touched_variables()
        return not any(name in touched for name in self.sources)

    def tier(self, result: RunResult) -> Tier | None:
        """Lowest tier among the sources a mechanism actually drives. ``None`` when inert.

        **Test the result with ``is None``, never with ``if tier:``.** ``Tier.SPECULATIVE``
        is the enum's zero and is therefore falsy, so a truthiness test silently reports the
        engine's least confident tier as "no tier at all" — the one confusion this whole
        return type exists to prevent.
        """
        if self.inert(result):
            return None
        touched = result.touched_variables()
        tiers = [
            result.traj.tier_map[n]
            for n in self.sources
            if n in touched and n in result.traj.tier_map
        ]
        return combine(tiers) if tiers else None

    def series(self, result: RunResult) -> FloatArray:
        return np.asarray(self.compute(result), dtype=np.float64)


def _state(name: str, slot: int | None = None) -> Callable[[RunResult], FloatArray]:
    def get(result: RunResult) -> FloatArray:
        block = np.atleast_2d(np.asarray(result.traj.series(name), dtype=np.float64))
        row: FloatArray = block[slot if slot is not None else 0]
        return row

    return get


def _sugar_total(result: RunResult) -> FloatArray:
    block = np.atleast_2d(np.asarray(result.traj.series("S"), dtype=float))
    return np.asarray(block.sum(axis=0), dtype=float)


def _scaled(name: str, factor: float) -> Callable[[RunResult], FloatArray]:
    def get(result: RunResult) -> FloatArray:
        return np.asarray(result.traj.series(name), dtype=float) * factor

    return get


def _analysis(
    fn: Callable[..., FloatArray], *, needs_params: bool = True, factor: float = 1.0
) -> Callable[[RunResult], FloatArray]:
    def get(result: RunResult) -> FloatArray:
        traj = result.traj.as_trajectory()
        out = fn(traj, result.param_values) if needs_params else fn(traj)
        return np.asarray(out, dtype=float) * factor

    return get


#: Trust the day-0 value and nothing after it. Taken from the engine's own note on this
#: series: the whole by-product pool is counted as if it were acid, so the line climbs where
#: a real ferment's would not.
_TA_CAVEAT = (
    "Trust the starting value, not the finish. Total acidity is what a lab measures by "
    "titrating a sample. This model counts everything the yeast makes as though it were "
    "acid, so the line climbs by roughly 3-4 g/L across a ferment. A real ferment holds "
    "steady or drifts slightly down. Day zero is sound; after that it reads too high and "
    "moves the wrong way."
)

_BOUND_SO2_CAVEAT = (
    "This is an under-count, and knowingly so. Several compounds in wine tie up sulfur "
    "dioxide; only one of them is modelled here. The real bound figure is higher than this "
    "line, which means the free SO2 line beside it is more optimistic than a lab would be."
)

_GAS_CAVEAT = (
    "This is the running total that has escaped into the headspace, not how much is left "
    "in the liquid."
)

#: Every slot that enters the charge balance pH is solved from — read from the engine's own
#: acid registry (both media's, filtered per run by ``Readout.available``) plus the three
#: non-acid charge carriers. Hard-coding this list would put the app one new acid behind the
#: engine and quietly mis-report pH's confidence from then on.
_CHARGE_ACTIVE: tuple[str, ...] = (
    *sorted(ALL_ACIDS),
    "Byp",
    "N",
    "cation_charge",
    "nitrogen_charge_excess",
)


READOUTS: tuple[Readout, ...] = (
    # -- the ferment itself ----------------------------------------------------------
    Readout("sugar", "Sugar (total)", "g/L", ("S",), _sugar_total),
    Readout("glucose", "Glucose", "g/L", ("S",), _state("S", 0), media=("beer",), slot=0),
    Readout("maltose", "Maltose", "g/L", ("S",), _state("S", 1), media=("beer",), slot=1),
    Readout("maltotriose", "Maltotriose", "g/L", ("S",), _state("S", 2), media=("beer",), slot=2),
    Readout("ethanol", "Ethanol", "g/L", ("E",), _state("E")),
    Readout(
        "abv",
        "Alcohol",
        "% ABV",
        ("E",),
        lambda r: np.array([abv_from_ethanol(v) for v in _state("E")(r)]),
    ),
    Readout(
        "brix",
        "Apparent sugar",
        "degrees Brix",
        ("S",),
        lambda r: np.array([sugar_gpl_to_brix(v) for v in _sugar_total(r)]),
        caveat="Worked out from the modelled sugar using a standard density, so it is roughly "
        "what a hydrometer would read — not a reading anyone took.",
    ),
    Readout("biomass", "Yeast (viable)", "g/L", ("X",), _state("X")),
    Readout("biomass_dead", "Yeast (dead)", "g/L", ("X_dead",), _state("X_dead")),
    Readout("co2", "CO2 evolved", "g/L", ("CO2",), _state("CO2")),
    Readout("glycerol", "Glycerol", "g/L", ("Gly",), _state("Gly")),
    Readout(
        "temperature",
        "Temperature",
        "degrees C",
        ("T",),
        lambda r: np.array([kelvin_to_celsius(v) for v in _state("T")(r)]),
    ),
    # -- nitrogen --------------------------------------------------------------------
    Readout("nitrogen", "Assimilable nitrogen", "mg/L", ("N",), _scaled("N", 1000.0)),
    Readout(
        "amino_acids",
        "Amino acids",
        "g/L",
        ("amino_acids",),
        _state("amino_acids"),
        media=("wine",),
    ),
    Readout(
        "stored_nitrogen",
        "Stored (intracellular) N",
        "g/L",
        ("stored_nitrogen",),
        _state("stored_nitrogen"),
        media=("wine",),
    ),
    # -- acids and pH ----------------------------------------------------------------
    Readout("ph", "pH", "pH", _CHARGE_ACTIVE, _analysis(analysis.ph_series)),
    Readout(
        "ta",
        "Titratable acidity",
        "g/L tartaric eq.",
        _CHARGE_ACTIVE,
        _analysis(analysis.titratable_acidity_series),
        caveat=_TA_CAVEAT,
        media=("wine",),
    ),
    Readout("malic", "Malic acid", "g/L", ("malic",), _state("malic")),
    Readout("lactic", "Lactic acid", "g/L", ("lactic",), _state("lactic")),
    Readout("tartaric", "Tartaric acid", "g/L", ("tartaric",), _state("tartaric"), media=("wine",)),
    # -- sulfur dioxide and oxygen ---------------------------------------------------
    Readout(
        "so2_free",
        "Free SO2",
        "mg/L",
        ("so2_total", "acetaldehyde"),
        _analysis(analysis.free_so2_series, factor=1000.0),
        media=("wine",),
    ),
    Readout(
        "so2_bound",
        "Bound SO2",
        "mg/L",
        ("so2_total", "acetaldehyde"),
        _analysis(analysis.bound_so2_series, factor=1000.0),
        caveat=_BOUND_SO2_CAVEAT,
        media=("wine",),
    ),
    Readout(
        "so2_molecular",
        "Molecular SO2",
        "mg/L",
        ("so2_total", "acetaldehyde"),
        _analysis(analysis.molecular_so2_series, factor=1000.0),
        media=("wine",),
    ),
    Readout("o2", "Dissolved oxygen", "mg/L", ("o2",), _scaled("o2", 1000.0)),
    Readout("quinone", "Quinones", "mg/L", ("quinone",), _scaled("quinone", 1000.0)),
    Readout(
        "acetaldehyde", "Acetaldehyde", "mg/L", ("acetaldehyde",), _scaled("acetaldehyde", 1000.0)
    ),
    # -- colour and phenolics ---------------------------------------------------------
    Readout("a420", "Browning (A420)", "absorbance", ("A420",), _state("A420")),
    Readout(
        "anthocyanin",
        "Free anthocyanin",
        "g/L",
        ("anthocyanin",),
        _state("anthocyanin"),
        media=("wine",),
    ),
    Readout(
        "polymeric_pigment",
        "Polymeric pigment",
        "g/L",
        ("polymeric_pigment",),
        _state("polymeric_pigment"),
        media=("wine",),
    ),
    Readout("tannin", "Tannin", "g/L", ("tannin",), _state("tannin"), media=("wine",)),
    # -- bitterness (beer) -------------------------------------------------------------
    Readout(
        "ibu",
        "Bitterness",
        "IBU",
        ("iso_alpha",),
        _analysis(analysis.ibu_series, needs_params=False),
        media=("beer",),
    ),
    # -- aroma -------------------------------------------------------------------------
    Readout("diacetyl", "Diacetyl", "mg/L", ("diacetyl",), _scaled("diacetyl", 1000.0)),
    Readout(
        "isoamyl_acetate",
        "Isoamyl acetate (banana)",
        "mg/L",
        ("isoamyl_acetate",),
        _scaled("isoamyl_acetate", 1000.0),
    ),
    Readout(
        "ethyl_hexanoate",
        "Ethyl hexanoate (apple)",
        "mg/L",
        ("ethyl_hexanoate",),
        _scaled("ethyl_hexanoate", 1000.0),
    ),
    Readout(
        "isoamyl_alcohol",
        "Isoamyl alcohol (fusel)",
        "mg/L",
        ("isoamyl_alcohol",),
        _scaled("isoamyl_alcohol", 1000.0),
    ),
    Readout("h2s", "Hydrogen sulfide", "ug/L", ("h2s",), _scaled("h2s", 1_000_000.0)),
    Readout(
        "h2s_gas",
        "H2S lost to headspace",
        "ug/L",
        ("h2s_gas",),
        _scaled("h2s_gas", 1_000_000.0),
        caveat=_GAS_CAVEAT,
    ),
    Readout(
        "ethylphenols",
        "Ethylphenols (Brett)",
        "ug/L",
        ("ethylphenols",),
        _scaled("ethylphenols", 1_000_000.0),
        media=("wine",),
    ),
    Readout(
        "vanillin",
        "Vanillin (oak)",
        "mg/L",
        ("vanillin",),
        _scaled("vanillin", 1000.0),
    ),
    Readout(
        "sotolon",
        "Sotolon (oxidative)",
        "ug/L",
        ("sotolon",),
        _scaled("sotolon", 1_000_000.0),
        media=("wine",),
    ),
)

BY_KEY: Mapping[str, Readout] = {r.key: r for r in READOUTS}


@dataclass(frozen=True)
class Group:
    """A named chart: a title, a shared axis, and the readouts drawn on it."""

    title: str
    keys: tuple[str, ...]
    #: Second y-axis for readouts whose scale would otherwise flatten the first.
    secondary: tuple[str, ...] = field(default_factory=tuple)
    blurb: str = ""


GROUPS: tuple[Group, ...] = (
    Group(
        "Fermentation",
        ("sugar", "ethanol", "biomass"),
        secondary=("biomass",),
        blurb="The three lines that say whether it fermented at all: sugar down, alcohol up, "
        "yeast growing and then settling out.",
    ),
    Group(
        "Sugars (beer)",
        ("glucose", "maltose", "maltotriose"),
        blurb="Yeast eats beer's three sugars in order rather than all at once. If the "
        "maltotriose line never comes down, the model is showing a beer that stops short "
        "of finishing.",
    ),
    Group(
        "Alcohol and gravity",
        ("abv", "brix"),
        blurb="The same run again, in the two numbers someone would actually measure on the tank.",
    ),
    Group(
        "Yeast",
        ("biomass", "biomass_dead", "nitrogen"),
        secondary=("nitrogen",),
        blurb="Yeast stops multiplying when it runs out of nitrogen, which usually happens long "
        "before the sugar is gone. Watch the nitrogen line hit the floor and the yeast line "
        "flatten just after.",
    ),
    Group(
        "Acids and pH",
        ("ph", "malic", "lactic", "tartaric", "ta"),
        secondary=("ph",),
        blurb="pH is not tracked directly. It is worked out from the acids present at each point "
        "in the run, so it only moves because they move.",
    ),
    Group(
        "Sulfur dioxide",
        ("so2_free", "so2_bound", "so2_molecular"),
        secondary=("so2_molecular",),
        blurb="Free sulfur dioxide dips early, then comes back. Nobody wrote that dip into the "
        "model: it happens because a compound the yeast makes early on ties the SO2 up, and "
        "releases it again as that compound is used up.",
    ),
    Group(
        "Oxygen and oxidation",
        ("o2", "quinone", "acetaldehyde", "a420"),
        secondary=("a420",),
        blurb="Oxygen going in at one end, browning coming out at the other.",
    ),
    Group(
        "Colour and phenolics",
        ("a420", "anthocyanin", "polymeric_pigment", "tannin"),
        secondary=("a420",),
        blurb="Loose colour pigments bind into stable ones over time. The wine keeps its colour "
        "and, from then on, cannot easily lose it.",
    ),
    Group(
        "Bitterness",
        ("ibu",),
        blurb="Bitterness drops during fermentation because the bittering compounds stick to the "
        "yeast and get carried out with it. Finished beer is less bitter than the kettle was.",
    ),
    Group(
        "Aroma",
        ("isoamyl_acetate", "ethyl_hexanoate", "isoamyl_alcohol", "diacetyl", "h2s"),
        secondary=("h2s",),
        blurb="The compounds a taster would put a word to - banana, apple, solventy, buttery, "
        "eggy. These are the least well grounded numbers on the whole page. Read whether a "
        "line rises or falls, not how far.",
    ),
    Group(
        "Oak and age",
        ("vanillin", "sotolon", "ethylphenols"),
        blurb="Nothing here moves until the run reaches an aging step. Before that these are all "
        "flat by construction.",
    ),
)


def groups_for(result: RunResult) -> list[tuple[Group, list[Readout]]]:
    """The groups that have something to draw for this run, with their live readouts."""
    medium = result.scenario.medium
    out: list[tuple[Group, list[Readout]]] = []
    for group in GROUPS:
        live = [
            BY_KEY[k]
            for k in group.keys
            if k in BY_KEY and medium in BY_KEY[k].media and BY_KEY[k].available(result)
        ]
        if live:
            out.append((group, live))
    return out


def headline_variables(result: RunResult) -> tuple[str, ...]:
    """State variables the convergence check compares. The ones a verdict rests on."""
    candidates = ("S", "E", "X", "N", "CO2", "Gly", "acetaldehyde")
    return tuple(n for n in candidates if n in result.schema)


def format_value(value: float) -> str:
    """Render a final value. Solver dust around zero is printed as zero, not as 1e-10."""
    if abs(value) < 1e-6:
        return "0"
    return f"{value:,.3g}"


def summary(result: RunResult) -> list[tuple[str, str, str, Tier | None]]:
    """Final-value table: ``(label, value, unit, tier)``. A ``None`` tier means inert."""
    medium = result.scenario.medium
    wanted: Sequence[str] = (
        ("sugar", "abv", "ph", "so2_free", "a420", "biomass")
        if medium == "wine"
        else ("sugar", "abv", "ph", "ibu", "diacetyl", "biomass")
    )
    rows: list[tuple[str, str, str, Tier | None]] = []
    for key in wanted:
        r = BY_KEY.get(key)
        if r is None or medium not in r.media or not r.available(result):
            continue
        value = float(r.series(result)[-1])
        rows.append((r.label, format_value(value), r.unit, r.tier(result)))
    return rows
