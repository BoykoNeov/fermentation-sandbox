"""Starter scenarios, and the metadata the scenario form is generated from.

The form is built from the engine's own tables — the allowed initial-composition keys per
medium and the intervention verb registry — rather than from a copy of them kept here. A
copy would drift the first time a verb is added, and the failure mode of a drifted form is
a traceback from the compile seam in front of the user. What *is* written here is
presentation only: which fields are common enough to show first, and what each one is
called in words rather than in variable names.

The gates at the bottom are the cross-cutting rules the compile seam enforces but that a
form can enforce earlier and more kindly: setting a target pH needs a starting pH to move
from, and sealing a bottle needs a closure named and an aging phase already begun.
"""

from __future__ import annotations

from dataclasses import dataclass

from fermentation.scenario import Intervention, Scenario, TemperaturePoint
from fermentation.scenario.compile import (
    _ALLOWED_KEYS,
    _INTERVENTION_VERBS,
    _SCENARIO_INTERVENTION_VERBS,
)

#: Media the app offers, in the order they appear.
MEDIA: tuple[str, ...] = ("wine", "beer")

#: Every verb the engine knows, read from its own registries so the form cannot drift.
VERBS: tuple[str, ...] = tuple(sorted({*_INTERVENTION_VERBS, *_SCENARIO_INTERVENTION_VERBS}))


def allowed_initial_keys(medium: str) -> tuple[str, ...]:
    """The initial-composition keys the compile seam accepts for this medium."""
    return tuple(sorted(_ALLOWED_KEYS.get(medium, frozenset())))


#: Shown in the main form. Everything else in ``allowed_initial_keys`` is offered under
#: "more inputs", so the common case is six fields and the full vocabulary is still reachable.
PRIMARY_INPUTS: dict[str, tuple[str, ...]] = {
    "wine": ("brix", "yan_mgl", "pitch_gpl", "initial_ph", "tartaric_gpl", "malic_gpl"),
    "beer": (
        "glucose_gpl",
        "maltose_gpl",
        "maltotriose_gpl",
        "yan_mgl",
        "pitch_gpl",
        "initial_ph",
    ),
}

#: Sensible starting values, so an empty form is a runnable scenario rather than an error.
DEFAULT_INPUTS: dict[str, dict[str, float]] = {
    "wine": {
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.5,
        "initial_ph": 3.4,
        "tartaric_gpl": 5.0,
        "malic_gpl": 3.0,
    },
    "beer": {
        "glucose_gpl": 15.0,
        "maltose_gpl": 70.0,
        "maltotriose_gpl": 20.0,
        "yan_mgl": 200.0,
        "pitch_gpl": 1.0,
        "initial_ph": 5.3,
    },
}

#: Plain-language names for the inputs. A key missing here falls back to a tidied-up
#: version of the variable name, so a new engine key still shows up readably.
INPUT_LABELS: dict[str, str] = {
    "brix": "Starting sugar",
    "yan_mgl": "Assimilable nitrogen (YAN)",
    "pitch_gpl": "Yeast pitched",
    "initial_ph": "Starting pH",
    "tartaric_gpl": "Tartaric acid",
    "malic_gpl": "Malic acid",
    "so2_total_mgl": "SO2 at the crusher",
    "mlf_pitch_gpl": "Malolactic culture pitched",
    "citrate_gpl": "Citric acid",
    "copper_gpl": "Copper",
    "glucose_gpl": "Glucose",
    "maltose_gpl": "Maltose",
    "maltotriose_gpl": "Maltotriose",
    "o2_mgl": "Dissolved oxygen at pitch",
    "peptide_buffer_gpl": "Peptide buffer",
    "amino_acids_gpl": "Amino acids",
    "anthocyanin_gpl": "Anthocyanin",
    "tannin_gpl": "Tannin",
    "carrying_capacity_gpl": "Biomass cap",
    "ethanol_gpl": "Ethanol already present",
}


#: One plain sentence per input, shown on the form's own help icon. These describe what the
#: quantity *is* and what a normal batch looks like; they are presentation, not constraint.
#: Nothing here narrows what the engine accepts — a range invented in the interface would be
#: a number with no source behind it, which is precisely what the parameter store forbids.
#: A key missing here simply gets no help icon, which is a visible gap rather than a wrong one.
INPUT_HELP: dict[str, str] = {
    "brix": "How sweet the juice is before fermenting, on the scale a refractometer reads. "
    "Most wine is picked between 20 and 26; higher means more alcohol at the end, and a "
    "harder job for the yeast.",
    "yan_mgl": "The part of the nitrogen in the juice that yeast can actually eat. It is the "
    "usual reason a ferment stalls: comfortable is 200-300, and under about 150 the yeast "
    "runs short before the sugar does.",
    "pitch_gpl": "How much dried yeast goes in per litre. Winemaking is usually 0.25-0.5, "
    "brewing rather more. More yeast starts faster but does not change where it finishes.",
    "initial_ph": "How acidic it starts out. Wine sits near 3.0-3.8, wort near 5.0-5.5. It "
    "affects how the yeast works and how much of the sulfur dioxide is doing anything.",
    "tartaric_gpl": "The main acid of grapes, and the one that decides most of the pH. "
    "Typically 4-8 g/L in juice.",
    "malic_gpl": "The sharp apple acid. Cool years leave more of it; a malolactic culture "
    "later converts it to the softer lactic acid.",
    "so2_total_mgl": "Sulfur dioxide added at the crusher, the standard antioxidant and "
    "antimicrobial. Commonly 30-60 mg/L. Only the small free fraction does the work.",
    "glucose_gpl": "The simple sugar in wort. Yeast eats this first, before the maltose.",
    "maltose_gpl": "The main sugar in wort, and most of the eventual alcohol. Usually the "
    "largest of the three by far.",
    "maltotriose_gpl": "The slowest of the three wort sugars, taken up last and often left "
    "partly behind — which is where a lot of the sweetness in a finished beer comes from.",
    "o2_mgl": "Oxygen dissolved into the wort when the yeast goes in. Brewers aerate on "
    "purpose, to around 8-10 mg/L; the yeast needs it to build cell membranes.",
    "mlf_pitch_gpl": "Malolactic bacteria, which convert the sharp malic acid to softer "
    "lactic acid. A very small dose does it.",
    "citrate_gpl": "Citric acid — a minor acid in grapes, sometimes added.",
    "copper_gpl": "Copper already in the juice. It is the fining for the rotten-egg smell, "
    "and also the metal that drives browning, so it cuts both ways.",
    "peptide_buffer_gpl": "Short protein fragments in wort that resist pH change. They are "
    "why a beer's pH falls less than the acid made would suggest.",
    "amino_acids_gpl": "The nitrogen the yeast eats, counted as the amino acids themselves "
    "rather than as one lump. Also the raw material for much of the aroma.",
    "anthocyanin_gpl": "Red grape pigment. A white must has none, which is why the colour "
    "chart is empty for one.",
    "tannin_gpl": "The compounds from skins, seeds and oak that taste astringent and that "
    "bind pigment as the wine ages.",
    "carrying_capacity_gpl": "A ceiling on how much yeast the batch can hold, whatever else "
    "is available. Left at zero, the model works it out instead.",
    "ethanol_gpl": "Alcohol already present at the start — for a batch being restarted, or "
    "one that has been fortified.",
}


def input_label(key: str) -> str:
    return INPUT_LABELS.get(key, key.replace("_gpl", "").replace("_mgl", "").replace("_", " "))


def input_help(key: str) -> str | None:
    """The help sentence for an input, or ``None`` where there is nothing to say yet."""
    return INPUT_HELP.get(key)


def input_unit(key: str) -> str:
    """Unit inferred from the key's own suffix — the engine's naming convention is the unit."""
    if key == "brix":
        return "degrees Brix"
    if key.endswith("_ph") or key == "initial_ph":
        return "pH"
    if key.endswith("_per_h"):
        return "1/h"
    if key.endswith("_ugl"):
        return "ug/L"
    if key.endswith("_mgl"):
        return "mg/L"
    if key.endswith("_gpl"):
        return "g/L"
    return ""


@dataclass(frozen=True)
class VerbSpec:
    """One intervention verb as the form needs to present it."""

    verb: str
    label: str
    #: ``(param name, label, unit, default)`` for each numeric argument.
    numbers: tuple[tuple[str, str, str, float], ...] = ()
    #: ``(param name, label, choices)`` for each categorical argument.
    choices: tuple[tuple[str, str, tuple[str, ...]], ...] = ()
    media: tuple[str, ...] = ("wine", "beer")
    note: str = ""


#: Presentation for each verb. The set of verbs itself comes from :data:`VERBS`; a verb the
#: engine grows that is missing here still appears, with its raw name and no arguments, which
#: is a visible gap rather than a silent omission.
VERB_SPECS: dict[str, VerbSpec] = {
    "add_dap": VerbSpec(
        "add_dap",
        "Add DAP (nutrient)",
        numbers=(("dap_gpl", "Diammonium phosphate", "g/L", 0.25),),
        note="Doses both ions — the ammonium the yeast eats and the phosphate it does not.",
    ),
    "add_so2": VerbSpec(
        "add_so2", "Add SO2", numbers=(("so2_mgl", "Sulfur dioxide", "mg/L", 50.0),)
    ),
    "add_copper": VerbSpec(
        "add_copper",
        "Add copper",
        numbers=(("copper_mgl", "Copper", "mg/L", 0.3),),
        note="The fining for hydrogen sulfide.",
    ),
    "add_acid": VerbSpec(
        "add_acid",
        "Acidify",
        numbers=(("gpl", "Amount", "g/L", 1.0),),
        choices=(("acid", "Acid", ("tartaric", "malic", "lactic", "citric")),),
    ),
    "set_ph": VerbSpec(
        "set_ph",
        "Set pH to a target",
        numbers=(("ph", "Target pH", "pH", 3.4),),
        note="Needs a starting pH on the scenario — it moves the anchor, it does not invent one.",
    ),
    "add_sugar": VerbSpec("add_sugar", "Add sugar", numbers=(("sugar_gpl", "Sugar", "g/L", 20.0),)),
    "add_oxygen": VerbSpec(
        "add_oxygen",
        "Expose to oxygen",
        numbers=(("o2_mgl", "Oxygen", "mg/L", 2.0),),
        note="A racking, a micro-oxygenation, a bottling line's pickup.",
    ),
    "add_ascorbate": VerbSpec(
        "add_ascorbate", "Add ascorbate", numbers=(("ascorbate_mgl", "Ascorbate", "mg/L", 50.0),)
    ),
    "add_oak": VerbSpec(
        "add_oak",
        "Add oak",
        numbers=(
            ("oak_gpl", "Oak", "g/L", 3.0),
            ("fill_number", "Barrel fill number", "", 1.0),
        ),
        choices=(
            ("toast", "Toast level", ("light", "medium", "heavy")),
            ("spirit", "Previous contents", ("none", "bourbon", "sherry")),
        ),
        media=("wine",),
    ),
    "rack": VerbSpec(
        "rack",
        "Rack off the lees",
        numbers=(("fraction", "Fraction of lees removed", "0-1", 0.9),),
    ),
    "pitch_mlf": VerbSpec(
        "pitch_mlf",
        "Pitch malolactic culture",
        numbers=(("pitch_gpl", "Culture", "g/L", 0.01),),
        media=("wine",),
    ),
    "pitch_brett": VerbSpec(
        "pitch_brett",
        "Pitch Brettanomyces",
        numbers=(("pitch_gpl", "Culture", "g/L", 0.001),),
        media=("wine",),
        note="A spoilage organism. Modelled because it happens, not because it is wanted.",
    ),
    "begin_aging": VerbSpec(
        "begin_aging",
        "Begin aging",
        note="Switches on the slow chemistry: oxidation, pigment polymerisation, ester "
        "hydrolysis. Nothing in the oak or age charts moves before this.",
    ),
    "seal_bottle": VerbSpec(
        "seal_bottle",
        "Seal the bottle",
        note="Doses the one-off oxygen burst trapped at bottling. Needs a closure chosen, and "
        "must not come before aging begins.",
        media=("wine",),
    ),
}

#: Closures the compile seam knows, each carrying its own sourced oxygen transmission rate.
CLOSURES: tuple[str, ...] = (
    "hermetic",
    "technical_cork",
    "screwcap",
    "natural_cork",
    "synthetic_nomacorc",
    "synthetic_supremecorq",
)


def gate_problems(scenario: Scenario) -> list[str]:
    """Cross-cutting rules a form can catch before the compile seam raises.

    These are the constraints that span two parts of a scenario, so no single field can
    validate them. Returning them as sentences lets the form disable the offending control
    and say why, which is a better answer than a traceback.
    """
    problems: list[str] = []
    actions = [iv.action for iv in scenario.interventions]

    if "set_ph" in actions and "initial_ph" not in scenario.initial:
        problems.append(
            "'Set pH to a target' needs a starting pH in the composition — the verb moves the "
            "acid-base anchor, and there is nothing to move without one."
        )
    if "seal_bottle" in actions:
        if scenario.closure is None:
            problems.append(
                "'Seal the bottle' needs a closure chosen. The oxygen burst it doses is looked "
                "up from the closure by name, so an unnamed closure has no dose."
            )
        seal_days = [iv.day for iv in scenario.interventions if iv.action == "seal_bottle"]
        aging_days = [iv.day for iv in scenario.interventions if iv.action == "begin_aging"]
        if not aging_days:
            problems.append(
                "'Seal the bottle' comes after 'Begin aging'. Add the aging step first."
            )
        elif min(seal_days) < min(aging_days):
            problems.append(
                f"'Seal the bottle' is scheduled on day {min(seal_days):g}, before aging begins "
                f"on day {min(aging_days):g}. A bottle cannot be sealed before it is filled."
            )
    if scenario.hops and scenario.batch_volume_liters is None:
        problems.append(
            "Hops are dosed by mass, so the batch needs a volume before grams can become a "
            "concentration."
        )
    for iv in scenario.interventions:
        if iv.day > scenario.duration_days:
            problems.append(
                f"'{VERB_SPECS.get(iv.action, VerbSpec(iv.action, iv.action)).label}' is "
                f"scheduled on day {iv.day:g}, after the run ends on day "
                f"{scenario.duration_days:g}. It will never fire."
            )
    return problems


# -- starter scenarios ---------------------------------------------------------------------


def _temps(*points: tuple[float, float]) -> list[TemperaturePoint]:
    return [TemperaturePoint(day=d, celsius=c) for d, c in points]


STARTERS: dict[str, Scenario] = {
    "White wine, cool and clean": Scenario(
        name="White wine, cool and clean",
        medium="wine",
        initial={
            "brix": 22.0,
            "yan_mgl": 220.0,
            "pitch_gpl": 0.5,
            "initial_ph": 3.3,
            "tartaric_gpl": 6.0,
            "malic_gpl": 3.5,
            "so2_total_mgl": 40.0,
        },
        temperature_schedule=_temps((0.0, 15.0), (14.0, 15.0)),
        duration_days=16.0,
    ),
    "Red wine, warm ferment then malolactic": Scenario(
        name="Red wine, warm ferment then malolactic",
        medium="wine",
        initial={
            "brix": 25.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.5,
            "initial_ph": 3.6,
            "tartaric_gpl": 5.5,
            "malic_gpl": 2.5,
            "anthocyanin_gpl": 0.5,
            "tannin_gpl": 1.5,
        },
        temperature_schedule=_temps((0.0, 24.0), (10.0, 26.0), (30.0, 18.0)),
        interventions=[
            Intervention(day=2.0, action="add_dap", params={"dap_gpl": 0.25}),
            Intervention(day=10.0, action="rack", params={"fraction": 0.9}),
            Intervention(day=11.0, action="pitch_mlf", params={"pitch_gpl": 0.01}),
        ],
        duration_days=40.0,
    ),
    "Nitrogen-starved must (a stuck ferment)": Scenario(
        name="Nitrogen-starved must (a stuck ferment)",
        medium="wine",
        initial={
            "brix": 26.0,
            "yan_mgl": 90.0,
            "pitch_gpl": 0.4,
            "initial_ph": 3.6,
            "tartaric_gpl": 4.0,
            "malic_gpl": 1.5,
        },
        temperature_schedule=_temps((0.0, 20.0), (21.0, 20.0)),
        duration_days=21.0,
    ),
    "Bottle-aged white, screwcap": Scenario(
        name="Bottle-aged white, screwcap",
        medium="wine",
        initial={
            "brix": 21.0,
            "yan_mgl": 230.0,
            "pitch_gpl": 0.5,
            "initial_ph": 3.2,
            "tartaric_gpl": 6.5,
            "malic_gpl": 3.0,
            "so2_total_mgl": 60.0,
        },
        temperature_schedule=_temps((0.0, 16.0), (14.0, 16.0), (15.0, 14.0), (365.0, 14.0)),
        interventions=[
            Intervention(day=14.0, action="begin_aging"),
            Intervention(day=14.0, action="seal_bottle"),
        ],
        closure="screwcap",
        duration_days=365.0,
    ),
    "Pale ale": Scenario(
        name="Pale ale",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 1.0,
            "initial_ph": 5.3,
            "o2_mgl": 8.0,
        },
        temperature_schedule=_temps((0.0, 19.0), (7.0, 20.0)),
        duration_days=9.0,
    ),
    "Lager, cold and slow": Scenario(
        name="Lager, cold and slow",
        medium="beer",
        initial={
            "glucose_gpl": 12.0,
            "maltose_gpl": 65.0,
            "maltotriose_gpl": 18.0,
            "yan_mgl": 180.0,
            "pitch_gpl": 1.5,
            "initial_ph": 5.4,
            "o2_mgl": 9.0,
        },
        temperature_schedule=_temps((0.0, 10.0), (10.0, 12.0), (14.0, 18.0), (21.0, 18.0)),
        duration_days=21.0,
    ),
}
