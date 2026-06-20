"""The scenario → core compile seam.

A :class:`~fermentation.scenario.schema.Scenario` is declarative and expressed in
industry units (degrees Brix, mg/L of nitrogen, degrees C, days).
:func:`compile_scenario` turns it into everything the runtime needs to integrate:

    * ``y0``         — the initial state vector in canonical units (g/L, K),
    * ``process_set`` — the medium's Processes, assembled against its schema,
    * ``parameters``  — the provenance-backed parameter set for that medium/strain.

This is the *only* place industry units cross into the canonical internal
representation (decision D-3); the core never sees a degree Brix. Physics does not
live here — it stays in the core's Processes — so this module is pure plumbing:
look up the medium, convert the initial composition, load the parameters, and
assemble the Process set.

The accepted ``Scenario.initial`` keys are validated here (the schema deliberately
leaves them as a free ``dict`` so the vocabulary can live at this boundary):

    wine: brix, yan_mgl, pitch_gpl, [ethanol_gpl]
    beer: glucose_gpl, maltose_gpl, maltotriose_gpl, yan_mgl, pitch_gpl, [ethanol_gpl]

Beer's three sugars are given explicitly rather than split from a single original
gravity: that wort spectrum is a provenance-backed parameter (Milestone 1's
sourcing task), not a magic constant to bury in the compile step.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from fermentation.core.media import get_medium
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.parameters.store import ParameterSet, default_data_dir, load_parameters
from fermentation.scenario.schema import Scenario
from fermentation.units.convert import (
    brix_to_sugar_gpl,
    celsius_to_kelvin,
    days_to_hours,
    mgl_to_gpl,
)

#: A name → value(s) mapping ready for :meth:`StateSchema.pack`.
_Initial = dict[str, float | list[float]]


@dataclass(frozen=True, eq=False)
class CompiledScenario:
    """Everything the runtime needs to integrate one scenario.

    Realises the documented compile seam ``(y0, ProcessSet, params)`` as a named
    record, plus the schema and time span that travel with them. ``param_values``
    is the plain ``{name: float}`` mapping ``simulate`` and ``Process.derivatives``
    consume; ``parameters`` retains the full provenance and tier information for
    honest downstream reporting.
    """

    scenario: Scenario
    schema: StateSchema
    y0: FloatArray
    process_set: ProcessSet
    parameters: ParameterSet
    t_span_h: tuple[float, float]

    @property
    def param_values(self) -> dict[str, float]:
        """Resolved ``{name: value}`` mapping for the integration hot loop."""
        return self.parameters.resolve()


# -- initial-composition vocabulary (the industry-unit boundary) --------------

#: Keys accepted in ``Scenario.initial`` per medium. Validated at compile time so
#: a typo ("brixx") fails loudly instead of being silently ignored.
_ALLOWED_KEYS: dict[str, frozenset[str]] = {
    "wine": frozenset({"brix", "yan_mgl", "pitch_gpl", "ethanol_gpl"}),
    "beer": frozenset(
        {"glucose_gpl", "maltose_gpl", "maltotriose_gpl", "yan_mgl", "pitch_gpl", "ethanol_gpl"}
    ),
}


def _nonneg(value: float, key: str) -> float:
    if value < 0.0:
        raise ValueError(f"scenario.initial[{key!r}] must be >= 0, got {value}")
    return value


def _require(values: Mapping[str, float], key: str, medium: str) -> float:
    if key not in values:
        raise ValueError(f"{medium} scenario.initial is missing required key {key!r}")
    return _nonneg(float(values[key]), key)


def _optional(values: Mapping[str, float], key: str, default: float) -> float:
    return _nonneg(float(values[key]), key) if key in values else default


def _wine_initial(values: Mapping[str, float], temperature_k: float) -> _Initial:
    return {
        "X": _require(values, "pitch_gpl", "wine"),
        "S": [brix_to_sugar_gpl(_require(values, "brix", "wine"))],
        "E": _optional(values, "ethanol_gpl", 0.0),
        "N": mgl_to_gpl(_require(values, "yan_mgl", "wine")),
        "T": temperature_k,
        "CO2": 0.0,
    }


def _beer_initial(values: Mapping[str, float], temperature_k: float) -> _Initial:
    return {
        "X": _require(values, "pitch_gpl", "beer"),
        "S": [
            _require(values, "glucose_gpl", "beer"),
            _require(values, "maltose_gpl", "beer"),
            _require(values, "maltotriose_gpl", "beer"),
        ],
        "E": _optional(values, "ethanol_gpl", 0.0),
        "N": mgl_to_gpl(_require(values, "yan_mgl", "beer")),
        "T": temperature_k,
        "CO2": 0.0,
    }


_INITIAL_BUILDERS: dict[str, Callable[[Mapping[str, float], float], _Initial]] = {
    "wine": _wine_initial,
    "beer": _beer_initial,
}


def _validate_initial_keys(scenario: Scenario) -> None:
    allowed = _ALLOWED_KEYS.get(scenario.medium)
    if allowed is None:
        raise ValueError(
            f"medium {scenario.medium!r} has no initial-composition vocabulary defined"
        )
    unknown = set(scenario.initial) - allowed
    if unknown:
        raise ValueError(
            f"scenario.initial has unknown key(s) {sorted(unknown)} for medium "
            f"{scenario.medium!r}; allowed: {sorted(allowed)}"
        )


def _initial_temperature_kelvin(scenario: Scenario) -> float:
    schedule = scenario.temperature_schedule
    if not schedule:
        raise ValueError(
            f"scenario {scenario.name!r}: temperature_schedule needs at least one point "
            "to seed the initial temperature"
        )
    earliest = min(schedule, key=lambda point: point.day)
    return celsius_to_kelvin(earliest.celsius)


def _load_parameters(
    scenario: Scenario,
    parameter_paths: Sequence[str | Path] | None,
    data_dir: str | Path | None,
) -> ParameterSet:
    if parameter_paths is not None:
        return load_parameters(*parameter_paths)
    base = Path(data_dir) if data_dir is not None else default_data_dir()
    path = base / f"{scenario.medium}_{scenario.strain}.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"no parameter file for medium={scenario.medium!r} strain={scenario.strain!r}: "
            f"expected {path}. Pass parameter_paths=... or add the YAML "
            "(see the Milestone 1 parameter-sourcing task)."
        )
    return load_parameters(path)


def compile_scenario(
    scenario: Scenario,
    *,
    parameter_paths: Sequence[str | Path] | None = None,
    data_dir: str | Path | None = None,
    strict: bool = False,
) -> CompiledScenario:
    """Compile a declarative scenario into an integrable :class:`CompiledScenario`.

    Industry units in ``scenario.initial`` are converted to canonical units here
    and nowhere else. ``parameter_paths`` overrides the default lookup of
    ``<medium>_<strain>.yaml`` under ``data_dir`` (or the packaged data dir);
    ``strict=True`` enables the Process ``touches`` contract on the returned set.

    Raises ``KeyError`` for an unknown medium, ``ValueError`` for an invalid
    initial composition or missing temperature, and ``FileNotFoundError`` when the
    medium/strain has no parameter file yet.
    """
    medium = get_medium(scenario.medium)
    _validate_initial_keys(scenario)

    builder = _INITIAL_BUILDERS.get(scenario.medium)
    if builder is None:
        raise ValueError(f"no initial-composition builder for medium {scenario.medium!r}")

    temperature_k = _initial_temperature_kelvin(scenario)
    y0 = medium.schema.pack(builder(scenario.initial, temperature_k))

    parameters = _load_parameters(scenario, parameter_paths, data_dir)
    process_set = medium.build_process_set(strict=strict)
    t_span_h = (0.0, days_to_hours(scenario.duration_days))

    return CompiledScenario(
        scenario=scenario,
        schema=medium.schema,
        y0=y0,
        process_set=process_set,
        parameters=parameters,
        t_span_h=t_span_h,
    )
