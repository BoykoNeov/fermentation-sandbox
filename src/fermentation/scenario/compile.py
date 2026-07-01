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

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from fermentation.core import acidbase
from fermentation.core.kinetics import MalolacticConversion
from fermentation.core.media import get_medium
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import combine
from fermentation.parameters.schema import Parameter, Provenance, Uncertainty
from fermentation.parameters.store import ParameterSet, default_data_dir, load_parameters
from fermentation.scenario.schema import Scenario
from fermentation.units.convert import (
    brix_to_sugar_gpl,
    celsius_to_kelvin,
    days_to_hours,
    mgl_to_gpl,
)

#: Coleman Y_X/N regression coefficients (decision D-14). Present iff a medium
#: ships the nitrogen-dependent biomass yield; gates the compile-time override.
_N_YIELD_COEFFS = ("biomass_N_yield_log_intercept", "biomass_N_yield_log_slope")

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
    "wine": frozenset(
        # tartaric_gpl/malic_gpl/initial_ph are the optional pH-solver inputs (D-18);
        # lactic is produced-only (MLF product) so it is not an input, and the
        # strong cation is back-solved from initial_ph, not given. so2_total_mgl is the
        # optional total-SO₂ dose for the free/bound + molecular-SO₂ readout (D-22/D-28);
        # mlf_pitch_gpl is the optional Oenococcus oeni dose driving malolactic conversion (D-23).
        {
            "brix",
            "yan_mgl",
            "pitch_gpl",
            "ethanol_gpl",
            "tartaric_gpl",
            "malic_gpl",
            "initial_ph",
            "so2_total_mgl",
            "mlf_pitch_gpl",
        }
    ),
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


def _wine_initial(
    values: Mapping[str, float], temperature_k: float, parameters: ParameterSet
) -> _Initial:
    # Brix measures *total* dissolved solids; only ~90-95% of ripe-must solids are
    # fermentable hexose (the rest is acids/minerals/phenolics). The sourced
    # must_fermentable_fraction corrects brix_to_sugar_gpl so a 24 Brix must loads
    # realistic fermentable sugar (~245 g/L, not 264) and the wine ABV is realistic
    # (decision D-16). Absent ⇒ 1.0 (no correction), so older parameter sets still
    # compile. Produced-only pools (X_dead, Gly, Byp, esters, fusels) default to 0
    # (see VarSpec) and so start empty at pitch.
    fermentable_fraction = (
        parameters["must_fermentable_fraction"].value
        if "must_fermentable_fraction" in parameters
        else 1.0
    )
    sugar_gpl = brix_to_sugar_gpl(_require(values, "brix", "wine")) * fermentable_fraction
    # pH-solver acid inputs (decision D-18), all optional so acid-free scenarios still
    # compile (slots default to 0, inert). tartaric/malic are must inputs in g/L;
    # lactic is produced-only (MLF product), 0 at pitch. The net strong cation is
    # back-solved from the measured initial_ph so the modelled pH reproduces it at t=0
    # (inverse anchoring): D-18 predicts pH *changes*, not absolute initial pH.
    tartaric = _optional(values, "tartaric_gpl", 0.0)
    malic = _optional(values, "malic_gpl", 0.0)
    initial: _Initial = {
        "X": _require(values, "pitch_gpl", "wine"),
        "S": [sugar_gpl],
        "E": _optional(values, "ethanol_gpl", 0.0),
        "N": mgl_to_gpl(_require(values, "yan_mgl", "wine")),
        "T": temperature_k,
        "CO2": 0.0,
        "X_dead": 0.0,  # no inactivated biomass at pitch
        "Gly": 0.0,  # no byproducts at pitch (decision D-16)
        "Byp": 0.0,
        "esters": 0.0,  # produced-only aroma pools, empty at pitch (decision D-19)
        "fusels": 0.0,
        "esters_gas": 0.0,  # volatilized-ester bookkeeping pool, empty at pitch (D-20)
        "tartaric": tartaric,
        "malic": malic,
        "lactic": 0.0,
        "cation_charge": 0.0,  # back-solved below iff initial_ph is given
        # Total-SO₂ dose for the free/bound + molecular-SO₂ readout (D-22/D-28); mg/L→g/L,
        # default 0 (no dose). Inert/conserved state (readout-only, not in the charge
        # balance), so it does NOT enter the cation back-solve below — SO₂'s minor bisulfite
        # charge is a scoped omission the inverse anchoring would absorb at t=0 anyway (D-22).
        # Free/bound are derived from this total + acetaldehyde at the solved pH (D-28).
        "so2_total": mgl_to_gpl(_optional(values, "so2_total_mgl", 0.0)),
        # Oenococcus oeni dose driving malolactic conversion (D-23); g/L, default 0 (no
        # MLF). Inert catalyst in v1 (no Process grows/kills it) and carbon-free, so an
        # undosed run is byte-for-byte the validated core; the compile step below disables
        # the MLF Process entirely when this is 0 (tier + perf isolability).
        "X_mlf": _optional(values, "mlf_pitch_gpl", 0.0),
    }
    if "initial_ph" in values:
        # Byp = 0 at pitch, so the anchoring cation reproduces initial_ph from the named
        # acids alone; as Byp accumulates during the ferment, pH drifts emergently.
        acid_gpl = {"tartaric": tartaric, "malic": malic, "lactic": 0.0}
        totals_molar = {n: g / acidbase.ACID_STATE[n].molar_mass for n, g in acid_gpl.items()}
        try:
            initial["cation_charge"] = acidbase.solve_cation_charge(
                totals_molar,
                byp_succinic_molar=0.0,
                pka_map=acidbase.build_pka_map(parameters.resolve()),
                target_ph=float(values["initial_ph"]),
            )
        except ValueError as exc:  # initial_ph below the acid load's intrinsic pH
            raise ValueError(f"wine scenario.initial['initial_ph'] is unphysical: {exc}") from exc
        except KeyError as exc:  # acidbase.yaml pKa parameters not loaded
            raise ValueError(
                "wine scenario gives 'initial_ph' but the pKa parameters are missing "
                f"({exc}); include acidbase.yaml in parameter_paths (the default lookup "
                "merges it automatically)."
            ) from exc
    return initial


def _beer_initial(
    values: Mapping[str, float], temperature_k: float, parameters: ParameterSet
) -> _Initial:
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
        "X_dead": 0.0,  # no inactivated biomass at pitch
        "Gly": 0.0,  # beer carries zero byproduct diversion in M1 (decision D-16)
        "Byp": 0.0,
        "esters": 0.0,  # produced-only aroma pools, empty at pitch (decision D-19)
        "fusels": 0.0,
        "esters_gas": 0.0,  # volatilized-ester bookkeeping pool, empty at pitch (D-20)
    }


_INITIAL_BUILDERS: dict[str, Callable[[Mapping[str, float], float, ParameterSet], _Initial]] = {
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
        # Caller-controlled override: a caller wanting the pH solver must include
        # acidbase.yaml in their paths (the pKa set the charge balance reads, D-18).
        return load_parameters(*parameter_paths)
    base = Path(data_dir) if data_dir is not None else default_data_dir()
    path = base / f"{scenario.medium}_{scenario.strain}.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"no parameter file for medium={scenario.medium!r} strain={scenario.strain!r}: "
            f"expected {path}. Pass parameter_paths=... or add the YAML "
            "(see the Milestone 1 parameter-sourcing task)."
        )
    # Merge the shared, medium-agnostic parameter files alongside the medium file so every
    # default-lookup scenario can compute pH (acidbase.yaml, decision D-18), run the diacetyl
    # pathway (vicinal_diketones.yaml, decision D-26 — the load-bearing decarb step is
    # non-enzymatic, so its constants are medium-agnostic), the acetaldehyde buffer
    # (acetaldehyde.yaml, decision D-27 — main-pathway yeast metabolism, likewise generic) and
    # H₂S production (hydrogen_sulfide.yaml, decision D-29 — the sulfate-reduction sequence,
    # generic yeast metabolism). The names are collision-free with the per-medium kinetic
    # parameters; load_parameters merges left-to-right.
    shared_files = [
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "hydrogen_sulfide.yaml",
    ]
    return load_parameters(path, *(f for f in shared_files if f.exists()))


def _apply_nitrogen_dependent_yield(scenario: Scenario, parameters: ParameterSet) -> ParameterSet:
    """Override ``biomass_N_fraction`` from Coleman's ``Y_X/N(N_init)`` regression.

    Coleman, Fish & Block (2007) found the cell-mass-per-nitrogen yield to depend
    on the *initial* nitrogen (Fig 4 / Table A2): ``ln(Y_X/N) = a0 + a1·YAN``
    (YAN in mg N/L). This is the one parameter that cannot be pre-evaluated into
    the YAML the way the temperature regressions are — the evaluation point is the
    scenario's nitrogen, not a fixed reference — so it is computed here at the
    compile boundary and nowhere else (decision D-14). Because every assimilated
    gram of nitrogen enters biomass in our model, ``Y_X/N = 1/f_N`` identically;
    setting ``biomass_N_fraction = 1/Y_X/N`` leaves the nitrogen balance exact (the
    ``total_nitrogen`` check reads this same per-run constant).

    Gated on the regression coefficients being present, so a medium without them
    (beer) keeps the static elemental ``biomass_N_fraction`` untouched.
    """
    if not all(name in parameters for name in _N_YIELD_COEFFS):
        return parameters
    yan_mgl = scenario.initial.get("yan_mgl")
    if yan_mgl is None:
        return parameters

    a0, a1 = (parameters[name] for name in _N_YIELD_COEFFS)
    y_xn = math.exp(a0.value + a1.value * float(yan_mgl))  # g cell / g N
    f_n = 1.0 / y_xn
    override = Parameter(
        name="biomass_N_fraction",
        value=f_n,
        unit="g/g",
        tier=combine((a0.tier, a1.tier)),
        uncertainty=Uncertainty(
            # Bracketing metadata, not a tuned value: f_N = 1/Y_X/N ranges
            # ~0.039-0.107 across Coleman's 70-350 mg N/L treatment span
            # (Y_X/N ~25.7 down to ~9.4); [0.03, 0.15] brackets that with margin.
            low=0.03,
            high=0.15,
            note="nitrogen-status-dependent; brackets f_N across Coleman's 70-350 mg N/L range",
        ),
        provenance=Provenance(
            source=a0.provenance.source,
            doi=a0.provenance.doi,
            conditions=(
                f"computed at compile from Coleman Y_X/N regression at YAN={float(yan_mgl):g} mg/L"
            ),
            notes=(
                f"Y_X/N = exp({a0.value} + {a1.value}*{float(yan_mgl):g}) = {y_xn:.2f} g cell/g N; "
                f"f_N = 1/Y_X/N = {f_n:.4f} g N/g cell. Overrides the static elemental "
                "biomass_N_fraction so a nitrogen-limited must builds realistically little "
                "biomass (decision D-14)."
            ),
        ),
    )
    return parameters.merge(ParameterSet([override]), override=True)


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

    # Parameters are loaded *before* y0 because the wine initial sugar applies a
    # sourced must_fermentable_fraction (decision D-16), mirroring how the
    # nitrogen-dependent yield (D-14) is also resolved at this boundary.
    parameters = _load_parameters(scenario, parameter_paths, data_dir)
    parameters = _apply_nitrogen_dependent_yield(scenario, parameters)

    y0 = medium.schema.pack(builder(scenario.initial, temperature_k, parameters))
    process_set = medium.build_process_set(strict=strict)

    # MLF isolability (decision D-23): the malolactic Process is wired into the wine medium
    # but contributes nothing until Oenococcus oeni is pitched. When it is not, DISABLE it
    # so (a) the inert ``malic``/``lactic`` slots keep their VALIDATED tier — an *enabled*
    # Process that touches them drops them to speculative even with a zero contribution,
    # since ``tier_of`` counts enabled, not nonzero, Processes — and (b) no per-RHS pH
    # ``brentq`` solve is paid on an undosed run. When pitched it is the first RHS consumer
    # of the D-18 pH solver and the D-22 molecular-SO₂ readout.
    if MalolacticConversion.name in process_set:
        mlf_pitch_gpl = float(scenario.initial.get("mlf_pitch_gpl", 0.0) or 0.0)
        if mlf_pitch_gpl <= 0.0:
            process_set.disable(MalolacticConversion.name)

    t_span_h = (0.0, days_to_hours(scenario.duration_days))

    return CompiledScenario(
        scenario=scenario,
        schema=medium.schema,
        y0=y0,
        process_set=process_set,
        parameters=parameters,
        t_span_h=t_span_h,
    )
