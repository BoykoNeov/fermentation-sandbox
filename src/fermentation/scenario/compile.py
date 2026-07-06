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
from dataclasses import dataclass, field
from pathlib import Path

from fermentation.core import acidbase
from fermentation.core.kinetics import (
    AminoAcidAssimilation,
    AutolyticHydrogenSulfide,
    AutolyticMercaptan,
    BiomassCarryingCapacity,
    BrettDeath,
    BrettDecarboxylation,
    BrettGrowth,
    BrettVinylphenolReduction,
    FuselAminoAcidReroute,
    MalolacticCitrateMetabolism,
    MalolacticConversion,
    MalolacticDeath,
    MalolacticGrowth,
    MalolacticSenescence,
    OenococcusDiacetylReduction,
    YeastAutolysis,
    YeastPOFDecarboxylation,
)
from fermentation.core.kinetics.temperature import RAMP_RATE
from fermentation.core.media import get_medium
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier, combine
from fermentation.parameters.schema import Parameter, Provenance, Uncertainty
from fermentation.parameters.store import ParameterSet, default_data_dir, load_parameters
from fermentation.runtime.ensemble import Ensemble, simulate_ensemble
from fermentation.runtime.schedule import ScheduledEvent, ScheduledTrajectory, simulate_scheduled
from fermentation.scenario.schema import Intervention, Scenario
from fermentation.units.convert import (
    brix_to_sugar_gpl,
    celsius_to_kelvin,
    days_to_hours,
    mgl_to_gpl,
)

#: Coleman Y_X/N regression coefficients (decision D-14). Present iff a medium
#: ships the nitrogen-dependent biomass yield; gates the compile-time override.
_N_YIELD_COEFFS = ("biomass_N_yield_log_intercept", "biomass_N_yield_log_slope")

#: The malolactic Processes gated on an *Oenococcus oeni* pitch (decisions D-23, D-31, D-39):
#: malate→lactate conversion, the citrate co-metabolism feeding the diacetyl reservoir, the
#: bacterial diacetyl reduction, and bacterial death/decay. They are wired into the wine medium
#: but contribute nothing until bacteria are present, so the compile step DISABLES them when
#: unpitched and the ``pitch_mlf`` intervention (decision D-36) re-enables *exactly* this set at
#: its breakpoint — a single source of truth so the compile-time gate and the mid-run pitch cannot
#: drift apart. :class:`MalolacticDeath` (D-39) and :class:`MalolacticSenescence` (MLF v2, D-41) are
#: pitch-gated too (bacteria die/age whether or not amino acids were dosed), unlike
#: :class:`MalolacticGrowth`, which is amino-acid-gated below.
_MLF_GATED_PROCESSES = (
    MalolacticConversion,
    MalolacticCitrateMetabolism,
    OenococcusDiacetylReduction,
    MalolacticDeath,
    MalolacticSenescence,
)

#: The *Brettanomyces* Processes gated on a Brett pitch (decision D-40): hydroxycinnamate
#: decarboxylation, vinylphenol reduction, and (D-40 pt3) the SO₂-driven :class:`BrettDeath`. Wired
#: into the wine medium but contributing nothing until Brett is present, so compile DISABLES them
#: (unpitched) and the ``pitch_brett`` intervention re-enables *exactly* this set at its
#: breakpoint — one source of truth so the compile-time gate and the mid-run pitch cannot drift
#: apart (the ``_MLF_GATED_PROCESSES`` pattern). :class:`BrettDeath` is pitch-gated too (Brett dies
#: whether or not amino acids were dosed), unlike :class:`BrettGrowth`, which is amino-acid-gated
#: below (the exact :class:`MalolacticDeath` vs :class:`MalolacticGrowth` split).
_BRETT_GATED_PROCESSES = (
    BrettDecarboxylation,
    BrettVinylphenolReduction,
    BrettDeath,
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
    #: Timed interventions compiled from the scenario, in canonical hours, ready to hand to
    #: :func:`fermentation.runtime.simulate_scheduled`: the temperature-schedule slope-change
    #: events (decision D-35) merged with the discrete winemaking verbs — ``add_dap`` / ``add_so2``
    #: / ``rack`` / ``pitch_mlf`` (decision D-36). Empty ⇒ an un-scheduled run (plain
    #: :func:`simulate` suffices).
    events: tuple[ScheduledEvent, ...] = field(default_factory=tuple)

    @property
    def param_values(self) -> dict[str, float]:
        """Resolved ``{name: value}`` mapping for the integration hot loop."""
        return self.parameters.resolve()

    def run(self, **kwargs: object) -> ScheduledTrajectory:
        """Integrate this scenario, **honouring its timed events** (decision D-35).

        The single "run a compiled scenario" entry point. It always dispatches through
        :func:`~fermentation.runtime.simulate_scheduled` with the compiled ``events``, so a
        temperature ramp (and, from D-36, a dosing/pitching schedule) is applied correctly —
        a multi-knot ramp changes slope at its breakpoints, a hold holds. With no events this
        is byte-for-byte a plain :func:`~fermentation.runtime.simulate` (an empty schedule is a
        single segment), so it is always the right call.

        This exists because a hand-wired ``simulate(cs.process_set, cs.param_values, cs.y0,
        cs.t_span_h)`` **silently ignores** ``events``: the injected ``temperature_ramp_rate``
        would then apply the *first* segment's slope for the whole run — correct only for a
        single-slope ramp. ``param_tiers`` defaults to the scenario's own tier map for honest
        D-1 reporting; ``t_eval``/solver kwargs pass straight through.

        The stochastic sibling is :meth:`run_ensemble`, which threads the same ``events`` into
        :func:`~fermentation.runtime.simulate_ensemble` (decision D-37).
        """
        kwargs.setdefault("param_tiers", self.parameters.tier_map())
        return simulate_scheduled(
            self.process_set,
            self.param_values,
            self.y0,
            self.t_span_h,
            events=self.events,
            **kwargs,  # type: ignore[arg-type]
        )

    def run_ensemble(self, **kwargs: object) -> Ensemble:
        """Run a stochastic ensemble of this scenario, **honouring its timed events** (D-37).

        The uncertainty-band counterpart to :meth:`run`: it hands the compiled ``events`` to
        :func:`~fermentation.runtime.simulate_ensemble` so every sampled member is integrated
        through the same schedule (temperature ramp + dosing/pitching), then reports the nominal
        run plus the median + spread over the parameters' provenance bands. Passes the full
        :class:`~fermentation.parameters.store.ParameterSet` (the ensemble needs the *bands*, not
        just resolved floats) and defaults ``param_tiers`` to the scenario's own tier map. Sampling
        scope, per-member Process-set isolation, and the per-member external-flow ledger are all
        handled by :func:`simulate_ensemble`; ``n_members``/``seed``/``sampler``/``t_eval`` and the
        solver kwargs pass straight through.
        """
        kwargs.setdefault("param_tiers", self.parameters.tier_map())
        return simulate_ensemble(
            self.process_set,
            self.parameters,
            self.y0,
            self.t_span_h,
            events=self.events,
            **kwargs,  # type: ignore[arg-type]
        )


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
        # carrying_capacity_gpl is the optional opt-in biomass cap K that enables the
        # residual-nitrogen floor (D-30); absent ⇒ the cap modifier is disabled (core untouched).
        # citrate_gpl is the optional citric-acid must input O. oeni co-metabolises into
        # MLF-derived diacetyl (D-31); absent ⇒ 0 (no citrate, the diacetyl branch is silent).
        # amino_acids_gpl is the optional assimilable amino-acid dose the AminoAcidAssimilation
        # swap funds biomass from (D-32); absent ⇒ 0 (the swap Process is disabled, core untouched).
        # autolysis_rate_per_h is the optional opt-in first-order autolysis rate (1/h) that both
        # enables YeastAutolysis and overrides its k_autolysis reference (D-34); absent/0 ⇒ the
        # autolysis Process is disabled (an undosed wine run is byte-for-byte the validated core).
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
            "carrying_capacity_gpl",
            "citrate_gpl",
            "amino_acids_gpl",
            "autolysis_rate_per_h",
            "hydroxycinnamic_gpl",
            "brett_pitch_gpl",
            # pof_positive is a binary POF+ strain opt-in (decision D-40 pt4): present/>0 enables
            # YeastPOFDecarboxylation (the yeast cinnamate decarboxylase filling vinylphenols during
            # AF), WHOLLY INDEPENDENT of brett_pitch_gpl. Absent/0 ⇒ a POF-negative wine — the
            # Process stays disabled and the run is byte-for-byte the validated core. Not a state
            # slot: it is a compile-time gate only, never packed into y0.
            "pof_positive",
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
        # the MLF Processes entirely when this is 0 (tier + perf isolability).
        "X_mlf": _optional(values, "mlf_pitch_gpl", 0.0),
        # Citric acid must input (decision D-31); g/L, default 0 (no citrate ⇒ no MLF-derived
        # diacetyl). O. oeni co-metabolises it into α-acetolactate feeding the shared VDK
        # reservoir. Carbon-active (weighted in total_carbon) but not charge-active (kept out
        # of the D-18 pH balance in v1); inert at 0, so an un-dosed run is unchanged.
        "citrate": _optional(values, "citrate_gpl", 0.0),
        # Assimilable amino-acid dose (decision D-32); g/L, default 0 (no amino-acid ledger).
        # Carbon- AND nitrogen-bearing (arginine); the AminoAcidAssimilation swap funds a
        # fraction of biomass from it, refunding sugar + ammonium N. Inert at 0 and the compile
        # step below disables the swap Process entirely when this is 0, so an undosed run is
        # byte-for-byte the validated core (tier + perf isolability, the MLF/carrying pattern).
        "amino_acids": _optional(values, "amino_acids_gpl", 0.0),
        # Non-assimilable cell-wall debris pool (decision D-34); produced-only, empty at pitch.
        # YeastAutolysis routes the carbon-rich remainder of autolysed dead biomass here after
        # releasing the nitrogen-rich amino acids; inert (weight 0) until autolysis is opted in.
        "debris": 0.0,
        # Lumped hydroxycinnamic-acid must precursors (decision D-40); g/L, default 0. Real must
        # carries ~10-200 mg/L (p-coumaric + ferulic); defaulted 0 for isolability (an undosed
        # run is byte-for-byte the validated core), so a Brett scenario doses it. Carbon-active
        # (weighted in total_carbon as p-coumaric); Brettanomyces decarboxylates it to vinylphenols.
        "hydroxycinnamics": _optional(values, "hydroxycinnamic_gpl", 0.0),
        "vinylphenols": 0.0,  # shared decarboxylase→reductase intermediate, empty at pitch (D-40)
        "ethylphenols": 0.0,  # terminal Brett volatile-phenol readout, empty at pitch (D-40)
        # Brettanomyces dose driving the volatile-phenol spoilage (decision D-40); g/L, default 0
        # (no Brett). Constant inert catalyst in pt1 and carbon-free, so an undosed run is
        # byte-for-byte the validated core; the compile step below disables the Brett Processes
        # entirely when this is 0 (tier + perf isolability, the mlf_pitch_gpl pattern).
        "X_brett": _optional(values, "brett_pitch_gpl", 0.0),
        "X_brett_dead": 0.0,  # non-viable Brett lees, empty until BrettDeath (D-40 pt3)
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
    # additions.yaml carries the industry-unit → canonical conversion constants the
    # discrete-intervention verb registry reads at this boundary (decision D-36); like the
    # others it is medium-agnostic and collision-free.
    shared_files = [
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "hydrogen_sulfide.yaml",
        base / "additions.yaml",
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


def _override_carrying_capacity(parameters: ParameterSet, cap_gpl: float) -> ParameterSet:
    """Override the reference ``biomass_carrying_capacity`` with a scenario opt-in value.

    Only reached when a wine scenario passes ``carrying_capacity_gpl > 0`` (decision D-30),
    so the biomass cap modifier is enabled and its cap ``K`` is the scenario's value rather
    than the YAML reference — letting a demonstration sweep the cap. Keeps the reference
    parameter's units/tier/uncertainty (the form and confidence are unchanged; only the
    operating point moves) and records the override in provenance for the audit trail.
    """
    base = parameters["biomass_carrying_capacity"]
    override = Parameter(
        name="biomass_carrying_capacity",
        value=cap_gpl,
        unit=base.unit,
        tier=base.tier,
        uncertainty=base.uncertainty,
        provenance=Provenance(
            source=base.provenance.source,
            doi=base.provenance.doi,
            conditions=(
                f"scenario opt-in override (decision D-30): carrying_capacity_gpl={cap_gpl:g} g/L, "
                f"replacing the {base.value:g} g/L YAML reference"
            ),
            notes=base.provenance.notes,
        ),
    )
    return parameters.merge(ParameterSet([override]), override=True)


def _override_autolysis_rate(parameters: ParameterSet, rate_per_h: float) -> ParameterSet:
    """Override the reference ``k_autolysis`` with a scenario opt-in value (decision D-34).

    Reached only when a wine scenario passes ``autolysis_rate_per_h > 0`` (which enables the
    autolysis Process), so the rate is the scenario's value rather than the YAML reference —
    letting a demonstration sweep the *sur lie* timescale. Keeps the reference parameter's
    units/tier/uncertainty (only the operating point moves) and records the override in
    provenance, mirroring :func:`_override_carrying_capacity`.
    """
    base = parameters["k_autolysis"]
    override = Parameter(
        name="k_autolysis",
        value=rate_per_h,
        unit=base.unit,
        tier=base.tier,
        uncertainty=base.uncertainty,
        provenance=Provenance(
            source=base.provenance.source,
            doi=base.provenance.doi,
            conditions=(
                f"scenario opt-in override (decision D-34): autolysis_rate_per_h={rate_per_h:g} "
                f"1/h, replacing the {base.value:g} 1/h YAML reference"
            ),
            notes=base.provenance.notes,
        ),
    )
    return parameters.merge(ParameterSet([override]), override=True)


def _temperature_ramp_schedule(
    scenario: Scenario, t_end_h: float
) -> tuple[float, tuple[ScheduledEvent, ...]]:
    """Compile ``temperature_schedule`` into an initial slope + slope-change events.

    A temperature schedule is a piecewise-*linear* ramp between its ``(day, celsius)``
    knots (decision D-35): between two knots ``dT/dt`` is a single constant, and the
    :class:`~fermentation.core.kinetics.temperature.TemperatureRamp` Process writes that
    constant into ``dT/dt``. This converts the knots (industry units) into canonical hours
    and Kelvin and returns ``(initial_slope, events)`` where ``events`` restart the
    integrator only at the interior knots where the slope **changes** — so collinear knots
    (a straight ramp given by three points) produce a single segment, and a flat or
    single-knot schedule produces no events and a zero initial slope (isothermal). ``T`` is
    held at the nearest knot's value outside the schedule's span (slope 0 before the first
    knot and after the last).
    """
    knots = sorted(scenario.temperature_schedule, key=lambda p: p.day)
    times = [days_to_hours(p.day) for p in knots]
    temps = [celsius_to_kelvin(p.celsius) for p in knots]

    def slope_after(t: float) -> float:
        # Slope of the segment starting at t. Held flat outside the schedule's span.
        if t < times[0] or t >= times[-1]:
            return 0.0
        for i in range(len(times) - 1):
            if times[i] <= t < times[i + 1]:
                dt = times[i + 1] - times[i]
                return (temps[i + 1] - temps[i]) / dt if dt > 0.0 else 0.0
        return 0.0

    initial_slope = slope_after(0.0)
    events: list[ScheduledEvent] = []
    prev = initial_slope
    for bt in sorted({t for t in times if 0.0 < t < t_end_h}):
        s = slope_after(bt)
        # Only a genuine slope change opens a new segment (collinear knots do not); isclose
        # absorbs float noise so a straight multi-point ramp stays one segment.
        if not math.isclose(s, prev, rel_tol=1e-12, abs_tol=1e-15):
            events.append(
                ScheduledEvent(
                    time_h=bt,
                    label=f"temperature_ramp@{bt / 24.0:g}d",
                    param_update={RAMP_RATE: s},
                )
            )
            prev = s
    return initial_slope, tuple(events)


def _inject_temperature_ramp_rate(parameters: ParameterSet, slope_k_per_h: float) -> ParameterSet:
    """Register the scenario's initial temperature-ramp slope as a provenance-backed parameter.

    Reached only when the schedule actually ramps (decision D-35). The rate is a
    scenario-*exact* set-point forcing, not an empirical kinetic constant, so it is
    VALIDATED with a zero-width band (never swept by the ensemble) — but prime directive #2
    still requires it to travel as a :class:`Parameter` with provenance, so it is minted
    here at the boundary (the D-14/D-30/D-34 pattern) rather than inlined. This value seeds
    the first segment; later segments' slopes are supplied by the ``simulate_scheduled``
    events. ``TemperatureRamp`` reads it with a ``0.0`` default, so an un-ramped scenario
    needs no such parameter at all.
    """
    param = Parameter(
        name=RAMP_RATE,
        value=slope_k_per_h,
        unit="K/h",
        tier=Tier.VALIDATED,
        uncertainty=Uncertainty(
            low=slope_k_per_h,
            high=slope_k_per_h,
            note="scenario-exact temperature set-point schedule; not an uncertain parameter",
        ),
        provenance=Provenance(
            source="scenario temperature_schedule",
            conditions=(
                f"initial piecewise-linear ramp slope {slope_k_per_h:g} K/h (decision D-35)"
            ),
            notes="later intervals' slopes are supplied by simulate_scheduled events",
        ),
    )
    return parameters.merge(ParameterSet([param]), override=True)


# -- discrete-intervention verb registry (the winemaking vocabulary boundary) -----------------
#
# ``scenario.interventions`` is a declarative timeline of winemaking verbs in industry units
# (decision D-36). This registry is where each verb's *meaning* lives — which canonical state
# slot a dose lands on, which unit conversion applies, which Processes a pitch enables — exactly
# the layer that owns the initial-composition vocabulary and the temperature-schedule compile
# above (decision D-3). The runtime driver (``simulate_scheduled``) stays verb-agnostic: a verb
# compiles to an opaque :class:`ScheduledEvent` and the driver just segments-and-restarts around
# it, booking each state jump as an :class:`~fermentation.runtime.schedule.ExternalFlow` for the
# conservation ledger. New verbs are added here and nowhere else.


def _iv_check_keys(iv: Intervention, allowed: frozenset[str], verb: str) -> None:
    unknown = set(iv.params) - allowed
    if unknown:
        raise ValueError(
            f"intervention {verb!r} at day {iv.day:g} has unknown param(s) {sorted(unknown)}; "
            f"allowed: {sorted(allowed)}"
        )


def _iv_float(iv: Intervention, key: str, verb: str) -> float:
    """Read a required numeric intervention param, non-negative (the ``_nonneg`` discipline)."""
    if key not in iv.params:
        raise ValueError(
            f"intervention {verb!r} at day {iv.day:g} is missing required param {key!r}"
        )
    try:
        value = float(iv.params[key])
    except (TypeError, ValueError):
        raise ValueError(
            f"intervention {verb!r} param {key!r} must be a number, got {iv.params[key]!r}"
        ) from None
    return _nonneg(value, key)


def _verb_add_dap(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_dap`` — dose diammonium phosphate, injecting assimilable nitrogen (decision D-36).

    Doses DAP by mass (``dap_gpl``, faithful to the commercial additive) and converts to the
    assimilable-nitrogen jump on the lumped ``N`` slot via the sourced ``dap_nitrogen_fraction``
    (exact (NH₄)₂HPO₄ stoichiometry, VALIDATED). Phosphate is dropped — the model tracks no
    phosphorus pool. The headline consequence is a *timing* effect the static D-29 lever could
    not produce: restoring N mid-ferment momentarily closes the inverse H₂S gate
    ``K_h2s_n/(K_h2s_n+N)`` while sugar (hence the flux the gate multiplies) is still present.
    """
    _iv_check_keys(iv, frozenset({"dap_gpl"}), "add_dap")
    dap_gpl = _iv_float(iv, "dap_gpl", "add_dap")
    try:
        n_fraction = parameters["dap_nitrogen_fraction"].value
    except KeyError as exc:  # additions.yaml not loaded (caller-supplied parameter_paths)
        raise ValueError(
            "intervention 'add_dap' needs 'dap_nitrogen_fraction' but it is missing "
            f"({exc}); include additions.yaml in parameter_paths (the default lookup merges "
            "it automatically)."
        ) from None
    added_n_gpl = dap_gpl * n_fraction
    n_slice = schema.slice("N")

    def mutate(_schema: StateSchema, y: FloatArray) -> FloatArray:
        out = y.copy()
        out[n_slice] += added_n_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_dap@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_add_so2(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_so2`` — dose total SO₂ onto the conserved ``so2_total`` slot (decision D-36).

    Doses total sulfur dioxide by the industry unit (``so2_mgl``, mg/L) and converts to the
    canonical g/L jump on ``so2_total`` — the same slot the initial ``so2_total_mgl`` addition
    lands on (D-22/D-28). Free/bound/molecular SO₂ are then re-derived at the solved pH from that
    total (D-28), so a mid-ferment addition raises the antimicrobial molecular fraction from that
    time forward. SO₂ carries neither carbon nor nitrogen, so this flow perturbs neither elemental
    ledger — the single-run carbon and nitrogen balances still close with no correction term.
    """
    _iv_check_keys(iv, frozenset({"so2_mgl"}), "add_so2")
    so2_mgl = _iv_float(iv, "so2_mgl", "add_so2")
    if "so2_total" not in schema:
        raise ValueError(
            f"intervention 'add_so2' at day {iv.day:g} needs a 'so2_total' slot, but medium "
            f"{schema!r} has none (SO₂ is a wine-only pool, decision D-22)"
        )
    added_gpl = mgl_to_gpl(so2_mgl)
    so2_slice = schema.slice("so2_total")

    def mutate(_schema: StateSchema, y: FloatArray) -> FloatArray:
        out = y.copy()
        out[so2_slice] += added_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_so2@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_add_copper(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_copper`` — copper-fine reductive sulfur (H₂S + mercaptans) out of the wine (D-44/D-45).

    The remediation half of the reductive-fault beat. Copper (Cu²⁺, dosed as copper sulfate)
    precipitates the dissolved reductive-sulfur compounds as insoluble copper salts that settle out
    with the lees — the standard fix for the sur-lie "reduction"
    :class:`~fermentation.core.kinetics.hydrogen_sulfide.AutolyticHydrogenSulfide` (D-44) and
    :class:`~fermentation.core.kinetics.mercaptans.AutolyticMercaptan` (D-45) build up un-stripped
    post-dryness. Doses copper by the industry unit (``copper_mgl``, mg/L Cu) and binds, in **order
    of affinity**:

    1. **H₂S first** — copper sulfide (Cu²⁺ + H₂S → CuS↓ + 2 H⁺, **1:1 mol**), CuS being far more
       insoluble (Ksp ~10⁻³⁶) than the mercaptide, so sulfide is bound preferentially. Capacity
       ``copper·copper_h2s_binding``; removes ``min(h2s, capacity)``.
    2. **Mercaptans with the leftover copper** — copper mercaptide (Cu²⁺ + 2 RSH → Cu(SR)₂↓ + 2 H⁺,
       **1:2 mol**, so a gram of Cu binds ~2.8× the thiol mass it does sulfide), capacity
       ``copper_left·copper_mercaptan_binding``; removes ``min(mercaptans, capacity)``.

    Copper in excess simply clears all dissolved reductive sulfur (the real outcome). **Ledger:**
    H₂S is carbon-free (D-29), so its removal perturbs neither elemental balance (the ``add_so2``
    precedent); **mercaptans carry carbon** (methanethiol, D-45), so removing them removes carbon
    from the wine as the precipitated mercaptide — a **negative external flow** the driver books
    (the racking-debris precedent), so the run-wide identity ``final == initial + Σ flows`` still
    holds. **On a default (autolysis-off) wine ``mercaptans ≡ 0``, so add_copper is carbon-neutral
    there** — the carbon flow appears only once the D-45 pool is non-empty. SCOPE (v1): the removal
    lever only. Residual copper (excess Cu left in the wine) is untracked; copper is imperfect on
    mercaptans and useless on the disulfides they oxidise to (see ``copper_mercaptan_binding``). The
    ``mercaptans`` slot is wine-only, so on a medium without it copper binds H₂S alone.
    """
    _iv_check_keys(iv, frozenset({"copper_mgl"}), "add_copper")
    copper_mgl = _iv_float(iv, "copper_mgl", "add_copper")
    if "h2s" not in schema:  # both current media carry h2s; guard for symmetry with add_so2
        raise ValueError(
            f"intervention 'add_copper' at day {iv.day:g} needs an 'h2s' slot, but medium "
            f"{schema!r} has none"
        )
    try:
        binding_h2s = parameters["copper_h2s_binding"].value
    except KeyError as exc:  # additions.yaml not loaded (caller-supplied parameter_paths)
        raise ValueError(
            "intervention 'add_copper' needs 'copper_h2s_binding' but it is missing "
            f"({exc}); include additions.yaml in parameter_paths (the default lookup merges "
            "it automatically)."
        ) from None
    copper_gpl = mgl_to_gpl(copper_mgl)  # g/L Cu dosed
    h2s_slice = schema.slice("h2s")
    # Mercaptans are wine-only; bind them with leftover copper iff the slot exists (D-45).
    has_mercaptans = "mercaptans" in schema
    binding_merc = parameters["copper_mercaptan_binding"].value if has_mercaptans else 0.0
    merc_slice = schema.slice("mercaptans") if has_mercaptans else None

    def mutate(_schema: StateSchema, y: FloatArray) -> FloatArray:
        out = y.copy()
        # 1. H₂S first (higher affinity). Clamp present ≥ 0 so a solver undershoot is not "removed".
        h2s_present = max(float(out[h2s_slice][0]), 0.0)
        removed_h2s = min(h2s_present, copper_gpl * binding_h2s)
        out[h2s_slice] = float(out[h2s_slice][0]) - removed_h2s
        # 2. Mercaptans with the copper left after binding H₂S (its stoichiometric share).
        if merc_slice is not None:
            copper_left = max(copper_gpl - removed_h2s / binding_h2s, 0.0)
            merc_present = max(float(out[merc_slice][0]), 0.0)
            removed_merc = min(merc_present, copper_left * binding_merc)
            out[merc_slice] = float(out[merc_slice][0]) - removed_merc
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_copper@{iv.day:g}d",
        mutate=mutate,
    )


#: The lees-associated pools racking removes: inactivated yeast biomass ``X_dead`` and (if autolysis
#: is opted in) the non-assimilable cell-wall ``debris`` (decision D-36); plus **both** *O. oeni*
#: pools — settled dead ``X_mlf_dead`` **and viable ``X_mlf``** (decision D-39). Racking viable
#: bacteria is the deliberate **asymmetry with yeast**: a rack leaves viable yeast ``X`` untouched
#: because it works in *suspension* during AF (racking gross lees leaves it fermenting), but
#: *O. oeni* carries out MLF *on the lees* and is drawn off with them — so racking removes the
#: bacteria that clear diacetyl, the physical twin of the SO₂ kill (the other half of the D-31
#: "rack/SO₂ locks in diacetyl" lever, D-39). Both bacterial pools carry biomass carbon and nitrogen
#: (weighted since D-38), so their removal books a negative C/N external flow like ``X_dead``.
#: Dissolved species (sugar, ethanol, acids, glycerol, byproducts, SO₂, YAN) stay with the
#: racked-off liquid — a concentration model has no volume change on racking, so touching them
#: would be physically wrong (decision D-36).
#:
#: **Both *Brettanomyces* pools** — viable ``X_brett`` and settled dead ``X_brett_dead`` — are
#: racked off too (decision D-40), the same lees-organism asymmetry as *O. oeni*: Brett colonises
#: the lees, so drawing the wine off them removes the spoilage catalyst and halts volatile-phenol
#: production — the physical twin of the SO₂ kill (:class:`~fermentation.core.kinetics.brett.\
#: BrettDeath`, pt3). ``X_brett`` is carbon-free in pt1 (constant catalyst, unweighted), so its
#: removal books no C/N flow yet; :class:`~fermentation.core.kinetics.brett.BrettGrowth` (pt2)
#: promotes both pools to weighted biomass, and the same ExternalFlow machinery then books their
#: removal like ``X_mlf`` (D-38/D-39).
_LEES_SLOTS = ("X_dead", "debris", "X_mlf", "X_mlf_dead", "X_brett", "X_brett_dead")


def _verb_rack(iv: Intervention, schema: StateSchema, parameters: ParameterSet) -> ScheduledEvent:
    """``rack`` — draw the wine off a fraction of its lees (decisions D-36, D-39).

    Removes ``fraction`` ∈ [0, 1] of each lees-associated pool (:data:`_LEES_SLOTS`: inactivated
    yeast ``X_dead`` and, when autolysis is opted in, the cell-wall ``debris``; plus both *O. oeni*
    pools — viable ``X_mlf`` and settled dead ``X_mlf_dead``), booking the negative jump as an
    :class:`~fermentation.runtime.schedule.ExternalFlow` (the ledger's removal side). Viable
    **yeast** ``X`` and every dissolved species are left untouched — a normal post-AF rack settles
    dead yeast, and a concentration model has no volume change on racking.

    **Racking removes viable *O. oeni* (the D-39 asymmetry, the D-31 lever's physical half).**
    Unlike viable yeast — which ferments in suspension, so a rack leaves it working — *O. oeni*
    carries out MLF on the lees and is drawn off with them. So racking removes the bacteria that
    clear diacetyl, the physical twin of the SO₂ kill (:class:`~fermentation.core.kinetics.\
    malolactic.MalolacticDeath`): the deferred D-31 "rack early ⇒ diacetyl locked in" case. Both
    bacterial pools carry biomass carbon and nitrogen (weighted since D-38), so — like ``X_dead`` —
    their removal shows up as a negative term in the run-wide carbon and nitrogen ledgers
    (``X_dead``/``X_mlf``/``X_mlf_dead`` all carry N; every racked pool carries C).
    """
    _iv_check_keys(iv, frozenset({"fraction"}), "rack")
    fraction = _iv_float(iv, "fraction", "rack")
    if fraction > 1.0:
        raise ValueError(
            f"intervention 'rack' at day {iv.day:g}: fraction must be in [0, 1], got {fraction:g}"
        )
    slices = [schema.slice(name) for name in _LEES_SLOTS if name in schema]
    retained = 1.0 - fraction

    def mutate(_schema: StateSchema, y: FloatArray) -> FloatArray:
        out = y.copy()
        for sl in slices:
            out[sl] *= retained
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"rack@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_pitch_mlf(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``pitch_mlf`` — inoculate *Oenococcus oeni* mid-run, enabling malolactic fermentation.

    The verb that exercises the driver's third effect (in-place reconfiguration): it both
    **mutates** ``X_mlf`` (adds the bacterial catalyst dose, ``pitch_gpl`` g/L) and
    **reconfigures** the Process set to enable :data:`_MLF_GATED_PROCESSES` — *exactly* the set an
    unpitched compile disables, so a sequential mid-run pitch is symmetric with an initial
    co-inoculation. Since D-38 ``X_mlf`` is real bacterial biomass carrying carbon and nitrogen (no
    longer the inert catalyst of D-23), so the pitch's state jump adds biomass C/N — booked as an
    :class:`~fermentation.runtime.schedule.ExternalFlow`, exactly like the ``add_dap`` dose, so the
    run-wide ledgers still close.

    Because the Processes are enabled only from the breakpoint onward, ``simulate_scheduled``
    min-combines the per-segment tier maps (D-35): the malate/lactate/citrate slots the enabled
    speculative Processes touch report speculative for the *whole* run, not just the back half.

    Honest scope (decisions D-23, D-31): a *post-AF* pitch lands past the Luong ethanol wall
    (~110 g/L; a 24-Brix wine finishes near ~135), so the environmental gate keeps conversion
    near zero — malolactic must be co-inoculated or pitched early to complete. This verb makes
    that timing a *scenario* choice; it does not change the kinetics.
    """
    _iv_check_keys(iv, frozenset({"pitch_gpl"}), "pitch_mlf")
    pitch_gpl = _iv_float(iv, "pitch_gpl", "pitch_mlf")
    if "X_mlf" not in schema:
        raise ValueError(
            f"intervention 'pitch_mlf' at day {iv.day:g} needs an 'X_mlf' slot, but medium "
            f"{schema!r} has none (malolactic fermentation is wine-only, decision D-23)"
        )
    x_mlf_slice = schema.slice("X_mlf")
    gated_names = tuple(p.name for p in _MLF_GATED_PROCESSES)

    def mutate(_schema: StateSchema, y: FloatArray) -> FloatArray:
        out = y.copy()
        out[x_mlf_slice] += pitch_gpl
        return out

    def reconfigure(ps: ProcessSet) -> None:
        for name in gated_names:
            if name in ps:
                ps.enable(name)

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"pitch_mlf@{iv.day:g}d",
        mutate=mutate,
        reconfigure=reconfigure,
    )


def _verb_pitch_brett(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``pitch_brett`` — inoculate *Brettanomyces* mid-run, enabling volatile-phenol spoilage.

    The Brett twin of :func:`_verb_pitch_mlf`: it both **mutates** ``X_brett`` (adds the spoilage
    dose, ``pitch_gpl`` g/L) and **reconfigures** the Process set to enable
    :data:`_BRETT_GATED_PROCESSES` — *exactly* the set an unpitched compile disables, so a mid-run
    contamination is symmetric with an initial co-inoculation. The realistic Brett scenario is a
    *post-AF* contamination in the cellar/barrel, which this verb expresses as a scenario choice.

    Because the Processes are enabled only from the breakpoint onward, ``simulate_scheduled``
    min-combines the per-segment tier maps (D-35): the ``vinylphenols``/``ethylphenols`` slots the
    enabled speculative Processes touch report speculative for the whole run. ``X_brett`` is a
    carbon-free catalyst in pt1, so the pitch's state jump adds no biomass carbon/nitrogen (no
    ExternalFlow needed); :class:`~fermentation.core.kinetics.brett.BrettGrowth` (pt2) makes it
    weighted biomass, at which point the pitch books an :class:`~fermentation.runtime.schedule.\
    ExternalFlow` like ``pitch_mlf`` (D-38).
    """
    _iv_check_keys(iv, frozenset({"pitch_gpl"}), "pitch_brett")
    pitch_gpl = _iv_float(iv, "pitch_gpl", "pitch_brett")
    if "X_brett" not in schema:
        raise ValueError(
            f"intervention 'pitch_brett' at day {iv.day:g} needs an 'X_brett' slot, but medium "
            f"{schema!r} has none (Brettanomyces spoilage is wine-only, decision D-40)"
        )
    x_brett_slice = schema.slice("X_brett")
    gated_names = tuple(p.name for p in _BRETT_GATED_PROCESSES)

    def mutate(_schema: StateSchema, y: FloatArray) -> FloatArray:
        out = y.copy()
        out[x_brett_slice] += pitch_gpl
        return out

    def reconfigure(ps: ProcessSet) -> None:
        for name in gated_names:
            if name in ps:
                ps.enable(name)

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"pitch_brett@{iv.day:g}d",
        mutate=mutate,
        reconfigure=reconfigure,
    )


#: action verb → compiler turning one :class:`Intervention` into a :class:`ScheduledEvent`.
_INTERVENTION_VERBS: dict[
    str, Callable[[Intervention, StateSchema, ParameterSet], ScheduledEvent]
] = {
    "add_dap": _verb_add_dap,
    "add_so2": _verb_add_so2,
    "add_copper": _verb_add_copper,
    "rack": _verb_rack,
    "pitch_mlf": _verb_pitch_mlf,
    "pitch_brett": _verb_pitch_brett,
}


def _compile_interventions(
    scenario: Scenario, schema: StateSchema, parameters: ParameterSet, t_end_h: float
) -> tuple[ScheduledEvent, ...]:
    """Compile ``scenario.interventions`` into timed :class:`ScheduledEvent`\\ s (decision D-36).

    Each verb is looked up in :data:`_INTERVENTION_VERBS`; an unknown action raises loudly (the
    ``_ALLOWED_KEYS`` discipline). An intervention at or after the run duration is rejected here
    with a scenario-level message rather than deferred to ``simulate_scheduled``'s window check,
    so the error names the scenario and the verb.
    """
    events: list[ScheduledEvent] = []
    for iv in scenario.interventions:
        verb = _INTERVENTION_VERBS.get(iv.action)
        if verb is None:
            raise ValueError(
                f"scenario {scenario.name!r}: unknown intervention action {iv.action!r}; "
                f"known verbs: {sorted(_INTERVENTION_VERBS)}"
            )
        if days_to_hours(iv.day) >= t_end_h:
            raise ValueError(
                f"scenario {scenario.name!r}: intervention {iv.action!r} at day {iv.day:g} is "
                f"at or beyond the run duration ({scenario.duration_days:g} d); interventions "
                "must fall within the run"
            )
        events.append(verb(iv, schema, parameters))
    return tuple(events)


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

    # MLF isolability (decisions D-23, D-31): the malolactic Processes are wired into the wine
    # medium but contribute nothing until Oenococcus oeni is pitched. When it is not, DISABLE
    # them all so (a) the inert ``malic``/``lactic``/``citrate`` slots keep their VALIDATED
    # tier — an *enabled* Process that touches them drops them to speculative even with a zero
    # contribution, since ``tier_of`` counts enabled, not nonzero, Processes — and (b) no
    # per-RHS pH ``brentq`` solve is paid on an undosed run. When pitched, MalolacticConversion
    # is the first RHS consumer of the D-18 pH solver / D-22 molecular-SO₂ readout, and the two
    # D-31 Processes co-metabolise citrate into diacetyl and reduce it on the lees.
    # An initial ``mlf_pitch_gpl`` co-inoculates at t0; a mid-run ``pitch_mlf`` intervention
    # (decision D-36) instead leaves this 0 and re-enables the same _MLF_GATED_PROCESSES at its
    # breakpoint. Either way, an unpitched compile disables them here.
    mlf_pitch_gpl = float(scenario.initial.get("mlf_pitch_gpl", 0.0) or 0.0)
    if mlf_pitch_gpl <= 0.0:
        for mlf_process in _MLF_GATED_PROCESSES:
            if mlf_process.name in process_set:
                process_set.disable(mlf_process.name)

    # Brett isolability (decision D-40): the volatile-phenol Processes are wired into the wine
    # medium but contribute nothing until Brettanomyces is pitched. When it is not, DISABLE them so
    # (a) the inert ``hydroxycinnamics``/``vinylphenols``/``ethylphenols`` slots keep their
    # VALIDATED tier (an *enabled* zero-contribution Process still drags them to speculative via
    # ``tier_of``) and (b) no per-RHS pH ``brentq`` is paid on an unpitched run. An initial
    # ``brett_pitch_gpl`` co-inoculates at t0; a mid-run ``pitch_brett`` intervention instead leaves
    # this 0 and re-enables the same _BRETT_GATED_PROCESSES at its breakpoint (the MLF pattern).
    brett_pitch_gpl = float(scenario.initial.get("brett_pitch_gpl", 0.0) or 0.0)
    if brett_pitch_gpl <= 0.0:
        for brett_process in _BRETT_GATED_PROCESSES:
            if brett_process.name in process_set:
                process_set.disable(brett_process.name)

    # POF+ yeast decarboxylase isolability (decision D-40 pt4): YeastPOFDecarboxylation is the yeast
    # cinnamate decarboxylase filling ``vinylphenols`` from must ``hydroxycinnamics`` during AF.
    # POF+ is a BINARY STRAIN TRAIT, gated on its own opt-in ``pof_positive`` and WHOLLY INDEPENDENT
    # of the Brett pitch (a POF+ ferment need not have Brett; a POF-negative wine must make no
    # vinylphenol). Absent/<=0 => DISABLE it so (a) the empty ``vinylphenols`` slot keeps its
    # VALIDATED tier - an *enabled* zero-contribution Process would drag it to speculative via
    # ``tier_of`` - and (b) no wasted flux/Monod recompute is paid on a POF- run, which is then
    # byte-for-byte the validated core (the Brett-unpitched pattern). Opted in => the Process runs;
    # ``vinylphenols`` honestly reports speculative, while ``ethylphenols`` stays VALIDATED at 0
    # unless Brett (the only reductase) is also present - the emergent stranding.
    pof_positive = float(scenario.initial.get("pof_positive", 0.0) or 0.0)
    if pof_positive <= 0.0 and YeastPOFDecarboxylation.name in process_set:
        process_set.disable(YeastPOFDecarboxylation.name)

    # Residual-nitrogen floor (decision D-30): the biomass carrying-capacity cap is a
    # deliberate DEPARTURE from the validated Coleman anchor (which caps nothing and strips
    # YAN to zero at every dose), so it ships OPT-IN. Absent ``carrying_capacity_gpl`` ⇒
    # DISABLE the modifier so (a) growth's whole contribution is unscaled (factor 1) and the
    # run is byte-for-byte the validated core, and (b) the enabled-but-inert modifier does not
    # drag growth's X/S/N outputs from PLAUSIBLE to speculative (``tier_of`` counts enabled,
    # not nonzero, modifiers — the exact MLF *tier* isolability argument above). Opted in ⇒
    # enable it and override the reference cap with the scenario's value so demonstrations can
    # sweep K; growth's outputs then honestly report speculative.
    if BiomassCarryingCapacity.name in process_set:
        raw_cap = scenario.initial.get("carrying_capacity_gpl")
        # A negative cap is a typo, not an intent — raise loudly like every other initial key
        # (the _nonneg gate), rather than silently disabling. Absent or 0 ⇒ opt out (disable).
        cap_gpl = _nonneg(float(raw_cap), "carrying_capacity_gpl") if raw_cap is not None else 0.0
        if cap_gpl <= 0.0:
            process_set.disable(BiomassCarryingCapacity.name)
        else:
            parameters = _override_carrying_capacity(parameters, cap_gpl)

    # Amino-acid ledger isolability (decisions D-32, D-33): the AminoAcidAssimilation swap and
    # the FuselAminoAcidReroute are wired into the wine medium but contribute nothing until amino
    # acids are dosed. When they are not, DISABLE them so (a) the empty ``amino_acids`` slot keeps
    # its VALIDATED tier — an *enabled* speculative Process touching ``S``/``N`` would drag those
    # outputs down even with a zero contribution (``tier_of`` counts enabled, not nonzero,
    # Processes) — and (b) no rate recompute is paid on an undosed run. Dosed, the swap funds a
    # fraction of biomass from amino acids (refunding sugar carbon and ammonium N, scaled alongside
    # growth by the wine Arrhenius/carrying modifiers so it never creates sugar, D-32), and the
    # re-route sources a fraction of Ehrlich fusel carbon from amino acids, deaminating the nitrogen
    # to ammonium (D-33). The re-route is paired with FuselAlcoholsEhrlich (it refunds sugar that
    # producer drew), which is always enabled in the wine set, so it is safe to enable here.
    amino_acids_gpl = float(scenario.initial.get("amino_acids_gpl", 0.0) or 0.0)
    if amino_acids_gpl <= 0.0:
        for aa_process in (AminoAcidAssimilation, FuselAminoAcidReroute):
            if aa_process.name in process_set:
                process_set.disable(aa_process.name)

    # MLF-growth isolability (decision D-38, the deferred growth beat). MalolacticGrowth builds
    # bacterial biomass X_mlf from the amino-acid pool, so the FEATURE it represents (amino-acid-fed
    # bacterial growth) is keyed on amino acids being dosed — the SAME gate as the swap/re-route
    # above. Disable it when amino_acids_gpl ≤ 0 so (a) the empty amino_acids slot / X_mlf keep
    # their tier — an enabled speculative Process touching them would drag tier_of even with a zero
    # contribution — and (b) no rate recompute is paid undosed. This alone prevents the D-23/D-31
    # tier regression: those tests pitch O. oeni but dose NO amino acids, so growth stays disabled.
    # NOT additionally gated on the pitch: "bacteria present" is runtime state the Process's own
    # ``X_mlf ≤ 0`` guard handles (zero until a co-inoculation dose or a mid-run pitch_mlf mutation
    # adds X_mlf), and — mirroring how MalolacticConversion trusts its ethanol gate rather than a
    # compile rule — whether post-pitch bacteria then GROW is left to the emergent environmental
    # gate (g_EtOH·γ(T)·…). So co-inoculation-dominance is emergent, not hard-coded: a high-ABV must
    # arrests growth via the ethanol wall, a normal-ABV sequential MLF can still grow (D-38).
    if MalolacticGrowth.name in process_set and amino_acids_gpl <= 0.0:
        process_set.disable(MalolacticGrowth.name)

    # Brett-growth isolability (decision D-40 pt2). BrettGrowth builds X_brett from the amino-acid
    # pool (drawing the carbon shortfall from ethanol, so Brett grows in a dry wine), so — exactly
    # like MalolacticGrowth — it is keyed on amino acids being dosed, NOT on the Brett pitch (the
    # Process's own ``X_brett ≤ 0`` guard handles "Brett present"; whether it then grows is left to
    # the emergent SO₂/temperature gate). Disable it when amino_acids_gpl ≤ 0 so the empty
    # amino_acids / X_brett slots keep their tier and no rate recompute is paid undosed. This keeps
    # every pitched-but-not-aa-dosed Brett run (e.g. the pt1 headline) with growth disabled.
    if BrettGrowth.name in process_set and amino_acids_gpl <= 0.0:
        process_set.disable(BrettGrowth.name)

    # Autolytic source (decisions D-34, D-44, D-45): YeastAutolysis refills the amino-acid pool from
    # dead biomass (X_dead) post-AF — the second MLF-with-growth prerequisite — AutolyticHydrogen-
    # Sulfide (D-44) feeds the shared h2s pool the sulfide those self-digesting cells release, and
    # AutolyticMercaptan (D-45) feeds the mercaptans pool their thiols (drawing carbon from
    # amino_acids, deaminating N) — all three yields on the SAME autolysis flux. They *consume/read*
    # core state gated on autolysis, so they ship OPT-IN TOGETHER: absent ``autolysis_rate_per_h`` ⇒
    # DISABLE all three so (a) the X_dead/amino_acids/debris/h2s/mercaptans/N columns are untouched
    # and the run is byte-for-byte the validated core, and (b) the inert
    # Processes do not drag those outputs to speculative (``tier_of`` counts enabled, not nonzero,
    # Processes — the MLF/carrying *tier* isolability argument). Opted in ⇒ enable them and override
    # k_autolysis with the scenario's rate so demonstrations can sweep the sur-lie timescale — the
    # override drives BOTH the peptide refill and the sulfide yield (they read one k_autolysis).
    if YeastAutolysis.name in process_set:
        raw_rate = scenario.initial.get("autolysis_rate_per_h")
        rate_per_h = (
            _nonneg(float(raw_rate), "autolysis_rate_per_h") if raw_rate is not None else 0.0
        )
        if rate_per_h <= 0.0:
            for autolysis_process in (
                YeastAutolysis,
                AutolyticHydrogenSulfide,
                AutolyticMercaptan,
            ):
                if autolysis_process.name in process_set:
                    process_set.disable(autolysis_process.name)
        else:
            parameters = _override_autolysis_rate(parameters, rate_per_h)

    t_span_h = (0.0, days_to_hours(scenario.duration_days))

    # Temperature schedule (decision D-35): compile the piecewise-linear ramp into the
    # TemperatureRamp's initial slope + slope-change events. Only when it actually ramps do
    # we mint the provenance-backed rate parameter and emit events — a flat/single-knot
    # schedule leaves ``temperature_ramp_rate`` absent, so the always-enabled TemperatureRamp
    # reads its 0.0 default and an isothermal run is byte-for-byte the pre-ramp core.
    initial_slope, ramp_events = _temperature_ramp_schedule(scenario, t_span_h[1])
    if initial_slope != 0.0 or ramp_events:
        parameters = _inject_temperature_ramp_rate(parameters, initial_slope)

    # Discrete winemaking interventions (decision D-36): compile the declarative timeline of
    # verbs into timed events and merge with the ramp's slope-change events into the single
    # ``events`` tuple ``simulate_scheduled`` sorts by time. Empty ⇒ the temp-only path is
    # unchanged (byte-for-byte core when there is no ramp either).
    intervention_events = _compile_interventions(scenario, medium.schema, parameters, t_span_h[1])
    events = (*ramp_events, *intervention_events)

    return CompiledScenario(
        scenario=scenario,
        schema=medium.schema,
        y0=y0,
        process_set=process_set,
        parameters=parameters,
        t_span_h=t_span_h,
        events=events,
    )
