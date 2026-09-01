"""Compile-and-run boundary. Everything the interface displays comes from here.

**The one rule this module exists to enforce: a compiled scenario is used for exactly one
run and then thrown away.** ``CompiledScenario.run()`` is not idempotent (decision D-206) —
an event's ``reconfigure`` mutates the Process set in place and nothing puts it back, so a
second ``run()`` on the same object starts with the first run's switches already live at
t = 0, silently and with no error. Under a UI framework that re-executes its script on every
widget change, and a cache that hands the *same object* back on a hit, that is not a
theoretical hazard: it is the default outcome. So :func:`run_once` compiles inside the cached
function, every single time, and the cache stores the finished, inert result.

Nothing here holds a live ``ProcessSet`` either. The mechanism metadata the provenance walk
needs (names, tiers, what each touches, what each reads) is static and is copied out at
compile, from the public ``active`` view taken on either side of the run. Comparing those
two views is the honest way to answer "was aging running?" — it names both the mechanisms
that ran from t = 0 and the ones an event switched on partway — without keeping a mutable
object alive past its run.
"""

from __future__ import annotations

import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from app.fidelity import Fidelity, tightened
from fermentation.core.process import Process, RateModifier
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import ParameterSet
from fermentation.runtime.ensemble import Ensemble
from fermentation.runtime.schedule import ScheduledTrajectory
from fermentation.scenario import Scenario, compile_scenario


@dataclass(frozen=True)
class MechanismInfo:
    """One Process or RateModifier, flattened to plain data at compile time."""

    name: str
    kind: str  # "process" | "modifier"
    tier: Tier
    touches: tuple[str, ...]
    reads: tuple[str, ...]
    modifies: tuple[str, ...] = ()
    enabled_at_start: bool = True
    enabled_at_end: bool = True

    @property
    def switched_on_mid_run(self) -> bool:
        return self.enabled_at_end and not self.enabled_at_start


@dataclass(frozen=True)
class RunResult:
    """A finished deterministic run, plus everything needed to explain it."""

    scenario: Scenario
    fidelity: Fidelity
    traj: ScheduledTrajectory
    parameters: ParameterSet
    param_tiers: Mapping[str, Tier]
    mechanisms: tuple[MechanismInfo, ...]
    wall_seconds: float
    compile_seconds: float
    #: Mechanisms present in the compiled set but switched off for the whole run.
    inactive_names: tuple[str, ...] = ()

    @property
    def schema(self) -> StateSchema:
        return self.traj.schema

    @property
    def days(self) -> FloatArray:
        return np.asarray(self.traj.t, dtype=float) / 24.0

    @property
    def param_values(self) -> dict[str, float]:
        return self.parameters.resolve()

    def touched_variables(self) -> frozenset[str]:
        """State variables some mechanism in this run actually writes to.

        Anything outside this set held whatever the compile seam seeded it with for the whole
        run. That matters for how it is *reported*, not just for what it looks like: the
        engine's tier combine returns the top tier for an empty input list, so a variable no
        mechanism touches comes back reading ``validated`` — the one word nothing in this
        engine has earned. An untouched variable is inert, not confirmed, and the interface
        must say so with a different word.
        """
        return frozenset(v for m in self.mechanisms for v in m.touches)

    def touching(self, variable: str) -> tuple[MechanismInfo, ...]:
        """Mechanisms that write ``variable`` (plus the modifiers scaling them)."""
        direct = tuple(m for m in self.mechanisms if variable in m.touches)
        names = {m.name for m in direct}
        scaling = tuple(
            m for m in self.mechanisms if m.kind == "modifier" and names & set(m.modifies)
        )
        return direct + scaling


@dataclass(frozen=True)
class EnsembleResult:
    """A finished uncertainty band. Carries the failure fraction, always."""

    scenario: Scenario
    fidelity: Fidelity
    ensemble: Ensemble
    parameters: ParameterSet
    param_tiers: Mapping[str, Tier]
    n_requested: int
    seed: int
    sampler: str
    wall_seconds: float
    #: Constants restricted to, if the sampling was narrowed. Empty ⇒ everything in scope.
    only: tuple[str, ...] = ()
    warnings: tuple[str, ...] = field(default_factory=tuple)

    @property
    def days(self) -> FloatArray:
        return np.asarray(self.ensemble.t, dtype=float) / 24.0

    @property
    def failure_fraction(self) -> float:
        return float(self.ensemble.failure_fraction)


def scenario_key(scenario: Scenario) -> str:
    """A stable cache key for a scenario. JSON, because the model holds lists and dicts."""
    return scenario.model_dump_json()


def _mechanisms(
    before: Sequence[Process | RateModifier],
    after: Sequence[Process | RateModifier],
    enabled: Mapping[str, bool],
) -> tuple[MechanismInfo, ...]:
    """Flatten the mechanisms live at either end of the run into plain, inert data.

    Built from ``ProcessSet.active`` / ``active_modifiers`` taken before and after the run
    rather than from the set's private registries, so this cannot break on an internal
    rename. A mechanism off at both ends contributed nothing and is reported by name only,
    through :attr:`RunResult.inactive_names`.
    """
    start_names = {m.name for m in before}
    seen: dict[str, MechanismInfo] = {}
    for source in (before, after):
        for m in source:
            if m.name in seen:
                continue
            modifier = isinstance(m, RateModifier)
            seen[m.name] = MechanismInfo(
                name=m.name,
                kind="modifier" if modifier else "process",
                tier=m.tier,
                touches=() if isinstance(m, RateModifier) else tuple(m.touches),
                reads=tuple(m.reads),
                modifies=tuple(m.modifies) if isinstance(m, RateModifier) else (),
                enabled_at_start=m.name in start_names,
                enabled_at_end=bool(enabled.get(m.name, False)),
            )
    return tuple(sorted(seen.values(), key=lambda m: (m.kind, m.name)))


def run_once(scenario: Scenario, fidelity: Fidelity) -> RunResult:
    """Compile this scenario **fresh** and integrate it once.

    Compiling inside this function is the D-206 guard described in the module docstring, not
    an inefficiency to optimise away: compile is a fraction of a second, and reusing a
    compiled object across two runs is a silent wrong answer.
    """
    t0 = time.perf_counter()
    compiled = compile_scenario(scenario, strict=fidelity.strict, oxidative=fidelity.oxidative)
    t1 = time.perf_counter()
    before = compiled.process_set.active + compiled.process_set.active_modifiers
    traj = compiled.run(**fidelity.solver_kwargs(compiled.t_span_h))
    t2 = time.perf_counter()
    after = compiled.process_set.active + compiled.process_set.active_modifiers
    final_enabled = compiled.process_set.enabled_snapshot()
    return RunResult(
        scenario=scenario,
        fidelity=fidelity,
        traj=traj,
        parameters=compiled.parameters,
        param_tiers=compiled.parameters.tier_map(),
        mechanisms=_mechanisms(before, after, final_enabled),
        inactive_names=tuple(sorted(n for n, on in final_enabled.items() if not on)),
        wall_seconds=t2 - t1,
        compile_seconds=t1 - t0,
    )


def varying_constants(result: RunResult, only: Sequence[str] | None = None) -> tuple[str, ...]:
    """Constants that would actually vary in an ensemble of this run.

    A parameter with a zero-width uncertainty range is drawn at its own value every time and
    explains no spread, so it does not count. The reason this is worth projecting *before*
    the ensemble runs: the spread ranking is a regression of the members on the drawn
    parameters, and it is underdetermined unless there are more members than varying
    parameters. On a wine that means ~90 constants and therefore ~90 runs — a minute of
    compute the person deserves to be told about before they press the button, not after.
    """
    names: set[str] = set()
    for m in result.mechanisms:
        names.update(m.reads)
    keep: list[str] = []
    for name in sorted(names):
        if only is not None and name not in only:
            continue
        if name not in result.parameters:
            continue
        u = result.parameters[name].uncertainty
        if u.high > u.low:
            keep.append(name)
    return tuple(keep)


def run_uncertainty(
    scenario: Scenario,
    fidelity: Fidelity,
    *,
    n_members: int = 48,
    seed: int = 0,
    sampler: str = "lhs",
    only: Sequence[str] | None = None,
) -> EnsembleResult:
    """Run a Monte-Carlo band over the parameters' published uncertainty ranges.

    ``n_members`` defaults to 48, **not** the engine's 200: the caller here is a person
    waiting at a screen, and 200 members of a two-week wine is well over a minute. The
    engine's default is right for a batch job and wrong for a button.

    ``max_failure_fraction`` is left at the engine's 0.5, which means an ensemble can lose
    up to half its members and still return a band. That is a legitimate result — some
    parameter draws genuinely do not integrate — but it must never be plotted silently, so
    the surviving-member count travels in the result and every renderer prints it.
    """
    compiled = compile_scenario(scenario, strict=fidelity.strict, oxidative=fidelity.oxidative)
    kwargs = fidelity.solver_kwargs(compiled.t_span_h)
    t0 = time.perf_counter()
    if only is not None:
        kwargs["only"] = list(only)
    ens = compiled.run_ensemble(n_members=n_members, seed=seed, sampler=sampler, **kwargs)
    t1 = time.perf_counter()

    warnings: list[str] = []
    frac = float(ens.failure_fraction)
    if frac > 0.0:
        warnings.append(
            f"{frac:.0%} of the {n_members} sampled runs failed to integrate. The band below "
            f"is computed from the {ens.n_succeeded} that survived, which is a selected sample, "
            "not the whole draw."
        )
    varying = len(ens.sampled_names)
    if n_members <= varying:
        warnings.append(
            f"{varying} constants varied across only {n_members} runs. The band itself is "
            "fine, but the spread ranking cannot be computed with fewer runs than varying "
            "constants — raise the member count above that, or narrow the sampling to the "
            "constants one readout depends on."
        )
    elif n_members < 50:
        warnings.append(
            "Below about 50 members the spread ranking is unstable. Treat the ordering as "
            "indicative until you re-run it larger."
        )
    return EnsembleResult(
        scenario=scenario,
        fidelity=fidelity,
        ensemble=ens,
        parameters=compiled.parameters,
        param_tiers=compiled.parameters.tier_map(),
        n_requested=n_members,
        seed=seed,
        sampler=sampler,
        only=tuple(only) if only is not None else (),
        wall_seconds=t1 - t0,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True)
class ConvergenceCheck:
    """Did the answer move when the solver was leaned on harder?

    ``worst`` is the largest relative change across the readouts checked, scaled by each
    series' own peak so a variable that is near zero everywhere cannot manufacture a huge
    relative difference out of numerical dust.
    """

    tighter: Fidelity
    per_variable: Mapping[str, float]
    worst: float
    worst_variable: str
    seconds: float

    @property
    def converged(self) -> bool:
        return self.worst < 1e-3

    @property
    def verdict(self) -> str:
        if self.worst < 1e-4:
            return "Converged. Tightening the solver changes nothing you could read off a chart."
        if self.worst < 1e-3:
            return "Converged for practical purposes — the answer moves by less than 0.1%."
        if self.worst < 1e-2:
            return "Marginal. The answer moves by up to 1%; quote it at the tighter setting."
        return "Not converged. The current setting is too loose for this scenario."


def check_convergence(result: RunResult, variables: tuple[str, ...]) -> ConvergenceCheck | None:
    """Re-run at the next precision up and report how far the answer moved.

    This is what turns the default from a claim into a demonstration. It costs one extra
    integration and is the only control in the fidelity panel that answers the question the
    panel actually raises: *was the default good enough here?*
    """
    tighter = tightened(result.fidelity)
    if tighter is None:
        return None
    t0 = time.perf_counter()
    ref = run_once(result.scenario, tighter)
    seconds = time.perf_counter() - t0

    per: dict[str, float] = {}
    for name in variables:
        if name not in result.schema or name not in ref.schema:
            continue
        a = np.atleast_2d(np.asarray(result.traj.series(name), dtype=float))
        b = np.atleast_2d(np.asarray(ref.traj.series(name), dtype=float))
        if a.shape != b.shape:
            continue
        scale = float(np.max(np.abs(b)))
        if scale <= 0.0:
            continue
        per[name] = float(np.max(np.abs(a - b)) / scale)
    if not per:
        return None
    worst_variable = max(per, key=lambda k: per[k])
    return ConvergenceCheck(
        tighter=tighter,
        per_variable=per,
        worst=per[worst_variable],
        worst_variable=worst_variable,
        seconds=seconds,
    )
