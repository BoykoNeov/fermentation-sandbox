"""Stochastic ensemble wrapper over the deterministic core (handoff §1.6, D-24).

The core is deterministic: given a state and a resolved ``{name: float}`` parameter
map it returns derivatives, and a single run is byte-for-byte reproducible. Real
ferments vary because their *parameters* are uncertain — every :class:`Parameter`
already carries a provenance-declared :class:`~fermentation.parameters.schema.Uncertainty`
band, and until now nothing at runtime read it. This module closes that loop: it
samples each parameter within its band, runs an ensemble of :func:`simulate_scheduled`
calls (a plain :func:`simulate` when no ``events`` are passed), and reports the
deterministic *nominal* trajectory alongside the ensemble *median* and inter-percentile
*spread*. Passing ``events`` makes it the uncertainty band of a *scheduled* scenario — a
temperature ramp, a DAP dose, a mid-run pitch — replayed under every draw (decision D-37).

**Where randomness lives.** All sampling happens *here*, in the runtime wrapper,
behind an explicit seed. The core stays pure (no ``Math.random`` in a Process), so a
single unsampled run remains reproducible and debuggable — exactly the split the
handoff (§1.6) and the architecture rule require. This wrapper takes the full
:class:`~fermentation.parameters.store.ParameterSet` (it needs the bands), which is
the natural seam distinguishing it from :func:`simulate` (which takes resolved
floats).

**Scope (decision D-24).**

* *Parameter* uncertainty only. Scenario/initial-condition uncertainty (Brix, YAN)
  is a separate axis and is **not** sampled here — ``y0`` is held fixed.
* Plain Monte Carlo (``sampler="mc"``, the method the handoff §1.6 names) by default.
  Latin-hypercube (``"lhs"``) and Sobol (``"sobol"``) low-discrepancy sequences are
  also available: they stratify the draws so a fixed member budget covers the band —
  and especially its tails — more evenly than i.i.d. Monte Carlo. They change only *how*
  the unit hypercube is drawn (via ``scipy.stats.qmc``), then map it through each
  parameter's inverse CDF; the ``only``/``exclude`` scoping and the failed-member /
  survivorship accounting are identical across samplers. Only the *varying* parameters
  take a hypercube dimension (a pinned zero-width band would waste one and degrade Sobol's
  balance). Sobol's balance property holds at powers of two, so a Sobol run requires
  ``n_members`` to be a power of two (it raises otherwise rather than silently returning
  an unbalanced sequence); LHS and MC take any count.
* Default distribution is **triangular** ``(low, mode=value, high)``: bounds plus a
  most-likely value is the textbook case for triangular, and ``value`` *is* the
  sourced, benchmarked most-likely estimate. ``uniform`` is available for a caller
  who wants to discard the point estimate. The reported band uses outer percentiles
  (P5/P95 by default), which keeps the full bracket visible and de-sensitises the
  result to the shape choice.
* By default only the parameters the **active** Process set ``reads`` are sampled
  (the rest are no-ops on the trajectory), so the spread means "sensitivity of *this*
  scenario". ``only`` overrides that set; ``exclude`` removes names from it.
* The ensemble **median is not the nominal trajectory** (the median of nonlinear
  trajectories is not the trajectory of median parameters). Both are reported: the
  nominal is the deterministic reference, the median+band is the uncertainty summary.

**Independence caveat.** Parameters are sampled independently, which ignores any
cross-parameter correlation or ordering constraint. The two live constraint groups
were checked against their actual bands (decision D-24) and found immaterial: the
realised-yield partition cannot exceed the theoretical split (glycerol/byproduct
carbon is *carved from* the Gay-Lussac flux, with a hard guard in the uptake
Process), and the load-bearing ``E_a > E_a_uptake`` byproduct ordering is preserved
for the overwhelming majority of triangular draws (the wine ester direction is
intentionally null; the fusel/uptake band overlap is a negligible tail). A caller who
wants a strictly-ordered ensemble can ``exclude`` the offending names to pin them.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from scipy.stats import qmc, triang

from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import ParameterSet
from fermentation.runtime.schedule import (
    ExternalFlow,
    ScheduledEvent,
    ScheduledTrajectory,
    simulate_scheduled,
)

#: Supported sampling distributions over a parameter's ``[low, high]`` band.
DISTRIBUTIONS = ("triangular", "uniform")

#: Supported ensemble sampling strategies. ``"mc"`` is i.i.d. Monte Carlo; ``"lhs"`` and
#: ``"sobol"`` are low-discrepancy sequences for better (tail) coverage per member.
SAMPLERS = ("mc", "lhs", "sobol")


class Band(NamedTuple):
    """A named variable's uncertainty band over the shared time grid.

    ``low``/``median``/``high`` each mirror :meth:`Trajectory.series`: a scalar
    variable comes back as shape ``(n_times,)``, a vector (multi-slot ``S``) as
    ``(n_slots, n_times)``. ``low``/``high`` are the requested percentiles across the
    surviving ensemble members; ``nominal`` is the deterministic run for reference.
    """

    low: FloatArray
    median: FloatArray
    high: FloatArray
    nominal: FloatArray


def sample_parameters(
    parameters: ParameterSet,
    rng: np.random.Generator,
    *,
    distribution: str = "triangular",
    names: Sequence[str] | None = None,
) -> dict[str, float]:
    """Draw one ``{name: float}`` sample from the parameters' uncertainty bands.

    Every parameter starts at its nominal ``value``; each name in ``names`` (or every
    parameter when ``names`` is ``None``) is then replaced by a draw from its band via
    ``distribution``. A zero-width band (``high <= low``) is pinned to ``value`` — a
    no-op draw that consumes no randomness — so a parameter with no stated spread never
    moves. ``names`` is consumed **in order** (not as a set) so the draw sequence is
    deterministic for a given seed regardless of hash randomisation.

    The :class:`Parameter` schema guarantees ``low <= value <= high``, which is exactly
    ``numpy``'s ``triangular(left, mode, right)`` precondition, so no clamping is needed.
    """
    if distribution not in DISTRIBUTIONS:
        raise ValueError(f"unknown distribution {distribution!r}; expected one of {DISTRIBUTIONS}")
    out = parameters.resolve()
    target = parameters.names if names is None else names
    for name in target:
        p = parameters[name]
        lo, hi, val = p.uncertainty.low, p.uncertainty.high, p.value
        if hi <= lo:  # no stated spread ⇒ pinned, and draw nothing (keeps RNG in step)
            out[name] = val
            continue
        if distribution == "triangular":
            out[name] = float(rng.triangular(lo, val, hi))
        else:  # uniform
            out[name] = float(rng.uniform(lo, hi))
    return out


def _inverse_cdf(q: float, lo: float, val: float, hi: float, distribution: str) -> float:
    """Map a unit-hypercube coordinate ``q ∈ [0, 1]`` to a value in ``[lo, hi]``.

    The inverse CDF (percent-point function) for the chosen band shape: triangular with
    mode ``val`` (``scipy.stats.triang`` parameterised by ``c = (val-lo)/(hi-lo)``), or
    uniform. This is the low-discrepancy counterpart to :func:`sample_parameters`'
    ``rng.triangular``/``rng.uniform`` draws — same distribution, but the coordinate comes
    from a stratified QMC engine instead of a PRNG. Callers pass only varying bands here
    (``hi > lo``), so the ``c`` denominator is never zero.
    """
    if distribution == "triangular":
        return float(triang.ppf(q, (val - lo) / (hi - lo), loc=lo, scale=hi - lo))
    return lo + q * (hi - lo)  # uniform


def _qmc_samples(
    parameters: ParameterSet,
    sampled_names: Sequence[str],
    n_members: int,
    *,
    seed: int,
    distribution: str,
    sampler: str,
) -> list[dict[str, float]]:
    """Build ``n_members`` parameter samples from a low-discrepancy (LHS/Sobol) sequence.

    Only the *varying* sampled parameters (``high > low``) take a hypercube dimension;
    pinned zero-width bands stay at their nominal ``value`` (giving them a QMC column would
    waste a dimension and, for Sobol, break its balance). Each member starts from the full
    nominal ``resolve()`` map (so unsampled and pinned names are at nominal, exactly as in
    the MC path) and has its varying names replaced by the inverse-CDF image of the QMC
    coordinate. The engine is seeded, so the whole matrix — and thus the ensemble — is
    reproducible; because the matrix is fixed up front, member ``i`` always uses row ``i``
    regardless of which members fail.
    """
    varying = [
        n for n in sampled_names if parameters[n].uncertainty.high > parameters[n].uncertainty.low
    ]
    base = parameters.resolve()
    if not varying:  # nothing to stratify ⇒ every member is the nominal run
        return [dict(base) for _ in range(n_members)]

    d = len(varying)
    if sampler == "lhs":
        engine: qmc.QMCEngine = qmc.LatinHypercube(d=d, seed=seed)
    else:  # sobol
        if n_members & (n_members - 1) != 0:  # not a power of two
            raise ValueError(
                f"sobol sampling needs n_members to be a power of two for its balance "
                f"property, got {n_members}; use a power of two, or sampler='lhs'/'mc'"
            )
        engine = qmc.Sobol(d=d, seed=seed)
    unit = engine.random(n_members)  # shape (n_members, d), coordinates in [0, 1)

    samples: list[dict[str, float]] = []
    for i in range(n_members):
        out = dict(base)
        for j, name in enumerate(varying):
            p = parameters[name]
            out[name] = _inverse_cdf(
                float(unit[i, j]), p.uncertainty.low, p.value, p.uncertainty.high, distribution
            )
        samples.append(out)
    return samples


def _build_samples(
    parameters: ParameterSet,
    sampled_names: Sequence[str],
    n_members: int,
    *,
    seed: int,
    distribution: str,
    sampler: str,
) -> list[dict[str, float]]:
    """The full list of ``n_members`` ``{name: float}`` samples for the chosen sampler.

    ``"mc"`` reproduces the original per-member PRNG draw sequence byte-for-byte (so a
    seed pins the same ensemble as before this sampler split existed); ``"lhs"``/``"sobol"``
    delegate to :func:`_qmc_samples`. Precomputing the whole list is equivalent to drawing
    inside the integration loop — the draw never depends on a member's integration outcome.
    """
    if sampler == "mc":
        rng = np.random.default_rng(seed)
        return [
            sample_parameters(parameters, rng, distribution=distribution, names=sampled_names)
            for _ in range(n_members)
        ]
    return _qmc_samples(
        parameters, sampled_names, n_members, seed=seed, distribution=distribution, sampler=sampler
    )


def _active_reads(process_set: ProcessSet) -> set[str]:
    """Union of parameter names the active Processes and modifiers declare reading."""
    reads: set[str] = set()
    for p in process_set.active:
        reads.update(p.reads)
    for m in process_set.active_modifiers:
        reads.update(m.reads)
    return reads


def _schedule_reads(process_set: ProcessSet, events: Sequence[ScheduledEvent]) -> set[str]:
    """Reads of every Process active at *any* point across a scheduled run.

    ``_active_reads`` alone reflects only the ``t0`` set — but a ``reconfigure`` event
    (a mid-run ``pitch_mlf``) enables Processes that were disabled at compile time, and
    those Processes' ``reads`` drive the *back half* of the run. Sampling only the ``t0``
    reads would under-sample exactly the parameters a pitched schedule's later segments
    depend on. So this replays every event's ``reconfigure`` onto a throwaway copy of the
    enable state and unions the reads at each configuration, then restores the set. Over-
    covering is safe (a sampled parameter no active Process reads is already a documented
    no-op, D-24); under-covering silently narrows the reported spread.
    """
    snapshot = process_set.enabled_snapshot()
    try:
        reads = _active_reads(process_set)
        for e in events:
            if e.reconfigure is not None:
                e.reconfigure(process_set)
                reads |= _active_reads(process_set)
    finally:
        process_set.restore_enabled(snapshot)
    return reads


def _resolve_sample_names(
    process_set: ProcessSet,
    parameters: ParameterSet,
    only: Iterable[str] | None,
    exclude: Iterable[str] | None,
    events: Sequence[ScheduledEvent] = (),
) -> tuple[str, ...]:
    """The sorted, deterministic set of parameter names to sample for this run.

    Defaults to the parameters the active Process set reads across the *whole schedule*
    (sampling anything else is a no-op on the trajectory and only dilutes the member
    count); ``only`` overrides that with an explicit set; ``exclude`` removes names from
    whichever set was chosen (the escape hatch for pinning, e.g. the pKa set that anchors
    the D-18 initial pH). Names not present as provenance parameters are dropped. Sorted
    so the per-member draw order — and thus the seeded reproducibility — does not depend
    on set ordering. With ``events=()`` the schedule union is just the ``t0`` reads, so an
    un-scheduled run samples exactly the names it did before scheduling existed.
    """
    chosen = set(only) if only is not None else _schedule_reads(process_set, events)
    chosen &= set(parameters.names)
    if exclude is not None:
        chosen -= set(exclude)
    return tuple(sorted(chosen))


@dataclass(frozen=True)
class Ensemble:
    """Result of a :func:`simulate_ensemble` run.

    ``t`` is the shared time grid (internal hours); ``nominal`` has shape
    ``(n_vars, n_times)`` and is the deterministic run on the parameters' nominal
    values; ``members`` has shape ``(n_succeeded, n_vars, n_times)`` — the surviving
    sampled trajectories, all on the same grid. ``tier_map`` is the derived
    per-variable confidence tier (a property of provenance, not of the spread, so it is
    the same as a single deterministic run). ``member_params[i]`` is the full resolved
    ``{name: float}`` map that produced ``members[i]`` — needed to audit a member with
    the accounting constants it actually used (e.g. a per-member conservation check must
    read that member's sampled ``biomass_C_fraction``, not the nominal one, or the
    genuine closure reads as drift). Sampling is *not silent*: ``n_requested``,
    ``n_succeeded`` and ``n_failed`` are all reported, so a caller can see how much of
    the ensemble survived and judge whether the spread is trustworthy.

    For a *scheduled* ensemble (``events`` passed to :func:`simulate_ensemble`, D-37),
    ``segment_bounds`` are the shared breakpoint times and ``member_flows[i]`` is member
    ``i``'s own across-jumps external-flow ledger (member-dependent — a ``rack`` removes a
    fraction of the *sampled* lees mass). :meth:`member_trajectory` /
    :meth:`nominal_trajectory` reconstruct full :class:`ScheduledTrajectory` objects from
    them so the ``final == initial + Σ flows`` identity is auditable per draw. All three
    are empty / ``(t0, t_end)`` for an un-scheduled ensemble.
    """

    schema: StateSchema
    t: FloatArray
    nominal: FloatArray
    members: FloatArray
    member_params: tuple[Mapping[str, float], ...]
    tier_map: Mapping[str, Tier]
    sampled_names: tuple[str, ...]
    distribution: str
    sampler: str
    seed: int
    n_requested: int
    n_succeeded: int
    n_failed: int
    failures: tuple[str, ...]
    #: Breakpoint times of the (optional) intervention schedule, ``t_span`` ends included.
    #: Scenario-fixed — the same for every member and the nominal — so stored once. Just
    #: ``(t0, t_end)`` for an un-scheduled ensemble (``events=()``).
    segment_bounds: tuple[float, ...] = ()
    #: Per-member external-flow ledger (``member_flows[i]`` produced ``members[i]``). The
    #: flows are **member-dependent**: a ``rack`` removes a *fraction of the settled lees*,
    #: whose mass at rack time depends on that member's sampled death/growth kinetics, so
    #: each member's removal delta differs. Storing them per member keeps the across-jumps
    #: conservation identity (``final == initial + Σ flows``) auditable for every draw, not
    #: just the nominal. Empty tuples for an un-scheduled ensemble.
    member_flows: tuple[tuple[ExternalFlow, ...], ...] = ()
    #: External-flow ledger of the deterministic nominal run (its own, since the flows
    #: depend on the trajectory). Lets :meth:`nominal_trajectory` audit the nominal identity.
    nominal_flows: tuple[ExternalFlow, ...] = ()

    # -- shape helpers --------------------------------------------------------

    def _rows(self, block: FloatArray, name: str) -> FloatArray:
        """Collapse a ``(..., n_slots, n_times)`` block to the :meth:`series` shape."""
        sl = self.schema.slice(name)
        sub = block[sl, :]
        return sub[0] if sub.shape[0] == 1 else sub

    def series_nominal(self, name: str) -> FloatArray:
        """Deterministic (nominal) time series of one variable — like ``Trajectory.series``."""
        return self._rows(self.nominal, name)

    # -- aggregates -----------------------------------------------------------

    def median(self) -> FloatArray:
        """Element-wise median over members, shape ``(n_vars, n_times)``."""
        return np.median(self.members, axis=0)

    def percentile(self, q: float) -> FloatArray:
        """Element-wise ``q``-th percentile over members, shape ``(n_vars, n_times)``."""
        return np.percentile(self.members, q, axis=0)

    def band(self, name: str, *, low: float = 5.0, high: float = 95.0) -> Band:
        """Low/median/high percentile band plus the nominal series for one variable.

        ``low``/``high`` default to P5/P95 (the outer bracket, per decision D-24). Each
        field follows :meth:`Trajectory.series` shape conventions.
        """
        if not 0.0 <= low <= high <= 100.0:
            raise ValueError(f"require 0 <= low <= high <= 100, got low={low}, high={high}")
        lo_block, hi_block = np.percentile(self.members, [low, high], axis=0)
        med_block = self.median()
        return Band(
            low=self._rows(lo_block, name),
            median=self._rows(med_block, name),
            high=self._rows(hi_block, name),
            nominal=self.series_nominal(name),
        )

    def member_trajectory(self, i: int) -> ScheduledTrajectory:
        """Reconstruct member ``i`` as a :class:`ScheduledTrajectory`.

        Lets a caller reuse the deterministic validation harness (``assert_conserved``,
        ``assert_nonnegative``, …) on any individual sampled member — the per-member
        conservation invariant is exactly how the project guards against a sampled
        parameter set silently breaking a balance. The reconstruction carries the member's
        own ``external_flows`` and the shared ``segment_bounds``, so a *scheduled* member's
        across-jumps identity (``final == initial + Σ flows``) is checkable too, not just
        the constant-total balance of an un-scheduled run (whose flows are empty).
        """
        return ScheduledTrajectory(
            schema=self.schema,
            t=self.t,
            y=self.members[i],
            success=True,
            message="ensemble member",
            tier_map=self.tier_map,
            segment_bounds=self.segment_bounds,
            external_flows=self.member_flows[i] if self.member_flows else (),
        )

    def nominal_trajectory(self) -> ScheduledTrajectory:
        """Reconstruct the deterministic nominal run as a :class:`ScheduledTrajectory`.

        Carries the nominal run's own ``external_flows`` and the ``segment_bounds``, so the
        nominal across-jumps conservation identity is auditable with the same harness as the
        members. For an un-scheduled ensemble the flows are empty and this is the plain
        nominal trajectory.
        """
        return ScheduledTrajectory(
            schema=self.schema,
            t=self.t,
            y=self.nominal,
            success=True,
            message="ensemble nominal",
            tier_map=self.tier_map,
            segment_bounds=self.segment_bounds,
            external_flows=self.nominal_flows,
        )

    @property
    def failure_fraction(self) -> float:
        return self.n_failed / self.n_requested if self.n_requested else 0.0

    def overall_tier(self) -> Tier:
        return min(self.tier_map.values(), default=Tier.VALIDATED)


def simulate_ensemble(
    process_set: ProcessSet,
    parameters: ParameterSet,
    y0: FloatArray,
    t_span: tuple[float, float],
    *,
    n_members: int = 200,
    seed: int = 0,
    t_eval: FloatArray | None = None,
    distribution: str = "triangular",
    sampler: str = "mc",
    only: Iterable[str] | None = None,
    exclude: Iterable[str] | None = None,
    param_tiers: Mapping[str, Tier] | None = None,
    max_failure_fraction: float = 0.5,
    events: Iterable[ScheduledEvent] = (),
    method: str = "BDF",
    rtol: float = 1e-6,
    atol: float = 1e-9,
    max_step: float = np.inf,
) -> Ensemble:
    """Run a Monte-Carlo ensemble of :func:`simulate_scheduled` over sampled parameters.

    Draws ``n_members`` parameter samples from ``parameters``' provenance bands (seeded
    by ``seed``, so the whole ensemble is reproducible), integrates each on the shared
    ``t_eval`` grid, and returns the deterministic nominal run plus the surviving
    members. All members use identical solver settings, ``y0``, and ``events``; only the
    sampled parameters differ.

    ``t_eval`` defaults to 200 evenly spaced points across ``t_span`` (a shared grid is
    required so members can be aggregated). ``param_tiers`` defaults to
    ``parameters.tier_map()`` so the reported tiers are honest (D-1) without the caller
    threading it. See the module docstring for ``distribution``/``sampler``/``only``/
    ``exclude``. ``sampler="mc"`` (default) is i.i.d. Monte Carlo; ``"lhs"``/``"sobol"``
    are low-discrepancy sequences with better per-member (tail) coverage (``"sobol"`` needs
    a power-of-two ``n_members``).

    **Scheduled ensembles (D-37).** ``events`` are timed interventions (a temperature
    ramp, a DAP dose, a mid-run ``pitch_mlf``) handed to :func:`simulate_scheduled`: every
    member is integrated through the *same* schedule, so the spread is the parameter-
    uncertainty band of a *scheduled* scenario. Three consequences follow, all handled
    here: (1) sampling scope is the union of reads across the whole schedule, so a Process
    a ``reconfigure`` enables mid-run has its parameters sampled too (not just the ``t0``
    set); (2) a ``reconfigure`` mutates ``process_set`` in place and is *not* self-
    restoring, so the set is reset to its pre-run enable state before every member (else
    one member's ``pitch`` leaks into the next member's pre-event segments); (3) each
    member's ``external_flows`` are stored per member (they are member-dependent — a
    ``rack`` removes a fraction of the *sampled* lees mass), keeping the across-jumps
    identity auditable draw by draw. With ``events=()`` this is byte-for-byte the previous
    un-scheduled ensemble (``simulate_scheduled`` with no events is a single ``simulate``
    segment). **Caveat:** a parameter that is *both* sampled and overwritten by an event's
    ``param_update`` uses its sampled value pre-event and the fixed compile-time value
    after — no current verb hits this (``temperature_ramp_rate`` declares no ``reads`` and
    is VALIDATED, so it is never sampled), but a future one could, and would need pinning
    via ``exclude``.

    A sampled parameter set can make a member fail — ``solve_ivp`` may report
    ``success=False``, or the right-hand side may *raise* (e.g. the uptake carbon-draw
    guard). Both are caught, recorded in ``failures``, and counted; the samples are drawn
    up front (one per member) so the same seed reproduces the same ensemble including its
    failures. If more than ``max_failure_fraction`` of members fail, this raises rather
    than return a survivorship-biased spread from the lucky survivors.
    """
    if n_members < 1:
        raise ValueError(f"n_members must be >= 1, got {n_members}")
    if distribution not in DISTRIBUTIONS:
        raise ValueError(f"unknown distribution {distribution!r}; expected one of {DISTRIBUTIONS}")
    if sampler not in SAMPLERS:
        raise ValueError(f"unknown sampler {sampler!r}; expected one of {SAMPLERS}")

    events = tuple(events)
    grid = np.linspace(t_span[0], t_span[1], 200) if t_eval is None else np.asarray(t_eval, float)
    tiers = parameters.tier_map() if param_tiers is None else param_tiers
    nominal_values = parameters.resolve()
    sampled_names = _resolve_sample_names(process_set, parameters, only, exclude, events)

    # A reconfigure event mutates process_set in place and does not restore it (D-35); the
    # pristine pre-run enable state is captured once and reset before every run so members
    # cannot leak an enable into one another (the isolation the shared-set ensemble needs).
    pristine = process_set.enabled_snapshot()

    def run(values: Mapping[str, float]) -> ScheduledTrajectory:
        process_set.restore_enabled(pristine)
        return simulate_scheduled(
            process_set,
            values,
            y0,
            t_span,
            events=events,
            param_tiers=tiers,
            t_eval=grid,
            method=method,
            rtol=rtol,
            atol=atol,
            max_step=max_step,
        )

    # Nominal = the deterministic baseline (D-24). Identical solver kwargs + grid + events,
    # so a caller reproduces it byte-for-byte with simulate_scheduled on parameters.resolve()
    # (and, with events=(), with a plain simulate()).
    nominal = run(nominal_values)
    if not nominal.success:
        raise RuntimeError(f"nominal (unsampled) run failed to integrate: {nominal.message}")

    # Draw every member's parameters up front (one row per member, in order). The draw
    # never depends on integration outcome, so this is equivalent to drawing inside the
    # loop — and it makes the same seed reproduce the same ensemble including its failures,
    # for MC and the QMC samplers alike.
    samples = _build_samples(
        parameters, sampled_names, n_members, seed=seed, distribution=distribution, sampler=sampler
    )
    members: list[FloatArray] = []
    member_params: list[Mapping[str, float]] = []
    member_flows: list[tuple[ExternalFlow, ...]] = []
    failures: list[str] = []
    for values in samples:
        try:
            traj = run(values)
        except Exception as exc:  # a sampled RHS can raise (uptake guard, etc.)
            failures.append(f"raised: {exc}")
            continue
        if not traj.success:
            failures.append(f"solver: {traj.message}")
            continue
        if traj.y.shape != nominal.y.shape:  # truncated grid ⇒ unusable for aggregation
            failures.append(f"shape {traj.y.shape} != nominal {nominal.y.shape}")
            continue
        members.append(traj.y)
        member_params.append(values)
        member_flows.append(traj.external_flows)

    # The ensemble is a batch that reset the shared set before every member; leave it in
    # the pristine pre-run state rather than the last member's reconfigured one, so it is
    # side-effect-free on process_set (unlike a single simulate_scheduled run, whose enable
    # is meant to persist — a distinction between a batch and a run).
    process_set.restore_enabled(pristine)

    n_succeeded = len(members)
    n_failed = n_members - n_succeeded
    if n_succeeded == 0:
        raise RuntimeError(
            f"all {n_members} ensemble members failed to integrate; first failures: {failures[:3]}"
        )
    if n_failed / n_members > max_failure_fraction:
        raise RuntimeError(
            f"{n_failed}/{n_members} ensemble members failed "
            f"(> {max_failure_fraction:.0%}); the surviving spread would be "
            f"survivorship-biased. First failures: {failures[:3]}"
        )

    return Ensemble(
        schema=process_set.schema,
        t=nominal.t,
        nominal=nominal.y,
        members=np.stack(members, axis=0),
        member_params=tuple(member_params),
        tier_map=nominal.tier_map,
        sampled_names=sampled_names,
        distribution=distribution,
        sampler=sampler,
        seed=int(seed),
        n_requested=n_members,
        n_succeeded=n_succeeded,
        n_failed=n_failed,
        failures=tuple(failures),
        segment_bounds=nominal.segment_bounds,
        member_flows=tuple(member_flows),
        nominal_flows=nominal.external_flows,
    )
