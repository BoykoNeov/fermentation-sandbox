"""Event-scheduled integration — discrete interventions and continuous forcing.

:func:`simulate` integrates a *fixed* Process set over a *fixed* parameter map. Real
ferments are not fixed: a winemaker doses DAP on day 2, adds SO₂ at pressing, racks
off the lees, pitches malolactic bacteria after dryness; and the cellar temperature
ramps between set-points rather than sitting at one value. Both are *time-driven*
departures from a single continuous ODE, and both need the same machinery: integrate
a segment, stop at a breakpoint, change something, restart.

This module is that machinery. :func:`simulate_scheduled` walks the run as a sequence
of segments separated by :class:`ScheduledEvent` breakpoints; between breakpoints it
calls :func:`simulate` unchanged (so the pure, deterministic core is never touched),
and at each breakpoint it can

* **mutate the state** — a dose adds mass to a slot, racking removes it (a jump in
  ``y``); and/or
* **reconfigure the Process set** — enable a Process that was dosed in (e.g. pitch
  malolactic bacteria mid-run), in place; and/or
* **update parameters** — change a value that takes effect from that time forward
  (this is how a piecewise-linear temperature schedule is driven: each segment gets a
  constant ``temperature_ramp_rate`` slope).

**Verb-agnostic by construction.** The driver knows nothing about DAP, SO₂, MLF, or
temperature. It takes opaque ``(time, mutate, reconfigure, param_update)`` events. The
*vocabulary* — which winemaking verb maps to which state mutation, which unit
conversion, which Processes to enable — lives at the scenario→core compile boundary
(``fermentation.scenario.compile``), exactly where the industry-unit conversion and the
Process disable-gates already sit (decision D-3). Runtime drives time; the boundary owns
vocabulary and units; the core stays pure physics.

**Why segment-and-restart, not ``solve_ivp(events=...)``.** SciPy's event mechanism
detects zero-crossings and can *terminate* integration, but it cannot mutate the state
and resume. A discrete dose is a genuine discontinuity in the trajectory, so the only
correct approach is to stop at it, apply the jump, and start a fresh ``solve_ivp`` from
the mutated state. The implicit BDF solver re-initialises its order at each restart;
that is correct behaviour at a discontinuity, not a performance bug.

**Conservation across a discontinuity (a prime directive, made a first-class output).**
A dose injects mass *from outside* the modelled system, so the single-run invariant
"total carbon at the end equals total carbon at the start" no longer holds — it becomes
``final == initial + Σ(external inputs) − Σ(external outputs)``. Each mutation's net
injection is recorded as an :class:`ExternalFlow` (the post-minus-pre state delta), so
the run-wide balance is auditable: the continuous ODE still closes *exactly within every
segment*, and the ledger is the correction term across the jumps. The driver books the
raw state delta and stays verb-agnostic; weighting it by carbon/nitrogen is the existing
``conservation`` helpers' job, so no per-verb chemistry lives here.

**Tier travels across reconfiguration.** When a mid-run pitch enables a speculative
Process for the back half of the run, the variables that Process touches must report
speculative for the *whole* concatenated trajectory — a run is only as trustworthy as
its least-trustworthy segment. So the per-segment ``tier_map`` is snapshotted and
:func:`~fermentation.core.tiers.combine`d (min) across segments (decision D-35).

**Isolability.** With no events and a single segment, :func:`simulate_scheduled` calls
:func:`simulate` once with the identical arguments, so its trajectory is byte-for-byte a
plain :func:`simulate` run — the scheduling layer adds nothing to an un-scheduled run
(the same discipline the stochastic ensemble wrapper and every speculative Process
follow).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field

import numpy as np

from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier, combine
from fermentation.runtime.integrate import Trajectory, simulate

#: A state mutation applied at a breakpoint: ``(schema, state) -> new_state``. Must
#: return a fresh array of the same shape (the driver books ``new − old`` as the
#: intervention's external flow), and must not mutate its argument in place.
StateMutation = Callable[[StateSchema, FloatArray], FloatArray]

#: An in-place reconfiguration of the Process set at a breakpoint, e.g.
#: ``lambda ps: ps.enable("malolactic_conversion")``. Applied to the *same* set that
#: subsequent segments integrate, so an enable persists for the rest of the run.
Reconfiguration = Callable[[ProcessSet], None]


@dataclass(frozen=True)
class ScheduledEvent:
    """One timed intervention: a breakpoint plus what changes there.

    Any combination of the three effects may be present (all optional); a
    param-update-only event (a temperature slope change) carries no ``mutate`` or
    ``reconfigure``. ``time_h`` is in canonical internal hours. Events at ``t_span[0]``
    are applied to the initial state/params *before* the first segment; interior events
    (``t0 < time_h < t_end``) each open a new segment. An event at or after ``t_span[1]``
    has nothing left to integrate and is rejected by the driver (the compile boundary
    decides whether a late scenario intervention is an error or a no-op).
    """

    time_h: float
    label: str = ""
    #: Parameters that take effect from ``time_h`` forward (merged into the running map).
    param_update: Mapping[str, float] = field(default_factory=dict)
    #: State jump applied at ``time_h`` (dose / racking); ``None`` ⇒ state is continuous.
    mutate: StateMutation | None = None
    #: In-place Process-set change at ``time_h`` (e.g. enable a pitched Process).
    reconfigure: Reconfiguration | None = None


@dataclass(frozen=True)
class ExternalFlow:
    """One intervention's net injection (+) or removal (−) of state mass.

    ``delta`` is the post-mutation minus pre-mutation state vector, so a DAP dose shows
    a positive entry on ``N`` and racking shows negative entries on the lees pools. The
    continuous ODE closes exactly within each segment, so the run-wide conservation
    identity is ``final == initial + Σ delta`` — weighted by
    ``fermentation.validation``'s ``total_carbon`` / ``total_nitrogen`` for the elemental
    balances. Booking the raw state delta keeps the ledger verb-agnostic (no per-dose
    chemistry here).
    """

    time_h: float
    label: str
    delta: FloatArray


@dataclass(frozen=True)
class ScheduledTrajectory:
    """Result of a :func:`simulate_scheduled` run.

    Mirrors :class:`~fermentation.runtime.integrate.Trajectory` (same ``t``/``y`` shapes
    and ``series``/``final``/``overall_tier`` helpers) so downstream code can treat a
    scheduled run like a plain one, and adds the two things scheduling introduces:
    ``external_flows`` (the conservation-across-jumps ledger) and ``segment_bounds`` (the
    breakpoint times used, including ``t_span`` ends). ``tier_map`` is min-combined across
    segments, so a Process enabled only for part of the run still drags its variables'
    reported tier down for the whole trajectory (decision D-35).
    """

    schema: StateSchema
    t: FloatArray
    y: FloatArray
    success: bool
    message: str
    tier_map: Mapping[str, Tier]
    segment_bounds: tuple[float, ...]
    external_flows: tuple[ExternalFlow, ...]

    def series(self, name: str) -> FloatArray:
        """Time series of one (scalar or vector) variable — like ``Trajectory.series``."""
        sl = self.schema.slice(name)
        block = self.y[sl, :]
        return block[0] if block.shape[0] == 1 else block

    def final(self) -> dict[str, float | FloatArray]:
        """The last state as a name -> value(s) mapping."""
        return self.schema.unpack(self.y[:, -1])

    def overall_tier(self) -> Tier:
        return min(self.tier_map.values(), default=Tier.VALIDATED)

    def as_trajectory(self) -> Trajectory:
        """Drop the scheduling extras, exposing the run as a plain :class:`Trajectory`."""
        return Trajectory(
            schema=self.schema,
            t=self.t,
            y=self.y,
            success=self.success,
            message=self.message,
            tier_map=self.tier_map,
        )


def simulate_scheduled(
    process_set: ProcessSet,
    params: Mapping[str, float],
    y0: FloatArray,
    t_span: tuple[float, float],
    *,
    events: Iterable[ScheduledEvent] = (),
    param_tiers: Mapping[str, Tier] | None = None,
    t_eval: FloatArray | None = None,
    method: str = "BDF",
    rtol: float = 1e-6,
    atol: float = 1e-9,
    max_step: float = np.inf,
) -> ScheduledTrajectory:
    """Integrate through timed interventions by segmenting and restarting :func:`simulate`.

    ``events`` are sorted by time (stably, so ties keep their input order — a
    deterministic same-instant ordering). Events at ``t_span[0]`` are applied to the
    initial state/params before the first segment; each interior event opens a new
    segment integrated with the mutated state, reconfigured Process set, and updated
    parameters in force. The solver settings mirror :func:`simulate`.

    Returns a :class:`ScheduledTrajectory` on the concatenated grid (breakpoint times are
    included and emitted **post**-mutation, so a dose appears as a clean jump and the time
    axis stays strictly monotone — no duplicate timestamps). With ``events=()`` this is a
    single :func:`simulate` call with identical arguments, hence byte-for-byte a plain run.

    Raises ``ValueError`` for a bad ``y0`` shape, a non-positive span, or an event outside
    ``[t_span[0], t_span[1])``.
    """
    schema = process_set.schema
    t0, t_end = float(t_span[0]), float(t_span[1])
    if t_end <= t0:
        raise ValueError(f"t_span must be increasing, got {t_span}")
    y0 = np.asarray(y0, dtype=np.float64)
    if y0.shape != (schema.size,):
        raise ValueError(f"y0 has shape {y0.shape}, expected ({schema.size},)")

    ordered = sorted(events, key=lambda e: e.time_h)  # stable ⇒ input order breaks ties
    for e in ordered:
        if e.time_h < t0 or e.time_h >= t_end:
            raise ValueError(
                f"event {e.label!r} at t={e.time_h} is outside the integrable window "
                f"[{t0}, {t_end}); interventions must fall within the run"
            )

    current_params: dict[str, float] = dict(params)
    current_y = np.array(y0, dtype=np.float64, copy=True)
    flows: list[ExternalFlow] = []

    def apply(event: ScheduledEvent) -> None:
        nonlocal current_y, current_params
        if event.mutate is not None:
            new_y = np.asarray(event.mutate(schema, current_y), dtype=np.float64)
            if new_y.shape != current_y.shape:
                raise ValueError(
                    f"mutation for event {event.label!r} returned shape {new_y.shape}, "
                    f"expected ({schema.size},)"
                )
            flows.append(ExternalFlow(event.time_h, event.label, new_y - current_y))
            current_y = new_y
        if event.reconfigure is not None:
            event.reconfigure(process_set)
        if event.param_update:
            current_params = {**current_params, **event.param_update}

    # Events exactly at t0 seed the run (the day-0 temperature slope, a pre-fermentation
    # addition) — applied before any integration so segment 0 starts from them.
    for e in ordered:
        if e.time_h == t0:
            apply(e)
    interior = [e for e in ordered if t0 < e.time_h < t_end]
    breakpoints = sorted({e.time_h for e in interior})
    bounds = [t0, *breakpoints, t_end]

    # Shared output grid, augmented with the breakpoint times so every intervention is
    # visible in the result. np.unique sorts and de-duplicates.
    grid = np.linspace(t0, t_end, 200) if t_eval is None else np.asarray(t_eval, dtype=np.float64)
    grid = np.unique(np.concatenate([grid, np.asarray(bounds, dtype=np.float64)]))
    grid = grid[(grid >= t0) & (grid <= t_end)]

    t_parts: list[FloatArray] = []
    y_parts: list[FloatArray] = []
    seg_tier_maps: list[Mapping[str, Tier]] = []

    n_seg = len(bounds) - 1
    for i in range(n_seg):
        a, b = bounds[i], bounds[i + 1]
        is_final = i == n_seg - 1
        # Evaluate at the grid points in [a, b]; force both endpoints present so the
        # segment's terminal state (at b) is available to seed the next segment / mutation.
        pts = grid[(grid >= a) & (grid <= b)]
        pts = np.unique(np.concatenate([pts, np.asarray([a, b], dtype=np.float64)]))
        seg = simulate(
            process_set,
            current_params,
            current_y,
            (a, b),
            param_tiers=param_tiers,
            t_eval=pts,
            method=method,
            rtol=rtol,
            atol=atol,
            max_step=max_step,
        )
        seg_tier_maps.append(process_set.tier_map(param_tiers))
        if not seg.success:
            t_done = np.concatenate(t_parts) if t_parts else np.asarray([], dtype=np.float64)
            y_done = (
                np.concatenate(y_parts, axis=1)
                if y_parts
                else np.empty((schema.size, 0), dtype=np.float64)
            )
            return ScheduledTrajectory(
                schema=schema,
                t=t_done,
                y=y_done,
                success=False,
                message=f"segment [{a}, {b}] failed to integrate: {seg.message}",
                tier_map={n: combine([m[n] for m in seg_tier_maps]) for n in schema.names},
                segment_bounds=tuple(bounds),
                external_flows=tuple(flows),
            )

        # The terminal state (at b) seeds the next segment. Emit [a, b) for non-final
        # segments so b is emitted once, post-mutation, as the next segment's first point;
        # the final segment emits its closed endpoint t_end.
        current_y = seg.y[:, -1]
        keep = (seg.t >= a) & (seg.t <= b) if is_final else (seg.t >= a) & (seg.t < b)
        t_parts.append(seg.t[keep])
        y_parts.append(seg.y[:, keep])

        if not is_final:
            for e in interior:
                if e.time_h == b:
                    apply(e)

    return ScheduledTrajectory(
        schema=schema,
        t=np.concatenate(t_parts),
        y=np.concatenate(y_parts, axis=1),
        success=True,
        message="ok",
        tier_map={n: combine([m[n] for m in seg_tier_maps]) for n in schema.names},
        segment_bounds=tuple(bounds),
        external_flows=tuple(flows),
    )
