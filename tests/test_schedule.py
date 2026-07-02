"""Event-scheduled integration driver (decision D-35).

These pin the *verb-agnostic* runtime contract of :func:`simulate_scheduled` with toy
Processes (no winemaking vocabulary — that lives at the compile boundary and is tested
separately): an un-scheduled run is byte-for-byte a plain :func:`simulate`; a
per-segment parameter change (the mechanism a temperature ramp rides on) takes effect at
its breakpoint and a constant-slope forcing integrates *exactly*; a state mutation books
its net injection in the external-flow ledger and appears as a clean, monotone jump; a
mid-run Process-set reconfiguration changes the dynamics from its breakpoint and drags
the reconfigured variable's tier down for the whole run; and out-of-window events are
rejected.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema, VarSpec
from fermentation.core.tiers import Tier
from fermentation.runtime import ScheduledEvent, simulate, simulate_scheduled

# -- toys ---------------------------------------------------------------------


class RampToy(Process):
    """``dA/dt = slope`` — a constant-slope forcing read from params (the shape a
    piecewise-linear temperature schedule drives). Linear, so BDF integrates it exactly."""

    name = "ramp_toy"
    tier = Tier.VALIDATED
    touches = ("A",)
    reads = ("slope",)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        d[schema.slice("A")] = params["slope"]
        return d


class ConstGrowthToy(Process):
    """``dA/dt = +1`` — a parameter-free validated flux (the always-on 'core' in the
    reconfigure test)."""

    name = "const_growth_toy"
    tier = Tier.VALIDATED
    touches = ("A",)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        d[schema.slice("A")] = 1.0
        return d


class SpecFluxToy(Process):
    """``dA/dt = +10`` — a *speculative* flux, disabled until an event enables it."""

    name = "spec_flux_toy"
    tier = Tier.SPECULATIVE
    touches = ("A",)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        d[schema.slice("A")] = 10.0
        return d


def _schema() -> StateSchema:
    return StateSchema([VarSpec("A", "g/L"), VarSpec("B", "g/L", default=0.0)])


def _grid(t0: float = 0.0, t_end: float = 10.0, n: int = 51) -> FloatArray:
    return np.linspace(t0, t_end, n)


# -- isolability: no events == a plain simulate -------------------------------


def test_no_events_is_byte_for_byte_a_plain_simulate():
    schema = _schema()
    ps = ProcessSet(schema, [RampToy()])
    y0 = schema.pack({"A": 20.0})
    params = {"slope": 0.5}
    grid = _grid()

    plain = simulate(ps, params, y0, (0.0, 10.0), t_eval=grid)
    sched = simulate_scheduled(ps, params, y0, (0.0, 10.0), t_eval=grid)

    assert sched.success
    assert np.array_equal(sched.t, plain.t)
    assert np.array_equal(sched.y, plain.y)
    assert sched.external_flows == ()
    assert sched.segment_bounds == (0.0, 10.0)


# -- per-segment parameter update (the temperature-ramp mechanism) ------------


def test_param_update_at_breakpoint_changes_the_slope_and_is_exact():
    schema = _schema()
    ps = ProcessSet(schema, [RampToy()])
    y0 = schema.pack({"A": 20.0})
    grid = _grid()
    # slope 1.0 for [0, 4], then 3.0 for [4, 10]: A(t) piecewise linear, continuous at t=4.
    sched = simulate_scheduled(
        ps,
        {"slope": 1.0},
        y0,
        (0.0, 10.0),
        events=[ScheduledEvent(time_h=4.0, label="steepen", param_update={"slope": 3.0})],
        t_eval=grid,
    )
    assert sched.success
    assert sched.segment_bounds == (0.0, 4.0, 10.0)
    a = sched.series("A")
    t = sched.t
    expected = np.where(t <= 4.0, 20.0 + 1.0 * t, 24.0 + 3.0 * (t - 4.0))
    # dA/dt piecewise-constant + segmenting at the slope change ⇒ integrated exactly.
    assert np.allclose(a, expected, atol=1e-10, rtol=0.0)
    # No state jump ⇒ no external flow, and A is continuous across the breakpoint.
    assert sched.external_flows == ()


def test_breakpoint_time_is_emitted_once_and_grid_is_monotone():
    schema = _schema()
    ps = ProcessSet(schema, [RampToy()])
    y0 = schema.pack({"A": 0.0})
    # t_eval deliberately omits the breakpoint (3.3); the driver must still emit it once.
    sched = simulate_scheduled(
        ps,
        {"slope": 1.0},
        y0,
        (0.0, 10.0),
        events=[ScheduledEvent(time_h=3.3, label="bp", param_update={"slope": 1.0})],
        t_eval=np.linspace(0.0, 10.0, 11),
    )
    assert sched.success
    assert np.all(np.diff(sched.t) > 0.0)  # strictly monotone: no duplicate timestamp
    assert np.count_nonzero(sched.t == 3.3) == 1


# -- state mutation + external-flow ledger ------------------------------------


def test_state_mutation_jumps_the_trajectory_and_books_the_flow():
    schema = _schema()
    ps = ProcessSet(schema, [RampToy()])
    y0 = schema.pack({"A": 10.0})
    grid = _grid(n=101)

    def add_five_to_a(sch: StateSchema, state: FloatArray) -> FloatArray:
        out = state.copy()
        out[sch.slice("A")] += 5.0
        return out

    sched = simulate_scheduled(
        ps,
        {"slope": 0.0},  # flat, so A only moves at the dose
        y0,
        (0.0, 10.0),
        events=[ScheduledEvent(time_h=5.0, label="dose_A", mutate=add_five_to_a)],
        t_eval=grid,
    )
    assert sched.success
    a = sched.series("A")
    t = sched.t
    # Flat before the dose at A0=10, flat after at 15 — a clean +5 step at t=5.
    assert np.allclose(a[t < 5.0], 10.0)
    assert np.allclose(a[t >= 5.0], 15.0)

    # The ledger books exactly the injected +5 on A (and nothing on B).
    assert len(sched.external_flows) == 1
    flow = sched.external_flows[0]
    assert flow.time_h == 5.0 and flow.label == "dose_A"
    assert flow.delta[schema.slice("A")][0] == 5.0
    assert flow.delta[schema.slice("B")][0] == 0.0

    # Run-wide balance: final == initial(from y0) + Σ flows. dA/dt = 0, so all of A's
    # change is the external injection.
    total_injected = sum(f.delta[schema.slice("A")][0] for f in sched.external_flows)
    assert np.isclose(a[-1] - 10.0, total_injected)


# -- mid-run reconfiguration + tier travel ------------------------------------


def _reconfig_ps() -> ProcessSet:
    schema = _schema()
    ps = ProcessSet(schema, [ConstGrowthToy(), SpecFluxToy()])
    ps.disable(SpecFluxToy.name)  # speculative flux off until pitched
    return ps


def _enable_spec_event() -> ScheduledEvent:
    return ScheduledEvent(
        time_h=5.0, label="enable_spec", reconfigure=lambda s: s.enable(SpecFluxToy.name)
    )


def test_reconfiguration_enables_a_process_from_its_breakpoint():
    ps = _reconfig_ps()
    schema = ps.schema
    y0 = schema.pack({"A": 0.0})
    grid = _grid(n=101)
    sched = simulate_scheduled(ps, {}, y0, (0.0, 10.0), events=[_enable_spec_event()], t_eval=grid)
    assert sched.success
    a = sched.series("A")
    t = sched.t
    # Before t=5: only +1/h (A == t). After: +11/h, continuous at t=5 (A(5)=5).
    assert np.allclose(a[t <= 5.0], t[t <= 5.0], atol=1e-9)
    late = t >= 5.0
    assert np.allclose(a[late], 5.0 + 11.0 * (t[late] - 5.0), atol=1e-9)


def test_tier_travels_across_the_run_when_a_speculative_process_is_enabled_late():
    ps = _reconfig_ps()
    schema = ps.schema
    y0 = schema.pack({"A": 0.0})
    grid = _grid(n=21)
    sched = simulate_scheduled(ps, {}, y0, (0.0, 10.0), events=[_enable_spec_event()], t_eval=grid)
    # A is VALIDATED for the first segment but SPECULATIVE for the second; the run-wide
    # reported tier is the least-trustworthy segment (decision D-35).
    assert sched.tier_map["A"] is Tier.SPECULATIVE
    assert sched.overall_tier() is Tier.SPECULATIVE

    # Control: without the reconfiguration, A stays VALIDATED throughout.
    ps2 = _reconfig_ps()
    baseline = simulate_scheduled(ps2, {}, y0, (0.0, 10.0), t_eval=grid)
    assert baseline.tier_map["A"] is Tier.VALIDATED


# -- day-0 and ordering edge cases --------------------------------------------


def test_event_at_t0_seeds_the_run_before_integration():
    schema = _schema()
    ps = ProcessSet(schema, [RampToy()])
    y0 = schema.pack({"A": 0.0})
    grid = _grid()
    # A day-0 param event sets the initial slope; no interior breakpoint is created.
    sched = simulate_scheduled(
        ps,
        {"slope": 0.0},
        y0,
        (0.0, 10.0),
        events=[ScheduledEvent(time_h=0.0, label="initial_slope", param_update={"slope": 2.0})],
        t_eval=grid,
    )
    assert sched.success
    assert sched.segment_bounds == (0.0, 10.0)  # single segment
    assert np.allclose(sched.series("A"), 2.0 * sched.t, atol=1e-10)


def test_same_time_events_apply_in_list_order():
    schema = _schema()
    ps = ProcessSet(schema, [RampToy()])
    y0 = schema.pack({"A": 0.0})
    grid = _grid()
    calls: list[str] = []
    sched = simulate_scheduled(
        ps,
        {"slope": 0.0},
        y0,
        (0.0, 10.0),
        events=[
            ScheduledEvent(time_h=5.0, label="a", reconfigure=lambda s: calls.append("first")),
            ScheduledEvent(time_h=5.0, label="b", reconfigure=lambda s: calls.append("second")),
        ],
        t_eval=grid,
    )
    assert sched.success
    assert calls == ["first", "second"]


def test_events_outside_the_window_raise():
    schema = _schema()
    ps = ProcessSet(schema, [RampToy()])
    y0 = schema.pack({"A": 0.0})
    for bad in (-1.0, 10.0, 11.0):  # before t0, at t_end, after t_end
        with pytest.raises(ValueError, match="outside the integrable window"):
            simulate_scheduled(
                ps,
                {"slope": 1.0},
                y0,
                (0.0, 10.0),
                events=[ScheduledEvent(time_h=bad, label="oops", param_update={"slope": 2.0})],
            )
