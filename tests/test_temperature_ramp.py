"""Temperature-schedule ramp — driven ``T`` state + the compile→event-loop path (D-35).

The temperature schedule was declared in the Scenario since Milestone 1 but inert (only
its earliest knot seeded the initial ``T``; ``T`` then sat constant). This activates it:
``TemperatureRamp`` drives ``dT/dt`` along the piecewise-linear schedule, the compile
boundary turns the knots into slope-change events for :func:`simulate_scheduled`, and the
Arrhenius kinetics see the true ramped temperature. These tests pin: the Process is an
isothermal no-op by default (byte-for-byte); a straight ramp integrates ``T`` *exactly*
and is segmented only where the slope changes; the rate is a VALIDATED, un-sampled
scenario input; and the ramp genuinely feeds the kinetics.
"""

import numpy as np
import pytest

from fermentation.core.kinetics import TemperatureRamp
from fermentation.core.kinetics.temperature import RAMP_RATE
from fermentation.core.media import get_medium
from fermentation.core.tiers import Tier
from fermentation.runtime import simulate, simulate_scheduled
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario.compile import _temperature_ramp_schedule
from fermentation.units.convert import celsius_to_kelvin, days_to_hours

# -- Process metadata + isothermal no-op --------------------------------------


def test_process_metadata_and_isothermal_default():
    p = TemperatureRamp()
    assert p.name == "temperature_ramp"
    assert p.tier is Tier.VALIDATED
    assert p.touches == ("T",)
    assert p.reads == ()  # scenario-exact forcing, deliberately not a D-1 read

    schema = get_medium("wine").schema
    y = schema.pack({"X": 0.25, "S": [200.0], "E": 0.0, "N": 0.1, "T": 293.15, "CO2": 0.0})
    # No rate in params ⇒ dT/dt = 0 (the isothermal no-op that keeps un-ramped runs exact).
    d = p.derivatives(0.0, y, schema, {})
    assert d[schema.slice("T")][0] == 0.0
    assert np.count_nonzero(d) == 0
    # With a rate, it writes exactly that slope onto T and nothing else.
    d = p.derivatives(0.0, y, schema, {RAMP_RATE: 0.05})
    assert d[schema.slice("T")][0] == 0.05
    assert np.count_nonzero(d) == 1


# -- compile: schedule → slope + events ---------------------------------------


def _wine(
    schedule: list[TemperaturePoint],
    *,
    days: float = 6.0,
    initial: dict[str, float] | None = None,
) -> Scenario:
    return Scenario(
        name="ramp-test",
        medium="wine",
        initial=initial or {"brix": 24.0, "yan_mgl": 150.0, "pitch_gpl": 0.25},
        temperature_schedule=schedule,
        duration_days=days,
    )


def test_single_knot_schedule_is_isothermal_no_events():
    cs = compile_scenario(_wine([TemperaturePoint(day=0.0, celsius=20.0)]))
    assert cs.events == ()
    assert RAMP_RATE not in cs.parameters  # no ramp ⇒ no injected rate parameter


def test_flat_multi_knot_schedule_is_isothermal_no_events():
    # Same temperature at every knot ⇒ zero slope everywhere ⇒ no ramp.
    cs = compile_scenario(
        _wine([TemperaturePoint(day=0.0, celsius=20.0), TemperaturePoint(day=6.0, celsius=20.0)])
    )
    assert cs.events == ()
    assert RAMP_RATE not in cs.parameters


def test_straight_ramp_has_initial_slope_and_no_interior_events():
    schedule = [TemperaturePoint(day=0.0, celsius=20.0), TemperaturePoint(day=6.0, celsius=26.0)]
    t_end = days_to_hours(6.0)
    slope, events = _temperature_ramp_schedule(_wine(schedule), t_end)
    expected = (celsius_to_kelvin(26.0) - celsius_to_kelvin(20.0)) / t_end
    assert slope == pytest.approx(expected)
    assert events == ()  # a single slope ⇒ one segment, no interior breakpoint


def test_collinear_knots_do_not_open_a_segment():
    # 20 → 23 → 26 at days 0/3/6 is one straight line: identical slope on both intervals.
    schedule = [
        TemperaturePoint(day=0.0, celsius=20.0),
        TemperaturePoint(day=3.0, celsius=23.0),
        TemperaturePoint(day=6.0, celsius=26.0),
    ]
    slope, events = _temperature_ramp_schedule(_wine(schedule), days_to_hours(6.0))
    assert slope > 0.0
    assert events == ()  # collinear ⇒ no slope change ⇒ no interior segment


def test_slope_change_opens_exactly_one_event():
    # Ramp 20 → 26 over days 0-3, then hold at 26 to day 6: slope changes to 0 at day 3.
    schedule = [
        TemperaturePoint(day=0.0, celsius=20.0),
        TemperaturePoint(day=3.0, celsius=26.0),
        TemperaturePoint(day=6.0, celsius=26.0),
    ]
    slope, events = _temperature_ramp_schedule(_wine(schedule), days_to_hours(6.0))
    assert slope > 0.0
    assert len(events) == 1
    assert events[0].time_h == pytest.approx(days_to_hours(3.0))
    assert events[0].param_update[RAMP_RATE] == pytest.approx(0.0)


def test_hold_before_first_and_after_last_knot():
    # First knot at day 2, last at day 4, run to day 6: slope 0 on [0,2], ramp on [2,4],
    # slope 0 on [4,6]. Two slope changes ⇒ two interior events; initial slope 0.
    schedule = [TemperaturePoint(day=2.0, celsius=18.0), TemperaturePoint(day=4.0, celsius=24.0)]
    slope, events = _temperature_ramp_schedule(_wine(schedule, days=6.0), days_to_hours(6.0))
    assert slope == 0.0  # held before the first knot
    times = [e.time_h for e in events]
    assert times == pytest.approx([days_to_hours(2.0), days_to_hours(4.0)])
    assert events[0].param_update[RAMP_RATE] > 0.0  # ramp starts at the first knot
    assert events[1].param_update[RAMP_RATE] == pytest.approx(0.0)  # held after the last


# -- exact integration + isolability through the runtime ----------------------


def test_ramp_integrates_temperature_exactly_to_the_analytic_line():
    schedule = [TemperaturePoint(day=0.0, celsius=20.0), TemperaturePoint(day=6.0, celsius=26.0)]
    cs = compile_scenario(_wine(schedule))
    grid = np.linspace(*cs.t_span_h, 200)
    traj = simulate_scheduled(
        cs.process_set, cs.param_values, cs.y0, cs.t_span_h, events=cs.events, t_eval=grid
    )
    assert traj.success
    slope = (celsius_to_kelvin(26.0) - celsius_to_kelvin(20.0)) / cs.t_span_h[1]
    expected = celsius_to_kelvin(20.0) + slope * traj.t
    # dT/dt constant + segmenting at slope changes ⇒ BDF integrates the line to round-off.
    assert np.allclose(traj.series("T"), expected, atol=1e-10, rtol=0.0)


def test_isothermal_scenario_scheduled_equals_plain_simulate():
    # No events ⇒ simulate_scheduled is a single simulate call with identical args.
    cs = compile_scenario(_wine([TemperaturePoint(day=0.0, celsius=20.0)]))
    grid = np.linspace(*cs.t_span_h, 150)
    plain = simulate(cs.process_set, cs.param_values, cs.y0, cs.t_span_h, t_eval=grid)
    sched = simulate_scheduled(
        cs.process_set, cs.param_values, cs.y0, cs.t_span_h, events=cs.events, t_eval=grid
    )
    assert np.array_equal(sched.t, plain.t)
    assert np.array_equal(sched.y, plain.y)


# -- tier + provenance of the rate --------------------------------------------


def test_ramp_rate_is_validated_and_unsampled_and_T_stays_validated():
    schedule = [TemperaturePoint(day=0.0, celsius=20.0), TemperaturePoint(day=6.0, celsius=26.0)]
    cs = compile_scenario(_wine(schedule))
    rate = cs.parameters[RAMP_RATE]
    assert rate.tier is Tier.VALIDATED
    assert rate.uncertainty.low == rate.uncertainty.high  # zero-width ⇒ never sampled
    # T is driven by a VALIDATED Process reading a VALIDATED (undeclared) forcing ⇒ VALIDATED.
    tiers = cs.parameters.tier_map()
    assert cs.process_set.tier_of("T", tiers) is Tier.VALIDATED


# -- the ramp actually feeds the kinetics -------------------------------------


def test_warming_ramp_ferments_between_the_two_isothermal_bounds():
    # A run that ramps 14 → 30 °C must finish with less residual sugar than a cold-held run
    # (warming speeds it) and more than a hot-held run (it spent its early biomass-building
    # phase cold) — proof that Arrhenius reads the true time-varying T, not a constant.
    days = 5.0

    def final_sugar(schedule: list[TemperaturePoint]) -> float:
        cs = compile_scenario(_wine(schedule, days=days))
        traj = simulate_scheduled(
            cs.process_set, cs.param_values, cs.y0, cs.t_span_h, events=cs.events
        )
        assert traj.success
        return float(traj.final()["S"])

    cold = final_sugar([TemperaturePoint(day=0.0, celsius=14.0)])
    hot = final_sugar([TemperaturePoint(day=0.0, celsius=30.0)])
    ramp = final_sugar(
        [TemperaturePoint(day=0.0, celsius=14.0), TemperaturePoint(day=days, celsius=30.0)]
    )
    assert hot < ramp < cold
