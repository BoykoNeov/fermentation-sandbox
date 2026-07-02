"""Discrete winemaking interventions — the verb registry at the compile boundary (D-36).

``Scenario.interventions`` was declared since Milestone 1 but never consumed: the compile
seam turned only the temperature schedule into events. This activates the timeline of
winemaking verbs. Each verb compiles at the vocabulary boundary (decision D-3) into an
opaque :class:`~fermentation.runtime.schedule.ScheduledEvent`, and the verb-agnostic driver
(``simulate_scheduled``, D-35) segments-and-restarts around it, booking each state jump as an
:class:`~fermentation.runtime.schedule.ExternalFlow` for the conservation ledger.

This module pins the first verb, ``add_dap``:

* the dose lands exactly on ``N`` and is booked as one external flow (nothing else moves);
* nitrogen conserves across the jump — ``final == initial + Σ flows`` (the crown-jewel ledger
  the D-35 driver built the ledger for), and carbon is untouched;
* the emergent H₂S headline the *static* D-29 lever could not produce — a mid-ferment dose
  momentarily closes the inverse H₂S gate while sugar (hence flux) is still present, so the
  production rate drops right after the dose and net cumulative H₂S falls;
* the vocabulary discipline: unknown verbs, out-of-window days, and bad params raise loudly;
* isolability — a scenario with no interventions (and no ramp) is byte-for-byte a plain run.
"""

import numpy as np
import pytest

from fermentation.core.tiers import Tier
from fermentation.runtime import simulate
from fermentation.runtime.schedule import ScheduledTrajectory
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.validation.conservation import total_carbon, total_nitrogen

_DAP_N_FRACTION = 0.2121  # exact (NH4)2HPO4 stoichiometry (additions.yaml)


def _wine(
    interventions: list[Intervention],
    *,
    yan_mgl: float = 80.0,
    days: float = 14.0,
    celsius: float = 20.0,
) -> Scenario:
    return Scenario(
        name="dap-test",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        interventions=interventions,
        duration_days=days,
    )


def _dap(day: float, dap_gpl: float) -> Intervention:
    return Intervention(day=day, action="add_dap", params={"dap_gpl": dap_gpl})


# -- the dose lands on N and is booked as one external flow --------------------


def test_add_dap_lands_on_nitrogen_and_books_one_flow():
    cs = compile_scenario(_wine([_dap(2.0, 0.4)]))
    traj = cs.run()

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "add_dap@2d"
    assert flow.time_h == pytest.approx(48.0)  # 2 days in canonical hours

    schema = cs.schema
    # The whole injection is on N (= dap_gpl * fraction); every other slot delta is zero.
    expected_n = 0.4 * _DAP_N_FRACTION
    assert flow.delta[schema.slice("N")][0] == pytest.approx(expected_n, rel=1e-12)
    others = np.delete(flow.delta, schema.slice("N").start)
    assert np.count_nonzero(others) == 0


def test_dap_nitrogen_fraction_is_validated_and_exact():
    cs = compile_scenario(_wine([_dap(2.0, 0.4)]))
    frac = cs.parameters["dap_nitrogen_fraction"]
    assert frac.tier is Tier.VALIDATED
    assert frac.value == pytest.approx(0.2121)
    # exact stoichiometry ⇒ zero-width band, never swept by the ensemble
    assert frac.uncertainty.low == frac.value == frac.uncertainty.high


# -- conservation across the jump: final == initial + Σ flows (crown jewel) ----


def test_add_dap_conserves_nitrogen_across_the_jump():
    cs = compile_scenario(_wine([_dap(3.0, 0.5)]))
    traj = cs.run()
    schema = cs.schema
    f_n = cs.param_values["biomass_N_fraction"]  # the value growth actually uses (D-14)
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)

    injected = sum(n_of(flow.delta) for flow in traj.external_flows)
    assert injected == pytest.approx(0.5 * _DAP_N_FRACTION, rel=1e-12)

    # The continuous ODE closes within each segment, so the run-wide identity is
    # final == initial + Σ(external injections). Nitrogen is conserved to machine precision.
    n_initial = n_of(cs.y0)
    n_final = n_of(traj.y[:, -1])
    assert n_final == pytest.approx(n_initial + injected, abs=1e-9)


def test_add_dap_leaves_carbon_untouched():
    cs = compile_scenario(_wine([_dap(3.0, 0.5)]))
    traj = cs.run()
    schema = cs.schema
    f_c = cs.param_values["biomass_C_fraction"]
    c_of = total_carbon(schema, biomass_carbon_fraction=f_c)

    # DAP injects nitrogen only (phosphate dropped, no carbon), so every flow is carbon-free
    # and the single-run carbon balance still closes with no ledger correction term.
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0), abs=1e-6)


# -- the emergent H2S headline (a timing effect a static dose cannot produce) --


def _run_pair() -> tuple[ScheduledTrajectory, ScheduledTrajectory]:
    grid = np.linspace(0.0, 14.0 * 24.0, 1400)
    undosed = compile_scenario(_wine([])).run(t_eval=grid)
    dosed = compile_scenario(_wine([_dap(2.0, 0.4)])).run(t_eval=grid)
    return undosed, dosed


def _rate_at(rate: np.ndarray, t_axis: np.ndarray, day: float) -> float:
    return float(rate[np.argmin(np.abs(t_axis - day * 24.0))])


def test_mid_ferment_dap_dose_drops_the_h2s_production_rate():
    undosed, dosed = _run_pair()
    # The dosed run has an extra grid point (the inserted day-2 breakpoint), so each run's
    # production rate is differenced against its OWN time axis.
    rate_u = np.gradient(undosed.series("h2s"), undosed.t)
    rate_d = np.gradient(dosed.series("h2s"), dosed.t)

    # Just before the day-2 dose the two runs are identical (same N history); just after, the
    # restored N re-represses the inverse gate K_h2s_n/(K_h2s_n+N) while sugar (the flux the
    # gate multiplies) is still present ⇒ the dosed production rate falls below the undosed one.
    assert _rate_at(rate_d, dosed.t, 1.8) == pytest.approx(
        _rate_at(rate_u, undosed.t, 1.8), rel=1e-6
    )  # identical pre-dose
    assert _rate_at(rate_d, dosed.t, 2.1) < 0.75 * _rate_at(rate_u, undosed.t, 2.1)  # gate closes


def test_dap_dose_lowers_net_cumulative_h2s():
    undosed, dosed = _run_pair()
    # Net over the whole ferment the gate closure dominates the competing extra-biomass flux
    # (more N ⇒ more growth ⇒ more flux later), so a DAP addition reduces total H₂S — the
    # realistic direction (DAP is the standard H₂S-management lever). This is EMERGENT, not
    # imposed: the model has no explicit "DAP lowers H₂S" term (decision D-36 / D-29).
    assert dosed.series("h2s")[-1] < undosed.series("h2s")[-1]


# -- vocabulary discipline: loud failures --------------------------------------


def test_unknown_action_raises():
    with pytest.raises(ValueError, match="unknown intervention action 'add_potion'"):
        compile_scenario(_wine([Intervention(day=1.0, action="add_potion")]))


def test_intervention_at_or_after_duration_raises():
    with pytest.raises(ValueError, match="at or beyond the run duration"):
        compile_scenario(_wine([_dap(14.0, 0.4)], days=14.0))


def test_add_dap_missing_required_param_raises():
    with pytest.raises(ValueError, match="missing required param 'dap_gpl'"):
        compile_scenario(_wine([Intervention(day=2.0, action="add_dap")]))


def test_add_dap_unknown_param_raises():
    bad = Intervention(day=2.0, action="add_dap", params={"dap_gpl": 0.4, "phosphate_gpl": 0.1})
    with pytest.raises(ValueError, match="unknown param"):
        compile_scenario(_wine([bad]))


def test_add_dap_negative_dose_raises():
    with pytest.raises(ValueError, match="must be >= 0"):
        compile_scenario(_wine([_dap(2.0, -0.4)]))


# -- isolability: no interventions ⇒ byte-for-byte a plain run -----------------


def test_no_interventions_is_byte_for_byte_plain_simulate():
    # A wine scenario with no interventions and a flat (single-knot) schedule compiles to an
    # empty event tuple, so CompiledScenario.run() is a single simulate() call with identical
    # arguments — the scheduling layer adds nothing (the D-35 isolability discipline).
    cs = compile_scenario(_wine([]))
    assert cs.events == ()

    grid = np.linspace(0.0, cs.t_span_h[1], 200)
    scheduled = cs.run(t_eval=grid)
    plain = simulate(
        cs.process_set,
        cs.param_values,
        cs.y0,
        cs.t_span_h,
        param_tiers=cs.parameters.tier_map(),
        t_eval=grid,
    )
    assert np.array_equal(scheduled.y, plain.y)
    assert scheduled.external_flows == ()
