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


# -- add_so2: lands on so2_total, rides neither elemental ledger ---------------


def _so2(day: float, so2_mgl: float) -> Intervention:
    return Intervention(day=day, action="add_so2", params={"so2_mgl": so2_mgl})


def test_add_so2_lands_on_so2_total_and_books_one_flow():
    cs = compile_scenario(_wine([_so2(5.0, 40.0)]))
    traj = cs.run()

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "add_so2@5d"
    schema = cs.schema
    # 40 mg/L → 0.040 g/L on so2_total; nothing else moves.
    assert flow.delta[schema.slice("so2_total")][0] == pytest.approx(0.040, rel=1e-12)
    others = np.delete(flow.delta, schema.slice("so2_total").start)
    assert np.count_nonzero(others) == 0


def test_add_so2_perturbs_neither_carbon_nor_nitrogen():
    cs = compile_scenario(_wine([_so2(5.0, 50.0)]))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])

    # SO₂ carries neither element, so the dose contributes nothing to either ledger and both
    # single-run balances still close with no correction term (unlike the DAP nitrogen jump).
    assert all(c_of(f.delta) == pytest.approx(0.0, abs=1e-15) for f in traj.external_flows)
    assert all(n_of(f.delta) == pytest.approx(0.0, abs=1e-15) for f in traj.external_flows)
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0), abs=1e-6)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0), abs=1e-9)


def test_add_so2_raises_the_molecular_so2_readout():
    from fermentation.analysis import molecular_so2_series

    cs = compile_scenario(_wine([_so2(5.0, 60.0)]))
    traj = cs.run()
    params = cs.param_values
    mol = molecular_so2_series(traj.as_trajectory(), params)

    # Before the day-5 dose there is no SO₂ (so2_total ≡ 0 ⇒ molecular ≡ 0); after it, the
    # readout is positive — the dose feeds the free/molecular partition from that time forward.
    before = mol[np.argmin(np.abs(traj.t - 4.5 * 24.0))]
    after = mol[np.argmin(np.abs(traj.t - 6.0 * 24.0))]
    assert before == pytest.approx(0.0, abs=1e-15)
    assert after > 0.0


def test_add_so2_requires_an_so2_slot():
    beer = Scenario(
        name="no-so2",
        medium="beer",
        initial={
            "glucose_gpl": 90.0,
            "maltose_gpl": 140.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 2.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[_so2(3.0, 40.0)],
        duration_days=10.0,
    )
    with pytest.raises(ValueError, match="needs a 'so2_total' slot"):
        compile_scenario(beer)


# -- rack: remove settled lees, leave the wine ---------------------------------


def _rack(day: float, fraction: float) -> Intervention:
    return Intervention(day=day, action="rack", params={"fraction": fraction})


def test_rack_removes_settled_lees_and_leaves_the_wine_untouched():
    # Autolysis opted in so debris (the second lees pool) is non-empty at racking.
    sc = Scenario(
        name="rack-test",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 150.0, "pitch_gpl": 0.25, "autolysis_rate_per_h": 0.002},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[_rack(15.0, 0.9)],
        duration_days=20.0,
    )
    cs = compile_scenario(sc)
    traj = cs.run()
    schema = cs.schema

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "rack@15d"

    # The settled pools drop (negative delta), and 90% of each is removed.
    for name in ("X_dead", "debris"):
        sl = schema.slice(name)
        assert flow.delta[sl][0] < 0.0
    # Everything that stays with the racked-off liquid is untouched (zero delta): viable biomass,
    # sugar, ethanol, YAN, glycerol, byproducts, acids, SO₂.
    for name in ("X", "S", "E", "N", "Gly", "Byp", "esters", "fusels", "tartaric", "so2_total"):
        assert np.all(flow.delta[schema.slice(name)] == 0.0)


def test_rack_fraction_out_of_range_raises():
    with pytest.raises(ValueError, match=r"fraction must be in \[0, 1\]"):
        compile_scenario(_wine([_rack(10.0, 1.5)]))


def test_dap_dose_and_rack_conserve_carbon_and_nitrogen_across_all_jumps():
    # The crown-jewel ledger the D-35 external-flow machinery was built for: a run with an
    # injection (DAP, +N) AND a removal (rack, −C/−N) still satisfies the run-wide identity
    # final == initial + Σ external_flows for BOTH elements, to machine precision.
    sc = Scenario(
        name="dap-and-rack",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 100.0, "pitch_gpl": 0.25, "autolysis_rate_per_h": 0.002},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[_dap(2.0, 0.4), _rack(15.0, 0.8)],
        duration_days=20.0,
    )
    cs = compile_scenario(sc)
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    assert len(traj.external_flows) == 2  # one dose, one rack

    for quantity, tol in ((c_of, 1e-6), (n_of, 1e-9)):
        injected = sum(quantity(f.delta) for f in traj.external_flows)
        initial = quantity(cs.y0)
        final = quantity(traj.y[:, -1])
        assert final == pytest.approx(initial + injected, abs=tol)
    # And the ledger is not trivially zero: the rack removed carbon, the DAP added nitrogen.
    assert sum(c_of(f.delta) for f in traj.external_flows) < 0.0
    assert sum(n_of(f.delta) for f in traj.external_flows) != pytest.approx(0.0, abs=1e-6)


# -- pitch_mlf: mutate the catalyst + reconfigure the Process set --------------


_MLF_NAMES = (
    "malolactic_conversion",
    "malolactic_citrate_metabolism",
    "oenococcus_diacetyl_reduction",
)


def _mlf_wine(interventions: list[Intervention]) -> Scenario:
    return Scenario(
        name="mlf-test",
        medium="wine",
        initial={
            "brix": 22.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.25,
            "malic_gpl": 3.0,
            "initial_ph": 3.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=22.0)],
        interventions=interventions,
        duration_days=25.0,
    )


def _pitch(day: float, pitch_gpl: float = 0.1) -> Intervention:
    return Intervention(day=day, action="pitch_mlf", params={"pitch_gpl": pitch_gpl})


def test_pitch_mlf_mutates_catalyst_and_enables_exactly_the_gated_set():
    cs = compile_scenario(_mlf_wine([_pitch(1.0, 0.1)]))
    # Unpitched at compile: the gated Processes are all disabled.
    assert all(not cs.process_set.is_enabled(name) for name in _MLF_NAMES)

    traj = cs.run()
    # The pitch enabled exactly the gated set (symmetric with the compile-time disable).
    assert all(cs.process_set.is_enabled(name) for name in _MLF_NAMES)

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "pitch_mlf@1d"
    schema = cs.schema
    # The whole jump is the X_mlf catalyst dose; nothing else moves.
    assert flow.delta[schema.slice("X_mlf")][0] == pytest.approx(0.1, rel=1e-12)
    others = np.delete(flow.delta, schema.slice("X_mlf").start)
    assert np.count_nonzero(others) == 0


def test_pitch_mlf_catalyst_perturbs_neither_ledger():
    cs = compile_scenario(_mlf_wine([_pitch(1.0, 0.1)]))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    # X_mlf is an inert carbon-/nitrogen-free catalyst (D-23): the pitch flow is on neither ledger.
    flow = traj.external_flows[0]
    assert c_of(flow.delta) == pytest.approx(0.0, abs=1e-15)
    assert n_of(flow.delta) == pytest.approx(0.0, abs=1e-15)


def test_early_pitch_converts_malic_but_post_af_pitch_stalls():
    # The honest headline (decisions D-23, D-31): pitch timing is now a scenario choice, and only
    # an early/co-inoculation pitch completes — a post-AF pitch lands past the Luong ethanol wall,
    # so the environmental gate keeps malate conversion near zero. EMERGENT from the gate, not
    # imposed: the same Process, only the pitch time differs.
    early = compile_scenario(_mlf_wine([_pitch(1.0)])).run()
    late = compile_scenario(_mlf_wine([_pitch(15.0)])).run()
    unpitched = compile_scenario(_mlf_wine([])).run()

    malic0 = 3.0
    early_malic = early.series("malic")[-1]
    late_malic = late.series("malic")[-1]
    assert unpitched.series("malic")[-1] == pytest.approx(malic0, abs=1e-9)  # inert, unchanged
    assert early_malic < 0.5 * malic0  # early pitch converts most of the malate
    assert late_malic > 0.9 * malic0  # post-AF pitch stalls at the ethanol wall
    assert early_malic < late_malic


def test_mid_run_pitch_drags_malolactic_tiers_for_the_whole_run():
    # tier travels (D-35): the Processes are enabled only from the breakpoint, but min-combining
    # the per-segment tier maps reports the touched slots speculative for the WHOLE trajectory.
    unpitched = compile_scenario(_mlf_wine([])).run()
    pitched = compile_scenario(_mlf_wine([_pitch(1.0)])).run()
    for name in ("malic", "lactic"):
        assert unpitched.tier_map[name] is Tier.VALIDATED  # disabled ⇒ inert ⇒ validated
        assert pitched.tier_map[name] is Tier.SPECULATIVE  # enabled mid-run ⇒ speculative run-wide


def test_pitch_mlf_requires_an_x_mlf_slot():
    beer = Scenario(
        name="no-mlf",
        medium="beer",
        initial={
            "glucose_gpl": 90.0,
            "maltose_gpl": 140.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 2.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[_pitch(3.0)],
        duration_days=10.0,
    )
    with pytest.raises(ValueError, match="needs an 'X_mlf' slot"):
        compile_scenario(beer)


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
