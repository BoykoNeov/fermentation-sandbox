"""Scenario-level aging-phase wiring — the ``begin_aging`` verb (decision D-70).

D-69 built :class:`~fermentation.core.kinetics.aging.EsterHydrolysis` (the first §4.1 aging
Process) and exercised it via a hand-built ``ProcessSet``. D-70 wires it into the *scenario*
pipeline end-to-end:

* ``EsterHydrolysis`` is wired into both media but DISABLED at the compile seam (aging is
  inherently post-ferment — there is no aging at t0), so an un-aged scenario is byte-for-byte
  the pre-aging core (prime directive #3);
* a ``begin_aging`` intervention — the ``pitch_mlf`` reconfigure pattern MINUS the state
  mutation — re-enables it from its breakpoint over a post-fermentation aging segment whose
  span is expressed by ``duration_days``;
* ``aging.yaml`` rides in every compiled scenario's ``ParameterSet`` (``compile.py``'s
  ``shared_files``) but is read by nothing until ``begin_aging`` fires;
* the §7 slow-phase integration comes for free from ``simulate_scheduled``'s segment restart:
  with the fermentative flux gone at dryness the solver takes large steps across the quiescent
  aging segment, and every other producer of ``esters``/``fusels``/``Byp`` is flux-gated and
  silent there — so the aging signal is unconfounded (Stance A).

These tests pin: the compile-seam enable/disable gate; the emergent aging headline (esters fade,
fusels + Byp rise) end-to-end vs an otherwise-identical un-aged run; end-to-end carbon closure
(a pure inter-pool transfer — no external flow, since ``begin_aging`` mutates no state); the
speculative tier floor min-combined across the whole run; and the loud errors (params, unknown
verb-params) at the vocabulary boundary.
"""

import numpy as np
import pytest

from fermentation.core.kinetics.aging import EsterHydrolysis
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.validation.conservation import assert_conserved, assert_nonnegative, total_carbon

# A short, dry-by-then ferment (14 d at 20 C takes a 24-Brix must to dryness), then a warm aging
# tail. begin_aging sits well past dryness so the flux-gated ester producers are quiescent and the
# aging effect is the ONLY thing moving esters/fusels/Byp (Stance A). Warm (25 C) aging + a
# multi-month tail makes the hydrolysis measurable within a fast test.
_FERMENT_DAYS = 30.0
_AGING_DAYS = 150.0


def _wine(interventions: list[Intervention], *, aging_celsius: float = 25.0) -> Scenario:
    return Scenario(
        name="aging-test",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25},
        # A gentle 20->aging_celsius ramp across the ferment window (the two knots compile to one
        # linear segment), held flat at the warmer aging temperature for the tail (schedule holds
        # the last knot). The ferment reaches dryness well before the begin_aging breakpoint, so
        # the aging segment runs at the warmer hold — enough to make the hydrolysis measurable.
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=_FERMENT_DAYS, celsius=aging_celsius),
        ],
        interventions=interventions,
        duration_days=_FERMENT_DAYS + _AGING_DAYS,
    )


def _beer(interventions: list[Intervention]) -> Scenario:
    # A minimal ~1.048-OG wort (glucose+maltose+maltotriose), fermented ~14 d then aged. Aging is
    # medium-agnostic (esters/fusels/Byp exist in the beer schema too), so begin_aging must drive
    # the beer scenario path just as it does wine — this is the beer-side smoke coverage.
    return Scenario(
        name="beer-aging-test",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 150.0,
            "pitch_gpl": 0.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=interventions,
        duration_days=14.0 + 120.0,
    )


def _begin_aging(day: float) -> Intervention:
    return Intervention(day=day, action="begin_aging")


# -- the compile-seam enable/disable gate -------------------------------------


def test_ester_hydrolysis_is_disabled_without_begin_aging():
    # Wired into the medium (present in the set) but DISABLED at compile — aging is inherently
    # post-ferment, so an un-aged scenario never activates it (byte-for-byte the pre-aging core).
    cs = compile_scenario(_wine([]))
    assert EsterHydrolysis.name in cs.process_set
    assert not cs.process_set.is_enabled(EsterHydrolysis.name)


def test_begin_aging_enables_ester_hydrolysis_at_the_breakpoint():
    # The reconfigure enables it. Compiling leaves it disabled (the pre-breakpoint segment); the
    # event's reconfigure callback turns it on, exactly as pitch_mlf enables the MLF Processes.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    assert not cs.process_set.is_enabled(EsterHydrolysis.name)  # still off at compile
    aging_events = [e for e in cs.events if e.label.startswith("begin_aging")]
    assert len(aging_events) == 1
    event = aging_events[0]
    # A pure phase switch: it reconfigures but mutates no state (aging inoculates nothing).
    assert event.mutate is None
    assert event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(EsterHydrolysis.name)


def test_aging_params_ride_in_every_compiled_scenario():
    # aging.yaml is in shared_files, so the params are present even on an un-aged scenario —
    # inert (read by nothing until begin_aging enables the Process), but loaded.
    cs = compile_scenario(_wine([]))
    for name in ("k_ester_hydrolysis", "E_a_ester_hydrolysis", "esters_eq"):
        assert name in cs.parameters


# -- the emergent aging headline (esters fade, fusels + Byp rise) -------------


def test_aging_fades_esters_and_raises_fusels_end_to_end():
    # The end-to-end payoff: an aged wine finishes with LOWER esters and HIGHER fusels/Byp than
    # the otherwise-identical un-aged wine. Both runs share the identical ferment (aging is off
    # until the breakpoint), so any difference at the end is exactly the aging Process's doing —
    # a clean A/B that also proves isolability (the un-aged run never activates the Process).
    aged = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    plain = compile_scenario(_wine([])).run()
    assert aged.success and plain.success

    esters_aged = float(aged.series("esters")[-1])
    esters_plain = float(plain.series("esters")[-1])
    fusels_aged = float(aged.series("fusels")[-1])
    fusels_plain = float(plain.series("fusels")[-1])
    byp_aged = float(aged.series("Byp")[-1])
    byp_plain = float(plain.series("Byp")[-1])

    # The wine actually made ester during the ferment, so there is something to hydrolyse.
    assert esters_plain > 0.0
    # Aging hydrolyses the fruity acetate esters: fewer esters, more fusels + Byp at the end.
    assert esters_aged < esters_plain
    assert fusels_aged > fusels_plain
    assert byp_aged > byp_plain


def test_aging_does_not_strip_esters_below_the_equilibrium_floor():
    # Net decay toward a LOWER floor, not decay-to-zero (D-68): even over a long warm aging tail
    # the esters pool relaxes toward, not past, esters_eq.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    traj = cs.run()
    assert traj.success
    assert float(traj.series("esters")[-1]) >= cs.param_values["esters_eq"]


# -- conservation + non-negativity end-to-end --------------------------------


def test_aged_run_closes_carbon_end_to_end():
    # begin_aging mutates no state (a pure reconfigure), so there is NO external flow: the
    # run-wide invariant is the plain final == initial carbon. The ferment routes carbon into the
    # aroma pools and the aging segment transfers esters -> fusels + Byp — both close to machine
    # precision, so total_carbon is flat across the whole ferment+aging trajectory.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    traj = cs.run()
    assert traj.success
    assert traj.external_flows == ()  # a pure phase switch injects/removes no mass

    f_c = cs.parameters.value("biomass_C_fraction")
    schema = cs.schema
    assert_conserved(
        traj.as_trajectory(), total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_nonnegative(traj.as_trajectory(), ("esters", "fusels", "Byp"), atol=1e-9)


def test_slow_phase_integration_succeeds_over_the_long_span():
    # The §7 multi-scale concern (do not integrate years at ferment resolution) is answered by
    # the segment restart: the ferment segment integrates at fine resolution, and after the
    # breakpoint the BDF solver re-inits and takes large steps across the quiescent aging segment
    # (flux gone at dryness). The full wine set is active throughout — verify the shipped config,
    # not the isolated D-69 case — over ~half a year.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], aging_celsius=25.0))
    traj = cs.run()
    assert traj.success, traj.message
    # The aging breakpoint appears in the segment bounds (ferment | aging).
    assert cs.events  # a ramp knot + the begin_aging event
    assert any(np.isclose(b, _FERMENT_DAYS * 24.0) for b in traj.segment_bounds)


# -- tier propagation (min-combined across the whole run) ---------------------


def test_aging_floors_touched_pools_at_speculative_for_the_whole_run():
    # EsterHydrolysis is enabled only for the aging back half, but simulate_scheduled min-combines
    # the per-segment tier maps (D-35): a run is only as trustworthy as its least-trustworthy
    # segment, so the speculative aging Process drags esters/fusels/Byp to speculative for the
    # WHOLE run — not just the aging segment.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    traj = cs.run()
    for pool in ("esters", "fusels", "Byp"):
        assert traj.tier_map[pool] is Tier.SPECULATIVE


# -- the vocabulary boundary (loud errors) ------------------------------------


def test_begin_aging_rejects_unknown_params():
    # A pure phase switch takes no params; a stray param is a typo, rejected loudly (the verb
    # registry discipline), not silently ignored.
    scenario = _wine([Intervention(day=_FERMENT_DAYS, action="begin_aging", params={"months": 6})])
    with pytest.raises(ValueError, match="unknown param"):
        compile_scenario(scenario)


def test_begin_aging_without_aging_params_fails_loudly_at_compile():
    # A caller-supplied parameter_paths WITHOUT aging.yaml would otherwise surface as a bare
    # KeyError deep in EsterHydrolysis.derivatives mid-integration; the verb guards the params at
    # compile (the add_dap/additions.yaml pattern) so the error names the scenario and the fix.
    data = default_data_dir()
    with pytest.raises(ValueError, match="aging.yaml"):
        compile_scenario(
            _wine([_begin_aging(_FERMENT_DAYS)]),
            parameter_paths=[data / "wine_generic.yaml", data / "acidbase.yaml"],
        )


# -- medium-agnostic: the beer scenario path too -------------------------------


def test_begin_aging_drives_the_beer_scenario_path():
    # Aging is medium-agnostic (D-70): esters/fusels/Byp exist in the beer schema, so begin_aging
    # must compile, enable EsterHydrolysis, and fade the beer esters exactly as it does for wine.
    # Beer-side smoke coverage — the D-69 Process math is beer-tested, but the beer scenario
    # compile -> begin_aging -> run path was otherwise unexercised.
    cs = compile_scenario(_beer([_begin_aging(14.0)]))
    assert not cs.process_set.is_enabled(EsterHydrolysis.name)  # disabled at compile
    aged = cs.run()
    plain = compile_scenario(_beer([])).run()
    assert aged.success and plain.success
    esters_plain = float(plain.series("esters")[-1])
    assert esters_plain > 0.0  # the beer ferment made ester to hydrolyse
    assert float(aged.series("esters")[-1]) < esters_plain  # aging fades it
    assert float(aged.series("fusels")[-1]) > float(plain.series("fusels")[-1])
