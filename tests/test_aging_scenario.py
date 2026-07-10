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

from fermentation.core.kinetics.aging import (
    EsterHydrolysis,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    StreckerDegradation,
    SulfiteOxidation,
)
from fermentation.core.media import get_medium
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.sensory import load_thresholds, oav_series
from fermentation.validation.conservation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)

# A short, dry-by-then ferment (14 d at 20 C takes a 24-Brix must to dryness), then a warm aging
# tail. begin_aging sits well past dryness so the flux-gated ester producers are quiescent and the
# aging effect is the ONLY thing moving esters/fusels/Byp (Stance A). Warm (25 C) aging + a
# multi-month tail makes the hydrolysis measurable within a fast test.
_FERMENT_DAYS = 30.0
_AGING_DAYS = 150.0


def _wine(
    interventions: list[Intervention],
    *,
    aging_celsius: float = 25.0,
    amino_acids_gpl: float = 0.0,
) -> Scenario:
    # amino_acids_gpl (default 0, byte-for-byte the pre-D-75 helper) doses the assimilable amino
    # must input the StreckerDegradation Process (D-75) draws its aldehyde carbon from; a residual
    # survives to the aging segment (AminoAcidAssimilation only draws it during active ferment).
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25}
    if amino_acids_gpl > 0.0:
        initial["amino_acids_gpl"] = amino_acids_gpl
    return Scenario(
        name="aging-test",
        medium="wine",
        initial=initial,
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


def _add_oxygen(day: float, o2_mgl: float) -> Intervention:
    return Intervention(day=day, action="add_oxygen", params={"o2_mgl": o2_mgl})


def _add_so2(day: float, so2_mgl: float) -> Intervention:
    return Intervention(day=day, action="add_so2", params={"so2_mgl": so2_mgl})


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


# -- oxidative aging: the O₂ substrate + add_oxygen (decision D-71) ------------


def test_oxidative_acetaldehyde_disabled_and_gated_with_ester_hydrolysis():
    # OxidativeAcetaldehyde rides the same aging tuple: wired into the medium, DISABLED at compile,
    # enabled by the SAME begin_aging reconfigure as EsterHydrolysis (one gate for the aging axis).
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    assert OxidativeAcetaldehyde.name in cs.process_set
    assert not cs.process_set.is_enabled(OxidativeAcetaldehyde.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(OxidativeAcetaldehyde.name)  # begin_aging turns it on too


def test_add_oxygen_doses_the_o2_slot():
    # add_oxygen is a pure carbon-free dose onto the o2 slot (the add_so2 pattern): it mutates o2,
    # reconfigures nothing, and — the o2 pool being off every ledger — books no external flow.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 40.0)]))
    dose_events = [e for e in cs.events if e.label.startswith("add_oxygen")]
    assert len(dose_events) == 1
    event = dose_events[0]
    assert event.reconfigure is None and event.mutate is not None
    before = cs.y0.copy()
    after = event.mutate(cs.schema, before)
    # 40 mg/L O₂ = 0.04 g/L lands on the o2 slot; nothing else moves.
    assert cs.schema.get(after, "o2") - cs.schema.get(before, "o2") == pytest.approx(0.04)


def test_oxidative_aging_raises_acetaldehyde_and_depletes_oxygen():
    # The end-to-end payoff: an oxygen-dosed aged wine finishes with HIGHER acetaldehyde (the
    # 'sherry'/oxidised note) than the otherwise-identical reductive (no-O₂) aged wine, and the
    # dosed O₂ is consumed over the aging tail — the saturating O₂-limited oxidation, end to end.
    o2_dose = 60.0  # mg/L — a generous cumulative aerobic exposure
    oxidative = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, o2_dose)])
    ).run()
    reductive = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()  # no O₂ dosed
    assert oxidative.success and reductive.success

    # Oxidation raises acetaldehyde above the reductive-aging baseline.
    assert float(oxidative.series("acetaldehyde")[-1]) > float(reductive.series("acetaldehyde")[-1])
    # The dosed O₂ is largely consumed by the end of the aging tail.
    assert float(oxidative.series("o2")[-1]) < 0.5 * (o2_dose / 1000.0)
    # The reductive run never accrued any O₂ (byte-for-byte no oxidation substrate).
    assert float(reductive.series("o2")[-1]) == 0.0


def test_reductive_aging_leaves_acetaldehyde_byte_for_byte():
    # Isolability (D-71): a begin_aging run WITHOUT add_oxygen is purely reductive aging — the
    # oxidation Process is inert at o2=0 — so acetaldehyde ends exactly where the un-aged run leaves
    # it (aging draws no acetaldehyde via ester hydrolysis, and viable-X-gated production/reduction
    # are quiescent post-dryness). The oxidation Process cannot move acetaldehyde without O₂.
    reductive = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    plain = compile_scenario(_wine([])).run()
    assert reductive.success and plain.success
    assert float(reductive.series("acetaldehyde")[-1]) == pytest.approx(
        float(plain.series("acetaldehyde")[-1]), rel=1e-9
    )


def test_oxygen_dosed_run_closes_carbon_end_to_end():
    # add_oxygen mutates the o2 slot, so the runtime books an external flow for it — but o2 is off
    # every ledger, so that flow is CARBON-FREE (the add_so2/add_dap idiom), and the oxidation
    # transfers E → acetaldehyde carbon-exactly. So total_carbon is flat (final == initial, no
    # ledger correction) across the whole ferment+aging+oxidation trajectory.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)]))
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    # The O₂ dose flow exists but injects zero carbon (o2 is off the carbon ledger).
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert_conserved(traj.as_trajectory(), c_of, label="carbon")
    assert_nonnegative(traj.as_trajectory(), ("o2", "acetaldehyde"), atol=1e-9)


def test_oxidative_aging_raises_the_acetaldehyde_oav():
    # Close the design loop through the STATED acceptance lens (milestone-3-plan: aging Processes
    # are "validated by the D-67 OAV lens"). The whole justification for this Process (D-71) was
    # that it moves an OAV the lens already reads — so assert the acetaldehyde OAV itself climbs on
    # the oxygen-dosed run vs the reductive baseline, not merely the raw pool. OAV = conc/threshold,
    # so it tracks the pool, but this ties the built behaviour to the reason it was built.
    thresholds = load_thresholds()
    oxidative = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)])
    ).run()
    reductive = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    assert oxidative.success and reductive.success

    oav_ox = float(oav_series(oxidative.as_trajectory(), thresholds, "acetaldehyde")[-1])
    oav_red = float(oav_series(reductive.as_trajectory(), thresholds, "acetaldehyde")[-1])
    # Oxidative aging lifts the acetaldehyde OAV (the 'oxidised' aroma the lens reads) to a positive
    # value above the reductive baseline — the sensory signature of the O₂-driven ethanol oxidation.
    # (The reductive baseline sits at ~0, a hair negative only from solver undershoot near an empty
    # pool, so compare against it directly rather than pinning its sign.)
    assert oav_ox > 0.0
    assert oav_ox > oav_red


def test_add_oxygen_rejects_unknown_params():
    # The verb takes exactly {o2_mgl}; a stray param is a typo, rejected loudly (verb registry
    # discipline).
    scenario = _wine(
        [Intervention(day=_FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": 40.0, "ppm": 5})]
    )
    with pytest.raises(ValueError, match="unknown param"):
        compile_scenario(scenario)


# -- oxidative aging: SO₂ scavenging on the shared O₂ budget (decision D-72) ---


def test_sulfite_oxidation_gated_by_begin_aging_wine_only():
    # SulfiteOxidation is WINE-ONLY (reads wine-only so2_total/pH slots) — present in the wine set,
    # absent from beer — and rides the SAME aging gate: disabled at compile, enabled by begin_aging.
    assert SulfiteOxidation.name in get_medium("wine").build_process_set()
    assert SulfiteOxidation.name not in get_medium("beer").build_process_set()

    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    assert SulfiteOxidation.name in cs.process_set
    assert not cs.process_set.is_enabled(SulfiteOxidation.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(SulfiteOxidation.name)  # begin_aging turns it on too


def test_so2_suppresses_oxidative_acetaldehyde_end_to_end():
    # THE HEADLINE end-to-end (D-72): dose the SAME O₂ charge on two aged wines, one also dosed with
    # SO₂ at the aging breakpoint. The SO₂ scavenges O₂ (bisulfite out-competes ethanol for it), so
    # the sulfited wine finishes with LOWER oxidative acetaldehyde and CONSUMES its SO₂ — the "SO₂
    # protects until exhausted" threshold, driven through the whole compile→schedule→run pipeline.
    o2_dose = 60.0
    day = _FERMENT_DAYS
    unprotected = compile_scenario(_wine([_begin_aging(day), _add_oxygen(day, o2_dose)])).run()
    protected = compile_scenario(
        _wine([_begin_aging(day), _add_oxygen(day, o2_dose), _add_so2(day, 200.0)])
    ).run()
    assert unprotected.success and protected.success

    # SO₂ present ⇒ less oxidative acetaldehyde than the identical wine without it.
    assert float(protected.series("acetaldehyde")[-1]) < float(
        unprotected.series("acetaldehyde")[-1]
    )
    # The protective SO₂ is genuinely consumed defending the wine (some of the 200 mg/L is burned).
    so2_end_gpl = float(protected.series("so2_total")[-1])
    assert so2_end_gpl < 0.200  # below the dosed 200 mg/L — SO₂ was spent oxidising
    # And it did its job: the dosed O₂ is still largely consumed (diverted to SO₂, not left over).
    assert float(protected.series("o2")[-1]) < 0.5 * (o2_dose / 1000.0)


def test_so2_dosed_oxidative_run_closes_carbon_end_to_end():
    # add_so2 + add_oxygen both mutate off-ledger slots (so2_total, o2), so the runtime books two
    # external flows — both CARBON-FREE. The SO₂-oxidation transfer moves nothing on the carbon
    # ledger (o2→? and so2→sulfate are both off it); only OxidativeAcetaldehyde's E→acetaldehyde is
    # on it, and it closes exactly. So total_carbon stays flat across the whole trajectory.
    day = _FERMENT_DAYS
    cs = compile_scenario(_wine([_begin_aging(day), _add_oxygen(day, 60.0), _add_so2(day, 150.0)]))
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    # Every dose flow (O₂ and SO₂) injects zero carbon.
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert_conserved(traj.as_trajectory(), c_of, label="carbon")
    assert_nonnegative(traj.as_trajectory(), ("o2", "so2_total", "acetaldehyde"), atol=1e-9)


# -- PhenolicBrowning (decision D-74) — the first ALWAYS-ON O₂ sink, end to end ------------


def test_phenolic_browning_disabled_and_gated_with_begin_aging():
    # PhenolicBrowning rides the same aging tuple: wired into the medium, DISABLED at compile,
    # enabled by the SAME begin_aging reconfigure as the other aging Processes (one gate).
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    assert PhenolicBrowning.name in cs.process_set
    assert not cs.process_set.is_enabled(PhenolicBrowning.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(PhenolicBrowning.name)  # begin_aging turns it on too


def test_oxidative_aging_browns_the_wine_end_to_end():
    # The D-74 payoff: an oxygen-dosed aged wine finishes BROWNED (A420 climbs from 0) — the
    # gold→amber→brown of oxidative aging — while the otherwise-identical reductive (no-O₂) aged
    # wine
    # never browns (A420 stays exactly 0, the byte-for-byte reductive isolability).
    o2_dose = 60.0  # mg/L cumulative aerobic exposure
    oxidative = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, o2_dose)])
    ).run()
    reductive = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()  # no O₂ dosed
    assert oxidative.success and reductive.success

    # Oxidative aging builds visible brown; reductive aging browns none.
    assert float(oxidative.series("A420")[-1]) > 0.0
    assert float(reductive.series("A420")[-1]) == 0.0


def test_reductive_aging_leaves_a420_at_zero_byte_for_byte():
    # Isolability (D-74): begin_aging WITHOUT add_oxygen is purely reductive — browning is inert at
    # o2=0, so A420 stays exactly 0, the same as an un-aged run (which never even enables browning).
    reductive = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    plain = compile_scenario(_wine([])).run()
    assert reductive.success and plain.success
    assert float(reductive.series("A420")[-1]) == 0.0
    assert float(plain.series("A420")[-1]) == 0.0


def test_browning_suppresses_oxidative_acetaldehyde_end_to_end():
    # The headline competition, end to end: browning is the DOMINANT always-on O₂ sink, so an
    # oxygen-dosed aged wine finishes with the O₂ SPLIT between brown pigment (A420 > 0) and
    # acetaldehyde — the acetaldehyde suppressed to its ~40% ethanol share (the balance browned).
    # Isolating browning off would raise acetaldehyde; here we assert the co-resident sinks coexist:
    # both observables are positive, and the O₂ is (largely) consumed by their sum.
    o2_dose = 60.0
    aged = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, o2_dose)])
    ).run()
    assert aged.success
    # Both oxidative products form (O₂ split between them) and the dosed O₂ is largely spent by the
    # two summed sinks (browning taking the majority).
    assert float(aged.series("A420")[-1]) > 0.0
    assert float(aged.series("acetaldehyde")[-1]) > 0.0
    assert float(aged.series("o2")[-1]) < 0.5 * (o2_dose / 1000.0)


def test_begin_aging_browns_the_beer_scenario():
    # MEDIUM-AGNOSTIC end to end (D-74, superseding D-73's provisional wine-only): beer carries
    # autoxidising polyphenols and browns too, and A420 exists in the beer schema — so an
    # oxygen-dosed aged beer browns (A420 climbs), the beer-side smoke coverage for browning.
    day = 14.0
    oxidative = compile_scenario(_beer([_begin_aging(day), _add_oxygen(day, 40.0)])).run()
    reductive = compile_scenario(_beer([_begin_aging(day)])).run()
    assert oxidative.success and reductive.success
    assert float(oxidative.series("A420")[-1]) > 0.0
    assert float(reductive.series("A420")[-1]) == 0.0


def test_browned_run_closes_carbon_end_to_end():
    # PhenolicBrowning touches only o2 + A420, BOTH off every ledger — so it moves nothing
    # conserved.
    # An oxygen-dosed aged run (browning + oxidative acetaldehyde + ester hydrolysis all active)
    # still
    # closes total_carbon exactly: the O₂ dose flow is carbon-free, browning adds no carbon term,
    # and
    # the only on-ledger moves (E→acetaldehyde, the ester 5:2 transfer) close to machine precision.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)]))
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert_conserved(traj.as_trajectory(), c_of, label="carbon")
    assert_nonnegative(traj.as_trajectory(), ("o2", "A420", "acetaldehyde"), atol=1e-9)


# -- StreckerDegradation (decision D-75) — the O₂/amino-acid Strecker aldehydes, end to end ----


def test_strecker_gated_by_begin_aging_wine_only():
    # StreckerDegradation is WINE-ONLY (reads wine-only amino_acids + deaminates to N) — present in
    # the wine set, absent from beer — and rides the aging gate: disabled at compile, then on by
    # begin_aging (the SulfiteOxidation pattern).
    assert StreckerDegradation.name in get_medium("wine").build_process_set()
    assert StreckerDegradation.name not in get_medium("beer").build_process_set()
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], amino_acids_gpl=0.5))
    assert StreckerDegradation.name in cs.process_set
    assert not cs.process_set.is_enabled(StreckerDegradation.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(StreckerDegradation.name)  # begin_aging turns it on


def test_strecker_produces_aldehydes_with_oxygen_and_amino_acids():
    # The end-to-end payoff: an oxygen-dosed, amino-acid-dosed aged wine finishes with BOTH Strecker
    # aldehydes accumulated (methional the cooked-potato off-note, phenylacetaldehyde the honey),
    # methional-dominant (f_methional = 0.6), and the dosed O₂ largely consumed over the aging tail.
    o2_dose = 60.0
    aged = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, o2_dose)],
            amino_acids_gpl=0.5,
        )
    ).run()
    assert aged.success
    methional = float(aged.series("methional")[-1])
    phenyl = float(aged.series("phenylacetaldehyde")[-1])
    # Both aldehydes form, at aroma-relevant (µg/L-scale) levels; methional dominates the split.
    assert methional > phenyl > 0.0
    # The dosed O₂ is largely consumed by the end of the aging tail (Strecker rides the shared O₂
    # alongside the dominant browning + ethanol-oxidation sinks).
    assert float(aged.series("o2")[-1]) < 0.5 * (o2_dose / 1000.0)


def test_strecker_silent_without_amino_acids():
    # Isolability (the D-75 substrate gate): an oxygen-dosed aged run with NO amino acids makes
    # NO Strecker aldehydes — the amino_acids ≤ 0 guard is exact — so a nutrient-free aging is
    # byte-for-byte the case without this Process (methional/phenylacetaldehyde stay 0).
    aged = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)])  # no amino_acids_gpl
    ).run()
    assert aged.success
    assert float(aged.series("methional")[-1]) == 0.0
    assert float(aged.series("phenylacetaldehyde")[-1]) == 0.0


def test_strecker_silent_reductive():
    # The other substrate gate: an amino-acid-dosed but REDUCTIVE (no add_oxygen) aging makes no
    # Strecker aldehydes — the o2 ≤ 0 guard is exact — so reductive aging is unchanged.
    aged = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], amino_acids_gpl=0.5)).run()
    assert aged.success
    assert float(aged.series("methional")[-1]) == 0.0
    assert float(aged.series("phenylacetaldehyde")[-1]) == 0.0


def test_strecker_closes_carbon_and_nitrogen_end_to_end():
    # StreckerDegradation draws carbon from amino_acids into methional + phenylacetaldehyde + CO₂,
    # deaminates the nitrogen to N — so BOTH ledgers must close end to end. The O₂ dose flow is
    # carbon- AND nitrogen-free (o2 off every ledger), the amino-acid dose is a t0 initial (not a
    # flow), so total_carbon and total_nitrogen are both flat (final == initial) across the whole
    # ferment + aging + Strecker trajectory.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)],
            amino_acids_gpl=0.5,
        )
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    n_of = total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n)
    # Every external flow (the O₂ dose) injects zero carbon AND zero nitrogen (o2 off every ledger).
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert all(n_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert_conserved(traj.as_trajectory(), c_of, label="carbon")
    assert_conserved(traj.as_trajectory(), n_of, label="nitrogen")
    assert_nonnegative(
        traj.as_trajectory(), ("o2", "amino_acids", "methional", "phenylacetaldehyde"), atol=1e-9
    )


def test_strecker_raises_the_strecker_oavs():
    # Close the design loop through the STATED acceptance lens (milestone-3-plan: aging Processes
    # "validated by the D-67 OAV lens"). The whole point of D-75 was to add the two Strecker aromas
    # the lens now reads — so assert BOTH OAVs climb positive on the oxygen+amino-acid-dosed run vs
    # the reductive baseline (which is byte-for-byte 0). methional's low (~0.5 µg/L) threshold makes
    # its OAV the larger of the two.
    thresholds = load_thresholds()
    aged = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)],
            amino_acids_gpl=0.5,
        )
    ).run()
    reductive = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], amino_acids_gpl=0.5)).run()
    assert aged.success and reductive.success

    for pool in ("methional", "phenylacetaldehyde"):
        oav_aged = float(oav_series(aged.as_trajectory(), thresholds, pool)[-1])
        oav_red = float(oav_series(reductive.as_trajectory(), thresholds, pool)[-1])
        assert oav_aged > 0.0
        assert oav_red == 0.0  # reductive aging raises neither aroma
        assert oav_aged > oav_red
