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

from fermentation.analysis import (
    astringency_series,
    color_series,
    observed_color_series,
    polymeric_pigment_series,
)
from fermentation.core.kinetics.aging import (
    AcetaldehydeBridgedCondensation,
    EllagitanninOxidation,
    EsterHydrolysis,
    OakExtraction,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    StreckerDegradation,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
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
    autolysis_rate_per_h: float = 0.0,
    anthocyanin_gpl: float = 0.0,
    tannin_gpl: float = 0.0,
) -> Scenario:
    # amino_acids_gpl (default 0, byte-for-byte the pre-D-75 helper) doses the assimilable amino
    # must input the StreckerDegradation Process (D-75) draws its aldehyde carbon from; a residual
    # survives to the aging segment (AminoAcidAssimilation only draws it during active ferment).
    # autolysis_rate_per_h (default 0, byte-for-byte the pre-D-76 helper) instead opts into lees
    # autolysis (D-34): dead biomass self-digests post-dryness, REFILLING amino_acids from the
    # physically-real nitrogen source — the emergent sur-lie → Strecker pathway (D-76), no dose.
    # anthocyanin_gpl/tannin_gpl (default 0, byte-for-byte the pre-D-79 helper) dose the grape must
    # inputs the TanninAnthocyaninCondensation Process (D-79) condenses into stable polymeric
    # pigment
    # during aging — a red wine (both > 0) softens + stabilizes colour; a white (both 0) is inert.
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25}
    if amino_acids_gpl > 0.0:
        initial["amino_acids_gpl"] = amino_acids_gpl
    if autolysis_rate_per_h > 0.0:
        initial["autolysis_rate_per_h"] = autolysis_rate_per_h
    if anthocyanin_gpl > 0.0:
        initial["anthocyanin_gpl"] = anthocyanin_gpl
    if tannin_gpl > 0.0:
        initial["tannin_gpl"] = tannin_gpl
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


def _add_oak(day: float, oak_gpl: float, toast: str) -> Intervention:
    return Intervention(day=day, action="add_oak", params={"oak_gpl": oak_gpl, "toast": toast})


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
    # phenylacetaldehyde-dominant (f_methional = 0.15 — phenylalanine is the more abundant must
    # precursor), and the dosed O₂ largely consumed over the aging tail.
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
    # Both aldehydes form, at aroma-relevant (µg/L-scale) levels; phenylacetaldehyde dominates the
    # split, methional the potent low-µg/L minority. Both land in the oxidised-white-wine range.
    assert phenyl > methional > 0.0
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
    # the reductive baseline (which is byte-for-byte 0). Both clear their thresholds; the split is
    # phenylacetaldehyde-dominant so both are perceptible without asserting an ordering here.
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


# -- D-76: the emergent sur-lie → Strecker pathway (autolysis refills amino_acids, no dose) -----
#
# D-75 exercised Strecker by DOSING amino_acids_gpl (an artificial nutrient add). D-76 closes the
# physically-real loop: opting into lees autolysis (D-34, autolysis_rate_per_h) lets dead biomass
# self-digest post-dryness and REFILL amino_acids — the very substrate the O₂/quinone Strecker
# route (D-75) draws — so Strecker is non-silent from the sur-lie nitrogen source, not a dose. No
# new physics: the refill Process (D-34) and the consumer (D-75) simply COMPOSE. The A-vs-B design
# fork (owner-decided: A) is recorded in DECISIONS D-76; the measurement behind it is that the
# active-ferment (pre-dryness) release is small (~15 mg/L) while the aging pool is dominated by
# legit post-dryness sur-lie autolysis, so autolysis-from-t0 needs no re-gating.


def test_sur_lie_autolysis_feeds_strecker_without_a_dose():
    # The D-76 headline: with NO amino_acids_gpl, opting into autolysis + O₂ + begin_aging still
    # makes BOTH Strecker aldehydes — dead lees autolyse, refilling amino_acids (D-34), which the
    # Strecker route (D-75) degrades. The nitrogen is the wine's own dead yeast, not a dose.
    aged = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)],
            autolysis_rate_per_h=1.0e-3,  # sur-lie refill; NO amino_acids_gpl dose
        )
    ).run()
    assert aged.success
    # Autolysis refilled the pool with no dose (measured ~0.83 g/L at the end of the sur-lie tail).
    assert float(aged.series("amino_acids")[-1]) > 0.1
    # Both Strecker aldehydes emerge from that autolytic nitrogen, phenylacetaldehyde-dominant. This
    # is a DIRECTIONAL assertion only (Tier-3 discipline): the absolute level (~154/1006 µg/L) runs
    # ~8× the D-75 dosed literature anchor because the autolytic pool (~0.4–0.8 g/L) floods the sink
    # vs the ~11.6 mg/L dosed case k_strecker was calibrated against — an order-of-magnitude figure,
    # not a prediction (the arginine-lump over-feed is a recorded open item, DECISIONS D-76).
    methional = float(aged.series("methional")[-1])
    phenyl = float(aged.series("phenylacetaldehyde")[-1])
    assert phenyl > methional > 0.0


def test_sur_lie_pathway_is_autolysis_driven_not_an_artifact():
    # The isolability contrast proving the pathway IS the autolysis refill: the same O₂-dosed aged
    # scenario with NEITHER autolysis NOR a dose keeps amino_acids exactly 0 (no producer) and
    # makes NO Strecker aldehydes — turning autolysis on is exactly what lights the pathway.
    plain = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)])  # no autolysis/dose
    ).run()
    assert plain.success
    assert float(plain.series("amino_acids")[-1]) == 0.0
    assert float(plain.series("methional")[-1]) == 0.0
    assert float(plain.series("phenylacetaldehyde")[-1]) == 0.0


def test_sur_lie_strecker_closes_carbon_and_nitrogen_end_to_end():
    # The NEW conservation combination D-76 introduces: autolysis REFILLS amino_acids (releasing
    # dead-cell N as arginine, routing the C-rich remainder to debris) WHILE StreckerDegradation
    # AND AutolyticMercaptan both DRAW it (arginine C → aldehydes/thiol + CO₂, N deaminated back to
    # N). D-34 and D-75 each pinned their half; this pins them COMPOSED. Both ledgers stay flat
    # (final == initial): the O₂ dose is the only external flow and carries neither element.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)],
            autolysis_rate_per_h=1.0e-3,
        )
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    n_of = total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n)
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert all(n_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert_conserved(traj.as_trajectory(), c_of, label="carbon")
    assert_conserved(traj.as_trajectory(), n_of, label="nitrogen")
    assert_nonnegative(
        traj.as_trajectory(),
        ("amino_acids", "debris", "methional", "phenylacetaldehyde", "mercaptans", "h2s"),
        atol=1e-9,
    )


# -- OakExtraction (decision D-77) — the NON-oxidative barrel/chip aroma axis, end to end -------

_OAK_EXTRACTIVES = ("whiskey_lactone", "vanillin", "guaiacol", "eugenol")
_OAK_CEILINGS = tuple(f"{c}_ceiling" for c in _OAK_EXTRACTIVES)


def test_oak_extraction_gated_by_begin_aging_both_media():
    # OakExtraction rides the same aging tuple: DISABLED at compile, enabled by the SAME begin_aging
    # reconfigure as the rest of the aging axis. Barrel-beer oak (D-86): the oak axis is a wood
    # property, so — unlike the wine-only SulfiteOxidation/StreckerDegradation — it is wired into
    # BOTH media and present (disabled-then-enabled) in each.
    for scenario in (_wine([_begin_aging(_FERMENT_DAYS)]), _beer([_begin_aging(14.0)])):
        cs = compile_scenario(scenario)
        assert OakExtraction.name in cs.process_set
        assert not cs.process_set.is_enabled(OakExtraction.name)  # off at compile (post-ferment)
        event = next(e for e in cs.events if e.label.startswith("begin_aging"))
        assert event.reconfigure is not None
        event.reconfigure(cs.process_set)
        assert cs.process_set.is_enabled(OakExtraction.name)  # begin_aging turns it on


def test_oak_params_ride_in_every_compiled_scenario():
    # oak.yaml is in shared_files, so the extraction constants are present even on an un-oaked
    # scenario — inert (read by nothing until add_oak + begin_aging), but loaded.
    cs = compile_scenario(_wine([]))
    for name in ("k_oak_extraction", "E_a_oak_extraction", "oak_yield_vanillin_medium"):
        assert name in cs.parameters


def test_add_oak_sets_the_ceilings_from_toast_yields():
    # add_oak is a pure dose onto the SET-AND-HOLD ceiling slots (the add_oxygen pattern): it
    # mutates the four ceilings to oak_gpl × oak_yield_<compound>_<toast>, reconfigures nothing, and
    # — the oak slots being off every ledger — perturbs no elemental balance. Medium toast at 4 g/L.
    cs = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium")])
    )
    event = next(e for e in cs.events if e.label.startswith("add_oak"))
    assert event.reconfigure is None and event.mutate is not None  # dose only (begin_aging enables)
    after = event.mutate(cs.schema, cs.y0.copy())
    for compound in _OAK_EXTRACTIVES:
        expected = 4.0 * cs.param_values[f"oak_yield_{compound}_medium"]
        assert cs.schema.get(after, f"{compound}_ceiling") == pytest.approx(expected)
        assert cs.schema.get(after, compound) == 0.0  # the extractive itself starts empty


def test_toast_selects_the_aroma_profile():
    # The load-bearing toast ORDERING (D-77): light toast is whiskey-lactone (coconut) dominant;
    # heavy toast is guaiacol (smoky) + eugenol (clove) dominant; vanillin peaks at medium. Compare
    # the finished ceilings across toast levels at the same oak dose — a discriminating check that
    # the categorical toast genuinely reshapes the profile, not just its scale.
    def ceilings(toast: str) -> dict[str, float]:
        cs = compile_scenario(
            _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, toast)])
        )
        event = next(e for e in cs.events if e.label.startswith("add_oak"))
        assert event.mutate is not None
        after = event.mutate(cs.schema, cs.y0.copy())
        return {c: float(cs.schema.get(after, f"{c}_ceiling")) for c in _OAK_EXTRACTIVES}

    light, medium, heavy = ceilings("light"), ceilings("medium"), ceilings("heavy")
    # Whiskey lactone (coconut) falls with toast; guaiacol + eugenol (smoky/clove) rise with toast.
    assert light["whiskey_lactone"] > medium["whiskey_lactone"] > heavy["whiskey_lactone"]
    assert heavy["guaiacol"] > medium["guaiacol"] > light["guaiacol"]
    assert heavy["eugenol"] > medium["eugenol"] > light["eugenol"]
    # Vanillin (vanilla) peaks at medium toast (lignin thermal release), not at the extremes.
    assert medium["vanillin"] > light["vanillin"] and medium["vanillin"] > heavy["vanillin"]


def test_oak_extraction_builds_the_extractives_end_to_end():
    # The end-to-end payoff: an oaked aged wine finishes with the four extractives risen from 0
    # toward their ceilings, while the otherwise-identical un-oaked aged wine has none — a clean A/B
    # (both share the identical ferment + aging; the only difference is the add_oak dose).
    oaked = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium")])
    ).run()
    plain = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    assert oaked.success and plain.success
    for compound in _OAK_EXTRACTIVES:
        oaked_end = float(oaked.series(compound)[-1])
        ceiling = float(oaked.series(f"{compound}_ceiling")[-1])
        assert 0.0 < oaked_end <= ceiling + 1e-15  # rose from 0, toward (not past) the ceiling
        assert float(plain.series(compound)[-1]) == 0.0  # un-oaked wine extracts nothing


def test_un_oaked_aging_leaves_the_oak_pools_zero_but_speculative():
    # The three-case isolability (D-77): begin_aging WITHOUT add_oak leaves every oak extractive and
    # ceiling identically 0 (byte-for-byte the case without oak) — BUT the pools still report
    # SPECULATIVE like the rest of the enabled aging axis (the Process is enabled, just un-dosed;
    # tier_of counts enabled, not nonzero, Processes). Zero value, speculative tier — both correct.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    traj = cs.run()
    assert traj.success
    for name in _OAK_EXTRACTIVES + _OAK_CEILINGS:
        assert np.max(np.abs(traj.series(name))) == 0.0  # identically zero — no oak dosed
    for compound in _OAK_EXTRACTIVES:
        assert traj.tier_map[compound] is Tier.SPECULATIVE  # enabled aging axis ⇒ speculative floor


def test_oaked_run_closes_every_ledger_end_to_end():
    # Oak extractives + ceilings are ALL off every ledger (exogenous wood-derived, the iso_alpha
    # precedent), so add_oak injects nothing conserved (its ExternalFlow carries no carbon/nitrogen)
    # and OakExtraction moves nothing conserved — the whole ferment+aging trajectory closes carbon
    # AND nitrogen to machine precision with the oak axis fully active.
    cs = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "heavy")])
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    n_of = total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n)
    # The add_oak ceiling jump carries neither element (off every ledger).
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert all(n_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert_conserved(traj.as_trajectory(), c_of, label="carbon")
    assert_conserved(traj.as_trajectory(), n_of, label="nitrogen")
    assert_nonnegative(traj.as_trajectory(), _OAK_EXTRACTIVES, atol=1e-15)


# -- EllagitanninOxidation (decision D-78) — oak PROTECTS the wine, end to end ------------------
# The user-facing path the unit tests in test_aging.py never drive: add_oak → compile →
# begin_aging + add_oxygen → run(). These pin the verb/compile wiring (add_oak sets the ellagitannin
# ceiling; begin_aging enables EllagitanninOxidation, wine-only) and — the SPINE — the
# oak-protection
# emergent through the compiled run(): an oaked+oxygenated wine browns less and makes less oxidative
# acetaldehyde than an un-oaked wine at the same O₂ dose.


def test_ellagitannin_oxidation_gated_by_begin_aging_both_media():
    # EllagitanninOxidation rides the same aging tuple: DISABLED at compile, enabled by the SAME
    # begin_aging reconfigure. Barrel-beer oak (D-86): the ellagitannin slot and o2 pool are both
    # medium-agnostic, so — unlike the wine-only Strecker/SulfiteOxidation — it is wired into BOTH
    # media and present (disabled-then-enabled) in each.
    for scenario in (_wine([_begin_aging(_FERMENT_DAYS)]), _beer([_begin_aging(14.0)])):
        cs = compile_scenario(scenario)
        assert EllagitanninOxidation.name in cs.process_set
        assert not cs.process_set.is_enabled(
            EllagitanninOxidation.name
        )  # off at compile (post-ferment)
        event = next(e for e in cs.events if e.label.startswith("begin_aging"))
        assert event.reconfigure is not None
        event.reconfigure(cs.process_set)
        assert cs.process_set.is_enabled(EllagitanninOxidation.name)  # begin_aging turns it on


def test_add_oak_sets_the_ellagitannin_ceiling_and_rate_params_ride_along():
    # D-78: the SAME add_oak dose that sets the four aroma ceilings now ALSO sets the ellagitannin
    # ceiling to oak_gpl × oak_yield_ellagitannin_<toast>, and the scavenging rate/E_a/yield params
    # ride in every compiled scenario (oak.yaml is in shared_files) — inert until oak is dosed.
    cs = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium")])
    )
    for name in ("k_ellagitannin_oxidation", "E_a_ellagitannin_oxidation", "y_ellag_per_o2"):
        assert name in cs.parameters
    event = next(e for e in cs.events if e.label.startswith("add_oak"))
    assert event.mutate is not None
    after = event.mutate(cs.schema, cs.y0.copy())
    expected = 4.0 * cs.param_values["oak_yield_ellagitannin_medium"]
    assert cs.schema.get(after, "ellagitannin_ceiling") == pytest.approx(expected)
    assert (
        cs.schema.get(after, "ellagitannin") == 0.0
    )  # the tannin itself starts empty (extracts in)


def test_ellagitannin_declines_with_toast():
    # The D-78 toast ORDERING: ellagitannin is THERMOLABILE — degraded by toasting — so it DECLINES
    # light > medium > heavy (heavy-toast barrels are rounder/less astringent), the same direction
    # as
    # whiskey lactone and opposite the guaiacol/eugenol pyrolysis phenols. A discriminating check
    # through the real add_oak verb at a fixed oak dose.
    def ceiling(toast: str) -> float:
        cs = compile_scenario(
            _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, toast)])
        )
        event = next(e for e in cs.events if e.label.startswith("add_oak"))
        assert event.mutate is not None
        after = event.mutate(cs.schema, cs.y0.copy())
        return float(cs.schema.get(after, "ellagitannin_ceiling"))

    assert ceiling("light") > ceiling("medium") > ceiling("heavy") > 0.0


def test_oak_protects_against_oxidation_end_to_end():
    # THE D-78 SPINE, through the compiled run() (not a hand-built ProcessSet): two identical
    # oxygenated aged wines — same 40 mg/L O₂ dose, same ferment + aging — differing ONLY in the
    # add_oak charge. The oaked wine's ellagitannin scavenges its share of the O₂, so it browns LESS
    # (lower A420) and accumulates LESS oxidative acetaldehyde. Suppression is PARTIAL (the oaked
    # wine
    # still shows SOME browning/acetaldehyde). And its astringency readout is positive (tannin
    # present).
    o2_dose = 40.0
    oaked = compile_scenario(
        _wine(
            [
                _begin_aging(_FERMENT_DAYS),
                _add_oak(
                    _FERMENT_DAYS, 6.0, "light"
                ),  # light toast ⇒ most ellagitannin (protection)
                _add_oxygen(_FERMENT_DAYS, o2_dose),
            ]
        )
    ).run()
    unoaked = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, o2_dose)])
    ).run()
    assert oaked.success and unoaked.success

    # PROTECTION: oak lowers BOTH the browning index and the oxidative acetaldehyde.
    assert float(oaked.series("A420")[-1]) < float(unoaked.series("A420")[-1])
    assert float(oaked.series("acetaldehyde")[-1]) < float(unoaked.series("acetaldehyde")[-1])
    # PARTIAL, not total — the oaked wine still browns and still makes some acetaldehyde.
    assert float(oaked.series("A420")[-1]) > 0.0
    # Astringency readout (mg/L ellagitannin, IBU-exact): positive on the oaked wine, zero un-oaked.
    astr_oaked = astringency_series(oaked.as_trajectory())
    assert np.allclose(astr_oaked, np.asarray(oaked.series("ellagitannin"), dtype=float) * 1000.0)
    assert astr_oaked[-1] > 0.0
    assert float(unoaked.series("ellagitannin")[-1]) == 0.0  # un-oaked wine has no tannin


def test_oak_extraction_raises_the_oak_oavs():
    # The sensory payoff: oak aging lifts the four oak-extractive OAVs from 0. Read the finished
    # OAV series through the D-67 lens — each rises above 0 with oak, stays 0 without (the readout
    # is a pure consumer of the trajectory; the ceiling slots are NOT aroma pools, so unread).
    thresholds = load_thresholds()
    oaked = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 6.0, "medium")])
    ).run()
    plain = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    for compound in _OAK_EXTRACTIVES:
        assert float(oav_series(oaked.as_trajectory(), thresholds, compound)[-1]) > 0.0
        assert float(oav_series(plain.as_trajectory(), thresholds, compound)[-1]) == 0.0


def test_add_oak_rejects_unknown_toast_and_accepts_beer():
    # The vocabulary boundary (loud errors): an unknown toast is a typo, rejected loudly.
    with pytest.raises(ValueError, match="unknown toast"):
        compile_scenario(_wine([_add_oak(_FERMENT_DAYS, 4.0, "charred")]))
    # Barrel-beer oak (D-86): oak is no longer wine-only — a beer scenario now ACCEPTS add_oak
    # (bourbon-barrel stouts / oak-aged sours), setting the same 5 ceilings from oak_gpl × toast
    # yield. (Was a "wine-only" rejection before D-86 — a legitimate expectation flip, not a
    # weakened test.) The unknown-toast rejection above still holds for beer too.
    cs = compile_scenario(_beer([_add_oak(14.0, 4.0, "medium")]))
    event = next(e for e in cs.events if e.label.startswith("add_oak"))
    assert event.mutate is not None  # dose only (a pure ceiling mutate, like the wine add_oak)
    after = event.mutate(cs.schema, cs.y0.copy())
    # The medium-toast dose sets each ceiling to oak_gpl × oak_yield_<compound>_medium (> 0).
    for compound in ("whiskey_lactone", "vanillin", "guaiacol", "eugenol", "ellagitannin"):
        assert float(cs.schema.get(after, f"{compound}_ceiling")) > 0.0


# -- D-86: barrel-beer oak end-to-end (the same wood axis, on beer) ------------------------------


def test_un_oaked_beer_aging_leaves_the_oak_pools_zero_and_closes_ledgers():
    # Isolability on beer (D-86, the D-77 three-case pattern): a begin_aging beer with NO add_oak
    # leaves every oak extractive + ceiling identically 0 — byte-for-byte the case without the oak
    # axis — and the whole ferment+aging trajectory still closes carbon AND nitrogen (the oak slots
    # are off every ledger, so wiring them into beer perturbs no elemental balance).
    cs = compile_scenario(_beer([_begin_aging(14.0)]))
    traj = cs.run()
    assert traj.success
    for name in _OAK_EXTRACTIVES + _OAK_CEILINGS + ("ellagitannin", "ellagitannin_ceiling"):
        assert np.max(np.abs(traj.series(name))) == 0.0  # identically zero — no oak dosed
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    assert_conserved(
        traj.as_trajectory(), total_carbon(cs.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj.as_trajectory(), total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n),
        label="nitrogen",
    )


def test_barrel_beer_oak_raises_the_oak_oavs_and_astringency():
    # The sensory payoff on beer (bourbon-barrel stout / oak-aged sour): oak aging lifts the four
    # oak-extractive OAVs from 0 (read through the D-67 lens against the beer-matrix thresholds,
    # D-86) and the ellagitannin astringency readout goes positive. An un-oaked beer reads 0.
    thresholds = load_thresholds()
    oaked = compile_scenario(
        _beer([_begin_aging(14.0), _add_oak(14.0, 6.0, "medium")])
    ).run()
    plain = compile_scenario(_beer([_begin_aging(14.0)])).run()
    assert oaked.success and plain.success
    for compound in _OAK_EXTRACTIVES:
        assert float(oav_series(oaked.as_trajectory(), thresholds, compound)[-1]) > 0.0
        assert float(oav_series(plain.as_trajectory(), thresholds, compound)[-1]) == 0.0
    # Astringency = oak ellagitannin alone (beer has no grape tannin slot — D-86): positive oaked.
    astr = astringency_series(oaked.as_trajectory())
    assert np.allclose(astr, np.asarray(oaked.series("ellagitannin"), dtype=float) * 1000.0)
    assert astr[-1] > 0.0
    assert float(plain.series("ellagitannin")[-1]) == 0.0


def test_barrel_beer_oak_protects_against_oxidation():
    # The D-78 protection spine, on beer (D-86): two identical oxygenated aged beers — same O₂ dose,
    # same ferment+aging — differing ONLY in the add_oak charge. The oaked beer's ellagitannin
    # scavenges its share of the O₂, so it browns LESS (lower A420) and makes LESS oxidative
    # acetaldehyde. Beer already runs the always-on O₂ sinks (OxidativeAcetaldehyde/PhenolicBrowning
    # in the medium-agnostic _AGING_PROCESSES), so the substrate-gated ellag sink adds on top.
    o2_dose = 40.0
    oaked = compile_scenario(
        _beer([_begin_aging(14.0), _add_oak(14.0, 6.0, "light"), _add_oxygen(14.0, o2_dose)])
    ).run()  # light toast ⇒ most ellagitannin (strongest protection)
    unoaked = compile_scenario(
        _beer([_begin_aging(14.0), _add_oxygen(14.0, o2_dose)])
    ).run()
    assert oaked.success and unoaked.success
    assert float(oaked.series("A420")[-1]) < float(unoaked.series("A420")[-1])
    assert float(oaked.series("acetaldehyde")[-1]) < float(unoaked.series("acetaldehyde")[-1])
    assert float(oaked.series("A420")[-1]) > 0.0  # PARTIAL, not total — still some browning


# -- D-79: tannin–anthocyanin condensation end-to-end (red-wine softening + colour stabilization) --


def test_polymerization_disabled_without_begin_aging():
    # Wired into the wine medium (present in the set) but DISABLED at compile — aging is inherently
    # post-ferment. Even a red wine (grape pools dosed) never activates it until begin_aging.
    cs = compile_scenario(_wine([], anthocyanin_gpl=0.3, tannin_gpl=2.0))
    assert TanninAnthocyaninCondensation.name in cs.process_set
    assert not cs.process_set.is_enabled(TanninAnthocyaninCondensation.name)


def test_begin_aging_enables_polymerization_at_the_breakpoint():
    # The reconfigure enables it alongside the other aging Processes (its name rides in
    # _AGING_GATED_PROCESSES) — a pure phase switch, mutating no state.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], anthocyanin_gpl=0.3, tannin_gpl=2.0))
    assert not cs.process_set.is_enabled(TanninAnthocyaninCondensation.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.mutate is None and event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(TanninAnthocyaninCondensation.name)


def test_polymerization_params_ride_in_every_compiled_scenario():
    # polymerization.yaml is in shared_files, so every compiled scenario carries the condensation
    # params — inert (read by nothing until begin_aging enables the Process on a grape-dosed wine).
    cs = compile_scenario(_wine([]))
    for name in ("k_polymerization", "E_a_polymerization", "y_tannin_per_anthocyanin"):
        assert name in cs.param_values


def test_white_wine_polymerization_is_byte_for_byte_inert():
    # Doubly substrate-gated: a white wine (no anthocyanin/tannin dose) aged with begin_aging is
    # byte-for-byte the case without this Process — the grape pools stay 0 and all three readouts
    # are identically zero (isolability #3).
    white = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    assert white.success
    assert float(white.series("anthocyanin")[-1]) == 0.0
    assert float(white.series("tannin")[-1]) == 0.0
    white_traj = white.as_trajectory()
    assert np.all(astringency_series(white_traj) == 0.0)
    assert np.all(polymeric_pigment_series(white_traj) == 0.0)
    assert np.all(color_series(white_traj) == 0.0)


def test_red_wine_softens_and_stabilizes_colour_end_to_end():
    # The D-79 spine through the compiled run() path (a red wine: both grape must inputs dosed,
    # begin_aging enabling the Process over the warm aging tail):
    #   (1) astringency SOFTENS — free tannin is drawn down from its young (post-ferment) value;
    #   (2) polymeric pigment RISES from 0 (the stable pigment = anthocyanin condensed);
    # (3) total colour is RETAINED (free anthocyanin declines, stable pigment rises, sum constant);
    #   (4) OAK-INDEPENDENT — no add_oak anywhere, yet it softens and stabilizes (a steel-tank red).
    red = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)], anthocyanin_gpl=0.3, tannin_gpl=2.0, aging_celsius=25.0
        )
    ).run()
    assert red.success
    red_traj = red.as_trajectory()

    astr = astringency_series(red_traj)
    # No oak, so astringency is exactly grape tannin × 1000; it softens over the aging tail.
    assert np.allclose(astr, np.asarray(red.series("tannin"), dtype=float) * 1000.0)
    assert astr[-1] < astr[0]

    pig = polymeric_pigment_series(red_traj)
    assert pig[0] == pytest.approx(0.0)
    assert pig[-1] > 0.0  # stable pigment formed

    # Colour largely RETAINED as it shifts form (free anthocyanin → stable pigment), MINUS a small
    # O₂-independent THERMAL fade (D-83: even this anaerobic red loses a few mg/L to the colourless
    # sink over the warm tail). The three-slot identity anthocyanin + pigment + faded ≡ antho0 holds
    # by construction, so color_series = (antho0 − faded) × 1000 — the stabilization physics
    # (pigment holds the colour condensation shifts) modulo the thermal loss; anthocyanin drawdown +
    # pigment rise are the load-bearing Process signals.
    col = color_series(red_traj)
    faded = np.asarray(red.series("faded_anthocyanin"), dtype=float)
    assert float(red.series("anthocyanin")[-1]) < 0.3  # the genuine dynamic
    assert 0.0 < faded[-1] < 0.01  # a SMALL thermal-only fade (anaerobic ⇒ no oxidative bleaching)
    assert np.allclose(col, (0.3 - faded) * 1000.0)  # retained modulo thermal fade (three-slot id.)


def test_red_wine_polymerization_off_ledger_end_to_end():
    # Both grape pools are off every ledger (grape-derived), so condensation moves nothing
    # conserved:
    # total_carbon and total_nitrogen close across the whole ferment+aging run (a pure reconfigure,
    # no external flow), exactly as the un-dosed aged run does (the OakExtraction precedent).
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], anthocyanin_gpl=0.3, tannin_gpl=2.0))
    traj = cs.run()
    assert traj.success
    assert traj.external_flows == ()  # grape must inputs are t0 initial conditions, not flows
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    schema = cs.schema
    assert_conserved(
        traj.as_trajectory(), total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj.as_trajectory(),
        total_nitrogen(schema, biomass_nitrogen_fraction=f_n),
        label="nitrogen",
    )
    assert_nonnegative(traj.as_trajectory(), ("anthocyanin", "tannin"), atol=1e-9)


def test_white_wine_tannin_softens_by_self_polymerization_end_to_end():
    # THE D-84 PAYOFF at scenario scale, retiring the D-80 "softening needs anthocyanin" note. A
    # tannin-dosed WHITE wine (no anthocyanin ⇒ both condensation routes inert) still SOFTENS over
    # the aging tail, purely by TanninSelfPolymerization drawing the free-tannin pool down. Warmer
    # softens more (the E_a lever), and it is oak-free (a steel-tank white). Colour stays zero
    # throughout (no anthocyanin), confirming the softener is decoupled from the colour axis.
    white_args = {"tannin_gpl": 3.0}  # tannin but NO anthocyanin
    cool = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS)], aging_celsius=12.0, **white_args)
    ).run()
    warm = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS)], aging_celsius=30.0, **white_args)
    ).run()
    assert cool.success and warm.success
    astr_cool = astringency_series(cool.as_trajectory())
    astr_warm = astringency_series(warm.as_trajectory())
    # Softens over aging (no anthocyanin, so this is self-polymerization ONLY — the retirement).
    assert astr_cool[-1] < astr_cool[0]
    assert astr_warm[-1] < astr_warm[0]
    # Warmer softens MORE (E_a > 0, the temperature lever).
    assert astr_warm[-1] < astr_cool[-1]
    # No anthocyanin ⇒ colour is identically zero throughout (softener decoupled from colour).
    assert np.all(color_series(warm.as_trajectory()) == 0.0)


# -- D-80: acetaldehyde-bridged condensation end-to-end (micro-oxygenation → colour, split ledger)
# --


def test_bridged_condensation_disabled_without_begin_aging():
    # Wired into the wine medium (present in the set) but DISABLED at compile — aging is
    # post-ferment.
    # Even a red wine (grape pools dosed) never activates the bridged route until begin_aging.
    cs = compile_scenario(_wine([], anthocyanin_gpl=0.3, tannin_gpl=2.0))
    assert AcetaldehydeBridgedCondensation.name in cs.process_set
    assert not cs.process_set.is_enabled(AcetaldehydeBridgedCondensation.name)


def test_begin_aging_enables_bridged_condensation_at_the_breakpoint():
    # The reconfigure enables it alongside the other aging Processes (its name rides in
    # _AGING_GATED_PROCESSES) — a pure phase switch, mutating no state.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], anthocyanin_gpl=0.3, tannin_gpl=2.0))
    assert not cs.process_set.is_enabled(AcetaldehydeBridgedCondensation.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.mutate is None and event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(AcetaldehydeBridgedCondensation.name)


def test_bridged_params_ride_in_every_compiled_scenario():
    # polymerization.yaml (D-80's params live with the D-79 ones) is in shared_files, so every
    # compiled scenario carries the bridged-route params — inert until begin_aging on a grape-dosed,
    # oxygenated wine.
    cs = compile_scenario(_wine([]))
    for name in (
        "k_acetaldehyde_bridge",
        "E_a_acetaldehyde_bridge",
        "y_acetaldehyde_per_anthocyanin",
    ):
        assert name in cs.param_values


def test_micro_oxygenation_drives_bridged_condensation_end_to_end():
    # THE D-80 SPINE — the first link from the oxidative sub-axis to red-wine colour, through the
    # compiled run(). ``ethyl_bridge`` (the acetaldehyde carbon captured on-ledger) is a PURE
    # micro-ox
    # signal: after fermentation the viable yeast has cleared acetaldehyde to ~0, so an anaerobic
    # aged
    # red bridges nothing (ethyl_bridge ≡ 0); dosing O₂ (add_oxygen) makes OxidativeAcetaldehyde
    # (D-71)
    # regenerate acetaldehyde, which the bridged route then consumes → ethyl_bridge accumulates. And
    # SO₂ DELAYS it (bound acetaldehyde can't bridge, D-47/D-80; plus SulfiteOxidation scavenges the
    # O₂, D-72) — so a sulfited oxygenated red bridges LESS than an unsulfited one. (Anthocyanin,
    # the
    # limiting reagent, fully condenses via the direct route in all runs, so the pigment ENDPOINT
    # saturates and is not the discriminator here — ``ethyl_bridge`` is.)
    red_args = {"anthocyanin_gpl": 0.3, "tannin_gpl": 2.0}
    noox = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], **red_args)).run()
    ox = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 40.0)], **red_args)
    ).run()
    ox_so2 = compile_scenario(
        _wine(
            [
                _begin_aging(_FERMENT_DAYS),
                _add_oxygen(_FERMENT_DAYS, 40.0),
                _add_so2(_FERMENT_DAYS, 100.0),
            ],
            **red_args,
        )
    ).run()
    assert noox.success and ox.success and ox_so2.success

    bridge_noox = float(np.asarray(noox.series("ethyl_bridge"), dtype=float)[-1])
    bridge_ox = float(np.asarray(ox.series("ethyl_bridge"), dtype=float)[-1])
    bridge_ox_so2 = float(np.asarray(ox_so2.series("ethyl_bridge"), dtype=float)[-1])
    # Anaerobic aged red: no acetaldehyde survives to aging ⇒ NO bridging (exactly zero).
    assert bridge_noox == 0.0
    # Micro-ox regenerates acetaldehyde ⇒ the bridge accumulates (the O₂ → colour link).
    assert bridge_ox > 1e-4
    # SO₂ delays it: the sulfited oxygenated red bridges strictly less (emergent, free-acetaldehyde
    # +
    # O₂-scavenging), but still some (partial protection).
    assert 0.0 < bridge_ox_so2 < bridge_ox

    # The direct route still (near-)exhausts anthocyanin in every run (colour endpoint saturates) —
    # so the bridged route is an ADDITIONAL pigment-formation pathway, not the only one. NOTE (D-84
    # coupling): TanninSelfPolymerization now competes for the tannin pool, drawing [tannin] down a
    # touch, so the anthocyanin-condensation rate (∝ [tannin]) is slightly lower and a tiny
    # anthocyanin residual survives (~1e-4 g/L, still >99.9 % consumed) — hence < 1e-3, not < 1e-6.
    assert float(ox.series("anthocyanin")[-1]) < 1e-3
    # Micro-ox sustains a residual acetaldehyde the anaerobic run lacks (O₂ keeps regenerating it).
    assert float(ox.series("acetaldehyde")[-1]) > float(noox.series("acetaldehyde")[-1])


def test_micro_oxygenation_now_fades_colour_end_to_end():
    # THE D-81 PAYOFF at scenario scale (that RETIRED the D-80 "colour O₂-invariant" pin): an
    # OXYGENATED red FADES more than an anaerobic one — dissolved O₂ adds an oxidative bleaching
    # on top. NOTE (D-83 retirement): the anaerobic red no longer holds PERFECTLY flat — it now
    # fades a SMALL THERMAL amount (ThermalAnthocyaninFade, O₂-independent), so this test pins the
    # O₂ CONTRAST (ox fades strictly MORE than no-ox) rather than "no-ox holds at antho₀ × 1000"
    # (that D-81 assertion is retired by D-83; see test_thermal_fade_adds_to_oxidative_fade for the
    # thermal decomposition). Both runs still satisfy the three-slot colour identity.
    red_args = {"anthocyanin_gpl": 0.3, "tannin_gpl": 2.0}
    noox = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], **red_args)).run()
    ox = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 40.0)], **red_args)
    ).run()
    assert noox.success and ox.success
    col_noox = color_series(noox.as_trajectory())
    col_ox = color_series(ox.as_trajectory())
    # Anaerobic red fades only THERMALLY (D-83) — a small background loss, well under the O₂ fade.
    faded_noox = float(np.asarray(noox.series("faded_anthocyanin"), dtype=float)[-1])
    assert 0.0 < faded_noox < 0.01  # small thermal-only fade (≈ a few mg/L over 150 d at 25 °C)
    # Oxygenated red FADES MORE: O₂ adds an oxidative bleaching sink on top of the thermal one.
    assert col_ox[-1] < col_noox[-1] - 5.0
    # Each run's lost colour is exactly its faded (colourless) fraction (the three-slot identity).
    faded_ox = float(np.asarray(ox.series("faded_anthocyanin"), dtype=float)[-1])
    assert faded_ox > faded_noox  # O₂ fade strictly exceeds thermal-only
    assert col_ox[-1] == pytest.approx((0.3 - faded_ox) * 1000.0, rel=1e-3)
    assert col_noox[-1] == pytest.approx((0.3 - faded_noox) * 1000.0, rel=1e-3)


def test_thermal_fade_adds_to_oxidative_fade_end_to_end():
    # THE D-83 PAYOFF at scenario scale: ThermalAnthocyaninFade fades a red even with NO O₂, and its
    # thermal loss is TEMPERATURE-driven (warm storage bleaches faster) and SO₂-UNPROTECTED — the
    # mirror of D-81's O₂-coupled, SO₂-protected fade. Two anaerobic reds, cellar-cool vs warm: the
    # warm one fades MORE, purely thermally (no O₂ in either), and a heavy SO₂ dose does NOT rescue
    # the warm red's colour (SO₂ protects only the O₂ route, D-81 — it can't touch a thermal one).
    red_args = {"anthocyanin_gpl": 0.3, "tannin_gpl": 2.0}
    cool = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS)], aging_celsius=12.0, **red_args)
    ).run()
    warm = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS)], aging_celsius=30.0, **red_args)
    ).run()
    warm_so2 = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_so2(_FERMENT_DAYS, 150.0)],
            aging_celsius=30.0,
            **red_args,
        )
    ).run()
    assert cool.success and warm.success and warm_so2.success
    faded_cool = float(np.asarray(cool.series("faded_anthocyanin"), dtype=float)[-1])
    faded_warm = float(np.asarray(warm.series("faded_anthocyanin"), dtype=float)[-1])
    faded_warm_so2 = float(np.asarray(warm_so2.series("faded_anthocyanin"), dtype=float)[-1])
    # Anaerobic reds still fade (thermally): both lose SOME colour with zero O₂ (the D-83 headline).
    assert faded_cool > 0.0 and faded_warm > 0.0
    # Warmer fades MORE (E_a > 0, the temperature lever) — 'warm cellars kill colour' even sealed.
    assert faded_warm > faded_cool
    # SO₂ does NOT protect the thermal route (the D-83 mirror of D-81): the sulfited warm red fades
    # essentially the SAME as the unsulfited warm red (no O₂ to scavenge ⇒ no emergent protection).
    assert faded_warm_so2 == pytest.approx(faded_warm, rel=1e-3)


# == D-82: the reversible SO₂/pH masking readout (observed_color_series), end to end ============
# The committed second half of D-81's "Both" fork: observed_color_series masks the free monomeric
# anthocyanin by χ(SO₂, pH) while counting the SO₂/pH-resistant polymeric pigment full — the Somers
# assay, distinct from color_series (intrinsic pigment CONTENT). A pure readout: no state slot.


def test_observed_colour_white_wine_is_zero():
    # Wine-only, doubly gated: a white wine (no anthocyanin/tannin dose) reads identically zero —
    # no monomeric to mask, no pigment to count (the color_series/astringency isolability contract).
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    white = cs.run()
    assert white.success
    assert np.all(observed_color_series(white.as_trajectory(), cs.param_values) == 0.0)


def test_observed_colour_is_masked_below_content_colour():
    # observed ≤ color_series ALWAYS (χ ≤ 1; the pigment term is identical in both). And it is
    # strictly BELOW at wine pH even with no SO₂, because only a minority of monomeric anthocyanin
    # is red (the flavylium fraction ~0.14 at pH 3.4) — the pH mask alone.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)], anthocyanin_gpl=0.3, tannin_gpl=2.0, aging_celsius=25.0
        )
    )
    red = cs.run()
    assert red.success
    traj = red.as_trajectory()
    obs = observed_color_series(traj, cs.param_values)
    col = color_series(traj)
    assert np.all(obs <= col + 1e-9)  # masked ≤ content, at every column
    assert obs[0] < col[0]  # strictly masked while colour is still monomeric (pH mask, no SO₂)


def test_so2_masks_observed_colour_opposite_sign_to_fade():
    # THE D-82 HEADLINE and the OPPOSITE-SIGN point vs D-81. Two reds identical but for an SO₂ dose,
    # aged ANAEROBICALLY (no O₂ ⇒ no fade ⇒ SO₂ is inert to color_series; SulfiteOxidation needs
    # O₂): color_series is IDENTICAL between them (content held fixed), yet observed_color_series is
    # LOWER for the sulfited wine — SO₂ reversibly BLEACHES the monomeric colour. This is the
    # opposite sign to D-81 (where SO₂ PROTECTS color_series by scavenging fade-driving O₂): a
    # different series and a different (reversible) mechanism, both real — do not "reconcile" them.
    red_args = {"anthocyanin_gpl": 0.3, "tannin_gpl": 2.0}
    plain_cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], **red_args))
    so2_cs = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_so2(_FERMENT_DAYS, 150.0)], **red_args)
    )
    plain, sulfited = plain_cs.run(), so2_cs.run()
    assert plain.success and sulfited.success
    plain_traj, sulf_traj = plain.as_trajectory(), sulfited.as_trajectory()
    # Content colour is identical (no O₂ ⇒ no fade; SO₂ does not touch anthocyanin/pigment content).
    assert np.allclose(color_series(plain_traj), color_series(sulf_traj))
    # …yet the sulfited wine shows LESS observed colour throughout — the bleaching mask.
    obs_plain = observed_color_series(plain_traj, plain_cs.param_values)
    obs_sulf = observed_color_series(sulf_traj, so2_cs.param_values)
    assert obs_sulf[-1] < obs_plain[-1]
    assert np.all(obs_sulf <= obs_plain + 1e-9)


def test_condensation_unmasks_observed_colour_while_content_flat():
    # The Somers "ageing shifts colour onto the SO₂/pH-resistant pigment" evolution: as monomeric
    # anthocyanin condenses to bleach-/pH-resistant polymeric pigment (counted FULL), observed
    # colour RISES over the aging tail — even though color_series (content) is NEARLY flat (this red
    # is anaerobic, so no oxidative fade; condensation conserves content and only a SMALL D-83
    # thermal fade nibbles it). observed_color_series RISES while color_series is flat-to-slightly-
    # declining — they trend OPPOSITELY here: the reason beat A was worth building alongside beat B.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)], anthocyanin_gpl=0.3, tannin_gpl=2.0, aging_celsius=25.0
        )
    )
    red = cs.run()
    assert red.success
    traj = red.as_trajectory()
    col = color_series(traj)
    obs = observed_color_series(traj, cs.param_values)
    # Content NEARLY flat: only a small thermal fade (anaerobic ⇒ no oxidative fade), well under
    # the strong observed rise below — the load-bearing contrast is the OPPOSITE trends.
    assert col[-1] > 0.3 * 1000.0 - 5.0  # ≈ flat (small thermal loss only)
    assert polymeric_pigment_series(traj)[-1] > 0.0  # pigment genuinely formed
    assert obs[-1] > obs[0] + 5.0  # observed RISES strongly: masked monomeric → unmasked pigment


def test_bridged_run_closes_carbon_end_to_end():
    # THE SPLIT-LEDGER CONSERVATION PROOF at scenario scale. The bridged route consumes ON-ledger
    # acetaldehyde (borrowed from E at D-71) and books its carbon into the on-ledger ethyl_bridge
    # slot, so total_carbon closes across the WHOLE E → acetaldehyde → ethyl_bridge chain (the D-71
    # borrow + the D-80 capture, both carbon-exact). Nitrogen closes too. total_mass is NOT
    # asserted:
    # the D-71 E → acetaldehyde transfer moves mass into an unweighted pool (the standing aging-axis
    # scope-out), and add_oxygen books an external o2 flow besides — carbon (o2 off-ledger) is the
    # invariant.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 40.0)],
            anthocyanin_gpl=0.3,
            tannin_gpl=2.0,
        )
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    schema = cs.schema
    assert_conserved(
        traj.as_trajectory(), total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj.as_trajectory(),
        total_nitrogen(schema, biomass_nitrogen_fraction=f_n),
        label="nitrogen",
    )
    # The ethyl_bridge slot genuinely accumulated (the closure is non-trivial, not vacuous).
    assert float(np.asarray(traj.series("ethyl_bridge"), dtype=float)[-1]) > 1e-4
    assert_nonnegative(traj.as_trajectory(), ("anthocyanin", "tannin", "ethyl_bridge"), atol=1e-9)


def test_micro_oxygenation_softens_white_tannin_via_ethyl_bridge_end_to_end():
    # THE D-85 PAYOFF at scenario scale, and the split-ledger proof for the tannin–tannin route. A
    # tannin-dosed WHITE wine (NO anthocyanin ⇒ both D-79/D-80 colour routes AND the anthocyanin
    # fade are inert) is micro-oxygenated: dissolved-O₂ acetaldehyde bridges two flavanols
    # (tannin–ethyl–tannin), so astringency softens MORE than the anaerobic control (where only the
    # D-84 direct self-polymerization acts). The bridge carbon is booked on-ledger (ethyl_bridge
    # accumulates), total_carbon closes across the E → acetaldehyde → ethyl_bridge chain, and colour
    # stays identically zero throughout (a colourless tannin–tannin polymer — the D-80 colour
    # difference).
    white_args = {"tannin_gpl": 3.0}  # tannin but NO anthocyanin
    noox = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], **white_args)).run()
    ox = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 40.0)], **white_args)
    ).run()
    assert noox.success and ox.success
    astr_noox = astringency_series(noox.as_trajectory())
    astr_ox = astringency_series(ox.as_trajectory())
    # Micro-ox softens MORE than the anaerobic control (the bridged route adds to the D-84 direct).
    assert astr_ox[-1] < astr_noox[-1] - 1.0
    # The bridge carbon accumulated on-ledger (the tannin–ethyl–tannin split-ledger capture).
    bridge_ox = float(np.asarray(ox.series("ethyl_bridge"), dtype=float)[-1])
    bridge_noox = float(np.asarray(noox.series("ethyl_bridge"), dtype=float)[-1])
    assert bridge_ox > 1e-4  # micro-ox regenerates acetaldehyde that bridges
    assert bridge_noox == pytest.approx(0.0, abs=1e-9)  # anaerobic: no acetaldehyde to bridge
    # Colour stays identically zero (no anthocyanin — the tannin–ethyl–tannin polymer is colourless)
    assert np.all(color_series(ox.as_trajectory()) == 0.0)
    # total_carbon closes across the whole E → acetaldehyde → ethyl_bridge chain (o2 off-ledger).
    cs = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 40.0)], **white_args)
    )
    traj = cs.run()
    f_c = cs.parameters.value("biomass_C_fraction")
    assert_conserved(
        traj.as_trajectory(), total_carbon(cs.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
