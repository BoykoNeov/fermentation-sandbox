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
  aging segment, and every other producer of ``esters``/``isoamyl_alcohol``/``Byp`` is
  flux-gated and
  silent there — so the aging signal is unconfounded (Stance A).

These tests pin: the compile-seam enable/disable gate; the emergent aging headline (esters fade,
isoamyl alcohol + Byp rise) end-to-end vs an otherwise-identical un-aged run; end-to-end
carbon closure
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
from fermentation.core.kinetics import FuselAminoAcidReroute
from fermentation.core.kinetics.aging import (
    AcetaldehydeBridgedCondensation,
    Caramelization,
    EllagitanninOxidation,
    EsterHydrolysis,
    MaillardStrecker,
    OakExtraction,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    StreckerDegradation,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
)
from fermentation.core.kinetics.amino_acid_pools import AMINO_ACID_SPECS
from fermentation.core.media import get_medium
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir
from fermentation.runtime import ScheduledTrajectory
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
# aging effect is the ONLY thing moving esters/isoamyl_alcohol/Byp (Stance A). Warm (25 C) aging + a
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
    brix: float = 24.0,
    duration_days: float | None = None,
) -> Scenario:
    # brix (default 24.0, byte-for-byte the pre-D-87 dry-fermenting helper) sets the must sugar; a
    # botrytis-level brix (~70) ferments to the ethanol-inactivation ceiling and ARRESTS with a
    # large
    # RESIDUAL sugar — a SWEET wine — the driver the non-oxidative thermal MaillardStrecker (D-87)
    # needs. (The modelled ABV at that brix runs high — a pre-existing EthanolInactivation
    # calibration limit, orthogonal to D-87: what matters here is the residual sugar it leaves.)
    # duration_days overrides the default ferment+aging span (a longer sweet-wine aging tail).
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
    initial: dict[str, float] = {"brix": brix, "yan_mgl": 250.0, "pitch_gpl": 0.25}
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
        duration_days=(_FERMENT_DAYS + _AGING_DAYS if duration_days is None else duration_days),
    )


def _beer(interventions: list[Intervention]) -> Scenario:
    # A minimal ~1.048-OG wort (glucose+maltose+maltotriose), fermented ~14 d then aged. Aging is
    # medium-agnostic (esters/isoamyl_alcohol/Byp exist in the beer schema too), so
    # begin_aging must drive
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


def _run_thermal_isolating_reroute(scenario: Scenario) -> ScheduledTrajectory:
    """Compile + run with the D-33 Ehrlich amino-acid reroute disabled (decision D-99).

    Any end-to-end thermal test that asserts a MAGNITUDE — thermal aldehydes at aroma-relevant
    levels, or "the route fired" — must isolate the reroute, or the assertion goes vacuous. The
    reroute and MaillardStrecker both draw the lumped `amino_acids` pool; since D-99 raised fusel
    production ~3.8x the reroute empties it during FERMENTATION (it needs flux, so it is inactive
    by the time begin_aging fires), leaving the aging-phase thermal routes ~nothing to work on.
    Isolating it restores the substrate WITHOUT weakening aging-phase conservation — the reroute
    is not active during aging anyway. The starvation itself is a real known limitation (D-100),
    pinned in test_the_ehrlich_reroute_starves_maillard_of_the_lumped_amino_acid_pool. A PURE
    conservation test (no magnitude claim) deliberately does NOT use this helper — closure must
    hold with every Process on, and the reroute makes that stress stronger, not weaker.
    """
    cs = compile_scenario(scenario)
    cs.process_set.disable(FuselAminoAcidReroute.name)
    return cs.run()


def _add_oxygen(day: float, o2_mgl: float) -> Intervention:
    return Intervention(day=day, action="add_oxygen", params={"o2_mgl": o2_mgl})


def _add_oak(
    day: float,
    oak_gpl: float,
    toast: str,
    fill_number: float | None = None,
    spirit: str | None = None,
) -> Intervention:
    params: dict[str, object] = {"oak_gpl": oak_gpl, "toast": toast}
    if fill_number is not None:  # D-91 barrel fill-number depletion (default omitted ⇒ fresh fill)
        params["fill_number"] = fill_number
    if spirit is not None:  # D-92 ex-spirit barrel soak-back (default omitted ⇒ no soak-back)
        params["spirit"] = spirit
    return Intervention(day=day, action="add_oak", params=params)


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
    for name in ("k_ester_hydrolysis", "E_a_ester_hydrolysis", "isoamyl_acetate_eq"):
        assert name in cs.parameters


# -- the emergent aging headline (esters fade, isoamyl alcohol + Byp rise) ----


def test_aging_fades_esters_and_raises_isoamyl_alcohol_end_to_end():
    # The end-to-end payoff: an aged wine finishes with LOWER esters and HIGHER
    # isoamyl_alcohol/Byp than
    # the otherwise-identical un-aged wine. Both runs share the identical ferment (aging is off
    # until the breakpoint), so any difference at the end is exactly the aging Process's doing —
    # a clean A/B that also proves isolability (the un-aged run never activates the Process).
    aged = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)])).run()
    plain = compile_scenario(_wine([])).run()
    assert aged.success and plain.success

    ester_aged = float(aged.series("isoamyl_acetate")[-1])
    ester_plain = float(plain.series("isoamyl_acetate")[-1])
    isoamyl_aged = float(aged.series("isoamyl_alcohol")[-1])
    isoamyl_plain = float(plain.series("isoamyl_alcohol")[-1])
    byp_aged = float(aged.series("Byp")[-1])
    byp_plain = float(plain.series("Byp")[-1])

    # The wine actually made ester during the ferment, so there is something to hydrolyse.
    assert ester_plain > 0.0
    # Aging hydrolyses the banana acetate ester: less ester, more isoamyl alcohol + Byp at the
    # end. Since D-99 the alcohol lands in the isoamyl pool SPECIFICALLY — hydrolysing isoamyl
    # acetate yields 3-methylbutan-1-ol and nothing else, so the old lump silently credited a
    # share of it to four other molecules.
    assert ester_aged < ester_plain
    assert isoamyl_aged > isoamyl_plain
    assert byp_aged > byp_plain


def test_aging_does_not_strip_esters_below_the_equilibrium_floor():
    # Net decay toward a LOWER floor, not decay-to-zero (D-68): even over a long warm aging tail
    # the esters pool relaxes toward, not past, isoamyl_acetate_eq.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    traj = cs.run()
    assert traj.success
    assert float(traj.series("isoamyl_acetate")[-1]) >= cs.param_values["isoamyl_acetate_eq"]


# -- conservation + non-negativity end-to-end --------------------------------


def test_aged_run_closes_carbon_end_to_end():
    # begin_aging mutates no state (a pure reconfigure), so there is NO external flow: the
    # run-wide invariant is the plain final == initial carbon. The ferment routes carbon into the
    # aroma pools and the aging segment transfers esters -> isoamyl_alcohol + Byp — both
    # close to machine
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
    assert_nonnegative(
        traj.as_trajectory(), ("isoamyl_acetate", "isoamyl_alcohol", "Byp"), atol=1e-9
    )


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
    # segment, so the speculative aging Process drags esters/isoamyl_alcohol/Byp to
    # speculative for the
    # WHOLE run — not just the aging segment.
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)]))
    traj = cs.run()
    for pool in ("isoamyl_acetate", "isoamyl_alcohol", "Byp"):
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
    # Aging is medium-agnostic (D-70): esters/isoamyl_alcohol/Byp exist in the beer schema,
    # so begin_aging
    # must compile, enable EsterHydrolysis, and fade the beer esters exactly as it does for wine.
    # Beer-side smoke coverage — the D-69 Process math is beer-tested, but the beer scenario
    # compile -> begin_aging -> run path was otherwise unexercised.
    cs = compile_scenario(_beer([_begin_aging(14.0)]))
    assert not cs.process_set.is_enabled(EsterHydrolysis.name)  # disabled at compile
    aged = cs.run()
    plain = compile_scenario(_beer([])).run()
    assert aged.success and plain.success
    ester_plain = float(plain.series("isoamyl_acetate")[-1])
    assert ester_plain > 0.0  # the beer ferment made ester to hydrolyse
    assert float(aged.series("isoamyl_acetate")[-1]) < ester_plain  # aging fades it
    assert float(aged.series("isoamyl_alcohol")[-1]) > float(plain.series("isoamyl_alcohol")[-1])


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
    #
    # D-90 supersession: the reductive (no-O₂) beer is no longer byte-for-byte A420 = 0. Caramel
    # ization is now medium-agnostic, so the beer's residual sugar (this scenario finishes at
    # S ≈ 5e-11 g/L — near-dry but NOT the exact ≤ 0 the wine scenarios hit, so the soft-sugar guard
    # does not fire) browns a NEGLIGIBLE thermal trace (~4e-8) over the warm aging tail. The
    # discriminating physics still holds: the O₂-driven PhenolicBrowning dominates by ~7 orders of
    # magnitude. (A meaningfully-brown reductive beer needs real residual sugar — see the
    # high-residual carbon-closure test in test_aging.py.)
    day = 14.0
    oxidative = compile_scenario(_beer([_begin_aging(day), _add_oxygen(day, 40.0)])).run()
    reductive = compile_scenario(_beer([_begin_aging(day)])).run()
    assert oxidative.success and reductive.success
    ox_a420 = float(oxidative.series("A420")[-1])
    red_a420 = float(reductive.series("A420")[-1])
    assert ox_a420 > 1e-2  # O₂-driven browning is substantial
    assert red_a420 < 1e-6  # reductive near-dry beer browns only a negligible thermal trace
    assert ox_a420 > 1e4 * red_a420  # oxidative browning dominates the caramelization trace


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
    #
    # NOW ON LEES (decision D-100). The phenyl-dominant ordering this test anchors is a property of
    # a REAL aged wine, and a real aged wine has had lees contact. Speciation revealed that the
    # ordering cannot hold without it: phenylalanine is 2-phenylethanol's Ehrlich precursor, so
    # fermentation strips it, while methionine (which no fusel eats) survives — leaving methional
    # dominant and phenylacetaldehyde silent, the inverse of the literature. Autolysis restores
    # phenylalanine and with it the ordering. See
    # test_thermal_aroma_from_drained_precursors_requires_autolysis for the mechanism, pinned.
    #
    # The lumped pool hid this: `f_methional` was standing in for precursor ABUNDANCE (its own
    # comment says so), which D-100 now models explicitly — so the old scenario only produced
    # phenylacetaldehyde by drawing on arginine, a molecule that makes none.
    o2_dose = 60.0
    aged = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, o2_dose)],
            amino_acids_gpl=0.5,
            autolysis_rate_per_h=1.0e-3,
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


# -- MaillardStrecker (decision D-87) — the NON-oxidative THERMAL Strecker axis, end to end -----
#
# The O₂-INDEPENDENT thermal mirror of StreckerDegradation (D-75): a SEALED (no add_oxygen) SWEET
# (residual sugar) amino-acid-dosed wine develops the sweet-wine/Madeira aldehyde suite from sugar +
# heat alone — where the O₂-only D-75 route, with no oxygen, makes NOTHING. A botrytis-level brix
# (~70) ferments to the ethanol-inactivation ceiling and arrests with a large residual sugar (the
# driver). These pin: the compile-seam gate (wine-only, disabled→begin_aging); the sealed-sweet
# aldehyde production through the FULL pipeline vs a dry-wine control (the soft sugar driver); and
# end-to-end carbon + nitrogen closure.

_MAILLARD_ALDEHYDES = (
    "methional",
    "phenylacetaldehyde",
    "2_methylbutanal",
    "3_methylbutanal",
    "2_methylpropanal",
    "sotolon",
)
_SWEET_BRIX = 70.0  # botrytis-level must → arrests with ~130 g/L residual sugar (a SWEET wine)
_SWEET_AGING_DAYS = 730.0  # a multi-year sweet-wine aging tail (thermal aging is slow)


def test_maillard_gated_by_begin_aging_wine_only():
    # MaillardStrecker is WINE-ONLY (reads wine-only amino_acids + deaminates to N) — present in the
    # wine set, absent from beer — and rides the aging gate: disabled at compile, then on by
    # begin_aging (the StreckerDegradation pattern).
    assert MaillardStrecker.name in get_medium("wine").build_process_set()
    assert MaillardStrecker.name not in get_medium("beer").build_process_set()
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], amino_acids_gpl=0.5))
    assert MaillardStrecker.name in cs.process_set
    assert not cs.process_set.is_enabled(MaillardStrecker.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(MaillardStrecker.name)  # begin_aging turns it on


def test_maillard_thermal_aldehydes_in_sealed_sweet_wine():
    # THE discriminating end-to-end through the FULL pipeline: a SEALED (NO add_oxygen) SWEET wine,
    # amino-acid-dosed, aged warm — develops ALL SIX thermal products from residual sugar + heat
    # with
    # no O₂ whatsoever. Compared against a DRY wine (brix 24 → ~0 residual sugar) under the
    # identical
    # sealed aging, which makes ~none (the soft sugar driver). This is the sugar-driven, O₂-free
    # route the D-75 oxidative one could never produce (no O₂ ⇒ D-75 is silent, per
    # test_strecker_silent_reductive above). Levels land aroma-relevant (sotolon ~µg/L, the
    # Sauternes/Madeira range) — anchored in thermal.yaml provenance, asserted directionally here.
    # The D-33 reroute is ISOLATED (D-99): it would drain the shared amino_acids pool to ~0 and
    # collapse these levels to noise, making the assertions below vacuous — see the helper.
    sweet = _run_thermal_isolating_reroute(
        _wine(
            [_begin_aging(_FERMENT_DAYS)],  # SEALED — no add_oxygen
            amino_acids_gpl=0.8,
            brix=_SWEET_BRIX,
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    )
    dry = _run_thermal_isolating_reroute(
        _wine(
            [_begin_aging(_FERMENT_DAYS)],
            amino_acids_gpl=0.8,
            brix=24.0,  # dry-fermenting control → residual sugar ≈ 0 ⇒ thermal route ~silent
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    )
    assert sweet.success and dry.success
    # The sweet wine really is sweet: a large residual sugar drives the route (no O₂ dosed
    # anywhere).
    assert float(sweet.series("S")[-1]) > 50.0
    assert float(dry.series("S")[-1]) < 1.0
    for pool in _MAILLARD_ALDEHYDES:
        sweet_end = float(sweet.series(pool)[-1])
        dry_end = float(dry.series(pool)[-1])
        assert sweet_end > 0.0  # the sealed sweet wine forms every thermal product
        assert sweet_end > 100.0 * dry_end  # orders more than the ~silent dry control


def test_maillard_raises_the_thermal_oavs():
    # Through the STATED acceptance lens (the D-67 OAV lens): the whole point of D-87 is the four
    # new
    # thermal aromas the lens now reads (+ the two shared with D-75). A sealed sweet wine raises
    # them
    # positive; a dry wine leaves them ~0. Sotolon (the curry/maple furanone) is the diagnostic one.
    #
    # THE D-33 REROUTE IS ISOLATED OUT (decision D-99), and the assertions below are UNCHANGED in
    # strength — this tests D-87 without an unrelated, known-limited D-33 interaction. Both
    # Processes draw on the ONE lumped `amino_acids` pool (arginine standing in for every amino
    # acid), and with honest D-99 fusel levels the Ehrlich reroute consumes all of it before
    # MaillardStrecker can, driving sotolon to zero. That interaction is a REAL known limitation
    # and is pinned by its own test below — see test_the_ehrlich_reroute_starves_maillard_of_the
    # _lumped_amino_acid_pool and D-100. Isolating it here is the `_isolate_*` pattern this suite
    # already uses; the pathology is captured in CI, not hidden from it.
    thresholds = load_thresholds()
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)],
            amino_acids_gpl=0.8,
            brix=_SWEET_BRIX,
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    )
    cs.process_set.disable(FuselAminoAcidReroute.name)
    sweet = cs.run()
    assert sweet.success
    for pool in _MAILLARD_ALDEHYDES:
        assert float(oav_series(sweet.as_trajectory(), thresholds, pool)[-1]) > 0.0
    # Sotolon and the honey phenylacetaldehyde clear their perception thresholds (OAV > 1) at these
    # aged-Sauternes-range levels — the sweet-wine thermal signature is genuinely perceptible.
    assert float(oav_series(sweet.as_trajectory(), thresholds, "sotolon")[-1]) > 1.0
    assert float(oav_series(sweet.as_trajectory(), thresholds, "phenylacetaldehyde")[-1]) > 1.0


def test_the_reroute_no_longer_starves_maillard_now_that_precursors_are_speciated():
    """THE D-100 TRIPWIRE, FLIPPED (decision D-100; supersedes the D-99 known-limitation pin).

    The deleted test pinned a pathology: `FuselAminoAcidReroute` (D-33) re-sourced Ehrlich fusel
    carbon from the LUMPED `amino_acids` pool — i.e. from ARGININE, which makes no higher alcohol —
    and at D-99's honest ~3.8x fusel rise it drained that pool to ~0 and drove sotolon from OAV
    ~1.2 to silence. Two speciated-scale consumers were sharing one lumped substrate.

    D-100 speciates the pool, and this test asserts the inverse of what its predecessor asserted:
    with the reroute ON, sotolon is now clearly perceptible. Kept as a POSITIVE assertion rather
    than deleted outright so the coverage the old test carried does not vanish with the pathology.
    """
    thresholds = load_thresholds()

    def sweet(disable_reroute: bool):
        cs = compile_scenario(
            _wine(
                [_begin_aging(_FERMENT_DAYS)],
                amino_acids_gpl=0.8,
                brix=_SWEET_BRIX,
                duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
            )
        )
        if disable_reroute:
            cs.process_set.disable(FuselAminoAcidReroute.name)
        res = cs.run()
        assert res.success
        return res

    with_reroute = sweet(disable_reroute=False)
    without = sweet(disable_reroute=True)

    # THE FIX, at the pool the pathology ran through: arginine now ends at essentially the same
    # level whether fusels are re-routed or not. The reroute cannot touch it — arginine makes no
    # higher alcohol, so `FuselAminoAcidReroute.touches` does not name it. What little difference
    # remains is second-order (the reroute spares sugar, which shifts the ferment slightly).
    arg_on = float(with_reroute.series("amino_acids")[-1])
    arg_off = float(without.series("amino_acids")[-1])
    assert arg_on > 0.5 * arg_off  # was ~1e-5 vs ~0.2 under the lump: five orders down
    assert arg_on > 1.0e-3

    # THE CONSEQUENCE, inverted: sotolon is perceptible WITH the reroute running. The old test
    # asserted OAV < 0.01 here and > 1.0 without; both routes now clear the threshold.
    sotolon_on = float(oav_series(with_reroute.as_trajectory(), thresholds, "sotolon")[-1])
    assert sotolon_on > 1.0
    assert float(oav_series(without.as_trajectory(), thresholds, "sotolon")[-1]) > 1.0

    # WHAT D-100 DOES NOT CLAIM TO FIX — the honest boundary. Threonine feeds BOTH propanol
    # (Ehrlich) and sotolon (D-87), so the reroute still costs sotolon something. That competition
    # is REAL CHEMISTRY over one molecule, unlike the retired fusels-vs-arginine competition, so
    # the model SHOULD show it: sotolon stays below its no-reroute value, just no longer silenced.
    assert sotolon_on < float(oav_series(without.as_trajectory(), thresholds, "sotolon")[-1])


def test_thermal_aroma_from_drained_precursors_requires_autolysis():
    """The D-100 EMERGENT consequence: aging precursors are dominantly autolysis-sourced.

    Not a wiring choice — a consequence. The branched-chain and aromatic amino acids are Ehrlich
    substrates, so a ferment consumes them (real must carries ~30-60 mg/L leucine while wine makes
    ~150-250 mg/L isoamyl alcohol — most higher alcohol is synthesised de novo from sugar, and the
    catabolic part exhausts the pool). Autolysis (D-34) is the only amino-acid SOURCE in the model,
    and since D-100 it releases the full must spectrum rather than pure arginine — so it restores
    the precursors, and with them the aroma. That lees contact enriches Maillard/Strecker character
    is the published sur-lie mechanism, arriving as an emergent consequence of speciation rather
    than as a modelled rule. The lumped model could not express it: arginine stood in for
    precursors that were long gone, so it made thermal aldehydes from a pool whose real precursor
    content was zero.

    THE DIRECTION IS SOURCED; THE MAGNITUDE ASSERTED HERE IS NOT. The `< 1e-9` bound below is
    **what the current speculative re-route fraction gives**, NOT a validated prediction.

    D-103 CORRECTS THE NUMBERS D-100 PUT IN THIS DOCSTRING, AND THIS TEST CORRECTLY DOES NOT MOVE.
    D-100 wrote that the catabolic fraction was "~0.5 via the shared `K_amino_acids` gate" against
    "a literature catabolic contribution nearer 20-50%". Measured exactly, it is **0.192** at this
    test's own dose and 0.21-0.33 at a realistic must; the 20-50% was **uncited**, and the sourced
    contribution is LOWER still (Rollero 2017: isoamyl 2-8%, isobutanol 5-15% by U-13C labelling)
    — so that band would have *acquitted* this model rather than convicting it.

    THE CITED CORROBORATION COULD NOT TEST THE CLAIM. The phenylacetaldehyde figure was taken from
    an ON-LEES run — and on lees autolysis refills phenylalanine to ~54 mg/L, so the pool ends
    FULL and nothing throttled. It measured the one run where the alleged over-drain is absent.
    (The direction is backwards too: over-draining would make that aldehyde LOW, and the evidence
    offered was that it is HIGH.)

    D-100 prophesied "if the re-route fraction is ever bounded, this test should move". **It is
    not bounded and this test has not moved** — D-103 turns no knob, because there is no single
    fraction to bound (the real defect is the gate's per-species SHAPE: ~8% for isoamyl vs ~82-93%
    for propanol, where reality is uniformly low) and `K_amino_acids` is one shared scalar that
    cannot reshape it. Do NOT read this test's unchanged bound as contradicting D-103 — an
    unchanged test IS the receipt for a documentation-only correction. What would move it is the
    missing anabolic/protein sink, or a re-formed gate; both are their own beats.
    """

    def sweet(autolysis: float | None):
        kw = {}
        if autolysis is not None:
            kw["autolysis_rate_per_h"] = autolysis
        cs = compile_scenario(
            _wine(
                [_begin_aging(_FERMENT_DAYS)],
                amino_acids_gpl=0.8,
                brix=_SWEET_BRIX,
                duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
                **kw,
            )
        )
        res = cs.run()
        assert res.success
        return res

    no_lees = sweet(None)
    on_lees = sweet(1.0e-3)

    # Fermentation strips the Ehrlich precursors to ~nothing without lees...
    for precursor in ("leucine", "isoleucine", "phenylalanine"):
        assert float(no_lees.series(precursor)[-1]) < 1.0e-4
        # ...and autolysis restores them by orders of magnitude.
        assert float(on_lees.series(precursor)[-1]) > 100.0 * max(
            float(no_lees.series(precursor)[-1]), 1e-12
        )

    # METHIONINE IS THE CONTROL: no fusel consumes it, so it survives the ferment on its own and
    # does NOT depend on lees for its presence. That the model distinguishes these two cases —
    # rather than treating "amino acids" as one thing — is the whole point of D-100.
    assert float(no_lees.series("methionine")[-1]) > 1.0e-3

    # The aroma follows the precursor: the leucine-derived malty aldehyde is silent without lees
    # and real with them.
    assert float(no_lees.series("3_methylbutanal")[-1]) < 1.0e-9
    assert float(on_lees.series("3_methylbutanal")[-1]) > 1.0e-6


def test_maillard_closes_carbon_and_nitrogen_end_to_end():
    # MaillardStrecker draws carbon from amino_acids into the six thermal products + CO₂, deaminates
    # the nitrogen to N — so BOTH ledgers must close end to end through the full ferment + sealed
    # sweet aging. No external flow at all (no O₂ dose; the amino-acid + sugar are t0 initials), so
    # total_carbon and total_nitrogen are both flat (final == initial).
    #
    # DELIBERATELY KEEPS THE D-33 REROUTE ON (unlike the magnitude tests above, D-99): this asserts
    # only CLOSURE and nonnegativity — invariants that must hold with EVERY Process active, and the
    # reroute (another amino_acids consumer) makes the ledger stress stronger, not weaker. It makes
    # no aroma-magnitude claim, so the reroute draining the pool cannot render it vacuous.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)],
            amino_acids_gpl=0.8,
            brix=_SWEET_BRIX,
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    assert_conserved(
        traj.as_trajectory(),
        total_carbon(cs.schema, biomass_carbon_fraction=f_c),
        label="carbon",
    )
    assert_conserved(
        traj.as_trajectory(),
        total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n),
        label="nitrogen",
    )
    # atol at 10x the solver's own atol (1e-9, runtime.integrate): N reaches a TRUE zero on this
    # scenario, and a state legitimately at zero undershoots by ~the integrator's absolute
    # tolerance — asserting tighter than the solver's own promise tests scipy, not the model. That
    # this only began mattering at D-100 is the fix working: the reroute's arginine draw
    # over-released nitrogen ~4x (D-33's documented lump), propping up the YAN that growth ate;
    # each precursor now releases the nitrogen it actually carries, so N reaches zero as it should.
    # Verified as noise, not drift: the excursion tracks the solver's atol ~1:1 (1e-9 -> -1.1e-9,
    # 1e-11 -> -1.3e-12, 1e-13 -> -1.1e-14) with the trajectory itself unchanged.
    assert_nonnegative(
        traj.as_trajectory(),
        (
            *(spec.pool for spec in AMINO_ACID_SPECS),  # every speciated pool (D-100)
            *_MAILLARD_ALDEHYDES,
            "N",
        ),
        atol=1e-8,
    )


# -- Caramelization (decision D-88) — the NON-oxidative THERMAL browning axis, end to end -------
#
# The O₂-INDEPENDENT thermal mirror of PhenolicBrowning (D-74): a SEALED SWEET wine (or, D-90,
# high-residual beer) browns thermally (residual sugar → melanoidin, raising shared A420) with no
# oxygen — where a DRY beverage (S ≈ 0 at the aging segment) is byte-for-byte inert. The FIRST aging
# Process to consume core S, so it carries the sugar carbon into the on-ledger melanoidin park
# (total_carbon closes). These pin: the compile-seam gate (medium-agnostic D-90, disabled→begin_
# aging); the sealed-sweet browning vs a dry control through the FULL pipeline; end-to-end carbon
# closure with core S consumed.


def test_caramelization_gated_by_begin_aging_medium_agnostic():
    # Caramelization is MEDIUM-AGNOSTIC (D-90: both media carry the melanoidin park + wire the
    # Process — beer's residual dextrins caramelize too) — present in BOTH sets, unlike wine-only
    # MaillardStrecker/MaillardBrowning — and rides the aging gate: disabled at compile, on by
    # begin_aging.
    assert Caramelization.name in get_medium("wine").build_process_set()
    assert Caramelization.name in get_medium("beer").build_process_set()
    cs = compile_scenario(_wine([_begin_aging(_FERMENT_DAYS)], brix=_SWEET_BRIX))
    assert Caramelization.name in cs.process_set
    assert not cs.process_set.is_enabled(Caramelization.name)  # off at compile
    event = next(e for e in cs.events if e.label.startswith("begin_aging"))
    assert event.reconfigure is not None
    event.reconfigure(cs.process_set)
    assert cs.process_set.is_enabled(Caramelization.name)  # begin_aging turns it on


def test_caramelization_browns_sealed_sweet_wine_not_dry():
    # THE discriminating end-to-end: a SEALED (no add_oxygen) SWEET wine browns thermally —
    # melanoidin
    # + the A420 browning index climb from 0 with no oxygen — where a DRY wine (S ≈ 0 at the aging
    # segment) leaves both at ~0 (byte-for-byte inert, the soft sugar gate). No O₂ anywhere, so this
    # is browning the O₂-driven PhenolicBrowning (D-74) could never produce.
    sweet = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)],  # SEALED — no add_oxygen
            brix=_SWEET_BRIX,
            aging_celsius=30.0,
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    ).run()
    dry = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)],
            brix=24.0,  # dry-fermenting control → residual sugar ≈ 0 ⇒ caramelization ~inert
            aging_celsius=30.0,
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    ).run()
    assert sweet.success and dry.success
    # The sweet wine browns: melanoidin forms and A420 rises, all from residual sugar + heat, no O₂.
    assert float(sweet.series("melanoidin")[-1]) > 0.0
    assert float(sweet.series("A420")[-1]) > 0.01  # a visible browning index
    # The dry wine (no residual sugar) makes essentially none.
    assert float(dry.series("melanoidin")[-1]) < 1e-6
    assert float(sweet.series("A420")[-1]) > 100.0 * float(dry.series("A420")[-1])
    # Residual sugar declines (consumed into melanoidin) but is not exhausted (browning is slow).
    assert float(sweet.series("S")[-1]) > 50.0


def test_caramelization_closes_carbon_end_to_end():
    # Caramelization CONSUMES core S into the melanoidin carbon-park — the first aging Process to
    # touch S — so total_carbon must still close end to end (the sugar carbon parks in the on-ledger
    # melanoidin pool, not destroyed). No external flow (no O₂ dose), so the ledger is flat.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS)],
            brix=_SWEET_BRIX,
            aging_celsius=30.0,
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.parameters.value("biomass_C_fraction")
    assert_conserved(
        traj.as_trajectory(),
        total_carbon(cs.schema, biomass_carbon_fraction=f_c),
        label="carbon",
    )
    assert_nonnegative(traj.as_trajectory(), ("S", "melanoidin", "A420"), atol=1e-9)


def test_thermal_and_oxidative_axes_coexist_and_close_end_to_end():
    # THE five-way interaction (the shared-pool stress test): a SWEET + OXYGENATED +
    # amino-acid-dosed
    # aged wine runs ALL FIVE amino-acid/browning aging routes at once — BOTH Strecker routes
    # (oxidative D-75 + thermal D-87) drawing the shared amino_acids, all THREE browning routes
    # writing the shared A420 (oxidative D-74 + sugar-only thermal Caramelization D-88 +
    # amino-acid-incorporating thermal MaillardBrowning D-89), with Caramelization AND
    # MaillardBrowning
    # BOTH consuming core S and MaillardBrowning ALSO drawing the shared amino_acids into the
    # N-bearing
    # maillard_melanoidin park. This is the stress case for the nitrogen ledger: MaillardBrowning is
    # the FIRST aging Process to PARK nitrogen in a product pool (D-87 deaminates its N back to N;
    # D-89 retains it), so with it live, three sinks pull amino_acids at once and total_nitrogen
    # must
    # still close through the maillard_melanoidin weighting. Each route sizes its OWN draw and
    # ProcessSet SUMS them, so by additivity carbon + nitrogen must still close and no shared pool
    # goes negative — the combination where a shared-pool interaction bug could hide, and exactly
    # where conservation would catch it. (The aa gate throttles all amino_acids draws to 0 as the
    # pool empties, so it never goes negative.)
    # The D-33 reroute is ISOLATED (D-99): the five AGING routes stressed here each need the
    # shared amino_acids pool, and the reroute would empty it during fermentation, so the
    # "all five fired" assertions below (sotolon > 0 etc.) would pass on ~noise. The reroute is
    # a SIXTH, fermentation-phase process — not one of the five aging routes — so isolating it
    # neither removes a stressed route nor weakens the aging-phase carbon/nitrogen closure (it
    # is inactive during aging regardless). See _run_thermal_isolating_reroute + D-100.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oxygen(_FERMENT_DAYS, 60.0)],  # OXYGENATED
            amino_acids_gpl=0.8,  # both Strecker routes have substrate
            brix=_SWEET_BRIX,  # SWEET — both browning routes + Caramelization's S driver
            aging_celsius=30.0,
            duration_days=_FERMENT_DAYS + _SWEET_AGING_DAYS,
        )
    )
    cs.process_set.disable(FuselAminoAcidReroute.name)
    traj = cs.run()
    assert traj.success
    tj = traj.as_trajectory()
    f_c = cs.parameters.value("biomass_C_fraction")
    f_n = cs.parameters.value("biomass_N_fraction")
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    n_of = total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n)
    # The O₂ dose is the only external flow and carries neither carbon nor nitrogen (o2 off every
    # ledger), so BOTH ledgers stay flat (final == initial) with all four axes live.
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert all(n_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert_conserved(tj, c_of, label="carbon")
    assert_conserved(tj, n_of, label="nitrogen")
    # No shared pool goes negative — each relative-depletion gate keeps its own amino acid ≥ 0
    # (D-100), and both carbon-parks (melanoidin + the N-bearing maillard_melanoidin) stay
    # nonnegative.
    #
    # atol at 10x the solver's own atol (1e-9, runtime.integrate): N reaches a TRUE zero here, and
    # a state legitimately at zero undershoots by ~the integrator's absolute tolerance — asserting
    # tighter than the solver's own promise tests scipy, not the model. That this only began
    # mattering at D-100 is the fix working: the reroute's arginine draw over-released nitrogen ~4x
    # (D-33's documented lump), propping up the YAN; each precursor now releases its REAL nitrogen.
    # Verified as noise, not drift: the excursion tracks the solver's atol ~1:1 (1e-9 -> -1.1e-9,
    # 1e-11 -> -1.3e-12, 1e-13 -> -1.1e-14) with the trajectory unchanged.
    assert_nonnegative(
        tj,
        (
            "amino_acids",
            "amino_acids_generic",
            *(spec.pool for spec in AMINO_ACID_SPECS),
            "S",
            "melanoidin",
            "maillard_melanoidin",
            "A420",
            "o2",
            *_MAILLARD_ALDEHYDES,
            "N",
        ),
        atol=1e-8,
    )
    # All five routes actually fired: the two thermal browning polymers (sugar-only melanoidin +
    # N-bearing maillard_melanoidin) + thermal aldehydes (sotolon, a thermal-route-only marker) are
    # live ALONGSIDE the oxidative browning (A420 gets all three browning routes) and the oxidative
    # Strecker (methional/phenylacetaldehyde, shared with the thermal route).
    assert float(traj.series("melanoidin")[-1]) > 0.0  # sugar-only thermal browning (D-88)
    assert float(traj.series("maillard_melanoidin")[-1]) > 0.0  # N-bearing thermal browning (D-89)
    assert float(traj.series("sotolon")[-1]) > 0.0  # thermal Strecker (O₂-independent marker)
    assert float(traj.series("A420")[-1]) > 0.0  # browning index (all three routes)
    assert float(traj.series("methional")[-1]) > 0.0  # Strecker (both routes)


# -- OakExtraction (decision D-77) — the NON-oxidative barrel/chip aroma axis, end to end -------

_OAK_EXTRACTIVES = ("whiskey_lactone", "vanillin", "guaiacol", "eugenol", "furaneol")
_OAK_CEILINGS = tuple(f"{c}_ceiling" for c in _OAK_EXTRACTIVES)
#: The bourbon aroma-soak-back subset (D-93/D-94) — mirrors compile._OAK_SPIRIT_AROMAS. An ex-
#: bourbon barrel BUMPS these four (vanilla/coconut/char — D-93 — + caramel furaneol — D-94);
#: eugenol (clove) + ellagitannin are untouched.
_OAK_SPIRIT_AROMAS = ("vanillin", "whiskey_lactone", "guaiacol", "furaneol")


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
    # Furaneol (caramel/toffee) RISES with toast (a thermal sugar-degradation furanone, co-varying
    # with guaiacol/eugenol — heavy/charred oak gives the most caramel; the D-94 sourced ordering).
    assert heavy["furaneol"] > medium["furaneol"] > light["furaneol"]
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
    # (bourbon-barrel stouts / oak-aged sours), setting the same 6 ceilings from oak_gpl × toast
    # yield. (Was a "wine-only" rejection before D-86 — a legitimate expectation flip, not a
    # weakened test.) The unknown-toast rejection above still holds for beer too.
    cs = compile_scenario(_beer([_add_oak(14.0, 4.0, "medium")]))
    event = next(e for e in cs.events if e.label.startswith("add_oak"))
    assert event.mutate is not None  # dose only (a pure ceiling mutate, like the wine add_oak)
    after = event.mutate(cs.schema, cs.y0.copy())
    # The medium-toast dose sets each ceiling to oak_gpl × oak_yield_<compound>_medium (> 0).
    for compound in (*_OAK_EXTRACTIVES, "ellagitannin"):
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
        traj.as_trajectory(),
        total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n),
        label="nitrogen",
    )


def test_barrel_beer_oak_raises_the_oak_oavs_and_astringency():
    # The sensory payoff on beer (bourbon-barrel stout / oak-aged sour): oak aging lifts the four
    # oak-extractive OAVs from 0 (read through the D-67 lens against the beer-matrix thresholds,
    # D-86) and the ellagitannin astringency readout goes positive. An un-oaked beer reads 0.
    thresholds = load_thresholds()
    oaked = compile_scenario(_beer([_begin_aging(14.0), _add_oak(14.0, 6.0, "medium")])).run()
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
    oaked_cs = compile_scenario(
        _beer([_begin_aging(14.0), _add_oak(14.0, 6.0, "light"), _add_oxygen(14.0, o2_dose)])
    )  # light toast ⇒ most ellagitannin (strongest protection)
    oaked = oaked_cs.run()
    unoaked = compile_scenario(_beer([_begin_aging(14.0), _add_oxygen(14.0, o2_dose)])).run()
    assert oaked.success and unoaked.success
    assert float(oaked.series("A420")[-1]) < float(unoaked.series("A420")[-1])
    assert float(oaked.series("acetaldehyde")[-1]) < float(unoaked.series("acetaldehyde")[-1])
    assert float(oaked.series("A420")[-1]) > 0.0  # PARTIAL, not total — still some browning
    # The oak axis + EllagitanninOxidation are all off every ledger (wood-derived), so even a
    # fully-active oaked+oxygenated beer closes carbon + nitrogen (mirrors the wine oaked-ledger
    # test; the un-oaked beer C+N check adds the off-ledger-slots-are-inert half).
    f_c = oaked_cs.parameters.value("biomass_C_fraction")
    f_n = oaked_cs.parameters.value("biomass_N_fraction")
    assert_conserved(
        oaked.as_trajectory(),
        total_carbon(oaked_cs.schema, biomass_carbon_fraction=f_c),
        label="carbon",
    )
    assert_conserved(
        oaked.as_trajectory(),
        total_nitrogen(oaked_cs.schema, biomass_nitrogen_fraction=f_n),
        label="nitrogen",
    )


# -- D-91: barrel fill-number depletion (a reused barrel extracts LESS — the add_oak dose lever) --
# fill_number is an ACROSS-FILL dose input, not new physics: it scales every add_oak ceiling by
# oak_fill_retention ** (fill_number - 1) at charge time. fill_number = 1 (a fresh first-fill
# barrel, the default) is UNSCALED — byte-for-byte the pre-D-91 dose — so the whole existing oak
# suite pins the first-fill case; these tests pin the DEPLETION. The signature lever is barrel-aged
# BEER, so a beer scenario carries the end-to-end weight.


def test_fill_number_defaults_to_first_fill_byte_for_byte():
    # The backward-compat anchor: omitting fill_number is identical to fill_number = 1 (r ** 0 ==
    # 1.0 exactly), so the ceilings match to the last bit — the guarantee that every pre-D-91 wine +
    # beer trajectory is untouched. Explicit and implicit first-fill produce identical ceilings.
    implicit = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium")])
    )
    explicit = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium", fill_number=1)])
    )
    ev_i = next(e for e in implicit.events if e.label.startswith("add_oak"))
    ev_e = next(e for e in explicit.events if e.label.startswith("add_oak"))
    assert ev_i.mutate is not None and ev_e.mutate is not None
    after_i = ev_i.mutate(implicit.schema, implicit.y0.copy())
    after_e = ev_e.mutate(explicit.schema, explicit.y0.copy())
    for compound in (*_OAK_EXTRACTIVES, "ellagitannin"):
        name = f"{compound}_ceiling"
        # byte-for-byte: explicit first-fill == implicit == the raw oak_gpl × yield (unscaled)
        assert implicit.schema.get(after_i, name) == explicit.schema.get(after_e, name)
        assert (
            implicit.schema.get(after_i, name)
            == 4.0 * implicit.param_values[f"oak_yield_{compound}_medium"]
        )


def test_higher_fill_number_geometrically_discounts_the_ceilings():
    # The D-91 core: each prior fill multiplies every ceiling by oak_fill_retention. A fresh
    # (fill 1), second-fill (fill 2) and fourth-fill (fill 4) barrel at the SAME oak_gpl/toast set
    # ceilings in the ratio 1 : r : r³ — a decreasing geometric sequence, all 5 extractives.
    def ceilings(fill: int) -> dict[str, float]:
        cs = compile_scenario(
            _wine(
                [
                    _begin_aging(_FERMENT_DAYS),
                    _add_oak(_FERMENT_DAYS, 4.0, "medium", fill_number=fill),
                ]
            )
        )
        ev = next(e for e in cs.events if e.label.startswith("add_oak"))
        assert ev.mutate is not None
        after = ev.mutate(cs.schema, cs.y0.copy())
        return {
            c: float(cs.schema.get(after, f"{c}_ceiling"))
            for c in (*_OAK_EXTRACTIVES, "ellagitannin")
        }

    r = compile_scenario(_wine([])).param_values["oak_fill_retention"]
    fresh, second, fourth = ceilings(1), ceilings(2), ceilings(4)
    for compound in (*_OAK_EXTRACTIVES, "ellagitannin"):
        assert fresh[compound] > second[compound] > fourth[compound] > 0.0  # monotone depletion
        assert second[compound] == pytest.approx(fresh[compound] * r)  # one prior fill ⇒ × r
        assert fourth[compound] == pytest.approx(fresh[compound] * r**3)  # three prior fills ⇒ × r³


def test_reused_barrel_beer_reads_lower_oak_oavs_and_astringency_end_to_end():
    # The motivating BEER payoff (D-91): a first-fill bourbon-barrel stout vs the SAME beer in a
    # fourth-fill (near-neutral) barrel — identical ferment + aging + oak_gpl/toast, differing ONLY
    # in fill_number. The reused barrel reads LOWER on every oak OAV and lower ellagitannin
    # astringency: a depleted barrel imparts less wood character (barrel-aged-beer practice).
    thresholds = load_thresholds()
    fresh = compile_scenario(
        _beer([_begin_aging(14.0), _add_oak(14.0, 6.0, "medium", fill_number=1)])
    ).run()
    reused = compile_scenario(
        _beer([_begin_aging(14.0), _add_oak(14.0, 6.0, "medium", fill_number=4)])
    ).run()
    assert fresh.success and reused.success
    for compound in _OAK_EXTRACTIVES:
        oav_fresh = float(oav_series(fresh.as_trajectory(), thresholds, compound)[-1])
        oav_reused = float(oav_series(reused.as_trajectory(), thresholds, compound)[-1])
        assert oav_fresh > oav_reused > 0.0  # both oaked, but the reused barrel imparts less
    # Astringency (oak ellagitannin alone on beer, D-86): reused barrel lower, still positive.
    astr_fresh = float(astringency_series(fresh.as_trajectory())[-1])
    astr_reused = float(astringency_series(reused.as_trajectory())[-1])
    assert astr_fresh > astr_reused > 0.0


def test_add_oak_rejects_a_zeroth_or_fractional_fill_number():
    # fill_number counts a physical barrel use — first/second/third — so it must be an integer ≥ 1.
    # A "zeroth fill" (< 1) and a fractional fill are both meaningless and rejected loudly at
    # compile (the toast-string rejection pattern), not silently coerced.
    for bad in (0, -1, 2.5):
        with pytest.raises(ValueError, match="fill_number must be an integer"):
            compile_scenario(
                _wine(
                    [
                        _begin_aging(_FERMENT_DAYS),
                        _add_oak(_FERMENT_DAYS, 4.0, "medium", fill_number=bad),
                    ]
                )
            )


# -- D-92: bourbon-barrel spirit soak-back — an ex-spirit barrel donates ETHANOL (raises ABV) ------
#
# A SEPARATE contribution from the wood extractives (D-77/78) and fill-number depletion (D-91): the
# ethanol comes from residual SPIRIT soaked into the staves, not the wood. add_oak {spirit: bourbon}
# adds a DISCRETE ethanol dose to the core E slot (the add_oxygen precedent), decoupled from oak_gpl
# and depleting with fill_number via its OWN steep spirit_soak_retention. spirit ABSENT ⇒ no ethanol
# ⇒ byte-for-byte the pre-D-92 charge (the existing oak suite pins no-spirit). ETHANOL (ABV) only
# this beat — the bourbon AROMA congeners are deferred.


def test_spirit_soak_back_absent_leaves_ethanol_untouched_byte_for_byte():
    # The backward-compat anchor: an add_oak with NO spirit does not touch E, and its ceilings are
    # exactly the raw oak_gpl × yield — so every pre-D-92 wine + beer trajectory is untouched. A
    # bourbon charge at the SAME dose sets IDENTICAL ceilings (soak-back is orthogonal to the wood).
    plain = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium")])
    )
    bourbon = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium", spirit="bourbon")]
        )
    )
    ev_p = next(e for e in plain.events if e.label.startswith("add_oak"))
    ev_b = next(e for e in bourbon.events if e.label.startswith("add_oak"))
    assert ev_p.mutate is not None and ev_b.mutate is not None
    after_p = ev_p.mutate(plain.schema, plain.y0.copy())
    after_b = ev_b.mutate(bourbon.schema, bourbon.y0.copy())
    # No spirit ⇒ E untouched by the dose (mutate only writes off-ledger ceilings).
    assert plain.schema.get(after_p, "E") == plain.schema.get(plain.y0, "E")
    # spirit does NOT touch the NON-soak ceilings (eugenol/clove + the ellagitannin taste tannin are
    # not bourbon congeners, D-93) — identical to the plain charge at the same oak_gpl/toast. The
    # three aroma-soak ceilings (vanillin/whiskey_lactone/guaiacol) ARE bumped — see the D-93 tests.
    untouched = set(_OAK_EXTRACTIVES + ("ellagitannin",)) - set(_OAK_SPIRIT_AROMAS)
    for compound in untouched:
        name = f"{compound}_ceiling"
        assert bourbon.schema.get(after_b, name) == plain.schema.get(after_p, name)


def test_bourbon_spirit_adds_the_full_first_fill_ethanol_bolus():
    # A first-fill (default fill_number = 1) ex-bourbon barrel donates the FULL soak-back ethanol
    # (spirit_soak_retention ** 0 == 1.0 exactly): the E slot rises by exactly spirit_soak_ethanol_
    # bourbon g/L (~8 g/L ≈ 1% ABV), the signature barrel-aged-stout ABV gain.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium", spirit="bourbon")]
        )
    )
    ev = next(e for e in cs.events if e.label.startswith("add_oak"))
    assert ev.mutate is not None
    after = ev.mutate(cs.schema, cs.y0.copy())
    delta_e = float(cs.schema.get(after, "E") - cs.schema.get(cs.y0, "E"))
    assert delta_e == pytest.approx(cs.param_values["spirit_soak_ethanol_bourbon"], rel=1e-12)
    assert delta_e > 0.0


def test_spirit_soak_back_depletes_geometrically_and_steeper_than_the_wood():
    # A reused ex-bourbon barrel donates LESS residual spirit: the E bolus at fills 1/2/3 is in the
    # ratio 1 : r_s : r_s² (r_s = spirit_soak_retention). And the spirit depletes STEEPER than the
    # wood extractables (r_s < oak_fill_retention) — "first-fill bourbon barrel" is the term of art
    # because a refill barrel is largely rinsed of spirit by its first use.
    def bolus(fill: int) -> float:
        cs = compile_scenario(
            _wine(
                [
                    _begin_aging(_FERMENT_DAYS),
                    _add_oak(_FERMENT_DAYS, 4.0, "medium", fill_number=fill, spirit="bourbon"),
                ]
            )
        )
        ev = next(e for e in cs.events if e.label.startswith("add_oak"))
        assert ev.mutate is not None
        after = ev.mutate(cs.schema, cs.y0.copy())
        return float(cs.schema.get(after, "E") - cs.schema.get(cs.y0, "E"))

    params = compile_scenario(_wine([])).param_values
    r_s, r_wood = params["spirit_soak_retention"], params["oak_fill_retention"]
    first, second, third = bolus(1), bolus(2), bolus(3)
    assert first > second > third > 0.0  # monotone depletion of the residual spirit
    assert second == pytest.approx(first * r_s)  # one prior fill ⇒ × r_s
    assert third == pytest.approx(first * r_s**2)  # two prior fills ⇒ × r_s²
    assert r_s < r_wood  # spirit rinses out faster than the wood extractables deplete


def test_spirit_soak_back_conserves_carbon_across_the_jump():
    # The crown-jewel ledger (D-92): ethanol is ON the carbon+mass ledger, so the soak-back dose
    # INJECTS carbon — but the scheduler books it as a POSITIVE external flow (the add_sugar
    # precedent), so the run-wide identity final == initial + Σ flows still closes to machine
    # precision. (Contrast the off-ledger o2/ceiling doses, whose flows are carbon-free.)
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium", spirit="bourbon")]
        )
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.param_values["biomass_C_fraction"]
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    injected = sum(c_of(flow.delta) for flow in traj.external_flows)
    # The soak-back flow carries POSITIVE carbon (ethanol's), unlike carbon-free o2/ceiling doses.
    assert injected > 0.0
    c_initial, c_final = c_of(cs.y0), c_of(traj.y[:, -1])
    assert c_final == pytest.approx(c_initial + injected, abs=1e-9)


def test_add_oak_rejects_an_unknown_spirit():
    # spirit is a categorical asserting an ex-spirit barrel; an unknown one is rejected loudly at
    # compile (the toast-string rejection pattern), never silently ignored.
    with pytest.raises(ValueError, match="unknown spirit"):
        compile_scenario(
            _wine(
                [
                    _begin_aging(_FERMENT_DAYS),
                    _add_oak(_FERMENT_DAYS, 4.0, "medium", spirit="tequila"),
                ]
            )
        )


def test_bourbon_barrel_stout_gains_abv_end_to_end():
    # The motivating BEER payoff (D-92): a bourbon-barrel imperial stout finishes at HIGHER ABV than
    # the SAME beer aged in an identical but spirit-free oak barrel — the residual-spirit soak-back.
    # A fourth-fill (near-rinsed) bourbon barrel donates LESS than a first-fill (steep depletion).
    def final_ethanol(iv: Intervention) -> float:
        traj = compile_scenario(_beer([_begin_aging(14.0), iv])).run()
        assert traj.success
        return float(traj.series("E")[-1])

    no_spirit = final_ethanol(_add_oak(14.0, 6.0, "medium"))
    first_fill = final_ethanol(_add_oak(14.0, 6.0, "medium", spirit="bourbon"))
    fourth_fill = final_ethanol(_add_oak(14.0, 6.0, "medium", fill_number=4, spirit="bourbon"))
    # First-fill bourbon barrel adds the most ABV; a near-neutral reused barrel adds little; a
    # spirit-free oak barrel adds none — the ordering a barrel-aged-beer program manages.
    assert first_fill > fourth_fill > no_spirit


# -- D-93: bourbon-barrel AROMA soak-back — an ex-spirit barrel donates residual CONGENERS -------
#
# The second half of the D-92 soak-back. Bourbon matures in CHARRED NEW OAK, so its residual spirit
# reads vanilla/coconut/char-forward: add_oak {spirit: bourbon} BUMPS the vanillin/whiskey_lactone/
# guaiacol ceilings (a DELIBERATE subset; not eugenol/clove or the ellagitannin tannin) by
# spirit_soak_<c>_bourbon × spirit_scale, and OakExtraction (D-77) leaches them in GRADUALLY on top
# of the wood diffusion. A CEILING bump — the ONLY wood + spirit ADDITIVE form (a bolus into the
# pool is erased by the extraction gate). Off the carbon/mass ledger (aroma ceilings, iso_alpha).
# spirit ABSENT ⇒ no bump ⇒ byte-for-byte the pre-D-92 charge on the aroma ceilings too.


def test_bourbon_aroma_bumps_the_signature_ceilings_by_exactly_the_spirit_soak():
    # A first-fill ex-bourbon barrel raises EACH of the three signature aroma ceilings
    # (vanillin/whiskey_lactone/guaiacol) by exactly its spirit_soak_<c>_bourbon g/L (retention**0
    # == 1.0), ON TOP of the wood ceiling — while eugenol (clove) and ellagitannin stay wood-only.
    plain = compile_scenario(
        _wine([_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium")])
    )
    bourbon = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium", spirit="bourbon")]
        )
    )
    ev_p = next(e for e in plain.events if e.label.startswith("add_oak"))
    ev_b = next(e for e in bourbon.events if e.label.startswith("add_oak"))
    assert ev_p.mutate is not None and ev_b.mutate is not None
    after_p = ev_p.mutate(plain.schema, plain.y0.copy())
    after_b = ev_b.mutate(bourbon.schema, bourbon.y0.copy())
    for compound in _OAK_SPIRIT_AROMAS:
        name = f"{compound}_ceiling"
        bump = bourbon.param_values[f"spirit_soak_{compound}_bourbon"]
        assert bump > 0.0
        # The bourbon ceiling is the wood ceiling PLUS the spirit bump — additive, not max/replace.
        assert bourbon.schema.get(after_b, name) == pytest.approx(
            plain.schema.get(after_p, name) + bump, rel=1e-12
        )
    # Non-signature extractives (clove + tannin) are NOT bumped — the deliberate subset.
    for compound in ("eugenol", "ellagitannin"):
        name = f"{compound}_ceiling"
        assert bourbon.schema.get(after_b, name) == plain.schema.get(after_p, name)


def test_bourbon_aroma_bump_depletes_geometrically_with_fill_number():
    # The residual-spirit congeners deplete with reuse by the SAME spirit_soak_retention as the
    # ethanol (one residual spirit, one depletion): the vanillin ceiling BUMP at fills 1/2/3 is in
    # the ratio 1 : r_s : r_s² — steeper than the wood's oak_fill_retention (discounts wood).
    def vanillin_bump(fill: int) -> float:
        # Isolate the spirit contribution: bourbon ceiling minus the same-fill plain (wood) ceiling.
        def ceil(spirit: str | None) -> float:
            cs = compile_scenario(
                _wine(
                    [
                        _begin_aging(_FERMENT_DAYS),
                        _add_oak(_FERMENT_DAYS, 4.0, "medium", fill_number=fill, spirit=spirit),
                    ]
                )
            )
            ev = next(e for e in cs.events if e.label.startswith("add_oak"))
            assert ev.mutate is not None
            after = ev.mutate(cs.schema, cs.y0.copy())
            return float(cs.schema.get(after, "vanillin_ceiling"))

        return ceil("bourbon") - ceil(None)

    r_s = compile_scenario(_wine([])).param_values["spirit_soak_retention"]
    first, second, third = vanillin_bump(1), vanillin_bump(2), vanillin_bump(3)
    assert first > second > third > 0.0
    assert second == pytest.approx(first * r_s)
    assert third == pytest.approx(first * r_s**2)


def test_bourbon_aroma_leaches_in_gradually_and_reads_forward_end_to_end():
    # The end-to-end payoff AND the additive proof at runtime: a bourbon-barrel beer finishes with
    # HIGHER extracted vanillin/whiskey_lactone/guaiacol than the SAME beer in an identical
    # spirit-free oak barrel — OakExtraction actually reaches the raised (wood + spirit) ceiling, so
    # the spirit congeners are NOT erased by the gate (the whole reason for the ceiling-bump
    # design). Their OAVs rise too (vanilla/coconut/char read FORWARD); eugenol (not bumped) is not.
    thresholds = load_thresholds()

    def finals(iv: Intervention) -> tuple[dict[str, float], dict[str, float]]:
        traj = compile_scenario(_beer([_begin_aging(14.0), iv])).run()
        assert traj.success
        pools = {c: float(traj.series(c)[-1]) for c in (*_OAK_SPIRIT_AROMAS, "eugenol")}
        # Exercise the OAV readout (Beat 1a) too, not just the raw pools — the test claims the aroma
        # "reads forward", so it must go through the sensory pathway, not merely the concentrations.
        tj = traj.as_trajectory()
        oav = {c: float(oav_series(tj, thresholds, c)[-1]) for c in _OAK_SPIRIT_AROMAS}
        return pools, oav

    (no_spirit, no_oav) = finals(_add_oak(14.0, 6.0, "medium"))
    (bourbon, bourbon_oav) = finals(_add_oak(14.0, 6.0, "medium", spirit="bourbon"))
    for compound in _OAK_SPIRIT_AROMAS:
        assert bourbon[compound] > no_spirit[compound]  # spirit congeners actually leach in
        # The sensory readout lifts too: each signature aroma reads MORE FORWARD (higher OAV) with
        # the bourbon soak-back, and clears its perception threshold (OAV > 1) on its own.
        assert bourbon_oav[compound] > no_oav[compound]
        assert bourbon_oav[compound] > 1.0
    # eugenol is not a bourbon congener ⇒ its extracted pool is identical either way.
    assert bourbon["eugenol"] == pytest.approx(no_spirit["eugenol"], rel=1e-9)


def test_bourbon_aroma_soak_back_is_off_ledger():
    # The aroma bump moves NO conserved quantity (aroma ceilings + extracted pools are off the
    # carbon/mass ledger, the iso_alpha precedent) — so a bourbon run's ONLY carbon injection is the
    # D-92 ethanol dose. Total injected carbon equals the ethanol-only amount; carbon still closes.
    cs = compile_scenario(
        _wine(
            [_begin_aging(_FERMENT_DAYS), _add_oak(_FERMENT_DAYS, 4.0, "medium", spirit="bourbon")]
        )
    )
    traj = cs.run()
    assert traj.success
    f_c = cs.param_values["biomass_C_fraction"]
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    injected = sum(c_of(flow.delta) for flow in traj.external_flows)
    # Ethanol carbon: the soak-back ethanol bolus is the ONLY carbon the whole add_oak dose injects
    # (the aroma bumps are carbon-free). Its magnitude is spirit_soak_ethanol_bourbon g/L ethanol.
    ethanol_gpl = cs.param_values["spirit_soak_ethanol_bourbon"]
    e_only = cs.schema.zeros()
    e_only[cs.schema.slice("E")] = ethanol_gpl
    ethanol_carbon = c_of(e_only)
    assert injected == pytest.approx(ethanol_carbon, rel=1e-9)
    c_initial, c_final = c_of(cs.y0), c_of(traj.y[:, -1])
    assert c_final == pytest.approx(c_initial + injected, abs=1e-9)


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
