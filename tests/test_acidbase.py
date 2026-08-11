"""The wine pH charge-balance solver and its derived pH/TA readout (decision D-18).

Ranked headline-first. The keystone's proof-of-purpose is
``test_headline_malic_to_lactic_raises_ph``: the malic→lactic deacidification (the
chemistry MLF performs) raises pH by 0.1–0.3, demonstrated *without* an MLF Process
built — the solver responds to acid dynamics on its own. The rest pin the balance
(residual ≈ 0, monotonicity, smoothness), the inverse anchoring (round-trip, physical
back-solved cation, the unphysical-initial_ph compile guard), that the new acid slots
leave carbon conservation intact, that TA lands in the textbook band, and that the
derived pH tier is computed explicitly as ``plausible`` (never ``validated``).
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.analysis import ph_series, titratable_acidity_series
from fermentation.core import acidbase
from fermentation.core.chemistry import M_LACTIC, M_MALIC, M_TARTARIC
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon

#: Potassium molar mass [g/mol] — used only to state test 6's expected cation range
#: from first principles (K⁺ 1–2 g/L is the physical wine range), NOT read by the solver.
M_K = 39.0983


@pytest.fixture
def pset():
    """Real wine kinetic params + the shared pKa set, as a ParameterSet."""
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    """Resolved ``{name: float}`` map the solver hot-loop signature consumes."""
    return pset.resolve()


@pytest.fixture
def pka(params):
    return acidbase.build_pka_map(params)


def _wine_state(schema: StateSchema, **acids: float) -> FloatArray:
    """A wine state vector with arbitrary (pH-irrelevant) bulk values + given acids."""
    base: dict[str, float | list[float]] = {
        "X": 0.5,
        "S": [240.0],
        "E": 0.0,
        "N": 0.2,
        "T": 293.15,
        "CO2": 0.0,
    }
    base.update(acids)
    return schema.pack(base)


def _anchor_cation(
    pka: Mapping[str, tuple[float, ...]],
    tartaric_gpl: float,
    malic_gpl: float,
    target_ph: float,
) -> float:
    totals = {
        "tartaric": tartaric_gpl / M_TARTARIC,
        "malic": malic_gpl / M_MALIC,
        "lactic": 0.0,
    }
    return acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, target_ph)


# -- 1. HEADLINE: malic→lactic deacidification raises pH 0.1–0.3 ---------------


def test_headline_malic_to_lactic_raises_ph(params, pka):
    # A malic-rich must (the case where MLF matters): tartaric 4 / malic 4 g/L, anchored
    # to a measured pH 3.4. Full MLF converts all malic → lactic mole-for-mole; pH must
    # rise into the deacidification band. Settled empirically — a tartaric-heavy must
    # would land below 0.1; the fix for an out-of-band number is a more malic-rich must,
    # NOT widening the band (CLAUDE.md forbids weakening benchmark tests).
    schema = wine_schema()
    cation = _anchor_cation(pka, 4.0, 4.0, 3.4)
    y0 = _wine_state(schema, tartaric=4.0, malic=4.0, lactic=0.0, cation_charge=cation)
    ph0 = acidbase.ph_of_state(y0, schema, params)
    assert ph0 == pytest.approx(3.4, abs=1e-3)  # anchoring exact at t=0

    # all malic → lactic, conserving moles (g/L scales by the molar-mass ratio)
    lactic_gpl = (4.0 / M_MALIC) * M_LACTIC
    y1 = _wine_state(schema, tartaric=4.0, malic=0.0, lactic=lactic_gpl, cation_charge=cation)
    ph1 = acidbase.ph_of_state(y1, schema, params)

    delta = ph1 - ph0
    assert 0.1 <= delta <= 0.3, f"malic→lactic ΔpH {delta:.3f} outside [0.1, 0.3]"


# -- 2. the balance actually balances -----------------------------------------


def test_charge_residual_zero_at_solved_ph(pka):
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    assert acidbase.charge_residual(ph, totals, cation, 0.0, 0.0, pka) == pytest.approx(
        0.0, abs=1e-9
    )


# -- 3. monotonicity ----------------------------------------------------------


def test_more_acid_lowers_ph_more_cation_raises_ph(pka):
    base = {"tartaric": 4.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(base, 0.0, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(base, cation, 0.0, 0.0, pka)

    more_tartaric = {**base, "tartaric": 6.0 / M_TARTARIC}
    assert acidbase.solve_ph(more_tartaric, cation, 0.0, 0.0, pka) < ph  # more acid → lower pH

    assert acidbase.solve_ph(base, cation * 1.2, 0.0, 0.0, pka) > ph  # more cation → higher pH


# -- 4. smoothness / C¹ (guards a future in-loop BDF consumer) ----------------


def test_ph_is_smooth_in_acid(pka):
    cation = _anchor_cation(pka, 5.0, 3.0, 3.4)
    tartaric = np.linspace(3.0, 7.0, 41) / M_TARTARIC
    ph = np.array(
        [
            acidbase.solve_ph(
                {"tartaric": t, "malic": 3.0 / M_MALIC, "lactic": 0.0}, cation, 0.0, 0.0, pka
            )
            for t in tartaric
        ]
    )
    d1 = np.diff(ph)  # first difference (∝ dpH/d tartaric)
    d2 = np.diff(d1)  # second difference — small & sign-stable ⇒ no kink
    assert np.all(d1 < 0.0)  # strictly monotone
    assert np.max(np.abs(d2)) < 1e-2  # no derivative jump on a fine grid


# -- 4b. totality: solve_ph clamps a non-physiological probe cation, never raises (D-46) --
# BDF's num_jac perturbs the ``cation_charge`` state slot far outside its ~0.03 mol/L range
# while building the Jacobian, which can push ``charge_residual`` positive (or negative) across
# the whole [0, 14] bracket. ``solve_ph`` must stay a TOTAL, bounded function and clamp to the
# window rather than let ``brentq`` throw "f(a) and f(b) must have different signs". The three
# Brett integration tests only catch this incidentally (a 120-day run happens to drive the
# probe there); these pin it at the function level so a refactor cannot silently un-total it.


def test_solve_ph_clamps_huge_probe_cation_to_14(pka):
    # A probe cation two orders of magnitude above the physical ~0.03 mol/L: no acid load can
    # neutralise it, so the electroneutral pH lies above the window ⇒ clamp to 14, not raise.
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    assert acidbase.solve_ph(totals, 3.81, 0.0, 0.0, pka) == 14.0


def test_solve_ph_clamps_negative_probe_cation_to_0(pka):
    # The mirror probe: a large acid load with a strongly NEGATIVE strong-cation charge is
    # net-negative even fully protonated ⇒ electroneutral pH below the window ⇒ clamp to 0.
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    assert acidbase.solve_ph(totals, -2.0, 0.0, 0.0, pka) == 0.0


def test_solve_ph_physiological_cation_falls_through_to_brentq(pka):
    # The untouched path: a physiological cation still returns an interior root, unclamped and
    # bit-for-bit identical to the brentq result (the clamp branches are never taken).
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    assert 0.0 < ph < 14.0
    assert ph == pytest.approx(3.4, abs=1e-6)  # inverts solve_cation_charge exactly


# -- 5. round-trip: a compiled scenario reproduces its measured initial_ph -----
# NB tautological w.r.t. the g/L→mol/L factor (solve_cation_charge / solve_ph are
# inverses applying the same conversion, so a unit bug cancels). Test 6 is the guard.


def _wine_scenario(**initial_extra: float) -> Scenario:
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5}
    initial.update(initial_extra)
    return Scenario(
        name="wine-ph",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=14.0,
    )


def test_compiled_scenario_reproduces_initial_ph():
    compiled = compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4))
    ph0 = acidbase.ph_of_state(compiled.y0, compiled.schema, compiled.param_values)
    assert ph0 == pytest.approx(3.4, abs=1e-3)


# -- 6. THE UNITS GUARD: back-solved cation is physical (K⁺ 1–2 g/L) -----------


def test_back_solved_cation_is_physical():
    # Independent of the solver's own arithmetic: a textbook must (TA ~6–9 g/L, pH 3.4)
    # is counter-charged by K⁺ ~1–2 g/L, i.e. ~25.6–51.2 meq/L (1–2 g/L ÷ 39.0983).
    # A g/L↔mol/L factor error (which the round-trip can't see) lands the cation orders
    # of magnitude outside this band, so this is the clean catch.
    compiled = compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4))
    cation = compiled.schema.get(compiled.y0, "cation_charge")  # mol⁺/L
    meq_per_l = cation * 1000.0
    lo, hi = 1.0 / M_K * 1000.0, 2.0 / M_K * 1000.0  # ~25.6 .. 51.2 meq/L
    assert lo <= meq_per_l <= hi, (
        f"cation {meq_per_l:.1f} meq/L outside physical K⁺ {lo:.1f}–{hi:.1f}"
    )


# -- 7. compile guard: initial_ph below the acid load's intrinsic pH raises ----


def test_unphysical_initial_ph_raises_at_compile():
    # pH 2.0 with 6 g/L tartaric + 3 g/L malic needs a NEGATIVE strong cation — unphysical.
    with pytest.raises(ValueError, match="unphysical"):
        compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=2.0))


def test_initial_ph_without_pka_params_raises_clearly():
    # The explicit parameter_paths override is caller-controlled: a caller asking for
    # initial_ph but NOT including acidbase.yaml must get a clear, actionable error
    # (the missing-pKa KeyError is caught and re-raised), not a raw KeyError.
    wine_only = default_data_dir() / "wine_generic.yaml"  # deliberately omit acidbase.yaml
    with pytest.raises(ValueError, match="acidbase.yaml"):
        compile_scenario(
            _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4),
            parameter_paths=[wine_only],
        )


# -- 8. carbon conservation unchanged by the new acid slots -------------------


def test_carbon_conserved_with_constant_acids():
    compiled = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4), strict=True
    )
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.parameters["biomass_C_fraction"].value
    )
    # Acids are inert (no Process touches them) — a constant offset that drifts 0; the
    # rest of the carbon ledger still closes to machine precision with them present.
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-9, label="total carbon (with acids)")


# -- 9. TA lands in the textbook band -----------------------------------------


def test_titratable_acidity_in_band(params):
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4)
    y = _wine_state(schema, tartaric=6.0, malic=3.0, lactic=0.0, cation_charge=cation)
    ta = acidbase.titratable_acidity(y, schema, params)
    assert 6.0 <= ta <= 9.0, f"TA {ta:.2f} g/L tartaric-equiv outside the 6–9 band"


# -- 10. tier is computed explicitly as plausible, never validated ------------


def test_ph_tier_is_plausible(pset):
    # Computed explicitly (not inherited): the lowest pKa tier floored at plausible. The
    # pKa params are all plausible, and pH is never validated however good the pKa source.
    #
    # SCOPED TO WINE since D-179. The schema argument is not decoration: beer's lumped
    # peptide buffer is honestly speculative, so an unscoped call combines it in and reports
    # SPECULATIVE — for an acid this wine does not carry. Asking about wine's pH now means
    # saying "wine".
    assert acidbase.ph_tier(pset.tier_map(), wine_schema()) is Tier.PLAUSIBLE


def test_ph_tier_unscoped_is_conservative_not_wine(pset):
    """The no-schema default reports the WORST tier across all media, and that is deliberate.

    A caller with no schema does not know which medium it holds. Defaulting to wine's set
    would be back-compatible and would quietly over-report a beer pH as plausible; defaulting
    to the union can only ever under-claim. Prime directive #1 says an output's tier is the
    LOWEST of its inputs, so the conservative default is the correct one — but it is a real
    behaviour change from D-178, so it is pinned rather than left implicit.
    """
    assert acidbase.ph_tier(pset.tier_map()) is Tier.SPECULATIVE
    assert acidbase.ph_tier(pset.tier_map(), beer_schema()) is Tier.SPECULATIVE


# -- analysis series + the emergent Byp pH drift (the second demonstration) ----


def test_ph_series_drifts_down_as_byp_accumulates():
    # The acid slots are constant, but pH is NOT flat: Byp (core realised-yield diversion)
    # grows over the ferment and the charge balance reads it (include-by-reading), so with
    # the cation frozen at pitch the pH series drifts mildly DOWN — emergent, unscripted.
    compiled = compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4))
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    ph = ph_series(traj, compiled.param_values)
    ta = titratable_acidity_series(traj, compiled.param_values)
    assert ph[0] == pytest.approx(3.4, abs=1e-3)  # anchored at pitch
    drift = ph[0] - ph[-1]
    assert 0.02 <= drift <= 0.15, f"Byp-driven pH drift {drift:.3f} outside the expected mild fall"

    # The MUST (t=0) TA is the fidelity-grade value, in the 6-9 g/L band. The TA SERIES
    # then RISES as Byp accumulates (whole pool read as titratable diprotic succinic) —
    # an ACKNOWLEDGED upstream artifact (D-16/D-19 pool booking), backwards to real wine
    # (TA flat-to-declining during ferment). Pinned here as known/directional, NOT
    # fidelity: see acidbase.titratable_acidity caveat. Don't "fix" it by changing D-18.
    assert ta.shape == traj.t.shape and np.all(ta > 0.0)
    assert 6.0 <= ta[0] <= 9.0  # must value is trustworthy
    assert ta[-1] > ta[0]  # documents the artifact rise (do not assert end-of-ferment band)


# ======================================================================================
# D-182 — dissolved CO2 as carbonic acid in the charge balance
# ======================================================================================


def test_carbonic_is_not_a_member_of_any_acid_registry():
    """The structural claim, asserted because getting it wrong FAILS SILENTLY.

    Dissolved CO2 has no state slot: it is derived from ``CO2`` and ``T``. The acid registries
    are lists of SLOTS, and ``_totals_molar`` skips any name not in the schema — so putting
    ``carbonic`` in ``WINE_ACIDS``/``BEER_ACIDS`` would not raise, it would resolve to nothing
    and the term would quietly vanish from the balance while looking present in the source.
    That is why it travels as its own argument, the way ``Byp`` does.
    """
    for registry in (acidbase.WINE_ACIDS, acidbase.BEER_ACIDS, acidbase.ALL_ACIDS):
        assert acidbase.CARBONIC_KEY not in registry
    assert acidbase.CARBONIC_KEY not in wine_schema()
    assert acidbase.CARBONIC_KEY not in beer_schema()
    # ...but it IS a key in the pKa map, exactly like Byp, or charge_residual would KeyError.
    data = default_data_dir()
    pmap = acidbase.build_pka_map(
        load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml").resolve()
    )
    assert acidbase.CARBONIC_KEY in pmap and acidbase.BYP_KEY in pmap


def test_the_carbonic_parameters_are_all_in_the_sampled_read_set():
    """The D-160 guard, in the one place the derivation could not reach.

    ``PKA_PARAM_NAMES`` derives itself from the acid registries so a new SLOT cannot be
    forgotten. Carbonic is not a registry member, so the derivation would have skipped it and
    ``pKa_carbonic_1`` would sit outside ``PH_SYSTEM_READS`` — undeclared, therefore never
    sampled, therefore narrowing the reported spread of every pH-dependent output below what
    its own provenance justifies. The three solubility parameters are not pKas at all and
    could not have arrived by that route under any derivation.
    """
    for name in ("pKa_carbonic_1", *acidbase.CO2_SOLUBILITY_PARAMS):
        assert name in acidbase.PH_SYSTEM_READS, f"{name} would never be sampled"
    # Still disjoint from the SO2 set — a pH-only caller must not declare the binding Ks.
    assert not set(acidbase.PH_SYSTEM_READS) & set(acidbase.SO2_BINDING_READS)


def test_carbonic_is_monoprotic_and_never_reaches_the_polyprotic_branch(params):
    """pKa2 (10.3) is omitted, not deferred: at any beverage pH the carbonate fraction is
    below 1e-9, so a second proton would be arithmetic on nothing. Asserted rather than
    asserted-in-prose, because "add the second pKa for completeness" is the obvious edit."""
    assert acidbase.CARBONIC_AS_CO2.protons == 1
    # What the omission is worth, computed with a hand-built diprotic set rather than quoted,
    # at the HIGHEST pH any medium here starts from (a 5.65 wort, taken at 5.7 to be generous
    # to the second proton). The number that matters is the RATIO to the first proton, not the
    # absolute: an absolute threshold would silently pass on a pool that had grown.
    diprotic = (params["pKa_carbonic_1"], 10.3)
    h = 10.0 ** (-5.7)
    mono = acidbase.mean_charge(h, (params["pKa_carbonic_1"],))
    di = acidbase.mean_charge(h, diprotic)
    assert (di - mono) / mono < 1e-4, (
        f"the second proton adds {(di - mono) / mono:.2e} of the first proton's charge "
        f"({di - mono:.2e} per mole against {mono:.3f}); D-182 measured 4.6e-05"
    )


def test_saturation_is_the_printed_constant_at_the_reference_and_falls_with_temperature(params):
    """Henry's law with a van 't Hoff transfer — exact at T_ref, monotone decreasing above it.

    At the reference temperature the transfer is exactly 1, so the model uses the PRINTED
    in-beer constant unscaled (the ``T_ref`` idiom of the Arrhenius modifiers). Away from it
    solubility falls as temperature rises, which is the direction every source states.
    """
    ref = params["T_ref_co2_solubility"]
    assert acidbase.co2_saturation_gpl(ref, params) == pytest.approx(
        params["H_co2_beverage"], rel=1e-12
    ), "the transfer must be exactly 1 at the temperature the constant is quoted at"
    temps = [275.0, 283.15, 288.15, 293.15, 303.15]
    sats = [acidbase.co2_saturation_gpl(t, params) for t in temps]
    assert all(a > b for a, b in zip(sats[:-1], sats[1:], strict=True)), (
        "CO2 solubility must FALL as temperature rises"
    )
    # A lager at 10 C holds meaningfully more than an ale at 20 C — the reason this ships
    # temperature-dependent rather than as one number.
    assert sats[1] / sats[3] > 1.25


def test_saturation_is_total_under_a_jacobian_probe(params):
    """Bounded like ``solve_ph``'s bracket (D-46): ``num_jac`` perturbs the ``T`` slot, and a
    non-positive temperature would make 1/T blow up. The probe must not be able to raise."""
    assert acidbase.co2_saturation_gpl(0.0, params) == 0.0
    assert acidbase.co2_saturation_gpl(-5.0, params) == 0.0
    assert acidbase.co2_saturation_gpl(1e-6, params) >= 0.0


def test_the_evolved_slot_is_capped_at_saturation_and_is_zero_before_fermentation(params):
    """``min(evolved, saturation)`` — the two ends of the term's whole behaviour.

    The ``CO2`` slot is a cumulative evolved-gas integral reaching ~40 g/L in a beer and
    ~100 g/L in a wine, against a saturation near 2. Reading it directly would put twenty to
    fifty times too much carbonic acid in the balance, which is the single most likely
    misreading of this module.
    """
    schema = wine_schema()
    m_co2 = acidbase.CARBONIC_AS_CO2.molar_mass
    sat = acidbase.co2_saturation_gpl(293.15, params)

    unfermented = _wine_state(schema, CO2=0.0, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(unfermented, schema, params) == 0.0

    finished = _wine_state(schema, CO2=100.0, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(finished, schema, params) * m_co2 == pytest.approx(sat)

    # Below saturation the liquid holds ALL of it — the mass-conservative half of the min.
    early = _wine_state(schema, CO2=0.5, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(early, schema, params) * m_co2 == pytest.approx(0.5)

    # Solver undershoot on the evolved pool cannot make the term negative.
    undershot = _wine_state(schema, CO2=-1e-9, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(undershot, schema, params) == 0.0


def test_the_inverse_anchor_cannot_absorb_the_carbonic_term(pka, params):
    """THE STRUCTURAL REASON THIS TERM MOVES pH AT ALL, and the contrast that proves it.

    D-178 rejected malt phosphate on the ground that a species of CONSTANT charge is absorbed
    outright by the back-solved cation and is a near no-op. Dissolved CO2 escapes that only
    because it is ~0 in the must/wort where the anchor is taken and saturated afterwards.

    So this asserts BOTH halves: anchoring CO2-free and then finishing saturated moves the
    finished pH, while anchoring a state that ALREADY carries the same CO2 gives it back —
    the phosphate result, reproduced on purpose as the counterfactual.
    """
    totals = {"tartaric": 5.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    sat_molar = acidbase.co2_saturation_gpl(293.15, params) / acidbase.CARBONIC_AS_CO2.molar_mass

    # Real ordering: anchor with no CO2, finish with it.
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    ph_free = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    ph_saturated = acidbase.solve_ph(totals, cation, 0.0, sat_molar, pka)
    assert ph_free == pytest.approx(3.4, abs=1e-9)
    assert ph_saturated < ph_free, "dissolved CO2 must acidify"

    # Counterfactual: the same CO2 present when the anchor is taken. The cation absorbs it and
    # the finished pH is the anchor again — a no-op, exactly as D-178 measured for phosphate.
    absorbed_cation = acidbase.solve_cation_charge(totals, 0.0, sat_molar, pka, 3.4)
    ph_absorbed = acidbase.solve_ph(totals, absorbed_cation, 0.0, sat_molar, pka)
    assert ph_absorbed == pytest.approx(3.4, abs=1e-9)
    assert absorbed_cation > cation, "absorbing the anion charge takes MORE counter-cation"


def test_titratable_acidity_models_a_degassed_sample(params):
    """A titration is run on a DEGASSED sample, so carbonic must not count — and the check has
    to be that TA is BITWISE unchanged by the CO2 in the state, not merely close.

    OIV titratable acidity explicitly excludes carbonic (and sulfurous) acid; brewers and
    winemakers degas before titrating for exactly this reason. The asymmetry with
    ``ph_of_state``, which DOES read it, is the measurement convention rather than an
    oversight — and it is the kind of asymmetry a later reader tidies away.
    """
    schema = wine_schema()
    flat = _wine_state(schema, CO2=0.0, tartaric=5.0, malic=3.0, cation_charge=0.05)
    carbonated = _wine_state(schema, CO2=100.0, tartaric=5.0, malic=3.0, cation_charge=0.05)

    assert acidbase.titratable_acidity(flat, schema, params) == acidbase.titratable_acidity(
        carbonated, schema, params
    ), "TA must be bit-for-bit identical — it models a degassed sample"
    # ...while pH is NOT, which is what makes the exclusion a choice rather than a no-op.
    assert acidbase.ph_of_state(carbonated, schema, params) < acidbase.ph_of_state(
        flat, schema, params
    )


def test_the_model_sits_below_the_printed_volumes_band_and_the_head_free_scope_is_why(params):
    """The one cross-check that is independent of every shipped number.

    "The Chemistry of Beer" sec 6.7.1 prints, as an OUTCOME rather than a constant, that beer
    contains 1.2-1.7 volumes of CO2 after a normal nonpressurized fermentation (1 volume =
    0.196 % w/w). The model gives ~1.02 volumes at 15 C. It lands BELOW that band, and the
    reason is the scope stated in the YAML: this term is fixed at 1 atm CO2 partial pressure,
    while a real tank carries metres of hydrostatic head. Pinned as DIRECTIONAL — model at or
    under the printed low edge — never as a target to hit, because the two are not the same
    quantity. Reproducing 1.7 volumes at 1 atm would need a fermentation at about -2 C.
    """
    gpl_per_volume = 1.96  # 1 volume = 0.196 % by weight, printed in the same sentence
    for celsius in (10.0, 15.0, 20.0):
        volumes = acidbase.co2_saturation_gpl(273.15 + celsius, params) / gpl_per_volume
        assert 0.8 < volumes <= 1.2, (
            f"at {celsius:g} C the model holds {volumes:.2f} volumes; it must stay at or below "
            "the printed 1.2-1.7 band's low edge, which the head-free 1 atm scope explains. "
            "Above it would mean the model had acquired vessel pressure from somewhere."
        )


def test_the_apparent_pka_is_the_right_one_and_the_true_one_would_overstate_the_charge(params):
    """APPARENT, and the word is load-bearing: the constant is referenced to TOTAL dissolved
    CO2, of which under 1 % is actually H2CO3. Our pool IS total dissolved CO2. Using the TRUE
    carbonic pKa (3.7) against that pool would overstate the anion charge by hundreds of
    times — the single most plausible wrong edit to this parameter, so its size is measured."""
    h = 10.0 ** (-4.9)  # a finished beer
    apparent = acidbase.mean_charge(h, (params["pKa_carbonic_1"],))
    true_h2co3 = acidbase.mean_charge(h, (3.7,))
    assert true_h2co3 / apparent == pytest.approx(32.7, rel=0.05), (
        f"substituting the TRUE carbonic pKa would multiply the anion charge by "
        f"{true_h2co3 / apparent:.1f}x at beer pH; D-182 measured 32.7x. The two constants "
        "differ by 2.7 pKa units because under 1 % of dissolved CO2 is actually H2CO3 — if "
        "this ratio has collapsed toward 1, the shipped value is no longer the apparent one."
    )


def test_a_wine_moves_but_barely_and_the_geometry_is_why(params, pka):
    """The measured cost of shipping this MEDIUM-AGNOSTIC rather than beer-scoped.

    Wine's pH sits ~3 units BELOW the apparent pKa, so under 0.1 % of its dissolved CO2 is
    dissociated, against ~4.5 % at beer's pH. That is the identical geometric argument D-178
    used to REJECT malt phosphate for beer, pointing the other way — which is what makes
    "medium-agnostic" a physics statement rather than a convenience.
    """
    totals = {"tartaric": 5.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    sat_molar = acidbase.co2_saturation_gpl(293.15, params) / acidbase.CARBONIC_AS_CO2.molar_mass
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    shift = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka) - acidbase.solve_ph(
        totals, cation, 0.0, sat_molar, pka
    )
    assert 0.0 < shift < 0.005, (
        f"the carbonic term moved a wine by {shift:.4f} pH; D-182 measured ~0.0007. It is "
        "small BY GEOMETRY, and a large value here means either the pKa or the saturation "
        "has moved somewhere it should not have."
    )
    # The dissociated fractions that produce that asymmetry, stated as the check itself.
    frac_wine = acidbase.mean_charge(10.0 ** (-3.33), (params["pKa_carbonic_1"],))
    frac_beer = acidbase.mean_charge(10.0 ** (-4.96), (params["pKa_carbonic_1"],))
    assert frac_wine < 0.001 < 0.02 < frac_beer < 0.10
