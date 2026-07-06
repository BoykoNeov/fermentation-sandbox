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
from fermentation.core.media import wine_schema
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
    return acidbase.solve_cation_charge(totals, 0.0, pka, target_ph)


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
    cation = acidbase.solve_cation_charge(totals, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(totals, cation, 0.0, pka)
    assert acidbase.charge_residual(ph, totals, cation, 0.0, pka) == pytest.approx(0.0, abs=1e-9)


# -- 3. monotonicity ----------------------------------------------------------


def test_more_acid_lowers_ph_more_cation_raises_ph(pka):
    base = {"tartaric": 4.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(base, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(base, cation, 0.0, pka)

    more_tartaric = {**base, "tartaric": 6.0 / M_TARTARIC}
    assert acidbase.solve_ph(more_tartaric, cation, 0.0, pka) < ph  # more acid → lower pH

    assert acidbase.solve_ph(base, cation * 1.2, 0.0, pka) > ph  # more cation → higher pH


# -- 4. smoothness / C¹ (guards a future in-loop BDF consumer) ----------------


def test_ph_is_smooth_in_acid(pka):
    cation = _anchor_cation(pka, 5.0, 3.0, 3.4)
    tartaric = np.linspace(3.0, 7.0, 41) / M_TARTARIC
    ph = np.array(
        [
            acidbase.solve_ph(
                {"tartaric": t, "malic": 3.0 / M_MALIC, "lactic": 0.0}, cation, 0.0, pka
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
    assert acidbase.solve_ph(totals, 3.81, 0.0, pka) == 14.0


def test_solve_ph_clamps_negative_probe_cation_to_0(pka):
    # The mirror probe: a large acid load with a strongly NEGATIVE strong-cation charge is
    # net-negative even fully protonated ⇒ electroneutral pH below the window ⇒ clamp to 0.
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    assert acidbase.solve_ph(totals, -2.0, 0.0, pka) == 0.0


def test_solve_ph_physiological_cation_falls_through_to_brentq(pka):
    # The untouched path: a physiological cation still returns an interior root, unclamped and
    # bit-for-bit identical to the brentq result (the clamp branches are never taken).
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(totals, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(totals, cation, 0.0, pka)
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
    assert acidbase.ph_tier(pset.tier_map()) is Tier.PLAUSIBLE


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
