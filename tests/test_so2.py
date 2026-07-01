"""SO₂ speciation: the pH-coupled molecular fraction (D-22) and free/bound split (D-28).

Ranked headline-first. The keystone payoff is ``test_headline_molecular_fraction_falls_with_ph``:
the antimicrobial *molecular* SO₂ fraction falls ~3× per 0.5 pH unit and lands on the
textbook ~6 % / 2 % / 0.6 % at pH 3.0 / 3.5 / 4.0 — the coupling the D-18 charge-balance
solver was built to make *emerge* ("dose SO₂ → speciation falls out of the current pH").
Sections 1–9 pin a real-world anchor (~0.8 mg/L molecular at 40 mg/L free, pH 3.5), the
neutral-fraction algebra, the compile/readout plumbing, the explicit ``plausible`` tier,
and — prime directive #3 — that SO₂ is **readout-only**: dosing it leaves pH and carbon
byte-for-byte unchanged (it is not in the charge balance and is carbon-free).

Section 10 adds the **D-28 free/bound split** now that acetaldehyde is real state (D-27):
the dosed slot is *total* SO₂ and free/bound are derived by the acetaldehyde-bisulfite
binding equilibrium. Its own headline is
``test_emergent_free_so2_dips_at_acetaldehyde_peak_then_recovers`` — the early acetaldehyde
peak transiently sequesters SO₂, crashing free/molecular, which recover as acetaldehyde is
reduced; and the regression anchor ``test_binding_recovers_d22_at_zero_acetaldehyde`` pins
that at acetaldehyde = 0 the split collapses to D-22 exactly (free == total).
"""

import numpy as np
import pytest

from fermentation.analysis import (
    bound_so2_series,
    free_so2_series,
    molecular_so2_series,
    ph_series,
)
from fermentation.core import acidbase
from fermentation.core.chemistry import M_MALIC, M_TARTARIC
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units.convert import gpl_to_mgl, mgl_to_gpl
from fermentation.validation import assert_conserved, total_carbon


@pytest.fixture
def pset():
    """Real wine kinetic params + the shared pKa set (incl. the sulfurous pair)."""
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    return pset.resolve()


@pytest.fixture
def so2_pka(params):
    """The sulfurous (pKa₁, pKa₂) tuple the molecular-SO₂ readout partitions on."""
    return tuple(params[n] for n in acidbase.SO2_PKA_PARAM_NAMES)


def _wine_state(schema: StateSchema, **slots: float) -> FloatArray:
    base: dict[str, float | list[float]] = {
        "X": 0.5, "S": [240.0], "E": 0.0, "N": 0.2, "T": 293.15, "CO2": 0.0,
    }  # fmt: skip
    base.update(slots)
    return schema.pack(base)


def _anchor_cation(pka, tartaric_gpl: float, malic_gpl: float, target_ph: float) -> float:
    totals = {
        "tartaric": tartaric_gpl / M_TARTARIC,
        "malic": malic_gpl / M_MALIC,
        "lactic": 0.0,
    }
    return acidbase.solve_cation_charge(totals, 0.0, pka, target_ph)


def _wine_scenario(**initial_extra: float) -> Scenario:
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5}
    initial.update(initial_extra)
    return Scenario(
        name="wine-so2",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=14.0,
    )


# -- 1. HEADLINE: molecular fraction falls ~3× per 0.5 pH, on the textbook curve ---


def test_headline_molecular_fraction_falls_with_ph(so2_pka):
    # The antimicrobial molecular SO₂ fraction at the canonical wine pH points. pKa₁≈1.81
    # gives the textbook ~6 % / 2 % / 0.6 % — the pH coupling the D-18 keystone exists for.
    f30 = acidbase.molecular_so2_fraction(3.0, so2_pka)
    f35 = acidbase.molecular_so2_fraction(3.5, so2_pka)
    f40 = acidbase.molecular_so2_fraction(4.0, so2_pka)
    assert f30 == pytest.approx(0.0606, abs=3e-3)
    assert f35 == pytest.approx(0.0200, abs=2e-3)
    assert f40 == pytest.approx(0.0064, abs=1e-3)
    # ~3× drop per 0.5 pH unit (the winemaker's rule of thumb), both intervals.
    assert 2.7 <= f30 / f35 <= 3.3
    assert 2.7 <= f35 / f40 <= 3.3


# -- 2. a real-world anchor beyond the bare fraction --------------------------


def test_microbial_stability_anchor(so2_pka):
    # The standard target is ~0.8 mg/L molecular SO₂; at pH 3.5 that needs ~40 mg/L free.
    free_gpl = mgl_to_gpl(40.0)
    molecular_gpl = free_gpl * acidbase.molecular_so2_fraction(3.5, so2_pka)
    assert gpl_to_mgl(molecular_gpl) == pytest.approx(0.8, abs=0.1)


# -- 3. monotonic: higher pH ⇒ less molecular (antimicrobial) SO₂ -------------


def test_molecular_fraction_strictly_decreasing_in_ph(so2_pka):
    ph = np.linspace(2.8, 4.2, 29)
    frac = np.array([acidbase.molecular_so2_fraction(p, so2_pka) for p in ph])
    assert np.all(np.diff(frac) < 0.0)  # strictly falls as pH rises
    assert np.all((frac > 0.0) & (frac < 1.0))


# -- 4. the neutral-fraction algebra (and that sulfite is negligible at wine pH) --


def test_neutral_fraction_diprotic_matches_monoprotic_at_wine_ph(so2_pka):
    # pKa₂≈7.2 is far above wine pH, so the diprotic h²/D collapses to the monoprotic
    # 1/(1+10^(pH−pKa₁)) there — carrying pKa₂ is for exactness, not a wine-pH effect.
    pka1 = so2_pka[0]
    for ph in (3.0, 3.4, 3.8):
        h = 10.0 ** (-ph)
        di = acidbase.neutral_fraction(h, so2_pka)
        mono = acidbase.neutral_fraction(h, (pka1,))
        # agree to ~1e-4: the gap IS the sulfite fraction (~6e-5 of the total), negligible
        # at wine pH — and di < mono (the extra Ka₁·Ka₂ denominator term) as it must.
        assert di == pytest.approx(mono, abs=1e-4)
        assert di < mono
        assert mono == pytest.approx(1.0 / (1.0 + 10.0 ** (ph - pka1)), abs=1e-12)


def test_neutral_fraction_rejects_triprotic():
    with pytest.raises(ValueError, match="mono- and diprotic"):
        acidbase.neutral_fraction(1e-3, (1.8, 7.2, 9.0))


# -- 5. molecular_so2 off a compiled scenario = free × fraction(solved pH) -----


def test_molecular_so2_pure_function(params):
    compiled = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.5, so2_total_mgl=40.0)
    )
    y0, schema = compiled.y0, compiled.schema
    # the dose was converted mg/L → g/L into the so2_total slot; at pitch acetaldehyde=0 so
    # free == total and molecular = total × fraction(pH) — the D-22 curve recovered exactly.
    assert schema.get(y0, "so2_total") == pytest.approx(mgl_to_gpl(40.0))

    ph0 = acidbase.ph_of_state(y0, schema, params)  # solved from the organic acids only
    so2_pka = tuple(params[n] for n in acidbase.SO2_PKA_PARAM_NAMES)
    expected = mgl_to_gpl(40.0) * acidbase.molecular_so2_fraction(ph0, so2_pka)
    assert acidbase.molecular_so2(y0, schema, compiled.param_values) == pytest.approx(expected)
    # ~0.8 mg/L molecular at the ~40 mg/L free / pH 3.5 stability target
    got_mgl = gpl_to_mgl(acidbase.molecular_so2(y0, schema, compiled.param_values))
    assert got_mgl == pytest.approx(0.8, abs=0.15)


# -- 6. zero when undosed or the slot is absent -------------------------------


def test_molecular_so2_zero_without_dose(params):
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4)
    y = _wine_state(schema, tartaric=6.0, malic=3.0, cation_charge=cation)  # so2_total → 0
    assert acidbase.molecular_so2(y, schema, params) == 0.0


def test_molecular_so2_zero_when_slot_absent(params):
    # Beer has no so2_total slot (D-22 is wine-only); the readout returns 0, not raises.
    schema = beer_schema()
    y = schema.pack(
        {"X": 0.5, "S": [60.0, 100.0, 20.0], "E": 0.0, "N": 0.2, "T": 293.15, "CO2": 0.0}
    )
    assert acidbase.molecular_so2(y, schema, params) == 0.0


# -- 7. tier is computed explicitly as plausible, never validated -------------


def test_molecular_so2_tier_is_plausible(pset):
    # Combines BOTH pKa sets (pH solver + sulfurous), floored at plausible: the readout
    # solves pH (every pH pKa) and partitions free SO₂ (the sulfurous pKa). All plausible.
    assert acidbase.molecular_so2_tier(pset.tier_map()) is Tier.PLAUSIBLE


# -- 8. PRIME DIRECTIVE #3: SO₂ is readout-only — leaves pH and carbon untouched ---


def test_dosing_so2_does_not_change_ph(params):
    # SO₂ is NOT in the charge balance (D-22/D-28), so dosing it cannot move the solved pH:
    # ph_of_state reads only the organic acids + cation + Byp, never so2_total.
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4)
    dry = _wine_state(schema, tartaric=6.0, malic=3.0, cation_charge=cation)
    dosed = _wine_state(
        schema, tartaric=6.0, malic=3.0, cation_charge=cation, so2_total=mgl_to_gpl(80.0)
    )
    assert acidbase.ph_of_state(dosed, schema, params) == acidbase.ph_of_state(dry, schema, params)


def test_so2_does_not_perturb_carbon_or_the_core_trajectory():
    # The strongest isolability claim: on a *shared* time grid, dosing SO₂ leaves every
    # other state column byte-identical and the pH series identical — and carbon still
    # closes (so2_free is carbon-free, weight 0). 236 → still green.
    t_eval = np.linspace(0.0, 14.0 * 24.0, 60)
    base = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4), strict=True
    )
    dosed = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4, so2_total_mgl=60.0),
        strict=True,
    )
    traj0 = simulate(base.process_set, base.param_values, base.y0, base.t_span_h, t_eval=t_eval)
    traj1 = simulate(dosed.process_set, dosed.param_values, dosed.y0, dosed.t_span_h, t_eval=t_eval)

    for name in base.schema.names:
        if name == "so2_total":
            continue  # the only slot that differs (constant 0 vs constant 0.06 g/L)
        assert np.allclose(traj0.series(name), traj1.series(name), rtol=1e-9, atol=1e-12), name
    assert np.allclose(
        ph_series(traj0, base.param_values), ph_series(traj1, dosed.param_values), atol=1e-12
    )
    carbon = total_carbon(
        dosed.schema, biomass_carbon_fraction=dosed.parameters["biomass_C_fraction"].value
    )
    assert_conserved(traj1, carbon, rtol=1e-6, atol=1e-9, label="total carbon (with SO₂)")


# -- 9. the analysis series tracks the (mildly drifting) pH with no scripting ---


def test_molecular_so2_series_tracks_ph_drift():
    # molecular SO₂ is recomputed off the solved pH at each column, so as Byp accrues and
    # pH drifts down (D-18 emergent), the molecular fraction drifts *up* — unscripted.
    compiled = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4, so2_total_mgl=50.0)
    )
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    mol = molecular_so2_series(traj, compiled.param_values)
    ph = ph_series(traj, compiled.param_values)
    assert mol.shape == traj.t.shape
    # per-column consistency with the scalar pure function
    assert mol[0] == pytest.approx(
        acidbase.molecular_so2(traj.y[:, 0], traj.schema, compiled.param_values)
    )
    # pH drifts down over the ferment ⇒ molecular SO₂ rises (more antimicrobial late).
    assert ph[-1] < ph[0]
    assert mol[-1] > mol[0]
    assert np.all(mol > 0.0)


# == 10. D-28: free/bound SO₂ split (acetaldehyde binding) ======================


def test_bisulfite_fraction_dominates_free_so2_at_wine_ph(so2_pka):
    # Bisulfite HSO₃⁻ is the reactive binder and the dominant free-SO₂ species at wine pH:
    # β ≈ 0.94–0.99, so the bisulfite-vs-total-free reference basis for K differs ≤ ~6%
    # (the D-28 provenance claim). neutral + bisulfite ≈ 1 (sulfite negligible below pKa₂).
    for ph in (3.0, 3.4, 3.5, 4.0):
        h = 10.0 ** (-ph)
        beta = acidbase.bisulfite_fraction(h, so2_pka)
        neutral = acidbase.neutral_fraction(h, so2_pka)
        assert 0.93 < beta < 0.995
        assert neutral + beta == pytest.approx(1.0, abs=1e-3)  # sulfite ~1e-4 at wine pH
    # monoprotic branch = α₁, and triprotic is rejected like the sibling fractions
    assert acidbase.bisulfite_fraction(1e-3, (1.81,)) == pytest.approx(
        10.0**-1.81 / (10.0**-1.81 + 1e-3)
    )
    with pytest.raises(ValueError, match="mono- and diprotic"):
        acidbase.bisulfite_fraction(1e-3, (1.8, 7.2, 9.0))


def test_binding_equilibrium_algebra_solves_and_conserves(so2_pka):
    # bound_so2_molar is pure algebra: solve (A−x)(C−x)β − Kx = 0 for the physical root.
    beta = acidbase.bisulfite_fraction(10.0**-3.4, so2_pka)
    k = 1.5e-6
    a, c = 9.0e-4, 7.8e-4  # ~40 mg/L acetaldehyde, ~50 mg/L SO₂ (mol/L)
    x = acidbase.bound_so2_molar(c, a, beta, k)
    assert 0.0 < x < min(a, c)  # cannot bind more than either pool holds
    # x actually satisfies the equilibrium it claims to solve
    assert (a - x) * (c - x) * beta - k * x == pytest.approx(0.0, abs=1e-14)
    # degenerate inputs → no binding (guards the brentq-free quadratic)
    assert acidbase.bound_so2_molar(0.0, a, beta, k) == 0.0
    assert acidbase.bound_so2_molar(c, 0.0, beta, k) == 0.0
    assert acidbase.bound_so2_molar(c, a, 0.0, k) == 0.0


def test_binding_recovers_d22_at_zero_acetaldehyde(params):
    # THE regression anchor: with no acetaldehyde, bound = 0 and free == total, so the whole
    # D-22 readout (molecular = total × fraction(pH)) is reproduced byte-for-byte — the
    # input-semantics change (free → total) is invisible at the dosing moment.
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.5)
    dose = mgl_to_gpl(40.0)
    y = _wine_state(schema, tartaric=6.0, malic=3.0, cation_charge=cation, so2_total=dose)
    spec = acidbase.speciate_so2(y, schema, params)
    assert spec.bound == 0.0
    assert spec.free == spec.total == pytest.approx(dose)
    ph = acidbase.ph_of_state(y, schema, params)
    pkas = tuple(params[n] for n in acidbase.SO2_PKA_PARAM_NAMES)
    assert spec.molecular == pytest.approx(dose * acidbase.molecular_so2_fraction(ph, pkas))


def test_acetaldehyde_sequesters_so2_near_stoichiometric(params):
    # With comparable molar acetaldehyde and SO₂, binding is near-stoichiometric: bound
    # approaches the smaller pool, free (and molecular) crash toward ~0 — the mechanism
    # behind the emergent dip. free + bound conserves total exactly.
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4)
    total = mgl_to_gpl(50.0)
    dry = _wine_state(schema, tartaric=6.0, malic=3.0, cation_charge=cation, so2_total=total)
    peak = _wine_state(
        schema, tartaric=6.0, malic=3.0, cation_charge=cation,
        so2_total=total, acetaldehyde=mgl_to_gpl(37.0),
    )  # fmt: skip
    s_dry = acidbase.speciate_so2(dry, schema, params)
    s_peak = acidbase.speciate_so2(peak, schema, params)
    assert s_dry.bound == 0.0 and s_dry.free == pytest.approx(total)
    assert s_peak.free < 0.05 * total  # free crashes at the peak (SO₂ ≈ limiting)
    assert s_peak.molecular < 0.05 * s_dry.molecular  # antimicrobial pool collapses
    # conservation: free + bound == total, at both states, exactly
    assert s_peak.free + s_peak.bound == pytest.approx(total)
    assert s_dry.free + s_dry.bound == pytest.approx(total)


def test_speciate_matches_scalar_wrappers(params):
    # speciate_so2 (one pH solve) and the scalar convenience wrappers agree.
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4)
    y = _wine_state(
        schema, tartaric=6.0, malic=3.0, cation_charge=cation,
        so2_total=mgl_to_gpl(60.0), acetaldehyde=mgl_to_gpl(30.0),
    )  # fmt: skip
    spec = acidbase.speciate_so2(y, schema, params)
    assert acidbase.bound_so2(y, schema, params) == pytest.approx(spec.bound)
    assert acidbase.free_so2(y, schema, params) == pytest.approx(spec.free)
    assert acidbase.molecular_so2(y, schema, params) == pytest.approx(spec.molecular)
    # molecular_so2_at_ph reuses a supplied pH and agrees with the full solve
    assert acidbase.molecular_so2_at_ph(y, schema, params, spec.ph) == pytest.approx(spec.molecular)


def test_emergent_free_so2_dips_at_acetaldehyde_peak_then_recovers():
    # THE D-28 headline: over a real ferment the early acetaldehyde peak transiently binds
    # SO₂ — free (analytically-measured) SO₂ crashes toward ~0 near the peak, then RECOVERS
    # to the dosed total as acetaldehyde is reduced back to ethanol (D-27). Unscripted: the
    # dosed total slot is constant; the dip emerges from the binding equilibrium tracking the
    # acetaldehyde state. free + bound conserves total at every column.
    compiled = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4, so2_total_mgl=50.0)
    )
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    free = free_so2_series(traj, compiled.param_values)
    bound = bound_so2_series(traj, compiled.param_values)
    acet = traj.series("acetaldehyde")
    total = mgl_to_gpl(50.0)

    assert free[0] == pytest.approx(total)  # no acetaldehyde at pitch ⇒ all free
    assert free.min() < 0.15 * total  # deep dip at the peak (SO₂ nearly all bound)
    assert free[-1] == pytest.approx(total, rel=1e-3)  # recovers as acetaldehyde clears
    # the free minimum coincides with the acetaldehyde maximum (the causal link)
    assert abs(int(np.argmin(free)) - int(np.argmax(acet))) <= 1
    assert bound.max() > 0.8 * total  # nearly all SO₂ is bound at the peak
    assert np.allclose(free + bound, total)  # conservation at every column


def test_speciation_tier_drops_with_binding_constant(pset, params):
    # The tier now folds in the binding constant K too (D-28): all plausible ⇒ plausible,
    # but a speculative K drags the whole speciation readout to speculative.
    assert acidbase.molecular_so2_tier(pset.tier_map()) is Tier.PLAUSIBLE
    tiers = dict(pset.tier_map())
    tiers[acidbase.SO2_BINDING_PARAM] = Tier.SPECULATIVE
    assert acidbase.molecular_so2_tier(tiers) is Tier.SPECULATIVE
