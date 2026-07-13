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
from fermentation.core.chemistry import M_ACETALDEHYDE, M_MALIC, M_SO2, M_TARTARIC
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


def test_so2_coupling_strands_acetaldehyde_but_spares_the_core_ferment():
    # D-47 RETIRED the "SO₂ is readout-only" invariant: dosing SO₂ now feeds back into the
    # acetaldehyde reduction (bound acetaldehyde is protected from ADH). This test pins the
    # *footprint* of that coupling on a shared time grid — it is confined to acetaldehyde:
    #   • acetaldehyde diverges order-unity — the undosed run clears it to ~0, the dosed run
    #     STRANDS a locked-in residual (~33 mg/L for a 60 mg/L dose);
    #   • every OTHER column (sugar, ethanol, biomass, CO₂, the acids) still agrees to a
    #     second-order ≤2e-3 — the only ripple is the borrowed-ethanol-carbon dip feeding the
    #     E→viability brake (the D-27 note), NOT a rewrite of the core ferment;
    #   • pH is unmoved to ~2e-6 — SO₂ is still NOT in the charge balance (D-22/D-28), it now
    #     couples only through acetaldehyde, and acetaldehyde is not a charge species;
    #   • carbon STILL closes (the one surviving invariant — the reduction only throttles the
    #     acetaldehyde→E transfer, it neither creates nor routes carbon).
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

    # acetaldehyde: the intended order-unity divergence (undosed clears; dosed strands).
    assert traj0.series("acetaldehyde")[-1] < 1e-6  # undosed: fully reduced back to ethanol
    assert traj1.series("acetaldehyde")[-1] > 0.02  # dosed: ~33 mg/L locked in by SO₂

    # every other column stays within the second-order E→viability ripple (was byte-identical
    # under the retired readout-only invariant; now ≤1e-3 of each column's own scale, driven by
    # the tiny borrowed-C dip). Compared as a fraction of the column scale, NOT elementwise: a
    # second-order shift in the timing of sugar exhaustion makes pointwise relative diffs blow up
    # on steep fronts (S, E) even when the trajectories are everywhere within a thousandth.
    for name in base.schema.names:
        if name in ("so2_total", "acetaldehyde"):
            continue
        a, b = traj0.series(name), traj1.series(name)
        scale = np.max(np.abs(b))
        assert np.max(np.abs(a - b)) <= 1e-3 * scale + 1e-6, name
    # pH essentially unmoved: SO₂ couples via acetaldehyde, which carries no charge.
    assert np.allclose(
        ph_series(traj0, base.param_values), ph_series(traj1, dosed.param_values), atol=1e-4
    )
    carbon = total_carbon(
        dosed.schema, biomass_carbon_fraction=dosed.parameters["biomass_C_fraction"].value
    )
    assert_conserved(traj1, carbon, rtol=1e-6, atol=1e-9, label="total carbon (with SO₂)")


# -- 9. the analysis series tracks the (mildly drifting) pH with no scripting ---


def test_molecular_so2_series_falls_as_stranded_acetaldehyde_binds_free_so2():
    # Two competing effects set the molecular-SO₂ trajectory, and D-47 flips which one wins:
    #   • the molecular *fraction* still rises as pH drifts down (the D-18/D-22 pH coupling —
    #     lower pH ⇒ more undissociated antimicrobial SO₂·H₂O), and
    #   • but the *free* SO₂ pool is now chronically depressed, because the acetaldehyde the
    #     yeast stranded (D-47: protected from ADH) stays bound to most of the dose.
    # Free falls ~4.5× while the fraction rises only ~1.16×, so molecular SO₂ (= free × fraction)
    # nets DOWN over the run — the opposite of the readout-only era, where free recovered to the
    # dosed total and molecular tracked the rising fraction upward. Unscripted: both effects fall
    # out of the state, not a schedule.
    compiled = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4, so2_total_mgl=50.0)
    )
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    mol = molecular_so2_series(traj, compiled.param_values)
    free = free_so2_series(traj, compiled.param_values)
    ph = ph_series(traj, compiled.param_values)
    so2_pka = tuple(compiled.param_values[n] for n in acidbase.SO2_PKA_PARAM_NAMES)
    assert mol.shape == traj.t.shape
    # per-column consistency with the scalar pure function
    assert mol[0] == pytest.approx(
        acidbase.molecular_so2(traj.y[:, 0], traj.schema, compiled.param_values)
    )
    assert ph[-1] < ph[0]  # pH still drifts down (Byp accrual, D-18)
    # the pH-driven molecular *fraction* still rises (the D-22 coupling survives) …
    assert acidbase.molecular_so2_fraction(ph[-1], so2_pka) > acidbase.molecular_so2_fraction(
        ph[0], so2_pka
    )
    # … but free SO₂ ends far below the dose (stranded acetaldehyde holds it bound) …
    assert free[-1] < 0.4 * free[0]
    # … so absolute molecular SO₂ nets DOWN — the free depression dominates (the D-47 flip).
    assert mol[-1] < mol[0]
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
    # bound_so2_molar is pure algebra: solve (A−x)(C−x)β − Kx = 0 for the physical root — the
    # n=1 case of the D-51 multi-carbonyl solver (a 1-tuple of (molar, K)).
    beta = acidbase.bisulfite_fraction(10.0**-3.4, so2_pka)
    k = 1.5e-6
    a, c = 9.0e-4, 7.8e-4  # ~40 mg/L acetaldehyde, ~50 mg/L SO₂ (mol/L)
    (x,) = acidbase.bound_so2_molar(c, ((a, k),), beta)
    assert 0.0 < x < min(a, c)  # cannot bind more than either pool holds
    # x actually satisfies the equilibrium it claims to solve
    assert (a - x) * (c - x) * beta - k * x == pytest.approx(0.0, abs=1e-14)
    # degenerate inputs → no binding (guards the brentq-free early-exit paths)
    assert acidbase.bound_so2_molar(0.0, ((a, k),), beta) == (0.0,)
    assert acidbase.bound_so2_molar(c, ((0.0, k),), beta) == (0.0,)
    assert acidbase.bound_so2_molar(c, ((a, k),), 0.0) == (0.0,)


# == D-51: the coupled multi-carbonyl equilibrium (pyruvate + α-ketoglutarate compete too) =====


def test_multi_carbonyl_reduces_exactly_to_the_single_carbonyl_form(so2_pka):
    # THE regression anchor: with only one carbonyl present (the other two molar amounts are 0),
    # the N-carbonyl solver must reproduce the ORIGINAL D-28 quadratic root exactly — a
    # keto-acid-pool-off run is byte-for-byte the pre-D-51 form.
    beta = acidbase.bisulfite_fraction(10.0**-3.4, so2_pka)
    c, a = 7.8e-4, 9.0e-4
    k_acet = 1.5e-6
    x_multi, x_pyr_zero, x_akg_zero = acidbase.bound_so2_molar(
        c, ((a, k_acet), (0.0, 5.55e-4), (0.0, 1.4e-4)), beta
    )
    assert x_pyr_zero == 0.0 and x_akg_zero == 0.0  # zero-molar carbonyls contribute nothing
    qb = -(beta * (a + c) + k_acet)
    disc = qb * qb - 4.0 * beta * beta * a * c
    x_single = (-qb - disc**0.5) / (2.0 * beta)
    assert x_multi == pytest.approx(min(max(x_single, 0.0), min(a, c)), abs=1e-12)


def test_multi_carbonyl_competition_conserves_and_binds_less_acetaldehyde(so2_pka):
    # With pyruvate/α-KG present alongside acetaldehyde, all three compete for ONE shared
    # bisulfite pool: acetaldehyde's bound share must be SMALLER than if it had the pool to
    # itself (some SO₂ goes to the competitors instead), and the three bound amounts + the
    # remaining free SO₂ must conserve the total exactly.
    beta = acidbase.bisulfite_fraction(10.0**-3.4, so2_pka)
    total = 7.8e-4  # ~50 mg/L SO₂ (mol/L)
    acet, pyr, akg = 5.8e-4, 3.4e-4, 1.4e-4  # ~ D-47/D-49/D-50 finished-wine molar levels
    k_acet, k_pyr, k_akg = 1.5e-6, 5.55e-4, 1.4e-4
    (x_acet_alone,) = acidbase.bound_so2_molar(total, ((acet, k_acet),), beta)
    x_acet, x_pyr, x_akg = acidbase.bound_so2_molar(
        total, ((acet, k_acet), (pyr, k_pyr), (akg, k_akg)), beta
    )
    assert x_acet < x_acet_alone  # competitors soak up some of what acetaldehyde would have taken
    assert x_pyr > 0.0 and x_akg > 0.0  # weaker binders still capture a real, non-zero share
    free = total - (x_acet + x_pyr + x_akg)
    assert free >= 0.0
    assert x_acet + x_pyr + x_akg + free == pytest.approx(total)  # conservation, exactly
    # each species still satisfies ITS OWN adduct equilibrium at the shared reactive bisulfite h
    h = beta * free
    for x, a, k in ((x_acet, acet, k_acet), (x_pyr, pyr, k_pyr), (x_akg, akg, k_akg)):
        assert (a - x) * h - k * x == pytest.approx(0.0, abs=1e-12)


def test_multi_carbonyl_order_independent_and_clamped(so2_pka):
    # The return tuple tracks input ORDER (not a fixed acetaldehyde/pyruvate/α-KG position), and
    # every entry stays within its own physical bound regardless of how the others are sized.
    beta = acidbase.bisulfite_fraction(10.0**-3.4, so2_pka)
    total = 5.0e-4
    carbonyls = ((3.0e-4, 5.55e-4), (1.0e-6, 1.5e-6), (2.0e-4, 1.4e-4))
    x_pyr, x_acet, x_akg = acidbase.bound_so2_molar(total, carbonyls, beta)
    reordered = ((1.0e-6, 1.5e-6), (2.0e-4, 1.4e-4), (3.0e-4, 5.55e-4))
    x_acet2, x_akg2, x_pyr2 = acidbase.bound_so2_molar(total, reordered, beta)
    assert (x_acet2, x_pyr2, x_akg2) == pytest.approx((x_acet, x_pyr, x_akg))
    for x, (a, _) in zip((x_pyr, x_acet, x_akg), carbonyls, strict=True):
        assert 0.0 <= x <= min(a, total)


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


def test_keto_acid_pools_widen_the_bound_so2_and_free_more_acetaldehyde(params):
    # THE D-51 wiring headline, at the state level (not just the pure-algebra tests above):
    # adding pyruvate/α-KG to an otherwise-identical state must (a) bind MORE total SO₂ than
    # acetaldehyde alone would, and (b) leave acetaldehyde itself LESS protected (free_acetaldehyde
    # rises) — competitors soak up part of the shared bisulfite pool acetaldehyde used to have to
    # itself. A keto-acid-pool-off state (the slots simply absent from ``slots``, defaulting to 0
    # via ``_wine_state``) is unaffected — this is the isolability the D-49/D-50 build promised.
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4)
    total = mgl_to_gpl(50.0)
    acet = mgl_to_gpl(60.0)
    acet_only = _wine_state(
        schema, tartaric=6.0, malic=3.0, cation_charge=cation, so2_total=total, acetaldehyde=acet
    )
    with_keto = _wine_state(
        schema, tartaric=6.0, malic=3.0, cation_charge=cation, so2_total=total, acetaldehyde=acet,
        pyruvate=mgl_to_gpl(30.0), alpha_ketoglutarate=mgl_to_gpl(20.0),
    )  # fmt: skip
    s_acet_only = acidbase.speciate_so2(acet_only, schema, params)
    s_with_keto = acidbase.speciate_so2(with_keto, schema, params)
    assert s_with_keto.bound > s_acet_only.bound  # keto acids add to the bound total
    assert s_with_keto.free < s_acet_only.free  # so free (and hence molecular) is further depressed
    assert s_with_keto.free + s_with_keto.bound == pytest.approx(total)  # still conserves exactly

    ph = acidbase.ph_of_state(with_keto, schema, params)
    free_acet_only = acidbase.free_acetaldehyde(acet_only, schema, params, ph)
    free_with_keto = acidbase.free_acetaldehyde(with_keto, schema, params, ph)
    assert free_with_keto > free_acet_only  # acetaldehyde itself is LESS bound (more ADH-reducible)
    assert free_with_keto <= mgl_to_gpl(60.0)  # never exceeds the total acetaldehyde present


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


def test_emergent_free_so2_dips_then_locks_in_a_stranded_acetaldehyde_residual():
    # THE D-47 headline (supersedes the D-28 "free recovers to total"): over a real ferment the
    # acetaldehyde peak binds SO₂ near-stoichiometrically — free (analytically-measured) SO₂
    # crashes toward ~0 near the peak — and because bound acetaldehyde is now PROTECTED FROM ADH
    # (D-47), the acetaldehyde never fully clears: it is LOCKED IN. Free SO₂ therefore only
    # *partially* recovers (the excess acetaldehyde above the SO₂ molar pool reduces away, but the
    # ~stoichiometric remainder stays bound), settling at ~0.22× the dose. The stranded
    # acetaldehyde ends at ~0.78 mol per mol SO₂ — the sub-to-near-stoichiometric field regime the
    # literature reports. Unscripted: the dosed total slot is constant; everything emerges from the
    # binding equilibrium tracking the acetaldehyde state. free + bound conserves total everywhere.
    compiled = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4, so2_total_mgl=50.0)
    )
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    free = free_so2_series(traj, compiled.param_values)
    bound = bound_so2_series(traj, compiled.param_values)
    acet = traj.series("acetaldehyde")
    total = mgl_to_gpl(50.0)

    assert free[0] == pytest.approx(total)  # no acetaldehyde at pitch ⇒ all free
    assert free.min() < 0.05 * total  # deep dip at the peak (SO₂ nearly all bound)
    # free PARTIALLY recovers from the dip (excess acetaldehyde clears) but locks in far below the
    # dose — the coupling's signature, NOT the readout-only-era recovery to the full total.
    assert free.min() < free[-1] < 0.4 * total
    assert free[-1] > 5.0 * free.min()  # materially higher than the peak dip …
    assert bound[-1] > 0.5 * total  # … yet most of the dose stays bound (locked in)
    # the free minimum coincides with the acetaldehyde maximum (the causal link)
    assert abs(int(np.argmin(free)) - int(np.argmax(acet))) <= 1
    assert bound.max() > 0.9 * total  # nearly all SO₂ is bound at the peak
    assert np.allclose(free + bound, total)  # conservation at every column
    # the stranded acetaldehyde ends near-stoichiometric with the SO₂ pool (literature: ~0.5–1:1)
    stranded_molar = acet[-1] / M_ACETALDEHYDE
    so2_molar = total / M_SO2
    assert 0.5 < stranded_molar / so2_molar < 1.0
    assert gpl_to_mgl(acet[-1]) > 20.0  # ~27 mg/L locked in (a sensorily-relevant residual)


def test_speciation_tier_drops_with_binding_constant(pset, params):
    # The tier now folds in the binding constant K too (D-28): all plausible ⇒ plausible,
    # but a speculative K drags the whole speciation readout to speculative.
    assert acidbase.molecular_so2_tier(pset.tier_map()) is Tier.PLAUSIBLE
    tiers = dict(pset.tier_map())
    tiers[acidbase.SO2_BINDING_PARAM] = Tier.SPECULATIVE
    assert acidbase.molecular_so2_tier(tiers) is Tier.SPECULATIVE


# == 11. D-82: the reversible SO₂/pH masking readout (anthocyanin coloured fraction) ============
# The pure competitive-denominator scalar behind analysis.observed_color_series — carbinol
# (hydration) and bisulfite adduct as PARALLEL drains of the flavylium pool, one denominator.


def test_coloured_fraction_competitive_denominator_form():
    # coloured = 1 / (1 + K_h/h + K·[HSO₃⁻]) — verified against the explicit expression, and shown
    # NOT to equal the product form [h/(h+K_h)]·[1/(1+K·B)] (which carries a spurious cross-term
    # physically implying bisulfite bleaches the colourless carbinol — it cannot).
    pk_h, k_bleach = 2.6, 2.5e4
    k_h = 10.0 ** (-pk_h)
    for ph, b in ((3.0, 0.0), (3.4, 3.0e-4), (3.8, 1.0e-4)):
        h = 10.0 ** (-ph)
        got = acidbase.anthocyanin_coloured_fraction(h, b, pk_h, k_bleach)
        assert got == pytest.approx(1.0 / (1.0 + k_h / h + k_bleach * b))
        product = (h / (h + k_h)) * (1.0 / (1.0 + k_bleach * b))
        if b > 0.0:  # the forms diverge only when BOTH drains are active
            assert got > product  # the competitive form retains more colour (no double-masking)


def test_coloured_fraction_pure_ph_limit_matches_neutral_fraction():
    # With no bisulfite (B=0) the bleaching term drops and the coloured fraction collapses to the
    # monoprotic flavylium⇌carbinol share h/(h+K_h) — i.e. neutral_fraction at the hydration pK,
    # the same textbook shape the molecular-SO₂ readout uses (a cross-check on the pH limb).
    pk_h = 2.6
    for ph in (3.0, 3.4, 3.8, 4.2):
        h = 10.0 ** (-ph)
        assert acidbase.anthocyanin_coloured_fraction(h, 0.0, pk_h, 2.5e4) == pytest.approx(
            acidbase.neutral_fraction(h, (pk_h,))
        )


def test_coloured_fraction_low_ph_is_redder_and_so2_bleaches():
    # The two sourced directions: LOWER pH ⇒ more red flavylium (acidify to brighten), and MORE free
    # bisulfite ⇒ less apparent colour (SO₂ bleaching). Both strictly monotone; bounded (0, 1].
    pk_h, k_bleach = 2.6, 2.5e4
    ph_grid = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0]
    fr = [
        acidbase.anthocyanin_coloured_fraction(10.0 ** (-p), 0.0, pk_h, k_bleach) for p in ph_grid
    ]
    assert all(a > b for a, b in zip(fr, fr[1:], strict=False))  # strictly falls as pH rises
    b_grid = [0.0, 1.0e-4, 3.0e-4, 1.0e-3]
    h = 10.0 ** (-3.4)
    fb = [acidbase.anthocyanin_coloured_fraction(h, b, pk_h, k_bleach) for b in b_grid]
    assert all(a > b for a, b in zip(fb, fb[1:], strict=False))  # strictly falls as bisulfite rises
    assert fb[0] == pytest.approx(acidbase.anthocyanin_coloured_fraction(h, 0.0, pk_h, k_bleach))
    for f in (*fr, *fb):
        assert 0.0 < f <= 1.0
    # At wine pH 3.4 only a minority of monomeric anthocyanin is red even with NO SO₂ (~0.14) — the
    # flavylium-minority textbook result the pH mask encodes.
    assert acidbase.anthocyanin_coloured_fraction(h, 0.0, pk_h, k_bleach) == pytest.approx(
        0.14, abs=0.02
    )
