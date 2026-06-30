"""Tests for the Tier-2 aroma-byproduct Processes (decision D-18/D-19).

:class:`EsterSynthesis` and :class:`FuselAlcoholsEhrlich` fill the ``esters``/
``fusels`` pools and, under **option a1 (decision D-19)**, route that carbon *out of
``S``* so the pools are real carbon-accounted state. These tests pin the closed-form
derivatives (including the exact sugar draw) and their guards, prove the properties
the beat requires of every byproduct Process — **monotone-increasing in temperature**,
**carbon-routed-from-sugar** (per-RHS the sugar carbon removed equals the carbon
deposited in the pool, so ``total_carbon`` closes), and **isolable** (the core is
built without them; enabling them perturbs only ``dS`` — never ``dX``/``dN``/``dE``/
``dCO2`` — by the trace sugar they consume) — and check the nitrogen gate on the
Ehrlich pathway plus tier propagation (the fusel form is speculative; the ester form
is plausible but its placeholder rate params cap its output at speculative).

The run-integrated "cleaner when colder" direction (lower T ⇒ fewer total
esters+fusels) is the *benchmark's* job (``test_lower_temperature_is_slower_but_\
cleaner``); these unit tests cover the per-Process mechanics it rests on.
"""

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics import (
    EsterSynthesis,
    FuselAlcoholsEhrlich,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
    arrhenius_factor,
)
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon

#: Representative species the pools book against (mirrors the Process constants).
_ESTER_C = carbon_mass_fraction("ethyl_acetate")
_FUSEL_C = carbon_mass_fraction("isoamyl_alcohol")
_GLUCOSE_C = carbon_mass_fraction("glucose")


@pytest.fixture
def store():
    # Real (sourced) wine parameters; now carry the byproduct placeholder constants.
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x: float = 2.0,
    s: float = 200.0,
    e: float = 0.0,
    n: float = 0.1,
    t: float = 293.15,
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": t, "CO2": 0.0})


# -- metadata -----------------------------------------------------------------


def test_ester_metadata():
    p = EsterSynthesis()
    assert p.name == "ester_synthesis"
    assert p.tier is Tier.PLAUSIBLE  # warmth-favoured, flux-coupled: standard form
    # Touches its own pool AND S — the ester carbon is routed from sugar (a1, D-19);
    # never E/CO2 (uptake's ethanol/CO2 production is left untouched).
    assert set(p.touches) == {"esters", "S"}
    assert set(p.reads) == {"k_ester", "K_sugar_uptake", "E_a_esters", "T_ref"}


def test_fusel_metadata():
    p = FuselAlcoholsEhrlich()
    assert p.name == "fusel_alcohols_ehrlich"
    # Speculative *form*: the nitrogen dependence is a knowingly-monotone simplification.
    assert p.tier is Tier.SPECULATIVE
    assert set(p.touches) == {"fusels", "S"}  # fusel carbon routed from sugar (a1, D-19)
    assert set(p.reads) == {"k_fusel", "K_sugar_uptake", "K_n", "E_a_fusels", "T_ref"}


# -- ester closed form & guards -----------------------------------------------


def test_ester_derivative_matches_closed_form(params):
    schema = wine_schema()
    x, s, t = 2.0, 200.0, 293.15
    y = _wine_y0(schema, x=x, s=s, t=t)
    d = EsterSynthesis().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    f_t = arrhenius_factor(t, params["E_a_esters"], params["T_ref"])
    rate = params["k_ester"] * flux * f_t
    assert schema.get(d, "esters") == pytest.approx(rate)
    # The ester carbon is routed from sugar (a1, D-19): one slot, so dS removes exactly
    # the ester carbon converted back to grams of glucose. Carbon balances per-RHS.
    assert schema.get(d, "S") == pytest.approx(-rate * _ESTER_C / _GLUCOSE_C)
    assert schema.get(d, "S") * _GLUCOSE_C == pytest.approx(-schema.get(d, "esters") * _ESTER_C)
    # Nothing else moves — ethanol/CO2 production is left to the uptake Process.
    for var in ("X", "E", "N", "CO2", "fusels"):
        assert schema.get(d, var) == 0.0


def test_ester_factor_is_one_at_reference_temperature(params):
    # At T_ref the embedded Arrhenius factor is exactly 1, so the rate is the bare
    # flux term times k_ester (the rate constants are anchored at T_ref).
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = _wine_y0(schema, x=x, s=s, t=params["T_ref"])
    d = EsterSynthesis().derivatives(0.0, y, schema, params)
    flux = x * (s / (params["K_sugar_uptake"] + s))
    assert schema.get(d, "esters") == pytest.approx(params["k_ester"] * flux)


def test_ester_rises_with_temperature(params):
    # The defining directional property: warmer ⇒ more esters per unit flux.
    schema = wine_schema()
    cold = EsterSynthesis().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = EsterSynthesis().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(warm, "esters") > schema.get(cold, "esters") > 0.0


def test_ester_scales_with_fermentative_flux(params):
    # Coupled to the biomass-catalysed sugar flux: more biomass ⇒ proportionally
    # more ester synthesis (the flux is linear in X).
    schema = wine_schema()
    r1 = EsterSynthesis().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = EsterSynthesis().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "esters") == pytest.approx(2.0 * schema.get(r1, "esters"))


def test_ester_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    assert EsterSynthesis().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)[
        schema.slice("esters")
    ] == pytest.approx(0.0)
    assert EsterSynthesis().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)[
        schema.slice("esters")
    ] == pytest.approx(0.0)


def test_ester_negative_excursion_does_not_produce(params):
    # A solver undershoot (S<0 or X<0) must not flip the clamp and create esters.
    schema = wine_schema()
    assert np.array_equal(
        EsterSynthesis().derivatives(0.0, _wine_y0(schema, s=-1e-6), schema, params),
        schema.zeros(),
    )
    assert np.array_equal(
        EsterSynthesis().derivatives(0.0, _wine_y0(schema, x=-1e-6), schema, params),
        schema.zeros(),
    )


# -- fusel closed form, nitrogen gate & guards --------------------------------


def test_fusel_derivative_matches_closed_form(params):
    schema = wine_schema()
    x, s, n, t = 2.0, 200.0, 0.1, 293.15
    y = _wine_y0(schema, x=x, s=s, n=n, t=t)
    d = FuselAlcoholsEhrlich().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    gate = n / (params["K_n"] + n)
    f_t = arrhenius_factor(t, params["E_a_fusels"], params["T_ref"])
    rate = params["k_fusel"] * flux * gate * f_t
    assert schema.get(d, "fusels") == pytest.approx(rate)
    # Fusel carbon routed from sugar (a1, D-19): dS removes exactly the fusel carbon
    # converted back to grams of glucose.
    assert schema.get(d, "S") == pytest.approx(-rate * _FUSEL_C / _GLUCOSE_C)
    assert schema.get(d, "S") * _GLUCOSE_C == pytest.approx(-schema.get(d, "fusels") * _FUSEL_C)
    for var in ("X", "E", "N", "CO2", "esters"):
        assert schema.get(d, var) == 0.0


def test_fusel_zero_without_nitrogen(params):
    # Ehrlich needs assimilable nitrogen (amino acids): no YAN, no fusels — the
    # mechanism that front-loads fusel formation into the early, N-replete ferment.
    schema = wine_schema()
    d = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, n=0.0), schema, params)
    assert schema.get(d, "fusels") == 0.0


def test_fusel_rises_with_nitrogen_monotone_branch(params):
    # v1 models the catabolic (monotone-increasing-in-N) branch: more YAN ⇒ more
    # Ehrlich fusels. (The real relationship is non-monotonic — the low-N
    # biosynthetic rise is the documented simplification, kept speculative.)
    schema = wine_schema()
    low = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, n=0.02), schema, params)
    high = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, n=0.2), schema, params)
    assert schema.get(high, "fusels") > schema.get(low, "fusels") > 0.0


def test_fusel_rises_with_temperature(params):
    schema = wine_schema()
    cold = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(warm, "fusels") > schema.get(cold, "fusels") > 0.0


def test_fusel_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    assert FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)[
        schema.slice("fusels")
    ] == pytest.approx(0.0)
    assert FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)[
        schema.slice("fusels")
    ] == pytest.approx(0.0)


# -- integration-level properties ---------------------------------------------


def test_both_run_strict_and_stay_nonnegative(params):
    # Build and run under the strict touches contract: each writes only its own pool
    # and the pools accumulate non-negatively over a full ferment.
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [
            GrowthNitrogenLimited(),
            SugarUptakeToEthanolCO2(),
            EsterSynthesis(),
            FuselAlcoholsEhrlich(),
        ],
        strict=True,
    )
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    # S is now also drawn by the byproduct Processes (a1) — it must stay non-negative
    # too: the proportional draw vanishes as a slot empties, so it cannot overshoot.
    assert_nonnegative(traj, ("esters", "fusels", "S"), atol=1e-12)
    # The aroma pools actually accumulate (the mechanisms are live).
    assert float(traj.series("esters")[-1]) > 0.0
    assert float(traj.series("fusels")[-1]) > 0.0
    # Trace, as expected — orders of magnitude below the g/L ethanol flux.
    assert float(traj.series("esters")[-1]) < 1.0
    assert float(traj.series("fusels")[-1]) < 1.0


def test_byproducts_perturb_only_sugar_and_close_carbon_per_rhs(params):
    # Prime directive #3 under a1 (D-19): the speculative byproduct Processes stay
    # isolable (the core is the ProcessSet built WITHOUT them), but they are no longer
    # byte-for-byte at the derivative level — they route carbon out of S. The exact
    # invariant: enabling them leaves dX/dN/dE/dCO2 byte-for-byte identical (they never
    # touch those), and the only change — on dS — removes EXACTLY the carbon they
    # deposit in esters+fusels. That per-RHS carbon balance is why total_carbon closes.
    schema = wine_schema()
    core = ProcessSet(schema, [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2()])
    with_byp = ProcessSet(
        schema,
        [
            GrowthNitrogenLimited(),
            SugarUptakeToEthanolCO2(),
            EsterSynthesis(),
            FuselAlcoholsEhrlich(),
        ],
    )
    # Evaluate at several representative states spanning the ferment.
    for s, e, n in ((245.0, 0.0, 0.08), (120.0, 60.0, 0.01), (5.0, 120.0, 0.0)):
        y = _wine_y0(schema, x=1.5, s=s, e=e, n=n)
        d_core = core.total_derivatives(0.0, y, params)
        d_byp = with_byp.total_derivatives(0.0, y, params)
        # The non-sugar core variables are byte-for-byte identical (untouched).
        for var in ("X", "E", "N", "CO2"):
            assert d_core[schema.slice(var)] == pytest.approx(d_byp[schema.slice(var)], abs=0.0)
        # dS gains only the byproduct draw, and the carbon balances exactly:
        #   Δ(dS)·c(glucose) + d[esters]·c(ester) + d[fusels]·c(fusel) == 0.
        delta_s = float(d_byp[schema.slice("S")][0] - d_core[schema.slice("S")][0])
        ester_rate = float(d_byp[schema.slice("esters")][0])
        fusel_rate = float(d_byp[schema.slice("fusels")][0])
        assert delta_s <= 0.0  # sugar is drawn down, never created
        carbon_residual = delta_s * _GLUCOSE_C + ester_rate * _ESTER_C + fusel_rate * _FUSEL_C
        assert carbon_residual == pytest.approx(0.0, abs=1e-12)


def test_byproducts_have_only_a_trace_effect_on_the_core(params):
    # Integrated counterpart: routing carbon from sugar (a1) does make the core
    # trajectory move when byproducts are on, but only by the TRACE sugar they consume
    # (~0.2 % of S0). Biomass/nitrogen are uncoupled (never touched) and stay ~identical
    # to solver tolerance; S/E/CO2 drift by well under the ~0.2 g/L the aroma pools draw.
    schema = wine_schema()
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    t_eval = np.linspace(0.0, 400.0, 201)

    core = ProcessSet(schema, [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2()])
    with_byp = ProcessSet(
        schema,
        [
            GrowthNitrogenLimited(),
            SugarUptakeToEthanolCO2(),
            EsterSynthesis(),
            FuselAlcoholsEhrlich(),
        ],
    )
    a = simulate(core, params=params, y0=y0, t_span=(0.0, 400.0), t_eval=t_eval)
    b = simulate(with_byp, params=params, y0=y0, t_span=(0.0, 400.0), t_eval=t_eval)
    assert a.success and b.success
    # Uncoupled — X and N differ only by solver step-adaptivity.
    for var in ("X", "N"):
        np.testing.assert_allclose(a.series(var), b.series(var), rtol=1e-3, atol=1e-3)
    # Coupled but trace — S/E/CO2 stay within the sugar the byproducts divert (< 0.5 g/L
    # on ~245/118/118 g/L pools), so the §2.2 ABV/CO2 bands are unmoved.
    for var in ("S", "E", "CO2"):
        np.testing.assert_allclose(a.series(var), b.series(var), rtol=1e-2, atol=0.5)
    # The ethanol shortfall is real, positive and trace: sugar diverted to aroma is
    # sugar not fermented to ethanol.
    assert 0.0 < float(a.series("E")[-1]) - float(b.series("E")[-1]) < 0.5


def test_total_carbon_closes_with_byproducts_on(params, store):
    # THE invariant option a1 buys (D-19): with esters/fusels carbon-routed from sugar
    # AND weighted in total_carbon, a full ferment conserves carbon to machine
    # precision — the aroma pools are real carbon-accounted state, not an unbooked leak.
    # (Under the interim (b) this assertion could not exist: weighting the pools then
    # would have double-counted the Byp higher-alcohol carbon.)
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [
            GrowthNitrogenLimited(),
            SugarUptakeToEthanolCO2(),
            EsterSynthesis(),
            FuselAlcoholsEhrlich(),
        ],
        strict=True,
    )
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    # Closure is non-trivial only because the pools actually accumulated.
    assert float(traj.series("esters")[-1]) > 0.0
    assert float(traj.series("fusels")[-1]) > 0.0


def test_total_carbon_closes_with_byproducts_on_beer_multislot():
    # The load-bearing multi-slot check (a1, D-19). For beer carbon is the SOLE invariant
    # (total_mass rejects the hydrolysing multi-component sugar), and the proportional
    # draw across three slots with DIFFERENT carbon fractions (glucose 0.40 / maltose
    # 0.42 / maltotriose 0.43) is the non-trivial logic the wine single-slot test cannot
    # exercise. The CO2-ratio benchmark's [0.95, 1.05] window would miss a sub-5% leak
    # here; machine-precision closure is what catches a bad slot distribution.
    store = load_parameters(default_data_dir() / "beer_generic.yaml")
    params = store.resolve()
    schema = beer_schema()
    ps = ProcessSet(
        schema,
        [
            GrowthNitrogenLimited(),
            SugarUptakeToEthanolCO2(),
            EsterSynthesis(),
            FuselAlcoholsEhrlich(),
        ],
        strict=True,
    )
    # The §2.2 ale wort spectrum (glucose/maltose/maltotriose), consumed sequentially.
    y0 = schema.pack(
        {"X": 0.6, "S": [13.2, 54.6, 20.2], "E": 0.0, "N": 0.2, "T": 293.15, "CO2": 0.0}
    )
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    # The proportional draw must vanish as each slot empties (glucose first under
    # sequential uptake), so no sugar slot is driven negative.
    assert_nonnegative(traj, ("S", "esters", "fusels"), atol=1e-9)
    assert float(traj.series("esters")[-1]) > 0.0
    assert float(traj.series("fusels")[-1]) > 0.0


# -- integrated falls-with-temperature property (the load-bearing constraint) --


def _wine_run_to_dryness(celsius: float, duration_days: float):
    """Compile + run the wine medium isothermally; return (reached_dryness, total
    esters+fusels at run end). End-of-run total ≈ total at dryness because byproduct
    production stops once the flux dies with the sugar."""
    sc = Scenario(
        name=f"wine-{celsius}C",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        duration_days=duration_days,
    )
    compiled = compile_scenario(sc, strict=True)
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    assert traj.success, traj.message
    reached_dryness = float(traj.series("S")[-1]) <= 2.0
    total_byproducts = float(traj.series("esters")[-1]) + float(traj.series("fusels")[-1])
    return reached_dryness, total_byproducts


def test_integrated_byproduct_total_falls_with_temperature():
    # THE load-bearing property and the regression guard for the E_a ordering: the
    # snapshot "rises with T" tests above pass for *any* positive E_a, but the
    # run-integrated total only falls with temperature when each byproduct E_a
    # exceeds E_a_uptake (the total scales as exp(-(ΔE_a/R)(1/T - 1/T_ref)); the flux
    # integral to dryness is fixed). If the sourcing step ever drops an E_a toward
    # E_a_uptake, this fails — *before* the formal benchmark is unskipped. Both runs
    # must reach dryness (else the comparison is meaningless), so the colder run gets
    # a generous duration. Mirrors test_lower_temperature_is_slower_but_cleaner.
    cold_dry, cold_total = _wine_run_to_dryness(14.0, 60.0)
    warm_dry, warm_total = _wine_run_to_dryness(25.0, 21.0)
    assert cold_dry and warm_dry, "both temperatures must reach dryness to compare"
    assert 0.0 < cold_total < warm_total, (
        f"colder ferment should be cleaner: cold {cold_total:.4f} vs warm {warm_total:.4f} g/L"
    )


# -- tier propagation ---------------------------------------------------------


def test_ester_tier_capped_by_placeholder_params(store):
    # The ester *form* is plausible, but its placeholder rate params are speculative,
    # so parameter-tier propagation (D-1) caps the esters output at speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterSynthesis()])
    # Structural (form-only) tier is plausible…
    assert ps.tier_of("esters") is Tier.PLAUSIBLE
    # …but folding in the real parameter tiers drops it to speculative.
    assert ps.tier_of("esters", store.tier_map()) is Tier.SPECULATIVE


def test_fusel_form_is_speculative_regardless_of_params():
    # The fusel form itself is speculative (monotone-N simplification), so the fusels
    # output is speculative even before any parameter caps it.
    schema = wine_schema()
    ps = ProcessSet(schema, [FuselAlcoholsEhrlich()])
    assert ps.tier_of("fusels") is Tier.SPECULATIVE
