"""Tests for the Tier-2 aroma-byproduct Processes (decision D-18/D-19).

:class:`EsterSynthesis` and :class:`FuselAlcoholsEhrlich` are additive,
produced-only Processes that fill the ``esters``/``fusels`` pools. These tests pin
the closed-form derivatives and their guards, prove the three properties the beat
requires of every byproduct Process — **produced-only** (touches only its own
pool), **monotone-increasing in temperature**, and **togglable-off ⇒ the validated
core is byte-for-byte unchanged** — and check the nitrogen gate on the Ehrlich
pathway plus tier propagation (the fusel form is speculative; the ester form is
plausible but its placeholder rate params cap its output at speculative).

The run-integrated "cleaner when colder" direction (lower T ⇒ fewer total
esters+fusels) is the *benchmark's* job (``test_lower_temperature_is_slower_but_\
cleaner``); these unit tests cover the per-Process mechanics it rests on.
"""

import numpy as np
import pytest

from fermentation.core.kinetics import (
    EsterSynthesis,
    FuselAlcoholsEhrlich,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
    arrhenius_factor,
)
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_nonnegative


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
    assert set(p.touches) == {"esters"}  # produced-only, touches nothing else
    assert set(p.reads) == {"k_ester", "K_sugar_uptake", "E_a_esters", "T_ref"}


def test_fusel_metadata():
    p = FuselAlcoholsEhrlich()
    assert p.name == "fusel_alcohols_ehrlich"
    # Speculative *form*: the nitrogen dependence is a knowingly-monotone simplification.
    assert p.tier is Tier.SPECULATIVE
    assert set(p.touches) == {"fusels"}
    assert set(p.reads) == {"k_fusel", "K_sugar_uptake", "K_n", "E_a_fusels", "T_ref"}


# -- ester closed form & guards -----------------------------------------------


def test_ester_derivative_matches_closed_form(params):
    schema = wine_schema()
    x, s, t = 2.0, 200.0, 293.15
    y = _wine_y0(schema, x=x, s=s, t=t)
    d = EsterSynthesis().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    f_t = arrhenius_factor(t, params["E_a_esters"], params["T_ref"])
    assert schema.get(d, "esters") == pytest.approx(params["k_ester"] * flux * f_t)
    # Produced-only: nothing else moves — not even the flux it reads from.
    for var in ("X", "S", "E", "N", "CO2", "fusels"):
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
    assert schema.get(d, "fusels") == pytest.approx(params["k_fusel"] * flux * gate * f_t)
    for var in ("X", "S", "E", "N", "CO2", "esters"):
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
    assert_nonnegative(traj, ("esters", "fusels"), atol=1e-12)
    # The aroma pools actually accumulate (the mechanisms are live).
    assert float(traj.series("esters")[-1]) > 0.0
    assert float(traj.series("fusels")[-1]) > 0.0
    # Trace, as expected — orders of magnitude below the g/L ethanol flux.
    assert float(traj.series("esters")[-1]) < 1.0
    assert float(traj.series("fusels")[-1]) < 1.0


def test_togglable_off_leaves_core_derivatives_byte_for_byte(params):
    # Prime directive #3: the speculative byproduct Processes must be isolable. The
    # exact invariant lives at the *derivative* level — they touch only esters/fusels,
    # so the total derivative on every core variable is byte-for-byte identical with
    # and without them. (The integrated trajectory differs only by solver step-size
    # adaptivity, since the error norm sees the new nonzero pool slots — that drift is
    # bounded separately below; the model-level guarantee is the exact one here.)
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
        for var in ("X", "S", "E", "N", "CO2"):
            assert d_core[schema.slice(var)] == pytest.approx(d_byp[schema.slice(var)], abs=0.0)


def test_byproducts_do_not_materially_drift_the_core(params):
    # The integrated core trajectory is unchanged to solver tolerance when the
    # byproduct Processes are added (the ~1e-6 step-adaptivity drift, not a coupling).
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
    for var in ("X", "S", "E", "N", "CO2"):
        np.testing.assert_allclose(a.series(var), b.series(var), rtol=1e-4, atol=1e-4)


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
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h
    )
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
