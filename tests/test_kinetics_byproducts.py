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

:class:`EsterVolatilization` (decision D-20) adds the CO2-stripping sink that moves
liquid ``esters`` into the ``esters_gas`` headspace pool — the physics behind wine's
*falling* liquid ester with warmth (Rollero 2014). Its tests pin the carbon-neutral
liquid→gas transfer (no sugar draw), the flux/temperature dependence, and the clamp.

The run-integrated "cleaner when colder" direction is the *benchmark's* job
(``test_lower_temperature_is_slower_but_cleaner``); these unit tests cover the
per-Process mechanics it rests on, plus the honest per-pool temperature directions
(fusels rise with T; wine liquid esters fall, beer liquid esters rise).
"""

import numpy as np
import pytest

from fermentation.core.chemistry import CARBON_ATOMS, carbon_mass_fraction
from fermentation.core.kinetics import (
    EsterSynthesis,
    EsterVolatilization,
    FuselAlcoholsEhrlich,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
    arrhenius_factor,
)
from fermentation.core.kinetics.carbon_routing import ESTER_SPECS
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
    # Touches ALL THREE ester pools AND S — each ester's carbon is routed from sugar (a1,
    # D-19); never E/CO2 (uptake's ethanol/CO2 production is left untouched). Since D-96 one
    # Process fills three single-molecule pools instead of one lump.
    assert set(p.touches) == {"ethyl_acetate", "isoamyl_acetate", "ethyl_hexanoate", "S"}
    # One INDEPENDENTLY-SOURCED rate per ester (D-96), plus the shared flux/temperature terms.
    # A single k split by a fitted ratio would have been exactly the fabricated-composition
    # constant the split exists to remove, so the plurality here is the point.
    assert set(p.reads) == {
        "k_ethyl_acetate",
        "k_isoamyl_acetate",
        "k_ethyl_hexanoate",
        "K_sugar_uptake",
        "E_a_esters",  # shared: one ATF1 enzyme for the acetates (D-96 documents the hexanoate)
        "T_ref",
    }


def test_fusel_metadata():
    p = FuselAlcoholsEhrlich()
    assert p.name == "fusel_alcohols_ehrlich"
    # Speculative *form*: the nitrogen dependence is a knowingly-monotone simplification.
    assert p.tier is Tier.SPECULATIVE
    assert set(p.touches) == {"fusels", "S"}  # fusel carbon routed from sugar (a1, D-19)
    assert set(p.reads) == {"k_fusel", "K_sugar_uptake", "K_n", "E_a_fusels", "T_ref"}


# -- ester closed form & guards -----------------------------------------------


def test_ester_derivative_matches_closed_form(params):
    """Each ester forms at its OWN rate and draws its OWN carbon — the D-96 ledger payoff.

    The pre-D-96 version of this test checked one lumped pool against one ``k_ester`` and one
    ethyl-acetate carbon fraction. The molecules differ (C4/C7/C8), so a single shared fraction
    would now mis-debit sugar for two of the three; this pins that each ester is weighted as
    itself.
    """
    schema = wine_schema()
    x, s, t = 2.0, 200.0, 293.15
    y = _wine_y0(schema, x=x, s=s, t=t)
    d = EsterSynthesis().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    f_t = arrhenius_factor(t, params["E_a_esters"], params["T_ref"])  # shared shape (D-96)

    total_ester_carbon = 0.0
    for spec in ESTER_SPECS:
        rate = params[spec.k_param] * flux * f_t
        assert schema.get(d, spec.pool) == pytest.approx(rate), spec.pool
        assert rate > 0.0, spec.pool  # a real contribution, not a vacuous 0 == 0
        total_ester_carbon += rate * carbon_mass_fraction(spec.species)

    # Every ester's carbon is routed from sugar (a1, D-19): one wine slot, so dS removes the
    # SUM of the three esters' carbon converted back to grams of glucose. Carbon balances
    # per-RHS — and it only balances because each ester is debited at its own fraction.
    assert schema.get(d, "S") == pytest.approx(-total_ester_carbon / _GLUCOSE_C)
    assert -schema.get(d, "S") * _GLUCOSE_C == pytest.approx(total_ester_carbon)
    # Nothing else moves — ethanol/CO2 production is left to the uptake Process.
    for var in ("X", "E", "N", "CO2", "fusels"):
        assert schema.get(d, var) == 0.0


def test_each_ester_is_carbon_weighted_as_its_own_molecule():
    """The three esters have genuinely different carbon fractions — C4 vs C7 vs C8 (D-96).

    Guards the split's premise. If a future edit collapsed these onto one representative
    species the ledger would silently mis-book two of the three, and every conservation test
    would still pass (the draw and the check would agree — they would just agree on the wrong
    molecule). That symmetry is exactly how the pre-D-96 seam survived so long, so the
    distinctness is asserted directly rather than inferred from closure.
    """
    fractions = {spec.species: carbon_mass_fraction(spec.species) for spec in ESTER_SPECS}
    assert len(set(fractions.values())) == 3, fractions
    assert {CARBON_ATOMS[spec.species] for spec in ESTER_SPECS} == {4, 7, 8}


def test_ester_factor_is_one_at_reference_temperature(params):
    # At T_ref the embedded Arrhenius factor is exactly 1, so the rate is the bare
    # flux term times k_ester (the rate constants are anchored at T_ref).
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = _wine_y0(schema, x=x, s=s, t=params["T_ref"])
    d = EsterSynthesis().derivatives(0.0, y, schema, params)
    flux = x * (s / (params["K_sugar_uptake"] + s))
    assert schema.get(d, "ethyl_acetate") == pytest.approx(params["k_ethyl_acetate"] * flux)


def test_ester_rises_with_temperature(params):
    # The defining directional property: warmer ⇒ more esters per unit flux.
    schema = wine_schema()
    cold = EsterSynthesis().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = EsterSynthesis().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(warm, "ethyl_acetate") > schema.get(cold, "ethyl_acetate") > 0.0


def test_ester_scales_with_fermentative_flux(params):
    # Coupled to the biomass-catalysed sugar flux: more biomass ⇒ proportionally
    # more ester synthesis (the flux is linear in X).
    schema = wine_schema()
    r1 = EsterSynthesis().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = EsterSynthesis().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "ethyl_acetate") == pytest.approx(2.0 * schema.get(r1, "ethyl_acetate"))


def test_ester_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    assert EsterSynthesis().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)[
        schema.slice("ethyl_acetate")
    ] == pytest.approx(0.0)
    assert EsterSynthesis().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)[
        schema.slice("ethyl_acetate")
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
    for var in ("X", "E", "N", "CO2", "ethyl_acetate"):
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


# -- ester volatilization (gas-stripping sink, decision D-20) -----------------


def _wine_y0_with_ester(schema: StateSchema, *, ester: float, **kw) -> FloatArray:
    """A wine state with the liquid ``esters`` pool pre-loaded (the sink needs ester to
    strip; the produced-only pool is 0 at pitch so it must be set explicitly here)."""
    y = _wine_y0(schema, **kw)
    y[schema.slice("ethyl_acetate")] = ester
    return y


def test_volatilization_metadata():
    p = EsterVolatilization()
    assert p.name == "ester_volatilization"
    # Plausible *form* (CO2 stripping is well-understood physics); speculative params cap it.
    assert p.tier is Tier.PLAUSIBLE
    # Pure liquid->gas transfers: touches each liquid ester pool and its OWN headspace twin,
    # never S/E/CO2 (it draws no fresh sugar, unlike synthesis). Since D-96 all three esters
    # are stripped, each into its own twin — a twin per ester is forced, since a pool and its
    # twin must share one molecule's carbon weight for the transfer to stay carbon-neutral.
    assert set(p.touches) == {
        *(spec.pool for spec in ESTER_SPECS),
        *(spec.gas_pool for spec in ESTER_SPECS),
    }
    # Physical Henry model (D-21): gas-flow rides E_a_uptake, partition rides the sourced
    # ethyl-acetate enthalpy dH_ester_volatil — NOT a fudged per-medium E_a_ester_volatil.
    assert set(p.reads) == {
        "k_ester_volatil",
        "K_sugar_uptake",
        "E_a_uptake",
        "dH_ester_volatil",
        "T_ref",
    }


def test_volatilization_derivative_matches_closed_form(params):
    schema = wine_schema()
    x, s, t, est = 2.0, 200.0, 298.15, 0.1  # off T_ref so both Arrhenius factors bite
    y = _wine_y0_with_ester(schema, ester=est, x=x, s=s, t=t)
    d = EsterVolatilization().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    f_gas = arrhenius_factor(t, params["E_a_uptake"], params["T_ref"])  # CO2 gas flow
    f_part = arrhenius_factor(t, params["dH_ester_volatil"], params["T_ref"])  # partition
    rate = params["k_ester_volatil"] * flux * f_gas * f_part * est
    # Liquid loses exactly what the headspace gains — a carbon-neutral transfer.
    assert schema.get(d, "ethyl_acetate") == pytest.approx(-rate)
    assert schema.get(d, "ethyl_acetate_gas") == pytest.approx(rate)
    # Both pools book as ethyl acetate, so the per-RHS carbon residual is exactly zero.
    carbon_residual = (
        schema.get(d, "ethyl_acetate") * _ESTER_C + schema.get(d, "ethyl_acetate_gas") * _ESTER_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=0.0)
    # No fresh sugar drawn, no ethanol/CO2 touched — unlike synthesis (which routes from S).
    for var in ("X", "S", "E", "N", "CO2", "fusels"):
        assert schema.get(d, var) == 0.0


def test_volatilization_zero_without_liquid_ester(params):
    # Nothing in the liquid pool to strip ⇒ no flux into the headspace.
    schema = wine_schema()
    d = EsterVolatilization().derivatives(0.0, _wine_y0(schema), schema, params)  # ester=0
    assert np.array_equal(d, schema.zeros())


def test_volatilization_zero_without_fermentative_flux(params):
    # Stripping rides the CO2-evolution (fermentative-flux) proxy: no flux (no biomass or
    # no sugar) ⇒ no stripping, so liquid esters freeze once the ferment dies (the
    # deliberate "no passive post-ferment evaporation" simplification).
    schema = wine_schema()
    assert np.array_equal(
        EsterVolatilization().derivatives(
            0.0, _wine_y0_with_ester(schema, ester=0.1, x=0.0), schema, params
        ),
        schema.zeros(),
    )
    assert np.array_equal(
        EsterVolatilization().derivatives(
            0.0, _wine_y0_with_ester(schema, ester=0.1, s=0.0), schema, params
        ),
        schema.zeros(),
    )


def test_volatilization_rises_with_temperature(params):
    # The stripping rate per unit liquid ester rises with temperature (volatility climbs):
    # the snapshot property behind the wine inversion (warmer strips more).
    schema = wine_schema()
    cold = EsterVolatilization().derivatives(
        0.0, _wine_y0_with_ester(schema, ester=0.1, t=283.15), schema, params
    )
    warm = EsterVolatilization().derivatives(
        0.0, _wine_y0_with_ester(schema, ester=0.1, t=303.15), schema, params
    )
    assert schema.get(warm, "ethyl_acetate_gas") > schema.get(cold, "ethyl_acetate_gas") > 0.0


def test_volatilization_negative_ester_does_not_strip(params):
    # A solver undershoot (esters < 0) must not flip the clamp into spurious gas creation.
    schema = wine_schema()
    d = EsterVolatilization().derivatives(
        0.0, _wine_y0_with_ester(schema, ester=-1e-6), schema, params
    )
    assert np.array_equal(d, schema.zeros())


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
            EsterVolatilization(),
        ],
        strict=True,
    )
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    # S is now also drawn by the byproduct Processes (a1) — it must stay non-negative
    # too: the proportional draw vanishes as a slot empties, so it cannot overshoot. The
    # liquid esters pool, drawn down by the volatilization sink, must also stay >= 0.
    assert_nonnegative(traj, ("ethyl_acetate", "fusels", "ethyl_acetate_gas", "S"), atol=1e-12)
    # The aroma pools actually accumulate (the mechanisms are live), and the sink fills
    # the headspace pool from the liquid esters it strips (D-20).
    assert float(traj.series("ethyl_acetate")[-1]) > 0.0
    assert float(traj.series("fusels")[-1]) > 0.0
    assert float(traj.series("ethyl_acetate_gas")[-1]) > 0.0
    # Trace, as expected — orders of magnitude below the g/L ethanol flux.
    assert float(traj.series("ethyl_acetate")[-1]) < 1.0
    assert float(traj.series("fusels")[-1]) < 1.0
    assert float(traj.series("ethyl_acetate_gas")[-1]) < 1.0


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
        #   Δ(dS)·c(glucose) + Σ_esters d[ester]·c(that ester) + d[fusels]·c(fusel) == 0.
        # Since D-96 the sum runs over three esters, EACH at its own molecule's fraction
        # (C4/C7/C8) — a single shared fraction would leave a residual here.
        delta_s = float(d_byp[schema.slice("S")][0] - d_core[schema.slice("S")][0])
        fusel_rate = float(d_byp[schema.slice("fusels")][0])
        assert delta_s <= 0.0  # sugar is drawn down, never created
        ester_carbon = sum(
            float(d_byp[schema.slice(spec.pool)][0]) * carbon_mass_fraction(spec.species)
            for spec in ESTER_SPECS
        )
        carbon_residual = delta_s * _GLUCOSE_C + ester_carbon + fusel_rate * _FUSEL_C
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
            EsterVolatilization(),
        ],
        strict=True,
    )
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    # Closure is non-trivial only because the pools actually accumulated — including the
    # volatilized-ester headspace pool, whose carbon must stay counted (D-20) or the
    # liquid→gas transfer would read as destroyed.
    assert float(traj.series("ethyl_acetate")[-1]) > 0.0
    assert float(traj.series("fusels")[-1]) > 0.0
    assert float(traj.series("ethyl_acetate_gas")[-1]) > 0.0


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
            EsterVolatilization(),
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
    # sequential uptake), so no sugar slot is driven negative; the volatilization sink
    # keeps the liquid esters pool non-negative as it strips into the headspace pool.
    assert_nonnegative(traj, ("S", "ethyl_acetate", "fusels", "ethyl_acetate_gas"), atol=1e-9)
    assert float(traj.series("ethyl_acetate")[-1]) > 0.0
    assert float(traj.series("fusels")[-1]) > 0.0
    assert float(traj.series("ethyl_acetate_gas")[-1]) > 0.0


# -- integrated temperature directions (the load-bearing E_a-ordering guards) --


def _wine_run_to_dryness(celsius: float, duration_days: float):
    """Compile + run the wine medium isothermally; return (reached_dryness, pools) where
    ``pools`` is a dict of end-of-run liquid esters / fusels / volatilized esters_gas.
    End-of-run ≈ at-dryness because byproduct production stops once the flux dies."""
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
    pools = {
        "ethyl_acetate": float(traj.series("ethyl_acetate")[-1]),
        "fusels": float(traj.series("fusels")[-1]),
        "ethyl_acetate_gas": float(traj.series("ethyl_acetate_gas")[-1]),
    }
    return reached_dryness, pools


def test_integrated_wine_aroma_temperature_directions():
    # THE load-bearing regression guard for the per-pool E_a ordering (decisions D-19 →
    # D-21). The snapshot "rises with T" tests pass for *any* positive E_a; only the
    # run-integrated pools encode the ordering that matters. The HONEST wine picture, not
    # a combined total that would hide the ester inversion:
    #   * FUSELS rise with T (E_a_fusels > E_a_uptake) — the "cleaner when colder"
    #     direction for the harsh higher alcohols; carries warmer⇒more-aroma for wine.
    #   * LIQUID esters FALL with T — the inversion: the physical Henry's-law stripping
    #     (sensitivity E_a_uptake + dH_ester_volatil ~ 100 kJ/mol) outruns wine's WEAK
    #     synthesis (E_a_esters ~ 15k, Mouret), so the warm ferment's esters end up in the
    #     gas, not the wine (Rollero 2014). The wine/beer split lives in synthesis (D-21).
    #   * VOLATILIZED esters_gas rises with T — the stripped fraction the headspace
    #     pool catches, and the proof the inversion is evaporation, not lost synthesis.
    # If the sourcing ever lifts wine E_a_esters above the stripping sensitivity, the
    # inversion fails here *before* the formal benchmark. Both runs must reach dryness.
    cold_dry, cold = _wine_run_to_dryness(14.0, 90.0)
    warm_dry, warm = _wine_run_to_dryness(25.0, 30.0)
    assert cold_dry and warm_dry, "both temperatures must reach dryness to compare"
    assert 0.0 < cold["fusels"] < warm["fusels"], (
        f"fusels should rise with T (cleaner cold): cold {cold['fusels']:.4f} vs "
        f"warm {warm['fusels']:.4f} g/L"
    )
    assert 0.0 < warm["ethyl_acetate"] < cold["ethyl_acetate"], (
        f"wine LIQUID esters should fall with T (volatilization inversion, D-20): "
        f"cold {cold['ethyl_acetate']:.4f} vs warm {warm['ethyl_acetate']:.4f} g/L"
    )
    assert 0.0 < cold["ethyl_acetate_gas"] < warm["ethyl_acetate_gas"], (
        f"volatilized esters_gas should rise with T (more stripping when warm): "
        f"cold {cold['ethyl_acetate_gas']:.4f} vs "
        f"warm {warm['ethyl_acetate_gas']:.4f} g/L"
    )
    # THE headline D-21 fidelity claim, locked (prime directive #2: enforced, not just
    # honoured): wine TOTAL ester production (liquid + volatilized) is ~FLAT in T. It is
    # the *consequence* of E_a_esters = E_a_uptake (the liquid-falls/gas-rises asserts
    # above can survive a drift that breaks flatness; this cannot). The 2% band is far
    # tighter than the tilt even a ~1 kJ/mol drift of either E_a would cause.
    cold_total = cold["ethyl_acetate"] + cold["ethyl_acetate_gas"]
    warm_total = warm["ethyl_acetate"] + warm["ethyl_acetate_gas"]
    assert cold_total == pytest.approx(warm_total, rel=0.02), (
        f"wine TOTAL ester production should be ~flat in T (D-21 mapping E_a_esters = "
        f"E_a_uptake): cold {cold_total:.4f} vs warm {warm_total:.4f} g/L"
    )


def test_wine_ester_synthesis_e_a_equals_uptake_for_flat_production(store):
    # D-21 mapping as a direct executable guard (the cause, paired with the consequence
    # asserted in ...temperature_directions). Run-integrated wine ester synthesis scales
    # as arrh(E_a_esters)/arrh(E_a_uptake) (the bare-flux integral to dryness is fixed by
    # total sugar), so it is T-INDEPENDENT iff E_a_esters = E_a_uptake — the Arrhenius form
    # of Mouret's flat/weak wine ester production. If a future M1 E_a-band review moves
    # E_a_uptake, this forces a deliberate re-decision rather than a silent tilt of the
    # flat-total fidelity claim.
    assert store.value("E_a_esters") == pytest.approx(store.value("E_a_uptake"))


# -- tier propagation ---------------------------------------------------------


def test_ester_tier_capped_by_placeholder_params(store):
    # The ester *form* is plausible, but its placeholder rate params are speculative,
    # so parameter-tier propagation (D-1) caps the esters output at speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterSynthesis()])
    # Structural (form-only) tier is plausible…
    assert ps.tier_of("ethyl_acetate") is Tier.PLAUSIBLE
    # …but folding in the real parameter tiers drops it to speculative.
    assert ps.tier_of("ethyl_acetate", store.tier_map()) is Tier.SPECULATIVE


def test_fusel_form_is_speculative_regardless_of_params():
    # The fusel form itself is speculative (monotone-N simplification), so the fusels
    # output is speculative even before any parameter caps it.
    schema = wine_schema()
    ps = ProcessSet(schema, [FuselAlcoholsEhrlich()])
    assert ps.tier_of("fusels") is Tier.SPECULATIVE


def test_volatilization_tier_capped_by_placeholder_params(store):
    # Like EsterSynthesis: the gas-stripping *form* is plausible, but its speculative
    # rate params cap the esters_gas output at speculative (parameter-tier propagation).
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterVolatilization()])
    assert ps.tier_of("ethyl_acetate_gas") is Tier.PLAUSIBLE
    assert ps.tier_of("ethyl_acetate_gas", store.tier_map()) is Tier.SPECULATIVE
