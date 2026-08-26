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

import math
import re
from pathlib import Path

import numpy as np
import pytest

from fermentation.core.chemistry import (
    CARBON_ATOMS,
    M_ISOAMYL_ACETATE,
    M_ISOAMYL_OH,
    carbon_mass_fraction,
    co2_yield,
)
from fermentation.core.kinetics import (
    EsterSynthesis,
    EsterVolatilization,
    FuselAlcoholsEhrlich,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
    arrhenius_factor,
)
from fermentation.core.kinetics.carbon_routing import (
    ESTER_SPECS,
    FUSEL_SPECS,
    fermentative_co2_rate,
    fermentative_flux_shape,
    realised_yield_carbon_diversion,
    realised_yield_scale,
)
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.sensory.descriptors import MaxRuleProjector
from fermentation.sensory.oav import sensory_profile
from fermentation.units import cells_per_ml_to_pitch_gpl
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon
from tests.conftest import BEER_COUNTED_PITCH_CELLS_PER_ML

#: Representative species the pools book against (mirrors the Process constants).
_ESTER_C = carbon_mass_fraction("ethyl_acetate")
_GLUCOSE_C = carbon_mass_fraction("glucose")

#: The five single-molecule higher-alcohol pools the lumped `fusels` pool became at D-99.
_FUSEL_POOLS = tuple(spec.pool for spec in FUSEL_SPECS)


def _total_fusels(traj) -> float:
    """End-of-run total higher alcohols [g/L] — what the pre-D-99 lumped pool held.

    Summed from the registry so a sixth alcohol needs no edit here. Tests asserting a
    DIRECTION on total fusels (warmer ⇒ more) use this; tests about the D-97 banana coupling
    or the D-69 hydrolysis must read `isoamyl_alcohol` SPECIFICALLY instead — those are facts
    about one molecule, not about the class.
    """
    return sum(float(traj.series(pool)[-1]) for pool in _FUSEL_POOLS)


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
    isoamyl_alcohol: float = 0.05,
) -> FloatArray:
    # `isoamyl_alcohol` defaults to a realistic mid-ferment ~50 mg/L rather than 0 because since
    # D-97 it is isoamyl acetate's PRECURSOR: at 0 the banana rate is identically 0 and any
    # assertion about it would be vacuously true. It names the MOLECULE since D-99 — the D-97
    # coupling is first-order in 3-methylbutan-1-ol specifically (Fujii 1998 measured ATF1's Km
    # for that alcohol), not in a lump, and not in its C5 isomer `active_amyl_alcohol`.
    return schema.pack(
        {"X": x, "S": [s], "E": e, "N": n, "T": t, "CO2": 0.0, "isoamyl_alcohol": isoamyl_alcohol}
    )


# -- metadata -----------------------------------------------------------------


def test_ester_metadata():
    p = EsterSynthesis()
    assert p.name == "ester_synthesis"
    assert p.tier is Tier.PLAUSIBLE  # warmth-favoured, flux-coupled: standard form
    # Touches ALL THREE ester pools AND S — each ester's carbon is routed from sugar (a1,
    # D-19); never E/CO2 (uptake's ethanol/CO2 production is left untouched). Since D-96 one
    # Process fills three single-molecule pools instead of one lump.
    # Since D-115 it ALSO debits the precursor alcohol (the 5:2-inverse re-route) and writes
    # both label tracers. The D-97-era set excluded ``isoamyl_alcohol`` on the "read, never
    # debited" scope call; D-115 retires that deliberately, so this widening is the beat's
    # content rather than a leak — see the re-route and label tests below for the real pins.
    assert set(p.touches) == {
        "ethyl_acetate",
        "isoamyl_acetate",
        "ethyl_hexanoate",
        "S",
        "isoamyl_alcohol",
        "isoamyl_alcohol_valine",
        "isoamyl_acetate_valine",
    }
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
    # The five single-molecule higher-alcohol pools (D-99) + S: fusel carbon is still routed
    # from sugar (a1, D-19), now at each species' OWN carbon fraction.
    assert set(p.touches) == {*_FUSEL_POOLS, "S"}
    # One k PER SPECIES, independently anchored (D-99) — but still ONE shared E_a and ONE
    # shared N-gate, which is the split's honest limit: the molecules are speciated, the
    # temperature/nitrogen response is not, so the composition is a fixed SPECTRUM.
    assert set(p.reads) == {
        *(spec.k_param for spec in FUSEL_SPECS),
        "K_sugar_uptake", "K_n", "E_a_fusels", "T_ref",
    }  # fmt: skip


# -- ester closed form & guards -----------------------------------------------


def test_ester_derivative_matches_closed_form(params):
    """Each ester forms at its OWN rate and draws its OWN carbon — the D-96 ledger payoff.

    The pre-D-96 version of this test checked one lumped pool against one ``k_ester`` and one
    ethyl-acetate carbon fraction. The molecules differ (C4/C7/C8), so a single shared fraction
    would now mis-debit sugar for two of the three; this pins that each ester is weighted as
    itself.

    Since **D-97** the closed form is no longer uniform across the three: isoamyl acetate
    carries a first-order ``[fusels]`` precursor factor (ATF1 is far from saturated in isoamyl
    alcohol), while the other two do not. That per-ester branch is driven by
    ``EsterSpec.precursor_pool``, so this reconstructs the rate from the registry rather than
    special-casing a name.
    """
    schema = wine_schema()
    x, s, t, isoamyl_alcohol = 2.0, 200.0, 293.15, 0.05
    y = _wine_y0(schema, x=x, s=s, t=t, isoamyl_alcohol=isoamyl_alcohol)
    d = EsterSynthesis().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    f_t = arrhenius_factor(t, params["E_a_esters"], params["T_ref"])  # shared shape (D-96)

    sugar_carbon = 0.0
    for spec in ESTER_SPECS:
        rate = params[spec.k_param] * flux * f_t
        if spec.precursor_pool is not None:
            rate *= schema.get(y, spec.precursor_pool)  # D-97: first-order in the precursor
        assert schema.get(d, spec.pool) == pytest.approx(rate), spec.pool
        assert rate > 0.0, spec.pool  # a real contribution, not a vacuous 0 == 0
        ester_carbon = rate * carbon_mass_fraction(spec.species)
        # Since D-115 the precursor-coupled ester funds only its C2 acetyl group from sugar —
        # the C5 comes off the alcohol pool (the 5:2-inverse re-route). The other two still
        # draw whole, so the closed form is now per-ester on the SOURCE as well as the rate.
        sugar_carbon += ester_carbon * (2.0 / 7.0 if spec.precursor_pool is not None else 1.0)

    # Carbon balances per-RHS — and it only balances because each ester is debited at its own
    # fraction, from its own source.
    assert schema.get(d, "S") == pytest.approx(-sugar_carbon / _GLUCOSE_C)
    assert -schema.get(d, "S") * _GLUCOSE_C == pytest.approx(sugar_carbon)
    # Nothing else moves — ethanol/CO2 production is left to the uptake Process. The isoamyl
    # pool is now excluded from this sweep because D-115 deliberately debits it; the other four
    # higher alcohols must still be untouched, which is what catches a re-route that leaked
    # onto the wrong C5 isomer (`active_amyl_alcohol` — the D-99 trap).
    untouched = tuple(p for p in _FUSEL_POOLS if p != "isoamyl_alcohol")
    assert len(untouched) == len(_FUSEL_POOLS) - 1, "the isoamyl pool must be in _FUSEL_POOLS"
    for var in ("X", "E", "N", "CO2", *untouched):
        assert schema.get(d, var) == 0.0
    assert schema.get(d, "isoamyl_alcohol") < 0.0, "D-115: the acetylation consumes its C5"


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
    shape = flux * gate * f_t  # shared by all five (D-99)

    # Each species gets its OWN k against the one shared shape — so the five rates stand in a
    # fixed ratio set by the k's, which are anchored to five INDEPENDENT measured concentrations
    # (never a ratio-split off one k_fusel; that is the D-96 rule this split obeys).
    for spec in FUSEL_SPECS:
        assert schema.get(d, spec.pool) == pytest.approx(params[spec.k_param] * shape), spec.pool

    # Fusel carbon routed from sugar (a1, D-19): dS removes exactly the summed fusel carbon,
    # each species weighted at ITS OWN molecule's fraction (D-99), converted back to grams of
    # glucose. Before the split this was one rate at isoamyl alcohol's fraction standing in for
    # all five — self-consistent, and therefore invisible to the carbon check.
    fusel_carbon = sum(
        schema.get(d, spec.pool) * carbon_mass_fraction(spec.species) for spec in FUSEL_SPECS
    )
    assert schema.get(d, "S") == pytest.approx(-fusel_carbon / _GLUCOSE_C)
    assert schema.get(d, "S") * _GLUCOSE_C == pytest.approx(-fusel_carbon)
    for var in ("X", "E", "N", "CO2", "ethyl_acetate"):
        assert schema.get(d, var) == 0.0


def test_fusel_zero_without_nitrogen(params):
    # Ehrlich needs assimilable nitrogen (amino acids): no YAN, no fusels — the
    # mechanism that front-loads fusel formation into the early, N-replete ferment.
    schema = wine_schema()
    d = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, n=0.0), schema, params)
    for pool in _FUSEL_POOLS:
        assert schema.get(d, pool) == 0.0, pool


def test_fusel_rises_with_nitrogen_monotone_branch(params):
    # v1 models the catabolic (monotone-increasing-in-N) branch: more YAN ⇒ more
    # Ehrlich fusels. (The real relationship is non-monotonic — the low-N
    # biosynthetic rise is the documented simplification, kept speculative.)
    schema = wine_schema()
    low = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, n=0.02), schema, params)
    high = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, n=0.2), schema, params)
    # All five rise together — they share the one N-gate (D-99). No ratio among them can move
    # with nitrogen, which is precisely the fixed-spectrum limitation the split does NOT retire.
    for pool in _FUSEL_POOLS:
        assert schema.get(high, pool) > schema.get(low, pool) > 0.0, pool


def test_fusel_rises_with_temperature(params):
    schema = wine_schema()
    cold = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = FuselAlcoholsEhrlich().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    # Again all five together: one shared E_a_fusels (D-99). A per-species E_a is the deferred
    # refinement that would make the SPECTRUM temperature-dependent — unsourced, so not built.
    for pool in _FUSEL_POOLS:
        assert schema.get(warm, pool) > schema.get(cold, pool) > 0.0, pool


def test_fusel_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    for starved in (_wine_y0(schema, x=0.0), _wine_y0(schema, s=0.0)):
        d = FuselAlcoholsEhrlich().derivatives(0.0, starved, schema, params)
        for pool in _FUSEL_POOLS:
            assert d[schema.slice(pool)] == pytest.approx(0.0), pool


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
    # Since D-115 the banana ester's label tracer joins them: stripping is non-fractionating, so
    # the tracer must fall with the pool or the enrichment would inflate for no physical reason.
    assert set(p.touches) == {
        *(spec.pool for spec in ESTER_SPECS),
        *(spec.gas_pool for spec in ESTER_SPECS),
        "isoamyl_acetate_valine",
    }
    # Physical Henry model (D-21): gas-flow rides E_a_uptake, partition rides the sourced
    # ethyl-acetate enthalpy dH_ester_volatil — NOT a fudged per-medium E_a_ester_volatil.
    #
    # D-227 GREW THIS SET, and the growth is the change rather than bookkeeping. The driver
    # stopped being a single-Monod stand-in for the CO2 stream and became the fermentative CO2
    # rate itself, so this Process now reads the whole uptake rate law — including the speed
    # knob `q_sugar_max`, which used to sit FOLDED INSIDE `k_ester_volatil` (that parameter's
    # own unit comment said so) where nothing could see it move. `reads` has two masters, tier
    # propagation AND sampler scope (D-160), so every one of them is declared.
    assert set(p.reads) == {
        "k_ester_volatil",
        "q_sugar_max",
        "K_sugar_uptake",
        "K_repression",
        "Y_glycerol_sugar",
        "Y_byproduct_sugar",
        "E_a_uptake",
        "dH_ester_volatil",
        "T_ref",
    }
    # The four that D-227 added are exactly the uptake Process's own parameter list. Asserted
    # as a relation rather than a second literal: if uptake's rate law gains a parameter and
    # this sink does not declare it, the sink's tier and its sampler scope both go quietly
    # wrong, and a hand-copied list is what would let that happen (D-225's rule).
    assert set(SugarUptakeToEthanolCO2().reads) <= set(p.reads)


def test_volatilization_derivative_matches_closed_form(params):
    schema = wine_schema()
    x, s, t, est = 2.0, 200.0, 298.15, 0.1  # off T_ref so both Arrhenius factors bite
    y = _wine_y0_with_ester(schema, ester=est, x=x, s=s, t=t)
    d = EsterVolatilization().derivatives(0.0, y, schema, params)

    # D-227: the driver is the CO2 the ferment actually evolves, so the closed form carries
    # `q_sugar_max` and `co2_yield` EXPLICITLY where they used to be folded into
    # `k_ester_volatil`. Wine is one sugar slot, so there is no repression term to write.
    monod = x * (s / (params["K_sugar_uptake"] + s))
    diverted_c = realised_yield_carbon_diversion(params)
    co2_rate = (
        params["q_sugar_max"]
        * monod
        * co2_yield("glucose")
        * realised_yield_scale("glucose", diverted_c)
    )
    f_gas = arrhenius_factor(t, params["E_a_uptake"], params["T_ref"])  # CO2 gas flow
    f_part = arrhenius_factor(t, params["dH_ester_volatil"], params["T_ref"])  # partition
    rate = params["k_ester_volatil"] * co2_rate * f_gas * f_part * est
    # The retired driver, kept as an explicit NON-assertion: the two differ by the factor the
    # re-anchoring absorbed, so a revert of the Process cannot pass this test by coincidence.
    assert rate != pytest.approx(params["k_ester_volatil"] * monod * f_gas * f_part * est, rel=1e-3)
    # Liquid loses exactly what the headspace gains — a carbon-neutral transfer.
    assert schema.get(d, "ethyl_acetate") == pytest.approx(-rate)
    assert schema.get(d, "ethyl_acetate_gas") == pytest.approx(rate)
    # Both pools book as ethyl acetate, so the per-RHS carbon residual is exactly zero.
    carbon_residual = (
        schema.get(d, "ethyl_acetate") * _ESTER_C + schema.get(d, "ethyl_acetate_gas") * _ESTER_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=0.0)
    # No fresh sugar drawn, no ethanol/CO2 touched — unlike synthesis (which routes from S).
    for var in ("X", "S", "E", "N", "CO2", *_FUSEL_POOLS):
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
    assert_nonnegative(traj, ("ethyl_acetate", *_FUSEL_POOLS, "ethyl_acetate_gas", "S"), atol=1e-12)
    # The aroma pools actually accumulate (the mechanisms are live), and the sink fills
    # the headspace pool from the liquid esters it strips (D-20).
    assert float(traj.series("ethyl_acetate")[-1]) > 0.0
    assert _total_fusels(traj) > 0.0
    assert float(traj.series("ethyl_acetate_gas")[-1]) > 0.0
    # Trace, as expected — orders of magnitude below the g/L ethanol flux.
    assert float(traj.series("ethyl_acetate")[-1]) < 1.0
    assert _total_fusels(traj) < 1.0
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
        #   Δ(dS)·c(glucose) + Σ_esters d[ester]·c(it) + Σ_fusels d[alcohol]·c(it) == 0.
        # Since D-96 the ester sum runs over three molecules and since D-99 the fusel sum runs
        # over five, EACH at its own fraction — a single shared fraction would leave a residual
        # here. NB it would NOT have before D-99: the lumped pool was drawn from S at the very
        # same stand-in fraction it was weighted by, so the error cancelled exactly and this
        # assertion passed on a wrong weight. Splitting the pool is what made it mean something.
        delta_s = float(d_byp[schema.slice("S")][0] - d_core[schema.slice("S")][0])
        assert delta_s <= 0.0  # sugar is drawn down, never created
        ester_carbon = sum(
            float(d_byp[schema.slice(spec.pool)][0]) * carbon_mass_fraction(spec.species)
            for spec in ESTER_SPECS
        )
        fusel_carbon = sum(
            float(d_byp[schema.slice(spec.pool)][0]) * carbon_mass_fraction(spec.species)
            for spec in FUSEL_SPECS
        )
        carbon_residual = delta_s * _GLUCOSE_C + ester_carbon + fusel_carbon
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
    assert _total_fusels(traj) > 0.0
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
    assert_nonnegative(traj, ("S", "ethyl_acetate", *_FUSEL_POOLS, "ethyl_acetate_gas"), atol=1e-9)
    assert float(traj.series("ethyl_acetate")[-1]) > 0.0
    assert _total_fusels(traj) > 0.0
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
        # Summed across the five single-molecule pools (D-99) — the same total quantity the
        # pre-D-99 lump held, so the "fusels rise with T" direction is tested on what it always
        # was. They share one E_a_fusels, so all five rise together and the sum is faithful.
        "fusels": _total_fusels(traj),
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
    #     synthesis (E_a_esters = 55,100 = E_a_uptake exactly, the condition for T-flat
    #     integrated production — Mouret-flat, D-21; NOT the "~15k" this comment carried
    #     before D-168), so the warm ferment's esters end up in the
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
    for pool in _FUSEL_POOLS:
        assert ps.tier_of(pool) is Tier.SPECULATIVE, pool


def test_volatilization_tier_capped_by_placeholder_params(store):
    # Like EsterSynthesis: the gas-stripping *form* is plausible, but its speculative
    # rate params cap the esters_gas output at speculative (parameter-tier propagation).
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterVolatilization()])
    assert ps.tier_of("ethyl_acetate_gas") is Tier.PLAUSIBLE
    assert ps.tier_of("ethyl_acetate_gas", store.tier_map()) is Tier.SPECULATIVE


# -- D-97: the ATF1 precursor coupling ----------------------------------------
#
# The banana ester's rate is FIRST-ORDER in its precursor alcohol (the `fusels` pool),
# because ATF1's measured Km for isoamyl alcohol (~29.8 mM; Fujii 1998, AEM 64:4076-4078)
# sits ~30-60x ABOVE the pool the sim actually carries (~0.5-1 mM) — the [S] << Km limit.
# Ethyl acetate gets no such term: its precursor is ethanol (~2 M), which saturates the
# same enzyme. Same enzyme, same rate law, opposite limits.


def _wine_isoamyl_run(yan_mgl: float):
    """A real solver run to dryness at a given YAN; returns the finished pools."""
    sc = Scenario(
        name=f"d97-yan-{yan_mgl}",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )
    compiled = compile_scenario(sc, strict=True)
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    assert traj.success, traj.message
    assert float(traj.series("S")[-1]) <= 2.0, "run must reach dryness to compare finished pools"
    return {
        "isoamyl_acetate": float(traj.series("isoamyl_acetate")[-1]),
        "ethyl_acetate": float(traj.series("ethyl_acetate")[-1]),
        # The D-97 PRECURSOR — isoamyl alcohol specifically since D-99, not the class. The
        # coupling is first-order in 3-methylbutan-1-ol because that is the alcohol Fujii 1998
        # measured ATF1's Km for; summing the five here (or reading the C5 isomer) would be a
        # different, unmeasured quantity.
        "isoamyl_alcohol": float(traj.series("isoamyl_alcohol")[-1]),
    }


def test_low_yan_makes_less_banana_than_high_yan_end_to_end():
    """THE D-97 OUTCOME — the observable the pre-D-97 model got WRONG.

    This is the test the beat exists for, and the one that pins the *result* rather than the
    mechanism (the D-96 done-call lesson: every other test here pins how it works, so a
    regression that silently unwired the coupling would slip past them all). Before D-97 the
    banana ester was **YAN-blind** — flat at 0.759/0.758/0.756 mg/L across YAN 40/80/250,
    because every ester shared one plain flux shape and nothing downstream of nitrogen
    reached it. That is physically wrong: ATF1 acetylates isoamyl alcohol, the Ehrlich
    pathway that builds isoamyl alcohol is nitrogen-gated, so a nitrogen-starved must has
    less precursor to acetylate and MUST make less banana.

    Deliberately asserted as a RATIO against the fusel swing, not as absolute values: the
    point is that the ester *tracks its precursor*, which is the derived consequence, whereas
    an absolute band would also pass for a k that had simply been retuned.
    """
    low = _wine_isoamyl_run(40.0)
    high = _wine_isoamyl_run(250.0)

    # The precursor pool must genuinely differ, or the test proves nothing about coupling.
    fusel_ratio = high["isoamyl_alcohol"] / low["isoamyl_alcohol"]
    assert fusel_ratio > 2.0, (
        f"precondition: YAN 40 -> 250 must move the fusel pool substantially for this test "
        f"to have teeth (got {fusel_ratio:.2f}x: {low['isoamyl_alcohol']:.4f} -> "
        f"{high['isoamyl_alcohol']:.4f} g/L)"
    )

    # THE CLAIM: the banana tracks its precursor. Pre-D-97 this ratio was 1.00 (0.756/0.759).
    ester_ratio = high["isoamyl_acetate"] / low["isoamyl_acetate"]
    assert ester_ratio > 2.0, (
        f"low-YAN must make substantially LESS isoamyl acetate than high-YAN (D-97): "
        f"{low['isoamyl_acetate'] * 1000:.3f} vs {high['isoamyl_acetate'] * 1000:.3f} mg/L "
        f"= {ester_ratio:.2f}x. Pre-D-97 this was 1.00x - the ester was YAN-blind."
    )
    # And it tracks it CLOSELY - first-order in the pool, so the ester ratio should land near
    # the fusel ratio rather than merely somewhere above 1. This is what distinguishes the
    # first-order coupling from any weaker or saturating dependence.
    assert ester_ratio == pytest.approx(fusel_ratio, rel=0.25), (
        f"first-order in the precursor => the ester swing ({ester_ratio:.2f}x) should track "
        f"the fusel swing ({fusel_ratio:.2f}x) closely"
    )


def test_ethyl_acetate_stays_yan_blind_while_the_banana_responds():
    """The ASYMMETRY, pinned: ethyl acetate must NOT inherit the precursor coupling.

    Both acetates are ATF1 products, so the naive "couple the acetate esters to their
    precursor" would gate BOTH - and this test is what fails if someone does. Ethanol runs
    ~2 M against the same mM-scale Km => ATF1 is saturated in it => zeroth-order => ethyl
    acetate cannot respond to its precursor's supply, and (having no nitrogen dependence
    anywhere upstream) stays YAN-blind. The asymmetry between two esters of the SAME enzyme
    is derived from the precursors' concentrations, not an exemption granted by hand.
    """
    low = _wine_isoamyl_run(40.0)
    high = _wine_isoamyl_run(250.0)

    # Ethyl acetate moves only via the small flux/biomass differences YAN causes - nowhere
    # near the banana's response.
    ea_ratio = high["ethyl_acetate"] / low["ethyl_acetate"]
    ia_ratio = high["isoamyl_acetate"] / low["isoamyl_acetate"]
    assert ea_ratio == pytest.approx(1.0, abs=0.15), (
        f"ethyl acetate should stay ~YAN-blind (its ethanol precursor saturates ATF1): "
        f"{ea_ratio:.3f}x"
    )
    assert ia_ratio > 2.0 * ea_ratio, (
        f"the banana must respond far more strongly to YAN than ethyl acetate does "
        f"({ia_ratio:.2f}x vs {ea_ratio:.2f}x) - that asymmetry IS the D-97 claim"
    )


def test_fusel_pool_stays_far_below_atf1_km_so_the_first_order_form_holds():
    """The FORM's precondition, checked on a real run - not left as an unchecked assumption.

    First-order-in-precursor is only correct in the ``[S] << Km`` limit. That is a claim about
    the concentrations the sim actually reaches, so it is checkable *here* rather than merely
    asserted in provenance: if the fusel pool ever climbed toward ATF1's Km, the linear form
    would silently become wrong (over-predicting the banana without saturating) and the
    saturable form with the measured Km would be required instead. Guarding at 10% of Km
    keeps the linear approximation within ~10%; the model in fact runs an order of magnitude
    below that.
    """
    # Fujii 1998, AEM 64:4076-4078 (citing Yoshioka & Hashimoto 1981): Km ~29.8 mM.
    km_mm = 29.8
    isoamyl_alcohol_molar_mass = 88.15  # g/mol, C5H12O (3-methylbutan-1-ol, CAS 123-51-3)

    # Reads the isoamyl_alcohol POOL since D-99, which is what makes this check exact rather
    # than indicative: Fujii's Km is for 3-methylbutan-1-ol, and before the split this compared
    # it against a LUMP of five molecules (~52% of which was actually isoamyl alcohol). The
    # comparison is now like-for-like — and the margin TIGHTENED, because honest anchoring
    # raised this pool ~2x. It still clears 10% of Km, so the [S] << Km limit still holds.
    worst = max(_wine_isoamyl_run(yan)["isoamyl_alcohol"] for yan in (40.0, 250.0))
    pool_mm = worst / isoamyl_alcohol_molar_mass * 1000.0
    assert pool_mm < 0.1 * km_mm, (
        f"the fusel pool ({pool_mm:.3f} mM) must stay well below ATF1's Km ({km_mm} mM) for "
        f"the first-order form to hold; it currently runs at {pool_mm / km_mm:.1%} of Km"
    )


def test_only_the_banana_declares_a_precursor_and_it_is_a_real_pool():
    """The registry contract: exactly one ester is precursor-coupled, and it names a real slot.

    Pins the D-96 promise that ``ESTER_SPECS`` drives every layer - a phantom precursor pool
    would otherwise fail only at runtime, deep inside the solver.
    """
    coupled = [spec for spec in ESTER_SPECS if spec.precursor_pool is not None]
    assert [spec.pool for spec in coupled] == ["isoamyl_acetate"], (
        "only isoamyl acetate is precursor-limited (D-97): ethyl acetate's ethanol precursor "
        "saturates ATF1, and ethyl hexanoate's hexanoyl-CoA precursor is not modelled at all"
    )
    # Names the MOLECULE since D-99, not the lump: Fujii 1998 measured ATF1's Km for
    # 3-methylbutan-1-ol specifically, and `active_amyl_alcohol` is a different substrate.
    assert coupled[0].precursor_pool == "isoamyl_alcohol"
    # No phantom: the named precursor must be a real state slot in BOTH media.
    for schema in (wine_schema(), beer_schema()):
        assert coupled[0].precursor_pool in schema.names


def test_the_acetylation_sources_its_c5_off_the_alcohol_and_only_its_c2_off_sugar(params):
    """THE 5:2-INVERSE RE-ROUTE (decision D-115) - what D-97 deferred, built.

    D-97 read the alcohol pool without debiting it and drew the whole ester from ``S``, on the
    correct reasoning that the re-route is **mass-negligible** (~0.5 mg/L of ester against an
    ~86 mg/L alcohol pool). D-114 established that this bounds the wrong quantity: the deferral
    was 100 % of the ester's valine enrichment, because an ester made wholly of sugar carbon is
    structurally incapable of carrying an amino-acid label. So the acetylation now takes the C5
    skeleton off ``isoamyl_alcohol`` and only the C2 acetyl group off ``S``.

    The predecessor of this test asserted the exact opposite (``pool not in proc.touches``,
    "reading a pool must not silently become debiting it"). That assertion is **retired by
    design, not broken** - it pinned a boundary D-97 chose and D-115 deliberately moves.

    **Asserted on the carbon split, not on ``touches``.** A ``touches`` assertion is a tautology
    about the declaration; the load-bearing facts are that the alcohol debit is *mole-for-mole*
    with the ester and that the sugar draw fell to exactly 2/7 of the ester's carbon.
    """
    schema = wine_schema()
    proc = EsterSynthesis()
    y = _wine_y0(schema, x=1.5, s=180.0, e=40.0, n=0.05)
    y[schema.slice("isoamyl_alcohol")] = 0.05
    d = proc.derivatives(0.0, y, schema, params)

    banana = float(d[schema.slice("isoamyl_acetate")][0])
    assert banana > 0.0, "vacuous: no ester is being synthesised at this state"

    # The alcohol is now DEBITED - and mole for mole, one alcohol molecule per ester molecule.
    alcohol = -float(d[schema.slice("isoamyl_alcohol")][0])
    assert alcohol > 0.0, "D-115: the acetylation must consume its precursor alcohol"
    assert alcohol / M_ISOAMYL_OH == pytest.approx(banana / M_ISOAMYL_ACETATE, rel=1e-12), (
        "one alcohol molecule per ester molecule - the C5 skeleton transfers as a unit"
    )

    # ...and the split is exactly 5:2 by carbon, the inverse of the D-69 hydrolysis.
    ester_carbon = banana * carbon_mass_fraction("isoamyl_acetate")
    alcohol_carbon = alcohol * carbon_mass_fraction("isoamyl_alcohol")
    assert alcohol_carbon == pytest.approx(ester_carbon * 5.0 / 7.0, rel=1e-12)

    # The OTHER two esters are untouched by this: they have no precursor pool, so they still
    # draw whole from S. Without this the test would pass on a build that re-routed everything.
    ungated_carbon = sum(
        float(d[schema.slice(spec.pool)][0]) * carbon_mass_fraction(spec.species)
        for spec in ESTER_SPECS
        if spec.precursor_pool is None
    )
    sugar_carbon = -float(d[schema.slice("S")][0]) * _GLUCOSE_C
    assert sugar_carbon == pytest.approx(ungated_carbon + ester_carbon * 2.0 / 7.0, rel=1e-12), (
        "sugar must now fund the ungated esters in full plus ONLY the C2 acetyl group of the "
        "banana ester"
    )


def test_the_alcohol_debit_cannot_drive_its_pool_negative(params):
    """The re-route is SELF-LIMITING, structurally - no clamp, and none needed (D-115).

    The debit is proportional to a rate that is itself first-order in the same pool (D-97), so
    it is first-order decay: it vanishes exactly as fast as the pool does. This matters because
    a clamp here would be actively harmful - it would break the ester's own mass balance rather
    than prevent anything that can happen - so the safety has to come from the algebra instead.

    Asserted as a *ratio* rather than at one pool value, because "the derivative is small at a
    small pool" is true of any decaying term; what makes it safe is that the ratio is CONSTANT.
    """
    schema = wine_schema()
    proc = EsterSynthesis()

    def debit_per_unit_pool(pool: float) -> float:
        y = _wine_y0(schema, x=1.5, s=180.0, e=40.0, n=0.05)
        y[schema.slice("isoamyl_alcohol")] = pool
        d = proc.derivatives(0.0, y, schema, params)
        return -float(d[schema.slice("isoamyl_alcohol")][0]) / pool

    reference = debit_per_unit_pool(0.05)
    assert reference > 0.0, "vacuous: nothing is being debited at all"
    for pool in (1e-9, 1e-6, 1e-3, 0.05, 0.5):
        assert debit_per_unit_pool(pool) == pytest.approx(reference, rel=1e-12), (
            "the debit must stay strictly proportional to the pool at every scale - that "
            "proportionality IS the guarantee the pool cannot go negative"
        )

    # And at exactly zero the whole branch is inert rather than dividing by zero.
    y = _wine_y0(schema, x=1.5, s=180.0, e=40.0, n=0.05)
    y[schema.slice("isoamyl_alcohol")] = 0.0
    d = proc.derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("isoamyl_alcohol")][0]) == 0.0
    assert float(d[schema.slice("isoamyl_acetate")][0]) == 0.0


def test_the_label_rides_across_the_acetylation_at_the_alcohol_pools_fraction(params):
    """D-115 - the ester now carries VALINE label, at its parent alcohol's current fraction.

    This replaces ``test_the_ester_carries_no_amino_acid_label_because_its_carbon_is_wholly_
    sugar_sourced``, whose own docstring predicted this: *"if a future beat builds the deferred
    D-69 5:2-inverse re-route ... this test FAILS, which is the finding's premise being
    correctly invalidated rather than silently going green on a stale claim."* It was, and this
    is the successor.

    **The mechanism.** Rollero 2017 defines enrichment as a **molecule** fraction (D-111 Finding
    3), and the C5 skeleton transfers through the acetylation as a unit - so an ester molecule
    is labelled exactly when its parent alcohol molecule was, and the fraction passes across the
    reaction unchanged. That identity is what this test pins, at three different pool fractions
    so it cannot pass on a build that hard-codes one.
    """
    schema = wine_schema()
    proc = EsterSynthesis()

    for fraction in (0.0, 0.25, 1.0):
        y = _wine_y0(schema, x=1.5, s=180.0, e=40.0, n=0.05)
        y[schema.slice("isoamyl_alcohol")] = 0.05
        y[schema.slice("isoamyl_alcohol_valine")] = 0.05 * fraction
        d = proc.derivatives(0.0, y, schema, params)

        banana = float(d[schema.slice("isoamyl_acetate")][0])
        assert banana > 0.0, "vacuous: no ester is being synthesised at this state"

        # The ester is credited with label at EXACTLY the alcohol pool's fraction.
        labelled = float(d[schema.slice("isoamyl_acetate_valine")][0])
        assert labelled == pytest.approx(banana * fraction, rel=1e-12), (
            "an ester molecule is labelled iff its parent alcohol was - the fraction must "
            "cross the acetylation unchanged"
        )
        # ...and the alcohol tracer is debited at the same fraction, so the ALCOHOL pool's own
        # enrichment is left exactly where it was. A re-route that moved mass without moving
        # label in step would silently enrich (or dilute) the source pool.
        alcohol = -float(d[schema.slice("isoamyl_alcohol")][0])
        alcohol_labelled = -float(d[schema.slice("isoamyl_alcohol_valine")][0])
        assert alcohol_labelled == pytest.approx(alcohol * fraction, rel=1e-12), (
            "the draw is NON-FRACTIONATING: it must remove labelled and unlabelled alcohol in "
            "the proportion the pool holds them"
        )


def test_stripping_is_non_fractionating_so_it_cannot_inflate_the_enrichment(params):
    """The classic tracer bug, pinned shut (D-115): mass leaves, label must leave with it.

    :class:`EsterVolatilization` strips liquid ester to the headspace. If it debited the pool
    but not the tracer, the enrichment ``tracer/bulk`` would climb toward 100 % purely because
    the denominator shrank - a rising "measurement" produced entirely by a bookkeeping omission,
    and one no conservation test could catch (a tracer slot carries no carbon weight, so the
    ledger is blind to it by construction - the D-89/D-90 family).

    Asserted as *the fraction is invariant*, which is the physical statement, rather than as
    "the tracer derivative is negative", which would pass at any wrong magnitude.
    """
    schema = wine_schema()
    proc = EsterVolatilization()
    y = _wine_y0(schema, x=1.5, s=180.0, e=40.0, n=0.05)
    y[schema.slice("isoamyl_acetate")] = 1.0e-3
    y[schema.slice("isoamyl_acetate_valine")] = 0.25e-3
    d = proc.derivatives(0.0, y, schema, params)

    stripped = -float(d[schema.slice("isoamyl_acetate")][0])
    assert stripped > 0.0, "vacuous: nothing is being stripped"
    stripped_label = -float(d[schema.slice("isoamyl_acetate_valine")][0])
    assert stripped_label == pytest.approx(stripped * 0.25, rel=1e-12), (
        "stripping must remove labelled and unlabelled molecules in the pool's own proportion"
    )


def test_banana_rate_is_first_order_in_the_fusel_pool(params):
    """The MECHANISM at the derivative level: double the precursor, double the banana rate -
    and leave the other two esters exactly where they were.

    The pure-function counterpart to the end-to-end outcome test: it isolates the coupling
    from every confound a full run carries (biomass, stripping, dryness).
    """
    schema = wine_schema()
    proc = EsterSynthesis()

    def rates(fusels: float):
        y = _wine_y0(schema, isoamyl_alcohol=fusels)
        d = proc.derivatives(0.0, y, schema, params)
        return {spec.pool: float(d[schema.slice(spec.pool)][0]) for spec in ESTER_SPECS}

    single = rates(0.05)
    double = rates(0.10)
    assert double["isoamyl_acetate"] == pytest.approx(2.0 * single["isoamyl_acetate"]), (
        "the banana rate is FIRST-ORDER in the fusel pool (D-97)"
    )
    # The ungated esters must be untouched by the precursor pool.
    assert double["ethyl_acetate"] == pytest.approx(single["ethyl_acetate"])
    assert double["ethyl_hexanoate"] == pytest.approx(single["ethyl_hexanoate"])
    # And with no precursor at all there is no banana - but the others still form.
    empty = rates(0.0)
    assert empty["isoamyl_acetate"] == 0.0, "no precursor alcohol => no acetylation => no banana"
    assert empty["ethyl_acetate"] > 0.0, "ethyl acetate does not depend on the fusel pool"


# -- D-224: the finished beer's aroma LEVELS, which nothing in the suite pinned ------------
#
# Seven of beer's rate constants are DEFINED by a landing level -- their `conditions:` fields
# say "k set to land finished X at ...". Nothing tested that they still do. Between D-211 and
# D-223 three beats moved beer's growth rate and its sugar-uptake rate, and all seven levels
# moved with them: the five Ehrlich higher alcohols by x2.87 and then back to x1.68, the two
# esters by x0.79. The whole 1842-test suite went red ONCE in the process, and only because
# D-223 had happened to pin a packaged ethyl-acetate number four commits earlier for an
# unrelated reason. These guards are the durable half of D-224 -- the re-anchoring is the
# cheap half.

#: The landing level each constant's own `conditions:` field names, mg/L. Four are
#: Wang/Frank/Steinhaus 2024 Table 1 beer means (n = 64-92 studies each); propan-1-ol is an
#: author estimate (that survey omits it); the two esters are mid-band picks inside their own
#: measured ale ranges (Meilgaard 1975: 10-30 and 0.5-3 mg/L).
_BEER_AROMA_TARGETS_MGL = {
    "propanol": 10.0,
    "isobutanol": 9.6,
    "active_amyl_alcohol": 10.3,
    "isoamyl_alcohol": 30.0,
    "2_phenylethanol": 25.7,
    "ethyl_acetate": 20.0,
    "isoamyl_acetate": 2.2,
}
#: How each target is SPELLED in its own parameter's `conditions:` field, so the numbers above
#: can be checked against the sentence that specifies them rather than standing as a second,
#: silently-divergent copy (the D-158 discipline: recompute from the seam, never re-transcribe).
_BEER_AROMA_TARGET_LITERALS = {
    "propanol": "10 mg/L",
    "isobutanol": "9.6 mg/L",
    "active_amyl_alcohol": "10.3 mg/L",
    "isoamyl_alcohol": "30.0 mg/L",
    "2_phenylethanol": "25.7 mg/L",
    "ethyl_acetate": "20.0 mg/L",
    "isoamyl_acetate": "2.2 mg/L",
}
#: The frame the constants are calibrated in, stated rather than implied. It is load-bearing:
#: `E_a_esters` is 200 kJ/mol, so the same run at 15 C lands ethyl acetate at 6.10 mg/L, below
#: its own 10 mg/L floor. 20 C is a typical ale ferment; the 15 C in the Sec 2.2 DURATION
#: criterion is Foster 2022's cool trial (D-221), a different frame and not a competing one.
_BEER_CALIBRATION_DAYS = 21.0
_BEER_CALIBRATION_CELSIUS = 20.0

#: The frame's INOCULUM, and it is a COUNT (D-228). It carried a flat 1.0 g/L from D-99 until
#: D-228 — 2.5x a counted ale pitch, and the same back-computed residual D-219 retired and
#: D-222 corrected in `TYRELL_SCENARIO`; this frame was simply never revisited. It is what made
#: "growth stops dead at day 0.92 with 81 % of the sugar unfermented" (D-226 §8) TRUE HERE and
#: nowhere the archive has counts: at 1.0 g/L the nitrogen ceiling is reached in 0.79 d and the
#: five higher alcohols are at 99.99 % of their finished level by day 1, where at this counted
#: pitch they are at 67 %. The eight calibrated LEVELS are invariant to the correction (6e-5 %,
#: measured across both registries below), which is why it is free — and why it owes a guard
#: that can see it, since no level guard can.
_BEER_CALIBRATION_PITCH_GPL = cells_per_ml_to_pitch_gpl(BEER_COUNTED_PITCH_CELLS_PER_ML)

#: What this frame carried before D-228, named rather than left in the history: the guards below
#: run BOTH arms, so a value that is only mentioned in prose could not be one of them.
_BEER_CALIBRATION_RETIRED_PITCH_GPL = 1.0

#: The five Ehrlich pools, spelled out rather than derived, because the claim under test is
#: exactly that these five move together and the two esters do not.
_BEER_EHRLICH_POOLS = (
    "propanol",
    "isobutanol",
    "active_amyl_alcohol",
    "isoamyl_alcohol",
    "2_phenylethanol",
)


def _beer_calibration_scenario(celsius: float = _BEER_CALIBRATION_CELSIUS) -> Scenario:
    return Scenario(
        name="beer-aroma-calibration",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 200.0,
            "pitch_gpl": _BEER_CALIBRATION_PITCH_GPL,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        duration_days=_BEER_CALIBRATION_DAYS,
    )


def _beer_aroma_levels_mgl(data_dir=None, celsius: float = _BEER_CALIBRATION_CELSIUS):
    """The seven finished levels, mg/L, from one compiled run at the calibration frame."""
    scenario = _beer_calibration_scenario(celsius)
    compiled = (
        compile_scenario(scenario, data_dir=data_dir)
        if data_dir is not None
        else compile_scenario(scenario)
    )
    traj = compiled.run(t_eval=np.array([0.0, _BEER_CALIBRATION_DAYS * 24.0]))
    return {pool: float(traj.series(pool)[-1]) * 1000.0 for pool in _BEER_AROMA_TARGETS_MGL}


def _beer_aroma_levels_at_pitch(pitch_gpl: float) -> dict[str, float]:
    """All eight pools' finished levels at one inoculum, everything else the shipped frame."""
    frame = _beer_calibration_scenario()
    scenario = Scenario(
        name=frame.name,
        medium=frame.medium,
        initial={**frame.initial, "pitch_gpl": pitch_gpl},
        temperature_schedule=frame.temperature_schedule,
        duration_days=frame.duration_days,
    )
    compiled = compile_scenario(scenario)
    traj = compiled.run(t_eval=np.array([0.0, _BEER_CALIBRATION_DAYS * 24.0]))
    pools = [spec.pool for spec in FUSEL_SPECS] + [spec.pool for spec in ESTER_SPECS]
    return {pool: float(traj.series(pool)[-1]) * 1000.0 for pool in pools}


def _beer_growth_gain_time(pitch_gpl: float, fraction: float) -> float:
    """Days for ``fraction`` of the biomass gain, on a FIXED one-minute grid."""
    frame = _beer_calibration_scenario()
    scenario = Scenario(
        name=frame.name,
        medium=frame.medium,
        initial={**frame.initial, "pitch_gpl": pitch_gpl},
        temperature_schedule=frame.temperature_schedule,
        duration_days=5.0,
    )
    compiled = compile_scenario(scenario)
    grid = np.linspace(0.0, 5.0 * 24.0, 5 * 24 * 60 + 1)
    x = np.asarray(compiled.run(t_eval=grid).y, dtype=float)[compiled.schema.slice("X").start, :]
    peak = int(np.argmax(x))
    gain = (x[: peak + 1] - x[0]) / (x[peak] - x[0])
    return float(np.interp(fraction, gain, grid[: peak + 1])) / 24.0


def _beer_all_aroma_levels_mgl(data_dir: Path | None = None) -> dict[str, float]:
    """All EIGHT beer aroma pools, enumerated from the REGISTRIES rather than a literal list.

    ``_beer_aroma_levels_mgl`` above returns the seven pools D-224 hand-listed, and that dict is
    deliberately frozen at seven so D-225's arm-A measurement stays reproducible. Anything that
    needs the whole set -- like D-226's invariance guard -- must read
    ``ESTER_SPECS``/``FUSEL_SPECS``, which is the rule D-225 exists to establish: a literal list
    cannot fail to be complete, a registry read can.
    """
    scenario = _beer_calibration_scenario()
    compiled = (
        compile_scenario(scenario, data_dir=data_dir)
        if data_dir is not None
        else compile_scenario(scenario)
    )
    traj = compiled.run(t_eval=np.array([0.0, _BEER_CALIBRATION_DAYS * 24.0]))
    pools = [spec.pool for spec in FUSEL_SPECS] + [spec.pool for spec in ESTER_SPECS]
    return {pool: float(traj.series(pool)[-1]) * 1000.0 for pool in pools}


def _beer_levels_at_celsius(celsius: float) -> dict[str, float]:
    """All eight aroma pools at one isothermal temperature, registry-enumerated."""
    traj = compile_scenario(_beer_calibration_scenario(celsius)).run(
        t_eval=np.array([0.0, _BEER_CALIBRATION_DAYS * 24.0])
    )
    pools = [spec.pool for spec in FUSEL_SPECS] + [spec.pool for spec in ESTER_SPECS]
    return {pool: float(traj.series(pool)[-1]) * 1000.0 for pool in pools}


def _beer_data_dir_at(tmp_path, allow_out_of_band: bool = False, **values: float):
    """A copied parameter dir with named ``beer_generic.yaml`` values overridden.

    A copied dir rather than a resolved-map patch, for the D-214/D-223 reason: the arm has to be
    present at LOAD time so the compiled ProcessSet carries it. Every override here is a band
    EDGE of the parameter it names, so the value stays inside its own ``uncertainty`` and the
    ``Parameter`` schema accepts it without the band having to move too.

    ``allow_out_of_band`` (D-226) widens the named parameter's own band to admit the value, and
    is for ONE purpose: replaying a parameter's own RETIRED value, which by definition sits
    outside the band that was narrowed after it was retired (``q_sugar_max`` 0.5 against today's
    0.634-0.818, ``mu_max`` 0.098 against 0.053-0.075). It is not a way to reach a value the
    file rejects -- the widened band travels with the arm and is discarded with it.
    """
    import re
    import shutil

    dest = tmp_path / (
        "beer_arm_"
        + "_".join(f"{k}{v}" for k, v in sorted(values.items()))
        + ("_oob" if allow_out_of_band else "")
    )
    if dest.exists():
        return dest
    shutil.copytree(default_data_dir(), dest)
    path = dest / "beer_generic.yaml"
    text = path.read_text(encoding="utf-8")
    for name, value in values.items():
        match = re.search(rf"(^{name}:\n  value: )([0-9.eE+-]+)", text, re.MULTILINE)
        assert match is not None, f"{name} not found in beer_generic.yaml"
        text = text[: match.start(2)] + repr(float(value)) + text[match.end(2) :]
        if allow_out_of_band:
            # Rebuild the edges by SLICE, never by str.replace on the parsed float: a YAML
            # literal like `6.93e-5` round-trips through float() as `6.93e-05`, so a replace
            # silently no-ops and the arm would load with its ORIGINAL band.
            band = re.search(
                rf"^{name}:\n(?:.*\n)*?  uncertainty: \{{ low: ([0-9.eE+-]+), high: ([0-9.eE+-]+)",
                text,
                re.MULTILINE,
            )
            assert band is not None, f"{name} has no uncertainty band to widen"
            low, high = float(band.group(1)), float(band.group(2))
            text = (
                text[: band.start(1)]
                + repr(min(low, value))
                + text[band.end(1) : band.start(2)]
                + repr(max(high, value))
                + text[band.end(2) :]
            )
    path.write_text(text, encoding="utf-8")
    return dest


def _band_edges(parameters, name: str, expected_low: float, expected_high: float):
    """The parameter's own drawn band edges, plus a loud check that they are the ones measured.

    Reading them keeps an "at the band edge" arm at the edge when the band moves; asserting the
    values keeps the RATIOS this file pins meaningful, since they were measured at these edges.
    A band that moves therefore fails here with its own name rather than quietly turning an edge
    test into an interior one.
    """
    band = parameters[name].uncertainty
    assert band is not None, name
    assert band.low == pytest.approx(expected_low), (
        f"{name}'s low edge is {band.low}, not the {expected_low} the ratios below were "
        f"measured at -- re-measure them before re-pinning."
    )
    assert band.high == pytest.approx(expected_high), (
        f"{name}'s high edge is {band.high}, not the {expected_high} the ratios below were "
        f"measured at -- re-measure them before re-pinning."
    )
    return band.low, band.high


def test_the_finished_beer_lands_the_aroma_levels_its_rate_constants_are_defined_by():
    """Every one of the seven lands the level its own provenance says the k was set for.

    This is the guard whose absence let three beats move seven calibrated levels in silence. It
    is deliberately an equality against the TARGET rather than a snapshot of whatever the model
    currently produces: a snapshot re-pinned each time the engine moves records the drift instead
    of catching it.
    """
    # First: the seven targets above are not a second copy of the spec, they ARE the spec --
    # each one appears verbatim in its own parameter's `conditions:` field, so a re-sourcing that
    # moves the sentence fails here instead of leaving two numbers to diverge.
    beer = load_parameters(default_data_dir() / "beer_generic.yaml")
    for pool, literal in _BEER_AROMA_TARGET_LITERALS.items():
        conditions = beer[f"k_{pool}"].provenance.conditions
        assert literal in conditions, (
            f"k_{pool}'s `conditions:` no longer says it is set to land {literal}; the target in "
            f"this file is now a transcription rather than a reading of the spec."
        )
        assert float(literal.split()[0]) == _BEER_AROMA_TARGETS_MGL[pool]

    levels = _beer_aroma_levels_mgl()
    for pool, target in _BEER_AROMA_TARGETS_MGL.items():
        assert levels[pool] == pytest.approx(target, rel=0.01), (
            f"{pool} finishes at {levels[pool]:.3f} mg/L against the {target} mg/L its "
            f"`conditions:` field says its k is set to land. Re-anchor the k (the rate is linear "
            f"in it) -- do NOT relax this tolerance, and read D-224 before deciding the level is "
            f"the thing that should move."
        )


def test_the_beer_calibration_frames_inoculum_is_a_counted_pitch():
    """D-222's boundary conversion, applied to the OTHER beer scenario (D-228).

    The frame the eight aroma constants are defined in pitched a flat 1.0 g/L from D-99 until
    D-228. That number is not a pitching rate anybody measured: D-219 showed it is a RESIDUAL,
    a dry-yeast DOSING convention (a gram of product per litre) back-computed into ~100 pg/cell
    against a settled 40. Every counted pitch in this repo goes through
    :func:`~fermentation.units.cells_per_ml_to_pitch_gpl`, and this one now does too.

    Asserted as a RELATION to the conversion rather than as the resulting float, so a change to
    the settled per-cell mass moves the frame instead of silently contradicting it.
    """
    scenario = _beer_calibration_scenario()
    assert scenario.initial["pitch_gpl"] == cells_per_ml_to_pitch_gpl(
        BEER_COUNTED_PITCH_CELLS_PER_ML
    ), (
        "the aroma calibration frame's pitch is no longer the counted inoculum converted at the "
        "boundary; if it is back to a flat gram-per-litre, read D-219 before re-introducing one"
    )
    assert scenario.initial["pitch_gpl"] != _BEER_CALIBRATION_RETIRED_PITCH_GPL


def test_the_pitch_correction_leaves_every_calibrated_aroma_level_where_it_was(tmp_path):
    """The measurement that LICENSES D-228, run rather than cited.

    The correction is free only because the eight levels do not move, and that is not obvious:
    the synthesis half is exactly pitch-invariant (``int(mu*X*f_growth) dt = YAN /
    biomass_N_fraction`` does not mention the inoculum) but the three esters are STRIPPED by a
    first-order sink integrating against a moving pool, and a lighter pitch spreads formation
    over 1.5x more time. D-226 measured this invariance when the sink read the flux SHAPE;
    D-227 moved the sink onto evolved CO2, so it is re-measured here on the shipped rate law.

    Enumerated from the REGISTRIES (D-225's rule), not from a literal list -- a ninth pool must
    join this comparison the day it is added.
    """
    # Without this the test is VACUOUS on a tree where the frame has been reverted: both arms
    # would be the same run and it would pass on identity rather than on measurement, which is
    # indistinguishable from passing on the real thing (feedback-verify-the-restore-between-
    # mutation-arms). D-228's arm A is exactly that tree, and it passed there for this reason.
    assert _BEER_CALIBRATION_PITCH_GPL != _BEER_CALIBRATION_RETIRED_PITCH_GPL, (
        "the shipped frame and the retired arm are the same pitch, so the comparison below "
        "measures nothing; the frame has been reverted -- see D-228"
    )
    shipped = _beer_all_aroma_levels_mgl()
    retired = _beer_aroma_levels_at_pitch(_BEER_CALIBRATION_RETIRED_PITCH_GPL)
    assert set(shipped) == set(retired)
    worst = max(abs(retired[pool] / shipped[pool] - 1.0) for pool in shipped)
    assert worst < 1e-5, (
        f"the retired 1.0 g/L pitch and the counted one now disagree by {worst:.3e} relative "
        f"across the eight pools; D-228 measured 6e-7 and shipped the correction BECAUSE it was "
        f"invisible to every level. If this is real, the eight constants need re-anchoring and "
        f"that is a beat, not a tolerance change."
    )


def test_the_calibration_frames_growth_window_is_no_longer_a_single_day():
    """The claim D-228 corrects, measured on the frame that produced it.

    D-226 §8 and D-227 §10 both name "growth stops dead at day 0.92 with 81 % of the sugar still
    to ferment" as the inherited limitation behind beer's aroma taper, and the memory ledger
    carries it as the next thing to fix. It is a statement about THIS frame's retired inoculum:
    nitrogen-limited growth reaches a ceiling that is fixed in ABSOLUTE terms, so a heavier pitch
    reaches it sooner. At the counted pitch the same wort takes 1.5x longer to get there.

    Read on a FIXED grid, because the quantity is a time on a steep limb
    (feedback-read-a-fast-curve-on-a-fixed-grid).
    """
    counted = _beer_growth_gain_time(_BEER_CALIBRATION_PITCH_GPL, 0.90)
    retired = _beer_growth_gain_time(_BEER_CALIBRATION_RETIRED_PITCH_GPL, 0.90)
    assert retired == pytest.approx(0.785, abs=0.01), (
        f"at the RETIRED 1.0 g/L pitch 90 % of the growth takes {retired:.3f} d; D-228 measured "
        f"0.785, which is the number D-226 §8 reported as 'day 0.92' (that is its 99 % point)"
    )
    assert counted > 1.15, (
        f"at the counted pitch 90 % of the growth takes {counted:.3f} d, back inside the single "
        f"day D-228 measured it out of (1.199). The window is what the correction buys -- the "
        f"levels are invariant, so this is the only assert that can see the change"
    )
    assert counted / retired == pytest.approx(1.528, abs=0.03)


def test_no_drawn_speed_knob_moves_the_five_higher_alcohols_and_barely_moves_the_esters(
    tmp_path,
):
    """The property D-226 exists to buy, and the exact size of the half it does NOT buy.

    This REPLACES D-224's ``..._moves_only_the_half_of_the_aroma_set_it_is_coupled_to``, whose
    whole premise was the defect: it asserted that ``mu_max`` moves the five Ehrlich alcohols
    (x1.094/x0.774) and ``q_sugar_max`` moves ethyl acetate (x1.111/x0.897), and it was RIGHT
    about the engine it was written for. Those sensitivities are what made four beats of
    ferment-speed re-anchoring silently move eight calibrated aroma levels. Under the
    growth-extent coupling neither knob reaches the synthesis term at all, because
    ``int(mu*X*f_growth) dt`` is ``YAN / biomass_N_fraction`` -- a conservation identity, not an
    approximation.

    **The five higher alcohols are EXACTLY invariant and the three esters are NOT, and the
    difference is a sink rather than a source.** The alcohols have no volatilisation Process, so
    nothing downstream re-introduces a speed dependence: they measure 1.00000 at both knobs' band
    edges AND at their retired values. The esters are stripped by
    :class:`~fermentation.core.kinetics.byproducts.EsterVolatilization`, which still rides the
    fermentative flux -- a slower beer evolves its CO2 over more biomass-hours and strips more
    away. That residue is REAL and is pinned here rather than rounded off, because a test named
    for invariance that quietly tolerated 4.5 % would be the same defect in the other direction.

    Both the retired values are included deliberately: they are what makes this a statement
    about the drift that actually happened rather than about a band nobody moved.
    """
    nominal = _beer_all_aroma_levels_mgl()
    # The edges are READ, never transcribed (D-224's hardening, kept).
    beer = load_parameters(default_data_dir() / "beer_generic.yaml")
    mu_lo, mu_hi = _band_edges(beer, "mu_max", 0.053, 0.075)
    q_lo, q_hi = _band_edges(beer, "q_sugar_max", 0.634, 0.818)

    arms = {
        "mu_max low": _beer_all_aroma_levels_mgl(_beer_data_dir_at(tmp_path, mu_max=mu_lo)),
        "mu_max high": _beer_all_aroma_levels_mgl(_beer_data_dir_at(tmp_path, mu_max=mu_hi)),
        "q_sugar_max low": _beer_all_aroma_levels_mgl(
            _beer_data_dir_at(tmp_path, q_sugar_max=q_lo)
        ),
        "q_sugar_max high": _beer_all_aroma_levels_mgl(
            _beer_data_dir_at(tmp_path, q_sugar_max=q_hi)
        ),
    }
    for arm, levels in arms.items():
        for pool in _BEER_EHRLICH_POOLS:
            assert levels[pool] / nominal[pool] == pytest.approx(1.0, abs=1e-4), (
                f"{pool} moved {levels[pool] / nominal[pool]:.6f}x at {arm}. The five Ehrlich "
                f"alcohols are coupled to growth EXTENT and have no stripping sink, so this "
                f"ratio is a conservation identity and should read 1.000000. If it does not, "
                f"the producer is reading a speed knob again -- read D-226 before re-pinning."
            )
        # The esters DO move, through the sink, and by how much is the claim. D-227 tightened
        # this from 5 % to 0.7 %: the sink's driver became the evolving CO2 itself, whose
        # integral to dryness is fixed by the sugar there was. The bound is set by the WORST
        # band-edge arm, which is `mu_max` high at 0.543 % -- not by `q_sugar_max`, whose two
        # edges now sit inside 0.33 %. That the surviving residue belongs to the GROWTH knob
        # rather than the uptake knob is the mechanism showing through, and it is asserted
        # pool-by-pool below rather than left to this blanket bound.
        assert levels["ethyl_acetate"] / nominal["ethyl_acetate"] == pytest.approx(
            1.0, abs=0.007
        ), arm

    # The residue, pinned per edge and per pool. Sizes RE-MEASURED at D-227; the SIGN is
    # still reversed from D-224's (a slower ferment packages LESS ester, because it strips for
    # longer), but the size is now a seventh of what D-226 shipped.
    for pool in ("ethyl_acetate", "isoamyl_acetate", "ethyl_hexanoate"):
        assert arms["q_sugar_max low"][pool] / nominal[pool] == pytest.approx(0.9970, abs=0.002), (
            f"{pool}'s residual q_sugar_max sensitivity is not the measured one. Under D-224's "
            f"flux coupling this edge read 1.111 (MORE ester from a slower beer); it now reads "
            f"0.997 (LESS), because synthesis is speed-invariant and the stripping sink now "
            f"rides the CO2 the ferment evolves, whose integral to dryness is fixed by the "
            f"sugar there was. D-226 shipped 0.955 here; if this reads that again, the sink "
            f"has gone back to the flux SHAPE."
        )
        assert arms["q_sugar_max high"][pool] / nominal[pool] == pytest.approx(1.0037, abs=0.002)
        # mu_max is nearly inert on the esters too -- it never reached them under either
        # coupling, and D-227 BARELY MOVES IT (1.003 -> 1.0025). That asymmetry is the point of
        # the pair: `q_sugar_max` acted through the driver's SIZE, which the CO2 integral now
        # fixes, while `mu_max` acts through its TIMING, which nothing here fixes. What is left
        # of the residue is a timing residue, and it is now the larger of the two.
        assert arms["mu_max low"][pool] / nominal[pool] == pytest.approx(1.0025, abs=0.002)
        assert arms["mu_max high"][pool] / nominal[pool] == pytest.approx(0.9937, abs=0.002)

    # And the historical arms: what D-223's and D-211's re-anchorings WOULD have done to these
    # levels under this coupling. This is the beat's headline as a number.
    retired_q = _beer_all_aroma_levels_mgl(
        _beer_data_dir_at(tmp_path, allow_out_of_band=True, q_sugar_max=0.5)
    )
    retired_mu = _beer_all_aroma_levels_mgl(
        _beer_data_dir_at(tmp_path, allow_out_of_band=True, mu_max=0.098)
    )
    for pool in _BEER_EHRLICH_POOLS:
        # D-211's mu_max change multiplied every one of these by 2.87 under the flux coupling.
        assert retired_mu[pool] / nominal[pool] == pytest.approx(1.0, abs=1e-4), pool
        assert retired_q[pool] / nominal[pool] == pytest.approx(1.0, abs=1e-4), pool
    # D-223's q_sugar_max 0.5 -> 0.72 moved packaged ethyl acetate 20.9 % under D-224's flux
    # coupling and 13.9 % under D-226's. Under D-227 it moves it 0.73 %: a 19x reduction on the
    # single largest historical drift this file records, and the reason it is nearly zero rather
    # than merely smaller is that total evolved CO2 does not depend on how fast the sugar went.
    assert retired_q["ethyl_acetate"] / nominal["ethyl_acetate"] == pytest.approx(0.9927, abs=0.005)
    # `mu_max`'s retired arm is UNIMPROVED (0.987 -> 0.990) and that is not a disappointment,
    # it is the mechanism showing through: growth rate re-times the CO2 stream relative to when
    # the ester is made, and a first-order sink integrated against a moving pool reads timing.
    assert retired_mu["ethyl_acetate"] / nominal["ethyl_acetate"] == pytest.approx(
        0.9901, abs=0.005
    )


def test_the_co2_rate_helper_is_bitwise_the_uptake_processs_own_co2(params):
    """The D-180 discipline, extended to the quantity D-227 added.

    ``EsterVolatilization`` strips on the CO2 the ferment evolves; ``SugarUptakeToEthanolCO2``
    books that same CO2 into the state. If the two ever computed it differently, the sink would
    be riding a gas stream the solver never ran -- the identical failure mode D-106 caught for
    the uptake RATES, with the identical absence of a symptom. So the helper is not "the same
    formula written twice", it is the Process's own arithmetic called from one place, and this
    is the test that says so.

    BITWISE rather than to a tolerance, for the reason the sibling test in ``test_organic_acids``
    gives: anything less than bit-equality means the arithmetic moved
    [[feedback-pin-tolerance-vs-solver-tolerance]].
    """
    uptake = SugarUptakeToEthanolCO2()
    for schema_fn, states in (
        (wine_schema, ([200.0], [40.0], [0.0])),
        (beer_schema, ([15.0, 70.0, 15.0], [0.0, 40.0, 12.0], [0.0, 0.0, 5.0], [0.0, 0.0, 0.0])),
    ):
        schema = schema_fn()
        for sugars in states:
            for x in (2.0, 0.0):
                y = schema.zeros()
                y[schema.slice("X")] = x
                y[schema.slice("S")] = np.asarray(sugars, dtype=float)
                y[schema.slice("T")] = 293.15
                d = uptake.derivatives(0.0, y, schema, params)
                assert fermentative_co2_rate(y, schema, params) == float(
                    d[schema.slice("CO2")][0]
                ), (
                    f"the shared CO2-rate helper and the uptake Process disagree for "
                    f"S={sugars}, X={x}: the ester sink is stripping on a gas stream the "
                    f"ferment is not evolving"
                )


def test_the_ester_sink_rides_evolved_co2_and_not_the_flux_shape_that_stood_in_for_it(params):
    """What D-227 changed, asserted as a PROPERTY of the two drivers rather than as a level.

    Every aroma LEVEL guard in this file is blind to this change by construction: the beat
    re-anchored ``k_ester_volatil`` so the calibration frame is unchanged, which is what makes
    its other consequences attributable -- and which means a revert of the Process alone would
    move the levels, but a revert of Process AND constant together would not. Only a test that
    names the driver can see that. [[feedback-prefer-the-variant-your-guards-can-see]]

    Two claims, and the second is the one that makes the change a mechanism rather than a rename:

    * the sink's derivative is proportional to ``fermentative_co2_rate``, not to
      ``fermentative_flux_shape``; and
    * in BEER the two are not proportional to each other, because three sugars under sequential
      catabolite repression evolve CO2 in a different shape from a single Monod on total sugar.
      In WINE they ARE proportional -- one sugar, no repression -- which is why wine's
      re-anchoring came out exactly equal to the analytic unfolding and beer's did not.
    """
    beer, wine = beer_schema(), wine_schema()
    sink = EsterVolatilization()

    # (a) the derivative tracks the CO2 rate, at a state where the two drivers DISAGREE.
    y = beer.zeros()
    y[beer.slice("X")] = 2.0
    y[beer.slice("S")] = np.asarray([5.0, 70.0, 15.0], dtype=float)  # glucose nearly gone
    y[beer.slice("T")] = 293.15
    y[beer.slice("ethyl_acetate")] = 0.02
    d = sink.derivatives(0.0, y, beer, params)
    f_gas = arrhenius_factor(293.15, params["E_a_uptake"], params["T_ref"])
    f_part = arrhenius_factor(293.15, params["dH_ester_volatil"], params["T_ref"])
    co2_rate = fermentative_co2_rate(y, beer, params)
    expected = params["k_ester_volatil"] * co2_rate * f_gas * f_part * 0.02
    assert float(d[beer.slice("ethyl_acetate")][0]) == pytest.approx(-expected, rel=1e-12)

    # (b) beer's two drivers are NOT proportional; wine's are. Measured as the spread of the
    # ratio across states, so it cannot be satisfied by a single lucky point.
    def ratios(schema: StateSchema, states: tuple[list[float], ...]) -> list[float]:
        out = []
        for sugars in states:
            yy = schema.zeros()
            yy[schema.slice("X")] = 2.0
            yy[schema.slice("S")] = np.asarray(sugars, dtype=float)
            yy[schema.slice("T")] = 293.15
            shape = fermentative_flux_shape(yy, schema, params["K_sugar_uptake"])
            out.append(fermentative_co2_rate(yy, schema, params) / shape)
        return out

    beer_r = ratios(beer, ([15.0, 70.0, 15.0], [0.5, 70.0, 15.0], [0.0, 5.0, 15.0]))
    wine_r = ratios(wine, ([200.0], [50.0], [5.0]))
    # 1.150 MEASURED across these three states, not a round number chosen to pass: the
    # per-species CO2 yields alone span 1.071 (0.4886 glucose / 0.5235 maltotriose, the
    # hydrolysis water) and sequential repression supplies the rest.
    assert max(beer_r) / min(beer_r) == pytest.approx(1.150, abs=0.01), beer_r
    assert max(beer_r) / min(beer_r) > 1.05, (
        f"beer's CO2 rate and the retired flux shape came out proportional (ratios {beer_r}). "
        f"They must not be: catabolite repression and the three sugars' different CO2 yields "
        f"are what make D-227 a change of mechanism in beer rather than the pure "
        f"reparameterisation it is in wine."
    )
    assert max(wine_r) / min(wine_r) == pytest.approx(1.0, abs=1e-12), (
        f"wine's two drivers stopped being proportional (ratios {wine_r}). Wine has ONE sugar "
        f"slot, so the ratio is the constant q_sugar_max*co2_yield*scale that D-227 unfolded "
        f"out of k_ester_volatil -- and wine's re-anchoring factor was DERIVED from it."
    )
    # And that constant is the number wine's k moved by, which is what makes the wine half a
    # derivation rather than a fit.
    assert wine_r[0] == pytest.approx(1.0 / 2.5252809, rel=1e-6)


def test_both_consumers_of_the_flux_helper_carry_the_uptake_arrhenius_exactly_once():
    """The hazard `fermentative_uptake_rates`' own docstring names, now that it has two callers.

    That helper returns UNMODIFIED rates, so anything built on it must acquire
    ``arrh(E_a_uptake)`` somehow or it books its yield against a flux the solver never ran --
    the D-32 coupling. Its two consumers acquire it by DIFFERENT routes:
    ``OrganicAcidExcretion`` is an extra target of ``ArrheniusTemperature.for_uptake`` and has
    its whole derivative scaled; ``EsterVolatilization`` is not a target and applies the factor
    itself, because it needs a second independent Arrhenius (the Henry partition) beside it and
    a modifier contributes only one.

    **This guard is for ATTRIBUTION, not detection, and the difference was measured rather than
    assumed.** The first draft of this docstring claimed the temperature guard would stay GREEN
    through a doubling. Both mutations were then actually run, and it does NOT: adding the sink
    to ``for_uptake`` while keeping ``f_gas`` takes the packaged 15/25 C ester span 6.8953 ->
    5.2784, and deleting ``f_gas`` takes it to 8.7751; ``test_beers_aroma_temperature_response_...``
    fails on both. What it CANNOT do is say why -- its message names ``E_a_esters`` and sends the
    reader to a parameter that did not move. This test reports "carries arrh(E_a_uptake) 2.0000
    times" and names the two edits that produce it. Verified RED in BOTH directions with those
    exact numbers [[feedback-verify-an-xfail-fails-for-its-stated-reason]].

    So the exponent is measured directly -- perturb ``E_a_uptake`` and read how many times the
    factor appears -- rather than inferred from the wiring, which is what makes the message
    specific enough to be worth having beside a guard that already goes red.
    """
    scenario = Scenario(
        name="exposure",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 1.0,
        },
        # Well off T_ref (20 C), or every Arrhenius factor is 1.0 and the exponent is undefined
        # -- the same blind spot D-226 arm C found in the level guards.
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=28.0)],
        duration_days=14.0,
    )
    compiled = compile_scenario(scenario)
    schema = compiled.schema
    y = compiled.y0.copy()
    y[schema.slice("X")] = 2.0
    y[schema.slice("S")] = np.asarray([10.0, 60.0, 14.0], dtype=float)
    y[schema.slice("T")] = 28.0 + 273.15
    y[schema.slice("ethyl_acetate")] = 0.02

    base = compiled.param_values["E_a_uptake"]
    t_ref = compiled.param_values["T_ref"]
    temp = 28.0 + 273.15

    def total_at(e_a: float) -> FloatArray:
        params = dict(compiled.param_values)
        params["E_a_uptake"] = e_a
        return compiled.process_set.total_derivatives(0.0, y, params)

    lo, hi = base * 0.9, base * 1.1
    d_lo, d_hi = total_at(lo), total_at(hi)
    one_factor = arrhenius_factor(temp, hi, t_ref) / arrhenius_factor(temp, lo, t_ref)
    assert one_factor > 1.05  # the perturbation has to bite, or the test is vacuous

    # `ethyl_acetate_gas` is the sink's OWN slot: the liquid pool is a net of synthesis and
    # stripping and would report a blend, which is why the twin is read instead.
    for label, slot in (
        ("the ester sink", "ethyl_acetate_gas"),
        ("the uptake Process itself", "CO2"),
    ):
        ratio = float(d_hi[schema.slice(slot)][0]) / float(d_lo[schema.slice(slot)][0])
        exponent = math.log(ratio) / math.log(one_factor)
        assert exponent == pytest.approx(1.0, abs=1e-3), (
            f"{label} carries arrh(E_a_uptake) {exponent:.4f} times, not once. 2.0 means the "
            f"sink was added to ArrheniusTemperature.for_uptake while keeping its own f_gas; "
            f"0.0 means f_gas was deleted. Read D-227 and the helper's docstring."
        )


def test_beers_aroma_temperature_response_is_the_one_its_activation_energies_were_solved_for():
    """Beer's packaged 15/25 C spans — the ONLY axis `E_a_esters` and `E_a_fusels` control.

    **This guard exists because a mutation came back GREEN.** D-226 arm C reverted both beer
    activation energies (values AND bands) to their pre-D-226 values, which more than DOUBLES the
    packaged ester temperature span, 6.90x -> 13.50x over 15-25 C. The full suite returned exactly
    ONE failure, and it was a band-overlap assert on a different pair. Every level guard is blind
    here by construction: the calibration frame is 20 C, which is `T_ref`, where both Arrhenius
    factors are exactly 1.0. Wine has had `test_integrated_wine_aroma_temperature_directions` for
    this since D-21; beer had nothing. A GREEN mutation is what owes a guard.

    The two numbers asserted are not snapshots -- they are the quantities the two constants were
    SOLVED to hold. D-226 changed what both parameters mean (under the retired flux coupling the
    observable's apparent activation energy was `E_a - E_a_uptake`; under growth-extent coupling
    the growth factor cancels against the nitrogen limit and it is `E_a` itself), and each value
    was then bisected to reproduce the span the retired shape delivered, so that the temperature
    axis is unchanged BY CONSTRUCTION and D-226's other consequences are attributable to the
    coupling alone. If either span moves, that construction has been broken.

    Read the fusel span off a pool the ester re-route does NOT debit. `isoamyl_alcohol` measures
    1.0395x rather than 1.2182x because D-115's re-route spends it on isoamyl acetate, whose own
    response is steep -- reading the span off that pool is how the first pass of this calibration
    solved `E_a_fusels` to 4,707 J/mol instead of 14,100.
    """
    warm = {}
    for celsius in (15.0, 20.0, 25.0):
        warm[celsius] = _beer_levels_at_celsius(celsius)

    def span(pool: str) -> float:
        return warm[25.0][pool] / warm[15.0][pool]

    # (1) DIRECTION -- the sourced claim (D-19/D-21): a warmer ale carries more of both families.
    for pool in _BEER_EHRLICH_POOLS + ("ethyl_acetate", "isoamyl_acetate", "ethyl_hexanoate"):
        assert warm[15.0][pool] < warm[20.0][pool], pool

    # (2) MAGNITUDE -- the numbers the two activation energies were solved to reproduce.
    # The four Ehrlich pools the ester re-route does not debit share one shape, so one number.
    for pool in ("propanol", "isobutanol", "active_amyl_alcohol", "2_phenylethanol"):
        assert span(pool) == pytest.approx(1.2182, rel=0.005), (
            f"{pool}'s 15/25 C span is {span(pool):.4f}, not the 1.2182 `E_a_fusels` was solved "
            f"to hold. That is the span the retired flux coupling delivered (E_a_fusels 70,000 "
            f"minus E_a_uptake 55,100); D-226 re-anchored the parameter to 14,100 to preserve it. "
            f"Read D-226 Sec 6 before re-pinning -- the parameter's own 'Q10 ~ 2.6' note describes "
            f"the RATE CONSTANT, not this observable."
        )
    # The two mechanically-locked esters share `E_a_esters` and the stripping constant exactly.
    assert span("ethyl_acetate") == pytest.approx(6.8953, rel=0.005), (
        f"ethyl acetate's 15/25 C span is {span('ethyl_acetate'):.4f}. `E_a_esters` was solved "
        f"(151,688, shipped 152,000) to hold the 6.8596 the retired biomass-hour coupling "
        f"delivered; the shipped 6.8953 is that plus the 0.52 % rounding cost. Reverting the "
        f"parameter to its retired 200,000 takes this to 13.50 -- which is precisely the mutation "
        f"the whole suite failed to notice before this test existed."
    )
    assert span("ethyl_hexanoate") == pytest.approx(span("ethyl_acetate"), rel=1e-4), (
        "the two esters are MECHANICALLY LOCKED (same rate law, same E_a_esters, same stripping "
        "constant), so their temperature spans are equal by construction, not by calibration."
    )
    # Isoamyl acetate is steeper than either: it multiplies the ester shape by a precursor pool
    # that is itself temperature-responsive.
    assert span("isoamyl_acetate") == pytest.approx(7.6002, rel=0.005)
    assert span("isoamyl_acetate") > span("ethyl_acetate")

    # (3) THE OTHER SOURCED ORDERING (D-19): warmth shifts the balance TOWARD esters, so the
    # ester/fusel ratio rises with temperature. Asserted as the ratio it is about, rather than as
    # the `E_a_fusels < E_a_esters` inequality, which is the same claim in the parameters'
    # coordinates and is guarded there.
    cold_ratio = warm[15.0]["ethyl_acetate"] / warm[15.0]["propanol"]
    hot_ratio = warm[25.0]["ethyl_acetate"] / warm[25.0]["propanol"]
    assert hot_ratio > cold_ratio, (
        f"the ester/fusel ratio FALLS with temperature ({cold_ratio:.4f} -> {hot_ratio:.4f}). "
        f"beer_generic.yaml states the opposite as a sourced direction."
    )


def test_beers_isoamyl_alcohol_stays_below_its_only_sourced_threshold_across_the_growth_band(
    tmp_path,
):
    """The SENSORY statement of what D-224 repaired, which is sharper than any of the means.

    ``isoamyl_alcohol`` is the one beer pool of the five with a sourced in-matrix threshold
    (Meilgaard 1975, ~50 mg/L, sensory.yaml), and beer_generic.yaml states the finished level
    must read OAV ~0.6 against it, because fusel character is a DEFECT in beer at elevated
    levels rather than a normal descriptor of it. The shipped model ran 48.261 mg/L -- OAV
    0.965, within 3.6 % of reporting a solventy note this file says an ale must not have -- and
    at ``mu_max``'s LOW edge it ran 52.809, i.e. over. That is a live descriptor defect, not a
    calibration nicety, and it is what makes the re-anchoring a fix rather than a preference.
    """
    # ug/L in sensory.yaml, mg/L here -- read from the file rather than written as a literal,
    # so a re-sourcing of the threshold moves this guard with it.
    threshold_mgl = (
        load_parameters(default_data_dir() / "sensory.yaml").resolve()[
            "threshold_isoamyl_alcohol_beer"
        ]
        / 1000.0
    )
    assert threshold_mgl == pytest.approx(50.0)
    for mu in (0.053, 0.058, 0.075):
        arm = _beer_aroma_levels_mgl(_beer_data_dir_at(tmp_path, mu_max=mu))
        oav = arm["isoamyl_alcohol"] / threshold_mgl
        assert oav < 1.0, (
            f"at mu_max={mu} beer finishes {arm['isoamyl_alcohol']:.3f} mg/L isoamyl alcohol, "
            f"OAV {oav:.3f} -- it is claiming a solventy/alcoholic note beer_generic.yaml says a "
            f"sound ale does not have."
        )
    # The nominal reads the 0.6 the parameter file prints, not merely something under 1.0.
    assert _beer_aroma_levels_mgl()["isoamyl_alcohol"] / threshold_mgl == pytest.approx(
        0.602, abs=0.01
    )


@pytest.mark.parametrize(
    ("medium_file", "stated"),
    [
        (
            "beer_generic.yaml",
            {
                "k_propanol": (0.2, 5.0),
                "k_isobutanol": (0.3, 3.0),
                "k_active_amyl_alcohol": (0.3, 3.0),
                "k_isoamyl_alcohol": (0.3, 3.0),
                "k_2_phenylethanol": (0.3, 3.0),
            },
        ),
        (
            "wine_generic.yaml",
            {
                "k_propanol": (0.2, 4.0),  # wine's propanol note says x0.2/x4, not beer's x0.2/x5
                "k_isobutanol": (0.3, 3.0),
                "k_active_amyl_alcohol": (0.3, 3.0),
                "k_isoamyl_alcohol": (0.3, 3.0),
                "k_2_phenylethanol": (0.3, 3.0),
            },
        ),
    ],
)
def test_the_ehrlich_bands_are_the_multiple_of_their_nominal_their_own_notes_state(
    medium_file, stated
):
    """D-99 shipped BEER's five bands as x0.3/x3 (x0.2/x5 for propanol) of a value 2.05x BELOW
    the nominal shipped beside them -- x0.145/x1.45 of the value actually in force. One
    draughting error across five entries: the nominal was fitted to land its target and the band
    was left where it was. D-224 corrects both halves; this is the arithmetic that would have
    caught it, and it is cheap enough that there was never a reason not to have had it.

    **Both media, because a guard is only as broad as the registry it names.** 955ebbc shipped
    wine's five in the same commit with the same note text, so the question "does wine carry the
    same error" is not answered by wine's LEVELS being right -- levels say nothing about
    multipliers. Measured at the follow-up: wine's five are CORRECT (x0.3/x3, and x0.2/x4 for its
    own propanol), so the draughting error is beer-only. That is a finding, and it is the reason
    this test is parametrized rather than left pointing at one file.
    """
    # rel=0.02: every value and edge in these files is rounded to 3 significant figures, so an
    # exactly-stated multiplier lands within ~1.5 % of it (wine's isobutanol is the worst at
    # 0.3037). The error this catches is 2.05x, two orders above the rounding.
    parameters = load_parameters(default_data_dir() / medium_file)
    for name, (low_mult, high_mult) in stated.items():
        param = parameters[name]
        assert param.uncertainty is not None
        assert param.uncertainty.low / param.value == pytest.approx(low_mult, rel=0.02), name
        assert param.uncertainty.high / param.value == pytest.approx(high_mult, rel=0.02), name


def test_the_ethyl_acetate_band_spans_its_sourced_ale_range_and_its_top_reaches_the_threshold(
    tmp_path,
):
    """The band is COMPUTED to span 10-30 mg/L, and the sensory consequence is pinned, not glossed.

    This is the guard D-224 did not have when it first re-anchored `k_ethyl_acetate`. Rescaling
    the pre-existing multipliers (x0.5/x1.5455) with the new nominal put the band top at 31.03
    mg/L -- over this molecule's own sourced 10-30 ale range AND over the ~30 mg/L Meilgaard beer
    threshold -- while the note beside it still made the OAV claim at the nominal ALONE. That is
    the exact defect this beat's headline is about, committed inside the fix for it.

    So the edges are computed from runs instead: they span 10.00-30.02 mg/L. The top therefore
    touches OAV 1.001, which is CORRECT rather than tolerated -- the source says a sound ale sits
    *at or below* the threshold, so the top of the molecule's own range should reach it. It is
    still a change: before D-224 the band spanned 7.93-24.52 mg/L (top OAV 0.817) and could not
    reach the threshold at all, because the nominal was 21 % low. Jointly with `q_sugar_max`'s
    low edge the reachable maximum is 33.36 mg/L, OAV 1.112.
    """
    beer = load_parameters(default_data_dir() / "beer_generic.yaml")
    band = beer["k_ethyl_acetate"].uncertainty
    assert band is not None
    threshold_mgl = (
        load_parameters(default_data_dir() / "sensory.yaml").resolve()[
            "threshold_ethyl_acetate_beer"
        ]
        / 1000.0
    )

    at_low = _beer_aroma_levels_mgl(_beer_data_dir_at(tmp_path, k_ethyl_acetate=band.low))
    at_high = _beer_aroma_levels_mgl(_beer_data_dir_at(tmp_path, k_ethyl_acetate=band.high))
    assert at_low["ethyl_acetate"] == pytest.approx(10.0, abs=0.05)
    assert at_high["ethyl_acetate"] == pytest.approx(30.0, abs=0.05)

    nominal = _beer_aroma_levels_mgl()
    assert nominal["ethyl_acetate"] / threshold_mgl == pytest.approx(0.669, abs=0.01)
    assert at_high["ethyl_acetate"] / threshold_mgl == pytest.approx(1.001, abs=0.01)

    # The joint corner with the OTHER drawn knob this pool reads. Stated as a number so a future
    # re-anchoring cannot quietly move it: a nominal-only OAV claim beside a drawn band is what
    # D-224 exists to stop.
    #
    # D-226 MOVED THIS CORNER AND REVERSED ITS DIRECTION, which is worth more than the number.
    # Under the flux coupling a SLOWER beer made MORE ester, so `q_sugar_max`'s LOW edge was the
    # bad corner and it reached OAV 1.112 -- above the threshold, and reachable in a drawn
    # ensemble. Under growth-extent coupling synthesis is speed-invariant and only the stripping
    # sink still reads the flux, so a slower beer now strips MORE and packages LESS: the low edge
    # falls and the worst corner is the HIGH edge instead. Both are asserted, because "the
    # reachable maximum" is the claim and D-226 changed which edge holds it.
    #
    # D-227 SHRANK THE WHOLE CORNER without moving its owner: the sink's driver became the
    # evolving CO2, whose integral to dryness `q_sugar_max` cannot change, so the two edges
    # close from 0.9555/1.0428 onto 0.9978/1.0039. The high edge still holds the maximum and it
    # is still (just) over the threshold, which is the honest reading -- the source says a sound
    # ale sits at or below it, and the top of this molecule's own range should reach it.
    q_lo, q_hi = _band_edges(beer, "q_sugar_max", 0.634, 0.818)
    slow = _beer_aroma_levels_mgl(
        _beer_data_dir_at(tmp_path, k_ethyl_acetate=band.high, q_sugar_max=q_lo)
    )
    fast = _beer_aroma_levels_mgl(
        _beer_data_dir_at(tmp_path, k_ethyl_acetate=band.high, q_sugar_max=q_hi)
    )
    assert slow["ethyl_acetate"] / threshold_mgl == pytest.approx(0.9978, abs=0.005)
    assert fast["ethyl_acetate"] / threshold_mgl == pytest.approx(1.0039, abs=0.005)
    assert max(slow["ethyl_acetate"], fast["ethyl_acetate"]) / threshold_mgl < 1.0428, (
        "D-226 REDUCED the reachable ethyl-acetate maximum from 1.112 to 1.043 and D-227 "
        "reduced it again to 1.004; a rise means the sink has stopped riding evolved CO2."
    )
    # The SPREAD between the two corners is the quantity D-227 acts on, and pinning it here
    # rather than only the two endpoints is deliberate: two endpoint pins with a 0.005 tolerance
    # would still pass if the band re-opened by 1 %, and the spread is what the beat claims.
    assert abs(fast["ethyl_acetate"] - slow["ethyl_acetate"]) / threshold_mgl < 0.010, (
        "the two q_sugar_max corners have re-opened; D-227 closed them from 0.087 to 0.006 "
        "of a threshold by making the sink read evolved CO2 rather than the flux shape."
    )


def test_the_isoamyl_acetate_band_is_rescaled_with_its_nominal_rather_than_recomputed():
    """Its band gets the OPPOSITE treatment from ethyl acetate's, and from the Ehrlich five.

    Three rules in one file, and keeping them straight is the point of this test and its two
    siblings:

    * the five Ehrlich `k` -- the stated MULTIPLIER (x0.3/x3) and the actual one (x0.145/x1.45)
      disagreed by 2.05x with nothing supporting the actual, so the multiplier is CORRECTED;
    * `k_ethyl_acetate` -- its band carries a sensory threshold at the top of its stated range,
      so the edges are COMPUTED from runs to span exactly that range;
    * `k_isoamyl_acetate` -- its stated span and its actual one agree, and its threshold is
      crossed at the NOMINAL by design (the banana note is an ale signature), so there is no
      nominal-below-threshold claim to protect. It is RESCALED with the nominal, the D-97/D-99
      convention, and the band top moved DOWN at D-224 (5.87 -> 4.47 mg/L) rather than up.
    """
    beer = load_parameters(default_data_dir() / "beer_generic.yaml")
    param = beer["k_isoamyl_acetate"]
    assert param.uncertainty is not None
    # The multipliers D-96/D-97/D-99 shipped, carried through every re-anchoring since.
    assert param.uncertainty.low / param.value == pytest.approx(0.3286, rel=0.005)
    assert param.uncertainty.high / param.value == pytest.approx(2.0896, rel=0.005)


# -- D-224 follow-up: the two surfaces the first pass of these guards could not see ---------


def test_wine_also_lands_the_five_ehrlich_levels_its_constants_are_defined_by():
    """The wine half of D-224's control, as a guard rather than a paragraph.

    D-224's argument rests on wine being untouched -- ``mu_max`` lives per medium, so the beer
    growth-rate changes that moved beer's higher alcohols by x2.87 could not reach wine, and wine
    still landing its own published means is what makes the beer drift a defect with a mechanism
    rather than a drift in the shared Process. That control was measured in two pools of five and
    written into a record; here it is all five, read from the same ``conditions:`` sentences.
    """
    wine = load_parameters(default_data_dir() / "wine_generic.yaml")
    # (value, the literal its own `conditions:` field spells it with)
    targets = {
        "propanol": (24.0, "24 mg/L"),
        "isobutanol": (33.0, "33.0 mg/L"),
        "active_amyl_alcohol": (70.1, "70.1 mg/L"),
        "isoamyl_alcohol": (172.0, "172 mg/L"),
        "2_phenylethanol": (28.7, "28.7 mg/L"),
    }
    for pool, (_, literal) in targets.items():
        conditions = wine[f"k_{pool}"].provenance.conditions
        assert literal in conditions, pool

    scenario = Scenario(
        name="wine-aroma-calibration",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )
    traj = compile_scenario(scenario).run(t_eval=np.array([0.0, 21.0 * 24.0]))
    for pool, (target, _) in targets.items():
        level = float(traj.series(pool)[-1]) * 1000.0
        assert level == pytest.approx(target, rel=0.02), (
            f"wine's {pool} finishes at {level:.3f} mg/L against the {target} its k is set to "
            f"land. If this fails alongside its beer twin, the drift is in the shared Process; "
            f"if it fails alone, read D-224 -- wine being untouched is that record's control."
        )


def test_beers_solventy_descriptor_axis_changed_owner_and_fell_at_d224():
    """The OUTPUT-level consequence of re-anchoring seven constants, which nothing else asserts.

    Under the D-95 MAX rule a descriptor reads its loudest contributing pool, so *which molecule
    owns an axis* is a modelled claim about the beer, not an internal detail. Re-anchoring moved
    isoamyl alcohol down and ethyl acetate up past each other, so beer's ``solventy`` axis
    **changed owner** -- ``isoamyl_alcohol`` at OAV 0.9652 before D-224, ``ethyl_acetate`` at
    0.6688 after -- and its magnitude fell 31 %. ``fruity`` keeps its owner and falls 21 %
    (2.3292 -> 1.8404). Both verdicts are unchanged and correct for a sound ale: solventy stays
    below threshold, fruity stays above it, which is the banana-forward ale D-96 anchored for.

    Pinned because a seven-value re-anchoring that silently swapped a descriptor's owner would
    otherwise be visible nowhere -- the same class of blind spot as the levels themselves.
    """
    compiled = compile_scenario(_beer_calibration_scenario())
    # simulate() rather than compiled.run(): this scenario has no scheduled events, so the two are
    # the same integration, and `sensory_profile` takes a plain Trajectory.
    traj = simulate(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
        param_tiers=compiled.parameters.tier_map(),
        t_eval=np.array([0.0, _BEER_CALIBRATION_DAYS * 24.0]),
    )
    profile = sensory_profile(traj, load_parameters(default_data_dir() / "sensory.yaml"))
    descriptors = MaxRuleProjector().project(profile)

    solventy = descriptors.readings["solventy"]
    assert solventy.dominant == "ethyl_acetate", (
        f"beer's solventy axis is owned by {solventy.dominant}; D-224 handed it from "
        f"isoamyl_alcohol to ethyl_acetate and the hand-over is the output-level statement of "
        f"that beat."
    )
    assert solventy.magnitude == pytest.approx(0.669, abs=0.01)
    assert solventy.magnitude < 1.0, "a sound ale must not read a solventy note"

    fruity = descriptors.readings["fruity"]
    assert fruity.dominant == "isoamyl_acetate"
    # D-227: 1.8333 -> 1.8294, the -0.22 % that beer's isoamyl acetate moved (it is the one
    # ester reading its own precursor pool, so the k solve that holds the other two EXACTLY
    # cannot hold it). Re-pinned to the value the engine returns, not to the 1.840 D-224 wrote:
    # that pin was 0.36 % away from the 1.8333 it was pinning and only its 0.01 tolerance hid
    # the gap [[feedback-pin-tolerance-vs-solver-tolerance]].
    assert fruity.magnitude == pytest.approx(1.8294, abs=0.005)
    assert fruity.magnitude > 1.0, "an ale's banana note is above threshold (D-96)"


# =====================================================================================
# D-225: the aroma level guard, driven by the REGISTRY rather than by a list
# =====================================================================================
#: The frame each medium's aroma constants are calibrated in, EXCEPT the YAN, which is read
#: per-parameter out of the `conditions:` sentence. Wine has TWO frames and that is the open
#: question D-225 records rather than resolves: its three ester k reproduce their stated levels
#: at YAN 80 and its five Ehrlich k at YAN 250, so a single wine frame would put one group or
#: the other out. Beer's clauses state no YAN and its whole set shares one frame.
_AROMA_FRAMES = {
    "beer": {"days": 21.0, "celsius": 20.0, "default_yan": 200.0},
    "wine": {"days": 21.0, "celsius": 20.0, "default_yan": 250.0},
}
#: "k set to land finished <molecule> at [~|the ]<n> mg/L". The target is PARSED, never
#: transcribed -- a second copy in this file would be free to diverge from the spec (D-158).
_LANDING_CLAUSE = re.compile(
    r"set to land finished [a-z0-9 ,\-]*?at (?:~|the )?([0-9.]+) mg/L", re.IGNORECASE
)
#: The frame's YAN, where the clause names one. Only wine's three esters do, and for exactly one
#: of them it is worth 1.65x.
_FRAME_YAN = re.compile(r"\bat YAN ([0-9.]+)", re.IGNORECASE)


def _aroma_k_params() -> tuple[str, ...]:
    """Every aroma rate constant, ENUMERATED FROM THE REGISTRY rather than listed.

    This function is the whole point of D-225. D-224 guarded seven aroma constants against a
    hand-written dict of seven, and beer's eighth -- ``k_ethyl_hexanoate``, the third entry in
    ``ESTER_SPECS`` -- was not in the dict, so it alone kept the full drift that record is about
    (0.787x of its own stated level) while the guard suite looked complete. A literal list cannot
    fail to be complete; a registry read can. Add a ninth ester tomorrow and it either carries a
    landing clause and gets levelled, or the completeness test below goes red naming it.
    """
    # The two spec types are iterated separately rather than unpacked into one tuple: they are
    # different dataclasses, so mypy joins them to `object` and loses `.k_param`.
    return tuple([spec.k_param for spec in ESTER_SPECS] + [spec.k_param for spec in FUSEL_SPECS])


def _aroma_spec(medium: str, k_param: str) -> tuple[float, float]:
    """(target mg/L, frame YAN) read out of the parameter's own ``conditions:`` sentence."""
    params = load_parameters(default_data_dir() / f"{medium}_generic.yaml")
    conditions = params[k_param].provenance.conditions
    match = _LANDING_CLAUSE.search(conditions)
    assert match is not None, (
        f"{medium}'s {k_param} states no landing level in its `conditions:` field. Every aroma "
        f"rate constant in ESTER_SPECS/FUSEL_SPECS is defined by the level it lands, so the "
        f"level belongs in the sentence that specifies the parameter -- not in `notes:`, where "
        f"no test can read it. That is exactly how beer's k_ethyl_hexanoate hid through D-224."
    )
    yan_match = _FRAME_YAN.search(conditions)
    yan = float(yan_match.group(1)) if yan_match else _AROMA_FRAMES[medium]["default_yan"]
    return float(match.group(1)), yan


def _aroma_level_mgl(medium: str, pool: str, yan: float) -> float:
    frame = _AROMA_FRAMES[medium]
    initial: dict[str, float] = {"yan_mgl": yan}
    if medium == "beer":
        initial.update(
            {"glucose_gpl": 15.0, "maltose_gpl": 70.0, "maltotriose_gpl": 15.0, "pitch_gpl": 1.0}
        )
    else:
        initial.update({"brix": 24.0, "pitch_gpl": 0.5})
    scenario = Scenario(
        name=f"{medium}-aroma-registry-census",
        medium=medium,
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=frame["celsius"])],
        duration_days=frame["days"],
    )
    traj = compile_scenario(scenario).run(t_eval=np.array([0.0, frame["days"] * 24.0]))
    return float(traj.series(pool)[-1]) * 1000.0


@pytest.mark.parametrize("medium", ["beer", "wine"])
def test_every_aroma_rate_constant_in_the_registry_declares_its_landing_level(medium):
    """The COMPLETENESS half: no aroma constant may be defined by a level it does not state.

    D-224's level guard could only check the constants someone had remembered to list. This one
    cannot be incomplete, because it asks the registry what the population is. It is deliberately
    separate from the level assertion below: a constant that states no target fails HERE, with a
    message about provenance, rather than being silently skipped by a census that then reports
    itself clean on a denominator it never measured.
    """
    for k_param in _aroma_k_params():
        target, yan = _aroma_spec(medium, k_param)
        assert target > 0.0, f"{medium}/{k_param} declares a non-positive target"
        assert yan > 0.0


@pytest.mark.parametrize("medium", ["beer", "wine"])
def test_every_aroma_rate_constant_lands_the_level_it_declares(medium):
    """The LEVEL half, over the registry: all eight constants in each medium, at their own frame.

    Supersedes the coverage of the two hand-listed level guards above without replacing them --
    those carry D-224's reasoning in their failure messages and pin its seven literals exactly.
    What this adds is that the EIGHTH is checked, and that a ninth would be.

    The frame's YAN is read from each constant's own `conditions:` sentence rather than fixed per
    medium, because wine genuinely has two: its esters reproduce at YAN 80 and its Ehrlich
    constants at YAN 250. That inconsistency is real and D-225 records it as open rather than
    papering over it -- but it must be VISIBLE in the data, so the test reads it instead of
    hard-coding one frame and quietly failing the other group.
    """
    specs: dict[str, str] = {spec.k_param: spec.pool for spec in ESTER_SPECS}
    specs.update({spec.k_param: spec.pool for spec in FUSEL_SPECS})
    for k_param in _aroma_k_params():
        target, yan = _aroma_spec(medium, k_param)
        level = _aroma_level_mgl(medium, specs[k_param], yan)
        assert level == pytest.approx(target, rel=0.02), (
            f"{medium}'s {specs[k_param]} finishes at {level:.4f} mg/L against the {target} mg/L "
            f"its own `conditions:` field says {k_param} is set to land (frame YAN {yan:.0f}). "
            f"The rate is linear in k, so re-anchor it -- do NOT relax this tolerance. Read "
            f"D-224 and D-225 first: both exist because a level defined by an upstream speed "
            f"drifts whenever anything upstream moves, and nothing noticed for three beats."
        )


def test_beers_apple_note_reaches_its_threshold_without_moving_any_descriptor_axis():
    """D-225's output-level statement: the pool crosses 1.0, the axis it feeds does not move.

    Re-anchoring ``k_ethyl_hexanoate`` x1.27 takes that pool's own OAV from 0.825 to 1.048, i.e.
    across the detection line. The crossing is INTENDED -- ``beer_generic.yaml`` says this
    molecule sits right at Meilgaard's ~0.21 mg/L threshold, and the pre-D-225 0.825 fell short
    only because the level had drifted 21 % low. Scored rather than asserted, because D-224's own
    follow-up exists for shipping seven values without scoring the axis they feed.

    The axis does NOT move: under the D-95 MAX rule ``fruity`` reads its loudest contributor, and
    isoamyl acetate at 1.840 outranks ethyl hexanoate at 1.048, so both the magnitude and the
    owner are unchanged and ``above_threshold()`` is ``['fruity']`` either way. That is the honest
    description -- the crossing is real at the pool and masked at the axis -- and pinning BOTH
    halves is what stops a future re-anchoring from flipping the owner unseen.
    """
    compiled = compile_scenario(_beer_calibration_scenario())
    traj = simulate(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
        param_tiers=compiled.parameters.tier_map(),
        t_eval=np.array([0.0, _BEER_CALIBRATION_DAYS * 24.0]),
    )
    profile = sensory_profile(traj, load_parameters(default_data_dir() / "sensory.yaml"))
    descriptors = MaxRuleProjector().project(profile)

    # The pool crosses its own threshold, and that is the point of the re-anchoring.
    assert profile.readings["ethyl_hexanoate"].oav == pytest.approx(1.048, abs=0.01)
    assert profile.readings["ethyl_hexanoate"].oav > 1.0, (
        "beer_generic.yaml says ethyl hexanoate sits right at its ~0.21 mg/L threshold; an OAV "
        "below 1.0 means the level has drifted low again (it read 0.825 before D-225)."
    )
    # ...and the axis is unmoved, because isoamyl acetate still shouts louder.
    fruity = descriptors.readings["fruity"]
    assert fruity.dominant == "isoamyl_acetate", (
        f"fruity is now owned by {fruity.dominant}. D-225 measured that ethyl hexanoate crosses "
        f"its threshold while isoamyl acetate keeps the axis at 1.840; a change of owner here is "
        f"a different beer and needs its own record."
    )
    # D-227: 1.8333 -> 1.8294, the -0.22 % that beer's isoamyl acetate moved (it is the one
    # ester reading its own precursor pool, so the k solve that holds the other two EXACTLY
    # cannot hold it). Re-pinned to the value the engine returns, not to the 1.840 D-224 wrote:
    # that pin was 0.36 % away from the 1.8333 it was pinning and only its 0.01 tolerance hid
    # the gap [[feedback-pin-tolerance-vs-solver-tolerance]].
    assert fruity.magnitude == pytest.approx(1.8294, abs=0.005)
    assert descriptors.above_threshold() == ["fruity"]


def test_the_ethyl_hexanoate_band_spans_its_sourced_ale_range_computed_at_the_edges(tmp_path):
    """D-225's band rule, in the COMPUTED form D-224 named the stronger of the two.

    Owed to a mutation that came back GREEN: reverting this band to its pre-D-225 edges while
    leaving the nominal alone broke nothing, because the level guards read the nominal and the
    two hand-listed band-rule tests name only the five Ehrlich k and the two acetate esters. A
    band nothing asserts is a band free to drift, and this one had drifted -- its stated span
    "~0.05-0.6 mg/L" was a gloss that measured 0.0433-0.5052 at the drifted nominal.

    Asserting landed CONCENTRATIONS rather than a multiplier is what makes this the stronger
    form (D-224 Sec 8): a multiplier is invariant to a joint rescale of value and band, so it
    stays green through exactly the re-anchoring this file keeps needing. These edges are
    computed to span the molecule's own sourced ~0.1-0.5 mg/L ale range, so a drawn beer can
    never sit outside the range its source reports -- which is why the rescale rule was
    REFUSED here: rescaled, the top would have reached 0.6416 mg/L, 28 % above that range.
    """
    beer = load_parameters(default_data_dir() / "beer_generic.yaml")
    lo, hi = _band_edges(beer, "k_ethyl_hexanoate", 7.76937e-5, 3.88467e-4)

    def level_at(k_value: float) -> float:
        # NOT `_beer_aroma_levels_mgl`: that helper returns only the SEVEN pools D-224 listed by
        # hand, so it has no `ethyl_hexanoate` key -- the same incompleteness this beat is about,
        # met inside the test written to fix it. That dict is deliberately left at seven so
        # D-225's arm-A measurement (its guard stays GREEN on this defect) stays reproducible;
        # coverage of all eight lives in the registry-driven tests below.
        scenario = _beer_calibration_scenario()
        compiled = compile_scenario(
            scenario, data_dir=_beer_data_dir_at(tmp_path, k_ethyl_hexanoate=k_value)
        )
        traj = compiled.run(t_eval=np.array([0.0, _BEER_CALIBRATION_DAYS * 24.0]))
        return float(traj.series("ethyl_hexanoate")[-1]) * 1000.0

    low_level = level_at(lo)
    high_level = level_at(hi)
    assert low_level == pytest.approx(0.100, abs=0.002), (
        f"the band's low edge lands {low_level:.4f} mg/L, not the 0.100 that is the floor of this "
        f"molecule's sourced 0.1-0.5 ale range. The edges are COMPUTED to span that range; "
        f"re-measure them rather than re-pinning this number."
    )
    assert high_level == pytest.approx(0.500, abs=0.005), (
        f"the band's high edge lands {high_level:.4f} mg/L, not the 0.500 top of the sourced ale "
        f"range. Above it, a drawn beer reports a concentration its own source rejects."
    )

    # The sensory reading at each edge, pinned because the nominal crossing IS the point here and
    # a band that crosses needs its reach stated rather than claimed away (D-224 Sec 7).
    threshold = (
        load_parameters(default_data_dir() / "sensory.yaml")["threshold_ethyl_hexanoate_beer"].value
        / 1000.0
    )
    assert low_level / threshold == pytest.approx(0.476, abs=0.01)
    assert high_level / threshold == pytest.approx(2.381, abs=0.02)
    # The nominal is above threshold by design; the band's FLOOR is the only place an ale reads
    # no apple note at all, and that is a beer at the bottom of the molecule's reported range.
    assert low_level / threshold < 1.0
    assert high_level / threshold > 1.0
