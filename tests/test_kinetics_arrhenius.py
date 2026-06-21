"""Tests for ArrheniusTemperature — the multiplicative temperature-dependence modifier.

ArrheniusTemperature is a RateModifier (like EthanolInhibition, reusing the D-10
hook): it scales a Process's flux by ``f = exp(-(E_a/R)(1/T - 1/T_ref))``,
normalised so ``f = 1`` at the reference temperature. It is *parameterised* — one
instance per rate, each reading its own activation energy (``E_a_growth`` for
growth, ``E_a_uptake`` for uptake), sharing ``T_ref`` (decision D-11).

These tests pin the closed-form factor and its reference-anchoring, prove the
factor scales each targeted Process's whole contribution (and toggles off cleanly),
that a speculative activation energy drags the scaled output's tier down (D-1), and
— the architectural point — that conservation survives *stacked* modifiers (uptake
scaled by ethanol inhibition AND Arrhenius is still one scalar on a conserving
vector). The behavioural check that matters for isothermal M1 runs: a warmer run
ferments faster than a cooler one.
"""

import math

import numpy as np
import pytest

from fermentation.core.kinetics import (
    ArrheniusTemperature,
    EthanolInhibition,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
)
from fermentation.core.kinetics.arrhenius import GAS_CONSTANT
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_mass,
    total_nitrogen,
)


@pytest.fixture
def store():
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x: float = 2.0,
    s: float = 264.0,
    e: float = 0.0,
    n: float = 0.0,
    temp: float = 293.15,
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": temp, "CO2": 0.0})


def test_metadata():
    growth = ArrheniusTemperature.for_growth()
    uptake = ArrheniusTemperature.for_uptake()
    # The Arrhenius law is textbook, so the mechanism is plausible (placeholder
    # E_a values are speculative; tier propagation handles that separately).
    assert growth.tier is Tier.PLAUSIBLE and uptake.tier is Tier.PLAUSIBLE
    # Per-instance names so the two coexist in one ProcessSet's shared name space.
    assert growth.name == "arrhenius_growth" and uptake.name == "arrhenius_uptake"
    assert growth.name != uptake.name
    # Each scales exactly its rate, referenced by that Process's own name.
    assert growth.modifies == (GrowthNitrogenLimited.name,)
    assert uptake.modifies == (SugarUptakeToEthanolCO2.name,)
    # Per-rate activation energy; shared reference temperature.
    assert set(growth.reads) == {"E_a_growth", "T_ref"}
    assert set(uptake.reads) == {"E_a_uptake", "T_ref"}


def test_factor_is_one_at_reference_temperature(params):
    # The defining property of the T_ref-normalised form: the measured rate
    # constant is used unscaled at its calibration temperature.
    schema = wine_schema()
    y = _wine_y0(schema, temp=params["T_ref"])
    assert ArrheniusTemperature.for_growth().factor(0.0, y, schema, params) == pytest.approx(1.0)
    assert ArrheniusTemperature.for_uptake().factor(0.0, y, schema, params) == pytest.approx(1.0)


def test_factor_matches_closed_form(params):
    schema = wine_schema()
    temp = 303.15  # 30 C, above the 20 C reference
    y = _wine_y0(schema, temp=temp)
    f = ArrheniusTemperature.for_growth().factor(0.0, y, schema, params)
    expected = math.exp(
        -(params["E_a_growth"] / GAS_CONSTANT) * (1.0 / temp - 1.0 / params["T_ref"])
    )
    assert f == pytest.approx(expected)


def test_factor_above_ref_speeds_below_ref_slows(params):
    schema = wine_schema()
    arr = ArrheniusTemperature.for_uptake()
    t_ref = params["T_ref"]
    warm = arr.factor(0.0, _wine_y0(schema, temp=t_ref + 10.0), schema, params)
    cool = arr.factor(0.0, _wine_y0(schema, temp=t_ref - 10.0), schema, params)
    assert warm > 1.0  # hotter than reference -> faster
    assert cool < 1.0  # colder than reference -> slower
    # Always positive (exp): no regime where the factor could flip a flux sign.
    assert cool > 0.0


def test_factor_monotonically_increasing_in_temperature(params):
    schema = wine_schema()
    arr = ArrheniusTemperature.for_growth()
    temps = np.linspace(283.15, 308.15, 25)  # 10 C .. 35 C
    factors = [arr.factor(0.0, _wine_y0(schema, temp=float(tt)), schema, params) for tt in temps]
    assert all(b > a for a, b in zip(factors, factors[1:], strict=False))
    assert all(f > 0.0 for f in factors)


def test_q10_is_biologically_sane(params):
    # With the placeholder E_a ~50-60 kJ/mol the rate roughly doubles per 10 C
    # (Q10 ~ 2), the standard biological range — a guard that the placeholders and
    # the gas constant compose into sensible physics, not just a passing formula.
    schema = wine_schema()
    arr = ArrheniusTemperature.for_growth()
    base, hot = 293.15, 303.15
    q10 = arr.factor(0.0, _wine_y0(schema, temp=hot), schema, params) / arr.factor(
        0.0, _wine_y0(schema, temp=base), schema, params
    )
    assert 1.7 < q10 < 2.6


def test_modifier_scales_whole_growth_contribution(params):
    # In a ProcessSet the modifier multiplies growth's entire (dX, dS, dN) vector
    # by the same scalar factor.
    schema = wine_schema()
    y = _wine_y0(schema, x=1.0, s=200.0, n=0.3, temp=303.15)
    growth_only = ProcessSet(schema, [GrowthNitrogenLimited()], strict=True)
    scaled = ProcessSet(
        schema,
        [GrowthNitrogenLimited()],
        modifiers=[ArrheniusTemperature.for_growth()],
        strict=True,
    )
    d_raw = growth_only.total_derivatives(0.0, y, params)
    d_scaled = scaled.total_derivatives(0.0, y, params)
    f = ArrheniusTemperature.for_growth().factor(0.0, y, schema, params)
    assert f > 1.0  # 30 C is above the 20 C reference
    for var in ("X", "S", "N"):
        assert schema.get(d_scaled, var) == pytest.approx(f * schema.get(d_raw, var))


def test_modifier_scales_whole_uptake_contribution(params):
    schema = wine_schema()
    y = _wine_y0(schema, x=2.0, s=200.0, temp=288.15)  # below reference -> slower
    uptake_only = ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True)
    scaled = ProcessSet(
        schema,
        [SugarUptakeToEthanolCO2()],
        modifiers=[ArrheniusTemperature.for_uptake()],
        strict=True,
    )
    d_raw = uptake_only.total_derivatives(0.0, y, params)
    d_scaled = scaled.total_derivatives(0.0, y, params)
    f = ArrheniusTemperature.for_uptake().factor(0.0, y, schema, params)
    assert 0.0 < f < 1.0
    for var in ("S", "E", "CO2"):
        assert schema.get(d_scaled, var) == pytest.approx(f * schema.get(d_raw, var))


def test_toggling_modifier_off_recovers_unscaled_rate(params):
    schema = wine_schema()
    y = _wine_y0(schema, x=2.0, s=200.0, temp=303.15)
    uptake_only = ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True)
    ps = ProcessSet(
        schema,
        [SugarUptakeToEthanolCO2()],
        modifiers=[ArrheniusTemperature.for_uptake()],
        strict=True,
    )
    ps.disable("arrhenius_uptake")
    assert np.allclose(
        ps.total_derivatives(0.0, y, params),
        uptake_only.total_derivatives(0.0, y, params),
    )
    # Re-enabling restores the (here, sped-up) scaled rate: consumption is faster.
    ps.enable("arrhenius_uptake")
    assert schema.get(ps.total_derivatives(0.0, y, params), "S") < schema.get(
        uptake_only.total_derivatives(0.0, y, params), "S"
    )  # more negative -> faster consumption above the reference temperature


def test_speculative_activation_energy_drags_output_tier_down(store):
    # Real path (D-1): the modifier mechanism is plausible, but its placeholder
    # E_a/T_ref are speculative, so via param-tier propagation the rate it scales
    # reports speculative. Disabling the modifier restores uptake's own tier.
    schema = wine_schema()
    param_tiers = store.tier_map()
    ps = ProcessSet(
        schema,
        [SugarUptakeToEthanolCO2()],
        modifiers=[ArrheniusTemperature.for_uptake()],
    )
    for var in ("S", "E", "CO2"):
        assert ps.tier_of(var, param_tiers) is Tier.SPECULATIVE
    ps.disable("arrhenius_uptake")
    # Uptake's own reads (q_sugar_max etc.) are themselves speculative placeholders,
    # so the structural (no param_tiers) tier is the clean check that the *modifier*
    # was what we removed: it drops from plausible-modifier-capped back to the
    # Process's plausible tier.
    assert ps.tier_of("S") is Tier.PLAUSIBLE


def test_arrhenius_alone_conserves_carbon_and_mass(params):
    # Scaling uptake by a single positive factor preserves its balances: a
    # temperature-scaled wine run still closes carbon and mass exactly.
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [SugarUptakeToEthanolCO2()],
        modifiers=[ArrheniusTemperature.for_uptake()],
        strict=True,
    )
    traj = simulate(ps, params=params, y0=_wine_y0(schema, temp=298.15), t_span=(0.0, 500.0))
    assert traj.success
    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=0.488),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_conserved(traj, total_mass(schema), rtol=1e-5, atol=1e-6, label="mass")
    assert_nonnegative(traj, ("S", "E", "CO2"), atol=1e-7)


def test_full_model_with_stacked_modifiers_conserves_carbon_and_nitrogen(params, store):
    # The architectural guarantee under STACKING: growth + uptake + ethanol
    # inhibition + Arrhenius (both rates). Uptake's vector is multiplied by two
    # factors (inhibition AND Arrhenius); that is still one combined scalar on a
    # conserving vector, so carbon and nitrogen close exactly end-to-end (D-11).
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2()],
        modifiers=[
            EthanolInhibition(),
            ArrheniusTemperature.for_growth(),
            ArrheniusTemperature.for_uptake(),
        ],
        strict=True,
    )
    y0 = schema.pack({"X": 0.1, "S": [264.0], "E": 0.0, "N": 0.3, "T": 298.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 500.0))
    assert traj.success

    f_c = store.value("biomass_C_fraction")
    f_n = store.value("biomass_N_fraction")
    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=f_c),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_conserved(
        traj,
        total_nitrogen(schema, biomass_nitrogen_fraction=f_n),
        rtol=1e-5,
        atol=1e-6,
        label="nitrogen",
    )
    assert_nonnegative(traj, ("X", "S", "N", "E", "CO2"), atol=1e-7)


def test_warmer_ferments_faster_than_cooler(params):
    # The behavioural point of the modifier on isothermal M1 runs: identical
    # configs differing only in temperature ferment at different speeds. We use
    # uptake-only (fixed catalyst, no inhibition) over a span the cooler run cannot
    # finish, so the comparison is unconfounded by the placeholder ethanol wall.
    schema = wine_schema()
    arr = ProcessSet(
        schema,
        [SugarUptakeToEthanolCO2()],
        modifiers=[ArrheniusTemperature.for_uptake()],
        strict=True,
    )
    warm = simulate(arr, params=params, y0=_wine_y0(schema, temp=298.15), t_span=(0.0, 100.0))
    cool = simulate(arr, params=params, y0=_wine_y0(schema, temp=288.15), t_span=(0.0, 100.0))
    assert warm.success and cool.success
    assert float(warm.series("S")[-1]) < float(cool.series("S")[-1])
