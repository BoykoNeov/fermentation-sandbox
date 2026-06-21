"""Tests for EthanolInhibition — the multiplicative ethanol-inhibition modifier.

EthanolInhibition is a RateModifier, not a summed Process: it scales the
sugar-uptake flux by ``f = (1 - E/E_max)**n``. These tests pin the closed-form
factor and its guards, prove the factor scales uptake's whole contribution (and
that toggling the modifier off recovers the uninhibited rate), and — the
architectural point — that an inhibited run still conserves carbon and mass,
because scaling a conserving flux by one scalar preserves its balances (D-10).

We deliberately do NOT assert dryness-under-inhibition: with the placeholder
``ethanol_tolerance`` (110 g/L, below a 24 Brix must's final ethanol) an inhibited
run stalls short of dryness. That is a parameter-tuning concern, not a property of
the form; the benchmark stays skipped until sourcing lands.
"""

import numpy as np
import pytest

from fermentation.core.kinetics import EthanolInhibition, SugarUptakeToEthanolCO2
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
)


@pytest.fixture
def store():
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema, *, x: float = 2.0, s: float = 264.0, e: float = 0.0
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": e, "N": 0.0, "T": 293.15, "CO2": 0.0})


def test_metadata():
    inh = EthanolInhibition()
    assert inh.name == "ethanol_inhibition"
    assert inh.tier is Tier.PLAUSIBLE
    # It scales the uptake Process, referenced by that Process's own name.
    assert inh.modifies == (SugarUptakeToEthanolCO2.name,)
    assert set(inh.reads) == {"ethanol_tolerance", "ethanol_inhibition_exponent"}


def test_factor_matches_closed_form(params):
    # f(E) = (1 - E/E_max)**n at a known ethanol level.
    schema = wine_schema()
    e = 55.0
    y = _wine_y0(schema, e=e)
    f = EthanolInhibition().factor(0.0, y, schema, params)
    e_max = params["ethanol_tolerance"]
    n = params["ethanol_inhibition_exponent"]
    assert f == pytest.approx((1.0 - e / e_max) ** n)


def test_factor_is_one_without_ethanol(params):
    schema = wine_schema()
    f = EthanolInhibition().factor(0.0, _wine_y0(schema, e=0.0), schema, params)
    assert f == pytest.approx(1.0)


def test_factor_monotonically_decreasing_and_in_unit_interval(params):
    schema = wine_schema()
    inh = EthanolInhibition()
    e_max = params["ethanol_tolerance"]
    levels = np.linspace(0.0, e_max, 25)
    factors = [inh.factor(0.0, _wine_y0(schema, e=float(e)), schema, params) for e in levels]
    assert factors[0] == pytest.approx(1.0)
    assert all(0.0 <= f <= 1.0 for f in factors)
    # Strictly decreasing until it reaches the wall (pairwise; lengths differ).
    assert all(b < a for a, b in zip(factors, factors[1:], strict=False))
    assert factors[-1] == pytest.approx(0.0)


def test_factor_clamps_to_zero_past_tolerance(params):
    # Above the tolerance the wall term goes negative; the flux is fully shut down,
    # never negative (which would let uptake *create* sugar).
    schema = wine_schema()
    e_max = params["ethanol_tolerance"]
    f = EthanolInhibition().factor(0.0, _wine_y0(schema, e=e_max + 20.0), schema, params)
    assert f == 0.0


def test_negative_ethanol_excursion_does_not_amplify(params):
    # A solver undershoot below zero must not read as *less* inhibition (factor>1).
    schema = wine_schema()
    f = EthanolInhibition().factor(0.0, _wine_y0(schema, e=-1e-6), schema, params)
    assert f == pytest.approx(1.0)


def test_modifier_scales_whole_uptake_contribution(params):
    # In a ProcessSet, the modifier multiplies uptake's entire (dS, dE, dCO2)
    # vector by the same scalar factor.
    schema = wine_schema()
    y = _wine_y0(schema, x=2.0, s=200.0, e=44.0)
    uptake_only = ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True)
    inhibited = ProcessSet(
        schema, [SugarUptakeToEthanolCO2()], modifiers=[EthanolInhibition()], strict=True
    )
    d_raw = uptake_only.total_derivatives(0.0, y, params)
    d_inh = inhibited.total_derivatives(0.0, y, params)
    f = EthanolInhibition().factor(0.0, y, schema, params)
    assert 0.0 < f < 1.0  # ethanol present -> partial inhibition
    for var in ("S", "E", "CO2"):
        assert schema.get(d_inh, var) == pytest.approx(f * schema.get(d_raw, var))


def test_toggling_modifier_off_recovers_uninhibited_rate(params):
    schema = wine_schema()
    y = _wine_y0(schema, x=2.0, s=200.0, e=44.0)
    ps = ProcessSet(
        schema, [SugarUptakeToEthanolCO2()], modifiers=[EthanolInhibition()], strict=True
    )
    uptake_only = ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True)

    ps.disable("ethanol_inhibition")
    assert np.allclose(
        ps.total_derivatives(0.0, y, params),
        uptake_only.total_derivatives(0.0, y, params),
    )
    # And re-enabling restores the scaled-down rate.
    ps.enable("ethanol_inhibition")
    assert schema.get(ps.total_derivatives(0.0, y, params), "S") > schema.get(
        uptake_only.total_derivatives(0.0, y, params), "S"
    )  # less negative -> slower consumption


def test_inhibition_drags_uptake_output_tier_down():
    # A speculative tolerance/exponent makes the modifier's effect speculative, so
    # the uptake outputs it scales (S, E, CO2) report speculative even though both
    # uptake and a vanilla modifier are plausible. Use a synthetic speculative
    # modifier to isolate the tier mechanism from the (currently plausible) real one.
    schema = wine_schema()

    class SpeculativeInhibition(EthanolInhibition):
        name = "spec_inhibition"
        tier = Tier.SPECULATIVE

    ps = ProcessSet(
        schema,
        [SugarUptakeToEthanolCO2()],
        modifiers=[SpeculativeInhibition()],
    )
    assert ps.tier_of("S") is Tier.SPECULATIVE
    assert ps.tier_of("E") is Tier.SPECULATIVE
    assert ps.tier_of("CO2") is Tier.SPECULATIVE
    # Disabling the speculative modifier restores uptake's own plausible tier.
    ps.disable("spec_inhibition")
    assert ps.tier_of("S") is Tier.PLAUSIBLE
    assert ps.overall_tier() is Tier.PLAUSIBLE


def test_inhibited_run_conserves_carbon_and_mass(params):
    # The architectural guarantee: scaling a conserving flux preserves its
    # balances. An inhibited wine run still closes carbon and mass exactly.
    schema = wine_schema()
    ps = ProcessSet(
        schema, [SugarUptakeToEthanolCO2()], modifiers=[EthanolInhibition()], strict=True
    )
    traj = simulate(ps, params=params, y0=_wine_y0(schema), t_span=(0.0, 500.0))
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
    # Ethanol never exceeds the tolerance wall (the flux stalls as it is approached).
    assert float(traj.series("E")[-1]) <= params["ethanol_tolerance"] + 1e-6


def test_inhibition_slows_fermentation(params):
    # Directional, parameter-robust: inhibition leaves strictly more residual sugar
    # over the same span than the uninhibited flux. (We do not assert an absolute
    # stall threshold — that depends on the placeholder tolerance.)
    schema = wine_schema()
    y0 = _wine_y0(schema)
    uninhibited = simulate(
        ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True),
        params=params,
        y0=y0,
        t_span=(0.0, 500.0),
    )
    inhibited = simulate(
        ProcessSet(
            schema, [SugarUptakeToEthanolCO2()], modifiers=[EthanolInhibition()], strict=True
        ),
        params=params,
        y0=y0,
        t_span=(0.0, 500.0),
    )
    assert uninhibited.success and inhibited.success
    assert float(inhibited.series("S")[-1]) > float(uninhibited.series("S")[-1])
