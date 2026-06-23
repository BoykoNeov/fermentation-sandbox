"""Tests for GrowthNitrogenLimited — the first validated-core kinetic Process.

Covers the closed-form derivative, the Monod shutoffs (no N / no sugar), and the
two atom balances the Process is built to conserve: nitrogen (free YAN + biomass
N) and carbon (sugar + biomass C). The carbon balance is exercised on beer too,
since that is the only path through the per-slot vector carbon draw.
"""

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics import GrowthNitrogenLimited
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)


@pytest.fixture
def store():
    # Real (sourced) wine parameters; the kinetics are medium-agnostic, so the
    # wine file suffices to exercise the Process mechanism (not beer-specific values).
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema, *, x: float = 0.1, s: float = 264.0, n: float = 0.3
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": 0.0, "N": n, "T": 293.15, "CO2": 0.0, "X_dead": 0.0})


def test_metadata():
    g = GrowthNitrogenLimited()
    assert g.name == "growth_nitrogen_limited"
    assert g.tier is Tier.PLAUSIBLE
    assert set(g.touches) == {"X", "S", "N"}
    # `reads` documents the params consumed (parameter-tier-propagation task).
    assert set(g.reads) == {"mu_max", "K_s", "K_n", "biomass_N_fraction", "biomass_C_fraction"}


def test_derivative_matches_closed_form(params):
    # Pin dX, dN, dS to the formula at a known wine state — no solver fuzz.
    schema = wine_schema()
    y = schema.pack(
        {"X": 1.0, "S": [200.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    d = GrowthNitrogenLimited().derivatives(0.0, y, schema, params)

    mu = params["mu_max"] * (200.0 / (params["K_s"] + 200.0)) * (0.3 / (params["K_n"] + 0.3))
    dx = mu * 1.0
    assert schema.get(d, "X") == pytest.approx(dx)
    assert schema.get(d, "N") == pytest.approx(-params["biomass_N_fraction"] * dx)
    # Carbon drawn from sugar closes exactly: c_frac_glucose * |dS| == f_C * dX.
    ds = schema.get(d, "S")
    assert ds < 0.0
    assert carbon_mass_fraction("glucose") * (-ds) == pytest.approx(
        params["biomass_C_fraction"] * dx
    )
    # Growth touches X, S, N only — never E or CO2.
    assert schema.get(d, "E") == 0.0
    assert schema.get(d, "CO2") == 0.0


def test_no_growth_without_nitrogen(params):
    schema = wine_schema()
    d = GrowthNitrogenLimited().derivatives(0.0, _wine_y0(schema, n=0.0), schema, params)
    assert np.array_equal(d, schema.zeros())


def test_no_growth_without_sugar(params):
    schema = wine_schema()
    d = GrowthNitrogenLimited().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)
    assert np.array_equal(d, schema.zeros())


def test_wine_run_conserves_carbon_and_nitrogen(params, store):
    schema = wine_schema()
    # strict=True enforces the touch contract on every solver step.
    ps = ProcessSet(schema, [GrowthNitrogenLimited()], strict=True)
    traj = simulate(ps, params=params, y0=_wine_y0(schema), t_span=(0.0, 300.0))
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
    assert_nonnegative(traj, ("X", "S", "N"), atol=1e-7)


def test_biomass_grows_then_caps_when_nitrogen_exhausted(params, store):
    schema = wine_schema()
    ps = ProcessSet(schema, [GrowthNitrogenLimited()], strict=True)
    x0, n0 = 0.1, 0.3
    traj = simulate(ps, params=params, y0=_wine_y0(schema, x=x0, n=n0), t_span=(0.0, 500.0))

    x = traj.series("X")
    n = traj.series("N")
    assert x[-1] > x[0]  # biomass grew (exponential phase)
    assert n[-1] < 1e-3  # YAN essentially exhausted (stationary phase)
    # Biomass is nitrogen-capped: N + f_N*X is conserved, so as N -> 0 the
    # biomass plateaus at X0 + N0/f_N. This IS the "stop dividing" mechanism.
    f_n = store.value("biomass_N_fraction")
    assert x[-1] == pytest.approx(x0 + n0 / f_n, rel=1e-3)


def test_beer_run_conserves_carbon_through_vector_draw(params, store):
    # Beer's three-slot sugar is the only path exercising the per-slot carbon
    # draw; carbon must still close to solver tolerance.
    schema = beer_schema()
    ps = ProcessSet(schema, [GrowthNitrogenLimited()], strict=True)
    y0 = schema.pack(
        {
            "X": 0.1,
            "S": [30.0, 60.0, 10.0],
            "E": 0.0,
            "N": 0.25,
            "T": 293.15,
            "CO2": 0.0,
            "X_dead": 0.0,
        }
    )
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 300.0))
    assert traj.success

    f_c = store.value("biomass_C_fraction")
    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=f_c),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_nonnegative(traj, ("X", "S", "N"), atol=1e-7)
