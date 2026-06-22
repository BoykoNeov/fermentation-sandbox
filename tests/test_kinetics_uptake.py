"""Tests for SugarUptakeToEthanolCO2 — the fermentative flux of the core.

Covers the closed-form derivative, the no-biomass shutoff, beer's sequential
(catabolite-repressed) uptake ordering, and conservation: carbon for wine and
beer, mass for wine (it does not close for beer — hydrolysis water — by D-8), and
carbon for a combined growth+uptake run (where mass is *not* expected to close,
since growth diverts sugar carbon into biomass).
"""

import numpy as np
import pytest

from fermentation.core.chemistry import co2_yield, ethanol_yield
from fermentation.core.kinetics import GrowthNitrogenLimited, SugarUptakeToEthanolCO2
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
    total_mass,
    total_nitrogen,
)


@pytest.fixture
def store():
    # Real (sourced) wine parameters; the kinetics are medium-agnostic, so the
    # wine file suffices to exercise the Process mechanism.
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema, *, x: float = 2.0, s: float = 264.0, e: float = 0.0
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": e, "N": 0.0, "T": 293.15, "CO2": 0.0, "X_dead": 0.0})


def _beer_y0(
    schema: StateSchema, *, x: float = 2.0, s: tuple[float, float, float] = (30.0, 60.0, 10.0)
) -> FloatArray:
    return schema.pack(
        {"X": x, "S": list(s), "E": 0.0, "N": 0.0, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )


def test_metadata():
    u = SugarUptakeToEthanolCO2()
    assert u.name == "sugar_uptake_to_ethanol_co2"
    assert u.tier is Tier.PLAUSIBLE
    assert set(u.touches) == {"S", "E", "CO2"}
    assert set(u.reads) == {"q_sugar_max", "K_sugar_uptake", "K_repression"}


def test_derivative_matches_closed_form(params):
    # Pin dS, dE, dCO2 to the formula at a known wine state — no solver fuzz.
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = schema.pack({"X": x, "S": [s], "E": 0.0, "N": 0.0, "T": 293.15, "CO2": 0.0, "X_dead": 0.0})
    d = SugarUptakeToEthanolCO2().derivatives(0.0, y, schema, params)

    r = params["q_sugar_max"] * x * (s / (params["K_sugar_uptake"] + s))
    assert schema.get(d, "S") == pytest.approx(-r)
    assert schema.get(d, "E") == pytest.approx(ethanol_yield("glucose") * r)
    assert schema.get(d, "CO2") == pytest.approx(co2_yield("glucose") * r)
    # Uptake touches S, E, CO2 only — never biomass or nitrogen.
    assert schema.get(d, "X") == 0.0
    assert schema.get(d, "N") == 0.0


def test_no_uptake_without_biomass(params):
    # Fermentation is biomass-catalysed: no yeast, no flux.
    schema = wine_schema()
    d = SugarUptakeToEthanolCO2().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)
    assert np.array_equal(d, schema.zeros())


def test_negative_sugar_excursion_does_not_create_sugar(params):
    # A clamp guards against solver undershoot below zero: rate must be >= 0, never
    # adding sugar back or driving E/CO2 negative.
    schema = wine_schema()
    y = schema.pack(
        {"X": 2.0, "S": [-1e-6], "E": 0.0, "N": 0.0, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    d = SugarUptakeToEthanolCO2().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_wine_run_conserves_carbon_and_mass(params):
    schema = wine_schema()
    ps = ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True)
    traj = simulate(ps, params=params, y0=_wine_y0(schema), t_span=(0.0, 500.0))
    assert traj.success

    # No biomass term active on S/E/CO2 here, so both balances close exactly.
    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=0.488),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_conserved(traj, total_mass(schema), rtol=1e-5, atol=1e-6, label="mass")
    assert_nonnegative(traj, ("S", "E", "CO2"), atol=1e-7)
    # Biomass-catalysed uptake runs the must to dryness without needing growth.
    assert traj.series("S")[-1] < 1.0


def test_beer_run_conserves_carbon(params):
    # Beer's di-/trisaccharides take up hydrolysis water, so mass does not close
    # (total_mass even rejects a multi-slot S); carbon is the cross-medium invariant.
    schema = beer_schema()
    ps = ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True)
    traj = simulate(ps, params=params, y0=_beer_y0(schema), t_span=(0.0, 500.0))
    assert traj.success

    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=0.488),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_nonnegative(traj, ("S", "E", "CO2"), atol=1e-7)


def test_beer_sequential_uptake_ordering(params):
    # Glucose > maltose > maltotriose: while glucose is plentiful, maltose (and
    # especially maltotriose) uptake is strongly repressed; once glucose is gone,
    # maltose proceeds. Evaluate the derivative directly at two states.
    schema = beer_schema()
    u = SugarUptakeToEthanolCO2()

    with_glucose = schema.pack(
        {
            "X": 2.0, "S": [20.0, 60.0, 10.0], "E": 0.0, "N": 0.0,
            "T": 293.15, "CO2": 0.0, "X_dead": 0.0,
        }
    )
    d1 = u.derivatives(0.0, with_glucose, schema, params)
    rates1 = -d1[schema.slice("S")]  # consumption rates, slot order
    r_glu, r_mal, r_tri = (float(v) for v in rates1)
    assert r_glu > 0.0
    # Repressed sugars trickle: maltose < 20% of glucose, maltotriose below maltose.
    assert r_mal < 0.2 * r_glu
    assert r_tri < r_mal

    no_glucose = schema.pack(
        {
            "X": 2.0, "S": [0.0, 60.0, 10.0], "E": 0.0, "N": 0.0,
            "T": 293.15, "CO2": 0.0, "X_dead": 0.0,
        }
    )
    d2 = u.derivatives(0.0, no_glucose, schema, params)
    r_mal_after = float(-d2[schema.slice("S")][1])
    # With glucose gone, maltose de-represses and accelerates sharply.
    assert r_mal_after > r_mal
    assert r_mal_after > 0.5 * r_glu


def test_beer_sequential_uptake_trajectory(params):
    # Emergent behaviour over a full run: the three wort sugars deplete in
    # preference order (glucose, then maltose, then maltotriose), not together --
    # the "sequential uptake" the task calls for, not just the snapshot mechanism.
    schema = beer_schema()
    ps = ProcessSet(schema, [SugarUptakeToEthanolCO2()], strict=True)
    s0 = (30.0, 60.0, 10.0)
    t_eval = np.linspace(0.0, 500.0, 1001)
    traj = simulate(
        ps, params=params, y0=_beer_y0(schema, s=s0), t_span=(0.0, 500.0), t_eval=t_eval
    )
    assert traj.success
    glucose, maltose, maltotriose = (traj.series("S")[i] for i in range(3))

    def first_time_below(series: FloatArray, frac: float, init: float) -> float:
        idx = int(np.argmax(series < frac * init))  # first True (or 0 if never)
        assert series[idx] < frac * init, "threshold never crossed"
        return float(traj.t[idx])

    # When glucose is essentially gone, maltose and maltotriose are still largely
    # intact: consumption is staggered, not simultaneous.
    g_idx = int(np.argmax(glucose < 0.05 * s0[0]))
    assert glucose[g_idx] < 0.05 * s0[0]
    assert maltose[g_idx] > 0.70 * s0[1]
    assert maltotriose[g_idx] > 0.70 * s0[2]

    # Depletion order: glucose before maltose before maltotriose.
    assert (
        first_time_below(glucose, 0.01, s0[0])
        < first_time_below(maltose, 0.01, s0[1])
        < first_time_below(maltotriose, 0.01, s0[2])
    )
    # And it all ferments out.
    assert float(traj.series("S")[:, -1].sum()) < 1.0


def test_combined_growth_and_uptake_conserves_carbon_and_nitrogen(params, store):
    # The two M1 processes run together on a wine must: growth draws sugar carbon
    # into biomass, uptake converts the rest to ethanol + CO2. Carbon and nitrogen
    # must still close. Mass is NOT expected to close — growth diverts sugar into
    # biomass, which total_mass{S,E,CO2} does not see (decision D-8).
    schema = wine_schema()
    ps = ProcessSet(schema, [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2()], strict=True)
    y0 = schema.pack(
        {"X": 0.1, "S": [264.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
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
