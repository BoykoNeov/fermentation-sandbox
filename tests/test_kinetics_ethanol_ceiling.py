"""Tests for EthanolToleranceDeath — the super-linear ethanol ceiling (decision D-129).

Coleman's ``k_d = k'_d·E`` (:class:`EthanolInactivation`) is linear and unbounded, so
past its 265-300 g/L validated envelope a high-sugar must over-ferments to an impossible
ABV instead of sticking. :class:`EthanolToleranceDeath` adds ``Φ(E) = k_d2·max(E − E_tol, 0)²``
extra death *within the same X→X_dead mechanism* (D-13-compatible, not an uptake wall).

These tests pin two things that make it honest: (1) it is **exactly zero** below the strain
tolerance, so an in-envelope ferment is byte-for-byte the Coleman-only core (structural
isolability, prime directive #3); and (2) once ethanol overshoots tolerance it bites
super-linearly, sticking the must while still closing carbon/nitrogen (a within-biomass
transfer, like :class:`EthanolInactivation`).
"""

import pytest

from fermentation.core.kinetics import (
    EthanolInactivation,
    EthanolToleranceDeath,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
)
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
    total_nitrogen,
)


@pytest.fixture
def store():
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema, *, x: float = 2.0, e: float = 150.0, xd: float = 0.0
) -> FloatArray:
    return schema.pack(
        {"X": x, "S": [50.0], "E": e, "N": 0.0, "T": 293.15, "CO2": 0.0, "X_dead": xd}
    )


def test_metadata():
    death = EthanolToleranceDeath()
    assert death.name == "ethanol_tolerance_death"
    assert death.tier is Tier.SPECULATIVE
    assert set(death.touches) == {"X", "X_dead"}
    assert set(death.reads) == {"k_d2_ethanol_tolerance_death", "ethanol_tolerance"}


def test_derivative_matches_closed_form(params):
    # Φ = k_d2·(E - E_tol)² ; r = Φ·X ; dX = -r, dX_dead = +r  (E above tolerance).
    schema = wine_schema()
    x, e = 2.0, 160.0
    tol = params["ethanol_tolerance"]
    k2 = params["k_d2_ethanol_tolerance_death"]
    y = _wine_y0(schema, x=x, e=e)
    d = EthanolToleranceDeath().derivatives(0.0, y, schema, params)
    r = k2 * (e - tol) ** 2 * x
    assert schema.get(d, "X") == pytest.approx(-r)
    assert schema.get(d, "X_dead") == pytest.approx(+r)
    for var in ("S", "E", "N", "T", "CO2"):
        assert schema.get(d, var) == pytest.approx(0.0)


def test_exactly_inert_at_and_below_tolerance(params):
    # THE isolability property: identically zero for E ≤ E_tol, so an in-envelope
    # ferment is byte-for-byte the Coleman-only core. Checked just below, exactly at,
    # and far below tolerance.
    schema = wine_schema()
    tol = params["ethanol_tolerance"]
    death = EthanolToleranceDeath()
    for e in (0.0, 80.0, tol - 1e-9, tol):
        d = death.derivatives(0.0, _wine_y0(schema, e=e), schema, params)
        assert schema.get(d, "X") == 0.0
        assert schema.get(d, "X_dead") == 0.0


def test_transfer_is_biomass_neutral(params):
    # Viable loss equals inactivated gain exactly — total biomass derivative is zero.
    schema = wine_schema()
    d = EthanolToleranceDeath().derivatives(0.0, _wine_y0(schema, e=160.0), schema, params)
    assert schema.get(d, "X") + schema.get(d, "X_dead") == pytest.approx(0.0)


def test_rate_is_quadratic_in_overshoot(params):
    # max(E - E_tol, 0)²: doubling the overshoot quadruples the death rate.
    schema = wine_schema()
    tol = params["ethanol_tolerance"]
    death = EthanolToleranceDeath()
    r1 = death.derivatives(0.0, _wine_y0(schema, e=tol + 10.0), schema, params)
    r2 = death.derivatives(0.0, _wine_y0(schema, e=tol + 20.0), schema, params)
    assert schema.get(r2, "X_dead") == pytest.approx(4.0 * schema.get(r1, "X_dead"))


def test_negative_excursions_do_not_resurrect(params):
    # A solver undershoot (E<0 or X<0) must not flip the sign and create/resurrect cells.
    schema = wine_schema()
    death = EthanolToleranceDeath()
    d_negE = death.derivatives(0.0, _wine_y0(schema, e=-1e-6), schema, params)
    d_negX = death.derivatives(0.0, _wine_y0(schema, x=-1e-6, e=160.0), schema, params)
    assert schema.get(d_negE, "X_dead") == 0.0
    assert schema.get(d_negX, "X_dead") == 0.0


def test_speculative_k_d2_drags_biomass_tier():
    # Parameter-tier propagation (D-1): the speculative ceiling makes the biomass pools
    # it writes report speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [EthanolToleranceDeath()])
    assert ps.tier_of("X") is Tier.SPECULATIVE
    assert ps.tier_of("X_dead") is Tier.SPECULATIVE


def test_in_envelope_run_is_byte_for_byte(params):
    # Numerical isolability: a must that dries below tolerance is identical with and
    # without the ceiling (the term is exactly zero the whole way).
    schema = wine_schema()
    y0 = schema.pack(
        {"X": 0.1, "S": [200.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    core = [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2(), EthanolInactivation()]
    t_core = simulate(ProcessSet(schema, core), params=params, y0=y0, t_span=(0.0, 800.0))
    t_ceil = simulate(
        ProcessSet(schema, [*core, EthanolToleranceDeath()]),
        params=params,
        y0=y0,
        t_span=(0.0, 800.0),
    )
    assert t_core.success and t_ceil.success
    # Finished below the 142 g/L tolerance, so the ceiling never bit.
    assert float(t_ceil.series("E")[-1]) < params["ethanol_tolerance"]
    assert float(t_ceil.series("E")[-1]) == pytest.approx(float(t_core.series("E")[-1]), rel=1e-9)
    assert float(t_ceil.series("S")[-1]) == pytest.approx(float(t_core.series("S")[-1]), abs=1e-9)


def test_high_sugar_sticks_and_conserves(params, store):
    # Past the envelope the ceiling arrests the ferment (residual sugar remains) while
    # still closing carbon AND nitrogen — the transfer is within-biomass, both pools at
    # identical f_C/f_N. Without the ceiling this must would over-ferment near-dry.
    schema = wine_schema()
    y0 = schema.pack(
        {"X": 0.1, "S": [340.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    core = [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2(), EthanolInactivation()]
    dry = simulate(ProcessSet(schema, core), params=params, y0=y0, t_span=(0.0, 2000.0))
    stuck = simulate(
        ProcessSet(schema, [*core, EthanolToleranceDeath()]),
        params=params,
        y0=y0,
        t_span=(0.0, 2000.0),
    )
    assert dry.success and stuck.success

    # The ceiling leaves substantial residual sugar the Coleman-only core burns through.
    assert float(stuck.series("S")[-1]) > 20.0
    assert float(stuck.series("S")[-1]) > float(dry.series("S")[-1]) + 20.0
    # Ethanol arrested above tolerance (it bit) but well short of the runaway dry ABV.
    assert float(stuck.series("E")[-1]) > params["ethanol_tolerance"]
    assert float(stuck.series("E")[-1]) < float(dry.series("E")[-1])

    f_c = store.value("biomass_C_fraction")
    f_n = store.value("biomass_N_fraction")
    assert_conserved(
        stuck,
        total_carbon(schema, biomass_carbon_fraction=f_c),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_conserved(
        stuck,
        total_nitrogen(schema, biomass_nitrogen_fraction=f_n),
        rtol=1e-5,
        atol=1e-6,
        label="nitrogen",
    )
    assert_nonnegative(stuck, ("X", "S", "N", "E", "CO2", "X_dead"), atol=1e-7)
