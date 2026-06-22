"""Tests for EthanolInactivation — the cumulative viability brake (decision D-13).

EthanolInactivation moves *viable* biomass ``X`` into the inactivated pool
``X_dead`` at the specific rate ``k_d = k'_d·E`` (Coleman 2007 eqs 2/7). These
tests pin the closed-form transfer and its guards, and prove the architectural
point: because the two biomass pools share one elemental composition, the
transfer is carbon- and nitrogen-neutral *by construction* — a full
growth+uptake+inactivation wine run still closes both atom balances exactly even
as cells die. That mass-neutrality is what lets the cumulative brake finish a
ferment where the reversible Luong wall stalls (see inhibition tests).
"""

import pytest

from fermentation.core.kinetics import (
    EthanolInactivation,
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
    schema: StateSchema, *, x: float = 2.0, e: float = 80.0, xd: float = 0.0
) -> FloatArray:
    return schema.pack(
        {"X": x, "S": [50.0], "E": e, "N": 0.0, "T": 293.15, "CO2": 0.0, "X_dead": xd}
    )


def test_metadata():
    inact = EthanolInactivation()
    assert inact.name == "ethanol_inactivation"
    assert inact.tier is Tier.PLAUSIBLE
    # Viable biomass leaves X and arrives in X_dead — both inside the contract.
    assert set(inact.touches) == {"X", "X_dead"}
    assert set(inact.reads) == {"k_prime_d"}


def test_derivative_matches_closed_form(params):
    # r = k'_d * E * X ; dX = -r, dX_dead = +r.
    schema = wine_schema()
    x, e = 2.0, 80.0
    y = _wine_y0(schema, x=x, e=e)
    d = EthanolInactivation().derivatives(0.0, y, schema, params)
    r = params["k_prime_d"] * e * x
    assert schema.get(d, "X") == pytest.approx(-r)
    assert schema.get(d, "X_dead") == pytest.approx(+r)
    # Nothing else moves.
    for var in ("S", "E", "N", "T", "CO2"):
        assert schema.get(d, var) == pytest.approx(0.0)


def test_transfer_is_biomass_neutral(params):
    # The defining property: viable loss equals inactivated gain, exactly. Their
    # sum (total biomass) has zero derivative from this Process alone.
    schema = wine_schema()
    d = EthanolInactivation().derivatives(0.0, _wine_y0(schema), schema, params)
    assert schema.get(d, "X") + schema.get(d, "X_dead") == pytest.approx(0.0)


def test_no_inactivation_without_ethanol(params):
    schema = wine_schema()
    d = EthanolInactivation().derivatives(0.0, _wine_y0(schema, e=0.0), schema, params)
    assert schema.get(d, "X") == 0.0
    assert schema.get(d, "X_dead") == 0.0


def test_no_inactivation_without_cells(params):
    schema = wine_schema()
    d = EthanolInactivation().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)
    assert schema.get(d, "X") == 0.0
    assert schema.get(d, "X_dead") == 0.0


def test_negative_excursions_do_not_resurrect(params):
    # A solver undershoot (E<0 or X<0) must not flip the sign and *create* viable
    # biomass / resurrect dead cells. Both clamp to zero rate.
    schema = wine_schema()
    d_negE = EthanolInactivation().derivatives(0.0, _wine_y0(schema, e=-1e-6), schema, params)
    d_negX = EthanolInactivation().derivatives(0.0, _wine_y0(schema, x=-1e-6), schema, params)
    assert schema.get(d_negE, "X_dead") == 0.0
    assert schema.get(d_negX, "X_dead") == 0.0


def test_rate_scales_with_ethanol(params):
    # k_d = k'_d·E is linear in ethanol: double the ethanol, double the rate.
    schema = wine_schema()
    inact = EthanolInactivation()
    r1 = inact.derivatives(0.0, _wine_y0(schema, e=40.0), schema, params)
    r2 = inact.derivatives(0.0, _wine_y0(schema, e=80.0), schema, params)
    assert schema.get(r2, "X_dead") == pytest.approx(2.0 * schema.get(r1, "X_dead"))


def test_strict_touches_contract(params):
    # Builds and runs under the strict touches check (writes only X, X_dead).
    schema = wine_schema()
    ps = ProcessSet(schema, [EthanolInactivation()], strict=True)
    y0 = _wine_y0(schema, x=2.0, e=80.0)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 200.0))
    assert traj.success
    # Viable biomass strictly declined; the inactivated pool absorbed the loss.
    assert float(traj.series("X")[-1]) < float(traj.series("X")[0])
    assert float(traj.series("X_dead")[-1]) > 0.0
    # Total biomass is conserved by this Process alone (pure transfer).
    total0 = float(traj.series("X")[0]) + float(traj.series("X_dead")[0])
    total1 = float(traj.series("X")[-1]) + float(traj.series("X_dead")[-1])
    assert total1 == pytest.approx(total0, rel=1e-5)


def test_full_model_conserves_carbon_and_nitrogen(params, store):
    # The architectural guarantee: growth + uptake + ethanol inactivation, run to
    # near-dryness, still closes carbon AND nitrogen exactly. Death is a within-
    # biomass transfer (X -> X_dead) and both pools carry identical f_C/f_N, so the
    # conservation checks (which weight both pools) are untouched by inactivation.
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2(), EthanolInactivation()],
        strict=True,
    )
    y0 = schema.pack(
        {"X": 0.1, "S": [264.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 600.0))
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
    assert_nonnegative(traj, ("X", "S", "N", "E", "CO2", "X_dead"), atol=1e-7)
    # Inactivation actually bit: a dead pool accumulated and viable biomass fell off
    # its peak by the end of the run.
    assert float(traj.series("X_dead")[-1]) > 0.5
    assert float(traj.series("X")[-1]) < float(traj.series("X").max())


def test_speculative_kprime_d_drags_biomass_tier():
    # Parameter-tier propagation (D-1): a speculative inactivation makes the biomass
    # pools it writes report speculative. Isolated with a synthetic speculative
    # subclass so the mechanism is tested independently of the real (plausible) tier.
    schema = wine_schema()

    class SpeculativeInactivation(EthanolInactivation):
        name = "spec_inactivation"
        tier = Tier.SPECULATIVE

    plausible = ProcessSet(schema, [EthanolInactivation()])
    assert plausible.tier_of("X") is Tier.PLAUSIBLE
    assert plausible.tier_of("X_dead") is Tier.PLAUSIBLE

    speculative = ProcessSet(schema, [SpeculativeInactivation()])
    assert speculative.tier_of("X") is Tier.SPECULATIVE
    assert speculative.tier_of("X_dead") is Tier.SPECULATIVE
