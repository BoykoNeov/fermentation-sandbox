"""Tests for the vicinal-diketone (VDK / diacetyl) pathway Processes (decision D-26).

The diacetyl beat models the *real* three-step pathway — α-acetolactate excretion →
spontaneous decarboxylation → yeast reduction — so the produce-then-reabsorb "diacetyl
rest" emerges rather than being scripted. This module grows one Process at a time:

* **Excretion** (:class:`AcetolactateExcretion`, this commit): fills the α-acetolactate
  reservoir from the fermentative flux, routing its C5 carbon out of ``S`` (option a1,
  D-19). Tested here for the closed-form derivative + exact sugar draw, temperature-
  flatness (a documented v1 simplification — the reservoir size is a weak lever), flux
  coupling, guards, per-RHS carbon closure and tier propagation.

Decarboxylation and reduction (and the emergent diacetyl-rest acceptance test) land in
the following commits.
"""

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics import (
    AcetolactateExcretion,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
)
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon

#: Carbon fractions the pools book against (mirror the Process/chemistry constants).
_ACETOLACTATE_C = carbon_mass_fraction("alpha_acetolactate")
_GLUCOSE_C = carbon_mass_fraction("glucose")


@pytest.fixture
def store():
    # Wine kinetics PLUS the shared VDK pathway constants (k_acetolactate, k_decarb, …).
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "vicinal_diketones.yaml",
    )


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


def test_excretion_metadata():
    p = AcetolactateExcretion()
    assert p.name == "acetolactate_excretion"
    # Speculative: rate magnitude is an order-of-magnitude estimate.
    assert p.tier is Tier.SPECULATIVE
    # Touches its own reservoir pool AND S — the C5 carbon is routed from sugar (a1, D-19);
    # never E/CO2 (the reservoir is the precursor, not a fermentation product yet).
    assert set(p.touches) == {"acetolactate", "S"}
    assert set(p.reads) == {"k_acetolactate", "K_sugar_uptake"}


# -- closed form & guards -----------------------------------------------------


def test_excretion_derivative_matches_closed_form(params):
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = _wine_y0(schema, x=x, s=s)
    d = AcetolactateExcretion().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    rate = params["k_acetolactate"] * flux
    assert schema.get(d, "acetolactate") == pytest.approx(rate)
    # The C5 carbon is routed from sugar (a1): one slot, so dS removes exactly the
    # acetolactate carbon converted back to grams of glucose. Carbon balances per-RHS.
    assert schema.get(d, "S") == pytest.approx(-rate * _ACETOLACTATE_C / _GLUCOSE_C)
    assert schema.get(d, "S") * _GLUCOSE_C == pytest.approx(
        -schema.get(d, "acetolactate") * _ACETOLACTATE_C
    )
    # Nothing else moves — the reservoir is a precursor; diacetyl/CO2 come later.
    for var in ("X", "E", "N", "CO2", "diacetyl", "butanediol"):
        assert schema.get(d, var) == 0.0


def test_excretion_is_temperature_flat(params):
    # Documented v1 simplification (D-26): excretion carries NO Arrhenius factor — the
    # reservoir size is a weak lever on the rest; temperature-criticality lives in the
    # decarboxylation. So the excretion rate is identical cold vs warm at equal flux.
    schema = wine_schema()
    cold = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(cold, "acetolactate") == pytest.approx(schema.get(warm, "acetolactate"))
    assert schema.get(cold, "acetolactate") > 0.0


def test_excretion_scales_with_fermentative_flux(params):
    # Coupled to the biomass-catalysed sugar flux (linear in X): twice the biomass ⇒
    # twice the excretion.
    schema = wine_schema()
    r1 = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "acetolactate") == pytest.approx(2.0 * schema.get(r1, "acetolactate"))


def test_excretion_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    assert AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)[
        schema.slice("acetolactate")
    ] == pytest.approx(0.0)
    assert AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)[
        schema.slice("acetolactate")
    ] == pytest.approx(0.0)


def test_excretion_negative_excursion_does_not_produce(params):
    # A solver undershoot (S<0 or X<0) must not flip the clamp and create acetolactate.
    schema = wine_schema()
    assert np.array_equal(
        AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, s=-1e-6), schema, params),
        schema.zeros(),
    )
    assert np.array_equal(
        AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=-1e-6), schema, params),
        schema.zeros(),
    )


# -- integration-level properties ---------------------------------------------


def test_excretion_runs_strict_accumulates_and_closes_carbon(params, store):
    # Under the strict touches contract: excretion writes only its reservoir pool (+ the
    # sugar draw), the reservoir accumulates non-negatively and stays TRACE, and — because
    # the C5 carbon is drawn from sugar AND acetolactate is weighted in total_carbon — a
    # full ferment conserves carbon to machine precision (the reservoir is real
    # carbon-accounted state, not an unbooked leak).
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2(), AcetolactateExcretion()],
        strict=True,
    )
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    assert_nonnegative(traj, ("acetolactate", "S"), atol=1e-12)
    # The reservoir actually fills (mechanism live) but stays trace (~mg/L), orders of
    # magnitude below the g/L ethanol flux.
    assert 0.0 < float(traj.series("acetolactate")[-1]) < 0.1
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_excretion_perturbs_only_sugar_per_rhs(params):
    # Isolability (prime directive #3) under a1: enabling excretion leaves dX/dN/dE/dCO2
    # byte-for-byte (it never touches them); the only change is on dS, which loses EXACTLY
    # the carbon deposited in the reservoir. That per-RHS balance is why total_carbon closes.
    schema = wine_schema()
    core = ProcessSet(schema, [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2()])
    with_vdk = ProcessSet(
        schema, [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2(), AcetolactateExcretion()]
    )
    for s, e, n in ((245.0, 0.0, 0.08), (120.0, 60.0, 0.01), (5.0, 120.0, 0.0)):
        y = _wine_y0(schema, x=1.5, s=s, e=e, n=n)
        d_core = core.total_derivatives(0.0, y, params)
        d_vdk = with_vdk.total_derivatives(0.0, y, params)
        for var in ("X", "E", "N", "CO2"):
            assert d_core[schema.slice(var)] == pytest.approx(d_vdk[schema.slice(var)], abs=0.0)
        delta_s = float(d_vdk[schema.slice("S")][0] - d_core[schema.slice("S")][0])
        ala_rate = float(d_vdk[schema.slice("acetolactate")][0])
        assert delta_s <= 0.0  # sugar drawn down, never created
        carbon_residual = delta_s * _GLUCOSE_C + ala_rate * _ACETOLACTATE_C
        assert carbon_residual == pytest.approx(0.0, abs=1e-12)


# -- tier propagation ---------------------------------------------------------


def test_excretion_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AcetolactateExcretion()])
    # Speculative form and speculative params ⇒ the reservoir output is speculative.
    assert ps.tier_of("acetolactate") is Tier.SPECULATIVE
    assert ps.tier_of("acetolactate", store.tier_map()) is Tier.SPECULATIVE
