"""Tests for the excreted keto-acid overflow pools — pyruvate (D-49) and alpha-ketoglutarate (D-50).

Both are the second- and third-strongest SO₂-binding carbonyls after acetaldehyde. Unlike the
acetaldehyde *buffer* (D-27, an ethanol-carbon borrow), overflow keto-acids are modelled as
**excreted side pools**: each is drawn out of sugar during active fermentation and re-assimilated
back to ethanol + CO₂. The load-bearing modelling choice is that re-assimilation is
**flux-linked (co-metabolic)**, NOT the no-flux/viable-``X`` ADH idiom — so at dryness both
terms die and the pool *freezes* at a persistent finished-wine residual (crash- and
duration-independent), rather than draining to ~0 over the long tail (which a viable-``X`` gate
would do, since a clean ferment ends with the yeast still viable).

* **Excretion** (:class:`PyruvateExcretion` / :class:`AlphaKetoglutarateExcretion`): flux-linked,
  temperature-flat; draws carbon out of ``S`` at the substrate's own carbon fraction.
* **Reassimilation** (:class:`PyruvateReassimilation` / :class:`AlphaKetoglutarateReassimilation`):
  flux-linked, carbon-closing to ethanol + CO₂; stops at dryness (freezing the residual), the
  *opposite* of the acetaldehyde reduction. Pyruvate's C3 → C2 + C1 is mole-for-mole (it happens
  to BE the Gay-Lussac 2:1 carbon ratio); α-KG's C5 does not divide evenly, so its reassimilation
  returns carbon at that SAME 2:1 ratio (5/3 mol ethanol + 5/3 mol CO2 per mole) instead of
  mole-for-mole — the one structural difference between the two Processes (D-50).

The unit tests pin each Process's closed form, the sugar-draw / release carbon accounting and
the guards; the acceptance section verifies the emergent persistent residual (in the real
finished-wine range), its crash- and duration-independence (the dryness freeze), machine-
precision carbon closure on a compiled run, and — the isolability guarantee — that the §2.2
ABV/CO₂ endpoints are preserved to ≪ 0.1 % with both pools wired in.
"""

import pytest

from fermentation.core.chemistry import (
    M_ALPHA_KETOGLUTARATE,
    M_CO2,
    M_ETHANOL,
    M_PYRUVATE,
    carbon_mass_fraction,
)
from fermentation.core.kinetics import (
    AlphaKetoglutarateExcretion,
    AlphaKetoglutarateReassimilation,
    PyruvateExcretion,
    PyruvateReassimilation,
)
from fermentation.core.media import beer_schema, get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units.convert import gpl_to_mgl
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon

#: Carbon fractions the draw/release books against (mirror the chemistry constants).
_PYRUVATE_C = carbon_mass_fraction("pyruvate")
_ALPHA_KG_C = carbon_mass_fraction("alpha_ketoglutarate")
_GLUCOSE_C = carbon_mass_fraction("glucose")  # wine's single sugar slot
_ETHANOL_C = carbon_mass_fraction("ethanol")
_CO2_C = carbon_mass_fraction("CO2")
#: The Gay-Lussac split unit count for alpha-KG's C5 (see AlphaKetoglutarateReassimilation): 5/3
#: mol ethanol + 5/3 mol CO2 per mole reassimilated, NOT mole-for-mole like pyruvate's C3.
_ALPHA_KG_SPLIT_UNITS = 5.0 / 3.0


@pytest.fixture
def store():
    # Wine kinetics PLUS the shared keto-acid constants (k_pyruvate_excretion, …).
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "keto_acids.yaml",
    )


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x: float = 2.0,
    s: float = 200.0,
    e: float = 50.0,
    n: float = 0.1,
    t: float = 293.15,
    pyruvate: float = 0.0,
    alpha_kg: float = 0.0,
) -> FloatArray:
    y = schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": t, "CO2": 0.0})
    y[schema.slice("pyruvate")] = pyruvate
    y[schema.slice("alpha_ketoglutarate")] = alpha_kg
    return y


# -- metadata -----------------------------------------------------------------


def test_excretion_metadata():
    p = PyruvateExcretion()
    assert p.name == "pyruvate_excretion"
    # Speculative: rate magnitude is an order-of-magnitude estimate.
    assert p.tier is Tier.SPECULATIVE
    # Draws carbon out of sugar into its own pool — touches pyruvate and S only (D-49).
    assert set(p.touches) == {"pyruvate", "S"}
    assert set(p.reads) == {"k_pyruvate_excretion", "K_sugar_uptake"}


def test_reassimilation_metadata():
    p = PyruvateReassimilation()
    assert p.name == "pyruvate_reassimilation"
    assert p.tier is Tier.SPECULATIVE
    # Returns carbon to ethanol + CO2 (C3 → C2 + C1) — touches pyruvate, E, CO2.
    assert set(p.touches) == {"pyruvate", "E", "CO2"}
    # FLUX-LINKED: reads K_sugar_uptake (co-metabolic), unlike the no-flux ADH reduction (D-27).
    assert set(p.reads) == {"k_pyruvate_reassimilation", "K_sugar_uptake"}


# -- excretion closed form & guards -------------------------------------------


def test_excretion_matches_closed_form(params):
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = _wine_y0(schema, x=x, s=s)
    d = PyruvateExcretion().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    rate = params["k_pyruvate_excretion"] * flux
    assert schema.get(d, "pyruvate") == pytest.approx(rate)
    # Carbon is drawn OUT OF SUGAR (booked at the C3 pyruvate fraction): the carbon leaving S
    # equals the carbon entering pyruvate, so the draw is carbon-exact per RHS.
    assert schema.get(d, "S") == pytest.approx(-rate * _PYRUVATE_C / _GLUCOSE_C)
    assert schema.get(d, "S") * _GLUCOSE_C == pytest.approx(
        -schema.get(d, "pyruvate") * _PYRUVATE_C
    )
    # Nothing else moves — no ethanol/CO2 borrow (that is the re-assimilation half).
    for var in ("X", "E", "N", "CO2"):
        assert schema.get(d, var) == 0.0


def test_excretion_is_temperature_flat(params):
    # Documented v1 simplification (D-49): excretion carries NO Arrhenius factor. Identical
    # cold vs warm at the same flux.
    schema = wine_schema()
    cold = PyruvateExcretion().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = PyruvateExcretion().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(cold, "pyruvate") == pytest.approx(schema.get(warm, "pyruvate"))
    assert schema.get(cold, "pyruvate") > 0.0


def test_excretion_scales_with_fermentative_flux(params):
    # Coupled to the biomass-catalysed sugar flux (linear in X): twice the biomass ⇒ twice the
    # excretion.
    schema = wine_schema()
    r1 = PyruvateExcretion().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = PyruvateExcretion().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "pyruvate") == pytest.approx(2.0 * schema.get(r1, "pyruvate"))


def test_excretion_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    no_x = PyruvateExcretion().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)
    no_s = PyruvateExcretion().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)
    assert schema.get(no_x, "pyruvate") == 0.0
    assert schema.get(no_s, "pyruvate") == 0.0
    # No excretion ⇒ sugar untouched too.
    assert schema.get(no_x, "S") == 0.0
    assert schema.get(no_s, "S") == 0.0


# -- reassimilation closed form & guards --------------------------------------


def test_reassimilation_matches_closed_form(params):
    schema = wine_schema()
    x, s, pyr = 2.0, 200.0, 0.03
    y = _wine_y0(schema, x=x, s=s, pyruvate=pyr)
    d = PyruvateReassimilation().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    loss = params["k_pyruvate_reassimilation"] * flux * pyr
    r = loss / M_PYRUVATE  # molar turnover: 1 pyruvate → 1 ethanol + 1 CO₂ (C3 → C2 + C1)
    assert schema.get(d, "pyruvate") == pytest.approx(-loss)
    assert schema.get(d, "E") == pytest.approx(r * M_ETHANOL)
    assert schema.get(d, "CO2") == pytest.approx(r * M_CO2)
    # Carbon-exact per RHS: the C3 lost from pyruvate = the C2 into ethanol + the C1 into CO2.
    moved = (
        schema.get(d, "pyruvate") * _PYRUVATE_C
        + schema.get(d, "E") * _ETHANOL_C
        + schema.get(d, "CO2") * _CO2_C
    )
    assert moved == pytest.approx(0.0, abs=1e-15)
    # Never touches sugar/nitrogen/biomass.
    for var in ("S", "N", "X"):
        assert schema.get(d, var) == 0.0


def test_reassimilation_stops_at_dryness(params):
    # THE load-bearing difference from the acetaldehyde reduction (D-27): re-assimilation is
    # FLUX-LINKED, so with sugar gone it STOPS even though pyruvate and viable yeast remain —
    # freezing the pool at its dryness value (the persistent residual). A no-flux gate would
    # instead keep draining it to ~0.
    schema = wine_schema()
    d = PyruvateReassimilation().derivatives(
        0.0, _wine_y0(schema, s=0.0, x=2.0, pyruvate=0.03), schema, params
    )
    assert schema.get(d, "pyruvate") == 0.0
    assert schema.get(d, "E") == 0.0
    assert schema.get(d, "CO2") == 0.0


def test_reassimilation_runs_during_active_ferment(params):
    # The flip side: while sugar remains (active ferment) it DOES clear pyruvate, returning it
    # to ethanol + CO2.
    schema = wine_schema()
    d = PyruvateReassimilation().derivatives(
        0.0, _wine_y0(schema, s=150.0, x=2.0, pyruvate=0.03), schema, params
    )
    assert schema.get(d, "pyruvate") < 0.0
    assert schema.get(d, "E") > 0.0
    assert schema.get(d, "CO2") > 0.0


def test_reassimilation_zero_without_biomass(params):
    # No viable yeast ⇒ flux is 0 ⇒ no re-assimilation (the pool is likewise stranded).
    schema = wine_schema()
    d = PyruvateReassimilation().derivatives(
        0.0, _wine_y0(schema, x=0.0, pyruvate=0.03), schema, params
    )
    assert schema.get(d, "pyruvate") == 0.0
    assert schema.get(d, "E") == 0.0


def test_reassimilation_zero_and_clamped_without_pyruvate(params):
    schema = wine_schema()
    empty = PyruvateReassimilation().derivatives(
        0.0, _wine_y0(schema, pyruvate=0.0), schema, params
    )
    negative = PyruvateReassimilation().derivatives(
        0.0, _wine_y0(schema, pyruvate=-1e-6), schema, params
    )
    assert schema.get(empty, "pyruvate") == 0.0
    # A solver undershoot below zero cannot manufacture ethanol from a negative pool.
    assert schema.get(negative, "pyruvate") == 0.0
    assert schema.get(negative, "E") == 0.0
    assert schema.get(negative, "CO2") == 0.0


# -- carbon: draw then release is carbon-neutral on the ledger ----------------


def test_excrete_then_reassimilate_is_carbon_neutral_on_the_ledger(params):
    # Each Process moves zero net carbon across the weighted ledger: excretion moves S → pyruvate
    # (carbon-exact), re-assimilation moves pyruvate → E + CO2 (carbon-exact). Together the pool is
    # a carbon no-op — carbon parked as pyruvate is exactly the carbon withheld from S.
    schema = wine_schema()
    y = _wine_y0(schema, pyruvate=0.03)
    exc = PyruvateExcretion().derivatives(0.0, y, schema, params)
    rea = PyruvateReassimilation().derivatives(0.0, y, schema, params)
    exc_moved = schema.get(exc, "pyruvate") * _PYRUVATE_C + schema.get(exc, "S") * _GLUCOSE_C
    rea_moved = (
        schema.get(rea, "pyruvate") * _PYRUVATE_C
        + schema.get(rea, "E") * _ETHANOL_C
        + schema.get(rea, "CO2") * _CO2_C
    )
    assert exc_moved == pytest.approx(0.0, abs=1e-15)
    assert rea_moved == pytest.approx(0.0, abs=1e-15)


# ==============================================================================
# Alpha-ketoglutarate (decision D-50) — same structure, Gay-Lussac carbon split
# ==============================================================================

# -- metadata -----------------------------------------------------------------


def test_alpha_kg_excretion_metadata():
    p = AlphaKetoglutarateExcretion()
    assert p.name == "alpha_kg_excretion"
    assert p.tier is Tier.SPECULATIVE
    assert set(p.touches) == {"alpha_ketoglutarate", "S"}
    assert set(p.reads) == {"k_alpha_kg_excretion", "K_sugar_uptake"}


def test_alpha_kg_reassimilation_metadata():
    p = AlphaKetoglutarateReassimilation()
    assert p.name == "alpha_kg_reassimilation"
    assert p.tier is Tier.SPECULATIVE
    assert set(p.touches) == {"alpha_ketoglutarate", "E", "CO2"}
    assert set(p.reads) == {"k_alpha_kg_reassimilation", "K_sugar_uptake"}


# -- excretion closed form & guards -------------------------------------------


def test_alpha_kg_excretion_matches_closed_form(params):
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = _wine_y0(schema, x=x, s=s)
    d = AlphaKetoglutarateExcretion().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    rate = params["k_alpha_kg_excretion"] * flux
    assert schema.get(d, "alpha_ketoglutarate") == pytest.approx(rate)
    assert schema.get(d, "S") == pytest.approx(-rate * _ALPHA_KG_C / _GLUCOSE_C)
    assert schema.get(d, "S") * _GLUCOSE_C == pytest.approx(
        -schema.get(d, "alpha_ketoglutarate") * _ALPHA_KG_C
    )
    for var in ("X", "E", "N", "CO2", "pyruvate"):
        assert schema.get(d, var) == 0.0


def test_alpha_kg_excretion_is_temperature_flat(params):
    schema = wine_schema()
    proc = AlphaKetoglutarateExcretion()
    cold = proc.derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = proc.derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(cold, "alpha_ketoglutarate") == pytest.approx(
        schema.get(warm, "alpha_ketoglutarate")
    )
    assert schema.get(cold, "alpha_ketoglutarate") > 0.0


def test_alpha_kg_excretion_scales_with_fermentative_flux(params):
    schema = wine_schema()
    r1 = AlphaKetoglutarateExcretion().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = AlphaKetoglutarateExcretion().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "alpha_ketoglutarate") == pytest.approx(
        2.0 * schema.get(r1, "alpha_ketoglutarate")
    )


def test_alpha_kg_excretion_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    no_x = AlphaKetoglutarateExcretion().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)
    no_s = AlphaKetoglutarateExcretion().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)
    assert schema.get(no_x, "alpha_ketoglutarate") == 0.0
    assert schema.get(no_s, "alpha_ketoglutarate") == 0.0
    assert schema.get(no_x, "S") == 0.0
    assert schema.get(no_s, "S") == 0.0


# -- reassimilation closed form & guards (the Gay-Lussac split) ---------------


def test_alpha_kg_reassimilation_matches_closed_form(params):
    # THE load-bearing difference from pyruvate: the carbon returns at the Gay-Lussac 2:1 ratio
    # (5/3 mol ethanol + 5/3 mol CO2 per mole alpha-KG), NOT mole-for-mole.
    schema = wine_schema()
    x, s, akg = 2.0, 200.0, 0.02
    y = _wine_y0(schema, x=x, s=s, alpha_kg=akg)
    d = AlphaKetoglutarateReassimilation().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    loss = params["k_alpha_kg_reassimilation"] * flux * akg
    molar_turnover = loss / M_ALPHA_KETOGLUTARATE
    units = molar_turnover * _ALPHA_KG_SPLIT_UNITS
    assert schema.get(d, "alpha_ketoglutarate") == pytest.approx(-loss)
    assert schema.get(d, "E") == pytest.approx(units * M_ETHANOL)
    assert schema.get(d, "CO2") == pytest.approx(units * M_CO2)
    # Carbon-exact per RHS: the C5 lost from alpha-KG = the ethanol-carbon + the CO2-carbon.
    moved = (
        schema.get(d, "alpha_ketoglutarate") * _ALPHA_KG_C
        + schema.get(d, "E") * _ETHANOL_C
        + schema.get(d, "CO2") * _CO2_C
    )
    assert moved == pytest.approx(0.0, abs=1e-15)
    for var in ("S", "N", "X", "pyruvate"):
        assert schema.get(d, var) == 0.0


def test_alpha_kg_reassimilation_is_not_mole_for_mole():
    # Guard against silently reverting to pyruvate's mole-for-mole form: for a C5 species that
    # would be a CO2-heavy 1+4 split, not the carbon-exact 2:1 Gay-Lussac ratio this Process uses.
    assert pytest.approx(5.0 / 3.0) == _ALPHA_KG_SPLIT_UNITS
    assert pytest.approx(1.0) != _ALPHA_KG_SPLIT_UNITS


def test_alpha_kg_reassimilation_stops_at_dryness(params):
    schema = wine_schema()
    d = AlphaKetoglutarateReassimilation().derivatives(
        0.0, _wine_y0(schema, s=0.0, x=2.0, alpha_kg=0.02), schema, params
    )
    assert schema.get(d, "alpha_ketoglutarate") == 0.0
    assert schema.get(d, "E") == 0.0
    assert schema.get(d, "CO2") == 0.0


def test_alpha_kg_reassimilation_runs_during_active_ferment(params):
    schema = wine_schema()
    d = AlphaKetoglutarateReassimilation().derivatives(
        0.0, _wine_y0(schema, s=150.0, x=2.0, alpha_kg=0.02), schema, params
    )
    assert schema.get(d, "alpha_ketoglutarate") < 0.0
    assert schema.get(d, "E") > 0.0
    assert schema.get(d, "CO2") > 0.0


def test_alpha_kg_reassimilation_zero_without_biomass(params):
    schema = wine_schema()
    d = AlphaKetoglutarateReassimilation().derivatives(
        0.0, _wine_y0(schema, x=0.0, alpha_kg=0.02), schema, params
    )
    assert schema.get(d, "alpha_ketoglutarate") == 0.0
    assert schema.get(d, "E") == 0.0


def test_alpha_kg_reassimilation_zero_and_clamped_without_substrate(params):
    schema = wine_schema()
    empty = AlphaKetoglutarateReassimilation().derivatives(
        0.0, _wine_y0(schema, alpha_kg=0.0), schema, params
    )
    negative = AlphaKetoglutarateReassimilation().derivatives(
        0.0, _wine_y0(schema, alpha_kg=-1e-6), schema, params
    )
    assert schema.get(empty, "alpha_ketoglutarate") == 0.0
    assert schema.get(negative, "alpha_ketoglutarate") == 0.0
    assert schema.get(negative, "E") == 0.0
    assert schema.get(negative, "CO2") == 0.0


# -- carbon: draw then release is carbon-neutral on the ledger ----------------


def test_alpha_kg_excrete_then_reassimilate_is_carbon_neutral_on_the_ledger(params):
    schema = wine_schema()
    y = _wine_y0(schema, alpha_kg=0.02)
    exc = AlphaKetoglutarateExcretion().derivatives(0.0, y, schema, params)
    rea = AlphaKetoglutarateReassimilation().derivatives(0.0, y, schema, params)
    exc_moved = (
        schema.get(exc, "alpha_ketoglutarate") * _ALPHA_KG_C + schema.get(exc, "S") * _GLUCOSE_C
    )
    rea_moved = (
        schema.get(rea, "alpha_ketoglutarate") * _ALPHA_KG_C
        + schema.get(rea, "E") * _ETHANOL_C
        + schema.get(rea, "CO2") * _CO2_C
    )
    assert exc_moved == pytest.approx(0.0, abs=1e-15)
    assert rea_moved == pytest.approx(0.0, abs=1e-15)


# -- wine-only wiring ---------------------------------------------------------


def test_pool_is_wine_only():
    # v1 scope (D-49/D-50): the SO₂-binding competition is a wine readout, no §2.2 beer benchmark
    # asserts a keto-acid level — so both pools are wired into wine only. Beer has neither slot
    # nor either keto-acid Process.
    assert "pyruvate" in wine_schema()
    assert "alpha_ketoglutarate" in wine_schema()
    assert "pyruvate" not in beer_schema()
    assert "alpha_ketoglutarate" not in beer_schema()
    wine_procs = {p().name for p in get_medium("wine").process_factories}
    beer_procs = {p().name for p in get_medium("beer").process_factories}
    keto_acid_procs = {
        "pyruvate_excretion",
        "pyruvate_reassimilation",
        "alpha_kg_excretion",
        "alpha_kg_reassimilation",
    }
    assert keto_acid_procs <= wine_procs
    assert beer_procs.isdisjoint(keto_acid_procs)


# -- acceptance: persistent residual, freeze, carbon, isolability -------------


def _wine_run(days: float = 21.0):
    scenario = Scenario(
        name="wine-pyruvate",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.5,
            "tartaric_gpl": 6.0,
            "malic_gpl": 3.0,
            "initial_ph": 3.4,
        },  # fmt: skip
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=days,
    )
    compiled = compile_scenario(scenario, strict=True)
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    assert traj.success, traj.message
    return traj, compiled


def test_pool_strands_a_persistent_finished_wine_residual():
    # The defining emergent behaviour (D-49): overflow pyruvate is drawn from sugar during the
    # ferment and FREEZES at dryness to a persistent residual in the real finished-wine range
    # (Jackowetz & Mira de Orduña 2013 — 10s-100s mg/L; nominal ratio targets ~30). Crucially the
    # residual persists WITH the yeast still viable (no crash needed) — the co-metabolic freeze.
    traj, _ = _wine_run(21.0)
    resid = gpl_to_mgl(float(traj.series("pyruvate")[-1]))
    assert 10.0 < resid < 60.0, f"finished-wine pyruvate {resid:.1f} mg/L off-range"
    # The yeast is STILL VIABLE at the end — the residual is a co-metabolic freeze, not a crash.
    assert float(traj.series("X")[-1]) > 0.1
    assert_nonnegative(traj, ("pyruvate",), atol=1e-12)


def test_alpha_kg_pool_strands_a_persistent_finished_wine_residual():
    # D-50 analogue: alpha-KG also freezes at dryness, targeted somewhat lower than pyruvate
    # (nominal ratio ~20 mg/L — alpha-KG is typically less abundant in finished wine).
    traj, _ = _wine_run(21.0)
    resid = gpl_to_mgl(float(traj.series("alpha_ketoglutarate")[-1]))
    assert 5.0 < resid < 50.0, f"finished-wine alpha-KG {resid:.1f} mg/L off-range"
    assert float(traj.series("X")[-1]) > 0.1
    assert_nonnegative(traj, ("alpha_ketoglutarate",), atol=1e-12)


def test_residual_is_duration_independent():
    # The dryness freeze makes the residual crash- AND duration-independent: extending the run
    # from 21 to 40 days leaves it unchanged (both terms are dead at dryness, so nothing drains
    # the frozen pool). A no-flux viable-X gate would instead keep draining it over the long tail.
    r21 = gpl_to_mgl(float(_wine_run(21.0)[0].series("pyruvate")[-1]))
    r40 = gpl_to_mgl(float(_wine_run(40.0)[0].series("pyruvate")[-1]))
    assert r21 == pytest.approx(r40, rel=1e-3)


def test_alpha_kg_residual_is_duration_independent():
    r21 = gpl_to_mgl(float(_wine_run(21.0)[0].series("alpha_ketoglutarate")[-1]))
    r40 = gpl_to_mgl(float(_wine_run(40.0)[0].series("alpha_ketoglutarate")[-1]))
    assert r21 == pytest.approx(r40, rel=1e-3)


def test_carbon_closes_on_a_compiled_run():
    # With both keto-acid pools wired into wine (D-49, D-50), a full compiled ferment conserves
    # carbon to machine precision — every draw/release is entirely on the weighted ledger, each
    # pool at its own carbon fraction (pyruvate C3, alpha-KG C5 at the Gay-Lussac split). Non-
    # trivial: both pools hold real sugar carbon at the end.
    traj, compiled = _wine_run(21.0)
    f_c = compiled.param_values["biomass_C_fraction"]
    assert gpl_to_mgl(float(traj.series("pyruvate")[-1])) > 10.0  # genuinely accumulates
    assert gpl_to_mgl(float(traj.series("alpha_ketoglutarate")[-1])) > 5.0  # genuinely accumulates
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )


def test_abv_co2_endpoints_preserved_by_the_pools():
    # ISOLABILITY (prime directive #3): both pools route a trace of sugar carbon on a detour to
    # ethanol (parking only the frozen residuals), so — unlike the byte-for-byte acetaldehyde
    # buffer — the ABV/CO2 endpoints are not bit-identical to the pools-off core. But the delta is
    # ≪ 0.1 %, so the §2.2 ABV / realised-yield / CO2 benchmarks are preserved to far below any
    # tolerance. Compared at the endpoint, not pointwise.
    on, _ = _wine_run(21.0)
    off_c = compile_scenario(
        Scenario(
            name="wine-keto-acids-off", medium="wine",
            initial={
                "brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5,
                "tartaric_gpl": 6.0, "malic_gpl": 3.0, "initial_ph": 3.4,
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)], duration_days=21.0,
        ),
        strict=True,
    )  # fmt: skip
    for n in (
        "pyruvate_excretion",
        "pyruvate_reassimilation",
        "alpha_kg_excretion",
        "alpha_kg_reassimilation",
    ):
        off_c.process_set.disable(n)
    off = simulate(off_c.process_set, off_c.param_values, off_c.y0, off_c.t_span_h)
    assert off.success, off.message
    for var in ("E", "CO2"):
        v_on = float(on.series(var)[-1])
        v_off = float(off.series(var)[-1])
        assert v_on == pytest.approx(v_off, rel=1e-3), f"{var} endpoint on={v_on} off={v_off}"


# -- tier propagation ---------------------------------------------------------


def test_excretion_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [PyruvateExcretion()])
    assert ps.tier_of("pyruvate") is Tier.SPECULATIVE
    assert ps.tier_of("pyruvate", store.tier_map()) is Tier.SPECULATIVE


def test_reassimilation_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [PyruvateReassimilation()])
    assert ps.tier_of("pyruvate") is Tier.SPECULATIVE
    assert ps.tier_of("pyruvate", store.tier_map()) is Tier.SPECULATIVE


def test_alpha_kg_excretion_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AlphaKetoglutarateExcretion()])
    assert ps.tier_of("alpha_ketoglutarate") is Tier.SPECULATIVE
    assert ps.tier_of("alpha_ketoglutarate", store.tier_map()) is Tier.SPECULATIVE


def test_alpha_kg_reassimilation_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AlphaKetoglutarateReassimilation()])
    assert ps.tier_of("alpha_ketoglutarate") is Tier.SPECULATIVE
    assert ps.tier_of("alpha_ketoglutarate", store.tier_map()) is Tier.SPECULATIVE
