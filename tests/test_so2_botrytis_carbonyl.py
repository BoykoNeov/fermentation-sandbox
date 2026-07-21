"""Botrytis 5-oxofructose SO₂ binding — the fourth competing carbonyl (decision D-130).

The SO₂-binding axis was already 3-carbonyl (acetaldehyde D-28, pyruvate/α-KG D-51). Barbe 2000
(via Handbook of Enology Vol 1 §8.4) adds the BOTRYTIS-specific carbonyls: 5-oxofructose is the
dominant one — a hexodiulose *Botrytis cinerea* makes by oxidising must fructose ON THE BERRY, so
it enters as a must-composition INPUT (the inert ``oxofructose`` slot, ``oxofructose_mgl`` scenario
key), is yeast-INERT ("not altered by alcoholic fermentation") and so persists into the bottle,
dominating a botrytized wine's high SO₂-combining power. Dihydroxyacetone (transient — reduced by
yeast during AF) and gluconolactone (Kd doesn't reconcile with its sourced binding point; the
gluconic⇌lactone split is unsourced) are DEFERRED, documented in ``K_5_oxofructose_so2`` notes.

Two headline gates (the D-129 discriminating-check pattern):
  * GATE 1 — ``test_gate1_isolable_byte_for_byte_at_zero_load``: at load 0 the binding is
    byte-for-byte the D-51 3-carbonyl form (structural, not param-luck), so every non-botrytis
    wine is unchanged.
  * GATE 2 — ``test_gate2_botrytis_load_collapses_molecular_so2``: at the sourced ~100 mg/L load
    5-oxofructose binds a Barbe-band share of SO₂ and collapses molecular (antimicrobial) SO₂ —
    opening the MLF/Brett gates, the reason a botrytized must needs a higher dose.
"""

import pytest

from fermentation.core import acidbase
from fermentation.core.acidbase import (
    _acetaldehyde_molar,
    _alpha_kg_molar,
    _oxofructose_molar,
    _pyruvate_molar,
    bisulfite_fraction,
    bound_so2_molar,
)
from fermentation.core.chemistry import (
    CARBON_ATOMS,
    M_5_OXOFRUCTOSE,
    M_MALIC,
    M_SO2,
    M_TARTARIC,
    NITROGEN_ATOMS,
    carbon_mass_fraction,
)
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units.convert import mgl_to_gpl


@pytest.fixture
def pset():
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    return pset.resolve()


def _anchor_cation(
    params: dict[str, float], tartaric_gpl: float, malic_gpl: float, target_ph: float
) -> float:
    totals = {
        "tartaric": tartaric_gpl / M_TARTARIC,
        "malic": malic_gpl / M_MALIC,
        "lactic": 0.0,
    }
    return acidbase.solve_cation_charge(totals, 0.0, acidbase.build_pka_map(params), target_ph)


def _wine_state(schema: StateSchema, **slots: float) -> FloatArray:
    # A finished-wine background: dry, ethanol present, typical carbonyl residuals.
    base: dict[str, float | list[float]] = {
        "X": 0.5, "S": [0.0], "E": 100.0, "N": 0.0, "T": 293.15, "CO2": 0.0,
    }  # fmt: skip
    base.update(slots)
    return schema.pack(base)


def _finished_wine(
    schema: StateSchema,
    params: dict[str, float],
    *,
    so2_mgl: float,
    oxofructose_mgl: float,
    ph: float = 3.5,
) -> FloatArray:
    """A dry wine at `ph` with typical acetaldehyde/pyruvate/α-KG residuals + the two doses."""
    cation = _anchor_cation(params, 6.0, 3.0, ph)
    return _wine_state(
        schema, tartaric=6.0, malic=3.0, cation_charge=cation,
        acetaldehyde=mgl_to_gpl(30.0), pyruvate=mgl_to_gpl(30.0),
        alpha_ketoglutarate=mgl_to_gpl(20.0),
        so2_total=mgl_to_gpl(so2_mgl), oxofructose=mgl_to_gpl(oxofructose_mgl),
    )  # fmt: skip


# -- 1. Provenance / metadata -------------------------------------------------


def test_kd_param_present_and_plausible(pset):
    p = pset["K_5_oxofructose_so2"]
    assert p.value == pytest.approx(3.3e-4)  # ~0.33 mM, B&S 1973b via Barbe 2000
    assert p.unit == "mol/L"
    assert p.tier is Tier.PLAUSIBLE  # measured Kd + measured load, apparent-in-wine


def test_chemistry_registry_6_carbon_nitrogen_free():
    # C6H10O6: mass ~178.14, six carbons, nitrogen-free (an oxidised sugar).
    assert pytest.approx(178.14, abs=0.05) == M_5_OXOFRUCTOSE
    assert CARBON_ATOMS["5_oxofructose"] == 6
    assert NITROGEN_ATOMS["5_oxofructose"] == 0
    assert carbon_mass_fraction("5_oxofructose") == pytest.approx(
        6 * 12.011 / M_5_OXOFRUCTOSE, rel=1e-6
    )


def test_wine_only_absent_from_beer():
    # A botrytis grape input has no meaning for beer — wine-only, like the acids/SO₂/keto-acids.
    assert "oxofructose" in wine_schema()
    assert "oxofructose" not in beer_schema()


# -- 2. GATE 1: structural isolability ---------------------------------------


def test_gate1_isolable_byte_for_byte_at_zero_load(params):
    # At load 0 the 4-carbonyl split is byte-for-byte the D-51 3-carbonyl form — the Langmuir
    # term is EXACTLY 0, so every non-botrytis wine is unchanged. Verified across a dose sweep
    # (structural, not param-luck) against an INDEPENDENT hand-computed 3-carbonyl bound.
    schema = wine_schema()
    max_gap = 0.0
    for dose_mgl in (20.0, 60.0, 120.0, 200.0):
        y = _finished_wine(schema, params, so2_mgl=dose_mgl, oxofructose_mgl=0.0)
        s = acidbase.speciate_so2(y, schema, params)
        pka = tuple(params[n] for n in acidbase.SO2_PKA_PARAM_NAMES)
        beta = bisulfite_fraction(10.0 ** (-s.ph), pka)
        acet, pyr, akg = bound_so2_molar(
            mgl_to_gpl(dose_mgl) / M_SO2,
            (
                (_acetaldehyde_molar(y, schema), params["K_acetaldehyde_so2"]),
                (_pyruvate_molar(y, schema), params["K_pyruvate_so2"]),
                (_alpha_kg_molar(y, schema), params["K_alpha_kg_so2"]),
            ),
            beta,
        )
        bound_d51 = (acet + pyr + akg) * M_SO2
        max_gap = max(max_gap, abs(s.bound - bound_d51))
    assert max_gap == 0.0  # exact — no fourth term contributes at load 0


# -- 3. GATE 2: the sourced botrytis outcome ---------------------------------


def test_gate2_botrytis_load_collapses_molecular_so2(params):
    # THE headline: adding a sourced ~100 mg/L botrytis 5-oxofructose load to an otherwise
    # identical finished wine sequesters SO₂, so molecular (antimicrobial) SO₂ collapses — the
    # emergent mechanism that OPENS the MLF/Brett antimicrobial gates (a botrytized must needs a
    # higher SO₂ dose for the same protection). Free SO₂ drops too; total is conserved.
    schema = wine_schema()
    dry = _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=0.0)
    botrytis = _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=100.0)
    s_dry = acidbase.speciate_so2(dry, schema, params)
    s_bot = acidbase.speciate_so2(botrytis, schema, params)
    assert s_bot.bound > s_dry.bound  # more SO₂ bound
    assert s_bot.free < s_dry.free  # less free
    assert s_bot.molecular < 0.85 * s_dry.molecular  # antimicrobial pool collapses (gate opens)
    # total conserved in both (free + bound == total)
    for s in (s_dry, s_bot):
        assert s.free + s.bound == pytest.approx(s.total)


def test_oxofructose_share_within_barbe_band(params):
    # At the sourced ~100 mg/L load, 5-oxofructose's OWN share of bound SO₂ lands inside Barbe's
    # measured "4-78% of combinations" band. Low end here is expected — this dry-wine background
    # carries heavy competing yeast-carbonyl residuals; an advanced-botrytis sweet must (higher
    # load, relatively less yeast carbonyl) pushes the share toward Barbe's cited ~60%.
    schema = wine_schema()
    y = _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=100.0)
    s = acidbase.speciate_so2(y, schema, params)
    pka = tuple(params[n] for n in acidbase.SO2_PKA_PARAM_NAMES)
    beta = bisulfite_fraction(10.0 ** (-s.ph), pka)
    _, _, _, oxo = acidbase._bound_molar_split(
        mgl_to_gpl(60.0) / M_SO2,
        _acetaldehyde_molar(y, schema),
        _pyruvate_molar(y, schema),
        _alpha_kg_molar(y, schema),
        _oxofructose_molar(y, schema),
        beta,
        params,
    )
    share = (oxo * M_SO2) / s.bound
    assert 0.04 <= share <= 0.78


def test_binding_monotone_in_load(params):
    # More botrytis carbonyl ⇒ monotonically more bound / less free / less molecular SO₂.
    schema = wine_schema()
    loads = [0.0, 50.0, 100.0, 150.0, 250.0]
    specs = [
        acidbase.speciate_so2(
            _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=L), schema, params
        )
        for L in loads
    ]
    assert all(a.bound < b.bound for a, b in zip(specs, specs[1:], strict=False))
    assert all(a.free > b.free for a, b in zip(specs, specs[1:], strict=False))
    assert all(a.molecular > b.molecular for a, b in zip(specs, specs[1:], strict=False))


def test_oxofructose_frees_acetaldehyde_for_adh(params):
    # A competing binder soaks up part of the shared bisulfite pool, so acetaldehyde itself is
    # LESS bound ⇒ free_acetaldehyde (the ADH-reducible share, D-47) RISES — the same competition
    # the keto-acids show (D-51), now driven by the botrytis input. Never exceeds total present.
    schema = wine_schema()
    ph = 3.5
    dry = _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=0.0, ph=ph)
    botrytis = _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=150.0, ph=ph)
    ph_dry = acidbase.ph_of_state(dry, schema, params)
    ph_bot = acidbase.ph_of_state(botrytis, schema, params)
    free_dry = acidbase.free_acetaldehyde(dry, schema, params, ph_dry)
    free_bot = acidbase.free_acetaldehyde(botrytis, schema, params, ph_bot)
    assert free_bot > free_dry
    assert free_bot <= mgl_to_gpl(30.0) + 1e-12  # never exceeds the acetaldehyde present


# -- 4. Conservation + wiring -------------------------------------------------


def test_inert_pool_is_carbon_conservation_neutral(params):
    # The oxofructose pool is weighted in total_carbon (C6) but INERT (no Process touches it), so
    # it is a constant term: two states differing only in the dosed load have total_carbon offset
    # by exactly the pool's carbon and neither changes over "time" (it has no derivative).
    from fermentation.validation import total_carbon

    schema = wine_schema()
    # biomass_carbon_fraction is required (schema has X) but cancels in the delta below — the two
    # states share X — so its exact value is irrelevant to what this test checks.
    tc = total_carbon(schema, biomass_carbon_fraction=0.48)
    dry = _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=0.0)
    botrytis = _finished_wine(schema, params, so2_mgl=60.0, oxofructose_mgl=100.0)
    delta = tc(botrytis) - tc(dry)
    expected = mgl_to_gpl(100.0) * carbon_mass_fraction("5_oxofructose")
    assert delta == pytest.approx(expected, rel=1e-9)


def test_scenario_wires_oxofructose_mgl(params):
    # oxofructose_mgl (mg/L) → the oxofructose slot (g/L); absent ⇒ 0 (a non-botrytis must).
    base = {"brix": 30.0, "pitch_gpl": 0.5, "yan_mgl": 250.0, "so2_total_mgl": 80.0}
    sweet = Scenario(
        name="sauternes",
        medium="wine",
        initial={**base, "oxofructose_mgl": 120.0},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        duration_days=1.0,
    )
    dry = Scenario(
        name="dry-white",
        medium="wine",
        initial=base,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        duration_days=1.0,
    )
    c_sweet = compile_scenario(sweet)
    c_dry = compile_scenario(dry)
    schema = c_sweet.schema
    sweet_oxo = float(c_sweet.y0[schema.slice("oxofructose")][0])
    assert sweet_oxo == pytest.approx(mgl_to_gpl(120.0))
    assert float(c_dry.y0[c_dry.schema.slice("oxofructose")][0]) == 0.0


def test_unknown_key_still_fails_loudly():
    # The allow-list guard is intact — a typo near the new key is rejected, not silently ignored.
    bad = Scenario(
        name="typo",
        medium="wine",
        initial={"brix": 24.0, "oxofructos_mgl": 100.0},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        duration_days=1.0,
    )
    with pytest.raises((ValueError, KeyError)):
        compile_scenario(bad)


def test_speciation_tier_stays_plausible(pset):
    # The new Kd is plausible, so the molecular-SO₂ readout tier is unmoved (still gated by the
    # pH-solver pKa tiers + the binding constants — all plausible).
    tier = acidbase.molecular_so2_tier(pset.tier_map())
    assert tier is Tier.PLAUSIBLE
