"""MLF medium-chain-fatty-acid inhibition — the yeast-secreted stuck-MLF cause (decision D-131).

The MLF environmental gate was ``g_pH·g_EtOH·g_SO₂·γ(T)`` (D-23→D-39). Lonvaud-Funel, Joyeux &
Desens (1988, *J. Sci. Food Agric.* 44:183) show *S. cerevisiae* secretes medium-chain fatty
acids (octanoic + decanoic — hexanoic is inert) during AF that inhibit the malolactic activity of
*O. oeni*, the mechanism behind the empirical fact that MLF fermentability of wines from the same
must varies with the AF yeast strain. This adds a fifth, **bacteriostatic** factor
``g_FA = exp(−[MCFA]/mcfa_inhib_mlf)`` reading a new inert wine-composition-at-MLF input slot
(``mcfa``, ``mcfa_mgl`` scenario key, octanoic-equivalent) — the D-130 oxofructose pattern.

Two headline gates (the D-129 discriminating-check pattern):
  * GATE 1 — ``test_gate1_isolable_byte_for_byte_at_zero_load``: at ``mcfa = 0`` ``g_FA = 1``
    EXACTLY, so the environmental gate is byte-for-byte the pre-D-131 ``toxicity·γ(T)`` and every
    no-MCFA MLF run is unchanged (structural, not param-luck).
  * GATE 2 — ``test_gate2_stressed_mcfa_collapses_conversion``: a stressed-AF MCFA load
    (~14.4 mg/L, Lonvaud-Funel Table 6 FA2) collapses the MLF conversion rate ~60% — the
    emergent sluggish/stuck MLF.

The gate lives in the conversion/growth gate and NOT the death term (the D-39 lesson — MCFA is a
"can't convert" signal, not a lethal one), pinned by ``test_mcfa_does_not_affect_death`` and its
structural companion ``test_death_reads_exclude_mcfa_param``.
"""

from collections.abc import Mapping

import pytest

from fermentation.core import acidbase
from fermentation.core.chemistry import (
    CARBON_ATOMS,
    M_MALIC,
    M_OCTANOIC_ACID,
    M_TARTARIC,
    NITROGEN_ATOMS,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.malolactic import (
    MalolacticConversion,
    MalolacticDeath,
    MalolacticGrowth,
    cardinal_temperature_factor,
    malolactic_environmental_gate,
    malolactic_fatty_acid_gate,
    malolactic_toxicity_gate,
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


def _octanoic_equiv_mgl(micromolar: float) -> float:
    """octanoic-equivalent mg/L for a given µmol/L C8+C10 total (the Table-6 unit)."""
    return micromolar * 1e-6 * M_OCTANOIC_ACID * 1e3


def _anchor_cation(params: Mapping[str, float], tartaric: float, malic: float, ph: float) -> float:
    totals = {"tartaric": tartaric / M_TARTARIC, "malic": malic / M_MALIC, "lactic": 0.0}
    return acidbase.solve_cation_charge(totals, 0.0, acidbase.build_pka_map(params), ph)


def _mlf_state(
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    mcfa_mgl: float,
    so2_mgl: float = 0.0,
    ph: float = 3.4,
) -> FloatArray:
    """A co-inoculated MLF state: O. oeni pitched, malic present, amino-acid fuel + sugar for
    growth, the given MCFA (octanoic-equiv) and optional SO₂ dose, cation anchoring ``ph``."""
    tartaric, malic = 4.0, 4.0
    cation = _anchor_cation(params, tartaric, malic, ph)
    base: dict[str, float | list[float]] = {
        "X": 1.0, "S": [80.0], "E": 40.0, "N": 0.05, "T": 293.15, "CO2": 0.0,
        "cation_charge": cation, "tartaric": tartaric, "malic": malic,
        "X_mlf": 0.2, "amino_acids": 1.0, "amino_acids_generic": 1.0,
        "mcfa": mgl_to_gpl(mcfa_mgl), "so2_total": mgl_to_gpl(so2_mgl),
    }  # fmt: skip
    return schema.pack(base)


def _malic_turnover(schema, params, **kw) -> float:
    d = MalolacticConversion().derivatives(0.0, _mlf_state(schema, params, **kw), schema, params)
    return float(-d[schema.slice("malic")][0])  # +malate turnover [g/L/h]


# -- 1. Provenance / metadata -------------------------------------------------


def test_mcfa_inhib_param_present_and_plausible(pset):
    p = pset["mcfa_inhib_mlf"]
    assert p.value == pytest.approx(1.09e-4)  # scale calibrated to Lonvaud-Funel Table 6 FA2
    assert p.unit == "mol/L"
    assert p.tier is Tier.PLAUSIBLE  # real-wine measured activity-vs-MCFA


def test_chemistry_registry_octanoic_c8_nitrogen_free():
    # C8H16O2: mass ~144.21, eight carbons, nitrogen-free (a fatty acid).
    assert pytest.approx(144.21, abs=0.05) == M_OCTANOIC_ACID
    assert CARBON_ATOMS["octanoic_acid"] == 8
    assert NITROGEN_ATOMS["octanoic_acid"] == 0
    assert carbon_mass_fraction("octanoic_acid") == pytest.approx(8 * 12.011 / M_OCTANOIC_ACID)


def test_wine_only_absent_from_beer():
    # Yeast MCFA inhibiting O. oeni is a wine-MLF concept — wine-only, like the acids/SO₂.
    assert "mcfa" in wine_schema()
    assert "mcfa" not in beer_schema()


# -- 2. GATE 1: structural isolability ----------------------------------------


def test_gate1_isolable_byte_for_byte_at_zero_load(params):
    # At mcfa = 0 the g_FA factor is EXACTLY 1, so the environmental gate is byte-for-byte the
    # pre-D-131 toxicity·γ(T). Verified across a sweep of pH/ethanol/SO₂ backgrounds (structural,
    # not param-luck) against an INDEPENDENT hand-computed FA-stripped gate.
    schema = wine_schema()
    max_gap = 0.0
    for ph, so2 in ((3.2, 0.0), (3.4, 20.0), (3.6, 40.0), (3.8, 0.0)):
        y = _mlf_state(schema, params, mcfa_mgl=0.0, so2_mgl=so2, ph=ph)
        assert malolactic_fatty_acid_gate(y, schema, params) == 1.0  # exact
        solved = acidbase.ph_of_state(y, schema, params)
        env = malolactic_environmental_gate(y, schema, params, solved)
        gamma = cardinal_temperature_factor(
            293.15, params["T_min_mlf"], params["T_opt_mlf"], params["T_max_mlf"]
        )
        env_no_fa = malolactic_toxicity_gate(y, schema, params, solved) * gamma
        max_gap = max(max_gap, abs(env - env_no_fa))
    assert max_gap == 0.0  # exact — no fifth factor contributes at load 0


# -- 3. GATE 2: the sourced stuck-MLF outcome ---------------------------------


def test_gate2_stressed_mcfa_collapses_conversion(params):
    # THE headline: a stressed-AF MCFA load (Lonvaud-Funel Table 6 FA2 = 100 µM octanoic+decanoic
    # ~ 14.4 mg/L octanoic-equiv) collapses the malolactic conversion rate ~60% — the emergent
    # sluggish/stuck MLF a high-MCFA yeast strain leaves behind.
    schema = wine_schema()
    r0 = _malic_turnover(schema, params, mcfa_mgl=0.0)
    r_stressed = _malic_turnover(schema, params, mcfa_mgl=_octanoic_equiv_mgl(100.0))
    assert r0 > 0.0
    assert r_stressed == pytest.approx(0.40 * r0, rel=0.05)  # ~60% activity loss (Table 6 FA2)


def test_gate_calibrated_to_table6_band(params):
    # g_FA at the three Lonvaud-Funel Table 6 mixtures reproduces the measured activity-loss band.
    schema = wine_schema()
    for micromolar, measured in ((20.0, 0.87), (100.0, 0.40), (150.0, 0.20)):
        y = _mlf_state(schema, params, mcfa_mgl=_octanoic_equiv_mgl(micromolar))
        g = malolactic_fatty_acid_gate(y, schema, params)
        assert g == pytest.approx(measured, abs=0.06)  # within biological scatter


def test_conversion_rate_monotone_in_mcfa(params):
    # More MCFA ⇒ monotonically slower conversion; g_FA ∈ (0, 1] never negative or > 1.
    schema = wine_schema()
    loads = [_octanoic_equiv_mgl(u) for u in (0.0, 20.0, 50.0, 100.0, 150.0, 220.0)]
    rates = [_malic_turnover(schema, params, mcfa_mgl=m) for m in loads]
    assert all(a > b for a, b in zip(rates, rates[1:], strict=False))
    for m in loads:
        g = malolactic_fatty_acid_gate(_mlf_state(schema, params, mcfa_mgl=m), schema, params)
        assert 0.0 < g <= 1.0


def test_mcfa_gates_growth_too(params):
    # The shared environmental gate throttles MalolacticGrowth as well as conversion — Lonvaud-Funel
    # Table 3 confirms the DIRECTION (the MCFA mixture inhibits bacterial growth too), but Table 3's
    # growth dose-response is NOT Table 6's resting-cell activity response: applying the same
    # activity-calibrated g_FA to growth is the v1 shared-gate MODELING choice (the pH/EtOH/SO₂
    # gates already act on both), not a sourced growth==activity equivalence. Here we pin that the
    # growth rate carries EXACTLY that g_FA factor (the shared-gate wiring), not some other path.
    schema = wine_schema()
    growth = MalolacticGrowth()
    dx0 = float(
        growth.derivatives(0.0, _mlf_state(schema, params, mcfa_mgl=0.0), schema, params)[
            schema.slice("X_mlf")
        ][0]
    )
    stressed = _octanoic_equiv_mgl(100.0)
    dxs = float(
        growth.derivatives(0.0, _mlf_state(schema, params, mcfa_mgl=stressed), schema, params)[
            schema.slice("X_mlf")
        ][0]
    )
    assert dx0 > 0.0
    g = malolactic_fatty_acid_gate(_mlf_state(schema, params, mcfa_mgl=stressed), schema, params)
    assert dxs == pytest.approx(g * dx0, rel=1e-9)  # exactly the g_FA factor, not some other path


# -- 4. The D-39 lesson: bacteriostatic, not a death driver -------------------


def test_mcfa_does_not_affect_death(params):
    # MCFA is bacteriostatic ("can't convert"), NOT lethal — MalolacticDeath (SO₂-driven, D-39)
    # is byte-for-byte independent of the MCFA load; death still fires on the SO₂ dose.
    schema = wine_schema()
    death = MalolacticDeath()
    rd0 = float(
        -death.derivatives(
            0.0, _mlf_state(schema, params, mcfa_mgl=0.0, so2_mgl=30.0), schema, params
        )[schema.slice("X_mlf")][0]
    )
    rdf = float(
        -death.derivatives(
            0.0,
            _mlf_state(schema, params, mcfa_mgl=_octanoic_equiv_mgl(150.0), so2_mgl=30.0),
            schema,
            params,
        )[schema.slice("X_mlf")][0]
    )
    assert rd0 > 0.0  # death is real (SO₂ dosed)
    assert rd0 == rdf  # exact — MCFA never touches the death rate


def test_death_reads_exclude_mcfa_param_but_conversion_includes_it():
    # Structural companion to the numerical death-independence: the fatty-acid parameter is in the
    # conversion/growth reads (it gates them) but NOT in the death/senescence reads (they can't be
    # accelerated by MCFA) — so a future rewire that leaked it into death would fail loudly here.
    assert "mcfa_inhib_mlf" in MalolacticConversion.reads
    assert "mcfa_inhib_mlf" in MalolacticGrowth.reads
    assert "mcfa_inhib_mlf" not in MalolacticDeath.reads


# -- 5. Conservation + wiring -------------------------------------------------


def test_inert_pool_is_carbon_conservation_neutral(params):
    # The mcfa pool is weighted in total_carbon (C8) but INERT (no Process touches it), so it is a
    # constant term: two states differing only in the dosed load differ in total_carbon by exactly
    # the pool's carbon, and it has no derivative (drifts 0).
    from fermentation.validation import total_carbon

    schema = wine_schema()
    tc = total_carbon(schema, biomass_carbon_fraction=0.48)  # cancels in the delta (shared X)
    lo = _mlf_state(schema, params, mcfa_mgl=0.0)
    hi = _mlf_state(schema, params, mcfa_mgl=15.0)
    delta = tc(hi) - tc(lo)
    expected = mgl_to_gpl(15.0) * carbon_mass_fraction("octanoic_acid")
    assert delta == pytest.approx(expected, rel=1e-9)


def test_scenario_wires_mcfa_mgl():
    # mcfa_mgl (mg/L) → the mcfa slot (g/L); absent ⇒ 0 (byte-for-byte the pre-D-131 MLF).
    base = {
        "brix": 24.0,
        "pitch_gpl": 0.5,
        "yan_mgl": 250.0,
        "malic_gpl": 4.0,
        "mlf_pitch_gpl": 0.2,
    }
    stressed = Scenario(
        name="stressed-af",
        medium="wine",
        initial={**base, "mcfa_mgl": 15.0},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=1.0,
    )
    clean = Scenario(
        name="clean-af",
        medium="wine",
        initial=base,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=1.0,
    )
    c_stressed = compile_scenario(stressed)
    c_clean = compile_scenario(clean)
    assert float(c_stressed.y0[c_stressed.schema.slice("mcfa")][0]) == pytest.approx(
        mgl_to_gpl(15.0)
    )
    assert float(c_clean.y0[c_clean.schema.slice("mcfa")][0]) == 0.0


def test_unknown_key_still_fails_loudly():
    # The allow-list guard is intact — a typo near the new key is rejected, not silently ignored.
    bad = Scenario(
        name="typo",
        medium="wine",
        initial={"brix": 24.0, "mcfa_mg": 15.0},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=1.0,
    )
    with pytest.raises((ValueError, KeyError)):
        compile_scenario(bad)


def test_mlf_output_tier_stays_speculative():
    # The MCFA scale is plausible, but it does NOT lift the MLF output tier: MalolacticConversion
    # is speculative (k_mlf), so malic/lactic stay speculative via Tier.combine (the advisor's
    # "tier won't move MLF"). The plausible param buys provenance, not a tier promotion.
    assert MalolacticConversion.tier is Tier.SPECULATIVE
    assert MalolacticGrowth.tier is Tier.SPECULATIVE
