"""Malolactic fermentation v1 — conversion-only (decision D-23).

Ranked headline-first. The keystone payoff is ``test_headline_mlf_raises_ph_emergently``:
the hand-built ``test_acidbase.test_headline_malic_to_lactic_raises_ph`` ΔpH ∈ [0.1, 0.3]
now *emerges* from the :class:`MalolacticConversion` Process on a malic-rich co-inoculated
must — measured as the no-MLF **control difference** ``pH_final(dosed) − pH_final(off)`` so
it isolates the malic→lactic swap from the (shared) Byp pH down-drift. The rest pin the
mole-for-mole stoichiometry + carbon closure, the gates (pH up ⇒ faster, ethanol/molecular
SO₂ ⇒ arrest — the first RHS consumers of D-18 pH and D-22 SO₂), the cardinal-temperature
optimum, the explicit ``speculative`` tier, and — prime directive #3 — isolability: an
undosed run is byte-for-byte the validated core *and* keeps ``malic``/``lactic`` at the
VALIDATED tier (the Process is disabled at compile when MLF is not pitched).
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.analysis import ph_series
from fermentation.core import acidbase
from fermentation.core.chemistry import (
    M_CO2,
    M_LACTIC,
    M_MALIC,
    M_TARTARIC,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.malolactic import (
    MalolacticConversion,
    cardinal_temperature_factor,
)
from fermentation.core.media import wine_schema
from fermentation.core.state import StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon


@pytest.fixture
def pset():
    """Real wine kinetic params (incl. the MLF set) + the shared pKa set."""
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    return pset.resolve()


def _anchor_cation(pka, tartaric_gpl: float, malic_gpl: float, lactic_gpl: float, target_ph: float):
    totals = {
        "tartaric": tartaric_gpl / M_TARTARIC,
        "malic": malic_gpl / M_MALIC,
        "lactic": lactic_gpl / M_LACTIC,
    }
    return acidbase.solve_cation_charge(totals, 0.0, pka, target_ph)


def _wine_state(schema: StateSchema, params: Mapping[str, float], *, target_ph: float, **slots):
    """A wine state with the bulk vars, given acid/MLF slots, and the cation that anchors
    ``target_ph`` from the named acids (so pH gates see a realistic pH)."""
    tartaric = slots.get("tartaric", 0.0)
    malic = slots.get("malic", 0.0)
    lactic = slots.get("lactic", 0.0)
    cation = _anchor_cation(acidbase.build_pka_map(params), tartaric, malic, lactic, target_ph)
    base: dict[str, float | list[float]] = {
        "X": 1.0, "S": [120.0], "E": 0.0, "N": 0.1, "T": 293.15, "CO2": 0.0,
        "cation_charge": cation,
    }  # fmt: skip
    base.update(slots)
    return schema.pack(base)


def _wine_scenario(**initial_extra) -> Scenario:
    initial: dict[str, float] = {
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.5,
        "tartaric_gpl": 4.0,
        "malic_gpl": 4.0,
        "initial_ph": 3.4,
    }
    initial.update(initial_extra)
    return Scenario(
        name="wine-mlf",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )


def _run(**initial_extra: float):
    compiled = compile_scenario(_wine_scenario(**initial_extra), strict=True)
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    return compiled, traj


# -- 1. HEADLINE: MLF deacidification emerges, ΔpH 0.1–0.3 (no-MLF control diff) ----


def test_headline_mlf_raises_ph_emergently():
    # The faithful emergent translation of the hand-built acidbase headline: dose O. oeni
    # on a malic-rich must and the pH rises into the deacidification band — measured as the
    # difference vs an otherwise-identical un-pitched run, so the (shared) Byp down-drift
    # cancels and the gap is purely the malic→lactic swap at the same final state.
    c_off, t_off = _run()
    c_on, t_on = _run(mlf_pitch_gpl=0.2)
    ph_off = ph_series(t_off, c_off.param_values)
    ph_on = ph_series(t_on, c_on.param_values)

    # MLF runs (co-inoculation) to near-completion in the early low-ethanol window
    assert t_off.series("malic")[-1] == pytest.approx(4.0)  # disabled when un-pitched
    assert t_on.series("malic")[-1] < 0.05  # malate essentially consumed
    assert t_on.series("lactic")[-1] > 2.0  # lactate produced

    delta = ph_on[-1] - ph_off[-1]
    assert 0.1 <= delta <= 0.3, f"emergent MLF ΔpH {delta:.3f} outside [0.1, 0.3]"


# -- 2. mole-for-mole stoichiometry + carbon closure at the RHS level --------------


def test_conversion_is_mole_for_mole_and_carbon_closing(params):
    schema = wine_schema()
    y = _wine_state(schema, params, target_ph=3.4, malic=4.0, tartaric=4.0, X_mlf=0.2)
    d = MalolacticConversion().derivatives(0.0, y, schema, params)
    dmalic = float(d[schema.slice("malic")][0])
    dlactic = float(d[schema.slice("lactic")][0])
    dco2 = float(d[schema.slice("CO2")][0])

    assert dmalic < 0.0 and dlactic > 0.0 and dco2 > 0.0
    # one mole of malate makes one mole of lactate and one of CO2
    r = -dmalic / M_MALIC
    assert dlactic / M_LACTIC == pytest.approx(r)
    assert dco2 / M_CO2 == pytest.approx(r)
    # carbon closes on the existing ledger (4 C = 3 C + 1 C): weighted sum is zero
    carbon_rate = (
        dmalic * carbon_mass_fraction("malic_acid")
        + dlactic * carbon_mass_fraction("lactic_acid")
        + dco2 * carbon_mass_fraction("CO2")
    )
    assert carbon_rate == pytest.approx(0.0, abs=1e-12)


def test_carbon_conserved_over_a_dosed_run():
    compiled, traj = _run(mlf_pitch_gpl=0.2)
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.parameters["biomass_C_fraction"].value
    )
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-9, label="total carbon (MLF on)")
    # the conversion self-limits (rate ∝ malate, guarded ≤0), so malate/lactate never go
    # meaningfully negative even as malate is driven to ~0
    assert_nonnegative(traj, ("malic", "lactic", "CO2"))


# -- 3. touches only malic/lactic/CO2 (strict contract) ---------------------------


def test_touches_only_malic_lactic_co2(params):
    schema = wine_schema()
    y = _wine_state(schema, params, target_ph=3.4, malic=4.0, tartaric=4.0, X_mlf=0.2)
    d = MalolacticConversion().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"malic", "lactic", "CO2"}


# -- 4. pH gate: higher pH ⇒ faster (the deacidification feedback) -----------------


def test_higher_ph_speeds_conversion(params):
    schema = wine_schema()
    # same malate, two anchored pHs: the rate rises with pH (O. oeni pH gate)
    y_lo = _wine_state(schema, params, target_ph=3.0, malic=4.0, tartaric=4.0, X_mlf=0.2)
    y_hi = _wine_state(schema, params, target_ph=3.6, malic=4.0, tartaric=4.0, X_mlf=0.2)
    rate_lo = -float(
        MalolacticConversion().derivatives(0.0, y_lo, schema, params)[schema.slice("malic")][0]
    )
    rate_hi = -float(
        MalolacticConversion().derivatives(0.0, y_hi, schema, params)[schema.slice("malic")][0]
    )
    assert 0.0 < rate_lo < rate_hi


# -- 5. ethanol gate: above O. oeni tolerance, MLF stops --------------------------


def test_ethanol_above_tolerance_arrests_mlf(params):
    schema = wine_schema()
    above = params["ethanol_tolerance_mlf"] + 5.0
    y = _wine_state(schema, params, target_ph=3.4, malic=4.0, tartaric=4.0, X_mlf=0.2, E=above)
    d = MalolacticConversion().derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("malic")][0]) == 0.0


# -- 6. molecular-SO₂ gate: dosing free SO₂ suppresses MLF (first RHS consumer) ----


def test_so2_dose_suppresses_mlf_unit(params):
    schema = wine_schema()
    from fermentation.units.convert import mgl_to_gpl

    common = {"target_ph": 3.4, "malic": 4.0, "tartaric": 4.0, "X_mlf": 0.2}
    y_clean = _wine_state(schema, params, **common)
    y_so2 = _wine_state(schema, params, so2_free=mgl_to_gpl(80.0), **common)
    rate_clean = -float(
        MalolacticConversion().derivatives(0.0, y_clean, schema, params)[schema.slice("malic")][0]
    )
    rate_so2 = -float(
        MalolacticConversion().derivatives(0.0, y_so2, schema, params)[schema.slice("malic")][0]
    )
    assert rate_clean > 0.0
    assert rate_so2 < 0.01 * rate_clean  # molecular SO₂ arrests the conversion


def test_so2_dose_suppresses_mlf_in_a_run():
    # The integration-level demonstration: SO₂ leaves pH/carbon untouched (D-22 readout-only)
    # yet, as the first RHS consumer of molecular SO₂, it gates MLF — so the dosed-SO₂ run
    # barely deacidifies relative to the same MLF dose without SO₂.
    c_off, t_off = _run()
    c_on, t_on = _run(mlf_pitch_gpl=0.2)
    c_so2, t_so2 = _run(mlf_pitch_gpl=0.2, so2_free_mgl=80.0)
    ph_off = ph_series(t_off, c_off.param_values)
    ph_on = ph_series(t_on, c_on.param_values)
    ph_so2 = ph_series(t_so2, c_so2.param_values)

    assert t_so2.series("malic")[-1] > 3.9  # MLF essentially blocked
    assert (ph_on[-1] - ph_off[-1]) > 0.1  # uninhibited MLF deacidifies
    assert (ph_so2[-1] - ph_off[-1]) == pytest.approx(0.0, abs=0.01)  # SO₂ arrests it


# -- 7. cardinal-temperature optimum (peaks at T_opt, declines warm, 0 outside) ----


def test_cardinal_temperature_factor_shape(params):
    t_min, t_opt, t_max = params["T_min_mlf"], params["T_opt_mlf"], params["T_max_mlf"]
    assert cardinal_temperature_factor(t_opt, t_min, t_opt, t_max) == pytest.approx(1.0)
    assert cardinal_temperature_factor(t_min, t_min, t_opt, t_max) == 0.0
    assert cardinal_temperature_factor(t_max, t_min, t_opt, t_max) == 0.0
    assert cardinal_temperature_factor(t_min - 5.0, t_min, t_opt, t_max) == 0.0
    assert cardinal_temperature_factor(t_max + 5.0, t_min, t_opt, t_max) == 0.0
    # peak at the optimum, in (0, 1) between cardinals, and DECLINES on the warm side
    # (the qualitative correctness a monotone Arrhenius factor lacks — decision D-23)
    warm = cardinal_temperature_factor(t_opt + 7.0, t_min, t_opt, t_max)
    cool = cardinal_temperature_factor(t_opt - 7.0, t_min, t_opt, t_max)
    assert 0.0 < warm < 1.0 and 0.0 < cool < 1.0
    grid = np.linspace(t_min + 0.1, t_max - 0.1, 200)
    vals = np.array([cardinal_temperature_factor(t, t_min, t_opt, t_max) for t in grid])
    assert vals.max() == pytest.approx(1.0, abs=1e-3)
    assert np.argmax(vals) == int(np.argmin(np.abs(grid - t_opt)))  # unimodal peak at T_opt


# -- 8. isolability (prime directive #3): undosed = validated core, byte + tier ----


def test_undosed_run_is_byte_for_byte_validated_core():
    # mlf_pitch_gpl absent vs explicitly 0 must integrate identically (the Process is
    # disabled both ways), so MLF adds nothing to an un-pitched wine ferment.
    t_eval = np.linspace(0.0, 21.0 * 24.0, 80)
    c_absent = compile_scenario(_wine_scenario(), strict=True)
    c_zero = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.0), strict=True)
    ta = simulate(
        c_absent.process_set, c_absent.param_values, c_absent.y0, c_absent.t_span_h, t_eval=t_eval
    )
    tz = simulate(
        c_zero.process_set, c_zero.param_values, c_zero.y0, c_zero.t_span_h, t_eval=t_eval
    )
    for name in c_absent.schema.names:
        assert np.allclose(ta.series(name), tz.series(name), rtol=1e-12, atol=1e-12), name
    assert not c_absent.process_set.is_enabled(MalolacticConversion.name)


def test_undosed_keeps_malic_lactic_validated_dosed_makes_speculative(pset):
    tm = pset.tier_map()
    off = compile_scenario(_wine_scenario()).process_set
    on = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.2)).process_set
    # un-pitched: nothing active touches the acid slots, so they stay VALIDATED
    assert off.tier_of("malic", tm) is Tier.VALIDATED
    assert off.tier_of("lactic", tm) is Tier.VALIDATED
    # pitched: the speculative MLF Process touches them ⇒ speculative (prime directive #1)
    assert on.tier_of("malic", tm) is Tier.SPECULATIVE
    assert on.tier_of("lactic", tm) is Tier.SPECULATIVE


def test_x_mlf_is_inert_over_a_run():
    # v1: no Process grows or kills the catalyst, so its concentration is constant.
    compiled, traj = _run(mlf_pitch_gpl=0.2)
    x_mlf = traj.series("X_mlf")
    assert np.allclose(x_mlf, x_mlf[0], rtol=0.0, atol=1e-12)
    assert x_mlf[0] == pytest.approx(0.2)


# -- 9. tier is speculative -------------------------------------------------------


def test_process_tier_is_speculative():
    assert MalolacticConversion.tier is Tier.SPECULATIVE
