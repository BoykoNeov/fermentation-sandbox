"""MLF-derived diacetyl (decision D-31): O. oeni citrate co-metabolism + bacterial reduction.

Headline-first. The payoff is ``test_headline_citrate_lifts_and_then_clears_diacetyl``: dosing
O. oeni *and* citrate makes wine diacetyl rise clearly above the yeast-only baseline, peak
*late* (post the early low-ethanol window, via the shared VDK reservoir's decarb lag), then fall
as reduction clears it — the buttery-then-cleaning-up MLF signature, emergent from the existing
D-26 machinery. The rest pin the lumped citrate → α-acetolactate + CO2 stoichiometry and its
carbon closure (6 = 5 + 1, the reason a *citrate* pool is needed — sugar carbon no-ops at
dryness), the shared O. oeni environmental gate, the bacterial reducer, the ``speculative`` tier,
and — prime directive #3 — isolability: an un-pitched (or citrate-free) run is byte-for-byte the
prior core and keeps ``citrate`` at the VALIDATED tier.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.chemistry import (
    M_ACETOLACTATE,
    M_CITRIC,
    M_CO2,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.malolactic import (
    MalolacticCitrateMetabolism,
    OenococcusDiacetylReduction,
)
from fermentation.core.media import wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon


@pytest.fixture
def pset():
    data = default_data_dir()
    return load_parameters(
        data / "wine_generic.yaml", data / "acidbase.yaml", data / "vicinal_diketones.yaml"
    )


@pytest.fixture
def params(pset):
    return pset.resolve()


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
        name="wine-mlf-diacetyl",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=30.0,
    )


def _run(**initial_extra: float):
    compiled = compile_scenario(_wine_scenario(**initial_extra), strict=True)
    t_eval = np.linspace(0.0, 30.0 * 24.0, 600)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    return compiled, traj


def _wine_state(schema: StateSchema, params: Mapping[str, float], **slots: float) -> FloatArray:
    """A wine state with sane bulk vars plus given acid/MLF/citrate slots (no pH anchoring
    needed — these RHS-level tests read the citrate/reservoir terms, and the pH gate is
    exercised via the run-level and gate-direction tests)."""
    base: dict[str, float | list[float]] = {
        "X": 1.0, "S": [120.0], "E": 0.0, "N": 0.1, "T": 293.15, "CO2": 0.0,
    }  # fmt: skip
    base.update(slots)
    return schema.pack(base)


# -- 1. HEADLINE: citrate lifts diacetyl above the yeast-only baseline, then clears -----


def test_headline_citrate_lifts_and_then_clears_diacetyl():
    # Yeast-only diacetyl (no O. oeni) vs O. oeni + citrate. The MLF-derived branch adds a
    # clearly-larger diacetyl pool that peaks LATE and then falls (reduction clears it).
    _, t_base = _run()  # MLF disabled at the compile seam (un-pitched)
    _, t_mlf = _run(mlf_pitch_gpl=0.2, citrate_gpl=0.3)

    dia_base = t_base.series("diacetyl") * 1000.0  # mg/L
    dia_mlf = t_mlf.series("diacetyl") * 1000.0

    peak_base = float(dia_base.max())
    peak_mlf = float(dia_mlf.max())
    peak_day = float(t_mlf.t[int(dia_mlf.argmax())] / 24.0)

    # MLF-derived diacetyl is a clear lever, not a rounding bump: at least ~2× the baseline
    assert peak_mlf > 2.0 * peak_base
    assert peak_mlf > 0.2  # lands in the buttery range (threshold ~0.2 mg/L wine)
    # …and it is a LATE peak (after the early low-ethanol conversion window), not an initial spike
    assert peak_day > 3.0
    # …then it clears: the final diacetyl is meaningfully below its own peak (reduction wins).
    # MEASURED ratio ~0.742 (D-53-corrected: X_mlf stays ~98.6% viable by day 30 — a real-wine
    # literature check found no support for D-52's faster senescence, so k_senescence_mlf dropped
    # ~50x — plenty of bacterial reductase remains on the lees to keep clearing diacetyl, closer to
    # D-41's original clean-clearing picture than D-52's transient 0.861 measurement).
    assert float(dia_mlf[-1]) < 0.80 * peak_mlf


def test_citrate_is_mostly_unconsumed_the_trace_branch():
    # The lumped stand-in is honest only because citrate is barely touched (the dominant
    # acetate branch is untracked, so full depletion would be a fiction). Assert <~20% used.
    _, t_mlf = _run(mlf_pitch_gpl=0.2, citrate_gpl=0.3)
    citrate = t_mlf.series("citrate")
    consumed_frac = (citrate[0] - citrate[-1]) / citrate[0]
    assert 0.0 < consumed_frac < 0.2


# -- 2. lumped stoichiometry + carbon closure (citrate C6 → acetolactate C5 + CO2 C1) ---


def test_citrate_metabolism_is_carbon_closing_at_the_rhs(params):
    schema = wine_schema()
    y = _wine_state(schema, params, citrate=0.3, X_mlf=0.2)
    d = MalolacticCitrateMetabolism().derivatives(0.0, y, schema, params)
    dcit = float(d[schema.slice("citrate")][0])
    dace = float(d[schema.slice("acetolactate")][0])
    dco2 = float(d[schema.slice("CO2")][0])

    assert dcit < 0.0 and dace > 0.0 and dco2 > 0.0
    # one mole of citrate makes one mole of α-acetolactate and one of CO2 (lumped stand-in)
    r = -dcit / M_CITRIC
    assert dace / M_ACETOLACTATE == pytest.approx(r)
    assert dco2 / M_CO2 == pytest.approx(r)
    # carbon closes on the existing ledger (6 C = 5 C + 1 C): weighted sum is zero
    carbon_rate = (
        dcit * carbon_mass_fraction("citric_acid")
        + dace * carbon_mass_fraction("alpha_acetolactate")
        + dco2 * carbon_mass_fraction("CO2")
    )
    assert carbon_rate == pytest.approx(0.0, abs=1e-12)


def test_carbon_conserved_over_a_dosed_run():
    compiled, traj = _run(mlf_pitch_gpl=0.2, citrate_gpl=0.3)
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.parameters["biomass_C_fraction"].value
    )
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-9, label="total carbon (MLF-diacetyl on)")
    assert_nonnegative(traj, ("citrate", "acetolactate", "diacetyl", "butanediol", "CO2"))


# -- 3. touches contracts (strict) ------------------------------------------------------


def test_citrate_metabolism_touches_only_citrate_acetolactate_co2(params):
    schema = wine_schema()
    y = _wine_state(schema, params, citrate=0.3, X_mlf=0.2)
    d = MalolacticCitrateMetabolism().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"citrate", "acetolactate", "CO2"}


def test_bacterial_reduction_touches_only_diacetyl_butanediol(params):
    schema = wine_schema()
    y = _wine_state(schema, params, diacetyl=0.0005, X_mlf=0.2)
    d = OenococcusDiacetylReduction().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"diacetyl", "butanediol"}


# -- 4. gates & catalyst dependence -----------------------------------------------------


def test_citrate_metabolism_needs_catalyst_and_substrate(params):
    schema = wine_schema()
    # no O. oeni ⇒ zero (structural value-isolability, no pH solve paid)
    no_cat = MalolacticCitrateMetabolism().derivatives(
        0.0, _wine_state(schema, params, citrate=0.3, X_mlf=0.0), schema, params
    )
    assert float(no_cat[schema.slice("citrate")][0]) == 0.0
    # O. oeni but no citrate ⇒ zero
    no_sub = MalolacticCitrateMetabolism().derivatives(
        0.0, _wine_state(schema, params, citrate=0.0, X_mlf=0.2), schema, params
    )
    assert float(no_sub[schema.slice("citrate")][0]) == 0.0


def test_ethanol_above_tolerance_arrests_citrate_metabolism(params):
    schema = wine_schema()
    above = params["ethanol_tolerance_mlf"] + 5.0
    y = _wine_state(schema, params, citrate=0.3, X_mlf=0.2, E=above)
    d = MalolacticCitrateMetabolism().derivatives(0.0, y, schema, params)
    # shares the O. oeni environmental gate with the malate conversion: ethanol wall ⇒ 0
    assert float(d[schema.slice("citrate")][0]) == 0.0


def test_bacterial_reduction_scales_with_catalyst_and_diacetyl(params):
    schema = wine_schema()
    proc = OenococcusDiacetylReduction()
    # no diacetyl ⇒ nothing to reduce; no catalyst ⇒ no bacterial reduction
    assert (
        float(
            proc.derivatives(
                0.0, _wine_state(schema, params, diacetyl=0.0, X_mlf=0.2), schema, params
            )[schema.slice("diacetyl")][0]
        )
        == 0.0
    )
    assert (
        float(
            proc.derivatives(
                0.0, _wine_state(schema, params, diacetyl=0.0005, X_mlf=0.0), schema, params
            )[schema.slice("diacetyl")][0]
        )
        == 0.0
    )
    # with both present the reduction is a loss (diacetyl falls, butanediol rises)
    d = proc.derivatives(
        0.0, _wine_state(schema, params, diacetyl=0.0005, X_mlf=0.2), schema, params
    )
    assert float(d[schema.slice("diacetyl")][0]) < 0.0
    assert float(d[schema.slice("butanediol")][0]) > 0.0


def test_more_oenococcus_clears_diacetyl_faster():
    # The bacterial reducer matters: a larger O. oeni dose leaves less residual diacetyl at
    # the end (it keeps clearing after the yeast is ethanol-inactivated — lees-contact clean-up).
    _, t_lo = _run(mlf_pitch_gpl=0.1, citrate_gpl=0.3)
    _, t_hi = _run(mlf_pitch_gpl=0.4, citrate_gpl=0.3)
    # normalise by the citrate actually consumed (more bacteria also make more diacetyl), so
    # this isolates the *reduction* effect: final/peak ratio is lower with more bacteria.
    ratio_lo = t_lo.series("diacetyl")[-1] / t_lo.series("diacetyl").max()
    ratio_hi = t_hi.series("diacetyl")[-1] / t_hi.series("diacetyl").max()
    assert ratio_hi < ratio_lo


# -- 5. isolability (prime directive #3) ------------------------------------------------


def test_undosed_run_is_byte_for_byte_prior_core():
    # No O. oeni and no citrate ⇒ the two D-31 Processes are disabled at the compile seam, so
    # the run is identical to one that names neither key (the yeast-pathway diacetyl core).
    t_eval = np.linspace(0.0, 30.0 * 24.0, 80)
    c_absent = compile_scenario(_wine_scenario(), strict=True)
    c_zero = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.0, citrate_gpl=0.0), strict=True)
    ta = simulate(
        c_absent.process_set, c_absent.param_values, c_absent.y0, c_absent.t_span_h, t_eval=t_eval
    )
    tz = simulate(
        c_zero.process_set, c_zero.param_values, c_zero.y0, c_zero.t_span_h, t_eval=t_eval
    )
    for name in c_absent.schema.names:
        assert np.allclose(ta.series(name), tz.series(name), rtol=1e-12, atol=1e-12), name
    assert not c_absent.process_set.is_enabled(MalolacticCitrateMetabolism.name)
    assert not c_absent.process_set.is_enabled(OenococcusDiacetylReduction.name)


def test_citrate_dosed_without_oenococcus_is_inert():
    # Citrate is a real carbon-bearing must input, but with no O. oeni nothing metabolises it,
    # so it sits constant and (via the disabled Processes) perturbs no other column. Diacetyl
    # matches the citrate-free yeast-only baseline.
    _, t_base = _run()
    _, t_cit = _run(citrate_gpl=0.3)  # citrate but no mlf pitch
    assert np.allclose(t_cit.series("citrate"), 0.3, atol=1e-12)
    assert np.allclose(t_cit.series("diacetyl"), t_base.series("diacetyl"), rtol=1e-9, atol=1e-12)


def test_undosed_keeps_citrate_validated_dosed_makes_speculative(pset):
    tm = pset.tier_map()
    off = compile_scenario(_wine_scenario(citrate_gpl=0.3)).process_set  # citrate, no O. oeni
    on = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.2, citrate_gpl=0.3)).process_set
    # un-pitched: the citrate Process is disabled, nothing active touches citrate ⇒ VALIDATED
    assert off.tier_of("citrate", tm) is Tier.VALIDATED
    # pitched: the speculative citrate Process touches it ⇒ speculative (prime directive #1)
    assert on.tier_of("citrate", tm) is Tier.SPECULATIVE


# -- 6. tiers ---------------------------------------------------------------------------


def test_processes_are_speculative():
    assert MalolacticCitrateMetabolism.tier is Tier.SPECULATIVE
    assert OenococcusDiacetylReduction.tier is Tier.SPECULATIVE
