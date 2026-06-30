"""Milestone 1 acceptance benchmarks (handoff section 2.2).

These are the §2.2 wine and beer acceptance criteria, encoded as executable tests
so the kinetics are test-driven against them. The wine, beer, and CO2 criteria are
live; the one still skipped is the Tier-2 byproduct-quality directional check,
which needs Processes outside Milestone 1's scope. Do not delete or weaken these to
make CI green — implement the model until they pass (handoff section 2.2: "encode
the benchmark, then iterate").
"""

import numpy as np
import pytest

from fermentation.core.chemistry import co2_yield, sugar_species
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import abv_from_ethanol, apparent_gravity, sg_to_plato
from fermentation.validation import BENCHMARKS, assert_conserved, total_carbon

pytestmark = pytest.mark.benchmark

KINETICS_PENDING = "Milestone 1: primary-fermentation kinetics not implemented yet"

#: Residual sugar [g/L] defining "dry". Below Coleman's stuck threshold (0.4 %
#: w/v ~ 4 g/L); 2 g/L is a solidly dry wine.
DRYNESS_GPL = 2.0

#: The §2.2 wine benchmark must (Coleman low-N anchor; see test_wine_*). Shared by the
#: dryness-window test and the Tier-2 temperature-direction benchmark.
_WINE_BENCH = {"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25}


def _days_to_dryness(scenario: Scenario) -> float:
    """Integrate ``scenario`` and return the day total sugar first falls to
    :data:`DRYNESS_GPL`. Returns ``inf`` if it never gets there (stuck)."""
    compiled = compile_scenario(scenario, strict=True)
    duration_h = compiled.t_span_h[1]
    t_eval = np.linspace(0.0, duration_h, int(duration_h) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message
    sugar = np.asarray(traj.series("S"))
    total = sugar if sugar.ndim == 1 else sugar.sum(axis=0)
    reached = np.where(total <= DRYNESS_GPL)[0]
    return float(traj.t[reached[0]] / 24.0) if reached.size else float("inf")


def test_wine_24brix_ferments_to_dryness_in_window():
    # A nitrogen-limited 24 Brix must at 20 C. Conditions are anchored to the
    # keystone source (Coleman low-N = 80 mg N/L; ~0.25 g/L = 25 g/hL pitch), NOT
    # tuned to the window — the validated core reproduces Coleman line-for-line
    # (decision D-14), and the window was re-anchored to that source.
    spec = BENCHMARKS["wine_dryness"]
    scenario = Scenario(
        name="wine-benchmark",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )
    days_to_dryness = _days_to_dryness(scenario)
    assert spec.passes(days_to_dryness), (
        f"days_to_dryness={days_to_dryness:.2f} outside [{spec.low}, {spec.high}] d"
    )


def test_wine_abv_and_glycerol_are_realistic():
    # Decision D-16: the realised-yield byproduct sink + must-fermentable-fraction
    # correction make a 24 Brix wine finish at a realistic ABV with realistic
    # glycerol — *without* tuning to a target. Each input (glycerol/byproduct yields,
    # fermentable fraction) is independently sourced; ABV and glycerol fall out. This
    # is a realism regression guard, not a §2.2 acceptance spec (no benchmark gates on
    # absolute ABV); the bands are the literature ranges with margin, not fitted.
    scenario = Scenario(
        name="wine-abv",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )
    compiled = compile_scenario(scenario, strict=True)
    duration_h = compiled.t_span_h[1]
    t_eval = np.linspace(0.0, duration_h, int(duration_h) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message

    sugar = np.asarray(traj.series("S"))
    s0 = float(sugar[0])
    s_final = float(sugar[-1])
    ethanol_final = float(traj.series("E")[-1])
    glycerol_final = float(traj.series("Gly")[-1])
    byproduct_final = float(traj.series("Byp")[-1])
    abv = abv_from_ethanol(ethanol_final)
    realised_yield = ethanol_final / (s0 - s_final)

    # Realistic 24 Brix potential alcohol (~14-15 %), not the 16.9 % the theoretical
    # split gave before D-16.
    assert 13.5 <= abv <= 15.5, f"wine ABV {abv:.2f}% outside realistic 13.5-15.5%"
    # Dry-wine glycerol is 4-10 g/L (Ribereau-Gayon); band carries margin.
    assert 5.0 <= glycerol_final <= 11.0, f"glycerol {glycerol_final:.2f} g/L outside 5-11"
    assert byproduct_final > 0.0  # the minor-byproduct lump accumulates too
    # Emergent realised yield sits in the literature 0.46-0.48 band (here ~0.48),
    # cross-checking Y_ethanol_sugar without being set to it.
    assert 0.46 <= realised_yield <= 0.50, f"realised Y_E {realised_yield:.4f} outside 0.46-0.50"

    # Carbon still closes to machine precision with byproducts tracked (D-16).
    f_c = compiled.parameters.value("biomass_C_fraction")
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )


#: Apparent (hydrometer) final gravity defining beer attenuation "done" (handoff
#: §2.2: ~1.010). It is an *apparent*, ethanol-depressed reading, not real extract.
TARGET_FG_SG = 1.010

#: A ~1.048 OG all-malt ale wort. The fermentable sugar spectrum (glucose/maltose/
#: maltotriose ≈ 15/62/23 % of fermentables) and the real degree of fermentation
#: (S0 ≈ 88 g/L of the ~125 g/L total extract at 1.048, RDF ~70 %) are sourced, NOT
#: tuned to the window: S0 ≈ 88 g/L is the initial fermentable sugar measured in
#: our beer source (Zamudio Lara et al. 2022), and the split is typical all-malt
#: (e.g. Briggs et al., "Brewing: Science and Practice"). The unfermentable extract
#: plus the ethanol-depressed apparent gravity make the wort finish near 1.007 —
#: comfortably below the 1.010 target — so the *endpoint* falls out of the wort and
#: the 1.010 crossing lands in the kinetic phase, not at a fragile asymptote (D-15).
_BEER_OG_SG = 1.048
_BEER_WORT: dict[str, float] = {
    "glucose_gpl": 13.2,
    "maltose_gpl": 54.6,
    "maltotriose_gpl": 20.2,
    "yan_mgl": 200.0,  # typical ale wort YAN (~150-250 mg/L)
    "pitch_gpl": 0.6,  # typical ale pitch, dry-cell-weight equivalent
}
_BEER_FERMENTABLE_S0 = (
    _BEER_WORT["glucose_gpl"] + _BEER_WORT["maltose_gpl"] + _BEER_WORT["maltotriose_gpl"]
)


def _beer_scenario(duration_days: float = 14.0) -> Scenario:
    return Scenario(
        name="beer-benchmark",
        medium="beer",
        initial=dict(_BEER_WORT),
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=duration_days,
    )


def _simulate_beer(scenario: Scenario):
    compiled = compile_scenario(scenario, strict=True)
    duration_h = compiled.t_span_h[1]
    t_eval = np.linspace(0.0, duration_h, int(duration_h) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message
    return compiled, traj


def _apparent_gravity_series(traj, og_sg: float, fermentable_s0_gpl: float):
    """Apparent (hydrometer) SG over a beer trajectory.

    Real extract = unfermentable extract + residual fermentable sugar; the
    unfermentable share is implicit in (OG extract − S0), so no extra state or
    parameter is needed. That real extract, depressed by the ethanol present
    (Balling/Tabarie), is the apparent gravity brewers quote as FG — which is why
    "1.010" is reachable by a realistic ~66 %-fermentable wort (real FG ~1.016)."""
    sugar = np.asarray(traj.series("S"))
    residual = sugar.sum(axis=0) if sugar.ndim == 2 else sugar
    oe_plato = sg_to_plato(og_sg)
    ferm_plato = fermentable_s0_gpl / (10.0 * og_sg)  # fermentable extract at OG
    unferm_plato = oe_plato - ferm_plato
    frac_remaining = np.clip(residual / fermentable_s0_gpl, 0.0, 1.0)
    real_extract_plato = unferm_plato + frac_remaining * ferm_plato
    return np.array([apparent_gravity(float(re), oe_plato) for re in real_extract_plato])


def test_beer_1048_og_attenuates_in_5_to_7_days():
    # A ~1.048 OG ale wort at 20 C must reach ~1.010 apparent gravity in 5-7 d.
    # The wort spectrum and fermentability are sourced (see _BEER_WORT), not swept
    # to fit; q_sugar_max was re-derived to a decoupled-equivalent rate (D-15), so
    # the 5-7 d window falls out rather than being dialed.
    spec = BENCHMARKS["beer_attenuation"]
    _, traj = _simulate_beer(_beer_scenario())
    apparent_sg = _apparent_gravity_series(traj, _BEER_OG_SG, _BEER_FERMENTABLE_S0)
    reached = np.where(apparent_sg <= TARGET_FG_SG)[0]
    days = float(traj.t[reached[0]] / 24.0) if reached.size else float("inf")
    assert spec.passes(days), (
        f"days_to_target_gravity={days:.2f} outside [{spec.low}, {spec.high}] d "
        f"(apparent FG floor ~{float(apparent_sg[-1]):.4f})"
    )


def test_co2_integral_tracks_sugar_consumed():
    # CO2 is the primary measurable validation channel: the evolved-CO2 integral
    # tracks the CO2 predicted from sugar consumed by Gay-Lussac stoichiometry,
    # summed over all three sugar slots (so the maltose 2x / maltotriose 3x hexose
    # factors are exercised; a flat single-hexose yield would mis-weight them).
    # The ratio is slightly BELOW 1: ~2-3% of sugar carbon is routed into biomass
    # by growth (no anabolic CO2 in M1), plus (since D-19) a trace routed into the
    # ester/fusel pools from sugar — together ~2.5% here (ratio 0.977 -> 0.975) — so a
    # touch less CO2 is evolved than total sugar consumed implies. The [0.95, 1.05]
    # window accommodates that diversion —
    # this is the measurable-channel check, not the machine-precision carbon audit
    # (that lives in the conservation tests). Run on beer to cover the 3-slot sum.
    spec = BENCHMARKS["co2_peak_then_tail"]
    compiled, traj = _simulate_beer(_beer_scenario())
    co2 = np.asarray(traj.series("CO2"))
    sugar = np.asarray(traj.series("S"))  # (3, n): glucose/maltose/maltotriose
    species = sugar_species(compiled.schema)
    consumed = sugar[:, [0]] - sugar  # per-species g/L consumed from each S0
    expected_co2 = sum(co2_yield(sp) * consumed[i] for i, sp in enumerate(species))
    ratio = float(co2[-1] / expected_co2[-1])
    assert spec.passes(ratio), (
        f"co2/sugar-consumed ratio={ratio:.4f} outside [{spec.low}, {spec.high}]"
    )
    # The spec's qualitative claim with real kinetic teeth: the CO2 *evolution
    # rate* rises to a peak then tails off as sugar depletes, so the cumulative
    # curve is sigmoidal — its gradient peaks in the interior, not at t=0.
    rate = np.gradient(co2, traj.t)
    peak = int(np.argmax(rate))
    assert 0 < peak < len(rate) - 1, "CO2 evolution rate should peak in the interior"
    assert rate[-1] < rate[peak], "CO2 evolution rate should tail off below its peak"


def _aroma_run(medium: str, initial: dict[str, float], celsius: float, duration_days: float):
    """Run a medium isothermally; return (days_to_dryness, liquid-pool dict).

    Reads the **liquid** ``esters``/``fusels`` pools only — ``esters_gas`` is the
    bookkeeping headspace pool (volatilized away), not aroma in the glass (D-20)."""
    sc = Scenario(
        name=f"{medium}-{celsius}C",
        medium=medium,
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        duration_days=duration_days,
    )
    days = _days_to_dryness(sc)
    compiled = compile_scenario(sc, strict=True)
    duration_h = compiled.t_span_h[1]
    t_eval = np.linspace(0.0, duration_h, int(duration_h) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message
    pools = {
        "esters": float(traj.series("esters")[-1]),
        "fusels": float(traj.series("fusels")[-1]),
    }
    return days, pools


def test_lower_temperature_is_slower_but_cleaner():
    # The §2.2 Tier-2 directional benchmark, honest per medium (decisions D-19 + D-20).
    # "Cleaner when colder" is real but NOT a single combined ester+fusel total: building
    # the volatilization sink (D-20) revealed that wine LIQUID esters *invert* — they
    # rise as T falls — because warm ferments strip esters into the headspace faster than
    # they synthesise them (Rollero 2014). So the honest, sourced directions are:
    #
    #   Both media:  colder ⇒ slower to dryness, AND fewer FUSELS (the harsh higher
    #                alcohols — the real "cleaner"; fusel synthesis rises with T).
    #   Beer:        colder ⇒ fewer liquid esters too (synthesis-dominated; warm ales are
    #                estery — de Andrés-Toro 1998).
    #   Wine:        colder ⇒ MORE liquid esters (stripping-dominated inversion — the
    #                warm ferment's esters end up volatilized, not in the wine).
    #
    # Reads the liquid pools only; the volatilized esters_gas headspace pool is not aroma
    # in the glass. Asserting a combined total here would hide the wine inversion that the
    # sink was built to surface — so we assert each pool's sourced direction explicitly.
    wine_cold_days, wine_cold = _aroma_run("wine", _WINE_BENCH, 14.0, 90.0)
    wine_warm_days, wine_warm = _aroma_run("wine", _WINE_BENCH, 25.0, 30.0)
    beer_cold_days, beer_cold = _aroma_run("beer", dict(_BEER_WORT), 14.0, 40.0)
    beer_warm_days, beer_warm = _aroma_run("beer", dict(_BEER_WORT), 25.0, 18.0)

    # Slower when colder (both media must actually reach dryness for the comparison).
    for label, cold_d, warm_d in (
        ("wine", wine_cold_days, wine_warm_days),
        ("beer", beer_cold_days, beer_warm_days),
    ):
        assert np.isfinite(cold_d) and np.isfinite(warm_d), f"{label}: both must reach dryness"
        assert cold_d > warm_d, f"{label}: colder should be slower ({cold_d:.1f} vs {warm_d:.1f} d)"

    # Cleaner when colder = fewer FUSELS, both media (the sourced "cleaner" direction).
    assert 0.0 < wine_cold["fusels"] < wine_warm["fusels"], "wine: colder ⇒ fewer fusels"
    assert 0.0 < beer_cold["fusels"] < beer_warm["fusels"], "beer: colder ⇒ fewer fusels"

    # Beer liquid esters: synthesis-dominated, so fewer when colder too.
    assert 0.0 < beer_cold["esters"] < beer_warm["esters"], "beer: colder ⇒ fewer liquid esters"

    # Wine liquid esters: the D-20 inversion — MORE when colder (warm strips them off).
    assert 0.0 < wine_warm["esters"] < wine_cold["esters"], (
        "wine: liquid esters should INVERT — more when colder (volatilization-dominated, "
        f"Rollero 2014): cold {wine_cold['esters']:.4f} vs warm {wine_warm['esters']:.4f} g/L"
    )
