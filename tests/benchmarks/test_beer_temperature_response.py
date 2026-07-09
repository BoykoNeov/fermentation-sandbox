"""D-63: beer temperature-response check — a cross-regime Arrhenius stress test.

**What this is, and (loudly) what it is not.** The project has a standing open
gap: no genuinely independent, in-regime (isothermal ale, ~1.048 OG, numeric
time-series) beer dataset is publicly accessible — the two richest candidates
are the model's own beer fit sources (circular; see D-59). The one independent
option found is an *off-regime* lager dataset. This file is the honest,
scope-limited thing that off-regime data supports — NOT a validation of absolute
lager kinetics.

**The data situation (D-63 investigation).** The only freely reconstructable
lager curve is Reid, Josey, MacIntosh, Maskell & Speers (2021), "Predicting
Fermentation Rates in Ale, Lager and Whisky," *Fermentation* 7(1):13
(doi:10.3390/fermentation7010013), Table 2 — an Australian lager at Original
Extract 14.1 °P, **single starting temperature 10 °C**, fit with a 3-parameter
ADF logistic (slope B = 0.06372 h⁻¹, midpoint M = 51.22 h ≈ 2.1 d). The
underlying fermentations are from Speers, Rogers & Smith (2003), "Non-linear
modelling of industrial brewing fermentations," *J. Inst. Brew.* 109(3):229–235
(doi:10.1002/j.2050-0416.2003.tb00163.x), which is where the *multi-temperature*
signal lives (starting temperature raises fermentation rate, p<0.01). Speers
2003 is Wiley-paywalled, and its temperature effect is a regression across many
*industrial* batches (brand/wort/pitch co-vary with temperature), so even if
obtained it may not be a clean controlled series.

**Why single-temperature lager data cannot be a validation band.** Comparing
"the engine's ale-yeast Arrhenius kinetics extrapolated to 10 °C" against "real
lager yeast (S. pastorianus) in a 14.1 °P industrial wort at 10 °C" conflates
the Arrhenius temperature law with the organism + wort + pitch-rate difference.
Concretely: the engine's low-pitch (0.6 g/L, homebrew-like) 10 °C run reaches
its attenuation midpoint at ~6.2 d, ~2.9× SLOWER than Speers' ~2.1 d industrial
midpoint — a gap dominated by pitch rate and organism, not the temperature
model. Guarding that gap as a regression band would guard a confound. So this
file does NOT assert the engine reproduces the 51 h midpoint; it deliberately
misses it.

**The confound-cancelling ratio test that WOULD have real signal is deferred.**
A rate *ratio* across two temperatures cancels the absolute-kinetics difference
between lager and ale yeast, isolating the temperature axis. Building it needs
Speers 2003's *controlled* temperature series (fitted rate/midpoint at ≥2
temperatures on one wort+yeast). Pending that source, the reusable helper
:func:`_apparent_activation_energy` below is the drop-in point: feed it two
(temperature, midpoint) pairs from Speers and compare the empirical E_a to the
engine's. See D-63 in docs/DECISIONS.md.

**What this file DOES assert (three claims, from the engine's own two-temperature
runs):**

  1. *Wiring / regression guard* (:func:`test_arrhenius_modifier_composes_over_
     full_beer_ferment`) — the whole-ferment apparent E_a recovered from the
     20 °C and 10 °C beer midpoints equals the input E_a's. This guards, on the
     BEER medium and over a full ferment composing BOTH growth and uptake, that
     the Arrhenius modifiers stay wired into fermentation timing — the D-57
     frozen-modifier bug class (an Arrhenius term left frozen at T_ref while
     others scale). Existing Arrhenius tests are directional only, uptake-only,
     and on wine; this is the quantitative, whole-ferment, beer extension.

  2. *Reality check* (:func:`test_beer_temperature_sensitivity_within_literature_
     range`) — the SAME recovered E_a sits inside the empirically observed range
     of apparent activation energies for S. cerevisiae alcoholic fermentation
     (~35–100 kJ/mol; free-cell values commonly ~62–97 kJ/mol, e.g. immobilised-
     cell kinetic studies, with lower ~35 kJ/mol also reported). This is the only
     reality-touching claim, and it is the one with teeth: it excludes the
     ~265 kJ/mol lumped-fit artifact the beer file explicitly rejected (de
     Andrés-Toro 1998; see beer_generic.yaml header) while staying humble about
     the lager/ale organism gap.

  3. *Cross-regime order-of-magnitude anchor* (:func:`test_cold_lager_completes_
     in_cross_regime_order_of_magnitude`) — the engine at 10 °C reaches near-
     completion within a loose "cold lager takes ~1–2 weeks" window. This is
     CONFOUNDED (organism + pitch) and deliberately loose: it only catches a
     temperature model that is order-of-magnitude wrong (hours, or months). The
     low-pitch assumption is why the engine sits at the slow end and does not
     match Speers' fast industrial timing.

**Firewall (prime directive 2): clean.** The engine's E_a's derive from the
Coleman 2007 wine fit; the reference data is Speers/Reid lager. Disjoint sources,
so comparing them is not self-confirming.
"""

import numpy as np
import pytest

from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

pytestmark = pytest.mark.benchmark

#: The §2.2 ~1.048 OG all-malt wort (glucose/maltose/maltotriose spectrum sourced
#: from Zamudio Lara et al. 2022; see tests/benchmarks/test_milestone1.py). Reused
#: unchanged so this file exercises the same beer core, only at a colder T.
_BEER_WORT: dict[str, float] = {
    "glucose_gpl": 13.2,
    "maltose_gpl": 54.6,
    "maltotriose_gpl": 20.2,
    "yan_mgl": 200.0,
    "pitch_gpl": 0.6,  # homebrew-like ale pitch; << industrial lager pitch (see docstring)
}
_S0 = _BEER_WORT["glucose_gpl"] + _BEER_WORT["maltose_gpl"] + _BEER_WORT["maltotriose_gpl"]

#: The engine's input Arrhenius activation energies (beer_generic.yaml). The
#: whole-ferment apparent E_a should recover these; fermentation pace is set by
#: uptake ("q_sugar_max sets beer's pace"), so the composed value tracks uptake.
_E_A_UPTAKE = 55100.0  # J/mol
_E_A_GROWTH = 55900.0  # J/mol

#: Independent literature range for S. cerevisiae alcoholic-fermentation apparent
#: activation energy (see docstring for sources). The engine's ~55 kJ/mol sits
#: inside; the rejected de Andrés-Toro ~265 kJ/mol artifact sits far outside.
_E_A_LITERATURE_LOW = 35000.0  # J/mol
_E_A_LITERATURE_HIGH = 100000.0  # J/mol

_R = 8.314  # J/(mol·K)
_KELVIN = 273.15


def _beer_scenario(celsius: float, duration_days: float) -> Scenario:
    return Scenario(
        name=f"beer-lager-{celsius:g}C",
        medium="beer",
        initial=dict(_BEER_WORT),
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        duration_days=duration_days,
    )


def _attenuation_series(scenario: Scenario) -> tuple[np.ndarray, np.ndarray]:
    """Integrate ``scenario``; return (hours, attenuated fraction) where the
    fraction is 1 − Σsugar/S0 over the three-slot beer sugar vector."""
    compiled = compile_scenario(scenario, strict=True)
    duration_h = compiled.t_span_h[1]
    t_eval = np.linspace(0.0, duration_h, int(duration_h) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message
    sugar = np.asarray(traj.series("S"))
    total = sugar.sum(axis=0) if sugar.ndim == 2 else sugar
    return np.asarray(traj.t), 1.0 - total / _S0


def _time_to_attenuation_h(scenario: Scenario, target_fraction: float) -> float:
    """Hours until ``target_fraction`` of the fermentable sugar is consumed
    (``inf`` if never)."""
    t_h, frac = _attenuation_series(scenario)
    reached = np.where(frac >= target_fraction)[0]
    return float(t_h[reached[0]]) if reached.size else float("inf")


def _apparent_activation_energy(
    midpoint_cold_h: float, t_cold_c: float, midpoint_warm_h: float, t_warm_c: float
) -> float:
    """Apparent Arrhenius activation energy [J/mol] from two (temperature,
    midpoint-time) pairs. Fermentation rate ∝ 1/midpoint, so

        E_a = R · ln(midpoint_cold / midpoint_warm) / (1/T_cold − 1/T_warm)

    with T in Kelvin. Reusable: the deferred Speers-based ratio test feeds this
    the real lager midpoints at two temperatures and compares to the engine's."""
    t_cold_k = t_cold_c + _KELVIN
    t_warm_k = t_warm_c + _KELVIN
    return float(_R * np.log(midpoint_cold_h / midpoint_warm_h) / (1.0 / t_cold_k - 1.0 / t_warm_k))


def _engine_apparent_e_a() -> float:
    """Recover the whole-ferment apparent E_a from the engine's own 20 °C and
    10 °C beer runs (attenuation midpoint = 50 % of fermentables consumed)."""
    warm_h = _time_to_attenuation_h(_beer_scenario(20.0, 21.0), 0.50)
    cold_h = _time_to_attenuation_h(_beer_scenario(10.0, 40.0), 0.50)
    return _apparent_activation_energy(cold_h, 10.0, warm_h, 20.0)


def test_arrhenius_modifier_composes_over_full_beer_ferment():
    # CLAIM 1 (wiring / regression guard). The apparent E_a recovered from the
    # engine's 20 C and 10 C beer midpoints must equal the input E_a's. This
    # guards -- on beer, over a full ferment composing growth AND uptake -- that
    # the Arrhenius modifiers stay wired into fermentation timing (the D-57
    # frozen-modifier bug class). Pace is set by uptake, so the composed value
    # tracks E_a_uptake (55.1 kJ/mol), pulled slightly by growth (55.9 kJ/mol);
    # the band allows for threshold-crossing discretization.
    e_a = _engine_apparent_e_a()
    assert 50000.0 <= e_a <= 60000.0, (
        f"apparent E_a = {e_a:.0f} J/mol strays from the input E_a_uptake="
        f"{_E_A_UPTAKE:.0f} / E_a_growth={_E_A_GROWTH:.0f} J/mol -- the Arrhenius "
        f"modifier may no longer compose into both beer rates over a full ferment "
        f"(D-57 frozen-modifier class)."
    )


def test_beer_temperature_sensitivity_within_literature_range():
    # CLAIM 2 (reality check -- the honest headline). The SAME recovered E_a must
    # sit inside the empirically observed S. cerevisiae fermentation E_a range
    # (~35-100 kJ/mol; see module docstring for sources). This is the only
    # reality-touching assertion: it excludes the ~265 kJ/mol de Andres-Toro
    # lumped-fit artifact the beer file rejects, while making no claim about
    # absolute lager kinetics (organism gap). Distinct from claim 1: if the input
    # E_a were retuned to another physical value, claim 1 flags the change but
    # this stays green -- it guards physical plausibility, not a specific value.
    e_a = _engine_apparent_e_a()
    assert _E_A_LITERATURE_LOW <= e_a <= _E_A_LITERATURE_HIGH, (
        f"engine apparent E_a = {e_a:.0f} J/mol is outside the literature range "
        f"[{_E_A_LITERATURE_LOW:.0f}, {_E_A_LITERATURE_HIGH:.0f}] J/mol for yeast "
        f"alcoholic fermentation -- the assumed temperature sensitivity is "
        f"physically implausible (cf. the rejected ~265 kJ/mol artifact)."
    )


def test_cold_lager_completes_in_cross_regime_order_of_magnitude():
    # CLAIM 3 (cross-regime order-of-magnitude anchor -- CONFOUNDED, loose). At
    # 10 C the engine must reach near-completion (90% attenuation) within a "cold
    # lager takes ~1-2 weeks" window. This is deliberately loose: the low-pitch
    # (0.6 g/L) homebrew-like wort ferments slower than Speers' high-pitch
    # industrial lager, which is precisely WHY the engine (~12 d here) does not
    # match Speers' fast ~51 h industrial midpoint and we do not assert it. The
    # band only catches an order-of-magnitude-wrong temperature model (hours, or
    # months), not the organism/pitch confound.
    completion_h = _time_to_attenuation_h(_beer_scenario(10.0, 40.0), 0.90)
    completion_d = completion_h / 24.0
    assert 5.0 <= completion_d <= 25.0, (
        f"10 C beer reaches 90% attenuation at {completion_d:.1f} d, outside the "
        f"cross-regime cold-lager order-of-magnitude window [5, 25] d -- the "
        f"temperature model may be order-of-magnitude wrong."
    )
