"""D-56/D-57: independent-data comparison against Varela, Pizarro & Agosin (2004).

Varela, C., Pizarro, F., & Agosin, E. (2004). "Biomass Content Governs Fermentation
Rate in Nitrogen-Deficient Wine Musts." Appl. Environ. Microbiol. 70(6):3392-3400.
doi:10.1128/AEM.70.6.3392-3400.2004.

This is the project's first check against a dataset genuinely independent of the
papers this model's own parameters were fit to (Coleman, Fish & Block 2007 and
friends) -- Pontificia Universidad Catolica de Chile, no author/lab overlap. Same
strain though (S. cerevisiae EC1118/Prise de Mousse), so this is a clean two-lab
comparison, not a strain confound.

**The model ran 2-4x too fast against Varela's measured endpoints (D-56); D-57
found and fixed a real bug that explains most of the extra N=50 gap, narrowing but
NOT closing the residual.** D-56 first attributed the uniform in-range gap to
q_sugar_max applying to "total" biomass with no active/inactive split -- that
diagnosis was WRONG (stale note: it described a pre-D-13 model; D-13 already gave
X/X_dead the exact split Coleman's own X_A/death-into-a-pool eqs 1-2 use, verified
line-for-line by test_coleman_reconstruction.py). The real bug D-57 found: Coleman's
own death-rate constant k_prime_d is the one QUADRATIC-in-temperature parameter in
his fit, and it shipped with NO temperature modifier at all ("M1 is isothermal at
20 C" -- true when written, stale once M2 added non-isothermal scenarios). So every
non-20 C run drove growth/uptake with Arrhenius scaling while leaving death frozen
at the 20 C rate -- inert on short/high-N runs, but compounding badly on the long,
nitrogen-limited runs this file exercises. D-57 wires in
``ColemanQuadraticDeathTemperature`` (Coleman's own regression, not an approximation)
and narrows the model's N50/N300 duration-ratio SHORTFALL against Varela's real
ratio (4.12x) from ~1.94x-too-small (2.12x modeled, pre-D-57) to ~1.17x-too-small
(3.53x modeled, post-D-57) -- real progress, NOT closure. These tests still do not
assert agreement with Varela's numbers -- they characterize the model's CURRENT gap
as a regression guard, so a future change that silently narrows OR widens it gets
caught either way. Do not widen these bands to make CI green without updating D-57
(or adding a new decision) to say why the underlying behaviour changed.

**Biomass comparison variable, corrected alongside D-57 (independent of the
k_prime_d fix itself):** Varela measures TOTAL dry cell weight (gravimetric
filtration, dead + viable cells combined -- confirmed against the paper's own
Materials and Methods), so the comparison variable is ``X + X_dead``, not viable
``X`` alone. Comparing viable-only was a standing D-56-era oversight, invisible
pre-D-57 because unscaled (too-weak) death kept ``X_dead`` small enough that
viable X approximated total. The corrected total-biomass comparison reproduces
D-56 finding 3 (the separate, already-documented Y_X/N cross-study yield gap)
almost exactly (~42% low at N=300, ~7% low at N=50) -- so these biomass
assertions now guard that growth-yield finding, cleanly separated from the
duration assertions above, which guard death/uptake timing.

**A related, NOT-fixed caveat worth flagging, not burying:** at N=50 the model's
own viable/dead SPLIT implies ~94-98% of biomass is "dead" by the time dryness
arrives, while Varela separately reports >97% viability throughout (LIVE/DEAD
membrane-integrity staining). This is not read as evidence k_prime_d's magnitude
is wrong: Coleman's own reference model shows the identical near-total X_A crash
at N=50 (see D-57), so the engine is faithfully reproducing Coleman here, not
diverging from him -- and ``X_dead`` is documented (inactivation.py) as a loss of
CATALYTIC (fermentative) capacity, the classic yeast *vitality* concept, which is
a different quantity from LIVE/DEAD *viability* (membrane integrity); the two are
not expected to agree. k_prime_d's magnitude is out of scope here (changing it
would break the Coleman line-for-line reconstruction) -- flagged for whoever next
touches death-rate calibration, not fixed by D-57.
"""

import numpy as np
import pytest
from scipy.optimize import brentq

from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import brix_to_sugar_gpl

pytestmark = pytest.mark.benchmark

#: Residual sugar [g/L] defining "dry" -- same threshold as test_milestone1.py.
DRYNESS_GPL = 2.0

#: Varela's synthetic MS300 must: 120 g/L glucose + 120 g/L fructose, both fully
#: fermentable. Scenario.initial only accepts "brix" for wine (brix_to_sugar_gpl *
#: must_fermentable_fraction), and must_fermentable_fraction=0.93 is a REAL-GRAPE-MUST
#: correction (dissolved solids that aren't fermentable hexose) that does not apply to
#: a pure-sugar synthetic must -- so it under-loads S0 by ~7% here. Documented, not
#: fixed: nowhere near enough to explain a 2-4x timing gap (D-56).
_TARGET_SUGAR_GPL = 240.0
_FERMENTABLE_FRACTION = 0.93


def _brix_for_sugar(target_gpl: float) -> float:
    def f(b: float) -> float:
        return brix_to_sugar_gpl(b) * _FERMENTABLE_FRACTION - target_gpl

    return float(brentq(f, 1.0, 40.0))


_BRIX = _brix_for_sugar(_TARGET_SUGAR_GPL)

#: Varela pitched 10^6 cells/mL (a research-lab dose, not a commercial pitch rate).
#: Converted via the standard ~18 pg/cell S. cerevisiae dry-weight figure (consistent
#: with the textbook ~1 OD600 ~ 3e7 cells/mL ~ 0.5 g/L DCW rule of thumb):
#'   10^6 cells/mL = 10^9 cells/L * ~18 pg/cell ~= 0.018 g/L
#: An order-of-magnitude conversion, not exact -- but the plausible range (~0.01-0.03
#: g/L) cannot explain a 2-4x duration gap on its own (D-56).
_PITCH_GPL_RESEARCH = 0.018

#: Varela's measured endpoints (Table 1; 3 independent experiments each; mean +/- SD).
_VARELA_HOURS_TO_DRYNESS = {300.0: 170.0, 50.0: 700.0}
_VARELA_FINAL_BIOMASS_GPL = {300.0: 5.8, 50.0: 1.5}


def _run_varela_condition(yan_mgl: float, duration_days: float) -> tuple[float, float]:
    """Run the wine scenario at Varela's conditions; return (hours_to_dryness,
    TOTAL biomass -- viable ``X`` + ethanol-inactivated ``X_dead`` -- at the moment
    dryness is first reached).

    Varela's Table 1 is dry cell weight by gravimetric filtration ("dried...to a
    constant weight at 85 C"), which does not distinguish live from dead cells --
    confirmed against the paper directly, not assumed. ``X`` alone (viable-only) is
    the WRONG comparison variable for it (a D-56-era test oversight, masked pre-D-57
    because weak, unscaled death made ``X_dead`` negligible so viable-only looked
    like a reasonable proxy for total). Total is also the more robust quantity to
    assert: EthanolInactivation only transfers mass between ``X``/``X_dead`` (D-13),
    so ``X + X_dead`` stops changing the moment nitrogen-limited growth ends
    (confirmed flat to 5 significant figures from ~40 h onward at both N levels
    tested here) -- unlike viable ``X`` alone, it carries no dependence on exactly
    when dryness is crossed."""
    scenario = Scenario(
        name=f"varela2004-{yan_mgl:g}mgN",
        medium="wine",
        initial={"brix": _BRIX, "yan_mgl": yan_mgl, "pitch_gpl": _PITCH_GPL_RESEARCH},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=28.0)],
        duration_days=duration_days,
    )
    compiled = compile_scenario(scenario, strict=True)
    duration_h = compiled.t_span_h[1]
    t_eval = np.linspace(0.0, duration_h, int(duration_h) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message
    sugar = np.asarray(traj.series("S"))
    total_sugar = sugar if sugar.ndim == 1 else sugar.sum(axis=0)
    reached = np.where(total_sugar <= DRYNESS_GPL)[0]
    total_biomass = np.asarray(traj.series("X")) + np.asarray(traj.series("X_dead"))
    if not reached.size:
        raise AssertionError(f"N={yan_mgl} mg/L never reached dryness within {duration_days} d")
    idx = reached[0]
    return float(traj.t[idx]), float(total_biomass[idx])


def test_varela2004_normal_n_gap_characterized():
    # 300 mg N/L (well-fed) is inside Coleman's fitted 70-350 mg N/L range and inside
    # its 11-35 C temperature range, so this is the CLEANEST comparison: same strain,
    # in-range conditions. The model still runs ~1.9x too fast. D-57 barely moves this
    # band (89 h vs pre-D-57's 83 h): at N=300 dryness arrives quickly enough that
    # cumulative ethanol-driven death is a minor contributor either way, so correcting
    # its temperature scaling has little to bite on here -- exactly as expected (the
    # fix mattered for the LONG N=50 run below, not this short one). This residual gap
    # is a genuine cross-study Coleman-vs-Varela difference (see D-57): Coleman's own
    # reference model (test_coleman_reconstruction.py's reconstruction, run at 28 C)
    # gives ~84.5 h for these exact inputs, matching the engine -- so the engine
    # faithfully reproduces Coleman even off his 20 C reference point, and the model
    # is not the source of the remaining gap to Varela's real 170 h. This is a
    # regression guard on that characterized gap, not a pass/fail against Varela.
    # Biomass here is TOTAL (X + X_dead, Varela's own dry-cell-weight measure; see
    # module docstring) and guards a SEPARATE finding from the duration assertion
    # below: the Y_X/N growth-yield cross-study gap (D-56 finding 3), not death/uptake
    # timing. It is also timing-insensitive (flat from ~40 h onward, confirmed), unlike
    # a viable-biomass reading would be.
    hours, biomass = _run_varela_condition(300.0, duration_days=15.0)
    assert 80.0 <= hours <= 100.0, (
        f"N=300: hours_to_dryness={hours:.1f} outside characterized [80, 100] h"
    )
    assert 3.0 <= biomass <= 3.7, (
        f"N=300: total_biomass_at_dryness={biomass:.2f} outside characterized [3.0, 3.7] g/L"
    )

    # The gap to Varela's real 170 +/- 12 h should stay in the characterized ~1.6-2.2x
    # band -- if this ever falls to ~1x, that is not a free win, it means something
    # shifted that should be recorded as a new decision; if it grows past ~2.2x, a
    # regression narrowed the model's fermentation rate.
    gap_ratio = _VARELA_HOURS_TO_DRYNESS[300.0] / hours
    assert 1.6 <= gap_ratio <= 2.2, (
        f"N=300: gap ratio to Varela {gap_ratio:.2f} outside characterized [1.6, 2.2]x"
    )


def test_varela2004_deficient_n_gap_characterized():
    # 50 mg N/L is BELOW Coleman's fitted 70-350 mg N/L floor -- genuine extrapolation,
    # not just a harder in-range case. Pre-D-57 the model ran ~4x too fast here (vs ~2x
    # at N=300); D-57's ColemanQuadraticDeathTemperature fix (see module docstring)
    # closes MOST of that extra gap -- this long, nitrogen-limited run is exactly where
    # a death rate silently frozen at the mild 20 C value (instead of the correct,
    # much steeper 28 C value) had the most hours to compound. Post-fix the model still
    # runs ~2.2x too fast here (down from ~4x): a smaller, genuine residual remains
    # (see test_varela2004_deficiency_widens_the_gap), left diagnosed not patched --
    # see D-57 before touching growth.py/uptake.py to "fix" this further, and read
    # D-57 before assuming a novel N-gated transporter mechanism (D-56's original
    # hypothesis) is still the right next step; most of what motivated it turned out
    # to be this simpler, sourced bug.
    #
    # Biomass here is TOTAL (X + X_dead, see module docstring) -- ~1.40 g/L, within 7%
    # of Varela's 1.5 (D-56 finding 3's already-documented near-match at this N level).
    # It is flat from ~40 h onward regardless of the death-rate fix (mass-neutral
    # X<->X_dead transfer), so this band is not knife-edge like a viable-only reading
    # would be. This is NOT a viability claim -- see module docstring's vitality/
    # viability caveat before reading anything into the model's internal X/X_dead split.
    hours, biomass = _run_varela_condition(50.0, duration_days=40.0)
    assert 280.0 <= hours <= 350.0, (
        f"N=50: hours_to_dryness={hours:.1f} outside characterized [280, 350] h"
    )
    assert 1.25 <= biomass <= 1.55, (
        f"N=50: total_biomass_at_dryness={biomass:.3f} outside characterized [1.25, 1.55] g/L"
    )

    gap_ratio = _VARELA_HOURS_TO_DRYNESS[50.0] / hours
    assert 1.9 <= gap_ratio <= 2.6, (
        f"N=50: gap ratio to Varela {gap_ratio:.2f} outside characterized [1.9, 2.6]x"
    )


def test_varela2004_deficiency_widens_the_gap():
    # The central, structural D-56 finding survives D-57's fix, just narrower: the
    # model's own N50/N300 slowdown ratio is still smaller than Varela's real one
    # (3.53x modeled vs Varela's 4.12x, post-D-57 -- was 2.12x vs 4.12x pre-D-57). If
    # this test starts failing because the model's ratio catches up to (or exceeds)
    # Varela's ~4.1x, that is a real change in the nitrogen-limitation kinetics worth
    # investigating -- not just noise.
    hours_300, _ = _run_varela_condition(300.0, duration_days=15.0)
    hours_50, _ = _run_varela_condition(50.0, duration_days=40.0)
    model_ratio = hours_50 / hours_300
    varela_ratio = _VARELA_HOURS_TO_DRYNESS[50.0] / _VARELA_HOURS_TO_DRYNESS[300.0]
    assert model_ratio < varela_ratio, (
        f"model N50/N300 duration ratio {model_ratio:.2f} should stay BELOW Varela's "
        f"real {varela_ratio:.2f} (D-56/D-57: the model still under-predicts how much "
        "severe nitrogen deficiency slows fermentation, relative to an in-range baseline, "
        "though D-57 narrowed the shortfall considerably)"
    )
