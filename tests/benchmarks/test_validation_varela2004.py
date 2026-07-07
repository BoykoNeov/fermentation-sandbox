"""D-56: independent-data comparison against Varela, Pizarro & Agosin (2004).

Varela, C., Pizarro, F., & Agosin, E. (2004). "Biomass Content Governs Fermentation
Rate in Nitrogen-Deficient Wine Musts." Appl. Environ. Microbiol. 70(6):3392-3400.
doi:10.1128/AEM.70.6.3392-3400.2004.

This is the project's first check against a dataset genuinely independent of the
papers this model's own parameters were fit to (Coleman, Fish & Block 2007 and
friends) -- Pontificia Universidad Catolica de Chile, no author/lab overlap. Same
strain though (S. cerevisiae EC1118/Prise de Mousse), so this is a clean two-lab
comparison, not a strain confound.

**The model runs 2-4x too fast against Varela's measured endpoints, and the gap is
diagnosed, not fixed** (see D-56 for the full investigation, including a prototyped
and DISPROVED single-term fix, and why fitting further would burn the project's only
independent dataset). These tests do not assert agreement with Varela's numbers --
they characterize the model's CURRENT gap as a regression guard, so a future change
that silently narrows OR widens the gap gets caught either way, instead of the
finding decaying into a stale doc comment. Do not widen these bands to make CI green
without updating D-56 to say why the underlying behaviour changed.
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
    biomass at the moment dryness is first reached -- NOT at the end of the full
    simulated window, which would include post-fermentation ethanol-driven cell death
    unrelated to what Varela measured)."""
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
    total = sugar if sugar.ndim == 1 else sugar.sum(axis=0)
    reached = np.where(total <= DRYNESS_GPL)[0]
    x = np.asarray(traj.series("X"))
    if not reached.size:
        raise AssertionError(f"N={yan_mgl} mg/L never reached dryness within {duration_days} d")
    idx = reached[0]
    return float(traj.t[idx]), float(x[idx])


def test_varela2004_normal_n_gap_characterized():
    # 300 mg N/L (well-fed) is inside Coleman's fitted 70-350 mg N/L range and inside
    # its 11-35 C temperature range, so this is the CLEANEST comparison: same strain,
    # in-range conditions. The model still runs ~2x too fast (D-56 finding 1) -- traced
    # to q_sugar_max's OWN documented caveat (wine_generic.yaml): the rate applies to
    # TOTAL biomass with no active/inactive split, unlike Coleman's own X_A, which
    # over-catalyses the tail. This is a regression guard on that characterized gap,
    # not a pass/fail against Varela's real 170 h.
    hours, biomass = _run_varela_condition(300.0, duration_days=15.0)
    assert 70.0 <= hours <= 95.0, (
        f"N=300: hours_to_dryness={hours:.1f} outside characterized [70, 95] h"
    )
    assert 2.5 <= biomass <= 3.2, (
        f"N=300: biomass_at_dryness={biomass:.2f} outside characterized [2.5, 3.2] g/L"
    )

    # The gap to Varela's real 170 +/- 12 h should stay in the characterized ~1.7-2.5x
    # band -- if this ever falls to ~1x, that is not a free win, it means something
    # shifted (e.g. an accidental active/inactive-biomass fix) that D-56 should record;
    # if it grows past ~2.5x, a regression narrowed the model's fermentation rate.
    gap_ratio = _VARELA_HOURS_TO_DRYNESS[300.0] / hours
    assert 1.7 <= gap_ratio <= 2.5, (
        f"N=300: gap ratio to Varela {gap_ratio:.2f} outside characterized [1.7, 2.5]x"
    )


def test_varela2004_deficient_n_gap_characterized():
    # 50 mg N/L is BELOW Coleman's fitted 70-350 mg N/L floor -- genuine extrapolation,
    # not just a harder in-range case. The model runs ~4x too fast here (vs ~2x at
    # N=300), and D-56 isolated the EXTRA ~2x (beyond the uniform gap above) to the
    # nitrogen-limited growth kinetics specifically, via a biomass-hours (integral of X
    # dt) argument: since K_sugar_uptake is negligible next to S for most of the run and
    # S0 is identical at both N levels, duration in this model is governed almost
    # entirely by how fast biomass X(t) builds. A single ethanol-driven, N-gated decline
    # term was prototyped and could NOT close both this gap and the N=300 gap
    # simultaneously (structural, not a sweep-resolution issue) -- so this is left
    # diagnosed, not patched. See D-56 before touching growth.py/uptake.py to "fix" this.
    hours, biomass = _run_varela_condition(50.0, duration_days=40.0)
    assert 155.0 <= hours <= 200.0, (
        f"N=50: hours_to_dryness={hours:.1f} outside characterized [155, 200] h"
    )
    assert 0.75 <= biomass <= 1.05, (
        f"N=50: biomass_at_dryness={biomass:.2f} outside characterized [0.75, 1.05] g/L"
    )

    gap_ratio = _VARELA_HOURS_TO_DRYNESS[50.0] / hours
    assert 3.2 <= gap_ratio <= 4.8, (
        f"N=50: gap ratio to Varela {gap_ratio:.2f} outside characterized [3.2, 4.8]x"
    )


def test_varela2004_deficiency_widens_the_gap():
    # The central, structural D-56 finding, independent of the exact band edges above:
    # the model's own N50/N300 slowdown ratio is smaller than Varela's real one. If
    # this test starts failing because the model's ratio catches up to (or exceeds)
    # Varela's ~4.1x, that is a real change in the nitrogen-limitation kinetics worth
    # investigating -- not just noise.
    hours_300, _ = _run_varela_condition(300.0, duration_days=15.0)
    hours_50, _ = _run_varela_condition(50.0, duration_days=40.0)
    model_ratio = hours_50 / hours_300
    varela_ratio = _VARELA_HOURS_TO_DRYNESS[50.0] / _VARELA_HOURS_TO_DRYNESS[300.0]
    assert model_ratio < varela_ratio, (
        f"model N50/N300 duration ratio {model_ratio:.2f} should stay BELOW Varela's "
        f"real {varela_ratio:.2f} (D-56: the model under-predicts how much severe "
        "nitrogen deficiency slows fermentation, relative to an in-range baseline)"
    )
