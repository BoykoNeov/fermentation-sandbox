"""D-59/D-60: independent-data comparison against Palma, Madeira, Mendes-Ferreira &
Sa-Correia (2012), "Impact of assimilable nitrogen availability in glucose uptake
kinetics in Saccharomyces cerevisiae during alcoholic fermentation", Microbial Cell
Factories 11:99. doi:10.1186/1475-2859-11-99, PMC3503800 (CC BY, open access).

This is the project's SECOND independent-data check (after Varela 2004) and the
first against a genuinely DIFFERENT strain: PYCC 4072, not the Coleman/Varela
Prise de Mousse lineage (D-59 Finding 0 -- the two datasets share a strain, so
Palma is the first strain-independent validation point this project has).

**Digitization method (D-59 Finding 3 flagged this as feasible, not yet built):**
Figure 1 panels C (glucose, g/L) and D (ethanol, %v/v) were fetched as the
CC-BY-licensed original image (PMC Open Access S3 mirror,
``PMC3503800.1/1475-2859-11-99-1.jpg``) and read off against a pixel grid
calibrated to the panels' own axis tick marks, at the paper's actual sampling
times (0, 6, 24, 48, 72, 80, 96, 144 h, confirmed from Methods). The paper
reports n=3 replicates with SD error bars on Figure 1, so digitization noise
(~2-5% of axis range, D-59) is well under the gaps this file characterizes.
Two of the paper's three conditions are digitized here -- CF (complete
fermentation, 320 mg N/L) and LF (nitrogen-limited, 90 mg N/L); the third, RF
(LF refed with DAP at 72 h), is a discrete mid-run intervention and is left for
whoever next wants to exercise ``add_dap`` timing fidelity specifically --
out of scope for a first glucose+ethanol pass (D-59's own scoping).

**Conditions (from Methods, not assumed):** synthetic grape-must medium,
glucose as the SOLE sugar (200 g/L, no fructose split -- this maps directly
onto the engine's single-slot wine ``S``, cleaner than Varela's glucose+
fructose must), pH 3.7, pitched at 10^6 CFU/mL (same research-lab dose and the
same order-of-magnitude g/L conversion already used and documented in
test_validation_varela2004.py: ~0.018 g/L), fermented at **20 C -- exactly the
engine's/Coleman's T_ref**, so unlike Varela (28 C) this comparison carries no
Arrhenius extrapolation uncertainty at all.

**Headline finding -- the CF/LF absolute timing gap flips direction from
Varela, and the two gaps (timing, yield) have DIFFERENT, deliberately NOT
conflated explanations -- neither is a fidelity signal:**
the engine reaches CF dryness at ~138 h against Palma's real ~72 h (text:
"consumed approximately after 3 days") -- i.e. the engine runs ~1.9x *slower*
here, the OPPOSITE direction from Varela (where the engine ran ~1.9x *faster*
at 300 mg N/L, 28 C). Cross-checked against Coleman's own reference model (the
same eqs-1-8 reconstruction test_coleman_reconstruction.py uses, re-run here at
Palma's exact S0=200 g/L, N0=320 mg/L, pitch=0.018 g/L, 20 C): it dries at
~140 h, ~1.5% from the engine's ~138 h -- **the engine is faithfully
reproducing Coleman at Palma's own inputs**, so the gap to Palma's real data is
a genuine Coleman-vs-Palma difference, not an engine bug (the exact D-57
argument, transplanted to a new dataset). **The timing gap's best-supported
explanation is strain, not protocol:** at 200 g/L glucose S. cerevisiae is
strongly Crabtree-repressed and ferments even under full aeration, so
respiratory carbon diversion cannot explain a 2x rate difference at this sugar
level -- ruling out "the flask was aerobic" as the timing story. PYCC 4072
(Palma) and Prise de Mousse (Coleman/Varela) are simply different strains with
different fermentation rates; this dataset's whole value is being the first
strain-independent check, so the gap is expected, not a red flag.
**Separately, the yield gap has its own, narrower explanation:** Palma's real
ethanol yield is only ~0.39-0.40 g ethanol / g glucose consumed at both N
levels (CF: 78.9 g/L ethanol on ~199 g/L consumed; LF: ~45.0 g/L on ~120 g/L
consumed) -- well below the anaerobic ~0.46-0.51 g/g range the engine itself
uses (~0.48 here) -- consistent with ethanol evaporating from a shaken,
cotton-stoppered 500 mL Erlenmeyer flask (120 rpm) over a multi-day shake.
Evaporation affects the reported ethanol *level*, not the glucose-consumption
*rate*, so it explains the yield gap only, not the timing gap above. A third,
weaker data point -- Varela's real CF (28 C, warmer) took LONGER (170 h) than
Palma's real CF (20 C, cooler, 72 h) -- shows the two "independent" reference
datasets disagree with each other by ~2.4x, at least as much as either
disagrees with the engine (the same "gap is at or below the reference data's
own discriminating power" shape D-59 reached for the SO2 overshoot); this is
NOT a clean temperature (anti-Arrhenius) comparison, since Varela and Palma are
also different strains -- strain is confounded with temperature here, so no
temperature-specific claim is made. **Absolute CF/LF duration and ethanol
level are therefore characterized here as regression guards (observed value +
margin), not asserted as agreement targets** -- do not tighten these into a
pass/fail against Palma's raw numbers.

**The regime-robust finding -- corroborated on an independent strain:** the
engine still under-predicts how much severe nitrogen limitation suppresses
fermentation PROGRESS, a relative comparison that cancels the yield/evaporation
confound (both conditions share the same flask protocol; only sugar consumed,
not ethanol produced, is compared -- comparing sugar rather than ethanol is
what insulates this finding from the evaporation confound above). Real Palma
LF is only ~60% through its
glucose by 144 h (residual ~80 g/L) and still visibly decelerating (122 to 80
g/L between 96 and 144 h) -- far from dry, NOT a clean plateau (careful:
this is weaker than "arrested", per D-58's overclaim lesson) -- while real
Palma CF is already fully dry by 72 h, giving Palma's own CF:LF
consumed-fraction ratio at 144 h of ~1.66. The engine's LF is much further
along (~79% consumed, residual ~41 g/L) against its own CF's ~99.7%, a ratio of
only ~1.26. **Same direction and same shape as D-56/D-57's Varela finding and
D-59's "model never reproduces arrest" framing, now independent of strain**:
this is the load-bearing, confound-independent signal in this dataset, not the
absolute timing.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import brix_to_sugar_gpl

pytestmark = pytest.mark.benchmark

#: Residual sugar [g/L] defining "dry" -- same threshold as test_milestone1.py
#: and test_validation_varela2004.py.
DRYNESS_GPL = 2.0

#: Palma's synthetic must: 200 g/L glucose as the SOLE sugar (Methods). Scenario.initial
#: only accepts "brix" for wine (brix_to_sugar_gpl * must_fermentable_fraction), and
#: must_fermentable_fraction=0.93 is a REAL-GRAPE-MUST correction that does not strictly
#: apply to a pure-sugar synthetic must -- the same documented, not-fixed approximation
#: test_validation_varela2004.py uses, kept identical here for methodological consistency
#: between the project's two independent-data benchmarks.
_TARGET_SUGAR_GPL = 200.0
_FERMENTABLE_FRACTION = 0.93


def _brix_for_sugar(target_gpl: float) -> float:
    def f(b: float) -> float:
        return brix_to_sugar_gpl(b) * _FERMENTABLE_FRACTION - target_gpl

    return float(brentq(f, 1.0, 40.0))


_BRIX = _brix_for_sugar(_TARGET_SUGAR_GPL)

#: Palma pitched 10^6 CFU/mL -- identical research-lab dose to Varela's, converted via
#: the same ~18 pg/cell S. cerevisiae dry-weight estimate documented in
#: test_validation_varela2004.py (an order-of-magnitude conversion, not exact).
_PITCH_GPL_RESEARCH = 0.018

#: Palma's fermentation temperature -- exactly the engine's/Coleman's T_ref, unlike
#: Varela's 28 C: this comparison carries zero Arrhenius extrapolation uncertainty.
_TEMPERATURE_C = 20.0

#: Digitized from Figure 1 panel C (glucose) -- see module docstring for method.
#: CF's dryness time is read from the text ("consumed approximately after 3 days"),
#: which is sharper than the figure alone can resolve at this axis scale.
_PALMA_CF_HOURS_TO_DRYNESS = 72.0

#: Digitized from Figure 1 panels C/D at t=144 h (the paper's last sampled point).
#: CF is already fully dry by 144 h; LF is far from it.
_PALMA_GLUCOSE_GPL_144H = {320.0: 1.0, 90.0: 80.0}


def _run_palma_condition(yan_mgl: float, duration_days: float) -> tuple[np.ndarray, np.ndarray]:
    """Run the wine scenario at Palma's conditions; return (t_hours, total_sugar_gpl)."""
    scenario = Scenario(
        name=f"palma2012-{yan_mgl:g}mgN",
        medium="wine",
        initial={"brix": _BRIX, "yan_mgl": yan_mgl, "pitch_gpl": _PITCH_GPL_RESEARCH},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=_TEMPERATURE_C)],
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
    return traj.t, total_sugar


def _hours_to_dryness(t_h: np.ndarray, sugar: np.ndarray, duration_days: float) -> float:
    reached = np.where(sugar <= DRYNESS_GPL)[0]
    if not reached.size:
        raise AssertionError(f"never reached dryness within {duration_days} d")
    return float(t_h[reached[0]])


def test_palma2012_cf_duration_characterized():
    # 320 mg N/L (complete fermentation, CF) at 20 C -- the engine's own T_ref, so this
    # is the cleanest possible comparison temperature-wise. The engine still runs ~1.9x
    # SLOWER than Palma's real ~72 h (opposite direction from Varela's ~1.9x too FAST at
    # 300 mg N/L, 28 C). Cross-checked against Coleman's own reference model (module
    # docstring): it dries at ~140 h here, ~1.5% from the engine's own ~138 h -- proof
    # this is a genuine Coleman-vs-Palma difference, not an engine defect -- best
    # explained by strain (PYCC 4072 vs Coleman/Varela's Prise de Mousse), since at this
    # 200 g/L glucose level S. cerevisiae is Crabtree-repressed and ferments even under
    # full aeration (see module docstring: the timing gap and the separate ~0.39 g/g
    # ethanol-yield gap have different causes, deliberately not conflated). This is a
    # regression guard on the engine's OWN characterized behaviour, not a pass/fail
    # against Palma.
    t_h, sugar = _run_palma_condition(320.0, duration_days=15.0)
    hours = _hours_to_dryness(t_h, sugar, 15.0)
    assert 125.0 <= hours <= 150.0, (
        f"CF: hours_to_dryness={hours:.1f} outside characterized [125, 150] h"
    )

    # The gap to Palma's real 72 h should stay in the characterized ~1.7-2.15x band --
    # if this ever falls toward ~1x, that is not a free win, it means something shifted
    # in the engine's fermentation-rate kinetics that should be recorded as a new
    # decision; if it grows past ~2.15x, a regression slowed the model down further.
    gap_ratio = hours / _PALMA_CF_HOURS_TO_DRYNESS
    assert 1.7 <= gap_ratio <= 2.15, (
        f"CF: gap ratio to Palma {gap_ratio:.2f} outside characterized [1.7, 2.15]x"
    )


def test_palma2012_lf_far_from_dry_at_144h_characterized():
    # 90 mg N/L (nitrogen-limited, LF) -- Palma's real LF is only ~60% through its
    # glucose by 144 h (residual ~80 g/L) and still visibly decelerating, i.e. far from
    # dry but NOT a clean plateau (weaker than "arrested" -- D-58's overclaim lesson).
    # The engine's LF is much further along at the same clock time: residual ~41 g/L,
    # ~79% consumed. This is a regression guard on that engine value, characterizing
    # (not closing) the same N-sensitivity shortfall test_palma2012_lf_vs_cf_progress_
    # ratio_understates_palma below measures directly.
    t_h, sugar = _run_palma_condition(90.0, duration_days=20.0)
    s_144 = float(np.interp(144.0, t_h, sugar))
    assert 35.0 <= s_144 <= 48.0, f"LF: S(144h)={s_144:.2f} g/L outside characterized [35, 48] g/L"
    # Sanity: the engine's LF should still be nowhere near Palma's LF residual (~80
    # g/L) -- if it ever rises to meet it, that's the N-sensitivity gap closing, worth
    # a new decision, not a silent re-band.
    assert s_144 < _PALMA_GLUCOSE_GPL_144H[90.0]


def test_palma2012_lf_vs_cf_progress_ratio_understates_palma():
    # The regime-robust, yield/evaporation-confound-independent finding: compare each
    # condition's OWN glucose-consumed fraction at 144 h (both CF and LF share the same
    # flask protocol, so a systematic yield/evaporation confound cancels in the ratio).
    # Real Palma: CF is ~99.5% consumed, LF only ~60% -- ratio ~1.66. The engine's CF is
    # ~99.7% consumed, LF ~79% -- ratio only ~1.26. Same direction as D-56/D-57's Varela
    # finding and D-59's framing, now corroborated on a genuinely independent strain
    # (PYCC 4072, not Coleman/Varela's Prise de Mousse -- D-59 Finding 0): the engine
    # still under-predicts how much severe nitrogen limitation suppresses fermentation
    # progress. If this ratio ever meets or exceeds Palma's real ~1.66, that is a real
    # change in the nitrogen-limitation kinetics worth investigating, not noise.
    t_cf, sugar_cf = _run_palma_condition(320.0, duration_days=15.0)
    t_lf, sugar_lf = _run_palma_condition(90.0, duration_days=20.0)
    consumed_cf = 1.0 - float(np.interp(144.0, t_cf, sugar_cf)) / _TARGET_SUGAR_GPL
    consumed_lf = 1.0 - float(np.interp(144.0, t_lf, sugar_lf)) / _TARGET_SUGAR_GPL
    model_ratio = consumed_cf / consumed_lf

    palma_consumed_cf = 1.0 - _PALMA_GLUCOSE_GPL_144H[320.0] / _TARGET_SUGAR_GPL
    palma_consumed_lf = 1.0 - _PALMA_GLUCOSE_GPL_144H[90.0] / _TARGET_SUGAR_GPL
    palma_ratio = palma_consumed_cf / palma_consumed_lf

    assert model_ratio < palma_ratio, (
        f"model CF:LF consumed-fraction ratio at 144h {model_ratio:.2f} should stay BELOW "
        f"Palma's real {palma_ratio:.2f} (the engine still under-predicts how much severe "
        "nitrogen limitation suppresses fermentation progress, now confirmed on an "
        "independent strain -- see D-56/D-57/D-59/D-60)"
    )
