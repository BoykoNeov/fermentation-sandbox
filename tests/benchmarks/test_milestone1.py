"""Milestone 1 acceptance benchmarks (handoff section 2.2).

These are the §2.2 wine and beer acceptance criteria, encoded as executable
tests *now* so the kinetics are test-driven against them later. They are skipped
until the validated-core Processes exist; the skip reason names what unblocks
each one. Do not delete or weaken these to make CI green — implement the model
until they pass (handoff section 2.2: "encode the benchmark, then iterate").
"""

import numpy as np
import pytest

from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import BENCHMARKS

pytestmark = pytest.mark.benchmark

KINETICS_PENDING = "Milestone 1: primary-fermentation kinetics not implemented yet"

#: Residual sugar [g/L] defining "dry". Below Coleman's stuck threshold (0.4 %
#: w/v ~ 4 g/L); 2 g/L is a solidly dry wine.
DRYNESS_GPL = 2.0


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


@pytest.mark.skip(reason=KINETICS_PENDING)
def test_beer_1048_og_attenuates_in_5_to_7_days():
    spec = BENCHMARKS["beer_attenuation"]
    days_to_target = ...  # noqa: F841
    assert spec.passes(days_to_target)  # type: ignore[arg-type]


@pytest.mark.skip(reason=KINETICS_PENDING)
def test_co2_integral_tracks_sugar_consumed():
    spec = BENCHMARKS["co2_peak_then_tail"]
    ratio = ...  # noqa: F841
    assert spec.passes(ratio)  # type: ignore[arg-type]


@pytest.mark.skip(reason=KINETICS_PENDING)
def test_lower_temperature_is_slower_but_cleaner():
    # Directional check: lower T -> longer time-to-dryness, fewer fusel/ester
    # byproducts. Qualitative once byproduct Processes exist (Tier 2).
    raise NotImplementedError
