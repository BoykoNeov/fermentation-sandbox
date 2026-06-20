"""Milestone 1 acceptance benchmarks (handoff section 2.2).

These are the §2.2 wine and beer acceptance criteria, encoded as executable
tests *now* so the kinetics are test-driven against them later. They are skipped
until the validated-core Processes exist; the skip reason names what unblocks
each one. Do not delete or weaken these to make CI green — implement the model
until they pass (handoff section 2.2: "encode the benchmark, then iterate").
"""

import pytest

from fermentation.validation import BENCHMARKS

pytestmark = pytest.mark.benchmark

KINETICS_PENDING = "Milestone 1: primary-fermentation kinetics not implemented yet"


@pytest.mark.skip(reason=KINETICS_PENDING)
def test_wine_24brix_ferments_to_dryness_in_10_to_14_days():
    # Build a ~240 g/L sugar, 20 C, nitrogen-limited wine scenario; integrate;
    # assert days_to_dryness lies in BENCHMARKS["wine_dryness"] window.
    spec = BENCHMARKS["wine_dryness"]
    days_to_dryness = ...  # noqa: F841  (filled in Milestone 1)
    assert spec.passes(days_to_dryness)  # type: ignore[arg-type]


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
