"""Tests for the conservation/sanity harness and the benchmark comparator."""

from collections.abc import Callable

import numpy as np
import pytest

from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.runtime import simulate
from fermentation.validation import (
    BENCHMARKS,
    BenchmarkSpec,
    ReferenceSeries,
    assert_conserved,
    assert_nonnegative,
    compare_series,
)


def _total(schema: StateSchema) -> Callable[[FloatArray], float]:
    sl_s, sl_e, sl_c = schema.slice("S"), schema.slice("E"), schema.slice("CO2")
    return lambda y: float(y[sl_s][0] + y[sl_e][0] + y[sl_c][0])


def test_assert_conserved_passes_for_conserving_model(toy_schema, toy_process):
    ps = ProcessSet(toy_schema, [toy_process])
    y0 = toy_schema.pack({"S": 100.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 100.0))
    assert_conserved(traj, _total(toy_schema), label="mass")


def test_assert_conserved_catches_a_carbon_creating_model(toy_schema):
    class CreatesMass(Process):
        name = "creates_mass"
        tier = Tier.SPECULATIVE
        touches = ("E",)

        def derivatives(self, t, y, schema, params):
            d = schema.zeros()
            d[schema.slice("E")] = 1.0  # ethanol from nowhere
            return d

    ps = ProcessSet(toy_schema, [CreatesMass()])
    y0 = toy_schema.pack({"S": 100.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 50.0))
    with pytest.raises(AssertionError, match="not conserved"):
        assert_conserved(traj, _total(toy_schema), label="mass")


def test_assert_nonnegative_catches_negative_excursion(toy_schema):
    class Drains(Process):
        name = "drains"
        tier = Tier.SPECULATIVE
        touches = ("S",)

        def derivatives(self, t, y, schema, params):
            d = schema.zeros()
            d[schema.slice("S")] = -1.0  # drives sugar below zero, unchecked
            return d

    ps = ProcessSet(toy_schema, [Drains()])
    y0 = toy_schema.pack({"S": 1.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 50.0))
    with pytest.raises(AssertionError, match="went negative"):
        assert_nonnegative(traj, ("S",))


def test_benchmarks_present_and_well_formed():
    assert {"wine_dryness", "beer_attenuation", "co2_peak_then_tail"} <= set(BENCHMARKS)
    for spec in BENCHMARKS.values():
        assert isinstance(spec, BenchmarkSpec)
        assert spec.low <= spec.high
        assert spec.source


def test_benchmark_passes_window():
    wine = BENCHMARKS["wine_dryness"]
    assert wine.passes(12.0)
    assert not wine.passes(20.0)
    assert not wine.passes(3.0)


def test_compare_series_perfect_fit_is_zero_error():
    t = np.linspace(0, 100, 50)
    v = 200.0 * np.exp(-t / 30.0)
    ref = ReferenceSeries(name="sugar", time_h=t, value=v, unit="g/L", source="synthetic")
    fit = compare_series(t, v, ref)
    assert fit.rmse == pytest.approx(0.0, abs=1e-9)
    assert fit.mae == pytest.approx(0.0, abs=1e-9)
    assert fit.n == 50


def test_compare_series_reports_error_for_offset_model():
    t = np.linspace(0, 100, 20)
    ref = ReferenceSeries(name="x", time_h=t, value=np.zeros_like(t), unit="g/L", source="s")
    fit = compare_series(t, np.full_like(t, 5.0), ref)
    assert fit.rmse == pytest.approx(5.0)
    assert fit.mae == pytest.approx(5.0)


def test_compare_series_requires_overlap():
    ref = ReferenceSeries(
        name="x",
        time_h=np.array([500.0, 600.0]),
        value=np.array([1.0, 2.0]),
        unit="g/L",
        source="s",
    )
    with pytest.raises(ValueError, match="No overlap"):
        compare_series(np.array([0.0, 10.0]), np.array([1.0, 1.0]), ref)


def test_reference_series_shape_mismatch_rejected():
    with pytest.raises(ValueError, match="equal shape"):
        ReferenceSeries(
            name="x",
            time_h=np.array([0.0, 1.0]),
            value=np.array([1.0]),
            unit="g/L",
            source="s",
        )
