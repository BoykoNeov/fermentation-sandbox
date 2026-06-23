"""Tests for the conservation/sanity harness and the benchmark comparator."""

from collections.abc import Callable

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema, VarSpec
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.validation import (
    BENCHMARKS,
    BenchmarkSpec,
    ReferenceSeries,
    assert_conserved,
    assert_nonnegative,
    compare_series,
    total_carbon,
    total_mass,
    total_nitrogen,
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


# -- chemistry-specific conserved-quantity builders ---------------------------


def test_total_mass_conserved_by_toy_fermentation(toy_schema, toy_process):
    ps = ProcessSet(toy_schema, [toy_process])
    y0 = toy_schema.pack({"S": 200.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 200.0))
    # S + E + CO2 is mass-balanced by the Gay-Lussac split (decision D-8).
    assert_conserved(traj, total_mass(toy_schema), rtol=1e-5, atol=1e-6, label="mass")


def test_total_carbon_conserved_by_toy_fermentation(toy_schema, toy_process):
    ps = ProcessSet(toy_schema, [toy_process])
    y0 = toy_schema.pack({"S": 200.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 200.0))
    # Carbon closes because the toy split derives from the same chemistry the
    # check uses — the single-source-of-truth point of decision D-8.
    assert_conserved(traj, total_carbon(toy_schema), rtol=1e-5, atol=1e-6, label="carbon")


def test_total_mass_rejects_multicomponent_sugar():
    # Beer's di-/trisaccharides hydrolyse, pulling solvent water into S+E+CO2, so
    # mass does not close — total_mass refuses rather than silently mislead (D-8).
    with pytest.raises(ValueError, match="hexose/wine check"):
        total_mass(beer_schema())


def test_total_mass_allows_single_sugar_wine():
    # Wine's single hexose slot is fine; the quantity is well-defined.
    mass = total_mass(wine_schema())
    y = wine_schema().pack(
        {"X": 1.0, "S": [200.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    assert mass(y) == pytest.approx(200.0)


def test_total_carbon_value_for_known_wine_state():
    schema = wine_schema()
    y = schema.pack(
        {"X": 2.0, "S": [100.0], "E": 50.0, "N": 0.3, "T": 293.15, "CO2": 20.0, "X_dead": 0.0}
    )
    carbon = total_carbon(schema, biomass_carbon_fraction=0.488)
    expected = (
        100.0 * carbon_mass_fraction("glucose")
        + 50.0 * carbon_mass_fraction("ethanol")
        + 20.0 * carbon_mass_fraction("CO2")
        + 2.0 * 0.488
    )
    assert carbon(y) == pytest.approx(expected)


def test_total_carbon_beer_uses_per_component_fractions():
    schema = beer_schema()
    y = schema.pack(
        {
            "X": 0.0,
            "S": [10.0, 20.0, 30.0],
            "E": 0.0,
            "N": 0.0,
            "T": 293.15,
            "CO2": 0.0,
            "X_dead": 0.0,
        }
    )
    carbon = total_carbon(schema, biomass_carbon_fraction=0.488)
    expected = (
        10.0 * carbon_mass_fraction("glucose")
        + 20.0 * carbon_mass_fraction("maltose")
        + 30.0 * carbon_mass_fraction("maltotriose")
    )
    assert carbon(y) == pytest.approx(expected)


def test_total_carbon_requires_biomass_fraction_when_X_present():
    # Omitting the fraction would silently under-count carbon -> false "conserved".
    with pytest.raises(ValueError, match="biomass_carbon_fraction"):
        total_carbon(wine_schema())


def test_total_carbon_catches_a_carbon_creating_model(toy_schema):
    class CreatesEthanol(Process):
        name = "creates_ethanol"
        tier = Tier.SPECULATIVE
        touches = ("E",)

        def derivatives(self, t, y, schema, params):
            d = schema.zeros()
            d[schema.slice("E")] = 1.0  # ethanol carbon from nowhere
            return d

    ps = ProcessSet(toy_schema, [CreatesEthanol()])
    y0 = toy_schema.pack({"S": 100.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 50.0))
    with pytest.raises(AssertionError, match="not conserved"):
        assert_conserved(traj, total_carbon(toy_schema), label="carbon")


def test_total_nitrogen_conserved_when_growth_uses_store_fraction():
    # The growth Process and the conservation check read the SAME biomass nitrogen
    # fraction from the parameter store, so nitrogen closes (decision D-8).
    params = load_parameters(default_data_dir() / "wine_generic.yaml")
    f_n = params.value("biomass_N_fraction")

    schema = StateSchema([VarSpec("N", "g/L"), VarSpec("X", "g/L")])

    class NitrogenLimitedGrowth(Process):
        name = "growth"
        tier = Tier.PLAUSIBLE
        touches = ("N", "X")

        def derivatives(self, t, y, schema, p):
            d = schema.zeros()
            n = schema.get(y, "N")
            growth = 0.05 * n if n > 0 else 0.0
            d[schema.slice("X")] = growth
            d[schema.slice("N")] = -f_n * growth  # nitrogen drawn into biomass
            return d

    ps = ProcessSet(schema, [NitrogenLimitedGrowth()], strict=True)
    y0 = schema.pack({"N": 0.3, "X": 0.1})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 100.0))

    assert traj.series("X")[-1] > traj.series("X")[0]  # biomass grew
    assert traj.series("N")[-1] < traj.series("N")[0]  # YAN was consumed
    nitrogen = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)
    assert_conserved(traj, nitrogen, rtol=1e-5, atol=1e-6, label="nitrogen")


def test_total_nitrogen_requires_biomass_fraction_when_X_present():
    schema = StateSchema([VarSpec("N", "g/L"), VarSpec("X", "g/L")])
    with pytest.raises(ValueError, match="biomass_nitrogen_fraction"):
        total_nitrogen(schema)
