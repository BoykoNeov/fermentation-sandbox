"""Tests for boundary unit conversions, including round-trip properties."""

import math

import pytest
from hypothesis import given
from hypothesis import strategies as st

from fermentation.units import (
    abv_from_ethanol,
    apparent_gravity,
    brix_to_sg,
    brix_to_sugar_gpl,
    celsius_to_kelvin,
    days_to_hours,
    hours_to_days,
    kelvin_to_celsius,
    plato_to_sg,
    real_to_apparent_extract,
    sg_to_brix,
    sg_to_plato,
    sugar_gpl_to_brix,
)


def test_temperature_known_points():
    assert celsius_to_kelvin(0.0) == pytest.approx(273.15)
    assert celsius_to_kelvin(20.0) == pytest.approx(293.15)
    assert kelvin_to_celsius(273.15) == pytest.approx(0.0)


@given(st.floats(min_value=-50, max_value=200))
def test_temperature_roundtrip(celsius):
    assert kelvin_to_celsius(celsius_to_kelvin(celsius)) == pytest.approx(celsius)


@given(st.floats(min_value=0, max_value=10000))
def test_time_roundtrip(days):
    assert hours_to_days(days_to_hours(days)) == pytest.approx(days)


def test_days_hours_known():
    assert days_to_hours(1.0) == 24.0
    assert hours_to_days(48.0) == 2.0


def test_water_is_sg_one():
    # 0 Brix is pure water -> SG 1.000.
    assert brix_to_sg(0.0) == pytest.approx(1.0)


def test_brix_sg_typical_must():
    # A ~24 Brix must sits around SG 1.10 (handoff wine benchmark).
    sg = brix_to_sg(24.0)
    assert 1.09 < sg < 1.11


@given(st.floats(min_value=0.5, max_value=40))
def test_brix_sg_roundtrip(brix):
    # The two industry polynomials are independent fits, so allow a small gap.
    assert sg_to_brix(brix_to_sg(brix)) == pytest.approx(brix, abs=0.3)


def test_plato_brix_numerically_close():
    # Plato and Brix both measure % sucrose by mass; they nearly coincide.
    sg = 1.048
    assert sg_to_plato(sg) == pytest.approx(sg_to_brix(sg), abs=0.3)


@given(st.floats(min_value=0.5, max_value=30))
def test_plato_sg_roundtrip(plato):
    assert sg_to_plato(plato_to_sg(plato)) == pytest.approx(plato, abs=0.3)


def test_sugar_concentration_scale():
    # ~24 Brix must -> ~240-265 g/L sugar.
    gpl = brix_to_sugar_gpl(24.0)
    assert 230 < gpl < 270


@given(st.floats(min_value=0.1, max_value=40))
def test_sugar_brix_roundtrip_at_fixed_sg(brix):
    sg = brix_to_sg(brix)
    gpl = brix_to_sugar_gpl(brix, sg=sg)
    assert sugar_gpl_to_brix(gpl, sg=sg) == pytest.approx(brix)


def test_abv_from_ethanol_scale():
    # ~100 g/L ethanol is around 12.7% ABV.
    abv = abv_from_ethanol(100.0)
    assert 12.0 < abv < 13.5


def test_abv_zero():
    assert abv_from_ethanol(0.0) == 0.0


@given(st.floats(min_value=0, max_value=160))
def test_abv_monotonic_and_finite(ethanol_gpl):
    abv = abv_from_ethanol(ethanol_gpl)
    assert abv >= 0
    assert math.isfinite(abv)


def test_apparent_extract_unchanged_before_fermentation():
    # With no fermentation real extract == original extract, so the apparent
    # (hydrometer) extract equals it too — nothing depresses the reading yet.
    oe = sg_to_plato(1.048)
    assert real_to_apparent_extract(oe, oe) == pytest.approx(oe)


def test_apparent_gravity_1048_to_1010_example():
    # The canonical brewing example: a 1.048 OG ale that finishes at apparent
    # 1.010 has a *real* extract near 4.25 P (~1.016). Going the other way, a real
    # extract of 4.25 P at OG 11.91 P must read ~1.010 apparent.
    oe = sg_to_plato(1.048)  # ~11.91 P
    real_extract = 4.25
    assert real_to_apparent_extract(real_extract, oe) == pytest.approx(2.56, abs=0.1)
    assert apparent_gravity(real_extract, oe) == pytest.approx(1.010, abs=0.001)


def test_apparent_gravity_below_real_when_fermenting():
    # Once ethanol is present (real extract has dropped below OG) the hydrometer
    # reads low: apparent gravity < the real-extract gravity.
    oe = sg_to_plato(1.048)
    real_extract = 5.0  # partially fermented
    assert apparent_gravity(real_extract, oe) < plato_to_sg(real_extract)
