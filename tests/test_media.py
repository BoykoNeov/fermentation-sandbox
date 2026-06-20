"""Medium state layouts and the medium registry."""

import dataclasses

import numpy as np
import pytest

from fermentation.core.media import MEDIA, Medium, beer_schema, get_medium, wine_schema

SHARED = ("X", "S", "E", "N", "T", "CO2")


def test_wine_schema_has_single_sugar_slot():
    schema = wine_schema()
    assert schema.names == SHARED
    assert schema.spec("S").size == 1
    # X, S(1), E, N, T, CO2
    assert schema.size == 6


def test_beer_schema_has_three_sequential_sugars():
    schema = beer_schema()
    assert schema.names == SHARED
    s = schema.spec("S")
    assert s.size == 3
    assert s.components == ("glucose", "maltose", "maltotriose")
    # X, S(3), E, N, T, CO2
    assert schema.size == 8


def test_shared_variable_units_are_canonical():
    schema = wine_schema()
    units = {spec.name: spec.unit for spec in schema.specs}
    assert units == {"X": "g/L", "S": "g/L", "E": "g/L", "N": "g/L", "T": "K", "CO2": "g/L"}


def test_registry_exposes_wine_and_beer():
    assert set(MEDIA) == {"wine", "beer"}
    assert get_medium("wine").schema.spec("S").size == 1
    assert get_medium("beer").schema.spec("S").size == 3


def test_unknown_medium_raises():
    with pytest.raises(KeyError, match="Unknown medium 'cider'"):
        get_medium("cider")


def test_medium_is_frozen():
    with pytest.raises(dataclasses.FrozenInstanceError):
        get_medium("wine").name = "beer"  # type: ignore[misc]


def test_empty_process_set_is_the_no_kinetics_baseline():
    # No Processes are wired yet, so the assembled set is empty and the total
    # derivative is identically zero — a valid constant-state baseline.
    medium = get_medium("wine")
    pset = medium.build_process_set()
    assert pset.active == ()
    y = medium.schema.zeros()
    deriv = pset.total_derivatives(0.0, y, {})
    assert np.array_equal(deriv, medium.schema.zeros())


def test_build_process_set_respects_strict_flag():
    pset = Medium(name="x", schema=wine_schema()).build_process_set(strict=True)
    assert pset.strict is True
