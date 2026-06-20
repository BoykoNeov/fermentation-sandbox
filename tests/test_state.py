"""Tests for the state schema and the ergonomic state-vector view."""

import numpy as np
import pytest

from fermentation.core.state import StateSchema, StateVector, VarSpec


def wine_schema() -> StateSchema:
    return StateSchema(
        [
            VarSpec("X", "g/L", description="viable biomass"),
            VarSpec("S", "g/L", description="fermentable sugar"),
            VarSpec("E", "g/L", description="ethanol"),
            VarSpec("N", "g/L", description="yeast assimilable nitrogen"),
            VarSpec("T", "K", description="temperature"),
            VarSpec("CO2", "g/L", description="evolved CO2"),
        ]
    )


def beer_schema() -> StateSchema:
    return StateSchema(
        [
            VarSpec("X", "g/L"),
            VarSpec(
                "S",
                "g/L",
                size=3,
                components=("glucose", "maltose", "maltotriose"),
                description="sequentially consumed sugars",
            ),
            VarSpec("E", "g/L"),
        ]
    )


def test_size_and_names():
    s = wine_schema()
    assert s.size == 6
    assert s.names == ("X", "S", "E", "N", "T", "CO2")


def test_vector_variable_occupies_contiguous_block():
    s = beer_schema()
    assert s.size == 1 + 3 + 1
    assert s.slice("S") == slice(1, 4)


def test_pack_unpack_roundtrip_scalar():
    s = wine_schema()
    values = {"X": 0.2, "S": 240.0, "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0}
    arr = s.pack(values)
    assert arr.shape == (6,)
    out = s.unpack(arr)
    assert out == pytest.approx(values)


def test_pack_unpack_roundtrip_vector():
    s = beer_schema()
    arr = s.pack({"X": 0.5, "S": [30.0, 60.0, 12.0], "E": 0.0})
    np.testing.assert_allclose(s.get(arr, "S"), [30.0, 60.0, 12.0])
    assert s.get(arr, "X") == 0.5


def test_pack_missing_variable_raises():
    s = wine_schema()
    with pytest.raises(ValueError, match="missing values"):
        s.pack({"X": 1.0})


def test_pack_unknown_variable_raises():
    s = wine_schema()
    full = {"X": 0.2, "S": 240.0, "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0}
    with pytest.raises(ValueError, match="unknown variables"):
        s.pack({**full, "bogus": 1.0})


def test_pack_wrong_vector_length_raises():
    s = beer_schema()
    with pytest.raises(ValueError, match="expects 3"):
        s.pack({"X": 0.5, "S": [1.0, 2.0], "E": 0.0})


def test_duplicate_names_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        StateSchema([VarSpec("X", "g/L"), VarSpec("X", "g/L")])


def test_empty_schema_rejected():
    with pytest.raises(ValueError, match="at least one"):
        StateSchema([])


def test_varspec_component_count_must_match_size():
    with pytest.raises(ValueError, match="component names"):
        VarSpec("S", "g/L", size=2, components=("a", "b", "c"))


def test_unknown_variable_lookup_raises():
    s = wine_schema()
    with pytest.raises(KeyError, match="ethXanol|Unknown"):
        s.slice("ethXanol")


def test_statevector_getset():
    s = wine_schema()
    sv = StateVector.from_values(
        s, {"X": 0.2, "S": 240.0, "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0}
    )
    assert sv["S"] == 240.0
    sv["S"] = 200.0
    assert sv["S"] == 200.0
    assert sv.as_dict()["S"] == 200.0


def test_statevector_defaults_to_zeros():
    s = wine_schema()
    sv = StateVector(schema=s)
    assert np.all(sv.array == 0.0)
    assert sv.array.shape == (6,)
