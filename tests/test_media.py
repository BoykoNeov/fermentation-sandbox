"""Medium state layouts and the medium registry."""

import dataclasses

import numpy as np
import pytest

from fermentation.core.media import MEDIA, Medium, beer_schema, get_medium, wine_schema

SHARED = ("X", "S", "E", "N", "T", "CO2", "X_dead")


def test_wine_schema_has_single_sugar_slot():
    schema = wine_schema()
    assert schema.names == SHARED
    assert schema.spec("S").size == 1
    # X, S(1), E, N, T, CO2, X_dead
    assert schema.size == 7


def test_beer_schema_has_three_sequential_sugars():
    schema = beer_schema()
    assert schema.names == SHARED
    s = schema.spec("S")
    assert s.size == 3
    assert s.components == ("glucose", "maltose", "maltotriose")
    # X, S(3), E, N, T, CO2, X_dead
    assert schema.size == 9


def test_shared_variable_units_are_canonical():
    schema = wine_schema()
    units = {spec.name: spec.unit for spec in schema.specs}
    assert units == {
        "X": "g/L",
        "S": "g/L",
        "E": "g/L",
        "N": "g/L",
        "T": "K",
        "CO2": "g/L",
        "X_dead": "g/L",
    }


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


def test_empty_medium_is_the_no_kinetics_baseline():
    # A Medium with no factories assembles an empty set whose total derivative is
    # identically zero — a valid constant-state baseline. The registered wine/beer
    # media now carry kinetics (see below); this is the property a *bare* Medium
    # still guarantees.
    medium = Medium(name="x", schema=wine_schema())
    pset = medium.build_process_set()
    assert pset.active == ()
    assert pset.active_modifiers == ()
    y = medium.schema.zeros()
    deriv = pset.total_derivatives(0.0, y, {})
    assert np.array_equal(deriv, medium.schema.zeros())


def test_build_process_set_respects_strict_flag():
    pset = Medium(name="x", schema=wine_schema()).build_process_set(strict=True)
    assert pset.strict is True


# -- the registered media now carry the validated-core primary-fermentation kinetics

# Growth + fermentative uptake + ethanol-driven cell inactivation; per-rate
# Arrhenius modifiers. The Luong ethanol wall (ethanol_inhibition) is retired from
# the default media in favour of the cumulative inactivation Process (decision D-13).
EXPECTED_PROCESSES = {
    "growth_nitrogen_limited",
    "sugar_uptake_to_ethanol_co2",
    "ethanol_inactivation",
}
EXPECTED_MODIFIERS = {"arrhenius_growth", "arrhenius_uptake"}


@pytest.mark.parametrize("medium", ["wine", "beer"])
def test_registered_media_wire_the_full_kinetic_set(medium):
    # Wine and beer share the same mechanism set — only the sugar vector differs,
    # and beer's sequential uptake lives inside the uptake Process, so no extra
    # Process is needed for it.
    pset = get_medium(medium).build_process_set(strict=True)
    assert {p.name for p in pset.active} == EXPECTED_PROCESSES
    assert {m.name for m in pset.active_modifiers} == EXPECTED_MODIFIERS


def test_each_build_returns_fresh_kinetic_instances():
    # Factories, not shared instances: two builds must not hand back the same
    # Process objects (a shared mutable Process across runs/media would be a bug).
    a = get_medium("wine").build_process_set()
    b = get_medium("wine").build_process_set()
    a_procs = {p.name: p for p in a.active}
    b_procs = {p.name: p for p in b.active}
    assert all(a_procs[n] is not b_procs[n] for n in a_procs)
