"""Medium state layouts and the medium registry."""

import dataclasses

import numpy as np
import pytest

from fermentation.core.media import MEDIA, Medium, beer_schema, get_medium, wine_schema

SHARED = (
    "X", "S", "E", "N", "T", "CO2", "X_dead", "Gly", "Byp", "esters", "fusels", "esters_gas",
)  # fmt: skip


def test_wine_schema_has_single_sugar_slot():
    schema = wine_schema()
    assert schema.names == SHARED
    assert schema.spec("S").size == 1
    # X, S(1), E, N, T, CO2, X_dead, Gly, Byp, esters, fusels, esters_gas (D-20)
    assert schema.size == 12


def test_beer_schema_has_three_sequential_sugars():
    schema = beer_schema()
    assert schema.names == SHARED
    s = schema.spec("S")
    assert s.size == 3
    assert s.components == ("glucose", "maltose", "maltotriose")
    # X, S(3), E, N, T, CO2, X_dead, Gly, Byp, esters, fusels, esters_gas (D-20)
    assert schema.size == 14


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
        "Gly": "g/L",
        "Byp": "g/L",
        "esters": "g/L",
        "fusels": "g/L",
        "esters_gas": "g/L",
    }


def test_produced_only_pools_default_to_zero_when_omitted():
    # X_dead/Gly/Byp/esters/fusels/esters_gas are produced-only pools (VarSpec.default=0),
    # so an initial state may omit them; substrate/condition vars stay required (test_state).
    schema = wine_schema()
    arr = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    assert schema.get(arr, "X_dead") == 0.0
    assert schema.get(arr, "Gly") == 0.0
    assert schema.get(arr, "Byp") == 0.0
    assert schema.get(arr, "esters") == 0.0
    assert schema.get(arr, "fusels") == 0.0
    assert schema.get(arr, "esters_gas") == 0.0


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


# -- the registered media now carry the validated-core kinetics + Tier-2 byproducts

# Validated core: growth + fermentative uptake + ethanol-driven cell inactivation,
# with per-rate Arrhenius modifiers. The Luong ethanol wall (ethanol_inhibition) is
# retired from the default media in favour of the cumulative inactivation Process
# (decision D-13).
CORE_PROCESSES = {
    "growth_nitrogen_limited",
    "sugar_uptake_to_ethanol_co2",
    "ethanol_inactivation",
}
# Tier-2 aroma byproducts (Milestone 2, decisions D-18/D-19/D-20): additive aroma
# Processes plus the ester gas-stripping sink (ester_volatilization, D-20: liquid
# esters → the esters_gas headspace pool). Wired in by default but isolable (prime
# directive #3) — disabling them leaves the validated core byte-for-byte.
BYPRODUCT_PROCESSES = {"ester_synthesis", "fusel_alcohols_ehrlich", "ester_volatilization"}
EXPECTED_PROCESSES = CORE_PROCESSES | BYPRODUCT_PROCESSES
EXPECTED_MODIFIERS = {"arrhenius_growth", "arrhenius_uptake"}


@pytest.mark.parametrize("medium", ["wine", "beer"])
def test_registered_media_wire_the_full_kinetic_set(medium):
    # Wine and beer share the same mechanism set — only the sugar vector differs,
    # and beer's sequential uptake lives inside the uptake Process, so no extra
    # Process is needed for it. Both media also carry the Tier-2 byproduct Processes
    # (esters/fusels), which are produced in both wine and beer.
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
