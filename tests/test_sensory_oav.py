"""Tests for the OAV sensory readout — the speculative aroma lens (decision D-67).

OAV = concentration / perception-threshold, a pure readout over aroma-active compounds the
chemistry already tracks. Coverage:

* **Plumbing / golden** — a known concentration at 2× its threshold reads OAV ≈ 2 (validates
  the arithmetic and the µg/L↔g/L unit crossing, NOT the threshold magnitude).
* **Definition of done** — OAV monotone in its pool; identically 0 when the pool is 0; the
  tier is the speculative floor *even for a validated input* (the non-vacuous, pure-function
  test); the reported compound set matches the medium.
* **Medium/matrix** — medium inferred from the schema signature; a wine-only pool on beer is a
  clear error; every threshold loads speculative with its measurement matrix recorded.
* **Isolation** — the readout adds no state and touches no chemistry (by construction; the
  full suite staying green is the end-to-end proof).
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.runtime.integrate import Trajectory
from fermentation.sensory.oav import (
    AROMA_COMPOUNDS,
    load_thresholds,
    medium_of,
    oav_series,
    oav_tier,
    sensory_profile,
)


@pytest.fixture(scope="module")
def thresholds():
    return load_thresholds()


def _traj(
    schema: StateSchema,
    pools: dict[str, float],
    *,
    tier_map: dict[str, Tier] | None = None,
    n: int = 4,
) -> Trajectory:
    """A synthetic constant-in-time trajectory with the named pools set (all else 0).

    Constructing the frozen :class:`Trajectory` directly (no solver) lets a test pin an input
    pool's tier — needed to prove the speculative floor caps even a VALIDATED input.
    """
    y: FloatArray = np.zeros((schema.size, n), dtype=np.float64)
    for pool, val in pools.items():
        y[schema.slice(pool), :] = val
    tm = tier_map if tier_map is not None else dict.fromkeys(schema.names, Tier.SPECULATIVE)
    return Trajectory(
        schema=schema,
        t=np.linspace(0.0, 1.0, n),
        y=y,
        success=True,
        message="",
        tier_map=tm,
    )


# -- golden / plumbing --------------------------------------------------------


def test_golden_diacetyl_plumbing_reads_twice_threshold(thresholds):
    """Diacetyl dosed at 2× the beer threshold reads OAV ≈ 2 — checks arithmetic + µg/L↔g/L.

    This validates the PLUMBING (ratio and the unit crossing), NOT the threshold number: OAV
    ≈ 2 here is true by construction. It catches a wrong ``ugl_to_gpl`` factor, nothing more.
    """
    thr_ugl = thresholds.value("threshold_diacetyl_beer")  # 100 µg/L
    conc_gpl = 2.0 * thr_ugl / 1_000_000.0  # 2× threshold, expressed back in g/L
    traj = _traj(beer_schema(), {"diacetyl": conc_gpl})
    series = oav_series(traj, thresholds, "diacetyl")
    assert series == pytest.approx(np.full_like(series, 2.0))


# -- definition of done -------------------------------------------------------


def test_oav_is_zero_when_pool_is_zero(thresholds):
    """A clean run raises no false aroma — every OAV is identically 0."""
    traj = _traj(beer_schema(), {})  # all pools 0
    for c in AROMA_COMPOUNDS["beer"]:
        assert np.all(oav_series(traj, thresholds, c.pool) == 0.0)


def test_oav_is_monotone_increasing_in_its_pool(thresholds):
    """More of a compound ⇒ strictly higher OAV (a positive scaling of the pool)."""
    lo = _traj(beer_schema(), {"esters": 1.0e-3})
    hi = _traj(beer_schema(), {"esters": 2.0e-3})
    o_lo = float(oav_series(lo, thresholds, "esters")[-1])
    o_hi = float(oav_series(hi, thresholds, "esters")[-1])
    assert 0.0 < o_lo < o_hi
    assert o_hi == pytest.approx(2.0 * o_lo, rel=1e-12)


def test_oav_tier_is_the_speculative_floor_even_for_a_validated_input():
    """The pure tier rule: a VALIDATED input + VALIDATED threshold still floors to speculative.

    This is the NON-VACUOUS floor test. Feeding a real trajectory would be a tautology (every
    aroma pool is produced by a speculative/plausible Process, so the input is never
    validated); testing the pure function with a validated input proves the mapping itself
    caps confidence.
    """
    assert oav_tier(Tier.VALIDATED, Tier.VALIDATED) is Tier.SPECULATIVE
    assert oav_tier(Tier.PLAUSIBLE, Tier.VALIDATED) is Tier.SPECULATIVE
    assert oav_tier(Tier.VALIDATED, Tier.SPECULATIVE) is Tier.SPECULATIVE


def test_profile_reading_tier_is_speculative_even_when_pool_is_validated(thresholds):
    """End-to-end floor: an *untouched* (VALIDATED) pool still reads speculative in the profile.

    Pins ``diacetyl``'s trajectory tier to VALIDATED (as it would be if no Process touched it),
    then confirms the profile reading is speculative — the floor applied through the real API,
    not just the pure function above.
    """
    schema = beer_schema()
    tm = dict.fromkeys(schema.names, Tier.SPECULATIVE)
    tm["diacetyl"] = Tier.VALIDATED
    traj = _traj(schema, {"diacetyl": 5.0e-4}, tier_map=tm)
    profile = sensory_profile(traj, thresholds)
    assert profile.readings["diacetyl"].tier is Tier.SPECULATIVE
    assert profile.tier() is Tier.SPECULATIVE


def test_profile_compound_set_matches_the_medium(thresholds):
    """Beer profiles the 5 common pools; wine adds the 5 wine-only pools — never mismatched."""
    beer = sensory_profile(_traj(beer_schema(), {}), thresholds)
    wine = sensory_profile(_traj(wine_schema(), {}), thresholds)
    assert set(beer.readings) == {"diacetyl", "acetaldehyde", "h2s", "esters", "fusels"}
    assert set(wine.readings) == set(beer.readings) | {
        "ethylphenols",
        "ethylguaiacols",
        "mercaptans",
        "methional",
        "phenylacetaldehyde",
    }
    # The wine-only pools are absent from the beer profile (not silently zero-filled).
    assert "ethylphenols" not in beer.readings
    assert beer.medium == "beer"
    assert wine.medium == "wine"


def test_above_threshold_flag_tracks_oav_one(thresholds):
    """`above_threshold` lists exactly the pools whose OAV exceeds 1 at the chosen time."""
    schema = beer_schema()
    # diacetyl at 3× its threshold (above); h2s far below its threshold.
    diacetyl_gpl = 3.0 * thresholds.value("threshold_diacetyl_beer") / 1_000_000.0
    h2s_gpl = 0.1 * thresholds.value("threshold_h2s_beer") / 1_000_000.0
    traj = _traj(schema, {"diacetyl": diacetyl_gpl, "h2s": h2s_gpl})
    profile = sensory_profile(traj, thresholds)
    assert profile.readings["diacetyl"].above_threshold is True
    assert profile.readings["diacetyl"].oav == pytest.approx(3.0)
    assert profile.readings["h2s"].above_threshold is False
    assert profile.above_threshold() == ["diacetyl"]


# -- medium inference & matrix ------------------------------------------------


def test_medium_inference_from_schema_signature():
    assert medium_of(beer_schema()) == "beer"
    assert medium_of(wine_schema()) == "wine"


def test_medium_inference_raises_on_unknown_schema():
    from fermentation.core.state import VarSpec

    bare = StateSchema([VarSpec("X", "g/L"), VarSpec("E", "g/L")])
    with pytest.raises(ValueError, match="cannot infer medium"):
        medium_of(bare)


def test_wine_only_pool_on_beer_is_a_clear_error(thresholds):
    """A beer schema has no ethylphenols slot at all — reject with a clear message."""
    traj = _traj(beer_schema(), {})
    with pytest.raises(ValueError, match="not an aroma-active pool for medium 'beer'"):
        oav_series(traj, thresholds, "ethylphenols")


# -- thresholds provenance ----------------------------------------------------


def test_every_threshold_is_speculative_with_a_recorded_matrix(thresholds):
    """All thresholds load speculative and record the matrix they were measured in (§4.2)."""
    keys = [f"threshold_{c.pool}_beer" for c in AROMA_COMPOUNDS["beer"]]
    keys += [f"threshold_{c.pool}_wine" for c in AROMA_COMPOUNDS["wine"]]
    for key in keys:
        p = thresholds[key]
        assert p.tier is Tier.SPECULATIVE, key
        assert p.unit == "ug/L", key
        assert p.provenance.conditions.strip(), key  # measurement matrix recorded


def test_lumped_thresholds_flag_the_fixed_composition_assumption(thresholds):
    """The three lumped reps carry the 'lump composition' honesty cost in their notes (D-66)."""
    for key in ("threshold_esters_beer", "threshold_fusels_beer", "threshold_mercaptans_wine"):
        notes = thresholds[key].provenance.notes.lower()
        assert "lump" in notes and "composition" in notes, key
