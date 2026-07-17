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
    lo = _traj(beer_schema(), {"isoamyl_acetate": 1.0e-3})
    hi = _traj(beer_schema(), {"isoamyl_acetate": 2.0e-3})
    o_lo = float(oav_series(lo, thresholds, "isoamyl_acetate")[-1])
    o_hi = float(oav_series(hi, thresholds, "isoamyl_acetate")[-1])
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
    """Beer profiles the 7 common + 5 oak pools (12); wine adds 2 fusel + 9 wine-only (23, D-99)."""
    beer = sensory_profile(_traj(beer_schema(), {}), thresholds)
    wine = sensory_profile(_traj(wine_schema(), {}), thresholds)
    # The 5 oak extractives (D-77 four + furaneol/caramel D-94) are SHARED by both media (barrel-
    # beer oak, D-86 — the oak axis is a wood property). The oak *ceiling* slots are NOT aroma
    # pools, so they must NOT appear.
    oak = {"whiskey_lactone", "vanillin", "guaiacol", "eugenol", "furaneol"}
    # The lumped `esters` pool became three single-molecule ester pools at D-96, taking the
    # common set from 5 to 7. D-99 split the lumped `fusels` pool into five single-molecule
    # higher alcohols but left this set at 7 and beer's total at 12: only `isoamyl_alcohol`
    # has a sourced in-matrix BEER threshold (Meilgaard 1975), so it alone replaces `fusels`
    # one-for-one here. Beer's other four higher alcohols are real modelled pools that carry
    # no OAV — present in `beer_schema()` and deliberately absent from this set.
    assert set(beer.readings) == {
        "diacetyl", "acetaldehyde", "h2s",
        "ethyl_acetate", "isoamyl_acetate", "ethyl_hexanoate", "isoamyl_alcohol",
    } | oak  # fmt: skip
    # The chemistry/sensory asymmetry D-99 introduced, asserted rather than described: beer
    # FERMENTS these four and the sim MODELS them, but no beer threshold is sourced for any of
    # them, so none may reach the OAV lens. This is the guard against a future edit "tidying"
    # the aroma set to match the schema — the mismatch is the honest position, not an oversight.
    for chemistry_only in ("propanol", "isobutanol", "active_amyl_alcohol", "2_phenylethanol"):
        assert chemistry_only in beer_schema(), chemistry_only
        assert chemistry_only not in beer.readings, chemistry_only
    assert set(wine.readings) == set(beer.readings) | {
        # The two higher alcohols whose MOLECULE is in both media but whose THRESHOLD is
        # sourced only for wine (Guth 1997, 10% aq. ethanol — the same paper and matrix as
        # isoamyl alcohol's). Beer ferments both; no beer threshold exists for either (D-99).
        "isobutanol",
        "2_phenylethanol",
        "ethylphenols",
        "ethylguaiacols",
        "methanethiol",
        "methional",
        "phenylacetaldehyde",
        # The four non-oxidative THERMAL Strecker aldehyde/sotolon pools (decision D-87).
        "2_methylbutanal",
        "3_methylbutanal",
        "2_methylpropanal",
        "sotolon",
        # The aged-wine truffle/black-olive odorant (decision D-102) — a single-molecule pool
        # (NOT lumped) whose precursor is grape-borne, so it accumulates on bottle aging with no
        # lees. Wine-only: beer's DMS is real but arrives by other routes (wort-boil cleavage
        # before pitch; yeast DMSO reduction), so no beer threshold is claimed here.
        "dms",
    }
    # The set-and-hold ceiling slots are read by OakExtraction but are not aroma compounds.
    assert "vanillin_ceiling" not in wine.readings
    assert "vanillin_ceiling" not in beer.readings
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


def test_no_pool_is_lumped_and_any_future_lump_must_declare_itself(thresholds):
    """The project carries NO lumped pool (D-110) — and the D-66 caveat still binds if one returns.

    **This test used to assert the opposite, and both halves are deliberate.** Through D-109 it
    opened ``assert lumped, "the lumped set should not be empty — mercaptans is still a lump"``,
    with a comment instructing that if a future decision split the last lump the honest response
    was to *delete this test*. D-110 retired that lump — so the instruction came due, and it is
    **not** followed, because D-101 (later, and reasoning about this exact retirement) prescribed
    replacing the tripwire with a *positive assertion* instead. Deleting would have thrown away a
    live guard to satisfy a stale note: the D-66 rule (*a lump must declare its
    fixed-composition assumption in its threshold provenance*) is **still true**, it simply has
    no instances today. A rule with no instances is dormant, not wrong.

    So the assertion inverts rather than vanishes: emptiness is now the **claim**, checked in
    both directions. ``assert not lumped`` fails the day a pool re-lumps — and if that day comes,
    the loop below is already waiting to demand its caveat.
    """
    lumped = {
        (c.pool, medium)
        for medium, compounds in AROMA_COMPOUNDS.items()
        for c in compounds
        if c.lumped
    }
    # The lump-composition risk class (D-66) is CLOSED: `esters` speciated at D-96, `fusels` at
    # D-99, and wine's `mercaptans` — the last one — retired at D-110, where the flag turned out
    # to be FALSE rather than merely retired: nothing in the model produces ethanethiol or any
    # other thiol, so the pool was one molecule under a plural name. Renaming it `methanethiol`
    # made the slot, the carbon weight, the threshold and the Stevens exponent name one compound.
    assert not lumped, (
        f"a pool re-lumped: {sorted(lumped)} — if that is deliberate, it must carry the D-66 "
        "fixed-composition caveat in its threshold provenance (the loop below enforces it), and "
        "this assertion should become the list of known lumps rather than be deleted"
    )
    # Dormant, not dead: the D-66 contract for any lump that ever returns. Vacuous today BY
    # CONSTRUCTION, which the assertion above is what makes safe — without it, this loop would
    # silently pass over an empty set and the caveat rule could rot unnoticed (the `0 == 0`
    # trap D-105 named, and D-109 hit again in its own tripwire).
    for pool, medium in lumped:  # pragma: no cover — no lumped pools exist since D-110
        notes = thresholds[f"threshold_{pool}_{medium}"].provenance.notes.strip().lower()
        # The convention is a LEADING declaration, checked with startswith rather than a bare
        # substring: prose elsewhere in a note may legitimately discuss lumping (a
        # single-molecule note may explain what it was split *from*), and only the opening
        # marker is the claim.
        assert notes.startswith("lumped pool"), (pool, medium)
        assert "composition" in notes, (pool, medium)


def test_single_molecule_thresholds_do_not_declare_a_lump_caveat(thresholds):
    """The converse of the D-66 flag, and the D-96 guard: the caveat only where it is TRUE.

    Before D-96 the ``esters`` pool was carbon-weighted as *ethyl acetate* but read against
    *isoamyl acetate*'s threshold, and the ``lumped`` flag was made to carry that split
    identity — which it could not honestly do: the resulting OAV was non-physical (~761 for a
    wine, implying ~23 mg/L isoamyl acetate against a real ceiling of ~1–3), not merely
    uncertain. A ``lumped`` flag can excuse a *coarse* reading; it cannot excuse reading a pool
    against a molecule the pool is not made of.

    Paired with the D-66 test above this pins the marker in **both** directions, so the two
    tests together make the flag and the provenance impossible to drift apart. If someone
    re-points a single-molecule pool at a different molecule's threshold and reaches for the
    lump caveat to excuse it, they must flip ``lumped`` — which trips the D-66 test unless the
    pool really is a lump. The honest fix is another pool, never another disclaimer.
    """
    for medium, compounds in AROMA_COMPOUNDS.items():
        for c in compounds:
            if c.lumped:
                continue
            notes = thresholds[f"threshold_{c.pool}_{medium}"].provenance.notes.strip().lower()
            assert not notes.startswith("lumped pool"), (
                f"{c.pool} ({medium}) is flagged lumped=False but its threshold declares a "
                "fixed-lump-composition assumption"
            )
