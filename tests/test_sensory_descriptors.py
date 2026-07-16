"""Tests for the descriptor-space projection — beat 1b slice 1 (decision D-95).

The projection groups beat 1a's per-compound OAVs under descriptor words and names the pool
driving each. Coverage:

* **Vocabulary integrity** — no orphan pool (every aroma pool feeds an axis) and no phantom
  pool (no axis names a pool that does not exist); the medium's axis set is *derived*.
* **The max rule** — the beat's thesis: a descriptor reads its loudest contributor, NEVER the
  sum, so the layer inherits beat 1a's refusal to assume perceptual additivity.
* **Definition of done** — monotone; identically silent on a clean run; the ``lumped`` honesty
  flag survives the projection; the tier is the speculative floor *even for a validated input*
  (the non-vacuous, pure-function test).
* **The seam** — an alternative projector satisfies the Protocol and swaps in, making the
  handoff §4.2 "replaceable by an ML model" requirement executable rather than aspirational.

The readout adds no state and touches no chemistry (by construction — no core file changes;
the full suite staying green is the end-to-end proof).
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.runtime.integrate import Trajectory
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.sensory.descriptors import (
    DESCRIPTOR_AXES,
    DescriptorProfile,
    DescriptorProjector,
    DescriptorReading,
    MaxRuleProjector,
    axes_for_medium,
    descriptor_tier,
)
from fermentation.sensory.oav import (
    AROMA_COMPOUNDS,
    SensoryProfile,
    load_thresholds,
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

    The beat-1a helper (``tests/test_sensory_oav.py``), duplicated rather than imported to keep
    the suites independent: constructing the frozen Trajectory directly (no solver) is what
    lets a test pin an input pool's tier.
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


def _at_oav(thresholds, pool: str, medium: str, oav: float) -> float:
    """The g/L concentration of ``pool`` that reads exactly ``oav`` in ``medium``."""
    return oav * float(thresholds.value(f"threshold_{pool}_{medium}")) / 1_000_000.0


def _project(thresholds, schema: StateSchema, pools: dict[str, float], **kw) -> DescriptorProfile:
    return MaxRuleProjector().project(sensory_profile(_traj(schema, pools, **kw), thresholds))


# -- vocabulary integrity -----------------------------------------------------


@pytest.mark.parametrize("medium", ["beer", "wine"])
def test_no_orphan_pool_every_aroma_pool_feeds_an_axis(medium):
    """Every aroma pool the OAV lens reads is projected somewhere — the coverage invariant.

    The guard for future beats: a new aroma pool (the axis has grown 13 Processes since D-67)
    fails loudly here until it is given a descriptor, rather than silently vanishing from the
    projection.
    """
    projected = {pool for axis in axes_for_medium(medium) for pool in axis.pools}
    assert projected == {c.pool for c in AROMA_COMPOUNDS[medium]}


def test_no_phantom_pool_every_axis_pool_is_a_real_aroma_pool():
    """No axis names a pool that does not exist in any medium (catches a typo in membership)."""
    known = {c.pool for m in AROMA_COMPOUNDS for c in AROMA_COMPOUNDS[m]}
    for axis in DESCRIPTOR_AXES:
        for pool in axis.pools:
            assert pool in known, f"{axis.name} names unknown pool {pool!r}"


def test_axis_set_is_derived_from_the_medium():
    """Beer's vocabulary is a strict subset of wine's, narrowed pool-wise — never declared.

    Wine → 15 axes over 20 pools; beer → 9 over 10. Beer cannot report the grape/Brett-only
    words by construction, and its shared axes are narrowed to the pools beer actually has.
    """
    beer = {a.name: a.pools for a in axes_for_medium("beer")}
    wine = {a.name: a.pools for a in axes_for_medium("wine")}
    assert len(wine) == 15
    assert len(beer) == 9
    assert set(beer) < set(wine)  # strict subset
    assert set(wine) - set(beer) == {
        "barnyard",  # ethylphenols — Brett, wine-only
        "floral_honey",  # phenylacetaldehyde — D-75
        "cooked_potato",  # methional — D-75
        "malty",  # the three D-87 thermal Strecker aldehydes
        "curry_maple",  # sotolon — D-87
        # dms — D-102. Its OWN axis rather than a member of `sulfidic`: DMS is a volatile sulfur
        # compound, but `sulfidic` means "rotten egg / drains / reductive" and DMS is the
        # truffle/black-olive note, POSITIVE near its threshold. Folding it in would repeat D-99's
        # finding exactly (2-phenylethanol booked `solventy` when it is ROSE) — a grouping that is
        # self-consistent and still wrong. Chemical family is not a perceptual axis.
        "truffle_olive",
    }
    # Shared axes are narrowed to the medium's pools: mercaptans + ethylguaiacols are wine-only.
    assert beer["sulfidic"] == ("h2s",)
    assert wine["sulfidic"] == ("h2s", "mercaptans")
    assert beer["smoky"] == ("guaiacol",)
    assert wine["smoky"] == ("guaiacol", "ethylguaiacols")


def test_ethylguaiacols_feeds_two_axes():
    """The many-to-many case: 4-EG genuinely smells both smoky and clove-spicy."""
    axes = {a.name: a.pools for a in axes_for_medium("wine")}
    assert "ethylguaiacols" in axes["smoky"]
    assert "ethylguaiacols" in axes["clove_spice"]


# -- the max rule (the beat's thesis) -----------------------------------------


def test_descriptor_reads_the_max_never_the_sum(thresholds):
    """`malty` with three contributors reads the LOUDEST (4.2), not their sum (5.6).

    THE load-bearing test. Summing OAVs assumes perceptual additivity — contested, and
    explicitly refused one layer down (D-67, ``SensoryProfile``); a projector that summed would
    reintroduce it silently. Max asserts nothing beyond "this compound is 4.2× over threshold".
    """
    profile = _project(
        thresholds,
        wine_schema(),
        {
            "3_methylbutanal": _at_oav(thresholds, "3_methylbutanal", "wine", 4.2),
            "2_methylbutanal": _at_oav(thresholds, "2_methylbutanal", "wine", 1.1),
            "2_methylpropanal": _at_oav(thresholds, "2_methylpropanal", "wine", 0.3),
        },
    )
    malty = profile.readings["malty"]
    assert malty.oav == pytest.approx(4.2)
    assert malty.oav != pytest.approx(5.6)  # the sum — explicitly NOT what we report
    assert malty.dominant == "3_methylbutanal"
    assert malty.contributors == pytest.approx(
        {"2_methylbutanal": 1.1, "3_methylbutanal": 4.2, "2_methylpropanal": 0.3}
    )


def test_sub_threshold_contributors_never_sum_into_a_perceived_descriptor(thresholds):
    """Three pools each at OAV 0.4 leave `malty` silent — a sum (1.2) would fake a smell.

    The consequence of the max rule that matters most: the projection cannot invent an
    above-threshold descriptor that no single compound justifies.
    """
    profile = _project(
        thresholds,
        wine_schema(),
        {
            p: _at_oav(thresholds, p, "wine", 0.4)
            for p in ("2_methylbutanal", "3_methylbutanal", "2_methylpropanal")
        },
    )
    assert profile.readings["malty"].oav == pytest.approx(0.4)
    assert profile.readings["malty"].above_threshold is False
    assert "malty" not in profile.above_threshold()


def test_dominant_names_the_argmax_and_tracks_it(thresholds):
    """`dominant` follows the loudest contributor as the balance shifts."""
    schema = wine_schema()
    h2s_loud = _project(
        thresholds,
        schema,
        {
            "h2s": _at_oav(thresholds, "h2s", "wine", 3.0),
            "mercaptans": _at_oav(thresholds, "mercaptans", "wine", 0.5),
        },
    )
    assert h2s_loud.readings["sulfidic"].dominant == "h2s"
    assert h2s_loud.readings["sulfidic"].oav == pytest.approx(3.0)

    merc_loud = _project(
        thresholds,
        schema,
        {
            "h2s": _at_oav(thresholds, "h2s", "wine", 0.5),
            "mercaptans": _at_oav(thresholds, "mercaptans", "wine", 3.0),
        },
    )
    assert merc_loud.readings["sulfidic"].dominant == "mercaptans"


def test_above_threshold_regroups_the_pool_level_flags(thresholds):
    """A descriptor clears iff one of its pools does — a regrouping, not a new claim.

    Pins the honest framing: under the max rule this layer adds vocabulary + attribution, not
    new above-threshold information (the D-80 "mechanism, not a behaviour change" precedent).
    """
    schema = wine_schema()
    traj = _traj(
        schema,
        {
            "ethylphenols": _at_oav(thresholds, "ethylphenols", "wine", 5.0),
            "diacetyl": _at_oav(thresholds, "diacetyl", "wine", 0.2),
        },
    )
    oav_profile = sensory_profile(traj, thresholds)
    desc_profile = MaxRuleProjector().project(oav_profile)

    assert oav_profile.above_threshold() == ["ethylphenols"]
    assert desc_profile.above_threshold() == ["barnyard"]
    # Every descriptor that clears is backed by a pool that clears, and vice versa.
    cleared_pools = set(oav_profile.above_threshold())
    for name in desc_profile.above_threshold():
        assert cleared_pools & set(desc_profile.readings[name].contributors)


# -- definition of done -------------------------------------------------------


def test_clean_run_raises_no_false_descriptor(thresholds):
    """All pools 0 ⇒ every axis reads 0 and nothing is perceived."""
    for schema in (beer_schema(), wine_schema()):
        profile = _project(thresholds, schema, {})
        assert profile.above_threshold() == []
        for reading in profile:
            assert reading.oav == 0.0
            assert reading.above_threshold is False


def test_descriptor_is_monotone_in_its_pools(thresholds):
    """Raising any pool never lowers its descriptor (max is monotone non-decreasing)."""
    schema = wine_schema()
    lo = _project(thresholds, schema, {"isoamyl_acetate": 1.0e-3})
    hi = _project(thresholds, schema, {"isoamyl_acetate": 2.0e-3})
    assert 0.0 < lo.readings["fruity"].oav < hi.readings["fruity"].oav

    # Raising a NON-dominant contributor leaves the descriptor unmoved (it is still not loudest)
    # — monotone non-decreasing, not strictly increasing. That is the max rule, working.
    base = _project(
        thresholds,
        schema,
        {
            "h2s": _at_oav(thresholds, "h2s", "wine", 3.0),
            "mercaptans": _at_oav(thresholds, "mercaptans", "wine", 0.5),
        },
    )
    bumped = _project(
        thresholds,
        schema,
        {
            "h2s": _at_oav(thresholds, "h2s", "wine", 3.0),
            "mercaptans": _at_oav(thresholds, "mercaptans", "wine", 1.5),
        },
    )
    assert bumped.readings["sulfidic"].oav == pytest.approx(base.readings["sulfidic"].oav)


def test_lumped_flag_propagates_from_the_dominant_contributor(thresholds):
    """D-66's fixed-lump-composition caveat survives the projection.

    `sulfidic` mixes h2s (clean, single-molecule) with mercaptans (a lumped pool read against
    methanethiol). The descriptor inherits the assumption exactly when the lumped pool is the
    one driving it — the honesty cost must not evaporate on crossing a layer.
    """
    schema = wine_schema()
    merc_loud = _project(
        thresholds,
        schema,
        {
            "h2s": _at_oav(thresholds, "h2s", "wine", 0.5),
            "mercaptans": _at_oav(thresholds, "mercaptans", "wine", 3.0),
        },
    )
    assert merc_loud.readings["sulfidic"].lumped is True

    h2s_loud = _project(
        thresholds,
        schema,
        {
            "h2s": _at_oav(thresholds, "h2s", "wine", 3.0),
            "mercaptans": _at_oav(thresholds, "mercaptans", "wine", 0.5),
        },
    )
    assert h2s_loud.readings["sulfidic"].lumped is False
    # A clean single-molecule axis never claims a lump assumption.
    assert h2s_loud.readings["buttery"].lumped is False
    # `fruity` read True here until D-96, when the lumped `esters` pool became two
    # single-molecule pools (isoamyl_acetate + ethyl_hexanoate). Both contributors are now real
    # molecules, so the axis carries no fixed-composition assumption at all — the caveat did not
    # move, it stopped existing.
    assert h2s_loud.readings["fruity"].lumped is False


def test_no_contributor_to_solventy_is_lumped_whoever_dominates(thresholds):
    """`solventy` carries NO lump caveat from any contributor — the D-99 converse (was D-96).

    THIS TEST ASSERTED THE OPPOSITE UNTIL D-99, and the flip is the decision's whole point.
    Between D-96 and D-99 `solventy` was the mixed case — `ethyl_acetate` (single-molecule,
    exact) beside `fusels` (a real lump read against isoamyl alcohol) — so it pinned that the
    axis inherits the lump assumption exactly when the lump is loudest. D-99 split `fusels`
    into five single-molecule pools, so every contributor to this axis is now its own molecule
    and there is no lump left to inherit. The assertion could not be repaired to say what it
    said; it is rewritten to pin what replaced it.

    The lumped-propagation MACHINERY is still tested — `sulfidic` (h2s + mercaptans) remains a
    genuinely mixed axis, because wine's `mercaptans` is now the last lump in the project. This
    test is that one's converse guard: it fires if anyone re-points a solventy pool at a
    molecule it is not made of and reaches for the `lumped` flag to excuse it.
    """
    schema = wine_schema()
    ea_loud = _project(
        thresholds,
        schema,
        {
            "ethyl_acetate": _at_oav(thresholds, "ethyl_acetate", "wine", 4.0),
            "isoamyl_alcohol": _at_oav(thresholds, "isoamyl_alcohol", "wine", 0.5),
            "isobutanol": _at_oav(thresholds, "isobutanol", "wine", 0.5),
        },
    )
    assert ea_loud.readings["solventy"].dominant == "ethyl_acetate"
    assert ea_loud.readings["solventy"].lumped is False

    # MAX still swaps the dominant — the axis is not merely un-lumped, it is still discriminating.
    for loud in ("isoamyl_alcohol", "isobutanol"):
        quiet = {
            p: _at_oav(thresholds, p, "wine", 0.5)
            for p in ("ethyl_acetate", "isoamyl_alcohol", "isobutanol")
        }
        profile = _project(
            thresholds, schema, {**quiet, loud: _at_oav(thresholds, loud, "wine", 4.0)}
        )
        assert profile.readings["solventy"].dominant == loud
        # The point: a higher-alcohol-dominated solventy no longer drags a lump caveat with it.
        assert profile.readings["solventy"].lumped is False


def test_descriptor_tier_is_the_speculative_floor_even_for_a_validated_input():
    """The pure tier rule: all-VALIDATED contributors still floor to speculative.

    The NON-VACUOUS floor test (D-67's advisor caught the vacuous form of exactly this: every
    aroma pool's OAV is already speculative, so asserting through a real profile is a
    tautology). Testing the pure function with validated inputs proves the PROJECTION itself
    caps confidence — grouping compounds under a word and naming one dominant is a perceptual
    claim no threshold measurement backs.
    """
    assert descriptor_tier([Tier.VALIDATED]) is Tier.SPECULATIVE
    assert descriptor_tier([Tier.VALIDATED, Tier.VALIDATED]) is Tier.SPECULATIVE
    assert descriptor_tier([Tier.PLAUSIBLE, Tier.VALIDATED]) is Tier.SPECULATIVE
    assert descriptor_tier([]) is Tier.SPECULATIVE


def test_profile_reading_tier_is_speculative_even_when_pool_is_validated(thresholds):
    """End-to-end floor: an untouched (VALIDATED) pool still reads speculative after projection."""
    schema = wine_schema()
    tm = dict.fromkeys(schema.names, Tier.SPECULATIVE)
    tm["diacetyl"] = Tier.VALIDATED
    profile = _project(thresholds, schema, {"diacetyl": 5.0e-4}, tier_map=tm)
    assert profile.readings["buttery"].tier is Tier.SPECULATIVE
    assert profile.tier() is Tier.SPECULATIVE


def test_projection_inherits_the_taste_exclusions(thresholds):
    """No descriptor can reach a TASTE pool — the odor/taste split is inherited from beat 1a.

    `iso_alpha`/IBU (bitterness, D-64) and `ellagitannin`/astringency (D-78) are tastes, not
    odors, and are absent from `SensoryProfile` — so consuming that profile makes it
    structurally impossible for a bitterness to leak into an aroma descriptor.
    """
    projected = {pool for axis in DESCRIPTOR_AXES for pool in axis.pools}
    for taste in ("iso_alpha", "ellagitannin"):
        assert taste not in projected
    # Nor the set-and-hold oak ceiling slots (read by OakExtraction; not aroma pools).
    assert not any(p.endswith("_ceiling") for p in projected)


def test_dominant_pools_view(thresholds):
    profile = _project(
        thresholds,
        wine_schema(),
        {"isoamyl_acetate": _at_oav(thresholds, "isoamyl_acetate", "wine", 2.0)},
    )
    assert profile.dominant_pools()["fruity"] == "isoamyl_acetate"


# -- the swappable seam (handoff §4.2) ----------------------------------------


def test_max_rule_projector_satisfies_the_seam():
    assert isinstance(MaxRuleProjector(), DescriptorProjector)


def test_an_alternative_projector_swaps_in(thresholds):
    """A stand-in 'panel model' satisfies the Protocol and is consumed identically.

    This is the handoff §4.2 requirement made executable: "a clean seam so it can later be
    replaced by an ML model trained on real sensory-panel data". The stub emits a per-descriptor
    INTENSITY (the thing a panel-trained model would, and a max cannot) without touching beat
    1a, the chemistry, or any caller — which is the whole point of the seam.
    """

    class SummingStub:
        """NOT a proposal — deliberately the additive rule this layer refuses, to prove that
        swapping the rule is a one-class change confined to this seam."""

        def project(self, profile: SensoryProfile) -> DescriptorProfile:
            readings: dict[str, DescriptorReading] = {}
            for axis in axes_for_medium(profile.medium):
                contributors = {p: profile.readings[p].oav for p in axis.pools}
                total = sum(contributors.values())
                dominant = max(contributors, key=lambda p: contributors[p])
                readings[axis.name] = DescriptorReading(
                    descriptor=axis.name,
                    contributors=contributors,
                    dominant=dominant,
                    oav=total,
                    above_threshold=total > 1.0,
                    lumped=profile.readings[dominant].lumped,
                    tier=descriptor_tier(profile.readings[p].tier for p in axis.pools),
                )
            return DescriptorProfile(
                medium=profile.medium, time_index=profile.time_index, readings=readings
            )

    stub: DescriptorProjector = SummingStub()
    assert isinstance(stub, DescriptorProjector)

    traj = _traj(
        wine_schema(),
        {
            p: _at_oav(thresholds, p, "wine", 0.4)
            for p in ("2_methylbutanal", "3_methylbutanal", "2_methylpropanal")
        },
    )
    oav_profile = sensory_profile(traj, thresholds)
    # Same input, different rule: the stub's additive malty clears threshold where v1's does not.
    assert MaxRuleProjector().project(oav_profile).readings["malty"].oav == pytest.approx(0.4)
    assert stub.project(oav_profile).readings["malty"].oav == pytest.approx(1.2)


# -- D-96 regression: the fruity OAV is PHYSICAL on a real run ----------------


def _finished_wine_profile(
    thresholds, **initial_extra: float
) -> tuple[SensoryProfile, DescriptorProfile]:
    """Compile + run a real wine to dryness and return its (SensoryProfile, DescriptorProfile).

    A REAL solver run, not a synthetic trajectory: the D-96 defect lived in the interaction
    between the carbon ledger, the parameter files and the OAV lens, so only an end-to-end run
    exercises the wiring that produced 761 in the first place.
    """
    sc = Scenario(
        name="d96-regression",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25, **initial_extra},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=40.0,
    )
    cs = compile_scenario(sc, strict=True)
    scheduled = cs.run(t_eval=np.linspace(0.0, 40.0 * 24.0, 200))
    assert scheduled.success, scheduled.message
    profile = sensory_profile(scheduled.as_trajectory(), thresholds)
    return profile, MaxRuleProjector().project(profile)


@pytest.mark.parametrize("brett", [False, True], ids=["clean-wine", "brett-wine"])
def test_fruity_oav_is_physical_on_a_real_wine_run(thresholds, brett):
    """THE D-96 regression: `fruity` reads tens, not 761 — pinned on a real finished wine.

    This is the test the whole decision exists to make possible, and the one thing every other
    D-96 test could not catch. The structural tests (5:2 exactness, per-ester carbon weighting,
    the lumped-marker pair) all pin MECHANISM; none of them pins the OUTCOME. Reintroduce a
    ledger/lens mismatch — carbon-weight a pool as one molecule while reading it against
    another's threshold — or fat-finger a threshold by two orders of magnitude, and every one of
    them still passes while `fruity` climbs back to a non-physical 761.

    The bound is the physical claim, not a golden number: a wine's fruity esters genuinely run
    OAV in the tens-to-low-hundreds (Guth 1997), so anything in [1, 200] is defensible and 761 —
    which implied ~23 mg/L of isoamyl acetate against a real ceiling of ~1-3 — is not. Run on
    BOTH a clean wine and the **Brett** wine whose 761 started this: Brett touches the phenol
    pools, not the esters, so its fruity reading must be equally physical. That the pre-D-96
    number came from the Brett scenario specifically is why it is parametrized here rather than
    checked on a generic wine alone.
    """
    extra = {"hydroxycinnamic_gpl": 0.15, "brett_pitch_gpl": 0.4} if brett else {}
    profile, descriptors = _finished_wine_profile(thresholds, **extra)

    fruity = descriptors.readings["fruity"]
    assert 1.0 < fruity.oav < 200.0, (
        f"fruity OAV {fruity.oav:.1f} is outside the physical band — the pre-D-96 lumped "
        f"reading was 761 (dominant {fruity.dominant})"
    )
    # The number must come from a REAL molecule, and MAX must be choosing between real
    # candidates — before D-96 `fruity` read one lumped pool, so `dominant` was vacuous.
    assert fruity.dominant in ("isoamyl_acetate", "ethyl_hexanoate")
    assert set(fruity.contributors) == {"isoamyl_acetate", "ethyl_hexanoate"}
    # No ester carries a fixed-lump-composition assumption any more (D-96), so the axis cannot
    # inherit one: the honesty flag must be off *because it is false*, not by omission.
    assert fruity.lumped is False

    # Each ester's own OAV is physical too — the axis max could hide one absurd contributor.
    for pool in ("isoamyl_acetate", "ethyl_hexanoate"):
        assert 0.0 < profile.readings[pool].oav < 200.0, pool
    # Ethyl acetate is the bulk ester by mass yet must read a MODEST OAV: it is the pool whose
    # mass the pre-D-96 lump borrowed while wearing isoamyl acetate's ~300x lower threshold —
    # the exact swap that manufactured 761. Reading single digits here is that seam staying shut.
    assert 0.0 < profile.readings["ethyl_acetate"].oav < 20.0
