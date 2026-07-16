"""Tests for Stevens compression — beat 1b slice 2 (decision D-98).

Slice 2 applies a **per-compound** compression curve (``I = OAV ** n``) before slice 1's max
rule. It is the first thing in the sensory layer built on numbers that are *entirely* author
estimates, and the suite is shaped by that: most of it exists to pin what the layer may
**not** claim.

* **The parameters are honest** — all 21 exponents are ``author estimate``/``speculative``,
  cover every aroma pool of both media exactly, and never leak into the chemistry.
* **Compression is threshold-preserving** — ``I > 1`` iff ``OAV > 1`` for any ``n``, so slice
  2 can neither invent nor silence a detectable smell. ``dominant`` is its *only* new
  observable.
* **A global exponent is a no-op** — the executable form of why the exponents are per-compound
  and why a single sourceable one was refused.
* **THE HEADLINE (`test_a_robust_dominance_flip_is_impossible_*`)** — at these bands no flip
  this layer produces can be trusted, on any trajectory, ever. Not a measurement: a
  consequence of the bands overlapping, which they do *because* the values are guesses.
* **The additivity through-line survives** — compression is per-compound, the combination rule
  is untouched, and slice 1 is byte-for-byte unaffected (isolability, prime directive #3).
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.runtime.integrate import Trajectory
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.sensory.compression import (
    FlipVerdict,
    StevensProjector,
    _axis_draws,
    _exponent_key,
    compressed_intensity,
    dominant_flip_sensitivity,
    load_exponents,
)
from fermentation.sensory.descriptors import (
    DescriptorProjector,
    MaxRuleProjector,
    axes_for_medium,
)
from fermentation.sensory.oav import AROMA_COMPOUNDS, load_thresholds, sensory_profile


@pytest.fixture(scope="module")
def thresholds():
    return load_thresholds()


@pytest.fixture(scope="module")
def exponents():
    return load_exponents()


def _traj(schema: StateSchema, pools: dict[str, float], n: int = 4) -> Trajectory:
    """A synthetic constant-in-time trajectory with the named pools set (all else 0)."""
    y: FloatArray = np.zeros((schema.size, n), dtype=np.float64)
    for pool, val in pools.items():
        y[schema.slice(pool), :] = val
    return Trajectory(
        schema=schema,
        t=np.linspace(0.0, 1.0, n),
        y=y,
        success=True,
        message="",
        tier_map=dict.fromkeys(schema.names, Tier.SPECULATIVE),
    )


def _at_oav(thresholds, pool: str, medium: str, oav: float) -> float:
    """The g/L concentration of ``pool`` that reads exactly ``oav`` in ``medium``."""
    return oav * float(thresholds.value(f"threshold_{pool}_{medium}")) / 1_000_000.0


# -- the parameters are honest ------------------------------------------------


def test_every_exponent_is_an_author_estimate_and_speculative(exponents):
    """No entry may claim a source it does not have — the §4.3 firewall at the file level.

    D-98's whole defensibility rests on these 23 numbers never being mistaken for measurements.
    Cain 1969 is cited in `notes` (it orders them and sets the spread's scale) and MUST NOT
    appear in `source`, which would let the values borrow a citation that measured none of
    these compounds, in a matrix Cain never used.

    The count tracks wine's aroma-pool set exactly (one exponent per pool the OAV lens can
    read), which is why D-99's fusel split moved it 21 → 23: `fusels` became `isoamyl_alcohol`
    one-for-one, and `isobutanol` / `2_phenylethanol` joined as wine-thresholded pools. The
    other two higher alcohols (`propanol`, `active_amyl_alcohol`) get NO entry here — they have
    no threshold, so they never reach this lens, and inventing exponents for them would be
    exactly the unforced speculation D-98 exists to forbid.

    D-99 is also where this file came closest to a real source and deliberately did not take
    it: Cain's measured class is n-aliphatic alcohols with the exponent falling by chain
    length, and four of the five new pools are aliphatic alcohols of differing length. That
    supplies an ORDERING, not values — so these stay author estimates with overlapping bands,
    and `test_a_robust_dominance_flip_is_impossible_at_these_bands` still holds.

    D-102 added the 24th, `stevens_n_dms`, when the DMS aroma pool joined the wine lens — the
    contract this count enforces (every aroma pool of every medium needs an exponent) is what
    caught its absence. Its solubility rank was CHECKED against a primary-citing source (PubChem:
    22 g/L at 25 C) rather than recalled, since D-101 shipped a wrong boiling point and D-102
    corrected a fabricated activation energy from exactly that habit.
    """
    assert len(exponents) == 24
    for p in exponents:
        assert p.provenance.source == "author estimate", p.name
        assert p.tier is Tier.SPECULATIVE, p.name
        assert p.uncertainty.low < p.value < p.uncertainty.high, p.name


@pytest.mark.parametrize("medium", ["beer", "wine"])
def test_exponent_coverage_is_exact(exponents, medium):
    """Every aroma pool has an exponent and every exponent has a pool.

    The coverage guard the next beat will need: add an aroma pool without an exponent and
    `StevensProjector` raises rather than silently defaulting to n=1 (which would make that one
    compound uncompressed — an invisible inconsistency).
    """
    pools = {c.pool for c in AROMA_COMPOUNDS[medium]}
    for pool in pools:
        assert _exponent_key(pool) in exponents, pool
    all_pools = {c.pool for m in AROMA_COMPOUNDS for c in AROMA_COMPOUNDS[m]}
    for name in exponents.names:
        assert name.removeprefix("stevens_n_") in all_pools, name


def test_exponents_are_ordered_by_the_documented_solubility_argument(exponents):
    """The file's ordering claim is executable, not just prose in a header.

    `psychophysics.yaml` justifies its ORDER by Cain 1969's solubility-exponent rank
    correlation. That is the only structure the citation supports, so it is the only structure
    pinned: the most water-soluble pool (acetaldehyde, miscible) must sit above the least
    (whiskey lactone, a lipophilic oak lactone), and the spread between them is Cain's reported
    ~2.5x soluble:insoluble ratio rather than a number invented here.
    """
    hi = exponents.value("stevens_n_acetaldehyde")
    lo = exponents.value("stevens_n_whiskey_lactone")
    assert hi > lo
    assert hi / lo == pytest.approx(2.5, abs=0.05)
    # The three esters must order by solubility: ethyl acetate (~80 g/L) steepest, ethyl
    # hexanoate (~0.5 g/L) flattest. This ordering decides the `fruity` axis, so it is the one
    # place the file's argument is load-bearing on an observable.
    assert (
        exponents.value("stevens_n_ethyl_acetate")
        > exponents.value("stevens_n_isoamyl_acetate")
        > exponents.value("stevens_n_ethyl_hexanoate")
    )


def test_exponents_never_reach_the_chemistry(exponents):
    """The §4.2 isolation firewall: no compiled scenario may see an exponent.

    `psychophysics.yaml` loads standalone (never via compile.py's `shared_files`), so the
    chemistry cannot read a perceptual guess even by accident — the same protection beat 1a
    gives the thresholds, and the reason a sensory file may hold 21 estimates at all.
    """
    sc = Scenario(
        name="d98-isolation",
        medium="wine",
        initial={"brix": 22.0, "yan_mgl": 200.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=1.0,
    )
    compiled = compile_scenario(sc, strict=True)
    for name in exponents.names:
        assert name not in compiled.parameters


# -- compression is threshold-preserving --------------------------------------


@pytest.mark.parametrize("n", [1.0, 0.8, 0.65, 0.5, 0.26, 0.2])
def test_intensity_is_one_at_threshold_for_every_exponent(n):
    """`I(1) == 1` for any n — the property that makes `above_threshold` rule-independent.

    Compression re-scales how loud a detectable smell is; it can never make an undetectable
    compound detectable or vice versa. This is why slice 2's ONLY new observable is `dominant`.
    """
    assert compressed_intensity(1.0, n) == pytest.approx(1.0)


@pytest.mark.parametrize("n", [0.8, 0.5, 0.26])
def test_intensity_is_monotone_zero_at_zero_and_compressive(n):
    """More compound never smells like less; nothing smells of an empty pool; ratios squash."""
    assert compressed_intensity(0.0, n) == 0.0
    xs = [0.1, 0.5, 1.0, 2.0, 10.0, 100.0]
    ys = [compressed_intensity(x, n) for x in xs]
    assert ys == sorted(ys)
    # Compressive for n < 1: pulled TOWARD 1 from both sides — the nose squashing ratios.
    assert compressed_intensity(100.0, n) < 100.0
    assert compressed_intensity(0.1, n) > 0.1


def test_above_threshold_is_identical_under_both_rules(thresholds, exponents):
    """Slice 2 flags exactly the same compounds as detectable as slice 1 — pinned, not assumed.

    `StevensProjector` computes `above_threshold` from the compressed intensity; if that ever
    diverged from beat 1a's `OAV > 1` the layer would be inventing or suppressing smells. The
    equivalence is exercised across the threshold in both directions.
    """
    schema = wine_schema()
    for oav in (0.4, 0.99, 1.01, 4.2):
        pools = {"diacetyl": _at_oav(thresholds, "diacetyl", "wine", oav)}
        profile = sensory_profile(_traj(schema, pools), thresholds)
        max_r = MaxRuleProjector().project(profile).readings["buttery"]
        st_r = StevensProjector(exponents).project(profile).readings["buttery"]
        assert st_r.above_threshold == max_r.above_threshold == (oav > 1.0)


def test_a_clean_run_raises_no_false_descriptor(thresholds, exponents):
    """Slice 1's silence invariant survives compression: 0**n == 0, nothing is claimed."""
    profile = sensory_profile(_traj(wine_schema(), {}), thresholds)
    assert StevensProjector(exponents).project(profile).above_threshold() == []


# -- why per-compound, and not global -----------------------------------------


def test_a_global_exponent_cannot_change_dominant(thresholds):
    """THE ARGUMENT FOR PER-COMPOUND EXPONENTS, EXECUTABLE.

    A single global exponent is sourceable — the wine literature's own stated workaround is to
    assign a global value of one — and it is a provable NO-OP: a monotone transform preserves
    argmax, so `dominant` never moves, and `I > 1` iff `OAV > 1`, so `above_threshold` never
    moves either. It would mint a parameter and change nothing, which is why D-98 refused it.
    Per-compound exponents are the ONLY version with an observable, and they are the version
    that cannot be sourced. That dilemma is the whole decision.
    """
    schema = wine_schema()
    pools = {
        "isoamyl_acetate": _at_oav(thresholds, "isoamyl_acetate", "wine", 42.0),
        "ethyl_hexanoate": _at_oav(thresholds, "ethyl_hexanoate", "wine", 79.0),
    }
    profile = sensory_profile(_traj(schema, pools), thresholds)
    for n in (1.0, 0.8, 0.5, 0.3, 0.2):
        contributors = {
            p: compressed_intensity(profile.readings[p].oav, n)
            for p in ("isoamyl_acetate", "ethyl_hexanoate")
        }
        assert max(contributors, key=lambda p: contributors[p]) == "ethyl_hexanoate"


def test_per_compound_exponents_can_flip_dominant(thresholds, exponents):
    """The payload EXISTS: the loudest-by-OAV compound need not be loudest-by-intensity.

    The documented critique of OAV, made executable — a compound detectable at a lower
    concentration can still be the weaker smell at realistic levels. Here ethyl hexanoate has
    ~1.9x isoamyl acetate's OAV and still loses, because its flatter exponent (the less
    water-soluble ester, per the file's ordering) compresses it harder.

    NB this test pins that the mechanism WORKS, not that the answer is TRUSTWORTHY. That the
    same flip is a coin-toss under honest uncertainty is the point of the impossibility test
    below, and the two must be read together.
    """
    pools = {
        "isoamyl_acetate": _at_oav(thresholds, "isoamyl_acetate", "wine", 42.0),
        "ethyl_hexanoate": _at_oav(thresholds, "ethyl_hexanoate", "wine", 79.0),
    }
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    assert MaxRuleProjector().project(profile).readings["fruity"].dominant == "ethyl_hexanoate"
    assert StevensProjector(exponents).project(profile).readings["fruity"].dominant == (
        "isoamyl_acetate"
    )


# -- THE HEADLINE: the payload is unreachable at honest uncertainty -----------


@pytest.mark.parametrize("medium", ["beer", "wine"])
def test_a_robust_dominance_flip_is_impossible_at_these_bands(exponents, medium):
    """D-98'S CENTRAL RESULT — and it is a THEOREM, not a measurement.

    A "flip" is compound j winning although OAV_j < OAV_i (both above threshold). Claim: if
    j and i's exponent bands OVERLAP, j cannot win every draw — so no flip is ever robust.

    Proof. Let v be any value in both bands. The draw n_i == n_j == v is admissible, and it is
    admissible under EITHER sampling scheme — including the order-preserving one, since equal
    exponents violate no strict rank. At that draw compression is a GLOBAL exponent, which is a
    provable no-op (a monotone transform preserves argmax — see
    `test_a_global_exponent_cannot_change_dominant`), so the higher-OAV compound i wins it. A
    neighbourhood of that draw has positive measure, so i wins with positive probability and j
    cannot be unanimous. Hence robustness requires DISJOINT bands. []

    This test asserts none are disjoint — in either medium, on any axis. Therefore **no
    dominance flip this layer can produce is robust, on any trajectory, for any drink, ever**:
    not because these particular guesses are bad, but because the bands are wide *because* they
    are guesses. Narrow enough bands to be disjoint would claim a precision an author estimate
    does not have. An honest band and a trustworthy flip from an estimate are mutually
    exclusive.

    NB this proof deliberately does NOT route through "j wins at min(n_j) against max(n_i)",
    which silently assumes the two are sampled INDEPENDENTLY. That was D-98's original argument
    and it was fragile: the default sampling is order-preserving (Cain's rank is a correlation),
    under which min(n_j) and max(n_i) are not jointly reachable. The equal-exponents argument
    above needs no independence and is why the conclusion survived the correction intact.

    The result is CONDITIONAL on the bands: a real measured exponent with a genuinely narrow
    band could produce a robust flip, and this test is what would then start failing — which is
    exactly the signal that slice 2 had become sourceable. It is written to fail loudly in that
    happy case rather than silently bless it.
    """
    for axis in axes_for_medium(medium):
        for i in axis.pools:
            for j in axis.pools:
                if i >= j:
                    continue
                bi = exponents[_exponent_key(i)].uncertainty
                bj = exponents[_exponent_key(j)].uncertainty
                disjoint = bj.low > bi.high or bi.low > bj.high
                assert not disjoint, (
                    f"{medium}/{axis.name}: {i} and {j} have DISJOINT exponent bands "
                    f"([{bi.low}, {bi.high}] vs [{bj.low}, {bj.high}]). A robust dominance "
                    f"flip is now possible — which is only honest if these stopped being "
                    f"author estimates. See D-98."
                )


def test_the_fruity_flip_is_never_trustworthy_under_its_own_uncertainty(thresholds, exponents):
    """The theorem above, instantiated on the pair that matters most.

    `fruity` is the one axis where compression changes wine's answer (apple -> banana). Under
    the honest order-preserving sampling banana wins a clear MAJORITY at high YAN (~78%) — but
    a majority is not robustness, and the verdict stays CONTESTED. The flip is therefore still a
    statement about `psychophysics.yaml`, not about wine: "cannot say which ester dominates".
    That is why `MaxRuleProjector` remains the default.
    """
    pools = {
        "isoamyl_acetate": _at_oav(thresholds, "isoamyl_acetate", "wine", 42.0),
        "ethyl_hexanoate": _at_oav(thresholds, "ethyl_hexanoate", "wine", 79.0),
    }
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    verdict = dominant_flip_sensitivity(profile, exponents, draws=4000, seed=0)["fruity"]
    assert verdict.contested
    assert not verdict.robust
    # A strong majority, and NOT unanimity — the whole distinction the beat turns on.
    assert 0.6 < verdict.share["isoamyl_acetate"] < 0.95
    assert "CONTESTED" in verdict.summary()


def test_sampling_respects_the_solubility_ordering_it_claims(exponents):
    """THE DONE-CALL CATCH, PINNED: the Monte Carlo may not contradict the file it samples.

    D-98 originally sampled every exponent INDEPENDENTLY while `psychophysics.yaml` asserts —
    and `test_exponents_are_ordered_by_the_documented_solubility_argument` pins — that the
    values are rank-ordered by solubility per Cain. ~28% of draws inverted the two fruity
    esters, i.e. the pass spent a quarter of its evidence on draws the file calls impossible.
    Since Cain's finding IS a rank correlation, the ordering is the best-supported structure the
    citation offers and the absolute values the least: independent sampling kept the weak part
    and discarded the strong one. It was not cosmetic — it moved wine's fruity contest from
    55/45 to 78/22 at YAN 250.
    """
    rng = np.random.default_rng(0)
    pools = ("isoamyl_acetate", "ethyl_hexanoate")  # nominal 0.36 > 0.28
    ordered = _axis_draws(pools, exponents, rng, 3000, preserve_order=True)
    assert np.all(ordered[:, 0] >= ordered[:, 1])
    # ...and the naive mode really does violate it, so the guard is not vacuous.
    naive = _axis_draws(pools, exponents, rng, 3000, preserve_order=False)
    assert np.any(naive[:, 0] < naive[:, 1])


def test_equal_exponents_hand_the_axis_to_the_higher_oav_compound(thresholds, exponents):
    """The theorem's engine, isolated: at n_i == n_j compression is global, hence a no-op.

    This draw is admissible under BOTH sampling schemes (equal values violate no strict rank),
    and it is why overlapping bands forbid a robust flip without any appeal to independence.
    """
    pools = {
        "isoamyl_acetate": _at_oav(thresholds, "isoamyl_acetate", "wine", 42.0),
        "ethyl_hexanoate": _at_oav(thresholds, "ethyl_hexanoate", "wine", 79.0),
    }
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    # 0.30 lies inside BOTH bands ([0.22, 0.54] and [0.20, 0.42]) — an admissible tie.
    for pool in pools:
        band = exponents[_exponent_key(pool)].uncertainty
        assert band.low <= 0.30 <= band.high
    intensities = {p: compressed_intensity(profile.readings[p].oav, 0.30) for p in pools}
    assert max(intensities, key=lambda p: intensities[p]) == "ethyl_hexanoate"


def test_a_wide_oav_gap_is_robust_precisely_because_compression_changed_nothing(
    thresholds, exponents
):
    """The converse, and the sharpest way to state D-98's result.

    Robustness is only ever available where the OAV gap is so wide that no exponent ratio in
    the bands can close it — i.e. exactly where compression did NOT change slice 1's answer.
    Every axis slice 2 leaves alone is trustworthy; every axis it moves is a coin toss. The
    layer is therefore informative only where it is redundant.
    """
    pools = {
        "isoamyl_acetate": _at_oav(thresholds, "isoamyl_acetate", "wine", 1.05),
        "ethyl_hexanoate": _at_oav(thresholds, "ethyl_hexanoate", "wine", 5000.0),
    }
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    verdict = dominant_flip_sensitivity(profile, exponents, draws=2000, seed=0)["fruity"]
    assert verdict.robust
    assert not verdict.contested
    # ...and the max rule already said so, unaided and without a single parameter.
    assert MaxRuleProjector().project(profile).readings["fruity"].dominant == "ethyl_hexanoate"
    assert verdict.nominal == "ethyl_hexanoate"


def test_an_absent_axis_is_reported_silent_not_robust(thresholds, exponents):
    """The guard against maximal confidence in an aroma that does not exist.

    An un-oaked wine has vanillin = whiskey lactone = 0. Since `0 ** n == 0`, every draw ties,
    the tie breaks to the first-listed pool, and a naive verdict reads "vanillin, wins every
    draw" — slice 1's "clean run raises no false descriptor" sin arriving one layer up wearing
    a statistic. `silent` catches it, and `robust` must be False.
    """
    profile = sensory_profile(_traj(wine_schema(), {}), thresholds)
    verdict = dominant_flip_sensitivity(profile, exponents, draws=200, seed=0)["vanilla_oak"]
    assert verdict.silent
    assert not verdict.robust
    assert "silent" in verdict.summary()


def test_sensitivity_is_deterministic_in_its_seed(thresholds, exponents):
    """A verdict that shifted run to run could not be quoted in a decision record."""
    pools = {"isoamyl_acetate": _at_oav(thresholds, "isoamyl_acetate", "wine", 42.0)}
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    a = dominant_flip_sensitivity(profile, exponents, draws=500, seed=7)["fruity"]
    b = dominant_flip_sensitivity(profile, exponents, draws=500, seed=7)["fruity"]
    assert a.share == b.share


# -- the seam, the rule tag, and isolability ----------------------------------


def test_stevens_projector_satisfies_the_seam(exponents):
    """§4.2's swappable-sensory-model requirement: slice 2 arrives THROUGH the seam.

    The best evidence D-95's Protocol was the right shape — the first real alternative
    projector needed no change to it, to slice 1, or to any caller.
    """
    assert isinstance(StevensProjector(exponents), DescriptorProjector)


def test_readings_self_identify_their_rule(thresholds, exponents):
    """`oav` is NOT the same quantity under both projectors, so a reading says which it is.

    Under the max rule the field holds a raw OAV; under Stevens it holds a compressed
    intensity. One field, two quantities is exactly D-96's category error (one pool, two
    molecular identities, in two layers), and D-96's rule was that the honest fix is structure,
    never a disclaimer. Hence `rule`.
    """
    pools = {"diacetyl": _at_oav(thresholds, "diacetyl", "wine", 4.0)}
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    assert MaxRuleProjector().project(profile).readings["buttery"].rule == "max"
    assert StevensProjector(exponents).project(profile).readings["buttery"].rule == "stevens"


def test_slice_one_is_unaffected_by_slice_two(thresholds):
    """Prime directive #3: slice 2 is togglable off, and slice 1 does not know it exists.

    `MaxRuleProjector` reads no exponent and imports nothing from `compression`; deleting
    `psychophysics.yaml` would leave beat 1a + slice 1 byte-for-byte identical. Pinned by
    reading the raw OAV straight back out of the max rule with slice 2 never constructed.
    """
    pools = {"diacetyl": _at_oav(thresholds, "diacetyl", "wine", 4.0)}
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    reading = MaxRuleProjector().project(profile).readings["buttery"]
    assert reading.oav == pytest.approx(4.0)  # uncompressed — the raw OAV, as before D-98
    assert reading.rule == "max"


def test_missing_exponent_raises_rather_than_defaulting(thresholds):
    """A pool without an exponent must fail loudly, never silently read as n=1.

    A silent default would leave one compound uncompressed among 20 compressed ones — an
    inconsistency invisible in the output and capable of manufacturing a `dominant` flip out of
    a missing YAML entry.
    """
    from fermentation.parameters.store import ParameterSet

    pools = {"diacetyl": _at_oav(thresholds, "diacetyl", "wine", 4.0)}
    profile = sensory_profile(_traj(wine_schema(), pools), thresholds)
    with pytest.raises(ValueError, match="no Stevens exponent"):
        StevensProjector(ParameterSet([])).project(profile)


def test_beer_medium_projects_too(thresholds, exponents):
    """Slice 2 is medium-agnostic, as slice 1 is: beer's 9 axes project unchanged."""
    pools = {"diacetyl": _at_oav(thresholds, "diacetyl", "beer", 3.0)}
    profile = sensory_profile(_traj(beer_schema(), pools), thresholds)
    projected = StevensProjector(exponents).project(profile)
    assert projected.readings["buttery"].above_threshold
    assert projected.tier() is Tier.SPECULATIVE


def test_flip_verdict_summary_shapes():
    """The one-line verdicts quoted in D-98 must say what they mean."""
    silent = FlipVerdict("x", "a", {"a": 1.0, "b": 0.0}, silent=True, _unanimous=True)
    assert not silent.robust and "silent" in silent.summary()
    robust = FlipVerdict("x", "a", {"a": 1.0, "b": 0.0}, silent=False, _unanimous=True)
    assert robust.robust and "robust" in robust.summary()
    split = FlipVerdict("x", "a", {"a": 0.6, "b": 0.4}, silent=False, _unanimous=False)
    assert not split.robust and split.contested and "CONTESTED" in split.summary()
