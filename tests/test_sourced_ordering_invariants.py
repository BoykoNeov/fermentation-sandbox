"""D-171: the sourced ordering invariants, guarded at the NOMINAL — and only there.

WHAT THIS FILE IS NOT. It does not close the breaches that motivated it. The beat was
"guard the five ordering breaches", and guarding a breach means asserting the ordering
*where it breaches* — at the joint band edge, the `E_a_decarb` shape:

    assert decarb.uncertainty.low >= reduction.uncertainty.high   # test_vicinal_diketones

That assertion is available there because the margin is non-negative (both edges meet at
60,000 J/mol). **Every margin in this file is negative**: each pair's bands overlap, so
`low >= high` is FALSE today and cannot be shipped without first moving a band edge. Nine
of the ten parameters involved are `author estimate` / `speculative`, and the D-171
pre-registration committed in advance that *an `author estimate` cannot license a band
narrowing* — narrowing a band so the author's own prose comes true is tuning the band to
the test, the mirror image of tuning the test to the band. So **zero edges moved**, and
the joint-draw inversions below remain live, unguarded, and unguardable until a citation
arrives that moves an edge on its own authority.

What each test here DOES forbid is stated per test, with the breach rate it does NOT
forbid (`feedback-name-guards-for-what-they-forbid` — a guard silent about its own scope
gets read as covering the breach, which is the false assurance D-170 §5 documented).

WHY THESE ARE NEW COVERAGE, NOT DECORATION
(`feedback-mutate-the-premise-before-building-the-guard`). Every ordering below was broken
first and the FULL suite watched, in four rounds:

  round 1  moved each nominal across its ordering.  Three arms — including the KNOWN-RED
           CONTROL — died at collection: the `Parameter` schema enforces
           `low <= value <= high`, so a crossing point outside the band is not a mutation,
           it is a schema rejection. The control "matched its prediction" having never
           executed the assertion it was built to exercise.
  round 2  pre-flighted every mutated tree through the store, and crossed each ordering
           with BOTH members moving onto their own band edges so every value stays legal.
  round 3  paired each RED with an ordering-PRESERVING move of the same size. Every round-2
           RED reproduced under a baseline: they were MAGNITUDE, not the inversion.
  round 4  re-ran the survivors at ASSERT granularity, because a node with four asserts has
           more than one cause.

Result: for every pair below the nominal moved across the ordering and the full 1518-test
suite stayed GREEN, or went red for a reason a baseline reproduced with the ordering
intact. Nothing in the suite forbids any of these inversions. That is what these tests add.

Both sides of every comparison are RECOMPUTED from the store — pinning a literal would
pass for the wrong reason if either band moved (D-154's `k_bound` precedent).
"""

from __future__ import annotations

import pytest

from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters


@pytest.fixture(scope="module")
def wine():
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture(scope="module")
def beer():
    return load_parameters(default_data_dir() / "beer_generic.yaml")


@pytest.fixture(scope="module")
def polymerization():
    return load_parameters(default_data_dir() / "polymerization.yaml")


@pytest.fixture(scope="module")
def thermal():
    return load_parameters(default_data_dir() / "thermal.yaml")


def _inverts(store, lo_name: str, hi_name: str) -> bool:
    """True iff some JOINT draw can invert the ordering — i.e. the two bands overlap.

    ARGUMENT ORDER IS ALWAYS (the one that should be SMALLER, the one that should be LARGER),
    never the order the neighbouring assert happens to be written in. For a `>` assert the
    two therefore look reversed, and that is deliberate:

        assert value(A) > value(B)          # A is the larger
        assert _inverts(store, "B", "A")    # so B is `lo_name` -- NOT a typo

    Do not "tidy" the args to match the assert: it would silently stop tracking the overlap
    and start asserting something else that is also True today. The inversion condition is
    `smaller.uncertainty.high > larger.uncertainty.low` — the sampler can draw the one that
    should be smaller near its top while the other is near its bottom.
    """
    return bool(store[lo_name].uncertainty.high > store[hi_name].uncertainty.low)


# -- the sourced directions, at the nominal ------------------------------------


def test_fusel_activation_energy_ordering_holds_in_each_mediums_own_coupling(wine, beer):
    """ "Warmer makes more higher alcohols" — asserted in the coordinates each medium ferments in.

    RENAMED AND RESCOPED AT D-226, and the rename is the finding. Until D-226 both media
    coupled fusel production to the fermentative flux, integrated production scaled as
    ``arrh(E_a_fusels)/arrh(E_a_uptake)``, and the physical claim "total fusels rise with
    temperature" was therefore the algebraic claim ``E_a_fusels > E_a_uptake``. That is what
    this test asserted, for both media, and it was correct.

    Beer now rides growth EXTENT. The growth Arrhenius factor cancels against the nitrogen
    limit (``int(mu*X*f_growth) dt = YAN / biomass_N_fraction``, a temperature-invariant
    identity), and the five higher alcohols have no stripping sink, so integrated production
    scales as ``arrh(E_a_fusels)`` alone and the SAME physical claim reads ``E_a_fusels > 0``.
    Beer's value re-anchored 70,000 -> 14,100 to hold the 1.2183x rise over 15-25 C the engine
    actually shipped, so it is now below ``E_a_uptake`` — and that is a change of coordinates,
    not a violation.

    **D-19's ordering constraint is a property of the RATE LAW, not of the physics.** Asserting
    the flux-form inequality against a medium that no longer uses the flux form would forbid a
    correct value; deleting it for beer would leave the physical claim unguarded. So each medium
    is asserted in its own coupling, and the coupling each one uses is asserted too — otherwise a
    future re-wiring of beer back to the flux form would leave this test quietly checking the
    wrong inequality.

    FORBIDS: nominal drift of either parameter, in either coupling.
    DOES NOT FORBID, for WINE: the joint-band inversion, which reaches **0.0051 %** of wine
    draws (2e6 triangular samples); the joint-edge margin is -3,000 J/mol.
    """
    from fermentation.core.media import MEDIA

    # The coupling is not assumed. If beer is ever re-wired to the flux producer, this fails
    # here by name rather than leaving the inequality below silently checking the wrong thing.
    beer_producers = {p.name for p in MEDIA["beer"].build_process_set().active}
    wine_producers = {p.name for p in MEDIA["wine"].build_process_set().active}
    assert "fusel_alcohols_ehrlich_growth_coupled" in beer_producers
    assert "fusel_alcohols_ehrlich" in wine_producers

    # WINE -- flux-coupled. Integrated production scales as arrh(E_a_fusels)/arrh(E_a_uptake).
    assert wine.value("E_a_fusels") > wine.value("E_a_uptake"), "wine"
    # The band-scoped half is NOT asserted -- it is false. Recorded so a future edge move is a
    # deliberate re-decision rather than a silent one.
    assert _inverts(wine, "E_a_uptake", "E_a_fusels"), (
        "wine: bands no longer overlap - a band edge moved. If that was deliberate and "
        "citation-backed, promote this pair to the joint-edge assertion "
        "(`E_a_fusels.uncertainty.low >= E_a_uptake.uncertainty.high`) and delete this line."
    )

    # BEER -- growth-extent coupled (D-226). "Rises with T" is E_a_fusels > 0, and it holds
    # across the WHOLE band rather than at the nominal alone, which the flux form never managed.
    assert beer.value("E_a_fusels") > 0.0, "beer"
    assert beer["E_a_fusels"].uncertainty.low > 0.0, "beer band floor"
    # The OTHER sourced ordering this parameter carries -- fusels below esters, so the
    # ester/fusel ratio rises with T -- is now band-scoped for beer where it was nominal-only:
    # E_a_fusels' high edge 50,400 sits below E_a_esters' LOW edge 87,000, so no joint draw can
    # invert it. Under the retired bands (250,000 against 120,000) it inverted on a real slice.
    assert beer.value("E_a_fusels") < beer.value("E_a_esters"), "beer nominal"
    assert not _inverts(beer, "E_a_fusels", "E_a_esters"), (
        "beer: E_a_fusels' band now reaches E_a_esters' - a joint draw can make the higher "
        "alcohols MORE temperature-sensitive than the esters, inverting the sourced claim that "
        "warmth shifts the balance toward esters. D-226 made this pair band-safe; keep it so."
    )


def test_brett_dies_slower_than_o_oeni_at_the_same_molecular_so2(wine):
    """`k_death_brett < k_death_mlf` — Brett is the more SO2-tolerant organism.

    The note declares the split outright: *"DIRECTION is sourced (molecular SO2
    kills/suppresses Brett …); the magnitude is the speculative modelling choice."*

    FORBIDS: nominal drift. Moving `k_death_brett` 0.03 → 0.06, above `k_death_mlf` 0.05,
    left the full suite GREEN (1518 passed) — nothing anywhere asserted this.
    DOES NOT FORBID: the joint-band inversion, **23.05 %** of draws. Joint-edge margin
    −0.06/h: `k_death_brett.high` 0.08 against `k_death_mlf.low` 0.02.
    """
    assert wine.value("k_death_brett") < wine.value("k_death_mlf")
    assert _inverts(wine, "k_death_brett", "k_death_mlf")


def test_mlf_bacteria_are_less_ethanol_tolerant_than_the_yeast(wine):
    """`ethanol_tolerance_mlf < ethanol_tolerance` — why high-alcohol musts challenge MLF.

    Direction sourced (Ribereau-Gayon, Handbook of Enology Vol. 1; G-Alegria 2004): O. oeni
    tolerates ~12-16 % ABV against the yeast's ~15-19 %.

    FORBIDS: nominal drift. Crossing the ordering in-band (mlf → 125, its own high edge;
    yeast → 120, its own low edge) went red on three MLF tests — but all three reproduce
    with IDENTICAL assertion values when `ethanol_tolerance_mlf` moves alone with the
    ordering intact, and the yeast-side move alone is green. The crossing itself is
    unasserted.
    DOES NOT FORBID: the joint-band inversion, **0.036 %** of draws. Joint-edge margin
    −5 g/L: `ethanol_tolerance_mlf.high` 125 against `ethanol_tolerance.low` 120.
    """
    assert wine.value("ethanol_tolerance_mlf") < wine.value("ethanol_tolerance")
    assert _inverts(wine, "ethanol_tolerance_mlf", "ethanol_tolerance")


def test_ethyl_bridging_costs_less_acetaldehyde_per_g_tannin_than_per_g_anthocyanin(
    polymerization,
):
    """`y_acetaldehyde_per_tannin < y_acetaldehyde_per_anthocyanin` — one bridge per PAIR.

    The direction is stoichiometric and cited (Timberlake & Bridle 1976; Es-Safi 1999): one
    acetaldehyde forms one ethylidene bridge per tannin-tannin adduct but per
    tannin-anthocyanin adduct, so per gram of tannin the acetaldehyde cost is roughly halved.

    CAVEAT THIS GUARD CANNOT DISCHARGE: that argument is MOLAR, and both yields ship as MASS
    ratios over lumped pools with, in the notes' own words, *"no clean molar mass"*. The
    bands exist to carry exactly the lumped-molar-mass uncertainty that can invert the
    ordering, so the sourced molar direction does not transfer to the banded mass ratios.

    FORBIDS: nominal drift. Doubling `y_acetaldehyde_per_tannin` 0.06 → 0.12, above 0.09,
    reds only `test_oxidative_cascade_guards::…reproduces_its_trajectory[1y/2y-acetaldehyde]`
    — and HALVING it to 0.03, ordering intact, reds the identical two nodes. Those are
    magnitude pins on a trajectory; nothing asserts the ordering.
    DOES NOT FORBID: the joint-band inversion, **34.21 %** of draws. Joint-edge margin −0.15.
    Carbon closes for any value (split-ledger transfer), so this is a fidelity claim only.
    """
    assert polymerization.value("y_acetaldehyde_per_tannin") < polymerization.value(
        "y_acetaldehyde_per_anthocyanin"
    )
    assert _inverts(polymerization, "y_acetaldehyde_per_tannin", "y_acetaldehyde_per_anthocyanin")


def test_maillard_melanoidins_are_browner_per_mass_than_caramelization_polymers(thermal):
    """`y_a420_per_maillard_melanoidin > y_a420_per_melanoidin` — the Hodge 1953 direction.

    Nitrogenous melanoidins are the dominant brown chromophores of thermally-aged
    sugar+amino-acid systems; the note sets the yield *"per the sourced browner-per-mass
    ordering"*.

    FORBIDS: nominal drift. Moving it 0.8 → 0.3, below the 0.4 sibling, left the full suite
    GREEN (1518 passed).
    DOES NOT FORBID: the joint-band inversion, **16.41 %** of draws. Joint-edge margin −1.3:
    `…maillard_melanoidin.low` 0.2 against `…melanoidin.high` 1.5.

    Both yields are lumped optical yields, not molar extinction coefficients, and this one
    is self-declared CALIBRATED with `k_maillard_browning` against a scenario that is not
    pinned (the D-89 problem) — so the ordering, not either magnitude, is the claim.
    """
    assert thermal.value("y_a420_per_maillard_melanoidin") > thermal.value("y_a420_per_melanoidin")
    # A `>` ordering, so the SIBLING is `lo_name` here: [0.1, 1.5] against [0.2, 3.0] means the
    # caramelan yield can be drawn at 1.5 while the Maillard one sits at 0.2. See _inverts.
    assert _inverts(thermal, "y_a420_per_melanoidin", "y_a420_per_maillard_melanoidin")


# -- deliberately NOT guarded --------------------------------------------------


def test_methanethiol_yield_ordering_is_recorded_as_unsourced_not_guarded(wine):
    """`y_methanethiol < y_h2s_autolysis` is NOT asserted here, on purpose.

    Its note says *"Set BELOW y_h2s_autolysis (2e-5) since reduction skews to H2S"* — but
    unlike the four orderings above, **no source is attached to the skew**. The provenance
    cites biomass methionine content and mercaptan sensory thresholds, neither of which
    fixes the H2S:MeSH branch ratio. Both are order-of-magnitude author estimates whose
    bands ([1e-6, 1e-4] and its sibling) overlap almost completely, inverting on **19.39 %**
    of joint draws.

    Asserting it would pin the author's own prose as if it were sourced — the same error as
    narrowing a band so the prose comes true. The mutation arm confirmed nothing asserts it
    (nominal moved 1e-5 → 3e-5, above the sibling; 1518 passed). It stays unasserted, and
    this test records WHY so the gap is not re-discovered as an oversight.

    What IS asserted: that both remain speculative author estimates. If either ever acquires
    a source that fixes the branch ratio, this test fails and the ordering can be promoted
    to a real guard — or to a joint-edge assertion, if the source also moves an edge.
    """
    for name in ("y_methanethiol", "y_h2s_autolysis"):
        assert wine[name].tier is Tier.SPECULATIVE, name
        assert "author estimate" in wine[name].provenance.source, name
    assert _inverts(wine, "y_methanethiol", "y_h2s_autolysis")
