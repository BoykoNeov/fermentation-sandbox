"""The n-protic speciation branch — added BESIDE the fast paths, never in place of them (D-178).

Beer's charge balance carries **citric acid**, which is triprotic, and the three speciation
helpers (:func:`mean_charge`, :func:`neutral_fraction`, :func:`bisulfite_fraction`) all raised
above two pKas. D-178 adds a general branch dispatched to only at ``len(pkas) >= 3``.

The load-bearing constraint is that **wine stays bitwise**. Wine's acids are all mono- or
diprotic, so the guarantee is structural rather than numerical: the two fast paths are still
the code that runs for them. These tests pin that structurally (by breaking the general branch
and showing wine's cases are unaffected), not by eyeballing a tolerance — a tolerance test
would pass just as happily if the fast paths had been deleted.

They also pin the thing that makes the "added beside" wording necessary in the first place:
the general branch is algebraically equal to the fast paths and **not** bitwise equal, because
it accumulates the same terms in a different order.
"""

import math

import pytest

from fermentation.core import acidbase
from fermentation.core.acidbase import (
    _polyprotic_terms,
    bisulfite_fraction,
    mean_charge,
    neutral_fraction,
)

# A pH grid spanning the bracket the solver actually uses, deliberately including values
# far outside any beverage window (the BDF Jacobian probe reaches them, D-46).
PH_GRID = [0.0, 1.0, 2.0, 3.0, 3.4, 4.0, 4.3, 4.76, 5.5, 7.0, 9.0, 11.0, 14.0]

# Wine's own acids, as shipped (values from acidbase.yaml's neighbourhood -- the exact
# numbers do not matter here, only that they are the 1- and 2-pKa shapes).
MONOPROTIC = (3.86,)
DIPROTIC = (3.04, 4.37)
# Citric acid, CRC 25 C / I=0. Triprotic, and TWO of its pKas sit inside beer's pH window.
CITRIC = (3.128, 4.761, 6.396)
# Phosphate, for the record: the acid this beat was opened on. Triprotic too, but both of
# its relevant pKas sit far OUTSIDE beer's window -- see test_phosphate_is_nearly_inert.
PHOSPHATE = (2.15, 7.20, 12.35)


def _general_mean_charge(h: float, pkas: tuple[float, ...]) -> float:
    terms, denom = _polyprotic_terms(h, pkas)
    return sum(k * term for k, term in enumerate(terms)) / denom


# -- the fast paths are STILL TAKEN (structural, not numerical) -----------------


@pytest.mark.parametrize("pkas", [MONOPROTIC, DIPROTIC])
@pytest.mark.parametrize("fn", [mean_charge, neutral_fraction, bisulfite_fraction])
def test_mono_and_diprotic_never_reach_the_general_branch(monkeypatch, fn, pkas):
    """Break the general branch; the 1- and 2-pKa cases must not notice.

    This is the wine-stays-bitwise guarantee stated as a structural fact. A tolerance
    comparison could not distinguish "the fast path ran" from "the general branch ran and
    agreed to 1e-15", which is precisely the distinction that matters here.
    """

    def _boom(*args, **kwargs):
        raise AssertionError("the general n-protic branch was reached for a mono/diprotic acid")

    monkeypatch.setattr(acidbase, "_polyprotic_terms", _boom)
    for ph in PH_GRID:
        fn(10.0 ** (-ph), pkas)  # must not raise


def test_the_general_branch_is_reached_for_three_pkas(monkeypatch):
    """The positive control for the test above — break it and the triprotic case DOES notice.

    Without this, `test_mono_and_diprotic_never_reach_the_general_branch` would still pass if
    the monkeypatch target were wrong (a self-sealing green,
    [[feedback-a-null-result-needs-a-positive-control]]).
    """

    def _boom(*args, **kwargs):
        raise AssertionError("reached")

    monkeypatch.setattr(acidbase, "_polyprotic_terms", _boom)
    with pytest.raises(AssertionError, match="reached"):
        mean_charge(10.0**-4.3, CITRIC)


# -- the general branch AGREES with the fast paths (algebraically, not bitwise) --


@pytest.mark.parametrize("pkas", [MONOPROTIC, DIPROTIC])
def test_the_general_branch_reproduces_the_fast_paths_to_tolerance(pkas):
    for ph in PH_GRID:
        h = 10.0 ** (-ph)
        terms, denom = _polyprotic_terms(h, pkas)
        assert _general_mean_charge(h, pkas) == pytest.approx(mean_charge(h, pkas), rel=1e-12)
        assert terms[0] / denom == pytest.approx(neutral_fraction(h, pkas), rel=1e-12)
        assert terms[1] / denom == pytest.approx(bisulfite_fraction(h, pkas), rel=1e-12)


def test_the_agreement_is_not_bitwise_and_that_is_why_the_fast_paths_stay():
    """At least one grid point must differ in the last bits.

    If this ever goes green-by-equality the "added beside, never a rewrite" justification has
    quietly become untestable — and someone will delete the fast paths as duplication. The
    assertion is deliberately weak (SOME point differs), because WHICH point differs is a
    property of the floating-point unit, not of this module.
    """
    differs = any(
        _general_mean_charge(10.0 ** (-ph), DIPROTIC) != mean_charge(10.0 ** (-ph), DIPROTIC)
        for ph in PH_GRID
    )
    assert differs, "general branch is bitwise-identical to the diprotic fast path"


# -- the triprotic case itself -------------------------------------------------


@pytest.mark.parametrize("pkas", [CITRIC, PHOSPHATE])
def test_species_fractions_are_a_partition(pkas):
    """Every species fraction in [0,1], summing to exactly 1 — the distribution is closed."""
    for ph in PH_GRID:
        terms, denom = _polyprotic_terms(10.0 ** (-ph), pkas)
        fractions = [t / denom for t in terms]
        assert all(0.0 <= f <= 1.0 for f in fractions)
        assert math.fsum(fractions) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("pkas", [CITRIC, PHOSPHATE])
def test_mean_charge_is_bounded_and_monotone_in_ph(pkas):
    """Bounded by [0, n] and rising with pH — what `solve_ph`'s single-root argument needs."""
    charges = [mean_charge(10.0 ** (-ph), pkas) for ph in PH_GRID]
    assert all(0.0 <= c <= len(pkas) for c in charges)
    assert charges == sorted(charges)
    assert charges[0] == pytest.approx(0.0, abs=1e-2)  # fully protonated at pH 0
    # NOT fully ionised at pH 14 in general: phosphate's pKa3 is 12.35, so only ~98 % of the
    # third proton is off at the top of the bracket. The bound is the honest one -- within one
    # proton's worth of the maximum -- rather than a tolerance tuned until citric passed.
    assert charges[-1] > len(pkas) - 1.0
    assert charges[-1] == pytest.approx(float(len(pkas)), abs=5e-2 if pkas[-1] < 10.0 else 3e-2)


def test_citric_at_its_own_pkas_splits_adjacent_species_evenly():
    """At pH = pKa_k the species either side of that step are equal — the textbook check.

    An independent property (no reference to the implementation's algebra), so it catches a
    transcription error in the cumulative products that a self-consistency test would not.
    """
    for k, pka in enumerate(CITRIC, start=1):
        terms, _ = _polyprotic_terms(10.0 ** (-pka), CITRIC)
        assert terms[k - 1] == pytest.approx(terms[k], rel=1e-12)


def test_phosphate_is_nearly_inert_across_beers_ph_window():
    """The measured reason this beat did NOT ship phosphate (decision D-178).

    Phosphate's pKas are 2.15 / 7.20 — beer sits at ~4.0–4.6, far from both — so it is ~99 %
    H₂PO₄⁻ throughout and its charge barely moves. Under this module's INVERSE-ANCHORED cation
    a species of constant charge is absorbed by the anchor entirely, so phosphate would be not
    merely a weak buffer but very nearly a no-op. Pinned so that "add malt phosphate" cannot be
    re-proposed without confronting the number.
    """
    charge_lo = mean_charge(10.0**-4.0, PHOSPHATE)
    charge_hi = mean_charge(10.0**-4.6, PHOSPHATE)
    assert charge_lo == pytest.approx(0.9867, abs=5e-4)
    assert charge_hi == pytest.approx(0.9990, abs=5e-4)
    # Swing across the whole window: 0.01227 charge units per mole.
    phosphate_swing = charge_hi - charge_lo
    assert phosphate_swing == pytest.approx(0.01227, abs=5e-5)
    # Citric, over the SAME window, swings 0.3602 -- 29.4x more. The ratio is RECOMPUTED here
    # rather than transcribed: the first draft of this test asserted "> 30x" from a guess and
    # went red at 29.4 [[feedback-compute-the-clean-fix-before-adopting-it]].
    citric_swing = mean_charge(10.0**-4.6, CITRIC) - mean_charge(10.0**-4.0, CITRIC)
    assert citric_swing == pytest.approx(0.3602, abs=5e-4)
    assert citric_swing / phosphate_swing == pytest.approx(29.4, abs=0.2)


def test_no_pkas_raises():
    with pytest.raises(ValueError, match="at least one pKa"):
        _polyprotic_terms(1e-4, ())


# -- the two identities D-178 rests its "unchanged" claim on ------------------
#
# The record does NOT claim a measured bitwise diff (the pre-registered before/after
# trajectory comparison was not run). It claims two structural identities, which are
# exhaustive where a sampled comparison would not be. They are checked here rather than
# asserted in prose, so that the day either stops holding, the suite says so.


def test_the_general_branch_is_live_and_exactly_one_acid_reaches_it():
    """D-179 fired the trigger D-178 armed: beer's citrate is the first ≥3-pKa acid shipped.

    The D-178 form of this test asserted that NO shipped acid reached the general branch, and
    read :data:`ACID_STATE` to check it. It is rewritten rather than retired because it went
    **falsely green** at D-179: ``ACID_STATE`` became wine's registry alone, so the assertion
    kept passing while citrate — three pKas — went live in beer's. A guard that names one
    medium's registry cannot police a claim about *every* shipped acid, so this iterates
    :data:`ACID_REGISTRIES` and would notice a new acid in either medium.
    """
    from fermentation.core.acidbase import ACID_REGISTRIES, BYP_AS_SUCCINIC

    protons = {
        f"{medium}.{name}": spec.protons
        for medium, registry in ACID_REGISTRIES.items()
        for name, spec in registry.items()
    }
    protons["Byp"] = BYP_AS_SUCCINIC.protons
    assert protons == {
        "wine.tartaric": 2,
        "wine.malic": 2,
        "wine.lactic": 1,
        "beer.lactic": 1,
        "beer.acetic": 1,
        "beer.citrate": 3,
        "beer.malic": 2,
        # The D-181 falling acids. Oxalic is DIPROTIC, so citrate is still the only shipped
        # acid that reaches the n-protic branch — the claim this test's name makes survives
        # three new acids, which is the interesting part rather than an accident.
        "beer.pyruvic": 1,
        "beer.formic": 1,
        "beer.oxalic": 2,
        "beer.succinic": 2,
        "beer.peptide_buffer": 1,
        "Byp": 2,
    }
    polyprotic = {name for name, n in protons.items() if n >= 3}
    assert polyprotic == {"beer.citrate"}, (
        "the set of acids reaching the n-protic branch changed; D-178/D-179 justify that "
        "branch on citrate alone, so a new member needs its own sourcing and its own "
        "before/after comparison"
    )
    # And WINE still reaches none of it — the structural half of "wine is unchanged".
    assert all(protons[f"wine.{n}"] <= 2 for n in ("tartaric", "malic", "lactic"))


def test_the_registries_are_scoped_so_wine_does_not_gain_citrate():
    """The load-bearing separation: both media carry a ``citrate`` SLOT, one charge balance.

    Wine's ``citrate`` is a carbon-active MLF substrate deliberately kept out of the pH
    balance (D-31); beer's is charge-active and triprotic (D-179). A single union registry
    would silently make wine's charge-active — which is why the registries are per-medium and
    selected by an explicit medium label rather than by sniffing for a slot.

    This is the assertion that fails first if anyone "simplifies" the two registries back
    into one.
    """
    from fermentation.core.acidbase import acid_registry
    from fermentation.core.media import beer_schema, wine_schema

    wine, beer = wine_schema(), beer_schema()
    assert "citrate" in wine and "citrate" in beer, "both media carry the slot"
    assert "citrate" not in acid_registry(wine), (
        "wine's citrate has entered the charge balance — that is D-31's stated v1 omission "
        "being reversed, which is a decision with its own before/after comparison to run, "
        "not a side effect of scoping beer's acids"
    )
    assert "citrate" in acid_registry(beer)
    assert set(acid_registry(wine)) == {"tartaric", "malic", "lactic"}


def test_beer_gained_the_pH_system_and_wine_kept_its_slots():
    """Beer now carries ``cation_charge`` WITHOUT ``tartaric`` — the asymmetry D-178 predicted.

    The D-178 form of this test asserted beer had *neither* slot, and said in its own failure
    message that the day beer carried one without the other, the gate rename would stop being
    a no-op and the trajectories would have to be compared. That day is D-179; the comparison
    was run (see the record), and this is the same invariant restated for the new layout.

    What matters now is the *direction* of the asymmetry: ``EsterHydrolysis`` gates on
    ``tartaric`` (grape tartrate catalysis, which does not transfer), so it stays OFF for beer
    — while ``EthylAcetateEsterification`` gates on ``cation_charge`` and therefore switches
    ON. Reversing either gate silently changes which medium ages how.
    """
    from fermentation.core.media import beer_schema, wine_schema

    wine, beer = wine_schema(), beer_schema()
    assert ("tartaric" in wine, "cation_charge" in wine) == (True, True)
    assert ("tartaric" in beer, "cation_charge" in beer) == (False, True), (
        "beer must carry the cation WITHOUT grape tartaric: tartrate catalysis is grape "
        "chemistry (D-125) and must stay off for beer, while the pH system itself is now on"
    )
