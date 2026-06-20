"""Tests for the confidence-tier lattice."""

from fermentation.core.tiers import Tier, combine


def test_ordering_most_to_least_trustworthy():
    assert Tier.VALIDATED > Tier.PLAUSIBLE > Tier.SPECULATIVE


def test_labels():
    assert Tier.VALIDATED.label == "validated"
    assert Tier.PLAUSIBLE.label == "plausible"
    assert Tier.SPECULATIVE.label == "speculative"


def test_combine_picks_least_trustworthy():
    assert combine([Tier.VALIDATED, Tier.SPECULATIVE]) is Tier.SPECULATIVE
    assert combine([Tier.VALIDATED, Tier.PLAUSIBLE]) is Tier.PLAUSIBLE
    assert combine([Tier.VALIDATED, Tier.VALIDATED]) is Tier.VALIDATED


def test_combine_empty_is_validated_identity():
    # The identity for min over this lattice: combining "nothing" changes nothing.
    assert combine([]) is Tier.VALIDATED


def test_speculative_dominates_any_mix():
    mix = [Tier.PLAUSIBLE, Tier.VALIDATED, Tier.SPECULATIVE, Tier.PLAUSIBLE]
    assert combine(mix) is Tier.SPECULATIVE
