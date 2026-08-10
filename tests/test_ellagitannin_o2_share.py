"""How much O2 the oak-tannin sink takes is a JOINT, not a point — decision **D-175**.

``k_ellagitannin_oxidation``'s note claimed the banding made oak's O2 protection *partial*:
"roughly a third-to-half of the O2 ... **banded so it never monopolizes the O2**". Both halves were
computed with this constant *and* the always-on total held at their **nominals**, varying only the
ellagitannin pool — the archive's recurring verified-at-a-point / sampled-over-a-band shape
(D-118/D-154/D-155/D-157/D-170). Over the two bands the sampler can actually reach, the measured
share of the aging-phase O2 budget spans **4.58 %–88.40 %**, so the sink *can* monopolize the O2 and
the banding named as the reason it could not is what makes it possible.

**What is pinned here, and what deliberately is not.**

* The **static corner shares are RECOMPUTED from the shipped band edges**, never transcribed
  (D-154/D-158). Move any of the three edges and these reds — it forbids **silence**, not a
  direction: re-sourcing an edge is allowed, doing it without revisiting the note is not
  (``feedback-name-guards-for-what-they-forbid``).
* The **withdrawn claim may not return**, and the test names the phrases it forbids. This is a
  prose assertion on purpose and it is *not* the D-171 class of pinning-my-own-prose: the phrases
  are withdrawn as **measured false**, not as a matter of wording.
* The **cascade set has no referent for the claim** — there ``EllagitanninOxidation`` does not touch
  ``o2`` at all and the same constant weights the *quinone* branch. Structural, so it is pinned
  structurally rather than by integrating.
* **No edge is pinned to a value and none moved.** The source is ``author estimate``, which cannot
  license a narrowing (D-171, refused seven times). The span is stated, not closed.

The *measured* trajectory shares (4.58/48.62/88.40 %) are **not** asserted here: each costs a
five-year integration. What is asserted is that the cheap static bound brackets them in the right
direction, which is the supply-limited relationship the measurement established (D-136 — the O2
budget is fixed by ingress, so a faster sink takes a larger share of the *same* total and the
ceiling-based static estimate is a tight upper bound).
"""

from __future__ import annotations

import pytest

from fermentation.core.kinetics.aging import EllagitanninOxidation
from fermentation.core.kinetics.oxidative_cascade import QuinoneEllagitanninOxidation
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

#: The dose the note's own arithmetic is anchored to. ``oak_gpl`` is a **scenario input** with no
#: band at all, and ``oak_yield_ellagitannin_light`` is drawn but unreachable (D-159 — ``add_oak``
#: bakes ``oak_gpl · yield`` into the event at compile time), so the pool is not a sampled axis and
#: the joint below is over the two bands that genuinely reach a run.
_OAK_GPL = 4.0

#: The measured aging-phase shares the note prints, at this dose, direct set, ``T_ref`` (D-175).
_MEASURED_LOW, _MEASURED_NOMINAL, _MEASURED_HIGH = 0.0458, 0.4862, 0.8840


@pytest.fixture(scope="module")
def params():
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "oak.yaml",
        default_data_dir() / "aging.yaml",  # k_o2_depletion_total, the denominator
    )


def _share(k_ellag: float, yield_light: float, k_total: float) -> float:
    """Share of the O2 between this sink and the always-on total, at ``T_ref``.

    ``R = k · [ellagitannin] / k_total`` with both Arrhenius factors exactly 1 at ``T_ref``; the
    share of a **fixed** (ingress-limited) O2 budget is then ``R/(1+R)``.
    """
    r = k_ellag * (_OAK_GPL * yield_light) / k_total
    return r / (1.0 + r)


def test_the_reachable_joint_is_recomputed_from_the_shipped_edges(params):
    k = params["k_ellagitannin_oxidation"].uncertainty
    y = params["oak_yield_ellagitannin_light"]
    total = params["k_o2_depletion_total"].uncertainty

    low = _share(k.low, y.uncertainty.low, total.high)
    nominal = _share(
        params.value("k_ellagitannin_oxidation"), y.value, params.value("k_o2_depletion_total")
    )
    high = _share(k.high, y.uncertainty.high, total.low)

    # These three follow from the six edges above and the 4 g/L dose. The tolerance is deliberately
    # 1e-6 rather than the 5e-4 this test first shipped with: at 5e-4 the expected values for `low`
    # and `high` could be — and were — copied from the adjacent R column instead of the share
    # (R = 0.0160 vs share 0.015748; R/(1+R) = 0.958904 vs a transcribed 0.9586), because R ~ share
    # when R << 1 and R/(1+R) ~ 1 when R >> 1. A tolerance that admits the neighbouring quantity is
    # not a recompute pin. See D-175 amendment (a) — the guard's own defect class, caught by review.
    assert low == pytest.approx(0.015748, abs=1e-6)
    assert nominal == pytest.approx(0.500000, abs=1e-6)
    assert high == pytest.approx(0.958904, abs=1e-6)

    # The claim the withdrawn parenthetical denied: at the high corner this sink takes the great
    # majority of the O2. Asserted as the INEQUALITY, so it survives a re-sourced edge that keeps
    # the conclusion; the equalities above are what catch a silent move.
    assert high > 0.5, "the withdrawn 'never monopolizes' would require this to stay below a half"

    # The ENSEMBLE-reachable corners hold the pool at the dose's nominal ceiling, because the yield
    # band cannot reach a run (D-159). These are the corners the measured figures come from.
    ens_low = _share(k.low, y.value, total.high)
    ens_high = _share(k.high, y.value, total.low)
    assert ens_low == pytest.approx(0.047619, abs=1e-6)
    assert ens_high == pytest.approx(0.892857, abs=1e-6)

    # Supply limitation makes the static estimate a tight UPPER bound on the measured share, in the
    # same direction at every corner (D-136 + D-175's pool/ceiling median of 0.96-0.998).
    for measured, static in (
        (_MEASURED_LOW, ens_low),
        (_MEASURED_NOMINAL, nominal),
        (_MEASURED_HIGH, ens_high),
    ):
        assert measured <= static
        assert static - measured < 0.02


def test_the_withdrawn_monopolize_claim_has_not_returned(params):
    note = params["k_ellagitannin_oxidation"].uncertainty.note
    assert note is not None

    # What it forbids, named: the two withdrawn phrases may appear only under their withdrawal.
    assert "never monopolizes the O2" in note and "WITHDRAWN" in note
    assert "NOMINAL-ONLY" in note, "the third-to-half figure must stay scoped to the nominal arm"

    # And the three scopes the old text did not carry (D-175 §3).
    for scope in ("DIRECT SET ONLY", "T_ref ONLY", "THE POOL IS NOT A BAND"):
        assert scope in note, f"the note lost its {scope!r} scope"

    # No edge moved: the source cannot license a narrowing (D-171).
    assert "author estimate" in params["k_ellagitannin_oxidation"].provenance.source
    assert "NO EDGE MOVED" in note


def _oaked(oxidative: str):
    return compile_scenario(
        Scenario(
            name="d175-scope",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 200.0, "pitch_gpl": 0.25},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=40.0,
            closure="synthetic_supremecorq",
        ),
        oxidative=oxidative,
    )


def test_the_share_claim_has_no_referent_in_the_cascade_set():
    """The measured span is DIRECT-set only, and that is a structural fact, not a calibration."""
    direct = _oaked("direct")
    cascade = _oaked("cascade")

    # The sets are REPLACEMENTS (D-141), not re-wirings: the direct Process is not merely disabled
    # in the cascade set, it is ABSENT, and a differently-named Process reads the same constant.
    assert EllagitanninOxidation.name in direct.process_set
    assert EllagitanninOxidation.name not in cascade.process_set
    assert QuinoneEllagitanninOxidation.name in cascade.process_set
    assert QuinoneEllagitanninOxidation.name not in direct.process_set

    # Same constant, different oxidant: `o2` in one, `quinone` in the other. So "this sink takes N %
    # of the O2" has no referent in the cascade set, which is why the scope has to be stated in
    # prose rather than being visible from the parameter name.
    assert "o2" in EllagitanninOxidation.touches
    assert "o2" not in QuinoneEllagitanninOxidation.touches
    assert "quinone" in QuinoneEllagitanninOxidation.touches
    for reads in (EllagitanninOxidation.reads, QuinoneEllagitanninOxidation.reads):
        assert "k_ellagitannin_oxidation" in reads
