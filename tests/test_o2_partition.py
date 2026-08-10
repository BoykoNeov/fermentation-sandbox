"""The always-on O2-depletion partition — decision **D-172**.

D-171 §6 flagged three parameters documented as one number and drawn as three:
``k_ethanol_oxidation`` + ``k_browning_base`` as a partition of a fixed 5.0e-4, and
``k_activation_floor`` as that same sum. D-172 replaced all three with
``k_o2_depletion_total`` and ``f_ethanol_o2_share``, so the identity holds by construction.

**What these tests are for, stated so they are not read as more than they are.** The sum
identity is now true *by construction* and asserting it proves nothing — that is exactly the
trap D-165 named ("median = nominal is an identity"). What is NOT free, and is what is
pinned here, is:

* that the partition is expressed as **two products** and not as a subtraction, which is an
  arithmetic requirement and the kind of thing a later tidy-up removes without noticing;
* that the three retired names stay retired, so the archive cannot silently regrow a fourth
  independent draw on the same number;
* that **both** oxidative sets read the same partition, since the pair the *cascade* drew
  (total against its own share) was the sharper of the two D-172 found;
* that the D-74 ordering claim is still **open** at the band edge — a guard that names the
  breach it does *not* forbid (``feedback-name-guards-for-what-they-forbid``).

Mutation-tested before it was written (``feedback-mutate-the-premise-before-building-the-
guard``): see D-172 §9 for which arms the shipped suite already caught and which it did not.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.kinetics.aging import OxidativeAcetaldehyde, PhenolicBrowning
from fermentation.core.kinetics.o2_partition import o2_depletion_shares
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario


def _compiled(medium: str):
    """A minimally-specified compiled scenario per medium — the real merge of parameter files."""
    if medium == "wine":
        # malic dosed so the named acid load can reach initial_ph without a negative cation.
        initial = {
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.25,
            "malic_gpl": 3.0,
            "initial_ph": 3.5,
        }
        celsius = 20.0
    else:
        initial = {
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 150.0,
            "pitch_gpl": 0.5,
        }
        celsius = 18.0
    return compile_scenario(
        Scenario(
            name=f"d172-{medium}",
            medium=medium,
            initial=initial,
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
            duration_days=10.0,
        )
    )


#: The names D-172 retired. Not a stylistic list — each was an independently sampled draw on
#: a quantity the other two also claimed to be, which is the defect this record repaired.
RETIRED = ("k_ethanol_oxidation", "k_browning_base", "k_activation_floor")


@pytest.fixture(scope="module")
def aging_params():
    d = default_data_dir()
    return load_parameters(d / "wine_generic.yaml", d / "aging.yaml", d / "acidbase.yaml")


def test_the_partition_reproduces_the_retired_nominals_bitwise(aging_params):
    # The whole repair rests on this being EXACT, not close: `k_browning_base` fed
    # `np.array_equal` isolability pins, so a one-ULP shift would red them for a reason with no
    # physical content. `==` on purpose, never pytest.approx — approx would pass on the very
    # error this test exists to forbid.
    params = aging_params.resolve()
    ethanol, browning = o2_depletion_shares(params)
    assert ethanol == 2.0e-4
    assert browning == 3.0e-4


def test_the_subtraction_form_is_why_this_is_written_as_two_products():
    # NOT a test of the code — a test of the ARITHMETIC FACT the code's shape depends on, so
    # that "simplify `total * (1 - f)` to `total - total * f`" fails here with the reason
    # attached instead of reddening four isolability pins somewhere else. The D-171 `_inverts`
    # precedent: make the argument un-tidyable by pinning it.
    total, f = 5.0e-4, 0.4
    assert total * (1.0 - f) == 3.0e-4
    assert total - total * f != 3.0e-4
    assert total - total * f == pytest.approx(3.0e-4)  # ... and approx cannot tell them apart


@pytest.mark.parametrize("medium", ["wine", "beer"])
def test_the_retired_names_are_gone_from_every_medium(medium):
    # Guards against the archive regrowing a fourth independent draw on the same number: if any
    # of these comes back as an entry, it is again a second opinion about `k_o2_depletion_total`.
    # Read through the COMPILE SEAM, not off a single YAML: a medium's ParameterSet is the merge
    # of several files, and a name reintroduced in any of them is a live draw.
    params = _compiled(medium).parameters
    present = [n for n in RETIRED if n in params.names]
    assert not present, f"{present} came back as YAML entries -- see D-172 before re-adding"
    assert "k_o2_depletion_total" in params.names
    assert "f_ethanol_o2_share" in params.names


def test_both_oxidative_sets_read_the_same_partition():
    # The pair the CASCADE drew was the sharper defect: the total against one of its own
    # shares, where the share's high edge was 8x the total's low. Both sets must now reach the
    # same two names, or that pairing can come back on one side only.
    direct = set(OxidativeAcetaldehyde.reads) | set(PhenolicBrowning.reads)
    assert {"k_o2_depletion_total", "f_ethanol_o2_share"} <= direct
    from fermentation.core.kinetics import oxidative_cascade as casc

    cascade_reads: set[str] = set()
    for obj in vars(casc).values():
        reads = getattr(obj, "reads", None)
        if isinstance(reads, tuple) and isinstance(obj, type):
            cascade_reads |= set(reads)
    assert "k_o2_depletion_total" in cascade_reads
    assert "f_ethanol_o2_share" in cascade_reads


def test_the_sum_is_the_total_across_the_band_interior_not_only_at_the_nominal(aging_params):
    # D-165: a distribution test at `x == mode` is vacuous, so this walks the band INTERIOR.
    # It is still a weak test by construction — the identity is mathematical — which is why the
    # assertion is on the RESIDUAL being at rounding scale, and why the docstring says the
    # falsifiable content lives in the ordering test below, not here.
    params = aging_params.resolve()
    p_total = aging_params["k_o2_depletion_total"].uncertainty
    p_f = aging_params["f_ethanol_o2_share"].uncertainty
    rng = np.random.default_rng(172)
    for total in np.linspace(p_total.low, p_total.high, 7):
        for f in np.linspace(p_f.low, p_f.high, 7):
            drawn = {**params, "k_o2_depletion_total": float(total), "f_ethanol_o2_share": float(f)}
            eth, brn = o2_depletion_shares(drawn)
            assert eth + brn == pytest.approx(total, rel=1e-15)
    for _ in range(50):  # off-grid draws too, so the grid cannot be the thing that passes
        total = float(rng.uniform(p_total.low, p_total.high))
        f = float(rng.uniform(p_f.low, p_f.high))
        eth, brn = o2_depletion_shares(
            {**params, "k_o2_depletion_total": total, "f_ethanol_o2_share": f}
        )
        assert eth + brn == pytest.approx(total, rel=1e-15)


def test_the_d74_ordering_is_still_open_at_the_split_fractions_high_edge(aging_params):
    # NAMES THE BREACH IT DOES NOT FORBID. D-74's load-bearing claim is that browning is the
    # DOMINANT share, i.e. f < 0.5. That holds at the nominal (0.4) and FAILS at the band's own
    # high edge (0.6), which its construction rule "~0.2-0.6x the total" licenses. D-172 did not
    # narrow that edge: both parents are `author estimate`, and narrowing a band so the author's
    # own prose comes true is the move D-171 refused seven times.
    #
    # So this test asserts the breach is REACHABLE, not that it is absent. It fails if someone
    # quietly moves the high edge to 0.5 — which would close the ordering, and which is a
    # decision that belongs in a record, not in a band edit.
    params = aging_params.resolve()
    high = aging_params["f_ethanol_o2_share"].uncertainty.high
    assert high > 0.5, "the D-74 ordering breach was closed by an edge move -- record it"

    at_high = {**params, "f_ethanol_o2_share": high}
    eth_hi, brn_hi = o2_depletion_shares(at_high)
    assert brn_hi < eth_hi  # the ordering INVERTS at the reachable edge

    eth_nom, brn_nom = o2_depletion_shares(params)
    assert brn_nom > eth_nom  # ... and holds at the nominal, which is all the guard covers

    # The breach rate is recomputed, never restated from the record (the D-154 rule): P(f >= 0.5)
    # for a triangular band is the closed-form upper tail.
    lo, mode = (
        aging_params["f_ethanol_o2_share"].uncertainty.low,
        aging_params["f_ethanol_o2_share"].value,
    )
    breach = (high - 0.5) ** 2 / ((high - lo) * (high - mode))
    assert breach == pytest.approx(0.125, abs=5e-4)
