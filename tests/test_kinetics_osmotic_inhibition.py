"""Tests for OsmoticSubstrateInhibition — the high-sugar brake (decision D-192).

Two things make this term honest and both are pinned here rather than argued:

1. **It is EXACTLY inert below its threshold** — the factor is the literal ``1.0`` returned
   by an early branch, not a small number, so every must whose sugar never reaches
   ``S_osmotic_threshold`` is *byte-for-byte* the pre-D-192 core. That covers the whole of
   Coleman's 265-300 g/L validated envelope, which is why a speculative term can default on.
2. **It is asymptotic, never zero** — so a very sweet must ferments glacially rather than
   being frozen in an absorbing state it could never leave.

The tests also pin what the term deliberately does NOT do (the Handbook's 200-vs-300 g/L
alcohol statement is out of reach at the shipped threshold) so a later reader cannot mistake
its absence for an oversight.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.kinetics import (
    AminoAcidAssimilation,
    EthanolInactivation,
    GrowthNitrogenLimited,
    OsmoticSubstrateInhibition,
    SugarUptakeToEthanolCO2,
)
from fermentation.core.media import get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import abv_from_ethanol
from fermentation.validation import assert_conserved, total_carbon, total_nitrogen

_MODIFIER = "osmotic_substrate_inhibition"


@pytest.fixture(scope="module")
def params():
    return load_parameters(default_data_dir() / "wine_generic.yaml").resolve()


def _y(schema, s_total: float):
    return schema.pack(
        {"X": 1.0, "S": [s_total], "E": 10.0, "N": 0.1, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )


def _factor(params, s_total: float, **over) -> float:
    schema = wine_schema()
    p = dict(params) | over
    return OsmoticSubstrateInhibition().factor(0.0, _y(schema, s_total), schema, p)


def _run(brix: float, *, osmotic: bool, days: float = 60.0, yan: float = 200.0, n_eval=2001):
    sc = Scenario(
        name="osm",
        medium="wine",
        initial={"brix": brix, "yan_mgl": yan, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=days,
    )
    compiled = compile_scenario(sc, strict=True)
    if not osmotic:
        compiled.process_set.disable(_MODIFIER)
    t_eval = np.linspace(0.0, compiled.t_span_h[1], n_eval)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message
    return traj


def _sugar(traj) -> np.ndarray:
    s = np.asarray(traj.series("S"))
    return s if s.ndim == 1 else s.sum(axis=0)


# --------------------------------------------------------------------- metadata


def test_metadata():
    m = OsmoticSubstrateInhibition()
    assert m.name == _MODIFIER
    assert m.tier is Tier.SPECULATIVE
    # Growth rides along because the Handbook's mechanism is growth-mediated; the amino-acid
    # swap MUST ride with growth or it would refund carbon/nitrogen against an unmodified
    # draw (decision D-32, the same coupling the growth Arrhenius carries).
    assert set(m.modifies) == {
        SugarUptakeToEthanolCO2.name,
        GrowthNitrogenLimited.name,
        AminoAcidAssimilation.name,
    }
    assert set(m.reads) == {
        "S_osmotic_threshold",
        "K_osmotic_inhibition",
        "n_osmotic_inhibition",
    }


def test_the_modifier_is_wired_into_wine_and_not_into_beer(params):
    """Beer's inertness is structural — it is not wired, so no gravity can fire it.

    That is a stronger claim than "the threshold is above any wort": a parameter value
    could be changed, an absent modifier cannot fire.
    """
    wine = {m.name for m in get_medium("wine").build_process_set().active_modifiers}
    beer = {m.name for m in get_medium("beer").build_process_set().active_modifiers}
    assert _MODIFIER in wine
    assert _MODIFIER not in beer


# ------------------------------------------------------------------- the factor


@pytest.mark.parametrize("s_total", [0.0, 1.0, 100.0, 245.31, 291.0, 299.999, 300.0])
def test_the_factor_is_exactly_one_at_or_below_the_threshold(params, s_total):
    """Not ``approx(1.0)`` — the literal float, because inertness here is structural.

    The Process early-returns ``1.0`` rather than evaluating a form that happens to be near
    one, which is what lets an in-envelope run be byte-for-byte rather than merely close.
    """
    assert _factor(params, s_total) == 1.0


def test_a_negative_sugar_excursion_reads_as_an_exhausted_must_not_as_less_inhibition(params):
    assert _factor(params, -5.0) == 1.0


@pytest.mark.parametrize("s_total", [300.001, 350.0, 450.0, 625.0, 900.0, 5000.0])
def test_the_factor_is_strictly_between_zero_and_one_above_the_threshold(params, s_total):
    f = _factor(params, s_total)
    assert 0.0 < f < 1.0


def test_the_factor_is_monotonically_decreasing_in_sugar(params):
    xs = np.linspace(250.0, 1200.0, 200)
    fs = [_factor(params, float(x)) for x in xs]
    assert all(b <= a for a, b in zip(fs, fs[1:], strict=False))


def test_the_far_anchor_is_hit_by_the_shipped_pair(params):
    """f(625) = 0.05 is the constructed reading of "practically unfermentable".

    It is the ONLY thing that ties ``K_osmotic_inhibition`` to ``n_osmotic_inhibition``; the
    two are a derived pair, not independent constants, so this is the invariant that a
    change to either must preserve.
    """
    # Two pins, because the loose one alone does not discriminate. rel=5e-3 on f(625) is the
    # SHIPPED ROUNDING, not slack (K is rounded 74.56 -> 74.6, landing f(625) at 0.05005), but
    # it leaves a window in which a hand-edited K would pass. The sharp pin is the DERIVATION
    # itself: K is not an independent constant, it is a function of n.
    assert _factor(params, 625.0) == pytest.approx(0.05, rel=5e-3)
    n = params["n_osmotic_inhibition"]
    derived = (625.0 - params["S_osmotic_threshold"]) / ((1.0 / 0.05 - 1.0) ** (1.0 / n))
    assert params["K_osmotic_inhibition"] == pytest.approx(derived, rel=1e-3)


@pytest.mark.parametrize("n", [2.0, 3.0, 4.0, 6.0])
def test_the_anchor_holds_along_the_whole_admissible_curve(params, n):
    """Every ``n`` in the band passes through the far anchor *with its own* ``K``.

    This is what makes an independent per-parameter sweep wrong: pairing ``n = 6`` with
    ``n = 2``'s ``K`` misses the anchor badly, and that combination is not a member of the
    band even though both values are individually inside their own.
    """
    k = 325.0 / ((1.0 / 0.05 - 1.0) ** (1.0 / n))
    on_curve = _factor(params, 625.0, K_osmotic_inhibition=k, n_osmotic_inhibition=n)
    assert on_curve == pytest.approx(0.05, rel=1e-6)
    off_curve = _factor(params, 625.0, n_osmotic_inhibition=n)  # keeps the shipped K
    if n != 2.0:
        assert abs(off_curve - 0.05) > 0.02


def test_an_exponent_below_the_smoothness_floor_is_refused_not_integrated(params):
    """A corner at the threshold is exactly the discontinuity a BDF step pays for (D-182),
    so it is rejected at read time rather than left to the solver to discover."""
    with pytest.raises(ValueError, match="derivative corner"):
        _factor(params, 400.0, n_osmotic_inhibition=1.0)


def test_the_factor_is_c1_smooth_where_the_brake_engages(params):
    """One-sided derivatives must both be ~0 at the threshold, not merely finite."""
    thr = params["S_osmotic_threshold"]
    h = 1e-4
    below = (_factor(params, thr) - _factor(params, thr - h)) / h
    above = (_factor(params, thr + h) - _factor(params, thr)) / h
    assert below == pytest.approx(0.0, abs=1e-9)
    assert above == pytest.approx(0.0, abs=1e-6)


# ------------------------------------------------------- structural inertness


@pytest.mark.parametrize("brix", [20.0, 24.0, 26.0, 28.0])
def test_a_must_inside_colemans_envelope_is_byte_for_byte_the_pre_d192_core(brix):
    """The whole point of putting the threshold at 300 g/L rather than the Handbook's 200.

    24/26/28 Brix load at 245/268/291 g/L — the entire Coleman envelope — and must come out
    *identical*, not close. ``== 0.0`` exactly, like D-129's Gate 1.
    """
    on, off = _run(brix, osmotic=True), _run(brix, osmotic=False)
    assert np.abs(_sugar(on) - _sugar(off)).max() == 0.0
    assert np.abs(np.asarray(on.series("E")) - np.asarray(off.series("E"))).max() == 0.0
    assert np.abs(np.asarray(on.series("X")) - np.asarray(off.series("X"))).max() == 0.0


def test_the_threshold_sits_above_the_whole_coleman_envelope(params):
    """A guard on the *reason* the test above can pass, not just on the fact that it does.

    If a future edit moved the threshold below 300 g/L the inertness test would start
    failing for a reason no one could read off it; this names the constraint directly.
    """
    assert params["S_osmotic_threshold"] >= 300.0


# ------------------------------------------------------------- the consequence


def test_a_concentrated_must_no_longer_ferments_to_an_impossible_abv():
    """The defect D-192 exists to fix.

    An 881 g/L must is Tokaji-Eszencia concentration — past the 600-650 g/L the Handbook
    calls "practically unfermentable". Before this term it fermented to ~19.8 % ABV.
    """
    on, off = _run(70.0, osmotic=True), _run(70.0, osmotic=False)
    abv_on = abv_from_ethanol(float(np.asarray(on.series("E"))[-1]))
    abv_off = abv_from_ethanol(float(np.asarray(off.series("E"))[-1]))
    assert abv_off > 19.0, "guard premise changed: the pre-D-192 behaviour was ~19.8 % ABV"
    assert abv_on < 6.0
    assert _sugar(on)[-1] > _sugar(off)[-1] + 200.0


def test_the_concentrated_must_is_genuinely_stuck_not_merely_slow():
    """The evidence behind D-192's ``Corrects: D-129`` marker, which claims more than "slow".

    D-129 argued substrate inhibition "lifts as S drops, so the ferment finishes late, not
    stuck; it cannot leave residual sugar". Run to 20 years the 881 g/L arm is unchanged from
    ~1 year on, because ``EthanolInactivation`` drives viable biomass to zero long before the
    brake lifts. It is arrested, not en route — so the marker is correct as written.
    """
    five, twenty = _run(70.0, osmotic=True, days=1825.0), _run(70.0, osmotic=True, days=7300.0)
    assert _sugar(twenty)[-1] == pytest.approx(_sugar(five)[-1], rel=1e-3)
    assert _sugar(twenty)[-1] > 700.0
    assert float(np.asarray(twenty.series("X"))[-1]) < 1e-6, "biomass must be spent, not idle"


def test_a_very_sweet_must_still_ferments_it_is_not_an_absorbing_state():
    """Asymptotic, not a wall.

    A hard zero would leave uptake at exactly 0 with growth nitrogen-capped and nothing able
    to remove sugar — the must could never ferment at all, ever. Real must at this
    concentration does ferment; Tokay Aszu takes "from two to five years or longer".
    """
    traj = _run(70.0, osmotic=True, days=365.0)
    sugar = _sugar(traj)
    assert sugar[-1] < sugar[0] - 20.0, "the must is frozen, not merely slow"
    assert float(np.asarray(traj.series("E"))[-1]) > 5.0


def test_across_sweet_wine_the_brake_changes_the_path_but_not_the_destination():
    """Sauternes/TBA/icewine (~341-438 g/L) — and the distinction here is the whole point.

    An earlier version of this test asserted only that final ABV moves "well under a tenth
    of a %" and called the term *negligible* across this range. That was false as stated:
    these musts sit above the 300 g/L threshold, so the brake IS engaged, and mid-run sugar
    differs by 15 / 94 / 241 g/L at 32.2 / 36 / 40 °Brix. What is negligible is the
    ENDPOINT — a supply-limited flux slowed on the way down still arrives (this is exactly
    D-129's "late, not stuck", which holds *here* and fails only at the extreme regime the
    test above covers).

    So both halves are pinned: the trajectory must move, and the endpoint must not. A pin on
    the endpoint alone would pass just as well if the modifier were dead.
    """
    for brix, floor in ((32.2, 5.0), (36.0, 20.0), (40.0, 50.0)):
        on, off = _run(brix, osmotic=True), _run(brix, osmotic=False)
        path = float(np.abs(_sugar(on) - _sugar(off)).max())
        endpoint = abs(
            abv_from_ethanol(float(np.asarray(on.series("E"))[-1]))
            - abv_from_ethanol(float(np.asarray(off.series("E"))[-1]))
        )
        assert path > floor, f"{brix} Brix: brake not engaged mid-run (max|dS| {path:.2f})"
        assert endpoint < 0.1, f"{brix} Brix moved the endpoint {endpoint:.4f} %ABV"


def test_the_handbooks_200_vs_300_alcohol_statement_is_out_of_reach_and_that_is_deliberate():
    """The forfeit, pinned so it reads as a decision and not as an oversight.

    The Handbook says alcohol production "can be lower in a must containing 300 g/l than in
    another containing only 200 g/l". At the shipped threshold that ordering does NOT hold —
    more sugar still means more alcohol — because reaching it needs a brake that bites at
    200 g/L too, which is a global rate cut refuted by the keystone (D-192 measured the
    Coleman RMSE at 170.9 against a 2.0 g/L threshold for the version that does reproduce it).
    """
    lo = abv_from_ethanol(float(np.asarray(_run(21.55, osmotic=True).series("E"))[-1]))
    hi = abv_from_ethanol(float(np.asarray(_run(32.20, osmotic=True).series("E"))[-1]))
    assert hi > lo


# ------------------------------------------------------------------ invariants


def test_scaling_a_conserving_flux_leaves_carbon_and_nitrogen_closed(params):
    """A uniformly slower carbon-neutral flux is still carbon-neutral — measured, not assumed.

    Deliberately built on the ISOLATED core rather than the full wine medium: over a long
    run the full medium legitimately moves nitrogen between pools this quantity function does
    not count (autolysis, the amino-acid ledger), so a full-medium check would fail for
    reasons that have nothing to do with this modifier and would tell us nothing about it.
    The claim under test is narrow — scaling three conserving Processes by one scalar keeps
    their balances closed — so the set is narrow too.
    """
    schema = wine_schema()
    y0 = schema.pack(
        {"X": 0.1, "S": [700.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    # ``AminoAcidAssimilation`` is present only because ``ProcessSet`` refuses a modifier
    # naming an absent target — which is itself worth noting: that check is what makes the
    # wine-only wiring structural rather than conventional, since beer has no such Process
    # and so could not carry this modifier even by mistake. The pool is empty here, so it
    # contributes nothing to the balances being measured.
    pset = ProcessSet(
        schema,
        [
            GrowthNitrogenLimited(),
            SugarUptakeToEthanolCO2(),
            EthanolInactivation(),
            AminoAcidAssimilation(),
        ],
        modifiers=[OsmoticSubstrateInhibition()],
    )
    traj = simulate(pset, params=params, y0=y0, t_span=(0.0, 3000.0))
    assert traj.success, traj.message
    # Non-vacuous: the brake is engaged for the whole run (S never falls to the threshold).
    assert float(traj.series("S")[-1]) > params["S_osmotic_threshold"]
    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=params["biomass_C_fraction"]),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_conserved(
        traj,
        total_nitrogen(schema, biomass_nitrogen_fraction=params["biomass_N_fraction"]),
        rtol=1e-5,
        atol=1e-6,
        label="nitrogen",
    )


def test_the_output_tier_of_sugar_is_not_made_worse_than_it_already_was():
    """The honest tier consequence, the D-26/D-27 parallel: a speculative modifier scaling a
    plausible Process drags that Process's outputs down.

    ``Tier`` is ordered SPECULATIVE < PLAUSIBLE < VALIDATED, so "no worse" is ``>=``. The
    assertion is deliberately one-directional: sugar was ALREADY speculative here through
    the ethanol ceiling's reads (D-129), so this modifier costs nothing that was not already
    spent — and pinning "it drops the tier" would be pinning something untrue.
    """
    pset = get_medium("wine").build_process_set()
    with_it = pset.tier_of("S")
    pset.disable(_MODIFIER)
    without = pset.tier_of("S")
    assert with_it is Tier.SPECULATIVE
    assert without >= with_it
