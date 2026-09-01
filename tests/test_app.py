"""Guards for the interface layer (``app/``).

These are not "does the page render" tests — nothing here imports Streamlit. They pin the
handful of claims the interface makes that would be *silently* wrong if broken, which is the
only kind worth a test:

* An untouched quantity must report as inert, never as validated. The engine's own tier
  combine returns the *top* tier for an empty list, so the naive reading of ``tier_map``
  hands a UI the word "validated" for anything no mechanism writes. That is the single most
  misleading thing this interface could display, and
  :func:`test_untouched_reads_inert_not_validated` is the test that stops it coming back.
* A scenario must be compiled fresh for every run (decision D-206). A cache that stored a
  compiled scenario and ran it twice would produce a believable wrong second answer with no
  error at all.
* pH must not get an uncertainty band. It is solved from every acid at once, so a band drawn
  on any single input and relabelled would be a different quantity wearing pH's name.
* The solver presets must not offer an explicit method. Fermentation is always stiff.
"""

from __future__ import annotations

import inspect
import math

import pytest

from app.fidelity import METHODS, PRECISION_ORDER, PRECISION_PRESETS, Fidelity, tightened
from app.library import STARTERS, gate_problems
from app.readouts import BY_KEY, groups_for, headline_variables, summary
from app.runner import run_once, run_uncertainty, varying_constants
from fermentation.analysis import attribute_spread
from fermentation.core.tiers import Tier
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Intervention, Scenario, TemperaturePoint

DRAFT = Fidelity.preset("draft", points=100)


@pytest.fixture(scope="module")
def white_wine():
    """A short wine run with no aging step — so the aging quantities stay untouched."""
    return run_once(STARTERS["White wine, cool and clean"], DRAFT)


@pytest.fixture(scope="module")
def aged_wine():
    """A wine that begins aging partway through, i.e. one that switches chemistry on mid-run."""
    return run_once(STARTERS["Bottle-aged white, screwcap"], DRAFT)


# -- every shipped starter has to actually work --------------------------------------------


@pytest.mark.parametrize("name", list(STARTERS))
def test_every_starter_runs(name):
    """A starter that does not compile is a broken front door, not a broken example."""
    result = run_once(STARTERS[name], DRAFT)
    assert result.traj.success, result.traj.message
    assert result.traj.t[-1] > 0
    assert groups_for(result), "a finished run must have something to draw"
    assert summary(result), "a finished run must have headline numbers"


# -- the honesty guards --------------------------------------------------------------------


def test_untouched_reads_inert_not_validated(white_wine):
    """An untouched quantity must NOT inherit the top tier from an empty combine.

    The assertion is deliberately two-sided. It is not enough that the readout says "inert";
    the test also pins that the raw engine value it is protecting against really is
    ``VALIDATED``, so if the engine ever changes that convention this test fails loudly
    rather than passing for a reason that no longer exists.
    """
    browning = BY_KEY["a420"]
    assert browning.available(white_wine)
    assert "A420" not in white_wine.touched_variables(), (
        "this guard needs a run in which nothing writes A420; the white-wine starter has no "
        "aging step, so if this fails the starter changed"
    )
    # What the engine reports for it, unfiltered — the trap.
    assert white_wine.traj.tier_map["A420"] is Tier.VALIDATED
    # What the interface is required to report instead.
    assert browning.inert(white_wine) is True
    assert browning.tier(white_wine) is None


def test_speculative_is_falsy_so_none_is_the_only_inert_signal():
    """Why every tier check in ``app/`` is ``is None`` and never a truthiness test."""
    assert not Tier.SPECULATIVE, "if this ever becomes truthy the `is None` discipline can relax"
    assert Tier.SPECULATIVE is not None


def test_tier_is_lowest_of_the_sources_that_are_actually_driven(white_wine):
    """A readout's confidence is the worst of its live inputs, ignoring the inert ones."""
    ethanol = BY_KEY["ethanol"]
    tier = ethanol.tier(white_wine)
    assert tier is not None
    assert tier == white_wine.traj.tier_map["E"]


def test_no_chemistry_in_a_default_run_claims_validated(white_wine):
    """The claim printed on every page: nothing computed here is checked against real data."""
    touched = white_wine.touched_variables()
    validated = {
        n for n, t in white_wine.traj.tier_map.items() if n in touched and t is Tier.VALIDATED
    }
    assert validated <= {"T"}, (
        f"chemistry now claiming validated: {sorted(validated - {'T'})}. If that is real, the "
        "interface copy saying no chemistry has met the bar must change with it."
    )


def test_known_weak_readouts_carry_their_warning():
    """The two readouts the engine documents as misleading must never lose their warning."""
    for key in ("ta", "so2_bound"):
        assert BY_KEY[key].caveat, f"{key} lost the warning that has to travel with it"
    assert "starting value" in BY_KEY["ta"].caveat
    assert "under-count" in BY_KEY["so2_bound"].caveat


# -- decision D-206: never reuse a compiled scenario ----------------------------------------


def test_each_run_starts_from_a_clean_scenario(aged_wine):
    """Running the same batch twice must give the same answer.

    The failure this guards against is not an exception. Running a compiled scenario leaves
    its switches flipped, so a cache that handed back the same compiled object would start
    the second run with the aging chemistry already on from t = 0 — a different, entirely
    believable curve. The aged starter is used because it is the one with a mid-run switch.
    """
    again = run_once(aged_wine.scenario, DRAFT)

    switched_first = {m.name for m in aged_wine.mechanisms if m.switched_on_mid_run}
    switched_again = {m.name for m in again.mechanisms if m.switched_on_mid_run}
    assert switched_first, "this guard needs a scenario that switches chemistry on mid-run"
    assert switched_first == switched_again, (
        "the second run began with the first run's chemistry already enabled — a compiled "
        "scenario is being reused (decision D-206)"
    )
    assert again.traj.y[:, -1] == pytest.approx(aged_wine.traj.y[:, -1], rel=1e-12, abs=1e-12)


# -- solver settings ------------------------------------------------------------------------


def test_standard_preset_is_the_engine_default():
    """The default setting must be byte-for-byte what calling the library plainly would do."""
    engine = inspect.signature(simulate).parameters
    assert PRECISION_PRESETS["standard"] == (
        engine["rtol"].default,
        engine["atol"].default,
    )
    assert Fidelity().method == engine["method"].default
    assert Fidelity().max_step == math.inf


def test_presets_never_choose_a_non_stiff_solver():
    """Fermentation always mixes fast and slow chemistry; RK45 belongs behind 'custom' only."""
    assert "RK45" in METHODS, "still offered as a diagnostic"
    for name in PRECISION_ORDER:
        assert Fidelity.preset(name).method == "BDF"


def test_tightening_walks_up_the_presets_and_stops():
    step = Fidelity.preset("draft")
    for expected in ("standard", "high"):
        nxt = tightened(step)
        assert nxt is not None and nxt.precision == expected
        assert nxt.rtol < step.rtol
        step = nxt
    assert tightened(step) is None, "nothing is tighter than the tightest preset"


def test_tightening_a_custom_setting_keeps_everything_but_the_tolerances():
    custom = Fidelity(precision="custom", method="Radau", rtol=1e-5, atol=1e-8, oxidative="cascade")
    tighter = tightened(custom)
    assert tighter is not None
    assert tighter.rtol == pytest.approx(1e-7)
    assert tighter.atol == pytest.approx(1e-10)
    # A convergence check must compare like with like: only the tolerances may move.
    assert tighter.method == "Radau"
    assert tighter.oxidative == "cascade"


def test_stored_points_reach_the_trajectory():
    """The output-resolution axis has to actually do the one thing it claims to do.

    The stored count can exceed the request by the number of scheduled breakpoints: a
    temperature knot or a dosing day is always added to the grid so the jump lands on a real
    stored point rather than being smeared between two. So this pins the floor and the
    ordering, not an exact count.
    """
    sparse = run_once(STARTERS["Pale ale"], Fidelity.preset("draft", points=100))
    denser = run_once(STARTERS["Pale ale"], Fidelity.preset("draft", points=500))
    breakpoints = len(sparse.traj.segment_bounds) - 2  # the two ends are not extra points
    assert 100 <= len(sparse.traj.t) <= 100 + breakpoints
    assert 500 <= len(denser.traj.t) <= 500 + breakpoints
    assert len(denser.traj.t) > len(sparse.traj.t)


def test_headline_variables_exist_in_the_schema(white_wine):
    for name in headline_variables(white_wine):
        assert name in white_wine.schema


# -- what may and may not carry an uncertainty band ------------------------------------------


def test_ph_refuses_a_band_and_ethanol_accepts_one(white_wine):
    """pH is not one tracked quantity scaled, so a band on it would be a different quantity."""
    plotly = pytest.importorskip("plotly")  # noqa: F841  (the ui dependency group)
    from app.render import bandable

    assert bandable(BY_KEY["ethanol"], white_wine) is True
    assert bandable(BY_KEY["abv"], white_wine) is True, "alcohol by volume is a plain scaling"
    assert bandable(BY_KEY["ph"], white_wine) is False, "pH comes from every acid at once"

    # Wine has one sugar, so its "total" is that one quantity and a band is honest. Beer has
    # three, so the total is a sum and a band drawn on the sum would not be one.
    assert bandable(BY_KEY["sugar"], white_wine) is True
    pale_ale = run_once(STARTERS["Pale ale"], DRAFT)
    assert bandable(BY_KEY["sugar"], pale_ale) is False


def test_the_projected_scope_is_the_engine_scope_minus_what_cannot_move():
    """The count the page prints has to be the engine's own answer, twice corrected.

    Re-deriving the scope by walking a finished run's mechanisms undercounts — it misses
    what the schedule reads and what the set-up step read to build the starting state
    (D-241). Measured before this was fixed: 12 short of 97 on a wine, 8 short of 89 on a
    beer, always low. And the engine's raw scope over-counts in the other direction, because
    it happily draws a number pinned to a single value, which has no variance and is dropped
    from the ranking.

    This runs a tiny ensemble purely to read back what it really sampled.
    """
    fidelity = Fidelity.preset("draft", points=40)
    for name in ("White wine, cool and clean", "Pale ale"):
        projected = set(varying_constants(STARTERS[name], fidelity))
        ens = run_uncertainty(STARTERS[name], fidelity, n_members=4, seed=0, sampler="mc")
        drawn = set(ens.ensemble.sampled_names)
        can_move = {
            n
            for n in drawn
            if ens.parameters[n].uncertainty.high > ens.parameters[n].uncertainty.low
        }
        assert projected <= drawn, f"{name}: projecting names the engine will not draw"
        assert projected == can_move, (
            f"{name}: the projection is out by {sorted(can_move ^ projected)} — the page "
            "would print a re-run count the ensemble does not agree with"
        )


def test_the_printed_re_run_count_is_the_one_the_ranking_actually_needs():
    """The promise, end to end, at a size the suite can afford.

    The page says the ranking needs more re-runs than there are varying numbers. That is only
    worth printing if it is exactly true at the boundary, so this pins both sides of it: one
    re-run per varying number is refused, one more than that succeeds. Run against a narrowed
    scope so the whole test is a handful of seconds rather than a minute.
    """
    fidelity = Fidelity.preset("draft", points=40)
    scenario = STARTERS["Pale ale"]

    # A deliberately small slice of the numbers, so the boundary is cheap to reach.
    candidates = varying_constants(scenario, fidelity)[:3]
    scope = varying_constants(scenario, fidelity, only=list(candidates))
    assert len(scope) == 3, "the narrowing did not survive into the engine's own resolution"

    def rank(n_members: int):
        ens = run_uncertainty(
            scenario, fidelity, n_members=n_members, seed=0, sampler="mc", only=list(scope)
        )
        return attribute_spread(ens.ensemble, "E", ens.param_tiers)

    with pytest.raises(ValueError, match="underdetermined"):
        rank(len(scope))
    attribution = rank(len(scope) + 1)
    assert set(attribution.per_param) <= set(scope)


def test_narrowing_the_scope_gives_a_smaller_set():
    """Focusing on one quantity must actually reduce what gets drawn."""
    fidelity = Fidelity.preset("draft", points=40)
    scenario = STARTERS["Pale ale"]
    everything = varying_constants(scenario, fidelity)
    result = run_once(scenario, fidelity)
    reads = sorted({n for m in result.touching("E") for n in m.reads})
    narrowed = varying_constants(scenario, fidelity, only=reads)
    assert 0 < len(narrowed) < len(everything)
    assert set(narrowed) <= set(reads)


def test_nothing_projected_is_pinned_to_a_single_value():
    """A number drawn at its own value every time explains nothing and must not be counted."""
    fidelity = Fidelity.preset("draft", points=40)
    scenario = STARTERS["White wine, cool and clean"]
    result = run_once(scenario, fidelity)
    names = varying_constants(scenario, fidelity)
    assert names
    pinned = [
        n
        for n in names
        if result.parameters[n].uncertainty.high <= result.parameters[n].uncertainty.low
    ]
    assert not pinned, f"drawn at a fixed value and so explaining nothing: {pinned}"


# -- the form's cross-cutting gates -----------------------------------------------------------


def _wine(**kw) -> Scenario:
    base: dict[str, object] = {
        "name": "gate",
        "medium": "wine",
        "initial": {"brix": 22.0, "yan_mgl": 200.0, "pitch_gpl": 0.5},
        "temperature_schedule": [TemperaturePoint(day=0.0, celsius=20.0)],
        "duration_days": 20.0,
    }
    base.update(kw)
    return Scenario(**base)


def test_setting_a_target_ph_needs_a_starting_ph():
    sc = _wine(interventions=[Intervention(day=2.0, action="set_ph", params={"ph": 3.4})])
    assert any("starting pH" in p for p in gate_problems(sc))


def test_sealing_a_bottle_needs_a_closure_and_an_aging_step():
    no_closure = _wine(
        interventions=[
            Intervention(day=5.0, action="begin_aging"),
            Intervention(day=5.0, action="seal_bottle"),
        ]
    )
    assert any("closure" in p for p in gate_problems(no_closure))

    before_aging = _wine(
        closure="screwcap",
        interventions=[
            Intervention(day=2.0, action="seal_bottle"),
            Intervention(day=5.0, action="begin_aging"),
        ],
    )
    assert any("before" in p for p in gate_problems(before_aging))


def test_an_intervention_after_the_run_ends_is_flagged():
    sc = _wine(
        duration_days=5.0,
        interventions=[Intervention(day=9.0, action="add_so2", params={"so2_mgl": 30.0})],
    )
    assert any("never fire" in p for p in gate_problems(sc))


def test_a_valid_scenario_has_no_gate_problems():
    assert gate_problems(STARTERS["Red wine, warm ferment then malolactic"]) == []


# -- the written report ------------------------------------------------------------------------


def test_report_is_one_self_contained_file(tmp_path, white_wine):
    pytest.importorskip("plotly")
    from app.report import write

    out = write(white_wine, tmp_path / "run.html", explain=("E",))
    text = out.read_text(encoding="utf-8")
    assert text.startswith("<!doctype html>")
    assert "plotly" in text.lower(), "the charts' javascript must be inlined, not linked"
    assert '<script src="http' not in text, "nothing may be fetched from the network"
    # The warnings have to survive the trip into the file.
    assert "Trust the starting value" in text or "under-count" in text
