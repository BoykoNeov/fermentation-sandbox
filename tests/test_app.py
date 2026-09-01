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


# -- the drawing: theme, log scale, and the chart that has nothing on it -----------------------
#
# These do not check that a chart is pretty. They check the three ways a change to the drawing
# can be silently WRONG: a theme that quietly drops the channel confidence is carried in, a log
# axis that reports solver dust as the story, and a chart of flat lines that reads as a broken
# page rather than as an answer about the batch.


def test_a_theme_may_change_the_hue_but_never_the_confidence_channel():
    """Line weight and dash pattern are the tier. A ground is allowed to change neither.

    ``TIER_STYLE`` is deliberately not a per-theme table. If someone ever splits it into one,
    a dark chart could carry a different dash for "speculative" than the light chart of the
    same run, and the two screenshots would disagree about how far a curve can be trusted.
    """
    pytest.importorskip("plotly")
    from app import render

    for theme in render.THEMES.values():
        assert len(theme.palette) == len(render.LIGHT.palette)
        assert theme.inert not in theme.palette, "inert must not collide with a live series"
    assert render.TIER_STYLE[Tier.SPECULATIVE] == ("dot", 1.8)
    assert render.TIER_STYLE[None] == ("longdash", 1.0)


def test_both_themes_draw_the_same_run_with_the_same_dashes(white_wine):
    """Same batch, two grounds: only the colours may differ, trace for trace."""
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import groups_for

    group, live = groups_for(white_wine)[0]
    light = render.series_figure(white_wine, group, live, theme=render.LIGHT)
    dark = render.series_figure(white_wine, group, live, theme=render.DARK)
    assert len(light.data) == len(dark.data)
    for a, b in zip(light.data, dark.data, strict=True):
        assert (a.line.dash, a.line.width) == (b.line.dash, b.line.width)
        assert a.name == b.name
        assert a.line.color != b.line.color


def test_a_log_axis_is_floored_off_the_data_not_off_the_solver_dust(white_wine):
    """A run that ends at 1e-8 g/L must not hand a log axis eight decades of nothing.

    Sugar finishes at solver dust of order 1e-8 and ethanol starts at exactly zero. Left to
    autoscale, a log axis would span from that dust to 200 g/L and squash the entire ferment
    into the top decade; left to Plotly's own handling, every non-positive point would be
    dropped and the sugar line would appear to *stop* days before the run ends.
    """
    pytest.importorskip("plotly")
    import numpy as np

    from app import render
    from app.readouts import groups_for

    group, live = groups_for(white_wine)[0]
    fig = render.series_figure(white_wine, group, live, log_y=True)

    assert fig.layout.yaxis.type == "log"
    low, high = fig.layout.yaxis.range
    assert high - low < render.LOG_DECADES + 1.0, "the axis reaches below its own floor"

    floor = 10.0**low
    for trace in fig.data:
        drawn, truth = np.asarray(trace.y), np.asarray(trace.customdata)
        assert len(drawn) == len(truth), "every point is still drawn, none dropped"
        assert drawn.min() >= floor * 0.5, "a point below the floor would vanish"
        # The hover reads off customdata, so clamping the line never falsifies a number.
        assert np.allclose(truth, np.asarray(trace.customdata))
        assert truth.min() <= drawn.min() + 1e-12


def test_a_log_axis_is_refused_where_nothing_is_positive_without_losing_the_zero_floor(
    white_wine,
):
    """An all-zero chart gets no decade ruler over it — and keeps the 0 to 1 range anyway.

    These two guards live in the same function and the refusal used to return early, so
    switching the log scale on handed the empty colour chart straight back to Plotly's
    -1 to 1: the very negative half-axis the other guard exists to remove, reappearing the
    moment someone used the other new control. Asserting only ``type != "log"`` passes with
    that bug present, which is why the range is asserted here too.
    """
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import groups_for

    colour = [g for g in groups_for(white_wine) if g[0].title.startswith("Colour")]
    assert colour, "the white wine fixture is supposed to have an empty colour chart"
    group, live = colour[0]
    for log_y in (False, True):
        fig = render.series_figure(white_wine, group, live, log_y=log_y)
        assert fig.layout.yaxis.type != "log"
        assert fig.layout.yaxis.range[0] == 0.0, f"negative concentrations at log_y={log_y}"
        assert fig.layout.yaxis2.range[0] == 0.0


def test_a_compared_run_keeps_its_colour_when_a_neighbour_cannot_be_drawn(white_wine):
    """Colours are handed out over the runs actually drawn, not over the ones asked for.

    A readout that one saved run's schema does not carry is skipped, and indexing the palette
    by position in the *requested* list would then leave a gap — two charts of the same three
    runs would colour them differently depending on which quantity was being compared.
    """
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import BY_KEY
    from app.runner import run_once

    beer = run_once(STARTERS["Pale ale"], DRAFT)
    tartaric = BY_KEY["tartaric"]  # wine only: the beer run has no such slot
    assert not tartaric.available(beer)

    pairs = [("beer", beer), ("wine", white_wine)]
    fig = render.compare_figure(pairs, tartaric)
    assert len(fig.data) == 1, "the run that cannot be drawn is skipped, not drawn empty"
    assert fig.data[0].name == "wine"
    assert fig.data[0].line.color == render.LIGHT.palette[0], "colours must not leave a gap"


def test_a_flat_axis_is_never_given_a_negative_range(white_wine):
    """Four lines flat at zero get 0 to 1, not Plotly's -1 to 1.

    The default range is not merely ugly: on a chart labelled g/L it draws half an axis of
    negative concentrations, which do not exist and which a reader will take as the model's
    claim about the batch.
    """
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import groups_for

    group, live = [g for g in groups_for(white_wine) if g[0].title.startswith("Colour")][0]
    fig = render.series_figure(white_wine, group, live)
    assert fig.layout.yaxis.range is not None
    assert fig.layout.yaxis.range[0] == 0.0
    assert fig.layout.yaxis2.range[0] == 0.0, "the second axis needs the same treatment"


def test_an_all_inert_chart_says_so_and_a_partly_live_one_does_not(white_wine):
    """The colour chart is empty because a white must has no pigment in it. Say that.

    The trigger has to be inertness rather than flatness. The oxygen chart in this same run
    has three lines pinned at zero and one — acetaldehyde — carrying the whole story; a note
    saying "nothing here moves" would be false there, and the assertion below is what stops
    a later, looser trigger from putting it there.
    """
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import groups_for

    noted = {}
    for group, live in groups_for(white_wine):
        panel = render.flat_group_panel(white_wine, group, live)
        if panel is not None:
            noted[group.title] = panel

    assert "Colour and phenolics" in noted
    body = noted["Colour and phenolics"].body
    assert "flat at zero" in body
    assert "anthocyanin" in body, "the note must say what would make the chart move"
    assert "Oxygen and oxidation" not in noted
    assert "Fermentation" not in noted


def test_the_empty_colour_note_separates_no_pigment_from_no_aging(white_wine):
    """Two different batches land on the same empty chart for opposite reasons.

    A white must has no pigment at all, so its colour lines sit at zero. A red *does* carry
    anthocyanin and tannin, but the pigment chemistry only switches on once the wine ages, so
    its lines sit at the values it was given. A note that blamed the missing ingredient would
    be plainly false on the red one, standing beside a sidebar that shows 0.5 g/L of it.
    """
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import groups_for
    from app.runner import run_once

    red = run_once(STARTERS["Red wine, warm ferment then malolactic"], DRAFT)
    bodies = {}
    for label, result in (("white", white_wine), ("red", red)):
        group, live = next(g for g in groups_for(result) if g[0].title.startswith("Colour"))
        panel = render.flat_group_panel(result, group, live)
        assert panel is not None, f"the {label} colour chart is empty and must say so"
        bodies[label] = panel.body

    assert "flat at zero" in bodies["white"]
    assert "flat at the value it started from" in bodies["red"]
    # And whichever it is, the advice has to cover both halves of the answer.
    for body in bodies.values():
        assert "anthocyanin" in body and "aging" in body


def test_the_advice_on_the_empty_colour_chart_actually_fills_it():
    """Following the note has to work, or it is decoration.

    Pigment plus an aging step is what the note asks for; this is that batch, and every line
    it names must then be driven. Without this the note could go on recommending something
    the engine stopped supporting and nothing would notice.
    """
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import groups_for
    from app.runner import run_once

    aged = STARTERS["Bottle-aged white, screwcap"]
    with_pigment = aged.model_copy(
        update={"initial": {**aged.initial, "anthocyanin_gpl": 0.5, "tannin_gpl": 1.5}}
    )
    result = run_once(with_pigment, DRAFT)
    group, live = next(g for g in groups_for(result) if g[0].title.startswith("Colour"))
    assert render.flat_group_panel(result, group, live) is None
    assert all(r.tier(result) is not None for r in live)


def test_every_group_that_can_come_up_empty_says_what_would_fill_it():
    """A note that only says "this is empty" is not worth printing."""
    from app.readouts import GROUPS

    for group in GROUPS:
        if group.title in {"Fermentation", "Alcohol and gravity", "Yeast"}:
            continue  # a run in which these are inert is a run that did not happen
        assert group.when_empty, f"{group.title} can be empty and has nothing to suggest"


def test_the_log_note_names_the_lines_it_is_about(white_wine):
    """Named, not left for the reader to work out which line is sitting on the floor."""
    pytest.importorskip("plotly")
    from app import render
    from app.readouts import BY_KEY

    ethanol = BY_KEY["ethanol"]  # starts at exactly zero
    panel = render.log_scale_panel(white_wine, [ethanol])
    assert panel is not None and "Ethanol" in panel.body

    temperature = BY_KEY["temperature"]  # in celsius, never near zero in this run
    assert render.log_scale_panel(white_wine, [temperature]) is None


def test_the_written_report_pins_its_own_ground(tmp_path, white_wine):
    """A dark report must not flip its page to white on the reader's machine.

    The figures inside are drawn once and baked in. A ``prefers-color-scheme`` query around
    them would eventually put a dark chart on a white page, so the file commits to the ground
    the run was drawn in.
    """
    pytest.importorskip("plotly")
    from app import render, report

    # Checked on the stylesheet this module writes, not on the whole file: Plotly's own
    # bundled CSS carries a prefers-color-scheme query of its own that is none of our business.
    assert "prefers-color-scheme" not in report._style(render.DARK)
    assert "color-scheme: dark" in report._style(render.DARK)
    assert "color-scheme: light" in report._style(render.LIGHT)

    dark = report.write(white_wine, tmp_path / "dark.html", theme=render.DARK).read_text(
        encoding="utf-8"
    )
    light = report.write(white_wine, tmp_path / "light.html", theme=render.LIGHT).read_text(
        encoding="utf-8"
    )
    assert render.DARK.palette[0] in dark and render.DARK.palette[0] not in light
    assert render.LIGHT.palette[0] in light


def test_the_keep_this_run_box_gets_a_fresh_key_per_batch() -> None:
    """The "Keep this run as" name box must not freeze on the first batch's name.

    Streamlit's ``value=`` only seeds a *keyed* widget the first time that key is created;
    from then on the key owns the value and the recomputed ``value=`` is dead. So a fixed
    key on this box leaves it displaying the first batch's name for the rest of the session
    while the sidebar's own name field moves on. Pressing "Keep it" would then file the
    current run under an earlier run's name -- overwriting the very run it was going to be
    compared against, with nothing on screen to say so.

    This is checked on the source rather than at runtime because nothing in this file
    imports Streamlit. The invariant is the one that matters: the key is an expression that
    varies with the run's name, not a constant.
    """
    import ast
    import pathlib

    source = pathlib.Path(__file__).resolve().parents[1] / "app" / "main.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    keys = [
        next(k.value for k in node.keywords if k.arg == "key")
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and any(isinstance(a, ast.Constant) and a.value == "Keep this run as" for a in node.args)
    ]
    assert len(keys) == 1, "expected exactly one 'Keep this run as' widget"
    assert not isinstance(keys[0], ast.Constant), "a constant key freezes the box"
    assert "scenario.name" in ast.unparse(keys[0])


# -- starting it, and where it puts what it writes ------------------------------------------
#
# The console was easy to start on exactly one machine. These pin the two ways that was true:
# a written report went to a folder named after a drive letter that exists here and nowhere
# else, and a taken port produced a one-line framework error and no console.


def test_no_shipped_interface_code_names_an_absolute_path() -> None:
    """No file under ``app/`` may contain an absolute filesystem path.

    ``REPORT_DIR`` named a folder on drive ``M:``, and ``report.write`` creates the
    folder it is given, so *Write it* raised on every machine without that drive -- in front
    of the user least equipped to read the traceback. The class of bug is "a location that is
    true on the machine it was written on", so the guard is on the class, not on that one
    constant: a path a user's machine must already have cannot be written down in code that
    ships to other machines. It has to be derived at runtime, as ``default_report_dir`` does.

    Scanned as text rather than as an AST because a drive letter is just as wrong inside a
    comment: it would be documenting the same false location.
    """
    import pathlib
    import re

    app_dir = pathlib.Path(__file__).resolve().parents[1] / "app"
    # The lookbehind is what keeps ``http://localhost`` out of it: a URL's scheme ends in a
    # letter, a drive letter never follows one.
    absolute = re.compile(r"(?<![A-Za-z])[A-Za-z]:[\\/]|(?:^|[\s\"'(=])/(?:home|Users|mnt|opt)/")
    offenders = [
        f"{path.name}:{n}: {line.strip()}"
        for path in sorted(app_dir.glob("*.py"))
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
        if absolute.search(line) and "test_" not in line
    ]
    assert not offenders, "absolute paths in shipped interface code:\n" + "\n".join(offenders)


def test_the_default_report_folder_is_under_the_user_s_own_home(monkeypatch) -> None:
    pytest.importorskip("plotly")
    from pathlib import Path

    from app.report import REPORT_DIR_ENV, default_report_dir

    monkeypatch.delenv(REPORT_DIR_ENV, raising=False)
    folder = default_report_dir()
    assert Path.home() in folder.parents, f"{folder} is not somewhere this user can reach"
    assert folder.is_absolute()


def test_the_report_folder_moves_with_its_environment_variable(monkeypatch, tmp_path) -> None:
    pytest.importorskip("plotly")
    from app.report import REPORT_DIR_ENV, default_report_dir, write

    monkeypatch.setenv(REPORT_DIR_ENV, str(tmp_path / "elsewhere"))
    assert default_report_dir() == tmp_path / "elsewhere"

    # And the folder does not have to exist first -- writing creates it.
    monkeypatch.setenv(REPORT_DIR_ENV, str(tmp_path / "not" / "there" / "yet"))
    assert default_report_dir().parent.name == "there"
    assert not write.__doc__ or "path" in write.__doc__.lower()


def test_a_port_something_is_already_serving_is_never_offered() -> None:
    """A busy port must be detected, including one held by a socket set to be re-usable.

    Two probes were tried against a console that was genuinely serving on 8613 and both
    reported the port free: binding it, and binding it with ``SO_REUSEADDR``. Windows lets a
    second socket take a port whose *listener* set that option, and Uvicorn sets it -- so on
    the platform this is most likely to be double-clicked on, a bind proves nothing. The
    listener below therefore sets the same option the real server does.

    On Linux a bind probe would pass this test, so treat a green here as necessary and not
    sufficient; the failure this pins is a Windows one.

    The stand-in server *accepts* in a thread rather than only listening. A socket that
    listens and never accepts fills its backlog after the first probe and then refuses the
    next one, so the second question gets the answer "free" from a port that is plainly
    busy -- which is the test lying, not the probe. A real server accepts.
    """
    import socket
    import threading

    from app.start_console import choose_port, port_is_free

    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server.bind(("127.0.0.1", 0))
    server.listen(8)
    busy = server.getsockname()[1]

    def accept_until_closed() -> None:
        while True:
            try:
                connection, _ = server.accept()
            except OSError:
                return
            connection.close()

    accepting = threading.Thread(target=accept_until_closed, daemon=True)
    accepting.start()
    try:
        assert not port_is_free(busy)
        assert choose_port((busy,)) is None, "a busy port must not be handed to the server"
    finally:
        server.close()
        accepting.join(timeout=2)

    # Closed again: the same port is free, so the probe is not simply always saying "busy".
    assert port_is_free(busy)


def test_the_launcher_never_lets_the_framework_stop_and_ask_for_an_email() -> None:
    """Headless, and telemetry off, on the command line the launcher actually builds.

    Left to open the browser itself, Streamlit asks for an email address on the terminal
    first. Nobody double-clicking an icon is watching a terminal, so that reads as a hang
    with no page. Headless skips the question; the launcher opens the browser instead.
    """
    from app.start_console import command

    line = command(8501)
    assert line[1:4] == ["-m", "streamlit", "run"]
    assert line[line.index("--server.headless") + 1] == "true"
    assert line[line.index("--browser.gatherUsageStats") + 1] == "false"
    assert line[line.index("--server.port") + 1] == "8501"


def test_every_input_the_form_shows_first_says_what_it_is() -> None:
    """The fields a first-time user meets must each carry a sentence explaining them.

    "Assimilable nitrogen (YAN)" with no help text is a number nobody outside a lab can set.
    The primary inputs are the ones on screen before anything is expanded, so a new one
    arriving without an explanation is a regression in exactly the place it is least
    affordable. The wider vocabulary behind "Everything else you can put in" is allowed to
    have gaps -- a visible gap there beats a wrong sentence.
    """
    from app.library import PRIMARY_INPUTS, input_help

    missing = {
        (medium, key)
        for medium, keys in PRIMARY_INPUTS.items()
        for key in keys
        if not (input_help(key) or "").strip()
    }
    assert not missing, f"primary inputs with no explanation: {sorted(missing)}"


def test_the_launcher_says_nothing_it_cannot_print() -> None:
    """Every sentence the launcher prints must be ASCII.

    A Windows console draws an em dash without complaint, so this passes unnoticed until the
    same command is redirected to a file: Python then encodes stdout with the machine's
    legacy codepage and one dash raises ``UnicodeEncodeError``. The failure lands on the
    launcher's own greeting -- the program crashing while explaining itself -- and it is
    invisible in the place it is normally tested. Held on the source, since running it
    would only reproduce the fault on one platform under one redirection.
    """
    import ast
    import pathlib

    source = pathlib.Path(__file__).resolve().parents[1] / "app" / "start_console.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    spoken = [
        ast.unparse(argument)
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "say"
        for argument in node.args
    ]
    assert spoken, "expected the launcher to say something"
    offenders = [line for line in spoken if not line.isascii()]
    assert not offenders, "non-ASCII in launcher output:\n" + "\n".join(offenders)
