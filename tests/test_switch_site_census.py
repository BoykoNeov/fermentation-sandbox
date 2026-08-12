"""Mechanism B — where a parameter's band straddles a code switch (decision D-166).

D-165 §4 found that ``k_d2_ethanol_tolerance_death``, the archive's worst-displaced band,
is **exactly inert** on a standard 24-Brix must, because ``EthanolToleranceDeath`` gates it
behind ``max(E - ethanol_tolerance, 0)**2`` and the ferment never reaches the tolerance. It
recorded the margin as **0.5 %** in ethanol concentration and left the general case open as
"Mechanism B, the switch-site census".

This file pins what that census found. Three things are worth stating up front.

**The enumeration is over parameters, not code sites.** D-165's own "Next" line named the
code's ``min``/``max``/``clip`` comparisons as the surface. That is wrong as an *enumerator*:
a ``params[`` -operand grep over ``src/`` returns four sites and does **not** contain
``ethanol_tolerance``, the instance the whole finding was built around, because it routes
through a local (``over = e - params["ethanol_tolerance"]``). The census therefore sweeps
every live band and pins each parameter at its band edges; code sites only *attribute* a
straddle once it is found.

**The classification is exact, not tolerance-based.** With a fixed ``t_eval`` grid, two runs
with equal parameter maps are bitwise identical, so "inert" means ``np.array_equal`` and
needs no epsilon:

    DEAD      bitwise identical to nominal at BOTH band edges
    ACTIVE    moves at both
    STRADDLE  moves at one edge and not the other -- the band decides whether the
              mechanism exists at all, which is a *qualitative* uncertainty, not a
              magnitude one

**The margin is much tighter in the coordinate a winemaker actually sets.** D-165's 0.5 % is
in ethanol. Re-expressed in Brix, ``ethanol_tolerance`` straddles over
**[24.56, 28.45] deg Brix** — and the reference must at 24.0 sits **0.56 Brix below** that
window. Picking half a Brix riper, well inside normal harvest variation and inside a
refractometer's read error, moves a speculative death term from "contributes nothing" to
"its band decides whether it fires". That is the guard this file exists for.

**What is deliberately NOT pinned here: a fix.** No band is edited, no gate is changed. D-166
measures and records, exactly as D-165 did with the log-triangular candidate.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario.compile import compile_scenario
from fermentation.scenario.schema import Scenario, TemperaturePoint

#: The pH-solver inputs ``_ALLOWED_KEYS`` documents (D-18). Run 1 of the D-166 census omitted
#: these and integrated a must with every acid pool at zero -- pH 7.0, i.e. water -- which
#: made every DEAD verdict unattributable between "a threshold gates it" and "the pool is
#: empty". Tests here use a populated must so that distinction is real.
ACID_MUST = {"tartaric_gpl": 6.0, "malic_gpl": 3.0, "initial_ph": 3.4, "so2_total_mgl": 50.0}

#: Pinned from ``wine_generic.yaml`` by ``test_the_reference_bands_are_the_shipped_ones``.
_TOL_LO, _TOL_VAL, _TOL_HI = 120.0, 142.0, 150.0

#: The straddle window measured by bisection to 0.1 Brix (D-166 §3).
_WINDOW_OPENS, _WINDOW_CLOSES = 24.56, 28.45
_REFERENCE_BRIX = 24.0


def _wine(brix: float = _REFERENCE_BRIX, days: float = 30.0, **knobs: float) -> Scenario:
    initial: dict[str, float] = {"brix": brix, "yan_mgl": 250.0, "pitch_gpl": 0.25, **ACID_MUST}
    initial.update(knobs)
    return Scenario(
        name=f"d166-{brix:g}",
        medium="wine",
        initial=initial,
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=14.0, celsius=25.0),
        ],
        interventions=[],
        duration_days=days,
    )


def _classify(scenario: Scenario, name: str, n_grid: int = 121) -> tuple[str, float]:
    """(class, peak ethanol) for one parameter: pin it at each band edge, compare bitwise.

    Everything else stays at nominal, so there is no sampling and therefore no
    draw-sequence confound to control for (contrast D-165 §0's ``exclude=``/``only=`` bug).
    """
    c = compile_scenario(scenario)
    events = tuple(c.events or ())
    grid = np.linspace(c.t_span_h[0], c.t_span_h[1], n_grid)
    base = c.parameters.resolve()

    def run(override: dict[str, float] | None = None):
        p = dict(base)
        if override:
            p.update(override)
        return simulate_scheduled(c.process_set, p, c.y0, c.t_span_h, events=events, t_eval=grid)

    nominal = run()
    assert nominal.success, nominal.message
    band = c.parameters[name].uncertainty
    moved_lo = not np.array_equal(run({name: band.low}).y, nominal.y)
    moved_hi = not np.array_equal(run({name: band.high}).y, nominal.y)
    if moved_lo and moved_hi:
        cls = "ACTIVE"
    elif moved_lo or moved_hi:
        cls = "STRADDLE"
    else:
        cls = "DEAD"
    return cls, float(np.max(nominal.series("E")))


def test_the_reference_bands_are_the_shipped_ones():
    """Non-vacuity. Every constant above is read back out of the real parameter store, so
    these tests cannot derive their expectations from themselves -- the D-108/D-109 shape
    that D-164 §5 and D-165 had each to name in their own test files."""
    c = compile_scenario(_wine())
    band = c.parameters["ethanol_tolerance"].uncertainty
    assert (band.low, c.parameters.value("ethanol_tolerance"), band.high) == (
        _TOL_LO,
        _TOL_VAL,
        _TOL_HI,
    )


def test_two_runs_with_equal_parameter_maps_are_bitwise_identical():
    """The premise the whole exact criterion rests on. If this ever fails, every DEAD
    verdict in this file degrades to "below some unstated noise floor" and the census
    would need a tolerance -- which is exactly what D-165's `approx`-on-either-bound
    lesson says not to reach for silently."""
    c = compile_scenario(_wine())
    events = tuple(c.events or ())
    grid = np.linspace(c.t_span_h[0], c.t_span_h[1], 121)
    base = c.parameters.resolve()
    a = simulate_scheduled(c.process_set, dict(base), c.y0, c.t_span_h, events=events, t_eval=grid)
    b = simulate_scheduled(c.process_set, dict(base), c.y0, c.t_span_h, events=events, t_eval=grid)
    assert np.array_equal(a.y, b.y)


def test_ethanol_tolerance_is_dead_at_the_reference_must():
    """D-165 §4 replicated by a different harness. Its measurement was an ensemble with
    ``only=[k_d2]``; this pins each band edge deterministically instead, on a must that
    actually contains acid. Same verdict."""
    cls, peak = _classify(_wine(), "ethanol_tolerance")
    assert cls == "DEAD", f"expected DEAD at {_REFERENCE_BRIX} Brix, got {cls} (peak E {peak:.4g})"
    assert peak < _TOL_LO


def test_the_curvature_constant_is_dead_wherever_its_gate_is_shut():
    """``k_d2_ethanol_tolerance_death`` is the archive's worst-displaced band (r = 300,
    sampled median 9.18x nominal, D-165 §2) and ``_active_reads`` says every wine ensemble
    draws it. It still contributes exactly nothing here, because it multiplies a term that
    is identically zero below tolerance. This is why a ``reads``-based reach count is an
    UPPER BOUND on consequential reach."""
    cls, _ = _classify(_wine(), "k_d2_ethanol_tolerance_death")
    assert cls == "DEAD"


def test_the_reference_must_sits_just_below_the_straddle_window():
    """**The guard.** D-165 stated the margin as 0.5 % of ethanol concentration. In Brix --
    the quantity a winemaker actually sets, and the one a refractometer reports -- the
    margin is **0.56 deg Brix**: a must picked at 24.56 instead of 24.0 puts
    ``ethanol_tolerance`` inside its straddle window, where the *band* decides whether a
    speculative death term fires at all.

    This fails if the window's lower edge drifts down onto the reference must. It is a
    guard on the gap, not on the gate (D-165's phrasing, and the same intent)."""
    assert _WINDOW_OPENS > _REFERENCE_BRIX
    cls_ref, peak_ref = _classify(_wine(_REFERENCE_BRIX), "ethanol_tolerance")
    cls_in, peak_in = _classify(_wine(26.0), "ethanol_tolerance")
    assert cls_ref == "DEAD", (
        f"the reference {_REFERENCE_BRIX}-Brix must has entered the straddle window "
        f"(peak E {peak_ref:.4g} vs band low {_TOL_LO}) -- ethanol_tolerance's band now "
        "decides whether EthanolToleranceDeath fires on a NORMAL wine. Re-measure D-165 §4 "
        "and D-166 §3 before relaxing this."
    )
    assert cls_in == "STRADDLE", (
        f"26 Brix was measured as inside the straddle window [{_WINDOW_OPENS}, "
        f"{_WINDOW_CLOSES}] Brix; got {cls_in} (peak E {peak_in:.4g})"
    )
    assert _TOL_LO < peak_in < _TOL_HI


def test_the_inertness_margin_is_held_open_by_the_carbon_diversion_bands():
    """**What actually keeps the gate shut at 24 Brix** -- and it is not the sugar load.

    The must holds 245.3 g/L sugar, so the Gay-Lussac ceiling is ``0.511 * 245.3 = 125.4``
    g/L ethanol, which is *above* ``ethanol_tolerance``'s band low of 120. The only reason a
    24-Brix ferment peaks at 117.3 instead is the D-16 realised-yield carbon diversion:
    ``Y_glycerol_sugar`` (0.035) and ``Y_byproduct_sugar`` (0.012) carve carbon out of the
    ethanol flux. Both are sampled in every wine ensemble.

    Measured deterministically at 24 Brix:

    ========================================  ==========  ==============
    diversion                                 peak E      margin to 120
    ========================================  ==========  ==============
    shipped nominal (0.035, 0.012)            117.280     +2.27 %
    both yield bands at their LOW edge        119.696     **+0.25 %**
    off entirely (0, 0) -- the Gay-Lussac core 122.959    **-2.47 %**
    ========================================  ==========  ==============

    So D-165 §4's "0.5 % of headroom" across a 64-member ensemble is not a coincidence of
    the draw: it is these two bands near their low ends. And with the diversion off -- which
    ``uptake.py`` documents as the togglable pure Gay-Lussac core -- the reference must is
    **already inside the straddle window**, with no change in Brix at all.

    This test pins the ordering, so a widened yield band or a raised diversion default
    cannot close the 0.30 g/L gap silently."""
    c = compile_scenario(_wine())
    events = tuple(c.events or ())
    grid = np.linspace(c.t_span_h[0], c.t_span_h[1], 121)
    base = c.parameters.resolve()

    def peak(overrides: dict[str, float]) -> float:
        p = dict(base)
        p.update(overrides)
        r = simulate_scheduled(c.process_set, p, c.y0, c.t_span_h, events=events, t_eval=grid)
        assert r.success, r.message
        return float(np.max(r.series("E")))

    gly, byp = c.parameters["Y_glycerol_sugar"], c.parameters["Y_byproduct_sugar"]
    assert gly.uncertainty.low == 0.02 and byp.uncertainty.low == 0.007

    peak_nominal = peak({})
    peak_worst = peak(
        {"Y_glycerol_sugar": gly.uncertainty.low, "Y_byproduct_sugar": byp.uncertainty.low}
    )
    peak_off = peak({"Y_glycerol_sugar": 0.0, "Y_byproduct_sugar": 0.0})

    assert peak_nominal < peak_worst < _TOL_LO < peak_off, (
        f"peak ethanol ordering changed: nominal {peak_nominal:.4g}, both yield bands low "
        f"{peak_worst:.4g}, diversion off {peak_off:.4g}, band low {_TOL_LO}"
    )
    # The gap that D-165 recorded as 0.5% across an ensemble is 0.25% deterministically here.
    assert (_TOL_LO - peak_worst) / _TOL_LO < 0.005, (
        f"the worst-case diversion draw now leaves only "
        f"{100 * (_TOL_LO - peak_worst) / _TOL_LO:.3f}% of headroom to ethanol_tolerance's "
        "band low -- re-measure D-165 §4 and D-166 before widening either yield band."
    )
    # And the Gay-Lussac ceiling really does sit above the gate, which is the reason this
    # margin is a diversion property rather than a sugar-load one.
    sugar0 = float(np.sum(c.y0[c.schema.slice("S")]))
    assert 0.511 * sugar0 > _TOL_LO


def test_far_past_the_window_the_gate_is_open_at_both_edges():
    """The window **closes**, and that is not obvious. A straddle needs the peak ethanol
    *inside* [120, 150]; by 36 Brix the must overshoots the band high, so both edges fire
    and the band no longer decides existence -- only magnitude. (Pre-registered as D-166 P6
    against a reviewer's prediction that 36 Brix would straddle; the harness adjudicated.)

    The window is narrow -- 3.89 Brix -- precisely because the mechanism it gates is
    self-limiting: the must sticks, so peak ethanol saturates near 156 g/L rather than
    climbing to the ~187 g/L the ``inactivation.py`` docstring quotes for the *uncapped*
    extrapolation.

    ``ACTIVE`` alone does **not** establish that both edges fire, and the D-166 mutation
    campaign is what showed it: ``_classify`` compares each edge against the *nominal* run,
    and past the window the nominal already has the gate open, so every edge differs from it
    -- including an edge the peak never reaches. Widening the band high to 200 g/L left this
    test green while the high-edge run had the gate physically **shut**. The assertion below
    therefore reads the high edge out of the store and requires the peak to clear it, which
    is the claim the docstring actually makes."""
    scenario = _wine(36.0, days=60.0)
    cls, peak = _classify(scenario, "ethanol_tolerance")
    assert cls == "ACTIVE", f"expected ACTIVE past the window, got {cls}"
    band_high = compile_scenario(scenario).parameters["ethanol_tolerance"].uncertainty.high
    assert peak > band_high, (
        f"peak ethanol {peak:.4g} no longer clears ethanol_tolerance's band high "
        f"{band_high}, so the high edge does not fire and the band still decides whether "
        "EthanolToleranceDeath exists -- the straddle window has not closed by 36 Brix"
    )
    assert peak > _TOL_HI


def test_a_dead_parameter_is_dead_across_its_whole_band_not_just_at_the_edges():
    """The three-point screen assumes monotone response. For these gates it is a *proof*
    (for ``q > theta``, inert at ``lo`` implies inert for every ``theta >= lo``; a
    multiplier on an identically-zero term is inert everywhere), but D-166 measured it at
    interior points rather than conceding the caveat in prose."""
    c = compile_scenario(_wine())
    events = tuple(c.events or ())
    grid = np.linspace(c.t_span_h[0], c.t_span_h[1], 121)
    base = c.parameters.resolve()
    nominal = simulate_scheduled(
        c.process_set, dict(base), c.y0, c.t_span_h, events=events, t_eval=grid
    )
    for name in ("ethanol_tolerance", "k_d2_ethanol_tolerance_death"):
        band = c.parameters[name].uncertainty
        for q in (0.0, 0.25, 0.5, 0.75, 1.0):
            p = dict(base)
            p[name] = band.low + q * (band.high - band.low)
            r = simulate_scheduled(c.process_set, p, c.y0, c.t_span_h, events=events, t_eval=grid)
            assert np.array_equal(r.y, nominal.y), (
                f"{name} is inert at both band edges but ACTIVE at {q:.0%} of the band -- "
                "the edge screen is not monotone here and the census undercounts."
            )


def _on_edge(scenario: Scenario) -> set[str]:
    """Live bands whose nominal coincides with one of their own band endpoints."""
    c = compile_scenario(scenario)
    return {
        n
        for n in c.parameters.names
        if c.parameters[n].uncertainty.high > c.parameters[n].uncertainty.low
        and c.parameters.value(n)
        in (c.parameters[n].uncertainty.low, c.parameters[n].uncertainty.high)
    }


def test_a_nominal_sitting_on_a_band_edge_cannot_be_classified_by_an_edge_screen():
    """A harness guard, not a model claim -- but a load-bearing one.

    For a parameter whose ``value`` equals a band endpoint, the run at that endpoint **is**
    the nominal run, so it returns bitwise identical *by construction* rather than by
    inertness. An edge screen then reports STRADDLE whenever the other edge moves, which is
    an artefact. The D-166 census classifies these over the band interior instead.

    Wine's three are all unread (hence DEAD anyway, so harmless). **Beer's are not**:
    ``Y_glycerol_sugar`` and ``Y_byproduct_sugar`` sit on the *low* edge at value 0.0, are
    read by an active Process, and are ACTIVE across the band interior -- exactly the case
    that would have manufactured two straddles. (These are the same two bands D-165 §2's
    ratio filter had to drop for a zero nominal.) Pinned so a new one cannot appear
    unnoticed.

    **``vant_hoff_co2_solubility`` joined at D-182, and its nominal sits on the low edge on
    purpose** — 2300 K is the source compilation's RECOMMENDED entry (a literature review),
    while 2400 K is what several independently measured entries carry, so the value has its
    own ground and the band's high edge has a different one. Moving the nominal to a
    midpoint purely to keep this set at three would be picking a number to satisfy a
    harness. What this test asks of it is met instead: anyone screening that parameter for
    inertness must classify over the band INTERIOR, because its low-edge run is the nominal
    run bitwise, by construction.

    NOTE this set now contains a READ parameter, which the old text's "wine's three are all
    unread (hence DEAD anyway, so harmless)" no longer covers: the CO2 coefficient is reached
    by every ``ph_of_state`` caller through ``PH_SYSTEM_READS``. So wine has acquired the
    same live case beer already had.

    **``bottling_burst_screwcap`` joined at D-187, and it is the cleanest example yet of why
    this set exists rather than being "fixed".** Its nominal sits on the HIGH edge because
    Lopes et al. 2007 prints ``<500 uL`` — a ceiling with no floor — so the entry ships the
    published bound and a CONSTRUCTED 0.0 low edge. The midpoint of ``[0, ceiling]`` was
    considered and rejected: it would invent a central estimate the source does not contain,
    purely to keep this set at four. The obligation this test names is met the same way as
    the CO2 coefficient's: anyone screening it must classify over the band interior. It is
    also not reached by any Process — the value seeds an event dose, so no sampler draws
    it."""
    wine_on_edge = _on_edge(_wine())
    assert wine_on_edge == {
        "f_non_ehrlich_phenylalanine",
        "copper_h2s_binding",
        "copper_mercaptan_binding",
        "vant_hoff_co2_solubility",
        "bottling_burst_screwcap",
    }, f"the wine nominal-on-edge set changed: {sorted(wine_on_edge)}"

    beer = Scenario(
        name="d166-beer",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 60.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[],
        duration_days=14.0,
    )
    beer_on_edge = _on_edge(beer)
    assert {"Y_glycerol_sugar", "Y_byproduct_sugar"} <= beer_on_edge, (
        "the two zero-nominal yield bands no longer sit on their band's low edge; the "
        f"beer nominal-on-edge set is {sorted(beer_on_edge)}"
    )
    cb = compile_scenario(beer)
    for n in ("Y_glycerol_sugar", "Y_byproduct_sugar"):
        assert cb.parameters.value(n) == cb.parameters[n].uncertainty.low == 0.0


@pytest.mark.parametrize("brix", [24.0, 26.0])
def test_peak_ethanol_is_the_quantity_the_window_is_measured_in(brix: float):
    """Ties the Brix window back to the ethanol margin it restates, so the two numbers
    cannot drift apart silently: below the window the peak is under the band low, inside it
    the peak is between the edges."""
    _, peak = _classify(_wine(brix), "ethanol_tolerance")
    if brix < _WINDOW_OPENS:
        assert peak < _TOL_LO
    else:
        assert _TOL_LO < peak < _TOL_HI
