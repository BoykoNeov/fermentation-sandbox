"""Beer's wort oxygen — the seed, the strip, and the inertness it is honest about (D-213).

Beer's ``o2`` slot existed from D-71 and **nothing seeded it**, so a beer ferment ran at exactly
0.000 mg/L dissolved oxygen from pitch to package — the one thing every brewery deliberately puts
into the wort. D-212 §3a found the gap while pricing an early-acetic candidate and could not use
it; this suite covers the Process that closes it.

**What these tests are FOR, since the beat ships an admittedly inert term.** The three O₂
consumers are aging-gated, so nothing reads this pool during fermentation. That makes the
interesting assertions the *negative* ones — this Process must not touch anything else, and the
pool must be gone before aging could start — plus the two positive ones that keep it from being
decoration: the seed is the sourced level at both band edges, and the oxygen genuinely leaves.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.kinetics.wort_oxygen import O2_SLOT, WortOxygenUptake
from fermentation.parameters import default_data_dir, load_parameters
from fermentation.scenario import (
    Intervention,
    Scenario,
    TemperaturePoint,
    compile_scenario,
)

#: The D-180/D-211 reference wort, so every number here is comparable with the beer acid suite.
TYRELL_SUGAR_GPL = 82.2388545
BEER_SCENARIO = {
    "glucose_gpl": 0.15 * TYRELL_SUGAR_GPL,
    "maltose_gpl": 0.70 * TYRELL_SUGAR_GPL,
    "maltotriose_gpl": 0.15 * TYRELL_SUGAR_GPL,
    "yan_mgl": 200.0,
    "pitch_gpl": 1.0,
    "initial_ph": 5.65,
}
#: BDF runs at rtol 1e-6; the mesh artifact below must sit inside it with room to spare.
SOLVER_RTOL = 1e-6


def _run(extra: dict[str, float] | None = None, days: float = 14.0, **solver: float):
    initial = dict(BEER_SCENARIO)
    if extra:
        initial.update(extra)
    compiled = compile_scenario(
        Scenario(
            name="d213",
            medium="beer",
            initial=initial,
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=days,
        )
    )
    res = compiled.run(**solver)
    return compiled, np.asarray(res.y, dtype=float), np.asarray(res.t, dtype=float)


def _at(schema, states, t_h, slot: str, hour: float) -> float:
    i = schema.slice(slot).start
    return float(np.interp(hour, t_h, states[i, :]))


@pytest.fixture(scope="module")
def beer_params():
    return load_parameters(default_data_dir() / "beer_generic.yaml")


def test_a_pitched_wort_starts_at_the_sourced_aeration_level(beer_params):
    """The seed IS the parameter, in mg/L, not a hard-coded number at the compile seam.

    Guards the units too: ``o2`` is a g/L slot and the source quotes mg/L, so a missing
    conversion would show up here as a 1000x error rather than as a plausible-looking curve.
    """
    compiled, states, t_h = _run(days=1.0)
    seeded_mgl = _at(compiled.schema, states, t_h, O2_SLOT, 0.0) * 1000.0
    assert seeded_mgl == pytest.approx(beer_params["o2_wort_aeration_beer"].value)
    assert seeded_mgl == pytest.approx(6.75)


def test_both_edges_of_the_aeration_band_seed_and_strip(beer_params):
    """BOTH printed edges, one assertion each — a nominal that passes says nothing about a band.

    ``o2_wort_aeration_beer`` is a printed range (5.5-8.0 mg/L) whose nominal is
    author-constructed, so the EDGES are the sourced content and are what a sampler will draw
    [[feedback-pin-the-band-not-the-nominal]].
    """
    band = beer_params["o2_wort_aeration_beer"].uncertainty
    assert (band.low, band.high) == (5.5, 8.0)
    for edge in (band.low, band.high):
        compiled, states, t_h = _run({"o2_mgl": edge}, days=2.0)
        assert _at(compiled.schema, states, t_h, O2_SLOT, 0.0) * 1000.0 == pytest.approx(edge)
        # Whatever the wort started with, the yeast strips it on the same lag-phase timescale.
        assert _at(compiled.schema, states, t_h, O2_SLOT, 24.0) * 1000.0 < 0.01


def test_the_yeast_strips_the_wort_oxygen_inside_the_lag_phase():
    """The sourced qualitative claim, as the only thing the rate constant has to deliver.

    *"The dissolved oxygen in beer rapidly disappears"* and pitched yeast takes *"several hours
    to adapt … before growth begins"* (The Chemistry of Beer). The assertion is deliberately a
    WINDOW, not a pinned curve: ``k_o2_uptake_beer`` is an author estimate and its whole band
    puts exhaustion inside this window, which is the only property anything depends on.
    """
    compiled, states, t_h = _run(days=2.0)
    schema = compiled.schema
    start = _at(schema, states, t_h, O2_SLOT, 0.0)
    assert _at(schema, states, t_h, O2_SLOT, 4.0) < 0.25 * start, "not stripped by 4 h"
    assert _at(schema, states, t_h, O2_SLOT, 12.0) < 0.01 * start, "still >1 % at 12 h"
    # Never MEANINGFULLY negative. The rate law is exponential decay and the Process clamps at
    # `o2 <= 0` besides, so the model cannot drive the pool negative — but the integrator can
    # still undershoot by an atol-sized amount on the way to zero, and it does (~1e-10 g/L
    # against atol 1e-9). Asserted against the seed rather than against 0, so this pins
    # "physically zero" and does not become a test of the integrator's rounding.
    assert states[schema.slice(O2_SLOT).start, :].min() > -1e-6 * start


def test_the_oxygen_is_gone_before_any_aging_could_begin():
    """The D-212 §7 hazard, defused and pinned rather than argued.

    A seed WITHOUT a sink would leave 6.75 mg/L lying in the pool, so a later ``begin_aging`` on
    a beer would oxidise against a phantom dose that no brewery's beer carries. This is the
    assertion that says the sink is why seeding is safe.
    """
    compiled, states, t_h = _run(days=14.0)
    assert _at(compiled.schema, states, t_h, O2_SLOT, 24.0) * 1000.0 < 1e-3


def test_the_process_writes_the_oxygen_slot_and_nothing_else():
    """Isolability at the DERIVATIVE, which is the claim that is exactly true.

    The integrated trajectory shifts by a solver-mesh amount (pinned in the next test); the RHS
    does not shift at all. Checked at several hours rather than one, because a Process that
    reads ``X`` could plausibly leak into it only once growth is running.
    """
    compiled, states, t_h = _run(days=7.0)
    schema, params = compiled.schema, compiled.parameters.resolve()
    process = WortOxygenUptake()
    assert process.touches == (O2_SLOT,)
    for hour in (0.0, 1.0, 4.0, 12.0, 24.0, 72.0):
        state = np.array([np.interp(hour, t_h, states[i, :]) for i in range(states.shape[0])])
        d = process.derivatives(hour, state, schema, params)
        written = [n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)]
        assert written in ([O2_SLOT], []), f"leaked into {written} at {hour} h"


def _worst_slot_difference(rtol: float, atol: float) -> tuple[float, str]:
    """Largest ABSOLUTE difference between the aerated and un-aerated beer, excluding ``o2``."""
    compiled_a, ya, ta = _run(days=14.0, rtol=rtol, atol=atol)
    _, yb, tb = _run({"o2_mgl": 0.0}, days=14.0, rtol=rtol, atol=atol)
    schema = compiled_a.schema
    grid = np.arange(0.0, 336.0, 1.0)
    worst, worst_slot = 0.0, ""
    for name in schema.names:
        if name == O2_SLOT:
            continue
        sl = schema.slice(name)
        for i in range(sl.start, sl.stop):
            a = np.interp(grid, ta, ya[i, :])
            b = np.interp(grid, tb, yb[i, :])
            gap = float(np.max(np.abs(a - b)))
            if gap > worst:
                worst, worst_slot = gap, name
    return worst, worst_slot


def test_seeding_the_oxygen_is_a_solver_MESH_artifact_and_not_a_coupling():
    """The inertness claim, separated from solver noise by CONVERGENCE rather than a bound.

    Seeding ``o2`` does shift every other column a little, and a fixed threshold could not tell
    you why: a real (small) coupling and an adaptive-mesh artifact both sit under any bound you
    pick generously enough to pass. **The two are distinguished by what happens when the solver
    tightens** — a physical pathway converges to a NON-zero difference, a mesh artifact converges
    to zero.

    Measured here across three decades of tolerance. At the shipped ``rtol`` 1e-6 the worst slot
    is ``h2s`` at ~1.1e-9 g/L — which looks alarming as a *relative* number (3.3e-4 of that
    pool's 3.4 µg/L peak) and is why this test does not use one. Tighten to 1e-9 and it falls to
    ~3.4e-12; to 1e-11 and it falls to ~1.1e-13. It tracks ``rtol`` essentially linearly, which
    is the signature of the mesh and not of a pathway.

    Consistent with the derivative test above: this Process writes ``o2`` only, and nothing reads
    ``o2`` while fermentation runs, so there is no channel for a coupling to exist in.
    [[feedback-pin-tolerance-vs-solver-tolerance]] and
    [[feedback-a-gate-is-a-discontinuity-the-solver-probes]].
    """
    loose, loose_slot = _worst_slot_difference(1e-6, 1e-9)
    tight, _ = _worst_slot_difference(1e-9, 1e-12)
    # A sanity CEILING, not the discriminator: 1e-4 g/L is 0.1 mg/L, far below anything
    # physically meaningful in any slot here, and comfortably below what a real O2->acid
    # coupling would move (D-212 priced the day-1 acetic requirement at 6-14 mg/L, i.e. 1e-2
    # g/L — two orders above this). The worst offender at the shipped tolerance is `S` at
    # ~1.4e-5 g/L, which is 1.7e-7 of an 82 g/L sugar charge.
    assert loose < 1e-4, f"{loose_slot} moved {loose:.2e} g/L even at the shipped tolerance"
    # The decisive assertion: two decades of rtol must buy at least one decade of agreement.
    # A genuine coupling would flatten out instead.
    assert tight < loose / 10.0, (
        f"difference did not converge with tolerance ({loose:.2e} -> {tight:.2e}); "
        "that is the signature of a real pathway, not the adaptive mesh"
    )


def test_removing_the_oxygen_seed_recovers_the_pre_d213_beer_exactly():
    """``o2_mgl = 0`` must reproduce the old anaerobic-from-t0 beer in the oxygen column itself.

    The isolability *gate* (prime directive 3), distinct from the test above: that one measures
    what the seed does to other slots, this one that the seed can be turned off at all.
    """
    compiled, states, t_h = _run({"o2_mgl": 0.0}, days=7.0)
    assert float(np.max(np.abs(states[compiled.schema.slice(O2_SLOT).start, :]))) == 0.0


def test_an_unpitched_wort_holds_its_oxygen():
    """The ``X <= 0`` branch — no yeast, no uptake. A positive control for the driver.

    Without it, "the oxygen disappears" would be consistent with a Process that simply decays the
    pool on a timer, and the biomass dependence would be untested decoration
    [[feedback-a-null-result-needs-a-positive-control]].
    """
    schema_carrier, states, t_h = _run(days=1.0)
    schema, params = schema_carrier.schema, schema_carrier.parameters.resolve()
    state = np.array([np.interp(0.0, t_h, states[i, :]) for i in range(states.shape[0])])
    state[schema.slice("X")] = 0.0
    d = WortOxygenUptake().derivatives(0.0, state, schema, params)
    assert float(d[schema.slice(O2_SLOT)][0]) == 0.0


def test_begin_aging_disables_the_wort_oxygen_uptake():
    """``begin_aging`` DISABLES this Process — the first thing that verb has ever switched off.

    Not a tidiness point: left enabled it competes with the aging oxidative sinks for a dosed
    ``add_oxygen`` and eats about half of it. That is a real defect the full suite caught —
    ``test_beer_depletes_its_packaging_oxygen`` saw 4.38 mg/L of an 8.0 mg/L dose survive to the
    post-dose sample — and packaging oxygen must be consumed by the sinks whose rates are
    calibrated against measured O₂ depletion, not by a fermentation-phase term carrying an
    author-estimate constant.

    Physically it is the same switch as the enables beside it: this Process is *pitched* yeast
    building membrane sterol in the lag phase, and past the ferment/aging breakpoint the yeast is
    settled and dormant.
    """
    compiled = compile_scenario(
        Scenario(
            name="d213-aging",
            medium="beer",
            initial=dict(BEER_SCENARIO),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=30.0,
            interventions=[Intervention(day=14.0, action="begin_aging", params={})],
        )
    )
    process_set = compiled.process_set
    assert process_set.is_enabled(WortOxygenUptake.name), "should run during fermentation"
    for event in compiled.events:
        if event.reconfigure is not None:
            event.reconfigure(process_set)
    assert not process_set.is_enabled(WortOxygenUptake.name), "begin_aging must disable it"


def test_wine_is_untouched_by_the_beer_wort_oxygen_process():
    """Beer-only: wine's O₂ is a post-ferment dose (D-71) and must not gain a fermentation sink.

    Checked through the registry rather than by running a wine, so it is a statement about
    WIRING and cannot pass because some wine scenario happened not to carry oxygen.
    """
    compiled = compile_scenario(
        Scenario(
            name="d213-wine",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    assert WortOxygenUptake.name not in compiled.process_set._processes
