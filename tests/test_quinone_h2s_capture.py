"""Guards for the quinone-H2S sulfide sink (decision D-201).

**What actually owes a guard here, and what does not.** D-200 established that a new quinone
consumer does *not* automatically owe an assertion: its ``MIAO_BAND[0] <= nominal`` benchmark
assert already detects anything that moves the branching, and a tighter pin would go brittle on
solver noise. That reasoning still holds — this Process takes 0.003 % of the quinone node and
moves the benchmark ratio by 0.0002 %, so **nothing here pins the branching**.

What does owe a guard is the part D-200 had no analogue for: this route's *own* reactant.

1. :func:`test_the_divalency_factor_cannot_move_the_sulfide` — **the beat's load-bearing design
   property.** The one unsourced number in the route (whether the divalent sulfide's adduct goes
   on to consume a *second* quinone) is written on the quinone side precisely so it cannot reach
   the reported quantity. That is a property of *how the rate law is written*, and the obvious
   "simplification" — making quinone primary and deriving the sulfide from it, which is how this
   beat's first probe wrote it — silently destroys it. Nothing else in the suite would notice.
2. :func:`test_the_nominal_sulfide_removal` — **a magnitude pin, and it has to be one.** The
   natural structural guard is the molar ratio between sulfide removed and quinone consumed, and
   that guard is *blind to the error that actually happened*: dropping the ``/ M_QUINONE`` unit
   bridge scales ``d(h2s)`` and ``d(quinone)`` by the same 108, so the ratio survives untouched
   while the sulfide draw is inflated a hundredfold. A common factor is invisible to a ratio.
   The inflated version reported 99.87 % removal and read as perfectly plausible prose.
3. :func:`test_the_sink_is_absent_from_the_default_build` /
   :func:`test_the_sink_does_not_run_before_begin_aging` — isolability (prime directive #3), and
   the second is not a formality: every cascade Process is disabled at the compile seam and
   re-enabled off a **fixed name tuple** in ``scenario.compile``. A Process omitted from that
   tuple is never disabled and therefore runs from t = 0, with nothing saying so.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.chemistry import M_H2S, M_QUINONE
from fermentation.core.kinetics import oxidative_cascade
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario.compile import compile_scenario
from fermentation.scenario.schema import Intervention, Scenario, TemperaturePoint

SINK = "quinone_h2s_capture"
FERM, SETTLE, AGE = 20.0, 10.0, 10.0
T0 = FERM + SETTLE


def _scenario() -> Scenario:
    """A finished white wine, sulfited, then a single oxygen dose under a sealed closure."""
    return Scenario(
        name="d201-quinone-h2s",
        medium="wine",
        strain="generic",
        initial={
            "brix": 21.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.3,
            "tannin_gpl": 0.3,
            "anthocyanin_gpl": 0.0,
            "amino_acids_gpl": 0.5,
        },
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=T0, celsius=20.0),
        ],
        closure="hermetic",
        duration_days=T0 + AGE,
        interventions=[
            Intervention(day=FERM - 1.0, action="add_so2", params={"so2_mgl": 80.0}),
            Intervention(day=FERM, action="begin_aging"),
            Intervention(day=T0, action="add_oxygen", params={"o2_mgl": 8.0}),
        ],
    )


def _run(oxidative: str = "cascade", **overrides: float):
    compiled = compile_scenario(_scenario(), oxidative=oxidative)
    params = dict(compiled.param_values)
    params.update(overrides)
    traj = simulate_scheduled(
        compiled.process_set,
        params,
        compiled.y0.copy(),
        (0.0, (T0 + AGE) * 24.0),
        events=compiled.events,
        t_eval=np.linspace(0.0, (T0 + AGE) * 24.0, 3000),
    )
    return compiled, traj


def _series(compiled, traj, name: str) -> np.ndarray:
    return np.asarray(traj.y[compiled.schema.slice(name), :][0], dtype=float)


@pytest.fixture(scope="module")
def cascade_run():
    return _run("cascade")


def test_the_sink_is_absent_from_the_default_build(cascade_run):
    """Cascade-only. The default (direct) build must not carry the Process at all."""
    compiled_cascade, _ = cascade_run
    assert SINK in compiled_cascade.process_set._processes

    compiled_direct = compile_scenario(_scenario(), oxidative="direct")
    assert SINK not in compiled_direct.process_set._processes


def test_the_sink_does_not_run_before_begin_aging(cascade_run):
    """Aging-gated: registered in the compile-seam tuple, so it is OFF until ``begin_aging``.

    Asserted on the trajectory rather than on the enabled flag, because the flag is what would
    be wrong if the registration were missed — checking it against itself proves nothing.
    Before ``begin_aging`` there is no quinone anyway, so the real content is the pairing: the
    sink is inert AND the pool is being driven by the fermentative Processes alone.
    """
    compiled, traj = cascade_run
    quinone = _series(compiled, traj, "quinone")
    pre = traj.t < FERM * 24.0
    assert np.all(quinone[pre] == 0.0), "quinone must be identically zero before begin_aging"

    proc = compiled.process_set._processes[SINK]
    params = dict(compiled.param_values)
    y_pre = traj.y[:, int(np.argmax(traj.t >= (FERM - 2.0) * 24.0))]
    d = proc.derivatives(0.0, y_pre, compiled.schema, params)
    assert float(d[compiled.schema.slice("h2s")][0]) == 0.0


def test_the_nominal_sulfide_removal(cascade_run):
    """Magnitude pin on the sulfide draw — a RATIO guard cannot see the error this catches.

    Dropping the ``/ M_QUINONE`` unit bridge scales the sulfide and quinone draws by the SAME
    factor, so every structural invariant survives it while the sulfide removal inflates ~108x.
    Only a magnitude assertion detects that, so this is a magnitude assertion.

    Banded loosely on purpose: it is guarding against a hundredfold error, not pinning a
    calibration, and a tight pin here would be brittle against the solver
    (``feedback-pin-tolerance-vs-solver-tolerance``).
    """
    compiled, traj = cascade_run
    h2s = _series(compiled, traj, "h2s")
    i_open = int(np.argmin(np.abs(traj.t - T0 * 24.0)))
    removed = 100.0 * (h2s[i_open] - h2s[-1]) / h2s[i_open]
    assert 1.0 < removed < 40.0, f"sulfide removal {removed:.4f} % is outside the sane envelope"


def test_the_divalency_factor_cannot_move_the_sulfide(monkeypatch):
    """**The guard this beat exists to leave behind.**

    ``_QUINONE_PER_H2S`` is the one number in the route the corpus does not give: H2S is
    divalent, so the mercapto-catechol adduct could in principle add to a second quinone, and
    *Understanding Wine Chemistry* is silent on it. The route is shippable anyway because the
    rate law is written **sulfide-first**, which confines that ambiguity to the quinone side —
    where it moves only a 0.003 % share of the node.

    So: doubling it must leave the sulfide **rate law** bitwise unchanged — and the *trajectory*
    unchanged only to within this Process's own share of the quinone node, which is the one
    indirect path by which the factor can reach the sulfide at all. Asserting bitwise equality on
    the trajectory would be asserting something false; the first version of this test did, and
    the 8.2e-07 it failed by is the real coupling, not noise.

    The mutation is also checked to LAND, so this cannot pass as a no-op confirming itself
    (``feedback-verify-the-restore-between-mutation-arms``).
    """
    compiled_1, traj_1 = _run("cascade")
    h2s_1 = _series(compiled_1, traj_1, "h2s")

    monkeypatch.setattr(oxidative_cascade, "_QUINONE_PER_H2S", 2.0)
    compiled_2, traj_2 = _run("cascade")
    h2s_2 = _series(compiled_2, traj_2, "h2s")

    # The mutation must LAND: the quinone draw at a fixed state doubles.
    proc = compiled_2.process_set._processes[SINK]
    params = dict(compiled_2.param_values)
    i = int(np.argmin(np.abs(traj_2.t - (T0 + 1.0) * 24.0)))
    y = traj_2.y[:, i].copy()
    d2 = proc.derivatives(0.0, y, compiled_2.schema, params)
    monkeypatch.setattr(oxidative_cascade, "_QUINONE_PER_H2S", 1.0)
    d1 = proc.derivatives(0.0, y, compiled_2.schema, params)
    q1 = float(d1[compiled_2.schema.slice("quinone")][0])
    q2 = float(d2[compiled_2.schema.slice("quinone")][0])
    assert q1 < 0.0, "the probe state must have a live quinone draw or the mutation is untested"
    assert q2 == pytest.approx(2.0 * q1, rel=1e-12), "the divalency factor did not land"

    # ...and having landed, it must not have reached the sulfide RATE LAW. This half IS exact:
    # at a fixed state the two arms compute the identical d(h2s)/dt, bit for bit.
    assert float(d1[compiled_2.schema.slice("h2s")][0]) == float(
        d2[compiled_2.schema.slice("h2s")][0]
    )

    # The TRAJECTORY, by contrast, is not bitwise identical, and asserting that it were would be
    # false. Doubling the quinone draw perturbs the quinone pool, and the sulfide rate is
    # bilinear in that pool — so the factor reaches the sulfide by exactly one indirect route,
    # bounded by this Process's own share of the quinone node. That is the same second-order
    # back-reaction measured at ``probe4``, and it is the honest bound to assert.
    #
    # Measured max relative difference: 8.2e-07 — below the node share (~3e-5) AND below the
    # solver's own rtol of 1e-6, i.e. not even resolvable. The tolerance is set an order above
    # the measurement and still two orders below the share, so the guard fails long before the
    # factor could matter, and does not fail on integrator noise.
    np.testing.assert_allclose(
        h2s_2,
        h2s_1,
        rtol=1e-5,
        atol=0.0,
        err_msg="the unsourced divalency factor reached the sulfide beyond its node share",
    )


def test_one_sulfide_per_event_is_the_sourced_half(cascade_run):
    """The sulfide side is 1:1 by *sourcing* — UWC assigns thiols to adduct formation, not redox.

    Pinned as the molar identity between the sulfide removed and the events implied by the
    quinone draw at ``_QUINONE_PER_H2S = 1``. This does NOT protect the unit bridge (see
    :func:`test_the_nominal_sulfide_removal` for why a ratio cannot) — it protects the
    stoichiometric coefficient itself.
    """
    compiled, traj = cascade_run
    proc = compiled.process_set._processes[SINK]
    params = dict(compiled.param_values)
    i = int(np.argmin(np.abs(traj.t - (T0 + 1.0) * 24.0)))
    d = proc.derivatives(0.0, traj.y[:, i], compiled.schema, params)
    h2s_mol = -float(d[compiled.schema.slice("h2s")][0]) / M_H2S
    quinone_mol = -float(d[compiled.schema.slice("quinone")][0]) / M_QUINONE
    assert h2s_mol > 0.0
    assert h2s_mol == pytest.approx(quinone_mol / oxidative_cascade._QUINONE_PER_H2S, rel=1e-12)
