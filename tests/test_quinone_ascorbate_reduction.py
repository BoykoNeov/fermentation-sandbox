"""Guards for the ascorbate route on the quinone node (decision D-202).

**What owes a guard here is different from every other cascade consumer, and the reason is the
default dose.** D-200 established that a new quinone consumer does not automatically owe an
assertion, because the benchmark's ``MIAO_BAND[0] <= nominal`` assert already detects anything
that moves the branching. That reasoning **fails here**, and it fails in the one direction that
matters: this Process is inert at its shipped default (``ascorbate = 0``), so *no shipped run
exercises it at all*. A defect anywhere in its rate law — the missing ``/ M_ASCORBIC`` unit
bridge is the obvious one, worth a factor of 176 — would leave the entire suite green.

So the guards below are the first in the cascade family whose job is simply **to run the thing**:

1. :func:`test_an_undosed_wine_is_exactly_inert` — the sourced default. UWC §24.4.3.2 says new
   wine has negligible ascorbic acid, which is *why* the slot defaults to 0; this asserts the
   model actually behaves that way, and it is the GREEN arm that makes the dosed arms below
   meaningful. It also pins the claim that a future beat must not quietly seed this pool the way
   D-134 seeded ``copper`` — that would move the benchmark by more than its whole margin.
2. :func:`test_the_dose_verb_lands` — ``add_ascorbate`` is the only way the pool becomes
   non-zero, so a broken verb silently disables the whole route.
3. :func:`test_the_dosed_node_share_and_pool_depletion` — **a magnitude pin, and it has to be
   one, for D-201's reason repeated**: dropping the molar bridge scales the quinone draw and the
   ascorbate draw by the same 176, so every structural invariant survives untouched while the
   route eats the node. A common factor is invisible to a ratio.
4. :func:`test_one_ascorbate_per_quinone` — the stoichiometry, and **unlike D-201's ratio guard
   this one is not blind**: ``_ASCORBATE_PER_QUINONE`` sits on the ascorbate side *only*, so
   mutating it does move the ratio. Stated as a contrast rather than copied as a formality.
5. :func:`test_the_route_is_absent_from_the_default_build` /
   :func:`test_the_route_does_not_run_before_begin_aging` — isolability (prime directive #3), the
   second because cascade Processes are re-enabled off a **fixed name tuple** in
   ``scenario.compile`` and one omitted from it is never disabled at all.
6. :func:`test_the_touches_and_reads_contracts` — ``touches`` is only enforced under
   ``strict=True``, and ``reads`` decides whether the new parameter is sampled by any ensemble.

**What is deliberately NOT guarded.** The beat's headline — that the model suppresses bisulfite's
share by only ~7 % where UWC reports ascorbate "completely prevent[s]" the reaction — is a
*finding about the branching structure*, not a property this Process should be pinned to. Pinning
it would freeze the very disagreement a later beat is meant to resolve.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.chemistry import M_ASCORBIC, M_QUINONE
from fermentation.core.kinetics import oxidative_cascade
from fermentation.core.process import ProcessSet
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario.compile import compile_scenario
from fermentation.scenario.schema import Intervention, Scenario, TemperaturePoint

ROUTE = "quinone_ascorbate_reduction"
FERM, SETTLE, AGE = 20.0, 10.0, 10.0
T0 = FERM + SETTLE
#: UWC §24.4.3.2's printed dose, and the one D-200/D-202 measure against. Corroborated as a
#: realistic magnitude by three other texts on disk: 50-100 mg/L at crushing (*Concepts in Wine
#: Chemistry*), 10 g/hL = 100 mg/L before bottling (*Handbook of Enology* Vol. 2), and 120 mg/L
#: alongside 30 mg/L SO2 (*Applied Wine Chemistry and Technology*).
DOSE_MGL = 60.0


def _scenario(*, dose: bool) -> Scenario:
    """A finished white wine, sulfited, oxygen-dosed under a sealed closure — with/without AA."""
    interventions = [
        Intervention(day=FERM - 1.0, action="add_so2", params={"so2_mgl": 80.0}),
        Intervention(day=FERM, action="begin_aging"),
        Intervention(day=T0, action="add_oxygen", params={"o2_mgl": 8.0}),
    ]
    if dose:
        # Dosed BEFORE aging begins, which is when a winemaker actually adds it.
        interventions.insert(
            1,
            Intervention(
                day=FERM - 1.0, action="add_ascorbate", params={"ascorbate_mgl": DOSE_MGL}
            ),
        )
    return Scenario(
        name=f"d202-ascorbate-{'dosed' if dose else 'undosed'}",
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
        interventions=interventions,
    )


def _run(*, dose: bool, oxidative: str = "cascade"):
    compiled = compile_scenario(_scenario(dose=dose), oxidative=oxidative)
    traj = simulate_scheduled(
        compiled.process_set,
        dict(compiled.param_values),
        compiled.y0.copy(),
        (0.0, (T0 + AGE) * 24.0),
        events=compiled.events,
        t_eval=np.linspace(0.0, (T0 + AGE) * 24.0, 3000),
    )
    return compiled, traj


def _series(compiled, traj, name: str) -> np.ndarray:
    return np.asarray(traj.y[compiled.schema.slice(name), :][0], dtype=float)


@pytest.fixture(scope="module")
def dosed():
    return _run(dose=True)


@pytest.fixture(scope="module")
def undosed():
    return _run(dose=False)


def test_the_route_is_absent_from_the_default_build(dosed):
    """Cascade-only. The default (direct) build must not carry the Process at all."""
    compiled_cascade, _ = dosed
    assert ROUTE in compiled_cascade.process_set._processes

    compiled_direct = compile_scenario(_scenario(dose=False), oxidative="direct")
    assert ROUTE not in compiled_direct.process_set._processes


def test_an_undosed_wine_is_exactly_inert(undosed):
    """The sourced default: new wine has negligible ascorbic acid, so nothing here runs.

    Two halves, because either alone is weak. The **pool** stays identically zero across a whole
    run — so nothing seeds it and nothing leaks into it — and the **derivative** the Process
    returns at a live oxidising state is an exact zero *array*, not merely a small number. The
    second half is what would catch a future beat wiring this route to some other reductant pool
    while leaving ``ascorbate`` at zero.
    """
    compiled, traj = undosed
    ascorbate = _series(compiled, traj, "ascorbate")
    assert np.all(ascorbate == 0.0), "an un-dosed wine grew ascorbate from somewhere"

    # A state where the rest of the cascade is demonstrably live.
    i = int(np.argmin(np.abs(traj.t - (T0 + 1.0) * 24.0)))
    y = traj.y[:, i]
    assert float(y[compiled.schema.slice("quinone")][0]) > 0.0, "no quinone: the arm is untested"

    proc = compiled.process_set._processes[ROUTE]
    d = proc.derivatives(0.0, y, compiled.schema, dict(compiled.param_values))
    assert not np.any(d), "the route is not exactly inert at the shipped default dose"


def test_the_dose_verb_lands(dosed, undosed):
    """``add_ascorbate`` is the only route into the pool, so a broken verb disables everything."""
    compiled, traj = dosed
    ascorbate = _series(compiled, traj, "ascorbate")
    pre = traj.t < (FERM - 1.0) * 24.0
    assert np.all(ascorbate[pre] == 0.0), "ascorbate appeared before its dose"

    i_dose = int(np.argmax(traj.t >= (FERM - 1.0) * 24.0))
    assert ascorbate[i_dose] == pytest.approx(DOSE_MGL * 1e-3, rel=1e-12)

    # ...and the dosed arm really is a different run from the un-dosed one, which is what makes
    # every comparison below meaningful (feedback-verify-the-restore-between-mutation-arms).
    compiled_u, traj_u = undosed
    assert _series(compiled, traj, "quinone")[-1] != _series(compiled_u, traj_u, "quinone")[-1]


def test_the_route_does_not_run_before_begin_aging(dosed):
    """Aging-gated, and the dose does not change that.

    Registered in the compile-seam tuple, so the Process is OFF until ``begin_aging`` even though
    the pool is already full — which is the case that tuple exists for: a wine dosed at crush must
    not have its ascorbate scavenged during fermentation.
    """
    compiled, traj = dosed
    ascorbate = _series(compiled, traj, "ascorbate")
    dosed_to_aging = (traj.t >= (FERM - 1.0) * 24.0) & (traj.t < FERM * 24.0)
    assert dosed_to_aging.any(), "the window between the dose and begin_aging is empty"
    assert np.all(ascorbate[dosed_to_aging] == pytest.approx(DOSE_MGL * 1e-3, rel=1e-12)), (
        "ascorbate was consumed before begin_aging"
    )


def test_the_dosed_node_share_and_pool_depletion(dosed):
    """Magnitude pin — a ratio guard cannot see the error this catches, and the envelope is
    sized against the parameter's own BAND rather than picked for looseness.

    Dropping the ``/ M_ASCORBIC`` molar bridge scales the quinone draw and the ascorbate draw by
    the SAME ~176, so the 1:1 stoichiometry below survives it untouched while this route takes
    essentially the whole quinone node. That is D-201's lesson repeated on a different constant,
    so the guard has to be a magnitude assertion. **Measured against the mutant: this test is the
    only one of the seven that goes red on it.**

    **The envelope is deliberately tight, and loosening it defeats the guard.** The dangerous
    error is not the gross one — it is the plausible near-miss of bridging with ``M_SO2``, the
    constant the sibling law legitimately uses, which is only **2.75x**. Measured on this exact
    scenario (probe5):

    ========================================  ==========  ============
    arm                                       consumed    node share
    ========================================  ==========  ============
    ``k_rel`` = 0.5 (declared band LOW)          2.582 %       3.990 %
    ``k_rel`` = 1.0 (nominal)                    4.915 %       7.651 %
    ``k_rel`` = 2.0 (declared band HIGH)         8.979 %      14.141 %
    bridged with ``M_SO2`` (2.75x)              11.605 %      18.396 %
    no bridge at all (176x)                     59.688 %      92.530 %
    ========================================  ==========  ============

    So the near-miss is separable from the whole declared band **only because the band (2x) is
    narrower than the mass ratio (2.75x)**, leaving a 1.29x window. The ceilings below sit inside
    it. A red here therefore means one of two things and **neither is fixed by relaxing the
    bound**: either the bridge is wrong, or ``k_so2_oxidation`` / the declared band has moved and
    the envelope must be **re-derived** from a fresh sweep (the D-154 idiom, and the same standing
    instruction the quinone-branching guard carries).
    """
    compiled, traj = dosed
    ascorbate = _series(compiled, traj, "ascorbate")
    i_open = int(np.argmin(np.abs(traj.t - T0 * 24.0)))
    consumed = 100.0 * (ascorbate[i_open] - ascorbate[-1]) / ascorbate[i_open]
    assert 1.0 < consumed < 10.5, f"ascorbate consumption {consumed:.4f} % is outside the envelope"

    # ...and the route's share of the quinone node, at a state where the node is live. The
    # missing-bridge version would read >90 % here.
    proc = compiled.process_set._processes[ROUTE]
    params = dict(compiled.param_values)
    i = int(np.argmin(np.abs(traj.t - (T0 + 1.0) * 24.0)))
    y = traj.y[:, i]
    mine = -float(
        proc.derivatives(0.0, y, compiled.schema, params)[compiled.schema.slice("quinone")][0]
    )
    total = 0.0
    for name, p in compiled.process_set._processes.items():
        if not compiled.process_set.is_enabled(name):
            continue
        draw = -float(
            p.derivatives(0.0, y, compiled.schema, params)[compiled.schema.slice("quinone")][0]
        )
        if draw > 0.0:
            total += draw
    share = 100.0 * mine / total
    assert 2.0 < share < 16.0, f"ascorbate takes {share:.4f} % of the quinone node"


def test_one_ascorbate_per_quinone(dosed):
    """1:1, sourced on BOTH sides — and unlike D-201's ratio guard, this one is not blind.

    ``_QUINONE_PER_H2S`` had to be kept off the reported quantity because the corpus does not
    settle it, which is exactly why a ratio could not police it. ``_ASCORBATE_PER_QUINONE`` is the
    opposite case: UWC §24.4.3.2 gives the fate as "reduction by ascorbic acid to regenerate the
    original o-diphenol", a two-electron reduction matched by ascorbate's two-electron oxidation
    to dehydroascorbate, so 1:1 is fixed on both sides at once. The coefficient therefore sits on
    the ascorbate side ALONE, a mutation of it does move this ratio, and the ratio is a real guard
    rather than a ceremonial one.
    """
    compiled, traj = dosed
    proc = compiled.process_set._processes[ROUTE]
    params = dict(compiled.param_values)
    i = int(np.argmin(np.abs(traj.t - (T0 + 1.0) * 24.0)))
    d = proc.derivatives(0.0, traj.y[:, i], compiled.schema, params)
    ascorbate_mol = -float(d[compiled.schema.slice("ascorbate")][0]) / M_ASCORBIC
    quinone_mol = -float(d[compiled.schema.slice("quinone")][0]) / M_QUINONE
    assert quinone_mol > 0.0, "the probe state has no quinone draw, so nothing is tested"
    assert ascorbate_mol == pytest.approx(
        quinone_mol * oxidative_cascade._ASCORBATE_PER_QUINONE, rel=1e-12
    )


def test_the_touches_and_reads_contracts(dosed):
    """The two contracts a new Process can satisfy *by assertion* and never have checked.

    ``touches`` is enforced only when the set is built with ``strict=True`` (CLAUDE.md), and
    ``compile_scenario`` builds a permissive set — so every other guard in this file would pass a
    Process writing outside its declaration. ``reads`` decides whether ``k_rel_ascorbate_quinone``
    is perturbed by any ensemble at all: a parameter no ``reads`` tuple claims is silently frozen
    (D-160, and the D-153→D-162 sampling work).
    """
    compiled, _ = dosed
    params = dict(compiled.param_values)
    schema = compiled.schema

    surface: set[str] = set()
    for p in compiled.process_set._processes.values():
        surface.update(p.reads)
    assert "k_rel_ascorbate_quinone" in surface, "the new parameter is frozen in every ensemble"

    proc = compiled.process_set._processes[ROUTE]
    assert not set(proc.reads) - set(params), "a declared read does not resolve"

    y = schema.zeros()
    y[schema.slice("T")] = 293.15
    y[schema.slice("quinone")] = 1.0e-3
    y[schema.slice("ascorbate")] = DOSE_MGL * 1e-3
    d = proc.derivatives(0.0, y, schema, params)
    written = {n for n in schema.names if float(np.atleast_1d(d[schema.slice(n)])[0]) != 0.0}
    assert written == {"quinone", "ascorbate"}, f"wrote outside its touches: {written}"

    # ...and the whole cascade integrated under the strict set, which is what enforces it.
    strict = ProcessSet(
        schema,
        list(compiled.process_set._processes.values()),
        modifiers=list(compiled.process_set._modifiers.values()),
        strict=True,
    )
    for name, on in compiled.process_set.enabled_snapshot().items():
        (strict.enable if on else strict.disable)(name)
    traj = simulate_scheduled(
        strict,
        params,
        compiled.y0.copy(),
        (0.0, (T0 + AGE) * 24.0),
        events=compiled.events,
        t_eval=np.linspace(0.0, (T0 + AGE) * 24.0, 600),
    )
    assert float(traj.y[schema.slice("ascorbate"), -1][0]) > 0.0
