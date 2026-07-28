"""Pre-build guards for the oxidative-cascade rebuild (decision D-139, Gate 3).

**Written BEFORE the cascade exists, while they still pass.** That ordering is the whole
point: D-139 closed the last of three gates with the instruction that three tests be
written first, because *"a guard written after the rebuild can only ratify whatever the
rebuild produced."* The rebuild re-homes every O2 sink behind an Fe(II)/O2 activation
node (D-137), adds a ``quinone`` slot and splits ``SulfiteOxidation`` in two (D-138).
Gate 3 (D-139) enumerated what must go red; these three tests cover the reds that would
otherwise go **silent** — passing while being wrong.

Baselines measured against the pre-build code by ``M:\\claud_projects\\temp\\d139-gate3\\``
(``baselines.py``, then ``measure_guards.py``, which re-took them without the
grid-dependent index anchoring of the first script).

**What each guard is really for**

1. :func:`test_beer_depletes_its_packaging_oxygen` — beer's O2 sink is not asserted by
   any other test in the suite. ``copper``/``burst_antioxidant`` are wine-only, no iron
   slot exists in either medium, and beer has no O2 *source* at all (``closure`` is
   rejected for beer), so beer's O2 is a one-shot dosed pool. If the activation node is
   made wine-only, beer's oxidation stops dead and **nothing else in the suite notices**.
   The signature is *"o2 stops declining"*, not *"flux goes to zero"* — hence the shape
   assertion alongside the magnitude.
2. :func:`test_old_oxidative_set_reproduces_its_trajectory` (+ the membership test) —
   prime directive #3. Isolability here is a **replacement**, not the additive
   ``_*_PROCESSES`` shape every other toggle uses: switching the cascade off must
   *restore* the direct sinks, not delete the oxidative axis. This pins what "restored"
   means, numerically, from before the cascade existed.
3. :func:`test_a420_baseline_survives_the_rebuild` — D-138 re-homes ``A420`` from an O2
   *yield* to a quinone *fate*. D-139 predicts the value survives while the meaning
   changes. That is the gate's own prediction, and this test is what makes it
   falsifiable.

**EDITING THESE TESTS IS THE THING THEY FORBID.** D-139's L6 already concedes that the
magnitude expectations in ``test_aging*.py`` and ``test_closure_ingress.py`` will be
re-derived by the build — which is exactly the moment they stop guarding anything. These
three are the ones that must survive un-re-derived. Exactly two seams exist for the build to
touch — :func:`_compile_with_old_oxidative_set` and :func:`_direct_oxidative_process_set` —
and only to point them at the ``_OXIDATIVE_DIRECT_PROCESSES`` alternative once it exists.
Every test here routes through one of them, so no test can be left with nowhere to be
pointed; the name lists, tolerances and pinned numbers are not the build's to adjust. If a
pin moves, that is a finding for the decision record, not a tolerance to widen.

**A correction this file records, found while taking the baselines.** D-138 and D-139
both speak of *"the seven O2 sinks"*, and D-139's isolability design specifies
``_OXIDATIVE_DIRECT_PROCESSES`` as *"the seven current O2 sinks, unchanged"*. The wired
count is **six**. :class:`~fermentation.core.kinetics.AntioxidantBurstOxidation` (D-133)
is defined, exported, unit-tested and listed in ``compile._AGING_GATED_PROCESSES`` — but
it is in **no medium's** ``process_factories``, so it is in no ``ProcessSet``, and the
``begin_aging`` enable loop skips it silently (it guards ``name in process_set``). Its
``burst_antioxidant`` slot is seeded from ``burst_antioxidant_initial`` and then consumed
by nothing. Verified by compiling wine *with* ``burst_antioxidant_gpl`` dosed and finding
the Process still absent. This is recorded rather than fixed because wiring it would move
every number pinned here, and D-138's constraint on where it belongs (a transient modifier
on the activation node, *not* a node of its own) is the rebuild's call to make.

**Superseded in part by D-147, and the part that changed matters to this file.** The call
was made: the burst is wired into a THIRD oxidative set, ``oxidative="direct_burst"``, and
is **not** in the default. So the paragraph above is still true of every build these 31
pins are taken on — the wired count in the DEFAULT set is still six — but "wired into no
medium at all" is no longer true of the package. The pins here did not move and were not
re-derived, which is precisely why the burst went opt-in: wiring it into the default would
have moved them by up to 37% (``o2`` -37.4%, ``A420`` -36.4% at 2 y), on the strength of a
joint calibration whose second constraint fails end-to-end. See
``tests/test_burst_oxidative_set.py`` for what wiring it does and what now forbids
re-describing it as a transient. The ``burst_antioxidant`` seed now follows its consumer,
so the "seeded and then consumed by nothing" defect is gone from the default build.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fermentation.core.media import get_medium
from fermentation.parameters import default_data_dir, load_parameters
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario.compile import CompiledScenario, compile_scenario
from fermentation.scenario.schema import Intervention, Scenario, TemperaturePoint

# ------------------------------------------------------------------------------------
# The oxidative axis as it is WIRED today — six sinks and one source.
# ------------------------------------------------------------------------------------

#: Every Process that DRAWS on the shared ``o2`` pool, by name, as wired into the wine
#: medium before the cascade. This list is the guard: after the rebuild the direct
#: alternative must present exactly these names again. ``antioxidant_burst_oxidation`` is
#: deliberately absent — it exists as a class but is wired into no medium (see the module
#: docstring); listing it here would make this test assert a thing the code has never done.
OLD_OXIDATIVE_SINKS = frozenset(
    {
        "oxidative_acetaldehyde",  # D-71 — O2 + ethanol -> acetaldehyde
        "phenolic_browning",  # D-74 — O2 + phenolics -> A420 pigment
        "sulfite_oxidation",  # D-72 — O2 + free SO2 -> sulfate (the one to be SPLIT)
        "strecker_degradation",  # D-75 — O2 + amino acids -> methional / phenylacetaldehyde
        "ellagitannin_oxidation",  # D-78 — O2 + oak ellagitannin (oak PROTECTION)
        "anthocyanin_fading",  # D-81 — O2 + anthocyanin -> faded_anthocyanin
    }
)

#: The only O2 *source* (D-136). Wine-only: ``closure`` is explicitly rejected for beer,
#: which is why beer's dosed O2 is a finite pool rather than a supplied flux.
O2_SOURCES = frozenset({"closure_oxygen_ingress"})

#: Beer wires the medium-agnostic subset only — no SO2/pH system (D-18), no amino acids
#: (D-32), no grape anthocyanin. Oak is available via ``add_oak``.
BEER_OXIDATIVE_SINKS = frozenset(
    {"oxidative_acetaldehyde", "phenolic_browning", "ellagitannin_oxidation"}
)


def _compile_with_old_oxidative_set(scenario: Scenario) -> CompiledScenario:
    """Compile ``scenario`` with the DIRECT (pre-cascade) oxidative sinks wired.

    **This function and its sibling below are the only seams the cascade build was expected
    to edit**, and only to select the direct alternative once it existed. D-141 built the
    cascade and made it the default wiring, so both seams now pass ``oxidative="direct"``.
    Everything else in this module — the name lists, the tolerances, the pinned numbers —
    is the guard, and none of it moved.
    """
    return compile_scenario(scenario, oxidative="direct")


def _direct_oxidative_process_set(medium: str):
    """The bare :class:`ProcessSet` for ``medium`` with the DIRECT oxidative sinks wired.

    The membership seam, sibling to :func:`_compile_with_old_oxidative_set`. It exists so
    that the membership test has somewhere to be pointed *at* after the rebuild: without
    it, once the cascade becomes a medium's default wiring, the test would go red with no
    seam to select, and the only move left would be editing the expectation — which is the
    single thing this module forbids. That is exactly what happened at D-141, and the seam
    did its job: selecting ``oxidative="direct"`` here was the whole fix.
    """
    return get_medium(medium, oxidative="direct").build_process_set(strict=True)


def _integrate(compiled: CompiledScenario, *, days: float, n: int = 4000):
    t_end = days * 24.0
    return simulate_scheduled(
        compiled.process_set,
        compiled.param_values,
        compiled.y0.copy(),
        (0.0, t_end),
        events=compiled.events,
        t_eval=np.linspace(0.0, t_end, n),
    )


def _at(trajectory, compiled: CompiledScenario, name: str, hours: float) -> float:
    """Read ``name`` at an explicit TIME, never a grid index.

    ``baselines.py`` anchored the post-dose reading at ``index + 2``, which silently moves
    with ``t_eval``; interpolating at a stated time keeps the pinned numbers meaningful if
    the grid is ever re-sized.
    """
    return float(np.interp(hours, trajectory.t, trajectory.y[compiled.schema.slice(name)][0]))


# ------------------------------------------------------------------------------------
# Guard 1 — beer's packaging oxygen (D-139 §2.1). The silent regression with no other
#           test in the suite watching it.
# ------------------------------------------------------------------------------------

_BEER_FERMENT_DAYS = 14.0
_BEER_AGE_DAYS = 120.0
_BEER_TOTAL_DAYS = _BEER_FERMENT_DAYS + _BEER_AGE_DAYS

#: ~air-saturation at packaging, the standard beer-staling entry point.
_O2_DOSE_MGL = 8.0

#: Measured on the pre-cascade code (``measure_guards.py``): plain beer consumes 5.7062
#: mg/L (71.36% of the dose) and oaked beer 6.3523 mg/L (79.44%). The ASSERTION below is
#: the floor D-139 states, not these values: re-homing the sinks behind an activation node
#: legitimately moves the exact figure, and pinning 5.7062 would manufacture a red the
#: build is right to argue with. What may NOT happen is beer's oxidation going quiet.
_BEER_MIN_CONSUMED_MGL = 5.0


def _beer_scenario(*, oak: bool) -> Scenario:
    interventions = [
        Intervention(day=_BEER_FERMENT_DAYS, action="begin_aging"),
        Intervention(day=_BEER_FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": _O2_DOSE_MGL}),
    ]
    if oak:
        interventions.insert(
            0,
            Intervention(
                day=_BEER_FERMENT_DAYS - 1.0,
                action="add_oak",
                params={"oak_gpl": 4.0, "toast": "medium"},
            ),
        )
    return Scenario(
        name="d139-guard-beer",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 150.0,
            "pitch_gpl": 0.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        duration_days=_BEER_TOTAL_DAYS,
        interventions=interventions,
    )


@pytest.fixture(scope="module")
def beer_runs():
    out = {}
    for oak in (False, True):
        compiled = _compile_with_old_oxidative_set(_beer_scenario(oak=oak))
        out[oak] = (compiled, _integrate(compiled, days=_BEER_TOTAL_DAYS))
    return out


@pytest.mark.parametrize("oak", [False, True], ids=["plain", "oaked"])
def test_beer_depletes_its_packaging_oxygen(beer_runs, oak):
    # Beer has no closure ingress and no other O2 source, so this is a one-shot pool: the
    # dose goes in at begin_aging and every mg/L that leaves it went through a sink. A
    # wine-only activation node would leave this flat at 8 mg/L for 120 days.
    compiled, trajectory = beer_runs[oak]
    t_dose = _BEER_FERMENT_DAYS * 24.0
    o2_post = _at(trajectory, compiled, "o2", t_dose + 1.0)
    o2_final = _at(trajectory, compiled, "o2", _BEER_TOTAL_DAYS * 24.0)

    # The dose landed: without this the depletion assertion could pass on an empty pool.
    assert o2_post == pytest.approx(_O2_DOSE_MGL * 1e-3, rel=1e-3)

    consumed_mgl = (o2_post - o2_final) * 1e3
    assert consumed_mgl >= _BEER_MIN_CONSUMED_MGL, (
        f"beer consumed only {consumed_mgl:.3f} mg/L of its {_O2_DOSE_MGL} mg/L packaging "
        f"O2 over {_BEER_AGE_DAYS:.0f} d (was 5.71 mg/L plain / 6.35 oaked before the "
        "cascade). If the activation node was made wine-only, beer has no metal to "
        "activate on and its oxidation has gone silent — see D-139 §2.1."
    )


@pytest.mark.parametrize("oak", [False, True], ids=["plain", "oaked"])
def test_beer_oxygen_only_ever_declines_after_the_dose(beer_runs, oak):
    # The SHAPE half of the guard. D-138 phrased beer's regression as "the flux goes to
    # zero"; because beer's O2 is a stock and not a flow, the observable signature is that
    # o2 stops declining. A rebuild that leaves a token sink alive could clear the
    # magnitude floor above while flat-lining partway; this catches monotonicity directly.
    compiled, trajectory = beer_runs[oak]
    t_dose = _BEER_FERMENT_DAYS * 24.0
    after = trajectory.y[compiled.schema.slice("o2")][0][trajectory.t >= t_dose + 1.0]
    assert after.size > 100
    # Strictly non-increasing, up to solver noise on a 1e-3 g/L pool.
    assert float(np.max(np.diff(after))) <= 1e-12


# ------------------------------------------------------------------------------------
# Guard 2 — the direct oxidative set, and the trajectory it must still produce
#           (D-139 §3, prime directive #3).
# ------------------------------------------------------------------------------------

_WINE_FERMENT_DAYS = 20.0
_WINE_YEARS = 2.0
_WINE_TOTAL_DAYS = _WINE_FERMENT_DAYS + 365.25 * _WINE_YEARS


def _wine_scenario(*, oak: bool, amino_acids_gpl: float) -> Scenario:
    """A red wine under natural cork — the closure supplies O2 continuously (D-136).

    ``oak`` and ``amino_acids_gpl`` are parameters because the *baseline* scenario D-139
    measured A420 on had neither, while the trajectory pin needs both: without amino acids
    ``strecker_degradation`` is a no-op and without oak ``ellagitannin_oxidation`` is, so
    two of the six sinks would be pinned at exactly 0.0 and would guard nothing.
    """
    initial: dict[str, float] = {
        "brix": 24.0,
        "yan_mgl": 200.0,
        "pitch_gpl": 0.25,
        "anthocyanin_gpl": 0.3,
        "tannin_gpl": 2.0,
    }
    if amino_acids_gpl:
        initial["amino_acids_gpl"] = amino_acids_gpl
    interventions = [
        Intervention(day=_WINE_FERMENT_DAYS - 1.0, action="add_so2", params={"so2_mgl": 60.0}),
        Intervention(day=_WINE_FERMENT_DAYS, action="begin_aging"),
    ]
    if oak:
        interventions.insert(
            0,
            Intervention(
                day=_WINE_FERMENT_DAYS - 1.0,
                action="add_oak",
                params={"oak_gpl": 4.0, "toast": "medium"},
            ),
        )
    return Scenario(
        name="d139-guard-wine",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=_WINE_TOTAL_DAYS,
        closure="natural_cork",
        interventions=interventions,
    )


#: Endpoints of the six wired sinks plus the shared pool they compete over, at 1 y and 2 y
#: of bottle aging, measured on the pre-cascade code. The per-sink O2 *split* is
#: deliberately NOT pinned — re-homing changes which Process draws the oxygen, which is the
#: rebuild's whole purpose; what must not change is the state the direct set arrives at.
_WINE_PINS: dict[str, tuple[float, float]] = {
    # slot: (value at 1 y, value at 2 y)
    "o2": (1.444070363339e-05, 1.576919862838e-05),
    "so2_total": (5.650150172989e-02, 5.286949820113e-02),
    "A420": (1.517312917092e-03, 2.586921779751e-03),
    "acetaldehyde": (4.856388639843e-05, 1.032235022293e-04),
    "methional": (2.394235591543e-07, 5.088022950883e-07),
    "phenylacetaldehyde": (1.084038787891e-09, 2.302542270316e-09),
    "anthocyanin": (4.364843106077e-05, 1.140357095547e-06),
    "faded_anthocyanin": (2.478027419280e-03, 2.479828995301e-03),
    "ellagitannin": (5.923487479412e-02, 5.997209704187e-02),
}

#: Tolerance chosen for what it must CATCH, and MEASURED rather than asserted.
#:
#: D-139 §3 asks the old set to reproduce the current trajectory "bit-for-bit". That phrasing
#: does not survive contact with the integrator: ``simulate_scheduled`` runs **BDF at
#: rtol=1e-6, atol=1e-9**, so a pin at 1e-6 asserts at exactly the solver's own error budget
#: and has no margin at all. The noise floor was measured (``probe_solver_noise.py``) by
#: re-integrating this scenario with the solver tightened four orders to rtol=1e-10 and
#: comparing — that difference IS the shipped run's error:
#:
#:   so2_total, o2                                  ~1e-8
#:   methional, ellagitannin, A420, acetaldehyde    1e-7 … 1e-6
#:   faded_anthocyanin                              6.1e-6
#:   phenylacetaldehyde                             6.7e-4     <- near-exhausted pool
#:   anthocyanin                                    2.9e-3     <- near-exhausted pool
#:
#: 1e-4 sits ~16x above the worst well-behaved slot and 100x below a 1% model change, which
#: is far smaller than anything the cascade does to these quantities.
#:
#: **The guard's resolution is bounded by the integrator, not by this constant** — tightening
#: it below ~1e-5 would pin solver noise, and adding the ``quinone`` slot alone perturbs step
#: selection (BDF's error norm is RMS-weighted over the state vector, so 93 -> 94 slots shifts
#: it) without any model change at all.
_PIN_RTOL = 1e-4

#: Two slots are near-exhausted pools where a RELATIVE tolerance is meaningless: ``anthocyanin``
#: has faded to 1.14e-6 g/L by 2 y and ``phenylacetaldehyde`` sits at ~2.3e-9 g/L, so the
#: solver's absolute error dominates their relative error. Each gets an absolute floor set
#: ~10x their measured absolute noise (3.3e-9 and 1.6e-12 respectively), which still leaves
#: the pin meaningful: 1e-7 is 0.2% of anthocyanin at 1 y, 1e-11 is 0.4% of
#: phenylacetaldehyde at 2 y. Every other slot is pinned on the relative tolerance alone.
_PIN_ATOL: dict[str, float] = {"anthocyanin": 1e-7, "phenylacetaldehyde": 1e-11}


@pytest.fixture(scope="module")
def wine_trajectory_run():
    compiled = _compile_with_old_oxidative_set(_wine_scenario(oak=True, amino_acids_gpl=0.5))
    return compiled, _integrate(compiled, days=_WINE_TOTAL_DAYS)


@pytest.mark.parametrize("medium", ["wine", "beer"])
def test_the_wired_oxidative_set_is_exactly_the_direct_sinks(medium):
    # Membership, asserted by NAME. This is what stops the cascade from being wired in
    # alongside the direct sinks and quietly double-counting the o2 pool: after the
    # rebuild, the direct alternative must present exactly these names and no others.
    process_set = _direct_oxidative_process_set(medium)
    touching = {p.name for p in process_set.active if "o2" in p.touches}
    expected = BEER_OXIDATIVE_SINKS if medium == "beer" else OLD_OXIDATIVE_SINKS | O2_SOURCES
    assert touching == expected, (
        f"{medium}'s O2-touching Process set changed. Under the cascade this test must be "
        "satisfied by selecting the DIRECT alternative in _direct_oxidative_process_set — "
        "not by editing this expectation."
    )


def test_the_burst_oxidation_process_is_wired_into_no_default_medium():
    # UPDATED AT D-147, which is what this test asked for: "if the cascade build wires it (or
    # deliberately retires it), this test is the one that should be updated, in the same commit
    # as the decision record that says which — so the change cannot pass silently."
    #
    # What changed is the STATUS of the absence, not the absence. D-133's Process is no longer
    # wired into *no medium at all* — `oxidative="direct_burst"` wires it (D-147), so it is
    # reachable, isolable and covered by tests/test_burst_oxidative_set.py. It stays out of the
    # DEFAULT build because only one of D-133's two joint constraints survives being run
    # end-to-end: the day-1 excess holds (0.93-0.95 vs ~1.0 mg/L/day), the self-exhaustion does
    # not, and under a natural cork the sink becomes a permanent ~37% tax on every oxidative fate
    # rather than the transient D-133 describes.
    #
    # So this remains a pin on a discrepancy — now a deliberate one — and it is what keeps the
    # 31 pins in this file meaningful: they are direct-set numbers, and they only stay direct-set
    # numbers while the burst is absent from the default.
    for medium in ("wine", "beer"):
        process_set = get_medium(medium).build_process_set(strict=True)
        assert "antioxidant_burst_oxidation" not in process_set
    # ...and the other half of the same fact, without which "absent from the default" would be
    # indistinguishable from "absent everywhere" — the state D-140 actually found.
    assert "antioxidant_burst_oxidation" in get_medium(
        "wine", oxidative="direct_burst"
    ).build_process_set(strict=True)


@pytest.mark.parametrize("medium", ["wine", "beer"])
def test_quinone_is_identically_zero_under_the_direct_set(medium):
    # D-139's ISOLABILITY CONDITION, made executable. The cascade's one new slot (D-141's
    # `quinone`, sized at D-138) is shared by both media and appended last to each, but while
    # the DIRECT oxidative alternative is selected nothing may produce it: "quinone must stay
    # identically 0 under the old set ... and if it does not, the sets are not isolable".
    #
    # Asserted STRUCTURALLY — no Process wired under the direct set declares `quinone` in its
    # `touches` — rather than by integrating and finding the pool still 0. A trajectory check
    # is the weaker claim: it can pass merely because the scenario never reached the producing
    # regime, and that false negative is self-sealing. `touches` is the contract ProcessSet
    # enforces under strict=True, so a Process that wrote quinone without declaring it would
    # already be a hard error. The trajectory form is asserted separately, on both media's
    # runs, by test_quinone_stays_zero_along_both_guard_trajectories.
    process_set = _direct_oxidative_process_set(medium)
    schema = get_medium(medium).schema
    assert "quinone" in schema, "the slot is shared — fork D1 — so BOTH media must carry it"
    touchers = {p.name for p in process_set.active if "quinone" in p.touches}
    assert touchers == set(), (
        f"{medium}: {sorted(touchers)} touch quinone while the DIRECT set is selected. The two "
        "oxidative sets are mutually exclusive (D-139) — this is a wiring error, not a "
        "tolerance to widen."
    )


def test_quinone_stays_zero_along_both_guard_trajectories(wine_trajectory_run, beer_runs):
    # The dynamic half of the isolability condition, over the two longest runs this file
    # already integrates: a 2 y corked red and a 120 d oaked beer, both with O2 flowing and
    # every direct sink active. Exact equality, not a tolerance — an untouched slot is not
    # merely small, it is bit-for-bit its seeded 0, and anything else means a Process is
    # writing a slot it did not declare.
    runs = [wine_trajectory_run, beer_runs[True]]
    for compiled, trajectory in runs:
        quinone = trajectory.y[compiled.schema.slice("quinone")][0]
        assert float(np.max(np.abs(quinone))) == 0.0


@pytest.mark.parametrize("slot", sorted(_WINE_PINS))
@pytest.mark.parametrize("years", [1.0, 2.0], ids=["1y", "2y"])
def test_old_oxidative_set_reproduces_its_trajectory(wine_trajectory_run, slot, years):
    # The strongest available evidence that the cascade did not quietly alter the validated
    # core: with the direct sinks selected, the model must still arrive at the state it
    # arrived at before the cascade was written.
    compiled, trajectory = wine_trajectory_run
    hours = (_WINE_FERMENT_DAYS + 365.25 * years) * 24.0
    expected = _WINE_PINS[slot][0 if years == 1.0 else 1]
    actual = _at(trajectory, compiled, slot, hours)
    assert actual == pytest.approx(expected, rel=_PIN_RTOL, abs=_PIN_ATOL.get(slot, 0.0)), (
        f"{slot} at {years:g} y moved from {expected:.6e} to {actual:.6e} under the direct "
        "oxidative set. The direct set is supposed to be untouched by the cascade — this is "
        "a finding for the decision record, not a tolerance to widen."
    )


# ------------------------------------------------------------------------------------
# Guard 3 — A420 keeps its value while its MEANING changes (D-139 §2.2).
# ------------------------------------------------------------------------------------

#: D-138 re-homes A420 from an O2 *yield* to a quinone *fate*, and D-139 predicts the
#: number survives. Measured on the pre-cascade code, on D-139's own baseline scenario
#: (no oak, no amino acids — oak competes for the same O2 and lowers A420 to 2.587e-3).
#:
#: The 2 y figure reproduces D-139's 0.002673 exactly. The 1 y figure does NOT reproduce
#: its 0.001522, and the discrepancy is an artefact of the receipts script rather than a
#: model difference: ``baselines.py`` read "1 y" as the midpoint INDEX of a grid spanning
#: the 20 d ferment plus 2 y, i.e. ~day 375 — about ten days short of one year after
#: begin_aging. A420 is still rising there, so its number is the smaller one. Interpolated
#: at the stated time, 1 y is 1.556718e-3.
_A420_PINS: dict[float, float] = {
    0.0: 0.0,
    1.0: 1.556717947222e-03,
    2.0: 2.672718096007e-03,
}


@pytest.fixture(scope="module")
def wine_a420_run():
    compiled = _compile_with_old_oxidative_set(_wine_scenario(oak=False, amino_acids_gpl=0.0))
    return compiled, _integrate(compiled, days=_WINE_TOTAL_DAYS)


@pytest.mark.parametrize("years", sorted(_A420_PINS), ids=["at_begin_aging", "1y", "2y"])
def test_a420_baseline_survives_the_rebuild(wine_a420_run, years):
    # This is the gate's OWN prediction under test. D-139 asserts A420's value survives the
    # re-homing while its meaning changes; nothing measured that, so it is asserted rather
    # than established. If the rebuild moves this number, the prediction was wrong and that
    # belongs in the decision record — along with the docstrings in analysis.py and media.py
    # that still describe A420 as an O2 yield, which D-139 §2.2 flags as the prose-vs-code
    # drift that recurs here by default unless someone acts.
    compiled, trajectory = wine_a420_run
    hours = (_WINE_FERMENT_DAYS + 365.25 * years) * 24.0
    actual = _at(trajectory, compiled, "A420", hours)
    assert actual == pytest.approx(_A420_PINS[years], rel=_PIN_RTOL, abs=1e-12)


# ------------------------------------------------------------------------------------
# Guard 4 — the copper multiplier's O2 BUDGET (D-139 §2.5, measured and closed at D-149).
# ------------------------------------------------------------------------------------

#: D-134 set ``k_copper_multiplier`` by ONE acceptance test: Ferreira 2015 caps the TOTAL
#: between-wine O2-consumption-rate spread from every compositional factor combined (Cu, Mn,
#: pH, TPI, phenolic acids…) at "never > factor 2.2", so copper ALONE must not over-spend it.
#: That is the arithmetic that forced the raw digitized 2000 L/g down to 600.
#:
#: ``aging.yaml``'s uncertainty note checks it on ``f_Cu`` **in isolation** — a ratio of the
#: multiplier, never of any run's O2 uptake, while Ferreira's 2.2x is a spread in MEASURED
#: TOTAL RATE. D-149 re-ran it in the frame that binds, and the isolated check turns out to be
#: exactly the CASCADE's number (copper multiplies 100% of that set's uptake) and a strict
#: upper bound for the direct sets (61% / 22% share, unsulfited / sulfited). So the guard is
#: written in the total-uptake form, which is the binding one under every wiring.
_FERREIRA_CU_LO_GPL = 1.68e-4  # Ferreira 2015 Table 1's own 15-red-wine copper range …
_FERREIRA_CU_HI_GPL = 6.79e-4  # … whose mean (2.61e-4) is `copper_typical`.
_FERREIRA_SPREAD_CEILING = 2.2

_COPPER_FERMENT_DAYS = 20.0
_COPPER_O2_DOSE_MGL = 8.0

#: The O2 draw ``k_copper_multiplier`` is allowed to move, per oxidative set. NOT an assertion
#: about which Processes *read* copper — it is measured by re-evaluating every active Process's
#: O2 contribution with the multiplier zeroed and collecting the ones whose draw moves, so a
#: Process that reached copper by a route nobody declared still shows up here.
_COPPER_MULTIPLIED_DRAWS: dict[str, frozenset[str]] = {
    "direct": frozenset({"phenolic_browning"}),  # D-134, on k_browning_eff
    "direct_burst": frozenset({"phenolic_browning"}),  # the burst sink is copper-free
    "cascade": frozenset({"oxygen_activation"}),  # D-141's re-home, onto the WHOLE node
}


def _copper_scenario(copper_gpl: float) -> Scenario:
    """A typical red (the D-132 phenolic anchor), unsulfited, dosed 8 mg/L O2 at aging.

    UNSULFITED deliberately: that is the worst case for the direct sets, whose copper-free
    ``sulfite_oxidation`` otherwise dilutes copper's share of total uptake (measured 61% -> 21%
    at 60 mg/L SO2). Guarding the sulfited arm would guard the better-behaved one.

    The acids are declared rather than defaulted. ``tartaric_gpl``/``malic_gpl`` both default to
    0.0, and a red with no titratable acidity solves to **pH 2.924** — outside any real red's
    3.4-3.8. It does not move this guard (unsulfited, so there is no pH-dependent term anywhere
    in the O2 draw, measured at exactly 1.0000x across pH 3.0-4.0), but calling that wine "a
    typical red" in the docstring of a guard about a real-wine calibration would be false, and
    the sulfited share it is contrasted against is genuinely pH-dependent.

    Runs only one day past the dose — the assertion is the INSTANTANEOUS draw at the dose, so
    integrating the aging tail would cost time and buy nothing.
    """
    return Scenario(
        name=f"d149-copper-{copper_gpl:.3e}",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.25,
            "anthocyanin_gpl": 0.3,
            "tannin_gpl": 2.0,
            "amino_acids_gpl": 0.5,
            "copper_gpl": copper_gpl,
            "tartaric_gpl": 6.0,
            "malic_gpl": 3.0,
            "initial_ph": 3.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=_COPPER_FERMENT_DAYS + 1.0,
        closure="hermetic",
        interventions=[
            Intervention(day=_COPPER_FERMENT_DAYS, action="begin_aging"),
            Intervention(
                day=_COPPER_FERMENT_DAYS,
                action="add_oxygen",
                params={"o2_mgl": _COPPER_O2_DOSE_MGL},
            ),
        ],
    )


def _o2_draws(
    compiled: CompiledScenario, params: dict[str, float], t: float, y: np.ndarray
) -> dict[str, float]:
    """Every active Process's own contribution to ``d[o2]/dt`` at ``(t, y)``, in g/L/h."""
    o2_slice = compiled.schema.slice("o2")
    draws: dict[str, float] = {}
    for process in compiled.process_set.active:
        contribution = np.asarray(process.derivatives(t, y, compiled.schema, params))
        value = float(contribution[o2_slice][0])
        if value != 0.0:
            draws[process.name] = value
    return draws


@pytest.fixture(scope="module")
def copper_dose_states():
    """``{(oxidative, copper_gpl): (compiled, params, t, y)}`` at the instant of the O2 dose.

    Both copper levels are integrated per set rather than evaluated off one shared state. The
    state at the dose IS copper-independent today (nothing reads ``copper`` until
    ``begin_aging`` enables the aging Processes, and the dose lands at that same instant) — but
    building that in would assume the very thing this guard exists to notice if it changes.
    """
    states: dict[tuple[str, float], tuple[CompiledScenario, dict[str, float], float, np.ndarray]]
    states = {}
    for oxidative in sorted(_COPPER_MULTIPLIED_DRAWS):
        for copper_gpl in (_FERREIRA_CU_LO_GPL, _FERREIRA_CU_HI_GPL):
            compiled = compile_scenario(_copper_scenario(copper_gpl), oxidative=oxidative)
            params = dict(compiled.param_values)  # a property: a FRESH dict per access (D-142)
            t_end = (_COPPER_FERMENT_DAYS + 1.0) * 24.0
            trajectory = simulate_scheduled(
                compiled.process_set,
                params,
                compiled.y0.copy(),
                (0.0, t_end),
                events=compiled.events,
                t_eval=np.linspace(0.0, t_end, 2001),
            )
            index = int(np.searchsorted(trajectory.t, _COPPER_FERMENT_DAYS * 24.0 + 1e-9))
            # trajectory.y is (n_states, n_times): a COLUMN is the state (D-147 amendment).
            states[(oxidative, copper_gpl)] = (
                compiled,
                params,
                float(trajectory.t[index]),
                np.asarray(trajectory.y[:, index], dtype=float),
            )
    return states


@pytest.mark.parametrize("oxidative", sorted(_COPPER_MULTIPLIED_DRAWS))
@pytest.mark.parametrize("at_band_high", [False, True], ids=["shipped", "band_high"])
def test_copper_may_not_over_spend_ferreiras_between_wine_spread(
    copper_dose_states, oxidative, at_band_high
):
    # WHAT THIS FORBIDS: raising `k_copper_multiplier` — or widening its band — until copper
    # ALONE accounts for more total-O2-uptake spread, across Ferreira's own 15-wine copper
    # range, than Ferreira measured for every compositional factor COMBINED.
    #
    # It is the live constraint, not a historical one. D-142 found Nguyen & Waterhouse 2021's
    # printed isolated-Cu table (Table 3.1, pH 3.5: 1.4e-4 -> 5.5e-4 /min over 0 -> 0.6355 mg/L
    # Cu, a 3.93x rise), and D-143 and D-148 both carried it forward as "a printed table arguing
    # k_copper_multiplier should go UP". D-149 converted it into the parameter: it implies
    # 2092 L/g, and at 2092 the cascade lands at 2.32x — over budget. Converted, the table
    # re-derives the value D-134 already rejected. It is not licence to move 600.
    #
    # Run at BOTH the shipped value and the band's declared HIGH edge, and the high edge is read
    # from the ParameterSet rather than hardcoded, so widening the band moves the guard with it
    # instead of leaving it pinned to a number the band no longer contains.
    compiled_lo, params_lo, t_lo, y_lo = copper_dose_states[(oxidative, _FERREIRA_CU_LO_GPL)]
    compiled_hi, params_hi, t_hi, y_hi = copper_dose_states[(oxidative, _FERREIRA_CU_HI_GPL)]
    if at_band_high:
        k_high = compiled_lo.parameters["k_copper_multiplier"].uncertainty.high
        params_lo = {**params_lo, "k_copper_multiplier": k_high}
        params_hi = {**params_hi, "k_copper_multiplier": k_high}

    rate_lo = -sum(_o2_draws(compiled_lo, params_lo, t_lo, y_lo).values())
    rate_hi = -sum(_o2_draws(compiled_hi, params_hi, t_hi, y_hi).values())
    assert rate_lo > 0.0 and rate_hi > 0.0, "no O2 is being consumed — the arm is not measuring"

    spread = rate_hi / rate_lo
    assert spread <= _FERREIRA_SPREAD_CEILING, (
        f"{oxidative}: copper alone moves total O2 uptake {spread:.4f}x across Ferreira's own "
        f"real-wine copper range ({_FERREIRA_CU_LO_GPL * 1e3:.3f}-{_FERREIRA_CU_HI_GPL * 1e3:.3f}"
        f" mg/L), over his {_FERREIRA_SPREAD_CEILING}x ceiling for EVERY compositional factor "
        "combined. That is the arithmetic that forced 2000 -> 600 at D-134, and a re-fit to "
        "Nguyen 2021's printed table lands at 2092 L/g and fails here. If this is red, the value "
        "or the band moved — that is a finding for the decision record, not a ceiling to raise."
    )


@pytest.mark.parametrize("oxidative", sorted(_COPPER_MULTIPLIED_DRAWS))
def test_only_the_declared_o2_draw_responds_to_the_copper_multiplier(copper_dose_states, oxidative):
    # WHAT THIS FORBIDS: silently re-homing copper onto another O2 sink, or adding a second one.
    # D-138 constraint 4 is that "a constant fitted to one structure does not survive being moved
    # to another" — D-141 moved copper from `phenolic_browning` (61% of unsulfited uptake) onto
    # `oxygen_activation` (100%), and the re-fit that implies went four decisions without being
    # run. This makes the next such move LOUD.
    #
    # It does NOT certify that the multiplied share is right for any set — the sibling test above
    # is what bounds that. It certifies only that the inventory is closed, and it is deliberately
    # weaker than that wording for exactly that reason.
    compiled, params, t, y = copper_dose_states[(oxidative, _FERREIRA_CU_HI_GPL)]
    draws = _o2_draws(compiled, params, t, y)
    zeroed = _o2_draws(compiled, {**params, "k_copper_multiplier": 0.0}, t, y)
    responders = {
        name
        for name, value in draws.items()
        if abs(value - zeroed.get(name, 0.0)) > 1e-18 * max(1.0, abs(value))
    }
    assert draws, "no active Process draws O2 — the arm is not measuring"
    assert responders == _COPPER_MULTIPLIED_DRAWS[oxidative], (
        f"{oxidative}: the O2 draws responding to k_copper_multiplier are {sorted(responders)}, "
        f"expected {sorted(_COPPER_MULTIPLIED_DRAWS[oxidative])}. Copper's calibration is scoped "
        "to what it multiplies (D-138 constraint 4) — re-homing or adding a draw invalidates it "
        "and needs a re-fit recorded, not this expectation edited."
    )


# ------------------------------------------------------------------------------------
# Guard 5 — the pH response's O2 BUDGET (D-150).
# ------------------------------------------------------------------------------------

#: Carrascon, Vallverdu-Queralt, Meudec, Sommerer, Fernandez-Zurbano & Ferreira (2018),
#: *Food Chemistry* 241:206-214, Tables 1 + 2 joined: eight commercial Spanish red wines
#: measured on the SAME repeated-air-saturation protocol Ferreira 2015 used (the protocol
#: D-132/D-133 are anchored on), as ``(code, pH, initial OCR, average OCR mg/L/day)``.
#:
#: This dataset is used rather than Ferreira's because it supplies the pH span and the
#: observed rate spread FROM THE SAME WINES — no mixing of one paper's span with another
#: paper's ceiling, which is the frame error D-149 was written to catch.
_CARRASCON_2018_REDS: tuple[tuple[str, float, float, float], ...] = (
    ("G1_09", 3.32, 3.70, 0.980),
    ("G2_13", 3.26, 4.90, 0.905),
    ("G3_14", 3.29, 1.95, 0.930),
    ("G4_14", 3.31, 4.29, 0.970),
    ("T1_11", 3.51, 4.40, 0.880),
    ("T2_12", 3.60, 5.13, 0.980),
    ("T3_10", 3.61, 6.83, 1.160),
    ("T4_14", 3.58, 5.29, 1.250),
)

#: The one statistic the paper prints against this exact pairing (Table 3): the correlation
#: between pH and the INITIAL OCR. It is the transcription's own check — see
#: ``test_the_carrascon_transcription_reproduces_its_papers_printed_correlation``. The "#"
#: it carries in the paper means p(t) < 0.1, i.e. NOT significant at 0.05: the paper's own
#: verdict is "only pH kept a non-significant positive correlation with initial OCRs".
_CARRASCON_PRINTED_R_PH_INITIAL = 0.689

#: Span and spread DERIVED from the transcription above, never hardcoded — editing a wine
#: moves both together, the way Guard 4 reads the copper band's high edge off the
#: ``ParameterSet`` instead of pinning a number the band may no longer contain.
_CARRASCON_PH_LO = min(w[1] for w in _CARRASCON_2018_REDS)
_CARRASCON_PH_HI = max(w[1] for w in _CARRASCON_2018_REDS)
#: **An OBSERVED spread, NOT a claimed ceiling.** Ferreira 2015's 2.2x is a stated maximum
#: over 15 wines ("never > factor 2.2"); this is merely what these 8 wines happened to show
#: on the STEADY (average) rate — the frame ``k_activation_*`` / ``k_browning_eff`` live in,
#: as opposed to the day-1 initial rate, which is ``burst_antioxidant``'s axis (D-133). It is
#: the right comparator because it comes from the same wines as the span, and it is a weaker
#: claim than Ferreira's; a green result here is NOT an endorsement of any pH response.
_CARRASCON_STEADY_SPREAD = max(w[3] for w in _CARRASCON_2018_REDS) / min(
    w[3] for w in _CARRASCON_2018_REDS
)

_PH_FERMENT_DAYS = 20.0
_PH_O2_DOSE_MGL = 8.0
#: Both arms are measured. Unsulfited is the worst case for the DIRECT sets (a pH term at the
#: ``f_Cu`` slot multiplies 62.4% of unsulfited uptake against 37.1% at 30 mg/L SO2, because
#: copper-free ``sulfite_oxidation`` dilutes it). The sulfited arm is the only one where the
#: model has ANY pH response today — ``k_activation_bisulfite * hso3`` through pH-dependent
#: ``free_bisulfite`` (D-149) — so dropping it would leave the one live route unguarded.
_PH_SO2_LEVELS = (0.0, 30.0)


def _ph_scenario(*, initial_ph: float, so2_mgl: float) -> Scenario:
    """The Guard 4 red, with ``initial_ph`` swept at a FIXED acid load.

    Only the back-solved strong-cation charge moves, so the sweep isolates pH rather than
    re-titrating the wine. The acids are declared for the reason D-149's amendment records:
    with ``tartaric_gpl``/``malic_gpl`` defaulted to 0 this wine solves to pH 2.924 and the
    solver then REFUSES a lower ``initial_ph`` outright.

    ``initial_ph`` anchors t=0, not the aging pH: byproducts accumulate through the ferment
    and pull it down, so 3.26/3.61 here arrive at the dose as 3.2084/3.5135. The guard asserts
    on the ratio between two runs of one sweep, so the offset cancels — but there is no way to
    SET an aging pH in this engine, and a future pH beat needs to know that.
    """
    interventions: list[Intervention] = [
        Intervention(day=_PH_FERMENT_DAYS, action="begin_aging"),
        Intervention(day=_PH_FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": _PH_O2_DOSE_MGL}),
    ]
    if so2_mgl > 0.0:
        interventions.insert(
            0,
            Intervention(day=_PH_FERMENT_DAYS - 1.0, action="add_so2", params={"so2_mgl": so2_mgl}),
        )
    return Scenario(
        name=f"d150-ph{initial_ph:g}-so2{so2_mgl:g}",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.25,
            "anthocyanin_gpl": 0.3,
            "tannin_gpl": 2.0,
            "amino_acids_gpl": 0.5,
            "tartaric_gpl": 6.0,
            "malic_gpl": 3.0,
            "initial_ph": initial_ph,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=_PH_FERMENT_DAYS + 1.0,
        closure="hermetic",
        interventions=interventions,
    )


@pytest.fixture(scope="module")
def ph_dose_states():
    """``{(oxidative, so2, initial_ph): (compiled, params, t, y)}`` at the instant of the dose."""
    states: dict[
        tuple[str, float, float], tuple[CompiledScenario, dict[str, float], float, np.ndarray]
    ]
    states = {}
    for oxidative in sorted(_COPPER_MULTIPLIED_DRAWS):
        for so2_mgl in _PH_SO2_LEVELS:
            for initial_ph in (_CARRASCON_PH_LO, _CARRASCON_PH_HI):
                compiled = compile_scenario(
                    _ph_scenario(initial_ph=initial_ph, so2_mgl=so2_mgl), oxidative=oxidative
                )
                params = dict(compiled.param_values)  # a property: a FRESH dict per access
                t_end = (_PH_FERMENT_DAYS + 1.0) * 24.0
                trajectory = simulate_scheduled(
                    compiled.process_set,
                    params,
                    compiled.y0.copy(),
                    (0.0, t_end),
                    events=compiled.events,
                    t_eval=np.linspace(0.0, t_end, 2001),
                )
                index = int(np.searchsorted(trajectory.t, _PH_FERMENT_DAYS * 24.0 + 1e-9))
                # trajectory.y is (n_states, n_times): a COLUMN is the state (D-147 amendment).
                states[(oxidative, so2_mgl, initial_ph)] = (
                    compiled,
                    params,
                    float(trajectory.t[index]),
                    np.asarray(trajectory.y[:, index], dtype=float),
                )
    return states


def test_the_carrascon_transcription_reproduces_its_papers_printed_correlation():
    # WHAT THIS FORBIDS: a silent transcription error in the eight wines above, which are the
    # only thing setting the span and the spread the sibling guard asserts on. The paper prints
    # exactly one statistic against this pairing (Table 3, r(pH, initial OCR) = 0.689), so
    # recomputing it checks the pH column and the initial-OCR column against the source rather
    # than trusting them — the D-149 precedent, where every (k1+k2)/2 entry in Nguyen's table
    # was reproduced from its own printed k1, k2 before the table was used.
    ph = [w[1] for w in _CARRASCON_2018_REDS]
    initial = [w[2] for w in _CARRASCON_2018_REDS]
    n = len(ph)
    mean_ph, mean_i = sum(ph) / n, sum(initial) / n
    cov = sum((a - mean_ph) * (b - mean_i) for a, b in zip(ph, initial, strict=True))
    var_ph = sum((a - mean_ph) ** 2 for a in ph)
    var_i = sum((b - mean_i) ** 2 for b in initial)
    r = cov / (var_ph * var_i) ** 0.5
    assert r == pytest.approx(_CARRASCON_PRINTED_R_PH_INITIAL, abs=0.005), (
        f"recomputed r(pH, initial OCR) = {r:.4f} against the paper's printed "
        f"{_CARRASCON_PRINTED_R_PH_INITIAL}. The transcription of Carrascon 2018 Tables 1+2 "
        "above is wrong, and it is what sets both the pH span and the steady-rate spread the "
        "sibling guard asserts on."
    )


@pytest.mark.parametrize("oxidative", sorted(_COPPER_MULTIPLIED_DRAWS))
@pytest.mark.parametrize("so2_mgl", _PH_SO2_LEVELS, ids=["unsulfited", "so2_30"])
def test_ph_may_not_out_spend_the_between_wine_steady_rate_spread(
    ph_dose_states, oxidative, so2_mgl
):
    # WHAT THIS FORBIDS: giving the O2 draw a pH dependence strong enough that pH ALONE moves
    # total uptake, across a real red's own pH span, further than every compositional factor
    # COMBINED moved the steady rate across the wines that span was measured on.
    #
    # It is green at 1.0000x (unsulfited) / ~1.006x (sulfited) today, because the model has no
    # pH term on the O2 draw at all — the D-140 "written while it still passes" pattern. What
    # makes it live is that the obvious source for such a term fails it: Nguyen & Waterhouse
    # 2021 Table 3.1's zero-copper column rises 3.87x per pH unit, which over this span is
    # 1.51x, and injected at the slot f_Cu occupies it lands at 1.5104x unsulfited / 1.5213x
    # sulfited in the CASCADE (where it would multiply ~100% of uptake) against an observed
    # 1.420x. In both direct sets the same exponent stays inside, at 1.32x / 1.15x. Identical
    # in shape to the copper failure Guard 4 forbids: one constant, over budget in one
    # structure and not another (D-138 constraint 4).
    #
    # It does NOT assert that the model's near-zero pH response is CORRECT — it is not; real
    # wine's steady rate rises ~1.62x per pH unit (D-150) and the model delivers ~1.006x. This
    # guard bounds the response from ABOVE only, and the under-response is recorded at D-150
    # as an open gap rather than guarded here, because nothing sources a within-wine pH rate
    # law to guard it against.
    compiled_lo, params_lo, t_lo, y_lo = ph_dose_states[(oxidative, so2_mgl, _CARRASCON_PH_LO)]
    compiled_hi, params_hi, t_hi, y_hi = ph_dose_states[(oxidative, so2_mgl, _CARRASCON_PH_HI)]
    rate_lo = -sum(_o2_draws(compiled_lo, params_lo, t_lo, y_lo).values())
    rate_hi = -sum(_o2_draws(compiled_hi, params_hi, t_hi, y_hi).values())
    assert rate_lo > 0.0 and rate_hi > 0.0, "no O2 is being consumed — the arm is not measuring"

    spread = rate_hi / rate_lo
    assert spread <= _CARRASCON_STEADY_SPREAD, (
        f"{oxidative} (SO2 {so2_mgl:g} mg/L): pH alone moves total O2 uptake {spread:.4f}x "
        f"across Carrascon 2018's own 8-red pH span ({_CARRASCON_PH_LO}-{_CARRASCON_PH_HI}), "
        f"past the {_CARRASCON_STEADY_SPREAD:.3f}x those same 8 wines showed on the STEADY "
        "rate from every compositional factor combined. Note this bound is an OBSERVED spread, "
        "not a claimed ceiling. If this is red, a pH term was added to the O2 draw: that needs "
        "a decision record saying which statistic it was sourced from (D-150 measured the "
        "steep published slopes to belong to the day-1 INITIAL rate, not the steady one) and "
        "how it separates from k_copper_multiplier -- which TWO independent experiments now say "
        "it does not: Nguyen 2021's own table (swing 2.011x per pH unit) and, since D-151, "
        "Carrasco-Quiroz 2022's copper-ORTHOGONAL L16 (swing 1.826x, F(1,38) = 5.95). See "
        "Guard 6."
    )


# ------------------------------------------------------------------------------------
# Guard 6 — pH x copper is NOT separable, corroborated in a copper-orthogonal design (D-151).
# ------------------------------------------------------------------------------------
#
# D-150 refused an ``f_pH`` term on the activation node on three legs. Leg 2 was the
# structural one: the architecture's rate modifiers are multiplicative and independent
# (D-10), so ``f_pH * f_Cu`` asserts the copper ratio is the same at every pH — and
# Nguyen & Waterhouse 2021 Table 3.1 shows it swinging 2.011x across one pH unit. The
# standing objection to that leg was that Nguyen's copper is dosed into the same wells the
# pH is set in, so the "interaction" could be an artefact of a non-orthogonal design.
#
# This section closes that objection with the unlock D-150's amendment named. Carrasco-Quiroz
# et al. 2022 crosses pH against copper in a genuine L16 orthogonal array, and the pH x Cu
# interaction is significant there too, in the same direction and at a comparable magnitude.
# **Leg 2 now rests on two independent experiments in two different media.**

#: Carrasco-Quiroz, Martinez-Gil, Nevares, Martinez-Martinez, Sanchez-Gomez &
#: del Alamo-Sanza (2022), *Foods* 11(13):1961, doi 10.3390/foods11131961 (PMC9266014).
#: Table 1, the L16(2^15) Taguchi design, as ``condition -> (pH, Fe, Cu, Mn, EtOH, AcH)``
#: in the paper's own units (mg/L except pH and %v/v). The MDPI host refuses automated
#: fetches; the PMC deposit and the EuropePMC ``fullTextXML`` both serve it and were
#: transcribed independently, agreeing cell-for-cell.
_CQ_DESIGN: dict[int, tuple[float, float, float, float, float, float]] = {
    1: (3.3, 1, 0.1, 4, 15, 30),
    2: (3.3, 8, 0.1, 1, 15, 30),
    3: (3.3, 8, 0.1, 4, 12, 10),
    4: (3.9, 1, 0.8, 1, 12, 10),
    5: (3.9, 8, 0.1, 1, 12, 30),
    6: (3.9, 1, 0.8, 4, 15, 30),
    7: (3.9, 1, 0.1, 1, 15, 10),
    8: (3.9, 8, 0.8, 4, 12, 10),
    9: (3.3, 8, 0.8, 1, 12, 30),
    10: (3.3, 8, 0.8, 4, 15, 10),
    11: (3.3, 1, 0.1, 1, 12, 10),
    12: (3.9, 1, 0.1, 4, 12, 30),
    13: (3.9, 8, 0.8, 1, 15, 30),
    14: (3.3, 1, 0.8, 4, 12, 30),
    15: (3.3, 1, 0.8, 1, 15, 10),
    16: (3.9, 8, 0.1, 4, 15, 10),
}
_CQ_FACTORS = ("pH", "Fe", "Cu", "Mn", "EtOH", "AcH")
#: (low, high) level of each factor exactly as printed in the Methods.
_CQ_LEVELS: dict[str, tuple[float, float]] = {
    "pH": (3.3, 3.9),
    "Fe": (1, 8),
    "Cu": (0.1, 0.8),
    "Mn": (1, 4),
    "EtOH": (12, 15),
    "AcH": (10, 30),
}
#: The three grape-extract wines the whole array was run on. They are a BLOCK, not a factor:
#: different extracts, same 16 conditions each.
_CQ_GEWS = ("A", "B", "C")

#: Table 2, ``R_max`` (hPa/h) — "maximum value of the oxygen consumption/rate curve".
#: A MAXIMUM rate on a single saturation, i.e. an initial-rate-class statistic in D-150 leg 3's
#: sorting, not the repeated-saturation steady rate the activation node is calibrated to (D-132).
_CQ_RMAX: dict[str, dict[int, float]] = {
    "A": {
        1: 5.6,
        2: 8.3,
        3: 8.3,
        4: 8.1,
        5: 10.0,
        6: 7.7,
        7: 19.7,
        8: 12.6,
        9: 7.6,
        10: 9.4,
        11: 5.6,
        12: 33.0,
        13: 12.5,
        14: 8.6,
        15: 6.7,
        16: 12.0,
    },
    "B": {
        1: 6.4,
        2: 7.9,
        3: 7.7,
        4: 7.2,
        5: 12.3,
        6: 8.1,
        7: 9.3,
        8: 13.0,
        9: 7.3,
        10: 8.8,
        11: 4.3,
        12: 30.7,
        13: 12.3,
        14: 7.9,
        15: 6.9,
        16: 13.2,
    },
    "C": {
        1: 4.6,
        2: 5.1,
        3: 5.3,
        4: 4.5,
        5: 8.1,
        6: 4.4,
        7: 6.2,
        8: 9.6,
        9: 5.0,
        10: 5.7,
        11: 3.0,
        12: 6.4,
        13: 7.4,
        14: 5.0,
        15: 4.4,
        16: 8.8,
    },
}
#: Table 2, ``delta O_90_10`` (hPa) and ``delta t_O_90_10`` (h). Carried ONLY to reproduce the
#: paper's printed significance verdict below. Their quotient is **not** a rate this archive
#: licenses: the numerator spans 1.49x (CV 0.074) across all 48 rows while the denominator
#: spans 4.79x (CV 0.398), so the quotient correlates with 1/dt at r = 0.988 — it is a
#: DURATION statistic wearing rate units, over one saturation of a dealcoholized grape-extract
#: reconstitution. See the verdict test's docstring.
_CQ_DO90: dict[str, dict[int, float]] = {
    "A": {
        1: 85.0,
        2: 96.9,
        3: 100.7,
        4: 101.6,
        5: 103.0,
        6: 112.1,
        7: 99.5,
        8: 106.1,
        9: 95.0,
        10: 95.8,
        11: 103.6,
        12: 97.6,
        13: 111.7,
        14: 102.7,
        15: 100.1,
        16: 107.8,
    },
    "B": {
        1: 101.6,
        2: 97.0,
        3: 100.7,
        4: 96.8,
        5: 112.0,
        6: 113.4,
        7: 97.0,
        8: 108.1,
        9: 93.7,
        10: 99.4,
        11: 94.3,
        12: 95.1,
        13: 112.7,
        14: 101.6,
        15: 104.8,
        16: 112.6,
    },
    "C": {
        1: 98.2,
        2: 94.4,
        3: 94.9,
        4: 94.3,
        5: 120.2,
        6: 102.9,
        7: 101.8,
        8: 110.1,
        9: 99.2,
        10: 95.5,
        11: 80.9,
        12: 101.1,
        13: 111.1,
        14: 96.8,
        15: 92.5,
        16: 104.8,
    },
}
_CQ_DT90: dict[str, dict[int, float]] = {
    "A": {
        1: 46.5,
        2: 29.0,
        3: 33.6,
        4: 36.8,
        5: 36.3,
        6: 35.5,
        7: 26.7,
        8: 17.8,
        9: 42.5,
        10: 27.3,
        11: 52.2,
        12: 17.4,
        13: 21.8,
        14: 31.9,
        15: 37.2,
        16: 21.1,
    },
    "B": {
        1: 42.4,
        2: 33.0,
        3: 35.1,
        4: 46.4,
        5: 39.4,
        6: 34.2,
        7: 35.1,
        8: 19.0,
        9: 43.2,
        10: 29.1,
        11: 59.9,
        12: 18.9,
        13: 20.8,
        14: 35.4,
        15: 42.8,
        16: 19.1,
    },
    "C": {
        1: 65.8,
        2: 52.1,
        3: 52.0,
        4: 68.5,
        5: 51.1,
        6: 77.3,
        7: 59.1,
        8: 28.8,
        9: 57.6,
        10: 47.7,
        11: 83.4,
        12: 44.7,
        13: 39.0,
        14: 61.5,
        15: 74.8,
        16: 31.8,
    },
}

#: The paper prints exactly two things this transcription can be checked against.
#: (1) Prose, on the pH main effect on ``R_max``: it rises "from 6.5 to 11.6 hPa/h".
_CQ_PRINTED_RMAX_PH_LEVEL_MEANS = (6.5, 11.6)
#: (2) A categorical verdict over all six factors: "pH, Fe2+ and Mn2+ being the significant
#: conditions" (ANOVA/LSD, p < 0.05).
_CQ_PRINTED_SIGNIFICANT = frozenset({"pH", "Fe", "Mn"})


def _cq_sign(cond: int, factor: str) -> int:
    """-1 at the factor's low level, +1 at its high level, as printed in Table 1."""
    lo, _hi = _CQ_LEVELS[factor]
    return -1 if _CQ_DESIGN[cond][_CQ_FACTORS.index(factor)] == lo else +1


def _cq_anova(values: dict[str, dict[int, float]]) -> dict[str, float]:
    """Six main effects + the pH x Cu interaction, over the 48 condition x GEw means.

    **This is NOT the paper's error term and no p-value from it is the paper's.** Table 2
    reports means of five replicates, so the replicate variance is not recoverable from the
    printed table; the residual here is that of the 48 condition x GEw means and absorbs every
    unmodelled interaction. That makes it a *conservative* error term for the effects it does
    fit, and it is used for exactly two things: reproducing the paper's categorical verdict as
    a transcription check, and estimating the pH x Cu interaction the paper never reports.
    """
    ys = [values[g][c] for g in _CQ_GEWS for c in _CQ_DESIGN]
    n = len(ys)
    grand = sum(ys) / n
    sst = sum((y - grand) ** 2 for y in ys)
    # GEw enters as a block: three different extracts, not a level of anything.
    ss_block = sum(
        len(_CQ_DESIGN) * (sum(values[g][c] for c in _CQ_DESIGN) / len(_CQ_DESIGN) - grand) ** 2
        for g in _CQ_GEWS
    )
    ss: dict[str, float] = {}
    for f in _CQ_FACTORS:
        lo = [values[g][c] for g in _CQ_GEWS for c in _CQ_DESIGN if _cq_sign(c, f) < 0]
        hi = [values[g][c] for g in _CQ_GEWS for c in _CQ_DESIGN if _cq_sign(c, f) > 0]
        ss[f] = (
            len(lo) * (sum(lo) / len(lo) - grand) ** 2 + len(hi) * (sum(hi) / len(hi) - grand) ** 2
        )
    ss_cells = 0.0
    for a in (-1, 1):
        for b in (-1, 1):
            cell = [
                values[g][c]
                for g in _CQ_GEWS
                for c in _CQ_DESIGN
                if _cq_sign(c, "pH") == a and _cq_sign(c, "Cu") == b
            ]
            ss_cells += len(cell) * (sum(cell) / len(cell) - grand) ** 2
    ss["pH_x_Cu"] = ss_cells - ss["pH"] - ss["Cu"]
    df_error = n - 1 - 2 - len(_CQ_FACTORS) - 1
    mse = (sst - ss_block - sum(ss[f] for f in _CQ_FACTORS) - ss["pH_x_Cu"]) / df_error
    out = {f"F_{k}": v / mse for k, v in ss.items()}
    out["df_error"] = float(df_error)
    return out


#: F(1, 38) at p = 0.05. Pinned rather than imported so the guard states its own threshold.
_CQ_F_CRIT = 4.0982


def test_the_carrasco_quiroz_l16_design_is_orthogonal_and_ph_x_cu_is_unaliased():
    # WHAT THIS FORBIDS: claiming a copper-ORTHOGONAL pH reading from a design that is not
    # orthogonal. The entire value of this dataset over Nguyen's is that copper varies
    # independently of pH; if the transcription of Table 1 were wrong, the pH x Cu interaction
    # the sibling guard measures could be any other factor's main effect wearing a disguise,
    # and leg 2's corroboration would be circular.
    for f in _CQ_FACTORS:
        assert sum(_cq_sign(c, f) for c in _CQ_DESIGN) == 0, (
            f"factor {f} is not balanced 8/8 across the 16 conditions — Table 1 is mis-transcribed"
        )
    for i, a in enumerate(_CQ_FACTORS):
        for b in _CQ_FACTORS[i + 1 :]:
            dot = sum(_cq_sign(c, a) * _cq_sign(c, b) for c in _CQ_DESIGN)
            assert dot == 0, f"{a} and {b} are not orthogonal (dot = {dot}); Table 1 is wrong"
    # The interaction column itself must be free of every studied main effect, or the
    # interaction estimate is confounded with one of them.
    for f in _CQ_FACTORS:
        dot = sum(_cq_sign(c, "pH") * _cq_sign(c, "Cu") * _cq_sign(c, f) for c in _CQ_DESIGN)
        assert dot == 0, (
            f"the pH x Cu interaction column is aliased with the {f} main effect (dot = {dot}). "
            "The interaction estimate below would then be unattributable, and leg 2's "
            "corroboration collapses."
        )


def test_the_carrasco_quiroz_transcription_reproduces_its_papers_printed_ph_level_means():
    # WHAT THIS FORBIDS: a silent transcription error in the R_max column, which is the only
    # column the load-bearing guard below reads. D-150 recorded this paper as "direction and
    # order only" precisely because its 6.5 -> 11.6 hPa/h comes from PROSE describing a figure.
    # Recomputing those two numbers from Table 2 is what converts it from prose into a sourced
    # table: the table reproduces the prose, so the prose was reporting the pH main effect.
    #
    # The match is 0.47% on the high level, NOT exact, and that is stated rather than rounded
    # away. No defensible aggregation lands on 11.6 -- arithmetic level means give 11.546,
    # geometric 10.116, mean-of-per-GEw-means 11.546. The residual is consistent with Table 2
    # being printed to one decimal place; it is not consistent with a different statistic.
    lo = [_CQ_RMAX[g][c] for g in _CQ_GEWS for c in _CQ_DESIGN if _cq_sign(c, "pH") < 0]
    hi = [_CQ_RMAX[g][c] for g in _CQ_GEWS for c in _CQ_DESIGN if _cq_sign(c, "pH") > 0]
    got = (sum(lo) / len(lo), sum(hi) / len(hi))
    for computed, printed in zip(got, _CQ_PRINTED_RMAX_PH_LEVEL_MEANS, strict=True):
        assert computed == pytest.approx(printed, rel=0.006), (
            f"recomputed R_max pH level means {got[0]:.4f} / {got[1]:.4f} against the paper's "
            f"printed {_CQ_PRINTED_RMAX_PH_LEVEL_MEANS}. The Table 2 R_max transcription is "
            "wrong, and it is what the separability guard below reads."
        )


def test_the_carrasco_quiroz_transcription_reproduces_its_papers_printed_significance_verdict():
    # WHAT THIS FORBIDS: a GROSS transcription error in Table 2's two duration columns. It
    # reproduces a CATEGORICAL verdict across all six factors -- "pH, Fe2+ and Mn2+ being the
    # significant conditions" -- which is why this dataset's transcription is trusted at all.
    #
    # ITS SENSITIVITY IS MEASURED, NOT ASSUMED, and it is coarse. The verdict's margins are
    # wide (pH 34.9, Fe 18.3, Mn 13.6 against Cu 0.23, EtOH 0.63, AcH 0.14 -- crit 4.10), so:
    #   - the smallest single-cell error it catches is 43% (dt B/13, 20.8 -> 11.9): a one-digit
    #     slip does NOT fire it;
    #   - a factor-aligned column distortion IS caught (dt scaled 1.5x on Fe-high drops Mn);
    #   - a GEw block permutation (A <-> C) and a halved condition-11 row are NOT caught.
    # These two columns feed nothing else -- the load-bearing R_max column is checked separately
    # above, and that check IS single-digit sensitive (33.0 -> 3.0 fires it). So this is a
    # coarse check on non-load-bearing columns, stated at the resolution it actually has.
    #
    # It runs on delta-O_90_10 / delta-t_O_90_10, which is the only statistic here that
    # reproduces the verdict. **That is a transcription check and nothing else.** That quotient
    # correlates with 1/dt at r = 0.988 (numerator CV 0.074, denominator CV 0.398): it is a
    # DURATION statistic in rate units, measured over ONE saturation of a dealcoholized
    # grape-extract reconstitution, and the activation node is calibrated to Ferreira's
    # repeated-saturation cycles-2-to-5 rate in real wine (D-132). Those are different physical
    # quantities, which is why Fe lights up here (F = 18.3) and vanishes on R_max (F = 0.017).
    #
    # DO NOT read a separability result off this statistic. Its pH x Cu interaction is
    # non-significant (F = 0.47), and D-151 records explicitly that this does NOT retire leg 2
    # for the node -- it is a limit on how far this dataset reaches, not a licence.
    rate = {g: {c: _CQ_DO90[g][c] / _CQ_DT90[g][c] for c in _CQ_DESIGN} for g in _CQ_GEWS}
    result = _cq_anova(rate)
    significant = frozenset(f for f in _CQ_FACTORS if result[f"F_{f}"] > _CQ_F_CRIT)
    f_by_factor = {f: round(result[f"F_{f}"], 2) for f in _CQ_FACTORS}
    assert significant == _CQ_PRINTED_SIGNIFICANT, (
        f"recomputed significant factors {sorted(significant)} against the paper's printed "
        f"{sorted(_CQ_PRINTED_SIGNIFICANT)} (F: {f_by_factor}, crit {_CQ_F_CRIT}). "
        "Table 2's delta-O_90_10 / delta-t_O_90_10 columns are mis-transcribed."
    )


def test_a_separable_f_ph_times_f_cu_is_rejected_by_the_copper_orthogonal_design():
    # WHAT THIS FORBIDS: reinstating the pH term D-150 refused by arguing that Nguyen's
    # pH/copper interaction is an artefact of his non-orthogonal dosing. It is not. Copper is
    # an independent factor here (proven by the orthogonality guard above), and the interaction
    # survives.
    #
    # The architecture's rate modifiers are multiplicative and independent (D-10), so any
    # f_pH * f_Cu asserts the pH ratio is IDENTICAL at every copper level. This measures that
    # assertion against the data: a separable model requires a swing of exactly 1.000x.
    #
    #   measured, R_max:  pH ratio 2.354x at Cu 0.1 mg/L vs 1.289x at Cu 0.8 mg/L
    #                     -> swing 1.826x over the 3.3-3.9 span, F(1,38) = 5.95, p = 0.0195
    #   Nguyen 2021:      swing 2.011x over one pH unit, opposite axis of the same table
    #
    # Same direction (pH's effect is weaker at high copper), comparable magnitude, two media,
    # two designs. Dropping condition 12 -- the one R_max outlier (33.0/30.7/6.4), which sits
    # in the pH-3.9/Cu-0.1 cell -- makes the interaction STRONGER, not weaker (F = 9.27,
    # p = 0.0044), so it is not an outlier artefact.
    #
    # Red-able on the real hazard: replacing the four cell means with their best multiplicative
    # -separable reconstruction drops F to 0.26, far below _CQ_F_CRIT. The statistic responds
    # to separability itself, not merely to the size of the numbers.
    result = _cq_anova(_CQ_RMAX)
    assert result["df_error"] == 38.0, "the error df moved; _CQ_F_CRIT is pinned at F(1, 38)"

    cells: dict[tuple[int, int], float] = {}
    for a in (-1, 1):
        for b in (-1, 1):
            xs = [
                _CQ_RMAX[g][c]
                for g in _CQ_GEWS
                for c in _CQ_DESIGN
                if _cq_sign(c, "pH") == a and _cq_sign(c, "Cu") == b
            ]
            cells[(a, b)] = sum(xs) / len(xs)
    swing = (cells[(1, -1)] / cells[(-1, -1)]) / (cells[(1, 1)] / cells[(-1, 1)])

    assert result["F_pH_x_Cu"] > _CQ_F_CRIT, (
        f"the pH x Cu interaction is F = {result['F_pH_x_Cu']:.3f} against F(1,38) crit "
        f"{_CQ_F_CRIT} (swing {swing:.4f}x). If this is green-turned-red by editing the table, "
        "the transcription is wrong. If it genuinely fell below the threshold, leg 2 of D-150 "
        "would rest on Nguyen alone again -- that is a finding for a decision record, not a "
        "threshold to lower."
    )
    assert swing > 1.0, (
        f"the pH effect is measured as STRONGER at high copper (swing {swing:.4f}x), reversing "
        "the direction Nguyen 2021 shows. Two experiments disagreeing in sign is a different "
        "finding from two agreeing, and needs a record."
    )


# ------------------------------------------------------------------------------------
# Guard 7 — the copper NULL turned into a BOUND, in the same design (D-152).
# ------------------------------------------------------------------------------------
#
# D-151 measured copper's own main effect in the only copper-orthogonal design in this
# archive's evidence set and found it non-significant: F = 0.23 on the duration statistic,
# F = 2.87 on R_max and pointing the WRONG way (more copper, LOWER max rate). It flagged
# that against D-134 and D-149 -- the two records that argued over whether
# ``k_copper_multiplier`` should go UP -- but refused to act, because *"it is a null, not a
# bound"*. Its Next named the fix: *"the SDs in Table 2 would support one, and that
# computation was not attempted here."*
#
# This section is that computation. A null says "no effect was detected"; a bound says "an
# effect this large is excluded", and only the second can be pointed at a constant.
#
# **THE FRAME IS THE SAME ONE THE CONSTANT WAS FITTED IN, WHICH IS WHY THIS TRANSFERS.**
# ``k_copper_multiplier`` was digitized from Danilewicz 2007 Figure 4 -- a MODEL WINE with
# FREE CuSO4 -- and then cut 2000 -> 600 with the explicit reason that *"free copper is more
# catalytically active, biasing the digitized slope high for real wine"*. This dataset is a
# dealcoholized grape-extract reconstitution with dosed CuCl2: the same class of medium, with
# the same free-copper bias, over a copper span (0.1-0.8 mg/L) that brackets Ferreira's own
# 15-red-wine range (0.168-0.679). So this is not a non-wine bound aimed at a wine value --
# it is a bound and a value in one frame, and the bound pushes FURTHER in the direction
# D-134's own correction already went.
#
# **AND IT IS LOOSE IN THE SAFE DIRECTION, TWICE.** ``R_max`` is a maximum rate on a single
# saturation, i.e. an initial-rate-class statistic in D-150 leg 3's sorting, while the
# multiplier scales a node calibrated to Ferreira's repeated-saturation cycles-2-to-5 steady
# rate (D-132). This archive's own pattern is that initial-rate statistics are MORE
# factor-sensitive than steady ones (Ferreira: ~15x initial between-wine spread against ~2.2x
# steady; D-150: every steep published pH slope turned out to be an initial-rate slope). That
# is an argument with its evidence, not a theorem -- but it points the same way as the
# free-copper bias, so a steady-rate bound in real wine should be TIGHTER than this one.

#: Table 2's ``+/- SD`` half, over FIVE replicates ("Five replicates were performed for each
#: one of the 16 experimental conditions and for each type of oxygen-saturated GEw").
#: NEW at D-152 -- D-151 transcribed the means only. Taken independently from the PMC deposit
#: and the EuropePMC ``fullTextXML``, which agree on all 48 cells and reproduce D-151's means
#: exactly; that agreement is the only provenance these SDs have, since the paper prints no
#: pooled statistic to check them against.
_CQ_RMAX_SD: dict[str, dict[int, float]] = {
    "A": {
        1: 0.8,
        2: 1.2,
        3: 0.6,
        4: 0.6,
        5: 0.5,
        6: 0.9,
        7: 1.0,
        8: 1.2,
        9: 0.6,
        10: 0.7,
        11: 0.7,
        12: 2.6,
        13: 1.0,
        14: 0.5,
        15: 0.2,
        16: 1.7,
    },
    "B": {
        1: 0.4,
        2: 1.2,
        3: 0.5,
        4: 0.8,
        5: 1.5,
        6: 1.0,
        7: 2.5,
        8: 0.9,
        9: 1.0,
        10: 0.4,
        11: 0.3,
        12: 0.7,
        13: 1.7,
        14: 0.3,
        15: 0.4,
        16: 1.1,
    },
    "C": {
        1: 0.4,
        2: 0.4,
        3: 0.6,
        4: 0.9,
        5: 1.1,
        6: 0.4,
        7: 1.8,
        8: 1.2,
        9: 0.8,
        10: 0.5,
        11: 0.4,
        12: 0.4,
        13: 0.7,
        14: 0.7,
        15: 0.4,
        16: 1.1,
    },
}
_CQ_N_REPLICATES = 5

#: Student's t at df = 38, TWO-SIDED 95%. Pinned rather than imported, like ``_CQ_F_CRIT``.
#: Two-sided deliberately, and this is the load-bearing methodological choice in the section:
#: copper's point estimate is NEGATIVE (ratio 0.789x marginal), so a one-sided 95% upper limit
#: lands at essentially zero effect and would "exclude" every candidate value down to and
#: including the shipped one -- purely because the noise happened to fall on that side of zero.
#: A bound that excludes everything because its point estimate went the other way is not a
#: measurement of anything.
_CQ_T_CRIT_38_TWO_SIDED = 2.0244

#: The two copper re-fits this archive has already rejected. Both were rejected on Ferreira's
#: 2.2x TOTAL between-wine spread budget -- an argument about how much room copper may take
#: from other factors. The bound below is an argument about copper's own measured size, from
#: a design in which copper varies independently. Two unrelated routes to one exclusion.
_D134_REJECTED_RAW_DANILEWICZ = 2000.0  # the raw Fig-4 digitization, cut to 600 at D-134
_D149_REJECTED_REFIT_NGUYEN = 2092.0  # what Nguyen 2021 Table 3.1 implies; D-149 closed it


def _cq_design_matrix(conds: list[int]) -> tuple[np.ndarray, list[str]]:
    """Intercept + 2 GEw block dummies + 6 main effects + pH x Cu, in +/-1 coding.

    Ten columns over 48 rows => df_error 38, i.e. **exactly the model behind the F(1, 38)
    Guard 6 reports**. On the full design this least-squares fit is algebraically identical to
    :func:`_cq_anova`'s level-mean sum-of-squares decomposition, and the guard below asserts
    that rather than assuming it. It is used here because a fit carries standard errors on its
    coefficients, which a sum-of-squares table does not.
    """
    names = ["intercept", "GEw_B", "GEw_C", *_CQ_FACTORS, "pH_x_Cu"]
    rows = [
        [
            1.0,
            1.0 if gew == "B" else 0.0,
            1.0 if gew == "C" else 0.0,
            *[float(_cq_sign(cond, f)) for f in _CQ_FACTORS],
            float(_cq_sign(cond, "pH") * _cq_sign(cond, "Cu")),
        ]
        for gew in _CQ_GEWS
        for cond in conds
    ]
    return np.asarray(rows), names


def _cq_fit(values: dict[str, dict[int, float]], conds: list[int]) -> dict[str, Any]:
    x, names = _cq_design_matrix(conds)
    y = np.asarray([values[g][c] for g in _CQ_GEWS for c in conds], dtype=float)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ beta
    df_error = len(y) - x.shape[1]
    return {
        "names": names,
        "beta": beta,
        "mse": float(resid @ resid) / df_error,
        "df_error": df_error,
        "xtx_inv": np.linalg.inv(x.T @ x),
    }


def _cq_copper_ratio_upper_bound(mse_override: float | None = None) -> dict[str, Any]:
    """The upper limit of a two-sided 95% interval on copper's effect on ``R_max``.

    **Taken at pH 3.3, the SIMPLE effect, not the pH-averaged main effect.** That is the
    conservative estimand and the reason is structural. D-151 measured a real pH x Cu
    interaction, so this design does not contain one copper ratio -- it contains two: 1.155x
    at pH 3.3 and 0.633x at pH 3.9. ``f_Cu`` is pH-independent (D-150 refused an ``f_pH``), so
    the model asserts a SINGLE ratio, and the data are internally inconsistent with any such
    ratio. **No single-ratio bound is fully coherent here, so the loosest slice is taken as the
    conservative resolution -- NOT because it is the better point estimate.** The pH-averaged
    marginal is not meaningless: real red wine's 3.4-3.8 straddles the design's levels, and
    interpolating to 3.6 gives ~0.894x, near the marginal's 0.789x. It is simply the arm that
    would let this dataset claim more than it can support -- it bounds k at 58 L/g, below the
    parameter's own band. D-152's amendment 1 prints all six arms rather than this one.

    Condition 12 -- the one ``R_max`` outlier (33.0/30.7/6.4) -- is KEPT, and the direction of
    that choice is the opposite of D-151's. There, dropping it made the interaction stronger.
    Here it sits in the pH-3.9/Cu-0.1 cell, so it does not touch the pH-3.3 point estimate at
    all (+0.9333 hPa/h either way) and only inflates the error term: dropping it would TIGHTEN
    this bound from 1.754x to 1.417x. Keeping it is both the wider interval and the model
    Guard 6 already ships.
    """
    conds = sorted(_CQ_DESIGN)
    model = _cq_fit(_CQ_RMAX, conds)
    names = model["names"]
    contrast = np.zeros(len(names))
    contrast[names.index("Cu")] = 2.0  # +/-1 coding: low -> high is two units
    contrast[names.index("pH_x_Cu")] = -2.0  # ... evaluated at the pH-3.3 level
    diff = float(contrast @ model["beta"])
    mse = model["mse"] if mse_override is None else mse_override
    se = float(np.sqrt(mse * (contrast @ model["xtx_inv"] @ contrast)))
    base = float(
        np.mean(
            [
                _CQ_RMAX[g][c]
                for g in _CQ_GEWS
                for c in conds
                if _cq_sign(c, "Cu") < 0 and _cq_sign(c, "pH") < 0
            ]
        )
    )
    return {
        "diff_hPa_h": diff,
        "se_hPa_h": se,
        "low_copper_reference_mean": base,
        "ratio_point": (base + diff) / base,
        "ratio_upper_bound": (base + diff + _CQ_T_CRIT_38_TWO_SIDED * se) / base,
        "mse": mse,
        "df_error": model["df_error"],
    }


def _cq_ratio_to_k(ratio: float, copper_typical_gpl: float) -> float:
    """The ``k_copper_multiplier`` implied by a given R_max ratio across 0.1 -> 0.8 mg/L Cu.

    ``f(Cu) = 1 + k * (Cu - copper_typical)``, so
    ``ratio = (1 + k*(Cu_hi - typ)) / (1 + k*(Cu_lo - typ))``, which inverts to the expression
    below. The conversion depends on where the mean-centering sits, so the guards run it
    across ``copper_typical``'s OWN declared band rather than at its point value alone.
    """
    d_hi = 0.8e-3 - copper_typical_gpl
    d_lo = 0.1e-3 - copper_typical_gpl
    return (ratio - 1.0) / (d_hi - ratio * d_lo)


@pytest.fixture(scope="module")
def aging_parameters():
    """``aging.yaml`` as a ``ParameterSet`` — the shipped copper value and band, from source."""
    return load_parameters(default_data_dir() / "aging.yaml")


def test_the_copper_orthogonal_bound_excludes_the_two_copper_refits_the_archive_rejected(
    aging_parameters,
):
    # WHAT THIS FORBIDS: re-opening `k_copper_multiplier` upward by arguing that the O2-budget
    # rejection of 2000/2092 was only ever a BUDGET argument -- that Ferreira's 2.2x ceiling is
    # about how much room copper may take from other factors, so a direct measurement of copper
    # could still license those values. It cannot. This bounds copper's own size in a design
    # where copper varies independently of everything else, and both rejected values sit
    # outside it. Guard 4 and this guard now reject them by two unrelated routes.
    #
    #   pH 3.3 simple effect on R_max:  +0.933 +/- 1.776 hPa/h on a 6.008 low-copper mean
    #   two-sided 95% upper limit:      ratio <= 1.754x across 0.1 -> 0.8 mg/L Cu
    #   converted:                      k_copper_multiplier <= 918 L/g
    #
    # The conversion is run across `copper_typical`'s whole declared band (0.168-0.679 mg/L),
    # over which the bound moves 1003 -> 663 L/g, and the exclusions must hold at the LOOSEST
    # end. Nothing is hardcoded: the bound is recomputed from the transcribed table and the
    # centring is read from the ParameterSet, so correcting either moves the guard with it.
    bound = _cq_copper_ratio_upper_bound()
    assert bound["df_error"] == 38, "the error df moved; _CQ_T_CRIT_38_TWO_SIDED is pinned at 38"

    # The fit must be the same model Guard 6's sum-of-squares route reports, or the two
    # readings of this one table have silently drifted apart.
    anova = _cq_anova(_CQ_RMAX)
    model = _cq_fit(_CQ_RMAX, sorted(_CQ_DESIGN))
    names = model["names"]
    for term in ("Cu", "pH_x_Cu"):
        i = names.index(term)
        f_from_fit = float(model["beta"][i] ** 2 / (model["mse"] * model["xtx_inv"][i, i]))
        assert f_from_fit == pytest.approx(anova[f"F_{term}"], rel=1e-9), (
            f"the least-squares fit and _cq_anova disagree on F_{term} ({f_from_fit:.6f} vs "
            f"{anova[f'F_{term}']:.6f}). They must be the same model -- Guard 6's F(1,38) and "
            "this bound's standard errors are taken from the same residual."
        )

    typical = aging_parameters["copper_typical"]
    k_bounds = {
        label: _cq_ratio_to_k(bound["ratio_upper_bound"], gpl)
        for label, gpl in (
            ("copper_typical band low", typical.uncertainty.low),
            ("copper_typical value", typical.value),
            ("copper_typical band high", typical.uncertainty.high),
        )
    }
    loosest = max(k_bounds.values())
    for name, rejected in (
        ("D-134's raw Danilewicz digitization", _D134_REJECTED_RAW_DANILEWICZ),
        ("D-149's re-fit from Nguyen 2021 Table 3.1", _D149_REJECTED_REFIT_NGUYEN),
    ):
        assert loosest < rejected, (
            f"{name} ({rejected:.0f} L/g) is NOT excluded by the copper-orthogonal bound: the "
            f"loosest bound across copper_typical's band is {loosest:.0f} L/g (ratio <= "
            f"{bound['ratio_upper_bound']:.4f}x; per-centring "
            f"{ {k: round(v) for k, v in k_bounds.items()} }). If this is red because the R_max "
            "transcription changed, the transcription is wrong. If the bound genuinely rose "
            "past a value D-134 and D-149 both rejected, that retires the only direct "
            "measurement of copper's size in this record -- a finding for a decision record, "
            "not a threshold to move."
        )


def test_the_shipped_copper_multiplier_is_not_excluded_by_the_copper_orthogonal_bound(
    aging_parameters,
):
    # WHAT THIS FORBIDS: raising `k_copper_multiplier` past the only direct, copper-orthogonal
    # measurement of copper's own size in this record. Guard 4 bounds copper by its O2 BUDGET
    # (how much of Ferreira's 2.2x between-wine spread copper alone may spend); this bounds it
    # by a measurement of copper itself. A value can clear the budget and still be excluded here.
    #
    # It asserts the shipped VALUE only. When written (D-152) that was forced: the band's HIGH
    # edge was 1500 L/g and already excluded, so asserting the band would have been red on
    # arrival. D-154 narrowed the edge to the bound, so the band is now assertable and the
    # sibling test below does exactly that. This one stays value-only so the two failure modes
    # stay distinguishable -- "the value moved" and "the edge moved" are different findings.
    param = aging_parameters["k_copper_multiplier"]
    typical = aging_parameters["copper_typical"]
    bound = _cq_copper_ratio_upper_bound()
    k_bound = min(
        _cq_ratio_to_k(bound["ratio_upper_bound"], gpl)
        for gpl in (typical.uncertainty.low, typical.value, typical.uncertainty.high)
    )
    assert param.value <= k_bound, (
        f"k_copper_multiplier = {param.value:.0f} L/g is EXCLUDED by the copper-orthogonal "
        f"bound of {k_bound:.0f} L/g (tightest across copper_typical's band; ratio <= "
        f"{bound['ratio_upper_bound']:.4f}x at pH 3.3, the loosest pH slice, two-sided 95%). "
        "Carrasco-Quiroz 2022's L16 is the only design in this record where copper varies "
        "independently of pH, and it is the same class of medium the value was digitized from. "
        "Raising the value past this needs a second copper-orthogonal experiment, not a "
        "re-reading of Nguyen's table (D-149). If THIS is red rather than its sibling below, "
        "the VALUE moved, not the band edge -- and the value is closed by D-149."
    )


def test_the_copper_bands_high_edge_is_not_excluded_by_its_own_bound(aging_parameters):
    # WHAT THIS FORBIDS: re-widening `k_copper_multiplier`'s band past the measurement that
    # set its edge -- the defect D-152 found and D-154 fixed, made unrepeatable.
    #
    # THE POINT OF THIS TEST IS THE SAMPLER, not the declaration. `runtime/ensemble.py` draws
    # `triangular(low, value, high)`, so a band edge is not documentation: every point up to it
    # is REACHABLE. Before D-154 the edge was 1500 L/g against a bound of 918, and ~29% of
    # ensemble draws landed on values this measurement excludes (D-152, flagged not fixed).
    # `rejected-values-must-be-unreachable`, live, in a field a green suite did not catch.
    #
    # WHY THE TIGHTEST CENTRING AND NOT THE SHIPPED 918. The ratio->k conversion depends on
    # where the mean-centering sits, `k_bound` is DECREASING in `copper_typical`, and
    # `copper_typical` is ITSELF a sampled band ([0.168, 0.679] mg/L). The two are drawn
    # INDEPENDENTLY, so the binding constraint is `copper_typical` at its MAXIMUM: 662.8 L/g,
    # not the 918 that holds only at the shipped centring. The sibling test above already took
    # `min()` across the centring band for this reason; this extends the same argument from the
    # value to the edge, which is the whole of D-154.
    #
    # Recomputed here, never trusted from the YAML note -- the D-118 breach-point template
    # (`test_the_de_novo_share_stays_above_its_analytic_breach_point`), because a claimed number
    # nobody recomputes is the D-96/D-109 defect class. It also catches the rounding trap that
    # would otherwise have shipped red: the exact bound is 662.802522 and D-152 PRINTED 663.
    param = aging_parameters["k_copper_multiplier"]
    typical = aging_parameters["copper_typical"]
    bound = _cq_copper_ratio_upper_bound()
    k_bound = min(
        _cq_ratio_to_k(bound["ratio_upper_bound"], gpl)
        for gpl in (typical.uncertainty.low, typical.value, typical.uncertainty.high)
    )
    # The binding centring must be copper_typical's HIGH edge -- if that ever stops being true
    # the argument above is wrong, and taking a min() would silently keep this test green.
    assert k_bound == _cq_ratio_to_k(bound["ratio_upper_bound"], typical.uncertainty.high), (
        "the tightest copper bound is no longer at copper_typical's band MAXIMUM, so k_bound is "
        "not monotone decreasing in the centring -- D-154's joint-sampling argument needs redoing"
    )
    assert param.uncertainty.high <= k_bound, (
        f"k_copper_multiplier's band HIGH edge ({param.uncertainty.high} L/g) is EXCLUDED by the "
        f"copper-orthogonal bound of {k_bound:.4f} L/g (at copper_typical's band maximum "
        f"{typical.uncertainty.high * 1e3:.3f} mg/L; ratio <= "
        f"{bound['ratio_upper_bound']:.4f}x at pH 3.3, the loosest pH slice, two-sided 95%). "
        "ensemble.py samples this band, so an excluded edge is REACHABLE, not merely declared. "
        "Re-widening needs a second copper-orthogonal experiment in real wine on a STEADY rate "
        "(D-152's named unlock) -- not a re-reading of Nguyen's table (D-149), and not this "
        "dataset read harder. NOTE the edge is deliberately truncated DOWN from the exact bound: "
        "D-152 printed 663 and the exact value is 662.802522, so 663 would fail here."
    )
    # ...and the band must still be a band: value strictly inside, low edge untouched by D-154.
    assert param.uncertainty.low < param.value < param.uncertainty.high


def test_the_printed_replicate_sds_cannot_carry_the_copper_bound():
    # WHAT THIS FORBIDS: rebuilding the bound above on the printed replicate SDs -- the obvious
    # thing to try, since D-151's own Next says "the SDs in Table 2 would support one", and
    # wrong by an order of magnitude in the DANGEROUS direction.
    #
    # A Table 2 cell is a mean of five replicates, so D-151's residual is
    # sigma^2_structure + sigma^2_rep/5. Decomposed:
    #
    #   pooled within-cell SD          1.006 hPa/h   (df 192)
    #   its share of a cell mean       0.202         (sigma^2_rep / 5)
    #   D-151's residual MSE          18.923
    #   => replicate noise is 1.1% of the residual; the rest is unmodelled structure
    #   => a replicate-only standard error understates by 9.7x
    #
    # So the SDs' real contribution is NEGATIVE: they prove they cannot carry the bound, which
    # is why the bound above uses the between-condition residual. That is not merely the
    # conservative choice, it is the only defensible one, and it is now measured rather than
    # asserted. A replicate-error bound reads ratio <= 0.983x -- i.e. it excludes the shipped
    # value, the whole band, and any positive copper effect at all -- on 1.1% of the variance.
    #
    # The reductio, on the SAME column and so independent of which column the paper's
    # categorical verdict came from: under a replicate-only error term copper's own F on R_max
    # is 269 against a critical value of 3.89. The paper reports copper as NOT significant.
    # Whatever error term produced that verdict, it was not these SDs.
    sds = [_CQ_RMAX_SD[g][c] for g in _CQ_GEWS for c in _CQ_DESIGN]
    assert len(sds) == 48, "the R_max SD transcription is not 16 conditions x 3 GEw"
    s2_rep_on_a_cell_mean = (sum(s * s for s in sds) / len(sds)) / _CQ_N_REPLICATES

    model = _cq_fit(_CQ_RMAX, sorted(_CQ_DESIGN))
    share = s2_rep_on_a_cell_mean / model["mse"]
    understatement = float(np.sqrt(model["mse"] / s2_rep_on_a_cell_mean))

    assert share < 0.10, (
        f"replicate noise is {share:.1%} of D-151's residual (sigma^2_rep/5 = "
        f"{s2_rep_on_a_cell_mean:.4f} against MSE {model['mse']:.4f}; a replicate-only standard "
        f"error would understate by {understatement:.2f}x). D-152's whole reason for "
        "using the between-condition residual is that this share is ~1%; if the two error terms "
        "have converged, the bound should be re-derived on the tighter one -- deliberately, in "
        "a record, not by editing this threshold."
    )
    # NOT a second threshold on the same quantity. `understatement` is sqrt(1/share), so any
    # assertion on it is implied by the one above and would be decoration in a file whose whole
    # discipline is that a guard forbids something. This asserts the independent fact instead:
    # under a replicate-only error term copper would be SIGNIFICANT on this column, while the
    # paper reports it as not. That is a statement about the SDs, not about their size.
    i_cu = model["names"].index("Cu")
    f_cu_on_replicates = float(
        model["beta"][i_cu] ** 2 / (s2_rep_on_a_cell_mean * model["xtx_inv"][i_cu, i_cu])
    )
    assert f_cu_on_replicates > _CQ_F_CRIT, (
        f"under a replicate-only error term copper's F on R_max is {f_cu_on_replicates:.1f}, "
        f"NOT past crit — so the reductio D-152 rests on has failed. The paper reports copper as "
        "non-significant; if the printed SDs can now reproduce that verdict, they may be the "
        "paper's own error term after all and the bound should be re-derived on them."
    )

    # And state by how much the choice of error term actually moves the thing that matters.
    on_residual = _cq_copper_ratio_upper_bound()["ratio_upper_bound"]
    on_replicates = _cq_copper_ratio_upper_bound(mse_override=s2_rep_on_a_cell_mean)[
        "ratio_upper_bound"
    ]
    assert on_replicates < on_residual, (
        f"the replicate-error bound ({on_replicates:.4f}x) is not tighter than the "
        f"between-condition one ({on_residual:.4f}x) -- the anti-conservatism this test exists "
        "to forbid has reversed sign, and D-152's choice of error term needs re-deriving."
    )
