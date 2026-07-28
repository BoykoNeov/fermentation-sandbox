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

import numpy as np
import pytest

from fermentation.core.media import get_medium
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
    ``sulfite_oxidation`` otherwise dilutes copper's share of total uptake (measured 61% -> 22%
    at 60 mg/L SO2). Guarding the sulfited arm would guard the better-behaved one.

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
def test_only_the_declared_o2_draw_responds_to_the_copper_multiplier(
    copper_dose_states, oxidative
):
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
