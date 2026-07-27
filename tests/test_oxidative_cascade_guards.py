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
three are the ones that must survive un-re-derived. The single seam the build is expected
to touch is :func:`_compile_with_old_oxidative_set`, and only to point it at the
``_OXIDATIVE_DIRECT_PROCESSES`` alternative once that alternative exists; the name list
and the pinned numbers below are not the build's to adjust. If a pin moves, that is a
finding for the decision record, not a tolerance to widen.

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

    **This function is the single seam the cascade build is expected to edit**, and only
    to select ``_OXIDATIVE_DIRECT_PROCESSES`` once that alternative exists. Today there is
    only one oxidative axis, so it is a plain compile. Everything else in this module —
    the name lists, the tolerances, the pinned numbers — is the guard, not the seam.
    """
    return compile_scenario(scenario)


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

#: Tolerance chosen for what it must CATCH, not for how tight it can be made. D-139 §3 asks
#: for the old set to reproduce the current trajectory "bit-for-bit"; that phrasing does not
#: survive contact with ``solve_ivp``, whose output is deterministic for identical inputs but
#: not reproducible to machine precision across scipy/BLAS builds — a 1e-12 pin becomes a
#: flake and then gets loosened, which is worse than no pin. The cascade moves these
#: quantities at the percent level, so 1e-6 detects it with four orders of margin.
_PIN_RTOL = 1e-6

#: ``phenylacetaldehyde`` sits at ~2e-9 g/L; a relative tolerance alone would be asserting
#: on solver noise for that slot.
_PIN_ATOL = 1e-15


@pytest.fixture(scope="module")
def wine_trajectory_run():
    compiled = _compile_with_old_oxidative_set(_wine_scenario(oak=True, amino_acids_gpl=0.5))
    return compiled, _integrate(compiled, days=_WINE_TOTAL_DAYS)


@pytest.mark.parametrize("medium", ["wine", "beer"])
def test_the_wired_oxidative_set_is_exactly_the_direct_sinks(medium):
    # Membership, asserted by NAME. This is what stops the cascade from being wired in
    # alongside the direct sinks and quietly double-counting the o2 pool: after the
    # rebuild, the direct alternative must present exactly these names and no others.
    process_set = get_medium(medium).build_process_set(strict=True)
    touching = {p.name for p in process_set.active if "o2" in p.touches}
    expected = BEER_OXIDATIVE_SINKS if medium == "beer" else OLD_OXIDATIVE_SINKS | O2_SOURCES
    assert touching == expected, (
        f"{medium}'s O2-touching Process set changed. Under the cascade this test must be "
        "satisfied by selecting the DIRECT alternative in "
        "_compile_with_old_oxidative_set — not by editing this expectation."
    )


def test_the_burst_oxidation_process_is_wired_into_no_medium():
    # NOT an endorsement — a pin on a documented discrepancy. D-138/D-139 describe "seven
    # O2 sinks"; six are wired. AntioxidantBurstOxidation (D-133) exists, is exported and
    # is listed in compile._AGING_GATED_PROCESSES, but no medium wires it, so its
    # burst_antioxidant pool is seeded and never drawn. If the cascade build wires it (or
    # deliberately retires it), this test is the one that should be updated, in the same
    # commit as the decision record that says which — so the change cannot pass silently.
    for medium in ("wine", "beer"):
        process_set = get_medium(medium).build_process_set(strict=True)
        assert "antioxidant_burst_oxidation" not in process_set


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
    assert actual == pytest.approx(expected, rel=_PIN_RTOL, abs=_PIN_ATOL), (
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
