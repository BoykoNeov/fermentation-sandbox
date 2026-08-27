"""The ``direct_burst`` oxidative set and the seed that follows it (decision D-147).

D-133 shipped :class:`~fermentation.core.kinetics.aging.AntioxidantBurstOxidation` — defined,
exported, unit-tested, listed in ``compile._AGING_GATED_PROCESSES`` — and wired into **no medium**.
It was therefore in no ``ProcessSet``, so the ``begin_aging`` enable loop skipped it silently (it
guards ``name in process_set``), and its ``burst_antioxidant`` pool was seeded from the sourced
``burst_antioxidant_initial`` and drawn by nothing. D-140 found that and pinned the absence rather
than fixing it, because *"wiring it would move every number pinned in this record, and which way it
goes is the rebuild's call."* D-147 made the call, having first measured what it costs.

**What the measurement found, because these tests exist to keep it true.** D-133's two constraints
JOINTLY pin ``k_burst_oxidation`` and ``burst_antioxidant_initial`` separately:

1. *The ~1.0 mg/L/day day-1 excess* — **HOLDS**. At Ferreira's own operating point (8 mg/L charge,
   hermetic, 20 C) the burst's own draw is 0.93-0.95 mg/L/day, -5 to -7% off target and well inside
   ``k_burst_oxidation``'s declared 0.5-5.0x band.
2. *The pool ~95% spent within one ~10-day saturation cycle* — **FAILS, structurally.** Spending
   the pool takes 3.3 mg/L of O2 through this route alone; the burst wins only ~35% of an 8 mg/L
   charge, so with SO2 present it PLATEAUS at 39.2% left because the oxygen runs out first. The
   isolated unit test that certifies this constraint
   (``test_burst_pool_exhausts_within_first_saturation_cycle``) passes only by setting
   ``o2 = 0.5 g/L`` — 500 mg/L, ~62x air saturation — and calls O2 running out a "confound". It is
   not a confound, it is the physics, and removing it removes the constraint's only content.

Under the DEFAULT operating point it is worse: a natural cork delivers 2.09 mg/L of O2 in two years
against a pool sized to absorb 3.3, so the *"fast, self-exhausting sink"* becomes a **permanent
~37% tax** on every oxidative fate. The self-exhausting SHAPE is a property of Ferreira's
saturation protocol, not of the sink — he delivers ~8 mg/L per 10-day cycle, a cork 2.09 mg/L per
two years, ~1400x apart.

So the burst ships **reachable, isolable and tested, but not default** — the D-141 cascade
treatment, for the same kind of reason. Measurements: ``M:\\claud_projects\\temp\\ferment\\
d147-burst-wiring\\`` (``FINDINGS.md`` plus the three scripts that produced it).

**These pins are the finding, not a convenience.** Two of them (the permanence guard and the
no-self-exhaustion guard) exist specifically to fail if someone re-describes the burst as a
transient without re-measuring it, which is the mistake D-133's own prose already made.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.media import get_medium
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario.compile import CompiledScenario, compile_scenario
from fermentation.scenario.schema import Intervention, Scenario, TemperaturePoint

BURST = "antioxidant_burst_oxidation"
OXIDATIVE_SETS = ("direct", "cascade", "direct_burst")

FERMENT_DAYS = 20.0
YEARS = 2.0
TOTAL_DAYS = FERMENT_DAYS + 365.25 * YEARS

# ======================================================================================
# WHY FOUR OF THE SIX NUMBERS BELOW MOVED AT D-182, AND WHAT WAS CHECKED FIRST
#
# D-182 put dissolved CO2 into the charge balance as carbonic acid. These tests state their
# own failure meaning — "the burst leaked into the default build" — so that had to be
# FALSIFIED before a single pin was touched, not assumed away:
#
#   * ``test_an_empty_burst_pool_reproduces_the_direct_trajectory_exactly`` — the bitwise
#     isolability guard — still PASSES. Direct and burst-with-an-empty-pool are byte-for-byte
#     the same model, which is the strongest available form of "the burst did not leak".
#   * The direct-vs-burst SEPARATION, computed from these very numbers, is unchanged: o2
#     37.2173 -> 37.2215 %, A420 35.8416 -> 35.8461 %, so2_total 2.2038 -> 2.2038 % at 1 y
#     (worst move 0.0045 percentage points). The signal these pins exist to protect is intact.
#
# WHAT ACTUALLY MOVED, and why it is small and confined: this scenario's wine supplies no
# ``initial_ph`` and no acids, so its charge balance carries only ``Byp`` against a zero
# cation — a state no real scenario produces, solving to a fictional pH ~2.92. The carbonic
# term shifts THAT by 0.0025 pH, which reaches these slots as ~2e-4 relative. An ANCHORED
# wine — the physical case — moves <=1e-5 across 400 days (D-182 measured it).
#
# ``so2_total`` is NOT re-pinned in either dict: it moved by 2.4e-6 / 4.7e-6, inside the
# tolerance, so its numbers are still D-140's and D-147's own. Only entries that actually
# broke are restated, each one recorded old -> new in the D-182 record.
# ======================================================================================

# ======================================================================================
# D-248 RE-PIN — o2 and A420 only, and the "did the burst leak?" reading is FALSIFIED first.
#
# This scenario doses ``amino_acids_gpl = 0.5``, so
# :class:`~fermentation.core.kinetics.amino_acids.AssimilableNitrogenUptake` is live in it: the
# yeast now consumes the must's assimilable nitrogen the way Crepin measures instead of stopping
# at what growth demands, which builds more biomass, ferments slightly differently, and reaches
# these aging slots. These tests state their own failure meaning — "the burst leaked into the
# default build" — so that had to be FALSIFIED before a pin was touched:
#
#   * ``test_an_empty_burst_pool_reproduces_the_direct_trajectory_exactly`` — the bitwise
#     isolability guard — still PASSES.
#   * The move is attributed by RE-RUNNING AT ``amino_acid_uptake_capacity_ratio = 0``, which
#     makes the new Process contribute exactly nothing while leaving the new state slot in the
#     schema. Every entry then sits within **7.5e-5** relative of its old pin — inside the 1e-4
#     tolerance, i.e. the extra schema dimension (which changes solve_ivp's RMS error norm and so
#     its step selection) does NOT on its own break anything. What breaks them is the Process,
#     and it moves them to **1.0-1.6e-4**. Small, one-sided, and confined to the two O2-driven
#     slots.
#   * ``so2_total`` is NOT re-pinned in either dict: it moved by 9e-6 to 2.1e-5, inside the
#     tolerance, so its numbers are still D-140's and D-147's own. Only entries that actually
#     broke are restated — the D-182 rule, applied again.
#
# The direct-vs-burst SEPARATION these pins exist to protect is intact: o2 37.22 -> 37.22 %,
# A420 35.85 -> 35.85 % at 1 y.
# ======================================================================================

#: D-140's own pins on the DIRECT set (``test_oxidative_cascade_guards._WINE_PINS``), restated here
#: as the baseline the burst set is measured AGAINST. Not re-derived — copied, so that if the
#: direct set ever moves, this file goes red for the same reason that one does.
DIRECT_PINS: dict[str, tuple[float, float]] = {
    # D-182: 1.444070363339e-05, 1.576919862838e-05 before the carbonic term.
    # D-248: 1.444363621566332e-05, 1.577233131562574e-05 before un-coupled nitrogen uptake.
    "o2": (1.4445900427430e-05, 1.5774539402510e-05),
    "so2_total": (5.650150172989e-02, 5.286949820113e-02),
    # D-182: 1.517312917092e-03, 2.586921779751e-03 before the carbonic term.
    # D-248: 1.5176242482004844e-03, 2.5874482824764554e-03 before un-coupled nitrogen uptake.
    "A420": (1.5178619521870e-03, 2.5878454466270e-03),
}

#: The same three slots under ``direct_burst``, measured at D-147 (``measure_burst.py``). The
#: SIGNAL here is 2-37%, and the measured solver noise floor is ~1e-8 (o2, so2_total) to ~7e-7
#: (A420) — taken by re-integrating the baseline at rtol 1e-8 on a 2x denser grid, because a pin
#: asserted at 1e-4 while BDF runs at 1e-6 is not self-evidently above its own noise.
BURST_PINS: dict[str, tuple[float, float]] = {
    # D-182: 9.066266268085e-06, 9.877702609069e-06 before the carbonic term.
    # D-248: 9.067499291264613e-06, 9.879107412627015e-06 before un-coupled nitrogen uptake.
    "o2": (9.0684747808570e-06, 9.8801657996040e-06),
    "so2_total": (5.774665845617e-02, 5.536926474690e-02),
    # D-182: 9.734830029596e-04, 1.644029205926e-03 before the carbonic term.
    # D-248: 9.73614572391898e-04, 1.6442540384438962e-03 before un-coupled nitrogen uptake.
    "A420": (9.7371569923880e-04, 1.6444273374770e-03),
}

_PIN_RTOL = 1e-4


def wine_scenario(*, burst_gpl: float | None = None) -> Scenario:
    """D-139/D-140's guard wine, so these numbers sit directly beside the pinned ones."""
    initial: dict[str, float] = {
        "brix": 24.0,
        "yan_mgl": 200.0,
        "pitch_gpl": 0.25,
        "anthocyanin_gpl": 0.3,
        "tannin_gpl": 2.0,
        "amino_acids_gpl": 0.5,
    }
    if burst_gpl is not None:
        initial["burst_antioxidant_gpl"] = burst_gpl
    return Scenario(
        name="d147-burst-set",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=TOTAL_DAYS,
        closure="natural_cork",
        interventions=[
            Intervention(
                day=FERMENT_DAYS - 1.0,
                action="add_oak",
                params={"oak_gpl": 4.0, "toast": "medium"},
            ),
            Intervention(day=FERMENT_DAYS - 1.0, action="add_so2", params={"so2_mgl": 60.0}),
            Intervention(day=FERMENT_DAYS, action="begin_aging"),
        ],
    )


def _integrate(compiled: CompiledScenario, *, n: int = 4000):
    t_end = TOTAL_DAYS * 24.0
    return simulate_scheduled(
        compiled.process_set,
        compiled.param_values,
        compiled.y0.copy(),
        (0.0, t_end),
        events=compiled.events,
        t_eval=np.linspace(0.0, t_end, n),
    )


def _at(traj, compiled: CompiledScenario, name: str, hours: float) -> float:
    """Read at an explicit TIME, never a grid index (the D-140 lesson)."""
    return float(np.interp(hours, traj.t, traj.y[compiled.schema.slice(name)][0]))


@pytest.fixture(scope="module")
def direct_run():
    compiled = compile_scenario(wine_scenario(), oxidative="direct")
    return compiled, _integrate(compiled)


@pytest.fixture(scope="module")
def burst_run():
    compiled = compile_scenario(wine_scenario(), oxidative="direct_burst")
    return compiled, _integrate(compiled)


# ------------------------------------------------------------------------------------
# Wiring — the thing D-133 believed it had done and had not.
# ------------------------------------------------------------------------------------


@pytest.mark.parametrize("oxidative", OXIDATIVE_SETS)
def test_only_the_burst_set_wires_the_burst_process(oxidative):
    # The whole point of D-147: the Process is now REACHABLE from a real registry, and reachable
    # from exactly one. Wine-only, because Ferreira's dataset is exclusively red wine and beer
    # carries no burst_antioxidant slot at all.
    wine = get_medium("wine", oxidative=oxidative).build_process_set(strict=True)
    beer = get_medium("beer", oxidative=oxidative).build_process_set(strict=True)
    assert (BURST in wine) is (oxidative == "direct_burst")
    assert BURST not in beer, "burst_antioxidant is wine-only; beer must never wire its consumer"


def test_the_burst_set_is_the_direct_set_plus_exactly_one_process():
    # An EXTENSION of the direct frame, not a third mechanism — which is why it is not an
    # alternative to the cascade. If it ever differs by more than this one name, the sets have
    # diverged and the "direct set stays byte-for-byte default" claim is no longer checkable.
    direct = {p.name for p in get_medium("wine", oxidative="direct").build_process_set().active}
    burst = {
        p.name for p in get_medium("wine", oxidative="direct_burst").build_process_set().active
    }
    assert burst - direct == {BURST}
    assert direct - burst == set()


def test_beer_burst_set_is_identical_to_beer_direct():
    # Asserted rather than assumed. The only thing this alternative adds is wine-only, so the two
    # beer builds are the SAME model — a real identity worth pinning, since a future beer-side
    # burst would otherwise arrive silently.
    direct = {p.name for p in get_medium("beer", oxidative="direct").build_process_set().active}
    burst = {
        p.name for p in get_medium("beer", oxidative="direct_burst").build_process_set().active
    }
    assert direct == burst


def test_the_burst_set_wires_no_cascade_node():
    # D-138 ruled the burst's node in the cascade frame "UNDETERMINED — do not force a node": it
    # was fitted to a day-1 O2 FLUX, and re-homed as an H2O2/quinone consumer the O2 spike vanishes
    # because only the activation node consumes O2 at all there. So there is deliberately no
    # `cascade_burst`, and this set must not have grown half of one.
    burst = {
        p.name for p in get_medium("wine", oxidative="direct_burst").build_process_set().active
    }
    cascade = {p.name for p in get_medium("wine", oxidative="cascade").build_process_set().active}
    direct = {p.name for p in get_medium("wine", oxidative="direct").build_process_set().active}
    assert burst & (cascade - direct) == set()


# ------------------------------------------------------------------------------------
# The seed follows its consumer (the M0 defect).
# ------------------------------------------------------------------------------------


@pytest.mark.parametrize("oxidative", OXIDATIVE_SETS)
def test_the_burst_pool_is_seeded_only_where_something_can_draw_it(oxidative):
    # D-133 seeded the pool from the sourced default on the D-45 "absent does not mean 0" argument.
    # That argument is sound only where a consumer exists: with none, the model ALREADY asserts the
    # burst never happens, and a non-zero seed additionally asserts an antioxidant that is present
    # and never spent. This is the assertion that makes the seed track the wiring instead of the
    # decision that introduced it.
    compiled = compile_scenario(wine_scenario(), oxidative=oxidative)
    seeded = float(compiled.y0[compiled.schema.slice("burst_antioxidant")][0])
    if oxidative == "direct_burst":
        assert seeded == pytest.approx(0.0033, rel=1e-12), "sourced burst_antioxidant_initial"
    else:
        assert seeded == 0.0


@pytest.mark.parametrize("oxidative", ["direct", "cascade"])
def test_dosing_the_pool_into_a_build_that_cannot_draw_it_raises(oxidative):
    # The hops-without-a-bitterness-model precedent: an input the wired model cannot consume is a
    # user error, not a silently-ignored field. Without this the dose would sit in the output
    # unspent, which is exactly the shape of the defect D-147 is fixing.
    with pytest.raises(ValueError, match="burst_antioxidant_gpl"):
        compile_scenario(wine_scenario(burst_gpl=0.005), oxidative=oxidative)


def test_the_dose_is_honoured_where_the_consumer_is_wired():
    compiled = compile_scenario(wine_scenario(burst_gpl=0.005), oxidative="direct_burst")
    assert float(compiled.y0[compiled.schema.slice("burst_antioxidant")][0]) == pytest.approx(0.005)


def test_the_default_build_emits_no_dangling_burst_pool(direct_run):
    # The M0 defect, stated as the thing it forbids: a 2 y default run used to emit
    # burst_antioxidant constant at 3.3e-3 with ptp == 0.0 exactly — a pool nothing could draw,
    # visible in every output column for the life of the wine. Exact equality, not a tolerance:
    # an untouched slot seeded 0 is bit-for-bit 0.
    compiled, traj = direct_run
    pool = traj.y[compiled.schema.slice("burst_antioxidant")][0]
    assert float(np.max(np.abs(pool))) == 0.0


# ------------------------------------------------------------------------------------
# Isolability (prime directive #3) — the default must be untouched.
# ------------------------------------------------------------------------------------


@pytest.mark.parametrize("slot", sorted(DIRECT_PINS))
@pytest.mark.parametrize("years", [1.0, 2.0], ids=["1y", "2y"])
def test_the_direct_set_is_unmoved_by_the_burst_set_existing(direct_run, slot, years):
    # Prime directive #3, and the whole justification for shipping the burst opt-in rather than
    # default: D-140's pins must still be D-140's pins. If this goes red the burst leaked into the
    # default build, and that is a finding for the decision record, not a tolerance to widen.
    compiled, traj = direct_run
    hours = (FERMENT_DAYS + 365.25 * years) * 24.0
    expected = DIRECT_PINS[slot][0 if years == 1.0 else 1]
    assert _at(traj, compiled, slot, hours) == pytest.approx(expected, rel=_PIN_RTOL)


def test_an_empty_burst_pool_reproduces_the_direct_trajectory_exactly():
    # The isolability condition in its strongest form: with the pool explicitly dosed to 0 the
    # burst set is not merely close to the direct set, it is the SAME model — the Process's
    # `burst <= 0` guard returns byte-for-byte zero, so adding it perturbs no RHS value, BDF
    # selects identical steps, and every slot must agree bitwise. An `approx` here would pass on a
    # set that had quietly started contributing.
    #
    # Compare the WHOLE (n_states, n_times) array. This first shipped as `.y[-1]`, which is not the
    # final state vector but the LAST SLOT's time series — `quinone`, identically 0.0 under both
    # sets by `test_quinone_is_identically_zero_under_the_direct_set`. It compared zeros to zeros
    # and would have passed on any divergence whatsoever: the exact defect class this decision
    # spends four sections indicting, committed inside the test that certifies it. Amended at D-147
    # rather than quietly rewritten, because "isolability is checked bitwise" is already in the
    # append-only archive.
    direct = compile_scenario(wine_scenario(), oxidative="direct")
    empty = compile_scenario(wine_scenario(burst_gpl=0.0), oxidative="direct_burst")
    direct_y = _integrate(direct, n=800).y
    empty_y = _integrate(empty, n=800).y
    assert direct_y.shape == empty_y.shape
    assert np.array_equal(direct_y, empty_y)


# ------------------------------------------------------------------------------------
# What wiring it actually does — the D-147 measurement, pinned.
# ------------------------------------------------------------------------------------


@pytest.mark.parametrize("slot", sorted(BURST_PINS))
@pytest.mark.parametrize("years", [1.0, 2.0], ids=["1y", "2y"])
def test_the_burst_set_reproduces_its_measured_trajectory(burst_run, slot, years):
    # The cost of the alternative, made checkable. These are the numbers the D-147 record quotes;
    # a green suite proves nothing about a quoted decimal unless something asserts it.
    compiled, traj = burst_run
    hours = (FERMENT_DAYS + 365.25 * years) * 24.0
    expected = BURST_PINS[slot][0 if years == 1.0 else 1]
    assert _at(traj, compiled, slot, hours) == pytest.approx(expected, rel=_PIN_RTOL)


def test_wiring_the_burst_suppresses_the_oxidative_fates_by_about_a_third(direct_run, burst_run):
    # The headline, as a BAND rather than a decimal: this is a structural consequence of adding a
    # sink to a supply-limited pool (o2 quasi-steady-states at closure_otr / sum(k), so a sink
    # holding ~37% of sum(k) scales every fate by ~0.63), not a tuned magnitude. Stated as a band
    # so it keeps meaning if a sibling rate is legitimately re-sourced.
    d_c, d_t = direct_run
    b_c, b_t = burst_run
    hours = (FERMENT_DAYS + 365.25 * 2.0) * 24.0
    for slot in ("o2", "A420"):
        ratio = _at(b_t, b_c, slot, hours) / _at(d_t, d_c, slot, hours)
        assert 0.55 <= ratio <= 0.70, f"{slot} suppression left the measured band at {ratio:.4f}"
    # SO2 survives LONGER, because it keeps a smaller share of a pool being drained faster.
    assert _at(b_t, b_c, "so2_total", hours) > _at(d_t, d_c, "so2_total", hours)


def test_the_burst_is_not_a_transient_under_a_cork(burst_run):
    """GUARD: forbids re-describing the burst as self-exhausting under the default closure.

    D-133's load-bearing claim is the SHAPE — *"a fast, self-exhausting sink... once exhausted,
    only the D-132 steady rate remains"* — with the magnitudes conceded as order-of-magnitude
    estimates. Under a natural cork that shape does not exist: the closure delivers 2.09 mg/L of
    O2 in two years against a pool sized to absorb 3.3 mg/L, so 77.7% of the pool is still there
    at 2 y and the sink's share of ``sum(k)`` never decays. Its effect GROWS monotonically to a
    plateau instead of fading, which is the opposite of a burst.

    This is named for what it forbids because the alternative is a guard that quietly ratifies the
    prose: asserting "the pool depletes" would pass at any depletion at all, including 22% in two
    years, and would read as confirmation of a claim the measurement refutes.
    """
    compiled, traj = burst_run
    t_age = FERMENT_DAYS * 24.0
    pool0 = float(compiled.y0[compiled.schema.slice("burst_antioxidant")][0])
    left_2y = _at(traj, compiled, "burst_antioxidant", t_age + 2 * 365.25 * 24.0) / pool0
    assert 0.70 <= left_2y <= 0.85, (
        f"{left_2y:.1%} of the burst pool remains at 2 y under a natural cork. D-147 measured "
        "77.7%. If this has become a small number the sink now genuinely self-exhausts, which "
        "would REINSTATE D-133's shape claim — re-measure and re-record before relaxing this."
    )
    # ...and the suppression it causes is still growing at 30 d, not decaying. A real transient
    # would peak early and fade; this is the signature that separates the two.
    o2_1d = _at(traj, compiled, "o2", t_age + 24.0)
    o2_30d = _at(traj, compiled, "o2", t_age + 720.0)
    o2_2y = _at(traj, compiled, "o2", t_age + 2 * 365.25 * 24.0)
    assert o2_1d < o2_30d < o2_2y, "o2 must still be RISING toward its suppressed steady state"


def test_the_pool_cannot_self_exhaust_within_one_saturation_cycle_at_a_real_o2_charge():
    """GUARD: forbids certifying D-133's constraint 2 at an unreachable O2 concentration.

    ``test_burst_pool_exhausts_within_first_saturation_cycle`` (test_aging.py) passes only because
    it sets ``o2 = 0.5 g/L`` — **500 mg/L, ~62x air saturation** — explicitly so *"the burst pool's
    decay is not confounded by O2 itself running out"*. Running out is not a confound; it is the
    binding constraint. Spending the pool takes 3.3 mg/L of O2 through this route alone, and at a
    REAL saturation charge (8 mg/L, the most any wine holds) the burst competes with five siblings
    for it and cannot win enough.

    Air saturation, hermetic, so the charge is a one-shot pool and every mg/L that leaves it went
    through a sink — Ferreira's own protocol, and D-133's own stated condition.
    """
    scen = Scenario(
        name="d147-ferreira-cycle",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.25,
            "anthocyanin_gpl": 0.3,
            "tannin_gpl": 2.0,
            "amino_acids_gpl": 0.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=FERMENT_DAYS + 60.0,
        closure="hermetic",
        interventions=[
            Intervention(day=FERMENT_DAYS - 1.0, action="add_so2", params={"so2_mgl": 60.0}),
            Intervention(day=FERMENT_DAYS, action="begin_aging"),
            Intervention(day=FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": 8.0}),
        ],
    )
    compiled = compile_scenario(scen, oxidative="direct_burst")
    t_end = (FERMENT_DAYS + 60.0) * 24.0
    traj = simulate_scheduled(
        compiled.process_set,
        compiled.param_values,
        compiled.y0.copy(),
        (0.0, t_end),
        events=compiled.events,
        t_eval=np.linspace(0.0, t_end, 6000),
    )
    pool0 = float(compiled.y0[compiled.schema.slice("burst_antioxidant")][0])
    t_dose = FERMENT_DAYS * 24.0
    left_10d = _at(traj, compiled, "burst_antioxidant", t_dose + 240.0) / pool0
    left_60d = _at(traj, compiled, "burst_antioxidant", t_dose + 1440.0) / pool0

    # D-133's constraint 2 wants <= 0.05 here. It measures 0.418.
    assert left_10d > 0.05, (
        "the burst pool self-exhausted within one saturation cycle at a REAL 8 mg/L charge. That "
        "would satisfy D-133's constraint 2 end-to-end for the first time and re-pin "
        "burst_antioxidant_initial — a finding for the record, not a test to update."
    )
    assert left_10d == pytest.approx(0.418, abs=0.01)
    # It PLATEAUS rather than continuing: by 60 d the 8 mg/L charge is spent, so the pool stops
    # where the oxygen ran out. That plateau is why no amount of extra TIME rescues the constraint.
    assert left_60d == pytest.approx(0.392, abs=0.01)
    assert left_10d - left_60d < 0.05, (
        "the pool must be arrested by O2 exhaustion, not still decaying"
    )
