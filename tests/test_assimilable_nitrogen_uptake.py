"""Assimilable-nitrogen uptake un-coupled from growth demand (decision D-248).

**What this closes.** D-246 §5 measured that the model leaves **40.8 %** of Crépin *et al.*
2017's assimilable nitrogen standing where its Data Set S1 measures **0.2 %**, and D-247 §7 left
that finding untouched and the repair the owner's call. The cause was arithmetic rather than
parametric: before this beat the only route from the speciated amino-acid pools into biomass
nitrogen was the D-32 swap at ``ψ·gate·f_N·base_dx``, strictly below growth's own ``f_N·base_dx``
draw, so the ammonium slot could only fall; when it reached zero growth's ``N/(K_n + N)`` Monod
shut growth off, the swap — proportional to ``base_dx`` — stopped with it, and the pools froze.
:class:`~fermentation.core.kinetics.amino_acids.AssimilableNitrogenUptake` is the flux that is
**not** proportional to the growth rate, so uptake outruns demand and the must is consumed.

**The load-bearing anchor is INTERNAL and was already shipped.**
:func:`~fermentation.scenario.compile._apply_nitrogen_dependent_yield` overrides
``biomass_N_fraction`` to ``1/Y_X/N(N_init)`` from Coleman, Fish & Block (2007) expressly so that
biomass comes out at ``Y_X/N × N_init`` — an identity that holds **only under complete
consumption**. Measured on Crépin's must, the speciated run reached **61.6 %** of that
prediction. Crépin's 0.2 % is the *independent* check and deliberately not the target; the
parameter is not fitted to either, and the tests below pin that by sweeping it.

**Two things this file is careful to keep separate.**

* The **extent** saturates. Residual, biomass and every fusel share are unmoved from
  ``r = 0.25`` to ``r = 50`` — a 200× sweep — so the shipped ``r = 1.0`` is a *bound* ("transport
  is not the bottleneck"), not a level anyone chose. The residual that remains is set by
  ``K_amino_acids``'s own asymptote, not by ``r``.
* The **timing** does not, and it misses. The model strips 90 % of the assimilable nitrogen by
  ~18.6 h against Crépin's measured N_T of 28 h. Pinned below so that "insensitive across 200×"
  can never be read as "the time course was checked".

**The carbon-park is the design's crux, and it has its own test.** The D-32 swap may refund its
drawn carbon to ``S`` only because its rate is proportional to growth's draw, which bounds the
refund. A flux that is not so proportional has no such bound and would create hexose at zero
growth rate. :func:`test_uptake_creates_no_sugar_at_a_state_where_growth_is_stopped` drives
exactly that state.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.kinetics import AssimilableNitrogenUptake, GrowthNitrogenLimited
from fermentation.core.kinetics.amino_acid_pools import AMINO_ACID_SPECS, ASSIMILABLE_SPECS
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS
from fermentation.core.kinetics.growth import biomass_growth_rate
from fermentation.core.media import get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.validation import assert_conserved, total_carbon, total_nitrogen
from tests.test_defined_media import (
    CREPIN_MEASURED_RESIDUAL_FRACTION,
    _assimilable_n_mgl,
    _run,
)
from tests.test_fusel_keto_acid_node import de_novo_share_of

UPTAKE = AssimilableNitrogenUptake.name
GROWTH = GrowthNitrogenLimited.name
SKELETON = "amino_acid_skeleton_carbon"
#: The D-250 intracellular store — where the nitrogen half of the same transfer goes.
STORED = "stored_nitrogen"

#: Crépin's own N_T — the time her Data Set S1 records the must's assimilable nitrogen exhausted
#: at, 28 h against an end of fermentation at 150 h.
CREPIN_EXHAUSTION_H = 28.0


@pytest.fixture(scope="module")
def full_params():
    base = default_data_dir()
    return load_parameters(
        base / "wine_generic.yaml",
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "keto_acids.yaml",
        base / "hydrogen_sulfide.yaml",
        base / "aging.yaml",
        base / "thermal.yaml",
    ).resolve()


def _y0(
    schema: StateSchema,
    *,
    x: float = 1.0,
    s: float = 200.0,
    n: float = 0.2,
    aa: float = 0.4,
    generic: float = 0.4,
    t_k: float = 293.15,
) -> FloatArray:
    return schema.pack(
        {
            "X": x,
            "S": [s],
            "E": 40.0,
            "N": n,
            "T": t_k,
            "CO2": 5.0,
            "amino_acids": aa,
            "amino_acids_generic": generic,
        }
    )


def _uptake_only(schema: StateSchema) -> ProcessSet:
    """The uptake Process alone, under the ``touches`` contract."""
    return ProcessSet(schema, [AssimilableNitrogenUptake()], strict=True)


# -- 1. The rate law is a CAPACITY, and that is the whole of the un-coupling ---------------


def test_the_uptake_rate_does_not_read_the_growth_rate(full_params):
    """Two states with the SAME cells and pools but very different growth rates draw alike.

    This is the claim the Process exists to make, tested as a claim rather than by reading the
    source: sugar is the only thing that differs between the two states below, so
    :func:`~fermentation.core.kinetics.growth.biomass_growth_rate` differs by an order of
    magnitude while the nitrogen drawn from the pools is identical. A rate that read ``base_dx``
    — which is what the D-32 swap reads, and what the pre-D-248 model had as its ONLY route —
    could not pass this.
    """
    schema = wine_schema()
    ps = _uptake_only(schema)
    fast = _y0(schema, s=200.0)
    slow = _y0(schema, s=0.5)  # deep into the sugar Monod, so growth nearly stops

    mu_fast = biomass_growth_rate(fast, schema, full_params)
    mu_slow = biomass_growth_rate(slow, schema, full_params)
    assert mu_fast > 10.0 * mu_slow > 0.0, (
        f"the two states no longer differ in growth rate ({mu_fast:.3e} vs {mu_slow:.3e}), so "
        "this test cannot distinguish a coupled rate from an un-coupled one"
    )

    d_fast = ps.total_derivatives(0.0, fast, full_params)
    d_slow = ps.total_derivatives(0.0, slow, full_params)
    # Read at the STORE since D-250: the nitrogen goes intracellular, not into the must's
    # ammonium. The claim under test is unchanged — it is about the rate, not the destination.
    n_fast = float(d_fast[schema.slice(STORED)][0])
    n_slow = float(d_slow[schema.slice(STORED)][0])
    assert n_fast > 0.0
    assert n_fast == pytest.approx(n_slow, rel=1e-12), (
        f"uptake drew {n_fast:.6e} g N/L/h at the fast-growing state and {n_slow:.6e} at the "
        "slow one — it has acquired a dependence on the instantaneous growth rate, which is "
        "exactly the coupling D-248 removed"
    )


def test_uptake_runs_at_a_state_where_growth_is_entirely_stopped(full_params):
    """The pre-D-248 deadlock, driven: no ammonium ⇒ no growth ⇒ (before) no assimilation.

    ``N = 0`` puts growth's Monod term at exactly zero, so ``base_dx`` is 0 and the D-32 swap
    contributes nothing. That was the trap: the swap could only refund ammonium *while growth
    ran*, and growth could only run *while ammonium lasted*, so the amino-acid pools froze with
    40.8 % of Crépin's nitrogen in them. Uptake must still run here or nothing is repaired.
    """
    schema = wine_schema()
    y = _y0(schema, n=0.0)
    assert biomass_growth_rate(y, schema, full_params) == 0.0

    d = _uptake_only(schema).total_derivatives(0.0, y, full_params)
    assert float(d[schema.slice(STORED)][0]) > 0.0, (
        "uptake contributes nothing at a stopped-growth state, so the assimilable nitrogen the "
        "cells can reach can still only fall and the deadlock D-248 names is back"
    )
    # ...and it must not do it through `N` any more (D-250): the must's ammonium is where
    # D-248 parked the surplus, which is what put nitrogen already inside the cell back into
    # the acid-base balance.
    assert float(d[schema.slice("N")][0]) == 0.0
    for spec in ASSIMILABLE_SPECS:
        assert float(d[schema.slice(spec.pool)][0]) < 0.0


# -- 2. The carbon-park, which is why this is not simply the swap with a different rate ----


def test_uptake_creates_no_sugar_at_a_state_where_growth_is_stopped(full_params):
    """The design's crux: an un-coupled flux must NOT refund its carbon to sugar.

    The D-32 swap's no-hexose guarantee is not a property of amino-acid chemistry — it is a
    property of the swap's *rate*, which is proportional to growth's own draw, bounding the
    refund at ``ψ·gate·(aa C:N)/(biomass C:N)`` of the carbon growth removed. Uptake has no such
    bound by construction, so at ``base_dx = 0`` a sugar refund would be pure gluconeogenesis —
    and no clamp fixes it without a C⁰ kink the stiff BDF solver catches on. The skeleton is
    parked in a carbon-only pool instead, and this asserts the parking rather than trusting it.
    """
    schema = wine_schema()
    y = _y0(schema, n=0.0)  # growth stopped, pools full
    d = _uptake_only(schema).total_derivatives(0.0, y, full_params)

    assert float(d[schema.slice("S")].sum()) == 0.0, (
        "uptake now writes the sugar slot. At this state growth draws NO carbon, so any positive "
        "write creates hexose outright and any negative one charges sugar for a skeleton the "
        "cells supplied themselves — the parking is not a stylistic choice (D-248)"
    )
    parked = float(d[schema.slice(SKELETON)][0])
    assert parked > 0.0, "the drawn skeleton went nowhere — carbon is being destroyed"

    # And the parked carbon is EXACTLY the carbon the debited mass carried, which is what makes
    # total_carbon close through the transfer rather than approximately close.
    from fermentation.core.chemistry import carbon_mass_fraction

    expected = sum(
        -float(d[schema.slice(spec.pool)][0]) * carbon_mass_fraction(spec.species)
        for spec in ASSIMILABLE_SPECS
    )
    assert parked == pytest.approx(expected, rel=1e-12)


def test_the_skeleton_park_is_on_the_carbon_ledger_and_off_the_nitrogen_one():
    """Weighted 1.0 on ``total_carbon`` (elemental carbon), absent from ``total_nitrogen``.

    The pool holds a blend of arginine and glutamine skeletons that no single molecule
    represents, so — unlike ``debris``/glucan — it is held as elemental carbon and weighted 1.0,
    the ``N``-slot idiom. Its nitrogen went to ``N``, so it must carry none of its own or the
    transfer would book the same nitrogen twice.
    """
    schema = wine_schema()
    probe = schema.zeros()
    probe[schema.slice(SKELETON)] = 1.0
    assert total_carbon(schema, biomass_carbon_fraction=0.45)(probe) == pytest.approx(1.0)
    assert total_nitrogen(schema, biomass_nitrogen_fraction=0.1)(probe) == pytest.approx(0.0)


def test_both_ledgers_close_on_a_full_crepin_run():
    """Carbon and nitrogen both close with every Process live and uptake doing real work.

    The transfer is aa → ``N`` (nitrogen) + skeleton (carbon), so neither ledger has any slack:
    a wrong carbon fraction on the park or a nitrogen double-book would show here.
    """
    traj, schema, params, cs = _run("crepin")
    parked = float(traj.y[schema.slice(SKELETON), -1][0])
    assert parked > 0.05, (
        f"only {parked:.4g} g C/L reached the skeleton park, so this run is not exercising "
        "uptake and the closure below proves nothing"
    )
    f_n = params["biomass_N_fraction"]
    f_c = params["biomass_C_fraction"]
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), rtol=1e-6)
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), rtol=1e-6)


# -- 3. Isolability (prime directive #3) ----------------------------------------------------


def test_an_empty_pool_is_an_exact_no_op(full_params):
    """The undosed guarantee: gate exactly 0 ⇒ every derivative exactly 0, not merely small."""
    schema = wine_schema()
    y = _y0(schema, aa=0.0, generic=0.0)
    d = _uptake_only(schema).total_derivatives(0.0, y, full_params)
    assert not np.any(d), f"an undosed run is no longer byte-for-byte the core: {d[np.nonzero(d)]}"


def test_the_compile_seam_disables_uptake_on_an_undosed_wine():
    """Wired into the wine medium, disabled when no amino acids are dosed (the D-32 pattern)."""
    from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

    def _scenario(**initial: float) -> Scenario:
        return Scenario(
            name="d248-gate",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25} | initial,
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            interventions=[],
            duration_days=2.0,
        )

    undosed = compile_scenario(_scenario())
    assert UPTAKE in undosed.process_set
    assert not undosed.process_set.is_enabled(UPTAKE)

    dosed = compile_scenario(_scenario(amino_acids_gpl=1.0))
    assert dosed.process_set.is_enabled(UPTAKE)


def test_uptake_is_not_a_rate_modifier_target():
    """Deliberately unscaled by the growth Arrhenius and the carrying cap (decision D-248).

    The swap is a target of both so its refunds track growth's *realised* draw — that is the
    D-32 correctness crux. This Process has no draw to track, and naming it would re-import the
    coupling it exists to remove. The cost is that uptake carries no temperature dependence,
    which the module docstring records as a named simplification rather than an oversight.
    """
    medium = get_medium("wine")
    modifiers = {f().name: f() for f in medium.modifier_factories}
    for name, modifier in modifiers.items():
        assert UPTAKE not in modifier.modifies, (
            f"{name} now scales {UPTAKE}. If that is deliberate the un-coupling claim needs "
            "re-deriving — a carrying cap in particular would stop uptake exactly where the "
            "residual-N regime makes it matter most"
        )
    assert "amino_acid_assimilation" in modifiers["arrhenius_growth"].modifies, (
        "the swap is no longer scaled by the growth Arrhenius, which is the D-32 crux this test "
        "is contrasting against — the contrast is vacuous if that has moved"
    )


def test_the_process_is_speculative_and_declares_what_it_touches():
    """Tier and contract, so the ``strict=True`` runs above mean something."""
    assert AssimilableNitrogenUptake.tier is Tier.SPECULATIVE
    assert set(AssimilableNitrogenUptake.touches) == {
        "amino_acids",
        "amino_acids_generic",
        STORED,
        SKELETON,
    }
    # `N` is NOT touched since D-250, and that is the repair, not a tidy-up: the surplus is
    # inside the cell, and `N` is the slot the acid-base charge balance reads (D-209/D-210).
    assert "N" not in AssimilableNitrogenUptake.touches
    # The six precursor pools are NOT touched, which is what keeps the fusel result attributable
    # to the biomass denominator alone (D-248) and keeps D-100's cross-subsystem starvation shut
    # against the Ehrlich re-route.
    precursors = {s.pool for s in AMINO_ACID_SPECS} - {s.pool for s in ASSIMILABLE_SPECS}
    assert precursors.isdisjoint(AssimilableNitrogenUptake.touches)
    assert "S" not in AssimilableNitrogenUptake.touches
    assert "X" not in AssimilableNitrogenUptake.touches


# -- 4. The measurement: extent saturates, timing does not ----------------------------------


@pytest.fixture(scope="module")
def sweep():
    """Crépin's must at five uptake capacities — computed ONCE (the D-245 pattern)."""
    out = {}
    for ratio in (0.0, 0.25, 1.0, 10.0, 50.0):
        traj, schema, params, cs = _run("crepin", scale={"amino_acid_uptake_capacity_ratio": ratio})
        out[ratio] = (traj, schema, params, cs)
    return out


def _residual(entry) -> float:
    traj, schema, _, _ = entry
    return _assimilable_n_mgl(traj, schema, -1) / _assimilable_n_mgl(traj, schema, 0)


def test_the_must_is_now_consumed_where_it_used_to_freeze(sweep):
    """40.8 % standing → 0.62 %, against Crépin's measured 0.2 % (decision D-248).

    ``r = 0`` reproduces the pre-D-248 model exactly (the Process contributes nothing), so the
    two rows below differ in one number and nothing else. The claim and the margin are on
    separate lines, D-245's own lesson: *the must is now consumed* is what a future beat must not
    lose; the exact 0.62 % is a pin that may legitimately move.
    """
    before = _residual(sweep[0.0])
    after = _residual(sweep[1.0])

    assert before == pytest.approx(0.408, abs=0.01), (
        f"the pre-D-248 residual reads {before:.4f}, not the 40.8 % D-246 §5 measured — the "
        "baseline this whole beat is a delta against has moved"
    )
    assert after < 0.05 * before, (
        f"uptake leaves {after:.4%} against the {before:.4%} it replaced; the repair has stopped "
        "working, which is a finding and not a reason to relax this bound"
    )
    assert 0.004 <= after <= 0.010, f"the residual left [0.4 %, 1.0 %]: {after:.5f}"
    # Still above Crépin's own, and the direction is recorded rather than smoothed over.
    assert after > CREPIN_MEASURED_RESIDUAL_FRACTION, (
        "the model now consumes MORE completely than Crépin measures, which would make the "
        "asymptote below a different claim than the one D-248 records"
    )


def test_the_residual_is_set_by_the_availability_GATE_not_by_the_new_parameter(sweep):
    """The single most load-bearing line in the beat: the outcome is not a function of ``r``.

    Across a **200×** sweep of ``amino_acid_uptake_capacity_ratio`` the residual does not move in
    the third decimal. That is what makes the shipped 1.0 a *bound* — "transport is not the
    bottleneck" — rather than a level fitted to Crépin's 0.2 %. What sets the residual instead is
    ``K_amino_acids``'s own asymptote: the shared D-100 depletion gate → 0 as the pool empties, so
    the last of the nitrogen is drawn ever more slowly whatever the capacity.
    """
    residuals = {r: _residual(sweep[r]) for r in (0.25, 1.0, 10.0, 50.0)}
    lo, hi = min(residuals.values()), max(residuals.values())
    assert hi - lo < 5e-4, (
        f"the residual now spans {lo:.5f}–{hi:.5f} across a 200x capacity sweep {residuals}; it "
        "has become a function of the new parameter, and the 'bound not midpoint' reading of the "
        "shipped value no longer holds"
    )


def test_biomass_now_reaches_the_coleman_yield_the_compile_seam_installs(sweep):
    """The INTERNAL anchor, and the one that is unfitted (decision D-248).

    ``_apply_nitrogen_dependent_yield`` sets ``biomass_N_fraction = 1/Y_X/N(N_init)`` from
    Coleman's regression expressly so peak biomass lands at ``Y_X/N × N_init``. That identity
    holds only if the nitrogen is all consumed, and before this beat the speciated path reached
    **61.6 %** of it — the model silently violating a regression it compiles in. Nothing here is
    fitted to Coleman: ``f_N`` is the seam's own override and the target is arithmetic from it.

    **The inoculum is read off the run, not written in (decision D-252).** This line carried
    ``0.25`` as a literal — correct at the pitch this fixture used to run, and brittle at any
    other, because at pitch 0.04 it subtracted an inoculum the run did not have and read 0.925
    instead of 0.985. D-249 §3 published that 0.925 as the cost of moving the fixture onto its
    only sourced pitch; it was the literal, and ``tests/test_inoculum_and_cell_nitrogen.py`` pins
    both readings.

    **D-253 then made that move, and the anchor is what proves it was free.** The fixture now
    pitches the sourced 0.04 g/L and the ``r = 1`` arm reads **0.9848** where it read 0.9844 —
    a 0.0004 drift across a 6.25× inoculum change, which is what "structurally insensitive"
    means. The ``r = 0`` arm is *not* insensitive and moves 0.616 → **0.591**: that arm is the
    pre-D-248 frozen path, where growth stops with the nitrogen still standing, so its peak is
    set by how much biomass was pitched rather than by how much nitrogen was consumed. The two
    arms moving differently is the finding, not noise — it is the same separation D-248 rests on.
    """
    for ratio, expected, tol in ((0.0, 0.591, 0.02), (1.0, 0.985, 0.02)):
        traj, schema, params, _ = sweep[ratio]
        initial = _assimilable_n_mgl(traj, schema, 0) / 1000.0
        x0 = float(traj.y[schema.slice("X"), :][0][0])
        predicted = x0 + initial / params["biomass_N_fraction"]
        peak = float(traj.y[schema.slice("X"), :][0].max())
        assert peak / predicted == pytest.approx(expected, abs=tol), (
            f"at r={ratio} peak biomass is {peak / predicted:.3f}x Coleman's own Y_X/N x N_init "
            f"({peak:.3f} vs {predicted:.3f}), not the {expected} D-248 measured"
        )


def test_the_exhaustion_TIME_now_LANDS_on_crepins_and_it_was_the_INOCULUM(sweep):
    """The extent saturates; the timing does not — and the timing has stopped missing (D-253).

    Crépin's Data Set S1 puts the must's assimilable nitrogen exhausted at N_T = 28 h, long
    before its end of fermentation at 150 h. D-248 §5 measured this model at **~18.6 h**, ~1.5×
    fast, and D-249 refused to repair it through ``amino_acid_uptake_capacity_ratio`` on the
    ground that the miss is unreachable that way — a sweep to 1000× saturates at ~16.4 h and
    never overtakes the run containing it. That refusal stands and is untouched.

    **What closed it was the fixture's inoculum, which was never Crépin's.** At the sourced pitch
    the same run reaches 90 % consumed at **~30.4 h** on this file's grid idiom, against her 28 h
    — 8.5 % slow rather than 50 % fast, and now on the *other* side of her measurement. Read on
    the extracellular quantity she actually sampled it is nearer still: 28.6 h, a 2 % miss
    (``tests/test_nitrogen_stored_intracellularly.py``). The attribution is measured at both
    pitches in ``tests/test_nitrogen_timing_attribution.py`` (nitrogen gap 1.59× at the house
    pitch, 0.94× at the sourced one) and is not re-derived here.

    **The purpose of this test is unchanged**, and it is why the band below is two-sided rather
    than a pass mark: "insensitive across 200×" is true of the *extent* and is the reason the
    parameter is not fitted, and it must never be read as "the time course was checked". The
    time course is a live observable of that same knob — the second half of this test drives it
    — while the residual it is measured beside is flat to the fifth decimal.
    """
    traj, schema, _, _ = sweep[1.0]
    series = np.array([_assimilable_n_mgl(traj, schema, i) for i in range(traj.y.shape[1])])
    below = np.nonzero(series <= 0.10 * series[0])[0]
    assert below.size, "the must is never 90 % consumed at the shipped capacity"
    t90 = float(traj.t[below[0]])

    assert 28.0 <= t90 <= 33.0, (
        f"the 90 %-consumed time reads {t90:.2f} h against Crépin's {CREPIN_EXHAUSTION_H} h. "
        "Below 28 h the model is FAST again, which is the pre-D-253 reading and means the "
        "fixture's pitch went back to the house 0.25; above 33 h the agreement D-253 bought "
        "has been spent somewhere else. Either way a re-decision, not a pin to widen"
    )
    assert t90 > CREPIN_EXHAUSTION_H, (
        f"the model exhausts the must at {t90:.2f} h, EARLIER than Crépin's "
        f"{CREPIN_EXHAUSTION_H} h. D-253 records the residual miss as slow-side; a flip to the "
        "fast side is the direction D-248 recorded returning"
    )

    # Non-vacuity, and the whole reason this test sits next to the residual one: the knob the
    # extent cannot see moves the timing by 10 h across the same sweep.
    times = []
    for ratio in sorted(sweep):
        if ratio == 0.0:
            continue  # the pre-D-248 path never reaches 90 % consumed at all
        arm_traj, arm_schema, _, _ = sweep[ratio]
        arm = np.array(
            [_assimilable_n_mgl(arm_traj, arm_schema, i) for i in range(arm_traj.y.shape[1])]
        )
        hit = np.nonzero(arm <= 0.10 * arm[0])[0]
        assert hit.size, f"r={ratio} never reaches 90 % consumed"
        times.append(float(arm_traj.t[hit[0]]))

    assert all(b <= a for a, b in zip(times, times[1:], strict=False)), (
        f"the 90 %-consumed time is no longer non-increasing in capacity: {times}"
    )
    assert times[0] - times[-1] > 5.0, (
        f"the timing moved only {times[0] - times[-1]:.1f} h across the sweep {times}. If it has "
        "gone flat too, then the extent's insensitivity and the timing's are the same fact and "
        "this test no longer distinguishes them"
    )


def test_uptake_touches_no_precursor_pool_in_a_driven_run(sweep):
    """Attributability, driven rather than declared (decision D-248).

    The ``touches`` contract already forbids it, but what the beat's fusel result rests on is the
    stronger claim that uptake moves those pools through **one** channel only — the biomass
    denominator. So the six precursor pools are compared between ``r = 0`` and ``r = 1`` at
    *pitch*, where nothing downstream has run yet: any difference there would be a direct draw.
    """
    for spec in AMINO_ACID_SPECS:
        if spec in ASSIMILABLE_SPECS:
            continue
        before = float(sweep[0.0][0].y[sweep[0.0][1].slice(spec.pool), 0][0])
        after = float(sweep[1.0][0].y[sweep[1.0][1].slice(spec.pool), 0][0])
        assert before == pytest.approx(after, rel=1e-12), (
            f"{spec.pool} is seeded differently under uptake, so the two columns are not one "
            "knob apart and no fusel move below can be attributed"
        )


def test_every_fusel_de_novo_share_rises_and_none_falls(sweep):
    """The one channel, measured on all five alcohols (decision D-248).

    Uptake draws no precursor, so the only route to a fusel's de-novo share is the denominator:
    more nitrogen consumed ⇒ more biomass ⇒ more de-novo alcohol off the sugar route, while the
    precursor-derived amount is capped by pools that were already exhausted. Every share must
    therefore rise, and the largest rise must be propanol's — the alcohol whose precursor supply
    is largest relative to its own carbon draw, i.e. the one the denominator moves most.
    """
    deltas = {}
    before_traj, before_schema, before_params, _ = sweep[0.0]
    after_traj, after_schema, after_params, _ = sweep[1.0]
    for spec in FUSEL_SPECS:
        a = de_novo_share_of(before_traj, before_schema, before_params, spec)
        b = de_novo_share_of(after_traj, after_schema, after_params, spec)
        deltas[spec.pool] = b - a
    assert all(v > 0.0 for v in deltas.values()), f"a de-novo share fell under uptake: {deltas}"
    assert max(deltas, key=lambda k: deltas[k]) == "propanol", (
        f"propanol is no longer the alcohol the biomass denominator moves most: {deltas}"
    )
