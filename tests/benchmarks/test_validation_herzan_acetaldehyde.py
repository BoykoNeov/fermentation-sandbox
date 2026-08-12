"""D-188: Herzan et al. 2020's SO2 ladder — the one dataset that grades every stage at once.

Source, read as the paper rather than as this archive's note (the D-187 lesson):

* Herzan, J., Prokes, K., Baron, M., Kumsta, M., Pavlousek, P. & Sochor, J. (2020). "Study of
  carbonyl compounds in white wine production." *Food Science & Nutrition* 8(11):5850-5859.
  doi:10.1002/fsn3.1855, **PMC7684598 — open access, full text read.** Table 1 gives total
  acetaldehyde (the method "is based on the hydrolysis of carbonyl compounds bound to SO2", so
  bound acetaldehyde is counted) for **seven** SO2 regimes, notated (must / tank / bottling) in
  mg/L, measured at the last sampling — day 182 from grape processing, after roughly three
  months in stainless steel and a bottle tail.

**D-108 carried FIVE of the seven rows.** The two it omitted, (0/60/0) and (60/60/0), are the
ones that isolate the tank column at double strength, and they are the reason the failure below
can be stated as a sign error rather than a magnitude error. Re-opening the source found them;
re-reading the note would not have.

**What "maintained" means, and why this file integrates in segments.** The tank column is not a
dose. The paper: *"Samples for measurement were taken once every 14 days. After each
determination of SO2, its volume in the wine was adjusted to the specified value."* So a (x/30/y)
wine was topped back up to 30 mg/L **free** SO2 every fortnight for three months. :func:`_run`
reproduces that — it integrates fortnight by fortnight and bisects for the dose that restores
free SO2 to target in the state it actually reaches. A single dose at the start of the tank phase
is a different treatment, and naming it after the paper's would be the error this file exists to
catch.

**THE HEADLINE IS NOT A RATIO, AND THE RATIO WOULD HAVE READ WELL.** The model's three
sulfited-must variants come out 26.256 / 26.507 / 26.604 — a 1.3 % spread against a published
span of 17.2 to 51.6, a factor of three. Quoting the middle one as "1.53x high" would be
reporting agreement from a model that emits **one constant against three published values**: with
that spread, one of the three ratios is guaranteed to look respectable. What is true is
structural, and it is what this file asserts:

    the model collapses a three-stage treatment onto the must stage alone.

Within sulfited must it is flat (26.256 -> 26.604) where the paper spans 3x. Within unsulfited
must it **declines** 11x (1.013 -> 0.089) where the paper **rises** 10x (2.7 -> 25.9).

**The attribution, measured rather than argued** (:func:`test_removing_the_sulfite_oxygen_
competition_flips_the_sign_but_not_the_size`). There are two candidate causes and they own
different halves of the failure:

1. *The sign* belongs to O2 competition. The model's only post-fermentation acetaldehyde source
   is :class:`OxidativeAcetaldehyde`, and :class:`SulfiteOxidation` draws on the same dissolved-O2
   pool, so dosing SO2 during maturation starves the source. Scale ``k_so2_oxidation`` to zero and
   the pair does flip to the published direction — the D-47 protection mechanism IS present and
   IS working.
2. *The size* belongs to a missing source. That flip is 1.0137 -> 1.0236, **+1.0 %**, against the
   paper's +298 %. Even with the competition entirely removed the model has no post-dryness
   acetaldehyde source at all: both terms of :class:`AcetaldehydeProduction` — the D-27 base
   borrow and the D-48 SO2-induced over-production — are gated on the fermentative flux
   ``X * S/(K+S)``, which is **exactly zero** from dryness onward. Herzan's mechanism for the tank
   column is precisely that regime: *"Acetaldehyde formation is a way to protect the yeast from
   the antiseptic effects of SO2"*, and *"the highest concentrations of acetaldehyde occur in the
   presence of free SO2 and active yeast"* — a young wine on lees, which this model does not have.

**Why nothing is built here.** A post-dryness SO2-driven production term has no magnitude anchor
except Table 1 itself, and calibrating against the table this file checks is a fit wearing a
benchmark's name. The gap is measured, attributed and pinned as ``xfail(strict=True)`` so it
cannot be closed silently — the D-142 idiom.

**The one thing that got BETTER without anyone aiming at it.** D-108 recorded the unsulfited
floor as 0.000 against Herzan's 2.7 and left it deliberately unpatched, because at the time a
sealed bottle admitted no oxygen at all. D-136's steady ingress and D-187's bottling burst gave
it one, and the floor now lands at 0.4668 (screwcap) to 2.2591 (SupremeCorq) with nothing fitted
— the published value sits just above the top of the closure menu rather than infinitely above a
zero. That span, not any single closure's number, is the claim.

**Scope, stated rather than conceded.**

* The model's aging axis has one phase, not two: :class:`ClosureOxygenIngress` runs from
  ``begin_aging``, so the three tank months are metered at the bottle's closure rate. A real
  stainless tank with weekly lees stirring takes up more oxygen than a corked bottle, so the tank
  phase here is oxygen-poor relative to the experiment. It does not rescue the sign — the
  competition probe removes the O2 term's effect entirely and still leaves +1 %.
* The paper's stated components (ferment + three months tank + a bottle tail) do not sum to its
  stated 182 days, so the interior split is not asserted. Traded off at a fixed 182-day total
  (tank 60/90/130 d), the sulfited-must rows move 26.25-26.99 and the unsulfited floor
  0.622-1.211: every claim here survives the ambiguity, and the receipts record the sweep.
* This is a white wine and the scenario seeds **no tannin and no anthocyanin**, so the
  acetaldehyde-bridging sinks (D-80 and siblings) are inert. That is asserted, not assumed —
  it is what makes "SO2 protection has almost nothing to protect" a measurement rather than a
  guess, and it is why this file must not be re-derived on a red wine.

Receipts: ``M:\\claud_projects\\temp\\ferment\\d188-herzan-ladder\\``.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.acidbase import free_so2
from fermentation.core.kinetics.carbon_routing import fermentative_flux_shape
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario.compile import compile_scenario
from fermentation.scenario.schema import Intervention, Scenario, TemperaturePoint

pytestmark = pytest.mark.benchmark

#: Herzan et al. 2020 Table 1, **total** acetaldehyde in mg/L, keyed by the paper's own
#: (must / tank / bottling) SO2 notation in mg/L. Transcribed from the paper's rendered table.
PUBLISHED: dict[tuple[float, float, float], float] = {
    (0.0, 0.0, 0.0): 2.7,
    (0.0, 0.0, 35.0): 6.5,
    (0.0, 30.0, 35.0): 25.9,
    (0.0, 60.0, 0.0): 9.6,
    (60.0, 0.0, 35.0): 17.2,
    (60.0, 30.0, 35.0): 51.6,
    (60.0, 60.0, 0.0): 22.6,
}

#: The three variants sharing a 60 mg/L must dose. Their published spread — 17.2 to 51.6 — is
#: entirely the tank and bottling columns, so it is the cleanest statement of what the model has
#: to reproduce beyond the must stage.
SULFITED_MUST = ((60.0, 0.0, 35.0), (60.0, 30.0, 35.0), (60.0, 60.0, 0.0))

#: "Every 14 days" — the paper's own re-measurement and top-up interval for the tank column.
TOP_UP_INTERVAL_DAYS = 14.0

_FERMENT_DAYS = 21.0
_TANK_DAYS = 90.0
_TOTAL_DAYS = 182.0

#: The closure menu, ascending in oxygen transmission (the ordering closure.yaml ships, with
#: technical cork deliberately below screwcap — D-136).
_ASCENDING_CLOSURES = (
    "hermetic",
    "technical_cork",
    "screwcap",
    "natural_cork",
    "synthetic_nomacorc",
    "synthetic_supremecorq",
)

#: The closure the ladder itself is run under. Any permeable choice tells the same structural
#: story; this one is the menu's mid-point and the one both closure primaries measure directly.
_LADDER_CLOSURE = "natural_cork"


def _scenario(must_so2: float, *, closure: str | None) -> Scenario:
    """The Herzan timeline WITHOUT the tank top-ups and the bottling dose (applied in :func:`_run`).

    A white must at 21 Brix and 18 C. The must dose is the one SO2 addition expressible as a
    plain intervention, because it is a single dose; everything downstream of dryness is a
    maintained level and has to be driven from the state.
    """
    interventions = []
    if must_so2:
        interventions.append(Intervention(day=0.0, action="add_so2", params={"so2_mgl": must_so2}))
    interventions.append(Intervention(day=_FERMENT_DAYS, action="begin_aging"))
    if closure:
        interventions.append(Intervention(day=_FERMENT_DAYS + _TANK_DAYS, action="seal_bottle"))
    return Scenario(
        name="herzan2020",
        medium="wine",
        initial={"brix": 21.0, "yan_mgl": 200.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        duration_days=_TOTAL_DAYS,
        closure=closure,
        interventions=interventions,
    )


def _dose_for_free_so2(y, schema, values, target_mgl: float) -> float:
    """The SO2 dose in mg/L that brings **free** SO2 to ``target_mgl`` in this exact state.

    Bisected on the shipped binding solver rather than assumed: how much total SO2 a target free
    level costs depends on the carbonyl load the wine has reached, which is the very quantity
    under test. A fixed top-up would silently be a different treatment in every variant.
    """
    lo, hi = 0.0, 500.0
    probe = np.asarray(y, dtype=np.float64).copy()
    so2 = schema.slice("so2_total")
    base = float(probe[so2][0])
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        probe[so2] = base + mid / 1000.0
        if float(free_so2(probe, schema, values)) * 1000.0 < target_mgl:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _run(
    must: float,
    tank: float,
    bottling: float,
    *,
    closure: str | None = _LADDER_CLOSURE,
    k_so2_scale: float = 1.0,
) -> dict[str, float]:
    """Integrate one Herzan variant, holding free SO2 at ``tank`` through the tank phase.

    The integration is restarted at each top-up. ``ScheduledEvent.reconfigure`` mutates the
    ``ProcessSet`` in place, so the ``begin_aging`` switch fired in an earlier segment stays in
    force across the restarts; each segment is handed only its own events so nothing fires twice.
    """
    compiled = compile_scenario(_scenario(must, closure=closure))
    schema = compiled.schema
    values = dict(compiled.param_values)
    values["k_so2_oxidation"] *= k_so2_scale

    bottling_day = _FERMENT_DAYS + _TANK_DAYS
    stops: list[float] = []
    if tank:
        day = _FERMENT_DAYS + TOP_UP_INTERVAL_DAYS
        while day < bottling_day - 1e-9:
            stops.append(day)
            day += TOP_UP_INTERVAL_DAYS
    stops.append(bottling_day)

    so2 = schema.slice("so2_total")
    y = compiled.y0.copy()
    t0 = 0.0
    for stop in [*stops, _TOTAL_DAYS]:
        segment_events = [e for e in compiled.events if t0 * 24.0 <= e.time_h < stop * 24.0 - 1e-9]
        traj = simulate_scheduled(
            compiled.process_set,
            values,
            y,
            (t0 * 24.0, stop * 24.0),
            events=segment_events,
            t_eval=np.array([t0 * 24.0, stop * 24.0]),
        )
        assert traj.success, traj.message
        y = traj.y[:, -1].copy()
        t0 = stop
        if tank and stop <= bottling_day + 1e-9:
            y[so2] += _dose_for_free_so2(y, schema, values, tank) / 1000.0
        if bottling and abs(stop - bottling_day) < 1e-9:
            y[so2] += bottling / 1000.0

    return {
        "acetaldehyde": float(y[schema.slice("acetaldehyde")][0]) * 1000.0,
        "so2_total": float(y[so2][0]) * 1000.0,
        "free_so2": float(free_so2(y, schema, values)) * 1000.0,
    }


@pytest.fixture(scope="module")
def ladder() -> dict[tuple[float, float, float], float]:
    """Total acetaldehyde in mg/L at day 182 for each of Herzan's seven variants."""
    return {variant: _run(*variant)["acetaldehyde"] for variant in PUBLISHED}


@pytest.fixture(scope="module")
def floor_by_closure() -> dict[str, float]:
    """The unsulfited (0/0/0) endpoint across the closure menu — D-108's floor, re-measured."""
    return {c: _run(0.0, 0.0, 0.0, closure=c)["acetaldehyde"] for c in _ASCENDING_CLOSURES}


# --------------------------------------------------------------------------------------------
# What the model gets right
# --------------------------------------------------------------------------------------------


def test_must_sulfiting_strands_acetaldehyde(ladder):
    """Sulfiting the must raises finished acetaldehyde — the one column the model has machinery for.

    This is the D-47 protection mechanism working on the stage it was built for: bound
    acetaldehyde is invisible to ADH, so the yeast cannot draw the ferment peak back down.

    **Asserted on Herzan's two MATCHED pairs**, where only the must column differs — grouping all
    sulfited-must variants against all unsulfited ones does not work, and the reason is the
    paper's point: (0/30/35) = 25.9 beats (60/0/35) = 17.2, so the tank column can out-do the
    must column. Any "sulfited must is higher" claim has to hold tank and bottling fixed.

    **And it is the ABSOLUTE lift that is asserted, not the ratio.** Published, the must dose adds
    10.7 and 25.7 mg/L; modelled it adds 25.95 and 26.41. As a ratio the model is wild (x86 and
    x275 against the paper's x2.6 and x2.0), because its unsulfited baseline is far too low — the
    ratio is reporting the denominator. The 26.41-vs-25.7 near-agreement on the second pair is
    the flat-model lottery the module docstring warns about, not corroboration: the model returns
    the same ~26 mg/L for every sulfited-must variant, so one published lift was going to be met.
    """
    for tank, bottling in ((0.0, 35.0), (30.0, 35.0)):
        published_lift = PUBLISHED[(60.0, tank, bottling)] - PUBLISHED[(0.0, tank, bottling)]
        modelled_lift = ladder[(60.0, tank, bottling)] - ladder[(0.0, tank, bottling)]
        assert published_lift > 0.0
        assert modelled_lift > 0.0  # the direction
        assert 0.5 < modelled_lift / published_lift < 3.0  # and the order of magnitude


def test_the_unsulfited_floor_is_no_longer_zero_and_the_closure_sets_it(floor_by_closure):
    """D-108's 0-vs-2.7 floor, re-measured after D-136 and D-187 gave a sealed bottle oxygen.

    D-108 measured 0.000 against Herzan's 2.7 and declined to patch it, on the ground that
    inventing a constant to hit 2.7 bought nothing observable. Nobody patched it; the closure
    oxygen axis closed most of it as a side effect. The claim is the SPAN, not any one closure:
    the published value sits just above the top of the menu.
    """
    # No oxygen, no source. The tolerance is the SOLVER's, not a hedge: acetaldehyde peaks near
    # 42 mg/L mid-ferment and is driven back to zero, so BDF at rtol=1e-6 leaves a residue of
    # order 1e-6 mg/L (measured: -7.9e-7). Pinning a hard zero here would be pinning the
    # integrator's luck, which is the D-108 "hard zero the chemistry cannot make" shape.
    assert floor_by_closure["hermetic"] == pytest.approx(0.0, abs=1e-5)
    permeable = {c: v for c, v in floor_by_closure.items() if c != "hermetic"}
    assert all(v > 0.0 for v in permeable.values())
    assert max(permeable.values()) < PUBLISHED[(0.0, 0.0, 0.0)]  # still under-predicts
    assert max(permeable.values()) > 0.8 * PUBLISHED[(0.0, 0.0, 0.0)]  # but only just
    # The most permeable closure is the one that gets closest — the floor is an oxygen supply
    # question now, which is exactly what D-108 said it was not able to be.
    assert max(permeable, key=lambda c: permeable[c]) == "synthetic_supremecorq"


def test_the_two_dosing_paths_are_the_same_mutation():
    """The must dose goes through the verb; the tank top-ups and bottling dose write the slot.

    Two paths in one experiment, and the finding is precisely that the must column behaves
    differently from the later columns — so the harness must not be the reason. ``_verb_add_so2``
    has to be *exactly* ``so2_total += mg/L / 1000``: no clamp, no companion slot, no parameter
    update. Asserted against the compiled event rather than read off the source, which is the
    D-187 ``seal_bottle`` ≡ ``add_oxygen`` idiom — a verb cannot quietly grow a second effect.
    """
    dose_mgl = 60.0
    compiled = compile_scenario(_scenario(dose_mgl, closure=_LADDER_CLOSURE))
    events = [e for e in compiled.events if e.label.startswith("add_so2@")]
    assert len(events) == 1
    event = events[0]
    assert event.param_update in (None, {}), "a dose must not move a parameter"

    before = compiled.y0.copy()
    assert event.mutate is not None
    after = np.asarray(event.mutate(compiled.schema, before))

    by_hand = before.copy()
    by_hand[compiled.schema.slice("so2_total")] += dose_mgl / 1000.0
    assert np.array_equal(after, by_hand)


def test_the_white_wine_has_no_acetaldehyde_bridging_substrate():
    """Tannin and anthocyanin are zero here, so the bridging sinks are inert.

    Load-bearing for the attribution: "SO2 protection has almost nothing left to protect" is a
    measurement about THIS scenario, not a property of the model. On a red wine the bridging
    routes (D-80 and the tannin/anthocyanin condensations) are live and SO2 protection would have
    real work to do, so this file's conclusion must not be re-derived there.
    """
    compiled = compile_scenario(_scenario(0.0, closure=_LADDER_CLOSURE))
    for name in ("tannin", "anthocyanin", "ethyl_bridge"):
        assert float(compiled.y0[compiled.schema.slice(name)][0]) == 0.0


def test_the_fermentative_flux_is_exactly_zero_once_the_wine_is_dry():
    """Both acetaldehyde PRODUCTION terms switch off at dryness — the missing-source half.

    ``AcetaldehydeProduction`` returns early on ``flux <= 0``, and its D-48 SO2-induced
    over-production term rides the same flux. So from dryness onward the model has no biological
    acetaldehyde source of any kind, whatever the SO2 level — which is the regime Herzan's tank
    column is entirely made of. Pinned here so the diagnosis cannot rot silently: if someone
    later gives the model a post-dryness source, this test is where the reasoning is recorded.
    """
    compiled = compile_scenario(_scenario(60.0, closure=_LADDER_CLOSURE))
    schema, values = compiled.schema, dict(compiled.param_values)
    traj = simulate_scheduled(
        compiled.process_set,
        values,
        compiled.y0,
        (0.0, _TOTAL_DAYS * 24.0),
        events=compiled.events,
        t_eval=np.linspace(0.0, _TOTAL_DAYS * 24.0, 400),
    )
    at_dryness = int(np.searchsorted(traj.t / 24.0, _FERMENT_DAYS))
    for i in range(at_dryness, traj.y.shape[1]):
        assert fermentative_flux_shape(traj.y[:, i], schema, values["K_sugar_uptake"]) == 0.0


def test_removing_the_sulfite_oxygen_competition_flips_the_sign_but_not_the_size(ladder):
    """The attribution, split between the two causes — and neither half is the whole failure.

    ``SulfiteOxidation`` and ``OxidativeAcetaldehyde`` draw on the same dissolved-O2 pool, so
    maturation SO2 starves the model's only post-ferment acetaldehyde source. Remove that
    competition (``k_so2_oxidation`` -> 0) and the published DIRECTION appears: the D-47
    protection is real and does raise acetaldehyde. But it appears at about +1 %, against the
    paper's +298 %. **The competition owns the sign; the missing post-dryness source owns the
    size**, and a fix aimed at only one of them will not move this benchmark.

    **The 5 % bound below is a bound on TODAY'S model and is EXPECTED to break** when a
    post-dryness source lands — unlike the ``xfail`` block further down, which asserts published
    directions only and never a published magnitude. It is here because it is what makes "the
    size belongs elsewhere" falsifiable rather than rhetorical; whoever closes the gap should
    delete it, not widen it.
    """
    with_competition = (ladder[(0.0, 0.0, 35.0)], ladder[(0.0, 30.0, 35.0)])
    assert with_competition[1] < with_competition[0]  # the shipped model runs the wrong way

    without = tuple(_run(0.0, tank, 35.0, k_so2_scale=0.0)["acetaldehyde"] for tank in (0.0, 30.0))
    assert without[1] > without[0]  # sign restored: protection alone runs the right way
    modelled_lift = without[1] / without[0] - 1.0
    published_lift = PUBLISHED[(0.0, 30.0, 35.0)] / PUBLISHED[(0.0, 0.0, 35.0)] - 1.0
    assert modelled_lift < 0.05  # ~1 %
    assert published_lift > 2.5  # ~298 %


# --------------------------------------------------------------------------------------------
# What the model gets wrong — pinned so it cannot be closed silently (the D-142 idiom)
# --------------------------------------------------------------------------------------------
#
# Every xfail below asserts the PUBLISHED DIRECTION and never a published magnitude. A magnitude
# here would become a fit target for whoever closes the gap, and the gap's only magnitude anchor
# is the very table this file checks against.


@pytest.mark.xfail(
    strict=True,
    reason="D-188: the model has no post-dryness acetaldehyde source, and maturation SO2 "
    "starves the oxidative one by competing for the same O2 — so it runs 0.303 -> 0.096 "
    "where Herzan runs 6.5 -> 25.9",
)
def test_maturation_sulfiting_raises_total_acetaldehyde(ladder):
    """Herzan's central result: holding free SO2 at 30 mg/L through the tank quadruples it."""
    assert PUBLISHED[(0.0, 30.0, 35.0)] > PUBLISHED[(0.0, 0.0, 35.0)]
    assert ladder[(0.0, 30.0, 35.0)] > ladder[(0.0, 0.0, 35.0)]


@pytest.mark.xfail(
    strict=True,
    reason="D-188: same mechanism at the bottling stage — the model runs 1.013 -> 0.303 where "
    "Herzan runs 2.7 -> 6.5",
)
def test_bottling_sulfiting_raises_total_acetaldehyde(ladder):
    """The smallest of the paper's three sulfiting contrasts, and the model still inverts it."""
    assert PUBLISHED[(0.0, 0.0, 35.0)] > PUBLISHED[(0.0, 0.0, 0.0)]
    assert ladder[(0.0, 0.0, 35.0)] > ladder[(0.0, 0.0, 0.0)]


@pytest.mark.xfail(
    strict=True,
    reason="D-188: the three sulfited-must variants come out within 1.3 % of each other "
    "(26.256 / 26.507 / 26.604) because the must dose alone sets the endpoint — the model "
    "collapses a three-stage treatment onto one stage",
)
def test_the_later_stages_separate_the_sulfited_must_variants(ladder):
    """Published, these three span a factor of three on the tank and bottling columns alone.

    Asserted as a spread rather than per-variant, because per-variant ratios against a model that
    emits one constant are a lottery — one of the three is bound to look good.
    """
    modelled = [ladder[v] for v in SULFITED_MUST]
    published = [PUBLISHED[v] for v in SULFITED_MUST]
    assert max(published) / min(published) > 2.5
    assert max(modelled) / min(modelled) > 2.0
