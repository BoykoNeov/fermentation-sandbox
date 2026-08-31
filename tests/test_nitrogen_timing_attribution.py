"""The nitrogen-exhaustion timing miss, ATTRIBUTED — it is the clock, not the nitrogen (D-249).

D-248 §5 recorded that the model strips 90 % of Crépin *et al.* 2017's assimilable nitrogen by
**18.6 h** against her measured ``N_T`` of **28 h** — ~1.5× fast — and §12 left that "open and
unattributed", the owner's call. This file attributes it and the attribution is a **refusal**:
there is no nitrogen-specific timing defect to fix.

**The load-bearing claim is pitch-invariant and time-free.** On the same run, the *whole
fermentation* is faster than Crépin's than the nitrogen channel is: 1.92× against 1.59× at the
shipped pitch, 1.54× against 0.94× at the sourced one. Expressed as a fraction of each run's own
duration — which removes the clock entirely — the model exhausts its nitrogen **later** in the
run than she does (0.225 and 0.308 against her 0.187). A channel that is slower than the run
containing it cannot be what makes the run fast.

**The confound pushes the other way.** The fixture ramps 28 → 18 °C over 14 d where Crépin held
28 °C constant, so the model is fast *despite* running cooler: isothermal at 28 °C it reaches
dryness in 73.0 h rather than 78.2 h.

**Corroboration, and it is only that.** ``tests/benchmarks/test_validation_varela2004.py`` already
pins this engine at **1.6–2.2× fast** against Varela 2004 at the same 28 °C and attributes it
there to a cross-study Coleman-vs-Varela difference rather than a model defect (D-56/D-57) —
Coleman's own reference model gives ~84.5 h for those inputs and the engine matches it. That band
is measured at Varela's **research pitch of 0.04 g/L**, so it is commensurate with Crépin's run
only at matched pitch; the 1.92× at the shipped 0.25 compares runs pitched 6.25× apart and is
**not** a band-membership claim. Matched (0.04, isothermal 28 °C) the run reads 89.2 h = 1.68×,
inside the band. §1's argument does not need this and does not rest on it.

**What the beat DOES identify is a different defect, and it is time-free.** Real yeast strip half
the must while holding a quarter of their final biomass (Crépin's own dry weights: 0.83, 1.39,
3.36, 3.39 g/L at 16, 20, 28, 150 h); this model needs over **half** of its peak to have eaten
the same half. Biomass-against-nitrogen contains no time, so no clock rescaling and no pitch
change can move it — and it barely moves across a 6.25× pitch change. The rate law
``ρ_N ∝ X`` is the wrong functional form early.

**That is NOT D-248 §10 re-opened.** §10 refused a storage-quota *cap* — a brake on uptake once
the cells' nitrogen quota fills. A late brake is the opposite sign and cannot raise the biomass
share at 50 % nitrogen consumed; §10's grounds are untouched and it stays refused. Replacing the
functional form needs a sourced per-cell transport capacity or an explicit storage compartment,
and neither is in hand — so it is **measured and priced, not built**.

**A source note, recorded and not resolved.** Crépin prints her four landmarks twice and
differently: Results gives "16, 20, 40, and 110 h", Methods gives "(16 h, 20 h, 28 h, and 150 h)".
:data:`~tests.test_assimilable_nitrogen_uptake.CREPIN_EXHAUSTION_H` follows Methods, the
operational statement of when the samples were drawn. Nothing here turns on the choice: 40 h makes
the miss larger.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.runtime import simulate_scheduled
from fermentation.scenario import compile_scenario
from fermentation.scenario.schema import TemperaturePoint
from tests.test_assimilable_nitrogen_uptake import CREPIN_EXHAUSTION_H
from tests.test_defined_media import (
    HOUSE_PITCH_GPL,
    SOURCED_PITCH_GPL,
    _assimilable_n_mgl,
    commensurate_scenario,
)
from tests.test_fusel_keto_acid_node import _OTHER_PRECURSOR_CONSUMERS

#: Total sugar at or below which the run counts as dry — the benchmarks' own threshold, so this
#: file's clock is the same clock ``test_validation_varela2004.py`` measures its gap on.
DRYNESS_GPL = 2.0

#: Crépin §Materials and Methods: samples at 16, 20, 28 and 150 h. EF is the one used as the run
#: length here; the Results section's disagreeing 110 h would make the gap *smaller* and the
#: 40 h it pairs with would make the nitrogen miss *larger*, so neither rescues the model.
CREPIN_EF_H = 150.0

#: Crépin's Data Set S1 panel B: dry weight at N½, N¾, N_T and EF. The last is the run's final
#: biomass and is the denominator of the time-free comparison below.
CREPIN_DRY_WEIGHT_GPL = {0.50: 0.83, 0.75: 1.39, 1.00: 3.36}
CREPIN_FINAL_BIOMASS_GPL = 3.39

#: Both pitches are defined in ``tests/test_defined_media.py``, the file that owns the fixture:
#: :data:`~tests.test_defined_media.SOURCED_PITCH_GPL` is what it now runs at, and
#: :data:`~tests.test_defined_media.HOUSE_PITCH_GPL` the 0.25 it ran at through D-252. Every
#: two-pitch comparison below is unchanged by D-253's move — each one passes its pitch in
#: explicitly — but which of the two is the *fixture's* has swapped, and the tests that read the
#: fixture's own default say so.


def _crepin_run(pitch_gpl: float, *, isothermal: bool = False, capacity_ratio: float | None = None):
    """Crépin's own must at a chosen pitch, other precursor consumers off (the D-245 pattern)."""
    scenario = commensurate_scenario("crepin")
    update: dict[str, object] = {"initial": scenario.initial | {"pitch_gpl": pitch_gpl}}
    if isothermal:
        update["temperature_schedule"] = [
            TemperaturePoint(day=0.0, celsius=28.0),
            TemperaturePoint(day=14.0, celsius=28.0),
        ]
    compiled = compile_scenario(scenario.model_copy(update=update))
    for name in _OTHER_PRECURSOR_CONSUMERS:
        compiled.process_set.disable(name)  # KeyErrors on a rename rather than silently no-op
    values = dict(compiled.param_values)
    if capacity_ratio is not None:
        assert "amino_acid_uptake_capacity_ratio" in values
        values["amino_acid_uptake_capacity_ratio"] = capacity_ratio
    traj = simulate_scheduled(
        compiled.process_set,
        values,
        compiled.y0,
        compiled.t_span_h,
        events=compiled.events,
        param_tiers=compiled.parameters.tier_map(),
    )
    assert traj.success, traj.message
    return traj, compiled.schema


class _Course:
    """The three time courses this file reads, and the crossings taken off them."""

    def __init__(self, traj, schema) -> None:
        self.t = traj.t
        self.biomass = traj.y[schema.slice("X"), :][0]
        self.sugar = traj.y[schema.slice("S"), :].sum(axis=0)
        nitrogen = np.array([_assimilable_n_mgl(traj, schema, i) for i in range(traj.y.shape[1])])
        self.consumed = 1.0 - nitrogen / nitrogen[0]

    def _cross(self, series: np.ndarray, level: float, *, rising: bool) -> float:
        hit = np.nonzero(series >= level if rising else series <= level)[0]
        assert hit.size, f"the run never {'reaches' if rising else 'falls to'} {level}"
        i = int(hit[0])
        if i == 0:
            return float(self.t[0])
        lo, hi = float(series[i - 1]), float(series[i])
        before, after = float(self.t[i - 1]), float(self.t[i])
        if lo > hi:  # np.interp needs an increasing x, and the falling series is the sugar one
            lo, hi, before, after = hi, lo, after, before
        return float(np.interp(level, [lo, hi], [before, after]))

    @property
    def hours_to_dryness(self) -> float:
        return self._cross(self.sugar, DRYNESS_GPL, rising=False)

    def hours_to_n_consumed(self, fraction: float) -> float:
        return self._cross(self.consumed, fraction, rising=True)

    def biomass_at(self, hours: float) -> float:
        return float(np.interp(hours, self.t, self.biomass))


@pytest.fixture(scope="module")
def courses():
    """The four runs this file needs, integrated ONCE (the D-245 pattern)."""
    return {
        (pitch, iso): _Course(*_crepin_run(pitch, isothermal=iso))
        for pitch in (HOUSE_PITCH_GPL, SOURCED_PITCH_GPL)
        for iso in (False, True)
    }


def test_the_nitrogen_channel_is_SLOWER_than_the_run_containing_it_at_both_pitches(courses):
    """The whole of D-249's refusal, in one comparison, and it does not depend on the pitch.

    D-248 §5 read the nitrogen as running ~1.5× fast. It does — but the fermentation it sits in
    runs faster still, so the nitrogen cannot be what makes the run fast. Both gaps are measured
    on the same trajectory, so nothing about the must, the medium or the parameters differs
    between them; only the observable does.

    **The four numbers moved at D-253 and the comparison did not.** D-249 measured them at
    ``r = 1``; the shipped capacity is now 2.6, so the nitrogen channel is faster in absolute
    terms at both pitches (1.59 → 1.65 house, 0.94 → 0.96 sourced) while the runs containing it
    barely move (1.92 → 1.93, 1.54 → 1.55). The margin narrowed and the sign did not, which is
    the honest way to state it: the refusal rests on the *ordering* below, and a beat that ever
    closes that ordering has reversed D-249 rather than tightened it.
    """
    for pitch, (clock, nitro) in (
        (HOUSE_PITCH_GPL, (1.93, 1.65)),
        (SOURCED_PITCH_GPL, (1.55, 0.96)),
    ):
        course = courses[(pitch, False)]
        clock_gap = CREPIN_EF_H / course.hours_to_dryness
        nitrogen_gap = CREPIN_EXHAUSTION_H / course.hours_to_n_consumed(0.90)

        assert clock_gap == pytest.approx(clock, abs=0.05), (
            f"pitch {pitch}: the whole-run gap reads {clock_gap:.2f}x, not the {clock}x D-249 "
            "measured — the baseline the attribution is a delta against has moved"
        )
        assert nitrogen_gap == pytest.approx(nitro, abs=0.05), (
            f"pitch {pitch}: the nitrogen gap reads {nitrogen_gap:.2f}x, not the {nitro}x "
            "D-249 measured"
        )
        assert nitrogen_gap < clock_gap, (
            f"pitch {pitch}: the nitrogen channel ({nitrogen_gap:.2f}x) is no longer slower than "
            f"the run containing it ({clock_gap:.2f}x). That REVERSES D-249's attribution and is "
            "a re-decision, not a pin to relax — a nitrogen-timing knob would become arguable "
            "again and D-248 §12 would re-open"
        )


def test_no_transport_CAPACITY_can_make_the_nitrogen_channel_outrun_its_own_run():
    """The refusal's strongest form: the miss is not reachable through the parameter at all.

    D-248 §5 pinned the *extent* as insensitive to ``amino_acid_uptake_capacity_ratio`` across a
    200× sweep and was careful to say the timing was a separate matter that misses. It is also
    unreachable. Taken to **1000×** the shipped value the exhaustion time saturates at ~16.4 h
    and the nitrogen gap tops out at ~1.71×, still under the run's own 1.94×.

    The reason is the rate law: ``ρ_N ∝ X`` makes uptake a slave to the biomass that performs it,
    so however large the capacity, nitrogen cannot be consumed faster than cells accumulate. That
    is the same functional-form defect
    :func:`test_the_model_needs_TWICE_the_biomass_share_to_have_eaten_the_same_nitrogen` measures
    time-free, reached from the other side — and it is why D-249 refuses a capacity knob rather
    than deferring one.

    **This guard is verified non-vacuous by the opposite arm** (D-249 probe 3): scaling
    ``q_sugar_max`` to 0.6 slows the run to 144.1 h while leaving the nitrogen at 17.59 h
    untouched, and the direction guard above — ``test_the_nitrogen_channel_is_SLOWER_than_the_``
    ``run_containing_it_at_both_pitches`` — goes RED. So a future change that fixed the clock
    without touching nitrogen *would* be caught, which is precisely the case in which this
    record's attribution reverses.
    """
    course = _Course(*_crepin_run(HOUSE_PITCH_GPL, capacity_ratio=1000.0))
    clock_gap = CREPIN_EF_H / course.hours_to_dryness
    nitrogen_gap = CREPIN_EXHAUSTION_H / course.hours_to_n_consumed(0.90)

    assert nitrogen_gap == pytest.approx(1.71, abs=0.05), (
        f"at 1000x capacity the nitrogen gap reads {nitrogen_gap:.2f}x, not the 1.71x D-249 "
        "measured — the saturation the refusal rests on has moved"
    )
    assert nitrogen_gap < clock_gap, (
        f"a 1000x transport capacity now makes the nitrogen channel ({nitrogen_gap:.2f}x) faster "
        f"than the run containing it ({clock_gap:.2f}x). The timing miss has become reachable "
        "through the parameter, which is the one thing that would make a capacity knob arguable "
        "and re-open what D-249 refused"
    )


def test_time_free_the_model_exhausts_its_nitrogen_LATER_in_its_run_than_crepin_does(courses):
    """The same claim with the clock divided out, so no rate calibration can be read into it.

    Crépin exhausts at 28 h of a 150 h fermentation — 18.7 % of the way in. Expressed the same
    way the model is at 21.8 % (house pitch) and 30.3 % (sourced). Both are *later*, which is
    the strongest form of "the nitrogen is not what is fast": it survives multiplying either
    run's time axis by any constant.
    """
    measured = CREPIN_EXHAUSTION_H / CREPIN_EF_H
    assert measured == pytest.approx(0.187, abs=0.002)

    for pitch, expected in ((HOUSE_PITCH_GPL, 0.218), (SOURCED_PITCH_GPL, 0.303)):
        course = courses[(pitch, False)]
        share = course.hours_to_n_consumed(0.90) / course.hours_to_dryness
        assert share == pytest.approx(expected, abs=0.01), (
            f"pitch {pitch}: nitrogen exhausts at {share:.3f} of the run, not the {expected} "
            "D-249 measured"
        )
        assert share > measured, (
            f"pitch {pitch}: the model now exhausts its nitrogen at {share:.3f} of its own run "
            f"against Crépin's {measured:.3f} — EARLIER rather than later, which reverses the "
            "time-free half of D-249's attribution"
        )


def test_the_temperature_confound_pushes_the_OTHER_way(courses):
    """The fixture is cooler than Crépin's run, so the model is fast *despite* the confound.

    Crépin held 28 °C throughout; :func:`~tests.test_defined_media.commensurate_scenario` ramps
    28 → 18 °C over its 14 days, because that ramp is the fixture's shared shape and not a claim
    about her protocol. Correcting it makes the run faster, so the fixture's schedule cannot be
    the explanation for the gap D-249 attributes — a finding that survives a confound pointing
    at it is worth more than one that needs the confound removed.
    """
    for pitch in (HOUSE_PITCH_GPL, SOURCED_PITCH_GPL):
        ramped = courses[(pitch, False)].hours_to_dryness
        isothermal = courses[(pitch, True)].hours_to_dryness
        assert isothermal < ramped, (
            f"pitch {pitch}: isothermal 28 °C now reaches dryness in {isothermal:.1f} h against "
            f"the ramp's {ramped:.1f} h — the confound has changed sign and D-249's "
            "'fast despite running cooler' reading no longer holds"
        )


def test_matched_to_varelas_own_research_pitch_the_clock_gap_lands_INSIDE_its_band(courses):
    """Corroboration only — and the pitch match is what makes it a comparison at all.

    ``test_validation_varela2004.py`` pins [1.6, 2.2]× at 28 °C **at 0.04 g/L**, its own
    ``_PITCH_GPL_RESEARCH``. Comparing the shipped fixture's 1.92× to that band would be
    comparing runs pitched 6.25× apart, which is why D-249 does not lead with it. Reconstructed
    at Crépin's stated temperature and the sourced pitch, the gap reads 1.69× — a different must
    and a different nitrogen level landing in the same characterised band, on the same engine.
    """
    course = courses[(SOURCED_PITCH_GPL, True)]
    gap = CREPIN_EF_H / course.hours_to_dryness
    assert gap == pytest.approx(1.69, abs=0.05), (
        f"the matched-pitch isothermal gap reads {gap:.2f}x, not the 1.69x D-253 measured"
    )
    assert 1.6 <= gap <= 2.2, (
        f"the matched gap {gap:.2f}x has left the band test_validation_varela2004.py pins for "
        "this engine at 28 °C; the corroboration is gone and D-249 §1 should be re-read on its "
        "own evidence rather than this one"
    )


def test_the_model_needs_TWICE_the_biomass_share_to_have_eaten_the_same_nitrogen(courses):
    """The residue D-249 identifies, and it contains no time (decision D-249 §3).

    Crépin's dry weights put her yeast at 24.5 % of final biomass with half the must's nitrogen
    already gone, and 41.0 % with three-quarters gone. The model needs 54.3 % and 77.7 %. Because
    both sides are ratios of a biomass to a biomass at a *nitrogen* level, the comparison has no
    time axis: rescaling either run's clock leaves it identical, and the near-equal columns at
    two pitches 6.25× apart show it is not an inoculum artefact either.

    **The frame these numbers are in, which D-249 did not label and D-251 does.** They are read
    through :func:`~tests.test_defined_media._assimilable_n_mgl`, which counts the D-250
    intracellular store as *unconsumed*. Crépin sampled the **medium**. In her frame the same
    shipped run reads **41.2 %** and **62.6 %**, so the gap is 1.68×/1.53× rather than the
    2.22×/1.90× below — see ``tests/test_nitrogen_uptake_shape_frame.py``. The numbers asserted
    here are correct and stay; what D-251 corrects is D-249 §5's *reading* of them, which
    concluded ``ρ_N ∝ X`` is the wrong functional form. In Crépin's own frame the shipped
    capacity knob crosses her landmark inside its declared band, so the form is not shown to be
    wrong — the level was never identifiable on an observable this repo could express before
    D-250 split the frames. Neither record re-opens the storage cap D-248 §10 refused, which is
    the opposite sign and could only lower these numbers late.
    """
    for pitch, expected in (
        (HOUSE_PITCH_GPL, {0.50: 0.543, 0.75: 0.777}),
        (SOURCED_PITCH_GPL, {0.50: 0.513, 0.75: 0.762}),
    ):
        course = courses[(pitch, False)]
        peak = float(course.biomass.max())
        for fraction, pinned in expected.items():
            hours = course.hours_to_n_consumed(fraction)
            share = course.biomass_at(hours) / peak
            measured = CREPIN_DRY_WEIGHT_GPL[fraction] / CREPIN_FINAL_BIOMASS_GPL

            assert share == pytest.approx(pinned, abs=0.02), (
                f"pitch {pitch}, {fraction:.0%} of the nitrogen gone: the model holds "
                f"{share:.3f} of its peak biomass, not the {pinned} D-249 measured"
            )
            assert share > 1.5 * measured, (
                f"pitch {pitch}, {fraction:.0%} of the nitrogen gone: the model now holds "
                f"{share:.3f} of peak against Crépin's {measured:.3f} — the front-loading gap "
                "has closed below 1.5x, which is the one thing that would make a replacement "
                "rate law unnecessary. That is a finding to record, not a bound to relax"
            )


def test_the_fixtures_pitch_is_sourced_and_the_house_default_is_6_25x_larger():
    """A declared-fact guard on the pitch, now that the move D-249 priced has been made (D-253).

    Crépin never states an inoculum — ``grep -c inocul`` over the full PMC text is **0** — so for
    her run the pitch is inherited by lineage: the same Bely, Sablayrolles & Barre 1990 medium out
    of the same lab as Minebois, who states 1 × 10⁶ cells mL⁻¹. For the Minebois fixture the same
    number is *directly* sourced, which is what makes the shared constant defensible rather than a
    second house default wearing a citation.

    **What the move bought and what it cost.** Nitrogen exhaustion goes 17.6 h → 29.9 h against
    her measured 28 h; peak biomass goes +0.8 % → −5.2 % of her measured 3.39 g/L. **The
    Coleman-anchor line of D-249's price list is deleted, not paid** (decision D-252): the
    0.984 → 0.925 was ``test_biomass_now_reaches_the_coleman_yield_the_compile_seam_installs``
    predicting against a hardcoded 0.25 g/L inoculum the run did not have. Against the run's own
    pitch the anchor reads 0.9848, inside the band, and is *structurally insensitive* to the
    inoculum — it drifts 0.0004 across the 6.25×. See ``tests/test_inoculum_and_cell_nitrogen.py``.

    **Minebois was the unpriced half and it does not regress**: her two in-study fusel shares read
    1.009×/0.977× at the house pitch and 1.014×/0.982× at the sourced one, inside the pin either
    way (``test_both_minebois_legs_land_on_her_own_measurement_once_uptake_is_uncoupled``).
    """
    assert abs(SOURCED_PITCH_GPL - 0.04) < 1e-6, (
        f"1e6 cells/mL now converts to {SOURCED_PITCH_GPL} g/L, not the 0.04 D-249 measured "
        "through the D-219 crossing"
    )
    ratio = HOUSE_PITCH_GPL / SOURCED_PITCH_GPL
    assert ratio == pytest.approx(6.25, abs=0.01), (
        "the house default is no longer 6.25x the sourced pitch; the two-pitch comparisons in "
        "this file are all quoted against that separation"
    )
    scenario = commensurate_scenario("crepin")
    assert scenario.initial["pitch_gpl"] == pytest.approx(SOURCED_PITCH_GPL), (
        "the shared Crépin fixture is no longer pitched at the sourced 1e6 cells/mL. If it went "
        "back to the house 0.25 that is a re-decision reversing D-253, not a default to restore "
        "quietly — every number in this file and in tests/test_assimilable_nitrogen_uptake.py is "
        "scored at the sourced pitch"
    )
    assert scenario.initial["pitch_gpl"] != pytest.approx(HOUSE_PITCH_GPL), (
        "the sourced and house pitches have collided, which makes every two-pitch comparison in "
        "this file vacuous rather than merely wrong"
    )
