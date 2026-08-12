"""Residual copper: fining leaves the metal in the wine, and the wine browns faster (D-191).

Copper fining does not remove copper. Understanding Wine Chemistry 2nd ed. §26.2.4.1 measured
white wines dosed to 1 mg/L Cu with equimolar H₂S keeping **>95 % of their copper** through
*filtering* or *5 days settling plus racking*, and concluded copper additions "are not
necessarily a 'fining' operation, as the copper remains in solution, albeit in a different
form"; Ch. 24 names the settles-with-the-lees account as an incorrect assumption of older
textbooks. Since the ``copper`` slot is :class:`PhenolicBrowning`'s mean-centred multiplier
input (D-134), crediting the dose is what finally connects the model's two coppers — the
disconnect D-149 recorded and ``add_copper``'s own docstring called unfixable for want of a
source.

**What is pinned here.** The retained share against its parameter, both band edges *with their
different status*, the consequence at a named site (A420 browning), and the two inertness
guarantees. What is deliberately NOT pinned: the ~29 % rate figure as a literal — it is
recomputed from the parameters, so a re-sourcing moves the test with the value.
"""

import numpy as np
import pytest

from fermentation.parameters import default_data_dir, load_parameters
from fermentation.runtime import simulate_ensemble
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario

_FERMENT_DAYS = 30.0
_AGING_DAYS = 150.0
_DOSE_MGL = 0.5


def _wine(*, fined: bool, oxygen_mgl: float = 60.0) -> Scenario:
    """A tannic wine aged warm, optionally copper-fined the day before aging begins."""
    interventions: list[Intervention] = []
    if fined:
        interventions.append(
            Intervention(
                day=_FERMENT_DAYS - 1.0,
                action="add_copper",
                params={"copper_mgl": _DOSE_MGL},
            )
        )
    interventions.append(Intervention(day=_FERMENT_DAYS, action="begin_aging"))
    if oxygen_mgl > 0.0:
        interventions.append(
            Intervention(day=_FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": oxygen_mgl})
        )
    return Scenario(
        name="d191-residual-copper",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25, "tannin_gpl": 1.5},
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=_FERMENT_DAYS, celsius=25.0),
        ],
        interventions=interventions,
        duration_days=_FERMENT_DAYS + _AGING_DAYS,
    )


# -- the sourced value and its two band edges ----------------------------------


def test_the_retained_fraction_carries_its_provenance():
    p = load_parameters(default_data_dir() / "additions.yaml")["copper_fining_residual_fraction"]
    assert p.value == 0.95
    assert p.tier.name.lower() == "plausible"
    # The source must be named, and the rendered-not-read flag must survive edits: the >95 %
    # is UWC's rendering of Clark et al. 2015, and that paper was NOT read (the D-190 rule).
    assert "Understanding Wine Chemistry" in p.provenance.source
    assert "NOT read" in p.provenance.source or "not read" in p.provenance.source


def test_the_band_edges_have_the_status_the_source_gives_them():
    """The two edges are NOT symmetric, and a sweep must not treat them as if they were.

    Low 0.95 is PRINTED — the only number the source states. High 1.0 is CONSTRUCTED — the
    physical ceiling, since a dose cannot be more than fully retained. The shipped value sits
    ON the printed edge deliberately, which means a band-edge sweep finds the low arm bitwise
    identical to nominal: that is arithmetic, not inertness (the D-164 trap).
    """
    p = load_parameters(default_data_dir() / "additions.yaml")["copper_fining_residual_fraction"]
    assert p.uncertainty is not None
    assert p.uncertainty.low == 0.95 == p.value  # nominal ON the printed edge
    assert p.uncertainty.high == 1.0  # the physical ceiling
    # A retained fraction above 1 would mean fining CREATES copper.
    assert p.uncertainty.high <= 1.0


# -- the credit itself ---------------------------------------------------------


def test_fining_leaves_the_retained_share_in_the_wine():
    cs = compile_scenario(_wine(fined=True))
    traj = cs.run()
    assert traj.success
    retained = cs.param_values["copper_fining_residual_fraction"]
    background = float(cs.y0[cs.schema.slice("copper")][0])
    expected = background + retained * (_DOSE_MGL / 1000.0)
    assert float(traj.series("copper")[-1]) == pytest.approx(expected, rel=1e-12)


def test_an_unfined_wine_sits_exactly_at_the_multiplier_neutral_point():
    """The control: with no fining, ``f_copper`` is exactly 1 and the D-132/D-133 rate stands."""
    cs = compile_scenario(_wine(fined=False))
    traj = cs.run()
    typical = cs.param_values["copper_typical"]
    assert float(traj.series("copper")[-1]) == pytest.approx(typical, rel=1e-12)
    f_copper = 1.0 + cs.param_values["k_copper_multiplier"] * (
        float(traj.series("copper")[-1]) - typical
    )
    assert f_copper == pytest.approx(1.0, abs=1e-15)


# -- the consequence, at a named site ------------------------------------------


def test_a_fined_wine_browns_faster_than_an_unfined_one():
    """The payoff: fining buys a reduction fix and pays for it in oxidation (D-191).

    Measured EARLY in the aging tail, while the dosed O₂ is still plentiful. The comparison
    interpolates both arms onto a common time in hours — they are separate adaptive solves
    with different step sequences, so a shared array INDEX is a different moment in each (the
    mistake probe 1 made, which reported the fined wine browning *slower*).
    """
    fined = compile_scenario(_wine(fined=True)).run()
    plain = compile_scenario(_wine(fined=False)).run()
    assert fined.success and plain.success

    # POSITIVE CONTROL: the browning term must actually be live, or "no difference" would be
    # a statement about the scenario rather than about copper.
    assert float(plain.series("A420")[-1]) > 0.0

    at_h = (_FERMENT_DAYS + 0.01 * _AGING_DAYS) * 24.0
    a_fined = float(np.interp(at_h, fined.t, np.asarray(fined.series("A420"), dtype=np.float64)))
    a_plain = float(np.interp(at_h, plain.t, np.asarray(plain.series("A420"), dtype=np.float64)))

    # The expected size is RECOMPUTED from the parameters, not pinned as a literal, so a
    # re-sourcing of either the retention or the multiplier moves this test with the value.
    cs = compile_scenario(_wine(fined=True))
    expected_f = 1.0 + cs.param_values["k_copper_multiplier"] * (
        cs.param_values["copper_fining_residual_fraction"] * (_DOSE_MGL / 1000.0)
    )
    assert expected_f > 1.0
    # Browning is faster, and by roughly the rate multiplier. The band is wide on purpose:
    # A420 is an integral of a rate that is already competing for O₂ with the other sinks,
    # so the ratio approaches f_copper from below rather than equalling it.
    assert a_fined > a_plain
    assert a_fined / a_plain == pytest.approx(expected_f, rel=0.10)


def test_the_finite_oxygen_dose_bounds_how_brown_the_wine_can_get():
    """Copper changes the RATE far more than the destination, and that is physics, not a bug.

    With a fixed O₂ dose the total browning is bounded by the oxygen supplied, so the fined
    wine's advantage decays as the stock is spent: it wins a larger SHARE of the same O₂
    rather than browning without limit. Pinning this keeps a future reader from "fixing" the
    small end-of-run difference.
    """
    fined = compile_scenario(_wine(fined=True)).run()
    plain = compile_scenario(_wine(fined=False)).run()
    early_h = (_FERMENT_DAYS + 0.01 * _AGING_DAYS) * 24.0

    early_fined = float(
        np.interp(early_h, fined.t, np.asarray(fined.series("A420"), dtype=np.float64))
    )
    early_plain = float(
        np.interp(early_h, plain.t, np.asarray(plain.series("A420"), dtype=np.float64))
    )
    early_ratio = early_fined / early_plain
    final_ratio = float(fined.series("A420")[-1]) / float(plain.series("A420")[-1])
    assert early_ratio > final_ratio > 1.0


# -- isolability (prime directive #3) ------------------------------------------


def test_a_wine_that_is_never_fined_has_no_event_that_could_write_the_slot():
    """No ``add_copper`` ⇒ no event writes ``copper`` ⇒ the D-132/D-133 rate is untouched.

    Checks the two things that can actually fail: no scheduled event mutates the slot, and the
    slot holds ``copper_typical`` end to end. An earlier draft compiled the same unfined
    scenario TWICE and asserted the trajectories matched — which is a solver-determinism check,
    not an isolability one: both arms carry this change, so it would have stayed green through
    the very regression its name claims to guard (a credit firing on unfined wines).
    """
    cs = compile_scenario(_wine(fined=False))
    traj = cs.run()
    assert not any(e.label.startswith("add_copper") for e in cs.events)

    typical = cs.param_values["copper_typical"]
    copper = np.asarray(traj.series("copper"), dtype=np.float64)
    # Constant at the multiplier's neutral point for the WHOLE run, not merely at the end.
    assert np.all(copper == typical)


def test_the_retention_band_is_not_propagated_by_the_sampler():
    """A scope limit, pinned rather than asserted in prose: the band's width is documentation.

    ``copper_fining_residual_fraction`` is **compile-consumed** — ``_verb_add_copper`` resolves
    ``.value`` when the scenario compiles and the ``mutate`` closure captures a plain float, so
    an ensemble that perturbs the parameter dict cannot move the credited copper. This is
    ``test_drawability_surface``'s class 2 (the ``_closure_otr`` shape), and it means the
    uncertainty band does NOT reach any output through the ensemble machinery.

    Two assertions, because either alone passes for the wrong reason: the draw must genuinely
    happen (distinct values), and a Process-read parameter must MOVE in the same harness, or
    "identical" cannot be told from a harness that perturbed nothing (the D-157 denominator).
    """
    compiled = compile_scenario(_wine(fined=True))
    forced = simulate_ensemble(
        compiled.process_set,
        compiled.parameters,
        compiled.y0,
        compiled.t_span_h,
        n_members=2,
        seed=0,
        only=["copper_fining_residual_fraction"],
        events=compiled.events,
    )
    assert "copper_fining_residual_fraction" in forced.sampled_names
    drawn = [float(m["copper_fining_residual_fraction"]) for m in forced.member_params]
    assert len(set(drawn)) == 2, f"the draw did not vary: {drawn}"
    assert np.array_equal(forced.members[0], forced.members[1])

    # POSITIVE CONTROL, same harness and scenario: the multiplier this parameter feeds DOES move.
    control = simulate_ensemble(
        compiled.process_set,
        compiled.parameters,
        compiled.y0,
        compiled.t_span_h,
        n_members=2,
        seed=0,
        only=["k_copper_multiplier"],
        events=compiled.events,
    )
    assert not np.array_equal(control.members[0], control.members[1])


def test_beer_has_no_copper_slot_so_fining_only_removes_sulfide():
    """The credit is guarded on the slot, exactly as the thiol half is guarded on its pool."""
    beer = Scenario(
        name="d191-beer",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 150.0,
            "pitch_gpl": 0.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[
            Intervention(day=10.0, action="add_copper", params={"copper_mgl": _DOSE_MGL})
        ],
        duration_days=20.0,
    )
    cs = compile_scenario(beer)
    assert "copper" not in cs.schema
    traj = cs.run()
    assert traj.success  # the verb still fines; it simply has no slot to credit
