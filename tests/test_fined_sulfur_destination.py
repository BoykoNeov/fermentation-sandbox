"""Copper fining complexes the sulfur, it does not destroy it (D-193).

D-191 credited the *metal* back to the wine and left the *sulfur* half measured but unbuilt, with
D-45's carbon booking flagged as resting on the same retracted precipitation mechanism. This beat
spends both. ``add_copper`` now transfers what it binds into the D-135 reservoirs —
``h2s -> bound_h2s`` and ``methanethiol -> bound_methanethiol`` — instead of annihilating it, so a
fined wine can turn reductive again in the bottle (UWC §26.2.4.1: the complexes regenerate H₂S in
storage, the "reappearance of off aromas" in that chapter's own consequences table) and the fining
stops booking thiol carbon out of the wine.

**What is pinned here.** The transfer is whole-mass and exact; the free pools are untouched at the
event so the odour fix is unchanged; the carbon ledger is neutral to machine precision; the aged
consequence at a named site; and the isolability guarantees. What is deliberately NOT pinned as a
literal: the 3-year free-H₂S figures, which are a property of this scenario rather than of the
mechanism — the ORDERING against the unfined control is pinned instead.
"""

import numpy as np
import pytest

from fermentation.parameters import default_data_dir, load_parameters
from fermentation.runtime.schedule import ExternalFlow, ScheduledTrajectory
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.validation.conservation import total_carbon

_FERMENT_DAYS = 30.0
_FINE_DAY = 25.0
_DOSE_MGL = 0.5
_AGING_DAYS = 365.0 * 3.0


def _reductive_wine(*, fined: bool, aged: bool = True, days: float = _AGING_DAYS) -> Scenario:
    """A wine with autolysis opted in, so BOTH reductive pools build before the fining.

    ``autolysis_rate_per_h`` is what makes this scenario able to answer the question at all: with
    it off, ``methanethiol`` is identically 0 and the thiol half of every assertion below would
    pass vacuously.
    """
    interventions: list[Intervention] = []
    if fined:
        interventions.append(
            Intervention(day=_FINE_DAY, action="add_copper", params={"copper_mgl": _DOSE_MGL})
        )
    if aged:
        interventions.append(Intervention(day=_FERMENT_DAYS, action="begin_aging"))
    return Scenario(
        name="d193-fined-sulfur",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 80.0,
            "pitch_gpl": 0.25,
            "autolysis_rate_per_h": 0.002,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=interventions,
        duration_days=_FERMENT_DAYS + days,
    )


def _fining_flow(traj: ScheduledTrajectory) -> ExternalFlow:
    return next(f for f in traj.external_flows if f.label.startswith("add_copper"))


# -- the transfer itself -------------------------------------------------------


def test_the_bound_sulfur_lands_in_the_reservoirs_mass_for_mass():
    """Whole mass, both species, no retention fraction applied — the D-193 shape.

    The transfer is 1:1 by construction (one molecule changing binding state), which is why no
    parameter mediates it. Both halves are asserted with a positive control that the removal was
    non-zero, since ``0 == -0`` would satisfy the equality on a wine with nothing to bind.
    """
    cs = compile_scenario(_reductive_wine(fined=True, aged=False, days=15.0))
    traj = cs.run()
    schema = cs.schema
    flow = _fining_flow(traj)

    for free, bound in (("h2s", "bound_h2s"), ("methanethiol", "bound_methanethiol")):
        removed = -float(flow.delta[schema.slice(free)][0])
        gained = float(flow.delta[schema.slice(bound)][0])
        assert removed > 0.0, f"{free} removal was vacuous — nothing to transfer"
        assert gained == pytest.approx(removed, rel=1e-15), free


def test_no_retention_fraction_is_applied_to_the_sulfur():
    """The copper keeps 0.95 of its dose; the sulfur keeps ALL of what was bound. Deliberate.

    0.95 is a PRINTED LOWER BOUND on copper retention measured after *filtering* or *racking* —
    operations ``add_copper`` does not model (the repo has a separate ``rack`` verb, and the
    parameter's own band note calls those cellar sinks "separate operations on longer timescales,
    not part of this event"). Scaling the sulfur by it was built as a probe and rejected: it turns
    a bound into a loss fraction and leaves an invented 5 % as the only carbon outflow. This test
    fails if a future edit "harmonises" the two shares.
    """
    cs = compile_scenario(_reductive_wine(fined=True, aged=False, days=15.0))
    traj = cs.run()
    schema = cs.schema
    flow = _fining_flow(traj)
    retained = cs.param_values["copper_fining_residual_fraction"]
    assert retained < 1.0  # or the distinction below is untestable

    # copper: the SCALED share of the dose
    assert float(flow.delta[schema.slice("copper")][0]) == pytest.approx(
        retained * (_DOSE_MGL / 1000.0), rel=1e-12
    )
    # sulfur: the WHOLE bound mass. Asserted as a RATIO against 1, and separately as NOT the
    # retained fraction. An earlier draft wrote `gained > retained * removed`, which is satisfied
    # by any fraction above 0.95 — including 0.96 — so it would have passed on the very thing this
    # test forbids while reading like coverage. That is this beat's own lesson turned on itself
    # (D-193 §8): an assert that passes on the forbidden case is not a guard.
    for free, bound in (("h2s", "bound_h2s"), ("methanethiol", "bound_methanethiol")):
        removed = -float(flow.delta[schema.slice(free)][0])
        gained = float(flow.delta[schema.slice(bound)][0])
        assert removed > 0.0, free  # positive control: 0/0 would make the ratio meaningless
        assert gained / removed == pytest.approx(1.0, rel=1e-15), free
        assert gained / removed != pytest.approx(retained, rel=1e-3), free


def test_the_free_pools_are_untouched_at_the_event():
    """The odour fix is unchanged: only the destination differs, never the removal.

    Compares the removal arithmetic against the same verb applied to a hand-built state, so a
    regression that "helpfully" left some sulfide behind to feed the reservoir would fail here.
    """
    cs = compile_scenario(_reductive_wine(fined=True, aged=False, days=15.0))
    ev = next(e for e in cs.events if e.label.startswith("add_copper"))
    assert ev.mutate is not None
    schema = cs.schema
    y = cs.y0.copy()
    y[schema.slice("h2s")] = 1.0e-3  # 1 mg/L, far above the dose's capacity
    y[schema.slice("methanethiol")] = 1.0e-3
    out = ev.mutate(schema, y, cs.param_values)

    capacity_h2s = (_DOSE_MGL / 1000.0) * cs.param_values["copper_h2s_binding"]
    assert float(out[schema.slice("h2s")][0]) == pytest.approx(1.0e-3 - capacity_h2s, rel=1e-12)
    # H₂S soaked up the whole dose, so no copper is left for the thiol — untouched, and its
    # reservoir therefore gains exactly nothing (the transfer cannot invent a bind).
    assert float(out[schema.slice("methanethiol")][0]) == pytest.approx(1.0e-3, rel=1e-15)
    assert float(out[schema.slice("bound_methanethiol")][0]) == pytest.approx(
        float(y[schema.slice("bound_methanethiol")][0]), rel=1e-15
    )


def test_a_negative_undershoot_transfers_nothing():
    """The ≥0 clamp must not let a solver undershoot CREATE reservoir out of a negative pool."""
    cs = compile_scenario(_reductive_wine(fined=True, aged=False, days=15.0))
    ev = next(e for e in cs.events if e.label.startswith("add_copper"))
    assert ev.mutate is not None
    schema = cs.schema
    y = cs.y0.copy()
    y[schema.slice("h2s")] = -1.0e-9
    before = float(y[schema.slice("bound_h2s")][0])
    out = ev.mutate(schema, y, cs.param_values)
    assert float(out[schema.slice("h2s")][0]) == pytest.approx(-1.0e-9, rel=1e-12)
    assert float(out[schema.slice("bound_h2s")][0]) == before


# -- the ledger ----------------------------------------------------------------


def test_the_fining_moves_no_carbon_now_that_the_mercaptide_stays():
    """D-45's carbon booking, corrected — the half D-191's ``Flags:`` marker identified.

    The two slots carry the identical carbon weight precisely so a bound↔free transfer closes
    exactly, so the fining's carbon flow is now zero to machine precision rather than the
    ~1e-06 g/L outflow it used to book.
    """
    cs = compile_scenario(_reductive_wine(fined=True, aged=False, days=15.0))
    traj = cs.run()
    c_of = total_carbon(cs.schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    flow = _fining_flow(traj)

    # positive control first: a thiol transfer really happened, so the zero is not vacuous
    assert float(flow.delta[cs.schema.slice("bound_methanethiol")][0]) > 0.0
    assert c_of(flow.delta) == pytest.approx(0.0, abs=1e-18)
    # and the state's own carbon no longer drops across the jump
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0), abs=1e-6)


# -- the consequence -----------------------------------------------------------


def test_a_fined_wine_gives_its_sulfide_back_in_the_bottle():
    """The behaviour that was previously impossible: fining is no longer permanent.

    Against an ARM-LEVEL control (the same wine, un-fined) rather than a literal, because the
    3-year figures belong to this scenario. Two orderings are pinned, and they point opposite
    ways on purpose: the fined wine ends up **less** reductive than the unfined one (the fining
    worked), but **more** reductive than the pre-D-193 model said (the sulfur came back).
    """
    fined = compile_scenario(_reductive_wine(fined=True)).run()
    unfined = compile_scenario(_reductive_wine(fined=False)).run()
    assert fined.success and unfined.success

    free_fined = float(fined.series("h2s")[-1])
    free_unfined = float(unfined.series("h2s")[-1])
    bound_fined = float(fined.series("bound_h2s")[-1])
    bound_unfined = float(unfined.series("bound_h2s")[-1])

    # the reservoir carries the fined sulfide, and it is a large addition, not a rounding one
    assert bound_fined > 1.5 * bound_unfined
    # release really ran: the end of the run sits BELOW the post-fining peak. (Compared against
    # the peak, not against `series[0]` — index 0 is the seeded level at pitch, which the fining
    # then jumps well above.)
    assert bound_fined < np.asarray(fined.series("bound_h2s"), dtype=np.float64).max()
    # fining still helps: the fined wine is the less reductive of the two at 3 years
    assert free_fined < free_unfined


def test_the_reservoir_is_inert_until_aging_begins():
    """Fining a wine that is never aged seeds the reservoir but changes no free pool.

    The release Process is compile-disabled without ``begin_aging``, so the routed sulfur sits
    there — which is what makes this change isolable rather than a new fermentation-phase term.
    """
    cs = compile_scenario(_reductive_wine(fined=True, aged=False, days=15.0))
    traj = cs.run()
    bound = np.asarray(traj.series("bound_h2s"), dtype=np.float64)
    seeded = float(cs.y0[cs.schema.slice("bound_h2s")][0])

    jumped = bound.max()
    assert jumped > seeded  # the fining did add to it
    # ...and once added, nothing draws it down again: no release outside the aging phase.
    after = bound[bound > seeded + 1e-15]
    assert np.allclose(after, jumped, rtol=0.0, atol=1e-18)


# -- isolability (prime directive #3) ------------------------------------------


def test_an_unfined_wine_is_untouched_by_this_change():
    """No ``add_copper`` ⇒ no event writes either reservoir ⇒ D-135's model stands unchanged."""
    cs = compile_scenario(_reductive_wine(fined=False))
    traj = cs.run()
    assert not any(e.label.startswith("add_copper") for e in cs.events)
    for bound, initial in (
        ("bound_h2s", "bound_h2s_initial"),
        ("bound_methanethiol", "bound_methanethiol_initial"),
    ):
        seeded = float(cs.y0[cs.schema.slice(bound)][0])
        assert seeded == pytest.approx(cs.param_values[initial] / 1.0e6, rel=1e-12)  # µg/L → g/L
        # a first-order release can only shrink it — nothing added to the pool
        assert float(traj.series(bound)[-1]) < seeded


def test_the_transfer_reads_no_parameter_so_it_cannot_be_tuned():
    """A structural guarantee worth pinning: no constant mediates the destination.

    The routing is 1:1 by construction. If a future edit introduced a "fined sulfide retention"
    parameter, this test would not catch it directly — but the two release rate constants are
    what D-135 refused to couple to copper, and those must stay copper-free. Checked at the
    Process's own declared ``reads``, which is the contract the ProcessSet enforces.
    """
    from fermentation.core.kinetics.aging import (
        BoundHydrogenSulfideRelease,
        BoundMethanethiolRelease,
    )

    for proc in (BoundHydrogenSulfideRelease(), BoundMethanethiolRelease()):
        assert not any("copper" in name for name in proc.reads), proc.name

    # And the reservoirs' seed levels are unchanged by this beat — the fining ADDS to them.
    params = load_parameters(default_data_dir() / "bound_sulfides.yaml")
    assert params["bound_h2s_initial"].value == 19.7
    assert params["bound_methanethiol_initial"].value == 1.4


def test_beer_keeps_the_removal_only_behaviour():
    """Guarded on the slots: beer has neither reservoir, so the verb still simply removes."""
    beer = Scenario(
        name="d193-beer",
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
    assert "bound_h2s" not in cs.schema
    assert "bound_methanethiol" not in cs.schema
    traj = cs.run()
    assert traj.success
