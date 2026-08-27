"""Wine's nitrogen budget: where the nitrogen enters, and where it ends up (decision D-243).

D-232 §10 entered *"wine's nitrogen budget was not audited the same way; only beer's was"* and
five records carried it forward untouched. Beer's audit (D-230) did two things: it **sourced the
numerator** of ``dX = YAN / f_N`` against a transcribed wort composition, and it **inverted the
identity** to ask what cell nitrogen content would reproduce a measured crop. This module is the
wine half, and its findings are not beer's.

**The numerator is sound, and the source did the hard part itself.** Varela, Pizarro & Agosin
2004 — the project's only independent wine dataset (D-56) — state their own proline subtraction:
total nitrogen was 380 / 65 mg/L, *"since proline is not assimilable by yeast under anaerobic
growth conditions, the ammonium salt and α-amino acids … were considered assimilable sources of
nitrogen. Therefore … the concentrations of assimilable nitrogen were 300 and 50 mg liter−1."*
That is the exact check that closed beer's first candidate on Peyer Table 16, and it lands the
same safe way: no subtraction is owed, the engine's ``yan_mgl`` is assimilable nitrogen.

**The inversion gives the OPPOSITE verdict to beer's, and the symmetry must not be written.**
D-230 refused beer's partition candidate because reproducing the counted crop demanded cells at
20.2-26.2 % nitrogen — 1.77x outside the band's high edge and roughly double anything a cell can
be. Wine's demand points the other way, and there is **no comparable floor**: ``wine_generic.yaml``
says 0.114 is the *N-replete* reference which "drops under nitrogen limitation", and the shipped
Coleman regression already puts f_N at 0.0362 at YAN=50. Varela's measured crop demands 0.0521 at
YAN=300 — *between* the engine's own two values. **A nitrogen-budget explanation for wine's gap is
ADMISSIBLE.** Do not cite D-230's inadmissibility argument on the wine side.

**The finding: assimilable nitrogen has TWO entry channels and the yield fit reads ONE.**
``yan_mgl`` seeds the ``N`` slot *and* is the point at which Coleman's ``Y_X/N(N_init)``
regression is evaluated (D-14). ``amino_acids_gpl`` seeds eight speciated pools that are **also**
on the nitrogen ledger (D-100, weighted in ``total_nitrogen``). They **add**; they do not
partition — and D-32's own text says *"amino acids **are** part of YAN"*, which is the premise the
seam does not implement. D-14 predates D-32 and its "in our model all assimilated nitrogen enters
biomass" was true when written; D-32 added the second store underneath it and neither record
states the consequence.

Nothing here is a conservation defect — the ledger closes to ~1e-14 relative on every arm
(:func:`test_the_nitrogen_ledger_closes_through_every_wine_sink`). It is a **declaration**
defect: a scenario says 250 and the run carries 362.7.

**These are characterization pins, in the Varela-benchmark style: they assert what the engine
CURRENTLY does, so a silent change is caught in either direction.** A repair that makes the two
channels partition SHOULD turn them red — update them and the decision together, do not widen
them. And note the naive repair does not work: summing the channels into the fit leaves Coleman's
own 70-350 mg N/L range at the suite's commonest dose and reaches f_N = 0.379 at 2 g/L, which is
three times anything physiological. The evaluation point is a design decision, not a beat's call.

Measurements: ``M:\\claud_projects\\temp\\ferment\\d243-wine-nitrogen-audit\\``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from fermentation.core.chemistry import nitrogen_mass_fraction
from fermentation.core.tiers import Tier
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario.compile import CompiledScenario
from fermentation.validation import assert_conserved, total_nitrogen

#: Every slot that carries nitrogen into a wine run through the amino-acid channel, with the
#: species ``total_nitrogen`` weights it at. Mirrors :func:`fermentation.validation.conservation
#: .total_nitrogen` — kept as a mapping here so a ninth pool shows up as a KeyError-free but
#: UNDERCOUNTED total rather than silently, which
#: :func:`test_the_amino_acid_channel_census_names_every_nitrogen_bearing_pool` forbids.
_AA_CHANNEL = {
    "amino_acids": "arginine",
    "amino_acids_generic": "glutamine",
    "leucine": "leucine",
    "isoleucine": "isoleucine",
    "valine": "valine",
    "threonine": "threonine",
    "phenylalanine": "phenylalanine",
    "methionine": "methionine",
}

_FERMENT_DAYS = 30.0
_AGING_DAYS = 150.0


def _wine(
    *,
    yan_mgl: float = 250.0,
    amino_acids_gpl: float = 0.0,
    autolysis_rate_per_h: float = 0.0,
    oxygen: bool = False,
) -> Scenario:
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if amino_acids_gpl:
        initial["amino_acids_gpl"] = amino_acids_gpl
    if autolysis_rate_per_h:
        initial["autolysis_rate_per_h"] = autolysis_rate_per_h
    interventions = [Intervention(day=_FERMENT_DAYS, action="begin_aging")]
    if oxygen:
        interventions.append(
            Intervention(day=_FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": 60.0})
        )
    return Scenario(
        name=f"n-budget-{yan_mgl:g}-{amino_acids_gpl:g}",
        medium="wine",
        initial=initial,
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=_FERMENT_DAYS, celsius=25.0),
        ],
        interventions=interventions,
        duration_days=_FERMENT_DAYS + _AGING_DAYS,
    )


def _channel_nitrogen(compiled: CompiledScenario) -> tuple[float, float]:
    """(ammonium channel, amino-acid channel) nitrogen at t=0, both in g N/L."""
    y0 = np.asarray(compiled.y0)
    schema = compiled.schema
    ammonium = float(y0[schema.slice("N")][0])
    amino = sum(
        float((y0[schema.slice(slot)] * nitrogen_mass_fraction(species)).sum())
        for slot, species in _AA_CHANNEL.items()
        if slot in schema
    )
    return ammonium, amino


def test_the_amino_acid_channel_census_names_every_nitrogen_bearing_pool():
    """:data:`_AA_CHANNEL` must stay the WHOLE amino-acid channel, or the tests below undercount.

    The failure this forbids is silent: a ninth speciated pool would be weighted in
    ``total_nitrogen`` and simply missing from :func:`_channel_nitrogen`, so the "declared vs
    actual" gap below would read SMALLER than it is and the pin would still pass. Scored against
    ``total_nitrogen``'s own weight vector rather than against a second hand-written list.
    """
    compiled = compile_scenario(_wine(amino_acids_gpl=0.5), strict=True)
    schema = compiled.schema
    f_n = compiled.param_values["biomass_N_fraction"]
    y0 = np.asarray(compiled.y0)

    biomass_slots = {"X", "X_dead", "X_mlf", "X_mlf_dead", "X_brett", "X_brett_dead"}
    total = float(np.asarray(total_nitrogen(schema, biomass_nitrogen_fraction=f_n)(y0)).sum())
    ammonium, amino = _channel_nitrogen(compiled)
    biomass = sum(float(y0[schema.slice(s)].sum()) * f_n for s in biomass_slots if s in schema)
    assert total == pytest.approx(ammonium + amino + biomass, rel=1e-12), (
        "t=0 nitrogen does not decompose into ammonium + amino-acid channel + biomass — a "
        f"nitrogen-bearing slot exists that `_AA_CHANNEL` does not name (total {total:.9f}, "
        f"parts {ammonium:.9f} + {amino:.9f} + {biomass:.9f})"
    )


def test_the_declared_yan_is_the_run_s_actual_assimilable_nitrogen_and_the_fit_point():
    """The D-243 gap, CLOSED at D-244 — and pinned so it cannot silently re-open.

    **This test is the inverse of the one it replaces.** D-243 found two entry channels and one
    evaluation point: ``yan_mgl`` seeded the ammonium slot *and* was where Coleman's regression
    was evaluated, while ``amino_acids_gpl`` seeded eight more nitrogen-bearing pools on top, so
    a wine declaring 250 mg N/L carried 362.7 and had its yield fitted for the poorer must. The
    predecessor pinned that gap at 1.45x. D-244 makes ``yan_mgl`` the must's TOTAL assimilable
    nitrogen and carves the amino-acid pools out of it, so all three numbers collapse onto each
    other: declared == actual == the point the fit is evaluated at.

    A RED here means the carve-out broke. The assertions fail for different causes and say so: the
    ammonium remainder is wrong (``_wine_ammonium_gpl``), the channel census misses a pool (see the
    census test), or the fit stopped reading the declared total.

    **The declared == actual assertion is an ARITHMETIC IDENTITY and is not the load-bearing one.**
    ``N`` is seeded as ``declared - amino`` and :func:`_channel_nitrogen` adds ``N`` back to the
    same pools at the same weights, so the sum closes whatever the must spectrum does. Measured:
    mutating ``must_aa_fraction_arginine`` 0.380 -> 0.470 (inside its own [0.25, 0.48] band, so the
    mutation RUNS rather than breaking ``load_parameters``) leaves it green and fires the two
    HALF-pins below. Those halves are the content; the identity is here because it is the sentence
    the repair claims, and a reader must not mistake it for the guard.
    """
    compiled = compile_scenario(_wine(amino_acids_gpl=0.5), strict=True)
    ammonium, amino = _channel_nitrogen(compiled)
    declared_mgl = 250.0
    actual_mgl = (ammonium + amino) * 1000.0

    assert actual_mgl == pytest.approx(declared_mgl, abs=1e-9), (
        f"a wine declaring yan_mgl=250 with amino_acids_gpl=0.5 carries {actual_mgl:.4f} mg N/L "
        "of assimilable nitrogen. Since D-244 the declaration IS the total, so any difference "
        "is the carve-out failing — not a channel to re-characterize"
    )
    # The carve-out is a SPLIT of one number, so pin both halves, not just the sum: a sum that
    # closes with both halves wrong is exactly what a one-sided read would miss.
    assert amino * 1000.0 == pytest.approx(112.7, abs=0.5), (
        f"the amino-acid channel carries {amino * 1000.0:.1f} mg N/L, not the characterized "
        "112.7 — the must amino-acid spectrum moved; re-measure, do not re-band"
    )
    assert ammonium * 1000.0 == pytest.approx(137.3, abs=0.5), (
        f"the ammonium remainder is {ammonium * 1000.0:.1f} mg N/L, not the characterized 137.3"
    )

    # The fit is evaluated at that same total. Asserted as an equality against a recompile of an
    # UNDOSED wine at the same declared YAN — not against a literal — so Coleman's coefficients
    # may move without touching the claim. Note this equality held before D-244 too and meant the
    # OPPOSITE thing: there it said the fit ignored the dose, here it says the dose is inside the
    # number the fit reads. The two halves above are what tell them apart.
    undosed = compile_scenario(_wine(amino_acids_gpl=0.0), strict=True)
    assert (
        compiled.param_values["biomass_N_fraction"] == undosed.param_values["biomass_N_fraction"]
    ), "biomass_N_fraction moved off the declared YAN — the D-14 evaluation point has shifted"

    # And the gap D-243 measured is GONE, stated as its own quantity so the repair cannot be
    # quietly reverted: the fit's value and the value at the run's own nitrogen are now one.
    a0 = compiled.parameters["biomass_N_yield_log_intercept"].value
    a1 = compiled.parameters["biomass_N_yield_log_slope"].value
    f_n_at_actual = 1.0 / math.exp(a0 + a1 * actual_mgl)
    assert f_n_at_actual / compiled.param_values["biomass_N_fraction"] == pytest.approx(1.0), (
        "the evaluation-point error is back — D-243 measured it at 1.502x and D-244 closed it"
    )


def test_the_dose_cannot_out_run_the_declaration_and_the_fit_is_held_at_colemans_edge():
    """The two refusals-turned-behaviours D-244 ships, pinned together because they are a pair.

    A must cannot hold more amino-acid nitrogen than it holds assimilable nitrogen, so the seam
    REFUSES rather than flooring the ammonium at zero (which would put the run straight back into
    the state D-243 named: carrying nitrogen it never declared). The refusal is what forces a
    pre-D-244 scenario to declare its real total — and that total is often above the 70-350 mg
    N/L span Coleman fitted, which before D-244 did not compile AT ALL: ``f_N`` left
    ``biomass_N_fraction``'s [0.03, 0.15] bracket at 444.0 mg N/L and raised an opaque pydantic
    band error naming neither nitrogen nor Coleman. Holding the fit at the fitted edge is what
    makes the migration possible, so the two ship together or neither does.

    **The hold is epistemic, and the tier is the claim.** Nothing here says Y_X/N saturates.
    """
    with pytest.raises(ValueError, match="more than the declared yan_mgl"):
        compile_scenario(_wine(yan_mgl=250.0, amino_acids_gpl=2.0), strict=True)

    # A must ABOVE Coleman's span compiles, holds, and says so in its tier.
    high = compile_scenario(_wine(yan_mgl=500.0, amino_acids_gpl=2.0), strict=True)
    edge = compile_scenario(_wine(yan_mgl=350.0), strict=True)
    assert high.param_values["biomass_N_fraction"] == pytest.approx(
        edge.param_values["biomass_N_fraction"]
    ), "the fit is not held at biomass_N_yield_fit_yan_max — extrapolation is back"
    assert high.parameters["biomass_N_fraction"].tier is Tier.SPECULATIVE, (
        "a held fit must carry the admission in its tier, or the hold is invisible downstream"
    )

    # Anti-vacuity: the hold must not fire INSIDE the span, and must not drag the tier down there.
    assert edge.parameters["biomass_N_fraction"].tier is not Tier.SPECULATIVE, (
        "the hold is firing at the edge itself — every wine in the suite would go speculative"
    )
    assert edge.param_values["biomass_N_fraction"] < 0.15, (
        "the held value left the bracket it exists to keep the fit inside"
    )

    # The LOW side is deliberately NOT held, and this is the pin that forbids adding one: it
    # would move Varela's 50 mg N/L arm, the project's only independent wine dataset (D-56).
    low = compile_scenario(_wine(yan_mgl=50.0), strict=True)
    assert low.param_values["biomass_N_fraction"] == pytest.approx(0.036171, abs=1e-6), (
        "the 50 mg N/L arm moved — a low-side hold was added, which re-anchors the Varela "
        "comparison on a number Coleman never fitted there"
    )
    assert low.parameters["biomass_N_fraction"].tier is not Tier.SPECULATIVE


def test_all_assimilable_nitrogen_reaches_biomass_whatever_channel_it_entered_by():
    """The D-14 identity, generalised — and the reason the amino-acid dose IS extra biomass.

    D-32 designed the swap to be nitrogen-neutral *through the transfer*, and it is. But the pool
    it debits drains to ``N``, and ``N`` drains to biomass, so cumulatively the amino-acid pool is
    a second nitrogen reservoir rather than a substitute for the first. On a bare run the identity
    is exact; with a dose it holds to a residue that is pool nitrogen the aging tail never
    re-assimilates.

    **THE RESIDUE GREW TWELVE-FOLD AT D-244 -- 0.6 % to 7.8 % -- and that was the carve-out
    working, not a defect.** Before D-244 a dosed wine's nitrogen was mostly ammonium (250 mg
    N/L of it, with the dose's 112.7 on top), and ammonium is what the growth Process reads
    directly, so almost all of it arrived. The same wine then held 137.3 mg N/L of ammonium and
    112.7 in the amino-acid pools, and pool nitrogen only reached biomass through the D-32 swap,
    which did not run to completion.

    **AND D-248 CLOSES IT: 7.75 % -> 0.09 %.** The swap's incompleteness was never a property of
    the swap's fraction psi -- it was that psi*gate*f_N*base_dx is strictly below growth's own
    f_N*base_dx draw, so ammonium could only fall and the pools froze once it hit zero.
    :class:`~fermentation.core.kinetics.amino_acids.AssimilableNitrogenUptake` un-couples uptake
    from the growth rate, and the amino-acid channel now delivers essentially all of its nitrogen
    to biomass. The identity D-14 states is a statement about BOTH channels again, and the
    "45 % must be assimilated rather than taken up" reading is retired: it is still assimilated
    through a second step, but the step now completes. Pinned two-sided as before -- a residue
    with only an upper bound cannot catch the channel slowing back down.

    **The bare arm's tolerance is 5e-9, not 1e-9, and the reason is the SOLVER rather than the
    model.** At ``aa = 0`` this Process is disabled at the compile seam and contributes exactly
    nothing -- measured, the arm is bit-identical at capacity 0 and capacity 1. What moved it by
    2.3e-9 is that D-248 adds a state slot: BDF's error norm is RMS-weighted over the state
    vector, so 97 -> 98 slots shifts step selection with no model change at all. The same effect
    is recorded at ``tests/test_oxidative_cascade_guards._PIN_RTOL`` for the ``quinone`` slot.
    """
    for aa, tol in ((0.0, 5e-9), (0.5, 8.0e-2)):
        compiled = compile_scenario(_wine(amino_acids_gpl=aa), strict=True)
        f_n = compiled.param_values["biomass_N_fraction"]
        ammonium, amino = _channel_nitrogen(compiled)
        x0 = float(np.asarray(compiled.y0)[compiled.schema.slice("X")][0])
        traj = compiled.run()
        assert traj.success, traj.message
        biomass = float(np.asarray(traj.series("X"))[-1] + np.asarray(traj.series("X_dead"))[-1])
        predicted = x0 + (ammonium + amino) / f_n
        assert biomass == pytest.approx(predicted, rel=tol), (
            f"aa={aa}: biomass {biomass:.6f} g/L is not X0 + (both channels)/f_N = "
            f"{predicted:.6f} — a nitrogen sink is holding mass the identity says reaches X"
        )
        if aa:
            # The residue as its own quantity, two-sided. A one-sided tolerance would stay
            # green if the swap got FASTER and the amino-acid channel stopped being
            # distinguishable from the ammonium one, which is what this pin exists to catch.
            assert 1.0 - biomass / predicted == pytest.approx(0.00092, abs=5e-4), (
                f"aa={aa}: the unassimilated residue is {1.0 - biomass / predicted:.5f}, not "
                "the D-248 0.00092 (D-244 characterized 0.0775, before assimilable-nitrogen "
                "uptake was un-coupled from growth demand) — the amino-acid channel's "
                "completeness moved, in one direction or the other; measure which first"
            )


def test_the_nitrogen_ledger_closes_through_every_wine_sink():
    """Closure, on the four arms that turn the optional sinks on. This is NOT the finding.

    Wine's code has twelve sites writing the ``N`` slot against beer's one, and the audit's first
    question is whether they balance. They do, to ~1e-14 relative. Recorded so the *declaration*
    finding above is not misread as a conservation break.
    """
    for label, scenario in (
        ("bare", _wine()),
        ("aa", _wine(amino_acids_gpl=0.5)),
        ("aa+autolysis", _wine(amino_acids_gpl=0.5, autolysis_rate_per_h=1e-4)),
        (
            "aa+autolysis+o2",
            _wine(amino_acids_gpl=0.5, autolysis_rate_per_h=1e-4, oxygen=True),
        ),
    ):
        compiled = compile_scenario(scenario, strict=True)
        traj = compiled.run()
        assert traj.success, traj.message
        assert_conserved(
            traj,
            total_nitrogen(
                compiled.schema,
                biomass_nitrogen_fraction=compiled.param_values["biomass_N_fraction"],
            ),
            label=f"nitrogen {label}",
        )


def test_autolysis_leaves_nearly_half_the_nitrogen_outside_biomass():
    """The sink census's one surprise, pinned — a pre-registered "no sink over 5 %" was wrong 5x.

    On a bare or merely aa-dosed run essentially all nitrogen ends in biomass (99.4 %+). Turn
    autolysis on and 49.7 % of it ends in the amino-acid pools instead: lees self-digestion
    returns biomass nitrogen to the pool faster than anything re-assimilates it once the ferment
    is dry. Pinned because it is the only configuration in which wine's nitrogen budget is *not*
    "biomass and rounding", and because a change to the autolysis/re-assimilation balance would
    otherwise move it silently.
    """
    compiled = compile_scenario(_wine(amino_acids_gpl=0.5, autolysis_rate_per_h=1e-4), strict=True)
    schema = compiled.schema
    f_n = compiled.param_values["biomass_N_fraction"]
    traj = compiled.run()
    assert traj.success, traj.message
    y0, y_end = np.asarray(compiled.y0), np.asarray(traj.y)[:, -1]
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)
    total0 = float(np.asarray(n_of(y0)).sum())

    in_pools = sum(
        float((y_end[schema.slice(slot)] * nitrogen_mass_fraction(species)).sum())
        for slot, species in _AA_CHANNEL.items()
        if slot in schema
    )
    # D-243 characterized this at 45.9 %. D-244's carve-out moved it to 49.7 %: the same must
    # started with 112.7 of its 250 mg N/L already IN the pools, and the swap did not clear all
    # of it before autolysis began refilling them. D-248 moves it to 38.6 %, and the direction is
    # the informative part: un-coupling assimilable-nitrogen uptake from growth demand clears the
    # MUST's pool nitrogen almost completely during fermentation, so what stands here afterwards
    # is much more nearly autolysis's own contribution rather than a dose the yeast never
    # finished. This is now the closest the number has been to measuring the thing its name says.
    assert in_pools / total0 == pytest.approx(0.386, abs=0.02), (
        f"autolysis now leaves {100 * in_pools / total0:.1f} % of the nitrogen in the amino-acid "
        "pools, not the D-248 ~38.6 % (D-244 characterized 49.7 %) — the autolysis/"
        "re-assimilation balance moved"
    )

    # The same run WITHOUT autolysis is the control: without it the 49.7 % is just a number, and
    # nothing shows it is autolysis that puts the nitrogen there.
    control = compile_scenario(_wine(amino_acids_gpl=0.5), strict=True)
    ctl = control.run()
    assert ctl.success, ctl.message
    ctl_end = np.asarray(ctl.y)[:, -1]
    ctl_pools = sum(
        float((ctl_end[control.schema.slice(slot)] * nitrogen_mass_fraction(species)).sum())
        for slot, species in _AA_CHANNEL.items()
        if slot in control.schema
    )
    # The control's floor moved with the carve-out at D-244 (an aa-dosed wine ended with the
    # 7.8 % swap residue in its pools even with autolysis OFF, where it had been under 1 % when
    # the dose was extra ammonium's worth of nitrogen) and moved BACK at D-248: measured
    # **0.0915 %**, because un-coupled uptake clears the must's own pool nitrogen during
    # fermentation. So the control is a near-zero floor again and the attribution is at its
    # strongest -- ~420x smaller than the arm above, against 6x at D-244. The bound is tightened
    # to match, because a 0.10 ceiling would now pass on a control that had regressed the whole
    # way back to D-244's residue and would attribute that to autolysis.
    assert ctl_pools / total0 < 0.005, (
        "the autolysis-off control now leaves nitrogen in the pools too, so the arm above no "
        f"longer attributes anything to autolysis ({100 * ctl_pools / total0:.3f} %, against "
        "D-248's measured 0.0915 %)"
    )
