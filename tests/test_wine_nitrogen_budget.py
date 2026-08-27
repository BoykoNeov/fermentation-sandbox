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


def test_the_yield_fit_is_evaluated_at_the_declared_yan_not_the_run_s_actual_nitrogen():
    """Forbid the two-channel gap being closed, widened, or restated without a decision.

    ``biomass_N_fraction`` is overridden at compile from Coleman's regression evaluated at
    ``yan_mgl`` alone (D-14). Dosing ``amino_acids_gpl`` adds assimilable nitrogen that the
    growth Process eventually reaches through the D-32 swap, so the run's nitrogen leaves the
    point its own yield was fitted at. The two errors compound in the same direction: more
    nitrogen, and a yield appropriate to a poorer must.

    Pinned at the suite's commonest dose. A RED here means one of three things, and the message
    should say which: the fit moved to total nitrogen (the repair — update D-243), the channels
    started partitioning (also the repair), or the must amino-acid spectrum moved and changed how
    much nitrogen 0.5 g/L carries (re-measure, do not re-band).
    """
    compiled = compile_scenario(_wine(amino_acids_gpl=0.5), strict=True)
    ammonium, amino = _channel_nitrogen(compiled)
    declared_mgl = ammonium * 1000.0
    actual_mgl = (ammonium + amino) * 1000.0

    assert declared_mgl == pytest.approx(250.0, abs=1e-9), (
        "the ammonium channel is no longer the declared yan_mgl"
    )
    assert actual_mgl == pytest.approx(362.7, abs=0.5), (
        f"a wine declaring yan_mgl=250 with amino_acids_gpl=0.5 now carries {actual_mgl:.1f} "
        "mg N/L of assimilable nitrogen, not the characterized 362.7 — the amino-acid channel's "
        "nitrogen content moved"
    )
    assert actual_mgl / declared_mgl == pytest.approx(1.451, abs=5e-3), (
        "the declared-vs-actual nitrogen ratio moved off its characterized 1.45x"
    )

    # The fit is evaluated at the DECLARED number. This is the defect, pinned as an equality
    # against a recompile at the SAME declared YAN with no dose — not against a literal, so the
    # regression's own coefficients may move without touching this claim.
    undosed = compile_scenario(_wine(amino_acids_gpl=0.0), strict=True)
    assert (
        compiled.param_values["biomass_N_fraction"] == undosed.param_values["biomass_N_fraction"]
    ), (
        "biomass_N_fraction now responds to the amino-acid dose — the yield fit has moved off "
        "the declared YAN. That is the repair D-243 declined to pick; update the decision."
    )

    # And what it would have been at the run's own nitrogen, so the size of the gap is on record
    # rather than inferred. Coleman's regression, evaluated where the run actually sits.
    a0 = compiled.parameters["biomass_N_yield_log_intercept"].value
    a1 = compiled.parameters["biomass_N_yield_log_slope"].value
    f_n_at_actual = 1.0 / math.exp(a0 + a1 * actual_mgl)
    assert f_n_at_actual / compiled.param_values["biomass_N_fraction"] == pytest.approx(
        1.502, abs=5e-3
    ), "the size of the evaluation-point error moved off its characterized 1.50x"

    # ...and 362.7 is OUTSIDE Coleman's own fitted range, which is why summing the channels into
    # the fit is not a drop-in repair. Asserted so the next reader cannot take that route as
    # obviously correct.
    assert actual_mgl > 350.0, (
        "the commonest dose no longer leaves Coleman's fitted 70-350 mg N/L range; the "
        "'evaluate at total nitrogen' repair may now be admissible — re-open D-243"
    )


def test_all_assimilable_nitrogen_reaches_biomass_whatever_channel_it_entered_by():
    """The D-14 identity, generalised — and the reason the amino-acid dose IS extra biomass.

    D-32 designed the swap to be nitrogen-neutral *through the transfer*, and it is. But the pool
    it debits drains to ``N``, and ``N`` drains to biomass, so cumulatively the amino-acid pool is
    a second nitrogen reservoir rather than a substitute for the first. On a bare run the identity
    is exact; with a dose it holds to ~0.6 %, the remainder being pool residue the aging tail
    never re-assimilates.
    """
    for aa, tol in ((0.0, 1e-9), (0.5, 7e-3)):
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
    autolysis on and 45.9 % of it ends in the amino-acid pools instead: lees self-digestion
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
    assert in_pools / total0 == pytest.approx(0.459, abs=0.02), (
        f"autolysis now leaves {100 * in_pools / total0:.1f} % of the nitrogen in the amino-acid "
        "pools, not the characterized ~45.9 % — the autolysis/re-assimilation balance moved"
    )

    # The same run WITHOUT autolysis is the control: without it the 45.9 % is just a number, and
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
    assert ctl_pools / total0 < 0.01, (
        "the autolysis-off control now leaves nitrogen in the pools too, so the arm above no "
        f"longer attributes anything to autolysis ({100 * ctl_pools / total0:.2f} %)"
    )
