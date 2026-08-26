"""The census D-233 §10 named: parameters read at COMPILE that are ALSO in the sampled set.

**The predicate, stated exactly.** A parameter is a *census member* iff (1) it is read during
``compile_scenario`` — through ``ParameterSet.__getitem__``/``.value`` or off the ``dict`` that
``.resolve()`` returns — and (2) it is in the set ``simulate_ensemble`` would draw for that same
compiled scenario, i.e. ``_resolve_sample_names(process_set, parameters, None, None, events)``.

A member is only *half* sampled: an ensemble re-draws the parameter map but compiles ``y0`` and
``events`` once, at the nominal, so the compile-time role stays pinned while the runtime role
moves. Where the two roles were designed to cancel, sampling breaks the cancellation in one
direction and the surviving term can carry the wrong sign
(``feedback-a-parameter-can-be-pinned-and-drawn``, D-206).

**This is NOT the drawability surface** (D-153/D-156/D-157/D-159, ``test_drawability_surface.py``).
That one asks whether a banded parameter is reachable by the sampler *at all*, and its answer for
the compile-consumed class is "never drawn". This asks the opposite overlap: the names the sampler
*does* reach whose compile-time half it cannot move. The two sets are disjoint by construction.

**Why the instrument is a recorder and not a grep.** A ``grep`` for ``parameters["..."]`` over
``scenario/compile.py`` finds 26 names and misses the 19 ``pKa_*``, because the pH anchor reads
them off ``parameters.resolve()`` inside ``acidbase.build_pka_map``'s generator expression. That
gap is why the census was owed at all, so the guard here re-derives membership by recording actual
reads rather than by matching source text (``feedback-grep-finds-claims-not-guards``). The
pre-registered estimate for this census was 12-20 names; the recorder found **32** over the
battery (wine 29, beer 20, wine-with-overrides 30), and 21 of them — the 19 ``pKa_*`` and the two
``nitrogen_uptake_charge_*`` — are exactly the block a grep cannot see.

**The four verdicts, all measured.**

* ``REPAIRED`` — the 19 ``pKa_*`` plus ``nitrogen_uptake_charge_<medium>`` in their *t=0 anchor*
  role. ``CompiledScenario.reanchor_for_member`` re-solves the cation slot per member (D-233), and
  :func:`test_set_ph_anchors_only_the_nominal_member` measures the result: t=0 spread 2.03e-11 pH.
* ``LIVE`` — the same 19 ``pKa_*`` in their *second* compile role, ``_verb_set_ph``, whose event
  closes over the compile-time resolved map. Members land up to 0.0790 pH from the target they
  asked for and span 0.1320 (24 members, D-186's setting), against 2.03e-11 at t=0.
* ``LIVE`` — ``copper_typical``, which seeds the ``copper`` state slot at compile *and* is
  ``PhenolicBrowning``'s mean-centering reference at runtime, so ``f(Cu) == 1`` — the D-134 design
  invariant — holds only at the nominal draw.
* ``CLOSED`` — the eight ``must_aa_fraction_*``, measured across all eight at D-206 and reproduced
  here to the digit by :func:`test_the_amino_acid_gate_table_still_reproduces_d206`.
* ``BY-DESIGN`` — ``biomass_carrying_capacity`` and ``k_autolysis`` under a scenario override. The
  override is the *mode* of the draw, not a value the sampler discards, and D-164's in-band bound
  exists precisely so ``triangular(low, value, high)`` stays valid. D-24's surviving exclusion
  (scenario inputs are never sampled) is not breached: what is sampled is the parameter the input
  moved, carrying the YAML's uncertainty about it.

**What is pinned here and what deliberately is not.** The census *count* is not pinned — freezing
61-name-style totals is the vacuity D-159 refused, and it would go red on every new parameter. What
is pinned is that **every member is classified**: a future compile-time read of a sampled parameter
fails :func:`test_every_census_member_is_classified` by name, which is the only thing that stops
this set growing silently again. The two ``LIVE`` rows additionally pin their defects on purpose,
in the D-233 idiom — **a RED there means the defect was repaired: delete the guard and say so in
the record, do not revert the repair.**
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from contextlib import contextmanager

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.core.kinetics.amino_acid_pools import AMINO_ACID_SPECS
from fermentation.parameters.schema import Parameter
from fermentation.parameters.store import ParameterSet
from fermentation.runtime import simulate_scheduled
from fermentation.runtime.ensemble import _resolve_sample_names
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario import compile as compile_mod

# -- the recorder --------------------------------------------------------------------------
#
# Patched at CLASS level, not on one instance: the seam builds several derived ParameterSets
# (`_apply_nitrogen_dependent_yield`, `_override_carrying_capacity`, `_override_autolysis_rate`,
# `_inject_temperature_ramp_rate`), and a wrapper around the loaded one would stop recording at
# the first of them. Restored in a `finally`, and never used around an integration — it walks no
# stack, but it does allocate a recording dict per `resolve()` call.

_ORIG = (ParameterSet.__getitem__, ParameterSet.value, ParameterSet.resolve)


class _RecordingMap(dict):  # type: ignore[type-arg]
    """A resolved ``{name: float}`` map that records which keys were actually read."""

    seen: set[str]

    def __getitem__(self, key):
        self.seen.add(key)
        return dict.__getitem__(self, key)

    def get(self, key, default=None):
        self.seen.add(key)
        return dict.get(self, key, default)


@contextmanager
def _recording() -> Iterator[set[str]]:
    seen: set[str] = set()

    def rec_getitem(self, name):
        seen.add(name)
        return _ORIG[0](self, name)

    def rec_value(self, name):
        seen.add(name)
        return _ORIG[1](self, name)

    def rec_resolve(self, names=None):
        out = _ORIG[2](self, names)
        if names is not None:
            # An explicit name list IS the read — the caller asked for exactly these.
            seen.update(out)
            return out
        recording = _RecordingMap(out)
        recording.seen = seen
        return recording

    ParameterSet.__getitem__ = rec_getitem  # type: ignore[method-assign]
    ParameterSet.value = rec_value  # type: ignore[method-assign]
    ParameterSet.resolve = rec_resolve  # type: ignore[method-assign]
    try:
        yield seen
    finally:
        ParameterSet.__getitem__, ParameterSet.value, ParameterSet.resolve = _ORIG  # type: ignore[method-assign]


def _census(scenario: Scenario) -> tuple[set[str], set[str]]:
    """``(compile-time reads, census members)`` for one scenario."""
    with _recording() as seen:
        compiled = compile_scenario(scenario)
    sampled = set(
        _resolve_sample_names(
            compiled.process_set, compiled.parameters, None, None, compiled.events
        )
    )
    return seen, seen & sampled


# -- the battery ---------------------------------------------------------------------------
#
# One scenario per medium, each exercising every verb that medium admits, because a dynamic
# branch nobody takes is a silently dropped denominator [[feedback-count-and-print-your-skips]].
# `begin_aging` is load-bearing for `copper_typical`: `phenolic_browning` is the only Process that
# reads it, and it is disabled until then, so a ferment-only wine would drop that member for an
# uninteresting reason.

FERMENT_D = 10.0
ANCHOR_PH = 3.5
SET_PH_TARGET = 3.4

#: Where `set_ph` fires in the measurement scenario. The miss grows with how far the ferment has
#: dragged pH away from the anchor, because that is what gives the drawn pKas leverage: at day 6
#: it is 0.0147, at day 10 0.0420, at day 12 0.0582 (12 members each). Day 12 is the archive's
#: own setting - D-186 sized this defect there - and at 24 members it reproduces that record's
#: figure exactly (worst 0.07896, spread 0.13202), which is what makes the 12-member run below a
#: cheaper view of the same number rather than a different one.
SET_PH_DAY = 12.0

WINE = Scenario(
    name="census-wine",
    medium="wine",
    initial={
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.25,
        "amino_acids_gpl": 0.8,
        "anthocyanin_gpl": 0.5,
        "tannin_gpl": 1.5,
        "initial_ph": ANCHOR_PH,
        "tartaric_gpl": 6.0,
        "malic_gpl": 3.0,
        "so2_total_mgl": 30.0,
    },
    temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
    closure="screwcap",
    interventions=[
        Intervention(day=1.0, action="add_dap", params={"dap_gpl": 0.4}),
        Intervention(day=2.0, action="add_so2", params={"so2_mgl": 20.0}),
        Intervention(day=3.0, action="add_sugar", params={"sugar_gpl": 10.0}),
        Intervention(day=4.0, action="add_acid", params={"acid": "tartaric", "gpl": 1.0}),
        Intervention(day=5.0, action="add_copper", params={"copper_mgl": 0.5}),
        Intervention(day=FERMENT_D, action="begin_aging"),
        Intervention(day=FERMENT_D, action="set_ph", params={"ph": SET_PH_TARGET}),
        Intervention(day=FERMENT_D, action="add_oak", params={"oak_gpl": 4.0, "toast": "medium"}),
        Intervention(day=FERMENT_D + 1, action="add_oxygen", params={"o2_mgl": 5.0}),
        Intervention(day=FERMENT_D + 2, action="add_ascorbate", params={"ascorbate_mgl": 50.0}),
        Intervention(day=FERMENT_D + 3, action="rack", params={"fraction": 0.8}),
        Intervention(day=FERMENT_D + 4, action="seal_bottle", params={}),
    ],
    duration_days=20.0,
)

BEER = Scenario(
    name="census-beer",
    medium="beer",
    initial={
        "maltose_gpl": 90.0,
        "glucose_gpl": 15.0,
        "maltotriose_gpl": 20.0,
        "yan_mgl": 200.0,
        "pitch_gpl": 0.5,
        "initial_ph": 5.65,
    },
    temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
    interventions=[
        Intervention(day=1.0, action="add_dap", params={"dap_gpl": 0.2}),
        Intervention(day=3.0, action="add_acid", params={"acid": "lactic", "gpl": 0.5}),
        Intervention(day=4.0, action="set_ph", params={"ph": 5.4}),
        Intervention(day=5.0, action="add_oxygen", params={"o2_mgl": 8.0}),
        Intervention(day=6.0, action="rack", params={"fraction": 0.9}),
    ],
    duration_days=12.0,
)

#: A wine that takes the two scenario-override knobs, which no other battery member reaches
#: (`_override_carrying_capacity` / `_override_autolysis_rate` are opt-in on `scenario.initial`).
WINE_OVERRIDES = Scenario(
    name="census-wine-overrides",
    medium="wine",
    initial={
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.25,
        "amino_acids_gpl": 0.8,
        "initial_ph": ANCHOR_PH,
        "tartaric_gpl": 6.0,
        "malic_gpl": 3.0,
        "carrying_capacity_gpl": 4.0,
        "autolysis_rate_per_h": 5.0e-3,
    },
    temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
    duration_days=14.0,
)

#: Measurement scenarios, deliberately LEANER than the enumeration battery. Enumeration wants
#: every verb reached; a measurement wants the shortest run that still exhibits the effect, and
#: paying the battery's twelve-event twenty-day integration three times over would put ~50 s of
#: wall clock into one xdist chunk for no extra evidence.
WINE_SET_PH = Scenario(
    name="census-set-ph",
    medium="wine",
    initial={
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.25,
        "initial_ph": ANCHOR_PH,
        "tartaric_gpl": 6.0,
        "malic_gpl": 3.0,
    },
    temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
    interventions=[
        Intervention(day=SET_PH_DAY, action="begin_aging"),
        Intervention(day=SET_PH_DAY, action="set_ph", params={"ph": SET_PH_TARGET}),
    ],
    duration_days=SET_PH_DAY + 1.0,
)

#: `phenolic_browning` is the only Process that reads `copper_typical`, it is disabled until
#: `begin_aging`, and it returns zero without oxygen — so all three are required for the arms to
#: be anything but vacuously equal.
#:
#: **Do not add `add_copper` here.** Nothing else writes the `copper` slot, which is what makes
#: `copper - copper_typical` bit-zero at every step of the control arm and `f_copper` exactly
#: 1.0. A dose would break that equality for a reason that has nothing to do with the defect.
WINE_BROWNING = Scenario(
    name="census-browning",
    medium="wine",
    initial={
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.25,
        "anthocyanin_gpl": 0.5,
        "tannin_gpl": 1.5,
    },
    temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
    interventions=[
        Intervention(day=6.0, action="begin_aging"),
        Intervention(day=6.0, action="add_oxygen", params={"o2_mgl": 5.0}),
    ],
    duration_days=10.0,
)


BATTERY = (WINE, BEER, WINE_OVERRIDES)

#: The pH map, read at compile in TWO roles: the t=0 anchor back-solve (repaired per member by
#: D-233's `reanchor_for_member`) and `_verb_set_ph`'s event, which still closes over the
#: compile-time nominal map. Same names, two verdicts — which is why the registry is keyed by
#: name and the verdict prose names the role.
_PKA_SPECS = (*acidbase.ALL_ACIDS.values(), acidbase.BYP_AS_SUCCINIC, acidbase.CARBONIC_AS_CO2)
_PKA_MEMBERS = tuple(sorted({name for spec in _PKA_SPECS for name in spec.pka_param_names}))

#: name → verdict. The guard is that this covers every member; the prose is why each is here.
CENSUS: Mapping[str, str] = {
    **dict.fromkeys(
        _PKA_MEMBERS,
        "REPAIRED at t=0 by reanchor_for_member (D-233); LIVE in _verb_set_ph's event, which "
        "closes over the compile-time resolved map",
    ),
    **dict.fromkeys(
        ("nitrogen_uptake_charge_wine", "nitrogen_uptake_charge_beer"),
        "REPAIRED — cation_charge_for_ph reads it from the member's map",
    ),
    **{
        spec.fraction_param: "CLOSED at D-206 — the gate's f_i cancels; all eight measured there"
        for spec in AMINO_ACID_SPECS
    },
    "copper_typical": "LIVE — seeds the `copper` slot at compile, mean-centers f(Cu) at runtime",
    "biomass_carrying_capacity": "BY-DESIGN — a scenario override is the MODE of the draw (D-164)",
    "k_autolysis": "BY-DESIGN — a scenario override is the MODE of the draw (D-164)",
}


# -- the instrument's own control ------------------------------------------------------------


def test_the_recorder_sees_a_read_no_grep_can():
    """Positive control: without it an empty ``seen`` set passes every census test vacuously.

    ``pKa_tartaric_1`` is never subscripted by name anywhere in ``scenario/compile.py`` — it is
    read inside ``acidbase.build_pka_map``'s generator expression off ``parameters.resolve()``.
    A recorder that only caught ``ParameterSet.__getitem__`` would miss it, and would then report
    a census of 13 rather than 32 while going green.
    """
    seen, _ = _census(WINE)
    assert "pKa_tartaric_1" in seen, "the recorder cannot see resolve()-mediated reads"
    assert "must_aa_fraction_methionine" in seen, "the recorder cannot see __getitem__ reads"
    # ...and it must not report the whole file: `resolve()` hands back every name, so a recorder
    # that logged the map rather than the keys read off it would call all ~280 of them members.
    assert len(seen) < len(compile_scenario(WINE).parameters.names) / 2


def test_the_recorder_restores_the_parameter_set():
    """The patch is class-level; a leaked patch would silently slow every later test."""
    with _recording():
        pass
    assert (ParameterSet.__getitem__, ParameterSet.value, ParameterSet.resolve) == _ORIG


# -- the guard ---------------------------------------------------------------------------------


@pytest.mark.parametrize("scenario", BATTERY, ids=lambda s: s.name)
def test_every_census_member_is_classified(scenario):
    """A compile-time read of a sampled parameter must be a KNOWN one.

    This is the whole point of the beat: ``must_aa_fraction_methionine`` (D-206) and the pH anchor
    (D-233) were two members of a set nobody had enumerated, each found by accident while chasing
    something else. A new one appears here by name, with the verdicts of its neighbours to compare
    against, instead of waiting for a third accident.
    """
    _, members = _census(scenario)
    unclassified = sorted(members - set(CENSUS))
    assert not unclassified, (
        f"{scenario.name}: {len(unclassified)} parameter(s) are read at compile AND sampled but "
        f"carry no verdict: {unclassified}. Each is half-sampled — the compile-time role stays at "
        "the nominal while the runtime role moves. Measure both arms (recompile vs patch-the-map) "
        "and add the verdict here."
    )


def test_no_classified_member_has_gone_stale():
    """The registry must not accumulate names that stopped being members.

    A rotted registry is the failure mode of every frozen list in this repo: it keeps passing
    while describing a tree that no longer exists. Scored over the UNION of the battery, because
    membership is scenario-conditional — ``copper_typical`` is only a member once ``begin_aging``
    enables ``phenolic_browning``, and the two override knobs only once a scenario opts in.
    """
    union: set[str] = set()
    for scenario in BATTERY:
        union |= _census(scenario)[1]
    stale = sorted(set(CENSUS) - union)
    assert not stale, (
        f"{len(stale)} classified name(s) are no longer read at compile, or no longer sampled: "
        f"{stale}. The battery stopped reaching them, or the seam stopped reading them — check "
        "which before deleting the row. Note a REPAIR does not land here: a repaired member is "
        "still compile-read and still sampled, only its verdict changes, and the two defect "
        "pins are what go red."
    )


def test_the_census_and_the_drawability_surface_are_disjoint():
    """A name cannot be both 'never drawn' (D-159) and 'drawn but half-pinned' (this census).

    The two audits share a predicate half, so a reader can easily take one for the other. The
    oak yields, ``otr_*`` and ``bottling_burst_*`` are compile-consumed AND unreachable by the
    sampler; every census member is compile-consumed and reachable. Pinning the disjointness is
    what keeps the distinction from collapsing the next time either set is edited.

    Scored PER SCENARIO, because that is the scope of the claim: membership is
    scenario-conditional, so a name that is compile-read here and sampled only in some *other*
    battery member is not a counterexample to anything — comparing this scenario's compile-only
    set against the battery-wide registry would false-fail on exactly that case.
    """
    for scenario in BATTERY:
        seen, members = _census(scenario)
        compile_only = seen - members
        assert not (compile_only & members)
    seen, members = _census(WINE)
    assert "otr_screwcap" in (seen - members)
    assert "oak_yield_vanillin_medium" in (seen - members)


# -- the two LIVE rows, pinned as defects ------------------------------------------------------


def test_set_ph_anchors_only_the_nominal_member():
    """PINS A DEFECT ON PURPOSE — a RED here means ``set_ph`` was made per-member.

    **If this goes red, delete it and say so in the record; do not revert the repair.**

    ``initial_ph`` and ``set_ph`` make the same promise — the beverage sits at the pH you named —
    and since D-233 only the first keeps it for a sampled member. D-186's own docstring forbade
    exactly this state: *"Do not 'fix' this one anchor alone — the two would then disagree about
    what a member's pH means, and the state-level one would be the odd anchor out."* D-233 fixed
    the other one. So the configuration D-186 named as wrong is the one that ships, and it got
    there by a beat that repaired the half it was looking at.

    Measured on ONE ensemble so the two numbers are a comparison and not two records' figures:
    at t=0 the members span 2.03e-11 pH about the anchor; 0.05 h after the ``set_ph`` event they
    span 0.111 and miss the target by up to 0.058 (12 members; 0.132 / 0.079 at 24, which is
    D-186's own setting). The nominal run is exact at both, which is why nothing else in the
    suite sees this.
    """
    compiled = compile_scenario(WINE_SET_PH)
    probe_h = 24.0 * SET_PH_DAY + 0.05
    ens = compiled.run_ensemble(n_members=12, seed=0, t_eval=np.array([0.0, probe_h]))
    # The driver inserts EVERY breakpoint into the returned grid, so the requested index is not
    # the returned index; look the time up rather than counting
    # [[feedback-read-a-fast-curve-on-a-fixed-grid]].
    post = int(np.flatnonzero(np.isclose(ens.t, probe_h))[0])

    def ph(member: int, t_index: int) -> float:
        return acidbase.ph_of_state(
            ens.members[member][:, t_index], ens.schema, ens.member_params[member]
        )

    n = ens.members.shape[0]
    assert n >= 2, "need at least two members for a spread"
    at_t0 = np.array([ph(i, 0) for i in range(n)])
    after = np.array([ph(i, post) for i in range(n)])

    # The repaired anchor, as the baseline that makes the second number mean something.
    assert np.abs(at_t0 - ANCHOR_PH).max() < 1e-9

    # The unrepaired one. Pinned as a floor, not a band: the exact spread depends on the member
    # count and the seed, and what is being asserted is that it is nowhere near the anchor's.
    assert np.abs(after - SET_PH_TARGET).max() > 0.02, (
        "set_ph now reproduces its target per member — the D-186 asymmetry is closed; delete "
        "this guard and record it"
    )
    assert (after.max() - after.min()) > 1e6 * (at_t0.max() - at_t0.min())

    # The nominal run must be exact at BOTH, or this is a different defect than the one claimed.
    nominal = compiled.param_values
    assert acidbase.ph_of_state(ens.nominal[:, 0], ens.schema, nominal) == pytest.approx(
        ANCHOR_PH, abs=1e-9
    )
    # 1e-6 rather than 1e-9 only because the probe sits 0.05 h AFTER the event and the wine
    # keeps fermenting through those three minutes; the jump itself is exact.
    assert acidbase.ph_of_state(ens.nominal[:, post], ens.schema, nominal) == pytest.approx(
        SET_PH_TARGET, abs=1e-6
    )


def _with_value(params: ParameterSet, name: str, value: float) -> ParameterSet:
    base = params[name]
    return params.merge(
        ParameterSet(
            [
                Parameter(
                    name=name,
                    value=value,
                    unit=base.unit,
                    tier=base.tier,
                    uncertainty=base.uncertainty,
                    provenance=base.provenance,
                )
            ]
        ),
        override=True,
    )


def _compiled_with(scenario: Scenario, name: str | None = None, value: float | None = None):
    """Compile ``scenario``, optionally with ``name`` moved BEFORE the seam reads it.

    The control arm of the D-206 pattern: a recompile moves every role of the parameter at once,
    where patching the resolved map moves only the runtime role — the same channel, one step
    earlier [[feedback-a-control-needs-mechanical-reach]].
    """
    if name is None:
        return compile_scenario(scenario)
    assert value is not None, "naming a parameter without a value moves nothing"
    moved = float(value)
    original = compile_mod._load_parameters

    def patched(sc, parameter_paths, data_dir):
        return _with_value(original(sc, parameter_paths, data_dir), name, moved)

    compile_mod._load_parameters = patched
    try:
        compiled = compile_scenario(scenario)
    finally:
        compile_mod._load_parameters = original
    assert compiled.parameters[name].value == moved, "the recompile did not move the parameter"
    return compiled


def _endpoint(compiled, slot: str, patch: dict[str, float] | None = None) -> float:
    params = compiled.param_values if patch is None else {**compiled.param_values, **patch}
    out = simulate_scheduled(
        compiled.process_set, params, compiled.y0, compiled.t_span_h, events=compiled.events
    )
    return float(np.atleast_1d(out.y[compiled.schema.slice(slot), -1]).sum())


def test_copper_typical_reaches_the_run_only_through_the_broken_cancellation():
    """PINS A DEFECT ON PURPOSE — a RED here means the ``copper`` seed was made per-member.

    **If this goes red, delete it and say so in the record; do not revert the repair.**

    ``copper_typical`` has two roles that D-134 made numerically identical on purpose: it seeds
    the ``copper`` state slot at the compile seam (for a wine that does not name ``copper_gpl``),
    and it is the mean-centering reference in ``PhenolicBrowning``'s multiplier
    ``f(Cu) = 1 + k·(copper − copper_typical)``, which is therefore exactly 1 for an
    un-overridden wine. Moving it coherently keeps ``f == 1``; a drawn member moves only the
    reference, so the invariant holds at the nominal draw and nowhere else.

    The consequence is the D-206 shape on a second member, and it is total rather than partial:
    the CONTROL arm is bit-identical to the baseline on ``A420``, so a coherent change to this
    parameter contributes **exactly nothing** to the reported band, while the sampler arm moves
    aged browning by ~14 % at the band top. Every bit of ``copper_typical``'s share of the spread
    is the broken cancellation.
    """
    baseline = _compiled_with(WINE_BROWNING)
    p = baseline.parameters["copper_typical"]
    seed = float(np.atleast_1d(baseline.y0[baseline.schema.slice("copper")])[0])
    assert seed == pytest.approx(p.value), (
        "the `copper` slot is no longer seeded from the reference"
    )

    a420_base = _endpoint(_compiled_with(WINE_BROWNING), "A420")
    a420_sampler = _endpoint(
        _compiled_with(WINE_BROWNING), "A420", patch={"copper_typical": p.uncertainty.high}
    )
    a420_control = _endpoint(
        _compiled_with(WINE_BROWNING, "copper_typical", p.uncertainty.high), "A420"
    )

    assert a420_base > 0.0, "no browning in this scenario — the arms would be vacuously equal"
    # The honest channel is silent. Exact, not toleranced: identical parameters into an identical
    # RHS give bit-identical BDF output, and f(Cu) is 1 in both runs by construction.
    assert a420_control == a420_base, (
        "a coherent copper_typical now moves browning — the mean-centering is no longer exact; "
        "re-derive this guard rather than loosening it"
    )
    # The sampler channel is not.
    assert abs(a420_sampler - a420_base) / a420_base > 0.05, (
        "the sampler no longer breaks the f(Cu) cancellation — the seed is per-member; delete "
        "this guard and record it"
    )


def test_the_arms_harness_can_see_a_difference_it_is_asked_for():
    """The positive control the copper NULL owes.

    ``feedback-a-null-result-needs-a-positive-control``.

    ``a420_control == a420_base`` above is an equality, and an equality is satisfiable by a
    harness that perturbed nothing. ``k_browning_phenolic`` is read only at runtime, so BOTH arms
    must move it and both must move it by the same amount — which is the strongest form of the
    control, because it also shows the recompile path and the patch path agree when there is no
    compile-time role to pin.
    """
    baseline = _compiled_with(WINE_BROWNING)
    p = baseline.parameters["k_browning_phenolic"]
    a420_base = _endpoint(_compiled_with(WINE_BROWNING), "A420")
    a420_sampler = _endpoint(
        _compiled_with(WINE_BROWNING), "A420", patch={"k_browning_phenolic": p.uncertainty.high}
    )
    a420_control = _endpoint(
        _compiled_with(WINE_BROWNING, "k_browning_phenolic", p.uncertainty.high), "A420"
    )

    assert abs(a420_sampler - a420_base) / a420_base > 0.05, "the harness perturbed nothing"
    assert a420_control == a420_sampler, "a runtime-only parameter must give identical arms"


# -- the two verdicts that are NOT defects, measured rather than asserted -----------------------


def test_the_amino_acid_gate_table_still_reproduces_d206():
    """The eight fractions are CLOSED, and this checks the record that closed them is current.

    D-206 enumerated all eight from ``AMINO_ACID_SPECS`` rather than arguing from methionine, and
    reported the gate each takes under a drawn fraction. Bands move between beats, so the census's
    job here is not to re-measure but to confirm the archive's table still describes what ships —
    ``feedback-reproduce-a-published-number-before-trusting-the-new-column``, applied to the
    repo's own record. Closed form: the pool is seeded ``dose·f_i/Σf`` at compile and the gate is
    ``aa_i/(K·f_i + aa_i)`` at runtime, so the fraction cancels exactly when both halves move.
    """
    compiled = compile_scenario(WINE)
    params = compiled.parameters
    k = params["K_amino_acids"].value
    dose = WINE.initial["amino_acids_gpl"]
    sum_f = sum(params[spec.fraction_param].value for spec in AMINO_ACID_SPECS)

    def gate(seeded: float, drawn: float) -> float:
        return seeded / (k * drawn + seeded)

    spans = {}
    for spec in AMINO_ACID_SPECS:
        p = params[spec.fraction_param]
        seeded = dose * p.value / sum_f  # pinned at the nominal by the compile seam
        # D-206: "every gate is 0.888889 at nominal — identical by construction".
        assert gate(seeded, p.value) == pytest.approx(0.888889, abs=5e-6)
        low, high = gate(seeded, p.uncertainty.low), gate(seeded, p.uncertainty.high)
        spans[spec.fraction_param] = (low, high)

    # D-206's two named extremes, to four decimals as it printed them.
    assert spans["must_aa_fraction_methionine"] == pytest.approx((0.9524, 0.8333), abs=5e-5)
    assert spans["must_aa_fraction_phenylalanine"] == pytest.approx((0.8924, 0.8854), abs=5e-5)
    widest = max(spans, key=lambda n: spans[n][0] - spans[n][1])
    narrowest = min(spans, key=lambda n: spans[n][0] - spans[n][1])
    assert widest == "must_aa_fraction_methionine"
    assert narrowest == "must_aa_fraction_phenylalanine"


def test_a_scenario_override_is_the_mode_of_the_draw_not_a_discarded_input():
    """``carrying_capacity_gpl`` / ``autolysis_rate_per_h`` are BY-DESIGN members, not defects.

    Both are read at compile and both are sampled, which is the census predicate — but the read
    builds a Parameter whose *value* is the scenario's number and whose *uncertainty* is the YAML
    reference's band, and the sampler draws ``triangular(low, value, high)``. So the stated
    operating point is the mode of the member distribution rather than a number the sampler throws
    away, which is also why D-164 requires the override to land inside the band.

    D-24's surviving exclusion is therefore intact: no *scenario input* is sampled here. What is
    sampled is the parameter the input moved, carrying the YAML's uncertainty about it — which is
    what an ensemble is for. Asserted on the compiled ParameterSet rather than on a run, because
    the claim is about the distribution the sampler is handed, not about a trajectory.
    """
    plain = compile_scenario(
        WINE_OVERRIDES.model_copy(
            update={
                "initial": {
                    k: v
                    for k, v in WINE_OVERRIDES.initial.items()
                    if k not in ("carrying_capacity_gpl", "autolysis_rate_per_h")
                }
            }
        )
    )
    overridden = compile_scenario(WINE_OVERRIDES)

    for name, key in (
        ("biomass_carrying_capacity", "carrying_capacity_gpl"),
        ("k_autolysis", "autolysis_rate_per_h"),
    ):
        asked = WINE_OVERRIDES.initial[key]
        assert overridden.parameters[name].value == asked, "the override is not the mode"
        assert overridden.parameters[name].value != plain.parameters[name].value
        # The band is the YAML reference's, unchanged — that is what makes the mode meaningful
        # and what D-164's in-band bound protects.
        assert overridden.parameters[name].uncertainty.low == plain.parameters[name].uncertainty.low
        assert (
            overridden.parameters[name].uncertainty.high == plain.parameters[name].uncertainty.high
        )
        assert (
            plain.parameters[name].uncertainty.low
            <= asked
            <= plain.parameters[name].uncertainty.high
        )
