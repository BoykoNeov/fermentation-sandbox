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
  role. ``CompiledScenario.reanchor_for_member`` re-solves the cation slot per member (D-233):
  t=0 spread 2.03e-11 pH.
* ``REPAIRED`` — ``pKa_peptide_buffer`` in its **third** compile role (D-238), the one D-233 named
  and declined: it is the pKa the beer peptide capacity was back-solved at OFFLINE, so a member
  carried a wort at BC 1.1161-1.180 rather than the 1.18 the constant exists to reproduce.
  ``y0_for_member``'s rule 3 re-roots the capacity on the member's own map and seeds it before the
  anchor reads it. The cost is the mirror image of the anchor role's and the pair is worth reading
  together: the anchor half is 100 % artefact at t=0 and washes out, this half is **exactly zero**
  at t=0 (the cation back-solve absorbs it) and grows to 0.0095 pH by day 14, converged.
* ``REPAIRED`` — the same 19 ``pKa_*`` in their *second* compile role, ``_verb_set_ph``. Its event
  closed over the compile-time resolved map until D-235 widened ``StateMutation`` to hand every
  mutation the **running** map; members landed up to 0.07896 pH from the target they asked for and
  spanned 0.13202 (24 members, D-186's setting), and now sit 1.98e-11 / 2.13e-11 — the t=0 anchor's
  own class. ``add_dap``'s z̄ read moved with it (worst member 0.0108 pH at the dose).
* ``REPAIRED`` — ``copper_typical``, which seeds the ``copper`` state slot at compile *and* is
  ``PhenolicBrowning``'s mean-centering reference at runtime, so ``f(Cu) == 1`` — the D-134 design
  invariant — held only at the nominal draw. D-236 re-seeds the slot per member, taking aged
  ``A420``'s spread from 16.65 % (all of it artefact: the coherent channel measured exactly zero)
  to bit-identical. A wine that NAMES ``copper_gpl`` is deliberately left alone — there the
  reference and the wine's copper are independent and drawing one is correct.
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
this set growing silently again. **There are no ``LIVE`` rows left**: D-235 repaired the pH pair
and D-236 the copper seed, and the two defect pins written in D-233's idiom were resolved in the
two ways that idiom allows. ``set_ph``'s went RED and was deleted for the positive form
(:func:`test_set_ph_reproduces_its_target_for_every_member`); copper's stayed GREEN, because its
arms drive ``simulate_scheduled`` by hand and never reach ``y0_for_member`` — so it was
**re-scoped**, not deleted, to pin the mechanism while
:func:`test_the_copper_seed_moves_with_the_member` pins the repair. A guard that keeps passing
across a repair is not evidence the repair was inert; read what path it drives.

**The census is WIRING-INVARIANT (decision D-237).** D-234 scored it over the ``direct`` oxidative
set only. Re-run against ``cascade`` and ``direct_burst`` on every battery member it is identical
in both directions, because the three sets differ in which *Processes* are wired and a Process is a
runtime consumer — compile-time reads come from the seam's seeding and the verbs, which do not
branch on the oxidative set. The copper repair reaches every wiring bit-identically, which D-236 §5
could only argue: PRE-repair 16.65 % (direct) / 15.38 % (cascade) / 16.73 % (burst), all to
0.0000 %. The cascade is the *mildest* arm, not the widest as D-234 §9 expected — it moves twice
the O₂ and less browning, because it multiplies one shared gate that every sink downstream splits.

**The one LIVE row here is not a census member at all** — it is the census's *complement*, found by
that widening: ``burst_antioxidant_initial`` is banded 50×, seeds a live pool, and no Process
declares it, so the sampler can never reach it
(:func:`test_the_burst_seed_carries_a_fifty_fold_band_the_ensemble_never_draws`). Three cells, and
only two of them had a guard before D-237: compile-read AND sampled (here), banded AND drawn but
inert (``test_drawability_surface.py``), banded and live and **never drawn** (nothing, until now).
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


#: The same wine, but STATING its copper. `copper_gpl` makes the wine's copper an independent
#: quantity from the reference `f(Cu)` is centred on, so drawing the reference alone is CORRECT
#: here and the seed must stay put — the arm that keeps D-236's repair from overwriting a
#: scenario input and breaching D-24's exclusion.
WINE_BROWNING_NAMED = Scenario(
    name="census-browning-named-copper",
    medium="wine",
    initial={**WINE_BROWNING.initial, "copper_gpl": 5.0e-4},
    temperature_schedule=WINE_BROWNING.temperature_schedule,
    interventions=WINE_BROWNING.interventions,
    duration_days=WINE_BROWNING.duration_days,
)


BATTERY = (WINE, BEER, WINE_OVERRIDES)

#: The pH map, read at compile in TWO roles: the t=0 anchor back-solve (repaired per member by
#: D-233's `reanchor_for_member`) and `_verb_set_ph`'s event (repaired by D-235, which hands the
#: mutation the running map). Both roles now move with the member; the registry stays keyed by
#: name with the role named in the prose, because that is what made the split legible when the
#: two verdicts differed.
_PKA_SPECS = (*acidbase.ALL_ACIDS.values(), acidbase.BYP_AS_SUCCINIC, acidbase.CARBONIC_AS_CO2)
_PKA_MEMBERS = tuple(sorted({name for spec in _PKA_SPECS for name in spec.pka_param_names}))

#: name → verdict. The guard is that this covers every member; the prose is why each is here.
CENSUS: Mapping[str, str] = {
    **dict.fromkeys(
        _PKA_MEMBERS,
        "REPAIRED in BOTH compile roles — the t=0 anchor by reanchor_for_member (D-233), and "
        "_verb_set_ph's event by the running map it is handed at the breakpoint (D-235)",
    ),
    **dict.fromkeys(
        ("nitrogen_uptake_charge_wine", "nitrogen_uptake_charge_beer"),
        "REPAIRED — cation_charge_for_ph reads it from the member's map, and since D-235 so do "
        "_verb_set_ph's and _verb_add_dap's mutations",
    ),
    **{
        spec.fraction_param: "CLOSED at D-206 — the gate's f_i cancels; all eight measured there"
        for spec in AMINO_ACID_SPECS
    },
    "copper_typical": (
        "REPAIRED — y0_for_member re-seeds the `copper` slot per member (D-236), restoring "
        "f(Cu) == 1; a scenario that NAMES copper_gpl is deliberately left alone"
    ),
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


# -- the pH pair (REPAIRED at D-235) and the one surviving LIVE row --------------------------


def test_set_ph_reproduces_its_target_for_every_member():
    """The repaired half of the pH-anchor pair (decision D-235), in its positive form.

    **This replaces a guard that pinned the defect on purpose.** Until D-235 ``_verb_set_ph``'s
    event closed over the compile-time resolved map, so a sampled member re-anchored with the
    nominal's pKa set and landed up to 0.07896 pH from the target it asked for (spread 0.13202,
    24 members, D-186's own setting) — against 2.03e-11 at t=0, which D-233 had just repaired.
    D-186's docstring had named that exact configuration as the thing not to produce. The
    mutation now reads the **running** map ``simulate_scheduled`` hands it, which under an
    ensemble is the member's own draw, so both anchors are per-member and agree again.

    Measured AT the event, not after it. The driver emits every breakpoint post-mutation, so the
    grid point at the event time IS the anchored state; the 0.05 h probe the old guard used
    carries three minutes of ferment on top and reads 7.04e-07 for that reason alone, not because
    the anchor is loose [[feedback-read-a-fast-curve-on-a-fixed-grid]].
    """
    compiled = compile_scenario(WINE_SET_PH)
    event_h = 24.0 * SET_PH_DAY
    ens = compiled.run_ensemble(n_members=12, seed=0, t_eval=np.array([0.0, event_h]))
    at_event = int(np.flatnonzero(np.isclose(ens.t, event_h))[0])

    def ph(member: int, t_index: int) -> float:
        return acidbase.ph_of_state(
            ens.members[member][:, t_index], ens.schema, ens.member_params[member]
        )

    n = ens.members.shape[0]
    assert n >= 2, "need at least two members for a spread"
    at_t0 = np.array([ph(i, 0) for i in range(n)])
    after = np.array([ph(i, at_event) for i in range(n)])

    # The t=0 anchor (D-233), as the baseline that says what "per-member" costs here.
    assert np.abs(at_t0 - ANCHOR_PH).max() < 1e-9

    # The second anchor, now in the same class. Pinned as a CEILING three orders of magnitude
    # below the 0.058 this scenario used to miss by at 12 members, and six above the root
    # finder's own residual — a band, not the digits [[feedback-pin-the-band-not-the-nominal]].
    assert np.abs(after - SET_PH_TARGET).max() < 1e-9, (
        "set_ph no longer reproduces its target per member — the D-186 pair has come apart again"
    )
    assert (after.max() - after.min()) < 1e-9

    # The nominal run stays exact at both, which is what says the widening moved nothing that was
    # already right: `run()` passes `param_values`, so the mutation reads the same numbers it read
    # when it closed over `resolve()`.
    nominal = compiled.param_values
    assert acidbase.ph_of_state(ens.nominal[:, 0], ens.schema, nominal) == pytest.approx(
        ANCHOR_PH, abs=1e-9
    )
    assert acidbase.ph_of_state(ens.nominal[:, at_event], ens.schema, nominal) == pytest.approx(
        SET_PH_TARGET, abs=1e-9
    )


def test_exactly_two_verbs_read_the_running_parameter_map():
    """The blast-radius claim D-235 had to prove, pinned rather than asserted.

    Widening :data:`~fermentation.runtime.schedule.StateMutation` handed every verb's mutation a
    third argument. Two of the thirteen read it — ``set_ph`` (the pH anchor) and ``add_dap`` (the
    z̄ a nitrogen pool's charge excess is measured relative to). The other eleven are byte-identical
    with or without it, because their jump is a scenario input in grams and not a parameter.

    Recorded rather than diffed: a numerical comparison under a perturbed map would call a verb a
    reader only where the read happens to MOVE something, and would call ``dosed_charge``'s
    compile-time read a non-read for the same reason. What is asserted is which names each
    mutation actually looks up [[feedback-grep-finds-claims-not-guards]].
    """
    compiled = compile_scenario(WINE)
    readers: dict[str, list[str]] = {}
    n_mutations = 0
    for event in compiled.events:
        if event.mutate is None:
            continue
        n_mutations += 1
        recording = _RecordingMap(compiled.param_values)
        recording.seen = set()
        event.mutate(compiled.schema, compiled.y0.copy(), recording)
        if recording.seen:
            readers[event.label.split("@")[0]] = sorted(recording.seen)

    assert n_mutations >= 11, f"the battery stopped reaching the verbs: {n_mutations} mutations"
    assert set(readers) == {"set_ph", "add_dap"}, (
        f"the set of verbs reading the running map changed: {sorted(readers)}. A NEW reader is "
        "not automatically wrong — it is a verb that became per-member — but it needs a measured "
        "record like D-235's, and a removed one means a member is back on the nominal."
    )
    # add_dap reads exactly the medium's z̄. `dap_nitrogen_charge` is deliberately NOT here: it is
    # a property of the compound, is in no Process's `reads`, and so has no member value.
    assert readers["add_dap"] == ["nitrogen_uptake_charge_wine"]
    # set_ph rebuilds the whole pKa map plus the CO2 solubility triple and the same z̄.
    assert set(_PKA_MEMBERS) <= set(readers["set_ph"])
    assert "nitrogen_uptake_charge_wine" in readers["set_ph"]


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


def test_a_hand_wired_run_on_a_patched_map_still_breaks_the_f_cu_cancellation():
    """The MECHANISM, which survives the repair — and is why the repair had to live where it does.

    **This was the defect pin, and it did NOT go red at D-236**, which is the thing to read
    carefully rather than as "the repair was inert". The arms below drive
    :func:`~fermentation.runtime.simulate_scheduled` directly on ``compiled.y0`` with a patched
    parameter map, so they never touch ``y0_for_member`` — the hook D-236's repair is made of.
    Patching a map after compile still moves the reference and not the seed, exactly as it always
    did; what changed is that an ENSEMBLE no longer does that, because it rebuilds the seed first.
    So this test is re-scoped rather than deleted: it pins the mechanism that made the defect
    possible, and :func:`test_the_copper_seed_moves_with_the_member` pins the repair.

    ``copper_typical`` has two roles that D-134 made numerically identical on purpose: it seeds
    the ``copper`` state slot at the compile seam (for a wine that does not name ``copper_gpl``),
    and it is the mean-centering reference in ``PhenolicBrowning``'s multiplier
    ``f(Cu) = 1 + k·(copper − copper_typical)``, which is therefore exactly 1 for an
    un-overridden wine. Moving it coherently keeps ``f == 1``; a drawn member moves only the
    reference, so the invariant holds at the nominal draw and nowhere else.

    The consequence is the D-206 shape on a second member, and it is total rather than partial:
    the CONTROL arm is bit-identical to the baseline on ``A420``, so a coherent change to this
    parameter contributes **exactly nothing** to the reported band, while the patched arm moves
    aged browning by ~14 % at the band top. Before D-236 that patched arm was what an ensemble
    member actually experienced, and every bit of ``copper_typical``'s share of the reported
    spread was the broken cancellation (16.65 % across 12 members, with the coherent channel at
    exactly zero). It is now what a hand-wired caller experiences, and nothing else.
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
    # The patched-map channel is not — and must stay so, because it is the mechanism whose
    # reach into ensembles D-236 closed. A RED here means patching the map stopped moving
    # browning at all, which would make the repair unfalsifiable rather than unnecessary.
    assert abs(a420_sampler - a420_base) / a420_base > 0.05, (
        "patching `copper_typical` after compile no longer breaks the f(Cu) cancellation — the "
        "seed moved with it, and the arms harness can no longer see the mechanism at all"
    )


def _seed_and_a420(compiled, n_members: int = 12):
    """One ``only=['copper_typical']`` ensemble, reported as (per-member seed, per-member A420).

    ``only=`` is the whole point: with a single name drawn, every member's RHS is identical to
    the nominal's *wherever the mean-centring is intact*, so the A420 comparison below can be
    exact equality rather than a tolerance — the same idiom as the control arm above
    [[feedback-the-setting-where-a-change-is-exact-is-the-control]].
    """
    ens = compiled.run_ensemble(n_members=n_members, seed=0, only=["copper_typical"])
    assert ens.sampled_names == ("copper_typical",), f"drew {ens.sampled_names}"
    assert ens.n_succeeded >= 2, "need at least two members for a spread"
    cu, a420 = ens.schema.slice("copper"), ens.schema.slice("A420")
    seeds = [float(ens.members[i][cu, 0][0]) for i in range(ens.n_succeeded)]
    drawn = [float(ens.member_params[i]["copper_typical"]) for i in range(ens.n_succeeded)]
    endpoints = [float(ens.members[i][a420, -1][0]) for i in range(ens.n_succeeded)]
    return seeds, drawn, endpoints, float(ens.nominal[a420, -1][0])


def test_the_copper_seed_moves_with_the_member():
    """The repair (decision D-236): a drawn ``copper_typical`` moves BOTH of its roles.

    ``copper_typical`` seeds the ``copper`` slot at compile *and* is the reference
    ``PhenolicBrowning`` mean-centres on, and D-134 made them numerically identical so that
    ``f(Cu)`` is **exactly** 1 for an un-overridden wine. An ensemble used to move the reference
    alone, so that invariant held at the nominal draw and nowhere else, and aged ``A420`` moved
    **16.65 %** across 12 members — against a coherent channel of **exactly zero**, so all of it
    was the broken cancellation, with the sign the parameter's name argues against.
    ``CompiledScenario.y0_for_member`` now re-seeds the slot from the member's own draw.

    The claim is exact equality, not a tolerance: with only this one name drawn and both roles
    moving together, every member integrates the identical RHS as the nominal run, so the whole
    aged endpoint is bit-identical. A tolerance here would pass on a repair that merely shrank
    the artefact.
    """
    seeds, drawn, endpoints, nominal = _seed_and_a420(compile_scenario(WINE_BROWNING))
    assert seeds == drawn, (
        "a member's `copper` seed is not its own drawn `copper_typical` — the compile-time role "
        "is pinned again and f(Cu) is 1 only at the nominal draw"
    )
    assert len(set(seeds)) == len(seeds), "the draw produced no spread; the arms are vacuous"
    assert nominal > 0.0, "no browning in this scenario — the members would agree vacuously"
    assert endpoints == [nominal] * len(endpoints), (
        "aged A420 is not bit-identical across members. With `copper_typical` the only name "
        "drawn and both of its roles moving, f(Cu) is exactly 1 for every member, so any spread "
        "at all is a surviving half-pinned role — re-derive it, do not loosen this to a rel="
    )


def test_a_scenario_that_names_its_copper_is_not_re_seeded():
    """The other half of D-236's branch — and D-24's exclusion, kept intact.

    When a scenario states ``copper_gpl`` the wine's copper is a scenario INPUT and the reference
    ``f(Cu)`` is centred on is a separate, sourced quantity. Drawing the reference alone is then
    *correct*: a wine at 0.5 mg/L compared against a lower typical browns faster, and that is
    physics rather than an artefact. Re-seeding here would overwrite an input with a parameter
    draw, which is exactly what D-24 excludes and what D-234 §7 went out of its way to certify
    intact — so the repair is conditional, and both halves of the condition are measured.

    Guarded rather than trusted, because the branch is invisible from the repaired side: every
    assertion in :func:`test_the_copper_seed_moves_with_the_member` would still pass if the rule
    fired unconditionally.
    """
    compiled = compile_scenario(WINE_BROWNING_NAMED)
    stated = float(WINE_BROWNING_NAMED.initial["copper_gpl"])
    seeds, drawn, endpoints, nominal = _seed_and_a420(compiled)

    assert float(compiled.y0[compiled.schema.slice("copper")][0]) == stated
    assert seeds == [stated] * len(seeds), (
        "the per-member builder overwrote a scenario input with a parameter draw — D-24's "
        "exclusion (scenario inputs are never sampled) is breached"
    )
    assert len(set(drawn)) == len(drawn), "the draw produced no spread; the arm is vacuous"
    # ...and the reference's own channel is genuinely live here, which is why the seed must not
    # follow it. A member browns differently because its wine really does sit off the reference.
    assert max(abs(e - nominal) for e in endpoints) / nominal > 0.05, (
        "a drawn reference no longer moves browning for a wine that stated its copper — the "
        "mean-centring has lost its physics, not just its artefact"
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


# -- the alternate oxidative wirings, measured rather than argued (decision D-237) --------------
#
# D-234 §9 named the `oxidative="cascade"` / `"direct_burst"` sets an unmeasured third battery
# member, "expected to widen §5's number". D-236 §5 then argued the copper repair reaches them for
# free — it moves the SEED, and every consumer mean-centres on the same reference — and said in as
# many words that this was "a mechanical argument, not a measurement". These are the measurement.

#: The two non-default wirings. `direct` is not here because the tests above already are it, and
#: beer is not here because `_OXIDATIVE_BURST_BEER == _OXIDATIVE_DIRECT_BEER` — a beer burst arm
#: would be a duplicate rather than a third number, asserted below rather than assumed.
_ALT_WIRINGS = ("cascade", "direct_burst")


def _seed_and_a420_wired(scenario, oxidative: str, *, repaired: bool, n_members: int = 8):
    """:func:`_seed_and_a420`'s arm, under a chosen wiring and with the repair on or off.

    ``y0_for_member=None`` reproduces the pre-D-233 fixed-``y0`` path exactly (``run_ensemble``
    uses ``kwargs.setdefault``), so the PRE arm here is the shipped code before the repair and
    not a re-implementation of it — which is what lets it serve as this test's positive control
    per wiring [[feedback-a-null-result-needs-a-positive-control]].

    8 members rather than 12: the POST claim is exact equality, which does not get stronger with
    more draws, and the PRE control clears 10 % at 8. The measurement in the record used 12.
    """
    compiled = compile_scenario(scenario, oxidative=oxidative)
    kwargs = {} if repaired else {"y0_for_member": None}
    ens = compiled.run_ensemble(n_members=n_members, seed=0, only=["copper_typical"], **kwargs)
    assert ens.sampled_names == ("copper_typical",), f"drew {ens.sampled_names}"
    assert ens.n_succeeded >= 2, "need at least two members for a spread"
    cu, a420 = ens.schema.slice("copper"), ens.schema.slice("A420")
    seeds = [float(ens.members[i][cu, 0][0]) for i in range(ens.n_succeeded)]
    drawn = [float(ens.member_params[i]["copper_typical"]) for i in range(ens.n_succeeded)]
    endpoints = [float(ens.members[i][a420, -1][0]) for i in range(ens.n_succeeded)]
    return seeds, drawn, endpoints, float(ens.nominal[a420, -1][0])


def test_beers_burst_set_is_beers_direct_set_so_it_is_not_a_third_number():
    """Why the wiring sweep is wine-only — the identity is asserted, not assumed.

    ``AntioxidantBurstOxidation`` is wine-only (Ferreira 2015's dataset is red wine and beer
    carries no ``burst_antioxidant`` slot), so ``get_medium("beer", oxidative="direct_burst")``
    is the direct beer build unchanged. Counting a beer burst arm as a third measurement would be
    reporting the same number twice — the denominator hazard
    [[feedback-count-and-print-your-skips]] in its quietest form, because nothing skips.
    """
    from fermentation.core.media import _OXIDATIVE_BURST_BEER, _OXIDATIVE_DIRECT_BEER

    assert _OXIDATIVE_BURST_BEER == _OXIDATIVE_DIRECT_BEER, (
        "beer's burst set has diverged from its direct set — a beer `direct_burst` arm is now a "
        "genuinely third wiring and this sweep no longer covers it"
    )


@pytest.mark.parametrize("oxidative", _ALT_WIRINGS)
def test_the_copper_repair_holds_under_every_oxidative_wiring(oxidative):
    """D-236 §5's mechanical argument, turned into a measurement — and it holds.

    Under every wiring ``copper_typical`` reaches the run only through a multiplier mean-centred
    on itself: ``PhenolicBrowning``'s ``f(Cu)`` under the direct sets, and
    :func:`~fermentation.core.kinetics.oxidative_cascade.activation_rate`'s re-homed copy of the
    same D-134 term under the cascade (where ``PhenolicBrowning`` is absent and ``A420`` comes
    from ``QuinonePolymerization`` instead). So re-seeding the slot restores ``f(Cu) == 1`` per
    member whatever is wired downstream, and the aged endpoint is bit-identical.

    The PRE arm is the control and it is not decoration: it is what makes the POST equality a
    null result rather than a harness that perturbed nothing. Measured at 12 members —
    cascade 15.38 %, direct_burst 16.73 % — both to 0.0000 %.
    """
    pre_seeds, drawn, pre_ends, pre_nominal = _seed_and_a420_wired(
        WINE_BROWNING, oxidative, repaired=False
    )
    assert len(set(drawn)) == len(drawn), "the draw produced no spread; the arms are vacuous"
    assert pre_nominal > 0.0, f"no browning under {oxidative} — the members would agree vacuously"
    assert len(set(pre_seeds)) == 1, "the PRE arm is meant to pin the seed at the nominal"
    pre_worst = max(abs(e - pre_nominal) for e in pre_ends) / pre_nominal
    assert pre_worst > 0.10, (
        f"under {oxidative} the un-repaired path moves aged A420 by only {pre_worst:.2%}; this "
        "control is what stops the POST equality below from being satisfied by a dead harness"
    )

    seeds, drawn, endpoints, nominal = _seed_and_a420_wired(WINE_BROWNING, oxidative, repaired=True)
    assert seeds == drawn, (
        f"under {oxidative} a member's `copper` seed is not its own drawn `copper_typical`"
    )
    assert nominal == pre_nominal, "the nominal member must not depend on the repair"
    assert endpoints == [nominal] * len(endpoints), (
        f"aged A420 is not bit-identical across members under {oxidative}. D-236 §5 argued the "
        "seed repair reaches every wiring because they all mean-centre on the same reference; a "
        "RED here means some consumer reads `copper_typical` OUTSIDE that difference, and the "
        "argument is false for this wiring — measure it, do not loosen this to a rel="
    )


@pytest.mark.parametrize("oxidative", _ALT_WIRINGS)
def test_a_named_copper_survives_the_repair_under_every_oxidative_wiring(oxidative):
    """The honest channel, checked under the alternate wirings too.

    A wine that states ``copper_gpl`` must keep its seed and must keep browning differently as
    the reference moves — that is physics, and D-24's exclusion says a parameter draw may never
    overwrite it. The condition is invisible from the repaired side, so it owes its own arm under
    every wiring the repaired side is claimed for. Measured at 12 members: cascade 13.05 %,
    direct_burst 14.42 %, and the PRE/POST numbers are identical because the rule does not fire.
    """
    stated = float(WINE_BROWNING_NAMED.initial["copper_gpl"])
    seeds, drawn, endpoints, nominal = _seed_and_a420_wired(
        WINE_BROWNING_NAMED, oxidative, repaired=True
    )
    assert seeds == [stated] * len(seeds), (
        f"under {oxidative} the per-member builder overwrote a scenario input with a parameter "
        "draw — D-24's exclusion is breached"
    )
    assert len(set(drawn)) == len(drawn), "the draw produced no spread; the arm is vacuous"
    assert max(abs(e - nominal) for e in endpoints) / nominal > 0.05, (
        f"under {oxidative} a drawn reference no longer moves browning for a wine that stated "
        "its copper — the mean-centring has lost its physics, not just its artefact"
    )


@pytest.mark.parametrize("scenario", BATTERY + (WINE_BROWNING,), ids=lambda s: s.name)
def test_the_census_itself_is_the_same_under_every_oxidative_wiring(scenario):
    """The census was scored over the DIRECT wiring only (D-234 §9). It is wiring-invariant.

    This is the widest reading of "the alternate wiring is unmeasured": not just whether the
    copper repair survives, but whether swapping the oxidative set puts any *new* name into the
    compile-read-AND-sampled overlap. It does not, in either direction, on any battery member —
    which is what makes D-234's 32-name census a statement about the engine rather than about
    one of its three wirings.

    The reason is structural and worth stating, because it is what a future beat would break:
    the three sets differ in which **Processes** are wired, and a Process is a runtime consumer.
    Compile-time reads come from the seam's seeding and the verbs, which do not branch on the
    oxidative set. A wiring that ever seeds its own state from a parameter would land here.
    """
    baseline = _census(scenario)[1]
    direct_names = set(compile_scenario(scenario).process_set.enabled_snapshot())
    for oxidative in _ALT_WIRINGS:
        with _recording() as seen:
            compiled = compile_scenario(scenario, oxidative=oxidative)
        members = seen & set(
            _resolve_sample_names(
                compiled.process_set, compiled.parameters, None, None, compiled.events
            )
        )
        # An equality is satisfiable by a harness that changed nothing, and `oxidative=` is
        # exactly the kind of argument that could be silently ignored. The wirings must really
        # differ before "the census is invariant across them" says anything at all. Membership
        # rather than `.active`: every oxidative Process is disabled until `begin_aging`, so at
        # compile time the three wirings have identical ACTIVE sets and this check written the
        # obvious way passes on all four scenarios while comparing nothing (it did). The one
        # exemption is named rather than hidden: beer's burst set IS beer's direct set, asserted
        # by `test_beers_burst_set_is_beers_direct_set_so_it_is_not_a_third_number`, so for beer
        # that arm is a genuine identity and the beer census is still compared against a real
        # difference by the cascade arm [[feedback-count-and-print-your-skips]].
        identical_by_design = scenario.medium == "beer" and oxidative == "direct_burst"
        assert (
            identical_by_design or set(compiled.process_set.enabled_snapshot()) != direct_names
        ), (
            f"`oxidative={oxidative!r}` produced the direct wiring's Process set on "
            f"{scenario.name} — this test is comparing the census against itself"
        )
        assert members == baseline, (
            f"the census under {oxidative} differs from the direct wiring's by "
            f"{sorted(members ^ baseline)}. Every name here is half-pinned by construction, so a "
            "new one is a new defect of D-206's class — classify it in `_VERDICTS` and measure it"
        )


def test_the_burst_seed_carries_a_fifty_fold_band_the_ensemble_never_draws():
    """A DEFECT is pinned here, on purpose — the census COMPLEMENT, found by widening to D-237.

    ``burst_antioxidant_initial`` is read at compile on every wiring (it seeds the
    ``burst_antioxidant`` pool under ``oxidative="direct_burst"`` and 0.0 otherwise) and its band
    is **0.0005 - 0.0033 - 0.025 g/L, a 50x span**. No Process declares it in ``reads``:
    ``AntioxidantBurstOxidation`` reads its rate constants, never its seed. So the sampler cannot
    reach it, and under the one wiring where that pool is the whole substrate its uncertainty is
    absent from every reported band.

    **This is not the census predicate and not the drawability surface.** The census is
    compile-read AND sampled — half-pinned, where the surviving half can carry a wrong sign. The
    drawability surface is banded AND drawn but unable to move the run (D-157's oak yield). This
    is the third cell: banded, able to move the run, and never drawn. Nothing guarded it, which is
    why it took a wiring sweep to surface.

    **A RED means it was FIXED** — some Process now declares the seed, or the band was retired.
    Do not revert that beat; delete this guard and say so in the record. It is pinned rather than
    repaired here because wiring a seed into ``reads`` changes what every burst ensemble reports,
    and shipping unmeasured movement to close a gap found in passing is the trap D-233 §6 named.
    """
    compiled = compile_scenario(WINE_BROWNING, oxidative="direct_burst")
    seed = compiled.parameters["burst_antioxidant_initial"]
    sampled = set(
        _resolve_sample_names(
            compiled.process_set, compiled.parameters, None, None, compiled.events
        )
    )

    assert seed.uncertainty.high / seed.uncertainty.low > 10.0, (
        "the band narrowed; re-measure what this gap now costs before trusting the old number"
    )
    assert "burst_antioxidant_initial" not in sampled, (
        "`burst_antioxidant_initial` is now drawn. If a beat wired the seed into a Process's "
        "`reads`, that is the fix — delete this guard and record the band it opened."
    )
    # The positive control the "not in" owes: `sampled` is a real set, not an empty one, and it
    # does contain a compile-read seed — so "never drawn" is this parameter's property and not
    # the harness's [[feedback-a-null-result-needs-a-positive-control]].
    assert "copper_typical" in sampled, "the sampler resolved nothing; the assertion above is void"
    # ...and the pool it seeds really is live under this wiring, so the gap is not academic.
    assert float(compiled.y0[compiled.schema.slice("burst_antioxidant")][0]) == seed.value
