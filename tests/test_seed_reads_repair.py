"""The D-45 fallback seeds are drawn, and the seed and the draw cannot come apart.

Decision D-241. D-240 §3 enumerated eight ``LIVE SEED, NEVER DRAWN`` rows — parameters the
compile seam reads to fill a state slot, that no Process reads at runtime — priced them, and
declined to repair them. Its reason was mechanical: drawing such a name does nothing on its own,
because the sampler re-draws the *parameter map* while the seed is already baked into ``y0``. A
repair therefore needs two coupled changes, and the second one looked like a declaration on some
Process's ``reads`` — a tier-propagation claim (``reads`` has two masters, D-160) that no
measurement supported. D-237 declined the same thing for one name and said why.

**The repair takes the other door.**
:attr:`~fermentation.scenario.compile.CompiledScenario.seed_reads` is a sampling-scope channel
that ``run_ensemble`` unions into the sampled set, and it is
**derived from the ``y0`` rules themselves**. So it carries no tier claim, and the half-repair
D-240 §10 warned about — *"the ensemble is drawing a name that cannot reach ``y0``, which is worse
than the gap"* — is unrepresentable rather than guarded. That derivation is the load-bearing
property of this beat, and it is asserted first, in
:func:`test_a_name_is_drawable_exactly_when_a_rule_reseeds_it`.

**What it is worth.** Paired on identical draws against a *shipped-before* arm — rules 1-3 live,
rule 4 off, so three earlier records' movement is not attributed to this one — the battery wine's
reported band widens **2.83x** for ``dms``, **2.07x** for ``methanethiol``, 1.36x for
``bound_h2s`` and 1.04x for ``E``; under ``direct_burst`` the ``burst_antioxidant`` band widens
**6.97x**. Beer's is unchanged to 1.000 on every slot, which is this module's null control and
not an oversight (:func:`test_beers_band_is_unchanged_and_the_rule_still_fired`).

**Two rows of D-240's eight are deliberately NOT repaired**, for a measured reason rather than a
cautious one: ``biomass_N_yield_log_intercept``/``_slope`` seed no slot at all — they derive the
``biomass_N_fraction`` override, and *that* parameter is itself sampled, over a band which
strictly contains the range the two coefficients imply and is 2.11x wider. Drawing them too would
put two bands on one physical quantity. They stay in :mod:`tests.test_banded_undrawn_census`'s
registry under a ``SUBSUMED`` verdict.

Measurements: ``M:\\claud_projects\\temp\\ferment\\d241-seed-repair\\``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping

import numpy as np
import pytest

from fermentation.core.state import FloatArray
from fermentation.runtime import simulate_ensemble
from fermentation.runtime.ensemble import _resolve_sample_names
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario.compile import _SEED_FALLBACKS, CompiledScenario
from tests.test_compile_sampled_census import BEER, WINE, WINE_BROWNING

#: Members enough to separate a band from its own sampling noise, few enough that the paired arms
#: stay affordable: every fixture here runs TWO ensembles, so the module pays 2x per scenario.
N_MEMBERS = 24

#: The battery, with each scenario under the wiring its own seeds live in. ``WINE_BROWNING`` is
#: scored under ``direct_burst`` for D-240 Arm C's reason: elsewhere D-147 zeroes the burst slot,
#: rule 4 correctly declines, and a guard written at the default would report a working repair as
#: a missing one [[feedback-a-guard-must-be-scored-where-its-subject-lives]].
_WIRED = ((WINE, "direct"), (BEER, "direct"), (WINE_BROWNING, "direct_burst"))


def _slot_index(compiled: CompiledScenario, slot: str) -> int:
    s = compiled.schema.slice(slot)
    return int(s.start) if isinstance(s, slice) else int(s)


def _spread(ens, compiled: CompiledScenario, slot: str) -> float:
    v = ens.members[:, _slot_index(compiled, slot), -1]
    return float(v.max() - v.min())


def _without_rule_4(
    compiled: CompiledScenario,
) -> Callable[[Mapping[str, float]], FloatArray] | None:
    """The **shipped-before** ``y0`` builder: rules 1-3 live, rule 4 removed.

    Deliberately not ``y0_for_member=None``. That would also switch off D-233's pH anchor,
    D-236's copper seed and D-238's peptide capacity, crediting three earlier records' movement
    to this one. A rule belongs to D-241 exactly when it declares a name — the same fact
    ``seed_reads`` is built from — so this filter cannot drift from what it is meant to isolate.
    """
    rules = [(names, fn) for names, fn in compiled._member_seed_rules() if not names]
    if not rules:
        return None
    base = compiled.y0

    def build(values: Mapping[str, float]) -> FloatArray:
        out = base.copy()
        for _, fn in rules:
            fn(out, values)
        return out

    return build


class _Paired:
    """Two ensembles over the SAME draws: shipped-before, and the repair."""

    def __init__(self, scenario: Scenario, oxidative: str = "direct") -> None:
        self.compiled = compile_scenario(scenario, oxidative=oxidative)
        grid = np.linspace(0.0, self.compiled.t_span_h[1], 60)
        # One `only` list for BOTH arms — the post-repair name set. Letting one arm default while
        # the other is narrowed would draw different members (a different name set reorders the
        # sample), and max-minus-min over 24 members is far too noisy to carry that comparison:
        # D-240 §7 got a 1.672x "widening" and a 0.42x "narrowing" that way, the same noise
        # wearing two signs, and had to throw both away
        # [[feedback-pair-the-arms-before-comparing-spreads]].
        self.names = list(
            _resolve_sample_names(
                self.compiled.process_set,
                self.compiled.parameters,
                None,
                None,
                self.compiled.events,
                self.compiled.seed_reads,
            )
        )
        common = {"n_members": N_MEMBERS, "seed": 0, "t_eval": grid, "only": self.names}
        self.before = self.compiled.run_ensemble(
            **common, y0_for_member=_without_rule_4(self.compiled)
        )
        self.after = self.compiled.run_ensemble(**common)

    def ratio(self, slot: str) -> float:
        b = _spread(self.before, self.compiled, slot)
        return _spread(self.after, self.compiled, slot) / b if b else float("nan")

    def worst_member_shift(self, slot: str) -> float:
        i = _slot_index(self.compiled, slot)
        return float(np.abs(self.after.members[:, i, -1] - self.before.members[:, i, -1]).max())


@pytest.fixture(scope="module")
def wine_pair() -> _Paired:
    return _Paired(WINE)


@pytest.fixture(scope="module")
def beer_pair() -> _Paired:
    return _Paired(BEER)


@pytest.fixture(scope="module")
def burst_pair() -> _Paired:
    return _Paired(WINE_BROWNING, oxidative="direct_burst")


# ==========================================================================================
# The invariant the whole repair rests on
# ==========================================================================================


@pytest.mark.parametrize(("scenario", "oxidative"), _WIRED, ids=["wine", "beer", "burst"])
def test_a_name_is_drawable_exactly_when_a_rule_reseeds_it(scenario, oxidative):
    """``seed_reads`` and the ``y0`` rules are one list, so neither half can ship alone.

    This is the guard that makes the rest of the repair safe, and it is an **iff** because the
    two failure modes are different sizes of wrong:

    * a name drawn with no rule to re-seed it is the state D-240 §10 calls *worse than the gap* —
      every member would draw a value that never reaches ``y0``, so the reported band would gain
      a dimension of pure noise;
    * a rule with no name in ``seed_reads`` is the silent one: the slot would be re-seeded from a
      value the sampler never varies, which is byte-for-byte the defect being repaired.

    Scored against the RULES rather than a hand-written list of names, because a hand-written
    list is exactly the thing that goes stale [[feedback-a-doc-rots-where-it-duplicates]].
    """
    compiled = compile_scenario(scenario, oxidative=oxidative)
    declared = {name for names, _ in compiled._member_seed_rules() for name in names}
    assert set(compiled.seed_reads) == declared
    assert len(compiled.seed_reads) == len(set(compiled.seed_reads)), "seed_reads must be a set"
    assert declared, "no rule declared anything, so the equality above is vacuous here"

    # Every declared name really moves `y0` when a member draws it — the mechanical reach the
    # equality cannot promise on its own [[feedback-a-control-needs-mechanical-reach]].
    build = compiled.y0_for_member()
    assert build is not None
    for name in sorted(declared):
        assert name in compiled.parameters, f"{name} is declared but is not a parameter"
        values = dict(compiled.parameters.resolve())
        values[name] = compiled.parameters[name].uncertainty.high
        assert not np.array_equal(build(values), compiled.y0), (
            f"{name} is in seed_reads but no rule moves y0 for it — the declaration half without "
            "the seed half is the state D-240 §10 calls worse than the gap"
        )

    # …and the union really did reach the sampler, which is the only reason `seed_reads` exists.
    sampled = set(
        _resolve_sample_names(
            compiled.process_set,
            compiled.parameters,
            None,
            None,
            compiled.events,
            compiled.seed_reads,
        )
    )
    assert declared <= sampled


def test_the_fallback_table_is_covered_and_every_row_is_reachable():
    """Each :data:`_SEED_FALLBACKS` row must fire somewhere in the battery.

    A row that fires nowhere is a rule describing a tree that is gone. It would keep passing the
    iff above — vacuously, since it declares nothing — while the seed it names went quietly back
    to being undrawn. The skip is counted rather than hidden
    [[feedback-count-and-print-your-skips]].
    """
    fired: set[str] = set()
    for scenario, oxidative in _WIRED:
        fired |= set(compile_scenario(scenario, oxidative=oxidative).seed_reads)
    never = sorted({row.param for row in _SEED_FALLBACKS} - fired)
    assert not never, (
        f"{len(never)} of {len(_SEED_FALLBACKS)} fallback rows never fire on any battery member: "
        f"{never}. Either the seam stopped seeding that slot from the parameter — delete the row "
        "and record the movement — or the battery no longer reaches it, and it must be widened."
    )


# ==========================================================================================
# The scope channel: what it deliberately does NOT change
# ==========================================================================================


def test_an_explicit_only_is_still_exactly_what_was_asked_for():
    """``seed_reads`` is unioned into the DEFAULT branch alone.

    ``only`` is the escape hatch every one-parameter sweep in this repo is built on — D-240's own
    pricing harness included. Widening it silently would put four extra draws into sweeps written
    to isolate one name, and every such measurement would quietly become a joint one
    [[feedback-a-one-parameter-sweep-is-not-the-band]].
    """
    compiled = compile_scenario(WINE)
    assert "dms_potential_initial" in compiled.seed_reads
    names = _resolve_sample_names(
        compiled.process_set,
        compiled.parameters,
        ["mu_max"],
        None,
        compiled.events,
        compiled.seed_reads,
    )
    assert names == ("mu_max",)


def test_a_seed_name_can_still_be_pinned_by_exclude():
    """The union happens BEFORE ``exclude``, so a caller can pin a seed back to its nominal.

    ``exclude`` is how this repo pins anything it needs held still (D-24's escape hatch). A seed
    that could not be pinned would be the only sampled name in the engine with no way off.
    """
    compiled = compile_scenario(WINE)
    names = set(
        _resolve_sample_names(
            compiled.process_set,
            compiled.parameters,
            None,
            ["dms_potential_initial"],
            compiled.events,
            compiled.seed_reads,
        )
    )
    assert "dms_potential_initial" not in names
    assert "bound_h2s_initial" in names, "exclude removed more than it was asked to"


def test_a_direct_simulate_ensemble_caller_is_byte_identical():
    """``seed_reads`` defaults to ``()``, so the low-level API is unchanged across this beat.

    The same compatibility contract D-233 made for ``y0_for_member``: the repair is opted into by
    ``CompiledScenario.run_ensemble`` and by nothing else, so a caller that assembled its own
    Process set and parameters sees no new draws.
    """
    compiled = compile_scenario(WINE)
    grid = np.linspace(0.0, compiled.t_span_h[1], 20)
    plain = simulate_ensemble(
        compiled.process_set,
        compiled.parameters,
        compiled.y0,
        compiled.t_span_h,
        n_members=4,
        seed=0,
        t_eval=grid,
        events=compiled.events,
    )
    assert compiled.seed_reads, "the scenario has no seeds, so this test proves nothing here"
    assert not set(compiled.seed_reads) & set(plain.sampled_names)


def test_a_scenario_that_states_the_level_is_not_re_seeded():
    """A stated initial level is a scenario INPUT and stays off the sampled axis (D-24).

    The same branch D-236 drew for ``copper_gpl``: a wine that names ``dms_potential_ugl`` is
    stating *this must's* DMS precursor, a different quantity from the sourced typical level the
    fallback stands in for. Re-seeding it from a draw would overwrite the recipe.
    """
    stated = WINE.model_copy(update={"initial": {**WINE.initial, "dms_potential_ugl": 40.0}})
    compiled = compile_scenario(stated)
    assert "dms_potential_initial" not in compiled.seed_reads
    # …and its neighbours in the same table are untouched, so the branch is per row rather than a
    # switch that turns the whole rule off.
    assert "bound_h2s_initial" in compiled.seed_reads

    build = compiled.y0_for_member()
    assert build is not None
    slot = compiled.schema.slice("dms_potential")
    values = dict(compiled.parameters.resolve())
    values["dms_potential_initial"] = compiled.parameters["dms_potential_initial"].uncertainty.high
    assert float(build(values)[slot][0]) == float(compiled.y0[slot][0])


def test_the_no_rule_branch_is_still_reachable():
    """``y0_for_member() is None`` still happens, and is still the byte-identical path.

    D-233 shipped that branch as its compatibility contract, and D-236's test reached it through
    a medium that lacked a ``copper`` slot. D-241 closed that door — an un-anchored beer now
    carries a wort-oxygen rule — so the branch is reached here by *stating* the seed instead of
    by relying on a medium to lack it. A branch nothing reaches is a branch nothing tests.
    """
    bare = compile_scenario(
        Scenario(
            name="d241-no-derived-seed",
            medium="beer",
            initial={
                "glucose_gpl": 15.0,
                "maltose_gpl": 60.0,
                "maltotriose_gpl": 10.0,
                "yan_mgl": 200.0,
                "pitch_gpl": 0.5,
                "o2_mgl": 8.0,  # states its own wort oxygen ⇒ the one beer rule declines
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
            duration_days=1.0,
        )
    )
    assert bare.seed_reads == ()
    assert bare.y0_for_member() is None, (
        "an un-anchored beer that states its own oxygen has no parameter-derived seed left; if a "
        "later beat gives it one, move this scenario rather than deleting the branch"
    )


# ==========================================================================================
# The prices — paired arms, and the null control that makes them readable
# ==========================================================================================


@pytest.mark.parametrize(
    ("slot", "floor"), [("dms", 2.0), ("methanethiol", 1.5), ("bound_h2s", 1.2)]
)
def test_the_repair_widens_the_bands_it_was_priced_to_widen(wine_pair, slot, floor):
    """The reported band really does grow, by about what D-240 predicted it would.

    D-240 §7 measured the *hidden* spread — what the sampler would add if it could reach these
    names — at 1.88 of the published band for ``methanethiol`` and 1.09 for ``dms``. Those are
    what a band GAINS, so a ratio above 2 for ``dms`` and above 1.5 for ``methanethiol`` is the
    same finding arriving from the other side. Floors rather than windows: the exact ratio is a
    property of this scenario and this member count, and pinning three digits would turn an
    honest change of either into a red.
    """
    assert wine_pair.ratio(slot) > floor, (
        f"the reported `{slot}` band no longer widens when the seeds are drawn — either rule 4 "
        "stopped firing, or something began absorbing the seed the way beer's t=0 anchor absorbs "
        "its wort acids (D-240 §5). Both are findings, not tolerances to lower."
    )


def test_the_widening_is_not_uniform_and_ethanol_is_the_small_one(wine_pair):
    """Band WIDTH does not order this class, which is the reading D-240 §4 asked for.

    ``must_fermentable_fraction`` is the narrowest band in the whole set — 1.06x — and it moves
    more of the engine than any other row, because it scales the sugar every downstream number is
    a fraction of. Yet it moves the *reported* ethanol band barely at all, because that band is
    already wide with kinetics. Both halves are asserted, since quoting either alone misdescribes
    the repair [[feedback-a-summary-statistic-is-not-the-curve]].
    """
    assert wine_pair.ratio("E") < 1.2, "the ethanol BAND is not where this shows up"
    assert wine_pair.worst_member_shift("E") > 1.0, (
        "…but an individual member's ethanol must move by grams per litre — measured 3.01 g/L "
        "worst case. A small band ratio with no member shift would mean rule 4 stopped firing."
    )
    assert wine_pair.ratio("dms") > wine_pair.ratio("E"), (
        "the ordering is the finding: a 1.06x band outranks a 1.53x one on the engine and not on "
        "the reported spread"
    )


def test_the_burst_wiring_is_where_the_fifty_fold_band_lands(burst_pair):
    """The widest band in D-240's set, and it exists only under the wiring that consumes it."""
    assert "burst_antioxidant_initial" in burst_pair.compiled.seed_reads
    assert burst_pair.ratio("burst_antioxidant") > 3.0
    assert burst_pair.ratio("A420") > 1.15, "the browning the burst pool protects moves with it"


def test_beers_band_is_unchanged_and_the_rule_still_fired(beer_pair):
    """The NULL CONTROL — and the second assertion is the half that makes it one.

    ``o2_wort_aeration_beer`` is repaired by the identical mechanism and is worth nothing: every
    ``o2`` consumer is aging-gated and ``begin_aging`` *disables* the Process that removes it
    (D-213, re-measured on the trajectory at D-240 §6). So beer's reported band must not move,
    which is what shows the repair does not manufacture spread out of the harness.

    **A zero on its own would also be produced by a rule that never ran** — the shape D-240's
    Arm B caught by making its pricer assert its own reach. Here the second assertion IS the
    reach: the member's ``y0`` really does differ, and the band still does not move.
    """
    assert "o2_wort_aeration_beer" in beer_pair.compiled.seed_reads
    for slot in ("X", "E", "o2", "acetic"):
        assert beer_pair.ratio(slot) == pytest.approx(1.0, abs=1e-6), (
            f"beer's reported `{slot}` band moved. Something wired O2 to growth, or stopped "
            "disabling WortOxygenUptake at begin_aging — D-213's scope decision, not a tolerance."
        )
    moved = max(beer_pair.worst_member_shift(s) for s in ("X", "E", "o2"))
    assert moved > 0.0, (
        "the beer seed rule reached nothing at all, so the four unchanged bands above are vacuous"
    )
    assert moved < 1e-6, "…and it must stay inert; a real movement here would contradict D-213"


def test_the_nominal_run_is_untouched_by_the_repair(wine_pair):
    """The deterministic baseline must not move — this repair is about members, not the nominal.

    D-24's byte-for-byte claim. Rule 4 re-seeds from a member's draw, and at the nominal draw
    that IS the compiled value, so both arms' nominals must be bit-identical to each other and
    the compiled ``y0`` unchanged from a fresh compile.
    """
    assert np.array_equal(wine_pair.before.nominal, wine_pair.after.nominal)
    assert np.array_equal(wine_pair.compiled.y0, compile_scenario(WINE).y0)
