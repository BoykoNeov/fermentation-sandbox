"""The banded-AND-compile-read-AND-**undrawn** census — the other half of D-234's set.

Decision D-240. Opened by the owner as D-237 §6's parked item — *"it implies a census of its
own — every banded parameter read at compile that no Process declares — which is not run
here"*.

**The predicate, and how it relates to its two neighbours.** A name is a member here iff it is
(1) read during ``compile_scenario`` and (2) drawn by ``_resolve_sample_names`` in **no** member
of the battery. That is exactly ``seen - members`` in :mod:`tests.test_compile_sampled_census`'s
terms, so this set and D-234's are the two halves of one partition of the compile reads and are
disjoint by construction — ``test_the_census_and_the_drawability_surface_are_disjoint`` over
there already names ``otr_screwcap`` and ``oak_yield_vanillin_medium`` as living here.

**It is NOT disjoint from D-153…D-159's drawability surface, and that is the point to read
first.** Eleven of the twenty-eight banded rows below — the six oak yields, ``otr_screwcap``,
``bottling_burst_screwcap`` and the three copper dose constants — are D-159's *structural* class
by another route. D-159 enumerated 61 such names and deliberately declined to pin the list,
pinning the **mechanism** on one exemplar per class instead. So enumeration is not the new thing
here. The new thing is the question D-159 never asked: **of the names the sampler can never
reach, which ones would have moved the run if it could** — i.e. whose band is uncertainty that
every reported spread silently omits. D-237 §4 named that cell (*"banded, able to move the run,
and never drawn"*) on one instance; this is the sweep.

**Why the battery is D-234's and is deliberately NOT widened.** Six closures, three toasts and
the bourbon soak would add ~40 names to the census and not one finding: ``oak_yield_vanillin_heavy``
carries ``_medium``'s verdict for ``_medium``'s reason. Variant enumeration inflates a count and
answers nothing [[feedback-count-and-print-your-skips]] — the skip is stated here rather than
hidden, and a battery that ever reaches another variant goes red on
:func:`test_every_banded_undrawn_name_is_classified` until its row is written.

**Nothing was repaired here, and D-241 repaired six of them — read the two together.** The
reason this record declined was mechanical: drawing one of these names would not, on its own,
do anything, because the sampler re-draws the parameter map while the seed is already baked into
``y0`` at compile. So a repair needs a ``y0_for_member`` rule (D-236's machinery) AND a way for
the name to enter the sampled set — and the second half looked like a declaration on a Process,
which is a tier-propagation claim (``reads`` has two masters, D-160) that no measurement here
supports. D-241 took the other route: ``CompiledScenario.seed_reads``, a sampling-scope channel
DERIVED from the ``y0`` rules, which carries no tier claim and cannot deliver either half of the
repair without the other. Six rows left this census as a result; the two that stayed are
:data:`_SUBSUMED_FIT`, and they stayed for a measured reason rather than the mechanical one.

**A guard in this module forbade nothing against that route, and that is the part worth reading.**
:func:`test_the_priced_names_are_still_undrawn` scored ``_resolve_sample_names`` without the new
argument, so it kept reporting the repaired names as undrawn and stayed GREEN through the very
event its error message describes. It was written against the one repair anyone had imagined —
and a repair that arrived by another door walked straight past it. That is D-240's own Arm C
lesson recurring one layer up [[feedback-a-guard-must-be-scored-where-its-subject-lives]]: the
instrument, not just the assertion, has to be scored where the subject lives.

Measurements: ``M:\\claud_projects\\temp\\ferment\\d240-banded-undeclared\\``.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.acidbase import ph_of_state
from fermentation.parameters.schema import Parameter
from fermentation.parameters.store import ParameterSet
from fermentation.runtime import simulate_ensemble, simulate_scheduled
from fermentation.runtime.ensemble import _resolve_sample_names
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario import compile as compile_mod
from fermentation.scenario.compile import CompiledScenario
from tests.test_compile_sampled_census import BATTERY, BEER, WINE, _census

# ==========================================================================================
# The registry: name -> verdict. Six classes, and the class is the finding.
# ==========================================================================================

#: The D-45 shape — a scenario that does not state a quantity falls back to a **sourced** level
#: rather than to 0, because "absent" does not mean "none". Every one of these seeds a live pool
#: and none is drawn, so its band is uncertainty no ensemble reports. These are the rows the
#: owner's "fix or pin" call is actually about.
#: **D-241 REPAIRED these six and they are no longer members of this census.** They are drawn
#: now — `CompiledScenario.seed_reads` unions them into the sampled set, and a `y0_for_member`
#: rule re-seeds the slot each one fills — so they moved to `test_compile_sampled_census.py`'s
#: CENSUS, which is where a repaired name lives. Kept here as a NAMED tuple rather than deleted,
#: because `test_the_six_repaired_seeds_really_are_drawn_now` is the mirror of the pin they used
#: to sit under: same teeth, opposite sign. A RED there means the repair was reverted.
_REPAIRED_AT_D241 = (
    "burst_antioxidant_initial",
    "dms_potential_initial",
    "bound_h2s_initial",
    "bound_methanethiol_initial",
    "must_fermentable_fraction",
    "o2_wort_aeration_beer",
)

#: The two rows D-240 filed under LIVE SEED that D-241 did NOT repair, because they seed no slot
#: at all: they derive the `biomass_N_fraction` override, and that parameter is itself sampled.
#: Measured at D-241 §2 rather than argued — at the battery wine's YAN the two coefficients'
#: own bands imply f_N in [0.051432, 0.108338], which sits STRICTLY INSIDE the override's own
#: [0.03, 0.15] and is 0.474 as wide. Drawing them as well would put two bands on one physical
#: quantity and report the wider one as narrower. Still undrawn, still pinned, now with a reason
#: that is a measurement instead of "different mechanism".
_SUBSUMED_FIT = (
    "biomass_N_yield_log_intercept",
    "biomass_N_yield_log_slope",
)

#: Beer's eight wort organic-acid levels. They seed the acid slots — but the t=0 cation anchor is
#: back-solved THROUGH those concentrations, so the anchor absorbs them at t=0 exactly (the D-238
#: §4 / D-239 §5 signature) and only the drift survives. Measured, not assumed: see
#: :func:`test_the_wort_acid_seeds_are_absorbed_by_the_t0_anchor`.
_WORT_ACIDS = (
    "acetic_typical_wort",
    "citric_typical_wort",
    "formic_typical_wort",
    "lactic_typical_wort",
    "malic_typical_wort",
    "oxalic_typical_wort",
    "pyruvic_typical_wort",
    "succinic_typical_wort",
)

#: The compile-consumed families D-157/D-159 already swept. Listed by the variant the battery
#: reaches, on purpose (see the module docstring).
_D159_STRUCTURAL = (
    "oak_yield_ellagitannin_medium",
    "oak_yield_eugenol_medium",
    "oak_yield_furaneol_medium",
    "oak_yield_guaiacol_medium",
    "oak_yield_vanillin_medium",
    "oak_yield_whiskey_lactone_medium",
    "otr_screwcap",
    "bottling_burst_screwcap",
)

#: Constants a DOSE VERB reads when it compiles its event — they scale a jump, not a seed.
_DOSE_VERB = (
    "copper_h2s_binding",
    "copper_mercaptan_binding",
    "copper_fining_residual_fraction",
)

#: Zero-width bands. A name here hides nothing *by construction*, which is the one verdict that
#: needs no price — and the one that would rot silently if a later beat banded it, so it is
#: asserted rather than trusted (:func:`test_the_zero_width_rows_really_are_zero_width`).
_ZERO_WIDTH = (
    "dap_nitrogen_charge",
    "dap_nitrogen_fraction",
    "dap_phosphate_fraction",
    "peptide_buffer_capacity_beer",
    "sucrose_inversion_mass_ratio",
)

VERDICTS: Mapping[str, str] = {
    # The one name that is in BOTH registries, and it is not a contradiction: membership is
    # scenario-conditional. `phenolic_browning` is the only Process that reads `copper_typical`
    # and it is disabled until `begin_aging`, so in a wine that never ages the seam still seeds
    # the `copper` slot from it and nothing draws it. D-236's repair re-seeds that slot from the
    # member's own draw — which is exactly a no-op where the member never draws it.
    "copper_typical": (
        "SCENARIO-CONDITIONAL — a D-234 census member in an AGED wine (REPAIRED there by "
        "y0_for_member's rule 2) and a member of THIS census in a wine that never ages, where no "
        "Process declares it. Nothing to repair: the slot it seeds is read by the Process whose "
        "absence is what put the name here"
    ),
    # The second wiring-conditional row, and the exact mirror of `copper_typical` above: that
    # name is a member HERE in the scenario where its reader is absent, and this one is a member
    # here in the WIRING where its reader is absent. Both prove membership is a property of the
    # configuration and not of the name.
    "burst_antioxidant_initial": (
        "WIRING-CONDITIONAL, and INERT in the wiring that puts it here — REPAIRED at D-241 under "
        '`oxidative="direct_burst"`, where it is drawn and re-seeded per member (there its band '
        "widens the reported `burst_antioxidant` spread 6.97x). Under the DEFAULT `direct` set "
        "and under `cascade` it stays a member of this census, and there is nothing to repair: "
        "D-147 zeroes the slot wherever AntioxidantBurstOxidation is not wired, so the seam reads "
        "the value and discards it. Measured, not argued — across the whole 50x band the ENTIRE "
        "y0 is bit-identical under both those wirings. The seed rule's equality guard declines "
        "for exactly that reason, so D-147's condition is honoured without the table knowing it "
        "exists"
    ),
    **dict.fromkeys(
        _SUBSUMED_FIT,
        "SUBSUMED — a live compile read that seeds NO slot: the two Coleman coefficients derive "
        "the `biomass_N_fraction` override, and that parameter IS sampled, over a band which "
        "strictly contains the range these two imply and is 2.11x wider (D-241 §2). So the "
        "uncertainty is not omitted from the reported spread — it is already in it, under "
        "another name and with margin. Drawing these too would double-count one quantity. This "
        "is the one verdict in the registry that says the gap is APPARENT rather than real",
    ),
    **dict.fromkeys(
        _WORT_ACIDS,
        "LIVE SEED, ABSORBED AT t=0 — it seeds a beer acid slot, but the cation anchor is "
        "back-solved through that concentration, so t=0 pH is untouched EXACTLY (D-238 §4's "
        "signature) and only the ferment drift carries the band. Its own class, not eight more "
        "members of the row above",
    ),
    **dict.fromkeys(
        _D159_STRUCTURAL,
        "COMPILE-CONSUMED — D-157/D-159's structural class, reached here by the other predicate. "
        "D-159 pinned the mechanism on one exemplar per class and declined to freeze the 61 "
        "names; nothing about that verdict changes here",
    ),
    **dict.fromkeys(
        ("copper_h2s_binding", "copper_mercaptan_binding"),
        "DOSE-VERB CONSTANT, and worth EXACTLY ZERO for a structural reason worth reading: "
        "_verb_add_copper removes min(pool_present, copper_gpl * binding), and at any dose a "
        "winemaker would use the copper is in ~1000x excess of the sulfide present — so the "
        "constant is the LOSING argument of that min and its band cannot reach the run at all. "
        "That is D-159's supply-limited warning inverted: not 'no substrate', but 'no shortage'",
    ),
    "copper_fining_residual_fraction": (
        "DOSE-VERB CONSTANT — it multiplies the copper credited to the slot, so unlike the two "
        "binding constants it does move: 2.5e-05 g/L across its band, 0.068 of the reported "
        "`copper` spread. Never drawn for the usual reason: no Process declares it"
    ),
    **dict.fromkeys(
        _ZERO_WIDTH,
        "ZERO-WIDTH — low == high, so there is no uncertainty for the sampler to be missing. "
        "peptide_buffer_capacity_beer is additionally re-derived per member by y0_for_member's "
        "rule 3 (D-238), so it is pinned twice over",
    ),
}


# ==========================================================================================
# The instrument — D-234's recorder, used through its own module so there is one of it
# ==========================================================================================


def _undrawn(scenario: Scenario) -> set[str]:
    """Names ``scenario`` reads at compile and that its own run would never draw."""
    seen, members = _census(scenario)
    return seen - members


def _banded_undrawn(scenario: Scenario) -> set[str]:
    compiled = compile_scenario(scenario)
    out = set()
    for name in _undrawn(scenario):
        if name not in compiled.parameters:
            continue
        u = compiled.parameters[name].uncertainty
        if u.low < u.high:
            out.add(name)
    return out


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


def _compiled_with(
    scenario: Scenario,
    name: str | None = None,
    value: float | None = None,
    *,
    oxidative: str = "direct",
):
    """Compile ``scenario`` with ``name`` moved BEFORE the seam reads it.

    A recompile is the only arm that reaches a seed: patching the resolved map moves the runtime
    role and leaves ``y0`` exactly where it was, which is the whole defect being priced.
    """
    if name is None:
        return compile_scenario(scenario, oxidative=oxidative)
    assert value is not None, "naming a parameter without a value moves nothing"
    original = compile_mod._load_parameters

    def patched(sc, parameter_paths, data_dir):
        return _with_value(original(sc, parameter_paths, data_dir), name, float(value))

    compile_mod._load_parameters = patched
    try:
        compiled = compile_scenario(scenario, oxidative=oxidative)
    finally:
        compile_mod._load_parameters = original
    assert compiled.parameters[name].value == float(value), "the recompile did not move it"
    return compiled


def _endpoint(compiled, slot: str) -> float:
    out = simulate_scheduled(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
        events=compiled.events,
    )
    return float(np.atleast_1d(out.y[compiled.schema.slice(slot), -1]).sum())


def _reported_spread(compiled, slot: str, n_members: int) -> float:
    """The spread the ensemble ALREADY publishes for ``slot`` — the denominator of every price.

    Without it a hidden span is a bare number: 0.02 g/L is either negligible or the whole
    reported band, and only the ratio says which.
    """
    ens = simulate_ensemble(
        compiled.process_set,
        compiled.parameters,
        compiled.y0,
        compiled.t_span_h,
        n_members=n_members,
        seed=0,
        events=compiled.events,
        y0_for_member=compiled.y0_for_member(),
    )
    vals = np.asarray(ens.members[:, compiled.schema.slice(slot), -1], float)
    vals = vals.reshape(len(ens.members), -1).sum(axis=1)
    return float(vals.max() - vals.min())


# ==========================================================================================
# Enumeration
# ==========================================================================================


@pytest.mark.parametrize("scenario", BATTERY, ids=lambda s: s.name)
def test_every_banded_undrawn_name_is_classified(scenario):
    """A banded parameter the seam reads and the sampler cannot reach must be a KNOWN one.

    This is the guard the beat exists for. ``burst_antioxidant_initial`` was found in passing
    while D-237 measured something else; before that, D-140 had found the same shape in the same
    parameter and pinned it. A third accident is not a method.
    """
    unclassified = sorted(_banded_undrawn(scenario) - set(VERDICTS))
    assert not unclassified, (
        f"{scenario.name}: {len(unclassified)} banded parameter(s) are read at compile and drawn "
        f"by nothing: {unclassified}. Each carries uncertainty that every reported spread omits. "
        "Price it against the spread the same ensemble publishes for the slot it reaches, then "
        "add the verdict here."
    )


def test_no_classified_name_has_gone_stale():
    """The registry must not accumulate names that stopped being members.

    Scored over the UNION of the battery, because membership is scenario-conditional: the beer
    rows are unreachable from a wine and vice versa.
    """
    banded: set[str] = set()
    undrawn: set[str] = set()
    for scenario in BATTERY:
        banded |= _banded_undrawn(scenario)
        undrawn |= _undrawn(scenario)
    # The zero-width five are members of the UNDRAWN set and not of the banded one, so they are
    # scored against the wider union rather than excluded from the check. Subtracting them
    # outright — the first version of this test — made their rows unfalsifiable: the seam could
    # stop reading one entirely and nothing here or in the zero-width test would notice, which is
    # the same shape D-240's Arm C found in the undrawn pin
    # [[feedback-a-guard-must-be-scored-where-its-subject-lives]].
    stale = sorted((set(VERDICTS) - set(_ZERO_WIDTH) - banded) | (set(_ZERO_WIDTH) - undrawn))
    assert not stale, (
        f"{len(stale)} classified name(s) are no longer banded, no longer compile-read, or have "
        f"become drawable: {stale}. A name that became DRAWABLE is a repair and belongs in "
        "tests/test_compile_sampled_census.py's CENSUS instead — move the row, do not delete it."
    )


def test_the_two_censuses_partition_the_compile_reads():
    """D-234's set and this one are the two halves of one partition, and neither is the whole.

    Stated as a guard because the two are one word apart in prose — *compile-read and sampled*
    against *compile-read and not* — and a reader who collapses them would conclude the compile
    seam had already been audited end to end.

    **The partition is per SCENARIO, and one name proves the distinction matters.** The registries
    are battery-wide unions, so a name can sit in both — ``copper_typical`` is a D-234 member in
    an aged wine and a member here in a wine that never ages, because the only Process that reads
    it is disabled until ``begin_aging``. Asserted in the positive direction rather than
    forbidden: a cross-registry name is the *expected* shape for anything an aging gate controls,
    and a test that outlawed it would go red on the first honest one.
    """
    from tests.test_compile_sampled_census import CENSUS

    for scenario in BATTERY:
        seen, members = _census(scenario)
        undrawn = seen - members
        assert undrawn | members == seen, "the two halves must cover every compile read"
        assert not (undrawn & members), "and must not overlap"
    assert "copper_typical" in CENSUS and "copper_typical" in VERDICTS
    assert "copper_typical" in _banded_undrawn(BATTERY[2]), "the un-aged wine's row"
    assert "copper_typical" in _census(BATTERY[0])[1], "the aged wine's row"


def test_the_zero_width_rows_really_are_zero_width():
    """The one verdict that needs no price is the one that rots silently if a beat bands it."""
    wine = compile_scenario(WINE)
    beer = compile_scenario(BEER)
    reached = set()
    for scenario in BATTERY:
        reached |= _undrawn(scenario)
    for name in _ZERO_WIDTH:
        assert name in reached, (
            f"{name} is no longer read at compile, or has become drawable — either way its "
            "zero-width row is describing a tree that no longer exists."
        )
        params = wine.parameters if name in wine.parameters else beer.parameters
        u = params[name].uncertainty
        assert u.low == u.high, (
            f"{name} now carries a band ({u.low} to {u.high}) and is drawn by nothing, so its "
            "uncertainty is absent from every reported spread. Price it and move its row."
        )


#: seed name -> (medium, the y0 slot it seeds, the oxidative wiring). The two biomass-yield
#: coefficients are deliberately absent: they seed no slot, they derive a *parameter*
#: (``_apply_nitrogen_dependent_yield``), and they get their own check below.
_SEEDED_SLOT = {
    "burst_antioxidant_initial": ("wine", "burst_antioxidant", "direct_burst"),
    "dms_potential_initial": ("wine", "dms_potential", "direct"),
    "bound_h2s_initial": ("wine", "bound_h2s", "direct"),
    "bound_methanethiol_initial": ("wine", "bound_methanethiol", "direct"),
    "must_fermentable_fraction": ("wine", "S", "direct"),
    "o2_wort_aeration_beer": ("beer", "o2", "direct"),
    "copper_typical": ("wine", "copper", "direct"),
    # The slot names are the acids' own, EXCEPT citric, whose slot is `citrate` — the one place
    # the parameter family and the state schema disagree, and a silent KeyError otherwise.
    **{
        f"{acid}_typical_wort": ("beer", slot, "direct")
        for acid, slot in (
            ("acetic", "acetic"),
            ("citric", "citrate"),
            ("formic", "formic"),
            ("lactic", "lactic"),
            ("malic", "malic"),
            ("oxalic", "oxalic"),
            ("pyruvic", "pyruvic"),
            ("succinic", "succinic"),
        )
    },
}


@pytest.mark.parametrize("name", sorted(_SEEDED_SLOT))
def test_a_band_edge_moves_the_slot_the_name_seeds(name):
    """The CONSEQUENTIAL half of D-159's two-assertion idiom, paid at compile.

    "Never drawn" is only a finding for a name that would have mattered. D-159 measured the
    other failure mode directly: five genuinely Process-read names froze at exactly zero movement
    in its harness for supply-limited reasons, so "the trajectory did not move" cannot tell
    *unreachable* from *reachable and zero here*. Here the reach is checked one step earlier and
    without an integration at all — move the parameter to each band edge, recompile, and read the
    slot out of ``y0``. A row that fails this is inert and belongs in another class.
    """
    medium, slot, oxidative = _SEEDED_SLOT[name]
    scenario = WINE if medium == "wine" else BEER
    base = _compiled_with(scenario, oxidative=oxidative)
    u = base.parameters[name].uncertainty

    def seeded(value: float) -> float:
        c = _compiled_with(scenario, name, value, oxidative=oxidative)
        return float(np.atleast_1d(c.y0[c.schema.slice(slot)]).sum())

    lo, hi = seeded(u.low), seeded(u.high)
    assert lo != hi, (
        f"{name} moved across its whole band ({u.low} to {u.high}) and the `{slot}` seed did not "
        "move at all — it is not the live row this registry says it is."
    )


def test_the_biomass_yield_coefficients_move_a_derived_parameter_not_a_slot():
    """The two coefficients that are consumed as a *fit*, not as a seed — the other live shape.

    ``_apply_nitrogen_dependent_yield`` reads both at compile and OVERRIDES
    ``biomass_N_fraction`` in the parameter map with a value fitted to the scenario's own YAN.
    So the band is hidden in a place ``y0`` cannot show it, and the check that catches it for the
    seeds above is structurally blind here.
    """
    derived = "biomass_N_fraction"
    base = compile_scenario(WINE)
    assert derived in base.param_values, "the nitrogen-dependent yield writes another name now"
    for name in ("biomass_N_yield_log_intercept", "biomass_N_yield_log_slope"):
        u = base.parameters[name].uncertainty
        values = {
            edge: _compiled_with(WINE, name, v).param_values[derived]
            for edge, v in (("low", u.low), ("high", u.high))
        }
        assert values["low"] != values["high"], (
            f"{name} no longer reaches {derived}; re-classify it before keeping the row."
        )
        assert min(values.values()) < base.param_values[derived] < max(values.values()), (
            "the nominal must sit strictly inside the two edges, or the band is one-sided"
        )


# ==========================================================================================
# The prices — what each undrawn band is worth against the spread that IS published
# ==========================================================================================


def _ph_course(compiled: CompiledScenario, hours: tuple[float, ...]) -> dict[float, float]:
    grid = np.linspace(compiled.t_span_h[0], compiled.t_span_h[1], 601)
    out = simulate_scheduled(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
        events=compiled.events,
        t_eval=grid,
    )
    return {
        h: ph_of_state(
            out.y[:, int(np.argmin(np.abs(out.t - h)))], compiled.schema, compiled.param_values
        )
        for h in hours
    }


#: `acetic` carries the largest day-1 excursion of the eight and `formic` one of the smallest,
#: which is the pair that shows the signature is about the ANCHOR rather than about size.
@pytest.mark.parametrize("name", ("acetic_typical_wort", "formic_typical_wort"))
def test_the_wort_acid_seeds_are_absorbed_by_the_t0_anchor(name):
    """t=0 is untouched EXACTLY; the band survives only as ferment drift (D-240 §arm A).

    The seam back-solves beer's ``cation_charge`` so the model reproduces ``initial_ph`` — and it
    solves it *through* these very concentrations. So moving one to a band edge moves the anchor
    by the same amount with the opposite sign, and t=0 pH does not move at all. That is the same
    signature D-238 §4 and D-239 §5 measured for two other re-partitions of a species present when
    the anchor is taken, and it is why these eight are their own verdict rather than eight more
    live seeds: what their band actually buys is a ~1e-3 pH wobble on day 1, not a different wort.

    Measured (12 d beer, both edges): acetic +-1.11e-3 at day 1 falling to +-9.8e-5 by day 7;
    formic +-9.9e-4 at day 1 falling to +-1.4e-6. The t=0 residual is 1e-13 or smaller for all
    sixteen arms — solver-free arithmetic, since it is read off ``y0``'s own pH.
    """
    base = compile_scenario(BEER)
    u = base.parameters[name].uncertainty
    hours = (0.0, 24.0, 168.0)
    nominal = _ph_course(base, hours)
    for edge in (u.low, u.high):
        moved = _ph_course(_compiled_with(BEER, name, edge), hours)
        assert abs(moved[0.0] - nominal[0.0]) < 1e-12, (
            f"{name} at {edge} moved t=0 pH by {moved[0.0] - nominal[0.0]:.3e}. The anchor is "
            "supposed to absorb it exactly; a residual here means the back-solve stopped reading "
            "this concentration and every beer now starts at a pH the scenario never asked for."
        )
        assert abs(moved[24.0] - nominal[24.0]) > 5e-5, (
            f"{name} at {edge} moved day-1 pH by only {moved[24.0] - nominal[24.0]:.3e} — the "
            "drift the anchor CANNOT absorb has gone, so the row is inert and not merely absorbed."
        )


def test_the_wort_aeration_seed_is_inert_in_beers_default_set():
    """The one banded, never-drawn, LIVE seed whose price is measurably zero (D-213).

    ``o2_wort_aeration_beer`` seeds beer's dissolved oxygen across a 1.45x band that no ensemble
    can draw — and it does not matter, because every Process that consumes ``o2`` is aging-gated
    and ``begin_aging`` *disables* the Process that removes it. D-213 argued the inertness on the
    derivative; this measures it on the trajectory, from the other side, as part of the sweep.
    Keeping the row is the point: "banded and undrawn" is not by itself a defect, and this is the
    row that proves the sweep can tell the difference.
    """
    base = compile_scenario(BEER)
    u = base.parameters["o2_wort_aeration_beer"].uncertainty
    for slot in ("X", "E"):
        lo = _endpoint(_compiled_with(BEER, "o2_wort_aeration_beer", u.low), slot)
        hi = _endpoint(_compiled_with(BEER, "o2_wort_aeration_beer", u.high), slot)
        assert abs(hi - lo) < 1e-6 * max(abs(lo), 1.0), (
            f"the wort-aeration seed now reaches `{slot}` ({lo!r} vs {hi!r}). Something wired O2 "
            "to growth or stopped disabling WortOxygenUptake at begin_aging — D-213's scope "
            "decision, not a tolerance to widen."
        )


#: The measurement wine: no interventions, no aging, 14 days — the shortest run that still
#: separates the two names below, because the ensemble in this test is paid per member.
_PRICE_WINE = Scenario(
    name="d240-price-wine",
    medium="wine",
    initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25},
    temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
    duration_days=14.0,
)


def test_an_undrawn_band_is_worth_as_much_as_a_drawn_one():
    """The headline, and the positive control that stops it being decoration.

    A hidden span is meaningless on its own — 6.43 g/L of ethanol is either noise or everything.
    The denominator is the spread the SAME ensemble already publishes for the same slot with the
    name undrawn, and the control is ``mu_max``: a name the sampler does draw, priced through the
    identical harness. Measured on a 14 d wine, 12 members, seed 0:

    * ``must_fermentable_fraction`` (undrawn) moves final ethanol 6.43 g/L across its band,
      against a reported spread of 11.19 — **0.574**.
    * ``mu_max`` (drawn) moves final biomass 0.2014 g/L against 0.4565 — **0.441**.

    So the undrawn name is not a rounding error next to the drawn one; it is the same order. That
    is the whole finding, and without the control it would read as an arbitrary fraction
    [[feedback-a-one-parameter-sweep-is-not-the-band]] — note both numbers are edge-to-edge
    excursions against a *sampled* spread, so both are upper bounds in the same way.
    """
    base = compile_scenario(_PRICE_WINE)
    prices = {}
    for name, slot in (("must_fermentable_fraction", "E"), ("mu_max", "X")):
        u = base.parameters[name].uncertainty
        lo = _endpoint(_compiled_with(_PRICE_WINE, name, u.low), slot)
        hi = _endpoint(_compiled_with(_PRICE_WINE, name, u.high), slot)
        prices[name] = abs(hi - lo) / _reported_spread(base, slot, n_members=12)

    assert prices["must_fermentable_fraction"] == pytest.approx(0.574, abs=0.02)
    assert prices["mu_max"] == pytest.approx(0.441, abs=0.02)
    assert prices["must_fermentable_fraction"] > prices["mu_max"], (
        "the undrawn band used to be worth MORE of the reported ethanol spread than mu_max is of "
        "the biomass one. If that ordering flipped, re-measure before quoting either."
    )


#: A wine that doses copper mid-ferment, where H2S is at its highest — the scenario most likely
#: to make the binding constant bite, chosen on purpose so the null is not a soft one.
_COPPER_DOSE_DAY = 8.0
_COPPER_DOSE_GPL = 5.0e-4
_COPPER_WINE = Scenario(
    name="d240-copper-wine",
    medium="wine",
    initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25, "so2_total_mgl": 30.0},
    temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
    interventions=[
        Intervention(day=_COPPER_DOSE_DAY, action="add_copper", params={"copper_mgl": 0.5}),
    ],
    duration_days=12.0,
)


def test_the_copper_binding_constants_lose_their_own_min_at_any_real_dose():
    """Why two banded, never-drawn dose constants are worth exactly zero — and it is structural.

    ``_verb_add_copper`` removes ``min(h2s_present, copper_gpl * copper_h2s_binding)``. A 0.5
    mg/L dose buys ~2.7e-04 g/L of binding capacity against a wine carrying ~1e-07 g/L of H2S,
    so the ``min`` takes the pool every time and the constant never appears in the arithmetic.
    Its 1.99x band is therefore unreachable twice over: the sampler cannot draw it, and drawing
    it would not matter.

    This is the case D-159's docstring warns about, with the sign reversed — there a parameter
    froze at zero because the substrate was absent, here because it is in excess — so the guard
    asserts the RATIO rather than the zero. A zero-only assertion would keep passing if a later
    beat changed the dose arithmetic and the term went live at some other magnitude.
    """
    compiled = compile_scenario(_COPPER_WINE)
    capacity = _COPPER_DOSE_GPL * compiled.parameters["copper_h2s_binding"].value
    out = simulate_scheduled(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
        events=compiled.events,
        t_eval=np.linspace(*compiled.t_span_h, 401),
    )
    at_dose = int(np.argmin(np.abs(out.t - _COPPER_DOSE_DAY * 24.0)))
    h2s_present = float(out.y[compiled.schema.slice("h2s"), at_dose].sum())
    assert capacity > 100.0 * h2s_present, (
        f"copper capacity {capacity:.3e} is no longer in large excess of the {h2s_present:.3e} "
        "g/L of H2S at the dose, so copper_h2s_binding has become reachable. Price it and move "
        "its row — it stopped being the losing side of the min."
    )
    u = compiled.parameters["copper_h2s_binding"].uncertainty
    lo = _endpoint(_compiled_with(_COPPER_WINE, "copper_h2s_binding", u.low), "h2s")
    hi = _endpoint(_compiled_with(_COPPER_WINE, "copper_h2s_binding", u.high), "h2s")
    assert lo == hi, "the whole band must be bit-identical while it loses the min"


@pytest.mark.parametrize("name", _SUBSUMED_FIT + _WORT_ACIDS + _DOSE_VERB)
def test_the_priced_names_are_still_undrawn(name):
    """A RED here means the gap was CLOSED — delete the row and say so in the record.

    D-233 §6's idiom, and D-237 shipped it for ``burst_antioxidant_initial`` alone. Wiring any of
    these into a Process's ``reads`` moves every ensemble that reaches it, so the change must
    arrive with its own measurement rather than as a green diff.

    **Scored under the name's OWN oxidative wiring, and the first version was not.** D-240's
    Arm C added ``burst_antioxidant_initial`` to ``AntioxidantBurstOxidation.reads`` — the exact
    repair this guard exists to catch — and the guard stayed **green**, because it compiled the
    default ``direct`` set, where that Process is not wired into the medium at all (D-147) and
    its ``reads`` never enter ``_schedule_reads``. A guard that scores the one configuration in
    which its subject cannot exist forbids nothing
    [[feedback-a-guard-must-be-scored-where-its-subject-lives]].
    """
    scenario = BEER if name.endswith("_beer") or name.endswith("_wort") else WINE
    oxidative = _SEEDED_SLOT.get(name, (None, None, "direct"))[2]
    compiled = compile_scenario(scenario, oxidative=oxidative)
    if name not in compiled.parameters:
        pytest.skip(f"{name} is not in {scenario.medium}'s parameter set")
    sampled = set(
        _resolve_sample_names(
            compiled.process_set, compiled.parameters, None, None, compiled.events
        )
    )
    assert name not in sampled, (
        f"{name} is now drawn. If a beat wired it into a Process's reads AND gave y0_for_member a "
        "rule for the slot it seeds, that is the repair D-240 priced and declined to ship: delete "
        "this row and record the movement. If only the declaration landed, the ensemble is now "
        "drawing a name that cannot reach y0, which is worse than the gap."
    )


@pytest.mark.parametrize("name", _REPAIRED_AT_D241)
def test_the_six_repaired_seeds_really_are_drawn_now(name):
    """The mirror of :func:`test_the_priced_names_are_still_undrawn` — same teeth, opposite sign.

    D-241 drew these six and gave each a ``y0_for_member`` rule. A RED here means one of them
    went back to being invisible to the sampler, which would silently restore the defect D-240
    priced: a banded seed whose uncertainty no reported spread contains.

    **Scored under the name's OWN wiring**, for D-240 Arm C's reason. ``burst_antioxidant_initial``
    is drawable only under ``oxidative="direct_burst"``: elsewhere
    :func:`_resolve_burst_antioxidant_seed`
    zeroes the slot, the seed rule's equality guard therefore declines, and the name correctly
    stays out of ``seed_reads``. That is not a gap — it is the D-147 condition being honoured for
    free, and a test that scored it under ``direct`` would report a repair that is working as a
    repair that vanished.
    """
    scenario = BEER if name.endswith("_beer") or name.endswith("_wort") else WINE
    oxidative = _SEEDED_SLOT.get(name, (None, None, "direct"))[2]
    compiled = compile_scenario(scenario, oxidative=oxidative)
    if name not in compiled.parameters:
        pytest.skip(f"{name} is not in {scenario.medium}'s parameter set")

    assert name in compiled.seed_reads, (
        f"{name} is no longer in seed_reads, so no ensemble can draw it and no member re-seeds "
        "the slot it fills. Both halves are derived from one list, so this is a reverted repair "
        "rather than a half-repair — restore the y0_for_member rule."
    )
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
    assert name in sampled, f"{name} is in seed_reads but the resolver dropped it"
    # …and the rule reaches y0, which is the half `seed_reads` alone cannot promise.
    build = compiled.y0_for_member()
    assert build is not None
    values = dict(compiled.parameters.resolve())
    values[name] = compiled.parameters[name].uncertainty.high
    assert not np.array_equal(build(values), compiled.y0), (
        f"{name} is drawn but a member at its band edge gets the compiled y0 — the declaration "
        "half landed without the seed half, which D-240 §10 calls worse than the gap."
    )


@pytest.mark.parametrize("oxidative", ("direct", "cascade"))
def test_the_burst_seed_is_inert_in_the_wirings_that_keep_it_in_this_census(oxidative):
    """The measurement behind the WIRING-CONDITIONAL verdict, and the reason it is not a gap.

    D-241 drew this seed under ``direct_burst`` and left it undrawn under the other two. That
    looks like a half-repair and is not: where :class:`AntioxidantBurstOxidation` is not wired,
    D-147 zeroes the slot after the pack, so the parameter is read at compile and its value is
    thrown away. There is no uncertainty for a sampler to be missing.

    Asserted on the WHOLE of ``y0`` rather than on the one slot, because the slot being 0 at both
    edges is also what a seam that stopped reading the name entirely would produce, and the
    stronger statement is the one that stays true: nothing about the compiled state depends on
    this parameter here. The 50x band makes it a real test — a leak of any size fails it.
    """
    u = compile_scenario(WINE).parameters["burst_antioxidant_initial"].uncertainty
    assert u.high / u.low > 10.0, (
        "the band narrowed; this null is only interesting while it is wide"
    )
    lo = _compiled_with(WINE, "burst_antioxidant_initial", u.low, oxidative=oxidative)
    hi = _compiled_with(WINE, "burst_antioxidant_initial", u.high, oxidative=oxidative)
    assert np.array_equal(lo.y0, hi.y0), (
        f"the burst seed now reaches y0 under `{oxidative}`. If a beat wired the burst Process "
        "into this set, the name belongs in seed_reads here too — add the row and price it."
    )
    assert "burst_antioxidant_initial" not in lo.seed_reads, (
        "the seed rule fired on a wiring whose slot D-147 zeroes, which would hand members a "
        "burst pool the compiled build cannot spend"
    )
    # The positive control the two nulls owe: the SAME harness, same wiring, does move the state
    # for a name that is genuinely live here [[feedback-a-null-result-needs-a-positive-control]].
    v = compile_scenario(WINE).parameters["dms_potential_initial"].uncertainty
    a = _compiled_with(WINE, "dms_potential_initial", v.low, oxidative=oxidative)
    b = _compiled_with(WINE, "dms_potential_initial", v.high, oxidative=oxidative)
    assert not np.array_equal(a.y0, b.y0), "the harness moved nothing at all; the nulls are void"
