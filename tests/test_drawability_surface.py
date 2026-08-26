"""Drawability: which banded seam parameters the ensemble sampler can actually reach.

Decision D-159, correcting D-157's Finding 1.

D-157 measured that a large slice of the compile seam's banded parameters is never drawn,
and named three consumption modes. Its table is hedged correctly ("not drawn by any
scenario **tried**"); its headline is not — it says 66 bands are *structurally*
undrawable. Measured against the sampler's own mechanism rather than a scenario census,
the structural count is **61**: five of the 66 (``E_a_pof``, ``k_pof_decarb``,
``k_pof_decarb_ferulic``, ``biomass_carrying_capacity``, ``k_iso_alpha_loss``) are declared
in some Process's ``reads`` and are merely *scenario*-inert — a 14th scenario draws them.

**What is pinned here, and what is deliberately not.**

The census is *not* pinned. Deriving "structural" from the declarations and comparing it to
what a run draws is vacuous: ``_resolve_sample_names`` narrows by ``_schedule_reads``, which
is itself the union of declared ``reads``, so ``structural ⊆ undrawn`` is a theorem about
set arithmetic rather than a fact about the code (the D-108/D-109 shape). Freezing 61 names
into the suite would also go red on every new parameter. Counts live in the D-159 record.

What *is* pinned is the **mechanism**, on one exemplar per structural class, and it takes
**two** assertions per exemplar, because either alone is satisfiable for the wrong reason:

1. **Consequential** — the value demonstrably reaches the model, at a named site. This is
   the half that a "the trajectory did not move" test cannot supply.
2. **Unreachable** — forcing the name into the sample via ``only=`` moves nothing.

**Why both halves are required, measured rather than argued.** The first draft of this test
was half 2 alone: force the name, assert the members are identical. It is unsound. Five
genuinely Process-read names froze at *exactly* ``max|dy| == 0.0`` in the same harness —
``k_copper_multiplier`` (no ``add_copper`` in *this* scenario — see the note below),
``ethanol_tolerance`` (``EthanolInactivation``'s term is clamped off below tolerance), and
``k_so2_oxidation`` / ``k_browning_base`` / ``k_ethanol_oxidation`` (supply-limited: no O₂
in the scenario). "Frozen" conflates *unreachable by the sampler* with *reachable but zero
here*, so half 2 alone would pass for the wrong reason on any supply-limited parameter —
D-155's decoration, and D-157's silent denominator, in a third costume.

:func:`test_the_same_harness_moves_a_process_read_parameter` is the positive control for
half 2, and it is not optional: without it "identical" cannot be told from "the harness
perturbed nothing". It runs at the *same* ``n_members``, span and settings as the exemplars.

**Two classes are deliberately absent as exemplars.**

* The **acidbase pKa/SO₂-binding set is NOT structural** and must never be added here.
  ``ph_of_state`` (``acidbase.py:390``) reads the pKa set at runtime for Processes in five
  kinetics modules, and ``free_acetaldehyde``/``bisulfite_so2_at_ph`` read the four
  carbonyl-bisulfite constants the same way. Until D-160 no Process declared any of them, so
  the sampler never perturbed them — that was a *defect* (flagged by D-159), not an example
  of one. **D-160 shipped the declarations** (``PH_SYSTEM_READS``/``SO2_BINDING_READS``,
  13 names over the two helper families), and every one is drawn by default now, so the
  class is drawable and adding it here would fail half 2 for the correct reason.
  Forcing it in isolation still moves 24 of 94 state slots — but D-161 measured what that
  is worth against a *full* sampled set and found that only one of those slots survives the
  competition, via ``ester_hydrolysis``'s ``pH_ref_ester_hydrolysis`` term. The per-slot
  numbers live in that record, deliberately: nothing here asserts them, and a measured
  figure in an unpinned docstring is the same thing that went stale above.
* ``ethanol_inhibition_exponent`` is undrawable for a fourth reason D-157 did not name —
  its ``EthanolInhibition`` is deliberately not wired into any medium (D-13: it would
  double-count ``EthanolInactivation``; ``media.py:1328`` retains the class "for
  optional/strain use"). Half 1 fails for it by construction, because nothing consumes it.
  Wiring it in is a legitimate future change that *should* draw the name.
"""

import numpy as np
import pytest

from fermentation.runtime import simulate_ensemble
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario.compile import _CLOSURES

#: Ferment/aging boundary and run end, in days. Short enough that five ensembles stay a
#: couple of seconds each, long enough that the positive control separates decisively.
FERMENT_DAYS = 12.0
RUN_DAYS = 60.0

#: The oak charge the verb exemplar doses. ``fill_number`` is left at its default 1 so the
#: ceiling is the UNSCALED ``oak_gpl · yield`` product (D-91 discounts reused fills).
OAK_GPL = 4.0
OAK_TOAST = "medium"


def _scenario(closure: str = "screwcap") -> Scenario:
    """A wine that ages under a named closure with an oak charge, so both exemplars compile.

    ``begin_aging`` is what enables ``OakExtraction`` and the closure-ingress Process, so a
    scenario without it would leave both exemplars inert for an uninteresting reason.
    """
    return Scenario(
        name="d159-drawability",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.25,
            "amino_acids_gpl": 1.0,
            "anthocyanin_gpl": 0.5,
            "tannin_gpl": 1.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        closure=closure,
        interventions=[
            Intervention(day=FERMENT_DAYS, action="begin_aging"),
            Intervention(
                day=FERMENT_DAYS,
                action="add_oak",
                params={"oak_gpl": OAK_GPL, "toast": OAK_TOAST},
            ),
        ],
        duration_days=RUN_DAYS,
    )


def _forced(compiled, name: str, n_members: int = 2):
    """Force ``name`` into the sample via ``only=`` and return the ensemble.

    ``only`` bypasses ``_schedule_reads`` entirely (``chosen = set(only)``), which is the
    whole point: it is the strongest reach the sampler API offers, so an invariance under
    it is a statement about consumption rather than about declaration.
    """
    return simulate_ensemble(
        compiled.process_set,
        compiled.parameters,
        compiled.y0,
        compiled.t_span_h,
        n_members=n_members,
        seed=0,
        only=[name],
        events=compiled.events,
    )


def _assert_the_draw_actually_happened(ensemble, name: str) -> None:
    """The guard against a vacuous invariance: the name was sampled, and the values differ.

    ``_resolve_sample_names`` intersects with ``parameters.names``, so a typo or a renamed
    parameter would silently yield an EMPTY sampled set — under which every member is the
    nominal run and "the trajectory did not move" is true for no interesting reason.
    """
    assert name in ensemble.sampled_names, f"{name} never entered the sample"
    drawn = [member[name] for member in ensemble.member_params]
    assert len(set(drawn)) == len(drawn), f"{name} drew the same value twice: {drawn}"


def _members_are_identical(ensemble) -> bool:
    """Exact equality, not a tolerance — the solver is deterministic.

    Identical parameters into an identical RHS give bit-identical BDF output, so a
    consumed-before-sampling parameter yields exactly 0.0 difference. A tolerance here
    would let a genuine small dependence (the pKa route moves some slots by ~1e-5) pass
    as invariance.
    """
    return bool(np.array_equal(ensemble.members[0], ensemble.members[1]))


# -- class 2: compile-consumed (``_closure_otr``) ------------------------------


@pytest.mark.parametrize("closure", _CLOSURES)
def test_the_closure_otr_value_is_consumed_into_state_before_any_sampling(closure):
    """Half 1 — ``otr_<closure>`` reaches the model as a STATE slot, fixed at compile.

    ``compile.py`` resolves ``scenario.closure`` once, through ``_closure_otr``, into
    ``y0``'s ``closure_otr`` slot. Asserting that identity says exactly where the value
    went, and — unlike a trajectory check — cannot be satisfied by a term that happens to
    be zero in this scenario. It is also the reason the sampler cannot reach it: the
    sampler draws *parameters*, and by the time it runs this number is *state*.

    Parametrized over the whole menu so a closure added to ``_CLOSURES`` without a matching
    ``otr_`` parameter, or wired to the wrong one, fails here.
    """
    compiled = compile_scenario(_scenario(closure))
    slot = compiled.y0[compiled.schema.slice("closure_otr")]
    assert slot.shape == (1,)
    assert float(slot[0]) == compiled.parameters.value(f"otr_{closure}")


def test_forcing_the_closure_otr_band_into_the_sample_cannot_move_the_run():
    """Half 2 — drawn across its band, ``otr_screwcap`` changes nothing.

    Paired with the compile-time identity above, this is the full claim: the parameter is
    consequential *and* the sampler cannot reach it. Neither half states it alone.

    This is also why D-157's Finding 2 — the closure menu's bands contradict D-136's
    ascending-OTR ordering on ~42 % of joint draws — cannot go live inside a run: one
    ``scenario.closure`` resolves to one state slot, and no ensemble perturbs it.
    """
    compiled = compile_scenario(_scenario("screwcap"))
    ensemble = _forced(compiled, "otr_screwcap")
    _assert_the_draw_actually_happened(ensemble, "otr_screwcap")
    assert _members_are_identical(ensemble)


# -- class 3: verb-consumed (``_verb_add_oak``) --------------------------------


def test_the_oak_yield_is_consumed_into_the_ceiling_slot_by_the_verb():
    """Half 1 — ``oak_yield_vanillin_medium`` reaches the model through the EVENT.

    ``add_oak`` computes ``oak_gpl · oak_yield_<compound>_<toast>`` and writes it to the
    compound's set-and-hold ceiling slot (D-77). Applying the compiled event's own
    ``mutate`` to ``y0`` reproduces that product exactly, with no integration — so this is
    a statement about the verb, not about whether ``OakExtraction`` happens to be running.

    The value travels as a state jump; ``params`` never carries it into a derivative, which
    is why sampling the parameter afterwards is a no-op.
    """
    compiled = compile_scenario(_scenario())
    (oak_event,) = [e for e in compiled.events if e.mutate is not None and "oak" in e.label]
    mutate = oak_event.mutate
    assert mutate is not None
    after = mutate(compiled.schema, compiled.y0, compiled.param_values)

    ceiling = float(after[compiled.schema.slice("vanillin_ceiling")][0])
    expected = OAK_GPL * compiled.parameters.value(f"oak_yield_vanillin_{OAK_TOAST}")
    assert ceiling == pytest.approx(expected, rel=1e-12)
    # ...and it is a real jump, not a slot that was already there.
    assert float(compiled.y0[compiled.schema.slice("vanillin_ceiling")][0]) == 0.0


def test_forcing_an_oak_yield_band_into_the_sample_cannot_move_the_run():
    """Half 2 — drawn across its band, ``oak_yield_vanillin_medium`` changes nothing.

    The measured half of D-157's cleanest evidence: firing the intervention verbs left the
    drawn count unchanged at 152 of 247. The verb reads ``.value`` off the ParameterSet at
    compile time, so the per-member draw never reaches it.

    **BOTH numbers are prose, both went stale, and D-182 caught them.** Re-measured at D-182
    on this module's own ``_scenario()``:

    * **DENOMINATOR — compile-seam distinct varying: 280.** 247 at D-157 → 246 at D-172 (three
      varying bands retired, two added) → 257 at D-179 (beer's five pKas) → 266 at D-180 (five
      ``*_typical_wort`` seeds, four ``Y_*_sugar_beer`` yields) → 277 at D-181 (three seeds,
      three floors, one rate, four pKas) → **280 at D-182** (``pKa_carbonic_1``,
      ``H_co2_beverage``, ``vant_hoff_co2_solubility``; ``T_ref_co2_solubility`` and
      ``peptide_buffer_capacity_beer`` are pinned zero-width and so are not counted).
    * **NUMERATOR — drawn by a default ensemble on this scenario: 185, NOT 152.** The claim
      that "the 152 is unaffected each time, because none of the added names is drawn by this
      scenario" **was false from D-179 onward** and nobody re-ran it. D-179 made
      ``PH_SYSTEM_READS`` the *union* across media **on purpose**, which puts beer's pKas into
      a WINE ensemble's sampled set — D-179 said so in its own record, and this line kept
      asserting the opposite in prose. D-182 added four more names to that same tuple.

    Reproduce both with ``ps[n].uncertainty.high > ps[n].uncertainty.low`` over
    ``compiled.parameters.names``, and ``len(simulate_ensemble(...).sampled_names)`` on an
    unforced ensemble.

    NOTE these are PROSE, not asserts, which is why retiring three parameters turned nothing
    red — recorded as a live gap in the surface guards at D-172 §8 and STILL OPEN. Two beats
    of silent drift on the numerator is what that gap costs, and it is the reason the standing
    instruction is **re-measure it, do not cite it**. Deliberately still not converted to an
    assertion — pinning it would fire on every future parameter addition regardless of
    relevance, which is a chore rather than a guard, and D-172 §8's gap is about the surface
    tests generally, not this one line.
    """
    compiled = compile_scenario(_scenario())
    ensemble = _forced(compiled, "oak_yield_vanillin_medium")
    _assert_the_draw_actually_happened(ensemble, "oak_yield_vanillin_medium")
    assert _members_are_identical(ensemble)


# -- the positive control ------------------------------------------------------


def test_the_same_harness_moves_a_process_read_parameter():
    """The control: an identical ``only=`` call on a Process-read name DOES move the run.

    Without this, both invariance tests above are indistinguishable from a harness that
    perturbs nothing — the failure mode that made D-157's first assertion sweep report
    "5 of 5 clean" on a denominator it never measured.

    ``mu_max`` is the pick because it drives the primary ferment from ``t0``, so it cannot
    be silently supply-limited the way the aging rate constants are: the same harness
    reports exactly 0.0 for ``k_so2_oxidation`` / ``k_browning_base`` /
    ``k_ethanol_oxidation`` (no O₂ here), for ``k_copper_multiplier`` (see below) and for
    ``ethanol_tolerance`` (its term is clamped off below tolerance) — every one of them
    genuinely Process-read. A control drawn from that group would have confirmed the wrong
    thing.

    **``k_copper_multiplier``'s freeze was SCENARIO-bound, and D-191 lifted it.** D-159
    attributed it to "the ``copper`` slot is never written", which was true then and is the
    right diagnosis: ``f_copper = 1 + k·(copper − copper_typical)`` is identically 1 for
    *any* k when copper sits at ``copper_typical``, so the multiplier cannot matter until
    something displaces the slot. Since D-191 the ``add_copper`` verb does. Measured in this
    same harness (`d191-residual-copper/probe3_thaw.py`), forcing the name gives
    ``max|dy|`` = **0.0** with no fining, **6.6e-4** with a 0.5 mg/L fining, and **1.3e-2**
    with fining plus dosed O₂ — the copper displacement is what thaws it and the oxidant
    only amplifies. This scenario has no ``add_copper``, so the 0.0 above still holds here;
    it is no longer a statement about the model.
    """
    compiled = compile_scenario(_scenario())
    ensemble = _forced(compiled, "mu_max")
    _assert_the_draw_actually_happened(ensemble, "mu_max")
    assert not _members_are_identical(ensemble)
