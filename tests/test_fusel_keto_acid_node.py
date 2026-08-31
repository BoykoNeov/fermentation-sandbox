"""The FUSEL side of the keto-acid node — the findings D-109 measured instead of building.

D-107 built the excreted ``alpha_ketobutyrate`` pool and left a work-list item that D-108 promoted
to "the largest open item": re-base propanol on that pool, so the *genuine propanol-vs-sotolon
competition* — over α-ketobutyrate, which propanol IS the decarboxylation of — becomes expressible.
**D-109 measured that item's premise before building it, and the premise is wrong twice over.**
This suite pins the measurements, because every one of them is a claim a future beat could quietly
invalidate.

**1. The pool D-107 chose FOR sotolon is the wrong pool for propanol — D-49's test, applied
symmetrically.** :class:`~fermentation.core.kinetics.aging.SotolonAldolCondensation` selects the
*excreted, extracellular residual* precisely because its aldol runs in a sealed bottle where no
intracellular pool can reach. Propanol is made **intracellularly, by living yeast, during active
fermentation** — it is exactly D-49's *flux intermediate*, the thing that module rejected as
acetaldehyde's precursor for the same reason. One pool cannot be both. The principle convicts on
its own (no literature needed — the sibling Processes already carry the argument), and
:func:`test_the_excreted_pool_cannot_supply_propanol` adds the arithmetic: the pool is not merely
the wrong *identity*, it is the wrong *size* by ~3×, so the re-base is **infeasible**, not just
mis-attributed.

**2. The work-list's "the competition is GONE" is CORRECT, and D-109 nearly "corrected" it.**
:class:`~fermentation.core.kinetics.keto_acids.AlphaKetobutyrateExcretion`'s rate is flux-only;
threonine's gate re-routes the *carbon source*, not the *rate*. So the competition over threonine
is **exactly zero, by construction and on purpose** (gating the rate on threonine would kill
sotolon in a threonine-free wine — the D-104 canary). D-109's first probe measured a 0.42%,
monotone, correctly-signed sotolon response to ``k_propanol`` and read it as "the competition is
present but small". It is the **sugar ledger**: propanol's de-novo carbon leaves ``S``, which every
flux-linked rate reads. See :func:`test_alpha_kb_production_is_exactly_threonine_independent` —
the mutation test is the only reason that write-up did not ship.

**3. The promised payoff is REAL, and the number that kills the shortcut is the number that says
so.** Propanol's molar demand is **1.36× the total α-KB the pool ever excretes** (D-245; it was
~2.8× until D-244's corrected yield halved the biomass propanol rides on, and the re-record is in
:func:`test_the_excreted_pool_cannot_supply_propanol`) — which is why it cannot be drawn from that
pool (finding 1), and equally why a *correctly placed* node would matter: propanol's 2-KB demand
still **exceeds the excretion flux itself**, so partitioning that flux honestly would couple
propanol and sotolon **materially**. The item is not dissolved by this beat;
it is **relocated** — off the excreted pool and onto the intracellular partition, where it needs
the milestone rather than a shortcut. (Rests on ``k_alpha_kb_excretion``, an author estimate ⇒ an
order-of-magnitude claim, not a calibration.)

**The claim is deliberately the RATIO, not "propanol is the node's dominant sink".** 2-KB's
committed anabolic route is **isoleucine biosynthesis** (ILV2 → KMV → isoleucine) — that is *why*
the cell makes 2-KB at all; propanol and excretion are both overflow off it. This model carries no
KMV, so propanol's share of *total* 2-KB synthesis is **unmeasured** and could well be the smaller
one. Propanol-vs-excretion is what the measurement supports and is all the argument needs; the
"dominant sink" phrasing was a scope error caught in review, in the beat whose subject is exactly
that (and it would have jarred against this file's own milestone scoping, which names KMV).

**Why the fusel-side node is a PARTITION and not a pool (the scoping result).** D-49's physics
says the intracellular keto acid is a vanishing pool carrying an enormous flux — i.e.
quasi-steady, ``synthesis == Σ consumption`` at every instant. A quasi-steady node is a **flux
partition**, not a state variable. So the two nodes are different in *kind*: the excreted
keto-acids are pools (they persist, they are measured, they bind SO₂), and the fusel-side node is
a partition of the sourcing. That is why the ``FuselAminoAcidReroute``/``PrecursorNonEhrlichFates``
sourcing layer — not the producer — is where the milestone belongs, and it is why nothing here
needs a sixth state slot.
"""

from __future__ import annotations

import pytest

from fermentation.core.chemistry import (
    CARBON_ATOMS,
    M_ALPHA_KETOBUTYRATE,
    M_PROPANOL,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS, SECONDARY_FUSEL_ROUTES
from fermentation.core.kinetics.precursor_fates import non_ehrlich_fraction_param
from fermentation.core.tiers import Tier
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import (
    Intervention,
    Scenario,
    TemperaturePoint,
    amino_acid_dose_nitrogen_mgl,
    compile_scenario,
)

#: Crépin *et al.* 2017's own must: 180 mg N/L, 28 °C. Every share below is quoted against that
#: paper's numbers, so the probe must run on a COMMENSURATE must — D-104's lesson, where a
#: ~470 mg N/L probe flattered the model against a 180 mg N/L source.
_CREPIN_YAN = 180.0
_CREPIN_TEMP = 28.0
_FERMENT_DAYS = 14.0
_AGING_DAYS = 720.0

#: Crépin's measured de-novo share of 2-ketobutyrate (the 81 of "19% exogenous / 81% newly
#: synthesised"), and Rollero's independent ">90% from the carbon central metabolism". The model
#: is asserted against the WEAKER of the two, so the test cannot be broken by the disagreement
#: between them — recorded as two bands, never averaged (D-103).
_SOURCED_DE_NOVO_FLOOR = 0.80

#: The precursor excluded from the floor above. **It is no longer excluded for being unsourced**
#: (D-117 sourced it: Minebois 2025, U-13C phenylalanine, ``f`` = 0.975, plausible) but because the
#: FLOOR's own citations — Crépin's 2-KB, Rollero's CCM — never measured **phenylpyruvate**. The
#: name is deliberately not ``_UNSOURCED_*`` any more: that would keep asserting a dead reason.
#: See :func:`test_2_phenylethanol_carries_no_sourced_de_novo_floor` for both legs.
_FLOOR_EXCLUDED_PRECURSOR = "phenylalanine"

#: The four alcohols whose keto acids the floor's own sources actually describe ([13C] leu/ile/val,
#: thr → Crépin's 2-KB, Rollero's CCM). Derived from the registry rather than hand-listed so a sixth
#: alcohol cannot silently inherit a floor no source covers (D-104: **a cited number binds only the
#: SET it describes** — and this test is where that rule is easiest to break, because parametrizing
#: over all of ``FUSEL_SPECS`` looks more thorough, not less). D-117 is the proof it is easy: 2-PE
#: now *passes* the floor by a wide margin and still may not be admitted to it on someone else's
#: citation.
_SOURCED_FUSEL_SPECS = tuple(
    s for s in FUSEL_SPECS if s.precursor_amino_acid != _FLOOR_EXCLUDED_PRECURSOR
)

#: Why the sourced de-novo FLOOR is missed since D-244 -- the cause MEASURED at D-245, where
#: D-244 could only assert it. One knob reproduces it: hold this fixture exactly as it is (same
#: must, same 1.0 g/L dose, same declared 405.4 mg N/L total) and sweep ONLY
#: ``biomass_N_fraction`` from the value the pre-D-244 evaluation point produced (0.0578, the fit
#: read at the declared 180) to the shipped one (0.1068, the fit read at the total the must
#: actually holds, held at Coleman's 350 edge). Every threshold crosses monotonically inside that
#: sweep, so the must's 2.25x nitrogen richness is a CONSTANT here and not a cause.
#:
#: **Which leg moved: the DENOMINATOR.** Consumed threonine (0.06700), valine (0.03700) and
#: leucine (0.03200) are identical to 5 dp at every f_N -- those pools exhaust whatever the
#: biomass, so the amino-acid draw is supply-limited and PINNED. The alcohol totals halve with
#: biomass, and the share rose entirely through its denominator.
#:
#: **And the same correction removed a ~2x OVER-production**, which is why this guard used to
#: pass: propanol 49.8 -> 27.0, isobutanol 68.4 -> 37.1, isoamyl 353.5 -> 191.5, 2-PE 59.5 ->
#: 32.3 mg/L, against per-molecule anchors of ~24 / 33.0 / 172 / 28.7 mg/L (those anchors describe
#: a typical UNDOSED must, so they do not strictly apply here -- the move is ~2x -> ~1.1x them,
#: not "now lands on the anchor"). The surplus was de-novo, sugar-sourced carbon sitting in this
#: ratio's denominator, so the floor was being cleared by an over-production rather than by the
#: supply structure the guard names.
#: **CLOSED AT D-248 — kept as the history of a xfail that ran for four records, not as a live
#: reason.** The gap was never a fusel parameter and never the availability gate: it was the
#: BIOMASS DENOMINATOR. Assimilation's only route into biomass nitrogen ran at
#: ``psi*gate*f_N*base_dx``, strictly below growth's own draw, so ammonium could only fall — and
#: when it reached zero, growth's Monod shut growth off and the swap stopped with it, freezing the
#: must with 40.8 % of its assimilable nitrogen unconsumed against Crepin's measured 0.2 %.
#: Un-coupling uptake from the growth rate (D-248) builds 98.4 % of the Coleman yield the compile
#: seam itself installs instead of 61.6 %, and propanol goes **0.7744 -> 0.8062** on this fixture
#: and 0.7963 -> 0.8784 on Crepin's own must. Nothing about threonine moved: this Process draws no
#: precursor at all. Retained verbatim below because the archive is append-only and because a
#: future beat reading "we lack the medium" or "the gate is mis-scaled" in an old record needs to
#: find where both were measured and set aside.
_D245_DE_NOVO_FLOOR_GAP_CLOSED_AT_D248 = (
    "D-245 (measuring D-244): propanol draws 77.4 % of its carbon de novo against the sourced "
    "80 % floor. MEASURED cause, not asserted: the corrected yield evaluation point roughly "
    "halves biomass, the alcohol total halves with it, and the amino-acid draw does not move at "
    "all (threonine exhausts either way), so the share rises through its denominator. The same "
    "correction removed a ~2x over-production that had been clearing this floor for the model. "
    "STRICT: the floor is a sourced target (Crepin 81 % newly-synthesised 2-KB; Rollero >90 % "
    "CCM) and is NOT to be lowered to fit the model. "
    "D-246 SOURCED THE MEDIUM AND IT DOES NOT CLOSE THIS: on Crepin's own must (its Data Set S1, "
    "measured, mean of 14) propanol reads 0.7963 against the 0.80 floor -- the commensurability "
    "violation D-244 section 6 recorded was 86 % of the gap and not the whole of it. "
    "D-247 MEASURED THE REST AND IT IS STILL UNATTRIBUTED: rescaling depletion_gate to the pool "
    "it actually gates -- the composition correction D-246 section 6 described, each share "
    "re-referenced to the declared must with the spectrum sum preserved -- moves propanol "
    "0.796275 -> 0.796017, i.e. 6.9 % of the remaining gap in the WRONG direction. D-246's probe "
    "cleared the floor only because scaling K_amino_acids uniformly is a LEVEL change, "
    "bit-identical to scaling all eight must_aa_fraction_* by that same factor, and nothing "
    "sources a smaller availability constant. The gate's reference is an unsourced modelling "
    "device either way (D-100 declined per-species Michaelis constants as the D-98 trap), so the "
    "repair is REFUSED and only literature half-saturations can settle it. See "
    "tests/test_defined_media.py, which pins all of that; do NOT re-cite 'we lack the medium'."
)

#: Why the two D-120 guards are STRICT xfails. Separate from the floor gap above because the
#: content is different: D-120 refused a de-novo cap for isoamyl on TWO measured legs, and D-245
#: measured both of them gone. Direction: every alcohol used to attribute LESS to amino acids
#: than Minebois measures, so a one-directional ceiling would have moved all three the wrong way;
#: isoamyl now reads 5.42 % against her 5.34 % and isobutanol 9.47 % against her 8.78 %.
#: Instrument: the cap was inert because phenylalanine exhausted with or without it, and it no
#: longer does (12.8 % survives), so the shipped cap now moves 2-PE's realised share by 12.7 %.
#: **BOTH LEGS RE-MEASURED AT D-248, AND THE DIRECTION ONE IS BACK.** The two guards this reason
#: xfailed are green again: isoamyl reads 0.873x and isobutanol 0.927x Minebois's in-study shares
#: on this fixture (1.01x / 0.98x on her own must), so nothing over-attributes and D-120's refusal
#: of a de-novo cap stands on its direction leg exactly as it originally did. The INSTRUMENT leg is
#: thinner rather than restored: the cap's bite on 2-PE's realised share shrinks 12.7 % -> 4.8 %
#: without reaching the inertness D-120 measured, so
#: ``test_the_de_novo_cap_is_inert_where_the_precursor_exhausts`` remains a strict xfail. Kept
#: verbatim as the history of what D-245 measured and D-248 undid.
_D245_D120_LEGS_GONE_DIRECTION_BACK_AT_D248 = (
    "D-245: D-120's refusal of a de-novo cap rested on two measured legs and the corrected yield "
    "evaluation point removed both. (1) DIRECTION: the model under-attributed to amino acids "
    "against Minebois's in-study shares; it now over-attributes -- isoamyl 5.42 % vs 5.34 %, "
    "isobutanol 9.47 % vs 8.78 %. (2) INSTRUMENT: the cap measured inert because phenylalanine "
    "exhausted either way; it no longer exhausts (12.8 % left) and the shipped cap now moves "
    "2-PE's realised share 1.603 % -> 1.400 %, a 12.7 % relative bite. STRICT: this is D-120's "
    "own tripwire firing as written, so the de-novo-entry build genuinely re-opens -- but the "
    "parameter it needs is unsourced for isoamyl, and deriving one from the model's own "
    "abundances is refused (D-206). "
    "D-246 SCORED BOTH LEGS ON MINEBOIS'S OWN MUST AND THEY GET WORSE, NOT BETTER: isoamyl "
    "0.0926 and isobutanol 0.1474 against her 0.0534 and 0.0878, i.e. ~1.7x each where the "
    "fixture read ~1.1x. That also retires D-245's own caveat that the isoamyl trip (1.5 % "
    "relative) sat inside this harness's cap-window systematic -- 73 % over is an order of "
    "magnitude outside it, so the over-attribution is a property of the model at her nitrogen "
    "rather than an artefact of a richer must. Pinned in tests/test_defined_media.py."
)

#: The routes that ALSO eat the speciated precursors. Disabled where the ``f : (1−f)`` split
#: invariant is used to size the Ehrlich draw from a state difference, so that invariant is
#: EXACT rather than approximately true (they measure inert in this regime anyway — but "measures
#: inert" is a reason to verify, not to assume: D-106).
#: ``strecker_degradation`` is the DIRECT oxidative set's name, and the direct set is the default
#: wiring (D-141), so this is the name ``disable`` must find. Under the cascade the same chemistry
#: is ``quinone_strecker_degradation`` — if the cascade ever becomes the default, this entry moves
#: with it.
#:
#: **Recorded because this module was a Gate 3 miss.** D-139 enumerated the expected reds by
#: subsystem (aging, media, closure) and never looked for a cross-domain consumer of an aging
#: Process *name*, so when D-141 briefly made the cascade the default, five tests here and five in
#: ``test_fusel_catabolic_shape.py`` went red unpredicted. They failed LOUDLY only because ``_run``
#: calls ``disable(name)`` rather than filtering — see the note there, now vindicated rather than
#: merely careful. A filter would have silently disabled nothing and quietly changed what the
#: ``f : (1-f)`` invariant was measuring.
_OTHER_PRECURSOR_CONSUMERS = (
    "alpha_kb_excretion",
    "maillard_strecker",
    "strecker_degradation",
    "autolytic_mercaptan",
)


def _scenario(*, o2_mgl: float = 20.0, aging: bool = True) -> Scenario:
    interventions = []
    duration = _FERMENT_DAYS
    if aging:
        interventions = [
            Intervention(day=_FERMENT_DAYS, action="begin_aging"),
            Intervention(day=_FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": o2_mgl}),
        ]
        duration = _FERMENT_DAYS + _AGING_DAYS
    initial: dict[str, float] = {
        "brix": 24.0,
        "yan_mgl": _CREPIN_YAN,
        "pitch_gpl": 0.25,
        "amino_acids_gpl": 1.0,
    }
    # D-244: ``yan_mgl`` is the must's TOTAL assimilable nitrogen and the amino-acid dose is
    # carved OUT of it. This fixture was authored when the two channels ADDED, so it declares
    # the sum -- which leaves its pitch state bit-for-bit and moves only the point Coleman's
    # yield fit is evaluated at, which is the defect D-243 found.
    initial["yan_mgl"] += amino_acid_dose_nitrogen_mgl(initial)
    # AND IT IS ALSO WHERE D-244 FOUND A LIVE COMMENSURABILITY VIOLATION, recorded and NOT
    # repaired here. The comment on the YAN constant above forbids scoring a richer must than
    # the paper's, citing D-104. Because the two nitrogen channels ADDED before D-244, this
    # probe was doing exactly that -- and the migration keeps it doing so.
    #
    # **THE REASON THIS USED TO GIVE IS DEAD (decision D-254).** It read: the alternative
    # "would substitute a GRAPE-must nitrogen partition for a DEFINED SYNTHETIC medium whose
    # composition is in the paper and not in this repo ... closing the gap needs the paper."
    # D-246 SOURCED that medium -- Bely, Sablayrolles & Barre 1990, the same one Crepin and
    # Minebois both run -- and put both musts in this repo as
    # ``tests.test_defined_media._MUSTS`` / ``commensurate_pools``. So the paper is here, and a
    # comment telling the next reader it is not was standing instruction to not even try. This
    # is the 8th time in this archive that "blocked on external sourcing" outlived its blocker
    # (D-191/196/199/208/209/211, and D-230's variant where the source was already on disk).
    #
    # **The fixture still does not migrate, and the reason is now a MEASURED one.** D-246 section 7
    # keeps this must as a legitimate characterization object; what was illegitimate was reading
    # paper numbers off it. The commensurate scoring lives beside it in
    # ``tests/test_defined_media.py``, on the papers' own musts. D-254 measured what migrating
    # the guards themselves would do and found it fails for one alcohol whichever estimator is
    # used -- see ``tests/test_fusel_provenance_estimator.py``. Migration is therefore blocked
    # on an OPEN MODEL QUESTION, not on a missing paper. Do not restore the old reason.
    return Scenario(
        name="d109-fusel-node",
        medium="wine",
        initial=initial,
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=_CREPIN_TEMP),
            TemperaturePoint(day=_FERMENT_DAYS, celsius=18.0),
        ],
        interventions=interventions,
        duration_days=duration,
    )


def _run(
    *,
    drop: tuple[str, ...] = (),
    scale: dict[str, float] | None = None,
    set_params: dict[str, float] | None = None,
    **kw,
):
    cs = compile_scenario(_scenario(**kw))
    for name in drop:
        cs.process_set.disable(name)  # KeyErrors on a renamed Process rather than silently no-op
    # ``cs.param_values`` is a PROPERTY that re-resolves from ``parameters`` on every access, so
    # ``cs.param_values[k] = v`` mutates a throwaway dict and ``cs.run()`` re-resolves a clean one
    # — a knob that looks connected and is not. Resolve ONCE and hand that dict to the runtime,
    # exactly as ``CompiledScenario.run`` does. (This bit D-109's own first probe: it reported a
    # perfect "+0.000000%" agreement with the answer being tested for.)
    pv = cs.param_values
    for key, factor in (scale or {}).items():
        assert key in pv, f"no such parameter {key!r}"
        pv[key] = pv[key] * factor
    for key, value in (set_params or {}).items():
        assert key in pv, f"no such parameter {key!r}"
        pv[key] = value
    traj = simulate_scheduled(
        cs.process_set,
        pv,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
    )
    assert traj.success, traj.message
    return traj, cs.process_set.schema


def _end(traj, schema, name: str) -> float:
    return float(traj.y[schema.slice(name), -1][0])


# -- finding 1: the excreted pool cannot supply propanol (design A is INFEASIBLE) --------------


def test_the_excreted_pool_cannot_supply_propanol():
    """Propanol's molar demand exceeds every gram of α-KB the pool ever carries (decision D-109).

    The prescribed "re-base propanol on ``alpha_ketobutyrate``" is rejected on **principle** —
    D-49's flux-intermediate-vs-excreted-residual test, applied symmetrically to the pool D-107
    selected for sotolon on exactly that test. This is the **arithmetic** that corroborates it, and
    it is the sharper statement: propanol is 1 mol 2-KB decarboxylated, so drawing it from this
    pool needs 1 mol of pool per mol of alcohol — and the pool never holds three quarters of that
    (74 % since D-245's re-record; 38 % at the pre-D-244 evaluation point).

    So the re-base does not merely mis-attribute the competition. It **cannot be built**: it would
    starve propanol (breaking its independently-anchored D-99 magnitude) *and* collapse sotolon's
    substrate — a large, dramatic, wrong result that would look exactly like the long-promised
    competition finally expressing itself.

    **THE SAME RATIO ARGUES BOTH WAYS, AND THAT IS THE BEAT'S RESULT.** Read as "can this pool
    supply propanol?" it is fatal to the shortcut. Read as "how big is propanol next to the
    excretion flux?" it says propanol's 2-KB demand exceeds that flux — so an honestly partitioned
    intracellular node would couple propanol and sotolon materially. The promised payoff is real;
    it is the *location* that was wrong. (Not "the node's dominant sink" — 2-KB's committed route
    is isoleucine biosynthesis, which this model does not carry; see the module docstring.)

    **RE-RECORDED AT D-245, OLD → NEW, BECAUSE THE MARGIN MOVED AND THE CLAIM DID NOT.** D-244's
    corrected yield halves the biomass; propanol is produced on the fermentative flux and halves
    with it, while α-KB excretion barely moves (throughput 0.3181 → 0.3304 mmol/L across the whole
    biomass sweep), so the ratio falls **2.60 → 1.358**. The finding — *this pool cannot supply
    propanol* — survives untouched and is asserted on its own line. What does not survive is the
    ``> 2.0`` pin, which was a margin, not the claim; it is replaced by a two-sided record of the
    margin the model now actually holds. **The honest downgrade is worth stating: the supply
    argument now rests on 36 % headroom over an author-estimated excretion rate** (D-109's own
    "order-of-magnitude claim, not a calibration"), where it used to rest on 160 %.

    TRIPWIRE, not a curiosity: if a future beat raises the excretion rate enough for the pool to
    supply propanol, the first assertion fails and the design question genuinely re-opens. The
    second fires on any move in either direction, which is the point of pinning it two-sided —
    a one-sided floor cannot catch the excretion flux *growing* toward the demand.
    """
    base, schema = _run(aging=False)
    # Reassimilation OFF ⇒ the pool accumulates EVERY gram excreted and nothing removes it ⇒ the
    # end state IS the total throughput, exactly. No rate integrated over interpolated states
    # (D-103's finding, where a nonlinear rate over linear interpolation overstated a draw 1.3–3.5×
    # and booked more carbon out of valine than valine ever held).
    no_reassim, _ = _run(aging=False, drop=("alpha_kb_reassimilation",))

    residual = _end(base, schema, "alpha_ketobutyrate")
    throughput = _end(no_reassim, schema, "alpha_ketobutyrate")
    propanol = _end(base, schema, "propanol")

    # Anti-vacuity: all three must be real, or the ratio below is a division of zeros that
    # "agrees" for reasons having nothing to do with the finding (D-106/D-108).
    assert residual > 0.0 and throughput > residual and propanol > 0.0

    demand_mmol = propanol / M_PROPANOL * 1e3
    supply_mmol = throughput / M_ALPHA_KETOBUTYRATE * 1e3
    ratio = demand_mmol / supply_mmol
    # THE CLAIM. Everything D-109 argues from this number needs only that it exceeds 1.
    assert ratio > 1.0, (
        f"propanol demand {demand_mmol:.4f} mmol/L vs total α-KB ever excreted "
        f"{supply_mmol:.4f} mmol/L — if the pool can now supply propanol, D-109's supply "
        "argument against re-basing it no longer holds and the design must be re-measured"
    )
    # THE MARGIN, re-recorded at D-245 (was > 2.0, measured 2.60, at the pre-D-244 evaluation
    # point) and AGAIN at D-248. Two-sided: a one-sided floor cannot catch the excretion flux
    # growing toward the demand, which is precisely the direction that would re-open the design.
    #
    # D-248 moves it BACK UP, 1.358 -> 1.599, and the direction is the good one. Propanol is made
    # on the fermentative flux, which scales with biomass; un-coupling nitrogen uptake restores
    # the biomass D-244's corrected yield had halved, while α-KB excretion barely moves (it is a
    # per-flux yield on a pool with no biomass term of its own). So the headroom this supply
    # argument rests on goes from 36 % back to 60 %. D-245's "honest downgrade" note above is
    # therefore partly retired — but only partly, and it stays on the page: 60 % over an
    # author-estimated excretion rate is still not the 160 % the claim originally had.
    assert 1.52 < ratio < 1.68, (
        f"propanol/α-KB demand ratio {ratio:.4f} left the D-248 band [1.52, 1.68] (measured "
        f"1.599; D-245 measured 1.358 before nitrogen uptake was un-coupled). The CLAIM above "
        "still holds — this is the margin it holds it by. Re-derive it rather than widening"
    )


# -- finding 2: the de-novo supply structure the node must preserve (SOURCED for four of five) --


def ehrlich_primary_share(params, precursor: str, *, f_override: float | None = None) -> float:
    """The fraction of CONSUMED ``precursor`` the model routes to its **primary** alcohol (D-245).

    ``1 − f_non_ehrlich_<precursor> − Σ secondary shares``. For every precursor but valine the sum
    is empty and this is the familiar ``(1 − f)``; **valine is the exception, and it is the whole
    reason this helper exists.** Since D-111 valine's consumption splits three ways — 0.15 to
    isobutanol (the residue this returns), 0.23 to isoamyl alcohol via α-ketoisocaproate, 0.62 to
    the non-Ehrlich lump — and :func:`~fermentation.core.kinetics.byproducts.ehrlich_draws` sizes
    the branches on exactly that residue (``consumed = primary.precursor_carbon / share_primary``).

    **This is a repair, and the defect it repairs was invisible for the usual reason (D-245).**
    :func:`_de_novo_share` used ``(1 − f) = 0.38`` for every precursor, which charges isobutanol
    with the valine carbon that actually becomes *isoamyl alcohol* — an over-count of **2.533×**,
    live from D-111 until D-245. It never failed a guard because isobutanol cleared the 80 % floor
    even mis-measured, until D-244's corrected yield pushed the mis-measured number under it and
    the miss read as a model defect. The sibling ``_amino_acid_share`` had the residue right and
    said so in a comment; **the two callers of "the same" quantity disagreed and only one was
    read** — the D-106 shared-helper failure exactly, one file over. Hence one helper, four
    callers, and it lives here because this is the module the sibling already imports from.

    **The residue is exact only outside the ``ehrlich_draws`` headroom-cap window, and the window
    is measured rather than assumed.** That cap refuses to source the secondary branch from more
    precursor than its alcohol is being made from; at t = 0 leucine's own gate is ≈1, leaving
    isoamyl no headroom, so the cap binds and valine's realised primary share starts at 0.2642
    instead of 0.15. Measured on a dense grid (0.024 h) it releases at t ≈ 1.35 h with **at most
    12.8 % of the valine pool consumed**, and the share falls monotonically to 0.15 across the
    window. So an isobutanol number computed from this residue is a **bound, not a point**: the
    worst-case blend 0.128 × 0.2642 + 0.872 × 0.15 = 0.1646 puts the de-novo share at ≥ 89.6 %
    against the 90.5 % this returns. Everything asserted against it clears both. Do not re-state
    the 90.5 % as exact, and do not "fix" the gap by integrating the branch — that is the
    quadrature D-103 forbids.
    """
    f = params[non_ehrlich_fraction_param(precursor)] if f_override is None else f_override
    secondary = sum(
        params[route.share_param]
        for route in SECONDARY_FUSEL_ROUTES
        if route.precursor == precursor
    )
    return float(1.0 - f - secondary)


def _de_novo_share(spec, *, f_override: float | None = None) -> float:
    """This alcohol's de-novo carbon share, from EXACT state differences (decision D-109).

    The D-104 split invariant (consumed precursor splits exactly ``f : (1−f)`` between the
    non-Ehrlich lump and the alcohol) sizes the Ehrlich draw with no quadrature — D-103, where a
    nonlinear rate integrated over linearly-interpolated states overstated a draw 1.3–3.5×. The
    other precursor consumers are disabled so that invariant is *exact* rather than nearly-true —
    verified rather than assumed at D-245: only ``fusel_amino_acid_reroute`` and
    ``precursor_non_ehrlich_fates`` touch a speciated precursor pool on this run.
    ``n_alc/(n_alc+1)`` removes the D-106 decarboxylation CO₂: the draw is a full mole of precursor
    per mole of alcohol, of which one carbon leaves as CO₂ instead of reaching the alcohol.

    **The Ehrlich side of that split is not one branch since D-111**, so the share going to *this*
    alcohol is :func:`ehrlich_primary_share`, not ``(1 − f)`` — see that helper for the 2.533×
    isobutanol over-count this line carried from D-111 to D-245, and for why its result is a bound
    while the ``ehrlich_draws`` headroom cap binds.

    **One helper, two callers** — the sourced-floor test and the 2-PE exclusion must compute this
    the *same* way, or the exclusion could be argued from arithmetic the floor never used (the
    D-33/D-99/D-106 shared-helper discipline; D-106 is the beat where two callers "recomputing
    exactly the same thing" agreed **by luck** until one of them changed).

    **This helper counts ``spec.precursor_amino_acid`` ONLY, and for isoamyl alcohol that omits a
    real branch** — valine → KIC, which ``_amino_acid_share`` does count. The omission is
    deliberate (it is the *primary*-route share) and is stated here rather than only in the
    sibling, because a floor asserted against a primary-only number is the same class of thing
    D-245 just repaired one alcohol over. Named so the next beat does not rediscover it as new.
    """
    param = non_ehrlich_fraction_param(spec.precursor_amino_acid)
    overrides = {} if f_override is None else {param: f_override}
    traj, schema = _run(aging=False, drop=_OTHER_PRECURSOR_CONSUMERS, set_params=overrides)
    params = compile_scenario(_scenario(aging=False)).param_values
    return de_novo_share_of(traj, schema, params, spec, f_override=f_override)


def de_novo_share_of(traj, schema, params, spec, *, f_override: float | None = None) -> float:
    """:func:`_de_novo_share`'s arithmetic, on a run the CALLER supplies (decision D-246).

    Split out so a probe scoring this share on a *different must* runs the identical arithmetic
    rather than a fourth copy of it — the D-106/D-245 shared-helper discipline, applied before
    the copy exists this time instead of after it has silently diverged for 134 records.
    :mod:`tests.test_defined_media` scores Crepin's and Minebois's own defined media through
    this helper, and the anchor that licenses the comparison is that doing so on the D-109
    fixture reproduces the shipped numbers to 4 dp.
    """
    precursor = spec.precursor_amino_acid
    consumed = float(traj.y[schema.slice(precursor), 0][0]) - _end(traj, schema, precursor)
    made = _end(traj, schema, spec.pool)
    assert consumed > 0.0 and made > 0.0, "vacuous: nothing consumed or nothing made"

    share = ehrlich_primary_share(params, precursor, f_override=f_override)
    n_alc = CARBON_ATOMS[spec.species]
    draw_carbon = share * consumed * carbon_mass_fraction(precursor)
    alcohol_carbon_from_precursor = draw_carbon * n_alc / (n_alc + 1.0)
    total_alcohol_carbon = made * carbon_mass_fraction(spec.species)
    return 1.0 - alcohol_carbon_from_precursor / total_alcohol_carbon


@pytest.mark.parametrize(
    "spec",
    [
        pytest.param(
            s,
            id=s.pool,
            # ISOBUTANOL IS NOT HERE, AND ITS REMOVAL IS NOT A WEAKENING (D-245). D-244 marked it
            # alongside propanol at a measured 76.0 %; that number was the harness over-charging
            # it 2.533× for valine carbon that becomes isoamyl alcohol (see
            # :func:`ehrlich_primary_share`). Measured on the residue the model actually routes,
            # isobutanol is 90.5 % de novo — clear of the floor by ten points, and ≥ 89.6 % even
            # on the worst-case cap-window bound. Propanol's miss is the real one: threonine has
            # no second branch, so its 77.4 % was never touched by the defect.
            #
            # AND PROPANOL'S XFAIL IS GONE AT D-248 — the parametrization now carries no marks at
            # all. It read 0.7744 here from D-245 to D-247 and reads **0.8062**; the cause was the
            # biomass denominator, not this alcohol's chemistry (see
            # ``_D245_DE_NOVO_FLOOR_GAP_CLOSED_AT_D248`` above). The floor itself is untouched:
            # ``_SOURCED_DE_NOVO_FLOOR`` is still 0.80 and was never lowered to meet the model,
            # which is what makes the close worth having.
        )
        for s in _SOURCED_FUSEL_SPECS
    ],
)
def test_every_sourced_fusel_is_de_novo_dominated(spec):
    """Every *sourced* Ehrlich alcohol draws ≥80% of its carbon de novo (decision D-109).

    The **sourced supply structure the fusel-side node has to preserve**, pinned here because the
    milestone will rewrite the sourcing layer that produces it: Crépin measures 2-ketobutyrate as
    81% newly synthesised, and Rollero independently puts >90% of the higher alcohols' carbon on
    the central carbon metabolism. This is the model reproducing that, and it is what makes the
    intracellular node *mostly de novo* — the property D-104 identified as the reason reality
    escapes its inverted split (isoamyl is built from KIC, so the leucine pool never faces the
    isoamyl demand).

    **It is NOT an argument that the competition is small** — D-109's first draft used it that way.
    Both arms of a keto-acid node would draw the *same de-novo 2-KB*, so de-novo dominance says
    where the carbon comes from, not whether the two consumers compete for it. The competition over
    threonine is zero for a **structural** reason instead (see
    :func:`test_alpha_kb_production_is_exactly_threonine_independent`).

    **2-phenylethanol is deliberately NOT in this parametrization** — see
    :func:`test_2_phenylethanol_carries_no_sourced_de_novo_floor`.

    Measured from EXACT state differences via the D-104 split invariant (consumed precursor splits
    exactly ``f : (1−f)`` between the non-Ehrlich lump and the alcohol), with the other precursor
    consumers disabled so that invariant is exact. The ``n_alc/(n_alc+1)`` factor removes the D-106
    decarboxylation CO₂: the Ehrlich draw is a full mole of precursor per mole of alcohol, of which
    one carbon leaves as CO₂ rather than reaching the alcohol.
    """
    de_novo_share = _de_novo_share(spec)
    assert de_novo_share >= _SOURCED_DE_NOVO_FLOOR, (
        f"{spec.pool} is only {de_novo_share:.1%} de novo, under the sourced floor "
        f"{_SOURCED_DE_NOVO_FLOOR:.0%} (Crépin 81% newly-synthesised 2-KB; Rollero >90% CCM). "
        "The fusel-side node is specified to PRESERVE this supply structure — if a beat moved it, "
        "the milestone's premise (the intracellular keto acid is mostly de novo) moved with it"
    )


def test_2_phenylethanol_carries_no_sourced_de_novo_floor():
    """Why 2-PE is *still* excluded from the floor above — on ONE leg now, not two (D-117).

    **This test was D-109's tripwire and it fired.** It pinned ``tier: speculative`` /
    ``source: "author estimate"`` so that *"a future beat that sources phenylalanine moves it
    deliberately"*. D-117 is that beat: **Minebois et al. 2025 labels U-13C phenylalanine** in
    synthetic must with *S. cerevisiae*, and ``f_non_ehrlich_phenylalanine`` is now **0.975,
    plausible**, off a stated consumed-precursor denominator. Both of the old exclusion's legs are
    re-examined below; **one died, one stands**, and keeping the exclusion for the surviving reason
    is the point of the test.

    **Leg 1 — "the precursor is unsourced" — IS DEAD.** It is sourced, per-compound, in the right
    organism and matrix, at the same growth stage Crépin quotes. This leg can never be cited again.

    **Leg 2 — "the floor is inside the parameter's own band" — IS ALSO DEAD, and measured so.**
    The old band's low end (0.38) gave **78.6%** de novo, *under* the 80% floor. D-117's band is
    **zero-width at the shipped bound 0.531**, where the de-novo share is **83.6%** — clear, and
    unsamplable besides. (It is zero-width on purpose: the honest interval reaches the measured
    lump 0.975, which breaks carbon conservation, so it is held in a test constant rather than in
    a field the ensemble draws from — see ``test_the_sourced_lump_breaks_the_carbon_refund_guard``.
    Across that *real* interval the de-novo share runs 83.6% → 99.1%, clear at both ends, so leg 2
    stays dead however the band is later restored.)

    **What survives is narrower and is purely a SCOPE claim: the floor's own sources do not
    describe this keto acid.** The floor is justified as *"Crépin 81% newly-synthesised 2-KB;
    Rollero >90% CCM"* — measurements of **2-ketobutyrate** and of bulk central-carbon contribution.
    Neither lab measured **phenylpyruvate**. Sourcing ``f_phe`` does not retroactively widen what
    Crépin's 2-KB number describes, so admitting 2-PE to a floor carrying that citation would still
    be **D-104's rule broken — a cited number binding a set it does not describe** — even though the
    model would now pass. **Passing is not the same as being covered, and this is the beat where
    that distinction was easiest to spend.**

    **The next beat is therefore sharp and small:** Minebois's Figure 6A splits each volatile into
    labelled/unlabelled (*"the unlabelled fraction ... represents the ... volatile compounds de novo
    synthesised from CCM precursors"*), which is **2-PE's de-novo share measured directly**. Extract
    that number and 2-PE joins ``_SOURCED_FUSEL_SPECS`` under its *own* citation — with 97.9% of
    headroom already measured. It was not extracted here rather than guessed.
    """
    cs = compile_scenario(_scenario(aging=False))
    param = non_ehrlich_fraction_param(_FLOOR_EXCLUDED_PRECURSOR)
    entry = cs.parameters[param]
    spec = next(s for s in FUSEL_SPECS if s.precursor_amino_acid == _FLOOR_EXCLUDED_PRECURSOR)

    # Leg 1 is dead: pin that it is sourced, so nobody re-argues the exclusion from "unsourced".
    assert entry.tier is Tier.PLAUSIBLE, f"{param} changed tier — re-derive BOTH legs, not just it"
    assert entry.provenance.source != "author estimate"
    assert "Minebois" in entry.provenance.source

    # Leg 2 is dead: the floor is OUTSIDE the band now, asserted with RUNS at both ends rather than
    # stated in prose. D-109's own lesson — the first draft of this test asserted
    # `low < value < high` under a comment claiming the floor was inside the band, which is
    # trivially true of every parameter and says nothing. The claim needs the runs.
    assert entry.uncertainty is not None
    for edge in (entry.uncertainty.low, entry.uncertainty.high):
        share = _de_novo_share(spec, f_override=edge)
        assert share > _SOURCED_DE_NOVO_FLOOR, (
            f"2-PE is {share:.1%} de novo at f_phe={edge} — back under the "
            f"{_SOURCED_DE_NOVO_FLOOR:.0%} floor, so leg 2 of the exclusion has revived "
            "and the docstring above is stale"
        )

    # ...and yet it stays excluded, on the SCOPE leg alone. If a beat sources phenylpyruvate's own
    # de-novo share (Minebois Fig. 6A), this is the assertion to move — deliberately, as D-109
    # intended and D-117 honoured.
    assert _FLOOR_EXCLUDED_PRECURSOR not in {s.precursor_amino_acid for s in _SOURCED_FUSEL_SPECS}


def test_alpha_kb_production_is_exactly_threonine_independent():
    """The competition over threonine is **structurally zero**, so "GONE" was right (D-109).

    D-109 set out to correct the work-list's "the propanol-vs-sotolon competition is GONE rather
    than wrong" and **the work-list was right**. :class:`AlphaKetobutyrateExcretion`'s rate is
    ``k · X · S/(K+S)`` — *flux-only*. Threonine's depletion gate re-routes the **carbon source**
    (threonine vs the sugar de-novo stand-in) and does not touch the **rate**. So draining
    threonine does not cost the α-KB pool one microgram, and propanol — which drains threonine —
    cannot reduce sotolon's substrate through it. Measured here at the derivative level, where it
    is exact rather than nearly-true.

    **D-107 chose this deliberately and the docstring says why**: gating the rate on threonine
    would empty the pool in a threonine-free wine and kill sotolon — the D-104 canary, one pool
    upstream. So this is not an oversight to fix; it is the reason the competition must be built
    at the **keto acid**, not at its grandparent.

    **WHY THIS TEST EXISTS, AND IT IS THE BEAT'S OWN NEAR-MISS.** D-109's first probe drove
    ``k_propanol`` from 0 to 10× (propanol ~0 → ~500 mg/L) and measured sotolon moving 0.42%,
    monotone and correctly signed. That was written up as "the competition is present, correctly
    signed, and negligible" — a *mechanism* that was never measured, attached to an assertion that
    measures something else. The channel is the **sugar ledger** (propanol's de-novo carbon comes
    out of ``S``, which every flux-linked rate reads), not threonine. The mutation test is what
    caught it: deleting the threonine draw entirely left that assertion **passing**. D-96/D-102/
    D-108's *"the sentence and the assertion are not the same claim"*, a fifth time — and this time
    it would have shipped as a correction to a claim that was **true**.
    """
    from fermentation.core.kinetics.keto_acids import AlphaKetobutyrateExcretion

    cs = compile_scenario(_scenario(aging=False))
    pv = cs.param_values
    schema = cs.process_set.schema

    def _derivative(threonine_gpl: float):
        y = schema.zeros()
        y[schema.slice("X")] = 2.0
        y[schema.slice("S")] = 100.0
        y[schema.slice("N")] = 0.1
        y[schema.slice("T")] = pv["T_ref"]
        y[schema.slice("threonine")] = threonine_gpl
        return AlphaKetobutyrateExcretion().derivatives(0.0, y, schema, pv)

    replete = _derivative(0.067)  # a real must's threonine
    starved = _derivative(0.0)  # none at all

    rate_replete = float(replete[schema.slice("alpha_ketobutyrate")][0])
    rate_starved = float(starved[schema.slice("alpha_ketobutyrate")][0])
    assert rate_replete > 0.0, "vacuous: the Process is not firing in this state at all"
    # EXACT equality, not approx: the rate does not read threonine, so this is not a small
    # dependence to bound — it is no dependence at all, and == is the only honest assertion.
    assert rate_starved == rate_replete

    # ...and the gate DID do its job — it moved the carbon source, which is the half that is real.
    # Without this the test above would also pass on a Process that had been silently disabled.
    assert float(replete[schema.slice("threonine")][0]) < 0.0  # drawn from threonine
    assert float(starved[schema.slice("threonine")][0]) == 0.0  # none to draw
    assert float(starved[schema.slice("S")][0]) < float(replete[schema.slice("S")][0])  # → sugar
