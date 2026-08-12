"""The oxidative cascade — every O2 sink re-homed behind one Fe(II)+O2 activation node (D-141).

This module is the **isolable alternative** to the direct oxidative sinks in :mod:`aging.py`.
The two sets are **mutually exclusive**: exactly one is wired into a given build (decision
D-139), and toggling the cascade off must *restore* the direct sinks rather than delete the
oxidative axis. That is a different shape from every other ``_*_PROCESSES`` tuple in
:mod:`~fermentation.core.media`, all of which are **additive** (off by default, switched on by
``begin_aging``, with an un-aged run byte-for-byte the pre-aging core).

**Why the rebuild exists.** Since D-71 this model has carried N independent sinks drawing on one
shared ``o2`` pool, splitting it by ``k_i / sum(k)`` — a structure that reads as competition and
sums correctly, but that asserts something Gate 1 (D-137) found to be false: that ethanol,
bisulfite, phenolics, amino acids, anthocyanin and oak tannin each react with **dissolved O2**.
In wine they do not. O2 reacts with essentially nothing except **Fe(II)**, and the two oxidants
that step produces — H2O2 and o-quinone — are what the six "sinks" actually consume::

    O2  + 2 Fe(II) + 2 H+  ->  H2O2 + 2 Fe(III)        the rate-determining ACTIVATION step
    2 Fe(III) + o-diphenol ->  o-quinone + 2 Fe(II)    catalyst turnover; makes the 2nd oxidant

So **one mol of O2 yields two oxidising equivalents**: one H2O2 and one quinone. That single
fact is the rebuild's whole falsifiable content, because it is what produces Danilewicz 2016's
measured sulfite series without any of it being fitted:

**One O2 more leaves the pool without yielding any oxidising equivalent at all** (D-196): the
1-hydroxyethyl radical the Fenton limb makes reacts with O2 to give acetaldehyde, so that limb
costs two O2 per acetaldehyde and only the first of them is an activation. See
:data:`_O2_PER_ACETALDEHYDE`. The series below is unaffected at its limits and only there,
because the draw is weighted by the ethanol branch's share of H2O2, which goes to zero exactly
where the limits are asserted.

* **1:1** — quinone route blocked (benzenesulfinic acid traps the quinone): only the H2O2 node
  can oxidise bisulfite, so one SO2 per O2.
* **1:2** — both routes open and bisulfite wins both: two SO2 per O2.
* **1:1.7** — real wine, where bisulfite wins the H2O2 node but only part of the quinone,
  the rest going to the competing nucleophiles and to polymerisation.

**1.7 must EMERGE from partial quinone capture. Do not fit it** (D-138 says this twice, and
this docstring says it a third time because this module is where it would be tempting).

**Sizing (D-138).** Only ``quinone`` becomes a state slot. ``H2O2`` does not — at a 17 ms
half-life against a five-year integration the differential and algebraic formulations agree to
~10 decimal places and only the algebraic one is non-stiff, so quasi-steady-state is the
numerically *correct* treatment here rather than a simplification of one. ``Fe(III)`` does not
either, on separate grounds: Nguyen's Table 3.1 measures it reduced 3.3-63x faster than O2 is
consumed (~18x at representative wine pH and copper), so it too is quasi-steady-state — with a
margin this module records as thin rather than comfortable.

**Consequence of QSS for the architecture.** ``ProcessSet`` sums independent Processes and no
Process sees another, so the H2O2 branch fractions cannot be computed by one Process and shared.
They do not need to be: the branch denominator is a **pure function of state and params**, so
each consumer computes it independently from ``y``. That is the existing shared-gate idiom —
:func:`~fermentation.core.kinetics.amino_acid_pools.depletion_gate` (D-100) and the malolactic
environmental gate (D-31/D-131) are both pure ``(y, schema, params, ...) -> float`` helpers
called independently by several Processes. :func:`activation_rate` and
:func:`h2o2_branch_fraction` here are the same shape.

The **quinone** branch needs no such helper at all, which is the quiet win of making it a slot:
its consumers each draw first-order in ``[quinone]`` from the shared pool, and the competition
between them emerges from ``ProcessSet`` summing exactly as the old ``o2`` competition did.

**What is re-homed vs what is new.** Every consumer keeps **its own existing rate constant**,
re-read against the oxidant that actually oxidises it instead of against ``o2``. No consumer
coefficient is invented. Genuinely new: the activation node's own parameters (three of which are
the pre-cascade constants re-expressed as one un-partitioned total — see ``aging.yaml``) and
:data:`k_quinone_polymerization`, the always-available fate that keeps quinone bounded and makes
``A420`` a quinone *fate* rather than an O2 *yield*.

**That economy does not make the re-home magnitude-neutral, and D-141 measured that it is not.**
Keeping ``k`` while swapping ``[o2]`` for ``[quinone]`` rescales every bilinear sink by
``[quinone]_ss / [o2]_ss`` — ~0.06 in an unsulfited red, and scenario-dependent because the two
pools are set by different balances. Measured against the direct set on the same scenario, the O2
*budget* reproduces (1.01x unsulfited, 0.88x at 60 mg/L SO2) while the *fates* move by up to 25x:
A420 1.4x/4.2x high, the Strecker aldehydes 20-25x low. **This module is therefore built, wired
and tested but NOT the default wiring** — see :data:`~fermentation.core.media._OXIDATIVE_SETS`
for why closing that gap is a sourcing problem (Nikolantonaki & Waterhouse 2012) and not a
tuning one.

**Medium-agnostic (fork D3), with a floor.** Beer carries no copper, no iron and none of the
reductant slots, so every read here is guarded and beer's activation sits at
``k_activation_floor`` — the same lumped anchor its O2 uptake already rode on. D-139 measured
what a wine-only activation node would have cost: beer consumes 71% of its packaging oxygen
today, and that would silently have become 0.00 with no existing test noticing. It is guarded
now (``tests/test_oxidative_cascade_guards.py``).

**Tier: speculative (fork D4)**, uniformly. Every sink this replaces is already speculative, so
``Tier.combine`` moves nothing downstream — verified, not assumed.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.acidbase import (
    PH_SYSTEM_READS,
    SO2_BINDING_READS,
    SO2_STATE_KEY,
    bisulfite_so2_at_ph,
    ph_of_state,
)
from fermentation.core.chemistry import (
    M_ACETALDEHYDE,
    M_CO2,
    M_ETHANOL,
    M_METHIONAL,
    M_O2,
    M_PHENYLACETALDEHYDE,
    M_QUINONE,
    M_SO2,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.amino_acid_pools import (
    SPEC_BY_SPECIES,
    depletion_gate,
    draw_precursor_carbon,
)
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.o2_partition import o2_depletion_shares
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The three reductant pools whose presence accelerates Fe(III) -> Fe(II) turnover, and hence the
#: activation rate itself. Read through ``k_activation_phenolic``; each is guarded, so a schema
#: without one contributes zero. ``hydroxycinnamics`` is the term the pre-cascade browning boost
#: did NOT have, and it is required: D-138's falsifier says that if activation reads only ``o2``
#: and a lumped ``k``, effective O2 uptake stays FLAT across five years where Nguyen predicts it
#: should DECLINE as the wine's reducing capacity is exhausted.
_PHENOLIC_REDUCTANT_POOLS = ("tannin", "anthocyanin", "hydroxycinnamics")

#: Moles of H2O2 produced per mole of O2 activated. Reaction stoichiometry (Danilewicz's
#: oxygen-reduction scheme: ``O2 + 2 Fe(II) + 2 H+ -> H2O2 + 2 Fe(III)``), NOT an uncertain YAML
#: parameter — the ``_SO2_PER_O2`` precedent in :mod:`aging.py`.
_H2O2_PER_O2 = 1.0

#: Moles of o-quinone produced per mole of O2 activated. The 2 Fe(III) that step makes are
#: returned to Fe(II) by one o-diphenol, which becomes one o-quinone
#: (``2 Fe(III) + o-diphenol -> o-quinone + 2 Fe(II)``). Together with :data:`_H2O2_PER_O2` this
#: is the **two oxidising equivalents per O2** that make Danilewicz's 1:2 limit emerge rather
#: than being asserted.
_QUINONE_PER_O2 = 1.0

#: Moles of bisulfite oxidised per mole of H2O2 (``HSO3- + H2O2 -> HSO4- + H2O``). Stoichiometry.
_SO2_PER_H2O2 = 1.0

#: Moles of bisulfite consumed per mole of quinone ~~sulfonated (one nucleophilic addition of
#: bisulfite to the quinone ring)~~ **consumed by bisulfite, by EITHER route (D-199)**.
#: Stoichiometry, and it is 1 both ways: the sulfonate route adds one bisulfite to the ring, and
#: the reduction route — which *Understanding Wine Chemistry* §24.4.3.2 gives as >90 % of the
#: reaction at wine pH — oxidises one bisulfite to sulfate while regenerating the o-diphenol.
#: So the product identity D-199 corrected on :class:`QuinoneSulfonation` moves nothing here.
_SO2_PER_QUINONE = 1.0

#: Moles of acetaldehyde per mole of H2O2 taking the ethanol route (Fenton: ``H2O2 + Fe(II) ->
#: HO. + OH-``; ``HO. + ethanol -> 1-hydroxyethyl radical -> acetaldehyde``). Stoichiometry, and
#: a deliberate change of KIND from the pre-cascade ``y_acetaldehyde_per_o2``, which was a
#: *fitted yield* per O2 because the pre-cascade structure had no mechanism to derive one from.
_ACETALDEHYDE_PER_H2O2 = 1.0

#: Ethanol debited per acetaldehyde formed — the carbon-exact C2 borrow (the D-27 reduction
#: reversed), identical to :mod:`aging.py`'s ``_ETHANOL_PER_ACETALDEHYDE``.
_ETHANOL_PER_ACETALDEHYDE = M_ETHANOL / M_ACETALDEHYDE

#: Moles of O2 consumed by the 1-hydroxyethyl radical per mole of acetaldehyde it yields —
#: the **second** O2 the Fenton limb costs, which this module did not book until D-196.
#:
#: Mechanism, printed: "reaction with the next highest concentration component of wine, ethanol,
#: generates the major detectable product of the Fenton reaction, the 1-hydroxyethyl radical.
#: **This radical goes on to react with oxygen and produce acetaldehyde** unless it can be
#: quenched by other processes, such as reaction with hydroxycinnamates" (Waterhouse, Sacks &
#: Jeffery, *Understanding Wine Chemistry* 2nd ed., 24.4.4.1). So one acetaldehyde made this way
#: costs TWO O2: one activated at the iron node, one taken by the radical.
#:
#: **Why it is pinned per ACETALDEHYDE and not per H2O2.** The source names a competitor for the
#: radical (hydroxycinnamate quenching) and gives no rate for it. That competition needs no
#: parameter here, because quenching removes the acetaldehyde and the second O2 *together* — a
#: quenched radical yields neither. Pinning the draw to the acetaldehyde this Process already
#: forms therefore absorbs the unquantified quench into an already-calibrated rate; pinning it to
#: the H2O2 flux would not, and would need the quench ratio nobody has published.
#:
#: ~~**This is an UPPER BOUND on the NET draw, deliberately.**~~ **STRUCK (D-198): it is a LOWER
#: bound.** The radical's other product is the hydroperoxyl radical, which in principle returns to
#: H2O2 via Fe(II) and hands an oxidising equivalent back. D-196 read that as a draw this constant
#: over-charges. Measured, it is the opposite: H2O2 is quasi-steady-state here, so a limb that
#: RETURNS H2O2 feeds its own production flux and closes as a geometric series, ``F = A/(1 - s)``
#: with ``s`` the ethanol branch share. That draws ``s^2/(1 - s)`` MORE O2 than this constant at
#: every ``s > 0`` — a RATE-LAW change that unseats the activation node as rate-determining step,
#: not the stoichiometric over-charge "upper bound" implies.
#:
#: D-196's OTHER claim here — that the limb pushes the SO2:O2 ratio the other way — is
#: **CONFIRMED** (the gap is ``1.5 s (1 - s)/(1 + s) >= 0``, zero only at the endpoints), and
#: D-198's own pre-registration predicted it false. So is this one: nothing sources **that limb's**
#: fate in wine (Danilewicz 2007's radical chain is *bisulfite* autoxidation, which he concludes
#: does not run in wine because polyphenols intercept it). What D-198 added is that the species has
#: THREE sourced fates in other contexts — HO2. -> H2O2, dismutation, Haber-Weiss — which disagree,
#: and that the series has a pole at ``s = 1``, which unsulfited wine reaches. Still named and NOT
#: built, now on four measured grounds rather than an inherited label. Recorded so a later
#: recycling term cannot be credited with a cancellation.
_O2_PER_ACETALDEHYDE = 1.0

#: Moles of CO2 released per Strecker aldehyde formed (the amino acid's own carboxyl carbon).
_CO2_PER_STRECKER_ALDEHYDE = 1.0

_METHIONAL_SPECIES = "methional"
_PHENYLACETALDEHYDE_SPECIES = "phenylacetaldehyde"
#: ``(product pool, precursor amino acid)`` — the D-100 speciated Strecker routes, each gated on
#: its OWN precursor so methional stops when methionine runs out and phenylacetaldehyde when
#: phenylalanine does, independently.
_STRECKER_ROUTES = (
    (_METHIONAL_SPECIES, "methionine"),
    (_PHENYLACETALDEHYDE_SPECIES, "phenylalanine"),
)


def _pool(y: FloatArray, schema: StateSchema, name: str) -> float:
    """A non-negative read of ``name``, or 0.0 if the schema has no such slot.

    The guarded-missing-slot idiom :class:`~fermentation.core.kinetics.aging.PhenolicBrowning`
    already uses for ``copper`` (D-134): a slot the medium does not carry reads as zero load
    rather than raising, which is what makes the whole cascade medium-agnostic. Clamping at 0
    also absorbs solver undershoot, so a momentarily negative pool cannot flip a rate's sign.
    """
    if name not in schema:
        return 0.0
    return max(float(y[schema.slice(name)][0]), 0.0)


def free_bisulfite(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Free bisulfite (HSO3-, g/L expressed as SO2), or 0.0 where SO2 is not tracked.

    The reactive antioxidant species is free **bisulfite**, not molecular SO2 (that is the
    antimicrobial form) and not total SO2 (bound bisulfite is already spent) — D-72's finding,
    carried over unchanged. Costs a pH solve, so every caller checks the cheap ``so2_total <= 0``
    guard first: an unsulfited or beer run must not pay a per-RHS pH solve for a zero result.
    """
    if SO2_STATE_KEY not in schema:
        return 0.0
    so2_total = float(y[schema.slice(SO2_STATE_KEY)][0])
    if so2_total <= 0.0:
        return 0.0
    return max(bisulfite_so2_at_ph(y, schema, params, ph_of_state(y, schema, params)), 0.0)


def activation_rate(
    y: FloatArray,
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    bisulfite: float | None = None,
) -> float:
    """The Fe(II)+O2 activation rate, g O2/L/h — the cascade's ONE draw on ``o2``.

    ::

        r_O2 = (k_activation_floor
                + k_activation_phenolic * ([tannin] + [anthocyanin] + [hydroxycinnamics])
                + k_activation_bisulfite * [HSO3-]) * f_Cu * f(T) * [o2]

    A pure function of state and params, so every downstream consumer can recompute it
    independently — the shared-gate idiom (D-100/D-31), and what lets the cascade work inside
    ``ProcessSet``'s sum-of-independent-Processes contract.

    **The reductant terms are not decoration.** They are the catalyst-turnover rate: Fe(III) must
    be reduced back to Fe(II) before the next O2 can be activated, and the pools that do the
    reducing are the o-diphenols and bisulfite. That is what makes Nguyen's predicted decline in
    O2 uptake across a long age *emergent* — as the reducing pools are exhausted, activation
    slows on its own. D-138 recorded the falsifier explicitly: a lumped constant here leaves
    uptake flat for five years, which would be wrong.

    ``f_Cu`` is D-134's mean-centered copper multiplier, re-homed from ``PhenolicBrowning``.
    Danilewicz's own measurement is what justifies the move: copper has no discernible effect on
    Fe(III) reduction at any pH, while O2 consumption rises with it — so copper accelerates the
    O2-activation step and nothing else. **This changes what the constant multiplies** (browning
    only, ~60% of the always-on rate, becomes the whole total including bisulfite), so it is a
    re-fit rather than a validation — a no-op at ``copper_typical``, +9.5% at 0.5 mg/L. Recorded
    at D-141, and it is the third time D-138's "a constant fitted to one structure does not
    survive being moved to another" has landed.

    Pass ``bisulfite`` when the caller has already paid for the pH solve, so a single RHS
    evaluation does not solve pH once per cascade Process.
    """
    o2 = _pool(y, schema, "o2")
    if o2 <= 0.0:
        return 0.0
    reductants = sum(_pool(y, schema, name) for name in _PHENOLIC_REDUCTANT_POOLS)
    hso3 = free_bisulfite(y, schema, params) if bisulfite is None else bisulfite
    k_eff = (
        params["k_o2_depletion_total"]
        + params["k_activation_phenolic"] * reductants
        + params["k_activation_bisulfite"] * hso3
    )
    # D-134's mean-centered copper multiplier, guarded: exactly 1.0 where copper is untracked
    # (beer) and exactly 1.0 at copper == copper_typical, so an un-overridden wine is unchanged.
    if "copper" in schema:
        k_eff *= 1.0 + params["k_copper_multiplier"] * (
            _pool(y, schema, "copper") - params["copper_typical"]
        )
    f_t = arrhenius_factor(
        float(y[schema.slice("T")][0]), params["E_a_activation"], params["T_ref"]
    )
    return max(k_eff * f_t * o2, 0.0)


def h2o2_branch_fraction(
    y: FloatArray,
    schema: StateSchema,
    params: Mapping[str, float],
    route: str,
    *,
    bisulfite: float | None = None,
) -> float:
    """``route``'s share of the H2O2 flux, in [0, 1] — the QSS branch fraction.

    H2O2 has no state slot (D-138: 17 ms half-life), so its production is partitioned among its
    consumers *instantaneously* by ``k_i[S_i] / sum_j k_j[S_j]``. Two consumers are tracked:

    * ``"ethanol"`` — weight ``k_ethanol_oxidation``, first-order (ethanol sits at ~100 g/L
      through aging and is effectively constant, exactly as D-71 argued for the pre-cascade form)
    * ``"bisulfite"`` — weight ``k_so2_oxidation * [HSO3-]``

    Both weights are the **pre-cascade constants**, re-used rather than re-fitted, so the
    ethanol:bisulfite competition on this node is the one D-72/D-73 already calibrated.

    The ethanol weight is a positive constant, so the denominator is never zero during aging and
    H2O2 always has somewhere to go — no floor term is needed here (unlike the quinone side,
    where :data:`k_quinone_polymerization` plays that role).
    """
    hso3 = free_bisulfite(y, schema, params) if bisulfite is None else bisulfite
    # D-172: the ethanol weight is still "the pre-cascade constant, re-used rather than re-fitted"
    # — that constant is now the ethanol half of the always-on total, formed rather than read.
    w_ethanol, _ = o2_depletion_shares(params)
    w_bisulfite = params["k_so2_oxidation"] * hso3
    total = w_ethanol + w_bisulfite
    if total <= 0.0:
        return 0.0
    return (w_ethanol if route == "ethanol" else w_bisulfite) / total


class OxygenActivation(Process):
    """The cascade's O2 *activation* node: Fe(II) + O2 -> H2O2 + Fe(III) (decision D-141).

    **No longer the sole O2 consumer** (D-196), and this docstring said it was in two places.
    :class:`PeroxideEthanolOxidation` now books the second O2 the 1-hydroxyethyl radical takes.
    Gate 1 (D-137) is **not** breached by that: its content is that ground-state wine molecules
    — ethanol, bisulfite, phenols — do not react with O2, so O2 goes through iron. A
    carbon-centred *radical* reacting with O2 at near-diffusion rate is not a counterexample to
    that claim; it is the standard next step of the chain Gate 1's own source describes. What is
    no longer true is the *bookkeeping* sentence, so it is corrected here rather than reinterpreted.

    Every gram of O2 that enters the **oxidant-producing** path still enters through here.
    Downstream Processes consume the **oxidants** this makes, never ``o2`` itself — with the one
    radical exception above, which consumes O2 without producing an oxidant, and that asymmetry
    is exactly why it moves the SO2:O2 ratio.

    Touches ``o2`` (consumed) and ``quinone`` (produced). H2O2 is *not* deposited anywhere — it
    is quasi-steady-state, so its consumers recompute the production flux from ``activation_rate``
    and take their branch share of it. That is why this Process's ``touches`` is only two slots
    even though it drives five downstream ones.

    Both slots are off every conservation ledger, so this moves nothing conserved: ``o2`` by the
    D-71 decision, ``quinone`` by fork D2 (its carbon comes from the untracked o-diphenol pool
    lumped into the rate constants, exactly as ``A420``'s pigment carbon does).
    """

    name = "oxygen_activation"
    tier = Tier.SPECULATIVE
    #: Consumes dissolved O2 and deposits the o-quinone oxidant — both off every ledger.
    touches = ("o2", "quinone")
    #: The activation node's own constants (aging.yaml, D-141) plus D-134's copper multiplier,
    #: re-homed here from PhenolicBrowning. The pKa/binding params read through ``acidbase`` to
    #: derive free bisulfite are omitted — all plausible, and this Process is already speculative,
    #: so they add no tier headline (the SulfiteOxidation/MalolacticConversion rule).
    reads: tuple[str, ...] = (
        "k_o2_depletion_total",
        "k_activation_phenolic",
        "k_activation_bisulfite",
        "E_a_activation",
        "k_copper_multiplier",
        "copper_typical",
        "T_ref",
        # Indirect, via the _free_so2 helper's ph_of_state + bisulfite_so2_at_ph (D-160).
        *PH_SYSTEM_READS,
        *SO2_BINDING_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if "o2" not in schema or "quinone" not in schema:
            return d
        r_o2 = activation_rate(y, schema, params)
        if r_o2 <= 0.0:
            return d
        d[schema.slice("o2")] = -r_o2
        # One o-quinone per O2 activated (the 2 Fe(III) returned to Fe(II) by one o-diphenol).
        d[schema.slice("quinone")] = _QUINONE_PER_O2 * (r_o2 / M_O2) * M_QUINONE
        return d


class PeroxideEthanolOxidation(Process):
    """H2O2 + ethanol -> acetaldehyde, via Fenton (decision D-141; re-homes D-71).

    The Fenton limb: ``H2O2 + Fe(II) -> HO.``, and the hydroxyl radical abstracts from ethanol to
    the 1-hydroxyethyl radical, which oxidises to acetaldehyde. D-138 recorded this as the one
    node the sources settle outright — Nguyen's Table 4.1 *is* its rate constant.

    **What changes from D-71, beyond the substrate.** The pre-cascade form produced
    ``y_acetaldehyde_per_o2`` (1.5) mol acetaldehyde per mol O2 — a *fitted yield*, because the
    old structure gave no mechanism to derive one from. Here the stoichiometry is mechanistic:
    one acetaldehyde per H2O2, and one H2O2 per O2. So the per-O2 acetaldehyde yield stops being
    a parameter and becomes a consequence of the branch fraction — it is 1.0 when ethanol wins
    the whole H2O2 node and less when bisulfite competes for it. That is a **fidelity gain and a
    magnitude change**, and the magnitude change is expected: it is why D-139 listed
    ``test_aging*.py`` magnitudes among the loud reds.

    Carbon closes exactly, unchanged from D-71: the acetaldehyde carbon is borrowed C2-for-C2
    from ``E`` (the D-27 reduction reversed).
    """

    name = "peroxide_ethanol_oxidation"
    tier = Tier.SPECULATIVE
    #: Books the oxidised carbon as ``acetaldehyde``, borrowed carbon-exactly from ``E``, and
    #: — since D-196 — the **second O2** the 1-hydroxyethyl radical takes. This comment used to
    #: read "does NOT touch ``o2`` — that is OxygenActivation's alone, which is the whole point
    #: of the rebuild", and that was a real claim, not a formality: it is corrected rather than
    #: quietly edited. The rebuild's point survives it (see :data:`_O2_PER_ACETALDEHYDE`).
    touches = ("o2", "acetaldehyde", "E")
    reads: tuple[str, ...] = (
        "k_o2_depletion_total",
        "k_activation_phenolic",
        "k_activation_bisulfite",
        "E_a_activation",
        "k_copper_multiplier",
        "copper_typical",
        "f_ethanol_o2_share",
        "k_so2_oxidation",
        "T_ref",
        # Indirect, via the _free_so2 helper's ph_of_state + bisulfite_so2_at_ph (D-160).
        *PH_SYSTEM_READS,
        *SO2_BINDING_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if "o2" not in schema or "acetaldehyde" not in schema:
            return d
        bisulfite = free_bisulfite(y, schema, params)
        r_o2 = activation_rate(y, schema, params, bisulfite=bisulfite)
        if r_o2 <= 0.0:
            return d
        share = h2o2_branch_fraction(y, schema, params, "ethanol", bisulfite=bisulfite)
        if share <= 0.0:
            return d
        # mol H2O2/L/h taking the ethanol branch -> mol acetaldehyde/L/h, 1:1.
        n_acet = _ACETALDEHYDE_PER_H2O2 * share * _H2O2_PER_O2 * (r_o2 / M_O2)
        rate = n_acet * M_ACETALDEHYDE
        d[schema.slice("acetaldehyde")] = rate
        d[schema.slice("E")] = -rate * _ETHANOL_PER_ACETALDEHYDE
        # The radical's own O2 (D-196). Weighted by ``share``, so it vanishes as bisulfite takes
        # the H2O2 node — which is what keeps Danilewicz's 1:1 and 1:2 asymptotes intact while a
        # same-sized FIXED fraction on the activation node breaks both (measured, D-196 3).
        d[schema.slice("o2")] = -_O2_PER_ACETALDEHYDE * n_acet * M_O2
        return d


class PeroxideSulfiteOxidation(Process):
    """H2O2 + HSO3- -> sulfate (decision D-141) — the FIRST half of D-72's split.

    D-138: ``SulfiteOxidation`` **splits in two**, and this is the node that survives when the
    quinone route is blocked. Danilewicz 2016 measured that limb at **1:1** — one SO2 per O2 —
    using benzenesulfinic acid to trap the quinones, and that is exactly what this Process alone
    produces (one H2O2 per O2, one SO2 per H2O2, when bisulfite wins the branch).

    **This split is the whole test.** One Process cannot produce Danilewicz's 1:1 / 1:2 / 1:1.7
    series; two can, and the 1.7 must *emerge* from partial quinone capture in
    :class:`QuinoneSulfonation`. A fitted 1.7 destroys the test — said at D-138, at D-139, and
    again here.

    Both slots off every ledger (SO2 is carbon-free, sulfate untracked), so nothing conserved
    moves. Distinct from the D-47 acetaldehyde-SO2 *binding*, which reversibly sequesters
    ``so2_total`` without removing it — the two do not double-count.
    """

    name = "peroxide_sulfite_oxidation"
    tier = Tier.SPECULATIVE
    touches = ("so2_total",)
    reads: tuple[str, ...] = (
        "k_o2_depletion_total",
        "k_activation_phenolic",
        "k_activation_bisulfite",
        "E_a_activation",
        "k_copper_multiplier",
        "copper_typical",
        "f_ethanol_o2_share",
        "k_so2_oxidation",
        "T_ref",
        # Indirect, via the _free_so2 helper's ph_of_state + bisulfite_so2_at_ph (D-160).
        *PH_SYSTEM_READS,
        *SO2_BINDING_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if SO2_STATE_KEY not in schema or "o2" not in schema:
            return d
        bisulfite = free_bisulfite(y, schema, params)
        if bisulfite <= 0.0:  # all free SO2 sequestered by carbonyls => nothing reactive left
            return d
        r_o2 = activation_rate(y, schema, params, bisulfite=bisulfite)
        if r_o2 <= 0.0:
            return d
        share = h2o2_branch_fraction(y, schema, params, "bisulfite", bisulfite=bisulfite)
        if share <= 0.0:
            return d
        n_so2 = _SO2_PER_H2O2 * share * _H2O2_PER_O2 * (r_o2 / M_O2)
        d[schema.slice(SO2_STATE_KEY)] = -n_so2 * M_SO2
        return d


class QuinoneSulfonation(Process):
    """o-quinone + bisulfite -> o-diphenol (major) / sulfonate (minor) — D-141, retitled D-199.

    Bisulfite consumes the quinone and is itself consumed. This is the node that lifts
    Danilewicz's series from 1:1 to 1:2 when it is open, and the node whose *partial* capture
    produces the real-wine 1:1.7 — because the quinone it does not take goes to the competing
    nucleophiles (Strecker, anthocyanin, ellagitannin) and to polymerisation.

    **The PRODUCT was named wrongly until D-199.** ~~"Bisulfite adds nucleophilically to the
    quinone ring"~~ describes the *minor* route. *Understanding Wine Chemistry* 2nd ed. §24.4.3.2,
    on an o-quinone **at wine pH**: *"most of the quinone (>90%) was reduced back to an o-diphenol,
    but the reaction yielded a small fraction of the sulfonate"*. (The ~25 %-sulfonate figure in
    the same passage is p-quinone at low pH, and at neutral pH the sulfonate becomes the only
    product — so the split is pH-dependent and wine sits at the reduction end.)

    **Nothing numerical moves, by construction**: either product consumes one bisulfite per
    quinone, so ``_SO2_PER_QUINONE`` stays 1.0 and the emergent 1:2 limit is untouched, and both
    products are off every ledger anyway. What changes is what the model is *claiming*.

    **And the correction runs in fork D2's favour.** :class:`OxygenActivation` never debits the
    o-diphenol that becomes the quinone (D-71/fork D2: "its carbon comes from the untracked
    o-diphenol pool"). Under reduction that pool is *regenerated*, so not debiting it is right for
    >90 % of this route's flux — the source corroborates the fork rather than breaking it. The
    residual sulfonate route is the part that genuinely removes a phenolic permanently and is not
    booked; recorded at D-199 rather than built, because no shipped pool is this route's substrate
    (``hydroxycinnamics``/``tannin`` are activation *drivers*, not the o-diphenol consumed here).

    **1.7 is an output of this module, never an input to it.** It is set by how much quinone this
    Process wins against its siblings, which is set by each sibling's own pre-cascade rate
    constant. Nothing here is tuned to reproduce it.

    This is also where SO2's *colour* protection becomes emergent: quinone captured here never
    reaches :class:`QuinonePolymerization`, so a sulfited wine browns less without any explicit
    "SO2 protects colour" term existing anywhere in the model.
    """

    name = "quinone_sulfonation"
    tier = Tier.SPECULATIVE
    touches = ("quinone", "so2_total")
    #: The pH-system and SO₂-binding sets are the indirect reads of this Process's own
    #: ``_free_so2`` call — ``ph_of_state`` + ``bisulfite_so2_at_ph`` (decision D-160).
    reads: tuple[str, ...] = (
        "k_so2_oxidation",
        "E_a_activation",
        "T_ref",
        *PH_SYSTEM_READS,
        *SO2_BINDING_READS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if SO2_STATE_KEY not in schema or "quinone" not in schema:
            return d
        quinone = _pool(y, schema, "quinone")
        if quinone <= 0.0:
            return d
        bisulfite = free_bisulfite(y, schema, params)
        if bisulfite <= 0.0:
            return d
        f_t = arrhenius_factor(
            float(y[schema.slice("T")][0]), params["E_a_activation"], params["T_ref"]
        )
        # Bilinear in the quinone pool and its nucleophile, at the SAME constant D-72 calibrated
        # against o2. NOTE, measured at D-141 and stated here because an earlier version of this
        # comment claimed the opposite: re-using the constant is NOT magnitude-neutral. The rate
        # is rescaled by [quinone]_ss / [o2]_ss (~0.06 unsulfited, scenario-dependent), which is
        # why the cascade is not the default wiring — see media._OXIDATIVE_SETS.
        r_q = params["k_so2_oxidation"] * f_t * quinone * bisulfite  # g quinone/L/h
        d[schema.slice("quinone")] = -r_q
        d[schema.slice(SO2_STATE_KEY)] = -_SO2_PER_QUINONE * (r_q / M_QUINONE) * M_SO2
        return d


class QuinoneStreckerDegradation(Process):
    """o-quinone oxidatively deaminates amino acids -> Strecker aldehydes (D-141; re-homes D-75).

    D-138 called this node **self-corroborated**: :class:`~aging.StreckerDegradation`'s own
    docstring already said "O2 (via quinones)" while its rate law drew directly on ``o2``. The
    sink knew its node; the rate law never followed. This is that correction.

    Everything downstream of the oxidant swap is D-75/D-100 unchanged: each aldehyde is throttled
    by its OWN precursor's relative-depletion gate, each draws its carbon from that precursor,
    and the nitrogen is deaminated to ``N``. Carbon and nitrogen close to machine precision.
    """

    name = "quinone_strecker_degradation"
    tier = Tier.SPECULATIVE
    touches = (
        "quinone",
        "methional",
        "phenylacetaldehyde",
        "CO2",
        "N",
        "methionine",
        "phenylalanine",
    )
    #: The last two entries are the D-100 relative-depletion gate's inputs, read indirectly
    #: through ``depletion_gate`` (decision D-160): the shared ``K_amino_acids`` half-saturation
    #: and each routed precursor's must-spectrum share. Derived from ``_STRECKER_ROUTES`` and the
    #: spec registry rather than re-listed, so adding a third Strecker route cannot leave the
    #: declaration behind.
    reads: tuple[str, ...] = (
        "k_strecker",
        "E_a_strecker",
        "y_strecker_per_o2",
        "f_methional",
        "T_ref",
        "K_amino_acids",
        *(SPEC_BY_SPECIES[precursor].fraction_param for _pool, precursor in _STRECKER_ROUTES),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if "methionine" not in schema or "quinone" not in schema:
            return d
        quinone = _pool(y, schema, "quinone")
        if quinone <= 0.0:
            return d
        f_t = arrhenius_factor(
            float(y[schema.slice("T")][0]), params["E_a_strecker"], params["T_ref"]
        )
        driver = params["k_strecker"] * f_t * quinone  # g quinone/L/h before the gates
        f_meth = params["f_methional"]
        shares = {_METHIONAL_SPECIES: f_meth, _PHENYLACETALDEHYDE_SPECIES: 1.0 - f_meth}

        r_q = 0.0
        co2_mol = 0.0
        product_rates: list[tuple[str, float]] = []
        precursor_carbon: dict[str, float] = {}
        for pool, precursor in _STRECKER_ROUTES:
            gate_i = depletion_gate(y, schema, params, (SPEC_BY_SPECIES[precursor],))
            if gate_i <= 0.0:
                continue
            r_q_i = shares[pool] * driver * gate_i
            n_i = params["y_strecker_per_o2"] * (r_q_i / M_QUINONE)
            if n_i <= 0.0:
                continue
            rate_i = n_i * (M_METHIONAL if pool == _METHIONAL_SPECIES else M_PHENYLACETALDEHYDE)
            co2_i = _CO2_PER_STRECKER_ALDEHYDE * n_i * M_CO2
            r_q += r_q_i
            co2_mol += _CO2_PER_STRECKER_ALDEHYDE * n_i
            product_rates.append((pool, rate_i))
            precursor_carbon[precursor] = rate_i * carbon_mass_fraction(pool) + co2_i * (
                carbon_mass_fraction("CO2")
            )
        if not product_rates:
            return d

        nitrogen = sum(
            draw_precursor_carbon(d, schema, precursor, carbon)
            for precursor, carbon in precursor_carbon.items()
        )
        d[schema.slice("quinone")] = -r_q
        for pool, rate_i in product_rates:
            d[schema.slice(pool)] = rate_i
        d[schema.slice("CO2")] = co2_mol * M_CO2
        d[schema.slice("N")] = nitrogen
        return d


class QuinoneAnthocyaninFading(Process):
    """o-quinone couples with anthocyanin -> colourless faded pigment (D-141; re-homes D-81).

    D-138 flagged this node as a **genuine unsourced fork**: ``AnthocyaninFading``'s docstring
    names *both* a quinone-coupled and an H2O2-bleach mechanism, so it does not resolve its own
    branch. The quinone limb is taken here, on two grounds recorded as a judgement rather than a
    citation: quinone-anthocyanin coupling is the route the wine-colour literature is written
    about, and the H2O2 node is already fully claimed by ethanol and bisulfite, whose constants
    are calibrated. **This is the cascade's weakest branch assignment after
    ``k_quinone_polymerization``'s magnitude, and it is reversible in one line here.**

    The pool transfer is D-81 unchanged: ``faded_anthocyanin`` gains exactly what ``anthocyanin``
    loses, so the colour identity closes by construction. Both off every ledger.
    """

    name = "quinone_anthocyanin_fading"
    tier = Tier.SPECULATIVE
    touches = ("quinone", "anthocyanin", "faded_anthocyanin")
    reads: tuple[str, ...] = (
        "k_anthocyanin_fade",
        "E_a_anthocyanin_fade",
        "y_anthocyanin_per_o2",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if "anthocyanin" not in schema or "quinone" not in schema:
            return d
        quinone = _pool(y, schema, "quinone")
        anthocyanin = _pool(y, schema, "anthocyanin")
        if quinone <= 0.0 or anthocyanin <= 0.0:
            return d
        f_t = arrhenius_factor(
            float(y[schema.slice("T")][0]), params["E_a_anthocyanin_fade"], params["T_ref"]
        )
        r_q = params["k_anthocyanin_fade"] * f_t * quinone * anthocyanin  # g quinone/L/h
        faded = params["y_anthocyanin_per_o2"] * r_q  # mass yield, the D-81 idiom
        d[schema.slice("quinone")] = -r_q
        d[schema.slice("anthocyanin")] = -faded
        d[schema.slice("faded_anthocyanin")] = faded
        return d


class QuinoneEllagitanninOxidation(Process):
    """o-quinone oxidises oak ellagitannin — the sacrificial oak sink (D-141; re-homes D-78).

    The other of D-138's two genuinely unsourced forks: ellagitannin may be oxidised via
    quinone/Fe(III) or may scavenge radicals directly. The quinone limb is taken here for
    consistency with the rest of the cascade — an ellagitannin that scavenged radicals directly
    would have to attach to the H2O2/Fenton node, and D-138 found no source placing it there
    either. **Recorded as a fork resolved by consistency, not by evidence.**

    Its *function* is unchanged and is what matters for the oak axis: oak tannin competes for the
    oxidant and so protects the wine, which is why an oaked run consumes more O2 in total (D-139
    measured 6.35 vs 5.71 mg/L in beer) while browning less per gram of O2.
    """

    name = "quinone_ellagitannin_oxidation"
    tier = Tier.SPECULATIVE
    touches = ("quinone", "ellagitannin")
    reads: tuple[str, ...] = (
        "k_ellagitannin_oxidation",
        "E_a_ellagitannin_oxidation",
        "y_ellag_per_o2",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if "ellagitannin" not in schema or "quinone" not in schema:
            return d
        quinone = _pool(y, schema, "quinone")
        ellag = _pool(y, schema, "ellagitannin")
        if quinone <= 0.0 or ellag <= 0.0:
            return d
        f_t = arrhenius_factor(
            float(y[schema.slice("T")][0]), params["E_a_ellagitannin_oxidation"], params["T_ref"]
        )
        r_q = params["k_ellagitannin_oxidation"] * f_t * quinone * ellag  # g quinone/L/h
        d[schema.slice("quinone")] = -r_q
        d[schema.slice("ellagitannin")] = -params["y_ellag_per_o2"] * r_q
        return d


class QuinonePolymerization(Process):
    """Uncaptured o-quinone polymerises to brown pigment, raising A420 (decision D-141).

    **This is where A420 stops being an O2 yield and becomes a quinone fate**, which D-139 listed
    as a silent red: the colour and sensory readouts keep reading a slot whose *meaning* has
    changed. The value changes too, and it should — under the cascade, pigment forms from the
    quinone that no nucleophile captured, so SO2, amino acids, anthocyanin and oak tannin all
    protect colour by competing for the same oxidant, with no protective term written anywhere.

    Structurally this is the **always-available** fate: the only quinone consumer that is not
    gated on a co-substrate. Without it, an unsulfited beverage carrying no nucleophiles would
    accumulate quinone without bound. Its rate constant is the cascade's weakest number —
    ``k_quinone_polymerization`` is anchored to an *asserted* hours-to-days quinone lifetime
    (D-138), ~~pending Nikolantonaki & Waterhouse 2012~~ **— STRUCK (D-199): that paper's key
    figure is in hand (Figure 24.12, reprinted with permission in *Understanding Wine Chemistry*
    2nd ed. p. 331, on disk since before D-141) and it does NOT settle this constant.** It ranks
    *nucleophiles*; this Process is by construction the fate of quinone that **no** nucleophile
    captured, so the figure has no entry for it and can have none. See
    :data:`k_quinone_polymerization`'s note for what the figure does and does not adjudicate.

    Both slots off every ledger: ``quinone`` by fork D2, ``A420`` because an absorbance carries
    no mass at all (D-74's argument, which the re-home leaves intact).
    """

    name = "quinone_polymerization"
    tier = Tier.SPECULATIVE
    touches = ("quinone", "A420")
    reads: tuple[str, ...] = (
        "k_quinone_polymerization",
        "E_a_activation",
        "y_a420_per_o2",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if "quinone" not in schema or "A420" not in schema:
            return d
        quinone = _pool(y, schema, "quinone")
        if quinone <= 0.0:
            return d
        f_t = arrhenius_factor(
            float(y[schema.slice("T")][0]), params["E_a_activation"], params["T_ref"]
        )
        r_q = params["k_quinone_polymerization"] * f_t * quinone  # g quinone/L/h
        d[schema.slice("quinone")] = -r_q
        # A420 accumulates per MOLE of quinone polymerised (the D-74 per-mole-of-oxidant idiom,
        # with quinone substituted for O2). Monotonic and irreversible, as D-74 requires.
        d[schema.slice("A420")] = params["y_a420_per_o2"] * (r_q / M_QUINONE)
        return d
