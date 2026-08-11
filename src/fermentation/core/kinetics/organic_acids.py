"""Beer's organic acids — what the ferment MAKES (D-180) and what it REMOVES (D-181).

Until D-180 beer's acid slots (D-179) were a **composition**: dosed at pitch, touched by no
Process, so a modelled beer's pH sat exactly where ``initial_ph`` put it for the whole run
while a real beer's falls by most of a pH unit. This module gives four of those five slots a
producer, which turns beer's pH from an input into a **prediction**.

The other half of that change lives at the compile seam, and the two are inseparable: the
seeds moved from finished-beer levels to **wort** levels (``*_typical_wort``,
``beer_acids.yaml``). A producer bolted onto finished-beer seeds would finish at roughly twice
the measured beer — verified at pH 4.26, below any real beer — so "produce the acids" and
"start from a wort" are one decision, not two.

**The source, and why it is unusually good for this engine.** Tyrell et al. (2013),
*BrewingScience* 66:75-84, fermented ONE 12 °Plato wort with FOUR strains in EBC tubes and
published, for the same ferments, both every organic acid's day-0-to-day-7 course (their
Figs 6-14) *and* the pH and extract curves (their Fig 4). So the yields are a measured
**difference on one wort**, not two studies subtracted — and the pH curve is then a **free
prediction**, because nothing in the parameter file is fitted to it.

**What the free prediction says (read the D-180, D-181 and D-182 records before quoting it),
and the scope it is true in.** Against a measured drop of **0.81** pH units — the mean of
Tyrell's extreme strains, which is the denominator the code divides by; D-180's prose quotes
0.8125, the four-strain mean, and the two must not be mixed — SINCE D-182:

* at **nominal**, across the sampled ``pKa_peptide_buffer`` band: **77.6-97.0 %** of it;
* over the **JOINT** band — **nine** drawn quantities, not one (the four yields, the three
  floors, the three seeds, the two acid pKas and D-182's carbonic pKa plus two solubility
  constants, moved together per dimension): **63.8-109.4 %**.

**BOTH of the omitted terms D-180 named are now built, and the order they were built in is the
point.** D-180 closed with its agreement held open by two terms of OPPOSITE sign and sized
both. D-181 built the larger one first — three wort acids a real ferment removes, which the
model could not lose because they were not state (:class:`WortAcidRemoval`) — which made the
prediction agree WORSE, dropping it to 42.7-62.2 % at nominal and 7.6-82.2 % band-wide, with
nothing in the reachable band covering the measurement. D-182 then built the other one,
dissolved CO₂ as carbonic acid in the charge balance
(:func:`~fermentation.core.acidbase.dissolved_co2_molar`), and the headline came back up.

**Building them the other way round would have produced the same final number and a false
story about it**: the flattering term would have arrived at 76-104 % agreement with D-181's
same-sized error still in place, and nobody would have had a reason to look for it.

Two things about the recovered agreement, so it is not over-read. It is still a SHORTFALL at
nominal (3-22 %), and the corner that now reaches the measurement is a corner of a
nine-dimensional hypercube rather than a draw anyone observed. And the two terms are **not
additive**: dissolved CO₂ dissociates more as pH rises, so it buffers AGAINST the falling
acids' removal — the same three acids that were worth +0.2094 pH at D-181 are worth +0.1128 pH
beside the CO₂ term. Adding their separately-measured sizes over-counts.

**Mechanism — an extra sliver off ``S``, never uptake's yields.** There are two shipped ways
to fund a byproduct's carbon, and the choice matters:

* D-16's ``Y_glycerol_sugar``/``Y_byproduct_sugar`` rescale the ethanol/CO₂ split *inside*
  uptake. Beer holds both at 0 so its sugar→ethanol stays theoretical and its CO₂-ratio
  benchmark is byte-for-byte validated. Routing acids that way would move that benchmark.
* D-19's byproduct idiom pulls an *additional* sliver of ``S`` into the producing pool and
  leaves ``dX``/``dN``/``dE``/``dCO2`` untouched. All four acids together divert
  2.435e-3 g/g — about a fifth of wine's single ``Y_byproduct_sugar`` and ~0.24 % of the
  sugar — so the drift is of the same order as the fusel/ester draws already accepted.

This module takes the second. ``Byp`` therefore stays empty for beer and beer's
``Y_byproduct_sugar`` stays 0 — which is now **load-bearing rather than incidental**: beer
carries its own ``succinic`` slot *and* the charge balance reads ``Byp`` as
succinic-equivalent (:data:`~fermentation.core.acidbase.BYP_AS_SUCCINIC`), so a non-zero beer
``Y_byproduct_sugar`` would count succinic twice in the pH solve. Pinned by
``tests/test_organic_acids.py::test_beer_byproduct_yield_stays_zero_or_succinic_double_counts``.

**Flux-linked, and the shapes that costs.** Production tracks the fermentative sugar-uptake
rate, so the run integral is a clean yield (``Y · ΔS``) on any beer, and it stops at dryness.
Against Tyrell's own curves that is right for succinic and wrong in two named ways:

* **acetic** really peaks around day 2 and falls back to 105-126 by day 7. **D-183 moved it off
  this rate law** (:class:`AceticAcidOverflow`) after mapping Fig 13 onto Fig 4 showed 86 % of
  its rise inside the first 15 % of the sugar flux. What that fixes is *when* the acid appears,
  not the peak — see below.
* **lactic** rises hardest *after* the extract curve goes flat (days 2-7). A flux-linked form
  makes it early instead. Tyrell name the candidate mechanism — end-of-ferment autolysis —
  and this engine has an autolysis Process, so a later beat has a route.

Reported, not tuned: one dataset cannot separate a late excretion from an autolytic release.

**THE MID-FERMENT SPIKE IS STILL NOT MODELLED, AND ITS FIX WAS REFUSED ON MEASUREMENT (D-183).**
D-180 §9 proposed the :mod:`~fermentation.core.kinetics.keto_acids` template — an excretion /
re-assimilation *pair*, both flux-linked. **The source's own figures falsify both halves of it**:

* the *excretion* half, because production is growth-phase-confined (the 86 %/15 % above);
* the *re-assimilation* half, because **half the fall happens at zero fermentative flux** —
  between days 5 and 7 the extract curve is flat while acetic drops 141.5 → 117.75 ppm. A
  co-metabolic term freezes the pool exactly there, which is the D-181 argument one figure over.

So the removal was built as a probe and refused, for three measured reasons rather than a
preference. Tyrell's decline **cannot discriminate** first-order from constant-rate (every fit
gap is under the ±3 ppm figure-read tolerance); it **cannot identify a floor** (the mean's best
fit sits on the bound at floor = 0, and the four strains split 2-2); and every law it *is*
consistent with makes the endpoint a function of the solver horizon — 108 / 65 / 20 / 0 ppm at
days 7 / 14 / 30 / 400 for the pure first-order form. A modelled beer's pH would then depend on
how long the run was left going, which is worse than a missing transient. **Do not re-propose
the pair without a dataset that fixes the removal law**; a second time course on a *different*
wort would do it, because the strain spread is what the floor is currently confounded with.

**What shipping the production half alone actually bought**, stated so it is not over-read: the
shape error against Tyrell's measured days 1-7 falls from **61.6 to 32.5 ppm RMSE** and the
endpoint lands on the four-strain mean by construction — but the curve is still **monotone**. A
plateau reached early is not a transient.

**Tier: plausible.** The yields are measured differences on a real wort with named strains at
both band edges, but the divisor passes through two derived steps (real-vs-apparent
attenuation, and a Plato→SG linearisation) and every level is a *figure read* rather than a
printed table. Parameter-tier propagation (D-1) does the rest.

**Isolable (prime directive 3).** Beer-only — no wine Process reads these yields, so a wine
ensemble never draws them — and disabled at the compile seam whenever ``initial_ph`` is
absent, exactly like the acid slots it fills. An un-anchored beer is therefore byte-for-byte
the pre-D-179 beer, and cannot acquire a charge balance mid-run by producing into empty slots.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics.carbon_routing import (
    draw_carbon_from_sugar as _draw_carbon_from_sugar,
)
from fermentation.core.kinetics.carbon_routing import (
    fermentative_uptake_rates as _fermentative_uptake_rates,
)
from fermentation.core.kinetics.growth import biomass_growth_rate as _biomass_growth_rate
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


@dataclass(frozen=True)
class OrganicAcidSpec:
    """One produced organic acid: its state slot, its yield parameter, its carbon species.

    ``species`` is the key :func:`~fermentation.core.chemistry.carbon_mass_fraction` weights
    the sugar draw by, so each acid's carbon is booked at **its own** molecule — the D-99
    correction applied here from the start rather than after a lumped stand-in ships.

    **No molar mass here, deliberately.** It would duplicate
    :data:`~fermentation.core.acidbase.ALL_ACIDS`, which is where every g/L -> mol/L conversion
    in the charge balance already comes from; a second copy is the single-source-of-truth
    failure this repo pins elsewhere. Callers needing mass read
    ``ALL_ACIDS[spec.slot].molar_mass``.
    """

    slot: str
    yield_param: str
    species: str


#: The acids beer's yeast PRODUCES, and the single source of truth every layer derives from
#: (the ``ESTER_SPECS``/``FUSEL_SPECS`` registry discipline, D-96/D-99): the Process's
#: ``touches`` and ``reads``, the carbon draw, and the tests all read this tuple, so adding a
#: fifth produced acid is one entry rather than a new code path.
#:
#: **``citrate`` is deliberately absent and must stay absent.** Tyrell's section 1.5 ("citric
#: … basically depends on concentrations of wort"), their Table 2 ("final concentration mostly
#: determined by wort concentration") and their Fig. 12 (four strains scattering ±20 ppm around
#: the wort value with no trend) all say the same thing three ways: beer's citrate is
#: malt-derived, not a fermentation product. It keeps a wort seed and stays inert.
#:
#: The three acids that *fall* over a real ferment — pyruvic, formic, oxalic — are absent from
#: THIS registry because they are not produced; they have their own, :data:`WORT_ACID_SINKS`,
#: and their own Process (:class:`WortAcidRemoval`, D-181). Two registries rather than one with
#: a sign, because the two halves share no parameter, no rate law and no ledger treatment.
#:
#: **``acetic`` LEFT THIS REGISTRY AT D-183 and must not be added back.** It is produced, but
#: not on this rate law: Tyrell's Fig 4 puts 86 % of its rise inside the first 15 % of the
#: fermentative flux, so it has its own growth-linked producer (:class:`AceticAcidOverflow`).
#: A fifth entry here would double-count it *and* restore the shape the source falsifies.
ORGANIC_ACID_SPECS: tuple[OrganicAcidSpec, ...] = (
    OrganicAcidSpec("lactic", "Y_lactic_sugar_beer", "lactic_acid"),
    OrganicAcidSpec("succinic", "Y_succinic_sugar_beer", "succinic_acid"),
    OrganicAcidSpec("malic", "Y_malic_sugar_beer", "malic_acid"),
)

#: The acid whose production is **growth**-linked rather than flux-linked (decision D-183), and
#: the molecule its carbon is booked at. One acid, so a tuple registry would be ceremony — but
#: the slot and the species are named here rather than inlined for the same single-source reason
#: :data:`ORGANIC_ACID_SPECS` exists: the Process's ``touches``, the carbon draw and the tests
#: all read these two names.
ACETIC_SLOT = "acetic"
ACETIC_SPECIES = "acetic_acid"


def organic_acid_rates(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float]
) -> list[tuple[OrganicAcidSpec, float]]:
    """Each produced acid's rate ``d(<slot>)/dt`` [g/L/h], paired with its spec.

    ``Y_<acid> · Σ_i r_i`` over the per-slot fermentative uptake rates — so the run integral
    is ``Y · ΔS``, a yield in the ordinary sense, on a wine-shaped single sugar slot or beer's
    three alike. Returns ``[]`` when nothing is fermenting, which is what stops production at
    dryness without a separate gate.
    """
    flux = sum(_fermentative_uptake_rates(y, schema, params))
    if flux <= 0.0:
        return []
    return [(spec, params[spec.yield_param] * flux) for spec in ORGANIC_ACID_SPECS]


def organic_acid_carbon_draw(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float]
) -> float:
    """Total carbon [g C/L/h] the producer books out of sugar, across all four acids.

    Each acid contributes ``rate · carbon_mass_fraction(species)`` at its **own** molecule's
    fraction. Factored out of the Process for the same reason the fusel draw is (D-33/D-99):
    a future consumer that refunds part of this carbon to another source must be able to
    compute the identical number rather than re-derive it.
    """
    return float(
        sum(
            rate * carbon_mass_fraction(spec.species)
            for spec, rate in organic_acid_rates(y, schema, params)
        )
    )


class OrganicAcidExcretion(Process):
    """Yeast-produced organic acids, flux-linked, carbon-routed out of sugar (decision D-180).

    ``d(<acid>)/dt = Y_<acid> · Σ_i r_i`` for each entry of :data:`ORGANIC_ACID_SPECS`, with
    the acids' carbon drawn out of ``S`` via
    :func:`~fermentation.core.kinetics.carbon_routing.draw_carbon_from_sugar` (option a1,
    D-19). So this touches the four acid slots and ``S`` — never ``E``/``CO2`` — and
    ``total_carbon`` closes exactly, because the carbon removed from sugar equals the carbon
    deposited in the pools by construction.

    The rates come from the shared
    :func:`~fermentation.core.kinetics.carbon_routing.fermentative_uptake_rates`, which the
    uptake Process also consumes. That sharing is what makes ``Y`` a yield rather than a
    coincidence, and it carries an obligation: **every RateModifier that scales uptake must
    also name this Process**, or the acid would be booked against a flux the solver never ran.
    Beer's medium definition does exactly that (``ArrheniusTemperature.for_uptake`` with this
    Process as an extra target) — the D-32 correctness coupling, in the form it takes here.

    **No CO₂ co-product**, matching the fusel sugar stand-in: routing sugar carbon into an
    acid pool asserts no named reaction, so no decarboxylation is claimed. Where this engine
    *does* assert a reaction (the D-106 Ehrlich branch) it charges the CO₂; here it does not.
    Real succinate formation via the reductive TCA branch and the GABA shunt draws on
    glutamate as much as on hexose, so the sugar source is a bookkeeping stand-in — the
    ledger closes, the metabolic claim is not made.
    """

    name = "organic_acid_excretion"
    tier = Tier.PLAUSIBLE
    #: The four produced acid slots plus ``S``, derived from the registry so a fifth acid
    #: cannot silently violate the ``touches`` contract.
    touches: tuple[str, ...] = (*(spec.slot for spec in ORGANIC_ACID_SPECS), "S")
    #: The four yields, plus the three uptake constants the shared rate helper reads —
    #: declared here too, because ``reads`` has two masters (tier propagation AND sampler
    #: scope, D-160): an undeclared read does not merely under-document a dependency, it
    #: narrows the reported ensemble spread below what the parameter's own provenance
    #: justifies.
    reads: tuple[str, ...] = (
        *(spec.yield_param for spec in ORGANIC_ACID_SPECS),
        "q_sugar_max",
        "K_sugar_uptake",
        "K_repression",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        rates = organic_acid_rates(y, schema, params)
        if not rates:
            return d
        for spec, rate in rates:
            d[schema.slice(spec.slot)] = rate
        # ONE draw for all four, each at its own molecule's carbon fraction.
        _draw_carbon_from_sugar(d, y, schema, organic_acid_carbon_draw(y, schema, params))
        return d


class AceticAcidOverflow(Process):
    """Acetic acid, produced with **growth** rather than with the sugar flux (decision D-183).

    ``d(acetic)/dt = Y_acetic_biomass_beer · dX/dt``, with the acid's carbon drawn out of ``S``
    at acetic acid's own C2 fraction, exactly as :class:`OrganicAcidExcretion` draws the other
    three. Touches ``acetic`` and ``S`` and nothing else — in particular **not** ``X`` or ``N``:
    it *reads* the growth rate, it does not add to it.

    **The rate law is the whole content of this Process, and it is sourced twice.** Until D-183
    acetic rode :data:`ORGANIC_ACID_SPECS` with the other three, i.e. a yield on the fermentative
    sugar flux. Mapping Tyrell's Fig 13 (the acid) onto their Fig 4 (the extract of the *same*
    ferments) with this file's own sugar divisor falsifies that:

    * by day 1 the wort is **15 %** attenuated and acetic has already made **86 %** of its whole
      rise (59 → 145 ppm of a 59 → 170.25 peak);
    * their Table 2 scores acetic ``++`` (strong increase) at **lower** Krausen, ``o`` at high
      Krausen and ``-`` at Krausen collapse — production over before the bulk of the sugar goes.

    Nothing proportional to the sugar flux can do that. The shipped flux-linked form put acetic
    at 65.7 ppm on day 1 and 74.4 on day 2 against a measured 145.0 and 170.25. Growth is the
    denominator the source itself reaches for elsewhere: their Fig 15 normalises acid production
    "related to grown cells" beside "related to fermented extract".

    **What this fixes, and what it explicitly does NOT.** Anchored on the day-7 level, a modelled
    beer now rises to 117.75 ppm by about day 1 and holds, which halves the shape error against
    Tyrell's own days 1-7 (RMSE 61.6 → 32.5 ppm). It is still **monotone**: the measured
    rise-then-fall is not reproduced, and this Process delivers **no mid-ferment spike**. The
    re-assimilation half that would deliver one was built as a probe, measured, and **refused**
    — see the module docstring. Do not read the improved RMSE as the transient being modelled.

    **The D-32 correctness coupling, in the form it takes here.** This recomputes the *base*
    growth rate from the shared :func:`~fermentation.core.kinetics.growth.biomass_growth_rate`,
    the same helper the growth Process builds biomass with. Growth's realised rate is scaled by
    ``ArrheniusTemperature.for_growth``, so this Process **must** be named as an extra target of
    that modifier or its yield would be booked against a growth flux the solver never ran — a
    cold beer would grow slowly while making acetic at the warm rate. Beer's medium does exactly
    that. Note the target moved: the flux-linked form belonged to ``for_uptake``'s extra targets,
    this one belongs to ``for_growth``'s, and leaving it on the old one would be silently wrong
    rather than loudly broken. Beer wires no carrying-capacity modifier (that is wine-only), so
    the growth Arrhenius is the complete list.

    Tier **speculative**, one step below :class:`OrganicAcidExcretion`'s plausible, and the
    reason is the denominator: the numerator is the same measured per-strain delta, but the
    biomass it is divided by is *this model's*, not Tyrell's. That drops ``acetic``'s output tier
    to speculative by parameter-tier propagation (D-1) — which changes no reported tier, because
    beer's ``ph_tier`` is already speculative through the peptide buffer (D-179).
    """

    name = "acetic_acid_overflow"
    tier = Tier.SPECULATIVE
    #: The acid slot and ``S``. NOT ``X``/``N`` — reading the growth rate is not contributing to
    #: it, and adding them would make ``ProcessSet(strict=True)`` permit a write that would
    #: double-count biomass.
    touches: tuple[str, ...] = (ACETIC_SLOT, "S")
    #: The yield plus the three constants
    #: :func:`~fermentation.core.kinetics.growth.biomass_growth_rate` reads. Declared here too
    #: because ``reads`` has two masters — tier propagation AND sampler scope (D-160): leaving
    #: the growth constants undeclared would narrow a beer ensemble's reported spread below what
    #: this Process's own dependence justifies.
    reads: tuple[str, ...] = ("Y_acetic_biomass_beer", "mu_max", "K_s", "K_n")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        growth = _biomass_growth_rate(y, schema, params)
        if growth <= 0.0:  # no growth ⇒ no overflow, which is what stops this at N exhaustion
            return d
        rate = params["Y_acetic_biomass_beer"] * growth  # [g acetic/L/h]
        d[schema.slice(ACETIC_SLOT)] = rate
        _draw_carbon_from_sugar(d, y, schema, rate * carbon_mass_fraction(ACETIC_SPECIES))
        return d


@dataclass(frozen=True)
class WortAcidSinkSpec:
    """One wort acid that is REMOVED during fermentation: its state slot and its floor.

    Deliberately a different dataclass from :class:`OrganicAcidSpec` rather than the same one
    with a sign: a produced acid needs a yield and a carbon ``species``, a removed one needs a
    floor and has no carbon destination at all (see :class:`WortAcidRemoval`). Sharing a spec
    would have made the two halves look interchangeable when their only common element is that
    both end up in beer's charge balance.
    """

    slot: str
    floor_param: str


#: The acids beer's wort carries and LOSES (decision D-181) — the sink half of beer's acid
#: model, and the single source of truth its Process, its ``touches``/``reads`` and its tests
#: all derive from.
#:
#: All three are OFF EVERY LEDGER, which is why no ``species`` appears here and why none of
#: them is in :data:`~fermentation.core.chemistry.MOLAR_MASS`. That is not an omission to fix
#: later: it is the structural form of "this carbon leaves the beer by a route the source does
#: not attribute" (the ``iso_alpha`` precedent, D-64). A future beat that gives one of these a
#: producer drawing on ``S`` would need to add the weight *and* the species, and until it does
#: the missing key raises instead of leaking carbon silently.
WORT_ACID_SINKS: tuple[WortAcidSinkSpec, ...] = (
    WortAcidSinkSpec("pyruvic", "pyruvic_floor_beer"),
    WortAcidSinkSpec("formic", "formic_floor_beer"),
    WortAcidSinkSpec("oxalic", "oxalic_floor_beer"),
)


class WortAcidRemoval(Process):
    """The three wort acids that FALL — beer's missing base, built (decision D-181).

    ``d(<acid>)/dt = −k_wort_acid_removal · max(<acid> − floor_<acid>, 0)`` for each entry of
    :data:`WORT_ACID_SINKS`. First-order relaxation toward a **measured** floor, temperature-
    flat, touching nothing but the three acid slots.

    **What it is for.** D-180 turned beer's pH into a prediction by producing the four acids
    that rise, and closed by naming two omitted terms of opposite sign that were holding its
    agreement open. This is the larger of the two: three acids a real wort carries and a real
    ferment removes (pyruvic 22 → ~1.3, formic 26 → ~4.75, oxalic 22 → ~5.6 mg/L), which the
    model could not lose because they were not state. Removing them removes anion charge, which
    raises the finished pH, which makes the predicted drop **smaller** — so this Process makes
    the headline number agree WORSE with the measurement, on purpose. Building the other omitted
    term (dissolved CO₂, the opposite sign) first would have moved the same number to near-exact
    agreement with an uncorrected error still in place.

    **D-182 has since built that other term, and it HALVES what this one appears to be worth**
    — not because these acids shrank, but because carbonic acid dissociates more the higher the
    pH, so it pushes back against exactly the rise this Process causes. In isolation these three
    acids move the finished pH by +0.2094; beside the CO₂ term, +0.1128.
    ``test_removing_the_falling_acids_raises_the_finished_ph_by_the_predicted_amount`` asserts
    BOTH numbers, so the buffering explanation stays falsifiable rather than merely asserted.

    **Why first-order-to-a-floor and NOT the flux-linked idiom this module's other Process
    uses.** Tyrell's Table 2 scores all three acids ``--`` (strong decrease) at *lower* Krausen
    and ``0`` at both high Krausen and Krausen collapse: the fall is over within about a day and
    nothing happens for the remaining six. A term riding
    :func:`~fermentation.core.kinetics.carbon_routing.fermentative_uptake_rates` is proportional
    to biomass and therefore peaks MID-ferment — it would put the removal exactly where the
    source says there is none. The archive's default is flux-linking; here the default is wrong,
    and this paragraph exists so it is not "restored".

    **No mechanism is asserted, and the ledger treatment is what says so.** Yeast uptake,
    calcium-oxalate precipitation and adsorption onto cell walls would all produce this curve;
    Tyrell distinguish none of them and explicitly report that pyruvate re-assimilation was *not*
    confirmed in their own lab-scale arm. So the acids' carbon is not routed to ``E``/``CO2``
    (which would claim metabolism) nor to a precipitate pool (which would claim precipitation):
    the three slots sit OFF every ledger, the ``iso_alpha`` treatment for exogenous mass that
    leaves the liquid by an unattributed route. ``total_carbon`` closes exactly as before, and
    **not** because this Process balances — because it touches nothing the ledger weighs. The
    price is ~19 mg C/L of malt carbon per litre untracked against ~33 g C/L on the ledger.

    **Beer-only and opt-in.** No wine Process reads these parameters and wine carries none of
    these slots. Like D-180's producer, the compile seam DISABLES this Process when a beer
    scenario supplies no ``initial_ph`` — not tidiness but the D-179 correctness gate: the acid
    slots are 0 without that opt-in, and a Process free to run on them would have nothing to
    remove but would still hold the empty slots' tier below VALIDATED.

    Tier **plausible** — unlike :class:`~fermentation.core.kinetics.hops.IsoAlphaAcidLoss`,
    whose structurally identical loss law is SPECULATIVE because its rate is an author estimate,
    every number here (both endpoints and the rate) is read off a measured time course on a real
    ferment. What is speculative is the *mechanism*, and the model does not claim one.
    """

    name = "wort_acid_removal"
    tier = Tier.PLAUSIBLE
    #: The three sink slots and nothing else — no ``S``, no ``E``, no ``CO2``. The shortest
    #: ``touches`` of any producing/consuming Process in the engine, and that brevity is the
    #: no-mechanism claim in machine-checkable form.
    touches: tuple[str, ...] = tuple(spec.slot for spec in WORT_ACID_SINKS)
    #: The shared rate constant plus the three floors. Derived from the registry so a fourth
    #: sink cannot be added without its floor entering the sampler's scope (``reads`` has two
    #: masters — tier propagation AND sampler scope, D-160).
    reads: tuple[str, ...] = (
        "k_wort_acid_removal",
        *(spec.floor_param for spec in WORT_ACID_SINKS),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        k = params["k_wort_acid_removal"]
        for spec in WORT_ACID_SINKS:
            pool = float(y[schema.slice(spec.slot)][0])
            # Clamped at the floor rather than at zero: an acid already at (or, through solver
            # undershoot, below) its floor has no removal term, so the pool cannot be driven
            # negative and cannot be REFILLED by a sign flip either.
            excess = pool - params[spec.floor_param]
            if excess <= 0.0:
                continue
            d[schema.slice(spec.slot)] = -k * excess
        return d
