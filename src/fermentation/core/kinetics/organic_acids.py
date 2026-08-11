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

**What the free prediction says (read the D-180 and D-181 records before quoting it), and the
scope it is true in.** Against a measured drop of 0.8125 pH units, SINCE D-181:

* at **nominal yields and floors**, across the sampled ``pKa_peptide_buffer`` band:
  **42.7-62.2 %** of it;
* over the **JOINT** band — the four yields and the three floors are sampled too, and widely:
  **8.7-81.4 %**.

**Those numbers were 63-92 % and 41-105 % at D-180, and they got worse on purpose.** D-180
closed with its agreement held open by two omitted terms of OPPOSITE sign, and D-181 built the
larger one — three wort acids a real ferment removes, which the model could not lose because
they were not state (:class:`WortAcidRemoval`). Removing them removes anion charge, so the
finished pH rises and the predicted drop shrinks. Nothing in the reachable band now covers the
measurement, where a corner used to. **That is the honest state**, and it is the reason the
D-180 agreement must never have been read as validation: it was propped up by an uncorrected
error pulling the other way.

One omitted term remains, and it pulls the opposite way again:

* dissolved CO₂ at end-of-fermentation saturation is not in the charge balance — worth about
  −0.3 pH. (The engine's ``CO2`` pool is *cumulative evolved gas*, not dissolved
  concentration, so it must NOT be read as carbonic acid; that needs a saturation model.)

It is a follow-up beat with a measured size, and it is not a reason to tune a yield. Note that
building it will move the headline back UP toward the measurement — which after D-181 would be
a real improvement rather than a compensation, and is precisely why the two were built in this
order rather than the other.

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

* **acetic** really peaks at ~200 ppm around day 2 (3.4× the wort level) and falls back to
  105-126 by day 7 as the yeast re-assimilate it. A monotone producer hits the endpoint and
  drops the peak — the same concession :mod:`~fermentation.core.kinetics.keto_acids` already
  makes for overflow pyruvate, and it needs the same excretion/re-assimilation *pair* to fix.
* **lactic** rises hardest *after* the extract curve goes flat (days 2-7). A flux-linked form
  makes it early instead. Tyrell name the candidate mechanism — end-of-ferment autolysis —
  and this engine has an autolysis Process, so a later beat has a route.

Reported, not tuned: one dataset cannot separate a late excretion from an autolytic release.

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
ORGANIC_ACID_SPECS: tuple[OrganicAcidSpec, ...] = (
    OrganicAcidSpec("acetic", "Y_acetic_sugar_beer", "acetic_acid"),
    OrganicAcidSpec("lactic", "Y_lactic_sugar_beer", "lactic_acid"),
    OrganicAcidSpec("succinic", "Y_succinic_sugar_beer", "succinic_acid"),
    OrganicAcidSpec("malic", "Y_malic_sugar_beer", "malic_acid"),
)


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
