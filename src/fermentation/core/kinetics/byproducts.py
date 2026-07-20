"""Temperature-/metabolism-driven aroma byproducts — esters and fusel alcohols.

The first Milestone-2 (Tier-2) beat (decision D-18 build order). Two *additive*
Processes that fill the produced-only ester and ``fusels`` pools the schema gained in the
byproducts beat. Since **decision D-96** "the ester pool" is three single-molecule pools —
``ethyl_acetate`` / ``isoamyl_acetate`` / ``ethyl_hexanoate``, registered in
:data:`~fermentation.core.kinetics.carbon_routing.ESTER_SPECS` — not one lump; see that
registry for why the lump could not survive contact with its own OAV readout. Both are
trace (mg/L–low-hundreds-mg/L) beside the g/L ethanol flux, and both **rise with
temperature** — the physics behind the
"warm ferments are estery/fusel-heavy" rule and the directional benchmark
``test_lower_temperature_is_slower_but_cleaner``.

**Additive, not multiplicative.** Unlike ethanol inhibition / Arrhenius (which
*scale* an existing flux and so are :class:`~fermentation.core.process.\
RateModifier` objects), these *produce* a compound, so they are ordinary summed
:class:`~fermentation.core.process.Process` objects (handoff §3.2). Each embeds its
own temperature sensitivity via :func:`~fermentation.core.kinetics.arrhenius.\
arrhenius_factor` (the shared Arrhenius shape, D-11) with its *own* activation
energy — kept in the Process rather than as a separate modifier so each byproduct
mechanism is one self-contained, independently-togglable object.

**Temperature directions — what rises, what falls, and why (decisions D-19 → D-21).**
Each *production* rate is the fermentative-flux Monod shape times its Arrhenius factor,
so over a run to dryness the total produced scales as ``f_byproduct(T) / f_flux(T)``
(the flux integral to dryness is fixed — the must ferments out either way), which rises
with T **iff** the byproduct's activation energy exceeds the fermentative flux's.

* **Fusels (both media)** must rise with T for the "cleaner-when-colder" benchmark, so
  ``E_a_fusels`` is held **above** ``E_a_uptake`` — a *sourced ordering* (Mouret 2015:
  higher alcohols rise with T; brewing/enology consensus). This is the half of the
  benchmark that genuinely carries "cleaner cold" for wine.
* **Esters** are *not* a simple production pool: liquid ester is the balance of
  synthesis (:class:`EsterSynthesis`, ``E_a_esters``) against the gas-stripping sink
  (:class:`EsterVolatilization`, stripping sensitivity ``E_a_uptake + dH_ester_volatil``
  ≈ 100 kJ/mol, the *same physical value in both media*). The wine/beer direction lives
  in **synthesis ``E_a_esters``, set per medium and sourced** — **beer** steep
  (de Andrés-Toro 1998, ester ride growth, apparent E_a ~200+ kJ/mol ⇒ synthesis wins ⇒
  liquid esters *rise* with T) and **wine** weak/~flat (Mouret 2015; Rollero 2014, wine
  ester synthesis is weak and non-monotonic ⇒ stripping wins ⇒ liquid esters *fall* with
  T, the Rollero evaporation inversion). See :class:`EsterVolatilization` for why a
  *sourced* (hence medium-independent) stripping forces the direction into synthesis.

Exact magnitudes stay speculative (directional check only, handoff §3.5); the *orderings*
are sourced. The integrated directions are **verified empirically** at 14/20/25 °C per
medium, not assumed — growth's sugar draw, the finite dryness cutoff, the inactivation
brake, and the stripping dynamics all perturb the clean cancellation.

**Carbon accounting — option (a)/a1 (decision D-19): carbon routed from sugar.**
Each Process draws its species' carbon *out of ``S``* and the pools are weighted in
``total_carbon`` (each ester by its OWN molecule's carbon fraction since D-96 — C4/C7/C8 —
and fusels by isoamyl alcohol's), so esters and fusels are **real carbon-accounted state**
under one consistent rule with ``Gly`` and ``Byp`` (D-16) — not diagnostic re-expressions.
This is the user's call (2026-06-29),
chosen over the interim option (b) (pools outside ``total_carbon``) and the
closure-neutral a2 variant (transfer from ``E``/``Byp`` without a sugar draw).

*The draw touches only ``S`` — never ``E``/``CO2``.* The uptake Process still ferments
``S`` to ethanol+CO2 unchanged; these Processes pull an *additional* sliver of ``S``
into their own pool. So at the derivative level only ``dS`` gains a (tiny, negative)
term — ``dX``/``dN``/``dE``/``dCO2`` are byte-for-byte identical with the byproducts
off. The integrated core therefore drifts only by the trace sugar these consume
(~0.2 % of ``S0``); the per-RHS carbon drawn from ``S`` exactly equals the carbon
deposited in the pool, so ``total_carbon`` closes to machine precision (verified by
the conservation tests). The carbon is split across sugar slots in proportion to each
slot's carbon content (:func:`_draw_carbon_from_sugar`), so the same code serves
wine's one slot and beer's three.

**The ``Byp`` double-count, resolved.** ``Byp`` formerly lumped "organic acids +
higher alcohols" (booked as succinic acid). Fusels *are* higher alcohols, so weighting
a separate ``fusels`` pool on top would book that carbon twice. Under D-19 ``Byp`` is
re-anchored to **organic acids / polyols only** (``Y_byproduct_sugar`` reduced to drop
the higher-alcohol share); the higher alcohols now live solely in the carbon-routed
``fusels`` pool. No overlap remains.

**Both carbon sources are bookkeeping stand-ins, not metabolic claims (D-19).**
(i) The Ehrlich pathway builds fusels from *amino-acid* skeletons, but ``N`` (YAN)
carries no carbon in :func:`total_carbon`, so fusel carbon is sourced from sugar.
(ii) An ester's ethanol moiety is carbon *already counted in ``E``*, so routing ester
carbon from sugar over-attributes fresh hexose. Both close the ledger exactly; neither
asserts where the carbon physically came from. The **sugar stand-in** carries **no CO2
co-product** (the Ehrlich decarboxylation is omitted) — a documented simplification that keeps
that draw a clean 1:1 sugar→pool carbon transfer.

**Since D-106 that omission is scoped to the stand-in, and the asymmetry is deliberate.**
:class:`FuselAminoAcidReroute` *does* emit the decarboxylation CO2, for the fraction it
re-sources onto real amino acids. The two are not inconsistent: a draw off an *amino acid* is a
chemical claim about a named reaction (leucine → isoamyl alcohol + CO2 + NH3), and it is wrong
unless it charges every one of leucine's six carbons — that was the D-105 finding. A draw off
*sugar* claims no reaction at all; it is a placeholder for de-novo synthesis whose real carbon
path (hexose → pyruvate → keto acid, decarboxylating on the way) the model does not trace. So
the CO2 appears exactly where a mechanism is asserted, and stays absent where one is not.

**The fusel sugar stand-in is now re-routable (decision D-33).** When the toggleable
``amino_acids`` pool (arginine; D-32) is dosed, :class:`FuselAminoAcidReroute` re-sources a
fraction of the fusel carbon off the sugar stand-in and onto that amino-acid pool — the
*physically-faithful* Ehrlich source — and **deaminates**, releasing the consumed amino
acids' nitrogen to the ammonium ``N`` pool. It is a separate wine-only *swap* Process (it
never touches ``fusels``; production stays here), sharing the one
:func:`fusel_production_rate` so its sugar refund matches this producer's draw exactly. That
deamination branch is the prerequisite the re-route was long deferred on (D-19/D-32); see
:class:`FuselAminoAcidReroute` for the closure algebra and the arginine N-over-release caveat.

**The gas-stripping sink (decisions D-20 → D-21).** The observed fall of wine *liquid*
ester with temperature is largely **evaporation** (Rollero 2014), not reduced synthesis.
That sink — logged as future work in D-19 — is built as :class:`EsterVolatilization`,
which strips each liquid ester into its own bookkeeping headspace pool on the
evolving-CO2 stream, with a **physical** Henry's-law partition (``dH_ester_volatil`` ≈
45 kJ/mol, sourced ethyl-acetate gas/liquid partition enthalpy; D-21 replaced D-20's
fudged per-medium ``E_a_ester_volatil``). Because a sourced stripping is medium-
independent, the wine/beer direction is carried by **per-medium sourced synthesis
``E_a_esters``** (beer steep / wine flat), not by the sink. Each transfer is carbon-neutral
(a pool and its twin share one molecule's weight), so ``total_carbon`` still closes to
machine precision. See that class for the full rationale, and for why v1 shares one sourced
ethyl-acetate partition enthalpy across all three esters (a documented approximation, D-96).

Tiers: :class:`EsterSynthesis` is **plausible** in form (warmth-favoured,
flux-coupled ester synthesis is the standard direction in the canonical *beer* model,
de Andrés-Toro 1998) with speculative rate parameters. *Wine caveat:* in wine the
warmer⇒more-ester direction is weaker and partly confounded — ester *synthesis* is
non-monotonic in temperature, so the warmer⇒more-aroma direction is carried by the
fusels while the **liquid** ester fall with warmth is the job of the volatilization sink
above (Mouret 2015; Rollero 2014). :class:`FuselAlcoholsEhrlich`
is **speculative** in form because its nitrogen dependence is knowingly simplified to a
single monotone branch (the real Ehrlich relationship is non-monotonic — handoff §3.2,
corroborated by Mouret 2015 / Rollero 2014: higher-alcohol synthesis is optimal at
~200–300 mg N/L). Parameter-tier propagation (D-1) caps the pool outputs at speculative
regardless. These earn no promotion past plausible: the benchmark is a *directional*
check (handoff §3.5).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from fermentation.core.chemistry import (
    CARBON_ATOMS,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics.amino_acid_pools import (
    SPEC_BY_SPECIES,
    AminoAcidSpec,
    depletion_gate,
    draw_precursor_carbon,
)
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import (
    ACETYL_CARBON_SHARE,
    ALCOHOL_CARBON_SHARE,
    CO2_PER_PRIMARY_EHRLICH_ALCOHOL,
    DE_NOVO_FUSEL_ROUTES,
    DE_NOVO_SHARE_BY_ALCOHOL,
    ESTER_SPECS,
    FUSEL_SPECS,
    ISOAMYL_ALCOHOL,
    LABELLED_PRECURSOR,
    SECONDARY_FUSEL_ROUTES,
    TRACER_BY_BULK,
    VALINE_LABEL_TRACERS,
    FuselSpec,
    labelled_fraction,
    non_ehrlich_fraction_param,
)
from fermentation.core.kinetics.carbon_routing import (
    draw_carbon_from_sugar as _draw_carbon_from_sugar,
)
from fermentation.core.kinetics.carbon_routing import (
    fermentative_flux_shape as _fermentative_flux_shape,
)
from fermentation.core.kinetics.carbon_routing import (
    refund_carbon_to_sugar as _refund_carbon_to_sugar,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: ``pool -> FuselSpec``, so the D-115 re-route can weight its debit at the precursor alcohol's
#: **own molecule**. Resolved through the registry rather than by assuming ``pool == species``:
#: that assumption happens to hold today, and a lookup that is right by coincidence is the
#: D-99 two-C5-isomers trap waiting to happen (``isoamyl_alcohol`` and ``active_amyl_alcohol``
#: share a formula and a carbon count, so a mis-resolution would be invisible to the ledger).
_FUSEL_SPEC_BY_POOL: dict[str, FuselSpec] = {spec.pool: spec for spec in FUSEL_SPECS}

#: Representative species whose formula carbon-accounts each aroma pool (D-19). The carbon
#: mass fraction of each weights both the sugar draw here and the pool in ``total_carbon``,
#: from the one chemistry source of truth — so the draw and the conservation check
#: can never disagree (cf. the ``Gly``/``Byp`` routing in the uptake Process).
#:
#: Neither the ester nor the fusel pools are lumped any more — each is its own molecule,
#: registered in ``ESTER_SPECS`` (decision D-96) / ``FUSEL_SPECS`` (decision D-99) and weighted
#: by itself, so there is no single ``_ESTER_SPECIES`` or ``_FUSEL_SPECIES`` to name here.
#: **No pool in the project is lumped any more (decision D-110).** ``mercaptans`` was the last —
#: and it turned out to be a FALSE lump rather than a real one, holding the single molecule it is
#: now named for (``methanethiol``). D-96 and D-99 split real mixtures; D-110 deleted a claim.

#: The Ehrlich decarboxylation releases exactly **1 mol CO₂ per mol higher alcohol** — the keto
#: acid's carboxyl carbon, the step that turns the transaminated amino acid into its aldehyde
#: (decision D-106). The twin of ``aging._CO2_PER_STRECKER_ALDEHYDE``, and for the same reason:
#: it is a real product on the carbon ledger, so charging it to the precursor's draw is what makes
#: that draw the molecule's actual stoichiometry.
#:
#: **"All five routes are C(precursor) == C(alcohol) + 1, so one constant covers the set" — that
#: note was retired at D-111**, which is exactly the ``+1`` it asserted. It was true of every
#: *primary* route and still is; it is false of the **secondary** valine → KIC → isoamyl route,
#: where the precursor is the *same* size as the alcohol (both C5), two carbons leave as CO₂, and
#: acetyl-CoA supplies the difference. Hence the alias: this is the PRIMARY route's count, and
#: :attr:`~fermentation.core.kinetics.carbon_routing.SecondaryFuselRoute.co2_per_alcohol` carries
#: the other. :func:`_branch` costs both from ``(n_p, n_a, k)`` alone, so neither is hand-written.
_CO2_PER_EHRLICH_ALCOHOL = float(CO2_PER_PRIMARY_EHRLICH_ALCOHOL)


def fusel_rate_shape(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """The Ehrlich rate shape SHARED by all five higher alcohols — everything but ``k`` [1/h].

    ``X·S_total/(K_sugar_uptake+S_total) · N/(K_n+N) · arrhenius(T, E_a_fusels, T_ref)``: the
    fermentative-flux Monod shape, gated on assimilable nitrogen, warmed by a
    steeper-than-uptake Arrhenius factor. Returns 0 under the guards
    :class:`FuselAlcoholsEhrlich` applies (no flux, no nitrogen).

    **That all five share this one shape is the honest limit of the D-99 split**, and is worth
    stating where it is implemented rather than only in the registry doc. The five *molecules*
    are now real — own carbon weight, own threshold, own independently-anchored ``k`` — but
    one N-gate and one ``E_a_fusels`` between them means the **spectrum is fixed**: warm the
    ferment or starve it of nitrogen and all five scale together, so no ratio among them can
    ever move. D-99 retires the wrong-molecule error, not the fixed-composition one; it
    downgrades it to a fixed *spectrum*. Per-species activation energies and per-amino-acid
    gates would make the composition dynamic and are deliberately NOT here: neither is sourced,
    and five author-estimated ``E_a``s would look like fidelity while adding only invention.
    """
    flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
    if flux <= 0.0:
        return 0.0
    n = max(float(y[schema.slice("N")][0]), 0.0)
    if n <= 0.0:  # Ehrlich needs assimilable nitrogen (amino acids)
        return 0.0
    nitrogen_gate = n / (params["K_n"] + n)
    temp = float(y[schema.slice("T")][0])
    f_t = arrhenius_factor(temp, params["E_a_fusels"], params["T_ref"])
    return float(flux * nitrogen_gate * f_t)


def fusel_production_rate(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], spec: FuselSpec
) -> float:
    """One higher alcohol's production rate ``d(<spec.pool>)/dt`` [g/L/h].

    ``k_<spec> · fusel_rate_shape(...)``. The per-species ``k`` is the ONLY thing that
    distinguishes the five (decision D-99), and each is anchored independently to its own
    molecule's measured concentration — never a share of a lumped ``k_fusel``.
    """
    shape = fusel_rate_shape(y, schema, params)
    if shape <= 0.0:
        return 0.0
    return float(params[spec.k_param] * shape)


def fusel_carbon_draw(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Total carbon [g C/L/h] the Ehrlich producer books out of sugar, across all five species.

    Each alcohol contributes ``rate_i · carbon_mass_fraction(species_i)`` — its own molecule's
    fraction, from the one chemistry source of truth (decision D-99). Before the split this was
    one rate at isoamyl alcohol's fraction standing in for all five.

    Factored out (decision D-33, generalised at D-99) as the single source of the fusel carbon
    draw so the *producer* (:class:`FuselAlcoholsEhrlich`, which draws it from sugar) and the
    *re-route* (:class:`FuselAminoAcidReroute`, which re-sources a fraction of it from the
    amino-acid pool and refunds the sugar) compute the **identical** number. Any divergence
    would break carbon closure, since the re-route refunds exactly what the producer drew —
    this helper makes that impossible (the shared ``biomass_growth_rate`` discipline of the
    D-32 swap, applied to fusels).
    """
    return float(sum(carbon for _, carbon in fusel_carbon_draw_by_species(y, schema, params)))


def fusel_carbon_draw_by_species(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float]
) -> list[tuple[FuselSpec, float]]:
    """Each higher alcohol's own carbon draw [g C/L/h], paired with its spec (decision D-100).

    The per-species decomposition of :func:`fusel_carbon_draw` (which now sums this, so the two
    can never disagree). D-33's re-route only ever needed the *total*, because it sourced every
    fusel's carbon from one lumped ``amino_acids`` pool. D-100 speciates that pool, so the
    re-route must know **which alcohol wants how much carbon** in order to debit that alcohol's
    own precursor — leucine for isoamyl alcohol, valine for isobutanol, and so on
    (``spec.precursor_amino_acid``, documentation-only since D-99, becomes load-bearing here).
    """
    shape = fusel_rate_shape(y, schema, params)
    if shape <= 0.0:
        return []
    return [
        (spec, params[spec.k_param] * shape * carbon_mass_fraction(spec.species))
        for spec in FUSEL_SPECS
    ]


def ehrlich_co2_carbon(spec: FuselSpec, alcohol_carbon: float) -> float:
    """The decarboxylation carbon [g C/L/h] charged to ``spec``'s precursor (decision D-106).

    ``alcohol_carbon`` is the **gated** carbon of the alcohol being re-sourced onto its amino acid.
    The Ehrlich pathway decarboxylates exactly once per alcohol, so the precursor must also give up
    one carbon beyond the alcohol's — and since ``alcohol_carbon`` is spread over
    ``CARBON_ATOMS[spec.species]`` carbons, that one carbon is simply the quotient. Charging it is
    what makes the draw 1 mol precursor per mol alcohol; without it the route consumes ``(n-1)/n``
    (the D-105 finding).

    **This exists as a shared helper because its absence is what let the two callers drift.**
    :class:`FuselAminoAcidReroute` draws it and
    :class:`~fermentation.core.kinetics.precursor_fates.PrecursorNonEhrlichFates` must scale its
    D-104 split against the *same* number — the sink's whole contract is that the realised split is
    exactly ``f : (1-f)`` of the precursor the re-route actually consumes. Before D-106 both
    recomputed ``gate x fusel_carbon`` independently and agreed by luck, because the two arithmetics
    were identical; adding the CO₂ to one silently moved the realised split in the other. One
    helper, two callers: they can no longer disagree (the D-33/D-99 shared-helper discipline, which
    this route had in name only).
    """
    return _CO2_PER_EHRLICH_ALCOHOL * alcohol_carbon / CARBON_ATOMS[spec.species]


@dataclass(frozen=True)
class EhrlichDraw:
    """One fully-costed Ehrlich sourcing branch: ``precursor → alcohol`` (decision D-111).

    **The shared helper the many-to-one map made necessary, and the D-106 discipline it enforces.**
    Before D-111 there was exactly one branch per alcohol, so the re-route and the D-104 sink could
    each recompute ``gate × fusel_carbon`` and agree. D-111 gives valine **two** branches
    (isobutanol via KIV, isoamyl alcohol via KIC), and the sink's contract is that consumed
    precursor splits exactly ``f : (1−f)`` between the non-Ehrlich lump and **everything Ehrlich** —
    so the sink must scale against the *sum* of a precursor's branches. Recomputing that in two
    places is the exact drift D-106 caught (two callers "computing the same thing" agreed by luck
    until one changed). One helper, two callers.

    Every field is carbon [g C/L/h], and **each branch is carbon-neutral by construction**:
    ``refund_carbon + co2_carbon − precursor_carbon == 0``. The alcohol itself is *not* credited
    here — the producer already made it and already drew sugar for it (D-109's constraint: the
    partition lives in the sourcing layer, never in the producer). A branch only moves *where that
    carbon came from*, which is what keeps
    :class:`FuselAlcoholsEhrlich` byte-for-byte when the precursor pools are undosed.
    """

    #: The amino-acid pool this branch debits.
    precursor: AminoAcidSpec
    #: The alcohol whose carbon is being re-sourced onto that precursor.
    alcohol: FuselSpec
    #: The alcohol carbon sourced on this branch — gated, so ≤ the alcohol's total draw.
    alcohol_carbon: float
    #: Carbon drawn out of the precursor pool: one whole molecule per alcohol.
    precursor_carbon: float
    #: Carbon leaving as CO₂ (one decarboxylation on the primary route, two via KIC).
    co2_carbon: float
    #: Carbon refunded to ``S``, undoing the producer's draw for the part not actually from sugar.
    refund_carbon: float

    @property
    def nitrogen(self) -> float:
        """Nitrogen [g N/L/h] deaminated out of the precursor — what the ``N`` pool gains."""
        return (
            self.precursor_carbon
            / carbon_mass_fraction(self.precursor.species)
            * nitrogen_mass_fraction(self.precursor.species)
        )


def _branch(
    precursor: AminoAcidSpec, alcohol: FuselSpec, alcohol_carbon: float, co2_per_alcohol: int
) -> EhrlichDraw:
    """Cost one branch from its stoichiometry alone — the general form of the D-106 charge.

    For a precursor of ``n_p`` carbons becoming an alcohol of ``n_a`` carbons while releasing
    ``k`` CO₂, per mole of alcohol: the precursor gives up ``n_p``, ``k`` carbons leave as CO₂,
    and the sugar refund is ``n_p − k`` — the producer drew ``n_a`` for the alcohol, the truth
    needs ``n_a + k − n_p`` from sugar (as acetyl-CoA), and the difference is what is handed back.

    **The refund is provenance-independent** (see :class:`SecondaryFuselRoute`): it falls out of
    the net sugar balance, not out of which carbons the decarboxylations remove — so D-109's
    ``{3,4,5}/5`` atom-assignment trap never reaches the ledger.

    **This reduces to the primary route exactly**, which is the check that D-111 moved nothing it
    did not mean to: there ``n_p = n_a + 1`` and ``k = 1``, giving ``precursor = alcohol·(n_a+1)/n_a
    = alcohol + ehrlich_co2_carbon``, ``co2 = alcohol/n_a`` and ``refund = alcohol·n_a/n_a =
    alcohol`` — the three numbers the pre-D-111 code wrote by hand.
    """
    n_a = CARBON_ATOMS[alcohol.species]
    n_p = CARBON_ATOMS[precursor.species]
    per_alcohol_carbon = alcohol_carbon / n_a
    return EhrlichDraw(
        precursor=precursor,
        alcohol=alcohol,
        alcohol_carbon=alcohol_carbon,
        precursor_carbon=per_alcohol_carbon * n_p,
        co2_carbon=per_alcohol_carbon * co2_per_alcohol,
        refund_carbon=per_alcohol_carbon * (n_p - co2_per_alcohol),
    )


def ehrlich_draws(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float]
) -> list[EhrlichDraw]:
    """Every Ehrlich sourcing branch, primary and secondary, fully costed (decision D-111).

    The single source of truth for *where each higher alcohol's carbon came from*, shared by
    :class:`FuselAminoAcidReroute` (which applies it) and
    :class:`~fermentation.core.kinetics.precursor_fates.PrecursorNonEhrlichFates` (which scales
    the non-Ehrlich lump against it). See :class:`EhrlichDraw` for why one helper is mandatory.

    **Primary branches** are each alcohol's own ``precursor_amino_acid``, gated on that pool's
    relative depletion (D-100) — unchanged from D-100/D-106.

    **Secondary branches**
    (:data:`~fermentation.core.kinetics.carbon_routing.SECONDARY_FUSEL_ROUTES`)
    are anchored *algebraically* against their precursor's primary branch rather than gated
    independently, because the sourced quantity is a **share of consumed precursor**: valine's
    consumption splits ``0.15 : 0.23 : 0.62`` between isobutanol, isoamyl alcohol and everything
    else, and the isobutanol branch is the one the model anchors (via ``k_isobutanol``). So
    ``consumed = primary_carbon / share_primary`` and the secondary branch takes
    ``share_secondary × consumed``. That is D-109's "algebraic partition in the sourcing layer":
    no new state, no second gate, and the primary gate already throttles the whole node as the
    precursor empties (so an undosed run is still a no-op, and the pool cannot go negative).
    """
    draws: list[EhrlichDraw] = []
    by_alcohol: dict[str, float] = {}
    primary_by_precursor: dict[str, EhrlichDraw] = {}
    for spec, fusel_carbon in fusel_carbon_draw_by_species(y, schema, params):
        if fusel_carbon <= 0.0:
            continue
        precursor = SPEC_BY_SPECIES[spec.precursor_amino_acid]
        gate = depletion_gate(y, schema, params, (precursor,))
        if gate <= 0.0:
            continue  # exhausted ⇒ this alcohol stays wholly on the sugar stand-in
        # D-118: an alcohol built mostly from central carbon metabolism may only be sourced from
        # its amino acid up to (1 − f_de_novo). The gate models AVAILABILITY; this models
        # PROVENANCE, and for 2-phenylethanol they differ by ~11× (see DeNovoFuselRoute). Applied
        # MULTIPLICATIVELY rather than as a min() cap so the gate's dynamics survive underneath —
        # the share scales what the gate would otherwise attribute, and an undosed run is still
        # byte-for-byte (gate 0 ⇒ 0, whatever the share).
        de_novo_param = DE_NOVO_SHARE_BY_ALCOHOL.get(spec.pool)
        if de_novo_param is not None:
            de_novo = params[de_novo_param]
            if not 0.0 <= de_novo < 1.0:
                raise ValueError(
                    f"{de_novo_param}={de_novo} outside [0, 1): it is the fraction of "
                    f"{spec.pool} built de novo rather than from consumed "
                    f"{spec.precursor_amino_acid}"
                )
            gate *= 1.0 - de_novo
            if gate <= 0.0:
                continue
        branch = _branch(precursor, spec, gate * fusel_carbon, CO2_PER_PRIMARY_EHRLICH_ALCOHOL)
        draws.append(branch)
        primary_by_precursor[spec.precursor_amino_acid] = branch
        by_alcohol[spec.pool] = branch.alcohol_carbon

    totals = dict(fusel_carbon_draw_by_species(y, schema, params))
    for route in SECONDARY_FUSEL_ROUTES:
        primary = primary_by_precursor.get(route.precursor)
        if primary is None:
            continue  # the precursor is exhausted (or makes no alcohol) ⇒ no anchor, no route
        f = params[non_ehrlich_fraction_param(route.precursor)]
        share = params[route.share_param]
        share_primary = 1.0 - f - share
        if share_primary <= 0.0:
            raise ValueError(
                f"{non_ehrlich_fraction_param(route.precursor)}={f} + {route.share_param}={share} "
                f"leaves no share for {primary.alcohol.pool}: the three fates of consumed "
                f"{route.precursor} must sum to 1 with a POSITIVE primary branch"
            )
        alcohol = next(s for s in FUSEL_SPECS if s.pool == route.alcohol_pool)
        precursor = SPEC_BY_SPECIES[route.precursor]
        consumed = primary.precursor_carbon / share_primary
        n_a, n_p = CARBON_ATOMS[alcohol.species], CARBON_ATOMS[precursor.species]
        alcohol_carbon = share * consumed * n_a / n_p
        # This alcohol cannot be sourced from more precursor than it is being made from. The
        # secondary branch is anchored off a DIFFERENT alcohol's draw, so nothing about its
        # arithmetic bounds it by this alcohol's own production — unlike a gated primary branch,
        # which is a fraction of exactly that. Measured headroom is large (valine ~1.8% + leucine
        # ~1.1% of isoamyl), so this never binds in practice; it is here because if it ever DID
        # bind, the refund would hand back more sugar than the producer drew for this alcohol and
        # the ledger would still close — the D-89/D-90 denominator-trap family, where conservation
        # is blind and only an explicit guard is not.
        headroom = totals.get(alcohol, 0.0) - by_alcohol.get(alcohol.pool, 0.0)
        if alcohol_carbon > headroom:
            alcohol_carbon = max(0.0, headroom)
        if alcohol_carbon <= 0.0:
            continue
        draws.append(_branch(precursor, alcohol, alcohol_carbon, route.co2_per_alcohol))
    return draws


class EsterSynthesis(Process):
    """Ester production, coupled to the fermentative flux and favoured by warmth.

    For each ester in :data:`~fermentation.core.kinetics.carbon_routing.ESTER_SPECS`,
    ``d(ester)/dt = k_<ester> · X · S_total/(K_sugar_uptake + S_total) · f(T)`` with
    ``f(T) = arrhenius_factor(T, E_a_esters, T_ref)``. Esters form alongside fermentation;
    tying synthesis to the biomass-catalysed sugar flux (sharing ``K_sugar_uptake``) couples
    them to that flux directly. ``E_a_esters`` is **sourced per medium** (decision D-21):
    **steep** for beer (de Andrés-Toro 1998, ester ride growth — synthesis outruns the
    stripping sink, so beer liquid esters *rise* with T) and **weak/~flat** for wine (Mouret
    2015 / Rollero 2014, wine ester synthesis is weak and non-monotonic — so the gas-stripping
    sink :class:`EsterVolatilization` wins and wine liquid esters *fall* with T). This is only
    the *synthesis* term; net liquid ester is synthesis minus stripping (see that class and
    the module docstring). Each ester's carbon is routed *out of ``S``* via
    :func:`_draw_carbon_from_sugar` **at its own molecule's carbon fraction** (option a1,
    D-19), so this touches the three ester pools and ``S`` — never ``E``/``CO2`` — and
    ``total_carbon`` closes exactly. See the module docstring for the ester carbon-source
    caveat.

    **Three esters, three independently-sourced rates (decision D-96).** Before D-96 a single
    ``k_ester`` filled one lumped pool. The pool is now split into three single-molecule pools
    and — the load-bearing call — each ``k`` is anchored to **its own molecule's** measured
    concentration range rather than to a share of the old total. Splitting one ``k_ester`` by a
    fitted ratio would have reproduced the same numbers while smuggling back the very
    fabricated-composition constant the split exists to remove; here the lump's composition and
    the total ester mass are *derived*, not targeted.

    **Shared temperature shape — a documented simplification (D-96).** All three read the one
    per-medium ``E_a_esters``. For the two acetates this is principled: same enzyme (ATF1), so
    no basis to differ. For ``ethyl_hexanoate`` (EEB1/EHT1 — a different enzyme, tied to
    fatty-acid synthesis) it is a genuine approximation, kept because no separate activation
    energy is sourced; a per-ester ``E_a`` is the named deferred refinement.

    **The ATF1 precursor coupling (decision D-97) — why isoamyl acetate alone carries an extra
    term.** ATF1 acetylates an *alcohol*, and Fujii 1998 (Appl. Environ. Microbiol.
    64:4076-4078) measures its ``Km`` for isoamyl alcohol at ~29.8 mM while stating outright
    that *"a major rate-limiting factor for isoamyl acetate production is the amount of isoamyl
    alcohol"*. The ``fusels`` pool runs ~0.5-1 mM — **~30-60x below that Km** — so the enzyme
    sits far down its linear stretch and the rate is **first-order in the pool**:
    ``d(isoamyl_acetate)/dt = k · fusels · X · S/(K_su+S) · f(T)``. Ethyl acetate gets no such
    term because *its* precursor is ethanol (~2 M, orders of magnitude above any mM-scale Km)
    ⇒ ATF1 is saturated in it ⇒ zeroth-order. Same enzyme, same rate law, opposite limits: the
    asymmetry is **derived from the two precursors' concentrations**, not an exemption.

    **Why first-order and not an explicit ``fusels/(Km+fusels)`` gate — identifiability, not
    parsimony (D-97).** In the ``[S] << Km`` limit the saturable form is numerically the linear
    one to within ~2 %, and only the *ratio* ``Vmax/Km`` is identifiable from any model output:
    scaling ``Km`` tenfold and refitting ``Vmax`` would give a byte-identical trajectory. Adding
    the measured ``Km`` as a parameter would therefore mint a sourced-looking constant that no
    model output could ever validate. ``k_isoamyl_acetate`` **is** that identifiable ratio, and
    the ``Km`` lives where it does its real work: in the parameter's provenance, as the sourced
    justification for this rate law's *form*.

    **What this buys — the banana note becomes YAN-responsive.** Reading the ``fusels`` pool
    makes isoamyl acetate inherit the nitrogen dependence of the Ehrlich pathway that builds its
    precursor: a low-YAN ferment makes less isoamyl alcohol, hence less banana. Before D-97 the
    ester was **YAN-blind** (flat ~0.758 mg/L across YAN 40-250 mg/L, where ``fusels`` swung
    2.9x) because every ester shared one plain flux shape. Note this couples to the *pool*, not
    to the fusel *production rate*: the alcohol **persists** after nitrogen is exhausted (in a
    wine run ``N`` empties around day 2 with ~75 % of the sugar still unfermented), and ATF1 goes
    on acetylating it for as long as the flux supplies acetyl-CoA. Coupling to the rate would
    have stopped the banana dead at day 2 with 51 mg/L of its substrate sitting in the vessel.

    **Read AND debited since D-115 — the 5:2-inverse re-route, built.** D-97 read the alcohol
    pool without debiting it and drew the whole ester from ``S``, naming the carbon re-route as
    a deferred refinement on the grounds that it is **mass-negligible** (~0.5 mg/L of ester
    against an ~86 mg/L alcohol pool). That reasoning was correct and remains correct; what
    D-114 established is that it does not bound what the re-route is *for*. Mass-negligible and
    **observable**-negligible are different quantities: the deferral was 100 % of the ester's
    valine enrichment, which was structurally pinned at zero for as long as every gram of ester
    carbon came from sugar. So the acetylation now takes the C5 skeleton *from* the alcohol and
    only the C2 acetyl group from ``S``, at
    :data:`~fermentation.core.kinetics.carbon_routing.ALCOHOL_CARBON_SHARE` — the **exact
    inverse** of the 5:2 split D-69's hydrolysis returns, reading the one shared ratio so the
    two directions of one reaction cannot drift apart.

    **The debit cannot empty the pool, and that is structural rather than guarded.** The rate is
    itself first-order in this pool (D-97), so the draw decays to zero exactly as the pool does.
    A clamp would be worse than useless here: it would silently break the ester's own mass
    balance rather than prevent anything that can happen.

    **What it buys — the ester carries label (D-115).** Because the C5 arrives as a unit from a
    pool with a known valine-derived fraction, an ester molecule is labelled exactly when its
    parent alcohol was, and the enrichment Rollero 2017 measures becomes a model output instead
    of a structural zero. Both tracer slots are written here: the alcohol's is debited and the
    ester's credited at the alcohol pool's *current* fraction. See
    :class:`~fermentation.core.kinetics.carbon_routing.LabelTracer` for why the ester needs a
    slot of its own rather than inheriting the alcohol's number at the analysis boundary.

    **A stated dependency, not a latent surprise (D-97).** Reading ``fusels`` means
    ``isoamyl_acetate`` synthesis **requires** :class:`FuselAlcoholsEhrlich` to be active — the
    two are co-wired in ``media._BYPRODUCT_PROCESSES`` (the only tuple either appears in), so
    every supported ProcessSet has both. A ProcessSet built with this Process alone is not
    broken, it is merely *precursor-free*: the banana rate is identically 0 (no alcohol to
    acetylate) while the other two esters form normally. Note the tier consequence is benign
    and needs no machinery: a **plausible**-form Process now reads a pool a **speculative**
    one fills, but ``k_isoamyl_acetate`` is speculative anyway, so the pool's *output* tier is
    already at the speculative floor and understates nothing. (Tier propagation runs through
    parameters and ``touches``, not through state dependencies — a pre-existing property of the
    tier system, not something this coupling introduces.)

    **Ensemble consequence (D-97).** Because the rate multiplies ``[fusels]``,
    ``isoamyl_acetate`` now **inherits the fusel pool's variance** (via ``k_fusel`` / ``K_n`` /
    ``E_a_fusels``) under a D-24 ensemble sweep. That is *more* correct — the precursor's
    uncertainty genuinely is the ester's — but it means the band widens beyond what
    ``k_isoamyl_acetate``'s own uncertainty implies. The parameter's "pins the same
    concentration span" note is a statement about the **point estimate** holding ``fusels``
    fixed; it is not a claim about the ensemble band.

    **Two inherited caveats, named not buried (D-97).** (i) ``fusels`` is a **lumped** pool
    (isoamyl alcohol representative, but really all the higher alcohols), so reading it whole as
    the ATF1 substrate over-states the true isoamyl-alcohol supply. The over-statement is
    absorbed into the re-anchored ``k``; what it costs is that the *YAN-response* assumes the
    lump's composition is fixed — the D-66 caveat, honest here because ``fusels`` is genuinely
    flagged lumped. (ii) Isoamyl acetate now inherits :class:`FuselAlcoholsEhrlich`'s
    **speculative monotone-in-N** shape, including its admitted omission of the low-YAN
    biosynthetic rise. Both sit inside this Process's existing plausible/speculative framing.
    """

    name = "ester_synthesis"
    tier = Tier.PLAUSIBLE
    #: Since D-115 this also **debits** the precursor alcohol (the 5:2-inverse re-route) and
    #: writes the two valine-label tracer slots, so all three join the three ester pools and
    #: ``S``. The D-97-era ``touches`` deliberately excluded ``isoamyl_alcohol`` — "read, never
    #: debited" — and a test asserted that exclusion; D-115 retires both, which is the
    #: *point* of the beat rather than a regression (see the class doc).
    touches = (
        *(spec.pool for spec in ESTER_SPECS),
        "S",
        ISOAMYL_ALCOHOL.pool,
        *(tracer.tracer_pool for tracer in VALINE_LABEL_TRACERS),
    )
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks;
    #: ``E_a_esters`` (sourced per medium — steep beer / flat wine, D-21) and ``T_ref``
    #: set the temperature shape, shared by all three esters (D-96). Each ester has its own
    #: independently-sourced rate constant. Their tiers cap the ester pools' output tier via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        *(spec.k_param for spec in ESTER_SPECS),
        "K_sugar_uptake",
        "E_a_esters",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_esters"], params["T_ref"])
        for spec in ESTER_SPECS:
            rate = params[spec.k_param] * flux * f_t
            if spec.precursor_pool is not None:
                # ATF1 is FAR from saturated in this alcohol ([S] ~ 0.5-1 mM vs Km ~29.8 mM),
                # so its supply limits the rate first-order (D-97).
                # Clamped >= 0 so a solver undershoot cannot flip synthesis negative.
                rate *= max(float(y[schema.slice(spec.precursor_pool)][0]), 0.0)
            d[schema.slice(spec.pool)] = rate
            # Each ester draws at ITS OWN carbon fraction (C4/C7/C8) — the D-96 split's
            # ledger payoff: no ester's carbon is booked through a stand-in molecule.
            ester_carbon = rate * carbon_mass_fraction(spec.species)
            if spec.precursor_pool is None:
                _draw_carbon_from_sugar(d, y, schema, ester_carbon)
                continue

            # --- the D-115 re-route: the C5 comes OFF the alcohol, only C2 off sugar ---------
            # The exact inverse of the D-69 hydrolysis split, using the same shared ratio. Two
            # sources, split 5:2 by carbon, which is mole-for-mole with the ester: the alcohol
            # gives up one whole molecule per ester molecule made.
            precursor = _FUSEL_SPEC_BY_POOL[spec.precursor_pool]
            alcohol_carbon = ester_carbon * ALCOHOL_CARBON_SHARE
            alcohol_rate = alcohol_carbon / carbon_mass_fraction(precursor.species)
            d[schema.slice(spec.precursor_pool)] -= alcohol_rate
            _draw_carbon_from_sugar(d, y, schema, ester_carbon * ACETYL_CARBON_SHARE)
            # The debit cannot drive the pool negative: `rate` is itself FIRST-ORDER in this
            # pool (D-97), so the draw decays to zero exactly as the pool does — first-order
            # decay, not a fixed subtraction. No clamp is needed and none is used, because a
            # clamp here would silently break the ester's own mass balance instead.

            # --- and the label rides across with it (D-115) ---------------------------------
            # Rollero's enrichment is a MOLECULE fraction and the C5 transfers as a unit, so an
            # ester molecule is labelled exactly when its parent alcohol was: both the alcohol
            # debit and the ester credit go at the alcohol pool's CURRENT fraction. This is the
            # only flow that moves label from the alcohol to the ester.
            f_alcohol = labelled_fraction(y, schema, TRACER_BY_BULK[spec.precursor_pool])
            if f_alcohol > 0.0:
                d[schema.slice(TRACER_BY_BULK[spec.precursor_pool].tracer_pool)] -= (
                    alcohol_rate * f_alcohol
                )
                d[schema.slice(TRACER_BY_BULK[spec.pool].tracer_pool)] += rate * f_alcohol
        return d


class FuselAlcoholsEhrlich(Process):
    """Fusel (higher) alcohols via the Ehrlich pathway — amino-acid-derived, warm.

    ``d(fusels)/dt = k_fusel · X · S_total/(K_sugar_uptake + S_total) ·
    N/(K_n + N) · f(T)`` with ``f(T) = arrhenius_factor(T, E_a_fusels, T_ref)``. The
    Ehrlich pathway makes higher alcohols by transamination/decarboxylation of
    amino acids, so production is gated on yeast-assimilable nitrogen availability
    (``N/(K_n + N)``, sharing the growth nitrogen half-saturation ``K_n``) on top of
    the fermentative flux — front-loading fusels into the nitrogen-replete early
    ferment, as observed. Steeper-than-uptake ``E_a_fusels`` gives the
    warmer-is-more direction. The fusel carbon (booked as isoamyl alcohol) is routed
    *out of ``S``* via :func:`_draw_carbon_from_sugar` (option a1, D-19), so it
    touches ``fusels`` and ``S`` — never ``E``/``CO2`` — and ``total_carbon`` (which
    now weights ``fusels``) closes exactly. The sugar source is a bookkeeping
    stand-in for the amino-acid skeleton (``N`` carries no carbon); see the module
    docstring.

    **Known simplification — monotone in nitrogen (why this Process is speculative).**
    The real fusel/nitrogen relationship is *non-monotonic*: higher alcohols rise
    again at very low YAN (the anabolic/biosynthetic route from sugar) as well as
    with amino-acid catabolism at higher YAN. v1 models only the catabolic,
    monotone-increasing-in-N branch and omits the low-N biosynthetic rise. Per the
    handoff (§3.2) this *models the pathway* (Ehrlich, N-gated) rather than fitting a
    slope, but the single-branch shape is the simplification that keeps the Process
    at the speculative tier until the non-monotonic form is sourced.
    """

    name = "fusel_alcohols_ehrlich"
    tier = Tier.SPECULATIVE
    #: The five single-molecule pools of ``FUSEL_SPECS`` plus ``S`` (decision D-99), derived
    #: from the registry so a sixth alcohol cannot silently violate the `touches` contract.
    touches: tuple[str, ...] = (*(spec.pool for spec in FUSEL_SPECS), "S")
    #: ``K_sugar_uptake``/``K_n`` are shared with the uptake/growth Processes;
    #: ``E_a_fusels`` (> ``E_a_uptake``) and ``T_ref`` set the temperature shape — one shared
    #: shape for all five (see :func:`fusel_rate_shape` for why that is the honest limit).
    #: The five ``k``s are per-species and independently anchored (D-99).
    reads: tuple[str, ...] = (
        *(spec.k_param for spec in FUSEL_SPECS),
        "K_sugar_uptake",
        "K_n",
        "E_a_fusels",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        shape = fusel_rate_shape(y, schema, params)
        if shape <= 0.0:
            return d
        for spec in FUSEL_SPECS:
            d[schema.slice(spec.pool)] = params[spec.k_param] * shape
        # ONE draw for all five, each at its own molecule's carbon fraction — shared verbatim
        # with the re-route (D-33/D-99) so the draw and its refund can never disagree.
        _draw_carbon_from_sugar(d, y, schema, fusel_carbon_draw(y, schema, params))
        return d


class FuselAminoAcidReroute(Process):
    """Ehrlich fusel carbon re-sourced from amino acids, with deamination (decision D-33).

    The physically-faithful other half of :class:`FuselAlcoholsEhrlich`. That producer books
    fusel carbon out of **sugar** — a documented bookkeeping stand-in, because the real Ehrlich
    pathway builds higher alcohols from *amino-acid* skeletons (transamination →
    decarboxylation → reduction), releasing the amino group as ammonium. Sugar is used only
    because ``N`` (YAN) carries no carbon in :func:`total_carbon` (D-19), so there was nowhere
    else to draw it from — *until* the toggleable ``amino_acids`` pool (arginine; D-32) gave
    the model a carbon- *and* nitrogen-bearing amino-acid source. This Process closes that gap:
    when amino acids are present it **re-routes** a fraction of the fusel carbon off the sugar
    stand-in and onto the amino-acid pool, and **deaminates** — releasing the nitrogen the
    consumed amino acids carried back to the ammonium ``N`` pool. That deamination branch is
    exactly the prerequisite the fusel re-route was deferred on (D-19/D-32); it also demonstrates
    the aa→N release path a later MLF-with-growth model needs.

    **A swap, not a producer — it never touches ``fusels``.** Like the D-32
    :class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation` swap, this leaves the
    *production* entirely to :class:`FuselAlcoholsEhrlich` (both call the one
    :func:`fusel_production_rate`); it only moves the carbon *source*. **Per alcohol** ``i``
    (decision D-100), for its own carbon draw ``F_i = rate_i·c_i`` and its own precursor's
    relative-depletion gate ``g_i = aa_i/(K_amino_acids·f_i + aa_i)``:

      * **refund sugar** by ``g_i·F_i`` (undoing the producer's draw for that fraction),
      * **debit that alcohol's precursor** by ``g_i·F_i / c_precursor``, and
      * **release ammonium** ``N`` by the nitrogen that mass carried (deamination).

    Carbon closes: the fusel gains ``F_i`` (from the producer), sourced now as ``(1−g_i)·F_i``
    from sugar + ``g_i·F_i`` from its precursor. Nitrogen closes: the precursors lose exactly the
    nitrogen the ``N`` pool gains. Net sugar is ``−Σ(1−g_i)·F_i ≤ 0`` for all ``g_i ≤ 1`` — the
    re-route never creates sugar (it only *spares* it), so the ABV bookkeeping caveat is the D-32
    one (spared sugar ferments to ethanol). **Wine-only** and **forced to be a separate Process**:
    declaring the precursor pools/``N`` in the both-media producer's ``touches`` would break
    beer's ProcessSet construction (beer tracks no amino acids).

    **The D-100 decoupling — this Process no longer touches ``amino_acids``.** Until D-100 it drew
    every alcohol's carbon from the lumped pool, i.e. **from arginine**, which does not make higher
    alcohols. That was not merely imprecise: at D-99's honest ~3.8× fusel rise the re-route drained
    the shared lump to ~0 and starved Maillard, MLF growth and Brett growth — three unrelated
    subsystems broken through one lumped substrate. Now each alcohol eats *its own* precursor and
    the identity-agnostic pools ({arginine, generic} — 81% of the must spectrum) are untouchable by
    fusel production, so that cross-subsystem starvation is structurally impossible rather than
    tuned away.

    **The anabolic/catabolic split became EMERGENT (the D-100 fidelity gain).** Real must carries
    ~30-60 mg/L leucine but wine makes ~150-250 mg/L isoamyl alcohol — most higher alcohol is
    synthesised **de novo from sugar**, not catabolised from amino acids. The lumped model could
    not represent that: it re-sourced a fixed gate-fraction of fusel carbon from a large arginine
    pool indefinitely. Now leucine's own gate throttles its re-route as leucine depletes and the
    remainder stays on the sugar stand-in — so the anabolic/catabolic ratio *falls out of* the must
    spectrum and the fusel demand instead of being a fitted fraction. The sugar stand-in it leaves
    behind is no longer an embarrassment but the **correct** book for de-novo synthesis.

    **[D-118 — TRUE IN SHAPE, WRONG BY AN ORDER OF MAGNITUDE FOR ONE MOLECULE. Read the paragraph
    above with this caveat attached.]** The gate does make the split emergent, but "emergent" is not
    "calibrated": measured at a wine-like must, this Process sourced **18.9%** of its
    2-phenylethanol from consumed phenylalanine against a derived **~1.7%** — an ~11×
    over-attribution the paragraph above would have you believe was already handled. **The gate
    encodes availability, not provenance.** Phenylalanine is available; it is simply not what most
    2-phenylethanol is made from, and a gate on the precursor's own concentration cannot express
    that — least of all here, where the must carries *fewer moles of phenylalanine than the wine
    makes of the alcohol*, so full sourcing is stoichiometrically impossible rather than merely
    generous. :data:`~fermentation.core.kinetics.carbon_routing.DE_NOVO_FUSEL_ROUTES` caps the
    branch for exactly the alcohols where the two quantities diverge. Only 2-phenylethanol is
    listed; **the other four are unmeasured, not verified innocent**, and isoamyl alcohol is the
    live suspect — D-104's inverted leucine split and D-113's "un-inverting leucine remains an
    unsourced build (de-novo-KIC relief …)" describe this same gap one precursor over, and this
    registry is the structure that would carry that fix once a leucine number is sourced.

    **The D-33 nitrogen over-release lump is RETIRED (decision D-100).** D-33 had to document that
    sourcing fusel carbon through N-rich arginine deaminated ``c_fusel/c_aa·y_N`` ≈ 0.78 g N per g
    fusel carbon — roughly **4× the real leucine→isoamyl-alcohol N:C** (leucine carries one amino
    group over six carbons). Drawing each alcohol's actual precursor releases exactly the nitrogen
    that molecule carries, so the ratio is now right by construction and the caveat is gone. Tier
    **speculative** (it inherits the fusel rate's speculative parameters and the spectrum/gate
    estimates).

    **Isolability (undosed-only, paired with the producer).** Every precursor gate → 0 at
    ``aa_i = 0``, so an undosed wine run is byte-for-byte the sugar-stand-in producer; the compile
    seam additionally *disables* this Process when ``amino_acids_gpl ≤ 0`` (tier isolability, the
    D-32 pattern). It is only valid while :class:`FuselAlcoholsEhrlich` is active — it refunds
    sugar that producer drew — so the two are kept paired (disabling the producer alone would let
    the re-route create sugar; the same acceptable swap↔producer coupling as D-32's swap↔growth).

    **The honest limit D-100 does NOT fix.** Speciation does not end the precursor competition — it
    *localises* it. Leucine feeds both isoamyl alcohol (here) and 3-methylbutanal (D-87); the same
    holds for valine/isoleucine/phenylalanine and their thermal aldehydes. Those are the same
    molecule in reality, so the competition is **real chemistry, not an artifact**, and the model
    should show it. What D-100 removes is the *false* competition — fusels vs bacterial growth over
    arginine, which shares no chemistry at all.

    **The threonine/sotolon example this paragraph used to give is RETIRED (D-107, measured D-109)
    — and it was wrong in BOTH halves.** It read "threonine still feeds both propanol (here) and
    sotolon (D-87)". Sotolon left D-87 at D-107 (it is not a Strecker product at all, but an aldol
    — :class:`~fermentation.core.kinetics.aging.SotolonAldolCondensation`), and its substrate is
    **α-ketobutyrate**, threonine's *child*. Threonine reaches sotolon only by supplying the
    **carbon source** of :class:`~fermentation.core.kinetics.keto_acids.AlphaKetobutyrateExcretion`,
    whose *rate* is flux-only — so draining threonine here does not cost that pool one microgram,
    and propanol and sotolon do **not** compete over threonine. D-109 measured it exactly (the two
    derivatives are bit-identical at threonine 67 mg/L and at 0). The real propanol-vs-sotolon
    competition is over α-ketobutyrate, is **not** expressible from this Process, and is the fusel
    side of the keto-acid node (D-109's scoping).
    """

    name = "fusel_amino_acid_reroute"
    tier = Tier.SPECULATIVE
    #: Refunds carbon to ``S``, debits **each alcohol's own precursor**, releases nitrogen to
    #: ``N``, and emits the decarboxylation ``CO2`` (D-106 — the term whose absence made the draw
    #: ``(n-1)/n``). Never touches ``fusels`` — production stays entirely in
    #: :class:`FuselAlcoholsEhrlich`
    #: — and, since D-100, **never touches ``amino_acids``/``amino_acids_generic``**: arginine does
    #: not make higher alcohols, so the re-route can no longer starve the consumers that live on
    #: those pools (yeast swap, MLF growth, Brett growth, Maillard browning). That absence is the
    #: whole D-100 decoupling and is pinned by a test.
    #: Since D-115 it also credits the **alcohol** label tracer: this is where a gram of isoamyl
    #: alcohol first becomes valine-derived, so it is where the label enters the model.
    touches = (
        "S",
        "N",
        "CO2",
        *(spec.precursor_amino_acid for spec in FUSEL_SPECS),
        VALINE_LABEL_TRACERS[0].tracer_pool,
    )
    #: Recomputes the fusel rate (so it reads the producer's parameters) plus ``K_amino_acids``
    #: and each precursor's ``must_aa_fraction_*`` for the relative-depletion gates (D-100). Since
    #: D-111 it also reads each secondary route's share **and its precursor's non-Ehrlich
    #: fraction** — not because this Process books the lump (it does not; that is the D-104 sink's
    #: job) but because the two together fix the *primary* branch's share as the residue
    #: ``1 − f_non_ehrlich − f_secondary``, which is what the secondary draw is anchored against.
    #: Their tiers cap the ``S``/precursor/``N`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        *(spec.k_param for spec in FUSEL_SPECS),
        "K_sugar_uptake",
        "K_n",
        "E_a_fusels",
        "T_ref",
        "K_amino_acids",
        *(SPEC_BY_SPECIES[spec.precursor_amino_acid].fraction_param for spec in FUSEL_SPECS),
        *(route.share_param for route in SECONDARY_FUSEL_ROUTES),
        *(non_ehrlich_fraction_param(route.precursor) for route in SECONDARY_FUSEL_ROUTES),
        # D-118: the de-novo share caps each listed alcohol's primary branch.
        *(route.share_param for route in DE_NOVO_FUSEL_ROUTES),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # EVERY sourcing branch, primary and secondary, costed from its own stoichiometry — the
        # one helper the D-104 sink also reads, so the split it imposes is against exactly the
        # precursor this Process consumes (D-106's shared-helper discipline, D-111's many-to-one).
        draws = ehrlich_draws(y, schema, params)
        if not draws:
            return d
        refund = 0.0  # total carbon re-sourced off sugar, across every branch
        nitrogen = 0.0  # total nitrogen deaminated out of the precursors
        co2_carbon = 0.0  # total decarboxylation carbon, precursor → CO2 (D-106; ×2 via KIC)
        for draw in draws:
            # Debit the precursor a whole molecule per alcohol. The gate that produced this
            # branch → 0 as the pool empties (D-100's relative-depletion rule), so the draw can
            # never drive it negative and an undosed run is a no-op.
            nitrogen += draw_precursor_carbon(
                d, schema, draw.precursor.species, draw.precursor_carbon
            )
            refund += draw.refund_carbon
            co2_carbon += draw.co2_carbon
            # D-115: this branch's alcohol carries VALINE's label iff valine is what it was
            # sourced from. Credited as grams of the labelled MOLECULE (Rollero's enrichment is
            # a molecule fraction), and only for the labelled precursor — the leucine → isoamyl
            # primary branch feeds the same pool and must NOT be credited here.
            tracer = TRACER_BY_BULK.get(draw.alcohol.pool)
            if tracer is not None and draw.precursor.species == LABELLED_PRECURSOR:
                d[schema.slice(tracer.tracer_pool)] += draw.alcohol_carbon / carbon_mass_fraction(
                    draw.alcohol.species
                )
        if refund <= 0.0:
            return d
        d[schema.slice("N")] = nitrogen  # DEAMINATION: precursor N → ammonium (the D-33 branch)
        # Precursor carbon that became CO₂, not alcohol: drawn above, emitted here, never refunded.
        d[schema.slice("CO2")] = co2_carbon / carbon_mass_fraction("CO2")
        # Refund the producer's sugar draw for the re-sourced fraction (the inverse of its draw),
        # so net sugar loss is only the un-rerouted remainder. rate > 0 ⇒ flux > 0 ⇒ sugar present,
        # so the refund always lands (no carbon leak).
        _refund_carbon_to_sugar(d, y, schema, refund)
        return d


class EsterVolatilization(Process):
    """CO2-stripping loss of liquid esters to the headspace — a physical Henry's-law sink.

    ``d(esters)/dt = -k · X·S_total/(K_sugar_uptake+S_total) · f_gas(T) · f_part(T) ·
    <ester>`` and the equal-and-opposite ``+`` into that ester's own ``<ester>_gas`` twin,
    for each ester in ``ESTER_SPECS`` (decision D-96), where

    * ``f_gas(T) = arrhenius_factor(T, E_a_uptake, T_ref)`` is the **gas-flow** factor: the
      stripping rides the evolving-CO2 stream, whose rate is the fermentative uptake flux
      scaled by the *same* Arrhenius factor the uptake Process carries (``E_a_uptake``).
      The bare ``X·S/(K+S)`` shape times ``f_gas`` is the CO2-evolution proxy; the
      ``q_sugar_max``/``co2_yield``/realised-yield constants fold into ``k_ester_volatil``.
    * ``f_part(T) = arrhenius_factor(T, dH_ester_volatil, T_ref)`` is the **gas/liquid
      partition** factor — a van't Hoff form for the ethyl-acetate Henry's-law constant,
      which *rises* with temperature (warmer ⇒ more volatile). ``dH_ester_volatil`` is the
      sourced partition enthalpy (~45 kJ/mol, ethyl acetate; NIST/Sander Henry data), a
      **physical Q10 ≈ 1.8**, *not* a fitted lever.

    So the total stripping temperature sensitivity is ``E_a_uptake + dH_ester_volatil`` ≈
    100 kJ/mol — the same physical value **in both media**, because the Henry's-law
    partition and the gas flow are properties of the molecule and the ferment, not the
    beverage (Morakul et al. 2011: the partition coefficient depends only on composition
    and temperature, not on which medium). The loss is **first-order in the liquid ester
    present** and **stops when fermentation stops** (``flux → 0`` at dryness) — a
    deliberate omission of slow passive evaporation after the cap goes on, keeping the
    sink a clean function of the gas stream (the ester analogue of Morakul's
    ``L = ∫ C_gas·Q_CO2 dt``).

    **Why a shared physical stripping is the faithful model (decisions D-20 → D-21).** D-20
    first built this sink with a *fudged per-medium* ``E_a_ester_volatil`` (wine above /
    beer below ``E_a_esters``) to make wine liquid esters fall while beer rises. But a
    *sourced* Henry's-law stripping is medium-independent, so it cannot push opposite
    directions by itself. The wine/beer divergence therefore lives where it is genuinely
    sourced — in ester **synthesis** ``E_a_esters``, set *per medium*: **beer** strongly
    T-sensitive (de Andrés-Toro 1998, ester ride the growth rate, apparent E_a ~200+
    kJ/mol) so synthesis outruns stripping and **beer liquid esters rise** with T; **wine**
    weak/~flat (Mouret 2015; Rollero 2014, wine ester synthesis is weak and non-monotonic
    in T) so stripping outruns synthesis and **wine liquid esters fall** with T (the
    Rollero evaporation inversion). Both directions now emerge from *physical + sourced*
    parameters, not a compensating constant. Verified empirically at 14/20/25 °C.

    **Carbon — three neutral liquid→gas transfers (no sugar draw).** Unlike
    :class:`EsterSynthesis`/:class:`FuselAlcoholsEhrlich`, this Process draws **no fresh
    sugar**: it only moves carbon already in each liquid ester pool into *that ester's own*
    headspace pool. So it touches the three ester pools and their three gas twins only —
    never ``S``/``E``/``CO2`` — and since ``total_carbon`` weights each pool and its twin at
    the *same* molecule's fraction, **every** transfer is carbon-neutral independently and
    closure stays at machine precision (a headspace pool is the ester analogue of evolved
    ``CO2``: carbon leaves the liquid but not the ledger). Each liquid pool is clamped ≥ 0 so
    a solver undershoot cannot strip a negative pool into spurious gas.

    Note the pairing is *why* each ester needs its own gas pool (decision D-96): a single
    shared headspace pool could carry only one carbon weight, so stripping a C7 ester into a
    C4-weighted pool would create or destroy carbon. Before D-96 one lumped pool and one gas
    twin were both booked as ethyl acetate, and the question could not arise.

    **Documented simplification — one partition enthalpy for three molecules (D-96).**
    ``dH_ester_volatil`` is sourced for **ethyl acetate** (~45 kJ/mol, NIST/Sander) and v1
    applies it to all three esters. That is a measured constant borrowed for two molecules it
    was not measured for — an approximation the ester van't Hoff enthalpies' natural clustering
    (~40–55 kJ/mol across this range) makes tolerable, but an approximation, and it is recorded
    here rather than left implicit. Its direction is known: isoamyl acetate and ethyl hexanoate
    are *less* volatile than ethyl acetate, so v1 slightly over-strips them. Per-ester sourced
    partition enthalpies are the named deferred refinement.

    **Documented simplification.** The full Morakul (2011) partition coefficient is also
    *ethanol-dependent* (``ln k_i = F1 + F2·E − (F3 + F4·E)·R·(1000/T − 1000/T_ref)``); we
    keep only the dominant temperature (van't Hoff) lever via ``dH_ester_volatil`` and omit
    the ethanol terms (the ``F`` coefficients are not openly available). Tier: **plausible**
    in form (gas stripping by the evolving CO2 is well-understood physics, the standard
    explanation for wine's liquid-ester temperature response), with speculative rate
    parameters that cap the pool outputs at speculative via parameter-tier propagation
    (D-1).
    """

    name = "ester_volatilization"
    tier = Tier.PLAUSIBLE
    #: Since D-115 the banana ester's label tracer joins them: stripping removes labelled and
    #: unlabelled molecules alike, so the tracer must fall with the pool or the enrichment would
    #: drift upward for no physical reason.
    touches = (
        *(spec.pool for spec in ESTER_SPECS),
        *(spec.gas_pool for spec in ESTER_SPECS),
        VALINE_LABEL_TRACERS[1].tracer_pool,
    )
    #: ``K_sugar_uptake``/``E_a_uptake`` are shared with the fermentative uptake whose CO2
    #: stream does the stripping (gas-flow factor); ``dH_ester_volatil`` is the sourced
    #: ethyl-acetate Henry's-law partition enthalpy (gas/liquid factor); ``T_ref`` anchors
    #: both. Their tiers cap the pool outputs via parameter-tier propagation (D-1). The
    #: wine/beer direction lives in ``E_a_esters`` (synthesis), not here — this stripping
    #: is the same physical mechanism in both media.
    reads: tuple[str, ...] = (
        "k_ester_volatil",
        "K_sugar_uptake",
        "E_a_uptake",
        "dH_ester_volatil",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_gas = arrhenius_factor(temp, params["E_a_uptake"], params["T_ref"])  # CO2 gas flow
        f_part = arrhenius_factor(temp, params["dH_ester_volatil"], params["T_ref"])  # partition
        for spec in ESTER_SPECS:
            liquid = max(float(y[schema.slice(spec.pool)][0]), 0.0)
            if liquid <= 0.0:  # nothing in this liquid pool to strip
                continue
            rate = params["k_ester_volatil"] * flux * f_gas * f_part * liquid
            # Liquid → its OWN headspace twin: the pair shares one molecule's carbon weight,
            # so each transfer is carbon-neutral independently (D-96).
            d[schema.slice(spec.pool)] = -rate
            d[schema.slice(spec.gas_pool)] = rate
            # D-115: stripping is NON-FRACTIONATING — the headspace takes labelled and unlabelled
            # molecules in the proportion the liquid holds them, so the tracer is debited at the
            # pool's own fraction and the *fraction* is left exactly where it was. Omitting this
            # would be the classic tracer bug: mass leaves, label does not, and the enrichment
            # climbs toward 100% purely because the pool shrank. (A real isotope effect on
            # volatility exists and is utterly negligible at these masses; not modelled, and it
            # would be a fractionation term here rather than an omission.)
            tracer = TRACER_BY_BULK.get(spec.pool)
            if tracer is not None:
                d[schema.slice(tracer.tracer_pool)] -= rate * labelled_fraction(y, schema, tracer)
        return d
