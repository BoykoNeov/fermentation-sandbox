"""Shared carbon-routing helpers for the metabolite-byproduct Processes.

Two small building blocks reused by every Process that produces a carbon-accounted
byproduct pool from the fermentative flux:

* :func:`fermentative_flux_shape` — the biomass-catalysed sugar Monod activity the
  byproduct Processes couple to, so their production tracks the *same* flux they are
  metabolically downstream of.
* :func:`draw_carbon_from_sugar` — routes a byproduct's carbon *out of* ``S`` (option a1,
  decision D-19) so the pool it fills is real carbon-accounted state and ``total_carbon``
  closes to machine precision.

Extracted from :mod:`fermentation.core.kinetics.byproducts` (decision D-26) so the ester/
fusel aroma Processes and the vicinal-diketone (diacetyl) Processes share one definition —
the same single-source-of-truth discipline the chemistry constants follow. The behaviour
is unchanged from the byproducts beat (D-19); only the home moved.

It also hosts :data:`ESTER_SPECS`, the **canonical ester registry** (decision D-96) — the
single source of truth every layer derives from (state slots, synthesis, stripping, the
carbon ledger, the OAV aroma set), so a fourth ester is *one entry*, not a new code path.
"""

from __future__ import annotations

from dataclasses import dataclass

from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.state import FloatArray, StateSchema


@dataclass(frozen=True)
class EsterSpec:
    """One aroma-active ester: its liquid pool, headspace twin, molecule and rate constant.

    ``pool`` is the liquid state variable, ``gas_pool`` its CO₂-stripped headspace
    bookkeeping twin (decision D-20), ``species`` the molecule **both** are carbon-weighted
    by in :func:`~fermentation.validation.conservation.total_carbon` — weighting a
    liquid/gas pair at the *same* fraction is what makes stripping carbon-neutral — and
    ``k_param`` the per-medium synthesis rate constant.

    **Each pool IS its molecule (decision D-96).** Before D-96 a single lumped ``esters``
    pool was carbon-weighted as *ethyl acetate* but read on the OAV lens against *isoamyl
    acetate*'s threshold — a split identity that made the fruity OAV non-physical (~761 for a
    wine, implying ~23 mg/L isoamyl acetate against a real ceiling of ~1–3 mg/L). Splitting
    the lump into single-molecule pools, each weighted **and** perceived as itself, removes
    that seam rather than papering over it: no ester pool carries a ``lumped`` composition
    assumption any more.

    ``precursor_pool`` (decision D-97) names the liquid state variable holding the ester's
    **precursor alcohol**, when — and only when — that alcohol is scarce enough for its
    supply to limit the rate. It is ``None`` for an ester whose precursor is a bulk pool that
    saturates the enzyme (see :data:`ESTER_SPECS` for the concentration argument that decides
    this per ester; it is *derived*, never asserted).
    """

    pool: str
    gas_pool: str
    species: str
    k_param: str
    #: Human-readable note for the liquid pool's :class:`~fermentation.core.state.VarSpec`
    #: description — kept here so the schema text cannot drift from the registry.
    note: str
    #: Liquid pool of the precursor alcohol this ester's synthesis is **first-order in**
    #: (decision D-97), or ``None`` when the precursor saturates the enzyme and the rate is
    #: zeroth-order in it. Read only — never debited (the carbon still comes from ``S``).
    precursor_pool: str | None = None


#: The THREE aroma-active esters the sim tracks (decision D-96) — each a **single-molecule**
#: pool, so none carries the fixed-lump-composition caveat (that now survives only on
#: ``fusels``/``mercaptans``, where it is true). Two families:
#:
#: * the **acetate esters** (ATF1/alcohol acetyltransferase) — ``ethyl_acetate``, the bulk
#:   solventy one (tens of mg/L), and ``isoamyl_acetate``, the trace potent banana one
#:   (~1–3 mg/L). Same enzyme ⇒ they share ``E_a_esters``, the per-medium sourced
#:   temperature shape (decision D-21).
#: * the **ethyl ester of a medium-chain fatty acid** — ``ethyl_hexanoate`` (apple/pineapple),
#:   the highest-OAV ester in wine. **Documented v1 simplification (D-96):** it is
#:   EEB1/EHT1-derived, a *different* enzyme from the acetates, but v1 shares
#:   ``E_a_esters`` because no separate activation energy is sourced. A per-ester ``E_a`` is
#:   the named deferred refinement.
#:
#: Each ``k`` is sourced **independently** to its own molecule's measured concentration range
#: (D-96, the load-bearing call): the lump's composition is *derived* from three
#: independently-anchored rates, never a single fitted ratio splitting one ``k_ester``. Adding
#: a fourth ester (ethyl octanoate, phenylethyl acetate) is one entry here plus its params.
#:
#: **The ATF1 precursor coupling — ONE enzyme, ONE rate law, TWO limits (decision D-97).**
#: The two acetate esters are made by the same enzyme from the same acetyl-CoA donor,
#: differing only in the alcohol it acetylates. Whether that alcohol *limits* the rate is
#: therefore not a modelling preference — it is decided by comparing each alcohol's
#: concentration to ATF1's measured ``Km`` for it (~29.8 mM for isoamyl alcohol; Fujii 1998,
#: Appl. Environ. Microbiol. 64:4076-4078, citing Yoshioka & Hashimoto 1981):
#:
#: * ``isoamyl_acetate`` — its precursor is *isoamyl alcohol*, the ``fusels`` pool, which runs
#:   ~0.5-1 mM: **~30-60x BELOW Km** ⇒ the enzyme sits far down its linear stretch ⇒ the rate
#:   is **first-order in the pool** (``precursor_pool="fusels"``). Fujii states the conclusion
#:   outright: *"a major rate-limiting factor for isoamyl acetate production is the amount of
#:   isoamyl alcohol in the sake mash"*. This is what makes the banana note **YAN-responsive**
#:   — it inherits the nitrogen dependence of the Ehrlich pathway that builds its precursor.
#: * ``ethyl_acetate`` — its precursor is *ethanol*, the bulk ``E`` pool, which runs ~2 M:
#:   **orders of magnitude ABOVE any mM-scale Km** ⇒ the enzyme is saturated in it ⇒ the rate
#:   is **zeroth-order** in the precursor ⇒ no term at all (``precursor_pool=None``). Ethanol
#:   is never scarce, so it can never limit.
#:
#: So the asymmetry between two esters of the *same* enzyme is **derived from their precursor
#: concentrations**, not an exemption granted to one of them. ``ethyl_hexanoate`` is ungated
#: for a different reason: a different enzyme (EEB1/EHT1) acylating from hexanoyl-CoA, a
#: precursor the sim does not model at all — so there is no pool to be first-order in.
ESTER_SPECS: tuple[EsterSpec, ...] = (
    EsterSpec(
        "ethyl_acetate",
        "ethyl_acetate_gas",
        "ethyl_acetate",
        "k_ethyl_acetate",
        note="ethyl acetate (C4H8O2) — the bulk solventy/nail-polish acetate ester (ATF1)",
        # No precursor term: ethanol (~2 M) saturates ATF1 ⇒ zeroth-order in it (D-97).
        precursor_pool=None,
    ),
    EsterSpec(
        "isoamyl_acetate",
        "isoamyl_acetate_gas",
        "isoamyl_acetate",
        "k_isoamyl_acetate",
        note="isoamyl acetate (C7H14O2) — the trace, potent BANANA acetate ester (ATF1); "
        "the only ester the D-69 aging hydrolysis fades",
        # First-order in isoamyl alcohol: the pool runs far below ATF1's Km ~29.8 mM
        # (Fujii 1998) ⇒ the [S] << Km limit. This is the D-97 coupling. Since D-99 it names
        # the isoamyl pool SPECIFICALLY rather than the old lump — which is what the Km was
        # always measured for. Fujii's Km is for 3-methylbutan-1-ol; pointing this at the
        # lump meant a mixture stood in for the molecule the enzyme actually sees, and
        # pointing it at `active_amyl_alcohol` (the C5 ISOMER) would be a different enzyme
        # substrate entirely. Hence the registry constant, not a literal.
        precursor_pool="isoamyl_alcohol",
    ),
    EsterSpec(
        "ethyl_hexanoate",
        "ethyl_hexanoate_gas",
        "ethyl_hexanoate",
        "k_ethyl_hexanoate",
        note="ethyl hexanoate (C8H16O2) — the APPLE/PINEAPPLE ethyl ester of a "
        "medium-chain fatty acid (EEB1/EHT1)",
    ),
)

#: The ester the D-69 aging hydrolysis acts on — the **banana** acetate. Its fade to
#: ``fusels`` + ``Byp`` is the whole point of that Process, and at D-96 the 5:2 carbon split
#: became **exact**: the debited molecule finally *is* isoamyl acetate (C7 → isoamyl alcohol
#: C5 + acetic acid C2), where D-69 had to debit ethyl acetate and split as if it were
#: isoamyl acetate. See :class:`~fermentation.core.kinetics.aging.EsterHydrolysis`.
HYDROLYSING_ESTER: EsterSpec = ESTER_SPECS[1]


@dataclass(frozen=True)
class FuselSpec:
    """One Ehrlich higher alcohol: its pool, molecule, rate constant and amino-acid precursor.

    The fusel twin of :class:`EsterSpec`, and deliberately simpler in one way: there is **no
    ``gas_pool``**. The esters carry a CO₂-stripped headspace twin (D-20) because they are
    volatile enough to be swept out during a vigorous ferment; the higher alcohols are not
    stripped in this model, so each is a liquid pool only.

    ``species`` is the molecule the pool is carbon-weighted by in
    :func:`~fermentation.validation.conservation.total_carbon`, and — since D-99 — **it is the
    pool's own molecule, not a stand-in**. ``precursor_amino_acid`` is documentation of the
    Ehrlich route, not a modelled pool: `N` (YAN) is a single lumped nitrogen number here, so
    the sim cannot know which amino acid was consumed. It is recorded because it is *why*
    these five and not others, and because it names exactly what a future speciated-YAN model
    would have to supply to make the composition dynamic (see :data:`FUSEL_SPECS`).
    """

    pool: str
    species: str
    k_param: str
    #: Human-readable note for the pool's :class:`~fermentation.core.state.VarSpec`
    #: description — kept here so the schema text cannot drift from the registry.
    note: str
    #: The Ehrlich-pathway amino acid whose skeleton becomes this alcohol. Documentation of
    #: provenance only — never read by any rate law (see the class doc).
    precursor_amino_acid: str


#: The FIVE Ehrlich higher alcohols the sim tracks (decision D-99) — each a **single-molecule**
#: pool, so none carries the fixed-lump-composition caveat. This is the D-96 ester split
#: applied one pool over, and for the same reason: until D-99 one lumped ``fusels`` pool was
#: carbon-weighted as isoamyl alcohol AND read on the OAV lens against isoamyl alcohol's
#: threshold — so unlike the pre-D-96 ester pool it had no *split identity*, but it made a
#: quieter and broader error: it asserted that **every** higher alcohol smells like isoamyl
#: alcohol and weighs like it. Four of the five are neither.
#:
#: **Each ``k`` is anchored INDEPENDENTLY to its own molecule's measured concentration** — the
#: load-bearing D-96 rule, and the reason this beat is worth anything. The lump's composition
#: is now *derived* from five independently-anchored rates; a ratio-split off one ``k_fusel``
#: would have smuggled back the fabricated composition constant the split exists to remove.
#: Concentrations come from Wang, Frank & Steinhaus 2024 (J. Agric. Food Chem. 72:22250-22257),
#: a meta-analysis reporting a per-compound MEAN over N independently published studies —
#: n=486 (isobutanol) / 128 (active amyl) / 555 (isoamyl) / 684 (2-phenylethanol) for wine.
#: The honest consequence is that the total RISES ~3.8× (wine ~86 → ~328 mg/L): the old lump
#: sat below even the *sum of the five species' low ends*, and ``k_fusel``'s own provenance
#: already admitted it ran at "the low end of the 150-400 mg/L wine higher-alcohol range".
#: The rise is **forced by honest per-molecule anchoring**, not chosen for an outcome.
#:
#: **The two C5 isomers are the trap.** ``isoamyl_alcohol`` (3-methylbutan-1-ol, CAS 123-51-3)
#: and ``active_amyl_alcohol`` (2-methylbutan-1-ol, CAS 137-32-6) share a formula, a molar mass
#: and a carbon count, and differ ~5.5× in odour potency. They coelute on common GC phases, so
#: routine OIV/BIPEA methods report a single combined "amyl alcohols" figure — but aroma
#: research resolves and quantifies them separately (n=128 wine / n=64 beer studies for active
#: amyl alone), which is what makes five independently-anchored pools sourceable rather than a
#: ratio-split. See :mod:`fermentation.core.chemistry` for the CAS warning; vendor literature
#: has been seen to invert these two names.
#:
#: **What this split does NOT buy — the honest limit (D-99).** All five share one N-gate and
#: one ``E_a_fusels``, so the *spectrum* is fixed even though the *molecules* are now real:
#: raise the temperature or starve the nitrogen and all five move together. Speciation retires
#: the wrong-molecule-threshold error; it replaces "fixed composition" with the weaker but
#: still-real "fixed spectrum". Making the spectrum dynamic needs per-species activation
#: energies and per-amino-acid gates — **deferred, because neither is sourced**, and an author
#: estimate per species would lower fidelity while looking like it raised it.
FUSEL_SPECS: tuple[FuselSpec, ...] = (
    FuselSpec(
        "propanol",
        "propanol",
        "k_propanol",
        note="propan-1-ol (C3H8O) — the shortest Ehrlich higher alcohol; no sourced odour "
        "threshold in either matrix, so it is chemistry-only (carries no OAV)",
        precursor_amino_acid="threonine",
    ),
    FuselSpec(
        "isobutanol",
        "isobutanol",
        "k_isobutanol",
        note="isobutanol / 2-methylpropan-1-ol (C4H10O) — fusel/alcoholic; aroma-active in "
        "wine (Guth threshold ~40 mg/L), chemistry-only in beer (no beer threshold sourced)",
        precursor_amino_acid="valine",
    ),
    FuselSpec(
        "active_amyl_alcohol",
        "active_amyl_alcohol",
        "k_active_amyl_alcohol",
        note="active amyl alcohol / 2-methylbutan-1-ol (C5H12O, CAS 137-32-6) — an ISOMER of "
        "isoamyl alcohol, not a synonym; no sourced odour threshold in either matrix, so it "
        "is chemistry-only (carries no OAV)",
        precursor_amino_acid="isoleucine",
    ),
    FuselSpec(
        "isoamyl_alcohol",
        "isoamyl_alcohol",
        "k_isoamyl_alcohol",
        note="isoamyl alcohol / 3-methylbutan-1-ol (C5H12O, CAS 123-51-3) — the DOMINANT "
        "higher alcohol of both media and the solventy/fusel note; the D-97 precursor of "
        "isoamyl acetate and the C5 half of the D-69 hydrolysis split",
        precursor_amino_acid="leucine",
    ),
    FuselSpec(
        "2_phenylethanol",
        "2_phenylethanol",
        "k_2_phenylethanol",
        note="2-phenylethanol (C8H10O) — the only AROMATIC higher alcohol: rose/honey, not "
        "solventy. Aroma-active in wine (Guth threshold ~10 mg/L vs ~28.7 mg/L typical), "
        "chemistry-only in beer (no beer threshold sourced)",
        precursor_amino_acid="phenylalanine",
    ),
)

#: The higher alcohol the D-97 ATF1 coupling and the D-69 aging hydrolysis both act on — the
#: **isoamyl** one, in both cases because the ester in question is isoamyl acetate. Named here
#: so neither Process has to spell the string, and so the D-99 split cannot silently re-point
#: either at the wrong C5 isomer (``active_amyl_alcohol`` is a different molecule).
ISOAMYL_ALCOHOL: FuselSpec = FUSEL_SPECS[3]


def draw_carbon_from_sugar(
    d: FloatArray, y: FloatArray, schema: StateSchema, carbon: float
) -> None:
    """Subtract ``carbon`` [g C/L/h] from ``S``, split across slots by carbon content.

    Routes a byproduct's carbon out of sugar (option a1, decision D-19) so the pool it
    fills becomes carbon-accounted state. Each slot ``i`` gives up carbon in proportion
    to the carbon it currently holds, ``s_i·c_i / Σ_j s_j·c_j``; converting that back to
    grams of sugar, ``d[S_i] -= carbon · s_i / Σ_j s_j·c_j``. By construction
    ``Σ_i (d[S_i]·c_i) = -carbon`` exactly, so the carbon removed from sugar equals the
    carbon the caller deposits in its pool and ``total_carbon`` closes to machine
    precision. Slots are clamped ≥ 0 (mirroring the flux/uptake guards) so a solver
    undershoot cannot flip a draw into sugar *creation*; with no sugar carbon present
    (``Σ s_j c_j ≤ 0``) nothing is drawn. This serves wine's single slot and beer's
    three identically (different carbon fractions per slot are handled exactly).
    """
    s_slice = schema.slice("S")
    species = sugar_species(schema)
    s = [max(float(y[s_slice.start + i]), 0.0) for i in range(len(species))]
    carbon_total = sum(s[i] * carbon_mass_fraction(sp) for i, sp in enumerate(species))
    if carbon_total <= 0.0:
        return
    for i in range(len(species)):
        if s[i] > 0.0:
            d[s_slice.start + i] -= carbon * s[i] / carbon_total


def refund_carbon_to_sugar(
    d: FloatArray, y: FloatArray, schema: StateSchema, carbon: float
) -> None:
    """Add ``carbon`` [g C/L/h] back to ``S``, split across slots by carbon content.

    The exact inverse of :func:`draw_carbon_from_sugar`: where a byproduct Process *drew*
    its carbon out of sugar, another Process may *refund* carbon to sugar (e.g. when a
    fraction of that carbon is instead sourced from the amino-acid pool — the D-33 fusel
    Ehrlich re-route — or when amino-acid-funded biomass spares sugar for ethanol — the
    D-32 swap). Each slot ``i`` receives carbon in proportion to the carbon it currently
    holds, so ``Σ_i (d[S_i]·c_i) = +carbon`` exactly and the refund matches whatever draw
    it undoes to machine precision. Uses ``+=`` so it composes with any draw the same
    Process makes. Slots are clamped ≥ 0 (mirroring the draw); with no sugar carbon present
    (``Σ s_j c_j ≤ 0``) nothing is refunded — callers must guarantee sugar is present (the
    fermentative flux they track is zero at ``S = 0``, so this edge is never hit in
    practice). Serves wine's single slot and beer's three identically.
    """
    s_slice = schema.slice("S")
    species = sugar_species(schema)
    s = [max(float(y[s_slice.start + i]), 0.0) for i in range(len(species))]
    carbon_total = sum(s[i] * carbon_mass_fraction(sp) for i, sp in enumerate(species))
    if carbon_total <= 0.0:
        return
    for i in range(len(species)):
        if s[i] > 0.0:
            d[s_slice.start + i] += carbon * s[i] / carbon_total


def fermentative_flux_shape(y: FloatArray, schema: StateSchema, k_sat: float) -> float:
    """Biomass-catalysed sugar Monod term ``X · S_total/(K + S_total)`` [g/L].

    The dimensionless-but-for-``X`` activity proxy the fermentative uptake Process
    runs on (``q_sugar_max·X·S/(K+S)``), reused by the byproduct Processes so their
    production tracks the *same* flux they are metabolically coupled to — which is what
    makes the run-integrated "total scales as f_byproduct/f_flux" cancellation clean and
    predictable (see the byproducts module docstring). Sugar is summed across slots
    (1 for wine, 3 for beer) and clamped ≥ 0 against solver undershoot, mirroring the
    guards in the uptake/growth Processes.
    """
    x = max(float(y[schema.slice("X")][0]), 0.0)
    s_total = max(float(y[schema.slice("S")].sum()), 0.0)
    if x <= 0.0 or s_total <= 0.0:
        return 0.0
    return x * (s_total / (k_sat + s_total))
