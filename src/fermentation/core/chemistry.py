"""Stoichiometric constants for the validated core.

Molar masses and carbon-atom counts of the species the core tracks. These are
*exact consequences of the chemical formulae* (and the standard atomic masses),
not empirical kinetic parameters — so, like the conversion factors in
:mod:`fermentation.units` (decision D-3), they live in code with citations rather
than in the provenance-backed parameter store.

This module is the **single source of truth** for fermentation stoichiometry: the
carbon/mass conservation checks (:mod:`fermentation.validation.conservation`) and
the sugar-uptake Process both derive their numbers here, so a conservation check
can never silently disagree with the kinetics it audits.

Empirical, strain-dependent quantities do **not** belong here. In particular the
*biomass elemental composition* (carbon/nitrogen content of dry yeast) is uncertain
and shared with the growth Process, so it is a :class:`~fermentation.parameters.\
schema.Parameter` (provenance store), passed into the conservation builders — see
decision D-8.

Standard atomic masses (IUPAC 2021, g/mol): C 12.011, H 1.008, O 15.999.
"""

from __future__ import annotations

from fermentation.core.state import StateSchema

# -- standard atomic masses (IUPAC 2021), g/mol -------------------------------
_M_C = 12.011
_M_H = 1.008
_M_O = 15.999
_M_S = 32.06
_M_N = 14.007

# -- molar masses of tracked species, g/mol (derived from their formulae) -----
#: Glucose / fructose, C6H12O6 — the lumped wine hexose and beer's first sugar.
M_GLUCOSE = 6 * _M_C + 12 * _M_H + 6 * _M_O
#: Maltose, C12H22O11.
M_MALTOSE = 12 * _M_C + 22 * _M_H + 11 * _M_O
#: Maltotriose, C18H32O16.
M_MALTOTRIOSE = 18 * _M_C + 32 * _M_H + 16 * _M_O
#: Ethanol, C2H6O.
M_ETHANOL = 2 * _M_C + 6 * _M_H + 1 * _M_O
#: Acetaldehyde (ethanal), C2H4O — the obligate intermediate on the main ethanol
#: pathway (pyruvate → acetaldehyde → ethanol), the carbonyl responsible for the early
#: "green apple" transient and the principal SO₂-binder (decision D-27). Same 2 carbons
#: as ethanol, so the yeast reduction acetaldehyde → ethanol is a mole-for-mole C2 → C2
#: transfer: modelling acetaldehyde as a transient buffer that *borrows* ethanol carbon
#: (production) and *returns* it (reduction) closes carbon to machine precision without
#: touching ``S`` or ``CO2`` — the faithful de-lumping of the uptake Process's single
#: sugar→ethanol step, chosen by the owner over a draw-from-sugar stand-in that would
#: double-count the main pathway and inflate ABV (decision D-27).
M_ACETALDEHYDE = 2 * _M_C + 4 * _M_H + 1 * _M_O
#: Carbon dioxide, CO2.
M_CO2 = 1 * _M_C + 2 * _M_O
#: Water, H2O (hydrolysis bookkeeping for di-/trisaccharide uptake).
M_WATER = 2 * _M_H + 1 * _M_O
#: Glycerol, C3H8O3 — the principal fermentation byproduct (realised-yield sink,
#: decision D-16).
M_GLYCEROL = 3 * _M_C + 8 * _M_H + 3 * _M_O
#: Succinic acid, C4H6O4 — the representative species for the lumped *minor*
#: byproduct pool. It carries the carbon of the ``Byp`` state variable so that
#: pool's carbon is accounted from a real formula rather than an ad-hoc fraction
#: (decision D-16). Under D-19 ``Byp`` is *organic acids / polyols only*: the higher
#: alcohols it formerly lumped now have their own carbon-routed ``fusels`` pool, so
#: there is no double-count between ``Byp`` (succinic) and ``fusels`` (isoamyl).
M_SUCCINIC = 4 * _M_C + 6 * _M_H + 4 * _M_O
#: Ethyl acetate, C4H8O2 — the representative species for the lumped ``esters`` aroma
#: pool, carbon-routed from sugar under decision D-19. BOOKKEEPING CAVEAT: a real
#: ester's ethanol moiety is carbon already counted in ``E``, so "route ester carbon
#: from sugar" over-attributes fresh hexose carbon — it closes the ledger exactly but
#: is an accounting stand-in, not a claim about the metabolic carbon origin (D-19).
M_ETHYL_ACETATE = 4 * _M_C + 8 * _M_H + 2 * _M_O
#: Isoamyl alcohol (3-methylbutan-1-ol), C5H12O — the representative species for the
#: lumped ``fusels`` higher-alcohol pool, carbon-routed from sugar under decision
#: D-19. BOOKKEEPING CAVEAT: the Ehrlich pathway builds fusels from amino-acid
#: skeletons, but ``N`` (YAN) carries no carbon in :func:`total_carbon`, so the
#: carbon is sourced from sugar as a stand-in — exact on the ledger, approximate on
#: the metabolism (D-19).
M_ISOAMYL_OH = 5 * _M_C + 12 * _M_H + 1 * _M_O
#: Tartaric acid, C4H6O6 — the dominant grape acid and the TA reference species
#: (equivalent weight M_TARTARIC/2 ≈ 75.04 g/eq). Diprotic; charge-active in the
#: wine pH solver (decision D-18).
M_TARTARIC = 4 * _M_C + 6 * _M_H + 6 * _M_O
#: L-malic acid, C4H6O5 — the second major grape acid and the MLF substrate.
#: Diprotic; a future MLF Process converts it to lactic + CO2 (4 = 3 + 1 carbons),
#: so these weights make that conversion carbon-closing (decision D-18).
M_MALIC = 4 * _M_C + 6 * _M_H + 5 * _M_O
#: L-lactic acid, C3H6O3 — the MLF product (produced-only). Monoprotic; the softer
#: acid that malic deacidifies *into*, the chemistry the pH solver must reproduce.
M_LACTIC = 3 * _M_C + 6 * _M_H + 3 * _M_O
#: Citric acid, C6H8O7 — the minor grape acid *Oenococcus oeni* co-metabolises during MLF,
#: the carbon source for **MLF-derived diacetyl** (decision D-31). A dosed must input
#: (~0.1–0.5 g/L), present *independent of sugar* so it can fund the diacetyl pool after the
#: wine is dry — the reason yeast-pathway sugar carbon (which no-ops at ``S=0``) cannot source
#: this beat. The v1 :class:`~fermentation.core.kinetics.malolactic.MalolacticCitrateMetabolism`
#: routes it as a lumped, carbon-closing ``citrate (C6) → α-acetolactate (C5) + CO2 (C1)``
#: stand-in feeding the shared VDK reservoir (6 = 5 + 1), so these weights make that
#: conversion carbon-closing on the existing ledger. BOOKKEEPING CAVEAT: real citrate
#: metabolism is ``citrate → acetate + oxaloacetate → pyruvate + CO2`` and takes ~2 citrate
#: per α-acetolactate, with acetate (a volatile-acidity contributor) as the *dominant*
#: co-product; the single-reaction stand-in balances carbon exactly but omits the acetate/
#: lactate branches and full citrate depletion (rate held low so citrate stays mostly
#: unconsumed — the trace diacetyl branch only, decision D-31). Triprotic; kept OUT of the
#: D-18 pH charge balance in v1 (a scoped omission the inverse anchoring absorbs at t=0, as
#: for SO₂'s bisulfite charge, D-22), so it is carbon-active but not charge-active.
M_CITRIC = 6 * _M_C + 8 * _M_H + 7 * _M_O
#: α-acetolactate (2-acetolactic acid), C5H8O4 — the vicinal-diketone (VDK) precursor
#: reservoir (decision D-26). Yeast excretes it during valine biosynthesis; it then
#: *spontaneously* (non-enzymatically) oxidatively decarboxylates to diacetyl + CO2,
#: the slow, temperature-critical step that makes the "diacetyl rest" a rest. The C5→C4
#: carbon (one carbon leaves as CO2) makes that decarboxylation carbon-closing on the
#: existing ledger, exactly as malic→lactic+CO2 (D-23). Better grounded than the
#: ester/fusel sugar stand-ins: α-acetolactate genuinely derives from pyruvate (sugar).
M_ACETOLACTATE = 5 * _M_C + 8 * _M_H + 4 * _M_O
#: Diacetyl (2,3-butanedione), C4H6O2 — the flavour-active vicinal diketone (buttery
#: off-note, the defining lager parameter). Produced by spontaneous decarboxylation of
#: α-acetolactate and reabsorbed by viable yeast, which reduces it to 2,3-butanediol —
#: the produce-then-reabsorb time course behind the diacetyl rest (decision D-26).
M_DIACETYL = 4 * _M_C + 6 * _M_H + 2 * _M_O
#: 2,3-Butanediol, C4H10O2 — the flavour-inactive terminal product of yeast diacetyl
#: reduction (via acetoin, lumped here into the diol; decision D-26). The real fate of
#: reabsorbed diacetyl, so tracking it as its own pool makes the reduction a genuine
#: carbon-conserving transfer (C4→C4, mole-for-mole) rather than a "returns-to-sugar"
#: bookkeeping stand-in — the fidelity the owner asked for over closure-only options.
M_BUTANEDIOL = 4 * _M_C + 10 * _M_H + 2 * _M_O
#: Sulfur dioxide, SO2 — the ``so2_free`` (free SO₂) state species (decision D-22).
#: The ONLY carbon-free tracked species, so it contributes nothing to ``total_carbon``
#: (registered with 0 carbon atoms below; cf. the charge-only ``cation_charge`` slot).
#: Free SO₂ is conventionally expressed *as SO₂* regardless of speciation, so the
#: pH-driven molecular/bisulfite/sulfite split is mass-preserving and the readout needs
#: no molar conversion; this molar mass is carried for completeness as the tracked
#: species' weight and for the deferred in-balance step (sulfurous-acid mol/L charge).
M_SO2 = 1 * _M_S + 2 * _M_O
#: Hydrogen sulfide, H2S — the "rotten egg" sulfidic off-aroma yeast releases when it
#: reduces sulfate faster than it can fix the sulfide onto nitrogen skeletons (the sulfate-
#: reduction sequence outruns the assimilation that needs O-acetylserine/-homoserine), so
#: production is de-repressed at low yeast-assimilable nitrogen (decision D-29). Like SO₂ it
#: is **carbon-free** (registered with 0 carbon atoms below), so it contributes nothing to
#: ``total_carbon`` and its produced-only pool sits on no conservation ledger — the sulfur it
#: carries is not tracked anywhere else (there is no sulfate/sulfur state), exactly as free
#: SO₂'s sulfur is not. This molar mass is carried for completeness (the tracked species'
#: weight) and for the deferred CO₂-stripping volatilization sink; the v1 production kinetics
#: work in g/L directly and need no molar conversion.
M_H2S = 2 * _M_H + 1 * _M_S
#: L-arginine, C6H14N4O2 — the representative species for the assimilable **amino-acid**
#: pool (decision D-32). Arginine is the *dominant* yeast-assimilable amino acid in grape
#: must (proline, though more abundant, is not assimilated anaerobically), so it is the
#: honest single-species stand-in for the lumped ``amino_acids`` pool — the succinic-for-
#: ``Byp`` / isoamyl-for-``fusels`` idiom (D-16/D-19). Crucially it is the FIRST tracked
#: species carrying **nitrogen** (four N per molecule; see ``NITROGEN_ATOMS`` below), so the
#: pool sits on *both* the carbon ledger (``total_carbon``) and the nitrogen ledger
#: (``total_nitrogen``) — the reason the amino-acid ledger needed a per-species nitrogen
#: accounting at all (nitrogen was previously tracked only as the elemental ``N`` slot plus
#: ``f_N·X``). Its mass C:N ratio (72.066 / 56.028 ≈ 1.29) is deliberately **N-rich** and
#: well below biomass's (``f_C/f_N`` ≈ 4.3): that is the load-bearing property that keeps the
#: :class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation` carbon refund
#: strictly below growth's sugar-carbon demand for any assimilation fraction ψ ≤ 1, so the
#: swap never creates hexose (gluconeogenesis) and needs no clamp — decision D-32.
M_ARGININE = 6 * _M_C + 14 * _M_H + 4 * _M_N + 2 * _M_O

#: Molar mass [g/mol] keyed by species name. ``fermentation.core.media`` sugar
#: component names ("glucose", "maltose", "maltotriose") are keys here.
MOLAR_MASS: dict[str, float] = {
    "glucose": M_GLUCOSE,
    "fructose": M_GLUCOSE,
    "maltose": M_MALTOSE,
    "maltotriose": M_MALTOTRIOSE,
    "ethanol": M_ETHANOL,
    "acetaldehyde": M_ACETALDEHYDE,
    "CO2": M_CO2,
    "glycerol": M_GLYCEROL,
    "succinic_acid": M_SUCCINIC,
    "ethyl_acetate": M_ETHYL_ACETATE,
    "isoamyl_alcohol": M_ISOAMYL_OH,
    "tartaric_acid": M_TARTARIC,
    "malic_acid": M_MALIC,
    "lactic_acid": M_LACTIC,
    "citric_acid": M_CITRIC,
    "sulfur_dioxide": M_SO2,
    "hydrogen_sulfide": M_H2S,
    "alpha_acetolactate": M_ACETOLACTATE,
    "diacetyl": M_DIACETYL,
    "butanediol": M_BUTANEDIOL,
    "arginine": M_ARGININE,
}

#: Carbon atoms per molecule, keyed by species name. The two sulfur species
#: ``sulfur_dioxide`` and ``hydrogen_sulfide`` are carried at **0** so
#: ``carbon_mass_fraction(...)`` returns 0.0 (not a KeyError) — the free-SO₂ pool (D-22) and
#: the H₂S pool (D-29) are correctly carbon-inert in any carbon sum.
CARBON_ATOMS: dict[str, int] = {
    "glucose": 6,
    "fructose": 6,
    "maltose": 12,
    "maltotriose": 18,
    "ethanol": 2,
    "acetaldehyde": 2,
    "CO2": 1,
    "glycerol": 3,
    "succinic_acid": 4,
    "ethyl_acetate": 4,
    "isoamyl_alcohol": 5,
    "tartaric_acid": 4,
    "malic_acid": 4,
    "lactic_acid": 3,
    "citric_acid": 6,
    "sulfur_dioxide": 0,
    "hydrogen_sulfide": 0,
    "alpha_acetolactate": 5,
    "diacetyl": 4,
    "butanediol": 4,
    "arginine": 6,
}

#: Nitrogen atoms per molecule, keyed by species name. Nitrogen was historically tracked
#: only as the elemental yeast-assimilable ``N`` slot (g N/L) plus the ``f_N·X`` bound in
#: biomass, so no species carried nitrogen — until the amino-acid pool (decision D-32).
#: Every carbon-tracked species is listed here at **0** except ``arginine`` (four N), exactly
#: mirroring ``CARBON_ATOMS`` so :func:`nitrogen_mass_fraction` returns 0.0 (not a KeyError)
#: for the carbon-only species and the check-vs-kinetics single-source discipline holds for
#: nitrogen as it does for carbon.
NITROGEN_ATOMS: dict[str, int] = {
    "glucose": 0,
    "fructose": 0,
    "maltose": 0,
    "maltotriose": 0,
    "ethanol": 0,
    "acetaldehyde": 0,
    "CO2": 0,
    "glycerol": 0,
    "succinic_acid": 0,
    "ethyl_acetate": 0,
    "isoamyl_alcohol": 0,
    "tartaric_acid": 0,
    "malic_acid": 0,
    "lactic_acid": 0,
    "citric_acid": 0,
    "sulfur_dioxide": 0,
    "hydrogen_sulfide": 0,
    "alpha_acetolactate": 0,
    "diacetyl": 0,
    "butanediol": 0,
    "arginine": 4,
}


def carbon_mass_fraction(species: str) -> float:
    """Grams of carbon per gram of ``species`` (exact from its formula).

    Used to weight each state variable when summing total carbon. Raises
    ``KeyError`` for an unknown species so a typo fails loudly rather than
    silently dropping a carbon-bearing term from a conservation check.
    """
    try:
        return CARBON_ATOMS[species] * _M_C / MOLAR_MASS[species]
    except KeyError:
        raise KeyError(f"unknown species {species!r}; known: {sorted(MOLAR_MASS)}") from None


def nitrogen_mass_fraction(species: str) -> float:
    """Grams of nitrogen per gram of ``species`` (exact from its formula).

    The nitrogen analogue of :func:`carbon_mass_fraction`, used to weight the
    amino-acid pool in :func:`~fermentation.validation.conservation.total_nitrogen`
    and to convert amino-acid mass to refunded ammonium ``N`` in the
    :class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation` swap
    (decision D-32). Returns 0.0 for the nitrogen-free species (all but arginine);
    raises ``KeyError`` for an unknown species so a typo fails loudly rather than
    silently dropping a nitrogen-bearing term from a conservation check.
    """
    try:
        return NITROGEN_ATOMS[species] * _M_N / MOLAR_MASS[species]
    except KeyError:
        raise KeyError(f"unknown species {species!r}; known: {sorted(MOLAR_MASS)}") from None


def sugar_species(schema: StateSchema) -> list[str]:
    """Map a schema's ``S`` slots to chemical species names, in slot order.

    Beer's ``S`` names its components (glucose/maltose/maltotriose); wine's single
    lumped slot is treated as a hexose (glucose). A multi-slot ``S`` without
    component names is an error — its carbon/mass weights are undefined without
    knowing which sugars occupy the slots.

    This is the single source of truth shared by every sugar-aware consumer —
    the carbon-conservation check (:mod:`fermentation.validation.conservation`)
    and the kinetic Processes that draw carbon from sugar — so a check can never
    disagree with the kinetics it audits (decision D-8). It lives here, in the
    core, because the validation layer may not be imported by the core that needs
    it (the one-directional dependency rule).
    """
    spec = schema.spec("S")
    if spec.components:
        return list(spec.components)
    if spec.size == 1:
        return ["glucose"]
    raise ValueError(
        f"sugar 'S' has {spec.size} slots but no component names; cannot assign carbon fractions"
    )


# -- Gay-Lussac stoichiometry of one hexose -----------------------------------
# Anaerobic alcoholic fermentation of a hexose:
#
#     C6H12O6 -> 2 C2H5OH + 2 CO2
#
# mass-balanced by atom count (2*M_ETHANOL + 2*M_CO2 == M_GLUCOSE to the third
# decimal). These are the *theoretical* maximum yields; the realised ethanol
# yield is a few percent lower because cells divert carbon to glycerol, organic
# acids and biomass. The sugar-uptake Process applies that realised-yield
# correction by scaling this theoretical split down and routing the diverted
# carbon into the ``Gly``/``Byp`` byproduct pools (decision D-16); the split is
# exposed here so the Process and the carbon-conservation check use one
# definition.

#: Grams of ethanol produced per gram of hexose consumed (theoretical, ~0.511).
ETHANOL_PER_HEXOSE = 2 * M_ETHANOL / M_GLUCOSE
#: Grams of CO2 evolved per gram of hexose consumed (theoretical, ~0.489).
CO2_PER_HEXOSE = 2 * M_CO2 / M_GLUCOSE

#: Hexose units released per molecule on complete hydrolysis. Glucose/fructose are
#: hexoses already; maltose -> 2, maltotriose -> 3. Each released hexose ferments
#: by Gay-Lussac (-> 2 ethanol + 2 CO2), so a sugar's per-gram ethanol/CO2 yield
#: scales with its hexose count. The di-/trisaccharide mass gain comes from
#: hydrolysis water pulled from the solvent (why beer's S+E+CO2 mass does not
#: close — see ``validation.total_mass`` and decision D-8).
HEXOSE_UNITS: dict[str, int] = {
    "glucose": 1,
    "fructose": 1,
    "maltose": 2,
    "maltotriose": 3,
}


def ethanol_yield(species: str) -> float:
    """Grams of ethanol per gram of ``species`` fermented (theoretical Gay-Lussac).

    Generalises :data:`ETHANOL_PER_HEXOSE` to the di-/trisaccharides: a sugar with
    ``n`` hexose units yields ``2n`` ethanol per molecule. M1 uses this theoretical
    split (not the realised ``Y_ethanol_sugar`` parameter) so carbon and mass close
    exactly — the realised yield, net of glycerol/biomass diversion, is a Tier-2
    concern (decision D-8). Raises ``KeyError`` for an unknown/unfermentable
    species so a typo fails loudly.
    """
    try:
        return 2 * HEXOSE_UNITS[species] * M_ETHANOL / MOLAR_MASS[species]
    except KeyError:
        raise KeyError(
            f"no fermentation yield for {species!r}; known: {sorted(HEXOSE_UNITS)}"
        ) from None


def co2_yield(species: str) -> float:
    """Grams of CO2 per gram of ``species`` fermented (theoretical Gay-Lussac).

    Companion to :func:`ethanol_yield`; for a hexose equals :data:`CO2_PER_HEXOSE`.
    ``ethanol_yield(s) + co2_yield(s)`` exceeds 1 for di-/trisaccharides by exactly
    the hydrolysis water taken up — mass closes only when that water is tracked,
    which M1 does not (decision D-8), so beer relies on the carbon balance.
    """
    try:
        return 2 * HEXOSE_UNITS[species] * M_CO2 / MOLAR_MASS[species]
    except KeyError:
        raise KeyError(
            f"no fermentation yield for {species!r}; known: {sorted(HEXOSE_UNITS)}"
        ) from None
