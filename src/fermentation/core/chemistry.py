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
#: Pyruvic acid, C3H4O3 — the terminal glycolytic keto-acid and the immediate precursor of
#: acetaldehyde (pyruvate → acetaldehyde + CO2 → ethanol). Yeast **excretes** overflow pyruvate
#: during active fermentation (an extracellular residual, 10s–100s mg/L, peaking mid-ferment and
#: slowly re-assimilated), and that persistent excreted pool is — after acetaldehyde — the second
#: strongest SO₂-binding carbonyl in wine (decision D-49). Modelled as an **excreted side pool**
#: (drawn from sugar like the D-19/D-26 byproducts, viable-yeast-gated re-assimilation returning
#: its carbon to ethanol + CO2 — a carbon-closing C3 → C2 + C1 step, like malic → lactic + CO2),
#: NOT as an on-pathway intermediate: the intracellular flux pyruvate never persists and never
#: binds SO₂, so conflating the two would be unphysical (the rejected D-27-rework, D-49).
M_PYRUVATE = 3 * _M_C + 4 * _M_H + 3 * _M_O
#: α-Ketoglutaric acid (2-oxoglutaric acid), C5H6O5 — the second excreted overflow keto-acid
#: SO₂-binder (decision D-50), after pyruvate (D-49). Real yeast physiology de-represses α-KG
#: dehydrogenase under fermentative/anaerobic conditions, so it (like pyruvate) overflows to an
#: extracellular excreted residual rather than cycling through the TCA cycle; its "real" forward
#: fate is largely glutamate synthesis (α-KG + NH4+ → glutamate, an N-coupled route this v1 does
#: not model), so — exactly like pyruvate's C3 → C2(ethanol) + C1(CO2) — the reassimilation carbon
#: destination is a carbon-closing lumped stand-in, not a metabolic claim (see the ``keto_acids``
#: module docstring). Modelled with the SAME excreted-side-pool structure as pyruvate.
M_ALPHA_KETOGLUTARATE = 5 * _M_C + 6 * _M_H + 5 * _M_O
#: α-Ketobutyric acid (2-oxobutanoic acid / 2-ketobutyrate — all three names are one molecule),
#: C4H6O3 — the **third** excreted overflow keto-acid (decision D-107), and the only one in the
#: model with a *consumer*: it is the C4 half of sotolon's aldol.
#:
#: **Why an excreted pool and not an intracellular flux intermediate** (the D-49 question, which
#: this molecule answers the other way): sotolon forms in a **sealed bottle with no living yeast**,
#: by a purely chemical aldol (Pham *et al.* 1995) — so the α-ketobutyrate it condenses **must** be
#: an extracellular residual present at bottling. Pons *et al.* 2010 measures exactly that: wine
#: yeast strains "released large quantities of 2-ketobutyric acid" during alcoholic fermentation,
#: "similar to those found in oxidized dry white wines". D-49 rejected routing acetaldehyde through
#: pyruvate because the persistent pool and the flux intermediate are physically distinct; that same
#: test here **selects** the excreted pool, because the reaction that consumes it happens where no
#: intracellular pool can reach.
#:
#: Its **sources** are threonine deamination (``L-threonine → 2-oxobutanoate + NH₃``, ILV1) and
#: methionine demethiolation (``methionine → methanethiol + 2-oxobutyrate + NH₃``, D-45) — the
#: second is why this pool unblocks the mercaptan's 5× under-draw: the four carbons D-45 had to
#: discard finally have somewhere to go.
M_ALPHA_KETOBUTYRATE = 4 * _M_C + 6 * _M_H + 3 * _M_O
#: 5-Oxo-D-fructose (5-keto-fructose), C6H10O6 — the dominant BOTRYTIS-specific carbonyl
#: SO₂-binder (decision D-130). A hexodiulose formed by *Botrytis cinerea* (and acetic-acid
#: bacteria) oxidising must fructose ON THE BERRY, pre-crush, so it is a **must-composition
#: input**, not a fermentation product — the model carries it as an inert dosed slot, never
#: produced or consumed by any Process (Barbe 2000, J. Agric. Food Chem. 48:3413; Ribéreau-Gayon
#: et al., Handbook of Enology Vol. 1 §8.4). It is **not metabolized by yeast** ("concentrations
#: not altered by alcoholic fermentation", Handbook via Barbe), so unlike the transient
#: dihydroxyacetone it persists unchanged into finished wine and accounts for the bulk of a
#: botrytized wine's high SO₂-combining power. Carbon-bearing (C6, weighted in ``total_carbon``);
#: as an inert pool it is a constant term that drifts 0, so conservation is untouched.
M_5_OXOFRUCTOSE = 6 * _M_C + 10 * _M_H + 6 * _M_O
#: Carbon dioxide, CO2.
M_CO2 = 1 * _M_C + 2 * _M_O
#: Water, H2O (hydrolysis bookkeeping for di-/trisaccharide uptake).
M_WATER = 2 * _M_H + 1 * _M_O
#: Dioxygen, O2 — the dissolved-oxygen aging substrate (decision D-71, the first §4.1
#: OXIDATIVE aging Process). Carried as a plain molar mass (like ``M_WATER``), NOT registered
#: in ``MOLAR_MASS``/``CARBON_ATOMS`` below: the ``o2`` state pool is **carbon-free and off
#: every conservation ledger** (``total_carbon``/``total_mass``/``total_nitrogen`` weight only
#: their explicitly-named pools), exactly like the untracked hydrolysis water — so O₂ needs no
#: species registration, only this weight to convert the g/L O₂ consumed into the moles that set
#: the acetaldehyde yield (``y_acetaldehyde_per_o2`` mol acetaldehyde per mol O₂).
M_O2 = 2 * _M_O
#: Glycerol, C3H8O3 — the principal fermentation byproduct (realised-yield sink,
#: decision D-16).
M_GLYCEROL = 3 * _M_C + 8 * _M_H + 3 * _M_O
#: Succinic acid, C4H6O4 — the representative species for the lumped *minor*
#: byproduct pool. It carries the carbon of the ``Byp`` state variable so that
#: pool's carbon is accounted from a real formula rather than an ad-hoc fraction
#: (decision D-16). Under D-19 ``Byp`` is *organic acids / polyols only*: the higher
#: alcohols it formerly lumped now have their own carbon-routed pools (the five
#: ``FUSEL_SPECS``, D-99), so there is no double-count between ``Byp`` (succinic) and those
#: higher-alcohol pools.
M_SUCCINIC = 4 * _M_C + 6 * _M_H + 4 * _M_O
#: Ethyl acetate, C4H8O2 — the ``ethyl_acetate`` pool's own species (decision D-96: the
#: pool *is* this molecule; before D-96 it stood in for a lumped ``esters`` pool whose OAV
#: was read as a *different* molecule — see the module note below). Carbon-routed from
#: sugar under decision D-19. BOOKKEEPING CAVEAT: a real ester's ethanol moiety is carbon
#: already counted in ``E``, so "route ester carbon from sugar" over-attributes fresh
#: hexose carbon — it closes the ledger exactly but is an accounting stand-in, not a claim
#: about the metabolic carbon origin (D-19).
M_ETHYL_ACETATE = 4 * _M_C + 8 * _M_H + 2 * _M_O
#: Isoamyl acetate (3-methylbutyl acetate), C7H14O2 — the ``isoamyl_acetate`` pool's own
#: species (decision D-96). The potent *banana* acetate ester: trace by mass but low-threshold,
#: so it dominates the fruity note that ethyl acetate's bulk mass cannot explain. Same D-19
#: sugar-carbon stand-in caveat as ``ethyl_acetate``.
M_ISOAMYL_ACETATE = 7 * _M_C + 14 * _M_H + 2 * _M_O
#: Ethyl hexanoate (ethyl caproate), C8H16O2 — the ``ethyl_hexanoate`` pool's own species
#: (decision D-96). The representative *ethyl ester of a medium-chain fatty acid* — the
#: apple/pineapple half of the fruity bouquet, and the highest-OAV ester in wine. Unlike the
#: two acetates (ATF1) it is EEB1/EHT1-derived; v1 shares the acetates' ``E_a_esters``
#: temperature shape, a documented simplification (D-96). Same D-19 carbon stand-in caveat.
M_ETHYL_HEXANOATE = 8 * _M_C + 16 * _M_H + 2 * _M_O
#: THE FIVE EHRLICH HIGHER ALCOHOLS (decision D-99) — each its own pool, each carbon-weighted
#: by its OWN molecule. Until D-99 a single lumped ``fusels`` pool stood for all five and was
#: weighted (and perceived) as isoamyl alcohol alone; the split is the D-96 ester pattern one
#: pool over. Each ``k`` is anchored **independently** to its own molecule's measured
#: concentration (Wang, Frank & Steinhaus 2024) — never a ratio-split off a lumped total.
#:
#: ⚠ **NOMENCLATURE — anchor on CAS, not on the trivial names.** ``isoamyl_alcohol`` and
#: ``active_amyl_alcohol`` are STRUCTURAL ISOMERS: same formula C5H12O, same molar mass,
#: ~5.5× apart in odour potency. Vendor literature is not reliable here — Restek's 2025
#: "Alcoholic Beverage Analysis by GC" note has the two synonyms **inverted**. The bindings
#: below are correct and are the ones to trust:
#:   * isoamyl alcohol      = 3-methylbutan-1-ol = **CAS 123-51-3**
#:   * active amyl alcohol  = 2-methylbutan-1-ol = **CAS 137-32-6**
#: (The physical cause of the confusion is real: the two coelute on common GC phases such as
#: Rtx-VMS "regardless of column dimensions", so routine OIV/BIPEA methods report one combined
#: "amyl alcohols" figure. That is a column artifact, not a fact about the molecules — aroma
#: research resolves them, and quantifies them separately.)
#:
#: BOOKKEEPING CAVEAT, shared by all five: the Ehrlich pathway builds higher alcohols from
#: amino-acid skeletons, but ``N`` (YAN) carries no carbon in :func:`total_carbon`, so the
#: carbon is sourced from sugar as a stand-in — exact on the ledger, approximate on the
#: metabolism (D-19).
#:
#: Propan-1-ol (CAS 71-23-8), C3H8O — the threonine-derived alcohol, shortest of the five.
M_PROPANOL = 3 * _M_C + 8 * _M_H + 1 * _M_O
#: Isobutanol (2-methylpropan-1-ol, CAS 78-83-1), C4H10O — the valine-derived alcohol.
M_ISOBUTANOL = 4 * _M_C + 10 * _M_H + 1 * _M_O
#: Active amyl alcohol (2-methylbutan-1-ol, **CAS 137-32-6**), C5H12O — the ISOLEUCINE-derived
#: alcohol. An isomer of isoamyl alcohol, NOT a synonym for it (see the CAS warning above).
M_ACTIVE_AMYL_OH = 5 * _M_C + 12 * _M_H + 1 * _M_O
#: Isoamyl alcohol (3-methylbutan-1-ol, **CAS 123-51-3**), C5H12O — the LEUCINE-derived
#: alcohol and the dominant higher alcohol of both media (wine mean ~172 mg/L, n=555 studies).
#: The precursor the D-97 ATF1 coupling makes isoamyl acetate first-order in, and the C5 half
#: of the D-69 hydrolysis split.
M_ISOAMYL_OH = 5 * _M_C + 12 * _M_H + 1 * _M_O
#: 2-Phenylethanol (CAS 60-12-8), C8H10O — the PHENYLALANINE-derived alcohol and the only
#: aromatic one of the five: the rose/honey note, not a solventy one. Its wine threshold
#: (~10 mg/L, Guth 1997) sits ~3× BELOW its typical wine concentration (~28.7 mg/L, n=684),
#: which is what gives wine a floral axis from fermentation alone (D-99).
M_2_PHENYLETHANOL = 8 * _M_C + 10 * _M_H + 1 * _M_O
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
#: Methanethiol (methyl mercaptan), CH3SH ≡ CH4S — the molecule the **``methanethiol``** pool IS
#: (decisions D-45/D-110). Methanethiol is the dominant "reduction" thiol (cooked-cabbage /
#: rotten, sensory threshold ~2–3 µg/L; its onion/rubber sibling, ethanethiol, is real chemistry
#: the model does NOT carry). Through D-109 the pool was the plural ``mercaptans`` and this note
#: called methanethiol its "representative species ... the honest single-species stand-in" — the
#: p-coumaric-for-``hydroxycinnamics`` idiom (D-40). **That was the wrong idiom**: a stand-in
#: names a mixture it approximates, and there was no mixture — no Process makes any thiol but this
#: one. D-110 renamed the pool to the compound, so this is not a representative of anything; it is
#: the pool's own molecule. **Unlike H₂S it carries
#: carbon** (one C), so — registered with 1 carbon below — the ``methanethiol`` pool sits on
#: ``total_carbon`` (contrast the carbon-free ``h2s``/SO₂). Nitrogen-free: methanethiol has no N,
#: so the autolytic-thiol Process (:class:`~fermentation.core.kinetics.mercaptans.\
#: AutolyticMercaptan`) **deaminates** the **methionine** nitrogen it draws back to the ``N`` pool
#: (the D-33 idiom; methionine since D-100). Carbon fraction 12.011/48.107 ≈ 0.2497.
M_METHANETHIOL = 1 * _M_C + 4 * _M_H + 1 * _M_S
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
#: The D-100 speciation of the lumped ``amino_acids`` pool into single-molecule amino acids.
#: The lump (arginine stand-in, D-32) could not serve the now-speciated-scale consumers drawing on
#: it (D-96/D-99 fusels + esters at honest levels): two speciated consumers cannot share one lumped
#: substrate (the D-99 done-call finding). Each amino acid below is a real molecule, carbon- AND
#: nitrogen-weighted by its own formula, so a per-precursor consumer (Ehrlich fusel, thermal/
#: oxidative Strecker, mercaptan) draws exactly the amino acid its chemistry names — the whole point
#: of the split. The five Ehrlich precursors map 1:1 onto the D-99 fusels (leucine→isoamyl,
#: isoleucine→active-amyl, valine→isobutanol, threonine→propanol, phenylalanine→2-phenylethanol) and
#: onto the D-87 thermal Strecker aldehydes (isoleucine→2-methylbutanal, leucine→3-methylbutanal,
#: valine→2-methylpropanal, threonine→sotolon via α-ketobutyrate); phenylalanine + methionine feed
#: the D-75 oxidative Strecker aldehydes (phenylacetaldehyde / methional) and methionine the D-45
#: mercaptan. All are N-rich (mass C:N well below biomass ≈ 4.3), so — like arginine — the D-32
#: yeast-swap carbon refund stays below growth's demand and never creates hexose (no clamp needed).
M_LEUCINE = 6 * _M_C + 13 * _M_H + 1 * _M_N + 2 * _M_O
#: L-isoleucine, C6H13NO2 — an isomer of leucine (identical formula, identical C/N weights); the two
#: differ only in odour-precursor identity (active-amyl vs isoamyl alcohol; 2- vs 3-methylbutanal).
M_ISOLEUCINE = 6 * _M_C + 13 * _M_H + 1 * _M_N + 2 * _M_O
M_VALINE = 5 * _M_C + 11 * _M_H + 1 * _M_N + 2 * _M_O
M_THREONINE = 4 * _M_C + 9 * _M_H + 1 * _M_N + 3 * _M_O
M_PHENYLALANINE = 9 * _M_C + 11 * _M_H + 1 * _M_N + 2 * _M_O
#: L-methionine, C5H11NO2S — sulfur-bearing (its thioether is the sulfur source of both methional,
#: D-75, and the autolytic mercaptan, D-45); the S is not tracked on any ledger (the D-45 idiom),
#: only its carbon and nitrogen are.
M_METHIONINE = 5 * _M_C + 11 * _M_H + 1 * _M_N + 2 * _M_O + 1 * _M_S
#: L-glutamine, C5H10N2O3 — the representative species for the ``amino_acids_generic`` bucket: every
#: assimilable amino acid without its own slot (glutamic acid / alanine / serine / aspartate /
#: histidine / lysine / glycine / tryptophan / tyrosine / cysteine / GABA), which the identity-
#: agnostic consumers (the D-32 yeast swap, D-38 MLF growth, D-40 Brett growth) draw alongside
#: arginine so arginine is not the sole generic source (D-100). Glutamine is the canonical
#: assimilable-N proxy and N-rich (mass C:N ≈ 2.14 < biomass ≈ 4.3), keeping the D-32 no-hexose
#: guarantee. **Proline is deliberately NOT tracked**: it is not assimilated anaerobically (excluded
#: from YAN by definition), so it is inert and absent from the pool — the ~48 %-of-must-AA proline
#: fraction simply never enters the model (the dose is an *assimilable* amino-acid dose).
M_GLUTAMINE = 5 * _M_C + 10 * _M_H + 2 * _M_N + 3 * _M_O
#: β-glucan / mannoprotein repeat unit, anhydroglucose C6H10O5 (glucose minus one water, the
#: polysaccharide monomer) — the representative species for the non-assimilable cell-wall
#: **debris** pool yeast autolysis leaves behind (decision D-34). Dead-cell biomass is C-rich
#: (mass C:N ≈ 4–11) while the assimilable amino acids it releases are N-rich (arginine C:N ≈
#: 1.29), so **most of the dead-cell carbon cannot leave as amino acids** — it stays as
#: cell-wall glucan/mannoprotein (the sur-lie lees). :class:`~fermentation.core.kinetics.\
#: autolysis.YeastAutolysis` routes that excess carbon here so autolysis conserves carbon *and*
#: nitrogen separately (nitrogen → ``amino_acids``, the C-rich remainder → ``debris``). Booked
#: **carbon-only** (registered with 0 nitrogen below): all the released nitrogen goes to the
#: amino-acid pool, so the glucan remainder is nitrogen-free — a documented simplification (real
#: mannoproteins retain some protein nitrogen). Carbon fraction 72.066/162.141 ≈ 0.4445.
M_GLUCAN = 6 * _M_C + 10 * _M_H + 5 * _M_O
#: p-Coumaric acid (4-hydroxycinnamic acid), C9H8O3 — the representative species for the lumped
#: **hydroxycinnamics** must-precursor pool (decision D-40, Brett volatile phenols). Grape must
#: carries free hydroxycinnamic acids (p-coumaric + ferulic, ~10–200 mg/L); *Brettanomyces*
#: (and POF+ *S. cerevisiae*) decarboxylate them to 4-vinylphenol/4-vinylguaiacol, which Brett
#: then reduces to the 4-ethylphenol/4-ethylguaiacol "barnyard"/"clove" off-aromas. The two
#: precursors and their two product chains are **lumped** (fork-2 choice, D-40): p-coumaric is the
#: honest single-species stand-in (the dominant 4-EP branch), the
#: succinic-for-``Byp`` idiom (D-16/D-32). Decarboxylation p-coumaric (9 C) → 4-vinylphenol (8 C)
#: + CO2 (1 C) is carbon-closing mole-for-mole on the existing ledger (9 = 8 + 1), exactly like
#: malic → lactic + CO2 (D-23). BOOKKEEPING CAVEAT: lumping p-coumaric for the (larger) ferulic
#: slightly under-counts the 4-EG branch's carbon; exact on the ledger, approximate on provenance.
M_P_COUMARIC = 9 * _M_C + 8 * _M_H + 3 * _M_O
#: 4-Vinylphenol, C8H8O — the representative species for the lumped **vinylphenols** shared
#: intermediate pool (4-vinylphenol + 4-vinylguaiacol; decision D-40). The decarboxylase product
#: and the reductase substrate: a *shared reservoir* POF+ yeast fills but cannot clear (it lacks
#: the reductase) and *Brettanomyces* drains — the emergent yeast/Brett coupling (the α-acetolactate
#: reservoir parallel, D-26/D-31). Reduction 4-vinylphenol (8 C) → 4-ethylphenol (8 C) is a
#: mole-for-mole C8 → C8 transfer (two H added from NADH, no carbon change), carbon-neutral between
#: two weighted pools like diacetyl → butanediol (D-26). Itself a medicinal off-aroma readout.
M_VINYLPHENOL = 8 * _M_C + 8 * _M_H + 1 * _M_O
#: 4-Ethylphenol, C8H10O — the representative species for the lumped **ethylphenols** end-product
#: pool (4-ethylphenol "horse-sweat/barnyard" + 4-ethylguaiacol "clove/smoky"; decision D-40). The
#: terminal Brett volatile-phenol readout, produced-only. Same 8 carbons as its 4-vinylphenol
#: precursor, so the vinylphenol → ethylphenol reduction is carbon-conserving (C8 → C8) and
#: ``total_carbon`` closes to machine precision through the whole precursor → intermediate →
#: product chain.
M_ETHYLPHENOL = 8 * _M_C + 10 * _M_H + 1 * _M_O
#: Ferulic acid (4-hydroxy-3-methoxycinnamic acid), C10H10O4 — the second Brett volatile-phenol
#: precursor, split out from the ``hydroxycinnamics`` (p-coumaric) lump (decision D-55, the D-40
#: pt4 "vinylguaiacol/vinylphenol split" deferral). Decarboxylation ferulic (10 C) → 4-vinylguaiacol
#: (9 C) + CO2 (1 C) closes carbon mole-for-mole (10 = 9 + 1), the same idiom as p-coumaric → 4-
#: vinylphenol + CO2 (9 = 8 + 1). A genuinely distinct molecule (not a fixed-fraction split of the
#: p-coumaric pool) — required because ferulic carries a different carbon count, so no fixed ratio
#: on the existing ``hydroxycinnamics`` output could stay carbon-exact.
M_FERULIC = 10 * _M_C + 10 * _M_H + 4 * _M_O
#: 4-Vinylguaiacol (2-methoxy-4-vinylphenol), C9H10O2 — ferulic acid's decarboxylation product, the
#: ferulic-branch counterpart to ``vinylphenol`` (decision D-55). Reduced by the same Brett
#: vinylphenol reductase (Tchobanov et al. 2008 confirm the enzyme acts on both 4-vinylguaiacol and
#: 4-vinylphenol) to 4-ethylguaiacol, a mole-for-mole C9 → C9 transfer (two H added from NADH, no
#: carbon change) — the "clove/smoky" volatile-phenol readout.
M_VINYLGUAIACOL = 9 * _M_C + 10 * _M_H + 2 * _M_O
#: 4-Ethylguaiacol (4-ethyl-2-methoxyphenol), C9H12O2 — the ferulic-branch terminal Brett
#: volatile-phenol product (decision D-55), the "clove/smoky" counterpart to ``ethylphenol``'s
#: "horse-sweat/barnyard". Same 9 carbons as its 4-vinylguaiacol precursor, so the reduction step is
#: carbon-conserving (C9 → C9) and ``total_carbon`` closes to machine precision.
M_ETHYLGUAIACOL = 9 * _M_C + 12 * _M_H + 2 * _M_O
#: Methional (3-(methylthio)propanal), C4H8OS — the Strecker aldehyde of **methionine**, the
#: potent "cooked-potato" oxidative off-note and one of the principal markers of oxidatively-aged
#: wine and stale beer (decision D-75). Formed when the o-quinones of phenol autoxidation
#: (the ``PhenolicBrowning`` cascade) oxidatively deaminate + decarboxylate methionine: the amino
#: acid loses its carboxyl carbon as CO₂, so the aldehyde carries **4** carbons to methionine's 5.
#: Sulfur-bearing (from methionine's thioether) like ``methanethiol``, and **nitrogen-free** — the
#: Strecker Process deaminates the drawn amino-acid nitrogen back to the ``N`` pool, exactly as
#: :class:`~fermentation.core.kinetics.mercaptans.AutolyticMercaptan` does (D-45).
M_METHIONAL = 4 * _M_C + 8 * _M_H + 1 * _M_O + 1 * _M_S
#: Phenylacetaldehyde (2-phenylacetaldehyde), C8H8O — the Strecker aldehyde of **phenylalanine**,
#: the "honey/floral" note of aged white and dessert wines (decision D-75). The pleasant-valence
#: counterpart to methional's off-note, produced by the same quinone-driven Strecker degradation
#: (phenylalanine's 9 carbons → an 8-carbon aldehyde + CO₂). Nitrogen-free (the amino-acid nitrogen
#: is deaminated to ``N``), so it sits on ``total_carbon`` but not ``total_nitrogen``.
M_PHENYLACETALDEHYDE = 8 * _M_C + 8 * _M_H + 1 * _M_O
#: Ethylidene bridge (—CH(CH₃)—), C2H4 — the representative species for the ``ethyl_bridge`` pool,
#: the acetaldehyde-derived linker of the **acetaldehyde-bridged (ethylidene) condensation** route
#: (decision D-80, the split-ledger colour beat deferred at D-79). When dissolved-O₂ acetaldehyde
#: (D-71) bridges a grape tannin to an anthocyanin — tannin—CH(CH₃)—anthocyanin — the acetaldehyde
#: **loses its carbonyl oxygen as water** (CH₃CHO + 2 Ar–H → Ar–CH(CH₃)–Ar + H₂O), so the retained
#: fragment is the two-carbon ethylidene C2H4. Registered (below) so ``carbon_mass_fraction``
#: weights
#: the ``ethyl_bridge`` slot: acetaldehyde's carbon (borrowed from ethanol at D-71, so **on** the
#: carbon ledger) does NOT vanish into the off-ledger grape-phenolic pigment — it is captured here,
#: on-ledger, and ``total_carbon`` closes to machine precision (the EsterHydrolysis carbon-exact
#: split: release at ``cf(acetaldehyde)``, redeposit at ``cf(ethylidene)`` — this is the "split
#: ledger", the grape bulk off-ledger, the acetaldehyde-derived bridge on it). Same 2 carbons as
#: acetaldehyde; the lost water O is the standing aging-axis mass gap (``total_mass`` weights only
#: ``{S, E, CO2}``, never asserted on an aging run — the D-71 E→acetaldehyde scope-out). Carbon
#: fraction 24.022/28.054 ≈ 0.8563. Nitrogen-free.
M_ETHYLIDENE = 2 * _M_C + 4 * _M_H
#: The non-oxidative THERMAL Strecker aldehydes (decision D-87) — the sweet-wine / Madeira /
#: baked-wine aroma suite the sugar-dicarbonyl (Maillard) route degrades amino acids into WITHOUT
#: O₂, the thermal mirror of the D-75 oxidative (quinone-driven) set. Three are branched-chain
#: Strecker aldehydes of the branched-chain amino acids, each one carbon short of its amino acid
#: (the
#: carboxyl lost as CO₂, exactly as methional/phenylacetaldehyde): **2-methylbutanal** (C5H10O, from
#: isoleucine, "malty/almond"), **3-methylbutanal** (C5H10O, from leucine, "malty/dark-chocolate")
#: and **2-methylpropanal** (isobutyraldehyde, C4H8O, from valine, "malty/grainy"). Nitrogen-free —
#: the amino-acid nitrogen is deaminated back to ``N`` (the D-45/D-75 idiom) — but carbon-bearing,
#: so
#: (like methional/phenylacetaldehyde) they sit on ``total_carbon`` and must be weighted below.
M_2_METHYLBUTANAL = 5 * _M_C + 10 * _M_H + 1 * _M_O
M_3_METHYLBUTANAL = 5 * _M_C + 10 * _M_H + 1 * _M_O
M_2_METHYLPROPANAL = 4 * _M_C + 8 * _M_H + 1 * _M_O
#: Sotolon (4,5-dimethyl-3-hydroxy-2(5H)-furanone), C6H8O3 — the "curry/fenugreek/maple/nutty"
#: marker of botrytized sweet wines (Sauternes), vin jaune, aged Port and Madeira (decision D-87).
#: Unlike the four Strecker aldehydes it is **NOT** a decarboxylation product: its true formation is
#: an aldol condensation of α-ketobutyrate (from threonine deamination) with acetaldehyde (and an
#: ascorbate/sugar-degradation branch), so it carries **no** CO₂ term. Its carbon is drawn from
#: **threonine** — the amino acid its α-ketobutyrate really comes from (decision D-100; before that
#: it was the ``amino_acids`` arginine lump, exact on the ledger but approximate on provenance).
#: Still approximate in one respect the split did not fix: 2 of its 6 carbons are really
#: acetaldehyde-derived, and the acetaldehyde-coupled sotolon route is deferred.
#: Nitrogen-free (deaminated), carbon-tracked.
M_SOTOLON = 6 * _M_C + 8 * _M_H + 3 * _M_O
#: Melanoidin — the brown thermal-browning polymer of :class:`~fermentation.core.kinetics.aging.\
#: Caramelization` (decision D-88), the non-oxidative sugar-only browning route (the O₂-independent
#: thermal mirror of ``PhenolicBrowning``, D-74). A heterogeneous caramelization polymer with no
#: clean molar mass; booked at a **caramelan stand-in** ``C12H18O9`` (two glucose − 3 water, the
#: canonical thermal-dehydration unit) so it carries a plausible carbon fraction (~0.47). It is a
#: **carbon-park** pool (the ``debris``/``glucan`` precedent): unlike the off-ledger oak/hop lumps,
#: it is ON ``total_carbon``, because :class:`Caramelization` *consumes core ``S``* to form it — so
#: the sugar carbon must land in a weighted pool or the transfer would read as carbon destroyed. The
#: transfer is carbon-exact (release at the sugar's fraction, redeposit at this one, the
#: ``EsterHydrolysis`` split idiom); the water lost on dehydration is the standing aging-axis mass
#: gap (``total_mass`` weights only ``{S, E, CO2}``, never asserted on an aging run). CO₂/volatile
#: evolution of real caramelization is lumped into this polymer (a documented v1 simplification).
#: Nitrogen-free — this is CARAMELIZATION (sugar-only), not amino-acid-incorporating Maillard
#: melanoidin (:data:`M_MAILLARD_MELANOIDIN` — the N-bearing route, D-89).
M_MELANOIDIN = 12 * _M_C + 18 * _M_H + 9 * _M_O
#: N-bearing Maillard melanoidin — the brown polymer of
#: :class:`~fermentation.core.kinetics.aging.MaillardBrowning` (decision D-89), the
#: **amino-acid-incorporating** thermal-browning route that D-88's sugar-only
#: :class:`~fermentation.core.kinetics.aging.Caramelization` deferred. Where caramelan
#: (:data:`M_MELANOIDIN`) is nitrogen-free, a true Maillard melanoidin **retains the amino-acid
#: nitrogen in the polymer** (that is what makes it *nitrogenous*), so this species is the FIRST
#: non-biomass, non-arginine pool on ``total_nitrogen``. Booked at a **glucose–glycine
#: model-melanoidin stand-in** ``C8H12O5N`` (a hexose + amino acid condensed and dehydrated — the
#: canonical glucose/glycine Maillard model system; Cämmerer & Kroh), a heterogeneous polymer with
#: no clean molar mass, giving molar C:N ≈ 8:1 and elemental ~47.5 % C / 6.9 % N / 39.6 % O —
#: squarely in the reported melanoidin range. Like caramelan it is an on-ledger **carbon-park** (the
#: ``debris``/``glucan`` precedent) — :class:`MaillardBrowning` consumes core ``S`` **and**
#: ``amino_acids`` to form it, so both the sugar carbon and the amino-acid carbon+nitrogen must
#: land in a weighted pool or the transfer would read as carbon/nitrogen destroyed. The draws are
#: sized to the melanoidin formed (its fixed C:N), so ``total_carbon`` AND ``total_nitrogen`` close
#: to machine precision (§ its class docstring); ALL drawn amino-acid nitrogen is retained here (the
#: deaminating branch is D-87's job). Water lost on dehydration is the standing aging-axis mass gap
#: (``total_mass`` = {S,E,CO2}, never asserted on an aging run). Tag speculative.
M_MAILLARD_MELANOIDIN = 8 * _M_C + 12 * _M_H + 5 * _M_O + 1 * _M_N

#: Molar mass [g/mol] keyed by species name. ``fermentation.core.media`` sugar
#: component names ("glucose", "maltose", "maltotriose") are keys here.
MOLAR_MASS: dict[str, float] = {
    "glucose": M_GLUCOSE,
    "fructose": M_GLUCOSE,
    "maltose": M_MALTOSE,
    "maltotriose": M_MALTOTRIOSE,
    "ethanol": M_ETHANOL,
    "acetaldehyde": M_ACETALDEHYDE,
    "pyruvate": M_PYRUVATE,
    "alpha_ketoglutarate": M_ALPHA_KETOGLUTARATE,
    "alpha_ketobutyrate": M_ALPHA_KETOBUTYRATE,
    "5_oxofructose": M_5_OXOFRUCTOSE,
    "CO2": M_CO2,
    "glycerol": M_GLYCEROL,
    "succinic_acid": M_SUCCINIC,
    "ethyl_acetate": M_ETHYL_ACETATE,
    "isoamyl_acetate": M_ISOAMYL_ACETATE,
    "ethyl_hexanoate": M_ETHYL_HEXANOATE,
    "propanol": M_PROPANOL,
    "isobutanol": M_ISOBUTANOL,
    "active_amyl_alcohol": M_ACTIVE_AMYL_OH,
    "isoamyl_alcohol": M_ISOAMYL_OH,
    "2_phenylethanol": M_2_PHENYLETHANOL,
    "tartaric_acid": M_TARTARIC,
    "malic_acid": M_MALIC,
    "lactic_acid": M_LACTIC,
    "citric_acid": M_CITRIC,
    "sulfur_dioxide": M_SO2,
    "hydrogen_sulfide": M_H2S,
    "methanethiol": M_METHANETHIOL,
    "alpha_acetolactate": M_ACETOLACTATE,
    "diacetyl": M_DIACETYL,
    "butanediol": M_BUTANEDIOL,
    "arginine": M_ARGININE,
    "leucine": M_LEUCINE,
    "isoleucine": M_ISOLEUCINE,
    "valine": M_VALINE,
    "threonine": M_THREONINE,
    "phenylalanine": M_PHENYLALANINE,
    "methionine": M_METHIONINE,
    "glutamine": M_GLUTAMINE,
    "glucan": M_GLUCAN,
    "p_coumaric_acid": M_P_COUMARIC,
    "vinylphenol": M_VINYLPHENOL,
    "ethylphenol": M_ETHYLPHENOL,
    "ferulic_acid": M_FERULIC,
    "vinylguaiacol": M_VINYLGUAIACOL,
    "ethylguaiacol": M_ETHYLGUAIACOL,
    "methional": M_METHIONAL,
    "phenylacetaldehyde": M_PHENYLACETALDEHYDE,
    "ethylidene": M_ETHYLIDENE,
    "2_methylbutanal": M_2_METHYLBUTANAL,
    "3_methylbutanal": M_3_METHYLBUTANAL,
    "2_methylpropanal": M_2_METHYLPROPANAL,
    "sotolon": M_SOTOLON,
    "melanoidin": M_MELANOIDIN,
    "maillard_melanoidin": M_MAILLARD_MELANOIDIN,
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
    #: Pyruvate (C3H4O3) carries three carbons — the excreted overflow-pyruvate SO₂-binder pool
    #: sits on ``total_carbon`` (decision D-49); its re-assimilation is a carbon-closing
    #: C3 → C2 (ethanol) + C1 (CO2) step.
    "pyruvate": 3,
    #: α-Ketoglutarate (C5H6O5) carries five carbons — the second excreted keto-acid SO₂-binder
    #: (decision D-50); its reassimilation is carbon-closing at the Gay-Lussac 2:1 split
    #: (5/3 mol ethanol + 5/3 mol CO2 per mole, C5 → C(10/3) ethanol-carbon + C(5/3) CO2-carbon).
    "alpha_ketoglutarate": 5,
    #: α-Ketobutyrate (C4H6O3) carries four carbons — the third excreted keto-acid (decision
    #: D-107). Two consumers read this count and **must** agree with it or their ledgers break:
    #: the mercaptan's ``5 = 1 + 4`` split (methionine → methanethiol + this) and sotolon's
    #: ``6 = 4 + 2`` aldol (this + acetaldehyde). Both are exact *because* the atom counts are.
    "alpha_ketobutyrate": 4,
    #: 5-Oxofructose (C6H10O6) carries six carbons — the botrytis-specific SO₂-binder must input
    #: (decision D-130). Weighted here for correctness of the absolute carbon total; as an INERT
    #: pool (no Process touches it) it is a constant term that drifts 0, so total_carbon closes
    #: exactly as before — the citrate/tartaric dosed-input precedent. Its grape-fructose carbon
    #: originates OUTSIDE the model (Botrytis oxidised it pre-crush), so it never enters or leaves
    #: the sugar→ethanol ledger; it only occupies the count so the sum is complete.
    "5_oxofructose": 6,
    "CO2": 1,
    "glycerol": 3,
    "succinic_acid": 4,
    "ethyl_acetate": 4,
    #: Isoamyl acetate (C7H14O2) carries SEVEN carbons and ethyl hexanoate (C8H16O2) EIGHT —
    #: each ester pool is weighted by its OWN molecule (decision D-96), so a liquid→gas strip
    #: or the D-69 hydrolysis moves that ester's real carbon, not a stand-in's.
    "isoamyl_acetate": 7,
    "ethyl_hexanoate": 8,
    #: The five Ehrlich alcohols each carry their OWN carbon count (decision D-99), which is
    #: the whole point of the split: the lumped pool weighted all of them as isoamyl's five.
    #: Per MOLE that is a 5-vs-3 error on propan-1-ol, but the ledger weights per GRAM, and
    #: the carbon *mass* fractions are far closer (0.5996 / 0.6482 / 0.6813 / 0.6813 / 0.7865)
    #: because the heavier molecules carry proportionally more hydrogen too — so the lump
    #: over-booked propan-1-ol's carbon by ~14% and under-booked 2-phenylethanol's by ~13%,
    #: not by the 67%/38% the atom counts alone suggest. The two C5 isomers agree exactly here:
    #: they differ in odour, not in carbon.
    "propanol": 3,
    "isobutanol": 4,
    "active_amyl_alcohol": 5,
    "isoamyl_alcohol": 5,
    "2_phenylethanol": 8,
    "tartaric_acid": 4,
    "malic_acid": 4,
    "lactic_acid": 3,
    "citric_acid": 6,
    "sulfur_dioxide": 0,
    "hydrogen_sulfide": 0,
    #: Methanethiol (CH3SH) carries ONE carbon — unlike the carbon-free H₂S/SO₂ — so the lumped
    #: ``mercaptans`` pool sits on ``total_carbon`` (decision D-45).
    "methanethiol": 1,
    "alpha_acetolactate": 5,
    "diacetyl": 4,
    "butanediol": 4,
    "arginine": 6,
    #: The D-100 speciated amino acids, each carrying its own carbon count (decision D-100):
    #: leucine/isoleucine C6 (isomers), valine C5, threonine C4, phenylalanine C9, methionine C5,
    #: glutamine C5 (the ``amino_acids_generic`` bucket representative). Each per-precursor consumer
    #: draws its own molecule, so the carbon it books is that amino acid's, not the arginine lump's.
    "leucine": 6,
    "isoleucine": 6,
    "valine": 5,
    "threonine": 4,
    "phenylalanine": 9,
    "methionine": 5,
    "glutamine": 5,
    "glucan": 6,
    "p_coumaric_acid": 9,
    "vinylphenol": 8,
    "ethylphenol": 8,
    #: Ferulic (10 C) -> vinylguaiacol (9 C) + CO2 (1 C); vinylguaiacol -> ethylguaiacol is a
    #: mole-for-mole C9 -> C9 transfer (decision D-55, the ferulic-acid branch of the Brett/POF
    #: volatile-phenol pathway, split out from the p-coumaric-only ``hydroxycinnamics`` lump).
    "ferulic_acid": 10,
    "vinylguaiacol": 9,
    "ethylguaiacol": 9,
    #: Methional (C4H8OS) carries FOUR carbons — methionine (5 C) loses its carboxyl carbon as CO₂
    #: in the Strecker decarboxylation, so the aldehyde is C4 (decision D-75). Sulfur-bearing but
    #: carbon-tracked (the ``methanethiol`` idiom): the ``methional`` pool sits on ``total_carbon``.
    "methional": 4,
    #: Phenylacetaldehyde (C8H8O) carries EIGHT carbons — phenylalanine (9 C) → an 8-carbon Strecker
    #: aldehyde + CO₂ (decision D-75), the honey-note counterpart to methional.
    "phenylacetaldehyde": 8,
    #: Ethylidene bridge (C2H4) carries TWO carbons — the same two as its acetaldehyde precursor
    #: (only the carbonyl O leaves, as water), so the acetaldehyde → ethyl_bridge transfer is a
    #: carbon-exact C2 → C2 step and the ``ethyl_bridge`` pool sits on ``total_carbon`` (decision
    #: D-80). This is what keeps the on-ledger acetaldehyde carbon from vanishing into the
    #: off-ledger grape-phenolic pigment — the split-ledger accounting.
    "ethylidene": 2,
    #: The non-oxidative THERMAL Strecker aldehydes (decision D-87). The three branched-chain
    #: aldehydes carry ONE fewer carbon than their branched-chain amino acid (the carboxyl lost as
    #: CO₂, the D-75 decarboxylation): 2-methylbutanal (C5, from isoleucine C6), 3-methylbutanal
    #: (C5, from leucine C6), 2-methylpropanal (C4, from valine C5). Sotolon (C6) is a furanone, NOT
    #: a decarboxylation product (no CO₂ term) — its carbon is booked from THREONINE (D-100). All
    #: carbon-tracked (they sit on ``total_carbon``, the methional/phenylacetaldehyde idiom).
    "2_methylbutanal": 5,
    "3_methylbutanal": 5,
    "2_methylpropanal": 4,
    "sotolon": 6,
    #: Melanoidin (caramelan stand-in C12H18O9) carries TWELVE carbons — the caramelization
    #: carbon-park pool Caramelization (D-88) forms from consumed sugar; ON total_carbon (it holds
    #: core-S carbon), so the sugar → melanoidin transfer closes to machine precision (decision
    #: D-88).
    "melanoidin": 12,
    #: N-bearing Maillard melanoidin (glucose–glycine stand-in C8H12O5N) carries EIGHT carbons —
    #: the carbon-park pool MaillardBrowning (D-89) forms from consumed sugar + amino acid; ON
    #: total_carbon (it holds core-S carbon AND amino-acid carbon), so the transfer closes to
    #: machine precision when the draws are sized to it (decision D-89).
    "maillard_melanoidin": 8,
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
    #: Pyruvate is nitrogen-free (a keto-acid), so the excreted-pyruvate pool is absent from
    #: total_nitrogen (decision D-49).
    "pyruvate": 0,
    #: α-Ketoglutarate is nitrogen-free as tracked here (a keto-acid); its real reassimilation
    #: fate is N-coupled (glutamate synthesis), but that coupling is not modelled in v1 (D-50).
    "alpha_ketoglutarate": 0,
    #: α-Ketobutyrate is nitrogen-free (a keto-acid) — and here that is **load-bearing rather than
    #: incidental** (decision D-107). Both routes that make it are *deaminations* (threonine →
    #: 2-oxobutyrate + NH₃; methionine → methanethiol + 2-oxobutyrate + NH₃), so every gram of
    #: precursor nitrogen must land in ``N``: the keto-acid keeps none of it. That is exactly what
    #: makes ``total_nitrogen`` close through both, and why sotolon — built from this pool — is
    #: the one aging product that releases **no** ammonium.
    "alpha_ketobutyrate": 0,
    #: 5-Oxofructose is nitrogen-free (an oxidised sugar), so the botrytis-carbonyl must-input
    #: pool is absent from total_nitrogen (decision D-130).
    "5_oxofructose": 0,
    "CO2": 0,
    "glycerol": 0,
    "succinic_acid": 0,
    "ethyl_acetate": 0,
    "isoamyl_acetate": 0,
    "ethyl_hexanoate": 0,
    #: All five Ehrlich alcohols are nitrogen-free: the pathway TRANSAMINATES the amino acid
    #: away before decarboxylating it, so the nitrogen leaves with the amino group and the
    #: alcohol that survives carries none. This is exactly why the D-33 reroute books the
    #: amino-acid skeleton's carbon while `N` accounting stays untouched (D-99).
    "propanol": 0,
    "isobutanol": 0,
    "active_amyl_alcohol": 0,
    "isoamyl_alcohol": 0,
    "2_phenylethanol": 0,
    "tartaric_acid": 0,
    "malic_acid": 0,
    "lactic_acid": 0,
    "citric_acid": 0,
    "sulfur_dioxide": 0,
    "hydrogen_sulfide": 0,
    #: Methanethiol is nitrogen-free: the autolytic-thiol Process deaminates the methionine
    #: nitrogen it draws back to the ``N`` pool, so the mercaptan itself carries none (D-45/D-100).
    "methanethiol": 0,
    "alpha_acetolactate": 0,
    "diacetyl": 0,
    "butanediol": 0,
    "arginine": 4,
    #: The D-100 speciated amino acids each carry nitrogen (decision D-100): the seven monoamino
    #: acids one N, glutamine (the generic-bucket stand-in) two. This is what keeps the pool on
    #: total_nitrogen after the split — each per-precursor draw deaminates its own amino acid's
    #: nitrogen back to ``N`` (the D-33/D-45/D-75 idiom), and the generic yeast/MLF/Brett swaps
    #: refund/consume each species' own nitrogen so total_nitrogen closes per species.
    "leucine": 1,
    "isoleucine": 1,
    "valine": 1,
    "threonine": 1,
    "phenylalanine": 1,
    "methionine": 1,
    "glutamine": 2,
    #: Cell-wall debris (glucan) is carbon-only: all dead-cell nitrogen is released as
    #: assimilable amino acids, so the remainder carries none (decision D-34).
    "glucan": 0,
    #: The Brett volatile-phenol species are nitrogen-free (hydroxycinnamic acids and their
    #: vinyl-/ethyl-phenol products carry no nitrogen; decision D-40).
    "p_coumaric_acid": 0,
    "vinylphenol": 0,
    "ethylphenol": 0,
    #: The ferulic-acid branch is likewise nitrogen-free (decision D-55).
    "ferulic_acid": 0,
    "vinylguaiacol": 0,
    "ethylguaiacol": 0,
    #: The Strecker aldehydes are nitrogen-free: the Process deaminates the drawn amino-acid
    #: nitrogen back to the ``N`` pool (the D-45 mercaptan idiom), so methional/phenylacetaldehyde
    #: themselves carry none — absent from ``total_nitrogen`` (decision D-75).
    "methional": 0,
    "phenylacetaldehyde": 0,
    #: The ethylidene bridge (C2H4) is nitrogen-free — the acetaldehyde-derived C2 linker carries no
    #: nitrogen, so the ``ethyl_bridge`` pool is absent from ``total_nitrogen`` (decision D-80).
    "ethylidene": 0,
    #: The non-oxidative thermal Strecker aldehydes / sotolon are nitrogen-free: the Process
    #: deaminates the drawn amino-acid nitrogen back to ``N`` (the D-45/D-75 idiom), so the
    #: aldehydes
    #: themselves carry none — absent from ``total_nitrogen`` (decision D-87).
    "2_methylbutanal": 0,
    "3_methylbutanal": 0,
    "2_methylpropanal": 0,
    "sotolon": 0,
    #: Melanoidin is nitrogen-free here: this is CARAMELIZATION (sugar-only, D-88), not amino-acid-
    #: incorporating Maillard melanoidin — so it carries no nitrogen (absent from total_nitrogen).
    "melanoidin": 0,
    #: N-bearing Maillard melanoidin (glucose–glycine stand-in C8H12O5N) carries ONE nitrogen — the
    #: amino-acid nitrogen RETAINED in the polymer (this is what makes a Maillard melanoidin
    #: nitrogenous, unlike sugar-only caramelan). The FIRST non-biomass, non-arginine species on
    #: total_nitrogen: MaillardBrowning (D-89) draws amino_acids sized so all its nitrogen lands
    #: here, so total_nitrogen closes to machine precision (the D-45/D-75 deamination idiom
    #: INVERTED — the nitrogen is parked in the product, not deaminated back to N).
    "maillard_melanoidin": 1,
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
