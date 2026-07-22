"""Kinetics-agnostic conservation and sanity checks over a trajectory.

The harness knows nothing about specific chemistry. A model supplies a
``quantity_fn(state) -> float`` that should be conserved (e.g. total carbon
across sugar, ethanol, CO2, and biomass), and these helpers assert it stays
constant along the trajectory. Encoding conservation as runtime/test assertions
is how we catch a model that quietly creates or destroys mass.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Protocol

import numpy as np

from fermentation.core.chemistry import (
    carbon_mass_fraction,
    nitrogen_mass_fraction,
    sugar_species,
)
from fermentation.core.kinetics.carbon_routing import ESTER_SPECS, FUSEL_SPECS
from fermentation.core.state import FloatArray, StateSchema

QuantityFn = Callable[[FloatArray], float]


class TrajectoryLike(Protocol):
    """The subset of a trajectory these kinetics-agnostic checks read.

    Structural so both a plain :class:`~fermentation.runtime.integrate.Trajectory` and a
    :class:`~fermentation.runtime.schedule.ScheduledTrajectory` (and an ensemble member
    reconstructed as either) satisfy it — the harness never needs the scheduling extras,
    only the state grid.
    """

    @property
    def schema(self) -> StateSchema: ...
    @property
    def t(self) -> FloatArray: ...
    @property
    def y(self) -> FloatArray: ...
    def series(self, name: str) -> FloatArray: ...


# -- conserved-quantity builders for the real chemistry -----------------------
#
# Each builder closes over a schema and returns a ``QuantityFn`` for the harness
# above. Stoichiometric weights come from ``fermentation.core.chemistry`` (a single
# source of truth shared with the kinetics), so the check cannot disagree with the
# model it audits. What each balance covers, and why, is decision D-8.


def _weighted_sum(weights: FloatArray) -> QuantityFn:
    """A ``QuantityFn`` returning ``weights @ y`` (weights precomputed once)."""

    def quantity(y: FloatArray) -> float:
        return float(weights @ y)

    return quantity


def total_carbon(
    schema: StateSchema, *, biomass_carbon_fraction: float | None = None
) -> QuantityFn:
    """Total carbon [g C / L] over sugars, ethanol, CO2 and biomass.

    Carbon is the rigorous conservation invariant (atoms cannot be created or
    destroyed): the model must route every gram of sugar carbon into ethanol,
    CO2 or biomass. Each ``S`` slot is weighted by its species' carbon fraction
    (single-slot ``S`` = hexose); ``E`` and ``CO2`` by theirs.

    Biomass ``X`` contributes ``biomass_carbon_fraction * X``. That fraction is a
    Parameter the growth Process also reads (decision D-8), so it is *passed in*
    here to keep the check consistent with the kinetics. If the schema has an
    ``X`` variable and no fraction is given, this raises rather than silently
    under-counting carbon (which would report a false conservation violation).
    """
    w = schema.zeros()
    s_slice = schema.slice("S")
    for offset, species in enumerate(sugar_species(schema)):
        w[s_slice.start + offset] = carbon_mass_fraction(species)
    if "E" in schema:
        w[schema.slice("E")] = carbon_mass_fraction("ethanol")
    if "CO2" in schema:
        w[schema.slice("CO2")] = carbon_mass_fraction("CO2")
    # Realised-yield byproduct sinks (decision D-16): the sugar-uptake Process
    # routes the carbon it diverts from ethanol/CO2 into these pools, so they must
    # be weighted or carbon would read as destroyed when glycerol is active. The
    # minor-byproduct lump is carbon-accounted as succinic acid (its representative
    # species). Both share the one stoichiometry source with the kinetics.
    if "Gly" in schema:
        w[schema.slice("Gly")] = carbon_mass_fraction("glycerol")
    if "Byp" in schema:
        w[schema.slice("Byp")] = carbon_mass_fraction("succinic_acid")
    # Aroma byproduct pools (decision D-19): the ester/fusel Processes route their
    # carbon out of sugar, so these pools must be weighted or carbon would read as
    # destroyed when they accumulate. Since D-96 (esters) and D-99 (higher alcohols) EVERY
    # one of these pools books as its own molecule — the same species the Processes draw
    # against, from the one chemistry source of truth, so the check matches the kinetics.
    # No overlap with Byp: under D-19 Byp is organic-acids/polyols only (higher alcohols
    # moved to their own pools), so the former Byp double-count is gone.
    # Each ester books as ITSELF (decision D-96): the three single-molecule pools that
    # replaced the lumped ``esters`` pool are weighted at their own molecule's carbon
    # fraction (ethyl acetate C4, isoamyl acetate C7, ethyl hexanoate C8), from the same
    # ``ESTER_SPECS`` registry the Processes draw against — so the draw and the check can
    # never disagree, and a fourth ester is weighted here automatically.
    #
    # Volatilized esters (decisions D-20/D-96): the EsterVolatilization sink moves carbon
    # from each liquid pool into THAT ester's headspace twin as CO2 strips it, so the twin
    # must be weighted at the *same* molecule's fraction or the liquid→gas transfer would
    # read as carbon destroyed. Pairing each pool with its own twin is exactly why one
    # shared gas pool is impossible now the esters differ — a C7 stripped into a
    # C4-weighted pool would create carbon. This mirrors how evolved ``CO2`` stays counted:
    # the carbon leaves the liquid but not the ledger, so closure holds to machine
    # precision while wine's liquid esters honestly fall with temperature.
    for spec in ESTER_SPECS:
        fraction = carbon_mass_fraction(spec.species)
        if spec.pool in schema:
            w[schema.slice(spec.pool)] = fraction
        if spec.gas_pool in schema:
            w[schema.slice(spec.gas_pool)] = fraction
    # Each higher alcohol books as ITSELF (decision D-99), from the same ``FUSEL_SPECS``
    # registry the producer draws against — the fusel twin of the ester rule above, and the
    # ledger half of the D-99 split. The lumped ``fusels`` pool that stood here weighted all
    # five at isoamyl alcohol's carbon mass fraction (0.6813), which over-booked propan-1-ol's
    # carbon by ~14% (its own is 0.5996) and under-booked 2-phenylethanol's by ~13% (0.7865).
    # Those are per-GRAM errors and are much milder than the per-MOLE 5-vs-3 carbon-atom gap
    # implies — the heavier alcohols carry proportionally more hydrogen.
    #
    # THE LEDGER STILL CLOSED UNDER THAT STAND-IN, and that is the instructive part: the
    # producer drew from S at the very same fraction the pool was weighted by, so the error
    # cancelled exactly. Closure alone could never have caught it — a self-consistent wrong
    # weight is invisible to a conservation test. It took splitting the pool to make each
    # number mean what it says.
    #
    # No gas twins here: higher alcohols are not stripped (contrast D-20's ester volatilization).
    for fusel in FUSEL_SPECS:
        if fusel.pool in schema:
            w[schema.slice(fusel.pool)] = carbon_mass_fraction(fusel.species)
    # Vicinal-diketone (VDK) pools (decision D-26): the diacetyl pathway routes carbon
    # sugar → α-acetolactate → diacetyl + CO2 → 2,3-butanediol, every step balanced on
    # this ledger. α-acetolactate is drawn from sugar (excretion) and booked at its own
    # C5 fraction; the spontaneous decarboxylation (C5 → C4 + CO2) is carbon-closing like
    # malic→lactic+CO2; the yeast reduction (diacetyl → butanediol) is a mole-for-mole
    # C4 → C4 transfer between two weighted pools (like esters → esters_gas). Weighting
    # all three at their own species' carbon fraction keeps total_carbon closed to machine
    # precision through the whole produce-then-reabsorb time course.
    if "acetolactate" in schema:
        w[schema.slice("acetolactate")] = carbon_mass_fraction("alpha_acetolactate")
    if "diacetyl" in schema:
        w[schema.slice("diacetyl")] = carbon_mass_fraction("diacetyl")
    if "butanediol" in schema:
        w[schema.slice("butanediol")] = carbon_mass_fraction("butanediol")
    # Acetaldehyde (decision D-27): a transient buffer on the main pathway. Production
    # borrows carbon from ethanol (E) and reduction returns it, both mole-for-mole C2 → C2,
    # so weighting acetaldehyde at its own carbon fraction keeps that borrow/return
    # carbon-neutral and total_carbon closed to machine precision through the produce-then-
    # reabsorb course. No sugar/CO2 term is involved — acetaldehyde de-lumps the uptake
    # Process's existing sugar→ethanol step rather than adding a parallel pathway.
    if "acetaldehyde" in schema:
        w[schema.slice("acetaldehyde")] = carbon_mass_fraction("acetaldehyde")
    # Excreted overflow pyruvate (decision D-49): PyruvateExcretion draws its carbon out of
    # sugar (C3, booked at pyruvate's fraction) and PyruvateReassimilation returns it to
    # ethanol + CO2 (C3 → C2 + C1, one mole each — carbon-closing like malic → lactic + CO2),
    # so the pool must be weighted here or the excrete-then-reassimilate course would read as
    # carbon created (during excretion) then destroyed (during re-assimilation). Weighted at
    # pyruvate's own carbon fraction — the SAME chemistry source the draw books against — so
    # the carbon into pyruvate equals the carbon out of S and total_carbon closes to machine
    # precision through the whole time course, leaving only the stranded residual as sugar
    # carbon parked as pyruvate. Re-assimilation returns to E/CO2 (not S) deliberately:
    # post-dryness S is 0, so a refund-to-sugar would be a no-op that destroys carbon. On a
    # keto-acid-pool-off run the pool is empty (constant 0 term). Nitrogen-free (a keto-acid),
    # so it is absent from total_nitrogen.
    if "pyruvate" in schema:
        w[schema.slice("pyruvate")] = carbon_mass_fraction("pyruvate")
    # Excreted overflow alpha-ketoglutarate (decision D-50): the SAME excreted-side-pool
    # structure as pyruvate above, weighted here at its own C5 carbon fraction so the
    # excrete-then-reassimilate course (S -> alpha_ketoglutarate -> E + CO2, at the Gay-Lussac
    # 2:1 split — NOT mole-for-mole, see AlphaKetoglutarateReassimilation) closes to machine
    # precision. Nitrogen-free, so absent from total_nitrogen.
    if "alpha_ketoglutarate" in schema:
        w[schema.slice("alpha_ketoglutarate")] = carbon_mass_fraction("alpha_ketoglutarate")
    # Excreted overflow alpha-ketobutyrate — THE KETO-ACID NODE (decision D-107): the same
    # excreted-side-pool structure again, weighted at its own C4 fraction. This pool carries
    # THREE flows across this ledger, not one, and all three close on ATOM COUNTS:
    #   in:  threonine (C4) -> alpha_ketobutyrate (C4)          [deamination, 4 == 4]
    #   in:  methionine (C5) -> methanethiol (C1) + this (C4)   [demethiolation, 5 == 1 + 4]
    #   out: this (C4) + acetaldehyde (C2) -> sotolon (C6)      [the aldol, 4 + 2 == 6]
    # plus the excretion's sugar de-novo share and the Gay-Lussac reassimilation it shares with
    # its two siblings. Weighting it here is what lets all of them close to machine precision —
    # and note that closure is NOT evidence any of them draws the right number of moles (D-105's
    # lesson: a draw sized to close the ledger can never violate it). Nitrogen-free, so absent
    # from total_nitrogen — which is itself load-bearing: both producing routes are deaminations,
    # so every gram of precursor nitrogen must land in N rather than ride along here.
    if "alpha_ketobutyrate" in schema:
        w[schema.slice("alpha_ketobutyrate")] = carbon_mass_fraction("alpha_ketobutyrate")
    # Botrytis 5-oxofructose SO₂-binder (decision D-130): a dosed, carbon-bearing (C6) must input,
    # INERT — no Process touches it, so it is a constant term that drifts 0 and total_carbon closes
    # exactly as before (the citrate/tartaric dosed-input precedent above). Its grape-fructose
    # carbon originates OUTSIDE the sugar→ethanol ledger (Botrytis oxidised it pre-crush), so it
    # never flows across this balance; it is weighted only so the absolute carbon total is complete.
    # Nitrogen-free (an oxidised sugar), so absent from total_nitrogen. Empty on a non-botrytis run
    # (constant 0 term) — byte-for-byte the pre-D-130 core.
    if "oxofructose" in schema:
        w[schema.slice("oxofructose")] = carbon_mass_fraction("5_oxofructose")
    # MCFA MLF-inhibitor pool (decision D-131): the aggregate yeast-secreted octanoic+decanoic
    # acids, carried octanoic-equivalent (C8), INERT — no Process touches it (v1 defers the
    # yeast-synthesis production layer), so it is a constant term that drifts 0 and total_carbon
    # closes exactly as before, like the oxofructose must-input above. Weighted only so the
    # absolute carbon total is complete; nitrogen-free, so absent from total_nitrogen. Empty on a
    # run with no MCFA dose (constant 0 term) — byte-for-byte the pre-D-131 core.
    if "mcfa" in schema:
        w[schema.slice("mcfa")] = carbon_mass_fraction("octanoic_acid")
    # Wine acid slots (decision D-18): the pH charge balance reads these, and the MLF
    # Process (decision D-23, when Oenococcus oeni is pitched) moves carbon malic (C4) ->
    # lactic (C3) + CO2 (C1) — balanced mole-for-mole — so they are weighted here for that
    # conversion to stay carbon-closing on this same ledger (no new conservation code). On
    # an undosed run the acids are inert (no active Process touches them, derivatives 0),
    # so this adds a constant term that drifts 0; ``cation_charge`` is a charge density,
    # not a carbon species, and stays weight 0 (``schema.zeros``). ``Byp`` weighting is
    # unchanged — the charge balance only *reads* it (include-by-reading), adding no carbon,
    # so ``total_carbon`` and the double-count are exactly as before.
    if "tartaric" in schema:
        w[schema.slice("tartaric")] = carbon_mass_fraction("tartaric_acid")
    if "malic" in schema:
        w[schema.slice("malic")] = carbon_mass_fraction("malic_acid")
    if "lactic" in schema:
        w[schema.slice("lactic")] = carbon_mass_fraction("lactic_acid")
    # Citrate (decision D-31): a dosed must input O. oeni co-metabolises during MLF. The
    # MalolacticCitrateMetabolism Process (when O. oeni is pitched) routes it as a lumped
    # citrate (C6) → α-acetolactate (C5) + CO2 (C1) conversion feeding the shared VDK
    # reservoir — carbon-closing mole-for-mole (6 = 5 + 1) on this same ledger, so weighting
    # citrate at its citric-acid fraction keeps total_carbon closed through that conversion.
    # On an undosed run citrate is inert (no active Process touches it, derivative 0), a
    # constant term that drifts 0. Charge-inactive (kept out of the D-18 pH balance in v1),
    # so it is a carbon term only — like malic/lactic, weighted here for the conversion.
    if "citrate" in schema:
        w[schema.slice("citrate")] = carbon_mass_fraction("citric_acid")
    # Amino-acid pool (decision D-32): a dosed, carbon- AND nitrogen-bearing wine pool. The
    # AminoAcidAssimilation swap debits it and refunds the displaced biomass carbon to sugar
    # and the displaced biomass nitrogen to ``N`` — a pure carbon/nitrogen-neutral transfer
    # (aa carbon → S, aa nitrogen → N), so weighting the pool at arginine's carbon fraction
    # (its representative species) keeps total_carbon closed through the swap. On an undosed
    # run the pool is empty and the Process is disabled (constant 0 term). Its NITROGEN side
    # is weighted in total_nitrogen below.
    if "amino_acids" in schema:
        w[schema.slice("amino_acids")] = carbon_mass_fraction("arginine")
    # The D-100 speciated amino-acid pools: the lumped `amino_acids` (arginine) split into seven
    # named amino acids + a generic bucket, each carbon-weighted at its OWN molecule's fraction —
    # the same per-species discipline as the D-96 esters / D-99 fusels. A per-precursor consumer
    # (Ehrlich reroute, thermal/oxidative Strecker, mercaptan) debits its species and deposits the
    # carbon in its aroma product, so weighting each at its own fraction keeps total_carbon closed
    # through that draw. The generic bucket is booked at glutamine (its stand-in). All empty on
    # an undosed / arginine-only run (constant 0 terms) — byte-for-byte the pre-D-100 core.
    for _aa in ("leucine", "isoleucine", "valine", "threonine", "phenylalanine", "methionine"):
        if _aa in schema:
            w[schema.slice(_aa)] = carbon_mass_fraction(_aa)
    if "amino_acids_generic" in schema:
        w[schema.slice("amino_acids_generic")] = carbon_mass_fraction("glutamine")
    # Cell-wall debris pool (decision D-34): yeast autolysis (YeastAutolysis) turns dead biomass
    # X_dead into amino acids + this non-assimilable carbon-rich remainder (booked as glucan). The
    # dead-cell carbon r·f_C splits into the amino acids' carbon and this debris carbon, weighted
    # at glucan's carbon fraction to keep total_carbon closed through autolysis (the excess carbon
    # leaves the metabolite pools but stays on the ledger — the esters_gas idiom). Debris is
    # carbon-only (released nitrogen all goes to amino_acids), so it is absent from
    # total_nitrogen. On an undosed / autolysis-off run the pool is empty (constant 0 term).
    if "debris" in schema:
        w[schema.slice("debris")] = carbon_mass_fraction("glucan")
    # Mercaptan (thiol) pool (decision D-45): AutolyticMercaptan draws carbon from the amino-acid
    # pool into this pool (deaminating the nitrogen to N), so — unlike the carbon-free h2s — it
    # carries carbon and must be weighted or that draw would read as carbon destroyed. Booked at
    # methanethiol's carbon fraction, so total_carbon closes to machine precision through the
    # transfer: the carbon into methanethiol equals the carbon out of methionine. Since D-110 the
    # slot name and the species name are ONE name (the pool was `mercaptans` through D-109, a false
    # plural), so this weighting is now the identity every single-molecule pool shows — the lookup
    # can no longer drift from the slot it weights. Nitrogen-free (the methionine N is deaminated
    # to the N pool), so it is absent from total_nitrogen. On an undosed / autolysis-off run the
    # pool is empty and the Process disabled (constant 0 term). NOTE: the add_copper verb (D-45)
    # removes the thiol as precipitated copper mercaptide — carbon that legitimately LEAVES the
    # wine, booked as a negative external flow (the racking-debris precedent), so the run-wide
    # identity final == initial + Σ flows still holds even though total_carbon(state) drops.
    if "methanethiol" in schema:
        w[schema.slice("methanethiol")] = carbon_mass_fraction("methanethiol")
    # The metal-complexed thiol reservoir (decision D-135) is weighted at the SAME carbon fraction
    # as the free pool it feeds, and that identity is the whole point: BoundMethanethiolRelease is a
    # 1:1 transfer of one molecule between two slots (complexed CH3SH -> free CH3SH — the ligand is
    # the same molecule, only its binding state changes), so equal weights make the release close
    # total_carbon to machine precision. Weighting the bonded form at 0 instead would read as carbon
    # CREATED out of nothing on every release step. Contrast bound_h2s, deliberately absent here:
    # H2S is carbon-free, so its reservoir sits off the ledger exactly as the free h2s pool does.
    if "bound_methanethiol" in schema:
        w[schema.slice("bound_methanethiol")] = carbon_mass_fraction("methanethiol")
    # Brett volatile-phenol pools (decision D-40): the decarboxylase routes carbon hydroxycinnamics
    # (p-coumaric, C9) → vinylphenols (C8) + CO2 (C1), carbon-closing mole-for-mole like malic →
    # lactic + CO2; the reductase moves vinylphenols (C8) → ethylphenols (C8), a mole-for-mole C8 →
    # C8 transfer between two weighted pools like diacetyl → butanediol. Weighting all three at
    # their representative species (p-coumaric / 4-vinylphenol / 4-ethylphenol) keeps total_carbon
    # closed through the whole precursor → intermediate → product chain. On an undosed / un-pitched
    # run the pools are empty and the Processes disabled (constant 0 terms). X_brett/X_brett_dead
    # are weighted as biomass in the X block below (decision D-40 pt2, the X_mlf D-38 precedent).
    if "hydroxycinnamics" in schema:
        w[schema.slice("hydroxycinnamics")] = carbon_mass_fraction("p_coumaric_acid")
    if "vinylphenols" in schema:
        w[schema.slice("vinylphenols")] = carbon_mass_fraction("vinylphenol")
    if "ethylphenols" in schema:
        w[schema.slice("ethylphenols")] = carbon_mass_fraction("ethylphenol")
    # The ferulic-acid branch (decision D-55): a genuinely distinct precursor/intermediate/product
    # chain from the p-coumaric branch above (different carbon counts: 10 C -> 9 C + CO2, vs 9 C ->
    # 8 C + CO2), weighted the same way at its own representative species so total_carbon closes
    # through this chain too. On an undosed run these pools are empty and the Processes disabled.
    if "ferulic_acid" in schema:
        w[schema.slice("ferulic_acid")] = carbon_mass_fraction("ferulic_acid")
    if "vinylguaiacols" in schema:
        w[schema.slice("vinylguaiacols")] = carbon_mass_fraction("vinylguaiacol")
    if "ethylguaiacols" in schema:
        w[schema.slice("ethylguaiacols")] = carbon_mass_fraction("ethylguaiacol")
    # Strecker aldehyde pools (decision D-75): StreckerDegradation draws carbon from the amino-acid
    # pool into methional (C4) + phenylacetaldehyde (C8) and releases the acid's carboxyl carbon
    # as CO2 (deaminating the nitrogen to N), so — like the mercaptans draw — these carry carbon and
    # must be weighted or the transfer would read as carbon destroyed. Booked at each aldehyde's own
    # carbon fraction (its representative species, the same one the Process deposits against), so
    # total_carbon closes to machine precision: the carbon into methional + phenylacetaldehyde + CO2
    # equals the carbon out of each aldehyde's own precursor. Nitrogen-free (that N is deaminated
    # to the N pool), so both are
    # absent from total_nitrogen. On a reductive / amino-acid-free run the pools are empty and the
    # Process contributes zero (constant 0 terms).
    if "methional" in schema:
        w[schema.slice("methional")] = carbon_mass_fraction("methional")
    if "phenylacetaldehyde" in schema:
        w[schema.slice("phenylacetaldehyde")] = carbon_mass_fraction("phenylacetaldehyde")
    # Non-oxidative THERMAL Strecker aldehydes + sotolon (decision D-87): MaillardStrecker draws
    # carbon from the amino-acid pool into the four new aroma pools (three branched-chain aldehydes
    # +
    # sotolon) and — for the five decarboxylating aldehydes only, NOT sotolon — releases the acid's
    # carboxyl carbon as CO2, deaminating the nitrogen to N. Like the D-75 pair they carry carbon
    # and
    # must be weighted or the transfer would read as carbon destroyed; booked at each product's own
    # carbon fraction (the same species the Process deposits against), so total_carbon closes to
    # machine precision: carbon into the aldehydes + sotolon + CO2 equals carbon out of each one's
    # own precursor. Nitrogen-free (that N deaminated to the N pool), so absent from total_nitrogen.
    # Empty (constant 0)
    # on any run where the Process is inert (no amino acids or no residual sugar).
    for _thermal_ald in ("2_methylbutanal", "3_methylbutanal", "2_methylpropanal", "sotolon"):
        if _thermal_ald in schema:
            w[schema.slice(_thermal_ald)] = carbon_mass_fraction(_thermal_ald)
    # Caramelization melanoidin (decision D-88; medium-agnostic D-90 — the guard fires for BOTH
    # schemas): the sugar-only thermal-browning carbon-park. Unlike every other aging-browning pool
    # (o2/A420 off-ledger), Caramelization CONSUMES core S to form melanoidin, so the sugar carbon
    # it draws must land in a weighted pool or the transfer would read as carbon destroyed. Booked
    # at melanoidin's own (caramelan stand-in) carbon fraction — the species the Process deposits
    # against — so total_carbon closes to machine precision: the carbon out of S (on beer, summed
    # over the vector at each sugar's own fraction — the D-90 vectorized draw) equals the
    # carbon into melanoidin (the EsterHydrolysis carbon-exact split). Nitrogen-free (sugar-only,
    # not amino-acid Maillard). Empty (constant 0) on any run where S ≈ 0 at aging.
    if "melanoidin" in schema:
        w[schema.slice("melanoidin")] = carbon_mass_fraction("melanoidin")
    # N-bearing Maillard melanoidin (decision D-89): the amino-acid-incorporating thermal-browning
    # carbon-park. Like caramelan (D-88) it holds core-S carbon, but MaillardBrowning ALSO draws
    # amino_acids, so this pool carries BOTH the sugar carbon and the amino-acid carbon. Booked at
    # its own (glucose–glycine stand-in) carbon fraction — the species the Process deposits against
    # — so total_carbon closes to machine precision: the carbon out of S + amino_acids equals the
    # carbon into maillard_melanoidin (the draws are sized to it; see the class docstring). It is
    # ALSO on total_nitrogen (the nitrogen it retains), the first non-biomass, non-arginine species
    # there. Empty (constant 0) on any run where S ≈ 0 or amino_acids == 0 at aging.
    if "maillard_melanoidin" in schema:
        w[schema.slice("maillard_melanoidin")] = carbon_mass_fraction("maillard_melanoidin")
    # Acetaldehyde-bridged condensation ethyl-bridge pool (decision D-80): the SPLIT-LEDGER capture.
    # AcetaldehydeBridgedCondensation consumes ON-ledger acetaldehyde (whose carbon is borrowed from
    # E at D-71) to form an ethylidene bridge —CH(CH₃)— in an OFF-ledger grape-phenolic pigment. If
    # that acetaldehyde carbon were not re-captured on the ledger it would vanish and fail the
    # carbon
    # balance (the trap D-79 named). So this ``ethyl_bridge`` slot is weighted here at ethylidene's
    # carbon fraction — the acetaldehyde-derived C2 bridge is ON the ledger while the grape bulk
    # (anthocyanin/tannin) stays OFF it (the "split ledger"). The Process re-deposits the consumed
    # acetaldehyde carbon here via the EsterHydrolysis carbon-exact split (release at
    # cf(acetaldehyde), redeposit at cf(ethylidene)), so total_carbon closes to machine precision:
    # acetaldehyde↓ exactly cancels ethyl_bridge↑ (a NON-trivial closure, unlike the direct D-79
    # route which moves nothing conserved). On a run without the bridged route the pool is 0
    # (constant
    # 0 term). Nitrogen-free (absent from total_nitrogen); the lost carbonyl O (→ water) is the
    # standing aging-axis mass gap (total_mass weights only {S,E,CO2}, never asserted on aging).
    if "ethyl_bridge" in schema:
        w[schema.slice("ethyl_bridge")] = carbon_mass_fraction("ethylidene")
    if "X" in schema:
        if biomass_carbon_fraction is None:
            raise ValueError(
                "schema has biomass 'X'; pass biomass_carbon_fraction (the value "
                "the growth Process uses) so the carbon check matches the kinetics"
            )
        w[schema.slice("X")] = biomass_carbon_fraction
        # Inactivated cells are still biomass of the same composition: ethanol
        # inactivation moves mass X -> X_dead, so counting both pools at the same
        # carbon fraction keeps that transfer carbon-neutral (decision D-13).
        if "X_dead" in schema:
            w[schema.slice("X_dead")] = biomass_carbon_fraction
        # Bacterial biomass X_mlf (decision D-38): once MalolacticGrowth is active it
        # builds O. oeni biomass from the amino-acid pool, drawing the carbon shortfall from
        # sugar, so X_mlf carries carbon and must be weighted or that growth would read as
        # carbon destroyed. Booked at the SAME biomass_carbon_fraction the growth stoichiometry
        # draws against (bacterial ≈ yeast elemental composition, a documented v1 simplification),
        # so carbon closes exactly. This supersedes the v1 "X_mlf is carbon-free" scoping (D-23):
        # on a conversion-only run X_mlf only ever moves into X_mlf_dead (death D-39 / senescence
        # D-41 — both to an equally-weighted pool), so the added term stays closure-neutral; a
        # co-inoculation dose / pitch_mlf flow now carries this bacterial-biomass carbon.
        if "X_mlf" in schema:
            w[schema.slice("X_mlf")] = biomass_carbon_fraction
        # Non-viable bacterial biomass X_mlf_dead (decisions D-39/D-41): the SO₂ kill AND the benign
        # senescence baseline both move X_mlf into it, so — like the yeast X → X_dead transfer
        # (D-13) — it must be weighted at the SAME biomass_carbon_fraction or that loss would read
        # as carbon destroyed. (Terminal sink: autolysis reads only yeast X_dead, not X_mlf_dead.)
        if "X_mlf_dead" in schema:
            w[schema.slice("X_mlf_dead")] = biomass_carbon_fraction
        # Brett biomass X_brett (decision D-40 pt2): BrettGrowth builds it from the amino-acid pool,
        # drawing the carbon shortfall from ETHANOL, so X_brett carries carbon and must be weighted
        # or that growth would read as carbon destroyed. Booked at the SAME biomass_carbon_fraction
        # the growth stoichiometry draws against (Brett ≈ yeast elemental composition, a v1
        # simplification), so carbon closes exactly. Supersedes the pt1 "X_brett is carbon-free"
        # scoping: on a no-growth run X_brett is constant, so the added term is a constant offset
        # that drifts 0; a co-inoculation dose / pitch_brett flow now carries this biomass carbon.
        if "X_brett" in schema:
            w[schema.slice("X_brett")] = biomass_carbon_fraction
        # Non-viable Brett biomass X_brett_dead (decision D-40 pt3): BrettDeath moves X_brett into
        # it under SO₂, so — like the X → X_dead / X_mlf → X_mlf_dead transfers — it is weighted at
        # the same fraction so that death reads as carbon-neutral.
        if "X_brett_dead" in schema:
            w[schema.slice("X_brett_dead")] = biomass_carbon_fraction
    return _weighted_sum(w)


def total_nitrogen(
    schema: StateSchema, *, biomass_nitrogen_fraction: float | None = None
) -> QuantityFn:
    """Total nitrogen [g N / L]: free YAN plus nitrogen bound in biomass.

    As cells grow they assimilate ``N`` (yeast-assimilable nitrogen) into biomass,
    so the invariant is ``N + biomass_nitrogen_fraction * X``. Like the biomass
    carbon fraction, ``biomass_nitrogen_fraction`` is a Parameter the growth
    Process reads and is passed in (decision D-8); required iff the schema has
    biomass. Conserved once the nitrogen-limited growth Process exists.
    """
    w = schema.zeros()
    if "N" in schema:
        w[schema.slice("N")] = 1.0
    # Amino-acid pool (decision D-32): the aa pool carries nitrogen (arginine, 4 N per
    # molecule) and the AminoAcidAssimilation swap moves that nitrogen aa → N (ammonium
    # refund) mole-for-mole, so weighting the pool at arginine's nitrogen fraction keeps
    # total_nitrogen closed through the swap. This is the first nitrogen-bearing tracked
    # species besides biomass; on an undosed run the pool is empty (constant 0 term).
    if "amino_acids" in schema:
        w[schema.slice("amino_acids")] = nitrogen_mass_fraction("arginine")
    # The D-100 speciated amino-acid pools carry nitrogen too (each monoamino acid one N, glutamine
    # two), weighted at their own fractions so total_nitrogen closes through every per-species draw:
    # a catabolic consumer deaminates the drawn amino acid's nitrogen back to `N` (the D-33/D-45
    # /D-75 idiom), and the generic yeast/MLF/Brett swaps refund/consume each species' own N. All
    # empty on an undosed / arginine-only run (constant 0 terms) — the pre-D-100 core is unchanged.
    for _aa in ("leucine", "isoleucine", "valine", "threonine", "phenylalanine", "methionine"):
        if _aa in schema:
            w[schema.slice(_aa)] = nitrogen_mass_fraction(_aa)
    if "amino_acids_generic" in schema:
        w[schema.slice("amino_acids_generic")] = nitrogen_mass_fraction("glutamine")
    # N-bearing Maillard melanoidin (decision D-89): the FIRST non-biomass, non-arginine species on
    # the nitrogen ledger. MaillardBrowning RETAINS the amino-acid nitrogen it draws in the
    # melanoidin polymer (the deaminating branch is D-87's job — this is the D-45/D-75 deamination
    # idiom INVERTED: the nitrogen is parked in the product, not refunded to N). Sizing the
    # amino_acids draw so all its nitrogen lands here means the pool loses exactly the nitrogen
    # maillard_melanoidin gains, so total_nitrogen closes to machine precision through the transfer.
    # Weighted at this species' own nitrogen fraction — the species the Process deposits against.
    # Empty (constant 0) on any run where S ≈ 0 or amino_acids == 0 at aging.
    if "maillard_melanoidin" in schema:
        w[schema.slice("maillard_melanoidin")] = nitrogen_mass_fraction("maillard_melanoidin")
    if "X" in schema:
        if biomass_nitrogen_fraction is None:
            raise ValueError(
                "schema has biomass 'X'; pass biomass_nitrogen_fraction (the value "
                "the growth Process uses) so the nitrogen check matches the kinetics"
            )
        w[schema.slice("X")] = biomass_nitrogen_fraction
        # Inactivated biomass retains its nitrogen: count X_dead so the X -> X_dead
        # inactivation transfer stays nitrogen-neutral (decision D-13).
        if "X_dead" in schema:
            w[schema.slice("X_dead")] = biomass_nitrogen_fraction
        # Bacterial biomass X_mlf (decision D-38): MalolacticGrowth builds it by
        # assimilating the amino-acid pool's nitrogen (arginine), so X_mlf carries nitrogen and
        # is weighted at the SAME biomass_nitrogen_fraction the growth draws against — the pool
        # loses exactly the nitrogen X_mlf gains, so total_nitrogen closes. On a conversion-only run
        # X_mlf only moves into X_mlf_dead (death/senescence), so closure holds; a co-inoculation
        # dose / pitch_mlf flow now carries it.
        if "X_mlf" in schema:
            w[schema.slice("X_mlf")] = biomass_nitrogen_fraction
        # Non-viable bacterial biomass X_mlf_dead (decisions D-39/D-41): retains its nitrogen, at
        # the same biomass_nitrogen_fraction so BOTH the X_mlf → X_mlf_dead death AND senescence
        # transfers are nitrogen-neutral (the yeast X → X_dead precedent, D-13).
        if "X_mlf_dead" in schema:
            w[schema.slice("X_mlf_dead")] = biomass_nitrogen_fraction
        # Brett biomass X_brett (decision D-40 pt2): BrettGrowth assimilates the amino-acid pool's
        # nitrogen (arginine) into biomass, so X_brett carries nitrogen, weighted at the same
        # biomass_nitrogen_fraction the growth draws against — the pool loses exactly the nitrogen
        # X_brett gains, so total_nitrogen closes. Constant (⇒ 0 drift) on a no-growth run.
        if "X_brett" in schema:
            w[schema.slice("X_brett")] = biomass_nitrogen_fraction
        # Non-viable Brett biomass X_brett_dead (decision D-40 pt3): retains its nitrogen, counted
        # at the same fraction so the X_brett → X_brett_dead death transfer is nitrogen-neutral.
        if "X_brett_dead" in schema:
            w[schema.slice("X_brett_dead")] = biomass_nitrogen_fraction
    return _weighted_sum(w)


def total_mass(schema: StateSchema) -> QuantityFn:
    """Total mass [g/L] of the abiotic conversion species ``{S, E, CO2}`` — wine only.

    Mass closes only for a *single hexose* (decision D-8): the wine reaction
    ``C6H12O6 -> 2 C2H5OH + 2 CO2`` is mass-balanced (180.156 = 92.138 + 88.018
    g/mol), so ``S + E + CO2`` is conserved to solver tolerance. It does **not**
    generalise, for the same untracked-solvent reason biomass is excluded:
    di-/trisaccharides *hydrolyse*, pulling water from the solvent into the product
    pool (maltose adds ~5.3% mass, maltotriose ~7.1%), so ``S + E + CO2`` is not a
    beer invariant. Carbon carries no such term (water has no carbon) and is the
    rigorous cross-medium invariant — so this builder **rejects a multi-component
    sugar** and beer relies on :func:`total_carbon`.

    It stays scoped to ``{S, E, CO2}`` for one further reason (decision D-16): the
    realised-yield byproduct sinks ``Gly``/``Byp`` are more reduced than the
    sugar→ethanol+CO2 route and draw redox H/O from the solvent like biomass, so
    ``{S, E, CO2}`` mass does not close once that diversion is active. This is thus
    the **validated-core** (byproduct-free) mass check — assert it only on a
    glycerol-off configuration; ``total_carbon`` (which weights ``Gly``/``Byp``)
    is the invariant when byproducts are on.
    """
    if "S" in schema and schema.spec("S").size > 1:
        raise ValueError(
            "total_mass is a hexose/wine check: a multi-component sugar hydrolyses "
            "(pulling solvent water into S+E+CO2), so mass does not close. Use "
            "total_carbon, which is conserved across media (decision D-8)."
        )
    w = schema.zeros()
    for name in ("S", "E", "CO2"):
        if name in schema:
            w[schema.slice(name)] = 1.0
    return _weighted_sum(w)


def _evaluate(traj: TrajectoryLike, quantity_fn: QuantityFn) -> FloatArray:
    return np.array([quantity_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])


def max_drift(traj: TrajectoryLike, quantity_fn: QuantityFn) -> float:
    """Maximum absolute deviation of ``quantity_fn`` from its initial value."""
    q = _evaluate(traj, quantity_fn)
    return float(np.max(np.abs(q - q[0]))) if q.size else 0.0


def assert_conserved(
    traj: TrajectoryLike,
    quantity_fn: QuantityFn,
    *,
    rtol: float = 1e-6,
    atol: float = 1e-9,
    label: str = "quantity",
) -> None:
    """Assert ``quantity_fn`` stays within tolerance of its initial value.

    Tolerance is relative to the initial magnitude plus an absolute floor, so it
    behaves sensibly when the conserved quantity is near zero. Raises
    ``AssertionError`` (so it reads naturally inside tests and runtime checks).
    """
    q = _evaluate(traj, quantity_fn)
    if q.size == 0:
        return
    tol = atol + rtol * abs(q[0])
    drift = np.max(np.abs(q - q[0]))
    if drift > tol:
        worst = int(np.argmax(np.abs(q - q[0])))
        raise AssertionError(
            f"{label} not conserved: drift {drift:.3e} > tol {tol:.3e} "
            f"(initial {q[0]:.6g}, at t={traj.t[worst]:.3g}h -> {q[worst]:.6g})"
        )


def assert_nonnegative(
    traj: TrajectoryLike, variables: tuple[str, ...], *, atol: float = 1e-9
) -> None:
    """Assert the named variables never go meaningfully negative.

    Concentrations and biomass are physical and must stay >= 0; a small negative
    excursion within ``atol`` is tolerated as solver noise.
    """
    for name in variables:
        series = traj.series(name)
        worst = float(np.min(series))
        if worst < -atol:
            idx = int(np.argmin(series))
            raise AssertionError(f"{name} went negative: {worst:.3e} at t={traj.t[idx]:.3g}h")
