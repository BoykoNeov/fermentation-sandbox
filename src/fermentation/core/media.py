"""Medium definitions — the named state layouts the validated core models.

A *medium* (wine, beer, …) fixes two things the rest of the engine builds on:

  * its :class:`~fermentation.core.state.StateSchema` — how many sugar slots, in
    what order, alongside biomass / ethanol / nitrogen / temperature / CO2; and
  * the Processes that act on that state (the kinetics).

Both are data here, not physics. This module declares *what a wine or beer state
looks like* and *which Processes apply*; the Processes themselves are ordinary
:class:`~fermentation.core.process.Process` subclasses elsewhere in the core, and
the industry-unit conversion boundary lives in ``fermentation.scenario.compile``.
Keeping the layout in the core gives the Processes (which reference variable names
like ``"S"`` and ``"N"``) and the scenario→core compile seam a single source of
truth to agree on.

The shared variables (decisions D-B / D-4):

    X      viable biomass        g/L (dry cell weight)
    S      sugar                 g/L — a *vector*: 1 slot for wine, 3 for beer
    E      ethanol               g/L
    N      yeast-assimilable N   g/L
    T      temperature           K
    CO2    evolved CO2           g/L
    X_dead ethanol-inactivated   g/L (non-viable biomass; carbon/nitrogen still
                                 counted, but no longer catalytic — decision D-13)
    Gly    glycerol              g/L (realised-yield byproduct sink — decision D-16)
    Byp    minor byproducts      g/L (lumped organic acids / higher alcohols,
                                 carbon-accounted as succinic acid — decision D-16)
    esters esters                g/L (aroma byproducts; lumped produced-only pool)
    fusels fusel/higher alcohols g/L (Ehrlich pathway; lumped produced-only pool)
    acetolactate α-acetolactate  g/L (vicinal-diketone precursor reservoir — decision D-26)
    diacetyl diacetyl (VDK)      g/L (buttery off-note; produced then reabsorbed — D-26)
    butanediol 2,3-butanediol    g/L (flavour-inactive diacetyl-reduction product — D-26)
    acetaldehyde acetaldehyde    g/L (main-pathway intermediate; transient ethanol-carbon
                                 buffer, produced then reduced back to ethanol — D-27)
    h2s      hydrogen sulfide     g/L (sulfidic "rotten egg" off-aroma; the *residual*
                                 dissolved pool, de-repressed at low nitrogen; carbon-free
                                 — decisions D-29 production / D-42 CO2-stripping sink)
    h2s_gas  H2S swept to gas     g/L (headspace bookkeeping; h2s + h2s_gas = cumulative
                                 produced; carbon-free, on no ledger — decision D-42)
    citrate  citric acid          g/L (wine-only must input; O. oeni co-metabolises it into
                                 MLF-derived diacetyl; carbon-active, not charge-active — D-31)

Sugar is always a vector so beer's sequential glucose → maltose → maltotriose
uptake needs no structural change to also support wine's single lumped sugar.
``X_dead``, ``Gly``, ``Byp``, ``esters``, ``fusels`` and the VDK pools
(``acetolactate``/``diacetyl``/``butanediol``) start at zero at pitch and are only
accumulated by the kinetics, so they declare a default initial of 0
(`VarSpec.default`) and need not be named at every initial-condition call site. The
``esters``/``fusels`` pools are filled by the Tier-2 byproduct Processes wired below;
the three VDK pools by the diacetyl-pathway Processes (decision D-26).
Under **decision D-19 (option a1)** those Processes route the aroma carbon *out of
``S``* and ``total_carbon`` weights the pools (as ethyl acetate / isoamyl alcohol), so
``esters``/``fusels`` are real carbon-accounted state alongside ``Gly``/``Byp`` — not
diagnostic re-expressions. The former ``Byp`` double-count (it once lumped higher
alcohols) is resolved by carving them out of ``Y_byproduct_sugar``; the draw touches
only ``S`` (never ``E``/``CO2``), so turning the byproducts on perturbs the core only
by the trace sugar they consume. See ``docs/plans/milestone-2-tasks.md`` and the
``kinetics.byproducts`` module docstring.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from fermentation.core.kinetics import (
    AcetaldehydeProduction,
    AcetaldehydeReduction,
    AcetolactateDecarboxylation,
    AcetolactateExcretion,
    AlphaKetoglutarateExcretion,
    AlphaKetoglutarateReassimilation,
    AminoAcidAssimilation,
    ArrheniusTemperature,
    AutolyticHydrogenSulfide,
    AutolyticMercaptan,
    BiomassCarryingCapacity,
    BrettDeath,
    BrettDecarboxylation,
    BrettEthanolToxicity,
    BrettGrowth,
    BrettVinylphenolReduction,
    ColemanQuadraticDeathTemperature,
    DiacetylReduction,
    EsterHydrolysis,
    EsterSynthesis,
    EsterVolatilization,
    EthanolInactivation,
    FuselAlcoholsEhrlich,
    FuselAminoAcidReroute,
    GrowthNitrogenLimited,
    HydrogenSulfideProduction,
    HydrogenSulfideVolatilization,
    IsoAlphaAcidLoss,
    MalolacticCitrateMetabolism,
    MalolacticConversion,
    MalolacticDeath,
    MalolacticGrowth,
    MalolacticSenescence,
    OenococcusDiacetylReduction,
    PyruvateExcretion,
    PyruvateReassimilation,
    SugarUptakeToEthanolCO2,
    TemperatureRamp,
    YeastAutolysis,
    YeastPOFDecarboxylation,
)
from fermentation.core.process import Process, ProcessSet, RateModifier
from fermentation.core.state import StateSchema, VarSpec


def _common_specs(sugar: VarSpec) -> list[VarSpec]:
    """The six state variables every medium shares, around its ``sugar`` spec.

    Order fixes the flat-array layout, so it is part of the contract: biomass
    first, then the (scalar or vector) sugar, then ethanol / nitrogen /
    temperature / evolved CO2.
    """
    return [
        VarSpec("X", "g/L", description="viable biomass (dry cell weight)"),
        sugar,
        VarSpec("E", "g/L", description="ethanol"),
        VarSpec("N", "g/L", description="yeast-assimilable nitrogen"),
        VarSpec("T", "K", description="temperature"),
        VarSpec("CO2", "g/L", description="evolved CO2"),
        VarSpec(
            "X_dead", "g/L", default=0.0, description="ethanol-inactivated (non-viable) biomass"
        ),
        VarSpec("Gly", "g/L", default=0.0, description="glycerol (realised-yield byproduct)"),
        VarSpec(
            "Byp",
            "g/L",
            default=0.0,
            description="minor byproducts (organic acids/higher alcohols; succinic-equivalent)",
        ),
        VarSpec(
            "esters",
            "g/L",
            default=0.0,
            description="esters (fermentation aroma; lumped produced-only pool)",
        ),
        VarSpec(
            "fusels",
            "g/L",
            default=0.0,
            description="fusel / higher alcohols (Ehrlich pathway; lumped produced-only pool)",
        ),
        VarSpec(
            "esters_gas",
            "g/L",
            default=0.0,
            description="esters lost to the headspace by CO2 stripping (volatilized; "
            "carbon-bookkeeping pool, decision D-20)",
        ),
        VarSpec(
            "acetolactate",
            "g/L",
            default=0.0,
            description="alpha-acetolactate — vicinal-diketone precursor reservoir "
            "(spontaneously decarboxylates to diacetyl; decision D-26)",
        ),
        VarSpec(
            "diacetyl",
            "g/L",
            default=0.0,
            description="diacetyl (2,3-butanedione) — buttery vicinal diketone; "
            "produced then yeast-reabsorbed (the diacetyl rest, decision D-26)",
        ),
        VarSpec(
            "butanediol",
            "g/L",
            default=0.0,
            description="2,3-butanediol — flavour-inactive terminal product of "
            "diacetyl reduction by viable yeast (decision D-26)",
        ),
        VarSpec(
            "acetaldehyde",
            "g/L",
            default=0.0,
            description="acetaldehyde (ethanal) — main-pathway intermediate; a transient "
            "ethanol-carbon buffer (produced then yeast-reduced back to ethanol; D-27)",
        ),
        VarSpec(
            "h2s",
            "g/L",
            default=0.0,
            description="hydrogen sulfide (H2S) — 'rotten egg' sulfidic off-aroma; the *residual* "
            "(dissolved) pool, de-repressed at low yeast-assimilable nitrogen; carbon-free "
            "(D-29 production; D-42 CO2-stripping sink makes this residual, not cumulative)",
        ),
        VarSpec(
            "h2s_gas",
            "g/L",
            default=0.0,
            description="hydrogen sulfide swept out of the liquid by the CO2 stream "
            "(headspace bookkeeping pool; carbon-free, on no ledger; h2s + h2s_gas is "
            "cumulative H2S produced — decision D-42)",
        ),
    ]


def wine_schema() -> StateSchema:
    """Wine state layout: a single lumped fermentable sugar slot, plus the wine-only
    charge-active acid + strong-cation slots the pH charge-balance solver reads
    (decision D-18), the free-SO₂ pool the molecular-SO₂ readout reads (decision D-22),
    the ``X_mlf`` malolactic-catalyst slot (decision D-23), the ``citrate`` slot
    *O. oeni* co-metabolises into MLF-derived diacetyl (decision D-31), the dosed
    ``amino_acids`` pool the amino-acid ledger swap funds biomass from (decision D-32), and the
    ``debris`` pool yeast autolysis routes non-assimilable cell-wall carbon into (decision D-34),
    and the ``X_mlf_dead`` pool bacterial death settles killed *O. oeni* biomass into (D-39).

    These ten slots are appended to ``wine_schema`` only (not ``_common_specs``), so
    ``beer_schema`` is untouched — beer's pH is a phosphate-buffered different acid
    system with no sourced data yet, explicitly deferred. ``default=0.0`` is
    load-bearing: existing wine scenarios/tests that name no acids still compile (all
    ten → 0), and with acids, cation, SO₂, ``X_mlf``, ``X_mlf_dead``, ``citrate``,
    ``amino_acids`` and ``debris`` at 0 the slots are inert — they
    contribute 0 to every conservation sum, so the validated core and its tests are
    untouched (prime directive #3). The acid/cation/SO₂ slots have no Process touching
    them in D-18/D-22; under D-23 :class:`~fermentation.core.kinetics.malolactic.\
    MalolacticConversion` depletes ``malic`` / grows ``lactic`` / evolves ``CO2`` *only
    when ``X_mlf`` is dosed* (and is disabled at the compile seam otherwise), so undosed
    wine runs keep a constant acid trajectory. Once *O. oeni* is pitched ``X_mlf`` is real
    biomass (weighted in ``total_carbon``/``total_nitrogen`` at the biomass fractions since the
    MLF-growth beat, D-38): :class:`~fermentation.core.kinetics.malolactic.MalolacticGrowth`
    grows it and :class:`~fermentation.core.kinetics.malolactic.MalolacticDeath` kills it into
    the ``X_mlf_dead`` lees (a carbon/nitrogen-neutral transfer, both pools weighted at the same
    fractions — decision D-39). On an un-pitched run both slots stay 0 (constant ⇒ 0 drift). pH is
    simply not meaningful for a no-acid scenario and is only *computed* when requested
    (``fermentation.analysis``). ``cation_charge`` is a charge density (mol⁺/L), not a
    mass concentration — state is already heterogeneous (``T`` in K) — back-solved from
    the scenario's measured ``initial_ph`` at compile and held constant (D-18).
    ``so2_total`` (g/L of SO₂-equivalent) is a dosed input read by ``acidbase.speciate_so2``
    to derive the free/bound split (acetaldehyde-bound vs free) and the antimicrobial
    molecular fraction at the solved pH; it is **not** in the charge balance (readout-only,
    D-22/D-28) and is carbon-free, so it leaves both pH and ``total_carbon`` unchanged.
    """
    specs = _common_specs(VarSpec("S", "g/L", description="fermentable sugar"))
    specs += [
        VarSpec("tartaric", "g/L", default=0.0, description="tartaric acid (must input; diprotic)"),
        VarSpec(
            "malic",
            "g/L",
            default=0.0,
            description="L-malic acid (must input; diprotic; MLF substrate)",
        ),
        VarSpec(
            "lactic", "g/L", default=0.0, description="L-lactic acid (produced-only; MLF product)"
        ),
        VarSpec(
            "citrate",
            "g/L",
            default=0.0,
            description="citric acid (must input; O. oeni co-metabolises it during MLF — the "
            "carbon source for MLF-derived diacetyl; carbon-active but not charge-active, D-31)",
        ),
        VarSpec(
            "cation_charge",
            "mol/L",
            default=0.0,
            description="net strong-cation charge (K+-dominant), constant; "
            "back-solved from initial_ph (D-18)",
        ),
        VarSpec(
            "so2_total",
            "g/L",
            default=0.0,
            description="total SO2 (as SO2); dosed input, inert/conserved; free/bound "
            "split + molecular-fraction readout derived at solved pH (D-22/D-28)",
        ),
        VarSpec(
            "X_mlf",
            "g/L",
            default=0.0,
            description="Oenococcus oeni viable biomass — the malolactic catalyst (scales the "
            "malolactic rate). Dosed at pitch; grown from amino acids (MalolacticGrowth, D-38); "
            "killed off into X_mlf_dead by molecular SO₂ (MalolacticDeath, D-39) and by benign "
            "baseline senescence (MalolacticSenescence, D-41)",
        ),
        VarSpec(
            "X_mlf_dead",
            "g/L",
            default=0.0,
            description="non-viable Oenococcus oeni biomass — the settled bacterial lees the SO₂ "
            "kill (MalolacticDeath, D-39) and benign senescence (MalolacticSenescence, D-41) move "
            "X_mlf into (carbon/nitrogen still counted at the biomass fractions, but no longer "
            "catalytic; racked off with the other lees)",
        ),
        VarSpec(
            "amino_acids",
            "g/L",
            default=0.0,
            description="assimilable amino-acid pool (dosed must input AND autolysis-refilled; "
            "represented as arginine). Carbon- AND nitrogen-bearing: the AminoAcidAssimilation "
            "swap funds a fraction of biomass from it, refunding sugar + ammonium N (D-32); "
            "YeastAutolysis refills it from dead biomass post-AF (decision D-34)",
        ),
        VarSpec(
            "debris",
            "g/L",
            default=0.0,
            description="non-assimilable cell-wall debris (glucan/mannoprotein; produced-only). "
            "The carbon-rich remainder yeast autolysis leaves behind after releasing the "
            "nitrogen-rich amino acids — carbon-accounted as glucan, nitrogen-free (D-34)",
        ),
        VarSpec(
            "hydroxycinnamics",
            "g/L",
            default=0.0,
            description="p-coumaric-acid must precursor (the p-coumaric branch of the Brett/POF "
            "volatile-phenol pathway; decision D-40). Decarboxylated to vinylphenols by "
            "Brettanomyces (and POF+ yeast). Split from ferulic_acid at decision D-55 — the two "
            "precursors are genuinely distinct molecules (9 C vs 10 C), not a fixed-ratio lump",
        ),
        VarSpec(
            "vinylphenols",
            "g/L",
            default=0.0,
            description="4-vinylphenol — the p-coumaric-branch decarboxylase→reductase "
            "intermediate reservoir (produced-only). POF+ yeast fills it but cannot clear it; "
            "Brettanomyces reduces it to ethylphenols (decision D-40; split from vinylguaiacols "
            "at D-55)",
        ),
        VarSpec(
            "ethylphenols",
            "g/L",
            default=0.0,
            description="4-ethylphenol — the p-coumaric-branch terminal Brett volatile-phenol "
            "off-aroma ('horse-sweat/barnyard'; produced-only readout, decision D-40; split from "
            "ethylguaiacols at D-55)",
        ),
        VarSpec(
            "ferulic_acid",
            "g/L",
            default=0.0,
            description="ferulic-acid must precursor — the second Brett/POF volatile-phenol "
            "branch, split out from the p-coumaric-only hydroxycinnamics pool because ferulic is "
            "a genuinely distinct molecule (10 C, vs p-coumaric's 9 C), not a fixed-ratio split "
            "(decision D-55). Decarboxylated to vinylguaiacols by Brettanomyces (and POF+ yeast), "
            "the same enzyme (Pad1/Fdc1) as the p-coumaric branch, at a literature-sourced "
            "slower relative rate (Edlin et al. 1998)",
        ),
        VarSpec(
            "vinylguaiacols",
            "g/L",
            default=0.0,
            description="4-vinylguaiacol — the ferulic-branch decarboxylase→reductase "
            "intermediate reservoir (produced-only; decision D-55), the counterpart to "
            "vinylphenols. POF+ yeast fills it but cannot clear it; Brettanomyces reduces it to "
            "ethylguaiacols (Tchobanov et al. 2008 confirm the same reductase acts on both "
            "vinylguaiacol and vinylphenol)",
        ),
        VarSpec(
            "ethylguaiacols",
            "g/L",
            default=0.0,
            description="4-ethylguaiacol — the ferulic-branch terminal Brett volatile-phenol "
            "off-aroma ('clove/smoky'; produced-only readout, decision D-55), the counterpart to "
            "ethylphenols",
        ),
        VarSpec(
            "X_brett",
            "g/L",
            default=0.0,
            description="Brettanomyces bruxellensis viable biomass — the spoilage catalyst scaling "
            "the decarboxylase/reductase rates. Dosed at pitch; grown (BrettGrowth, D-40 pt2) and "
            "killed off into X_brett_dead by SO₂ (BrettDeath, D-40 pt3) or high ethanol "
            "(BrettEthanolToxicity, D-58)",
        ),
        VarSpec(
            "X_brett_dead",
            "g/L",
            default=0.0,
            description="non-viable Brettanomyces biomass — the settled lees BrettDeath/"
            "BrettEthanolToxicity move X_brett into (carbon/nitrogen still counted at the biomass "
            "fractions, no longer catalytic; racked off with the other lees, decisions D-40/D-58)",
        ),
        VarSpec(
            "mercaptans",
            "g/L",
            default=0.0,
            description="lumped volatile thiols (methanethiol stand-in) — the carbon-bearing "
            "reductive off-aroma. AutolyticMercaptan fills it as a yield on the autolysis flux, "
            "drawing carbon from amino_acids and deaminating the nitrogen to N (Option A, D-45); "
            "carbon-accounted as methanethiol, nitrogen-free. Copper-fined out by add_copper",
        ),
        VarSpec(
            "pyruvate",
            "g/L",
            default=0.0,
            description="excreted overflow pyruvate (C3 keto-acid; excreted-then-reassimilated). "
            "PyruvateExcretion draws it from sugar during active ferment; the flux-linked "
            "(co-metabolic) PyruvateReassimilation returns it to ethanol+CO2 and stops at dryness, "
            "freezing a persistent finished-wine residual — the second-strongest SO2-binding "
            "carbonyl after acetaldehyde (D-49)",
        ),
        VarSpec(
            "alpha_ketoglutarate",
            "g/L",
            default=0.0,
            description="excreted overflow alpha-ketoglutarate (C5 keto-acid; excreted-then-"
            "reassimilated, same structure as pyruvate). AlphaKetoglutarateExcretion draws it "
            "from sugar during active ferment; the flux-linked (co-metabolic) "
            "AlphaKetoglutarateReassimilation returns it to ethanol+CO2 at the Gay-Lussac 2:1 "
            "carbon split and stops at dryness, freezing a persistent finished-wine residual — "
            "the third SO2-binding carbonyl, after acetaldehyde and pyruvate (D-50)",
        ),
    ]
    return StateSchema(specs)


def beer_schema() -> StateSchema:
    """Beer state layout: three sugars consumed sequentially.

    Glucose is taken up first, then maltose, then maltotriose — the order the
    ``components`` tuple records and the sugar-uptake Process will honour.
    """
    specs = _common_specs(
        VarSpec(
            "S",
            "g/L",
            size=3,
            description="fermentable sugars (sequential uptake)",
            components=("glucose", "maltose", "maltotriose"),
        )
    )
    # Iso-alpha-acids (isohumulones) — the bitter compounds (decision D-64). Made in the boil
    # by thermal isomerization of hop alpha-acids (computed at the compile seam and wired at
    # t=0, like initial_ph), then lost during fermentation by yeast adsorption (IsoAlphaAcidLoss).
    # BEER-ONLY (appended here, not in _common_specs, so wine_schema is untouched). Off the carbon
    # ledger (exogenous hop-derived mass, like dosed SO2), so it never perturbs total_carbon.
    # Default 0 ⇒ an unhopped beer carries no bitterness and the loss Process is inert/disabled.
    # 1 IBU ≈ 1 mg/L iso-alpha, so the ibu_series readout is this slot × 1000.
    specs.append(
        VarSpec("iso_alpha", "g/L", default=0.0, description="iso-alpha-acids (bitterness)")
    )
    return StateSchema(specs)


@dataclass(frozen=True)
class Medium:
    """A named beverage family: its state schema plus the kinetics that act on it.

    ``process_factories`` are zero-argument callables that each build one additive
    :class:`Process`; ``modifier_factories`` likewise build the multiplicative
    :class:`RateModifier` objects (ethanol inhibition, Arrhenius temperature
    dependence) that scale those Processes. Both are *factories* rather than shared
    instances so every ``build_process_set`` call gets fresh objects — two media (or
    two runs) never share a mutable Process/modifier. Kinetics read their parameters
    at ``derivatives``/``factor`` time, not construction time, so the factories need
    no arguments.

    An empty pair of tuples integrates to a constant trajectory — the honest
    "no kinetics" baseline a bare :class:`Medium` still provides.
    """

    name: str
    schema: StateSchema
    process_factories: tuple[Callable[[], Process], ...] = ()
    modifier_factories: tuple[Callable[[], RateModifier], ...] = field(default=())

    def build_process_set(self, *, strict: bool = False) -> ProcessSet:
        """Assemble this medium's Processes and modifiers into a :class:`ProcessSet`."""
        return ProcessSet(
            self.schema,
            [factory() for factory in self.process_factories],
            modifiers=[factory() for factory in self.modifier_factories],
            strict=strict,
        )


#: The validated-core primary-fermentation kinetics, as zero-argument factories.
#: Wine and beer share the *same* mechanism set — biomass growth, fermentative
#: sugar uptake, and ethanol-driven cell inactivation (the cumulative viability
#: brake that sets the fermentation timescale, Coleman 2007), with per-rate
#: temperature dependence scaling all three — Arrhenius for growth/uptake,
#: Coleman's own quadratic regression for death (``ColemanQuadraticDeathTemperature``,
#: decision D-57 — a single Arrhenius E_a cannot reproduce that curvature). The
#: only structural difference between the two media is the sugar vector (1 slot
#: vs 3): beer's
#: sequential glucose→maltose→maltotriose uptake is handled *inside*
#: :class:`~fermentation.core.kinetics.uptake.SugarUptakeToEthanolCO2` via catabolite
#: repression, so it needs no extra Process here.
#:
#: The instantaneous Luong ethanol wall (``EthanolInhibition``) is **not** wired in:
#: the cumulative inactivation Process is the mechanistically-correct ethanol brake,
#: and stacking an instantaneous wall on top would double-count ethanol toxicity
#: (decision D-13). The class is retained for optional/strain use.
_PRIMARY_FERMENTATION_PROCESSES: tuple[Callable[[], Process], ...] = (
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
    EthanolInactivation,
)
_PRIMARY_FERMENTATION_MODIFIERS: tuple[Callable[[], RateModifier], ...] = (
    ArrheniusTemperature.for_growth,
    ArrheniusTemperature.for_uptake,
    ColemanQuadraticDeathTemperature,
)

#: Tier-2 temperature-/metabolism-driven aroma byproducts (Milestone 2, decision
#: D-18/D-19): ester synthesis and Ehrlich-pathway fusel alcohols. Kept as a
#: *separate* tuple from the validated-core primary set so the speculative beat stays
#: **isolable** (prime directive #3): building a ProcessSet without this tuple is the
#: pure validated core. Under D-19 (option a1) they route aroma carbon out of ``S``
#: and ``total_carbon`` weights the ``esters``/``fusels`` pools, so they no longer
#: leave the core byte-for-byte when enabled — turning them on draws a *trace* of
#: sugar (~0.2 % of ``S0``), perturbing only ``dS`` (never ``dE``/``dCO2``). Carbon
#: still closes to machine precision with them on, and the §2.2 trio stays in band.
#: See D-19 / milestone-2-tasks.md.
#:
#: :class:`EsterVolatilization` (decision D-20) is the gas-stripping sink that moves
#: liquid ``esters`` into the bookkeeping ``esters_gas`` headspace pool as CO2 sparges
#: the must — the physics behind wine's "warmer ⇒ *less* liquid ester" (Rollero 2014):
#: with ``E_a_ester_volatil`` set *per medium* it is held **above** ``E_a_esters`` for
#: wine (stripping outruns synthesis, liquid esters fall with T) and **below** it for
#: beer (synthesis dominates, esters rise with T — de Andrés-Toro). The transfer is
#: carbon-neutral (``esters`` → ``esters_gas``, both booked as ethyl acetate), so it is
#: in this isolable tuple too and ``total_carbon`` still closes to machine precision.
_BYPRODUCT_PROCESSES: tuple[Callable[[], Process], ...] = (
    EsterSynthesis,
    FuselAlcoholsEhrlich,
    EsterVolatilization,
)

#: Vicinal-diketone (VDK / diacetyl) pathway (Milestone 2, decision D-26): the three-step
#: sugar → α-acetolactate → diacetyl + CO2 → 2,3-butanediol chain that makes the "diacetyl
#: rest" emerge. Kept as its own isolable tuple (prime directive #3): a ProcessSet built
#: without it is the prior core. Diacetyl is *intrinsic yeast metabolism* (not a dosed
#: organism like MLF), so — unlike ``_MLF_PROCESSES`` — it is wired into BOTH media and runs
#: on every default fermentation, like the ester/fusel byproducts. Turning it on draws only
#: a *trace* of sugar into the reservoir (α-acetolactate peaks ~mg/L, roughly an order of
#: magnitude below the ester draw), so it leaves ``dX``/``dE``/``dCO2``/``dN`` byte-for-byte
#: until the decarb/reduction move that carbon on; ``total_carbon`` closes to machine
#: precision throughout (each step is on the weighted ledger). One honest tier consequence
#: (D-26, the D-19 ``S`` parallel): the always-on speculative decarboxylation touches the
#: shared ``CO2`` slot, so the *structural* ``tier_of("CO2")`` drops PLAUSIBLE→SPECULATIVE —
#: but the param-aware tier users see was *already* SPECULATIVE (the uptake Process reads
#: speculative params), so there is no headline change, and the CO2 pool genuinely does hold
#: a speculative decarb trace. Excretion is temperature-flat; the temperature-criticality of
#: the rest lives in the spontaneous, non-yeast-gated decarboxylation (``E_a_decarb`` >
#: ``E_a_reduction``); reduction is gated on VIABLE ``X`` with no flux term, so a warm rest
#: with live yeast clears diacetyl fast while an early crash strands a rising diacetyl.
#: SCOPE (v1): yeast valine-pathway diacetyl only — MLF/citrate diacetyl is deferred, so wine
#: yeast-pathway diacetyl understates real wine diacetyl. VDK params live in the shared
#: ``vicinal_diketones.yaml`` (the load-bearing decarb step is non-enzymatic, medium-agnostic).
_VDK_PROCESSES: tuple[Callable[[], Process], ...] = (
    AcetolactateExcretion,
    AcetolactateDecarboxylation,
    DiacetylReduction,
)

#: Acetaldehyde pathway (Milestone 2, decision D-27): the main-pathway intermediate as a
#: transient ethanol-carbon *buffer* — flux-linked production that borrows carbon from ``E``
#: and viable-``X``-gated reduction that returns it. Kept as its own isolable tuple (prime
#: directive #3): a ProcessSet built without it is the prior core. Like the ester/VDK pools
#: (and unlike the *dosed* MLF organism), acetaldehyde is intrinsic yeast metabolism, so it
#: is wired into BOTH media and runs on every default ferment. It touches only
#: ``acetaldehyde`` and ``E`` (never ``S``/``CO2``/``N``/``X``) at the derivative level, so
#: ``dS``/``dCO2``/``dN`` are byte-for-byte given the same state; the only integrated coupling
#: is second-order (``E`` feeds the inactivation viability brake, a ~1e-4 relative path
#: perturbation). The ``E`` endpoint reconverges to the buffer-off core to relative ~1e-8 (the
#: pool fully reduces back), so the §2.2 ABV / realised-yield / CO2 benchmarks are preserved to
#: far below any tolerance. This is the owner's buffer choice (D-27) over a draw-from-sugar
#: stand-in, which would double-count the uptake
#: Process's already-complete sugar→ethanol conversion and inflate ABV with net-new ethanol.
#: One honest tier consequence (the D-26 ``CO2`` parallel): the always-on speculative
#: production is the first such Process to *write* ``E``, so the *structural* ``tier_of("E")``
#: drops PLAUSIBLE→SPECULATIVE — but the param-aware tier users see was already SPECULATIVE
#: (the uptake Process reads speculative params), so there is no headline change. Params live
#: in the shared, medium-agnostic ``acetaldehyde.yaml`` (main-pathway yeast metabolism, not a
#: beverage property). SCOPE (v1): the SO₂-binding free/bound split is a separate readout beat.
_ACETALDEHYDE_PROCESSES: tuple[Callable[[], Process], ...] = (
    AcetaldehydeProduction,
    AcetaldehydeReduction,
)

#: Hop bittering (BEER-ONLY, decision D-64): the §3.3 additive beat. The boil isomerization of
#: alpha-acids to iso-alpha-acids is a wort-side compile-seam calc (``iso_alpha_from_boil``,
#: wired into ``iso_alpha`` at t=0 like ``initial_ph``), NOT a Process — the only *dynamic*
#: content is :class:`IsoAlphaAcidLoss`, the fermentation-time adsorption of iso-alpha onto
#: viable yeast (the ~5-20% wort-to-beer bitterness drop). Kept in its own isolable tuple (prime
#: directive #3) and wired into the BEER medium only (wine has no ``iso_alpha`` slot). It touches
#: ``iso_alpha`` alone — OFF the carbon ledger (exogenous hop-derived mass, like dosed SO2), so
#: the whole beat leaves ``total_carbon`` byte-for-byte unchanged. On an unhopped beer
#: ``iso_alpha`` starts 0 and the term is inert; the compile seam additionally DISABLES it when
#: no hops are scheduled, so the empty ``iso_alpha`` slot keeps its VALIDATED tier and no flux is
#: paid (the MLF/Brett isolability pattern). Params live in the shared ``hops.yaml``.
_HOPS_PROCESSES: tuple[Callable[[], Process], ...] = (IsoAlphaAcidLoss,)

#: Aging chemistry — the slow post-fermentation "years" axis (Milestone 3 / Tier-3, decisions
#: D-68/D-69/D-70): :class:`EsterHydrolysis`, the first §4.1 aging Process (young fruity acetate
#: esters hydrolyse back toward equilibrium with age, releasing carbon 5:2 into ``fusels`` +
#: ``Byp``). MEDIUM-AGNOSTIC — acid-catalysed ester hydrolysis is a property of the molecule and
#: the wine/beer pH, not the biology (the ``vicinal_diketones.yaml`` / shared-file pattern), and
#: ``esters``/``fusels``/``Byp`` exist in both schemas — so it is wired into BOTH media. Kept in
#: its OWN isolable tuple (prime directive #3): a ProcessSet built without it is the pre-aging
#: model. Unlike the always-on intrinsic aroma pools, aging is INHERENTLY post-ferment (there is
#: no aging at t0), so the compile seam DISABLES it unconditionally and a ``begin_aging``
#: intervention (decision D-70, the ``pitch_mlf`` reconfigure pattern MINUS the state mutation)
#: re-enables it over a post-fermentation aging segment — off during the ferment, on during aging.
#: An un-aged run is thus byte-for-byte the pre-aging core (disabled ⇒ skipped by ``active`` /
#: ``tier_of`` / the strict ``touches`` check). During a post-dryness aging segment every OTHER
#: producer of ``esters``/``fusels``/``Byp`` (``ester_synthesis``, ``ester_volatilization``,
#: ``fusel_alcohols_ehrlich``, and the ``Byp`` uptake routing) is fermentative-flux-gated and
#: quiescent at ``S ≈ 0`` (``fermentative_flux_shape`` is 0 when sugar OR biomass is 0), so the
#: aging ester/fusel signal is UNCONFOUNDED — only :class:`EsterHydrolysis` moves those pools
#: (Stance A, D-70: leave the ferment set on, the aging effect emerges). Params live in the
#: shared, medium-agnostic ``aging.yaml``.
_AGING_PROCESSES: tuple[Callable[[], Process], ...] = (EsterHydrolysis,)

#: Excreted keto-acid overflow pool (wine-only, decision D-49): pyruvate as the
#: second-strongest SO₂-binding carbonyl after acetaldehyde. :class:`PyruvateExcretion`
#: draws carbon *out of ``S``* into the ``pyruvate`` pool on the fermentative flux (so it
#: fills during active ferment and stops at dryness); :class:`PyruvateReassimilation` returns
#: it to ``E``/``CO2``, *also* flux-linked (co-metabolic — NOT the no-flux ADH idiom), so both
#: terms die at dryness and the pool **freezes** at the quasi-steady plateau
#: ``k_pyruvate_excretion / k_pyruvate_reassimilation`` as a **persistent finished-wine
#: residual** — crash- and duration-independent (a no-flux viable-``X`` gate would instead
#: drain it to ~0 over the long tail, since a clean ferment ends with the yeast still viable).
#: That stranded residual is the carbonyl that will share dosed SO₂ with acetaldehyde in the
#: D-51 multi-carbonyl binding equilibrium (decision D-49, option A). Modelled as
#: an *excreted side pool*, NOT acetaldehyde's on-pathway precursor (the intracellular flux
#: pyruvate never persists and never binds SO₂ — see the ``keto_acids`` module docstring for
#: why the "route acetaldehyde through pyruvate" rework was rejected as unphysical), so
#: acetaldehyde / D-27 / D-47 / D-48 stay untouched. Kept in its own isolable tuple (prime
#: directive #3): a ProcessSet built without it is the prior core. Unlike the byte-for-byte-
#: isolable acetaldehyde buffer, excretion touches ``S`` and re-assimilation touches
#: ``E``/``CO2``, so turning it on routes a *trace* slice of sugar carbon on a detour to
#: ethanol; the only endpoint difference from the pool-off core is the stranded residual
#: (a few tens of mg/L of sugar carbon parked as pyruvate rather than fermented on), ≪ 0.1 %
#: of ABV, so the §2.2 CO2/ABV/realised-yield benchmarks are preserved far below tolerance.
#: WINE-ONLY (v1): the SO₂-binding competition it exists for is a wine readout and no §2.2
#: beer benchmark asserts a keto-acid level — beer overflow pyruvate/α-KG is deferred. Params
#: live in the shared, medium-agnostic ``keto_acids.yaml`` (overflow-keto-acid metabolism is
#: generic yeast, not a beverage property).
#:
#: :class:`AlphaKetoglutarateExcretion` / :class:`AlphaKetoglutarateReassimilation` (decision
#: D-50) add the third SO₂-binding carbonyl with the SAME structure: excretion draws the C5
#: pool from ``S``, flux-linked co-metabolic reassimilation returns it and freezes a lower
#: (~20 mg/L nominal) persistent residual at dryness. The one load-bearing difference from
#: pyruvate: the reassimilation carbon split. Pyruvate's C3 → C2(ethanol) + C1(CO2) is
#: mole-for-mole *because* 3 carbons is exactly one Gay-Lussac fermentation unit (2 carbon to
#: ethanol : 1 carbon to CO2) — the coincidence that keeps its detour stoichiometrically
#: identical to the main pathway. α-KG's C5 does not divide evenly 1:1, so its reassimilation
#: returns carbon at the SAME 2:1 ratio instead (5/3 mol ethanol + 5/3 mol CO2 per mole), not
#: mole-for-mole — copying pyruvate's form naively would have diverted reassimilation
#: *throughput* (not just the residual, ~10–20× larger) away from ethanol, large enough to
#: threaten the §2.2 ABV/CO₂ benchmarks. See the ``keto_acids`` module docstring.
_KETO_ACID_PROCESSES: tuple[Callable[[], Process], ...] = (
    PyruvateExcretion,
    PyruvateReassimilation,
    AlphaKetoglutarateExcretion,
    AlphaKetoglutarateReassimilation,
)

#: Hydrogen-sulfide production + CO₂-stripping (Milestone 2, decisions D-29 / D-42): the
#: low-nitrogen "rotten egg" off-aroma. :class:`HydrogenSulfideProduction` (D-29) is one
#: flux-linked producer gated by an *inverse*-nitrogen term; :class:`HydrogenSulfideVolatilization`
#: (D-42) is the CO₂-stripping sink that sweeps the volatile H₂S out of the liquid ``h2s`` pool
#: into the ``h2s_gas`` headspace pool as the ferment sparges CO₂ — so ``h2s`` is now the
#: *residual* (dissolved, µg/L) pool and ``h2s + h2s_gas`` is cumulative produced (the ester
#: D-19→D-20→D-21 precedent, but carbon-free, so *simpler*: neither pool is on any ledger).
#: Kept as their own isolable tuple (prime directive #3): a ProcessSet built without it is the
#: prior core, and dropping *just* the sink recovers the D-29 produced-only ``h2s`` byte-for-byte
#: (``h2s_gas`` stays 0). Like the ester/VDK/acetaldehyde pools (and unlike the *dosed* MLF
#: organism), H₂S is intrinsic yeast metabolism, so both run on every default ferment in BOTH
#: media. This is the most isolable beat in the model: H₂S is CARBON-FREE (on no conservation
#: ledger) and the Processes touch ONLY ``h2s``/``h2s_gas`` while merely *reading* ``X``/``S``/
#: ``N``/``T`` — so disabling them leaves the RHS of every other column byte-for-byte identical
#: (nothing reads ``h2s``/``h2s_gas`` to feed anything back); the integrated trajectory then
#: differs only by a ~1e-7 adaptive-solver mesh artifact, cleaner than the acetaldehyde buffer's
#: *genuine* second-order E→viability coupling (D-27). No tier headline either: they write pools
#: nothing reads, so no other column's structural tier drops (contrast the D-26 ``CO2`` / D-27
#: ``E`` cases). Params live in the shared, medium-agnostic ``hydrogen_sulfide.yaml`` (both
#: sulfate-reduction and the Henry's-law stripping are generic, medium-agnostic physics).
_H2S_PROCESSES: tuple[Callable[[], Process], ...] = (
    HydrogenSulfideProduction,
    HydrogenSulfideVolatilization,
)

#: Malolactic fermentation (wine-only, decision D-23): the *Oenococcus oeni* malate →
#: lactate + CO2 conversion, the first RHS consumer of the D-18 pH solver and the D-22
#: molecular-SO₂ readout. Kept as its own tuple so it stays **isolable** (prime directive
#: #3): the conversion contributes zero before the pH solve whenever ``X_mlf`` is undosed
#: (structural *value* isolability), and the compile seam *disables* it when MLF is not
#: pitched so the inert ``malic``/``lactic`` slots keep their VALIDATED tier rather than
#: being dragged to speculative by an enabled-but-zero Process (*tier* isolability —
#: ``ProcessSet.tier_of`` counts enabled, not nonzero, Processes). Wine-only: beer has no
#: ``malic``/``lactic`` slots, so it is never wired there.
#:
#: MLF-derived diacetyl (decision D-31) adds two more *O. oeni* Processes to this same dosed,
#: isolable tuple: :class:`MalolacticCitrateMetabolism` co-metabolises the dosed ``citrate`` must
#: input into α-acetolactate + CO2 (feeding the shared VDK reservoir, so diacetyl emerges from
#: the always-on D-26 decarboxylation), and :class:`OenococcusDiacetylReduction` clears diacetyl
#: on the lees (``X_mlf``-gated). Both are disabled at the compile seam with the malate Process
#: when *O. oeni* is un-pitched, so an un-pitched wine run stays byte-for-byte the validated core
#: and the ``citrate`` slot keeps its VALIDATED tier (like ``malic``/``lactic``). Citrate — not
#: sugar — sources this carbon because MLF-diacetyl is a late/post-dryness phenomenon and the
#: sugar-draw helper no-ops at ``S=0`` (decision D-31; see the malolactic module docstring).
#: :class:`MalolacticDeath` (decision D-39) rides in this same pitch-gated tuple: it moves viable
#: ``X_mlf`` into ``X_mlf_dead`` under **molecular SO₂** (``1 − g_SO₂``) with its own Arrhenius
#: temperature factor, so bacteria die off when SO₂ is dosed — the mechanism that lets an SO₂
#: addition (or a rack removing the bacteria) *lock in* MLF-derived diacetyl by halting
#: :class:`OenococcusDiacetylReduction`. :class:`MalolacticSenescence` (MLF v2, decision D-41) rides
#: alongside it: the *benign baseline* mortality (``k_senescence_mlf · X_mlf · arrhenius(T)``, no
#: SO₂/pH/ethanol term) that lifts the v1 "unsulfited bacteria never die" tradeoff — over
#: weeks-to-months a pitched, untreated culture slowly declines into the same ``X_mlf_dead`` pool.
#: Both are pitch-gated (not amino-acid-gated like growth): bacteria age and die whether or not they
#: were growing, so they belong with the conversion set, disabled at the compile seam on an
#: un-pitched run. Both transfers are carbon/nitrogen-neutral (both pools weighted at the biomass
#: fractions since D-38/D-39), so they add no conservation code.
_MLF_PROCESSES: tuple[Callable[[], Process], ...] = (
    MalolacticConversion,
    MalolacticCitrateMetabolism,
    OenococcusDiacetylReduction,
    MalolacticDeath,
    MalolacticSenescence,
)

#: Malolactic *growth* (wine-only, the deferred MLF-growth beat, decision D-38). Makes ``X_mlf``
#: dynamic: :class:`MalolacticGrowth` builds O. oeni biomass from the ``amino_acids`` pool (D-32,
#: autolysis-refilled D-34), which — since :class:`MalolacticConversion` is linear in ``X_mlf`` —
#: accelerates deacidification autocatalytically. Kept in its OWN tuple, DELIBERATELY SEPARATE
#: from ``_MLF_PROCESSES`` because it is gated on a different feature: amino-acid fuel, NOT the
#: pitch. The compile seam disables it when ``amino_acids_gpl ≤ 0`` (the swap/re-route gate), which
#: alone prevents the tier-isolability regression — every pitched-but-not-aa-dosed D-23/D-31 run
#: keeps it disabled, so it never drags the ``amino_acids``/``S``/``X_mlf`` tier via ``tier_of``.
#: It is NOT additionally gated on the pitch: the Process's own ``X_mlf ≤ 0`` guard keeps it inert
#: until bacteria are present, and whether post-pitch bacteria GROW is left to the emergent
#: environmental gate (the ethanol wall etc.), mirroring how conversion trusts its gate rather than
#: a compile rule — so co-inoculation dominance is emergent, not hard-coded (D-38). Wine-only.
_MLF_GROWTH_PROCESSES: tuple[Callable[[], Process], ...] = (MalolacticGrowth,)

#: *Brettanomyces* volatile-phenol spoilage (wine-only, decision D-40): the mixed-culture beat that
#: closes Milestone 2. :class:`BrettDecarboxylation` takes must ``hydroxycinnamics`` →
#: ``vinylphenols`` + CO2 and :class:`BrettVinylphenolReduction` reduces ``vinylphenols`` →
#: ``ethylphenols`` — Brett carries BOTH enzymes, so a dosed culture spoils POF-negative wine
#: unaided (the canonical funk mechanism). Kept in its own tuple so it stays **isolable** (prime
#: directive #3), mirroring the *dosed* MLF organism (and unlike the always-on intrinsic aroma
#: pools): the Processes contribute zero before any pH work when ``X_brett`` is undosed, and the
#: compile seam DISABLES them when Brett is not pitched so the inert ``hydroxycinnamics``/
#: ``vinylphenols``/``ethylphenols`` slots keep their VALIDATED tier (``tier_of`` counts enabled,
#: not nonzero, Processes — the D-23 MLF pattern). :class:`BrettGrowth` (D-40 pt2) is amino-acid-
#: gated in its own tuple below; :class:`BrettDeath` (D-40 pt3, the SO₂ lever) and
#: :class:`BrettEthanolToxicity` (D-58, the ethanol-toxicity lever — needs no SO₂) both ride in THIS
#: pitch-gated tuple — Brett dies whether or not it was growing, so they belong with the phenol
#: Processes, disabled at the compile seam on an unpitched run (mirroring how
#: :class:`~fermentation.core.kinetics.malolactic.MalolacticDeath` sits in ``_MLF_PROCESSES``, not
#: the amino-acid-gated growth tuple). The ``X_brett → X_brett_dead`` transfer is carbon/nitrogen-
#: neutral (both pools weighted at the biomass fractions since pt2), so it adds no new ledger code.
#: Wine-only: beer has no ``hydroxycinnamics``/phenol slots, so Brett is never wired there.
_BRETT_PROCESSES: tuple[Callable[[], Process], ...] = (
    BrettDecarboxylation,
    BrettVinylphenolReduction,
    BrettDeath,
    BrettEthanolToxicity,
)

#: *Brettanomyces* growth (wine-only, decision D-40 pt2). Makes ``X_brett`` dynamic:
#: :class:`BrettGrowth` builds Brett biomass from the ``amino_acids`` pool (D-32, autolysis-refilled
#: D-34) — but draws its carbon shortfall from **ethanol**, not sugar, so Brett grows in a *dry*
#: finished wine (its post-AF/barrel niche), and the phenol spoilage then *accelerates* as the
#: population multiplies (decarboxylase/reductase are linear in ``X_brett``). Kept in its OWN tuple,
#: DELIBERATELY SEPARATE from ``_BRETT_PROCESSES`` because it is gated on a different feature:
#: amino-acid fuel, NOT the Brett pitch (the exact ``_MLF_GROWTH_PROCESSES`` split). The compile
#: seam disables it when ``amino_acids_gpl ≤ 0``; the Process's own ``X_brett ≤ 0`` guard keeps it
#: inert until Brett is present. Wine-only.
_BRETT_GROWTH_PROCESSES: tuple[Callable[[], Process], ...] = (BrettGrowth,)

#: POF+ yeast decarboxylase (wine-only, decision D-40 pt4). :class:`YeastPOFDecarboxylation` is the
#: *yeast* half of the volatile-phenol story: a POF+ (phenolic-off-flavour-positive) primary strain
#: carries the cinnamate decarboxylase - the same reaction as :class:`BrettDecarboxylation`, drawing
#: the same ``hydroxycinnamics`` pool 9 = 8 + 1 into ``vinylphenols`` + CO2 - but **not** the
#: reductase, so it fills the shared reservoir it cannot drain (``vinylphenols`` strand with no
#: Brett; a later Brett gets a head start on the pre-filled reservoir - the emergent yeast/Brett
#: coupling, the D-26/D-31 parallel). Flux-coupled to the yeast's fermentative activity (catalyst =
#: viable ``X``, not ``X_brett``) so it runs during AF and stops at dryness; temperature-flat (the
#: :class:`~fermentation.core.kinetics.vicinal_diketones.AcetolactateExcretion` precedent). Kept in
#: its OWN tuple, DELIBERATELY SEPARATE from ``_BRETT_PROCESSES``, because POF+ is a *binary strain
#: trait* gated on its own opt-in (``pof_positive``), WHOLLY INDEPENDENT of the Brett pitch: a POF+
#: ferment need not have Brett, and a POF-negative default wine must make no vinylphenol. The
#: compile seam DISABLES it unless the strain is opted in, so a default (POF-) run is byte-for-byte
#: validated core and the phenol slots keep their VALIDATED tier (``tier_of`` counts enabled, not
#: nonzero, Processes - the Brett-unpitched pattern). Wine-only (beer has no phenol slots).
_POF_PROCESSES: tuple[Callable[[], Process], ...] = (YeastPOFDecarboxylation,)

#: Biomass carrying-capacity cap (wine-only, decision D-30): the opt-in residual-nitrogen
#: floor. A logistic ``(1 - X/K)`` RateModifier on growth that saturates biomass below the
#: nitrogen ceiling, leaving a dose-dependent residual of yeast-assimilable nitrogen — which
#: restores the D-29 cross-must H₂S inverse-N lever (muted in the core because growth strips
#: YAN to ~0 at every dose). A residual-N floor is a deliberate DEPARTURE from the validated
#: Coleman 2007 anchor (which caps nothing; ``test_coleman_reconstruction`` pins the match at
#: 80 *and* 330 mg N/L), so — like the *dosed* MLF organism and unlike the always-on intrinsic
#: aroma pools — it is kept in its own tuple and the compile seam DISABLES it unless a scenario
#: opts in via ``carrying_capacity_gpl``. Disabled ⇒ factor 1 *and* excluded from tier
#: derivation, so an undosed wine run is byte-for-byte the validated core and growth stays
#: PLAUSIBLE. Wine-only in v1 (the H₂S lever and the prospective MLF-with-growth model are wine
#: concerns), mirroring the wine-only MLF wiring; beer carrying capacity is deferred.
_CARRYING_CAPACITY_MODIFIERS: tuple[Callable[[], RateModifier], ...] = (BiomassCarryingCapacity,)

#: Amino-acid ledger (wine-only, decision D-32): the toggleable ``amino_acids`` pool the
#: :class:`AminoAcidAssimilation` swap funds a fraction of biomass from — refunding sugar
#: carbon and ammonium nitrogen so the pool sits on *both* conservation ledgers. Kept in its
#: own isolable tuple (prime directive #3): like the *dosed* MLF organism (and unlike the
#: always-on intrinsic aroma pools), it contributes only when amino acids are dosed, so the
#: compile seam DISABLES it when ``amino_acids_gpl`` ≤ 0 — an undosed wine run is byte-for-byte
#: the validated core and the empty ``amino_acids`` slot keeps its VALIDATED tier. Dosed, the
#: swap correctly perturbs the run (refunded N/S act like supplementary YAN) and its speculative
#: tier drops growth's ``S``/``N`` outputs to speculative. CORRECTNESS COUPLING (decision D-32):
#: the swap's refund must track growth's *realised* (post-modifier) draw, so the wine growth
#: Arrhenius (:data:`_WINE_FERMENTATION_MODIFIERS`) and the carrying-capacity modifier both name
#: it in their ``modifies`` — otherwise a cold ferment or a near-saturation carrying cap (M < 1)
#: would let the base-rate refund exceed the scaled draw and create sugar. Wine-only; beer
#: deferred with the wine-only nitrogen model (D-30).
#:
#: :class:`FuselAminoAcidReroute` (decision D-33) rides in this same dosed, wine-only tuple: it
#: re-sources a fraction of Ehrlich fusel carbon off its sugar stand-in and onto the ``amino_acids``
#: pool, **deaminating** the consumed amino acids' nitrogen to ammonium ``N`` — the deamination
#: branch the fusel re-route was deferred on (D-19/D-32). Unlike the swap it is NOT scaled by the
#: growth Arrhenius / carrying-cap modifiers: it recomputes the *fusel* production rate (which
#: carries its own ``E_a_fusels`` Arrhenius and is scaled by no RateModifier), so to refund exactly
#: what :class:`FuselAlcoholsEhrlich` drew it must stay unscaled too — the producer and re-route
#: share :func:`~fermentation.core.kinetics.byproducts.fusel_production_rate` and neither is a
#: modifier target. Disabled with the swap at the compile seam when amino acids are un-dosed.
_AMINO_ACID_PROCESSES: tuple[Callable[[], Process], ...] = (
    AminoAcidAssimilation,
    FuselAminoAcidReroute,
)

#: Yeast autolysis (wine-only, decisions D-34, D-44): the autolytic-peptide source that refills the
#: ``amino_acids`` pool from dead biomass (``X_dead``) post-AF — the second prerequisite (after the
#: D-33 fusel re-route) the deferred MLF-with-growth beat needs, since the pool is empty at the MLF
#: pitch point (D-23). The first consumer of ``X_dead``: it liberates the dead-cell nitrogen as
#: amino acids and routes the carbon-rich remainder to the ``debris`` pool (carbon + nitrogen close
#: separately). Like the *dosed* MLF organism / carrying cap and UNLIKE the always-on intrinsic
#: aroma pools, it *consumes* core state (``X_dead``), so it is kept isolable and the compile
#: seam DISABLES it unless a scenario opts in via ``autolysis_rate_per_h`` — an undosed wine run is
#: then byte-for-byte the validated core. Wine-only (mirrors the wine-only ``amino_acids`` pool and
#: nitrogen model, D-30/D-32); beer deferred.
#:
#: :class:`AutolyticHydrogenSulfide` (decision D-44) rides in this same opt-in tuple: it feeds the
#: shared ``h2s`` pool a **yield on the autolysis flux** (``y_h2s_autolysis·k_autolysis·f_T·
#: X_dead``) — the sulfide dead cells release as they self-digest. Sharing the gate keeps peptide
#: and sulfide release on one clock (both read the ``autolysis_rate_per_h`` override), and its
#: **non-flux-linked** form is the point: the D-42 CO₂-stripping sink gates off at dryness, so this
#: autolytic H₂S accumulates un-stripped as *residual* — the sur-lie "reduction" fault. Carbon-free,
#: touches only ``h2s`` (nothing reads it back), so like the D-34 refill it stays isolable and drops
#: to the validated core when autolysis is un-opted.
#:
#: :class:`AutolyticMercaptan` (decision D-45) rides here too — the *carbon-bearing* twin: it fills
#: the ``mercaptans`` (thiol) pool on the same autolysis flux, but draws the mercaptan carbon from
#: ``amino_acids`` and **deaminates** the nitrogen to ``N`` (Option A, the D-33 idiom — methanethiol
#: carries carbon, unlike H₂S, so it cannot draw from nothing). Also non-flux-linked ⇒ accumulates
#: un-stripped post-dryness. It is the **first autolysis-gated ``N``-writer**, so an autolysis-on
#: run drops the structural ``tier_of("N")`` to speculative (the D-27 ``E`` parallel). All three
#: Processes are disabled together at the compile seam absent ``autolysis_rate_per_h``.
_AUTOLYSIS_PROCESSES: tuple[Callable[[], Process], ...] = (
    YeastAutolysis,
    AutolyticHydrogenSulfide,
    AutolyticMercaptan,
)

#: Temperature-schedule ramp (decision D-35): the single Process that drives ``T`` along a
#: piecewise-linear temperature schedule (``dT/dt = temperature_ramp_rate``). Medium-agnostic —
#: cellar temperature is not a beverage property — so wired into BOTH media, and (unlike the
#: dosed/opt-in gated tuples) kept ALWAYS ENABLED: it reads the slope with a ``0.0`` isothermal
#: default, so an un-ramped run contributes exactly ``0.0`` to ``dT/dt`` and is byte-for-byte the
#: pre-ramp core (``T`` stays VALIDATED). The per-segment slope is supplied by the runtime event
#: loop (``simulate_scheduled``); the scenario compile boundary injects the provenance-backed
#: ``temperature_ramp_rate`` parameter and the slope-change events only when the schedule ramps.
_TEMPERATURE_PROCESSES: tuple[Callable[[], Process], ...] = (TemperatureRamp,)

#: Wine growth/uptake Arrhenius modifiers (decision D-32). Identical to
#: :data:`_PRIMARY_FERMENTATION_MODIFIERS` except the growth Arrhenius *also* scales the
#: amino-acid swap (``for_growth`` extra target), so the swap's carbon/nitrogen refunds carry
#: the same temperature factor as growth's draw (see the correctness coupling above). Beer keeps
#: the plain :data:`_PRIMARY_FERMENTATION_MODIFIERS` (no amino-acid pool). ``for_uptake`` is
#: unchanged — the swap tracks *growth*, not the fermentative sugar-uptake flux.
_WINE_FERMENTATION_MODIFIERS: tuple[Callable[[], RateModifier], ...] = (
    lambda: ArrheniusTemperature.for_growth(AminoAcidAssimilation.name),
    ArrheniusTemperature.for_uptake,
    ColemanQuadraticDeathTemperature,
)


#: The registry of known media. Adding a beverage family = adding an entry here
#: (and, at the I/O boundary, an initial-composition vocabulary in
#: ``fermentation.scenario.compile``).
MEDIA: dict[str, Medium] = {
    "wine": Medium(
        name="wine",
        schema=wine_schema(),
        process_factories=(
            _PRIMARY_FERMENTATION_PROCESSES
            + _TEMPERATURE_PROCESSES
            + _BYPRODUCT_PROCESSES
            + _VDK_PROCESSES
            + _ACETALDEHYDE_PROCESSES
            + _KETO_ACID_PROCESSES
            + _H2S_PROCESSES
            + _MLF_PROCESSES
            + _MLF_GROWTH_PROCESSES
            + _BRETT_PROCESSES
            + _BRETT_GROWTH_PROCESSES
            + _POF_PROCESSES
            + _AMINO_ACID_PROCESSES
            + _AUTOLYSIS_PROCESSES
            + _AGING_PROCESSES
        ),
        modifier_factories=_WINE_FERMENTATION_MODIFIERS + _CARRYING_CAPACITY_MODIFIERS,
    ),
    "beer": Medium(
        name="beer",
        schema=beer_schema(),
        process_factories=(
            _PRIMARY_FERMENTATION_PROCESSES
            + _TEMPERATURE_PROCESSES
            + _BYPRODUCT_PROCESSES
            + _VDK_PROCESSES
            + _ACETALDEHYDE_PROCESSES
            + _H2S_PROCESSES
            + _HOPS_PROCESSES
            + _AGING_PROCESSES
        ),
        modifier_factories=_PRIMARY_FERMENTATION_MODIFIERS,
    ),
}


def get_medium(name: str) -> Medium:
    """Look up a registered :class:`Medium` by name."""
    try:
        return MEDIA[name]
    except KeyError:
        raise KeyError(f"Unknown medium {name!r}; known media: {sorted(MEDIA)}") from None
