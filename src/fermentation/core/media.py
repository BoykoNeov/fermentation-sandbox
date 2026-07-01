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
    h2s      hydrogen sulfide     g/L (sulfidic "rotten egg" off-aroma; produced-only,
                                 de-repressed at low nitrogen; carbon-free — decision D-29)
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
    AminoAcidAssimilation,
    ArrheniusTemperature,
    BiomassCarryingCapacity,
    DiacetylReduction,
    EsterSynthesis,
    EsterVolatilization,
    EthanolInactivation,
    FuselAlcoholsEhrlich,
    GrowthNitrogenLimited,
    HydrogenSulfideProduction,
    MalolacticCitrateMetabolism,
    MalolacticConversion,
    OenococcusDiacetylReduction,
    SugarUptakeToEthanolCO2,
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
            description="hydrogen sulfide (H2S) — 'rotten egg' sulfidic off-aroma; produced-only "
            "pool, de-repressed at low yeast-assimilable nitrogen; carbon-free (D-29)",
        ),
    ]


def wine_schema() -> StateSchema:
    """Wine state layout: a single lumped fermentable sugar slot, plus the wine-only
    charge-active acid + strong-cation slots the pH charge-balance solver reads
    (decision D-18), the free-SO₂ pool the molecular-SO₂ readout reads (decision D-22),
    the ``X_mlf`` malolactic-catalyst slot (decision D-23), the ``citrate`` slot
    *O. oeni* co-metabolises into MLF-derived diacetyl (decision D-31), and the dosed
    ``amino_acids`` pool the amino-acid ledger swap funds biomass from (decision D-32).

    These eight slots are appended to ``wine_schema`` only (not ``_common_specs``), so
    ``beer_schema`` is untouched — beer's pH is a phosphate-buffered different acid
    system with no sourced data yet, explicitly deferred. ``default=0.0`` is
    load-bearing: existing wine scenarios/tests that name no acids still compile (all
    eight → 0), and with acids, cation, SO₂, ``X_mlf``, ``citrate`` and ``amino_acids`` at 0
    the slots are inert — they
    contribute 0 to every conservation sum, so the validated core and its tests are
    untouched (prime directive #3). The acid/cation/SO₂ slots have no Process touching
    them in D-18/D-22; under D-23 :class:`~fermentation.core.kinetics.malolactic.\
    MalolacticConversion` depletes ``malic`` / grows ``lactic`` / evolves ``CO2`` *only
    when ``X_mlf`` is dosed* (and is disabled at the compile seam otherwise), so undosed
    wine runs keep a constant acid trajectory. ``X_mlf`` itself is inert in v1 (no Process
    grows or kills it) and carbon-free in ``total_carbon`` (constant ⇒ 0 drift); it enters
    the carbon ledger only when the later MLF-growth beat lands. pH is
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
            description="Oenococcus oeni biomass — dosed-but-inert MLF catalyst "
            "(scales the malolactic rate; no growth/death in v1, decision D-23)",
        ),
        VarSpec(
            "amino_acids",
            "g/L",
            default=0.0,
            description="assimilable amino-acid pool (dosed must input; represented as "
            "arginine). Carbon- AND nitrogen-bearing: the AminoAcidAssimilation swap funds "
            "a fraction of biomass from it, refunding sugar + ammonium N (decision D-32)",
        ),
    ]
    return StateSchema(specs)


def beer_schema() -> StateSchema:
    """Beer state layout: three sugars consumed sequentially.

    Glucose is taken up first, then maltose, then maltotriose — the order the
    ``components`` tuple records and the sugar-uptake Process will honour.
    """
    return StateSchema(
        _common_specs(
            VarSpec(
                "S",
                "g/L",
                size=3,
                description="fermentable sugars (sequential uptake)",
                components=("glucose", "maltose", "maltotriose"),
            )
        )
    )


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
#: Arrhenius temperature dependence scaling growth and uptake. The only structural
#: difference between the two media is the sugar vector (1 slot vs 3): beer's
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

#: Hydrogen-sulfide production (Milestone 2, decision D-29): the low-nitrogen "rotten egg"
#: off-aroma, one flux-linked producer gated by an *inverse*-nitrogen term. Kept as its own
#: isolable tuple (prime directive #3): a ProcessSet built without it is the prior core. Like
#: the ester/VDK/acetaldehyde pools (and unlike the *dosed* MLF organism), H₂S is intrinsic
#: yeast metabolism, so it is wired into BOTH media and runs on every default ferment. This is
#: the most isolable beat in the model: H₂S is CARBON-FREE (on no conservation ledger) and the
#: Process touches ONLY ``h2s`` while merely *reading* ``X``/``S``/``N`` — so disabling it
#: leaves the RHS of every other column byte-for-byte identical (no ``h2s`` consumer exists to
#: feed anything back); the integrated trajectory then differs only by a ~1e-7 adaptive-solver
#: mesh artifact, cleaner than the acetaldehyde buffer's *genuine* second-order E→viability
#: coupling (D-27).
#: No tier headline either: it writes a fresh pool nothing reads, so no other column's
#: structural tier drops (contrast the D-26 ``CO2`` / D-27 ``E`` cases). Params live in the
#: shared, medium-agnostic ``hydrogen_sulfide.yaml`` (sulfate-reduction is generic yeast
#: metabolism). SCOPE (v1): produced-only (the CO₂-stripping sink is the deferred follow-up,
#: the ester D-19→D-20 precedent), so ``h2s`` is cumulative-produced (overstates residual);
#: and the cross-must YAN lever is muted by the upstream N→0 stripping gap (decision D-29).
_H2S_PROCESSES: tuple[Callable[[], Process], ...] = (HydrogenSulfideProduction,)

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
_MLF_PROCESSES: tuple[Callable[[], Process], ...] = (
    MalolacticConversion,
    MalolacticCitrateMetabolism,
    OenococcusDiacetylReduction,
)

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
_AMINO_ACID_PROCESSES: tuple[Callable[[], Process], ...] = (AminoAcidAssimilation,)

#: Wine growth/uptake Arrhenius modifiers (decision D-32). Identical to
#: :data:`_PRIMARY_FERMENTATION_MODIFIERS` except the growth Arrhenius *also* scales the
#: amino-acid swap (``for_growth`` extra target), so the swap's carbon/nitrogen refunds carry
#: the same temperature factor as growth's draw (see the correctness coupling above). Beer keeps
#: the plain :data:`_PRIMARY_FERMENTATION_MODIFIERS` (no amino-acid pool). ``for_uptake`` is
#: unchanged — the swap tracks *growth*, not the fermentative sugar-uptake flux.
_WINE_FERMENTATION_MODIFIERS: tuple[Callable[[], RateModifier], ...] = (
    lambda: ArrheniusTemperature.for_growth(AminoAcidAssimilation.name),
    ArrheniusTemperature.for_uptake,
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
            + _BYPRODUCT_PROCESSES
            + _VDK_PROCESSES
            + _ACETALDEHYDE_PROCESSES
            + _H2S_PROCESSES
            + _MLF_PROCESSES
            + _AMINO_ACID_PROCESSES
        ),
        modifier_factories=_WINE_FERMENTATION_MODIFIERS + _CARRYING_CAPACITY_MODIFIERS,
    ),
    "beer": Medium(
        name="beer",
        schema=beer_schema(),
        process_factories=(
            _PRIMARY_FERMENTATION_PROCESSES
            + _BYPRODUCT_PROCESSES
            + _VDK_PROCESSES
            + _ACETALDEHYDE_PROCESSES
            + _H2S_PROCESSES
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
