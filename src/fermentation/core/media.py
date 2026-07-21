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
    ethyl_acetate    ethyl acetate    g/L (solventy/nail-polish acetate ester — D-96)
    isoamyl_acetate  isoamyl acetate  g/L (banana acetate ester — D-96)
    ethyl_hexanoate  ethyl hexanoate  g/L (apple/pineapple fatty-acid ethyl ester — D-96)
    isoamyl_alcohol / active_amyl_alcohol / isobutanol / propanol / 2_phenylethanol
                                 g/L (the five Ehrlich higher alcohols, produced-only; each its
                                 own molecule since D-99 — the split of the former lumped
                                 ``fusels`` pool, which weighted and perceived all five as
                                 isoamyl alcohol)
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
``X_dead``, ``Gly``, ``Byp``, the three ester pools, the five higher-alcohol pools and
the VDK pools (``acetolactate``/``diacetyl``/``butanediol``) start at zero at pitch and
are only accumulated by the kinetics, so they declare a default initial of 0
(`VarSpec.default`) and need not be named at every initial-condition call site. The
ester and higher-alcohol pools are filled by the Tier-2 byproduct Processes wired below;
the three VDK pools by the diacetyl-pathway Processes (decision D-26).
Under **decision D-19 (option a1)** those Processes route the aroma carbon *out of
``S``* and ``total_carbon`` weights the pools (each ester as ITSELF since D-96, each higher
alcohol as ITSELF since D-99), so they are real carbon-accounted state alongside
``Gly``/``Byp`` — not diagnostic re-expressions. The former ``Byp`` double-count (it once
lumped higher alcohols) is resolved by carving them out of ``Y_byproduct_sugar``; the draw
touches
only ``S`` (never ``E``/``CO2``), so turning the byproducts on perturbs the core only
by the trace sugar they consume. See ``docs/plans/milestone-2-tasks.md`` and the
``kinetics.byproducts`` module docstring.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from fermentation.core.kinetics import (
    AcetaldehydeBridgedCondensation,
    AcetaldehydeProduction,
    AcetaldehydeReduction,
    AcetolactateDecarboxylation,
    AcetolactateExcretion,
    AlphaKetobutyrateExcretion,
    AlphaKetobutyrateReassimilation,
    AlphaKetoglutarateExcretion,
    AlphaKetoglutarateReassimilation,
    AminoAcidAssimilation,
    AnthocyaninFading,
    ArrheniusTemperature,
    AutolyticHydrogenSulfide,
    AutolyticMercaptan,
    BiomassCarryingCapacity,
    BrettDeath,
    BrettDecarboxylation,
    BrettEthanolToxicity,
    BrettGrowth,
    BrettVinylphenolReduction,
    Caramelization,
    ColemanQuadraticDeathTemperature,
    DiacetylReduction,
    EllagitanninOxidation,
    EsterHydrolysis,
    EsterSynthesis,
    EsterVolatilization,
    EthanolInactivation,
    EthylAcetateEsterification,
    EthylHexanoateHydrolysis,
    FuselAlcoholsEhrlich,
    FuselAminoAcidReroute,
    GrowthNitrogenLimited,
    HydrogenSulfideProduction,
    HydrogenSulfideVolatilization,
    IsoAlphaAcidLoss,
    MaillardBrowning,
    MaillardStrecker,
    MalolacticCitrateMetabolism,
    MalolacticConversion,
    MalolacticDeath,
    MalolacticGrowth,
    MalolacticSenescence,
    OakExtraction,
    OenococcusDiacetylReduction,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    PrecursorNonEhrlichFates,
    PyruvateExcretion,
    PyruvateReassimilation,
    SMMHydrolysis,
    SotolonAldolCondensation,
    StreckerDegradation,
    SugarUptakeToEthanolCO2,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
    TanninEthylTanninCondensation,
    TanninSelfPolymerization,
    TemperatureRamp,
    ThermalAnthocyaninFade,
    YeastAutolysis,
    YeastPOFDecarboxylation,
)
from fermentation.core.kinetics.carbon_routing import (
    ESTER_SPECS,
    FUSEL_SPECS,
    VALINE_LABEL_TRACERS,
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
        #: The THREE single-molecule ester pools (decision D-96), derived from the canonical
        #: ``ESTER_SPECS`` registry so a fourth ester is one entry there — never a hand-edit
        #: here that could drift from the carbon ledger or the OAV aroma set. Each replaced a
        #: share of the pre-D-96 lumped ``esters`` pool, which was weighted as ethyl acetate
        #: but *perceived* as isoamyl acetate.
        *(VarSpec(spec.pool, "g/L", default=0.0, description=spec.note) for spec in ESTER_SPECS),
        #: The FIVE single-molecule Ehrlich higher-alcohol pools (decision D-99), derived from
        #: the canonical ``FUSEL_SPECS`` registry for the same reason the esters are: a sixth
        #: alcohol is one entry there, never a hand-edit here that could drift from the carbon
        #: ledger or the OAV aroma set. Together they replace the pre-D-99 lumped ``fusels``
        #: pool, which weighted AND perceived all five as isoamyl alcohol. Unlike the esters
        #: these have no headspace twin — they are not stripped in this model.
        *(VarSpec(spec.pool, "g/L", default=0.0, description=spec.note) for spec in FUSEL_SPECS),
        #: The TWO valine-label tracer slots (decision D-115): the valine-derived part of the
        #: isoamyl alcohol pool and of the isoamyl acetate pool, in g/L of the labelled molecule,
        #: so Rollero's enrichment is ``tracer / bulk`` directly. **Not** the provenance metadata
        #: D-1 excludes from the state floats — a ¹³C isotopologue concentration is a conserved
        #: extensive quantity that flows and integrates like any other pool (see
        #: :class:`~fermentation.core.kinetics.carbon_routing.LabelTracer`). They carry carbon
        #: weight **zero** in ``total_carbon``, because each is a sub-quantity of a pool already
        #: weighted there; weighting them would double-count every labelled gram.
        *(
            VarSpec(
                tracer.tracer_pool,
                "g/L",
                default=0.0,
                description=(
                    f"valine-derived {tracer.bulk_pool} (D-115 label tracer) — a SUB-QUANTITY "
                    f"of {tracer.bulk_pool}, not additional mass; carbon-ledger weight 0"
                ),
            )
            for tracer in VALINE_LABEL_TRACERS
        ),
        #: Each ester's headspace twin (decision D-20, generalised per-ester at D-96): the
        #: CO2-stripping sink moves liquid ester carbon here. A pool and its twin share ONE
        #: molecule's carbon weight, which is what makes the strip carbon-neutral — so a
        #: SINGLE shared gas pool is impossible once the esters differ (C4/C7/C8).
        *(
            VarSpec(
                spec.gas_pool,
                "g/L",
                default=0.0,
                description=f"{spec.pool} lost to the headspace by CO2 stripping "
                "(volatilized; carbon-bookkeeping pool, decisions D-20/D-96)",
            )
            for spec in ESTER_SPECS
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
        VarSpec(
            "o2",
            "g/L",
            default=0.0,
            description="dissolved oxygen — the OXIDATIVE-aging substrate (decision D-71). Dosed "
            "post-ferment by add_oxygen (bottle ingress / micro-oxygenation); drawn down by the "
            "always-on O₂ sinks OxidativeAcetaldehyde (→ acetaldehyde, the 'sherry'/oxidised note) "
            "and PhenolicBrowning (→ brown pigment, D-74), plus wine's SulfiteOxidation. "
            "Carbon-free and on NO ledger (like h2s/iso_alpha). Default 0 ⇒ an un-oxygenated "
            "(reductive) aging is byte-for-byte the ester-hydrolysis-only case",
        ),
        VarSpec(
            "A420",
            "AU",
            default=0.0,
            description="oxidative-browning index — absorbance at 420 nm, the standard measure of "
            "wine/beer browning (decision D-74). Accumulated by PhenolicBrowning as dissolved O₂ "
            "oxidises phenolics to brown quinone/melanoidin pigment (the gold→amber→brown of aged "
            "white wine; oxidative darkening in beer). An OPTICAL INDEX (dimensionless AU, 1 cm "
            "path), NOT a pigment mass — so carbon-free and on NO ledger (like o2/iso_alpha), and "
            "its carbon (from untracked phenols) is sidestepped by construction. Cumulative and "
            "monotonic (d(A420)/dt ≥ 0). Default 0 ⇒ a reductive/un-oxygenated aging is byte-for-"
            "byte the case without browning",
        ),
    ]


def _oak_specs() -> list[VarSpec]:
    """The oak-extraction axis slots — SHARED by wine and barrel-beer (decisions D-77/D-78/D-86).

    The barrel/chip aroma-extractive aging axis (decision D-77), plus the ellagitannin BRIDGE
    (decision D-78) and the ``furaneol`` caramel furanone (decision D-94). FIVE extracted AROMA
    pools + the ellagitannin TASTE pool (six rising toward their ceilings) + six SET-AND-HOLD
    ceiling slots (the ``cation_charge`` idiom — state written ONLY by the ``add_oak`` verb, never
    by a Process). The aroma five are a SEPARATE, non-oxidative axis (draw no O₂); ellagitannin
    bridges to the ``o2`` sub-axis — :class:`EllagitanninOxidation` (D-78) consumes it to scavenge
    O₂ (protecting the beverage). All twelve slots are OFF EVERY LEDGER
    (exogenous WOOD-derived mass, like the hop-derived ``iso_alpha``, D-64), so they never perturb
    ``total_carbon``/``total_mass``/``total_nitrogen``. ``default=0`` ⇒ an un-oaked beverage carries
    no ceiling and both oak Processes are byte-for-byte inert (the ceiling ≤ 0 guard).

    **Medium-agnostic (D-86 — barrel-beer oak).** D-77/D-78 wired these into ``wine_schema`` ONLY
    (the barrel-beer extension was explicitly deferred as trivial); D-86 factored them into this
    helper so both :func:`wine_schema` (unchanged position — wine stays byte-for-byte) and
    :func:`beer_schema` carry the identical axis. Oak extraction is a WOOD property (``oak.yaml``
    yields are g-extractive-per-g-oak, matrix-independent); only the OAV perception thresholds are
    matrix-specific (``threshold_<compound>_beer`` vs ``_wine``, ``sensory.yaml``). The GRAPE colour
    axis (anthocyanin/tannin condensation + fade) stays wine-only — it is grape chemistry, not oak.
    """
    return [
        VarSpec(
            "whiskey_lactone",
            "g/L",
            default=0.0,
            description="whiskey lactone (β-methyl-γ-octalactone, cis+trans lumped) — the "
            "'coconut' oak-lactone note (decision D-77), LIGHT-toast dominant. Produced-only: "
            "OakExtraction rises it toward whiskey_lactone_ceiling (oak diffusion). Off every "
            "ledger (wood-derived). Read by the OAV lens (threshold_whiskey_lactone_<medium>)",
        ),
        VarSpec(
            "vanillin",
            "g/L",
            default=0.0,
            description="vanillin — the 'vanilla' oak extractive (decision D-77), MEDIUM-toast "
            "peak (lignin thermal release). Produced-only: OakExtraction rises it toward "
            "vanillin_ceiling. Off every ledger (wood-derived). OAV lens "
            "(threshold_vanillin_<medium>)",
        ),
        VarSpec(
            "guaiacol",
            "g/L",
            default=0.0,
            description="guaiacol — the 'smoky/toasty' oak extractive (decision D-77), HEAVY-toast "
            "dominant (lignin pyrolysis). DISTINCT from the Brett 4-ethylguaiacol (D-55). "
            "Produced-only: OakExtraction rises it toward guaiacol_ceiling. Off every ledger "
            "(exogenous wood-derived). Read by the OAV lens (threshold_guaiacol_<medium>)",
        ),
        VarSpec(
            "eugenol",
            "g/L",
            default=0.0,
            description="eugenol — the 'clove/spice' oak extractive (decision D-77), HEAVY-toast "
            "(co-varies with guaiacol). Produced-only: OakExtraction rises it toward "
            "eugenol_ceiling. Off every ledger (wood-derived). OAV lens "
            "(threshold_eugenol_<medium>)",
        ),
        VarSpec(
            "furaneol",
            "g/L",
            default=0.0,
            description="furaneol (HDMF, 4-hydroxy-2,5-dimethyl-3(2H)-furanone) — the "
            "'caramel/toffee' oak/bourbon furanone (decision D-94), a thermal sugar-degradation "
            "product RISING with toast (co-varies with guaiacol). Produced-only: OakExtraction "
            "rises it toward furaneol_ceiling (wood diffusion + the D-93 bourbon spirit soak-back "
            "bump). Off every ledger (wood/spirit-derived, so no collision with the on-ledger "
            "D-88 caramelization melanoidin). Read by the OAV lens (threshold_furaneol_<medium>)",
        ),
        VarSpec(
            "whiskey_lactone_ceiling",
            "g/L",
            default=0.0,
            description="SET-AND-HOLD saturation ceiling for whiskey_lactone (decision D-77): "
            "oak_gpl × oak_yield_whiskey_lactone_<toast>, written ONLY by the add_oak verb "
            "(constant state no Process touches, the cation_charge idiom). OakExtraction reads it. "
            "Off every ledger. Default 0 ⇒ no oak ⇒ inert",
        ),
        VarSpec(
            "vanillin_ceiling",
            "g/L",
            default=0.0,
            description="SET-AND-HOLD saturation ceiling for vanillin (decision D-77): "
            "oak_gpl × oak_yield_vanillin_<toast>, written ONLY by add_oak. Off every ledger. "
            "Default 0 ⇒ inert",
        ),
        VarSpec(
            "guaiacol_ceiling",
            "g/L",
            default=0.0,
            description="SET-AND-HOLD saturation ceiling for guaiacol (decision D-77): "
            "oak_gpl × oak_yield_guaiacol_<toast>, written ONLY by add_oak. Off every ledger. "
            "Default 0 ⇒ inert",
        ),
        VarSpec(
            "eugenol_ceiling",
            "g/L",
            default=0.0,
            description="SET-AND-HOLD saturation ceiling for eugenol (decision D-77): "
            "oak_gpl × oak_yield_eugenol_<toast>, written ONLY by add_oak. Off every ledger. "
            "Default 0 ⇒ inert",
        ),
        VarSpec(
            "furaneol_ceiling",
            "g/L",
            default=0.0,
            description="SET-AND-HOLD saturation ceiling for furaneol (decision D-94): "
            "oak_gpl × oak_yield_furaneol_<toast> (+ the D-93 spirit_soak_furaneol_<spirit> bump "
            "for an ex-bourbon barrel), written ONLY by add_oak. Off every ledger. "
            "Default 0 ⇒ inert",
        ),
        # Ellagitannin — the BRIDGE extractive (decision D-78). Unlike the four aroma extractives
        # above (pure diffusion axis, O₂-orthogonal), ellagitannin is DYNAMIC: OakExtraction rises
        # it toward its ceiling (diffusion in), AND EllagitanninOxidation draws its share of the
        # shared o2 budget and CONSUMES it (the oak O₂-scavenging PROTECTION, the D-78 spine). It is
        # a TASTE extractive — astringency, read out by analysis.astringency_series (the
        # iso_alpha/IBU precedent), NOT the D-67 OAV aroma lens. Both slots are OFF EVERY LEDGER
        # (wood-derived, the iso_alpha precedent), so neither Process perturbs
        # total_carbon/total_mass/total_nitrogen. default=0 ⇒ an un-oaked beverage carries no
        # ceiling and both Processes are byte-for-byte inert.
        VarSpec(
            "ellagitannin",
            "g/L",
            default=0.0,
            description="ellagitannin — oak's hydrolysable TANNIN (decision D-78), the ASTRINGENCY "
            "extractive AND an O₂ scavenger, LIGHT-toast dominant (thermolabile — degraded by "
            "toasting). Dynamic: OakExtraction rises it toward ellagitannin_ceiling (diffusion), "
            "EllagitanninOxidation consumes it as it scavenges dissolved O₂ (protecting the "
            "beverage). Off every ledger (wood-derived). A TASTE — read by "
            "analysis.astringency_series (mg/L tannin), NOT the OAV odor lens (the iso_alpha/IBU "
            "exclusion)",
        ),
        VarSpec(
            "ellagitannin_ceiling",
            "g/L",
            default=0.0,
            description="SET-AND-HOLD saturation ceiling for ellagitannin (decision D-78): "
            "oak_gpl × oak_yield_ellagitannin_<toast>, written ONLY by the add_oak verb (the "
            "cation_charge idiom). OakExtraction reads it (never written by a Process). Off every "
            "ledger. Default 0 ⇒ no oak ⇒ inert",
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
            "methanethiol",
            "g/L",
            default=0.0,
            description="methanethiol (CH3SH) — the carbon-bearing reductive off-aroma. Named for "
            "the one molecule it contains since D-110: it was `mercaptans` through D-109, a plural "
            "asserting a thiol MIXTURE the mass balance never held (nothing in the model makes "
            "ethanethiol or any other thiol). AutolyticMercaptan fills it as a yield on the "
            "autolysis flux, drawing carbon from methionine (D-100 — the actual sulfur-bearing "
            "precursor, not the arginine lump) at demethiolation's real 1:1 stoichiometry, booking "
            "the C4 co-product to alpha_ketobutyrate and deaminating the nitrogen to N (D-107); "
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
        VarSpec(
            "alpha_ketobutyrate",
            "g/L",
            default=0.0,
            description="excreted overflow alpha-ketobutyrate (2-oxobutyrate; C4 keto-acid; "
            "excreted-then-reassimilated, the pyruvate/alpha-KG structure). THE KETO-ACID NODE "
            "(D-107): unlike its two siblings this pool has a CONSUMER — it is the C4 half of "
            "sotolon's aldol (SotolonAldolCondensation) — and a second PRODUCER, the D-45 "
            "mercaptan, whose 2-oxobutyrate by-product had nowhere to go until this slot existed. "
            "AlphaKetobutyrateExcretion draws it threonine:sugar by threonine's depletion gate "
            "(Crepin's 19% exogenous / 81% de novo, reproduced not asserted); the flux-linked "
            "reassimilation freezes a persistent finished-wine residual at dryness, which is what "
            "the bottle-aging aldol then eats",
        ),
        VarSpec(
            "methional",
            "g/L",
            default=0.0,
            description="methional — the 'cooked-potato' Strecker aldehyde "
            "of methionine, the principal OXIDATIVE off-note of aged wine/stale beer (decision "
            "D-75). Produced-only: StreckerDegradation forms it as dissolved O₂ (via the phenol "
            "autoxidation quinones) oxidatively deaminates + decarboxylates amino acids; carbon "
            "drawn from METHIONINE, its real precursor (D-100), the nitrogen deaminated to N, "
            "one CO₂ "
            "released per aldehyde. Read by the D-67 OAV lens (threshold_methional_wine)",
        ),
        VarSpec(
            "phenylacetaldehyde",
            "g/L",
            default=0.0,
            description="phenylacetaldehyde — the 'honey/floral' Strecker aldehyde of "
            "phenylalanine (decision D-75), the pleasant-valence counterpart to methional from the "
            "SAME quinone-driven Strecker route. Produced-only, carbon from amino_acids + CO₂, "
            "nitrogen deaminated to N. Read by the OAV lens (threshold_phenylacetaldehyde_wine)",
        ),
    ]
    # Oak extraction axis (D-77 aroma four + D-78 ellagitannin) — SHARED with barrel-beer (D-86),
    # factored into _oak_specs() and inserted here at its ORIGINAL wine position so wine stays
    # byte-for-byte identical. See _oak_specs for the full axis story; beer_schema appends the
    # same helper. The GRAPE colour axis below stays wine-only (grape chemistry, not oak).
    specs += _oak_specs()
    specs += [
        # Tannin–anthocyanin condensation — the red-wine colour-stabilization +
        # astringency-softening
        # aging axis (decision D-79). Two GRAPE-derived must-input pools (the hydroxycinnamic_gpl
        # precedent): free monomeric anthocyanin (bleachable red pigment) + condensed grape tannin
        # (harsh young astringency). TanninAnthocyaninCondensation (bilinear) consumes BOTH into a
        # stable polymeric pigment — a SEPARATE, non-oxidative GRAPE axis: it draws NO o2 (unlike
        # every oxidative sink) and reads NO oak pool (grape condensed tannin ≠ oak hydrolysable
        # ellagitannin). The polymeric pigment is a POST-HOC readout (anthocyanin₀ − anthocyanin),
        # NOT
        # a slot (the A420 discriminator — anthocyanin's single fate makes it reconstructible). Both
        # slots are OFF EVERY LEDGER (grape-derived, the iso_alpha/ellagitannin precedent), so the
        # Process perturbs nothing conserved. default=0 ⇒ a white / no-tannin wine carries neither
        # and
        # the Process is byte-for-byte inert (doubly substrate-gated). Wine. Read as TASTE/COLOUR by
        # analysis.astringency_series / polymeric_pigment_series / color_series, NOT the OAV odor
        # lens.
        VarSpec(
            "anthocyanin",
            "g/L",
            default=0.0,
            description="free monomeric anthocyanin — the bright, bleachable purple-red grape "
            "pigment (decision D-79). GRAPE must input (default 0 ⇒ white wine). Has TWO fates as "
            "the wine ages: condensed into stable polymeric_pigment (TanninAnthocyaninCondensation "
            "D-79 / AcetaldehydeBridgedCondensation D-80 — the young purple → aged brick-red "
            "evolution) AND oxidatively faded to colourless faded_anthocyanin (AnthocyaninFading "
            "D-81 — the irreversible bleaching loss). Off every ledger (grape-derived, the "
            "iso_alpha/ellagitannin precedent). Read as COLOUR by analysis.color_series (free "
            "anthocyanin + polymeric_pigment; faded is colourless) and, SO₂/pH-masked, by "
            "analysis.observed_color_series (the reversible Somers bleaching readout, D-82), NOT "
            "the OAV odor lens (colour is not an aroma)",
        ),
        VarSpec(
            "tannin",
            "g/L",
            default=0.0,
            description="condensed grape (flavan-3-ol, skin/seed) tannin — the harsh young-red "
            "astringency (decision D-79). GRAPE must input (default 0). A DIFFERENT molecule from "
            "oak's hydrolysable ellagitannin (D-78): this is the grape `tannin` the D-78 namespace "
            "note left free. Consumed by TanninAnthocyaninCondensation (with anthocyanin) into "
            "soft polymeric pigment, so astringency SOFTENS. Off every ledger (grape-derived). "
            "Read "
            "as TASTE by analysis.astringency_series (mg/L, summed WITH oak ellagitannin — both "
            "harsh), NOT the OAV odor lens (astringency is a taste, the iso_alpha/IBU exclusion)",
        ),
        # Acetaldehyde-bridged (ethylidene) condensation — the SPLIT-LEDGER colour beat (decision
        # D-80), the D-79-deferred second pigment-formation pathway. AcetaldehydeBridgedCondensation
        # (trilinear in acetaldehyde × anthocyanin × tannin) bridges grape tannin to anthocyanin
        # with an acetaldehyde-derived ethylidene linker (tannin–ethyl–anthocyanin). Unlike the D-79
        # direct route (moves nothing conserved), the bridged route consumes ON-LEDGER acetaldehyde
        # (its carbon borrowed from E at D-71), so this `ethyl_bridge` slot CAPTURES that carbon
        # on-ledger — weighted at cf(ethylidene) in total_carbon — instead of letting it vanish into
        # the off-ledger grape pigment (the "split ledger": grape bulk off-ledger, acetaldehyde-
        # derived bridge on it). The FIRST aging colour slot ON the carbon ledger. Filled BY the
        # Process (starts at 0, no must input — it accumulates the bridged acetaldehyde carbon);
        # default=0 ⇒ inert until the Process fires (needs anthocyanin + tannin + acetaldehyde all
        # present + begin_aging). Wine.
        VarSpec(
            "ethyl_bridge",
            "g/L",
            default=0.0,
            description="acetaldehyde-derived ethylidene bridge carbon (—CH(CH₃)—) locked into "
            "polymeric pigment by AcetaldehydeBridgedCondensation (decision D-80). ON the carbon "
            "ledger (weighted at cf(ethylidene) in total_carbon): it captures the acetaldehyde "
            "carbon the bridged route consumes — borrowed from ethanol at D-71 — so it does NOT "
            "vanish into the off-ledger grape pigment (the SPLIT-LEDGER accounting). Filled by the "
            "Process (no must input, starts 0); an integrated slot, not a readout, because "
            "acetaldehyde has competing fates (the A420 discriminator). NOT read by any sensory "
            "lens (colour is captured via anthocyanin drawdown; this is carbon bookkeeping)",
        ),
        # Polymeric pigment PROMOTED to an integrated slot + the colourless fade sink — the SO₂/pH
        # anthocyanin-bleaching beat (decision D-81). D-79/D-80 kept the stable pigment a POST-HOC
        # readout (anthocyanin₀ − anthocyanin) because condensation was anthocyanin's SOLE fate.
        # D-81's AnthocyaninFading gives anthocyanin a SECOND, irreversible fate (oxidative
        # degradation → colourless), so that reconstruction identity breaks (it would wrongly count
        # faded anthocyanin as pigment) and the pigment MUST become a real slot (the A420
        # discriminator, D-74). `polymeric_pigment` is now filled by BOTH condensation routes
        # (direct D-79 + bridged D-80, d/dt = +r each); `faded_anthocyanin` is filled by
        # AnthocyaninFading. Both OFF EVERY LEDGER (grape-derived colour-equivalents, the
        # anthocyanin/tannin precedent), both filled BY their Processes (no must input, start 0),
        # both wine-only. Together they close the three-slot colour identity anthocyanin +
        # polymeric_pigment + faded_anthocyanin ≡ anthocyanin₀ (holds by construction — the d/dt
        # terms sum to zero — NOT via assert_conserved, whose weights are 0 for these off-ledger
        # slots).
        VarSpec(
            "polymeric_pigment",
            "g/L",
            default=0.0,
            description="stable polymeric pigment (tannin–anthocyanin condensate) — the "
            "SO₂/pH-STABLE aged red colour form (decision D-81, promoted from the D-79 post-hoc "
            "readout). Filled by BOTH TanninAnthocyaninCondensation (direct, D-79) and "
            "AcetaldehydeBridgedCondensation (bridged, D-80), each writing +r in "
            "anthocyanin-equivalents. An integrated SLOT, not a readout, because D-81's "
            "AnthocyaninFading gives anthocyanin a second fate → the anthocyanin₀ − anthocyanin "
            "reconstruction no longer isolates the pigment (the A420 discriminator, D-74). Off "
            "every ledger (grape-derived colour-equivalent, the anthocyanin/tannin precedent); "
            "starts 0, no must input. Read as COLOUR by analysis.color_series / "
            "polymeric_pigment_series and counted FULL (SO₂/pH-resistant) by "
            "analysis.observed_color_series (the bleach-RESISTANT fraction — the colour-stability "
            "payoff, D-82), NOT the OAV odor lens (colour is not an aroma)",
        ),
        VarSpec(
            "faded_anthocyanin",
            "g/L",
            default=0.0,
            description="colourless anthocyanin-degradation products — the IRREVERSIBLE fade sink "
            "(decisions D-81/D-83). Filled by TWO routes into one pool: AnthocyaninFading "
            "(O₂-coupled oxidative bleaching, D-81: r = k_fade·f(T)·o2·[anthocyanin], drawing the "
            "shared o2 pool) AND ThermalAnthocyaninFade (O₂-INDEPENDENT thermal/hydrolytic fade, "
            "D-83: r = k_thermal·f(T)·[anthocyanin], no o2 — so a sealed/sulfited/anaerobic red "
            "still fades, SO₂ giving no protection). Both capture the free monomeric anthocyanin "
            "lost to fading so it is NOT double-counted as pigment. Together they are the second "
            "anthocyanin fate that makes analysis.color_series GENUINELY decline (young bleachable "
            "colour is lost; the stable polymeric_pigment survives — the colour-stability payoff). "
            "Off every ledger (grape-derived, the anthocyanin precedent); starts 0, no must input, "
            "wine-only. NOT read by any sensory lens — it is colourless (the whole point), so it "
            "carries no colour and no odor",
        ),
    ]
    # Non-oxidative THERMAL Strecker aldehydes + sotolon (decision D-87), appended last: the four
    # NEW aroma pools MaillardStrecker produces from residual sugar + amino acids + heat with NO O₂
    # (the sweet-wine / Madeira / baked-wine suite). methional + phenylacetaldehyde (D-75, above)
    # are SHARED with this route — same molecules — so only these four are new. Carbon-bearing
    # (booked from each aldehyde's OWN precursor, D-100, deaminated to N), so on total_carbon like
    # the D-75 pair.
    # Wine-only. Read by the D-67 OAV lens against their own thresholds (threshold_<pool>_wine).
    specs += [
        VarSpec(
            "2_methylbutanal",
            "g/L",
            default=0.0,
            description="2-methylbutanal — the 'malty/almond' branched-chain Strecker aldehyde of "
            "isoleucine (decision D-87). Produced-only by MaillardStrecker: residual sugar forms "
            "α-dicarbonyls that deaminate + decarboxylate isoleucine WITH NO O₂ (the thermal "
            "mirror of the D-75 oxidative route); carbon from ISOLEUCINE, its real precursor "
            "(D-100), "
            "nitrogen deaminated to N, one CO₂ released. Read by the OAV lens "
            "(threshold_2_methylbutanal_wine)",
        ),
        VarSpec(
            "3_methylbutanal",
            "g/L",
            default=0.0,
            description="3-methylbutanal — the 'malty/dark-chocolate' branched-chain Strecker "
            "aldehyde of leucine (decision D-87), typically the most prominent thermal/staling "
            "branched-chain aldehyde. Produced-only by MaillardStrecker (sugar+heat, no O₂); "
            "carbon from amino_acids + CO₂, nitrogen deaminated to N. Read by the OAV lens "
            "(threshold_3_methylbutanal_wine)",
        ),
        VarSpec(
            "2_methylpropanal",
            "g/L",
            default=0.0,
            description="2-methylpropanal (isobutyraldehyde) — the 'malty/grainy' Strecker "
            "aldehyde of valine (decision D-87). Produced-only by MaillardStrecker (sugar+heat, "
            "no O₂); "
            "carbon from amino_acids + CO₂, nitrogen deaminated to N. Read by the OAV lens "
            "(threshold_2_methylpropanal_wine)",
        ),
        VarSpec(
            "sotolon",
            "g/L",
            default=0.0,
            description="sotolon (4,5-dimethyl-3-hydroxy-2(5H)-furanone) — the 'curry/maple/nutty' "
            "furanone marker of botrytized sweet wine (Sauternes), vin jaune, aged Port and "
            "Madeira (decision D-87). Produced-only by MaillardStrecker but NOT a decarboxylation "
            "Strecker aldehyde (a threonine/acetaldehyde aldol furanone), so it carries NO CO₂ "
            "term; carbon booked from THREONINE, its real precursor (D-100; its 2 "
            "acetaldehyde-derived "
            "carbons lumped in). Trace by mass but potent. Read by the OAV lens "
            "(threshold_sotolon_wine)",
        ),
    ]
    # Caramelization melanoidin carbon-park (decision D-88; medium-agnostic D-90), appended: the
    # brown thermal-browning polymer Caramelization forms by consuming residual sugar (the
    # O₂-independent mirror of PhenolicBrowning D-74). The FIRST aging pool that holds consumed
    # core-S carbon, so — unlike the off-ledger oak/colour lumps — it is ON total_carbon
    # (sugar → melanoidin closes exactly). It raises the SAME A420 index D-74 accumulates (not a
    # new observable). Sugar-only (nitrogen-free — caramelization, not Maillard). Now in BOTH media
    # (D-90: beer thermal browning — the beer_schema counterpart is appended in beer_schema()).
    specs.append(
        VarSpec(
            "melanoidin",
            "g/L",
            default=0.0,
            description="melanoidin — the brown caramelization polymer (caramelan stand-in "
            "C12H18O9) Caramelization forms from residual sugar by HEAT with NO O₂ (decision D-88, "
            "the O₂-independent thermal mirror of PhenolicBrowning D-74). A carbon-park pool (the "
            "debris/glucan precedent): ON total_carbon (it holds the consumed core-S carbon, so "
            "sugar → melanoidin closes exactly), unlike the off-ledger oak/colour lumps. "
            "Sugar-only (nitrogen-free — caramelization, not amino-acid Maillard). Raises the "
            "shared A420 "
            "browning index (read by analysis.a420); the melanoidin MASS itself is not a sensory "
            "pool. In BOTH media (D-90: beer thermal browning)",
        )
    )
    # N-bearing Maillard melanoidin carbon+nitrogen-park (decision D-89), appended last: the brown
    # amino-acid-incorporating thermal-browning polymer MaillardBrowning forms by consuming residual
    # sugar AND amino acids (the N-incorporating browning branch D-88 deferred). It holds consumed
    # core-S carbon AND amino-acid carbon+nitrogen, so — unlike the off-ledger oak/colour lumps and
    # unlike the nitrogen-free caramelan — it is ON total_carbon AND total_nitrogen (the FIRST
    # non-biomass, non-arginine species on the nitrogen ledger; sugar + amino_acids → this pool
    # closes both exactly). It raises the SAME A420 index D-74/D-88 accumulate. Wine-only v1.
    specs.append(
        VarSpec(
            "maillard_melanoidin",
            "g/L",
            default=0.0,
            description="maillard_melanoidin — the brown N-bearing Maillard melanoidin polymer "
            "(glucose–glycine stand-in C8H12O5N, molar C:N ≈ 8:1) MaillardBrowning forms from "
            "residual sugar + amino acids by HEAT with NO O₂ (decision D-89, the "
            "amino-acid-incorporating browning branch D-88's sugar-only Caramelization deferred). "
            "A carbon+nitrogen-park pool: ON total_carbon AND total_nitrogen (the FIRST "
            "non-biomass, non-arginine species on the nitrogen ledger — it RETAINS the amino-acid "
            "nitrogen, what makes a Maillard melanoidin nitrogenous; sugar + amino_acids → this "
            "pool closes both ledgers exactly). Raises the shared A420 browning index (read by "
            "analysis.a420); the MASS itself is not a sensory pool. Wine-only",
        )
    )
    # The D-100 speciation of the lumped `amino_acids` pool (arginine stand-in, D-32) into
    # single-molecule amino acids, appended last so existing wine slot indices are unchanged. The
    # existing `amino_acids` slot ABOVE is retained as the ARGININE pool (its representative species
    # was always arginine); these seven siblings give every other consumed amino acid its own slot,
    # each carbon- AND nitrogen-weighted by its own molecule (chemistry.py). The five Ehrlich
    # precursors (leucine/isoleucine/valine/threonine/phenylalanine) map onto the D-99 fusels and
    # the D-87 thermal Strecker aldehydes; phe + met feed the D-75 oxidative Strecker aldehydes and
    # methionine the D-45 mercaptan; `amino_acids_generic` (glutamine stand-in) is the bucket of
    # every other assimilable amino acid the identity-agnostic yeast/MLF/Brett swaps draw alongside
    # arginine (D-100). Proline is NOT tracked — not assimilated anaerobically (excluded from YAN).
    # All default=0: an undosed / arginine-only run leaves them inert (0 on every ledger), byte-for-
    # byte the pre-D-100 core. The consumers are rewired to draw them at the atomic flip (D-100
    # commit 2); until then this is pure inert scaffolding.
    specs += [
        VarSpec(
            "leucine",
            "g/L",
            default=0.0,
            description="assimilable L-leucine (dosed must input, split from amino_acids_gpl by "
            "the must spectrum; autolysis-refilled). Ehrlich precursor of isoamyl_alcohol (D-99) "
            "and the thermal Strecker precursor of 3-methylbutanal (D-87). Carbon- AND "
            "nitrogen-bearing (deaminated to N when catabolised). Decision D-100",
        ),
        VarSpec(
            "isoleucine",
            "g/L",
            default=0.0,
            description="assimilable L-isoleucine (dosed must input + autolysis-refilled). The "
            "Ehrlich precursor of active_amyl_alcohol (D-99) and the thermal Strecker precursor of "
            "2-methylbutanal (D-87). An isomer of leucine (same C/N weights). Decision D-100",
        ),
        VarSpec(
            "valine",
            "g/L",
            default=0.0,
            description="assimilable L-valine (dosed must input + autolysis-refilled). The Ehrlich "
            "precursor of isobutanol (D-99) and the thermal Strecker precursor of 2-methylpropanal "
            "(D-87). Carbon- AND nitrogen-bearing. Decision D-100",
        ),
        VarSpec(
            "threonine",
            "g/L",
            default=0.0,
            description="assimilable L-threonine (dosed must input + autolysis-refilled). The "
            "Ehrlich precursor of propanol (D-99) AND of sotolon (via α-ketobutyrate, D-87) — "
            "the ONE amino acid shared between a fusel and a thermal-aging product, so a real "
            "(not artifact) propanol-vs-sotolon competition survives the split. Decision D-100",
        ),
        VarSpec(
            "phenylalanine",
            "g/L",
            default=0.0,
            description="assimilable L-phenylalanine (dosed must input + autolysis-refilled). The "
            "Ehrlich precursor of 2_phenylethanol (D-99) AND the oxidative Strecker precursor of "
            "phenylacetaldehyde (D-75) — shared between a fusel and an oxidative-aging aldehyde. "
            "Carbon- AND nitrogen-bearing. Decision D-100",
        ),
        VarSpec(
            "methionine",
            "g/L",
            default=0.0,
            description="assimilable L-methionine (dosed must input + autolysis-refilled). The "
            "sulfur-bearing precursor of the oxidative Strecker aldehyde methional (D-75) and the "
            "autolytic mercaptan methanethiol (D-45). Its sulfur is untracked (the D-45 idiom); "
            "carbon- AND nitrogen-bearing. Trace in must. Decision D-100",
        ),
        VarSpec(
            "amino_acids_generic",
            "g/L",
            default=0.0,
            description="the generic assimilable amino-acid bucket (glutamine stand-in) — every "
            "assimilable amino acid without its own slot (glutamic acid/alanine/serine/aspartate/"
            "histidine/lysine/glycine/tryptophan/tyrosine/cysteine/GABA). Drawn by the identity-"
            "agnostic yeast swap (D-32), MLF growth (D-38) and Brett growth (D-40) ALONGSIDE "
            "arginine so it is not the sole generic source (D-100). Carbon- AND nitrogen-bearing "
            "(glutamine, N-rich). Dosed must input + autolysis-refilled. Decision D-100",
        ),
        # The DMS axis (decision D-102), appended LAST so existing wine slot indices are unchanged
        # (the D-100 convention). Both slots are OFF every ledger — the pair's carbon comes from
        # untracked SMM and lands in untracked DMS, so SMMHydrolysis moves nothing conserved (the
        # D-74 A420 argument), and at µg/L it is ~1e-6 of the carbon ledger regardless.
        VarSpec(
            "dms_potential",
            "g/L",
            default=0.0,
            description="DMS potential (DMSp) — the grape-borne precursor pool that hydrolyses to "
            "dimethyl sulfide during bottle aging (decision D-102), chiefly S-methylmethionine "
            "(SMM). Booked in DMS-EQUIVALENTS (g of the DMS it can release, NOT g of SMM) — the "
            "unit the wine literature reports DMSp in, which makes the SMMHydrolysis conversion "
            "1:1 by construction and sidesteps SMM's molar-mass/iodide-salt-form ambiguity. A "
            "GRAPE property, not a winemaking dose, so _wine_initial seeds it from the sourced "
            "dms_potential_initial rather than 0 (a 0 default would assert aged wine makes no DMS "
            "— the D-45 hard-zero defect). Carbon-negligible at µg/L and on NO ledger (the "
            "h2s/o2/A420 pattern). Decision D-102",
        ),
        VarSpec(
            "dms",
            "g/L",
            default=0.0,
            description="dimethyl sulfide (DMS) — the aged-wine 'truffle / black olive / cooked "
            "corn' odorant (decision D-102), accumulated by SMMHydrolysis as a first-order "
            "Arrhenius decay of dms_potential. Cumulative and monotonic (d(dms)/dt >= 0): unlike "
            "the D-42 h2s residual there is no CO2 stream to strip it (aging is post-dryness), so "
            "it accumulates rather than settling to a residual. Off every ledger (with "
            "dms_potential — see that slot). Default 0 ⇒ no DMS at pitch. Decision D-102",
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
    # Oak extraction axis — barrel-beer oak (decision D-86). The SAME wood-extractive axis wine
    # carries (D-77 aroma four + D-78 ellagitannin), appended here so barrel/foeder-aged beer
    # (bourbon-barrel stouts, oak-aged sours) extracts oak aroma + tannin. Off every ledger, so an
    # un-oaked beer is byte-for-byte unchanged (every ceiling 0 ⇒ the oak Processes are inert).
    # Extraction is a WOOD property (medium-agnostic oak.yaml yields); only the OAV perception
    # thresholds are matrix-specific (threshold_<compound>_beer, sensory.yaml).
    specs += _oak_specs()
    # Caramelization melanoidin carbon-park (decision D-88; medium-agnostic D-90). The SAME brown
    # thermal-browning polymer Caramelization forms in wine — beer's residual dextrins (unfermented
    # maltose/maltotriose) caramelize by HEAT with NO O₂, browning an aged/warm-stored beer and
    # raising the SAME A420 index D-74/D-86 accumulate. Appended here (not in _common_specs, so
    # wine_schema keeps its single append) — a carbon-park ON total_carbon: the consumed core-S
    # carbon (all three beer sugars, apportioned at each component's own carbon fraction) redeposits
    # here, so sugar → melanoidin closes exactly. Sugar-only (nitrogen-free — caramelization, not
    # Maillard; the N-incorporating MaillardBrowning D-89 stays wine-only, since beer's amino-acid
    # pool is untracked, D-32). Default 0 ⇒ a dry-finished beer (S ≈ 0 at begin_aging) browns
    # negligibly; a high-residual/under-attenuated beer browns meaningfully.
    specs.append(
        VarSpec(
            "melanoidin",
            "g/L",
            default=0.0,
            description="melanoidin — the brown caramelization polymer (caramelan stand-in "
            "C12H18O9) Caramelization forms from residual sugar by HEAT with NO O₂ (decision D-88, "
            "the O₂-independent thermal mirror of PhenolicBrowning D-74; medium-agnostic D-90 — "
            "beer's residual dextrins caramelize too). A carbon-park pool (the debris/glucan "
            "precedent): ON total_carbon (it holds the consumed core-S carbon — all three beer "
            "sugars at their own carbon fractions — so sugar → melanoidin closes exactly), unlike "
            "the off-ledger oak lumps. Sugar-only (nitrogen-free — caramelization, not amino-acid "
            "Maillard). Raises the shared A420 browning index (read by analysis.a420); the "
            "melanoidin MASS itself is not a sensory pool. In BOTH media (D-90)",
        )
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
#: and ``total_carbon`` weights the ester/``fusels`` pools, so they no longer
#: leave the core byte-for-byte when enabled — turning them on draws a *trace* of
#: sugar (~0.2 % of ``S0``), perturbing only ``dS`` (never ``dE``/``dCO2``). Carbon
#: still closes to machine precision with them on, and the §2.2 trio stays in band.
#: See D-19 / milestone-2-tasks.md.
#:
#: :class:`EsterVolatilization` (decision D-20) is the gas-stripping sink that moves
#: each liquid ester into its own bookkeeping headspace twin as CO2 sparges
#: the must — the physics behind wine's "warmer ⇒ *less* liquid ester" (Rollero 2014):
#: with ``E_a_ester_volatil`` set *per medium* it is held **above** ``E_a_esters`` for
#: wine (stripping outruns synthesis, liquid esters fall with T) and **below** it for
#: beer (synthesis dominates, esters rise with T — de Andrés-Toro). The transfer is
#: carbon-neutral (each pool and its twin book as the SAME one of the three D-96 ester
#: molecules), so it is in this isolable tuple too and ``total_carbon`` still closes to
#: machine precision.
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
#: D-68..D-74). Three medium-agnostic Processes: :class:`EsterHydrolysis` (D-69, the first §4.1
#: Process — young fruity acetate esters hydrolyse back toward equilibrium with age, releasing
#: carbon 5:2 into ``fusels`` + ``Byp``), :class:`OxidativeAcetaldehyde` (D-71, the first OXIDATIVE
#: Process — dissolved O₂ oxidises ethanol → acetaldehyde, the 'sherry'/oxidised note, saturating as
#: the ``o2`` charge is spent) and :class:`PhenolicBrowning` (D-74, the second always-on O₂ sink —
#: O₂ oxidises phenolics to brown pigment, accumulating the ``A420`` browning index; the dominant
#: O₂ consumer, so it diverts O₂ from — and suppresses — oxidative acetaldehyde). ALL
#: MEDIUM-AGNOSTIC
#: — hydrolysis and oxidation are properties of the molecules and the wine/beer pH, not the biology
#: (the ``vicinal_diketones.yaml`` / shared-file pattern); the ester/``fusels``/``Byp``/
#: ``acetaldehyde``/``o2``/``A420`` exist in both schemas, and both wine and beer carry autoxidising
#: polyphenols that consume O₂ and brown (D-74) — so all three are wired into BOTH media. Kept in
#: their OWN isolable tuple (prime directive #3): a
#: ProcessSet built without it is the pre-aging model. Unlike the always-on intrinsic aroma pools,
#: aging is INHERENTLY post-ferment (there is no aging at t0), so the compile seam DISABLES the
#: whole tuple unconditionally and a ``begin_aging`` intervention (decision D-70, the ``pitch_mlf``
#: reconfigure pattern MINUS the state mutation) re-enables it over a post-fermentation aging
#: segment — off during the ferment, on during aging. An un-aged run is thus byte-for-byte the
#: pre-aging core (disabled ⇒ skipped by ``active`` / ``tier_of`` / the strict ``touches`` check).
#: During a post-dryness aging segment every OTHER producer of the ester/``fusels``/``Byp`` /
#: ``acetaldehyde`` (``ester_synthesis``, ``ester_volatilization``, ``fusel_alcohols_ehrlich``, the
#: ``Byp`` uptake routing, and ``acetaldehyde_production``/``_reduction``) is fermentative-flux- or
#: viable-``X``-gated and quiescent at ``S ≈ 0`` / ``X = 0``, so the aging signal is UNCONFOUNDED —
#: only the aging Processes move those pools (Stance A, D-70). :class:`OxidativeAcetaldehyde` and
#: :class:`PhenolicBrowning` add a further gate: both are inert unless O₂ is dosed (``add_oxygen``),
#: so a ``begin_aging`` run with no oxygen is purely *reductive* aging — byte-for-byte the
#: ester-hydrolysis-only case (D-71/D-74; ``o2 ≤ 0`` ⇒ both contribute zero, ``A420`` stays 0).
#: Params live in the shared, medium-agnostic ``aging.yaml``.
_AGING_PROCESSES: tuple[Callable[[], Process], ...] = (
    EsterHydrolysis,
    EthylHexanoateHydrolysis,
    EthylAcetateEsterification,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
)

#: WINE-ONLY oxidative-aging Processes that draw on wine-only state (decision D-72). Unlike the
#: medium-agnostic ``_AGING_PROCESSES`` above, :class:`SulfiteOxidation` reads ``so2_total`` and the
#: acid/cation pH slots — all wine-only (beer's pH/SO₂ system is deferred, D-18) — so it is wired
#: into the *wine* medium only, exactly like ``_MLF_PROCESSES``/``_BRETT_PROCESSES``. It is the
#: first sink to claim its share of the shared ``o2`` budget opened by ``OxidativeAcetaldehyde``
#: (D-71): dissolved O₂ oxidises free bisulfite → sulfate, so while free SO₂ lasts O₂ is diverted
#: from ethanol oxidation and oxidative acetaldehyde is suppressed — the "SO₂ protects until
#: exhausted" threshold, emergent from the two Processes summing over ``o2``. Like the rest of the
#: aging axis it is DISABLED at the compile seam and re-enabled by ``begin_aging`` (its name rides
#: in :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in ``aging.yaml``.
_OXIDATIVE_SO2_PROCESSES: tuple[Callable[[], Process], ...] = (SulfiteOxidation,)

#: WINE-ONLY Strecker-degradation aging Process (decision D-75). Like ``_OXIDATIVE_SO2_PROCESSES``,
#: :class:`StreckerDegradation` reads wine-only state (``amino_acids`` + the ``N``-deamination),
#: so it is wired into the *wine* medium only. It is the third oxidative sibling on shared ``o2``
#: budget (after ``OxidativeAcetaldehyde``/``PhenolicBrowning``): dissolved O₂ — via the phenol-
#: oxidation quinones — degrades amino acids to the Strecker aldehydes ``methional`` (cooked-potato)
#: and ``phenylacetaldehyde`` (honey), drawing carbon from ``amino_acids`` and deaminating the
#: nitrogen to ``N`` (the D-45 mercaptan idiom + a CO₂ decarboxylation term). DOUBLY substrate-gated
#: (on ``o2`` AND ``amino_acids``), so — like ``SulfiteOxidation`` — it adds on top of the O₂ budget
#: WITHOUT re-baselining the anchor (superseding the D-71→D-74 forward-guess; see the Process
#: docstring and D-75). Kept in its OWN tuple (isolable, directive #3): DISABLED at the compile
#: seam and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in ``aging.yaml``.
_STRECKER_PROCESSES: tuple[Callable[[], Process], ...] = (StreckerDegradation,)

#: WINE-ONLY non-oxidative THERMAL Strecker aging Process (decision D-87) — the O₂-INDEPENDENT
#: thermal mirror of ``_STRECKER_PROCESSES``. :class:`MaillardStrecker` reads wine-only state
#: (``amino_acids`` + the ``N``-deamination) and the residual-sugar driver, so it is wired into the
#: *wine* medium only. Unlike the D-75 oxidative route it draws NO ``o2``: residual sugar forms
#: α-dicarbonyls (Maillard) that deaminate + decarboxylate amino acids to the sweet-wine / Madeira
#: aldehyde suite — methional + phenylacetaldehyde (shared with D-75), the three branched-chain
#: malty aldehydes (``2_methylbutanal`` / ``3_methylbutanal`` / ``2_methylpropanal``) and
#: ``sotolon`` (the curry/maple furanone) — with NO oxygen, so a sealed sweet wine ages thermally.
#: DOUBLY substrate-gated (on residual sugar AND ``amino_acids``): like ``_STRECKER_PROCESSES`` it
#: adds on top of the aging trajectory WITHOUT re-baselining, and — sharing the ``amino_acids``
#: limiting reagent with the D-75 route via ``ProcessSet`` summing — the two Strecker routes are
#: additive over that pool (the o2-sharing pattern applied to amino_acids). Kept in its OWN tuple
#: (isolable, directive #3): DISABLED at the compile seam and re-enabled by ``begin_aging`` (its
#: name
#: rides in :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in
#: ``thermal.yaml``.
_MAILLARD_STRECKER_PROCESSES: tuple[Callable[[], Process], ...] = (
    MaillardStrecker,
    SotolonAldolCondensation,
)

#: MEDIUM-AGNOSTIC (WINE + BEER, decision D-88; extended to beer D-90) non-oxidative THERMAL
#: browning Process — the O₂-INDEPENDENT thermal mirror of :class:`PhenolicBrowning` (D-74) and the
#: browning half of the thermal axis :class:`MaillardStrecker` (D-87) opened.
#: :class:`Caramelization`
#: browns **residual sugar** to ``melanoidin`` by heat alone (no ``o2``), raising the SAME ``A420``
#: index D-74/D-86 accumulate — so a sealed sweet wine *or* high-residual beer still darkens with
#: age. It is the **first aging Process to consume core** ``S``, so its carbon lands in the
#: on-ledger
#: ``melanoidin`` carbon-park (``total_carbon`` closes); the D-90 vectorized draw apportions the
#: debit across beer's 3-slot ``S`` at each component's own carbon fraction, so both media close.
#: SUGAR-ONLY (nitrogen-free — caramelization, not Maillard; the N-incorporating
#: :class:`MaillardBrowning` D-89 stays wine-only since beer's ``amino_acids`` are untracked, D-32).
#: Wired into BOTH media (D-90: the ``melanoidin`` slot is appended to both ``wine_schema`` and
#: ``beer_schema``, the D-86 oak-to-beer pattern), unlike the wine-only
#: ``_MAILLARD_STRECKER_PROCESSES`` /
#: ``_MAILLARD_BROWNING_PROCESSES``. Kept in its OWN tuple (isolable, directive #3): DISABLED at the
#: compile seam and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). With ``S ≈ 0`` at the aging
#: segment (a standard dry aging run) it is byte-for-byte (wine) / numerically (beer) inert (the
#: per-component clamp + ``S ≤ 0`` guard). Params live in ``thermal.yaml``.
_CARAMELIZATION_PROCESSES: tuple[Callable[[], Process], ...] = (Caramelization,)

#: WINE-ONLY non-oxidative amino-acid-incorporating THERMAL browning Process (decision D-89) — the
#: N-bearing browning branch D-88's sugar-only :class:`Caramelization` deferred. :class:`\
#: MaillardBrowning` browns **residual sugar + amino acids** to a nitrogen-bearing
#: ``maillard_melanoidin`` polymer by heat alone (no ``o2``), raising the SAME ``A420`` index
#: D-74/D-88 accumulate. It consumes core ``S`` **and** ``amino_acids`` and RETAINS the amino-acid
#: nitrogen in the polymer, so ``maillard_melanoidin`` is the FIRST non-biomass, non-arginine
#: species on ``total_nitrogen`` (the draws are sized to it, so both carbon and nitrogen close
#: exactly). Draws the shared ``amino_acids`` reagent with :class:`MaillardStrecker` (D-87) and the
#: shared ``S`` with :class:`Caramelization` (D-88) — ``ProcessSet`` sums the three thermal branches
#: over those reagents (the o2-sharing pattern, no double-count). Wine-only for v1 (beer thermal
#: browning deferred, the D-86 oak-to-beer pattern). Kept in its OWN tuple (isolable, directive #3):
#: DISABLED at the compile seam and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Isolability on the
#: ``amino_acids`` HARD gate (undosed ⇒ byte-for-byte inert). Params live in ``thermal.yaml``.
_MAILLARD_BROWNING_PROCESSES: tuple[Callable[[], Process], ...] = (MaillardBrowning,)

#: DMS-via-SMM-hydrolysis aging Process (decision D-102) — the aged-wine truffle/black-olive
#: odorant. WINE-ONLY. :class:`SMMHydrolysis` is a **distinct route**, and that is the point: every
#: other sulfur pool here is autolysis-gated (D-44 ``h2s``, D-45 ``mercaptans``), whereas DMS
#: accumulates by spontaneous hydrolysis of the grape-borne precursor during bottle aging — lees or
#: no lees — so it carries its **own anchor** instead of ratio-splitting a shared autolytic yield
#: (the D-96 linchpin ``mercaptans`` could not satisfy, D-101). Booked in **DMS-equivalents**, so
#: ``dms_potential`` → ``dms`` is 1:1 with no yield parameter and no molar-mass conversion. OFF
#: EVERY LEDGER (both slots — the D-74 ``A420`` argument: untracked precursor → untracked product,
#: so it moves nothing conserved; at µg/L it is ~1e-6 of the carbon ledger regardless), so it needs
#: no ``chemistry.py`` species registration. Wine-only because the CONSTANTS are wine-anchored:
#: beer's DMS is real and better-studied, but arrives by other routes entirely (wort-boil SMM
#: cleavage *before* pitch; yeast DMSO reduction during ferment), and transferring these constants
#: to it would be the exact wort→wine mechanism error D-102 rejects Scheuren for. Kept in its OWN
#: tuple (isolable, directive #3): DISABLED at the compile seam and re-enabled by ``begin_aging``
#: (its name rides in :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). With no
#: precursor seeded ⇒ byte-for-byte inert (the ``dms_potential <= 0`` guard). Params: ``dms.yaml``.
_DMS_PROCESSES: tuple[Callable[[], Process], ...] = (SMMHydrolysis,)

#: Oak-extraction aging Process (decision D-77) — the barrel/chip extractive axis. WINE + BARREL-
#: BEER (D-86: wired into BOTH media — the oak axis is a wood property, not a grape one).
#: :class:`OakExtraction` is the first **non-oxidative** aging Process: it draws NO O₂, so it takes
#: no share of the shared ``o2`` budget. As a finished beverage sits in oak, four AROMA extractives
#: — ``whiskey_lactone`` (coconut), ``vanillin`` (vanilla), ``guaiacol`` (smoky) and ``eugenol``
#: (clove) — PLUS the ``ellagitannin`` TASTE extractive (D-78) diffuse in and rise toward a
#: per-compound saturation ceiling (first-order approach from below, the inverse of
#: :class:`EsterHydrolysis`). The aroma four are a **separate axis**, O₂-orthogonal;
#: ``ellagitannin`` bridges to the O₂ sub-axis (see :data:`_ELLAGITANNIN_PROCESSES`) but its
#: *extraction* is the same pure diffusion this Process performs. The ceilings are SET-AND-HOLD
#: state slots the ``add_oak`` verb writes (``oak_gpl`` × toast-specific yield); this
#: Process reads them and rises the extracted pools toward them. Wired into BOTH media (D-86: the
#: oak slots are appended to both ``wine_schema`` and ``beer_schema`` via :func:`_oak_specs`), so —
#: unlike ``_OXIDATIVE_SO2_PROCESSES``/``_STRECKER_PROCESSES`` (wine-only, grape/pH-coupled) — it
#: is medium-agnostic (like ``_AGING_PROCESSES``). OFF EVERY LEDGER (exogenous wood-derived
#: mass, the ``iso_alpha`` precedent), so it moves nothing conserved and — a pure g/L transfer —
#: needs no ``chemistry.py`` species registration. Kept in its OWN tuple (isolable, directive #3):
#: DISABLED at the compile seam and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). With no oak dosed every ceiling
#: is 0 ⇒ byte-for-byte inert (the ceiling ≤ 0 guard). Params live in ``oak.yaml``.
_OAK_PROCESSES: tuple[Callable[[], Process], ...] = (OakExtraction,)

#: Ellagitannin O₂-scavenging aging Process (decision D-78) — the BRIDGE from the oak extractive
#: axis to the O₂ sub-axis. WINE + BARREL-BEER (D-86: wired into BOTH media alongside
#: ``_OAK_PROCESSES`` — the ``ellagitannin`` slots and the ``o2`` pool are both medium-agnostic).
#: :class:`EllagitanninOxidation` is the fourth oxidative
#: sibling to claim a share of the shared ``o2`` budget (after
#: ``OxidativeAcetaldehyde``/``PhenolicBrowning``/ ``SulfiteOxidation``): oak's hydrolysable tannin
#: (the ``ellagitannin`` pool that ``OakExtraction`` fills) is a sacrificial antioxidant — dissolved
#: O₂ oxidises it (bilinear ``[o2]·[ellagitannin]``, the :class:`SulfiteOxidation` form), CONSUMING
#: the tannin as it scavenges. So an oaked + oxygenated beverage browns LESS and accumulates LESS
#: oxidative acetaldehyde than an un-oaked one at the same O₂ dose — the oak-PROTECTION emergent
#: (the D-78 spine, the D-72 "SO₂ protects" threshold with a *renewable* buffer: the wood
#: re-supplies tannin below the ceiling). SUBSTRATE-GATED on the ``ellagitannin`` pool ⇒ zero unless
#: oak is dosed ⇒ adds on top of the anchor with NO re-baseline (the D-72/D-75 rule;
#: ``k_ethanol_oxidation + k_browning = 5.0e-4`` untouched). Wired into BOTH media (D-86), unlike
#: wine-only ``_STRECKER_PROCESSES``. OFF EVERY LEDGER (both ``o2`` and ``ellagitannin``
#: are unweighted), so — like ``SulfiteOxidation`` — it moves nothing conserved. Kept in its OWN
#: tuple (isolable, directive #3): DISABLED at compile and re-enabled by ``begin_aging`` (its name
#: rides in :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). With no oak dosed the
#: ``ellagitannin`` pool is 0 ⇒ byte-for-byte inert. Params live in ``oak.yaml`` (with the tannin's
#: extraction yields — all ellagitannin data together).
_ELLAGITANNIN_PROCESSES: tuple[Callable[[], Process], ...] = (EllagitanninOxidation,)

#: WINE-ONLY tannin–anthocyanin condensation aging Process (decision D-79) — the red-wine
#: colour-stabilization + astringency-softening axis, the DOMINANT softening mechanism D-77/D-78
#: deferred. :class:`TanninAnthocyaninCondensation` is the second **non-oxidative** aging Process
#: (after :class:`OakExtraction`) and a **third separate axis**: as a finished red wine ages, free
#: grape ``anthocyanin`` and condensed ``tannin`` combine (bilinear ``[anthocyanin]·[tannin]``,
#: the :class:`SulfiteOxidation` form) into a stable polymeric pigment — softening the astringency
#: and stabilizing the colour. **OAK- AND O₂-INDEPENDENT** (the D-79 crux): it draws NO share of the
#: shared ``o2`` budget (unlike every D-71..D-78 oxidative sink) and reads NO oak pool — grape
#: condensed ``tannin`` differs from oak hydrolysable ``ellagitannin`` (D-78), so a
#: steel-tank red with no oak and no oxygen still polymerizes (a reused-ellagitannin design would
#: wrongly require ``add_oak``). Wired into the *wine* medium only (both grape slots are wine-only),
#: like ``_OAK_PROCESSES``/``_ELLAGITANNIN_PROCESSES``. OFF EVERY LEDGER (both grape pools are
#: unweighted, the ``iso_alpha``/``ellagitannin`` precedent), so — like :class:`OakExtraction` — it
#: moves nothing conserved. DOUBLY substrate-gated on ``anthocyanin`` AND ``tannin`` ⇒ zero unless
#: BOTH are dosed (a white / no-tannin wine is byte-for-byte inert) ⇒ adds ON TOP, NO re-baseline
#: (and trivially so — no ``o2`` term, so it never touches the ``k_ethanol_oxidation + k_browning``
#: anchor). The polymeric-pigment product is a POST-HOC readout
#: (:func:`~fermentation.analysis.polymeric_pigment_series` = ``anthocyanin₀ − anthocyanin``), NOT a
#: state slot (the A420 discriminator — anthocyanin's single fate makes it reconstructible). Kept in
#: its OWN isolable tuple (directive #3): DISABLED at compile and re-enabled by ``begin_aging`` (its
#: name rides in :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in
#: ``polymerization.yaml``.
_POLYMERIZATION_PROCESSES: tuple[Callable[[], Process], ...] = (TanninAnthocyaninCondensation,)

#: WINE-ONLY acetaldehyde-bridged condensation aging Process (decision D-80) — the SPLIT-LEDGER
#: colour beat D-79 deferred, and the second pigment-formation pathway (after
#: :class:`TanninAnthocyaninCondensation`). :class:`AcetaldehydeBridgedCondensation` is the third
#: **non-oxidative** aging Process and the FIRST aging colour Process on the **carbon ledger**: as a
#: finished red wine takes up O₂ (micro-oxygenation), the dissolved-O₂ acetaldehyde
#: (:class:`OxidativeAcetaldehyde`, D-71) forms an ethylidene bridge ``—CH(CH₃)—`` linking grape
#: ``tannin`` to ``anthocyanin`` (trilinear ``[acetaldehyde]·[anthocyanin]·[tannin]``), stabilizing
#: colour and softening astringency — **the first link from the oxidative sub-axis to red-wine
#: colour** (the "controlled micro-ox stabilizes colour" payoff D-79 named). The grape bulk stays
#: OFF
#: every ledger (the D-79 precedent), but acetaldehyde's carbon is ON the ledger (borrowed from
#: ``E``
#: at D-71), so a new on-ledger ``ethyl_bridge`` slot captures it (weighted at ``cf(ethylidene)`` in
#: ``total_carbon``) — the SPLIT LEDGER that keeps carbon from vanishing (the trap D-79 named).
#: **Reads FREE acetaldehyde** under SO₂ (bound acetaldehyde can't bridge — the D-47 precedent), so
#: SO₂ *delays* colour stabilization (emergent). TRIPLY substrate-gated on ``acetaldehyde`` AND
#: ``anthocyanin`` AND ``tannin`` ⇒ zero unless all present (a white / no-tannin / no-acetaldehyde
#: wine is byte-for-byte inert). Wine-only (the grape/bridge slots are wine-only), like
#: ``_POLYMERIZATION_PROCESSES``. Kept in its OWN isolable tuple (directive #3): DISABLED at compile
#: and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in
#: ``polymerization.yaml`` (with the direct route's — all condensation data together).
_ACETALDEHYDE_BRIDGE_PROCESSES: tuple[Callable[[], Process], ...] = (
    AcetaldehydeBridgedCondensation,
)

#: WINE-ONLY oxidative anthocyanin-fading aging Process (decision D-81) — the O₂-coupled bleaching
#: loss that finally makes :func:`~fermentation.analysis.color_series` genuinely DECLINE.
#: :class:`AnthocyaninFading` is the sixth **oxidative** aging sink on the shared ``o2`` pool (after
#: :class:`OxidativeAcetaldehyde`/:class:`SulfiteOxidation`/:class:`PhenolicBrowning`/\
#: :class:`StreckerDegradation`/:class:`EllagitanninOxidation`): dissolved O₂ degrades free
#: ``anthocyanin`` to the colourless ``faded_anthocyanin`` slot (bilinear ``[o2]·[anthocyanin]``,
#: the :class:`EllagitanninOxidation` form), a pure off-ledger transfer. Because it draws the SHARED
#: o2 budget, **SO₂ protection is emergent** (SO₂ scavenges o2 via D-72, leaving less to fade the
#: colour) — nothing scripted. It is the second ``anthocyanin`` fate that forced promoting the
#: pigment to a slot (the A420 discriminator). Doubly substrate-gated on ``o2`` AND ``anthocyanin``
#: ⇒
#: zero unless a red is dosed AND oxygenated (a white / reductive / all-beer run is byte-for-byte
#: inert) ⇒ adds ON TOP of the oxidative sub-axis. Wine-only (the grape slots are wine-only), like
#: ``_POLYMERIZATION_PROCESSES``. Kept in its OWN isolable tuple (directive #3): DISABLED at compile
#: and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in
#: ``polymerization.yaml`` (with the condensation data — all colour-axis data together).
_ANTHOCYANIN_FADING_PROCESSES: tuple[Callable[[], Process], ...] = (AnthocyaninFading,)

#: WINE-ONLY O₂-INDEPENDENT thermal anthocyanin-fade aging Process (decision D-83) — the second,
#: non-oxidative fate that fades free ``anthocyanin`` to colourless, the pathway D-81 deferred.
#: :class:`ThermalAnthocyaninFade` degrades free monomeric ``anthocyanin`` to the SAME colourless
#: ``faded_anthocyanin`` slot the D-81 oxidative fade fills, but by a **thermal/hydrolytic** route
#: needing **no oxygen** (first-order ``[anthocyanin]``, the :class:`EsterHydrolysis` form, NOT the
#: D-81 bilinear ``[o2]·[anthocyanin]``). Because it touches **no ``o2``**, **SO₂ does NOT protect**
#: it (the mirror of D-81's emergent SO₂ protection): a sealed, sulfited, anaerobic red still fades
#: thermally, and only cold storage (``E_a > 0``) slows it — so a **reductive** red, flat under D-81
#: alone, now genuinely declines (retiring the D-81 "anaerobic sealed red holds its colour" note). A
#: pure OFF-LEDGER transfer (the D-81 colour identity still closes by construction), **no yield**
#: (the rate is already g anthocyanin/L/h). Wine-only (the grape slots are wine-only), like
#: ``_ANTHOCYANIN_FADING_PROCESSES``. Kept in its OWN isolable tuple (directive #3): DISABLED at
#: compile and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in
#: ``polymerization.yaml`` (with the condensation/fade data — all colour-axis data together).
_THERMAL_FADE_PROCESSES: tuple[Callable[[], Process], ...] = (ThermalAnthocyaninFade,)

#: WINE-ONLY tannin self-polymerization aging Process (decision D-84) — the first of the
#: tannin–tannin axis the D-79/D-80 condensation beats deferred. :class:`TanninSelfPolymerization`
#: condenses grape ``tannin`` WITH ITSELF (bimolecular ``[tannin]²``, a true self-reaction) into a
#: softer polymer, drawing the free-tannin pool down as a **pure off-ledger sink** (the soft polymer
#: goes to no slot — the D-79/D-80 tannin-is-a-pure-sink precedent, since no ledger reads tannin
#: mass). So astringency (:func:`~fermentation.analysis.astringency_series`) softens **WITHOUT
#: needing anthocyanin** — a white / tannin-only wine now softens, retiring the D-80
#: "one-directional-per-pool" honesty note. OAK- AND O₂-INDEPENDENT (grape condensed tannin, not oak
#: ``ellagitannin``; not an oxidation, draws no ``o2``) and acetaldehyde-free (the DIRECT route —
#: the bridged tannin–ethyl–tannin is D-85), so a steel-tank red still self-polymerizes. OFF EVERY
#: LEDGER (grape-derived), so it moves nothing
#: conserved. Substrate-gated on ``tannin`` ⇒ a no-tannin run is byte-for-byte inert. Wine-only (the
#: grape slot is wine-only), like ``_POLYMERIZATION_PROCESSES``. Kept in its OWN isolable tuple
#: (directive #3): DISABLED at compile and re-enabled by ``begin_aging`` (its name rides in
#: :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in
#: ``polymerization.yaml`` (with the condensation/fade data — all colour/tannin-axis data together).
_TANNIN_SELF_POLYMERIZATION_PROCESSES: tuple[Callable[[], Process], ...] = (
    TanninSelfPolymerization,
)

#: WINE-ONLY acetaldehyde-bridged tannin–ethyl–tannin aging Process (decision D-85) — the second of
#: the tannin–tannin axis, the acetaldehyde-bridged sibling of :class:`TanninSelfPolymerization`
#: (exactly as :class:`AcetaldehydeBridgedCondensation` D-80 is of
#: :class:`TanninAnthocyaninCondensation` D-79). :class:`TanninEthylTanninCondensation` bridges
#: **two** grape ``tannin`` flavanols with a
#: dissolved-O₂ acetaldehyde ethylidene linker (trilinear ``[acetaldehyde]·[tannin]²``), softening
#: astringency — so micro-oxygenation softens even an **anthocyanin-free** tannin pool. Like D-80 it
#: consumes ON-ledger ``acetaldehyde`` and captures its carbon in the **shared** on-ledger
#: ``ethyl_bridge`` slot (the split-ledger carbon-exact transfer), but with its **own**
#: ``y_acetaldehyde_per_tannin`` (one acetaldehyde per two flavanols) and — unlike D-80 — deposits
#: **no** ``polymeric_pigment`` (a colourless tannin–tannin polymer; the tannin sink goes to no
#: slot,
#: the D-84 precedent). **Reads FREE acetaldehyde** under SO₂ (bound can't bridge — the D-47/D-80
#: precedent), so SO₂ *delays* the softening (emergent). TRIPLY substrate-gated on ``acetaldehyde``
#: AND ``tannin`` ⇒ a no-tannin / no-acetaldehyde wine is byte-for-byte inert. Wine-only (the
#: grape/bridge slots are wine-only), like ``_ACETALDEHYDE_BRIDGE_PROCESSES``. Kept in its OWN
#: isolable tuple (directive #3): DISABLED at compile and re-enabled by ``begin_aging`` (its name
#: rides in :data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`). Params live in
#: ``polymerization.yaml`` (with the condensation/fade data — all colour/tannin-axis data together).
_TANNIN_ETHYL_TANNIN_PROCESSES: tuple[Callable[[], Process], ...] = (TanninEthylTanninCondensation,)

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
    AlphaKetobutyrateExcretion,
    AlphaKetobutyrateReassimilation,
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
#:
#: :class:`PrecursorNonEhrlichFates` (decision D-104) rides here too, and for the same reason the
#: re-route does — it is **not** a modifier target. It scales the re-route's own per-species draw
#: by ``f/(1−f)``, so it inherits that draw's temperature shape exactly; scaling it by the growth
#: modifiers (as the swap is) would break the ratio it exists to impose. It closes the sink D-100
#: documented and left out: before it, the re-route was each precursor's ONLY consumer, so 100% of
#: consumed leucine was attributed to isoamyl alcohol where real yeast send 77–86% of it to
#: protein. Disabled with the swap and the re-route when amino acids are un-dosed.
_AMINO_ACID_PROCESSES: tuple[Callable[[], Process], ...] = (
    AminoAcidAssimilation,
    FuselAminoAcidReroute,
    PrecursorNonEhrlichFates,
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
            + _OXIDATIVE_SO2_PROCESSES
            + _STRECKER_PROCESSES
            + _MAILLARD_STRECKER_PROCESSES
            + _CARAMELIZATION_PROCESSES
            + _MAILLARD_BROWNING_PROCESSES
            + _OAK_PROCESSES
            + _ELLAGITANNIN_PROCESSES
            + _POLYMERIZATION_PROCESSES
            + _ACETALDEHYDE_BRIDGE_PROCESSES
            + _ANTHOCYANIN_FADING_PROCESSES
            + _THERMAL_FADE_PROCESSES
            + _TANNIN_SELF_POLYMERIZATION_PROCESSES
            + _TANNIN_ETHYL_TANNIN_PROCESSES
            + _DMS_PROCESSES
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
            + _CARAMELIZATION_PROCESSES
            + _OAK_PROCESSES
            + _ELLAGITANNIN_PROCESSES
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
