"""Medium state layouts and the medium registry."""

import dataclasses

import numpy as np
import pytest

from fermentation.core.media import MEDIA, Medium, beer_schema, get_medium, wine_schema

SHARED = (
    "X", "S", "E", "N", "T", "CO2", "X_dead", "Gly", "Byp", "esters", "fusels", "esters_gas",
    "acetolactate", "diacetyl", "butanediol", "acetaldehyde", "h2s", "h2s_gas", "o2", "A420",
)  # fmt: skip


#: Wine appends three charge-active acid slots + the citric-acid input + the strong-cation
#: slot to the shared set (the pH charge-balance solver reads all but citrate, D-18; citrate
#: is carbon-active-not-charge-active, D-31), then the free-SO₂ pool for the molecular-SO₂
#: readout (decision D-22), then the malolactic-catalyst slot (decision D-23); beer does
#: not (its acid system, SO₂ and MLF are deferred).
WINE_ACID_SLOTS = ("tartaric", "malic", "lactic", "citrate", "cation_charge")
WINE_SO2_SLOTS = ("so2_total",)
WINE_MLF_SLOTS = ("X_mlf", "X_mlf_dead")
WINE_AMINO_ACID_SLOTS = ("amino_acids",)
# The non-assimilable cell-wall debris pool yeast autolysis fills (D-34).
WINE_DEBRIS_SLOTS = ("debris",)
# Brettanomyces volatile-phenol slots (decision D-40), appended last: the p-coumaric-branch
# precursor/intermediate/readout (hydroxycinnamics/vinylphenols/ethylphenols), the ferulic-branch
# precursor/intermediate/readout split out at decision D-55 (ferulic_acid/vinylguaiacols/
# ethylguaiacols — a genuinely distinct molecule, not a fixed-ratio split of the p-coumaric pool),
# and the viable/dead Brett biomass pools (X_brett a constant catalyst in pt1; X_brett_dead filled
# by BrettDeath, pt3).
WINE_BRETT_SLOTS = (
    "hydroxycinnamics", "vinylphenols", "ethylphenols",
    "ferulic_acid", "vinylguaiacols", "ethylguaiacols",
    "X_brett", "X_brett_dead",
)  # fmt: skip
# The carbon-bearing volatile-thiol pool AutolyticMercaptan fills (decision D-45), appended last.
WINE_MERCAPTAN_SLOTS = ("mercaptans",)
# The excreted keto-acid overflow pools (decisions D-49, D-50), appended last: pyruvate then
# alpha-ketoglutarate, the second- and third-strongest SO₂-binding carbonyls after acetaldehyde.
WINE_KETO_ACID_SLOTS = ("pyruvate", "alpha_ketoglutarate")
# The two Strecker-aldehyde aroma pools (decision D-75), appended last: methional (potato) and
# phenylacetaldehyde (honey), the oxidative-aging markers StreckerDegradation produces from amino
# acids. Wine-only (the Process reads wine-only amino_acids and deaminates to N).
WINE_STRECKER_SLOTS = ("methional", "phenylacetaldehyde")
# The oak-extraction axis (decisions D-77/D-78): 4 extracted aroma pools
# (whiskey_lactone/coconut, vanillin/vanilla, guaiacol/smoky, eugenol/clove) + their 4 SET-AND-HOLD
# ceiling slots, then the ellagitannin TASTE/O₂-scavenging pool + its ceiling (D-78, the bridge to
# the O₂ sub-axis). The add_oak verb writes all 5 ceilings (oak_gpl × toast-specific yield); all off
# every ledger (exogenous wood-derived mass, the iso_alpha precedent). SHARED by wine and barrel-
# beer (decision D-86 — the oak axis is a wood property; both media carry it via _oak_specs).
OAK_SLOTS = (
    "whiskey_lactone", "vanillin", "guaiacol", "eugenol",
    "whiskey_lactone_ceiling", "vanillin_ceiling", "guaiacol_ceiling", "eugenol_ceiling",
    "ellagitannin", "ellagitannin_ceiling",
)  # fmt: skip
# The tannin–anthocyanin condensation axis (decision D-79), appended: the two GRAPE must-input
# pools TanninAnthocyaninCondensation condenses into stable polymeric pigment — anthocyanin (the
# bleachable red pigment) + condensed grape tannin (the harsh young astringency). Both off every
# ledger (grape-derived, the iso_alpha/ellagitannin precedent); grape `tannin` is distinct from oak
# `ellagitannin`. Wine-only. The polymeric-pigment product is a post-hoc readout, not a slot. Then
# the ethyl_bridge slot (decision D-80, appended last): the acetaldehyde-derived ethylidene-bridge
# carbon AcetaldehydeBridgedCondensation captures ON the carbon ledger (the split-ledger accounting)
# so the on-ledger acetaldehyde carbon does not vanish into the off-ledger grape pigment. The FIRST
# aging colour slot on the carbon ledger; filled by the Process (no must input). Wine-only.
# Then the D-81 colour-form pair (appended last): polymeric_pigment — the stable pigment PROMOTED
# from the D-79 post-hoc readout to an integrated slot (filled by both condensation routes) — and
# faded_anthocyanin — the colourless, irreversible oxidative-fade sink (filled by AnthocyaninFading,
# O₂-coupled). Both off every ledger (grape-derived colour-equivalents), both wine-only, both filled
# by their Processes (no must input). Together they let color_series genuinely decline and close the
# identity anthocyanin + polymeric_pigment + faded_anthocyanin ≡ anthocyanin₀.
WINE_POLYMERIZATION_SLOTS = (
    "anthocyanin",
    "tannin",
    "ethyl_bridge",
    "polymeric_pigment",
    "faded_anthocyanin",
)
# The four non-oxidative THERMAL Strecker aldehyde/sotolon aroma pools (decision D-87), appended
# last: the sweet-wine / Madeira suite MaillardStrecker produces from residual sugar + amino acids +
# heat with NO O₂. methional + phenylacetaldehyde (WINE_STRECKER_SLOTS, above) are SHARED with this
# route (same molecules), so only these four are new. Carbon-bearing (booked from amino_acids as
# arginine, deaminated to N), on total_carbon like the D-75 pair. Wine-only.
WINE_MAILLARD_SLOTS = (
    "2_methylbutanal",
    "3_methylbutanal",
    "2_methylpropanal",
    "sotolon",
)
# The caramelization melanoidin carbon-park (decision D-88; D-90 medium-agnostic): brown thermal-
# browning polymer Caramelization forms by consuming residual sugar (the O₂-independent mirror of
# PhenolicBrowning D-74). The FIRST aging pool holding consumed core-S carbon, so — unlike the
# off-ledger oak/colour lumps — it is ON total_carbon (sugar → melanoidin closes exactly).
# Sugar-only
# (nitrogen-free — caramelization, not Maillard). In BOTH media (D-90: beer's residual dextrins
# caramelize too, appended to both schemas). Raises the shared A420 index (D-74).
CARAMELIZATION_SLOTS = ("melanoidin",)
# The N-bearing Maillard melanoidin carbon+nitrogen-park (decision D-89), appended last: the brown
# amino-acid-incorporating thermal-browning polymer MaillardBrowning forms by consuming residual
# sugar AND amino acids (the N-incorporating branch D-88's sugar-only Caramelization deferred). ON
# total_carbon AND total_nitrogen — the FIRST non-biomass, non-arginine species on the nitrogen
# ledger (it RETAINS the amino-acid nitrogen). Wine-only. Raises the shared A420 index (D-74/D-88).
WINE_MAILLARD_BROWNING_SLOTS = ("maillard_melanoidin",)

# Beer appends the iso-alpha-acid (bitterness) slot to the shared set — the boil-derived,
# fermentation-lost hop bitterness (decision D-64). Beer-only, exactly as wine's acid/MLF/Brett
# slots are wine-only; off the carbon ledger (exogenous hop-derived mass).
BEER_HOP_SLOTS = ("iso_alpha",)


def test_wine_schema_has_single_sugar_slot():
    schema = wine_schema()
    assert schema.names == (
        SHARED
        + WINE_ACID_SLOTS
        + WINE_SO2_SLOTS
        + WINE_MLF_SLOTS
        + WINE_AMINO_ACID_SLOTS
        + WINE_DEBRIS_SLOTS
        + WINE_BRETT_SLOTS
        + WINE_MERCAPTAN_SLOTS
        + WINE_KETO_ACID_SLOTS
        + WINE_STRECKER_SLOTS
        + OAK_SLOTS
        + WINE_POLYMERIZATION_SLOTS
        + WINE_MAILLARD_SLOTS
        + CARAMELIZATION_SLOTS
        + WINE_MAILLARD_BROWNING_SLOTS
    )
    assert schema.spec("S").size == 1
    # 20 shared (X, S(1), E, N, T, CO2, X_dead, Gly, Byp, esters, fusels, esters_gas,
    # acetolactate, diacetyl, butanediol — the VDK pathway, D-26 — acetaldehyde, D-27,
    # h2s + h2s_gas, D-29 production / D-42 CO2-stripping sink, o2 — the dissolved-oxygen
    # aging substrate, D-71 — and A420 — the oxidative-browning index, D-74)
    # + 3 wine-only acid slots
    # + citrate (D-31) + cation_charge (D-18) + 1 free-SO₂ slot (D-22) + X_mlf + X_mlf_dead
    # slots (D-23 catalyst / D-39 bacterial lees) + 1 amino_acids slot (D-32) + 1 debris slot
    # (D-34) + 8 Brett slots (hydroxycinnamics, vinylphenols, ethylphenols — the p-coumaric
    # branch, D-40; ferulic_acid, vinylguaiacols, ethylguaiacols — the ferulic branch, D-55;
    # X_brett, X_brett_dead) + 1 mercaptans slot (D-45) + 2 keto-acid slots (pyruvate D-49,
    # alpha_ketoglutarate D-50) + 2 Strecker-aldehyde slots (methional, phenylacetaldehyde — D-75)
    # + 10 oak slots (4 aroma extractives whiskey_lactone/vanillin/guaiacol/eugenol + 4 set-and-hold
    # ceilings — the non-oxidative barrel/chip axis, D-77 — plus the ellagitannin
    # TASTE/O₂-scavenging pool + its ceiling, the bridge to the O₂ sub-axis, D-78)
    # + 2 grape polymerization slots (anthocyanin + condensed tannin — the red-wine colour-
    # stabilization + astringency-softening axis, D-79)
    # + 1 ethyl_bridge slot (the acetaldehyde-bridged / split-ledger colour beat, D-80: the first
    # aging colour slot ON the carbon ledger, capturing the acetaldehyde carbon the bridged route
    # consumes so it does not vanish into the off-ledger grape pigment)
    # + 2 D-81 colour-form slots (polymeric_pigment — the stable pigment PROMOTED from readout to an
    # integrated slot — and faded_anthocyanin — the colourless oxidative-fade sink; both off-ledger,
    # letting color_series genuinely decline)
    # + 4 D-87 non-oxidative THERMAL Strecker aldehyde/sotolon slots (2-/3-methylbutanal,
    # 2-methylpropanal, sotolon — the sweet-wine/Madeira suite MaillardStrecker produces from
    # residual sugar + amino acids + heat with NO O₂; methional + phenylacetaldehyde are shared with
    # the D-75 oxidative route, so only these four are new)
    # + 1 D-88 caramelization melanoidin carbon-park slot (the O₂-independent thermal browning of
    # residual sugar → melanoidin, raising the shared A420; the first aging pool on total_carbon
    # that
    # holds consumed core-S carbon)
    # + 1 D-89 N-bearing Maillard melanoidin carbon+nitrogen-park slot (the amino-acid-incorporating
    # thermal browning of residual sugar + amino acids → maillard_melanoidin, raising the same A420;
    # the FIRST non-biomass, non-arginine species on total_nitrogen — it retains the amino-acid N)
    assert schema.size == 64


def test_beer_schema_has_three_sequential_sugars():
    schema = beer_schema()
    assert schema.names == SHARED + BEER_HOP_SLOTS + OAK_SLOTS + CARAMELIZATION_SLOTS
    s = schema.spec("S")
    assert s.size == 3
    assert s.components == ("glucose", "maltose", "maltotriose")
    # 20 shared (X, S(3), E, N, T, CO2, X_dead, Gly, Byp, esters, fusels, esters_gas (D-20),
    # acetolactate, diacetyl, butanediol (VDK pathway, D-26), acetaldehyde (D-27),
    # h2s (D-29) + h2s_gas (D-42 CO2-stripping sink), o2 (D-71 aging substrate), A420 (D-74
    # oxidative-browning index)) + 1 beer-only iso_alpha bitterness slot (D-64) — S occupies 3
    # slots, so 20 shared names span 22 slots, + iso_alpha = 23
    # + 10 oak-axis slots (4 aroma extractives + 4 ceilings — the barrel/chip axis, D-77 — plus the
    # ellagitannin TASTE/O₂-scavenging pool + its ceiling, D-78) SHARED with wine by barrel-beer oak
    # (decision D-86): the oak axis is a wood property, so both media carry it (via _oak_specs) = 33
    # + 1 caramelization melanoidin carbon-park slot (D-90: sugar-only thermal browning is medium-
    # agnostic — beer's residual dextrins caramelize too, so melanoidin is appended to beer_schema
    # too, ON total_carbon; the N-incorporating maillard_melanoidin stays wine-only, D-32) = 34
    assert schema.size == 34


def test_shared_variable_units_are_canonical():
    # Beer carries the shared variable set plus its one beer-only iso_alpha slot (D-64), so it
    # pins the canonical units of the shared layout (wine adds the D-18 acid/cation slots on top
    # — checked below).
    schema = beer_schema()
    units = {spec.name: spec.unit for spec in schema.specs}
    assert units == {
        "X": "g/L",
        "S": "g/L",
        "E": "g/L",
        "N": "g/L",
        "T": "K",
        "CO2": "g/L",
        "X_dead": "g/L",
        "Gly": "g/L",
        "Byp": "g/L",
        "esters": "g/L",
        "fusels": "g/L",
        "esters_gas": "g/L",
        "acetolactate": "g/L",
        "diacetyl": "g/L",
        "butanediol": "g/L",
        "acetaldehyde": "g/L",
        "h2s": "g/L",
        "h2s_gas": "g/L",
        "o2": "g/L",  # dissolved-oxygen aging substrate (D-71)
        "A420": "AU",  # oxidative-browning index — absorbance at 420 nm, dimensionless (D-74)
        "iso_alpha": "g/L",  # beer-only bitterness slot (D-64)
        # Oak-extraction axis — SHARED with wine by barrel-beer oak (D-86); all off-ledger g/L.
        "whiskey_lactone": "g/L",  # coconut oak-lactone (D-77)
        "vanillin": "g/L",  # vanilla (D-77)
        "guaiacol": "g/L",  # oak smoky/toasty (D-77)
        "eugenol": "g/L",  # clove/spice (D-77)
        "whiskey_lactone_ceiling": "g/L",  # set-and-hold ceiling (D-77)
        "vanillin_ceiling": "g/L",
        "guaiacol_ceiling": "g/L",
        "eugenol_ceiling": "g/L",
        "ellagitannin": "g/L",  # oak hydrolysable tannin — taste + O₂ scavenger (D-78)
        "ellagitannin_ceiling": "g/L",
        "melanoidin": "g/L",  # caramelization carbon-park (D-88; medium-agnostic D-90)
    }


def test_wine_acid_slot_units_are_canonical():
    # The D-18 pH-solver slots: acids in g/L (mass concentration, like every other
    # species), the net strong cation as a charge density in mol/L. The D-22 free-SO₂
    # pool is g/L of SO₂-equivalent (mass concentration).
    units = {spec.name: spec.unit for spec in wine_schema().specs}
    assert units["tartaric"] == "g/L"
    assert units["malic"] == "g/L"
    assert units["lactic"] == "g/L"
    assert units["citrate"] == "g/L"  # citric-acid must input (D-31)
    assert units["cation_charge"] == "mol/L"
    assert units["so2_total"] == "g/L"
    assert units["X_mlf"] == "g/L"
    assert units["X_mlf_dead"] == "g/L"  # bacterial lees (D-39)


def test_produced_only_pools_default_to_zero_when_omitted():
    # X_dead/Gly/Byp/esters/fusels/esters_gas are produced-only pools (VarSpec.default=0),
    # so an initial state may omit them; substrate/condition vars stay required (test_state).
    schema = wine_schema()
    arr = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    assert schema.get(arr, "X_dead") == 0.0
    assert schema.get(arr, "Gly") == 0.0
    assert schema.get(arr, "Byp") == 0.0
    assert schema.get(arr, "esters") == 0.0
    assert schema.get(arr, "fusels") == 0.0
    assert schema.get(arr, "esters_gas") == 0.0
    assert schema.get(arr, "acetolactate") == 0.0
    assert schema.get(arr, "diacetyl") == 0.0
    assert schema.get(arr, "butanediol") == 0.0
    assert schema.get(arr, "acetaldehyde") == 0.0


def test_registry_exposes_wine_and_beer():
    assert set(MEDIA) == {"wine", "beer"}
    assert get_medium("wine").schema.spec("S").size == 1
    assert get_medium("beer").schema.spec("S").size == 3


def test_unknown_medium_raises():
    with pytest.raises(KeyError, match="Unknown medium 'cider'"):
        get_medium("cider")


def test_medium_is_frozen():
    with pytest.raises(dataclasses.FrozenInstanceError):
        get_medium("wine").name = "beer"  # type: ignore[misc]


def test_empty_medium_is_the_no_kinetics_baseline():
    # A Medium with no factories assembles an empty set whose total derivative is
    # identically zero — a valid constant-state baseline. The registered wine/beer
    # media now carry kinetics (see below); this is the property a *bare* Medium
    # still guarantees.
    medium = Medium(name="x", schema=wine_schema())
    pset = medium.build_process_set()
    assert pset.active == ()
    assert pset.active_modifiers == ()
    y = medium.schema.zeros()
    deriv = pset.total_derivatives(0.0, y, {})
    assert np.array_equal(deriv, medium.schema.zeros())


def test_build_process_set_respects_strict_flag():
    pset = Medium(name="x", schema=wine_schema()).build_process_set(strict=True)
    assert pset.strict is True


# -- the registered media now carry the validated-core kinetics + Tier-2 byproducts

# Validated core: growth + fermentative uptake + ethanol-driven cell inactivation,
# with per-rate Arrhenius modifiers. The Luong ethanol wall (ethanol_inhibition) is
# retired from the default media in favour of the cumulative inactivation Process
# (decision D-13).
CORE_PROCESSES = {
    "growth_nitrogen_limited",
    "sugar_uptake_to_ethanol_co2",
    "ethanol_inactivation",
}
# Tier-2 aroma byproducts (Milestone 2, decisions D-18/D-19/D-20): additive aroma
# Processes plus the ester gas-stripping sink (ester_volatilization, D-20: liquid
# esters → the esters_gas headspace pool). Wired in by default but isolable (prime
# directive #3) — disabling them leaves the validated core byte-for-byte.
BYPRODUCT_PROCESSES = {"ester_synthesis", "fusel_alcohols_ehrlich", "ester_volatilization"}
# Vicinal-diketone (VDK / diacetyl) pathway (decision D-26): the three-step sugar →
# α-acetolactate → diacetyl + CO2 → 2,3-butanediol chain. Diacetyl is intrinsic yeast
# metabolism, so — unlike MLF — it is wired into BOTH media by default (isolable but always
# on, like the ester/fusel byproducts).
VDK_PROCESSES = {
    "acetolactate_excretion",
    "acetolactate_decarboxylation",
    "diacetyl_reduction",
}
ACETALDEHYDE_PROCESSES = {"acetaldehyde_production", "acetaldehyde_reduction"}
# Hydrogen-sulfide production + CO₂-stripping (decisions D-29 / D-42): the low-nitrogen sulfidic
# off-aroma. A flux-linked producer with an inverse-N gate (D-29), plus the Henry's-law stripping
# sink that sweeps volatile H₂S into the h2s_gas headspace pool (D-42, the ester D-20/D-21
# precedent). Intrinsic yeast metabolism, so — like the ester/VDK/acetaldehyde pools — both are
# wired into BOTH media by default (isolable but always on).
H2S_PROCESSES = {"hydrogen_sulfide_production", "hydrogen_sulfide_volatilization"}
# Temperature-schedule ramp (decision D-35): drives T along a piecewise-linear schedule.
# Medium-agnostic, so wired into BOTH media, and always enabled (it reads the slope with a
# 0.0 isothermal default, so an un-ramped run is byte-for-byte the pre-ramp core).
TEMPERATURE_PROCESSES = {"temperature_ramp"}
# Malolactic fermentation (decisions D-23, D-31) is wired into the WINE medium only (beer has
# no malic/lactic/citrate slots); enabled in a bare build_process_set and disabled at the
# compile seam when O. oeni is not pitched (so undosed wine runs keep malic/lactic/citrate at
# VALIDATED). D-31 adds the citrate → diacetyl co-metabolism and the bacterial diacetyl reducer.
# D-39 adds bacterial death/decay to this same pitch-gated tuple (moves X_mlf → X_mlf_dead under
# molecular SO₂). D-41 (MLF v2) adds the benign senescence baseline — a slow, always-on-when-pitched
# mortality into the same X_mlf_dead pool. Both are pitch-gated, not amino-acid-gated like growth,
# since bacteria die/age whether or not they were growing.
MLF_PROCESSES = {
    "malolactic_conversion",
    "malolactic_citrate_metabolism",
    "oenococcus_diacetyl_reduction",
    "malolactic_death",
    "malolactic_senescence",
}
# Amino-acid ledger (decision D-32) is wired into the WINE medium only (beer has no amino_acids
# slot); enabled in a bare build_process_set and disabled at the compile seam when amino acids
# are not dosed (so undosed wine runs keep the empty amino_acids slot at VALIDATED and are
# byte-for-byte the core).
# The fusel Ehrlich re-route (D-33) rides in the same dosed, wine-only tuple as the swap: it
# sources a fraction of fusel carbon from amino acids and deaminates, disabled at the compile
# seam when amino acids are undosed.
AMINO_ACID_PROCESSES = {"amino_acid_assimilation", "fusel_amino_acid_reroute"}
# Malolactic growth (decision D-23, the deferred growth beat) is wired into the WINE medium only;
# enabled in a bare build_process_set and disabled at the compile seam unless a scenario BOTH
# co-inoculates O. oeni AND doses amino acids (mlf_pitch_gpl>0 AND amino_acids_gpl>0) — a stricter
# gate than conversion, so it is a SEPARATE tuple from MLF_PROCESSES (avoids dragging amino_acids/S
# tier on pitched-but-not-aa-dosed runs).
MLF_GROWTH_PROCESSES = {"malolactic_growth"}
# Brettanomyces volatile phenols (decision D-40): wine-only, pitch-gated — present in a bare build,
# disabled at the compile seam unless a scenario pitches Brett (brett_pitch_gpl>0 or a pitch_brett
# intervention), so an unpitched wine run is byte-for-byte the validated core (the MLF pattern).
# brett_death (D-40 pt3) and brett_ethanol_toxicity (D-58) are pitch-gated too (Brett dies whether
# or not it grew), so they ride in this same tuple — the MalolacticDeath-in-_MLF_PROCESSES pattern.
BRETT_PROCESSES = {
    "brett_decarboxylation", "brett_vinylphenol_reduction", "brett_death", "brett_ethanol_toxicity",
}  # fmt: skip
# Brett growth (decision D-40 pt2, makes X_brett dynamic) is wired into the WINE medium only;
# enabled in a bare build and disabled at the compile seam unless a scenario BOTH pitches Brett AND
# doses amino acids — a stricter gate than the phenol Processes, so it is a SEPARATE tuple (avoids
# dragging the amino_acids/E tier on pitched-but-not-aa-dosed runs), mirroring MLF_GROWTH_PROCESSES.
BRETT_GROWTH_PROCESSES = {"brett_growth"}
# POF+ yeast decarboxylase (decision D-40 pt4) is wired into the WINE medium only; enabled in a bare
# build and disabled at the compile seam unless a scenario opts in via pof_positive. A SEPARATE
# tuple from the Brett Processes: it is gated on the (binary) POF+ strain trait, WHOLLY INDEPENDENT
# of the Brett pitch — a POF+ ferment need not have Brett, and a POF- wine makes no vinylphenol.
POF_PROCESSES = {"yeast_pof_decarboxylation"}
# Yeast autolysis (D-34) + autolytic H₂S source (D-44): wine-only, opt-in — both present in a bare
# build, disabled together at the compile seam unless a scenario passes autolysis_rate_per_h (the
# carrying-capacity opt-in pattern). AutolyticHydrogenSulfide feeds the shared h2s pool the sulfide
# dead cells release, on the same autolysis flux (decision D-44).
AUTOLYSIS_PROCESSES = {"yeast_autolysis", "autolytic_hydrogen_sulfide", "autolytic_mercaptan"}
# Excreted keto-acid overflow pool (decision D-49): wine-only, always-on intrinsic yeast
# metabolism (like the acetaldehyde/VDK/H₂S pools, not a dosed organism) — but wine-only because
# the SO₂-binding competition it exists for is a wine readout (no §2.2 beer benchmark asserts a
# keto-acid level). Excretion draws pyruvate from sugar; the flux-linked reassimilation returns
# it to ethanol+CO₂ and freezes the finished-wine residual at dryness.
KETO_ACID_PROCESSES = {
    "pyruvate_excretion", "pyruvate_reassimilation",
    "alpha_kg_excretion", "alpha_kg_reassimilation",
}  # fmt: skip
# Hop bittering (decision D-64) is wired into the BEER medium only (wine has no iso_alpha slot).
# The boil isomerization is a compile-seam calc, not a Process; the only dynamic member is the
# fermentation-time iso-alpha loss. It is present in a bare build_process_set and disabled at the
# compile seam when no hops are scheduled (the MLF/Brett isolability pattern).
HOP_PROCESSES = {"iso_alpha_acid_loss"}
# Aging chemistry (Milestone 3 / Tier-3, decisions D-68..D-74): the medium-agnostic §4.1 aging
# Processes, ester_hydrolysis (D-69), oxidative_acetaldehyde (D-71, the O₂-driven ethanol oxidation)
# and phenolic_browning (D-74, the O₂-driven browning that accumulates the A420 index). All
# medium-agnostic (hydrolysis/oxidation are molecule/pH properties, and esters/fusels/Byp/
# acetaldehyde/o2/A420 exist in both schemas; both media carry autoxidising polyphenols), so present
# in a bare build_process_set for BOTH media — but DISABLED unconditionally at the compile seam
# (aging is inherently post-ferment, no aging at t0), re-enabled only by a begin_aging intervention.
# An un-aged run is byte-for-byte the pre-aging core (the MLF/Brett isolability pattern, but with no
# t0 co-inoculation path).
AGING_PROCESSES = {"ester_hydrolysis", "oxidative_acetaldehyde", "phenolic_browning"}
# WINE-ONLY aging: sulfite_oxidation (D-72, the O₂-driven SO₂ scavenging) reads wine-only so2_total
# acid-pH slots (beer's pH/SO₂ system is deferred, D-18); strecker_degradation (D-75, the O₂/amino-
# acid-driven Strecker aldehydes methional + phenylacetaldehyde) reads wine-only amino_acids and
# deaminates to N — so — like the MLF/Brett Processes — both are wired into the wine medium only.
# Same compile-seam disable / begin_aging re-enable as the rest.
WINE_AGING_PROCESSES = {"sulfite_oxidation", "strecker_degradation"}
# NON-oxidative oak aging, SHARED by wine and barrel-beer (decision D-86): oak_extraction (D-77, the
# barrel/chip aroma-extractive axis) reads the oak ceiling/extractive slots and draws no O₂ — a
# SEPARATE axis from the oxidative siblings. ellagitannin_oxidation (D-78) is the BRIDGE: it draws
# the O₂-scavenging share of the shared o2 budget as the tannin OakExtraction fills is oxidised (oak
# protection), reading the ellagitannin pool. Both wired into BOTH media (the oak axis is a wood
# property; both schemas carry the slots via _oak_specs), same compile-seam disable / begin_aging
# re-enable.
OAK_PROCESSES = {"oak_extraction", "ellagitannin_oxidation"}
# WINE-ONLY, NON-oxidative aging (D-79): tannin_anthocyanin_condensation condenses the two GRAPE
# pools (anthocyanin + condensed tannin) into stable polymeric pigment — the red-wine colour-
# stabilization + astringency-softening axis. A THIRD separate axis: it draws no O₂ (unlike every
# oxidative sink) and reads no oak pool (grape tannin ≠ oak ellagitannin), so a steel-tank red still
# polymerizes. Wired into the wine medium only, same compile-seam disable / begin_aging re-enable.
WINE_POLYMERIZATION_PROCESSES = {"tannin_anthocyanin_condensation"}
# WINE-ONLY, NON-oxidative aging (D-80): acetaldehyde_bridged_condensation — the SPLIT-LEDGER beat,
# the second pigment-formation pathway and the first link from the oxidative sub-axis to red-wine
# colour. Dissolved-O₂ acetaldehyde bridges grape tannin to anthocyanin (trilinear), capturing the
# on-ledger acetaldehyde carbon in the ethyl_bridge slot. Wine-only, same compile-seam disable /
# begin_aging re-enable.
WINE_BRIDGE_PROCESSES = {"acetaldehyde_bridged_condensation"}
# WINE-ONLY, OXIDATIVE aging (D-81): anthocyanin_fading — the O₂-coupled bleaching loss that makes
# color_series genuinely decline. Dissolved O₂ fades free grape anthocyanin to the colourless
# faded_anthocyanin slot (bilinear [o2]·[anthocyanin], the ellagitannin_oxidation form), drawing the
# shared o2 pool so SO₂ protection is emergent. Wine-only, same compile-seam disable / begin_aging
# re-enable.
WINE_FADING_PROCESSES = {"anthocyanin_fading"}
# WINE-ONLY O₂-INDEPENDENT thermal anthocyanin-fade Process (decision D-83) — the second, non-
# oxidative fate that fades free anthocyanin to the SAME colourless faded_anthocyanin slot, but by a
# thermal/hydrolytic route needing NO oxygen (first-order [anthocyanin]). Touching no o2, SO₂ does
# protect it (the mirror of D-81), so a sealed/sulfited/anaerobic red still fades. Wine-only, same
# compile-seam disable / begin_aging re-enable.
WINE_THERMAL_FADE_PROCESSES = {"thermal_anthocyanin_fade"}
# WINE-ONLY tannin self-polymerization Process (decision D-84) — the first of the tannin–tannin axis
# the D-79/D-80 condensation beats deferred. Condenses grape tannin WITH ITSELF (bimolecular
# [tannin]²) into a soft polymer, a pure off-ledger tannin sink, so astringency softens WITHOUT
# needing anthocyanin (retiring the "one-directional-per-pool" note). Wine-only, same compile-seam
# disable / begin_aging re-enable.
WINE_TANNIN_SELF_POLY_PROCESSES = {"tannin_self_polymerization"}
# WINE-ONLY acetaldehyde-bridged tannin–ethyl–tannin Process (decision D-85) — the second of the
# tannin–tannin axis. Bridges TWO grape tannin flavanols with a dissolved-O₂ acetaldehyde ethylidene
# linker (trilinear [acetaldehyde]·[tannin]²), an O₂-driven softener that reuses the D-80
# split-ledger carbon capture (acetaldehyde → shared ethyl_bridge) but deposits NO pigment
# (colourless). Wine-only,
# same compile-seam disable / begin_aging re-enable.
WINE_TANNIN_ETHYL_TANNIN_PROCESSES = {"tannin_ethyl_tannin_condensation"}
# WINE-ONLY, NON-oxidative THERMAL aging (decision D-87): maillard_strecker — the O₂-independent
# thermal mirror of strecker_degradation. Residual sugar + heat (α-dicarbonyls, no O₂) degrade amino
# acids to the sweet-wine/Madeira aldehyde suite; wine-only (reads amino_acids, deaminates to N),
# additive with the D-75 route over the shared amino_acids limiting reagent. Same compile-seam
# disable / begin_aging re-enable.
WINE_MAILLARD_PROCESSES = {"maillard_strecker"}
# MEDIUM-AGNOSTIC (WINE + BEER), NON-oxidative THERMAL browning (decision D-88; extended to beer
# D-90): caramelization — the O₂-independent thermal mirror of phenolic_browning (D-74). Residual
# sugar browns to the on-ledger melanoidin carbon-park by heat with no O₂, raising the shared A420
# index; the first aging Process to consume core S. In BOTH media (D-90: beer's residual dextrins
# caramelize too, the vectorized draw apportions across beer 3-slot S). Same compile-seam disable/
# begin_aging re-enable.
CARAMELIZATION_PROCESSES = {"caramelization"}
# WINE-ONLY, NON-oxidative amino-acid-incorporating THERMAL browning (decision D-89):
# maillard_browning — the N-bearing browning branch D-88's sugar-only caramelization deferred.
# Residual sugar + amino acids brown to the on-ledger N-bearing maillard_melanoidin park by heat
# with
# no O₂, raising the same A420; the first aging Process on the nitrogen ledger. Draws the shared
# amino_acids (with D-87) and shared S (with D-88). Wine-only v1. Same disable / begin_aging
# re-enable.
WINE_MAILLARD_BROWNING_PROCESSES = {"maillard_browning"}
EXPECTED_PROCESSES = {
    "wine": (
        CORE_PROCESSES
        | TEMPERATURE_PROCESSES
        | BYPRODUCT_PROCESSES
        | VDK_PROCESSES
        | ACETALDEHYDE_PROCESSES
        | H2S_PROCESSES
        | MLF_PROCESSES
        | MLF_GROWTH_PROCESSES
        | BRETT_PROCESSES
        | BRETT_GROWTH_PROCESSES
        | POF_PROCESSES
        | AMINO_ACID_PROCESSES
        | AUTOLYSIS_PROCESSES
        | KETO_ACID_PROCESSES
        | AGING_PROCESSES
        | WINE_AGING_PROCESSES
        | WINE_MAILLARD_PROCESSES
        | CARAMELIZATION_PROCESSES
        | WINE_MAILLARD_BROWNING_PROCESSES
        | OAK_PROCESSES
        | WINE_POLYMERIZATION_PROCESSES
        | WINE_BRIDGE_PROCESSES
        | WINE_FADING_PROCESSES
        | WINE_THERMAL_FADE_PROCESSES
        | WINE_TANNIN_SELF_POLY_PROCESSES
        | WINE_TANNIN_ETHYL_TANNIN_PROCESSES
    ),
    "beer": (
        CORE_PROCESSES
        | TEMPERATURE_PROCESSES
        | BYPRODUCT_PROCESSES
        | VDK_PROCESSES
        | ACETALDEHYDE_PROCESSES
        | H2S_PROCESSES
        | HOP_PROCESSES
        | AGING_PROCESSES
        | CARAMELIZATION_PROCESSES
        | OAK_PROCESSES
    ),
}
# Per-rate Arrhenius modifiers wire into both media. Wine additionally carries the opt-in
# biomass carrying-capacity cap (decision D-30) — enabled in a bare build, but disabled at the
# compile seam unless a scenario passes carrying_capacity_gpl, so an undosed wine run is
# byte-for-byte the validated core (the wine-only MLF *modifier* parallel).
CARRYING_CAPACITY_MODIFIER = "biomass_carrying_capacity"
EXPECTED_MODIFIERS = {
    "wine": {
        "arrhenius_growth",
        "arrhenius_uptake",
        "coleman_death_temperature",
        CARRYING_CAPACITY_MODIFIER,
    },
    "beer": {"arrhenius_growth", "arrhenius_uptake", "coleman_death_temperature"},
}


@pytest.mark.parametrize("medium", ["wine", "beer"])
def test_registered_media_wire_the_full_kinetic_set(medium):
    # Wine and beer share the validated core + Tier-2 aroma byproducts (only the sugar
    # vector differs, and beer's sequential uptake lives inside the uptake Process). Wine
    # additionally carries the malolactic Process (D-23) and the opt-in carrying-capacity
    # modifier (D-30); beer does not (no malic/lactic; carrying cap is wine-only in v1).
    pset = get_medium(medium).build_process_set(strict=True)
    assert {p.name for p in pset.active} == EXPECTED_PROCESSES[medium]
    assert {m.name for m in pset.active_modifiers} == EXPECTED_MODIFIERS[medium]


def test_wine_growth_arrhenius_scales_the_amino_acid_swap_but_uptake_does_not():
    # The amino-acid swap (D-32) anchors to growth's base rate, so the GROWTH Arrhenius must scale
    # it too (else its refunds outrun the realised draw at M<1 and create sugar); the UPTAKE
    # Arrhenius must NOT (the swap tracks growth, not the fermentative sugar-uptake flux).
    mods = {m.name: m for m in get_medium("wine").build_process_set().active_modifiers}
    assert "amino_acid_assimilation" in mods["arrhenius_growth"].modifies
    assert "amino_acid_assimilation" not in mods["arrhenius_uptake"].modifies


def test_each_build_returns_fresh_kinetic_instances():
    # Factories, not shared instances: two builds must not hand back the same
    # Process objects (a shared mutable Process across runs/media would be a bug).
    a = get_medium("wine").build_process_set()
    b = get_medium("wine").build_process_set()
    a_procs = {p.name: p for p in a.active}
    b_procs = {p.name: p for p in b.active}
    assert all(a_procs[n] is not b_procs[n] for n in a_procs)
