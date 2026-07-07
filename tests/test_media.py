"""Medium state layouts and the medium registry."""

import dataclasses

import numpy as np
import pytest

from fermentation.core.media import MEDIA, Medium, beer_schema, get_medium, wine_schema

SHARED = (
    "X", "S", "E", "N", "T", "CO2", "X_dead", "Gly", "Byp", "esters", "fusels", "esters_gas",
    "acetolactate", "diacetyl", "butanediol", "acetaldehyde", "h2s", "h2s_gas",
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
# Brettanomyces volatile-phenol slots (decision D-40), appended last: the lumped hydroxycinnamic
# precursor, the shared vinylphenol intermediate reservoir, the ethylphenol readout, and the
# viable/dead Brett biomass pools (X_brett a constant catalyst in pt1; X_brett_dead filled by
# BrettDeath, pt3).
WINE_BRETT_SLOTS = (
    "hydroxycinnamics", "vinylphenols", "ethylphenols", "X_brett", "X_brett_dead",
)  # fmt: skip
# The carbon-bearing volatile-thiol pool AutolyticMercaptan fills (decision D-45), appended last.
WINE_MERCAPTAN_SLOTS = ("mercaptans",)
# The excreted keto-acid overflow pools (decisions D-49, D-50), appended last: pyruvate then
# alpha-ketoglutarate, the second- and third-strongest SO₂-binding carbonyls after acetaldehyde.
WINE_KETO_ACID_SLOTS = ("pyruvate", "alpha_ketoglutarate")


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
    )
    assert schema.spec("S").size == 1
    # 18 shared (X, S(1), E, N, T, CO2, X_dead, Gly, Byp, esters, fusels, esters_gas,
    # acetolactate, diacetyl, butanediol — the VDK pathway, D-26 — acetaldehyde, D-27, and
    # h2s + h2s_gas, D-29 production / D-42 CO2-stripping sink) + 3 wine-only acid slots
    # + citrate (D-31) + cation_charge (D-18) + 1 free-SO₂ slot (D-22) + X_mlf + X_mlf_dead
    # slots (D-23 catalyst / D-39 bacterial lees) + 1 amino_acids slot (D-32) + 1 debris slot
    # (D-34) + 5 Brett slots (hydroxycinnamics, vinylphenols, ethylphenols, X_brett,
    # X_brett_dead — decision D-40) + 1 mercaptans slot (D-45) + 2 keto-acid slots (pyruvate
    # D-49, alpha_ketoglutarate D-50)
    assert schema.size == 36


def test_beer_schema_has_three_sequential_sugars():
    schema = beer_schema()
    assert schema.names == SHARED
    s = schema.spec("S")
    assert s.size == 3
    assert s.components == ("glucose", "maltose", "maltotriose")
    # X, S(3), E, N, T, CO2, X_dead, Gly, Byp, esters, fusels, esters_gas (D-20),
    # acetolactate, diacetyl, butanediol (VDK pathway, D-26), acetaldehyde (D-27),
    # h2s (D-29) + h2s_gas (D-42 CO2-stripping sink)
    assert schema.size == 20


def test_shared_variable_units_are_canonical():
    # Beer carries exactly the shared variable set, so it pins the canonical units of
    # the shared layout (wine adds the D-18 acid/cation slots on top — checked below).
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
# brett_death (D-40 pt3) is pitch-gated too (Brett dies whether or not it grew), so it rides in this
# same tuple — the MalolacticDeath-in-_MLF_PROCESSES precedent.
BRETT_PROCESSES = {"brett_decarboxylation", "brett_vinylphenol_reduction", "brett_death"}
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
    ),
    "beer": (
        CORE_PROCESSES
        | TEMPERATURE_PROCESSES
        | BYPRODUCT_PROCESSES
        | VDK_PROCESSES
        | ACETALDEHYDE_PROCESSES
        | H2S_PROCESSES
    ),
}
# Per-rate Arrhenius modifiers wire into both media. Wine additionally carries the opt-in
# biomass carrying-capacity cap (decision D-30) — enabled in a bare build, but disabled at the
# compile seam unless a scenario passes carrying_capacity_gpl, so an undosed wine run is
# byte-for-byte the validated core (the wine-only MLF *modifier* parallel).
CARRYING_CAPACITY_MODIFIER = "biomass_carrying_capacity"
EXPECTED_MODIFIERS = {
    "wine": {"arrhenius_growth", "arrhenius_uptake", CARRYING_CAPACITY_MODIFIER},
    "beer": {"arrhenius_growth", "arrhenius_uptake"},
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
