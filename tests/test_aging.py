"""Tests for the Tier-3 aging Process :class:`EsterHydrolysis` (decision D-69).

The first §4.1 aging Process: a first-order **net decay** of the lumped ``esters`` pool
toward a lower equilibrium floor ``isoamyl_acetate_eq`` (young fruity acetate esters hydrolyse and
fade with age), warmed by an Arrhenius factor (warmer ages faster), routing the released
ester carbon **5:2** into ``isoamyl_alcohol`` (isoamyl alcohol, the alcohol product) and ``Byp``
(succinic-stand-in acetic acid, the acid product). These tests pin the closed-form
derivative and the exact 5:2 split; prove the properties the aging axis requires — **net
decay toward equilibrium** (zero at/below ``isoamyl_acetate_eq``, not decay-to-zero), **warmer ⇒
faster**, and an **on-ledger carbon transfer that closes ``total_carbon`` to machine
precision** (the D-68 "conservation is back in force" requirement, unlike the D-67 readout);
check the solver-undershoot guards; and confirm the tier floors at speculative and the
Process touches only ``esters``/``isoamyl_alcohol``/``Byp`` (no ``S``/``E``/``CO2`` — aging draws no
sugar). The scenario-level aging-phase wiring (the ``age N months`` verb + the reconfigure
enable) is D-70; here the Process is exercised directly via a hand-built ``ProcessSet`` (the
D-64 loss-Process pattern), off the fermentation ProcessSet so isolability is preserved.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.analysis import astringency_series, color_series, polymeric_pigment_series
from fermentation.core.acidbase import (
    ACID_STATE,
    bisulfite_fraction,
    bisulfite_so2_at_ph,
    build_pka_map,
    free_acetaldehyde,
    neutral_fraction,
    ph_of_state,
    solve_cation_charge,
)
from fermentation.core.chemistry import (
    CARBON_ATOMS,
    M_2_METHYLBUTANAL,
    M_2_METHYLPROPANAL,
    M_3_METHYLBUTANAL,
    M_ACETALDEHYDE,
    M_ALPHA_KETOBUTYRATE,
    M_CO2,
    M_ETHANOL,
    M_ISOAMYL_ACETATE,
    M_ISOAMYL_OH,
    M_MALIC,
    M_METHIONAL,
    M_O2,
    M_PHENYLACETALDEHYDE,
    M_SO2,
    M_SOTOLON,
    M_TARTARIC,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics import (
    AcetaldehydeBridgedCondensation,
    AnthocyaninFading,
    Caramelization,
    EllagitanninOxidation,
    EsterHydrolysis,
    EthylAcetateEsterification,
    EthylHexanoateHydrolysis,
    MaillardBrowning,
    MaillardStrecker,
    OakExtraction,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    SotolonAldolCondensation,
    StreckerDegradation,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
    TanninEthylTanninCondensation,
    TanninSelfPolymerization,
    ThermalAnthocyaninFade,
    arrhenius_factor,
)
from fermentation.core.kinetics.aging import (
    _ACETIC_ACID_CARBONS,
    _BYP_CARBON_SHARE,
    _FUSEL_CARBON_SHARE,
    _ISOAMYL_ALCOHOL_CARBONS,
    _MAILLARD_PRODUCTS,
    _SO2_PER_O2,
    _tartrate_hydrolysis_backbone,
)
from fermentation.core.kinetics.amino_acid_pools import (
    AMINO_ACID_SPECS,
    ASSIMILABLE_SPECS,
    GENERIC_SPECIES,
    assimilable_carbon_per_nitrogen,
)
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.media import beer_schema, get_medium, wine_schema
from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import Trajectory, simulate
from fermentation.units.convert import mgl_to_gpl
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_mass,
    total_nitrogen,
)
from tests.conftest import seed_amino_acids

#: Carbon fractions of the three pools the transfer touches (mirror the Process constants).
#: D-96: the hydrolysis debits isoamyl acetate ITSELF, not an ethyl-acetate stand-in.
_ESTER_C = carbon_mass_fraction("isoamyl_acetate")
_FUSEL_C = carbon_mass_fraction("isoamyl_alcohol")
_BYP_C = carbon_mass_fraction("succinic_acid")
#: The 5:2 split, from the isoamyl-acetate stand-in reaction (isoamyl alcohol 5 C : acetic 2 C).
_FUSEL_SHARE = 5.0 / 7.0
_BYP_SHARE = 2.0 / 7.0
#: Carbon fractions of the two pools the oxidation transfer moves carbon between (E → acetaldehyde).
_ETHANOL_C = carbon_mass_fraction("ethanol")
_ACET_C = carbon_mass_fraction("acetaldehyde")
#: Carbon fraction of ethyl acetate (C4) — the D-127 bidirectional esterification pool.
_ETHYL_ACETATE_C = carbon_mass_fraction("ethyl_acetate")
#: Carbon fraction of the ethylidene bridge (C2H4) — the D-80 split-ledger on-ledger capture
#: species.
_ETHYLIDENE_C = carbon_mass_fraction("ethylidene")
#: Carbon fractions of the Strecker pools + the amino-acid source (D-75), for the closure checks.
_METHIONAL_C = carbon_mass_fraction("methional")
_PHENYLACET_C = carbon_mass_fraction("phenylacetaldehyde")
_CO2_C = carbon_mass_fraction("CO2")
_AA_C = carbon_mass_fraction(AMINO_ACID_SPECIES)
_AA_N = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
#: The generic bucket's stand-in (glutamine, D-100) — the other half of every identity-agnostic
#: draw, and the C-richer end of the blend the D-89 denominator must clear.
_GENERIC_C = carbon_mass_fraction(GENERIC_SPECIES)
_GENERIC_N = nitrogen_mass_fraction(GENERIC_SPECIES)
# The four non-oxidative THERMAL Strecker aldehyde/sotolon carbon fractions (decision D-87).
_2MB_C = carbon_mass_fraction("2_methylbutanal")
_3MB_C = carbon_mass_fraction("3_methylbutanal")
_2MP_C = carbon_mass_fraction("2_methylpropanal")
_SOTOLON_C = carbon_mass_fraction("sotolon")
# Caramelization (decision D-88): the sugar (glucose) and melanoidin-carbon-park carbon fractions.
_GLUCOSE_C = carbon_mass_fraction("glucose")
_MELANOIDIN_C = carbon_mass_fraction("melanoidin")
# Beer's other two sugars (decision D-90, medium-agnostic caramelization): each caramelizes at its
# OWN carbon fraction (glucose/maltose/maltotriose differ), so the vectorized draw must weight the
# melanoidin transfer per component — these pin the beer carbon-closure test.
_MALTOSE_C = carbon_mass_fraction("maltose")
_MALTOTRIOSE_C = carbon_mass_fraction("maltotriose")
# MaillardBrowning (decision D-89): the N-bearing melanoidin park's carbon + nitrogen fractions, and
# arginine's nitrogen fraction — for the sized-draw closed form + carbon/nitrogen closure checks.
_MAILLARD_MELANOIDIN_C = carbon_mass_fraction("maillard_melanoidin")
_MAILLARD_MELANOIDIN_N = nitrogen_mass_fraction("maillard_melanoidin")


@pytest.fixture
def store():
    # Real wine parameters (T_ref, biomass_C_fraction, ...) MERGED with the aging.yaml
    # hydrolysis constants — the shared, medium-agnostic aging file (D-69). This mirrors the
    # compile seam D-70 wires; here the test loads it directly. acidbase.yaml (the pH-system
    # pKa constants) is now REQUIRED: since D-124 EsterHydrolysis solves the wine pH for its
    # first-order [H+] factor, exactly as a compiled wine scenario carries acidbase.yaml (D-18).
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "aging.yaml",
        default_data_dir() / "acidbase.yaml",
    )


@pytest.fixture
def params(store):
    return store.resolve()


def _aged_wine(schema: StateSchema, *, ester: float = 0.1, t: float = 293.15, **kw) -> FloatArray:
    """A finished, racked wine at the start of aging: yeast gone (X=0), dry (S=0), with the
    liquid ``esters`` pool pre-loaded (nothing produces it during aging — the Process only
    decays it). ``isoamyl_alcohol``/``Byp`` default to 0 so their aging gains are unambiguous.

    Carries a real acid load (``tartaric`` + ``cation_charge`` ⇒ pH ≈ 3.3, the ``_sulfited_wine``
    values) so ``ph_of_state`` solves into the wine range — REQUIRED since D-124 made the
    hydrolysis rate pH-dependent (an acid-free state would solve to ~pH 7 and crush the rate). A
    test wanting an exact target pH overrides ``cation_charge`` (see :func:`_wine_at_ph`)."""
    y = schema.pack({"X": 0.0, "S": [0.0], "E": 100.0, "N": 0.0, "T": t, "CO2": 0.0})
    y[schema.slice("isoamyl_acetate")] = ester
    if "cation_charge" in schema:  # wine only — the beer schema has no acid/pH slots (D-18)
        y[schema.slice("tartaric")] = 4.0
        y[schema.slice("cation_charge")] = 0.012
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def _wine_at_ph(
    schema: StateSchema,
    params: Mapping[str, float],
    target_ph: float,
    *,
    ester: float = 0.1,
    t: float = 293.15,
    tartaric: float | None = None,
) -> FloatArray:
    """An aged wine whose acid load + back-solved ``cation_charge`` reproduce ``target_ph`` exactly
    (the ``solve_cation_charge`` inverse-anchoring idiom, D-18) — for pinning the pH-dependent
    hydrolysis rate at a known pH (D-124/D-125).

    ``tartaric`` defaults to the sourced reference ``tartaric_ref_ester_hydrolysis`` (R&O's 7.5 g/L
    model solution) so that at ``target_ph = pH_ref`` the D-125 multi-species factor is exactly 1.0
    and the anchor is byte-for-byte. A test probing tartrate-dependence overrides it (D-125)."""
    if tartaric is None:
        tartaric = params["tartaric_ref_ester_hydrolysis"]
    totals = {"tartaric": tartaric / M_TARTARIC}
    cation = solve_cation_charge(totals, 0.0, build_pka_map(params), target_ph)
    return _aged_wine(schema, ester=ester, t=t, tartaric=tartaric, cation_charge=cation)


# -- metadata -----------------------------------------------------------------


def test_metadata():
    p = EsterHydrolysis()
    assert p.name == "ester_hydrolysis"
    # Speculative: the aging axis is the Tier-3 frontier (form sourced, magnitudes estimated).
    assert p.tier is Tier.SPECULATIVE
    # An on-ledger inter-pool transfer: decays esters, routes carbon to the alcohol
    # (isoamyl_alcohol)
    # and acid (Byp) products — never S/E/CO2 (aging draws no sugar, unlike the M2 producers).
    # Since D-115 the two label tracers ride along: the transfer carries its VALINE label as
    # well as its carbon, or an aging segment would dilute the alcohol pool's enrichment with
    # returned molecules silently booked as unlabelled. They carry no carbon weight, so the
    # "never S/E/CO2" statement above is untouched — the ledger still sees three pools.
    assert set(p.touches) == {
        "isoamyl_acetate",
        "isoamyl_alcohol",
        "Byp",
        "isoamyl_acetate_valine",
        "isoamyl_alcohol_valine",
    }
    assert set(p.reads) == {
        "k_ester_hydrolysis",
        "E_a_ester_hydrolysis",
        "isoamyl_acetate_eq",
        "pH_ref_ester_hydrolysis",  # the multi-species factor's reference pH (D-124)
        "r_h2t_ester_hydrolysis",  # k_H2T/k_H+, the undissociated-tartaric ratio (D-125)
        "r_ht_ester_hydrolysis",  # k_HT-/k_H+, the bitartrate ratio (D-125)
        "tartaric_ref_ester_hydrolysis",  # the reference tartrate the factor normalizes at (D-125)
        "T_ref",
    }


# -- closed form & the 5:2 split ----------------------------------------------


def test_derivative_matches_closed_form(params):
    schema = wine_schema()
    ester, t = 0.1, 298.15  # off T_ref so the Arrhenius factor bites
    y = _aged_wine(schema, ester=ester, t=t)
    d = EsterHydrolysis().derivatives(0.0, y, schema, params)

    f_t = arrhenius_factor(t, params["E_a_ester_hydrolysis"], params["T_ref"])
    # D-124/D-125: the rate carries the multi-species acid-catalysis factor. Recompute it INLINE
    # from R&O's backbone (the SulfiteOxidation-test idiom of rebuilding the driver from the solved
    # pH + the state) rather than calling the Process helper, so this is an independent check:
    # N(pH, T) = [H+] + r_h2t*[H2T] + r_ht*[HT-], normalized at (pH_ref, tartaric_ref).
    pkas = tuple(params[n] for n in ACID_STATE["tartaric"].pka_param_names)
    m_tart = ACID_STATE["tartaric"].molar_mass

    def _backbone(ph: float, tart_gpl: float) -> float:
        h = float(10.0**-ph)
        t_mol = tart_gpl / m_tart
        return float(
            h
            + params["r_h2t_ester_hydrolysis"] * neutral_fraction(h, pkas) * t_mol
            + params["r_ht_ester_hydrolysis"] * bisulfite_fraction(h, pkas) * t_mol
        )

    tart_wine = float(schema.get(y, "tartaric"))
    h_factor = _backbone(ph_of_state(y, schema, params), tart_wine) / _backbone(
        params["pH_ref_ester_hydrolysis"], params["tartaric_ref_ester_hydrolysis"]
    )
    rate = params["k_ester_hydrolysis"] * f_t * h_factor * (ester - params["isoamyl_acetate_eq"])
    carbon_released = rate * _ESTER_C

    assert schema.get(d, "isoamyl_acetate") == pytest.approx(-rate)
    # 5:2 split of the released carbon, re-deposited via each product's own carbon fraction.
    assert schema.get(d, "isoamyl_alcohol") == pytest.approx(
        _FUSEL_SHARE * carbon_released / _FUSEL_C
    )
    assert schema.get(d, "Byp") == pytest.approx(_BYP_SHARE * carbon_released / _BYP_C)
    # Aging touches nothing else — no sugar draw, no ethanol/CO2, no biomass.
    for var in ("X", "S", "E", "N", "CO2"):
        assert schema.get(d, var) == 0.0


def test_carbon_closes_per_rhs(params):
    # THE D-68 "conservation is back in force" invariant: the ester carbon lost equals the
    # fusel + Byp carbon gained, to machine precision — a pure on-ledger inter-pool transfer
    # (no S involvement), so total_carbon closes for ANY split summing to 1 (here 5:2).
    schema = wine_schema()
    d = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, ester=0.1, t=298.15), schema, params)
    carbon_residual = (
        schema.get(d, "isoamyl_acetate") * _ESTER_C
        + schema.get(d, "isoamyl_alcohol") * _FUSEL_C
        + schema.get(d, "Byp") * _BYP_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-15)


def test_five_to_two_split_exactly_partitions_the_debited_molecule():
    """The 5:2 split accounts for EVERY carbon of the ester it debits — impossible before D-96.

    This is the test that could not have been written at D-69. Back then the Process debited
    the lumped ``esters`` pool at its ledger-fixed **ethyl acetate** (C4) weighting, while
    splitting the released carbon 5:2 as though the molecule were **isoamyl acetate** — the
    stand-in reaction the split ratio actually comes from. Debited molecule ≠ split molecule, so
    ``5 + 2 = 7`` could not equal the debited molecule's 4 carbons, and D-69 could only
    *document* the mismatch. Crucially, every conservation test still passed: closure holds for
    ANY split summing to 1, so the ledger could not see the seam.

    D-96 split the lump, so the Process now debits isoamyl acetate itself and 5 + 2 = 7 = its
    real carbon count. The stand-in became the thing. (``aging.py`` also asserts this at import;
    pinned here too, since the seam it closes was invisible to conservation for 27 decisions.)
    """
    assert _ISOAMYL_ALCOHOL_CARBONS + _ACETIC_ACID_CARBONS == CARBON_ATOMS["isoamyl_acetate"] == 7
    # And the ratio really is read off those carbon counts, not a free parameter.
    assert pytest.approx(1.0) == _FUSEL_CARBON_SHARE + _BYP_CARBON_SHARE
    assert pytest.approx(5.0 / 2.0) == _FUSEL_CARBON_SHARE / _BYP_CARBON_SHARE


def test_split_is_five_to_two_by_carbon(params):
    # The carbon (not mass) partition is exactly 5:2 — isoamyl acetate's real hydrolysis
    # stoichiometry (isoamyl alcohol 5 C : acetic acid 2 C), the advisor-settled crux (D-69),
    # exact rather than a stand-in since D-96. Verified as carbon so it is independent of the
    # pools' differing mass weightings.
    schema = wine_schema()
    d = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, ester=0.1), schema, params)
    fusel_carbon = schema.get(d, "isoamyl_alcohol") * _FUSEL_C
    byp_carbon = schema.get(d, "Byp") * _BYP_C
    assert fusel_carbon / byp_carbon == pytest.approx(5.0 / 2.0)
    # Both products gain, and the fusel share is the larger (5/7) — the stronger fusel-OAV
    # rise the owner asked for (D-68 fork 2), plus the VA/pH-drifting Byp acid product.
    assert fusel_carbon > byp_carbon > 0.0


# -- net decay toward equilibrium (not to zero) -------------------------------


def test_zero_at_and_below_equilibrium(params):
    # Net decay toward a LOWER floor, not decay-to-zero (D-68): at or below isoamyl_acetate_eq the
    # rate is zero (the reverse ester-formation half is the deferred bidirectional term).
    schema = wine_schema()
    eq = params["isoamyl_acetate_eq"]
    at_eq = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, ester=eq), schema, params)
    below = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, ester=eq * 0.5), schema, params)
    assert np.array_equal(at_eq, schema.zeros())
    assert np.array_equal(below, schema.zeros())


def test_decays_only_the_excess_above_equilibrium(params):
    # The rate is proportional to (ester - isoamyl_acetate_eq), so a pool twice as far above it
    # decays twice as fast — the linear approach to equilibrium (Ramey & Ough first-order form).
    schema = wine_schema()
    eq = params["isoamyl_acetate_eq"]
    near = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, ester=eq + 0.02), schema, params)
    far = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, ester=eq + 0.04), schema, params)
    assert schema.get(far, "isoamyl_acetate") == pytest.approx(
        2.0 * schema.get(near, "isoamyl_acetate")
    )
    assert (
        schema.get(far, "isoamyl_acetate") < schema.get(near, "isoamyl_acetate") < 0.0
    )  # both decaying


def test_solver_undershoot_does_not_create_pools(params):
    # A solver undershoot (esters < 0) must not flip max(0, ...) into spurious production of
    # isoamyl_alcohol/Byp (or negative decay). A floor > 0 makes the excess negative ⇒ clamped to 0.
    schema = wine_schema()
    d = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, ester=-1e-6), schema, params)
    assert np.array_equal(d, schema.zeros())


# -- temperature direction (warmer ages faster) -------------------------------


def test_rises_with_temperature(params):
    # The sourced ordering (E_a_ester_hydrolysis > 0): warmer storage hydrolyses the esters
    # faster — why warm cellars age wine faster and cold storage preserves fruity esters.
    schema = wine_schema()
    cold = EsterHydrolysis().derivatives(
        0.0, _aged_wine(schema, ester=0.1, t=283.15), schema, params
    )
    warm = EsterHydrolysis().derivatives(
        0.0, _aged_wine(schema, ester=0.1, t=303.15), schema, params
    )
    # Faster decay (more negative) and a correspondingly larger fusel/Byp gain when warm.
    assert schema.get(warm, "isoamyl_acetate") < schema.get(cold, "isoamyl_acetate") < 0.0
    assert schema.get(warm, "isoamyl_alcohol") > schema.get(cold, "isoamyl_alcohol") > 0.0


def test_factor_is_one_at_reference_temperature(params):
    # At T_ref the Arrhenius factor is exactly 1, so the rate is the first-order term times only
    # the acid-catalysis factor (D-124/D-125). Build the wine at the REFERENCE pH AND reference
    # tartrate (the _wine_at_ph default) so h(pH, tartrate) is also exactly 1 and the rate is the
    # bare k*(ester - eq) — the byte-for-byte anchor at (T_ref, pH_ref, tartaric_ref).
    schema = wine_schema()
    ester = 0.1
    y = _wine_at_ph(
        schema, params, params["pH_ref_ester_hydrolysis"], ester=ester, t=params["T_ref"]
    )
    d = EsterHydrolysis().derivatives(0.0, y, schema, params)
    expected = params["k_ester_hydrolysis"] * (ester - params["isoamyl_acetate_eq"])
    assert schema.get(d, "isoamyl_acetate") == pytest.approx(-expected)


# -- Ramey & Ough 1980 re-anchor (decision D-123) -----------------------------
# R&O measured the isoamyl-acetate hydrolysis rate in REAL wine (open scanned PDF; the ACS copy is
# paywalled, which is why D-121 recorded it as blocked). Table IX (Pinot noir, pH 3.36): pseudo-
# first-order k_obsd = 54.72e-9 /s at 21.1 C, which shifts to ~1.80e-4 /h at T_ref = 20 C via E_a.
# Table X: E_a = 14.1 kcal/mol (Pinot) = 59.0 kJ/mol. These now anchor k_ester_hydrolysis /
# E_a_ester_hydrolysis, which were author estimates at D-69.
_RAMEY_OUGH_PINOT_KOBS_PER_H_20C = 1.80e-4  # R&O Table IX, Pinot noir, shifted to T_ref (/h)


def test_reanchored_params_sit_in_the_ramey_ough_measured_ranges(params):
    # Provenance teeth: the loaded k / E_a fall in R&O's measured WINE ranges, so a revert to the
    # D-69 author estimates would fail — notably k = 1.0e-4, which coincided with R&O's MODEL-
    # solution value (~1.06e-4) and sits below the real-wine band this re-anchor adopts.
    assert 55000.0 <= params["E_a_ester_hydrolysis"] <= 69000.0  # R&O wine 59-64, model 69 kJ/mol
    # D-124 NARROWED this band: the wine-to-wine pH spread that widened it to 4.0e-4 is now carried
    # by the explicit 10^(pH_ref - pH) factor, so k is the at-reference-pH rate (Pinot, grafted).
    assert 1.6e-4 <= params["k_ester_hydrolysis"] <= 2.5e-4  # R&O Pinot at pH 3.36, floor-grafted


def test_reanchored_rate_reproduces_ramey_ough_young_wine_fade(params):
    # THE floor-graft's design claim (D-123): EsterHydrolysis decays toward isoamyl_acetate_eq, but
    # R&O saw pure log-linear decay (NO floor over 200 d), so k was inflated (~x1.25) so that at a
    # representative young isoamyl_acetate level the sim's rate k_sim*(ester - eq) reproduces R&O's
    # measured disappearance k_obsd*[ester]. Check at 1.0 mg/L, T_ref (the level it was calibrated
    # at). With the old k = 1.0e-4 the sim rate would be ~44% of R&O's, so this pins the re-anchor.
    schema = wine_schema()
    ester = 1.0e-3  # g/L = 1.0 mg/L, a representative young-wine isoamyl acetate
    # R&O's k_obsd is the PINOT NOIR (pH 3.36) value, so evaluate at the reference pH where the
    # first-order [H+] factor is exactly 1.0 (D-124) — the fade is then the pinned k*(ester - eq).
    y = _wine_at_ph(
        schema, params, params["pH_ref_ester_hydrolysis"], ester=ester, t=params["T_ref"]
    )
    d = EsterHydrolysis().derivatives(0.0, y, schema, params)
    sim_fade = -float(schema.get(d, "isoamyl_acetate"))  # g/L/h, the modelled disappearance rate
    ramey_ough_fade = _RAMEY_OUGH_PINOT_KOBS_PER_H_20C * ester  # R&O's measured k_obsd*[ester]
    assert sim_fade == pytest.approx(ramey_ough_fade, rel=0.05)


# -- pH dependence: the multi-species acid-catalysis law, acetate fade tracks wine pH ---------
# (decisions D-124 [H+] backbone, D-125 tartrate terms). R&O 1980's full solved law is k_obsd =
# k_H+[H+] + k_H2T[H2T] + k_HT-[HT-] (Table VII). EsterHydrolysis scales its rate by the normalized
# factor h = N(pH, tartaric_wine)/N(pH_ref, tartaric_ref), N = [H+] + r_h2t[H2T] + r_ht[HT-]. A
# lower-pH wine fades its banana ester faster (the [H+] backbone dominates); h == 1 at the reference
# (pH_ref AND tartaric_ref, byte-for-byte the D-123 anchor) and in beer (no pH system, D-18).


def test_ph_ref_is_in_reads(params):
    # The factor reads pH_ref_ester_hydrolysis, so it is declared (tier propagation, D-1). The
    # plausible pH-system params (pKa_*, cation_charge) read inside ph_of_state / the tartaric
    # speciation are NOT declared — the Process is already speculative (the SulfiteOxidation
    # convention).
    assert "pH_ref_ester_hydrolysis" in EsterHydrolysis().reads
    assert "pH_ref_ester_hydrolysis" in params  # and it is actually loaded (aging.yaml)


def test_tartrate_law_params_in_reads(params):
    # D-125's three tartrate-law params are declared (so their speculative tiers propagate to the
    # output pools, D-1) and actually loaded (aging.yaml).
    for name in (
        "r_h2t_ester_hydrolysis",
        "r_ht_ester_hydrolysis",
        "tartaric_ref_ester_hydrolysis",
    ):
        assert name in EsterHydrolysis().reads
        assert name in params


def test_reference_gives_unit_factor_byte_for_byte(params):
    # At the reference (pH_ref AND tartaric_ref) the acid-catalysis factor is exactly 1.0, so the
    # rate is the pre-D-124 k*f_t*excess — the D-123 anchor and every reference trajectory are
    # preserved byte-for-byte (prime directive #3). _wine_at_ph defaults tartaric to tartaric_ref.
    schema = wine_schema()
    ester, t = 0.1, 298.15  # off T_ref so f_t bites and the check is non-vacuous
    y = _wine_at_ph(schema, params, params["pH_ref_ester_hydrolysis"], ester=ester, t=t)
    assert ph_of_state(y, schema, params) == pytest.approx(params["pH_ref_ester_hydrolysis"])
    assert float(schema.get(y, "tartaric")) == pytest.approx(
        params["tartaric_ref_ester_hydrolysis"]
    )
    d = EsterHydrolysis().derivatives(0.0, y, schema, params)
    f_t = arrhenius_factor(t, params["E_a_ester_hydrolysis"], params["T_ref"])
    expected = params["k_ester_hydrolysis"] * f_t * (ester - params["isoamyl_acetate_eq"])
    assert schema.get(d, "isoamyl_acetate") == pytest.approx(-expected)


def test_without_tartrate_the_backbone_is_first_order_in_hydrogen_ion(params):
    # D-125's [H+] BACKBONE: with no tartrate the law reduces to R&O's pure first-order form (the
    # D-124 behaviour). N(pH, 0) = [H+], so dropping pH by log10(2) doubles [H+] and doubles the
    # backbone exactly. (With tartrate present the law is DELIBERATELY not pure first-order — the
    # speciation adds curvature; that is what the high-pH ratio test below measures.) Tested at the
    # backbone helper because an acid-free wine cannot be brought to an acidic pH
    # (solve_cation_charge
    # rejects the negative strong-cation charge it would need).
    pkas = tuple(params[n] for n in ACID_STATE["tartaric"].pka_param_names)
    r_h2t, r_ht = params["r_h2t_ester_hydrolysis"], params["r_ht_ester_hydrolysis"]
    hi_ph = 3.6
    lo_ph = hi_ph - float(np.log10(2.0))  # [H+] doubled
    n_hi = _tartrate_hydrolysis_backbone(hi_ph, 0.0, r_h2t, r_ht, pkas)
    n_lo = _tartrate_hydrolysis_backbone(lo_ph, 0.0, r_h2t, r_ht, pkas)
    assert n_hi == pytest.approx(10.0**-hi_ph)  # no tartrate ⇒ backbone is [H+] alone
    assert n_lo == pytest.approx(2.0 * n_hi)  # first-order in [H+]


def test_high_ph_rate_ratio_matches_ramey_ough(params):
    # THE D-125 deliverable + its NON-CIRCULAR validation (not absolute-k, which would need R&O's
    # own
    # 12%-ethanol speciation): the model's k(pH 4.10)/k(pH 3.58) reproduces R&O's MEASURED ratio,
    # which the pure-[H+] law (D-124) could not. R&O Table V (isoamyl acetate, model solution):
    # k_obsd(4.10)/k_obsd(3.58) = 14.10/32.60 = 0.433. Pure-[H+] predicts [H+](4.10)/[H+](3.58) =
    # 10^(3.58-4.10) = 0.302. The multi-species law closes most of that gap (bitartrate catalysis at
    # high pH). Evaluated at R&O's model-solution tartrate (the default), the apples-to-apples
    # matrix.
    schema = wine_schema()
    ester, t = 0.1, 293.15

    def rate(ph: float) -> float:
        d = EsterHydrolysis().derivatives(
            0.0, _wine_at_ph(schema, params, ph, ester=ester, t=t), schema, params
        )
        return -float(schema.get(d, "isoamyl_acetate"))

    model_ratio = rate(4.10) / rate(3.58)
    ramey_ough_ratio = 14.10 / 32.60  # Table V measured k_obsd, isoamyl acetate
    pure_h_ratio = 10.0 ** (3.58 - 4.10)  # what D-124's [H+]-only law would give
    # Much closer to R&O than pure-[H+], and within ~5% of the measured ratio (the residual is the
    # sim's aqueous pKa vs R&O's 12%-ethanol speciation — the documented, ratio-washed mismatch).
    assert model_ratio == pytest.approx(ramey_ough_ratio, rel=0.06)
    assert abs(model_ratio - ramey_ough_ratio) < abs(pure_h_ratio - ramey_ough_ratio)


def test_lower_ph_hydrolyses_faster(params):
    # THE headline observable (the pH-blind model could not express it): a lower-pH wine fades its
    # banana ester faster (R&O: "pH is far more important than total acidity"; the [H+] backbone
    # dominates over the wine range). Spans the wine range at the reference tartrate.
    schema = wine_schema()
    low = EsterHydrolysis().derivatives(
        0.0, _wine_at_ph(schema, params, 3.0, ester=0.1), schema, params
    )
    high = EsterHydrolysis().derivatives(
        0.0, _wine_at_ph(schema, params, 3.8, ester=0.1), schema, params
    )
    assert -schema.get(low, "isoamyl_acetate") > -schema.get(high, "isoamyl_acetate") > 0.0
    # The fusel/Byp products rise correspondingly more in the faster (lower-pH) wine.
    assert schema.get(low, "isoamyl_alcohol") > schema.get(high, "isoamyl_alcohol") > 0.0


def test_off_reference_tartrate_changes_the_rate(params):
    # D-125's NEW state-dependence (D-124 was pH-only, tartrate-blind): at a FIXED high pH — where
    # bitartrate (HT-) is the dominant tartrate species — MORE tartaric acid hydrolyses FASTER.
    # This pins that the factor reads the `tartaric` state slot, and it is the honest supersession
    # of
    # D-124's "byte-for-byte at pH 3.36 regardless of tartrate": off the reference tartrate the rate
    # diverges. (At LOW pH the negative k_H2T flips this; tested here at pH 3.8 where the sign is
    # clean.)
    schema = wine_schema()
    ph = 3.8

    def rate(tart: float) -> float:
        d = EsterHydrolysis().derivatives(
            0.0, _wine_at_ph(schema, params, ph, ester=0.1, tartaric=tart), schema, params
        )
        return -float(schema.get(d, "isoamyl_acetate"))

    # Ascending tartrate ⇒ ascending rate at high pH (bitartrate catalysis grows with [tartrate]).
    assert 0.0 < rate(2.0) < rate(params["tartaric_ref_ester_hydrolysis"]) < rate(10.0)


def test_beer_ester_hydrolysis_is_ph_blind_byte_for_byte(params):
    # Beer has no pH system (D-18), so the [H+] factor is held at 1.0 (the cation_charge gate) and
    # the beer rate is the pH_ref-anchored k*f_t*excess — byte-for-byte the pre-D-124 beer behaviour
    # (prime directive #3: the pH build must not silently perturb the other medium).
    beer = beer_schema()
    ester, t = 0.08, 298.15
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": t, "CO2": 0.0})
    yb[beer.slice("isoamyl_acetate")] = ester
    d = EsterHydrolysis().derivatives(0.0, yb, beer, params)
    f_t = arrhenius_factor(t, params["E_a_ester_hydrolysis"], params["T_ref"])
    expected = params["k_ester_hydrolysis"] * f_t * (ester - params["isoamyl_acetate_eq"])
    assert beer.get(d, "isoamyl_acetate") == pytest.approx(-expected)


def test_integrated_lower_ph_wine_fades_more_over_a_year(params):
    # END-TO-END (the deliverable): over a year of aging, a lower-pH wine loses MORE of its banana
    # ester than a higher-pH wine at the same temperature — acetate fade tracks wine pH (D-124).
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()], strict=True)
    ester0, span = 0.1, (0.0, 24.0 * 365.0)
    low = simulate(
        ps, params=params, y0=_wine_at_ph(schema, params, 3.0, ester=ester0), t_span=span
    )
    high = simulate(
        ps, params=params, y0=_wine_at_ph(schema, params, 3.8, ester=ester0), t_span=span
    )
    assert low.success and high.success, (low.message, high.message)
    low_end = float(low.series("isoamyl_acetate")[-1])
    high_end = float(high.series("isoamyl_acetate")[-1])
    assert low_end < high_end < ester0  # both fade; the acidic wine fades further


# -- integrated aging segment (conservation + direction) ----------------------


def test_integrated_aging_closes_carbon_and_fades_esters(params, store):
    # Run a long aging segment (a racked, dry wine — X=0, S=0) with ONLY EsterHydrolysis
    # active, under the strict touches contract. Over the aging span the esters pool fades,
    # the isoamyl_alcohol and Byp pools rise, and total_carbon closes to machine precision —
    # the pure
    # on-ledger transfer of the D-68 aging axis (no sugar drawn, no ethanol touched).
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()], strict=True)
    esters0 = 0.1
    y0 = _aged_wine(schema, ester=esters0, t=293.15)
    # ~1 year of aging (large steps are fine — no ferment stiffness; the §7 slow-phase point).
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message

    # The fruity esters fade (toward, not past, the equilibrium floor) and the products rise.
    ester_end = float(traj.series("isoamyl_acetate")[-1])
    assert params["isoamyl_acetate_eq"] <= ester_end < esters0
    assert float(traj.series("isoamyl_alcohol")[-1]) > 0.0
    assert float(traj.series("Byp")[-1]) > 0.0
    # Non-negative pools and machine-precision carbon closure (X=0 throughout, so the biomass
    # term is inert; the invariant is the esters → isoamyl_alcohol + Byp inter-pool transfer).
    assert_nonnegative(traj, ("isoamyl_acetate", "isoamyl_alcohol", "Byp"), atol=1e-12)
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    # Carbon is the invariant; mass carries a small (~4.5%) documented stand-in gap (aging.yaml):
    # splitting a carbon-exact budget across pools with heterogeneous fixed weightings is not
    # mass-conserving (real hydrolysis consumes untracked water — the D-8/D-16/D-26 precedent).
    # total_mass weights only {S, E, CO2}, NONE of which this Process touches, so the gap is
    # scoped OUT by construction: the validated-core mass check stays flat through the aging span.
    assert_conserved(traj, total_mass(schema), rtol=1e-9, atol=1e-9, label="mass")


def test_isolable_from_the_core_when_below_equilibrium(params):
    # Isolability corner (prime directive #3): with the ester pool below its floor the Process
    # contributes exactly nothing, so an aging segment on an ester-poor wine is byte-for-byte
    # the no-aging state — the Process cannot create aroma out of an empty pool.
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()], strict=True)
    y = _aged_wine(schema, ester=params["isoamyl_acetate_eq"] * 0.5)
    assert np.array_equal(ps.total_derivatives(0.0, y, params), schema.zeros())


def test_integrated_aging_closes_carbon_beer_multislot(store):
    # The multi-slot (beer) counterpart: the aging transfer is sugar-free so the 3-slot S
    # vector is irrelevant to it, but running the strict ProcessSet on beer_schema proves the
    # Process is medium-agnostic (esters/isoamyl_alcohol/Byp exist in both) and closes carbon
    # there too.
    beer = load_parameters(
        default_data_dir() / "beer_generic.yaml", default_data_dir() / "aging.yaml"
    )
    params = beer.resolve()
    schema = beer_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()], strict=True)
    y0 = schema.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    y0[schema.slice("isoamyl_acetate")] = 0.08
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 180.0))
    assert traj.success, traj.message
    assert float(traj.series("isoamyl_acetate")[-1]) < 0.08
    f_c = beer.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


# ==============================================================================================
# ETHYL HEXANOATE HYDROLYSIS (decision D-126) — the apple/pineapple ethyl ester fades on aging, the
# sibling of EsterHydrolysis. Real-wine kinetics from Makhotkina & Kilmartin 2012 (PMID 22868118).
# Floored + grafted (the isoamyl D-123 idiom); NO pH factor (deferred, D-126). Routes the released
# carbon 2:6 → ethanol (core E slot, the OxidativeAcetaldehyde precedent) + hexanoic acid (Byp).
# ==============================================================================================

#: Makhotkina & Kilmartin 2012 ethyl-hexanoate hydrolysis k_obs at 20 C, /h. The Table 2
#: interpolation (18↔28 C) and the Table 3 Arrhenius fit (E_a=68 kJ/mol, A=4e4/s) agree at
#: ~1.1e-4 /h. The grafted k_sim reproduces this k_obs·[ester] at the young level (D-126).
_MAKHOTKINA_ETHYL_HEXANOATE_KOBS_PER_H_20C = 1.10e-4
#: The representative young ethyl-hexanoate level the graft is calibrated at (~0.4 mg/L, the sim's
#: calibrated young value; sensory.yaml/wine_generic.yaml). Mirrors isoamyl's 1.0 mg/L point.
_ETHYL_HEXANOATE_YOUNG_REF = 4.0e-4  # g/L


def test_ethyl_hexanoate_metadata():
    p = EthylHexanoateHydrolysis()
    assert p.name == "ethyl_hexanoate_hydrolysis"
    # Speculative: the aging axis is the Tier-3 frontier; here an honest floor (poor Arrhenius fit).
    assert p.tier is Tier.SPECULATIVE
    # An on-ledger inter-pool transfer: decays ethyl_hexanoate, routes carbon to ETHANOL (the core
    # E slot) and hexanoic acid (Byp). Unlike EsterHydrolysis it touches E — a total_mass{S,E,CO2}
    # pool (the OxidativeAcetaldehyde precedent) — and carries no fusel/label pools (ethyl hexanoate
    # has no valine label and its alcohol product is bulk ethanol, not a fusel).
    assert set(p.touches) == {"ethyl_hexanoate", "E", "Byp"}


def test_ethyl_hexanoate_reads_declared_and_loaded(params):
    p = EthylHexanoateHydrolysis()
    for name in (
        "k_ethyl_hexanoate_hydrolysis",
        "E_a_ethyl_hexanoate_hydrolysis",
        "ethyl_hexanoate_eq",
        "T_ref",
    ):
        assert name in p.reads  # declared (tier propagation, D-1)
        assert name in params  # and actually loaded (aging.yaml)
    # NO pH-system params are read — the pH/tartrate catalysis is DEFERRED (D-126), unlike
    # EsterHydrolysis which reads pH_ref_ester_hydrolysis + the D-125 tartrate ratios.
    assert not any(r.startswith(("pH_ref", "r_h2t", "r_ht", "tartaric_ref")) for r in p.reads)


def test_ethyl_hexanoate_params_in_makhotkina_ranges(params):
    # Provenance teeth: k / E_a fall in Makhotkina & Kilmartin 2012's measured ranges (k grafted).
    assert 9.0e-5 <= params["k_ethyl_hexanoate_hydrolysis"] <= 2.5e-4
    assert 38000.0 <= params["E_a_ethyl_hexanoate_hydrolysis"] <= 98000.0  # 68 ± 30 kJ/mol
    # The floor sits below the young apple level so most of it fades (a residuum remains).
    assert 0.0 < params["ethyl_hexanoate_eq"] < _ETHYL_HEXANOATE_YOUNG_REF


def test_ethyl_hexanoate_rate_reproduces_makhotkina_young_wine_fade(params):
    # THE floor-graft's design claim (D-126, the isoamyl D-123 idiom): EthylHexanoateHydrolysis
    # decays toward ethyl_hexanoate_eq, but Makhotkina's k_obs is a floor-less disappearance
    # constant, so k was grafted (~×1.33) so that at the ~0.4 mg/L young level the sim's rate
    # k_sim·(ester − eq) reproduces the observed k_obs·[ester]. No pH factor, so a bare wine does.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, t=params["T_ref"], ethyl_hexanoate=_ETHYL_HEXANOATE_YOUNG_REF)
    d = EthylHexanoateHydrolysis().derivatives(0.0, y, schema, params)
    sim_fade = -float(schema.get(d, "ethyl_hexanoate"))  # g/L/h, the modelled disappearance rate
    makhotkina_fade = _MAKHOTKINA_ETHYL_HEXANOATE_KOBS_PER_H_20C * _ETHYL_HEXANOATE_YOUNG_REF
    assert sim_fade == pytest.approx(makhotkina_fade, rel=0.05)


def test_ethyl_hexanoate_carbon_split_is_two_to_six(params):
    # The released C8 carbon splits 2:6 → ethanol (E, 1/4) + hexanoic acid (Byp, 3/4), re-deposited
    # through each pool's own carbon fraction so the transfer is carbon-exact for any split summing
    # to 1. Pin the exact routing (mirror of the isoamyl 5:2 derivative test).
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, t=params["T_ref"], ethyl_hexanoate=_ETHYL_HEXANOATE_YOUNG_REF)
    d = EthylHexanoateHydrolysis().derivatives(0.0, y, schema, params)
    rate = -float(schema.get(d, "ethyl_hexanoate"))
    carbon_released = rate * carbon_mass_fraction("ethyl_hexanoate")
    assert schema.get(d, "E") == pytest.approx(0.25 * carbon_released / _ETHANOL_C)
    assert schema.get(d, "Byp") == pytest.approx(0.75 * carbon_released / _BYP_C)
    # Carbon-exact: carbon out of the ester equals carbon into E + Byp.
    c_in = schema.get(d, "E") * _ETHANOL_C + schema.get(d, "Byp") * _BYP_C
    assert c_in == pytest.approx(carbon_released, rel=1e-12)


def test_ethyl_hexanoate_no_ph_factor(params):
    # D-126 minimal choice: NO pH/tartrate catalysis (deferred). Unlike EsterHydrolysis (whose rate
    # rises at low pH, D-124), this rate is pH-INDEPENDENT — the fade at pH 3.0 equals the fade at
    # pH 4.0 for the same pool. Anchored at Makhotkina's wine pH; the pH refinement is the named
    # follow-on. (Contrast the isoamyl low-pH-fades-faster behaviour.)
    schema = wine_schema()
    proc = EthylHexanoateHydrolysis()
    low = _wine_at_ph(schema, params, 3.0, ester=0.0, t=params["T_ref"])
    low[schema.slice("ethyl_hexanoate")] = _ETHYL_HEXANOATE_YOUNG_REF
    high = _wine_at_ph(schema, params, 4.0, ester=0.0, t=params["T_ref"])
    high[schema.slice("ethyl_hexanoate")] = _ETHYL_HEXANOATE_YOUNG_REF
    fade_low = -float(schema.get(proc.derivatives(0.0, low, schema, params), "ethyl_hexanoate"))
    fade_high = -float(schema.get(proc.derivatives(0.0, high, schema, params), "ethyl_hexanoate"))
    assert fade_low == pytest.approx(fade_high, rel=1e-12)


def test_ethyl_hexanoate_warmer_ages_faster(params):
    # The load-bearing DIRECTION (E_a > 0): a warmer wine fades its apple ester faster (cold storage
    # preserves it). Makhotkina's headline, and the only claim the poor Arrhenius fit (r²=0.572)
    # robustly supports.
    schema = wine_schema()
    proc = EthylHexanoateHydrolysis()
    cold = _aged_wine(schema, ester=0.0, t=283.15, ethyl_hexanoate=_ETHYL_HEXANOATE_YOUNG_REF)
    warm = _aged_wine(schema, ester=0.0, t=303.15, ethyl_hexanoate=_ETHYL_HEXANOATE_YOUNG_REF)
    fade_cold = -float(schema.get(proc.derivatives(0.0, cold, schema, params), "ethyl_hexanoate"))
    fade_warm = -float(schema.get(proc.derivatives(0.0, warm, schema, params), "ethyl_hexanoate"))
    assert fade_warm > fade_cold > 0.0


def test_ethyl_hexanoate_isolable_below_equilibrium(params):
    # Isolability (prime directive #3): below the floor the Process contributes exactly nothing —
    # an aging segment on an apple-ester-poor wine is byte-for-byte the no-aging state.
    schema = wine_schema()
    ps = ProcessSet(schema, [EthylHexanoateHydrolysis()], strict=True)
    y = _aged_wine(schema, ester=0.0, ethyl_hexanoate=params["ethyl_hexanoate_eq"] * 0.5)
    assert np.array_equal(ps.total_derivatives(0.0, y, params), schema.zeros())


def test_ethyl_hexanoate_aging_closes_carbon_and_fades_apple(params, store):
    # Integrated aging segment (racked dry wine, X=0, S=0) with ONLY EthylHexanoateHydrolysis active
    # under the strict touches contract. The apple ester fades toward (not past) its floor, E and
    # Byp rise, and total_CARBON closes to machine precision — the on-ledger transfer
    # ethyl_hexanoate → ethanol (E) + hexanoic acid (Byp).
    schema = wine_schema()
    ps = ProcessSet(schema, [EthylHexanoateHydrolysis()], strict=True)
    eh0 = _ETHYL_HEXANOATE_YOUNG_REF
    y0 = _aged_wine(schema, ester=0.0, t=293.15, ethyl_hexanoate=eh0)
    e0 = float(schema.get(y0, "E"))
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    eh_end = float(traj.series("ethyl_hexanoate")[-1])
    assert params["ethyl_hexanoate_eq"] <= eh_end < eh0  # fades toward, not past, the floor
    assert float(traj.series("E")[-1]) > e0  # ethanol released into the core E slot
    assert float(traj.series("Byp")[-1]) > 0.0  # hexanoic acid → Byp (succinic stand-in)
    assert_nonnegative(traj, ("ethyl_hexanoate", "Byp"), atol=1e-12)
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    # total_mass is DELIBERATELY NOT asserted flat here (contrast the isoamyl
    # test_integrated_aging_closes_carbon_and_fades_esters, which does): E is in the
    # total_mass{S,E,CO2} sub-ledger, and this Process ADDS ethanol to E from ethyl_hexanoate (a
    # pool OUTSIDE the sub-ledger), so total_mass gains a hair — exactly the OxidativeAcetaldehyde
    # precedent (which loses E-mass to acetaldehyde). The gain corresponds to real untracked
    # hydrolysis water (the D-8/D-16/D-26 stand-in gap) and is ~1e-6 g/L even over a year; the
    # invariant this Process closes exactly is total_CARBON.


def test_ethyl_hexanoate_wired_into_both_media():
    # The Process is wired into the medium-agnostic _AGING_PROCESSES (like EsterHydrolysis), so both
    # a compiled wine and beer carry it — the ethyl_hexanoate pool exists in both (D-96).
    for medium in ("wine", "beer"):
        names = {p().name for p in get_medium(medium).process_factories}
        assert "ethyl_hexanoate_hydrolysis" in names


# =====================================================================================
# EthylAcetateEsterification (decision D-127) — the THIRD ester Process and the ONLY bidirectional
# one: ethyl acetate sits ~AT its esterification equilibrium, so it FADES toward the floor from
# above (net hydrolysis → ethanol + acetic acid) and FORMS toward it from below (net esterification
# ← ethanol + acetic acid). Model-derived speculative rate + equilibrium (Shinohara 1979 approach
# time + ~10% E-rate; R&O 1980 acetate-cluster cross-check). These tests pin the metadata, the
# signed carbon closure for BOTH directions, the sign flip about eq, the integrated fade and form,
# and the both-media wiring.


def test_ethyl_acetate_esterification_metadata():
    p = EthylAcetateEsterification()
    assert p.name == "ethyl_acetate_esterification"
    assert p.tier is Tier.SPECULATIVE
    # Signed inter-pool transfer: ethyl_acetate <=> ethanol (E) + acetic acid (Byp, succinic).
    assert set(p.touches) == {"ethyl_acetate", "E", "Byp"}
    assert set(p.reads) == {
        "k_ethyl_acetate_esterification",
        "E_a_ethyl_acetate_esterification",
        "ethyl_acetate_eq",
        "pH_ref_ethyl_acetate_esterification",
        "T_ref",
    }


def test_ethyl_acetate_eq_sits_near_the_young_level(params):
    # Unlike the two hydrolysis floors (~25% of their young ester, far ABOVE equilibrium), ethyl
    # acetate's equilibrium is ~AT the sim's calibrated ~50 mg/L young level — that is the whole
    # reason this Process is bidirectional (a sound wine barely moves; it acts on off-equilibrium
    # wines). Pin that the floor is in the sound-wine 30–80 mg/L band.
    assert mgl_to_gpl(30.0) < params["ethyl_acetate_eq"] < mgl_to_gpl(80.0)


@pytest.mark.parametrize("side", ["above", "below"])
def test_ethyl_acetate_carbon_closes_per_rhs_both_directions(params, side):
    # The signed C4 ⇌ C2 (ethanol) + C2 (acetic) transfer closes total_carbon to machine precision
    # for EITHER flux sign — fading (ester>eq) and forming (ester<eq). This is the D-127 crux: the
    # SAME split releases/deposits, so closure holds regardless of direction.
    schema = wine_schema()
    eq = params["ethyl_acetate_eq"]
    ester = eq * 2.0 if side == "above" else eq * 0.5
    # Seed a realistic Byp so a forming flux (which DEBITS Byp) has acetic to draw from.
    y = _aged_wine(schema, ester=0.0, t=298.15, ethyl_acetate=ester, Byp=1.0)
    d = EthylAcetateEsterification().derivatives(0.0, y, schema, params)
    # Carbon leaving ethyl_acetate == carbon entering E + Byp (signed), to machine precision.
    residual = (
        schema.get(d, "ethyl_acetate") * _ETHYL_ACETATE_C
        + schema.get(d, "E") * _ETHANOL_C
        + schema.get(d, "Byp") * _BYP_C
    )
    assert residual == pytest.approx(0.0, abs=1e-15)
    # Direction check: above eq fades (ester↓, E↑, Byp↑); below eq forms (ester↑, E↓, Byp↓).
    if side == "above":
        assert schema.get(d, "ethyl_acetate") < 0.0
        assert schema.get(d, "E") > 0.0 and schema.get(d, "Byp") > 0.0
    else:
        assert schema.get(d, "ethyl_acetate") > 0.0
        assert schema.get(d, "E") < 0.0 and schema.get(d, "Byp") < 0.0
    # Touches nothing else — no sugar, CO2, biomass, or the other ester pools.
    for var in ("X", "S", "N", "CO2", "isoamyl_acetate", "isoamyl_alcohol", "ethyl_hexanoate"):
        assert schema.get(d, var) == 0.0


def test_ethyl_acetate_inert_exactly_at_equilibrium(params):
    # At ethyl_acetate == eq the signed gap is exactly zero, so the Process contributes
    # byte-for-byte nothing — the pivot of the bidirectional relaxation.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, t=298.15, ethyl_acetate=params["ethyl_acetate_eq"])
    d = EthylAcetateEsterification().derivatives(0.0, y, schema, params)
    for var in ("ethyl_acetate", "E", "Byp"):
        assert schema.get(d, var) == 0.0


def test_ethyl_acetate_aging_fades_high_va_wine_and_closes_carbon(params, store):
    # Integrated aging of a HIGH-EtOAc wine (above eq): the solventy note fades toward — not past —
    # its equilibrium (Shinohara's observed EtOAc decrease in stored wine), E and Byp rise, and
    # total_CARBON closes to machine precision. Strict touches contract, only this Process active.
    schema = wine_schema()
    ea0 = params["ethyl_acetate_eq"] * 2.0  # a high-VA / high-EtOAc wine, above equilibrium
    y0 = _aged_wine(schema, ester=0.0, t=293.15, ethyl_acetate=ea0, Byp=1.0)
    e0 = float(schema.get(y0, "E"))
    ps = ProcessSet(schema, [EthylAcetateEsterification()], strict=True)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    ea_end = float(traj.series("ethyl_acetate")[-1])
    assert params["ethyl_acetate_eq"] <= ea_end < ea0  # fades toward, not past, equilibrium
    assert float(traj.series("E")[-1]) > e0  # ethanol released to core E
    assert float(traj.series("Byp")[-1]) > 1.0  # acetic acid → Byp (succinic stand-in)
    assert_nonnegative(traj, ("ethyl_acetate", "Byp"), atol=1e-12)
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_ethyl_acetate_aging_forms_in_below_equilibrium_wine(params, store):
    # The FORMATION half — the one the sim models for this ester alone (D-127). A below-equilibrium
    # wine's ethyl acetate RISES toward eq, consuming ethanol (E↓) and acetic acid (Byp↓), and
    # total_CARBON still closes. Byp is seeded realistically so the tiny acetic draw won't go under.
    schema = wine_schema()
    ea0 = params["ethyl_acetate_eq"] * 0.5  # below equilibrium → forms
    y0 = _aged_wine(schema, ester=0.0, t=293.15, ethyl_acetate=ea0, Byp=1.0)
    e0, byp0 = float(schema.get(y0, "E")), float(schema.get(y0, "Byp"))
    ps = ProcessSet(schema, [EthylAcetateEsterification()], strict=True)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    ea_end = float(traj.series("ethyl_acetate")[-1])
    assert ea0 < ea_end <= params["ethyl_acetate_eq"]  # forms toward, not past, equilibrium
    assert float(traj.series("E")[-1]) < e0  # ethanol consumed by esterification
    assert float(traj.series("Byp")[-1]) < byp0  # acetic acid consumed
    assert_nonnegative(traj, ("ethyl_acetate", "Byp"), atol=1e-12)
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_ethyl_acetate_wired_into_both_media():
    # Wired into the medium-agnostic _AGING_PROCESSES (like the two hydrolysis siblings); the
    # ethyl_acetate pool exists in both media, and h(pH) is held at 1 in beer (no pH system, D-18).
    for medium in ("wine", "beer"):
        names = {p().name for p in get_medium(medium).process_factories}
        assert "ethyl_acetate_esterification" in names


# -- tier propagation ---------------------------------------------------------


def test_tier_floored_at_speculative(store):
    # The aging Process is speculative in FORM (Tier-3 frontier), so every pool it writes is
    # speculative even before parameters cap it — and folding in the (speculative) aging
    # parameter tiers keeps it there. Non-vacuous: esters/isoamyl_alcohol/Byp are all speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()])
    for pool in ("isoamyl_acetate", "isoamyl_alcohol", "Byp"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# OxidativeAcetaldehyde (decision D-71) — the first OXIDATIVE aging Process: dissolved O₂
# oxidises ethanol → acetaldehyde on the ``o2`` substrate pool. O₂-limited (first-order in the
# finite o2 pool ⇒ SATURATING, not the unbounded ethanol-first alternative), Arrhenius warmer-
# faster, a molar yield ``y_acetaldehyde_per_o2`` of the consumed O₂ becoming acetaldehyde and the
# oxidised carbon borrowed carbon-exactly from ``E`` (the D-27 reduction reversed). These tests pin
# the closed form, the O₂-off-ledger carbon closure (E → acetaldehyde), the reductive-aging
# isolability (inert at ``o2 = 0``), the warmer-faster ordering, the saturating depletion over an
# integrated segment, and the speculative tier floor.


def test_oxidation_metadata():
    p = OxidativeAcetaldehyde()
    assert p.name == "oxidative_acetaldehyde"
    assert p.tier is Tier.SPECULATIVE
    # Consumes the O₂ substrate, books the oxidised carbon as acetaldehyde borrowed from E.
    assert set(p.touches) == {"o2", "acetaldehyde", "E"}
    assert set(p.reads) == {
        "k_ethanol_oxidation",
        "E_a_ethanol_oxidation",
        "y_acetaldehyde_per_o2",
        "T_ref",
    }


def test_oxidation_matches_closed_form(params):
    schema = wine_schema()
    o2, t = 0.03, 298.15  # off T_ref so the Arrhenius factor bites
    y = _aged_wine(schema, ester=0.0, t=t, o2=o2)  # esters=0 so only oxidation moves anything
    d = OxidativeAcetaldehyde().derivatives(0.0, y, schema, params)

    f_t = arrhenius_factor(t, params["E_a_ethanol_oxidation"], params["T_ref"])
    r_o2 = params["k_ethanol_oxidation"] * f_t * o2
    acet_rate = params["y_acetaldehyde_per_o2"] * (r_o2 / M_O2) * M_ACETALDEHYDE

    assert schema.get(d, "o2") == pytest.approx(-r_o2)
    assert schema.get(d, "acetaldehyde") == pytest.approx(acet_rate)
    # Carbon-exact C2 borrow from ethanol (the reduction reversed).
    assert schema.get(d, "E") == pytest.approx(-acet_rate * M_ETHANOL / M_ACETALDEHYDE)
    # Oxidation touches nothing else — no sugar, no CO2, no esters/isoamyl_alcohol/Byp, no biomass.
    for var in ("X", "S", "N", "CO2", "isoamyl_acetate", "isoamyl_alcohol", "Byp"):
        assert schema.get(d, var) == 0.0


def test_oxidation_carbon_closes_per_rhs(params):
    # O₂ is OFF every ledger, so the only on-ledger movement is E → acetaldehyde, both C2 — the
    # carbon lost from ethanol equals the carbon gained as acetaldehyde, to machine precision.
    schema = wine_schema()
    d = OxidativeAcetaldehyde().derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=298.15, o2=0.03), schema, params
    )
    carbon_residual = schema.get(d, "E") * _ETHANOL_C + schema.get(d, "acetaldehyde") * _ACET_C
    assert carbon_residual == pytest.approx(0.0, abs=1e-15)


def test_oxidation_inert_without_oxygen(params):
    # Reductive aging (screwcap/inert) + the exact isolability guard: with no dissolved O₂ the
    # Process contributes byte-for-byte zero, so a begin_aging run with no add_oxygen is the
    # ester-hydrolysis-only case. Cannot oxidise ethanol out of an empty O₂ pool.
    schema = wine_schema()
    ps = ProcessSet(schema, [OxidativeAcetaldehyde()], strict=True)
    y = _aged_wine(schema, ester=0.0, o2=0.0)
    assert np.array_equal(ps.total_derivatives(0.0, y, params), schema.zeros())


def test_oxidation_solver_undershoot_does_not_create_acetaldehyde(params):
    # A solver undershoot (o2 < 0) must not flip into spurious acetaldehyde production: the
    # ``o2 <= 0`` guard returns zeros (no oxidant ⇒ no oxidation).
    schema = wine_schema()
    d = OxidativeAcetaldehyde().derivatives(0.0, _aged_wine(schema, o2=-1e-6), schema, params)
    assert np.array_equal(d, schema.zeros())


def test_oxidation_rises_with_temperature(params):
    # The sourced ordering (E_a_ethanol_oxidation > 0): warmer storage oxidises (maderises) faster
    # — more O₂ consumed and more acetaldehyde made per hour when warm.
    schema = wine_schema()
    cold = OxidativeAcetaldehyde().derivatives(
        0.0, _aged_wine(schema, o2=0.03, t=283.15), schema, params
    )
    warm = OxidativeAcetaldehyde().derivatives(
        0.0, _aged_wine(schema, o2=0.03, t=303.15), schema, params
    )
    # Warmer ⇒ faster O₂ depletion (more negative) and a larger acetaldehyde gain.
    assert schema.get(warm, "o2") < schema.get(cold, "o2") < 0.0
    assert schema.get(warm, "acetaldehyde") > schema.get(cold, "acetaldehyde") > 0.0


def test_oxidation_is_first_order_in_oxygen(params):
    # First-order in the O₂ pool (the D-71 crux — O₂, not ethanol, is the rate-limiter): twice the
    # dissolved O₂ ⇒ twice the instantaneous oxidation rate. This linearity is what makes the pool
    # SATURATE (the rate falls as O₂ is spent), unlike a constant ethanol-first rate.
    schema = wine_schema()
    lo = OxidativeAcetaldehyde().derivatives(0.0, _aged_wine(schema, o2=0.02), schema, params)
    hi = OxidativeAcetaldehyde().derivatives(0.0, _aged_wine(schema, o2=0.04), schema, params)
    assert schema.get(hi, "acetaldehyde") == pytest.approx(2.0 * schema.get(lo, "acetaldehyde"))


def test_integrated_oxidation_saturates_and_closes_carbon(params, store):
    # Run a long aging segment (racked, dry wine — X=0, S=0) with ONLY OxidativeAcetaldehyde and a
    # dosed O₂ charge. Over the span the O₂ is consumed (depletes toward 0), acetaldehyde rises to a
    # PLATEAU (saturating, not unbounded — the whole point of the O₂-limited form), and total_carbon
    # closes to machine precision (E → acetaldehyde, O₂ off the ledger).
    schema = wine_schema()
    ps = ProcessSet(schema, [OxidativeAcetaldehyde()], strict=True)
    o2_0 = 0.04  # ~40 mg/L cumulative O₂ exposure
    y0 = _aged_wine(schema, ester=0.0, t=298.15, o2=o2_0)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))  # ~1 year
    assert traj.success, traj.message

    o2_end = float(traj.series("o2")[-1])
    acet_end = float(traj.series("acetaldehyde")[-1])
    # The O₂ charge is largely consumed and acetaldehyde has risen from nothing.
    assert o2_end < 0.1 * o2_0
    assert acet_end > 0.0
    # Saturating bound: acetaldehyde cannot exceed the yield-limited ceiling y·(o2_0/M_O2)·M_acet.
    ceiling = params["y_acetaldehyde_per_o2"] * (o2_0 / M_O2) * M_ACETALDEHYDE
    assert acet_end <= ceiling + 1e-12
    # And it lands in the oxidised-wine ballpark (tens–hundreds of mg/L) — a sanity anchor, not a
    # pinned magnitude (the yield is speculative): ≥ half the ceiling once O₂ is ~spent.
    assert acet_end >= 0.5 * ceiling
    assert_nonnegative(traj, ("o2", "acetaldehyde"), atol=1e-12)
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_oxidation_tier_floored_at_speculative(store):
    # Speculative in FORM (Tier-3 frontier), so the pools it writes are speculative even before the
    # (speculative) aging parameter tiers cap them. Non-vacuous: o2/acetaldehyde/E all speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [OxidativeAcetaldehyde()])
    for pool in ("o2", "acetaldehyde", "E"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# SulfiteOxidation (decision D-72) — the first oxidative sub-axis SINK to claim its share of the
# shared ``o2`` budget (D-71): dissolved O₂ oxidises free BISULFITE (HSO₃⁻ — the antioxidant
# nucleophile, NOT molecular SO₂ the antimicrobial form) → sulfate, consuming protective SO₂ at the
# Danilewicz 2:1 mol SO₂:O₂ stoichiometry. Bilinear in [o2]·[HSO₃⁻], Arrhenius warmer-faster. Both
# ``o2`` and ``so2_total`` are off every ledger (no sulfur ledger), so nothing conserved moves.
# WINE-ONLY (so2_total + the acid-pH slots are wine-only, D-18). These tests pin the closed form and
# the 2:1 stoichiometry, prove the SO₂-protection THRESHOLD (with the sibling OxidativeAcetaldehyde,
# SO₂ suppresses oxidative acetaldehyde until it is spent, then acetaldehyde climbs), the
# double-substrate isolability (inert without O₂ or without SO₂), the wine-only no-op on beer, the
# warmer-faster ordering, and the speculative tier floor.


@pytest.fixture
def so2_store():
    # Wine + the acidbase/acetaldehyde/keto-acid pKa + SO₂-binding params SulfiteOxidation's
    # pH/bisulfite readout reads (D-72), merged with aging.yaml — the shared_files set the D-72
    # compile seam wires. (The plain ``store`` fixture omits them: EsterHydrolysis/
    # OxidativeAcetaldehyde never solve pH, but this Process does.)
    d = default_data_dir()
    return load_parameters(
        d / "wine_generic.yaml",
        d / "acidbase.yaml",
        d / "acetaldehyde.yaml",
        d / "keto_acids.yaml",
        d / "aging.yaml",
    )


@pytest.fixture
def so2_params(so2_store):
    return so2_store.resolve()


def _sulfited_wine(
    schema: StateSchema, *, so2: float = 0.03, o2: float = 0.03, t: float = 293.15, **kw
) -> FloatArray:
    """A finished, racked wine at the start of aging with a real acid load (so pH solves into the
    wine range), dosed SO₂ and O₂. ``tartaric`` + a ``cation_charge`` set an acidic ~pH 3.3; the
    dosed ``so2_total``/``o2`` are the substrates. ``acetaldehyde`` defaults to 0 so free SO₂ ==
    total (no binding) and the bisulfite driver is unambiguous unless a test sets otherwise."""
    y = _aged_wine(schema, ester=0.0, t=t, so2_total=so2, o2=o2, tartaric=4.0, cation_charge=0.012)
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def test_sulfite_oxidation_metadata():
    p = SulfiteOxidation()
    assert p.name == "sulfite_oxidation"
    assert p.tier is Tier.SPECULATIVE
    # Consumes the O₂ substrate and oxidises the free-bisulfite share of so2_total — both off every
    # ledger, so nothing conserved moves. Touches those two and nothing else.
    assert set(p.touches) == {"o2", "so2_total"}
    # Only its OWN aging params + shared T_ref; the plausible pKa/binding params read via acidbase
    # are omitted (Process already speculative — the MalolacticConversion/brett convention).
    assert set(p.reads) == {"k_so2_oxidation", "E_a_so2_oxidation", "T_ref"}


def test_sulfite_oxidation_matches_closed_form(so2_params):
    schema = wine_schema()
    so2, o2, t = 0.03, 0.03, 298.15  # off T_ref so the Arrhenius factor bites
    y = _sulfited_wine(schema, so2=so2, o2=o2, t=t)
    d = SulfiteOxidation().derivatives(0.0, y, schema, so2_params)

    ph = ph_of_state(y, schema, so2_params)
    bisulfite = bisulfite_so2_at_ph(y, schema, so2_params, ph)  # the reactive HSO₃⁻ driver, g/L
    f_t = arrhenius_factor(t, so2_params["E_a_so2_oxidation"], so2_params["T_ref"])
    r_o2 = so2_params["k_so2_oxidation"] * f_t * o2 * bisulfite  # bilinear g O₂/L/h

    assert bisulfite > 0.0  # the driver is live (guards against a vacuous pass)
    assert schema.get(d, "o2") == pytest.approx(-r_o2)
    assert schema.get(d, "so2_total") == pytest.approx(-_SO2_PER_O2 * (r_o2 / M_O2) * M_SO2)
    # Touches nothing else — not ethanol/acetaldehyde/esters, no sugar, no CO2, no biomass.
    for var in (
        "X",
        "S",
        "E",
        "N",
        "CO2",
        "acetaldehyde",
        "isoamyl_acetate",
        "isoamyl_alcohol",
        "Byp",
    ):
        assert schema.get(d, var) == 0.0


def test_sulfite_oxidation_two_to_one_stoichiometry(so2_params):
    # 2 mol SO₂ oxidised per mol O₂ consumed (the Danilewicz quinone-reduction + peroxide-scavenging
    # mechanism), which is exactly the classic ~4:1 SO₂:O₂ MASS rule (2·M_SO2/M_O2 = 4). A pure code
    # constant, not a fitted parameter — verified in both molar and mass form.
    schema = wine_schema()
    d = SulfiteOxidation().derivatives(0.0, _sulfited_wine(schema), schema, so2_params)
    moles_o2 = schema.get(d, "o2") / M_O2
    moles_so2 = schema.get(d, "so2_total") / M_SO2
    assert moles_so2 / moles_o2 == pytest.approx(_SO2_PER_O2)  # 2 mol SO₂ per mol O₂
    assert schema.get(d, "so2_total") / schema.get(d, "o2") == pytest.approx(2.0 * M_SO2 / M_O2)
    assert schema.get(d, "so2_total") / schema.get(d, "o2") == pytest.approx(4.0, rel=1e-3)


def test_sulfite_oxidation_bilinear_in_both_substrates(so2_params):
    # Bilinear: with acetaldehyde = 0 (free SO₂ == total, no binding), doubling either O₂ or SO₂
    # doubles the instantaneous rate. This is what makes SO₂ out-compete ethanol for O₂ in
    # proportion to how much free SO₂ remains — the mechanism behind the depletion threshold.
    schema = wine_schema()
    base = SulfiteOxidation().derivatives(
        0.0, _sulfited_wine(schema, so2=0.03, o2=0.03), schema, so2_params
    )
    dbl_o2 = SulfiteOxidation().derivatives(
        0.0, _sulfited_wine(schema, so2=0.03, o2=0.06), schema, so2_params
    )
    dbl_so2 = SulfiteOxidation().derivatives(
        0.0, _sulfited_wine(schema, so2=0.06, o2=0.03), schema, so2_params
    )
    assert schema.get(dbl_o2, "o2") == pytest.approx(2.0 * schema.get(base, "o2"))
    # Doubling total SO₂ doubles free bisulfite (acetaldehyde = 0 ⇒ free = total), hence the rate.
    assert schema.get(dbl_so2, "o2") == pytest.approx(2.0 * schema.get(base, "o2"))


def test_sulfite_oxidation_inert_without_oxygen(so2_params):
    # No oxidant ⇒ no scavenging (and no wasted pH solve): a reductive begin_aging (no add_oxygen)
    # is byte-for-byte the case without this Process. Also the o2<0 solver-undershoot guard.
    schema = wine_schema()
    ps = ProcessSet(schema, [SulfiteOxidation()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _sulfited_wine(schema, o2=0.0), so2_params), schema.zeros()
    )
    assert np.array_equal(
        SulfiteOxidation().derivatives(0.0, _sulfited_wine(schema, o2=-1e-6), schema, so2_params),
        schema.zeros(),
    )


def test_sulfite_oxidation_inert_without_so2(so2_params):
    # No SO₂ ⇒ no scavenging: an unsulfited aging is byte-for-byte the case without this Process
    # (only OxidativeAcetaldehyde then acts on the o2 pool). Also the so2<0 undershoot guard.
    schema = wine_schema()
    ps = ProcessSet(schema, [SulfiteOxidation()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _sulfited_wine(schema, so2=0.0), so2_params), schema.zeros()
    )
    assert np.array_equal(
        SulfiteOxidation().derivatives(0.0, _sulfited_wine(schema, so2=-1e-6), schema, so2_params),
        schema.zeros(),
    )


def test_sulfite_oxidation_rises_with_temperature(so2_params):
    # The sourced ordering (E_a_so2_oxidation > 0): warmer storage oxidises the protective SO₂
    # faster — warm cellars burn through SO₂ (and lose oxidative protection) faster.
    schema = wine_schema()
    cold = SulfiteOxidation().derivatives(0.0, _sulfited_wine(schema, t=283.15), schema, so2_params)
    warm = SulfiteOxidation().derivatives(0.0, _sulfited_wine(schema, t=303.15), schema, so2_params)
    assert schema.get(warm, "so2_total") < schema.get(cold, "so2_total") < 0.0
    assert schema.get(warm, "o2") < schema.get(cold, "o2") < 0.0


def test_sulfite_oxidation_is_wine_only_noop_on_beer(so2_params):
    # WINE-ONLY: so2_total + the acid-pH slots are wine-only (beer's pH/SO₂ system deferred, D-18),
    # so the SO2_STATE_KEY-absent guard makes this a hard no-op on beer even if it were wired there.
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    yb[beer.slice("o2")] = 0.03  # o2 exists in both media, but so2_total does not
    assert np.array_equal(SulfiteOxidation().derivatives(0.0, yb, beer, so2_params), beer.zeros())


def test_sulfite_oxidation_protects_until_exhausted(so2_params, so2_store):
    # THE HEADLINE (D-72): run the two oxidative Processes together over a ~1-year aging segment
    # with a fixed O₂ charge. With no SO₂ the O₂ makes ~full oxidative acetaldehyde; SO₂ competes
    # for that O₂ and SUPPRESSES acetaldehyde in proportion to the dose, and enough SO₂ (> ~4× O₂)
    # nearly fully protects — the celebrated "SO₂ protects until exhausted" threshold, emergent from
    # the two Processes summing over the shared o2 pool. Off-ledger o2/so2, so carbon still closes.
    schema = wine_schema()
    o2_0 = 0.04  # ~40 mg/L O₂ charge

    def run(so2_0: float) -> Trajectory:
        ps = ProcessSet(schema, [OxidativeAcetaldehyde(), SulfiteOxidation()], strict=True)
        y0 = _sulfited_wine(schema, so2=so2_0, o2=o2_0, t=298.15)
        traj = simulate(ps, params=so2_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
        assert traj.success, traj.message
        return traj

    none = run(0.0)
    modest = run(0.05)  # ~50 mg/L: partial protection, SO₂ largely consumed
    ample = run(0.30)  # ~300 mg/L (> 4·40): near-full protection, SO₂ left over

    acet_none = float(none.series("acetaldehyde")[-1])
    acet_modest = float(modest.series("acetaldehyde")[-1])
    acet_ample = float(ample.series("acetaldehyde")[-1])

    # Monotone protection: more SO₂ ⇒ less oxidative acetaldehyde.
    assert acet_none > acet_modest > acet_ample
    # Ample SO₂ nearly abolishes the oxidative note (well under a tenth of the unprotected level).
    assert acet_ample < 0.1 * acet_none
    # SO₂ is genuinely consumed where it is the limiting reactant (modest dose ~spent), and remains
    # where it is in excess (ample dose retains a protective reserve) — the depletion threshold.
    assert float(modest.series("so2_total")[-1]) < 0.5 * 0.05
    assert float(ample.series("so2_total")[-1]) > 0.30 - 4.1 * o2_0  # ≥ dose − stoichiometric burn
    # Non-negative pools; carbon closes (E → acetaldehyde is the only on-ledger move; o2/so2 off).
    assert_nonnegative(ample, ("o2", "so2_total", "acetaldehyde"), atol=1e-9)
    f_c = so2_store.value("biomass_C_fraction")
    assert_conserved(ample, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_sulfite_oxidation_tier_floored_at_speculative(so2_store):
    # Speculative in FORM (Tier-3 frontier): the pools it writes are speculative even before the
    # (speculative) aging parameter tiers cap them. Non-vacuous: o2/so2_total both speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [SulfiteOxidation()])
    for pool in ("o2", "so2_total"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, so2_store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# PhenolicBrowning (decision D-74) — the first ALWAYS-ON sink on the shared ``o2`` budget (D-71):
# dissolved O₂ oxidises phenolics → brown pigment, accumulating the ``A420`` browning index (an
# optical absorbance, dimensionless AU — NOT a mass). First-order in [o2] (its OWN, DOMINANT share
# ``k_browning_base`` > ``k_ethanol_oxidation``), Arrhenius warmer-faster. Touches only ``o2`` +
# ``A420``,
# BOTH off every ledger — so it moves NOTHING conserved (the cleanest aging Process; not even a
# carbon borrow). MEDIUM-AGNOSTIC (both media brown; ``A420`` exists in both schemas). These tests
# pin
# the closed form, the first-order-in-O₂ linearity, the monotonic A420 accumulation + saturation,
# the
# medium-agnostic run on beer, the reductive-aging isolability (inert without O₂), the
# off-every-ledger
# invariance (carbon AND mass both flat), the headline O₂-diversion (browning suppresses oxidative
# acetaldehyde — the always-on analogue of the D-72 SO₂ threshold), the warmer-faster ordering, and
# the speculative tier floor.


def test_browning_metadata():
    p = PhenolicBrowning()
    assert p.name == "phenolic_browning"
    assert p.tier is Tier.SPECULATIVE
    # Consumes its O₂ share and books the oxidised phenol as the A420 browning index — both off
    # every
    # ledger, so nothing conserved moves. Touches those two and nothing else (not even a carbon
    # borrow). tannin/anthocyanin are READ (D-132) but never written, so — like T elsewhere — they
    # are not part of touches.
    assert set(p.touches) == {"o2", "A420"}
    assert set(p.reads) == {
        "k_browning_base",
        "k_browning_phenolic",
        "E_a_browning",
        "y_a420_per_o2",
        "T_ref",
    }


def test_browning_matches_closed_form(params):
    schema = wine_schema()
    o2, t = 0.03, 298.15  # off T_ref so the Arrhenius factor bites
    y = _aged_wine(schema, ester=0.0, t=t, o2=o2)
    d = PhenolicBrowning().derivatives(0.0, y, schema, params)
    f_t = arrhenius_factor(t, params["E_a_browning"], params["T_ref"])
    # No tannin/anthocyanin dosed ⇒ k_browning_eff is byte-for-byte k_browning_base (D-132
    # isolability at zero grape phenolics).
    r_o2 = params["k_browning_base"] * f_t * o2
    assert schema.get(d, "o2") == pytest.approx(-r_o2)
    assert schema.get(d, "A420") == pytest.approx(params["y_a420_per_o2"] * (r_o2 / M_O2))
    # Touches ONLY o2 + A420 — nothing else moves (not even E/acetaldehyde: browning borrows no
    # carbon).
    for var in schema.names:
        if var not in ("o2", "A420"):
            assert schema.get(d, var) == 0.0


def test_browning_is_dominant_share_over_ethanol_oxidation(params):
    # The load-bearing D-74 ordering: browning is the DOMINANT always-on O₂ sink, so at the same
    # [o2]
    # it draws a larger O₂ rate than ethanol oxidation (k_browning_base > k_ethanol_oxidation), and
    # the
    # two BASELINE shares sum to the calibrated always-on total (5.0e-4) that holds the
    # O₂-depletion
    # timescale (D-132 leaves this floor untouched; a real red's effective rate is higher still).
    assert params["k_browning_base"] > params["k_ethanol_oxidation"]
    assert params["k_browning_base"] + params["k_ethanol_oxidation"] == pytest.approx(5.0e-4)
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, o2=0.03)
    brown = PhenolicBrowning().derivatives(0.0, y, schema, params)
    ethanol = OxidativeAcetaldehyde().derivatives(0.0, y, schema, params)
    assert -schema.get(brown, "o2") > -schema.get(ethanol, "o2") > 0.0


def test_browning_phenolic_boost_matches_closed_form(params):
    # D-132: dosed tannin + anthocyanin lift k_browning_eff above the baseline. Closed form:
    # k_browning_eff = k_browning_base + k_browning_phenolic * (tannin + anthocyanin).
    schema = wine_schema()
    o2, t = 0.03, 298.15
    tannin, anthocyanin = 2.0, 0.3  # the "typical red" anchor (polymerization.yaml D-79/D-81/D-84)
    y = _aged_wine(schema, ester=0.0, t=t, o2=o2, tannin=tannin, anthocyanin=anthocyanin)
    d = PhenolicBrowning().derivatives(0.0, y, schema, params)
    f_t = arrhenius_factor(t, params["E_a_browning"], params["T_ref"])
    k_eff = params["k_browning_base"] + params["k_browning_phenolic"] * (tannin + anthocyanin)
    assert k_eff > params["k_browning_base"]  # the boost strictly raises the rate
    r_o2 = k_eff * f_t * o2
    assert schema.get(d, "o2") == pytest.approx(-r_o2)
    assert schema.get(d, "A420") == pytest.approx(params["y_a420_per_o2"] * (r_o2 / M_O2))


def test_browning_phenolic_boost_isolable_at_zero_phenolics(params):
    # D-132 isolability (the D-129/D-131 GATE-1 pattern): explicitly dosing zero tannin AND zero
    # anthocyanin is byte-for-byte the case without the phenolic term at all.
    schema = wine_schema()
    o2, t = 0.03, 298.15
    zero = PhenolicBrowning().derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=t, o2=o2, tannin=0.0, anthocyanin=0.0), schema, params
    )
    baseline = PhenolicBrowning().derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=t, o2=o2), schema, params
    )
    assert np.array_equal(zero, baseline)


def test_browning_phenolic_boost_rises_with_phenolic_load(params):
    # Monotone in combined grape phenolics: more tannin+anthocyanin ⇒ faster O₂ uptake / A420 rise
    # (the D-132 sourced ordering — more phenolics, more oxidation).
    schema = wine_schema()
    o2, t = 0.03, 298.15
    low = PhenolicBrowning().derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=t, o2=o2, tannin=0.5, anthocyanin=0.05), schema, params
    )
    high = PhenolicBrowning().derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=t, o2=o2, tannin=2.0, anthocyanin=0.3), schema, params
    )
    assert -schema.get(high, "o2") > -schema.get(low, "o2") > 0.0
    assert schema.get(high, "A420") > schema.get(low, "A420") > 0.0


def test_browning_phenolic_boost_absent_on_beer(params):
    # tannin/anthocyanin are wine-only must-input slots (D-79), absent from beer's schema. The
    # Process guards their absence rather than gating the whole run (unlike a wine-only Process):
    # beer still browns, at the unboosted k_browning_base rate.
    beer = beer_schema()
    assert "tannin" not in beer and "anthocyanin" not in beer
    o2, t = 0.03, 298.15
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": t, "CO2": 0.0})
    yb[beer.slice("o2")] = o2
    d = PhenolicBrowning().derivatives(0.0, yb, beer, params)
    f_t = arrhenius_factor(t, params["E_a_browning"], params["T_ref"])
    r_o2 = params["k_browning_base"] * f_t * o2
    assert beer.get(d, "o2") == pytest.approx(-r_o2)


def test_browning_typical_red_lands_in_ferreira_band(params):
    # THE D-132 HEADLINE CALIBRATION: at a fresh ~8 mg/L O₂ saturation charge and a typical red's
    # grape phenolic load, the TOTAL O₂-depletion rate (ethanol oxidation + browning) lands in
    # Ferreira 2015's measured real-wine average of 0.5-0.7 mg/L/day — the previous medium-agnostic
    # rate under-predicted this by ~6-8x (see aging.yaml's k_browning_phenolic provenance).
    schema = wine_schema()
    o2 = 0.008  # ~8 mg/L, T_ref (f_t = 1, no Arrhenius correction needed)
    tannin, anthocyanin = 2.0, 0.3
    y = _aged_wine(schema, ester=0.0, o2=o2, tannin=tannin, anthocyanin=anthocyanin)
    brown = PhenolicBrowning().derivatives(0.0, y, schema, params)
    ethanol = OxidativeAcetaldehyde().derivatives(0.0, y, schema, params)
    total_rate_g_l_h = -schema.get(brown, "o2") - schema.get(ethanol, "o2")
    total_rate_mg_l_day = total_rate_g_l_h * 1000.0 * 24.0
    assert 0.5 <= total_rate_mg_l_day <= 0.7


def test_browning_is_first_order_in_oxygen(params):
    # First-order in the O₂ pool (its own share of the shared substrate): twice the dissolved O₂ ⇒
    # twice the instantaneous browning rate. This linearity is what makes A420 SATURATE as O₂ is
    # spent.
    schema = wine_schema()
    lo = PhenolicBrowning().derivatives(0.0, _aged_wine(schema, o2=0.02), schema, params)
    hi = PhenolicBrowning().derivatives(0.0, _aged_wine(schema, o2=0.04), schema, params)
    assert schema.get(hi, "A420") == pytest.approx(2.0 * schema.get(lo, "A420"))


def test_browning_inert_without_oxygen(params):
    # No oxidant ⇒ no browning: a reductive/un-oxygenated aging is byte-for-byte the case without
    # this
    # Process (A420 stays 0). The o2 ≤ 0 guard also absorbs a solver undershoot (o2 < 0 ⇒ no
    # spurious
    # browning), keeping d(A420)/dt ≥ 0 (A420 monotonic, never reversed).
    schema = wine_schema()
    assert np.array_equal(
        PhenolicBrowning().derivatives(0.0, _aged_wine(schema, o2=0.0), schema, params),
        schema.zeros(),
    )
    assert np.array_equal(
        PhenolicBrowning().derivatives(0.0, _aged_wine(schema, o2=-1e-6), schema, params),
        schema.zeros(),
    )


def test_browning_rises_with_temperature(params):
    # The sourced ordering (E_a_browning > 0): warmer storage browns (maderises) faster — more O₂
    # consumed and more A420 built per hour when warm.
    schema = wine_schema()
    cold = PhenolicBrowning().derivatives(
        0.0, _aged_wine(schema, o2=0.03, t=283.15), schema, params
    )
    warm = PhenolicBrowning().derivatives(
        0.0, _aged_wine(schema, o2=0.03, t=303.15), schema, params
    )
    assert schema.get(warm, "o2") < schema.get(cold, "o2") < 0.0
    assert schema.get(warm, "A420") > schema.get(cold, "A420") > 0.0


def test_browning_is_medium_agnostic_on_beer(params):
    # MEDIUM-AGNOSTIC (D-74, superseding D-73's provisional "wine-only"): both media carry
    # autoxidising
    # polyphenols and brown, and A420 exists in both schemas — so browning runs on beer too, at the
    # same closed form (o2/A420/T are shared slots; the shared aging.yaml params apply to both
    # media).
    beer = beer_schema()
    o2, t = 0.03, 298.15
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": t, "CO2": 0.0})
    yb[beer.slice("o2")] = o2
    d = PhenolicBrowning().derivatives(0.0, yb, beer, params)
    f_t = arrhenius_factor(t, params["E_a_browning"], params["T_ref"])
    r_o2 = params["k_browning_base"] * f_t * o2
    assert beer.get(d, "o2") == pytest.approx(-r_o2)
    assert beer.get(d, "A420") == pytest.approx(params["y_a420_per_o2"] * (r_o2 / M_O2))


def test_integrated_browning_accumulates_a420_and_saturates(params, store):
    # Run a long aging segment (racked, dry wine) with ONLY PhenolicBrowning and a dosed O₂ charge.
    # Over the span the O₂ is consumed (depletes toward 0), A420 rises MONOTONICALLY to a PLATEAU
    # (saturating as the O₂ charge is spent — the browning of an aged wine), and — because o2 + A420
    # are BOTH off every ledger — this Process moves NOTHING conserved: total_carbon AND total_mass
    # are both exactly flat (unusually strong: not even the E→acetaldehyde borrow
    # OxidativeAcetaldehyde
    # carries). A420's ceiling is y_a420_per_o2·(o2_0/M_O2) (every mol O₂ builds that much
    # absorbance).
    schema = wine_schema()
    ps = ProcessSet(schema, [PhenolicBrowning()], strict=True)
    o2_0 = 0.04  # ~40 mg/L cumulative O₂ exposure
    y0 = _aged_wine(schema, ester=0.0, t=298.15, o2=o2_0)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))  # ~1 year
    assert traj.success, traj.message

    a420 = traj.series("A420")
    o2_end = float(traj.series("o2")[-1])
    a420_end = float(a420[-1])
    # A420 rose from 0 and the O₂ charge is largely spent.
    assert o2_end < 0.1 * o2_0
    assert a420_end > 0.0
    # Monotone accumulation (d(A420)/dt ≥ 0 everywhere — no reversal).
    assert np.all(np.diff(np.asarray(a420, dtype=float)) >= -1e-15)
    # Saturating bound: A420 cannot exceed the yield-limited ceiling y·(o2_0/M_O2), and once O₂ is
    # ~spent it has reached most of it (a browned wine, a sanity anchor not a pinned magnitude).
    ceiling = params["y_a420_per_o2"] * (o2_0 / M_O2)
    assert a420_end <= ceiling + 1e-12
    assert a420_end >= 0.5 * ceiling
    assert_nonnegative(traj, ("o2", "A420"), atol=1e-12)
    # BOTH ledgers flat — browning moves nothing conserved (the cleanest aging Process).
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")


def test_browning_diverts_o2_and_suppresses_acetaldehyde(params, store):
    # THE HEADLINE (D-74): browning is a co-resident always-on O₂ sink, so putting it alongside
    # OxidativeAcetaldehyde over the same fixed O₂ charge DIVERTS most of the O₂ to brown pigment
    # and
    # correspondingly SUPPRESSES the oxidative acetaldehyde — the always-on analogue of the D-72 SO₂
    # threshold (but permanent: browning is not spent). Same k_ethanol_oxidation in both runs; the
    # only difference is whether browning competes for the O₂.
    schema = wine_schema()
    o2_0 = 0.04  # ~40 mg/L O₂ charge
    # A LONG span so BOTH runs fully plateau (the O₂ charge is spent): the ethanol-only run depletes
    # O₂ at the slower k_ethanol_oxidation alone, so it needs the long tail to reach its ceiling —
    # only once both have plateaued does the clean partition ratio
    # k_ethanol/(k_ethanol+k_browning_base) hold (in finite time the slower run lags its ceiling and
    # the ratio reads high).
    span = (0.0, 24.0 * 365.0 * 5.0)

    def run(processes: list[Process]) -> Trajectory:
        ps = ProcessSet(schema, processes, strict=True)
        y0 = _aged_wine(schema, ester=0.0, t=298.15, o2=o2_0)
        traj = simulate(ps, params=params, y0=y0, t_span=span)
        assert traj.success, traj.message
        return traj

    ethanol_only = run([OxidativeAcetaldehyde()])
    with_browning = run([OxidativeAcetaldehyde(), PhenolicBrowning()])

    acet_alone = float(ethanol_only.series("acetaldehyde")[-1])
    acet_diverted = float(with_browning.series("acetaldehyde")[-1])
    # Browning diverts O₂ ⇒ LESS oxidative acetaldehyde, and it builds visible brown (A420 > 0)
    # where
    # the ethanol-only run browns none. The suppression tracks the share: browning takes
    # ~k_browning_base /
    # (k_browning_base + k_ethanol) of the O₂, so acetaldehyde falls toward the ethanol share ~40%
    # (no tannin/anthocyanin dosed here ⇒ k_browning_eff is byte-for-byte k_browning_base, D-132).
    assert acet_diverted < acet_alone
    assert float(with_browning.series("A420")[-1]) > 0.0
    assert float(ethanol_only.series("A420")[-1]) == 0.0
    share_ethanol = params["k_ethanol_oxidation"] / (
        params["k_ethanol_oxidation"] + params["k_browning_base"]
    )
    assert acet_diverted == pytest.approx(share_ethanol * acet_alone, rel=0.05)
    # Carbon still closes (E → acetaldehyde the only on-ledger move; o2/A420 off every ledger).
    f_c = store.value("biomass_C_fraction")
    assert_conserved(
        with_browning, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon"
    )


def test_browning_tier_floored_at_speculative(store):
    # Speculative in FORM (Tier-3 frontier): the pools it writes are speculative even before the
    # (speculative) aging parameter tiers cap them. Non-vacuous: o2/A420 both speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [PhenolicBrowning()])
    for pool in ("o2", "A420"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# StreckerDegradation (decision D-75) — the WINE-ONLY, DOUBLY substrate-gated oxidative aging sink
# that produces the Strecker aldehydes methional (cooked-potato off-note) + phenylacetaldehyde
# (honey). Dissolved O2 — via the phenol-oxidation quinones — degrades amino acids: the carbon is
# drawn from ``amino_acids`` (arginine stand-in), the nitrogen deaminated to ``N``, one CO2 released
# per aldehyde (the D-45 mercaptan idiom + a decarboxylation term). Gated on BOTH ``o2`` AND
# ``amino_acids`` (the O2 draw itself carries the availability gate), so — like SulfiteOxidation —
# it ADDS ON TOP of the shared o2 budget WITHOUT re-baselining (superseding the D-71..D-74
# "reduce k_ethanol again" forward-guess). These tests pin the closed form, carbon AND nitrogen
# closure per-RHS, the double-substrate isolability (inert without O2 or without amino acids), the
# first-order-in-O2 linearity + aa-availability gate, the methional:phenylacetaldehyde split, the
# wine-only no-op on beer, the warmer-faster ordering, the integrated saturation + closure, and the
# speculative tier floor (incl. the structural N-write).

#: Each aldehyde now draws ITS OWN precursor (D-100): methional really is methionine minus its
#: carboxyl, phenylacetaldehyde really is phenylalanine minus its. The lumped arginine draw —
#: which booked a molecule containing no sulfur as methional's source — is retired.
_STRECKER_TOUCHES = {
    "o2",
    "methional",
    "phenylacetaldehyde",
    "CO2",
    "N",
    "methionine",
    "phenylalanine",
}
#: (precursor, product) — the D-100 map the two Strecker routes and the thermal route share.
_STRECKER_PRECURSORS = {"methional": "methionine", "phenylacetaldehyde": "phenylalanine"}


def _strecker_wine(
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    aa: float = 0.05,
    o2: float = 0.03,
    t: float = 293.15,
    **kw,
) -> FloatArray:
    """A finished, racked wine at the start of aging with dosed amino acids + O2 — the two Strecker
    substrates. Amino acids are seeded at MUST-SPECTRUM composition (D-100), the state in which
    every per-precursor gate provably equals the pre-split lumped gate, so the closed form below
    asserts the same numbers the lumped suite did. ``esters`` defaults to 0 (irrelevant here)."""
    y = _aged_wine(schema, ester=0.0, t=t, o2=o2)
    seed_amino_acids(y, schema, params, aa)
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def _strecker_closed_form(
    schema: StateSchema, params: dict[str, float], y: FloatArray, t: float
) -> dict[str, float]:
    """The Process's own algebra, recomputed independently for the closed-form assertions."""
    o2 = float(y[schema.slice("o2")][0])
    # At must-spectrum composition every per-precursor gate aa_i/(K·f_i + aa_i) collapses to this
    # one lumped value (the D-100 reduction property), so the pre-split algebra is EXACT here —
    # which is the property this closed form now doubles as a proof of.
    aa = sum(float(y[schema.slice(spec.pool)][0]) for spec in AMINO_ACID_SPECS)
    gate = aa / (params["K_amino_acids"] + aa)
    f_t = arrhenius_factor(t, params["E_a_strecker"], params["T_ref"])
    r_o2 = params["k_strecker"] * f_t * o2 * gate
    n_ald = params["y_strecker_per_o2"] * (r_o2 / M_O2)
    f_meth = params["f_methional"]
    return {
        "o2": -r_o2,
        "methional": f_meth * n_ald * M_METHIONAL,
        "phenylacetaldehyde": (1.0 - f_meth) * n_ald * M_PHENYLACETALDEHYDE,
        "CO2": n_ald * M_CO2,
    }


def test_strecker_metadata():
    p = StreckerDegradation()
    assert p.name == "strecker_degradation"
    # Speculative: the aging axis is the Tier-3 frontier (form sourced, magnitudes estimated).
    assert p.tier is Tier.SPECULATIVE
    # Consumes its aa-gated O2 share; books the two aldehydes + the decarboxylation CO2, drawing
    # carbon from amino_acids and deaminating the nitrogen to N. Touches those six and nothing else.
    assert set(p.touches) == _STRECKER_TOUCHES
    assert "amino_acids" not in p.touches  # the lumped draw is RETIRED (D-100)
    assert set(p.reads) == {
        "k_strecker",
        "E_a_strecker",
        "y_strecker_per_o2",
        "f_methional",
        "K_amino_acids",
        "T_ref",
        "must_aa_fraction_methionine",
        "must_aa_fraction_phenylalanine",
    }


def test_strecker_matches_closed_form(params):
    schema = wine_schema()
    aa, o2, t = 0.05, 0.03, 298.15  # off T_ref so the Arrhenius factor bites
    y = _strecker_wine(schema, params, aa=aa, o2=o2, t=t)
    d = StreckerDegradation().derivatives(0.0, y, schema, params)
    cf = _strecker_closed_form(schema, params, y, t)

    assert cf["methional"] > 0.0  # the products are live (guards against a vacuous pass)
    assert schema.get(d, "o2") == pytest.approx(cf["o2"])
    assert schema.get(d, "methional") == pytest.approx(cf["methional"])
    assert schema.get(d, "phenylacetaldehyde") == pytest.approx(cf["phenylacetaldehyde"])
    assert schema.get(d, "CO2") == pytest.approx(cf["CO2"])
    # EACH precursor is drawn sized to the carbon of the product IT made, plus the CO2 that
    # product's own carboxyl released (D-100). The CO2 attribution is load-bearing and invisible to
    # conservation — closure holds for any split — so it is asserted per-precursor here.
    n_ald = cf["CO2"] / M_CO2  # total mol aldehyde; 1 CO2 each (both routes decarboxylate)
    f_meth = params["f_methional"]
    per_precursor = {
        "methionine": cf["methional"] * _METHIONAL_C + f_meth * n_ald * M_CO2 * _CO2_C,
        "phenylalanine": (
            cf["phenylacetaldehyde"] * _PHENYLACET_C + (1.0 - f_meth) * n_ald * M_CO2 * _CO2_C
        ),
    }
    expected_n = 0.0
    for precursor, carbon in per_precursor.items():
        mass = carbon / carbon_mass_fraction(precursor)
        assert schema.get(d, precursor) == pytest.approx(-mass)
        expected_n += mass * nitrogen_mass_fraction(precursor)
    # ...and the deamination releases the nitrogen THOSE molecules carried, not arginine's.
    assert schema.get(d, "N") == pytest.approx(expected_n)
    assert schema.get(d, "amino_acids") == 0.0  # the lumped pool is untouched (D-100)
    # Touches nothing else — no ethanol/esters/higher alcohols/acetaldehyde, no sugar, no biomass.
    for var in ("X", "S", "E", "isoamyl_acetate", "isoamyl_alcohol", "Byp", "acetaldehyde"):
        assert schema.get(d, var) == 0.0


def test_strecker_carbon_closes_per_rhs(params):
    # CARBON closes to machine precision: the carbon leaving EACH precursor, at ITS OWN carbon
    # fraction (D-100), equals the carbon entering methional + phenylacetaldehyde + CO2 — a pure
    # on-ledger transfer, off-ledger o2 aside.
    schema = wine_schema()
    d = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, t=298.15), schema, params
    )
    carbon_residual = (
        schema.get(d, "methional") * _METHIONAL_C
        + schema.get(d, "phenylacetaldehyde") * _PHENYLACET_C
        + schema.get(d, "CO2") * _CO2_C
        + sum(schema.get(d, aa) * carbon_mass_fraction(aa) for aa in _STRECKER_PRECURSORS.values())
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)


def test_strecker_nitrogen_closes_per_rhs(params):
    # NITROGEN closes: all the nitrogen leaving the precursors, at THEIR OWN fractions (D-100),
    # lands in the N pool (the aldehydes are nitrogen-free — the deamination, the D-45 idiom).
    schema = wine_schema()
    d = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, t=298.15), schema, params
    )
    nitrogen_residual = (
        sum(schema.get(d, aa) * nitrogen_mass_fraction(aa) for aa in _STRECKER_PRECURSORS.values())
        + schema.get(d, "N") * 1.0
    )
    assert nitrogen_residual == pytest.approx(0.0, abs=1e-18)


def test_strecker_inert_without_oxygen(params):
    # No oxidant ⇒ no Strecker: a reductive begin_aging (no add_oxygen) is byte-for-byte the case
    # without this Process. Also the o2<0 solver-undershoot guard.
    schema = wine_schema()
    ps = ProcessSet(schema, [StreckerDegradation()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _strecker_wine(schema, params, o2=0.0), params), schema.zeros()
    )
    assert np.array_equal(
        StreckerDegradation().derivatives(
            0.0, _strecker_wine(schema, params, o2=-1e-6), schema, params
        ),
        schema.zeros(),
    )


def test_strecker_inert_without_amino_acids(params):
    # No amino acids ⇒ no Strecker: an amino-acid-free aging is byte-for-byte the case without this
    # Process (the substrate gate that makes it ADD ON TOP with no re-baseline, D-75). Also aa<0.
    schema = wine_schema()
    ps = ProcessSet(schema, [StreckerDegradation()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _strecker_wine(schema, params, aa=0.0), params), schema.zeros()
    )
    assert np.array_equal(
        StreckerDegradation().derivatives(
            0.0, _strecker_wine(schema, params, aa=-1e-6), schema, params
        ),
        schema.zeros(),
    )


def test_strecker_first_order_in_oxygen(params):
    # First-order in the O2 pool (at fixed amino acids, so the availability gate is held constant):
    # doubling [o2] doubles the instantaneous O2 draw and every product rate.
    schema = wine_schema()
    base = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, o2=0.02), schema, params
    )
    dbl = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, o2=0.04), schema, params
    )
    for pool in (
        "o2",
        "methional",
        "phenylacetaldehyde",
        "CO2",
        "methionine",
        "phenylalanine",
        "N",
    ):
        assert schema.get(dbl, pool) == pytest.approx(2.0 * schema.get(base, pool))


def test_strecker_availability_gate_saturates(params):
    # The amino-acid availability gate aa/(K+aa) throttles the draw at low aa and SATURATES toward a
    # ceiling at high aa (the smooth swap/reroute gate, D-33). At aa >> K the rate approaches the
    # ungated k*f*[o2]; at aa == K it is ~half that — a monotone, saturating aa dependence.
    schema = wine_schema()
    k = params["K_amino_acids"]
    low = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, aa=0.1 * k), schema, params
    )
    mid = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, aa=k), schema, params
    )
    high = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, aa=100.0 * k), schema, params
    )
    # Monotone increasing O2 draw magnitude with amino acids, but saturating (not linear).
    assert 0.0 < -schema.get(low, "o2") < -schema.get(mid, "o2") < -schema.get(high, "o2")
    # mid (aa = K) is ~half the high-aa ceiling (gate = 0.5 vs -> 1); low (aa = 0.1K) far below.
    assert -schema.get(mid, "o2") == pytest.approx(0.5 * -schema.get(high, "o2"), rel=0.02)


def test_strecker_split_phenylacetaldehyde_dominant(params):
    # The mol split between the two aldehydes is f_methional : (1 - f_methional); with the default
    # f_methional = 0.15 PHENYLACETALDEHYDE dominates (the split is production flux — phenylalanine
    # is far more abundant in must than methionine; potency lives in the OAV threshold, not here).
    schema = wine_schema()
    d = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, params), schema, params)
    meth_mol = schema.get(d, "methional") / M_METHIONAL
    phenyl_mol = schema.get(d, "phenylacetaldehyde") / M_PHENYLACETALDEHYDE
    f_meth = params["f_methional"]
    assert meth_mol / phenyl_mol == pytest.approx(f_meth / (1.0 - f_meth))
    assert phenyl_mol > meth_mol > 0.0  # phenylacetaldehyde-dominant (the more abundant precursor)
    # One CO2 per aldehyde (the decarboxylation) — total aldehyde mol equals CO2 mol.
    assert (meth_mol + phenyl_mol) == pytest.approx(schema.get(d, "CO2") / M_CO2)


def test_strecker_rises_with_temperature(params):
    # The sourced ordering (E_a_strecker > 0): warmer storage forms Strecker aldehydes faster — the
    # canonical warm-storage staling/oxidation direction.
    schema = wine_schema()
    cold = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, t=283.15), schema, params
    )
    warm = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, params, t=303.15), schema, params
    )
    assert schema.get(warm, "methional") > schema.get(cold, "methional") > 0.0
    assert schema.get(warm, "o2") < schema.get(cold, "o2") < 0.0  # more O2 drawn when warm


def test_strecker_is_wine_only_noop_on_beer(params):
    # WINE-ONLY: amino_acids + the N-deamination read wine-only slots (beer's amino-acid pool is not
    # tracked, D-32), so the "amino_acids not in schema" guard makes this a hard no-op on beer even
    # though o2 exists in both media.
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    yb[beer.slice("o2")] = 0.03
    assert np.array_equal(StreckerDegradation().derivatives(0.0, yb, beer, params), beer.zeros())


def test_integrated_strecker_saturates_and_closes(params, store):
    # Integrate Strecker ALONGSIDE the dominant always-on O2 sinks (OxidativeAcetaldehyde +
    # PhenolicBrowning — its real co-residents, which deplete the shared o2 charge on the ~weeks-
    # months timescale) over a warm ~1-year aging segment with a fixed O2 + amino-acid charge. The
    # Strecker aldehydes ACCUMULATE and SATURATE as the O2 is spent — a bounded, substrate-limited
    # climb, not the unbounded rise a rate first-order in ethanol would give. This also exercises
    # D-75 headline: Strecker ADDS ON TOP of the shared budget (its own aa-gated share) without
    # perturbing the sibling sinks. Carbon AND nitrogen both close to machine precision.
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [OxidativeAcetaldehyde(), PhenolicBrowning(), StreckerDegradation()],
        strict=True,
    )
    y0 = _strecker_wine(schema, params, aa=0.05, o2=0.04, t=298.15)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message

    meth = traj.series("methional")
    phenyl = traj.series("phenylacetaldehyde")
    # Both aldehydes accumulate (monotone, produced-only) and end well above zero.
    assert meth[-1] > meth[0] == 0.0
    assert phenyl[-1] > phenyl[0] == 0.0
    # Saturation (by TIME, robust to a non-uniform solver mesh): the second-half gain is a small
    # fraction of the first-half gain — the O2 charge is spent, so production has plateaued rather
    # than climbing linearly (the D-71 saturating-vs-unbounded distinction, inherited via the shared
    # o2 pool).
    mid = int(np.searchsorted(traj.t, 0.5 * traj.t[-1]))
    first_half = meth[mid] - meth[0]
    second_half = meth[-1] - meth[mid]
    assert 0.0 < second_half < 0.2 * first_half
    # Non-negative pools; carbon + nitrogen close (off-ledger o2/A420 aside; acetaldehyde from E).
    assert_nonnegative(
        traj, ("o2", "amino_acids", "methional", "phenylacetaldehyde", "N"), atol=1e-9
    )
    f_c = store.value("biomass_C_fraction")
    f_n = store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_strecker_tier_floored_at_speculative(store):
    # Speculative in FORM (Tier-3 frontier): every pool it writes is speculative even before the
    # (speculative) aging parameter tiers cap them. Non-vacuous across all six touched pools —
    # including the structural N-write (the first aging Process to write N, the D-45 note).
    schema = wine_schema()
    ps = ProcessSet(schema, [StreckerDegradation()])
    for pool in _STRECKER_TOUCHES:
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# MaillardStrecker (decision D-87) — the WINE-ONLY, NON-oxidative THERMAL Strecker route: the
# O₂-INDEPENDENT thermal mirror of StreckerDegradation (D-75). Residual SUGAR + HEAT (α-dicarbonyls,
# NO O₂) degrade amino acids to the sweet-wine / Madeira suite: methional + phenylacetaldehyde
# (SHARED with D-75) + three branched-chain malty aldehydes (2-/3-methylbutanal, 2-methylpropanal) +
# sotolon (the curry/maple furanone). Carbon from ``amino_acids`` (arginine), nitrogen deaminated to
# ``N``, 1 CO2 per DECARBOXYLATING aldehyde (the five Strecker aldehydes) but NONE for sotolon (a
# furanone). S is a read-only DRIVER (not consumed here); NO ``o2`` term. These tests pin the closed
# form, carbon AND nitrogen closure, the aa HARD gate (isolability) + the sugar SOFT gate, the
# O₂-INDEPENDENCE (identical with/without O₂ — the discriminating mirror of D-75), the first-order-
# in-sugar linearity + aa-availability saturation, the composition-split normalization + the
# sotolon-no-CO2 flag, the warmer-faster ordering, the wine-only no-op on beer, the integrated
# sealed-sweet accumulation + closure, the discriminating contrast vs the O₂-only D-75 route, and
# the
# speculative tier floor (incl. the structural N-write).

#: Each product draws ITS OWN precursor (D-100) — five distinct amino acids over six products.
#: sotolon<-threonine is the load-bearing entry: threonine is ALSO propanol's Ehrlich precursor, so
#: the propanol-vs-sotolon competition is real chemistry over one molecule (unlike the retired
#: fusels-vs-arginine competition, which was an artifact of the lump).
_MAILLARD_PRECURSORS = {
    "methional": "methionine",
    "phenylacetaldehyde": "phenylalanine",
    "2_methylbutanal": "isoleucine",
    "3_methylbutanal": "leucine",
    "2_methylpropanal": "valine",
}
# SOTOLON LEFT THIS MAP AT D-107. It was never a Strecker degradation of threonine: it is an aldol
# of alpha-ketobutyrate + acetaldehyde (Pham et al. 1995), and threonine is only its GRANDparent
# (threonine -> alpha-ketobutyrate -> sotolon). See SotolonAldolCondensation and its tests below.
_MAILLARD_TOUCHES = {
    "methional",
    "phenylacetaldehyde",
    "2_methylbutanal",
    "3_methylbutanal",
    "2_methylpropanal",
    "CO2",
    "N",
    # S is READ-ONLY again since D-107 (it was a WRITE at D-104, for sotolon's de-novo sugar
    # stand-in). A Strecker degradation has no business drawing sugar carbon, and the only row that
    # did was the row that was not a Strecker degradation.
    *_MAILLARD_PRECURSORS.values(),
}
# Per-product carbon fraction, keyed by pool (for the closed-form carbon accounting).
_MAILLARD_C = {
    "methional": _METHIONAL_C,
    "phenylacetaldehyde": _PHENYLACET_C,
    "2_methylbutanal": _2MB_C,
    "3_methylbutanal": _3MB_C,
    "2_methylpropanal": _2MP_C,
}


@pytest.fixture
def maillard_store():
    # Wine params + the thermal.yaml constants (k_maillard_strecker, E_a_maillard_strecker, the six
    # w_maillard_* weights) and — for the shared K_amino_acids gate — wine_generic.yaml, the
    # shared_files the D-87 compile seam wires. (The plain ``store`` fixture omits thermal.yaml.)
    return load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "thermal.yaml"
    )


@pytest.fixture
def maillard_params(maillard_store):
    return maillard_store.resolve()


def _maillard_wine(
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    aa: float = 0.3,
    s: float = 80.0,
    t: float = 298.15,
    **kw,
) -> FloatArray:
    """A finished, SEALED (o2 = 0 — the whole point) SWEET wine at the start of aging: residual
    sugar ``s`` (the dicarbonyl driver) + dosed amino acids at MUST-SPECTRUM composition (D-100 —
    the state where every per-precursor gate provably equals the pre-split lumped gate), warm."""
    y = _aged_wine(schema, ester=0.0, t=t)
    seed_amino_acids(y, schema, params, aa)
    y[schema.slice("S")] = s
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def _maillard_closed_form(
    schema: StateSchema, params: dict[str, float], y: FloatArray, t: float
) -> dict[str, float]:
    """The Process's own algebra, recomputed independently for the closed-form assertions."""
    # At must-spectrum composition every per-precursor gate collapses to this one lumped value
    # (the D-100 reduction property), so the pre-split algebra is EXACT here.
    aa = sum(float(y[schema.slice(spec.pool)][0]) for spec in AMINO_ACID_SPECS)
    s_total = float(y[schema.slice("S")].sum())
    gate = aa / (params["K_amino_acids"] + aa)
    f_t = arrhenius_factor(t, params["E_a_maillard_strecker"], params["T_ref"])
    driver = params["k_maillard_strecker"] * f_t * s_total  # PRE-gate since D-104
    weights = [params[wname] for (_, _, wname, _) in _MAILLARD_PRODUCTS]
    w_sum = sum(weights)
    out: dict[str, float] = {}
    co2_mol = 0.0
    for (pool, m_i, _wn, _prec), w_i in zip(_MAILLARD_PRODUCTS, weights, strict=True):
        # Every product here is a true Strecker degradation of its own amino acid, so every one is
        # gated on its precursor (no leucine, no 3-methylbutanal) and every one decarboxylates.
        # The two exception branches this loop used to carry — `de_novo` (rate ungated, carbon off
        # sugar) and `decarboxylates=False` — were sotolon's alone, and sotolon left at D-107.
        n_i = (w_i / w_sum) * driver * gate
        out[pool] = n_i * m_i
        co2_mol += n_i
    out["CO2"] = co2_mol * M_CO2
    return out


def test_maillard_metadata():
    p = MaillardStrecker()
    assert p.name == "maillard_strecker"
    # Speculative: the aging axis is the Tier-3 frontier (form sourced, magnitudes estimated).
    assert p.tier is Tier.SPECULATIVE
    # Books the FIVE thermal Strecker aldehydes + the decarboxylation CO2, drawing carbon from each
    # product's own precursor and deaminating the nitrogen to N. Still NO o2: the thermal axis is
    # the
    # oxygen-free one, which is the whole point of the D-87 split from D-75.
    assert set(p.touches) == _MAILLARD_TOUCHES
    assert "o2" not in p.touches
    # D-107: S is a READ-ONLY driver again (it was a WRITE at D-104 for sotolon's de-novo sugar
    # share), and sotolon is not this Process's product at all any more.
    assert "S" not in p.touches
    assert "sotolon" not in p.touches
    assert "amino_acids" not in p.touches  # the lumped draw is RETIRED (D-100)
    assert set(p.reads) == {
        "k_maillard_strecker",
        "E_a_maillard_strecker",
        "w_maillard_methional",
        "w_maillard_phenylacetaldehyde",
        "w_maillard_2_methylbutanal",
        "w_maillard_3_methylbutanal",
        "w_maillard_2_methylpropanal",
        "K_amino_acids",
        "T_ref",
        *(f"must_aa_fraction_{aa}" for aa in set(_MAILLARD_PRECURSORS.values())),
    }
    # w_maillard_sotolon is RETIRED (D-107): sotolon is not one of six co-produced Strecker
    # products competing for a shared flux, so it has no relative composition weight. It has its
    # own second-order rate constant instead (k_sotolon_aldol).
    assert not any(r.startswith("w_maillard_sotolon") for r in p.reads)


def test_maillard_matches_closed_form(maillard_params):
    schema = wine_schema()
    y = _maillard_wine(
        schema, maillard_params, aa=0.3, s=80.0, t=298.15
    )  # off T_ref so the Arrhenius factor bites
    d = MaillardStrecker().derivatives(0.0, y, schema, maillard_params)
    cf = _maillard_closed_form(schema, maillard_params, y, 298.15)

    assert cf["methional"] > 0.0  # products are live (guards against a vacuous pass)
    for pool in _MAILLARD_C:
        assert schema.get(d, pool) == pytest.approx(cf[pool])
    assert schema.get(d, "CO2") == pytest.approx(cf["CO2"])
    # EACH precursor drawn sized to the carbon of the product IT made, plus the CO2 that product's
    # own carboxyl released (D-100). The CO2 attribution is invisible to conservation, so it is
    # pinned here — and charging it is exactly what makes each carbon-sized draw land on the true
    # 1 mol precursor per mol product (D-105's signature).
    per_precursor: dict[str, float] = {}
    for pool, m_i, _wname, precursor in _MAILLARD_PRODUCTS:
        carbon = cf[pool] * _MAILLARD_C[pool] + (cf[pool] / m_i) * M_CO2 * _CO2_C
        per_precursor[precursor] = per_precursor.get(precursor, 0.0) + carbon
    expected_n = 0.0
    for precursor, carbon in per_precursor.items():
        mass = carbon / carbon_mass_fraction(precursor)
        assert schema.get(d, precursor) == pytest.approx(-mass)
        expected_n += mass * nitrogen_mass_fraction(precursor)
    assert schema.get(d, "N") == pytest.approx(expected_n)
    assert schema.get(d, "amino_acids") == 0.0  # the lumped pool is untouched (D-100)
    # S IS NOT CONSUMED since D-107: it is the dicarbonyl DRIVER, read and never written. The only
    # row that drew sugar carbon was sotolon's de-novo share (D-104), and sotolon now takes its
    # carbon from the tracked alpha_ketobutyrate pool in its own Process.
    assert float(d[schema.slice("S")].sum()) == 0.0
    assert schema.get(d, "o2") == 0.0
    for var in ("X", "E", "isoamyl_acetate", "isoamyl_alcohol", "Byp", "acetaldehyde"):
        assert schema.get(d, var) == 0.0


def test_maillard_carbon_closes_per_rhs(maillard_params):
    # CARBON closes to machine precision: the carbon leaving EACH precursor, at ITS OWN carbon
    # fraction (D-100), equals the carbon entering all six products + CO2 — a pure on-ledger
    # transfer. Five distinct molecules now, each weighted as itself rather than as arginine.
    schema = wine_schema()
    d = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params), schema, maillard_params
    )
    carbon_residual = (
        sum(schema.get(d, pool) * _MAILLARD_C[pool] for pool in _MAILLARD_C)
        + schema.get(d, "CO2") * _CO2_C
        + sum(
            schema.get(d, aa) * carbon_mass_fraction(aa)
            for aa in set(_MAILLARD_PRECURSORS.values())
        )
        # The D-104 sugar leg is gone (D-107): sotolon's de-novo share left with sotolon, so this
        # Process is a pure amino-acid -> aldehyde + CO2 transfer again. `test_maillard_matches_
        # closed_form` asserts dS == 0 head-on, so dropping the term here is not hiding a leg.
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)


def test_maillard_nitrogen_closes_per_rhs(maillard_params):
    # NITROGEN closes: all the nitrogen leaving the five precursors, at THEIR OWN fractions
    # (D-100), lands in the N pool (every product is nitrogen-free — the deamination, the
    # D-45/D-75 idiom).
    schema = wine_schema()
    d = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params), schema, maillard_params
    )
    nitrogen_residual = (
        sum(
            schema.get(d, aa) * nitrogen_mass_fraction(aa)
            for aa in set(_MAILLARD_PRECURSORS.values())
        )
        + schema.get(d, "N") * 1.0
    )
    assert nitrogen_residual == pytest.approx(0.0, abs=1e-18)


def test_maillard_is_wholly_inert_without_amino_acids(maillard_params):
    """The unit-level isolability D-104 had to give up, RESTORED at D-107.

    D-87/D-100 asserted the whole Process was inert at aa=0. D-104 broke that -- correctly, given
    where sotolon then lived: sotolon is de-novo-dominated (Crepin: "19% consumed threonine ... 81%
    newly synthesized"), so gating it on the must pool made a threonine-free wine produce no sotolon
    at all, false for every aged Sauternes. ``de_novo=True`` un-gated its rate and sourced the
    shortfall off sugar, which meant this Process fired at aa=0 and prime directive #3 fell back
    entirely onto the compile seam.

    D-107 gets the property back **without** re-breaking sotolon, because the exception was never
    about this Process: sotolon is an aldol, not a Strecker degradation, and it now lives in
    :class:`SotolonAldolCondensation` with its carbon coming from a tracked keto-acid pool. What is
    left here is five true degradations OF an amino acid -- no leucine, no 3-methylbutanal -- so
    aa=0 silences every one of them exactly, as it always should have. Pinned per-product so a
    future change cannot silently un-gate one.
    """
    schema = wine_schema()
    for aa in (0.0, -1e-6):  # 0 and the solver-undershoot guard
        d = MaillardStrecker().derivatives(
            0.0, _maillard_wine(schema, maillard_params, aa=aa), schema, maillard_params
        )
        for pool, _m, _w, _prec in _MAILLARD_PRODUCTS:
            assert schema.get(d, pool) == 0.0, f"{pool} fired at aa={aa}"
        for precursor in set(_MAILLARD_PRECURSORS.values()):
            assert schema.get(d, precursor) == 0.0
        assert schema.get(d, "N") == 0.0  # nothing deaminated
        # …and no sugar drawn: the D-104 de-novo leg left with sotolon.
        assert float(d[schema.slice("S")].sum()) == 0.0
        assert np.array_equal(d, schema.zeros())  # byte-for-byte the no-op


def test_scenario_isolability_disables_the_thermal_route_when_undosed():
    # The compile seam is the OUTER isolability guarantee (prime directive #3). D-104 made it the
    # ONLY one -- the Process fired at aa=0 via sotolon's de-novo sugar draw, so the unit-level
    # belt-and-braces was gone. D-107 restored the inner one (see
    # test_maillard_is_wholly_inert_without_amino_acids), so this is belt-and-braces again rather
    # than the sole line of defence. Kept either way: the two guarantees are independent, and a
    # seam regression is not a unit-level one.
    from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

    compiled = compile_scenario(
        Scenario(
            name="undosed",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25},  # NO amino_acids_gpl
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        ),
        strict=True,
    )
    assert not compiled.process_set.is_enabled(MaillardStrecker.name)


def test_maillard_inert_without_sugar(maillard_params):
    # No residual sugar ⇒ no dicarbonyls ⇒ nothing: a dry wine (S = 0) makes ~none. A SOFT driver
    # (unlike the aa hard gate) — physically a dry wine still makes a trace, but the S=0 guard is a
    # clean no-op here (and absorbs the S<0 solver undershoot).
    schema = wine_schema()
    assert np.array_equal(
        MaillardStrecker().derivatives(
            0.0, _maillard_wine(schema, maillard_params, s=0.0), schema, maillard_params
        ),
        schema.zeros(),
    )
    assert np.array_equal(
        MaillardStrecker().derivatives(
            0.0, _maillard_wine(schema, maillard_params, s=-1e-6), schema, maillard_params
        ),
        schema.zeros(),
    )


def test_maillard_is_oxygen_independent(maillard_params):
    # THE discriminating property (the mirror of D-75): this route needs NO O₂, so a SEALED wine
    # (o2 = 0) produces EXACTLY the same as an oxygenated one. Contrast StreckerDegradation, which
    # is
    # zero at o2 = 0. Same sweet state, o2 swept 0 → high: byte-for-byte identical derivative.
    schema = wine_schema()
    sealed = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, o2=0.0), schema, maillard_params
    )
    oxygenated = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, o2=0.05), schema, maillard_params
    )
    assert np.array_equal(sealed, oxygenated)
    assert schema.get(sealed, "methional") > 0.0  # non-vacuous: it IS producing while sealed


def test_maillard_first_order_in_sugar(maillard_params):
    # The residual sugar is the dicarbonyl driver: doubling S doubles the production (first-order),
    # holding the aa gate fixed. The bounded-vs-unbounded concern is on amino_acids (the limiting
    # reagent), not sugar — sugar drives the RATE.
    schema = wine_schema()
    base = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, s=40.0), schema, maillard_params
    )
    dbl = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, s=80.0), schema, maillard_params
    )
    assert schema.get(dbl, "methional") == pytest.approx(2.0 * schema.get(base, "methional"))


def test_maillard_availability_gate_saturates(maillard_params):
    # The aa availability gate aa/(K+aa) saturates: production per unit sugar rises then plateaus as
    # amino_acids climbs (the same smooth-Monod shape D-75 uses). Below/at/well-above K.
    schema = wine_schema()
    k = maillard_params["K_amino_acids"]
    # Equal-RATIO aa steps (k → 10k → 100k), so a saturating gate shows strictly diminishing gains.
    low = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, aa=k), schema, maillard_params
    )
    mid = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, aa=10.0 * k), schema, maillard_params
    )
    high = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, aa=100.0 * k), schema, maillard_params
    )
    # D-104 had to read this on a specifically-chosen GATED product, because sotolon was `de_novo`
    # and therefore deliberately flat in aa. Since D-107 every product here is gated, so the choice
    # no longer carries that caveat -- any of the five would do.
    lo, md, hi = (schema.get(x, "3_methylbutanal") for x in (low, mid, high))
    assert lo < md < hi  # monotone increasing in aa
    assert (md - lo) > (hi - md)  # but saturating (diminishing returns) — the gate flattens


def test_maillard_split_normalizes_and_every_product_charges_one_co2(maillard_params):
    # The five composition weights NORMALIZE to fractions summing to 1 (the split-hygiene the
    # advisor flagged), and EVERY product charges exactly 1 CO2 -- so the CO2 mole rate equals the
    # TOTAL product mole rate. Through D-106 this test asserted the opposite: that sotolon (a
    # furanone, NOT a decarboxylation product) was EXCLUDED from the CO2 sum. D-107 moved sotolon
    # out of this Process, so the exception it guarded no longer exists and the assertion inverts
    # from a strict subset to an equality. It stays load-bearing for the reason it always was:
    # conservation closes for ANY CO2 attribution, so only an explicit check pins the keying -- and
    # charging that CO2 is exactly what makes each carbon-sized draw land on the true 1:1 (D-105).
    schema = wine_schema()
    y = _maillard_wine(schema, maillard_params)
    d = MaillardStrecker().derivatives(0.0, y, schema, maillard_params)
    # mole rate of each product = mass rate / molar mass
    masses = {
        "methional": M_METHIONAL,
        "phenylacetaldehyde": M_PHENYLACETALDEHYDE,
        "2_methylbutanal": M_2_METHYLBUTANAL,
        "3_methylbutanal": M_3_METHYLBUTANAL,
        "2_methylpropanal": M_2_METHYLPROPANAL,
    }
    n = {pool: schema.get(d, pool) / masses[pool] for pool in masses}
    n_total = sum(n.values())
    assert n_total > 0.0  # non-vacuous
    # The normalized split sums to 1 (fractions = n_i / n_total).
    assert sum(n[pool] / n_total for pool in n) == pytest.approx(1.0)
    # CO2 == EVERY product's mole sum: all five decarboxylate (D-107), so D-87's strict subset is
    # now an equality.
    assert schema.get(d, "CO2") / M_CO2 == pytest.approx(n_total)
    # And sotolon is not this Process's business at all any more.
    assert schema.get(d, "sotolon") == 0.0


def test_maillard_rises_with_temperature(maillard_params):
    # Warmer ages faster — the sourced direction (thermal Strecker is strongly temperature-driven).
    schema = wine_schema()
    cold = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, t=283.15), schema, maillard_params
    )
    warm = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, maillard_params, t=303.15), schema, maillard_params
    )
    assert schema.get(warm, "methional") > schema.get(cold, "methional") > 0.0


def test_maillard_is_more_thermally_sensitive_than_oxidative(maillard_params, params):
    # The sourced ORDERING (D-87): Maillard/caramelization out-accelerates the oxidative aging
    # reactions with temperature, so E_a_maillard_strecker > the oxidative E_a's. The Q10 (ratio of
    # a 10 K temperature step) is strictly LARGER for the thermal route than for the oxidative one.
    def q10(e_a: float) -> float:
        return arrhenius_factor(303.15, e_a, maillard_params["T_ref"]) / arrhenius_factor(
            293.15, e_a, maillard_params["T_ref"]
        )

    assert q10(maillard_params["E_a_maillard_strecker"]) > q10(params["E_a_strecker"])


def test_maillard_is_wine_only_noop_on_beer(maillard_params):
    # amino_acids is wine-only (D-32), so — like StreckerDegradation — MaillardStrecker is a hard
    # no-op on beer (the "amino_acids"/"sotolon" not in schema guard), even with residual wort
    # sugar.
    beer = beer_schema()
    yb = beer.zeros()
    yb[beer.slice("S")] = 50.0
    yb[beer.slice("T")] = 298.15
    assert np.array_equal(
        MaillardStrecker().derivatives(0.0, yb, beer, maillard_params), beer.zeros()
    )


def test_maillard_sealed_sweet_accumulates_where_oxidative_gives_zero(
    maillard_store, maillard_params
):
    # THE discriminating end-to-end (ProcessSet integration): a SEALED (o2 = 0) SWEET (residual S)
    # amino-acid-dosed wine accumulates the thermal aldehydes over a warm aging year — WHERE the
    # O₂-only StreckerDegradation, on the identical sealed state, produces exactly ZERO (no O₂). The
    # single contrast that proves the route does something the oxidative sink cannot. Carbon AND
    # nitrogen close to machine precision through the integrated segment.
    schema = wine_schema()
    y0 = _maillard_wine(
        schema, maillard_params, aa=0.4, s=80.0, t=301.15
    )  # 28 °C, sealed (o2 defaults to 0)
    # The oxidative route on the SAME sealed state: identically zero (its o2 <= 0 guard short-
    # circuits before it even reads a param, so maillard_params suffices — nothing is computed).
    d_oxid = StreckerDegradation().derivatives(0.0, y0, schema, maillard_params)
    assert np.array_equal(d_oxid, schema.zeros())

    ps = ProcessSet(schema, [MaillardStrecker()], strict=True)
    traj = simulate(ps, params=maillard_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    for pool in ("methional", "phenylacetaldehyde", "2_methylpropanal", "3_methylbutanal"):
        series = traj.series(pool)
        assert series[-1] > series[0] == 0.0  # produced-only, monotone accumulation from 0
    assert_nonnegative(traj, ("amino_acids", "methional", "N"), atol=1e-9)
    f_c = maillard_store.value("biomass_C_fraction")
    f_n = maillard_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_maillard_tier_floored_at_speculative(maillard_store):
    # Speculative in FORM (Tier-3 frontier): every pool it writes is speculative. Non-vacuous across
    # all nine touched pools — including the structural N-write (deamination, the D-45/D-75 note).
    schema = wine_schema()
    ps = ProcessSet(schema, [MaillardStrecker()])
    for pool in _MAILLARD_TOUCHES:
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, maillard_store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# SotolonAldolCondensation (decision D-107) — the keto-acid NODE's consumer: the purely chemical
# aldol of alpha-ketobutyrate + acetaldehyde. NOT a Strecker route: no sugar driver, no CO2, no N.
# =====================================================================================


def _aldol_wine(
    schema: StateSchema,
    *,
    keto: float = 2.0e-3,  # ~2 mg/L — the D-107 excreted residual
    acetaldehyde: float = 33.0e-3,  # ~33 mg/L — what a stuck/sweet ferment leaves
    s: float = 80.0,
    t: float = 298.15,
) -> FloatArray:
    y = schema.zeros()
    y[schema.slice("T")] = t
    y[schema.slice("S")] = s
    y[schema.slice("alpha_ketobutyrate")] = keto
    y[schema.slice("acetaldehyde")] = acetaldehyde
    return y


def test_sotolon_aldol_metadata():
    p = SotolonAldolCondensation()
    assert p.name == "sotolon_aldol_condensation"
    assert p.tier is Tier.SPECULATIVE
    # Its two real substrates and its one product. The three ABSENCES are the content of D-107:
    # no S (an aldol needs substrates, not a dicarbonyl driver), no CO2 (nothing decarboxylates —
    # there is no carboxyl to lose), no N (both substrates are nitrogen-free).
    assert set(p.touches) == {"sotolon", "alpha_ketobutyrate", "acetaldehyde"}
    assert "S" not in p.touches
    assert "CO2" not in p.touches
    assert "N" not in p.touches
    assert "threonine" not in p.touches  # threonine is sotolon's GRANDparent, not its precursor
    assert set(p.reads) == {"k_sotolon_aldol", "E_a_maillard_strecker", "T_ref"}


def test_the_sotolon_aldol_draws_one_mole_of_each_substrate_when_driven(maillard_params):
    """The D-105 signature, satisfied by construction rather than by sizing (decision D-107).

    Sotolon's entry on ``_KNOWN_NON_STOICHIOMETRIC`` said it OVER-drew threonine 1.5×: all six of
    its carbons were taken from a C4 amino acid, because its two acetaldehyde-derived carbons were
    lumped into that draw (the D-87 scope note). That ratio was never fixable by re-sizing — it was
    the wrong question. Sotolon is an aldol of two molecules, so it consumes **1 mol of each**, and
    the carbon closes because ``4 + 2 == 6`` on the atom counts, not because a mass was chosen to
    make it close.

    That distinction is the point: a route with no carbon-sized draw **has no D-105 blind spot**.
    This test reads the stoichiometry off ``dy/dt``, so it is the code layer, not a declaration.
    """
    schema = wine_schema()
    d = SotolonAldolCondensation().derivatives(0.0, _aldol_wine(schema), schema, maillard_params)
    n_sot = schema.get(d, "sotolon") / M_SOTOLON
    n_keto = -schema.get(d, "alpha_ketobutyrate") / M_ALPHA_KETOBUTYRATE
    n_acet = -schema.get(d, "acetaldehyde") / M_ACETALDEHYDE
    assert n_sot > 0.0  # non-vacuous
    assert n_keto / n_sot == pytest.approx(1.0, abs=1e-12)
    assert n_acet / n_sot == pytest.approx(1.0, abs=1e-12)


def test_sotolon_aldol_carbon_closes_on_atom_counts(maillard_params):
    # C6 out == C4 + C2 in. This is the strongest closure in the tree: it holds because the
    # chemistry balances, not because the draw was sized to make it hold.
    schema = wine_schema()
    d = SotolonAldolCondensation().derivatives(0.0, _aldol_wine(schema), schema, maillard_params)
    residual = (
        schema.get(d, "sotolon") * carbon_mass_fraction("sotolon")
        + schema.get(d, "alpha_ketobutyrate") * carbon_mass_fraction("alpha_ketobutyrate")
        + schema.get(d, "acetaldehyde") * carbon_mass_fraction("acetaldehyde")
    )
    assert residual == pytest.approx(0.0, abs=1e-18)
    # And the atom counts really are what balances it (the identity this route rests on).
    assert CARBON_ATOMS["sotolon"] == (
        CARBON_ATOMS["alpha_ketobutyrate"] + CARBON_ATOMS["acetaldehyde"]
    )


def test_sotolon_aldol_is_bimolecular_in_both_substrates(maillard_params):
    """Mass-action in BOTH substrates — the rate law Pham et al. 1995 measured (D-107).

    Doubling either substrate doubles the rate; doubling both quadruples it. This is what makes the
    route sugar-independent (an aldol needs its two substrates, not a dicarbonyl), and it is the
    MUTATION-SENSITIVE property: a rate that was first-order in one substrate, or that kept the
    D-87 sugar driver, fails here.
    """
    schema = wine_schema()
    p = SotolonAldolCondensation()

    def rate(**kw: float) -> float:
        d = p.derivatives(0.0, _aldol_wine(schema, **kw), schema, maillard_params)
        return float(schema.get(d, "sotolon"))

    base = rate()
    assert base > 0.0
    assert rate(keto=4.0e-3) == pytest.approx(2.0 * base)
    assert rate(acetaldehyde=66.0e-3) == pytest.approx(2.0 * base)
    assert rate(keto=4.0e-3, acetaldehyde=66.0e-3) == pytest.approx(4.0 * base)


def test_sotolon_aldol_is_sugar_independent(maillard_params):
    """The empirical correction D-107 makes, stated as an assertion (Pons et al. 2010).

    Inside MaillardStrecker sotolon's rate was ``k · f(T) · S`` — pseudo-first-order in residual
    sugar, which made it a SWEETNESS marker. Pons *et al.* 2010 identified this same aldol as the
    sotolon pathway in prematurely aged DRY white wines, where it is the premature-oxidation
    marker. Sugar is absent from the rate law because it is absent from the reaction, so a bone-dry
    state and a botrytis-sweet one give byte-for-byte the same derivative.
    """
    schema = wine_schema()
    p = SotolonAldolCondensation()
    dry = p.derivatives(0.0, _aldol_wine(schema, s=0.0), schema, maillard_params)
    sweet = p.derivatives(0.0, _aldol_wine(schema, s=150.0), schema, maillard_params)
    assert np.array_equal(dry, sweet)
    assert schema.get(dry, "sotolon") > 0.0  # non-vacuous: it fires with NO sugar at all


def test_sotolon_aldol_is_exactly_isolable_without_the_keto_acid_pool(maillard_params):
    """Isolability is EXACT and free, which is WHY the rate is mass-action (decision D-107).

    A ProcessSet built without ``_KETO_ACID_PROCESSES`` leaves ``alpha_ketobutyrate`` at 0, and a
    rate that is the product of its substrates is then **exactly** 0 — no clamp, no epsilon, no
    availability constant. The alternative (keep the sugar driver, gate the draws on the pools)
    would have needed a FABRICATED half-saturation per substrate: the faithful rate law is the one
    with FEWER invented numbers. Both substrates are checked, plus the solver-undershoot guard.
    """
    schema = wine_schema()
    p = SotolonAldolCondensation()
    for kw in ({"keto": 0.0}, {"acetaldehyde": 0.0}, {"keto": -1e-9}, {"acetaldehyde": -1e-9}):
        d = p.derivatives(0.0, _aldol_wine(schema, **kw), schema, maillard_params)
        assert np.array_equal(d, schema.zeros()), kw


def test_sotolon_aldol_rises_with_temperature(maillard_params):
    # Pham et al. 1995: "the formation of sotolon increases by increasing temperature". Only the
    # DIRECTION is sourced — the magnitude rides E_a_maillard_strecker as a labelled carry-over
    # (sotolon already rode that constant from D-87 to D-106), because inventing a sotolon-specific
    # activation energy to look precise is the E_a D-101 fabricated and D-102 had to retract.
    schema = wine_schema()
    p = SotolonAldolCondensation()
    cold = p.derivatives(0.0, _aldol_wine(schema, t=283.15), schema, maillard_params)
    warm = p.derivatives(0.0, _aldol_wine(schema, t=303.15), schema, maillard_params)
    assert schema.get(warm, "sotolon") > schema.get(cold, "sotolon") > 0.0


def test_sotolon_aldol_is_wine_only_noop_on_beer(maillard_params):
    # sotolon/alpha_ketobutyrate are wine-only slots, so this is a hard no-op on beer even with
    # residual wort sugar and acetaldehyde present.
    beer = beer_schema()
    yb = beer.zeros()
    yb[beer.slice("S")] = 50.0
    yb[beer.slice("T")] = 298.15
    yb[beer.slice("acetaldehyde")] = 30.0e-3
    assert np.array_equal(
        SotolonAldolCondensation().derivatives(0.0, yb, beer, maillard_params), beer.zeros()
    )


def test_sotolon_aldol_tier_floored_at_speculative(maillard_store):
    schema = wine_schema()
    ps = ProcessSet(schema, [SotolonAldolCondensation()])
    for pool in ("sotolon", "alpha_ketobutyrate", "acetaldehyde"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, maillard_store.tier_map()) is Tier.SPECULATIVE


# -------------------------------------------------------------------------------------
# D-108: the aldol reads FREE acetaldehyde, not total. SO₂-bound acetaldehyde is the
# bisulfite adduct — its carbonyl is BLOCKED, and an aldol condensation IS a nucleophilic
# attack on that carbonyl, so the adduct cannot condense. Identical to the argument
# AcetaldehydeBridging (D-80) and the tannin polymerization already make for the ethylidene
# bridge, and AcetaldehydeReduction (D-47) for ADH. D-107 read the TOTAL slot while its
# docstring claimed it "reads the pool the binding depletes" — but `free = total − bound` is
# a derived overlay and the binding depletes nothing, so SO₂ came out RAISING sotolon.
# NOTE these are the FIRST tests of sotolon under SO₂ at all: the pre-D-108 code moved
# sulfited sotolon by 85× and all 1152 tests stayed green (the D-105 tripwire lesson).
# -------------------------------------------------------------------------------------


@pytest.fixture
def aldol_so2_store():
    # thermal.yaml (k_sotolon_aldol / E_a) + acidbase.yaml (pKas + the D-51 binding constants
    # free_acetaldehyde reads) + wine_generic.yaml. The plain maillard_store omits acidbase.
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "thermal.yaml",
        default_data_dir() / "acidbase.yaml",
    )


@pytest.fixture
def aldol_so2_params(aldol_so2_store):
    return aldol_so2_store.resolve()


def _aldol_sulfitable(
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    so2_mgl: float,
    acetaldehyde_mgl: float = 50.0,
) -> FloatArray:
    """An aging wine at pH ~3.4 with α-ketobutyrate, acetaldehyde and (maybe) SO₂ dosed."""
    tartaric, malic = 6.0, 3.0
    totals = {"tartaric": tartaric / M_TARTARIC, "malic": malic / M_MALIC, "lactic": 0.0}
    cation = solve_cation_charge(totals, 0.0, build_pka_map(params), 3.4)
    y = schema.zeros()
    y[schema.slice("T")] = 298.15
    y[schema.slice("tartaric")] = tartaric
    y[schema.slice("malic")] = malic
    y[schema.slice("cation_charge")] = cation
    y[schema.slice("alpha_ketobutyrate")] = 2.0e-3
    y[schema.slice("acetaldehyde")] = mgl_to_gpl(acetaldehyde_mgl)
    if so2_mgl > 0.0:
        y[schema.slice("so2_total")] = mgl_to_gpl(so2_mgl)
    return y


def _sotolon_rate(schema: StateSchema, params: Mapping[str, float], y: FloatArray) -> float:
    d = SotolonAldolCondensation().derivatives(0.0, y, schema, params)
    return float(d[schema.slice("sotolon")][0])


def test_unsulfited_aldol_is_byte_for_byte_the_total_acetaldehyde_form(aldol_so2_params):
    # The `so2_total > 0` guard is EXACT: with no SO₂ the rate reads the TOTAL pool (free ==
    # total), no per-RHS pH brentq is paid, and the Process is byte-for-byte the D-107 code.
    # This is why every output D-107 measured — all of them unsulfited — is unmoved by D-108.
    schema = wine_schema()
    y = _aldol_sulfitable(schema, aldol_so2_params, so2_mgl=0.0)
    ph = ph_of_state(y, schema, aldol_so2_params)
    total = float(y[schema.slice("acetaldehyde")][0])
    assert free_acetaldehyde(y, schema, aldol_so2_params, ph) == pytest.approx(total)
    # Written from the CHEMISTRY (Pham's bimolecular aldol), not re-derived from the code —
    # D-107 lesson (iv): the same arithmetic twice cannot find an error in itself.
    f_t = arrhenius_factor(
        298.15, aldol_so2_params["E_a_maillard_strecker"], aldol_so2_params["T_ref"]
    )
    n_expected = (
        aldol_so2_params["k_sotolon_aldol"]
        * f_t
        * (2.0e-3 / M_ALPHA_KETOBUTYRATE)
        * (total / M_ACETALDEHYDE)
    )
    assert _sotolon_rate(schema, aldol_so2_params, y) == pytest.approx(n_expected * M_SOTOLON)


def test_so2_throttles_the_aldol_to_the_free_share(aldol_so2_params):
    # THE D-108 FIX. Bound acetaldehyde has no free carbonyl ⇒ it cannot condense ⇒ dosing SO₂
    # must LOWER the sotolon rate, in proportion to free/total. Before D-108 this ran backwards:
    # SO₂ strands acetaldehyde (D-47 protects it from ADH), the rate read the swollen TOTAL, and
    # a sulfited dry wine made MORE sotolon — against Pons, for whom low free SO₂ is the prémox
    # RISK factor.
    schema = wine_schema()
    y_clean = _aldol_sulfitable(schema, aldol_so2_params, so2_mgl=0.0)
    y_comparable = _aldol_sulfitable(schema, aldol_so2_params, so2_mgl=60.0)
    y_excess = _aldol_sulfitable(schema, aldol_so2_params, so2_mgl=400.0)
    r_clean = _sotolon_rate(schema, aldol_so2_params, y_clean)
    r_comparable = _sotolon_rate(schema, aldol_so2_params, y_comparable)
    r_excess = _sotolon_rate(schema, aldol_so2_params, y_excess)
    assert r_clean > 0.0
    assert r_comparable < r_clean  # the DIRECTION is the whole point: SO₂ protects
    assert r_excess < r_comparable  # more SO₂ ⇒ less sotolon, monotone
    assert r_excess < 0.02 * r_clean  # SO₂ ≫ acetaldehyde ⇒ ~no free carbonyl left to condense
    # Pin the rate to the free-acetaldehyde readout exactly (the D-47 assertion shape).
    ph = ph_of_state(y_comparable, schema, aldol_so2_params)
    free = free_acetaldehyde(y_comparable, schema, aldol_so2_params, ph)
    f_t = arrhenius_factor(
        298.15, aldol_so2_params["E_a_maillard_strecker"], aldol_so2_params["T_ref"]
    )
    n_expected = (
        aldol_so2_params["k_sotolon_aldol"]
        * f_t
        * (2.0e-3 / M_ALPHA_KETOBUTYRATE)
        * (free / M_ACETALDEHYDE)
    )
    assert r_comparable == pytest.approx(n_expected * M_SOTOLON)


def test_the_aldol_draw_debits_TOTAL_acetaldehyde_even_though_the_rate_reads_free(aldol_so2_params):
    """The D-47 idiom, and the half that keeps carbon closing (decision D-108).

    Only the RATE reads the free share. The DRAW still debits the ``acetaldehyde`` slot — which
    holds the TOTAL — because consuming free acetaldehyde removes it from the total and the
    binding equilibrium then re-splits what is left. Booking the draw against a derived "free
    pool" instead would debit a quantity that is not a state slot, and carbon would stop closing
    on the atom counts. So the 1:1:1 mole relation must survive SO₂ untouched.
    """
    schema = wine_schema()
    y = _aldol_sulfitable(schema, aldol_so2_params, so2_mgl=60.0)
    d = SotolonAldolCondensation().derivatives(0.0, y, schema, aldol_so2_params)
    n_sot = float(d[schema.slice("sotolon")][0]) / M_SOTOLON
    assert n_sot > 0.0
    # 1 mol of each substrate per mol of sotolon — unchanged under SO₂ (4 + 2 == 6 still closes).
    assert float(d[schema.slice("alpha_ketobutyrate")][0]) == pytest.approx(
        -n_sot * M_ALPHA_KETOBUTYRATE
    )
    assert float(d[schema.slice("acetaldehyde")][0]) == pytest.approx(-n_sot * M_ACETALDEHYDE)


def test_the_aldol_asymptotes_toward_zero_under_excess_so2_but_never_reaches_it(aldol_so2_params):
    """The end of the ladder is an ASYMPTOTE, not a hard zero — and this test first claimed a zero.

    The binding is an equilibrium (``free = total − bound`` off the D-51 split), so ``bound`` stays
    below ``total`` and a vanishing free share survives however much SO₂ is dosed: at **5000 mg/L
    SO₂ against
    1 mg/L acetaldehyde**, ``free_acetaldehyde`` still returns ~2e-8 g/L. This test originally
    asserted ``free <= 0`` and an exactly-zero derivative — a hard zero the chemistry cannot make.
    So the ``acetaldehyde <= 0.0`` early-return in the Process is a **defensive mirror** of the D-47
    idiom rather than a reachable branch, and is documented as such rather than claimed as coverage.

    What IS true, and is what a sulfited wine actually does, is asserted instead: the rate is driven
    arbitrarily close to zero, and stays strictly positive.
    """
    schema = wine_schema()
    y = _aldol_sulfitable(schema, aldol_so2_params, so2_mgl=5000.0, acetaldehyde_mgl=1.0)
    ph = ph_of_state(y, schema, aldol_so2_params)
    free = free_acetaldehyde(y, schema, aldol_so2_params, ph)
    assert 0.0 < free < 1.0e-6  # vanishing, but an equilibrium never gives it up entirely
    y_clean = _aldol_sulfitable(schema, aldol_so2_params, so2_mgl=0.0, acetaldehyde_mgl=1.0)
    r_excess = _sotolon_rate(schema, aldol_so2_params, y)
    assert 0.0 < r_excess < 1.0e-4 * _sotolon_rate(schema, aldol_so2_params, y_clean)


# =====================================================================================
# Caramelization (decision D-88; MEDIUM-AGNOSTIC D-90) — the NON-oxidative THERMAL browning: the
# O₂-INDEPENDENT thermal mirror of PhenolicBrowning (D-74). Residual SUGAR browns to melanoidin by
# HEAT (no O₂), raising the SAME A420 index D-74 accumulates — so a sealed sweet wine *or*
# high-residual beer still darkens. The FIRST aging Process to consume core S: the sugar carbon
# lands in the on-ledger melanoidin carbon-park (the debris/glucan precedent), so carbon closes
# exactly (release at each sugar's own fraction, redeposit at melanoidin's — the D-90 vectorized
# draw apportions across beer's 3-slot S). SUGAR-ONLY (nitrogen-free — caramelization, not Maillard;
# N-incorporating MaillardBrowning D-89 stays wine-only, beer's amino_acids untracked, D-32). These
# tests pin the closed form, carbon closure per-RHS, the sugar SOFT gate (inert at S ≈ 0 /
# undershoot), the O₂-independence (no o2 term at all), the first-order-in-sugar linearity, the
# monotone A420 rise, the warmer-faster ordering, the medium-agnostic beer + per-component carbon
# closure (D-90), the integrated sweet-wine + residual-beer browning + closure, and the speculative
# tier floor.

_CARAMEL_TOUCHES = {"S", "melanoidin", "A420"}


@pytest.fixture
def caramel_store():
    # Wine params + the thermal.yaml constants (k_caramelization, E_a_caramelization,
    # y_a420_per_melanoidin, T_ref) — the shared_files the D-88 compile seam wires.
    return load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "thermal.yaml"
    )


@pytest.fixture
def caramel_params(caramel_store):
    return caramel_store.resolve()


def _caramel_wine(schema: StateSchema, *, s: float = 100.0, t: float = 298.15, **kw) -> FloatArray:
    """A finished SWEET wine at the start of aging: residual sugar ``s``, warm, sealed (o2 = 0)."""
    y = _aged_wine(schema, ester=0.0, t=t)
    y[schema.slice("S")] = s
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def _caramel_closed_form(
    schema: StateSchema, params: dict[str, float], y: FloatArray, t: float
) -> dict[str, float]:
    """The Process's own algebra, recomputed independently."""
    s_total = float(y[schema.slice("S")].sum())
    f_t = arrhenius_factor(t, params["E_a_caramelization"], params["T_ref"])
    r = params["k_caramelization"] * f_t * s_total
    mel_rate = r * _GLUCOSE_C / _MELANOIDIN_C
    return {
        "S": -r,
        "melanoidin": mel_rate,
        "A420": params["y_a420_per_melanoidin"] * mel_rate,
    }


def test_caramelization_metadata():
    p = Caramelization()
    assert p.name == "caramelization"
    assert p.tier is Tier.SPECULATIVE
    # Consumes core S into the on-ledger melanoidin carbon-park + raises the shared A420. Touches
    # those three and nothing else — NO o2 (the whole point), no amino acids (sugar-only).
    assert set(p.touches) == _CARAMEL_TOUCHES
    assert "o2" not in p.touches
    assert set(p.reads) == {
        "k_caramelization",
        "E_a_caramelization",
        "y_a420_per_melanoidin",
        "T_ref",
    }


def test_caramelization_matches_closed_form(caramel_params):
    schema = wine_schema()
    y = _caramel_wine(schema, s=100.0, t=298.15)  # off T_ref so the Arrhenius factor bites
    d = Caramelization().derivatives(0.0, y, schema, caramel_params)
    cf = _caramel_closed_form(schema, caramel_params, y, 298.15)
    assert cf["melanoidin"] > 0.0  # live (guards against a vacuous pass)
    assert schema.get(d, "S") == pytest.approx(cf["S"])
    assert schema.get(d, "melanoidin") == pytest.approx(cf["melanoidin"])
    assert schema.get(d, "A420") == pytest.approx(cf["A420"])
    # Touches nothing else — no o2, no aroma pools, no amino acids, no E.
    for var in ("o2", "amino_acids", "E", "isoamyl_acetate", "acetaldehyde", "sotolon", "N"):
        assert schema.get(d, var) == 0.0


def test_caramelization_carbon_closes_per_rhs(caramel_params):
    # CARBON closes to machine precision: the sugar carbon leaving S equals the carbon entering the
    # melanoidin carbon-park (release at the sugar fraction, redeposit at melanoidin's). A420 is an
    # optical index (off every ledger), so it carries none.
    schema = wine_schema()
    d = Caramelization().derivatives(0.0, _caramel_wine(schema), schema, caramel_params)
    carbon_residual = schema.get(d, "S") * _GLUCOSE_C + schema.get(d, "melanoidin") * _MELANOIDIN_C
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)


def test_caramelization_inert_without_sugar(caramel_params):
    # No residual sugar ⇒ no browning: a dry wine (S = 0) is byte-for-byte inert (the SOFT sugar
    # gate — every standard dry aging run ferments to S ≈ 0 before begin_aging, so caramelization
    # leaves it unchanged). Also absorbs the S < 0 solver undershoot.
    schema = wine_schema()
    ps = ProcessSet(schema, [Caramelization()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _caramel_wine(schema, s=0.0), caramel_params), schema.zeros()
    )
    assert np.array_equal(
        Caramelization().derivatives(0.0, _caramel_wine(schema, s=-1e-6), schema, caramel_params),
        schema.zeros(),
    )


def test_caramelization_is_oxygen_independent(caramel_params):
    # NO o2 term at all: a sealed wine (o2 = 0) browns exactly as an oxygenated one — the whole
    # point
    # (the O₂-independent thermal mirror of the O₂-driven PhenolicBrowning). o2 is untouched.
    schema = wine_schema()
    sealed = Caramelization().derivatives(
        0.0, _caramel_wine(schema, o2=0.0), schema, caramel_params
    )
    oxygenated = Caramelization().derivatives(
        0.0, _caramel_wine(schema, o2=0.05), schema, caramel_params
    )
    assert np.array_equal(sealed, oxygenated)
    assert schema.get(sealed, "o2") == 0.0  # never draws o2
    assert schema.get(sealed, "melanoidin") > 0.0  # non-vacuous: browning while sealed


def test_caramelization_first_order_in_sugar(caramel_params):
    # First-order in residual sugar: doubling S doubles the browning rate (and the A420 rise).
    schema = wine_schema()
    base = Caramelization().derivatives(0.0, _caramel_wine(schema, s=50.0), schema, caramel_params)
    dbl = Caramelization().derivatives(0.0, _caramel_wine(schema, s=100.0), schema, caramel_params)
    assert schema.get(dbl, "melanoidin") == pytest.approx(2.0 * schema.get(base, "melanoidin"))
    assert schema.get(dbl, "A420") == pytest.approx(2.0 * schema.get(base, "A420"))


def test_caramelization_a420_rises_monotone(caramel_params):
    # A420 accumulates (produced-only, the D-74 optical-index idiom): d(A420)/dt ≥ 0 whenever sugar
    # is present, so the browning index is monotone (never reversed). Melanoidin likewise.
    schema = wine_schema()
    d = Caramelization().derivatives(0.0, _caramel_wine(schema), schema, caramel_params)
    assert schema.get(d, "A420") > 0.0
    assert schema.get(d, "melanoidin") > 0.0
    assert schema.get(d, "S") < 0.0  # sugar is consumed


def test_caramelization_rises_with_temperature(caramel_params):
    # Warmer browns faster — the sourced direction (caramelization is strongly thermal, why Madeira
    # estufagem browns so much faster than cellar aging).
    schema = wine_schema()
    cold = Caramelization().derivatives(
        0.0, _caramel_wine(schema, t=283.15), schema, caramel_params
    )
    warm = Caramelization().derivatives(
        0.0, _caramel_wine(schema, t=303.15), schema, caramel_params
    )
    assert schema.get(warm, "melanoidin") > schema.get(cold, "melanoidin") > 0.0


def test_caramelization_runs_on_beer_and_closes_carbon_per_component(caramel_params):
    # MEDIUM-AGNOSTIC (D-90 supersedes D-88's wine-only v1): beer's residual dextrins caramelize.
    # The vectorized draw apportions the sugar debit across the 3-slot S vector and releases each
    # component's carbon at its OWN fraction (glucose/maltose/maltotriose differ), so total_carbon
    # closes per-RHS on beer. A warm beer with ALL THREE sugars present — so every per-component
    # carbon fraction (glucose slot 0 included) is load-bearing in the closure — browns, and the
    # per-component carbon balance is exact.
    beer = beer_schema()
    yb = beer.zeros()
    s_vec = [10.0, 40.0, 20.0]  # glucose + maltose + maltotriose all residual
    yb[beer.slice("S")] = s_vec
    yb[beer.slice("T")] = 303.15
    d = Caramelization().derivatives(0.0, yb, beer, caramel_params)
    # It RUNS (not the old no-op): melanoidin + A420 climb, sugar is consumed.
    mel_rate = float(d[beer.slice("melanoidin")][0])
    assert mel_rate > 0.0
    assert float(d[beer.slice("A420")][0]) > 0.0
    dS = d[beer.slice("S")]
    # Every slot is debited (per-component apportionment, not a broadcast onto one slot).
    assert dS[0] < 0.0 and dS[1] < 0.0 and dS[2] < 0.0
    # The draw is apportioned by share: maltose (40) loses twice maltotriose (20), 4× glucose (10).
    assert float(dS[1]) == pytest.approx(2.0 * float(dS[2]))
    assert float(dS[1]) == pytest.approx(4.0 * float(dS[0]))
    # CARBON closes per-RHS: carbon leaving each S slot (at its own fraction) == carbon into
    # melanoidin.
    carbon_residual = (
        float(dS[0]) * _GLUCOSE_C
        + float(dS[1]) * _MALTOSE_C
        + float(dS[2]) * _MALTOTRIOSE_C
        + mel_rate * _MELANOIDIN_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)


def test_caramelization_browns_a_residual_beer_and_closes_carbon(caramel_store):
    # Integrated beer counterpart of the sweet-wine browning test (D-90): a warm HIGH-RESIDUAL beer
    # (under-attenuated big stout — maltose + maltotriose left) browns over an aging year through
    # STRICT ProcessSet, with total_carbon closing to machine precision across the multi-slot S
    # vector (the per-component fractions redeposit correctly). The load-bearing beer-closure test.
    beer = load_parameters(
        default_data_dir() / "beer_generic.yaml", default_data_dir() / "thermal.yaml"
    )
    params = beer.resolve()
    schema = beer_schema()
    y0 = schema.zeros()
    y0[schema.slice("S")] = [0.0, 45.0, 25.0]  # residual maltose + maltotriose
    y0[schema.slice("E")] = 40.0
    y0[schema.slice("T")] = 303.15  # 30 °C warm store
    ps = ProcessSet(schema, [Caramelization()], strict=True)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    assert float(traj.series("melanoidin")[-1]) > 0.0  # browns from 0
    assert float(traj.series("A420")[-1]) > 0.0  # A420 climbs from 0
    assert float(traj.series("S")[:, -1].sum()) < 70.0  # residual sugar declined
    assert_nonnegative(traj, ("melanoidin", "A420"), atol=1e-9)
    f_c = beer.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_caramelization_browns_a_sweet_wine_and_closes_carbon(caramel_store, caramel_params):
    # Integrated: a warm SWEET wine (residual sugar) browns over an aging year — melanoidin + A420
    # climb from 0, S declines — with carbon closing to machine precision through the segment (the
    # sugar carbon parks in melanoidin). Sealed (o2 defaults to 0) — the browning needs no oxygen.
    schema = wine_schema()
    y0 = _caramel_wine(schema, s=120.0, t=303.15)  # 30 °C sweet wine
    ps = ProcessSet(schema, [Caramelization()], strict=True)
    traj = simulate(ps, params=caramel_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    mel = traj.series("melanoidin")
    a420 = traj.series("A420")
    s = traj.series("S")
    assert mel[-1] > mel[0] == 0.0  # melanoidin accumulates from 0
    assert a420[-1] > a420[0] == 0.0  # A420 browning index climbs from 0
    assert s[-1] < s[0]  # residual sugar declines (consumed into melanoidin)
    assert_nonnegative(traj, ("S", "melanoidin", "A420"), atol=1e-9)
    f_c = caramel_store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_caramelization_tier_floored_at_speculative(caramel_store):
    # Speculative in FORM (Tier-3 frontier): every pool it writes is speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [Caramelization()])
    for pool in _CARAMEL_TOUCHES:
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, caramel_store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# MaillardBrowning (decision D-89) — the WINE-ONLY, NON-oxidative amino-acid-incorporating THERMAL
# browning: the N-bearing browning branch D-88's sugar-only Caramelization deferred. Residual SUGAR
# +
# AMINO ACID brown by HEAT (no O₂) to a NITROGEN-bearing maillard_melanoidin polymer, raising the
# SAME
# A420 index D-74/D-88 accumulate. It consumes core S AND amino_acids and RETAINS the amino-acid
# nitrogen in the polymer (the deaminating branch is D-87), so maillard_melanoidin is the FIRST
# non-biomass, non-arginine species on total_nitrogen. Closure is by SIZING both draws to the
# melanoidin formed (its fixed C:N stand-in), so carbon AND nitrogen close for any formula given the
# sign of the denominator (c_m − n_m·c(arg)/n(arg) > 0). These tests pin the closed form, carbon AND
# nitrogen closure per-RHS, the load-bearing denominator sign, the amino-acid HARD gate
# (isolability),
# the sugar SOFT gate, O₂-independence, the availability-gate saturation, warmer-faster, the
# wine-only
# no-op on beer, the integrated sweet browning + dual-ledger closure, and the speculative tier
# floor.

#: Browning retains amino-acid NITROGEN without caring which molecule carried it, so it draws the
#: IDENTITY-AGNOSTIC pools (D-100) — never a precursor. That is what makes it immune to the Ehrlich
#: re-route, which is the starvation channel D-100 exists to close.
_MAILLARD_BROWNING_TOUCHES = {
    "S",
    "amino_acids",
    "amino_acids_generic",
    "maillard_melanoidin",
    "A420",
}
#: The two extremes of the drawn blend's mass C:N (D-100). The realised ratio is a
#: nitrogen-weighted average of these, so ANY pool composition lands between them.
_BLEND_CN_BOUNDS = (_AA_C / _AA_N, _GENERIC_C / _GENERIC_N)


def _mb_denom(ratio: float) -> float:
    """The sized-draw denominator ``c_m − n_m·R`` at a given blend C:N ``R``.

    Must be > 0 or ``r_m`` flips sign and the Process silently CREATES sugar — and no conservation
    test would catch it, because closure holds for either sign (the advisor's must-check trap,
    D-89). D-100 makes ``R`` state-dependent (it was arginine's fixed 1.29), so the check below
    asserts the sign across the whole reachable range rather than at one species.
    """
    return _MAILLARD_MELANOIDIN_C - _MAILLARD_MELANOIDIN_N * ratio


def _maillard_browning_wine(
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    aa: float = 0.3,
    s: float = 80.0,
    t: float = 298.15,
    **kw,
) -> FloatArray:
    """A finished, SEALED (o2 = 0) SWEET wine at the start of aging: residual sugar ``s`` + dosed
    amino acids at MUST-SPECTRUM composition (D-100), warm — the regime where amino-acid Maillard
    browning runs. Extra pools via kw."""
    y = _aged_wine(schema, ester=0.0, t=t)
    seed_amino_acids(y, schema, params, aa)
    y[schema.slice("S")] = s
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def _maillard_browning_closed_form(
    schema: StateSchema, params: dict[str, float], y: FloatArray, t: float
) -> dict[str, float]:
    """The Process's own sized-draw algebra, recomputed independently for the assertions."""
    # The identity-agnostic substrate, gated on ITS combined must-spectrum share (D-100). At
    # must-spectrum composition this equals the pre-split lumped gate exactly — the reduction
    # property — so the D-89 algebra is reproduced here rather than reinvented.
    assimilable = sum(float(y[schema.slice(spec.pool)][0]) for spec in ASSIMILABLE_SPECS)
    share = sum(params[spec.fraction_param] for spec in ASSIMILABLE_SPECS)
    s_total = float(y[schema.slice("S")].sum())
    gate = assimilable / (params["K_amino_acids"] * share + assimilable)
    f_t = arrhenius_factor(t, params["E_a_maillard_browning"], params["T_ref"])
    r_sugar = params["k_maillard_browning"] * f_t * s_total * gate
    ratio = assimilable_carbon_per_nitrogen(y, schema)  # the drawn blend's mass C:N
    mel_rate = r_sugar * _GLUCOSE_C / _mb_denom(ratio)
    nitrogen = mel_rate * _MAILLARD_MELANOIDIN_N
    # The nitrogen demand splits across {arginine, generic} by the nitrogen each pool holds.
    held = {
        spec.pool: float(y[schema.slice(spec.pool)][0]) * nitrogen_mass_fraction(spec.species)
        for spec in ASSIMILABLE_SPECS
    }
    total_n = sum(held.values())
    out = {
        "S": -r_sugar,
        "maillard_melanoidin": mel_rate,
        "A420": params["y_a420_per_maillard_melanoidin"] * mel_rate,
    }
    for spec in ASSIMILABLE_SPECS:
        out[spec.pool] = (
            -(held[spec.pool] / total_n) * nitrogen / nitrogen_mass_fraction(spec.species)
        )
    return out


def test_maillard_browning_metadata():
    p = MaillardBrowning()
    assert p.name == "maillard_browning"
    assert p.tier is Tier.SPECULATIVE
    # Consumes core S + shared amino_acids into the on-ledger N-bearing maillard_melanoidin park +
    # raises the shared A420. Touches those four and nothing else — NO o2 (the whole point), and NO
    # CO2/N (all carbon+nitrogen retained in the polymer; the deaminating/decarboxylating branch is
    # MaillardStrecker, D-87).
    assert set(p.touches) == _MAILLARD_BROWNING_TOUCHES
    assert "o2" not in p.touches
    assert "CO2" not in p.touches
    assert "N" not in p.touches
    assert set(p.reads) == {
        "k_maillard_browning",
        "E_a_maillard_browning",
        "y_a420_per_maillard_melanoidin",
        "K_amino_acids",
        "T_ref",
        *(spec.fraction_param for spec in ASSIMILABLE_SPECS),
    }


def test_maillard_browning_denominator_is_positive_with_margin():
    # THE load-bearing sign check (advisor's must-verify): the sized-draw denominator must be
    # comfortably positive, else r_m flips sign and the Process CREATES sugar with no conservation
    # test catching it (closure holds either sign). D-100 makes the drawn blend's C:N
    # state-dependent, so the check is now over the WHOLE reachable range — arginine-only (1.29) to
    # generic-only (2.14) — which is strictly stronger than the old single-species assertion. The
    # C-rich melanoidin (C:N ≈ 8:1) clears even the worst case by ~3.7×.
    for ratio in _BLEND_CN_BOUNDS:
        assert _mb_denom(ratio) > 0.0
        assert ratio < _MAILLARD_MELANOIDIN_C / _MAILLARD_MELANOIDIN_N
    # The margin narrows with D-100 (glutamine is less N-rich than arginine) but stays healthy:
    # 0.81·c_m at arginine's ratio, 0.69·c_m at glutamine's — the worst case the blend can reach.
    assert 0.65 < _mb_denom(max(_BLEND_CN_BOUNDS)) / _MAILLARD_MELANOIDIN_C < 0.72
    assert 0.77 < _mb_denom(min(_BLEND_CN_BOUNDS)) / _MAILLARD_MELANOIDIN_C < 0.85


def test_maillard_browning_matches_closed_form(caramel_params):
    schema = wine_schema()
    y = _maillard_browning_wine(
        schema, caramel_params, aa=0.3, s=80.0, t=298.15
    )  # off T_ref so Arrhenius bites
    d = MaillardBrowning().derivatives(0.0, y, schema, caramel_params)
    cf = _maillard_browning_closed_form(schema, caramel_params, y, 298.15)
    assert cf["maillard_melanoidin"] > 0.0  # live (guards against a vacuous pass)
    assert schema.get(d, "S") == pytest.approx(cf["S"])
    assert schema.get(d, "amino_acids") == pytest.approx(cf["amino_acids"])
    assert schema.get(d, "maillard_melanoidin") == pytest.approx(cf["maillard_melanoidin"])
    assert schema.get(d, "A420") == pytest.approx(cf["A420"])
    # Touches nothing else — no o2, no CO2, no N, no aroma pools, no E.
    for var in ("o2", "CO2", "N", "E", "isoamyl_acetate", "acetaldehyde", "sotolon", "melanoidin"):
        assert schema.get(d, var) == 0.0


def test_maillard_browning_carbon_closes_per_rhs(caramel_params):
    # CARBON closes to machine precision: the carbon leaving S (sugar) + BOTH identity-agnostic
    # pools, each at its own carbon fraction (D-100), equals the carbon entering the
    # maillard_melanoidin park. A420 is an optical index (off every ledger), so it carries none.
    # This closes only because the melanoidin was sized with the SAME blend C:N the draw realises
    # — the reason `assimilable_carbon_per_nitrogen` and `draw_assimilable_nitrogen` share a
    # split rule rather than each computing their own.
    schema = wine_schema()
    d = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params), schema, caramel_params
    )
    carbon_residual = (
        schema.get(d, "S") * _GLUCOSE_C
        + sum(
            schema.get(d, spec.pool) * carbon_mass_fraction(spec.species)
            for spec in ASSIMILABLE_SPECS
        )
        + schema.get(d, "maillard_melanoidin") * _MAILLARD_MELANOIDIN_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)


def test_maillard_browning_nitrogen_closes_per_rhs(caramel_params):
    # NITROGEN closes to machine precision — the D-89 novelty: the nitrogen leaving the
    # identity-agnostic pools is RETAINED in the maillard_melanoidin polymer (NOT deaminated to N —
    # that is D-87's branch), so the drawn nitrogen exactly cancels maillard_melanoidin↑·n(mm).
    # maillard_melanoidin is the FIRST non-biomass, non-amino-acid species on total_nitrogen.
    schema = wine_schema()
    d = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params), schema, caramel_params
    )
    nitrogen_residual = (
        sum(
            schema.get(d, spec.pool) * nitrogen_mass_fraction(spec.species)
            for spec in ASSIMILABLE_SPECS
        )
        + schema.get(d, "maillard_melanoidin") * _MAILLARD_MELANOIDIN_N
    )
    assert nitrogen_residual == pytest.approx(0.0, abs=1e-18)
    assert schema.get(d, "N") == 0.0  # nitrogen RETAINED in the polymer, not refunded to N


def test_maillard_browning_inert_without_amino_acids(caramel_params):
    # HARD amino-acid gate — the isolability guarantee: an undosed wine (aa = 0) is byte-for-byte
    # the
    # case without this Process, even with abundant residual sugar (unlike the SOFT sugar gate).
    # Also
    # absorbs the aa < 0 solver undershoot.
    schema = wine_schema()
    ps = ProcessSet(schema, [MaillardBrowning()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(
            0.0, _maillard_browning_wine(schema, caramel_params, aa=0.0, s=120.0), caramel_params
        ),
        schema.zeros(),
    )
    assert np.array_equal(
        MaillardBrowning().derivatives(
            0.0,
            _maillard_browning_wine(schema, caramel_params, aa=-1e-6, s=120.0),
            schema,
            caramel_params,
        ),
        schema.zeros(),
    )


def test_maillard_browning_soft_sugar_gate(caramel_params):
    # SOFT sugar driver: a dry wine (S = 0) makes ~none (every standard dry aging run ferments to
    # S ≈ 0 before begin_aging), and the S < 0 solver undershoot is absorbed — but this is a driver,
    # not the isolability gate (that is the amino-acid HARD gate above).
    schema = wine_schema()
    assert np.array_equal(
        MaillardBrowning().derivatives(
            0.0, _maillard_browning_wine(schema, caramel_params, s=0.0), schema, caramel_params
        ),
        schema.zeros(),
    )
    assert np.array_equal(
        MaillardBrowning().derivatives(
            0.0, _maillard_browning_wine(schema, caramel_params, s=-1e-6), schema, caramel_params
        ),
        schema.zeros(),
    )


def test_maillard_browning_is_oxygen_independent(caramel_params):
    # NO o2 term at all: a sealed wine (o2 = 0) browns exactly as an oxygenated one — the whole
    # point
    # (the O₂-independent thermal mirror). o2 is untouched.
    schema = wine_schema()
    sealed = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params, o2=0.0), schema, caramel_params
    )
    oxygenated = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params, o2=0.05), schema, caramel_params
    )
    assert np.array_equal(sealed, oxygenated)
    assert schema.get(sealed, "o2") == 0.0  # never draws o2
    assert schema.get(sealed, "maillard_melanoidin") > 0.0  # non-vacuous: browning while sealed


def test_maillard_browning_availability_gate_saturates(caramel_params):
    # The amino-acid availability gate aa/(K+aa) is the same smooth-Monod shape as the Strecker
    # routes: at aa ≫ K the rate saturates (≈ proportional to sugar alone), at aa = K it is half.
    schema = wine_schema()
    k = caramel_params["K_amino_acids"]
    half = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params, aa=k), schema, caramel_params
    )
    saturated = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params, aa=1.0e6 * k), schema, caramel_params
    )
    # At aa = K the melanoidin rate is half its saturated value (gate = 0.5 vs ≈ 1.0).
    assert schema.get(half, "maillard_melanoidin") == pytest.approx(
        0.5 * schema.get(saturated, "maillard_melanoidin"), rel=1e-3
    )


def test_maillard_browning_rises_with_temperature(caramel_params):
    # Warmer browns faster — the sourced direction (Maillard browning is strongly thermal, why
    # Madeira estufagem develops nitrogenous-melanoidin colour so much faster than cellar aging).
    schema = wine_schema()
    cold = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params, t=283.15), schema, caramel_params
    )
    warm = MaillardBrowning().derivatives(
        0.0, _maillard_browning_wine(schema, caramel_params, t=303.15), schema, caramel_params
    )
    assert schema.get(warm, "maillard_melanoidin") > schema.get(cold, "maillard_melanoidin") > 0.0


def test_maillard_browning_is_wine_only_noop_on_beer(caramel_params):
    # Wine-only v1 (amino_acids + the maillard_melanoidin park are wine slots): a hard no-op on beer
    # even with residual wort sugar (the "maillard_melanoidin"/"amino_acids" not in schema guard).
    beer = beer_schema()
    yb = beer.zeros()
    yb[beer.slice("S")] = 60.0
    yb[beer.slice("T")] = 303.15
    assert np.array_equal(
        MaillardBrowning().derivatives(0.0, yb, beer, caramel_params), beer.zeros()
    )


def test_maillard_browning_browns_sweet_wine_closes_carbon_and_nitrogen(
    caramel_store, caramel_params
):
    # Integrated: a warm SWEET + amino-acid-dosed wine browns over an aging year —
    # maillard_melanoidin
    # + A420 climb from 0, S + amino_acids decline — with BOTH carbon AND nitrogen closing to
    # machine
    # precision through the segment (the sugar + amino-acid carbon and the amino-acid nitrogen park
    # in
    # maillard_melanoidin). Sealed (o2 defaults to 0) — the browning needs no oxygen.
    schema = wine_schema()
    y0 = _maillard_browning_wine(
        schema, caramel_params, aa=0.6, s=120.0, t=303.15
    )  # 30 °C sweet + amino acids
    ps = ProcessSet(schema, [MaillardBrowning()], strict=True)
    traj = simulate(ps, params=caramel_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    mel = traj.series("maillard_melanoidin")
    a420 = traj.series("A420")
    s = traj.series("S")
    aa = traj.series("amino_acids")
    assert mel[-1] > mel[0] == 0.0  # N-melanoidin accumulates from 0
    assert a420[-1] > a420[0] == 0.0  # A420 browning index climbs from 0
    assert s[-1] < s[0]  # residual sugar declines (consumed into melanoidin)
    assert aa[-1] < aa[0]  # amino acids consumed (retained in the N-bearing polymer)
    assert_nonnegative(traj, ("S", "amino_acids", "maillard_melanoidin", "A420"), atol=1e-9)
    f_c = caramel_store.value("biomass_C_fraction")
    f_n = caramel_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_maillard_browning_tier_floored_at_speculative(caramel_store):
    # Speculative in FORM (Tier-3 frontier): every pool it writes is speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [MaillardBrowning()])
    for pool in _MAILLARD_BROWNING_TOUCHES:
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, caramel_store.tier_map()) is Tier.SPECULATIVE


# =====================================================================================
# OakExtraction (decision D-77) — the WINE-ONLY, NON-oxidative barrel/chip aroma-extraction axis, a
# SEPARATE axis from the oxidative sub-axis (draws NO O₂). Four wood extractives — whiskey_lactone
# (coconut), vanillin (vanilla), guaiacol (smoky), eugenol (clove) — diffuse into the wine and rise
# toward per-compound SET-AND-HOLD ceilings (the cation_charge idiom, set by the add_oak verb):
# d(C_i)/dt = k_oak_extraction·f(T)·max(0, ceiling_i − C_i), the inverse of EsterHydrolysis's decay.
# OFF EVERY LEDGER (exogenous wood-derived, the iso_alpha precedent) — touches only the extracted
# slots, moves nothing conserved, needs no chemistry.py registration. These tests pin the closed
# form, the first-order-in-gap linearity, the ceiling ≤ 0 guard (undosed AND undershoot — the floor
# is 0, so max() alone is insufficient), the monotone approach + saturation over an integrated
# segment, the warmer-faster ordering, the wine-only no-op on beer, the off-every-ledger invariance
# (carbon AND mass AND nitrogen flat), and the speculative tier floor.

_OAK_COMPOUNDS = ("whiskey_lactone", "vanillin", "guaiacol", "eugenol", "furaneol")


@pytest.fixture
def oak_store():
    # Wine params + the oak.yaml extraction constants (k_oak_extraction, the weak diffusion
    # E_a_oak_extraction, T_ref), the shared_files the D-77 compile seam wires. (The plain ``store``
    # fixture omits oak.yaml; OakExtraction reads its own rate/E_a from here.)
    return load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "oak.yaml"
    )


@pytest.fixture
def oak_params(oak_store):
    return oak_store.resolve()


def _oak_wine(
    schema: StateSchema,
    *,
    ceilings: dict[str, float] | None = None,
    t: float = 293.15,
    **kw,
) -> FloatArray:
    """A finished, racked wine at the start of oak aging: the four ceiling slots pre-set (as the
    add_oak verb would, oak_gpl × toast yield), the extracted pools starting at 0. ``ceilings`` maps
    a compound name to its ceiling g/L; any extra pool via kwargs."""
    y = _aged_wine(schema, ester=0.0, t=t)
    for compound, ceiling in (ceilings or {}).items():
        y[schema.slice(f"{compound}_ceiling")] = ceiling
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def test_oak_metadata():
    p = OakExtraction()
    assert p.name == "oak_extraction"
    assert p.tier is Tier.SPECULATIVE
    # Writes ONLY the six extracted-compound slots — the five aroma extractives (four D-77 +
    # furaneol/caramel D-94) plus the ellagitannin taste extractive (D-78). The ceilings are read,
    # never written (a set-and-hold constant the add_oak verb owns). Off every ledger, none moves.
    assert set(p.touches) == {
        "whiskey_lactone",
        "vanillin",
        "guaiacol",
        "eugenol",
        "furaneol",
        "ellagitannin",
    }
    # Only its own shared rate/E_a + T_ref; the per-compound ceilings ride in STATE, not params.
    assert set(p.reads) == {"k_oak_extraction", "E_a_oak_extraction", "T_ref"}


def test_oak_matches_closed_form(oak_params):
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    # Distinct ceilings and a partial pre-fill so each compound's gap is unambiguous.
    ceilings = {
        "whiskey_lactone": 8.0e-5,
        "vanillin": 2.0e-4,
        "guaiacol": 1.5e-5,
        "eugenol": 8.0e-6,
    }
    y = _oak_wine(schema, ceilings=ceilings, t=t)
    y[schema.slice("vanillin")] = 5.0e-5  # partly extracted already ⇒ gap = ceiling − conc
    d = OakExtraction().derivatives(0.0, y, schema, oak_params)

    f_t = arrhenius_factor(t, oak_params["E_a_oak_extraction"], oak_params["T_ref"])
    k = oak_params["k_oak_extraction"]
    for compound, ceiling in ceilings.items():
        conc = 5.0e-5 if compound == "vanillin" else 0.0
        assert schema.get(d, compound) == pytest.approx(k * f_t * (ceiling - conc))
    # Extraction touches nothing else — no O₂ (a separate axis), no sugar/ethanol/CO2, no ceilings.
    for var in ("X", "S", "E", "N", "CO2", "o2", "A420", "acetaldehyde"):
        assert schema.get(d, var) == 0.0
    for compound in _OAK_COMPOUNDS:
        assert (
            schema.get(d, f"{compound}_ceiling") == 0.0
        )  # ceilings are set-and-hold (never moved)


def test_oak_is_first_order_in_the_gap(oak_params):
    # First-order approach from below: twice the gap (ceiling − conc) ⇒ twice the instantaneous
    # extraction rate. This linearity is what makes each compound SATURATE at its ceiling.
    schema = wine_schema()
    near = OakExtraction().derivatives(
        0.0, _oak_wine(schema, ceilings={"vanillin": 2.0e-4}, vanillin=1.6e-4), schema, oak_params
    )
    far = OakExtraction().derivatives(
        0.0, _oak_wine(schema, ceilings={"vanillin": 2.0e-4}, vanillin=1.2e-4), schema, oak_params
    )
    # gap_far = 8e-5 = 2 × gap_near (4e-5), so the rate is doubled.
    assert schema.get(far, "vanillin") == pytest.approx(2.0 * schema.get(near, "vanillin"))
    assert schema.get(far, "vanillin") > schema.get(near, "vanillin") > 0.0  # both extracting


def test_oak_inert_without_dose_and_on_undershoot(oak_params):
    # No oak dosed (every ceiling 0) ⇒ byte-for-byte inert: a begin_aging run with no add_oak is the
    # case without oak. And the EXPLICIT ceiling ≤ 0 guard is load-bearing — the floor is 0 (unlike
    # the ester floor is > 0), so a solver undershoot conc = −ε must NOT flip
    # max(0, ceiling − conc) into
    # spurious extraction. Both the undosed and the undershoot case return byte-for-byte zero.
    schema = wine_schema()
    ps = ProcessSet(schema, [OakExtraction()], strict=True)
    assert np.array_equal(ps.total_derivatives(0.0, _oak_wine(schema), oak_params), schema.zeros())
    # Undershoot on an UNDOSED compound (ceiling 0, conc = −ε): the ceiling ≤ 0 guard blocks it.
    undershoot = _oak_wine(schema)
    undershoot[schema.slice("guaiacol")] = -1e-9
    assert np.array_equal(
        OakExtraction().derivatives(0.0, undershoot, schema, oak_params), schema.zeros()
    )


def test_oak_stops_at_the_ceiling(oak_params):
    # At/above the ceiling the gap ≤ 0 ⇒ no further extraction (monotone rise, never overshoots).
    schema = wine_schema()
    at_ceiling = _oak_wine(schema, ceilings={"vanillin": 2.0e-4}, vanillin=2.0e-4)
    d = OakExtraction().derivatives(0.0, at_ceiling, schema, oak_params)
    assert schema.get(d, "vanillin") == 0.0


def test_oak_rises_with_temperature(oak_params):
    # The sourced (but deliberately WEAK — diffusion-limited) ordering E_a_oak_extraction > 0: a
    # warmer barrel extracts a little faster. Still monotone in T even at the low diffusion E_a.
    schema = wine_schema()
    cold = OakExtraction().derivatives(
        0.0, _oak_wine(schema, ceilings={"vanillin": 2.0e-4}, t=283.15), schema, oak_params
    )
    warm = OakExtraction().derivatives(
        0.0, _oak_wine(schema, ceilings={"vanillin": 2.0e-4}, t=303.15), schema, oak_params
    )
    assert schema.get(warm, "vanillin") > schema.get(cold, "vanillin") > 0.0


def test_oak_is_wine_only_noop_on_beer(oak_params):
    # WINE-ONLY: the oak slots are appended to wine_schema (beer has none), so the
    # "whiskey_lactone" not in schema guard makes this a hard no-op on beer even if wired there.
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    assert np.array_equal(OakExtraction().derivatives(0.0, yb, beer, oak_params), beer.zeros())


def test_integrated_oak_saturates_and_moves_nothing_conserved(oak_params, oak_store):
    # Run a long oak-aging segment (racked, dry wine) with ONLY OakExtraction and the four ceilings
    # pre-set. Over the span every extractive rises MONOTONICALLY to a PLATEAU at its ceiling
    # (saturating — the barrel aroma builds then levels), and — because the extractives + ceilings
    # are ALL off every ledger — this Process moves NOTHING conserved: total_carbon, total_mass AND
    # total_nitrogen are all exactly flat (the cleanest possible aging Process — the iso_alpha/A420
    # off-ledger invariance, extended to carbon AND mass AND nitrogen together).
    schema = wine_schema()
    ps = ProcessSet(schema, [OakExtraction()], strict=True)
    ceilings = {
        "whiskey_lactone": 8.0e-5,
        "vanillin": 2.0e-4,
        "guaiacol": 6.0e-5,
        "eugenol": 2.5e-5,
        "furaneol": 4.8e-5,  # D-94: the caramel furanone extracts + is off every ledger too
    }
    y0 = _oak_wine(schema, ceilings=ceilings, t=298.15)
    traj = simulate(ps, params=oak_params, y0=y0, t_span=(0.0, 24.0 * 365.0))  # ~1 year
    assert traj.success, traj.message

    for compound, ceiling in ceilings.items():
        series = np.asarray(traj.series(compound), dtype=float)
        # Rose from 0, monotonically, toward (never past) the ceiling.
        assert series[0] == 0.0
        assert np.all(np.diff(series) >= -1e-18)  # monotone rise
        assert series[-1] <= ceiling + 1e-15  # never overshoots
        assert (
            series[-1] >= 0.5 * ceiling
        )  # reached most of the ceiling over a year (sanity anchor)
        # The set-and-hold ceiling never moved.
        assert float(traj.series(f"{compound}_ceiling")[-1]) == pytest.approx(ceiling)
    assert_nonnegative(traj, _OAK_COMPOUNDS, atol=1e-15)
    # ALL THREE ledgers flat — oak extraction moves nothing conserved (exogenous wood-derived mass,
    # off every ledger). X=0 throughout, so the biomass terms are inert.
    f_c = oak_store.value("biomass_C_fraction")
    f_n = oak_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_oak_tier_floored_at_speculative(oak_store):
    # Speculative in FORM (Tier-3 frontier): each extracted pool it writes is speculative even pre-
    # the (speculative) oak parameter tiers cap it. Non-vacuous across all four extractives.
    schema = wine_schema()
    ps = ProcessSet(schema, [OakExtraction()])
    for pool in _OAK_COMPOUNDS:
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, oak_store.tier_map()) is Tier.SPECULATIVE


def test_oak_also_extracts_ellagitannin(oak_params):
    # D-78: OakExtraction now extracts a FIFTH pool, ellagitannin (the taste/O₂-scavenging tannin),
    # by the IDENTICAL diffusion-to-a-ceiling form — d(ellag)/dt = k·f(T)·max(0, ceiling − ellag).
    schema = wine_schema()
    t = 298.15
    y = _oak_wine(schema, ceilings={"ellagitannin": 0.1}, t=t)
    y[schema.slice("ellagitannin")] = 0.04  # partly extracted ⇒ gap = 0.1 − 0.04
    d = OakExtraction().derivatives(0.0, y, schema, oak_params)
    f_t = arrhenius_factor(t, oak_params["E_a_oak_extraction"], oak_params["T_ref"])
    k = oak_params["k_oak_extraction"]
    assert schema.get(d, "ellagitannin") == pytest.approx(k * f_t * (0.1 - 0.04))
    # Extraction draws NO O₂ (a diffusion process) and never moves the ceiling.
    assert schema.get(d, "o2") == 0.0
    assert schema.get(d, "ellagitannin_ceiling") == 0.0


def test_oak_extracts_furaneol_the_caramel_furanone(oak_params):
    # D-94: OakExtraction extracts a SIXTH pool, furaneol (the caramel/toffee furanone — the caramel
    # note D-93 deferred), by the IDENTICAL diffusion-to-a-ceiling form as the other aroma four.
    schema = wine_schema()
    t = 298.15
    y = _oak_wine(schema, ceilings={"furaneol": 1.0e-4}, t=t)
    y[schema.slice("furaneol")] = 3.0e-5  # partly extracted ⇒ gap = 1.0e-4 − 3.0e-5
    d = OakExtraction().derivatives(0.0, y, schema, oak_params)
    f_t = arrhenius_factor(t, oak_params["E_a_oak_extraction"], oak_params["T_ref"])
    k = oak_params["k_oak_extraction"]
    assert schema.get(d, "furaneol") == pytest.approx(k * f_t * (1.0e-4 - 3.0e-5))
    # Off every ledger like the other aroma pools: draws no O₂ and never moves the ceiling.
    assert schema.get(d, "o2") == 0.0
    assert schema.get(d, "furaneol_ceiling") == 0.0


@pytest.fixture
def caramel_oak_store():
    # Wine + thermal.yaml (Caramelization: k/E_a/y_a420) + oak.yaml (OakExtraction rate/E_a/yields):
    # the two Processes D-94's collision thesis is about run from one combined store.
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "thermal.yaml",
        default_data_dir() / "oak.yaml",
    )


def test_furaneol_and_caramelization_coexist_without_collision(caramel_oak_store):
    # D-94's load-bearing thesis: the caramel AROMA (furaneol, on the oak axis — off every ledger)
    # and the caramel COLOUR (D-88 melanoidin, ON total_carbon) are the SAME browning chemistry read
    # two ways, and they DO NOT collide. Run a SWEET wine (residual sugar) aging in oak with BOTH
    # Caramelization and OakExtraction active: melanoidin forms (on-ledger, from core S) AND the
    # furaneol pool extracts toward its ceiling (off-ledger, from the oak/spirit ceiling). Carbon
    # EXACTLY — the sugar→melanoidin transfer is carbon-exact and furaneol adds NOTHING to it,
    # because the two never share a conserved pool (the collision the D-93 deferral feared is
    # dissolved by putting caramel AROMA on the off-ledger oak axis, not an on-ledger co-product).
    params = caramel_oak_store.resolve()
    schema = wine_schema()
    y0 = _caramel_wine(schema, s=100.0, t=308.15)  # warm sweet wine ⇒ caramelization runs (Madeira)
    y0[schema.slice("furaneol_ceiling")] = 1.0e-4  # ~100 µg/L ceiling (as add_oak wood + bourbon)
    ps = ProcessSet(schema, [Caramelization(), OakExtraction()], strict=True)
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))  # ~1 year
    assert traj.success, traj.message

    # BOTH ran: melanoidin (on-ledger caramelization COLOUR) AND furaneol (off-ledger AROMA) rose.
    assert float(traj.series("melanoidin")[-1]) > 0.0  # caramelization browned the sugar
    furaneol_end = float(traj.series("furaneol")[-1])
    assert 0.0 < furaneol_end <= 1.0e-4 + 1e-15  # extracted from 0 toward (not past) its ceiling
    assert float(traj.series("S")[-1]) < 100.0  # residual sugar was consumed into melanoidin

    # NO COLLISION: total_carbon closes to machine precision though both Processes ran. Furaneol's
    # extraction moves nothing conserved (off every ledger, the iso_alpha precedent), so it cannot
    # perturb the D-88 sugar→melanoidin carbon closure — the whole reason caramel AROMA lives on the
    # off-ledger oak axis rather than as an on-ledger caramelization co-product (D-94).
    f_c = caramel_oak_store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


# =====================================================================================
# EllagitanninOxidation (decision D-78) — the WINE-ONLY oak-tannin O₂-scavenging sink, the BRIDGE
# from the D-77 oak extractive axis to the O₂ sub-axis. Oak's hydrolysable tannin (the ellagitannin
# pool OakExtraction fills) is a sacrificial antioxidant: dissolved O₂ oxidises it (bilinear
# [o2]·[ellagitannin], the SulfiteOxidation form), CONSUMING it at a mass-based yield y_ellag_per_o2
# (g ellag / g O₂ — no fake molar mass for the lumped macromolecule). The EMERGENT SPINE is oak
# PROTECTION: an oaked + oxygenated wine browns LESS (lower A420) and makes LESS oxidative
# acetaldehyde than an un-oaked wine at the same O₂ dose (the D-72 "SO₂ protects" threshold with a
# RENEWABLE buffer). Substrate-gated on ellagitannin ⇒ adds on top with NO re-baseline of the
# k_ethanol_oxidation + k_browning_base = 5.0e-4 anchor. Off every ledger (both slots unweighted),
# so it moves nothing conserved. These tests pin the closed form, the bilinearity, the
# reaction-scale temperature ordering, the doubly-substrate-gated inertness (KeyError-safe without
# oak.yaml), the wine-only no-op on beer, THE PROTECTION SPINE (partial, not total), the
# sacrificial-consumption softening (astringency_series), the off-every-ledger invariance, and the
# speculative tier floor.


@pytest.fixture
def ellag_store():
    # Wine params + aging.yaml (the O₂ sinks browning/ethanol oxidation) + oak.yaml (the
    # ellagitannin extraction yields AND the D-78 scavenging rate/E_a/yield) — the shared_files the
    # compile seam wires. EllagitanninOxidation + OakExtraction + PhenolicBrowning +
    # OxidativeAcetaldehyde all resolve from this combined store.
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "aging.yaml",
        default_data_dir() / "oak.yaml",
    )


@pytest.fixture
def ellag_params(ellag_store):
    return ellag_store.resolve()


def test_ellagitannin_oxidation_metadata():
    p = EllagitanninOxidation()
    assert p.name == "ellagitannin_oxidation"
    # Speculative: the aging axis is the Tier-3 frontier.
    assert p.tier is Tier.SPECULATIVE
    # Touches ONLY o2 (consumed) and ellagitannin (oxidised) — both off every ledger, so nothing
    # conserved moves (the SulfiteOxidation precedent). No carbon borrow.
    assert set(p.touches) == {"o2", "ellagitannin"}
    assert set(p.reads) == {
        "k_ellagitannin_oxidation",
        "E_a_ellagitannin_oxidation",
        "y_ellag_per_o2",
        "T_ref",
    }


def test_ellagitannin_oxidation_closed_form(ellag_params):
    # d(o2)/dt = −k·f(T)·[o2]·[ellag] (bilinear, the SulfiteOxidation form); d(ellag)/dt = −y·r_o2
    # (mass-based consumption). Verify both exactly, and that nothing else moves.
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    o2, ellag = 0.03, 0.08
    y = _aged_wine(schema, ester=0.0, t=t, o2=o2, ellagitannin=ellag)
    d = EllagitanninOxidation().derivatives(0.0, y, schema, ellag_params)

    f_t = arrhenius_factor(t, ellag_params["E_a_ellagitannin_oxidation"], ellag_params["T_ref"])
    k = ellag_params["k_ellagitannin_oxidation"]
    y_ellag = ellag_params["y_ellag_per_o2"]
    r_o2 = k * f_t * o2 * ellag
    assert schema.get(d, "o2") == pytest.approx(-r_o2)
    assert schema.get(d, "ellagitannin") == pytest.approx(-y_ellag * r_o2)
    # Off every ledger — no carbon borrow (unlike OxidativeAcetaldehyde's E→acetaldehyde), no A420,
    # no other pool touched.
    for var in ("X", "S", "E", "N", "CO2", "A420", "acetaldehyde", "so2_total"):
        assert schema.get(d, var) == 0.0


def test_ellagitannin_oxidation_is_bilinear(ellag_params):
    # Bilinear in BOTH drivers: doubling o2 OR doubling ellagitannin doubles the O₂-scavenging rate.
    schema = wine_schema()
    base = EllagitanninOxidation().derivatives(
        0.0, _aged_wine(schema, ester=0.0, o2=0.02, ellagitannin=0.05), schema, ellag_params
    )
    twice_o2 = EllagitanninOxidation().derivatives(
        0.0, _aged_wine(schema, ester=0.0, o2=0.04, ellagitannin=0.05), schema, ellag_params
    )
    twice_ellag = EllagitanninOxidation().derivatives(
        0.0, _aged_wine(schema, ester=0.0, o2=0.02, ellagitannin=0.10), schema, ellag_params
    )
    assert schema.get(twice_o2, "o2") == pytest.approx(2.0 * schema.get(base, "o2"))
    assert schema.get(twice_ellag, "o2") == pytest.approx(2.0 * schema.get(base, "o2"))
    assert schema.get(base, "o2") < 0.0  # actually scavenging


def test_ellagitannin_oxidation_inert_without_o2_or_tannin(ellag_params):
    # Doubly substrate-gated: no O₂ OR no ellagitannin ⇒ byte-for-byte zero. A reductive (no
    # add_oxygen) or an un-oaked aging is exactly the case without this Process (isolability #3).
    schema = wine_schema()
    p = EllagitanninOxidation()
    no_o2 = _aged_wine(schema, ester=0.0, o2=0.0, ellagitannin=0.08)
    no_tannin = _aged_wine(schema, ester=0.0, o2=0.03, ellagitannin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_o2, schema, ellag_params), schema.zeros())
    assert np.array_equal(p.derivatives(0.0, no_tannin, schema, ellag_params), schema.zeros())
    # <= 0 also absorbs solver undershoot (a spurious −ε in either driver ⇒ no draw).
    undershoot = _aged_wine(schema, ester=0.0, o2=-1e-9, ellagitannin=0.08)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, ellag_params), schema.zeros())


def test_ellagitannin_oxidation_gate_before_params_is_keyerror_safe(params):
    # Gate on the ellagitannin STATE before reading any oak param, so an enabled-but-undosed Process
    # never KeyErrors when oak.yaml is ABSENT (the ``params`` fixture is wine+aging only, no
    # k_ellagitannin_oxidation). An un-oaked wine (ellag=0) returns zero without touching oak
    # params.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, o2=0.03, ellagitannin=0.0)
    assert np.array_equal(
        EllagitanninOxidation().derivatives(0.0, y, schema, params), schema.zeros()
    )


def test_ellagitannin_oxidation_rises_with_temperature(ellag_params):
    # REACTION-scale E_a_ellagitannin_oxidation > 0 (its OWN param, ~50 kJ/mol — distinct from the
    # WEAK diffusion E_a_oak_extraction that governs the tannin's extraction): warmer scavenges
    # faster.
    schema = wine_schema()
    p = EllagitanninOxidation()
    cold = p.derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=283.15, o2=0.03, ellagitannin=0.08),
        schema,
        ellag_params,
    )
    warm = p.derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=303.15, o2=0.03, ellagitannin=0.08),
        schema,
        ellag_params,
    )
    # More negative (faster scavenging) when warmer.
    assert schema.get(warm, "o2") < schema.get(cold, "o2") < 0.0


def test_ellagitannin_oxidation_wine_only_noop_on_beer(ellag_params):
    # WINE-ONLY: the ellagitannin slot is appended to wine_schema (beer has none), so the
    # "ellagitannin" not in schema guard makes this a hard no-op on beer even if wired there.
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    yb[beer.slice("o2")] = 0.03  # o2 present, but no ellagitannin slot exists
    assert np.array_equal(
        EllagitanninOxidation().derivatives(0.0, yb, beer, ellag_params), beer.zeros()
    )


def test_oak_protects_against_oxidation(ellag_params):
    # THE D-78 SPINE. Two identical oxygenated aging runs — same O₂ dose, full oxidative process set
    # — differing ONLY in whether the wine is oaked (ellagitannin present). The oaked wine browns
    # LESS (lower A420) and makes LESS oxidative acetaldehyde, because the tannin scavenges its
    # share of the O₂ (oak protection, the SO₂-protection analogue with a renewable buffer).
    # Suppression is PARTIAL (the oaked wine still shows SOME browning/acetaldehyde — banded so
    # ellagitannin never monopolizes the O₂). And the tannin is CONSUMED (its pool declines) — the
    # sacrificial softening.
    schema = wine_schema()
    processes = [
        OakExtraction(),
        EllagitanninOxidation(),
        PhenolicBrowning(),
        OxidativeAcetaldehyde(),
    ]
    ps = ProcessSet(schema, processes, strict=True)
    o2_dose = 0.04  # ~40 mg/L cumulative O₂ exposure
    span = (0.0, 24.0 * 365.0)  # ~1 year

    # Oaked: an ellagitannin charge at its ceiling (renewable — OakExtraction re-supplies below it).
    oaked0 = _aged_wine(schema, ester=0.0, o2=o2_dose, ellagitannin=0.1, ellagitannin_ceiling=0.1)
    oaked = simulate(ps, params=ellag_params, y0=oaked0, t_span=span)
    # Un-oaked: identical but no tannin (and no ceiling), same O₂ dose.
    unoaked0 = _aged_wine(schema, ester=0.0, o2=o2_dose, ellagitannin=0.0, ellagitannin_ceiling=0.0)
    unoaked = simulate(ps, params=ellag_params, y0=unoaked0, t_span=span)
    assert oaked.success and unoaked.success

    a420_oaked = float(oaked.series("A420")[-1])
    a420_unoaked = float(unoaked.series("A420")[-1])
    acet_oaked = float(oaked.series("acetaldehyde")[-1])
    acet_unoaked = float(unoaked.series("acetaldehyde")[-1])

    # PROTECTION: oak lowers BOTH the browning index and the oxidative acetaldehyde.
    assert a420_oaked < a420_unoaked
    assert acet_oaked < acet_unoaked
    # PARTIAL, not total (ellagitannin does not monopolize the O₂ — banded so oaked wines still
    # oxidise SOME): the oaked wine still browns and still makes acetaldehyde.
    assert a420_oaked > 0.0 and acet_oaked > 0.0
    # SACRIFICIAL: the tannin IS drawn on (its pool sits at/below the ceiling — the buffer is being
    # spent). NB with the wood re-supply this is only a weak witness (OakExtraction refills the pool
    # back toward the ceiling once the O₂ is gone, so it ends only fractionally below 0.1) — the
    # GENUINE monotone softening is isolated in test_ellagitannin_consumed_softens_astringency (no
    # re-supply). Here the physical truth is a BUFFERED tannin, not permanent softening.
    assert float(oaked.series("ellagitannin")[-1]) <= 0.1


def test_ellagitannin_consumed_softens_astringency(ellag_params):
    # WITHOUT the renewable re-supply (EllagitanninOxidation alone, no OakExtraction), the
    # sacrificial tannin is drawn down MONOTONICALLY under O₂ — so astringency_series (mg/L
    # ellagitannin) SOFTENS over the run. One directional contributor (oxidative consumption);
    # polymerization deferred.
    schema = wine_schema()
    ps = ProcessSet(schema, [EllagitanninOxidation()], strict=True)
    y0 = _aged_wine(schema, ester=0.0, o2=0.05, ellagitannin=0.1)
    traj = simulate(ps, params=ellag_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    astr = astringency_series(traj)
    # astringency_series is exactly ellagitannin (g/L) × 1000 (mg/L), IBU-exact (reads no
    # threshold).
    assert np.allclose(astr, np.asarray(traj.series("ellagitannin"), dtype=float) * 1000.0)
    # Softens: monotone non-increasing, and ends below the start (some tannin spent on the O₂).
    assert np.all(np.diff(astr) <= 1e-15)
    assert astr[-1] < astr[0]
    assert astr[-1] >= 0.0  # never negative (the o2 charge is finite, tannin only partly spent)


def test_ellagitannin_oxidation_moves_nothing_conserved(ellag_store, ellag_params):
    # Both o2 and ellagitannin are off every ledger (wood-derived / carbon-free), so oxidising the
    # tannin moves NOTHING conserved — total_carbon, total_mass AND total_nitrogen all exactly flat
    # (the SulfiteOxidation off-every-ledger invariance). X=0 throughout.
    schema = wine_schema()
    ps = ProcessSet(schema, [EllagitanninOxidation()], strict=True)
    y0 = _aged_wine(schema, ester=0.0, o2=0.05, ellagitannin=0.1)
    traj = simulate(ps, params=ellag_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    f_c = ellag_store.value("biomass_C_fraction")
    f_n = ellag_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_ellagitannin_oxidation_tier_floored_at_speculative(ellag_store):
    # Speculative in FORM (Tier-3 frontier) and capped by its speculative oak params (D-1). The
    # ellagitannin pool it writes is speculative both pre- and post-parameter-tier propagation.
    schema = wine_schema()
    ps = ProcessSet(schema, [EllagitanninOxidation()])
    assert ps.tier_of("ellagitannin") is Tier.SPECULATIVE
    assert ps.tier_of("ellagitannin", ellag_store.tier_map()) is Tier.SPECULATIVE


# -- D-79: TanninAnthocyaninCondensation — red-wine colour stabilization + astringency softening --
#
# The eighth aging Process, the SECOND non-oxidative one (after OakExtraction) and a THIRD separate
# axis: grape anthocyanin + grape tannin condense (bilinear) into stable polymeric pigment — the
# DOMINANT softening + colour-evolution mechanism D-77/D-78 deferred. These tests pin the
# closed-form
# derivative, the bilinearity, the OFF-EVERY-LEDGER invariance, the doubly-substrate-gated
# isolability
# (a white / no-tannin wine is byte-for-byte inert; no o2 term ⇒ oak- AND O₂-independent), the
# warmer-faster ordering, the wine-only no-op on beer, the speculative tier floor, and the three
# readouts (astringency softens, polymeric pigment rises, colour is retained).


@pytest.fixture
def poly_store():
    # Wine params + polymerization.yaml (k_polymerization, E_a_polymerization,
    # y_tannin_per_anthocyanin,
    # T_ref) — the shared_files the D-79 compile seam wires. TanninAnthocyaninCondensation reads its
    # own rate/E_a/yield from here. (The plain ``store`` fixture omits polymerization.yaml.)
    return load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "polymerization.yaml"
    )


@pytest.fixture
def poly_params(poly_store):
    return poly_store.resolve()


def test_polymerization_metadata():
    p = TanninAnthocyaninCondensation()
    assert p.name == "tannin_anthocyanin_condensation"
    # Speculative: the aging axis is the Tier-3 frontier.
    assert p.tier is Tier.SPECULATIVE
    # Touches the two grape pools it condenses PLUS the polymeric_pigment slot it fills (D-81
    # promotion) — all three off every ledger, so nothing conserved moves (the OakExtraction/
    # EllagitanninOxidation precedent). NO o2 (a non-oxidative grape axis).
    assert set(p.touches) == {"anthocyanin", "tannin", "polymeric_pigment"}
    assert "o2" not in p.touches  # oak- AND O₂-independent (the D-79 crux)
    assert set(p.reads) == {
        "k_polymerization",
        "E_a_polymerization",
        "y_tannin_per_anthocyanin",
        "T_ref",
    }


def test_polymerization_closed_form(poly_params):
    # r = k·f(T)·[anthocyanin]·[tannin] (bilinear); d(antho)/dt = −r; d(tannin)/dt = −y·r
    # (mass-based
    # consumption). Verify both exactly, and that nothing else moves (no o2, no carbon borrow).
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    antho, tannin = 0.3, 2.0
    y = _aged_wine(schema, ester=0.0, t=t, anthocyanin=antho, tannin=tannin)
    d = TanninAnthocyaninCondensation().derivatives(0.0, y, schema, poly_params)

    f_t = arrhenius_factor(t, poly_params["E_a_polymerization"], poly_params["T_ref"])
    k = poly_params["k_polymerization"]
    y_tannin = poly_params["y_tannin_per_anthocyanin"]
    r = k * f_t * antho * tannin
    assert schema.get(d, "anthocyanin") == pytest.approx(-r)
    assert schema.get(d, "tannin") == pytest.approx(-y_tannin * r)
    # The condensed anthocyanin is deposited into the polymeric_pigment slot (D-81), +r exactly
    # balancing the anthocyanin drawdown (anthocyanin-equivalents).
    assert schema.get(d, "polymeric_pigment") == pytest.approx(r)
    # A separate, non-oxidative grape axis: it touches NO o2 and borrows no carbon — nothing else.
    for var in ("X", "S", "E", "N", "CO2", "o2", "A420", "acetaldehyde", "ellagitannin"):
        assert schema.get(d, var) == 0.0


def test_polymerization_is_bilinear(poly_params):
    # Bilinear in BOTH grape drivers: doubling anthocyanin OR doubling tannin doubles the rate.
    schema = wine_schema()
    base = TanninAnthocyaninCondensation().derivatives(
        0.0, _aged_wine(schema, ester=0.0, anthocyanin=0.2, tannin=1.5), schema, poly_params
    )
    twice_antho = TanninAnthocyaninCondensation().derivatives(
        0.0, _aged_wine(schema, ester=0.0, anthocyanin=0.4, tannin=1.5), schema, poly_params
    )
    twice_tannin = TanninAnthocyaninCondensation().derivatives(
        0.0, _aged_wine(schema, ester=0.0, anthocyanin=0.2, tannin=3.0), schema, poly_params
    )
    assert schema.get(twice_antho, "anthocyanin") == pytest.approx(
        2.0 * schema.get(base, "anthocyanin")
    )
    assert schema.get(twice_tannin, "anthocyanin") == pytest.approx(
        2.0 * schema.get(base, "anthocyanin")
    )
    assert schema.get(base, "anthocyanin") < 0.0  # actually condensing


def test_polymerization_inert_without_anthocyanin_or_tannin(poly_params):
    # Doubly substrate-gated: no anthocyanin OR no tannin ⇒ byte-for-byte zero. A white wine (no
    # anthocyanin) or a no-tannin wine is exactly the case without this Process (isolability #3).
    schema = wine_schema()
    p = TanninAnthocyaninCondensation()
    no_antho = _aged_wine(schema, ester=0.0, anthocyanin=0.0, tannin=2.0)
    no_tannin = _aged_wine(schema, ester=0.0, anthocyanin=0.3, tannin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_antho, schema, poly_params), schema.zeros())
    assert np.array_equal(p.derivatives(0.0, no_tannin, schema, poly_params), schema.zeros())
    # Solver undershoot (negative) is likewise absorbed.
    undershoot = _aged_wine(schema, ester=0.0, anthocyanin=-1e-9, tannin=2.0)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_polymerization_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks k_polymerization,
    # yet a white wine (anthocyanin 0) returns zero — the gate-on-STATE-before-params discipline.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, anthocyanin=0.0, tannin=0.0)
    d = TanninAnthocyaninCondensation().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_polymerization_rises_with_temperature(poly_params):
    # Warmer condenses faster (E_a > 0, reaction-scale): the |anthocyanin| draw grows with T.
    schema = wine_schema()
    cold = TanninAnthocyaninCondensation().derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=283.15, anthocyanin=0.3, tannin=2.0),
        schema,
        poly_params,
    )
    warm = TanninAnthocyaninCondensation().derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=303.15, anthocyanin=0.3, tannin=2.0),
        schema,
        poly_params,
    )
    assert abs(float(schema.get(warm, "anthocyanin"))) > abs(float(schema.get(cold, "anthocyanin")))


def test_polymerization_wine_only_noop_on_beer(poly_params):
    # Wine-only (anthocyanin/tannin are appended to wine_schema): a hard no-op on beer (no slots).
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 100.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    assert np.array_equal(
        TanninAnthocyaninCondensation().derivatives(0.0, yb, beer, poly_params), beer.zeros()
    )


def test_polymerization_moves_nothing_conserved(poly_store, poly_params):
    # Both anthocyanin and tannin are off every ledger (grape-derived / carbon-free), so condensing
    # them moves NOTHING conserved — total_carbon, total_mass AND total_nitrogen all exactly flat
    # (the OakExtraction/EllagitanninOxidation off-every-ledger invariance). X=0 throughout.
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninAnthocyaninCondensation()], strict=True)
    y0 = _aged_wine(schema, ester=0.0, anthocyanin=0.3, tannin=2.0)
    traj = simulate(ps, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    f_c = poly_store.value("biomass_C_fraction")
    f_n = poly_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_polymerization_softens_and_stabilizes_colour(poly_params):
    # The D-79 spine, exercised over a long aging span (a red wine: both grape pools dosed):
    #   (1) astringency_series (free tannin) SOFTENS — tannin is drawn down monotonically as it
    #       condenses (here ellagitannin ≡ 0, so astringency == tannin × 1000);
    #   (2) polymeric_pigment_series RISES from 0 (the stable pigment = anthocyanin condensed);
    # (3) color_series is RETAINED (== initial anthocyanin × 1000, conserved in v1) — the observable
    #       dynamic is the monomeric → polymeric shift, not colour loss (bleaching deferred).
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninAnthocyaninCondensation()], strict=True)
    antho0, tannin0 = 0.3, 2.0
    y0 = _aged_wine(schema, ester=0.0, anthocyanin=antho0, tannin=tannin0)
    traj = simulate(ps, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0 * 2.0))
    assert traj.success

    # (1) astringency softens — monotone non-increasing, ends well below the start (tol >> the ~1e-7
    # mg/L adaptive-solver jitter at the tail plateau, << the ~900 mg/L signal). Tannin asymptotes
    # to tannin0 − y·antho0 as the limiting anthocyanin depletes (~1100 mg/L here — a MODEST draw,
    #     the "one directional contributor" scope: self-polymerization deferred).
    astr = astringency_series(traj)
    assert np.all(np.diff(astr) <= 1e-4)
    assert astr[-1] < astr[0] - 100.0
    # With no oak, astringency is exactly grape tannin × 1000.
    assert np.allclose(astr, np.asarray(traj.series("tannin"), dtype=float) * 1000.0)

    # (2) polymeric pigment rises from 0, read from the integrated polymeric_pigment slot (D-81).
    # With condensation as anthocyanin's ONLY fate here (no AnthocyaninFading in this ProcessSet),
    # the slot equals the anthocyanin condensed (antho0 − antho) × 1000 — a cross-check that the
    # promoted slot reproduces the old reconstruction exactly (the promotion preserves behaviour).
    pig = polymeric_pigment_series(traj)
    assert pig[0] == pytest.approx(0.0)
    assert pig[-1] > 0.0
    assert np.all(np.diff(pig) >= -1e-4)  # monotone non-decreasing (bar solver jitter)
    expected_pig = (antho0 - np.asarray(traj.series("anthocyanin"), dtype=float)) * 1000.0
    assert np.allclose(pig, expected_pig)

    # (3) total colour is RETAINED — free anthocyanin declines but polymeric pigment rises, so the
    # sum holds at the initial anthocyanin (condensation loses no colour, it stabilizes it). NOTE:
    # colour genuinely DECLINES only under AnthocyaninFading (D-81, the second anthocyanin fate → a
    # colourless sink) — see test_fading_makes_colour_decline. With ONLY condensation here,
    # color_series is flat at antho0 × 1000 (pigment gained == free anthocyanin lost), documenting
    # the stabilization physics; the genuine Process signal is the anthocyanin drawdown next.
    col = color_series(traj)
    assert np.allclose(col, antho0 * 1000.0)
    assert float(traj.series("anthocyanin")[-1]) < antho0  # free anthocyanin genuinely declined


def test_colour_form_identity_holds_by_construction(poly_params):
    # The D-81 three-slot colour identity: anthocyanin + polymeric_pigment + faded_anthocyanin ≡
    # anthocyanin₀ at ALL times. It holds BY CONSTRUCTION — condensation moves anthocyanin →
    # polymeric_pigment (−r, +r) and fading moves anthocyanin → faded_anthocyanin, so the three
    # d/dt terms sum to zero for any rate law. This CANNOT go through assert_conserved: all three
    # slots are off every ledger (weight 0), so it is a direct three-slot sum check (advisor note).
    # With only condensation in this ProcessSet, faded_anthocyanin stays ≡ 0 and the identity
    # reduces to anthocyanin + polymeric_pigment ≡ antho0 — the promotion's colour-conservation
    # proof.
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninAnthocyaninCondensation()], strict=True)
    antho0 = 0.3
    y0 = _aged_wine(schema, ester=0.0, anthocyanin=antho0, tannin=2.0)
    traj = simulate(ps, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0 * 2.0))
    assert traj.success
    antho = np.asarray(traj.series("anthocyanin"), dtype=float)
    pigment = np.asarray(traj.series("polymeric_pigment"), dtype=float)
    faded = np.asarray(traj.series("faded_anthocyanin"), dtype=float)
    assert np.allclose(antho + pigment + faded, antho0, atol=1e-9)
    assert np.allclose(faded, 0.0)  # no fade Process here ⇒ the colourless sink stays empty
    assert pigment[-1] > 0.0  # pigment genuinely accumulated in its own slot


def test_polymerization_tier_floored_at_speculative(poly_store):
    # Speculative in FORM (Tier-3 frontier) and capped by its speculative polymerization params
    # (D-1).
    # Both grape pools it writes are speculative pre- and post-parameter-tier propagation.
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninAnthocyaninCondensation()])
    for pool in ("anthocyanin", "tannin"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, poly_store.tier_map()) is Tier.SPECULATIVE


# -- D-80: AcetaldehydeBridgedCondensation — the acetaldehyde-bridged (ethylidene) / split-ledger --
#
# The ninth aging Process, the THIRD non-oxidative one and the D-79-deferred split-ledger colour
# beat:
# dissolved-O₂ acetaldehyde (D-71) bridges grape tannin to anthocyanin (tannin–ethyl–anthocyanin),
# the
# FIRST link from the oxidative sub-axis to red-wine colour. Unlike the D-79 direct route (moves
# nothing conserved), the bridged route consumes ON-ledger acetaldehyde (carbon borrowed from E at
# D-71), so a new on-ledger ``ethyl_bridge`` slot captures that carbon — the SPLIT LEDGER (grape
# bulk
# off-ledger, acetaldehyde-derived bridge on it). These tests pin: the trilinear closed form + the
# carbon-exact acetaldehyde→ethyl_bridge deposit, the trilinearity, the NON-TRIVIAL carbon closure
# (acetaldehyde↓ exactly cancels ethyl_bridge↑ — the split-ledger proof), the FREE-acetaldehyde read
# under SO₂ (bound acetaldehyde can't bridge — SO₂ delays colour stabilization, D-47), the exact
# unsulfited/undosed isolability, the warmer-faster ordering, the wine-only no-op on beer, and
# tiers.


@pytest.fixture
def bridge_store():
    # Wine params + acidbase.yaml (the pKa params free_acetaldehyde/ph_of_state read under SO₂) +
    # polymerization.yaml (k_acetaldehyde_bridge, E_a_acetaldehyde_bridge, y_acetaldehyde_per_
    # anthocyanin, y_tannin_per_anthocyanin, T_ref). acidbase is needed for the SO₂ path AND for any
    # simulate: the BDF numerical Jacobian perturbs so2_total off exact-zero, tripping the pH branch
    # (the SulfiteOxidation precedent — real compiled wine always loads acidbase.yaml).
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "acidbase.yaml",
        default_data_dir() / "polymerization.yaml",
    )


@pytest.fixture
def bridge_params(bridge_store):
    return bridge_store.resolve()


def test_bridge_metadata():
    p = AcetaldehydeBridgedCondensation()
    assert p.name == "acetaldehyde_bridged_condensation"
    assert p.tier is Tier.SPECULATIVE
    # Touches the two grape pools (off every ledger) PLUS the on-ledger acetaldehyde/ethyl_bridge
    # pair — the split ledger — PLUS the polymeric_pigment slot it fills (D-81, shared with the
    # direct route). The FIRST aging colour Process to touch the carbon ledger.
    assert set(p.touches) == {
        "acetaldehyde",
        "ethyl_bridge",
        "anthocyanin",
        "tannin",
        "polymeric_pigment",
    }
    assert set(p.reads) == {
        "k_acetaldehyde_bridge",
        "E_a_acetaldehyde_bridge",
        "y_acetaldehyde_per_anthocyanin",
        "y_tannin_per_anthocyanin",
        "T_ref",
    }


def test_bridge_closed_form(poly_params):
    # r = k·f(T)·[acetaldehyde]·[anthocyanin]·[tannin] (trilinear); anchored on anthocyanin. Verify
    # all four derivatives exactly, INCLUDING the carbon-exact acetaldehyde→ethyl_bridge split
    # (release at cf(acetaldehyde), redeposit at cf(ethylidene)). No SO₂ ⇒ reads total acetaldehyde.
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    acet, antho, tannin = 0.05, 0.3, 2.0
    y = _aged_wine(schema, ester=0.0, t=t, acetaldehyde=acet, anthocyanin=antho, tannin=tannin)
    d = AcetaldehydeBridgedCondensation().derivatives(0.0, y, schema, poly_params)

    f_t = arrhenius_factor(t, poly_params["E_a_acetaldehyde_bridge"], poly_params["T_ref"])
    k = poly_params["k_acetaldehyde_bridge"]
    y_tannin = poly_params["y_tannin_per_anthocyanin"]
    y_acet = poly_params["y_acetaldehyde_per_anthocyanin"]
    r = k * f_t * acet * antho * tannin
    acet_consumed = y_acet * r
    assert schema.get(d, "anthocyanin") == pytest.approx(-r)
    assert schema.get(d, "tannin") == pytest.approx(-y_tannin * r)
    assert schema.get(d, "acetaldehyde") == pytest.approx(-acet_consumed)
    # The split-ledger capture: the acetaldehyde carbon released is re-deposited into ethyl_bridge
    # at
    # the ethylidene fraction, so carbon released == carbon deposited (machine precision).
    expected_bridge = acet_consumed * _ACET_C / _ETHYLIDENE_C
    assert schema.get(d, "ethyl_bridge") == pytest.approx(expected_bridge)
    assert acet_consumed * _ACET_C == pytest.approx(schema.get(d, "ethyl_bridge") * _ETHYLIDENE_C)
    # Deposits the bridged anthocyanin into the SHARED polymeric_pigment slot (D-81), +r — same pool
    # the direct route fills (off-ledger colour-equivalent, distinct from the on-ledger ethyl_bridge
    # carbon; no double-count).
    assert schema.get(d, "polymeric_pigment") == pytest.approx(r)
    # Nothing else moves (no o2, no E borrow, no A420, no oak pools).
    for var in ("X", "S", "E", "N", "CO2", "o2", "A420", "ellagitannin"):
        assert schema.get(d, var) == 0.0


def test_bridge_is_trilinear(poly_params):
    # Trilinear: doubling acetaldehyde OR anthocyanin OR tannin each doubles the rate.
    schema = wine_schema()
    p = AcetaldehydeBridgedCondensation()

    def rate(acet: float, antho: float, tannin: float) -> float:
        y = _aged_wine(schema, ester=0.0, acetaldehyde=acet, anthocyanin=antho, tannin=tannin)
        return float(schema.get(p.derivatives(0.0, y, schema, poly_params), "anthocyanin"))

    base = rate(0.04, 0.2, 1.5)
    assert rate(0.08, 0.2, 1.5) == pytest.approx(2.0 * base)  # 2× acetaldehyde
    assert rate(0.04, 0.4, 1.5) == pytest.approx(2.0 * base)  # 2× anthocyanin
    assert rate(0.04, 0.2, 3.0) == pytest.approx(2.0 * base)  # 2× tannin
    assert base < 0.0  # actually bridging


def test_bridge_inert_without_any_substrate(poly_params):
    # Triply substrate-gated: no acetaldehyde OR no anthocyanin OR no tannin ⇒ byte-for-byte zero,
    # so
    # a white / no-tannin / no-acetaldehyde (reductive, un-oxygenated) wine is exactly the case
    # without this Process (isolability #3). Solver undershoot (negative) is likewise absorbed.
    schema = wine_schema()
    p = AcetaldehydeBridgedCondensation()
    for kw in (
        {"acetaldehyde": 0.0, "anthocyanin": 0.3, "tannin": 2.0},
        {"acetaldehyde": 0.05, "anthocyanin": 0.0, "tannin": 2.0},
        {"acetaldehyde": 0.05, "anthocyanin": 0.3, "tannin": 0.0},
        {"acetaldehyde": -1e-9, "anthocyanin": 0.3, "tannin": 2.0},
    ):
        y = _aged_wine(schema, ester=0.0, **kw)
        assert np.array_equal(p.derivatives(0.0, y, schema, poly_params), schema.zeros())


def test_bridge_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks
    # k_acetaldehyde_bridge, yet a white wine (anthocyanin 0) returns zero — the
    # gate-on-STATE-before-params discipline.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, acetaldehyde=0.05, anthocyanin=0.0, tannin=0.0)
    d = AcetaldehydeBridgedCondensation().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_bridge_reads_free_acetaldehyde_under_so2(bridge_params):
    # SO₂-bound acetaldehyde is the bisulfite adduct — its carbonyl is blocked, so it CANNOT bridge
    # (the D-47 precedent). At the SAME total acetaldehyde, a sulfited wine bridges SLOWER than an
    # unsulfited one (some acetaldehyde is bound and unavailable) — SO₂ DELAYS colour stabilization.
    schema = wine_schema()
    p = AcetaldehydeBridgedCondensation()
    acids = {"tartaric": 4.0, "cation_charge": 0.012}  # a real acid state so pH solves
    unsulfited = _aged_wine(
        schema, ester=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.0, **acids
    )
    sulfited = _aged_wine(
        schema, ester=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.05, **acids
    )
    r_unsulfited = float(
        schema.get(p.derivatives(0.0, unsulfited, schema, bridge_params), "anthocyanin")
    )
    r_sulfited = float(
        schema.get(p.derivatives(0.0, sulfited, schema, bridge_params), "anthocyanin")
    )
    # Both bridge (negative), but the sulfited rate is throttled toward zero (less free
    # acetaldehyde).
    assert r_sulfited < 0.0
    assert r_unsulfited < 0.0
    assert abs(r_sulfited) < abs(r_unsulfited)


def test_bridge_unsulfited_is_exactly_total_acetaldehyde(bridge_params):
    # The so2_total > 0 guard is EXACT: at so2_total = 0 the rate uses TOTAL acetaldehyde (no pH
    # solve), so it equals the trilinear closed form on total acetaldehyde — byte-for-byte the
    # no-SO₂-branch case. (This is why an unsulfited aging run pays no per-RHS brentq.)
    schema = wine_schema()
    acet, antho, tannin, t = 0.05, 0.3, 2.0, 298.15
    y = _aged_wine(
        schema, ester=0.0, t=t, acetaldehyde=acet, anthocyanin=antho, tannin=tannin, so2_total=0.0
    )
    d = AcetaldehydeBridgedCondensation().derivatives(0.0, y, schema, bridge_params)
    f_t = arrhenius_factor(t, bridge_params["E_a_acetaldehyde_bridge"], bridge_params["T_ref"])
    r = bridge_params["k_acetaldehyde_bridge"] * f_t * acet * antho * tannin
    assert schema.get(d, "anthocyanin") == pytest.approx(-r)


def test_bridge_rises_with_temperature(poly_params):
    # Warmer bridges faster (E_a > 0, reaction-scale): the |anthocyanin| draw grows with T.
    schema = wine_schema()
    p = AcetaldehydeBridgedCondensation()
    cold = p.derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=283.15, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0),
        schema,
        poly_params,
    )
    warm = p.derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=303.15, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0),
        schema,
        poly_params,
    )
    assert abs(float(schema.get(warm, "anthocyanin"))) > abs(float(schema.get(cold, "anthocyanin")))


def test_bridge_wine_only_noop_on_beer(poly_params):
    # Wine-only (anthocyanin/tannin/ethyl_bridge are appended to wine_schema): a hard no-op on beer.
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 100.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    yb[beer.slice("acetaldehyde")] = 0.05  # acetaldehyde is medium-agnostic, but no grape pools
    assert np.array_equal(
        AcetaldehydeBridgedCondensation().derivatives(0.0, yb, beer, poly_params), beer.zeros()
    )


def test_bridge_carbon_closes_nontrivially(bridge_store, bridge_params):
    # THE SPLIT-LEDGER PROOF. Unlike the D-79 direct route (moves nothing conserved), the bridged
    # route genuinely MOVES carbon: it consumes on-ledger acetaldehyde and books its carbon into the
    # on-ledger ethyl_bridge slot. total_carbon is flat NON-TRIVIALLY — acetaldehyde↓ exactly
    # cancels
    # ethyl_bridge↑ (that cancellation IS the split-ledger accounting). total_mass/total_nitrogen
    # are
    # flat too, but trivially (the Process touches no {S,E,CO2} or N species). X=0 throughout.
    schema = wine_schema()
    ps = ProcessSet(schema, [AcetaldehydeBridgedCondensation()], strict=True)
    # A real acid state so the BDF Jacobian's so2_total perturbation solves pH cleanly (so2_total=0,
    # so the exact RHS reads total acetaldehyde — the acids are inert here, just pH-well-posed).
    y0 = _aged_wine(
        schema,
        ester=0.0,
        acetaldehyde=0.05,
        anthocyanin=0.3,
        tannin=2.0,
        tartaric=4.0,
        cation_charge=0.012,
    )
    traj = simulate(ps, params=bridge_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    f_c = bridge_store.value("biomass_C_fraction")
    f_n = bridge_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")
    # The closure is NON-trivial: acetaldehyde genuinely fell and ethyl_bridge genuinely rose, and
    # the
    # carbon lost from one equals the carbon gained by the other (machine precision).
    acet = np.asarray(traj.series("acetaldehyde"), dtype=float)
    bridge = np.asarray(traj.series("ethyl_bridge"), dtype=float)
    assert acet[-1] < acet[0] - 1e-3  # acetaldehyde consumed
    assert bridge[-1] > bridge[0] + 1e-4  # ethyl_bridge accumulated
    dC_acet = (acet[-1] - acet[0]) * _ACET_C
    dC_bridge = (bridge[-1] - bridge[0]) * _ETHYLIDENE_C
    assert dC_acet + dC_bridge == pytest.approx(0.0, abs=1e-12)


def test_bridge_tier_floored_at_speculative(bridge_store):
    # Speculative in FORM (Tier-3 frontier) and capped by its speculative polymerization params
    # (D-1),
    # for every pool it writes — including the on-ledger acetaldehyde/ethyl_bridge pair.
    schema = wine_schema()
    ps = ProcessSet(schema, [AcetaldehydeBridgedCondensation()])
    for pool in ("anthocyanin", "tannin", "acetaldehyde", "ethyl_bridge"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, bridge_store.tier_map()) is Tier.SPECULATIVE


# -- D-81: AnthocyaninFading — the O₂-coupled oxidative bleaching loss --------------------------
#
# The tenth aging Process and the beat that makes color_series genuinely DECLINE: dissolved O₂ fades
# free monomeric anthocyanin to the COLOURLESS faded_anthocyanin slot (a second anthocyanin fate
# besides condensation into stable pigment). Bilinear O₂ sink (the EllagitanninOxidation form), a
# pure off-ledger transfer, drawing the SHARED o2 pool — so SO₂ protection is EMERGENT (SO₂
# o2 via D-72, leaving less to fade). These tests pin the closed form, the bilinearity, the gates,
# the wine-only no-op, off-every-ledger invariance, the genuine colour DECLINE + the identity,
# the emergent SO₂ protection, and tiers.


@pytest.fixture
def fade_store():
    # Wine + polymerization.yaml (k_anthocyanin_fade, E_a_anthocyanin_fade, y_anthocyanin_per_o2 +
    # the condensation params) + aging.yaml (k_so2_oxidation, for the emergent-SO₂-protection test's
    # SulfiteOxidation) + the acidbase/acetaldehyde/keto-acid pKa params SulfiteOxidation's
    # pH/bisulfite readout reads — the full oxidative-axis + fade parameter set.
    d = default_data_dir()
    return load_parameters(
        d / "wine_generic.yaml",
        d / "acidbase.yaml",
        d / "acetaldehyde.yaml",
        d / "keto_acids.yaml",
        d / "aging.yaml",
        d / "polymerization.yaml",
    )


@pytest.fixture
def fade_params(fade_store):
    return fade_store.resolve()


def test_fading_metadata():
    p = AnthocyaninFading()
    assert p.name == "anthocyanin_fading"
    assert p.tier is Tier.SPECULATIVE
    # Draws the SHARED o2 pool (so SO₂ protection is emergent) and TRANSFERS anthocyanin into the
    # colourless faded_anthocyanin slot — all three off every ledger, so nothing conserved moves.
    assert set(p.touches) == {"o2", "anthocyanin", "faded_anthocyanin"}
    assert "o2" in p.touches  # O₂-COUPLED (the D-81 crux — SO₂ protection is emergent)
    assert set(p.reads) == {
        "k_anthocyanin_fade",
        "E_a_anthocyanin_fade",
        "y_anthocyanin_per_o2",
        "T_ref",
    }


def test_fading_closed_form(poly_params):
    # r_o2 = k·f(T)·[o2]·[anthocyanin] (bilinear); d(o2)/dt = −r_o2; anthocyanin → faded at a
    # mass-based yield (a pure transfer). Verify all three derivatives exactly and that nothing else
    # moves — in particular NO pigment (fading is colourless), NO tannin, NO carbon borrow.
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    o2, antho = 0.03, 0.3
    y = _aged_wine(schema, ester=0.0, t=t, o2=o2, anthocyanin=antho)
    d = AnthocyaninFading().derivatives(0.0, y, schema, poly_params)

    f_t = arrhenius_factor(t, poly_params["E_a_anthocyanin_fade"], poly_params["T_ref"])
    k = poly_params["k_anthocyanin_fade"]
    y_fade = poly_params["y_anthocyanin_per_o2"]
    r_o2 = k * f_t * o2 * antho
    faded = y_fade * r_o2
    assert schema.get(d, "o2") == pytest.approx(-r_o2)
    assert schema.get(d, "anthocyanin") == pytest.approx(-faded)
    assert schema.get(d, "faded_anthocyanin") == pytest.approx(faded)
    # Pure transfer: anthocyanin lost == faded gained (the colour identity closes by construction).
    assert schema.get(d, "anthocyanin") == pytest.approx(-schema.get(d, "faded_anthocyanin"))
    # Colourless fade adds NO pigment, and nothing else moves (no tannin, no E/CO2, no A420).
    for var in ("X", "S", "E", "N", "CO2", "A420", "tannin", "polymeric_pigment", "acetaldehyde"):
        assert schema.get(d, var) == 0.0


def test_fading_is_bilinear(poly_params):
    # Bilinear in BOTH drivers: doubling o2 OR doubling anthocyanin doubles the fade rate.
    schema = wine_schema()
    p = AnthocyaninFading()

    def fade_rate(o2: float, antho: float) -> float:
        y = _aged_wine(schema, ester=0.0, o2=o2, anthocyanin=antho)
        return float(schema.get(p.derivatives(0.0, y, schema, poly_params), "anthocyanin"))

    base = fade_rate(0.02, 0.2)
    assert fade_rate(0.04, 0.2) == pytest.approx(2.0 * base)  # 2× o2
    assert fade_rate(0.02, 0.4) == pytest.approx(2.0 * base)  # 2× anthocyanin
    assert base < 0.0  # actually fading


def test_fading_inert_without_o2_or_anthocyanin(poly_params):
    # Doubly substrate-gated: no O₂ (reductive) OR no anthocyanin (white) ⇒ byte-for-byte zero, so
    # such a run is exactly the case without this Process (isolability #3). Undershoot absorbed too.
    schema = wine_schema()
    p = AnthocyaninFading()
    no_o2 = _aged_wine(schema, ester=0.0, o2=0.0, anthocyanin=0.3)
    no_antho = _aged_wine(schema, ester=0.0, o2=0.03, anthocyanin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_o2, schema, poly_params), schema.zeros())
    assert np.array_equal(p.derivatives(0.0, no_antho, schema, poly_params), schema.zeros())
    undershoot = _aged_wine(schema, ester=0.0, o2=-1e-9, anthocyanin=0.3)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_fading_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks the fade rate,
    # yet a white wine (anthocyanin 0, even with O₂ dosed) returns zero — gate-before-params.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, o2=0.03, anthocyanin=0.0)
    d = AnthocyaninFading().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_fading_rises_with_temperature(poly_params):
    # Warmer fades faster (E_a > 0, reaction-scale): the |anthocyanin| draw grows with T.
    schema = wine_schema()
    p = AnthocyaninFading()
    cold = p.derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=283.15, o2=0.03, anthocyanin=0.3), schema, poly_params
    )
    warm = p.derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=303.15, o2=0.03, anthocyanin=0.3), schema, poly_params
    )
    assert abs(float(schema.get(warm, "anthocyanin"))) > abs(float(schema.get(cold, "anthocyanin")))


def test_fading_wine_only_noop_on_beer(poly_params):
    # Wine-only (anthocyanin/faded_anthocyanin are appended to wine_schema): a hard no-op on beer
    # (o2 exists in both media, but the grape colour slots do not).
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 100.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    yb[beer.slice("o2")] = 0.03
    assert np.array_equal(AnthocyaninFading().derivatives(0.0, yb, beer, poly_params), beer.zeros())


def test_fading_moves_nothing_conserved(poly_store, poly_params):
    # o2 (carbon-free) and both anthocyanin/faded_anthocyanin (grape-derived) are off every ledger,
    # so fading anthocyanin to colourless products moves NOTHING conserved — carbon, mass
    # AND total_nitrogen all exactly flat (the EllagitanninOxidation off-every-ledger invariance).
    schema = wine_schema()
    ps = ProcessSet(schema, [AnthocyaninFading()], strict=True)
    y0 = _aged_wine(schema, ester=0.0, o2=0.04, anthocyanin=0.3)
    traj = simulate(ps, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    f_c = poly_store.value("biomass_C_fraction")
    f_n = poly_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_fading_makes_colour_decline(poly_params):
    # THE D-81 HEADLINE: with condensation AND fading acting together on an oxygenated red, colour
    # GENUINELY declines (unlike the D-79/D-80 flat line). Free anthocyanin leaves partly to STABLE
    # pigment (colour retained) and partly to COLOURLESS faded (colour lost), so color_series falls
    # by exactly the faded amount while the polymeric pigment survives — the stability payoff.
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninAnthocyaninCondensation(), AnthocyaninFading()], strict=True)
    antho0 = 0.3
    y0 = _aged_wine(schema, ester=0.0, o2=0.05, anthocyanin=antho0, tannin=2.0)
    traj = simulate(ps, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0 * 2.0))
    assert traj.success

    antho = np.asarray(traj.series("anthocyanin"), dtype=float)
    pigment = np.asarray(traj.series("polymeric_pigment"), dtype=float)
    faded = np.asarray(traj.series("faded_anthocyanin"), dtype=float)
    col = color_series(traj)

    # (1) colour GENUINELY declines (the whole point of D-81) — the pigment slot makes this real.
    assert col[-1] < col[0] - 10.0  # well beyond solver jitter; col[0] == antho0 × 1000 == 300
    assert np.all(np.diff(col) <= 1e-6)  # monotone non-increasing
    # (2) it falls by exactly the faded amount: color ≡ (antho0 − faded) × 1000.
    assert np.allclose(col, (antho0 - faded) * 1000.0)
    # (3) the stable pigment SURVIVED (condensation still ran) — the colour-stability payoff.
    assert pigment[-1] > 0.0
    assert faded[-1] > 0.0  # some anthocyanin genuinely bleached to colourless
    # (4) the three-slot colour identity holds by construction, now with faded > 0 (non-trivially).
    assert np.allclose(antho + pigment + faded, antho0, atol=1e-9)


# -- D-83: ThermalAnthocyaninFade — the O₂-INDEPENDENT thermal/hydrolytic bleaching loss --------
# ThermalAnthocyaninFade is the second, O₂-independent fate that fades free monomeric anthocyanin to
# the SAME colourless faded_anthocyanin slot the D-81 oxidative fade fills — but by a thermal route
# needing NO oxygen (first-order [anthocyanin], the EsterHydrolysis form, NOT the D-81 bilinear
# o2·anthocyanin). Its crux is the MIRROR of D-81: touching no o2, SO₂ does NOT protect it, so a
# sealed/sulfited/anaerobic red still fades and only cold storage slows it. These tests pin the
# closed form, the first-order (o2-free) shape, the single gate, the wine-only no-op, off-every-
# ledger invariance, that it fades a REDUCTIVE red (retiring the D-81 anaerobic-holds note), that
# SO₂ does NOT protect it (vs the D-81 contrast), the shared faded sink + identity, and tiers.


def test_thermal_fade_metadata():
    p = ThermalAnthocyaninFade()
    assert p.name == "thermal_anthocyanin_fade"
    assert p.tier is Tier.SPECULATIVE
    # TRANSFERS anthocyanin into the colourless faded_anthocyanin slot — both off every ledger, so
    # nothing conserved moves. It touches NO o2 (the D-83 crux — O₂-INDEPENDENT, so SO₂ can't
    # protect it), unlike its D-81 oxidative sibling.
    assert set(p.touches) == {"anthocyanin", "faded_anthocyanin"}
    assert "o2" not in p.touches  # O₂-INDEPENDENT (the D-83 crux, the mirror of D-81)
    # No yield (contrast D-81): the rate is already g anthocyanin/L/h — a direct −r/+r transfer.
    assert set(p.reads) == {
        "k_anthocyanin_thermal_fade",
        "E_a_anthocyanin_thermal_fade",
        "T_ref",
    }


def test_thermal_fade_closed_form(poly_params):
    # r = k·f(T)·[anthocyanin] (FIRST-ORDER, no o2 term); anthocyanin → faded 1:1 (a pure transfer,
    # no yield). Verify both derivatives exactly and that nothing else moves — in particular NO o2
    # draw (this is the O₂-free route), NO pigment (fading is colourless), NO tannin, NO carbon.
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    antho = 0.3
    y = _aged_wine(schema, ester=0.0, t=t, o2=0.03, anthocyanin=antho)
    d = ThermalAnthocyaninFade().derivatives(0.0, y, schema, poly_params)

    f_t = arrhenius_factor(t, poly_params["E_a_anthocyanin_thermal_fade"], poly_params["T_ref"])
    r = poly_params["k_anthocyanin_thermal_fade"] * f_t * antho
    assert schema.get(d, "anthocyanin") == pytest.approx(-r)
    assert schema.get(d, "faded_anthocyanin") == pytest.approx(r)
    # Pure transfer: anthocyanin lost == faded gained (the colour identity closes by construction).
    assert schema.get(d, "anthocyanin") == pytest.approx(-schema.get(d, "faded_anthocyanin"))
    # O₂-FREE (the whole point): despite o2 present in state, this route draws NONE. And nothing
    # else moves — no pigment, no tannin, no E/CO2, no A420.
    for var in ("o2", "X", "S", "E", "N", "CO2", "A420", "tannin", "polymeric_pigment"):
        assert schema.get(d, var) == 0.0


def test_thermal_fade_is_first_order_and_o2_independent(poly_params):
    # FIRST-ORDER in anthocyanin (doubling it doubles the rate) and INDEPENDENT of o2 (doubling —
    # or zeroing — o2 leaves the rate unchanged, the D-83 crux vs the D-81 bilinear o2 sink).
    schema = wine_schema()
    p = ThermalAnthocyaninFade()

    def fade_rate(o2: float, antho: float) -> float:
        y = _aged_wine(schema, ester=0.0, o2=o2, anthocyanin=antho)
        return float(schema.get(p.derivatives(0.0, y, schema, poly_params), "anthocyanin"))

    base = fade_rate(0.02, 0.2)
    assert fade_rate(0.02, 0.4) == pytest.approx(2.0 * base)  # 2× anthocyanin ⇒ 2× rate
    assert fade_rate(0.04, 0.2) == pytest.approx(base)  # 2× o2 ⇒ SAME rate (o2-independent)
    assert fade_rate(0.0, 0.2) == pytest.approx(base)  # ZERO o2 ⇒ SAME rate (fades anaerobically!)
    assert base < 0.0  # actually fading


def test_thermal_fade_inert_without_anthocyanin(poly_params):
    # Singly substrate-gated on anthocyanin only (NO o2 gate — it fades even with zero o2). A white
    # wine (no anthocyanin) ⇒ byte-for-byte zero, so it is exactly the case without this Process.
    schema = wine_schema()
    p = ThermalAnthocyaninFade()
    no_antho = _aged_wine(schema, ester=0.0, o2=0.03, anthocyanin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_antho, schema, poly_params), schema.zeros())
    # Undershoot absorbed too.
    undershoot = _aged_wine(schema, ester=0.0, o2=0.03, anthocyanin=-1e-9)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_thermal_fade_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks the thermal-fade
    # rate, yet a white wine (anthocyanin 0) returns zero — gate-before-params.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, o2=0.03, anthocyanin=0.0)
    d = ThermalAnthocyaninFade().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_thermal_fade_rises_with_temperature(poly_params):
    # Warmer fades faster (E_a > 0, reaction-scale): the |anthocyanin| draw grows with T — the
    # 'warm storage kills colour even anaerobically' temperature lever, the only lever this route
    # has (no o2/SO₂ coupling).
    schema = wine_schema()
    p = ThermalAnthocyaninFade()
    cold = p.derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=283.15, o2=0.0, anthocyanin=0.3), schema, poly_params
    )
    warm = p.derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=303.15, o2=0.0, anthocyanin=0.3), schema, poly_params
    )
    assert abs(float(schema.get(warm, "anthocyanin"))) > abs(float(schema.get(cold, "anthocyanin")))


def test_thermal_fade_wine_only_noop_on_beer(poly_params):
    # Wine-only (anthocyanin/faded_anthocyanin are appended to wine_schema): a hard no-op on beer.
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 100.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    d = ThermalAnthocyaninFade().derivatives(0.0, yb, beer, poly_params)
    assert np.array_equal(d, beer.zeros())


def test_thermal_fade_moves_nothing_conserved(poly_store, poly_params):
    # Both anthocyanin/faded_anthocyanin (grape-derived) are off every ledger, so thermal fading
    # anthocyanin to colourless products moves NOTHING conserved — carbon, mass AND nitrogen all
    # exactly flat (the AnthocyaninFading off-every-ledger invariance, now with no o2 involved).
    schema = wine_schema()
    ps = ProcessSet(schema, [ThermalAnthocyaninFade()], strict=True)
    y0 = _aged_wine(schema, ester=0.0, o2=0.0, anthocyanin=0.3)  # o2=0: fades anyway
    traj = simulate(ps, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    f_c = poly_store.value("biomass_C_fraction")
    f_n = poly_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_thermal_fade_bleaches_a_reductive_red(poly_params):
    # THE D-83 HEADLINE, retiring the D-81 "anaerobic sealed red holds its colour" note. With ONLY
    # ThermalAnthocyaninFade on a fully REDUCTIVE (o2 = 0) red — which is byte-for-byte FLAT under
    # AnthocyaninFading alone (D-81 fades only via o2) — colour now GENUINELY declines, purely
    # thermally, into the colourless faded slot. Contrast pinned: the D-81 oxidative fade does
    # nothing here (no o2).
    schema = wine_schema()
    antho0 = 0.3
    y0 = _aged_wine(schema, ester=0.0, o2=0.0, anthocyanin=antho0)  # anaerobic, sealed red

    # D-81 alone on a reductive red: byte-for-byte flat (fades only via o2).
    ps_oxid = ProcessSet(schema, [AnthocyaninFading()], strict=True)
    traj_oxid = simulate(ps_oxid, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0 * 2.0))
    assert traj_oxid.success
    assert color_series(traj_oxid)[-1] == pytest.approx(antho0 * 1000.0)  # unchanged, D-81 inert

    # D-83 alone on the SAME reductive red: colour genuinely declines (the retirement).
    ps_therm = ProcessSet(schema, [ThermalAnthocyaninFade()], strict=True)
    traj = simulate(ps_therm, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0 * 2.0))
    assert traj.success
    antho = np.asarray(traj.series("anthocyanin"), dtype=float)
    faded = np.asarray(traj.series("faded_anthocyanin"), dtype=float)
    col = color_series(traj)
    assert col[-1] < col[0] - 10.0  # genuinely declines (col[0] == antho0 × 1000 == 300)
    assert np.all(np.diff(col) <= 1e-6)  # monotone non-increasing
    assert col[-1] == pytest.approx((antho0 - faded[-1]) * 1000.0)  # falls by exactly the faded amt
    # Two-slot identity (no pigment/o2 route here): anthocyanin + faded ≡ anthocyanin₀.
    assert np.allclose(antho + faded, antho0, atol=1e-9)


def test_thermal_fade_unprotected_by_so2(fade_params):
    # THE D-83 MIRROR OF D-81: SO₂ does NOT protect the thermal route. Because it draws no o2, a
    # heavily-sulfited red fades thermally EXACTLY as an unsulfited one — the physically-honest
    # split vs D-81, where SO₂ protects the o2-coupled fade emergently. Two identical reds, one
    # dosed with SO₂: the thermal anthocyanin draw is byte-for-byte the same.
    schema = wine_schema()
    p = ThermalAnthocyaninFade()
    no_so2 = _aged_wine(schema, ester=0.0, o2=0.0, anthocyanin=0.3)
    with_so2 = _aged_wine(schema, ester=0.0, o2=0.0, anthocyanin=0.3, so2_total=0.05)
    d_no = p.derivatives(0.0, no_so2, schema, fade_params)
    d_so2 = p.derivatives(0.0, with_so2, schema, fade_params)
    # SO₂ present or not, the thermal fade rate is identical (no o2 to scavenge ⇒ no protection).
    assert schema.get(d_so2, "anthocyanin") == pytest.approx(schema.get(d_no, "anthocyanin"))
    assert schema.get(d_so2, "anthocyanin") < 0.0  # and it IS fading despite the SO₂ dose


def test_thermal_fade_tier_is_speculative(poly_store):
    # Parameter-tier propagation (D-1): both fade params are speculative, so the anthocyanin /
    # faded_anthocyanin outputs report speculative (the aging-axis frontier).
    schema = wine_schema()
    ps = ProcessSet(schema, [ThermalAnthocyaninFade()], strict=True)
    tier_map = poly_store.tier_map()
    assert ps.tier_of("anthocyanin", tier_map) is Tier.SPECULATIVE
    assert ps.tier_of("faded_anthocyanin", tier_map) is Tier.SPECULATIVE


# -- D-84: TanninSelfPolymerization — the direct, non-oxidative tannin–tannin softener ----------
# TanninSelfPolymerization condenses grape tannin WITH ITSELF (bimolecular [tannin]², a true self-
# reaction) into a soft polymer, drawing the free-tannin pool down as a PURE off-ledger sink (the
# soft polymer goes to no slot — the D-79/D-80 tannin-is-a-pure-sink precedent). So astringency
# softens WITHOUT needing anthocyanin (retiring the D-80 "one-directional-per-pool" note). These
# tests pin the closed form, the [tannin]² second-order shape, the single gate, the wine-only no-op,
# off-every-ledger invariance, the anthocyanin-free softening headline, o2/anthocyanin-independence,
# and tiers.


def test_tannin_self_poly_metadata():
    p = TanninSelfPolymerization()
    assert p.name == "tannin_self_polymerization"
    assert p.tier is Tier.SPECULATIVE
    # Consumes the single grape tannin pool as a PURE off-ledger sink (the soft polymer goes to no
    # slot), so nothing conserved moves. No o2 (not an oxidation), no acetaldehyde (the DIRECT route
    # — bridged tannin–ethyl–tannin is D-85), no anthocyanin (a self-reaction).
    assert set(p.touches) == {"tannin"}
    assert "o2" not in p.touches and "anthocyanin" not in p.touches
    # No yield (a self-reaction: one pool, no second reactant).
    assert set(p.reads) == {
        "k_tannin_self_polymerization",
        "E_a_tannin_self_polymerization",
        "T_ref",
    }


def test_tannin_self_poly_closed_form(poly_params):
    # r = k·f(T)·[tannin]² (BIMOLECULAR self-reaction); d(tannin)/dt = −r; NOTHING else moves (pure
    # off-ledger sink — no destination slot, no o2, no anthocyanin, no acetaldehyde/carbon).
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    tannin = 2.0
    y = _aged_wine(schema, ester=0.0, t=t, anthocyanin=0.3, tannin=tannin)
    d = TanninSelfPolymerization().derivatives(0.0, y, schema, poly_params)

    f_t = arrhenius_factor(t, poly_params["E_a_tannin_self_polymerization"], poly_params["T_ref"])
    r = poly_params["k_tannin_self_polymerization"] * f_t * tannin * tannin
    assert schema.get(d, "tannin") == pytest.approx(-r)
    assert r > 0.0
    # Pure sink: nothing else moves — no o2 draw, no anthocyanin, no pigment, no carbon, no faded.
    for var in (
        "o2",
        "anthocyanin",
        "polymeric_pigment",
        "faded_anthocyanin",
        "acetaldehyde",
        "ethyl_bridge",
        "E",
        "CO2",
        "A420",
    ):
        assert schema.get(d, var) == 0.0


def test_tannin_self_poly_is_second_order(poly_params):
    # BIMOLECULAR in tannin: doubling tannin QUADRUPLES the rate (2² — the self-reaction signature,
    # distinct from the D-79 bilinear two-pool form which would only double).
    schema = wine_schema()
    p = TanninSelfPolymerization()

    def rate(tannin: float) -> float:
        y = _aged_wine(schema, ester=0.0, anthocyanin=0.3, tannin=tannin)
        return float(schema.get(p.derivatives(0.0, y, schema, poly_params), "tannin"))

    base = rate(1.0)
    assert rate(2.0) == pytest.approx(4.0 * base)  # 2× tannin ⇒ 4× rate (second-order)
    assert base < 0.0  # actually consuming tannin


def test_tannin_self_poly_independent_of_anthocyanin_and_o2(poly_params):
    # It needs NEITHER anthocyanin NOR o2 (the D-84 crux): the tannin draw is identical whether the
    # wine has anthocyanin/o2 or not — a self-reaction on the tannin pool alone.
    schema = wine_schema()
    p = TanninSelfPolymerization()

    def rate(anthocyanin: float, o2: float) -> float:
        y = _aged_wine(schema, ester=0.0, anthocyanin=anthocyanin, o2=o2, tannin=2.0)
        return float(schema.get(p.derivatives(0.0, y, schema, poly_params), "tannin"))

    base = rate(0.3, 0.03)
    assert rate(0.0, 0.0) == pytest.approx(base)  # no anthocyanin, no o2 ⇒ SAME rate
    assert rate(0.6, 0.06) == pytest.approx(base)  # more anthocyanin/o2 ⇒ SAME rate
    assert base < 0.0


def test_tannin_self_poly_inert_without_tannin(poly_params):
    # Singly substrate-gated on tannin: a no-tannin wine ⇒ byte-for-byte zero, so it is exactly the
    # case without this Process. Undershoot absorbed too.
    schema = wine_schema()
    p = TanninSelfPolymerization()
    no_tannin = _aged_wine(schema, ester=0.0, anthocyanin=0.3, tannin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_tannin, schema, poly_params), schema.zeros())
    undershoot = _aged_wine(schema, ester=0.0, anthocyanin=0.3, tannin=-1e-9)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_tannin_self_poly_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks the rate, yet a
    # no-tannin wine returns zero — gate-before-params.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, anthocyanin=0.3, tannin=0.0)
    d = TanninSelfPolymerization().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_tannin_self_poly_rises_with_temperature(poly_params):
    # Warmer polymerizes (softens) faster (E_a > 0, reaction-scale): the |tannin| draw grows with T.
    schema = wine_schema()
    p = TanninSelfPolymerization()
    cold = p.derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=283.15, tannin=2.0), schema, poly_params
    )
    warm = p.derivatives(
        0.0, _aged_wine(schema, ester=0.0, t=303.15, tannin=2.0), schema, poly_params
    )
    assert abs(float(schema.get(warm, "tannin"))) > abs(float(schema.get(cold, "tannin")))


def test_tannin_self_poly_wine_only_noop_on_beer(poly_params):
    # Wine-only (tannin is appended to wine_schema): a hard no-op on beer (no tannin slot).
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 100.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    d = TanninSelfPolymerization().derivatives(0.0, yb, beer, poly_params)
    assert np.array_equal(d, beer.zeros())


def test_tannin_self_poly_moves_nothing_conserved(poly_store, poly_params):
    # tannin (grape-derived) is off every ledger, so self-polymerizing it into a soft polymer moves
    # NOTHING conserved — carbon, mass AND nitrogen all exactly flat (the pure off-ledger sink).
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninSelfPolymerization()], strict=True)
    y0 = _aged_wine(schema, ester=0.0, tannin=3.0)  # no anthocyanin — softens anyway
    traj = simulate(ps, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    f_c = poly_store.value("biomass_C_fraction")
    f_n = poly_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")


def test_tannin_self_poly_softens_without_anthocyanin(poly_params):
    # THE D-84 HEADLINE, retiring the D-80 "softening needs anthocyanin" honesty note. On a wine
    # with tannin but ZERO anthocyanin (white, or exhausted red) — where BOTH condensation
    # routes are inert (they need anthocyanin) — astringency now GENUINELY softens, purely by tannin
    # self-polymerization. Contrast pinned: the D-79 direct condensation does nothing here.
    schema = wine_schema()
    tannin0 = 3.0
    y0 = _aged_wine(schema, ester=0.0, anthocyanin=0.0, tannin=tannin0)  # white / no anthocyanin

    # D-79 alone with no anthocyanin: byte-for-byte flat (condensation needs anthocyanin).
    ps_cond = ProcessSet(schema, [TanninAnthocyaninCondensation()], strict=True)
    traj_cond = simulate(ps_cond, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0 * 2.0))
    assert traj_cond.success
    assert astringency_series(traj_cond)[-1] == pytest.approx(tannin0 * 1000.0)  # unchanged

    # D-84 alone on the SAME no-anthocyanin wine: astringency genuinely softens (the retirement).
    ps_self = ProcessSet(schema, [TanninSelfPolymerization()], strict=True)
    traj = simulate(ps_self, params=poly_params, y0=y0, t_span=(0.0, 24.0 * 365.0 * 2.0))
    assert traj.success
    astr = astringency_series(traj)
    assert astr[-1] < astr[0] - 10.0  # genuinely softens (astr[0] == tannin0 × 1000 == 3000)
    assert np.all(np.diff(astr) <= 1e-6)  # monotone non-increasing (tannin only declines)


def test_tannin_self_poly_tier_is_speculative(poly_store):
    # Parameter-tier propagation (D-1): both params are speculative, so the tannin output reports
    # speculative (the aging-axis frontier).
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninSelfPolymerization()], strict=True)
    assert ps.tier_of("tannin", poly_store.tier_map()) is Tier.SPECULATIVE


# -- D-85: TanninEthylTanninCondensation — the acetaldehyde-bridged tannin–ethyl–tannin softener ---
# TanninEthylTanninCondensation bridges TWO grape tannin flavanols with a dissolved-O₂ acetaldehyde
# ethylidene linker (trilinear [acetaldehyde]·[tannin]², the D-84 self form + the D-80 acetaldehyde
# factor), softening astringency WITHOUT anthocyanin. It reuses the D-80 split-ledger carbon capture
# (acetaldehyde → shared ethyl_bridge slot, own y_acetaldehyde_per_tannin) but deposits NO pigment
# (colourless tannin–tannin polymer). These tests pin the closed form + carbon-exact split, the
# [tannin]² shape, the free-acetaldehyde SO₂ read, the gates, the wine-only no-op, NON-TRIVIAL
# carbon closure end-to-end, the anthocyanin-free (colourless) softening, and tiers.


def test_tannin_ethyl_metadata():
    p = TanninEthylTanninCondensation()
    assert p.name == "tannin_ethyl_tannin_condensation"
    assert p.tier is Tier.SPECULATIVE
    # Consumes off-ledger tannin (pure sink, no destination slot) + on-ledger acetaldehyde, whose
    # carbon it books into the shared on-ledger ethyl_bridge slot. Touches NO anthocyanin and NO
    # polymeric_pigment (a colourless tannin–tannin polymer — the colour difference from D-80).
    assert set(p.touches) == {"acetaldehyde", "ethyl_bridge", "tannin"}
    assert "anthocyanin" not in p.touches and "polymeric_pigment" not in p.touches
    # Its OWN acetaldehyde yield (distinct from D-80's y_acetaldehyde_per_anthocyanin).
    assert set(p.reads) == {
        "k_tannin_ethyl_tannin",
        "E_a_tannin_ethyl_tannin",
        "y_acetaldehyde_per_tannin",
        "T_ref",
    }


def test_tannin_ethyl_closed_form(poly_params):
    # r = k·f(T)·[acetaldehyde]·[tannin]² (trilinear, second-order in tannin); anchored on tannin.
    # Verify the tannin draw, the acetaldehyde drawdown at the OWN yield, and the carbon-exact
    # acetaldehyde→ethyl_bridge split (release at cf(acetaldehyde), redeposit at cf(ethylidene)).
    # No SO₂ ⇒ reads total acetaldehyde.
    schema = wine_schema()
    t = 298.15  # off T_ref so the Arrhenius factor bites
    acet, tannin = 0.05, 2.0
    y = _aged_wine(schema, ester=0.0, t=t, acetaldehyde=acet, anthocyanin=0.3, tannin=tannin)
    d = TanninEthylTanninCondensation().derivatives(0.0, y, schema, poly_params)

    f_t = arrhenius_factor(t, poly_params["E_a_tannin_ethyl_tannin"], poly_params["T_ref"])
    k = poly_params["k_tannin_ethyl_tannin"]
    y_acet = poly_params["y_acetaldehyde_per_tannin"]
    r = k * f_t * acet * tannin * tannin
    acet_consumed = y_acet * r
    assert schema.get(d, "tannin") == pytest.approx(-r)
    assert schema.get(d, "acetaldehyde") == pytest.approx(-acet_consumed)
    # The split-ledger capture: acetaldehyde carbon released == carbon deposited into ethyl_bridge.
    expected_bridge = acet_consumed * _ACET_C / _ETHYLIDENE_C
    assert schema.get(d, "ethyl_bridge") == pytest.approx(expected_bridge)
    assert acet_consumed * _ACET_C == pytest.approx(schema.get(d, "ethyl_bridge") * _ETHYLIDENE_C)
    # NO pigment, NO anthocyanin, NO colour (the D-80 colour difference), and nothing else moves.
    for var in (
        "anthocyanin",
        "polymeric_pigment",
        "faded_anthocyanin",
        "o2",
        "X",
        "S",
        "E",
        "N",
        "CO2",
        "A420",
        "ellagitannin",
    ):
        assert schema.get(d, var) == 0.0


def test_tannin_ethyl_is_second_order_in_tannin_and_needs_acetaldehyde(poly_params):
    # BIMOLECULAR in tannin (doubling tannin QUADRUPLES the rate) and FIRST-order in acetaldehyde
    # (doubling acetaldehyde doubles it) — the D-84 self form accelerated by the D-80 factor.
    schema = wine_schema()
    p = TanninEthylTanninCondensation()

    def rate(acet: float, tannin: float) -> float:
        y = _aged_wine(schema, ester=0.0, acetaldehyde=acet, anthocyanin=0.3, tannin=tannin)
        return float(schema.get(p.derivatives(0.0, y, schema, poly_params), "tannin"))

    base = rate(0.04, 1.0)
    assert rate(0.04, 2.0) == pytest.approx(4.0 * base)  # 2× tannin ⇒ 4× rate (second-order)
    assert rate(0.08, 1.0) == pytest.approx(2.0 * base)  # 2× acetaldehyde ⇒ 2× rate (first-order)
    assert base < 0.0  # actually bridging


def test_tannin_ethyl_inert_without_tannin_or_acetaldehyde(poly_params):
    # Doubly substrate-gated: no tannin OR no acetaldehyde ⇒ byte-for-byte zero (a no-tannin /
    # reductive un-oxygenated wine is exactly the case without this Process). Undershoot absorbed.
    schema = wine_schema()
    p = TanninEthylTanninCondensation()
    for kw in (
        {"acetaldehyde": 0.0, "tannin": 2.0},
        {"acetaldehyde": 0.05, "tannin": 0.0},
        {"acetaldehyde": 0.05, "tannin": -1e-9},
    ):
        y = _aged_wine(schema, ester=0.0, anthocyanin=0.3, **kw)
        assert np.array_equal(p.derivatives(0.0, y, schema, poly_params), schema.zeros())


def test_tannin_ethyl_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture lacks k_tannin_ethyl_tannin, yet a no-tannin wine returns zero.
    schema = wine_schema()
    y = _aged_wine(schema, ester=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=0.0)
    d = TanninEthylTanninCondensation().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_tannin_ethyl_reads_free_acetaldehyde_under_so2(bridge_params):
    # SO₂-bound acetaldehyde can't bridge (the D-47/D-80 precedent): at the SAME total acetaldehyde
    # a sulfited wine bridges SLOWER than an unsulfited one — SO₂ DELAYS the tannin softening.
    schema = wine_schema()
    p = TanninEthylTanninCondensation()
    acids = {"tartaric": 4.0, "cation_charge": 0.012}  # a real acid state so pH solves
    unsulfited = _aged_wine(
        schema, ester=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.0, **acids
    )
    sulfited = _aged_wine(
        schema, ester=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.05, **acids
    )
    r_unsulf = float(schema.get(p.derivatives(0.0, unsulfited, schema, bridge_params), "tannin"))
    r_sulf = float(schema.get(p.derivatives(0.0, sulfited, schema, bridge_params), "tannin"))
    assert r_sulf < 0.0 and r_unsulf < 0.0
    assert abs(r_sulf) < abs(r_unsulf)  # SO₂ throttles the free-acetaldehyde-driven rate


def test_tannin_ethyl_unsulfited_is_exactly_total_acetaldehyde(bridge_params):
    # The so2_total > 0 guard is EXACT: at so2_total = 0 the rate uses TOTAL acetaldehyde (no pH
    # solve), byte-for-byte the trilinear closed form — an unsulfited run pays no per-RHS brentq.
    schema = wine_schema()
    acet, tannin, t = 0.05, 2.0, 298.15
    y = _aged_wine(
        schema, ester=0.0, t=t, acetaldehyde=acet, anthocyanin=0.3, tannin=tannin, so2_total=0.0
    )
    d = TanninEthylTanninCondensation().derivatives(0.0, y, schema, bridge_params)
    f_t = arrhenius_factor(t, bridge_params["E_a_tannin_ethyl_tannin"], bridge_params["T_ref"])
    r = bridge_params["k_tannin_ethyl_tannin"] * f_t * acet * tannin * tannin
    assert schema.get(d, "tannin") == pytest.approx(-r)


def test_tannin_ethyl_rises_with_temperature(poly_params):
    # Warmer bridges faster (E_a > 0, reaction-scale): the |tannin| draw grows with T.
    schema = wine_schema()
    p = TanninEthylTanninCondensation()
    cold = p.derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=283.15, acetaldehyde=0.05, tannin=2.0),
        schema,
        poly_params,
    )
    warm = p.derivatives(
        0.0,
        _aged_wine(schema, ester=0.0, t=303.15, acetaldehyde=0.05, tannin=2.0),
        schema,
        poly_params,
    )
    assert abs(float(schema.get(warm, "tannin"))) > abs(float(schema.get(cold, "tannin")))


def test_tannin_ethyl_wine_only_noop_on_beer(poly_params):
    # Wine-only (tannin/ethyl_bridge are appended to wine_schema): a hard no-op on beer (aldehyde
    # exists in both media, but the grape/bridge slots do not).
    beer = beer_schema()
    yb = beer.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 100.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    yb[beer.slice("acetaldehyde")] = 0.05
    assert np.array_equal(
        TanninEthylTanninCondensation().derivatives(0.0, yb, beer, poly_params), beer.zeros()
    )


def test_tannin_ethyl_carbon_closes_nontrivially(bridge_store, bridge_params):
    # THE D-85 SPLIT-LEDGER PROOF (the ledger-touching beat). It consumes on-ledger acetaldehyde and
    # books its carbon into the on-ledger ethyl_bridge slot, so total_carbon is flat NON-TRIVIALLY —
    # acetaldehyde↓ exactly cancels ethyl_bridge↑. mass/nitrogen flat too (touches no {S,E,CO2}/N).
    # No anthocyanin here — the softening is purely tannin–ethyl–tannin (colourless).
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninEthylTanninCondensation()], strict=True)
    y0 = _aged_wine(
        schema,
        ester=0.0,
        acetaldehyde=0.05,
        anthocyanin=0.0,  # NO anthocyanin — a pure tannin–ethyl–tannin run
        tannin=3.0,
        tartaric=4.0,
        cation_charge=0.012,
    )
    traj = simulate(ps, params=bridge_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success
    f_c = bridge_store.value("biomass_C_fraction")
    f_n = bridge_store.value("biomass_N_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_conserved(traj, total_mass(schema), label="mass")
    assert_conserved(traj, total_nitrogen(schema, biomass_nitrogen_fraction=f_n), label="nitrogen")
    # NON-trivial: acetaldehyde genuinely fell, ethyl_bridge genuinely rose, carbon lost == gained.
    acet = np.asarray(traj.series("acetaldehyde"), dtype=float)
    bridge = np.asarray(traj.series("ethyl_bridge"), dtype=float)
    tannin = np.asarray(traj.series("tannin"), dtype=float)
    assert acet[-1] < acet[0] - 1e-3  # acetaldehyde consumed
    assert bridge[-1] > bridge[0] + 1e-4  # ethyl_bridge accumulated
    assert tannin[-1] < tannin[0]  # tannin softened (the astringency payoff)
    dC_acet = (acet[-1] - acet[0]) * _ACET_C
    dC_bridge = (bridge[-1] - bridge[0]) * _ETHYLIDENE_C
    assert dC_acet + dC_bridge == pytest.approx(0.0, abs=1e-12)


def test_tannin_ethyl_tier_floored_at_speculative(bridge_store):
    # Speculative in FORM + capped by its speculative params (D-1), for every pool it writes —
    # including the on-ledger acetaldehyde/ethyl_bridge pair.
    schema = wine_schema()
    ps = ProcessSet(schema, [TanninEthylTanninCondensation()])
    for pool in ("tannin", "acetaldehyde", "ethyl_bridge"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, bridge_store.tier_map()) is Tier.SPECULATIVE


def test_fading_so2_protects_colour_emergently(fade_params):
    # THE EMERGENT PAYOFF: SO₂ protects the colour with NOTHING scripted. Two oxygenated reds
    # (SulfiteOxidation + AnthocyaninFading on the shared o2 pool), identical but for the SO₂ dose:
    # the sulfited wine scavenges o2 via D-72, leaving LESS o2 to fade the anthocyanin, so it keeps
    # MORE colour (less faded). This falls out of the shared o2 split — no SO₂ term in the rate.
    schema = wine_schema()
    ps = ProcessSet(schema, [SulfiteOxidation(), AnthocyaninFading()], strict=True)

    def run(so2_0: float) -> Trajectory:
        y0 = _sulfited_wine(schema, so2=so2_0, o2=0.05, t=298.15, anthocyanin=0.3)
        traj = simulate(ps, params=fade_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
        assert traj.success
        return traj

    protected = run(0.08)  # a real SO₂ dose (~80 mg/L)
    unprotected = run(0.0)  # no SO₂ — o2 goes entirely to fading (and ethanol oxidation)
    faded_protected = float(protected.series("faded_anthocyanin")[-1])
    faded_unprotected = float(unprotected.series("faded_anthocyanin")[-1])
    # Emergent: SO₂ diverts o2 to bisulfite oxidation, so LESS anthocyanin fades — more colour kept.
    assert faded_protected < faded_unprotected
    # And the surviving colour is correspondingly higher with SO₂ (the colour-stability payoff).
    assert float(color_series(protected)[-1]) > float(color_series(unprotected)[-1])


def test_fading_tier_floored_at_speculative(poly_store):
    # Speculative in FORM (Tier-3 frontier) and capped by its speculative fade params (D-1), for
    # every pool it writes — the shared o2 pool and the grape anthocyanin/faded pair.
    schema = wine_schema()
    ps = ProcessSet(schema, [AnthocyaninFading()])
    for pool in ("o2", "anthocyanin", "faded_anthocyanin"):
        assert ps.tier_of(pool) is Tier.SPECULATIVE
        assert ps.tier_of(pool, poly_store.tier_map()) is Tier.SPECULATIVE


def test_hydrolysis_returns_the_label_with_the_c5_so_aging_cannot_dilute_the_enrichment(params):
    """D-115 - the label comes back with the carbon, at the ESTER's own fraction.

    Hydrolysis is the exact reverse of the acetylation, so an ester molecule that was
    valine-derived returns a valine-derived alcohol molecule, mole for mole. Two things have to
    hold and neither is implied by the other:

    * **the ester side is non-fractionating** - the tracer falls at the pool's own fraction, so
      decaying the pool leaves its enrichment where it was;
    * **the alcohol side carries that fraction across** - the returned alcohol is credited as
      labelled in the ester's proportion. **This is the half the wiring exists for.** Omit it and
      an aging segment silently dilutes the alcohol pool's enrichment with returned molecules
      booked as unlabelled - a drift that no conservation test could ever catch, because a tracer
      slot carries no carbon weight by construction (the D-89/D-90 family).

    Deliberately debited at ``f_ester`` and **not** at the alcohol's fraction: using the alcohol's
    would assume the ester tracks it, which is the very thing the two-slot design exists to
    measure rather than assert (the D-98/D-108 vacuity trap, relocated into the RHS).

    Covered here because nothing else can: this Process is inert through every fermentation-phase
    run in the suite, so the enrichment tests in ``test_fusel_reroute.py`` never exercise it.
    """
    schema = wine_schema()
    ester = 0.1
    for fraction in (0.0, 0.25, 1.0):
        y = _aged_wine(schema, ester=ester, t=298.15)
        y[schema.slice("isoamyl_acetate_valine")] = ester * fraction
        d = EsterHydrolysis().derivatives(0.0, y, schema, params)

        decayed = -float(d[schema.slice("isoamyl_acetate")][0])
        assert decayed > 0.0, "vacuous: no hydrolysis at this state"

        # (a) the ester side is non-fractionating.
        assert -float(d[schema.slice("isoamyl_acetate_valine")][0]) == pytest.approx(
            decayed * fraction, rel=1e-12
        ), "hydrolysis must consume labelled and unlabelled ester in the pool's own proportion"

        # (b) the returned ALCOHOL carries exactly the ester's fraction - the half that keeps an
        # aging segment from diluting the alcohol pool's enrichment.
        returned = float(d[schema.slice("isoamyl_alcohol")][0])
        returned_label = float(d[schema.slice("isoamyl_alcohol_valine")][0])
        assert returned > 0.0, "vacuous: no alcohol is being returned"
        assert returned_label == pytest.approx(returned * fraction, rel=1e-12), (
            "the returned alcohol's labelled fraction must equal the ester's - the C5 skeleton "
            "comes back as a unit, so the molecule fraction crosses the reaction unchanged"
        )
        # ...and mole for mole with the ester consumed, which is what makes the two independent.
        assert returned / M_ISOAMYL_OH == pytest.approx(decayed / M_ISOAMYL_ACETATE, rel=1e-12)
