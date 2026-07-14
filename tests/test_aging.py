"""Tests for the Tier-3 aging Process :class:`EsterHydrolysis` (decision D-69).

The first §4.1 aging Process: a first-order **net decay** of the lumped ``esters`` pool
toward a lower equilibrium floor ``esters_eq`` (young fruity acetate esters hydrolyse and
fade with age), warmed by an Arrhenius factor (warmer ages faster), routing the released
ester carbon **5:2** into ``fusels`` (isoamyl alcohol, the alcohol product) and ``Byp``
(succinic-stand-in acetic acid, the acid product). These tests pin the closed-form
derivative and the exact 5:2 split; prove the properties the aging axis requires — **net
decay toward equilibrium** (zero at/below ``esters_eq``, not decay-to-zero), **warmer ⇒
faster**, and an **on-ledger carbon transfer that closes ``total_carbon`` to machine
precision** (the D-68 "conservation is back in force" requirement, unlike the D-67 readout);
check the solver-undershoot guards; and confirm the tier floors at speculative and the
Process touches only ``esters``/``fusels``/``Byp`` (no ``S``/``E``/``CO2`` — aging draws no
sugar). The scenario-level aging-phase wiring (the ``age N months`` verb + the reconfigure
enable) is D-70; here the Process is exercised directly via a hand-built ``ProcessSet`` (the
D-64 loss-Process pattern), off the fermentation ProcessSet so isolability is preserved.
"""

import numpy as np
import pytest

from fermentation.analysis import astringency_series, color_series, polymeric_pigment_series
from fermentation.core.acidbase import bisulfite_so2_at_ph, ph_of_state
from fermentation.core.chemistry import (
    M_2_METHYLBUTANAL,
    M_2_METHYLPROPANAL,
    M_3_METHYLBUTANAL,
    M_ACETALDEHYDE,
    M_CO2,
    M_ETHANOL,
    M_METHIONAL,
    M_O2,
    M_PHENYLACETALDEHYDE,
    M_SO2,
    M_SOTOLON,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics import (
    AcetaldehydeBridgedCondensation,
    AnthocyaninFading,
    Caramelization,
    EllagitanninOxidation,
    EsterHydrolysis,
    MaillardStrecker,
    OakExtraction,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    StreckerDegradation,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
    TanninEthylTanninCondensation,
    TanninSelfPolymerization,
    ThermalAnthocyaninFade,
    arrhenius_factor,
)
from fermentation.core.kinetics.aging import _MAILLARD_PRODUCTS, _SO2_PER_O2
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import Trajectory, simulate
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_mass,
    total_nitrogen,
)

#: Carbon fractions of the three pools the transfer touches (mirror the Process constants).
_ESTER_C = carbon_mass_fraction("ethyl_acetate")
_FUSEL_C = carbon_mass_fraction("isoamyl_alcohol")
_BYP_C = carbon_mass_fraction("succinic_acid")
#: The 5:2 split, from the isoamyl-acetate stand-in reaction (isoamyl alcohol 5 C : acetic 2 C).
_FUSEL_SHARE = 5.0 / 7.0
_BYP_SHARE = 2.0 / 7.0
#: Carbon fractions of the two pools the oxidation transfer moves carbon between (E → acetaldehyde).
_ETHANOL_C = carbon_mass_fraction("ethanol")
_ACET_C = carbon_mass_fraction("acetaldehyde")
#: Carbon fraction of the ethylidene bridge (C2H4) — the D-80 split-ledger on-ledger capture
#: species.
_ETHYLIDENE_C = carbon_mass_fraction("ethylidene")
#: Carbon fractions of the Strecker pools + the amino-acid source (D-75), for the closure checks.
_METHIONAL_C = carbon_mass_fraction("methional")
_PHENYLACET_C = carbon_mass_fraction("phenylacetaldehyde")
_CO2_C = carbon_mass_fraction("CO2")
_AA_C = carbon_mass_fraction(AMINO_ACID_SPECIES)
_AA_N = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
# The four non-oxidative THERMAL Strecker aldehyde/sotolon carbon fractions (decision D-87).
_2MB_C = carbon_mass_fraction("2_methylbutanal")
_3MB_C = carbon_mass_fraction("3_methylbutanal")
_2MP_C = carbon_mass_fraction("2_methylpropanal")
_SOTOLON_C = carbon_mass_fraction("sotolon")
# Caramelization (decision D-88): the sugar (glucose) and melanoidin-carbon-park carbon fractions.
_GLUCOSE_C = carbon_mass_fraction("glucose")
_MELANOIDIN_C = carbon_mass_fraction("melanoidin")


@pytest.fixture
def store():
    # Real wine parameters (T_ref, biomass_C_fraction, ...) MERGED with the aging.yaml
    # hydrolysis constants — the shared, medium-agnostic aging file (D-69). This mirrors the
    # compile seam D-70 will wire; here the test loads it directly.
    return load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "aging.yaml"
    )


@pytest.fixture
def params(store):
    return store.resolve()


def _aged_wine(schema: StateSchema, *, esters: float = 0.1, t: float = 293.15, **kw) -> FloatArray:
    """A finished, racked wine at the start of aging: yeast gone (X=0), dry (S=0), with the
    liquid ``esters`` pool pre-loaded (nothing produces it during aging — the Process only
    decays it). ``fusels``/``Byp`` default to 0 so their aging gains are unambiguous."""
    y = schema.pack({"X": 0.0, "S": [0.0], "E": 100.0, "N": 0.0, "T": t, "CO2": 0.0})
    y[schema.slice("esters")] = esters
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


# -- metadata -----------------------------------------------------------------


def test_metadata():
    p = EsterHydrolysis()
    assert p.name == "ester_hydrolysis"
    # Speculative: the aging axis is the Tier-3 frontier (form sourced, magnitudes estimated).
    assert p.tier is Tier.SPECULATIVE
    # An on-ledger inter-pool transfer: decays esters, routes carbon to the alcohol (fusels)
    # and acid (Byp) products — never S/E/CO2 (aging draws no sugar, unlike the M2 producers).
    assert set(p.touches) == {"esters", "fusels", "Byp"}
    assert set(p.reads) == {
        "k_ester_hydrolysis",
        "E_a_ester_hydrolysis",
        "esters_eq",
        "T_ref",
    }


# -- closed form & the 5:2 split ----------------------------------------------


def test_derivative_matches_closed_form(params):
    schema = wine_schema()
    esters, t = 0.1, 298.15  # off T_ref so the Arrhenius factor bites
    y = _aged_wine(schema, esters=esters, t=t)
    d = EsterHydrolysis().derivatives(0.0, y, schema, params)

    f_t = arrhenius_factor(t, params["E_a_ester_hydrolysis"], params["T_ref"])
    rate = params["k_ester_hydrolysis"] * f_t * (esters - params["esters_eq"])
    carbon_released = rate * _ESTER_C

    assert schema.get(d, "esters") == pytest.approx(-rate)
    # 5:2 split of the released carbon, re-deposited via each product's own carbon fraction.
    assert schema.get(d, "fusels") == pytest.approx(_FUSEL_SHARE * carbon_released / _FUSEL_C)
    assert schema.get(d, "Byp") == pytest.approx(_BYP_SHARE * carbon_released / _BYP_C)
    # Aging touches nothing else — no sugar draw, no ethanol/CO2, no biomass.
    for var in ("X", "S", "E", "N", "CO2"):
        assert schema.get(d, var) == 0.0


def test_carbon_closes_per_rhs(params):
    # THE D-68 "conservation is back in force" invariant: the ester carbon lost equals the
    # fusel + Byp carbon gained, to machine precision — a pure on-ledger inter-pool transfer
    # (no S involvement), so total_carbon closes for ANY split summing to 1 (here 5:2).
    schema = wine_schema()
    d = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, esters=0.1, t=298.15), schema, params)
    carbon_residual = (
        schema.get(d, "esters") * _ESTER_C
        + schema.get(d, "fusels") * _FUSEL_C
        + schema.get(d, "Byp") * _BYP_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-15)


def test_split_is_five_to_two_by_carbon(params):
    # The carbon (not mass) partition is exactly 5:2 — the isoamyl-acetate stand-in ratio
    # (isoamyl alcohol 5 C : acetic acid 2 C), the advisor-settled crux (D-69). Verified as
    # carbon so it is independent of the pools' differing mass weightings.
    schema = wine_schema()
    d = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, esters=0.1), schema, params)
    fusel_carbon = schema.get(d, "fusels") * _FUSEL_C
    byp_carbon = schema.get(d, "Byp") * _BYP_C
    assert fusel_carbon / byp_carbon == pytest.approx(5.0 / 2.0)
    # Both products gain, and the fusel share is the larger (5/7) — the stronger fusel-OAV
    # rise the owner asked for (D-68 fork 2), plus the VA/pH-drifting Byp acid product.
    assert fusel_carbon > byp_carbon > 0.0


# -- net decay toward equilibrium (not to zero) -------------------------------


def test_zero_at_and_below_equilibrium(params):
    # Net decay toward a LOWER floor, not decay-to-zero (D-68): at or below esters_eq the
    # rate is zero (the reverse ester-formation half is the deferred bidirectional term).
    schema = wine_schema()
    eq = params["esters_eq"]
    at_eq = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, esters=eq), schema, params)
    below = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, esters=eq * 0.5), schema, params)
    assert np.array_equal(at_eq, schema.zeros())
    assert np.array_equal(below, schema.zeros())


def test_decays_only_the_excess_above_equilibrium(params):
    # The rate is proportional to (esters - esters_eq), so a pool twice as far above the floor
    # decays twice as fast — the linear approach to equilibrium (Ramey & Ough first-order form).
    schema = wine_schema()
    eq = params["esters_eq"]
    near = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, esters=eq + 0.02), schema, params)
    far = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, esters=eq + 0.04), schema, params)
    assert schema.get(far, "esters") == pytest.approx(2.0 * schema.get(near, "esters"))
    assert schema.get(far, "esters") < schema.get(near, "esters") < 0.0  # both decaying


def test_solver_undershoot_does_not_create_pools(params):
    # A solver undershoot (esters < 0) must not flip max(0, ...) into spurious production of
    # fusels/Byp (or negative decay). esters_eq > 0 makes the excess negative ⇒ clamped to 0.
    schema = wine_schema()
    d = EsterHydrolysis().derivatives(0.0, _aged_wine(schema, esters=-1e-6), schema, params)
    assert np.array_equal(d, schema.zeros())


# -- temperature direction (warmer ages faster) -------------------------------


def test_rises_with_temperature(params):
    # The sourced ordering (E_a_ester_hydrolysis > 0): warmer storage hydrolyses the esters
    # faster — why warm cellars age wine faster and cold storage preserves fruity esters.
    schema = wine_schema()
    cold = EsterHydrolysis().derivatives(
        0.0, _aged_wine(schema, esters=0.1, t=283.15), schema, params
    )
    warm = EsterHydrolysis().derivatives(
        0.0, _aged_wine(schema, esters=0.1, t=303.15), schema, params
    )
    # Faster decay (more negative) and a correspondingly larger fusel/Byp gain when warm.
    assert schema.get(warm, "esters") < schema.get(cold, "esters") < 0.0
    assert schema.get(warm, "fusels") > schema.get(cold, "fusels") > 0.0


def test_factor_is_one_at_reference_temperature(params):
    # At T_ref the Arrhenius factor is exactly 1, so the rate is the bare first-order term.
    schema = wine_schema()
    esters = 0.1
    d = EsterHydrolysis().derivatives(
        0.0, _aged_wine(schema, esters=esters, t=params["T_ref"]), schema, params
    )
    expected = params["k_ester_hydrolysis"] * (esters - params["esters_eq"])
    assert schema.get(d, "esters") == pytest.approx(-expected)


# -- integrated aging segment (conservation + direction) ----------------------


def test_integrated_aging_closes_carbon_and_fades_esters(params, store):
    # Run a long aging segment (a racked, dry wine — X=0, S=0) with ONLY EsterHydrolysis
    # active, under the strict touches contract. Over the aging span the esters pool fades,
    # the fusels and Byp pools rise, and total_carbon closes to machine precision — the pure
    # on-ledger transfer of the D-68 aging axis (no sugar drawn, no ethanol touched).
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()], strict=True)
    esters0 = 0.1
    y0 = _aged_wine(schema, esters=esters0, t=293.15)
    # ~1 year of aging (large steps are fine — no ferment stiffness; the §7 slow-phase point).
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message

    # The fruity esters fade (toward, not past, the equilibrium floor) and the products rise.
    esters_end = float(traj.series("esters")[-1])
    assert params["esters_eq"] <= esters_end < esters0
    assert float(traj.series("fusels")[-1]) > 0.0
    assert float(traj.series("Byp")[-1]) > 0.0
    # Non-negative pools and machine-precision carbon closure (X=0 throughout, so the biomass
    # term is inert; the invariant is the esters → fusels + Byp inter-pool transfer).
    assert_nonnegative(traj, ("esters", "fusels", "Byp"), atol=1e-12)
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    # Carbon is the invariant; mass carries a small (~4.5%) documented stand-in gap (aging.yaml):
    # splitting a carbon-exact budget across pools with heterogeneous fixed weightings is not
    # mass-conserving (real hydrolysis consumes untracked water — the D-8/D-16/D-26 precedent).
    # total_mass weights only {S, E, CO2}, NONE of which this Process touches, so the gap is
    # scoped OUT by construction: the validated-core mass check stays flat through the aging span.
    assert_conserved(traj, total_mass(schema), rtol=1e-9, atol=1e-9, label="mass")


def test_isolable_from_the_core_when_below_equilibrium(params):
    # Isolability corner (prime directive #3): with the esters pool below esters_eq the Process
    # contributes exactly nothing, so an aging segment on an ester-poor wine is byte-for-byte
    # the no-aging state — the Process cannot create aroma out of an empty pool.
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()], strict=True)
    y = _aged_wine(schema, esters=params["esters_eq"] * 0.5)
    assert np.array_equal(ps.total_derivatives(0.0, y, params), schema.zeros())


def test_integrated_aging_closes_carbon_beer_multislot(store):
    # The multi-slot (beer) counterpart: the aging transfer is sugar-free so the 3-slot S
    # vector is irrelevant to it, but running the strict ProcessSet on beer_schema proves the
    # Process is medium-agnostic (esters/fusels/Byp exist in both) and closes carbon there too.
    beer = load_parameters(
        default_data_dir() / "beer_generic.yaml", default_data_dir() / "aging.yaml"
    )
    params = beer.resolve()
    schema = beer_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()], strict=True)
    y0 = schema.pack({"X": 0.0, "S": [0.0, 0.0, 0.0], "E": 40.0, "N": 0.0, "T": 293.15, "CO2": 0.0})
    y0[schema.slice("esters")] = 0.08
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 24.0 * 180.0))
    assert traj.success, traj.message
    assert float(traj.series("esters")[-1]) < 0.08
    f_c = beer.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


# -- tier propagation ---------------------------------------------------------


def test_tier_floored_at_speculative(store):
    # The aging Process is speculative in FORM (Tier-3 frontier), so every pool it writes is
    # speculative even before parameters cap it — and folding in the (speculative) aging
    # parameter tiers keeps it there. Non-vacuous: esters/fusels/Byp are all speculative.
    schema = wine_schema()
    ps = ProcessSet(schema, [EsterHydrolysis()])
    for pool in ("esters", "fusels", "Byp"):
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
    y = _aged_wine(schema, esters=0.0, t=t, o2=o2)  # esters=0 so only oxidation moves anything
    d = OxidativeAcetaldehyde().derivatives(0.0, y, schema, params)

    f_t = arrhenius_factor(t, params["E_a_ethanol_oxidation"], params["T_ref"])
    r_o2 = params["k_ethanol_oxidation"] * f_t * o2
    acet_rate = params["y_acetaldehyde_per_o2"] * (r_o2 / M_O2) * M_ACETALDEHYDE

    assert schema.get(d, "o2") == pytest.approx(-r_o2)
    assert schema.get(d, "acetaldehyde") == pytest.approx(acet_rate)
    # Carbon-exact C2 borrow from ethanol (the reduction reversed).
    assert schema.get(d, "E") == pytest.approx(-acet_rate * M_ETHANOL / M_ACETALDEHYDE)
    # Oxidation touches nothing else — no sugar, no CO2, no esters/fusels/Byp, no biomass.
    for var in ("X", "S", "N", "CO2", "esters", "fusels", "Byp"):
        assert schema.get(d, var) == 0.0


def test_oxidation_carbon_closes_per_rhs(params):
    # O₂ is OFF every ledger, so the only on-ledger movement is E → acetaldehyde, both C2 — the
    # carbon lost from ethanol equals the carbon gained as acetaldehyde, to machine precision.
    schema = wine_schema()
    d = OxidativeAcetaldehyde().derivatives(
        0.0, _aged_wine(schema, esters=0.0, t=298.15, o2=0.03), schema, params
    )
    carbon_residual = schema.get(d, "E") * _ETHANOL_C + schema.get(d, "acetaldehyde") * _ACET_C
    assert carbon_residual == pytest.approx(0.0, abs=1e-15)


def test_oxidation_inert_without_oxygen(params):
    # Reductive aging (screwcap/inert) + the exact isolability guard: with no dissolved O₂ the
    # Process contributes byte-for-byte zero, so a begin_aging run with no add_oxygen is the
    # ester-hydrolysis-only case. Cannot oxidise ethanol out of an empty O₂ pool.
    schema = wine_schema()
    ps = ProcessSet(schema, [OxidativeAcetaldehyde()], strict=True)
    y = _aged_wine(schema, esters=0.0, o2=0.0)
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
    y0 = _aged_wine(schema, esters=0.0, t=298.15, o2=o2_0)
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
    y = _aged_wine(schema, esters=0.0, t=t, so2_total=so2, o2=o2, tartaric=4.0, cation_charge=0.012)
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
    for var in ("X", "S", "E", "N", "CO2", "acetaldehyde", "esters", "fusels", "Byp"):
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
# ``k_browning`` > ``k_ethanol_oxidation``), Arrhenius warmer-faster. Touches only ``o2`` +
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
    # borrow).
    assert set(p.touches) == {"o2", "A420"}
    assert set(p.reads) == {"k_browning", "E_a_browning", "y_a420_per_o2", "T_ref"}


def test_browning_matches_closed_form(params):
    schema = wine_schema()
    o2, t = 0.03, 298.15  # off T_ref so the Arrhenius factor bites
    y = _aged_wine(schema, esters=0.0, t=t, o2=o2)
    d = PhenolicBrowning().derivatives(0.0, y, schema, params)
    f_t = arrhenius_factor(t, params["E_a_browning"], params["T_ref"])
    r_o2 = params["k_browning"] * f_t * o2
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
    # it draws a larger O₂ rate than ethanol oxidation (k_browning > k_ethanol_oxidation), and the
    # two
    # shares sum to the calibrated always-on total (5.0e-4) that holds the O₂-depletion timescale.
    assert params["k_browning"] > params["k_ethanol_oxidation"]
    assert params["k_browning"] + params["k_ethanol_oxidation"] == pytest.approx(5.0e-4)
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, o2=0.03)
    brown = PhenolicBrowning().derivatives(0.0, y, schema, params)
    ethanol = OxidativeAcetaldehyde().derivatives(0.0, y, schema, params)
    assert -schema.get(brown, "o2") > -schema.get(ethanol, "o2") > 0.0


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
    r_o2 = params["k_browning"] * f_t * o2
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
    y0 = _aged_wine(schema, esters=0.0, t=298.15, o2=o2_0)
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
    # only once both have plateaued does the clean partition ratio k_ethanol/(k_ethanol+k_browning)
    # hold (in finite time the slower run lags its ceiling and the ratio reads high).
    span = (0.0, 24.0 * 365.0 * 5.0)

    def run(processes: list[Process]) -> Trajectory:
        ps = ProcessSet(schema, processes, strict=True)
        y0 = _aged_wine(schema, esters=0.0, t=298.15, o2=o2_0)
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
    # ~k_browning /
    # (k_browning + k_ethanol) of the O₂, so acetaldehyde falls toward the ethanol share ~40%.
    assert acet_diverted < acet_alone
    assert float(with_browning.series("A420")[-1]) > 0.0
    assert float(ethanol_only.series("A420")[-1]) == 0.0
    share_ethanol = params["k_ethanol_oxidation"] / (
        params["k_ethanol_oxidation"] + params["k_browning"]
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

_STRECKER_TOUCHES = {"o2", "methional", "phenylacetaldehyde", "CO2", "amino_acids", "N"}


def _strecker_wine(
    schema: StateSchema, *, aa: float = 0.05, o2: float = 0.03, t: float = 293.15, **kw
) -> FloatArray:
    """A finished, racked wine at the start of aging with dosed amino acids + O2 — the two Strecker
    substrates. ``esters`` defaults to 0 (irrelevant here); any extra pool via kwargs."""
    y = _aged_wine(schema, esters=0.0, t=t, o2=o2, amino_acids=aa)
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def _strecker_closed_form(
    schema: StateSchema, params: dict[str, float], y: FloatArray, t: float
) -> dict[str, float]:
    """The Process's own algebra, recomputed independently for the closed-form assertions."""
    o2 = float(y[schema.slice("o2")][0])
    aa = float(y[schema.slice("amino_acids")][0])
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
    assert set(p.reads) == {
        "k_strecker",
        "E_a_strecker",
        "y_strecker_per_o2",
        "f_methional",
        "K_amino_acids",
        "T_ref",
    }


def test_strecker_matches_closed_form(params):
    schema = wine_schema()
    aa, o2, t = 0.05, 0.03, 298.15  # off T_ref so the Arrhenius factor bites
    y = _strecker_wine(schema, aa=aa, o2=o2, t=t)
    d = StreckerDegradation().derivatives(0.0, y, schema, params)
    cf = _strecker_closed_form(schema, params, y, t)

    assert cf["methional"] > 0.0  # the products are live (guards against a vacuous pass)
    assert schema.get(d, "o2") == pytest.approx(cf["o2"])
    assert schema.get(d, "methional") == pytest.approx(cf["methional"])
    assert schema.get(d, "phenylacetaldehyde") == pytest.approx(cf["phenylacetaldehyde"])
    assert schema.get(d, "CO2") == pytest.approx(cf["CO2"])
    # amino_acids drawn sized to the product carbon; nitrogen deaminated to N.
    product_carbon = (
        cf["methional"] * _METHIONAL_C
        + cf["phenylacetaldehyde"] * _PHENYLACET_C
        + cf["CO2"] * _CO2_C
    )
    aa_mass = product_carbon / _AA_C
    assert schema.get(d, "amino_acids") == pytest.approx(-aa_mass)
    assert schema.get(d, "N") == pytest.approx(aa_mass * _AA_N)
    # Touches nothing else — no ethanol/esters/fusels/acetaldehyde, no sugar, no biomass.
    for var in ("X", "S", "E", "esters", "fusels", "Byp", "acetaldehyde"):
        assert schema.get(d, var) == 0.0


def test_strecker_carbon_closes_per_rhs(params):
    # CARBON closes to machine precision: the arginine carbon leaving amino_acids equals the carbon
    # entering methional + phenylacetaldehyde + CO2 (the draw is sized to match) — a pure on-ledger
    # transfer, off-ledger o2 aside.
    schema = wine_schema()
    d = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, t=298.15), schema, params)
    carbon_residual = (
        schema.get(d, "methional") * _METHIONAL_C
        + schema.get(d, "phenylacetaldehyde") * _PHENYLACET_C
        + schema.get(d, "CO2") * _CO2_C
        + schema.get(d, "amino_acids") * _AA_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)


def test_strecker_nitrogen_closes_per_rhs(params):
    # NITROGEN closes: all the arginine nitrogen leaving amino_acids lands in the N pool (the
    # aldehydes are nitrogen-free — the deamination, the D-45 mercaptan idiom).
    schema = wine_schema()
    d = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, t=298.15), schema, params)
    nitrogen_residual = schema.get(d, "amino_acids") * _AA_N + schema.get(d, "N") * 1.0
    assert nitrogen_residual == pytest.approx(0.0, abs=1e-18)


def test_strecker_inert_without_oxygen(params):
    # No oxidant ⇒ no Strecker: a reductive begin_aging (no add_oxygen) is byte-for-byte the case
    # without this Process. Also the o2<0 solver-undershoot guard.
    schema = wine_schema()
    ps = ProcessSet(schema, [StreckerDegradation()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _strecker_wine(schema, o2=0.0), params), schema.zeros()
    )
    assert np.array_equal(
        StreckerDegradation().derivatives(0.0, _strecker_wine(schema, o2=-1e-6), schema, params),
        schema.zeros(),
    )


def test_strecker_inert_without_amino_acids(params):
    # No amino acids ⇒ no Strecker: an amino-acid-free aging is byte-for-byte the case without this
    # Process (the substrate gate that makes it ADD ON TOP with no re-baseline, D-75). Also aa<0.
    schema = wine_schema()
    ps = ProcessSet(schema, [StreckerDegradation()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _strecker_wine(schema, aa=0.0), params), schema.zeros()
    )
    assert np.array_equal(
        StreckerDegradation().derivatives(0.0, _strecker_wine(schema, aa=-1e-6), schema, params),
        schema.zeros(),
    )


def test_strecker_first_order_in_oxygen(params):
    # First-order in the O2 pool (at fixed amino acids, so the availability gate is held constant):
    # doubling [o2] doubles the instantaneous O2 draw and every product rate.
    schema = wine_schema()
    base = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, o2=0.02), schema, params)
    dbl = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, o2=0.04), schema, params)
    for pool in ("o2", "methional", "phenylacetaldehyde", "CO2", "amino_acids", "N"):
        assert schema.get(dbl, pool) == pytest.approx(2.0 * schema.get(base, pool))


def test_strecker_availability_gate_saturates(params):
    # The amino-acid availability gate aa/(K+aa) throttles the draw at low aa and SATURATES toward a
    # ceiling at high aa (the smooth swap/reroute gate, D-33). At aa >> K the rate approaches the
    # ungated k*f*[o2]; at aa == K it is ~half that — a monotone, saturating aa dependence.
    schema = wine_schema()
    k = params["K_amino_acids"]
    low = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, aa=0.1 * k), schema, params)
    mid = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, aa=k), schema, params)
    high = StreckerDegradation().derivatives(
        0.0, _strecker_wine(schema, aa=100.0 * k), schema, params
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
    d = StreckerDegradation().derivatives(0.0, _strecker_wine(schema), schema, params)
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
    cold = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, t=283.15), schema, params)
    warm = StreckerDegradation().derivatives(0.0, _strecker_wine(schema, t=303.15), schema, params)
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
    y0 = _strecker_wine(schema, aa=0.05, o2=0.04, t=298.15)
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

_MAILLARD_TOUCHES = {
    "methional",
    "phenylacetaldehyde",
    "2_methylbutanal",
    "3_methylbutanal",
    "2_methylpropanal",
    "sotolon",
    "CO2",
    "amino_acids",
    "N",
}
# Per-product carbon fraction, keyed by pool (for the closed-form carbon accounting).
_MAILLARD_C = {
    "methional": _METHIONAL_C,
    "phenylacetaldehyde": _PHENYLACET_C,
    "2_methylbutanal": _2MB_C,
    "3_methylbutanal": _3MB_C,
    "2_methylpropanal": _2MP_C,
    "sotolon": _SOTOLON_C,
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
    schema: StateSchema, *, aa: float = 0.3, s: float = 80.0, t: float = 298.15, **kw
) -> FloatArray:
    """A finished, SEALED (o2 = 0 — the whole point) SWEET wine at the start of aging: residual
    sugar ``s`` (the dicarbonyl driver) + dosed ``amino_acids``, warm. Any extra pool via kwargs."""
    y = _aged_wine(schema, esters=0.0, t=t, amino_acids=aa)
    y[schema.slice("S")] = s
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def _maillard_closed_form(
    schema: StateSchema, params: dict[str, float], y: FloatArray, t: float
) -> dict[str, float]:
    """The Process's own algebra, recomputed independently for the closed-form assertions."""
    aa = float(y[schema.slice("amino_acids")][0])
    s_total = float(y[schema.slice("S")].sum())
    gate = aa / (params["K_amino_acids"] + aa)
    f_t = arrhenius_factor(t, params["E_a_maillard_strecker"], params["T_ref"])
    n_ald = params["k_maillard_strecker"] * f_t * s_total * gate
    weights = [params[wname] for (_, _, wname, _) in _MAILLARD_PRODUCTS]
    w_sum = sum(weights)
    out: dict[str, float] = {}
    co2_mol = 0.0
    for (pool, m_i, _wn, decarb), w_i in zip(_MAILLARD_PRODUCTS, weights, strict=True):
        n_i = (w_i / w_sum) * n_ald
        out[pool] = n_i * m_i
        if decarb:
            co2_mol += n_i
    out["CO2"] = co2_mol * M_CO2
    return out


def test_maillard_metadata():
    p = MaillardStrecker()
    assert p.name == "maillard_strecker"
    # Speculative: the aging axis is the Tier-3 frontier (form sourced, magnitudes estimated).
    assert p.tier is Tier.SPECULATIVE
    # Books the six thermal products + the decarboxylation CO2, drawing carbon from amino_acids and
    # deaminating the nitrogen to N. Touches those nine and nothing else — NO o2, and S is a
    # read-only driver (not in touches).
    assert set(p.touches) == _MAILLARD_TOUCHES
    assert "o2" not in p.touches
    assert "S" not in p.touches
    assert set(p.reads) == {
        "k_maillard_strecker",
        "E_a_maillard_strecker",
        "w_maillard_methional",
        "w_maillard_phenylacetaldehyde",
        "w_maillard_2_methylbutanal",
        "w_maillard_3_methylbutanal",
        "w_maillard_2_methylpropanal",
        "w_maillard_sotolon",
        "K_amino_acids",
        "T_ref",
    }


def test_maillard_matches_closed_form(maillard_params):
    schema = wine_schema()
    y = _maillard_wine(schema, aa=0.3, s=80.0, t=298.15)  # off T_ref so the Arrhenius factor bites
    d = MaillardStrecker().derivatives(0.0, y, schema, maillard_params)
    cf = _maillard_closed_form(schema, maillard_params, y, 298.15)

    assert cf["sotolon"] > 0.0  # products are live (guards against a vacuous pass)
    for pool in _MAILLARD_C:
        assert schema.get(d, pool) == pytest.approx(cf[pool])
    assert schema.get(d, "CO2") == pytest.approx(cf["CO2"])
    # amino_acids drawn sized to the TOTAL product carbon (all six products + CO2); N deaminated.
    product_carbon = sum(cf[pool] * _MAILLARD_C[pool] for pool in _MAILLARD_C) + cf["CO2"] * _CO2_C
    aa_mass = product_carbon / _AA_C
    assert schema.get(d, "amino_acids") == pytest.approx(-aa_mass)
    assert schema.get(d, "N") == pytest.approx(aa_mass * _AA_N)
    # S is a read-only driver — NOT consumed here (its draw is booked by D-88). And NO o2 term.
    assert schema.get(d, "S") == 0.0
    assert schema.get(d, "o2") == 0.0
    for var in ("X", "E", "esters", "fusels", "Byp", "acetaldehyde"):
        assert schema.get(d, var) == 0.0


def test_maillard_carbon_closes_per_rhs(maillard_params):
    # CARBON closes to machine precision: the arginine carbon leaving amino_acids equals the carbon
    # entering all six products + CO2 (the draw is sized to match) — a pure on-ledger transfer.
    schema = wine_schema()
    d = MaillardStrecker().derivatives(0.0, _maillard_wine(schema), schema, maillard_params)
    carbon_residual = (
        sum(schema.get(d, pool) * _MAILLARD_C[pool] for pool in _MAILLARD_C)
        + schema.get(d, "CO2") * _CO2_C
        + schema.get(d, "amino_acids") * _AA_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)


def test_maillard_nitrogen_closes_per_rhs(maillard_params):
    # NITROGEN closes: all the arginine nitrogen leaving amino_acids lands in the N pool (every
    # product is nitrogen-free — the deamination, the D-45/D-75 idiom).
    schema = wine_schema()
    d = MaillardStrecker().derivatives(0.0, _maillard_wine(schema), schema, maillard_params)
    nitrogen_residual = schema.get(d, "amino_acids") * _AA_N + schema.get(d, "N") * 1.0
    assert nitrogen_residual == pytest.approx(0.0, abs=1e-18)


def test_maillard_inert_without_amino_acids(maillard_params):
    # No amino acids ⇒ no thermal Strecker: the HARD gate that IS the isolability guarantee — an
    # undosed wine is byte-for-byte the case without this Process. Also the aa<0 undershoot guard.
    schema = wine_schema()
    ps = ProcessSet(schema, [MaillardStrecker()], strict=True)
    assert np.array_equal(
        ps.total_derivatives(0.0, _maillard_wine(schema, aa=0.0), maillard_params), schema.zeros()
    )
    assert np.array_equal(
        MaillardStrecker().derivatives(
            0.0, _maillard_wine(schema, aa=-1e-6), schema, maillard_params
        ),
        schema.zeros(),
    )


def test_maillard_inert_without_sugar(maillard_params):
    # No residual sugar ⇒ no dicarbonyls ⇒ nothing: a dry wine (S = 0) makes ~none. A SOFT driver
    # (unlike the aa hard gate) — physically a dry wine still makes a trace, but the S=0 guard is a
    # clean no-op here (and absorbs the S<0 solver undershoot).
    schema = wine_schema()
    assert np.array_equal(
        MaillardStrecker().derivatives(0.0, _maillard_wine(schema, s=0.0), schema, maillard_params),
        schema.zeros(),
    )
    assert np.array_equal(
        MaillardStrecker().derivatives(
            0.0, _maillard_wine(schema, s=-1e-6), schema, maillard_params
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
        0.0, _maillard_wine(schema, o2=0.0), schema, maillard_params
    )
    oxygenated = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, o2=0.05), schema, maillard_params
    )
    assert np.array_equal(sealed, oxygenated)
    assert schema.get(sealed, "sotolon") > 0.0  # non-vacuous: it IS producing while sealed


def test_maillard_first_order_in_sugar(maillard_params):
    # The residual sugar is the dicarbonyl driver: doubling S doubles the production (first-order),
    # holding the aa gate fixed. The bounded-vs-unbounded concern is on amino_acids (the limiting
    # reagent), not sugar — sugar drives the RATE.
    schema = wine_schema()
    base = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, s=40.0), schema, maillard_params
    )
    dbl = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, s=80.0), schema, maillard_params
    )
    assert schema.get(dbl, "sotolon") == pytest.approx(2.0 * schema.get(base, "sotolon"))


def test_maillard_availability_gate_saturates(maillard_params):
    # The aa availability gate aa/(K+aa) saturates: production per unit sugar rises then plateaus as
    # amino_acids climbs (the same smooth-Monod shape D-75 uses). Below/at/well-above K.
    schema = wine_schema()
    k = maillard_params["K_amino_acids"]
    # Equal-RATIO aa steps (k → 10k → 100k), so a saturating gate shows strictly diminishing gains.
    low = MaillardStrecker().derivatives(0.0, _maillard_wine(schema, aa=k), schema, maillard_params)
    mid = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, aa=10.0 * k), schema, maillard_params
    )
    high = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, aa=100.0 * k), schema, maillard_params
    )
    lo, md, hi = (schema.get(x, "sotolon") for x in (low, mid, high))
    assert lo < md < hi  # monotone increasing in aa
    assert (md - lo) > (hi - md)  # but saturating (diminishing returns) — the gate flattens


def test_maillard_split_normalizes_and_sotolon_has_no_co2(maillard_params):
    # The six composition weights NORMALIZE to fractions summing to 1 (the split-hygiene the advisor
    # flagged), and — the load-bearing flag — sotolon (a furanone, NOT a decarboxylation product)
    # contributes NO CO2, while the CO2 exactly matches the FIVE Strecker aldehydes' mole sum.
    schema = wine_schema()
    y = _maillard_wine(schema)
    d = MaillardStrecker().derivatives(0.0, y, schema, maillard_params)
    # mole rate of each product = mass rate / molar mass
    masses = {
        "methional": M_METHIONAL,
        "phenylacetaldehyde": M_PHENYLACETALDEHYDE,
        "2_methylbutanal": M_2_METHYLBUTANAL,
        "3_methylbutanal": M_3_METHYLBUTANAL,
        "2_methylpropanal": M_2_METHYLPROPANAL,
        "sotolon": M_SOTOLON,
    }
    n = {pool: schema.get(d, pool) / masses[pool] for pool in masses}
    n_total = sum(n.values())
    # The normalized split sums to 1 (fractions = n_i / n_total).
    assert sum(n[pool] / n_total for pool in n) == pytest.approx(1.0)
    # CO2 == the FIVE decarboxylating aldehydes' mole sum (sotolon EXCLUDED) — the flag is load-
    # bearing (conservation closes for ANY CO2 attribution, so only this asserts it is keyed right).
    decarb_mol = sum(n[pool] for pool in n if pool != "sotolon")
    assert schema.get(d, "CO2") / M_CO2 == pytest.approx(decarb_mol)
    # Concretely: NOT the total (sotolon really is excluded — a non-vacuous check).
    assert decarb_mol < n_total


def test_maillard_rises_with_temperature(maillard_params):
    # Warmer ages faster — the sourced direction (thermal Strecker is strongly temperature-driven).
    schema = wine_schema()
    cold = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, t=283.15), schema, maillard_params
    )
    warm = MaillardStrecker().derivatives(
        0.0, _maillard_wine(schema, t=303.15), schema, maillard_params
    )
    assert schema.get(warm, "sotolon") > schema.get(cold, "sotolon") > 0.0


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
    y0 = _maillard_wine(schema, aa=0.4, s=80.0, t=301.15)  # 28 °C, sealed (o2 defaults to 0)
    # The oxidative route on the SAME sealed state: identically zero (its o2 <= 0 guard short-
    # circuits before it even reads a param, so maillard_params suffices — nothing is computed).
    d_oxid = StreckerDegradation().derivatives(0.0, y0, schema, maillard_params)
    assert np.array_equal(d_oxid, schema.zeros())

    ps = ProcessSet(schema, [MaillardStrecker()], strict=True)
    traj = simulate(ps, params=maillard_params, y0=y0, t_span=(0.0, 24.0 * 365.0))
    assert traj.success, traj.message
    for pool in ("methional", "phenylacetaldehyde", "sotolon", "3_methylbutanal"):
        series = traj.series(pool)
        assert series[-1] > series[0] == 0.0  # produced-only, monotone accumulation from 0
    assert_nonnegative(traj, ("amino_acids", "sotolon", "methional", "N"), atol=1e-9)
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
# Caramelization (decision D-88) — the WINE-ONLY, NON-oxidative THERMAL browning: the O₂-INDEPENDENT
# thermal mirror of PhenolicBrowning (D-74). Residual SUGAR browns to melanoidin by HEAT (no O₂),
# raising the SAME A420 index D-74 accumulates — so a sealed sweet wine still darkens. The FIRST
# aging Process to consume core S: the sugar carbon lands in the on-ledger melanoidin carbon-park
# (the debris/glucan precedent), so total_carbon closes exactly (release at the sugar fraction,
# redeposit at melanoidin's). SUGAR-ONLY (nitrogen-free — caramelization, not Maillard). These tests
# pin the closed form, carbon closure per-RHS, the sugar SOFT gate (inert at S ≈ 0 / undershoot),
# the O₂-independence (no o2 term at all), the first-order-in-sugar linearity, the monotone A420
# rise, the warmer-faster ordering, the wine-only no-op on beer, the integrated sweet browning +
# closure, and the speculative tier floor.

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
    y = _aged_wine(schema, esters=0.0, t=t)
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
    for var in ("o2", "amino_acids", "E", "esters", "acetaldehyde", "sotolon", "N"):
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


def test_caramelization_is_wine_only_noop_on_beer(caramel_params):
    # Wine-only v1 (the melanoidin carbon-park is a wine slot): a hard no-op on beer even with
    # residual wort sugar (the "melanoidin" not in schema guard).
    beer = beer_schema()
    yb = beer.zeros()
    yb[beer.slice("S")] = 60.0
    yb[beer.slice("T")] = 303.15
    assert np.array_equal(Caramelization().derivatives(0.0, yb, beer, caramel_params), beer.zeros())


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

_OAK_COMPOUNDS = ("whiskey_lactone", "vanillin", "guaiacol", "eugenol")


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
    y = _aged_wine(schema, esters=0.0, t=t)
    for compound, ceiling in (ceilings or {}).items():
        y[schema.slice(f"{compound}_ceiling")] = ceiling
    for name, val in kw.items():
        y[schema.slice(name)] = val
    return y


def test_oak_metadata():
    p = OakExtraction()
    assert p.name == "oak_extraction"
    assert p.tier is Tier.SPECULATIVE
    # Writes ONLY the five extracted-compound slots — the four aroma extractives (D-77) plus the
    # ellagitannin taste extractive (D-78). The ceilings are read, never written (a set-and-hold
    # constant the add_oak verb owns). Off every ledger, so nothing conserved moves.
    assert set(p.touches) == {"whiskey_lactone", "vanillin", "guaiacol", "eugenol", "ellagitannin"}
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
    # esters_eq > 0), so a solver undershoot conc = −ε must NOT flip max(0, ceiling − conc) into
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


# =====================================================================================
# EllagitanninOxidation (decision D-78) — the WINE-ONLY oak-tannin O₂-scavenging sink, the BRIDGE
# from the D-77 oak extractive axis to the O₂ sub-axis. Oak's hydrolysable tannin (the ellagitannin
# pool OakExtraction fills) is a sacrificial antioxidant: dissolved O₂ oxidises it (bilinear
# [o2]·[ellagitannin], the SulfiteOxidation form), CONSUMING it at a mass-based yield y_ellag_per_o2
# (g ellag / g O₂ — no fake molar mass for the lumped macromolecule). The EMERGENT SPINE is oak
# PROTECTION: an oaked + oxygenated wine browns LESS (lower A420) and makes LESS oxidative
# acetaldehyde than an un-oaked wine at the same O₂ dose (the D-72 "SO₂ protects" threshold with a
# RENEWABLE buffer). Substrate-gated on ellagitannin ⇒ adds on top with NO re-baseline of the
# k_ethanol_oxidation + k_browning = 5.0e-4 anchor. Off every ledger (both slots unweighted), so it
# moves nothing conserved. These tests pin the closed form, the bilinearity, the reaction-scale
# temperature ordering, the doubly-substrate-gated inertness (KeyError-safe without oak.yaml), the
# wine-only no-op on beer, THE PROTECTION SPINE (partial, not total), the sacrificial-consumption
# softening (astringency_series), the off-every-ledger invariance, and the speculative tier floor.


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
    y = _aged_wine(schema, esters=0.0, t=t, o2=o2, ellagitannin=ellag)
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
        0.0, _aged_wine(schema, esters=0.0, o2=0.02, ellagitannin=0.05), schema, ellag_params
    )
    twice_o2 = EllagitanninOxidation().derivatives(
        0.0, _aged_wine(schema, esters=0.0, o2=0.04, ellagitannin=0.05), schema, ellag_params
    )
    twice_ellag = EllagitanninOxidation().derivatives(
        0.0, _aged_wine(schema, esters=0.0, o2=0.02, ellagitannin=0.10), schema, ellag_params
    )
    assert schema.get(twice_o2, "o2") == pytest.approx(2.0 * schema.get(base, "o2"))
    assert schema.get(twice_ellag, "o2") == pytest.approx(2.0 * schema.get(base, "o2"))
    assert schema.get(base, "o2") < 0.0  # actually scavenging


def test_ellagitannin_oxidation_inert_without_o2_or_tannin(ellag_params):
    # Doubly substrate-gated: no O₂ OR no ellagitannin ⇒ byte-for-byte zero. A reductive (no
    # add_oxygen) or an un-oaked aging is exactly the case without this Process (isolability #3).
    schema = wine_schema()
    p = EllagitanninOxidation()
    no_o2 = _aged_wine(schema, esters=0.0, o2=0.0, ellagitannin=0.08)
    no_tannin = _aged_wine(schema, esters=0.0, o2=0.03, ellagitannin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_o2, schema, ellag_params), schema.zeros())
    assert np.array_equal(p.derivatives(0.0, no_tannin, schema, ellag_params), schema.zeros())
    # <= 0 also absorbs solver undershoot (a spurious −ε in either driver ⇒ no draw).
    undershoot = _aged_wine(schema, esters=0.0, o2=-1e-9, ellagitannin=0.08)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, ellag_params), schema.zeros())


def test_ellagitannin_oxidation_gate_before_params_is_keyerror_safe(params):
    # Gate on the ellagitannin STATE before reading any oak param, so an enabled-but-undosed Process
    # never KeyErrors when oak.yaml is ABSENT (the ``params`` fixture is wine+aging only, no
    # k_ellagitannin_oxidation). An un-oaked wine (ellag=0) returns zero without touching oak
    # params.
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, o2=0.03, ellagitannin=0.0)
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
        _aged_wine(schema, esters=0.0, t=283.15, o2=0.03, ellagitannin=0.08),
        schema,
        ellag_params,
    )
    warm = p.derivatives(
        0.0,
        _aged_wine(schema, esters=0.0, t=303.15, o2=0.03, ellagitannin=0.08),
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
    oaked0 = _aged_wine(schema, esters=0.0, o2=o2_dose, ellagitannin=0.1, ellagitannin_ceiling=0.1)
    oaked = simulate(ps, params=ellag_params, y0=oaked0, t_span=span)
    # Un-oaked: identical but no tannin (and no ceiling), same O₂ dose.
    unoaked0 = _aged_wine(
        schema, esters=0.0, o2=o2_dose, ellagitannin=0.0, ellagitannin_ceiling=0.0
    )
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
    y0 = _aged_wine(schema, esters=0.0, o2=0.05, ellagitannin=0.1)
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
    y0 = _aged_wine(schema, esters=0.0, o2=0.05, ellagitannin=0.1)
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
    y = _aged_wine(schema, esters=0.0, t=t, anthocyanin=antho, tannin=tannin)
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
        0.0, _aged_wine(schema, esters=0.0, anthocyanin=0.2, tannin=1.5), schema, poly_params
    )
    twice_antho = TanninAnthocyaninCondensation().derivatives(
        0.0, _aged_wine(schema, esters=0.0, anthocyanin=0.4, tannin=1.5), schema, poly_params
    )
    twice_tannin = TanninAnthocyaninCondensation().derivatives(
        0.0, _aged_wine(schema, esters=0.0, anthocyanin=0.2, tannin=3.0), schema, poly_params
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
    no_antho = _aged_wine(schema, esters=0.0, anthocyanin=0.0, tannin=2.0)
    no_tannin = _aged_wine(schema, esters=0.0, anthocyanin=0.3, tannin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_antho, schema, poly_params), schema.zeros())
    assert np.array_equal(p.derivatives(0.0, no_tannin, schema, poly_params), schema.zeros())
    # Solver undershoot (negative) is likewise absorbed.
    undershoot = _aged_wine(schema, esters=0.0, anthocyanin=-1e-9, tannin=2.0)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_polymerization_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks k_polymerization,
    # yet a white wine (anthocyanin 0) returns zero — the gate-on-STATE-before-params discipline.
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, anthocyanin=0.0, tannin=0.0)
    d = TanninAnthocyaninCondensation().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_polymerization_rises_with_temperature(poly_params):
    # Warmer condenses faster (E_a > 0, reaction-scale): the |anthocyanin| draw grows with T.
    schema = wine_schema()
    cold = TanninAnthocyaninCondensation().derivatives(
        0.0,
        _aged_wine(schema, esters=0.0, t=283.15, anthocyanin=0.3, tannin=2.0),
        schema,
        poly_params,
    )
    warm = TanninAnthocyaninCondensation().derivatives(
        0.0,
        _aged_wine(schema, esters=0.0, t=303.15, anthocyanin=0.3, tannin=2.0),
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
    y0 = _aged_wine(schema, esters=0.0, anthocyanin=0.3, tannin=2.0)
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
    y0 = _aged_wine(schema, esters=0.0, anthocyanin=antho0, tannin=tannin0)
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
    y0 = _aged_wine(schema, esters=0.0, anthocyanin=antho0, tannin=2.0)
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
    y = _aged_wine(schema, esters=0.0, t=t, acetaldehyde=acet, anthocyanin=antho, tannin=tannin)
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
        y = _aged_wine(schema, esters=0.0, acetaldehyde=acet, anthocyanin=antho, tannin=tannin)
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
        y = _aged_wine(schema, esters=0.0, **kw)
        assert np.array_equal(p.derivatives(0.0, y, schema, poly_params), schema.zeros())


def test_bridge_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks
    # k_acetaldehyde_bridge, yet a white wine (anthocyanin 0) returns zero — the
    # gate-on-STATE-before-params discipline.
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, acetaldehyde=0.05, anthocyanin=0.0, tannin=0.0)
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
        schema, esters=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.0, **acids
    )
    sulfited = _aged_wine(
        schema, esters=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.05, **acids
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
        schema, esters=0.0, t=t, acetaldehyde=acet, anthocyanin=antho, tannin=tannin, so2_total=0.0
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
        _aged_wine(schema, esters=0.0, t=283.15, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0),
        schema,
        poly_params,
    )
    warm = p.derivatives(
        0.0,
        _aged_wine(schema, esters=0.0, t=303.15, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0),
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
        esters=0.0,
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
    y = _aged_wine(schema, esters=0.0, t=t, o2=o2, anthocyanin=antho)
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
        y = _aged_wine(schema, esters=0.0, o2=o2, anthocyanin=antho)
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
    no_o2 = _aged_wine(schema, esters=0.0, o2=0.0, anthocyanin=0.3)
    no_antho = _aged_wine(schema, esters=0.0, o2=0.03, anthocyanin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_o2, schema, poly_params), schema.zeros())
    assert np.array_equal(p.derivatives(0.0, no_antho, schema, poly_params), schema.zeros())
    undershoot = _aged_wine(schema, esters=0.0, o2=-1e-9, anthocyanin=0.3)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_fading_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks the fade rate,
    # yet a white wine (anthocyanin 0, even with O₂ dosed) returns zero — gate-before-params.
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, o2=0.03, anthocyanin=0.0)
    d = AnthocyaninFading().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_fading_rises_with_temperature(poly_params):
    # Warmer fades faster (E_a > 0, reaction-scale): the |anthocyanin| draw grows with T.
    schema = wine_schema()
    p = AnthocyaninFading()
    cold = p.derivatives(
        0.0, _aged_wine(schema, esters=0.0, t=283.15, o2=0.03, anthocyanin=0.3), schema, poly_params
    )
    warm = p.derivatives(
        0.0, _aged_wine(schema, esters=0.0, t=303.15, o2=0.03, anthocyanin=0.3), schema, poly_params
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
    y0 = _aged_wine(schema, esters=0.0, o2=0.04, anthocyanin=0.3)
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
    y0 = _aged_wine(schema, esters=0.0, o2=0.05, anthocyanin=antho0, tannin=2.0)
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
    y = _aged_wine(schema, esters=0.0, t=t, o2=0.03, anthocyanin=antho)
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
        y = _aged_wine(schema, esters=0.0, o2=o2, anthocyanin=antho)
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
    no_antho = _aged_wine(schema, esters=0.0, o2=0.03, anthocyanin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_antho, schema, poly_params), schema.zeros())
    # Undershoot absorbed too.
    undershoot = _aged_wine(schema, esters=0.0, o2=0.03, anthocyanin=-1e-9)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_thermal_fade_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks the thermal-fade
    # rate, yet a white wine (anthocyanin 0) returns zero — gate-before-params.
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, o2=0.03, anthocyanin=0.0)
    d = ThermalAnthocyaninFade().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_thermal_fade_rises_with_temperature(poly_params):
    # Warmer fades faster (E_a > 0, reaction-scale): the |anthocyanin| draw grows with T — the
    # 'warm storage kills colour even anaerobically' temperature lever, the only lever this route
    # has (no o2/SO₂ coupling).
    schema = wine_schema()
    p = ThermalAnthocyaninFade()
    cold = p.derivatives(
        0.0, _aged_wine(schema, esters=0.0, t=283.15, o2=0.0, anthocyanin=0.3), schema, poly_params
    )
    warm = p.derivatives(
        0.0, _aged_wine(schema, esters=0.0, t=303.15, o2=0.0, anthocyanin=0.3), schema, poly_params
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
    y0 = _aged_wine(schema, esters=0.0, o2=0.0, anthocyanin=0.3)  # o2=0: fades anyway
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
    y0 = _aged_wine(schema, esters=0.0, o2=0.0, anthocyanin=antho0)  # anaerobic, sealed red

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
    no_so2 = _aged_wine(schema, esters=0.0, o2=0.0, anthocyanin=0.3)
    with_so2 = _aged_wine(schema, esters=0.0, o2=0.0, anthocyanin=0.3, so2_total=0.05)
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
    y = _aged_wine(schema, esters=0.0, t=t, anthocyanin=0.3, tannin=tannin)
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
        y = _aged_wine(schema, esters=0.0, anthocyanin=0.3, tannin=tannin)
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
        y = _aged_wine(schema, esters=0.0, anthocyanin=anthocyanin, o2=o2, tannin=2.0)
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
    no_tannin = _aged_wine(schema, esters=0.0, anthocyanin=0.3, tannin=0.0)
    assert np.array_equal(p.derivatives(0.0, no_tannin, schema, poly_params), schema.zeros())
    undershoot = _aged_wine(schema, esters=0.0, anthocyanin=0.3, tannin=-1e-9)
    assert np.array_equal(p.derivatives(0.0, undershoot, schema, poly_params), schema.zeros())


def test_tannin_self_poly_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture (wine_generic + aging.yaml, NO polymerization.yaml) lacks the rate, yet a
    # no-tannin wine returns zero — gate-before-params.
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, anthocyanin=0.3, tannin=0.0)
    d = TanninSelfPolymerization().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_tannin_self_poly_rises_with_temperature(poly_params):
    # Warmer polymerizes (softens) faster (E_a > 0, reaction-scale): the |tannin| draw grows with T.
    schema = wine_schema()
    p = TanninSelfPolymerization()
    cold = p.derivatives(
        0.0, _aged_wine(schema, esters=0.0, t=283.15, tannin=2.0), schema, poly_params
    )
    warm = p.derivatives(
        0.0, _aged_wine(schema, esters=0.0, t=303.15, tannin=2.0), schema, poly_params
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
    y0 = _aged_wine(schema, esters=0.0, tannin=3.0)  # no anthocyanin — softens anyway
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
    y0 = _aged_wine(schema, esters=0.0, anthocyanin=0.0, tannin=tannin0)  # white / no anthocyanin

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
    y = _aged_wine(schema, esters=0.0, t=t, acetaldehyde=acet, anthocyanin=0.3, tannin=tannin)
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
        y = _aged_wine(schema, esters=0.0, acetaldehyde=acet, anthocyanin=0.3, tannin=tannin)
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
        y = _aged_wine(schema, esters=0.0, anthocyanin=0.3, **kw)
        assert np.array_equal(p.derivatives(0.0, y, schema, poly_params), schema.zeros())


def test_tannin_ethyl_gate_before_params_is_keyerror_safe(params):
    # An enabled-but-undosed Process must not KeyError when polymerization.yaml is absent: the
    # ``params`` fixture lacks k_tannin_ethyl_tannin, yet a no-tannin wine returns zero.
    schema = wine_schema()
    y = _aged_wine(schema, esters=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=0.0)
    d = TanninEthylTanninCondensation().derivatives(0.0, y, schema, params)
    assert np.array_equal(d, schema.zeros())


def test_tannin_ethyl_reads_free_acetaldehyde_under_so2(bridge_params):
    # SO₂-bound acetaldehyde can't bridge (the D-47/D-80 precedent): at the SAME total acetaldehyde
    # a sulfited wine bridges SLOWER than an unsulfited one — SO₂ DELAYS the tannin softening.
    schema = wine_schema()
    p = TanninEthylTanninCondensation()
    acids = {"tartaric": 4.0, "cation_charge": 0.012}  # a real acid state so pH solves
    unsulfited = _aged_wine(
        schema, esters=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.0, **acids
    )
    sulfited = _aged_wine(
        schema, esters=0.0, acetaldehyde=0.05, anthocyanin=0.3, tannin=2.0, so2_total=0.05, **acids
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
        schema, esters=0.0, t=t, acetaldehyde=acet, anthocyanin=0.3, tannin=tannin, so2_total=0.0
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
        _aged_wine(schema, esters=0.0, t=283.15, acetaldehyde=0.05, tannin=2.0),
        schema,
        poly_params,
    )
    warm = p.derivatives(
        0.0,
        _aged_wine(schema, esters=0.0, t=303.15, acetaldehyde=0.05, tannin=2.0),
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
        esters=0.0,
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
