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

from fermentation.core.acidbase import bisulfite_so2_at_ph, ph_of_state
from fermentation.core.chemistry import (
    M_ACETALDEHYDE,
    M_CO2,
    M_ETHANOL,
    M_METHIONAL,
    M_O2,
    M_PHENYLACETALDEHYDE,
    M_SO2,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics import (
    EsterHydrolysis,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    StreckerDegradation,
    SulfiteOxidation,
    arrhenius_factor,
)
from fermentation.core.kinetics.aging import _SO2_PER_O2
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
#: Carbon fractions of the Strecker pools + the amino-acid source (D-75), for the closure checks.
_METHIONAL_C = carbon_mass_fraction("methional")
_PHENYLACET_C = carbon_mass_fraction("phenylacetaldehyde")
_CO2_C = carbon_mass_fraction("CO2")
_AA_C = carbon_mass_fraction(AMINO_ACID_SPECIES)
_AA_N = nitrogen_mass_fraction(AMINO_ACID_SPECIES)


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


def test_strecker_split_methional_dominant(params):
    # The mol split between the two aldehydes is f_methional : (1 - f_methional); with the default
    # f_methional = 0.6 methional dominates. Verified as a MOLAR ratio (independent of the two
    # differing molar masses).
    schema = wine_schema()
    d = StreckerDegradation().derivatives(0.0, _strecker_wine(schema), schema, params)
    meth_mol = schema.get(d, "methional") / M_METHIONAL
    phenyl_mol = schema.get(d, "phenylacetaldehyde") / M_PHENYLACETALDEHYDE
    f_meth = params["f_methional"]
    assert meth_mol / phenyl_mol == pytest.approx(f_meth / (1.0 - f_meth))
    assert meth_mol > phenyl_mol > 0.0  # methional-dominant (the cooked-potato oxidative marker)
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
