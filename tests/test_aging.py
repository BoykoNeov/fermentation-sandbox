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

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics import EsterHydrolysis, arrhenius_factor
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_mass,
)

#: Carbon fractions of the three pools the transfer touches (mirror the Process constants).
_ESTER_C = carbon_mass_fraction("ethyl_acetate")
_FUSEL_C = carbon_mass_fraction("isoamyl_alcohol")
_BYP_C = carbon_mass_fraction("succinic_acid")
#: The 5:2 split, from the isoamyl-acetate stand-in reaction (isoamyl alcohol 5 C : acetic 2 C).
_FUSEL_SHARE = 5.0 / 7.0
_BYP_SHARE = 2.0 / 7.0


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
