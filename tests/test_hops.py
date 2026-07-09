"""Tests for hop bittering — boil isomerization + fermentation loss (decision D-64).

Bitterness comes from *iso*-alpha-acids made in the boil by thermal isomerization of the
hop's alpha-acids (Malowicki & Shellhammer 2005, sourced kinetics), then partly lost during
fermentation by adsorption onto viable yeast. The engine treats the boil as a wort-side
compile-seam calc (closed-form consecutive first-order reaction, wired into ``iso_alpha`` at
t=0 like ``initial_ph``); only :class:`IsoAlphaAcidLoss` is a dynamic Process.

Coverage:

* **Boil kinetics** — the Malowicki rate constants at 100 C; the closed-form intermediate
  matches a numerically-integrated A->B->C reaction; ~48% of alpha isomerizes at 60 min/100 C
  on the RISING limb (peak ~3 h); monotone in alpha dose; hotter boil isomerizes faster; a
  zero-length boil gives zero.
* **The loss Process** — iso-alpha declines during fermentation, gated on viable yeast (no
  yeast => frozen bitterness), inert when unhopped, first-order, touches only ``iso_alpha``.
* **Isolability / off-ledger** — a hopped beer's ``total_carbon`` is BYTE-FOR-BYTE the unhopped
  beer's (iso-alpha is exogenous, off the conservation ledger); unhopped runs disable the loss
  Process and keep ``iso_alpha`` VALIDATED; hops on a non-beer medium and hops without a volume
  are loud errors.
* **Tinseth cross-check** (fit-vs-fit, NOT validation, per §3.5) — canonical recipes land in
  the ballpark of the Tinseth utilization model.
* **Tier** — the finished IBU is capped speculative by the loss/efficiency inputs even though
  the boil kinetics are sourced/plausible (D-1 propagation).
"""

import math
from collections.abc import Sequence

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from fermentation.analysis import ibu_series
from fermentation.core.kinetics import IsoAlphaAcidLoss, boil_rate_constants, iso_alpha_fraction
from fermentation.core.media import beer_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import Trajectory, simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario.schema import HopAddition
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon

BOIL_K = 373.15  # 100 C reference boil


def _simulate(sc: Scenario) -> Trajectory:
    """Compile and integrate (no events on these single-knot scenarios, so ``simulate`` ==
    ``run`` here, and it returns the plain ``Trajectory`` the series readouts type against)."""
    cs = compile_scenario(sc, strict=True)
    return simulate(cs.process_set, cs.param_values, cs.y0, cs.t_span_h)


def _ibu(sc: Scenario) -> FloatArray:
    return ibu_series(_simulate(sc))


@pytest.fixture
def params():
    return load_parameters(default_data_dir() / "hops.yaml").resolve()


# -- boil isomerization kinetics (pure closed form) ---------------------------


def test_rate_constants_match_malowicki_at_100c(params):
    """k1, k2 at 100 C reproduce Malowicki & Shellhammer 2005 (~0.0125, ~0.0031 /min)."""
    k1, k2 = boil_rate_constants(BOIL_K, params)
    assert k1 == pytest.approx(0.0125, abs=2e-4)
    assert k2 == pytest.approx(0.0031, abs=2e-4)
    # k2 < k1 at boil temperature ⇒ iso-alpha accumulates on the rising limb of a normal boil.
    assert k2 < k1


def test_activation_energies_are_the_reported_values(params):
    """Ea/R equals the reported Arrhenius exponents (98.6 / 108.0 kJ/mol), internal-consistency."""
    R = 8.314462618
    assert params["Ea_iso"] / R == pytest.approx(11858.0, rel=1e-3)
    assert params["Ea_iso_degradation"] / R == pytest.approx(12994.0, rel=1e-3)
    # degradation is more temperature-sensitive than isomerization (Malowicki's ordering)
    assert params["Ea_iso_degradation"] > params["Ea_iso"]


def test_closed_form_matches_numerical_consecutive_reaction(params):
    """The analytic intermediate equals a numerically-integrated A -k1-> B -k2-> C."""
    k1, k2 = boil_rate_constants(BOIL_K, params)

    def rhs(_t, y):
        a, b = y
        return [-k1 * a, k1 * a - k2 * b]

    for minutes in (15.0, 45.0, 90.0, 180.0):
        sol = solve_ivp(rhs, (0.0, minutes), [1.0, 0.0], rtol=1e-10, atol=1e-12)
        numeric_b = float(sol.y[1, -1])
        assert iso_alpha_fraction(minutes, BOIL_K, params) == pytest.approx(numeric_b, rel=1e-6)


def test_about_half_isomerizes_at_60min_on_the_rising_limb(params):
    """~48% of alpha is iso-alpha at 60 min/100 C, and it is still rising (peak ~3 h)."""
    f60 = iso_alpha_fraction(60.0, BOIL_K, params)
    assert 0.40 <= f60 <= 0.55  # the advisor's sourced sanity band
    # Rising limb: monotone increasing through a normal boil, peaks well past 90 min.
    fr1 = [iso_alpha_fraction(t, BOIL_K, params) for t in (15, 30, 60, 90)]
    assert fr1 == sorted(fr1)
    assert iso_alpha_fraction(90.0, BOIL_K, params) > f60
    # A very long boil eventually turns over (degradation wins) — the peak is finite.
    assert iso_alpha_fraction(240.0, BOIL_K, params) < iso_alpha_fraction(180.0, BOIL_K, params)


def test_zero_boil_gives_no_isomerization(params):
    assert iso_alpha_fraction(0.0, BOIL_K, params) == 0.0


def test_hotter_boil_isomerizes_faster(params):
    """A hotter boil converts more alpha at a fixed short time (Arrhenius)."""
    hot = iso_alpha_fraction(20.0, 373.15 + 10.0, params)
    cool = iso_alpha_fraction(20.0, 373.15 - 20.0, params)  # e.g. a whirlpool
    assert hot > cool


# -- the fermentation loss Process --------------------------------------------


def _beer_state(
    schema: StateSchema, *, x: float = 2.0, iso: float = 0.03, t_k: float = 293.15
) -> FloatArray:
    y = schema.zeros()
    y[schema.slice("X")] = x
    y[schema.slice("iso_alpha")] = iso
    y[schema.slice("T")] = t_k
    return y


def test_loss_is_first_order_in_iso_and_gated_on_viable_yeast(params):
    schema = beer_schema()
    proc = IsoAlphaAcidLoss()
    # First-order in iso-alpha: doubling iso doubles the loss rate.
    d1 = proc.derivatives(0.0, _beer_state(schema, iso=0.02), schema, params)
    d2 = proc.derivatives(0.0, _beer_state(schema, iso=0.04), schema, params)
    r1 = float(d1[schema.slice("iso_alpha")][0])
    r2 = float(d2[schema.slice("iso_alpha")][0])
    assert r1 < 0.0  # it is a loss
    assert r2 == pytest.approx(2.0 * r1, rel=1e-12)
    # First-order in viable biomass: no yeast ⇒ no loss (frozen bitterness).
    d0 = proc.derivatives(0.0, _beer_state(schema, x=0.0), schema, params)
    assert float(d0[schema.slice("iso_alpha")][0]) == 0.0


def test_loss_is_inert_when_unhopped(params):
    """iso_alpha = 0 ⇒ zero contribution (an unhopped beer is untouched)."""
    schema = beer_schema()
    d = IsoAlphaAcidLoss().derivatives(0.0, _beer_state(schema, iso=0.0), schema, params)
    assert not np.any(d)


def test_loss_touches_only_iso_alpha(params):
    """Off the carbon ledger: the Process contributes to no column but iso_alpha (strict)."""
    schema = beer_schema()
    ps = ProcessSet(schema, [IsoAlphaAcidLoss()], strict=True)
    d = ps.total_derivatives(0.0, _beer_state(schema), params)
    for name in schema.names:
        if name != "iso_alpha":
            assert float(np.sum(np.abs(d[schema.slice(name)]))) == 0.0
    assert float(d[schema.slice("iso_alpha")][0]) < 0.0


# -- compile + integrate ------------------------------------------------------


def _beer(
    hops: Sequence[HopAddition] = (),
    volume: float | None = None,
    boil_c: float = 100.0,
    days: float = 10.0,
) -> Scenario:
    return Scenario(
        name="test-beer",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 60.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 1.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        hops=list(hops),
        batch_volume_liters=volume,
        boil_celsius=boil_c,
        duration_days=days,
    )


def test_hopped_beer_seeds_iso_alpha_and_it_declines(params):
    """The boil iso-alpha is wired at t=0 and falls during fermentation to a positive finish."""
    sc = _beer(
        hops=[HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=60)], volume=18.93
    )
    ibu = _ibu(sc)
    assert ibu[0] == pytest.approx(19.5, abs=1.0)  # end-of-boil, matches the hand calc
    assert 0.0 < ibu[-1] < ibu[0]  # declines but does not vanish
    loss_frac = (ibu[0] - ibu[-1]) / ibu[0]
    assert 0.05 <= loss_frac <= 0.20  # the ~5-20% wort-to-beer bitterness drop


def test_more_hops_and_longer_boil_give_more_bitterness(params):
    def finish(hops: Sequence[HopAddition]) -> float:
        return float(_ibu(_beer(hops=hops, volume=18.93))[-1])

    base = [HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=60)]
    more_mass = [HopAddition(alpha_acid_percent=5.0, grams=56.70, boil_minutes=60)]
    longer = [HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=90)]
    assert finish(more_mass) > finish(base)
    assert finish(longer) > finish(base)


def test_multiple_additions_sum(params):
    one = [HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=60)]
    two = [
        HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=60),
        HopAddition(alpha_acid_percent=8.0, grams=14.0, boil_minutes=15),
    ]
    f1 = _ibu(_beer(hops=one, volume=18.93))[0]
    f2 = _ibu(_beer(hops=two, volume=18.93))[0]
    assert f2 > f1


def test_unhopped_beer_has_zero_bitterness_and_disables_the_loss(params):
    cs = compile_scenario(_beer(), strict=True)
    assert not cs.process_set.is_enabled(IsoAlphaAcidLoss.name)
    ibu = ibu_series(simulate(cs.process_set, cs.param_values, cs.y0, cs.t_span_h))
    assert np.all(ibu == 0.0)
    # With the loss Process disabled, the empty iso_alpha slot keeps its VALIDATED tier.
    tier = cs.process_set.tier_of("iso_alpha", cs.parameters.tier_map())
    assert tier is Tier.VALIDATED


def test_hops_on_wine_is_an_error():
    sc = Scenario(
        name="bad",
        medium="wine",
        initial={"brix": 22.0, "yan_mgl": 200.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        hops=[HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=60)],
        batch_volume_liters=20.0,
    )
    with pytest.raises(ValueError, match="beer-only"):
        compile_scenario(sc)


def test_hops_without_volume_is_rejected_by_the_schema():
    with pytest.raises(ValueError, match="batch_volume_liters"):
        _beer(hops=[HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=60)], volume=None)


# -- carbon invariant: iso-alpha is off the ledger ----------------------------


def test_hopping_leaves_total_carbon_byte_for_byte(params):
    """A hopped beer's carbon trajectory is IDENTICAL to the unhopped beer's (off-ledger)."""
    unhopped = compile_scenario(_beer(), strict=True)
    hopped = compile_scenario(
        _beer(
            hops=[HopAddition(alpha_acid_percent=6.0, grams=40.0, boil_minutes=75)], volume=18.93
        ),
        strict=True,
    )
    f_c = unhopped.param_values["biomass_C_fraction"]
    tu = unhopped.run()
    th = hopped.run()
    # Both conserve carbon ...
    assert_conserved(th, total_carbon(th.schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_nonnegative(th, ("iso_alpha",))
    # ... and the carbon series are identical to machine precision (hops perturb nothing on the
    # ledger; iso_alpha has weight 0 in total_carbon and the loss Process touches nothing else).
    carbon = total_carbon(th.schema, biomass_carbon_fraction=f_c)
    cu = np.array([carbon(tu.y[:, i]) for i in range(tu.y.shape[1])])
    ch = np.array([carbon(th.y[:, i]) for i in range(th.y.shape[1])])
    assert np.allclose(cu, ch, rtol=0, atol=1e-9)


# -- Tinseth cross-check (fit-vs-fit, not validation — §3.5) ------------------


def _tinseth_ibu(
    alpha_pct: float, grams: float, volume_l: float, boil_min: float, sg: float
) -> float:
    """Reference finished IBU from the Tinseth utilization model (an independent empirical fit)."""
    bigness = 1.65 * 0.000125 ** (sg - 1.0)
    boil_factor = (1.0 - math.exp(-0.04 * boil_min)) / 4.15
    util = bigness * boil_factor
    conc = (alpha_pct / 100.0) * grams * 1000.0 / volume_l  # mg/L if fully utilized
    return float(conc * util)


@pytest.mark.parametrize(
    "alpha,grams,boil",
    [(5.0, 28.35, 60.0), (8.0, 28.35, 60.0), (5.0, 56.7, 30.0)],
)
def test_finished_ibu_is_in_the_tinseth_ballpark(params, alpha, grams, boil):
    """Finished IBU tracks Tinseth within ~30% across recipes — fit-vs-fit, not validation."""
    sc = _beer(
        hops=[HopAddition(alpha_acid_percent=alpha, grams=grams, boil_minutes=boil)], volume=18.93
    )
    model_ibu = float(_ibu(sc)[-1])
    tinseth = _tinseth_ibu(alpha, grams, 18.93, boil, sg=1.050)
    assert model_ibu == pytest.approx(tinseth, rel=0.30)


# -- tier propagation ---------------------------------------------------------


def test_finished_ibu_tier_is_capped_speculative(params):
    """Sourced boil kinetics are plausible, but the loss/efficiency inputs cap finished IBU."""
    cs = compile_scenario(
        _beer(
            hops=[HopAddition(alpha_acid_percent=5.0, grams=28.35, boil_minutes=60)], volume=18.93
        ),
        strict=True,
    )
    tier = cs.process_set.tier_of("iso_alpha", cs.parameters.tier_map())
    assert tier is Tier.SPECULATIVE
    # The sourced boil constants themselves are plausible (not speculative, not validated).
    assert cs.parameters["A_iso"].tier is Tier.PLAUSIBLE
    assert cs.parameters["Ea_iso"].tier is Tier.PLAUSIBLE
    # ... and the loss/utilization inputs are what drag the finished readout down.
    assert cs.parameters["k_iso_alpha_loss"].tier is Tier.SPECULATIVE
    assert cs.parameters["hop_utilization_efficiency"].tier is Tier.SPECULATIVE
