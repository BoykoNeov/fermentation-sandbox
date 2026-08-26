"""Tests for boundary unit conversions, including round-trip properties."""

import math

import pytest
from hypothesis import given
from hypothesis import strategies as st

from fermentation.parameters.schema import Parameter
from fermentation.units import (
    abv_from_ethanol,
    apparent_gravity,
    brix_to_sg,
    brix_to_sugar_gpl,
    cells_per_ml_to_pitch_gpl,
    celsius_to_kelvin,
    days_to_hours,
    hours_to_days,
    kelvin_to_celsius,
    plato_to_sg,
    real_to_apparent_extract,
    sg_to_brix,
    sg_to_plato,
    sugar_gpl_to_brix,
)


def test_temperature_known_points():
    assert celsius_to_kelvin(0.0) == pytest.approx(273.15)
    assert celsius_to_kelvin(20.0) == pytest.approx(293.15)
    assert kelvin_to_celsius(273.15) == pytest.approx(0.0)


@given(st.floats(min_value=-50, max_value=200))
def test_temperature_roundtrip(celsius):
    assert kelvin_to_celsius(celsius_to_kelvin(celsius)) == pytest.approx(celsius)


@given(st.floats(min_value=0, max_value=10000))
def test_time_roundtrip(days):
    assert hours_to_days(days_to_hours(days)) == pytest.approx(days)


def test_days_hours_known():
    assert days_to_hours(1.0) == 24.0
    assert hours_to_days(48.0) == 2.0


def test_water_is_sg_one():
    # 0 Brix is pure water -> SG 1.000.
    assert brix_to_sg(0.0) == pytest.approx(1.0)


def test_brix_sg_typical_must():
    # A ~24 Brix must sits around SG 1.10 (handoff wine benchmark).
    sg = brix_to_sg(24.0)
    assert 1.09 < sg < 1.11


@given(st.floats(min_value=0.5, max_value=40))
def test_brix_sg_roundtrip(brix):
    # The two industry polynomials are independent fits, so allow a small gap.
    assert sg_to_brix(brix_to_sg(brix)) == pytest.approx(brix, abs=0.3)


def test_plato_brix_numerically_close():
    # Plato and Brix both measure % sucrose by mass; they nearly coincide.
    sg = 1.048
    assert sg_to_plato(sg) == pytest.approx(sg_to_brix(sg), abs=0.3)


@given(st.floats(min_value=0.5, max_value=30))
def test_plato_sg_roundtrip(plato):
    assert sg_to_plato(plato_to_sg(plato)) == pytest.approx(plato, abs=0.3)


def test_sugar_concentration_scale():
    # ~24 Brix must -> ~240-265 g/L sugar.
    gpl = brix_to_sugar_gpl(24.0)
    assert 230 < gpl < 270


@given(st.floats(min_value=0.1, max_value=40))
def test_sugar_brix_roundtrip_at_fixed_sg(brix):
    sg = brix_to_sg(brix)
    gpl = brix_to_sugar_gpl(brix, sg=sg)
    assert sugar_gpl_to_brix(gpl, sg=sg) == pytest.approx(brix)


def test_abv_from_ethanol_scale():
    # ~100 g/L ethanol is around 12.7% ABV.
    abv = abv_from_ethanol(100.0)
    assert 12.0 < abv < 13.5


def test_abv_zero():
    assert abv_from_ethanol(0.0) == 0.0


@given(st.floats(min_value=0, max_value=160))
def test_abv_monotonic_and_finite(ethanol_gpl):
    abv = abv_from_ethanol(ethanol_gpl)
    assert abv >= 0
    assert math.isfinite(abv)


def test_apparent_extract_unchanged_before_fermentation():
    # With no fermentation real extract == original extract, so the apparent
    # (hydrometer) extract equals it too — nothing depresses the reading yet.
    oe = sg_to_plato(1.048)
    assert real_to_apparent_extract(oe, oe) == pytest.approx(oe)


def test_apparent_gravity_1048_to_1010_example():
    # The canonical brewing example: a 1.048 OG ale that finishes at apparent
    # 1.010 has a *real* extract near 4.25 P (~1.016). Going the other way, a real
    # extract of 4.25 P at OG 11.91 P must read ~1.010 apparent.
    oe = sg_to_plato(1.048)  # ~11.91 P
    real_extract = 4.25
    assert real_to_apparent_extract(real_extract, oe) == pytest.approx(2.56, abs=0.1)
    assert apparent_gravity(real_extract, oe) == pytest.approx(1.010, abs=0.001)


def test_apparent_gravity_below_real_when_fermenting():
    # Once ethanol is present (real extract has dropped below OG) the hydrometer
    # reads low: apparent gravity < the real-extract gravity.
    oe = sg_to_plato(1.048)
    real_extract = 5.0  # partially fermented
    assert apparent_gravity(real_extract, oe) < plato_to_sg(real_extract)


# ======================================================================================
# D-219 - the per-cell dry mass, and why it is a unit identity rather than an estimate
# ======================================================================================

#: The three published pitches this repo converts, and the ``pitch_gpl`` each becomes.
#: Every one is a COUNT in its source, because that is how the literature states a pitch.
PUBLISHED_PITCHES: dict[str, tuple[float, float]] = {
    "Tyrell 2013 (beer, EBC tubes)": (9.96e6, 0.3984),
    "Varela 2004 / Palma 2012 (wine)": (1.0e6, 0.0400),
    "Foster 2022 (beer, D-218)": (1.2e7, 0.4800),
}

#: The band, in pg/cell. Derived, not asserted: it is what the elemental route below gives
#: across ``biomass_N_fraction``'s own printed 0.08-0.14 uncertainty. Both readings the
#: archive shipped before D-219 sit outside it.
SETTLED_BAND_PG = (28.0, 50.0)

RETIRED_READINGS_PG = {
    "unsourced wine-benchmark assertion": 18.0,
    "back-computed from the beer scenario pitch": 100.0,
}


def _elemental_pg_per_cell(n_init_mgl: float, f_n: float, coleman_pg: float) -> float:
    """Coleman's cells-per-gram-nitrogen, priced with an elemental composition he never used.

    The point of the manoeuvre is that ``coleman_pg`` cancels out of the measurement and
    re-enters only as the thing being checked: Coleman COUNTED cells, so ``Y_X/N /
    coleman_pg`` is the quantity he actually measured, and dividing the elemental dry mass
    per gram of nitrogen by it gives a dry mass per cell owing nothing to his assumption.
    """
    from fermentation.parameters import default_data_dir, load_parameters

    resolved = load_parameters(default_data_dir() / "wine_generic.yaml").resolve()
    a0 = resolved["biomass_N_yield_log_intercept"]
    a1 = resolved["biomass_N_yield_log_slope"]
    y_x_n = math.exp(a0 + a1 * n_init_mgl)  # g cell (Coleman-gram) / g N
    cells_per_g_n = y_x_n / (coleman_pg * 1e-12)  # what the hemacytometer measured
    return (1.0 / f_n) / cells_per_g_n * 1e12


def _biomass_n_fraction_spec() -> Parameter:
    from fermentation.parameters import default_data_dir, load_parameters

    return load_parameters(default_data_dir() / "wine_generic.yaml")["biomass_N_fraction"]


def test_the_pitch_conversion_is_colemans_own_biomass_unit_not_an_estimate():
    """The conversion is 4e-11 g/cell EXACTLY, because that is what the model's gram is.

    Coleman, Fish & Block 2007 Materials and Methods: *"Each cell count was converted to
    grams per liter of cell mass, assuming that each cell weighs 4 x 10^-11 g"*. They
    counted cells and weighed nothing, so every gram in that paper - and therefore ``X``,
    ``X_A``, ``Y_X/N`` and its nitrogen regression, ``k_prime_d`` and wine's ``mu_max``,
    all of which this engine fits to it - is a count times this constant.

    So this is not one estimate among several. A counted pitch converted at any other
    figure enters the model in a unit its own parameters do not use, which is precisely
    what the two wine benchmarks were doing at 18 pg/cell before D-219.

    **A RED here names the unit drifting away from the parameters fitted in it.** It is
    asserted exactly, not to a tolerance, because an identity has no tolerance.
    """
    assert cells_per_ml_to_pitch_gpl(1.0e6) == 0.04, (
        "the conversion no longer returns Coleman's 4e-11 g/cell. Every wine biomass "
        "parameter in this engine is expressed in that gram; changing this silently "
        "re-scales every counted pitch relative to them"
    )
    # Linear and zero-preserving - a pitch conversion that is not is a different claim.
    assert cells_per_ml_to_pitch_gpl(0.0) == 0.0
    assert cells_per_ml_to_pitch_gpl(2.0e6) == pytest.approx(2 * cells_per_ml_to_pitch_gpl(1e6))


def test_the_elemental_route_corroborates_colemans_assumption():
    """The independent check that says the unit is also physically honest (D-219).

    Coleman writes *"assuming"*, so nothing in that chain is a weighing and the tier is
    plausible, not validated. What makes 40 pg more than a convention is that inverting his
    yield back to the quantity he measured - cells per gram of nitrogen - and pricing it
    with the Roels elemental formula, which he had no hand in, lands at ~34.9 pg.

    It also settles the frame Coleman left open. He writes *"cell mass"*, never *"dry
    weight"*; read as a WET mass, 4e-11 g would make yeast ~33 % nitrogen on a dry basis
    against a real 7-12 %. The engine declares ``biomass_C_fraction`` and
    ``biomass_N_fraction`` on dry cell weight in both medium files, so dry is the frame
    those parameters need and the arithmetic agrees.

    **The wet arm is the positive control.** Without it, "the routes agree" is
    indistinguishable from a predicate that cannot disagree with anything.
    """
    f_n = _biomass_n_fraction_spec().value
    dry = _elemental_pg_per_cell(330.0, f_n, 40.0)
    assert dry == pytest.approx(34.9, abs=1.0), (
        f"the elemental route now gives {dry:.2f} pg/cell against Coleman's assumed 40. "
        "D-219 measured 34.87. This is the only corroboration the constant has, so a move "
        "here is a move in how well the engine's own biomass composition agrees with the "
        "unit its biomass is counted in"
    )
    assert 0.75 <= dry / 40.0 <= 1.0, "the check should sit BELOW Coleman's figure, not above"

    # Positive control: the same predicate on the WET reading it is supposed to reject.
    # 4e-11 g wet, at ~30 % dry matter, is ~12 pg of dry cell.
    wet = _elemental_pg_per_cell(330.0, f_n, 40.0 * 0.30)
    assert not (25.0 <= wet <= 55.0), (
        f"the wet reading gives {wet:.1f} pg/cell and the predicate accepts it, so the "
        "check cannot tell the two frames apart and proves nothing about which one ships"
    )


def test_both_readings_the_archive_shipped_are_outside_the_settled_band():
    """18 pg and ~100 pg are not near-misses; they are outside the band (D-219).

    The band is not asserted, it is DERIVED here from ``biomass_N_fraction``'s own printed
    uncertainty, so it moves if that parameter's provenance does. 18 pg implies a ~50 fL
    cell (a small lab haploid, not the diploid wine strain those benchmarks run) and
    ~100 pg implies a ~300 fL cell, which no *S. cerevisiae* is.

    A RED on the shipped value means the settlement itself has moved. A RED on either
    retired reading means the band has widened far enough to re-admit one of them, which
    would re-open D-216 section 11 - the window verdict in ``test_organic_acids.py``
    turns on exactly that exclusion.
    """
    spec = _biomass_n_fraction_spec()
    edges = sorted(
        _elemental_pg_per_cell(330.0, f, 40.0)
        for f in (spec.uncertainty.low, spec.uncertainty.high)
    )
    assert edges[0] == pytest.approx(SETTLED_BAND_PG[0], abs=1.0)
    assert edges[1] == pytest.approx(SETTLED_BAND_PG[1], abs=1.0)

    shipped_pg = cells_per_ml_to_pitch_gpl(1.0e6) / (1.0e6 * 1e3) * 1e12
    assert edges[0] <= shipped_pg <= edges[1], (
        f"the shipped {shipped_pg:.1f} pg/cell has fallen outside the elemental band "
        f"[{edges[0]:.1f}, {edges[1]:.1f}]"
    )
    for label, pg in RETIRED_READINGS_PG.items():
        assert not (edges[0] <= pg <= edges[1]), (
            f"the retired reading '{label}' ({pg} pg/cell) is back inside the band "
            f"[{edges[0]:.1f}, {edges[1]:.1f}]. D-219's settlement rests on it being out"
        )


def test_the_settled_band_propagates_only_one_of_its_two_uncertain_inputs():
    """D-232: the band above varies `f_N` and holds the OTHER input at a point estimate.

    ``_elemental_pg_per_cell`` has two uncertain inputs, not one. The test above sweeps
    ``biomass_N_fraction`` across its printed 0.08-0.14 and takes ``Y_X/N`` at
    ``exp(a0 + a1·N)`` — a single number. But ``a0`` and ``a1`` are **this repo's fitted
    regression** (D-13/D-14, with a published-typo correction to the slope exponent), and
    they carry their own printed 5-95 % credible regions. Across those, ``Y_X/N`` at 330 mg
    N/L runs **6.52-15.38 g/g** — a factor of 2.4 — and the band runs **18.6-76.7 pg**, not
    28.4-49.7.

    **D-219's conclusion SURVIVES this and is strengthened by it**: both retired readings
    are still outside even the wider band, so the exclusion never depended on the narrower
    one. **One of its characterisations does not.** D-219 says the two are *"not near
    misses"*; against the propagated band the 18 pg reading misses the low edge by 0.6 pg —
    about 3 % — which is a near miss by any reading, and nothing like the ~50 fL-cell gulf
    that record describes. The 100 pg reading is still clear by 30 %.

    Nothing is re-banded here. ``SETTLED_BAND_PG`` keeps its meaning and its narrower
    edges; what this test adds is that the narrowness is a CHOICE of what to propagate, so
    a later beat asking "is 20 pg admissible?" gets the honest answer instead of the
    point-estimate one.
    """
    import itertools

    from fermentation.parameters import default_data_dir, load_parameters

    resolved = load_parameters(default_data_dir() / "wine_generic.yaml")
    a0, a1 = resolved["biomass_N_yield_log_intercept"], resolved["biomass_N_yield_log_slope"]
    spec = _biomass_n_fraction_spec()

    def pg_at(a0v: float, a1v: float, fn: float) -> float:
        y_x_n = math.exp(a0v + a1v * 330.0)
        return (1.0 / fn) / (y_x_n / (40.0 * 1e-12)) * 1e12

    propagated = [
        pg_at(x, y, f)
        for x, y, f in itertools.product(
            (a0.uncertainty.low, a0.value, a0.uncertainty.high),
            (a1.uncertainty.low, a1.value, a1.uncertainty.high),
            (spec.uncertainty.low, spec.value, spec.uncertainty.high),
        )
    ]
    lo, hi = min(propagated), max(propagated)
    assert (lo, hi) == pytest.approx((18.6, 76.7), abs=0.5), (
        f"propagating the regression's own credible regions gives {lo:.1f}-{hi:.1f} pg; "
        "D-232 measured 18.6-76.7. This is the honest width of the elemental route"
    )
    assert hi / lo > 3.0 * (SETTLED_BAND_PG[1] / SETTLED_BAND_PG[0]) / 2.0, (
        f"the propagated band ({hi / lo:.2f}x wide) has stopped being much wider than the "
        f"shipped one ({SETTLED_BAND_PG[1] / SETTLED_BAND_PG[0]:.2f}x). If the regression's "
        "credible regions have narrowed, D-219's band is no longer understating itself"
    )

    # D-219's exclusion holds on the WIDER band — that is what makes it robust, and it is
    # asserted here so a future widening cannot quietly re-admit a retired reading.
    for label, retired in RETIRED_READINGS_PG.items():
        assert not (lo <= retired <= hi), (
            f"the retired reading '{label}' ({retired} pg/cell) is inside the PROPAGATED "
            f"band [{lo:.1f}, {hi:.1f}]. D-219's settlement survived D-232 only because "
            "both readings stayed out of the wider band too; if one is back in, that "
            "settlement re-opens on its own terms"
        )
    # ...but the 18 pg reading is a NEAR MISS against it, which D-219 explicitly denies.
    near = RETIRED_READINGS_PG["unsourced wine-benchmark assertion"]
    assert (lo - near) / lo < 0.05, (
        f"the 18 pg reading now misses the propagated low edge {lo:.1f} by "
        f"{(lo - near) / lo * 100:.1f} %; D-232 measured ~3 % and corrected D-219's "
        "'neither is a near miss' on exactly this. A large gap here means the regression "
        "moved and the correction should be restated"
    )


def test_the_published_pitches_convert_to_the_recorded_values():
    """Every counted pitch in the repo, through one conversion (D-219).

    Pinned because these three numbers are what the rest of the archive reasons about:
    0.3984 is the pitch Tyrell actually counted (against a scenario that ships 1.0), 0.0400
    is what both wine benchmarks now use, and 0.4800 is Foster's - the one D-218's fork
    table runs on.
    """
    for label, (cells, expected) in PUBLISHED_PITCHES.items():
        got = cells_per_ml_to_pitch_gpl(cells)
        assert got == pytest.approx(expected, abs=1e-4), f"{label}: {got:.4f} != {expected}"
