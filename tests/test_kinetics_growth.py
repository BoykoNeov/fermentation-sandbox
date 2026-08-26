"""Tests for GrowthNitrogenLimited — the first validated-core kinetic Process.

Covers the closed-form derivative, the Monod shutoffs (no N / no sugar), and the
two atom balances the Process is built to conserve: nitrogen (free YAN + biomass
N) and carbon (sugar + biomass C). The carbon balance is exercised on beer too,
since that is the only path through the per-slot vector carbon draw.
"""

from pathlib import Path

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics import GrowthNitrogenLimited
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import CompiledScenario
from fermentation.units import cells_per_ml_to_pitch_gpl
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)


@pytest.fixture
def store():
    # Real (sourced) wine parameters; the kinetics are medium-agnostic, so the
    # wine file suffices to exercise the Process mechanism (not beer-specific values).
    return load_parameters(default_data_dir() / "wine_generic.yaml")


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema, *, x: float = 0.1, s: float = 264.0, n: float = 0.3
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": 0.0, "N": n, "T": 293.15, "CO2": 0.0, "X_dead": 0.0})


def test_metadata():
    g = GrowthNitrogenLimited()
    assert g.name == "growth_nitrogen_limited"
    assert g.tier is Tier.PLAUSIBLE
    assert set(g.touches) == {"X", "S", "N"}
    # `reads` documents the params consumed (parameter-tier-propagation task).
    assert set(g.reads) == {"mu_max", "K_s", "K_n", "biomass_N_fraction", "biomass_C_fraction"}


def test_derivative_matches_closed_form(params):
    # Pin dX, dN, dS to the formula at a known wine state — no solver fuzz.
    schema = wine_schema()
    y = schema.pack(
        {"X": 1.0, "S": [200.0], "E": 0.0, "N": 0.3, "T": 293.15, "CO2": 0.0, "X_dead": 0.0}
    )
    d = GrowthNitrogenLimited().derivatives(0.0, y, schema, params)

    mu = params["mu_max"] * (200.0 / (params["K_s"] + 200.0)) * (0.3 / (params["K_n"] + 0.3))
    dx = mu * 1.0
    assert schema.get(d, "X") == pytest.approx(dx)
    assert schema.get(d, "N") == pytest.approx(-params["biomass_N_fraction"] * dx)
    # Carbon drawn from sugar closes exactly: c_frac_glucose * |dS| == f_C * dX.
    ds = schema.get(d, "S")
    assert ds < 0.0
    assert carbon_mass_fraction("glucose") * (-ds) == pytest.approx(
        params["biomass_C_fraction"] * dx
    )
    # Growth touches X, S, N only — never E or CO2.
    assert schema.get(d, "E") == 0.0
    assert schema.get(d, "CO2") == 0.0


def test_no_growth_without_nitrogen(params):
    schema = wine_schema()
    d = GrowthNitrogenLimited().derivatives(0.0, _wine_y0(schema, n=0.0), schema, params)
    assert np.array_equal(d, schema.zeros())


def test_no_growth_without_sugar(params):
    schema = wine_schema()
    d = GrowthNitrogenLimited().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)
    assert np.array_equal(d, schema.zeros())


def test_wine_run_conserves_carbon_and_nitrogen(params, store):
    schema = wine_schema()
    # strict=True enforces the touch contract on every solver step.
    ps = ProcessSet(schema, [GrowthNitrogenLimited()], strict=True)
    traj = simulate(ps, params=params, y0=_wine_y0(schema), t_span=(0.0, 300.0))
    assert traj.success

    f_c = store.value("biomass_C_fraction")
    f_n = store.value("biomass_N_fraction")
    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=f_c),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_conserved(
        traj,
        total_nitrogen(schema, biomass_nitrogen_fraction=f_n),
        rtol=1e-5,
        atol=1e-6,
        label="nitrogen",
    )
    assert_nonnegative(traj, ("X", "S", "N"), atol=1e-7)


def test_biomass_grows_then_caps_when_nitrogen_exhausted(params, store):
    schema = wine_schema()
    ps = ProcessSet(schema, [GrowthNitrogenLimited()], strict=True)
    x0, n0 = 0.1, 0.3
    traj = simulate(ps, params=params, y0=_wine_y0(schema, x=x0, n=n0), t_span=(0.0, 500.0))

    x = traj.series("X")
    n = traj.series("N")
    assert x[-1] > x[0]  # biomass grew (exponential phase)
    assert n[-1] < 1e-3  # YAN essentially exhausted (stationary phase)
    # Biomass is nitrogen-capped: N + f_N*X is conserved, so as N -> 0 the
    # biomass plateaus at X0 + N0/f_N. This IS the "stop dividing" mechanism.
    f_n = store.value("biomass_N_fraction")
    assert x[-1] == pytest.approx(x0 + n0 / f_n, rel=1e-3)


def test_beer_run_conserves_carbon_through_vector_draw(params, store):
    # Beer's three-slot sugar is the only path exercising the per-slot carbon
    # draw; carbon must still close to solver tolerance.
    schema = beer_schema()
    ps = ProcessSet(schema, [GrowthNitrogenLimited()], strict=True)
    y0 = schema.pack(
        {
            "X": 0.1,
            "S": [30.0, 60.0, 10.0],
            "E": 0.0,
            "N": 0.25,
            "T": 293.15,
            "CO2": 0.0,
            "X_dead": 0.0,
        }
    )
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 300.0))
    assert traj.success

    f_c = store.value("biomass_C_fraction")
    assert_conserved(
        traj,
        total_carbon(schema, biomass_carbon_fraction=f_c),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_nonnegative(traj, ("X", "S", "N"), atol=1e-7)


# ---------------------------------------------------------------------------
# D-211: beer's growth RATE, validated against Tyrell 2013 Fig. 4's cell counts
# ---------------------------------------------------------------------------
#
# Fig. 4 has three panels from one EBC-tube trial at 15 C: extract, pH and total
# cell count. D-180 used the extract panel and D-207 the pH panel; the cell-count
# panel went untranscribed for four beats even though it is the only measured
# growth curve in reach. D-211 transcribed it (pixel calibration on the printed
# gridlines, markers found by width, the unfitted 20-gridline predicted to
# 0.018e6/mL) and re-derived ``mu_max`` from it.
#
# The comparison is in NORMALISED growth fraction, never in cells: cells/mL and
# g/L differ by a cell-mass conversion nothing here sources, and the
# normalisation cancels it exactly. Uptake is growth-coupled and
# ``growth_nitrogen_limited`` is the ONLY active beer Process touching ``N``
# (checked below rather than assumed), so this curve is also the nitrogen
# drawdown curve that sets beer's pH timing.

#: Tyrell 2013 Fig. 4, cell-count panel: four-strain envelope in 10^6 cells/mL.
TYRELL_CELL_COUNT = {
    0: (9.98, 9.94),
    1: (15.10, 19.77),
    2: (27.07, 32.13),
    3: (29.06, 34.69),
    4: (15.96, 26.91),  # the settling limb: counts FALL, which this model cannot do
}
TYRELL_PITCH_CELLS = 9.96
TYRELL_PEAK_DAY = 3

#: The wort nitrogen Tyrell's trial is modelled at, mg N/L. **An ASSUMPTION** - the paper prints
#: none (the wort is a dilution of Bavarian Pilsener malt extract and the strains' amino-acid
#: profile is "data not shown"). Named rather than left as a literal at D-230, so the guard that
#: scores it against a sourced malt-wort composition reads the same object the run does. The
#: aroma calibration frame in ``test_kinetics_byproducts`` assumes the same 200 for the same
#: reason; D-230 checked the assumption and did NOT move it, in either frame.
TYRELL_ASSUMED_YAN_MGL = 200.0
TYRELL_CELL_SUGAR_GPL = 82.2388545


def _tyrell_growth_fraction_spread(day: int) -> tuple[float, float]:
    """Measured growth fraction at ``day``, normalised on the peak-day MIDPOINT.

    Both edges use the same denominator, so the spread this returns is the
    strain-to-strain scatter and nothing else.
    """
    mid = {d: (lo + hi) / 2.0 for d, (lo, hi) in TYRELL_CELL_COUNT.items()}
    span = mid[TYRELL_PEAK_DAY] - TYRELL_PITCH_CELLS
    lo, hi = TYRELL_CELL_COUNT[day]
    return (lo - TYRELL_PITCH_CELLS) / span, (hi - TYRELL_PITCH_CELLS) / span


def _tyrell_beer_growth(
    data_dir: Path | None = None,
) -> tuple["CompiledScenario", FloatArray, FloatArray, FloatArray]:
    from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

    compiled = compile_scenario(
        Scenario(
            name="tyrell-growth",
            medium="beer",
            initial={
                "glucose_gpl": 0.15 * TYRELL_CELL_SUGAR_GPL,
                "maltose_gpl": 0.70 * TYRELL_CELL_SUGAR_GPL,
                "maltotriose_gpl": 0.15 * TYRELL_CELL_SUGAR_GPL,
                "yan_mgl": TYRELL_ASSUMED_YAN_MGL,
                # Tyrell's OWN counted pitch (D-222). It was a flat 1.0 g/L — 2.51x the
                # biomass Tyrell pitched — until D-222 converted the paper's count through
                # the settled 40 pg/cell. `mu_max` is fitted ON this curve, so the fit frame
                # has to be the trial's own inoculum or the rate absorbs the pitch error.
                "pitch_gpl": cells_per_ml_to_pitch_gpl(TYRELL_PITCH_CELLS * 1e6),
                "initial_ph": 5.65,
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=14.0,
        ),
        data_dir=data_dir,
    )
    res = compiled.run()
    t_h = np.asarray(res.t, dtype=float)
    x = np.asarray(res.y, dtype=float)[compiled.schema.slice("X").start, :]
    return compiled, t_h, x, (x - x[0]) / (x.max() - x[0])


def _beer_params_with_mu(dest: Path, value: float, band: tuple[float, float] | None = None) -> Path:
    """Copy the packaged parameter dir into ``dest`` with beer ``mu_max`` replaced.

    The value and its band move together because the Parameter schema rejects a
    value outside its own uncertainty band — an arm that moved only one would not
    load, which is a loud failure rather than a silent half-change.
    """
    import re
    import shutil

    src = default_data_dir()
    shutil.copytree(src, dest)
    lo, hi = band if band is not None else (value, value)
    text = (dest / "beer_generic.yaml").read_text(encoding="utf-8")
    # Matched by pattern, not by formatting the shipped float back into a string:
    # the YAML writes 0.040 and Python renders the same number as 0.04, so an exact
    # string replace silently matches nothing and the arm runs as the baseline.
    text, n = re.subn(r"^(mu_max:\n  value: )[-0-9.eE]+$", rf"\g<1>{value!r}", text, flags=re.M)
    assert n == 1, f"mu_max value pattern matched {n} times"
    text, n = re.subn(
        r"^(mu_max:\n(?:  (?!uncertainty)[^\n]*\n)*  uncertainty: )\{[^\n]*\}$",
        rf'\g<1>{{ low: {lo!r}, high: {hi!r}, note: "test arm" }}',
        text,
        flags=re.M,
    )
    assert n == 1, f"mu_max uncertainty pattern matched {n} times"
    (dest / "beer_generic.yaml").write_text(text, encoding="utf-8")
    loaded = load_parameters(dest / "beer_generic.yaml")["mu_max"]
    assert loaded.value == value, f"mu_max override did not take: {loaded.value}"
    return dest


def test_only_growth_draws_beer_nitrogen():
    # The whole comparison rests on the cell-count curve BEING the nitrogen
    # drawdown curve. That is a claim about the registry, not about one Process,
    # so it is checked against the active set rather than read off a docstring.
    compiled, _, _, _ = _tyrell_beer_growth()
    drawers = [p.name for p in compiled.process_set.active if "N" in p.touches]
    assert drawers == ["growth_nitrogen_limited"], (
        f"beer's `N` is touched by {drawers}; if another Process joined, the cell-count "
        "curve is no longer the uptake curve and D-211's calibration needs re-deriving"
    )


def test_beer_growth_rate_matches_tyrells_measured_cell_counts():
    # The calibration D-211 shipped, scored on the observable it was fitted to.
    # Days 1 and 2 are the only informative points: day 0 is identically 0 and
    # day 3 identically 1 under peak normalisation, and day 4 is the settling
    # limb the model does not represent.
    _, t_h, _, frac = _tyrell_beer_growth()
    for day in (1, 2):
        lo, hi = _tyrell_growth_fraction_spread(day)
        model = float(np.interp(day * 24.0, t_h, frac))
        assert lo <= model <= hi, (
            f"day {day}: model growth fraction {model:.4f} outside Tyrell's four-strain "
            f"spread {lo:.4f}-{hi:.4f}. `mu_max` is calibrated ON this curve (D-211), so a "
            "failure here means the growth path moved, not that the data did"
        )


@pytest.mark.parametrize("edge", ["low", "high"])
def test_both_mu_max_band_edges_stay_inside_the_measured_spread(edge, tmp_path):
    # The band is not decoration: D-211 CONSTRUCTED its edges as the admissible
    # set of `mu_max` on this very spread, so both edges owe the same check the
    # nominal gets. Pinning only the nominal would leave the edges free to drift
    # to values the data rules out.
    shipped = load_parameters(default_data_dir() / "beer_generic.yaml")["mu_max"]
    val = shipped.uncertainty.low if edge == "low" else shipped.uncertainty.high
    dest = _beer_params_with_mu(tmp_path / "params", val)

    _, t_h, _, frac = _tyrell_beer_growth(data_dir=dest)
    for day in (1, 2):
        lo, hi = _tyrell_growth_fraction_spread(day)
        model = float(np.interp(day * 24.0, t_h, frac))
        assert lo <= model <= hi, (
            f"{edge} band edge mu_max={val}: day {day} growth fraction {model:.4f} outside "
            f"the measured spread {lo:.4f}-{hi:.4f} the band was built from"
        )


def test_the_retired_growth_rate_is_ruled_out_by_the_counts(tmp_path):
    # THE MUTATION THE CLAIM NAMES. "The prior 0.098/h was ~3x too fast" is only
    # worth writing down if running it fails, so this runs it. It is also the guard
    # that stops the value drifting back toward Zamudio's Droop-form magnitude.
    dest = _beer_params_with_mu(tmp_path / "params", 0.098)
    _, t_h, _, frac = _tyrell_beer_growth(data_dir=dest)
    lo, hi = _tyrell_growth_fraction_spread(1)
    day1 = float(np.interp(24.0, t_h, frac))
    assert day1 > hi, (
        f"at the retired 0.098/h the day-1 growth fraction is {day1:.4f}, which the counts' "
        f"{lo:.4f}-{hi:.4f} would ADMIT. D-211's re-derivation rests on it being ruled out"
    )
    # RE-PINNED at D-222. "The crop is at its ceiling by day 1" was a pitch-1.0 statement:
    # at Tyrell's own counted pitch the same retired rate has 69 % of the crop grown by day 1,
    # because the ceiling is unchanged in ABSOLUTE terms while the fold is 2.51x larger. The
    # ruling-out above is unaffected — 0.690 is still well above the measured 0.448 — and it is
    # still most of the nitrogen charge step arriving before the day-1 pH reading (D-209 §7).
    assert day1 == pytest.approx(0.690, abs=0.02), (
        f"at the retired 0.098/h the day-1 growth fraction is {day1:.4f}; D-222 measured 0.690 "
        "at Tyrell's counted pitch (it was 1.000, the ceiling, at the retired 1.0 g/L pitch)"
    )


def test_the_measured_multiplication_is_not_reproduced_and_that_is_a_separate_defect():
    """A RECORDED DEVIATION, pinned so it is not mistaken for something a rate fixed.

    **D-222 INVERTED THIS TEST'S SIGN, and that inversion is the beat.** D-211 measured the
    model's multiplication at 2.75x against Tyrell's 2.92-3.48x and called the extent "right,
    ~6 % low". That was measured at a 1.0 g/L pitch — 2.51x the biomass Tyrell counted — and
    the fold is precisely what a pitch SETS. Beer's growth is nitrogen-limited, so the gain is
    fixed in ABSOLUTE terms (``dX = YAN / biomass_N_fraction`` ~= 1.75 g/L at either pitch);
    dividing the same gain by an inoculum 2.51x lighter gives 5.39x, **1.55x ABOVE** the
    measured envelope's high edge. The old agreement was an artefact of the excess pitch.

    **No growth RATE can repair it** — ``mu_max`` sets WHEN the ceiling is reached, not WHERE
    it is; the test below runs that claim rather than asserting it. That is why D-222 refits
    the rate and REPORTS the extent instead of tuning it, and why the timing comparison above
    is normalised on each curve's own peak.

    **What it would take, stated as a quantity rather than a mood.** Tyrell's counted growth is
    (29.06-34.69) - 9.96 = 19.1-24.7 x10^6 cells/mL, i.e. 0.76-0.99 g/L at the settled 40 pg,
    so at this file's ``biomass_N_fraction`` only 87-113 mg/L of nitrogen reached suspended
    biomass — against the 200 mg/L the scenario assumes.

    **D-230 SETTLES WHICH CANDIDATE THAT IS, and it is neither of the two D-222 named.** Tyrell
    still prints no FAN, but the repo already carries a sourced 10-12 degP malt-wort composition
    (Peyer 2017, the source ``nitrogen_uptake_charge_beer`` is derived from) and it puts a malt
    wort's ASSIMILABLE nitrogen at 189-194 mg N/L. The assumed 200 is 3-6 % above that envelope,
    not 1.55x, so the nitrogen assumption cannot carry the gap and was NOT moved. The partition
    candidate is refused by arithmetic rather than by sourcing: see the two tests below, which
    show a nitrogen-budget repair would need cell nitrogen at 20-26 % of dry weight. The residue
    is a frame ambiguity with both branches priced, not a defect in this Process.
    """
    _, _, x, _ = _tyrell_beer_growth()
    model_fold = float(x.max() / x[0])
    meas_lo = TYRELL_CELL_COUNT[TYRELL_PEAK_DAY][0] / TYRELL_PITCH_CELLS
    meas_hi = TYRELL_CELL_COUNT[TYRELL_PEAK_DAY][1] / TYRELL_PITCH_CELLS
    assert model_fold > meas_hi, (
        f"model multiplication {model_fold:.3f}x is no longer ABOVE Tyrell's measured "
        f"{meas_lo:.3f}-{meas_hi:.3f}x. D-222's finding is that the extent agreement D-211 "
        "recorded belonged to the 2.51x pitch excess; if it is back, say which repair did it"
    )
    assert model_fold / meas_hi == pytest.approx(1.546, abs=0.02), (
        f"the extent overshoot is {model_fold / meas_hi:.3f}x the measured high edge; D-222 "
        "measured 1.546"
    )


def test_no_growth_rate_in_the_band_can_repair_the_extent(tmp_path):
    """The mutation the claim above names: a RATE cannot move a nitrogen-limited CEILING.

    Without this, "no growth rate can repair it" is a mechanism assertion resting on the
    reader's agreement rather than on a run. Both re-constructed band edges are integrated and
    the fold is read off each; a RED means growth has stopped being nitrogen-limited, which
    would re-open D-211 §1 — the identity that makes the cell-count panel a nitrogen curve.
    """
    shipped = load_parameters(default_data_dir() / "beer_generic.yaml")["mu_max"]
    meas_hi = TYRELL_CELL_COUNT[TYRELL_PEAK_DAY][1] / TYRELL_PITCH_CELLS
    folds = []
    for edge, val in (("low", shipped.uncertainty.low), ("high", shipped.uncertainty.high)):
        dest = _beer_params_with_mu(tmp_path / f"params-{edge}", val)
        _, _, x, _ = _tyrell_beer_growth(data_dir=dest)
        folds.append(float(x.max() / x[0]))
        assert folds[-1] > meas_hi, (
            f"at the {edge} band edge mu_max={val} the multiplication is {folds[-1]:.3f}x, "
            f"inside or below Tyrell's measured envelope (high edge {meas_hi:.3f}x). A rate "
            "would then be a lever on the extent after all"
        )
    assert abs(folds[1] - folds[0]) < 0.05, (
        f"the fold moves {abs(folds[1] - folds[0]):.3f} across the whole band ({folds[0]:.3f} "
        f"to {folds[1]:.3f}); D-222 measured 0.010. A rate that moves the extent means growth "
        "is no longer nitrogen-limited"
    )


def test_the_assumed_wort_nitrogen_is_inside_a_sourced_malt_worts_envelope():
    """D-222's FIRST extent candidate, scored and CLOSED (D-230).

    ``yan_mgl = 200`` has been the beer scenarios' wort nitrogen since D-178 and had never been
    checked against anything — Tyrell prints no FAN, and "nothing sources a repair" is what the
    test above said for eight records. That was true about *Tyrell's* wort and false about *a*
    wort: this repo transcribed a 10-12 degP malt wort's full free-amino-acid composition at
    D-209 to derive ``nitrogen_uptake_charge_beer``, and the same table sizes the pool that
    parameter is a charge per mole OF. Nobody had summed it.

    **Why the sum is the assimilable pool with no correction owed.** Peyer's Table 16 has 18
    amino acids and proline is not one of them, which is exactly the right scope: proline is
    Jones & Pierce Group D, brewing yeast does not assimilate it, and the ``N`` slot is
    assimilable nitrogen by definition. A FAN figure would have needed a proline subtraction and
    this does not — that was the fork this test was written to resolve, because a FAN number
    sitting in an assimilable-N slot would have been a real scope error worth perhaps 20-30 %.

    **The verdict is that the assumption is sound, so nothing moves.** 200 sits 3.1-5.8 % above
    the sourced envelope. Moving it would re-anchor all eight beer aroma constants (their
    synthesis integral IS ``YAN / biomass_N_fraction``, D-226/D-228), ``Y_acetic_biomass_beer``
    (whose denominator is literally the biomass this wort's nitrogen builds), and the pH course
    (the ``N`` slot is the charge-balance cation seed) — to swap one wort's assumption for a
    DIFFERENT wort's measurement. ``acidbase.yaml`` already records Peyer's wort as "a DIFFERENT
    wort from Tyrell 2013's ... an accepted deviation, stated not tuned", and that is the right
    call here too. **Do not re-propose this as a free win.**

    **And the cascade is FIVE-way, not four — measured, and it was NOT predicted.** D-230's
    falsification arm A set this constant to the 113 mg/L that would land Tyrell's counts and
    pre-registered four REDs; it produced eight. The four unpredicted ones are every rate and
    timing test in this file (``test_beer_growth_rate_matches_tyrells_measured_cell_counts``,
    ``test_both_mu_max_band_edges_stay_inside_the_measured_spread``,
    ``test_the_retired_growth_rate_is_ruled_out_by_the_counts``,
    ``test_no_growth_rate_in_the_band_can_repair_the_extent``). The reason is structural: the
    growth fraction ``mu_max`` is fitted on is normalised on the peak, and the peak IS the
    nitrogen-limited ceiling, so wort nitrogen is not separable from the rate fit. Moving this
    constant re-opens D-222's refit as well — **the wort-nitrogen candidate was never an
    isolable knob**, which the arm establishes and the prediction missed.
    """
    from tests.conftest import (
        PEYER_WORT_AMMONIUM_MG_N_PER_L,
        peyer_wort_assimilable_nitrogen_mg_per_l,
    )

    low, high = peyer_wort_assimilable_nitrogen_mg_per_l()
    assert (low, high) == pytest.approx((189.0, 194.0), abs=0.5), (
        f"the sourced malt-wort envelope is {low:.1f}-{high:.1f} mg N/L; D-230 measured "
        "189-194. Both edges are the AMMONIUM range alone (Peyer Table 2's printed 25-30) — "
        "the amino-acid half is one composition column with no band of its own"
    )
    # The band is the ammonium row's width and nothing else, so it must be exactly that wide.
    assert high - low == pytest.approx(
        PEYER_WORT_AMMONIUM_MG_N_PER_L[1] - PEYER_WORT_AMMONIUM_MG_N_PER_L[0], abs=1e-9
    )

    excess = TYRELL_ASSUMED_YAN_MGL / high, TYRELL_ASSUMED_YAN_MGL / low
    assert max(excess) < 1.10, (
        f"the assumed {TYRELL_ASSUMED_YAN_MGL:.0f} mg N/L is {max(excess):.3f}x a sourced malt "
        f"wort's assimilable nitrogen ({low:.1f}-{high:.1f}). D-230 measured 1.031-1.058 and "
        "closed the wort-nitrogen candidate on it; above ~1.10 that closure stops holding and "
        "the extent overshoot would have a nitrogen component again"
    )
    assert min(excess) > 1.0, (
        "the assumption no longer sits ABOVE the sourced envelope. D-230's closure is that it "
        "overstates by 3-6 %, which is the WRONG SIGN to be hiding under an extent overshoot "
        "of 1.55x but the right sign to be honest about"
    )


def test_the_extent_overshoot_cannot_be_a_nitrogen_budget_error():
    """D-222's SECOND candidate, REFUSED by arithmetic rather than parked (D-230).

    The remaining candidate was "the model books one lumped N pool wholly into suspended
    biomass". Stated as a partition fraction it is unbuildable — it invents the one constant
    that sets the curve it would be fitted against, D-213 section 7's refusal. But it does not
    need to be built to be scored, because the growth law is an identity: ``dX = YAN / f_N``,
    with both inputs sourced. Invert it against Tyrell's counted crop and it prices ITSELF.

    **At the settled 40 pg/cell, Tyrell's cells would have to be 20-26 % nitrogen by dry
    weight.** Real yeast is 7-12 % and this engine's own ``biomass_N_fraction`` band tops out at
    14 %. So no admissible partition of the nitrogen budget reproduces the counted crop: the
    overshoot is not a nitrogen error, and a term draining nitrogen away from suspended biomass
    would have to drain about half of it to a destination nothing in the corpus names.

    The comparison is a RELATION to the shipped ``biomass_N_fraction`` and its band, not a
    hardcoded 0.14, so a later change to that parameter moves this guard instead of silently
    contradicting it (D-228's idiom for the pitch frame).
    """
    f_n = load_parameters(default_data_dir() / "beer_generic.yaml")["biomass_N_fraction"]
    pitch_gpl = cells_per_ml_to_pitch_gpl(TYRELL_PITCH_CELLS * 1e6)
    per_cell_g = pitch_gpl / (TYRELL_PITCH_CELLS * 1e6 * 1e3)  # the engine's gram, D-219

    yan_gpl = TYRELL_ASSUMED_YAN_MGL / 1000.0
    implied = []
    for peak in TYRELL_CELL_COUNT[TYRELL_PEAK_DAY]:
        counted_new_per_l = (peak - TYRELL_PITCH_CELLS) * 1e6 * 1e3
        dx_counted = counted_new_per_l * per_cell_g  # g/L of new biomass, at the engine's gram
        implied.append(yan_gpl / dx_counted)  # g N per g cell this would demand

    assert min(implied) > f_n.uncertainty.high, (
        f"reproducing Tyrell's counted crop at the engine's own per-cell mass demands cell "
        f"nitrogen of {min(implied):.3f}-{max(implied):.3f} g N/g, and the shipped band's high "
        f"edge is {f_n.uncertainty.high}. If this passes no longer, a nitrogen-budget repair "
        "has become admissible and D-230's refusal of the partition candidate re-opens"
    )
    assert (min(implied), max(implied)) == pytest.approx((0.202, 0.262), abs=0.005), (
        f"D-230 measured 0.202-0.262 g N/g cell; this run gives {min(implied):.3f}-"
        f"{max(implied):.3f}. Both edges are Tyrell's peak-day strain spread at one assumed YAN"
    )
    assert min(implied) / f_n.value > 1.7, (
        "the demanded nitrogen fraction is within 1.7x of the shipped one. D-230's refusal "
        "rests on the demand being physically impossible (real yeast 7-12 % N), not merely high"
    )


def test_the_extent_residue_is_a_two_way_frame_ambiguity_and_both_branches_are_priced():
    """What is LEFT once both nitrogen candidates are closed — and it is not one thing (D-230).

    With ``YAN`` sourced and ``f_N`` pinned by chemistry, the only free quantity in
    ``fold = (X0 + YAN/f_N) / X0`` is the count-to-gram conversion, and it enters twice: once
    in ``X0`` and once in reading Tyrell's crop. Two DIFFERENT readings close the gap and this
    beat cannot separate them, so both are pinned rather than one being chosen:

    * **Tyrell's cells are heavier than the engine's gram** — the sourced nitrogen budget
      demands 71-92 pg per counted new cell against the settled 40 (1.8-2.3x). Suggestive
      rather than decisive: it is an INDEPENDENT third estimate (Peyer's wort nitrogen times the
      Roels elemental fraction against Tyrell's own counts, with no Coleman input), and it lands
      near the ~100 pg D-219 retired as a back-computed residual. That agreement is not evidence
      the residual was right — it is a reason not to call this branch settled either.
    * **Or 44-56 % of the crop had already left suspension at the day-3 peak**, so the counted
      peak is not the cells made. No cell-mass change needed at all. This is consistent with
      Tyrell's OWN next reading: their counts fall 22-45 % over the single following day
      (``TYRELL_CELL_COUNT[4]``), a settling rate fast enough to have removed about half the
      crop by day 3 if it began, as flocculation does, when the sugar ran down. That 22-45 % is
      an ENVELOPE-EDGE fall (day 4's low edge against day 3's low edge, high against high), not
      a per-strain one — Tyrell's figure does not identify which strain holds which edge on
      which day, so it bounds the fall rather than measuring one strain's.

    **Neither branch is built here.** The first would move ``cells_per_ml_to_pitch_gpl``, which
    D-219 settled as the DEFINITION of this engine's biomass gram, on an inference rather than
    on the count-plus-weighing that record says is what would settle it. The second is the
    missing settling Process, whose rate this corpus does not print — the five beer texts
    describe flocculation only qualitatively (trigger and strain-dependence, no constant).

    What this test protects is the AMBIGUITY: a later beat that closes the extent gap must say
    which branch it closed, because a repair that lands the fold without naming one has assumed
    the other away.
    """
    f_n = load_parameters(default_data_dir() / "beer_generic.yaml")["biomass_N_fraction"].value
    pitch_gpl = cells_per_ml_to_pitch_gpl(TYRELL_PITCH_CELLS * 1e6)
    per_cell_g = pitch_gpl / (TYRELL_PITCH_CELLS * 1e6 * 1e3)
    dx_gpl = (TYRELL_ASSUMED_YAN_MGL / 1000.0) / f_n  # what the nitrogen builds, in grams

    counted_new = [
        (peak - TYRELL_PITCH_CELLS) * 1e6 * 1e3 for peak in TYRELL_CELL_COUNT[TYRELL_PEAK_DAY]
    ]
    per_cell_demanded = sorted(dx_gpl / n * 1e12 for n in counted_new)  # pg
    assert per_cell_demanded == pytest.approx([70.9, 91.9], abs=1.0), (
        f"branch 1: the sourced nitrogen budget demands {per_cell_demanded[0]:.1f}-"
        f"{per_cell_demanded[1]:.1f} pg per counted new cell; D-230 measured 70.9-91.9"
    )
    assert per_cell_demanded[0] / (per_cell_g * 1e12) > 1.5, (
        "branch 1 has collapsed onto the engine's own gram, which would mean the extent gap is "
        "gone. If a repair did that, D-230's ambiguity is resolved and this test should say how"
    )

    made = dx_gpl / per_cell_g  # cells/L the nitrogen builds AT the engine's gram
    settled_share = sorted(1.0 - n / made for n in counted_new)
    assert settled_share == pytest.approx([0.436, 0.565], abs=0.01), (
        f"branch 2: closing the gap by settling alone needs {settled_share[0] * 100:.1f}-"
        f"{settled_share[1] * 100:.1f} % of the crop already out of suspension at the peak; "
        "D-230 measured 43.6-56.5 %"
    )
    # The consistency check that makes branch 2 a live alternative rather than an escape hatch:
    # Tyrell's own day-3 -> day-4 fall must be of the same order as the share branch 2 needs.
    day3, day4 = TYRELL_CELL_COUNT[3], TYRELL_CELL_COUNT[4]
    one_day_fall = sorted((hi - lo) / hi for hi, lo in zip(day3, day4, strict=True))
    assert one_day_fall == pytest.approx([0.224, 0.451], abs=0.01), (
        f"Tyrell's measured one-day settling fall is {one_day_fall[0] * 100:.1f}-"
        f"{one_day_fall[1] * 100:.1f} %; D-230 read 22.4-45.1 % off TYRELL_CELL_COUNT[4]"
    )
    assert max(one_day_fall) > min(settled_share) * 0.5, (
        "the measured one-day fall is far too slow to have removed branch 2's share by day 3, "
        "which would make the settling branch untenable and leave the cell-mass branch alone. "
        "D-230's two-way ambiguity would then be a one-way finding"
    )
