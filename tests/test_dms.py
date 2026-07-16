"""Tests for the DMS-via-SMM-hydrolysis aging Process :class:`SMMHydrolysis` (decision D-102).

The aged-wine "truffle / black olive / cooked corn" odorant, and the first sulfur pool in the model
that is **not** autolysis-gated: DMS accumulates by spontaneous hydrolysis of the grape-borne
precursor during bottle aging, lees or no lees, so it carries its own anchor instead of
ratio-splitting a shared autolytic yield (the D-96 linchpin ``mercaptans`` could not satisfy —
D-101).

These tests pin the closed-form derivative; prove the properties the beat rests on — **1:1 transfer
in DMS-equivalents** (``dms_potential + dms`` invariant to machine precision), **warmer ⇒ faster**
(the load-bearing sourced direction), **monotonic accumulation** (no D-42-style stripping sink),
**O₂- and pH-independence** (a hydrolysis, and a *sourced null* below pH 5) — and confirm the
un-seeded no-op, the off-ledger isolability, and the speculative tier.

**The validation tests are the point of the file.** ``test_reproduces_the_amarone_45c_anchor`` and
``test_reproduces_the_35c_model_wine_anchor`` check the shipped constants against the two in-matrix
measurements they were fitted through; ``test_cellar_halflife_lands_in_years`` checks the
extrapolation the model actually uses against the third, independent observation (six real Syrah
wines cellared at 18 °C for 24 months) — the check that would have caught the D-101 activation
energy, and does catch Scheuren's.

**What the aged-wine comparison may and may not claim (an overclaim, corrected).**
``test_predicted_aging_dms_is_the_right_order_against_syrah_at_the_MATCHING_age`` first read
*"brackets the Syrah observation"* — and did not: it compared the model at **5 years** against data
at **t24 = 2 years**, and asserted only that two ranges **overlap**, which nearly any positive
conversion satisfies. **The assertion was weaker than the sentence around it** — the D-96 failure
mode (prove the mechanism, never pin the observable) wearing a new costume. Two confounds make a
quantitative bracket unavailable at *any* age: the model reports DMS **formed by aging** (it seeds
``dms = 0``) while those wines already held **29.9–314.9 µg/L at bottling**, and their totals
**fall** t12→t24 as DMS permeates out through the closure — a shape a monotonic-from-zero model
cannot produce. So the comparison claims **order of magnitude and direction**, nothing sharper.
``test_the_amarone_miss_is_recorded_not_tuned`` carries the same caveat and the same honest
strength: a *recorded finding*, not a defect (see D-102 and ``dms.yaml``).
"""

import math
from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.kinetics import SMMHydrolysis
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.units.convert import gpl_to_ugl, ugl_to_gpl

#: The shipped constants are READ from dms.yaml rather than restated here — so a re-anchoring of
#: the parameters cannot leave these tests asserting the old numbers (the D-100 lesson: a test that
#: hard-codes the value it should be reading is a test of itself).
_GAS_CONSTANT = 8.314  # J/(mol·K), for the closed-form Arrhenius the tests re-derive independently


@pytest.fixture
def params() -> dict[str, float]:
    # dms.yaml MERGED with wine_generic.yaml, exactly as the compile seam merges them — T_ref is
    # the shared Arrhenius reference every rate anchors on, not dms.yaml's to own.
    store = load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "dms.yaml"
    )
    return store.resolve()


def _state(schema: StateSchema, *, dms_potential_ugl: float, temp_k: float) -> FloatArray:
    y = schema.zeros()
    y[schema.slice("dms_potential")] = ugl_to_gpl(dms_potential_ugl)
    y[schema.slice("T")] = temp_k
    return y


def _rate(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> FloatArray:
    return SMMHydrolysis().derivatives(0.0, y, schema, params)


def _k_at(params: Mapping[str, float], temp_k: float) -> float:
    """The first-order constant at ``temp_k`` — re-derived here, not read from the Process."""
    return params["k_smm_hydrolysis"] * math.exp(
        -params["E_a_smm_hydrolysis"] / _GAS_CONSTANT * (1.0 / temp_k - 1.0 / params["T_ref"])
    )


def _converted_fraction(params: Mapping[str, float], temp_k: float, hours: float) -> float:
    """Closed-form first-order conversion — the quantity both anchor papers report."""
    return 1.0 - math.exp(-_k_at(params, temp_k) * hours)


# --------------------------------------------------------------------------------------
# The Process contract
# --------------------------------------------------------------------------------------


def test_touches_only_the_two_off_ledger_dms_slots():
    # The isolability claim in one assertion: SMMHydrolysis cannot move anything conserved,
    # anything oxidative, or anything another Process owns. No o2 (a hydrolysis, not an
    # oxidation), no carbon/nitrogen pool (both slots are off every ledger — the D-74 A420
    # argument), no S/E/CO2 (aging draws no sugar).
    assert SMMHydrolysis.touches == ("dms_potential", "dms")
    assert SMMHydrolysis.tier is Tier.SPECULATIVE


def test_strict_process_set_accepts_the_touches_contract(params):
    schema = wine_schema()
    pset = ProcessSet(schema, [SMMHydrolysis()], strict=True)
    y = _state(schema, dms_potential_ugl=338.0, temp_k=291.15)
    d = pset.total_derivatives(0.0, y, params)
    assert d[schema.slice("dms")][0] > 0.0


def test_tier_of_dms_is_speculative():
    schema = wine_schema()
    pset = ProcessSet(schema, [SMMHydrolysis()])
    assert pset.tier_of("dms") is Tier.SPECULATIVE


# --------------------------------------------------------------------------------------
# The closed form, and the properties the beat rests on
# --------------------------------------------------------------------------------------


def test_derivative_matches_the_closed_form_first_order_arrhenius_decay(params):
    schema = wine_schema()
    temp = 291.15
    y = _state(schema, dms_potential_ugl=338.0, temp_k=temp)
    d = _rate(y, schema, params)
    expected = _k_at(params, temp) * ugl_to_gpl(338.0)
    assert d[schema.slice("dms_potential")][0] == pytest.approx(-expected, rel=1e-12)
    assert d[schema.slice("dms")][0] == pytest.approx(expected, rel=1e-12)


def test_transfer_is_exactly_one_to_one_in_dms_equivalents(params):
    # THE conservation law of this beat, and the reason the pool is denominated in DMS-equivalents
    # rather than in SMM: what leaves the precursor enters the product with no yield parameter and
    # no molar-mass conversion, so dms_potential + dms is invariant to machine precision. Booking
    # the pool in grams of SMM instead would need a stoichiometric yield AND SMM's molar mass —
    # and SMM's salt form (iodide, in both source papers) makes that ambiguous.
    schema = wine_schema()
    for dmsp in (1.0, 119.0, 338.0, 958.4):
        for temp in (277.15, 291.15, 318.15):
            d = _rate(_state(schema, dms_potential_ugl=dmsp, temp_k=temp), schema, params)
            total = d[schema.slice("dms_potential")][0] + d[schema.slice("dms")][0]
            assert total == pytest.approx(0.0, abs=1e-20)


def test_warmer_is_faster_the_load_bearing_sourced_direction(params):
    # E_a > 0 is what the two in-matrix anchors establish and what every cellar-temperature claim
    # rides on. It is also the discriminator that kills the Scheuren transfer: not the SIGN (which
    # both agree on) but the MAGNITUDE — see test_cellar_halflife_lands_in_years.
    schema = wine_schema()
    rates = [
        _rate(_state(schema, dms_potential_ugl=338.0, temp_k=t), schema, params)[
            schema.slice("dms")
        ][0]
        for t in (277.15, 283.15, 291.15, 298.15, 308.15, 318.15)
    ]
    assert all(a < b for a, b in zip(rates[:-1], rates[1:], strict=True))
    assert params["E_a_smm_hydrolysis"] > 0.0


def test_accumulation_is_monotonic_no_stripping_sink(params):
    # The D-42 contrast: h2s has a CO2-stripping sink that makes it a RESIDUAL. Aging is
    # post-dryness, so there is no CO2 stream — dms only rises. d(dms)/dt >= 0 always, with no
    # clamp: it is a decay of a non-negative pool by a positive rate.
    schema = wine_schema()
    for dmsp in (0.0, 1e-6, 1.0, 958.4):
        d = _rate(_state(schema, dms_potential_ugl=dmsp, temp_k=291.15), schema, params)
        assert d[schema.slice("dms")][0] >= 0.0
        assert d[schema.slice("dms_potential")][0] <= 0.0


def test_no_precursor_is_a_byte_for_byte_no_op(params):
    # Prime directive #3 at the slot level: an un-seeded wine (or one whose precursor is spent)
    # contributes exactly zero, so the Process cannot perturb a run that has no DMS story.
    schema = wine_schema()
    d = _rate(_state(schema, dms_potential_ugl=0.0, temp_k=291.15), schema, params)
    assert not np.any(d)


def test_solver_undershoot_below_zero_is_absorbed(params):
    schema = wine_schema()
    y = _state(schema, dms_potential_ugl=0.0, temp_k=291.15)
    y[schema.slice("dms_potential")] = -1e-18  # a solver undershoot, not real physics
    assert not np.any(_rate(y, schema, params))


# --------------------------------------------------------------------------------------
# The compile seam — the must-seeding call, and the gap driving the model end-to-end found
# --------------------------------------------------------------------------------------


def test_dms_potential_is_seeded_from_the_sourced_level_not_zero():
    # THE anti-hard-zero call (D-102), pinned. `dms_potential` is the ONLY optional wine slot that
    # does not default to 0: so2_total/oak/anthocyanin are winemaking DOSES (0 is a true statement
    # about a scenario that added nothing), but DMSp is a property of the GRAPE that every must
    # carries — a 0 default would assert that aged wine develops no DMS, which is the D-45
    # hard-zero defect.
    from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

    cs = compile_scenario(
        Scenario(
            name="dms-seed",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    seeded_ugl = gpl_to_ugl(float(cs.y0[cs.schema.slice("dms_potential")][0]))
    assert seeded_ugl == pytest.approx(cs.parameters["dms_potential_initial"].value)
    assert seeded_ugl > 0.0  # the whole point
    assert float(cs.y0[cs.schema.slice("dms")][0]) == 0.0  # produced-only: no DMS at pitch


@pytest.mark.parametrize("override_ugl", [0.0, 119.0, 958.4])
def test_scenario_can_override_the_precursor_level(override_ugl):
    # REGRESSION. dms.yaml, the VarSpec and the compile seam all document `dms_potential_ugl` as
    # the per-scenario override — and it was NOT in the medium's allowed-initial-keys set, so every
    # scenario that tried to use it raised. All 16 unit tests passed regardless, because they build
    # state by hand and never cross the compile seam. Driving a real scenario end-to-end is what
    # found it. Explicit 0 must still be honoured (it makes the Process byte-for-byte inert), which
    # is why it is parametrized here.
    from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

    cs = compile_scenario(
        Scenario(
            name="dms-override",
            medium="wine",
            initial={
                "brix": 24.0,
                "yan_mgl": 250.0,
                "pitch_gpl": 0.25,
                "dms_potential_ugl": override_ugl,
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    assert gpl_to_ugl(float(cs.y0[cs.schema.slice("dms_potential")][0])) == pytest.approx(
        override_ugl
    )


def test_smm_hydrolysis_is_disabled_without_begin_aging():
    # The aging-axis gate (D-70): wired into the medium but DISABLED at compile, since aging is
    # inherently post-ferment. So a seeded-but-un-aged wine never converts any precursor.
    from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

    cs = compile_scenario(
        Scenario(
            name="dms-ungated",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    assert SMMHydrolysis.name in cs.process_set
    assert not cs.process_set.is_enabled(SMMHydrolysis.name)


def test_rate_is_independent_of_oxygen_and_ph(params):
    # TWO sourced nulls, deliberately asserted rather than left implicit.
    # (1) O2: this is a hydrolysis, not an oxidation — a sealed, inerted, sulfited bottle makes DMS
    #     exactly as fast (the D-83 thermal/O2-independent relationship).
    # (2) pH: De La Burgade et al. measured SMM degradation at pH 2.8 vs 3.8 and found NO
    #     significant difference — below pH 5 the nucleophilic-substitution mechanism's rate is
    #     pH-insensitive. Wine's whole range sits inside that regime, so a pH term would be a knob
    #     with a measured value of zero. THIS NULL IS ALSO WHY SCHEUREN IS UNUSABLE: his wort at
    #     pH 5.2 is in the OTHER mechanism's regime (see D-102).
    schema = wine_schema()
    base = _rate(_state(schema, dms_potential_ugl=338.0, temp_k=291.15), schema, params)[
        schema.slice("dms")
    ][0]
    for slot, value in (("o2", 5e-3), ("so2_total", 1e-1), ("tartaric", 6.0), ("malic", 3.0)):
        y = _state(schema, dms_potential_ugl=338.0, temp_k=291.15)
        y[schema.slice(slot)] = value
        assert _rate(y, schema, params)[schema.slice("dms")][0] == pytest.approx(base, rel=1e-12)


# --------------------------------------------------------------------------------------
# VALIDATION — the shipped constants against the measurements they came from
# --------------------------------------------------------------------------------------


def test_reproduces_the_amarone_45c_anchor(params):
    # Ugliano et al. 2023: real Amarone spiked with SMM converted 49-54% in ONE MONTH at 45 C.
    # This is one of the two points E_a was fitted through, so agreement is expected — the test
    # exists so that a future re-anchoring cannot silently break the fit it claims to have.
    assert 0.49 <= _converted_fraction(params, 318.15, 730.0) <= 0.54


def test_reproduces_the_35c_model_wine_anchor(params):
    # De La Burgade et al. 2025: model wine (pH 2.8 and 3.8) spiked with SMM degraded 45-49%
    # (crown-capped) at 35 C. The SOURCE CONTRADICTS ITSELF on the duration — methods say 3 months,
    # Table 1's caption and section 4.2 say 4 — which is exactly where E_a's 101-125 kJ/mol band
    # comes from (see dms.yaml). So the honest assertion is that the shipped band STRADDLES the
    # observation across both readings, not that it hits a single number.
    at_3_months = _converted_fraction(params, 308.15, 3 * 730.0)
    at_4_months = _converted_fraction(params, 308.15, 4 * 730.0)
    assert at_3_months < 0.49  # the fast-E_a end under-converts by 4 months...
    assert at_4_months > 0.45  # ...and the slow-E_a end over-converts by 3
    assert at_3_months <= 0.49 and at_4_months >= 0.45  # the observed 45-49% is bracketed


def test_cellar_halflife_lands_in_years_and_rejects_the_scheuren_transfer(params):
    # THE test that discriminates. Every rival activation energy agrees on the SIGN; they disagree
    # by orders of magnitude on the only question the model asks — how fast DMS appears in a
    # cellared bottle. Extrapolating the shipped constants to the Syrah study's 18 C cellar:
    half_life_yr = math.log(2) / _k_at(params, 291.15) / 8766.0
    assert 2.0 < half_life_yr < 7.0  # years — the bottle-aging timescale

    # For contrast, the rejected activation energies — compared FAIRLY, on D-101's own proposed
    # design: anchor k in-matrix on the Amarone 45 C measurement (which every candidate shares,
    # so it cannot be the difference) and take ONLY E_a from elsewhere. The rivals then differ in
    # exactly one thing, and it changes the answer by an order of magnitude.
    k_45 = _k_at(params, 318.15)  # the shared in-wine anchor

    def cellar_half_life_yr_with(e_a: float) -> float:
        k_18 = k_45 * math.exp(-e_a / _GAS_CONSTANT * (1.0 / 291.15 - 1.0 / 318.15))
        return math.log(2) / k_18 / 8766.0

    # Scheuren's real value: wort at pH 5.2 — the WRONG mechanism for wine, per his own discussion
    # ("a high pH value affects the intramolecular substitution ... a low pH value promotes the
    # nucleophile substitution"). It predicts a cellared wine essentially never develops DMS.
    assert cellar_half_life_yr_with(186_000.0) > 40.0  # decades — refuted by the Syrah observation
    # D-101's banked 128 +/- 37 was a SEARCH-SUMMARY FABRICATION of that same paper's number. It
    # lands nearer the truth than the value it misquoted — the right order for the wrong reason,
    # from a source that (read properly) forbids the transfer entirely.
    assert 5.0 < cellar_half_life_yr_with(128_000.0) < 10.0
    # The shipped band, ordered and both ends in years.
    assert 2.0 < cellar_half_life_yr_with(125_000.0) < 7.0
    assert 2.0 < cellar_half_life_yr_with(101_000.0) < 7.0
    assert cellar_half_life_yr_with(101_000.0) < cellar_half_life_yr_with(125_000.0)


def test_dmsp_decay_direction_matches_the_syrah_cellar_observation(params):
    # The independent check (NOT a fit point — the two are different observables, which is why the
    # 3-point Arrhenius is non-linear and why forcing them onto one line would give a spurious
    # ~61 kJ/mol; see D-102). Six real Syrah wines cellared at 18 C lost 61% of their DMSp in 24
    # months => a DMSp half-life of ~1.47 yr. DMSp is the TOTAL precursor pool (SMM is only 21-74%
    # of it) and is broader and more labile, so it must decay FASTER than SMM alone. The model's
    # SMM half-life must therefore sit ABOVE 1.47 yr — but still in years, not decades.
    smm_half_life_yr = math.log(2) / _k_at(params, 291.15) / 8766.0
    observed_dmsp_half_life_yr = math.log(2) / (-math.log(1 - 0.61) / (24 * 730.0)) / 8766.0
    assert observed_dmsp_half_life_yr == pytest.approx(1.47, abs=0.05)
    assert smm_half_life_yr > observed_dmsp_half_life_yr  # the expected direction
    assert smm_half_life_yr < 10.0  # but not so slow that bottle aging could not show it


def test_predicted_aging_dms_is_the_right_order_against_syrah_at_the_MATCHING_age(params):
    # WHAT THE MODEL SAYS, WITH NOTHING TUNED — and what that CANNOT be compared to.
    #
    # THIS TEST IS THE CORRECTED VERSION OF AN OVERCLAIM (D-102). It first read "brackets the Syrah
    # observation", comparing the model at 5 YEARS against Syrah data at t24 = 2 YEARS, and its
    # assertion (`lo < 399.5 and hi > 68.2`) merely checked that two ranges OVERLAP — which almost
    # any positive conversion satisfies. The assertion was weaker than the sentence around it: the
    # D-96 failure mode (prove the mechanism, never pin the observable) in a new costume.
    #
    # TWO CONFOUNDS make a quantitative bracket unavailable at ANY age, and they are the point:
    #  1. The model reports DMS *FORMED BY AGING* (it seeds dms = 0). The Syrah wines already held
    #     29.9-314.9 ug/L DMS *AT BOTTLING* ("DMS was already present in all wines at bottling"),
    #     so the observed t24 totals are dominated by DMS the model structurally starts without.
    #  2. The observed totals FALL from t12 to t24 (LR1: 539.8 -> 399.5) because DMS PERMEATES OUT
    #     through the closure — a shape a monotonic-from-zero model cannot produce at all.
    # Formed-DMS and total-DMS are not commensurate. Only the ORDER OF MAGNITUDE is checkable, and
    # that is therefore all this test asserts.
    conv = _converted_fraction(params, 291.15, 2 * 8766.0)  # t24 — the MATCHING age
    lo, hi = 119.0 * conv, 958.4 * conv
    assert 0.25 < conv < 0.35  # ~28.5% of the precursor converts in the Syrah study's 2 years

    # Right order of magnitude against the observed t24 totals (68.2-399.5): tens-to-hundreds of
    # ug/L, not single digits and not thousands. A real check (a 10x-wrong rate fails it) but a
    # WEAKER one than "brackets" — which is the honest strength given the confounds above.
    assert 10.0 < lo < 100.0
    assert 100.0 < hi < 1000.0

    # The one comparison the paper makes PAIRABLE by naming its wines, recorded but not asserted as
    # validation: LR1 has the max DMSp (958.4) and went 314.9 -> 539.8 (t12) -> 399.5 (t24), i.e. it
    # GROSSLY formed >= 224.9 before permeation removed ~140. The model forms ~273 from that
    # precursor in 2 years — the right size against gross formation, and ~3x the NET increase of
    # 84.6. Both readings are consistent with a model that omits permeation; neither is a bracket.
    assert 200.0 < hi < 350.0


def test_the_amarone_miss_is_recorded_not_tuned(params):
    # The model OVER-PREDICTS Amarone (2.9-64.3 ug/L, mean 27.9), whose Corvina fruit evidently
    # carries far less precursor than Syrah — the sourced default. THE MISS IS RECORDED, NOT TUNED
    # (D-102): backing dms_potential_initial out of Amarone's observed DMS would make OAV ~= 1 land
    # by construction rather than by evidence, and OAV ~= 1 *is* the low-precursor corner, so the
    # temptation is specific and real. The brief's "OAV ~= 1 in aged wine" is Amarone's MEAN
    # specifically; reality has no single value, and this pins that the model does not pretend to.
    #
    # CAVEATED THE SAME WAY as the Syrah comparison above: Amarone's at-bottling DMS is UNSTATED, so
    # "over-predicts" also compares formed-DMS to a total. The DIRECTION is the honest read (Corvina
    # carries less precursor than Syrah); the factor is softer than it looks.
    conv = _converted_fraction(params, 291.15, 5 * 8766.0)
    assert 119.0 * conv > 64.3  # even the lowest-precursor wine exceeds Amarone's observed maximum


def test_dms_oav_is_reported_and_crosses_threshold_in_aged_wine(params):
    # The claim the beat actually supports, stated as the test: not "OAV = 1", but "DMS crosses
    # its perceptual threshold during bottle aging, across a wide band set by the grape".
    from fermentation.sensory.oav import load_thresholds

    threshold_ugl = load_thresholds()["threshold_dms_wine"].value
    assert threshold_ugl == pytest.approx(25.0)

    # NB this is the model's own OAV — DMS FORMED BY AGING over its threshold. It is not offered as
    # a reproduction of any measured aged-wine OAV: real wines also carry at-bottling DMS (which
    # would raise it) and lose DMS through the closure (which would lower it), and this Process
    # models neither. What is claimed is the SHAPE: crosses threshold, and does so across a band
    # ~8x wide because the grape's precursor level is ~8x wide.
    conv = _converted_fraction(params, 291.15, 5 * 8766.0)
    oav_lo, oav_hi = 119.0 * conv / threshold_ugl, 958.4 * conv / threshold_ugl
    assert oav_lo > 1.0  # even the lowest-precursor wine crosses threshold by 5 years
    assert oav_hi > 10.0  # and the band is WIDE — an 8x spread within one variety
    assert oav_hi / oav_lo == pytest.approx(958.4 / 119.0)  # the spread IS the precursor's, exactly
