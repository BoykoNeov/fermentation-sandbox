"""What Crépin's own biomass data can and cannot settle about the inoculum (decision D-252).

D-249 §3 left the fixture's unsourced 0.25 g/L pitch as "the owner's call" and priced the move
onto the only sourced value, 0.04 g/L, as *"it buys the timing and sells the anchor"* — the
anchor being ``test_biomass_now_reaches_the_coleman_yield_the_compile_seam_installs``, whose
ratio it measured going 0.984 → 0.925, outside that guard's own ±0.02. D-251 then priced the
front-loading calibration against ``biomass_N_fraction``'s declared band and declined to take it
because the cells would carry 0.137–0.139 g N/g DW against a sourced 0.114.

**Both prices are void, and neither voiding moves the pitch.**

**The anchor was never a quantity this move could sell.** That guard computed
``predicted = 0.25 + initial / f_N`` with the fixture's pitch as a *literal*, so at any other
pitch it subtracted an inoculum the run did not have. ``f_N`` compiles from ``N_init`` alone and
is identical at both pitches, so the biomass *increment* is too, and the ratio
``(X₀ + ΔX)/(X₀ + N₀/f_N)`` is **structurally insensitive** to X₀ — not invariant, it drifts
toward 1 as X₀ grows, by a measured **0.0004** across a 6.25× change. Read against the run's own
pitch the ratio is 0.9848 at 0.04 against 0.9844 at 0.25. The 0.925 was the literal.

**The cell-nitrogen price was charged against the wrong kind of band.** ``biomass_N_fraction``'s
[0.08, 0.14] is an *end-state elemental composition*. Crépin's own end state falls **below the
whole band** — her measured 3.39 g/L final dry weight and her must's nitrogen give < 0.0786
(model frame) / < 0.0703 (paper frame), from ``X₀ < 0.83`` alone. A band whose low edge sits
above the measured whole-run value is not a ceiling on a transient store. And at the one moment
she actually measured, the landing capacity **undershoots** her cells rather than overloading
them: 0.092–0.097 g N/g DW where her landmark implies 0.108–0.156.

**What her data cannot do is settle the inoculum.** ``X₀ = X_final − YAN/f_N`` is a lever with the
wrong gain: across the whole candidate range 0 → 0.25 g/L the implied ``f_N`` moves only 8.0 %,
and separating 0.04 from 0.25 would need ``f_N`` known to 6.7 % where its own declared band spans
a factor of 1.75. The point estimate is frame-dependent too (0.168 g/L model frame, 0.507 paper
frame) and the archive records no assay frame for Coleman's own YAN, so neither end is usable.

**Two conditionals, stated rather than absorbed.** (1) The transient-density reading is an
*internal-consistency* argument on the model's own accounting convention —
``_apply_nitrogen_dependent_yield``'s "every assimilated gram of nitrogen enters biomass, so
Y_X/N = 1/f_N identically". Consumed is not the same as in-cells if yeast excrete nitrogen, and
nothing here sources that either way; this is not an independent source-side observation of
storage. (2) The moments are named, not matched: the model's total cell nitrogen peaks at
**20.5 %** (landing capacity) and **10.6 %** (sourced pitch) of the nitrogen gone, well before
Crépin's first landmark at 50 %, so the peak is compared to no measurement at all — which is
what D-251 §4 did.

**What this file forbids.** Re-charging the Coleman anchor against an inoculum change; re-pricing
a transient nitrogen store against ``biomass_N_fraction``'s end-state band; and reading "her data
cannot settle the inoculum" as "her data are silent" — they cap the end-state fraction below the
whole band, in both frames.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import compile_scenario
from tests.test_defined_media import (
    CREPIN_MUST_MM,
    _assimilable_n_mgl,
    commensurate_scenario,
    model_frame_mgn,
    paper_frame_mgn,
)
from tests.test_fusel_keto_acid_node import _OTHER_PRECURSOR_CONSUMERS
from tests.test_nitrogen_timing_attribution import (
    CREPIN_FINAL_BIOMASS_GPL,
    SHIPPED_PITCH_GPL,
    SOURCED_PITCH_GPL,
)

STORED = "stored_nitrogen"

#: Crépin Data Set S1 panel B, her FIRST sampled dry weight — at 16 h, her 50 %-nitrogen landmark.
#: Every bound below rests on this plus "no net biomass loss before the first sample", so the
#: inoculum cannot have exceeded it.
CREPIN_FIRST_DRY_WEIGHT_GPL = 0.83
CREPIN_FIRST_LANDMARK_FRACTION = 0.50

#: The capacity that lands her half-nitrogen landmark at each pitch (D-251). NOT shipped.
LANDING_R = {SHIPPED_PITCH_GPL: 3.9, SOURCED_PITCH_GPL: 2.65}


def _band() -> tuple[float, float]:
    """``biomass_N_fraction``'s declared band, read from the parameter file itself."""
    params = load_parameters(default_data_dir() / "wine_generic.yaml")
    unc = params["biomass_N_fraction"].uncertainty
    assert unc is not None, "biomass_N_fraction no longer declares a band"
    return float(unc.low), float(unc.high)


def _run(pitch: float, capacity_ratio: float) -> dict[str, Any]:
    """Crépin's must at one pitch and one capacity, with the courses this file reads."""
    scenario = commensurate_scenario("crepin")
    compiled = compile_scenario(
        scenario.model_copy(update={"initial": scenario.initial | {"pitch_gpl": pitch}})
    )
    for name in _OTHER_PRECURSOR_CONSUMERS:
        compiled.process_set.disable(name)  # KeyErrors on a rename rather than silently no-op
    values = dict(compiled.param_values)
    values["amino_acid_uptake_capacity_ratio"] = capacity_ratio
    traj = simulate_scheduled(
        compiled.process_set,
        values,
        compiled.y0,
        compiled.t_span_h,
        events=compiled.events,
        param_tiers=compiled.parameters.tier_map(),
    )
    assert traj.success, traj.message

    schema = compiled.schema
    biomass = traj.y[schema.slice("X"), :][0]
    store = traj.y[schema.slice(STORED), :][0]
    nitrogen = np.array([_assimilable_n_mgl(traj, schema, i) for i in range(traj.y.shape[1])])
    consumed = 1.0 - nitrogen / nitrogen[0]
    f_n = float(values["biomass_N_fraction"])
    # Total cell nitrogen per gram of ALL biomass -- the basis D-251 §4 priced on. The store is
    # the transient part; ``f_n`` is the structural part every gram carries by construction.
    cell_n = (store + biomass * f_n) / biomass
    half = int(np.nonzero(consumed >= 0.50)[0][0])
    top = int(np.argmax(cell_n))
    return {
        "pitch": pitch,
        "r": capacity_ratio,
        "biomass_N_fraction": f_n,
        "initial_n_gpl": nitrogen[0] / 1000.0,
        "peak_biomass": float(biomass.max()),
        "cell_n_max": float(cell_n[top]),
        "consumed_at_cell_n_max": float(consumed[top]),
        "cell_n_at_half": float(cell_n[half]),
    }


@pytest.fixture(scope="module")
def runs() -> dict[tuple[float, float], dict[str, Any]]:
    """The four runs this file reads — computed ONCE (the D-245 pattern)."""
    wanted = (
        (SHIPPED_PITCH_GPL, 1.0),
        (SOURCED_PITCH_GPL, 1.0),
        (SHIPPED_PITCH_GPL, LANDING_R[SHIPPED_PITCH_GPL]),
        (SOURCED_PITCH_GPL, LANDING_R[SOURCED_PITCH_GPL]),
    )
    return {key: _run(*key) for key in wanted}


def _crepin_cell_n(fraction: float, biomass_gpl: float, yan_gpl: float, x0: float, f_inoc: float):
    """Her cells' nitrogen per gram of ALL biomass at one landmark, on the model's own convention.

    ``(nitrogen she has consumed by then) + (what the inoculum itself brought)``, over the dry
    weight she measured. The inoculum's own nitrogen must be credited or the basis is the biomass
    *increment*, which is not what D-251 §4 priced.
    """
    return (fraction * yan_gpl + x0 * f_inoc) / biomass_gpl


# --------------------------------------------------------------------------------------------
# 1. The anchor D-249 charged
# --------------------------------------------------------------------------------------------


def test_the_coleman_anchor_is_STRUCTURALLY_INSENSITIVE_and_d249s_price_was_a_LITERAL(runs):
    """D-249 §3 priced the pitch move at "sells the anchor". It sells nothing (decision D-252).

    ``test_biomass_now_reaches_the_coleman_yield_the_compile_seam_installs`` compared peak biomass
    against ``0.25 + initial/f_N`` with the pitch written in as a literal. ``f_N`` is compiled from
    ``N_init`` alone — asserted below — so the increment is pitch-independent and both terms of the
    ratio shift by the same X₀. It is *structurally insensitive*, not invariant: the ratio drifts
    toward 1 as X₀ grows, and this pins the drift so "insensitive" can never be read as "constant".
    """
    shipped, sourced = runs[(SHIPPED_PITCH_GPL, 1.0)], runs[(SOURCED_PITCH_GPL, 1.0)]

    same_f_n = pytest.approx(sourced["biomass_N_fraction"], rel=1e-12)
    assert shipped["biomass_N_fraction"] == same_f_n, (
        "biomass_N_fraction now differs between the two pitches, so the compile seam has started "
        "reading something other than N_init; the whole insensitivity argument below rests on it "
        "reading N_init alone"
    )

    own, literal = {}, {}
    for pitch, run in ((SHIPPED_PITCH_GPL, shipped), (SOURCED_PITCH_GPL, sourced)):
        predicted_own = pitch + run["initial_n_gpl"] / run["biomass_N_fraction"]
        predicted_literal = SHIPPED_PITCH_GPL + run["initial_n_gpl"] / run["biomass_N_fraction"]
        own[pitch] = run["peak_biomass"] / predicted_own
        literal[pitch] = run["peak_biomass"] / predicted_literal

    assert own[SHIPPED_PITCH_GPL] == pytest.approx(0.9844, abs=0.002), (
        f"the anchor reads {own[SHIPPED_PITCH_GPL]:.4f} at the fixture's own pitch, not the "
        "0.9844 D-248 measured and this record re-derived"
    )
    assert own[SOURCED_PITCH_GPL] == pytest.approx(0.9848, abs=0.002), (
        f"the anchor reads {own[SOURCED_PITCH_GPL]:.4f} at the sourced pitch when the prediction "
        "uses that pitch; D-252 measured 0.9848"
    )

    drift = abs(own[SOURCED_PITCH_GPL] - own[SHIPPED_PITCH_GPL])
    assert drift == pytest.approx(0.0004, abs=0.0006), (
        f"the anchor moves {drift:.5f} across a 6.25x inoculum change, not the 0.0004 D-252 "
        "measured -- 'structurally insensitive' is a measured smallness, not an identity, and "
        "this arm is what stops it being quoted as one"
    )

    # The literal is what produced D-249 §3's price, and it is outside the guard's own +/-0.02.
    assert literal[SOURCED_PITCH_GPL] == pytest.approx(0.9252, abs=0.002), (
        f"the guard's hardcoded 0.25 now yields {literal[SOURCED_PITCH_GPL]:.4f} at pitch 0.04, "
        "not the 0.925 D-249 §3 published as the cost of moving the fixture"
    )
    assert abs(literal[SOURCED_PITCH_GPL] - 0.984) > 0.02 > abs(own[SOURCED_PITCH_GPL] - 0.984), (
        "the literal no longer breaks the guard's band where the run's own pitch keeps it inside; "
        "that contrast IS the correction and without it this file pins nothing"
    )


# --------------------------------------------------------------------------------------------
# 2. What her data DO settle
# --------------------------------------------------------------------------------------------


def test_crepins_own_data_cap_the_cell_nitrogen_fraction_BELOW_the_whole_declared_band():
    """Her end state sits under [0.08, 0.14] in BOTH frames (decision D-252).

    On the model's own convention that every assimilated gram of nitrogen enters biomass, her
    must's nitrogen and her measured final dry weight fix the whole-run yield up to the inoculum:
    ``f_N = YAN / (X_final − X₀)``. Her first sampled dry weight is 0.83 g/L and the curve rises,
    so ``X₀ < 0.83`` and ``f_N`` is capped. **The cap is below the band's own low edge**, which is
    what makes the band the wrong object to price a transient store against.

    Frame-robust by construction: the model frame counts every nitrogen atom, the papers' frame
    counts what yeast can release, and the cap holds under both.
    """
    low, high = _band()
    assert (low, high) == (0.08, 0.14), (
        f"biomass_N_fraction's declared band is now [{low}, {high}]; every claim in this file "
        "about the band being an end-state object is scored against [0.08, 0.14]"
    )

    caps = {}
    for frame, yan_mgl in (
        ("model", model_frame_mgn(CREPIN_MUST_MM)),
        ("paper", paper_frame_mgn(CREPIN_MUST_MM)),
    ):
        cap = (yan_mgl / 1000.0) / (CREPIN_FINAL_BIOMASS_GPL - CREPIN_FIRST_DRY_WEIGHT_GPL)
        caps[frame] = cap
        assert cap < low, (
            f"in the {frame} frame Crépin's own data cap the cell-nitrogen fraction at {cap:.4f}, "
            f"which is no longer below the declared band's low edge {low}; the argument that the "
            "band cannot be a ceiling on a transient store rests on exactly this"
        )

    assert caps["model"] == pytest.approx(0.0786, abs=0.001), caps
    assert caps["paper"] == pytest.approx(0.0703, abs=0.001), caps


def test_her_data_CANNOT_settle_the_inoculum_and_the_gain_required_says_why():
    """The mass balance is a lever with the wrong gain (decision D-252).

    ``X₀ = X_final − YAN/f_N`` is exact, and useless: the inoculum is a few per cent of the final
    biomass either way, so the whole candidate range barely moves the implied ``f_N``. This pins
    the *gain* rather than the conclusion, because "her data cannot settle it" is only worth
    recording alongside how far off the required precision is.
    """
    low, high = _band()
    yan_gpl = model_frame_mgn(CREPIN_MUST_MM) / 1000.0

    implied = {x0: yan_gpl / (CREPIN_FINAL_BIOMASS_GPL - x0) for x0 in (0.0, 0.04, 0.25)}
    spread = implied[0.25] / implied[0.0] - 1.0
    assert spread == pytest.approx(0.080, abs=0.005), (
        f"the implied cell-nitrogen fraction now moves {spread:.1%} across the whole candidate "
        "inoculum range 0 -> 0.25 g/L, not the 8.0 % D-252 measured"
    )

    separation = abs(implied[0.25] / implied[0.04] - 1.0)
    assert separation == pytest.approx(0.067, abs=0.005), (
        f"separating the sourced 0.04 from the fixture's 0.25 now needs f_N to {separation:.1%}, "
        "not the 6.7 % D-252 measured"
    )
    assert separation < (high / low - 1.0) / 5.0, (
        f"the precision the separation needs ({separation:.1%}) is no longer an order of magnitude "
        f"inside f_N's own declared band (a factor {high / low:.2f}); if the band ever narrows "
        "that far, her data DO start to discriminate and this null result reopens"
    )

    # The point estimate is frame-dependent, which is the other half of why it is not usable.
    point = {
        frame: CREPIN_FINAL_BIOMASS_GPL - (yan / 1000.0) / 0.06241
        for frame, yan in (
            ("model", model_frame_mgn(CREPIN_MUST_MM)),
            ("paper", paper_frame_mgn(CREPIN_MUST_MM)),
        )
    }
    assert point["model"] == pytest.approx(0.168, abs=0.01), point
    assert point["paper"] == pytest.approx(0.507, abs=0.01), point
    assert point["model"] < SHIPPED_PITCH_GPL < point["paper"], (
        f"the two frames' point estimates {point} no longer straddle the fixture's own pitch; "
        "the reason neither end is usable is that the frame decides which side of 0.25 it lands"
    )


# --------------------------------------------------------------------------------------------
# 3. The price D-251 charged
# --------------------------------------------------------------------------------------------


def test_at_HER_landmark_the_landing_capacity_UNDERSHOOTS_her_cell_nitrogen(runs):
    """D-251 §4 compared the model's peak against a band, and neither against her (decision D-252).

    Read at the one moment Crépin sampled — half the assimilable nitrogen gone — the model at the
    capacity that lands her biomass share carries **less** nitrogen per gram of cell than her own
    numbers imply, at every candidate inoculum and in both frames. The calibration does not
    overload the cells relative to the source; it still under-loads them.

    Her side is computed on the same total-cell basis as the model's, inoculum nitrogen credited,
    and it is an internal-consistency reading rather than a measurement of storage — see the
    module docstring's first conditional.
    """
    her = []
    for _frame, yan_mgl in (
        ("model", model_frame_mgn(CREPIN_MUST_MM)),
        ("paper", paper_frame_mgn(CREPIN_MUST_MM)),
    ):
        for x0 in (0.0, SOURCED_PITCH_GPL, SHIPPED_PITCH_GPL):
            for f_inoc in (0.06241, 0.114):
                her.append(
                    _crepin_cell_n(
                        CREPIN_FIRST_LANDMARK_FRACTION,
                        CREPIN_FIRST_DRY_WEIGHT_GPL,
                        yan_mgl / 1000.0,
                        x0,
                        f_inoc,
                    )
                )
    her_low, her_high = min(her), max(her)
    assert her_low == pytest.approx(0.108, abs=0.003), (
        f"the floor of Crépin's implied cell nitrogen at her 50 % landmark is now {her_low:.4f}, "
        "not the 0.108 D-252 measured across both frames, three inocula and two credits"
    )
    assert her_high == pytest.approx(0.156, abs=0.003), (
        f"its ceiling is now {her_high:.4f}, not the 0.156 D-252 measured"
    )

    for pitch, expected in ((SHIPPED_PITCH_GPL, 0.0920), (SOURCED_PITCH_GPL, 0.0971)):
        run = runs[(pitch, LANDING_R[pitch])]
        at_half = run["cell_n_at_half"]
        assert at_half == pytest.approx(expected, abs=0.003), (
            f"at pitch {pitch} and the landing capacity the model holds {at_half:.4f} g N/g DW "
            f"when half the nitrogen is gone, not the {expected} D-252 measured"
        )
        assert at_half < her_low, (
            f"the model now holds {at_half:.4f} at her own landmark, no longer below the "
            f"{her_low:.4f} floor of what her data imply -- D-251 §4's 'the calibration overloads "
            "the cells' would start to be true and its refusal would need re-deriving"
        )


def test_the_models_cell_nitrogen_PEAK_is_at_a_moment_crepin_never_sampled(runs):
    """The number D-251 §4 priced is a peak, and it falls before her first sample (decision D-252).

    This is the second conditional made checkable: 0.1391 is not "the model at her landmark", it
    is the model's own transient maximum, at 20.5 % of the nitrogen gone where her first landmark
    is at 50 %. A price charged at a moment with no measurement in it is a price against nothing.
    """
    for pitch, top, consumed in (
        (SHIPPED_PITCH_GPL, 0.1391, 0.205),
        (SOURCED_PITCH_GPL, 0.1372, 0.106),
    ):
        run = runs[(pitch, LANDING_R[pitch])]
        assert run["cell_n_max"] == pytest.approx(top, abs=0.002), (
            f"at pitch {pitch} the landing capacity's peak cell nitrogen is "
            f"{run['cell_n_max']:.4f}, not the {top} D-251 §4 priced"
        )
        assert run["consumed_at_cell_n_max"] == pytest.approx(consumed, abs=0.03), (
            f"that peak now falls at {run['consumed_at_cell_n_max']:.1%} of the nitrogen gone, "
            f"not the {consumed:.1%} D-252 measured"
        )
        assert run["consumed_at_cell_n_max"] < CREPIN_FIRST_LANDMARK_FRACTION, (
            "the peak has moved to or past Crépin's first landmark, so it is no longer a moment "
            "she never sampled and the second conditional in this file's docstring is stale"
        )


def test_the_shipped_capacity_is_the_CONTROL_and_undershoots_her_further(runs):
    """Without the control the undershoot could be a property of the observable (decision D-252).

    At the shipped capacity the model carries less again — so the gap against her landmark is a
    thing the capacity moves, in the direction the calibration moves it, and the comparison in
    :func:`test_at_HER_landmark_the_landing_capacity_UNDERSHOOTS_her_cell_nitrogen` is not
    measuring something inert.
    """
    for pitch in (SHIPPED_PITCH_GPL, SOURCED_PITCH_GPL):
        shipped = runs[(pitch, 1.0)]
        landing = runs[(pitch, LANDING_R[pitch])]
        assert shipped["cell_n_at_half"] < landing["cell_n_at_half"], (
            f"at pitch {pitch} the shipped capacity now holds at least as much cell nitrogen at "
            f"her landmark as the landing capacity ({shipped['cell_n_at_half']:.4f} vs "
            f"{landing['cell_n_at_half']:.4f}); the capacity no longer moves this observable and "
            "the undershoot reading needs re-deriving"
        )
        assert shipped["cell_n_max"] < landing["cell_n_max"], (
            f"at pitch {pitch} the shipped capacity's peak cell nitrogen no longer sits below the "
            f"landing capacity's -- {shipped['cell_n_max']:.4f} vs {landing['cell_n_max']:.4f}"
        )

    assert runs[(SHIPPED_PITCH_GPL, 1.0)]["cell_n_max"] == pytest.approx(0.0833, abs=0.002), (
        "the shipped capacity's peak cell nitrogen is no longer the 0.0833 D-251 §4 measured"
    )
