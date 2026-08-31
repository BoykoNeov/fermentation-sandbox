"""D-249's front-loading residue, read in the frame Crépin actually sampled (decision D-251).

D-249 §5 measured that the model needs **54.3 %** of its peak biomass to have eaten half the
must's assimilable nitrogen where Crépin's yeast need **24.5 %**, and attributed it: "``ρ_N ∝ X``
is the wrong functional form early". Two things about that measurement are corrected here, and
neither needs a line of ``src/``.

**The frame.** ``_assimilable_n_mgl`` counts the D-250 intracellular store as unconsumed — on
purpose, because D-248's "40.8 % → 0.62 % residual" was measured when the surplus sat in ``N`` and
dropping the store would redefine that number rather than migrate it. Crépin sampled the
**medium**. Read on the narrower extracellular quantity D-250 made separable for the first time,
the same shipped run holds **41.2 %** and **62.6 %** at her two landmarks rather than 54.3 % and
77.7 % — two fifths of the front-loading gap is the frame, not the model.

**The reachability.** D-249 §4 swept ``amino_acid_uptake_capacity_ratio`` across three orders of
magnitude and found the exhaustion *time* saturating, which is true and is why the timing miss is
refused. It never scored the sweep on the **shape** observable. In the total frame that sweep is
flat to three decimals (0.543 → 0.542 across 1000×), which is exactly why nothing was seen; in
Crépin's frame the same knob walks the half-nitrogen share 0.412 → 0.091 and crosses her 0.245
at ``r ≈ 3.9``, **inside the parameter's own declared [0.5, 10] band**. It lands her
three-quarters landmark at the same capacity, and the landmark *ratio* — which one knob is not
free to choose independently — lands with it (1.648 against her 1.673).

So the form is not shown to be wrong: the level was never identifiable on any observable this repo
could express before D-250 split the two frames. What is **not** settled is whether the capacity
that lands is affordable. Front-loading is paid for in stored nitrogen, and at the landing
capacity the cells carry a total of **0.133–0.139 g N per g dry weight** against the repo's own
sourced N-replete elemental reference of **0.114** (band [0.08, 0.14]) — it spends essentially the
whole band, at both pitches. That, and the fact that the landing capacity is not pitch-invariant
(3.9 at the fixture's unsourced 0.25 g/L, ~2.65 at the only sourced pitch), is why this file
measures the calibration and does **not** take it. Moving ``r`` is the owner's trade.

**What this file forbids.** Re-proposing "the uptake rate law's functional form is wrong early" as
a settled finding; re-running the capacity sweep against the timing (D-249 §4 did it, it
saturates, it is refused); and taking the front-loading calibration without pricing the cell
nitrogen it costs.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from fermentation.core.state import FloatArray
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import compile_scenario
from tests.test_defined_media import (
    HOUSE_PITCH_GPL,
    SOURCED_PITCH_GPL,
    _assimilable_n_mgl,
    commensurate_scenario,
)
from tests.test_fusel_keto_acid_node import _OTHER_PRECURSOR_CONSUMERS
from tests.test_nitrogen_timing_attribution import (
    CREPIN_DRY_WEIGHT_GPL,
    CREPIN_FINAL_BIOMASS_GPL,
    _crepin_run,
)

STORED = "stored_nitrogen"

#: Crépin's own dry weights as a share of her final biomass, at the nitrogen landmarks her Data
#: Set S1 panel B reports them at. The first two are the discriminating pair; the third is
#: matched by this model at **every** capacity including the shipped one, so it separates nothing
#: — see :func:`test_the_third_landmark_discriminates_NOTHING_and_is_labelled_so`.
CREPIN_SHARE = {f: CREPIN_DRY_WEIGHT_GPL[f] / CREPIN_FINAL_BIOMASS_GPL for f in (0.50, 0.75, 1.00)}

#: The capacity that lands Crépin's half-nitrogen landmark, per pitch. NOT shipped — the shipped
#: value is 1.0 and this file does not move it. Measured on the grid in :func:`sweep`.
LANDING_R = {HOUSE_PITCH_GPL: 3.9, SOURCED_PITCH_GPL: 2.65}

#: The grid each pitch is read on. 1000 is D-249 §4's own saturation point, kept so the
#: total-frame flatness claim is made against the same extreme that record used.
GRID = {HOUSE_PITCH_GPL: (1.0, 2.0, 3.9, 5.0, 1000.0), SOURCED_PITCH_GPL: (1.0, 2.65)}


def _cross(t: FloatArray, series: FloatArray, level: float) -> float:
    """Interpolated first crossing of a RISING series — D-249's idiom, restricted.

    That file's ``_Course._cross`` also handles a falling series (it reads sugar); this one
    only ever sees a consumed-fraction, so the falling branch is dropped rather than carried
    untested. The two agree on the rising case, which is the one both files compare hours on.
    """
    hit = np.nonzero(series >= level)[0]
    assert hit.size, f"the run never reaches {level}"
    i = int(hit[0])
    if i == 0:
        return float(t[0])
    return float(np.interp(level, [series[i - 1], series[i]], [t[i - 1], t[i]]))


def _read(pitch: float, r: float) -> dict[str, Any]:
    """Both nitrogen frames off one trajectory, plus the store expressed as a per-cell quota."""
    traj, schema = _crepin_run(pitch, capacity_ratio=r)
    t = traj.t
    biomass = traj.y[schema.slice("X"), :][0]
    total = np.array([_assimilable_n_mgl(traj, schema, i) for i in range(traj.y.shape[1])])
    store = np.maximum(traj.y[schema.slice(STORED), :][0], 0.0) * 1000.0  # mg N/L
    peak = float(biomass.max())

    frames = {
        "total": 1.0 - total / total[0],
        # What a sample of Crépin's medium would have contained: the store is inside the
        # cells. ``store[0]`` is zero on every scenario this repo compiles — nothing seeds a
        # pre-loaded store — but it is subtracted explicitly so the denominator stays the
        # medium's own initial nitrogen if one ever is, rather than changing meaning silently.
        "medium": 1.0 - (total - store) / (total[0] - store[0]),
    }
    shares = {
        frame: {
            f: float(np.interp(_cross(t, consumed, f), t, biomass)) / peak
            for f in (0.50, 0.75, 0.99)
        }
        for frame, consumed in frames.items()
    }
    # g stored N per g dry cell. Below 0.05 g/L of cells the ratio is the pitch instant's
    # rounding error rather than a quota, so it is excluded rather than allowed to set the max.
    live = biomass > 0.05
    quota = store / np.maximum(biomass, 1e-12) / 1000.0
    return {
        "shares": shares,
        "quota_max": float(quota[live].max()),
        "store_share": float(store.max() / total[0]),
        "peak_biomass": peak,
    }


@pytest.fixture(scope="module")
def sweep() -> dict[tuple[float, float], dict[str, Any]]:
    """Seven integrations, paid once (the D-245 pattern)."""
    return {(pitch, r): _read(pitch, r) for pitch, grid in GRID.items() for r in grid}


@pytest.fixture(scope="module")
def wine_params():
    return load_parameters(default_data_dir() / "wine_generic.yaml")


def test_read_in_the_frame_crepin_SAMPLED_two_fifths_of_the_front_loading_gap_is_the_FRAME(sweep):
    """D-249's 2.2× front-loading gap is 1.7× in the medium, with nothing in ``src/`` changed.

    Crépin's Data Set S1 measures what is left in the must. The model's own quantity counts the
    D-250 store as still-unconsumed, so its nitrogen clock runs late against hers and the biomass
    that has accumulated by the crossing is correspondingly larger. Both readings come off the
    **same trajectory** here — only the observable differs — so this is a frame statement and
    cannot be a model one.

    The total-frame numbers this asserts are D-249's own pinned 0.543/0.777, re-derived rather
    than quoted, so the two files cannot drift apart on the baseline.
    """
    shares = sweep[(HOUSE_PITCH_GPL, 1.0)]["shares"]

    assert shares["total"][0.50] == pytest.approx(0.543, abs=0.01)
    assert shares["total"][0.75] == pytest.approx(0.777, abs=0.01)
    assert shares["medium"][0.50] == pytest.approx(0.412, abs=0.01)
    assert shares["medium"][0.75] == pytest.approx(0.626, abs=0.01)

    for landmark, expected_total, expected_medium in ((0.50, 2.22, 1.68), (0.75, 1.90, 1.53)):
        measured = CREPIN_SHARE[landmark]
        gap_total = shares["total"][landmark] / measured
        gap_medium = shares["medium"][landmark] / measured
        assert gap_total == pytest.approx(expected_total, abs=0.05)
        assert gap_medium == pytest.approx(expected_medium, abs=0.05), (
            f"{landmark:.0%} of the nitrogen gone: the medium-frame gap reads {gap_medium:.2f}× "
            f"against D-251's {expected_medium}× — the frame correction has moved"
        )
        assert gap_medium < gap_total, (
            "the medium frame must read a SMALLER gap than the total frame; if it does not, the "
            "store is empty or is being counted on the wrong side and this whole file is void"
        )

    assert sweep[(HOUSE_PITCH_GPL, 1.0)]["store_share"] > 0.10, (
        "the two frames can only differ while the store holds something — at 16.3 % of the "
        "must's nitrogen it does, and this is the non-vacuity arm for every comparison above"
    )


def test_in_the_TOTAL_frame_the_capacity_knob_is_INVISIBLE_which_is_why_d249_saw_nothing(sweep):
    """The mechanism of the earlier miss, pinned so it cannot be re-made (decision D-251).

    D-249 §4's sweep was scored on exhaustion *time*, and §5's shape table was computed in the
    total frame at ``r = 1`` only. Had the sweep been scored on the shape in that frame it would
    still have found nothing: across a **1000×** capacity change the total-frame shares move in
    the third decimal. The knob moves nitrogen between the medium and the store, and the total
    frame is blind to that transfer by construction — it counts both sides.
    """
    lo = sweep[(HOUSE_PITCH_GPL, 1.0)]["shares"]["total"]
    hi = sweep[(HOUSE_PITCH_GPL, 1000.0)]["shares"]["total"]

    for landmark in (0.50, 0.75):
        assert abs(hi[landmark] - lo[landmark]) < 0.005, (
            f"{landmark:.0%}: the total-frame share moved {abs(hi[landmark] - lo[landmark]):.4f} "
            "across 1000× of capacity. It used to be flat, and that flatness is the stated "
            "reason D-249's sweep could not see the knob"
        )

    medium_lo = sweep[(HOUSE_PITCH_GPL, 1.0)]["shares"]["medium"][0.50]
    medium_hi = sweep[(HOUSE_PITCH_GPL, 1000.0)]["shares"]["medium"][0.50]
    assert medium_lo - medium_hi > 0.30, (
        "the same 1000× in the MEDIUM frame must move the half-nitrogen share by a lot "
        f"(0.412 → 0.091 when measured); it moved {medium_lo - medium_hi:.3f}. Without this arm "
        "the flatness above would be evidence of a dead knob rather than of a blind observable"
    )


def test_the_shipped_knob_REACHES_crepins_landmark_inside_its_own_declared_band(sweep, wine_params):
    """D-249 §5's "wrong functional form" does not survive the frame correction (decision D-251).

    In Crépin's frame the half-nitrogen share falls monotonically with capacity and crosses her
    0.245 at ``r ≈ 3.9`` — inside the ``[0.5, 10]`` uncertainty the parameter already declares, so
    reaching her landmark needs no new parameter, no widened band and no new rate law. What this
    does **not** establish is that the crossing capacity is affordable; that is the next test.
    """
    band = wine_params["amino_acid_uptake_capacity_ratio"].uncertainty
    assert band is not None and (band.low, band.high) == (0.5, 10.0), (
        "the declared band moved; the claim below is that the crossing lies INSIDE it, so the "
        "band's edges are load-bearing literals here"
    )

    for pitch in (HOUSE_PITCH_GPL, SOURCED_PITCH_GPL):
        grid = sorted(GRID[pitch])
        shares = [sweep[(pitch, r)]["shares"]["medium"][0.50] for r in grid]
        assert all(b < a for a, b in zip(shares, shares[1:], strict=False)), (
            f"pitch {pitch}: the half-nitrogen share must fall monotonically with capacity; "
            f"it read {shares} across {grid}"
        )
        assert shares[0] > CREPIN_SHARE[0.50] > shares[-1], (
            f"pitch {pitch}: Crépin's {CREPIN_SHARE[0.50]:.3f} is no longer BRACKETED by the "
            f"grid ({shares[0]:.3f} … {shares[-1]:.3f}), so the knob no longer crosses her "
            "landmark and D-249 §5's unreachability verdict would be back"
        )
        landing = LANDING_R[pitch]
        assert band.low <= landing <= band.high, (
            f"pitch {pitch}: the landing capacity {landing} left the parameter's declared band"
        )
        assert sweep[(pitch, landing)]["shares"]["medium"][0.50] == pytest.approx(
            CREPIN_SHARE[0.50], abs=0.01
        ), (
            f"pitch {pitch}: capacity {landing} no longer lands the half-nitrogen landmark, so "
            "LANDING_R is stale and every price quoted against it is quoted at the wrong point"
        )


def test_one_capacity_lands_BOTH_discriminating_landmarks_and_their_RATIO(sweep):
    """Two landmarks and one knob — which is what separates this from fitting a number to itself.

    Crépin's three-quarters landmark is not free once the half is matched: a single capacity sets
    both, and their ratio (1.673 measured) is a statement about the *shape* of her nitrogen-versus-
    biomass curve that one knob cannot tune independently. At the landing capacity the model reads
    1.648 at the fixture's pitch and 1.701 at the sourced one, against her 1.673 — so the shape
    lands at both, with different capacities.

    This is the anti-D-98 arm: the record is entitled to say the form reproduces her landmarks
    only while the *second* landmark and the ratio come along uninvited.
    """
    for pitch, expected in ((HOUSE_PITCH_GPL, 1.648), (SOURCED_PITCH_GPL, 1.701)):
        shares = sweep[(pitch, LANDING_R[pitch])]["shares"]["medium"]

        assert shares[0.75] == pytest.approx(CREPIN_SHARE[0.75], abs=0.015), (
            f"pitch {pitch}: with the half-nitrogen landmark landed the three-quarters one reads "
            f"{shares[0.75]:.3f} against Crépin's {CREPIN_SHARE[0.75]:.3f}. It is not free — if "
            "it stops coming along, one knob is fitting one number and the record's shape claim "
            "is void"
        )
        ratio = shares[0.75] / shares[0.50]
        assert ratio == pytest.approx(expected, abs=0.02)
        assert ratio == pytest.approx(CREPIN_SHARE[0.75] / CREPIN_SHARE[0.50], abs=0.05), (
            f"pitch {pitch}: the landmark ratio reads {ratio:.3f} against Crépin's "
            f"{CREPIN_SHARE[0.75] / CREPIN_SHARE[0.50]:.3f}"
        )


def test_the_third_landmark_discriminates_NOTHING_and_is_labelled_so(sweep):
    """Crépin's N_T landmark is matched at every capacity, the shipped one included.

    Her yeast hold 99.1 % of final biomass when the must's nitrogen is gone; so does this model,
    at ``r = 1`` and at ``r = 1000`` alike. That is structural — near-complete consumption and
    near-peak biomass are the same event in a nitrogen-limited run — so the landmark carries no
    information about the capacity and must never be counted as a third point of agreement.
    """
    for (pitch, r), reading in sweep.items():
        near_complete = reading["shares"]["medium"][0.99]
        assert near_complete == pytest.approx(CREPIN_SHARE[1.00], abs=0.005), (
            f"pitch {pitch}, r={r}: the near-complete landmark reads {near_complete:.3f} against "
            f"Crépin's {CREPIN_SHARE[1.00]:.3f}"
        )

    spread = max(reading["shares"]["medium"][0.99] for reading in sweep.values()) - min(
        reading["shares"]["medium"][0.99] for reading in sweep.values()
    )
    assert spread < 0.01, (
        f"the near-complete landmark spread {spread:.4f} across the whole grid. If it ever "
        "becomes capacity-sensitive it turns into a real third constraint — a finding to record, "
        "not a bound to relax"
    )


def test_the_front_loading_is_PAID_FOR_in_stored_nitrogen_at_the_cell_N_ceiling(sweep, wine_params):
    """The price of landing Crépin's landmarks, against the repo's OWN sourced ceiling (D-251).

    Nitrogen taken up ahead of anabolic need has to sit somewhere, and D-250's store is where.
    ``biomass_N_fraction``'s shipped 0.114 g N/g dry weight is the canonical N-**replete** yeast
    elemental reference (Roels 1983 / Heijnen, CH1.8O0.5N0.2) — *total* cell nitrogen, band
    [0.08, 0.14]. The wine seam overrides it downward per Coleman's ``Y_X/N(N_init)``, to 0.0624
    on this must, so the headroom a stored pool may occupy is what is left below the ceiling.

    At the shipped capacity the store is comfortably inside it. At the landing capacity the cells
    carry 0.133–0.139 g N/g DW: past the reference value by 17–22 % and within a hair of the
    band's top edge, at **both** pitches. The calibration is therefore not free, and this test
    exists so that a beat which takes it has to say so.
    """
    ceiling = wine_params["biomass_N_fraction"]
    assert ceiling.value == pytest.approx(0.114) and ceiling.uncertainty is not None
    top = ceiling.uncertainty.high
    assert top == pytest.approx(0.14), "the ceiling's high edge is a load-bearing literal here"

    structural = float(
        compile_scenario(commensurate_scenario("crepin")).param_values["biomass_N_fraction"]
    )
    assert structural == pytest.approx(0.0624, abs=0.001), (
        "the compiled structural nitrogen fraction moved, so every headroom figure below is "
        "measured against a different cell"
    )

    shipped = structural + sweep[(HOUSE_PITCH_GPL, 1.0)]["quota_max"]
    assert shipped == pytest.approx(0.083, abs=0.005)
    assert shipped < ceiling.value, (
        "at the SHIPPED capacity the cells must stay below the N-replete reference — this is the "
        "control that makes the ceiling pressure below a property of the calibration"
    )

    for pitch, expected in ((HOUSE_PITCH_GPL, 0.1391), (SOURCED_PITCH_GPL, 0.1371)):
        reading = sweep[(pitch, LANDING_R[pitch])]
        loaded = structural + reading["quota_max"]
        assert loaded == pytest.approx(expected, abs=0.005), (
            f"pitch {pitch}: at the landing capacity the cells carry {loaded:.4f} g N/g DW, not "
            f"the {expected} D-251 priced"
        )
        assert loaded > ceiling.value, (
            f"pitch {pitch}: landing Crépin's landmarks no longer pushes the cells past the "
            f"N-replete reference ({loaded:.4f} ≤ {ceiling.value}). That would make the "
            "calibration cheap, and it is the reason D-251 left it to the owner — re-price it "
            "rather than deleting this"
        )
        assert loaded <= top, (
            f"pitch {pitch}: {loaded:.4f} g N/g DW is outside the sourced band's own high edge "
            f"{top}, which would make the calibration inadmissible rather than expensive"
        )
        assert reading["store_share"] > 0.30, (
            f"pitch {pitch}: the store holds {reading['store_share']:.1%} of the must's nitrogen "
            "at its peak — the transient the price above is charged for"
        )


def test_the_stores_nonnegativity_dip_is_SOLVER_NOISE_and_is_there_at_the_shipped_capacity():
    """Why the −2.3e-9 dip in ``stored_nitrogen`` is not an over-draw (decision D-251).

    Growth's draw is split between ``N`` and the store in proportion to what each holds, so the
    store's own share of the draw goes to zero as the store empties and no over-draw is possible
    by construction. The falsifiable consequence is that the dip must scale with **solver
    tolerance** and not with the uptake capacity. It does: at the shipped capacity on this must
    the store reaches −2.3e-9 at ``rtol=1e-6`` and −4e-20 at ``rtol=1e-10``, a collapse of eleven
    orders for four of tolerance.

    It is recorded because :func:`~fermentation.validation.assert_nonnegative`'s 1e-9 default sits
    *below* that excursion: the guard is already within a factor of two of firing on shipped
    parameters, and no shipped test happens to integrate a trajectory where it does. That is a
    tolerance calibration, not a physics defect, and the arms below are what license saying so.
    """
    minima = {}
    for rtol, atol in ((1e-6, 1e-9), (1e-10, 1e-13)):
        compiled = compile_scenario(commensurate_scenario("crepin"))
        for name in _OTHER_PRECURSOR_CONSUMERS:
            compiled.process_set.disable(name)
        traj = simulate_scheduled(
            compiled.process_set,
            compiled.param_values,
            compiled.y0,
            compiled.t_span_h,
            events=compiled.events,
            param_tiers=compiled.parameters.tier_map(),
            rtol=rtol,
            atol=atol,
        )
        assert traj.success, traj.message
        minima[rtol] = float(traj.y[compiled.schema.slice(STORED), :][0].min())

    assert -1e-8 < minima[1e-6] < -1e-9, (
        f"the default-tolerance dip reads {minima[1e-6]:.3e}; D-251 measured −2.3e-9, i.e. past "
        "assert_nonnegative's 1e-9 default. If it is now inside that default the guard's "
        "brittleness has gone away and the record should say so"
    )
    assert abs(minima[1e-10]) < abs(minima[1e-6]) / 1e4, (
        f"tightening the solver took the dip from {minima[1e-6]:.3e} only to "
        f"{minima[1e-10]:.3e}. A dip that does NOT collapse with tolerance is an over-draw in "
        "the proportional split — a physics defect, and a different finding from this one"
    )
