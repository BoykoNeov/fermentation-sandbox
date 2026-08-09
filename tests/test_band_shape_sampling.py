"""The band-SHAPE assumption in the sampler, pinned (decision D-165).

``ensemble.py``'s module docstring claimed that reporting outer percentiles
"de-sensitises the result to the shape choice". D-165 measured that claim. It is true
across the two shapes the module *offers* — both linear in the parameter — and false for
the shape that 84 of the archive's wide bands are actually built to:

=========================  ==================  ==================  ===============
shape (band [1e-4, 1e-3, 1e-2], r = 100)       P5 / nominal        P95 / nominal
=========================  ==================  ==================  ===============
triangular (the default)                       0.7675x             7.8893x
uniform (also offered)                         0.5950x             9.5050x
log-triangular (NOT offered)                   0.2071x             4.8281x
=========================  ==================  ==================  ===============

The two offered shapes agree to within 1.30x at both ends; the log-scale shape sits
3.71x from the triangular P5. So the de-sensitising claim was true as far as it was ever
exercised and silent about the case the wide bands need. The docstring now says so, and
these tests are what keep it honest.

**Why this matters beyond a docstring.** 123 of 337 live bands span ``hi/lo >= 10``, and
84 of those put the nominal at the *geometric* centre — for which a linear triangular has
mean ``(lo+m+hi)/3`` and is therefore dominated by ``hi``. The worst,
``k_d2_ethanol_tolerance_death`` (``hi/lo = 300``), draws a median 9.18x its own stated
nominal, in every wine and beer ensemble.

**What is deliberately NOT pinned here: a fix.** No shape field is shipped. Which bands
were stated multiplicatively is not in the schema — :class:`Uncertainty` carries endpoints
only — so selecting a shape per parameter would be a magic constant in a new place
(prime directive 2). D-165 measured the log-triangular candidate and recorded it as
flagged, not adopted, exactly as D-164 §6 did with the admissible-range field.

**Non-vacuity.** ``test_reference_band_is_the_shipped_one`` reads the band from the real
parameter store and asserts it equals the ``_LO``/``_VAL``/``_HI`` these tests compute
from. Without it every other test here would derive its expectation from its own
constants and pass against any YAML whatsoever — the D-108/D-109 self-derivation shape.
Each test below was checked against its own mutation arm (D-155 standard); the arms are
named in D-165 §5.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from fermentation.parameters.store import ParameterSet, default_data_dir, load_parameters
from fermentation.runtime.ensemble import _inverse_cdf, sample_parameters, simulate_ensemble
from fermentation.scenario.compile import compile_scenario
from fermentation.scenario.schema import Scenario, TemperaturePoint

#: The reference band: ``k_autolysis``, log-symmetric over a full two decades. PRINTED
#: from ``wine_generic.yaml``, and pinned against it by the first test below.
_LO, _VAL, _HI = 1.0e-4, 1.0e-3, 1.0e-2

#: How close the two OFFERED (linear) shapes come to each other at the outer percentiles.
#: This is the quantity the docstring's "de-sensitises" claim asserts is small.
_LINEAR_PAIR_AGREE_WITHIN = 1.30

#: How far the NOT-offered log-scale shape sits from the default at P5. The claim's gap.
_LOG_SHAPE_DIVERGES_BY_AT_LEAST = 3.5


def _wine() -> ParameterSet:
    """The band under test lives in wine_generic.yaml."""
    return load_parameters(default_data_dir() / "wine_generic.yaml")


def _log_triangular_quantile(q: float, lo: float, val: float, hi: float) -> float:
    """The shape the module does not offer: triangular in log-space."""
    return math.exp(_inverse_cdf(q, math.log(lo), math.log(val), math.log(hi), "triangular"))


def test_reference_band_is_the_shipped_one():
    """Non-vacuity: the constants above ARE ``k_autolysis``'s band, not this file's idea
    of it. If the YAML band moves, every other test in this file is measuring a band the
    archive no longer has, and must be re-derived rather than re-pinned."""
    p = _wine()["k_autolysis"]
    assert (p.uncertainty.low, p.value, p.uncertainty.high) == (_LO, _VAL, _HI)


def test_triangular_mean_is_displaced_by_the_band_width():
    """The mechanism, in closed form: a linear triangular's mean is ``(lo+m+hi)/3``, so a
    band whose ``hi`` is 10x the nominal is mean-dominated by ``hi``. 3.70x, not ~1x."""
    assert pytest.approx(3.70, abs=5e-3) == (_LO + _VAL + _HI) / (3 * _VAL)


def test_the_displacement_survives_in_the_statistic_actually_reported():
    """``Band`` reports the MEDIAN, never the mean — so a displacement that lived only in
    the mean would be a fair objection to this whole record. It does not: the median is
    3.33x the nominal. Mutation arm: assert on the mean here and the test still passes,
    which is why it is written against the median."""
    median = _inverse_cdf(0.5, _LO, _VAL, _HI, "triangular")
    assert median / _VAL == pytest.approx(3.3254, rel=1e-3)


def test_the_two_offered_shapes_do_agree_at_the_outer_percentiles():
    """The docstring's claim, in the scope where it holds. Both offered shapes are linear
    in the parameter, and at the outer percentiles they agree to within 1.30x."""
    for q in (0.05, 0.95):
        tri = _inverse_cdf(q, _LO, _VAL, _HI, "triangular")
        uni = _inverse_cdf(q, _LO, _VAL, _HI, "uniform")
        ratio = max(tri, uni) / min(tri, uni)
        assert ratio < _LINEAR_PAIR_AGREE_WITHIN, (
            f"the two OFFERED shapes disagree by {ratio:.3f}x at q={q}; the docstring's "
            "de-sensitising claim is scoped to exactly this pair and would need rewording"
        )


def test_but_they_do_not_de_sensitise_against_a_log_scale_shape():
    """The claim's gap, and the reason the docstring was corrected. A log-triangular over
    the same band puts P5 at 0.207x where the linear pair put it at 0.60-0.77x. Outer
    percentiles do not rescue a linear shape assumption on a band stated over decades."""
    tri_p5 = _inverse_cdf(0.05, _LO, _VAL, _HI, "triangular")
    log_p5 = _log_triangular_quantile(0.05, _LO, _VAL, _HI)
    assert tri_p5 / log_p5 > _LOG_SHAPE_DIVERGES_BY_AT_LEAST
    # ... and the log-scale median lands exactly on the nominal, because this band is
    # log-symmetric. That is an ALGEBRAIC IDENTITY, not evidence the shape is right --
    # it is asserted here only so the contrast above cannot be read as a shape endorsement.
    assert _log_triangular_quantile(0.5, _LO, _VAL, _HI) == pytest.approx(_VAL, rel=1e-12)


def test_a_narrow_band_is_not_displaced_at_all():
    """DESIGNED NULL. Without this the file would be measuring "sampling moves things"
    rather than "band WIDTH moves things". A +/-20% band displaces the mean by 1.1%."""
    lo, val, hi = 1.0 / 1.2, 1.0, 1.2
    assert (lo + val + hi) / (3 * val) == pytest.approx(1.0111, abs=5e-4)
    assert _inverse_cdf(0.5, lo, val, hi, "triangular") / val == pytest.approx(1.0, abs=0.01)


def test_the_worst_displaced_band_is_inert_at_24_brix_and_by_how_little():
    """The reason D-165 is a *sampling* finding and not a trajectory catastrophe.

    ``k_d2_ethanol_tolerance_death`` is the archive's worst-displaced band (r=300, median
    9.18x nominal) and ``_active_reads`` says it is drawn in every wine ensemble. It is
    nonetheless **exactly inert** on a normal must: ``EthanolToleranceDeath`` multiplies it
    by ``max(E - ethanol_tolerance, 0)**2``, and a 24-Brix wine peaks at ~117 g/L ethanol.
    So a ``reads``-based reach count is an UPPER BOUND on consequential reach, and any
    audit that treats "drawn" as "in play" overstates its own surface.

    The margin is the point. ``ethanol_tolerance`` is itself sampled, over [120, 150] --
    so the comparison is between two draws, and the peak across the whole ensemble came in
    at 119.4 g/L against a band low of 120.0. **0.5%.** This test fails if that margin
    closes, because a speculative death term switching on for part of an ensemble is not
    something that should happen silently. It is a guard on the gap, not on the gate."""
    scen = Scenario(
        name="d165-margin",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25},
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=14.0, celsius=25.0),
        ],
        interventions=[],
        duration_days=30.0,
    )
    c = compile_scenario(scen)
    t_end = 30.0 * 24.0
    ens = simulate_ensemble(
        c.process_set,
        c.parameters,
        c.y0,
        (0.0, t_end),
        n_members=16,
        seed=0,
        t_eval=np.linspace(0.0, t_end, 60),
        sampler="lhs",
    )
    tol_low = c.parameters["ethanol_tolerance"].uncertainty.low
    peak = float(np.max(ens.members[:, c.schema.slice("E"), :]))
    assert peak < tol_low, (
        f"peak ensemble ethanol {peak:.4g} g/L has reached the ethanol_tolerance band low "
        f"{tol_low:.4g} g/L -- EthanolToleranceDeath is no longer inert for every member of "
        "a standard 24-Brix ensemble, so k_d2's 300-fold band is now live on a normal must. "
        "Re-measure D-165's reach numbers before relaxing this."
    )
    # The band low is the binding edge, not the nominal 142 -- pin that it is what we read.
    assert tol_low == 120.0 and c.parameters.value("ethanol_tolerance") == 142.0


def test_the_sampler_really_draws_this_way_not_just_the_inverse_cdf():
    """The tests above exercise ``_inverse_cdf`` (the LHS/Sobol path). Pin that the plain
    Monte-Carlo path -- ``sample_parameters``, the default -- shows the same displacement,
    so neither path is quietly exempt."""
    pset = _wine()
    rng = np.random.default_rng(0)
    draws = np.array(
        [sample_parameters(pset, rng, names=["k_autolysis"])["k_autolysis"] for _ in range(20_000)]
    )
    assert draws.mean() / _VAL == pytest.approx(3.70, rel=0.02)
    assert np.median(draws) / _VAL == pytest.approx(3.33, rel=0.03)
