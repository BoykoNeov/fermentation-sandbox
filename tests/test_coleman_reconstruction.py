"""Faithful reconstruction of Coleman, Fish & Block (2007) as a tracked regression.

This is the evidence base for decision **D-14**, which overturned the earlier
"uptake-speed gap" reading. An independent re-implementation of Coleman's own
comprehensive model (their eqs 1-8 with the Table A2 regressions evaluated at
20 C, which the paper validates against the measured Fig 6c sugar curves) is
integrated with pure scipy -- no engine code -- to generate a reference sugar
trajectory. ``compare_series`` then scores our engine against it through the
standard ``ReferenceSeries`` seam. The claim it guards: the validated core
reproduces Coleman line-for-line across the nitrogen range (low-N 80 mg/L and
normal-N 330 mg/L), so a wine that ferments ~9 d at 20 C is *correct*, not too
fast -- the 10-14 d expectation was a generic heuristic (see docs/DECISIONS.md).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation.benchmarks import ReferenceSeries, compare_series

_T_C = 20.0  # the wine benchmark temperature; Coleman's regressions are read here


def _coleman_params(n0_mgL: float) -> dict[str, float]:
    """Coleman Table A2 regressions. ``Y_X/N`` regresses on INITIAL nitrogen
    (mg/L); the rest on temperature (deg C). The ``k'_d`` and ``Y_X/N`` a1
    exponents are typo-corrected (decisions D-13, D-14)."""
    return {
        "mu_max": math.exp(-3.92 + 0.0782 * _T_C),
        "K_N": math.exp(-4.73),
        "beta_max": math.exp(-2.30 + 0.0771 * _T_C),
        "Y_ES": math.exp(-0.598),
        "K_S": math.exp(2.33),
        "k_prime_d": math.exp(-9.81 - 0.108 * _T_C + 0.00478 * _T_C**2),
        "Y_XN": math.exp(3.50 - 3.61e-3 * n0_mgL),
    }


def _coleman_reference(
    n0_mgL: float, *, s0: float = 264.0, pitch: float = 0.25, hours: float = 720.0
) -> ReferenceSeries:
    """Integrate Coleman's eqs 1-8 and return its sugar(t) as a ``ReferenceSeries``."""
    p = _coleman_params(n0_mgL)

    def rhs(t: float, y: np.ndarray) -> list[float]:
        x, x_a, n, e, s = y
        mu = p["mu_max"] * n / (p["K_N"] + n) if n > 0.0 else 0.0
        beta = p["beta_max"] * s / (p["K_S"] + s) if s > 0.0 else 0.0
        k_d = p["k_prime_d"] * e
        return [
            mu * x_a,
            mu * x_a - k_d * x_a,
            -mu * x_a / p["Y_XN"],
            beta * x_a,
            -beta * x_a / p["Y_ES"],
        ]

    t_eval = np.linspace(0.0, hours, int(hours // 6) + 1)
    sol = solve_ivp(
        rhs,
        (0.0, hours),
        [pitch, pitch, n0_mgL / 1000.0, 0.0, s0],
        method="BDF",
        t_eval=t_eval,
        rtol=1e-9,
        atol=1e-11,
    )
    assert sol.success, sol.message
    return ReferenceSeries(
        name=f"coleman-sugar-{n0_mgL:.0f}mgN",
        time_h=sol.t,
        value=sol.y[4],
        unit="g/L",
        source="Coleman, Fish & Block 2007 comprehensive model (eqs 1-8, Table A2 @ 20 C)",
        tier="reconstructed from published model",
    )


def _our_sugar(n0_mgL: float, *, pitch: float = 0.25) -> tuple[np.ndarray, np.ndarray]:
    scenario = Scenario(
        name="coleman-recon",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": n0_mgL, "pitch_gpl": pitch},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=_T_C)],
        duration_days=30.0,
    )
    compiled = compile_scenario(scenario, strict=True)
    t_eval = np.linspace(0.0, 30.0 * 24.0, 4001)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success
    return traj.t, np.asarray(traj.series("S"))


@pytest.mark.parametrize("n0_mgL", [80.0, 330.0])
def test_engine_reproduces_coleman_sugar_curve(n0_mgL: float) -> None:
    # Line-for-line agreement on the sugar trajectory at both the low- and
    # normal-nitrogen ends. This is a model-vs-model comparison of the *kinetics*,
    # so both sides must start from the same initial sugar. Our engine now loads the
    # must_fermentable_fraction-corrected sugar (~245 g/L at 24 Brix, decision D-16),
    # so we feed that same S0 to the Coleman reference rather than its raw ~264 g/L
    # default; the glycerol diversion leaves dS untouched, so the curves still track.
    # Observed RMSE is ~1.3 g/L (~0.5 % of S0) at both N levels; the 2.0 g/L
    # threshold is observed + ~50 % margin (NOT a loose pass).
    t_h, sugar = _our_sugar(n0_mgL)
    ref = _coleman_reference(n0_mgL, s0=float(sugar[0]))
    fit = compare_series(t_h, sugar, ref)
    assert fit.rmse < 2.0, f"RMSE {fit.rmse:.2f} g/L vs Coleman at {n0_mgL:.0f} mg N/L"
    # Both finish (the low-N case is the one that used to stick before D-14).
    assert float(sugar[-1]) < 2.0


def test_low_nitrogen_completes_after_dynamic_yield() -> None:
    # Regression guard for the D-14 fix proper: at 80 mg N/L the fixed-yield model
    # built too little biomass and STUCK; the nitrogen-dependent yield lets it
    # finish, matching Coleman.
    _, sugar = _our_sugar(80.0)
    assert float(sugar[-1]) < 1.0
