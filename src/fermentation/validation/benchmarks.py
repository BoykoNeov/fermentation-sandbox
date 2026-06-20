"""Benchmark specifications (handoff section 2.2) and a measured-series comparator.

The §2.2 acceptance criteria for Milestone 1 are encoded here as data, before the
model is tuned, so the model is test-driven against them. Each :class:`BenchmarkSpec`
states what to measure, the acceptable window, and the literature basis.

:class:`ReferenceSeries` + :func:`compare_series` are the seam for *real* datasets:
when measured Brix/gravity/temperature time-series become available they drop in
here and the same comparator (RMSE/MAE) scores the model against them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class BenchmarkSpec:
    """A declarative acceptance criterion for the validated core.

    ``metric`` names the observable (e.g. ``"days_to_dryness"``); ``low``/``high``
    bound the acceptable value; ``conditions`` and ``source`` record the setup and
    its basis. Whether a run *passes* is evaluated by Milestone 1 code once the
    kinetics exist — the spec itself is stable data.
    """

    key: str
    description: str
    metric: str
    low: float
    high: float
    unit: str
    conditions: str
    source: str

    def passes(self, value: float) -> bool:
        return self.low <= value <= self.high


@dataclass(frozen=True)
class ReferenceSeries:
    """A measured time-series to validate against (real data, when we have it).

    ``time_h`` and ``value`` are equal-length arrays in canonical internal units.
    ``tier`` is informational text (e.g. 'measured', 'digitized from figure').
    """

    name: str
    time_h: NDArray[np.float64]
    value: NDArray[np.float64]
    unit: str
    source: str
    tier: str = "measured"
    meta: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        t = np.asarray(self.time_h, dtype=np.float64)
        v = np.asarray(self.value, dtype=np.float64)
        if t.shape != v.shape:
            raise ValueError(
                f"ReferenceSeries {self.name!r}: time {t.shape} and value {v.shape} "
                "must have equal shape"
            )
        if t.ndim != 1:
            raise ValueError(f"ReferenceSeries {self.name!r}: expected 1-D arrays")
        object.__setattr__(self, "time_h", t)
        object.__setattr__(self, "value", v)


@dataclass(frozen=True)
class FitResult:
    rmse: float
    mae: float
    n: int


def compare_series(
    model_time_h: NDArray[np.float64],
    model_value: NDArray[np.float64],
    reference: ReferenceSeries,
) -> FitResult:
    """Score a model trajectory against a reference series.

    The model is linearly interpolated onto the reference's time points (within
    the model's covered span) and RMSE/MAE are computed there. Reference points
    outside the model's time span are ignored. Raises if there is no overlap.
    """
    mt = np.asarray(model_time_h, dtype=np.float64)
    mv = np.asarray(model_value, dtype=np.float64)
    if mt.ndim != 1 or mt.shape != mv.shape:
        raise ValueError("model_time_h and model_value must be equal-length 1-D arrays")

    in_span = (reference.time_h >= mt.min()) & (reference.time_h <= mt.max())
    if not in_span.any():
        raise ValueError(
            f"No overlap between model span [{mt.min():.3g}, {mt.max():.3g}] h "
            f"and reference {reference.name!r}"
        )
    ref_t = reference.time_h[in_span]
    ref_v = reference.value[in_span]
    model_at_ref = np.interp(ref_t, mt, mv)
    residual = model_at_ref - ref_v
    rmse = float(np.sqrt(np.mean(residual**2)))
    mae = float(np.mean(np.abs(residual)))
    return FitResult(rmse=rmse, mae=mae, n=int(in_span.sum()))


# -- The Milestone 1 acceptance criteria (handoff section 2.2) ----------------

BENCHMARKS: dict[str, BenchmarkSpec] = {
    "wine_dryness": BenchmarkSpec(
        key="wine_dryness",
        description=(
            "A ~24 Brix must (~264 g/L sugar) at 20 C ferments to dryness with a "
            "visible lag -> exponential -> stationary biomass trajectory."
        ),
        metric="days_to_dryness",
        low=10.0,
        high=14.0,
        unit="day",
        conditions="~264 g/L initial sugar, 20 C, single wine strain, nitrogen-limited",
        source="handoff section 2.2 (wine benchmark)",
    ),
    "beer_attenuation": BenchmarkSpec(
        key="beer_attenuation",
        description="A ~1.048 OG ale wort at 20 C attenuates to roughly 1.010.",
        metric="days_to_target_gravity",
        low=5.0,
        high=7.0,
        unit="day",
        conditions="OG 1.048 -> ~1.010, 20 C, single ale strain",
        source="handoff section 2.2 (beer benchmark)",
    ),
    "co2_peak_then_tail": BenchmarkSpec(
        key="co2_peak_then_tail",
        description=(
            "CO2 evolution rate rises to a peak then tails off; its integral "
            "tracks sugar consumed (primary measurable validation channel)."
        ),
        metric="co2_integral_vs_sugar_consumed_ratio",
        low=0.95,
        high=1.05,
        unit="dimensionless",
        conditions="carbon balance between evolved CO2 and sugar consumed",
        source="handoff sections 1.2, 2.2",
    ),
}
"""Acceptance criteria for the validated core. Evaluation lands in Milestone 1;
the specs are stable data now so the model can be test-driven against them."""
