"""End-to-end: the architecture must be runnable before any real kinetics exist."""

import numpy as np

from fermentation.core.process import ProcessSet
from fermentation.runtime import simulate
from fermentation.validation import assert_conserved, assert_nonnegative, max_drift


def test_simulate_runs_and_conserves_mass(toy_schema, toy_process):
    ps = ProcessSet(toy_schema, [toy_process], strict=True)
    y0 = toy_schema.pack({"S": 200.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 200.0))

    assert traj.success
    # Sugar consumed, ethanol and CO2 produced.
    assert traj.series("S")[-1] < 1.0
    assert traj.series("E")[-1] > 90.0
    assert traj.series("CO2")[-1] > 90.0

    # Total mass S+E+CO2 is conserved by construction.
    def total(y):
        return float(sum(y[toy_schema.slice(v)][0] for v in ("S", "E", "CO2")))

    assert max_drift(traj, total) < 1e-3
    assert_conserved(traj, total, rtol=1e-5, atol=1e-6, label="total mass")
    assert_nonnegative(traj, ("S", "E", "CO2"))


def test_simulate_respects_schema_size(toy_schema, toy_process):
    ps = ProcessSet(toy_schema, [toy_process])
    bad_y0 = np.zeros(99)
    try:
        simulate(ps, params={}, y0=bad_y0, t_span=(0.0, 1.0))
    except ValueError as e:
        assert "shape" in str(e)
    else:  # pragma: no cover
        raise AssertionError("expected ValueError on wrong y0 shape")


def test_trajectory_tier_map_propagates(toy_schema, toy_process):
    ps = ProcessSet(toy_schema, [toy_process])
    y0 = toy_schema.pack({"S": 50.0, "E": 0.0, "CO2": 0.0})
    traj = simulate(ps, params={}, y0=y0, t_span=(0.0, 50.0))
    assert traj.tier_map["S"].label == "validated"
    assert traj.overall_tier().label == "validated"
