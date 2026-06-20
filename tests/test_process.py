"""Tests for the Process abstraction and ProcessSet aggregation/tiering."""

from collections.abc import Mapping

import pytest

from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema, VarSpec
from fermentation.core.tiers import Tier


def schema() -> StateSchema:
    return StateSchema([VarSpec("S", "g/L"), VarSpec("E", "g/L"), VarSpec("X", "g/L")])


class ToyFermentation(Process):
    """Consumes sugar, makes ethanol at a fixed yield — a stand-in for kinetics."""

    name = "toy_fermentation"
    tier = Tier.VALIDATED
    touches = ("S", "E")

    def __init__(self, rate: float = 1.0):
        self.rate = rate

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        s = schema.get(y, "S")
        yield_ = params["Y_ethanol_sugar"]
        consume = self.rate if s > 0 else 0.0
        d[schema.slice("S")] = -consume
        d[schema.slice("E")] = consume * yield_
        return d


class SpeculativeAging(Process):
    name = "spec_aging"
    tier = Tier.SPECULATIVE
    touches = ("E",)

    def derivatives(self, t, y, schema, params):
        d = schema.zeros()
        d[schema.slice("E")] = -1e-3
        return d


class Leaky(Process):
    """Declares it only touches X but actually writes to S — contract violation."""

    name = "leaky"
    tier = Tier.PLAUSIBLE
    touches = ("X",)

    def derivatives(self, t, y, schema, params):
        d = schema.zeros()
        d[schema.slice("S")] = 1.0  # not declared!
        return d


def test_total_derivatives_sums_processes():
    s = schema()
    ps = ProcessSet(s, [ToyFermentation(rate=2.0)])
    y = s.pack({"S": 100.0, "E": 0.0, "X": 1.0})
    d = ps.total_derivatives(0.0, y, {"Y_ethanol_sugar": 0.47})
    assert s.get(d, "S") == pytest.approx(-2.0)
    assert s.get(d, "E") == pytest.approx(2.0 * 0.47)
    assert s.get(d, "X") == 0.0


def test_disable_process_drops_its_contribution():
    s = schema()
    ps = ProcessSet(s, [ToyFermentation(), SpeculativeAging()])
    y = s.pack({"S": 100.0, "E": 0.0, "X": 1.0})
    params = {"Y_ethanol_sugar": 0.47}

    with_spec = ps.total_derivatives(0.0, y, params).copy()
    ps.disable("spec_aging")
    without_spec = ps.total_derivatives(0.0, y, params)
    # The speculative process removed exactly its -1e-3 ethanol drain.
    assert s.get(with_spec, "E") == pytest.approx(s.get(without_spec, "E") - 1e-3)


def test_tier_of_variable_is_lowest_contributor():
    s = schema()
    ps = ProcessSet(s, [ToyFermentation(), SpeculativeAging()])
    # E is touched by a VALIDATED and a SPECULATIVE process -> speculative wins.
    assert ps.tier_of("E") is Tier.SPECULATIVE
    # S only by the validated process.
    assert ps.tier_of("S") is Tier.VALIDATED
    # X untouched -> trivially validated.
    assert ps.tier_of("X") is Tier.VALIDATED


def test_tier_map_and_overall():
    s = schema()
    ps = ProcessSet(s, [ToyFermentation(), SpeculativeAging()])
    assert ps.tier_map() == {"S": Tier.VALIDATED, "E": Tier.SPECULATIVE, "X": Tier.VALIDATED}
    assert ps.overall_tier() is Tier.SPECULATIVE
    ps.disable("spec_aging")
    assert ps.overall_tier() is Tier.VALIDATED
    assert ps.tier_of("E") is Tier.VALIDATED


def test_strict_mode_catches_touch_contract_violation():
    s = schema()
    ps = ProcessSet(s, [Leaky()], strict=True)
    y = s.pack({"S": 1.0, "E": 0.0, "X": 1.0})
    with pytest.raises(ValueError, match="undeclared variables"):
        ps.total_derivatives(0.0, y, {})


def test_non_strict_mode_does_not_check():
    s = schema()
    ps = ProcessSet(s, [Leaky()], strict=False)
    y = s.pack({"S": 1.0, "E": 0.0, "X": 1.0})
    # No exception; the leak is summed (this is the speed/safety tradeoff).
    d = ps.total_derivatives(0.0, y, {})
    assert s.get(d, "S") == 1.0


def test_duplicate_process_names_rejected():
    s = schema()
    with pytest.raises(ValueError, match="Duplicate process names"):
        ProcessSet(s, [ToyFermentation(), ToyFermentation()])


def test_process_touching_unknown_variable_rejected():
    s = schema()

    class Bad(Process):
        name = "bad"
        tier = Tier.PLAUSIBLE
        touches = ("nonexistent",)

        def derivatives(self, t, y, schema, params):
            return schema.zeros()

    with pytest.raises(ValueError, match="unknown variables"):
        ProcessSet(s, [Bad()])


def test_enable_disable_unknown_raises():
    s = schema()
    ps = ProcessSet(s, [ToyFermentation()])
    with pytest.raises(KeyError):
        ps.disable("nope")
