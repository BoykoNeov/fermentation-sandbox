"""Tests for the Process abstraction and ProcessSet aggregation/tiering."""

from collections.abc import Mapping

import pytest

from fermentation.core.process import Process, ProcessSet, RateModifier
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


# -- rate modifiers (the multiplicative path) --------------------------------


class HalfRate(RateModifier):
    """Scales its target Process's whole contribution by 0.5."""

    name = "half_rate"
    tier = Tier.PLAUSIBLE
    modifies = ("toy_fermentation",)

    def factor(self, t, y, schema, params):
        return 0.5


class SpecDamper(RateModifier):
    """A speculative half-rate modifier, for tier-propagation tests."""

    name = "spec_damper"
    tier = Tier.SPECULATIVE
    modifies = ("toy_fermentation",)

    def factor(self, t, y, schema, params):
        return 0.5


def test_modifier_scales_only_its_target_process():
    s = schema()
    # ToyFermentation (S, E) is halved; SpeculativeAging (E only) is untouched.
    ps = ProcessSet(s, [ToyFermentation(rate=2.0), SpeculativeAging()], modifiers=[HalfRate()])
    y = s.pack({"S": 100.0, "E": 0.0, "X": 1.0})
    d = ps.total_derivatives(0.0, y, {"Y_ethanol_sugar": 0.47})
    # Sugar consumption halved: -2.0 -> -1.0.
    assert s.get(d, "S") == pytest.approx(-1.0)
    # Ethanol = halved toy yield (2*0.47*0.5) + the *unscaled* aging drain (-1e-3).
    assert s.get(d, "E") == pytest.approx(2.0 * 0.47 * 0.5 - 1e-3)


def test_modifiers_compose_multiplicatively():
    s = schema()
    ps = ProcessSet(s, [ToyFermentation(rate=2.0)], modifiers=[HalfRate(), SpecDamper()])
    y = s.pack({"S": 100.0, "E": 0.0, "X": 1.0})
    d = ps.total_derivatives(0.0, y, {"Y_ethanol_sugar": 0.47})
    # 0.5 * 0.5 = 0.25 applied to -2.0.
    assert s.get(d, "S") == pytest.approx(-0.5)


def test_disabled_modifier_contributes_unit_factor():
    s = schema()
    ps = ProcessSet(s, [ToyFermentation(rate=2.0)], modifiers=[HalfRate()])
    y = s.pack({"S": 100.0, "E": 0.0, "X": 1.0})
    params = {"Y_ethanol_sugar": 0.47}
    ps.disable("half_rate")
    d = ps.total_derivatives(0.0, y, params)
    assert s.get(d, "S") == pytest.approx(-2.0)  # full, unscaled rate restored


def test_modifier_drags_target_variable_tier_down():
    s = schema()
    # ToyFermentation is VALIDATED; a SPECULATIVE modifier scaling it makes the
    # variables it touches (S, E) report speculative.
    ps = ProcessSet(s, [ToyFermentation()], modifiers=[SpecDamper()])
    assert ps.tier_of("S") is Tier.SPECULATIVE
    assert ps.tier_of("E") is Tier.SPECULATIVE
    assert ps.overall_tier() is Tier.SPECULATIVE
    # X is untouched by the process, so the modifier cannot reach it.
    assert ps.tier_of("X") is Tier.VALIDATED
    # Disabling the modifier restores the process's own validated tier.
    ps.disable("spec_damper")
    assert ps.tier_of("S") is Tier.VALIDATED
    assert ps.overall_tier() is Tier.VALIDATED


def test_modifier_preserves_strict_touch_contract():
    s = schema()
    # Scaling a contribution by a scalar keeps undeclared slots at zero, so strict
    # mode still passes for a well-behaved process under a modifier.
    ps = ProcessSet(s, [ToyFermentation(rate=2.0)], modifiers=[HalfRate()], strict=True)
    y = s.pack({"S": 100.0, "E": 0.0, "X": 1.0})
    d = ps.total_derivatives(0.0, y, {"Y_ethanol_sugar": 0.47})
    assert s.get(d, "X") == 0.0


def test_modifier_name_clashing_with_process_rejected():
    s = schema()

    class Clash(RateModifier):
        name = "toy_fermentation"  # same as the process
        tier = Tier.PLAUSIBLE
        modifies = ("toy_fermentation",)

        def factor(self, t, y, schema, params):
            return 1.0

    with pytest.raises(ValueError, match="clash with process names"):
        ProcessSet(s, [ToyFermentation()], modifiers=[Clash()])


def test_modifier_targeting_unknown_process_rejected():
    s = schema()

    class BadTarget(RateModifier):
        name = "bad_target"
        tier = Tier.PLAUSIBLE
        modifies = ("nonexistent_process",)

        def factor(self, t, y, schema, params):
            return 1.0

    with pytest.raises(ValueError, match="unknown process"):
        ProcessSet(s, [ToyFermentation()], modifiers=[BadTarget()])


def test_duplicate_modifier_names_rejected():
    s = schema()
    with pytest.raises(ValueError, match="Duplicate modifier names"):
        ProcessSet(s, [ToyFermentation()], modifiers=[HalfRate(), HalfRate()])


# -- parameter-tier propagation (D-1) ----------------------------------------


class ParamReadingProcess(Process):
    """A VALIDATED process that reads one parameter — exercises D-1 tier capping."""

    name = "param_reading"
    tier = Tier.VALIDATED
    touches = ("S", "E")
    reads = ("k_rate",)

    def derivatives(self, t, y, schema, params):
        return schema.zeros()


class ParamReadingModifier(RateModifier):
    """A VALIDATED modifier reading one parameter, scaling ``param_reading``."""

    name = "param_reading_mod"
    tier = Tier.VALIDATED
    modifies = ("param_reading",)
    reads = ("k_mod",)

    def factor(self, t, y, schema, params):
        return 1.0


def test_speculative_param_caps_validated_process_output():
    s = schema()
    ps = ProcessSet(s, [ParamReadingProcess()])
    # The process is VALIDATED, but it reads a speculative parameter, so the
    # variables it touches must report speculative — no credibility borrowing (D-1).
    spec = {"k_rate": Tier.SPECULATIVE}
    assert ps.tier_of("S", spec) is Tier.SPECULATIVE
    assert ps.tier_of("E", spec) is Tier.SPECULATIVE
    assert ps.overall_tier(spec) is Tier.SPECULATIVE
    assert ps.tier_map(spec) == {"S": Tier.SPECULATIVE, "E": Tier.SPECULATIVE, "X": Tier.VALIDATED}
    # A variable the process does not touch is unaffected by its reads.
    assert ps.tier_of("X", spec) is Tier.VALIDATED


def test_validated_param_leaves_process_tier_unchanged():
    s = schema()
    ps = ProcessSet(s, [ParamReadingProcess()])
    # Flip the same parameter to VALIDATED and the output stays VALIDATED — the cap
    # is exactly the lowest input tier, nothing lower invented.
    val = {"k_rate": Tier.VALIDATED}
    assert ps.tier_of("S", val) is Tier.VALIDATED
    assert ps.overall_tier(val) is Tier.VALIDATED


def test_param_tiers_omitted_is_structural_only():
    s = schema()
    ps = ProcessSet(s, [ParamReadingProcess()])
    # Without param_tiers the result is the process's own (validated) tier — the
    # narrower, pre-D-1 structural answer that ignores its reads.
    assert ps.tier_of("S") is Tier.VALIDATED
    assert ps.tier_map() == {"S": Tier.VALIDATED, "E": Tier.VALIDATED, "X": Tier.VALIDATED}


def test_speculative_param_read_by_modifier_drags_target_down():
    s = schema()
    ps = ProcessSet(s, [ParamReadingProcess()], modifiers=[ParamReadingModifier()])
    # Both process and modifier are VALIDATED, but the modifier reads a speculative
    # parameter -> every variable the modified process touches goes speculative.
    tiers = {"k_rate": Tier.VALIDATED, "k_mod": Tier.SPECULATIVE}
    assert ps.tier_of("S", tiers) is Tier.SPECULATIVE
    assert ps.tier_of("E", tiers) is Tier.SPECULATIVE
    # Disabling the modifier removes its speculative read from the variable's tier.
    ps.disable("param_reading_mod")
    assert ps.tier_of("S", tiers) is Tier.VALIDATED


def test_declared_read_missing_from_param_tiers_raises():
    s = schema()
    ps = ProcessSet(s, [ParamReadingProcess()])
    # k_rate is a declared read but absent from the supplied map: fail loudly rather
    # than silently treat an unknown-provenance input as validated.
    with pytest.raises(KeyError, match="k_rate"):
        ps.tier_of("S", {"unrelated": Tier.VALIDATED})
