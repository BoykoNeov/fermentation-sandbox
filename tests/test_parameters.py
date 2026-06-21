"""Tests for the provenance-enforced parameter store."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import StateSchema, VarSpec
from fermentation.core.tiers import Tier
from fermentation.parameters import Parameter, ParameterSet, Uncertainty
from fermentation.parameters.store import load_parameters

DATA = Path(__file__).resolve().parents[1] / "src" / "fermentation" / "parameters" / "data"


def good_param(**overrides: object) -> dict[str, object]:
    base = {
        "name": "mu_max",
        "value": 0.33,
        "unit": "1/h",
        "tier": "plausible",
        "uncertainty": {"low": 0.2, "high": 0.5},
        "provenance": {"source": "Coleman et al. 2007", "conditions": "wine must, 20 C"},
    }
    base.update(overrides)
    return base


def test_valid_parameter_constructs():
    p = Parameter(**good_param())
    assert p.value == 0.33
    assert p.tier is Tier.PLAUSIBLE
    assert p.provenance.source.startswith("Coleman")


def test_tier_accepts_label_or_int():
    assert Parameter(**good_param(tier="speculative")).tier is Tier.SPECULATIVE
    assert Parameter(**good_param(tier=2)).tier is Tier.VALIDATED


def test_unknown_tier_label_rejected():
    with pytest.raises(ValidationError, match="unknown tier"):
        Parameter(**good_param(tier="guesstimate"))


def test_missing_provenance_rejected():
    bad = good_param()
    del bad["provenance"]
    with pytest.raises(ValidationError):
        Parameter(**bad)


def test_missing_uncertainty_rejected():
    bad = good_param()
    del bad["uncertainty"]
    with pytest.raises(ValidationError):
        Parameter(**bad)


def test_empty_source_rejected():
    with pytest.raises(ValidationError):
        Parameter(**good_param(provenance={"source": "  ", "conditions": "x"}))


def test_value_outside_uncertainty_rejected():
    with pytest.raises(ValidationError, match="outside"):
        Parameter(**good_param(value=99.0))


def test_uncertainty_low_above_high_rejected():
    with pytest.raises(ValidationError, match="low.*high|high"):
        Uncertainty(low=5.0, high=1.0)


def test_extra_fields_forbidden():
    with pytest.raises(ValidationError):
        Parameter(**good_param(typo_field=1))


def test_parameterset_resolve_and_tier():
    p1 = Parameter(**good_param(name="a", tier="validated"))
    p2 = Parameter(**good_param(name="b", tier="speculative"))
    ps = ParameterSet([p1, p2])
    assert ps.resolve() == {"a": 0.33, "b": 0.33}
    assert ps.resolve(["a"]) == {"a": 0.33}
    assert ps.lowest_tier() is Tier.SPECULATIVE
    assert ps.lowest_tier(["a"]) is Tier.VALIDATED


def test_parameterset_tier_map_bridges_to_process_tier_propagation():
    """The production seam D-1 closes: a parameter's YAML tier ->
    ``ParameterSet.tier_map()`` -> ``ProcessSet.tier_of`` caps a validated process.
    Every other propagation test hand-builds the ``{name: Tier}`` dict; this one
    runs it through the real ``ParameterSet`` so the bridge itself is covered."""
    params = ParameterSet(
        [
            Parameter(**good_param(name="k_rate", tier="speculative")),
            Parameter(**good_param(name="k_other", tier="validated")),
        ]
    )
    assert params.tier_map() == {"k_rate": Tier.SPECULATIVE, "k_other": Tier.VALIDATED}

    class Reader(Process):
        name = "reader"
        tier = Tier.VALIDATED
        touches = ("S",)
        reads = ("k_rate",)

        def derivatives(self, t, y, schema, params):
            return schema.zeros()

    pset = ProcessSet(StateSchema([VarSpec("S", "g/L")]), [Reader()])
    # Structural-only over-reports; the real param tier map drags it to speculative.
    assert pset.tier_of("S") is Tier.VALIDATED
    assert pset.tier_of("S", params.tier_map()) is Tier.SPECULATIVE


def test_parameterset_duplicate_rejected():
    p = Parameter(**good_param(name="a"))
    with pytest.raises(ValueError, match="Duplicate"):
        ParameterSet([p, p])


def test_parameterset_merge_override():
    base = ParameterSet([Parameter(**good_param(name="a", value=0.3))])
    overlay = ParameterSet([Parameter(**good_param(name="a", value=0.4))])
    with pytest.raises(ValueError, match="both"):
        base.merge(overlay)
    merged = base.merge(overlay, override=True)
    assert merged.value("a") == 0.4


def test_missing_parameter_lookup_message():
    ps = ParameterSet([Parameter(**good_param(name="a"))])
    with pytest.raises(KeyError, match="No parameter"):
        _ = ps["zzz"]


def test_load_shipped_wine_parameters():
    ps = load_parameters(DATA / "wine_generic.yaml")
    assert "mu_max" in ps
    assert "Y_ethanol_sugar" in ps
    # The realised ethanol yield is literature-grounded.
    assert ps.tier_of("Y_ethanol_sugar") is Tier.PLAUSIBLE
    # mu_max is now sourced (Coleman 2007), so it is promoted out of speculative.
    assert ps.tier_of("mu_max") is Tier.PLAUSIBLE
    # K_s has no literature analogue (Coleman growth is Monod-on-N only) - still speculative.
    assert ps.tier_of("K_s") is Tier.SPECULATIVE
    assert 0.46 <= ps.value("Y_ethanol_sugar") <= 0.48


def test_load_shipped_beer_parameters():
    ps = load_parameters(DATA / "beer_generic.yaml")
    # Defines every parameter the medium-agnostic kinetics read.
    read_params = {
        "mu_max", "K_s", "K_n", "q_sugar_max", "K_sugar_uptake", "K_repression",
        "ethanol_tolerance", "ethanol_inhibition_exponent", "E_a_growth",
        "E_a_uptake", "T_ref", "biomass_C_fraction", "biomass_N_fraction",
    }
    assert read_params <= set(ps.names)
    # Sourced from Zamudio Lara et al. 2022 (open access).
    assert ps.tier_of("mu_max") is Tier.PLAUSIBLE
    assert ps.value("mu_max") == pytest.approx(0.098)
    assert ps.tier_of("K_sugar_uptake") is Tier.PLAUSIBLE
    # Honestly thinner than wine: transferred/derived values stay speculative.
    assert ps.tier_of("K_n") is Tier.SPECULATIVE  # transferred from the wine fit
    assert ps.tier_of("q_sugar_max") is Tier.SPECULATIVE  # derived, growth-coupled origin
    assert ps.tier_of("E_a_growth") is Tier.SPECULATIVE  # de Andres-Toro primary not read


def test_load_rejects_param_missing_provenance(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "naked_constant:\n  value: 1.0\n  unit: '1/h'\n  tier: plausible\n"
        "  uncertainty: {low: 0.0, high: 2.0}\n",
        encoding="utf-8",
    )
    with pytest.raises(ValidationError):
        load_parameters(bad)


def test_load_rejects_key_name_mismatch(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "alpha:\n  name: beta\n  value: 1.0\n  unit: '1/h'\n  tier: plausible\n"
        "  uncertainty: {low: 0, high: 2}\n"
        "  provenance: {source: x, conditions: y}\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="disagrees"):
        load_parameters(bad)
