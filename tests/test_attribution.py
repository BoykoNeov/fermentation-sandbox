"""Ensemble spread attribution — first-order sensitivity (analysis layer, D-24 follow-up).

The wrapper (D-24) reports an ensemble's *spread*; this layer answers *which parameters
drive it* and *how that spread partitions across confidence tiers*. It is a post-hoc
standardized-regression (SRC) decomposition over one ensemble's stored draws — no extra
integrations. These tests pin the contract: the SRC² budget plus the unexplained
(nonlinear/interaction) remainder accounts for the spread, a dominant parameter is
identified with the right sign, widening a band shifts its share, pinned parameters are
excluded, and the degenerate/ill-posed cases are handled honestly rather than silently.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.analysis import attribute_spread
from fermentation.core.chemistry import CO2_PER_HEXOSE, ETHANOL_PER_HEXOSE
from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema, VarSpec
from fermentation.core.tiers import Tier
from fermentation.parameters.schema import Parameter, Provenance, Uncertainty
from fermentation.parameters.store import ParameterSet
from fermentation.runtime import simulate_ensemble

# -- a toy whose output is (early) near-linear in a sampled parameter ----------
#
# Monod uptake: over an early window (before S depletes) the ethanol produced is
# essentially linear in ``vmax`` and nearly flat in ``ks``, so SRC should assign almost
# all of the variance to ``vmax`` with R^2 ~ 1 — a clean anchor for the budget contract.


class ParamFerm(Process):
    name = "pf"
    tier = Tier.VALIDATED
    touches = ("S", "E", "CO2")
    reads = ("vmax", "ks")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        s = float(y[schema.slice("S")][0])
        if s <= 0.0:
            return d
        r = params["vmax"] * s / (params["ks"] + s)
        d[schema.slice("S")] = -r
        d[schema.slice("E")] = r * ETHANOL_PER_HEXOSE
        d[schema.slice("CO2")] = r * CO2_PER_HEXOSE
        return d


def _param(name, value, low, high, tier=Tier.PLAUSIBLE) -> Parameter:
    return Parameter(
        name=name,
        value=value,
        unit="1/h",
        tier=tier,
        uncertainty=Uncertainty(low=low, high=high),
        provenance=Provenance(source="author estimate", conditions="toy"),
    )


def _schema() -> StateSchema:
    return StateSchema([VarSpec("S", "g/L"), VarSpec("E", "g/L"), VarSpec("CO2", "g/L")])


def _pset(
    *,
    vmax_band: tuple[float, float] = (3.0, 7.0),
    ks_band: tuple[float, float] = (4.0, 6.0),
) -> ParameterSet:
    # ks is tagged SPECULATIVE, vmax PLAUSIBLE, so the per-tier rollup is exercised.
    return ParameterSet(
        [
            _param("vmax", 5.0, *vmax_band),
            _param("ks", 5.0, *ks_band, tier=Tier.SPECULATIVE),
        ]
    )


def _ensemble(pset: ParameterSet, *, n=150, seed=3):
    schema = _schema()
    ps = ProcessSet(schema, [ParamFerm()])
    y0 = schema.pack({"S": 200.0, "E": 0.0, "CO2": 0.0})
    grid = np.linspace(0.0, 100.0, 25)
    return simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=n, seed=seed, t_eval=grid), pset


# time index into the 25-point grid where the response is still ~linear in vmax
_EARLY = 6


def test_budget_accounts_for_spread_and_names_the_driver():
    ens, pset = _ensemble(_pset())
    attr = attribute_spread(ens, "E", pset.tier_map(), time_index=_EARLY)

    # Near-linear early window: R^2 close to 1, so almost nothing is left unexplained.
    assert attr.r_squared > 0.9
    # The SRC^2 shares plus the unexplained remainder form a variance budget summing to
    # ~1 (exactly under a linear, independent-input model; approximate otherwise).
    assert sum(attr.per_param.values()) + attr.unexplained == pytest.approx(1.0, abs=0.1)
    # vmax is the driver, and raising it raises ethanol (positive SRC).
    assert attr.per_param["vmax"] > attr.per_param["ks"]
    assert attr.per_param_signed["vmax"] > 0.0
    assert attr.ranked()[0][0] == "vmax"
    # The per-tier rollup is just the per-parameter shares grouped by tier.
    assert attr.per_tier[Tier.PLAUSIBLE] == pytest.approx(attr.per_param["vmax"])
    assert attr.per_tier[Tier.SPECULATIVE] == pytest.approx(attr.per_param["ks"])
    assert sum(attr.per_tier.values()) == pytest.approx(sum(attr.per_param.values()))


def test_widening_a_band_raises_that_parameters_share():
    # Give ks the dominant relative spread and vmax almost none: the driver flips to ks.
    ens, pset = _ensemble(_pset(vmax_band=(4.9, 5.1), ks_band=(1.0, 9.0)))
    attr = attribute_spread(ens, "E", pset.tier_map(), time_index=_EARLY)
    assert attr.per_param["ks"] > attr.per_param["vmax"]
    assert attr.ranked()[0][0] == "ks"


def test_srrc_method_runs_and_unknown_method_raises():
    ens, pset = _ensemble(_pset())
    src = attribute_spread(ens, "E", pset.tier_map(), time_index=_EARLY, method="src")
    srrc = attribute_spread(ens, "E", pset.tier_map(), time_index=_EARLY, method="srrc")
    assert src.method == "src" and srrc.method == "srrc"
    # Both should agree that vmax is the driver on this monotone response.
    assert srrc.ranked()[0][0] == "vmax"
    with pytest.raises(ValueError, match="unknown method"):
        attribute_spread(ens, "E", pset.tier_map(), method="anova")


def test_pinned_parameter_is_excluded_from_the_budget():
    # A zero-width band for ks: it is in sampled_names (ParamFerm reads it) but drawn
    # constant, so it explains no variance and must be dropped, not divided-by-zero.
    ens, pset = _ensemble(_pset(ks_band=(5.0, 5.0)))
    assert "ks" in ens.sampled_names
    attr = attribute_spread(ens, "E", pset.tier_map(), time_index=_EARLY)
    assert "ks" not in attr.per_param
    assert "vmax" in attr.per_param


def test_no_spread_is_a_degenerate_budget():
    # Sampling nothing -> every member is the nominal run -> the output has no spread, so
    # there is nothing to attribute (empty budget, zero unexplained, not a crash).
    ens, pset = _ensemble(_pset())
    ens_flat = simulate_ensemble(
        ProcessSet(_schema(), [ParamFerm()]),
        pset,
        _schema().pack({"S": 200.0, "E": 0.0, "CO2": 0.0}),
        (0.0, 100.0),
        n_members=6,
        seed=0,
        t_eval=np.linspace(0.0, 100.0, 25),
        only=[],
    )
    attr = attribute_spread(ens_flat, "E", pset.tier_map())
    assert attr.per_param == {}
    assert attr.per_tier == {}
    assert attr.r_squared == 0.0
    assert attr.unexplained == 0.0


def test_underdetermined_fit_raises():
    # Fewer members than varying parameters -> the OLS fit is underdetermined.
    ens, pset = _ensemble(_pset(), n=2)
    with pytest.raises(ValueError, match="underdetermined"):
        attribute_spread(ens, "E", pset.tier_map(), time_index=_EARLY)


def test_bad_variable_and_slot_raise():
    ens, pset = _ensemble(_pset())
    with pytest.raises(KeyError, match="no variable"):
        attribute_spread(ens, "nope", pset.tier_map())
    with pytest.raises(ValueError, match="slot"):
        attribute_spread(ens, "E", pset.tier_map(), slot=3)
