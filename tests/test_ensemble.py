"""Stochastic ensemble wrapper (handoff §1.6, decision D-24).

The wrapper's contract: sample each parameter within its provenance uncertainty band,
run a Monte-Carlo ensemble of the *pure* deterministic core, and report the nominal
run plus a median + spread — without letting randomness leak into the core. These
tests pin that contract: a single unsampled run stays byte-for-byte reproducible, the
ensemble is seed-reproducible, the spread tracks the input band width, every sampled
member still conserves mass, and failed members are counted (never silently dropped
into a survivorship-biased band).
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.chemistry import CO2_PER_HEXOSE, ETHANOL_PER_HEXOSE
from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import FloatArray, StateSchema, VarSpec
from fermentation.core.tiers import Tier
from fermentation.parameters.schema import Parameter, Provenance, Uncertainty
from fermentation.parameters.store import ParameterSet
from fermentation.runtime import Ensemble, sample_parameters, simulate, simulate_ensemble
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon, total_mass

# -- a toy that actually READS a sampled parameter ----------------------------
#
# The conftest toy takes vmax as an __init__ arg, so sampling would not move it. This
# one reads ``vmax``/``ks`` from the resolved param map, so a sampled draw changes the
# trajectory — exactly what the ensemble must propagate. Mass S+E+CO2 is conserved by
# the Gay-Lussac split regardless of the sampled rate, which is the per-member invariant.


class ParamFermentation(Process):
    name = "param_fermentation"
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


class PickyFermentation(ParamFermentation):
    """Same conserving kinetics, but raises above a vmax threshold — to exercise the
    failed-member accounting (a sampled RHS that *raises*, like the uptake carbon guard)."""

    name = "picky_fermentation"

    def __init__(self, max_vmax: float) -> None:
        self.max_vmax = max_vmax

    def derivatives(self, t, y, schema, params):
        if params["vmax"] > self.max_vmax:
            raise ValueError(f"vmax {params['vmax']} exceeds {self.max_vmax}")
        return super().derivatives(t, y, schema, params)


def _toy_schema() -> StateSchema:
    return StateSchema([VarSpec("S", "g/L"), VarSpec("E", "g/L"), VarSpec("CO2", "g/L")])


def _param(
    name: str, value: float, low: float, high: float, tier: Tier = Tier.PLAUSIBLE
) -> Parameter:
    return Parameter(
        name=name,
        value=value,
        unit="1/h",
        tier=tier,
        uncertainty=Uncertainty(low=low, high=high),
        provenance=Provenance(source="author estimate", conditions="toy test"),
    )


def _toy_pset(*, vmax_band=(3.0, 7.0), ks_band=(4.0, 6.0)) -> ParameterSet:
    return ParameterSet(
        [
            _param("vmax", 5.0, *vmax_band),
            _param("ks", 5.0, *ks_band),
            # a zero-width band: must be pinned, never drawn
            _param("inert", 1.0, 1.0, 1.0),
        ]
    )


def _toy_ps(process: Process | None = None) -> ProcessSet:
    schema = _toy_schema()
    return ProcessSet(schema, [process or ParamFermentation()])


def _toy_y0(schema: StateSchema) -> FloatArray:
    return schema.pack({"S": 200.0, "E": 0.0, "CO2": 0.0})


# -- sample_parameters --------------------------------------------------------


def test_sample_pins_zero_width_and_leaves_unnamed_at_nominal():
    pset = _toy_pset()
    rng = np.random.default_rng(0)
    # names=None samples everything; the zero-width 'inert' band pins to its value.
    full = sample_parameters(pset, rng, names=None)
    assert full["inert"] == 1.0
    # Restricting names leaves the others at their nominal resolved value.
    rng = np.random.default_rng(0)
    just_vmax = sample_parameters(pset, rng, names=["vmax"])
    assert just_vmax["ks"] == 5.0  # untouched nominal
    assert just_vmax["inert"] == 1.0
    assert 3.0 <= just_vmax["vmax"] <= 7.0  # the one named param was drawn from its band


def test_sample_stays_within_bands_and_uniform_available():
    pset = _toy_pset(vmax_band=(3.0, 7.0))
    for dist in ("triangular", "uniform"):
        rng = np.random.default_rng(1)
        draws = [
            sample_parameters(pset, rng, distribution=dist, names=["vmax"])["vmax"]
            for _ in range(500)
        ]
        assert min(draws) >= 3.0 and max(draws) <= 7.0
        # a real spread was produced (not all pinned)
        assert np.std(draws) > 0.1
    with pytest.raises(ValueError, match="unknown distribution"):
        sample_parameters(pset, np.random.default_rng(0), distribution="normal")


def test_sample_draw_order_is_deterministic():
    # Same seed -> identical sample, independent of set/hash ordering (names sorted).
    pset = _toy_pset()
    a = sample_parameters(pset, np.random.default_rng(42), names=["vmax", "ks"])
    b = sample_parameters(pset, np.random.default_rng(42), names=["vmax", "ks"])
    assert a == b


# -- simulate_ensemble: determinism & reproducibility -------------------------


def _short_grid() -> FloatArray:
    return np.linspace(0.0, 100.0, 25)


def test_nominal_equals_deterministic_simulate():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=5, seed=0, t_eval=grid)
    # The nominal member must reproduce a plain simulate() on the resolved values,
    # byte-for-byte, with identical solver kwargs + grid (decision D-24).
    ref = simulate(ps, pset.resolve(), y0, (0.0, 100.0), param_tiers=pset.tier_map(), t_eval=grid)
    assert np.array_equal(ens.nominal, ref.y)
    assert ens.tier_map == ref.tier_map


def test_reproducible_by_seed():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    a = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=8, seed=7, t_eval=grid)
    b = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=8, seed=7, t_eval=grid)
    c = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=8, seed=8, t_eval=grid)
    assert np.array_equal(a.members, b.members)  # same seed -> identical ensemble
    assert not np.array_equal(a.members, c.members)  # different seed -> different draws


def test_only_empty_is_a_degenerate_ensemble():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=6, seed=0, t_eval=grid, only=[])
    assert ens.sampled_names == ()
    # Nothing sampled -> every member is the nominal run; median coincides with nominal.
    for i in range(ens.n_succeeded):
        assert np.array_equal(ens.members[i], ens.nominal)
    assert np.array_equal(ens.median(), ens.nominal)


# -- spread tracks the input band width ---------------------------------------


def _final_e_band_width(pset: ParameterSet) -> float:
    ps = _toy_ps()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=120, seed=3, t_eval=grid)
    band = ens.band("E")
    return float(band.high[-1] - band.low[-1])


def test_spread_tracks_band_width():
    narrow = _final_e_band_width(_toy_pset(vmax_band=(4.5, 5.5)))
    wide = _final_e_band_width(_toy_pset(vmax_band=(2.0, 8.0)))
    assert narrow > 0.0  # a real, nonzero spread
    assert wide > narrow  # widening the input band widens the output band


def test_band_brackets_are_ordered_and_contain_median():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=100, seed=2, t_eval=grid)
    band = ens.band("E", low=5.0, high=95.0)
    assert np.all(band.low <= band.median + 1e-9)
    assert np.all(band.median <= band.high + 1e-9)
    with pytest.raises(ValueError, match="0 <= low <= high <= 100"):
        ens.band("E", low=90.0, high=10.0)


# -- per-member conservation (the crown-jewel invariant) ----------------------


def test_every_member_conserves_mass():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=30, seed=5, t_eval=grid)
    conserved = total_mass(ps.schema)
    for i in range(ens.n_succeeded):
        # member_trajectory lets the deterministic harness audit each sampled draw.
        assert_conserved(
            ens.member_trajectory(i), conserved, rtol=1e-5, atol=1e-6, label=f"mass[{i}]"
        )


# -- failed-member accounting -------------------------------------------------


def test_partial_failures_are_counted_not_dropped():
    # PickyFermentation raises when a sampled vmax exceeds 6.0; the band [3,7] means some
    # draws fail. They must be counted, and n_succeeded + n_failed == n_requested.
    ps = _toy_ps(PickyFermentation(max_vmax=6.0))
    pset = _toy_pset(vmax_band=(3.0, 7.0))
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(
        ps, pset, y0, (0.0, 100.0), n_members=60, seed=4, t_eval=grid, max_failure_fraction=1.0
    )
    assert ens.n_failed > 0
    assert ens.n_succeeded > 0
    assert ens.n_succeeded + ens.n_failed == ens.n_requested == 60
    assert len(ens.failures) == ens.n_failed
    assert all("raised" in f for f in ens.failures)


def test_survivorship_threshold_raises():
    # Threshold at the nominal vmax (5.0) with a right-skewed [3,9] band (mode 5): the
    # nominal run passes (5.0 is not > 5.0) but ~2/3 of draws exceed it and fail, so the
    # default 50% survivorship guard must refuse to return a biased spread.
    ps = _toy_ps(PickyFermentation(max_vmax=5.0))
    pset = _toy_pset(vmax_band=(3.0, 9.0))
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    with pytest.raises(RuntimeError, match="survivorship"):
        simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=40, seed=1, t_eval=grid)


# -- end-to-end on a real provenance-backed wine scenario ---------------------


def _wine_scenario(**overrides) -> Scenario:
    kwargs: dict[str, object] = {
        "name": "wine-ensemble",
        "medium": "wine",
        "initial": {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5},
        "temperature_schedule": [TemperaturePoint(day=0.0, celsius=20.0)],
        "duration_days": 14.0,
    }
    kwargs.update(overrides)
    return Scenario(**kwargs)


def test_wine_ensemble_scopes_to_active_reads_and_conserves_carbon():
    compiled = compile_scenario(_wine_scenario())
    grid = np.linspace(0.0, compiled.t_span_h[1], 20)
    ens = simulate_ensemble(
        compiled.process_set,
        compiled.parameters,
        compiled.y0,
        compiled.t_span_h,
        n_members=8,
        seed=0,
        t_eval=grid,
    )
    assert isinstance(ens, Ensemble)

    # Sampling is scoped to what the ACTIVE Process set reads (decision D-24): a subset
    # of the provenance params, and never a parameter no active Process consumes.
    active_reads: set[str] = set()
    for p in compiled.process_set.active:
        active_reads.update(p.reads)
    for m in compiled.process_set.active_modifiers:
        active_reads.update(m.reads)
    assert set(ens.sampled_names) <= active_reads
    assert "q_sugar_max" in ens.sampled_names  # a kinetic driver is in scope
    # MLF is undosed here, so the pKa set (its reads) is NOT active and NOT sampled —
    # the D-24 no-op-avoidance and the pH-anchor safety both fall out of the scoping.

    # Every sampled member conserves carbon to solver tolerance — a sampled parameter
    # set must never break a balance (prime directive; the structural draw-from-S
    # accounting closes regardless of the rate constants drawn). The check must use the
    # member's OWN biomass_C_fraction (an accounting constant the growth Process draws
    # sugar carbon against, and which is itself sampled) or genuine closure reads as drift.
    nominal_fc = compiled.parameters.value("biomass_C_fraction")
    for i in range(ens.n_succeeded):
        fc = ens.member_params[i].get("biomass_C_fraction", nominal_fc)
        carbon = total_carbon(compiled.schema, biomass_carbon_fraction=fc)
        assert_conserved(
            ens.member_trajectory(i), carbon, rtol=1e-6, atol=1e-6, label=f"carbon[{i}]"
        )

    # A real, ordered ethanol band that brackets a plausible dry-wine ABV region.
    e_band = ens.band("E")
    assert np.all(e_band.low <= e_band.high + 1e-9)
    assert e_band.median[-1] > 100.0  # ~130-150 g/L ethanol by day 14
