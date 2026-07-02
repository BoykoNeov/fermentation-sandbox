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
from fermentation.runtime import (
    Ensemble,
    ScheduledEvent,
    sample_parameters,
    simulate,
    simulate_ensemble,
    simulate_scheduled,
)
from fermentation.scenario import (
    Intervention,
    Scenario,
    TemperaturePoint,
    compile_scenario,
)
from fermentation.validation import assert_conserved, total_carbon, total_mass, total_nitrogen

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


# -- LHS / Sobol low-discrepancy samplers -------------------------------------


def _mean_e_at(pset: ParameterSet, sampler: str, *, n: int, seed: int, idx: int) -> float:
    ps = _toy_ps()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(
        ps, pset, y0, (0.0, 100.0), n_members=n, seed=seed, t_eval=grid, sampler=sampler
    )
    return float(np.mean(ens.members[:, ps.schema.slice("E").start, idx]))


def test_qmc_samplers_run_conserve_and_stay_in_band():
    for sampler in ("lhs", "sobol"):
        ps, pset = _toy_ps(), _toy_pset()
        y0 = _toy_y0(ps.schema)
        grid = _short_grid()
        ens = simulate_ensemble(
            ps, pset, y0, (0.0, 100.0), n_members=32, seed=0, t_eval=grid, sampler=sampler
        )
        assert ens.sampler == sampler
        assert ens.n_succeeded == 32
        conserved = total_mass(ps.schema)
        for i in range(ens.n_succeeded):
            # The crown-jewel per-member invariant is sampler-agnostic.
            assert_conserved(ens.member_trajectory(i), conserved, rtol=1e-5, atol=1e-6)
            # Draws respect the provenance bands, same as MC.
            assert 3.0 <= ens.member_params[i]["vmax"] <= 7.0
            assert 4.0 <= ens.member_params[i]["ks"] <= 6.0


def test_qmc_center_matches_mc():
    # A low-discrepancy sequence should not shift the center — only cover the band more
    # evenly. The median trajectory must agree with MC to a loose tolerance.
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    mc = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=128, seed=1, t_eval=grid)
    for sampler in ("lhs", "sobol"):
        q = simulate_ensemble(
            ps, pset, y0, (0.0, 100.0), n_members=128, seed=1, t_eval=grid, sampler=sampler
        )
        e_mc = mc.median()[ps.schema.slice("E").start]
        e_q = q.median()[ps.schema.slice("E").start]
        assert np.allclose(e_mc, e_q, atol=1.0)


def test_qmc_reduces_estimator_variance():
    # The point of LHS/Sobol: at a fixed member budget the estimator (here the ensemble
    # mean of E) is far more stable seed-to-seed than i.i.d. MC. Observed ~8x; assert a
    # comfortable margin so this is not flaky.
    pset = _toy_pset()
    idx = 12
    seeds = range(12)
    mc_std = np.std([_mean_e_at(pset, "mc", n=24, seed=s, idx=idx) for s in seeds])
    lhs_std = np.std([_mean_e_at(pset, "lhs", n=24, seed=s, idx=idx) for s in seeds])
    assert mc_std > 0.0
    assert lhs_std < 0.5 * mc_std


def test_qmc_reproducible_by_seed():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    for sampler in ("lhs", "sobol"):
        a = simulate_ensemble(
            ps, pset, y0, (0.0, 100.0), n_members=16, seed=7, t_eval=grid, sampler=sampler
        )
        b = simulate_ensemble(
            ps, pset, y0, (0.0, 100.0), n_members=16, seed=7, t_eval=grid, sampler=sampler
        )
        c = simulate_ensemble(
            ps, pset, y0, (0.0, 100.0), n_members=16, seed=8, t_eval=grid, sampler=sampler
        )
        assert np.array_equal(a.members, b.members)
        assert not np.array_equal(a.members, c.members)


def test_qmc_pins_zero_width_bands():
    # A zero-width read parameter is in sampled_names but must NOT take a QMC dimension —
    # it stays at its nominal value in every member while the real band still varies.
    ps = _toy_ps()
    pset = _toy_pset(ks_band=(5.0, 5.0))
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(
        ps, pset, y0, (0.0, 100.0), n_members=16, seed=0, t_eval=grid, sampler="lhs"
    )
    assert "ks" in ens.sampled_names  # read, so in scope
    assert all(mp["ks"] == 5.0 for mp in ens.member_params)  # but pinned, never drawn
    assert np.std([mp["vmax"] for mp in ens.member_params]) > 0.0  # the real band varied


def test_sobol_requires_power_of_two():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    with pytest.raises(ValueError, match="power of two"):
        simulate_ensemble(
            ps, pset, y0, (0.0, 100.0), n_members=100, seed=0, t_eval=grid, sampler="sobol"
        )


def test_unknown_sampler_raises():
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    with pytest.raises(ValueError, match="unknown sampler"):
        simulate_ensemble(
            ps, pset, y0, (0.0, 100.0), n_members=8, seed=0, t_eval=grid, sampler="halton"
        )


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


@pytest.mark.parametrize("sampler", ["mc", "lhs", "sobol"])
def test_wine_ensemble_scopes_to_active_reads_and_conserves_carbon(sampler):
    # The crown-jewel per-member conservation is enforced on EVERY sampler path (not just
    # the MC one it was born on): LHS/Sobol produce member_params via a different mechanism
    # (inverse-CDF, varying/pinned split), so carbon and nitrogen must be shown to close
    # there too. n_members=8 is a power of two, valid for Sobol.
    compiled = compile_scenario(_wine_scenario())
    grid = np.linspace(0.0, compiled.t_span_h[1], 20)
    ens = simulate_ensemble(
        compiled.process_set,
        compiled.parameters,
        compiled.y0,
        compiled.t_span_h,
        sampler=sampler,
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

    # The nitrogen ledger closes per-member too. As cells grow they assimilate free
    # ``N`` into biomass, so ``N + biomass_N_fraction * X`` is the invariant. Like carbon,
    # the check must read the member's OWN sampled ``biomass_N_fraction`` (the growth
    # Process draws N against it, and it is itself sampled) or genuine closure reads as
    # drift. There is no aa-ledger yet (deferred, D-23) and fusels route carbon — not N —
    # from sugar, so nothing sinks N except biomass and the balance closes to ~1e-12.
    nominal_fn = compiled.parameters.value("biomass_N_fraction")
    for i in range(ens.n_succeeded):
        fn = ens.member_params[i].get("biomass_N_fraction", nominal_fn)
        nitrogen = total_nitrogen(compiled.schema, biomass_nitrogen_fraction=fn)
        assert_conserved(
            ens.member_trajectory(i), nitrogen, rtol=1e-6, atol=1e-6, label=f"nitrogen[{i}]"
        )

    # A real, ordered ethanol band that brackets a plausible dry-wine ABV region.
    e_band = ens.band("E")
    assert np.all(e_band.low <= e_band.high + 1e-9)
    assert e_band.median[-1] > 100.0  # ~130-150 g/L ethanol by day 14


# -- scheduled ensembles: an ensemble over a multi-segment schedule (D-37) -----
#
# ``events`` threads the D-35/D-36 intervention machinery through the stochastic wrapper:
# every sampled member is integrated through the SAME schedule, so the spread is the
# parameter-uncertainty band of a scheduled scenario. These tests pin the three things
# scheduling adds — schedule-union sampling scope, per-member Process-set isolation, and the
# per-member across-jumps conservation ledger — plus the un-scheduled isolability.


def test_unscheduled_ensemble_has_empty_ledger_and_default_bounds():
    # Isolability: with no events the ensemble is byte-for-byte the pre-scheduling one, and the
    # new fields degenerate cleanly — empty flows, a single [t0, t_end] segment.
    ps, pset = _toy_ps(), _toy_pset()
    y0 = _toy_y0(ps.schema)
    grid = _short_grid()
    ens = simulate_ensemble(ps, pset, y0, (0.0, 100.0), n_members=5, seed=0, t_eval=grid)
    assert ens.segment_bounds == (0.0, 100.0)
    assert ens.member_flows == ((),) * ens.n_succeeded
    assert ens.nominal_flows == ()
    # member_trajectory / nominal_trajectory still audit (empty ledger ⇒ constant total).
    conserved = total_mass(ps.schema)
    for i in range(ens.n_succeeded):
        assert_conserved(ens.member_trajectory(i), conserved, rtol=1e-5, atol=1e-6)
    assert np.array_equal(ens.nominal_trajectory().y, ens.nominal)


def _dap_rack_scenario() -> Scenario:
    # An injection (DAP, +N at day 2) AND a removal (rack, −C/−N at day 15) on one run — the
    # only scenario that populates both sides of the external-flow ledger. Autolysis is opted
    # in so the rack has a non-empty debris pool to draw off.
    return Scenario(
        name="dap-rack-ensemble",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 100.0, "pitch_gpl": 0.25, "autolysis_rate_per_h": 0.002},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[
            Intervention(day=2.0, action="add_dap", params={"dap_gpl": 0.4}),
            Intervention(day=15.0, action="rack", params={"fraction": 0.8}),
        ],
        duration_days=20.0,
    )


def test_scheduled_ensemble_conserves_across_jumps_per_member():
    # The crown-jewel ledger (D-36) extended to the ensemble: EVERY sampled member satisfies
    # the run-wide identity final == initial + Σ external_flows for carbon AND nitrogen, read
    # with that member's OWN sampled accounting fractions. The flows are member-dependent (a
    # rack removes a fraction of the sampled lees mass), so they must be stored per member.
    cs = compile_scenario(_dap_rack_scenario())
    grid = np.linspace(0.0, cs.t_span_h[1], 30)
    ens = cs.run_ensemble(n_members=8, seed=0, t_eval=grid)

    # Both breakpoints (DAP@2d, rack@15d) are in the shared segment bounds; two flows per member.
    assert len(ens.segment_bounds) == 4  # t0, 48h, 360h, t_end
    assert all(len(flows) == 2 for flows in ens.member_flows)  # one dose + one rack each

    nominal_fc = cs.parameters.value("biomass_C_fraction")
    nominal_fn = cs.parameters.value("biomass_N_fraction")
    rack_carbon = []
    for i in range(ens.n_succeeded):
        fc = ens.member_params[i].get("biomass_C_fraction", nominal_fc)
        fn = ens.member_params[i].get("biomass_N_fraction", nominal_fn)
        c_of = total_carbon(cs.schema, biomass_carbon_fraction=fc)
        n_of = total_nitrogen(cs.schema, biomass_nitrogen_fraction=fn)
        traj = ens.member_trajectory(i)  # a ScheduledTrajectory carrying THIS member's flows
        for quantity, tol in ((c_of, 1e-6), (n_of, 1e-9)):
            injected = sum(quantity(f.delta) for f in traj.external_flows)
            initial = quantity(cs.y0)
            final = quantity(traj.y[:, -1])
            assert final == pytest.approx(initial + injected, abs=tol)
        rack_carbon.append(c_of(traj.external_flows[1].delta))  # flows[1] is the rack

    # The rack removal is genuinely member-dependent — the lees mass at rack time varies with
    # the sampled death/growth kinetics, so storing one nominal ledger would have been wrong.
    assert all(rc < 0.0 for rc in rack_carbon)  # racking removes carbon
    assert np.std(rack_carbon) > 0.0  # …and by a member-varying amount


def _pitch_wine_scenario() -> Scenario:
    return Scenario(
        name="pitch-ensemble",
        medium="wine",
        initial={
            "brix": 22.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.25,
            "malic_gpl": 3.0,
            "initial_ph": 3.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=22.0)],
        interventions=[Intervention(day=1.0, action="pitch_mlf", params={"pitch_gpl": 0.1})],
        duration_days=25.0,
    )


_MLF_NAMES = (
    "malolactic_conversion",
    "malolactic_citrate_metabolism",
    "oenococcus_diacetyl_reduction",
)


def test_scheduled_ensemble_samples_mid_run_enabled_reads():
    # A pitch_mlf schedule enables the malolactic Processes only from the breakpoint. Their
    # kinetics drive the back half of the run, but they are DISABLED at t0 — so sampling only
    # the t0-active reads would silently miss exactly the parameters the pitched half depends
    # on. The ensemble unions reads across the whole schedule (D-37): k_mlf is in scope even
    # though no t0-active Process reads it.
    cs = compile_scenario(_pitch_wine_scenario())
    t0_reads: set[str] = set()
    for p in cs.process_set.active:
        t0_reads.update(p.reads)
    assert "k_mlf" not in t0_reads  # disabled at t0 ⇒ absent from the naive scope
    ens = cs.run_ensemble(n_members=6, seed=0, t_eval=np.linspace(0.0, cs.t_span_h[1], 20))
    assert "k_mlf" in ens.sampled_names  # …but the schedule-union scope picks it up


def test_scheduled_ensemble_restores_process_set_and_travels_tier():
    cs = compile_scenario(_pitch_wine_scenario())
    assert all(not cs.process_set.is_enabled(n) for n in _MLF_NAMES)  # unpitched at compile
    ens = cs.run_ensemble(n_members=6, seed=1, t_eval=np.linspace(0.0, cs.t_span_h[1], 20))
    # The ensemble reset the shared set before every member and leaves it PRISTINE — unlike a
    # single cs.run(), whose enable persists. This isolation is what stops one member's pitch
    # from leaking into the next member's pre-pitch segments.
    assert all(not cs.process_set.is_enabled(n) for n in _MLF_NAMES)
    # Tier travels across the mid-run reconfiguration (D-35): the malate/lactate slots the
    # pitched speculative Processes touch report speculative for the WHOLE ensemble.
    for name in ("malic", "lactic"):
        assert ens.tier_map[name] is Tier.SPECULATIVE


class _UngatedExtraFlux(Process):
    """A speculative Process, disabled at t0, that adds a constant ungated flux once enabled.

    Unlike ``pitch_mlf`` — whose enabled Processes are gated by the ``X_mlf`` catalyst and so
    contribute nothing until the pitch *mutation* — this one acts the instant it is enabled,
    with no state gate. That makes a leaked enable numerically *visible*, which is exactly what
    the isolation invariant needs to be tested against.
    """

    name = "ungated_extra_flux"
    tier = Tier.SPECULATIVE
    touches = ("E",)

    def derivatives(self, t, y, schema, params):
        d = schema.zeros()
        d[schema.slice("E")] = 1.0  # constant, ungated
        return d


def _iso_ps() -> ProcessSet:
    """A fresh toy set with the ungated extra flux present but disabled (the t0 configuration)."""
    schema = _toy_schema()
    ps = ProcessSet(schema, [ParamFermentation(), _UngatedExtraFlux()])
    ps.disable("ungated_extra_flux")
    return ps


def test_scheduled_ensemble_isolates_members_from_each_others_reconfigure():
    # The load-bearing invariant: the reconfigure is reset before EVERY member, so no member's
    # enable leaks into the next member's pre-event segments. Tested with an UNGATED reconfigure
    # (see _UngatedExtraFlux) so a leak is numerically visible — pitch_mlf's leak is inert
    # (X_mlf-gated), so per-member conservation cannot catch this and a byte-for-byte comparison
    # to independent fresh-set runs is the only discriminating check. Fails without the reset.
    event = ScheduledEvent(
        time_h=50.0,
        label="enable_extra",
        reconfigure=lambda ps: ps.enable("ungated_extra_flux"),
    )
    pset = _toy_pset()
    y0 = _toy_y0(_toy_schema())
    grid = _short_grid()
    ens = simulate_ensemble(
        _iso_ps(), pset, y0, (0.0, 100.0), n_members=4, seed=0, t_eval=grid, events=[event]
    )
    assert ens.n_succeeded == 4
    for i in range(ens.n_succeeded):
        ref = simulate_scheduled(
            _iso_ps(),  # a pristine set: extra disabled until the day-50 reconfigure
            ens.member_params[i],
            y0,
            (0.0, 100.0),
            events=[event],
            param_tiers=pset.tier_map(),
            t_eval=grid,
        )
        assert np.array_equal(ens.members[i], ref.y)  # a leaked prior-member enable breaks this
