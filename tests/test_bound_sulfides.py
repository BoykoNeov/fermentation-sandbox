"""Tests for the bottle-reduction sulfide-release Processes (decision D-135).

:class:`BoundHydrogenSulfideRelease` and :class:`BoundMethanethiolRelease` empty the
metal-complexed ("bonded") H2S / methanethiol reservoirs into the free pools during anaerobic
bottle aging — the beat D-101 parked as unbuildable and
:mod:`fermentation.core.kinetics.mercaptans` still names as its own missing route.

**What these tests are really guarding.** The beat's whole claim to fidelity is that its rate
constants and reservoir levels are not invented: they come from Franco-Luesma & Ferreira 2016
(*Food Chemistry* 199:42-50), which reports the same physics through three mutually independent
measurements. So the tests below are mostly *cross-checks between parts of that paper that the
paper itself never compared*, re-derived here from the shipped YAML:

1. ``k x reservoir`` must reproduce the paper's separately-regressed **free-form** accumulation
   (+0.38 ug/L/yr H2S, +0.23 ug/L/yr MeSH). For H2S it lands at 98 %.
2. Where it does NOT reproduce it — MeSH, at 49 % — the shortfall must equal the paper's
   independently-measured **release share** (47.5 %), because the missing half is the de novo
   route this beat deliberately does not model. **The under-claim is asserted, not apologised
   for**: a test that let the MeSH coverage drift up to 100 % would be hiding the gap.

Every shipped constant is READ from ``bound_sulfides.yaml`` rather than restated (the D-100
lesson: a test that hard-codes the value it should be reading is a test of itself). The numbers
that ARE literal here are the *paper's observations* — the things the model is checked against.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.kinetics import BoundHydrogenSulfideRelease, BoundMethanethiolRelease
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.units.convert import gpl_to_ugl, ugl_to_gpl

#: Hours per Julian year — the conversion used to turn the paper's %/yr regression slopes into the
#: shipped 1/h constants, re-derived here rather than imported so the test cannot inherit an error
#: from the code it checks.
_HOURS_PER_YEAR = 365.25 * 24.0


@pytest.fixture
def params() -> dict[str, float]:
    store = load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "bound_sulfides.yaml"
    )
    return store.resolve()


def _state(
    schema: StateSchema, *, bound_h2s_ugl: float = 0.0, bound_mesh_ugl: float = 0.0
) -> FloatArray:
    y = schema.zeros()
    y[schema.slice("bound_h2s")] = ugl_to_gpl(bound_h2s_ugl)
    y[schema.slice("bound_methanethiol")] = ugl_to_gpl(bound_mesh_ugl)
    return y


def _released_ugl(params: Mapping[str, float], k_name: str, reservoir_ugl: float, hours: float):
    """Closed-form first-order release — the free-form gain the paper regresses."""
    return reservoir_ugl * (1.0 - math.exp(-params[k_name] * hours))


# --------------------------------------------------------------------------------------
# The Process contract
# --------------------------------------------------------------------------------------


def test_touches_only_each_reservoir_and_its_own_free_pool():
    # The isolability claim in two assertions. No o2 (complex dissociation is not an oxidation —
    # this is the ANOXIC route, measured under strict argon), no sugar/ethanol, and crucially each
    # Process touches only ITS OWN species: a lumped sulfide Process could not carry the measured
    # 4.3x rate asymmetry between the two.
    assert BoundHydrogenSulfideRelease.touches == ("bound_h2s", "h2s")
    assert BoundMethanethiolRelease.touches == ("bound_methanethiol", "methanethiol")
    assert BoundHydrogenSulfideRelease.tier is Tier.SPECULATIVE
    assert BoundMethanethiolRelease.tier is Tier.SPECULATIVE


def test_strict_process_set_accepts_the_touches_contract(params):
    schema = wine_schema()
    pset = ProcessSet(
        schema, [BoundHydrogenSulfideRelease(), BoundMethanethiolRelease()], strict=True
    )
    y = _state(schema, bound_h2s_ugl=19.7, bound_mesh_ugl=1.4)
    d = pset.total_derivatives(0.0, y, params)
    assert d[schema.slice("h2s")][0] > 0.0
    assert d[schema.slice("methanethiol")][0] > 0.0


def test_tier_of_released_pools_is_speculative(params):
    schema = wine_schema()
    pset = ProcessSet(schema, [BoundHydrogenSulfideRelease(), BoundMethanethiolRelease()])
    assert pset.tier_of("h2s") is Tier.SPECULATIVE
    assert pset.tier_of("methanethiol") is Tier.SPECULATIVE


def test_unseeded_reservoirs_make_both_processes_byte_for_byte_inert(params):
    # Prime directive 3: an explicit 0 reservoir (or an older ParameterSet that seeds neither slot)
    # must leave the derivative of EVERY column exactly zero — not small, zero.
    schema = wine_schema()
    pset = ProcessSet(schema, [BoundHydrogenSulfideRelease(), BoundMethanethiolRelease()])
    d = pset.total_derivatives(0.0, _state(schema), params)
    assert np.array_equal(d, np.zeros_like(d))


def test_release_is_temperature_flat_by_decision(params):
    # NOT an omission — the paper never states the storage temperature behind its ambient
    # regression, so no E_a is fitted and none is shipped (see bound_sulfides.yaml's header). This
    # test pins that as a DELIBERATE property: the rate must be identical at cellar and at cellar
    # +20 K, and neither Process may read T_ref or an activation energy.
    schema = wine_schema()
    cold = _state(schema, bound_h2s_ugl=19.7, bound_mesh_ugl=1.4)
    hot = cold.copy()
    cold[schema.slice("T")] = 285.15
    hot[schema.slice("T")] = 305.15
    pset = ProcessSet(schema, [BoundHydrogenSulfideRelease(), BoundMethanethiolRelease()])
    assert np.array_equal(
        pset.total_derivatives(0.0, cold, params), pset.total_derivatives(0.0, hot, params)
    )
    for proc in (BoundHydrogenSulfideRelease, BoundMethanethiolRelease):
        assert not any("E_a" in name or name == "T_ref" for name in proc.reads)


# --------------------------------------------------------------------------------------
# 1:1 transfer, and the ledger split between the two species
# --------------------------------------------------------------------------------------


def test_each_release_conserves_its_species_exactly(params):
    # The ligand is the same molecule either way — only its binding state changes — so what leaves
    # the reservoir must enter the free pool to machine precision, with no yield parameter and no
    # molar-mass conversion (the dms_potential idiom).
    schema = wine_schema()
    y = _state(schema, bound_h2s_ugl=19.7, bound_mesh_ugl=1.4)
    d = ProcessSet(
        schema, [BoundHydrogenSulfideRelease(), BoundMethanethiolRelease()]
    ).total_derivatives(0.0, y, params)
    assert d[schema.slice("bound_h2s")][0] == -d[schema.slice("h2s")][0]
    assert d[schema.slice("bound_methanethiol")][0] == -d[schema.slice("methanethiol")][0]


def test_thiol_release_is_carbon_neutral_but_h2s_release_is_off_ledger(params):
    # The deliberate asymmetry (D-135): methanethiol carries one carbon and D-45 weights the free
    # pool, so the BONDED thiol is weighted identically and the transfer closes total_carbon
    # exactly. Weighting it at 0 instead would read as carbon created from nothing on every step —
    # which is precisely what this test would catch. H2S is carbon-free, so BOTH its slots are off
    # the ledger and its release moves no carbon at all.
    from fermentation.validation.conservation import total_carbon

    schema = wine_schema()
    carbon = total_carbon(schema, biomass_carbon_fraction=0.45)

    # Probe the ledger weights directly by pushing 1 g/L through one slot at a time.
    def weight_of(slot: str) -> float:
        probe = schema.zeros()
        probe[schema.slice(slot)] = 1.0
        return carbon(probe)

    assert weight_of("bound_methanethiol") == weight_of("methanethiol") > 0.0
    assert weight_of("bound_h2s") == weight_of("h2s") == 0.0

    # The release itself therefore moves zero net carbon: the thiol's carbon leaves one weighted
    # slot and enters another at the identical weight, and H2S's slots carry none.
    y = _state(schema, bound_h2s_ugl=19.7, bound_mesh_ugl=1.4)
    d = ProcessSet(
        schema, [BoundHydrogenSulfideRelease(), BoundMethanethiolRelease()]
    ).total_derivatives(0.0, y, params)
    assert carbon(d) == pytest.approx(0.0, abs=1e-18)


# --------------------------------------------------------------------------------------
# THE CROSS-CHECKS — the paper's own numbers, compared against each other
# --------------------------------------------------------------------------------------


def test_shipped_rates_reproduce_the_papers_percent_per_year_release(params):
    # The regression the constants were derived FROM: 1.9 %/yr of bonded H2S and 8.1 %/yr of bonded
    # MeSH released in real bottles (Franco-Luesma & Ferreira 2016, 16 reds across vintages).
    h2s_pct = 1.0 - math.exp(-params["k_bound_h2s_release"] * _HOURS_PER_YEAR)
    mesh_pct = 1.0 - math.exp(-params["k_bound_methanethiol_release"] * _HOURS_PER_YEAR)
    assert h2s_pct == pytest.approx(0.019, abs=5e-4)
    assert mesh_pct == pytest.approx(0.081, abs=5e-4)


def test_species_asymmetry_survives_in_the_shipped_constants(params):
    # Bonded MeSH releases ~4.3x faster than bonded H2S in the ambient regression. The ordering is
    # the load-bearing part (it also holds, at 1.9x, in the independent 50 C experiment) and is the
    # reason the beat ships two pools and two constants rather than one lumped sulfide reservoir.
    ratio = params["k_bound_methanethiol_release"] / params["k_bound_h2s_release"]
    assert ratio == pytest.approx(4.3, abs=0.2)


def test_h2s_release_reproduces_the_independently_regressed_free_form_slope(params):
    # THE CROSS-CHECK THE PAPER NEVER RAN ON ITSELF. Its bonded-fraction regression (1.9 %/yr) and
    # its free-form regression (+0.38 +/- 0.11 ug/L/yr) are independent observables. k x reservoir
    # must land inside the measured band — and it does, at 98 % of the central value. Because the
    # reservoir enters multiplicatively, this simultaneously corroborates bound_h2s_initial.
    predicted = _released_ugl(
        params, "k_bound_h2s_release", params["bound_h2s_initial"], _HOURS_PER_YEAR
    )
    assert predicted == pytest.approx(0.38, abs=0.11)


def test_thiol_release_under_predicts_by_exactly_the_missing_de_novo_route(params):
    # THE UNDER-CLAIM, ASSERTED RATHER THAN CONFESSED. Release-only recovers just ~49 % of the
    # paper's measured +0.23 ug/L/yr free-MeSH accumulation. That is not a calibration failure to
    # be tuned away: the missing 51 % is de novo formation, whose mechanism the authors themselves
    # list as an open question, and 49 % agrees with the INDEPENDENTLY measured 47.5 % release
    # share from the 50 C mass balance (Table 4) to 1.5 points.
    #
    # If this test ever starts failing high, someone has quietly given the model a de novo route.
    measured_free_gain_ugl_per_yr = 0.23
    predicted = _released_ugl(
        params,
        "k_bound_methanethiol_release",
        params["bound_methanethiol_initial"],
        _HOURS_PER_YEAR,
    )
    coverage = predicted / measured_free_gain_ugl_per_yr
    assert coverage == pytest.approx(params["release_share_methanethiol_red"], abs=0.05)
    assert coverage < 0.6, "release-only must remain a LOWER BOUND on aged-wine methanethiol"


def test_recorded_release_shares_rank_the_coverage_gap_correctly(params):
    # The four mass-balance shares are carried as data so the gap is checkable, not prose. Their
    # ORDERING is the fact that matters: this model is nearly complete for red-wine H2S and barely
    # a quarter of the story for white-wine MeSH.
    red_h2s = params["release_share_h2s_red"]
    white_h2s = params["release_share_h2s_white"]
    red_mesh = params["release_share_methanethiol_red"]
    white_mesh = params["release_share_methanethiol_white"]
    assert red_h2s > white_h2s > red_mesh > white_mesh
    assert red_h2s == pytest.approx(0.903)
    assert white_mesh == pytest.approx(0.241)
    # The worst case is bad enough to be worth a named assertion: a white wine's modelled thiol is
    # a ~4x under-prediction, and nothing downstream should present it as a total.
    assert 1.0 / white_mesh == pytest.approx(4.15, abs=0.1)


# --------------------------------------------------------------------------------------
# Reservoirs, seeding, and the aging gate
# --------------------------------------------------------------------------------------


def test_seeded_reservoirs_match_the_measured_bonded_fractions(params):
    # Table 1's red-wine averages: total 20.8 / free 1.1 H2S (94 % bonded) and total 2.3 / free 0.9
    # MeSH (62 % bonded). The reservoir IS total minus free, so re-deriving the bonded fraction
    # from the shipped seed and the paper's totals must return the paper's percentages.
    assert params["bound_h2s_initial"] / 20.8 == pytest.approx(0.94, abs=0.01)
    # 1.4/2.3 = 60.9 %, against Table 1's printed 62 %: the paper averages the per-wine bonded
    # fractions, which is not the same as the ratio of the averaged totals. A 1-point gap between
    # those two orderings of the same arithmetic is expected, and the tolerance says so.
    assert params["bound_methanethiol_initial"] / 2.3 == pytest.approx(0.62, abs=0.02)
    # And the H2S reservoir dwarfs the thiol one ~14x, which is why bottle reduction reads as
    # H2S-dominated by mass despite the thiol moving faster.
    assert params["bound_h2s_initial"] / params["bound_methanethiol_initial"] == pytest.approx(
        14.0, abs=0.5
    )


def test_a_decade_of_cellaring_moves_free_sulfides_into_the_sensible_range(params):
    # An integration-free sanity check on the beat's headline behaviour: over 10 years a red wine
    # should pick up single-digit ug/L of free H2S and well under 1 ug/L of free MeSH — a wine
    # drifting reductive, not one turning into a stink bomb. (H2S sensory threshold ~1-2 ug/L.)
    hours = 10.0 * _HOURS_PER_YEAR
    h2s = _released_ugl(params, "k_bound_h2s_release", params["bound_h2s_initial"], hours)
    mesh = _released_ugl(
        params, "k_bound_methanethiol_release", params["bound_methanethiol_initial"], hours
    )
    assert 3.0 < h2s < 4.5
    assert 0.6 < mesh < 1.2
    # The reservoirs are barely dented: bottle reduction is slow, and the model must say so.
    assert h2s / params["bound_h2s_initial"] < 0.25


def test_wine_compile_seeds_the_reservoirs_and_gates_the_processes():
    # End-to-end: the sourced levels reach y0 through the compile seam (absent scenario keys mean
    # the SOURCED level, not 0 — the D-45 hard-zero defect), and both Processes are present but
    # DISABLED until begin_aging.
    from fermentation.scenario.compile import compile_scenario
    from fermentation.scenario.schema import Scenario, TemperaturePoint

    cs = compile_scenario(
        Scenario(
            name="d135-seed",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 200.0, "pitch_gpl": 0.25},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    schema = cs.schema
    assert gpl_to_ugl(float(cs.y0[schema.slice("bound_h2s")][0])) == pytest.approx(19.7, abs=0.01)
    assert gpl_to_ugl(float(cs.y0[schema.slice("bound_methanethiol")][0])) == pytest.approx(
        1.4, abs=0.01
    )
    for name in ("bound_h2s_release", "bound_methanethiol_release"):
        assert name in cs.process_set
        assert not cs.process_set.is_enabled(name)


def test_explicit_zero_reservoirs_are_honoured_and_inert():
    # An explicit 0 is a documented deviation, not a silent one, and must produce the pre-D-135
    # model exactly.
    from fermentation.scenario.compile import compile_scenario
    from fermentation.scenario.schema import Scenario, TemperaturePoint

    cs = compile_scenario(
        Scenario(
            name="d135-zero",
            medium="wine",
            initial={
                "brix": 24.0,
                "yan_mgl": 200.0,
                "pitch_gpl": 0.25,
                "bound_h2s_ugl": 0.0,
                "bound_methanethiol_ugl": 0.0,
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    schema = cs.schema
    assert float(cs.y0[schema.slice("bound_h2s")][0]) == 0.0
    assert float(cs.y0[schema.slice("bound_methanethiol")][0]) == 0.0
