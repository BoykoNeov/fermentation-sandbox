"""Discrete winemaking interventions — the verb registry at the compile boundary (D-36).

``Scenario.interventions`` was declared since Milestone 1 but never consumed: the compile
seam turned only the temperature schedule into events. This activates the timeline of
winemaking verbs. Each verb compiles at the vocabulary boundary (decision D-3) into an
opaque :class:`~fermentation.runtime.schedule.ScheduledEvent`, and the verb-agnostic driver
(``simulate_scheduled``, D-35) segments-and-restarts around it, booking each state jump as an
:class:`~fermentation.runtime.schedule.ExternalFlow` for the conservation ledger.

This module pins the first verb, ``add_dap``:

* the dose lands exactly on ``N`` and is booked as one external flow (nothing else moves);
* nitrogen conserves across the jump — ``final == initial + Σ flows`` (the crown-jewel ledger
  the D-35 driver built the ledger for), and carbon is untouched;
* the emergent H₂S headline the *static* D-29 lever could not produce — a mid-ferment dose
  momentarily closes the inverse H₂S gate while sugar (hence flux) is still present, so the
  production rate drops right after the dose and net cumulative H₂S falls;
* the vocabulary discipline: unknown verbs, out-of-window days, and bad params raise loudly;
* isolability — a scenario with no interventions (and no ramp) is byte-for-byte a plain run.
"""

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics.carbon_routing import ESTER_SPECS, FUSEL_SPECS
from fermentation.core.tiers import Tier
from fermentation.runtime import simulate
from fermentation.runtime.schedule import ScheduledTrajectory
from fermentation.scenario import (
    Intervention,
    Scenario,
    TemperaturePoint,
    amino_acid_dose_nitrogen_mgl,
    compile_scenario,
)
from fermentation.units.convert import mgl_to_gpl
from fermentation.validation.conservation import total_carbon, total_nitrogen

_DAP_N_FRACTION = 0.2121  # exact (NH4)2HPO4 stoichiometry (additions.yaml)
#: g phosphoric-acid-equivalent phosphate per g DAP — M(H3PO4)/M((NH4)2HPO4), D-210.
_DAP_PHOSPHATE_FRACTION = 0.74206
_COPPER_H2S_BINDING = 0.536  # g H2S bound / g Cu, stoichiometric CuS 1:1 (additions.yaml, D-44)


def _wine(
    interventions: list[Intervention],
    *,
    yan_mgl: float = 80.0,
    days: float = 14.0,
    celsius: float = 20.0,
    anchored: bool = False,
) -> Scenario:
    # anchored (default False, byte-for-byte the pre-D-210 helper) seeds acids and an `initial_ph`,
    # which is what opts the charge balance in (D-179). It matters to `add_dap` since D-210: the
    # dose's two CHARGE writes ride that gate, so an unanchored wine gets the nitrogen jump D-36
    # shipped for the H2S gate and nothing in the balance.
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if anchored:
        initial |= {"tartaric_gpl": 6.0, "malic_gpl": 3.0, "initial_ph": 3.40}
    return Scenario(
        name="dap-test",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        interventions=interventions,
        duration_days=days,
    )


def _dap(day: float, dap_gpl: float) -> Intervention:
    return Intervention(day=day, action="add_dap", params={"dap_gpl": dap_gpl})


# -- the dose lands on N and is booked as one external flow --------------------


def test_add_dap_writes_both_ions_and_nothing_else():
    """The dose is a SALT: nitrogen, its phosphate counter-anion, and the pool's charge (D-210).

    This test used to assert ``count_nonzero(everything except N) == 0`` — "the whole injection
    is on N". That assert was correct about the code and wrong about the compound, and widening
    what the verb writes is exactly the kind of change it existed to catch. Now it pins all three
    writes individually, so the same protection survives at the wider scope.

    On an ANCHORED wine, because the two charge writes ride D-179's gate — see
    ``test_an_unanchored_dose_puts_nothing_in_the_charge_balance`` for the other side.
    """
    cs = compile_scenario(_wine([_dap(2.0, 0.4)], anchored=True))
    traj = cs.run()

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "add_dap@2d"
    assert flow.time_h == pytest.approx(48.0)  # 2 days in canonical hours

    schema = cs.schema
    expected_n = 0.4 * _DAP_N_FRACTION
    assert flow.delta[schema.slice("N")][0] == pytest.approx(expected_n, rel=1e-12)
    # The phosphate arrives with it, as phosphoric-acid-equivalent mass (D-210).
    assert flow.delta[schema.slice("phosphate")][0] == pytest.approx(
        0.4 * _DAP_PHOSPHATE_FRACTION, rel=1e-12
    )
    # And the pool's mean charge is re-mixed toward pure ammonium.
    excess = flow.delta[schema.slice(acidbase.NITROGEN_CHARGE_EXCESS_KEY)][0]
    assert excess > 0.0
    written = {"N", "phosphate", acidbase.NITROGEN_CHARGE_EXCESS_KEY}
    others = np.delete(flow.delta, [schema.slice(name).start for name in written])
    assert np.count_nonzero(others) == 0


def test_an_unanchored_dose_puts_nothing_in_the_charge_balance():
    """A nutrient addition is not pH information, and must not become it (D-210, D-179).

    Found in review after the first green suite, and it was two faults rather than one. The
    nitrogen half rides ``charge_balance_is_populated``; ``phosphate`` is an ordinary acid slot
    that ``_totals_molar`` reads unconditionally. So on a wine with no ``initial_ph`` and no dosed
    acids — an empty balance, and the shape of the Palma 2012 benchmark — the dose booked the
    salt's ANION with its CATION suppressed, *and* opened the gate by writing an acid slot, which
    switched the whole D-209 nitrogen term on for a beverage whose pH no scenario supplied. That
    was worth pH 3.10 → **4.53** at the dose instant and **−0.647 pH** at the end.

    Both writes are now atomic and gated together, so this wine's pH is what it was before D-210
    and a dose cannot open the gate.
    """
    schema = compile_scenario(_wine([_dap(2.0, 1.1)])).schema
    undosed_cs = compile_scenario(_wine([]))
    dosed_cs = compile_scenario(_wine([_dap(2.0, 1.1)]))
    undosed = np.asarray(undosed_cs.run().y, float)[:, -1]
    dosed = np.asarray(dosed_cs.run().y, float)[:, -1]

    # The gate is shut and the dose has not opened it.
    assert not acidbase.charge_balance_is_populated(undosed, schema)
    assert not acidbase.charge_balance_is_populated(dosed, schema)
    assert float(dosed[schema.slice("phosphate")][0]) == 0.0
    assert float(dosed[schema.slice(acidbase.NITROGEN_CHARGE_EXCESS_KEY)][0]) == 0.0
    assert acidbase.nitrogen_charge_molar(dosed, schema, dosed_cs.param_values) == 0.0

    # The nitrogen jump D-36 shipped for the H2S gate is untouched, so pH differs only by the
    # kinetic consequences of the extra nitrogen — not by ~0.6 pH of half-booked salt.
    ph_undosed = acidbase.ph_of_state(undosed, schema, undosed_cs.param_values)
    ph_dosed = acidbase.ph_of_state(dosed, schema, dosed_cs.param_values)
    assert abs(ph_dosed - ph_undosed) < 0.02, (
        f"an unanchored wine's pH moved {ph_dosed - ph_undosed:+.4f} on a DAP dose; the charge "
        "writes must ride D-179's gate, or a nutrient addition becomes pH information"
    )


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("dap_nitrogen_fraction", 0.2121),
        ("dap_phosphate_fraction", _DAP_PHOSPHATE_FRACTION),
        ("dap_nitrogen_charge", 1.0),
    ],
)
def test_the_dap_stoichiometry_params_are_validated_and_exact(name: str, value: float):
    """All three are formula, not measurement — VALIDATED with a zero-width band (D-36, D-210).

    ``dap_nitrogen_charge`` = 1.0 belongs here rather than in ``acidbase.yaml`` beside the two
    composition averages precisely because of this contrast: wine's 0.3245 and beer's 0.1772 are
    DERIVED from composition tables and carry real bands, while the charge of an ammonium ion is
    the compound's formula.
    """
    cs = compile_scenario(_wine([_dap(2.0, 0.4)]))
    param = cs.parameters[name]
    assert param.tier is Tier.VALIDATED
    assert param.value == pytest.approx(value)
    # exact stoichiometry ⇒ zero-width band, never swept by the ensemble
    assert param.uncertainty.low == param.value == param.uncertainty.high


# -- conservation across the jump: final == initial + Σ flows (crown jewel) ----


def test_add_dap_conserves_nitrogen_across_the_jump():
    cs = compile_scenario(_wine([_dap(3.0, 0.5)]))
    traj = cs.run()
    schema = cs.schema
    f_n = cs.param_values["biomass_N_fraction"]  # the value growth actually uses (D-14)
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)

    injected = sum(n_of(flow.delta) for flow in traj.external_flows)
    assert injected == pytest.approx(0.5 * _DAP_N_FRACTION, rel=1e-12)

    # The continuous ODE closes within each segment, so the run-wide identity is
    # final == initial + Σ(external injections). Nitrogen is conserved to machine precision.
    n_initial = n_of(cs.y0)
    n_final = n_of(traj.y[:, -1])
    assert n_final == pytest.approx(n_initial + injected, abs=1e-9)


def test_add_dap_leaves_carbon_untouched():
    cs = compile_scenario(_wine([_dap(3.0, 0.5)]))
    traj = cs.run()
    schema = cs.schema
    f_c = cs.param_values["biomass_C_fraction"]
    c_of = total_carbon(schema, biomass_carbon_fraction=f_c)

    # DAP is nitrogen plus phosphate and neither carries carbon (D-210 put the phosphate in the
    # charge balance, and phosphoric acid has no carbon atom — it is absent from MOLAR_MASS /
    # CARBON_ATOMS for that reason), so every flow is still carbon-free and the single-run
    # carbon balance closes with no ledger correction term.
    assert all(c_of(flow.delta) == pytest.approx(0.0, abs=1e-15) for flow in traj.external_flows)
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0), abs=1e-6)


# -- the emergent H2S headline (a timing effect a static dose cannot produce) --


def _run_pair() -> tuple[ScheduledTrajectory, ScheduledTrajectory]:
    grid = np.linspace(0.0, 14.0 * 24.0, 1400)
    undosed = compile_scenario(_wine([])).run(t_eval=grid)
    dosed = compile_scenario(_wine([_dap(2.0, 0.4)])).run(t_eval=grid)
    return undosed, dosed


def _rate_at(rate: np.ndarray, t_axis: np.ndarray, day: float) -> float:
    return float(rate[np.argmin(np.abs(t_axis - day * 24.0))])


def _produced_h2s(traj: ScheduledTrajectory) -> np.ndarray:
    # Cumulative H₂S produced = residual (h2s) + swept-to-gas (h2s_gas). The D-42 CO2-stripping
    # sink holds the residual pool at a quasi-steady µg/L level, so the PRODUCTION RATE is the
    # gradient of the produced SUM, not of the residual pool alone.
    return np.asarray(traj.series("h2s")) + np.asarray(traj.series("h2s_gas"))


def test_mid_ferment_dap_dose_drops_the_h2s_production_rate():
    undosed, dosed = _run_pair()
    # The dosed run has an extra grid point (the inserted day-2 breakpoint), so each run's
    # production rate is differenced against its OWN time axis.
    rate_u = np.gradient(_produced_h2s(undosed), undosed.t)
    rate_d = np.gradient(_produced_h2s(dosed), dosed.t)

    # Just before the day-2 dose the two runs are identical (same N history); just after, the
    # restored N re-represses the inverse gate K_h2s_n/(K_h2s_n+N) while sugar (the flux the
    # gate multiplies) is still present ⇒ the dosed production rate falls below the undosed one.
    # Pre-dose the two runs are the same physics; the tolerance is loose because the dosed run
    # restarts BDF at day 2, so its adaptive step sequence on [0, 48h] is not byte-identical to
    # the undosed run's single [0, 336h] integration (a mesh artifact, not a physical difference).
    assert _rate_at(rate_d, dosed.t, 1.8) == pytest.approx(
        _rate_at(rate_u, undosed.t, 1.8), rel=1e-3
    )  # ~identical pre-dose
    assert _rate_at(rate_d, dosed.t, 2.1) < 0.75 * _rate_at(rate_u, undosed.t, 2.1)  # gate closes


def test_dap_dose_lowers_net_cumulative_h2s():
    undosed, dosed = _run_pair()
    # Net over the whole ferment the gate closure dominates the competing extra-biomass flux
    # (more N ⇒ more growth ⇒ more flux later), so a DAP addition reduces total H₂S — the
    # realistic direction (DAP is the standard H₂S-management lever). This is EMERGENT, not
    # imposed: the model has no explicit "DAP lowers H₂S" term (decision D-36 / D-29). Read
    # cumulative PRODUCED (h2s + h2s_gas), since the D-42 stripping sink now splits produced H₂S
    # between the residual and headspace pools.
    assert _produced_h2s(dosed)[-1] < _produced_h2s(undosed)[-1]


# -- add_so2: lands on so2_total, rides neither elemental ledger ---------------


def _so2(day: float, so2_mgl: float) -> Intervention:
    return Intervention(day=day, action="add_so2", params={"so2_mgl": so2_mgl})


def test_add_so2_lands_on_so2_total_and_books_one_flow():
    cs = compile_scenario(_wine([_so2(5.0, 40.0)]))
    traj = cs.run()

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "add_so2@5d"
    schema = cs.schema
    # 40 mg/L → 0.040 g/L on so2_total; nothing else moves.
    assert flow.delta[schema.slice("so2_total")][0] == pytest.approx(0.040, rel=1e-12)
    others = np.delete(flow.delta, schema.slice("so2_total").start)
    assert np.count_nonzero(others) == 0


def test_add_so2_perturbs_neither_carbon_nor_nitrogen():
    cs = compile_scenario(_wine([_so2(5.0, 50.0)]))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])

    # SO₂ carries neither element, so the dose contributes nothing to either ledger and both
    # single-run balances still close with no correction term (unlike the DAP nitrogen jump).
    assert all(c_of(f.delta) == pytest.approx(0.0, abs=1e-15) for f in traj.external_flows)
    assert all(n_of(f.delta) == pytest.approx(0.0, abs=1e-15) for f in traj.external_flows)
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0), abs=1e-6)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0), abs=1e-9)


def test_add_so2_raises_the_molecular_so2_readout():
    from fermentation.analysis import molecular_so2_series

    cs = compile_scenario(_wine([_so2(5.0, 60.0)]))
    traj = cs.run()
    params = cs.param_values
    mol = molecular_so2_series(traj.as_trajectory(), params)

    # Before the day-5 dose there is no SO₂ (so2_total ≡ 0 ⇒ molecular ≡ 0); after it, the
    # readout is positive — the dose feeds the free/molecular partition from that time forward.
    before = mol[np.argmin(np.abs(traj.t - 4.5 * 24.0))]
    after = mol[np.argmin(np.abs(traj.t - 6.0 * 24.0))]
    assert before == pytest.approx(0.0, abs=1e-15)
    assert after > 0.0


def test_add_so2_requires_an_so2_slot():
    beer = Scenario(
        name="no-so2",
        medium="beer",
        initial={
            "glucose_gpl": 90.0,
            "maltose_gpl": 140.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 2.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[_so2(3.0, 40.0)],
        duration_days=10.0,
    )
    with pytest.raises(ValueError, match="needs a 'so2_total' slot"):
        compile_scenario(beer)


# -- rack: remove settled lees, leave the wine ---------------------------------


def _rack(day: float, fraction: float) -> Intervention:
    return Intervention(day=day, action="rack", params={"fraction": fraction})


def test_rack_removes_settled_lees_and_leaves_the_wine_untouched():
    # Autolysis opted in so debris (the second lees pool) is non-empty at racking.
    sc = Scenario(
        name="rack-test",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 150.0, "pitch_gpl": 0.25, "autolysis_rate_per_h": 0.002},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[_rack(15.0, 0.9)],
        duration_days=20.0,
    )
    cs = compile_scenario(sc)
    traj = cs.run()
    schema = cs.schema

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "rack@15d"

    # The settled pools drop (negative delta), and 90% of each is removed.
    for name in ("X_dead", "debris"):
        sl = schema.slice(name)
        assert flow.delta[sl][0] < 0.0
    # Everything that stays with the racked-off liquid is untouched (zero delta): viable biomass,
    # sugar, ethanol, YAN, glycerol, byproducts, acids, SO₂.
    for name in (
        "X", "S", "E", "N", "Gly", "Byp",
        *(spec.pool for spec in ESTER_SPECS),  # the three D-96 ester pools
        *(spec.pool for spec in FUSEL_SPECS),  # the five D-99 higher-alcohol pools
        "tartaric", "so2_total",
    ):  # fmt: skip
        assert np.all(flow.delta[schema.slice(name)] == 0.0)


def test_rack_fraction_out_of_range_raises():
    with pytest.raises(ValueError, match=r"fraction must be in \[0, 1\]"):
        compile_scenario(_wine([_rack(10.0, 1.5)]))


def test_dap_dose_and_rack_conserve_carbon_and_nitrogen_across_all_jumps():
    # The crown-jewel ledger the D-35 external-flow machinery was built for: a run with an
    # injection (DAP, +N) AND a removal (rack, −C/−N) still satisfies the run-wide identity
    # final == initial + Σ external_flows for BOTH elements, to machine precision.
    sc = Scenario(
        name="dap-and-rack",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 100.0, "pitch_gpl": 0.25, "autolysis_rate_per_h": 0.002},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[_dap(2.0, 0.4), _rack(15.0, 0.8)],
        duration_days=20.0,
    )
    cs = compile_scenario(sc)
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    assert len(traj.external_flows) == 2  # one dose, one rack

    for quantity, tol in ((c_of, 1e-6), (n_of, 1e-9)):
        injected = sum(quantity(f.delta) for f in traj.external_flows)
        initial = quantity(cs.y0)
        final = quantity(traj.y[:, -1])
        assert final == pytest.approx(initial + injected, abs=tol)
    # And the ledger is not trivially zero: the rack removed carbon, the DAP added nitrogen.
    assert sum(c_of(f.delta) for f in traj.external_flows) < 0.0
    assert sum(n_of(f.delta) for f in traj.external_flows) != pytest.approx(0.0, abs=1e-6)


# -- rack removes O. oeni off the lees — both viable X_mlf and settled X_mlf_dead (D-39) --


def _mlf_so2_rack(rack_fraction: float) -> Scenario:
    # Co-inoculate O. oeni with amino acids (so X_mlf grows), SO₂ at day 5 to kill some bacteria
    # (so X_mlf_dead accumulates while viable X_mlf remains), then a day-7 rack draws both off.
    initial: dict[str, float] = {
        "brix": 22.0,
        "yan_mgl": 200.0,
        "pitch_gpl": 0.25,
        "malic_gpl": 3.0,
        "initial_ph": 3.5,
        "mlf_pitch_gpl": 0.2,
        "amino_acids_gpl": 1.0,
    }  # fmt: skip
    # D-244: ``yan_mgl`` is the must's TOTAL assimilable nitrogen and the amino-acid dose is
    # carved OUT of it. This fixture was authored when the two channels ADDED, so it declares
    # the sum -- which leaves its pitch state bit-for-bit and moves only the point Coleman's
    # yield fit is evaluated at, which is the defect D-243 found.
    initial["yan_mgl"] += amino_acid_dose_nitrogen_mgl(initial)
    return Scenario(
        name="mlf-rack",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=22.0)],
        interventions=[
            Intervention(day=5.0, action="add_so2", params={"so2_mgl": 40.0}),
            _rack(7.0, rack_fraction),
        ],
        duration_days=14.0,
    )


def test_rack_removes_both_oenococcus_pools_and_leaves_viable_yeast():
    # decision D-39: racking draws O. oeni off the lees — BOTH viable X_mlf and settled X_mlf_dead —
    # the physical twin of the SO₂ kill (the deferred D-31 "rack early ⇒ diacetyl locked in" lever).
    # Viable YEAST X is left untouched: it ferments in suspension, so a rack leaves it working.
    cs = compile_scenario(_mlf_so2_rack(0.8))
    traj = cs.run()
    schema = cs.schema
    rack_flow = next(f for f in traj.external_flows if f.label == "rack@7d")

    # both bacterial pools are non-empty at the rack (viable remnant + SO₂-killed dead) and drop
    x_mlf, x_mlf_dead = traj.series("X_mlf"), traj.series("X_mlf_dead")
    i_pre = int(np.searchsorted(traj.t, 7.0 * 24.0)) - 1
    assert x_mlf[i_pre] > 0.0 and x_mlf_dead[i_pre] > 0.0
    for name in ("X_mlf", "X_mlf_dead"):
        assert rack_flow.delta[schema.slice(name)][0] < 0.0  # bacteria drawn off with the lees
    # viable yeast and every dissolved species stay with the racked-off liquid (untouched)
    for name in ("X", "S", "E", "N", "malic", "so2_total"):
        assert np.all(rack_flow.delta[schema.slice(name)] == 0.0)


def test_rack_on_an_mlf_run_conserves_carbon_and_nitrogen():
    # Racking O. oeni off the lees removes biomass carbon AND nitrogen (both X_mlf pools weighted
    # since D-38), booked as a negative external flow — so the run-wide identity still closes to
    # machine precision for BOTH elements (SO₂ carries neither, so only the rack moves the ledger).
    cs = compile_scenario(_mlf_so2_rack(0.8))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])

    for quantity, tol in ((c_of, 1e-6), (n_of, 1e-9)):
        injected = sum(quantity(f.delta) for f in traj.external_flows)
        assert quantity(traj.y[:, -1]) == pytest.approx(quantity(cs.y0) + injected, abs=tol)
    # the rack removed biomass carbon AND nitrogen — a non-trivial ledger move, not a wash
    rack_flow = next(f for f in traj.external_flows if f.label == "rack@7d")
    assert c_of(rack_flow.delta) < 0.0
    assert n_of(rack_flow.delta) < 0.0


# -- pitch_mlf: mutate the catalyst + reconfigure the Process set --------------


_MLF_NAMES = (
    "malolactic_conversion",
    "malolactic_citrate_metabolism",
    "oenococcus_diacetyl_reduction",
)


def _mlf_wine(interventions: list[Intervention]) -> Scenario:
    return Scenario(
        name="mlf-test",
        medium="wine",
        initial={
            "brix": 22.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.25,
            "malic_gpl": 3.0,
            "initial_ph": 3.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=22.0)],
        interventions=interventions,
        duration_days=25.0,
    )


def _pitch(day: float, pitch_gpl: float = 0.1) -> Intervention:
    return Intervention(day=day, action="pitch_mlf", params={"pitch_gpl": pitch_gpl})


def test_pitch_mlf_mutates_catalyst_and_enables_exactly_the_gated_set():
    cs = compile_scenario(_mlf_wine([_pitch(1.0, 0.1)]))
    # Unpitched at compile: the gated Processes are all disabled.
    assert all(not cs.process_set.is_enabled(name) for name in _MLF_NAMES)

    traj = cs.run()
    # The pitch enabled exactly the gated set (symmetric with the compile-time disable).
    assert all(cs.process_set.is_enabled(name) for name in _MLF_NAMES)

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "pitch_mlf@1d"
    schema = cs.schema
    # The whole jump is the X_mlf catalyst dose; nothing else moves.
    assert flow.delta[schema.slice("X_mlf")][0] == pytest.approx(0.1, rel=1e-12)
    others = np.delete(flow.delta, schema.slice("X_mlf").start)
    assert np.count_nonzero(others) == 0


def test_pitch_mlf_flow_carries_bacterial_biomass_carbon_and_nitrogen():
    # Since the MLF-growth beat (decision D-38) promoted X_mlf from an inert carbon-/
    # nitrogen-free catalyst to real bacterial biomass, X_mlf is weighted in both ledgers, so a
    # pitch_mlf dose adds bacterial-biomass carbon AND nitrogen — booked on the external flow (the
    # run-wide ledger final == initial + Σ flows still closes; that is checked in test_schedule).
    cs = compile_scenario(_mlf_wine([_pitch(1.0, 0.1)]))
    traj = cs.run()
    schema = cs.schema
    f_c = cs.param_values["biomass_C_fraction"]
    f_n = cs.param_values["biomass_N_fraction"]
    c_of = total_carbon(schema, biomass_carbon_fraction=f_c)
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)
    flow = traj.external_flows[0]
    # The 0.1 g/L X_mlf dose carries exactly its biomass carbon/nitrogen and nothing else.
    assert c_of(flow.delta) == pytest.approx(0.1 * f_c, rel=1e-12)
    assert n_of(flow.delta) == pytest.approx(0.1 * f_n, rel=1e-12)


def test_early_pitch_converts_malic_but_post_af_pitch_stalls():
    # The honest headline (decisions D-23, D-31): pitch timing is now a scenario choice, and only
    # an early/co-inoculation pitch completes — a post-AF pitch lands past the Luong ethanol wall,
    # so the environmental gate keeps malate conversion near zero. EMERGENT from the gate, not
    # imposed: the same Process, only the pitch time differs.
    early = compile_scenario(_mlf_wine([_pitch(1.0)])).run()
    late = compile_scenario(_mlf_wine([_pitch(15.0)])).run()
    unpitched = compile_scenario(_mlf_wine([])).run()

    malic0 = 3.0
    early_malic = early.series("malic")[-1]
    late_malic = late.series("malic")[-1]
    assert unpitched.series("malic")[-1] == pytest.approx(malic0, abs=1e-9)  # inert, unchanged
    assert early_malic < 0.5 * malic0  # early pitch converts most of the malate
    assert late_malic > 0.9 * malic0  # post-AF pitch stalls at the ethanol wall
    assert early_malic < late_malic


def test_mid_run_pitch_drags_malolactic_tiers_for_the_whole_run():
    # tier travels (D-35): the Processes are enabled only from the breakpoint, but min-combining
    # the per-segment tier maps reports the touched slots speculative for the WHOLE trajectory.
    unpitched = compile_scenario(_mlf_wine([])).run()
    pitched = compile_scenario(_mlf_wine([_pitch(1.0)])).run()
    for name in ("malic", "lactic"):
        assert unpitched.tier_map[name] is Tier.VALIDATED  # disabled ⇒ inert ⇒ validated
        assert pitched.tier_map[name] is Tier.SPECULATIVE  # enabled mid-run ⇒ speculative run-wide


def test_pitch_mlf_requires_an_x_mlf_slot():
    beer = Scenario(
        name="no-mlf",
        medium="beer",
        initial={
            "glucose_gpl": 90.0,
            "maltose_gpl": 140.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 2.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[_pitch(3.0)],
        duration_days=10.0,
    )
    with pytest.raises(ValueError, match="needs an 'X_mlf' slot"):
        compile_scenario(beer)


# -- vocabulary discipline: loud failures --------------------------------------


def test_unknown_action_raises():
    with pytest.raises(ValueError, match="unknown intervention action 'add_potion'"):
        compile_scenario(_wine([Intervention(day=1.0, action="add_potion")]))


def test_intervention_at_or_after_duration_raises():
    with pytest.raises(ValueError, match="at or beyond the run duration"):
        compile_scenario(_wine([_dap(14.0, 0.4)], days=14.0))


def test_add_dap_missing_required_param_raises():
    with pytest.raises(ValueError, match="missing required param 'dap_gpl'"):
        compile_scenario(_wine([Intervention(day=2.0, action="add_dap")]))


def test_add_dap_unknown_param_raises():
    bad = Intervention(day=2.0, action="add_dap", params={"dap_gpl": 0.4, "phosphate_gpl": 0.1})
    with pytest.raises(ValueError, match="unknown param"):
        compile_scenario(_wine([bad]))


def test_add_dap_negative_dose_raises():
    with pytest.raises(ValueError, match="must be >= 0"):
        compile_scenario(_wine([_dap(2.0, -0.4)]))


# -- the merge: a temperature ramp AND an intervention on the one driver -------


def test_ramp_events_and_intervention_events_merge_on_one_driver():
    # The realistic winemaking scenario — a fermentation temperature schedule AND a dosing
    # timeline — and the literal embodiment of the D-35→D-36 claim that both ride the same
    # driver. Every other test populates only one side of ``events``; this exercises the merge.
    from fermentation.units.convert import celsius_to_kelvin

    sc = Scenario(
        name="ramp-and-dose",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 100.0, "pitch_gpl": 0.25},
        # 16 → 24 °C over days 0–4, then down to 18 °C by day 12: two interior slope changes.
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=16.0),
            TemperaturePoint(day=4.0, celsius=24.0),
            TemperaturePoint(day=12.0, celsius=18.0),
        ],
        interventions=[_dap(2.0, 0.4)],
        duration_days=14.0,
    )
    cs = compile_scenario(sc)
    # Both sides are genuinely populated: ramp slope-change events AND the DAP dose.
    labels = [e.label for e in cs.events]
    assert any(lbl.startswith("temperature_ramp@") for lbl in labels)
    assert "add_dap@2d" in labels

    traj = cs.run()
    assert traj.success
    schema = cs.schema

    # The dose (which touches only N) did not perturb the temperature path: after the day-2 dose,
    # still on the first ramp leg (t < day 4), T follows the analytic line 16 + slope·t exactly
    # (BDF integrates a constant slope to round-off, D-35). Assert at the actual grid time so the
    # coarse default grid does not turn the ramp offset into a spurious failure.
    idx = np.argmin(np.abs(traj.t - 3.0 * 24.0))
    t_h = traj.t[idx]
    assert 0.0 <= t_h < 4.0 * 24.0
    expected_t = celsius_to_kelvin(16.0 + (24.0 - 16.0) / (4.0 * 24.0) * t_h)
    assert traj.series("T")[idx] == pytest.approx(expected_t, abs=1e-6)

    # The DAP flow still lands on N, and the nitrogen ledger still closes across the jump.
    dap_flows = [f for f in traj.external_flows if f.label == "add_dap@2d"]
    assert len(dap_flows) == 1
    dap_n = dap_flows[0].delta[schema.slice("N")][0]
    assert dap_n == pytest.approx(0.4 * _DAP_N_FRACTION, rel=1e-12)
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    injected = sum(n_of(f.delta) for f in traj.external_flows)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0) + injected, abs=1e-9)


# -- isolability: no interventions ⇒ byte-for-byte a plain run -----------------


def test_no_interventions_is_byte_for_byte_plain_simulate():
    # A wine scenario with no interventions and a flat (single-knot) schedule compiles to an
    # empty event tuple, so CompiledScenario.run() is a single simulate() call with identical
    # arguments — the scheduling layer adds nothing (the D-35 isolability discipline).
    cs = compile_scenario(_wine([]))
    assert cs.events == ()

    grid = np.linspace(0.0, cs.t_span_h[1], 200)
    scheduled = cs.run(t_eval=grid)
    plain = simulate(
        cs.process_set,
        cs.param_values,
        cs.y0,
        cs.t_span_h,
        param_tiers=cs.parameters.tier_map(),
        t_eval=grid,
    )
    assert np.array_equal(scheduled.y, plain.y)
    assert scheduled.external_flows == ()


# -- add_copper: copper-fine dissolved H₂S out of the wine (decision D-44) -----
#
# The remediation half of the reductive-fault beat. Copper (Cu²⁺) precipitates dissolved sulfide
# as insoluble CuS (1:1 mol), so add_copper removes min(h2s_present, capacity) from the h2s pool,
# where capacity = (copper g/L) · copper_h2s_binding. H₂S is carbon-free, so — like add_so2 — the
# removal rides neither elemental ledger. This is the standard fix for the un-stripped autolytic
# reduction the D-44 AutolyticHydrogenSulfide source builds up post-dryness.


def _copper(day: float, copper_mgl: float) -> Intervention:
    return Intervention(day=day, action="add_copper", params={"copper_mgl": copper_mgl})


def _copper_event(cs, day: float):
    return next(e for e in cs.events if e.label == f"add_copper@{day:g}d")


def test_add_copper_removes_stoichiometric_h2s():
    # Verb-level, applied to a hand-built state so the removal is exact and not confounded by
    # concurrent production: copper binds H₂S 1:1 as CuS, so a dose removes min(present, capacity).
    cs = compile_scenario(_wine([_copper(5.0, 0.05)]))
    ev = _copper_event(cs, 5.0)
    schema = cs.schema
    h2s = schema.slice("h2s")
    capacity = mgl_to_gpl(0.05) * cs.param_values["copper_h2s_binding"]  # ~26.8 µg/L
    assert capacity == pytest.approx(mgl_to_gpl(0.05) * _COPPER_H2S_BINDING, rel=1e-12)

    # present ≫ capacity ⇒ partial removal of EXACTLY the capacity (copper is the limiting reagent)
    high = cs.y0.copy()
    high[h2s] = 1.0e-3  # 1 mg/L dissolved, far above capacity
    assert float(ev.mutate(schema, high, cs.param_values)[h2s][0]) == pytest.approx(
        1.0e-3 - capacity, rel=1e-12
    )

    # present < capacity ⇒ ALL of it removed, and no more (the pool floors at 0, not negative)
    low = cs.y0.copy()
    low[h2s] = 1.0e-5  # 10 µg/L, below capacity
    assert float(ev.mutate(schema, low, cs.param_values)[h2s][0]) == pytest.approx(0.0, abs=1e-18)


def test_add_copper_does_not_strip_a_negative_undershoot():
    # Guard symmetry with the kinetics' ≥0 clamps: if a solver undershoot left h2s < 0, copper must
    # not "remove" a negative amount (which would ADD H₂S). The clamp leaves it untouched.
    cs = compile_scenario(_wine([_copper(5.0, 0.5)]))
    ev = _copper_event(cs, 5.0)
    schema = cs.schema
    h2s = schema.slice("h2s")
    y = cs.y0.copy()
    y[h2s] = -1.0e-9
    assert float(ev.mutate(schema, y, cs.param_values)[h2s][0]) == pytest.approx(-1.0e-9, rel=1e-12)


def test_add_copper_lands_on_h2s_and_copper_and_books_one_flow():
    cs = compile_scenario(_wine([_copper(5.0, 0.5)]))
    traj = cs.run()
    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "add_copper@5d"
    schema = cs.schema
    # The h2s delta is ≤ 0 (copper never ADDS H₂S), the copper delta is > 0 (D-191: the dose stays
    # in the wine), and the bound_h2s delta is the EXACT NEGATIVE of the h2s one (D-193: the
    # sulfide is complexed, not destroyed). Those are the only three slots this verb may touch on a
    # default wine — methanethiol is empty with autolysis off, so its removal, and hence its own
    # transfer, is exactly 0 here.
    assert flow.delta[schema.slice("h2s")][0] <= 0.0
    assert flow.delta[schema.slice("copper")][0] > 0.0
    assert flow.delta[schema.slice("bound_h2s")][0] == pytest.approx(
        -flow.delta[schema.slice("h2s")][0], rel=1e-15
    )
    # ...and the transfer is not vacuous: some H₂S really was bound on this wine.
    assert flow.delta[schema.slice("bound_h2s")][0] > 0.0
    others = np.delete(
        flow.delta,
        [
            schema.slice("h2s").start,
            schema.slice("copper").start,
            schema.slice("bound_h2s").start,
        ],
    )
    assert np.count_nonzero(others) == 0


def test_add_copper_leaves_the_sourced_residual_share_in_the_wine():
    """D-191: fining does not remove the copper — >95 % of the dose stays dissolved.

    Pins the arithmetic against the PARAMETER rather than a literal, so a re-sourcing of
    ``copper_fining_residual_fraction`` moves this test with it; the separate band test below is
    what pins the sourced edge itself.
    """
    cs = compile_scenario(_wine([_copper(5.0, 0.5)]))
    traj = cs.run()
    schema = cs.schema
    retained = cs.param_values["copper_fining_residual_fraction"]
    background = cs.y0[schema.slice("copper")][0]
    dosed_gpl = 0.5e-3  # 0.5 mg/L Cu in canonical g/L

    # The slot ACCUMULATES onto the grape/must background rather than being set to the dose.
    assert traj.y[schema.slice("copper"), -1][0] == pytest.approx(
        background + retained * dosed_gpl, rel=1e-12
    )
    # And the background is copper_typical, so before the fining f_copper was exactly 1.
    assert background == pytest.approx(cs.param_values["copper_typical"], rel=1e-12)


def test_repeated_finings_accumulate_copper():
    """Two doses leave twice the copper: the verb credits, it does not set."""
    once = compile_scenario(_wine([_copper(5.0, 0.5)]))
    twice = compile_scenario(_wine([_copper(5.0, 0.5), _copper(6.0, 0.5)]))
    schema = once.schema
    background = once.y0[schema.slice("copper")][0]
    gained_once = once.run().y[schema.slice("copper"), -1][0] - background
    gained_twice = twice.run().y[schema.slice("copper"), -1][0] - background
    assert gained_twice == pytest.approx(2.0 * gained_once, rel=1e-12)


def test_add_copper_perturbs_neither_carbon_nor_nitrogen():
    cs = compile_scenario(_wine([_copper(5.0, 0.5)]))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    # On this DEFAULT wine (autolysis off) mercaptans ≡ 0, so copper binds only the carbon-free H₂S:
    # the removal contributes nothing to either ledger and both balances close with no correction
    # term (the add_so2 case). Since D-193 that holds on a THIOL-BEARING wine too — the mercaptide
    # stays in the wine, so its carbon transfers rather than leaves
    # (test_add_copper_mercaptan_removal_moves_no_carbon_out_of_the_wine).
    assert all(c_of(f.delta) == pytest.approx(0.0, abs=1e-15) for f in traj.external_flows)
    assert all(n_of(f.delta) == pytest.approx(0.0, abs=1e-15) for f in traj.external_flows)
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0), abs=1e-6)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0), abs=1e-9)


def _autolysis_wine(interventions: list[Intervention]) -> Scenario:
    # A wine with autolysis opted in, so AutolyticHydrogenSulfide builds an un-stripped reductive
    # H₂S residual post-dryness (decision D-44) — the fault the copper fining below remediates.
    return Scenario(
        name="reduction-cu",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 80.0,
            "pitch_gpl": 0.25,
            "autolysis_rate_per_h": 0.002,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=interventions,
        duration_days=40.0,
    )


def test_add_copper_clears_the_autolytic_reduction():
    # THE remediation headline (D-44): a copper fining late in a reductive (autolysing) wine
    # precipitates the accumulated H₂S residual and drops the h2s pool toward zero — copper clears
    # reduction. Compared against the same wine with NO copper (identical up to the dose).
    grid = np.linspace(0.0, 40.0 * 24.0, 401)
    no_cu = compile_scenario(_autolysis_wine([])).run(t_eval=grid)
    fined = compile_scenario(_autolysis_wine([_copper(38.0, 0.5)])).run(t_eval=grid)
    h2s_nocu = np.asarray(no_cu.series("h2s"))
    h2s_cu = np.asarray(fined.series("h2s"))

    # the reductive residual really did build up (not a µg/L trace) — this is the fault to fix
    assert h2s_nocu[-1] > 1.0e-5  # > 10 µg/L, well above the stripped-default residual
    # before the day-38 fining the two runs share the same H₂S history
    i_before = int(np.argmin(np.abs(fined.t / 24.0 - 37.0)))
    assert h2s_cu[i_before] == pytest.approx(h2s_nocu[i_before], rel=1e-2)
    # after it, the fined run's residual collapses — copper removed ~all the dissolved H₂S
    assert h2s_cu[-1] < h2s_cu[i_before]
    assert h2s_cu[-1] < 0.25 * h2s_nocu[-1]


def test_add_copper_missing_param_raises():
    with pytest.raises(ValueError, match="missing required param 'copper_mgl'"):
        compile_scenario(_wine([Intervention(day=5.0, action="add_copper", params={})]))


# -- add_copper also binds mercaptans (H₂S-first partition, carbon-bearing removal, D-45) -------

_COPPER_MERC_BINDING = 1.514  # g MeSH / g Cu, stoichiometric Cu(SR)₂ 1:2 (additions.yaml, D-45)


def test_add_copper_binds_h2s_first_then_mercaptans():
    # Verb-level partition (D-45): copper binds H₂S first (CuS far more insoluble), then binds
    # mercaptans with the LEFTOVER copper (Cu(SR)₂, 1 Cu : 2 thiol). Applied to a hand-built state
    # so both removals are exact.
    cs = compile_scenario(_wine([_copper(5.0, 0.5)]))
    ev = _copper_event(cs, 5.0)
    schema = cs.schema
    h2s, merc = schema.slice("h2s"), schema.slice("methanethiol")
    copper_gpl = mgl_to_gpl(0.5)
    b_h2s = cs.param_values["copper_h2s_binding"]
    b_merc = cs.param_values["copper_mercaptan_binding"]
    assert b_merc == pytest.approx(_COPPER_MERC_BINDING, rel=1e-12)

    # H₂S below its share (fully removed), leaving copper over for mercaptans
    y = cs.y0.copy()
    y[h2s] = 1.0e-5  # 10 µg/L H₂S, well under the ~268 µg/L H₂S capacity
    y[merc] = 1.0e-3  # 1 mg/L mercaptans, above the leftover-copper capacity
    out = ev.mutate(schema, y, cs.param_values)
    assert float(out[h2s][0]) == pytest.approx(0.0, abs=1e-18)  # all H₂S bound first
    copper_left = copper_gpl - 1.0e-5 / b_h2s  # copper after the (tiny) H₂S bind
    expected_merc_removed = copper_left * b_merc
    assert float(out[merc][0]) == pytest.approx(1.0e-3 - expected_merc_removed, rel=1e-9)


def test_add_copper_h2s_consumes_copper_before_mercaptans():
    # If H₂S alone soaks up ALL the copper, none is left for mercaptans (the affinity order).
    cs = compile_scenario(_wine([_copper(5.0, 0.001)]))  # tiny dose: ~0.5 µg/L H₂S capacity
    ev = _copper_event(cs, 5.0)
    schema = cs.schema
    h2s, merc = schema.slice("h2s"), schema.slice("methanethiol")
    y = cs.y0.copy()
    y[h2s] = 1.0e-3  # 1 mg/L H₂S ≫ capacity ⇒ soaks all copper
    y[merc] = 1.0e-4
    out = ev.mutate(schema, y, cs.param_values)
    assert float(out[h2s][0]) < 1.0e-3  # some H₂S removed
    assert float(out[merc][0]) == pytest.approx(1.0e-4, abs=1e-18)  # mercaptans untouched


def test_add_copper_mercaptan_removal_moves_no_carbon_out_of_the_wine():
    # D-193 INVERTS THIS TEST'S CLAIM. Through D-192 it asserted `c_of(flow.delta) < 0.0` —
    # mercaptans carry carbon, so binding them was booked as carbon LEAVING the wine. That rested
    # on the precipitation mechanism D-191 found retracted: the mercaptide stays dispersed, so the
    # carbon is transferred to bound_methanethiol at the identical weight and the fining is
    # carbon-neutral to machine precision.
    #
    # WHY THE ASSERT IS A PIN AT ZERO AND NOT A SIGN CHECK. The old `< 0.0` does not survive as a
    # guard in either direction: after the change the flow is float noise (~1e-23) whose SIGN is
    # arbitrary — measured +6.7e-23 on this scenario and −1.2e-23 on the D-193 probe's — so it
    # would pass or fail on rounding, not on behaviour. Hence a tolerance pin plus a positive
    # control that the transfer actually happened.
    sc = Scenario(
        name="reduction-cu-carbon",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25, "autolysis_rate_per_h": 0.002},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[_copper(38.0, 0.5)],
        duration_days=40.0,
    )
    cs = compile_scenario(sc)
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    # POSITIVE CONTROL: thiol really was bound on this wine, and it landed in the bonded pool
    # mole-for-mole. Without this the carbon zero below would also hold on a wine with no thiol.
    removed = -flow.delta[schema.slice("methanethiol")][0]
    assert removed > 0.0
    assert flow.delta[schema.slice("bound_methanethiol")][0] == pytest.approx(removed, rel=1e-15)
    # The fining moves NEITHER element: the thiol's carbon changed binding state inside the wine,
    # and thiols are nitrogen-free as before.
    assert c_of(flow.delta) == pytest.approx(0.0, abs=1e-18)
    assert n_of(flow.delta) == pytest.approx(0.0, abs=1e-15)
    # crown-jewel: carbon closes across the jump — now with a zero correction term, not an excused
    # one, which is the strongest form of the identity rather than a weaker one.
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0) + c_of(flow.delta), abs=1e-6)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0), abs=1e-9)


def test_add_copper_clears_both_h2s_and_mercaptans():
    # THE full remediation headline (D-45): a copper fining of a reductive wine drops BOTH the
    # accumulated H₂S residual and the mercaptans toward zero — copper clears reduction.
    grid = np.linspace(0.0, 40.0 * 24.0, 401)
    no_cu = compile_scenario(_autolysis_wine([])).run(t_eval=grid)
    fined = compile_scenario(_autolysis_wine([_copper(38.0, 0.5)])).run(t_eval=grid)
    for pool in ("h2s", "methanethiol"):
        built = np.asarray(no_cu.series(pool))
        cleared = np.asarray(fined.series(pool))
        assert built[-1] > 1.0e-6, pool  # the fault really built up
        assert cleared[-1] < 0.25 * built[-1], pool  # copper cleared ~all of it


# -- add_acid: dose a charge-active organic acid; pH drops, TA rises (decision D-65) -----------
#
# The §3.3 acidulation verb, general over the D-18 charge-active acids (tartaric/malic/lactic).
# The dose is the PURE acid (own protons, no counter-cation), so it lands on the acid slot but NOT
# on cation_charge — the D-18 charge balance re-solves the same fixed strong cation against more
# anion, so pH DROPS and TA RISES *emergently* (potassium bitartrate, which adds a counter-cation,
# would be different). Each acid carries carbon, so the dose is a POSITIVE carbon external flow
# (the add_dap +N precedent, opposite sign to the copper mercaptan −C removal) and nitrogen-free.


def _acid(day: float, acid: str, gpl: float) -> Intervention:
    return Intervention(day=day, action="add_acid", params={"acid": acid, "gpl": gpl})


def test_add_acid_lands_on_the_named_acid_slot_and_leaves_cation_untouched():
    cs = compile_scenario(_wine([_acid(3.0, "tartaric", 1.5)]))
    traj = cs.run()
    schema = cs.schema

    assert len(traj.external_flows) == 1
    flow = traj.external_flows[0]
    assert flow.label == "add_acid@3d"
    # The whole dose lands on the tartaric slot; every other slot delta is zero — crucially
    # cation_charge is UNTOUCHED (the load-bearing modelling choice that makes pH drop).
    assert flow.delta[schema.slice("tartaric")][0] == pytest.approx(1.5, rel=1e-12)
    others = np.delete(flow.delta, schema.slice("tartaric").start)
    assert np.count_nonzero(others) == 0
    assert flow.delta[schema.slice("cation_charge")][0] == 0.0


def test_add_acid_is_general_over_the_charge_active_acids():
    # Owner chose the general add_acid {acid, gpl} over a tartaric-only verb (decision D-65): any
    # D-18 charge-active acid slot can be dosed. Malic starts at 0 in this must but its slot exists.
    cs = compile_scenario(_wine([_acid(3.0, "malic", 2.0)]))
    traj = cs.run()
    schema = cs.schema
    flow = traj.external_flows[0]
    assert flow.delta[schema.slice("malic")][0] == pytest.approx(2.0, rel=1e-12)
    others = np.delete(flow.delta, schema.slice("malic").start)
    assert np.count_nonzero(others) == 0


def test_add_acid_books_a_positive_carbon_flow_and_no_nitrogen():
    # Tartaric (C4H6O6) is carbon-weighted in total_carbon (D-18), so the dose ADDS carbon — a
    # positive external flow (the mirror of the copper mercaptan −C removal). Nitrogen-free. The
    # crown-jewel identity final == initial + Σ flows still closes for both elements.
    cs = compile_scenario(_wine([_acid(3.0, "tartaric", 1.5)]))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    flow = traj.external_flows[0]

    # the flow carries exactly the dosed acid's carbon (1.5 g/L × tartaric carbon fraction) and no N
    assert c_of(flow.delta) == pytest.approx(1.5 * carbon_mass_fraction("tartaric_acid"), rel=1e-12)
    assert c_of(flow.delta) > 0.0
    assert n_of(flow.delta) == pytest.approx(0.0, abs=1e-15)
    # crown-jewel: carbon closes across the jump once the positive external flow is counted
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0) + c_of(flow.delta), abs=1e-6)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0), abs=1e-9)


def _ph_ta(traj, params, schema) -> tuple[float, float]:
    y_final = traj.y[:, -1]
    return (
        acidbase.ph_of_state(y_final, schema, params),
        acidbase.titratable_acidity(y_final, schema, params),
    )


def test_add_acid_lowers_ph_and_raises_ta():
    # THE headline (D-65): a tartaric addition acidifies — pH DROPS and TA RISES vs an otherwise
    # identical undosed wine. Emergent from the D-18 charge balance (more anion, same cation), NOT
    # scripted. Direction + a loose band only: acidbase.py claims directional/slope fidelity for the
    # concentration-based apparent pKa, so a tight pH-delta would over-claim. pH does not feed back
    # into the yeast kinetics (D-18), so S/E/X are identical dosed-vs-undosed, only tartaric moves.
    undosed = compile_scenario(_wine([], days=14.0)).run()
    dosed = compile_scenario(_wine([_acid(1.0, "tartaric", 2.0)], days=14.0)).run()

    base = compile_scenario(_wine([]))
    params, schema = base.param_values, base.schema
    ph_u, ta_u = _ph_ta(undosed, params, schema)
    ph_d, ta_d = _ph_ta(dosed, params, schema)

    assert ph_d < ph_u  # acidification lowers pH
    assert ph_u - ph_d < 1.0  # but within a sane band (a ~2 g/L tartaric bump, not a collapse)
    assert ta_d > ta_u  # and raises titratable acidity


def test_add_acid_moves_no_tier():
    # Unlike pitch_mlf, add_acid enables no Processes and touches an inert slot, so no tier moves:
    # the acid slots stay VALIDATED (no Process touches them) and the derived pH tier is unchanged.
    undosed = compile_scenario(_wine([])).run()
    dosed = compile_scenario(_wine([_acid(1.0, "tartaric", 2.0)])).run()
    for name in ("tartaric", "malic", "S", "E"):
        assert dosed.tier_map[name] is undosed.tier_map[name]


def test_add_acid_unknown_acid_raises():
    with pytest.raises(ValueError, match="unknown acid 'citric'"):
        compile_scenario(_wine([_acid(3.0, "citric", 1.0)]))


def test_add_acid_on_beer_raises_wine_only():
    beer = Scenario(
        name="no-acid",
        medium="beer",
        initial={
            "glucose_gpl": 90.0,
            "maltose_gpl": 140.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 2.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[_acid(3.0, "tartaric", 1.0)],
        duration_days=10.0,
    )
    with pytest.raises(ValueError, match="needs a 'tartaric' slot"):
        compile_scenario(beer)


def test_add_acid_missing_and_negative_params_raise():
    with pytest.raises(ValueError, match="missing required param 'acid'"):
        compile_scenario(_wine([Intervention(day=3.0, action="add_acid", params={"gpl": 1.0})]))
    with pytest.raises(ValueError, match="missing required param 'gpl'"):
        compile_scenario(
            _wine([Intervention(day=3.0, action="add_acid", params={"acid": "tartaric"})])
        )
    with pytest.raises(ValueError, match="must be >= 0"):
        compile_scenario(_wine([_acid(3.0, "tartaric", -1.0)]))


# -- add_sugar: chaptalize (dose sucrose, inverted to fermentable hexose) (decision D-65) -------
#
# Sucrose is dosed by mass and inverted AT THE DOSE (invertase is fast vs the ferment) to hexose-
# equivalent via the exact sucrose_inversion_mass_ratio (~1.0526; the +5.26 % is hydrolysis water).
# The hexose lands on the fermentable sugar slot (wine's lumped S; beer's glucose specifically).
# More sugar ⇒ higher finished ethanol once it ferments out (emergent). Carbon is conserved through
# inversion, so the flow books exactly the sucrose carbon (a positive flow, add_dap precedent).

_SUCROSE_INVERSION_RATIO = 1.0526  # g hexose / g sucrose (additions.yaml, D-65)


def _sugar(day: float, sugar_gpl: float) -> Intervention:
    return Intervention(day=day, action="add_sugar", params={"sugar_gpl": sugar_gpl})


def test_add_sugar_inverts_sucrose_onto_the_hexose_slot():
    cs = compile_scenario(_wine([_sugar(2.0, 20.0)]))
    frac = cs.parameters["sucrose_inversion_mass_ratio"]
    assert frac.tier is Tier.VALIDATED  # exact stoichiometry
    assert frac.uncertainty.low == frac.value == frac.uncertainty.high  # zero-width band

    traj = cs.run()
    schema = cs.schema
    flow = traj.external_flows[0]
    assert flow.label == "add_sugar@2d"
    # 20 g/L sucrose → 20 × 1.0526 g/L hexose on the single lumped S slot; nothing else moves.
    expected_hexose = 20.0 * _SUCROSE_INVERSION_RATIO
    assert flow.delta[schema.slice("S")][0] == pytest.approx(expected_hexose, rel=1e-4)
    others = np.delete(flow.delta, schema.slice("S").start)
    assert np.count_nonzero(others) == 0


def test_add_sugar_books_positive_carbon_and_no_nitrogen():
    # The inverted hexose carries carbon (weighted at the glucose fraction), so the dose is a
    # positive carbon external flow; sugar is nitrogen-free. Carbon closes across the jump. Because
    # inversion is carbon-conserving (water carries none), the booked carbon equals the sucrose's.
    cs = compile_scenario(_wine([_sugar(2.0, 20.0)]))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])
    flow = traj.external_flows[0]

    # the flow books EXACTLY the inverted hexose's carbon (20 g/L sucrose × ratio × glucose
    # fraction) — which equals the dosed sucrose's carbon, since inversion is carbon-conserving
    # (hydrolysis water carries none): the docstring's "books exactly the sucrose carbon" claim.
    expected_c = 20.0 * _SUCROSE_INVERSION_RATIO * carbon_mass_fraction("glucose")
    assert c_of(flow.delta) == pytest.approx(expected_c, rel=1e-4)
    assert c_of(flow.delta) > 0.0
    assert n_of(flow.delta) == pytest.approx(0.0, abs=1e-15)
    # crown-jewel: carbon closes across the jump once the positive external flow is counted
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0) + c_of(flow.delta), abs=1e-6)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0), abs=1e-9)


def test_add_sugar_raises_finished_ethanol():
    # THE headline: chaptalizing lifts the finished ethanol once the added sugar ferments out.
    # Dosed early in a modest must and run long, so both runs finish dry; the dosed run then has
    # strictly more ethanol (the extra hexose fermented) — emergent, no explicit ABV term.
    grid = np.linspace(0.0, 20.0 * 24.0, 401)
    undosed = compile_scenario(_wine([], yan_mgl=200.0, days=20.0)).run(t_eval=grid)
    dosed = compile_scenario(_wine([_sugar(1.0, 30.0)], yan_mgl=200.0, days=20.0)).run(t_eval=grid)

    # both finished dry (added sugar fermented out, not stranded), and the dosed run made more EtOH
    assert dosed.series("S")[-1] < 5.0  # g/L residual — dry
    assert dosed.series("E")[-1] > undosed.series("E")[-1]


def test_add_sugar_on_beer_targets_glucose_only():
    # Owner chose wine + beer scope (decision D-65). Beer's S is a 3-vector (glucose/maltose/
    # maltotriose); add_sugar must land the inverted hexose on the GLUCOSE slot alone, never
    # broadcast — fructose from the inversion lumps as glucose-equivalent (isomers, exact carbon).
    beer = Scenario(
        name="prime",
        medium="beer",
        initial={
            "glucose_gpl": 90.0,
            "maltose_gpl": 140.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 2.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=[_sugar(1.0, 10.0)],
        duration_days=10.0,
    )
    cs = compile_scenario(beer)
    traj = cs.run()
    schema = cs.schema
    flow = traj.external_flows[0]
    s_start = schema.slice("S").start
    # the inverted hexose lands on the first (glucose) sugar slot; maltose/maltotriose are untouched
    assert flow.delta[s_start] == pytest.approx(10.0 * _SUCROSE_INVERSION_RATIO, rel=1e-4)
    assert flow.delta[s_start + 1] == 0.0
    assert flow.delta[s_start + 2] == 0.0


def test_add_sugar_missing_and_negative_params_raise():
    with pytest.raises(ValueError, match="missing required param 'sugar_gpl'"):
        compile_scenario(_wine([Intervention(day=2.0, action="add_sugar", params={})]))
    with pytest.raises(ValueError, match="must be >= 0"):
        compile_scenario(_wine([_sugar(2.0, -10.0)]))


# -- set_ph: the pH a beverage AGES at (decision D-186, closing D-150's open item) ------------

# The ferment is long enough to reach dryness, so the aging breakpoint sits where every aging
# Process actually runs and Byp has stopped moving; the tail is short because nothing here
# asserts an aging MAGNITUDE — only which pH the aging segment runs at.
_PH_FERMENT_DAYS = 20.0
_PH_AGING_DAYS = 20.0
_ANCHOR_PH = 3.60


def _ph_wine(
    interventions: list[Intervention],
    *,
    initial_ph: float | None = _ANCHOR_PH,
) -> Scenario:
    # tartaric/malic are seeded either way (they are their own scenario keys); initial_ph is what
    # opts the pH system in, and omitting it is how the gate below is exercised.
    initial: dict[str, float] = {
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.25,
        "tartaric_gpl": 6.0,
        "malic_gpl": 3.0,
    }
    if initial_ph is not None:
        initial["initial_ph"] = initial_ph
    return Scenario(
        name="set-ph-test",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=interventions,
        duration_days=_PH_FERMENT_DAYS + _PH_AGING_DAYS,
    )


def _ph_beer(interventions: list[Intervention], *, initial_ph: float | None = 5.50) -> Scenario:
    initial: dict[str, float] = {
        "glucose_gpl": 15.0,
        "maltose_gpl": 70.0,
        "maltotriose_gpl": 15.0,
        "yan_mgl": 150.0,
        "pitch_gpl": 0.5,
    }
    if initial_ph is not None:
        initial["initial_ph"] = initial_ph
    return Scenario(
        name="set-ph-beer-test",
        medium="beer",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=18.0)],
        interventions=interventions,
        duration_days=_PH_FERMENT_DAYS + _PH_AGING_DAYS,
    )


def _set_ph(day: float, ph: float) -> Intervention:
    return Intervention(day=day, action="set_ph", params={"ph": ph})


def _aging_at(ph: float | None) -> list[Intervention]:
    """``begin_aging`` at the ferment/aging boundary, optionally re-anchored to ``ph`` there."""
    events = [Intervention(day=_PH_FERMENT_DAYS, action="begin_aging")]
    if ph is not None:
        events.append(_set_ph(_PH_FERMENT_DAYS, ph))
    return events


#: A grid that straddles the aging breakpoint by a hair on each side, so "the pH the aging
#: segment runs at" is read from a real sample point rather than inferred. The breakpoint hour
#: itself is deliberately NOT probed: an event time is a segment BOUNDARY, and which side of the
#: mutation a sample there reports is an implementation detail of the driver, not a claim this
#: verb should be pinned to. The driver also augments any ``t_eval`` with its own breakpoints,
#: so the probes are located by TIME below rather than by index.
_BOUNDARY_H = _PH_FERMENT_DAYS * 24.0
_END_H = (_PH_FERMENT_DAYS + _PH_AGING_DAYS) * 24.0
_PROBE_H = np.array([0.0, _BOUNDARY_H - 0.01, _BOUNDARY_H + 0.01, _END_H])


def _ph_after_boundary(cs, traj) -> float:
    """pH at the first sample strictly after the aging breakpoint — the pH aging runs at."""
    i = int(np.argmax(traj.t > _BOUNDARY_H))
    return acidbase.ph_of_state(traj.y[:, i], cs.schema, cs.param_values)


def _ph_before_boundary(cs, traj) -> float:
    """pH at the last sample strictly before the aging breakpoint (the untouched half)."""
    i = int(np.max(np.flatnonzero(traj.t < _BOUNDARY_H)))
    return acidbase.ph_of_state(traj.y[:, i], cs.schema, cs.param_values)


def test_cation_charge_for_ph_inverts_ph_of_state_on_a_fermented_state():
    # The core claim, at the state level: the new inverse is exact against the forward solver on a
    # state that has FERMENTED — Byp accumulated and the dissolved-CO2 term saturated, the two
    # things the compile-seam anchor gets to assume away. Round trip is closed-form on the cation
    # side, so it is pinned far tighter than solve_ph's own 1e-10 brentq tolerance would need.
    cs = compile_scenario(_ph_wine([]))
    traj = cs.run()
    y_dry = traj.y[:, -1]
    schema, params = cs.schema, cs.param_values
    assert traj.series("Byp")[-1] > 0.0  # the state really is fermented
    assert acidbase.dissolved_co2_molar(y_dry, schema, params) > 0.0

    cations = []
    for target in (3.20, 3.40, 3.60, 3.80):
        y = y_dry.copy()
        y[schema.slice("cation_charge")] = acidbase.cation_charge_for_ph(y, schema, params, target)
        assert acidbase.ph_of_state(y, schema, params) == pytest.approx(target, abs=1e-9)
        cations.append(float(y[schema.slice("cation_charge")][0]))
    # more base ⇒ higher pH: the inverse is monotone, which is what makes it a real anchor
    assert cations == sorted(cations)


def test_set_ph_to_the_states_own_ph_leaves_the_cation_where_it_was():
    # Inertness stated as a round trip rather than as a byte-for-byte run: asking for the pH the
    # state already has must return the cation already there. A one-sided inverse (say, one that
    # silently dropped the carbonic term the forward solver reads) would still pass every
    # "pH moved the right way" test and fail HERE.
    cs = compile_scenario(_ph_wine([]))
    traj = cs.run()
    schema, params = cs.schema, cs.param_values
    y_dry = traj.y[:, -1]
    own_ph = acidbase.ph_of_state(y_dry, schema, params)
    before = float(y_dry[schema.slice("cation_charge")][0])
    after = acidbase.cation_charge_for_ph(y_dry, schema, params, own_ph)
    assert after == pytest.approx(before, rel=1e-9)


def test_the_aging_ph_is_the_one_asked_for_not_the_one_the_ferment_left():
    # THE headline (D-186). Without the verb, the aging segment runs at whatever pH the ferment
    # happened to leave — measurably BELOW the anchor, because Byp accumulates against a cation
    # frozen at pitch (D-150's finding). With it, the aging segment runs at the stated pH.
    drift_cs = compile_scenario(_ph_wine(_aging_at(None)))
    drift = drift_cs.run(t_eval=_PROBE_H)
    drifted_ph = _ph_after_boundary(drift_cs, drift)
    assert drifted_ph < _ANCHOR_PH - 0.01, "the drift this verb exists to fix is not present"

    set_cs = compile_scenario(_ph_wine(_aging_at(_ANCHOR_PH)))
    set_traj = set_cs.run(t_eval=_PROBE_H)
    # just before the breakpoint the run is untouched; just after it, the wine sits where asked
    assert _ph_before_boundary(set_cs, set_traj) == pytest.approx(
        _ph_before_boundary(drift_cs, drift), abs=1e-6
    )
    assert _ph_after_boundary(set_cs, set_traj) == pytest.approx(_ANCHOR_PH, abs=1e-5)


def test_two_aging_arms_span_exactly_what_they_asked_for_where_initial_ph_falls_short():
    # The defect, paired with its baseline (the D-150 measurement, re-taken): asking for a 0.35-unit
    # pH SPAN via initial_ph delivers less than 0.35 at the aging boundary, because each arm drifts
    # by its own amount. Asking via set_ph delivers 0.35. Both arms are measured the same way, so
    # the comparison is of the two mechanisms and not of two protocols.
    low, high = 3.30, 3.65
    span_asked = high - low

    def boundary_ph(scenario) -> float:
        cs = compile_scenario(scenario)
        return _ph_after_boundary(cs, cs.run(t_eval=_PROBE_H))

    anchored = [boundary_ph(_ph_wine(_aging_at(None), initial_ph=p)) for p in (low, high)]
    assert anchored[1] - anchored[0] < span_asked - 0.005, "initial_ph delivered the span asked for"

    at_ph = [boundary_ph(_ph_wine(_aging_at(p))) for p in (low, high)]
    assert at_ph[0] == pytest.approx(low, abs=1e-5)
    assert at_ph[1] == pytest.approx(high, abs=1e-5)
    assert at_ph[1] - at_ph[0] == pytest.approx(span_asked, abs=1e-5)


def test_set_ph_books_a_flow_that_is_carbon_and_nitrogen_free():
    # Unlike add_acid (+C) and add_copper (−C), a strong cation is charge, not matter the ledgers
    # weigh: cation_charge carries weight 0 in total_carbon and is absent from total_nitrogen. So
    # the jump is booked like every other mutation and contributes exactly nothing to either
    # element — final == initial + Σ flows closes with a zero term, not with an excused one.
    cs = compile_scenario(_ph_wine(_aging_at(3.40)))
    traj = cs.run()
    schema = cs.schema
    c_of = total_carbon(schema, biomass_carbon_fraction=cs.param_values["biomass_C_fraction"])
    n_of = total_nitrogen(schema, biomass_nitrogen_fraction=cs.param_values["biomass_N_fraction"])

    flow = next(f for f in traj.external_flows if f.label.startswith("set_ph"))
    # the cation really moved — otherwise the zero below would be vacuous
    assert flow.delta[schema.slice("cation_charge")][0] != 0.0
    assert c_of(flow.delta) == pytest.approx(0.0, abs=1e-15)
    assert n_of(flow.delta) == pytest.approx(0.0, abs=1e-15)
    # and nothing else in the state was touched
    others = np.delete(flow.delta, schema.slice("cation_charge").start)
    assert np.count_nonzero(others) == 0
    assert c_of(traj.y[:, -1]) == pytest.approx(c_of(cs.y0), abs=1e-6)
    assert n_of(traj.y[:, -1]) == pytest.approx(n_of(cs.y0), abs=1e-9)


def test_set_ph_moves_no_tier():
    # pH's tier is already the PLAUSIBLE-floored pKa tier (D-18) and this writes a STATE slot, not
    # a parameter — so the reported tier map must be identical to the same run without the verb.
    plain = compile_scenario(_ph_wine(_aging_at(None))).run()
    with_set = compile_scenario(_ph_wine(_aging_at(3.40))).run()
    assert dict(with_set.tier_map) == dict(plain.tier_map)


def test_a_target_below_the_acid_loads_own_ph_raises_and_names_the_floor():
    # The one check that CANNOT be made at compile (the reachable range is a property of the state
    # at the moment of the adjustment). It surfaces from the event application with the verb's
    # label and the achievable floor, not as a bare traceback out of the solver.
    cs = compile_scenario(_ph_wine(_aging_at(1.50)))
    with pytest.raises(ValueError, match=r"set_ph@20d.*intrinsic pH"):
        cs.run()


def test_set_ph_needs_the_scenario_to_have_opted_into_the_ph_system():
    # Both media, because the reason differs: beer without initial_ph has an EMPTY acid load, wine
    # has a populated-but-unanchored one. The gate is the same and it fires at compile.
    with pytest.raises(ValueError, match="needs the pH system"):
        compile_scenario(_ph_wine(_aging_at(3.40), initial_ph=None))
    with pytest.raises(ValueError, match="needs the pH system"):
        compile_scenario(_ph_beer(_aging_at(4.20), initial_ph=None))


def test_beer_can_set_its_aging_ph_too():
    # Medium-agnostic by construction — it reads the medium's own acid registry (D-179) — so the
    # beer path is exercised rather than assumed. Beer's finished pH is a PREDICTION (D-180), and
    # this does not change that: it re-anchors from the day it fires, downstream of the prediction.
    cs = compile_scenario(_ph_beer(_aging_at(4.20)))
    traj = cs.run(t_eval=_PROBE_H)
    assert _ph_after_boundary(cs, traj) == pytest.approx(4.20, abs=1e-5)


def test_set_ph_vocabulary_errors_are_loud():
    with pytest.raises(ValueError, match="missing required param 'ph'"):
        compile_scenario(_ph_wine([Intervention(day=2.0, action="set_ph", params={})]))
    with pytest.raises(ValueError, match=r"unknown param\(s\) \['target'\]"):
        compile_scenario(
            _ph_wine([Intervention(day=2.0, action="set_ph", params={"ph": 3.4, "target": 3.4})])
        )
    with pytest.raises(ValueError, match="must be >= 0"):
        compile_scenario(_ph_wine([_set_ph(2.0, -1.0)]))
    with pytest.raises(ValueError, match=r"strictly inside \(0, 14\)"):
        compile_scenario(_ph_wine([_set_ph(2.0, 14.0)]))
    with pytest.raises(ValueError, match=r"strictly inside \(0, 14\)"):
        compile_scenario(_ph_wine([_set_ph(2.0, 0.0)]))


# -- D-210: DAP is a SALT, and both of its ions are in the charge balance -------


_D210_DOSE_DAY = 3.0
_D210_DAYS = 14.0


def _dap_charge_wine(dap_gpl: float, *, dose_day: float = _D210_DOSE_DAY) -> Scenario:
    """An ANCHORED wine (the D-210 arm), optionally dosed.

    ``initial_ph`` is what opts the charge balance in (D-179), so an unanchored wine would show
    none of this — nitrogen is not pH information on its own.
    """
    interventions = (
        [Intervention(day=dose_day, action="add_dap", params={"dap_gpl": dap_gpl})]
        if dap_gpl
        else []
    )
    return Scenario(
        name="d210-dap-charge",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.25,
            "tartaric_gpl": 6.0,
            "malic_gpl": 3.0,
            "initial_ph": 3.40,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=interventions,
        duration_days=_D210_DAYS,
    )


def _end_state(dap_gpl: float, *, dose_day: float = _D210_DOSE_DAY):
    cs = compile_scenario(_dap_charge_wine(dap_gpl, dose_day=dose_day))
    traj = cs.run()
    return cs, np.asarray(traj.y, float)[:, -1]


def _ph_with(cs, state, *, ammonium: bool = True, phosphate: bool = True) -> float:
    """pH of one state with either D-210 half switched off IN THE STATE.

    Each half is isolated by zeroing its OWN slot, which is what makes the attribution honest:
    the two have opposite signs at the dose instant, so a single "with/without the beat" number
    would credit each with the other's cancellation.
    """
    probe = state.copy()
    if not ammonium:
        probe[cs.schema.slice(acidbase.NITROGEN_CHARGE_EXCESS_KEY)] = 0.0
    if not phosphate:
        probe[cs.schema.slice("phosphate")] = 0.0
    return acidbase.ph_of_state(probe, cs.schema, cs.param_values)


def test_a_dap_dose_permanently_acidifies_a_wine_through_its_phosphate():
    """The endpoint claim, re-RUN and not re-scored (D-210).

    ``ph_of_state`` is read inside ``EsterHydrolysis``, ``EthylAcetateEsterification``,
    ``AcetaldehydeReduction``, ``PhenolicBrowning`` and both Strecker limbs, so a static
    re-score of a fixed trajectory is an instrument and only a re-run is the claim (D-209 §6).
    """
    cs_undosed, undosed = _end_state(0.0)
    ph_undosed = acidbase.ph_of_state(undosed, cs_undosed.schema, cs_undosed.param_values)

    cs, dosed = _end_state(0.3)
    ph_dosed = acidbase.ph_of_state(dosed, cs.schema, cs.param_values)

    # A cellar dose leaves the wine measurably more acidic than the same wine undosed. Before
    # D-210 this difference was ~2e-4 pH — the dose's kinetic side-effects only.
    assert ph_dosed < ph_undosed - 0.03
    assert ph_dosed == pytest.approx(3.1683, abs=2e-3)
    assert ph_undosed == pytest.approx(3.2123, abs=2e-3)


def test_the_phosphate_half_owns_the_endpoint_and_the_ammonium_half_the_excursion():
    """Attribution at the granularity that can be wrong — and it INVERTS the pre-registration.

    D-209 §8c predicted the dosed nitrogen's charge would matter because the acidification was
    "understated by ~3x on the dosed fraction". True of the *excursion*, false of the endpoint:
    on a DRY wine the dosed nitrogen is fully consumed either way, so the pool's charge returns
    to ~0 whatever it was and the ammonium half is worth exactly ZERO at the end. Everything
    permanent is the phosphate. Shipping the ammonium half alone would have credited it with
    +0.235 pH of dose-time rise that the phosphate cancels 71 % of.
    """
    cs, end = _end_state(1.1)  # Palma 2012's refeed dose
    both = acidbase.ph_of_state(end, cs.schema, cs.param_values)

    # The dosed nitrogen is gone, so its charge cannot be in the endpoint either way.
    assert float(end[cs.schema.slice("N")][0]) == pytest.approx(0.0, abs=1e-6)
    assert _ph_with(cs, end, ammonium=False) == pytest.approx(both, abs=1e-9)
    # The phosphate is still there and is the whole of it.
    assert _ph_with(cs, end, phosphate=False) - both == pytest.approx(0.162, abs=5e-3)


def test_the_ammonium_half_reaches_the_endpoint_when_nitrogen_is_left_standing():
    """The ammonium half is not endpoint-free in general — only on a wine that finishes dry.

    Dose late enough that the run's horizon ends with the supplement still in the pool (a stuck
    or sluggish ferment, a sweet wine, or simply a short horizon) and standing NH4+ is real
    cation charge: the wine reads 0.23 pH HIGHER than the pool-average charge said. That is the
    case a "the ammonium half does nothing" reading of the test above would get wrong.
    """
    cs, end = _end_state(1.1, dose_day=12.0)
    standing = float(end[cs.schema.slice("N")][0])
    assert standing > 0.2  # the supplement is still there at the horizon
    both = acidbase.ph_of_state(end, cs.schema, cs.param_values)
    assert both - _ph_with(cs, end, ammonium=False) == pytest.approx(0.228, abs=5e-3)


def test_the_dose_time_ph_rise_is_emergent_and_the_two_halves_oppose_there():
    """At the dose instant the salt is electroneutral as it ARRIVES (+2 ammonium, -2 hydrogen
    phosphate) and the balance gains ~+1.05 net cation per mole — which IS the ~1.05 protons the
    dosed HPO4(2-) takes up on its way to H2PO4(-) at wine pH. So the pH rise falls out of the
    balance, and nothing may add a proton-consumption term on top of it.

    Pinned at 1.1 g/L, where the halves are largest: ammonium +0.235, phosphate -0.166, net
    +0.070. A single with/without number would report only that net and hide both.
    """
    cs = compile_scenario(_dap_charge_wine(1.1))
    traj = cs.run(t_eval=np.linspace(0.0, _D210_DAYS * 24.0, 2001))
    y = np.asarray(traj.y, float)
    k = int(np.searchsorted(np.asarray(traj.t, float), _D210_DOSE_DAY * 24.0 + 1e-9))
    state = y[:, k]

    both = acidbase.ph_of_state(state, cs.schema, cs.param_values)
    ammonium_worth = both - _ph_with(cs, state, ammonium=False)
    phosphate_worth = both - _ph_with(cs, state, phosphate=False)
    assert ammonium_worth == pytest.approx(0.235, abs=8e-3)
    assert phosphate_worth == pytest.approx(-0.166, abs=8e-3)
    # Opposite signs, and the phosphate cancels most of the ammonium's bump.
    assert ammonium_worth > 0.0 > phosphate_worth
    assert abs(phosphate_worth) / ammonium_worth == pytest.approx(0.705, abs=0.05)


def test_an_undosed_wine_carries_the_new_slots_at_exactly_zero():
    """The isolability claim, in the only form checkable from inside (prime directive #3).

    Both D-210 slots default to 0.0 and only ``add_dap`` writes them, so an undosed wine's charge
    balance gains the exact float 0.0 — not a small number. ``nitrogen_charge_excess`` stores the
    EXCESS over the medium average rather than the mean charge itself precisely so that this is
    true without a sentinel: 0.0 is both the default and the correct undosed value.
    """
    cs, end = _end_state(0.0)
    assert float(end[cs.schema.slice("phosphate")][0]) == 0.0
    assert float(end[cs.schema.slice(acidbase.NITROGEN_CHARGE_EXCESS_KEY)][0]) == 0.0
    # And the pool's charge is exactly the medium average, i.e. the D-209 value unchanged.
    n_gpl = float(end[cs.schema.slice("N")][0])
    assert acidbase.nitrogen_charge_molar(end, cs.schema, cs.param_values) == pytest.approx(
        acidbase.nitrogen_charge_from_gpl(n_gpl, "wine", cs.param_values), rel=1e-12
    )


def _setph_after_dose(target_ph: float, *, set_day: float = 1.02) -> Scenario:
    return Scenario(
        name="d210-setph-after-dose",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 300.0,
            "pitch_gpl": 0.25,
            "tartaric_gpl": 2.0,
            "initial_ph": 3.40,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[
            Intervention(day=1.0, action="add_dap", params={"dap_gpl": 3.0}),
            Intervention(day=set_day, action="set_ph", params={"ph": target_ph}),
        ],
        duration_days=_D210_DAYS,
    )


def test_a_set_ph_after_a_big_dap_dose_names_the_nitrogen_and_not_the_acid_load():
    """The second way into ``cation_slot_after_nitrogen``, and a mis-diagnosis it exposed (D-210).

    Two separate faults, both found here. First, D-209's guard was reachable only from the compile
    seam, so its message said "Reduce yan_mgl, add acid, or raise initial_ph" — advice a mid-run
    ``set_ph`` caller cannot act on, against a scenario key they did not use. After a dose the pool
    carries up to 3x the average charge per mole N, so a target that was fine before the dose can
    demand a negative remainder: 3 g/L of DAP into a low-acid wine puts 0.047 mol+/L of nitrogen
    charge against the 0.045 a pH-4.5 target demands.

    Second and worse, ``_verb_set_ph`` caught bare ``ValueError`` and rewrote it as "target pH is
    below this state's intrinsic pH — the acid load alone holds it there". Every configuration
    scanned while writing this test reported that message, and in every one of them the acid load
    was NOT the cause. Hence ``NitrogenExceedsCationDemandError``: the two failures have different
    remedies, so they need different messages.
    """
    with pytest.raises(acidbase.NitrogenExceedsCationDemandError, match=r"set_ph=4\.5 needs"):
        compile_scenario(_setph_after_dose(4.5)).run()
    with pytest.raises(ValueError, match="pure ammonium"):
        compile_scenario(_setph_after_dose(4.5)).run()
    # It is still a ValueError, so a caller catching the broad type is unaffected.
    with pytest.raises(ValueError):
        compile_scenario(_setph_after_dose(4.5)).run()


def test_the_set_ph_floor_message_still_fires_where_the_acid_load_really_is_the_cause():
    """The floor branch is narrowed by D-210, not removed — the positive control for the test above.

    Once the supplement has been assimilated the pool's charge is ~0 and an unreachable target is
    unreachable for the ordinary D-186 reason. Both branches must stay live, or "it names the
    nitrogen" would be true by having only one branch left.

    A normally-acidified wine, dosed at a cellar rate and re-anchored on day 6 when its nitrogen
    is long gone: the floor is pH 2.1475 and a target below it reports the floor, not nitrogen.
    """
    scenario = Scenario(
        name="d210-setph-floor",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.25,
            "tartaric_gpl": 6.0,
            "malic_gpl": 3.0,
            "initial_ph": 3.40,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[
            Intervention(day=3.0, action="add_dap", params={"dap_gpl": 0.3}),
            Intervention(day=6.0, action="set_ph", params={"ph": 2.0}),
        ],
        duration_days=_D210_DAYS,
    )
    with pytest.raises(ValueError, match="is below this state's intrinsic pH") as caught:
        compile_scenario(scenario).run()
    assert not isinstance(caught.value, acidbase.NitrogenExceedsCationDemandError)


def test_the_compile_seam_guard_keeps_its_seed_shaped_remedy():
    """At compile no dose can have fired, so the remedy is the scenario's own keys (D-209/D-210)."""
    with pytest.raises(ValueError, match="Reduce yan_mgl"):
        compile_scenario(
            Scenario(
                name="d210-unphysical-seed",
                medium="wine",
                initial={
                    "brix": 24.0,
                    "yan_mgl": 400.0,
                    "pitch_gpl": 0.25,
                    "tartaric_gpl": 2.0,
                    "initial_ph": 2.90,
                },
                temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
                duration_days=_D210_DAYS,
            )
        )
