"""Tests for the acetaldehyde pathway Processes (decision D-27).

Acetaldehyde is the obligate intermediate on the *main* ethanol pathway, so — unlike the
side-pool aroma byproducts — it is modelled as a transient ethanol-carbon **buffer**, not a
draw from sugar: production borrows a C2 slice of ethanol and reduction returns it. That
de-lumps the uptake Process's already-complete sugar→ethanol step instead of adding a
parallel pathway (which would double-count and inflate ABV).

* **Production** (:class:`AcetaldehydeProduction`): flux-linked, temperature-flat; borrows
  carbon from ``E`` (``d[acetaldehyde] += r``, ``d[E] -= r·M_eth/M_acet``).
* **Reduction** (:class:`AcetaldehydeReduction`): enzymatic C2 → C2 to ethanol, gated on
  *viable* yeast with no flux term — so it clears acetaldehyde during the rest but stops when
  the yeast is crashed (stranding it, with the borrowed ethanol carbon un-returned).

The unit tests pin each Process's closed form, the borrow/return carbon accounting and the
guards; the acceptance section verifies the *emergent* early peak (produced during active
fermentation, then reduced to a low residual), the machine-precision carbon closure on a
compiled run, and — the buffer's headline guarantee — that the ``E`` endpoint reconverges to
the no-acetaldehyde core (to relative ~1e-8), so the §2.2 benchmarks are preserved to far
below any tolerance (a tiny ~1e-4 second-order path drift via the E→viability brake aside).
"""

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.core.chemistry import (
    M_ACETALDEHYDE,
    M_ETHANOL,
    M_MALIC,
    M_TARTARIC,
    carbon_mass_fraction,
)
from fermentation.core.kinetics import (
    AcetaldehydeProduction,
    AcetaldehydeReduction,
    SugarUptakeToEthanolCO2,
    arrhenius_factor,
)
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.units.convert import gpl_to_mgl, mgl_to_gpl
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon

#: Carbon fractions the borrow/return books against (mirror the chemistry constants).
_ACETALDEHYDE_C = carbon_mass_fraction("acetaldehyde")
_ETHANOL_C = carbon_mass_fraction("ethanol")
#: The mole-for-mole ethanol⇄acetaldehyde mass ratio used by both Processes.
_ETH_PER_ACET = M_ETHANOL / M_ACETALDEHYDE


@pytest.fixture
def store():
    # Wine kinetics PLUS the shared acetaldehyde constants (k_acetaldehyde, k_acet_reduction…).
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "acetaldehyde.yaml",
    )


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x: float = 2.0,
    s: float = 200.0,
    e: float = 50.0,
    n: float = 0.1,
    t: float = 293.15,
    acetaldehyde: float = 0.0,
) -> FloatArray:
    y = schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": t, "CO2": 0.0})
    y[schema.slice("acetaldehyde")] = acetaldehyde
    return y


# -- metadata -----------------------------------------------------------------


def test_production_metadata():
    p = AcetaldehydeProduction()
    assert p.name == "acetaldehyde_production"
    # Speculative: rate magnitude is an order-of-magnitude estimate.
    assert p.tier is Tier.SPECULATIVE
    # Touches ONLY its own pool and ethanol — the buffer borrows C2 from E, never S/CO2 (D-27).
    assert set(p.touches) == {"acetaldehyde", "E"}
    # ``k_acet_so2_induced`` is the D-48 SO₂-induced over-production coefficient (total-SO₂ driven).
    assert set(p.reads) == {"k_acetaldehyde", "K_sugar_uptake", "k_acet_so2_induced"}


def test_reduction_metadata():
    p = AcetaldehydeReduction()
    assert p.name == "acetaldehyde_reduction"
    assert p.tier is Tier.SPECULATIVE
    assert set(p.touches) == {"acetaldehyde", "E"}
    assert set(p.reads) == {"k_acet_reduction", "E_a_acet_reduction", "T_ref"}


# -- production closed form & guards ------------------------------------------


def test_production_matches_closed_form(params):
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = _wine_y0(schema, x=x, s=s)
    d = AcetaldehydeProduction().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    rate = params["k_acetaldehyde"] * flux
    assert schema.get(d, "acetaldehyde") == pytest.approx(rate)
    # Ethanol is decremented mole-for-mole (C2→C2): the carbon leaving E equals the carbon
    # entering acetaldehyde, so the borrow is carbon-exact per RHS.
    assert schema.get(d, "E") == pytest.approx(-rate * _ETH_PER_ACET)
    assert schema.get(d, "E") * _ETHANOL_C == pytest.approx(
        -schema.get(d, "acetaldehyde") * _ACETALDEHYDE_C
    )
    # Nothing else moves — no sugar draw, no CO2 (the whole point of the buffer model).
    for var in ("X", "S", "N", "CO2"):
        assert schema.get(d, var) == 0.0


def test_production_is_temperature_flat(params):
    # Documented v1 simplification (D-27): production carries NO Arrhenius factor — the
    # temperature dependence lives in the enzymatic reduction. Identical cold vs warm at flux.
    schema = wine_schema()
    cold = AcetaldehydeProduction().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = AcetaldehydeProduction().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(cold, "acetaldehyde") == pytest.approx(schema.get(warm, "acetaldehyde"))
    assert schema.get(cold, "acetaldehyde") > 0.0


def test_production_scales_with_fermentative_flux(params):
    # Coupled to the biomass-catalysed sugar flux (linear in X): twice the biomass ⇒ twice
    # the borrow.
    schema = wine_schema()
    r1 = AcetaldehydeProduction().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = AcetaldehydeProduction().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "acetaldehyde") == pytest.approx(2.0 * schema.get(r1, "acetaldehyde"))


def test_production_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    no_x = AcetaldehydeProduction().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)
    no_s = AcetaldehydeProduction().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)
    assert schema.get(no_x, "acetaldehyde") == 0.0
    assert schema.get(no_s, "acetaldehyde") == 0.0
    # No borrow ⇒ ethanol untouched too.
    assert schema.get(no_x, "E") == 0.0
    assert schema.get(no_s, "E") == 0.0


# -- reduction closed form & guards ------------------------------------------


def test_reduction_matches_closed_form(params):
    schema = wine_schema()
    x, acet = 2.0, 0.04
    y = _wine_y0(schema, x=x, acetaldehyde=acet)
    d = AcetaldehydeReduction().derivatives(0.0, y, schema, params)

    f_t = arrhenius_factor(293.15, params["E_a_acet_reduction"], params["T_ref"])
    loss = params["k_acet_reduction"] * x * f_t * acet
    assert schema.get(d, "acetaldehyde") == pytest.approx(-loss)
    # Ethanol gains mole-for-mole (C2→C2): the carbon returned to E equals the carbon lost
    # from acetaldehyde, so the return is carbon-exact per RHS.
    assert schema.get(d, "E") == pytest.approx(loss * _ETH_PER_ACET)
    assert schema.get(d, "E") * _ETHANOL_C == pytest.approx(
        -schema.get(d, "acetaldehyde") * _ACETALDEHYDE_C
    )
    # Never touches sugar/CO2/nitrogen.
    for var in ("S", "CO2", "N"):
        assert schema.get(d, var) == 0.0


def test_reduction_is_temperature_dependent(params):
    # Enzymatic ⇒ Arrhenius in T: warmer reduces faster (unlike production, which is flat).
    schema = wine_schema()
    cold = AcetaldehydeReduction().derivatives(
        0.0, _wine_y0(schema, t=283.15, acetaldehyde=0.04), schema, params
    )
    warm = AcetaldehydeReduction().derivatives(
        0.0, _wine_y0(schema, t=303.15, acetaldehyde=0.04), schema, params
    )
    warm_rate = abs(float(schema.get(warm, "acetaldehyde")))
    cold_rate = abs(float(schema.get(cold, "acetaldehyde")))
    assert warm_rate > cold_rate > 0.0


def test_reduction_gated_on_viable_yeast(params):
    # Gated on VIABLE X (not X_dead), no flux term: with no live yeast the acetaldehyde is
    # STRANDED — the borrowed ethanol carbon is never returned (the crash-at-peak behaviour).
    schema = wine_schema()
    stranded = AcetaldehydeReduction().derivatives(
        0.0, _wine_y0(schema, x=0.0, acetaldehyde=0.04), schema, params
    )
    assert schema.get(stranded, "acetaldehyde") == 0.0
    assert schema.get(stranded, "E") == 0.0


def test_reduction_runs_without_sugar(params):
    # NO flux term (unlike production): reduction must run during the rest, after sugar is
    # gone. With S=0 but live yeast and acetaldehyde present, it still clears.
    schema = wine_schema()
    d = AcetaldehydeReduction().derivatives(
        0.0, _wine_y0(schema, s=0.0, acetaldehyde=0.04), schema, params
    )
    assert schema.get(d, "acetaldehyde") < 0.0
    assert schema.get(d, "E") > 0.0


def test_reduction_zero_and_clamped_without_acetaldehyde(params):
    schema = wine_schema()
    empty = AcetaldehydeReduction().derivatives(
        0.0, _wine_y0(schema, acetaldehyde=0.0), schema, params
    )
    negative = AcetaldehydeReduction().derivatives(
        0.0, _wine_y0(schema, acetaldehyde=-1e-6), schema, params
    )
    assert schema.get(empty, "acetaldehyde") == 0.0
    # A solver undershoot below zero cannot manufacture ethanol from a negative pool.
    assert schema.get(negative, "acetaldehyde") == 0.0
    assert schema.get(negative, "E") == 0.0


# -- carbon: borrow then return is exactly neutral ----------------------------


def test_production_then_reduction_is_carbon_neutral(params):
    # The two Processes together move zero net carbon at equal magnitude: what production
    # borrows from E, reduction can return — the whole buffer is a carbon no-op on E+pool.
    schema = wine_schema()
    y = _wine_y0(schema, acetaldehyde=0.04)
    prod = AcetaldehydeProduction().derivatives(0.0, y, schema, params)
    red = AcetaldehydeReduction().derivatives(0.0, y, schema, params)
    for d in (prod, red):
        moved = schema.get(d, "acetaldehyde") * _ACETALDEHYDE_C + schema.get(d, "E") * _ETHANOL_C
        assert moved == pytest.approx(0.0, abs=1e-15)


# -- acceptance: the emergent early peak, carbon closure, E preservation ------


def _run(medium: str, celsius: float, days: float):
    """Compile+integrate a default ferment; return trajectory, compiled scenario, and the
    same run with the acetaldehyde Processes disabled (the buffer-off control)."""
    if medium == "wine":
        initial = {"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25}
    else:
        initial = {
            "glucose_gpl": 12.0,
            "maltose_gpl": 66.0,
            "maltotriose_gpl": 12.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 1.0,
        }
    scenario = Scenario(
        name=f"{medium}-acet",
        medium=medium,
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        duration_days=days,
    )
    compiled = compile_scenario(scenario, strict=True)
    dur = compiled.t_span_h[1]
    t_eval = np.linspace(0.0, dur, int(dur) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message

    off = compile_scenario(scenario, strict=True)
    for n in ("acetaldehyde_production", "acetaldehyde_reduction"):
        off.process_set.disable(n)
    traj_off = simulate(off.process_set, off.param_values, off.y0, off.t_span_h, t_eval=t_eval)
    assert traj_off.success, traj_off.message
    return traj, traj_off, compiled


@pytest.mark.parametrize(
    ("medium", "celsius", "days"),
    [("wine", 20.0, 21.0), ("beer", 18.0, 14.0)],
)
def test_acetaldehyde_produce_then_reabsorb_early_peak(medium, celsius, days):
    # The defining emergent behaviour (verified empirically before this assertion, D-27):
    # acetaldehyde rises to a peak DURING active fermentation and is then reduced to a low
    # residual — the produce-then-reabsorb "green apple" transient, not a monotone pool.
    traj, _, _ = _run(medium, celsius, days)
    t_days = traj.t / 24.0
    acet_mgl = np.asarray(traj.series("acetaldehyde")) * 1000.0  # g/L → mg/L
    peak_i = int(np.argmax(acet_mgl))
    peak = float(acet_mgl[peak_i])
    final = float(acet_mgl[-1])

    # Peak is a real, non-trivial excursion in the observed range (wine ~30-80, beer peaks
    # ~20-40 mg/L; threshold ~10-25 mg/L green apple) — the band carries margin.
    assert 10.0 < peak < 100.0, f"{medium} peak {peak:.1f} mg/L off-range"
    # The peak is EARLY — during active fermentation, well before the end.
    assert float(t_days[peak_i]) < 0.4 * days, (
        f"{medium} peak at day {t_days[peak_i]:.1f} not early"
    )
    # Reabsorbed: the final residual is a small fraction of the peak (yeast reduced it back).
    assert final < 0.1 * peak, f"{medium} final {final:.2f} not reabsorbed below peak {peak:.1f}"
    assert_nonnegative(traj, ("acetaldehyde",), atol=1e-12)


def test_carbon_closes_on_a_compiled_run():
    # With the acetaldehyde buffer wired into the default medium (D-27), a full compiled
    # ferment still conserves carbon to machine precision — the borrow/return is entirely on
    # the weighted E+acetaldehyde ledger (no sugar/CO2 term). Non-trivial because the pool
    # actually accumulates to ~40 mg/L mid-ferment.
    traj, _, compiled = _run("wine", 20.0, 21.0)
    f_c = compiled.param_values["biomass_C_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )


# == D-47: SO₂-bound acetaldehyde is protected from ADH (the free/bound RHS coupling) ==========
#
# The D-28 free/bound SO₂ split now feeds back into the reduction: alcohol dehydrogenase reduces
# only UNBOUND acetaldehyde (literature: "acetaldehyde bound to SO₂ could not be metabolized by
# yeast during fermentation; only free acetaldehyde could impact metabolism"). So dosed SO₂ locks
# acetaldehyde in. This retires the D-22/D-28 "SO₂ is readout-only" invariant for sulfited runs;
# an unsulfited run is byte-for-byte the D-27 core (the ``so2_total > 0`` guard is exact).


@pytest.fixture
def params_so2():
    # Reduction's SO₂ protection needs the sulfurous pKas + binding K + the pH-solver pKas, which
    # live in acidbase.yaml — layer it over the wine + acetaldehyde kinetics.
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "acetaldehyde.yaml",
        default_data_dir() / "acidbase.yaml",
    ).resolve()


def _sulfitable_state(
    schema: StateSchema, params, *, acetaldehyde_mgl: float, so2_mgl: float, x: float = 1.0
) -> FloatArray:
    """A mid-ferment wine state at pH ~3.4 with viable yeast, dosed acetaldehyde and (maybe) SO₂."""
    tartaric, malic = 6.0, 3.0
    totals = {"tartaric": tartaric / M_TARTARIC, "malic": malic / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(totals, 0.0, acidbase.build_pka_map(params), 3.4)
    y = schema.pack(
        {
            "X": x, "S": [150.0], "E": 40.0, "N": 0.1, "T": 293.15, "CO2": 0.0,
            "tartaric": tartaric, "malic": malic, "cation_charge": cation,
        }
    )  # fmt: skip
    y[schema.slice("acetaldehyde")] = mgl_to_gpl(acetaldehyde_mgl)
    if so2_mgl > 0.0:
        y[schema.slice("so2_total")] = mgl_to_gpl(so2_mgl)
    return y


def _reduction_rate(schema: StateSchema, params, y: FloatArray) -> float:
    """Mass rate of acetaldehyde reduction (g/L/h) at ``y`` — positive = clearing."""
    d = AcetaldehydeReduction().derivatives(0.0, y, schema, params)
    return -float(d[schema.slice("acetaldehyde")][0])


def test_unsulfited_reduction_is_byte_for_byte_the_closed_form(params_so2):
    # The guard is EXACT: with no SO₂ dosed, the reduction reads the *total* acetaldehyde — the
    # D-27 closed form — so the coupling is inert (no pH brentq, no protection). free == total.
    schema = wine_schema()
    y = _sulfitable_state(schema, params_so2, acetaldehyde_mgl=50.0, so2_mgl=0.0)
    acet = float(y[schema.slice("acetaldehyde")][0])
    f_t = arrhenius_factor(293.15, params_so2["E_a_acet_reduction"], params_so2["T_ref"])
    closed_form = params_so2["k_acet_reduction"] * 1.0 * f_t * acet
    assert _reduction_rate(schema, params_so2, y) == pytest.approx(closed_form)
    # free_acetaldehyde returns the whole pool when SO₂ is absent (the readout side of the guard)
    ph = acidbase.ph_of_state(y, schema, params_so2)
    assert acidbase.free_acetaldehyde(y, schema, params_so2, ph) == pytest.approx(acet)


def test_so2_throttles_the_reduction_to_the_free_share(params_so2):
    # Dosed SO₂ binds acetaldehyde and protects it: the reduction rate falls to (free/total) of
    # the unsulfited rate. Near-stoichiometric SO₂ throttles it hard; a large excess ~arrests it.
    schema = wine_schema()
    y_clean = _sulfitable_state(schema, params_so2, acetaldehyde_mgl=50.0, so2_mgl=0.0)
    y_comparable = _sulfitable_state(schema, params_so2, acetaldehyde_mgl=50.0, so2_mgl=60.0)
    y_excess = _sulfitable_state(schema, params_so2, acetaldehyde_mgl=50.0, so2_mgl=400.0)
    r_clean = _reduction_rate(schema, params_so2, y_clean)
    r_comparable = _reduction_rate(schema, params_so2, y_comparable)
    r_excess = _reduction_rate(schema, params_so2, y_excess)
    assert r_clean > 0.0
    assert r_comparable < 0.5 * r_clean  # comparable molar SO₂ ⇒ most acetaldehyde bound
    assert r_excess < 0.02 * r_clean  # SO₂ ≫ acetaldehyde ⇒ reduction ~fully arrested
    # the rate is exactly k·X·f(T)·free — pin it against the free-acetaldehyde readout
    ph = acidbase.ph_of_state(y_comparable, schema, params_so2)
    free = acidbase.free_acetaldehyde(y_comparable, schema, params_so2, ph)
    f_t = arrhenius_factor(293.15, params_so2["E_a_acet_reduction"], params_so2["T_ref"])
    assert r_comparable == pytest.approx(params_so2["k_acet_reduction"] * 1.0 * f_t * free)


def _so2_wine_run(so2_at_pitch: float = 0.0, so2_late: tuple[float, float] | None = None):
    """A default wine ferment, optionally dosed SO₂ at pitch and/or via a late ``add_so2``."""
    initial: dict[str, float] = {
        "brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5,
        "tartaric_gpl": 6.0, "malic_gpl": 3.0, "initial_ph": 3.4,
    }  # fmt: skip
    if so2_at_pitch > 0.0:
        initial["so2_total_mgl"] = so2_at_pitch
    interventions = []
    if so2_late is not None:
        day, mgl = so2_late
        interventions.append(Intervention(day=day, action="add_so2", params={"so2_mgl": mgl}))
    scenario = Scenario(
        name="wine-so2-acet", medium="wine", initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0, interventions=interventions,
    )  # fmt: skip
    c = compile_scenario(scenario, strict=True)
    traj = c.run()  # the intervention-aware runner (segment-and-restart at each add_so2)
    return traj, c


def test_post_af_so2_dose_strands_far_less_than_a_pitch_dose():
    # SCENARIO REALISM (the literature scoping): a 50 mg/L SO₂ dose at PITCH is sequestered by the
    # acetaldehyde peak and locks in a large residual; the SAME dose added POST-AF (day 16, after
    # the yeast has cleared acetaldehyde) finds little to bind, so it strands almost nothing and
    # its free SO₂ stays near the full dose. This is why SO₂ timing matters in the cellar.
    pitch, c_p = _so2_wine_run(so2_at_pitch=50.0)
    late, c_l = _so2_wine_run(so2_late=(16.0, 50.0))
    acet_pitch = gpl_to_mgl(pitch.series("acetaldehyde")[-1])
    acet_late = gpl_to_mgl(late.series("acetaldehyde")[-1])
    assert acet_pitch > 20.0  # pitch dose locks in a sensorily-relevant residual
    assert acet_late < 1.0  # post-AF dose strands ~nothing (acetaldehyde already gone)
    assert acet_pitch > 20.0 * acet_late
    # free SO₂ endpoint: depressed for the pitch dose, ~full for the late dose
    total = mgl_to_gpl(50.0)
    free_pitch = acidbase.free_so2(pitch.y[:, -1], pitch.schema, c_p.param_values)
    free_late = acidbase.free_so2(late.y[:, -1], late.schema, c_l.param_values)
    assert free_pitch < 0.4 * total
    assert free_late > 0.9 * total


def test_carbon_closes_on_a_sulfited_run():
    # Carbon is the SURVIVING invariant (D-47 retired only the trajectory-isolation one): the
    # reduction throttle only slows the acetaldehyde→E transfer, it neither creates nor routes
    # carbon, so a sulfited ferment that strands acetaldehyde still closes to machine precision.
    traj, c = _so2_wine_run(so2_at_pitch=50.0)
    assert gpl_to_mgl(traj.series("acetaldehyde")[-1]) > 20.0  # genuinely strands (non-trivial)
    f_c = c.param_values["biomass_C_fraction"]
    assert_conserved(traj, total_carbon(c.schema, biomass_carbon_fraction=f_c), label="carbon")


@pytest.mark.parametrize("method", ["RK45", "LSODA"])
def test_sulfited_reduction_agrees_across_solvers(method):
    # The reduction rate is now nonlinear in acetaldehyde/SO₂ (the bound_so2_molar quadratic root
    # + clamps) on an always-on RHS, and it runs under BDF's num_jac probe. Pin that the stranded
    # acetaldehyde endpoint agrees across BDF (default) and the RK45/LSODA references — no stiff
    # artefact, no probe-induced divergence (the D-40 pt2 idiom).
    bdf, _ = _so2_wine_run(so2_at_pitch=50.0)
    ref_c = compile_scenario(
        Scenario(
            name="wine-so2-ref", medium="wine",
            initial={
                "brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5, "tartaric_gpl": 6.0,
                "malic_gpl": 3.0, "initial_ph": 3.4, "so2_total_mgl": 50.0,
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)], duration_days=21.0,
        ),
        strict=True,
    )  # fmt: skip
    ref = simulate(
        ref_c.process_set, ref_c.param_values, ref_c.y0, ref_c.t_span_h, method=method
    )
    assert ref.success, ref.message
    assert gpl_to_mgl(bdf.series("acetaldehyde")[-1]) == pytest.approx(
        gpl_to_mgl(ref.series("acetaldehyde")[-1]), rel=5e-3, abs=0.5
    )


@pytest.mark.parametrize(
    ("medium", "celsius", "days"), [("wine", 20.0, 21.0), ("beer", 18.0, 14.0)]
)
def test_ethanol_endpoint_is_preserved_by_the_buffer(medium, celsius, days):
    # THE HEADLINE GUARANTEE (D-27): because acetaldehyde is a transient buffer that borrows
    # ethanol carbon and returns it, the final ethanol matches the buffer-off core to a
    # relative ~1e-8 (the pool fully reduces back), so the §2.2 ABV / realised-yield
    # benchmarks are preserved to far below any tolerance — the failure mode a draw-from-sugar
    # model does NOT avoid (it adds net-new ethanol scaling with pool turnover). Not bit-exact
    # only because ethanol feeds the ethanol-inactivation viability brake, so the transient E
    # dip perturbs the *path* at second order (see below); the endpoint reconverges once the
    # pool is reduced.
    traj, traj_off, _ = _run(medium, celsius, days)
    e_on = float(traj.series("E")[-1])
    e_off = float(traj_off.series("E")[-1])
    assert e_on == pytest.approx(e_off, rel=1e-6), f"{medium} E endpoint on={e_on} off={e_off}"

    # The buffer touches ONLY acetaldehyde and E at the derivative level (the closed-form unit
    # tests pin dS=dCO2=dN=0). So the other core outputs move only through the second-order
    # E→viability feedback above — a tiny path difference, not a direct draw — and the ferment
    # still completes to the same endpoints (compared here, not pointwise: a sub-hour timing
    # shift during the steep sugar crash makes a fixed-t path comparison spuriously large).
    def final_total(t: object, name: str) -> float:
        arr = np.asarray(t.series(name))  # type: ignore[attr-defined]
        arr = arr if arr.ndim == 1 else arr.sum(axis=0)
        return float(arr[-1])

    for var in ("S", "CO2", "N"):
        on = final_total(traj, var)
        off = final_total(traj_off, var)
        assert on == pytest.approx(off, rel=1e-3, abs=1e-6), (
            f"{medium} {var} endpoint on={on} off={off}"
        )


def test_ethanol_carries_a_transient_dip_while_buffered():
    # The honest flip side of the endpoint guarantee: DURING the ferment ethanol sits a hair
    # below the buffer-off core (carbon held as acetaldehyde), recovering as reduction
    # completes. Tiny (~0.04 g/L against ~118 g/L) and in the safe direction.
    traj, traj_off, _ = _run("wine", 20.0, 21.0)
    e_on = np.asarray(traj.series("E"))
    e_off = np.asarray(traj_off.series("E"))
    dip = e_off - e_on
    assert dip.max() > 0.0  # there IS a dip while acetaldehyde is elevated
    assert dip.max() < 0.5  # but it is negligible beside the ~118 g/L ethanol


# -- tier propagation ---------------------------------------------------------


def test_production_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AcetaldehydeProduction()])
    assert ps.tier_of("acetaldehyde") is Tier.SPECULATIVE
    assert ps.tier_of("acetaldehyde", store.tier_map()) is Tier.SPECULATIVE


def test_reduction_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AcetaldehydeReduction()])
    assert ps.tier_of("acetaldehyde") is Tier.SPECULATIVE
    assert ps.tier_of("acetaldehyde", store.tier_map()) is Tier.SPECULATIVE


def test_ethanol_tier_reflects_the_speculative_buffer(store):
    # The always-on speculative production is the FIRST byproduct Process to WRITE ethanol E
    # (esters/fusels/VDK touch S/CO2, MLF is disabled unpitched), so it drops the *structural*
    # tier_of("E") PLAUSIBLE→SPECULATIVE — the exact D-26 CO2 parallel, pinned here so it can
    # never regress silently (prime directive #1). Honest: E now genuinely holds a speculative
    # buffered slice. Crucially the param-aware tier users SEE is UNCHANGED — already
    # speculative, because the uptake Process itself reads speculative params (E_a_uptake,
    # realised-yield).
    schema = wine_schema()
    tm = store.tier_map()
    core = ProcessSet(schema, [SugarUptakeToEthanolCO2()])  # the validated ethanol producer
    with_buffer = ProcessSet(schema, [SugarUptakeToEthanolCO2(), AcetaldehydeProduction()])
    # Structural tier: PLAUSIBLE for the core alone, dropping to SPECULATIVE with the buffer.
    assert core.tier_of("E") is Tier.PLAUSIBLE
    assert with_buffer.tier_of("E") is Tier.SPECULATIVE
    # Param-aware tier (what users see): SPECULATIVE either way — no headline change.
    assert core.tier_of("E", tm) is Tier.SPECULATIVE
    assert with_buffer.tier_of("E", tm) is Tier.SPECULATIVE


# == D-48: SO₂-induced over-production — the transient-peak half of the elevation ==============
#
# D-47 protection shields *existing* acetaldehyde from ADH; D-48 adds the other half — trapping the
# terminal electron acceptor makes the yeast over-excrete acetaldehyde (the glyceropyruvic redox
# pull, Han 2020). It is scoped to the transient PEAK, not the end state: the finished-wine residual
# is capped by the SO₂-binding equilibrium (D-28), and D-47 protection ALONE already meets the field
# ~0.39 mg/mg total-acetaldehyde-vs-total-SO₂ slope, so an additive end-state bump would overshoot —
# D-48 leaves the stranded residual unchanged. Driver is TOTAL SO₂ (free SO₂ collapses to ~0 at the
# peak, so it is empirically inert there); the term is flux-gated (no runaway) and a carbon-exact
# borrow from E. The ``so2_total > 0`` guard is EXACT — an unsulfited run is byte-for-byte the core.


def _induced_run(so2_mgl: float, k_induced: float | None = None) -> np.ndarray:
    """A default wine ferment at a pitch SO₂ dose. ``k_induced=None`` uses the SHIPPED
    ``k_acet_so2_induced`` (so the tests track the YAML value); a float overrides it."""
    initial: dict[str, float] = {
        "brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5,
        "tartaric_gpl": 6.0, "malic_gpl": 3.0, "initial_ph": 3.4,
    }  # fmt: skip
    if so2_mgl > 0.0:
        initial["so2_total_mgl"] = so2_mgl
    c = compile_scenario(
        Scenario(
            name="wine-so2-induced", medium="wine", initial=initial,
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)], duration_days=21.0,
        ),
        strict=True,
    )  # fmt: skip
    pv = dict(c.param_values)  # the property returns a fresh resolve() — mutate a copy
    if k_induced is not None:
        pv["k_acet_so2_induced"] = k_induced
    traj = simulate(c.process_set, pv, c.y0, c.t_span_h)
    return np.asarray(traj.series("acetaldehyde"))


def test_induced_term_is_exact_when_undosed(params):
    # The so2_total>0 guard is EXACT: with no SO₂ the induced term is skipped entirely and the
    # production derivative is byte-for-byte the D-27 base borrow — no dependence on
    # k_acet_so2_induced whatsoever. Pin with == (not approx): the branch is simply not taken.
    schema = wine_schema()
    y = _wine_y0(schema, x=2.0, s=200.0)  # so2_total slot defaults to 0
    d = AcetaldehydeProduction().derivatives(0.0, y, schema, params)
    flux = 2.0 * (200.0 / (params["K_sugar_uptake"] + 200.0))
    assert schema.get(d, "acetaldehyde") == params["k_acetaldehyde"] * flux  # exact base only


def test_induced_production_closed_form(params):
    # With SO₂ dosed, production = base borrow + k_acet_so2_induced·flux·so2_total, all borrowed
    # from E mole-for-mole (carbon-exact). SO₂ is read-only; nothing but acetaldehyde and E moves.
    schema = wine_schema()
    x, s, so2 = 2.0, 200.0, mgl_to_gpl(100.0)
    y = _wine_y0(schema, x=x, s=s)
    y[schema.slice("so2_total")] = so2
    d = AcetaldehydeProduction().derivatives(0.0, y, schema, params)
    flux = x * (s / (params["K_sugar_uptake"] + s))
    rate = params["k_acetaldehyde"] * flux + params["k_acet_so2_induced"] * flux * so2
    assert schema.get(d, "acetaldehyde") == pytest.approx(rate)
    assert schema.get(d, "E") == pytest.approx(-rate * _ETH_PER_ACET)  # C2 borrow from ethanol
    # carbon-exact: the carbon leaving E equals the carbon entering acetaldehyde, per RHS
    assert schema.get(d, "E") * _ETHANOL_C == pytest.approx(
        -schema.get(d, "acetaldehyde") * _ACETALDEHYDE_C
    )
    for var in ("X", "S", "N", "CO2", "so2_total"):
        assert schema.get(d, var) == 0.0  # no sugar/CO2 draw; SO₂ is not consumed here


def test_induced_over_production_lifts_the_peak_not_the_end_state():
    # THE D-48 HEADLINE: the SO₂-induced bump raises the transient acetaldehyde PEAK but leaves the
    # finished-wine residual unchanged — the end state is capped by the SO₂-binding equilibrium
    # (an over-produced slice is reduced back; only the bound pool survives), so this is a genuine
    # peak-only effect, not a double-count of the end state D-47 already delivers.
    dose = 50.0
    off = _induced_run(dose, 0.0)  # D-47 protection only
    on = _induced_run(dose)  # + D-48 induced over-production (shipped k)
    peak_off, peak_on = gpl_to_mgl(off.max()), gpl_to_mgl(on.max())
    end_off, end_on = gpl_to_mgl(off[-1]), gpl_to_mgl(on[-1])
    assert peak_on > peak_off + 2.0  # a real, non-trivial peak lift (~+3.8 mg/L at 50 mg/L)
    assert end_on == pytest.approx(end_off, abs=0.5)  # end state / stranded residual unchanged


def test_induced_peak_lift_scales_with_dose():
    # Driven by TOTAL SO₂, so a larger dose lifts the peak more (near dose-proportional). Confirms
    # the driver is the dosed SO₂ state, not a saturating free-SO₂ readout (which is inert here).
    # Uses the SHIPPED k (None) so it tracks the YAML value, not a hardcoded coefficient.
    lift_50 = gpl_to_mgl(_induced_run(50.0).max()) - gpl_to_mgl(_induced_run(50.0, 0.0).max())
    lift_200 = gpl_to_mgl(_induced_run(200.0).max()) - gpl_to_mgl(_induced_run(200.0, 0.0).max())
    assert lift_50 > 1.0
    assert lift_200 > 2.5 * lift_50  # 4× the dose ⇒ a substantially larger peak lift
