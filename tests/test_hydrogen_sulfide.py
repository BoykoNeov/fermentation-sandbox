"""Tests for the hydrogen-sulfide (H₂S) Process (decision D-29).

H₂S is the "rotten egg" sulfidic off-aroma yeast releases when it reduces sulfate faster than
it can fix the sulfide onto nitrogen skeletons — so production is **de-repressed at low
nitrogen** (the inverse of the Ehrlich fusel gate). It is the accounting-easiest beat: H₂S is
carbon-free, on no conservation ledger, and its one Process (:class:`HydrogenSulfideProduction`)
touches only ``h2s`` — reading ``X``/``S``/``N`` without writing them.

The unit tests pin the closed form, the inverse-nitrogen gate direction/limits (the mechanism),
the temperature-flat/flux-linked shape and the guards. The acceptance section verifies:

* the **derivative-level isolability** is *exact* — disabling H₂S leaves every other column's
  RHS byte-for-byte (nothing reads ``h2s``); the integrated trajectory then differs only by a
  ~1e-7 adaptive-solver mesh artifact, cleaner than the acetaldehyde buffer's real second-order
  coupling (D-27);
* carbon still closes to machine precision on a compiled run (H₂S is carbon-free);
* the **emergent gate signal**: a low-YAN must produces markedly more H₂S *early* (its gate
  opens sooner) — even though it has *less* biomass — while the cumulative endpoint lever is
  **muted**, because the nitrogen model strips ``N`` to ~0 early regardless of dose (the
  documented upstream gap; D-29).
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics import HydrogenSulfideProduction
from fermentation.core.media import get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon


@pytest.fixture
def store():
    # Wine kinetics PLUS the shared H₂S constants (k_h2s, K_h2s_n).
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "hydrogen_sulfide.yaml",
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
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": t, "CO2": 0.0})


def _gate(params: Mapping[str, float], n: float) -> float:
    k = float(params["K_h2s_n"])
    return k / (k + n)


# -- metadata -----------------------------------------------------------------


def test_metadata():
    p = HydrogenSulfideProduction()
    assert p.name == "hydrogen_sulfide_production"
    # Speculative: rate magnitude is an order-of-magnitude estimate.
    assert p.tier is Tier.SPECULATIVE
    # Touches ONLY its own carbon-free pool — never S/CO2/E/N/X (reads them, writes none).
    assert set(p.touches) == {"h2s"}
    assert set(p.reads) == {"k_h2s", "K_sugar_uptake", "K_h2s_n"}


def test_h2s_is_carbon_free():
    # Registered with 0 carbon (like SO₂), so it is inert in every carbon sum — the reason
    # this beat needs no conservation code and closure is trivial.
    assert carbon_mass_fraction("hydrogen_sulfide") == 0.0


# -- closed form, gate, guards ------------------------------------------------


def test_matches_closed_form(params):
    schema = wine_schema()
    x, s, n = 2.0, 200.0, 0.1
    y = _wine_y0(schema, x=x, s=s, n=n)
    d = HydrogenSulfideProduction().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    rate = params["k_h2s"] * flux * _gate(params, n)
    assert schema.get(d, "h2s") == pytest.approx(rate)
    assert schema.get(d, "h2s") > 0.0
    # Carbon-free produced pool: nothing else moves — no sugar draw, no CO2, no ethanol/N/X.
    for var in ("X", "S", "E", "N", "CO2"):
        assert schema.get(d, var) == 0.0


def test_inverse_nitrogen_gate_direction(params):
    # THE mechanism, isolated at the derivative level (flux fixed, only N varies): H₂S
    # production is HIGHER at low nitrogen and LOWER at high nitrogen — de-repression, the
    # opposite of the Ehrlich fusel gate. This is the clean gate test (the integrated-run
    # cross-must version below is confounded by biomass, so this pins the direction).
    schema = wine_schema()
    low_n = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, n=0.01), schema, params)
    high_n = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, n=0.3), schema, params)
    assert schema.get(low_n, "h2s") > schema.get(high_n, "h2s") > 0.0


def test_gate_limits(params):
    # N → 0: gate → 1, rate is the full flux-linked production. N large: gate → 0, rate → 0.
    schema = wine_schema()
    x, s = 2.0, 200.0
    flux = x * (s / (params["K_sugar_uptake"] + s))
    depleted = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, n=0.0), schema, params)
    assert schema.get(depleted, "h2s") == pytest.approx(params["k_h2s"] * flux)  # gate == 1
    flooded = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, n=1e6), schema, params)
    assert schema.get(flooded, "h2s") == pytest.approx(0.0, abs=1e-12)  # gate → 0


def test_temperature_flat(params):
    # Documented v1 simplification (D-29): production carries NO Arrhenius factor. Identical
    # cold vs warm at the same flux/nitrogen state.
    schema = wine_schema()
    cold = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(cold, "h2s") == pytest.approx(schema.get(warm, "h2s"))
    assert schema.get(cold, "h2s") > 0.0


def test_scales_with_fermentative_flux(params):
    # Coupled to the biomass-catalysed sugar flux (linear in X): twice the biomass ⇒ twice
    # the release.
    schema = wine_schema()
    r1 = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "h2s") == pytest.approx(2.0 * schema.get(r1, "h2s"))


def test_zero_without_biomass_or_sugar(params):
    # Flux-linked ⇒ no biomass or no sugar ⇒ no production (stops at dryness).
    schema = wine_schema()
    no_x = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)
    no_s = HydrogenSulfideProduction().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)
    assert schema.get(no_x, "h2s") == 0.0
    assert schema.get(no_s, "h2s") == 0.0


# -- tier propagation ---------------------------------------------------------


def test_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [HydrogenSulfideProduction()])
    assert ps.tier_of("h2s") is Tier.SPECULATIVE
    assert ps.tier_of("h2s", store.tier_map()) is Tier.SPECULATIVE


# -- acceptance: isolability, carbon closure, emergent gate signal ------------


def _run(yan_mgl: float, *, celsius: float = 20.0, days: float = 21.0):
    """Compile+integrate a default wine ferment at the given YAN on a daily grid; return the
    trajectory and compiled scenario."""
    scenario = Scenario(
        name=f"wine-h2s-{yan_mgl:.0f}",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25},
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
    return traj, compiled


def test_isolable_at_derivative_level():
    # THE strong isolability claim (D-29): H₂S writes only its own pool, so removing it cannot
    # change any other column's RHS. Verified EXACTLY (0.0) across several states — stronger
    # than the acetaldehyde buffer, whose E write forces a second-order coupling (D-27). Uses
    # the FULL default wine ProcessSet (all shared kinetics), so it needs the merged params.
    base = default_data_dir()
    full = load_parameters(
        base / "wine_generic.yaml",
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "hydrogen_sulfide.yaml",
    ).resolve()
    m = get_medium("wine")
    schema = m.schema
    on = m.build_process_set()
    off = m.build_process_set()
    off.disable("hydrogen_sulfide_production")
    h2s_i = schema.slice("h2s").start
    for x, s, n in [(2.0, 200.0, 0.1), (3.0, 50.0, 0.0), (1.0, 150.0, 0.05)]:
        y = schema.pack({"X": x, "S": [s], "E": 60.0, "N": n, "T": 293.15, "CO2": 10.0})
        diff = on.total_derivatives(0.0, y, full) - off.total_derivatives(0.0, y, full)
        other = np.delete(diff, h2s_i)
        assert np.max(np.abs(other)) == 0.0  # byte-for-byte on every non-h2s column
        assert diff[h2s_i] > 0.0  # ...and H₂S itself is the only thing that moved


def test_trajectory_isolability_is_solver_mesh_only():
    # At the integrated level the non-h2s columns are NOT bit-identical, but only by a ~1e-7
    # adaptive-solver mesh artifact (adding the h2s equation shifts step selection) — there is
    # no physical pathway, since nothing reads h2s. Pinned loosely to document it as numerical.
    traj_on, compiled = _run(80.0)
    off = compile_scenario(
        Scenario(
            name="wine-h2s-off",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=21.0,
        ),
        strict=True,
    )
    off.process_set.disable("hydrogen_sulfide_production")
    t_eval = traj_on.t
    traj_off = simulate(off.process_set, off.param_values, off.y0, off.t_span_h, t_eval=t_eval)
    assert traj_off.success
    for name in compiled.schema.names:
        if name == "h2s":
            continue
        a = np.asarray(traj_on.series(name))
        b = np.asarray(traj_off.series(name))
        scale = max(float(np.max(np.abs(b))), 1e-9)
        assert float(np.max(np.abs(a - b))) / scale < 1e-5, f"{name} drifts more than mesh noise"


def test_carbon_closes_on_a_compiled_run():
    # H₂S is carbon-free, so a full compiled ferment with it wired in (D-29) still conserves
    # carbon to machine precision — H₂S is simply absent from total_carbon. Non-trivial because
    # the pool actually accumulates (~0.5 mg/L) yet contributes nothing to the ledger.
    traj, compiled = _run(80.0)
    f_c = compiled.param_values["biomass_C_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert float(traj.series("h2s")[-1]) > 0.0  # it really did produce H₂S


def test_h2s_is_produced_only_and_plateaus_at_dryness():
    # Produced-only pool: monotonically non-decreasing (no reduction/stripping in v1), rising
    # during active fermentation and flat once sugar is gone. Honest caveat: this is cumulative
    # PRODUCED, which overstates residual (real H₂S is largely CO2-stripped).
    traj, _ = _run(80.0)
    h2s = np.asarray(traj.series("h2s"))
    assert_nonnegative(traj, ("h2s",), atol=1e-12)
    # non-decreasing (allow tiny solver noise)
    assert np.min(np.diff(h2s)) > -1e-9
    # a real, non-trivial excursion
    assert h2s[-1] > 1e-4
    # plateaus: the last-quarter growth is a tiny fraction of the total (production has stopped)
    q = len(h2s) // 4
    late_growth = h2s[-1] - h2s[-q]
    assert late_growth < 0.02 * h2s[-1]


def test_low_yan_produces_more_h2s_early():
    # THE emergent gate signal, integrated (D-29): a low-YAN must's inverse-N gate opens sooner,
    # so it produces markedly more H₂S EARLY — even though it grows LESS biomass (so this is the
    # gate winning over flux, not a flux artifact). Sampled at day 1, while N still differs
    # between the musts (both deplete to ~0 by ~day 1.3).
    lo, _ = _run(80.0)
    hi, _ = _run(300.0)
    day = 1.0
    i_lo = int(np.argmin(np.abs(lo.t / 24.0 - day)))
    i_hi = int(np.argmin(np.abs(hi.t / 24.0 - day)))
    # low-YAN grows less biomass but makes more H₂S — the gate dominates the higher flux.
    assert float(hi.series("X")[i_hi]) > float(lo.series("X")[i_lo])
    assert float(lo.series("h2s")[i_lo]) > 1.4 * float(hi.series("h2s")[i_hi])


def test_cross_must_endpoint_lever_is_muted():
    # THE honest caveat (D-29): the cumulative-endpoint YAN lever is MUTED — a low- and a
    # high-YAN must end within a few percent — because the nitrogen model strips N to ~0 early
    # regardless of dose (no residual-N floor), so the inverse gate is ~1 for the rest of the
    # ferment for BOTH. Direction is preserved (low ≥ high) but the gap is small; the lever
    # becomes real only once the deferred residual-N floor lands. NOT a hollow assertion — this
    # pins the gap as small on purpose, documenting the gap rather than papering over it.
    lo, _ = _run(80.0)
    hi, _ = _run(300.0)
    final_lo = float(lo.series("h2s")[-1])
    final_hi = float(hi.series("h2s")[-1])
    assert final_lo >= final_hi  # direction preserved (low YAN ≥ high YAN)
    assert (final_lo - final_hi) / final_lo < 0.15  # ...but muted: within ~15%
