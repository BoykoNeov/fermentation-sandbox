"""Tests for the vicinal-diketone (VDK / diacetyl) pathway Processes (decision D-26).

The diacetyl beat models the *real* three-step pathway — α-acetolactate excretion →
spontaneous decarboxylation → yeast reduction — so the produce-then-reabsorb "diacetyl
rest" emerges rather than being scripted:

* **Excretion** (:class:`AcetolactateExcretion`): fills the α-acetolactate reservoir from
  the fermentative flux, routing its C5 carbon out of ``S`` (option a1, D-19). Temperature-
  flat (a documented v1 simplification — the reservoir size is a weak lever).
* **Decarboxylation** (:class:`AcetolactateDecarboxylation`): spontaneous, non-yeast-gated,
  strongly temperature-dependent C5 → C4 + CO2 (carbon-closing like malic → lactic + CO2).
  The rate-limiting, temperature-critical step.
* **Reduction** (:class:`DiacetylReduction`): enzymatic C4 → C4 to 2,3-butanediol, gated on
  *viable* yeast with no flux term — so it clears diacetyl during the rest but stops when
  the yeast is crashed.

The unit tests pin each Process's closed form, carbon accounting and guards; the acceptance
section verifies the *emergent* diacetyl rest — warmer ⇒ cleaner, peak-then-fall when warm,
and cold stranding with an unconverted reservoir — plus carbon closure on a compiled run.
"""

import numpy as np
import pytest

from fermentation.core.chemistry import (
    M_ACETOLACTATE,
    M_BUTANEDIOL,
    M_CO2,
    M_DIACETYL,
    carbon_mass_fraction,
)
from fermentation.core.kinetics import (
    AcetolactateDecarboxylation,
    AcetolactateExcretion,
    DiacetylReduction,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
    arrhenius_factor,
)
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon

#: Carbon fractions the pools book against (mirror the Process/chemistry constants).
_ACETOLACTATE_C = carbon_mass_fraction("alpha_acetolactate")
_DIACETYL_C = carbon_mass_fraction("diacetyl")
_CO2_C = carbon_mass_fraction("CO2")
_GLUCOSE_C = carbon_mass_fraction("glucose")


@pytest.fixture
def store():
    # Wine kinetics PLUS the shared VDK pathway constants (k_acetolactate, k_decarb, …).
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "vicinal_diketones.yaml",
    )


@pytest.fixture
def params(store):
    return store.resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x: float = 2.0,
    s: float = 200.0,
    e: float = 0.0,
    n: float = 0.1,
    t: float = 293.15,
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": t, "CO2": 0.0})


# -- metadata -----------------------------------------------------------------


def test_excretion_metadata():
    p = AcetolactateExcretion()
    assert p.name == "acetolactate_excretion"
    # Speculative: rate magnitude is an order-of-magnitude estimate.
    assert p.tier is Tier.SPECULATIVE
    # Touches its own reservoir pool AND S — the C5 carbon is routed from sugar (a1, D-19);
    # never E/CO2 (the reservoir is the precursor, not a fermentation product yet).
    assert set(p.touches) == {"acetolactate", "S"}
    assert set(p.reads) == {"k_acetolactate", "K_sugar_uptake"}


# -- closed form & guards -----------------------------------------------------


def test_excretion_derivative_matches_closed_form(params):
    schema = wine_schema()
    x, s = 2.0, 200.0
    y = _wine_y0(schema, x=x, s=s)
    d = AcetolactateExcretion().derivatives(0.0, y, schema, params)

    flux = x * (s / (params["K_sugar_uptake"] + s))
    rate = params["k_acetolactate"] * flux
    assert schema.get(d, "acetolactate") == pytest.approx(rate)
    # The C5 carbon is routed from sugar (a1): one slot, so dS removes exactly the
    # acetolactate carbon converted back to grams of glucose. Carbon balances per-RHS.
    assert schema.get(d, "S") == pytest.approx(-rate * _ACETOLACTATE_C / _GLUCOSE_C)
    assert schema.get(d, "S") * _GLUCOSE_C == pytest.approx(
        -schema.get(d, "acetolactate") * _ACETOLACTATE_C
    )
    # Nothing else moves — the reservoir is a precursor; diacetyl/CO2 come later.
    for var in ("X", "E", "N", "CO2", "diacetyl", "butanediol"):
        assert schema.get(d, var) == 0.0


def test_excretion_is_temperature_flat(params):
    # Documented v1 simplification (D-26): excretion carries NO Arrhenius factor — the
    # reservoir size is a weak lever on the rest; temperature-criticality lives in the
    # decarboxylation. So the excretion rate is identical cold vs warm at equal flux.
    schema = wine_schema()
    cold = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, t=283.15), schema, params)
    warm = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, t=303.15), schema, params)
    assert schema.get(cold, "acetolactate") == pytest.approx(schema.get(warm, "acetolactate"))
    assert schema.get(cold, "acetolactate") > 0.0


def test_excretion_scales_with_fermentative_flux(params):
    # Coupled to the biomass-catalysed sugar flux (linear in X): twice the biomass ⇒
    # twice the excretion.
    schema = wine_schema()
    r1 = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=1.0), schema, params)
    r2 = AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params)
    assert schema.get(r2, "acetolactate") == pytest.approx(2.0 * schema.get(r1, "acetolactate"))


def test_excretion_zero_without_biomass_or_sugar(params):
    schema = wine_schema()
    assert AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=0.0), schema, params)[
        schema.slice("acetolactate")
    ] == pytest.approx(0.0)
    assert AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, s=0.0), schema, params)[
        schema.slice("acetolactate")
    ] == pytest.approx(0.0)


def test_excretion_negative_excursion_does_not_produce(params):
    # A solver undershoot (S<0 or X<0) must not flip the clamp and create acetolactate.
    schema = wine_schema()
    assert np.array_equal(
        AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, s=-1e-6), schema, params),
        schema.zeros(),
    )
    assert np.array_equal(
        AcetolactateExcretion().derivatives(0.0, _wine_y0(schema, x=-1e-6), schema, params),
        schema.zeros(),
    )


# -- integration-level properties ---------------------------------------------


def test_excretion_runs_strict_accumulates_and_closes_carbon(params, store):
    # Under the strict touches contract: excretion writes only its reservoir pool (+ the
    # sugar draw), the reservoir accumulates non-negatively and stays TRACE, and — because
    # the C5 carbon is drawn from sugar AND acetolactate is weighted in total_carbon — a
    # full ferment conserves carbon to machine precision (the reservoir is real
    # carbon-accounted state, not an unbooked leak).
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2(), AcetolactateExcretion()],
        strict=True,
    )
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    assert_nonnegative(traj, ("acetolactate", "S"), atol=1e-12)
    # The reservoir actually fills (mechanism live) but stays trace (~mg/L), orders of
    # magnitude below the g/L ethanol flux.
    assert 0.0 < float(traj.series("acetolactate")[-1]) < 0.1
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


def test_excretion_perturbs_only_sugar_per_rhs(params):
    # Isolability (prime directive #3) under a1: enabling excretion leaves dX/dN/dE/dCO2
    # byte-for-byte (it never touches them); the only change is on dS, which loses EXACTLY
    # the carbon deposited in the reservoir. That per-RHS balance is why total_carbon closes.
    schema = wine_schema()
    core = ProcessSet(schema, [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2()])
    with_vdk = ProcessSet(
        schema, [GrowthNitrogenLimited(), SugarUptakeToEthanolCO2(), AcetolactateExcretion()]
    )
    for s, e, n in ((245.0, 0.0, 0.08), (120.0, 60.0, 0.01), (5.0, 120.0, 0.0)):
        y = _wine_y0(schema, x=1.5, s=s, e=e, n=n)
        d_core = core.total_derivatives(0.0, y, params)
        d_vdk = with_vdk.total_derivatives(0.0, y, params)
        for var in ("X", "E", "N", "CO2"):
            assert d_core[schema.slice(var)] == pytest.approx(d_vdk[schema.slice(var)], abs=0.0)
        delta_s = float(d_vdk[schema.slice("S")][0] - d_core[schema.slice("S")][0])
        ala_rate = float(d_vdk[schema.slice("acetolactate")][0])
        assert delta_s <= 0.0  # sugar drawn down, never created
        carbon_residual = delta_s * _GLUCOSE_C + ala_rate * _ACETOLACTATE_C
        assert carbon_residual == pytest.approx(0.0, abs=1e-12)


# -- decarboxylation: α-acetolactate -> diacetyl + CO2 (spontaneous, decision D-26) --


def _wine_y0_with_acetolactate(schema: StateSchema, *, acetolactate: float, **kw) -> FloatArray:
    """A wine state with the α-acetolactate reservoir pre-loaded (the decarb needs a
    reservoir to convert; the produced-only pool is 0 at pitch so it must be set here)."""
    y = _wine_y0(schema, **kw)
    y[schema.slice("acetolactate")] = acetolactate
    return y


def test_decarb_metadata():
    p = AcetolactateDecarboxylation()
    assert p.name == "acetolactate_decarboxylation"
    assert p.tier is Tier.SPECULATIVE
    # Carbon-closing decarboxylation: reservoir -> diacetyl + CO2, no sugar draw.
    assert set(p.touches) == {"acetolactate", "diacetyl", "CO2"}
    assert set(p.reads) == {"k_decarb", "E_a_decarb", "T_ref"}


def test_decarb_derivative_matches_closed_form_and_closes_carbon(params):
    schema = wine_schema()
    ala, t = 0.02, 298.15  # off T_ref so the Arrhenius factor bites
    y = _wine_y0_with_acetolactate(schema, acetolactate=ala, t=t)
    d = AcetolactateDecarboxylation().derivatives(0.0, y, schema, params)

    f_t = arrhenius_factor(t, params["E_a_decarb"], params["T_ref"])
    r = params["k_decarb"] * f_t * ala / M_ACETOLACTATE  # molar turnover
    assert schema.get(d, "acetolactate") == pytest.approx(-r * M_ACETOLACTATE)
    assert schema.get(d, "diacetyl") == pytest.approx(r * M_DIACETYL)
    assert schema.get(d, "CO2") == pytest.approx(r * M_CO2)
    # Carbon closes mole-for-mole (C5 -> C4 + C1), like malic -> lactic + CO2 (D-23).
    carbon_residual = (
        schema.get(d, "acetolactate") * _ACETOLACTATE_C
        + schema.get(d, "diacetyl") * _DIACETYL_C
        + schema.get(d, "CO2") * _CO2_C
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-15)
    # No sugar drawn (unlike excretion), no ethanol/biomass touched.
    for var in ("X", "S", "E", "N", "butanediol"):
        assert schema.get(d, var) == 0.0


def test_decarb_factor_is_one_at_reference_temperature(params):
    schema = wine_schema()
    ala = 0.02
    y = _wine_y0_with_acetolactate(schema, acetolactate=ala, t=params["T_ref"])
    d = AcetolactateDecarboxylation().derivatives(0.0, y, schema, params)
    r = params["k_decarb"] * ala / M_ACETOLACTATE
    assert schema.get(d, "diacetyl") == pytest.approx(r * M_DIACETYL)


def test_decarb_rises_with_temperature(params):
    # The load-bearing property (D-26): the spontaneous conversion accelerates with
    # temperature — this is what makes a warm rest clear diacetyl faster.
    schema = wine_schema()
    cold = AcetolactateDecarboxylation().derivatives(
        0.0, _wine_y0_with_acetolactate(schema, acetolactate=0.02, t=283.15), schema, params
    )
    warm = AcetolactateDecarboxylation().derivatives(
        0.0, _wine_y0_with_acetolactate(schema, acetolactate=0.02, t=303.15), schema, params
    )
    assert schema.get(warm, "diacetyl") > schema.get(cold, "diacetyl") > 0.0


def test_decarb_is_not_yeast_gated(params):
    # THE mechanistic point (D-26): the reaction is non-enzymatic, so it proceeds with NO
    # viable yeast (X=0) — the reason the reservoir keeps making diacetyl after a crash.
    schema = wine_schema()
    d = AcetolactateDecarboxylation().derivatives(
        0.0, _wine_y0_with_acetolactate(schema, acetolactate=0.02, x=0.0), schema, params
    )
    assert schema.get(d, "diacetyl") > 0.0


def test_decarb_first_order_in_acetolactate(params):
    # First-order: twice the reservoir ⇒ twice the conversion rate.
    schema = wine_schema()
    r1 = AcetolactateDecarboxylation().derivatives(
        0.0, _wine_y0_with_acetolactate(schema, acetolactate=0.01), schema, params
    )
    r2 = AcetolactateDecarboxylation().derivatives(
        0.0, _wine_y0_with_acetolactate(schema, acetolactate=0.02), schema, params
    )
    assert schema.get(r2, "diacetyl") == pytest.approx(2.0 * schema.get(r1, "diacetyl"))


def test_decarb_zero_and_clamped_without_reservoir(params):
    schema = wine_schema()
    # No reservoir ⇒ nothing to convert.
    assert np.array_equal(
        AcetolactateDecarboxylation().derivatives(0.0, _wine_y0(schema), schema, params),
        schema.zeros(),
    )
    # Solver undershoot (acetolactate < 0) must not manufacture diacetyl.
    assert np.array_equal(
        AcetolactateDecarboxylation().derivatives(
            0.0, _wine_y0_with_acetolactate(schema, acetolactate=-1e-9), schema, params
        ),
        schema.zeros(),
    )


def test_excretion_plus_decarb_accumulates_diacetyl_and_closes_carbon(params, store):
    # Two steps wired together: the reservoir fills (excretion) and converts to diacetyl +
    # CO2 (decarb). Diacetyl accumulates (nothing reduces it yet — reduction is commit 3),
    # stays trace, and carbon closes to machine precision across sugar -> reservoir ->
    # diacetyl + CO2 (every step is on the weighted ledger).
    schema = wine_schema()
    ps = ProcessSet(
        schema,
        [
            GrowthNitrogenLimited(),
            SugarUptakeToEthanolCO2(),
            AcetolactateExcretion(),
            AcetolactateDecarboxylation(),
        ],
        strict=True,
    )
    y0 = schema.pack({"X": 0.25, "S": [245.0], "E": 0.0, "N": 0.08, "T": 293.15, "CO2": 0.0})
    traj = simulate(ps, params=params, y0=y0, t_span=(0.0, 400.0))
    assert traj.success
    assert_nonnegative(traj, ("acetolactate", "diacetyl", "S"), atol=1e-12)
    assert 0.0 < float(traj.series("diacetyl")[-1]) < 0.1  # trace, buttery off-note
    f_c = store.value("biomass_C_fraction")
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")


# -- reduction: diacetyl -> 2,3-butanediol (viable-yeast-gated, decision D-26) --


def _wine_y0_with_diacetyl(schema: StateSchema, *, diacetyl: float, **kw) -> FloatArray:
    """A wine state with the ``diacetyl`` pool pre-loaded (the reduction needs diacetyl to
    reduce; the produced-only pool is 0 at pitch so it must be set here)."""
    y = _wine_y0(schema, **kw)
    y[schema.slice("diacetyl")] = diacetyl
    return y


def test_reduction_metadata():
    p = DiacetylReduction()
    assert p.name == "diacetyl_reduction"
    assert p.tier is Tier.SPECULATIVE
    # Carbon-neutral C4 -> C4 transfer to the terminal diol pool; no sugar/CO2 touched.
    assert set(p.touches) == {"diacetyl", "butanediol"}
    assert set(p.reads) == {"k_reduction", "E_a_reduction", "T_ref"}


def test_reduction_derivative_matches_closed_form_and_conserves_carbon(params):
    schema = wine_schema()
    dia, x, t = 0.001, 2.0, 298.15  # off T_ref so the Arrhenius factor bites
    y = _wine_y0_with_diacetyl(schema, diacetyl=dia, x=x, t=t)
    d = DiacetylReduction().derivatives(0.0, y, schema, params)

    f_t = arrhenius_factor(t, params["E_a_reduction"], params["T_ref"])
    loss = params["k_reduction"] * x * f_t * dia
    assert schema.get(d, "diacetyl") == pytest.approx(-loss)
    assert schema.get(d, "butanediol") == pytest.approx(loss * M_BUTANEDIOL / M_DIACETYL)
    # Mole-for-mole C4 -> C4: carbon leaving diacetyl equals carbon entering butanediol.
    carbon_residual = schema.get(d, "diacetyl") * _DIACETYL_C + schema.get(d, "butanediol") * (
        carbon_mass_fraction("butanediol")
    )
    assert carbon_residual == pytest.approx(0.0, abs=1e-18)
    for var in ("X", "S", "E", "N", "CO2", "acetolactate"):
        assert schema.get(d, var) == 0.0


def test_reduction_requires_viable_yeast(params):
    # THE stranding mechanism (D-26): reduction is enzymatic, so with NO viable yeast
    # (X=0) it stops dead — the reservoir keeps making diacetyl (decarb is not yeast-gated)
    # but nothing reduces it, so diacetyl is stranded. This is why packaging/crashing too
    # early leaves diacetyl high.
    schema = wine_schema()
    d = DiacetylReduction().derivatives(
        0.0, _wine_y0_with_diacetyl(schema, diacetyl=0.001, x=0.0), schema, params
    )
    assert np.array_equal(d, schema.zeros())


def test_reduction_has_no_flux_term(params):
    # Must run during the REST (after sugar is gone, flux ≈ 0). Reduction depends only on
    # viable X and diacetyl, NOT on the fermentative flux — so a sugar-free state (S=0)
    # with live yeast still reduces diacetyl (contrast excretion, which needs the flux).
    schema = wine_schema()
    d = DiacetylReduction().derivatives(
        0.0, _wine_y0_with_diacetyl(schema, diacetyl=0.001, x=2.0, s=0.0), schema, params
    )
    assert schema.get(d, "diacetyl") < 0.0
    assert schema.get(d, "butanediol") > 0.0


def test_reduction_scales_with_viable_biomass_and_diacetyl(params):
    schema = wine_schema()
    base = DiacetylReduction().derivatives(
        0.0, _wine_y0_with_diacetyl(schema, diacetyl=0.001, x=1.0), schema, params
    )
    more_x = DiacetylReduction().derivatives(
        0.0, _wine_y0_with_diacetyl(schema, diacetyl=0.001, x=2.0), schema, params
    )
    more_dia = DiacetylReduction().derivatives(
        0.0, _wine_y0_with_diacetyl(schema, diacetyl=0.002, x=1.0), schema, params
    )
    assert schema.get(more_x, "diacetyl") == pytest.approx(2.0 * schema.get(base, "diacetyl"))
    assert schema.get(more_dia, "diacetyl") == pytest.approx(2.0 * schema.get(base, "diacetyl"))


def test_reduction_zero_and_clamped_without_diacetyl(params):
    schema = wine_schema()
    assert np.array_equal(
        DiacetylReduction().derivatives(0.0, _wine_y0(schema, x=2.0), schema, params),
        schema.zeros(),
    )
    assert np.array_equal(
        DiacetylReduction().derivatives(
            0.0, _wine_y0_with_diacetyl(schema, diacetyl=-1e-12, x=2.0), schema, params
        ),
        schema.zeros(),
    )


# -- the emergent diacetyl rest (the D-26 acceptance gate) ---------------------


def _run_diacetyl(medium: str, celsius: float, days: float):
    """Compile + run a medium isothermally with the full VDK pathway wired in; return a
    dict of the diacetyl / α-acetolactate trajectories (mg/L) and the compiled run so the
    carbon check can read the resolved biomass fraction. The three-step pathway makes the
    diacetyl rest *emerge* — nothing here scripts a peak or a rest."""
    if medium == "wine":
        init = {"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25}
    else:
        init = {
            "glucose_gpl": 13.2,
            "maltose_gpl": 54.6,
            "maltotriose_gpl": 20.2,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.6,
        }
    sc = Scenario(
        name=f"{medium}-{celsius}C",
        medium=medium,
        initial=init,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        duration_days=days,
    )
    compiled = compile_scenario(sc, strict=True)
    t_eval = np.linspace(0.0, days * 24.0, 400)
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h,
                    t_eval=t_eval)  # fmt: skip
    assert traj.success, traj.message
    dia = traj.series("diacetyl") * 1000.0  # mg/L
    return {
        "diacetyl_mgl": dia,
        "diacetyl_peak": float(dia.max()),
        "diacetyl_final": float(dia[-1]),
        "acetolactate_final_mgl": float(traj.series("acetolactate")[-1] * 1000.0),
        "butanediol_final_mgl": float(traj.series("butanediol")[-1] * 1000.0),
        "traj": traj,
        "compiled": compiled,
    }


@pytest.mark.parametrize(
    ("medium", "temps", "days"),
    [
        ("beer", (10.0, 18.0, 25.0), 30.0),
        ("wine", (14.0, 20.0, 28.0), 45.0),
    ],
)
def test_warmer_ferment_is_cleaner_the_diacetyl_rest(medium, temps, days):
    # THE D-26 headline (the money shot): final diacetyl falls monotonically as fermentation
    # temperature rises — "warm it up to clear the diacetyl faster." This EMERGES from the
    # three-step pathway (temperature-critical spontaneous decarb emptying the reservoir
    # faster + faster enzymatic reduction while yeast is viable), not from any scripted rest.
    finals = [_run_diacetyl(medium, c, days)["diacetyl_final"] for c in temps]
    assert finals[0] > finals[1] > finals[2] > 0.0, (
        f"{medium} final diacetyl should fall with T (warm rest cleaner): "
        f"{dict(zip(temps, finals, strict=True))} mg/L"
    )
    # The cold run leaves a perceptible buttery note (above the ~0.1 mg/L lager threshold).
    assert finals[0] > 0.1
    # D-57 gave EthanolInactivation its (previously missing) real Coleman quadratic
    # temperature scaling, so warm ferments now correctly lose viable/reductase-capable
    # biomass faster than before (when death was silently frozen at the 20 C rate at
    # EVERY temperature). Wine's 28 C run no longer clears diacetyl below the ~0.1 mg/L
    # perceptibility threshold on isolated yeast reductase alone (0.162 mg/L measured,
    # was ~0.03 pre-D-57) -- a real, sourced physics change, not a loosened pass. The
    # practical "warm it up" claim still holds as a large, measured reduction (~3x cold
    # vs warm here); asserting sub-perceptibility would need additional reductase
    # capacity this isolated-yeast scenario doesn't model (e.g. MLF bacteria). Beer's
    # much lower ethanol tolerance keeps k_prime_d inert (see beer_generic.yaml), so its
    # warm run still clears comfortably below threshold — no relative-margin needed there.
    if medium == "wine":
        assert finals[0] / finals[-1] > 2.5, (
            f"warm rest should still cut diacetyl by a large factor: {finals}"
        )
    else:
        assert finals[-1] < 0.1


def test_warm_rest_shows_peak_then_fall_and_clears():
    # Produce-then-reabsorb: at a warm temperature diacetyl RISES to an interior peak (the
    # reservoir decarboxylates faster than young yeast clears it) then FALLS as the reservoir
    # empties and viable yeast reduces the diacetyl — clearing it below threshold. The
    # interior peak (peak >> final) is the signature a monotone-accumulate pool cannot show.
    r = _run_diacetyl("beer", 25.0, 30.0)
    assert r["diacetyl_peak"] > 5.0 * r["diacetyl_final"]  # a real interior peak, then cleared
    assert r["diacetyl_final"] < 0.1  # cleared below the lager flavour threshold
    # The interior peak is genuinely interior (not the last sample): the series comes back down.
    peak_idx = int(np.argmax(r["diacetyl_mgl"]))
    assert peak_idx < len(r["diacetyl_mgl"]) - 1


def test_cold_strands_diacetyl_with_an_unconverted_reservoir():
    # The other half of the rest: run cold and the diacetyl is STRANDED — it rises to its
    # peak and never comes down (final ≈ peak), because the spontaneous decarb is too slow
    # to empty the α-acetolactate reservoir before the yeast inactivates. The tell-tale is a
    # large UNCONVERTED reservoir left behind, which the warm run consumes. This is the
    # emergent "crash/package too early ⇒ diacetyl stays high" that motivated the 3-pool model.
    cold = _run_diacetyl("beer", 10.0, 30.0)
    warm = _run_diacetyl("beer", 25.0, 30.0)
    # Cold: stranded — final is essentially the peak (never reabsorbed).
    assert cold["diacetyl_final"] == pytest.approx(cold["diacetyl_peak"], rel=0.05)
    # Cold leaves a big unconverted reservoir; warm has consumed almost all of it.
    assert cold["acetolactate_final_mgl"] > 1.0
    assert cold["acetolactate_final_mgl"] > 20.0 * warm["acetolactate_final_mgl"]
    # And the stranded cold diacetyl is far above the cleared warm one.
    assert cold["diacetyl_final"] > 10.0 * warm["diacetyl_final"]


def test_vdk_pathway_conserves_carbon_on_a_default_compiled_run():
    # With the VDK pathway wired into the default medium (D-26), a full compiled ferment
    # still conserves carbon to machine precision — sugar → α-acetolactate → diacetyl + CO2
    # → 2,3-butanediol is fully on the weighted ledger. Non-trivial because the pools
    # actually accumulate.
    r = _run_diacetyl("wine", 20.0, 45.0)
    traj, compiled = r["traj"], r["compiled"]
    schema = compiled.process_set.schema
    f_c = compiled.param_values["biomass_C_fraction"]
    assert_conserved(traj, total_carbon(schema, biomass_carbon_fraction=f_c), label="carbon")
    assert_nonnegative(traj, ("acetolactate", "diacetyl", "butanediol"), atol=1e-12)
    assert r["butanediol_final_mgl"] > 0.0  # the terminal reduction product accumulated


# -- tier propagation ---------------------------------------------------------


def test_excretion_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AcetolactateExcretion()])
    # Speculative form and speculative params ⇒ the reservoir output is speculative.
    assert ps.tier_of("acetolactate") is Tier.SPECULATIVE
    assert ps.tier_of("acetolactate", store.tier_map()) is Tier.SPECULATIVE


def test_decarb_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AcetolactateDecarboxylation()])
    assert ps.tier_of("diacetyl") is Tier.SPECULATIVE
    assert ps.tier_of("diacetyl", store.tier_map()) is Tier.SPECULATIVE


def test_reduction_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [DiacetylReduction()])
    assert ps.tier_of("butanediol") is Tier.SPECULATIVE
    assert ps.tier_of("butanediol", store.tier_map()) is Tier.SPECULATIVE


def test_co2_tier_reflects_the_speculative_decarb_trace(store):
    # The always-on speculative decarboxylation is the FIRST byproduct Process to write the
    # shared CO2 slot (esters/fusels touch S; MLF is disabled unpitched), so it drops the
    # *structural* tier_of("CO2") PLAUSIBLE→SPECULATIVE — the exact D-19 `S` parallel, made
    # explicit here so it can never regress silently (prime directive #1). This is honest:
    # the CO2 pool genuinely holds a speculative decarb trace (real evolved CO2, so it belongs
    # there). Crucially the param-aware tier users SEE is UNCHANGED — already speculative,
    # because the uptake Process itself reads speculative params (E_a_uptake, realised-yield).
    schema = wine_schema()
    tm = store.tier_map()
    core = ProcessSet(schema, [SugarUptakeToEthanolCO2()])  # the validated CO2 producer
    with_decarb = ProcessSet(schema, [SugarUptakeToEthanolCO2(), AcetolactateDecarboxylation()])
    # Structural tier: PLAUSIBLE for the core alone, dropping to SPECULATIVE with the decarb.
    assert core.tier_of("CO2") is Tier.PLAUSIBLE
    assert with_decarb.tier_of("CO2") is Tier.SPECULATIVE
    # Param-aware tier (what users see): SPECULATIVE either way — no headline change from VDK.
    assert core.tier_of("CO2", tm) is Tier.SPECULATIVE
    assert with_decarb.tier_of("CO2", tm) is Tier.SPECULATIVE
