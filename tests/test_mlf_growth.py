"""Malolactic *growth* — the deferred MLF-growth beat (decision D-23).

Ranked headline-first. The acceptance payoff is ``test_growth_accelerates_conversion``: a
co-inoculated, amino-acid-fed must deacidifies **faster** than an otherwise-identical
fixed-``X_mlf`` control, because :class:`MalolacticGrowth` multiplies the bacterial catalyst
and :class:`MalolacticConversion` is linear in it — autocatalysis that *emerges*, and that
**vanishes** the moment the growth Process is removed (the fail-first control). The rest pin
the nitrogen-anchored / carbon-shortfall-from-sugar stoichiometry and its exact carbon+nitrogen
closure (with ``X_mlf`` now weighted as real bacterial biomass — decision D-38), the
``touches`` contract, the compile-seam gate (growth needs BOTH a pitch AND amino acids), the
guards (no catalyst / no fuel / no sugar ⇒ zero), that it never creates sugar, the ethanol-wall
arrest that makes growth a co-inoculation phenomenon, and the explicit ``speculative`` tier.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.kinetics.malolactic import MalolacticGrowth
from fermentation.core.media import wine_schema
from fermentation.core.state import StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)


@pytest.fixture
def pset():
    """Real wine kinetic params (incl. the MLF + growth set) + the shared pKa set."""
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    return pset.resolve()


def _wine_scenario(*, days: float = 14.0, **initial_extra) -> Scenario:
    initial: dict[str, float] = {
        "brix": 24.0,
        "yan_mgl": 300.0,
        "pitch_gpl": 0.25,
        "tartaric_gpl": 3.0,
        "malic_gpl": 3.0,
        "initial_ph": 3.4,
    }
    initial.update(initial_extra)
    return Scenario(
        name="wine-mlf-growth",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=days,
    )


def _run(*, days: float = 14.0, n_eval: int = 400, **initial_extra: float):
    compiled = compile_scenario(_wine_scenario(days=days, **initial_extra), strict=True)
    t_eval = np.linspace(0.0, days * 24.0, n_eval)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    return compiled, traj


def _malic_at_day(traj, day: float) -> float:
    return float(np.interp(day * 24.0, traj.t, traj.series("malic")))


# -- 1. ACCEPTANCE (fail-first): growth accelerates malolactic conversion ----------


def test_growth_accelerates_conversion():
    # Fail-first / apples-to-apples: the SAME co-inoculated, amino-acid-fed scenario, run with
    # MalolacticGrowth enabled vs the identical run with the growth Process disabled (the
    # fixed-X_mlf control). Both dose amino acids, so the yeast-side D-32/D-33 perturbation is
    # identical and the ONLY difference is bacterial growth — isolating the autocatalysis.
    grow = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.05, amino_acids_gpl=2.0), strict=True)
    assert grow.process_set.is_enabled("malolactic_growth")
    # Disable the benign senescence baseline (MLF v2, D-41) in BOTH runs: it erodes X_mlf on its own
    # slow timescale regardless of growth, which would confound this GROWTH-isolation contrast (the
    # fixed-X_mlf control below could no longer be exactly constant). Senescence is exercised on its
    # own in test_malolactic; here we isolate growth's autocatalysis (the control disables growth).
    grow.process_set.disable("malolactic_senescence")
    t_eval = np.linspace(0.0, 14.0 * 24.0, 400)
    t_on = simulate(grow.process_set, grow.param_values, grow.y0, grow.t_span_h, t_eval=t_eval)

    fixed = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.05, amino_acids_gpl=2.0), strict=True)
    fixed.process_set.disable("malolactic_growth")  # fixed-X_mlf control
    fixed.process_set.disable("malolactic_senescence")  # ...held truly fixed (no baseline decay)
    t_off = simulate(fixed.process_set, fixed.param_values, fixed.y0, fixed.t_span_h, t_eval=t_eval)

    # X_mlf multiplies several-fold under growth; constant in the control.
    assert t_on.series("X_mlf")[-1] > 2.0 * t_on.series("X_mlf")[0]
    assert t_off.series("X_mlf")[-1] == pytest.approx(t_off.series("X_mlf")[0], rel=1e-9)

    # The autocatalytic payoff: at day 3 the growing must has converted far more malate — and the
    # gap vanishes when the growth Process is removed (fail-first).
    assert _malic_at_day(t_on, 3.0) < 0.5 * _malic_at_day(t_off, 3.0)


# -- 2. CONSERVATION: carbon AND nitrogen close over a growing run -----------------


def test_carbon_and_nitrogen_conserved_over_a_growing_run():
    compiled, traj = _run(mlf_pitch_gpl=0.05, amino_acids_gpl=2.0)
    fc = compiled.param_values["biomass_C_fraction"]
    fn = compiled.param_values["biomass_N_fraction"]
    carbon = total_carbon(compiled.schema, biomass_carbon_fraction=fc)
    nitrogen = total_nitrogen(compiled.schema, biomass_nitrogen_fraction=fn)
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-6, label="total carbon (MLF growth)")
    assert_conserved(traj, nitrogen, rtol=1e-6, atol=1e-6, label="total nitrogen (MLF growth)")
    assert_nonnegative(traj, ("X_mlf", "amino_acids", "malic", "lactic"))


# -- 3. touches (X_mlf, amino_acids, S) only; never N ------------------------------


def test_touches_only_x_mlf_amino_acids_sugar(params):
    proc = MalolacticGrowth()
    assert set(proc.touches) == {"X_mlf", "amino_acids", "S"}
    assert "N" not in proc.touches  # nitrogen comes from amino acids, not the ammonium pool

    schema = wine_schema()
    y = schema.pack(
        {"X": 1.0, "S": [120.0], "E": 20.0, "N": 0.05, "T": 293.15, "CO2": 0.0,
         "X_mlf": 0.1, "amino_acids": 2.0, "malic": 3.0, "tartaric": 3.0}
    )  # fmt: skip
    d = proc.derivatives(0.0, y, schema, params)
    moved = {name for name in schema.names if float(d[schema.slice(name)][0]) != 0.0}
    assert moved == {"X_mlf", "amino_acids", "S"}


# -- 4. compile-seam gate: keyed on the amino-acid feature, NOT the pitch ----------


@pytest.mark.parametrize(
    "mlf,aa,expected",
    [(0.0, 0.0, False), (0.05, 0.0, False), (0.0, 2.0, True), (0.05, 2.0, True)],
)
def test_growth_enabled_iff_amino_acids_dosed(mlf, aa, expected):
    # The gate keys on the FEATURE (amino-acid-fed bacterial growth ⇒ amino_acids dosed), the same
    # gate as the swap/re-route — NOT on the pitch. Pitched-but-not-aa (every D-23/D-31 test) stays
    # disabled (the tier-isolability guarantee); aa-dosed enables it regardless of the pitch, and
    # the X_mlf<=0 guard keeps it inert until bacteria are actually present (decision D-38).
    extra: dict[str, float] = {}
    if mlf:
        extra["mlf_pitch_gpl"] = mlf
    if aa:
        extra["amino_acids_gpl"] = aa
    cs = compile_scenario(_wine_scenario(**extra), strict=True)
    assert cs.process_set.is_enabled("malolactic_growth") is expected


def test_mid_run_pitch_growth_is_emergently_gated_by_ethanol():
    # Growth is NOT compile-gated on co-inoculation. With amino acids dosed it is enabled from t0,
    # inert until a mid-run pitch_mlf adds X_mlf (the X_mlf<=0 guard), and whether the pitched
    # bacteria then GROW is left to the emergent environmental gate — as MalolacticConversion trusts
    # its ethanol gate. So an early pitch (low ethanol) grows; a late post-AF pitch (past the
    # ~110 g/L O. oeni wall) does not. Co-inoculation dominance is emergent, not hard-coded (D-38).
    from fermentation.scenario.schema import Intervention

    # Isolate the GROWTH signal as a control DIFFERENCE (growth-on minus growth-off) at each pitch
    # day. Since MLF v2 (D-41) the benign senescence baseline also moves X_mlf, and it is re-enabled
    # by the mid-run pitch_mlf (part of _MLF_GATED_PROCESSES), so it can't simply be toggled off for
    # the run — but it acts in BOTH the growth-on and growth-off arms, so it cancels in their
    # difference. Growth is amino-acid-gated (NOT in _MLF_GATED_PROCESSES), so disabling it survives
    # the pitch. The residual = growth's net contribution, cleanly separated from senescence.
    def _end_x_mlf(day: float, *, growth: bool) -> float:
        base = _wine_scenario(days=21.0, amino_acids_gpl=2.0)
        sc = Scenario(
            name=base.name, medium=base.medium, strain=base.strain, initial=base.initial,
            temperature_schedule=base.temperature_schedule, duration_days=base.duration_days,
            interventions=[Intervention(day=day, action="pitch_mlf", params={"pitch_gpl": 0.1})],
        )  # fmt: skip
        cs = compile_scenario(sc, strict=True)
        assert cs.process_set.is_enabled("malolactic_growth")  # enabled by aa, regardless of pitch
        if not growth:
            cs.process_set.disable("malolactic_growth")  # survives pitch_mlf (aa-gated, not pitch)
        tr = cs.run(t_eval=np.linspace(0.0, 21.0 * 24.0, 500))
        return float(tr.series("X_mlf")[-1])

    early_gain = _end_x_mlf(2.0, growth=True) - _end_x_mlf(
        2.0, growth=False
    )  # E ~64 g/L (below wall)
    late_gain = _end_x_mlf(8.0, growth=True) - _end_x_mlf(
        8.0, growth=False
    )  # E ~115 g/L (past wall)
    # Early pitch: growth measurably ADDS biomass over the no-growth baseline (well above solver
    # noise). Late pitch: ethanol-arrested, so growth adds nothing. The gap is the emergent
    # co-inoculation advantage, not a hard-coded rule — senescence-independent (it cancels out).
    assert early_gain > 1e-4
    assert late_gain == pytest.approx(0.0, abs=1e-5)
    assert early_gain > late_gain  # the emergent early-vs-late growth advantage


# -- 5. derivative-level guards + never-creates-sugar ------------------------------


def _growing_state(schema: StateSchema, params: Mapping[str, float], **overrides) -> np.ndarray:
    base: dict[str, float | list[float]] = {
        "X": 1.0, "S": [120.0], "E": 10.0, "N": 0.05, "T": 293.15, "CO2": 0.0,
        "X_mlf": 0.1, "amino_acids": 2.0, "malic": 3.0, "tartaric": 3.0,
    }  # fmt: skip
    base.update(overrides)
    return schema.pack(base)


@pytest.mark.parametrize("zeroed", ["X_mlf", "amino_acids", "S"])
def test_no_catalyst_no_fuel_no_sugar_gives_zero(params, zeroed):
    # Each guard independently zeroes the contribution BEFORE the pH solve (no catalyst, no
    # amino-acid fuel, or no sugar — the last so the carbon shortfall never targets an empty S).
    schema = wine_schema()
    override = {zeroed: [0.0] if zeroed == "S" else 0.0}
    y = _growing_state(schema, params, **override)
    d = MalolacticGrowth().derivatives(0.0, y, schema, params)
    assert not np.any(d)


def test_growth_never_creates_sugar(params):
    # The carbon shortfall is DRAWN from sugar, so growth's contribution to dS is <= 0 always
    # (the structurally-positive shortfall coefficient f_C - f_N*y_C/y_N; decision D-38).
    schema = wine_schema()
    y = _growing_state(schema, params)
    d = MalolacticGrowth().derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("S")][0]) < 0.0
    # The CARBON drawn from sugar (mass rate × glucose carbon fraction) is a positive shortfall
    # strictly below the biomass carbon built — arginine's carbon covers the rest.
    fc = params["biomass_C_fraction"]
    dx = float(d[schema.slice("X_mlf")][0])
    sugar_carbon_drawn = -float(d[schema.slice("S")][0]) * carbon_mass_fraction("glucose")
    assert 0.0 < sugar_carbon_drawn < fc * dx


def test_stoichiometry_closes_at_the_derivative_level(params):
    # Carbon: X_mlf gains f_C*dX = amino-acid carbon + sugar shortfall. Nitrogen: X_mlf gains
    # f_N*dX = amino-acid nitrogen. Checked directly on the derivatives (the closure the run test
    # integrates), using arginine's fractions and the biomass fractions X_mlf is weighted at.
    schema = wine_schema()
    y = _growing_state(schema, params)
    d = MalolacticGrowth().derivatives(0.0, y, schema, params)
    fc, fn = params["biomass_C_fraction"], params["biomass_N_fraction"]
    y_c = carbon_mass_fraction(AMINO_ACID_SPECIES)
    y_n = nitrogen_mass_fraction(AMINO_ACID_SPECIES)

    dx = float(d[schema.slice("X_mlf")][0])
    d_aa = float(d[schema.slice("amino_acids")][0])
    d_s = float(d[schema.slice("S")][0])
    # nitrogen: biomass gain == amino-acid nitrogen lost
    assert fn * dx == pytest.approx(-d_aa * y_n, rel=1e-12)
    # carbon: biomass gain == amino-acid carbon lost + sugar carbon lost
    assert fc * dx == pytest.approx(-d_aa * y_c - d_s * carbon_mass_fraction("glucose"), rel=1e-12)


# -- 6. the ethanol wall makes growth a co-inoculation phenomenon ------------------


def test_ethanol_wall_arrests_growth(params):
    # High ethanol (a post-AF must) collapses the environmental gate, so bacteria cannot build up
    # even with amino acids present — the emergent reason MLF-growth is a co-inoculation mode.
    schema = wine_schema()
    low_e = MalolacticGrowth().derivatives(
        0.0, _growing_state(schema, params, E=10.0), schema, params
    )
    high_e = MalolacticGrowth().derivatives(
        0.0, _growing_state(schema, params, E=130.0), schema, params
    )
    assert float(low_e[schema.slice("X_mlf")][0]) > 0.0
    assert float(high_e[schema.slice("X_mlf")][0]) == pytest.approx(0.0, abs=1e-12)


# -- 7. tier: growth is speculative and caps its outputs ---------------------------


def test_growth_is_speculative_and_caps_x_mlf_tier():
    assert MalolacticGrowth.tier is Tier.SPECULATIVE
    # Discriminating: an amino-acid-dosed but UN-pitched run disables MalolacticConversion (no
    # pitch), so growth is the ONLY enabled Process touching X_mlf. Its speculative tier therefore
    # is what drops X_mlf's output tier — proven by toggling just the growth Process.
    cs = compile_scenario(_wine_scenario(amino_acids_gpl=2.0), strict=True)
    assert not cs.process_set.is_enabled("malolactic_conversion")
    assert cs.process_set.is_enabled("malolactic_growth")
    tm = cs.parameters.tier_map()
    assert cs.process_set.tier_of("X_mlf", tm) is Tier.SPECULATIVE  # growth caps it
    cs.process_set.disable("malolactic_growth")
    assert cs.process_set.tier_of("X_mlf", tm) is Tier.VALIDATED  # nothing enabled touches it
