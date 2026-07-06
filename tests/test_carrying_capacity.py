"""Tests for the biomass carrying-capacity cap — the opt-in residual-nitrogen floor (D-30).

In the validated core, :class:`GrowthNitrogenLimited` is the *sole* nitrogen sink and its
only shutoff is a tiny-``K_n`` Monod term, so a wine ferment strips yeast-assimilable nitrogen
(YAN) to ~0 by day ~1.3 *regardless of dose*. That mutes every downstream low-N signal — most
visibly the D-29 H₂S inverse-nitrogen lever, which reads ``N→0`` for every must.

:class:`BiomassCarryingCapacity` is a logistic ``(1 - X/K)`` RateModifier on growth that
saturates biomass below the nitrogen ceiling, leaving a dose-dependent residual of YAN. Because
it scales growth's *whole* contribution by one scalar, ``dN = -f_N·dX`` and the carbon skeleton
draw stay proportional, so nitrogen and carbon still close with the cap on — the nitrogen simply
stays in the ``N`` pool. A residual-N floor is a deliberate DEPARTURE from the validated Coleman
anchor (which caps nothing; ``test_coleman_reconstruction`` pins the zero-residual match at 80
*and* 330 mg N/L), so the cap ships **isolable and disabled by default**: the compile seam
enables it only when a scenario passes ``carrying_capacity_gpl``.

These tests pin: the factor form; that a default (undosed) compile disables it and keeps growth
PLAUSIBLE and the RHS byte-for-byte; that opting in enables it, overrides the cap value, drops
growth's structural tier to speculative, and still conserves carbon + nitrogen; and the emergent
payoff — the H₂S cross-must lever restored and a dose-dependent residual YAN left behind.
"""

import numpy as np
import pytest

from fermentation.core.kinetics import (
    AminoAcidAssimilation,
    BiomassCarryingCapacity,
    EthanolInactivation,
    GrowthNitrogenLimited,
)
from fermentation.core.media import get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon, total_nitrogen

CAP = BiomassCarryingCapacity.name


@pytest.fixture
def full_params():
    # The full default-wine parameter surface: every YAML the wired wine Processes read, so a
    # bare-built wine ProcessSet can evaluate its RHS (mirrors the H₂S isolability fixture).
    base = default_data_dir()
    return load_parameters(
        base / "wine_generic.yaml",
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "keto_acids.yaml",
        base / "hydrogen_sulfide.yaml",
    ).resolve()


def _wine_y0(schema: StateSchema, *, x: float, s: float = 200.0, n: float = 0.1) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": 40.0, "N": n, "T": 293.15, "CO2": 5.0})


# -- metadata + factor form ---------------------------------------------------


def test_metadata():
    m = BiomassCarryingCapacity()
    assert m.name == "biomass_carrying_capacity"
    assert m.tier is Tier.SPECULATIVE
    # Scales growth AND the amino-acid swap (decision D-32): the swap's refunds are anchored to
    # growth's base rate, so the cap must throttle it alongside growth or a near-saturation cap
    # would let the refund outrun the realised draw and create sugar.
    assert m.modifies == (GrowthNitrogenLimited.name, AminoAcidAssimilation.name)
    assert m.reads == ("biomass_carrying_capacity",)


@pytest.mark.parametrize(
    ("x", "expected"),
    [
        (0.0, 1.0),  # empty: no cap
        (1.25, 0.5),  # halfway to K=2.5
        (2.5, 0.0),  # at the cap: growth fully shut
        (3.0, 0.0),  # overshoot clamped to 0 (never negative -> never a source)
        (-1.0, 1.0),  # negative solver excursion clamped to X=0 -> factor 1, not >1
    ],
)
def test_factor_is_clamped_logistic(x, expected):
    schema = wine_schema()
    m = BiomassCarryingCapacity()
    y = _wine_y0(schema, x=x)
    assert m.factor(0.0, y, schema, {"biomass_carrying_capacity": 2.5}) == pytest.approx(expected)


# -- tier: structural drop only when enabled ----------------------------------


def test_cap_drops_growth_output_tier_structurally_only_when_enabled():
    # Growth (PLAUSIBLE) and inactivation (PLAUSIBLE) both touch X, so structural tier_of("X")
    # is PLAUSIBLE. A speculative modifier ON growth folds into X (tier_of extends the
    # least-trustworthy rule to the multiplicative path), so enabling the cap drops X to
    # SPECULATIVE — but ONLY when enabled: a disabled modifier is excluded from tier derivation
    # (the wine-only MLF *tier* isolability argument). Undosed wine keeps growth PLAUSIBLE.
    schema = wine_schema()
    # The cap now also targets the amino-acid swap (D-32), so any set carrying it must include
    # that Process; disable it here (it touches S/N, not X, and is undosed) so this X-tier test
    # is unchanged.
    procs = [GrowthNitrogenLimited(), EthanolInactivation(), AminoAcidAssimilation()]
    off = ProcessSet(schema, procs, modifiers=[BiomassCarryingCapacity()])
    off.disable(CAP)
    off.disable(AminoAcidAssimilation.name)
    on = ProcessSet(schema, procs, modifiers=[BiomassCarryingCapacity()])
    on.disable(AminoAcidAssimilation.name)
    assert off.tier_of("X") is Tier.PLAUSIBLE
    assert on.tier_of("X") is Tier.SPECULATIVE


def test_param_aware_tier_of_x_is_speculative_either_way():
    # Honest no-headline note (the D-26 CO2 / D-27 E parallel): growth already reads the
    # speculative K_s sugar-co-limitation guard, so the param-aware tier of X is SPECULATIVE
    # with or without the cap. The cap's structural drop above is the isolability signal; the
    # user-facing param-aware tier does not change.
    tier_map = load_parameters(default_data_dir() / "wine_generic.yaml").tier_map()
    schema = wine_schema()
    # Include the amino-acid swap (the cap now targets it, D-32) but disable it — it touches
    # S/N, not X, so the X tier is unaffected either way.
    procs = [GrowthNitrogenLimited(), AminoAcidAssimilation()]
    off = ProcessSet(schema, procs, modifiers=[BiomassCarryingCapacity()])
    off.disable(CAP)
    off.disable(AminoAcidAssimilation.name)
    on = ProcessSet(schema, procs, modifiers=[BiomassCarryingCapacity()])
    on.disable(AminoAcidAssimilation.name)
    assert off.tier_of("X", tier_map) is Tier.SPECULATIVE
    assert on.tier_of("X", tier_map) is Tier.SPECULATIVE


# -- isolability at the derivative level --------------------------------------


def test_disabled_cap_equals_the_uncapped_rhs_exactly(full_params):
    # THE byte-for-byte claim: a disabled cap contributes nothing. Compare a wine set with the
    # cap DISABLED against the same set with the cap ENABLED at K=∞ (factor = 1 - X/∞ = 1.0
    # exactly in float). Both reduce to the uncapped core RHS, so the difference is exactly 0.0
    # on every column and every state — stronger than a tolerance (the compile seam gives every
    # undosed run this disabled path).
    schema = get_medium("wine").schema
    off = get_medium("wine").build_process_set()
    off.disable(CAP)
    on_inf = get_medium("wine").build_process_set()  # enabled, K=∞ -> factor 1.0
    params_inf = {**full_params, "biomass_carrying_capacity": np.inf}
    for x, s, n in [(0.5, 240.0, 0.15), (2.0, 120.0, 0.02), (3.5, 40.0, 0.0)]:
        y = _wine_y0(schema, x=x, s=s, n=n)
        diff = on_inf.total_derivatives(0.0, y, params_inf) - off.total_derivatives(
            0.0, y, full_params
        )
        assert np.max(np.abs(diff)) == 0.0


def test_enabled_cap_scales_growth_by_the_logistic_factor(full_params):
    # When enabled at a biting K, the cap scales growth's WHOLE contribution by (1 - X/K).
    # Isolate growth (the full-set X column also carries the inactivation X→X_dead flux, which
    # the cap does NOT scale): on a growth-only set, every column growth touches — X, N, and the
    # sugar draw — scales by exactly the same factor, so nitrogen/carbon proportions are kept.
    schema = wine_schema()
    k = 2.5
    # Include the amino-acid swap (the cap now targets it, D-32) but disable it — undosed, it
    # contributes nothing, so growth's scaling is isolated exactly as before.
    procs = [GrowthNitrogenLimited(), AminoAcidAssimilation()]
    off = ProcessSet(schema, procs, modifiers=[BiomassCarryingCapacity()])
    off.disable(CAP)
    off.disable(AminoAcidAssimilation.name)
    on = ProcessSet(schema, procs, modifiers=[BiomassCarryingCapacity()])
    on.disable(AminoAcidAssimilation.name)
    params_on = {**full_params, "biomass_carrying_capacity": k}
    x = 1.5
    y = _wine_y0(schema, x=x, s=200.0, n=0.2)
    d_off = off.total_derivatives(0.0, y, full_params)
    d_on = on.total_derivatives(0.0, y, params_on)
    assert d_off[schema.slice("X").start] > 0.0  # growth is running
    for var in ("X", "N"):
        i = schema.slice(var).start
        assert d_on[i] == pytest.approx(d_off[i] * (1.0 - x / k))


# -- behaviour through the compile seam ---------------------------------------


def _run(yan_mgl: float, *, cap_gpl: float | None = None, days: float = 21.0):
    """Compile+integrate a wine ferment at the given YAN; opt into the carrying-capacity cap
    when ``cap_gpl`` is set. Returns (trajectory, compiled)."""
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if cap_gpl is not None:
        initial["carrying_capacity_gpl"] = cap_gpl
    scenario = Scenario(
        name=f"wine-cap-{yan_mgl:.0f}-{cap_gpl or 0:.1f}",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
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


def _final_mgl(traj, name: str) -> float:
    return float(traj.series(name)[-1]) * 1000.0


def _final_produced_h2s_mgl(traj) -> float:
    # Cumulative H₂S produced = residual (h2s) + swept-to-gas (h2s_gas); the D-42 CO2-stripping
    # sink splits produced between the two pools, so the cross-must LEVER (a produced quantity)
    # reads the sum, not the µg/L residual pool alone.
    return float(traj.series("h2s")[-1] + traj.series("h2s_gas")[-1]) * 1000.0


def test_default_compile_disables_the_cap():
    # No carrying_capacity_gpl ⇒ the modifier is present (wired into wine) but DISABLED, so an
    # undosed run is the validated core.
    _, compiled = _run(80.0)
    assert CAP in compiled.process_set
    assert not compiled.process_set.is_enabled(CAP)


def test_opt_in_enables_and_overrides_the_cap_value():
    _, compiled = _run(80.0, cap_gpl=2.2)
    assert compiled.process_set.is_enabled(CAP)
    # The scenario value overrides the 2.5 g/L YAML reference (lets a demonstration sweep K).
    assert compiled.param_values["biomass_carrying_capacity"] == pytest.approx(2.2)


def test_negative_cap_raises_like_the_other_initial_keys():
    # A negative cap is a typo, not "opt out" — it must raise loudly (the _nonneg gate every
    # other initial key gets), not silently disable.
    scenario = Scenario(
        name="wine-cap-neg",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25, "carrying_capacity_gpl": -1.0},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )
    with pytest.raises(ValueError, match="carrying_capacity_gpl"):
        compile_scenario(scenario, strict=True)


def test_carbon_and_nitrogen_close_with_the_cap_on():
    # The crux: scaling growth's whole contribution keeps dN = -f_N·dX and the carbon-skeleton
    # draw proportional, so BOTH atom balances still close to solver tolerance with the cap
    # biting — the residual nitrogen is left in the N pool, not created or destroyed.
    traj, compiled = _run(300.0, cap_gpl=2.5)
    f_c = compiled.param_values["biomass_C_fraction"]
    f_n = compiled.param_values["biomass_N_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj, total_nitrogen(compiled.schema, biomass_nitrogen_fraction=f_n), label="nitrogen"
    )


def test_cap_leaves_dose_dependent_residual_nitrogen():
    # Core (cap off) strips YAN to ~0 at EVERY dose — the documented gap. With the cap on, a
    # dose-dependent residual survives: a high-YAN must ends well above a low-YAN one, and the
    # low-YAN must still bottoms out near zero (the correct clinical picture).
    core_lo = _final_mgl(_run(80.0)[0], "N")
    core_hi = _final_mgl(_run(300.0)[0], "N")
    assert core_lo < 1.0 and core_hi < 1.0  # muted: both stripped to ~0

    cap_lo = _final_mgl(_run(80.0, cap_gpl=2.5)[0], "N")
    cap_hi = _final_mgl(_run(300.0, cap_gpl=2.5)[0], "N")
    assert cap_lo < 10.0  # low-YAN must still (nearly) exhausts YAN
    assert cap_hi > cap_lo + 10.0  # high-YAN must leaves a real residual


def test_cap_restores_the_h2s_cross_must_lever():
    # The headline payoff (D-29 → D-30): with N stripped to ~0 at every dose, the core H₂S
    # inverse-N gate barely separates musts (muted lever). The cap leaves dose-dependent
    # residual N, so the gate stays partly shut at high YAN ⇒ H₂S is monotonically ordered by
    # dose and its span widens materially versus the core. Asserted as ordering + ratio (not
    # brittle absolute values — speculative-on-speculative).
    core = {yan: _final_produced_h2s_mgl(_run(yan)[0]) for yan in (80.0, 150.0, 300.0)}
    cap = {yan: _final_produced_h2s_mgl(_run(yan, cap_gpl=2.5)[0]) for yan in (80.0, 150.0, 300.0)}

    # Monotone in dose: less nitrogen ⇒ more sulfide.
    assert cap[80.0] > cap[150.0] > cap[300.0]

    core_span = core[80.0] / core[300.0]
    cap_span = cap[80.0] / cap[300.0]
    assert cap_span > 1.4  # a real cross-must lever...
    assert cap_span > core_span  # ...materially wider than the muted core


def test_opt_in_wine_still_reaches_dryness():
    # Nice-to-have (NOT a §2.2 gate — the benchmark tests the default disabled core): a capped
    # wine still ferments to dryness within the window. Less biomass slows the tail, but the
    # per-cell uptake keeps going, so an opted-in must does not stall short.
    traj, _ = _run(80.0, cap_gpl=2.5)
    s = np.asarray(traj.series("S"))
    total_sugar = s if s.ndim == 1 else s.sum(axis=0)
    assert float(total_sugar[-1]) < 2.0
