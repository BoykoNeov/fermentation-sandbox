"""*Brettanomyces* volatile-phenol spoilage — the mixed-culture beat (decision D-40, pt1).

Ranked headline-first. The acceptance payoff is
``test_headline_brett_raises_ethylphenols_emergently``: a **POF-negative** wine (no yeast
decarboxylase) dosed with hydroxycinnamic precursors accumulates the ``ethylphenols`` off-aroma
**only when Brett is pitched** — Brett carries both the decarboxylase and the reductase, so it
spoils normal wine unaided, and the 4-EP signal **vanishes** in the otherwise-identical no-Brett
control (the fail-first control). ``test_so2_suppresses`` and ``test_rack_suppresses`` are the two
winemaking levers (molecular SO₂ arrests the metabolism; racking draws Brett off the lees). The
rest pin the decarboxylase / reductase stoichiometry and its exact carbon closure, the ``touches``
contracts, the compile-seam pitch gate + tier isolability (an unpitched run is byte-for-byte the
validated core and the phenol slots keep their VALIDATED tier), the guards, the warmer-than-MLF
temperature optimum, and the explicit ``speculative`` tier.
"""

import numpy as np
import pytest

from fermentation.core.chemistry import (
    M_ETHYLGUAIACOL,
    M_ETHYLPHENOL,
    M_VINYLGUAIACOL,
    M_VINYLPHENOL,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.brett import (
    BrettDeath,
    BrettDecarboxylation,
    BrettEthanolToxicity,
    BrettGrowth,
    BrettVinylphenolReduction,
    YeastPOFDecarboxylation,
    brett_environmental_gate,
    brett_ethanol_survival_factor,
)
from fermentation.core.kinetics.malolactic import cardinal_temperature_factor
from fermentation.core.media import wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario.schema import Intervention
from fermentation.units.convert import mgl_to_gpl
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)


@pytest.fixture
def pset():
    """Real wine kinetic params (incl. the Brett set) + the shared pKa set."""
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    return pset.resolve()


@pytest.fixture
def schema():
    return wine_schema()


def _state(schema: StateSchema, **overrides: float | list[float]) -> FloatArray:
    """A wine state with sane defaults, overridable per-slot (T defaults to 20 °C)."""
    values: dict[str, float | list[float]] = {
        "X": 0.1,
        "S": [200.0],
        "E": 90.0,
        "N": 0.0,
        "T": 293.15,
        "CO2": 0.0,
        "tartaric": 3.0,
        "malic": 2.0,
    }
    values.update(overrides)
    return schema.pack(values)


# -- scenario helpers ---------------------------------------------------------


def _wine_scenario(
    *, days: float = 40.0, interventions: list[Intervention] | None = None, **initial_extra: float
) -> Scenario:
    initial: dict[str, float] = {
        "brix": 22.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.2,
        "tartaric_gpl": 3.0,
        "malic_gpl": 2.0,
        "initial_ph": 3.5,
    }
    initial.update(initial_extra)
    return Scenario(
        name="wine-brett",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=interventions or [],
        duration_days=days,
    )


def _run(
    *,
    days: float = 40.0,
    n_eval: int = 300,
    interventions: list[Intervention] | None = None,
    **initial_extra: float,
):
    compiled = compile_scenario(
        _wine_scenario(days=days, interventions=interventions, **initial_extra), strict=True
    )
    t_eval = np.linspace(0.0, days * 24.0, n_eval)
    traj = compiled.run(t_eval=t_eval)
    return compiled, traj


# -- 1. HEADLINE: emergent, Brett-gated volatile phenols ----------------------


def test_headline_brett_raises_ethylphenols_emergently():
    """POF- wine + hydroxycinnamics: 4-EP rises **only** when Brett is pitched.

    The canonical funk mechanism. A dosed Brett culture (both enzymes) takes must
    hydroxycinnamics all the way to ethylphenols; the identical must with **no Brett** keeps
    ethylphenols at exactly 0 — no yeast decarboxylase in the POF-negative default (the pt4
    opt-in POF+ strain is off). The signal is emergent from the pitch, not scripted.
    """
    _, brett = _run(hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.3)
    _, none = _run(hydroxycinnamic_gpl=0.1)  # POF- wine, no Brett

    ep_brett = brett.series("ethylphenols")
    ep_none = none.series("ethylphenols")

    assert ep_brett[-1] > 10.0 * (ep_brett[0] + 1e-9)  # a clear rise
    assert float(np.max(np.abs(ep_none))) == 0.0  # no Brett ⇒ no volatile phenols at all
    # And the precursor is genuinely consumed into the pathway (not merely diluted).
    hc = brett.series("hydroxycinnamics")
    assert hc[-1] < 0.9 * hc[0]


def test_pitch_brett_post_af_at_high_ethanol():
    """The ``pitch_brett`` verb: Brett contaminates a *finished*, full-strength wine.

    The honest Brett framing — a post-AF/barrel spoiler, not a t0 co-inoculant. Exercises the
    intervention verb end to end (mutate ``X_brett`` + reconfigure the Process set at the
    breakpoint) and, in the same stroke, asserts Brett's **ethanol tolerance** at the integration
    level: ethylphenols is *exactly* 0 before the pitch (the verb gates correctly) and rises after,
    at ~14 % ABV ethanol — the property that would silently die if anyone re-added an MLF-style
    ethanol wall to the Brett gate.
    """
    pitch_day = 25.0
    _, traj = _run(
        days=60.0,
        hydroxycinnamic_gpl=0.1,  # precursor present in the must from t0
        interventions=[
            Intervention(day=pitch_day, action="pitch_brett", params={"pitch_gpl": 0.3})
        ],
    )
    t_days = traj.t / 24.0
    ep = traj.series("ethylphenols")

    # Fermentation has finished by the pitch: ethanol is at full strength (~14 % ABV ≈ 110 g/L).
    e_at_pitch = float(np.interp(pitch_day * 24.0, traj.t, traj.series("E")))
    assert e_at_pitch > 100.0

    before = ep[t_days <= pitch_day]
    assert float(np.max(np.abs(before))) == 0.0  # verb gates: nothing until the pitch
    assert ep[-1] > 1e-3  # Brett makes phenols at full-strength-wine ethanol (no ethanol wall)


def test_so2_suppresses_ethylphenols():
    """A molecular-SO₂ dose arrests Brett metabolism ⇒ far less 4-EP (the winemaker's lever)."""
    _, base = _run(hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.3)
    _, sulfited = _run(hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.3, so2_total_mgl=60.0)
    assert sulfited.series("ethylphenols")[-1] < 0.1 * base.series("ethylphenols")[-1]


def test_rack_suppresses_ethylphenols():
    """Racking draws Brett off the lees (X_brett → 0) ⇒ phenol production halts at the rack."""
    rack = [Intervention(day=5.0, action="rack", params={"fraction": 1.0})]
    _, base = _run(hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.3)
    _, racked = _run(hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.3, interventions=rack)
    assert racked.series("X_brett")[-1] == pytest.approx(0.0, abs=1e-12)
    # Only what accrued in the first 5 days survives — well below the un-racked endpoint.
    assert racked.series("ethylphenols")[-1] < 0.2 * base.series("ethylphenols")[-1]


# -- 2. conservation ----------------------------------------------------------


def test_carbon_closes_through_the_chain():
    """total_carbon closes to tolerance through precursor → intermediate → product + CO₂."""
    compiled, traj = _run(hydroxycinnamic_gpl=0.15, brett_pitch_gpl=0.4)
    fn = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.param_values["biomass_C_fraction"]
    )
    assert_conserved(traj, fn, label="carbon")
    assert_nonnegative(traj, ("hydroxycinnamics", "vinylphenols", "ethylphenols"))


def test_carbon_closes_through_the_ferulic_branch():
    """The D-55 ferulic branch, end-to-end: total_carbon closes with BOTH branches active.

    Dosing ferulic_acid_gpl alongside hydroxycinnamic_gpl exercises the full pipeline (scenario
    dosing -> decarboxylation -> reduction) for the second precursor, not just the per-Process
    derivatives() calls above — a genuine wiring check the unit-level tests can't catch (e.g. a
    typo in the compile-seam initial-condition key, or a slot the reduction step forgot to drain).
    """
    compiled, traj = _run(hydroxycinnamic_gpl=0.15, ferulic_acid_gpl=0.08, brett_pitch_gpl=0.4)
    fn = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.param_values["biomass_C_fraction"]
    )
    assert_conserved(traj, fn, label="carbon (ferulic branch)")
    assert_nonnegative(
        traj,
        (
            "hydroxycinnamics",
            "vinylphenols",
            "ethylphenols",
            "ferulic_acid",
            "vinylguaiacols",
            "ethylguaiacols",
        ),
    )
    # Both branches actually ran and reduced through to their terminal readouts.
    assert float(traj.series("ethylphenols")[-1]) > 0.0
    assert float(traj.series("ethylguaiacols")[-1]) > 0.0
    assert float(traj.series("ferulic_acid")[-1]) < 0.08  # genuinely consumed


# -- 3. per-Process stoichiometry + touches -----------------------------------


def test_decarboxylation_stoichiometry_and_touches(schema, params):
    """hydroxycinnamics down; vinylphenols + CO₂ up; carbon flux sums to 0; touches honoured."""
    y = _state(schema, X_brett=0.3, hydroxycinnamics=0.1)
    d = BrettDecarboxylation().derivatives(0.0, y, schema, params)

    assert d[schema.slice("hydroxycinnamics")][0] < 0.0
    assert d[schema.slice("vinylphenols")][0] > 0.0
    assert d[schema.slice("CO2")][0] > 0.0
    assert d[schema.slice("ferulic_acid")][0] == 0.0  # undosed branch is exactly inert
    assert d[schema.slice("vinylguaiacols")][0] == 0.0
    # Carbon flux across the three touched slots must cancel (9 C = 8 C + 1 C).
    c_flux = (
        d[schema.slice("hydroxycinnamics")][0] * carbon_mass_fraction("p_coumaric_acid")
        + d[schema.slice("vinylphenols")][0] * carbon_mass_fraction("vinylphenol")
        + d[schema.slice("CO2")][0] * carbon_mass_fraction("CO2")
    )
    assert c_flux == pytest.approx(0.0, abs=1e-12)
    assert set(BrettDecarboxylation.touches) == {
        "hydroxycinnamics",
        "vinylphenols",
        "ferulic_acid",
        "vinylguaiacols",
        "CO2",
    }


def test_decarboxylation_ferulic_branch_stoichiometry(schema, params):
    """The D-55 ferulic branch: ferulic_acid down, vinylguaiacols + CO₂ up, carbon flux sums to 0.

    A genuinely distinct precursor (10 C) from p-coumaric's ``hydroxycinnamics`` (9 C), so this
    checks its own carbon-closing reaction (10 = 9 + 1), both alone and composed with the
    p-coumaric branch running simultaneously (both draw the same catalyst/gate, independently).
    """
    y = _state(schema, X_brett=0.3, ferulic_acid=0.05)  # no hydroxycinnamics dosed
    d = BrettDecarboxylation().derivatives(0.0, y, schema, params)

    assert d[schema.slice("ferulic_acid")][0] < 0.0
    assert d[schema.slice("vinylguaiacols")][0] > 0.0
    assert d[schema.slice("CO2")][0] > 0.0
    assert d[schema.slice("hydroxycinnamics")][0] == 0.0  # undosed branch is exactly inert
    assert d[schema.slice("vinylphenols")][0] == 0.0
    c_flux = (
        d[schema.slice("ferulic_acid")][0] * carbon_mass_fraction("ferulic_acid")
        + d[schema.slice("vinylguaiacols")][0] * carbon_mass_fraction("vinylguaiacol")
        + d[schema.slice("CO2")][0] * carbon_mass_fraction("CO2")
    )
    assert c_flux == pytest.approx(0.0, abs=1e-12)

    # Both branches active simultaneously: each is independent (same catalyst/gate, distinct pools).
    y_pcoumaric_only = _state(schema, X_brett=0.3, hydroxycinnamics=0.1)
    d_pcoumaric_only = BrettDecarboxylation().derivatives(0.0, y_pcoumaric_only, schema, params)

    y_both = _state(schema, X_brett=0.3, hydroxycinnamics=0.1, ferulic_acid=0.05)
    d_both = BrettDecarboxylation().derivatives(0.0, y_both, schema, params)
    assert d_both[schema.slice("hydroxycinnamics")][0] < 0.0
    assert d_both[schema.slice("ferulic_acid")][0] < 0.0
    assert d_both[schema.slice("vinylphenols")][0] > 0.0
    assert d_both[schema.slice("vinylguaiacols")][0] > 0.0
    # The p-coumaric branch's own rate is unaffected by the ferulic branch running alongside it
    # (independent Monod terms sharing only the catalyst/gate `activity`, not each other's pool).
    assert d_both[schema.slice("hydroxycinnamics")][0] == pytest.approx(
        d_pcoumaric_only[schema.slice("hydroxycinnamics")][0], rel=1e-9
    )


def test_reduction_is_carbon_neutral_transfer(schema, params):
    """vinylphenol → ethylphenol is mole-for-mole (C8 → C8): moles produced == moles consumed."""
    y = _state(schema, X_brett=0.3, vinylphenols=0.02)
    d = BrettVinylphenolReduction().derivatives(0.0, y, schema, params)

    dvp = d[schema.slice("vinylphenols")][0]
    dep = d[schema.slice("ethylphenols")][0]
    assert dvp < 0.0 and dep > 0.0
    moles_consumed = -dvp / M_VINYLPHENOL
    moles_produced = dep / M_ETHYLPHENOL
    assert moles_produced == pytest.approx(moles_consumed, rel=1e-12)
    assert d[schema.slice("vinylguaiacols")][0] == 0.0  # undosed branch is exactly inert
    assert d[schema.slice("ethylguaiacols")][0] == 0.0
    assert set(BrettVinylphenolReduction.touches) == {
        "vinylphenols",
        "ethylphenols",
        "vinylguaiacols",
        "ethylguaiacols",
    }


def test_reduction_ferulic_branch_is_carbon_neutral_transfer(schema, params):
    """vinylguaiacol → ethylguaiacol is mole-for-mole (C9 → C9) — the D-55 ferulic-branch analogue.

    Tchobanov et al. 2008 directly confirm Brett's vinylphenol reductase acts on 4-vinylguaiacol
    too, reusing the same k_brett_reduction (no differential rate sourced between the two
    substrates — a documented simplification, unlike the decarboxylase branches' sourced ratio).
    """
    y = _state(schema, X_brett=0.3, vinylguaiacols=0.02)
    d = BrettVinylphenolReduction().derivatives(0.0, y, schema, params)

    dvg = d[schema.slice("vinylguaiacols")][0]
    deg = d[schema.slice("ethylguaiacols")][0]
    assert dvg < 0.0 and deg > 0.0
    moles_consumed = -dvg / M_VINYLGUAIACOL
    moles_produced = deg / M_ETHYLGUAIACOL
    assert moles_produced == pytest.approx(moles_consumed, rel=1e-12)
    assert d[schema.slice("vinylphenols")][0] == 0.0
    assert d[schema.slice("ethylphenols")][0] == 0.0

    # Both branches active simultaneously, independently (same catalyst/gate, distinct pools).
    y_vp_only = _state(schema, X_brett=0.3, vinylphenols=0.02)
    d_vp_only = BrettVinylphenolReduction().derivatives(0.0, y_vp_only, schema, params)

    y_both = _state(schema, X_brett=0.3, vinylphenols=0.02, vinylguaiacols=0.02)
    d_both = BrettVinylphenolReduction().derivatives(0.0, y_both, schema, params)
    assert d_both[schema.slice("vinylphenols")][0] == pytest.approx(
        d_vp_only[schema.slice("vinylphenols")][0], rel=1e-9
    )
    assert d_both[schema.slice("vinylguaiacols")][0] == pytest.approx(dvg, rel=1e-9)


# -- 4. guards ----------------------------------------------------------------


def test_guards_zero_without_catalyst_or_substrate(schema, params):
    """No X_brett ⇒ zero from both; no substrate ⇒ zero from the relevant Process."""
    no_brett = _state(schema, X_brett=0.0, hydroxycinnamics=0.1, vinylphenols=0.02)
    assert not np.any(BrettDecarboxylation().derivatives(0.0, no_brett, schema, params))
    assert not np.any(BrettVinylphenolReduction().derivatives(0.0, no_brett, schema, params))

    no_hc = _state(schema, X_brett=0.3, hydroxycinnamics=0.0)
    assert not np.any(BrettDecarboxylation().derivatives(0.0, no_hc, schema, params))
    no_vp = _state(schema, X_brett=0.3, vinylphenols=0.0)
    assert not np.any(BrettVinylphenolReduction().derivatives(0.0, no_vp, schema, params))


# -- 5. isolability: unpitched run is inert + keeps VALIDATED tier -------------


def test_unpitched_run_is_inert(schema):
    """No Brett pitch ⇒ phenol slots stay exactly 0 and keep their VALIDATED tier.

    The precursor may be present (dosed hydroxycinnamics), but with no Brett the decarboxylase
    Process is disabled at the compile seam, so nothing touches the phenol slots — byte-for-byte
    the validated core (prime directive #3).
    """
    compiled, traj = _run(hydroxycinnamic_gpl=0.1)  # precursor dosed, no Brett
    assert float(np.max(np.abs(traj.series("vinylphenols")))) == 0.0
    assert float(np.max(np.abs(traj.series("ethylphenols")))) == 0.0
    # tier isolability: an *enabled* zero Process would still drag these to speculative.
    assert compiled.process_set.tier_of("ethylphenols") is Tier.VALIDATED
    assert compiled.process_set.tier_of("vinylphenols") is Tier.VALIDATED


# -- 6. temperature optimum (warmer than O. oeni) -----------------------------


def test_temperature_optimum_warmer_than_mlf(params):
    """γ(T) peaks warm: 32 °C (Brett optimum) beats a cool 12 °C cellar; both beat the extremes."""
    schema = wine_schema()

    def gate_at(celsius: float) -> float:
        y = _state(schema, X_brett=0.3, T=celsius + 273.15)
        return brett_environmental_gate(y, schema, params, ph=0.0)

    assert gate_at(32.0) > gate_at(12.0) > 0.0  # warm-optimum, but active in a cool cellar
    assert gate_at(32.0) == pytest.approx(1.0, abs=1e-9)  # γ=1 at the optimum
    assert gate_at(5.0) == 0.0 and gate_at(50.0) == 0.0  # 0 outside the cardinals


# -- 7. tier ------------------------------------------------------------------


def test_processes_are_speculative():
    assert BrettDecarboxylation.tier is Tier.SPECULATIVE
    assert BrettVinylphenolReduction.tier is Tier.SPECULATIVE


# =============================================================================
# pt2 — BrettGrowth: X_brett becomes dynamic (decision D-40 pt2)
# =============================================================================
#
# Ranked headline-first. The payoff is the *autocatalytic* spoilage — a growing Brett population
# makes 4-EP accelerate over a barrel's months (``test_growth_accelerates_phenols``). The crown
# regression is ``test_growth_bounded_and_amino_acids_nonneg_under_bdf`` + its solver cross-check:
# an early build drove X_brett to 23 g/L and amino_acids to −4.5 g/L because BrettGrowth's ``E ≤ 0``
# guard had no smooth *shadow* (unlike the amino-acid Monod and the (1−X/K) brake), so the BDF
# Jacobian straddled an on/off step at E=0 during primary AF and the autocatalytic mode blew up.
# The ethanol Monod ``E/(K_E_brett+E)`` is that shadow; these tests pin the fix under the *default
# BDF solver* (RK45/LSODA never saw the bug — they build no Jacobian).


def _run_dry_post_af(*, method: str = "BDF", days: float = 200.0, **initial_extra):
    """A *finished* wine (no sugar, full-strength ethanol) pitched with Brett — Brett's real niche.

    ``brix=0``/``ethanol_gpl≈90`` puts the ethanol Monod at ≈1 from t0 (no near-zero cover), so a
    dosed low amino-acid pool is the tight resource — the case that discriminates "fixed the whole
    BDF blow-up class" from "fixed only the E-crossing instance". Solver ``method`` is explicit so
    the regression can cross-check BDF against the non-stiff solvers.
    """
    initial: dict[str, float] = {
        "brix": 0.0,
        "ethanol_gpl": 90.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.2,
        "brett_pitch_gpl": 0.05,
        "tartaric_gpl": 3.0,
        "malic_gpl": 2.0,
        "initial_ph": 3.5,
    }
    initial.update(initial_extra)
    sc = Scenario(
        name="brett-dry-post-af",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=25.0)],
        interventions=[],
        duration_days=days,
    )
    compiled = compile_scenario(sc, strict=True)
    traj = compiled.run(t_eval=np.linspace(0.0, days * 24.0, 500), method=method)
    return compiled, traj


# -- 8. HEADLINE: growth makes spoilage autocatalytic --------------------------


def test_growth_accelerates_phenols():
    """A growing Brett population (amino acids dosed) spoils *more* than a constant catalyst.

    Same Brett pitch and same precursor; the only difference is the amino-acid fuel that lets
    ``X_brett`` grow. The dynamic-population run must end with **more** ethylphenols *and* a higher
    ``X_brett`` than the constant-catalyst control — the "it gets worse the longer the barrel sits"
    dynamic a fixed catalyst (pt1) cannot produce.
    """
    _, dynamic = _run(
        days=120.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05, amino_acids_gpl=1.0
    )
    _, constant = _run(
        days=120.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05
    )  # no aa ⇒ no growth

    assert dynamic.series("X_brett")[-1] > 2.0 * constant.series("X_brett")[-1]  # population grew
    assert dynamic.series("ethylphenols")[-1] > constant.series("ethylphenols")[-1]  # more spoilage
    # The constant control's X_brett is genuinely constant (growth disabled without amino acids).
    xb_const = constant.series("X_brett")
    assert float(np.max(xb_const) - np.min(xb_const)) == pytest.approx(0.0, abs=1e-12)


# -- 9. CROWN REGRESSION: no runaway; amino_acids >= 0 under the default BDF ---


def test_growth_bounded_and_amino_acids_nonneg_under_bdf():
    """The X_brett=23 / amino_acids=−4.5 blow-up must not recur — under the DEFAULT BDF solver.

    A finished wine (ethanol Monod ≈1, no near-zero cover) with a small amino-acid pool: the
    carrying-capacity brake must bound ``X_brett`` and the pool must never go meaningfully negative.
    ``assert_nonnegative`` is the real gate (it is what caught the original bug); ``atol=1e-8`` sits
    an order of magnitude above the integrator's own ``atol=1e-9`` solver noise.
    """
    compiled, traj = _run_dry_post_af(method="BDF", amino_acids_gpl=0.05, hydroxycinnamic_gpl=0.1)
    k = compiled.param_values["brett_carrying_capacity"]
    xb = traj.series("X_brett")
    assert float(np.max(xb)) <= k * 1.01  # brake holds: bounded by the carrying capacity
    assert xb[-1] > xb[0]  # but it did grow (the process is active, not just inert)
    assert_nonnegative(traj, ("X_brett", "amino_acids", "E"), atol=1e-8)


def test_growth_bdf_matches_nonstiff_solvers():
    """BDF must agree with RK45 and LSODA — the bug was a BDF-only Jacobian artefact.

    The original runaway was invisible to every non-stiff solver (they build no Jacobian, so never
    straddled the E=0 step). Encoding "all three solvers agree" is the direct regression on the
    *class* of bug: if the ethanol-Monod shadow ever regresses, BDF diverges from RK45/LSODA here.
    """
    finals = {}
    for method in ("BDF", "RK45", "LSODA"):
        _, traj = _run_dry_post_af(method=method, amino_acids_gpl=0.05, hydroxycinnamic_gpl=0.1)
        finals[method] = traj.series("X_brett")[-1]
    assert finals["BDF"] == pytest.approx(finals["RK45"], rel=1e-3)
    assert finals["BDF"] == pytest.approx(finals["LSODA"], rel=1e-3)


# -- 10. conservation with growth active --------------------------------------


def test_growth_conserves_carbon_and_nitrogen():
    """total_carbon and total_nitrogen close while growth builds X_brett from amino acids + E."""
    compiled, traj = _run_dry_post_af(method="BDF", amino_acids_gpl=0.1, hydroxycinnamic_gpl=0.1)
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.param_values["biomass_C_fraction"]
    )
    nitrogen = total_nitrogen(
        compiled.schema, biomass_nitrogen_fraction=compiled.param_values["biomass_N_fraction"]
    )
    assert_conserved(traj, carbon, label="carbon")
    assert_conserved(traj, nitrogen, label="nitrogen")


# -- 11. stoichiometry, touches, and the ethanol carbon source ----------------


def test_growth_draws_ethanol_not_sugar_or_ammonium(schema, params):
    """Growth builds X_brett from amino acids + ETHANOL; touches neither sugar S nor ammonium N."""
    y = _state(schema, X_brett=0.1, amino_acids=0.5, E=90.0)
    d = BrettGrowth().derivatives(0.0, y, schema, params)

    assert d[schema.slice("X_brett")][0] > 0.0  # biomass grows
    assert d[schema.slice("amino_acids")][0] < 0.0  # nitrogen fuel consumed
    assert d[schema.slice("E")][0] < 0.0  # carbon shortfall drawn from ethanol
    assert float(d[schema.slice("S")].sum()) == 0.0  # NOT sugar (Brett's dry-wine niche)
    assert d[schema.slice("N")][0] == 0.0  # NOT ammonium (nitrogen-anchored on amino acids)
    assert set(BrettGrowth.touches) == {"X_brett", "amino_acids", "E"}


def test_growth_needs_ethanol_the_carbon_source(schema, params):
    """The ethanol Monod: no ethanol ⇒ ≈no growth; ample ethanol ⇒ growth — Brett feeds on ethanol.

    Also the smooth shadow of the ``E ≤ 0`` guard: growth eases to zero as E → 0 (rather than
    switching off at a step), which is what keeps the BDF Jacobian well-conditioned.
    """
    dry = BrettGrowth().derivatives(
        0.0, _state(schema, X_brett=0.1, amino_acids=0.5, E=0.0), schema, params
    )
    wet = BrettGrowth().derivatives(
        0.0, _state(schema, X_brett=0.1, amino_acids=0.5, E=90.0), schema, params
    )
    assert float(dry[schema.slice("X_brett")][0]) == 0.0  # no ethanol ⇒ no growth
    assert float(wet[schema.slice("X_brett")][0]) > 0.0  # ethanol present ⇒ growth

    # Smoothness of the shadow: a tiny ethanol level gives a tiny (but continuous, nonzero) rate —
    # no on/off step for the finite-difference Jacobian to straddle.
    trace = BrettGrowth().derivatives(
        0.0, _state(schema, X_brett=0.1, amino_acids=0.5, E=0.05), schema, params
    )
    r_trace = float(trace[schema.slice("X_brett")][0])
    r_wet = float(wet[schema.slice("X_brett")][0])
    assert 0.0 < r_trace < r_wet  # ramps up with ethanol availability, not a cliff


# -- 12. carrying-capacity brake ----------------------------------------------


def test_growth_brake_shuts_off_at_carrying_capacity(schema, params):
    """(1 − X/K): growth eases to 0 as X_brett → K and is exactly 0 at/above it (no runaway)."""
    k = params["brett_carrying_capacity"]
    below = BrettGrowth().derivatives(
        0.0, _state(schema, X_brett=0.5 * k, amino_acids=0.5, E=90.0), schema, params
    )
    at_cap = BrettGrowth().derivatives(
        0.0, _state(schema, X_brett=k, amino_acids=0.5, E=90.0), schema, params
    )
    over = BrettGrowth().derivatives(
        0.0, _state(schema, X_brett=1.5 * k, amino_acids=0.5, E=90.0), schema, params
    )
    assert float(below[schema.slice("X_brett")][0]) > 0.0  # growing below the cap
    assert float(at_cap[schema.slice("X_brett")][0]) == pytest.approx(0.0, abs=1e-15)  # shut at K
    assert float(over[schema.slice("X_brett")][0]) == 0.0  # clamped, never a biomass source


# -- 13. guards ---------------------------------------------------------------


def test_growth_guards_zero_without_catalyst_fuel_or_ethanol(schema, params):
    """Zero contribution when there is no Brett, no amino-acid fuel, or no ethanol."""
    no_brett = _state(schema, X_brett=0.0, amino_acids=0.5, E=90.0)
    no_aa = _state(schema, X_brett=0.1, amino_acids=0.0, E=90.0)
    no_e = _state(schema, X_brett=0.1, amino_acids=0.5, E=0.0)
    assert not np.any(BrettGrowth().derivatives(0.0, no_brett, schema, params))
    assert not np.any(BrettGrowth().derivatives(0.0, no_aa, schema, params))
    assert not np.any(BrettGrowth().derivatives(0.0, no_e, schema, params))


# -- 14. isolability + tier ---------------------------------------------------


def test_growth_disabled_without_amino_acids(schema):
    """Pitched Brett, no amino-acid dose ⇒ growth disabled at the compile seam ⇒ X_brett constant.

    The phenol pathway still runs (constant catalyst), but with no fuel the growth Process is
    disabled, so ``X_brett`` never changes and the ``amino_acids``/``E`` slots keep their tier — the
    stricter growth gate is isolable from the pt1 phenol gate (mirrors MLF growth vs conversion).
    """
    compiled, traj = _run(hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.3)  # pitched, no amino acids
    xb = traj.series("X_brett")
    assert float(np.max(xb) - np.min(xb)) == pytest.approx(0.0, abs=1e-12)
    assert "brett_growth" not in {p.name for p in compiled.process_set.active}


def test_growth_is_speculative():
    assert BrettGrowth.tier is Tier.SPECULATIVE


# =============================================================================
# pt3 — BrettDeath: the SO₂-driven Brett kill (decision D-40 pt3)
# =============================================================================
#
# Ranked headline-first. The payoff is ``test_so2_crashes_growing_brett_population``: a molecular
# SO₂ dose does not merely PAUSE a growing Brett population (its growth gate's g_SO₂ already does
# that) — it KILLS it, so ``X_brett`` DECLINES into ``X_brett_dead`` and the volatile-phenol accrual
# is curtailed vs the un-sulfited control. The rest mirror the ``MalolacticDeath`` RHS suite (D-39):
# death is exactly 0 without SO₂, is a carbon/nitrogen-neutral X_brett→X_brett_dead transfer, is
# monotone in the SO₂ dose, carries an ARRHENIUS temperature factor (not the cardinal γ(T), so cold
# preserves rather than kills), closes both ledgers at the integration level, and is speculative.


def _death_state(
    schema: StateSchema,
    params,
    *,
    so2_mgl: float = 0.0,
    temp_k: float = 293.15,
    x_brett: float = 0.2,
    x_brett_dead: float = 0.0,
) -> FloatArray:
    """A pitched wine state for exercising BrettDeath at the RHS level (the MLF `_death_state`)."""
    overrides: dict[str, float | list[float]] = {
        "X_brett": x_brett,
        "X_brett_dead": x_brett_dead,
        "T": temp_k,
    }
    if so2_mgl > 0.0:
        overrides["so2_total"] = mgl_to_gpl(so2_mgl)
    return _state(schema, **overrides)


def test_death_is_exactly_zero_without_so2(schema, params):
    """Death is driven by molecular SO₂ ALONE, so an unsulfited pitched run never kills Brett.

    The v1 tradeoff, enforced at the RHS: without SO₂ the population is inert (Brett persists in the
    barrel) — no ethanol/pH decay term, unlike a naive copy of an ethanol-tolerance organism.
    """
    y = _death_state(schema, params, so2_mgl=0.0, x_brett=0.2)
    d = BrettDeath().derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("X_brett")][0]) == 0.0
    assert float(d[schema.slice("X_brett_dead")][0]) == 0.0


def test_so2_drives_death_as_a_neutral_transfer(schema, params):
    """With SO₂ dosed, viable X_brett leaves and the SAME mass enters X_brett_dead (the D-13 idiom).

    ``d[X_brett] = −d[X_brett_dead]`` exactly, so — both weighted at the biomass fractions since
    pt2 — the move is carbon- and nitrogen-neutral by construction.
    """
    y = _death_state(schema, params, so2_mgl=80.0, x_brett=0.2)
    d = BrettDeath().derivatives(0.0, y, schema, params)
    dx = float(d[schema.slice("X_brett")][0])
    dxd = float(d[schema.slice("X_brett_dead")][0])
    assert dx < 0.0 and dxd > 0.0  # Brett dies
    assert dxd == pytest.approx(-dx)  # mass-conserving transfer (neutral in both ledgers)


def test_death_touches_only_the_x_brett_pools(schema, params):
    y = _death_state(schema, params, so2_mgl=80.0, x_brett=0.2)
    d = BrettDeath().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"X_brett", "X_brett_dead"}
    assert set(BrettDeath.touches) == {"X_brett", "X_brett_dead"}


def test_more_so2_kills_faster(schema, params):
    """Monotone in the antimicrobial dose: more SO₂ ⇒ higher molecular SO₂ ⇒ larger 1 − g_SO₂."""
    rate_lo = -float(
        BrettDeath().derivatives(0.0, _death_state(schema, params, so2_mgl=20.0), schema, params)[
            schema.slice("X_brett")
        ][0]
    )
    rate_hi = -float(
        BrettDeath().derivatives(0.0, _death_state(schema, params, so2_mgl=60.0), schema, params)[
            schema.slice("X_brett")
        ][0]
    )
    assert 0.0 < rate_lo < rate_hi


def test_cold_preserves_brett_via_arrhenius_not_gamma(schema, params):
    """Death carries its OWN Arrhenius factor, not the cardinal γ(T) — so cold preserves, not kills.

    Below ``T_min_brett`` the cardinal γ(T) = 0; if death reused it, cold would spuriously HALT the
    SO₂ kill. Arrhenius instead merely SLOWS it: a sulfited Brett culture in a cold cellar still
    dies (just slower), and warmer ⇒ faster — the "cold *preserves* a Brett infection" direction
    γ(T) cannot supply. So Brett survives cold storage; it is cleared by SO₂, not by chilling.
    """
    t_below_min = 278.15  # 5 °C, below T_min_brett (10 °C) ⇒ cardinal γ(T) = 0
    assert (
        cardinal_temperature_factor(
            t_below_min, params["T_min_brett"], params["T_opt_brett"], params["T_max_brett"]
        )
        == 0.0
    )
    rate_cold = -float(
        BrettDeath().derivatives(
            0.0, _death_state(schema, params, so2_mgl=80.0, temp_k=t_below_min), schema, params
        )[schema.slice("X_brett")][0]
    )
    rate_warm = -float(
        BrettDeath().derivatives(
            0.0, _death_state(schema, params, so2_mgl=80.0, temp_k=298.15), schema, params
        )[schema.slice("X_brett")][0]
    )
    assert rate_cold > 0.0  # still dying below T_min — proves Arrhenius, not γ(T)
    assert rate_warm > rate_cold  # warm accelerates the kill (Arrhenius direction)


# -- HEADLINE: SO₂ crashes a growing Brett population (not just arrests it) ----


def test_so2_crashes_growing_brett_population():
    """A molecular-SO₂ dose KILLS a growing Brett population — it does not merely pause it.

    Brett is dosed with amino acids so ``X_brett`` grows autocatalytically (pt2). A mid-run SO₂
    addition then both arrests growth (g_SO₂ in the growth gate) AND kills it off (BrettDeath).
    The unambiguous *death* signal — distinct from growth-arrest, which the g_SO₂ gate alone would
    give — is that ``X_brett_dead`` **accumulates** and ``X_brett`` **falls below its value at the
    dose**; the winemaking payoff is that ethylphenols end **lower** than the un-sulfited control.
    """
    dose_day = 40.0
    so2 = [Intervention(day=dose_day, action="add_so2", params={"so2_mgl": 80.0})]
    _, sulfited = _run(
        days=140.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05, amino_acids_gpl=1.0,
        interventions=so2,
    )  # fmt: skip
    _, control = _run(
        days=140.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05, amino_acids_gpl=1.0
    )  # no SO₂: Brett keeps growing and spoiling

    xb = sulfited.series("X_brett")
    xbd = sulfited.series("X_brett_dead")
    xb_at_dose = float(np.interp(dose_day * 24.0, sulfited.t, xb))

    # Death, not mere arrest: the dead pool fills and the viable pool declines past the dose point.
    assert xbd[-1] > 0.0
    assert xb[-1] < xb_at_dose
    # The control (no SO₂) keeps growing — the population was genuinely rising when SO₂ hit it.
    assert control.series("X_brett")[-1] > xb_at_dose
    # Winemaking payoff: killing Brett curtails the volatile-phenol spoilage.
    assert sulfited.series("ethylphenols")[-1] < control.series("ethylphenols")[-1]


def test_death_run_conserves_carbon_and_nitrogen():
    """total_carbon and total_nitrogen close with death ACTIVE (SO₂ dosed mid-run, X_brett → dead).

    The X_brett → X_brett_dead transfer is carbon/nitrogen-neutral (both weighted since pt2), SO₂
    carries neither element, and growth draws from the weighted amino-acid/ethanol pools — so both
    run-wide ledgers close across the sulfite jump (final == initial + Σ external flows).
    """
    so2 = [Intervention(day=40.0, action="add_so2", params={"so2_mgl": 80.0})]
    compiled, traj = _run(
        days=120.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05, amino_acids_gpl=1.0,
        interventions=so2,
    )  # fmt: skip
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.param_values["biomass_C_fraction"]
    )
    nitrogen = total_nitrogen(
        compiled.schema, biomass_nitrogen_fraction=compiled.param_values["biomass_N_fraction"]
    )
    assert_conserved(traj, carbon, label="carbon (Brett death on)")
    assert_conserved(traj, nitrogen, label="nitrogen (Brett death on)")
    assert_nonnegative(traj, ("X_brett", "X_brett_dead"))


def test_death_is_speculative():
    assert BrettDeath.tier is Tier.SPECULATIVE


# =============================================================================
# D-58 — BrettEthanolToxicity: the sourced ethanol-toxicity kill (no SO₂ needed)
# =============================================================================
#
# D-52 declined a generic age-based "BrettSenescence"; D-58's research re-confirmed that but
# surfaced a real, DIFFERENT, sourced mechanism (Barata et al. 2008): Brett grows normally to ~14%
# v/v ethanol and is fully arrested by ~14.5-15%, no SO₂ required. Ranked headline-first. The payoff
# is ``test_headline_high_ethanol_crashes_unsulfited_brett_population``: an UNSULFITED, high-ethanol
# (~13% ABV) wine crashes a growing Brett population purely on ethanol, contrasting with a normal-
# strength (~11% ABV) control that keeps growing — the "no SO₂ needed" headline. The rest pin the
# exact-zero guard at/below the onset (ordinary wine strength is unaffected), the neutral transfer,
# monotonicity between onset and ceiling, the reused Arrhenius temperature shape, the growth wall's
# reconciliation with BrettGrowth's existing ethanol-as-carbon-source Monod, the survival-factor
# helper's own boundary values, touches/reads, and the speculative tier.


def _ethanol_toxicity_state(
    schema: StateSchema,
    *,
    e_gpl: float,
    temp_k: float = 293.15,
    x_brett: float = 0.2,
    x_brett_dead: float = 0.0,
) -> FloatArray:
    """A pitched wine state for exercising BrettEthanolToxicity at the RHS level."""
    return _state(schema, X_brett=x_brett, X_brett_dead=x_brett_dead, E=e_gpl, T=temp_k)


def test_ethanol_toxicity_is_exactly_zero_at_or_below_onset(schema, params):
    """Ordinary wine strength (E <= onset, ~14% v/v) sees NO contribution — a genuine zero, not a
    small one, so a typical ~11-13% ABV wine is completely unaffected by this Process."""
    onset = params["brett_ethanol_toxicity_onset"]
    for e in (90.0, onset - 1.0, onset):
        y = _ethanol_toxicity_state(schema, e_gpl=e)
        d = BrettEthanolToxicity().derivatives(0.0, y, schema, params)
        assert float(d[schema.slice("X_brett")][0]) == 0.0
        assert float(d[schema.slice("X_brett_dead")][0]) == 0.0


def test_ethanol_toxicity_kills_above_onset_as_a_neutral_transfer(schema, params):
    """Above the onset, viable X_brett leaves and the SAME mass enters X_brett_dead (the BrettDeath
    idiom) — carbon/nitrogen-neutral by construction, no SO₂ dosed."""
    onset = params["brett_ethanol_toxicity_onset"]
    y = _ethanol_toxicity_state(schema, e_gpl=onset + 5.0)
    d = BrettEthanolToxicity().derivatives(0.0, y, schema, params)
    dx = float(d[schema.slice("X_brett")][0])
    dxd = float(d[schema.slice("X_brett_dead")][0])
    assert dx < 0.0 and dxd > 0.0
    assert dxd == pytest.approx(-dx)


def test_ethanol_toxicity_touches_only_the_x_brett_pools(schema, params):
    onset = params["brett_ethanol_toxicity_onset"]
    y = _ethanol_toxicity_state(schema, e_gpl=onset + 5.0)
    d = BrettEthanolToxicity().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"X_brett", "X_brett_dead"}
    assert set(BrettEthanolToxicity.touches) == {"X_brett", "X_brett_dead"}


def test_higher_ethanol_kills_faster_up_to_the_ceiling(schema, params):
    """Monotone toxicity between the onset and the ceiling: more ethanol => more death."""
    onset = params["brett_ethanol_toxicity_onset"]
    ceiling = params["brett_ethanol_toxicity_ceiling"]
    mid = (onset + ceiling) / 2.0
    rate_lo = -float(
        BrettEthanolToxicity().derivatives(
            0.0, _ethanol_toxicity_state(schema, e_gpl=mid), schema, params
        )[schema.slice("X_brett")][0]
    )
    rate_hi = -float(
        BrettEthanolToxicity().derivatives(
            0.0, _ethanol_toxicity_state(schema, e_gpl=ceiling), schema, params
        )[schema.slice("X_brett")][0]
    )
    assert 0.0 < rate_lo < rate_hi
    # At (or above) the ceiling, toxicity is at its full k_death_brett-scaled rate (survival = 0).
    assert rate_hi == pytest.approx(params["k_death_brett"] * 0.2, rel=1e-9)


def test_ethanol_toxicity_warmer_kills_faster(schema, params):
    """Reuses BrettDeath's Arrhenius factor (D-58): warm accelerates, mirroring the SO₂ kill."""
    ceiling = params["brett_ethanol_toxicity_ceiling"]
    rate_cold = -float(
        BrettEthanolToxicity().derivatives(
            0.0, _ethanol_toxicity_state(schema, e_gpl=ceiling, temp_k=283.15), schema, params
        )[schema.slice("X_brett")][0]
    )
    rate_warm = -float(
        BrettEthanolToxicity().derivatives(
            0.0, _ethanol_toxicity_state(schema, e_gpl=ceiling, temp_k=303.15), schema, params
        )[schema.slice("X_brett")][0]
    )
    assert 0.0 < rate_cold < rate_warm


def test_ethanol_toxicity_needs_no_catalyst_or_ph_solve(schema, params):
    """No X_brett ⇒ zero contribution; the Process never solves pH (no SO₂ term at all)."""
    y = _ethanol_toxicity_state(schema, e_gpl=130.0, x_brett=0.0)
    d = BrettEthanolToxicity().derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("X_brett")][0]) == 0.0
    assert float(d[schema.slice("X_brett_dead")][0]) == 0.0


def test_ethanol_toxicity_is_speculative():
    assert BrettEthanolToxicity.tier is Tier.SPECULATIVE


def test_survival_factor_boundary_values(params):
    """Direct unit tests of brett_ethanol_survival_factor: 1 at/below onset, 0 at/above ceiling."""
    onset = params["brett_ethanol_toxicity_onset"]
    ceiling = params["brett_ethanol_toxicity_ceiling"]
    assert brett_ethanol_survival_factor(0.0, params) == 1.0
    assert brett_ethanol_survival_factor(onset, params) == 1.0
    assert brett_ethanol_survival_factor(ceiling, params) == 0.0
    assert brett_ethanol_survival_factor(ceiling + 50.0, params) == 0.0  # clamped, no overshoot
    mid_value = brett_ethanol_survival_factor((onset + ceiling) / 2.0, params)
    assert 0.0 < mid_value < 1.0


def test_growth_wall_leaves_normal_wine_strength_unaffected(schema, params):
    """BrettGrowth at ordinary wine strength (E <= onset) is UNCHANGED by the D-58 wall — the
    growth-rate direction test this codebase's own `test_pitch_brett_post_af_at_high_ethanol`
    integration test depends on (no wall at full-strength wine ethanol, ~106 g/L < onset)."""
    onset = params["brett_ethanol_toxicity_onset"]
    y_below = _state(schema, X_brett=0.1, amino_acids=1.0, E=onset - 20.0)
    y_at_onset = _state(schema, X_brett=0.1, amino_acids=1.0, E=onset)
    d_below = BrettGrowth().derivatives(0.0, y_below, schema, params)
    d_at_onset = BrettGrowth().derivatives(0.0, y_at_onset, schema, params)
    assert float(d_below[schema.slice("X_brett")][0]) > 0.0
    assert float(d_at_onset[schema.slice("X_brett")][0]) > 0.0


def test_growth_wall_arrests_growth_near_the_ceiling(schema, params):
    """Above the onset, BrettGrowth eases toward 0 as E approaches the ceiling (decision D-58) —
    the reconciliation of ethanol-as-carbon-source (low E) and ethanol-as-toxin (high E) on the
    SAME state variable."""
    ceiling = params["brett_ethanol_toxicity_ceiling"]
    onset = params["brett_ethanol_toxicity_onset"]
    y_mid = _state(schema, X_brett=0.1, amino_acids=1.0, E=(onset + ceiling) / 2.0)
    y_ceiling = _state(schema, X_brett=0.1, amino_acids=1.0, E=ceiling)
    rate_mid = float(
        BrettGrowth().derivatives(0.0, y_mid, schema, params)[schema.slice("X_brett")][0]
    )
    rate_ceiling = float(
        BrettGrowth().derivatives(0.0, y_ceiling, schema, params)[schema.slice("X_brett")][0]
    )
    assert rate_mid > 0.0
    assert rate_ceiling == 0.0  # survival factor is exactly 0 at the ceiling: growth fully arrested


# -- HEADLINE: high ethanol crashes an unsulfited, growing Brett population ---


def test_headline_high_ethanol_crashes_unsulfited_brett_population():
    """A HIGH-ethanol wine (~13% ABV, above brett_ethanol_toxicity_onset) crashes a growing,
    UNSULFITED Brett population — the D-58 headline: no SO₂ needed, unlike BrettDeath.

    Contrasts with a normal-strength (~11% ABV, brix 22, below onset) control that keeps growing —
    the same population, the same amino-acid dose, differing only in must sugar/ethanol.
    """
    _, high_ethanol = _run(
        days=140.0, brix=26.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05, amino_acids_gpl=1.0
    )
    _, normal = _run(
        days=140.0, brix=22.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05, amino_acids_gpl=1.0
    )

    xb_high = high_ethanol.series("X_brett")
    xbd_high = high_ethanol.series("X_brett_dead")

    assert xbd_high[-1] > 0.0  # the dead pool fills, purely from ethanol toxicity (no SO₂ dosed)
    assert xb_high[-1] < float(np.max(xb_high))  # viable population declines from its peak
    # The normal-strength control (below onset) keeps growing rather than crashing.
    assert normal.series("X_brett")[-1] >= float(np.max(normal.series("X_brett"))) * 0.999


def test_ethanol_toxicity_run_conserves_carbon_and_nitrogen():
    """total_carbon/total_nitrogen close with the ethanol-toxicity kill ACTIVE (no SO₂ dosed)."""
    compiled, traj = _run(
        days=140.0, brix=26.0, hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.05, amino_acids_gpl=1.0
    )
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.param_values["biomass_C_fraction"]
    )
    nitrogen = total_nitrogen(
        compiled.schema, biomass_nitrogen_fraction=compiled.param_values["biomass_N_fraction"]
    )
    assert_conserved(traj, carbon, label="carbon (Brett ethanol toxicity on)")
    assert_conserved(traj, nitrogen, label="nitrogen (Brett ethanol toxicity on)")
    assert_nonnegative(traj, ("X_brett", "X_brett_dead"))


# =============================================================================
# pt4 — YeastPOFDecarboxylation: the POF+ yeast opt-in + emergent reservoir (decision D-40 pt4)
# =============================================================================
#
# Ranked headline-first. The PRIMARY payoff is the STRANDING test
# (``test_pof_strands_vinylphenols_without_brett``): a POF+ strain (opted in, NO Brett)
# decarboxylates must hydroxycinnamics into ``vinylphenols`` DURING AF, but - lacking the
# reductase - cannot take them to ethylphenols, so the reservoir *strands* (``ethylphenols`` stays
# exactly 0 and VALIDATED). It is the timing-independent control-difference, the parallel of the pt1
# headline. The SECONDARY payoff is the emergent HEAD START (``test_pof_gives_brett_a_head_start``):
# a Brett contamination arriving after a POF+ AF finds a pre-filled vinylphenol reservoir and
# reaches a given 4-EP level SOONER than into a POF-negative wine. That is deliberately an
# EARLY-TIME (kinetic) claim, NOT an endpoint one: with the same total hydroxycinnamics in both
# arms, conservation forces the *asymptotic* ethylphenols EQUAL (all hc -> ep eventually), so the
# difference is only that POF+ has vinylphenol pre-made and ready to reduce. The rest pin the
# decarboxylase stoichiometry / carbon closure / ``touches``, the flux-coupled guards, the
# POF-independent-of-Brett compile gate (a POF+ ferment need not have Brett), the default (POF-)
# isolability, and the ``speculative`` tier.


# -- 15. HEADLINE: POF+ fills the reservoir but strands it (no Brett reductase) ---


def test_pof_strands_vinylphenols_without_brett():
    """POF+ yeast makes ``vinylphenols`` during AF; with no Brett they **strand** (ep stays 0).

    The emergent yeast/Brett coupling the 3-pool design was chosen for (the α-acetolactate-reservoir
    parallel, D-26/D-31). A POF+ strain carries the decarboxylase but not the reductase, so it fills
    the shared reservoir it cannot drain: ``vinylphenols`` rise and remain, while ``ethylphenols``
    stays **exactly 0** — nothing reduces vinylphenol without Brett. Tier is the honest consequence:
    ``vinylphenols`` reports speculative (the enabled POF Process touches it) while ``ethylphenols``
    stays **VALIDATED at 0** (no enabled Process touches it — the reductase is Brett's, and Brett is
    absent).
    """
    compiled, traj = _run(hydroxycinnamic_gpl=0.1, pof_positive=1.0)  # POF+, no Brett

    vp = traj.series("vinylphenols")
    ep = traj.series("ethylphenols")
    hc = traj.series("hydroxycinnamics")

    assert float(np.max(vp)) > 1e-3  # a real reservoir accumulates during AF
    assert vp[-1] > 0.5 * float(np.max(vp))  # and it STRANDS (no reductase to drain it)
    assert float(np.max(np.abs(ep))) == 0.0  # no Brett ⇒ no ethylphenols at all
    assert hc[-1] < 0.9 * hc[0]  # the precursor is genuinely consumed into vinylphenol

    # Tier isolability of the stranding: vinylphenol is touched (speculative); ethylphenol is not.
    assert compiled.process_set.tier_of("vinylphenols") is Tier.SPECULATIVE
    assert compiled.process_set.tier_of("ethylphenols") is Tier.VALIDATED


def test_pof_strands_vinylguaiacols_too_without_brett():
    """The D-55 ferulic branch strands identically to the p-coumaric branch (decision D-55).

    Same mechanism as ``test_pof_strands_vinylphenols_without_brett``: POF+ yeast decarboxylates
    ferulic_acid to vinylguaiacols during AF but has no reductase, so with no Brett present
    ethylguaiacols stays exactly 0 too — the reduction step is entirely Brett's, for both branches.
    """
    compiled, traj = _run(ferulic_acid_gpl=0.06, pof_positive=1.0)  # POF+, no Brett, no hc dosed

    vg = traj.series("vinylguaiacols")
    eg = traj.series("ethylguaiacols")
    fa = traj.series("ferulic_acid")

    assert float(np.max(vg)) > 1e-3  # a real reservoir accumulates during AF
    assert vg[-1] > 0.5 * float(np.max(vg))  # and it STRANDS (no reductase to drain it)
    assert float(np.max(np.abs(eg))) == 0.0  # no Brett ⇒ no ethylguaiacols at all
    assert fa[-1] < 0.9 * fa[0]  # the precursor is genuinely consumed into vinylguaiacol

    assert compiled.process_set.tier_of("vinylguaiacols") is Tier.SPECULATIVE
    assert compiled.process_set.tier_of("ethylguaiacols") is Tier.VALIDATED


# -- 16. SECONDARY: a POF+ AF gives a later Brett a head start (early-time claim) --


def _run_brett_pitched_post_af(*, pof: bool, days: float = 120.0, pitch_day: float = 10.0):
    """Ferment (POF+ or POF-), then pitch Brett post-AF - the head-start comparison's two arms."""
    interventions = [Intervention(day=pitch_day, action="pitch_brett", params={"pitch_gpl": 0.3})]
    if pof:
        return _run(
            days=days, interventions=interventions, hydroxycinnamic_gpl=0.1, pof_positive=1.0
        )
    return _run(days=days, interventions=interventions, hydroxycinnamic_gpl=0.1)


def test_pof_gives_brett_a_head_start():
    """Brett into a POF+ (pre-filled reservoir) wine reaches 4-EP SOONER than into a POF− wine.

    The emergent coupling's winemaking meaning: a POF+ primary ferment hands a subsequent Brett
    contamination a running start, because the vinylphenol is already made and only needs reducing.

    This is asserted as an **early-time / time-to-threshold** claim, NOT an endpoint one. With the
    same total hydroxycinnamics in both arms, conservation forces the asymptotic ethylphenols EQUAL
    (all hc → ep eventually); the POF+ advantage is purely kinetic (vinylphenol pre-made). Asserting
    higher *final* ep would be wrong — the arms converge.
    """
    pitch_day = 10.0
    _, pos = _run_brett_pitched_post_af(pof=True, pitch_day=pitch_day)
    _, neg = _run_brett_pitched_post_af(pof=False, pitch_day=pitch_day)

    t = pos.t
    ep_pos = pos.series("ethylphenols")
    ep_neg = neg.series("ethylphenols")

    # Early-time: shortly after the pitch the POF+ arm is far ahead (vinylphenol ready to reduce).
    probe_h = (pitch_day + 3.0) * 24.0
    ep_pos_early = float(np.interp(probe_h, t, ep_pos))
    ep_neg_early = float(np.interp(probe_h, t, ep_neg))
    assert ep_pos_early > 1e-3  # POF+ is already producing 4-EP from the pre-made reservoir
    assert ep_pos_early > 5.0 * ep_neg_early  # a clear head start over the from-scratch POF− arm

    # Time-to-threshold (half the POF− endpoint): POF+ crosses it comfortably sooner.
    thr = 0.5 * ep_neg[-1]
    t_pos = t[np.argmax(ep_pos >= thr)]
    t_neg = t[np.argmax(ep_neg >= thr)]
    assert t_pos < t_neg - 5.0 * 24.0  # at least ~5 days sooner


# -- 17. per-Process stoichiometry + touches ----------------------------------


def test_pof_decarboxylation_stoichiometry_and_touches(schema, params):
    """Same reaction as Brett's decarboxylase: hc down, vinylphenols + CO₂ up, carbon flux sums 0.

    Catalyst is viable yeast via the fermentative flux (``X``/``S`` present in ``_state``), not
    ``X_brett`` — so the RHS is nonzero here with no Brett dosed at all.
    """
    y = _state(schema, hydroxycinnamics=0.1)  # X=0.1, S=[200] ⇒ fermentative flux > 0, no X_brett
    d = YeastPOFDecarboxylation().derivatives(0.0, y, schema, params)

    assert d[schema.slice("hydroxycinnamics")][0] < 0.0
    assert d[schema.slice("vinylphenols")][0] > 0.0
    assert d[schema.slice("CO2")][0] > 0.0
    assert d[schema.slice("ferulic_acid")][0] == 0.0  # undosed branch is exactly inert
    assert d[schema.slice("vinylguaiacols")][0] == 0.0
    # Carbon flux across the three touched slots must cancel (9 C = 8 C + 1 C).
    c_flux = (
        d[schema.slice("hydroxycinnamics")][0] * carbon_mass_fraction("p_coumaric_acid")
        + d[schema.slice("vinylphenols")][0] * carbon_mass_fraction("vinylphenol")
        + d[schema.slice("CO2")][0] * carbon_mass_fraction("CO2")
    )
    assert c_flux == pytest.approx(0.0, abs=1e-12)
    assert set(YeastPOFDecarboxylation.touches) == {
        "hydroxycinnamics",
        "vinylphenols",
        "ferulic_acid",
        "vinylguaiacols",
        "CO2",
    }


def test_pof_decarboxylation_ferulic_branch_stoichiometry(schema, params):
    """The D-55 ferulic branch for POF+ yeast: ferulic_acid down, vinylguaiacols + CO₂ up.

    Same carbon-closing form (10 = 9 + 1) as the Brett ferulic branch, but flux-coupled (no
    X_brett needed) — the POF+ yeast decarboxylase analogue of
    ``test_decarboxylation_ferulic_branch_stoichiometry``.
    """
    y = _state(schema, ferulic_acid=0.05)  # no hydroxycinnamics dosed
    d = YeastPOFDecarboxylation().derivatives(0.0, y, schema, params)

    assert d[schema.slice("ferulic_acid")][0] < 0.0
    assert d[schema.slice("vinylguaiacols")][0] > 0.0
    assert d[schema.slice("CO2")][0] > 0.0
    assert d[schema.slice("hydroxycinnamics")][0] == 0.0
    assert d[schema.slice("vinylphenols")][0] == 0.0
    c_flux = (
        d[schema.slice("ferulic_acid")][0] * carbon_mass_fraction("ferulic_acid")
        + d[schema.slice("vinylguaiacols")][0] * carbon_mass_fraction("vinylguaiacol")
        + d[schema.slice("CO2")][0] * carbon_mass_fraction("CO2")
    )
    assert c_flux == pytest.approx(0.0, abs=1e-12)


# -- 18. guards: flux-coupled, no X_brett needed ------------------------------


def test_pof_guards_zero_without_precursor_or_flux(schema, params):
    """Zero without precursor, and zero post-AF (no fermentative flux: S=0 or dead yeast)."""
    no_hc = _state(schema, hydroxycinnamics=0.0)
    assert not np.any(YeastPOFDecarboxylation().derivatives(0.0, no_hc, schema, params))

    no_sugar = _state(schema, hydroxycinnamics=0.1, S=[0.0])  # dryness ⇒ flux 0 ⇒ POF stops
    assert not np.any(YeastPOFDecarboxylation().derivatives(0.0, no_sugar, schema, params))

    no_yeast = _state(schema, hydroxycinnamics=0.1, X=0.0)  # crashed yeast ⇒ flux 0
    assert not np.any(YeastPOFDecarboxylation().derivatives(0.0, no_yeast, schema, params))


# -- 19. carbon closes with POF active (alone, and alongside Brett) -----------


def test_pof_carbon_closes(schema):
    """total_carbon closes with POF decarboxylation active — alone and composed with Brett."""
    compiled, traj = _run(hydroxycinnamic_gpl=0.12, pof_positive=1.0)  # POF+ only
    fn = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.param_values["biomass_C_fraction"]
    )
    assert_conserved(traj, fn, label="carbon (POF+)")
    assert_nonnegative(traj, ("hydroxycinnamics", "vinylphenols", "ethylphenols"))

    # Both decarboxylases active: POF+ yeast and pitched Brett draw the same hydroxycinnamic pool.
    compiled2, traj2 = _run(hydroxycinnamic_gpl=0.12, pof_positive=1.0, brett_pitch_gpl=0.3)
    fn2 = total_carbon(
        compiled2.schema, biomass_carbon_fraction=compiled2.param_values["biomass_C_fraction"]
    )
    assert_conserved(traj2, fn2, label="carbon (POF+ and Brett)")


# -- 20. isolability: POF- default is inert + keeps VALIDATED tier ------------


def test_pof_negative_default_is_inert(schema):
    """No POF opt-in ⇒ phenol slots stay exactly 0 and keep VALIDATED — byte-for-byte the core.

    A POF-negative wine (the default) must make no vinylphenol even with hydroxycinnamics dosed: the
    POF Process is disabled at the compile seam, so nothing touches the phenol slots. This is a
    SEPARATE gate from the Brett pitch — here Brett is also absent, but the POF gate alone suffices.
    """
    compiled, traj = _run(hydroxycinnamic_gpl=0.1)  # precursor dosed, POF off, no Brett
    assert float(np.max(np.abs(traj.series("vinylphenols")))) == 0.0
    assert float(np.max(np.abs(traj.series("ethylphenols")))) == 0.0
    assert compiled.process_set.tier_of("vinylphenols") is Tier.VALIDATED
    assert "yeast_pof_decarboxylation" not in {p.name for p in compiled.process_set.active}


def test_pof_gate_is_independent_of_the_brett_pitch(schema):
    """POF+ enables its decarboxylase with NO Brett pitch; Brett-only leaves POF disabled.

    The two gates are orthogonal (decision D-40 pt4): ``pof_positive`` enables the yeast
    decarboxylase, ``brett_pitch_gpl`` enables the Brett Processes — neither implies the other.
    """
    pof_only, _ = _run(hydroxycinnamic_gpl=0.1, pof_positive=1.0)
    brett_only, _ = _run(hydroxycinnamic_gpl=0.1, brett_pitch_gpl=0.3)

    pof_active = {p.name for p in pof_only.process_set.active}
    brett_active = {p.name for p in brett_only.process_set.active}
    assert "yeast_pof_decarboxylation" in pof_active
    assert "brett_decarboxylation" not in pof_active  # POF+ does not enable Brett
    assert "yeast_pof_decarboxylation" not in brett_active  # a Brett pitch does not make yeast POF+
    assert "brett_decarboxylation" in brett_active


def test_pof_is_speculative():
    assert YeastPOFDecarboxylation.tier is Tier.SPECULATIVE


# -- 21. E_a_pof: intrinsic rate rises with warmth, but NET conversion falls (decision D-54) --


def test_pof_own_rate_rises_with_warmth(schema, params):
    """The decarboxylase's OWN rate (fixed flux/precursor) is faster warm than cold — E_a_pof > 0.

    Isolates the raw Arrhenius direction from the flux-coupling effect: at the SAME flux and
    precursor concentration, only T differs, so this pins that ``E_a_pof`` itself is a genuine
    positive (enzyme-accelerates-with-warmth) term — not the net finished-wine direction, which the
    flux-coupling flips (see ``test_pof_net_conversion_falls_with_warmer_fermentation``).
    """
    y_cold = _state(schema, hydroxycinnamics=0.1, T=285.15)  # 12 C
    y_warm = _state(schema, hydroxycinnamics=0.1, T=301.15)  # 28 C
    r_cold = float(
        YeastPOFDecarboxylation().derivatives(0.0, y_cold, schema, params)[schema.slice("CO2")][0]
    )
    r_warm = float(
        YeastPOFDecarboxylation().derivatives(0.0, y_warm, schema, params)[schema.slice("CO2")][0]
    )
    assert r_warm > r_cold > 0.0


def _pof_run_at_temperature(
    *, celsius: float, days: float = 60.0, hydroxycinnamic_gpl: float = 0.1
):
    """A POF+, no-Brett wine ferment held isothermal at ``celsius`` (net-conversion-vs-T probe)."""
    initial: dict[str, float] = {
        "brix": 22.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.2,
        "tartaric_gpl": 3.0,
        "malic_gpl": 2.0,
        "initial_ph": 3.5,
        "hydroxycinnamic_gpl": hydroxycinnamic_gpl,
        "pof_positive": 1.0,
    }
    sc = Scenario(
        name="wine-pof-temperature",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
        interventions=[],
        duration_days=days,
    )
    compiled = compile_scenario(sc, strict=True)
    traj = compiled.run(t_eval=np.linspace(0.0, days * 24.0, 400))
    return compiled, traj


def test_pof_net_conversion_falls_with_warmer_fermentation():
    """The NET (finished-must) vinylphenol conversion is HIGHER from a cooler ferment (D-54).

    POF decarboxylation is flux-coupled, so a warmer ferment finishes faster and gives the
    decarboxylase a shorter time window to act — and that shrinking window outweighs the enzyme's
    own (smaller) Arrhenius acceleration (``E_a_pof < E_a_uptake``, by design). Net result: a cool
    ferment strands MORE vinylphenol than a warm one, matching the sourced real-world direction
    (cooler wheat-beer fermentation retains more clove/4-vinylguaiacol character — the same Pad1/
    Fdc1 enzyme). Both runs reach dryness comfortably within the shared 60-day window, so this
    compares two frozen (post-dryness) totals, not two different-length integration windows.
    """
    _, cool = _pof_run_at_temperature(celsius=12.0)
    _, warm = _pof_run_at_temperature(celsius=28.0)

    vp_cool = float(cool.series("vinylphenols")[-1])
    vp_warm = float(warm.series("vinylphenols")[-1])
    assert vp_cool > vp_warm > 0.0
