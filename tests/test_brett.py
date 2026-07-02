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
    M_ETHYLPHENOL,
    M_VINYLPHENOL,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.brett import (
    BrettDecarboxylation,
    BrettVinylphenolReduction,
    brett_environmental_gate,
)
from fermentation.core.media import wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.scenario.schema import Intervention
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon


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


# -- 3. per-Process stoichiometry + touches -----------------------------------


def test_decarboxylation_stoichiometry_and_touches(schema, params):
    """hydroxycinnamics down; vinylphenols + CO₂ up; carbon flux sums to 0; touches honoured."""
    y = _state(schema, X_brett=0.3, hydroxycinnamics=0.1)
    d = BrettDecarboxylation().derivatives(0.0, y, schema, params)

    assert d[schema.slice("hydroxycinnamics")][0] < 0.0
    assert d[schema.slice("vinylphenols")][0] > 0.0
    assert d[schema.slice("CO2")][0] > 0.0
    # Carbon flux across the three touched slots must cancel (9 C = 8 C + 1 C).
    c_flux = (
        d[schema.slice("hydroxycinnamics")][0] * carbon_mass_fraction("p_coumaric_acid")
        + d[schema.slice("vinylphenols")][0] * carbon_mass_fraction("vinylphenol")
        + d[schema.slice("CO2")][0] * carbon_mass_fraction("CO2")
    )
    assert c_flux == pytest.approx(0.0, abs=1e-12)
    assert set(BrettDecarboxylation.touches) == {"hydroxycinnamics", "vinylphenols", "CO2"}


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
    assert set(BrettVinylphenolReduction.touches) == {"vinylphenols", "ethylphenols"}


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
