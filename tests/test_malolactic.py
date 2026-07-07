"""Malolactic fermentation v1 — conversion-only (decision D-23).

Ranked headline-first. The keystone payoff is ``test_headline_mlf_raises_ph_emergently``:
the hand-built ``test_acidbase.test_headline_malic_to_lactic_raises_ph`` ΔpH ∈ [0.1, 0.3]
now *emerges* from the :class:`MalolacticConversion` Process on a malic-rich co-inoculated
must — measured as the no-MLF **control difference** ``pH_final(dosed) − pH_final(off)`` so
it isolates the malic→lactic swap from the (shared) Byp pH down-drift. The rest pin the
mole-for-mole stoichiometry + carbon closure, the gates (pH up ⇒ faster, ethanol/molecular
SO₂ ⇒ arrest — the first RHS consumers of D-18 pH and D-22 SO₂), the cardinal-temperature
optimum, the explicit ``speculative`` tier, and — prime directive #3 — isolability: an
undosed run is byte-for-byte the validated core *and* keeps ``malic``/``lactic`` at the
VALIDATED tier (the Process is disabled at compile when MLF is not pitched).
"""

import math
from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.analysis import ph_series
from fermentation.core import acidbase
from fermentation.core.chemistry import (
    M_CO2,
    M_LACTIC,
    M_MALIC,
    M_TARTARIC,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.malolactic import (
    MalolacticConversion,
    MalolacticDeath,
    MalolacticSenescence,
    cardinal_temperature_factor,
    malolactic_environmental_gate,
    malolactic_toxicity_gate,
)
from fermentation.core.media import wine_schema
from fermentation.core.state import StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario
from fermentation.units.convert import mgl_to_gpl
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)


@pytest.fixture
def pset():
    """Real wine kinetic params (incl. the MLF set) + the shared pKa set."""
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    return pset.resolve()


def _anchor_cation(pka, tartaric_gpl: float, malic_gpl: float, lactic_gpl: float, target_ph: float):
    totals = {
        "tartaric": tartaric_gpl / M_TARTARIC,
        "malic": malic_gpl / M_MALIC,
        "lactic": lactic_gpl / M_LACTIC,
    }
    return acidbase.solve_cation_charge(totals, 0.0, pka, target_ph)


def _wine_state(schema: StateSchema, params: Mapping[str, float], *, target_ph: float, **slots):
    """A wine state with the bulk vars, given acid/MLF slots, and the cation that anchors
    ``target_ph`` from the named acids (so pH gates see a realistic pH)."""
    tartaric = slots.get("tartaric", 0.0)
    malic = slots.get("malic", 0.0)
    lactic = slots.get("lactic", 0.0)
    cation = _anchor_cation(acidbase.build_pka_map(params), tartaric, malic, lactic, target_ph)
    base: dict[str, float | list[float]] = {
        "X": 1.0, "S": [120.0], "E": 0.0, "N": 0.1, "T": 293.15, "CO2": 0.0,
        "cation_charge": cation,
    }  # fmt: skip
    base.update(slots)
    return schema.pack(base)


def _wine_scenario(**initial_extra) -> Scenario:
    initial: dict[str, float] = {
        "brix": 24.0,
        "yan_mgl": 250.0,
        "pitch_gpl": 0.5,
        "tartaric_gpl": 4.0,
        "malic_gpl": 4.0,
        "initial_ph": 3.4,
    }
    initial.update(initial_extra)
    return Scenario(
        name="wine-mlf",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )


def _run(**initial_extra: float):
    compiled = compile_scenario(_wine_scenario(**initial_extra), strict=True)
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    return compiled, traj


# -- 1. HEADLINE: MLF deacidification emerges, ΔpH 0.1–0.3 (no-MLF control diff) ----


def test_headline_mlf_raises_ph_emergently():
    # The faithful emergent translation of the hand-built acidbase headline: dose O. oeni
    # on a malic-rich must and the pH rises into the deacidification band — measured as the
    # difference vs an otherwise-identical un-pitched run, so the (shared) Byp down-drift
    # cancels and the gap is purely the malic→lactic swap at the same final state.
    c_off, t_off = _run()
    c_on, t_on = _run(mlf_pitch_gpl=0.2)
    ph_off = ph_series(t_off, c_off.param_values)
    ph_on = ph_series(t_on, c_on.param_values)

    # MLF runs (co-inoculation) to near-completion in the early low-ethanol window
    assert t_off.series("malic")[-1] == pytest.approx(4.0)  # disabled when un-pitched
    assert t_on.series("malic")[-1] < 0.05  # malate essentially consumed
    assert t_on.series("lactic")[-1] > 2.0  # lactate produced

    delta = ph_on[-1] - ph_off[-1]
    assert 0.1 <= delta <= 0.3, f"emergent MLF ΔpH {delta:.3f} outside [0.1, 0.3]"


# -- 2. mole-for-mole stoichiometry + carbon closure at the RHS level --------------


def test_conversion_is_mole_for_mole_and_carbon_closing(params):
    schema = wine_schema()
    y = _wine_state(schema, params, target_ph=3.4, malic=4.0, tartaric=4.0, X_mlf=0.2)
    d = MalolacticConversion().derivatives(0.0, y, schema, params)
    dmalic = float(d[schema.slice("malic")][0])
    dlactic = float(d[schema.slice("lactic")][0])
    dco2 = float(d[schema.slice("CO2")][0])

    assert dmalic < 0.0 and dlactic > 0.0 and dco2 > 0.0
    # one mole of malate makes one mole of lactate and one of CO2
    r = -dmalic / M_MALIC
    assert dlactic / M_LACTIC == pytest.approx(r)
    assert dco2 / M_CO2 == pytest.approx(r)
    # carbon closes on the existing ledger (4 C = 3 C + 1 C): weighted sum is zero
    carbon_rate = (
        dmalic * carbon_mass_fraction("malic_acid")
        + dlactic * carbon_mass_fraction("lactic_acid")
        + dco2 * carbon_mass_fraction("CO2")
    )
    assert carbon_rate == pytest.approx(0.0, abs=1e-12)


def test_carbon_conserved_over_a_dosed_run():
    compiled, traj = _run(mlf_pitch_gpl=0.2)
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.parameters["biomass_C_fraction"].value
    )
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-9, label="total carbon (MLF on)")
    # the conversion self-limits (rate ∝ malate, guarded ≤0), so malate/lactate never go
    # meaningfully negative even as malate is driven to ~0
    assert_nonnegative(traj, ("malic", "lactic", "CO2"))


# -- 3. touches only malic/lactic/CO2 (strict contract) ---------------------------


def test_touches_only_malic_lactic_co2(params):
    schema = wine_schema()
    y = _wine_state(schema, params, target_ph=3.4, malic=4.0, tartaric=4.0, X_mlf=0.2)
    d = MalolacticConversion().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"malic", "lactic", "CO2"}


# -- 4. pH gate: higher pH ⇒ faster (the deacidification feedback) -----------------


def test_higher_ph_speeds_conversion(params):
    schema = wine_schema()
    # same malate, two anchored pHs: the rate rises with pH (O. oeni pH gate)
    y_lo = _wine_state(schema, params, target_ph=3.0, malic=4.0, tartaric=4.0, X_mlf=0.2)
    y_hi = _wine_state(schema, params, target_ph=3.6, malic=4.0, tartaric=4.0, X_mlf=0.2)
    rate_lo = -float(
        MalolacticConversion().derivatives(0.0, y_lo, schema, params)[schema.slice("malic")][0]
    )
    rate_hi = -float(
        MalolacticConversion().derivatives(0.0, y_hi, schema, params)[schema.slice("malic")][0]
    )
    assert 0.0 < rate_lo < rate_hi


# -- 5. ethanol gate: above O. oeni tolerance, MLF stops --------------------------


def test_ethanol_above_tolerance_arrests_mlf(params):
    schema = wine_schema()
    above = params["ethanol_tolerance_mlf"] + 5.0
    y = _wine_state(schema, params, target_ph=3.4, malic=4.0, tartaric=4.0, X_mlf=0.2, E=above)
    d = MalolacticConversion().derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("malic")][0]) == 0.0


# -- 6. molecular-SO₂ gate: dosing free SO₂ suppresses MLF (first RHS consumer) ----


def test_so2_dose_suppresses_mlf_unit(params):
    schema = wine_schema()
    common = {"target_ph": 3.4, "malic": 4.0, "tartaric": 4.0, "X_mlf": 0.2}
    y_clean = _wine_state(schema, params, **common)
    y_so2 = _wine_state(schema, params, so2_total=mgl_to_gpl(80.0), **common)
    rate_clean = -float(
        MalolacticConversion().derivatives(0.0, y_clean, schema, params)[schema.slice("malic")][0]
    )
    rate_so2 = -float(
        MalolacticConversion().derivatives(0.0, y_so2, schema, params)[schema.slice("malic")][0]
    )
    assert rate_clean > 0.0
    assert rate_so2 < 0.01 * rate_clean  # molecular SO₂ arrests the conversion


def test_so2_dose_suppresses_mlf_in_a_run():
    # The integration-level demonstration: SO₂ gates MLF as the first RHS consumer of molecular
    # SO₂ — but only *partially*, because of the emergent acetaldehyde–SO₂ competition made
    # dynamic by D-47. When SO₂ is dosed at pitch, the yeast's acetaldehyde binds it near-
    # stoichiometrically AND is then locked in (D-47: bound acetaldehyde is protected from ADH),
    # so free SO₂ stays chronically depressed for the whole run rather than recovering, and a
    # *substantial* slice of MLF proceeds, vs the uninhibited run consuming all of it and the
    # un-pitched run converting none. This is the faithful "bound SO₂ is not antimicrobial"
    # consequence — dosing SO₂ into acetaldehyde-rich fermenting must is a far weaker MLF brake
    # than the same dose in a cleared wine.
    #
    # POST-D-51 CAVEAT (real physics, not a regression): the multi-carbonyl equilibrium now also
    # lets the always-on pyruvate/α-KG pools compete for bisulfite, which ADDS binding capacity
    # beyond acetaldehyde alone (D-51 is a shared pool, not just a diversion) — so overall free
    # SO₂ ends lower than pre-D-51 (~15 % of the 80 mg/L dose here, down from ~21 %), the brake is
    # a little weaker still, and MLF now edges just past the halfway point (~51 % converted,
    # ~1.95 g/L malic remains) rather than retaining a majority.
    c_off, t_off = _run()
    c_on, t_on = _run(mlf_pitch_gpl=0.2)
    c_so2, t_so2 = _run(mlf_pitch_gpl=0.2, so2_total_mgl=80.0)
    ph_off = ph_series(t_off, c_off.param_values)
    ph_on = ph_series(t_on, c_on.param_values)
    ph_so2 = ph_series(t_so2, c_so2.param_values)

    assert t_off.series("malic")[-1] == pytest.approx(4.0)  # no bacteria ⇒ no conversion
    assert t_on.series("malic")[-1] < 0.5  # uninhibited MLF consumes ~all of it
    # SO₂ still suppresses substantially — but PARTIALLY, because stranded acetaldehyde plus the
    # competing keto-acid pools (D-51) sequester most of the antimicrobial pool (~1.7-2.2 g/L
    # remains, ~1.95 g/L nominal — just under half the malic converts, D-48's induced
    # over-production and the D-51 competition both weaken the brake a little further; the shipped
    # k's are sized to KEEP this "roughly half retained" regime, so the band is unchanged):
    assert 1.7 < t_so2.series("malic")[-1] < 2.2  # ~1.95 g/L remains (~51 % converted)
    assert t_so2.series("malic")[-1] > 5.0 * t_on.series("malic")[-1]  # …yet well above uninhibited
    assert (ph_on[-1] - ph_off[-1]) > 0.15  # uninhibited MLF deacidifies fully
    assert 0.03 < (ph_so2[-1] - ph_off[-1]) < 0.12  # SO₂ run deacidifies only partially


# -- 7. cardinal-temperature optimum (peaks at T_opt, declines warm, 0 outside) ----


def test_cardinal_temperature_factor_shape(params):
    t_min, t_opt, t_max = params["T_min_mlf"], params["T_opt_mlf"], params["T_max_mlf"]
    assert cardinal_temperature_factor(t_opt, t_min, t_opt, t_max) == pytest.approx(1.0)
    assert cardinal_temperature_factor(t_min, t_min, t_opt, t_max) == 0.0
    assert cardinal_temperature_factor(t_max, t_min, t_opt, t_max) == 0.0
    assert cardinal_temperature_factor(t_min - 5.0, t_min, t_opt, t_max) == 0.0
    assert cardinal_temperature_factor(t_max + 5.0, t_min, t_opt, t_max) == 0.0
    # peak at the optimum, in (0, 1) between cardinals, and DECLINES on the warm side
    # (the qualitative correctness a monotone Arrhenius factor lacks — decision D-23)
    warm = cardinal_temperature_factor(t_opt + 7.0, t_min, t_opt, t_max)
    cool = cardinal_temperature_factor(t_opt - 7.0, t_min, t_opt, t_max)
    assert 0.0 < warm < 1.0 and 0.0 < cool < 1.0
    grid = np.linspace(t_min + 0.1, t_max - 0.1, 200)
    vals = np.array([cardinal_temperature_factor(t, t_min, t_opt, t_max) for t in grid])
    assert vals.max() == pytest.approx(1.0, abs=1e-3)
    assert np.argmax(vals) == int(np.argmin(np.abs(grid - t_opt)))  # unimodal peak at T_opt


# -- 8. isolability (prime directive #3): undosed = validated core, byte + tier ----


def test_undosed_run_is_byte_for_byte_validated_core():
    # mlf_pitch_gpl absent vs explicitly 0 must integrate identically (the Process is
    # disabled both ways), so MLF adds nothing to an un-pitched wine ferment.
    t_eval = np.linspace(0.0, 21.0 * 24.0, 80)
    c_absent = compile_scenario(_wine_scenario(), strict=True)
    c_zero = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.0), strict=True)
    ta = simulate(
        c_absent.process_set, c_absent.param_values, c_absent.y0, c_absent.t_span_h, t_eval=t_eval
    )
    tz = simulate(
        c_zero.process_set, c_zero.param_values, c_zero.y0, c_zero.t_span_h, t_eval=t_eval
    )
    for name in c_absent.schema.names:
        assert np.allclose(ta.series(name), tz.series(name), rtol=1e-12, atol=1e-12), name
    assert not c_absent.process_set.is_enabled(MalolacticConversion.name)


def test_undosed_keeps_malic_lactic_validated_dosed_makes_speculative(pset):
    tm = pset.tier_map()
    off = compile_scenario(_wine_scenario()).process_set
    on = compile_scenario(_wine_scenario(mlf_pitch_gpl=0.2)).process_set
    # un-pitched: nothing active touches the acid slots, so they stay VALIDATED
    assert off.tier_of("malic", tm) is Tier.VALIDATED
    assert off.tier_of("lactic", tm) is Tier.VALIDATED
    # pitched: the speculative MLF Process touches them ⇒ speculative (prime directive #1)
    assert on.tier_of("malic", tm) is Tier.SPECULATIVE
    assert on.tier_of("lactic", tm) is Tier.SPECULATIVE


def test_so2_crashes_bacteria_over_the_slow_senescence_baseline():
    # The v2 (D-41) picture, superseding the v1 "no-SO₂ pitched run is byte-for-byte inert" premise:
    # on a pitched, amino-acid-free run (growth disabled) the two movers of X_mlf are the benign
    # senescence baseline (MalolacticSenescence, always on when pitched — since D-52 modulated by a
    # BOUNDED ethanol/starvation stress factor, up to 2.5× as fermentation ethanol rises and the
    # amino-acid pool empties, D-52) and the SO₂-driven kill (MalolacticDeath, 0 without SO₂). So
    # WITHOUT SO₂ the catalyst declines SLOWLY (weeks-scale senescence, faster than the D-41
    # environment-free baseline but still no longer inert), and an ``add_so2`` dose CRASHES it
    # toward zero within ~1–3 days (the D-31 lever) — still a rate orders of magnitude above even
    # the fully-stressed baseline, so the two timescales stay cleanly separated.
    # (a) no SO₂ ⇒ X_mlf declines slowly & monotonically (senescence), retaining most over the run.
    _c, t_clean = _run(mlf_pitch_gpl=0.2)
    x_clean = t_clean.series("X_mlf")
    assert x_clean[0] == pytest.approx(0.2)  # the pitched dose
    assert np.all(np.diff(x_clean) <= 1e-12)  # monotone non-increasing — senescence only, no SO₂
    assert x_clean[-1] < 0.98 * x_clean[0]  # NOT inert: it has measurably declined (v1 superseded)
    # MEASURED ratio ~0.608 (D-52: stress ramps 1.5x -> ~2.0x as ethanol rises then holds, vs the
    # D-41 environment-free ~0.71 at the same 21 d) — banded around the measured value, not just
    # threshold-checked, so a regression in the stress-multiplier magnitude is caught either way.
    assert 0.55 * x_clean[0] < x_clean[-1] < 0.65 * x_clean[0]

    # (b) add SO₂ mid-run ⇒ gentle-but-D-52-stressed senescence until the dose, then a sharp crash.
    dosed = Scenario(
        name="wine-mlf-so2",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.5,
            "tartaric_gpl": 4.0,
            "malic_gpl": 4.0,
            "initial_ph": 3.4,
            "mlf_pitch_gpl": 0.2,
        },  # fmt: skip
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[Intervention(day=6.0, action="add_so2", params={"so2_mgl": 50.0})],
        duration_days=21.0,
    )
    traj = compile_scenario(dosed, strict=True).run()
    t_h, x_mlf = traj.t, traj.series("X_mlf")
    i6 = int(np.searchsorted(t_h, 6.0 * 24.0))  # index just at/after the day-6 SO₂ dose
    # MEASURED ratio ~0.875 at day 6 (D-52 stress is already ramping as AF ethanol rises and
    # amino acids empty — faster than D-41's environment-free ~0.95, still far above the crash).
    assert 0.8 * 0.2 < x_mlf[i6] < 0.95 * 0.2
    post = x_mlf[i6:]
    assert np.all(np.diff(post) <= 1e-12)  # monotone non-increasing once SO₂ is present
    assert x_mlf[-1] < 0.1 * x_mlf[i6]  # SO₂ crashed the population (fast kill dominates baseline)


# -- 8b. MalolacticDeath — the SO₂-driven bacterial kill (decision D-39) -----------


def _death_state(
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    so2_mgl: float = 0.0,
    temp_k: float = 293.15,
    x_mlf: float = 0.2,
    x_mlf_dead: float = 0.0,
):
    """A pitched wine state for exercising MalolacticDeath at the RHS level."""
    slots = {"X_mlf": x_mlf, "X_mlf_dead": x_mlf_dead, "T": temp_k, "malic": 2.0, "tartaric": 4.0}
    if so2_mgl > 0.0:
        slots["so2_total"] = mgl_to_gpl(so2_mgl)
    return _wine_state(schema, params, target_ph=3.4, **slots)


def test_death_is_exactly_zero_without_so2(params):
    # The v1 tradeoff, enforced: death is driven by molecular SO₂ ALONE, so an unsulfited pitched
    # run has an identically-zero death contribution (O. oeni persists) — no ethanol/pH decay.
    schema = wine_schema()
    y = _death_state(schema, params, so2_mgl=0.0, x_mlf=0.2)
    d = MalolacticDeath().derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("X_mlf")][0]) == 0.0
    assert float(d[schema.slice("X_mlf_dead")][0]) == 0.0


def test_so2_drives_death_as_a_neutral_transfer(params):
    # With SO₂ dosed, viable X_mlf leaves and the SAME mass enters X_mlf_dead — the D-13 X→X_dead
    # transfer. d[X_mlf] = −d[X_mlf_dead] exactly, so (both weighted at the biomass fractions since
    # D-38) the move is carbon- and nitrogen-neutral by construction.
    schema = wine_schema()
    y = _death_state(schema, params, so2_mgl=80.0, x_mlf=0.2)
    d = MalolacticDeath().derivatives(0.0, y, schema, params)
    dx = float(d[schema.slice("X_mlf")][0])
    dxd = float(d[schema.slice("X_mlf_dead")][0])
    assert dx < 0.0 and dxd > 0.0  # bacteria die
    assert dxd == pytest.approx(-dx)  # mass-conserving transfer (neutral in both ledgers)


def test_death_touches_only_the_x_mlf_pools(params):
    schema = wine_schema()
    y = _death_state(schema, params, so2_mgl=80.0, x_mlf=0.2)
    d = MalolacticDeath().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"X_mlf", "X_mlf_dead"}


def test_more_so2_kills_faster(params):
    # Monotone in the antimicrobial dose: a larger SO₂ addition ⇒ higher molecular SO₂ ⇒ a larger
    # 1 − g_SO₂ ⇒ a faster kill. (The D-31 lever's strength scales with the sulfite dose.)
    schema = wine_schema()
    rate_lo = -float(
        MalolacticDeath().derivatives(
            0.0, _death_state(schema, params, so2_mgl=20.0), schema, params
        )[schema.slice("X_mlf")][0]
    )
    rate_hi = -float(
        MalolacticDeath().derivatives(
            0.0, _death_state(schema, params, so2_mgl=60.0), schema, params
        )[schema.slice("X_mlf")][0]
    )
    assert 0.0 < rate_lo < rate_hi


def test_cold_preserves_bacteria_via_arrhenius_not_gamma(params):
    # The load-bearing D-39 choice: death carries its OWN Arrhenius factor, not the cardinal γ(T).
    # Below T_min_mlf (8 °C) γ(T) = 0 — if death reused γ(T), cold would spuriously HALT the kill.
    # Instead Arrhenius merely SLOWS it: a sulfited culture in the cold still dies (just slower),
    # and warmer ⇒ faster. Cold preserving bacteria is the correct direction γ(T) cannot supply.
    schema = wine_schema()
    t_below_min = 278.15  # 5 °C, below T_min_mlf = 281.15 K (8 °C) ⇒ cardinal γ(T) = 0
    assert cardinal_temperature_factor(t_below_min, 281.15, 296.15, 310.15) == 0.0
    rate_cold = -float(
        MalolacticDeath().derivatives(
            0.0, _death_state(schema, params, so2_mgl=80.0, temp_k=t_below_min), schema, params
        )[schema.slice("X_mlf")][0]
    )
    rate_warm = -float(
        MalolacticDeath().derivatives(
            0.0, _death_state(schema, params, so2_mgl=80.0, temp_k=298.15), schema, params
        )[schema.slice("X_mlf")][0]
    )
    assert rate_cold > 0.0  # still dying below T_min — proves Arrhenius, not γ(T)
    assert rate_warm > rate_cold  # warm accelerates the kill (Arrhenius direction)


def test_death_run_conserves_carbon_and_nitrogen():
    # Integration-level closure with death ACTIVE: pitch O. oeni, add SO₂ mid-run so X_mlf → dead.
    # The X_mlf → X_mlf_dead transfer is carbon/nitrogen-neutral (both weighted, D-38), SO₂ carries
    # neither element, and malate→lactate+CO₂ is carbon-closing — so both ledgers close to machine
    # precision across the sulfite jump (final == initial + Σ external flows; SO₂ adds no C/N).
    scen = Scenario(
        name="wine-mlf-death",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.5,
            "tartaric_gpl": 4.0,
            "malic_gpl": 4.0,
            "initial_ph": 3.4,
            "mlf_pitch_gpl": 0.2,
        },  # fmt: skip
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[Intervention(day=6.0, action="add_so2", params={"so2_mgl": 60.0})],
        duration_days=21.0,
    )
    cs = compile_scenario(scen, strict=True)
    traj = cs.run()
    f_c = cs.param_values["biomass_C_fraction"]
    f_n = cs.param_values["biomass_N_fraction"]
    carbon = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    nitrogen = total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n)
    # SO₂ is carbon/nitrogen-free, so the run-wide balances close with no external-flow correction.
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-9, label="total carbon (MLF death on)")
    assert_conserved(traj, nitrogen, rtol=1e-6, atol=1e-9, label="total nitrogen (MLF death on)")
    assert_nonnegative(traj, ("X_mlf", "X_mlf_dead"))


def test_death_tier_is_speculative():
    assert MalolacticDeath.tier is Tier.SPECULATIVE


# -- 8c. MalolacticSenescence — the benign baseline mortality (MLF v2, decision D-41) --


def test_senescence_is_a_neutral_transfer_without_so2(params):
    # The v2 headline, DISTINGUISHING senescence from the SO₂ kill: even with NO SO₂, viable X_mlf
    # ages into X_mlf_dead — a nonzero contribution where MalolacticDeath is identically zero. It is
    # the same D-13 neutral transfer: d[X_mlf] = −d[X_mlf_dead] exactly (both weighted at the
    # biomass fractions ⇒ carbon- and nitrogen-neutral), so it adds no conservation code.
    schema = wine_schema()
    y = _death_state(schema, params, so2_mgl=0.0, x_mlf=0.2)
    d = MalolacticSenescence().derivatives(0.0, y, schema, params)
    dx = float(d[schema.slice("X_mlf")][0])
    dxd = float(d[schema.slice("X_mlf_dead")][0])
    assert dx < 0.0 and dxd > 0.0  # bacteria age/die even with no SO₂ (v1 "never die" superseded)
    assert dxd == pytest.approx(-dx)  # mass-conserving transfer (neutral in both ledgers)


def test_senescence_is_zero_without_bacteria(params):
    schema = wine_schema()
    y = _death_state(schema, params, x_mlf=0.0)
    d = MalolacticSenescence().derivatives(0.0, y, schema, params)
    assert float(d[schema.slice("X_mlf")][0]) == 0.0
    assert float(d[schema.slice("X_mlf_dead")][0]) == 0.0


def test_senescence_touches_only_the_x_mlf_pools(params):
    schema = wine_schema()
    y = _death_state(schema, params, x_mlf=0.2)
    d = MalolacticSenescence().derivatives(0.0, y, schema, params)
    touched = {n for n in schema.names if np.any(d[schema.slice(n)] != 0.0)}
    assert touched == {"X_mlf", "X_mlf_dead"}


def test_senescence_is_so2_independent(params):
    # Retained from D-41's environment-free property (the D-39 crux reused): the benign baseline
    # still carries NO SO₂ term at all — only the acute :class:`MalolacticDeath` kill sees SO₂. D-52
    # adds ethanol/starvation stress (see below) but SO₂ independence is untouched.
    schema = wine_schema()

    def rate(y) -> float:
        d = MalolacticSenescence().derivatives(0.0, y, schema, params)
        return -float(d[schema.slice("X_mlf")][0])

    r_no_so2 = rate(_death_state(schema, params, so2_mgl=0.0, x_mlf=0.2))
    r_hi_so2 = rate(_death_state(schema, params, so2_mgl=80.0, x_mlf=0.2))
    assert r_no_so2 > 0.0
    assert r_hi_so2 == pytest.approx(r_no_so2, rel=1e-12)


def test_senescence_ethanol_stress_is_bounded(params):
    # D-52 lifts the D-41 "environment-free" deferral: senescence now RISES with ethanol stress
    # (a real effect the D-41 baseline knowingly omitted), but the rise is a smooth, BOUNDED
    # Monod-type factor — not the Luong wall's near-binary collapse that caused the D-39 wipeout. So
    # higher ethanol must raise the rate, but the ratio against the zero-stress floor can never
    # exceed the design ceiling (1 + k_senescence_ethanol_scale + k_senescence_starvation_scale).
    schema = wine_schema()

    def rate(e: float, aa: float) -> float:
        y = _wine_state(
            schema, params, target_ph=3.4, X_mlf=0.2, T=293.15, E=e, malic=2.0, amino_acids=aa
        )
        d = MalolacticSenescence().derivatives(0.0, y, schema, params)
        return -float(d[schema.slice("X_mlf")][0])

    floor = rate(0.0, 10.0)  # no ethanol, amino acids abundant ⇒ stress == 1 (pure baseline)
    lo = rate(10.0, 10.0)
    hi = rate(120.0, 10.0)  # at/above ethanol_tolerance_mlf (~110 g/L) — the O. oeni ethanol wall
    assert floor > 0.0
    assert floor < lo < hi  # monotone rising with ethanol stress (the D-52 lift)
    ceiling = 1.0 + params["k_senescence_ethanol_scale"] + params["k_senescence_starvation_scale"]
    assert hi / floor < ceiling  # BOUNDED — never approaches the full ceiling at finite E


def test_senescence_starvation_stress_tracks_amino_acid_depletion(params):
    # D-52's second stress term reuses the growth fuel pool: senescence is faster once amino acids
    # are exhausted (the ~1.3 d post-pitch norm, D-23) than while the pool is replete — the term
    # tracks the real nutrient dynamic rather than acting as a flat multiplier.
    schema = wine_schema()

    def rate(aa: float) -> float:
        y = _wine_state(
            schema, params, target_ph=3.4, X_mlf=0.2, T=293.15, E=0.0, malic=2.0, amino_acids=aa
        )
        d = MalolacticSenescence().derivatives(0.0, y, schema, params)
        return -float(d[schema.slice("X_mlf")][0])

    r_starved = rate(0.0)
    r_replete = rate(10.0)  # well above K_aa_mlf (0.05 g/L) ⇒ starvation term ≈ 0
    assert r_starved > r_replete > 0.0


def test_senescence_no_wipeout_at_worst_case_stress_at_benchmark_temperature(params):
    # The D-52 wipeout guard, verified empirically (not just threshold-checked, per the D-51/D-48
    # project discipline): at the wine benchmark T_ref (20 °C, Arrhenius factor 1) — the scope this
    # bound covers — even the WORST case simultaneously (saturating ethanol AND amino acids fully
    # exhausted) keeps the half-life in the weeks range, nowhere near the ~1-week D-39 wipeout
    # regime that deferred ethanol coupling to v2 in the first place. Temperature is a SEPARATE
    # stress axis (a warm cellar legitimately shortens this further via Arrhenius — see
    # ``test_senescence_warm_worst_case_stress_stays_far_below_the_so2_kill`` for that axis).
    schema = wine_schema()
    y = _wine_state(
        schema, params, target_ph=3.4, X_mlf=0.2, T=293.15, E=1.0e4, malic=2.0, amino_acids=0.0
    )
    d = MalolacticSenescence().derivatives(0.0, y, schema, params)
    r_worst = -float(d[schema.slice("X_mlf")][0])
    specific_rate = r_worst / 0.2  # [1/h], independent of the X_mlf dose used above
    half_life_h = math.log(2.0) / specific_rate
    assert half_life_h > 14.0 * 24.0  # > 2 weeks at T_ref — clear of the ~1-week wipeout regime


def test_senescence_warm_worst_case_stress_stays_far_below_the_so2_kill(params):
    # Temperature is a separate stress axis from ethanol/starvation (Arrhenius, the D-39/D-41
    # precedent): a warm cellar legitimately shortens the worst-case combined-stress half-life
    # further than the T_ref bound above — that is correct physics (warm accelerates senescence),
    # NOT a re-emergence of the D-39 ethanol-at-normal-temperature wipeout bug. The invariant that
    # actually holds at ANY temperature is the one D-41 established (test_senescence_is_slow_
    # relative_to_the_so2_kill): the chronic senescence baseline, even fully stressed, stays far
    # below the acute SO₂-driven kill — both share the same arrhenius(T, E_a_death_mlf, T_ref)
    # factor, so the ratio is temperature-invariant by construction, verified here at a warm 30 °C.
    schema = wine_schema()
    warm = 303.15  # 30 °C
    y_sen = _wine_state(
        schema, params, target_ph=3.4, X_mlf=0.2, T=warm, E=1.0e4, malic=2.0, amino_acids=0.0
    )
    r_worst_warm = -float(
        MalolacticSenescence().derivatives(0.0, y_sen, schema, params)[schema.slice("X_mlf")][0]
    )
    y_kill = _death_state(schema, params, so2_mgl=80.0, x_mlf=0.2, temp_k=warm)
    r_kill_warm = -float(
        MalolacticDeath().derivatives(0.0, y_kill, schema, params)[schema.slice("X_mlf")][0]
    )
    assert 0.0 < r_worst_warm < 0.1 * r_kill_warm  # chronic stays far below the kill, even warm


def test_senescence_warm_accelerates_cold_preserves_via_arrhenius(params):
    # Like MalolacticDeath, senescence carries the Arrhenius factor, NOT the cardinal γ(T): below
    # T_min_mlf (γ(T) = 0) it still proceeds (cold merely SLOWS it — dormancy, not immortality),
    # and warm accelerates it. γ(T) would spuriously make senescence PEAK at the 23 °C growth
    # optimum and switch OFF in the warm — backwards for a decline (the D-41 rationale).
    schema = wine_schema()

    def rate(tk: float) -> float:
        y = _death_state(schema, params, temp_k=tk, x_mlf=0.2)
        d = MalolacticSenescence().derivatives(0.0, y, schema, params)
        return -float(d[schema.slice("X_mlf")][0])

    t_below_min = 278.15  # 5 °C, below T_min_mlf = 8 °C ⇒ cardinal γ(T) = 0
    assert cardinal_temperature_factor(t_below_min, 281.15, 296.15, 310.15) == 0.0
    assert rate(t_below_min) > 0.0  # still aging below T_min — proves Arrhenius, not γ(T)
    assert rate(298.15) > rate(t_below_min)  # warm accelerates senescence (Arrhenius direction)


def test_senescence_is_slow_relative_to_the_so2_kill(params):
    # The two timescales are cleanly separated by ~100×: the benign baseline is far below a full SO₂
    # kill, so a stabilizing SO₂ dose dominates the trajectory (crash in days) while the baseline is
    # the weeks-to-months persistence limit. Compared at the same X_mlf and temperature.
    schema = wine_schema()
    r_sen = -float(
        MalolacticSenescence().derivatives(
            0.0, _death_state(schema, params, so2_mgl=0.0, x_mlf=0.2), schema, params
        )[schema.slice("X_mlf")][0]
    )
    r_kill = -float(
        MalolacticDeath().derivatives(
            0.0, _death_state(schema, params, so2_mgl=80.0, x_mlf=0.2), schema, params
        )[schema.slice("X_mlf")][0]
    )
    assert 0.0 < r_sen < 0.1 * r_kill  # baseline is a small fraction of the acute SO₂ kill


def test_senescence_run_conserves_carbon_and_nitrogen():
    # Integration-level closure with senescence ACTIVE and NO SO₂: pitch O. oeni, leave it untreated
    # so the ONLY X_mlf mover is the senescence baseline (growth off — no amino acids; death 0 — no
    # SO₂). The X_mlf → X_mlf_dead transfer is carbon/nitrogen-neutral (both weighted, D-38), so
    # both ledgers close to machine precision with no external flow (nothing is dosed after t0).
    scen = Scenario(
        name="wine-mlf-senescence",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 250.0,
            "pitch_gpl": 0.5,
            "tartaric_gpl": 4.0,
            "malic_gpl": 4.0,
            "initial_ph": 3.4,
            "mlf_pitch_gpl": 0.2,
        },  # fmt: skip
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=21.0,
    )
    cs = compile_scenario(scen, strict=True)
    assert cs.process_set.is_enabled("malolactic_senescence")  # pitch-gated, enabled by the pitch
    traj = cs.run()
    x_mlf = traj.series("X_mlf")
    assert x_mlf[-1] < x_mlf[0]  # senescence has moved viable biomass into the dead pool
    f_c = cs.param_values["biomass_C_fraction"]
    f_n = cs.param_values["biomass_N_fraction"]
    carbon = total_carbon(cs.schema, biomass_carbon_fraction=f_c)
    nitrogen = total_nitrogen(cs.schema, biomass_nitrogen_fraction=f_n)
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-9, label="total carbon (MLF senescence on)")
    assert_conserved(traj, nitrogen, rtol=1e-6, atol=1e-9, label="total N (MLF senescence on)")
    assert_nonnegative(traj, ("X_mlf", "X_mlf_dead"))


def test_senescence_needs_no_ph_solve(params):
    # Performance/isolability: unlike MalolacticDeath (which reads molecular SO₂ at the solved pH),
    # senescence reads no SO₂ and no pH — even with the D-52 stress terms, E/amino_acids are read
    # directly off state — so its reads tuple excludes every acidbase/SO₂ parameter and it never
    # triggers a brentq. Pin the declared reads so a future SO₂/pH coupling can't slip in.
    assert set(MalolacticSenescence.reads) == {
        "k_senescence_mlf",
        "E_a_death_mlf",
        "T_ref",
        "k_senescence_ethanol_scale",
        "ethanol_tolerance_mlf",
        "k_senescence_starvation_scale",
        "K_aa_mlf",
    }


def test_senescence_tier_is_speculative():
    assert MalolacticSenescence.tier is Tier.SPECULATIVE


def test_environmental_gate_is_toxicity_times_gamma(params):
    # The D-39 refactor invariant: splitting the gate into toxicity × γ(T) left the growth/
    # conversion consumers byte-for-byte — the full gate is exactly the product of its two halves.
    schema = wine_schema()
    y = _wine_state(schema, params, target_ph=3.4, malic=2.0, tartaric=4.0, X_mlf=0.2, E=40.0)
    ph = acidbase.ph_of_state(y, schema, params)
    toxicity = malolactic_toxicity_gate(y, schema, params, ph)
    gamma = cardinal_temperature_factor(
        float(y[schema.slice("T")][0]),
        params["T_min_mlf"],
        params["T_opt_mlf"],
        params["T_max_mlf"],
    )
    full = malolactic_environmental_gate(y, schema, params, ph)
    assert full == toxicity * gamma  # exact — same multiplication grouping as before the split


# -- 9. tier is speculative -------------------------------------------------------


def test_process_tier_is_speculative():
    assert MalolacticConversion.tier is Tier.SPECULATIVE
