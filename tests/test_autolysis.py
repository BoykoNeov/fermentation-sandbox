"""Tests for yeast autolysis — the autolytic-peptide source (decision D-34).

MLF with bacterial growth is blocked because the ``amino_acids`` pool is empty at the MLF pitch
point (D-23). :class:`YeastAutolysis` refills it: the first consumer of the ``X_dead`` pool (dead
biomass, from D-13 ethanol inactivation), it liberates the dead-cell nitrogen as assimilable amino
acids (arginine) and routes the carbon-rich remainder to a non-assimilable ``debris`` pool (glucan).
Because dead biomass is carbon-rich (C:N ≈ 4–11) and arginine is nitrogen-rich (C:N ≈ 1.29), most of
the dead-cell carbon cannot leave as amino acids — it is cell-wall debris (the *sur lie* lees), so
carbon and nitrogen close *separately*. This suite pins that closure, the debris split, the
X_dead-guard, opt-in isolability, and the headline emergent behaviour: amino acids rise post-AF.
"""

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics import YeastAutolysis
from fermentation.core.media import get_medium
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)

AUTOLYSIS = YeastAutolysis.name
_AA_SPECIES = "arginine"
_DEBRIS_SPECIES = "glucan"


@pytest.fixture
def full_params():
    base = default_data_dir()
    return load_parameters(
        base / "wine_generic.yaml",
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "hydrogen_sulfide.yaml",
    ).resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x_dead: float = 1.5,
    t_k: float = 293.15,
) -> FloatArray:
    # A post-AF-ish state: some dead biomass to autolyse, ethanol high, sugar low, N exhausted.
    return schema.pack(
        {"X": 0.5, "S": [5.0], "E": 100.0, "N": 0.0, "T": t_k, "CO2": 100.0, "X_dead": x_dead}
    )


def _isolate_autolysis(full_params: dict[str, float]) -> ProcessSet:
    """A wine ProcessSet with only autolysis active, so the X_dead/amino_acids/debris columns
    reflect exactly its release + debris split."""
    ps = get_medium("wine").build_process_set()
    for p in ps.active:
        if p.name != AUTOLYSIS:
            ps.disable(p.name)
    return ps


# -- metadata -----------------------------------------------------------------


def test_metadata():
    p = YeastAutolysis()
    assert p.name == "yeast_autolysis"
    assert p.tier is Tier.SPECULATIVE
    assert set(p.touches) == {"X_dead", "amino_acids", "debris"}
    for r in ("k_autolysis", "E_a_autolysis", "biomass_N_fraction", "biomass_C_fraction"):
        assert r in p.reads


# -- carbon + nitrogen close separately; debris carries the excess carbon -----


def test_autolysis_contribution_is_carbon_and_nitrogen_neutral(full_params):
    # The heart of the design: autolysis moves carbon X_dead → amino_acids + debris and nitrogen
    # X_dead → amino_acids, so its OWN contribution changes neither total carbon nor total nitrogen.
    ps = _isolate_autolysis(full_params)
    schema = ps.schema
    f_c = full_params["biomass_C_fraction"]
    f_n = full_params["biomass_N_fraction"]
    carbon = total_carbon(schema, biomass_carbon_fraction=f_c)
    nitrogen = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)
    for x_dead in (0.3, 1.5, 4.0):
        y = _wine_y0(schema, x_dead=x_dead)
        d = ps.total_derivatives(0.0, y, full_params)
        assert abs(carbon(d)) < 1e-14
        assert abs(nitrogen(d)) < 1e-14


def test_debris_carries_the_excess_carbon_and_is_positive(full_params):
    # Nitrogen-anchored: amino acids gain exactly the dead-cell nitrogen; the C-rich remainder goes
    # to debris. Verify the algebra and that the excess is strictly positive (biomass is always more
    # carbon-rich than arginine, so the split never flips — no clamp needed).
    ps = _isolate_autolysis(full_params)
    schema = ps.schema
    y = _wine_y0(schema, x_dead=1.5)
    d = ps.total_derivatives(0.0, y, full_params)
    f_c = full_params["biomass_C_fraction"]
    f_n = full_params["biomass_N_fraction"]
    y_n = nitrogen_mass_fraction(_AA_SPECIES)
    y_c = carbon_mass_fraction(_AA_SPECIES)
    r = -d[schema.slice("X_dead")][0]  # X_dead consumed
    assert r > 0.0
    # amino acids carry exactly the dead-cell nitrogen:
    assert d[schema.slice("amino_acids")][0] == pytest.approx(r * f_n / y_n, rel=1e-12)
    assert d[schema.slice("amino_acids")][0] * y_n == pytest.approx(
        r * f_n, rel=1e-12
    )  # N released
    # debris carbon = dead-cell carbon − amino-acid carbon, and it is positive:
    debris_carbon = d[schema.slice("debris")][0] * carbon_mass_fraction(_DEBRIS_SPECIES)
    assert debris_carbon > 0.0
    assert debris_carbon == pytest.approx(r * (f_c - f_n * y_c / y_n), rel=1e-12)
    # the dominant fate really is debris (most dead-cell carbon is non-assimilable cell wall):
    aa_carbon = d[schema.slice("amino_acids")][0] * y_c
    assert debris_carbon > aa_carbon


def test_no_dead_cells_is_a_noop(full_params):
    # First consumer of X_dead: with no dead biomass there is nothing to autolyse (also guards the
    # clamp against a negative solver excursion).
    ps = _isolate_autolysis(full_params)
    schema = ps.schema
    y = _wine_y0(schema, x_dead=0.0)
    d = ps.total_derivatives(0.0, y, full_params)
    assert np.max(np.abs(d)) == 0.0


def test_temperature_accelerates_autolysis(full_params):
    # Enzymatic ⇒ warmer lees autolyse faster (E_a_autolysis > 0): X_dead consumption rises with T.
    ps = _isolate_autolysis(full_params)
    schema = ps.schema
    cold = ps.total_derivatives(0.0, _wine_y0(schema, x_dead=1.5, t_k=283.15), full_params)
    warm = ps.total_derivatives(0.0, _wine_y0(schema, x_dead=1.5, t_k=303.15), full_params)
    assert -warm[schema.slice("X_dead")][0] > -cold[schema.slice("X_dead")][0]


# -- isolability (opt-in) -----------------------------------------------------


def test_disabled_autolysis_is_byte_for_byte_core(full_params):
    # A bare wine build with autolysis DISABLED equals one where it is simply absent from the sum,
    # exactly, on every column — the opt-in guarantee (the compile seam disables it when undosed).
    on = get_medium("wine").build_process_set()
    off = get_medium("wine").build_process_set()
    off.disable(AUTOLYSIS)
    schema = on.schema
    # A live-ferment state (X_dead present) so the difference would be nonzero if it leaked:
    y = _wine_y0(schema, x_dead=1.5)
    on.disable(AUTOLYSIS)
    diff = on.total_derivatives(0.0, y, full_params) - off.total_derivatives(0.0, y, full_params)
    assert np.max(np.abs(diff)) == 0.0


def test_compile_seam_toggles_and_overrides():
    _, undosed = _run(150.0)
    assert AUTOLYSIS in undosed.process_set
    assert not undosed.process_set.is_enabled(AUTOLYSIS)
    _, dosed = _run(150.0, autolysis_rate_per_h=1e-2)
    assert dosed.process_set.is_enabled(AUTOLYSIS)
    assert dosed.param_values["k_autolysis"] == pytest.approx(1e-2)  # scenario overrides the YAML


def test_negative_autolysis_rate_raises():
    with pytest.raises(ValueError, match="autolysis_rate_per_h"):
        _run(150.0, autolysis_rate_per_h=-1.0)


# -- tier: structural drop only when enabled ---------------------------------


def test_autolysis_drops_output_tiers_only_when_enabled(full_params):
    # The speculative Process touches X_dead/amino_acids/debris; enabling it drops their output tier
    # to SPECULATIVE, but only when enabled (a disabled Process is excluded from tier derivation).
    schema = get_medium("wine").schema
    procs = [YeastAutolysis()]
    on = ProcessSet(schema, procs)
    off = ProcessSet(schema, procs)
    off.disable(AUTOLYSIS)
    for var in ("X_dead", "amino_acids", "debris"):
        assert on.tier_of(var) is Tier.SPECULATIVE
        assert off.tier_of(var) is Tier.VALIDATED  # nothing enabled touches it


# -- behaviour through the compile seam --------------------------------------


def _run(
    yan_mgl: float,
    *,
    autolysis_rate_per_h: float | None = None,
    amino_acids_gpl: float | None = None,
    days: float = 30.0,
):
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if autolysis_rate_per_h is not None:
        initial["autolysis_rate_per_h"] = autolysis_rate_per_h
    if amino_acids_gpl is not None:
        initial["amino_acids_gpl"] = amino_acids_gpl
    scenario = Scenario(
        name=f"wine-autolysis-{yan_mgl:.0f}",
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


def test_headline_amino_acids_rise_post_af_from_autolysis():
    # THE payoff: with autolysis on, dead biomass (X_dead) accumulates as the ferment ends and
    # then feeds the amino-acid pool, which rises from empty — the pool a later MLF-with-growth
    # model will draw on. Undosed amino acids, so nothing consumes it (consumers disabled): a rise.
    traj, _ = _run(150.0, autolysis_rate_per_h=1e-2)
    aa = np.asarray(traj.series("amino_acids"))
    x_dead = np.asarray(traj.series("X_dead"))
    debris = np.asarray(traj.series("debris"))
    assert aa[0] == 0.0
    assert x_dead.max() > 0.1  # dead biomass really accumulates (the source)
    assert aa[-1] > 1e-3  # materially refilled from empty
    assert aa[-1] > aa[len(aa) // 2]  # still rising in the post-AF tail
    assert debris[-1] > aa[-1]  # most autolysed carbon is non-assimilable debris (C-rich biomass)
    assert_nonnegative(traj, ("amino_acids", "debris", "X_dead"))


def test_carbon_and_nitrogen_close_with_autolysis_on():
    # Crown jewel: over a full autolysis-on run, carbon (X_dead → amino_acids + debris, weighted)
    # and nitrogen (X_dead → amino_acids) both close to solver tolerance.
    traj, compiled = _run(150.0, autolysis_rate_per_h=1e-2)
    f_c = compiled.param_values["biomass_C_fraction"]
    f_n = compiled.param_values["biomass_N_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj, total_nitrogen(compiled.schema, biomass_nitrogen_fraction=f_n), label="nitrogen"
    )


def test_autolysis_feeds_the_pool_while_swap_and_reroute_drain_it():
    # The three-way composition these two beats (D-33 + D-34) exist to enable — the actual
    # MLF-growth-prerequisite configuration, which every other test isolates apart. Dose BOTH the
    # amino-acid ledger (enabling the AminoAcidAssimilation swap AND the FuselAminoAcidReroute) AND
    # autolysis (refilling the pool from X_dead): all three touch ``amino_acids`` at once —
    # autolysis *feeds* it while the swap + re-route *drain* it. Carbon and nitrogen must still
    # close over the full run and the pools stay non-negative. Conservation as a TEST of the
    # composition, not the
    # assumption that three individually-neutral contributions sum to neutral (conservation is a
    # TEST, not an assumption — the CLAUDE.md conservation-laws-are-tests convention).
    traj, compiled = _run(150.0, autolysis_rate_per_h=1e-2, amino_acids_gpl=2.0)
    ps = compiled.process_set
    assert ps.is_enabled("amino_acid_assimilation")
    assert ps.is_enabled("fusel_amino_acid_reroute")
    assert ps.is_enabled(AUTOLYSIS)
    f_c = compiled.param_values["biomass_C_fraction"]
    f_n = compiled.param_values["biomass_N_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj, total_nitrogen(compiled.schema, biomass_nitrogen_fraction=f_n), label="nitrogen"
    )
    assert_nonnegative(traj, ("amino_acids", "debris", "X_dead"))
