"""Tests for the fusel Ehrlich re-route — the deamination branch (decision D-33).

:class:`FuselAlcoholsEhrlich` books fusel carbon out of *sugar* (a stand-in — the real Ehrlich
pathway builds higher alcohols from amino-acid skeletons and releases the amino nitrogen). Once
the toggleable ``amino_acids`` pool (arginine; D-32) exists, :class:`FuselAminoAcidReroute`
re-sources a fraction ``g = aa/(K_amino_acids+aa)`` of that carbon off sugar and onto amino acids
and **deaminates** — releasing the consumed amino acids' nitrogen to the ammonium ``N`` pool. It
is a separate wine-only *swap* (never touches ``fusels``; production stays in the producer), and
the two share one :func:`fusel_production_rate` so the re-route's sugar refund matches the
producer's draw exactly. This suite pins the closure algebra, the deamination direction, that
production is untouched at the derivative level, and undosed isolability.
"""

import numpy as np
import pytest

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics import (
    FuselAlcoholsEhrlich,
    FuselAminoAcidReroute,
    fusel_carbon_draw,
)
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS
from fermentation.core.media import get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon, total_nitrogen

REROUTE = FuselAminoAcidReroute.name
PRODUCER = FuselAlcoholsEhrlich.name
# The five single-molecule higher-alcohol pools the lumped `fusels` pool became at D-99.
_FUSEL_POOLS = tuple(spec.pool for spec in FUSEL_SPECS)
_AA_SPECIES = "arginine"


@pytest.fixture
def full_params():
    base = default_data_dir()
    return load_parameters(
        base / "wine_generic.yaml",
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "keto_acids.yaml",
        base / "hydrogen_sulfide.yaml",
        base / "aging.yaml",  # the bare wine set now carries the (default-off) aging Process (D-70)
        base / "thermal.yaml",  # D-87/D-88: Caramelization fires at S>0 (must), so the bare
        # wine set needs its params (the first aging Process not inert at a fermentation state)
    ).resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x: float = 1.0,
    s: float = 200.0,
    n: float = 0.15,
    aa: float = 0.0,
    e: float = 40.0,
    t_k: float = 293.15,
) -> FloatArray:
    return schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": t_k, "CO2": 5.0, "amino_acids": aa})


def _isolate_fusel(full_params: dict[str, float]) -> ProcessSet:
    """A wine ProcessSet with only the fusel producer + re-route active, so the S/amino_acids/N
    columns reflect exactly the producer's draw and the re-route's refund/debit/deamination."""
    ps = get_medium("wine").build_process_set()
    for p in ps.active:
        if p.name not in (PRODUCER, REROUTE):
            ps.disable(p.name)
    return ps


# -- metadata -----------------------------------------------------------------


def test_metadata():
    p = FuselAminoAcidReroute()
    assert p.name == "fusel_amino_acid_reroute"
    assert p.tier is Tier.SPECULATIVE
    # A swap: it moves the carbon SOURCE and releases nitrogen — it never produces fusels.
    assert set(p.touches) == {"S", "amino_acids", "N"}
    for pool in _FUSEL_POOLS:
        assert pool not in p.touches
    # Reads all FIVE per-species k's since D-99 — it must reproduce the producer's TOTAL
    # draw, which is the sum over five rates at five different carbon fractions.
    for r in (*(spec.k_param for spec in FUSEL_SPECS), "E_a_fusels", "K_amino_acids"):
        assert r in p.reads


# -- the re-route's own contribution is carbon- and nitrogen-neutral ----------


def test_reroute_contribution_is_carbon_and_nitrogen_neutral(full_params):
    # The heart of the swap: the re-route only moves carbon aa → S and nitrogen aa → N, so its
    # OWN contribution changes neither total carbon nor total nitrogen (fusels untouched). Isolate
    # it by enable/disable differencing and apply the conservation weights to that derivative.
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    f_c = full_params["biomass_C_fraction"]
    f_n = full_params["biomass_N_fraction"]
    carbon = total_carbon(schema, biomass_carbon_fraction=f_c)
    nitrogen = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)
    for x, s, n, aa in [(1.0, 200.0, 0.15, 5.0), (0.5, 240.0, 0.1, 1.0), (2.0, 80.0, 0.05, 0.3)]:
        y = _wine_y0(schema, x=x, s=s, n=n, aa=aa)
        ps.enable(REROUTE)
        d_both = ps.total_derivatives(0.0, y, full_params)
        ps.disable(REROUTE)
        d_prod = ps.total_derivatives(0.0, y, full_params)
        ps.enable(REROUTE)
        d_reroute = d_both - d_prod
        assert abs(carbon(d_reroute)) < 1e-14
        assert abs(nitrogen(d_reroute)) < 1e-14


# -- production is untouched; only the carbon source moves --------------------


def test_reroute_never_touches_any_fusel_pool_or_other_columns(full_params):
    # Derivative-level guard (the warm=more-fusel benchmark is protected): enabling the re-route
    # changes ONLY S/amino_acids/N — never the five fusel pools/E/CO2/X (D-99). Production
    # stays in the producer.
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    y = _wine_y0(schema, x=1.0, s=200.0, n=0.15, aa=5.0)
    ps.enable(REROUTE)
    d_both = ps.total_derivatives(0.0, y, full_params)
    ps.disable(REROUTE)
    d_prod = ps.total_derivatives(0.0, y, full_params)
    ps.enable(REROUTE)
    d_reroute = d_both - d_prod
    for col in (*_FUSEL_POOLS, "E", "CO2", "X"):
        assert d_reroute[schema.slice(col)][0] == 0.0
    # It DOES move exactly these three:
    assert d_reroute[schema.slice("amino_acids")][0] < 0.0  # debited
    assert d_reroute[schema.slice("N")][0] > 0.0  # deaminated (ammonium released)


def test_reroute_matches_the_producer_draw_exactly(full_params):
    # The single-helper guarantee, made concrete: the re-route refunds sugar carbon equal to
    # ``g`` × the carbon the producer drew, so the amino-acid debit (as carbon) equals the sugar
    # refund (as carbon) — i.e. exactly the fraction of the producer's draw is re-sourced.
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    aa_val = 0.3
    y = _wine_y0(schema, x=1.5, s=180.0, n=0.12, aa=aa_val)
    g = aa_val / (full_params["K_amino_acids"] + aa_val)
    # The producer's TOTAL draw across all five species, each at its OWN carbon fraction
    # (D-99). Before the split this was one rate times one stand-in fraction.
    fusel_carbon = fusel_carbon_draw(y, schema, full_params)

    ps.enable(REROUTE)
    d_both = ps.total_derivatives(0.0, y, full_params)
    ps.disable(REROUTE)
    d_prod = ps.total_derivatives(0.0, y, full_params)
    ps.enable(REROUTE)
    d_reroute = d_both - d_prod

    # sugar carbon refunded by the re-route (single wine slot):
    sugar_refund_c = d_reroute[schema.slice("S")][0] * carbon_mass_fraction("glucose")
    assert sugar_refund_c == pytest.approx(g * fusel_carbon, rel=1e-12)
    # amino-acid carbon debited equals that same carbon:
    aa_debit_c = -d_reroute[schema.slice("amino_acids")][0] * carbon_mass_fraction(_AA_SPECIES)
    assert aa_debit_c == pytest.approx(g * fusel_carbon, rel=1e-12)
    # and the deamination releases exactly that amino-acid mass's nitrogen:
    aa_mass = g * fusel_carbon / carbon_mass_fraction(_AA_SPECIES)
    assert d_reroute[schema.slice("N")][0] == pytest.approx(
        aa_mass * nitrogen_mass_fraction(_AA_SPECIES), rel=1e-12
    )


def test_reroute_never_creates_sugar(full_params):
    # Net sugar (producer draw + re-route refund) is −(1−g)·F_c ≤ 0 for any g ≤ 1: the re-route
    # only spares sugar, never creates it. Also the producer+re-route net draw is *smaller* than
    # the producer's alone (some carbon now comes from amino acids).
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    for aa in (0.1, 1.0, 10.0):
        y = _wine_y0(schema, x=1.5, s=180.0, n=0.12, aa=aa)
        ps.enable(REROUTE)
        d_both = ps.total_derivatives(0.0, y, full_params)
        ps.disable(REROUTE)
        d_prod = ps.total_derivatives(0.0, y, full_params)
        ps.enable(REROUTE)
        assert d_both[schema.slice("S")][0] <= 0.0  # never creates sugar
        assert (
            d_both[schema.slice("S")][0] > d_prod[schema.slice("S")][0]
        )  # smaller draw than producer


# -- isolability (undosed-only) ----------------------------------------------


def test_empty_pool_reroute_is_a_noop(full_params):
    # With the pool empty the availability gate → 0, so the full wine RHS with the re-route ENABLED
    # equals the RHS with it DISABLED, exactly, on every column (byte-for-byte the sugar stand-in).
    on = get_medium("wine").build_process_set()
    off = get_medium("wine").build_process_set()
    off.disable(REROUTE)
    schema = on.schema
    for x, s, n in [(0.5, 240.0, 0.15), (2.0, 120.0, 0.05), (1.0, 200.0, 0.12)]:
        y = _wine_y0(schema, x=x, s=s, n=n, aa=0.0)
        diff = on.total_derivatives(0.0, y, full_params) - off.total_derivatives(
            0.0, y, full_params
        )
        assert np.max(np.abs(diff)) == 0.0


def test_compile_seam_toggles_the_reroute():
    _, undosed = _run(150.0)
    assert REROUTE in undosed.process_set
    assert not undosed.process_set.is_enabled(REROUTE)
    _, dosed = _run(150.0, amino_acids_gpl=2.0)
    assert dosed.process_set.is_enabled(REROUTE)


# -- tier: structural drop only when enabled ---------------------------------


def test_reroute_drops_s_and_n_output_tier_only_when_enabled():
    # The speculative re-route touches S and N; enabling it drops tier_of("S")/("N") to
    # SPECULATIVE, but only when enabled (a disabled Process is excluded from tier derivation).
    schema = wine_schema()
    procs = [FuselAlcoholsEhrlich(), FuselAminoAcidReroute()]
    off = ProcessSet(schema, procs)
    off.disable(REROUTE)
    on = ProcessSet(schema, procs)
    # The producer is itself speculative and touches S, so S is speculative either way; N is only
    # touched by the re-route here, so it is the clean witness of the enabled-only drop.
    assert off.tier_of("N") is Tier.VALIDATED  # nothing enabled touches N
    assert on.tier_of("N") is Tier.SPECULATIVE


# -- behaviour through the compile seam --------------------------------------


def _run(yan_mgl: float, *, amino_acids_gpl: float | None = None, days: float = 14.0):
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if amino_acids_gpl is not None:
        initial["amino_acids_gpl"] = amino_acids_gpl
    scenario = Scenario(
        name=f"wine-fusel-reroute-{yan_mgl:.0f}",
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


def test_carbon_and_nitrogen_close_with_the_reroute_on():
    # Crown jewel: over a full dosed run the fusels pool is fed partly from amino acids (the
    # re-route) and partly from sugar (the producer), and the deamination refunds N — carbon and
    # nitrogen both close to solver tolerance (atoms only move between weighted pools).
    traj, compiled = _run(150.0, amino_acids_gpl=2.0)
    f_c = compiled.param_values["biomass_C_fraction"]
    f_n = compiled.param_values["biomass_N_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj, total_nitrogen(compiled.schema, biomass_nitrogen_fraction=f_n), label="nitrogen"
    )
