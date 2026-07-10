"""Tests for the amino-acid ledger — the toggleable amino-acid swap (decision D-32).

Yeast build biomass mostly from amino acids, but the validated core sources all biomass
carbon from sugar and all biomass nitrogen from the lumped ammonium ``N`` pool. The
:class:`AminoAcidAssimilation` swap restores carbon/nitrogen honesty: dosed a ``default=0``
``amino_acids`` pool (represented as arginine), it funds a fraction of biomass from amino
acids by refunding sugar carbon and ammonium nitrogen and debiting the pool — a carbon- AND
nitrogen-neutral transfer, biomass untouched. It is nitrogen-anchored
(``ρ = ψ·gate(aa)·f_N·base_dx/y_N``) and, with an N-rich representative amino acid, its carbon
refund stays strictly below growth's sugar-carbon draw for any ψ ≤ 1, so it never creates
hexose.

**The correctness crux (decision D-32).** The safety uses growth's *pre-modifier* rate, but
growth's realised biomass is ``base_dx·M`` where ``M`` is the Arrhenius × (opt-in) carrying
modifiers ``ProcessSet`` applies. If the swap refunded at ``base_dx`` while growth drew at
``M·base_dx``, then at ``M < 1`` (cold ferment, or the carrying cap near saturation) the refund
could outrun the draw and CREATE sugar. The fix is that the wine growth Arrhenius and the
carrying-capacity modifier scale the swap too. The guard tests below are written to FAIL with an
unscaled swap (they exercise ``M < 1`` states through the full ``ProcessSet`` — at ``T_ref`` the
mismatch never fires) and to pass once the scaling lands.
"""

import numpy as np
import pytest

from fermentation.core.kinetics import (
    AminoAcidAssimilation,
    BiomassCarryingCapacity,
    GrowthNitrogenLimited,
    arrhenius_factor,
)
from fermentation.core.media import get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon, total_nitrogen

SWAP = AminoAcidAssimilation.name
GROWTH = GrowthNitrogenLimited.name


@pytest.fixture
def full_params():
    # The full default-wine parameter surface: every YAML the wired wine Processes read, so a
    # bare-built wine ProcessSet can evaluate its RHS (mirrors the carrying-capacity fixture).
    base = default_data_dir()
    return load_parameters(
        base / "wine_generic.yaml",
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "keto_acids.yaml",
        base / "hydrogen_sulfide.yaml",
        base / "aging.yaml",  # the bare wine set now carries the (default-off) aging Process (D-70)
    ).resolve()


def _wine_y0(
    schema: StateSchema,
    *,
    x: float = 1.0,
    s: float = 200.0,
    n: float = 0.2,
    aa: float = 0.0,
    t_k: float = 293.15,
) -> FloatArray:
    return schema.pack(
        {"X": x, "S": [s], "E": 40.0, "N": n, "T": t_k, "CO2": 5.0, "amino_acids": aa}
    )


def _isolate_growth_and_swap(full_params: dict[str, float]) -> ProcessSet:
    """A wine ProcessSet with every Process but growth + the swap disabled, so the ``S`` and
    ``N`` columns of the RHS reflect only growth's draw and the swap's refund. Modifiers stay
    as the wine medium wires them (the growth Arrhenius + carrying cap), which is the whole
    point: the guard tests read the *realised* scaling, not the raw derivatives."""
    ps = get_medium("wine").build_process_set()
    for p in ps.active:
        if p.name not in (GROWTH, SWAP):
            ps.disable(p.name)
    return ps


# -- metadata -----------------------------------------------------------------


def test_metadata():
    p = AminoAcidAssimilation()
    assert p.name == "amino_acid_assimilation"
    assert p.tier is Tier.SPECULATIVE
    # Refunds carbon to S and nitrogen to N; debits the aa pool. Does NOT touch X (growth
    # builds biomass; the swap only re-sources its atoms).
    assert set(p.touches) == {"amino_acids", "N", "S"}
    assert "X" not in p.touches
    assert "amino_acid_assimilation_fraction" in p.reads
    assert "K_amino_acids" in p.reads


# -- guard tests: the swap must never create sugar (decision D-32) -------------


def test_swap_never_creates_sugar_under_carrying_saturation(full_params):
    # Correctness guard (fail-first). With the carrying cap near saturation (X→K) the growth
    # draw is throttled to nearly zero while nitrogen is still available, so an UNSCALED swap
    # refunds sugar carbon faster than growth removed it → net sugar creation. Scaling the swap
    # by the same carrying factor closes it. Isolated to growth + swap so the S/N columns are
    # exactly draw + refund (T_ref, so Arrhenius = 1: this probes the CARRYING coupling).
    ps = _isolate_growth_and_swap(full_params)
    schema = ps.schema
    params = {**full_params, "biomass_carrying_capacity": 2.5}
    y = _wine_y0(schema, x=2.4, s=200.0, n=0.2, aa=5.0, t_k=293.15)  # X just below K=2.5
    d = ps.total_derivatives(0.0, y, params)
    assert d[schema.slice("S").start] <= 1e-12  # no net sugar creation
    assert d[schema.slice("N").start] <= 1e-12  # no net ammonium creation


def test_swap_refund_carries_the_growth_arrhenius_factor(full_params):
    # Correctness guard (fail-first) for the temperature coupling. With arginine's margin a
    # realistic cold ferment does not by itself create sugar, so the honest non-vacuous test is
    # that the swap's refund SCALES with growth's Arrhenius factor (else its aa-funded share of
    # biomass would wrongly drift with temperature, and at extreme cold outrun the draw). Isolate
    # the swap's own contribution by differencing enable/disable, with the carrying modifier off
    # so only the growth Arrhenius scales.
    ps = _isolate_growth_and_swap(full_params)
    ps.disable(BiomassCarryingCapacity.name)  # isolate the Arrhenius factor
    schema = ps.schema
    n_slot = schema.slice("N").start

    def swap_refund_n(t_k: float) -> float:
        y = _wine_y0(schema, x=1.0, s=200.0, n=0.2, aa=5.0, t_k=t_k)
        ps.enable(SWAP)
        d_both = ps.total_derivatives(0.0, y, full_params)
        ps.disable(SWAP)
        d_growth = ps.total_derivatives(0.0, y, full_params)
        ps.enable(SWAP)
        return float(d_both[n_slot] - d_growth[n_slot])  # the swap's N refund alone

    t_cold = 283.15
    ratio = swap_refund_n(t_cold) / swap_refund_n(293.15)
    expected = arrhenius_factor(t_cold, full_params["E_a_growth"], full_params["T_ref"])
    assert ratio == pytest.approx(expected, rel=1e-9)


# -- the swap is carbon- and nitrogen-neutral by construction -----------------


def test_swap_contribution_is_carbon_and_nitrogen_neutral(full_params):
    # The heart of the design: the swap only moves carbon aa → S and nitrogen aa → N, so its
    # OWN contribution changes neither total carbon nor total nitrogen — for any state, exactly
    # (biomass X is untouched). Isolate the swap's contribution (enable/disable difference) and
    # apply the conservation weight vectors to that derivative: both must be ~0.
    ps = _isolate_growth_and_swap(full_params)
    schema = ps.schema
    f_c = full_params["biomass_C_fraction"]
    f_n = full_params["biomass_N_fraction"]
    carbon = total_carbon(schema, biomass_carbon_fraction=f_c)
    nitrogen = total_nitrogen(schema, biomass_nitrogen_fraction=f_n)
    for x, s, n, aa in [(1.0, 200.0, 0.2, 5.0), (0.5, 240.0, 0.1, 1.0), (2.0, 80.0, 0.05, 0.3)]:
        y = _wine_y0(schema, x=x, s=s, n=n, aa=aa)
        ps.enable(SWAP)
        d_both = ps.total_derivatives(0.0, y, full_params)
        ps.disable(SWAP)
        d_growth = ps.total_derivatives(0.0, y, full_params)
        ps.enable(SWAP)
        d_swap = d_both - d_growth  # the swap's own (scaled) contribution
        assert abs(carbon(d_swap)) < 1e-12
        assert abs(nitrogen(d_swap)) < 1e-12


# -- isolability (undosed-only) ----------------------------------------------


def test_empty_pool_swap_is_a_noop(full_params):
    # THE byte-for-byte claim: with the pool empty the swap contributes nothing, so the full
    # wine RHS with the swap ENABLED equals the RHS with it DISABLED, exactly, on every column.
    # This also pins that folding the swap into the growth Arrhenius / carrying ``modifies``
    # tuples is transparent when the swap is inactive (a modifier naming a zero-contribution
    # Process is a no-op).
    on = get_medium("wine").build_process_set()
    off = get_medium("wine").build_process_set()
    off.disable(SWAP)
    schema = on.schema
    for x, s, n in [(0.5, 240.0, 0.15), (2.0, 120.0, 0.02), (3.0, 40.0, 0.0)]:
        y = _wine_y0(schema, x=x, s=s, n=n, aa=0.0)  # empty aa pool
        diff = on.total_derivatives(0.0, y, full_params) - off.total_derivatives(
            0.0, y, full_params
        )
        assert np.max(np.abs(diff)) == 0.0


def test_default_compile_disables_the_swap():
    # No amino_acids_gpl ⇒ the swap is present (wired into wine) but DISABLED, so an undosed run
    # is the validated core.
    _, compiled = _run(80.0)
    assert SWAP in compiled.process_set
    assert not compiled.process_set.is_enabled(SWAP)


def test_dose_enables_the_swap():
    _, compiled = _run(80.0, amino_acids_gpl=2.0)
    assert compiled.process_set.is_enabled(SWAP)
    assert compiled.y0[compiled.schema.slice("amino_acids").start] == pytest.approx(2.0)


# -- tier: structural drop only when enabled ---------------------------------


def test_swap_drops_s_and_n_output_tier_structurally_only_when_enabled():
    # Growth (PLAUSIBLE) touches S and N. The speculative swap also touches them, so enabling it
    # drops tier_of("S")/("N") to SPECULATIVE — but ONLY when enabled: a disabled Process is
    # excluded from tier derivation (the wine-only MLF/carrying *tier* isolability argument).
    # Undosed wine keeps them PLAUSIBLE. Minimal [growth, swap] set so growth is the only other
    # toucher (the full wine set's speculative byproducts already make S speculative regardless).
    schema = wine_schema()
    procs = [GrowthNitrogenLimited(), AminoAcidAssimilation()]
    off = ProcessSet(schema, procs)
    off.disable(SWAP)
    on = ProcessSet(schema, procs)
    for var in ("S", "N"):
        assert off.tier_of(var) is Tier.PLAUSIBLE
        assert on.tier_of(var) is Tier.SPECULATIVE


# -- behaviour through the compile seam --------------------------------------


def _run(
    yan_mgl: float,
    *,
    amino_acids_gpl: float | None = None,
    carrying_capacity_gpl: float | None = None,
    days: float = 14.0,
):
    """Compile + integrate a wine ferment at the given YAN; dose the amino-acid ledger when
    ``amino_acids_gpl`` is set and opt into the carrying cap when ``carrying_capacity_gpl`` is.
    Returns (trajectory, compiled)."""
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if amino_acids_gpl is not None:
        initial["amino_acids_gpl"] = amino_acids_gpl
    if carrying_capacity_gpl is not None:
        initial["carrying_capacity_gpl"] = carrying_capacity_gpl
    scenario = Scenario(
        name=f"wine-aa-{yan_mgl:.0f}-{amino_acids_gpl or 0:.1f}",
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


def test_carbon_and_nitrogen_close_with_the_swap_on():
    # The crown-jewel invariant with the aa pool weighted in BOTH ledgers: over a full dosed
    # run, carbon and nitrogen close to solver tolerance (the swap moves atoms between pools, it
    # neither creates nor destroys them).
    traj, compiled = _run(150.0, amino_acids_gpl=2.0)
    f_c = compiled.param_values["biomass_C_fraction"]
    f_n = compiled.param_values["biomass_N_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj, total_nitrogen(compiled.schema, biomass_nitrogen_fraction=f_n), label="nitrogen"
    )


def test_amino_acid_pool_depletes_and_stays_nonnegative():
    from fermentation.validation import assert_nonnegative

    traj, _ = _run(150.0, amino_acids_gpl=2.0)
    aa = np.asarray(traj.series("amino_acids"))
    assert aa[0] == pytest.approx(2.0)
    assert aa[-1] < aa[0]  # consumed
    assert_nonnegative(traj, ("amino_acids",))


def test_both_opt_ins_compose_through_the_compile_seam():
    # The end-to-end verification of the modifier fix's target regime: dose the amino-acid ledger
    # AND opt into the carrying cap (near-saturation with residual N is exactly where an unscaled
    # swap would create sugar). Through the real compile seam both modifiers end up scaling the
    # enabled swap, so carbon + nitrogen close over the whole run and sugar reaches dryness with
    # no creation — the integrated counterpart to the isolated-RHS guard tests.
    traj, compiled = _run(300.0, amino_acids_gpl=2.0, carrying_capacity_gpl=2.5)
    assert compiled.process_set.is_enabled(SWAP)
    assert compiled.process_set.is_enabled(BiomassCarryingCapacity.name)
    f_c = compiled.param_values["biomass_C_fraction"]
    f_n = compiled.param_values["biomass_N_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj, total_nitrogen(compiled.schema, biomass_nitrogen_fraction=f_n), label="nitrogen"
    )
    s = np.asarray(traj.series("S"))
    total_sugar = s if s.ndim == 1 else s.sum(axis=0)
    assert float(total_sugar[-1]) < 2.0  # reaches dryness, no sugar created along the way


def test_dose_behaves_like_supplementary_yan():
    # The emergent (second-order) effect: growth's derivatives are untouched, but the swap
    # refunds ammonium N, so the pool growth reads is replenished on the next step — dosing amino
    # acids acts like extra YAN, building more biomass. Directional, not a fit value.
    peak_undosed = float(np.max(_run(80.0)[0].series("X")))
    peak_dosed = float(np.max(_run(80.0, amino_acids_gpl=2.0)[0].series("X")))
    assert peak_dosed > peak_undosed * 1.05  # materially more biomass with the aa dose
