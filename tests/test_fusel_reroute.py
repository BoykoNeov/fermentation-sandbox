"""Tests for the fusel Ehrlich re-route — the deamination branch (decision D-33).

:class:`FuselAlcoholsEhrlich` books fusel carbon out of *sugar* (a stand-in — the real Ehrlich
pathway builds higher alcohols from amino-acid skeletons and releases the amino nitrogen). Once
the amino-acid pools exist (D-32, SPECIATED at D-100), :class:`FuselAminoAcidReroute` re-sources
a fraction ``g_i = aa_i/(K_amino_acids·f_i + aa_i)`` of **each alcohol's** carbon off sugar and onto
**that alcohol's own precursor** — leucine for isoamyl alcohol, valine for isobutanol, and so on —
and **deaminates**, releasing the consumed amino acids' nitrogen to the ammonium ``N`` pool. It
is a separate wine-only *swap* (never touches ``fusels``; production stays in the producer), and
the two share one :func:`fusel_production_rate` so the re-route's sugar refund matches the
producer's draw exactly. This suite pins the closure algebra, the deamination direction, that
production is untouched at the derivative level, and undosed isolability.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.chemistry import CARBON_ATOMS, carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics import (
    FuselAlcoholsEhrlich,
    FuselAminoAcidReroute,
    fusel_carbon_draw,
)
from fermentation.core.kinetics.byproducts import ehrlich_draws, fusel_carbon_draw_by_species
from fermentation.core.kinetics.carbon_routing import DE_NOVO_FUSEL_ROUTES, FUSEL_SPECS
from fermentation.core.media import get_medium, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon, total_nitrogen
from tests.conftest import seed_amino_acids

REROUTE = FuselAminoAcidReroute.name
PRODUCER = FuselAlcoholsEhrlich.name
# The five single-molecule higher-alcohol pools the lumped `fusels` pool became at D-99.
_FUSEL_POOLS = tuple(spec.pool for spec in FUSEL_SPECS)
#: Each alcohol's precursor amino acid — the pools the re-route actually debits since D-100.
#: `amino_acids` (arginine) is NOT among them, which is the whole decoupling.
_PRECURSORS = tuple(spec.precursor_amino_acid for spec in FUSEL_SPECS)


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
    params: Mapping[str, float],
    *,
    x: float = 1.0,
    s: float = 200.0,
    n: float = 0.15,
    aa: float = 0.0,
    e: float = 40.0,
    t_k: float = 293.15,
) -> FloatArray:
    # Amino acids seeded at MUST-SPECTRUM composition (D-100) — the state in which every
    # per-precursor gate aa_i/(K·f_i + aa_i) provably equals the pre-split lumped gate aa/(K + aa),
    # so the closed forms below assert the same numbers the lumped suite did.
    y = schema.pack({"X": x, "S": [s], "E": e, "N": n, "T": t_k, "CO2": 5.0})
    return seed_amino_acids(y, schema, params, aa)


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
    # ``CO2`` joined at D-106: the Ehrlich decarboxylation the re-sourced fraction really performs.
    # ``isoamyl_alcohol_valine`` joined at D-115: this is where a gram of isoamyl alcohol first
    # becomes valine-derived, so it is where the label enters the model. It is a SUB-QUANTITY of
    # the alcohol pool carrying zero carbon weight — so the "never produces fusels" claim below
    # still holds exactly, and is still checked against the real pools.
    assert set(p.touches) == {"S", "N", "CO2", *_PRECURSORS, "isoamyl_alcohol_valine"}
    for pool in _FUSEL_POOLS:
        assert pool not in p.touches
    # THE D-100 DECOUPLING, pinned: the re-route no longer touches the identity-agnostic pools.
    # Arginine does not make higher alcohols, so fusel production can no longer drain the substrate
    # the yeast swap / MLF growth / Brett growth / Maillard browning live on. This single assertion
    # is what makes D-99's ~3.8x fusel rise unable to starve three unrelated subsystems again.
    assert "amino_acids" not in p.touches
    assert "amino_acids_generic" not in p.touches
    # Reads all FIVE per-species k's since D-99 — it must reproduce the producer's TOTAL
    # draw, which is the sum over five rates at five different carbon fractions — plus each
    # precursor's must-spectrum share, which scales its own gate (D-100).
    for r in (*(spec.k_param for spec in FUSEL_SPECS), "E_a_fusels", "K_amino_acids"):
        assert r in p.reads
    for precursor in _PRECURSORS:
        assert f"must_aa_fraction_{precursor}" in p.reads


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
        y = _wine_y0(schema, full_params, x=x, s=s, n=n, aa=aa)
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
    # changes ONLY S/precursors/N/CO2 — never the five fusel pools/E/X (D-99). Production
    # stays in the producer. ``CO2`` left this exclusion list at D-106 and is asserted POSITIVE
    # below instead: the re-sourced fraction is a real decarboxylation, so a zero there is now the
    # bug rather than the contract.
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    y = _wine_y0(schema, full_params, x=1.0, s=200.0, n=0.15, aa=5.0)
    ps.enable(REROUTE)
    d_both = ps.total_derivatives(0.0, y, full_params)
    ps.disable(REROUTE)
    d_prod = ps.total_derivatives(0.0, y, full_params)
    ps.enable(REROUTE)
    d_reroute = d_both - d_prod
    for col in (*_FUSEL_POOLS, "E", "X"):
        assert d_reroute[schema.slice(col)][0] == 0.0
    # It DOES debit every precursor, deaminate, and decarboxylate:
    for precursor in _PRECURSORS:
        assert d_reroute[schema.slice(precursor)][0] < 0.0  # debited
    assert d_reroute[schema.slice("N")][0] > 0.0  # deaminated (ammonium released)
    assert d_reroute[schema.slice("CO2")][0] > 0.0  # decarboxylated (D-106)
    # ...and it leaves the identity-agnostic pools ALONE (the D-100 decoupling, at the
    # derivative level rather than the metadata level).
    assert d_reroute[schema.slice("amino_acids")][0] == 0.0
    assert d_reroute[schema.slice("amino_acids_generic")][0] == 0.0


def test_reroute_matches_the_producer_draw_exactly(full_params):
    # The single-helper guarantee, made concrete: the re-route refunds sugar carbon equal to
    # ``g`` × the carbon the producer drew — i.e. exactly the fraction of the producer's draw is
    # re-sourced, and the REFUND is the invariant the two Processes share.
    #
    # **The debit is no longer that same number (decision D-106).** Until D-106 the precursor debit
    # equalled the sugar refund exactly, because the draw was sized to the alcohol's carbon alone —
    # which is precisely what made it (n-1)/n of a mole instead of 1. The draw now also carries the
    # decarboxylation CO2, so the identity that holds is the three-way one:
    #
    #     precursor carbon OUT == sugar carbon REFUNDED + CO2 carbon EMITTED
    #
    # **And the refund is no longer ``g × the producer's draw`` either (decision D-111).** That
    # identity WAS the one-alcohol-one-precursor world: every branch refunded the whole alcohol
    # carbon it re-sourced, so the total was just the gated fraction of what the producer drew.
    # Isoamyl alcohol now has a SECOND precursor — valine via alpha-ketoisocaproate — which
    # re-sources FURTHER isoamyl carbon off sugar, so the refund is strictly larger. The extra
    # term is derived below from the SOURCED shares and the route's stoichiometry, never from the
    # code's own output: fitting the expectation to the implementation would assert nothing (the
    # D-96/D-102/D-109 "the sentence and the assertion are not the same claim" trap).
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    # 0.05 rather than the 0.3 this test used pre-D-111, and the reason is the finding rather than
    # convenience: at 0.3 the KIC branch's headroom clamp BINDS (leucine's gate already claims 75%
    # of isoamyl's carbon and valine wants another 26%), so the closed form below — which is the
    # UNCLAMPED algebra — is not what the code should compute there. This test asserts the sourcing
    # algebra, so it is run in the regime where the algebra applies; the clamp is pinned on its own
    # in `test_the_kic_branch_cannot_source_more_isoamyl_than_is_being_made`, at a dose where it
    # really does bind. Splitting them keeps each assertion able to fail for exactly one reason.
    aa_val = 0.05
    y = _wine_y0(schema, full_params, x=1.5, s=180.0, n=0.12, aa=aa_val)
    # AT MUST-SPECTRUM COMPOSITION every per-precursor gate collapses to this one lumped value
    # (the D-100 reduction property), so the aggregate below is exact, not approximate.
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
    # The D-111 valine -> KIC -> isoamyl branch's extra refund, from the parameter file's shares
    # and the route's stoichiometry ALONE:
    #   consumed valine splits share_but : share_iso : f  (isobutanol : isoamyl : everything else),
    #   the isobutanol branch is the anchor and draws (n+1)/n of its alcohol carbon (D-106),
    #   and the KIC branch refunds (n_p - k)/n_a = (5-2)/5 of the isoamyl carbon it re-sources
    #   -- valine C5 + 2 C acetyl-CoA -> isoamyl C5 + 2 CO2, of which the producer already drew
    #   all 5 C from sugar while the truth needs only 2, so 3 come back.
    f_val = full_params["f_non_ehrlich_valine"]
    share_iso = full_params["f_valine_to_isoamyl"]
    share_but = 1.0 - f_val - share_iso
    but_carbon = next(
        c
        for spec, c in fusel_carbon_draw_by_species(y, schema, full_params)
        if spec.pool == "isobutanol"
    )
    n_but = CARBON_ATOMS["isobutanol"]
    consumed_valine_c = (g * but_carbon * (n_but + 1) / n_but) / share_but
    kic_isoamyl_c = share_iso * consumed_valine_c  # valine C5 -> isoamyl C5, so 1:1 in carbon
    extra_refund = kic_isoamyl_c * 3.0 / 5.0
    assert extra_refund > 0.0, "vacuous: the D-111 KIC branch contributed nothing"
    # D-118: a de-novo-dominated alcohol's primary branch is capped at (1 - f_de_novo), so the
    # producer's draw for the de-novo share is deliberately NOT refunded — that carbon genuinely
    # did come from sugar (the shikimate route to phenylpyruvate), which is the whole content of
    # the route. Derived from the registry rather than hard-coded, so a second de-novo entry
    # cannot silently invalidate this closed form.
    de_novo_forgone = 0.0
    for route in DE_NOVO_FUSEL_ROUTES:
        alcohol_c = next(
            c
            for spec, c in fusel_carbon_draw_by_species(y, schema, full_params)
            if spec.pool == route.alcohol_pool
        )
        de_novo_forgone += g * full_params[route.share_param] * alcohol_c
    assert de_novo_forgone > 0.0, "vacuous: the D-118 de-novo cap withheld nothing"
    assert sugar_refund_c == pytest.approx(
        g * fusel_carbon - de_novo_forgone + extra_refund, rel=1e-12
    )
    # The carbon debited ACROSS THE FIVE PRECURSORS — each at its OWN carbon fraction. Before
    # D-100 this was one debit at arginine's fraction; the sum is what must now match, because each
    # alcohol eats a different molecule. Since D-106 it matches refund + CO2, not refund alone.
    aa_debit_c = sum(
        -d_reroute[schema.slice(precursor)][0] * carbon_mass_fraction(precursor)
        for precursor in _PRECURSORS
    )
    co2_c = d_reroute[schema.slice("CO2")][0] * carbon_mass_fraction("CO2")
    # The same D-111 term appears here too, and it MUST: the three-way identity is over whatever
    # the precursors actually gave up, and valine now gives up a second molecule for the KIC route.
    # Its sugar refund is 3/5 of that branch (not the whole of it, as on the primary routes), which
    # is why the extra term is `extra_refund` and not `kic_isoamyl_c` — the missing 2/5 is the two
    # CO2, already inside `co2_c`.
    # D-118's cap enters here for the same reason it entered the refund identity above: the
    # de-novo share is never sourced from phenylalanine, so the precursors never give it up.
    assert aa_debit_c == pytest.approx(
        g * fusel_carbon - de_novo_forgone + extra_refund + co2_c, rel=1e-12
    )
    # The CO2 is a real share of the draw, not a rounding term: one carbon per alcohol means the
    # precursors give up ~1/n MORE carbon than they did before D-106. Pin that it is neither zero
    # nor the whole draw, so a silently-dropped term cannot pass as "approximately equal".
    assert 0.1 * aa_debit_c < co2_c < 0.35 * aa_debit_c
    # ...and the deamination releases exactly the nitrogen THOSE molecules carried. This is the
    # D-33 over-release lump being retired, not restated: arginine (4 N over 6 C) would have
    # released ~4x the real leucine->isoamyl N:C, and now each precursor releases its own.
    expected_n = sum(
        -d_reroute[schema.slice(precursor)][0] * nitrogen_mass_fraction(precursor)
        for precursor in _PRECURSORS
    )
    assert d_reroute[schema.slice("N")][0] == pytest.approx(expected_n, rel=1e-12)


def test_reroute_never_creates_sugar(full_params):
    # Net sugar (producer draw + re-route refund) is −(1−g)·F_c ≤ 0 for any g ≤ 1: the re-route
    # only spares sugar, never creates it. Also the producer+re-route net draw is *smaller* than
    # the producer's alone (some carbon now comes from amino acids).
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    for aa in (0.1, 1.0, 10.0):
        y = _wine_y0(schema, full_params, x=1.5, s=180.0, n=0.12, aa=aa)
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
        y = _wine_y0(schema, full_params, x=x, s=s, n=n, aa=0.0)
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


# -- D-111: the many-to-one map (valine -> KIC -> isoamyl alcohol) ------------


def test_valine_feeds_isoamyl_alcohol_as_well_as_isobutanol(full_params):
    """The D-111 route exists at all: valine sources TWO alcohols, not one.

    D-104 named this route as the model's missing one ("the gap is a missing route (valine -> KIC
    -> isoamyl)") and declined to build it; D-109 scoped it as the milestone's real content.
    Before D-111 ``FuselSpec`` pinned one alcohol to one precursor, so this could not be expressed.
    """
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    y = _wine_y0(schema, full_params, x=1.5, s=180.0, n=0.12, aa=0.05)
    by_alcohol = {
        d.alcohol.pool: d
        for d in ehrlich_draws(y, schema, full_params)
        if d.precursor.species == "valine"
    }
    assert set(by_alcohol) == {"isobutanol", "isoamyl_alcohol"}
    # The KIC route's stoichiometry is what makes it a different KIND of branch, not just another
    # one: valine C5 + 2 C acetyl-CoA -> isoamyl C5 + 2 CO2, against the primary route's
    # C(precursor) == C(alcohol) + 1 and ONE CO2. Both are asserted, because the D-106 helper's
    # "one constant covers the set" is exactly what this route retires.
    kic = by_alcohol["isoamyl_alcohol"]
    assert kic.co2_carbon == pytest.approx(kic.alcohol_carbon * 2.0 / 5.0, rel=1e-12)
    assert kic.precursor_carbon == pytest.approx(kic.alcohol_carbon, rel=1e-12)  # C5 -> C5
    assert kic.refund_carbon == pytest.approx(kic.alcohol_carbon * 3.0 / 5.0, rel=1e-12)
    primary = by_alcohol["isobutanol"]
    assert primary.co2_carbon == pytest.approx(primary.alcohol_carbon / 4.0, rel=1e-12)


def test_every_ehrlich_branch_is_carbon_neutral_in_the_sourcing_layer(full_params):
    """Each branch only moves carbon's SOURCE: refund + CO2 == precursor drawn, exactly.

    The invariant that lets D-111 add a route with two CO2 and an acetyl-CoA co-substrate without
    touching ``total_carbon``: the producer already made the alcohol and already drew sugar for it,
    so a branch that credited the alcohol would double-count. Holds per branch, so no branch can
    be wrong in a way another branch's error hides.
    """
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    for aa in (0.01, 0.05, 0.3, 1.0):
        y = _wine_y0(schema, full_params, x=1.5, s=180.0, n=0.12, aa=aa)
        draws = ehrlich_draws(y, schema, full_params)
        assert draws, "vacuous: no branches at all"
        for d in draws:
            assert d.refund_carbon + d.co2_carbon == pytest.approx(d.precursor_carbon, rel=1e-12)


def test_the_kic_branch_cannot_source_more_isoamyl_than_is_being_made(full_params):
    """The headroom clamp, and the D-103 gate-shape defect it makes VISIBLE.

    The KIC branch is anchored off ISOBUTANOL's draw (the share is a fraction of consumed valine),
    so nothing in its arithmetic bounds it by isoamyl's own production — unlike a gated primary
    branch, which is a fraction of exactly that. At a realistic must (aa = 1.0 g/L => leucine
    ~32 mg/L, valine ~37 mg/L) leucine's gate already claims **90.9%** of isoamyl's instantaneous
    carbon and the KIC branch wants a further **31.8%**: two independently *sourced* claims summing
    to **122.7%** of one alcohol.

    **That sum is D-103's finding becoming visible rather than a new defect.** D-103 measured the
    gate's shape as the real problem ("which no scalar can fix") — the model sources far too much
    of each alcohol from its precursor early, and only survives integrated because the precursor
    pools are tiny and empty fast. One sourcing claim can be that wrong silently; TWO cannot, and
    the many-to-one map is what turns the error into an arithmetic impossibility.

    Without the clamp the refund would hand back more sugar than the producer ever drew for this
    alcohol, and **conservation would not notice** — the ledger closes either way (the D-89/D-90
    denominator-trap family, where only an explicit guard is not blind).
    """
    ps = _isolate_fusel(full_params)
    schema = ps.schema
    y = _wine_y0(schema, full_params, x=1.5, s=180.0, n=0.12, aa=1.0)
    totals = dict(fusel_carbon_draw_by_species(y, schema, full_params))
    f_iso = next(c for spec, c in totals.items() if spec.pool == "isoamyl_alcohol")
    by_precursor = {
        d.precursor.species: d
        for d in ehrlich_draws(y, schema, full_params)
        if d.alcohol.pool == "isoamyl_alcohol"
    }
    sourced = sum(d.alcohol_carbon for d in by_precursor.values())
    # The clamp is real at this dose: it must actually be binding, or this test is vacuous and
    # would pass on a model that never needed a guard.
    unclamped = full_params["f_valine_to_isoamyl"] * (
        next(
            d.precursor_carbon
            for d in ehrlich_draws(y, schema, full_params)
            if d.alcohol.pool == "isobutanol"
        )
        / (1.0 - full_params["f_non_ehrlich_valine"] - full_params["f_valine_to_isoamyl"])
    )
    assert unclamped + by_precursor["leucine"].alcohol_carbon > f_iso, (
        "vacuous: the clamp is not binding here, so this asserts nothing"
    )
    # ...and it holds: isoamyl is never sourced from more precursor than it is made from.
    assert sourced <= f_iso * (1.0 + 1e-12)


# -- D-115: the label reaches the ESTER, and what it is measured against ------


#: Rollero *et al.* 2017's own SM250 condition - 250 mg N/L, 24 C. The enrichments below are
#: quoted against that paper, so the probe must run on a COMMENSURATE must (D-104's lesson,
#: where a ~470 mg N/L probe flattered the model against a 180 mg N/L source).
_ROLLERO_YAN = 250.0
_ROLLERO_TEMP = 24.0

#: What Rollero measures, and the ONLY thing this beat may be scored against.
#:
#: The **alcohol** band is the model's existing, already-validated D-111 target (2.1-7.5%). The
#: **ester** figure is D-114's corrected reading: ~4% at the well-measured SM250/SM425 columns,
#: **not** the bare 0.0-19.7% full spread. That correction is load-bearing in the sceptical
#: direction - the four high values all sit in the low-nitrogen SM70 blocks and breach the
#: physical bound that an acetate's enrichment cannot exceed its parent alcohol's (by 3.1-5.3x),
#: so they are excluded on chemistry rather than on a noise argument. Quoted loosely as
#: "0-19.7%" the band would swallow almost any model output, including D-114's structural zero.
_ROLLERO_ALCOHOL_BAND = (0.021, 0.075)
_ROLLERO_ESTER = 0.04


def _rollero_run(amino_acids_gpl: float, *, days: float = 14.0):
    scenario = Scenario(
        name=f"d115-enrichment-{amino_acids_gpl}",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": _ROLLERO_YAN,
            "pitch_gpl": 0.25,
            "amino_acids_gpl": amino_acids_gpl,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=_ROLLERO_TEMP)],
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


def _enrichment(traj, schema, bulk: str, tracer: str) -> float:
    b = float(traj.y[schema.slice(bulk), -1][0])
    assert b > 0.0, f"vacuous: {bulk} is empty, so its enrichment is a division of zeros"
    return float(traj.y[schema.slice(tracer), -1][0]) / b


def test_the_ester_carries_valine_label_at_its_parent_alcohols_enrichment():
    """D-115's headline: the ester's enrichment is a model output, and it TRACKS THE ALCOHOL.

    **What is being claimed, and what is deliberately NOT.** Rollero measures isoamyl acetate at
    ~4% valine-labelled against its parent alcohol's 3.4-7.5% - the ester sits *at or marginally
    below* its alcohol throughout, which is what the chemistry predicts (the C5 skeleton
    transfers as a unit, and enrichment is a molecule fraction). **That RATIO is the structural
    claim this beat delivers.** The model's absolute enrichment is ~1.9%, roughly half the
    measured ~4%, and this test asserts *only the same order of magnitude* - because the
    shortfall is **inherited**, not new: the alcohol itself runs at ~1.8% against a measured
    2.1-7.5%, an already-documented D-111 gap.

    **Why the DoD is not "hit 4%".** Closing the residual gap here would mean raising
    ``f_valine_to_isoamyl``, which is independently sourced twice (Crepin 2017's 23%, sitting
    inside Rollero's own 42-45% KIC flux) and whose denominator is *consumed valine* while this
    test's is *total ester*. Tuning it to land this number is precisely the D-104 error - "the
    gap is a missing route ... NOT a value to tune" - with the route now built and the
    temptation moved one compound downstream. The parameter is untouched by this beat.
    """
    traj, compiled = _rollero_run(1.0)
    schema = compiled.schema
    alcohol = _enrichment(traj, schema, "isoamyl_alcohol", "isoamyl_alcohol_valine")
    ester = _enrichment(traj, schema, "isoamyl_acetate", "isoamyl_acetate_valine")

    # Anti-vacuity FIRST: a structural zero would satisfy "same order of magnitude as zero".
    assert ester > 0.0, (
        "the ester carries NO label - the D-115 re-route is not delivering, and this is exactly "
        "the structural zero D-114 measured"
    )
    assert alcohol > 0.0

    # THE STRUCTURAL CLAIM: the ester tracks its parent alcohol. Rollero's ester/alcohol ratio
    # is ~1 (at or marginally below); the model's is marginally ABOVE 1, because the alcohol
    # pool's fraction FALLS over the run and the ester integrates some of the earlier, richer
    # pool. That direction is a real, reportable difference from the source, not a fit.
    assert ester / alcohol == pytest.approx(1.0, abs=0.15), (
        "the ester must carry label at essentially its parent alcohol's enrichment - the C5 "
        "transfers as a unit, so the two fractions cannot diverge"
    )

    # ORDER OF MAGNITUDE ONLY against the source, and the band is deliberately generous in the
    # direction of the known shortfall. Tightening this toward Rollero's ~4% would convert an
    # honest inherited gap into pressure on a sourced parameter.
    assert 0.25 * _ROLLERO_ESTER < ester < 2.0 * _ROLLERO_ESTER, (
        f"model ester enrichment {ester:.4f} is not even the right order against Rollero's "
        f"~{_ROLLERO_ESTER:.2f} - something structural, not a calibration drift"
    )
    # And the alcohol's own gap is UNCHANGED by this beat: still just under D-111's band.
    lo, _hi = _ROLLERO_ALCOHOL_BAND
    assert 0.5 * lo < alcohol < lo, (
        f"alcohol enrichment {alcohol:.4f} moved out of its known just-under-band position - "
        "the re-route was supposed to cost the alcohol pool almost nothing (mass-negligible)"
    )


def test_no_amino_acid_dose_means_no_label_anywhere():
    """Isolability, and the anti-vacuity floor for every enrichment number above (D-115).

    With no amino-acid dose there is no valine, so no branch of ``ehrlich_draws`` sources
    isoamyl alcohol from it and **both** tracers must stay identically zero - while the pools
    themselves fill normally. Without this, an implementation that credited the tracer from the
    *leucine* branch (the same pool's primary precursor) would produce healthy-looking
    enrichments that no source supports, and nothing else in the suite would notice: a tracer
    slot carries no carbon weight, so conservation is blind to it by construction.
    """
    traj, compiled = _rollero_run(0.0)
    schema = compiled.schema

    # The pools are real - this is a working ferment, not an empty one.
    assert float(traj.y[schema.slice("isoamyl_alcohol"), -1][0]) > 0.0
    assert float(traj.y[schema.slice("isoamyl_acetate"), -1][0]) > 0.0
    # ...and carry no label whatsoever.
    assert float(traj.y[schema.slice("isoamyl_alcohol_valine"), -1][0]) == 0.0
    assert float(traj.y[schema.slice("isoamyl_acetate_valine"), -1][0]) == 0.0


def test_the_ester_needed_its_own_tracer_slot_because_the_alcohol_fraction_is_not_flat():
    """THE RECEIPT for the two-slot design (D-115) - a measurement, not an argument.

    The cheap design was one slot: track the alcohol's valine fraction and let the ester inherit
    it, justified by D-114's observation that ~93% of the ester forms after valine is exhausted,
    so the fraction should be *frozen* by then. **That premise is false, and this test is where
    it is falsified.** The alcohol's fraction is strongly time-varying across the window in which
    the ester forms - it starts high (a tiny young pool is dominated by the valine route) and
    falls by more than an order of magnitude as de-novo synthesis takes over.

    So the ester's enrichment is a flux-weighted average of a moving fraction, and inheriting a
    single number gets it wrong in a measurable way. Reconstructing that average after the fact
    is quadrature over interpolated states - the D-103 defect that overstated a draw 1.3-3.5x -
    which is why the ester carries its own integrated slot instead.

    **A tripwire, not a curiosity.** If a future beat flattens the alcohol's fraction (a
    different gate shape, a de-novo ceiling), the second slot stops earning its keep and this
    test fails - at which point the design question genuinely re-opens.
    """
    traj, compiled = _rollero_run(1.0)
    schema = compiled.schema
    alcohol = traj.y[schema.slice("isoamyl_alcohol")][0]
    alcohol_label = traj.y[schema.slice("isoamyl_alcohol_valine")][0]
    ester = traj.y[schema.slice("isoamyl_acetate")][0]

    made = np.clip(np.diff(ester, prepend=ester[0]), 0.0, None)
    assert made.sum() > 0.0, "vacuous: no ester formed over the run"
    fraction = np.where(alcohol > 0.0, alcohol_label / np.where(alcohol > 0.0, alcohol, 1.0), 0.0)

    forming = fraction[made > 0.0]
    assert forming.max() > 5.0 * forming.min(), (
        "the alcohol's valine fraction is nearly FLAT across the ester-forming window, so the "
        "ester could inherit it and the second tracer slot is redundant - re-open the D-115 "
        "design decision rather than carrying a slot that buys nothing"
    )

    # And the error the one-slot design would actually have made: the flux-weighted alcohol
    # fraction is NOT the ester's own integrated enrichment.
    inherited = float((fraction * made).sum() / made.sum())
    actual = _enrichment(traj, schema, "isoamyl_acetate", "isoamyl_acetate_valine")
    assert abs(inherited - actual) / actual > 0.05, (
        f"one-slot inheritance would have reported {inherited:.5f} against the tracer's "
        f"{actual:.5f} - if that error is now negligible the second slot is not paying for "
        "itself"
    )
