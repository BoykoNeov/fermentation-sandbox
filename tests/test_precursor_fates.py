"""Tests for the precursors' non-Ehrlich fates — the sink D-100 left out (decision D-104).

Before D-104 the Ehrlich re-route was each precursor's ONLY consumer, so the model attributed
**100% of consumed leucine to isoamyl alcohol**; Crépin *et al.* 2017 measures 77-86% of it going
to protein. :class:`PrecursorNonEhrlichFates` draws ``f/(1-f)`` times the re-route's own
per-species draw, so consumed precursor splits exactly ``f : (1-f)`` between every non-Ehrlich
fate and the alcohol.

This suite pins: the closure algebra (carbon + nitrogen), that the imposed split is *exactly*
``f`` at the ProcessSet level, that production is untouched, the joint-nitrogen-budget guard the
D-32 swap's ``psi`` no longer covers alone, and undosed isolability.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.chemistry import (
    carbon_mass_fraction,
    sugar_species,
)
from fermentation.core.kinetics import (
    AminoAcidAssimilation,
    FuselAlcoholsEhrlich,
    FuselAminoAcidReroute,
    PrecursorNonEhrlichFates,
    non_ehrlich_fraction_param,
)
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon, total_nitrogen

SINK = PrecursorNonEhrlichFates.name
REROUTE = FuselAminoAcidReroute.name
PRODUCER = FuselAlcoholsEhrlich.name

#: Every precursor the re-route draws — methionine is NOT among them (no Ehrlich alcohol), which
#: is why it has no ``f_non_ehrlich_*`` entry.
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
        base / "aging.yaml",
        base / "thermal.yaml",
    ).resolve()


#: Minebois 2025's measured PROTEIN share of consumed phenylalanine (Sc, Fig. 6A). A hard LOWER
#: BOUND on the non-Ehrlich lump (protein is one of the lump's fates, so the lump cannot be
#: smaller). D-117 shipped it as the *value*; since D-118 it is the **floor of the band**.
_PHE_PROTEIN_SHARE_BOUND = 0.531

#: The measured non-Ehrlich LUMP, 1 − 0.025 (Minebois 2025: 2.5% of consumed phenylalanine reaches
#: 2-phenylethanol). **[D-118: THE MODEL CARRIES THIS NOW — it is the shipped value.]** Until the
#: de-novo phenylpyruvate route landed it could not: at 0.975 the D-104 sink's joint carbon refund
#: reached 1.125× growth's own draw, past the sparing credit's ceiling, inventing extracellular
#: sugar in ``S``. It lived here as an unsampleable constant for exactly that reason. With
#: 2-phenylethanol's phenylalanine branch capped by ``f_de_novo_2_phenylethanol`` the same value
#: measures 0.584×, and the YAML band reaches it safely;
#: see :func:`test_the_de_novo_route_is_what_makes_the_sourced_lump_shippable`.
_PHE_MEASURED_LUMP = 0.975


def _run(*, amino_acids_gpl: float | None, days: float = 14.0, yan_mgl: float = 250.0):
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if amino_acids_gpl is not None:
        initial["amino_acids_gpl"] = amino_acids_gpl
    scenario = Scenario(
        name="wine-precursor-fates",
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


# -- the parameter contract ---------------------------------------------------


def test_every_reroute_precursor_has_a_sourced_non_ehrlich_fraction(full_params):
    # A sixth Ehrlich alcohol added to FUSEL_SPECS must fail loudly here rather than silently
    # acquiring no sink (which would restore the D-100 defect for that one molecule).
    for species in _PRECURSORS:
        key = non_ehrlich_fraction_param(species)
        assert key in full_params, f"{species} has no {key}"
        assert 0.0 <= full_params[key] < 1.0


def test_methionine_has_no_non_ehrlich_fraction(full_params):
    # Deliberate absence (D-104): methionine has no Ehrlich alcohol, so the sink has no draw to
    # scale for it and the parameter would never be read. Pins the reasoning against a
    # well-meaning future edit that "completes the set".
    assert "f_non_ehrlich_methionine" not in full_params
    assert "methionine" not in _PRECURSORS


# -- the split algebra --------------------------------------------------------


def _state(schema: StateSchema, params: Mapping[str, float], aa: float = 0.05) -> FloatArray:
    y = schema.zeros()
    y[schema.slice("X")] = 2.0
    y[schema.slice("S")] = 100.0
    y[schema.slice("N")] = 0.1
    y[schema.slice("T")] = params["T_ref"]
    for species in _PRECURSORS:
        y[schema.slice(species)] = aa
    return y


def test_the_realised_split_is_exactly_the_sourced_fraction(full_params):
    # THE CRUX (D-104). The sink must impose f : (1-f) on the CONSUMED precursor exactly — not
    # approximately, and at every instant, so the split holds on any trajectory. Compared at the
    # derivative level against the re-route's own draw, per species.
    schema = wine_schema()
    y = _state(schema, full_params)
    reroute = FuselAminoAcidReroute().derivatives(0.0, y, schema, full_params)
    sink = PrecursorNonEhrlichFates().derivatives(0.0, y, schema, full_params)
    for species in _PRECURSORS:
        ehrlich = -float(reroute[schema.slice(species)][0])
        lump = -float(sink[schema.slice(species)][0])
        assert ehrlich > 0.0 and lump > 0.0
        f = full_params[non_ehrlich_fraction_param(species)]
        # consumed = ehrlich + lump; the lump's share of it must be exactly f
        assert lump / (ehrlich + lump) == pytest.approx(f, rel=1e-12), species


def test_the_sink_does_not_touch_production_or_the_identity_agnostic_pools(full_params):
    # The D-100 decoupling, preserved: arginine does not make higher alcohols, so this sink must
    # never touch the pools the yeast swap / MLF / Brett / Maillard consumers live on. If it did,
    # fusel production could starve bacterial growth again — the exact pathology D-100 fixed.
    schema = wine_schema()
    y = _state(schema, full_params)
    d = PrecursorNonEhrlichFates().derivatives(0.0, y, schema, full_params)
    for pool in ("amino_acids", "amino_acids_generic", *(s.pool for s in FUSEL_SPECS)):
        assert float(d[schema.slice(pool)][0]) == 0.0, pool


def test_a_fraction_at_or_above_one_raises_rather_than_returning_inf(full_params):
    # f → 1 demands an infinite draw against a finite alcohol. An ensemble sampling the
    # uncertainty band is the realistic way this is ever reached; a silent inf would poison the
    # solver instead of failing.
    schema = wine_schema()
    params = dict(full_params)
    params[non_ehrlich_fraction_param("leucine")] = 1.0
    with pytest.raises(ValueError, match="outside"):
        PrecursorNonEhrlichFates().derivatives(0.0, _state(schema, params), schema, params)


# -- conservation -------------------------------------------------------------


def test_carbon_and_nitrogen_close_over_a_dosed_run(full_params):
    # Crown jewel: the sink moves precursor carbon to sugar and precursor nitrogen to ammonium.
    # Atoms only move between weighted pools, so both ledgers close to solver tolerance. (What is
    # a stand-in is the DESTINATION — Crépin's 20% unrecovered is booked as biomass — not the
    # balance; see the module docstring.)
    traj, compiled = _run(amino_acids_gpl=1.0)
    schema = compiled.process_set.schema
    pv = compiled.param_values
    assert_conserved(
        traj, total_carbon(schema, biomass_carbon_fraction=pv["biomass_C_fraction"]), rtol=1e-6
    )
    assert_conserved(
        traj, total_nitrogen(schema, biomass_nitrogen_fraction=pv["biomass_N_fraction"]), rtol=1e-6
    )


def test_the_sink_never_drives_a_precursor_negative(full_params):
    # Both consumers ride the SAME gate, which → 0 as the pool empties, so the combined draw is
    # self-limiting however large f/(1-f) gets. Near an empty pool the gate is LINEAR in it
    # (aa/(K·f_i + aa) → aa/(K·f_i)), so the draw decays exponentially and cannot reach zero in
    # exact arithmetic, let alone cross it. Any negative here is therefore BDF undershoot, and the
    # claim is structural — not an epsilon.
    #
    # **Stated scale-relatively since D-106**, which is what it always meant. The old absolute
    # -1e-9 was calibrated to a pre-D-106 undershoot (phenylalanine -4.8e-10); charging the Ehrlich
    # CO2 grew the draw and phenylalanine's undershoot went to -2.1e-9 — through a bound that was
    # never a physical statement, only that run's noise floor.
    #
    # That undershoot is solver noise, not the bigger draw, and the ORDERING proves it:
    # phenylalanine has the SMALLEST draw increase of the five (+12.5%, C9→C8) yet grew its
    # undershoot the most (4.3x), while valine and isoleucine — drawing 25%/20% more — improved to
    # ~1e-15. If the draw drove the undershoot that ordering would have to reverse.
    traj, compiled = _run(amino_acids_gpl=1.0)
    schema = compiled.process_set.schema
    for species in _PRECURSORS:
        pool = traj.y[schema.slice(species)][0]
        # Relative to the pool the run actually started with: a structural failure (a draw that
        # outruns its own gate) lands orders of magnitude above this, while numerical undershoot
        # sits far below it. The absolute ceiling keeps a tiny seeded pool from making it vacuous.
        assert pool.min() > -1e-6 * float(pool[0]), species
        assert pool.min() > -1e-7, species


def _worst_joint_refund(traj, compiled) -> tuple[float, float]:
    """Peak joint (swap + sink) N and C refund, as a multiple of growth's own draw.

    Only states with real growth are considered — at ``base_dx`` → 0 both quantities vanish and
    the ratio is a meaningless 0/0.
    """
    from fermentation.core.kinetics.growth import biomass_growth_rate

    schema = compiled.process_set.schema
    ps, pv = compiled.process_set, compiled.param_values
    worst_n = worst_c = 0.0
    for i in range(traj.y.shape[1]):
        y = traj.y[:, i]
        base_dx = biomass_growth_rate(y, schema, pv)
        if base_dx <= 1e-6:
            continue
        n = c = 0.0
        for proc in (AminoAcidAssimilation(), PrecursorNonEhrlichFates()):
            if not ps.is_enabled(proc.name):
                continue
            d = proc.derivatives(float(traj.t[i]), y, schema, pv)
            n += float(d[schema.slice("N")][0])
            # Weight each sugar slot by its own species' carbon fraction — the ledger's own rule,
            # so this serves wine's single slot and beer's three identically.
            s_slice = schema.slice("S")
            for offset, species in enumerate(sugar_species(schema)):
                c += float(d[s_slice.start + offset]) * carbon_mass_fraction(species)
        worst_n = max(worst_n, n / (pv["biomass_N_fraction"] * base_dx))
        worst_c = max(worst_c, c / (pv["biomass_C_fraction"] * base_dx))
    return worst_n, worst_c


def test_the_joint_carbon_refund_never_creates_sugar(full_params):
    # THE GUARD that matters (D-104, the owner's separate-Process-plus-guard call). Two Processes
    # now refund biomass carbon — the D-32 swap (on {arginine, generic}) and this sink (on the
    # C-RICH precursors: leucine C:N 5.5, phenylalanine 7.7, both ABOVE biomass's 4.3). Their
    # joint guarantee is NOT structural the way D-32's was for the swap alone: nothing bounds
    # f/(1-f) against growth's draw. The refund is a SPARING CREDIT and you cannot spare more sugar
    # than growth was charged: past 1.0 it CREATES extracellular sugar in `S`. That is a
    # mass-balance violation, not a metabolic pathway (D-117) — nothing puts glucose back into the
    # must, and gluconeogenesis would not either (intracellular G6P, never secreted). Pin it.
    traj, compiled = _run(amino_acids_gpl=1.0)
    _, worst_c = _worst_joint_refund(traj, compiled)
    assert worst_c < 1.0, f"joint C refund reached {worst_c:.2f}x growth's draw — creates sugar"


def test_the_net_sugar_derivative_is_never_positive(full_params):
    # The same guarantee at the ProcessSet level, where it actually has to hold: whatever the
    # individual Processes do, the SUMMED right-hand side must never make sugar appear.
    traj, compiled = _run(amino_acids_gpl=1.0)
    schema, ps, pv = compiled.process_set.schema, compiled.process_set, compiled.param_values
    worst = max(
        float(ps.total_derivatives(float(traj.t[i]), traj.y[:, i], pv)[schema.slice("S")].sum())
        for i in range(traj.y.shape[1])
    )
    assert worst <= 0.0, f"net dS/dt reached {worst:+.3e} g/L/h — sugar created"


def test_the_joint_nitrogen_refund_exceeds_growths_draw_at_pitch_and_that_is_deamination():
    # DOCUMENTED, NOT A BUG — and it falsifies a claim D-32 makes about itself.
    #
    # D-32's docstring argues its N refund is "≤ f_N·base_dx (growth's nitrogen draw) for all
    # ψ·gate ≤ 1 — never over-refunds, so no deamination branch is needed in v1". That holds for
    # the swap ALONE. With this sink on it is FALSE: at a vigorously growing state (NOT a
    # degenerate tail) the joint refund reaches **1.090x** (measured at D-248; it was 1.040x from
    # D-104 to D-106 and **1.171x** from D-106 to D-248).
    #
    # **D-248 moved it 1.171 -> 1.090, and the direction matters more than the size.** Neither the
    # swap nor the sink changed; what changed is the trajectory this peak is scanned over. Uptake
    # keeps ammonium non-zero far longer, so growth runs at states it never used to reach and the
    # worst-case ratio lands elsewhere. Note WHICH edge that consumes: the load-bearing bound here
    # is the LOWER one (`> 1.0` is the claim — the joint refund really does exceed growth's draw,
    # which IS the deamination), and its margin has HALVED, 0.171 -> 0.090. The ceiling is authored
    # tripwire margin; the floor is the finding, and it is now the tighter of the two. The joint
    # CARBON refund is unmoved at 0.5838x (measured at D-248 against the recorded 0.584).
    #
    # It is physical. The refund is always the drawn amino acid's own nitrogen; whether the NET
    # is negative (aa nitrogen spares ammonium growth would have drawn) or positive (the excess is
    # deaminated and released) falls out of the arithmetic rather than needing its own branch. So
    # the over-refund IS the deamination — no branch required, just a claim to correct.
    #
    # **D-106 moved this 1.040 → 1.171 and the ceiling was re-authored, which needs its reason
    # stated.** Charging the Ehrlich decarboxylation CO2 made the re-route draw a FULL mole of
    # precursor per alcohol instead of (n-1)/n, and a full mole carries a full mole of nitrogen —
    # so the deamination rose by the same ~12.6% the consumption did. The band moved because the
    # model got MORE right, and the CO2 charge does not rest on this band: it is fixed by atom
    # counts and pinned by a mutation-tested driven test. That is the difference between this and
    # the D-103 trap, where a band was nearly used to ACQUIT a model whose correctness was not
    # independently established.
    #
    # **The qualitative call the tripwire exists to force**: is "slight deamination at pitch" still
    # fair at 1.090x? Yes — more comfortably than at the 1.171x it was asked of at D-106. The net
    # ammonium release is 9% of growth's draw (17% at D-106, 4% at D-104); the direction, the
    # mechanism, and the conservation are unchanged, and 9% is nowhere near inverting the nitrogen
    # story (which would need the
    # precursors to become a dominant N SOURCE, i.e. multiples of growth's draw). The ceiling below
    # is AUTHORED, not sourced — its width is tripwire margin over the measured value, not physics.
    #
    # It is bounded and it conserves: carbon stays at ~0.55x (the test above), and total nitrogen
    # closes to 1e-14 (the conservation test) because the nitrogen is TRANSFERRED from the
    # precursor pools, never created. Pinned so a ψ/dose/fraction change that pushes this far
    # higher — where "slight deamination at pitch" would stop being a fair description — fails
    # here instead of quietly inverting the nitrogen story.
    traj, compiled = _run(amino_acids_gpl=1.0)
    worst_n, _ = _worst_joint_refund(traj, compiled)
    assert 1.0 < worst_n < 1.20, (
        f"joint N refund {worst_n:.3f}x — outside the documented band. D-248 measured 1.090; the "
        "LOWER edge is the one carrying the claim, so a value at or under 1.0 means the joint "
        "refund no longer exceeds growth's draw and D-32's corrected docstring goes back"
    )


def test_the_de_novo_route_is_what_makes_the_sourced_lump_shippable(full_params):
    """THE D-117 BLOCKER, RESOLVED — and this test is the inverse of the one it replaces (D-118).

    Its predecessor, ``test_the_sourced_lump_breaks_the_carbon_refund_guard``, was **designed to
    fail** when the de-novo phenylpyruvate route landed, and instructed its own deletion. Deleting
    it would have thrown away the executable knowledge, so it is inverted instead: the same three
    quantities are measured, and each now asserts the *opposite* outcome **plus the counterfactual
    that proves the route is the cause**.

    **What changed.** D-117 shipped ``f_non_ehrlich_phenylalanine`` at Minebois's *protein* share
    (0.531, a lower bound) because the measured **lump** (0.975) drove the D-104 sink's joint
    carbon refund to **1.125x growth's own draw** — past the sparing credit's ceiling, inventing
    extracellular sugar in ``S``. The cause was never the parameter: the model charged *all* of its
    ``k``-calibrated 2-phenylethanol to consumed phenylalanine, while the molecule is
    overwhelmingly built de novo. D-118's
    :data:`~fermentation.core.kinetics.carbon_routing.DE_NOVO_FUSEL_ROUTES` caps that branch, and
    the measured lump now ships.

    **The counterfactual is the load-bearing half.** Asserting only that 0.975 is now safe would
    not distinguish "the route fixed it" from "something else drifted and happens to mask it" —
    the D-108/D-109 vacuity trap, where a test measures an outcome it also assumes. So this
    switches the route **off** at the same parameter values and pins that the identical 0.975
    breaches again.
    """
    traj, compiled = _run(amino_acids_gpl=1.0)
    entry = compiled.parameters["f_non_ehrlich_phenylalanine"]

    # The MEASURED lump ships now — not the bound.
    assert entry.value == pytest.approx(_PHE_MEASURED_LUMP)

    # ...and its band is real rather than pinned shut, which was D-117's whole compromise.
    assert entry.uncertainty is not None
    assert entry.uncertainty.high > entry.uncertainty.low, (
        "f_non_ehrlich_phenylalanine's band is zero-width again — the D-117 pin was a workaround "
        "for a breach the de-novo route removed, not a permanent shape"
    )
    assert entry.uncertainty.high == pytest.approx(_PHE_MEASURED_LUMP)
    assert entry.uncertainty.low == pytest.approx(_PHE_PROTEIN_SHARE_BOUND)

    # THE BAND'S WHOLE SPAN MUST BE SAFE, not merely its mode — the D-117 follow-up lesson
    # ("a rejected value is not rejected until it is unreachable"). The refund scales
    # monotonically in f, so the top of the band is the worst case and pinning it suffices.
    worst_n, worst_c = _worst_joint_refund(traj, compiled)
    assert worst_c < 1.0, f"joint C refund {worst_c:.3f}x at the shipped lump — creates sugar"
    assert 1.0 < worst_n < 1.20, (
        f"joint N refund {worst_n:.3f}x — outside the documented band. D-248 measured 1.090; the "
        "LOWER edge is the one carrying the claim, so a value at or under 1.0 means the joint "
        "refund no longer exceeds growth's draw and D-32's corrected docstring goes back"
    )

    # THE COUNTERFACTUAL: switch the de-novo route off and the SAME lump breaches again. This is
    # what makes the assertions above attributable to the route rather than to drift.
    pv = dict(compiled.param_values)
    pv["f_de_novo_2_phenylethanol"] = 0.0
    dur = compiled.t_span_h[1]
    without = simulate(
        compiled.process_set,
        pv,
        compiled.y0,
        compiled.t_span_h,
        t_eval=np.linspace(0.0, dur, int(dur) + 1),
    )
    assert without.success, without.message

    class _Shim:
        process_set = compiled.process_set
        param_values = pv

    _, worst_c_off = _worst_joint_refund(without, _Shim())
    assert worst_c_off > 1.0, (
        f"with the de-novo route OFF the sourced lump no longer breaches (joint C refund "
        f"{worst_c_off:.3f}x) — either the sink's refund destination changed or the guard moved. "
        "The route is supposed to be the ONLY thing standing between 0.975 and a breach; if it is "
        "not, this beat's causal claim is wrong and D-118 needs re-deriving"
    )
    # ...and the route must be doing real work, not rounding: an order of magnitude between them.
    assert worst_c_off > 1.5 * worst_c


def test_the_de_novo_share_stays_above_its_analytic_breach_point():
    """The band's floor is a MODEL limit, and it is derived here rather than trusted (D-118).

    ``f_de_novo_2_phenylethanol``'s lower bound is not a measurement spread — it is the point
    where the phenylalanine branch grows enough to break the carbon guard at the shipped lump.
    The breach point is analytic: the sink's draw scales ``f/(1-f)``, so the de-novo-capped phe
    refund is ``(1 - f_de_novo) * f/(1-f)``, and it reaches the pre-route (and safely-shipping)
    0.531 configuration's ``0.531/0.469`` when ``(1 - f_de_novo) * 39 == 1.132``.

    Pinned because the YAML note *claims* this number, and a claimed number nobody recomputes is
    the D-96/D-109 class of defect this project keeps re-learning.

    **Evaluated at phe's band MAXIMUM, not its nominal (D-155).** Both parameters are sampled and
    drawn INDEPENDENTLY, and ``breach()`` is strictly increasing in ``f`` — so the constraint that
    must hold for every joint draw is the one at ``f``'s top, not at its point value. Today those
    coincide (D-118 put the mode on the high edge because 0.531 is a hard measured protein floor),
    which is why the nominal form was safe; but it was safe by a coincidence enforced in a
    *different* test, ``test_the_de_novo_route_is_what_makes_the_sourced_lump_shippable``, whose
    subject is the D-117 band shape and which has no idea it is holding a conservation property
    (it pins ``uncertainty.high == _PHE_MEASURED_LUMP``). Reading the edge
    here makes this test self-sufficient. Measured, not assumed: widening phe's high edge alone
    left this test GREEN before the change and turns it red after (D-155's mutation A).
    """
    _, compiled = _run(amino_acids_gpl=1.0, days=1.0)
    entry = compiled.parameters["f_de_novo_2_phenylethanol"]
    phe = compiled.parameters["f_non_ehrlich_phenylalanine"]
    assert phe.uncertainty is not None
    f = phe.uncertainty.high  # the joint worst case; see the docstring

    breach_point = 1.0 - (_PHE_PROTEIN_SHARE_BOUND / (1.0 - _PHE_PROTEIN_SHARE_BOUND)) / (
        f / (1.0 - f)
    )
    # SAFETY FIRST, sanity second — deliberately this order. The pin below has an abs=5e-4
    # tolerance against a margin of 3.07e-5, i.e. ~16x looser than the headroom it sits beside,
    # so a drift that consumes the entire margin passes it. That is not a hole (the assert here
    # is exact and binding) but it does mean the pin must not be the first thing to fire, or an
    # unsafe band reports itself as "the breach point isn't 0.971".
    assert entry.uncertainty is not None
    assert entry.uncertainty.low >= breach_point, (
        f"f_de_novo_2_phenylethanol's band floor {entry.uncertainty.low} is below the analytic "
        f"breach point {breach_point:.6f} — an ensemble draw can now break the carbon guard. "
        f"Computed at f_non_ehrlich_phenylalanine's band HIGH edge ({f}), which is the joint "
        "worst case because breach() increases in f and the two are drawn independently. So this "
        "is red if EITHER band moved: phe's top edge up, or the de-novo floor down."
    )
    # ...and the monotonicity that makes the band's top the worst case, checked not assumed.
    lower = 1.0 - (_PHE_PROTEIN_SHARE_BOUND / (1.0 - _PHE_PROTEIN_SHARE_BOUND)) / (
        phe.value / (1.0 - phe.value)
    )
    assert lower <= breach_point, (
        "breach() is no longer increasing in f_non_ehrlich_phenylalanine, so its band's TOP is "
        "not the joint worst case and this test is checking the wrong end of the band"
    )
    assert breach_point == pytest.approx(0.971, abs=5e-4), breach_point


# -- isolability --------------------------------------------------------------


def test_an_undosed_run_is_byte_for_byte_the_validated_core():
    # Prime directive #3. Undosed, the compile seam disables the sink outright; even enabled every
    # gate is exactly 0 at aa=0 and the re-route's draw it scales is 0 too.
    dosed_off, _ = _run(amino_acids_gpl=None)
    schema = wine_schema()
    for spec in FUSEL_SPECS:
        assert float(dosed_off.y[schema.slice(spec.pool)][0][-1]) > 0.0  # alcohols still made
    for species in _PRECURSORS:
        assert float(dosed_off.y[schema.slice(species)][0][-1]) == 0.0  # pools stay empty


def test_the_compile_seam_disables_the_sink_undosed():
    for dose, expected in ((None, False), (1.0, True)):
        _, compiled = _run(amino_acids_gpl=dose, days=1.0)
        assert compiled.process_set.is_enabled(SINK) is expected


def test_the_sink_is_speculative_and_only_taints_tiers_when_enabled():
    schema = wine_schema()
    procs = [FuselAlcoholsEhrlich(), FuselAminoAcidReroute(), PrecursorNonEhrlichFates()]
    off = ProcessSet(schema, procs)
    off.disable(SINK)
    off.disable(REROUTE)
    on = ProcessSet(schema, procs)
    assert off.tier_of("N") is Tier.VALIDATED
    assert on.tier_of("N") is Tier.SPECULATIVE
    assert PrecursorNonEhrlichFates.tier is Tier.SPECULATIVE


# -- D-257: why the temporal fusel repair is BLOCKED HERE, measured rather than asserted -------

#: The uniform factor a level-preserving timing repair must apply to all five higher-alcohol
#: rate constants. Derived at D-257, not chosen: adding a nitrogen-independent de-novo term to
#: ``fusel_rate_shape`` and fitting its weight to Rollero's NT/EF fraction (48.2 % on his SM250)
#: lands finished isoamyl alcohol at 426 mg/L on the D-112 anchor must, so holding the Wang 2024
#: 172 mg/L anchor costs exactly this factor. The five share one rate shape, so one factor holds
#: every species' mean at once and no D-99 ratio moves.
_D257_LEVEL_PRESERVING_K_FACTOR = 0.4033


def _run_with_scaled_fusel_k(factor: float):
    """The D-104 fixture again, with every higher-alcohol rate constant scaled by ``factor``."""
    scenario = Scenario(
        name="d257-fusel-k-sweep",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25, "amino_acids_gpl": 1.0},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=14.0,
    )
    compiled = compile_scenario(scenario, strict=True)
    params = dict(compiled.param_values)
    for spec in FUSEL_SPECS:
        params[spec.k_param] = params[spec.k_param] * factor
    dur = compiled.t_span_h[1]
    traj = simulate(
        compiled.process_set,
        params,
        compiled.y0,
        compiled.t_span_h,
        t_eval=np.linspace(0.0, dur, int(dur) + 1),
    )
    assert traj.success, traj.message
    return traj, compiled


def _surviving_fraction(traj, compiled, pool: str) -> float:
    series = traj.y[compiled.process_set.schema.slice(pool)][0]
    assert float(series[0]) > 0.0, f"vacuous: the must carries no {pool}"
    return float(series[-1]) / float(series[0])


def test_precursor_consumption_rides_the_ALCOHOL_rate_and_that_is_what_blocks_the_d257_repair():
    r"""THE BLOCKER, pinned inside one model and one pair of runs (decision D-257).

    :class:`PrecursorNonEhrlichFates` consumes precursor at ``f/(1-f)`` times the Ehrlich
    re-route's own draw, and that draw is proportional to the higher-alcohol PRODUCTION RATE.
    So the model has an amino acid's **dominant** fate - protein, 77-86 % of consumed leucine by
    Crepin's measurement, and the whole reason this sink exists - riding how much fusel alcohol
    is being made, rather than riding growth. Nothing in the must's biology works that way.

    **Why this is the finding of D-257 rather than a note.** D-256 measured a large temporal
    defect: the model makes 100.6 % of its isoamyl alcohol before its own nitrogen gate shuts,
    at peak fermentative flux, where Rollero measures 42-54 % in by then across six
    fermentations. Every repair for that defect must slow production in the growth window and
    make up the rest afterwards, and holding the sourced level anchor then forces the five ``k``
    down by :data:`_D257_LEVEL_PRESERVING_K_FACTOR`. This test measures what that costs: the
    must's phenylalanine stops being consumed. **The repair was built, measured and reverted on
    exactly this** - the patch is kept at
    ``M:\claud_projects\temp\ferment\d257-fusel-de-novo\attempted-repair.patch``.

    **Verified by mutation, and the mutation is sharper than the test.** Toggling the sink OFF
    and repeating the same sweep, the must's phenylalanine stops caring about the higher-alcohol
    rate entirely: 97.8 % -> 99.1 % survives (a 1.3-point swing) against 20.3 % -> 65.8 % with
    the sink on (45.5 points). So this sink is not merely *a* consumer of phenylalanine, it is
    effectively the ONLY one - :class:`AminoAcidAssimilation` leaves 97.8 % of it in the must -
    and what it consumes is set by how much 2-phenylethanol is being made. Probe kept at
    ``M:\claud_projects\temp\ferment\d257-fusel-de-novo\mutation.py``.

    **The anti-vacuity arm is the load-bearing half**, and it is what makes this an attribution
    rather than a restatement: the same scaling leaves biomass and the run's TOTAL nitrogen
    consumption essentially where they were. If those moved too, "consumption rides the alcohol
    rate" would just be "the whole ferment moved". They do not - so the coupling is specifically
    between higher-alcohol production and the amino-acid pools, which is the inversion.
    """
    base_traj, base_compiled = _run_with_scaled_fusel_k(1.0)
    slow_traj, slow_compiled = _run_with_scaled_fusel_k(_D257_LEVEL_PRESERVING_K_FACTOR)

    base_phe = _surviving_fraction(base_traj, base_compiled, "phenylalanine")
    slow_phe = _surviving_fraction(slow_traj, slow_compiled, "phenylalanine")

    assert base_phe < 0.30, (
        f"the shipped model leaves {base_phe:.1%} of the must's phenylalanine unconsumed; D-257 "
        "measured 20.3 %. If this has risen the baseline for the comparison below has moved and "
        "the coupling must be re-measured rather than re-pinned"
    )
    assert slow_phe > 0.55, (
        f"scaling the five higher-alcohol rate constants by "
        f"{_D257_LEVEL_PRESERVING_K_FACTOR} leaves {slow_phe:.1%} of the must's phenylalanine "
        f"unconsumed against {base_phe:.1%} unscaled. D-257 measured 65.8 %. If this has fallen, "
        "precursor consumption no longer rides the alcohol rate - which would UNBLOCK the "
        "temporal fusel repair D-257 refused, so re-open that record rather than re-pinning here"
    )

    # -- the anti-vacuity arm: the rest of the ferment did NOT move with it ---------------------
    schema = base_compiled.process_set.schema
    base_x = float(base_traj.y[schema.slice("X"), -1][0])
    slow_x = float(slow_traj.y[slow_compiled.process_set.schema.slice("X"), -1][0])
    assert base_x > 0.0, "vacuous: no biomass was built"
    assert abs(slow_x - base_x) < 0.05 * base_x, (
        f"final biomass moved {base_x:.4f} -> {slow_x:.4f} g/L under the k scaling. This test "
        "attributes the phenylalanine change to the ALCOHOL rate specifically; if the whole "
        "ferment moved, that attribution is not established and the finding must be re-measured"
    )

    base_n = _surviving_fraction(base_traj, base_compiled, "N")
    slow_n = _surviving_fraction(slow_traj, slow_compiled, "N")
    assert abs(slow_n - base_n) < 0.02, (
        f"the extracellular nitrogen slot moved {base_n:.4f} -> {slow_n:.4f} of its initial "
        "value under the k scaling - the same objection as the biomass arm above"
    )
