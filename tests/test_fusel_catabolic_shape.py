"""The fusel catabolic SHAPE, measured — the receipts for D-112.

D-111 left "the leucine-derived isoamyl shortfall (1.12% vs Rollero 3.4-17.3%)" as the sharpest
open item, framed as D-103's gate *shape* that the keto-acid node would fix. **D-112 measured it
and retired that framing.** This suite pins the three measurements so a future beat cannot quietly
re-inherit the stale story:

1. **D-103's gate-shape SPREAD is absorbed by the D-104 non-Ehrlich sink** — it is large with the
   sink off (isoamyl ~6% vs propanol ~67%) and compresses to a uniform low band with it on. So the
   "minor alcohols are wildly over-attributed" defect is gone, and only isoamyl (UNDER) survives.
   **The "uniform low band" half was an xfail from D-244 to D-247** — propanol read 0.2256
   sink-on against a 0.20 band that IS the sourced 80% floor written upside down; see
   ``test_every_sink_on_share_sits_inside_rollers_low_band``, split out at D-245 so that failure
   stopped burying the compression assertion beside it. **D-248 closed it at 0.1938 without
   touching the 0.20**: un-coupling assimilable-nitrogen uptake from growth demand consumes the
   whole must, which restores the biomass denominator every de-novo share is measured against.
   The compression itself held throughout.
2. **Isoamyl sits on its ``(1-f)`` mass-conservation ceiling**, which no sourcing-layer change —
   the keto-acid node included — can lift: a gate cap (the "obvious" fix) does not move it, because
   leucine is too scarce to persist under any draw rate. **The ``(1-f)`` here is leucine's and is
   still literally right** — D-245 corrected that notation for VALINE alone, whose Ehrlich carbon
   splits again between isobutanol and isoamyl (see
   :func:`~tests.test_fusel_keto_acid_node.ehrlich_primary_share`). Finding 2 was audited against
   that repair and does not move; the note is here because this is the line a reader will grep.
3. **``k_isoamyl_alcohol`` is correctly calibrated** to the Wang 2024 172 mg/L anchor at typical
   must nitrogen with no amino-acid dose — so the ~2x isoamyl over-production in the D-109
   characterization must is the ``amino_acids_gpl=1.0`` dose (Finding 4), not a mis-set ``k``.
   **D-245 corrects the MECHANISM inside that finding, and only the mechanism.** The dose is still
   the cause and ``k`` is still right (the undosed anchor test below is untouched by any of this),
   but the route was never "the dose's deamination-N sustaining the fusel gate": the dose is what
   made a scenario's declared YAN differ from the total it held, so Coleman's yield fit was read
   at a poorer must than the one simulated and the model built ~2x the biomass. Measured on ONE
   knob at D-245 — sweep ``biomass_N_fraction`` alone across this fixture and isoamyl runs 353.5
   → 191.5 mg/L, the whole of the 2x. The over-production is gone as of D-244.

All shares are measured from EXACT state differences (no quadrature — D-103's trapezoid error),
reusing the vetted D-109 harness where possible.

**D-113 adds the fourth receipt: the node vs the inversion.** D-112 named D-104's inverted split
(leucine 20.9% to protein vs Crépin's 77–86%) as the node's one live motivation and left "whether
D-111's valine-side route touches it" untested. ``test_the_valine_route_does_not_touch_leucines_
anabolic_split`` measures that it does not — leucine's Ehrlich draw is *bit-invariant* under the
route toggle (a headroom-fill never re-gates leucine) and total biomass with it, so leucine's
emergent protein share is route-invariant for any biomass composition. The route touches only
valine, the least-inverted species. Un-inverting leucine needs a de-novo-KIC route that RELIEVES
leucine's isoamyl demand — D-111 built a valine drain, not a leucine relief.
"""

from __future__ import annotations

import pytest

from fermentation.core.chemistry import CARBON_ATOMS, carbon_mass_fraction
from fermentation.core.kinetics.amino_acid_pools import depletion_gate as _orig_gate
from fermentation.core.kinetics.byproducts import ehrlich_draws
from fermentation.core.kinetics.carbon_routing import (
    FUSEL_SPECS,
    ISOAMYL_ALCOHOL,
    SECONDARY_FUSEL_ROUTES,
)
from fermentation.core.state import FloatArray, StateSchema
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from tests.test_fusel_keto_acid_node import (
    _D245_D120_LEGS_GONE_DIRECTION_BACK_AT_D248,
    _FERMENT_DAYS,
    _OTHER_PRECURSOR_CONSUMERS,
    _de_novo_share,
    _end,
    _run,
    _scenario,
    ehrlich_primary_share,
)

_SINK = "precursor_non_ehrlich_fates"
_PROPANOL = next(s for s in FUSEL_SPECS if s.pool == "propanol")


def _direct_catabolic_share_sink_off(spec) -> float:
    """Precursor-derived carbon fraction of ``spec``'s alcohol with the D-104 sink DISABLED.

    With the sink off the re-route is the precursor's ONLY consumer, so consumed precursor carbon
    times ``n/(n+1)`` (removing the D-106 decarboxylation CO2) is exactly the alcohol carbon sourced
    from it — an exact state difference, no quadrature. This is ~D-103's gate value; contrast with
    :func:`_de_novo_share` (sink on), which yields the share compressed by the precursor's
    *primary* Ehrlich residue (``(1-f)`` for every precursor but valine — see
    :func:`~tests.test_fusel_keto_acid_node.ehrlich_primary_share`).
    """
    prec = spec.precursor_amino_acid
    traj, schema = _run(aging=False, drop=(*_OTHER_PRECURSOR_CONSUMERS, _SINK))
    consumed = float(traj.y[schema.slice(prec), 0][0]) - _end(traj, schema, prec)
    made = _end(traj, schema, spec.pool)
    n_alc = CARBON_ATOMS[spec.species]
    alc_from_prec = consumed * carbon_mass_fraction(prec) * n_alc / (n_alc + 1.0)
    return alc_from_prec / (made * carbon_mass_fraction(spec.species))


def test_the_d104_sink_absorbs_the_d103_gate_shape_spread():
    """D-103's 11x catabolic spread is absorbed by the D-104 sink (decision D-112, finding 1).

    With the sink OFF the gate over-attributes minor alcohols exactly as D-103 diagnosed — propanol
    (small carbon draw against an abundant precursor) reads many times isoamyl (large draw, which
    throttles its own gate). With the sink ON the ``(1-f)`` multiplier compresses propanol out of
    the tens of percent. So "minor alcohols are wildly over-attributed" is no longer true.

    **SPLIT AT D-245, AND THE SPLIT IS A REPAIR RATHER THAN A CONVENIENCE.** D-244 marked this
    whole test ``xfail`` for its *middle* assertion (every alcohol inside Rollero's low band),
    which left the third one — propanol's compression, the finding's own payoff — **dead behind
    the failure**, never reached and never run. The band assertion now lives in
    :func:`test_every_sink_on_share_sits_inside_rollers_low_band` carrying the xfail alone, and
    the two assertions that still hold are guards again. Nothing is weakened: every assertion in
    the D-112 original still runs somewhere, on the same numbers.
    """
    iso_off = _direct_catabolic_share_sink_off(ISOAMYL_ALCOHOL)
    prop_off = _direct_catabolic_share_sink_off(_PROPANOL)
    # D-103's spread: sink OFF, propanol reads several times isoamyl (a minor alcohol over-read).
    assert prop_off > 4.0 * iso_off, (
        f"sink-OFF spread collapsed (propanol {prop_off:.1%} vs isoamyl {iso_off:.1%}) — D-103's "
        "gate-shape spread is the thing D-112 says the sink absorbs; if it is gone with the sink "
        "off too, finding 1's baseline has moved"
    )
    # And the compression is real: propanol specifically drops from tens-of-% to the low band.
    # This is the assertion D-244's xfail buried; it passes, and it did throughout (D-245).
    on_propanol = 1.0 - _de_novo_share(_PROPANOL)
    assert on_propanol < 0.5 * prop_off, (
        f"propanol's sink-ON share {on_propanol:.1%} is no longer half its sink-OFF "
        f"{prop_off:.1%} — the compression IS finding 1, so this is the finding moving"
    )


def test_every_sink_on_share_sits_inside_rollers_low_band():
    """Sink ON, every alcohol compresses into Rollero's uniform low band (D-112 finding 1).

    Split out of :func:`test_the_d104_sink_absorbs_the_d103_gate_shape_spread` at D-245 so its
    failure stops burying that test's third assertion.

    **The 0.20 here is NOT an independent threshold — it is the sourced 80 % de-novo floor written
    the other way up**, over all five alcohols instead of the four the floor's sources describe.
    That is why it failed for the same reason and on the same number as
    ``test_every_sourced_fusel_is_de_novo_dominated[propanol]``: propanol read 0.2256 sink-on, and
    raising 0.20 to admit it would have been lowering the sourced floor under a different name.

    **XFAIL REMOVED AT D-248, and the threshold was never touched.** Propanol reads **0.1938**
    sink-on, inside a 0.20 that is still the sourced floor inverted. What moved is the biomass
    denominator: un-coupling assimilable-nitrogen uptake from growth demand consumes the whole
    must, so more alcohol is made de novo off the sugar route against a precursor draw that was
    already supply-limited. Rollero's own uniform range tops out at 17.3 %, so the model is now
    inside the guard but still above the measurement it was drawn from — worth keeping in view.

    **Isobutanol is no longer the maximum and that is a D-245 repair, not a drift**: it read
    0.2400 under the harness defect ``ehrlich_primary_share`` fixes and reads 0.0814 measured on
    the residue the model actually routes.
    """
    on = {s.pool: 1.0 - _de_novo_share(s) for s in FUSEL_SPECS}
    assert max(on.values()) < 0.20, f"sink-on shares no longer compressed to the low band: {on}"


def test_isoamyl_sits_on_the_one_minus_f_mass_conservation_ceiling(monkeypatch):
    """Isoamyl's leucine share is a ``(1-f)`` ceiling a gate cap cannot move (D-112, finding 2).

    Leucine's only two AF fates are the non-Ehrlich lump ``f`` and isoamyl ``(1-f)``, so its share
    of isoamyl is ``(1-f) x leucine_C/isoamyl_C`` — a mass-conservation ceiling the model sits on.
    The keto-acid node reallocates HOW leucine reaches isoamyl (via KIC) but not HOW MUCH, and
    Crépin's ``f`` already prices in every non-isoamyl fate. The proof the ceiling binds:
    a gate cap — the intuitive "stop the gate over-claiming" fix — does NOT move the realised share,
    because leucine (~32 mg/L) is too scarce to persist under any draw rate (the advisor predicted
    the opposite; the probe refuted it). If a cap ever DOES move it, leucine has stopped exhausting
    and the ceiling argument must be re-derived.
    """
    base = 1.0 - _de_novo_share(ISOAMYL_ALCOHOL)
    assert 0.005 < base < 0.03, f"isoamyl leucine share {base:.4f} left its measured ~1.1% regime"

    # Cap the availability gate at 0.10 (the D-104 sink rides ehrlich_draws->depletion_gate, so both
    # the re-route and the sink shrink together — a faithful "catabolic cap with the sink intact").
    # String target so mypy does not need depletion_gate re-exported from byproducts.
    monkeypatch.setattr(
        "fermentation.core.kinetics.byproducts.depletion_gate",
        lambda *a, **k: min(_orig_gate(*a, **k), 0.10),
    )
    capped = 1.0 - _de_novo_share(ISOAMYL_ALCOHOL)
    assert abs(capped - base) < 0.005, (
        f"the gate cap moved isoamyl leucine share {base:.4f} -> {capped:.4f}: it should be inert "
        "(leucine exhausts under any draw rate), which is why no cap and no node lifts the ceiling"
    )


def _isoamyl_no_dose(yan_mgl: float, celsius: float = 20.0) -> float:
    """Finished isoamyl (mg/L) with NO amino-acid dose, at ``yan_mgl`` and ``celsius`` (D-112)."""
    scenario = Scenario(
        name="d112-anchor",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25, "amino_acids_gpl": 0.0},
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=celsius),
            TemperaturePoint(day=_FERMENT_DAYS, celsius=celsius),
        ],
        duration_days=_FERMENT_DAYS,
    )
    cs = compile_scenario(scenario)
    traj = simulate_scheduled(
        cs.process_set,
        cs.param_values,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
    )
    assert traj.success, traj.message
    return _end(traj, cs.process_set.schema, ISOAMYL_ALCOHOL.pool) * 1e3


@pytest.mark.parametrize("yan", [250.0, 300.0])
def test_k_isoamyl_alcohol_lands_the_wang_2024_anchor_at_typical_must_n(yan):
    """``k_isoamyl_alcohol`` lands finished isoamyl on its 172 mg/L anchor (D-112, finding 4).

    The provenance sets ``k`` to land 3-methylbutan-1-ol at the Wang 2024 mean **172 mg/L**. With
    NO amino-acid dose, at typical must nitrogen (250-300 mgN/L) and ``T_ref = 20 °C``, the model
    reproduces that mean — so the ~2x isoamyl over-production in the D-109 characterization must
    was the ``amino_acids_gpl=1.0`` dose, NOT a mis-set ``k``. This is what makes the isoamyl
    catabolic DENOMINATOR a probe artifact, not a calibration bug: the ceiling comparison in
    finding 3 rests on the ``k`` being right here.

    **THIS TEST IS THE CONTROL THAT SURVIVED D-244 UNMOVED, AND THAT IS WHY IT IS WORTH SAYING SO.**
    It declares no amino-acid dose, so its declared YAN already WAS the total the must held and
    D-244's re-reading of the yield fit changes nothing here — bit-for-bit. The dosed fixture next
    door moved 2x on the same commit. The pair is the cleanest statement of what D-244 did.
    D-245 corrects only the mechanism D-112 named for the dosed over-production ("the dose's
    deamination-N sustaining the fusel gate"); see the module docstring, finding 3.
    """
    isoamyl = _isoamyl_no_dose(yan)
    assert 140.0 < isoamyl < 205.0, (
        f"isoamyl {isoamyl:.1f} mg/L at YAN {yan:.0f} (no aa dose, 20 °C) left the 172 mg/L Wang "
        "anchor band — if k moved, D-112 finding 4 (the over-production is the aa dose, not the k) "
        "must be re-measured"
    )


# -- D-113: the node vs the inversion -- the receipt for "still untested against the inversion" ---


def _dosed_midferment_state(schema: StateSchema, params) -> FloatArray:
    """A mid-ferment dosed state with every precursor present and headroom on isoamyl.

    Built explicitly (not run) so the branch comparison below is an EXACT derivative-level read,
    free of solver noise — the shape of ``test_alpha_kb_production_is_exactly_threonine_indep``
    (D-109). ``N`` is left ample so leucine's isoamyl gate does not saturate to 1, leaving the
    headroom the valine branch fills (if it did saturate, the valine branch would clamp to 0 and the
    anti-vacuity half of the test below would be vacuous — which is the D-111 Finding 5 regime).
    """
    y = schema.zeros()
    y[schema.slice("X")] = 2.0
    y[schema.slice("S")] = 100.0
    y[schema.slice("N")] = 0.5
    y[schema.slice("T")] = params["T_ref"]
    for precursor in ("leucine", "isoleucine", "valine", "threonine", "phenylalanine"):
        y[schema.slice(precursor)] = 0.05
    return y


def _branch_carbon(schema: StateSchema, params, f_valine: float) -> dict[tuple[str, str], float]:
    """``{(precursor, alcohol): alcohol_carbon}`` from ``ehrlich_draws`` at this ``f_valine``."""
    p = dict(params)
    p["f_valine_to_isoamyl"] = f_valine
    y = _dosed_midferment_state(schema, p)
    return {
        (d.precursor.species, d.alcohol.pool): d.alcohol_carbon for d in ehrlich_draws(y, schema, p)
    }


def test_the_valine_route_does_not_touch_leucines_anabolic_split():
    """The keto-acid node's valine route leaves D-104's INVERTED species untouched (decision D-113).

    D-112 retired the leucine *shortfall* as the node's motivation and left one live reason: D-104's
    inverted anabolic split — model leucine **20.9%** to protein against Crépin's **77–86%**, order
    ``leu<ile<val<thr`` exactly reversed. D-111 built the valine → KIC → isoamyl route (the
    mechanism D-104 said the model lacked) on the valine side only, and D-112 left *"whether that
    touches the inversion"* explicitly untested. **It does not, and the reason is structural.**

    The inversion is a property of leucine's EMERGENT protein share, ``D_bio,leu / (D_bio,leu +
    D_ehrlich,leu)`` under a demand-anchored sink. That sink is the *rejected* D-100 prescription,
    not the shipped model (which IMPOSES the split via static ``f_non_ehrlich`` — measuring the
    shipped split would just read Crépin back out, the D-108 vacuity trap). So the node moves the
    inversion **iff** it moves one of those two inputs — and it moves NEITHER, for ANY biomass
    composition ``w_leu`` (no invented yeast-protein spectrum needed — the D-98 trap that
    reconstructing the sink would have sprung):

    * **``D_ehrlich,leu`` is bit-invariant (derivative level, EXACT).** ``ehrlich_draws`` gates
      leucine's isoamyl branch at ``gate_leu · fusel_carbon_isoamyl`` and clamps the valine branch
      to the headroom **above** it (D-111 Finding 5: the 122.7% over-claim clamp cut the *KIC*
      branch 31.8→9.1%, never leucine's 90.9%). So the route relieves leucine of **0%** of isoamyl
      demand — it adds a valine drain, it does not lift a leucine one.
    * **``D_bio,leu`` (∝ total biomass built) is invariant end-to-end**, and leucine consumed (∝
      ``D_ehrlich,leu`` via the sink's exact ``f:(1−f)`` split, ``f_non_ehrlich_leucine``
      route-invariant) with it.

    **Anti-vacuity — the route IS live and DOES move a species, just the wrong one.** The valine →
    isoamyl branch appears with the route and vanishes without it, so this is no disabled-Process
    no-op (D-106/D-108). The point is precisely that it touches **valine** — the *least*-inverted
    species (model 45.8 vs Crépin 41) — and never leucine, the most-inverted. Un-inverting leucine
    needs a route that RELIEVES leucine's isoamyl demand (isoamyl from de-novo KIC), which the model
    still lacks: D-111 built a valine drain, not a leucine relief. If a future beat claims the node
    fixed the inversion, this fails — the split it must move is leucine's, and leucine is invariant.
    """
    cs = compile_scenario(_scenario(aging=False))
    schema, pv = cs.process_set.schema, cs.param_values

    # (1) DERIVATIVE LEVEL, EXACT: the route leaves every branch bit-identical and only ADDS
    # valine → isoamyl. Leucine's isoamyl branch is the load-bearing one — it must be *identical*,
    # bit-for-bit, because the route is a headroom-fill that never re-gates leucine.
    on = _branch_carbon(schema, pv, 0.23)
    off = _branch_carbon(schema, pv, 0.0)
    leu, val = ("leucine", "isoamyl_alcohol"), ("valine", "isoamyl_alcohol")

    assert on[leu] == off[leu], (
        f"leucine's isoamyl branch moved with the valine route ({off[leu]!r} → {on[leu]!r}): the "
        "route is specified as a headroom-fill that never re-gates leucine — if it now reduces "
        "leucine's Ehrlich draw, it has begun relieving leucine's isoamyl demand and D-113's core "
        "claim (0% relief, inversion untouched) must be re-measured"
    )
    # Anti-vacuity: the route is live and it is valine it touches.
    assert on.get(val, 0.0) > 0.0 and val not in off, (
        f"the valine → isoamyl branch is not live (ON {on.get(val)!r}, OFF {off.get(val)!r}): the "
        "invariance above would then be a no-op on a disabled route, not a finding"
    )
    # Every OTHER branch is untouched too — the route adds exactly one branch, changes none.
    assert {k: v for k, v in on.items() if k != val} == off, (
        "the valine route changed a branch other than its own valine → isoamyl: it is specified to "
        "add exactly one branch and re-gate no precursor"
    )

    # (2) END TO END: leucine's two emergent-share inputs are invariant under the toggle, so the
    # share itself is — for any w_leu. Solver-noise tolerance (the route's real perturbation is via
    # valine's N/CO₂, which does not reach leucine's exhaustion; measured ~2e-6 relative).
    on_traj, sch = _run(aging=False, set_params={"f_valine_to_isoamyl": 0.23})
    off_traj, _ = _run(aging=False, set_params={"f_valine_to_isoamyl": 0.0})

    def leucine_consumed(traj) -> float:
        return float(traj.y[sch.slice("leucine"), 0][0]) - _end(traj, sch, "leucine")

    def biomass_built(traj) -> float:
        return _end(traj, sch, "X") + _end(traj, sch, "X_dead")

    leu_on, leu_off = leucine_consumed(on_traj), leucine_consumed(off_traj)
    bio_on, bio_off = biomass_built(on_traj), biomass_built(off_traj)
    d_ehrlich = abs(leu_on - leu_off) / leu_off
    d_bio = abs(bio_on - bio_off) / bio_off
    assert d_ehrlich < 1e-4, f"leucine consumed (∝ D_ehrlich,leu) moved {d_ehrlich:.1e} with route"
    assert d_bio < 1e-4, f"total biomass (∝ D_bio,leu) moved {d_bio:.1e} with the route"


# -- D-120: the de-novo CAP is the wrong instrument for isoamyl (and for the class) -------------


#: Minebois 2025 Fig. 6A, Sc (T73) at T4 — ``labelled µM / total µM`` per alcohol, i.e. the share
#: of that alcohol built FROM ITS AMINO ACID. Bar semantics are per-bar and MIXED (D-119): 2-PE's
#: large number is the total, isoamyl's and isobutanol's are the unlabelled segment.
_MINEBOIS_AMINO_ACID_SHARE = {
    "2_phenylethanol": 4.05 / 109.0,  # 3.72 %
    "isoamyl_alcohol": (77.0 + 104.0) / 3392.0,  # 5.34 % — leucine AND valine labels
    "isobutanol": 39.0 / 444.0,  # 8.78 %
}


def _amino_acid_share(traj, schema, params, pool: str) -> float:
    """This alcohol's carbon share sourced from amino acids — EVERY branch, not just the primary.

    :func:`_de_novo_share` counts an alcohol's ``precursor_amino_acid`` only. Isoamyl uniquely has
    a second (valine → KIC, :data:`SECONDARY_FUSEL_ROUTES`), and Minebois's de-novo share is
    ``1 − I.E`` over BOTH its labels — so comparing her number to a primary-only share would
    understate the model against its own target. Exact state differences throughout (D-103).
    """
    spec = next(s for s in FUSEL_SPECS if s.pool == pool)
    alcohol_carbon = _end(traj, schema, pool) * carbon_mass_fraction(spec.species)
    n_alc = CARBON_ATOMS[spec.species]
    route = next(r for r in SECONDARY_FUSEL_ROUTES if r.precursor == "valine")
    share_val_iso = params[route.share_param]

    def consumed(species: str) -> float:
        return float(traj.y[schema.slice(species), 0][0]) - _end(traj, schema, species)

    precursor = spec.precursor_amino_acid
    # valine's (1−f) splits again between isobutanol (primary) and isoamyl (secondary), so the
    # primary share is NOT (1−f) for it — using (1−f) would double-count the isoamyl branch.
    # THIS FILE HAD IT RIGHT AND ITS SIBLING DID NOT, from D-111 to D-245: the same rule is now
    # ONE helper both call, because two callers of one quantity disagreeing in prose while only
    # one of them is read is the D-106 failure, and it cost isobutanol a false floor miss.
    primary_share = ehrlich_primary_share(params, precursor)
    # ×n/(n+1) removes the D-106 decarboxylation CO₂: a full mole of precursor per mole of alcohol.
    total = (
        primary_share
        * consumed(precursor)
        * carbon_mass_fraction(precursor)
        * n_alc
        / (n_alc + 1.0)
    )
    if pool == "isoamyl_alcohol":
        total += (
            share_val_iso
            * consumed("valine")
            * carbon_mass_fraction("valine")
            * n_alc
            / CARBON_ATOMS["valine"]
        )
    return float(total / alcohol_carbon)


@pytest.fixture(scope="module")
def shipped_run():
    """The dosed characterization run, other precursor consumers off — computed ONCE (D-245).

    Three tests below read the same trajectory. Before D-245 they each paid their own
    ``solve_ivp``; sharing it is what makes the per-alcohol split below cost nothing extra.
    Module-scoped rather than session-scoped because ``--dist load`` hands a worker a contiguous
    chunk, so one module's fixture is paid once per run (see ``CLAUDE.md`` on why not worksteal).
    """
    return _run(aging=False, drop=_OTHER_PRECURSOR_CONSUMERS)


#: **D-120'S DIRECTION LEG IS BACK, AND THE MARKS ARE GONE (decision D-248).** D-245 measured
#: isoamyl and isobutanol crossing from under- to over-attribution (1.015x and 1.079x Minebois's
#: in-study shares) and xfailed them strictly; D-246 scored both on her OWN must and they roughly
#: doubled, to 1.73x and 1.68x. Un-coupling assimilable-nitrogen uptake from growth demand puts
#: them at **0.873x and 0.927x here, and 1.01x / 0.98x on her own must** — under the measurement
#: on this fixture, and essentially ON it where the comparison is commensurate. So every alcohol
#: is once again at or below Minebois, which is precisely the condition D-120 refused a de-novo
#: cap on: a ``(1-f_de_novo)`` ceiling can only REDUCE amino-acid sourcing and would move all
#: three the wrong way. Nothing was fitted to her: the uptake parameter sits at the bound where
#: transport stops limiting and every number here is unmoved across a 200x sweep of it.
#:
#: What is NOT fully restored is D-120's second (instrument) leg — see
#: :func:`test_the_de_novo_cap_now_bites_because_phenylalanine_no_longer_exhausts`, where the
#: cap's bite shrinks 12.7 % -> 4.8 % without reaching the inertness D-120 measured. One leg
#: back and one leg thinner is the honest statement, and it is why the de-novo-entry build is
#: not simply re-refused here.
#:
#: The ``_OVER_ATTRIBUTING`` tuple that used to select those marks is **deleted, not emptied**. An
#: empty one would read to the next person as a live switch that currently selects nothing — the
#: guards-that-quietly-forbid-nothing shape this suite exists to avoid — where the truth is that
#: no alcohol needs a mark at all.


@pytest.mark.parametrize(
    ("pool", "measured"),
    [
        pytest.param(pool, measured, id=pool)
        for pool, measured in _MINEBOIS_AMINO_ACID_SHARE.items()
    ],
)
def test_no_alcohol_over_attributes_to_amino_acids_so_no_de_novo_cap_is_warranted(
    pool, measured, shipped_run
):
    """D-118's "class of error" hypothesis, tested across all three — and REFUTED (decision D-120).

    **SPLIT PER ALCOHOL AT D-245, AND THE SPLIT IS THE POINT.** As one test with a loop this read
    ``2-PE → isoamyl → isobutanol`` and died on isoamyl, so the two alcohols after it were never
    evaluated — including **isobutanol, which is the leg D-245's record calls the real flip**
    (9.47 % against Minebois's 8.78 %, and the one that survives the cap-window correction). The
    leg that fired was the one the same record calls noise (isoamyl 5.42 vs 5.34 %, inside the
    harness's own systematic). A finding sitting behind a failing assert is exactly what this beat
    repaired twice elsewhere; it was live a third time in its own output
    [[feedback-an-xfail-buries-the-asserts-after-it]]. Now each alcohol carries its own mark:
    2-phenylethanol is a **green** guard again, and the two that trip say so independently.

    D-118 capped 2-phenylethanol's amino-acid sourcing and left "isoamyl's de-novo entry" as the
    next build, on the premise that all three alcohols are de-novo dominated in-study so the error
    is a *class*. **Measured, none of the three warrants the cap, and the sign is the reason.**

    :class:`~fermentation.core.kinetics.carbon_routing.DeNovoFuselRoute` is applied as
    ``gate *= (1 − f_de_novo)``: a ONE-DIRECTIONAL ceiling that can only ever REDUCE amino-acid
    sourcing. It is therefore warranted **iff the model over-attributes** to amino acids. Against
    Minebois's own in-study shares every alcohol sits on the *other* side — the model attributes
    LESS to amino acids than she measures — so the cap would drive all three further from the
    measurement it would be sourced from.

    This is the tripwire for the refused build: if a future beat makes any alcohol over-attribute,
    the de-novo-entry question genuinely re-opens for that alcohol and this fails.

    **IT FIRED, AS WRITTEN, AT D-244 — AND D-245 MEASURED WHICH ALCOHOLS AND BY HOW MUCH.**
    Isoamyl reads 5.42 % against Minebois's 5.34 % and isobutanol 9.47 % against her 8.78 %; 2-PE
    is still well under. Two things must be said about the size before anyone builds on it. The
    **isoamyl trip is inside this harness's own systematic**: the ``ehrlich_draws`` headroom cap
    binds over the first ~1.35 h (at most 12.8 % of the valine pool), during which valine's
    isoamyl branch is roughly half its designed 0.23 share, so the honest estimate is ~5.30 % —
    *under* the target, on a 1.5 % relative trip read off figure bars with mixed per-bar semantics
    (D-119). **Isobutanol's is the one that survives that correction** — the cap moves its share
    the other way, to ~10.4 % worst case, so the direction flip is real for it. Both are scored on
    a must carrying 1.0 g/L of amino acids against Minebois's own stock, which is the same
    commensurability caveat D-244 section 6 recorded and declined to repair.

    **THIS IS A CLAIM ABOUT THE D-109 FIXTURE, NOT ABOUT MINEBOIS'S MUST — and D-254 measured the
    difference rather than leaving it as a caveat.** Migrating this guard onto her own must is
    the move D-246 §7 left to the owner. It cannot be made: there, one of the two valine-derived
    alcohols over-attributes whichever estimator is used, and the two estimators name *different*
    alcohols — ``_amino_acid_share`` as written says isoamyl at 1.014×, and the same quantity
    measured from the draws the run actually applied says isobutanol at 1.049×. Relaxing the
    assertion to fit either one would be fitting a threshold to an outcome.

    On THIS fixture every alcohol is under both ways (isoamyl 0.873× / 0.854×, isobutanol
    0.927× / 0.977×), so the guard is not living on a measurement error and correcting the
    estimator in place would leave it green. See ``tests/test_fusel_provenance_estimator.py``,
    which pins both halves and holds the reasoning. **Do not migrate this to a commensurate must
    without re-deciding D-120's refusal**; and do not read the caveat above as merely unrepaired,
    because what is unrepaired is now identified.
    """
    traj, schema = shipped_run
    params = compile_scenario(_scenario(aging=False)).param_values

    model = _amino_acid_share(traj, schema, params, pool)
    assert 0.0 < model < measured, (
        f"{pool}: model sources {model:.2%} of it from amino acids against Minebois's "
        f"{measured:.2%}. D-120 refused a de-novo cap for isoamyl because a (1−f_de_novo) "
        "ceiling can only REDUCE amino-acid sourcing and every alcohol was already under. If "
        "this one now over-attributes, that refusal must be re-derived for it"
    )


@pytest.mark.xfail(strict=True, reason=_D245_D120_LEGS_GONE_DIRECTION_BACK_AT_D248)
def test_the_de_novo_cap_is_inert_where_the_precursor_exhausts(shipped_run):
    """Why the cap is the wrong INSTRUMENT, not merely the wrong direction (decision D-120).

    The cap multiplies the availability gate, i.e. a RATE. Where a precursor is fully consumed the
    total drawn is fixed by SUPPLY, so the realised amino-acid share sits on the D-112
    ``(1−f) × pool / alcohol`` mass-conservation ceiling and no rate knob can move it. D-112 found
    exactly this for leucine ("a gate cap does not move the realised share"); it generalises.

    **Measured on the SHIPPED parameter:** at the characterization must, phenylalanine exhausts
    with or without the cap, so ``f_de_novo_2_phenylethanol`` moves 2-PE's realised share by ~0.

    This does NOT make D-118 wrong, and the distinction is the point. The route's *measured* job is
    the instantaneous carbon-refund guard — pinned by its own counterfactual in
    ``test_the_de_novo_route_is_what_makes_the_sourced_lump_shippable`` — which is a different
    quantity from the integrated realised share. What does not survive is reading the route as the
    fix for the 18.9% over-attribution: at this dose that correction was delivered by
    ``f_non_ehrlich_phenylalanine`` (0.53 → 0.975), and the route is inert beside it.

    **THE PREMISE IS NOW FALSE — MEASURED AT D-245, NOT INFERRED — AND THIS TEST IS NOW THAT
    PREMISE ALONE.** D-244's corrected yield halves the biomass, so 2-PE's draw falls far enough
    that phenylalanine stops exhausting: **12.8 % of the pool survives**. The supply-limited regime
    this test names is simply gone, so its claim — "a rate knob cannot move a supply-limited
    share" — is true and no longer *applies*, which is why this is an xfail and not a re-pin.

    **Its second assertion has MOVED rather than died**, into
    :func:`test_the_de_novo_cap_now_bites_because_phenylalanine_no_longer_exhausts`, because what
    it now measures is a green fact about the shipped model and the receipt for D-245's
    ``Flags: D-120``. Left here it would have sat behind this failing premise, never running —
    which is the defect D-245 repaired twice elsewhere and shipped a third time in its own output
    [[feedback-an-xfail-buries-the-asserts-after-it]].
    """
    shipped, schema = shipped_run
    initial = float(shipped.y[schema.slice("phenylalanine"), 0][0])
    left = _end(shipped, schema, "phenylalanine") / initial
    assert left < 0.01, (
        f"phenylalanine no longer exhausts ({left:.1%} left) — the supply-limited premise of this "
        "test has moved, and the cap may now bite on the realised share"
    )


def test_the_de_novo_cap_now_bites_because_phenylalanine_no_longer_exhausts(shipped_run):
    """The cap MOVES 2-PE's realised share — D-120's instrument leg, inverted and pinned (D-245).

    **This is a re-record, not a new claim.** D-120 measured ``|with − without| < 1e-3 × without``
    and read it as proof that a de-novo cap is the wrong *instrument*: it multiplies a rate, and
    where a precursor is fully consumed the realised share sits on a mass-conservation ceiling no
    rate can move. That assertion lived inside
    :func:`test_the_de_novo_cap_is_inert_where_the_precursor_exhausts`, one line below an
    exhaustion premise that D-244 falsified — so from D-244 until D-245 **it never ran**, and the
    inversion below went unmeasured while the record above it was being written.

    Measured at D-245: phenylalanine survived at 12.8 % with the cap and exhausted without it, and
    the cap moved the realised amino-acid share **1.603 % → 1.400 %, 12.7 % relative** — four
    orders of magnitude past the old inertness threshold. Pinned two-sided, because both directions
    matter: shrinking toward 0 means the supply-limited regime has returned and D-120's instrument
    argument revives; growing means the cap is taking over the sourcing this fixture measures.

    **D-248 MOVED IT MOST OF THE WAY BACK, AND THE BAND IS RE-RECORDED RATHER THAN WIDENED.**
    Un-coupling assimilable-nitrogen uptake from growth demand restores the biomass D-244's
    corrected yield had halved, so 2-PE's draw rises again and phenylalanine goes from 12.8 % left
    to **4.85 %** — much nearer exhaustion, though not at it. The cap's bite falls with it:
    **1.3790 % → 1.3125 %, 4.82 % relative**. So the direction is *toward* D-120's inertness
    without reaching it, which is why
    :func:`test_the_de_novo_cap_is_inert_where_the_precursor_exhausts` is still a strict xfail and
    this test is still a live guard rather than being retired into it.

    **This does not license building a cap for another alcohol, and D-248 makes that stronger
    rather than weaker.** D-245 removed both of D-120's measured legs; D-248 gives the DIRECTION
    leg back outright — every alcohol is at or below Minebois again — and leaves this INSTRUMENT
    leg thinner than D-120 measured but no longer four orders of magnitude away from it. The build
    would still need a sourced ``f_de_novo_isoamyl``, which this repo does not hold (D-206).
    """
    shipped, schema = shipped_run
    uncapped, _ = _run(
        aging=False,
        drop=_OTHER_PRECURSOR_CONSUMERS,
        set_params={"f_de_novo_2_phenylethanol": 0.0},
    )
    params = compile_scenario(_scenario(aging=False)).param_values

    # Anti-vacuity, and the mechanism in one line: the cap is what leaves phenylalanine standing.
    initial = float(shipped.y[schema.slice("phenylalanine"), 0][0])
    left_uncapped = _end(uncapped, schema, "phenylalanine") / initial
    assert left_uncapped < 0.01, (
        f"phenylalanine no longer exhausts even with the cap OFF ({left_uncapped:.2%} left) — the "
        "counterfactual this test rests on has moved, so the bite below is not the cap's doing"
    )

    with_cap = _amino_acid_share(shipped, schema, params, "2_phenylethanol")
    without = _amino_acid_share(
        uncapped, schema, {**params, "f_de_novo_2_phenylethanol": 0.0}, "2_phenylethanol"
    )
    move = abs(with_cap - without) / without
    assert 0.030 < move < 0.070, (
        f"the de-novo cap moves 2-PE's realised amino-acid share {without:.4%} → {with_cap:.4%} "
        f"({move:.2%} relative), outside the D-248 band [3.0 %, 7.0 %] (measured 4.82 %; D-245 "
        "measured 12.7 % before nitrogen uptake was un-coupled). Falling toward 0 means the "
        "supply-limited regime has fully returned, D-120's instrument argument is whole again, "
        "and the xfail above should be lifted; rising means the biomass denominator has moved "
        "back. Either way this is a re-derivation, not a band to widen"
    )
