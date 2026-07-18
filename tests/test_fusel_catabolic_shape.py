"""The fusel catabolic SHAPE, measured — the receipts for D-112.

D-111 left "the leucine-derived isoamyl shortfall (1.12% vs Rollero 3.4-17.3%)" as the sharpest
open item, framed as D-103's gate *shape* that the keto-acid node would fix. **D-112 measured it
and retired that framing.** This suite pins the three measurements so a future beat cannot quietly
re-inherit the stale story:

1. **D-103's gate-shape SPREAD is absorbed by the D-104 non-Ehrlich sink** — it is large with the
   sink off (isoamyl ~6% vs propanol ~67%) and compresses to a uniform low band with it on. So the
   "minor alcohols are wildly over-attributed" defect is gone, and only isoamyl (UNDER) survives.
2. **Isoamyl sits on its ``(1-f)`` mass-conservation ceiling**, which no sourcing-layer change —
   the keto-acid node included — can lift: a gate cap (the "obvious" fix) does not move it, because
   leucine is too scarce to persist under any draw rate.
3. **``k_isoamyl_alcohol`` is correctly calibrated** to the Wang 2024 172 mg/L anchor at typical
   must nitrogen with no amino-acid dose — so the ~2x isoamyl over-production in the D-109
   characterization must is the ``amino_acids_gpl=1.0`` dose (Finding 4), not a mis-set ``k``.

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
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS, ISOAMYL_ALCOHOL
from fermentation.core.state import FloatArray, StateSchema
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from tests.test_fusel_keto_acid_node import (
    _FERMENT_DAYS,
    _OTHER_PRECURSOR_CONSUMERS,
    _de_novo_share,
    _end,
    _run,
    _scenario,
)

_SINK = "precursor_non_ehrlich_fates"
_PROPANOL = next(s for s in FUSEL_SPECS if s.pool == "propanol")


def _direct_catabolic_share_sink_off(spec) -> float:
    """Precursor-derived carbon fraction of ``spec``'s alcohol with the D-104 sink DISABLED.

    With the sink off the re-route is the precursor's ONLY consumer, so consumed precursor carbon
    times ``n/(n+1)`` (removing the D-106 decarboxylation CO2) is exactly the alcohol carbon sourced
    from it — an exact state difference, no quadrature. This is ~D-103's gate value; contrast with
    :func:`_de_novo_share` (sink on), which yields the compressed ``(1-f)``-scaled share.
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
    throttles its own gate). With the sink ON the ``(1-f)`` multiplier compresses every alcohol into
    Rollero's uniform low band. So "minor alcohols are wildly over-attributed" is no longer true,
    and the only survivor is isoamyl — which is UNDER, the opposite direction. If a future beat
    re-widens the sink-on spread, the D-112 premise (the gate-shape defect is retired) has moved.
    """
    iso_off = _direct_catabolic_share_sink_off(ISOAMYL_ALCOHOL)
    prop_off = _direct_catabolic_share_sink_off(_PROPANOL)
    # D-103's spread: sink OFF, propanol reads several times isoamyl (a minor alcohol over-read).
    assert prop_off > 4.0 * iso_off, (
        f"sink-OFF spread collapsed (propanol {prop_off:.1%} vs isoamyl {iso_off:.1%}) — D-103's "
        "gate-shape spread is the thing D-112 says the sink absorbs; if it is gone with the sink "
        "off too, finding 1's baseline has moved"
    )
    # Sink ON: every alcohol compresses to the low band (all five <= ~17%, Rollero's uniform range).
    on = {s.pool: 1.0 - _de_novo_share(s) for s in FUSEL_SPECS}
    assert max(on.values()) < 0.20, f"sink-on shares no longer compressed to the low band: {on}"
    # And the compression is real: propanol specifically drops from tens-of-% to the low band.
    assert on["propanol"] < 0.5 * prop_off


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
    reproduce that mean — so the ~2x isoamyl over-production in the D-109 characterization must
    (307 mg/L) is the ``amino_acids_gpl=1.0`` dose's deamination-N sustaining the fusel gate, NOT a
    mis-set ``k``. This is what makes the isoamyl catabolic DENOMINATOR a probe artifact, not a
    calibration bug: the ceiling comparison in finding 3 rests on the ``k`` being right here.
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
