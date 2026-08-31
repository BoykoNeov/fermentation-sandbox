"""What ``_amino_acid_share`` actually measures, and the fate split a backstop is setting.

**The beat (decision D-254).** D-246 §7 left "migrate the fusel guards onto the papers' own
musts" as the owner's call, and D-249/D-250 carried the note forward unchanged. Taking it turned
up two things that are not about test fixtures at all.

**1. The reason the migration was blocked had expired.** ``test_fusel_keto_acid_node._scenario``
and ``test_fusel_reroute._rollero_run`` carried the same verbatim comment: the mismatch is
"recorded and NOT repaired" because the paper's medium is "not in this repo … closing the gap
needs the paper". D-246 sourced Crépin's and Minebois's medium (Bely, Sablayrolles & Barre 1990)
and put both musts in this repo. Rollero's is genuinely still missing and its comment says so
now; Crépin's was a standing instruction not to try. Eighth instance in this archive of a
sourcing blocker outliving its blocker.

**2. ``_amino_acid_share`` is a BIASED estimator, in opposite directions for the two alcohols
built from valine — and the bias is a conservation backstop doing physiology.**

The helper credits each valine branch at its *designed* share of consumed valine:
``f_valine_to_isoamyl`` (0.23) for the secondary branch, ``1 − f − share`` (0.15) for
isobutanol's primary. :func:`~fermentation.core.kinetics.byproducts.ehrlich_draws` truncates the
secondary branch at ``headroom`` — an alcohol cannot be sourced from more precursor than it is
being made from — so the *realised* split is not the designed one. Because the non-Ehrlich lump
scales against the Ehrlich draw, the TOTAL Ehrlich draw on valine is pinned at ``1 − f``
whatever the truncation does, so every gram it takes off isoamyl is absorbed by isobutanol.

**The truncation binding is not the new part, and this file must not be read as discovering it.**
D-111's Finding 5 measured exactly this when the clamp was added — "realised 0.2233 against a
sourced 0.23, isobutanol 0.1567 against 0.15" — and credited it with isoamyl's 1.80 → 1.74. What
is new is that ``_amino_acid_share``, written 130-odd records later at D-245, credits the
*designed* shares and so silently re-assumes what D-111 had already falsified; and that on the
papers' own musts, sourceable only since D-246, the resulting bias decides which alcohol
over-attributes. The comment on the truncation meanwhile said it "never binds in practice",
contradicting the record that shipped it.

**Why this decides the migration.** Against Minebois's own in-study shares, on her own must:

======================  ==========  ==============
quantity                 isoamyl     isobutanol
======================  ==========  ==============
harness (designed)       1.0143×      0.9819×
realised (measured)      0.9889×      **1.0489×**
======================  ==========  ==============

``test_fusel_catabolic_shape.test_no_alcohol_over_attributes_to_amino_acids_so_no_de_novo_cap_
is_warranted`` asserts ``0 < model < measured`` and is D-120's tripwire: a ``(1 − f_de_novo)``
cap can only REDUCE amino-acid sourcing, so it is warranted iff the model over-attributes.
**On the commensurate must one of the two valine alcohols over-attributes whichever estimator
is used** — the harness names isoamyl, the realised measurement names isobutanol and by more.
So the guard cannot be migrated, and it must not be relaxed to fit; it stays on the D-109
fixture as a claim about that fixture (where every alcohol is under, both ways), and the
commensurate statement stays with ``test_defined_media``'s two-sided band, now carrying an
explicit flag that its band contains an over-attribution.

**What is NOT claimed here.** The truncation is not itself wrong: early in the run isoamyl's
production is genuinely small while the secondary branch — anchored *algebraically* off
isobutanol's draw rather than gated on its own — asks for a share of a draw isoamyl is not yet
making. The incommensurate part is that anchoring. Identifying a replacement is a rate-law
question and is **refused here as unidentified**, not deferred as easy.

**The no-cap counterfactual is a probe, not a test**, because reproducing it needs a copy of
``ehrlich_draws`` with the truncation deleted, and a copy in a test file drifts from the
original. It is at ``M:\\claud_projects\\temp\\ferment\\d254-commensurate-guards\\drive_nocap.py``
(patch verified live: realised/credited becomes 1.000003 against the shipped 0.95536). Its
result is the row above: remove the truncation and the over-attribution moves back to isoamyl.
What licenses reading the harness value as that counterfactual is
:func:`test_the_cap_can_move_provenance_but_never_state`.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.chemistry import CARBON_ATOMS, carbon_mass_fraction
from fermentation.core.kinetics.byproducts import ehrlich_draws
from fermentation.core.kinetics.carbon_routing import (
    FUSEL_SPECS,
    SECONDARY_FUSEL_ROUTES,
    non_ehrlich_fraction_param,
)
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import compile_scenario
from tests.test_defined_media import commensurate_scenario
from tests.test_fusel_catabolic_shape import _MINEBOIS_AMINO_ACID_SHARE, _amino_acid_share
from tests.test_fusel_keto_acid_node import _OTHER_PRECURSOR_CONSUMERS, ehrlich_primary_share

ROUTE = next(r for r in SECONDARY_FUSEL_ROUTES if r.precursor == "valine")

#: Points in the ``t_eval`` grid the realised draw is integrated over. **Load-bearing, in the
#: D-253 ``SHIPPED_R`` sense**: the pinned biases below are quadrature results, so the grid is
#: part of what they mean and cannot be changed without re-measuring them.
#: :func:`test_the_integrated_bias_is_grid_converged` is the receipt that it is dense enough.
GRID = 800
#: The refinement arm. 4× the points; the bias must agree to :data:`GRID_TOL`.
FINE_GRID = 3200
GRID_TOL = 5e-4


def _run(which: str, n: int):
    cs = compile_scenario(commensurate_scenario(which))
    for name in _OTHER_PRECURSOR_CONSUMERS:
        cs.process_set.disable(name)  # KeyErrors on a rename rather than silently no-op
    traj = simulate_scheduled(
        cs.process_set,
        cs.param_values,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
        t_eval=np.linspace(cs.t_span_h[0], cs.t_span_h[1], n),
    )
    assert traj.success, traj.message
    return traj, cs.schema, cs.param_values


@pytest.fixture(scope="module")
def minebois():
    return _run("minebois", GRID)


@pytest.fixture(scope="module")
def crepin():
    return _run("crepin", GRID)


def _consumed(traj, schema, species: str) -> float:
    return float(traj.y[schema.slice(species), 0][0]) - float(traj.y[schema.slice(species), -1][0])


def _designed_share(params, pool: str) -> float:
    share = float(params[ROUTE.share_param])
    if pool == "isoamyl_alcohol":
        return share
    return 1.0 - float(params[non_ehrlich_fraction_param("valine")]) - share


def realised_branch_carbon(traj, schema, params, pool: str) -> float:
    """Alcohol carbon this run actually drew from **valine** for ``pool``, integrated.

    Quadrature, and it has to be: valine is fully consumed either way, so the total is fixed by
    supply (D-112) and the realised *split* between the two branches lives in no state slot. The
    counterfactual route — toggle ``f_valine_to_isoamyl`` to 0 and diff the state — measures
    nothing for the same reason, which is why it is not used. The quadrature is checked two
    ways: against a refined grid (:func:`test_the_integrated_bias_is_grid_converged`) and against
    an exact state difference (:func:`test_the_quadrature_reproduces_an_exact_state_difference`).
    """
    rates = [
        sum(
            d.alcohol_carbon
            for d in ehrlich_draws(traj.y[:, i], schema, params)
            if d.precursor.species == "valine" and d.alcohol.pool == pool
        )
        for i in range(traj.t.size)
    ]
    return float(np.trapezoid(rates, traj.t))


def credited_branch_carbon(traj, schema, params, pool: str) -> float:
    """What ``_amino_acid_share`` credits the same branch: the DESIGNED share of consumed valine."""
    spec = next(s for s in FUSEL_SPECS if s.pool == pool)
    valine_carbon = _consumed(traj, schema, "valine") * carbon_mass_fraction("valine")
    return (
        _designed_share(params, pool)
        * valine_carbon
        * CARBON_ATOMS[spec.species]
        / CARBON_ATOMS["valine"]
    )


def corrected_share(traj, schema, params, pool: str) -> float:
    """``_amino_acid_share`` with the valine branch replaced by what was actually drawn.

    Leucine's primary branch rides along unchanged: it is gated, not anchored, so nothing
    truncates it and the designed-share credit is exact for it.
    """
    spec = next(s for s in FUSEL_SPECS if s.pool == pool)
    n_alc = CARBON_ATOMS[spec.species]
    alcohol_carbon = float(traj.y[schema.slice(pool), -1][0]) * carbon_mass_fraction(spec.species)
    total = realised_branch_carbon(traj, schema, params, pool)
    if pool == "isoamyl_alcohol":
        total += (
            ehrlich_primary_share(params, "leucine")
            * _consumed(traj, schema, "leucine")
            * carbon_mass_fraction("leucine")
            * n_alc
            / (n_alc + 1.0)
        )
    return total / alcohol_carbon


# -- the measurement machinery, before anything is read off it ---------------------------------


def test_the_quadrature_reproduces_an_exact_state_difference(minebois):
    """The integration is faithful: integrated Ehrlich valine carbon ÷ exact consumed = ``1 − f``.

    **This checks the QUADRATURE, not the truncation.** The equality is an identity in the model
    — ``PrecursorNonEhrlichFates`` scales the non-Ehrlich lump against the Ehrlich draw, so the
    Ehrlich share of consumed valine is ``1 − f`` however the branches split underneath. That is
    exactly what makes it a good check of the numerics and a useless one for the bias: the
    denominator is an exact state difference and the numerator is the quadrature, so if the
    integration were coarse the ratio would miss ``1 − f``. It does not.

    It is also the reason the bias has to be equal and opposite between the two alcohols rather
    than a net loss: the total is pinned, so the truncation can only MOVE carbon between them.
    """
    traj, schema, params = minebois
    integrated = float(
        np.trapezoid(
            [
                sum(
                    d.precursor_carbon
                    for d in ehrlich_draws(traj.y[:, i], schema, params)
                    if d.precursor.species == "valine"
                )
                for i in range(traj.t.size)
            ],
            traj.t,
        )
    )
    exact = _consumed(traj, schema, "valine") * carbon_mass_fraction("valine")
    expected = 1.0 - params[non_ehrlich_fraction_param("valine")]
    assert integrated / exact == pytest.approx(expected, rel=1e-3), (
        f"integrated Ehrlich valine carbon is {integrated / exact:.5f} of the exact consumed "
        f"total against the identity's {expected:.5f} — the quadrature in this file is too "
        f"coarse to be read, and every bias below is suspect"
    )


def test_the_integrated_bias_is_grid_converged():
    """The pinned biases do not move on a 4× finer grid.

    Without this the numbers below would be quadrature artefacts wearing a decimal point.
    ``GRID`` is load-bearing and this is what says so.
    """
    c_traj, c_schema, c_params = _run("minebois", GRID)
    f_traj, f_schema, f_params = _run("minebois", FINE_GRID)
    for pool in ("isoamyl_alcohol", "isobutanol"):
        c = realised_branch_carbon(c_traj, c_schema, c_params, pool) / credited_branch_carbon(
            c_traj, c_schema, c_params, pool
        )
        f = realised_branch_carbon(f_traj, f_schema, f_params, pool) / credited_branch_carbon(
            f_traj, f_schema, f_params, pool
        )
        assert abs(c - f) <= GRID_TOL, (
            f"{pool}: bias reads {c:.5f} at n={GRID} and {f:.5f} at n={FINE_GRID}; the grid is "
            f"not dense enough for the pins in this file"
        )


def test_the_cap_can_move_provenance_but_never_state(minebois):
    """Valine exhausts and each alcohol's amount is its own rate law's — so state is cap-blind.

    **This is what licenses reading ``_amino_acid_share`` as the no-cap counterfactual.** The
    truncation changes how much valine the secondary branch draws *per instant*; it cannot change
    how much valine is consumed in total (the pool goes to zero either way — supply-limited,
    D-112's own finding) and it cannot change how much alcohol is made (the reroute sets carbon
    PROVENANCE; the rate law sets the amount). Driven here by zeroing the secondary route
    entirely, which is a far larger perturbation than the truncation ever applies: nothing in the
    state moves.

    So the cap is invisible to every guard in this repo that reads a concentration, and only a
    provenance measurement can see it at all. That is why it survived to D-254 with a comment
    saying it never fires.
    """
    on_traj, on_schema, _ = minebois
    cs = compile_scenario(commensurate_scenario("minebois"))
    for name in _OTHER_PRECURSOR_CONSUMERS:
        cs.process_set.disable(name)
    params = dict(cs.param_values)
    params[ROUTE.share_param] = 0.0
    off_traj = simulate_scheduled(
        cs.process_set,
        params,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
        t_eval=np.linspace(cs.t_span_h[0], cs.t_span_h[1], GRID),
    )
    assert off_traj.success, off_traj.message
    off_schema = cs.schema

    assert float(on_traj.y[on_schema.slice("valine"), -1][0]) < 1e-9, (
        "valine no longer exhausts on Minebois's must — the supply-limited premise of this "
        "file's whole argument has moved and every bias in it needs re-deriving"
    )
    # Valine is compared as CONSUMED, not as residual: it exhausts, so its residual is numerical
    # zero (~1e-17 here, ~2e-14 with the route off) and a relative bound on that compares two
    # noise figures rather than the quantity the argument uses. Consumed is ~0.031 g/L and is
    # what "supply-limited" is a claim about.
    quantities = {
        "valine consumed": (
            _consumed(on_traj, on_schema, "valine"),
            _consumed(off_traj, off_schema, "valine"),
        ),
        "isoamyl_alcohol": (
            float(on_traj.y[on_schema.slice("isoamyl_alcohol"), -1][0]),
            float(off_traj.y[off_schema.slice("isoamyl_alcohol"), -1][0]),
        ),
        "isobutanol": (
            float(on_traj.y[on_schema.slice("isobutanol"), -1][0]),
            float(off_traj.y[off_schema.slice("isobutanol"), -1][0]),
        ),
    }
    # **Not exactly zero, and the size is the claim.** The reroute REFUNDS the sourced carbon to
    # sugar, so deleting the route perturbs ``S`` a little and that feeds every rate: isoamyl
    # moves 3.4e-6 relative. What matters is the ratio to the provenance effect it must not be
    # confused with — 4.5 % — which is four orders larger. Asserting exact invariance here would
    # be false; asserting the separation is what the no-cap reading actually needs.
    for name, (a, b) in quantities.items():
        assert a > 1e-6, f"{name} is at noise scale ({a:.3e}); this comparison would be vacuous"
        moved = abs(a - b) / a
        assert moved <= 1e-5, (
            f"{name} moved {a:.8e} -> {b:.8e} ({moved:.2e} relative) when the secondary route "
            "was switched off. The state is meant to be blind to this route's share to within "
            "the sugar refund's own footprint; a larger move means the truncation has a real "
            "state footprint and this file's no-cap reading is void"
        )
        assert moved < 0.045 / 100.0, (
            f"{name}'s state footprint ({moved:.2e}) is within a hundredth of the 4.5 % "
            "provenance shift it is supposed to be negligible against"
        )


# -- the bias itself ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("which", "pool", "expected"),
    [
        ("minebois", "isoamyl_alcohol", 0.95536),
        ("minebois", "isobutanol", 1.06827),
        ("crepin", "isoamyl_alcohol", 0.97332),
        ("crepin", "isobutanol", 1.04069),
    ],
    ids=lambda v: v if isinstance(v, str) else f"{v}",
)
def test_the_estimator_is_biased_in_opposite_directions_for_the_two_valine_alcohols(
    which, pool, expected, request
):
    """Realised ÷ credited, per branch, on both papers' own musts (decision D-254).

    Under 1 for isoamyl and over 1 for isobutanol, on both musts, because the total is pinned
    (see the quadrature test) — the truncation moves carbon from one to the other rather than
    losing it. Minebois's must is the sharper case: 4.5 % of the secondary branch.

    **The direction is the whole point, the third decimal is not.** If a future beat changes the
    anchoring these numbers move, and that is fine; what must not happen silently is the bias
    going to zero (the truncation stopped binding — then say so and delete this file's premise)
    or reversing sign (isoamyl gaining at isobutanol's expense, which no branch here can do).
    """
    traj, schema, params = request.getfixturevalue(which)
    bias = realised_branch_carbon(traj, schema, params, pool) / credited_branch_carbon(
        traj, schema, params, pool
    )
    assert bias == pytest.approx(expected, abs=2e-3), (
        f"{which}/{pool}: `_amino_acid_share` credits this branch {1 / bias:.4f}x what the run "
        f"actually drew (bias {bias:.5f}, was {expected:.5f})"
    )
    if pool == "isoamyl_alcohol":
        assert bias < 1.0, "the truncation can only ever REDUCE the secondary branch"
    else:
        assert bias > 1.0, "isobutanol's primary must absorb what the truncation takes"


def test_the_headroom_cap_binds_where_its_comment_said_it_never_would(minebois):
    """Falsifies the prose that stood in ``ehrlich_draws`` from D-111 to D-254.

    It read: "Measured headroom is large (valine ~1.8% + leucine ~1.1% of isoamyl), so this never
    binds in practice". Asserted here on the consequence rather than by re-deriving ``headroom``
    in a test: if the truncation never bound, realised would equal credited and the bias above
    would be exactly 1. It is 0.955.

    Kept separate from the parametrized bias test because it is a different claim — that one
    pins a magnitude, this one refuses a sentence — and a magnitude pin drifting into a tolerance
    would silently take the refusal with it.
    """
    traj, schema, params = minebois
    bias = realised_branch_carbon(traj, schema, params, "isoamyl_alcohol") / (
        credited_branch_carbon(traj, schema, params, "isoamyl_alcohol")
    )
    assert bias < 0.99, (
        f"the headroom truncation appears not to bind (bias {bias:.5f}). If that is real, the "
        "D-254 finding has reversed: say so, and re-derive both the estimator correction and "
        "the fate split it was setting — do not just widen this bound"
    )


# -- what it means for the guard that wanted to migrate -----------------------------------------


def test_one_valine_alcohol_over_attributes_on_mineboiss_own_must_either_way(minebois):
    """Why D-120's tripwire cannot migrate onto the commensurate must (decision D-254).

    The tripwire asserts ``0 < model < measured`` for every alcohol, because a ``(1 − f_de_novo)``
    cap can only REDUCE amino-acid sourcing and is therefore warranted iff the model
    over-attributes. On the D-109 fixture every alcohol is under, both ways, and it passes there.

    On Minebois's own must it fails whichever estimator is used, and they disagree about which
    alcohol: the harness names isoamyl at 1.014×, the realised measurement names isobutanol at
    1.049×. **That disagreement is why this is a finding and not a re-pin.** Migrating the guard
    and relaxing the assertion to fit would be fitting a threshold to an outcome; migrating it
    without relaxing it would xfail a test whose stated premise ("no alcohol over-attributes") is
    simply not true on that must.

    So the guard stays where it is, as a claim about the D-109 fixture, and this records what the
    migration would have said. **The de-novo cap question is therefore OPEN on the commensurate
    must** — for isobutanol on the measurement this file argues is the right one, which has no
    de-novo route today (``DE_NOVO_FUSEL_ROUTES`` carries ``f_de_novo_2_phenylethanol`` alone).
    It is not built here: D-120's refusal was reasoned on the *direction* of a one-sided cap, and
    re-opening it needs a sourced ``f_de_novo`` for isobutanol, which nothing in this repo has.
    """
    traj, schema, params = minebois
    harness = {
        p: _amino_acid_share(traj, schema, params, p) / _MINEBOIS_AMINO_ACID_SHARE[p]
        for p in ("isoamyl_alcohol", "isobutanol")
    }
    realised = {
        p: corrected_share(traj, schema, params, p) / _MINEBOIS_AMINO_ACID_SHARE[p]
        for p in ("isoamyl_alcohol", "isobutanol")
    }
    assert harness["isoamyl_alcohol"] == pytest.approx(1.0143, abs=5e-3)
    assert harness["isobutanol"] == pytest.approx(0.9819, abs=5e-3)
    assert realised["isoamyl_alcohol"] == pytest.approx(0.9889, abs=5e-3)
    assert realised["isobutanol"] == pytest.approx(1.0489, abs=5e-3)

    over_harness = [p for p, r in harness.items() if r > 1.0]
    over_realised = [p for p, r in realised.items() if r > 1.0]
    assert over_harness and over_realised, (
        f"no alcohol over-attributes on Minebois's own must any more (harness {harness}, "
        f"realised {realised}). That would REMOVE the reason D-120's tripwire cannot migrate, "
        "and the migration D-246 §7 left open becomes takeable — a re-decision, not a fix"
    )
    assert over_harness != over_realised, (
        f"the two estimators now agree on which alcohol over-attributes ({over_harness}). The "
        "D-254 finding is that they disagree, which is what makes the de-novo cap question "
        "un-actionable rather than merely open; if they have converged, it is actionable"
    )


def test_the_d109_fixture_is_under_on_both_estimators_so_the_tripwire_is_safe_there():
    """The other half of the same claim: migrating is what breaks it, not the estimator fix.

    Every alcohol on the D-109 characterization must is under Minebois both ways (isoamyl 0.873×
    harness / 0.854× realised, isobutanol 0.927× / 0.977×). So the tripwire is not living on a
    measurement error — correcting the estimator in place would leave it green — and the record
    cannot be read as "the guard was wrong all along". It is the must that decides.
    """
    from tests.test_fusel_keto_acid_node import _scenario

    cs = compile_scenario(_scenario(aging=False))
    for name in _OTHER_PRECURSOR_CONSUMERS:
        cs.process_set.disable(name)
    traj = simulate_scheduled(
        cs.process_set,
        cs.param_values,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
        t_eval=np.linspace(cs.t_span_h[0], cs.t_span_h[1], GRID),
    )
    assert traj.success, traj.message
    for pool in ("isoamyl_alcohol", "isobutanol"):
        measured = _MINEBOIS_AMINO_ACID_SHARE[pool]
        assert _amino_acid_share(traj, cs.schema, cs.param_values, pool) < measured
        assert corrected_share(traj, cs.schema, cs.param_values, pool) < measured, (
            f"{pool} over-attributes on the D-109 fixture under the corrected estimator. That "
            "moves the tripwire's failure off the commensurate must and onto the fixture it "
            "actually runs on, which is a live defect rather than a recorded one"
        )
