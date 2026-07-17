"""The FUSEL side of the keto-acid node — the findings D-109 measured instead of building.

D-107 built the excreted ``alpha_ketobutyrate`` pool and left a work-list item that D-108 promoted
to "the largest open item": re-base propanol on that pool, so the *genuine propanol-vs-sotolon
competition* — over α-ketobutyrate, which propanol IS the decarboxylation of — becomes expressible.
**D-109 measured that item's premise before building it, and the premise is wrong twice over.**
This suite pins the measurements, because every one of them is a claim a future beat could quietly
invalidate.

**1. The pool D-107 chose FOR sotolon is the wrong pool for propanol — D-49's test, applied
symmetrically.** :class:`~fermentation.core.kinetics.aging.SotolonAldolCondensation` selects the
*excreted, extracellular residual* precisely because its aldol runs in a sealed bottle where no
intracellular pool can reach. Propanol is made **intracellularly, by living yeast, during active
fermentation** — it is exactly D-49's *flux intermediate*, the thing that module rejected as
acetaldehyde's precursor for the same reason. One pool cannot be both. The principle convicts on
its own (no literature needed — the sibling Processes already carry the argument), and
:func:`test_the_excreted_pool_cannot_supply_propanol` adds the arithmetic: the pool is not merely
the wrong *identity*, it is the wrong *size* by ~3×, so the re-base is **infeasible**, not just
mis-attributed.

**2. The work-list's "the competition is GONE" is CORRECT, and D-109 nearly "corrected" it.**
:class:`~fermentation.core.kinetics.keto_acids.AlphaKetobutyrateExcretion`'s rate is flux-only;
threonine's gate re-routes the *carbon source*, not the *rate*. So the competition over threonine
is **exactly zero, by construction and on purpose** (gating the rate on threonine would kill
sotolon in a threonine-free wine — the D-104 canary). D-109's first probe measured a 0.42%,
monotone, correctly-signed sotolon response to ``k_propanol`` and read it as "the competition is
present but small". It is the **sugar ledger**: propanol's de-novo carbon leaves ``S``, which every
flux-linked rate reads. See :func:`test_alpha_kb_production_is_exactly_threonine_independent` —
the mutation test is the only reason that write-up did not ship.

**3. The promised payoff is REAL, and the number that kills the shortcut is the number that says
so.** Propanol's molar demand is **~2.8× the total α-KB the pool ever excretes** — which is why it
cannot be drawn from that pool (finding 1), and equally why a *correctly placed* node would matter:
propanol is ~73% of the 2-ketobutyrate flux, so it is the node's dominant sink, and partitioning
that flux honestly would couple propanol and sotolon **materially**. The item is not dissolved by
this beat; it is **relocated** — off the excreted pool and onto the intracellular partition, where
it needs the milestone rather than a shortcut. (The 73% rests on ``k_alpha_kb_excretion``, an
author estimate — so it is an order-of-magnitude claim, not a calibration.)

**Why the fusel-side node is a PARTITION and not a pool (the scoping result).** D-49's physics
says the intracellular keto acid is a vanishing pool carrying an enormous flux — i.e.
quasi-steady, ``synthesis == Σ consumption`` at every instant. A quasi-steady node is a **flux
partition**, not a state variable. So the two nodes are different in *kind*: the excreted
keto-acids are pools (they persist, they are measured, they bind SO₂), and the fusel-side node is
a partition of the sourcing. That is why the ``FuselAminoAcidReroute``/``PrecursorNonEhrlichFates``
sourcing layer — not the producer — is where the milestone belongs, and it is why nothing here
needs a sixth state slot.
"""

from __future__ import annotations

import pytest

from fermentation.core.chemistry import (
    CARBON_ATOMS,
    M_ALPHA_KETOBUTYRATE,
    M_PROPANOL,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS
from fermentation.core.kinetics.precursor_fates import non_ehrlich_fraction_param
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import Intervention, Scenario, TemperaturePoint, compile_scenario

#: Crépin *et al.* 2017's own must: 180 mg N/L, 28 °C. Every share below is quoted against that
#: paper's numbers, so the probe must run on a COMMENSURATE must — D-104's lesson, where a
#: ~470 mg N/L probe flattered the model against a 180 mg N/L source.
_CREPIN_YAN = 180.0
_CREPIN_TEMP = 28.0
_FERMENT_DAYS = 14.0
_AGING_DAYS = 720.0

#: Crépin's measured de-novo share of 2-ketobutyrate (the 81 of "19% exogenous / 81% newly
#: synthesised"), and Rollero's independent ">90% from the carbon central metabolism". The model
#: is asserted against the WEAKER of the two, so the test cannot be broken by the disagreement
#: between them — recorded as two bands, never averaged (D-103).
_SOURCED_DE_NOVO_FLOOR = 0.80

#: The routes that ALSO eat the speciated precursors. Disabled where the ``f : (1−f)`` split
#: invariant is used to size the Ehrlich draw from a state difference, so that invariant is
#: EXACT rather than approximately true (they measure inert in this regime anyway — but "measures
#: inert" is a reason to verify, not to assume: D-106).
_OTHER_PRECURSOR_CONSUMERS = (
    "alpha_kb_excretion",
    "maillard_strecker",
    "strecker_degradation",
    "autolytic_mercaptan",
)


def _scenario(*, o2_mgl: float = 20.0, aging: bool = True) -> Scenario:
    interventions = []
    duration = _FERMENT_DAYS
    if aging:
        interventions = [
            Intervention(day=_FERMENT_DAYS, action="begin_aging"),
            Intervention(day=_FERMENT_DAYS, action="add_oxygen", params={"o2_mgl": o2_mgl}),
        ]
        duration = _FERMENT_DAYS + _AGING_DAYS
    return Scenario(
        name="d109-fusel-node",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": _CREPIN_YAN,
            "pitch_gpl": 0.25,
            "amino_acids_gpl": 1.0,
        },
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=_CREPIN_TEMP),
            TemperaturePoint(day=_FERMENT_DAYS, celsius=18.0),
        ],
        interventions=interventions,
        duration_days=duration,
    )


def _run(*, drop: tuple[str, ...] = (), scale: dict[str, float] | None = None, **kw):
    cs = compile_scenario(_scenario(**kw))
    for name in drop:
        cs.process_set.disable(name)  # KeyErrors on a renamed Process rather than silently no-op
    # ``cs.param_values`` is a PROPERTY that re-resolves from ``parameters`` on every access, so
    # ``cs.param_values[k] = v`` mutates a throwaway dict and ``cs.run()`` re-resolves a clean one
    # — a knob that looks connected and is not. Resolve ONCE and hand that dict to the runtime,
    # exactly as ``CompiledScenario.run`` does. (This bit D-109's own first probe: it reported a
    # perfect "+0.000000%" agreement with the answer being tested for.)
    pv = cs.param_values
    for key, factor in (scale or {}).items():
        assert key in pv, f"no such parameter {key!r}"
        pv[key] = pv[key] * factor
    traj = simulate_scheduled(
        cs.process_set,
        pv,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
    )
    assert traj.success, traj.message
    return traj, cs.process_set.schema


def _end(traj, schema, name: str) -> float:
    return float(traj.y[schema.slice(name), -1][0])


# -- finding 1: the excreted pool cannot supply propanol (design A is INFEASIBLE) --------------


def test_the_excreted_pool_cannot_supply_propanol():
    """Propanol's molar demand exceeds every gram of α-KB the pool ever carries (decision D-109).

    The prescribed "re-base propanol on ``alpha_ketobutyrate``" is rejected on **principle** —
    D-49's flux-intermediate-vs-excreted-residual test, applied symmetrically to the pool D-107
    selected for sotolon on exactly that test. This is the **arithmetic** that corroborates it, and
    it is the sharper statement: propanol is 1 mol 2-KB decarboxylated, so drawing it from this
    pool needs 1 mol of pool per mol of alcohol — and the pool never holds a third of that.

    So the re-base does not merely mis-attribute the competition. It **cannot be built**: it would
    starve propanol (breaking its independently-anchored D-99 magnitude) *and* collapse sotolon's
    substrate — a large, dramatic, wrong result that would look exactly like the long-promised
    competition finally expressing itself.

    **THE SAME RATIO ARGUES BOTH WAYS, AND THAT IS THE BEAT'S RESULT.** Read as "can this pool
    supply propanol?" it is fatal to the shortcut. Read as "how big is propanol inside the node?"
    it says propanol is ~73% of the 2-ketobutyrate flux — the **dominant sink** — so an honestly
    partitioned intracellular node would couple propanol and sotolon materially. The promised
    payoff is real; it is the *location* that was wrong.

    TRIPWIRE, not a curiosity: if a future beat raises the excretion rate enough for the pool to
    supply propanol, this fails, and the design question genuinely re-opens.
    """
    base, schema = _run(aging=False)
    # Reassimilation OFF ⇒ the pool accumulates EVERY gram excreted and nothing removes it ⇒ the
    # end state IS the total throughput, exactly. No rate integrated over interpolated states
    # (D-103's finding, where a nonlinear rate over linear interpolation overstated a draw 1.3–3.5×
    # and booked more carbon out of valine than valine ever held).
    no_reassim, _ = _run(aging=False, drop=("alpha_kb_reassimilation",))

    residual = _end(base, schema, "alpha_ketobutyrate")
    throughput = _end(no_reassim, schema, "alpha_ketobutyrate")
    propanol = _end(base, schema, "propanol")

    # Anti-vacuity: all three must be real, or the ratio below is a division of zeros that
    # "agrees" for reasons having nothing to do with the finding (D-106/D-108).
    assert residual > 0.0 and throughput > residual and propanol > 0.0

    demand_mmol = propanol / M_PROPANOL * 1e3
    supply_mmol = throughput / M_ALPHA_KETOBUTYRATE * 1e3
    assert demand_mmol / supply_mmol > 2.0, (
        f"propanol demand {demand_mmol:.4f} mmol/L vs total α-KB ever excreted "
        f"{supply_mmol:.4f} mmol/L — if the pool can now supply propanol, D-109's supply "
        "argument against re-basing it no longer holds and the design must be re-measured"
    )


# -- finding 2: the de-novo dominance that makes the competition small, and is SOURCED ---------


@pytest.mark.parametrize("spec", FUSEL_SPECS, ids=lambda s: s.pool)
def test_every_fusel_is_de_novo_dominated(spec):
    """Every Ehrlich alcohol draws ≥80% of its carbon de novo, not from its precursor (D-109).

    The **sourced supply structure the fusel-side node has to preserve**, pinned here because the
    milestone will rewrite the sourcing layer that produces it: Crépin measures 2-ketobutyrate as
    81% newly synthesised, and Rollero independently puts >90% of the higher alcohols' carbon on
    the central carbon metabolism. This is the model reproducing that, and it is what makes the
    intracellular node *mostly de novo* — the property D-104 identified as the reason reality
    escapes its inverted split (isoamyl is built from KIC, so the leucine pool never faces the
    isoamyl demand).

    **It is NOT an argument that the competition is small** — D-109's first draft used it that way.
    Both arms of a keto-acid node would draw the *same de-novo 2-KB*, so de-novo dominance says
    where the carbon comes from, not whether the two consumers compete for it. The competition over
    threonine is zero for a **structural** reason instead (see
    :func:`test_alpha_kb_production_is_exactly_threonine_independent`).

    Measured from EXACT state differences via the D-104 split invariant (consumed precursor splits
    exactly ``f : (1−f)`` between the non-Ehrlich lump and the alcohol), with the other precursor
    consumers disabled so that invariant is exact. The ``n_alc/(n_alc+1)`` factor removes the D-106
    decarboxylation CO₂: the Ehrlich draw is a full mole of precursor per mole of alcohol, of which
    one carbon leaves as CO₂ rather than reaching the alcohol.
    """
    traj, schema = _run(aging=False, drop=_OTHER_PRECURSOR_CONSUMERS)
    cs = compile_scenario(_scenario(aging=False))
    pv = cs.param_values

    precursor = spec.precursor_amino_acid
    consumed = float(traj.y[schema.slice(precursor), 0][0]) - _end(traj, schema, precursor)
    made = _end(traj, schema, spec.pool)
    assert consumed > 0.0 and made > 0.0, "vacuous: nothing consumed or nothing made"

    f = pv[non_ehrlich_fraction_param(precursor)]
    n_alc = CARBON_ATOMS[spec.species]
    draw_carbon = (1.0 - f) * consumed * carbon_mass_fraction(precursor)
    alcohol_carbon_from_precursor = draw_carbon * n_alc / (n_alc + 1.0)
    total_alcohol_carbon = made * carbon_mass_fraction(spec.species)

    de_novo_share = 1.0 - alcohol_carbon_from_precursor / total_alcohol_carbon
    assert de_novo_share >= _SOURCED_DE_NOVO_FLOOR, (
        f"{spec.pool} is only {de_novo_share:.1%} de novo; D-109's finding that the keto-acid "
        "competition is negligible rests on the sourced de-novo dominance (Crépin 81%, "
        "Rollero >90% CCM)"
    )


# -- finding 3: the competition is present, correctly signed, and negligible -------------------


def test_alpha_kb_production_is_exactly_threonine_independent():
    """The competition over threonine is **structurally zero**, so "GONE" was right (D-109).

    D-109 set out to correct the work-list's "the propanol-vs-sotolon competition is GONE rather
    than wrong" and **the work-list was right**. :class:`AlphaKetobutyrateExcretion`'s rate is
    ``k · X · S/(K+S)`` — *flux-only*. Threonine's depletion gate re-routes the **carbon source**
    (threonine vs the sugar de-novo stand-in) and does not touch the **rate**. So draining
    threonine does not cost the α-KB pool one microgram, and propanol — which drains threonine —
    cannot reduce sotolon's substrate through it. Measured here at the derivative level, where it
    is exact rather than nearly-true.

    **D-107 chose this deliberately and the docstring says why**: gating the rate on threonine
    would empty the pool in a threonine-free wine and kill sotolon — the D-104 canary, one pool
    upstream. So this is not an oversight to fix; it is the reason the competition must be built
    at the **keto acid**, not at its grandparent.

    **WHY THIS TEST EXISTS, AND IT IS THE BEAT'S OWN NEAR-MISS.** D-109's first probe drove
    ``k_propanol`` from 0 to 10× (propanol ~0 → ~500 mg/L) and measured sotolon moving 0.42%,
    monotone and correctly signed. That was written up as "the competition is present, correctly
    signed, and negligible" — a *mechanism* that was never measured, attached to an assertion that
    measures something else. The channel is the **sugar ledger** (propanol's de-novo carbon comes
    out of ``S``, which every flux-linked rate reads), not threonine. The mutation test is what
    caught it: deleting the threonine draw entirely left that assertion **passing**. D-96/D-102/
    D-108's *"the sentence and the assertion are not the same claim"*, a fifth time — and this time
    it would have shipped as a correction to a claim that was **true**.
    """
    from fermentation.core.kinetics.keto_acids import AlphaKetobutyrateExcretion

    cs = compile_scenario(_scenario(aging=False))
    pv = cs.param_values
    schema = cs.process_set.schema

    def _derivative(threonine_gpl: float):
        y = schema.zeros()
        y[schema.slice("X")] = 2.0
        y[schema.slice("S")] = 100.0
        y[schema.slice("N")] = 0.1
        y[schema.slice("T")] = pv["T_ref"]
        y[schema.slice("threonine")] = threonine_gpl
        return AlphaKetobutyrateExcretion().derivatives(0.0, y, schema, pv)

    replete = _derivative(0.067)  # a real must's threonine
    starved = _derivative(0.0)  # none at all

    rate_replete = float(replete[schema.slice("alpha_ketobutyrate")][0])
    rate_starved = float(starved[schema.slice("alpha_ketobutyrate")][0])
    assert rate_replete > 0.0, "vacuous: the Process is not firing in this state at all"
    # EXACT equality, not approx: the rate does not read threonine, so this is not a small
    # dependence to bound — it is no dependence at all, and == is the only honest assertion.
    assert rate_starved == rate_replete

    # ...and the gate DID do its job — it moved the carbon source, which is the half that is real.
    # Without this the test above would also pass on a Process that had been silently disabled.
    assert float(replete[schema.slice("threonine")][0]) < 0.0  # drawn from threonine
    assert float(starved[schema.slice("threonine")][0]) == 0.0  # none to draw
    assert float(starved[schema.slice("S")][0]) < float(replete[schema.slice("S")][0])  # → sugar
