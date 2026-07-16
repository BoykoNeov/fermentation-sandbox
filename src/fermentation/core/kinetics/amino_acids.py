"""Amino-acid assimilation — the toggleable amino-acid ledger (decision D-32).

**What this closes (decision D-23 → D-32).** Yeast build biomass mostly from
amino acids, not hexose, but the validated core sources *all* biomass carbon from
sugar and *all* biomass nitrogen from the lumped ammonium ``N`` pool, and ``N`` is
deliberately carbon-free in ``total_carbon`` (decision D-19). Making amino acids a
carbon source is therefore a change to the protected carbon *and* nitrogen ledgers.
The owner's toggleable **amino-acid ledger** (D-23) restores isolability: a
``default=0`` ``amino_acids`` pool that, when dosed, funds a fraction of biomass
from amino acids instead of sugar+ammonium. Following the advisor's refinement, it
is a **separate isolable Process** — a pure *swap* — rather than a branch inside the
core's hottest kinetic, so growth (and the Coleman reconstruction) stay byte-for-byte
and isolability is structural, not a tested coincidence.

**The swap.** For biomass built at the growth rate ``dX/dt`` (the shared
:func:`~fermentation.core.kinetics.growth.biomass_growth_rate`), this Process
consumes amino acids at rate ``ρ`` and:

  * **debits** the amino-acid pool by ``ρ`` (``d[amino_acids] = -ρ``),
  * **refunds ammonium** ``N`` by the nitrogen that ``ρ`` carries
    (``d[N] = +ρ·y_N``), and
  * **refunds sugar** by the carbon that ``ρ`` carries (``d[S] += +ρ·y_C``),

where ``y_N``/``y_C`` are arginine's nitrogen/carbon mass fractions (the
representative amino acid, D-32). Biomass ``X`` is untouched — growth still builds
it — so the swap is a pure transfer ``aa → S`` (carbon) and ``aa → N`` (nitrogen):
**carbon- and nitrogen-neutral by construction** for any ``ρ`` (``total_carbon`` and
``total_nitrogen`` close, the aa pool now weighted in both). The physical reading:
using amino acids for biomass *spares* sugar for ethanol, so the sugar and ammonium
that growth's stoichiometry charged are credited back.

BOOKKEEPING CAVEAT (the D-19/D-31 stand-in discipline): mechanically the aa carbon is
refunded to **sugar**, biomass carbon still comes from growth's sugar draw, and the
spared sugar then ferments to ethanol — arginine's carbon skeleton is booked as spared
hexose, not tracked through arginine catabolism. This is carbon-closing and physically
defensible (aa-fed biomass really does spare sugar for ethanol), but it is a stand-in,
not a claim about arginine's metabolic fate. One consequence: dosing amino acids nudges
ethanol *up* by ~0.15–0.3 % of sugar (the spared carbon), tiny and — since the §2.2
benchmarks run undosed — leaving them untouched.

**The rate is nitrogen-anchored (decision D-32).** Amino acids *are* part of
yeast-assimilable nitrogen, so ``ρ`` is tied to the fraction of biomass **nitrogen**
sourced from the pool::

    ρ = ψ · gate(aa) · f_N · base_dx / y_N          [g aa / L / h]
    gate(aa) = aa / (K_amino_acids + aa)            (smooth availability, → 0 as aa → 0)

with ``ψ = amino_acid_assimilation_fraction ∈ [0, 1]`` the max aa-funded share of
biomass nitrogen. The nitrogen refund is then ``ρ·y_N = ψ·gate·f_N·base_dx ≤
f_N·base_dx`` (growth's nitrogen draw) for all ``ψ·gate ≤ 1`` — **this Process** never
over-refunds, so it needs no deamination branch of its own.

.. warning::

   **That bound is this swap's alone, and it is NOT the system's (decision D-104).** It was
   written when the swap was the only Process refunding biomass nitrogen. Since D-104
   :class:`~fermentation.core.kinetics.precursor_fates.PrecursorNonEhrlichFates` refunds the
   *precursors'* nitrogen too, and nothing bounds the pair against ``f_N·base_dx``: measured at
   the shipped ``ψ = 0.5`` with a 1 g/L dose, the **joint** refund reaches **1.171× growth's draw
   at pitch** (D-104 measured 1.04×; **D-106 raised it to 1.171×** by charging the Ehrlich
   decarboxylation CO₂, which makes the re-route consume a full mole of precursor per alcohol
   instead of ``(n-1)/n`` — and a full mole carries a full mole of nitrogen to deaminate). The
   excess is *net deamination* — physical, and needing no branch, because the
   refund is always the drawn molecule's own nitrogen and the sign of the net falls out of the
   arithmetic. Nitrogen still closes exactly (it is transferred from the precursor pools, never
   created). The guarantee that *does* still bind the pair is the **carbon** one below — no
   sugar creation — and it is pinned jointly in ``tests/test_precursor_fates.py``.

The carbon refund is
``ρ·y_C = ψ·gate·base_dx·f_N·(y_C/y_N)``; dividing by growth's carbon draw
``f_C·base_dx`` gives ``ψ·gate·(f_N/f_C)·(y_C/y_N) = ψ·gate·(aa C:N)/(biomass C:N)``.
With arginine (mass C:N ≈ 1.29) and biomass (``f_C/f_N`` ≈ 4.3) this is
``≈ 0.30·ψ·gate ≤ 0.30`` — the carbon refund is *strictly* below growth's demand for
any ``ψ ≤ 1``, so the swap **never creates hexose** (gluconeogenesis to sugar, which
fermenting yeast do not do) and needs no clamp. That N-rich representative is the
load-bearing modelling choice; a carbon-rich amino acid would force a clamp (a C⁰
kink the stiff BDF solver catches on) or leak sugar.

**Why the modifier scaling matters — the correctness crux (decision D-32).** The
guarantee above uses ``base_dx``, growth's *pre-modifier* rate. Growth's realised
biomass is ``base_dx · M`` where ``M`` is the product of the Arrhenius and (opt-in)
carrying-capacity :class:`~fermentation.core.process.RateModifier` factors applied by
:class:`~fermentation.core.process.ProcessSet`. If the swap refunded at ``base_dx``
while growth drew at ``M·base_dx``, then at ``M < 0.30`` (cold ferment, or the
carrying cap near saturation with nitrogen still available — the D-30 residual-N
regime) the refund would exceed the draw and **create sugar**. The fix: the growth
Arrhenius and carrying-capacity modifiers scale *this Process too* (they name it in
their ``modifies``), so refund and draw carry the *same* ``M``::

    net dS = M·f_C·base_dx·(0.30·ψ·gate − 1) ≤ 0
    net dN = M·f_N·base_dx·(ψ·gate − 1)      ≤ 0

for all ``ψ·gate ≤ 1`` — never creates sugar, never deaminates. This is verified at
the *ProcessSet* level (not the raw derivatives) at ``M < 1`` states in the tests,
because at the reference temperature ``M = 1`` and the mismatch would never fire.

**Isolability (undosed-only).** When ``amino_acids`` is empty the compile seam
*disables* this Process (so its speculative tier does not drag growth's ``S``/``N``
outputs down and no work is done), and even enabled the availability gate → 0 at
``aa = 0`` — an undosed wine run is byte-for-byte the validated core. **Dosed**, the
swap *correctly* perturbs the trajectory: refunded ``N``/``S`` raise the pools growth
reads on the next step, so dosing amino acids behaves like supplementary YAN
(nitrogen lasts longer ⇒ more biomass / more sugar consumed) — a second-order
feedback, not a first-order growth edit (growth's derivatives are untouched).

Tier: **speculative** — the swap form is sound but ``ψ`` and ``K_amino_acids`` are
author estimates and the single-representative-amino-acid lumping is a simplification.

SCOPE (v1): the swap only (primary-fermentation yeast carbon/nitrogen honesty). The
D-19 fusel Ehrlich re-route (drawing fusel carbon from this pool instead of sugar) is
its natural later home but needs the deamination branch, so it is deferred (D-23).
Wine-only, mirroring the wine-only nitrogen-model wiring (D-30); beer is deferred.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.kinetics.amino_acid_pools import (
    AMINO_ACID_SPECIES,
    ASSIMILABLE_SPECS,
    depletion_gate,
    draw_assimilable_nitrogen,
)
from fermentation.core.kinetics.carbon_routing import refund_carbon_to_sugar
from fermentation.core.kinetics.growth import biomass_growth_rate
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Re-exported for the consumers that named this module as the home of the representative
#: species before D-100 moved the pool registry to
#: :mod:`~fermentation.core.kinetics.amino_acid_pools`. The species itself is unchanged.
__all__ = ["AMINO_ACID_SPECIES", "AminoAcidAssimilation"]


class AminoAcidAssimilation(Process):
    """Nitrogen-anchored amino-acid → biomass swap (decision D-32).

    Consumes amino acids at ``ρ = ψ·gate(aa)·f_N·base_dx/y_N`` and refunds the carbon
    to sugar and the nitrogen to ammonium ``N``, leaving biomass untouched — a
    carbon- and nitrogen-neutral transfer ``aa → S`` / ``aa → N`` (module docstring).
    The growth Arrhenius/carrying-capacity modifiers scale this Process too, so its
    refunds track growth's *realised* draw and never create sugar (decision D-32).
    """

    name = "amino_acid_assimilation"
    tier = Tier.SPECULATIVE
    #: Refunds carbon to ``S`` and nitrogen to ``N``; debits the two identity-agnostic pools
    #: ``amino_acids`` (arginine) + ``amino_acids_generic`` (D-100). Does NOT touch ``X``
    #: (growth builds biomass; this only re-sources its atoms), and does not touch the six
    #: precursor pools — yeast build biomass from any assimilable amino acid, but leucine's
    #: fate in this model is the Ehrlich pathway (D-33/D-99), not protein.
    touches = (*(spec.pool for spec in ASSIMILABLE_SPECS), "N", "S")
    #: ``mu_max``/``K_s``/``K_n``/``biomass_N_fraction`` reach it through the shared
    #: growth-rate helper (the swap anchors to the same base rate); ``ψ``, ``K_amino_acids``
    #: and the two assimilable ``must_aa_fraction_*`` shares (which scale the D-100
    #: relative-depletion gate) are its own. ``ProcessSet.tier_of`` folds these into the
    #: swap's output tier for ``S``/``N``/the assimilable pools when it is enabled (D-1).
    reads: tuple[str, ...] = (
        "mu_max",
        "K_s",
        "K_n",
        "biomass_N_fraction",
        "amino_acid_assimilation_fraction",
        "K_amino_acids",
        *(spec.fraction_param for spec in ASSIMILABLE_SPECS),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        base_dx = biomass_growth_rate(y, schema, params)
        if base_dx <= 0.0:
            return d  # no growth ⇒ no biomass to re-source ⇒ nothing to swap
        # The identity-agnostic substrate's relative-depletion gate (D-100): {arginine, generic}
        # with K scaled by their combined must-spectrum share. → 0 on an empty pool, so this is
        # still the undosed no-op, and the draw can never drive either pool negative.
        gate = depletion_gate(y, schema, params, ASSIMILABLE_SPECS)
        if gate <= 0.0:
            return d  # empty pool ⇒ nothing to assimilate (also the undosed no-op)

        # The nitrogen demand [g N/L/h]: the fraction ψ·gate of the biomass nitrogen growth is
        # drawing. ≤ f_N·base_dx for all ψ·gate ≤ 1 (module docstring), so the N refund never
        # exceeds growth's draw and no deamination branch is needed.
        nitrogen = (
            params["amino_acid_assimilation_fraction"]
            * gate
            * params["biomass_N_fraction"]
            * base_dx
        )
        # Split that demand across {arginine, generic} in proportion to the nitrogen each holds,
        # debiting both (D-100). Returns the carbon those amino acids carry — a blend of
        # arginine's C:N ≈ 1.29 and glutamine's ≈ 2.14, both far below biomass's ≈ 4.3, so the
        # carbon refund stays strictly below growth's demand and the swap still cannot create
        # hexose for any ψ ≤ 1 (the guarantee is structural for the blend, not just for arginine).
        carbon = draw_assimilable_nitrogen(d, y, schema, nitrogen)
        d[schema.slice("N")] = nitrogen  # refund displaced biomass nitrogen to ammonium

        # Refund the displaced biomass carbon to sugar — the inverse of growth's draw,
        # distributed across sugar slots by their current carbon content so that
        # Σ_i (d[S_i]·c_i) equals it exactly. base_dx > 0 guarantees s_total > 0, so the
        # refund always has somewhere to go (no silent carbon leak). Shares the single
        # carbon-routing helper with the fusel re-route (D-33) so the draw and its inverse
        # can never drift apart (single source of truth, decision D-8).
        refund_carbon_to_sugar(d, y, schema, carbon)  # [g C/L/h]
        return d
