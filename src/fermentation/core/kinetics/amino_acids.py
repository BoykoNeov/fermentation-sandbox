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
from fermentation.core.kinetics.growth import (
    STORED_NITROGEN_KEY,
    add_assimilable_nitrogen,
    biomass_growth_rate,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Re-exported for the consumers that named this module as the home of the representative
#: species before D-100 moved the pool registry to
#: :mod:`~fermentation.core.kinetics.amino_acid_pools`. The species itself is unchanged.
__all__ = ["AMINO_ACID_SPECIES", "AminoAcidAssimilation", "AssimilableNitrogenUptake"]


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
    touches = (*(spec.pool for spec in ASSIMILABLE_SPECS), "N", STORED_NITROGEN_KEY, "S")
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
        # Refund the displaced biomass nitrogen on the SAME split growth drew it on (D-250).
        # This is a correctness coupling, not tidiness: growth's draw is now shared between the
        # must pool and the intracellular store, so a refund booked wholly to `N` would credit a
        # pool growth only partly debited and `N` would drift UP -- the D-248 charge defect
        # returning through a second door.
        add_assimilable_nitrogen(d, y, schema, nitrogen)

        # Refund the displaced biomass carbon to sugar — the inverse of growth's draw,
        # distributed across sugar slots by their current carbon content so that
        # Σ_i (d[S_i]·c_i) equals it exactly. base_dx > 0 guarantees s_total > 0, so the
        # refund always has somewhere to go (no silent carbon leak). Shares the single
        # carbon-routing helper with the fusel re-route (D-33) so the draw and its inverse
        # can never drift apart (single source of truth, decision D-8).
        refund_carbon_to_sugar(d, y, schema, carbon)  # [g C/L/h]
        return d


class AssimilableNitrogenUptake(Process):
    """Amino-acid nitrogen taken up *beyond* growth's anabolic demand (decision D-248).

    **The gap this closes (D-246 §5 → D-247 §7).** Crépin *et al.* 2017's Data Set S1 measures
    the residual as well as the initial amino acids of its synthetic must, and its yeast consumes
    essentially all of the assimilable nitrogen — the end-of-fermentation block reads 0.0000 for
    ammonium, arginine, glutamine and every precursor, with the YAN column falling 9.636 → ~0.02
    mM (**0.2 % left**), exhausted at N_T ≈ 28 h. Reconstructed on that same must this model left
    **40.8 %** standing, and held at 30 and 60 days instead of 14 the residual did not move in the
    fourth decimal. It was never slowness.

    **The arithmetic of the defect, which is what makes the repair structural rather than a knob.**
    Before this Process the *only* route from the speciated pools into biomass nitrogen was
    :class:`AminoAcidAssimilation`, whose rate is ``ψ·gate·f_N·base_dx`` — strictly **below**
    growth's own draw ``f_N·base_dx`` for every ``ψ·gate ≤ 1``. So the net ammonium derivative
    ``f_N·base_dx·(ψ·gate − 1)`` is negative at every state that has ever existed: ``N`` can only
    fall. When it reaches zero growth's Monod term ``N/(K_n + N)`` shuts growth off, ``base_dx``
    goes to zero, and the swap — being proportional to ``base_dx`` — stops with it. The amino-acid
    pools then freeze wherever they stood. **Uptake could only ever consume what growth demanded,
    and growth's nitrogen signal is the ammonium slot alone**, which the speciated pools are not
    in. Real yeast over-accumulate assimilable nitrogen far past immediate anabolic need, storing
    the surplus (vacuolar amino-acid pools) — the flux this model had no way to express.

    **Two published anchors, and the load-bearing one is INTERNAL.** The compile seam already
    overrides ``biomass_N_fraction`` to ``1/Y_X/N(N_init)`` from Coleman, Fish & Block (2007), for
    the express purpose that biomass comes out at ``Y_X/N × N_init``
    (:func:`~fermentation.scenario.compile._apply_nitrogen_dependent_yield`). **That identity holds
    only under complete consumption.** Measured on Crépin's must: with all the nitrogen in the
    ammonium slot the run reached 3.445 g/L against the fit's own 3.472 (0.8 %), and with the same
    nitrogen speciated into pools it reached 2.12 — **61 % of a regression the model itself
    compiles in**. That anchor is internal, already shipped and fitted to nothing here. Crépin's
    0.2 % residual and 28 h exhaustion are the *independent* check, deliberately not the target.

    **The rate, and why it is un-coupled.** ::

        ρ_N = r · mu_max · f_N · X · gate(aa)              [g N / L / h]
        gate(aa) = Σaa / (K_amino_acids·Σf + Σaa)          (the shared D-100 depletion gate)

    ``mu_max · f_N`` is the *maximum specific* nitrogen demand of growth — nitrogen per gram of
    cell per hour at the top of the growth curve — so ``r``
    (``amino_acid_uptake_capacity_ratio``) is a dimensionless statement about transport
    **capacity** relative to peak anabolic demand, and no new dimensional constant enters. The
    un-coupling is that the rate reads ``mu_max``, a constant, and **not**
    :func:`~fermentation.core.kinetics.growth.biomass_growth_rate`: transport scales with the
    cells present, not with what they happen to need this instant, so it keeps running after
    growth has stopped. That is the entire mechanism, stated as a rate law.

    **The carbon does NOT go back to sugar, and that is the design's crux.** The D-32 swap may
    refund its drawn carbon to ``S`` only because its rate is proportional to growth's own draw,
    which is what bounds the refund at ``ψ·gate·(aa C:N)/(biomass C:N) ≈ 0.30·ψ·gate`` of the
    carbon growth removed. A flux that is *not* proportional to ``base_dx`` has no such bound —
    at ``base_dx = 0`` it would refund carbon against a draw of zero and **create hexose**, the
    gluconeogenesis the whole D-32 design forbids and which no clamp fixes without a C⁰ kink for
    the BDF solver. So the nitrogen is refunded to ``N`` and the skeleton is parked in
    ``amino_acid_skeleton_carbon``, a carbon-only pool weighted 1.0 on ``total_carbon``.

    That parking is a documented **stand-in**, in the D-19/D-31 tradition, and it is the
    conservative direction: real yeast do route some of an assimilated amino acid's skeleton into
    biomass and central metabolism, so booking none of it back to ``S`` slightly *over*-charges
    sugar for biomass carbon. The magnitude is ~0.15 g C/L against ~96 g C/L of must carbon
    (~0.16 %), i.e. immaterial to every carbon-ratio benchmark, and the error's sign cannot create
    sugar or ethanol. What it must not be read as is a claim about where arginine's carbon goes.

    **Scope: the identity-agnostic pair ONLY** (``amino_acids``/arginine and
    ``amino_acids_generic``/glutamine), which is a measurement discipline and not a shortcut.
    Those two are where the residual sits — measured at 61.4 % of their initial mass each while
    every precursor pool was already at zero — so nothing is lost by the narrow scope. What the
    scope *buys* is attributability: the six precursor pools are untouched, so this Process cannot
    move a fusel's de-novo share through precursor starvation, and the only channel left to those
    shares is the biomass denominator. A wider draw would move both at once and no result would be
    nameable. It also keeps D-100's cross-subsystem starvation shut: the Ehrlich re-route still
    touches no member of this pair.

    **Isolability (prime directive #3).** The gate is exactly 0 on an empty pool and the compile
    seam disables this Process along with the rest of the amino-acid ledger when
    ``amino_acids_gpl <= 0``, so an undosed wine run is byte-for-byte the validated core and the
    Coleman reconstruction is untouched. It is **not** a target of the growth Arrhenius or
    carrying-capacity modifiers: the swap is scaled by them so its refunds track growth's realised
    draw, and this Process has no draw to track — importing that scaling would re-introduce the
    coupling the Process exists to remove. The consequence is that uptake carries **no temperature
    dependence** in v1, a named simplification: a sourced activation energy for amino-acid
    permease transport is not in hand, and borrowing ``E_a_growth`` would be an unsourced claim
    dressed as fidelity (the D-98 trap).

    Tier: **speculative** — the form is standard (biomass-proportional, saturable transport) but
    ``r`` is an author estimate, the storage destination is a stand-in, and the surplus nitrogen
    is booked as extracellular ammonium rather than an intracellular reserve.
    """

    name = "assimilable_nitrogen_uptake"
    tier = Tier.SPECULATIVE
    #: Debits the two identity-agnostic pools, refunds their nitrogen to ``N`` and parks their
    #: carbon in ``amino_acid_skeleton_carbon``. Touches neither ``X`` (it funds no biomass of its
    #: own — growth does that, out of the ammonium this fills) nor ``S`` (see the docstring: a
    #: sugar refund here would create hexose) nor any of the six precursor pools.
    touches = (
        *(spec.pool for spec in ASSIMILABLE_SPECS),
        STORED_NITROGEN_KEY,
        "amino_acid_skeleton_carbon",
    )
    #: ``mu_max`` and ``biomass_N_fraction`` set the *capacity scale* (peak specific nitrogen
    #: demand) — read as constants, never as the instantaneous growth rate, which is the whole
    #: point. ``K_amino_acids`` and the pair's two ``must_aa_fraction_*`` shares are the shared
    #: D-100 depletion gate. ``ProcessSet.tier_of`` folds these into the output tier of ``N``,
    #: the two pools and the skeleton pool when this Process is enabled (D-1).
    reads: tuple[str, ...] = (
        "mu_max",
        "biomass_N_fraction",
        "amino_acid_uptake_capacity_ratio",
        "K_amino_acids",
        *(spec.fraction_param for spec in ASSIMILABLE_SPECS),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if "amino_acid_skeleton_carbon" not in schema:
            return d  # beer: the amino-acid ledger is wine-only (D-30/D-32)
        x = float(y[schema.slice("X")][0])
        if x <= 0.0:
            return d  # no cells ⇒ no transport capacity
        gate = depletion_gate(y, schema, params, ASSIMILABLE_SPECS)
        if gate <= 0.0:
            return d  # empty pool ⇒ nothing to take up (also the undosed no-op)

        # The capacity, NOT the demand: mu_max·f_N is the maximum specific nitrogen draw of
        # growth, so this is transport sized against peak anabolic need and scaled by the cells
        # present. `biomass_growth_rate` is deliberately not called — reading it would restore
        # exactly the coupling this Process removes.
        nitrogen = (
            params["amino_acid_uptake_capacity_ratio"]
            * params["mu_max"]
            * params["biomass_N_fraction"]
            * x
            * gate
        )
        # Split across {arginine, generic} in proportion to the nitrogen each holds (the shared
        # D-100 idiom, so this and the swap can never drift apart on the arithmetic), and take
        # the carbon that mass carried.
        carbon = draw_assimilable_nitrogen(d, y, schema, nitrogen)
        # The surplus goes INTO THE CELL (D-250), not into the must's ammonium. D-248 booked it
        # in `N`, which the acid-base balance reads at the must's mean nitrogen charge -- so
        # nitrogen the yeast had already transported in went on titrating the must, worth up to
        # +0.215 pH mid-run on a 2 g/L dosed must. `stored_nitrogen` is in no charge balance.
        d[schema.slice(STORED_NITROGEN_KEY)] = nitrogen
        # Park the skeleton. NOT refunded to sugar — see the docstring: unbounded by growth's
        # draw, a sugar refund would create hexose whenever growth is stopped.
        d[schema.slice("amino_acid_skeleton_carbon")] = carbon
        return d
