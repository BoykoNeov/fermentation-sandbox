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
f_N·base_dx`` (growth's nitrogen draw) for all ``ψ·gate ≤ 1`` — never over-refunds,
so no deamination branch is needed in v1 (excess-aa deamination is deferred with the
fusel Ehrlich re-route, D-23/D-32). The carbon refund is
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

from fermentation.core.chemistry import (
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics.carbon_routing import refund_carbon_to_sugar
from fermentation.core.kinetics.growth import biomass_growth_rate
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The representative species for the lumped assimilable amino-acid pool (decision D-32).
#: Arginine is the dominant yeast-assimilable amino acid in grape must and is N-rich (mass
#: C:N ≈ 1.29 ≪ biomass ≈ 4.3), the property that keeps the carbon refund below growth's
#: demand for any ψ ≤ 1 (see the module docstring). Kept as a module constant so the swap
#: and the conservation weighting name one species. Public so the D-33 fusel Ehrlich
#: re-route (which sources fusel carbon from this same pool) references the identical species.
AMINO_ACID_SPECIES = "arginine"


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
    #: Refunds carbon to ``S`` and nitrogen to ``N``; debits the ``amino_acids`` pool.
    #: Does NOT touch ``X`` (growth builds biomass; this only re-sources its atoms).
    touches = ("amino_acids", "N", "S")
    #: ``mu_max``/``K_s``/``K_n``/``biomass_N_fraction`` reach it through the shared
    #: growth-rate helper (the swap anchors to the same base rate); ``ψ`` and
    #: ``K_amino_acids`` are its own. ``ProcessSet.tier_of`` folds these into the
    #: swap's output tier for ``S``/``N``/``amino_acids`` when it is enabled (D-1).
    reads: tuple[str, ...] = (
        "mu_max",
        "K_s",
        "K_n",
        "biomass_N_fraction",
        "amino_acid_assimilation_fraction",
        "K_amino_acids",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        base_dx = biomass_growth_rate(y, schema, params)
        if base_dx <= 0.0:
            return d  # no growth ⇒ no biomass to re-source ⇒ nothing to swap
        aa = max(float(y[schema.slice("amino_acids")][0]), 0.0)
        if aa <= 0.0:
            return d  # empty pool ⇒ nothing to assimilate (also the undosed no-op)

        gate = aa / (params["K_amino_acids"] + aa)  # smooth availability, in [0, 1)
        y_n = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
        y_c = carbon_mass_fraction(AMINO_ACID_SPECIES)
        # ρ [g aa/L/h]: aa consumption anchored to the fraction ψ·gate of biomass
        # nitrogen sourced from the pool. ρ·y_N ≤ f_N·base_dx and ρ·y_C < f_C·base_dx
        # for all ψ·gate ≤ 1 (module docstring), so neither refund exceeds growth's draw.
        rho = (
            params["amino_acid_assimilation_fraction"]
            * gate
            * params["biomass_N_fraction"]
            * base_dx
            / y_n
        )
        d[schema.slice("amino_acids")] = -rho
        d[schema.slice("N")] = rho * y_n  # refund displaced biomass nitrogen to ammonium

        # Refund the displaced biomass carbon to sugar — the inverse of growth's draw,
        # distributed across sugar slots by their current carbon content so that
        # Σ_i (d[S_i]·c_i) = ρ·y_C exactly. base_dx > 0 guarantees s_total > 0, so the
        # refund always has somewhere to go (no silent carbon leak). Shares the single
        # carbon-routing helper with the fusel re-route (D-33) so the draw and its inverse
        # can never drift apart (single source of truth, decision D-8).
        refund_carbon_to_sugar(d, y, schema, rho * y_c)  # ρ·y_C [g C/L/h]
        return d
