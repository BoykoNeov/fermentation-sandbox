"""Excreted keto-acid overflow pools — the non-acetaldehyde SO₂-binding carbonyls (D-49).

Acetaldehyde is the *principal* SO₂-binder in wine, but it is not the only one: yeast
excretes **overflow keto-acids** — pyruvate first, α-ketoglutarate next — during active
fermentation, and those persist in the finished wine (10s–100s mg/L) as real, measured
SO₂-binding carbonyls (Jackowetz & Mira de Orduña 2013; Burroughs & Sparks 1973). This
module builds the first, pyruvate (α-ketoglutarate follows in D-50), as an **excreted side
pool**, so the free/bound-SO₂ split (D-28) and the acetaldehyde protection (D-47) can share
SO₂ across *all* the carbonyls that really compete for it (the D-51 multi-carbonyl binding).

**Why a side pool, not an on-pathway intermediate (the load-bearing modelling choice, D-49).**
The tempting "faithful" model — route acetaldehyde's carbon *through* pyruvate (its real
metabolic precursor: pyruvate → acetaldehyde + CO₂ → ethanol) — was designed and **rejected**.
It conflates two physically distinct pools: acetaldehyde's precursor is the *intracellular*
flux intermediate (enormous flux, vanishing pool, never persists, never measured), whereas
the SO₂-binding pyruvate is the *extracellular* excreted overflow residual (small flux,
persistent, measured). One pool cannot be both, and the persistence mechanism the rework
needed — "SO₂ protects pyruvate from pyruvate decarboxylase" — is not real: pyruvate
decarboxylase is intracellular, the excreted residual never meets it. Worse, that shielding
would make dosed SO₂ *sequester acetaldehyde's precursor* and **suppress** acetaldehyde — the
exact opposite of the SO₂-induced over-production the model just shipped (D-48). So the
excreted-overflow side pool is the *more* faithful structure for the quantity that matters
here, and acetaldehyde / D-27 / D-47 / D-48 stay untouched.

**Carbon — drawn from sugar, returned to ethanol + CO₂ (closing on the existing ledger).**
Following the D-19/D-26 byproduct idiom:

* :class:`PyruvateExcretion` draws pyruvate's carbon *out of ``S``* (via
  :func:`~fermentation.core.kinetics.carbon_routing.draw_carbon_from_sugar`, booked at the C3
  pyruvate fraction) while the yeast ferments — flux-linked, so it stops at dryness, leaving a
  full overflow pool.
* :class:`PyruvateReassimilation` returns it: pyruvate re-enters metabolism and is oxidised to
  ethanol + CO₂ (``C3 → C2 + C1``, one mole each — carbon-closing on the ledger exactly like
  malic → lactic + CO₂, D-23). **Flux-linked (co-metabolic)**, sharing excretion's
  ``X·S/(K_sugar_uptake+S)`` shape, so it clears pyruvate while the yeast ferments and **stops at
  dryness** — which *freezes* the pool at its dryness value as a **persistent finished-wine
  residual** for SO₂ to bind. This is deliberately the *opposite* of the acetaldehyde-reduction
  template (D-27, no-flux, viable-``X``-gated, runs through the post-ferment rest): ADH keeps
  reducing acetaldehyde after sugar is gone, but a normal ferment finishes with the yeast still
  viable, so a no-flux gate would drain overflow pyruvate to ~0 over the long tail — a
  crash-dependent residual a clean ferment never strands. Co-metabolic flux-linking pegs the
  residual to *end-of-fermentation*, not yeast death, so it is crash- and duration-independent
  (decision D-49, option A). Held temperature-flat (a documented v1 simplification, like the
  acetaldehyde production and α-acetolactate excretion).

Because both terms ride the same flux shape, the pool rises monotonically to the quasi-steady
plateau ``k_pyruvate_excretion / k_pyruvate_reassimilation`` and freezes there — a v1
simplification that drops the real mid-ferment *peak-then-decline* transient (overflow pyruvate's
exponential-growth excretion). Nothing reads the peak; D-51 reads only the residual, so the growth-
coupled excretion that would restore the transient (option B) is deferred. Re-assimilation returns
carbon to ``E``/``CO2`` (not ``S``) because that is pyruvate's real metabolic fate (forward to
ethanol, not back to sugar) — and a refund-to-sugar would be a no-op at post-dryness ``S = 0``
anyway. ``total_carbon`` (which weights ``pyruvate`` at its own C3 fraction) closes to **machine
precision** through the whole excrete-then-reassimilate course; ``total_mass`` carries the usual
small gap (the oxidation moves untracked NAD(P)H / O) — carbon is the invariant, as for the
diacetyl rest (D-26) and the acetaldehyde buffer (D-27).

**Isolability (prime directive #3).** The pool lives in its own ``_KETO_ACID_PROCESSES`` tuple
(:mod:`fermentation.core.media`) and is **wine-only** (v1): the SO₂-binding competition it
exists for is a wine readout, and no §2.2 beer benchmark asserts a keto-acid level — beer
overflow pyruvate is deferred. A ProcessSet built without it is the prior core. Unlike the
byte-for-byte-isolable acetaldehyde buffer, excretion touches ``S`` and re-assimilation touches
``E``/``CO2``, so the pool routes a genuine (small) slice of sugar carbon on a detour to
ethanol; the endpoint difference from the pool-off core is only the **stranded residual** (a
few tens of mg/L of sugar carbon parked as pyruvate rather than fermented on to ethanol), which
is ≪ 0.1 % of ABV — the §2.2 CO₂/ABV/realised-yield benchmarks are preserved far below any
tolerance (verified at endpoint, not just transiently).

**Tiers.** Both rate constants are order-of-magnitude estimates, so both Processes are
**speculative**; parameter-tier propagation (D-1) caps the ``pyruvate`` output tier at
speculative regardless. The excreted-overflow *mechanism* (yeast excretes pyruvate during
active fermentation, re-assimilating it late) is textbook; only the RATE magnitudes are the
author's estimates, tuned so the pool peaks in the real ~100s mg/L range and settles to a
~10–40 mg/L finished-wine residual.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import (
    M_CO2,
    M_ETHANOL,
    M_PYRUVATE,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.carbon_routing import (
    draw_carbon_from_sugar,
    fermentative_flux_shape,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The species whose C3 formula carbon-accounts the ``pyruvate`` pool (decision D-49). Its
#: carbon mass fraction weights both the sugar draw here and the pool in ``total_carbon`` —
#: one chemistry source of truth, so the draw and the conservation check cannot disagree.
_PYRUVATE_SPECIES = "pyruvate"


class PyruvateExcretion(Process):
    """Overflow-pyruvate excretion — fills the excreted keto-acid pool during fermentation.

    ``d(pyruvate)/dt = k_pyruvate_excretion · X · S_total/(K_sugar_uptake + S_total)`` with the
    carbon drawn *out of ``S``* (booked as pyruvate, C3). Yeast overflows pyruvate while it
    ferments, so production is tied to the biomass-catalysed sugar flux (sharing
    ``K_sugar_uptake`` with the uptake Process) and **stops at dryness**. The equally flux-linked
    re-assimilation draws the pool toward the quasi-steady plateau
    ``k_pyruvate_excretion / k_pyruvate_reassimilation`` during the ferment and then freezes it at
    dryness — the persistent finished-wine residual (see :class:`PyruvateReassimilation`).

    Held **temperature-flat** (no explicit Arrhenius factor) as a documented v1 simplification,
    like the α-acetolactate excretion (D-26) and the acetaldehyde production (D-27): no §2.2
    benchmark asserts a pyruvate temperature direction in v1. Touches ``pyruvate`` and ``S``
    only, so with the keto-acid pool off the core is byte-for-byte and with it on only ``dS``
    gains a negative excretion term. Tier **speculative** (rate magnitude estimate).
    """

    name = "pyruvate_excretion"
    tier = Tier.SPECULATIVE
    touches = ("pyruvate", "S")
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks;
    #: ``k_pyruvate_excretion`` sets the excretion magnitude. Their tiers cap ``pyruvate``'s
    #: output tier via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_pyruvate_excretion", "K_sugar_uptake")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        rate = params["k_pyruvate_excretion"] * flux
        d[schema.slice("pyruvate")] = rate
        draw_carbon_from_sugar(d, y, schema, rate * carbon_mass_fraction(_PYRUVATE_SPECIES))
        return d


class PyruvateReassimilation(Process):
    """Co-metabolic re-assimilation of overflow pyruvate → ethanol + CO₂ — freezes the residual.

    ``d(pyruvate)/dt = −L`` with ``L = k_pyruvate_reassimilation · X · S/(K_sugar_uptake+S) ·
    [pyruvate]`` (mass loss), returned to metabolism as ``d(E)/dt = +r·M_ethanol`` and
    ``d(CO2)/dt = +r·M_CO2`` with the molar turnover ``r = L/M_pyruvate``. Pyruvate is oxidised to
    ethanol + CO₂ (``C3 → C2 + C1``, one mole each), carbon-closing mole-for-mole on the existing
    ledger like malic → lactic + CO₂ (D-23); the carbon returns to ``E``/``CO2`` rather than ``S``
    because post-dryness ``S`` is zero and a refund there would be a no-op that destroys carbon.

    **FLUX-LINKED (co-metabolic), sharing excretion's ``X · S/(K_sugar_uptake+S)`` shape — NOT the
    no-flux ADH idiom.** Overflow-pyruvate re-assimilation tracks *active fermentation*, so it
    clears pyruvate while the yeast ferments and **stops at dryness** (``S → 0`` ⇒ flux ``→ 0`` ⇒
    both this term and excretion die together). That is the load-bearing choice, and it is the
    *opposite* of the acetaldehyde reduction template (D-27, no-flux, runs through the post-ferment
    rest): ADH genuinely keeps reducing acetaldehyde after sugar is gone, but overflow pyruvate is
    co-metabolic — a normal ferment finishes with the yeast still viable, so a viable-``X`` no-flux
    gate would drain the pool to ~0 over the long tail (a crash-dependent residual that a clean
    ferment never strands). Flux-linking instead **freezes the pool at its dryness value**: with
    excretion also flux-linked the pool rises toward the quasi-steady plateau
    ``k_pyruvate_excretion / k_pyruvate_reassimilation`` and holds there at dryness — a
    **persistent finished-wine residual** that is crash-*independent* and run-duration-*independent*
    (decision D-49, owner's finish-and-wrap-up call: option A). That residual — set by the
    excretion/re-assimilation *ratio* here, NOT by SO₂ — is what binds SO₂ in the D-51
    multi-carbonyl equilibrium.

    V1 SIMPLIFICATION (documented, in-idiom with the temperature-flat calls): because both terms
    ride the same flux shape, the pool rises monotonically to the plateau rather than showing the
    real mid-ferment *peak-then-decline* (the descriptive transient overflow pyruvate really shows
    during exponential growth). Nothing reads the peak — D-51 reads only the finished-wine residual
    — so the transient is dropped in v1; the growth-coupled excretion that would restore it (option
    B) is deferred. Held temperature-flat (v1). ``pyruvate`` is clamped ≥ 0 and the shared
    ``fermentative_flux_shape`` clamps ``X``/``S`` against solver undershoot. Mass carries a small
    gap (the oxidation moves untracked NAD(P)H) — carbon is the invariant. Tier **speculative**.
    """

    name = "pyruvate_reassimilation"
    tier = Tier.SPECULATIVE
    touches = ("pyruvate", "E", "CO2")
    #: ``k_pyruvate_reassimilation`` sets the re-assimilation magnitude; with
    #: ``k_pyruvate_excretion`` its *ratio* sets the frozen residual. ``K_sugar_uptake`` is shared
    #: with the fermentative-uptake flux this co-metabolically tracks (so it dies at dryness). Their
    #: tiers cap the outputs via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_pyruvate_reassimilation", "K_sugar_uptake")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        pyruvate = max(float(y[schema.slice("pyruvate")][0]), 0.0)
        if pyruvate <= 0.0:  # nothing to re-assimilate
            return d
        flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:  # dryness / no viable yeast ⇒ re-assimilation stops, pool is frozen
            return d
        loss = params["k_pyruvate_reassimilation"] * flux * pyruvate  # mass loss of pyruvate
        r = loss / M_PYRUVATE  # molar turnover: 1 pyruvate → 1 ethanol + 1 CO₂ (C3 → C2 + C1)
        d[schema.slice("pyruvate")] = -loss
        d[schema.slice("E")] = r * M_ETHANOL
        d[schema.slice("CO2")] = r * M_CO2
        return d
