"""Acetaldehyde — the main-pathway intermediate as a transient ethanol-carbon buffer.

The second §3.2 aroma beat after diacetyl (decision D-27). Acetaldehyde (ethanal,
CH₃CHO) is the obligate intermediate on the *main* alcoholic-fermentation pathway

    sugar → … → pyruvate → acetaldehyde → ethanol
                              (produced)   (reduced)

and, like diacetyl (D-26) but on the trunk rather than a side branch, it is
**produced then reabsorbed** — it accumulates to an early peak during vigorous
fermentation (the "green apple" note) and is then reduced to ethanol by viable
yeast. It is also the principal carbonyl that binds SO₂, so building it as real
state unlocks the deferred free/bound-SO₂ split (a separate readout beat, D-22).

**Carbon — a transient buffer that borrows ethanol carbon, NOT a draw from sugar
(the load-bearing modelling choice, decision D-27).** This is where acetaldehyde
differs from every prior byproduct. Esters, fusels and the VDK pools land their
carbon in *side* pools genuinely removed from ethanol, so routing their carbon out
of ``S`` (option a1, D-19) is benchmark-neutral. Acetaldehyde's product is
**ethanol ``E`` itself** — and the uptake Process *already* performs the complete
(lumped) sugar → ethanol + CO₂ conversion, which implicitly includes this
intermediate. Drawing *fresh* sugar → acetaldehyde → *new* ethanol would therefore
be a **second, parallel** sugar→ethanol pathway: net-new ethanol that inflates ABV
and the realised yield by an amount that *scales with pool turnover* (cumulative
acetaldehyde produced, not its peak). That is a real double-count, not a trace.

The faithful model de-lumps the existing step instead of duplicating it. Because
acetaldehyde and ethanol are **both two-carbon**, the reduction acetaldehyde →
ethanol is a mole-for-mole C2 → C2 transfer. So:

* :class:`AcetaldehydeProduction` **holds back** a transient slice of the ethanol
  the uptake Process just made — reclassifying it as the true intermediate:
  ``d(acetaldehyde)/dt = +r`` with the equal-carbon ``d(E)/dt = −r·M_ethanol/
  M_acetaldehyde``. No fresh sugar, no CO₂.
* :class:`AcetaldehydeReduction` **returns** it: ``d(acetaldehyde)/dt = −L`` with
  ``d(E)/dt = +L·M_ethanol/M_acetaldehyde``.

``total_carbon`` (which now weights ``acetaldehyde``, see
:mod:`fermentation.validation.conservation`) closes to **machine precision** through
the whole produce-then-reabsorb course, touching **neither ``S`` nor ``CO2``**. The
``E`` **endpoint reconverges** to the no-acetaldehyde core to **relative ~1e-8** once
the pool is reduced back — so the §2.2 ABV / realised-yield / CO₂ benchmarks are
preserved to far below any tolerance; during the ferment ``E`` carries a tiny
transient dip (carbon held as acetaldehyde) that recovers as reduction completes, and
that dip feeds a *second-order* ~1e-4 relative path perturbation of the rest of the
core through the ``E``→viability coupling (see Isolability below). Mass carries a small gap (the
oxidative production and the reduction move NAD(P)H that the model does not track) —
carbon is the invariant, exactly as for the diacetyl reduction (D-26) and beer's
hydrolysis water (D-8).

This is the owner's call (decision D-27), chosen over the draw-from-sugar stand-in
the D-26 forward note anticipated (that note applied the ester/fusel template before
anyone noticed acetaldehyde's product is ``E``, not a side pool). It is *more*
faithful, not merely benchmark-safe: acetaldehyde genuinely **is** obligate in-transit
ethanol carbon, so borrowing from ``E`` asserts exactly the right provenance, whereas a
sugar draw would assert a parallel pathway that does not exist.

**Why the early peak emerges.** Production tracks the fermentative flux (shared
``K_sugar_uptake``), so it is strong during active fermentation and stops at dryness.
Reduction is gated on **viable** ``X`` (not ``X_dead``) with **no flux term**, so it
keeps clearing acetaldehyde during the rest and after sugar is gone, but stops dead
once the yeast is crashed / racked / ethanol-inactivated — the same structural pair as
the diacetyl decarb/reduction (D-26). Acetaldehyde thus rises while the flux outruns
the still-building reductive capacity, peaks, and is drawn back down to a low residual
as fermentation slows; a crash before it is cleared **strands** it (and leaves the
borrowed ethanol carbon parked as acetaldehyde rather than returned to ``E``).

**Isolability (prime directive #3).** Both Processes live in their own
``_ACETALDEHYDE_PROCESSES`` tuple (:mod:`fermentation.core.media`), so a ProcessSet
built without them is the prior core. Acetaldehyde is intrinsic yeast metabolism (not
a dosed organism like MLF), so it is wired into *both* media and runs on every default
ferment — like esters and the VDK pools. It touches **only ``acetaldehyde`` and ``E``**
(never ``S``/``CO2``/``N``/``X``) at the derivative level, so ``dS``/``dCO2``/``dN`` are
byte-for-byte identical given the same state — the same derivative-level isolability the
D-19 sugar draws claim. The one integrated coupling is second-order: because ``E`` feeds
the ethanol-inactivation viability brake, the transient ``E`` dip perturbs the *path* by a
tiny ~1e-4 relative amount, and the ``E`` **endpoint** reconverges to the buffer-off core
to relative ~1e-8 (the pool fully reduces back), so the §2.2 ABV / realised-yield / CO₂
benchmarks are preserved to far below any tolerance.

**Tiers.** Every rate constant is an order-of-magnitude estimate, so both Processes are
**speculative**; parameter-tier propagation (D-1) caps the ``acetaldehyde`` output tier
at speculative regardless. Because :class:`AcetaldehydeProduction` is the first
always-on speculative Process to *write* ``E``, the *structural* ``tier_of("E")`` drops
PLAUSIBLE → SPECULATIVE — but the param-aware tier users see was *already* speculative
(the uptake Process reads speculative params), so there is no headline change, exactly
the honest ``CO2`` consequence recorded for D-26. Production is held temperature-flat (a
documented simplification, like the acetolactate excretion, D-26); the enzymatic reduction
carries an Arrhenius factor.

**SO₂ protection — the reduction reads free, not total (decision D-47).** The D-28 free/bound
SO₂ split is fed back into this RHS: the acetaldehyde-bisulfite adduct is protected from
alcohol dehydrogenase (literature: bound acetaldehyde "could not be metabolized by yeast during
fermentation; only free acetaldehyde could impact metabolism"), so :class:`AcetaldehydeReduction`
reduces only the unbound share (:func:`fermentation.core.acidbase.free_acetaldehyde`). Dosed SO₂
therefore **locks in** acetaldehyde — a sulfited wine strands a residual acetaldehyde pool
(≈ the SO₂ molar amount, capped at what is present) and its free SO₂ stays depressed, both
*emergent* from the binding equilibrium tracking the acetaldehyde state. This intentionally retires
the D-22/D-28 "SO₂ is readout-only" invariant *for sulfited runs*; an un-dosed run is byte-for-byte
the D-27 core (the ``so2_total > 0`` guard is exact) and no §2.2 benchmark doses SO₂. CAVEAT
(speculative): bound acetaldehyde is treated inert-to-ADH — real adduct slowly dissociates and
degrades over months, so the stranding is an *upper bound* on persistence; at field (sub-
stoichiometric) doses the model reproduces the observed ~0.76× degradation-rate slowdown.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.acidbase import (
    SO2_STATE_KEY,
    free_acetaldehyde,
    ph_of_state,
)
from fermentation.core.chemistry import M_ACETALDEHYDE, M_ETHANOL
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import fermentative_flux_shape
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Ethanol ⇄ acetaldehyde are both C2, so the reduction (and its inverse, the
#: production "borrow") is a mole-for-mole transfer weighted by this molar-mass ratio —
#: exactly the ``M_BUTANEDIOL/M_DIACETYL`` shape of the diacetyl reduction (D-26). Carbon
#: closes because ``M_acet·cf_acet == M_eth·cf_eth·(M_acet/M_eth) == 2·M_C`` per mole.
_ETHANOL_PER_ACETALDEHYDE = M_ETHANOL / M_ACETALDEHYDE


class AcetaldehydeProduction(Process):
    """Acetaldehyde formation — a transient slice of ethanol carbon held back as the
    true main-pathway intermediate.

    ``d(acetaldehyde)/dt = +r`` and ``d(E)/dt = −r·M_ethanol/M_acetaldehyde`` with
    ``r = k_acetaldehyde · X · S_total/(K_sugar_uptake + S_total)`` (the shared
    fermentative-flux shape). It does **not** make acetaldehyde from sugar: it borrows
    carbon from the ethanol the uptake Process is depositing on the same flux and
    reclassifies it as acetaldehyde (decision D-27). Because ethanol and acetaldehyde are
    both two-carbon, the borrow is carbon-exact and touches only ``acetaldehyde`` and
    ``E`` — never ``S``/``CO2`` — so ``dS``/``dCO2``/``dN`` are byte-for-byte given the same
    state; the only perturbation to the rest of the core is the transient ``E`` dip the
    reduction then repays (and its tiny second-order feedback through viability).

    Coupled to the flux (so production stops at dryness) and held **temperature-flat**
    (a documented v1 simplification, like the α-acetolactate excretion, D-26): the
    interesting temperature dependence lives in the enzymatic reduction, and no benchmark
    asserts an acetaldehyde temperature direction in v1. ``E`` is not clamped here — net
    ``dE`` stays positive throughout because the uptake ethanol rate rides the *same*
    flux with a far larger coefficient — but the borrow itself cannot drive ``E`` below
    the uptake's own deposit, so no guard is needed. Tier **speculative** (rate magnitude
    estimate).
    """

    name = "acetaldehyde_production"
    tier = Tier.SPECULATIVE
    touches = ("acetaldehyde", "E")
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks;
    #: ``k_acetaldehyde`` sets the borrow magnitude (and, with ``k_acet_reduction``, the
    #: quasi-steady peak). Their tiers cap ``acetaldehyde``'s output tier via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_acetaldehyde", "K_sugar_uptake")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        rate = params["k_acetaldehyde"] * flux  # mass rate of acetaldehyde formation
        d[schema.slice("acetaldehyde")] = rate
        d[schema.slice("E")] = -rate * _ETHANOL_PER_ACETALDEHYDE  # borrow C2 from ethanol
        return d


class AcetaldehydeReduction(Process):
    """Enzymatic acetaldehyde → ethanol by *viable* yeast — returns the borrowed carbon.

    ``d(acetaldehyde)/dt = −L`` and ``d(E)/dt = +L·M_ethanol/M_acetaldehyde`` with the
    mass loss ``L = k_acet_reduction · X · f(T) · [acetaldehyde]`` and ``f(T) =
    arrhenius_factor(T, E_a_acet_reduction, T_ref)``. Alcohol dehydrogenase reduces
    acetaldehyde to ethanol (the last step of the main pathway), a mole-for-mole C2 → C2
    transfer, so weighting the ethanol gain by ``M_ethanol/M_acetaldehyde`` returns exactly
    the carbon :class:`AcetaldehydeProduction` borrowed — closing ``total_carbon`` to
    machine precision. Touches only ``acetaldehyde`` and ``E`` (never ``S``/``CO2``). Mass
    carries a small gap (the reduction consumes untracked NAD(P)H) — carbon is the invariant.

    **Gated on VIABLE biomass ``X`` (not ``X_dead``), with NO fermentative-flux term
    (decision D-27), mirroring the diacetyl reduction (D-26).** Live-yeast gating makes the
    reduction stop the moment the yeast is crashed / racked / ethanol-inactivated, so an
    early crash **strands** acetaldehyde with its borrowed ethanol carbon un-returned; the
    absence of a flux term lets the reduction keep clearing acetaldehyde *after* sugar is
    gone (the produce-then-reabsorb tail). Together with the flux-linked production, that
    pair makes acetaldehyde rise to an early peak and then be drawn back down as
    fermentation slows. ``acetaldehyde`` and ``X`` are clamped ≥ 0 against solver
    undershoot. Tier **speculative** (rate magnitude estimate).

    **SO₂-bound acetaldehyde is protected from ADH (decision D-47).** When SO₂ is dosed
    (``so2_total > 0``) the loss reads the *free* (unbound) acetaldehyde
    (:func:`fermentation.core.acidbase.free_acetaldehyde`) rather than the total — the
    hydroxysulphonate adduct cannot be reduced — so SO₂ locks acetaldehyde in. The guard is
    exact: an unsulfited run pays no per-RHS pH ``brentq`` and is byte-for-byte the D-27 core.
    ``touches`` is unchanged (still only ``acetaldehyde``/``E``); ``reads`` is unchanged too —
    the SO₂/pH params are read *inside* :func:`free_acetaldehyde`/:func:`ph_of_state` and the
    output is already speculative, so declaring them would not move any tier (the MLF-gate
    precedent, D-39).
    """

    name = "acetaldehyde_reduction"
    tier = Tier.SPECULATIVE
    touches = ("acetaldehyde", "E")
    #: ``k_acet_reduction`` sets the enzymatic reduction magnitude; ``E_a_acet_reduction``
    #: and ``T_ref`` set the temperature shape. Their tiers cap the outputs via
    #: parameter-tier propagation (D-1). Reads viable ``X`` from state.
    reads: tuple[str, ...] = ("k_acet_reduction", "E_a_acet_reduction", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        acetaldehyde = max(float(y[schema.slice("acetaldehyde")][0]), 0.0)
        if acetaldehyde <= 0.0:  # nothing to reduce
            return d
        x_viable = max(float(y[schema.slice("X")][0]), 0.0)
        if x_viable <= 0.0:  # no viable yeast ⇒ no reduction (acetaldehyde is stranded)
            return d
        # SO₂-bound acetaldehyde is protected from ADH (decision D-47): reduce only the *free*
        # share. The ``so2_total > 0`` guard is EXACT — an unsulfited run pays no per-RHS pH
        # ``brentq`` and its contribution is byte-for-byte the D-27 core (the MLF/Brett SO₂-gate
        # isolability idiom). When dosed, SO₂ binds acetaldehyde near-stoichiometrically and
        # *locks it in*: the reducible pool shrinks toward the excess of acetaldehyde over SO₂.
        reducible = acetaldehyde
        if SO2_STATE_KEY in schema and float(y[schema.slice(SO2_STATE_KEY)][0]) > 0.0:
            reducible = free_acetaldehyde(y, schema, params, ph_of_state(y, schema, params))
            if reducible <= 0.0:  # all acetaldehyde bound ⇒ nothing ADH can reach this step
                return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_acet_reduction"], params["T_ref"])
        loss = params["k_acet_reduction"] * x_viable * f_t * reducible  # mass loss (free only)
        d[schema.slice("acetaldehyde")] = -loss
        d[schema.slice("E")] = loss * _ETHANOL_PER_ACETALDEHYDE  # return C2 to ethanol
        return d
