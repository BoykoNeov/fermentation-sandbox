"""Hop bittering — alpha-acid isomerization in the boil, and iso-alpha loss in the ferment.

The §3.3 additive-with-a-clear-mechanism beat (decision D-64). Bitterness comes from
*iso*-alpha-acids (isohumulones), which do not pre-exist in the hop: they are made in the
boil by thermal isomerization of the hop's alpha-acids (humulones), and then partly lost
again during fermentation. Two physically distinct regimes, handled in the two places they
belong:

    1. THE BOIL (~373 K, ~60-90 min, PRE-fermentation, no yeast). A CONSECUTIVE first-order
       reaction, both rate constants measured by Malowicki & Shellhammer 2005 (J. Agric.
       Food Chem. 53(11):4434-4439, doi:10.1021/jf0481296) over 90-130 C in a model wort:

           alpha --k1--> iso-alpha --k2--> uncharacterized degradation products

       This module gives the CLOSED-FORM solution (:func:`iso_alpha_fraction`), which the
       compile seam evaluates once per hop addition and sums — NOT a boil ODE phase. Running
       the boil through the fermentation integrator would drive the (yeast-free, sugar-full)
       wort at 373 K, which is meaningless. The boil is a wort-side input, exactly like the
       measured ``initial_ph`` the compile seam wires at t=0; only its *result* (the
       iso-alpha delivered to the fermenter) enters the state.

    2. FERMENTATION (the engine's native regime). :class:`IsoAlphaAcidLoss` removes iso-alpha
       by adsorption onto viable yeast and scavenging into the krausen/foam — the standard
       account of the ~5-20% wort-to-beer bitterness drop. This is the dynamic content of the
       beat and the reason hops touch the ODE at all.

**Off the carbon ledger (decision D-64).** Iso-alpha-acids are exogenous (they arrive via
hops, at the mg/L scale) and never touch ``S``/``E``/``CO2``/``N``. Like dosed SO2 (D-22)
they are absent from ``total_carbon``/``total_nitrogen`` (an unreferenced state slot gets
weight 0 in :mod:`fermentation.validation.conservation`), so the whole beat leaves the carbon
invariant BYTE-FOR-BYTE unchanged. The fermentation loss is adsorptive removal of hop-derived
mass, not a conversion within the fermentation carbon budget.

**The closed form.** For a consecutive first-order reaction A -k1-> B -k2-> C with only A
present initially, the intermediate B (iso-alpha) obeys::

    [B](t) / [A]_0 = k1 / (k2 - k1) * (exp(-k1 t) - exp(-k2 t))     (k1 != k2)

(the k1 == k2 degenerate limit, [A]_0 * k1 t exp(-k1 t), is guarded for completeness though
hop k1 != k2 at every temperature). At 100 C / 60 min this is ~0.48 — about 48% of the
dissolved alpha is present as iso-alpha, still on the RISING limb (k2 < k1, so the iso-alpha
peak is ~2-3 h out; a 60-90 min boil has not reached it), matching brewing practice.

**Tiers (decision D-64, split; let D-1 derive it — do not assert one tier for the readout).**
The boil rate constants are SOURCED/measured (Malowicki), so the end-of-boil iso-alpha
fraction is PLAUSIBLE. But the finished-beer IBU also depends on the SPECULATIVE
``hop_utilization_efficiency`` (the kettle->fermenter non-fermentation losses, applied at the
compile seam) and this SPECULATIVE :class:`IsoAlphaAcidLoss`, so parameter-tier propagation
caps the finished ``iso_alpha`` output at speculative. The Malowicki reaction is measured but
its mapping to real wort (extraction, gravity, hop form) is an honest-mapping step, so even
the boil kinetics stay plausible rather than validated. Utilization's gravity-dependence,
dry-hop/whirlpool bitterness, and hop-form effects are documented v1 deferrals (see hops.yaml).
"""

from __future__ import annotations

import math
from collections.abc import Mapping

from fermentation.core.kinetics.arrhenius import GAS_CONSTANT
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


def boil_rate_constants(boil_temp_k: float, params: Mapping[str, float]) -> tuple[float, float]:
    """Malowicki 2005 isomerization/degradation rate constants ``(k1, k2)`` [1/min] at ``T``.

    Absolute Arrhenius ``k = A * exp(-E_a / (R T))`` (T in Kelvin) — the boil kinetics have a
    measured pre-exponential, so this uses the absolute form, not the reference-anchored
    :func:`~fermentation.core.kinetics.arrhenius.arrhenius_factor` the fermentation-rate
    modifiers use (D-11). ``k1`` isomerizes alpha -> iso-alpha; ``k2`` degrades iso-alpha.
    """
    k1 = params["A_iso"] * math.exp(-params["Ea_iso"] / (GAS_CONSTANT * boil_temp_k))
    k2 = params["A_iso_degradation"] * math.exp(
        -params["Ea_iso_degradation"] / (GAS_CONSTANT * boil_temp_k)
    )
    return k1, k2


def iso_alpha_fraction(
    boil_minutes: float, boil_temp_k: float, params: Mapping[str, float]
) -> float:
    """Fraction of dissolved alpha-acid present as iso-alpha at the end of a boil.

    The closed-form consecutive-first-order intermediate ``k1/(k2-k1)·(e^{-k1 t}-e^{-k2 t})``
    (the degenerate ``k1==k2`` limit ``k1 t e^{-k1 t}`` is guarded for completeness). Pure and
    parameter-driven; the compile seam multiplies it by the dissolved alpha concentration and
    ``hop_utilization_efficiency`` and sums over the hop schedule. A zero-length boil returns 0
    (no isomerization), never negative.
    """
    if boil_minutes <= 0.0:
        return 0.0
    k1, k2 = boil_rate_constants(boil_temp_k, params)
    t = boil_minutes
    if abs(k2 - k1) < 1e-12:  # degenerate limit (not reached by real hop kinetics)
        return float(k1 * t * math.exp(-k1 * t))
    frac = k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))
    return float(max(frac, 0.0))


class IsoAlphaAcidLoss(Process):
    """Fermentation-time loss of iso-alpha-acids to yeast adsorption and krausen scavenging.

    ``d(iso_alpha)/dt = -k_iso_alpha_loss · X_viable · iso_alpha`` — first-order in both the
    iso-alpha present and the VIABLE biomass ``X`` (not ``X_dead``), the standard mechanistic
    account of the ~5-20% wort-to-beer bitterness drop: iso-alpha adsorbs onto growing yeast
    cell walls and is scavenged into the foam. Gating on live ``X`` means a crashed or racked
    ferment stops losing bitterness, and an unhopped beer (``iso_alpha`` = 0) contributes
    exactly zero — so with no hops the term is inert (the compile seam additionally DISABLES
    this Process when no hops are scheduled, keeping the empty ``iso_alpha`` slot's tier
    VALIDATED and paying no flux, the MLF/Brett isolability pattern).

    Touches ``iso_alpha`` only — no ``S``/``E``/``CO2``/``N`` coupling — so it is OFF the
    carbon ledger and cannot perturb ``total_carbon`` (iso-alpha is exogenous hop-derived
    mass, adsorptively removed, not fermentation carbon). Both ``iso_alpha`` and ``X`` are
    clamped >= 0 against solver undershoot. Tier **speculative** (the loss-rate magnitude is
    an author estimate; the yeast-adsorption mechanism is the sourced-standard account).
    Beer-only (wine has no ``iso_alpha`` slot).
    """

    name = "iso_alpha_acid_loss"
    tier = Tier.SPECULATIVE
    touches = ("iso_alpha",)
    reads: tuple[str, ...] = ("k_iso_alpha_loss",)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        iso = max(float(y[schema.slice("iso_alpha")][0]), 0.0)
        if iso <= 0.0:  # unhopped (or fully stripped) — nothing to lose
            return d
        x_viable = max(float(y[schema.slice("X")][0]), 0.0)
        if x_viable <= 0.0:  # no viable yeast ⇒ no adsorption (bitterness frozen)
            return d
        d[schema.slice("iso_alpha")] = -params["k_iso_alpha_loss"] * x_viable * iso
        return d
