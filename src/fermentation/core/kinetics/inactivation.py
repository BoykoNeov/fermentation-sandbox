"""Ethanol-driven cell inactivation — the cumulative viability brake on the tail.

The mechanism that actually sets a wine fermentation's *timescale*. As ethanol
accumulates it does not merely slow living cells (a reversible, instantaneous
effect) — it progressively *kills* them, and dead cells never ferment again. This
irreversible, cumulative loss of catalytic biomass is what decelerates the back
half of a ferment and lets a 24 Brix must finish in ~10-14 days rather than
either racing to dryness in ~6 days (no brake) or stalling forever against an
instantaneous wall (see decision D-13 and :mod:`~fermentation.core.kinetics.\
inhibition`).

**Sourced from the keystone wine model (Coleman, Fish & Block 2007).** Their model
splits biomass into active cells ``X_A`` and tracks an ethanol-proportional
inactivation rate (their eqs 2 and 7)::

    dX_A/dt = mu*X_A - k_d*X_A          (eq 2)
    k_d     = k'_d * E                  (eq 7)

so cells inactivate at a specific rate proportional to the ethanol they sit in.
``k'_d`` is one of only three parameters Coleman found strongly temperature
dependent (its Table A2 regression is *quadratic* in T, far steeper than growth or
uptake — which is exactly why high temperature drives stuck fermentations).

**Two-pool representation (decision D-13).** We keep ``X`` as the *viable*
biomass it has always been (growth and uptake are catalysed by ``X``), and add an
inactivated pool ``X_dead``. Inactivation moves mass from one pool to the other at
equal rate::

    r        = k'_d * E * X             [g/L/h]
    dX/dt   += -r
    dX_dead += +r

Because the two pools have identical elemental composition, this transfer is
carbon- and nitrogen-neutral *by construction* — a gram leaving ``X`` arrives in
``X_dead`` with the same ``f_C``/``f_N``, so ``total_carbon`` and ``total_nitrogen``
(which weight both pools by those fractions) are untouched by death. The viability
loss is genuinely *stateful*: it is the integral of ethanol damage over time, held
in the ``X``/``X_dead`` split, not a function of the instantaneous state — which is
why it finishes where the reversible Luong wall stalled.

This Process *replaces* the Luong ethanol wall as the validated core's
ethanol-brake (D-13): with a mechanistic, cumulative inactivation in hand, also
multiplying uptake by an instantaneous ``(1 - E/E_max)`` factor would double-count
ethanol toxicity. :class:`~fermentation.core.kinetics.inhibition.EthanolInhibition`
is retained for optional/strain use but is no longer wired into the default media.

Tier: **plausible** — a sound, literature-sourced mechanism (Coleman 2007, the same
fit our growth/uptake constants come from). It sets the wine timescale the §2.2
dryness benchmark now passes; that confirms the plausible tier is earned but does
deliberately **not** promote it to validated, which awaits validation against
independent *measured* time-series (handoff D-C: none exist yet; decision D-17).
``k'_d`` itself is sourced (Coleman Table A2, evaluated at the 20 C T_ref).
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


class EthanolInactivation(Process):
    """Ethanol-proportional inactivation of viable biomass (Coleman 2007).

    Moves viable biomass ``X`` into the inactivated pool ``X_dead`` at the specific
    rate ``k_d = k'_d·E`` (``r = k'_d·E·X`` g/L/h). Carbon/nitrogen neutral because
    the two pools share one composition; see the module docstring and decision D-13.
    """

    name = "ethanol_inactivation"
    tier = Tier.PLAUSIBLE
    #: Viable biomass leaves ``X``; the same mass enters the inactivated pool
    #: ``X_dead``. Declaring both keeps the transfer inside the ``touches`` contract.
    touches = ("X", "X_dead")
    #: The ethanol-sensitivity constant (Coleman eq 7). Its tier caps the tier of
    #: ``X``/``X_dead`` via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_prime_d",)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        x = max(float(y[schema.slice("X")][0]), 0.0)
        e = max(float(y[schema.slice("E")][0]), 0.0)
        # Clamp X and E to >= 0 before the product: a negative solver excursion
        # would otherwise flip the inactivation sign (resurrecting dead cells /
        # creating viable biomass). Mirrors the guards in the other Processes.
        if x <= 0.0 or e <= 0.0:  # no ethanol or no cells -> no inactivation
            return d
        rate = params["k_prime_d"] * e * x  # [g/L/h] viable biomass inactivated
        d[schema.slice("X")] = -rate
        d[schema.slice("X_dead")] = rate
        return d
