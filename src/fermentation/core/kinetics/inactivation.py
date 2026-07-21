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

**The tolerance supplement (decision D-129).** Coleman's ``k_d = k'_d·E`` is *linear
and unbounded* in ethanol by construction — verified in the source (Appl. Environ.
Microbiol. 73(18), eq 7), which was fit and validated only over 265-300 g/L initial
sugar (~26.5-30 °Brix). Inside that envelope the brake sets the timescale correctly
and real high-tolerance strains *do* finish (a 28 °Brix must dries at ~140 g/L EtOH,
below EC-1118's rated 142 g/L). *Extrapolated past it*, the linear death has no
ceiling: a 36 °Brix must marches to ~187 g/L EtOH (~23.7% ABV) instead of *sticking*
sweet like a real dessert/icewine must. :class:`EthanolToleranceDeath` supplies the
missing ceiling as a *super-linear* addition to the *same death mechanism* (staying
within D-13's "cumulative death, not an instantaneous uptake wall" ruling, so no
double-count): ``Φ(E) = k_d2·max(E − E_tol, 0)²`` extra specific death, so viability
collapses once ethanol overshoots the strain tolerance ``E_tol``
(``ethanol_tolerance``). Anchored to that *sourced* tolerance for **where** the
ceiling sits and to Schenk et al. 2014 (arXiv:1412.6068) for the tolerance-gated
quadratic **form**; the curvature ``k_d2`` (speculative) sets only **how sharply** it
arrests — the D-129 sweep shows the sticking threshold and arrest-ABV plateau are
governed by the sourced ``E_tol`` and shape-insensitive across 2.5 orders of ``k_d2``,
with only the residual-sugar *amount* speculative. The ``max(·, 0)²`` form is exactly
zero below ``E_tol`` (so an in-envelope ferment is **byte-for-byte** the Coleman-only
core — structural isolability, prime directive #3) and C¹-smooth at ``E_tol`` (no BDF
derivative kink). It is tier **speculative** and kept in its own isolable tuple.

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


class EthanolToleranceDeath(Process):
    """Super-linear viability collapse once ethanol overshoots the strain tolerance.

    The ceiling the linear Coleman brake lacks (decision D-129). Extra specific death
    ``Φ(E) = k_d2·max(E − E_tol, 0)²`` moved from viable ``X`` into ``X_dead`` at
    ``r = Φ·X`` g/L/h, with ``E_tol = ethanol_tolerance`` (the *sourced* strain rating)
    and ``k_d2 = k_d2_ethanol_tolerance_death`` (speculative curvature). This is an
    addition to the *same* cumulative death mechanism as
    :class:`EthanolInactivation` — not an instantaneous uptake wall — so it does not
    double-count ethanol toxicity (D-13). Carbon/nitrogen neutral by construction: the
    two pools share one elemental composition, so a gram leaving ``X`` arrives in
    ``X_dead`` at the same ``f_C``/``f_N`` (see the module docstring).

    The ``max(·, 0)²`` form is **exactly zero** for ``E ≤ E_tol``: any ferment that
    finishes within Coleman's validated envelope (E stays below tolerance) is
    byte-for-byte the Coleman-only core, and the term is C¹-smooth at ``E_tol`` (no BDF
    kink). Only when ethanol overshoots the rated tolerance does viability collapse
    super-linearly and the must *stick* with residual sugar. Tier **speculative** — the
    quadratic form is sourced (Schenk et al. 2014) and ``E_tol`` is sourced, but ``k_d2``
    (the arrest sharpness) is an author estimate; via parameter-tier propagation (D-1)
    it caps the tier of the biomass pools this writes when the Process is active.
    """

    name = "ethanol_tolerance_death"
    tier = Tier.SPECULATIVE
    #: Viable biomass leaves ``X``; the same mass enters ``X_dead`` (same transfer,
    #: same neutrality, as :class:`EthanolInactivation`).
    touches = ("X", "X_dead")
    #: The speculative curvature and the sourced strain tolerance. Their tiers cap the
    #: tier of ``X``/``X_dead`` when this Process is active (D-1).
    reads: tuple[str, ...] = ("k_d2_ethanol_tolerance_death", "ethanol_tolerance")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        x = max(float(y[schema.slice("X")][0]), 0.0)
        e = max(float(y[schema.slice("E")][0]), 0.0)
        # Clamp X and E to >= 0 before the term (mirrors EthanolInactivation): a
        # negative solver excursion must not flip the death sign or resurrect cells.
        over = e - params["ethanol_tolerance"]
        if x <= 0.0 or over <= 0.0:  # at/below tolerance the ceiling is inert
            return d
        phi = params["k_d2_ethanol_tolerance_death"] * over * over  # [1/h] specific rate
        rate = phi * x  # [g/L/h] viable biomass killed past tolerance
        d[schema.slice("X")] = -rate
        d[schema.slice("X_dead")] = rate
        return d
