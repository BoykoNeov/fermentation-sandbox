"""Arrhenius temperature dependence — a multiplicative rate modifier.

Every kinetic rate constant rises with temperature. The standard description is
the Arrhenius law ``k(T) = A·exp(-E_a / R T)``: a higher activation energy
``E_a`` makes the rate more temperature-sensitive. Like ethanol inhibition this
mechanism *scales* an existing flux rather than *adding* one, so it is a
:class:`~fermentation.core.process.RateModifier` (reusing the hook D-10 built for
:class:`~fermentation.core.kinetics.inhibition.EthanolInhibition`), not a summed
:class:`~fermentation.core.process.Process`. Decision **D-11**.

This module also holds :class:`ColemanQuadraticDeathTemperature` (decision D-57),
which scales :class:`~fermentation.core.kinetics.inactivation.EthanolInactivation`.
It is *not* Arrhenius-shaped (Coleman's own ``k'_d`` regression is quadratic in
°C, not exponential in 1/T) but lives here because its role — a per-rate,
reference-anchored temperature factor on a primary-fermentation Process — is the
same concern this module owns.

**Reference-anchored form (no separate pre-exponential).** Rather than carry an
absolute pre-exponential ``A``, the factor is normalised to a reference
temperature ``T_ref`` — the temperature the rate constants were measured at::

    f(T) = exp( -(E_a / R) · (1/T - 1/T_ref) )

At ``T = T_ref`` the factor is exactly 1, so the *measured* rate constant
(``mu_max`` / ``q_sugar_max``) is used unscaled; above ``T_ref`` the factor
exceeds 1 (faster), below it the factor is < 1 (slower). This matters because the
measured constants *already* encode ``A·exp(-E_a / R T_ref)``: adding a standalone
``A`` would double-book the pre-exponential and could silently disagree with the
rate constant it multiplies. So only ``E_a`` and ``T_ref`` are parameters; ``A``
is not (D-11).

**Per-rate, not one shared E_a.** Growth and fermentation differ in their
temperature sensitivity, so this modifier is *parameterised*: each instance names
the Process it scales and the activation-energy Parameter it reads. The wine
configuration uses two — :meth:`for_growth` (``E_a_growth`` on
:class:`~fermentation.core.kinetics.growth.GrowthNitrogenLimited`) and
:meth:`for_uptake` (``E_a_uptake`` on
:class:`~fermentation.core.kinetics.uptake.SugarUptakeToEthanolCO2`) — sharing one
``T_ref``. (Per-rate activation energies are the documented intent:
``milestone-1-context.md`` "Arrhenius A + E_a *per rate*"; D-11.)

**Conservation is automatic, no clamp needed.** ``exp`` is always positive, so the
factor scales each targeted Process's whole contribution vector by a single
positive scalar — every balance that Process respects is preserved (a uniformly
faster carbon-neutral flux is still carbon-neutral). Unlike the wall-type
inhibition form there is no regime where the factor could turn negative, so no
clamp is required (a defensive one would be inconsistent noise).

**Isothermal in M1.** No M1 Process drives ``T``, so within a run the factor is
constant; its role is to make *different-temperature* runs differ in rate (the
directional "lower T → slower" check). It reads ``T`` from the state vector
(Kelvin, canonical per D-3) rather than from params, so it is already correct for
the non-isothermal temperature dynamics of a later tier.

Tier: **plausible** — the Arrhenius law is textbook, literature-standard kinetics.
It stays plausible after the §2.2 promotion sweep precisely because the M1
benchmarks are isothermal at ``T_ref`` (``f = 1``), so they never exercise this
modifier — an untested mechanism cannot be promoted on their strength (decision
D-17). The placeholder ``E_a`` values are *speculative*; parameter-tier propagation
(D-1) caps the scaled outputs accordingly.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

from fermentation.core.kinetics.growth import GrowthNitrogenLimited
from fermentation.core.kinetics.inactivation import EthanolInactivation
from fermentation.core.kinetics.uptake import SugarUptakeToEthanolCO2
from fermentation.core.process import RateModifier
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Universal gas constant [J/(mol·K)]. Exact by SI definition (2019 redefinition:
#: ``R = N_A · k_B`` with both defining constants exact); CODATA value
#: 8.314462618…. A universal physical constant, not a stoichiometric one, so it
#: lives here local to the modifier rather than in ``core.chemistry`` (whose scope
#: is molar masses / carbon counts) or the provenance store (which is for
#: empirical, uncertain quantities) — same code-with-citation rule as D-3/D-8.
GAS_CONSTANT = 8.314462618


def arrhenius_factor(temp: float, e_a: float, t_ref: float) -> float:
    """Reference-anchored Arrhenius factor ``exp(-(E_a/R)·(1/T - 1/T_ref))``.

    The dimensionless temperature scaling that is exactly 1 at ``T = T_ref`` (so a
    rate constant measured at ``T_ref`` is used unscaled), > 1 above and < 1 below.
    ``temp``/``t_ref`` are in Kelvin, ``e_a`` in J/mol. Always positive (``exp``),
    so multiplying a conserving flux by it cannot flip a sign or break a balance.

    Shared single source of truth for the Arrhenius shape: the multiplicative
    :class:`ArrheniusTemperature` modifier (which scales growth/uptake) and the
    additive Milestone-2 byproduct Processes (which *embed* their own, steeper
    temperature sensitivity) both call this, so the law is written once (D-11).
    """
    return float(math.exp(-(e_a / GAS_CONSTANT) * (1.0 / temp - 1.0 / t_ref)))


class ArrheniusTemperature(RateModifier):
    """Arrhenius temperature scaling of one kinetic Process's flux.

    ``factor = exp(-(E_a / R)·(1/T - 1/T_ref))`` (``E_a`` the per-instance
    activation-energy parameter, ``T_ref`` the shared reference temperature, ``T``
    read from the state vector in Kelvin). Multiplied onto the named Process's
    whole contribution by :class:`~fermentation.core.process.ProcessSet`.

    Parameterised because growth and fermentation differ in temperature
    sensitivity: construct one per rate (see :meth:`for_growth` / :meth:`for_uptake`)
    rather than sharing a single ``E_a``. ``name`` is set per instance (not a class
    attribute) so two instances do not collide in a ``ProcessSet``'s shared
    Process/modifier name space.
    """

    tier = Tier.PLAUSIBLE

    def __init__(
        self, *, name: str, modifies: str | tuple[str, ...], activation_energy: str
    ) -> None:
        #: Per-instance so multiple Arrhenius modifiers coexist in one ProcessSet.
        self.name = name
        #: One or more Processes scaled by this temperature factor. A single name is
        #: the common case (one flux); a tuple lets one Arrhenius factor scale several
        #: Processes that share a temperature sensitivity — e.g. growth *and* the
        #: amino-acid swap, so the swap's refunds track growth's realised draw (D-32).
        self.modifies = (modifies,) if isinstance(modifies, str) else tuple(modifies)
        #: The activation-energy Parameter varies per instance; ``T_ref`` is shared.
        self._e_a_param = activation_energy
        self.reads = (activation_energy, "T_ref")

    @classmethod
    def for_growth(cls, *also_scales: str) -> ArrheniusTemperature:
        """Arrhenius scaling of biomass growth (reads ``E_a_growth``).

        ``also_scales`` names extra Processes to scale by the *same* growth factor —
        used by wine to also scale the amino-acid assimilation swap, so its
        carbon/nitrogen refunds carry growth's temperature factor and cannot outrun
        the realised draw (decision D-32). Beer passes none (growth only).
        """
        return cls(
            name="arrhenius_growth",
            modifies=(GrowthNitrogenLimited.name, *also_scales),
            activation_energy="E_a_growth",
        )

    @classmethod
    def for_uptake(cls) -> ArrheniusTemperature:
        """Arrhenius scaling of the fermentative sugar-uptake flux (reads ``E_a_uptake``)."""
        return cls(
            name="arrhenius_uptake",
            modifies=SugarUptakeToEthanolCO2.name,
            activation_energy="E_a_uptake",
        )

    def factor(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> float:
        temp = float(y[schema.slice("T")][0])  # K, read from state (not params)
        # Always positive (exp), so no clamp: a single positive scalar on a
        # conserving vector cannot break a balance or flip a flux sign (D-11).
        return arrhenius_factor(temp, params[self._e_a_param], params["T_ref"])


#: Kelvin/Celsius offset (exact, SI-defined). Local to this module for the same
#: reason ``GAS_CONSTANT`` is (D-3/D-8's "code-with-citation" rule): the
#: :class:`ColemanQuadraticDeathTemperature` regression below is Coleman's own,
#: written directly in degrees C (their Table A2), so it needs this conversion —
#: nothing else in ``core`` does, by design (canonical internal temperature is
#: Kelvin; Celsius only appears at I/O edges per ``CLAUDE.md``, except here where
#: the *source regression itself* is in Celsius).
_KELVIN_OFFSET = 273.15


class ColemanQuadraticDeathTemperature(RateModifier):
    """Coleman's own quadratic temperature regression for ethanol-driven death.

    ``k'_d`` (:class:`~fermentation.core.kinetics.inactivation.EthanolInactivation`'s
    rate constant) is the one Coleman 2007 Table A2 parameter that is **quadratic**
    in temperature (°C), not log-linear like ``mu_max``/``beta_max`` — "the only
    QUADRATIC Coleman parameter and the steepest T-dependence" (wine_generic.yaml's
    ``k_prime_d`` provenance note). A single-``E_a`` Arrhenius tangent (the
    :class:`ArrheniusTemperature` form above) cannot reproduce a quadratic's
    curvature, so this modifier implements the regression directly rather than
    approximating it — decision **D-57**.

    ``ln(k'_d(T)) = a0 + a1·T_C + a2·T_C²`` (Coleman Table A2); the intercept
    ``a0`` cancels when normalised to ``T_ref`` (the temperature ``k_prime_d``
    itself is evaluated at, 20 C), leaving::

        factor(T) = exp( a1·(T_C - T_ref_C) + a2·(T_C² - T_ref_C²) )

    which is exactly 1 at ``T = T_ref`` (the measured ``k_prime_d`` used unscaled),
    matching the reference-anchored pattern :class:`ArrheniusTemperature` uses (D-11).

    **Why this was missing until D-57, and why it matters.** ``k_prime_d`` shipped
    with *no* temperature modifier at all (D-11/D-12): "M1 is isothermal at 20 C so
    no Arrhenius modifier is attached." That was correct scoping *for M1* — but M2
    added non-isothermal scenarios (temperature ramps, D-35/36) without anyone
    revisiting this, so every non-20 C wine/beer run since has driven growth and
    uptake with Arrhenius scaling while leaving death frozen at the 20 C rate. That
    asymmetry is inert on short/high-nitrogen runs (death is a minor contributor by
    the time dryness arrives) but compounds badly on long/nitrogen-limited runs at
    T != T_ref, exactly the regime D-56's Varela 2004 comparison (28 C) exposed:
    most of what D-56 read as a missing nitrogen-gated capacity-decline mechanism
    was this simpler, sourced T-scaling gap (see D-57 for the measured before/after).

    **Unphysical extrapolation guard.** The quadratic's vertex sits at ~11.3 C
    (``-a1/(2a2)``); below it the curve turns around and predicts *more* death as
    it gets *colder*, which is backwards. Coleman's own fitted range floors at
    11 C, so ``T_C`` is clamped to ``k_prime_d_t_floor`` before the quadratic is
    evaluated — a cold-cellar/lagering scenario gets the floor's (still gentle)
    death rate instead of a runaway one. No ceiling clamp: Coleman's fit runs to
    35 C and the quadratic keeps accelerating in the physically-correct direction
    above that (steeper death at higher T is exactly the "why heat causes stuck
    fermentations" mechanism this Process exists for).
    """

    name = "coleman_death_temperature"
    tier = Tier.PLAUSIBLE
    modifies = (EthanolInactivation.name,)
    #: ``T_ref`` is shared with :class:`ArrheniusTemperature`; the other three are
    #: this modifier's own (parameter-tier propagation caps ``X``/``X_dead`` by
    #: whichever tier is lowest, so wine's PLAUSIBLE/beer's SPECULATIVE `k_prime_d_a*`
    #: sourcing correctly floors each medium's output tier, D-1).
    reads = ("k_prime_d_a1", "k_prime_d_a2", "k_prime_d_t_floor", "T_ref")

    def factor(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> float:
        temp_c = float(y[schema.slice("T")][0]) - _KELVIN_OFFSET
        t_ref_c = params["T_ref"] - _KELVIN_OFFSET
        # Clamp only the low end (see class docstring); Coleman's fit has no
        # physically-motivated ceiling within fermentation-relevant temperatures.
        temp_c = max(temp_c, params["k_prime_d_t_floor"])
        a1 = params["k_prime_d_a1"]
        a2 = params["k_prime_d_a2"]
        # Always positive (exp): a single positive scalar on a conserving
        # X -> X_dead transfer cannot break the carbon/nitrogen balance (D-11).
        return float(math.exp(a1 * (temp_c - t_ref_c) + a2 * (temp_c**2 - t_ref_c**2)))
