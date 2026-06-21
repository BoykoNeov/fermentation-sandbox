"""Ethanol inhibition of the fermentative flux — a multiplicative rate modifier.

Ethanol is the yeast's own product and its own poison: as it accumulates it slows
sugar uptake and, past a strain-specific tolerance, collapses viability. Unlike
growth and uptake this mechanism does not *add* a flux — it *scales* an existing
one — so it cannot be a summed :class:`~fermentation.core.process.Process` in the
additive :class:`~fermentation.core.process.ProcessSet`. It is a
:class:`~fermentation.core.process.RateModifier` instead: it returns a factor that
``ProcessSet`` multiplies onto the *whole* contribution of the sugar-uptake
Process. (Decision **D-10** — building the modifier hook here, one task ahead of
the Arrhenius temperature modifier that reuses it.)

**Functional form (Levenspiel / Luong "toxic power").** With ethanol ``E`` and the
strain tolerance ``E_max`` (``ethanol_tolerance``)::

    f(E) = (1 - E/E_max)**n   for 0 <= E < E_max,   else 0

``E_max`` is a *wall*: the flux reaches zero there, matching the parameter's
"viability collapses past tolerance" semantics. The exponent ``n`` (>1) keeps the
touchdown smooth — ``f'(E_max) = 0``, so the rate eases to zero with no derivative
kink for the BDF solver to choke on (a raw ``n=1`` linear form has that kink). At
``E=0`` the factor is 1 (no inhibition); it decreases monotonically to 0.

**Conservation is automatic.** The factor scales uptake's entire ``(dS, dE, dCO2)``
contribution by one scalar, so the Gay-Lussac carbon/mass balance uptake already
respects is preserved — a slower carbon-neutral flux is still carbon-neutral.

**Only uptake is targeted in M1.** Ethanol also inhibits growth, but growth is shut
off by nitrogen limitation long before ethanol climbs high (see
:mod:`~fermentation.core.kinetics.growth`), so inhibiting it would change nothing
the benchmarks probe. The Arrhenius modifier (next task) will target both rates.

Tier: **plausible** — a standard, literature-backed inhibition form (Levenspiel
1980; Luong 1985), not yet validated against the §2.2 benchmark curves. With the
*placeholder* ``ethanol_tolerance`` (110 g/L, below a 24 °Brix must's ~124-135 g/L
final ethanol) an inhibited wine run stalls short of dryness; that is a
parameter-sourcing/tuning concern (a high-alcohol must implies a high-tolerance
strain), not a flaw in the form — see D-10 and milestone-1-tasks.md.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.kinetics.uptake import SugarUptakeToEthanolCO2
from fermentation.core.process import RateModifier
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


class EthanolInhibition(RateModifier):
    """Levenspiel/Luong ethanol inhibition scaling the fermentative uptake flux.

    ``factor = (1 - E/E_max)**n`` clamped to ``[0, 1]`` (``E_max =
    ethanol_tolerance``, ``n = ethanol_inhibition_exponent``). Multiplied onto
    :class:`~fermentation.core.kinetics.uptake.SugarUptakeToEthanolCO2`'s whole
    contribution by :class:`~fermentation.core.process.ProcessSet`.
    """

    name = "ethanol_inhibition"
    tier = Tier.PLAUSIBLE
    #: Reference the uptake Process by its ``name`` (rename-safe) rather than a
    #: bare string literal.
    modifies = (SugarUptakeToEthanolCO2.name,)
    reads: tuple[str, ...] = ("ethanol_tolerance", "ethanol_inhibition_exponent")

    def factor(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> float:
        # Clamp E >= 0 before it enters the wall term: a negative solver excursion
        # would otherwise read as *less* inhibition (factor > 1), speeding uptake.
        e = max(float(y[schema.slice("E")][0]), 0.0)
        remaining = 1.0 - e / params["ethanol_tolerance"]
        if remaining <= 0.0:  # at/above tolerance the flux is fully shut down
            return 0.0
        return float(remaining ** params["ethanol_inhibition_exponent"])
