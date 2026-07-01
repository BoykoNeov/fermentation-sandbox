"""Biomass carrying-capacity cap on growth — the opt-in residual-nitrogen floor.

**The gap this closes (decision D-30, §3.2).** In the validated core,
:class:`~fermentation.core.kinetics.growth.GrowthNitrogenLimited` is the *sole*
nitrogen sink, and its only shutoff is the Monod term ``N/(K_n + N)`` with a tiny
``K_n`` (~0.009 g/L). Nothing else caps biomass, so a wine ferment builds
``X ≈ X0 + N0/f_N`` and strips yeast-assimilable nitrogen (YAN) to ~0 by day ~1.3
*regardless of dose* — 80 and 300 mg N/L both bottom out at zero. That erases every
downstream low-N signal: the D-29 H₂S ``K_h2s_n/(K_h2s_n + N)`` inverse gate reads
``N→0`` for every must, so the cross-must "high-YAN suppresses sulfide" lever is
muted to a few percent, and a residual-YAN readout (for a future MLF-with-growth
model) is impossible.

**The mechanism.** Real yeast populations saturate below the nitrogen ceiling —
oxygen/sterol limitation, cell-density quorum effects, membrane crowding — leaving
assimilable nitrogen unconsumed at a dose-dependent floor. The textbook lumped form
is a logistic carrying capacity: growth slows as biomass ``X`` approaches a cap
``K`` (``biomass_carrying_capacity``) and stops at it. Because this *scales* an
existing flux rather than adding one, it is a
:class:`~fermentation.core.process.RateModifier`, not a summed Process — it returns

    factor(X) = clamp(1 - X/K, 0, 1)

which :class:`~fermentation.core.process.ProcessSet` multiplies onto the *whole*
contribution of the growth Process. Linear ``1 - X/K`` (not the smoothed
``(1-·)**n`` of the ethanol wall) is deliberate: unlike ethanol, ``X`` self-limits
— growth → 0 as ``X → K``, so the state never gets driven past the wall and there
is no derivative kink for the BDF solver to catch on. The ``[0, 1]`` clamp still
guards a solver overshoot ``X > K`` from flipping the factor negative (which would
turn growth into a biomass/nitrogen *source*).

**Conservation is automatic.** Growth removes nitrogen and sugar-carbon in fixed
stoichiometric proportion to ``dX`` (``dN = -f_N·dX``, carbon skeleton drawn from
``S``). Scaling that whole contribution by one scalar preserves every proportion,
so ``total_nitrogen`` and ``total_carbon`` still close to solver tolerance with the
cap on — the nitrogen simply stays in the ``N`` pool once growth saturates, instead
of being forced into biomass. This is the crux that makes a carrying-capacity cap
the right vehicle: less biomass, exact balances, residual N left behind.

**Departure from Coleman — why this is opt-in, not default (decision D-30).**
Coleman, Fish & Block (2007), the keystone wine model, has *no* biomass cap: it
consumes all YAN and builds full nitrogen-proportional biomass at every dose, and
``tests/test_coleman_reconstruction.py`` confirms our core reproduces that
line-for-line at 80 *and* 330 mg N/L. A residual-N floor necessarily departs from
that curve (leaving assimilable N ⇒ less biomass ⇒ a slower uptake tail). The two
are genuinely incompatible: restoring the H₂S lever *requires* residual assimilable
N that differs by dose, which *means* not matching Coleman's zero-residual biomass.
So this modifier ships **isolable and disabled by default** (prime directive #3):
it is wired into the wine medium but the compile seam *disables* it unless a
scenario opts in via ``carrying_capacity_gpl``. Disabled ⇒ it contributes a factor
of 1 **and** is excluded from tier derivation (``ProcessSet`` counts enabled, not
nonzero, modifiers), so an undosed wine run is byte-for-byte the validated core and
keeps growth at PLAUSIBLE — the Coleman reconstruction, §2.2 dryness/ABV, and the
fusel/ester benchmarks are untouched. Opt in and growth's ``X``/``S``/``N`` outputs
drop to SPECULATIVE, honestly flagging the residual-N departure.

Tier: **speculative** — a standard functional form (logistic carrying capacity) but
the cap value is an author estimate, and the mechanism is a deliberate departure
from the validated Coleman anchor. The cross-must H₂S lever it restores and the
residual YAN it leaves are both speculative-on-speculative.

SCOPE (v1): wine-only (the H₂S lever and the prospective MLF-with-growth model are
wine concerns), mirroring the wine-only MLF wiring; beer carrying capacity is a
deferred follow-up. The MLF *unblock* is **prospective**: MLF v1 is conversion-only
with pH/ethanol/molecular-SO₂/cardinal-T gates and *no* nitrogen gate, so residual N
does not change current MLF behaviour — it enables a future MLF-with-growth model.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.kinetics.growth import GrowthNitrogenLimited
from fermentation.core.process import RateModifier
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


class BiomassCarryingCapacity(RateModifier):
    """Logistic carrying-capacity cap scaling the nitrogen-limited growth flux.

    ``factor = clamp(1 - X/K, 0, 1)`` with ``K = biomass_carrying_capacity``.
    Multiplied onto :class:`~fermentation.core.kinetics.growth.GrowthNitrogenLimited`'s
    whole contribution by :class:`~fermentation.core.process.ProcessSet`, so nitrogen
    and carbon balances are preserved and growth eases to zero as ``X → K``, leaving a
    dose-dependent residual of yeast-assimilable nitrogen behind (decision D-30).
    """

    name = "biomass_carrying_capacity"
    tier = Tier.SPECULATIVE
    #: Reference the growth Process by its ``name`` (rename-safe) rather than a bare
    #: string literal.
    modifies = (GrowthNitrogenLimited.name,)
    reads: tuple[str, ...] = ("biomass_carrying_capacity",)

    def factor(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> float:
        # Clamp X >= 0 before the wall term so a negative solver excursion cannot read
        # as factor > 1 (which would *speed* growth past the uncapped rate).
        x = max(float(y[schema.slice("X")][0]), 0.0)
        remaining = 1.0 - x / params["biomass_carrying_capacity"]
        if remaining <= 0.0:  # at/above the cap, growth is fully shut down
            return 0.0
        return float(remaining)
