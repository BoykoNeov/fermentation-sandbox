"""Nitrogen-limited yeast growth — the first validated-core kinetic Process.

Models the textbook stuck/sluggish mechanism: cells divide while
yeast-assimilable nitrogen (YAN) lasts, then stop dividing once it runs out
(Monod shutoff on ``N``), while the separate sugar-uptake Process keeps
fermenting at a declining rate. Growth is also Monod-limited on total sugar so it
cannot run on an empty sugar pool, giving the lag→exponential→stationary biomass
shape the milestone targets.

Conservation (decision D-8) is built into the stoichiometry, not bolted on, so
this Process conserves both rigorous atom balances *on its own*:

  * **Nitrogen.** New biomass assimilates YAN, so ``dN/dt = -f_N · dX/dt`` and
    ``N + f_N · X`` (``total_nitrogen``) is conserved to solver tolerance.
  * **Carbon.** The biomass carbon skeleton is drawn from sugar, so ``dS`` removes
    exactly ``f_C · dX/dt`` grams of carbon and ``carbon(S) + f_C · X``
    (``total_carbon``) is conserved to solver tolerance.

Both ``f_N`` (``biomass_N_fraction``) and ``f_C`` (``biomass_C_fraction``) are the
same provenance-backed Parameters the conservation checks read, which is why the
balances close exactly rather than approximately.

**M1 simplification — no anabolic CO2.** Every gram of sugar carbon this Process
removes goes into biomass; respiratory/anabolic CO2 is not modelled. That makes
the biomass yield carbon-cheap (~0.82 g biomass / g sugar in isolation), which is
unrealistic on its own but immaterial for M1: nitrogen caps biomass near
``X0 + N0/f_N`` (~2-3 g/L for wine), so only ~1-2 % of sugar is ever diverted
here, and no M1 benchmark probes biomass yield. (Revisit when tightening the
``co2_peak_then_tail`` ±5 % carbon-ratio benchmark — the missing anabolic CO2
eats into that budget.) ``total_mass`` over ``{S, E, CO2}`` therefore does **not**
close under growth (sugar leaves without matching E/CO2); carbon is the invariant
to check, exactly as decision D-8 scopes it.

Tier: **plausible** — a sound, literature-standard mechanism (Monod kinetics,
nitrogen limitation), not yet validated against the benchmark curves. Promote to
validated once it reproduces the §2.2 lag→exp→stationary shape.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


class GrowthNitrogenLimited(Process):
    """Monod biomass growth, co-limited by sugar and yeast-assimilable nitrogen.

    ``mu = mu_max · S_total/(K_s + S_total) · N/(K_n + N)`` and ``dX/dt = mu·X``.
    Nitrogen and the biomass carbon skeleton are drawn from ``N`` and ``S`` so the
    run conserves both atom balances (see the module docstring and decision D-8).
    """

    name = "growth_nitrogen_limited"
    tier = Tier.PLAUSIBLE
    touches = ("X", "S", "N")
    #: Parameters this Process reads. Declared so the parameter-tier-propagation
    #: task can cap the Process's output tier by its inputs' tiers (D-1); not yet
    #: consumed by ``ProcessSet.tier_of`` — see ``docs/plans/milestone-1-tasks.md``.
    reads: tuple[str, ...] = ("mu_max", "K_s", "K_n", "biomass_N_fraction", "biomass_C_fraction")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        x = float(y[schema.slice("X")][0])
        n = max(float(y[schema.slice("N")][0]), 0.0)
        s_slice = schema.slice("S")
        s_block = y[s_slice]
        s_total = max(float(s_block.sum()), 0.0)
        # Monod shutoff: no biomass, no sugar, or no nitrogen -> no growth. The
        # guards also keep the rate non-negative against small solver excursions.
        if x <= 0.0 or s_total <= 0.0 or n <= 0.0:
            return d

        mu = (
            params["mu_max"]
            * (s_total / (params["K_s"] + s_total))
            * (n / (params["K_n"] + n))
        )
        dx = mu * x  # biomass growth rate [g/L/h]
        d[schema.slice("X")] = dx
        d[schema.slice("N")] = -params["biomass_N_fraction"] * dx  # YAN into biomass

        # Carbon skeleton for the new biomass is drawn from sugar. Distribute the
        # carbon demand across the available sugar slots in proportion to their
        # current mass, converting each to a sugar-mass rate by *that slot's own*
        # carbon fraction. Because the mass weights sum to 1, the carbon removed
        # equals the demand exactly for any split, so carbon closes for wine (one
        # slot) and beer (three) alike (decision D-8).
        carbon_demand = params["biomass_C_fraction"] * dx  # [g C / L / h]
        for offset, species in enumerate(sugar_species(schema)):
            s_i = float(s_block[offset])
            if s_i <= 0.0:
                continue
            carbon_i = (s_i / s_total) * carbon_demand
            d[s_slice.start + offset] = -carbon_i / carbon_mass_fraction(species)
        return d
