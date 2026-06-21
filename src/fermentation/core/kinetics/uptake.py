"""Sugar uptake to ethanol + CO2 — the fermentative flux of the validated core.

The second M1 kinetic Process, and the one that does the bulk work: it converts
sugar to ethanol and CO2 by the theoretical Gay-Lussac stoichiometry. The defining
choice (design decision recorded with the milestone) is that uptake is **catalysed
by biomass, not coupled to growth**::

    r = q_sugar_max · X · S/(K_sugar_uptake + S)

so the rate depends on how much yeast is present, not on whether it is dividing.
Fermentation therefore continues at a declining rate after yeast-assimilable
nitrogen runs out and :class:`~fermentation.core.kinetics.growth.GrowthNitrogen\
Limited` has shut off — the textbook "non-growing cells finish the ferment"
mechanism, and what carries a wine must from ~24 Brix to dryness. (A
growth-coupled uptake would stall at high residual sugar the moment nitrogen
ran out.)

**Theoretical Gay-Lussac yields (decision D-8).** Ethanol and CO2 yields are not
parameters: they come from :func:`~fermentation.core.chemistry.ethanol_yield` /
:func:`~fermentation.core.chemistry.co2_yield`, the same stoichiometry the
conservation checks use. M1 uses the *theoretical* split (0.511 / 0.489 g per g
hexose) rather than the realised ``Y_ethanol_sugar`` (~0.47), so carbon and mass
conserve to machine precision. The realised yield — the few percent of carbon
cells divert to glycerol and biomass — is a Tier-2 concern; ``Y_ethanol_sugar``
is its hook and is deliberately unread in M1.

**Beer's sequential uptake (smooth catabolite repression).** Wort sugars are
consumed in preference order (glucose, then maltose, then maltotriose). Each
higher sugar is suppressed while a more-preferred one remains, via a smooth
repression factor ``K_repression / (K_repression + S_preferred)`` per
more-preferred slot. Smooth rather than a hard switch so the derivative stays
continuous for the BDF solver. This relies on the schema's ``S`` slot order being
the preference order, which :func:`~fermentation.core.media.beer_schema` defines.
Wine has a single sugar slot, so repression never applies.

**Conservation.** Every gram of sugar carbon removed reappears as ethanol or CO2
carbon (``ethanol_yield``/``co2_yield`` are carbon-balanced per hexose unit), so
``total_carbon`` is conserved for wine and beer alike. For a hexose the two yields
sum to exactly 1, so ``total_mass`` over ``{S, E, CO2}`` is also conserved for
wine; beer's di-/trisaccharides take up hydrolysis water, so mass does not close
there and beer relies on the carbon balance (decision D-8).

Tier: **plausible** — a sound, literature-standard mechanism (saturating
biomass-catalysed uptake, catabolite repression), not yet validated against the
§2.2 benchmark curves.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import co2_yield, ethanol_yield, sugar_species
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


class SugarUptakeToEthanolCO2(Process):
    """Biomass-catalysed sugar uptake fermented to ethanol + CO2 (Gay-Lussac).

    Per sugar slot ``i`` (in schema/preference order)::

        r_i = q_sugar_max · X · S_i/(K_sugar_uptake + S_i) · repression_i
        repression_i = Π_{j<i} K_repression / (K_repression + S_j)

    with ``dS_i = -r_i``, ``dE = Σ ethanol_yield(species_i)·r_i`` and
    ``dCO2 = Σ co2_yield(species_i)·r_i``. See the module docstring for why this is
    decoupled from growth and conserves carbon (and mass, for wine).
    """

    name = "sugar_uptake_to_ethanol_co2"
    tier = Tier.PLAUSIBLE
    touches = ("S", "E", "CO2")
    #: Parameters this Process reads (for parameter-tier propagation, D-1).
    #: ``K_repression`` is read only on the multi-sugar (beer) path, but is
    #: declared unconditionally — under D-1 it will conservatively cap wine's
    #: E/CO2 tier even though wine never represses. See milestone-1-tasks.md.
    reads: tuple[str, ...] = ("q_sugar_max", "K_sugar_uptake", "K_repression")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        x = float(y[schema.slice("X")][0])
        if x <= 0.0:  # no catalyst, no fermentation
            return d

        s_slice = schema.slice("S")
        s_block = y[s_slice]
        species = sugar_species(schema)
        # Clamp each slot to >= 0 *before* it enters a Monod term or a repression
        # denominator: a small negative solver excursion would otherwise flip the
        # uptake sign (creating sugar, driving E/CO2 negative). Mirrors the guards
        # in GrowthNitrogenLimited.
        s = [max(float(s_block[i]), 0.0) for i in range(len(species))]

        q_max = params["q_sugar_max"]
        k_su = params["K_sugar_uptake"]
        k_rep = params["K_repression"]

        e_slice = schema.slice("E")
        co2_slice = schema.slice("CO2")

        repression = 1.0  # the most-preferred sugar (slot 0) is never repressed
        for i, sp in enumerate(species):
            if i > 0:
                # Suppress this sugar while the next-more-preferred one remains.
                repression *= k_rep / (k_rep + s[i - 1])
            s_i = s[i]
            if s_i <= 0.0:
                continue
            r_i = q_max * x * (s_i / (k_su + s_i)) * repression  # g sugar_i /L/h
            d[s_slice.start + i] = -r_i
            d[e_slice] += ethanol_yield(sp) * r_i
            d[co2_slice] += co2_yield(sp) * r_i
        return d
