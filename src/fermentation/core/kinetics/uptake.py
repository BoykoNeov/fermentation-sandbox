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

**Gay-Lussac yields with a realised-yield correction (decisions D-8, D-16).**
Ethanol and CO2 yields come from :func:`~fermentation.core.chemistry.ethanol_yield`
/ :func:`~fermentation.core.chemistry.co2_yield`, the same *theoretical* split
(0.511 / 0.489 g per g hexose) the conservation checks use. Real ferments divert a
few percent of sugar carbon to glycerol, organic acids and higher alcohols, so the
realised ethanol yield is ~0.47-0.49. We apply that correction *here* rather than
swapping in a lower ``Y_ethanol_sugar`` constant: the theoretical split is scaled
by ``(1 - f)`` and the diverted carbon ``f`` is deposited into the ``Gly``
(glycerol) and ``Byp`` (minor byproducts, booked as succinic acid) state pools.
Because the deposited carbon exactly equals the carbon removed from ethanol+CO2,
``total_carbon`` still closes to machine precision *with the byproducts tracked*
(decision D-16). The diversion is set by two parameters, ``Y_glycerol_sugar`` and
``Y_byproduct_sugar`` (g byproduct / g sugar consumed); **both default to 0**, and
at 0 this Process is exactly the theoretical Gay-Lussac core (so it stays togglable
off and the validated-core mass balance is intact). Wine sources them from
glycerol/byproduct data; beer carries 0 (its byproducts are out of M1 scope), so
its CO2 benchmark is untouched. Glycerol/byproducts are more reduced than the
ethanol route, drawing redox H/O from the solvent — so wine ``{S,E,CO2}`` mass no
longer closes once the diversion is on (carbon, which counts ``Gly``/``Byp``, is
the invariant then — exactly as for biomass).

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
biomass-catalysed uptake, catabolite repression). The §2.2 benchmarks now pass,
which confirms the plausible tier is earned; it is deliberately **not** promoted to
validated, which is reserved for validation against independent *measured*
time-series (handoff D-C: none exist yet; decision D-17). Output tiers stay capped
at speculative regardless via the parameters this reads (``K_repression`` and
``Y_byproduct_sugar`` are speculative — D-1 propagation), so the cap is on the
inputs, not the form.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import (
    carbon_mass_fraction,
    co2_yield,
    ethanol_yield,
    sugar_species,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Representative species whose formula carbon-accounts each byproduct pool. The
#: minor-byproduct lump (``Byp``) is booked as succinic acid (decision D-16).
_GLYCEROL = "glycerol"
_BYPRODUCT = "succinic_acid"


class SugarUptakeToEthanolCO2(Process):
    """Biomass-catalysed sugar uptake fermented to ethanol + CO2 (Gay-Lussac).

    Per sugar slot ``i`` (in schema/preference order)::

        r_i = q_sugar_max · X · S_i/(K_sugar_uptake + S_i) · repression_i
        repression_i = Π_{j<i} K_repression / (K_repression + S_j)

    with ``dS_i = -r_i``. The carbon of ``r_i`` splits into ethanol/CO2 and the
    realised-yield byproduct pools (decision D-16)::

        f_C   = Y_glycerol_sugar·c(glycerol) + Y_byproduct_sugar·c(succinic)
        scale = 1 - f_C / c(species_i)            # share still going to ethanol+CO2
        dE   = Σ ethanol_yield(species_i)·scale·r_i
        dCO2 = Σ co2_yield(species_i)·scale·r_i
        dGly = Σ Y_glycerol_sugar·r_i
        dByp = Σ Y_byproduct_sugar·r_i

    where ``c(·)`` is the species' carbon mass fraction. The carbon deposited in
    ``Gly``/``Byp`` exactly equals the carbon scaled out of ethanol+CO2, so
    ``total_carbon`` closes for any yields. With both yields 0, ``scale = 1`` and
    this reduces to the theoretical Gay-Lussac core. See the module docstring for
    why uptake is decoupled from growth and how this preserves conservation.
    """

    name = "sugar_uptake_to_ethanol_co2"
    tier = Tier.PLAUSIBLE
    touches = ("S", "E", "CO2", "Gly", "Byp")
    #: Parameters this Process reads; :meth:`ProcessSet.tier_of` folds their tiers
    #: into the output tier of ``S``/``E``/``CO2``/``Gly``/``Byp`` (parameter-tier
    #: propagation, D-1). ``K_repression`` is read only on the multi-sugar (beer)
    #: path, and ``Y_glycerol_sugar``/``Y_byproduct_sugar`` are 0 for beer, but all
    #: are declared unconditionally so they conservatively cap the output tier.
    reads: tuple[str, ...] = (
        "q_sugar_max",
        "K_sugar_uptake",
        "K_repression",
        "Y_glycerol_sugar",
        "Y_byproduct_sugar",
    )

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
        # Realised-yield byproduct diversion (decision D-16). Default 0 ⇒ the
        # theoretical Gay-Lussac core (togglable off). ``.get`` keeps hand-built
        # test param maps that predate these knobs working as the pure core.
        y_gly = params.get("Y_glycerol_sugar", 0.0)
        y_byp = params.get("Y_byproduct_sugar", 0.0)
        # Carbon [g C] diverted to byproducts per g sugar consumed — independent of
        # which sugar, since it is booked against each pool's own carbon fraction.
        diverted_c = y_gly * carbon_mass_fraction(_GLYCEROL) + y_byp * carbon_mass_fraction(
            _BYPRODUCT
        )

        e_slice = schema.slice("E")
        co2_slice = schema.slice("CO2")
        divert = diverted_c > 0.0
        if divert:
            if "Gly" not in schema or "Byp" not in schema:
                raise ValueError(
                    "byproduct diversion is on (Y_glycerol_sugar/Y_byproduct_sugar > 0) "
                    "but the schema has no 'Gly'/'Byp' pool to receive the carbon"
                )
            gly_slice = schema.slice("Gly")
            byp_slice = schema.slice("Byp")

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
            if divert:
                # Share of this sugar's carbon still fermented to ethanol+CO2.
                scale = 1.0 - diverted_c / carbon_mass_fraction(sp)
                if scale < 0.0:
                    raise ValueError(
                        f"byproduct yields divert more carbon ({diverted_c:.4g} g C/g) than "
                        f"{sp} carries ({carbon_mass_fraction(sp):.4g} g C/g); reduce "
                        "Y_glycerol_sugar/Y_byproduct_sugar"
                    )
                d[gly_slice] += y_gly * r_i
                d[byp_slice] += y_byp * r_i
            else:
                scale = 1.0
            d[e_slice] += ethanol_yield(sp) * scale * r_i
            d[co2_slice] += co2_yield(sp) * scale * r_i
        return d
