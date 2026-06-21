"""Stoichiometric constants for the validated core.

Molar masses and carbon-atom counts of the species the core tracks. These are
*exact consequences of the chemical formulae* (and the standard atomic masses),
not empirical kinetic parameters — so, like the conversion factors in
:mod:`fermentation.units` (decision D-3), they live in code with citations rather
than in the provenance-backed parameter store.

This module is the **single source of truth** for fermentation stoichiometry: the
carbon/mass conservation checks (:mod:`fermentation.validation.conservation`) and
the sugar-uptake Process both derive their numbers here, so a conservation check
can never silently disagree with the kinetics it audits.

Empirical, strain-dependent quantities do **not** belong here. In particular the
*biomass elemental composition* (carbon/nitrogen content of dry yeast) is uncertain
and shared with the growth Process, so it is a :class:`~fermentation.parameters.\
schema.Parameter` (provenance store), passed into the conservation builders — see
decision D-8.

Standard atomic masses (IUPAC 2021, g/mol): C 12.011, H 1.008, O 15.999.
"""

from __future__ import annotations

from fermentation.core.state import StateSchema

# -- standard atomic masses (IUPAC 2021), g/mol -------------------------------
_M_C = 12.011
_M_H = 1.008
_M_O = 15.999

# -- molar masses of tracked species, g/mol (derived from their formulae) -----
#: Glucose / fructose, C6H12O6 — the lumped wine hexose and beer's first sugar.
M_GLUCOSE = 6 * _M_C + 12 * _M_H + 6 * _M_O
#: Maltose, C12H22O11.
M_MALTOSE = 12 * _M_C + 22 * _M_H + 11 * _M_O
#: Maltotriose, C18H32O16.
M_MALTOTRIOSE = 18 * _M_C + 32 * _M_H + 16 * _M_O
#: Ethanol, C2H6O.
M_ETHANOL = 2 * _M_C + 6 * _M_H + 1 * _M_O
#: Carbon dioxide, CO2.
M_CO2 = 1 * _M_C + 2 * _M_O
#: Water, H2O (hydrolysis bookkeeping for di-/trisaccharide uptake).
M_WATER = 2 * _M_H + 1 * _M_O

#: Molar mass [g/mol] keyed by species name. ``fermentation.core.media`` sugar
#: component names ("glucose", "maltose", "maltotriose") are keys here.
MOLAR_MASS: dict[str, float] = {
    "glucose": M_GLUCOSE,
    "fructose": M_GLUCOSE,
    "maltose": M_MALTOSE,
    "maltotriose": M_MALTOTRIOSE,
    "ethanol": M_ETHANOL,
    "CO2": M_CO2,
}

#: Carbon atoms per molecule, keyed by species name.
CARBON_ATOMS: dict[str, int] = {
    "glucose": 6,
    "fructose": 6,
    "maltose": 12,
    "maltotriose": 18,
    "ethanol": 2,
    "CO2": 1,
}


def carbon_mass_fraction(species: str) -> float:
    """Grams of carbon per gram of ``species`` (exact from its formula).

    Used to weight each state variable when summing total carbon. Raises
    ``KeyError`` for an unknown species so a typo fails loudly rather than
    silently dropping a carbon-bearing term from a conservation check.
    """
    try:
        return CARBON_ATOMS[species] * _M_C / MOLAR_MASS[species]
    except KeyError:
        raise KeyError(f"unknown species {species!r}; known: {sorted(MOLAR_MASS)}") from None


def sugar_species(schema: StateSchema) -> list[str]:
    """Map a schema's ``S`` slots to chemical species names, in slot order.

    Beer's ``S`` names its components (glucose/maltose/maltotriose); wine's single
    lumped slot is treated as a hexose (glucose). A multi-slot ``S`` without
    component names is an error — its carbon/mass weights are undefined without
    knowing which sugars occupy the slots.

    This is the single source of truth shared by every sugar-aware consumer —
    the carbon-conservation check (:mod:`fermentation.validation.conservation`)
    and the kinetic Processes that draw carbon from sugar — so a check can never
    disagree with the kinetics it audits (decision D-8). It lives here, in the
    core, because the validation layer may not be imported by the core that needs
    it (the one-directional dependency rule).
    """
    spec = schema.spec("S")
    if spec.components:
        return list(spec.components)
    if spec.size == 1:
        return ["glucose"]
    raise ValueError(
        f"sugar 'S' has {spec.size} slots but no component names; cannot assign carbon fractions"
    )


# -- Gay-Lussac stoichiometry of one hexose -----------------------------------
# Anaerobic alcoholic fermentation of a hexose:
#
#     C6H12O6 -> 2 C2H5OH + 2 CO2
#
# mass-balanced by atom count (2*M_ETHANOL + 2*M_CO2 == M_GLUCOSE to the third
# decimal). These are the *theoretical* maximum yields; the realised ethanol
# yield is a few percent lower because cells divert carbon to glycerol, organic
# acids and biomass (the realised-yield Parameter, and the Tier-2 glycerol
# Process — see decision D-8). The split is exposed here so the sugar-uptake
# Process and the carbon-conservation check use one definition.

#: Grams of ethanol produced per gram of hexose consumed (theoretical, ~0.511).
ETHANOL_PER_HEXOSE = 2 * M_ETHANOL / M_GLUCOSE
#: Grams of CO2 evolved per gram of hexose consumed (theoretical, ~0.489).
CO2_PER_HEXOSE = 2 * M_CO2 / M_GLUCOSE

#: Hexose units released per molecule on complete hydrolysis. Glucose/fructose are
#: hexoses already; maltose -> 2, maltotriose -> 3. Each released hexose ferments
#: by Gay-Lussac (-> 2 ethanol + 2 CO2), so a sugar's per-gram ethanol/CO2 yield
#: scales with its hexose count. The di-/trisaccharide mass gain comes from
#: hydrolysis water pulled from the solvent (why beer's S+E+CO2 mass does not
#: close — see ``validation.total_mass`` and decision D-8).
HEXOSE_UNITS: dict[str, int] = {
    "glucose": 1,
    "fructose": 1,
    "maltose": 2,
    "maltotriose": 3,
}


def ethanol_yield(species: str) -> float:
    """Grams of ethanol per gram of ``species`` fermented (theoretical Gay-Lussac).

    Generalises :data:`ETHANOL_PER_HEXOSE` to the di-/trisaccharides: a sugar with
    ``n`` hexose units yields ``2n`` ethanol per molecule. M1 uses this theoretical
    split (not the realised ``Y_ethanol_sugar`` parameter) so carbon and mass close
    exactly — the realised yield, net of glycerol/biomass diversion, is a Tier-2
    concern (decision D-8). Raises ``KeyError`` for an unknown/unfermentable
    species so a typo fails loudly.
    """
    try:
        return 2 * HEXOSE_UNITS[species] * M_ETHANOL / MOLAR_MASS[species]
    except KeyError:
        raise KeyError(
            f"no fermentation yield for {species!r}; known: {sorted(HEXOSE_UNITS)}"
        ) from None


def co2_yield(species: str) -> float:
    """Grams of CO2 per gram of ``species`` fermented (theoretical Gay-Lussac).

    Companion to :func:`ethanol_yield`; for a hexose equals :data:`CO2_PER_HEXOSE`.
    ``ethanol_yield(s) + co2_yield(s)`` exceeds 1 for di-/trisaccharides by exactly
    the hydrolysis water taken up — mass closes only when that water is tracked,
    which M1 does not (decision D-8), so beer relies on the carbon balance.
    """
    try:
        return 2 * HEXOSE_UNITS[species] * M_CO2 / MOLAR_MASS[species]
    except KeyError:
        raise KeyError(
            f"no fermentation yield for {species!r}; known: {sorted(HEXOSE_UNITS)}"
        ) from None
