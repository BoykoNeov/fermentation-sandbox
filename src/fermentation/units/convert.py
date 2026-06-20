"""Boundary unit conversions.

Formulas are standard enology/brewing approximations; each is cited inline. They
are convenience conversions, not kinetic parameters, so they live in code rather
than the provenance store — but the citations make their assumptions auditable.
"""

from __future__ import annotations

# Density of pure ethanol at 20 C, g/mL. (CRC Handbook.)
_ETHANOL_DENSITY_GPML = 0.78924

# Density of pure water at 20 C, g/mL.
_WATER_DENSITY_GPML = 0.99820


def celsius_to_kelvin(celsius: float) -> float:
    return celsius + 273.15


def kelvin_to_celsius(kelvin: float) -> float:
    return kelvin - 273.15


def days_to_hours(days: float) -> float:
    return days * 24.0


def hours_to_days(hours: float) -> float:
    return hours / 24.0


def brix_to_sg(brix: float) -> float:
    """Degrees Brix -> specific gravity (20/20 C).

    Inverse of the standard cubic Brix(SG) polynomial; expressed here as the
    widely used closed form (e.g. Wikipedia "Brix", attributed to the ASBC
    extract tables):

        SG = 1 + brix / (258.6 - (brix / 258.2) * 227.1)
    """
    return 1.0 + brix / (258.6 - (brix / 258.2) * 227.1)


def sg_to_brix(sg: float) -> float:
    """Specific gravity -> degrees Brix (ASBC cubic, accurate 0-40 Brix).

    Brix = -668.962 + 1262.45*SG - 776.43*SG^2 + 182.94*SG^3
    """
    return -668.962 + 1262.45 * sg - 776.43 * sg**2 + 182.94 * sg**3


def sg_to_plato(sg: float) -> float:
    """Specific gravity -> degrees Plato (ASBC cubic).

        P = -616.868 + 1111.14*SG - 630.272*SG^2 + 135.997*SG^3

    Degrees Plato (brewing) and degrees Brix (enology) both measure % sucrose by
    mass and are numerically near-identical; the small difference comes from the
    fitting polynomial used by each industry.
    """
    return -616.868 + 1111.14 * sg - 630.272 * sg**2 + 135.997 * sg**3


def plato_to_sg(plato: float) -> float:
    """Degrees Plato -> specific gravity.

    SG = 1 + plato / (258.6 - (plato / 258.2) * 227.1)
    """
    return 1.0 + plato / (258.6 - (plato / 258.2) * 227.1)


def brix_to_sugar_gpl(brix: float, sg: float | None = None) -> float:
    """Degrees Brix -> dissolved sugar concentration in g/L.

    Brix is grams sucrose per 100 g solution, so the volumetric concentration is

        sugar [g/L] = brix [g/100 g] * density [g/mL] * 10

    If ``sg`` is not supplied it is derived from ``brix`` via :func:`brix_to_sg`
    (density of the solution ~= SG * water density, and SG is dimensionless
    relative to water at 20/20 C so density [g/mL] ~= SG * 0.9982).
    """
    if sg is None:
        sg = brix_to_sg(brix)
    density_gpml = sg * _WATER_DENSITY_GPML
    return brix * density_gpml * 10.0


def sugar_gpl_to_brix(sugar_gpl: float, sg: float = 1.0) -> float:
    """Dissolved sugar [g/L] -> degrees Brix at a known/assumed ``sg``.

    Inverse of :func:`brix_to_sugar_gpl` for a *given* solution density. Because
    Brix and density are mutually dependent, callers that need high accuracy
    should pass the measured ``sg``; the default ``sg=1.0`` is a dilute-solution
    approximation.
    """
    density_gpml = sg * _WATER_DENSITY_GPML
    return sugar_gpl / (density_gpml * 10.0)


def abv_from_ethanol(ethanol_gpl: float) -> float:
    """Ethanol concentration [g/L] -> alcohol by volume [% v/v].

        ABV [%] = (ethanol [g/L] / ethanol_density [g/mL]) / 10

    i.e. volume of ethanol per volume of solution, in percent. Uses pure-ethanol
    density at 20 C; this ignores volume contraction on mixing and so is a close
    approximation rather than an exact figure.
    """
    ethanol_ml_per_l = ethanol_gpl / _ETHANOL_DENSITY_GPML
    return ethanol_ml_per_l / 10.0
