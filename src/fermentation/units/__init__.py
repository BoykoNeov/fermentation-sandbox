"""Unit conversions — applied **only at I/O boundaries**.

Canonical internal units (a single representation; conversions live here and are
called only when crossing the engine's edges):

    quantity         internal unit        rationale
    ---------------  -------------------  -------------------------------------
    concentration    g/L  (== kg/m3)      fermentation-literature convention,
                                          numerically identical to SI kg/m3
    temperature      K (kelvin)           Arrhenius needs absolute temperature
    time             h (hours)            dominant unit for kinetic constants
                                          (mu_max /h); benchmarks quoted in days
    volume           L (litre)
    mass             g (gram)

Industry units (degrees Brix, specific gravity, degrees Plato, %ABV, degrees C,
days) appear *only* on the far side of these functions. The hot integration loop
sees nothing but the canonical floats above.
"""

from fermentation.units.convert import (
    abv_from_ethanol,
    apparent_gravity,
    brix_to_sg,
    brix_to_sugar_gpl,
    celsius_to_kelvin,
    days_to_hours,
    hours_to_days,
    kelvin_to_celsius,
    plato_to_sg,
    real_to_apparent_extract,
    sg_to_brix,
    sg_to_plato,
    sugar_gpl_to_brix,
)

__all__ = [
    "abv_from_ethanol",
    "apparent_gravity",
    "brix_to_sg",
    "brix_to_sugar_gpl",
    "celsius_to_kelvin",
    "days_to_hours",
    "hours_to_days",
    "kelvin_to_celsius",
    "plato_to_sg",
    "real_to_apparent_extract",
    "sg_to_brix",
    "sg_to_plato",
    "sugar_gpl_to_brix",
]
