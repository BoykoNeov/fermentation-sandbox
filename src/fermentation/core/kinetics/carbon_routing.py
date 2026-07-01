"""Shared carbon-routing helpers for the metabolite-byproduct Processes.

Two small building blocks reused by every Process that produces a carbon-accounted
byproduct pool from the fermentative flux:

* :func:`fermentative_flux_shape` — the biomass-catalysed sugar Monod activity the
  byproduct Processes couple to, so their production tracks the *same* flux they are
  metabolically downstream of.
* :func:`draw_carbon_from_sugar` — routes a byproduct's carbon *out of* ``S`` (option a1,
  decision D-19) so the pool it fills is real carbon-accounted state and ``total_carbon``
  closes to machine precision.

Extracted from :mod:`fermentation.core.kinetics.byproducts` (decision D-26) so the ester/
fusel aroma Processes and the vicinal-diketone (diacetyl) Processes share one definition —
the same single-source-of-truth discipline the chemistry constants follow. The behaviour
is unchanged from the byproducts beat (D-19); only the home moved.
"""

from __future__ import annotations

from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.state import FloatArray, StateSchema


def draw_carbon_from_sugar(
    d: FloatArray, y: FloatArray, schema: StateSchema, carbon: float
) -> None:
    """Subtract ``carbon`` [g C/L/h] from ``S``, split across slots by carbon content.

    Routes a byproduct's carbon out of sugar (option a1, decision D-19) so the pool it
    fills becomes carbon-accounted state. Each slot ``i`` gives up carbon in proportion
    to the carbon it currently holds, ``s_i·c_i / Σ_j s_j·c_j``; converting that back to
    grams of sugar, ``d[S_i] -= carbon · s_i / Σ_j s_j·c_j``. By construction
    ``Σ_i (d[S_i]·c_i) = -carbon`` exactly, so the carbon removed from sugar equals the
    carbon the caller deposits in its pool and ``total_carbon`` closes to machine
    precision. Slots are clamped ≥ 0 (mirroring the flux/uptake guards) so a solver
    undershoot cannot flip a draw into sugar *creation*; with no sugar carbon present
    (``Σ s_j c_j ≤ 0``) nothing is drawn. This serves wine's single slot and beer's
    three identically (different carbon fractions per slot are handled exactly).
    """
    s_slice = schema.slice("S")
    species = sugar_species(schema)
    s = [max(float(y[s_slice.start + i]), 0.0) for i in range(len(species))]
    carbon_total = sum(s[i] * carbon_mass_fraction(sp) for i, sp in enumerate(species))
    if carbon_total <= 0.0:
        return
    for i in range(len(species)):
        if s[i] > 0.0:
            d[s_slice.start + i] -= carbon * s[i] / carbon_total


def fermentative_flux_shape(y: FloatArray, schema: StateSchema, k_sat: float) -> float:
    """Biomass-catalysed sugar Monod term ``X · S_total/(K + S_total)`` [g/L].

    The dimensionless-but-for-``X`` activity proxy the fermentative uptake Process
    runs on (``q_sugar_max·X·S/(K+S)``), reused by the byproduct Processes so their
    production tracks the *same* flux they are metabolically coupled to — which is what
    makes the run-integrated "total scales as f_byproduct/f_flux" cancellation clean and
    predictable (see the byproducts module docstring). Sugar is summed across slots
    (1 for wine, 3 for beer) and clamped ≥ 0 against solver undershoot, mirroring the
    guards in the uptake/growth Processes.
    """
    x = max(float(y[schema.slice("X")][0]), 0.0)
    s_total = max(float(y[schema.slice("S")].sum()), 0.0)
    if x <= 0.0 or s_total <= 0.0:
        return 0.0
    return x * (s_total / (k_sat + s_total))
