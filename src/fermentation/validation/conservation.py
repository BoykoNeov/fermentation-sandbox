"""Kinetics-agnostic conservation and sanity checks over a trajectory.

The harness knows nothing about specific chemistry. A model supplies a
``quantity_fn(state) -> float`` that should be conserved (e.g. total carbon
across sugar, ethanol, CO2, and biomass), and these helpers assert it stays
constant along the trajectory. Encoding conservation as runtime/test assertions
is how we catch a model that quietly creates or destroys mass.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.state import FloatArray, StateSchema
from fermentation.runtime.integrate import Trajectory

QuantityFn = Callable[[FloatArray], float]


# -- conserved-quantity builders for the real chemistry -----------------------
#
# Each builder closes over a schema and returns a ``QuantityFn`` for the harness
# above. Stoichiometric weights come from ``fermentation.core.chemistry`` (a single
# source of truth shared with the kinetics), so the check cannot disagree with the
# model it audits. What each balance covers, and why, is decision D-8.


def _weighted_sum(weights: FloatArray) -> QuantityFn:
    """A ``QuantityFn`` returning ``weights @ y`` (weights precomputed once)."""

    def quantity(y: FloatArray) -> float:
        return float(weights @ y)

    return quantity


def total_carbon(
    schema: StateSchema, *, biomass_carbon_fraction: float | None = None
) -> QuantityFn:
    """Total carbon [g C / L] over sugars, ethanol, CO2 and biomass.

    Carbon is the rigorous conservation invariant (atoms cannot be created or
    destroyed): the model must route every gram of sugar carbon into ethanol,
    CO2 or biomass. Each ``S`` slot is weighted by its species' carbon fraction
    (single-slot ``S`` = hexose); ``E`` and ``CO2`` by theirs.

    Biomass ``X`` contributes ``biomass_carbon_fraction * X``. That fraction is a
    Parameter the growth Process also reads (decision D-8), so it is *passed in*
    here to keep the check consistent with the kinetics. If the schema has an
    ``X`` variable and no fraction is given, this raises rather than silently
    under-counting carbon (which would report a false conservation violation).
    """
    w = schema.zeros()
    s_slice = schema.slice("S")
    for offset, species in enumerate(sugar_species(schema)):
        w[s_slice.start + offset] = carbon_mass_fraction(species)
    if "E" in schema:
        w[schema.slice("E")] = carbon_mass_fraction("ethanol")
    if "CO2" in schema:
        w[schema.slice("CO2")] = carbon_mass_fraction("CO2")
    # Realised-yield byproduct sinks (decision D-16): the sugar-uptake Process
    # routes the carbon it diverts from ethanol/CO2 into these pools, so they must
    # be weighted or carbon would read as destroyed when glycerol is active. The
    # minor-byproduct lump is carbon-accounted as succinic acid (its representative
    # species). Both share the one stoichiometry source with the kinetics.
    if "Gly" in schema:
        w[schema.slice("Gly")] = carbon_mass_fraction("glycerol")
    if "Byp" in schema:
        w[schema.slice("Byp")] = carbon_mass_fraction("succinic_acid")
    # Aroma byproduct pools (decision D-19): the ester/fusel Processes route their
    # carbon out of sugar, so these pools must be weighted or carbon would read as
    # destroyed when they accumulate. Esters book as ethyl acetate, fusels as isoamyl
    # alcohol — the same representative species the Processes draw against, from the
    # one chemistry source of truth, so the check matches the kinetics. No overlap
    # with Byp: under D-19 Byp is organic-acids/polyols only (higher alcohols moved
    # to the fusels pool), so the former Byp double-count is gone.
    if "esters" in schema:
        w[schema.slice("esters")] = carbon_mass_fraction("ethyl_acetate")
    if "fusels" in schema:
        w[schema.slice("fusels")] = carbon_mass_fraction("isoamyl_alcohol")
    # Volatilized esters (decision D-20): the EsterVolatilization sink moves carbon
    # from the liquid ``esters`` pool into this headspace pool as CO2 strips it, so it
    # must be weighted at the *same* ethyl-acetate fraction or that liquid→gas transfer
    # would read as carbon destroyed. This mirrors how evolved ``CO2`` stays counted —
    # the carbon leaves the liquid but not the ledger, so closure holds to machine
    # precision while wine's liquid esters honestly fall with temperature.
    if "esters_gas" in schema:
        w[schema.slice("esters_gas")] = carbon_mass_fraction("ethyl_acetate")
    if "X" in schema:
        if biomass_carbon_fraction is None:
            raise ValueError(
                "schema has biomass 'X'; pass biomass_carbon_fraction (the value "
                "the growth Process uses) so the carbon check matches the kinetics"
            )
        w[schema.slice("X")] = biomass_carbon_fraction
        # Inactivated cells are still biomass of the same composition: ethanol
        # inactivation moves mass X -> X_dead, so counting both pools at the same
        # carbon fraction keeps that transfer carbon-neutral (decision D-13).
        if "X_dead" in schema:
            w[schema.slice("X_dead")] = biomass_carbon_fraction
    return _weighted_sum(w)


def total_nitrogen(
    schema: StateSchema, *, biomass_nitrogen_fraction: float | None = None
) -> QuantityFn:
    """Total nitrogen [g N / L]: free YAN plus nitrogen bound in biomass.

    As cells grow they assimilate ``N`` (yeast-assimilable nitrogen) into biomass,
    so the invariant is ``N + biomass_nitrogen_fraction * X``. Like the biomass
    carbon fraction, ``biomass_nitrogen_fraction`` is a Parameter the growth
    Process reads and is passed in (decision D-8); required iff the schema has
    biomass. Conserved once the nitrogen-limited growth Process exists.
    """
    w = schema.zeros()
    if "N" in schema:
        w[schema.slice("N")] = 1.0
    if "X" in schema:
        if biomass_nitrogen_fraction is None:
            raise ValueError(
                "schema has biomass 'X'; pass biomass_nitrogen_fraction (the value "
                "the growth Process uses) so the nitrogen check matches the kinetics"
            )
        w[schema.slice("X")] = biomass_nitrogen_fraction
        # Inactivated biomass retains its nitrogen: count X_dead so the X -> X_dead
        # inactivation transfer stays nitrogen-neutral (decision D-13).
        if "X_dead" in schema:
            w[schema.slice("X_dead")] = biomass_nitrogen_fraction
    return _weighted_sum(w)


def total_mass(schema: StateSchema) -> QuantityFn:
    """Total mass [g/L] of the abiotic conversion species ``{S, E, CO2}`` — wine only.

    Mass closes only for a *single hexose* (decision D-8): the wine reaction
    ``C6H12O6 -> 2 C2H5OH + 2 CO2`` is mass-balanced (180.156 = 92.138 + 88.018
    g/mol), so ``S + E + CO2`` is conserved to solver tolerance. It does **not**
    generalise, for the same untracked-solvent reason biomass is excluded:
    di-/trisaccharides *hydrolyse*, pulling water from the solvent into the product
    pool (maltose adds ~5.3% mass, maltotriose ~7.1%), so ``S + E + CO2`` is not a
    beer invariant. Carbon carries no such term (water has no carbon) and is the
    rigorous cross-medium invariant — so this builder **rejects a multi-component
    sugar** and beer relies on :func:`total_carbon`.

    It stays scoped to ``{S, E, CO2}`` for one further reason (decision D-16): the
    realised-yield byproduct sinks ``Gly``/``Byp`` are more reduced than the
    sugar→ethanol+CO2 route and draw redox H/O from the solvent like biomass, so
    ``{S, E, CO2}`` mass does not close once that diversion is active. This is thus
    the **validated-core** (byproduct-free) mass check — assert it only on a
    glycerol-off configuration; ``total_carbon`` (which weights ``Gly``/``Byp``)
    is the invariant when byproducts are on.
    """
    if "S" in schema and schema.spec("S").size > 1:
        raise ValueError(
            "total_mass is a hexose/wine check: a multi-component sugar hydrolyses "
            "(pulling solvent water into S+E+CO2), so mass does not close. Use "
            "total_carbon, which is conserved across media (decision D-8)."
        )
    w = schema.zeros()
    for name in ("S", "E", "CO2"):
        if name in schema:
            w[schema.slice(name)] = 1.0
    return _weighted_sum(w)


def _evaluate(traj: Trajectory, quantity_fn: QuantityFn) -> FloatArray:
    return np.array([quantity_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])


def max_drift(traj: Trajectory, quantity_fn: QuantityFn) -> float:
    """Maximum absolute deviation of ``quantity_fn`` from its initial value."""
    q = _evaluate(traj, quantity_fn)
    return float(np.max(np.abs(q - q[0]))) if q.size else 0.0


def assert_conserved(
    traj: Trajectory,
    quantity_fn: QuantityFn,
    *,
    rtol: float = 1e-6,
    atol: float = 1e-9,
    label: str = "quantity",
) -> None:
    """Assert ``quantity_fn`` stays within tolerance of its initial value.

    Tolerance is relative to the initial magnitude plus an absolute floor, so it
    behaves sensibly when the conserved quantity is near zero. Raises
    ``AssertionError`` (so it reads naturally inside tests and runtime checks).
    """
    q = _evaluate(traj, quantity_fn)
    if q.size == 0:
        return
    tol = atol + rtol * abs(q[0])
    drift = np.max(np.abs(q - q[0]))
    if drift > tol:
        worst = int(np.argmax(np.abs(q - q[0])))
        raise AssertionError(
            f"{label} not conserved: drift {drift:.3e} > tol {tol:.3e} "
            f"(initial {q[0]:.6g}, at t={traj.t[worst]:.3g}h -> {q[worst]:.6g})"
        )


def assert_nonnegative(traj: Trajectory, variables: tuple[str, ...], *, atol: float = 1e-9) -> None:
    """Assert the named variables never go meaningfully negative.

    Concentrations and biomass are physical and must stay >= 0; a small negative
    excursion within ``atol`` is tolerated as solver noise.
    """
    for name in variables:
        series = traj.series(name)
        worst = float(np.min(series))
        if worst < -atol:
            idx = int(np.argmin(series))
            raise AssertionError(f"{name} went negative: {worst:.3e} at t={traj.t[idx]:.3g}h")
