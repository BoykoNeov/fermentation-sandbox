"""Kinetics-agnostic conservation and sanity checks over a trajectory.

The harness knows nothing about specific chemistry. A model supplies a
``quantity_fn(state) -> float`` that should be conserved (e.g. total carbon
across sugar, ethanol, CO2, and biomass), and these helpers assert it stays
constant along the trajectory. Encoding conservation as runtime/test assertions
is how we catch a model that quietly creates or destroys mass.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Protocol

import numpy as np

from fermentation.core.chemistry import (
    carbon_mass_fraction,
    nitrogen_mass_fraction,
    sugar_species,
)
from fermentation.core.state import FloatArray, StateSchema

QuantityFn = Callable[[FloatArray], float]


class TrajectoryLike(Protocol):
    """The subset of a trajectory these kinetics-agnostic checks read.

    Structural so both a plain :class:`~fermentation.runtime.integrate.Trajectory` and a
    :class:`~fermentation.runtime.schedule.ScheduledTrajectory` (and an ensemble member
    reconstructed as either) satisfy it — the harness never needs the scheduling extras,
    only the state grid.
    """

    @property
    def schema(self) -> StateSchema: ...
    @property
    def t(self) -> FloatArray: ...
    @property
    def y(self) -> FloatArray: ...
    def series(self, name: str) -> FloatArray: ...


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
    # Vicinal-diketone (VDK) pools (decision D-26): the diacetyl pathway routes carbon
    # sugar → α-acetolactate → diacetyl + CO2 → 2,3-butanediol, every step balanced on
    # this ledger. α-acetolactate is drawn from sugar (excretion) and booked at its own
    # C5 fraction; the spontaneous decarboxylation (C5 → C4 + CO2) is carbon-closing like
    # malic→lactic+CO2; the yeast reduction (diacetyl → butanediol) is a mole-for-mole
    # C4 → C4 transfer between two weighted pools (like esters → esters_gas). Weighting
    # all three at their own species' carbon fraction keeps total_carbon closed to machine
    # precision through the whole produce-then-reabsorb time course.
    if "acetolactate" in schema:
        w[schema.slice("acetolactate")] = carbon_mass_fraction("alpha_acetolactate")
    if "diacetyl" in schema:
        w[schema.slice("diacetyl")] = carbon_mass_fraction("diacetyl")
    if "butanediol" in schema:
        w[schema.slice("butanediol")] = carbon_mass_fraction("butanediol")
    # Acetaldehyde (decision D-27): a transient buffer on the main pathway. Production
    # borrows carbon from ethanol (E) and reduction returns it, both mole-for-mole C2 → C2,
    # so weighting acetaldehyde at its own carbon fraction keeps that borrow/return
    # carbon-neutral and total_carbon closed to machine precision through the produce-then-
    # reabsorb course. No sugar/CO2 term is involved — acetaldehyde de-lumps the uptake
    # Process's existing sugar→ethanol step rather than adding a parallel pathway.
    if "acetaldehyde" in schema:
        w[schema.slice("acetaldehyde")] = carbon_mass_fraction("acetaldehyde")
    # Excreted overflow pyruvate (decision D-49): PyruvateExcretion draws its carbon out of
    # sugar (C3, booked at pyruvate's fraction) and PyruvateReassimilation returns it to
    # ethanol + CO2 (C3 → C2 + C1, one mole each — carbon-closing like malic → lactic + CO2),
    # so the pool must be weighted here or the excrete-then-reassimilate course would read as
    # carbon created (during excretion) then destroyed (during re-assimilation). Weighted at
    # pyruvate's own carbon fraction — the SAME chemistry source the draw books against — so
    # the carbon into pyruvate equals the carbon out of S and total_carbon closes to machine
    # precision through the whole time course, leaving only the stranded residual as sugar
    # carbon parked as pyruvate. Re-assimilation returns to E/CO2 (not S) deliberately:
    # post-dryness S is 0, so a refund-to-sugar would be a no-op that destroys carbon. On a
    # keto-acid-pool-off run the pool is empty (constant 0 term). Nitrogen-free (a keto-acid),
    # so it is absent from total_nitrogen.
    if "pyruvate" in schema:
        w[schema.slice("pyruvate")] = carbon_mass_fraction("pyruvate")
    # Wine acid slots (decision D-18): the pH charge balance reads these, and the MLF
    # Process (decision D-23, when Oenococcus oeni is pitched) moves carbon malic (C4) ->
    # lactic (C3) + CO2 (C1) — balanced mole-for-mole — so they are weighted here for that
    # conversion to stay carbon-closing on this same ledger (no new conservation code). On
    # an undosed run the acids are inert (no active Process touches them, derivatives 0),
    # so this adds a constant term that drifts 0; ``cation_charge`` is a charge density,
    # not a carbon species, and stays weight 0 (``schema.zeros``). ``Byp`` weighting is
    # unchanged — the charge balance only *reads* it (include-by-reading), adding no carbon,
    # so ``total_carbon`` and the double-count are exactly as before.
    if "tartaric" in schema:
        w[schema.slice("tartaric")] = carbon_mass_fraction("tartaric_acid")
    if "malic" in schema:
        w[schema.slice("malic")] = carbon_mass_fraction("malic_acid")
    if "lactic" in schema:
        w[schema.slice("lactic")] = carbon_mass_fraction("lactic_acid")
    # Citrate (decision D-31): a dosed must input O. oeni co-metabolises during MLF. The
    # MalolacticCitrateMetabolism Process (when O. oeni is pitched) routes it as a lumped
    # citrate (C6) → α-acetolactate (C5) + CO2 (C1) conversion feeding the shared VDK
    # reservoir — carbon-closing mole-for-mole (6 = 5 + 1) on this same ledger, so weighting
    # citrate at its citric-acid fraction keeps total_carbon closed through that conversion.
    # On an undosed run citrate is inert (no active Process touches it, derivative 0), a
    # constant term that drifts 0. Charge-inactive (kept out of the D-18 pH balance in v1),
    # so it is a carbon term only — like malic/lactic, weighted here for the conversion.
    if "citrate" in schema:
        w[schema.slice("citrate")] = carbon_mass_fraction("citric_acid")
    # Amino-acid pool (decision D-32): a dosed, carbon- AND nitrogen-bearing wine pool. The
    # AminoAcidAssimilation swap debits it and refunds the displaced biomass carbon to sugar
    # and the displaced biomass nitrogen to ``N`` — a pure carbon/nitrogen-neutral transfer
    # (aa carbon → S, aa nitrogen → N), so weighting the pool at arginine's carbon fraction
    # (its representative species) keeps total_carbon closed through the swap. On an undosed
    # run the pool is empty and the Process is disabled (constant 0 term). Its NITROGEN side
    # is weighted in total_nitrogen below.
    if "amino_acids" in schema:
        w[schema.slice("amino_acids")] = carbon_mass_fraction("arginine")
    # Cell-wall debris pool (decision D-34): yeast autolysis (YeastAutolysis) turns dead biomass
    # X_dead into amino acids + this non-assimilable carbon-rich remainder (booked as glucan). The
    # dead-cell carbon r·f_C splits into the amino acids' carbon and this debris carbon, weighted
    # at glucan's carbon fraction to keep total_carbon closed through autolysis (the excess carbon
    # leaves the metabolite pools but stays on the ledger — the esters_gas idiom). Debris is
    # carbon-only (released nitrogen all goes to amino_acids), so it is absent from
    # total_nitrogen. On an undosed / autolysis-off run the pool is empty (constant 0 term).
    if "debris" in schema:
        w[schema.slice("debris")] = carbon_mass_fraction("glucan")
    # Mercaptan (thiol) pool (decision D-45): AutolyticMercaptan draws carbon from the amino-acid
    # pool into this pool (deaminating the nitrogen to N), so — unlike the carbon-free h2s — it
    # carries carbon and must be weighted or that draw would read as carbon destroyed. Booked at
    # methanethiol's carbon fraction (its representative species, the same one the Process draws
    # against), so total_carbon closes to machine precision through the transfer: the carbon into
    # mercaptans equals the carbon out of amino_acids. Nitrogen-free (the arginine N is deaminated
    # to the N pool), so it is absent from total_nitrogen. On an undosed / autolysis-off run the
    # pool is empty and the Process disabled (constant 0 term). NOTE: the add_copper verb (D-45)
    # removes mercaptans as precipitated copper mercaptide — carbon that legitimately LEAVES the
    # wine, booked as a negative external flow (the racking-debris precedent), so the run-wide
    # identity final == initial + Σ flows still holds even though total_carbon(state) drops.
    if "mercaptans" in schema:
        w[schema.slice("mercaptans")] = carbon_mass_fraction("methanethiol")
    # Brett volatile-phenol pools (decision D-40): the decarboxylase routes carbon hydroxycinnamics
    # (p-coumaric, C9) → vinylphenols (C8) + CO2 (C1), carbon-closing mole-for-mole like malic →
    # lactic + CO2; the reductase moves vinylphenols (C8) → ethylphenols (C8), a mole-for-mole C8 →
    # C8 transfer between two weighted pools like diacetyl → butanediol. Weighting all three at
    # their representative species (p-coumaric / 4-vinylphenol / 4-ethylphenol) keeps total_carbon
    # closed through the whole precursor → intermediate → product chain. On an undosed / un-pitched
    # run the pools are empty and the Processes disabled (constant 0 terms). X_brett/X_brett_dead
    # are weighted as biomass in the X block below (decision D-40 pt2, the X_mlf D-38 precedent).
    if "hydroxycinnamics" in schema:
        w[schema.slice("hydroxycinnamics")] = carbon_mass_fraction("p_coumaric_acid")
    if "vinylphenols" in schema:
        w[schema.slice("vinylphenols")] = carbon_mass_fraction("vinylphenol")
    if "ethylphenols" in schema:
        w[schema.slice("ethylphenols")] = carbon_mass_fraction("ethylphenol")
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
        # Bacterial biomass X_mlf (decision D-38): once MalolacticGrowth is active it
        # builds O. oeni biomass from the amino-acid pool, drawing the carbon shortfall from
        # sugar, so X_mlf carries carbon and must be weighted or that growth would read as
        # carbon destroyed. Booked at the SAME biomass_carbon_fraction the growth stoichiometry
        # draws against (bacterial ≈ yeast elemental composition, a documented v1 simplification),
        # so carbon closes exactly. This supersedes the v1 "X_mlf is carbon-free" scoping (D-23):
        # on a conversion-only run X_mlf only ever moves into X_mlf_dead (death D-39 / senescence
        # D-41 — both to an equally-weighted pool), so the added term stays closure-neutral; a
        # co-inoculation dose / pitch_mlf flow now carries this bacterial-biomass carbon.
        if "X_mlf" in schema:
            w[schema.slice("X_mlf")] = biomass_carbon_fraction
        # Non-viable bacterial biomass X_mlf_dead (decisions D-39/D-41): the SO₂ kill AND the benign
        # senescence baseline both move X_mlf into it, so — like the yeast X → X_dead transfer
        # (D-13) — it must be weighted at the SAME biomass_carbon_fraction or that loss would read
        # as carbon destroyed. (Terminal sink: autolysis reads only yeast X_dead, not X_mlf_dead.)
        if "X_mlf_dead" in schema:
            w[schema.slice("X_mlf_dead")] = biomass_carbon_fraction
        # Brett biomass X_brett (decision D-40 pt2): BrettGrowth builds it from the amino-acid pool,
        # drawing the carbon shortfall from ETHANOL, so X_brett carries carbon and must be weighted
        # or that growth would read as carbon destroyed. Booked at the SAME biomass_carbon_fraction
        # the growth stoichiometry draws against (Brett ≈ yeast elemental composition, a v1
        # simplification), so carbon closes exactly. Supersedes the pt1 "X_brett is carbon-free"
        # scoping: on a no-growth run X_brett is constant, so the added term is a constant offset
        # that drifts 0; a co-inoculation dose / pitch_brett flow now carries this biomass carbon.
        if "X_brett" in schema:
            w[schema.slice("X_brett")] = biomass_carbon_fraction
        # Non-viable Brett biomass X_brett_dead (decision D-40 pt3): BrettDeath moves X_brett into
        # it under SO₂, so — like the X → X_dead / X_mlf → X_mlf_dead transfers — it is weighted at
        # the same fraction so that death reads as carbon-neutral.
        if "X_brett_dead" in schema:
            w[schema.slice("X_brett_dead")] = biomass_carbon_fraction
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
    # Amino-acid pool (decision D-32): the aa pool carries nitrogen (arginine, 4 N per
    # molecule) and the AminoAcidAssimilation swap moves that nitrogen aa → N (ammonium
    # refund) mole-for-mole, so weighting the pool at arginine's nitrogen fraction keeps
    # total_nitrogen closed through the swap. This is the first nitrogen-bearing tracked
    # species besides biomass; on an undosed run the pool is empty (constant 0 term).
    if "amino_acids" in schema:
        w[schema.slice("amino_acids")] = nitrogen_mass_fraction("arginine")
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
        # Bacterial biomass X_mlf (decision D-38): MalolacticGrowth builds it by
        # assimilating the amino-acid pool's nitrogen (arginine), so X_mlf carries nitrogen and
        # is weighted at the SAME biomass_nitrogen_fraction the growth draws against — the pool
        # loses exactly the nitrogen X_mlf gains, so total_nitrogen closes. On a conversion-only run
        # X_mlf only moves into X_mlf_dead (death/senescence), so closure holds; a co-inoculation
        # dose / pitch_mlf flow now carries it.
        if "X_mlf" in schema:
            w[schema.slice("X_mlf")] = biomass_nitrogen_fraction
        # Non-viable bacterial biomass X_mlf_dead (decisions D-39/D-41): retains its nitrogen, at
        # the same biomass_nitrogen_fraction so BOTH the X_mlf → X_mlf_dead death AND senescence
        # transfers are nitrogen-neutral (the yeast X → X_dead precedent, D-13).
        if "X_mlf_dead" in schema:
            w[schema.slice("X_mlf_dead")] = biomass_nitrogen_fraction
        # Brett biomass X_brett (decision D-40 pt2): BrettGrowth assimilates the amino-acid pool's
        # nitrogen (arginine) into biomass, so X_brett carries nitrogen, weighted at the same
        # biomass_nitrogen_fraction the growth draws against — the pool loses exactly the nitrogen
        # X_brett gains, so total_nitrogen closes. Constant (⇒ 0 drift) on a no-growth run.
        if "X_brett" in schema:
            w[schema.slice("X_brett")] = biomass_nitrogen_fraction
        # Non-viable Brett biomass X_brett_dead (decision D-40 pt3): retains its nitrogen, counted
        # at the same fraction so the X_brett → X_brett_dead death transfer is nitrogen-neutral.
        if "X_brett_dead" in schema:
            w[schema.slice("X_brett_dead")] = biomass_nitrogen_fraction
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


def _evaluate(traj: TrajectoryLike, quantity_fn: QuantityFn) -> FloatArray:
    return np.array([quantity_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])


def max_drift(traj: TrajectoryLike, quantity_fn: QuantityFn) -> float:
    """Maximum absolute deviation of ``quantity_fn`` from its initial value."""
    q = _evaluate(traj, quantity_fn)
    return float(np.max(np.abs(q - q[0]))) if q.size else 0.0


def assert_conserved(
    traj: TrajectoryLike,
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


def assert_nonnegative(
    traj: TrajectoryLike, variables: tuple[str, ...], *, atol: float = 1e-9
) -> None:
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
