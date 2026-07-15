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

It also hosts :data:`ESTER_SPECS`, the **canonical ester registry** (decision D-96) — the
single source of truth every layer derives from (state slots, synthesis, stripping, the
carbon ledger, the OAV aroma set), so a fourth ester is *one entry*, not a new code path.
"""

from __future__ import annotations

from dataclasses import dataclass

from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.state import FloatArray, StateSchema


@dataclass(frozen=True)
class EsterSpec:
    """One aroma-active ester: its liquid pool, headspace twin, molecule and rate constant.

    ``pool`` is the liquid state variable, ``gas_pool`` its CO₂-stripped headspace
    bookkeeping twin (decision D-20), ``species`` the molecule **both** are carbon-weighted
    by in :func:`~fermentation.validation.conservation.total_carbon` — weighting a
    liquid/gas pair at the *same* fraction is what makes stripping carbon-neutral — and
    ``k_param`` the per-medium synthesis rate constant.

    **Each pool IS its molecule (decision D-96).** Before D-96 a single lumped ``esters``
    pool was carbon-weighted as *ethyl acetate* but read on the OAV lens against *isoamyl
    acetate*'s threshold — a split identity that made the fruity OAV non-physical (~761 for a
    wine, implying ~23 mg/L isoamyl acetate against a real ceiling of ~1–3 mg/L). Splitting
    the lump into single-molecule pools, each weighted **and** perceived as itself, removes
    that seam rather than papering over it: no ester pool carries a ``lumped`` composition
    assumption any more.
    """

    pool: str
    gas_pool: str
    species: str
    k_param: str
    #: Human-readable note for the liquid pool's :class:`~fermentation.core.state.VarSpec`
    #: description — kept here so the schema text cannot drift from the registry.
    note: str


#: The THREE aroma-active esters the sim tracks (decision D-96) — each a **single-molecule**
#: pool, so none carries the fixed-lump-composition caveat (that now survives only on
#: ``fusels``/``mercaptans``, where it is true). Two families:
#:
#: * the **acetate esters** (ATF1/alcohol acetyltransferase) — ``ethyl_acetate``, the bulk
#:   solventy one (tens of mg/L), and ``isoamyl_acetate``, the trace potent banana one
#:   (~1–3 mg/L). Same enzyme ⇒ they share ``E_a_esters``, the per-medium sourced
#:   temperature shape (decision D-21).
#: * the **ethyl ester of a medium-chain fatty acid** — ``ethyl_hexanoate`` (apple/pineapple),
#:   the highest-OAV ester in wine. **Documented v1 simplification (D-96):** it is
#:   EEB1/EHT1-derived, a *different* enzyme from the acetates, but v1 shares
#:   ``E_a_esters`` because no separate activation energy is sourced. A per-ester ``E_a`` is
#:   the named deferred refinement.
#:
#: Each ``k`` is sourced **independently** to its own molecule's measured concentration range
#: (D-96, the load-bearing call): the lump's composition is *derived* from three
#: independently-anchored rates, never a single fitted ratio splitting one ``k_ester``. Adding
#: a fourth ester (ethyl octanoate, phenylethyl acetate) is one entry here plus its params.
ESTER_SPECS: tuple[EsterSpec, ...] = (
    EsterSpec(
        "ethyl_acetate",
        "ethyl_acetate_gas",
        "ethyl_acetate",
        "k_ethyl_acetate",
        note="ethyl acetate (C4H8O2) — the bulk solventy/nail-polish acetate ester (ATF1)",
    ),
    EsterSpec(
        "isoamyl_acetate",
        "isoamyl_acetate_gas",
        "isoamyl_acetate",
        "k_isoamyl_acetate",
        note="isoamyl acetate (C7H14O2) — the trace, potent BANANA acetate ester (ATF1); "
        "the only ester the D-69 aging hydrolysis fades",
    ),
    EsterSpec(
        "ethyl_hexanoate",
        "ethyl_hexanoate_gas",
        "ethyl_hexanoate",
        "k_ethyl_hexanoate",
        note="ethyl hexanoate (C8H16O2) — the APPLE/PINEAPPLE ethyl ester of a "
        "medium-chain fatty acid (EEB1/EHT1)",
    ),
)

#: The ester the D-69 aging hydrolysis acts on — the **banana** acetate. Its fade to
#: ``fusels`` + ``Byp`` is the whole point of that Process, and at D-96 the 5:2 carbon split
#: became **exact**: the debited molecule finally *is* isoamyl acetate (C7 → isoamyl alcohol
#: C5 + acetic acid C2), where D-69 had to debit ethyl acetate and split as if it were
#: isoamyl acetate. See :class:`~fermentation.core.kinetics.aging.EsterHydrolysis`.
HYDROLYSING_ESTER: EsterSpec = ESTER_SPECS[1]


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


def refund_carbon_to_sugar(
    d: FloatArray, y: FloatArray, schema: StateSchema, carbon: float
) -> None:
    """Add ``carbon`` [g C/L/h] back to ``S``, split across slots by carbon content.

    The exact inverse of :func:`draw_carbon_from_sugar`: where a byproduct Process *drew*
    its carbon out of sugar, another Process may *refund* carbon to sugar (e.g. when a
    fraction of that carbon is instead sourced from the amino-acid pool — the D-33 fusel
    Ehrlich re-route — or when amino-acid-funded biomass spares sugar for ethanol — the
    D-32 swap). Each slot ``i`` receives carbon in proportion to the carbon it currently
    holds, so ``Σ_i (d[S_i]·c_i) = +carbon`` exactly and the refund matches whatever draw
    it undoes to machine precision. Uses ``+=`` so it composes with any draw the same
    Process makes. Slots are clamped ≥ 0 (mirroring the draw); with no sugar carbon present
    (``Σ s_j c_j ≤ 0``) nothing is refunded — callers must guarantee sugar is present (the
    fermentative flux they track is zero at ``S = 0``, so this edge is never hit in
    practice). Serves wine's single slot and beer's three identically.
    """
    s_slice = schema.slice("S")
    species = sugar_species(schema)
    s = [max(float(y[s_slice.start + i]), 0.0) for i in range(len(species))]
    carbon_total = sum(s[i] * carbon_mass_fraction(sp) for i, sp in enumerate(species))
    if carbon_total <= 0.0:
        return
    for i in range(len(species)):
        if s[i] > 0.0:
            d[s_slice.start + i] += carbon * s[i] / carbon_total


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
