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
nitrogen limitation). The §2.2 wine benchmark now passes and this Process
reproduces Coleman 2007 line-for-line (``tests/test_coleman_reconstruction.py``),
which *confirms* the plausible tier is earned — but it is deliberately **not**
promoted to validated. VALIDATED is reserved for validation against independent
*measured* time-series (handoff D-C: none exist yet); the §2.2 wine window was
re-anchored to Coleman, the same source these constants come from, so clearing it
is a faithful-implementation cross-check, not independent validation. Promote when
real curves drop into the data-ready harness (decision D-17).
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The intracellular assimilable-nitrogen store (decision D-250) — wine schemas only. ``N``
#: keeps its D-209 meaning, EXTRACELLULAR assimilable nitrogen, and is the only one of the two
#: the acid-base charge balance reads; this one titrates nothing because it is inside a cell.
STORED_NITROGEN_KEY = "stored_nitrogen"


def assimilable_nitrogen_pools(y: FloatArray, schema: StateSchema) -> tuple[float, float]:
    """``(extracellular N, intracellular store)`` [g N/L], each clamped >= 0 (decision D-250).

    Yeast growth is limited by, and draws on, *both*: nitrogen the cell has already transported
    in is exactly as available to anabolism as nitrogen still in the must. The store is absent
    from beer's schema (the amino-acid ledger is wine-only, D-30/D-32), which reads as 0.0 and
    makes every caller reduce to its pre-D-250 form term-for-term.
    """
    n = max(float(y[schema.slice("N")][0]), 0.0)
    if STORED_NITROGEN_KEY not in schema:
        return n, 0.0
    return n, max(float(y[schema.slice(STORED_NITROGEN_KEY)][0]), 0.0)


def add_assimilable_nitrogen(
    d: FloatArray, y: FloatArray, schema: StateSchema, flux: float
) -> None:
    """Book ``flux`` [g N/L/h] across ``N`` and the store **in proportion to what each holds**.

    The single source of truth for how anabolic nitrogen moves between the two pools
    (decision D-250, the D-8/D-32 idiom): growth's draw and the
    :class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation` refund of that same
    draw must use the *identical* split, or the swap refunds to a pool growth did not debit and
    ``N`` drifts upward — which is the D-248 charge defect returning through a second door.

    **Proportional rather than store-first**, deliberately: a preferential draw puts a C0 kink
    exactly where the store empties, and the BDF Jacobian probe straddles precisely that kind of
    gate ([[feedback-a-gate-is-a-discontinuity-the-solver-probes]]). Proportional splitting is
    smooth, keeps the two pools' *ratio* invariant under drawdown — which is what lets the SUM
    ``N + store`` follow exactly the trajectory ``N`` alone followed before D-250 — and cannot
    drive either pool negative faster than it drives the total.

    On an empty total the whole flux is booked to ``N``; both callers guard on a positive growth
    rate, which requires a positive total, so that branch is unreachable and exists only so the
    helper is total-function.
    """
    n, stored = assimilable_nitrogen_pools(y, schema)
    total = n + stored
    if total <= 0.0 or stored <= 0.0:
        d[schema.slice("N")] += flux
        return
    d[schema.slice("N")] += flux * (n / total)
    d[schema.slice(STORED_NITROGEN_KEY)] += flux * (stored / total)


def biomass_growth_rate(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Base (pre-modifier) biomass growth rate ``dX/dt = mu·X`` [g/L/h].

    ``mu = mu_max · S_total/(K_s + S_total) · N/(K_n + N)``. The Monod shutoff guards
    (no biomass, no sugar, or no nitrogen ⇒ 0, clamped ≥ 0 against solver excursions)
    are the same ones :class:`GrowthNitrogenLimited` applies. Factored out (decision
    D-32) as the single source of the growth rate so the amino-acid assimilation swap
    (:class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation`) anchors to
    the *identical* rate the growth Process builds biomass at — a divergence here would
    let the swap refund more sugar carbon than growth drew. This returns the **base**
    rate; the Arrhenius/carrying-capacity :class:`~fermentation.core.process.RateModifier`
    scaling is applied by :class:`~fermentation.core.process.ProcessSet` to growth *and*
    the swap alike, so both track the realised rate (decision D-32).
    """
    x = float(y[schema.slice("X")][0])
    # Both nitrogen pools limit growth (D-250): what is still in the must and what the cells
    # have already transported in. Pre-D-250 this read `N` alone, which was the same number
    # because uptake's surplus was booked there -- see `assimilable_nitrogen_pools`.
    n = sum(assimilable_nitrogen_pools(y, schema))
    s_total = max(float(y[schema.slice("S")].sum()), 0.0)
    if x <= 0.0 or s_total <= 0.0 or n <= 0.0:
        return 0.0
    mu = params["mu_max"] * (s_total / (params["K_s"] + s_total)) * (n / (params["K_n"] + n))
    return mu * x


class GrowthNitrogenLimited(Process):
    """Monod biomass growth, co-limited by sugar and yeast-assimilable nitrogen.

    ``mu = mu_max · S_total/(K_s + S_total) · N/(K_n + N)`` and ``dX/dt = mu·X``.
    Nitrogen and the biomass carbon skeleton are drawn from ``N`` and ``S`` so the
    run conserves both atom balances (see the module docstring and decision D-8).
    """

    name = "growth_nitrogen_limited"
    tier = Tier.PLAUSIBLE
    touches = ("X", "S", "N")
    #: ``stored_nitrogen`` is medium-conditional (D-250): growth draws its nitrogen from the
    #: extracellular pool and the intracellular store together, split by what each holds, but
    #: the store is wine-only (it mirrors the wine-only amino-acid ledger, D-30/D-32). Growth is
    #: the one primary-fermentation Process wired into BOTH media, so it is the case
    #: :attr:`Process.touches_where_present` exists for.
    touches_where_present = (STORED_NITROGEN_KEY,)
    #: Parameters this Process reads. :meth:`ProcessSet.tier_of` folds their tiers
    #: into the output tier of ``X``/``S``/``N`` (parameter-tier propagation, D-1),
    #: so a speculative parameter here caps those outputs at speculative.
    reads: tuple[str, ...] = ("mu_max", "K_s", "K_n", "biomass_N_fraction", "biomass_C_fraction")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        s_slice = schema.slice("S")
        s_block = y[s_slice]
        s_total = max(float(s_block.sum()), 0.0)
        # Monod shutoff (no biomass/sugar/nitrogen -> no growth) lives in the shared
        # rate helper, which also clamps the rate non-negative against solver excursions.
        dx = biomass_growth_rate(y, schema, params)  # biomass growth rate [g/L/h]
        if dx <= 0.0:
            return d
        d[schema.slice("X")] = dx
        # YAN into biomass, split across the must pool and the intracellular store in
        # proportion to what each holds (D-250). The shared helper is what keeps this draw and
        # the D-32 swap's refund of it on the same split.
        add_assimilable_nitrogen(d, y, schema, -params["biomass_N_fraction"] * dx)

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
