"""Temperature-/metabolism-driven aroma byproducts — esters and fusel alcohols.

The first Milestone-2 (Tier-2) beat (decision D-18 build order). Two *additive*
Processes that fill the produced-only ``esters`` and ``fusels`` pools the schema
gained in the byproducts beat. Both are trace (mg/L–low-hundreds-mg/L) beside the
g/L ethanol flux, and both **rise with temperature** — the physics behind the
"warm ferments are estery/fusel-heavy" rule and the directional benchmark
``test_lower_temperature_is_slower_but_cleaner``.

**Additive, not multiplicative.** Unlike ethanol inhibition / Arrhenius (which
*scale* an existing flux and so are :class:`~fermentation.core.process.\
RateModifier` objects), these *produce* a compound, so they are ordinary summed
:class:`~fermentation.core.process.Process` objects (handoff §3.2). Each embeds its
own temperature sensitivity via :func:`~fermentation.core.kinetics.arrhenius.\
arrhenius_factor` (the shared Arrhenius shape, D-11) with its *own* activation
energy — kept in the Process rather than as a separate modifier so each byproduct
mechanism is one self-contained, independently-togglable object.

**Why the benchmark needs steeper-than-flux activation energies.** Each byproduct
rate is the fermentative-flux Monod shape times its Arrhenius factor, so over a run
to dryness the *total* produced scales as ``f_byproduct(T) / f_flux(T) =
exp(-(ΔE_a/R)·(1/T - 1/T_ref))`` (the flux integral to dryness is fixed — the must
ferments out either way). That total falls with temperature **iff** the byproduct's
activation energy exceeds the fermentative flux's (``ΔE_a > 0``): a byproduct that
shared the flux's temperature sensitivity would integrate to a T-independent total
and the "cleaner-when-colder" half of the benchmark would not hold. So
``E_a_esters`` / ``E_a_fusels`` are sourced **above** ``E_a_uptake`` — which is also
the physically-correct ordering (ester/fusel-forming steps are more
temperature-sensitive than the fermentative flux). The cancellation is exact only
for pure uptake to exact dryness; growth's sugar draw, the finite dryness cutoff,
and the (un-Arrhenius-scaled) inactivation brake perturb it, so the integrated
direction is **verified empirically**, not assumed.

**Carbon accounting — option (a)/a1 (decision D-19): carbon routed from sugar.**
Each Process draws its species' carbon *out of ``S``* and the pools are weighted in
``total_carbon`` (by ethyl-acetate / isoamyl-alcohol carbon fractions), so esters and
fusels are **real carbon-accounted state** under one consistent rule with ``Gly`` and
``Byp`` (D-16) — not diagnostic re-expressions. This is the user's call (2026-06-29),
chosen over the interim option (b) (pools outside ``total_carbon``) and the
closure-neutral a2 variant (transfer from ``E``/``Byp`` without a sugar draw).

*The draw touches only ``S`` — never ``E``/``CO2``.* The uptake Process still ferments
``S`` to ethanol+CO2 unchanged; these Processes pull an *additional* sliver of ``S``
into their own pool. So at the derivative level only ``dS`` gains a (tiny, negative)
term — ``dX``/``dN``/``dE``/``dCO2`` are byte-for-byte identical with the byproducts
off. The integrated core therefore drifts only by the trace sugar these consume
(~0.2 % of ``S0``); the per-RHS carbon drawn from ``S`` exactly equals the carbon
deposited in the pool, so ``total_carbon`` closes to machine precision (verified by
the conservation tests). The carbon is split across sugar slots in proportion to each
slot's carbon content (:func:`_draw_carbon_from_sugar`), so the same code serves
wine's one slot and beer's three.

**The ``Byp`` double-count, resolved.** ``Byp`` formerly lumped "organic acids +
higher alcohols" (booked as succinic acid). Fusels *are* higher alcohols, so weighting
a separate ``fusels`` pool on top would book that carbon twice. Under D-19 ``Byp`` is
re-anchored to **organic acids / polyols only** (``Y_byproduct_sugar`` reduced to drop
the higher-alcohol share); the higher alcohols now live solely in the carbon-routed
``fusels`` pool. No overlap remains.

**Both carbon sources are bookkeeping stand-ins, not metabolic claims (D-19).**
(i) The Ehrlich pathway builds fusels from *amino-acid* skeletons, but ``N`` (YAN)
carries no carbon in :func:`total_carbon`, so fusel carbon is sourced from sugar.
(ii) An ester's ethanol moiety is carbon *already counted in ``E``*, so routing ester
carbon from sugar over-attributes fresh hexose. Both close the ledger exactly; neither
asserts where the carbon physically came from. Fusels carry **no CO2 co-product**
(the Ehrlich decarboxylation is omitted) — a documented simplification that keeps the
draw a clean 1:1 sugar→pool carbon transfer.

Tiers: :class:`EsterSynthesis` is **plausible** in form (warmth-favoured,
flux-coupled ester synthesis is the literature-standard direction) with speculative
rate parameters; :class:`FuselAlcoholsEhrlich` is **speculative** in form because
its nitrogen dependence is knowingly simplified to a single monotone branch (the
real Ehrlich relationship is non-monotonic — handoff §3.2). Parameter-tier
propagation (D-1) caps the pool outputs at speculative regardless. These earn no
promotion past plausible: the benchmark is a *directional* check (handoff §3.5).
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Representative species whose formula carbon-accounts each aroma pool (D-19). Ester
#: ⇒ ethyl acetate (C4H8O2); fusel ⇒ isoamyl alcohol (C5H12O). The carbon mass
#: fraction of each weights both the sugar draw here and the pool in ``total_carbon``,
#: from the one chemistry source of truth — so the draw and the conservation check
#: can never disagree (cf. the ``Gly``/``Byp`` routing in the uptake Process).
_ESTER_SPECIES = "ethyl_acetate"
_FUSEL_SPECIES = "isoamyl_alcohol"


def _draw_carbon_from_sugar(
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


def _fermentative_flux_shape(y: FloatArray, schema: StateSchema, k_sat: float) -> float:
    """Biomass-catalysed sugar Monod term ``X · S_total/(K + S_total)`` [g/L].

    The dimensionless-but-for-``X`` activity proxy the fermentative uptake Process
    runs on (``q_sugar_max·X·S/(K+S)``), reused here so byproduct production tracks
    the *same* flux it is metabolically coupled to — which is what makes the
    run-integrated "total scales as f_byproduct/f_flux" cancellation clean and
    predictable (see the module docstring). Sugar is summed across slots (1 for
    wine, 3 for beer) and clamped ≥ 0 against solver undershoot, mirroring the
    guards in the uptake/growth Processes.
    """
    x = max(float(y[schema.slice("X")][0]), 0.0)
    s_total = max(float(y[schema.slice("S")].sum()), 0.0)
    if x <= 0.0 or s_total <= 0.0:
        return 0.0
    return x * (s_total / (k_sat + s_total))


class EsterSynthesis(Process):
    """Ester production, coupled to the fermentative flux and favoured by warmth.

    ``d(esters)/dt = k_ester · X · S_total/(K_sugar_uptake + S_total) · f(T)`` with
    ``f(T) = arrhenius_factor(T, E_a_esters, T_ref)``. Esters (isoamyl acetate,
    ethyl esters, ethyl acetate) form alongside fermentation; tying synthesis to the
    biomass-catalysed sugar flux (sharing ``K_sugar_uptake``) couples them to that
    flux directly, and the steeper-than-uptake ``E_a_esters`` makes the
    run-integrated total fall with temperature (module docstring). The ester carbon
    (booked as ethyl acetate) is routed *out of ``S``* via
    :func:`_draw_carbon_from_sugar` (option a1, D-19), so it touches ``esters`` and
    ``S`` — never ``E``/``CO2`` — and ``total_carbon`` (which now weights ``esters``)
    closes exactly. See the module docstring for the ester carbon-source caveat.
    """

    name = "ester_synthesis"
    tier = Tier.PLAUSIBLE
    touches = ("esters", "S")
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks;
    #: ``E_a_esters`` (> ``E_a_uptake``) and ``T_ref`` set the temperature shape.
    #: Their tiers cap ``esters``'s output tier via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_ester", "K_sugar_uptake", "E_a_esters", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_esters"], params["T_ref"])
        rate = params["k_ester"] * flux * f_t
        d[schema.slice("esters")] = rate
        _draw_carbon_from_sugar(d, y, schema, rate * carbon_mass_fraction(_ESTER_SPECIES))
        return d


class FuselAlcoholsEhrlich(Process):
    """Fusel (higher) alcohols via the Ehrlich pathway — amino-acid-derived, warm.

    ``d(fusels)/dt = k_fusel · X · S_total/(K_sugar_uptake + S_total) ·
    N/(K_n + N) · f(T)`` with ``f(T) = arrhenius_factor(T, E_a_fusels, T_ref)``. The
    Ehrlich pathway makes higher alcohols by transamination/decarboxylation of
    amino acids, so production is gated on yeast-assimilable nitrogen availability
    (``N/(K_n + N)``, sharing the growth nitrogen half-saturation ``K_n``) on top of
    the fermentative flux — front-loading fusels into the nitrogen-replete early
    ferment, as observed. Steeper-than-uptake ``E_a_fusels`` gives the
    warmer-is-more direction. The fusel carbon (booked as isoamyl alcohol) is routed
    *out of ``S``* via :func:`_draw_carbon_from_sugar` (option a1, D-19), so it
    touches ``fusels`` and ``S`` — never ``E``/``CO2`` — and ``total_carbon`` (which
    now weights ``fusels``) closes exactly. The sugar source is a bookkeeping
    stand-in for the amino-acid skeleton (``N`` carries no carbon); see the module
    docstring.

    **Known simplification — monotone in nitrogen (why this Process is speculative).**
    The real fusel/nitrogen relationship is *non-monotonic*: higher alcohols rise
    again at very low YAN (the anabolic/biosynthetic route from sugar) as well as
    with amino-acid catabolism at higher YAN. v1 models only the catabolic,
    monotone-increasing-in-N branch and omits the low-N biosynthetic rise. Per the
    handoff (§3.2) this *models the pathway* (Ehrlich, N-gated) rather than fitting a
    slope, but the single-branch shape is the simplification that keeps the Process
    at the speculative tier until the non-monotonic form is sourced.
    """

    name = "fusel_alcohols_ehrlich"
    tier = Tier.SPECULATIVE
    touches = ("fusels", "S")
    #: ``K_sugar_uptake``/``K_n`` are shared with the uptake/growth Processes;
    #: ``E_a_fusels`` (> ``E_a_uptake``) and ``T_ref`` set the temperature shape.
    reads: tuple[str, ...] = ("k_fusel", "K_sugar_uptake", "K_n", "E_a_fusels", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        n = max(float(y[schema.slice("N")][0]), 0.0)
        if n <= 0.0:  # Ehrlich needs assimilable nitrogen (amino acids)
            return d
        nitrogen_gate = n / (params["K_n"] + n)
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_fusels"], params["T_ref"])
        rate = params["k_fusel"] * flux * nitrogen_gate * f_t
        d[schema.slice("fusels")] = rate
        _draw_carbon_from_sugar(d, y, schema, rate * carbon_mass_fraction(_FUSEL_SPECIES))
        return d
