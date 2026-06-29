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

**Carbon accounting — produced-only, NOT carbon-routed (decision D-19).** The beat
plan floated routing these pools' carbon from sugar D-16-style for machine-precision
closure. We deliberately do **not**, because the carbon is *already in the ledger*:
the ``Byp`` minor-byproduct lump is explicitly "organic acids / higher alcohols,
carbon-accounted as succinic acid" — and higher alcohols *are* the fusels — while
ester carbon merely re-expresses ethanol (``E``) and acid (⊂ ``Byp``) carbon already
counted. Adding ``esters``/``fusels`` to ``total_carbon`` would therefore
**double-count**, not close a leak; and routing them from sugar would force carving
higher alcohols back out of ``Byp`` (re-anchoring the sourced ``Y_byproduct_sugar``
and risking the realised-yield / ABV realism guard). So these Processes touch *only*
their own pools, ``total_carbon`` is left untouched, and carbon closure stays
byte-for-byte what it was. (A future beat that wants rigorous per-species carbon
routing must first resolve the ``Byp`` double-count, then decrement ``E``/``Byp``/
sugar accordingly — see D-19.)

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

from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


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
    run-integrated total fall with temperature (module docstring). Produced-only:
    touches ``esters`` alone, carbon left to the existing ledger (D-19).
    """

    name = "ester_synthesis"
    tier = Tier.PLAUSIBLE
    touches = ("esters",)
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
        d[schema.slice("esters")] = params["k_ester"] * flux * f_t
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
    warmer-is-more direction. Produced-only: touches ``fusels`` alone (D-19).

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
    touches = ("fusels",)
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
        d[schema.slice("fusels")] = params["k_fusel"] * flux * nitrogen_gate * f_t
        return d
