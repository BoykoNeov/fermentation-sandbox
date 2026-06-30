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

**Temperature directions — what rises, what falls, and why (decisions D-19 → D-21).**
Each *production* rate is the fermentative-flux Monod shape times its Arrhenius factor,
so over a run to dryness the total produced scales as ``f_byproduct(T) / f_flux(T)``
(the flux integral to dryness is fixed — the must ferments out either way), which rises
with T **iff** the byproduct's activation energy exceeds the fermentative flux's.

* **Fusels (both media)** must rise with T for the "cleaner-when-colder" benchmark, so
  ``E_a_fusels`` is held **above** ``E_a_uptake`` — a *sourced ordering* (Mouret 2015:
  higher alcohols rise with T; brewing/enology consensus). This is the half of the
  benchmark that genuinely carries "cleaner cold" for wine.
* **Esters** are *not* a simple production pool: liquid ester is the balance of
  synthesis (:class:`EsterSynthesis`, ``E_a_esters``) against the gas-stripping sink
  (:class:`EsterVolatilization`, stripping sensitivity ``E_a_uptake + dH_ester_volatil``
  ≈ 100 kJ/mol, the *same physical value in both media*). The wine/beer direction lives
  in **synthesis ``E_a_esters``, set per medium and sourced** — **beer** steep
  (de Andrés-Toro 1998, ester ride growth, apparent E_a ~200+ kJ/mol ⇒ synthesis wins ⇒
  liquid esters *rise* with T) and **wine** weak/~flat (Mouret 2015; Rollero 2014, wine
  ester synthesis is weak and non-monotonic ⇒ stripping wins ⇒ liquid esters *fall* with
  T, the Rollero evaporation inversion). See :class:`EsterVolatilization` for why a
  *sourced* (hence medium-independent) stripping forces the direction into synthesis.

Exact magnitudes stay speculative (directional check only, handoff §3.5); the *orderings*
are sourced. The integrated directions are **verified empirically** at 14/20/25 °C per
medium, not assumed — growth's sugar draw, the finite dryness cutoff, the inactivation
brake, and the stripping dynamics all perturb the clean cancellation.

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

**The gas-stripping sink (decisions D-20 → D-21).** The observed fall of wine *liquid*
ester with temperature is largely **evaporation** (Rollero 2014), not reduced synthesis.
That sink — logged as future work in D-19 — is built as :class:`EsterVolatilization`,
which strips liquid ``esters`` into the bookkeeping ``esters_gas`` headspace pool on the
evolving-CO2 stream, with a **physical** Henry's-law partition (``dH_ester_volatil`` ≈
45 kJ/mol, sourced ethyl-acetate gas/liquid partition enthalpy; D-21 replaced D-20's
fudged per-medium ``E_a_ester_volatil``). Because a sourced stripping is medium-
independent, the wine/beer direction is carried by **per-medium sourced synthesis
``E_a_esters``** (beer steep / wine flat), not by the sink. The transfer is carbon-neutral
(``esters`` → ``esters_gas``, both ethyl acetate), so ``total_carbon`` still closes to
machine precision. See that class for the full rationale.

Tiers: :class:`EsterSynthesis` is **plausible** in form (warmth-favoured,
flux-coupled ester synthesis is the standard direction in the canonical *beer* model,
de Andrés-Toro 1998) with speculative rate parameters. *Wine caveat:* in wine the
warmer⇒more-ester direction is weaker and partly confounded — ester *synthesis* is
non-monotonic in temperature, so the warmer⇒more-aroma direction is carried by the
fusels while the **liquid** ester fall with warmth is the job of the volatilization sink
above (Mouret 2015; Rollero 2014). :class:`FuselAlcoholsEhrlich`
is **speculative** in form because its nitrogen dependence is knowingly simplified to a
single monotone branch (the real Ehrlich relationship is non-monotonic — handoff §3.2,
corroborated by Mouret 2015 / Rollero 2014: higher-alcohol synthesis is optimal at
~200–300 mg N/L). Parameter-tier propagation (D-1) caps the pool outputs at speculative
regardless. These earn no promotion past plausible: the benchmark is a *directional*
check (handoff §3.5).
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
    flux directly. ``E_a_esters`` is **sourced per medium** (decision D-21): **steep**
    for beer (de Andrés-Toro 1998, ester ride growth — synthesis outruns the stripping
    sink, so beer liquid esters *rise* with T) and **weak/~flat** for wine (Mouret 2015 /
    Rollero 2014, wine ester synthesis is weak and non-monotonic — so the gas-stripping
    sink :class:`EsterVolatilization` wins and wine liquid esters *fall* with T). This is
    only the *synthesis* term; net liquid ester is synthesis minus stripping (see that
    class and the module docstring). The ester carbon (booked as ethyl acetate) is routed
    *out of ``S``* via :func:`_draw_carbon_from_sugar` (option a1, D-19), so it touches
    ``esters`` and ``S`` — never ``E``/``CO2`` — and ``total_carbon`` (which now weights
    ``esters``) closes exactly. See the module docstring for the ester carbon-source caveat.
    """

    name = "ester_synthesis"
    tier = Tier.PLAUSIBLE
    touches = ("esters", "S")
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks;
    #: ``E_a_esters`` (sourced per medium — steep beer / flat wine, D-21) and ``T_ref``
    #: set the temperature shape.
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


class EsterVolatilization(Process):
    """CO2-stripping loss of liquid esters to the headspace — a physical Henry's-law sink.

    ``d(esters)/dt = -k · X·S_total/(K_sugar_uptake+S_total) · f_gas(T) · f_part(T) ·
    esters`` and the equal-and-opposite ``+`` into ``esters_gas``, where

    * ``f_gas(T) = arrhenius_factor(T, E_a_uptake, T_ref)`` is the **gas-flow** factor: the
      stripping rides the evolving-CO2 stream, whose rate is the fermentative uptake flux
      scaled by the *same* Arrhenius factor the uptake Process carries (``E_a_uptake``).
      The bare ``X·S/(K+S)`` shape times ``f_gas`` is the CO2-evolution proxy; the
      ``q_sugar_max``/``co2_yield``/realised-yield constants fold into ``k_ester_volatil``.
    * ``f_part(T) = arrhenius_factor(T, dH_ester_volatil, T_ref)`` is the **gas/liquid
      partition** factor — a van't Hoff form for the ethyl-acetate Henry's-law constant,
      which *rises* with temperature (warmer ⇒ more volatile). ``dH_ester_volatil`` is the
      sourced partition enthalpy (~45 kJ/mol, ethyl acetate; NIST/Sander Henry data), a
      **physical Q10 ≈ 1.8**, *not* a fitted lever.

    So the total stripping temperature sensitivity is ``E_a_uptake + dH_ester_volatil`` ≈
    100 kJ/mol — the same physical value **in both media**, because the Henry's-law
    partition and the gas flow are properties of the molecule and the ferment, not the
    beverage (Morakul et al. 2011: the partition coefficient depends only on composition
    and temperature, not on which medium). The loss is **first-order in the liquid ester
    present** and **stops when fermentation stops** (``flux → 0`` at dryness) — a
    deliberate omission of slow passive evaporation after the cap goes on, keeping the
    sink a clean function of the gas stream (the ester analogue of Morakul's
    ``L = ∫ C_gas·Q_CO2 dt``).

    **Why a shared physical stripping is the faithful model (decisions D-20 → D-21).** D-20
    first built this sink with a *fudged per-medium* ``E_a_ester_volatil`` (wine above /
    beer below ``E_a_esters``) to make wine liquid esters fall while beer rises. But a
    *sourced* Henry's-law stripping is medium-independent, so it cannot push opposite
    directions by itself. The wine/beer divergence therefore lives where it is genuinely
    sourced — in ester **synthesis** ``E_a_esters``, set *per medium*: **beer** strongly
    T-sensitive (de Andrés-Toro 1998, ester ride the growth rate, apparent E_a ~200+
    kJ/mol) so synthesis outruns stripping and **beer liquid esters rise** with T; **wine**
    weak/~flat (Mouret 2015; Rollero 2014, wine ester synthesis is weak and non-monotonic
    in T) so stripping outruns synthesis and **wine liquid esters fall** with T (the
    Rollero evaporation inversion). Both directions now emerge from *physical + sourced*
    parameters, not a compensating constant. Verified empirically at 14/20/25 °C.

    **Carbon — a neutral liquid→gas transfer (no sugar draw).** Unlike
    :class:`EsterSynthesis`/:class:`FuselAlcoholsEhrlich`, this Process draws **no fresh
    sugar**: it only moves carbon already in the liquid ``esters`` pool into the
    ``esters_gas`` headspace pool, both booked as ethyl acetate. So it touches
    ``esters`` and ``esters_gas`` only — never ``S``/``E``/``CO2`` — and since
    ``total_carbon`` weights both pools at the *same* ethyl-acetate fraction, the transfer
    is carbon-neutral and closure stays at machine precision (the headspace pool is the
    ester analogue of evolved ``CO2``: carbon leaves the liquid but not the ledger).
    ``esters`` is clamped ≥ 0 so a solver undershoot cannot strip a negative pool into
    spurious gas.

    **Documented simplification.** The full Morakul (2011) partition coefficient is also
    *ethanol-dependent* (``ln k_i = F1 + F2·E − (F3 + F4·E)·R·(1000/T − 1000/T_ref)``); we
    keep only the dominant temperature (van't Hoff) lever via ``dH_ester_volatil`` and omit
    the ethanol terms (the ``F`` coefficients are not openly available). Tier: **plausible**
    in form (gas stripping by the evolving CO2 is well-understood physics, the standard
    explanation for wine's liquid-ester temperature response), with speculative rate
    parameters that cap the pool outputs at speculative via parameter-tier propagation
    (D-1).
    """

    name = "ester_volatilization"
    tier = Tier.PLAUSIBLE
    touches = ("esters", "esters_gas")
    #: ``K_sugar_uptake``/``E_a_uptake`` are shared with the fermentative uptake whose CO2
    #: stream does the stripping (gas-flow factor); ``dH_ester_volatil`` is the sourced
    #: ethyl-acetate Henry's-law partition enthalpy (gas/liquid factor); ``T_ref`` anchors
    #: both. Their tiers cap the pool outputs via parameter-tier propagation (D-1). The
    #: wine/beer direction lives in ``E_a_esters`` (synthesis), not here — this stripping
    #: is the same physical mechanism in both media.
    reads: tuple[str, ...] = (
        "k_ester_volatil",
        "K_sugar_uptake",
        "E_a_uptake",
        "dH_ester_volatil",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        esters_liquid = max(float(y[schema.slice("esters")][0]), 0.0)
        if esters_liquid <= 0.0:  # nothing in the liquid pool to strip
            return d
        temp = float(y[schema.slice("T")][0])
        f_gas = arrhenius_factor(temp, params["E_a_uptake"], params["T_ref"])  # CO2 gas flow
        f_part = arrhenius_factor(temp, params["dH_ester_volatil"], params["T_ref"])  # partition
        rate = params["k_ester_volatil"] * flux * f_gas * f_part * esters_liquid
        d[schema.slice("esters")] = -rate
        d[schema.slice("esters_gas")] = rate
        return d
