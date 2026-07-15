"""Temperature-/metabolism-driven aroma byproducts — esters and fusel alcohols.

The first Milestone-2 (Tier-2) beat (decision D-18 build order). Two *additive*
Processes that fill the produced-only ester and ``fusels`` pools the schema gained in the
byproducts beat. Since **decision D-96** "the ester pool" is three single-molecule pools —
``ethyl_acetate`` / ``isoamyl_acetate`` / ``ethyl_hexanoate``, registered in
:data:`~fermentation.core.kinetics.carbon_routing.ESTER_SPECS` — not one lump; see that
registry for why the lump could not survive contact with its own OAV readout. Both are
trace (mg/L–low-hundreds-mg/L) beside the g/L ethanol flux, and both **rise with
temperature** — the physics behind the
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
``total_carbon`` (each ester by its OWN molecule's carbon fraction since D-96 — C4/C7/C8 —
and fusels by isoamyl alcohol's), so esters and fusels are **real carbon-accounted state**
under one consistent rule with ``Gly`` and ``Byp`` (D-16) — not diagnostic re-expressions.
This is the user's call (2026-06-29),
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

**The fusel sugar stand-in is now re-routable (decision D-33).** When the toggleable
``amino_acids`` pool (arginine; D-32) is dosed, :class:`FuselAminoAcidReroute` re-sources a
fraction of the fusel carbon off the sugar stand-in and onto that amino-acid pool — the
*physically-faithful* Ehrlich source — and **deaminates**, releasing the consumed amino
acids' nitrogen to the ammonium ``N`` pool. It is a separate wine-only *swap* Process (it
never touches ``fusels``; production stays here), sharing the one
:func:`fusel_production_rate` so its sugar refund matches this producer's draw exactly. That
deamination branch is the prerequisite the re-route was long deferred on (D-19/D-32); see
:class:`FuselAminoAcidReroute` for the closure algebra and the arginine N-over-release caveat.

**The gas-stripping sink (decisions D-20 → D-21).** The observed fall of wine *liquid*
ester with temperature is largely **evaporation** (Rollero 2014), not reduced synthesis.
That sink — logged as future work in D-19 — is built as :class:`EsterVolatilization`,
which strips each liquid ester into its own bookkeeping headspace pool on the
evolving-CO2 stream, with a **physical** Henry's-law partition (``dH_ester_volatil`` ≈
45 kJ/mol, sourced ethyl-acetate gas/liquid partition enthalpy; D-21 replaced D-20's
fudged per-medium ``E_a_ester_volatil``). Because a sourced stripping is medium-
independent, the wine/beer direction is carried by **per-medium sourced synthesis
``E_a_esters``** (beer steep / wine flat), not by the sink. Each transfer is carbon-neutral
(a pool and its twin share one molecule's weight), so ``total_carbon`` still closes to
machine precision. See that class for the full rationale, and for why v1 shares one sourced
ethyl-acetate partition enthalpy across all three esters (a documented approximation, D-96).

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

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import ESTER_SPECS
from fermentation.core.kinetics.carbon_routing import (
    draw_carbon_from_sugar as _draw_carbon_from_sugar,
)
from fermentation.core.kinetics.carbon_routing import (
    fermentative_flux_shape as _fermentative_flux_shape,
)
from fermentation.core.kinetics.carbon_routing import (
    refund_carbon_to_sugar as _refund_carbon_to_sugar,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Representative species whose formula carbon-accounts each aroma pool (D-19). The carbon
#: mass fraction of each weights both the sugar draw here and the pool in ``total_carbon``,
#: from the one chemistry source of truth — so the draw and the conservation check
#: can never disagree (cf. the ``Gly``/``Byp`` routing in the uptake Process).
#:
#: The fusel pool is still a genuinely LUMPED pool standing in for the higher alcohols;
#: the ester pools no longer are — each is its own molecule, registered in ``ESTER_SPECS``
#: and weighted by itself (decision D-96), so there is no single ``_ESTER_SPECIES`` any more.
_FUSEL_SPECIES = "isoamyl_alcohol"


def fusel_production_rate(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Ehrlich fusel-alcohol production rate ``d(fusels)/dt`` [g/L/h] — the shared rate.

    ``k_fusel · X·S_total/(K_sugar_uptake+S_total) · N/(K_n+N) · arrhenius(T, E_a_fusels,
    T_ref)`` (the fermentative-flux Monod shape, gated on assimilable nitrogen, warmed by a
    steeper-than-uptake Arrhenius factor). Returns 0 under the same guards
    :class:`FuselAlcoholsEhrlich` applies (no flux, no nitrogen).

    Factored out (decision D-33) as the single source of the fusel rate so the *producer*
    (:class:`FuselAlcoholsEhrlich`, which deposits it in the ``fusels`` pool and draws the
    carbon from sugar) and the *re-route* (:class:`FuselAminoAcidReroute`, which re-sources a
    fraction of that carbon from the amino-acid pool) compute the **identical** rate. Any
    divergence between the two would break carbon closure, since the re-route refunds exactly
    the sugar carbon the producer drew — this helper makes that impossible (the shared
    ``biomass_growth_rate`` discipline of the D-32 swap, applied to fusels).
    """
    flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
    if flux <= 0.0:
        return 0.0
    n = max(float(y[schema.slice("N")][0]), 0.0)
    if n <= 0.0:  # Ehrlich needs assimilable nitrogen (amino acids)
        return 0.0
    nitrogen_gate = n / (params["K_n"] + n)
    temp = float(y[schema.slice("T")][0])
    f_t = arrhenius_factor(temp, params["E_a_fusels"], params["T_ref"])
    return float(params["k_fusel"] * flux * nitrogen_gate * f_t)


class EsterSynthesis(Process):
    """Ester production, coupled to the fermentative flux and favoured by warmth.

    For each ester in :data:`~fermentation.core.kinetics.carbon_routing.ESTER_SPECS`,
    ``d(ester)/dt = k_<ester> · X · S_total/(K_sugar_uptake + S_total) · f(T)`` with
    ``f(T) = arrhenius_factor(T, E_a_esters, T_ref)``. Esters form alongside fermentation;
    tying synthesis to the biomass-catalysed sugar flux (sharing ``K_sugar_uptake``) couples
    them to that flux directly. ``E_a_esters`` is **sourced per medium** (decision D-21):
    **steep** for beer (de Andrés-Toro 1998, ester ride growth — synthesis outruns the
    stripping sink, so beer liquid esters *rise* with T) and **weak/~flat** for wine (Mouret
    2015 / Rollero 2014, wine ester synthesis is weak and non-monotonic — so the gas-stripping
    sink :class:`EsterVolatilization` wins and wine liquid esters *fall* with T). This is only
    the *synthesis* term; net liquid ester is synthesis minus stripping (see that class and
    the module docstring). Each ester's carbon is routed *out of ``S``* via
    :func:`_draw_carbon_from_sugar` **at its own molecule's carbon fraction** (option a1,
    D-19), so this touches the three ester pools and ``S`` — never ``E``/``CO2`` — and
    ``total_carbon`` closes exactly. See the module docstring for the ester carbon-source
    caveat.

    **Three esters, three independently-sourced rates (decision D-96).** Before D-96 a single
    ``k_ester`` filled one lumped pool. The pool is now split into three single-molecule pools
    and — the load-bearing call — each ``k`` is anchored to **its own molecule's** measured
    concentration range rather than to a share of the old total. Splitting one ``k_ester`` by a
    fitted ratio would have reproduced the same numbers while smuggling back the very
    fabricated-composition constant the split exists to remove; here the lump's composition and
    the total ester mass are *derived*, not targeted.

    **Shared temperature shape — a documented simplification (D-96).** All three read the one
    per-medium ``E_a_esters``. For the two acetates this is principled: same enzyme (ATF1), so
    no basis to differ. For ``ethyl_hexanoate`` (EEB1/EHT1 — a different enzyme, tied to
    fatty-acid synthesis) it is a genuine approximation, kept because no separate activation
    energy is sourced; a per-ester ``E_a`` is the named deferred refinement.

    **The ATF1 precursor coupling (decision D-97) — why isoamyl acetate alone carries an extra
    term.** ATF1 acetylates an *alcohol*, and Fujii 1998 (Appl. Environ. Microbiol.
    64:4076-4078) measures its ``Km`` for isoamyl alcohol at ~29.8 mM while stating outright
    that *"a major rate-limiting factor for isoamyl acetate production is the amount of isoamyl
    alcohol"*. The ``fusels`` pool runs ~0.5-1 mM — **~30-60x below that Km** — so the enzyme
    sits far down its linear stretch and the rate is **first-order in the pool**:
    ``d(isoamyl_acetate)/dt = k · fusels · X · S/(K_su+S) · f(T)``. Ethyl acetate gets no such
    term because *its* precursor is ethanol (~2 M, orders of magnitude above any mM-scale Km)
    ⇒ ATF1 is saturated in it ⇒ zeroth-order. Same enzyme, same rate law, opposite limits: the
    asymmetry is **derived from the two precursors' concentrations**, not an exemption.

    **Why first-order and not an explicit ``fusels/(Km+fusels)`` gate — identifiability, not
    parsimony (D-97).** In the ``[S] << Km`` limit the saturable form is numerically the linear
    one to within ~2 %, and only the *ratio* ``Vmax/Km`` is identifiable from any model output:
    scaling ``Km`` tenfold and refitting ``Vmax`` would give a byte-identical trajectory. Adding
    the measured ``Km`` as a parameter would therefore mint a sourced-looking constant that no
    model output could ever validate. ``k_isoamyl_acetate`` **is** that identifiable ratio, and
    the ``Km`` lives where it does its real work: in the parameter's provenance, as the sourced
    justification for this rate law's *form*.

    **What this buys — the banana note becomes YAN-responsive.** Reading the ``fusels`` pool
    makes isoamyl acetate inherit the nitrogen dependence of the Ehrlich pathway that builds its
    precursor: a low-YAN ferment makes less isoamyl alcohol, hence less banana. Before D-97 the
    ester was **YAN-blind** (flat ~0.758 mg/L across YAN 40-250 mg/L, where ``fusels`` swung
    2.9x) because every ester shared one plain flux shape. Note this couples to the *pool*, not
    to the fusel *production rate*: the alcohol **persists** after nitrogen is exhausted (in a
    wine run ``N`` empties around day 2 with ~75 % of the sugar still unfermented), and ATF1 goes
    on acetylating it for as long as the flux supplies acetyl-CoA. Coupling to the rate would
    have stopped the banana dead at day 2 with 51 mg/L of its substrate sitting in the vessel.

    **Read, never debited — a deliberate v1 scope call (D-97).** The rate reads ``fusels``; the
    ester's carbon is still drawn from ``S``, so ``touches`` is unchanged and ``total_carbon``
    is untouched. Physically the acetylation takes the C5 skeleton *from* the alcohol (and C2
    from acetyl-CoA) — the exact inverse of the 5:2 split D-69's hydrolysis returns to
    ``fusels`` — so a carbon re-route off ``fusels`` is the honest end state and is the named
    deferred refinement. It is mass-negligible here (~0.5 mg/L of ester against an ~86 mg/L
    fusel pool), and the sugar draw is the same documented stand-in
    :class:`FuselAlcoholsEhrlich` already uses for the amino-acid skeleton.

    **Two inherited caveats, named not buried (D-97).** (i) ``fusels`` is a **lumped** pool
    (isoamyl alcohol representative, but really all the higher alcohols), so reading it whole as
    the ATF1 substrate over-states the true isoamyl-alcohol supply. The over-statement is
    absorbed into the re-anchored ``k``; what it costs is that the *YAN-response* assumes the
    lump's composition is fixed — the D-66 caveat, honest here because ``fusels`` is genuinely
    flagged lumped. (ii) Isoamyl acetate now inherits :class:`FuselAlcoholsEhrlich`'s
    **speculative monotone-in-N** shape, including its admitted omission of the low-YAN
    biosynthetic rise. Both sit inside this Process's existing plausible/speculative framing.
    """

    name = "ester_synthesis"
    tier = Tier.PLAUSIBLE
    touches = (*(spec.pool for spec in ESTER_SPECS), "S")
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks;
    #: ``E_a_esters`` (sourced per medium — steep beer / flat wine, D-21) and ``T_ref``
    #: set the temperature shape, shared by all three esters (D-96). Each ester has its own
    #: independently-sourced rate constant. Their tiers cap the ester pools' output tier via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        *(spec.k_param for spec in ESTER_SPECS),
        "K_sugar_uptake",
        "E_a_esters",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = _fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_esters"], params["T_ref"])
        for spec in ESTER_SPECS:
            rate = params[spec.k_param] * flux * f_t
            if spec.precursor_pool is not None:
                # ATF1 is FAR from saturated in this alcohol ([S] ~ 0.5-1 mM vs Km ~29.8 mM),
                # so its supply limits the rate first-order (D-97). READ ONLY — the pool is
                # never debited; the ester's carbon still comes from S (see the class doc).
                # Clamped >= 0 so a solver undershoot cannot flip synthesis negative.
                rate *= max(float(y[schema.slice(spec.precursor_pool)][0]), 0.0)
            d[schema.slice(spec.pool)] = rate
            # Each ester draws at ITS OWN carbon fraction (C4/C7/C8) — the D-96 split's
            # ledger payoff: no ester's carbon is booked through a stand-in molecule.
            _draw_carbon_from_sugar(d, y, schema, rate * carbon_mass_fraction(spec.species))
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
        rate = fusel_production_rate(y, schema, params)  # shared with the re-route (D-33)
        if rate <= 0.0:
            return d
        d[schema.slice("fusels")] = rate
        _draw_carbon_from_sugar(d, y, schema, rate * carbon_mass_fraction(_FUSEL_SPECIES))
        return d


class FuselAminoAcidReroute(Process):
    """Ehrlich fusel carbon re-sourced from amino acids, with deamination (decision D-33).

    The physically-faithful other half of :class:`FuselAlcoholsEhrlich`. That producer books
    fusel carbon out of **sugar** — a documented bookkeeping stand-in, because the real Ehrlich
    pathway builds higher alcohols from *amino-acid* skeletons (transamination →
    decarboxylation → reduction), releasing the amino group as ammonium. Sugar is used only
    because ``N`` (YAN) carries no carbon in :func:`total_carbon` (D-19), so there was nowhere
    else to draw it from — *until* the toggleable ``amino_acids`` pool (arginine; D-32) gave
    the model a carbon- *and* nitrogen-bearing amino-acid source. This Process closes that gap:
    when amino acids are present it **re-routes** a fraction of the fusel carbon off the sugar
    stand-in and onto the amino-acid pool, and **deaminates** — releasing the nitrogen the
    consumed amino acids carried back to the ammonium ``N`` pool. That deamination branch is
    exactly the prerequisite the fusel re-route was deferred on (D-19/D-32); it also demonstrates
    the aa→N release path a later MLF-with-growth model needs.

    **A swap, not a producer — it never touches ``fusels``.** Like the D-32
    :class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation` swap, this leaves the
    *production* entirely to :class:`FuselAlcoholsEhrlich` (both call the one
    :func:`fusel_production_rate`); it only moves the carbon *source*. For the amino-acid-sourced
    fraction ``g = aa/(K_amino_acids + aa)`` (the same smooth availability gate the swap uses,
    → 0 as the pool empties) of the fusel carbon ``F_c = rate·c_fusel``:

      * **refund sugar** by ``g·F_c`` (undoing the producer's draw for that fraction),
      * **debit amino acids** by ``g·F_c / c_aa`` (the arginine mass carrying that carbon), and
      * **release ammonium** ``N`` by ``(g·F_c/c_aa)·y_N`` (deamination).

    Carbon closes: the fusel gains ``F_c`` (from the producer), sourced now as ``(1−g)·F_c`` from
    sugar + ``g·F_c`` from amino acids. Nitrogen closes: the amino acids lose exactly the nitrogen
    the ``N`` pool gains. Net sugar is ``−(1−g)·F_c ≤ 0`` for all ``g ≤ 1`` — the re-route never
    creates sugar (it only *spares* it), so the ABV bookkeeping caveat is the D-32 one (spared
    sugar ferments to ethanol). **Wine-only** and **forced to be a separate Process**: declaring
    ``amino_acids``/``N`` in the both-media producer's ``touches`` would break beer's ProcessSet
    construction (beer has no ``amino_acids`` slot).

    **Isolability (undosed-only, paired with the producer).** The availability gate → 0 at
    ``aa = 0``, so an undosed wine run is byte-for-byte the sugar-stand-in producer; the compile
    seam additionally *disables* this Process when ``amino_acids_gpl ≤ 0`` (tier isolability, the
    D-32 pattern). It is only valid while :class:`FuselAlcoholsEhrlich` is active — it refunds
    sugar that producer drew — so the two are kept paired (disabling the producer alone would let
    the re-route create sugar; the same acceptable swap↔producer coupling as D-32's swap↔growth).

    **Documented lump — arginine over-releases nitrogen.** Sourcing fusel carbon through the
    N-rich representative amino acid deaminates ``c_fusel/c_aa · y_N`` ≈ 0.78 g N per g fusel
    carbon — roughly **4× the real leucine→isoamyl-alcohol N:C** (leucine carries one amino group
    over six carbons). This is conservation-exact but a forced consequence of the single-species
    ``amino_acids`` lump (arginine, chosen N-rich for the D-32 swap), the same class of stand-in
    as the sugar-carbon fiction it replaces. The released nitrogen feeds back as supplementary YAN,
    but fusels are trace so the effect is second-order and tiny. Tier **speculative** (it inherits
    the fusel rate's speculative parameters and the ``amino_acids`` gate estimate).
    """

    name = "fusel_amino_acid_reroute"
    tier = Tier.SPECULATIVE
    #: Refunds carbon to ``S``, debits ``amino_acids``, releases nitrogen to ``N``. Never
    #: touches ``fusels`` — production stays entirely in :class:`FuselAlcoholsEhrlich`.
    touches = ("S", "amino_acids", "N")
    #: Recomputes the fusel rate (so it reads the producer's parameters) plus ``K_amino_acids``
    #: for the availability gate. Their tiers cap the ``S``/``amino_acids``/``N`` output tier via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_fusel",
        "K_sugar_uptake",
        "K_n",
        "E_a_fusels",
        "T_ref",
        "K_amino_acids",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        rate = fusel_production_rate(y, schema, params)  # identical to the producer's (D-33)
        if rate <= 0.0:
            return d
        aa = max(float(y[schema.slice("amino_acids")][0]), 0.0)
        if aa <= 0.0:
            return d  # empty pool ⇒ producer sources all fusel carbon from sugar (undosed no-op)

        gate = aa / (params["K_amino_acids"] + aa)  # smooth availability, in [0, 1)
        fusel_carbon = rate * carbon_mass_fraction(_FUSEL_SPECIES)  # what the producer drew from S
        aa_carbon = gate * fusel_carbon  # the fraction re-sourced from amino acids
        c_aa = carbon_mass_fraction(AMINO_ACID_SPECIES)
        y_n = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
        aa_mass = aa_carbon / c_aa  # arginine mass consumed to supply that carbon

        d[schema.slice("amino_acids")] = -aa_mass
        d[schema.slice("N")] = (
            aa_mass * y_n
        )  # DEAMINATION: aa nitrogen → ammonium (the D-33 branch)
        # Refund the producer's sugar draw for the re-sourced fraction (the inverse of its draw),
        # so net sugar loss is only the (1−g) fraction still taken from sugar. rate > 0 ⇒ flux > 0
        # ⇒ sugar present, so the refund always lands (no carbon leak).
        _refund_carbon_to_sugar(d, y, schema, aa_carbon)
        return d


class EsterVolatilization(Process):
    """CO2-stripping loss of liquid esters to the headspace — a physical Henry's-law sink.

    ``d(esters)/dt = -k · X·S_total/(K_sugar_uptake+S_total) · f_gas(T) · f_part(T) ·
    <ester>`` and the equal-and-opposite ``+`` into that ester's own ``<ester>_gas`` twin,
    for each ester in ``ESTER_SPECS`` (decision D-96), where

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

    **Carbon — three neutral liquid→gas transfers (no sugar draw).** Unlike
    :class:`EsterSynthesis`/:class:`FuselAlcoholsEhrlich`, this Process draws **no fresh
    sugar**: it only moves carbon already in each liquid ester pool into *that ester's own*
    headspace pool. So it touches the three ester pools and their three gas twins only —
    never ``S``/``E``/``CO2`` — and since ``total_carbon`` weights each pool and its twin at
    the *same* molecule's fraction, **every** transfer is carbon-neutral independently and
    closure stays at machine precision (a headspace pool is the ester analogue of evolved
    ``CO2``: carbon leaves the liquid but not the ledger). Each liquid pool is clamped ≥ 0 so
    a solver undershoot cannot strip a negative pool into spurious gas.

    Note the pairing is *why* each ester needs its own gas pool (decision D-96): a single
    shared headspace pool could carry only one carbon weight, so stripping a C7 ester into a
    C4-weighted pool would create or destroy carbon. Before D-96 one lumped pool and one gas
    twin were both booked as ethyl acetate, and the question could not arise.

    **Documented simplification — one partition enthalpy for three molecules (D-96).**
    ``dH_ester_volatil`` is sourced for **ethyl acetate** (~45 kJ/mol, NIST/Sander) and v1
    applies it to all three esters. That is a measured constant borrowed for two molecules it
    was not measured for — an approximation the ester van't Hoff enthalpies' natural clustering
    (~40–55 kJ/mol across this range) makes tolerable, but an approximation, and it is recorded
    here rather than left implicit. Its direction is known: isoamyl acetate and ethyl hexanoate
    are *less* volatile than ethyl acetate, so v1 slightly over-strips them. Per-ester sourced
    partition enthalpies are the named deferred refinement.

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
    touches = (*(spec.pool for spec in ESTER_SPECS), *(spec.gas_pool for spec in ESTER_SPECS))
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
        temp = float(y[schema.slice("T")][0])
        f_gas = arrhenius_factor(temp, params["E_a_uptake"], params["T_ref"])  # CO2 gas flow
        f_part = arrhenius_factor(temp, params["dH_ester_volatil"], params["T_ref"])  # partition
        for spec in ESTER_SPECS:
            liquid = max(float(y[schema.slice(spec.pool)][0]), 0.0)
            if liquid <= 0.0:  # nothing in this liquid pool to strip
                continue
            rate = params["k_ester_volatil"] * flux * f_gas * f_part * liquid
            # Liquid → its OWN headspace twin: the pair shares one molecule's carbon weight,
            # so each transfer is carbon-neutral independently (D-96).
            d[schema.slice(spec.pool)] = -rate
            d[schema.slice(spec.gas_pool)] = rate
        return d
