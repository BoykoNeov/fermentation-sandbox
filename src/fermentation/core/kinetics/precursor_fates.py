"""The precursors' non-Ehrlich fates — the sink the speciated pools never had (decision D-104).

**What this closes (decision D-100 → D-103 → D-104).** D-100 speciated the amino-acid lump into
eight single-molecule pools so each consumer draws the molecule it actually eats. It left one
sink out, and said so: the D-33/D-99 Ehrlich re-route is each precursor's **only** consumer, so
the model attributed **100% of consumed leucine to isoamyl alcohol**. Real yeast send most of it
to **protein**. Every precursor therefore disappeared by the one route the model owned, whatever
reality's split. This Process is that missing sink.

**It is NOT "protein incorporation", and the difference is the whole beat.** ``f_non_ehrlich_*``
is a **lump of every fate that is not this alcohol** — protein incorporation, downstream
amino-acid biosynthesis, and unrecovered losses. Using the *protein share* instead is the trap
this beat fell into and measured its way out of: Crépin's threonine has **four** fates (38%
protein, ~18% propanol, ~17% proteinogenic isoleucine, 20% unrecovered) and this model has
**two**, so booking 0.38 dumps the other 44% onto propanol. Measured on a Crépin-commensurate
must: the protein share gives propanol **55.7%** catabolic against a sourced **19%** — *worse
than the crude estimate it replaced*; the lump gives **16.4%**. The parameter answers exactly one
question — *how much does not become this alcohol* — and reading it as "protein" will
silently re-break the model.

**The form.** For each Ehrlich alcohol the re-route sources from precursor ``i``::

    ehrlich_draw_i = gate_i · fusel_carbon_i          (the re-route's own draw, unchanged)
    sink_draw_i    = ehrlich_draw_i · f_i / (1 − f_i)

so the consumed precursor splits **exactly** ``f_i : (1 − f_i)`` between "everything else" and
the alcohol, for any gate and any trajectory. The sink debits the precursor, **refunds its carbon
to sugar** and **its nitrogen to ammonium** — the D-32 swap shape, one molecule at a time: the
physical reading is that a precursor spent on biomass *spares* the sugar and ammonium growth's
stoichiometry already charged.

**Why it rides the re-route's draw rather than growth.** Anchoring to growth (``w_i·dX/dt`` at
biomass composition) is the obvious form and it is **measurably wrong**: it makes each precursor's
split fall out of *demand*, and the model's Ehrlich demand pulls the **amino acid**, so the
largest alcohol drains its own precursor hardest. Measured against Crépin, that split is
**monotonically inverted** — model leucine 20.9% to protein against a measured 77–86%, model
order leu<ile<val<thr against Crépin's thr<val<ile<leu, exactly reversed. No biomass composition
repairs it: for leucine to reach 77–86% by demand its protein demand would need ~4× the isoamyl
demand — roughly half the biomass being leucine. Reality escapes the inversion by building
isoamyl from **KIC**, a keto acid that is mostly *de novo*, so the leucine pool never faces that
demand at all. This model has no keto-acid node, so it cannot reproduce that mechanism — it can
only impose the outcome the mechanism produces. Riding the re-route's draw is what makes the
imposed ratio *exact*; it is a bookkeeping anchor, not a claim that Ehrlich catabolism *causes*
protein synthesis.

**The honest consequence (decision D-104): the split is STATIC.** ``f_i`` are measured *fates* at
one must composition and one time point, not rate constants, so they do not respond to must
composition — where D-100's gate, wrong as it was, at least moved. This **retires D-100's "the
anabolic/catabolic ratio is emergent, not a fitted fraction" for this split**, and the reason is
measured rather than preferred: the emergent alternative is inverted, and *an emergent wrong
answer is worse than a sourced static one*. What survives of D-100's claim is the part it was
really about — the **catabolic fraction** of each alcohol is still emergent, because it is
``(1−f_i)`` times a supply/demand ratio the model computes.

**What it earns, measured on a Crépin-commensurate must (180 mg N/L, 28 °C):**

===================  ==========  =========  ====================================
alcohol              no sink     with sink  sourced
===================  ==========  =========  ====================================
propanol             72.9%       **16.4%**  19% (Crépin, [13C]threonine)
isobutanol           43.6%       **6.5%**   5–15% (Rollero, U-13C valine)
isoamyl alcohol      7.4%        **1.35%**  2–8% (Rollero, U-13C leucine) — UNDER
===================  ==========  =========  ====================================

**Isoamyl's miss is sourced, not a loose end — do not tune it.** Crépin measures **23% of
consumed valine reaching KIC → isoamyl alcohol**. Reality feeds isoamyl exogenous carbon from
*both* leucine and valine; this model routes only leucine, so it **must** under-count. The gap is
a missing route (valine → KIC → isoamyl), and closing it means a keto-acid node — named in D-104,
not built here. ``active_amyl_alcohol`` is over-attributed by an unknown margin for the mirror
reason: Crépin does not resolve isoleucine's active-amyl share, so ``f_non_ehrlich_isoleucine``
falls back to the protein share, a *lower bound* on the true lump.

**At fermentation this is an ATTRIBUTION fix, not an aroma fix.** Every alcohol moves <1%
(+0.627%, all five identically, via the nitrogen refund feeding the shared ``N/(K_n+N)`` term).
What changes is what the carbon and nitrogen ledgers *mean*: leucine now goes mostly to biomass,
as it does in a real ferment. Where it should move an output is **aging** — it drives every
precursor to ~0 during AF, so sur-lie Strecker/sotolon become autolysis-dependent.

BOOKKEEPING CAVEAT (the D-19/D-31 stand-in discipline): the sink books the whole lump as if it
became **biomass** — carbon refunded to sugar, nitrogen to ammonium. Crépin's "20% unrecovered"
did not become biomass, so that share is booked to the wrong home. Carbon and nitrogen still
close to machine precision (atoms only move between weighted pools); what is approximate is the
**destination**, not the balance.

**Isolability (prime directive #3, undosed-only).** Every gate → 0 at ``aa_i = 0`` and the draw
is proportional to the re-route's, which is itself 0 undosed — so an undosed wine run is
byte-for-byte the validated core. The compile seam additionally *disables* this Process when
``amino_acids_gpl ≤ 0``, the D-32 tier-isolability pattern. It is only valid while the re-route
is active (it scales that Process's draw), so the two are kept paired.

Tier: **speculative** — it inherits the re-route's speculative rate parameters, and two of its
five fractions (phenylalanine, isoleucine's fallback) are author estimates or bounds.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acid_pools import SPEC_BY_SPECIES, depletion_gate
from fermentation.core.kinetics.byproducts import fusel_carbon_draw_by_species
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS, refund_carbon_to_sugar
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

__all__ = ["NON_EHRLICH_FRACTION_PARAMS", "PrecursorNonEhrlichFates", "non_ehrlich_fraction_param"]


def non_ehrlich_fraction_param(species: str) -> str:
    """The ``f_non_ehrlich_<species>`` parameter naming rule — the one place it is spelled.

    Derived from the species rather than hand-listed so a sixth Ehrlich alcohol added to
    :data:`~fermentation.core.kinetics.carbon_routing.FUSEL_SPECS` cannot silently acquire a
    fraction that no parameter file defines (it fails at load, loudly).
    """
    return f"f_non_ehrlich_{species}"


#: The fraction parameter for every precursor the re-route actually draws — derived from the
#: canonical fusel registry, so this can never drift from the set of alcohols that exist.
#: Methionine is deliberately absent: it has no Ehrlich alcohol, so it has no draw to scale
#: (see the parameter file's note where its entry would be).
NON_EHRLICH_FRACTION_PARAMS: tuple[str, ...] = tuple(
    non_ehrlich_fraction_param(spec.precursor_amino_acid) for spec in FUSEL_SPECS
)


class PrecursorNonEhrlichFates(Process):
    """Every fate of a consumed precursor except its own higher alcohol (decision D-104).

    Draws ``f_i/(1−f_i)`` times the Ehrlich re-route's own per-species draw, so consumed
    precursor splits exactly ``f_i : (1−f_i)`` between the non-Ehrlich lump and the alcohol.
    Debits the precursor, refunds its carbon to ``S`` and its nitrogen to ``N`` (the D-32 swap
    shape, per molecule). See the module docstring for why ``f_i`` is **not** the protein share
    and why the draw rides the re-route rather than growth.
    """

    name = "precursor_non_ehrlich_fates"
    tier = Tier.SPECULATIVE
    #: Debits **each alcohol's own precursor**, refunds carbon to ``S`` and nitrogen to ``N``.
    #: Never touches ``fusels`` (production stays in the producer) and never touches
    #: ``amino_acids``/``amino_acids_generic`` — the D-100 decoupling: arginine does not make
    #: higher alcohols, so this sink cannot starve the yeast swap / MLF / Brett / Maillard
    #: consumers that live on those pools. Same ``touches`` set as the re-route it scales.
    touches = ("S", "N", *(spec.precursor_amino_acid for spec in FUSEL_SPECS))
    #: Recomputes the re-route's per-species draw (so it reads the fusel producer's parameters
    #: plus ``K_amino_acids`` and the ``must_aa_fraction_*`` gates), and its own five
    #: ``f_non_ehrlich_*`` fractions. Their tiers cap the ``S``/precursor/``N`` output tiers via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        *(spec.k_param for spec in FUSEL_SPECS),
        "K_sugar_uptake",
        "K_n",
        "E_a_fusels",
        "T_ref",
        "K_amino_acids",
        *(SPEC_BY_SPECIES[spec.precursor_amino_acid].fraction_param for spec in FUSEL_SPECS),
        *NON_EHRLICH_FRACTION_PARAMS,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # The re-route's own per-species carbon draw, from the identical helper it uses — so the
        # ratio this Process imposes is against the exact flux the re-route will book, and the two
        # can never drift apart (the D-33/D-99 shared-helper discipline).
        draws = fusel_carbon_draw_by_species(y, schema, params)
        if not draws:
            return d
        carbon = 0.0  # total precursor carbon re-sourced off sugar into the non-Ehrlich lump
        nitrogen = 0.0  # total nitrogen it carries, refunded to ammonium
        for spec, fusel_carbon in draws:
            if fusel_carbon <= 0.0:
                continue
            precursor = SPEC_BY_SPECIES[spec.precursor_amino_acid]
            # This precursor's OWN relative-depletion gate — the SAME gate the re-route applies
            # (D-100), so both consumers are throttled together and the imposed split holds at
            # every instant, not just in the integral. → 0 as the pool empties, so the combined
            # draw can never drive it negative, and an undosed run is a no-op.
            gate = depletion_gate(y, schema, params, (precursor,))
            if gate <= 0.0:
                continue
            f = params[non_ehrlich_fraction_param(spec.precursor_amino_acid)]
            # f is a fraction of the consumed precursor, so f ∈ [0, 1); f → 1 would demand an
            # infinite draw against a finite alcohol. The parameter file bounds every entry well
            # below this; the guard is here because a bad ENSEMBLE draw off the uncertainty band
            # is the realistic way it would ever be reached, and a silent inf would poison the
            # solver rather than fail.
            if not 0.0 <= f < 1.0:
                raise ValueError(
                    f"{non_ehrlich_fraction_param(spec.precursor_amino_acid)}={f} outside [0, 1): "
                    "it is the fraction of consumed precursor NOT becoming its alcohol"
                )
            ehrlich_carbon = gate * fusel_carbon  # exactly what the re-route books
            lump_carbon = ehrlich_carbon * f / (1.0 - f)  # ⇒ realised split is exactly f : (1−f)
            mass = lump_carbon / carbon_mass_fraction(precursor.species)
            d[schema.slice(precursor.pool)] -= mass
            nitrogen += mass * nitrogen_mass_fraction(precursor.species)
            carbon += lump_carbon
        if carbon <= 0.0:
            return d
        # The precursor spent on biomass spares the ammonium and sugar growth's stoichiometry
        # already charged — the D-32 swap's physical reading, one molecule at a time. Carbon and
        # nitrogen close exactly; only the DESTINATION is a stand-in (module docstring).
        d[schema.slice("N")] = nitrogen
        refund_carbon_to_sugar(d, y, schema, carbon)
        return d
