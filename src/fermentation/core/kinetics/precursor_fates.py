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

**Isoamyl's miss is sourced, not a loose end — do not tune it. [D-111 built the valine route; D-112
MEASURED that it does not close the leucine gap, and neither does the node.]** The valine → KIC →
isoamyl route D-104 named as missing is real and was built at D-111 (0% → 1.74%, matching Rollero's
valine tracer). But the *leucine*-derived shortfall (1.12% vs Rollero's *leucine* tracer 3.4–17.3%)
is a **different quantity** a valine route cannot touch, and D-112 measured that no sourcing-layer
change touches it either: leucine's only AF fates are this lump ``f`` and isoamyl ``(1−f)``, so
``isoamyl_leucine = (1−f)·leucine_C/isoamyl_C`` is a **mass-conservation ceiling** the model already
sits on, and Crépin's ``f`` already prices in every non-isoamyl fate. Most of the residual gap is an
incommensurate isoamyl **denominator** (the probe must's ``aa`` dose inflates isoamyl ~2×; at
Rollero's isoamyl the share is ~2.9%, ≈ Rollero's floor) plus a raw-enrichment-vs-net-carbon
mismatch. So the keto-acid node's motivation is D-104's **inverted split**, never this gap.
**[D-113 measured whether D-111's valine route touches that inverted split: it does not. Leucine's
Ehrlich branch is bit-invariant under the route — it is a headroom-fill that relieves leucine of 0%
of its isoamyl demand (a valine *drain*, not the leucine *relief* the inversion needs). Un-inverting
leucine remains an unsourced build (de-novo-KIC relief + kinetically-limited transamination).]**
``active_amyl_alcohol`` is over-attributed by an unknown margin for the mirror
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

Tier: **speculative** — it inherits the re-route's speculative rate parameters. **[D-117: all five
fractions now rest on measurements**, phenylalanine's last of all (Minebois 2025, U-13C
phenylalanine). Two of the five still ship as *bounds* rather than values — isoleucine's because its
true lump is unknown, phenylalanine's because its true lump is **known and unshippable** (below).
The Process tier does not move: it is set by the re-route's rates, not by these fractions.**]**

**[D-117 — phenylalanine's fraction was WRONG, not merely unsourced, and the error is this
docstring's own trap.** D-104 shipped ``f_non_ehrlich_phenylalanine = 0.53``. Minebois 2025 measures
Sc sending **2.5%** of consumed phenylalanine to 2-phenylethanol and **53.1%** to protein — so the
lump is **0.975**, and the old 0.53 sat almost exactly on the *protein share*, i.e. on the trap the
paragraph above says this beat "measured its way out of". It fell in again for the one precursor it
could not check, and the coincidence made the error look derived. The stated derivation ("the mean
of Crépin's four measured splits") did not even yield 0.53 — that mean is 0.749. **A value whose
stated derivation does not reproduce it is a defect even while it is only "speculative".**

**AND THE MEASURED LUMP CANNOT BE SHIPPED — this Process is why.** The draw scales ``f/(1−f)``,
which goes **1.13 → 39** between 0.53 and 0.975. At 0.975 the joint (D-32 swap + this sink) carbon
refund reaches **1.125× growth's own draw**: it refunds more carbon than growth was ever charged,
which is gluconeogenesis, and it trips the hard ``< 1.0`` guard in
``test_the_joint_carbon_refund_never_creates_sugar``. The cause is **the missing de-novo route, not
the parameter**: this model charges *all* of its ``k``-calibrated 2-phenylethanol to consumed
phenylalanine, while reality builds ~97% of 2-PE from **de-novo phenylpyruvate**. So the sink must
eat phenylalanine at ~40× the Ehrlich draw to feed an alcohol reality mostly makes from sugar —
**the same shape as the isoamyl/KIC problem this docstring already describes, one precursor over.**
The shipped value is therefore Minebois's protein share, **0.531, an explicit lower bound**, and the
blocker is pinned by ``test_the_sourced_lump_breaks_the_carbon_refund_guard``, which is *designed to
fail* when the phenylpyruvate route lands. See DECISIONS D-117.]**
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acid_pools import SPEC_BY_SPECIES
from fermentation.core.kinetics.byproducts import ehrlich_draws
from fermentation.core.kinetics.carbon_routing import (
    FUSEL_SPECS,
    SECONDARY_FUSEL_ROUTES,
    non_ehrlich_fraction_param,
    refund_carbon_to_sugar,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Re-exported: ``non_ehrlich_fraction_param`` moved to
#: :mod:`~fermentation.core.kinetics.carbon_routing` at D-111 (the fusel sourcing layer needs it
#: to derive a secondary route's primary-branch share, and cannot import this module — that would
#: be a cycle). Kept in ``__all__`` so every existing importer, including the tests, is unaffected.
__all__ = ["NON_EHRLICH_FRACTION_PARAMS", "PrecursorNonEhrlichFates", "non_ehrlich_fraction_param"]


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

    **D-106 changed what "the re-route's own draw" means, and the split survived it exactly.**
    Charging the Ehrlich decarboxylation CO₂ made that draw a full mole of precursor per alcohol
    instead of ``(n-1)/n``, so this Process now scales against alcohol carbon **+ CO₂ carbon**
    (via the shared :func:`~fermentation.core.kinetics.byproducts.ehrlich_co2_carbon`). Scaling
    against the alcohol alone would have realised a *lower* ``f`` than the file's sourced one —
    measured, threonine 0.82 → 0.774 — which is what the crux test caught.

    **The sourced fraction needed no recalibration, and the reason is structural**: ``f_i`` is a
    *ratio*, and D-106 scales the Ehrlich branch and this lump **equally** (the lump is defined off
    the Ehrlich draw), so the realised split is still exactly the sourced 77–86%. What D-106 does
    move is the **absolute** consumption, up ~12.6% — because the Ehrlich branch is anchored to
    *alcohol production* and this sink to the *split ratio*, so correcting the stoichiometry
    necessarily scales both rather than holding total consumption fixed and shifting the split.
    That is a modelling choice, not an accident: the alternative — hold consumption, move the split
    — would have silently overridden a sourced number with a stoichiometric correction.
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
        # D-111: the shared draw helper resolves every secondary route's share, so this Process
        # reads it too — the lump it books is scaled against those branches as well.
        *(route.share_param for route in SECONDARY_FUSEL_ROUTES),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # EXACTLY the branches the re-route will book, from the identical helper it uses — so the
        # ratio this Process imposes is against the precursor carbon actually consumed, and the two
        # can never drift apart (the D-33/D-99 shared-helper discipline, generalised at D-111).
        draws = ehrlich_draws(y, schema, params)
        if not draws:
            return d
        # EVERY branch of a precursor must reach this sum — the D-111 correctness pin. `f` is the
        # fraction of consumed precursor going anywhere except **this model's Ehrlich alcohols**,
        # plural since D-111: valine feeds isobutanol AND isoamyl alcohol. Drop one branch
        # (accumulate with `=` rather than `+=`, so the second silently wins) and valine's realised
        # non-Ehrlich share lands at **0.497 against the file's sourced 0.62** — a sourced number
        # overridden by an arithmetic slip that conservation cannot see, since both books still
        # close. D-106's catch one route further on: there the omission was the CO₂, here the second
        # branch. `test_the_realised_split_is_exactly_the_sourced_fraction` fails on it; measured.
        #
        # What is NOT a hazard, recorded because this comment first claimed it was: applying f/(1−f)
        # per branch instead of to the sum is a **no-op**. `f` is per-SPECIES and both valine
        # branches share the species, so f/(1−f) is one constant and `Σ(cᵢ)·k == Σ(cᵢ·k)` is
        # linearity — measured identical at 0.000e+00 relative difference. The mutation test caught
        # the false claim; the *code* was right and its stated reason was not, which is this
        # project's oldest lesson (D-96/D-102/D-108/D-109) landing on prose in the beat citing it.
        ehrlich_by_precursor: dict[str, float] = {}
        for draw in draws:
            ehrlich_by_precursor[draw.precursor.species] = (
                ehrlich_by_precursor.get(draw.precursor.species, 0.0) + draw.precursor_carbon
            )
        carbon = 0.0  # total precursor carbon re-sourced off sugar into the non-Ehrlich lump
        nitrogen = 0.0  # total nitrogen it carries, refunded to ammonium
        for species, ehrlich_carbon in ehrlich_by_precursor.items():
            if ehrlich_carbon <= 0.0:
                continue
            precursor = SPEC_BY_SPECIES[species]
            f = params[non_ehrlich_fraction_param(species)]
            # f is a fraction of the consumed precursor, so f ∈ [0, 1); f → 1 would demand an
            # infinite draw against a finite alcohol. The parameter file bounds every entry well
            # below this; the guard is here because a bad ENSEMBLE draw off the uncertainty band
            # is the realistic way it would ever be reached, and a silent inf would poison the
            # solver rather than fail.
            if not 0.0 <= f < 1.0:
                raise ValueError(
                    f"{non_ehrlich_fraction_param(species)}={f} outside [0, 1): "
                    "it is the fraction of consumed precursor NOT becoming an Ehrlich alcohol"
                )
            # `ehrlich_carbon` is the SUM over this precursor's branches of the whole-molecule
            # draw the re-route books — the alcohol's carbon **plus the decarboxylation CO₂
            # charged to the same precursor** (D-106; two CO₂ on the D-111 KIC route). The CO₂
            # term is not optional: the split is over *consumed precursor*, and since D-106 an
            # Ehrlich branch consumes a full mole per alcohol rather than (n-1)/n. Scaling against
            # alcohol carbon alone would silently realise a LOWER f than the file's sourced one
            # (threonine: 0.82 → 0.77).
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
