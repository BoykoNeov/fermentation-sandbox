"""The speciated amino-acid pool registry — the single source of truth for D-100.

**What this closes (decision D-99 → D-100).** The D-99 fusel split raised the higher alcohols
~3.8× (an honest, forced rise) and broke three tests in unrelated subsystems — Maillard, Brett
growth, MLF growth — all through one channel: the D-33 Ehrlich re-route drained the single lumped
``amino_acids`` pool to ~0, starving every other consumer. A pre-D-99 emulation proved D-99 did
not *create* that pathology (the re-route already ate ~96.5% of the pool; sotolon passed only on
the old lump's under-production). The finding was structural: **two speciated-scale consumers
cannot share one lumped substrate.** This module is the fix — the lump becomes eight
single-molecule pools, and every consumer draws the molecule it actually eats.

**The eight pools (:data:`AMINO_ACID_SPECS`).** Six are *precursors* — each is the specific amino
acid whose skeleton becomes a specific product, so drawing it is a chemical claim, not a
bookkeeping stand-in::

    leucine       → isoamyl_alcohol (D-99)        → 3-methylbutanal   (thermal Strecker, D-87)
    isoleucine    → active_amyl_alcohol           → 2-methylbutanal
    valine        → isobutanol                    → 2-methylpropanal
    threonine     → propanol                      → sotolon (via alpha-ketobutyrate)
    phenylalanine → 2_phenylethanol               → phenylacetaldehyde (oxidative Strecker, D-75)
    methionine    → (no fusel)                    → methional (D-75/D-87) + methanethiol (D-45)

Two are *identity-agnostic*: ``arginine`` (kept in the original ``amino_acids`` slot — the D-32
representative, renaming it would touch every consumer twice for no fidelity gain) and
``amino_acids_generic`` (every assimilable amino acid without its own slot, lumped as glutamine).
These two are what the yeast swap (D-32), MLF growth (D-38), Brett growth (D-40) and Maillard
browning (D-89) draw: those consumers build *biomass* or *melanoidin* from amino-acid nitrogen and
genuinely do not care which molecule supplied it. **The re-route never touches them** — which is
the decoupling that resolves D-100's cross-subsystem starvation.

**The two draw idioms.** Precursor consumers call :func:`draw_precursor_carbon` (carbon out of one
named species, its nitrogen deaminated to ``N``). Identity-agnostic consumers call
:func:`draw_assimilable_nitrogen` (a nitrogen demand met from {arginine, generic}, split by the
nitrogen each pool holds). Both return what the caller needs to close its own ledger, and both are
here so no two consumers can drift apart on the arithmetic (the D-33/D-99 shared-helper
discipline).

**The relative-depletion gate — the one uniform rule (the D-100 gate choice).** Every consumer
gates on its own substrate with the *shared* ``K_amino_acids`` scaled by that substrate's
must-spectrum share::

    gate_i = aa_i / (K_amino_acids · f_i + aa_i)      f_i = must_aa_fraction_<i>

with ``f`` summed over the substrate for a multi-pool draw (so the generic gate uses
``f_arg + f_generic`` ≈ 0.81). **This introduces no new parameters** — it is derived from
``K_amino_acids`` (already sourced, already shared by the swap/re-route/D-75/D-87/D-89/D-45 gates)
and the D-100 must spectrum (recorded from literature before any wiring). Two properties make it
the honest choice:

  * **It reduces.** At must-spectrum composition ``aa_i = f_i·aa_total``, so
    ``gate_i = f_i·aa_total/(f_i·K + f_i·aa_total) = aa_total/(K + aa_total)`` — the same gate for
    **every species and every subset**, so the split does not silently move the dosed baseline.
    Pinned by ``test_every_per_species_gate_reduces_to_the_lumped_gate_at_spectrum_composition``.

    **Two claims, with different standing — do not conflate them.** That all subsets agree is
    *structural*: it falls out of the algebra for any fractions. That the common value equals the
    **pre-split lumped gate** is *contingent on ``Σf = 1``*, which the sourced spectrum happens to
    satisfy exactly today. Since a dose is apportioned as ``f_i·D/Σf`` (:func:`normalized
    <fermentation.scenario.compile._wine_amino_acids>`), a spectrum summing to ``F ≠ 1`` seeds
    ``D/F`` and gives ``(D/F)/(K + D/F)`` — not ``D/(K+D)``. So an **ensemble sampling the
    fractions' uncertainty bands shifts the baseline slightly** on nearly every draw (acceptable at
    this tier; the fractions are speculative), and a future re-source that breaks the sum would move
    it permanently. The reduction test pins *current* values and will catch that.
  * **It bites where the physics is.** When one species is preferentially drained — exactly the
    D-100 pathology — *its own* gate falls while the others' do not. The lumped gate could not
    express that: it read a pool 38% arginine and concluded leucine was abundant.

What it is **not** is a Michaelis constant. ``K_amino_acids`` never was one (it is an availability
proxy calibrated at pool scale); scaling it by the spectrum keeps it an availability proxy at
species scale — a *relative-depletion* measure, honest because it is labelled as one. Per-species
Michaelis constants would be eight unsourced numbers wearing the costume of fidelity — the D-98
trap, declined here for the same reason D-99 declined per-species activation energies.

**The emergent anabolic/catabolic split (the D-100 finding, SOURCED at D-103).** Real must carries
~30-60 mg/L leucine, but wine makes ~150-250 mg/L isoamyl alcohol: **most higher alcohol is NOT from
amino-acid catabolism** — it is synthesised de novo from sugar via the valine/leucine biosynthetic
pathway. The lumped model could not know this; it let the re-route draw a fixed fraction of fusel
carbon from a big arginine pool forever, which is nonsense (arginine does not make isoamyl alcohol).
With the split, leucine's own gate throttles the re-route as leucine depletes, and the remaining
fusel carbon stays on the sugar stand-in. **The anabolic/catabolic ratio is therefore emergent** — a
consequence of the must spectrum and the fusel demand, not a fitted fraction. That is a fidelity
gain the split delivers for free, and the mechanism by which the D-100 tripwire flips.

D-100 made that argument by hand and **it holds**: it reasoned "leucine supplies only ~7% of the
isoamyl carbon", the model gives **7.15%** at D-100's dose, and D-103 found the source D-100 never
had — Rollero *et al.* 2017 measures **2-8%** by U-13C leucine labelling. The two independent
isotope methods agree *qualitatively* (de novo dominates) and **disagree quantitatively**: Rollero's
">90% ... from the carbon central metabolism" implies <10% catabolic, while the [U-13C]-glucose
tracer's ">75% hexose derived" implies <25%. **Recorded as two bands, never averaged into one** —
this model sits over the first and inside the second.

**A SINK THIS MODULE DOES NOT HAVE (a real gap, deliberately not patched here).** The re-route is
each precursor's *only* consumer, so **100% of consumed leucine is attributed to Ehrlich**. Real
yeast incorporate much of it into **protein** — the anabolic sink. Every precursor here therefore
disappears by the one route the model owns, whatever reality's split. Building it needs a sourced
biomass amino-acid composition and touches the conservation ledger ⇒ its own beat, not a footnote.

**Isolability (prime directive #3).** Every speciated slot defaults to 0, so an undosed run has
every gate at exactly 0 and is byte-for-byte the validated core — the D-32 guarantee, unchanged.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.state import FloatArray, StateSchema

#: The representative species for the identity-agnostic assimilable pool (decision D-32).
#: Arginine is the dominant yeast-assimilable amino acid in grape must and is N-rich (mass
#: C:N ≈ 1.29 ≪ biomass ≈ 4.3), the property that keeps the carbon refund below growth's
#: demand for any ψ ≤ 1 (:mod:`~fermentation.core.kinetics.amino_acids`). It keeps the
#: original ``amino_acids`` slot name at D-100 (the slot IS the arginine pool): renaming it
#: would touch every consumer a second time and buy no fidelity.
AMINO_ACID_SPECIES = "arginine"

#: The state slot holding :data:`AMINO_ACID_SPECIES`.
ARGININE_POOL = "amino_acids"

#: The species standing in for :data:`GENERIC_POOL` — every assimilable amino acid without its own
#: slot (glu/ala/ser/asp/his/lys/gly/trp/tyr/cys/GABA), lumped as **glutamine**: the canonical
#: N-rich assimilable-nitrogen proxy. Mass C:N ≈ 2.14 — still far below biomass's ≈ 4.3, so the
#: {arginine, generic} blend keeps the D-32 no-sugar-creation guarantee structurally (a blend of
#: two species both under the biomass ratio is itself under it, for any split).
GENERIC_SPECIES = "glutamine"

#: The state slot holding the identity-agnostic remainder (decision D-100).
GENERIC_POOL = "amino_acids_generic"


@dataclass(frozen=True)
class AminoAcidSpec:
    """One speciated amino-acid pool: its state slot, its molecule, its must-spectrum share.

    ``pool`` is the state slot name and ``species`` the :mod:`~fermentation.core.chemistry` key
    that weights it in ``total_carbon``/``total_nitrogen``. They differ only for the two
    identity-agnostic pools (``amino_acids``→arginine, ``amino_acids_generic``→glutamine); every
    precursor slot is named for its own molecule, so the slot IS the species.

    ``fraction_param`` names the ``must_aa_fraction_*`` parameter that both splits the dose at the
    compile seam and scales this species' depletion gate (module docstring).
    """

    pool: str
    species: str
    fraction_param: str


#: The eight speciated assimilable amino-acid pools (decision D-100). Ordered arginine, generic,
#: then the six precursors — the identity-agnostic pair first, mirroring how the compile seam
#: splits the dose. ``proline`` is deliberately absent: it is not assimilated anaerobically (it is
#: excluded from YAN by definition), so the ~48%-of-must proline fraction never enters the model
#: and ``amino_acids_gpl`` is honestly an *assimilable* amino-acid dose.
AMINO_ACID_SPECS: tuple[AminoAcidSpec, ...] = (
    AminoAcidSpec(ARGININE_POOL, AMINO_ACID_SPECIES, "must_aa_fraction_arginine"),
    AminoAcidSpec(GENERIC_POOL, GENERIC_SPECIES, "must_aa_fraction_generic"),
    AminoAcidSpec("leucine", "leucine", "must_aa_fraction_leucine"),
    AminoAcidSpec("isoleucine", "isoleucine", "must_aa_fraction_isoleucine"),
    AminoAcidSpec("valine", "valine", "must_aa_fraction_valine"),
    AminoAcidSpec("threonine", "threonine", "must_aa_fraction_threonine"),
    AminoAcidSpec("phenylalanine", "phenylalanine", "must_aa_fraction_phenylalanine"),
    AminoAcidSpec("methionine", "methionine", "must_aa_fraction_methionine"),
)

#: The identity-agnostic subset — the ONLY pools the yeast swap (D-32), MLF growth (D-38), Brett
#: growth (D-40) and Maillard browning (D-89) draw. Keeping them a named pair (rather than "the
#: first two of :data:`AMINO_ACID_SPECS`") is what makes the D-100 decoupling explicit: the
#: Ehrlich re-route touches no member of this tuple, so no amount of fusel production can starve
#: bacterial growth again.
ASSIMILABLE_SPECS: tuple[AminoAcidSpec, ...] = AMINO_ACID_SPECS[:2]

#: Slot lookup by molecule, for the precursor consumers (which know their chemistry, not the
#: schema layout). Every precursor's slot is its own species name, so this is an identity map
#: except for the two identity-agnostic pools.
SPEC_BY_SPECIES: dict[str, AminoAcidSpec] = {spec.species: spec for spec in AMINO_ACID_SPECS}


def amino_acid_pool(y: FloatArray, schema: StateSchema, pool: str) -> float:
    """The clamped mass [g/L] in one amino-acid slot; 0.0 when the schema lacks it.

    The ``max(..., 0.0)`` absorbs solver undershoot, and the ``not in schema`` fallback makes
    every amino-acid consumer a hard no-op on beer (whose amino-acid pools are not tracked, D-32).
    """
    if pool not in schema:
        return 0.0
    return max(float(y[schema.slice(pool)][0]), 0.0)


def depletion_gate(
    y: FloatArray,
    schema: StateSchema,
    params: Mapping[str, float],
    specs: Sequence[AminoAcidSpec],
    k_param: str = "K_amino_acids",
) -> float:
    """The relative-depletion availability gate over ``specs`` — in [0, 1) (decision D-100).

    ``Σaa / (K·Σf + Σaa)``: the caller's half-saturation scaled by the substrate's must-spectrum
    share, so the gate measures depletion *relative to typical abundance* rather than against a
    pool-scale constant a single species could never reach. Reduces exactly to the old lumped
    ``aa/(K + aa)`` at must-spectrum composition, for any ``specs`` (module docstring).

    ``k_param`` defaults to the shared ``K_amino_acids`` — the constant the yeast swap, the
    re-route and the D-45/D-75/D-87/D-89 gates all read. The two bacterial growth Processes pass
    their **own** half-saturation instead (``K_aa_mlf``, ``K_aa_brett``): those were never the
    shared availability constant, and D-100 scales whatever constant a consumer already used
    rather than quietly re-pointing it at a different one.

    Returns 0.0 on an empty (or untracked) substrate — the isolability guarantee every caller
    leans on, and the guard that stops any draw from driving a pool negative.
    """
    total = sum(amino_acid_pool(y, schema, spec.pool) for spec in specs)
    if total <= 0.0:
        return 0.0
    k_scaled = params[k_param] * sum(params[spec.fraction_param] for spec in specs)
    return float(total / (k_scaled + total))


def draw_precursor_carbon(
    d: FloatArray,
    schema: StateSchema,
    species: str,
    carbon: float,
) -> float:
    """Debit ``species`` for the mass carrying ``carbon`` [g C/L/h]; return the nitrogen released.

    The per-precursor idiom (decision D-100): the amino-acid mass is sized so the carbon leaving
    the pool **equals** the carbon entering the caller's products, and the nitrogen that mass
    carried is returned for the caller to deaminate into ``N``. Carbon and nitrogen therefore
    close to machine precision at every call site, exactly as the D-45/D-75 lumped idiom did — but
    now at the *right molecule's* carbon fraction, so the draw is a chemical claim rather than a
    stand-in.

    Accumulates into ``d`` with ``+=``: several products can share one precursor (methional and
    3-methylbutanal both come off :class:`MaillardStrecker`'s single pass; threonine feeds both
    propanol and sotolon), and a plain ``=`` would silently drop all but the last.
    """
    spec = SPEC_BY_SPECIES[species]
    mass = carbon / carbon_mass_fraction(spec.species)  # g amino acid/L/h supplying that carbon
    d[schema.slice(spec.pool)] -= mass
    return float(mass * nitrogen_mass_fraction(spec.species))


def spectrum_carbon_per_nitrogen(params: Mapping[str, float]) -> float:
    """The must spectrum's mass C:N — g carbon per g nitrogen, if released at must composition.

    ``Σ(f_i·c_i) / Σ(f_i·n_i)`` ≈ **1.9** at the D-100 spectrum, sitting between arginine's 1.29
    and the C-rich precursors' 3.4–7.7. This is what :func:`release_spectrum_nitrogen` will
    realise, and :class:`~fermentation.core.kinetics.autolysis.YeastAutolysis` needs it *before*
    releasing in order to size its debris remainder.

    **The load-bearing ordering (decision D-34, preserved at D-100):** dead biomass's C:N is 4–11
    across Coleman's whole nitrogen range, so it stays above this ratio and the debris carbon
    ``r·(f_C − f_N·R)`` is **structurally non-negative** — no clamp, no C⁰ kink for the BDF solver.
    D-100 narrows that margin (1.29 → 1.9 against a floor of 4) but does not threaten it: the
    spectrum would have to reach C:N 4 — i.e. become mostly phenylalanine — to flip the sign.
    """
    numerator = sum(
        params[spec.fraction_param] * carbon_mass_fraction(spec.species)
        for spec in AMINO_ACID_SPECS
    )
    denominator = sum(
        params[spec.fraction_param] * nitrogen_mass_fraction(spec.species)
        for spec in AMINO_ACID_SPECS
    )
    return float(numerator / denominator)


def release_spectrum_nitrogen(
    d: FloatArray,
    schema: StateSchema,
    params: Mapping[str, float],
    nitrogen: float,
) -> float:
    """Deposit ``nitrogen`` [g N/L/h] across all eight pools at must-spectrum composition.

    The inverse of the draw idioms — the one *source* of amino acids in the model
    (:class:`~fermentation.core.kinetics.autolysis.YeastAutolysis`, D-34). Masses are apportioned
    so the released **mass** vector matches the must spectrum while the released **nitrogen** sums
    to exactly ``nitrogen``; the carbon it carries is returned so the caller can close its ledger
    (``nitrogen · spectrum_carbon_per_nitrogen(params)``, by construction).

    **This is what makes the refill matter at D-100.** Under the lump, autolysis released pure
    arginine — a molecule that feeds no Ehrlich alcohol and no Strecker aldehyde — so the pool it
    refilled could not restore the aroma routes that had drained. Releasing the spectrum puts
    leucine, methionine and threonine back, which is how a real *sur lie* wine regenerates the
    precursors fermentation consumed. Aging precursors are therefore now **dominantly
    autolysis-sourced**, and thermal/oxidative Strecker aroma is **strongly lees-dependent** — the
    published sur-lie mechanism, arriving as a consequence of speciation rather than a modelled
    rule.

    **The DIRECTION is sourced; the MAGNITUDE is not — do not read the extreme as a prediction.**
    The model currently says a no-lees aged wine makes *essentially zero* branched-chain Strecker
    aldehyde. **D-103 corrects what D-100 wrote here.** D-100 claimed the re-route's catabolic
    fraction was "a lumped estimate (~0.5 via the shared ``K_amino_acids`` gate)" against "a
    literature contribution nearer 20-50%". Both numbers were wrong. Measured exactly (the
    re-route is each precursor's only consumer, so the fraction is a state difference, not a
    quadrature): it is **0.192** at D-100's own dose, and **0.21-0.33** at a realistic must
    (leucine 30-60 mg/L) — never ~0.5, which needs ~5 g/L amino acids. And the "20-50%" was
    **uncited**; the sourced contribution is **lower**, so that band would have *acquitted* a
    model that should be convicted (see D-103 for why an uncited number can no more acquit than
    convict).

    **There is no single fraction to bound — the defect is this gate's SHAPE (decision D-103).**
    Rollero *et al.* 2017 (Microb. Biotechnol. 10:1649-1662, U-13C leucine/valine in synthetic
    must) measures the catabolic contribution as **uniformly low across species** — isoamyl
    alcohol **2-8%**, isobutanol **5-15%**, ">90% of the acids and higher alcohols ... derived
    from intermediates produced by the carbon central metabolism". This gate instead spans
    **~8% (isoamyl) to ~82-93% (propanol)** — an 11x spread set by each alcohol's *demand/supply*
    ratio, so it over-attributes precisely the MINOR alcohols whose precursor is abundant relative
    to their small carbon draw. ``K_amino_acids`` is one shared scalar (assimilation/MLF/Brett
    read it too) and **cannot** reshape a per-species spread ⇒ there was never a knob here.

    **What the gate DOES earn, and it is not nothing:** both Rollero trends are reproduced — the
    fraction *rises with nitrogen* ("increased with the initial nitrogen content") and *falls as
    precursors deplete* ("enrichments decreased throughout the fermentation process ... as the
    exogenous amino acids were depleted"). And on the compound carrying ~53% of the fusel carbon
    the model is close: isoamyl **8.0-13.4%** vs a measured **2-8%**. That is why the *aggregate*
    over-attribution looks mild (14-23.5% over Rollero's measured pair) — isoamyl's weight masks
    the per-species defect. **The aggregate is the misleading statistic here; the shape is the
    finding.**

    "Silent without lees" is **still not validated, and D-103 does not settle it either way**: it
    is an output of this same gate (phenylalanine exhausts *because* the overstated 2-PE draw
    takes all its carbon), so a correct lower draw might leave residual precursor. In reality
    phenylalanine does exhaust — but via protein synthesis, a sink this model does not have
    (below). Right outcome, wrong route, unresolved. The lees-enrichment direction stands alone.

    **DOCUMENTED LIMITATION (decision D-100):** the release uses the **must** spectrum, because a
    sourced yeast-**autolysate** amino-acid spectrum is not in hand. Autolysate is protein
    hydrolysate, not must: it is richer in the branched-chain and sulfur amino acids and much
    poorer in arginine, so this understates exactly the precursors (leucine, methionine) the
    aroma routes want. The error is conservative — it under-produces autolytic aroma rather than
    inventing it — and re-sourcing it is a one-parameter-file change, not a structural one.
    """
    ratio = spectrum_carbon_per_nitrogen(params)
    denominator = sum(
        params[spec.fraction_param] * nitrogen_mass_fraction(spec.species)
        for spec in AMINO_ACID_SPECS
    )
    # alpha scales the whole spectrum vector so its nitrogen totals exactly `nitrogen`; each pool
    # then gets alpha·f_i grams, preserving the spectrum's mass ratios by construction.
    alpha = nitrogen / denominator
    for spec in AMINO_ACID_SPECS:
        d[schema.slice(spec.pool)] += alpha * params[spec.fraction_param]
    return float(nitrogen * ratio)


def _assimilable_nitrogen_shares(y: FloatArray, schema: StateSchema) -> list[float]:
    """Each assimilable pool's share of the nitrogen currently held in {arginine, generic}.

    The one place the identity-agnostic split rule lives, so :func:`draw_assimilable_nitrogen` and
    :func:`assimilable_carbon_per_nitrogen` cannot disagree about it — a disagreement would break
    carbon closure in :class:`~fermentation.core.kinetics.aging.MaillardBrowning`, which must size
    its melanoidin from the ratio *before* it knows the draw. Returns an empty list on an empty
    substrate (every caller gates to 0 there anyway).
    """
    held = [
        amino_acid_pool(y, schema, spec.pool) * nitrogen_mass_fraction(spec.species)
        for spec in ASSIMILABLE_SPECS
    ]
    total = sum(held)
    if total <= 0.0:
        return []
    return [pool_n / total for pool_n in held]


def assimilable_carbon_per_nitrogen(y: FloatArray, schema: StateSchema) -> float:
    """The {arginine, generic} blend's mass C:N — g carbon per g nitrogen drawn (decision D-100).

    Exactly the ratio :func:`draw_assimilable_nitrogen` will realise (both read
    :func:`_assimilable_nitrogen_shares`), which is what
    :class:`~fermentation.core.kinetics.aging.MaillardBrowning` needs: it sizes the melanoidin it
    forms from this ratio and only then draws, so a mismatch would silently break its carbon
    ledger. Lies between arginine's ≈ 1.29 and glutamine's ≈ 2.14 for any pool composition, hence
    always below biomass's ≈ 4.3 and below melanoidin's ≈ 8 — the two orderings the D-32
    no-sugar-creation guarantee and the D-89 positive-denominator guarantee rest on.

    Returns 0.0 on an empty substrate (the callers' gates are 0 there).
    """
    shares = _assimilable_nitrogen_shares(y, schema)
    if not shares:
        return 0.0
    return float(
        sum(
            share * carbon_mass_fraction(spec.species) / nitrogen_mass_fraction(spec.species)
            for spec, share in zip(ASSIMILABLE_SPECS, shares, strict=True)
        )
    )


def draw_assimilable_nitrogen(
    d: FloatArray,
    y: FloatArray,
    schema: StateSchema,
    nitrogen: float,
) -> float:
    """Debit {arginine, generic} for ``nitrogen`` [g N/L/h]; return the carbon that mass carries.

    The identity-agnostic idiom (decision D-100): the yeast swap, MLF growth, Brett growth and
    Maillard browning all state their demand in **nitrogen** (``ρ = f_N·dX/y_N``) because they
    build biomass or melanoidin from amino-acid nitrogen without regard to which molecule supplied
    it. So the demand is split **in proportion to the nitrogen each pool holds** — the
    dimensionally natural rule for a nitrogen-anchored demand, and the one that keeps both pools
    draining to zero *together* (a mass-proportional split would empty the N-poorer pool first and
    could drive it negative while the shared gate still read "available").

    This is a modelling choice, not a sourced claim: real yeast strongly prefer some amino acids
    over others (nitrogen catabolite repression takes ammonium and glutamine early, arginine
    late). Modelling that preference needs per-species uptake constants that are not sourced — the
    D-98 trap — so non-preferential uptake is the honest v1, and the generic bucket exists
    precisely because these consumers make no identity claim.

    Returns the carbon [g C/L/h] the drawn mass carries, which the caller either refunds to sugar
    (D-32) or books into biomass (D-38/D-40). Because both species sit far below biomass's C:N
    (arginine ≈ 1.29, glutamine ≈ 2.14, biomass ≈ 4.3), any split of the demand carries less
    carbon than the biomass it funds — so the **no-sugar-creation guarantee is structural for the
    blend**, not just for arginine (decision D-32's load-bearing property, preserved).
    """
    shares = _assimilable_nitrogen_shares(y, schema)
    if not shares:
        return 0.0  # empty substrate — the caller's gate is 0 here too, so this is belt-and-braces
    carbon = 0.0
    for spec, share in zip(ASSIMILABLE_SPECS, shares, strict=True):
        mass = (share * nitrogen) / nitrogen_mass_fraction(spec.species)
        d[schema.slice(spec.pool)] -= mass
        carbon += mass * carbon_mass_fraction(spec.species)
    return float(carbon)
