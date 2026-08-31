"""The two papers' OWN defined media, sourced — and what they do to D-245's five xfails.

D-245 §7 left one item: the fusel guards are scored on a must carrying 2.25x Crépin's nitrogen,
and "whether 77.4 % against Crépin's 19 %-exogenous is a model defect or the expected reading of
a richer must cannot be settled without that paper's medium. Nothing here composes one." The
memory note sharpened it to **one sourcing ask gating four of the five remaining failures**.

**The ask was answerable, and answering it closes none of the four** (decision D-246).

1. **Both papers are ONE recipe.** Crépin and Minebois both use Bely, Sablayrolles & Barre 1990.
   Minebois prints the stock table verbatim (§2.2). Crépin does not tabulate it — it points at
   its supplemental **Data Set S1**, which is better than a recipe: *measured* initial
   concentrations, mM, all 20 species, mean of 14 independent fermentations. Species by species
   the two musts sit at a near-constant ratio (1.40–1.86, median 1.55), which is what "same
   recipe, different nitrogen" looks like.
2. **Neither paper's headline nitrogen number is this model's ``yan_mgl``**, and the difference is
   not rounding. See :func:`test_the_papers_own_nitrogen_convention_reproduces_their_own_number`.
   **D-248 measured that conflation and REFUSED the repair** — the carve-out and release frames
   cancel exactly, so the shipped seam is outcome-correct and fixing only the declaration would
   cost 6.1 % of it. Do not re-propose it as an open item.
3. **Scored on the media their own sources used, the four guards do not recover.** Propanol's
   floor closes 86 % of its gap and still misses; both Minebois legs get *materially worse*.
4. **What the sourcing did buy is a different, upstream finding.** Crépin's Data Set S1 also
   reports the *residual* nitrogen, and it is ~0.2 % of the must's YAN at end of fermentation.
   This model leaves **40.8 %** standing, at any duration. That divergence sits underneath every
   number in this file — biomass is the denominator of every de-novo share. See the last three
   tests.
5. **D-247 corrects §6 of that record.** D-246 read the last 0.46 % of propanol's miss as fitting
   inside the availability gate's own commensurability defect. The correction it *described* —
   re-referencing each pool's half-saturation from the must-spectrum share to the share the
   declared must really holds — was measured here and moved propanol **0.796275 → 0.796017**,
   6.9 % of the gap in the wrong direction. Its probe cleared the floor only because scaling the
   shared constant uniformly is a **level** change, bit-identical to scaling all eight spectrum
   shares by that same factor. The residual miss was therefore left unattributed, and the repair
   measured-and-refused rather than owed.
6. **D-248 builds item 4 and it moves every column in this file.** The 40.8 % was never a
   parameter: assimilation's only route into biomass nitrogen ran at ``ψ·gate·f_N·base_dx``,
   strictly below growth's own draw, so ammonium could only fall — and when it hit zero growth's
   Monod shut growth off and the swap, being proportional to the growth rate, stopped with it.
   :class:`~fermentation.core.kinetics.amino_acids.AssimilableNitrogenUptake` un-couples uptake
   from that rate. The must is now consumed to **0.62 %** against Crépin's 0.2 %, peak biomass
   goes from **61.6 %** to **98.4 %** of the Coleman yield the compile seam itself installs, and
   through that one denominator **propanol clears its floor (0.7963 → 0.8784), the band clears
   (0.2037 → 0.1216), and both Minebois legs land ON her published shares** (1.73× → 1.01×,
   1.68× → 0.98×). Every test below is rewritten to the new behaviour with the old number kept
   in its docstring; none was relaxed.
7. **And D-247's refusal is REINFORCED rather than overturned.** Both gate rescalings it measured
   are now inert — the uniform probe is worth +1.4e-5 and the composition correction +4.3e-6,
   against the +0.0065 and −0.00026 they were. Once the pool is drawn to the gate's asymptote
   whatever its scale, how the gate is scaled stops being what decides the residual.

**Nothing here is tuned.** The knobs turned below are turned in probes; D-248's own parameter is
shipped at the bound where transport stops being limiting, and its insensitivity across a 200×
sweep is pinned in ``tests/test_assimilable_nitrogen_uptake.py``, not here.

**The comparison is ANCHORED, not merely plausible.** Every share below is computed by
:func:`~tests.test_fusel_keto_acid_node.de_novo_share_of` and
:func:`~tests.test_fusel_catabolic_shape._amino_acid_share` — the shipped helpers, not a
re-implementation. Run on the D-109 fixture they reproduce D-245's published 0.7744 / 0.2256 /
0.0542 / 0.0947 to 4 dp, which is what licenses reading the columns against each other.
"""

from __future__ import annotations

from collections.abc import Mapping

import pytest

from fermentation.core.chemistry import MOLAR_MASS, NITROGEN_ATOMS, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acid_pools import (
    AMINO_ACID_SPECS,
    ARGININE_POOL,
    ASSIMILABLE_SPECS,
    GENERIC_POOL,
    AminoAcidSpec,
)
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import brix_to_sugar_gpl
from fermentation.units.convert import cells_per_ml_to_pitch_gpl
from tests.test_fusel_catabolic_shape import _MINEBOIS_AMINO_ACID_SHARE, _amino_acid_share
from tests.test_fusel_keto_acid_node import (
    _OTHER_PRECURSOR_CONSUMERS,
    _SOURCED_DE_NOVO_FLOOR,
    de_novo_share_of,
)

_M_N = 14.007

#: Molar mass [g/mol] and N atoms for the species this model does **not** track. The seven it does
#: track are read from :mod:`~fermentation.core.chemistry` instead, so a transcription here can
#: never quietly disagree with the constants the run itself uses.
_EXTRA_MOLAR_MASS = {
    "alanine": 89.09, "aspartate": 133.10, "cysteine": 121.16, "glutamate": 147.13,
    "glycine": 75.07, "histidine": 155.15, "lysine": 146.19, "proline": 115.13,
    "serine": 105.09, "tryptophan": 204.23, "tyrosine": 181.19, "ammonium": 18.04,
}  # fmt: skip
_EXTRA_N_ATOMS = {
    "alanine": 1, "aspartate": 1, "cysteine": 1, "glutamate": 1, "glycine": 1,
    "histidine": 3, "lysine": 2, "proline": 1, "serine": 1, "tryptophan": 2, "tyrosine": 1,
    "ammonium": 1,
}  # fmt: skip


def _mw(species: str) -> float:
    return MOLAR_MASS[species] if species in MOLAR_MASS else _EXTRA_MOLAR_MASS[species]


def _n_atoms(species: str) -> int:
    """N atoms per molecule — the TOTAL, which is the frame ``nitrogen_mass_fraction`` uses."""
    return NITROGEN_ATOMS[species] if species in NITROGEN_ATOMS else _EXTRA_N_ATOMS[species]


#: Crépin *et al.* 2017, Appl. Environ. Microbiol. 83:e02617-16, supplemental "Dataset S1: Raw
#: data: residual amino acids (mM)", row ``Mean`` / "Initial concentrations (from 14 replicates)".
#:
#: **Transcription hazard, recorded so a re-read does not repeat it.** The mean row is split across
#: three physical lines by the PDF's column layout, and its middle block (Gly…Pro) lands on the
#: line *labelled* ``SD``. Taking the row as printed would silently swap nine means for nine
#: standard deviations. Cross-checked against the 14 individual ``T0`` replicates above it
#: (e.g. leucine's 14 values mean to 0.1733 against the 0.173 transcribed here).
CREPIN_MUST_MM = {
    "alanine": 0.750, "arginine": 0.882, "aspartate": 0.127, "cysteine": 0.023,
    "glutamine": 1.397, "glutamate": 0.375, "glycine": 0.114, "histidine": 0.094,
    "isoleucine": 0.116, "leucine": 0.173, "lysine": 0.048, "methionine": 0.088,
    "ammonium": 4.040, "phenylalanine": 0.107, "proline": 2.392, "serine": 0.335,
    "threonine": 0.298, "tryptophan": 0.444, "tyrosine": 0.045, "valine": 0.179,
}  # fmt: skip

#: The paper's own ``Nass`` row at N_T, mg N/L — the number its Methods prints as "180 mg
#: nitrogen · liter⁻¹". It is **consumed** assimilable nitrogen, not the must's content.
CREPIN_NASS_MGN = 179.64
#: The paper's own ``YAN`` column: the molar sum of assimilable species, proline excluded [mM].
CREPIN_YAN_MOLAR_MM = 9.636
#: Data Set S1's ``EF`` block: every assimilable species reads 0.0000 but for Ala/Gly/Cys traces,
#: and the YAN column lands at ~0.02 mM of the 9.636 it started at.
CREPIN_MEASURED_RESIDUAL_FRACTION = 0.02 / 9.636

#: Minebois *et al.* 2025, Microb. Biotechnol. 18:e70087 §2.2, VERBATIM: "The composition of the
#: stock solution of amino acids was (in g L−1): …", dosed at 9.24 mL per litre of SM.
#:
#: **The stock as printed ALSO lists "ammonium chloride (46.0)" and the same sentence then adds
#: "0.325 g of ammonium chloride … per litre of SM".** Reading both as real puts the must at
#: 435.6 mg N/L against a stated 300, so they cannot both be the ammonium that was added. Only
#: the 0.325 g/L is used here, which lands the must at 290.2 in the paper's own frame against its
#: printed 300 — a 3 % miss attributable to the stock table's 2-significant-figure entries.
#: Recorded rather than reconciled: this is the paper's ambiguity, not a transcription error.
MINEBOIS_STOCK_GPL = {
    "tyrosine": 1.4, "tryptophan": 13.7, "isoleucine": 2.5, "aspartate": 3.4,
    "glutamate": 9.2, "arginine": 28.6, "leucine": 3.7, "threonine": 5.8,
    "glycine": 1.4, "glutamine": 38.6, "alanine": 11.1, "valine": 3.4,
    "methionine": 2.4, "phenylalanine": 2.9, "serine": 6.0, "histidine": 2.5,
    "lysine": 1.3, "cysteine": 1.0, "proline": 46.8,
}  # fmt: skip
MINEBOIS_STOCK_ML_PER_L = 9.24
MINEBOIS_NH4CL_GPL = 0.325
_M_NH4CL = 53.49

#: Arginine at **3** N, tryptophan at **1**, histidine at **1**; everything else at its total.
#: **Derived from Crépin's own numbers, not assumed**: its per-species consumed mg N/L at N_T
#: divided by the initial mM above recovers exactly these counts, for all 19 species, to 3 s.f.
#: Chemically it is the assimilable subset — the indole and imidazole nitrogens are not released.
_PAPER_N_ATOMS = {"arginine": 3, "tryptophan": 1, "histidine": 1}

#: The six pools named for their own molecule. Arginine and the generic lump are handled apart.
_NAMED_POOLS = ("leucine", "isoleucine", "valine", "threonine", "phenylalanine", "methionine")
#: The eleven assimilable species with no slot of their own, which lump into ``GENERIC_POOL``.
_LUMPED = (
    "alanine", "aspartate", "cysteine", "glutamine", "glutamate", "glycine",
    "histidine", "lysine", "serine", "tryptophan", "tyrosine",
)  # fmt: skip

_MUST_FERMENTABLE = 0.930  # ``must_fermentable_fraction``; asserted against the file below


def minebois_must_mm() -> dict[str, float]:
    """Minebois's must [mM] — the stock at 9.24 mL/L, plus the separately-added ammonium."""
    mm = {s: g * MINEBOIS_STOCK_ML_PER_L / _mw(s) for s, g in MINEBOIS_STOCK_GPL.items()}
    mm["ammonium"] = MINEBOIS_NH4CL_GPL * 1000.0 / _M_NH4CL
    return mm


def paper_frame_mgn(mm: dict[str, float]) -> float:
    """Assimilable N [mg N/L] in the PAPERS' frame — Arg 3, Trp 1, His 1, proline excluded."""
    return sum(
        v * _PAPER_N_ATOMS.get(s, _n_atoms(s)) * _M_N for s, v in mm.items() if s != "proline"
    )


def model_frame_mgn(mm: dict[str, float]) -> float:
    """Assimilable N [mg N/L] in THIS MODEL's frame — every N atom, proline excluded.

    ``nitrogen_mass_fraction`` weights each pool at its molecule's full formula, which is right
    for the conservation ledger (mass balance must count every atom) and wrong for "assimilable",
    since yeast releases neither tryptophan's indole nor histidine's imidazole nitrogen. The two
    meanings share one field, ``yan_mgl``, and that conflation is what the 201.11-vs-179.91 gap
    below is made of.

    **MEASURED AND REFUSED AT D-248 — no longer "named, not repaired".** The conflation is real
    and it is one species wide on this model's own registry (arginine, 4 N against 3 assimilable;
    tryptophan and histidine are not tracked as pools). But the compile seam carves the pools out
    of ``yan_mgl`` at their TOTAL nitrogen and every in-run deamination releases that same total,
    so the two errors are the same error with opposite signs and the nitrogen the run makes
    available equals the number declared, **exactly, for any dose**. Repairing the declaration
    alone leaves more ammonium behind and the unchanged release frame still delivers all of it:
    measured, +15.28 mg N/L at a 0.5 g/L dose (6.1 % of the declaration) and +30.55 at 1 g/L, all
    of it now reaching biomass because D-248's uncoupled uptake consumes the must. So the halves
    are one repair, the shipped seam is outcome-correct, and the complete version is priced rather
    than built. See ``tests/test_wine_nitrogen_budget.py``'s D-248 block.
    """
    return sum(v * _n_atoms(s) * _M_N for s, v in mm.items() if s != "proline")


def _brix_for(sugar_gpl: float) -> float:
    """The Brix a must must declare to load ``sugar_gpl`` of fermentable sugar."""
    lo, hi = 1.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if brix_to_sugar_gpl(mid) * _MUST_FERMENTABLE < sugar_gpl:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def commensurate_pools(mm: dict[str, float]) -> tuple[dict[str, float], float]:
    """(the eight per-species overrides [g/L], the ``yan_mgl`` that must be declared with them).

    **The generic bucket conserves NITROGEN, not mass**, and the choice is load-bearing enough to
    state rather than bury: eleven species with no slot of their own collapse into one glutamine
    pool, and glutamine is not their average molecule. Preserving their model-frame nitrogen keeps
    the ledger every guard here reads exact; preserving mass would over-state that nitrogen by
    ~20 %. Both were measured (see the record) and the miss below survives either, but the choice
    argued for on principle is also the one that flatters the model most — so the miss is not an
    artefact of a pessimistic reconstruction.

    ``yan_mgl`` is then the model-frame total, NOT the paper's headline: the pools carry real,
    assimilable nitrogen and the ``N`` slot must be left holding exactly the must's own ammonium.
    Declaring the paper's 180 instead would strand 21 mg N/L of it — the D-244 §6 error inverted.
    """
    pools = {f"{s}_gpl": mm[s] * _mw(s) / 1000.0 for s in _NAMED_POOLS}
    pools["arginine_gpl"] = mm["arginine"] * _mw("arginine") / 1000.0
    lump_n_mm = sum(mm[s] * _n_atoms(s) for s in _LUMPED)
    pools[f"{GENERIC_POOL}_gpl"] = (lump_n_mm / _n_atoms("glutamine")) * _mw("glutamine") / 1000.0
    return pools, model_frame_mgn(mm)


def _override_key(spec: AminoAcidSpec) -> str:
    """The ``initial`` key overriding one pool, by the compile seam's own rule (D-100)."""
    return "arginine_gpl" if spec.pool == ARGININE_POOL else f"{spec.pool}_gpl"


def commensurate_gate_scales(mm: dict[str, float], values: Mapping[str, float]) -> dict[str, float]:
    """Per-species multipliers on ``K·f_i`` that re-reference the gate to the declared must (D-247).

    :func:`~fermentation.core.kinetics.amino_acid_pools.depletion_gate` scales its
    half-saturation by ``f_i``, the **must-spectrum** share, so that at spectrum composition
    ``gate_i`` collapses to the pre-split lumped ``aa/(K + aa)`` (D-100's reduction property).
    A per-species override breaks the premise: the pool no longer holds ``D·f_i/Σf``, so the
    scale belongs to a composition the run does not have. This returns the multipliers that
    restore the property for the composition it *does* have —

        ``f_eff_i = (held_i / Σheld) · Σf``   ⇒   ``multiplier_i = f_eff_i / f_i``

    — i.e. the realised shares renormalised to the spectrum's own sum. **Σf is preserved on
    purpose**: it is what keeps the correction a statement about composition only. Under it every
    gate reads ``Σheld/(K·Σf + Σheld)`` at t=0, the lumped gate at the must's *real* total, so
    the run's nitrogen level still reaches the gate through the numerator — where it always did.

    Each multiplier is therefore ``(realised share)/(spectrum share)``, a pure ratio of two
    compositions, and the weighted set averages to 1 by construction. That is what distinguishes
    it from the probe D-246 ran: see
    :func:`test_the_uniform_rescaling_that_clears_propanol_is_a_LEVEL_change_not_a_composition_one`.
    """
    pools, _ = commensurate_pools(mm)
    held = {spec: pools[_override_key(spec)] for spec in AMINO_ACID_SPECS}
    fractions = {spec: values[spec.fraction_param] for spec in AMINO_ACID_SPECS}
    sum_held, sum_f = sum(held.values()), sum(fractions.values())
    return {
        spec.fraction_param: (held[spec] / sum_held) / (fractions[spec] / sum_f)
        for spec in AMINO_ACID_SPECS
    }


_MUSTS = {
    # Crépin §Materials and Methods: EC1118, 28 °C, pH 3.5, "240 g glucose · liter−1".
    "crepin": (CREPIN_MUST_MM, 240.0, 28.0),
    # Minebois §2.2: T73 (and two non-cerevisiae), 24 °C, pH 3.3, "100 g L−1 of glucose,
    # 100 g L−1 of fructose".
    "minebois": (minebois_must_mm(), 200.0, 24.0),
}
_FERMENT_DAYS = 14.0

#: The pitch both fixtures run at, and it is **sourced** (decision D-253). Minebois §Materials and
#: Methods — the same Bely, Sablayrolles & Barre 1990 medium and the same lab lineage as Crépin —
#: inoculates at 1 × 10⁶ cells mL⁻¹, converted through the engine's own D-219 crossing. So the
#: value is *directly* sourced for the Minebois fixture and inherited by lineage for Crépin, who
#: states no inoculum at all (``grep -c`` over her text = 0).
SOURCED_PITCH_GPL = cells_per_ml_to_pitch_gpl(1.0e6)

#: What these fixtures used to pitch: the wine benchmark's house default, 6.25× larger and
#: **unsourced by either paper**. Retained as the contrast arm of D-249/D-251/D-252's measurements,
#: never as a fixture default again — D-253 moved it. Do not "simplify" the two-pitch comparisons
#: in ``test_nitrogen_timing_attribution.py`` down to one by deleting this.
HOUSE_PITCH_GPL = 0.25


def commensurate_scenario(which: str, *, days: float = _FERMENT_DAYS) -> Scenario:
    """The named paper's own must, as a scenario.

    ``amino_acids_gpl`` is held POSITIVE and then overridden away pool by pool. That is not
    decoration: the compile seam disables the fusel re-route and the D-104 sink outright when the
    dose is ``<= 0`` (the D-32 isolability gate), so a pure-override fixture would silently switch
    off the Processes every test in this file measures and the guards would then be answering a
    question nobody asked. :func:`test_the_commensurate_musts_carry_the_papers_own_pools` asserts
    both halves — the Processes are live AND no pool kept the dose's spectrum share.

    **The pitch is :data:`SOURCED_PITCH_GPL`, not the house 0.25 it used to be (decision D-253).**
    Both papers run the same medium out of the same lab lineage and Minebois states the inoculum;
    the 0.25 was the wine benchmark's default, inherited by accident and unsourced by either
    paper. Moving it buys Crépin's nitrogen-exhaustion clock (29.9 h against her measured 28 h,
    from 17.6 h) and costs peak biomass (3.21 g/L against her 3.39, from 3.42). Every number in
    this file is measured at the new pitch; :func:`~tests.test_nitrogen_timing_attribution.
    test_the_fixtures_pitch_is_sourced_and_the_house_default_is_6_25x_larger` pins the move.
    """
    mm, sugar, celsius = _MUSTS[which]
    pools, yan = commensurate_pools(mm)
    initial = {"brix": _brix_for(sugar), "yan_mgl": yan, "pitch_gpl": SOURCED_PITCH_GPL} | pools
    initial["amino_acids_gpl"] = 1.0
    return Scenario(
        name=f"d246-{which}",
        medium="wine",
        initial=initial,
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=celsius),
            TemperaturePoint(day=days, celsius=18.0),
        ],
        interventions=[],
        duration_days=days,
    )


def _run(
    which: str,
    *,
    days: float = _FERMENT_DAYS,
    scale: dict[str, float] | None = None,
    override: dict[str, float] | None = None,
):
    """``scale`` MULTIPLIES the compiled value; ``override`` REPLACES it.

    Both exist because callers mean different things. A commensurability probe means "this
    constant, times that ratio" and must follow the shipped value wherever it goes. A capacity
    sweep means "r = 2.6", full stop — and if it is written as a multiplier it silently becomes
    a different sweep the moment the shipped value moves, which is exactly what D-253 would have
    done to ``test_assimilable_nitrogen_uptake``'s grid (every ``r = 10`` in its prose would have
    started meaning 26). Use ``override`` whenever the docstring names an absolute value.
    """
    cs = compile_scenario(commensurate_scenario(which, days=days))
    for name in _OTHER_PRECURSOR_CONSUMERS:
        cs.process_set.disable(name)  # KeyErrors on a rename rather than silently no-op
    params = cs.param_values
    for key, factor in (scale or {}).items():
        assert key in params, f"no such parameter {key!r}"
        params[key] = params[key] * factor
    for key, value in (override or {}).items():
        assert key in params, f"no such parameter {key!r}"
        params[key] = value
    traj = simulate_scheduled(
        cs.process_set,
        params,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
    )
    assert traj.success, traj.message
    return traj, cs.schema, params, cs


@pytest.fixture(scope="module")
def crepin_run():
    """Crépin's own must, other precursor consumers off — computed ONCE (the D-245 pattern)."""
    return _run("crepin")


@pytest.fixture(scope="module")
def minebois_run():
    return _run("minebois")


def _assimilable_n_mgl(traj, schema, index: int) -> float:
    """Assimilable N [mg N/L] the yeast has NOT yet built into biomass, at ``index``.

    The ``N`` slot, the eight pools, and — since D-250 — the intracellular store. The store is
    part of this total because that is what the quantity has always meant here: D-248's
    "40.8 % -> 0.62 % residual" was measured when uptake's surplus sat in ``N``, so dropping the
    store would silently redefine the number rather than migrate it.

    **D-250 makes a second, narrower reading possible for the first time**, and it is the one
    Crépin actually measures: assimilable N still in the MEDIUM, i.e. this total minus the store.
    The two were indistinguishable before D-250 because the model held them in one slot. See
    ``test_nitrogen_stored_intracellularly.py`` for that reading and its number; D-249's
    conclusion is unchanged under it.
    """
    total = float(traj.y[schema.slice("N"), index][0])
    if "stored_nitrogen" in schema:
        total += max(float(traj.y[schema.slice("stored_nitrogen"), index][0]), 0.0)
    for spec in AMINO_ACID_SPECS:
        if spec.pool in schema:
            pool = max(float(traj.y[schema.slice(spec.pool), index][0]), 0.0)
            total += pool * nitrogen_mass_fraction(spec.species)
    return total * 1000.0


def test_the_papers_own_nitrogen_convention_reproduces_their_own_number():
    """Both headline N figures are recovered — in a frame that is NOT this model's (D-246).

    The convention is **derived, not assumed**. Crépin's Data Set S1 prints consumed nitrogen per
    species in mg N/L alongside the initial mM; dividing one by the other recovers an integer N
    count for every species, and those counts are Arg **3**, Trp **1**, His **1**, everything else
    its formula's total. Under them the must reads 179.91 against the paper's printed 179.64 —
    reproduced to 0.15 %, which is what makes it the paper's frame rather than a fitted one.

    The same convention lands Minebois at 290.2 against her printed 300.

    **And this model reads the same must 12 % heavier**, because ``nitrogen_mass_fraction`` counts
    tryptophan's indole and histidine's imidazole nitrogen. Both readings are correct for what
    they are — one is assimilable N, the other is total N — and ``yan_mgl`` is asked to be both.
    **That conflation was measured and the repair REFUSED at D-248**, on the grounds that the
    carve-out and release frames cancel exactly and repairing only the first would make the
    outcome worse by 6.1 % of the declaration; see :func:`model_frame_mgn` and
    ``tests/test_wine_nitrogen_budget.py``. It is not an open item and it is not unmeasured.
    """
    crepin_paper = paper_frame_mgn(CREPIN_MUST_MM)
    assert crepin_paper == pytest.approx(CREPIN_NASS_MGN, rel=0.005), (
        f"Crépin's must reads {crepin_paper:.2f} mg N/L in the convention derived from its own "
        f"consumed-vs-initial table, against the {CREPIN_NASS_MGN} its Nass row prints — if this "
        "drifts, the derivation of that convention no longer holds and every column below is "
        "being read in a frame the paper does not use"
    )
    assert paper_frame_mgn(minebois_must_mm()) == pytest.approx(300.0, rel=0.05)

    # The molar sum the paper's own YAN column reports, as an independent transcription check.
    molar = sum(v for s, v in CREPIN_MUST_MM.items() if s != "proline")
    assert molar == pytest.approx(CREPIN_YAN_MOLAR_MM, abs=0.002)

    # And the gap that makes the paper's headline unusable as ``yan_mgl``.
    crepin_model = model_frame_mgn(CREPIN_MUST_MM)
    assert crepin_model == pytest.approx(201.11, abs=0.5)
    assert crepin_model - crepin_paper == pytest.approx(21.2, abs=1.0), (
        "the two frames differ by ~21 mg N/L on Crépin's must; typing the paper's 180 into "
        "yan_mgl strands that much real assimilable nitrogen, which is D-244 §6 inverted"
    )


def test_the_commensurate_musts_carry_the_papers_own_pools(crepin_run):
    """Anti-vacuity: the fixture is the paper's must AND the Processes it measures are live.

    Two ways this file could pass while measuring nothing, both closed here. The pools could still
    be carrying the ``amino_acids_gpl`` dose's must-spectrum shares (threonine would read ~67 mg/L
    rather than Crépin's 35.5); or the dose could have been zeroed to "clean up" the override, at
    which point the compile seam disables the re-route and the D-104 sink and every share below
    becomes a division by an alcohol nothing routed.
    """
    traj, schema, _, cs = crepin_run
    assert "fusel_amino_acid_reroute" in cs.process_set, (
        "the fusel re-route is not in the compiled set — the compile seam disables it when "
        "amino_acids_gpl <= 0, so this fixture would be measuring the core, not the model"
    )
    assert "precursor_non_ehrlich_fates" in cs.process_set

    for species in ("threonine", "leucine", "valine", "phenylalanine"):
        loaded = float(traj.y[schema.slice(species), 0][0]) * 1000.0
        expected = CREPIN_MUST_MM[species] * _mw(species)
        assert loaded == pytest.approx(expected, rel=1e-6), (
            f"{species} loaded at {loaded:.2f} mg/L, not Crépin's measured {expected:.2f}"
        )

    # The ``N`` slot holds exactly the must's own ammonium, which is the whole point of declaring
    # yan_mgl in the model's frame rather than the paper's.
    ammonium = float(traj.y[schema.slice("N"), 0][0]) * 1000.0
    assert ammonium == pytest.approx(CREPIN_MUST_MM["ammonium"] * _M_N, rel=1e-3)


def test_the_de_novo_floor_is_cleared_on_crepins_own_must_once_uptake_is_uncoupled(crepin_run):
    """Propanol clears the sourced floor — and the cause is the BIOMASS DENOMINATOR (D-248).

    **The history, because the claim reversed and a future reader must know on what.** D-244 §6
    recorded the commensurability violation and declined to repair it. D-245 §7 named the paper's
    medium as the unlock; D-246 sourced it and found 86 % of the gap was the must — 0.7744 on the
    D-109 fixture, **0.7963** here — and still a miss against 0.80. D-247 measured the
    availability gate's own rescaling and found it worth −6.9 % of what was left, the wrong way,
    leaving the residual miss *unattributed*.

    **D-248 attributes it, and it was never a fusel parameter.** Un-coupling assimilable-nitrogen
    uptake from growth demand carries propanol to **0.8784**. Nothing about propanol's chemistry
    moved: uptake draws no threonine (its ``touches`` forbids it, and
    ``test_uptake_touches_no_precursor_pool_in_a_driven_run`` drives the claim), so the only
    channel is the denominator. Consuming the whole must builds 98.4 % of Coleman's own yield
    instead of 61.6 %, which is more de-novo alcohol off the sugar route against a
    threonine-derived amount that was already capped by an exhausted pool.

    The claim and the margin are on separate lines (D-245's own lesson): *it now clears* is what a
    future beat must not lose; the exact 0.8784 is a pin that may legitimately move.
    """
    traj, schema, params, _ = crepin_run
    propanol = next(s for s in FUSEL_SPECS if s.pool == "propanol")
    share = de_novo_share_of(traj, schema, params, propanol)

    assert share >= _SOURCED_DE_NOVO_FLOOR, (
        f"propanol reads {share:.4f} de novo on Crépin's OWN must and no longer clears the "
        f"sourced {_SOURCED_DE_NOVO_FLOOR:.0%} floor. D-248's central finding is that "
        "un-coupling uptake carries it over; if that has reversed, the xfail on "
        "test_every_sourced_fusel_is_de_novo_dominated[propanol] comes back and this is a "
        "re-decision — not an assertion to relax"
    )
    assert 0.870 <= share <= 0.886, (
        f"propanol's commensurate de-novo share left the [0.870, 0.886] pin: {share:.4f}"
    )


def test_the_low_band_is_cleared_on_crepins_own_must_once_uptake_is_uncoupled(crepin_run):
    """The band restating the floor moves with it, and clears in the same place (D-248).

    ``max(amino-acid share) < 0.20`` is the 80 % floor written upside down, so it is **not
    independent evidence** — it is the same measurement, and it is here to keep the pair moving
    together rather than to corroborate. 0.2256 on the D-109 fixture, 0.2037 under D-246, and
    **0.1216** once uptake is un-coupled. That the two move together is the point; if they ever
    stop, one of them has acquired a second cause.
    """
    traj, schema, params, _ = crepin_run
    shares = {s.pool: 1.0 - de_novo_share_of(traj, schema, params, s) for s in FUSEL_SPECS}
    worst = max(shares, key=lambda k: shares[k])
    assert worst == "propanol", f"the band's worst alcohol is now {worst}, not propanol: {shares}"
    assert shares[worst] < 0.20, (
        f"the low band is missed again on Crépin's must: {shares} — this is the floor above "
        "written upside down, so it and that test must move together or one has a second cause"
    )
    assert 0.114 <= shares[worst] <= 0.130


def test_both_minebois_legs_land_on_her_own_measurement_once_uptake_is_uncoupled(minebois_run):
    """The two D-120 legs were a nitrogen defect all along, and they close (decision D-248).

    **The history.** D-245 §5 measured the model over-attributing to amino acids against
    Minebois's in-study shares — isoamyl 5.42 % vs her 5.34 %, isobutanol 9.47 % vs her 8.78 % —
    and warned the isoamyl trip sat inside this harness's cap-window systematic. D-246 scored both
    on her own must and they roughly doubled: **1.73×** and **1.68×** her measurement, an order of
    magnitude outside that systematic, which retired the shield and made the over-attribution a
    real property of the model at her nitrogen.

    **Un-coupling uptake removes it: 1.01× and 0.98×.** This is the beat's strongest external
    result and the discipline behind it matters more than the number. Nothing was fitted to
    Minebois — ``amino_acid_uptake_capacity_ratio`` is set at the bound where transport stops
    limiting and every outcome here is unmoved across a 200× sweep of it
    (``tests/test_assimilable_nitrogen_uptake.py``). Nothing touched her precursors either: uptake
    draws only the two identity-agnostic pools, so valine and leucine reach the Ehrlich route
    exactly as before and the whole move is the denominator.

    **The pin is deliberately two-sided and tight.** A rising ratio is the D-246 defect returning;
    a falling one puts the model UNDER her measurement, which re-opens D-120's direction leg from
    the other side and would mean the denominator has over-shot.

    **FLAG (decision D-254): this band CONTAINS an over-attribution, and that is not an accident
    of width.** ``_amino_acid_share`` credits each valine branch at its designed share, but
    ``ehrlich_draws`` truncates the secondary branch at ``headroom`` and the total Ehrlich draw
    on valine is pinned at ``1 − f``, so what the truncation takes off isoamyl isobutanol absorbs.
    Measured from the draws the run actually applied, the two legs here are **0.989× and 1.049×**
    rather than the 1.014× and 0.982× this test reads. Either way one alcohol sits above
    Minebois, which by D-120's own logic — a ``(1 − f_de_novo)`` cap is warranted iff the model
    over-attributes — re-opens the de-novo cap question on this must. It is NOT actionable,
    because the two estimators disagree about which alcohol, and because isobutanol has no
    sourced ``f_de_novo`` in this repo.

    **The band is deliberately left as it is.** It was set at D-248 on this estimator and moving
    it to the corrected one would re-pin a passing guard to make a flag look tidy. The estimator
    correction lives in ``tests/test_fusel_provenance_estimator.py``, with the direction of the
    bias pinned per alcohol; this docstring exists so nobody reads 0.90–1.12 as evidence that
    every alcohol is under her measurement. It is not.
    """
    traj, schema, params, _ = minebois_run
    for pool in ("isoamyl_alcohol", "isobutanol"):
        model = _amino_acid_share(traj, schema, params, pool)
        measured = _MINEBOIS_AMINO_ACID_SHARE[pool]
        assert 0.90 <= model / measured <= 1.12, (
            f"{pool} reads {model / measured:.3f}x Minebois's {measured:.4f} on her own must "
            f"({model:.5f}); D-248 measured 1.01x and 0.98x. Above the band is D-246's "
            "over-attribution returning; below it the biomass denominator has over-shot and "
            "D-120's direction leg re-opens from the other side. Either way a re-decision"
        )

    # 2-PE stays the green leg it became at D-245, on her must as on the fixture — and it is the
    # one leg un-coupling barely moves (0.01505 -> 0.01400), because its share is set by the
    # sourced f_de_novo cap rather than by precursor availability.
    two_pe = _amino_acid_share(traj, schema, params, "2_phenylethanol")
    assert 0.0 < two_pe < _MINEBOIS_AMINO_ACID_SHARE["2_phenylethanol"]


def test_the_model_now_consumes_crepins_must_where_it_used_to_leave_two_fifths(crepin_run):
    """The divergence the sourcing bought, and its repair (decisions D-246 → D-248).

    Data Set S1 reports residuals as well as initials. Crépin's yeast consumes essentially **all**
    of the must's assimilable nitrogen: the EF block reads 0.0000 for ammonium, arginine,
    glutamine and every precursor, and the YAN column falls 9.636 -> ~0.02 mM (0.2 % left), with
    exhaustion reached at N_T ≈ 28 h against an EF at 150 h.

    **This model used to leave 40.8 % standing, at any duration.** It was never slowness: held at
    30 and 60 days instead of 14 the residual did not move in the fourth decimal. Growth had
    stopped with the nitrogen still there, because the only route into biomass nitrogen ran at
    ``ψ·gate·f_N·base_dx``, strictly below growth's own draw — so ammonium could only fall, and
    when it reached zero growth's Monod shut growth off and the swap stopped with it.

    **D-248 un-couples uptake from that rate and the must is consumed: 0.62 %.** Still ~3× Crépin's
    own, and the direction is recorded rather than smoothed over. The residual that remains is set
    by ``K_amino_acids``'s asymptote and **not** by the new parameter — it reads 0.622 % at every
    capacity from 0.25 to 50, which is what makes the shipped value a bound rather than a fit
    (pinned in ``tests/test_assimilable_nitrogen_uptake.py``).

    This sits UNDERNEATH every other number in this file: biomass is the denominator of every
    de-novo share, and the model now builds 3.42 g/L where it built 2.14.
    """
    traj, schema, _, _ = crepin_run
    initial = _assimilable_n_mgl(traj, schema, 0)
    left = _assimilable_n_mgl(traj, schema, -1)
    fraction = left / initial

    assert fraction < 0.05, (
        f"the model leaves {fraction:.1%} of Crépin's assimilable nitrogen standing. D-248's "
        "un-coupled uptake is what consumes it, so this is that repair having stopped working — "
        "and every column in this file rests on it"
    )
    assert fraction > CREPIN_MEASURED_RESIDUAL_FRACTION, (
        f"the model now consumes MORE completely ({fraction:.4%}) than Crépin measures "
        f"({CREPIN_MEASURED_RESIDUAL_FRACTION:.4%}); the direction D-248 recorded has flipped"
    )
    assert 0.004 <= fraction <= 0.010, f"unconsumed nitrogen left [0.4 %, 1.0 %]: {fraction:.5f}"

    # The duration control is KEPT and its meaning inverted. Under D-246 it proved the residual
    # was not slowness; here it proves the repair is not slowness either — the must is consumed
    # inside the 14 days, not merely on its way there by day 60.
    long_traj, long_schema, _, _ = _run("crepin", days=60.0)
    long_fraction = _assimilable_n_mgl(long_traj, long_schema, -1) / _assimilable_n_mgl(
        long_traj, long_schema, 0
    )
    assert long_fraction == pytest.approx(fraction, abs=1e-3), (
        f"at 60 days the residual is {long_fraction:.4f} against {fraction:.4f} at 14 — the "
        "consumption is supposed to be complete well inside the window, so a duration dependence "
        "here means uptake has become rate-limited and the 'bound not midpoint' reading is gone"
    )


def test_the_uniform_rescaling_that_clears_propanol_is_a_LEVEL_change_not_a_composition_one(
    crepin_run,
):
    """The probe that carries propanol over the floor moves the LEVEL, not the composition (D-247).

    **Kept, renamed, and its claim corrected.** D-246 §6 ran this probe and read the result as
    "the residual propanol miss fits inside the gate's own commensurability defect". The
    arithmetic below is unchanged and the number still stands — 0.7963 shipped, 0.8028 rescaled —
    but the reading was wrong, and the assertion that now closes this test is what shows it:
    scaling ``K_amino_acids`` by 0.7155 is **indistinguishable from scaling all eight**
    ``must_aa_fraction_*`` **shares by 0.7155**, because
    :func:`~fermentation.core.kinetics.amino_acid_pools.depletion_gate` only ever reads the
    product ``K·f_i``. A uniform scaling of every share changes no share *relative to any other*.
    It is a statement that the availability constant is too big, which nothing here sources —
    not a statement about Crépin's must having a different composition from a typical one.

    The ratio it uses (0.5796 held against the 0.810 a 1.0 g/L dose at spectrum composition would
    seed) compares a real pool mass against a *declared dose* the fixture only carries to keep the
    D-32 isolability gate open — so it folds the run's nitrogen level into what was presented as a
    composition correction. The composition-only correction is measured in the next test and is
    worth −0.0003. **What survives of D-246 §6 is the number, not the attribution.**

    The two runs are compared at ``rel=1e-12`` rather than bit-for-bit: they differ only in
    whether the solver sees ``(K·a)·f_i`` or ``K·(f_i·a)``, and float multiplication is not
    associative, so a bit-for-bit pin would be a platform pin — exactly D-238's scar.
    """
    _, _, shipped_params, cs = crepin_run
    spectrum = sum(cs.parameters[s.fraction_param].value for s in ASSIMILABLE_SPECS)
    pools, _ = commensurate_pools(CREPIN_MUST_MM)
    held = pools["arginine_gpl"] + pools[f"{GENERIC_POOL}_gpl"]
    ratio = held / spectrum
    assert 0.70 <= ratio <= 0.73, f"the gate's commensurate rescaling moved: {ratio:.4f}"

    propanol = next(s for s in FUSEL_SPECS if s.pool == "propanol")
    shipped_traj, shipped_schema = crepin_run[0], crepin_run[1]
    shipped = de_novo_share_of(shipped_traj, shipped_schema, shipped_params, propanol)
    traj, schema, params, _ = _run("crepin", scale={"K_amino_acids": ratio})
    rescaled = de_novo_share_of(traj, schema, params, propanol)

    # D-248 UPDATE, and it reinforces this record rather than overturning it. Under D-246 this
    # probe was worth +0.0065 (it carried propanol across the floor and was read as the
    # attribution). Once uptake is un-coupled the pools are drawn down to the gate's asymptote
    # whatever the gate's scale, so the probe is worth **+1.4e-5** — inert. The level was never
    # the mechanism; what the gate's scale bought was a slightly lower asymptote on a residual
    # that is now 0.6 % rather than 40.8 %.
    assert rescaled == pytest.approx(shipped, abs=5e-4), (
        f"the uniform rescaling moves propanol {shipped:.6f} -> {rescaled:.6f}, which is no "
        "longer the null D-248 measured. If the gate's LEVEL has regained purchase, the residual "
        "it acts on has grown back and D-248's uptake repair is the thing to check first"
    )
    assert rescaled > shipped, (
        "the probe no longer points the way it did — a lower K still means a lower asymptote and "
        "so slightly MORE nitrogen consumed; a reversal here is a re-decision"
    )
    assert params["K_amino_acids"] < shipped_params["K_amino_acids"], (
        "the probe must scale K_amino_acids DOWN; if this stops holding the two runs are the "
        "same run and the comparison is vacuous"
    )

    # The correction of D-246 §6, as an identity rather than an argument: the same run comes back
    # from scaling every spectrum share instead of the constant, so no share moved against any
    # other and nothing about composition was tested.
    uniform = {spec.fraction_param: ratio for spec in AMINO_ACID_SPECS}
    by_shares, share_schema, share_params, _ = _run("crepin", scale=uniform)
    assert de_novo_share_of(by_shares, share_schema, share_params, propanol) == pytest.approx(
        rescaled, rel=1e-12
    ), (
        "scaling all eight must_aa_fraction_* by the same factor no longer reproduces the "
        "K_amino_acids probe — if the gate has stopped reading only the product K·f_i, D-247's "
        "reason for calling this probe a level change needs re-deriving"
    )


def test_the_gates_OWN_commensurate_rescaling_leaves_propanol_where_it_was(crepin_run):
    """Re-referencing the gate to the must actually declared is worth −0.0003 (D-247).

    :func:`commensurate_gate_scales` is the correction D-246 §6 *described* — each pool's
    half-saturation re-scaled from the spectrum's share to the share the declared must really
    holds, with ``Σf`` preserved so nothing about the level moves. It is derived from two
    compositions and fitted to nothing.

    **It does not close propanol and it does not even point that way.** 0.796275 shipped,
    0.796017 rescaled: 6.9 % of the remaining gap to the floor, in the wrong direction, against
    D-246 §6's "it spans the whole of what is left of propanol's miss". The move is stable to six
    decimals across rtol/atol from 1e-6/1e-9 to 1e-10/1e-13, so it is a real (tiny) result and not
    solver noise; measured with every Process live it moves no fusel by more than 0.5 %.

    **So the defect is real and the repair is refused** — measured, not inherited, the D-120 /
    D-189 / D-195 pattern. What the gate scales by is an unsourced modelling device either way:
    D-100 justified the spectrum scaling as dynamic range with *zero new parameters* and
    explicitly declined per-species Michaelis constants as "eight unsourced numbers wearing the
    costume of fidelity" (the D-98 trap). Re-referencing to the run's own must does not source
    anything either — it swaps one unsourced reference for another, and buys nothing measurable.
    The reference is the open question; **only per-species half-saturations from literature can
    settle it**, and that is a sourcing ask, not a core change.

    The three multipliers pinned below are the shape of Crépin's must against the shipped
    spectrum, and they are why the effect cancels: arginine is poorer than typical, the generic
    remainder richer, and the tiny methionine pool nearly 4x its spectrum share.
    """
    traj_shipped, schema_shipped, shipped_params, _ = crepin_run
    propanol = next(s for s in FUSEL_SPECS if s.pool == "propanol")
    shipped = de_novo_share_of(traj_shipped, schema_shipped, shipped_params, propanol)
    scales = commensurate_gate_scales(CREPIN_MUST_MM, shipped_params)

    assert scales["must_aa_fraction_arginine"] == pytest.approx(0.574, abs=0.01)
    assert scales["must_aa_fraction_generic"] == pytest.approx(1.406, abs=0.01)
    assert scales["must_aa_fraction_methionine"] == pytest.approx(3.726, abs=0.05)
    assert scales["must_aa_fraction_arginine"] < 1.0 < scales["must_aa_fraction_generic"], (
        "the two identity-agnostic pools no longer move in opposite directions — the "
        "reconstruction's shape against the spectrum has changed and every number here with it"
    )

    traj, schema, params, _ = _run("crepin", scale=scales)
    rescaled = de_novo_share_of(traj, schema, params, propanol)

    # D-248 UPDATE: the correction has gone from tiny-and-wrong-way to INERT. It was worth
    # −0.00026 (−6.9 % of the gap to the floor) when 40.8 % of the must's nitrogen was still
    # standing; with uptake un-coupled it is worth **+4.3e-6**. The refusal therefore stands on
    # firmer ground than when it was made: the repair D-246 §6 described buys nothing measurable
    # at all now, and the reference question it turns on (D-247 §4) is still a sourcing ask.
    assert rescaled == pytest.approx(shipped, abs=5e-5), (
        f"the composition-only rescaling moves propanol {shipped:.6f} -> {rescaled:.6f}, which "
        "is no longer the null D-248 measured. D-247 refused this repair on a measured payoff of "
        "nil; if it has regained purchase that refusal needs re-deciding, not this pin relaxing"
    )
    # The "share of the gap to the floor" statistic D-247 reported is RETIRED rather than
    # re-pinned: propanol now sits ABOVE the floor, so that denominator is negative and the ratio
    # is not a quantity any more. What replaces it is the absolute move, which is the thing the
    # refusal actually rests on.
    assert abs(rescaled - shipped) < 5e-5, (
        f"the composition correction moves propanol by {rescaled - shipped:+.2e}; D-247 measured "
        "-2.6e-4 against a 40.8 % residual and D-248 measures it inert against a 0.6 % one"
    )
    assert 0.870 <= rescaled <= 0.886
