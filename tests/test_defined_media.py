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
3. **Scored on the media their own sources used, the four guards do not recover.** Propanol's
   floor closes 86 % of its gap and still misses; both Minebois legs get *materially worse*.
4. **What the sourcing did buy is a different, upstream finding.** Crépin's Data Set S1 also
   reports the *residual* nitrogen, and it is ~0.2 % of the must's YAN at end of fermentation.
   This model leaves **40.8 %** standing, at any duration. That divergence sits underneath every
   number in this file — biomass is the denominator of every de-novo share — and the residual
   propanol miss is smaller than the span of one unfixed commensurability defect in the
   availability gate. See the last two tests.

**Nothing here is tuned.** No parameter file and no ``src/`` file changed for this record; the one
knob turned below is turned in a probe, and its shipped value is left exactly where it was.

**The comparison is ANCHORED, not merely plausible.** Every share below is computed by
:func:`~tests.test_fusel_keto_acid_node.de_novo_share_of` and
:func:`~tests.test_fusel_catabolic_shape._amino_acid_share` — the shipped helpers, not a
re-implementation. Run on the D-109 fixture they reproduce D-245's published 0.7744 / 0.2256 /
0.0542 / 0.0947 to 4 dp, which is what licenses reading the columns against each other.
"""

from __future__ import annotations

import pytest

from fermentation.core.chemistry import MOLAR_MASS, NITROGEN_ATOMS, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acid_pools import (
    AMINO_ACID_SPECS,
    ASSIMILABLE_SPECS,
    GENERIC_POOL,
)
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import brix_to_sugar_gpl
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
    below is made of. Named, not repaired — the fix is a core change and is the owner's call.
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


_MUSTS = {
    # Crépin §Materials and Methods: EC1118, 28 °C, pH 3.5, "240 g glucose · liter−1".
    "crepin": (CREPIN_MUST_MM, 240.0, 28.0),
    # Minebois §2.2: T73 (and two non-cerevisiae), 24 °C, pH 3.3, "100 g L−1 of glucose,
    # 100 g L−1 of fructose".
    "minebois": (minebois_must_mm(), 200.0, 24.0),
}
_FERMENT_DAYS = 14.0


def commensurate_scenario(which: str, *, days: float = _FERMENT_DAYS) -> Scenario:
    """The named paper's own must, as a scenario.

    ``amino_acids_gpl`` is held POSITIVE and then overridden away pool by pool. That is not
    decoration: the compile seam disables the fusel re-route and the D-104 sink outright when the
    dose is ``<= 0`` (the D-32 isolability gate), so a pure-override fixture would silently switch
    off the Processes every test in this file measures and the guards would then be answering a
    question nobody asked. :func:`test_the_commensurate_musts_carry_the_papers_own_pools` asserts
    both halves — the Processes are live AND no pool kept the dose's spectrum share.
    """
    mm, sugar, celsius = _MUSTS[which]
    pools, yan = commensurate_pools(mm)
    initial = {"brix": _brix_for(sugar), "yan_mgl": yan, "pitch_gpl": 0.25} | pools
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


def _run(which: str, *, days: float = _FERMENT_DAYS, scale: dict[str, float] | None = None):
    cs = compile_scenario(commensurate_scenario(which, days=days))
    for name in _OTHER_PRECURSOR_CONSUMERS:
        cs.process_set.disable(name)  # KeyErrors on a rename rather than silently no-op
    params = cs.param_values
    for key, factor in (scale or {}).items():
        assert key in params, f"no such parameter {key!r}"
        params[key] = params[key] * factor
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
    """Assimilable N [mg N/L] standing in the ``N`` slot plus the eight pools, at ``index``."""
    total = float(traj.y[schema.slice("N"), index][0])
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


def test_the_de_novo_floor_is_still_missed_on_crepins_own_must(crepin_run):
    """The sourcing beat's headline: the commensurate must does NOT close propanol (D-246).

    D-244 §6 recorded the commensurability violation and declined to repair it; D-245 §7 named the
    paper's medium as the unlock. The medium is now in hand and **86 % of the gap was the must**:
    0.7744 on the D-109 fixture, 0.7963 here, against a floor of 0.80. The claim and the margin
    are on separate lines (D-245's own lesson), because the claim — *it still misses* — is what a
    future beat must not lose, while the exact 0.7963 is a pin that may legitimately move.
    """
    traj, schema, params, _ = crepin_run
    propanol = next(s for s in FUSEL_SPECS if s.pool == "propanol")
    share = de_novo_share_of(traj, schema, params, propanol)

    assert share < _SOURCED_DE_NOVO_FLOOR, (
        f"propanol reads {share:.4f} de novo on Crépin's OWN must and now CLEARS the sourced "
        f"{_SOURCED_DE_NOVO_FLOOR:.0%} floor. That is this record's central finding reversing, "
        "so the xfail on test_every_sourced_fusel_is_de_novo_dominated[propanol] and the reason "
        "strings citing D-246 both need re-deciding — not this assertion relaxing"
    )
    assert 0.790 <= share <= 0.800, (
        f"propanol's commensurate de-novo share left the [0.790, 0.800] pin: {share:.4f}"
    )


def test_the_low_band_is_still_missed_on_crepins_own_must(crepin_run):
    """The band restating the floor moves with it, and lands in the same place (D-246).

    ``max(amino-acid share) < 0.20`` is the 80 % floor written upside down, so it is not
    independent evidence — it is the same measurement and it is here to keep the pair moving
    together. 0.2256 on the D-109 fixture, 0.2037 here.
    """
    traj, schema, params, _ = crepin_run
    shares = {s.pool: 1.0 - de_novo_share_of(traj, schema, params, s) for s in FUSEL_SPECS}
    worst = max(shares, key=lambda k: shares[k])
    assert worst == "propanol", f"the band's worst alcohol is now {worst}, not propanol: {shares}"
    assert shares[worst] >= 0.20, f"the low band is now cleared on Crépin's must: {shares}"
    assert 0.200 <= shares[worst] <= 0.210


def test_both_minebois_legs_get_WORSE_on_her_own_must(minebois_run):
    """The two D-120 legs are not medium-limited — the commensurate must moves them away (D-246).

    D-245 §5 measured the model over-attributing to amino acids against Minebois's in-study
    shares: isoamyl 5.42 % vs her 5.34 %, isobutanol 9.47 % vs her 8.78 % — trips of 1.5 % and
    7.9 % relative, and it explicitly warned that the isoamyl one sits inside this harness's own
    cap-window systematic and "must not be built on".

    **Scored on her own must both roughly double, and that caveat no longer shields either.**
    Isoamyl reads ~73 % over her measurement and isobutanol ~68 % over: an order of magnitude
    outside the systematic D-245 invoked. So the over-attribution is a real property of the model
    at her nitrogen, not an artefact of scoring it on a must 25 % richer than hers.
    """
    traj, schema, params, _ = minebois_run
    for pool in ("isoamyl_alcohol", "isobutanol"):
        model = _amino_acid_share(traj, schema, params, pool)
        measured = _MINEBOIS_AMINO_ACID_SHARE[pool]
        assert model > measured, (
            f"{pool} no longer over-attributes on Minebois's own must ({model:.4f} vs her "
            f"{measured:.4f}) — D-120's direction leg would be back, which is a re-decision"
        )
        assert 1.55 <= model / measured <= 1.90, (
            f"{pool} reads {model / measured:.2f}x Minebois's {measured:.2f} on her own must; "
            "the pin is two-sided because a FALLING ratio re-opens D-120 and a rising one is a "
            "worsening neither D-245 nor this record predicted"
        )

    # 2-PE stays the green leg it became at D-245, on her must as on the fixture.
    two_pe = _amino_acid_share(traj, schema, params, "2_phenylethanol")
    assert 0.0 < two_pe < _MINEBOIS_AMINO_ACID_SHARE["2_phenylethanol"]


def test_the_model_leaves_two_fifths_of_crepins_nitrogen_standing(crepin_run):
    """What the sourcing actually bought — an upstream divergence, measured (decision D-246).

    Data Set S1 reports residuals as well as initials. Crépin's yeast consumes essentially **all**
    of the must's assimilable nitrogen: the EF block reads 0.0000 for ammonium, arginine,
    glutamine and every precursor, and the YAN column falls 9.636 -> ~0.02 mM. This model leaves
    **40.8 %** of it standing, almost all in the two identity-agnostic pools.

    **It is not slowness.** Held at 60 days instead of 14 the residual does not move in the fourth
    decimal — growth has stopped with the nitrogen still there, because assimilation here is
    growth-coupled and can consume only what growth demands. Real yeast over-accumulates well past
    immediate anabolic need, and Crépin's is exhausted at N_T (28 h), long before its EF at 150 h.

    This sits UNDERNEATH every other number in this file: biomass is the denominator of every
    de-novo share, and the model builds 2.14 g/L where consuming the whole must at its own
    ``biomass_N_fraction`` would build ~3.2.
    """
    traj, schema, _, _ = crepin_run
    initial = _assimilable_n_mgl(traj, schema, 0)
    left = _assimilable_n_mgl(traj, schema, -1)
    fraction = left / initial

    assert fraction > 10.0 * CREPIN_MEASURED_RESIDUAL_FRACTION, (
        f"the model now leaves {fraction:.1%} of Crépin's assimilable nitrogen against her "
        f"measured {CREPIN_MEASURED_RESIDUAL_FRACTION:.1%} — the divergence this record named "
        "has closed, which is a finding and not a reason to delete the guard"
    )
    assert 0.38 <= fraction <= 0.43, f"unconsumed nitrogen left [38 %, 43 %]: {fraction:.3f}"

    long_traj, long_schema, _, _ = _run("crepin", days=60.0)
    long_fraction = _assimilable_n_mgl(long_traj, long_schema, -1) / _assimilable_n_mgl(
        long_traj, long_schema, 0
    )
    assert long_fraction == pytest.approx(fraction, abs=1e-3), (
        f"at 60 days the residual is {long_fraction:.4f} against {fraction:.4f} at 14 — if time "
        "now consumes it, the cause is slow uptake rather than growth-coupled uptake and the "
        "mechanism this test names is wrong"
    )


def test_the_residual_propanol_miss_fits_inside_the_gates_own_commensurability_defect(crepin_run):
    """The 0.46 % that is left is smaller than one unfixed defect one layer down (D-246).

    :func:`~fermentation.core.kinetics.amino_acid_pools.depletion_gate` scales its half-saturation
    by ``Σ params[must_aa_fraction_*]`` — the **must-spectrum** shares, i.e. the composition a
    1.0 g/L dose would produce. A per-species override bypasses that spectrum, so on a real must
    the gate is scaled for a pool it is not gating: Crépin's two identity-agnostic pools hold
    0.5796 g/L where the spectrum implies 0.810.

    Rescaling the constant by that ratio — **derived from the two numbers, not fitted to an
    outcome** — carries propanol across the floor. The knob is turned in this probe only; its
    shipped value is untouched, and this test asserts the shipped run still misses so the two can
    never be confused.

    So the honest verdict on propanol is neither pass nor defect: on its source's own must the
    model lands ON the sourced floor to within the span of a commensurability defect that has not
    been repaired, with the sign unfavourable. Repairing that gate is a core change and the
    owner's call; it is named here so the next beat does not re-derive it.
    """
    _, _, shipped_params, cs = crepin_run
    spectrum = sum(cs.parameters[s.fraction_param].value for s in ASSIMILABLE_SPECS)
    pools, _ = commensurate_pools(CREPIN_MUST_MM)
    held = pools["arginine_gpl"] + pools[f"{GENERIC_POOL}_gpl"]
    ratio = held / spectrum
    assert 0.70 <= ratio <= 0.73, f"the gate's commensurate rescaling moved: {ratio:.4f}"

    propanol = next(s for s in FUSEL_SPECS if s.pool == "propanol")
    traj, schema, params, _ = _run("crepin", scale={"K_amino_acids": ratio})
    rescaled = de_novo_share_of(traj, schema, params, propanol)

    assert rescaled >= _SOURCED_DE_NOVO_FLOOR, (
        f"rescaling the gate to the pool it actually gates leaves propanol at {rescaled:.4f}, "
        f"still under the {_SOURCED_DE_NOVO_FLOOR:.0%} floor — this record's claim that the "
        "residual miss fits inside that defect would then be false"
    )
    assert 0.800 <= rescaled <= 0.810
    assert params["K_amino_acids"] < shipped_params["K_amino_acids"], (
        "the probe must scale K_amino_acids DOWN; if this stops holding the two runs are the "
        "same run and the comparison is vacuous"
    )
