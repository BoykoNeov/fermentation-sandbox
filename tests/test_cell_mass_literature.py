"""What the literature actually offers for a yeast cell's dry mass (decision D-271).

D-219 left one item open and D-232 restated it unsupplied: **a source pairing a yeast CELL
COUNT with a GRAVIMETRIC DRY WEIGHT in one ferment at one timepoint.** Divide the one by the
other and D-230's branch 1 — Tyrell's counted crop read as a per-cell dry mass of
``70.9-91.9 pg`` — is settled or refuted outright.

D-271 hunted it and it does not exist as published. The negative is not a count of files
searched: every source reports the ONE currency its own endpoint needs, and never both. The
four instances are in :data:`_PAIRING_CANDIDATES` with the reason each fails, two of them
papers this repo already reads for other purposes.

What the literature offers *instead* is never a measurement of a cell's mass. It is either a
**geometric estimate** (volume x density x dry fraction) or a **rule of thumb**. Both classes
are represented here, and the finding that matters is the ordering: **every figure in reach
lands BELOW branch 1's demand**, including the engine's own settled gram and D-270 §7's
re-pricing of the elemental estimate. That is a direction, not a settlement — none of them is
a pairing, and none is an ale cell in wort.

The size route cannot rescue it either, and that is pinned rather than left for a future beat
to re-discover: the two sourced size series **disagree by more than a factor of four in volume
on the same laboratory strains**, so branch 1's demanded cell reads as an ordinary haploid in
one frame and as three times a diploid in the other. A route whose answer is chosen by which
source you open cannot adjudicate anything.

One thing does ship: D-219's geometric cross-check used ``rho ~ 1.11`` and an **unsourced**
"~30 % dry matter". Klis, de Koster & Brul 2014 sources both constants, and its 0.34 is the
value D-219's own printed upper edge already implies.

**Nothing here moves the engine's gram.** 4e-11 g is D-219's *definition* of the unit Coleman's
fit is counted in, not an estimate competing with these; a literature pg/cell cannot move it and
none is scored as though it could.

Receipts: ``docs/receipts/d271-cell-mass-pairing-hunt/`` -- ``PREREGISTER.md`` carries the
admissibility rule this beat fixed BEFORE looking, which is why a figure landing inside the
demand could not be read as support for it.
"""

from __future__ import annotations

import math
from typing import TypedDict

import pytest

from fermentation.units import cells_per_ml_to_pitch_gpl
from tests.test_biomass_nitrogen_frame import _D230_BRANCH_ONE_PG_PER_CELL


# ---------------------------------------------------------------------------
# The two-sided negative, as data
# ---------------------------------------------------------------------------
#
#: Every source reachable from this repo that could have carried the pairing, and what it
#: reports instead. ``counts`` and ``weighs`` are about the SAME ferment at the SAME timepoint —
#: which is the whole difficulty, and why Varela is False on both despite doing each once.
#:
#: This table is the machine-readable form of D-271's negative. A future beat that finds a
#: genuine pairing adds a row with both flags True and
#: :func:`test_no_source_in_reach_pairs_a_count_with_a_weighing` turns RED, which forces the
#: finding into a record instead of letting it land silently.
class _Candidate(TypedDict):
    """One source that could have carried the pairing, and why it does not."""

    counts: bool
    weighs: bool
    why: str


_PAIRING_CANDIDATES: dict[str, _Candidate] = {
    "Crepin 2017 Data Set S1": {
        "counts": False,
        "weighs": True,
        "why": (
            "gravimetric dry weight (filtered, washed, 105 C / 48 h to constant weight) with "
            "no cell count anywhere -- the closest anyone gets, and on disk since D-246"
        ),
    },
    "Varela 2004 (AEM 70:3392)": {
        "counts": False,
        "weighs": True,
        "why": (
            "counts (Neubauer) only the 1e6 cells/mL INOCULUM and weighs only the crop, so the "
            "count and the weighing are at different timepoints -- not a pairing"
        ),
    },
    "Foster 2022 (Front. Microbiol. 13:747546)": {
        "counts": True,
        "weighs": False,
        "why": (
            "counts cells by haemocytometer; dry cell weight appears ONLY as a trehalose "
            "normaliser (mg/g DCW), never at a fermentation timepoint -- already on record at "
            "D-219 §4"
        ),
    },
    "Tyrell 2013 (BrewingScience 66:76-83)": {
        "counts": True,
        "weighs": False,
        "why": (
            "the counted crop branch 1 is read off; pitched by count (MEBAK III 10.4.3) and "
            "prints no pitch by weight -- already on record at D-219 §4"
        ),
    },
}

# ---------------------------------------------------------------------------
# Klis, de Koster & Brul 2014, Eukaryot. Cell 13(1):2-9 (PMC3910951)
# ---------------------------------------------------------------------------
#
#: The stated route, verbatim: *"The biomass (dry weight) was calculated by multiplying the
#: volume with the density (1.11) to obtain the biomass (wet weight) of the cell and multiplying
#: the obtained value with the dry weight fraction (0.34) of the wet weight."*
#:
#: A COMPUTED figure, not a weighing — which is exactly why it is here and not treated as the
#: pairing. Its value to this repo is that it **sources the two constants D-219's own geometric
#: cross-check used**, one of which D-219 carried unsourced.
_KLIS_DENSITY_G_PER_ML = 1.11
_KLIS_DRY_FRACTION = 0.34

#: Klis Table 1, exponential-phase cells in rich medium at 30 C: volume [fL] and the dry mass
#: [pg] the paper prints for it. Both are transcribed so
#: :func:`test_klis_route_reproduces_its_own_printed_masses` can check the transcription against
#: the route rather than trusting either alone.
_KLIS_CELLS = {
    "haploid": {"volume_fl": 44.0, "printed_pg": 16.5},
    "diploid": {"volume_fl": 83.0, "printed_pg": 31.2},
}

#: BNID 101795, from *Physical Biology of the Cell* Table 1.1 — a stated rule of thumb
#: ("estimates ... carried out using a stick in the sand"), with no method and no condition.
#: Carried because it is the only figure in reach ABOVE the engine's gram, so excluding it
#: would flatter :func:`test_every_figure_in_reach_lands_below_branch_one`.
_RULE_OF_THUMB_PG = 60.0

#: Okada et al. 2023, Sci. Rep. 13:1722 (PMC9883461) Table 2: **apparent** diameters [um] from
#: flow-cytometry forward scatter calibrated on size beads. A different measurement frame from
#: Klis's microscopy volumes, which is the entire point of holding both.
_OKADA_APPARENT_DIAMETER_UM = {
    "haploid BY4741L": 7.3,
    "diploid BY4743L": 9.4,
    "brewing diploid K7A": 12.6,
}

# ---------------------------------------------------------------------------
# D-219's own numbers, and D-270's re-pricing
# ---------------------------------------------------------------------------
#
#: D-219's geometric cross-check: its assumed volume range for a wine/ale cell, and the band it
#: printed off an unsourced "~30 % dry matter".
_D219_ASSUMED_VOLUME_FL = (100.0, 150.0)
_D219_PRINTED_CROSSCHECK_PG = (30.0, 57.0)

#: D-219's settled per-cell dry mass and the band the elemental route gives across
#: ``biomass_N_fraction``'s own 0.08-0.14 uncertainty. **Not this module's to move.**
_D219_SETTLED_PG = 40.0
_D219_SETTLED_BAND_PG = (28.0, 50.0)

#: D-270 §7's re-pricing of the ELEMENTAL estimate onto the sourced composition. It is the
#: engine-side estimate, **not** branch 1's demand: D-270 narrowed the GAP between the two
#: (2.03-2.64x to 1.14-1.93x) and left the demand where D-230 put it. Pinned here because the
#: two were conflated once already, in D-271's own first draft.
_D270_ELEMENTAL_REPRICED_PG = (47.71, 62.12)


def _volume_fl_to_dry_pg(volume_fl: float, dry_fraction: float = _KLIS_DRY_FRACTION) -> float:
    """Cell volume [fL] -> dry mass [pg] by Klis's stated route (1 fL = 1e-12 mL)."""
    return volume_fl * 1e-12 * _KLIS_DENSITY_G_PER_ML * dry_fraction * 1e12


def _dry_pg_to_volume_fl(dry_pg: float, dry_fraction: float = _KLIS_DRY_FRACTION) -> float:
    """Inverse of :func:`_volume_fl_to_dry_pg`."""
    return dry_pg * 1e-12 / dry_fraction / _KLIS_DENSITY_G_PER_ML * 1e12


def _sphere_volume_fl(diameter_um: float) -> float:
    """Equivalent-sphere volume [fL] for a diameter [um] (1 fL == 1 um^3)."""
    return math.pi / 6.0 * diameter_um**3


# ---------------------------------------------------------------------------
# The negative itself
# ---------------------------------------------------------------------------


def test_no_source_in_reach_pairs_a_count_with_a_weighing() -> None:
    """D-271's finding: nobody reports both currencies for one ferment at one timepoint.

    The forcing device described on :data:`_PAIRING_CANDIDATES`: this is the guard a future
    beat trips by finding the pairing, so the find cannot land without a record.
    """
    paired = [
        name for name, entry in _PAIRING_CANDIDATES.items() if entry["counts"] and entry["weighs"]
    ]
    assert paired == [], (
        f"a source now pairs a count with a weighing: {paired}. That settles or refutes D-230 "
        "branch 1 outright and D-219's open item with it -- append a record, do not delete "
        "this guard"
    )
    # Both halves exist in the corpus; it is only their conjunction that does not.
    assert any(entry["weighs"] for entry in _PAIRING_CANDIDATES.values())
    assert any(entry["counts"] for entry in _PAIRING_CANDIDATES.values())


def test_every_candidate_carries_the_reason_it_fails() -> None:
    """A registry of names would let a later reader assume the wrong failure mode."""
    for name, entry in _PAIRING_CANDIDATES.items():
        assert entry["why"].strip(), name
        assert not (entry["counts"] and entry["weighs"]), name


# ---------------------------------------------------------------------------
# The geometric route, and D-219's cross-check made sourced
# ---------------------------------------------------------------------------


def test_klis_route_reproduces_its_own_printed_masses() -> None:
    """The transcription check: volume x 1.11 x 0.34 must return Klis's own printed pg."""
    for name, cell in _KLIS_CELLS.items():
        derived = _volume_fl_to_dry_pg(cell["volume_fl"])
        assert derived == pytest.approx(cell["printed_pg"], rel=0.01), (
            f"{name}: the stated route gives {derived:.2f} pg against a printed "
            f"{cell['printed_pg']} pg -- one of the two is transcribed wrong"
        )


def test_the_sourced_dry_fraction_is_what_d219s_own_upper_edge_implied() -> None:
    """D-219's "~30 %" was silently doing duty as a range, and 0.34 is its top.

    The repair is a *sourcing* one and this guard says so: the printed edges imply
    0.270 and 0.342, so adopting Klis's sourced 0.34 does not contradict D-219 — it names the
    value the upper edge was already computed at and removes the unsourced spread below it.
    """
    implied = tuple(
        pg * 1e-12 / (fl * 1e-12 * _KLIS_DENSITY_G_PER_ML)
        for pg, fl in zip(_D219_PRINTED_CROSSCHECK_PG, _D219_ASSUMED_VOLUME_FL, strict=True)
    )
    assert implied[0] == pytest.approx(0.270, abs=0.005)
    assert implied[1] == pytest.approx(0.342, abs=0.005)
    assert implied[0] < _KLIS_DRY_FRACTION <= implied[1] + 1e-9, (
        "the sourced dry fraction must sit at the TOP of the spread D-219's printed band "
        "implies -- if it falls outside, the cross-check is not the same calculation"
    )


def test_the_sourced_crosscheck_narrows_but_does_not_move_the_settled_band() -> None:
    """37.7-56.6 pg on the sourced fraction, and the settlement is untouched.

    The settled 28-50 pg band comes from the ELEMENTAL route across
    ``biomass_N_fraction``'s uncertainty, not from this cross-check. A cross-check that moved
    the settlement would be a fit, so the guard pins the overlap rather than a replacement.
    """
    sourced = tuple(_volume_fl_to_dry_pg(fl) for fl in _D219_ASSUMED_VOLUME_FL)
    assert sourced[0] == pytest.approx(37.74, abs=0.05)
    assert sourced[1] == pytest.approx(56.61, abs=0.05)
    # narrower on the low side than the printed band, because the unsourced 0.27 edge goes
    assert sourced[0] > _D219_PRINTED_CROSSCHECK_PG[0]
    # and it still contains the settled gram, so nothing about D-219's settlement changes
    lo, hi = _D219_SETTLED_BAND_PG
    assert sourced[0] < _D219_SETTLED_PG < sourced[1]
    assert sourced[0] < hi and sourced[1] > lo, "the two routes must still overlap"


# ---------------------------------------------------------------------------
# THE HEADLINE: the direction every figure in reach points
# ---------------------------------------------------------------------------


def test_every_figure_in_reach_lands_below_branch_one() -> None:
    """Not one per-cell mass available anywhere reaches branch 1's 70.9-91.9 pg.

    Eight figures of three different classes — a geometric estimate, a rule of thumb, the
    engine's settled gram, the elemental estimate at both its old and its re-priced edges, and
    the sourced cross-check — and every one of them is below the demand. **This is a
    direction, not a settlement**: none is a pairing and none is an ale cell in wort, so the
    uniformity says the literature has nothing as heavy as branch 1 needs, not that branch 1
    is false.
    """
    demand_lo, _demand_hi = _D230_BRANCH_ONE_PG_PER_CELL
    figures = {
        "Klis haploid": _KLIS_CELLS["haploid"]["printed_pg"],
        "Klis diploid": _KLIS_CELLS["diploid"]["printed_pg"],
        "Physical Biology of the Cell rule of thumb": _RULE_OF_THUMB_PG,
        "D-219 settled gram": _D219_SETTLED_PG,
        "D-219 settled band, high edge": _D219_SETTLED_BAND_PG[1],
        "D-270 §7 elemental re-pricing, low": _D270_ELEMENTAL_REPRICED_PG[0],
        "D-270 §7 elemental re-pricing, high": _D270_ELEMENTAL_REPRICED_PG[1],
        "D-219 cross-check on the sourced fraction, high": _volume_fl_to_dry_pg(
            _D219_ASSUMED_VOLUME_FL[1]
        ),
    }
    above = {name: pg for name, pg in figures.items() if pg >= demand_lo}
    assert above == {}, (
        f"a figure now reaches branch 1's demand: {above}. D-271's headline is that none does; "
        "if one now can, the ordering has changed and owes a record"
    )


def test_d270s_repricing_is_the_engine_side_estimate_not_branch_ones_demand() -> None:
    """The conflation guard: 47.71-62.12 pg is the estimate, 70.9-91.9 pg is the demand.

    D-271's own first draft read D-270 §7 as having *re-priced the demand*, which would have
    made the hunt look nearly settled. What D-270 moved was the engine-side elemental estimate,
    narrowing the GAP to 1.14-1.93x while leaving D-230's demand exactly where it was.
    """
    demand_lo, demand_hi = _D230_BRANCH_ONE_PG_PER_CELL
    estimate_lo, estimate_hi = _D270_ELEMENTAL_REPRICED_PG
    assert estimate_hi < demand_lo, "the re-priced estimate must still sit below the demand"
    assert demand_hi / estimate_lo == pytest.approx(1.93, abs=0.01)
    assert demand_lo / estimate_hi == pytest.approx(1.14, abs=0.01)


# ---------------------------------------------------------------------------
# Why the size route cannot rescue it
# ---------------------------------------------------------------------------


def test_the_size_route_is_frame_broken_on_the_same_strains() -> None:
    """Two sourced size series disagree by >4x in VOLUME on the same laboratory strains.

    Klis's microscopy volumes make a haploid 44 fL (4.4 um equivalent sphere); Okada's
    flow-cytometry *apparent* diameter makes the same ploidy 7.3 um (204 fL). Combining them —
    reading a mass off one and a size off the other — is the frame error this guard exists to
    forbid.
    """
    for ploidy, klis_key, okada_key in (
        ("haploid", "haploid", "haploid BY4741L"),
        ("diploid", "diploid", "diploid BY4743L"),
    ):
        klis_fl = _KLIS_CELLS[klis_key]["volume_fl"]
        okada_fl = _sphere_volume_fl(_OKADA_APPARENT_DIAMETER_UM[okada_key])
        assert okada_fl / klis_fl > 4.0, (
            f"{ploidy}: the two size frames now agree to within 4x "
            f"({okada_fl / klis_fl:.2f}x) -- the size route may be re-openable, which owes a "
            "record rather than a silent re-use"
        )


def test_branch_ones_demanded_cell_is_chosen_by_which_size_source_you_open() -> None:
    """The demand reads as an ordinary haploid in one frame and 2-3x a diploid in the other.

    That is the operational form of the frame break: the size route does not merely fail to
    settle branch 1, it answers *both ways*, so no verdict may be taken from it.
    """
    demanded_fl = tuple(_dry_pg_to_volume_fl(pg) for pg in _D230_BRANCH_ONE_PG_PER_CELL)
    assert demanded_fl[0] == pytest.approx(187.9, abs=0.5)
    assert demanded_fl[1] == pytest.approx(243.5, abs=0.5)

    # In Klis's frame: 2.3-2.9x a diploid, and above D-219's own assumed wine/ale cell.
    klis_diploid_fl = _KLIS_CELLS["diploid"]["volume_fl"]
    assert demanded_fl[0] / klis_diploid_fl == pytest.approx(2.26, abs=0.02)
    assert demanded_fl[1] / klis_diploid_fl == pytest.approx(2.93, abs=0.02)
    assert demanded_fl[0] > _D219_ASSUMED_VOLUME_FL[1]

    # In Okada's frame: the SAME demand is about one haploid.
    okada_haploid_fl = _sphere_volume_fl(_OKADA_APPARENT_DIAMETER_UM["haploid BY4741L"])
    assert 0.9 < demanded_fl[0] / okada_haploid_fl < 1.3
    assert 0.9 < demanded_fl[1] / okada_haploid_fl < 1.3


# ---------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------


def test_none_of_this_moves_the_engines_gram() -> None:
    """4e-11 g is D-219's DEFINITION of the unit Coleman's fit is counted in.

    A literature pg/cell scores branch 1; it cannot move the gram, and the tripwire says so by
    pinning the conversion these figures are all read against.
    """
    gram_per_cell = cells_per_ml_to_pitch_gpl(1.0) / 1e3
    assert gram_per_cell == pytest.approx(4.0e-11, rel=1e-12), (
        "the engine's gram moved. D-271's figures are literature estimates scored AGAINST it "
        "and none of them licenses moving it -- see D-219"
    )


def test_the_literature_figures_are_scoring_targets_not_parameters() -> None:
    """Prime directive 2 governs numbers the MODEL reads, and ``src/`` reads none of these.

    Same standing as ``TYRELL_CELL_COUNT`` and ``CHEMISTRY_OF_BEER_GROWTH_FOLD``: transcribed
    beside the guards that score them. If a Process is ever built to one, it moves to YAML then.
    """
    from fermentation.parameters import default_data_dir

    haystack = "\n".join(
        path.read_text(encoding="utf-8") for path in sorted(default_data_dir().glob("*.yaml"))
    )
    # Deliberately NOT a search for the bare literals: 0.34 and 1.11 are ordinary numbers and
    # the suite already ships an unrelated `value: 0.34`, so a literal grep fires on
    # coincidence and forbids nothing. A provenance string cannot collide by accident.
    for citation in ("Klis", "Okada", "Physical Biology of the Cell"):
        assert citation not in haystack, (
            f"a parameter now cites {citation}: one of D-271's literature figures has become a "
            "number the MODEL reads. That is a real change of standing and owes a record -- "
            "these are scoring targets, and none of them is a pairing"
        )
