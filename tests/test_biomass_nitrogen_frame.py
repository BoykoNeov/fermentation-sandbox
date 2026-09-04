"""What ``biomass_N_fraction`` actually is in each medium, and what it is not (decision D-270).

D-267 §6 recorded a mismatch: *"the engine ships ``biomass_N_fraction`` = 0.114 … 1.6-1.8x the
nitrogen those two wine-yeast statements imply"*. The comparison was made against the YAML
literal, and **wine never runs on the YAML literal** — :func:`_apply_nitrogen_dependent_yield`
overrides it at the compile boundary from Coleman's ``Y_X/N`` regression, so wine's value is a
function of the must's declared assimilable nitrogen. Across the suite's own musts it spans
**0.0362 to 0.1068**, a factor 2.95, and it *straddles* the sourced range rather than sitting
above it.

That is the frame this module pins. The quantity is not a measured composition in wine at all:
it is defined as ``1/Y_X/N``, a yield residual, and
:func:`test_the_compiled_fraction_hands_colemans_own_gram_back_exactly` proves it by feeding the
compiled value through ``convert.py``'s independent per-cell check, which then returns **exactly**
the 4e-11 g the yield was built on. The check is informative only when fed a nitrogen fraction
Coleman had no hand in, which is why ``convert.py`` uses the static elemental 0.114 — correct as
it stands, and not this module's to change.

What survives as a real mismatch is **beer's**, because beer is gated off the override and keeps
the static value. Beer's ensemble draws over ``[0.08, 0.14]`` while every sourced estimate lies in
``[0.0640, 0.0833]``, so 94 % of the band is above all of them. **Moving it is refused** — the
growth law is the identity ``dX = YAN/f_N``, so the sourced range raises beer's biomass ceiling
1.37-1.78x and takes the growth extent from 5.378x to 7.36-9.58x against a counted 2.918-3.483x
and a printed 4-5x. Refused, flagged, not built.

The sourcing **strengthens** D-230's arithmetic refusal of the nitrogen-partition candidate
rather than weakening it, and that is asserted here so a reader cannot take it the other way.

Receipts: ``M:\\claud_projects\\temp\\ferment\\d270-nitrogen-frame\\``.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from fermentation.parameters import Parameter, default_data_dir, load_parameters
from fermentation.scenario import (
    Scenario,
    TemperaturePoint,
    amino_acid_dose_nitrogen_mgl,
    compile_scenario,
)
from fermentation.units import cells_per_ml_to_pitch_gpl
from tests.test_precursor_fates import _D267_PROTEIN_FRACTION_ANCHORS

#: The two conversion conventions the sourced protein statements carry with them. Neither is a
#: choice made here: the *Concise Encyclopedia* states its own figure is "based on N x 6.25", and
#: the *Understanding Wine Chemistry* note states protein is one-sixth nitrogen. A crude-protein
#: factor asserts TOTAL cell nitrogen, not protein nitrogen, which is what makes these statements
#: comparable to ``biomass_N_fraction`` at all.
_CRUDE_PROTEIN_FACTOR = 6.25
_UWC_NITROGEN_SHARE_OF_PROTEIN = 1.0 / 6.0

#: Total cell nitrogen [g N / g dry weight] the sourced wine-yeast protein statements imply.
#: DERIVED from D-267 §2's anchors through each statement's own convention — never a literal, so
#: a re-reading of either source moves this and every guard below with it.
_SOURCED_N_FRACTION = {
    "Concise Encyclopedia, 40 % crude protein": (
        _D267_PROTEIN_FRACTION_ANCHORS["Concise Encyclopedia, wine yeast, low edge"]
        / _CRUDE_PROTEIN_FACTOR
    ),
    "Concise Encyclopedia, 45 % crude protein": (
        _D267_PROTEIN_FRACTION_ANCHORS["Concise Encyclopedia, wine yeast, high edge"]
        / _CRUDE_PROTEIN_FACTOR
    ),
    "Understanding Wine Chemistry, 50 % protein at 1/6 N": (
        _D267_PROTEIN_FRACTION_ANCHORS["Understanding Wine Chemistry 2nd ed, chapter note"]
        * _UWC_NITROGEN_SHARE_OF_PROTEIN
    ),
}
_SOURCED_LO = min(_SOURCED_N_FRACTION.values())
_SOURCED_HI = max(_SOURCED_N_FRACTION.values())

#: Every declared assimilable nitrogen [mg N/L] the suite's wine scenarios evaluate the yield fit
#: at, low to high. The bare values; the three dosed arms below carry the D-244 migration sum.
_WINE_DECLARED_YAN = (50.0, 80.0, 100.0, 150.0, 250.0, 300.0, 330.0, 350.0, 400.0, 500.0)
_WINE_DOSED_ARMS = ((250.0, 2.0), (250.0, 4.0), (80.0, 2.0))

#: Coleman, Fish & Block 2007 Fig. 4, and the gram his counts were converted at (D-219).
_COLEMAN_A0, _COLEMAN_A1 = 3.50, -3.61e-3
_COLEMAN_GRAM_PER_CELL = 4.0e-11
_COLEMAN_REFERENCE_YAN = 330.0

#: Tyrell 2013's counted crop, as D-230 priced it: the cell nitrogen that crop would demand at
#: the engine's own gram, and the per-cell dry mass its branch-1 reading implies instead.
_D230_DEMANDED_N_FRACTION = (0.202, 0.262)
_D230_BRANCH_ONE_PG_PER_CELL = (70.9, 91.9)
_BEER_SHIPPED_EXTENT_FOLD = 5.378
_TYRELL_COUNTED_FOLD = (2.918, 3.483)
_CHEMISTRY_OF_BEER_PRINTED_FOLD = (4.0, 5.0)


def _compiled_f_n(medium: str, initial: dict[str, float]) -> Parameter:
    """The ``biomass_N_fraction`` a COMPILED scenario actually carries — never the YAML literal."""
    scenario = Scenario(
        name=f"d270-{medium}",
        medium=medium,
        initial=initial,
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0 if medium == "wine" else 15.0)
        ],
        duration_days=10.0,
    )
    return compile_scenario(scenario, strict=True).parameters["biomass_N_fraction"]


def _wine_f_n(yan_mgl: float, amino_acids_gpl: float | None = None) -> float:
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if amino_acids_gpl is not None:
        initial["amino_acids_gpl"] = amino_acids_gpl
        # D-244 migration form: a dosed fixture declares the sum, and the sum is what the fit
        # is evaluated at. This is the arm that reaches the hold.
        initial["yan_mgl"] += amino_acid_dose_nitrogen_mgl(initial)
    return _compiled_f_n("wine", initial).value


def _beer_f_n() -> float:
    return _compiled_f_n(
        "beer",
        {
            "glucose_gpl": 15.0,
            "maltose_gpl": 60.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.3984,
        },
    ).value


def _elemental_pg_per_cell(f_n: float, yan_mgl: float = _COLEMAN_REFERENCE_YAN) -> float:
    """``convert.py``'s independent route: Coleman's cells per gram N, priced at a composition."""
    y_xn = math.exp(_COLEMAN_A0 + _COLEMAN_A1 * yan_mgl)
    cells_per_g_n = y_xn / _COLEMAN_GRAM_PER_CELL
    return (1.0 / f_n) / cells_per_g_n * 1e12


def test_the_sourced_nitrogen_range_is_derived_from_d267s_anchors_and_not_transcribed():
    """The derivation itself, so a re-read of either source cannot leave this module stale.

    Each anchor goes through ITS OWN convention. Mixing them — applying ``N x 6.25`` to the
    *Understanding Wine Chemistry* figure, say — would silently invent a fourth number, and the
    two conventions are not interchangeable: 6.25 and 6.0 differ by 4 %.
    """
    assert _SOURCED_N_FRACTION["Concise Encyclopedia, 40 % crude protein"] == pytest.approx(0.0640)
    assert _SOURCED_N_FRACTION["Concise Encyclopedia, 45 % crude protein"] == pytest.approx(0.0720)
    assert _SOURCED_N_FRACTION[
        "Understanding Wine Chemistry, 50 % protein at 1/6 N"
    ] == pytest.approx(0.08333, abs=1e-5)
    assert pytest.approx((0.0640, 0.08333), abs=1e-5) == (_SOURCED_LO, _SOURCED_HI)

    # The fourth anchor D-267 §2 lists is deliberately NOT in the range: van Gulik & Heijnen's
    # 0.42 is a protein MASS from a chemostat balance, with no stated nitrogen convention, so
    # converting it would mean choosing one on the source's behalf.
    assert "van Gulik & Heijnen 1995 Table I (after Verduyn)" in _D267_PROTEIN_FRACTION_ANCHORS
    assert not any("van Gulik" in key for key in _SOURCED_N_FRACTION)


def test_wines_nitrogen_fraction_is_a_function_of_the_must_and_not_a_composition():
    """The frame correction: wine has no single ``biomass_N_fraction`` to compare to a source.

    ``_apply_nitrogen_dependent_yield`` sets it to ``1/Y_X/N(YAN)``, and ``ln(Y_X/N)`` falls
    linearly in YAN, so the fraction RISES with the must's nitrogen — monotonically, by a factor
    2.95 across the musts this suite actually runs, until the D-244 hold pins it at the fit's top.
    A cell composition cannot depend on the medium's nitrogen by a factor of three; a yield
    residual can, and must.
    """
    values = [_wine_f_n(yan) for yan in _WINE_DECLARED_YAN]

    assert all(later >= earlier for earlier, later in zip(values, values[1:], strict=False)), (
        f"f_N must be monotone non-decreasing in declared YAN; got {values}"
    )
    assert values[0] == pytest.approx(0.0362, abs=5e-4), "50 mg N/L (Varela's low arm)"
    assert max(values) == pytest.approx(0.1068, abs=5e-4), "held at the fit's 350 mg N/L top"
    assert max(values) / min(values) == pytest.approx(2.95, abs=0.05), (
        "the span across the suite's own musts is the whole point: this is not a constant"
    )
    # The hold: everything above the fit's domain shares one value, and that is D-244's
    # epistemic rule, not saturation of the yield.
    held = [_wine_f_n(yan) for yan in (350.0, 400.0, 500.0)]
    assert held[0] == pytest.approx(held[1]) == pytest.approx(held[2])


def test_the_compiled_fraction_hands_colemans_own_gram_back_exactly():
    """Why the compiled value cannot corroborate anything: the check closes on itself.

    ``convert.py`` prices Coleman's cells-per-gram-nitrogen with an elemental composition *he had
    no hand in* and gets 34.9 pg against his assumed 40 — a 13 % agreement between independent
    routes. Feed it the COMPILED fraction instead and the agreement becomes exact to machine
    precision, because ``f_N = 1/Y_X/N`` makes the two sides the same statement:

        (1 / f_N) / (Y_X/N / 4e-11) = Y_X/N * 4e-11 / Y_X/N = 4e-11

    That is the proof that wine's ``biomass_N_fraction`` carries no compositional information,
    and it is also why ``convert.py``'s use of the static 0.114 is right and must stay.
    """
    compiled = _wine_f_n(_COLEMAN_REFERENCE_YAN)
    assert _elemental_pg_per_cell(compiled) == pytest.approx(
        _COLEMAN_GRAM_PER_CELL * 1e12, rel=1e-12
    ), "the compiled fraction is Coleman's own inverse; this identity is exact, not approximate"

    static = load_parameters(default_data_dir() / "beer_generic.yaml")["biomass_N_fraction"].value
    assert _elemental_pg_per_cell(static) == pytest.approx(34.9, abs=0.1), (
        "the independent route, which is only independent because 0.114 is not Coleman's"
    )


def test_the_suites_musts_straddle_the_sourced_range_rather_than_sitting_above_it():
    """D-267 §6's "1.6-1.8x above" is a statement about a literal wine does not read.

    Thirteen wine arms: **four below** the sourced range, **one inside**, **eight above** — and
    the eight include every dosed fixture, which reaches the hold by declaring the D-244 sum.
    The low arms miss by more than the high ones: 0.0362 against a 0.0640 floor is 1.77x low,
    while the held 0.1068 against a 0.0833 ceiling is 1.28x high.
    """
    values = [_wine_f_n(yan) for yan in _WINE_DECLARED_YAN]
    values += [_wine_f_n(yan, dose) for yan, dose in _WINE_DOSED_ARMS]

    below = [v for v in values if v < _SOURCED_LO]
    inside = [v for v in values if _SOURCED_LO <= v <= _SOURCED_HI]
    above = [v for v in values if v > _SOURCED_HI]
    assert (len(below), len(inside), len(above)) == (4, 1, 8), (
        f"census over {len(values)} wine arms: {len(below)} below / {len(inside)} inside / "
        f"{len(above)} above the sourced {_SOURCED_LO:.4f}-{_SOURCED_HI:.4f}"
    )
    assert _SOURCED_LO / min(below) == pytest.approx(1.77, abs=0.03)
    assert max(above) / _SOURCED_HI == pytest.approx(1.28, abs=0.03)


def test_beers_band_lies_almost_entirely_above_every_sourced_estimate():
    """The mismatch that survives the frame correction, stated as a share of the band.

    Beer is gated OFF the override — ``_apply_nitrogen_dependent_yield`` returns unchanged when
    the regression coefficients are absent — so beer runs the static elemental value and its
    ensemble draws over the static band. The overlap with the sourced range is the sliver
    ``[0.0800, 0.0833]``: **5.6 % of the band's width**. Not "never draws a supported value", but
    close enough that the nominal itself is above every sourced estimate.
    """
    spec = load_parameters(default_data_dir() / "beer_generic.yaml")["biomass_N_fraction"]
    assert _beer_f_n() == pytest.approx(spec.value), "beer keeps the static value; no override"

    low, high = spec.uncertainty.low, spec.uncertainty.high
    overlap = max(0.0, min(high, _SOURCED_HI) - max(low, _SOURCED_LO))
    assert overlap / (high - low) == pytest.approx(0.056, abs=0.005), (
        f"band [{low}, {high}] overlaps the sourced [{_SOURCED_LO:.4f}, {_SOURCED_HI:.4f}] over "
        f"{overlap:.4f}, i.e. {overlap / (high - low):.1%} of its width"
    )
    assert spec.value > _SOURCED_HI
    assert spec.value / _SOURCED_HI == pytest.approx(1.368, abs=0.005)


def test_the_sourced_composition_is_refused_for_beer_by_beers_own_growth_extent():
    """Why the mismatch is flagged and not repaired: the repair is priced and it is unaffordable.

    Beer's growth is nitrogen-limited, so the ceiling is the identity ``dX = YAN/f_N`` — lowering
    the nitrogen fraction builds proportionally MORE biomass out of the same wort nitrogen. The
    extent is already over both of its targets; the sourced range takes it further out of both.

    This is the same arithmetic that forbids harmonising beer onto wine's regression, which at
    200 mg N/L gives 0.062 and would take the overshoot past 9x. Nothing here is a new
    measurement; it is the identity applied to a sourced number for the first time.
    """
    static = load_parameters(default_data_dir() / "beer_generic.yaml")["biomass_N_fraction"].value
    folds = {
        name: _BEER_SHIPPED_EXTENT_FOLD * static / f for name, f in _SOURCED_N_FRACTION.items()
    }

    assert min(folds.values()) == pytest.approx(7.36, abs=0.05)
    assert max(folds.values()) == pytest.approx(9.58, abs=0.05)
    assert min(folds.values()) > max(_CHEMISTRY_OF_BEER_PRINTED_FOLD), (
        "even the mildest sourced composition puts the extent above the PRINTED 4-5x fold, "
        "which is the target the model was already closest to (D-258)"
    )
    assert min(folds.values()) / max(_TYRELL_COUNTED_FOLD) > 2.0, (
        "and more than 2x Tyrell's counted crop, against 1.54x for the shipped value"
    )


def test_the_sourcing_strengthens_d230s_arithmetic_refusal_it_does_not_weaken_it():
    """A reader could take a LOWER real composition as making the partition candidate easier.

    It makes it harder. D-230 refused "not all assimilated N reaches suspended biomass" because
    reproducing Tyrell's counted crop demands cells at 20-26 % nitrogen, outside the shipped
    band's 14 % top by 1.44-1.87x. Measured against the sourced range instead, the same demand is
    outside by **2.42-4.09x**. The refusal's ground moves further from admissible, not closer.
    """
    spec = load_parameters(default_data_dir() / "beer_generic.yaml")["biomass_N_fraction"]
    lo, hi = _D230_DEMANDED_N_FRACTION

    against_band = (lo / spec.uncertainty.high, hi / spec.uncertainty.high)
    against_source = (lo / _SOURCED_HI, hi / _SOURCED_LO)
    assert against_band == pytest.approx((1.44, 1.87), abs=0.01)
    assert against_source == pytest.approx((2.42, 4.09), abs=0.01)
    assert against_source[0] > against_band[0] and against_source[1] > against_band[1], (
        "sourcing the composition must move the demanded fraction FURTHER outside, or D-230's "
        "refusal has been weakened rather than strengthened and the record says the opposite"
    )


def test_repricing_the_independent_per_cell_check_moves_it_toward_the_counted_branch():
    """Reported, flagged, and NOT acted on: D-219's gram is a unit definition, not an estimate.

    ``convert.py``'s corroboration of the 4e-11 g gram runs Coleman's cells-per-gram-nitrogen
    through an elemental composition. Swap Roels' generic 0.114 for the sourced wine-yeast range
    and the same route gives **47.7-62.1 pg/cell** instead of 34.9 — the low edge still inside
    D-219's 28-50 pg band, the high edge outside it.

    The direction matters because D-230's branch 1 reads Tyrell's counts as 70.9-91.9 pg/cell and
    D-232 could only widen that ambiguity. The sourced composition narrows the gap from
    2.03-2.63x to 1.14-1.93x. That is evidence on a three-way residue, not a settlement: the
    gram is the unit wine's fitted parameters live in, and re-opening it is D-219's to re-open.
    """
    repriced = sorted(_elemental_pg_per_cell(f) for f in _SOURCED_N_FRACTION.values())
    assert repriced[0] == pytest.approx(47.71, abs=0.05)
    assert repriced[-1] == pytest.approx(62.12, abs=0.05)

    branch_lo, branch_hi = _D230_BRANCH_ONE_PG_PER_CELL
    static = load_parameters(default_data_dir() / "beer_generic.yaml")["biomass_N_fraction"].value
    shipped_pg = _elemental_pg_per_cell(static)
    assert (branch_lo / shipped_pg, branch_hi / shipped_pg) == pytest.approx((2.03, 2.64), abs=0.01)
    assert (branch_lo / repriced[-1], branch_hi / repriced[0]) == pytest.approx(
        (1.14, 1.93), abs=0.01
    )
    assert branch_lo / repriced[-1] > 1.0, (
        "the count-derived branch must still sit ABOVE the elemental route, or the sourcing has "
        "closed a residue D-232 left three-way and that needs its own record"
    )


def test_the_unsourced_seven_to_twelve_percent_range_no_longer_stands_unattributed():
    """The prose repair this record ships, pinned so it cannot silently revert.

    ``convert.py``'s wet-vs-dry argument leaned on "a real 7-12 %" — an unsourced range doing
    sanity-check duty in three files. It is replaced by the sourced 6.4-8.3 %, which points the
    same way and points further: 33 % nitrogen on a wet reading is more absurd against a lower
    ceiling, not less.
    """
    text = (
        Path(__file__).parent.parent / "src" / "fermentation" / "units" / "convert.py"
    ).read_text(encoding="utf-8")
    assert "7-12" not in text, "the unsourced range is back in convert.py"
    assert "6.4-8.3" in text, "convert.py must cite the sourced range instead"


def test_the_engines_gram_is_untouched_by_this_record():
    """A tripwire, because the record's most quotable number is one it deliberately did not use.

    47.7-62.1 pg is not a proposal. If ``cells_per_ml_to_pitch_gpl`` ever moves off Coleman's
    4e-11 g, that is D-219's settlement being re-opened and it owes its own record — this guard
    fires first.
    """
    assert cells_per_ml_to_pitch_gpl(1e6) == pytest.approx(1e6 * 1e3 * _COLEMAN_GRAM_PER_CELL)
