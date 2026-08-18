"""Beer's yeast-produced organic acids — the producer, and the pH it predicts (decision D-180).

D-179 gave beer charge-active acid slots and left them **inert**: dosed at pitch, touched by
nothing, so a modelled beer's pH sat where ``initial_ph`` put it for the whole run. D-180 gives
four of them a producer, which forces the seeds to become **wort** levels and turns beer's
finished pH from an input into a prediction.

The source is unusually good for this engine. Tyrell et al. (2013) fermented ONE 12 °P wort
with FOUR strains and published, for the same ferments, every organic acid's day-0-to-day-7
course (Figs 6-14) *and* the pH and extract curves (Fig 4). So the yields are a measured
difference on one wort, and the pH is a FREE prediction — nothing in the parameter file is
fitted to it. The acceptance section below is therefore a real external comparison rather than
a round-trip, with one exception which says so in its own name.

**D-181 adds the sink half** (section 5): the three wort acids Tyrell measure FALLING, which
beer's model previously could not lose. That is why this file's headline numbers moved DOWN.
**D-182 then adds dissolved CO₂** to the charge balance — the second of the two terms D-180
named as omitted, pulling the opposite way — which is why they moved back UP.

**THE PH COMPARISON RUNS IN THE DEGASSED FRAME SINCE D-208, AND IT HALVES.** Tyrell's pH is a
*decarbonated* reading (their cited MEBAK II 2.14 is "pH (EBC)"; Analytica-EBC 9.35's scope line
is "the determination of pH at 20 °C of DECARBONATED beer"), while ``ph_of_state`` reports the pH
inside the vessel — worth 0.29 pH at day 7 and 0.50 at day 1. Against a measured drop of **0.81**
pH — the mean of the extreme strains, which is what ``measured_drop`` below computes; D-180's
prose quotes the four-strain mean 0.8125 and the two must not be mixed — the model gives
**43.2-62.9 %** at nominal across the sampled ``pKa_peptide_buffer`` band and **8.3-82.7 % over
the joint band** of all TEN drawn quantities, so **nothing in the band reaches the measurement**.
The numbers this file quoted for four beats — **77.8-97.3 %** at nominal, **64.0-109.7 %** joint,
with a corner reaching — are the IN-VESSEL quantity; they are still pinned, as a property of the
model rather than as an agreement with anything published.

The history of those in-vessel numbers is most of the story of this axis: 63-92 % and 41-105 % at
D-180 (with a corner reaching), 42.7-62.2 % and 7.6-82.2 % at D-181 (with nothing reaching),
77.6-97.0 % and 63.8-109.4 % at D-182 — because D-181 removed an error that was propping the
agreement up and D-182 supplied the term that was genuinely missing — and then **D-183's +0.2 pp,
which was pre-registered before a line of ``src/`` changed**. That last move is deliberately too
small to be a result: acetic's endpoint went from 116.06 to 117.75 mg/L at 0.1233 pp of agreement
per ppm, so a rate-law change that reshapes the whole acetic curve is **headline-neutral by
construction** and must not be credited here.

**D-208 is the rest of that story and it is not another move along the same line**: the D-182 rise
was real chemistry scored in the wrong frame, and re-scoring it returns the degassed comparison to
within ~0.7 pp of the pre-carbonic D-181 numbers — close, but not identical, because D-183's
acetic rate law and the seeds moved in between. Nothing was un-built to get there: the term stays
in ``ph_of_state``, where every pH-reading Process needs it.

**D-183 also removes one of the two shape failures from the "unmodelled" list and hardens the
other.** Acetic's producer moved from the sugar flux to **growth** (``AceticAcidOverflow``),
because mapping Tyrell's Fig 13 onto their Fig 4 puts 86 % of acetic's rise inside the first
15 % of the fermentative flux. That halves the shape error against their measured days 1-7
(RMSE 61.6 → 32.5 ppm) and fixes *when* the acid appears — but the modelled curve is still
**MONOTONE**. The mid-ferment spike is **still not modelled**, and the excretion/re-assimilation
pair D-180 §9 proposed for it was built as a probe and **REFUSED on measurement**: the same
figures falsify its re-assimilation half too (half the measured fall happens at zero
fermentative flux), the decline cannot discriminate first-order from constant-rate at the ±3 ppm
read tolerance, no floor is identifiable, and every law the data admits makes the endpoint a
function of the solver horizon. Lactic's late rise is untouched.

**This still is not validation, and since D-208 it is a named failure rather than a shortfall.**
In the measured frame the nominal falls short by 37-57 %, the day-7 pH is above the four-strain
envelope for every retained-CO₂ fraction in [0, 1], and the acceptance claim is carried by a
``strict=True`` xfail (``test_the_model_reaches_tyrells_measured_beer_ph``) so it cannot be closed
silently — the D-188 idiom. The acetic spike and lactic's late rise remain unmodelled and are
charge-balance non-terms either way. No test here is named or phrased as validating the produced
acids alone.
"""

import re
import shutil

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.core.acidbase import charge_balance_is_populated
from fermentation.core.chemistry import M_NITROGEN, carbon_mass_fraction, sugar_species
from fermentation.core.kinetics import (
    ACETIC_SLOT,
    ORGANIC_ACID_SPECS,
    WORT_ACID_SINKS,
    AceticAcidOverflow,
    OrganicAcidExcretion,
    WortAcidRemoval,
)
from fermentation.core.kinetics.carbon_routing import fermentative_uptake_rates
from fermentation.core.kinetics.uptake import SugarUptakeToEthanolCO2
from fermentation.core.media import beer_schema, get_medium, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.parameters import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import cells_per_ml_to_pitch_gpl

# --------------------------------------------------------------------------------------
# Tyrell et al. (2013), BrewingScience 66:75-84 — the matched dataset. FIGURE READS; the
# rendered crops and the read tolerances live with the D-180 working files. These are the
# numbers the acceptance section compares against and they are written here ONCE.
# --------------------------------------------------------------------------------------

#: EBC-tube trial, day 7, four strains: the measured range each acid finished in [mg/L].
TYRELL_BEER_PPM = {
    "acetic": (105.0, 126.0),
    "lactic": (87.0, 142.0),
    "succinic": (32.0, 76.0),
    "malic": (81.0, 116.0),
}
#: Acetic's growth-linked yield parameter (D-183). Named once, here, because three tests and
#: the isolability check all have to mean the same name.
ACETIC_YIELD_PARAM = "Y_acetic_biomass_beer"

#: Fig. 13, the acetic course of those same ferments, as the **four-strain mean per day** —
#: the *peak of the mean* (170.25 ppm at day 2), NOT the mean of the four strains' own peaks
#: (182.75 ppm), which is a different quantity and must not be substituted. Both are
#: AUTHOR-CONSTRUCTED means over FIGURE READS at ±3 ppm. Day 0 and day 7 are cross-checked
#: against D-180's independently recorded seed and per-strain deltas and agree exactly.
TYRELL_ACETIC_MEAN_PPM = {
    0: 59.00,
    1: 145.00,
    2: 170.25,
    3: 161.75,
    4: 153.25,
    5: 141.50,
    6: 130.00,
    7: 117.75,
}

#: Fig. 4, same ferments: wort pH and the four strains' day-7 beer pH.
TYRELL_WORT_PH = 5.65
TYRELL_BEER_PH = (4.78, 4.90)

#: The temperature every Tyrell comparison in this file is scored at — and its provenance,
#: written here ONCE because it was an inference for nine records (D-217).
#:
#: Tyrell ran **two** trials. §2.4.1's laboratory flasks are the ones with a printed
#: fermentation temperature — *"Fermentation temperature was isothermal at 20 °C"* — and they
#: are NOT the source of Fig. 4. Fig. 4 is §2.4.2's 2 L EBC tubes, the trial whose *"yeast
#: concentration, pH and extract development was checked daily"*, and the only temperature the
#: paper prints for it is the FILL: *"cooled down to 15 °C"*. D-178 transcribed that correctly
#: as "pitched at 15 °C"; from D-207 on, the archive said "Tyrell ferments at 15 °C".
#:
#: The source does settle the direction, in §3.2's comparison of the two scales: *"higher oxygen
#: intake, **higher fermentation temperature** and lower hydrostatic pressure of **lab scale**
#: fermentation"*. The lab scale is the 20 °C one, so the tube trial ran BELOW 20 °C. 15 °C is
#: the only value the paper offers for it, and it is now sourced in direction rather than
#: assumed — but the exact value is still not printed, which is why the two tests in section 11
#: pin what depends on it.
TYRELL_TRIAL_CELSIUS = 15.0

#: Fig. 4's pH panel as a COURSE — the four-strain envelope ``(low, high)`` per day, in the
#: paper's own daily sampling ("yeast concentration, pH and extract development was checked
#: daily", §2.4.2). D-207 transcribed it; D-180 had read this same panel for its two endpoints
#: only, which is why the pH agreement was an endpoint number for four beats.
#:
#: **An envelope, NOT a mean.** The four strains are resolved individually and no mean is
#: constructed, because a merged marker blob's centroid is not either marker's position — the
#: trap D-180's ``TYRELL_ACETIC_MEAN_PPM`` note names for the acetic course ("peak of the mean",
#: never "mean of the peaks"). Day 7's envelope is also the same *shape* of quantity as
#: ``TYRELL_BEER_PH`` above, so the two are directly comparable.
#:
#: **Provenance: a PIXEL read, calibrated and cross-checked, not an eye read.** Least-squares
#: y-calibration through the four printed gridlines (worst residual 0.0028 pH); the pH-5.0 line
#: was withheld from the fit and is predicted at the exact midpoint of its neighbours; and day 0
#: — where all four strains MUST coincide, being one wort — reads as an envelope of span 0.002
#: pH centred on 5.6507 against the 5.65 D-180 recorded independently. The paper's prose is a
#: fourth check: "only strain 7 shows a slightly higher pH-value throughout fermentation", and
#: strain 7 is the top trace here.
#:
#: **READ TOLERANCE 0.024 pH**, and it is the *disagreement between two independent reads of
#: this figure*, not the extraction precision (0.0028). At day 7 this read gives 4.804-4.916
#: where D-180's eye read gave 4.78-4.90. Nothing is re-anchored on that: moving the shipped
#: ``measured_drop`` from 0.81 to 0.792 would shift the headline fraction ~2 %, and that is a
#: decision for a beat that means to take it, not a side effect of transcribing a curve.
#:
#: **Only day 7 is read by an assert, and only inside an ``xfail``** (D-208). D-207 shipped the
#: whole course as data no assert read, on the ground that a pin on the day-1 SHAPE would encode
#: the defect it measured and a correct fix would have to delete it. That still holds for days
#: 1-6. Day 7 is different: it is a LEVEL the model misses in every measurement frame, so an
#: expected-fail assert on it names the gap instead of protecting it, and a fix turns it green.
#: **Day 0's pair is DEGENERATE and is stored ``(low, high)`` like every other day.** The four
#: strains share one wort, so their true span is zero; the extraction returned 5.652 and 5.651,
#: i.e. zero to within the 0.0028 pH extraction precision. It is ordered here rather than left
#: as the raw ``(5.652, 5.651)`` the probe printed, so that ``lo <= x <= hi`` cannot fail on day
#: 0 for a correct value — the span is noise, not a band.
TYRELL_PH_COURSE = {
    0: (5.651, 5.652),
    1: (5.258, 5.377),
    2: (4.871, 5.063),
    3: (4.804, 5.036),
    4: (4.764, 4.956),
    5: (4.790, 4.936),
    6: (4.804, 4.922),
    7: (4.804, 4.916),
}
#: The tolerance above, named so a future assert cannot quietly pick a friendlier one.
TYRELL_PH_READ_TOL = 0.024
#: Fermentable sugar of that wort, derived in beer_acids.yaml's group-2 header from the
#: printed extract curve (apparent → real attenuation → g/L). The whole yield anchor.
TYRELL_SUGAR_GPL = 82.2388545

#: Figs 9 (malic), 10 (lactic) and 14 (succinic) as COURSES — the four-strain mean per day, the
#: same shape of quantity as :data:`TYRELL_ACETIC_MEAN_PPM` above and read at the same ±2 ppm
#: figure tolerance ``beer_acids.yaml`` states for them.
#:
#: **These figures' interiors were never read until D-215.** D-180 mined Figs 6-14 for two points
#: each (day 0, day 7) to build the yields; D-183 transcribed only Fig 13's interior, which is what
#: retired acetic's flux-linked rate law. The other three carried a rate law nothing had ever
#: scored against the days it claims to predict — the same endpoint-scored/trajectory-unscored gap
#: D-207 found in the pH panel of this very figure set.
#:
#: **Anchored at BOTH ends against numbers recorded in earlier beats**, which is what makes the
#: transcription checkable rather than merely careful (the argument D-183 used for Fig 13):
#: lactic day 0 reads 47.75 against a recorded seed of 48 — the 0.25 is strain 7 reading 47 where
#: the other three read 48, i.e. a per-strain read on one wort, NOT a disagreement with the seed;
#: do not "correct" it to 48.0, that would be editing a transcription. Its day-7 span reads 88-143
#: against a recorded 87-142, with strain 15 at +40 and strain 6 at +95 against
#: ``Y_lactic_sugar_beer``'s
#: recorded band edges of +39 and +94; malic reads 78 / 81-116 against a recorded 78 / 81-116,
#: with strain 15 at +3 ("produces essentially none") and strain 31 at 1.57× the mean ("half as
#: much again"); succinic reads 15 / 32-76 against a recorded 15 / 32-76.
#:
#: **Succinic is the CONTROL, not a third defendant** — see
#: ``test_the_three_flux_linked_acid_courses_are_mistimed``. Its measured shape is the one closest
#: to the shipped rate law, so a test that failed on all three alike would not be measuring timing.
TYRELL_ACID_COURSE_PPM: dict[str, dict[int, float]] = {
    "lactic": {0: 47.75, 1: 54.75, 2: 58.25, 3: 76.5, 4: 103.25, 5: 104.0, 6: 113.0, 7: 119.5},
    "malic": {0: 78.0, 1: 75.75, 2: 77.0, 3: 76.75, 4: 98.25, 5: 93.0, 6: 96.75, 7: 102.25},
    "succinic": {0: 15.0, 1: 23.75, 2: 36.25, 3: 48.75, 4: 59.25, 5: 56.25, 6: 61.25, 7: 61.25},
}
#: The ±2 ppm figure read tolerance, named so a later assert cannot quietly pick a friendlier one.
TYRELL_ACID_COURSE_READ_TOL = 2.0

#: Fig. 4's EXTRACT panel as a course — the fraction of the fermentable extract consumed by each
#: day. Transcribed at D-183 §2 (apparent extract → real attenuation → g/L over this file's
#: :data:`TYRELL_SUGAR_GPL` divisor) and used there to falsify acetic's flux-linkage.
#:
#: **It has never been scored against the model's own attenuation, and that is the gap D-215
#: found.** Fig. 4 has three panels and each has now been read by a different beat for a different
#: purpose — D-180 the extract endpoints, D-207 the pH panel, D-211 the cell-count panel — but no
#: test has ever asked whether this engine ferments Tyrell's wort on Tyrell's schedule. It does
#: not: see ``test_the_model_ferments_tyrells_wort_on_tyrells_schedule``.
TYRELL_FLUX_FRACTION: dict[int, float] = {0: 0.0, 1: 0.150, 2: 0.594, 5: 1.003, 7: 0.997}

#: slot -> the wort-seed parameter feeding it, for the three acids D-181 drains. Their bands
#: are drawn, so the joint-band arm varies them; measured, they are worth ~1.5 % of the drop.
WORT_SEED_PARAMS = {
    "pyruvic": "pyruvic_typical_wort",
    "formic": "formic_typical_wort",
    "oxalic": "oxalic_typical_wort",
}

#: Tyrell's wort as this engine states it: the derived fermentable extract split in the
#: customary all-malt proportions, pitched at their 15 °C, anchored at their measured wort pH.
TYRELL_SCENARIO = {
    "glucose_gpl": 0.15 * TYRELL_SUGAR_GPL,
    "maltose_gpl": 0.70 * TYRELL_SUGAR_GPL,
    "maltotriose_gpl": 0.15 * TYRELL_SUGAR_GPL,
    "yan_mgl": 200.0,
    "pitch_gpl": 1.0,
    "initial_ph": TYRELL_WORT_PH,
}


@pytest.fixture(scope="module")
def beer_params():
    data = default_data_dir()
    return load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    )


def _run(initial: dict[str, float], *, days: float = 14.0, celsius: float = 15.0):
    compiled = compile_scenario(
        Scenario(
            name="d180",
            medium="beer",
            initial=initial,
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
            duration_days=days,
        )
    )
    return compiled, compiled.run()


def _beer_state(schema: StateSchema, **over: float | list[float]) -> FloatArray:
    base: dict[str, float | list[float]] = {
        "X": 1.5,
        "S": [10.0, 40.0, 8.0],
        "E": 20.0,
        "N": 0.1,
        "T": 288.15,
        "CO2": 20.0,
    }
    base.update(over)
    return schema.pack(base)


# ======================================================================================
# 1. The shared rate helper — the thing that makes a yield a yield
# ======================================================================================


def test_uptake_rates_helper_is_bitwise_the_uptake_process(beer_params):
    """The extracted helper and the Process that ferments must agree to the last bit.

    ``OrganicAcidExcretion`` books ``Y · Σ r_i`` against the SAME per-slot rates
    ``SugarUptakeToEthanolCO2`` consumes. If the two ever computed the flux differently, ``Y``
    would silently stop being grams-of-acid-per-gram-of-sugar and become grams per gram of
    something-nearly-that — a failure with no symptom, which is exactly the drift D-106 caught
    when two callers "computing the same thing" agreed by luck until one changed.

    Pinned BITWISE rather than to a tolerance, deliberately: the helper was extracted from the
    Process verbatim and in its original operation order, so anything less than bit-equality
    means the arithmetic moved [[feedback-pin-tolerance-vs-solver-tolerance]].
    """
    schema = beer_schema()
    params = beer_params.resolve()
    uptake = SugarUptakeToEthanolCO2()
    s_slice = schema.slice("S")

    # Across the states that actually exercise the branches: full wort, mid-ferment with
    # glucose gone (repression released), a single sugar left, dryness, and no biomass.
    for sugars in ([15.0, 70.0, 15.0], [0.0, 40.0, 12.0], [0.0, 0.0, 5.0], [0.0, 0.0, 0.0]):
        for x in (2.0, 0.0):
            y = _beer_state(schema, S=sugars, X=x)
            rates = fermentative_uptake_rates(y, schema, params)
            d = uptake.derivatives(0.0, y, schema, params)
            for i in range(len(sugar_species(schema))):
                assert -float(d[s_slice.start + i]) == rates[i], (
                    f"uptake's dS and the shared helper disagree at slot {i} for "
                    f"S={sugars}, X={x}: the yield's denominator has drifted from the flux"
                )


def test_a_negative_sugar_excursion_cannot_make_the_producer_create_acid(beer_params):
    """The clamp moved into the helper with the arithmetic it guards — check it came along.

    A small negative ``S`` excursion is normal solver behaviour. In the uptake Process that
    clamp stops a sign flip from *creating* sugar; here it must additionally stop the producer
    booking a negative rate, which would drain an acid pool that nothing consumes.
    """
    schema = beer_schema()
    params = beer_params.resolve()
    y = _beer_state(schema, S=[-1e-9, -1e-12, -1e-10])
    assert fermentative_uptake_rates(y, schema, params) == [0.0, 0.0, 0.0]
    d = OrganicAcidExcretion().derivatives(0.0, y, schema, params)
    for spec in ORGANIC_ACID_SPECS:
        assert float(d[schema.slice(spec.slot)][0]) == 0.0


# ======================================================================================
# 2. Isolability — wine untouched, un-anchored beer untouched
# ======================================================================================


def test_wine_never_wires_or_reads_the_beer_acid_producer():
    """Beer-only, structurally — not "wine's numbers happen not to move".

    Two distinct claims, because they fail independently: the Process is absent from wine's
    set (so no wine trajectory can carry it), AND no active wine Process declares the yields
    in ``reads`` (so no wine ENSEMBLE draws them — ``reads`` has two masters since D-160, and
    an extra declared name shifts a wine ensemble's draw sequence even when nominals do not
    move).
    """
    wine = get_medium("wine").build_process_set(strict=True)
    assert OrganicAcidExcretion.name not in {p.name for p in wine.active}
    assert AceticAcidOverflow.name not in {p.name for p in wine.active}

    yields = {spec.yield_param for spec in ORGANIC_ACID_SPECS} | {ACETIC_YIELD_PARAM}
    declared = {name for p in wine.active for name in p.reads}
    assert not (yields & declared), (
        f"a wine Process declares beer's organic-acid yields {sorted(yields & declared)}; "
        "that would put them in wine's sampled set and shift its draw sequence"
    )
    # SCOPED TO THE YIELD ON PURPOSE, and this is the one place the D-183 Process differs in
    # kind from D-180's. `AceticAcidOverflow.reads` also names `mu_max`/`K_s`/`K_n`, which wine's
    # own growth Process declares and MUST keep declaring — so the disjointness claim above
    # cannot be "none of its reads appear in wine" without being false for a correct model.
    # What has to be beer-only is the yield, because that is the only name whose presence in
    # wine's sampled set would be a defect rather than a shared dependency.
    assert set(AceticAcidOverflow.reads) & declared == {"mu_max", "K_s", "K_n"}, (
        "the growth-linked producer's overlap with wine's declared reads changed; it must be "
        "exactly the three growth constants it recomputes the growth rate from"
    )


def test_an_unanchored_beer_keeps_the_producer_disabled_and_its_acid_slots_empty():
    """The correctness gate, not a tidiness one — and the reason it is a gate.

    Without ``initial_ph`` every acid slot and the counter-cation are 0: an EMPTY charge
    balance, which ``charge_balance_is_populated`` reports as "this beverage does not claim a
    pH" so the aging rate laws hold their factor at 1 instead of aging against water's 7.0.
    An ENABLED producer would fill those slots from sugar as the ferment ran, so a beer that
    supplied no pH would ACQUIRE a populated charge balance mid-run — the D-179 defect wearing
    a producer's hat, arriving silently on a green suite.
    """
    unanchored = {k: v for k, v in TYRELL_SCENARIO.items() if k != "initial_ph"}
    compiled, res = _run(unanchored)

    assert OrganicAcidExcretion.name not in {p.name for p in compiled.process_set.active}
    assert AceticAcidOverflow.name not in {p.name for p in compiled.process_set.active}
    for slot in (*(spec.slot for spec in ORGANIC_ACID_SPECS), ACETIC_SLOT):
        assert float(res.series(slot)[-1]) == 0.0, (
            f"{slot} grew on a beer that supplied no initial_ph"
        )
    # Sugar was fully fermented all the same — the run is a normal beer, just pH-less.
    assert float(res.series("S")[:, -1].sum()) < 0.1
    assert not charge_balance_is_populated(res.y[:, 0], compiled.schema)


def test_the_producer_honours_its_touches_contract():
    """``strict=True`` is what turns ``touches`` from documentation into a contract."""
    pset = get_medium("beer").build_process_set(strict=True)
    assert OrganicAcidExcretion.name in pset
    schema = beer_schema()
    data = default_data_dir()
    params = load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    ).resolve()
    d = OrganicAcidExcretion().derivatives(0.0, _beer_state(schema), schema, params)
    touched = {n for n in schema.names if float(abs(d[schema.slice(n)]).sum()) > 0.0}
    assert touched <= set(OrganicAcidExcretion.touches)
    # And it really does move all three acids plus sugar — an inert Process would pass the
    # subset check above vacuously.
    assert touched == {*(spec.slot for spec in ORGANIC_ACID_SPECS), "S"}
    assert ACETIC_SLOT not in touched, (
        "acetic left ORGANIC_ACID_SPECS at D-183; if the flux-linked producer moves it again "
        "the acid is being made twice, once on each rate law"
    )


def test_the_growth_linked_producer_honours_its_touches_contract():
    """The D-183 sibling — and the ``X``/``N`` exclusion is the claim worth pinning.

    ``AceticAcidOverflow`` *reads* the growth rate and must not *contribute* to it. Declaring
    ``X``/``N`` in ``touches`` would let a future edit add biomass or consume nitrogen here
    while ``strict=True`` stayed green, and the second copy of growth's stoichiometry is exactly
    the drift D-106 caught one module over.
    """
    pset = get_medium("beer").build_process_set(strict=True)
    assert AceticAcidOverflow.name in pset
    schema = beer_schema()
    data = default_data_dir()
    params = load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    ).resolve()
    d = AceticAcidOverflow().derivatives(0.0, _beer_state(schema), schema, params)
    touched = {n for n in schema.names if float(abs(d[schema.slice(n)]).sum()) > 0.0}
    assert touched <= set(AceticAcidOverflow.touches)
    assert touched == {ACETIC_SLOT, "S"}
    assert not ({"X", "N"} & set(AceticAcidOverflow.touches))


# ======================================================================================
# 3. The Byp double-count, and the constraint that prevents it
# ======================================================================================


def test_beer_byproduct_yield_stays_zero_or_succinic_double_counts(beer_params):
    """Beer's ``Y_byproduct_sugar = 0`` was incidental before D-180. It is load-bearing now.

    Beer carries its own ``succinic`` state slot (D-179) *and* the charge balance reads ``Byp``
    as succinic-equivalent (``BYP_AS_SUCCINIC``, the D-18 include-by-reading coupling). While
    beer produced no ``Byp`` those two could not overlap. Give beer a non-zero
    ``Y_byproduct_sugar`` now and the same acid is counted twice — once in its own slot, once
    through the lump — in the pH solve AND on the carbon ledger.

    The zero is asserted, and then the double-count is DEMONSTRATED rather than described, so
    that a future edit to the yield fails here with the reason attached rather than shifting a
    pH by a few hundredths somewhere downstream.
    """
    assert beer_params["Y_byproduct_sugar"].value == 0.0
    assert beer_params["Y_glycerol_sugar"].value == 0.0

    schema = beer_schema()
    params = beer_params.resolve()
    pka = acidbase.build_pka_map(params)
    totals = {"succinic": 0.102 / acidbase.ALL_ACIDS["succinic"].molar_mass}
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 4.4)

    honest = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    # The same succinic, half of it also booked through the Byp lump.
    doubled = acidbase.solve_ph(
        {"succinic": totals["succinic"] * 0.5}, cation, totals["succinic"] * 0.5, 0.0, pka
    )
    assert honest == pytest.approx(doubled, abs=1e-9), (
        "sanity: splitting succinic between its slot and Byp must be pH-neutral"
    )
    # ...which is precisely why ADDING a Byp on top of a full succinic slot is not.
    inflated = acidbase.solve_ph(totals, cation, totals["succinic"], 0.0, pka)
    assert inflated < honest - 0.05, (
        "a non-zero beer Y_byproduct_sugar would count succinic twice in the charge balance"
    )
    assert schema.slice("succinic") is not None and "Byp" in schema


def test_carbon_closes_with_the_producer_running():
    """The acids' carbon comes OUT of sugar, so the ledger must close to machine precision.

    This is what makes ``acetic``/``succinic`` weighted terms in ``total_carbon`` (D-180) and
    ``acetic_acid`` an entry in the chemistry registries — before D-180 the slots were inert
    dosed constants, invisible to a conservation *difference*, so neither was needed.
    """
    from fermentation.validation.conservation import total_carbon

    compiled, res = _run(dict(TYRELL_SCENARIO))
    weights = total_carbon(
        compiled.schema,
        biomass_carbon_fraction=compiled.parameters["biomass_C_fraction"].value,
    )
    start, end = weights(res.y[:, 0]), weights(res.y[:, -1])
    assert end == pytest.approx(start, rel=1e-9), (
        f"carbon drifted {end - start:.6e} g C/L with the organic-acid producer active"
    )


# ======================================================================================
# 4. Acceptance — the matched Tyrell comparison
# ======================================================================================


def test_produced_acids_land_in_tyrells_measured_beer_bands():
    """A WIRING check, not a validation — and the distinction is worth naming.

    The yields were derived as ``(beer − wort)/ΔS`` from these very curves, so a fully
    attenuated ferment of this wort reproduces these levels close to by construction. What it
    is NOT free of is the wiring: a wrong divisor, a mg/L-vs-g/L slip, a yield read against the
    wrong flux, or a modifier that scales uptake but not the producer would all show up here.
    That is why it is worth having and why it is not evidence the model is right.
    """
    _, res = _run(dict(TYRELL_SCENARIO))
    for slot, (lo, hi) in TYRELL_BEER_PPM.items():
        got = float(res.series(slot)[-1]) * 1000.0
        assert lo <= got <= hi, (
            f"{slot} finished at {got:.1f} mg/L, outside Tyrell's measured {lo}-{hi} mg/L"
        )


def _daily_ppm(res, slot: str, days: int = 7) -> dict[int, float]:
    """A slot's value at each whole day [mg/L], read on a FIXED hourly grid.

    Never index the solver's own output for a value at a named time: BDF places its steps
    adaptively and the nearest one can sit 20+ minutes off the hour, which on a steep limb is a
    double-digit percentage error that is stable enough to look measured
    [[feedback-read-a-fast-curve-on-a-fixed-grid]].
    """
    grid = np.linspace(0.0, days * 24.0, days * 24 + 1)
    series = np.interp(grid, res.t, np.asarray(res.series(slot), dtype=float)) * 1000.0
    return {d: float(series[d * 24]) for d in range(days + 1)}


@pytest.mark.xfail(strict=True, reason="D-215: the model ferments this wort ~2.8x too slowly")
def test_the_model_ferments_tyrells_wort_on_tyrells_schedule():
    """The extract panel of Fig. 4, scored against the model for the first time (D-215).

    **Why this is an xfail and not a pin.** It states the thing that is true of the source and
    false of the model, so a correct fix turns it GREEN and nothing has to be deleted — the D-208
    idiom. A plain assert on the model's *current* fraction would encode the defect instead, which
    is exactly why D-207 shipped the pH course as data no assert read.

    The gap is not subtle. Tyrell's wort is **59.4 %** fermented by day 2; this engine books
    **~21 %**, and does not reach dryness until about day 10 where their extract curve is flat by
    day 5. The tolerance below is ±0.10 of the fermentable — roughly five times any plausible read
    error on the extract panel — so this cannot be argued down to a transcription quibble.

    **This is upstream of every acid course.** The acids are produced as ``Y · ΔS``, so an extract
    curve that is too gradual makes every flux-linked acid too gradual with it. Fix this and the
    acid courses move without any rate law changing; fix the acid rate laws while this stands and
    they will be fitted to compensate for it [[feedback-a-margin-can-be-borrowed-from-a-defect]].

    **SCOPE: this is measured on ONE scenario** — Tyrell's wort, 15 °C, YAN 200, pitch 1.0 g/L. It
    is a statement about how this engine ferments *that* wort, NOT an established property of beer's
    kinetics in general; nothing here separates a scenario artefact (their EBC-tube trial, their
    pitch rate) from a general rate defect. D-211 measured total attenuation at 6.08 d inside §2.2's
    5-7 d window, which is consistent with either. Do not generalise this without measuring it.

    **What it does NOT claim.** It does not say the day-1 pH miss is caused by this. D-215 probed
    that directly and could not attribute it: re-scoring the pH matched on extent rather than on
    the calendar OVERSHOOTS, moving day 1 from 0.070 too alkaline to 0.063 too acidic. Both a
    too-slow ferment and a too-generous acid yield are live, and they partly cancel on a
    time-matched read. Do not cite this test as the diagnosis for D-211 §9's brief.
    """
    compiled, res = _run(dict(TYRELL_SCENARIO), days=7.0)
    grid = np.linspace(0.0, 7 * 24.0, 7 * 24 + 1)
    sugar = np.asarray(res.y, dtype=float)[compiled.schema.slice("S"), :]
    total = np.vstack([np.interp(grid, res.t, row) for row in sugar]).sum(axis=0)
    initial = float(total[0])
    assert initial == pytest.approx(TYRELL_SUGAR_GPL, rel=1e-9), (
        "the scenario is not carrying Tyrell's fermentable divisor, so the fractions below "
        "would be measured against the wrong wort"
    )
    for day, measured in TYRELL_FLUX_FRACTION.items():
        if day > 7:
            continue
        got = (initial - float(total[day * 24])) / initial
        assert got == pytest.approx(measured, abs=0.10), (
            f"day {day}: the model has fermented {got:.1%} of the wort against Tyrell's measured "
            f"{measured:.1%}. Their extract panel is the driver of every flux-linked acid course "
            "in this file; see TYRELL_FLUX_FRACTION."
        )


@pytest.mark.xfail(strict=True, reason="D-215: lactic/malic/succinic courses are mistimed")
def test_the_three_flux_linked_acid_courses_are_mistimed():
    """The three ``Y·ΔS`` acids against the days they claim to predict (D-215).

    D-183 scored acetic against Fig. 13's interior and the flux-linked law lost. The other three
    were never scored the same way — only their day-7 endpoints were, and those agree close to by
    construction because the yields are ``(day7 − day0)/ΔS`` off these very curves
    (``test_produced_acids_land_in_tyrells_measured_beer_bands`` says so in its own docstring).

    Tolerance is **three times** the ±2 ppm figure read tolerance, so a failure here is a shape
    disagreement and not a quarrel about pixels.

    **The three do not fail the same way, and that is the finding**, not a nuisance:

    * **succinic** is the closest to the shipped law and acts as the control — measured 45.9 % of
      its rise done by day 2 against a modelled 20.5 %, i.e. the model is *late*;
    * **malic** errs the other way, ~0 % measured against 20.5 % modelled, i.e. the model is
      *early*;
    * **lactic** sits between them and is nearly right at day 2, then falls behind.

    So there is **no single timing correction that helps all three**, and the two large errors have
    OPPOSITE sign — which is why the shared flux-linked shape survived four beats of endpoint
    checks. A fix that made all three match by moving one rate law would be fitting a compromise,
    not a mechanism [[feedback-a-hit-can-be-two-errors-cancelling]].
    """
    _, res = _run(dict(TYRELL_SCENARIO), days=7.0)
    tol = 3.0 * TYRELL_ACID_COURSE_READ_TOL
    got = {slot: _daily_ppm(res, slot) for slot in TYRELL_ACID_COURSE_PPM}

    # DAY 2 FIRST, and across all three — because day 2 is where the opposing-sign finding
    # lives and this test is what carries it. Iterating acid-by-acid instead would die on
    # lactic's day 4 (a real but different miss) and the RED would never name the claim the
    # record leads with [[feedback-grep-finds-claims-not-guards]].
    for slot, course in TYRELL_ACID_COURSE_PPM.items():
        assert got[slot][2] == pytest.approx(course[2], abs=tol), (
            f"{slot} day 2: model {got[slot][2]:.1f} vs Tyrell's measured {course[2]:.1f} mg/L. "
            "Day 2 is where the three errors have OPPOSING signs — succinic late, malic early, "
            f"lactic nearly right (tolerance ±{tol:.0f}, i.e. 3× the figure read tolerance)"
        )
    for slot, course in TYRELL_ACID_COURSE_PPM.items():
        for day, measured in course.items():
            assert got[slot][day] == pytest.approx(measured, abs=tol), (
                f"{slot} day {day}: model {got[slot][day]:.1f} vs Tyrell's measured "
                f"{measured:.1f} mg/L (tolerance ±{tol:.0f}, i.e. 3× the read tolerance)"
            )


def test_the_acid_courses_are_anchored_to_numbers_recorded_before_them():
    """The transcription's own check — and it is NOT an xfail, because it must hold today.

    ``TYRELL_ACID_COURSE_PPM`` was read off the figures at D-215, years of beats after the seeds
    and yields derived from the same figures were recorded. If the new reads did not reproduce
    those older numbers, the transcription would be worthless and every conclusion resting on it
    would be a fresh error rather than a measurement
    [[feedback-a-hit-can-be-two-errors-cancelling]].

    Day 0 against the shipped wort seed, and day 7 against the shipped measured beer band.
    """
    data = default_data_dir()
    params = load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    ).resolve()
    for slot, course in TYRELL_ACID_COURSE_PPM.items():
        seed_ppm = params[f"{slot}_typical_wort"] * 1000.0
        assert course[0] == pytest.approx(seed_ppm, abs=TYRELL_ACID_COURSE_READ_TOL), (
            f"{slot}'s day-0 read {course[0]:.2f} disagrees with the shipped wort seed "
            f"{seed_ppm:.2f} mg/L by more than the figure read tolerance"
        )
        lo, hi = TYRELL_BEER_PPM[slot]
        assert lo - TYRELL_ACID_COURSE_READ_TOL <= course[7] <= hi + TYRELL_ACID_COURSE_READ_TOL, (
            f"{slot}'s day-7 four-strain mean {course[7]:.2f} falls outside the recorded measured "
            f"band {lo}-{hi} mg/L even allowing the read tolerance"
        )


def test_the_yield_is_a_yield_at_every_temperature():
    """The D-32 correctness coupling in the form it takes here (the uptake Arrhenius).

    ``OrganicAcidExcretion`` recomputes the *unmodified* uptake rates, while uptake's own
    contribution is scaled by the uptake Arrhenius. Unless the producer carries the identical
    factor, a cold beer would ferment slowly and make acid at the warm rate — g-acid-per-g-sugar
    would stop being constant, silently. Two temperatures, same wort, both run to dryness: the
    finished acid levels must agree.
    """
    _, cold = _run(dict(TYRELL_SCENARIO), celsius=10.0, days=40.0)
    _, warm = _run(dict(TYRELL_SCENARIO), celsius=22.0, days=40.0)
    for spec in ORGANIC_ACID_SPECS:
        c = float(cold.series(spec.slot)[-1])
        w = float(warm.series(spec.slot)[-1])
        assert c == pytest.approx(w, rel=2e-3), (
            f"{spec.slot} differs by {abs(c - w) / w:.2%} between a 10 °C and a 22 °C ferment "
            "of the same wort; the producer is not carrying the uptake Arrhenius factor"
        )
    # Acetic rides the SAME coupling against a DIFFERENT modifier (D-183). Its yield is per gram
    # of biomass, and biomass is nitrogen-capped rather than temperature-capped, so the finished
    # level must be temperature-invariant for the same reason — but only if the Process is named
    # by the GROWTH Arrhenius. Left on the uptake one it would still integrate and still pass
    # every other test in this file; this assert is what makes that mis-wiring loud.
    cold_acetic = float(cold.series(ACETIC_SLOT)[-1])
    warm_acetic = float(warm.series(ACETIC_SLOT)[-1])
    assert cold_acetic == pytest.approx(warm_acetic, rel=2e-3), (
        f"acetic differs by {abs(cold_acetic - warm_acetic) / warm_acetic:.2%} between a 10 °C "
        "and a 22 °C ferment; the growth-linked producer is not carrying the GROWTH Arrhenius "
        "factor (it is a `for_growth` extra target, not a `for_uptake` one)"
    )


def test_a_flux_linked_acetic_yield_puts_the_acid_where_the_source_says_it_is_not():
    """THE COUNTERFACTUAL THAT RETIRED A RATE LAW (D-183) — recomputed, never hard-coded.

    ``Y_acetic_sugar_beer`` is read by no Process since D-183. It is kept for the same reason
    group 1's finished-beer levels are: it is the only thing that keeps the retirement
    falsifiable instead of merely written down.

    The retirement rests on mapping Tyrell's Fig 13 onto their Fig 4. By day 1 that wort is
    **15 % attenuated** while acetic has already made **77 %** of its whole rise. So a yield on
    the sugar flux — which is what ``Y_acetic_sugar_beer`` is, a measured day7−day0 difference
    divided by a sugar divisor — must put the acid far too late. This recomputes exactly how
    late, from the flux the shipped run actually booked (recovered through *lactic*, which still
    rides that flux), rather than trusting a number in a comment.

    (D-183 and this docstring both said **86 %**, which is a units slip corrected at D-211: 86.00
    is acetic's day-1 rise in **mg/L**, and its share of the rise to the measured peak is
    111.25 mg/L, i.e. 77.3 %. The disproportion the retirement rests on is unaffected — 77 %
    against 15 % of the flux — so no shipped value moves. It is corrected because the number is
    now asserted below rather than quoted.)

    It is a claim about SHAPE, not size: the same yield reproduces the day-7 endpoint by
    construction, which is precisely why the error was invisible until the figure interiors were
    read.

    **D-211 narrowed this margin and the reason is worth stating.** When D-183 chose the
    growth-linked producer, growth in this model finished inside ~20 h, so a producer tied to it
    delivered ~100 % of its acid by day 1 and appeared to explain Tyrell's early rise. D-211
    measured the growth rate against Tyrell's OWN cell counts and found it 2.88x too fast: the real
    crop has made only ~36 % of its growth by day 1. So acetic rises faster than growth AND faster
    than flux, and the old rate was flattering the attribution. The retirement still holds —
    growth-linked beats flux-linked by every measure below — but it now wins while booking
    **0.360** of acetic's rise by day 1 against a measured **0.773**, a **2.15x** shortfall it
    previously had no way to show. The residual is pinned at the foot of this test.
    """
    compiled, res = _run(dict(TYRELL_SCENARIO))
    params = compiled.parameters.resolve()
    t_h = np.asarray(res.t, dtype=float)
    # The fermentative flux integral, recovered from an acid that still rides it.
    lactic = np.asarray(res.series("lactic"), dtype=float)
    phi = (lactic - lactic[0]) / params["Y_lactic_sugar_beer"]
    seed = float(res.series(ACETIC_SLOT)[0]) * 1000.0

    for day in (1, 2):
        flux_linked = (
            seed + params["Y_acetic_sugar_beer"] * float(np.interp(day * 24.0, t_h, phi)) * 1000.0
        )
        shipped = float(np.interp(day * 24.0, t_h, np.asarray(res.series(ACETIC_SLOT)))) * 1000.0
        measured = TYRELL_ACETIC_MEAN_PPM[day]
        assert flux_linked < measured - 60.0, (
            f"day {day}: the retired flux-linked yield gives {flux_linked:.1f} mg/L against "
            f"Tyrell's measured {measured:.2f}. If this gap has closed, the sugar-uptake "
            "kinetics moved, and the whole basis for the D-183 rate-law change needs re-reading"
        )
        assert abs(shipped - measured) < abs(flux_linked - measured), (
            f"day {day}: the shipped growth-linked producer ({shipped:.1f} mg/L) is no closer to "
            f"the measured {measured:.2f} than the retired flux-linked one ({flux_linked:.1f}); "
            "the rate-law change bought nothing"
        )

    # And the size claim: the shipped form halves the shape error over Tyrell's measured days.
    days = [d for d in TYRELL_ACETIC_MEAN_PPM if d >= 1]
    shipped_rmse = float(
        np.sqrt(
            np.mean(
                [
                    (
                        float(np.interp(d * 24.0, t_h, np.asarray(res.series(ACETIC_SLOT))))
                        * 1000.0
                        - TYRELL_ACETIC_MEAN_PPM[d]
                    )
                    ** 2
                    for d in days
                ]
            )
        )
    )
    flux_rmse = float(
        np.sqrt(
            np.mean(
                [
                    (
                        seed
                        + params["Y_acetic_sugar_beer"]
                        * float(np.interp(d * 24.0, t_h, phi))
                        * 1000.0
                        - TYRELL_ACETIC_MEAN_PPM[d]
                    )
                    ** 2
                    for d in days
                ]
            )
        )
    )
    assert shipped_rmse < flux_rmse, (
        f"the growth-linked form scores {shipped_rmse:.1f} ppm RMSE against Tyrell's days 1-7 "
        f"and the retired flux-linked one {flux_rmse:.1f}; D-183 measured 32.5 vs 61.6 and "
        "D-211's slower growth moved them to 40.7 vs 65.3. The RETIREMENT is what this pins — "
        "if the shipped form ever scores worse, the rate-law change has stopped paying"
    )
    # The ratio is pinned SEPARATELY from the ordering, because the two say different things and
    # a single threshold would let one hide the other: the ordering is D-183's verdict, the ratio
    # is how much of it survives. D-183 measured 0.528; D-211's re-derived growth rate cost it
    # ~0.10, which is the price this beat paid on a neighbouring claim and is recorded, not tuned.
    assert shipped_rmse / flux_rmse == pytest.approx(0.624, abs=0.03), (
        f"the shape-error ratio is {shipped_rmse / flux_rmse:.3f}; D-183 measured 0.528 before "
        "D-211 slowed growth to Tyrell's measured rate and 0.624 after. A move here means the "
        "growth path changed again"
    )
    # THE RESIDUAL D-211 EXPOSED, pinned where it was found. Acetic makes 86 % of its measured
    # rise by day 1; the model's growth-linked producer can only make as much as the crop has
    # grown, which the cell counts put at ~34 %. Neither driver in the file explains the early
    # rise — the old growth rate merely concealed that by finishing inside 20 h.
    shipped_rise = {
        d: float(np.interp(d * 24.0, t_h, np.asarray(res.series(ACETIC_SLOT)))) * 1000.0 - seed
        for d in (1, 7)
    }
    model_share_day1 = shipped_rise[1] / shipped_rise[7]
    measured_share_day1 = (TYRELL_ACETIC_MEAN_PPM[1] - TYRELL_ACETIC_MEAN_PPM[0]) / (
        max(TYRELL_ACETIC_MEAN_PPM.values()) - TYRELL_ACETIC_MEAN_PPM[0]
    )
    assert measured_share_day1 == pytest.approx(0.774, abs=0.01), (
        "the measured day-1 share of acetic's rise is a transcription, not a model output"
    )
    assert measured_share_day1 / model_share_day1 == pytest.approx(2.15, abs=0.15), (
        f"the model books {model_share_day1:.3f} of acetic's rise by day 1 against Tyrell's "
        f"{measured_share_day1:.3f}, a {measured_share_day1 / model_share_day1:.2f}x shortfall; "
        "D-211 measured 2.15x once the growth rate was corrected. If it has closed, a producer "
        "faster than growth has been built and D-183's growth-linked choice should be re-read"
    )
    # THE HONEST CEILING, asserted so it cannot quietly be read as a transient: the shipped
    # curve is MONOTONE. It reaches its endpoint early and holds; it does not peak and fall.
    # Delivering the peak needs the re-assimilation half, which D-183 measured and REFUSED.
    shipped_series = np.asarray(res.series(ACETIC_SLOT), dtype=float)
    # NON-DECREASING, not "ends at its maximum". The weaker form passes for any curve that
    # overshoots and comes back to exactly the endpoint — which is precisely the re-assimilated
    # transient this asserts the absence of — and it passes for the retired producer too.
    #
    # THE TOLERANCE IS A MEASURED NOISE FLOOR, not a round number: this pool is mathematically
    # non-decreasing and the shipped BDF run still takes 38 negative steps out of 199, the worst
    # −2.0e−10 g/L (−1.7e−9 relative to the pool). 1e−7 sits ~58× above that and ~1000× below
    # the step a real re-assimilation term produces, so it separates the two rather than
    # splitting the difference ([[feedback-pin-tolerance-vs-solver-tolerance]]).
    steps = np.diff(shipped_series)
    assert float(steps.min()) >= -1e-7 * float(shipped_series.max()), (
        f"acetic decreases somewhere (worst step {float(steps.min()):.3e} g/L) — something "
        "restored the re-assimilation half that D-183 refused on duration-dependence, and the "
        "D-182 pH headline is now a function of the run's duration. Read the D-183 record "
        "before keeping this"
    )
    assert shipped_rmse > 25.0, (
        "the shape error against Tyrell's own days collapsed; a MONOTONE curve cannot fit a "
        "rise-then-fall well, so this passing would mean the comparison stopped being real"
    )


def test_citrate_is_seeded_and_stays_inert():
    """Sourced three ways, so pin it three-ways-worth.

    Tyrell's §1.5 ("citric … basically depends on concentrations of wort"), their Table 2
    ("final concentration mostly determined by wort concentration") and their Fig. 12 (four
    strains scattering ±20 ppm with no trend) all say beer's citrate is malt-derived. It is the
    one shipped beer acid with a wort seed and NO yield, and it must stay that way — giving it
    one would be a mechanism claim this file's own provenance contradicts.
    """
    assert "citrate" not in {spec.slot for spec in ORGANIC_ACID_SPECS}
    _, res = _run(dict(TYRELL_SCENARIO))
    series = res.series("citrate")
    assert float(series[0]) == pytest.approx(0.205, abs=1e-12)
    assert float(series[-1]) == pytest.approx(float(series[0]), abs=1e-12)


def test_the_predicted_ph_drop_over_the_joint_yield_and_pka_band(beer_params):
    """THE FREE PREDICTION — over the JOINT band of every quantity the sampler actually draws.

    Nothing in ``beer_acids.yaml`` is fitted to pH, so comparing the modelled trajectory with
    Tyrell's Fig. 4 is a genuine external test.

    **THE COMPARISON RUNS IN THE DEGASSED FRAME SINCE D-208, AND THAT IS THE WHOLE HEADLINE.**
    Tyrell cite *"MEBAK, Band II, 4th edition … 2.14 (pH)"*; MEBAK's pH method is *"pH (EBC)"*
    and instructs that carbonated beverages be decarbonated before measurement, and
    Analytica-EBC 9.35 states its scope as *"the determination of pH at 20 °C of DECARBONATED
    beer using a pH meter"*. So the published number is a **decarbonated** reading, while
    ``ph_of_state`` reports the pH **inside the vessel** — the two differ by 0.29 pH at day 7
    and 0.50 at day 1 (D-207 §5). For four beats this test scored the in-vessel pH against the
    decarbonated measurement, and the agreement it reported was the frame difference:

    * **at nominal, across ``pKa_peptide_buffer``'s band, DEGASSED — the comparison with
      Tyrell: 43.2-62.9 %** of the measured 0.81 pH drop.
    * **the same members IN-VESSEL: 77.8-97.3 %** — pinned here as a property of the model, and
      **not** a comparison with anything published. This is the number D-180 → D-183 quote
      (63-92 % at D-180, 42.7-62.2 % at D-181, 77.6-97.0 % at D-182), so their history is a
      history of the in-vessel quantity.
    * **over the JOINT band — TEN drawn quantities, not one: 8.3-82.7 % degassed,
      63.8-109.7 % in-vessel.** In the measured frame **nothing in the band reaches the
      measurement**, which is where D-181 left this axis before the carbonic term was built.

    The degassed frame is a **bound, not a point**: a real decarbonation leaves a residue, so the
    honest object is the one-parameter family "sample retained fraction ``s`` of saturation",
    whose ends are the two columns above. D-208 walked it and **no member reproduces the
    measured course** — the ``s`` that fits day 1 (0.150) leaves day 7 at 5.17 against a measured
    4.804-4.916, and the day-7 pH sits above the four-strain envelope for every ``s`` in [0, 1].
    That is why the acceptance claim now lives in the ``xfail`` below rather than in a floor here.

    Two scopes, and conflating them is the other trap this test exists to avoid:

    **The joint band has grown three times, and the middle one is the instructive one.**
    D-180's amendment added the four ``Y_*_sugar_beer`` after the first version asserted "the
    model must fall short" on the nominal alone — a constraint verified at a POINT where the
    sampler reads a BAND. D-181 then shipped the three FLOORS as a third dimension and
    **committed the identical mistake one level out**, holding ``pKa_oxalic_2``,
    ``pKa_pyruvic`` and the three wort SEEDS at nominal while calling the result band-wide.
    D-182 adds the three the carbonic term brought (``pKa_carbonic_1``, ``H_co2_beverage``,
    ``vant_hoff_co2_solubility``) in the same beat that ships them, rather than in the record
    that documents the beat — which is what the two previous instances failed to do. All nine
    are drawn; all nine are varied here, and the corner COUNT is asserted so a future
    dimension cannot be added to the registry without being added here.

    **What each dimension is worth, measured** — IN-VESSEL (all-nominal fraction 0.8559, whole
    band each). In the degassed frame the last three are **algebraically inert**: the term they
    parameterise is multiplied by zero, so the degassed joint band holds 3⁷ distinct values in
    3¹⁰ corners. That is inertness by construction, not the measured kind
    [[feedback-a-nominal-on-a-band-edge-is-not-inertness]], and the corner count is asserted in
    both frames so a new dimension still has to be added here:
    the peptide pKa moves it 0.776-0.970; the yields and floors comparably; the three seeds
    ~0.015; ``pKa_oxalic_2``/``pKa_pyruvic`` **0.0003**; and of D-182's three, the carbonic
    pKa is worth ~0.01 and the two solubility parameters ~0.01 between them — the CO2 term's
    SIZE is consequential but its band is not, because both edges sit within 10 % of a
    nominal that is itself a printed in-beer measurement.

    **A CORNER REACHES THE MEASUREMENT AGAIN — AND SINCE D-208 ONLY IN THE WRONG FRAME.** D-180
    had one at 104.5 %, and that reach belonged to the falling acids' absence rather than to the
    model; D-181 removed it and wrote that "a future change that makes one reach again is a
    signal to find out which omitted term arrived"
    [[feedback-a-margin-is-a-claim-about-what-holds-it-open]]. One did: dissolved CO2, the last
    of the two terms D-180 named, built at D-182 and pre-registered at 76-104 %. D-208 then
    settled which frame the measurement lives in, and **every one of those reaching corners is
    in-vessel**: in the degassed frame the joint high corner is 82.7 %, so nothing reaches. The
    signal D-181 asked for was real but it pointed at a **comparison** defect, not at an
    arrived term — the term is right where it is and the scoring was not
    [[feedback-a-summary-statistic-is-not-the-curve]].

    **The arm RE-ANCHORS the cation per member, and getting that wrong is the other trap.**
    Re-reading the shipped trajectory's pH at a different pKa while holding the
    ``cation_charge`` the compile back-solved at the NOMINAL pKa reports 72-80 % and moves the
    START pH off the 5.65 the scenario supplied — which no ensemble member ever does. A member
    draws its pKa and the compile back-solves ITS cation to hit ``initial_ph``, so every member
    starts at 5.65 and they differ only in where they end.
    """
    compiled, res = _run(dict(TYRELL_SCENARIO))
    params = compiled.parameters.resolve()
    pka_map = acidbase.build_pka_map(params)
    slots = (
        "lactic",
        "acetic",
        "citrate",
        "malic",
        "succinic",
        "pyruvic",
        "formic",
        "oxalic",
        "peptide_buffer",
    )
    molar = {s: acidbase.ALL_ACIDS[s].molar_mass for s in slots}
    start_molar = {s: float(res.series(s)[0]) / molar[s] for s in slots}
    measured_drop = TYRELL_WORT_PH - sum(TYRELL_BEER_PH) / 2.0

    # The fermentative flux the shipped run actually booked its yields against — recovered
    # from what it produced, so this cannot drift from the run. Recovered through LACTIC since
    # D-183: acetic no longer rides this flux, and dividing its production by a yield it does
    # not use would silently mis-scale every produced acid in this arm.
    flux = (float(res.series("lactic")[-1]) - float(res.series("lactic")[0])) / params[
        "Y_lactic_sugar_beer"
    ]
    # Acetic's own denominator (D-183): the gross biomass this run formed, N being growth's only
    # consumer in beer. Recovered from the run for the same reason the flux is.
    biomass_formed = (float(res.series("N")[0]) - max(float(res.series("N")[-1]), 0.0)) / params[
        "biomass_N_fraction"
    ]
    # D-209: the nitrogen the run started with and ended on, in g N/L. Recovered from the run for
    # the same reason as the two above — the charge this arm removes must be the charge the
    # solver's own nitrogen carried, not a nominal 200 mg/L that a scenario override could
    # silently diverge from.
    nitrogen_gpl_start = float(res.series("N")[0])
    nitrogen_gpl_end = max(float(res.series("N")[-1]), 0.0)

    def edge(param: str, pick: str) -> float:
        p = beer_params[param]
        return float({"lo": p.uncertainty.low, "nom": p.value, "hi": p.uncertainty.high}[pick])

    def fraction(
        pka: float,
        pick: str,
        floor_pick: str,
        ox2: str,
        pyr: str,
        seed: str,
        carbonic: str,
        henry: str,
        vant_hoff: str,
        acetic_pick: str,
        nitro: str,
    ) -> tuple[float, float]:
        """``(degassed, in_vessel)`` fraction of the measured drop for one band member.

        Both frames come out of ONE member so they cannot drift apart: everything up to the
        endpoint solve is shared and only the carbonic term differs, which also keeps the cost
        of the second frame to one extra root-find per corner rather than a second sweep.
        """
        # Every pKa this beat added is drawn too (PH_SYSTEM_READS is the union, D-179), so
        # pinning them here would be the same point-vs-band mistake one level further out.
        member = {
            **pka_map,
            "peptide_buffer": (pka,),
            "oxalic": (pka_map["oxalic"][0], edge("pKa_oxalic_2", ox2)),
            "pyruvic": (edge("pKa_pyruvic", pyr),),
            "carbonic": (edge("pKa_carbonic_1", carbonic),),
        }
        # ...and so are the three wort seeds, which set how much charge there is to lose.
        seeded = dict(start_molar)
        for seed_slot, seed_param in WORT_SEED_PARAMS.items():
            seeded[seed_slot] = edge(seed_param, seed) / molar[seed_slot]
        # D-209: the wort's assimilable nitrogen is itself on the cation side, so the anchored
        # total splits into a frozen slot plus a term that MOVES as the yeast takes the nitrogen
        # up. The anchor is unaffected (the start uses the total, exactly as before), which is
        # why `start` still lands on the supplied wort pH bit-for-bit; what changes is the END,
        # where the pool is gone. Written as `cation` (the total at t=0) minus the charge the
        # consumed nitrogen carried, so the arm mirrors the shipped code's re-allocation rather
        # than re-deriving it.
        zbar = edge("nitrogen_uptake_charge_beer", nitro)
        nitrogen_charge_lost = zbar * (nitrogen_gpl_start - nitrogen_gpl_end) / M_NITROGEN
        cation = acidbase.solve_cation_charge(seeded, 0.0, 0.0, member, TYRELL_WORT_PH)
        start = acidbase.solve_ph(seeded, cation, 0.0, 0.0, member)
        assert start == pytest.approx(TYRELL_WORT_PH, abs=1e-6), (
            "every member must start at the supplied wort pH — that is what anchoring means"
        )
        cation_end = cation - nitrogen_charge_lost
        end = dict(seeded)
        for spec in ORGANIC_ACID_SPECS:
            # Produced acids build on the RUN's seed, not the varied one: a yield is a measured
            # production and does not depend on how much of a DIFFERENT acid the wort carried.
            end[spec.slot] = (
                start_molar[spec.slot] + (edge(spec.yield_param, pick) * flux) / molar[spec.slot]
            )
        # Acetic is the same construction against its OWN denominator (D-183). It is a separate
        # dimension rather than a fourth member of `pick` because its band is a spread over a
        # DIFFERENT quantity — g acid per g biomass, not per g sugar — so tying it to the other
        # three would assert a correlation between two unrelated strain rankings.
        end[ACETIC_SLOT] = (
            start_molar[ACETIC_SLOT]
            + (edge(ACETIC_YIELD_PARAM, acetic_pick) * biomass_formed) / molar[ACETIC_SLOT]
        )
        # The D-181 sinks land ON their floors: every k in the band clears the fall inside the
        # first days of this 14-day run, which the endpoint-insensitivity test pins separately.
        for sink in WORT_ACID_SINKS:
            end[sink.slot] = edge(sink.floor_param, floor_pick) / molar[sink.slot]
        # D-182's carbonic term, at THIS member's own solubility parameters. It is 0 in the
        # wort arm above and non-zero here, which is the whole shape of the term: the anchor
        # cannot absorb what is not present when the anchor is taken.
        member_params = {
            **params,
            "H_co2_beverage": edge("H_co2_beverage", henry),
            "vant_hoff_co2_solubility": edge("vant_hoff_co2_solubility", vant_hoff),
        }
        sat = acidbase.co2_saturation_gpl(float(res.series("T")[-1]), member_params)
        evolved = float(res.series("CO2")[-1])
        carbonic_molar = min(evolved, sat) / acidbase.CARBONIC_AS_CO2.molar_mass
        in_vessel = acidbase.solve_ph(end, cation_end, 0.0, carbonic_molar, member)
        # The frame Tyrell measured in (D-208): the same member with the sample decarbonated.
        degassed = acidbase.solve_ph(end, cation_end, 0.0, 0.0, member)
        return ((start - degassed) / measured_drop, (start - in_vessel) / measured_drop)

    pka_band = (
        beer_params["pKa_peptide_buffer"].uncertainty.low,
        beer_params["pKa_peptide_buffer"].value,
        beer_params["pKa_peptide_buffer"].uncertainty.high,
    )

    # Scope 1 — everything nominal but the peptide pKa, in BOTH frames. The comparison with
    # Tyrell is the degassed column; the in-vessel column is a property of the model and is
    # pinned beside it so the two can never be silently re-conflated (D-208).
    at_nominal = [
        fraction(pka, "nom", "nom", "nom", "nom", "nom", "nom", "nom", "nom", "nom", "nom")
        for pka in pka_band
    ]
    degassed_nominal = [f[0] for f in at_nominal]
    vessel_nominal = [f[1] for f in at_nominal]
    # BOTH edges pinned, not just the one a floor would guard: swapping a band edge for its
    # neighbour has passed a one-sided pin before [[feedback-pin-the-band-not-the-nominal]].
    assert min(degassed_nominal) == pytest.approx(0.914, abs=0.005), (
        f"the degassed prediction moved to {min(degassed_nominal):.1%} of Tyrell's measured "
        "0.81 pH drop at nominal yields; D-209 measures 91.4-127.1 % across the pKa band, "
        "against D-208's 43.2-62.9 % for the same arm before the nitrogen pool's charge was in "
        "the balance. This is the number that compares with a published beer pH"
    )
    assert max(degassed_nominal) == pytest.approx(1.271, abs=0.005), (
        f"the degassed prediction's high edge moved to {max(degassed_nominal):.1%}; D-209 "
        "measures 127.1 %, so at this arm's DAY-14 endpoint the high edge of the peptide pKa "
        "band now OVERSHOOTS the measured drop. Two things keep that honest and neither is a "
        "tuning knob: Tyrell's measurement stops at day 7 while this arm reads day 14, where "
        "the model is still producing acid (the day-7 comparison is the acceptance test below, "
        "101.7-107.1 %); and z-bar is DERIVED from published wort composition, never fitted"
    )
    assert min(vessel_nominal) == pytest.approx(1.068, abs=0.005), (
        f"the IN-VESSEL fraction moved to {min(vessel_nominal):.1%}; D-209 measures 106.8 % "
        "(D-183's 77.8 %, D-182's 77.6 %, and D-181's 42.7 % was this same model with no "
        "dissolved CO2 in its charge balance). Pinned as a model property: no published beer pH "
        "is measured in this frame, so a change here is a change to the vessel's chemistry, not "
        "to an agreement"
    )
    assert max(vessel_nominal) == pytest.approx(1.381, abs=0.005), (
        f"the IN-VESSEL fraction's high edge moved to {max(vessel_nominal):.1%}; D-209 measures "
        "138.1 % against D-183's 97.3 %. It now exceeds 1: the in-vessel pH carries a term the "
        "measurement excludes AND the nitrogen charge, so it is expected to run above the "
        "measured drop. The previous version of this pin asserted the opposite would be "
        "surprising — that reasoning belonged to a model with no nitrogen term"
    )

    # Scope 2 — the joint band the sampler can actually reach. No upper bound is asserted on the
    # IN-VESSEL band, because a corner of it legitimately exceeds the measured drop; what is
    # pinned there is the SPAN, so that a change which narrows or shifts it has to be looked at.
    # The DEGASSED band does carry one, since in the measured frame nothing reaches. D-181 adds the
    # three FLOORS as a third band dimension: their edges are named strains' own day-7 values,
    # so leaving them at nominal here would repeat the point-vs-band mistake one level out —
    # which is exactly what D-180's amendment had to correct in this very test.
    picks = ("lo", "nom", "hi")
    joint = [
        fraction(
            pka, pick, floor_pick, ox2, pyr, seed, carbonic, henry, vant_hoff, acetic_pick, nitro
        )
        for pka in pka_band
        for pick in picks
        for floor_pick in picks
        for ox2 in picks
        for pyr in picks
        for seed in picks
        for carbonic in picks
        for henry in picks
        for vant_hoff in picks
        for acetic_pick in picks
        for nitro in picks
    ]
    # ELEVEN dimensions since D-209, and the eleventh went in the SAME COMMIT that shipped the
    # band it varies. That ordering is the archive's fix for its most-repeated shape — a
    # constraint verified at a POINT where the sampler reads a BAND — which landed six times,
    # twice inside the record documenting the previous instance. At D-183 the tenth was
    # `Y_acetic_biomass_beer` replacing `Y_acetic_sugar_beer`, so the count did not move then;
    # `nitrogen_uptake_charge_beer` is a genuine addition, and it is drawn because it is in
    # `PH_SYSTEM_READS` (D-160's sampler-scope master).
    assert len(joint) == 3**11, "every drawn dimension must be varied, not a subset of them"
    degassed_joint = [f[0] for f in joint]
    vessel_joint = [f[1] for f in joint]
    # The measured frame first, because it is the one that compares with Tyrell (D-208).
    assert min(degassed_joint) == pytest.approx(0.712, abs=0.02), (
        f"the degassed joint low corner moved to {min(degassed_joint):.1%}; D-209 measures "
        "71.2 %, against D-208's 8.3 % and D-181's pre-carbonic 7.6 %. The low corner moved "
        "further than the high one for the same reason it did at D-182: a member predicting "
        "little acidification finishes at a higher pH, and the nitrogen term is a fixed charge "
        "removal, so it buys the most where the acids buy the least"
    )
    assert max(degassed_joint) == pytest.approx(1.401, abs=0.02), (
        f"the degassed joint high corner moved to {max(degassed_joint):.1%}; D-209 measures "
        "140.1 % (D-208's 82.7 %, D-181's pre-carbonic 82.2 %)"
    )
    # D-208 asserted `max(degassed_joint) < 1.0` here and wrote, in that assert's own failure
    # message, that if it ever fired "an acidification term genuinely arrived and the xfail
    # below should be re-scored before anything else is concluded". At D-209 it fired, the term
    # had arrived, and the xfail flipped green. So the claim is INVERTED rather than deleted: a
    # corner of the band must now reach, because a band whose top no longer covers the
    # measurement would mean the term had been lost again.
    assert max(degassed_joint) > 1.0, (
        f"no corner of the joint band reaches the measured drop in the degassed frame "
        f"({max(degassed_joint):.1%}). Since D-209 the nitrogen pool's charge is in the balance "
        "and the band is expected to straddle 100 %; if this fires, that term has gone missing "
        "or been re-scoped, and `test_the_model_reaches_tyrells_measured_beer_ph` should be the "
        "next thing looked at"
    )
    assert min(vessel_joint) == pytest.approx(0.928, abs=0.02), (
        f"the joint low corner moved to {min(vessel_joint):.1%}; D-209 measures 92.8 % over "
        "ELEVEN dimensions, D-183 measured 64.0 % over TEN "
        "dimensions (D-182 measured 63.8 % over nine; D-181's was 7.6 %) — "
        "yields at their low edge, peptide pKa HIGH, floors at their LOW edge (the "
        "strains that clear the most wort acid) and the seeds HIGH. A LOW floor means MORE "
        "acid removed and therefore a SMALLER net drop, which is the opposite of the "
        "intuition that a lower residue means a more acidic beer: what moves pH is charge "
        "lost, not acid left. "
        "THE LOW CORNER MOVED EIGHT TIMES MORE THAN THE HIGH ONE (7.6 -> 63.8 % against "
        "82.2 -> 109.4 %), and the asymmetry is the carbonic term's own geometry rather than "
        "a mistake: a member that predicts LITTLE acidification finishes at a HIGHER pH, and "
        "carbonic acid dissociates more the higher the pH, so the CO2 contribution is largest "
        "exactly where the acids contribute least. The joint band is therefore COMPRESSED, "
        "from a 74.6-point span to a 45.6-point one — this term is a stabiliser of the "
        "prediction, not just an offset to it."
    )
    assert max(vessel_joint) == pytest.approx(1.495, abs=0.02), (
        f"the joint high corner moved to {max(vessel_joint):.1%}; D-209 measures 149.5 % over "
        "ELEVEN dimensions, D-183 measured 109.7 % over TEN "
        "dimensions (D-182 measured 109.4 % over nine; D-181's "
        "was 82.2 %) — yields high, peptide pKa low, floors high, seeds low, and the three "
        "CO2 parameters at the edges that dissolve the most and dissociate it hardest. NB "
        "this is a CORNER of a 10-D hypercube, not a member any ensemble was seen to draw. "
        "**A corner REACHES the measured drop again, and it was pre-registered.** D-180 had "
        "one at 104.5 %; D-181 removed it and wrote that a future change restoring it would "
        "signal that an omitted term had arrived. One has — dissolved CO2, D-180's own arm C, "
        "predicted at 76-104 % before D-182 was written. Reaching is therefore the EXPECTED "
        "outcome here and not evidence the model is right — and since D-208 it is not even "
        "evidence about the measurement: this frame is the pH INSIDE the vessel, and the "
        "published number is a decarbonated reading, where the same corner reaches 82.7 %."
    )


def test_the_model_reaches_tyrells_measured_beer_ph():
    """THE ACCEPTANCE CLAIM, and since D-209 it **passes** — this was an ``xfail(strict)``.

    The history is the point, so it is kept rather than summarised. D-207 read Fig. 4's pH panel
    as a course and found the shape wrong while the endpoint metric passed. D-208 settled the
    measurement frame (a published beer pH is a decarbonated reading), whereupon the endpoint
    metric did not pass either, and this test was written as an expected failure on the day-7
    LEVEL: the model reached 43.2-62.9 % of the measured 0.81 pH drop and finished day 7 at
    5.2474 against a four-strain envelope of 4.804-4.916, for every retained-CO₂ fraction in
    [0, 1]. D-208 §5 chose the level over a day-1 shape pin precisely because *a fix flips it
    green, which is what strict=True is for*. D-209 is that fix, and it flipped: the assimilable
    nitrogen pool carries net positive charge, so the yeast taking it up acidifies the liquid,
    and that term was missing from the charge balance rather than from the chemistry.

    **What passing here does and does not license.** It says day 7 lands inside the envelope. It
    does NOT say the course is reproduced: day 1 is still outside, and it is outside by MORE than
    before D-209 (0.29-0.34 pH against 0.19), because the model empties its nitrogen inside 20 h
    and so delivers the whole charge step before the day-1 reading. That is a statement about
    nitrogen-uptake TIMING, not about the charge arithmetic, and it is pinned separately in
    :func:`test_the_day_1_pH_is_still_missed_and_the_miss_is_uptake_timing`. Nor is the term the
    whole of nitrogen uptake's effect: it is the charge half only, so it is a lower bound (see
    ``acidbase.nitrogen_charge_molar``), and two further same-sign sources — K⁺/H⁺ antiport and
    trub precipitation — are identified and unbuilt.

    ``TYRELL_PH_COURSE`` was shipped at D-207 as data no assert read. This is the assert that
    reads it, and it reads day 7 ONLY.
    """
    compiled, res = _run(dict(TYRELL_SCENARIO))
    params = compiled.parameters.resolve()
    states = np.asarray(res.y, dtype=float)
    t_h = np.asarray(res.t, dtype=float)
    # Day 7 lands between solver steps, so interpolate the state slot-wise rather than taking
    # the nearest sample: the pH is read off a state, not off a series of pHs.
    y = np.array([np.interp(7 * 24.0, t_h, states[i, :]) for i in range(states.shape[0])])
    modelled = acidbase.degassed_ph_of_state(y, compiled.schema, params)
    lo, hi = TYRELL_PH_COURSE[7]
    assert lo - TYRELL_PH_READ_TOL <= modelled <= hi + TYRELL_PH_READ_TOL, (
        f"modelled degassed day-7 pH {modelled:.4f} is outside Tyrell's four-strain envelope "
        f"{lo:.3f}-{hi:.3f} widened by the {TYRELL_PH_READ_TOL} pH read tolerance"
    )


def test_a_beer_seeded_at_finished_beer_levels_would_overshoot(beer_params):
    """Why the seeds HAD to move — kept as a test so the reason cannot be forgotten.

    D-179 dosed the ``*_typical_beer`` levels at pitch. Bolting a producer onto those seeds
    would add a whole ferment's production on top of an already-finished composition. This
    recomputes that counterfactual on the charge balance: it lands below any real beer, which
    is what makes "produce the acids" and "start from a wort" one decision rather than two.
    """
    params = beer_params.resolve()
    pka = acidbase.build_pka_map(params)
    # ALL FIVE finished-beer levels, citrate included. The counterfactual is the whole D-179
    # seed set, so leaving citrate out would both understate it and leave
    # ``citric_typical_beer`` read by nothing at all — which the YAML group-1 header claims is
    # not the case. A prose claim its guard does not back is the shape
    # [[feedback-grep-finds-claims-not-guards]] exists for.
    slots = ("lactic", "acetic", "citrate", "malic", "succinic")
    molar = {s: acidbase.ALL_ACIDS[s].molar_mass for s in slots}
    seeds = {
        s: beer_params[f"{'citric' if s == 'citrate' else s}_typical_beer"].value for s in slots
    }

    totals = {s: g / molar[s] for s, g in seeds.items()}
    totals["peptide_buffer"] = (
        beer_params["peptide_buffer_capacity_beer"].value
        / acidbase.ALL_ACIDS["peptide_buffer"].molar_mass
    )
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 4.4)
    produced = dict(totals)
    for spec in ORGANIC_ACID_SPECS:  # citrate is seeded but has no yield
        produced[spec.slot] = (
            produced[spec.slot] + (params[spec.yield_param] * TYRELL_SUGAR_GPL) / molar[spec.slot]
        )
    assert acidbase.solve_ph(produced, cation, 0.0, 0.0, pka) < 4.3, (
        "producing on top of finished-beer seeds should overshoot below any real beer; "
        "if it no longer does, the seeds or the yields have moved"
    )


def test_wine_and_beer_schemas_disagree_about_which_slots_are_produced():
    """The registry is beer's, and every slot it names must exist in beer and be absent or
    differently-sourced in wine — the structural version of "wine is unchanged"."""
    beer, wine = beer_schema(), wine_schema()
    for spec in ORGANIC_ACID_SPECS:
        assert spec.slot in beer
        assert carbon_mass_fraction(spec.species) > 0.0
    assert "acetic" not in wine and "succinic" not in wine


# ======================================================================================
# 5. The three acids that FALL — beer's missing base (decision D-181)
#
# D-180 closed with its agreement explicitly held open by two omitted terms of opposite sign
# and sized both. This section is the larger one built: three malt-derived wort acids a real
# ferment removes, which the model could not lose because they were not state. Every test here
# is written so that the beat's own direction — the headline agrees WORSE — is asserted rather
# than tolerated.
# ======================================================================================

#: Tyrell's EBC-tube day-0 and day-7 values for the three, in mg/L (Figs 6, 11, 7 — FIGURE
#: READS). The wort seed, then the four-strain day-7 range the floor is the mean of.
TYRELL_FALLING_PPM = {
    "pyruvic": (22.0, (0.2, 2.0)),
    "formic": (26.0, (3.5, 6.0)),
    "oxalic": (22.0, (5.0, 6.0)),
}


def test_the_falling_acids_reach_their_measured_floors_and_stop():
    """Acceptance: each acid falls from its wort seed and settles ON its measured floor.

    Named a WIRING check like its D-180 counterpart, not a validation: the floors were read
    off the same figures this compares against, so it is close to a round-trip. What it would
    catch is a wrong seed, a sign error, an mg/L slip, or a removal that overshoots the floor.

    The stronger claim is the second one — the pools do not merely *approach* the floor, they
    are AT it well before the run ends, which is what makes the finished pH insensitive to the
    rate constant (pinned separately below).
    """
    compiled, res = _run(dict(TYRELL_SCENARIO))
    params = compiled.parameters.resolve()
    for sink in WORT_ACID_SINKS:
        seed_ppm, (lo_ppm, hi_ppm) = TYRELL_FALLING_PPM[sink.slot]
        series = res.series(sink.slot)
        floor = params[sink.floor_param]
        assert float(series[0]) == pytest.approx(seed_ppm / 1000.0, rel=1e-9), (
            f"{sink.slot} was not seeded at Tyrell's measured wort level"
        )
        assert float(series[-1]) == pytest.approx(floor, rel=1e-6), (
            f"{sink.slot} finished at {float(series[-1]) * 1000:.2f} mg/L, not its measured "
            f"floor of {floor * 1000:.2f}"
        )
        # The floor is the four-strain mean, so it must sit inside the measured spread.
        assert lo_ppm / 1000.0 <= floor <= hi_ppm / 1000.0
        # Never below the floor — checked on the trajectory rather than on the rate law, so a
        # solver undershoot would show up here.
        assert float(series.min()) >= floor - 1e-12, f"{sink.slot} was driven below its floor"


def _finished_ph(over: dict[str, float], param_over: dict[str, float] | None = None):
    """Finished pH of a Tyrell beer, with optional initial- and parameter-level overrides.

    Goes through ``simulate_scheduled`` with ``compiled.events`` rather than a bare
    ``simulate`` — the D-35 trap: a hand-wired call silently drops the schedule. ``param_values``
    is a property returning a FRESH dict each access (a D-147 trap), so mutating the copy is
    safe and mutating it in place would have been a no-op.
    """
    from fermentation.runtime import simulate_scheduled

    compiled = compile_scenario(
        Scenario(
            name="d181",
            medium="beer",
            initial={**TYRELL_SCENARIO, **over},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=14.0,
        )
    )
    params = compiled.param_values
    params.update(param_over or {})
    res = simulate_scheduled(
        compiled.process_set, params, compiled.y0, compiled.t_span_h, events=compiled.events
    )
    # ``ph_of_state`` on the final column rather than the whole series: only the endpoint is
    # read here, and it keeps the helper working on a ScheduledTrajectory (which duck-types
    # Trajectory but is not one).
    ph = acidbase.ph_of_state(res.y[:, -1], compiled.schema, params)
    return float(ph), compiled, res


def test_the_removal_rate_does_not_move_the_finished_ph(beer_params):
    """P4, pre-registered: the rate sets the SHAPE of the first days, not the outcome.

    Every value in ``k_wort_acid_removal``'s band clears the fall inside the first days of a
    14-day run, so the finished pH must be insensitive to it. This matters because the rate is
    the weakest number in the beat — inverted from a single day-1 figure read, with two of the
    twelve strain-acid fits unrepresentable — and it is what licenses shipping ONE shared
    constant instead of three: if the endpoint moved with it, the parsimony would be buying
    error rather than honesty.

    A one-sided insensitivity claim is worthless without a positive control, so the same
    comparison is run on a quantity that MUST move: the early pool, which is what the rate
    actually sets [[feedback-a-null-result-needs-a-positive-control]].
    """
    base = beer_params["k_wort_acid_removal"].value

    def arm(scale: float):
        ph, _compiled, res = _finished_ph({}, {"k_wort_acid_removal": base * scale})
        # The pool ~24 h in — where the rate constant actually lives.
        day1 = int(res.t.searchsorted(24.0))
        return ph, float(res.series("formic")[day1])

    slow_ph, slow_day1 = arm(1.0 / 3.0)
    fast_ph, fast_day1 = arm(3.0)

    assert abs(fast_ph - slow_ph) < 0.005, (
        f"tripling and thirding the removal rate moved the finished pH by "
        f"{abs(fast_ph - slow_ph):.4f}. The floor is no longer being reached inside the run, "
        "so the endpoint has started depending on a number that only sets the transient."
    )
    assert slow_day1 > fast_day1 * 1.5, (
        f"the positive control did not move: day-1 formic was {slow_day1 * 1000:.3f} mg/L at "
        f"1/3x and {fast_day1 * 1000:.3f} at 3x. If the rate moves neither the transient nor "
        "the endpoint, the Process is not running and the insensitivity above is vacuous."
    )


def test_the_falling_acids_are_off_every_ledger_by_construction():
    """P5: closure is untouched, and NOT because the removal is balanced.

    The three slots carry real malt carbon out of the beer with no destination pool, which is
    only legitimate because they are weighted nowhere — the ``iso_alpha`` treatment for
    exogenous mass removed by an unattributed route. Two independent checks, because they fail
    apart: the ledger weight really is zero at those slots, and the species really are absent
    from the chemistry registries, so a future producer drawing one of them out of ``S`` raises
    instead of silently leaking carbon.
    """
    from fermentation.core.chemistry import CARBON_ATOMS, MOLAR_MASS
    from fermentation.validation.conservation import total_carbon

    compiled, res = _run(dict(TYRELL_SCENARIO))
    weights = total_carbon(
        compiled.schema,
        biomass_carbon_fraction=compiled.parameters["biomass_C_fraction"].value,
    )
    probe = compiled.schema.zeros()
    for sink in WORT_ACID_SINKS:
        probe[compiled.schema.slice(sink.slot)] = 1.0
        assert weights(probe) == 0.0, (
            f"{sink.slot} carries a carbon weight, so WortAcidRemoval destroys carbon: either "
            "weight it AND give the removal a destination pool, or leave it off the ledger"
        )
        probe[compiled.schema.slice(sink.slot)] = 0.0
        for name in (sink.slot, f"{sink.slot}_acid"):
            assert name not in MOLAR_MASS and name not in CARBON_ATOMS

    start, end = weights(res.y[:, 0]), weights(res.y[:, -1])
    assert end == pytest.approx(start, rel=1e-9), (
        f"carbon drifted {end - start:.6e} g C/L with the wort-acid sink active"
    )


def test_removing_the_falling_acids_raises_the_finished_ph_by_the_predicted_amount():
    """THE BEAT'S OWN DIRECTION, asserted: this makes the model agree WORSE, on purpose.

    Run against the pre-D-181 configuration — the same engine with the three acids dosed to
    zero, which is exactly the one-sided acidification D-180 shipped — the finished pH must be
    HIGHER and the predicted drop SMALLER. D-180 sized the omission at roughly +0.2-0.3 pH from
    a static charge balance; this checks the built version lands in that window rather than
    merely moving the right way.

    Both arms are the same scenario with one override, so the ONLY difference is the three
    acids [[feedback-pair-the-red-with-an-ordering-preserving-baseline]] — in particular the
    peptide capacity is the same re-anchored value in both, so this measures the acids and not
    the re-anchor that shipped beside them.

    **THIS NUMBER HAS NOW MOVED THREE TIMES, AND ONLY THE FIRST MOVE WAS ABOUT THE ACIDS.**
    D-181 built them and measured +0.2094 pH with no carbonic term in the balance. D-182 added
    dissolved CO2 and the shipped figure became +0.1128 — halved by buffering, not by the acids
    shrinking, which is what the CO2-free arm below existed to prove. D-209 then put the
    assimilable nitrogen pool's charge into the balance, which drops the finished pH by ~0.35,
    and BOTH figures moved again: the CO2-free one fell to **+0.1417** while the shipped one
    barely budged, to **+0.1111**.

    The direction is the acids' own speciation. Pyruvic (pKa 2.39), formic (3.75) and oxalic
    (1.25/4.14) are LESS dissociated at a lower pH, so removing them removes less charge — the
    CO2-free measurement of the removal shrinks. Carbonic meanwhile dissociates less at the
    lower pH too, so it buffers less against the removal, which pushes the shipped figure the
    other way. The two therefore CONVERGE, from 0.097 apart to 0.031, and that convergence is
    asserted below. The lesson generalises past this test: **none of these separately-measured
    term sizes is a constant of the model — each is quoted at the pH the model reached when it
    was measured**, so adding them to predict a total over-counts, and re-quoting an old one
    without re-running it is a stale number [[feedback-a-summary-statistic-is-not-the-curve]].

    D-182's own paragraph, kept because its reasoning is still correct at its own pH: D-181 measured
    this at +0.2094 pH with NO carbonic term in the balance. Dissolved CO2 (D-182) does not
    merely add its own shift on top: it **buffers against this one**. Carbonic's dissociated
    fraction RISES with pH (pKa 6.43 sits above the beer, so the model is on the steep side of
    its curve), so removing the falling acids — which pushes pH up — makes carbonic give up
    more charge, which pushes back. The same three acids are therefore worth +0.1128 pH in the
    shipped model. **The two omitted terms D-180 named are not additive**, and anyone adding
    their separately-measured sizes to predict a total will over-count.

    So this asserts BOTH numbers: the shipped one, and D-181's own recomputed CO2-free on the
    identical pair of end states. The second is what keeps D-181's claim falsifiable instead
    of quietly superseded — and it is only possible because ``solve_ph`` takes the carbonic
    term as an explicit argument rather than reading it implicitly.
    """
    with_sinks, c1, r1 = _finished_ph({})
    without, c2, r2 = _finished_ph({"pyruvic_gpl": 0.0, "formic_gpl": 0.0, "oxalic_gpl": 0.0})

    assert with_sinks > without, (
        "removing the falling acids must RAISE the finished pH — that is the missing base "
        f"D-180 sized. Got {with_sinks:.4f} with them and {without:.4f} without."
    )
    assert with_sinks - without == pytest.approx(0.1128, abs=0.01), (
        f"the missing base is worth {with_sinks - without:.4f} pH in the shipped model; "
        "D-182 measured 0.1128. It is SMALLER than D-181's CO2-free 0.2094 because the "
        "carbonic term buffers against the removal, not because the acids shrank — the "
        "CO2-free assertion below is what tells those two explanations apart."
    )

    # The same two end states, re-solved with the carbonic term switched off. Until D-209 this
    # reproduced D-181's +0.2094 pH, which was what told the CO2-buffer explanation apart from
    # "something in the seeds, floors or charge arithmetic actually moved". It no longer does,
    # and the reason is the third movement of this same number (see the docstring).
    def co2_free(compiled, res) -> float:
        params = compiled.parameters.resolve()
        y = res.y[:, -1]
        return acidbase.solve_ph(
            acidbase._totals_molar(y, compiled.schema),
            # ``params`` became required at D-209 (the nitrogen pool joined the cation side);
            # this is exactly ``degassed_ph_of_state``, kept spelled out because the point of
            # the arm is that only the carbonic argument differs from the shipped solve.
            acidbase._cation(y, compiled.schema, params),
            acidbase._byp_succinic_molar(y, compiled.schema),
            0.0,
            acidbase.build_pka_map(params),
        )

    co2_free_gap = co2_free(c1, r1) - co2_free(c2, r2)
    assert co2_free_gap == pytest.approx(0.1417, abs=0.01), (
        f"with the carbonic term off the missing base is worth {co2_free_gap:.4f} pH; D-209 "
        "measures 0.1417, against D-181's 0.2094 on the SAME arm. The fall is the acids' own "
        "speciation at a lower pH, not a change to the acids: see the docstring. If this moves "
        "back toward 0.21 the nitrogen term has stopped acting on the finished pH."
    )
    # The two numbers have CONVERGED, and pinning the gap between them is what makes the
    # explanation falsifiable rather than merely plausible: if the fall were caused by something
    # other than the finished pH moving down, there is no reason the carbonic-buffered and
    # CO2-free measurements of the same removal would end up 0.03 apart instead of 0.10.
    assert co2_free_gap - (with_sinks - without) == pytest.approx(0.031, abs=0.01), (
        f"the CO2-free and shipped measurements of the same removal differ by "
        f"{co2_free_gap - (with_sinks - without):.4f} pH; D-209 measures 0.031, against D-182's "
        "0.097. Carbonic buffers against the removal in proportion to how dissociated it is, "
        "and at D-209's lower finished pH it is less dissociated, so it pushes back less."
    )


def test_the_peptide_capacity_still_reproduces_peyers_published_wort_bc():
    """The re-anchor is a round-trip, so the round-trip is a test.

    ``peptide_buffer_capacity_beer`` is back-solved so that this engine's charge balance
    reproduces Peyer's published wort BC = 1.18 on top of the wort acid table. D-181 puts three
    more acids in that table, so the capacity had to move (1.6708 -> 1.5481 g/L, -7.34 %) — and
    a back-solved constant with no test is a number that silently stops meaning anything the
    next time the table changes. This recomputes the titration from the SHIPPED parameters.

    Two of the three stacked mismatches the capacity's own block admits to are untouched by
    this check and stay open: it is fitted at wort pH and applied across a traverse, and 1.18
    is a wort measurement with no published finished-beer counterpart.
    """
    import math

    v_in, v_acid, c_acid = 25.0, 0.375, 1.0
    v_fin = v_in + v_acid
    data = default_data_dir()
    params = load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    )
    pka = acidbase.build_pka_map(params.resolve())
    seeds = {
        "lactic": "lactic_typical_wort",
        "acetic": "acetic_typical_wort",
        "citrate": "citric_typical_wort",
        "malic": "malic_typical_wort",
        "succinic": "succinic_typical_wort",
        "pyruvic": "pyruvic_typical_wort",
        "formic": "formic_typical_wort",
        "oxalic": "oxalic_typical_wort",
    }
    totals = {
        slot: params[p].value / acidbase.ALL_ACIDS[slot].molar_mass for slot, p in seeds.items()
    }
    totals["peptide_buffer"] = (
        params["peptide_buffer_capacity_beer"].value
        / acidbase.ALL_ACIDS["peptide_buffer"].molar_mass
    )
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 5.5)
    ph_in = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    f = v_in / v_fin
    diluted = {k: v * f for k, v in totals.items()}
    ph_fin = acidbase.solve_ph(diluted, cation * f - (v_acid * c_acid) / v_fin, 0.0, 0.0, pka)
    h_inc = v_fin * 10.0 ** (-ph_fin) - v_in * 10.0 ** (-ph_in)
    bc = math.log10((v_acid * c_acid) / h_inc)

    assert bc == pytest.approx(1.18, abs=1e-9), (
        f"the shipped peptide capacity reproduces BC = {bc:.6f}, not Peyer's published 1.18. "
        "The wort acid table changed without the back-solve being re-run."
    )


def test_no_process_touches_the_peptide_buffer_pool():
    """The peptide pool is ANCHOR-TIME state, and a Process that drained it would break beer.

    Trub settling was carried as an unbuilt same-sign acidification from D-209 sec 8 until D-214
    measured it. It is not buildable as a fermentation-phase Process, and the reason is structural
    rather than a missing rate: ``peptide_buffer`` rides ``_BEER_ACID_SEEDS``, so the t=0 cation
    back-solve is fitted WITH this pool present as the counter-anion. Removing it after the anchor
    leaves that cation with nothing to balance against, and the charge balance answers with
    hydroxide. Measured at D-214: cutting 20 % of the pool at 6 h takes day 1 from pH 5.448 to
    **7.085**, and cutting all of it takes day 7 to **11.66**. That is not a small same-sign term,
    it is a charge-balance violation.

    So this is a guard on a hazard the archive nearly walked into, not on a shipped number. The
    mechanism is real — protein does drop out — but it drops out as hot trub in the BOIL and as
    cold break during CHILLING, both before pitching, and *Chemistry of Beer* sec 2.9 says the
    coagulate is "removed before the wort is fermented". The model already carries it there: the
    capacity is back-solved from Peyer's 1.18 **control wort**, an already-boiled figure.

    The pool being real state is exactly what makes this reachable — expressible is not buildable,
    which is D-205's lesson arriving from the other direction.
    """
    beer = get_medium("beer").build_process_set(strict=True)
    touching = sorted(p.name for p in beer.active if "peptide_buffer" in p.touches)
    assert not touching, (
        f"{touching} declares `peptide_buffer` in `touches`. That pool is the t=0 anchor's "
        "counter-anion, not a fermentation-phase reservoir: draining it post-anchor drives the "
        "charge balance alkaline (D-214 measures pH 7.08 for a 20 % cut at 6 h, 11.66 for a full "
        "one). If a beat really means to model less buffering protein, it belongs PRE-anchor, in "
        "`peptide_buffer_capacity_beer` or the `peptide_buffer_gpl` scenario key, where the "
        "cation back-solve sees it."
    )


def _tyrell_ph_with_peptide_loss(data_dir, loss: float, day: float) -> float:
    """Beer's degassed pH at ``day`` with ``loss`` of the wort's buffering protein gone PRE-pitch.

    The loss is applied through the ``peptide_buffer_gpl`` scenario key, which is seeded before
    the cation back-solve — the placement the sources actually describe (boil + chill), and the
    only placement that is not a charge-balance violation (see the test above).
    """
    params = load_parameters(data_dir / "beer_acids.yaml")
    capacity = params["peptide_buffer_capacity_beer"].value
    compiled = compile_scenario(
        Scenario(
            name="d214-peptide-loss",
            medium="beer",
            initial={**TYRELL_SCENARIO, "peptide_buffer_gpl": capacity * (1.0 - loss)},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=14.0,
        ),
        data_dir=data_dir,
    )
    res = compiled.run()
    resolved = compiled.parameters.resolve()
    states = np.asarray(res.y, dtype=float)
    t_h = np.asarray(res.t, dtype=float)
    y = np.array([np.interp(day * 24.0, t_h, states[i, :]) for i in range(states.shape[0])])
    return float(acidbase.degassed_ph_of_state(y, compiled.schema, resolved))


def test_losing_wort_protein_acidifies_late_not_early():
    """The shape claim that refuses trub settling — asserted, because it is the whole argument.

    D-211 sec 9 left a brief for the next beer-pH beat: beer wants **acidification early and none
    late**, because at the high ``nitrogen_uptake_charge_beer`` edge day 1 sits above its ceiling
    while day 7 has almost no headroom above its floor. Less buffering protein acidifies — that
    part of D-209 sec 8 is right — but it does so with the OPPOSITE profile, because removing
    buffer amplifies acid production, and acid production is cumulative. So the effect grows with
    time instead of fading, and it is disqualified by its shape rather than by its size.

    Pinning the ratio and not just the two deltas is the point: a beat that later proposes any
    buffer-removal term is answered by this test, whatever magnitude it picks.
    """
    data_dir = default_data_dir()
    base_1 = _tyrell_ph_with_peptide_loss(data_dir, 0.0, 1.0)
    base_7 = _tyrell_ph_with_peptide_loss(data_dir, 0.0, 7.0)
    lost_1 = _tyrell_ph_with_peptide_loss(data_dir, 0.20, 1.0)
    lost_7 = _tyrell_ph_with_peptide_loss(data_dir, 0.20, 7.0)

    # Tolerances are a solver-noise allowance, not a band: both runs are deterministic, so the
    # measured 0.017252 / 0.058699 should reproduce exactly unless the model moved. A looser
    # tolerance here would let the two values drift far enough to break the ratio below while
    # both still "passed" [[feedback-pin-tolerance-vs-solver-tolerance]].
    early, late = base_1 - lost_1, base_7 - lost_7
    assert early == pytest.approx(0.01725, abs=0.0005), (
        f"a 20 % pre-pitch protein loss moves day 1 by {early:.6f} pH; D-214 measures 0.017252"
    )
    assert late == pytest.approx(0.05870, abs=0.0015), (
        f"a 20 % pre-pitch protein loss moves day 7 by {late:.6f} pH; D-214 measures 0.058699"
    )
    assert late > 3.0 * early, (
        f"buffer removal is supposed to be LATE-weighted (D-214 measures day 7 at 3.4x day 1); "
        f"here day 7 is {late / early:.2f}x day 1. If this ratio has fallen below 3, the shape "
        "argument that refuses trub settling as an answer to D-211 sec 9's brief no longer holds "
        "and the refusal needs re-measuring, not re-asserting."
    )


def test_the_trub_window_is_empty_at_the_edge_that_parks_it(tmp_path, beer_params):
    """The refusal as arithmetic, on the arm the parking question actually lives on.

    D-210 parked trub settling on the HIGH ``nitrogen_uptake_charge_beer`` edge and D-211 sec 9
    re-priced it there, so the nominal cannot answer it [[feedback-pin-the-band-not-the-nominal]].
    At that edge D-214 measures: day 1 needs a loss of **>= 27.6 %** to come inside its ceiling,
    while day 7 can afford **<= 3.1 %** before falling through its floor. The window is empty by
    about ninefold, so no protein-loss fraction satisfies both ends — which is why the refusal is
    arithmetic rather than a judgement about how much trub a brewery carries over.

    This asserts the tight half: a 5 % loss, well below what day 1 would need, already takes day 7
    out of the envelope.

    **Both margins here are small, and that is stated rather than hidden.** The control passes by
    **+0.0086** pH and the assertion fails the floor by **−0.0054**. The archive has reported that
    control quantity two ways — D-211 §9 said 0.0082, D-214 §7 says 0.0086 — and the difference was
    a read artefact, not a model change (`argmin` over the solver's adaptive output against
    `np.interp` onto the exact hour; this helper uses the latter, which is why it should read
    0.0086). So if this test ever fails on the CONTROL line, check the read before the model: the
    pass margin is comparable to the size of that old disagreement
    [[feedback-read-a-fast-curve-on-a-fixed-grid]].
    """
    hi = beer_params["nitrogen_uptake_charge_beer"].uncertainty.high
    data_dir = _beer_data_dir_with_nitrogen_charge(tmp_path, hi)
    floor_7 = 4.804 - 0.024

    assert _tyrell_ph_with_peptide_loss(data_dir, 0.0, 7.0) > floor_7, (
        "the unmodified high edge should still be inside the day-7 envelope (D-211 sec 9 "
        "measures 0.0086 pH of headroom on a fixed-grid read); if it is not, this test's "
        "premise moved and the trub arithmetic below is scored against the wrong baseline"
    )
    with_loss_7 = _tyrell_ph_with_peptide_loss(data_dir, 0.05, 7.0)
    assert with_loss_7 < floor_7, (
        f"a 5 % pre-pitch protein loss puts day 7 at {with_loss_7:.4f}, which should be BELOW the "
        f"admissible floor {floor_7:.4f}. D-214's refusal rests on day 7 affording only ~3 % while "
        "day 1 needs ~28 %; if 5 % now fits, that window has opened and the refusal is stale."
    )


def test_wine_never_wires_or_reads_the_wort_acid_sink():
    """Beer-only, structurally — the D-180 isolability claim, restated for the sink half."""
    wine = get_medium("wine").build_process_set(strict=True)
    assert WortAcidRemoval.name not in {p.name for p in wine.active}
    declared = {name for p in wine.active for name in p.reads}
    assert not (set(WortAcidRemoval.reads) & declared)
    schema = wine_schema()
    for sink in WORT_ACID_SINKS:
        assert sink.slot not in schema
    # Wine carries the same MOLECULE as `pyruvic` under a different slot name, and that pool
    # must stay charge-inactive: a union registry would make it charge-active behind this
    # beat's back, which is the D-179 argument for per-medium registries, one acid later.
    assert "pyruvate" in schema
    assert "pyruvate" not in acidbase.WINE_ACIDS


def test_an_unanchored_beer_keeps_the_sink_disabled_and_its_slots_empty():
    """The same opt-in gate as the producer, for the same correctness reason (D-179/D-180)."""
    unanchored = {k: v for k, v in TYRELL_SCENARIO.items() if k != "initial_ph"}
    compiled, res = _run(unanchored)
    assert WortAcidRemoval.name not in {p.name for p in compiled.process_set.active}
    for sink in WORT_ACID_SINKS:
        assert float(res.series(sink.slot)[0]) == 0.0
        assert float(res.series(sink.slot)[-1]) == 0.0


def test_the_sink_honours_its_touches_contract():
    """``strict=True`` turns ``touches`` into a contract — and here it carries the claim.

    ``WortAcidRemoval`` asserts no mechanism, and the machine-checkable form of that is its
    ``touches``: three acid slots and nothing else. If ``S``, ``E`` or ``CO2`` ever appears
    here, the Process has started claiming where the carbon went.
    """
    pset = get_medium("beer").build_process_set(strict=True)
    assert WortAcidRemoval.name in pset
    schema = beer_schema()
    data = default_data_dir()
    params = load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    ).resolve()
    state = _beer_state(schema, pyruvic=0.022, formic=0.026, oxalic=0.022)
    d = WortAcidRemoval().derivatives(0.0, state, schema, params)
    touched = {n for n in schema.names if float(abs(d[schema.slice(n)]).sum()) > 0.0}
    assert touched == {sink.slot for sink in WORT_ACID_SINKS}
    for sink in WORT_ACID_SINKS:
        assert float(d[schema.slice(sink.slot)][0]) < 0.0, f"{sink.slot} is not being removed"


def test_an_acid_at_or_below_its_floor_is_frozen_not_refilled():
    """The clamp is at the FLOOR, not at zero, and it must not run backwards.

    A pool already at its floor has no removal term; a pool below it (solver undershoot, or a
    scenario that doses less than the floor) must be frozen rather than topped back up. The
    second case is the one a naive ``-k * (pool - floor)`` gets wrong: without the guard the
    sign flips and the Process becomes a producer of an acid nothing produces.
    """
    schema = beer_schema()
    data = default_data_dir()
    params = load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    ).resolve()
    for sink in WORT_ACID_SINKS:
        floor = params[sink.floor_param]
        for level in (floor, floor * 0.5, 0.0, -1e-12):
            state = _beer_state(schema, **{sink.slot: level})
            d = WortAcidRemoval().derivatives(0.0, state, schema, params)
            assert float(d[schema.slice(sink.slot)][0]) == 0.0, (
                f"{sink.slot} at {level} g/L (floor {floor}) produced a nonzero rate"
            )


def _beer_data_dir_with_nitrogen_charge(tmp_path, value: float):
    """A copy of the packaged parameter dir with ``nitrogen_uptake_charge_beer`` set to ``value``.

    Heavier than a ``param_values`` patch, and deliberately so: this parameter is read at COMPILE
    time (both anchors subtract the nitrogen charge off the fitted cation) *and* at runtime (every
    rate that reads pH). Patching the resolved map after compile would move only the runtime half
    and leave the cation slot carrying the shipped value's subtraction — a silently wrong arm
    rather than a broken one [[feedback-a-parameter-can-be-pinned-and-drawn]]. Copying the dir is
    what makes both halves move together, so this helper is also the guard on that coupling.
    """
    dest = tmp_path / f"data_{value}"
    shutil.copytree(default_data_dir(), dest)
    path = dest / "acidbase.yaml"
    text = path.read_text(encoding="utf-8")
    text, n_value = re.subn(
        r"^(nitrogen_uptake_charge_beer:\n  value: )[-0-9.eE]+$",
        r"\g<1>" + repr(value),
        text,
        flags=re.MULTILINE,
    )
    zero_width = "{ low: " + repr(value) + ", high: " + repr(value) + ', note: "band arm" }'
    band_pattern = (
        r"^(nitrogen_uptake_charge_beer:\n(?:  (?!uncertainty)[^\n]*\n)*"
        r"  uncertainty: )\{[^\n]*\}$"
    )
    text, n_band = re.subn(band_pattern, r"\g<1>" + zero_width, text, flags=re.MULTILINE)
    assert (n_value, n_band) == (1, 1), "the parameter file's shape moved; fix this helper"
    path.write_text(text, encoding="utf-8")
    # The Parameter schema rejects a value outside its own band, so a loadable arm has to move
    # both — checked here rather than discovered as a confusing compile error downstream.
    assert load_parameters(path)["nitrogen_uptake_charge_beer"].value == value
    return dest


def _tyrell_degassed_ph_at_day(data_dir, day: float) -> float:
    compiled = compile_scenario(
        Scenario(
            name="d209-edge",
            medium="beer",
            initial=dict(TYRELL_SCENARIO),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=14.0,
        ),
        data_dir=data_dir,
    )
    res = compiled.run()
    params = compiled.parameters.resolve()
    states = np.asarray(res.y, dtype=float)
    t_h = np.asarray(res.t, dtype=float)
    y = np.array([np.interp(day * 24.0, t_h, states[i, :]) for i in range(states.shape[0])])
    return float(acidbase.degassed_ph_of_state(y, compiled.schema, params))


def test_both_edges_of_the_nitrogen_charge_band_keep_day_7_in_the_envelope(tmp_path, beer_params):
    """BOTH band edges, one threshold each — not a shared floor the tight edge rides for free.

    ``test_the_model_reaches_tyrells_measured_beer_ph`` scores the NOMINAL value only, and a
    nominal that passes says nothing about a band whose edges the sampler actually draws
    [[feedback-pin-the-band-not-the-nominal]]. This runs the full model at the low, nominal and
    high edges of ``nitrogen_uptake_charge_beer`` and pins each day-7 pH separately.

    **The finding this test exists to hold visible: the HIGH edge very nearly overshoots.** With
    the envelope at 4.804-4.916 and a 0.024 pH read tolerance, the admissible window is
    [4.780, 4.940]; the high edge lands at 4.783, a margin of 0.003 pH. So the derived band does
    not merely reach the measurement, it straddles the point of going past it — which is the
    honest statement about a term whose size was derived from published wort composition and
    never fitted, and which is a LOWER bound besides (the buffer-removal half of nitrogen uptake
    is inexpressible here, and it pushes the same way). A future beat that adds any further
    acidification has to re-examine this edge rather than celebrate it.
    """
    param = beer_params["nitrogen_uptake_charge_beer"]
    lo_ph = _tyrell_degassed_ph_at_day(
        _beer_data_dir_with_nitrogen_charge(tmp_path, param.uncertainty.low), 7.0
    )
    nom_ph = _tyrell_degassed_ph_at_day(
        _beer_data_dir_with_nitrogen_charge(tmp_path, param.value), 7.0
    )
    hi_ph = _tyrell_degassed_ph_at_day(
        _beer_data_dir_with_nitrogen_charge(tmp_path, param.uncertainty.high), 7.0
    )
    assert lo_ph > nom_ph > hi_ph, "more charge per mole N must mean a more acidic beer"

    lo_bound, hi_bound = TYRELL_PH_COURSE[7]
    window = (lo_bound - TYRELL_PH_READ_TOL, hi_bound + TYRELL_PH_READ_TOL)
    for label, value in (("low", lo_ph), ("nominal", nom_ph), ("high", hi_ph)):
        assert window[0] <= value <= window[1], (
            f"the {label} edge finishes day 7 at {value:.4f}, outside Tyrell's envelope "
            f"{lo_bound:.3f}-{hi_bound:.3f} widened by the {TYRELL_PH_READ_TOL} read tolerance"
        )
    # Each edge pinned on its own, so a shift that moved them together could not hide inside a
    # single containment check.
    assert lo_ph == pytest.approx(4.826, abs=0.01)
    assert nom_ph == pytest.approx(4.804, abs=0.01)
    assert hi_ph == pytest.approx(4.783, abs=0.01)
    # ...and the high edge's margin, asserted as the small number it is rather than described.
    assert 0.0 < hi_ph - window[0] < 0.01, (
        f"the high edge's margin to the bottom of the admissible window is {hi_ph - window[0]:.4f}"
        " pH; D-209 measures 0.003. If this has grown, something reduced the model's "
        "acidification and the day-7 agreement is no longer as tight as D-209 recorded"
    )


def test_the_day_1_pH_miss_survives_the_timing_fix_and_has_CHANGED_SIDES():
    """THE REMAINING DEFECT, re-pinned after D-211 removed its stated cause.

    D-209 left day 1 **0.315 pH too acidic** and attributed it to uptake TIMING: the whole
    nitrogen charge step landed inside ~20 h, before the day-1 reading. D-211 measured that
    timing against Tyrell's own cell-count panel and re-derived ``mu_max`` from it, and the
    charge step now takes ~2.5 days like the measured crop does.

    What is left is **0.046 pH on the ALKALINE side** — the miss shrank ~6x and reversed. That
    reversal is the point of this test. A residual that is now too alkaline is exactly the
    direction of the half D-209 section 8 recorded as unbuilt: assimilation removes the
    nitrogen pool's CHARGE but not its BUFFERING, and buffer removal pushes pH down. So the
    remaining day-1 gap is evidence for a term already identified, not for the timing being
    wrong again.

    **This test must not be "fixed" by moving ``mu_max``.** The value is calibrated on the cell
    counts (``tests/test_kinetics_growth.py``), which is the observable a growth rate IS; the pH
    course is three mechanisms downstream with a known-missing term. Tuning the rate to close
    this gap would consume the headroom that term needs and mis-attribute a missing buffer to a
    growth rate. The 60-hour ramp counterfactual D-209 used to locate the defect is retired with
    the defect: it priced a gap that no longer exists.
    """
    compiled, res = _run(dict(TYRELL_SCENARIO))
    params = compiled.parameters.resolve()
    schema = compiled.schema
    states = np.asarray(res.y, dtype=float)
    t_h = np.asarray(res.t, dtype=float)

    def state_at(hours: float):
        return np.array([np.interp(hours, t_h, states[i, :]) for i in range(states.shape[0])])

    lo, hi = TYRELL_PH_COURSE[1]
    day1 = float(acidbase.degassed_ph_of_state(state_at(24.0), schema, params))
    assert day1 > hi + TYRELL_PH_READ_TOL, (
        f"day 1 reads {day1:.4f}, inside or below Tyrell's {lo:.3f}-{hi:.3f}. D-211 leaves it "
        "ABOVE — too alkaline. If it is now inside, the last identified beer-pH gap has closed "
        "and the unbuilt buffer-removal half (D-209 section 8) needs re-pricing"
    )
    assert day1 - hi == pytest.approx(0.070, abs=0.02), (
        f"the day-1 miss is {day1 - hi:.4f} pH above the envelope; D-211 measures 0.070, where "
        f"D-209 measured 0.315 BELOW and D-208 0.186 above. All three matter: the size and the "
        "SIDE are what say whether a beat moved the timing or the arithmetic"
    )

    # Every other measured day is inside. The claim is that ONE day is missed, not that the
    # course is roughly right — a count is what stops a later beat trading days for days.
    inside = [
        band_lo - TYRELL_PH_READ_TOL
        <= float(acidbase.degassed_ph_of_state(state_at(day * 24.0), schema, params))
        <= band_hi + TYRELL_PH_READ_TOL
        for day, (band_lo, band_hi) in TYRELL_PH_COURSE.items()
    ]
    missed = [d for d, ok in zip(TYRELL_PH_COURSE, inside, strict=True) if not ok]
    assert sum(inside) == 7 and not inside[1], (
        f"{sum(inside)}/8 days inside with misses at {missed}; D-211 leaves day 1 alone"
    )

    # The attribution, measured rather than asserted: the charge step is NO LONGER complete by
    # the day-1 reading. This is the exact quantity D-209 pinned at >99 %, inverted.
    nitrogen = states[schema.slice("N").start, :]
    consumed_by_24h = float(np.interp(24.0, t_h, nitrogen[0] - nitrogen))
    fraction = consumed_by_24h / (nitrogen[0] - min(nitrogen))
    assert fraction == pytest.approx(0.363, abs=0.03), (
        f"{fraction:.3f} of the nitrogen is drawn by 24 h; D-209 measured >0.99 and D-211's "
        "re-derived rate puts it at 0.363, inside Tyrell's measured 0.234-0.448 cell-count "
        "spread. That agreement is what makes the timing MEASURED rather than fitted"
    )


# ---------------------------------------------------------------------------------------
# D-216 — the two anchors on beer's fermentation SPEED, and what they forbid
# ---------------------------------------------------------------------------------------

#: The ``q_sugar_max`` that reproduces Tyrell's day-2 extract fraction exactly (D-216 §3),
#: found by bisection on the shipped model. It is INSIDE the parameter's own printed band
#: (0.3-1.5), which is why this needs a test rather than a sentence: a future beat looking at
#: :func:`test_the_model_ferments_tyrells_wort_on_tyrells_schedule` will find that the one knob
#: that closes it is not even out of band. What stops it is the OTHER anchor, below.
Q_SUGAR_MAX_MATCHING_TYRELL = 1.397

#: ``K_repression`` removed entirely — not a candidate value, the unbounded LIMIT of the term
#: that owns 79 % of the lag (D-216 §5). Used to make the refusal two-tier.
K_REPRESSION_REMOVED = 1.0e6


def _beer_days_to_target_gravity(q_sugar_max: float | None = None, **overrides: float) -> float:
    """Days for §2.2's 1.048 ale wort at 20 °C to reach 1.010 apparent, optionally re-rated.

    Deliberately reuses the benchmark's own wort and gravity construction rather than
    restating them, so this test cannot drift away from the anchor it is about.
    """
    from tests.benchmarks.test_milestone1 import (
        _BEER_FERMENTABLE_S0,
        _BEER_OG_SG,
        TARGET_FG_SG,
        _apparent_gravity_series,
        _beer_scenario,
    )

    compiled = compile_scenario(_beer_scenario(duration_days=14.0))
    params = dict(compiled.param_values)
    if q_sugar_max is not None:
        params["q_sugar_max"] = q_sugar_max
    for name, value in overrides.items():
        assert name in params, f"{name} is not a compiled parameter"
        params[name] = value
    grid = np.linspace(0.0, compiled.t_span_h[1], int(compiled.t_span_h[1]) + 1)
    traj = simulate(compiled.process_set, params, compiled.y0, compiled.t_span_h, t_eval=grid)
    assert traj.success, traj.message
    apparent = _apparent_gravity_series(traj, _BEER_OG_SG, _BEER_FERMENTABLE_S0)
    reached = np.where(apparent <= TARGET_FG_SG)[0]
    return float(traj.t[reached[0]] / 24.0) if reached.size else float("inf")


def _tyrell_flux_fraction(
    days: int = 7, *, celsius: float = TYRELL_TRIAL_CELSIUS, **overrides: float
) -> dict[int, float]:
    """Fraction of Tyrell's fermentable consumed per day, on a FIXED grid (D-214's lesson)."""
    compiled = compile_scenario(
        Scenario(
            name="d216",
            medium="beer",
            initial=dict(TYRELL_SCENARIO),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
            duration_days=float(days),
        )
    )
    params = dict(compiled.param_values)
    for name, value in overrides.items():
        assert name in params, f"{name} is not a compiled parameter"
        params[name] = value
    grid = np.linspace(0.0, days * 24.0, days * 24 + 1)
    traj = simulate(compiled.process_set, params, compiled.y0, compiled.t_span_h, t_eval=grid)
    assert traj.success, traj.message
    total = np.asarray(traj.y, dtype=float)[compiled.schema.slice("S"), :].sum(axis=0)
    s0 = float(total[0])
    return {d: float((s0 - total[d * 24]) / s0) for d in range(days + 1)}


def test_matching_tyrells_extract_schedule_breaks_the_attenuation_benchmark():
    """Why D-215's extract xfail cannot be closed on the uptake rate (D-216 §4).

    The two anchors on beer's fermentation speed are not compatible under this rate law:

    * **Tyrell's measured extract course** — their wort is 59.4 % fermented by day 2, where the
      model books 21.2 % (``test_the_model_ferments_tyrells_wort_on_tyrells_schedule``);
    * **§2.2's acceptance criterion** — a 1.048 ale wort at 20 °C reaching 1.010 apparent in
      5-7 days (``test_beer_1048_og_attenuates_in_5_to_7_days``).

    Beer's uptake is ``q_sugar_max · X · Monod(S)``, so one constant scales both. The value that
    lands Tyrell's day 2 is 1.397 — **inside** the printed 0.3-1.5 band — and it takes the
    benchmark to 2.71 d. The window is already violated at q ≈ 0.6, having closed under a fifth
    of the gap.

    **The baseline is asserted first and deliberately.** A test that only showed the re-rated
    arm failing would not distinguish "this override breaks the benchmark" from "the benchmark
    is broken" [[feedback-pair-the-red-with-an-ordering-preserving-baseline]].
    """
    shipped = _beer_days_to_target_gravity()
    assert 5.0 <= shipped <= 7.0, (
        f"the benchmark wort attenuates in {shipped:.2f} d at the shipped q_sugar_max, outside "
        "§2.2's 5-7 d window. The control failed, so nothing below is attributable to the "
        "override — fix this before reading the arm"
    )

    matched = _beer_days_to_target_gravity(Q_SUGAR_MAX_MATCHING_TYRELL)
    assert matched < 5.0, (
        f"re-rated to the q_sugar_max that reproduces Tyrell's day-2 extract "
        f"({Q_SUGAR_MAX_MATCHING_TYRELL}), the benchmark wort attenuates in {matched:.2f} d, "
        "which is INSIDE §2.2's window. D-216 measured 2.71 d. If this is now inside, the two "
        "anchors no longer conflict and D-215's extract xfail is closable on this knob — which "
        "is a result, not a test failure: re-open D-216 §4"
    )


def test_removing_catabolite_repression_entirely_still_misses_tyrells_schedule():
    """The refusal's second tier: not even the unbounded limit of the dominant term (D-216 §5).

    ``K_repression`` = 2.0 g/L is tier ``speculative``, source *"author estimate"* — the
    functional form is Gee & Ramirez's but ``beer_generic.yaml`` records that their numeric
    constants were not accessible in-source. With glucose at 12.3 g/L it holds maltose at 14 %
    of its rate on day 0 and maltotriose at 0.5 %, so the model's day 1 is essentially
    glucose-only. It is by far the largest single contributor to the early-limb lag: removing it
    takes day 2 from 0.212 to 0.514, **79 % of the gap**, and it has the right SHAPE (a brake
    that vanishes as glucose clears, matching a lag that peaks at day 2 and closes by day 7).

    So the obvious next move is to re-source that constant. **This test says it would not be
    enough.** Removed ENTIRELY — not a candidate value, the limit — the model still falls short
    of Tyrell's day 2, while putting the benchmark at 3.42 d. The refusal in D-216 §6 is
    therefore not "no in-band point works" but the stronger "not even the unbounded limit of the
    term that owns most of it".
    """
    unrepressed = _tyrell_flux_fraction(K_repression=K_REPRESSION_REMOVED)
    measured_day2 = TYRELL_FLUX_FRACTION[2]
    assert unrepressed[2] < measured_day2 - 0.05, (
        f"with catabolite repression removed entirely the model ferments {unrepressed[2]:.1%} of "
        f"Tyrell's wort by day 2 against their measured {measured_day2:.1%}. D-216 measured "
        "51.4 %. If the limit now reaches the measurement, the lag IS the repression term and "
        "re-sourcing K_repression becomes the beat — see D-216 §5"
    )

    shipped = _tyrell_flux_fraction()
    closed = (unrepressed[2] - shipped[2]) / (measured_day2 - shipped[2])
    assert closed == pytest.approx(0.79, abs=0.05), (
        f"removing repression closes {closed:.0%} of the day-2 gap; D-216 measured 79 %. This "
        "share is the reason the term is named as the dominant contributor, so it is pinned "
        "rather than left as prose"
    )


def test_the_beer_ph_agreement_is_conditional_on_the_scenario_pitch():
    """The pH course's 7/8 is scored on a pitch nothing independently sources (D-216 §8).

    D-207/208/209/211 all score beer's pH on ``TYRELL_SCENARIO``, which pitches **1.0 g/L**
    against Tyrell's stated 9.96e6 cells/mL — ``beer_generic.yaml`` already records that this
    implies ~100 pg dry weight per cell, about **2× the textbook 40-60 pg**. That note reads as
    a caveat about a scenario detail. It is not: two published results are conditioned on it.

    At a textbook-honest 0.5 g/L the day-1 miss D-211 pinned at **0.070 becomes 0.354** — five
    times worse — day 2 falls out of its envelope as well, and D-211's measured attribution (the
    nitrogen fraction drawn by 24 h, 0.363, "inside Tyrell's measured 0.234-0.448 spread") drops
    to 0.181, **outside** that spread.

    **This is not a claim that the pitch is wrong.** Read the other way, 1.0 g/L is the value at
    which the model's biomass reproduces Tyrell's measured growth *timing* in the units that
    drive the rate, and two independent observables endorse it against the per-cell arithmetic.
    The extract lag gets **worse**, not better, at the honest pitch — 2.81× → **5.31×** at
    0.5 g/L. (D-216 §7 printed 3.51× there; that is wrong and D-219 corrects it. It is prose
    only and was never asserted, and D-218 §5's own day-2 fraction at this pitch, 0.112,
    gives 0.594/0.112 = 5.31×. At Tyrell's counted pitch of 0.398 g/L it is 6.58×.)
    What this test forbids is reading D-211's 0.070 as unconditional.
    """
    compiled = compile_scenario(
        Scenario(
            name="d216-pitch",
            medium="beer",
            initial={**TYRELL_SCENARIO, "pitch_gpl": 0.5},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=14.0,
        )
    )
    res = compiled.run()
    params = compiled.parameters.resolve()
    states = np.asarray(res.y, dtype=float)
    t_h = np.asarray(res.t, dtype=float)

    def ph_at(day: int) -> float:
        y = np.array([np.interp(day * 24.0, t_h, states[i, :]) for i in range(states.shape[0])])
        return float(acidbase.degassed_ph_of_state(y, compiled.schema, params))

    day1_excess = ph_at(1) - TYRELL_PH_COURSE[1][1]
    assert day1_excess == pytest.approx(0.354, abs=0.03), (
        f"at a 0.5 g/L pitch day 1 reads {day1_excess:.4f} pH above Tyrell's envelope; D-216 "
        "measured 0.354 against the 0.070 the shipped 1.0 g/L pitch gives. If these have "
        "converged, the pH result no longer depends on the pitch and D-216 §8 is spent"
    )

    inside = [
        band_lo - TYRELL_PH_READ_TOL <= ph_at(day) <= band_hi + TYRELL_PH_READ_TOL
        for day, (band_lo, band_hi) in TYRELL_PH_COURSE.items()
    ]
    assert sum(inside) == 6, (
        f"{sum(inside)}/8 days inside at the honest pitch; D-216 measured 6, against 7 at the "
        "shipped pitch. The count is what says the agreement is conditional"
    )


# ======================================================================================
# 11. The temperature sensitivity of uptake — the lever D-216 named, and what holds it
#     open (decision D-217)
# ======================================================================================


def test_the_uptake_activation_energy_is_inert_at_the_attenuation_benchmark():
    """§2.2's benchmark cannot see `E_a_uptake` at all, and that is exact (D-217 §1).

    D-216 §6 named `E_a_uptake` as the only lever that decouples beer's two speed anchors,
    on the grounds that the benchmark runs at exactly ``T_ref`` = 20 °C so its Arrhenius
    factor is 1.0 by construction. That was inherited from a note in ``beer_generic.yaml``
    and spot-checked at one band edge. This asserts it across the whole band and well
    outside it, on both signs, because the argument the archive now rests on is not "the
    benchmark barely moves" but "the benchmark does not move".

    The arm at −97,000 J/mol is de Andres-Toro's fitted sugar-uptake term, the one figure
    that would nearly close Tyrell's gap; it is swept here to show that even a value of the
    opposite SIGN leaves this anchor untouched. Nothing here endorses it.

    **A RED names a change to the frame, not to a rate.** The only way this fails is if the
    benchmark stops running at ``T_ref``, or uptake stops being Arrhenius in it — either of
    which silently converts D-216's refusal from "two anchors conflict" into "one knob
    moves both", and would have to be priced before any of section 10 is read again.
    """
    baseline = _beer_days_to_target_gravity()
    assert 5.0 <= baseline <= 7.0, (
        f"the benchmark wort attenuates in {baseline:.2f} d at the shipped parameters, outside "
        "§2.2's 5-7 d window. The control failed, so nothing below is attributable"
    )

    for e_a in (-97000.0, 0.0, 30000.0, 55100.0, 63000.0, 80000.0):
        got = _beer_days_to_target_gravity(E_a_uptake=e_a)
        assert got == baseline, (
            f"E_a_uptake = {e_a:.0f} J/mol moves the benchmark to {got:.4f} d against a "
            f"baseline of {baseline:.4f} d. D-216 §6's decoupling argument requires this to "
            "be EXACTLY zero — the benchmark is supposed to sit at T_ref, where the Arrhenius "
            "factor is 1.0 by construction"
        )


def test_the_uptake_activation_energy_is_a_lever_only_because_the_trial_ran_cool():
    """The lever's whole size is Tyrell's distance from ``T_ref`` (D-217 §4).

    `E_a_uptake` is free at the benchmark (test above), so it is the one parameter that can
    move Tyrell without moving §2.2. That is true only because Tyrell's tube trial ran
    *below* 20 °C. Had it run at ``T_ref``, both anchors would sit at an Arrhenius factor of
    1.0 and the lever would not exist — sweeping the entire printed band would move nothing.

    Both arms are asserted because one alone would mislead. The 15 °C arm on its own reads
    as "there is a lever"; the 20 °C arm on its own reads as "there is no lever". Together
    they say what is true: the lever is worth exactly what the frame is worth, and the frame
    rests on §3.2's comparative sentence (see :data:`TYRELL_TRIAL_CELSIUS`), not on a printed
    fermentation temperature.

    The span is small either way. Across the whole 30,000-63,000 band the day-2 fraction moves
    0.045, taking D-215's 2.81× gap to 2.41× — under a fifth of it — which is why D-216
    refused the low edge on provenance rather than on power.
    """

    def day2_span(celsius: float) -> tuple[float, float]:
        lo = _tyrell_flux_fraction(celsius=celsius, E_a_uptake=30000.0)[2]
        hi = _tyrell_flux_fraction(celsius=celsius, E_a_uptake=63000.0)[2]
        return lo, hi

    cool_lo, cool_hi = day2_span(TYRELL_TRIAL_CELSIUS)
    span_cool = cool_lo - cool_hi
    assert span_cool == pytest.approx(0.0449, abs=0.005), (
        f"at {TYRELL_TRIAL_CELSIUS:.0f} °C the printed E_a_uptake band moves Tyrell's day-2 "
        f"fraction by {span_cool:.4f}; D-217 measured 0.0449. If this has collapsed, the lever "
        "D-216 §6 named is gone and its refusal needs re-reading"
    )

    ref_lo, ref_hi = day2_span(20.0)
    assert ref_lo == ref_hi, (
        f"at T_ref the same band gives {ref_lo:.9f} and {ref_hi:.9f}. These must be identical: "
        "at 20 °C the Arrhenius factor is 1.0 whatever E_a is, so the lever exists ONLY because "
        "Tyrell's trial ran cooler — which the paper states comparatively (§3.2), never as a "
        "printed temperature"
    )


# ======================================================================================
# 12. The third-party yardstick for beer's SPEED, and the unsourced conversion that
#     decides between it and the brief (decision D-218)
# ======================================================================================

#: Foster et al. 2022, *Kveik Brewing Yeasts Demonstrate Wide Flexibility in Beer
#: Fermentation Temperature Tolerance and Exhibit Enhanced Trehalose Accumulation*,
#: Front. Microbiol. 13:747546, doi:10.3389/fmicb.2022.747546 (PMC8966892).
#:
#: **Why this trial and not another.** D-216 §11 asked which of beer's two speed anchors to
#: calibrate against — §2.2's 5-7 d acceptance criterion or Tyrell's measured course — and
#: D-217 §8 left it open as the owner's call. This is the first third-party source the
#: archive has that carries the whole tuple the question needs: a wort within 0.003 SG of
#: §2.2's own, a **counted** pitch, two temperatures for the **same three strains**, and the
#: **same** 1.010 target gravity (the paper takes it from Parker 2008, so the target §2.2
#: asserts is independently sourced — that much of the brief is corroborated).
#:
#: *"The hopped wort was prepared using Canadian 2-row malt to an original gravity of
#: 12.5°Plato (1.045 specific gravity)"*; *"inoculated at a rate of 1.2 x 10^7 cells/mL into
#: 200 mL of sterilized wort in 250 mL glass bottles fitted with airlocks"*.
FOSTER_OG_SG = 1.045

#: Foster's fermentable extract in this file's units. 12.5 °P is ~128 g/L total extract at
#: the same all-malt apparent-attenuable fraction §2.2's own wort is built on, which is what
#: makes the two comparable at all.
FOSTER_FERMENTABLE_GPL = 92.0

#: The counted pitch. This is the number the whole beat turns on, because the engine works
#: in g/L and the paper works in cells/mL — see :data:`PER_CELL_DRY_MASS_PG`.
FOSTER_PITCH_CELLS_PER_ML = 1.2e7

#: Gravity was sampled at these hours only (Fig. 2 panels A-D). **This is why Foster's
#: "3 days" is a CEILING and not a duration:** the prose reads *"The Beer 1 control strains
#: reached SG < 1.01 after only 3 days ... (Figure 2C)"* and Figure 2C **is** the 72 h panel.
#: The strains were below target when someone next looked, not at 72 h exactly. If they were
#: still above at 48 h — which the word *"only"* implies but does not state — the true
#: endpoint is somewhere in (48, 72] h. Every use of it below is one-sided and says so.
FOSTER_SAMPLE_HOURS = (12.0, 48.0, 72.0, 120.0)

#: *"reached SG < 1.01 after only 3 days"* — at 22 °C, for the three **Beer 1** controls
#: (Cali Ale, Vermont Ale, Kölsch; all *S. cerevisiae* ale strains). St. Lucifer is the
#: Beer 2 control and *"never reached an FG < 1.01, regardless of temperature"*, so it is
#: excluded from this figure and from the one below — the pair really is same-strain.
FOSTER_DAYS_TO_FG_AT_22C = 3.0

#: *"Hornindal1 was the only kveik strain that completed a 12°C fermentation within 10 days
#: along with the Beer 1 controls"*. 10 days is the **incubation length** (*"Incubations were
#: for 5 days, except for the 12°C fermentations (10 days)"*), so this is a ceiling by
#: construction, not a measurement.
FOSTER_CEILING_DAYS_AT_12C = 10.0

#: The tightest reading of Foster's within-source duration ratio, and the loosest.
#:
#: Both endpoints above are ceilings, and that matters for the DIRECTION of the bound. An
#: upper bound on ``d12 / d22`` needs an upper bound on ``d12`` (10 d, have it) and a LOWER
#: bound on ``d22`` (do not have it). Taking 3 d as exact gives 10/3; taking the open end of
#: the sampling interval, ``d22 > 2``, gives 10/2. The truth is in between and nothing in
#: the paper narrows it further.
FOSTER_RATIO_BOUND_TIGHTEST = FOSTER_CEILING_DAYS_AT_12C / FOSTER_DAYS_TO_FG_AT_22C
FOSTER_RATIO_BOUND_LOOSEST = FOSTER_CEILING_DAYS_AT_12C / 2.0

#: Yeast dry mass per cell, pg — the conversion that turns Foster's counted pitch into a
#: ``pitch_gpl``, and (D-218 §3) the quantity that decides beer's whole speed question.
#:
#: **SETTLED AT D-219 at 40 pg/cell, and neither reading the archive shipped survives.** The
#: value is not a literature pick. Coleman, Fish & Block 2007 — the paper this engine's wine
#: growth, death and yield parameters are fitted to — state verbatim that *"each cell count
#: was converted to grams per liter of cell mass, assuming that each cell weighs 4 x 10^-11
#: g"*. Every gram in that paper is a hemacytometer count times that constant, so 40 pg/cell
#: is the DEFINITION of the gram this engine's biomass state variable is expressed in, and a
#: pitch converted at anything else is in a unit the model's own parameters do not use. It
#: lives in :func:`fermentation.units.cells_per_ml_to_pitch_gpl`, which carries the quote,
#: the independent elemental corroboration (34.9 pg) and the band (28-50 pg).
#:
#: The dict is kept as the BRACKET the sweep below runs over — the readings are what make the
#: verdict pattern legible — not as a live fork. 18 pg was an unsourced assertion in the two
#: wine benchmarks (corrected at D-219); ~100 pg was never chosen by anyone, being
#: back-computed from :data:`TYRELL_SCENARIO`'s 1.0 g/L against Tyrell's counted 9.96e6
#: cells/mL, so it is a residual absorbing the true cell mass *and* every error in the
#: model's per-gram uptake rate. Both sit outside the settled 28-50 pg band.
PER_CELL_DRY_MASS_PG: dict[str, float] = {
    "retired: unsourced wine-benchmark assertion": 18.0,
    "SETTLED (Coleman 2007, D-219)": 40.0,
    "textbook high": 60.0,
    "retired: back-computed from the beer scenario pitch": 100.0,
}

#: ``q_sugar_max`` that lands the model exactly on Foster's 72 h endpoint at 22 °C, per
#: reading of the per-cell dry mass. Bisected on the shipped model on the hourly grid, and
#: verified live by the test below rather than trusted.
#:
#: **18 pg/cell is deliberately NOT in here** — see
#: :data:`Q_SUGAR_MAX_UNREACHABLE_AT_18_PG`. Every value in this dict is a crossing; a
#: saturation in the same container would be a value a later beat could iterate over and
#: read as fitted, which is D-177's lesson at the level of the data structure rather than
#: the prose.
Q_SUGAR_MAX_REACHING_FOSTER: dict[float, float] = {40.0: 0.9242, 50.0: 0.8176, 100.0: 0.5602}

#: At 18 pg/cell — the repo's own wine conversion — **no in-band value reaches Foster's
#: endpoint at all.** The printed band's CEILING of 1.5 still needs 3.04 d against a
#: published ≤3, missing by 1.4 %. This is a saturation, not a crossing.
Q_SUGAR_MAX_UNREACHABLE_AT_18_PG = 1.5

#: The 2 d arm on the 100 pg row: the ``q_sugar_max`` reaching Foster's target at the OPEN
#: end of its sampling interval (:data:`FOSTER_SAMPLE_HOURS`) rather than at 72 h. It is the
#: whole difference between "one corner of the bracket survives" and "none does", so it is
#: named here rather than left as a literal in the one test that reads it.
Q_SUGAR_MAX_REACHING_FOSTER_AT_2D_100PG = 0.8971

#: What §2.2's benchmark then reads, per reading. Only one lands inside 5-7 d.
HANDOFF_DAYS_AT_FOSTER_RATE: dict[float, float] = {
    18.0: 2.5833,
    40.0: 3.6250,
    50.0: 4.0,
    100.0: 5.5,
}


def _foster_pitch_gpl(pg_per_cell: float) -> float:
    """Foster's counted pitch as a ``pitch_gpl``, at one reading of the per-cell dry mass."""
    return FOSTER_PITCH_CELLS_PER_ML * pg_per_cell * 1e-9


def _foster_days_to_target_gravity(
    pitch_gpl: float,
    celsius: float,
    *,
    days: float = 30.0,
    per_hour: int = 1,
    **overrides: float,
) -> float:
    """Days for Foster's 12.5 °P wort to reach 1.010 apparent, on a FIXED grid.

    Reuses the benchmark's own gravity construction and target, exactly as
    :func:`_beer_days_to_target_gravity` does, so the two are commensurable.

    ``per_hour`` sets the grid's resolution, and it is load-bearing rather than cosmetic
    (D-214's lesson, one turn on). The reported duration is quantised to ``1/(24*per_hour)``
    days, so a RATIO of two durations carries a quantum of roughly
    ``(1 + ratio) / (24 * per_hour * d_hot)``. On the hourly grid at the heavy pitch that is
    ~0.04 — larger than several of the effects §12 measures. Each caller below states the
    resolution its claim needs.
    """
    from tests.benchmarks.test_milestone1 import TARGET_FG_SG, _apparent_gravity_series

    s0 = FOSTER_FERMENTABLE_GPL
    compiled = compile_scenario(
        Scenario(
            name="d218-foster",
            medium="beer",
            initial={
                "glucose_gpl": s0 * 0.132 / 0.88,
                "maltose_gpl": s0 * 0.546 / 0.88,
                "maltotriose_gpl": s0 * 0.202 / 0.88,
                "yan_mgl": 200.0,
                "pitch_gpl": pitch_gpl,
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=celsius)],
            duration_days=days,
        )
    )
    params = dict(compiled.param_values)
    for name, value in overrides.items():
        assert name in params, f"{name} is not a compiled parameter"
        params[name] = value
    steps = int(compiled.t_span_h[1] * per_hour)
    grid = np.linspace(0.0, compiled.t_span_h[1], steps + 1)
    traj = simulate(compiled.process_set, params, compiled.y0, compiled.t_span_h, t_eval=grid)
    assert traj.success, traj.message
    apparent = _apparent_gravity_series(traj, FOSTER_OG_SG, s0)
    reached = np.where(apparent <= TARGET_FG_SG)[0]
    return float(traj.t[reached[0]] / 24.0) if reached.size else float("inf")


#: The grid the temperature guards run on. Six minutes, not an hour, and the difference
#: decides whether two of the three claims below are measurements or noise: at the heavy
#: pitch the hourly grid's quantum on a duration RATIO is ~0.04, and the residual's distance
#: from 1.0 at the band's high edge is 0.002.
FOSTER_GRID_PER_HOUR = 10


@pytest.fixture(scope="module")
def foster_temperature_sweep() -> dict[tuple[float, float], float]:
    """(pg/cell, ``E_a_uptake``) -> the model's 12 °C / 22 °C duration ratio on Foster's wort.

    Module-scoped because three tests read it and each entry costs two integrations.

    **Both per-cell readings are swept, and which one each claim is read on matters.** The
    bound in §3 is cleared by ~0.9 at the heavy pitch, hundreds of times the grid quantum, so
    it is safe anywhere. The residual crossing is not: at 1.2 g/L the residual only reaches
    0.998 at the band's high edge — a margin of 0.002 against a 6-minute quantum of 0.004, so
    the crossing is **not resolved at that pitch**. At 0.216 g/L the same edge reads 0.987
    against a quantum of 0.002, and it is. The crossing test therefore runs on the light
    pitch; it is the pitch where the claim can be measured rather than the pitch that is
    cheapest.
    """
    out = {}
    for pg in (18.0, 100.0):
        pitch = _foster_pitch_gpl(pg)
        for e_a in (30000.0, 55100.0, 63000.0, 90000.0):
            if pg == 18.0 and e_a == 90000.0:
                continue  # only the bound test reads the control arm, and it runs heavy
            hot = _foster_days_to_target_gravity(
                pitch, 22.0, days=12.0, per_hour=FOSTER_GRID_PER_HOUR, E_a_uptake=e_a
            )
            cold = _foster_days_to_target_gravity(
                pitch, 12.0, days=25.0, per_hour=FOSTER_GRID_PER_HOUR, E_a_uptake=e_a
            )
            out[(pg, e_a)] = cold / hot
    return out


def test_fosters_temperature_pair_cannot_discriminate_anywhere_in_the_uptake_band(
    foster_temperature_sweep,
):
    """Foster's two temperatures do NOT test this model's temperature response (D-218 §2).

    This is the guard that stops the beat's most quotable line from being read as a pass.
    Foster gives 3 d at 22 °C and <=10 d at 12 °C for the same three strains, the same wort
    and the same pitch — a within-source ratio, immune to both the per-cell-mass fork and any
    species confound, and the only such pair the archive has ever had. The model's ratio is
    2.21, comfortably inside it, and an earlier draft of D-218 reported exactly that.

    **It means nothing, and this test is why.** The bound is cleared at *every* printed value
    of ``E_a_uptake``, from 1.59 at the low edge to 2.44 at the high — and cleared against the
    TIGHTEST reading of the bound, the one that treats Foster's 3 d as exact. A test no
    in-band configuration can fail is not evidence about the configuration that ships; D-216
    §10 and D-217 §6 both record that an out-of-band arm tests nothing, and this is the same
    fact stated about the *literature* rather than about a mutation.

    The 90,000 J/mol arm is the **positive control**: it is out of band (1.43x the high edge)
    and it DOES fire, which is what proves the predicate can distinguish anything at all.
    Without it, "everything clears" would be indistinguishable from a broken predicate.

    A RED on the in-band arms names the temperature response gaining reach — the bound
    becoming informative, which would be news and would license re-opening D-217's refusal.
    A RED on the control names the opposite: the predicate has gone slack.
    """
    baseline = _beer_days_to_target_gravity()
    assert 5.0 <= baseline <= 7.0, (
        f"the benchmark wort attenuates in {baseline:.2f} d at the shipped parameters, outside "
        "§2.2's 5-7 d window. The control failed, so nothing below is attributable"
    )

    for e_a in (30000.0, 55100.0, 63000.0):
        ratio = foster_temperature_sweep[(100.0, e_a)]
        assert ratio < FOSTER_RATIO_BOUND_TIGHTEST, (
            f"E_a_uptake = {e_a:.0f} J/mol gives a 12/22 °C duration ratio of {ratio:.3f}, which "
            f"now EXCEEDS Foster's tightest bound of {FOSTER_RATIO_BOUND_TIGHTEST:.2f}. That "
            "would make the pair discriminating for the first time — D-218 §2's whole point is "
            "that it is not, and D-217's E_a refusal was left standing on that basis"
        )

    fires = foster_temperature_sweep[(100.0, 90000.0)]
    assert fires > FOSTER_RATIO_BOUND_TIGHTEST, (
        f"the out-of-band control at 90,000 J/mol gives {fires:.3f}, which does NOT exceed "
        f"{FOSTER_RATIO_BOUND_TIGHTEST:.2f}. The predicate can no longer distinguish anything, "
        "so the 'every in-band value clears' result above is vacuous rather than measured"
    )
    assert fires < FOSTER_RATIO_BOUND_LOOSEST, (
        f"{fires:.3f} exceeds even the LOOSEST reading of the bound "
        f"({FOSTER_RATIO_BOUND_LOOSEST:.2f}, from d22 > 2 d). Recorded because it bounds what "
        "the sampling-grid ambiguity in Foster's 3 d could ever be worth: not enough to make an "
        "in-band value fire"
    )


def test_the_apparent_arrhenius_identity_in_beers_temperature_response_is_a_crossing(
    foster_temperature_sweep,
):
    """The ratio equals the bare Arrhenius factor at the nominal BY COINCIDENCE (D-218 §2).

    At the shipped ``E_a_uptake`` the model's cross-temperature duration ratio (2.212) sits
    within 0.7 % of ``exp[(E_a/R)(1/T_12 - 1/T_22)]`` = 2.198 — the Arrhenius factor with no
    model in it at all. That looks like a structural fact: duration inversely proportional to
    uptake rate, the ferment uptake-limited end to end, nothing else with a temperature
    dependence in the loop.

    **It is not.** The residual ``ratio / arrhenius`` runs 1.108 -> 1.012 -> 0.987 across the
    printed band and crosses 1.0 inside it. The nominal sits short of that crossing, so the
    near-identity is a property of one value, not of the rate law, and a beat that pinned it
    as a law would be pinning a coincidence — the same trap D-217 §2 named around -90,000
    J/mol.

    **Read at the LIGHT pitch, and that is not a free choice.** The crossing's location moves
    with the pitch (~58 kJ/mol at 0.216 g/L, ~62 kJ/mol at 1.2 g/L), and at the heavy pitch
    the residual only reaches 0.998 at the band's high edge — inside the grid's own quantum,
    so the crossing is unresolvable there however the assert is written. See
    :func:`foster_temperature_sweep`.

    Asserting the two ENDS is what makes the middle mean something. Without them,
    ``|residual - 1| < 0.02`` reads as "the model is a bare Arrhenius response"; with them it
    reads "the model is not, and this is where the two curves cross".
    """
    gas_constant = 8.314462618
    t_cold, t_hot = 273.15 + 12.0, 273.15 + 22.0

    def residual(e_a: float) -> float:
        arrhenius = float(np.exp((e_a / gas_constant) * (1.0 / t_cold - 1.0 / t_hot)))
        return float(foster_temperature_sweep[(18.0, e_a)]) / arrhenius

    low, nominal, high = residual(30000.0), residual(55100.0), residual(63000.0)

    assert low > 1.0 and high < 1.0, (
        f"the residual is {low:.4f} at the band's low edge and {high:.4f} at its high edge. "
        "D-218 §3 measured 1.108 and 0.987 — a sweep that CROSSES 1.0 inside the band. If both "
        "are now on one side there is no crossing, and the near-identity at the nominal below "
        "would be a structural property of the rate law rather than a coincidence"
    )
    assert min(low - 1.0, 1.0 - high) > 0.005, (
        f"the crossing's margins are {low - 1.0:.4f} and {1.0 - high:.4f}; D-218 measured 0.108 "
        "and 0.013. The 6-minute grid's quantum on this ratio is ~0.002, so a margin under "
        "0.005 is not resolved and the "
        "assert above would be reading the grid rather than the model"
    )
    assert abs(nominal - 1.0) < 0.02, (
        f"at the shipped E_a_uptake the residual is {nominal:.4f}; D-218 measured 1.0119. The "
        "nominal sitting within 1 % of the bare Arrhenius factor is what makes the model's "
        "temperature ratio LOOK like an identity, and §3's finding is that it is a crossing"
    )


def test_the_beer_temperature_ratio_is_pitch_invariant_only_near_the_nominal(
    foster_temperature_sweep,
):
    """Pitch-invariance of the temperature ratio is local, not structural (D-218 §2).

    Across a 5.6x pitch span — the full width of the per-cell-mass fork, 0.216 to 1.2 g/L —
    the 12/22 °C duration ratio moves by 0.007 at the shipped ``E_a_uptake`` and by 0.10 at
    the band's low edge. So the ferment is uptake-limited end to end *near the nominal* and
    measurably not at the low edge, where the growth phase gets long enough to matter.

    **This is recorded as a defect signature, not a pass.** Real ferments do show a
    pitch-dependent temperature response; a model whose response is pitch-blind at its shipped
    value is agreeing with the data for a reason the data does not have. It also bounds what
    §3's fork can do: the fork moves the pitch 5.6x, and this says that move cannot rescue the
    temperature response, only the durations.

    **The hourly grid could not have measured this and an earlier draft tried.** On it the
    nominal spread reads -0.0002, which is not a small number but a number below the quantum;
    on the 6-minute grid this fixture uses it reads +0.007. A guard asserting "under 0.02" on
    the hourly grid would have stayed green while the true spread grew past its own threshold.
    """
    spans = {
        e_a: abs(foster_temperature_sweep[(18.0, e_a)] - foster_temperature_sweep[(100.0, e_a)])
        for e_a in (30000.0, 55100.0)
    }

    assert spans[55100.0] < 0.02, (
        f"at the shipped E_a_uptake a 5.6x pitch change moves the temperature ratio by "
        f"{spans[55100.0]:.4f}; D-218 measured 0.007. If this has grown, the response has "
        "acquired a pitch dependence and §4's fork can no longer be read as moving durations "
        "only"
    )
    assert spans[30000.0] > 0.05, (
        f"the low band edge moves it by only {spans[30000.0]:.4f}; D-218 measured 0.10, against "
        f"{spans[55100.0]:.4f} at the nominal. That contrast is what makes the invariance LOCAL. "
        "If the edge has gone quiet too, the invariance is structural after all and the "
        "defect-signature reading in this docstring is wrong"
    )


def test_fosters_endpoint_breaks_the_handoff_window_at_every_reading_but_one():
    """§2.2's 5-7 d window survives in ONE corner, and it is the doubtful one (D-218 §3).

    D-216 §11 asked which anchor to calibrate beer's speed against. Foster answers it on the
    terms the brief itself sets — same target gravity, a wort within 0.003 SG — and the answer
    is *faster than the brief*, at every reading. What stops that from being a verdict is that
    the model needs Foster's cells/mL as a g/L, and the conversion is unsourced and internally
    inconsistent (:data:`PER_CELL_DRY_MASS_PG`).

    Each row: the ``q_sugar_max`` that puts the model on Foster's 72 h endpoint, and what
    §2.2's benchmark then reads.

    ======================  =========  ===============  ==========
    per-cell mass           pitch g/L  ``q_sugar_max``  §2.2 reads
    ======================  =========  ===============  ==========
    18 pg (retired)             0.216  1.5 = CEILING       2.58 d
    40 pg (SETTLED, D-219)      0.480  0.924             **3.63 d**
    50 pg                       0.600  0.818               4.00 d
    100 pg (retired)            1.200  0.560               5.50 d
    ======================  =========  ===============  ==========

    Only the last survives 5-7 d, and **D-219 retired it**: the settled conversion is 40
    pg/cell, the row above it, so §2.2's window does not survive. When D-218 wrote this the
    survivor was merely the doubtful corner; it is now the one reading known to be wrong.

    **And it survives only because Foster's 3 d is read as exact.** Take the open
    end of the sampling interval instead (d22 -> 2 d, still consistent with every word in the
    paper) and the 100 pg row goes to 3.71 d: the window then breaks in all eight bracket
    cells. That arm is asserted here because it is the whole difference between "one corner
    survives" and "nothing does".

    **What this test does NOT do is retire the window.** §2.2's benchmark passes today at
    6.08 d and is untouched. Retiring it means moving ``q_sugar_max``, which D-216 §4 already
    priced: no value of it reproduces Tyrell's shape, so a rate that satisfies Foster still
    misses the trajectory. A RED here names one of the pins moving, and the first thing to
    check is whether the per-cell mass finally got sourced.
    """
    from tests.benchmarks.test_milestone1 import TARGET_FG_SG

    assert TARGET_FG_SG == 1.010, (
        f"the benchmark's target gravity is {TARGET_FG_SG}; Foster's is 1.01 (they cite Parker "
        "2008). The comparison in this test is only sound while the two agree"
    )

    # The saturated row first: it is the one that is NOT a crossing, and keeping it out of
    # the loop below is what keeps that dict free of values a later beat could read as fitted.
    at_ceiling = _foster_days_to_target_gravity(
        _foster_pitch_gpl(18.0), 22.0, days=12.0, q_sugar_max=Q_SUGAR_MAX_UNREACHABLE_AT_18_PG
    )
    assert at_ceiling > FOSTER_DAYS_TO_FG_AT_22C, (
        f"at 18 pg/cell the band CEILING q = {Q_SUGAR_MAX_UNREACHABLE_AT_18_PG} now reaches "
        f"Foster's endpoint in {at_ceiling:.4f} d. D-218 measured 3.0417 — a 1.4 % miss, which "
        "is what makes this reading unreachable in band rather than merely expensive. If it is "
        "reachable, 1.5 stops being a saturation and the row needs bisecting for a real crossing"
    )
    verdicts = {
        18.0: 5.0
        <= _beer_days_to_target_gravity(q_sugar_max=Q_SUGAR_MAX_UNREACHABLE_AT_18_PG)
        <= 7.0
    }

    for pg, q in Q_SUGAR_MAX_REACHING_FOSTER.items():
        pitch = _foster_pitch_gpl(pg)
        got = _foster_days_to_target_gravity(pitch, 22.0, days=12.0, q_sugar_max=q)
        assert got == pytest.approx(FOSTER_DAYS_TO_FG_AT_22C, abs=0.05), (
            f"q = {q} at {pg:.0f} pg/cell reaches Foster's endpoint in {got:.4f} d, not "
            f"{FOSTER_DAYS_TO_FG_AT_22C}. The pinned crossing has moved, so every §2.2 "
            "reading below is against the wrong rate"
        )

        benchmark = _beer_days_to_target_gravity(q_sugar_max=q)
        assert benchmark == pytest.approx(HANDOFF_DAYS_AT_FOSTER_RATE[pg], abs=0.05), (
            f"at {pg:.0f} pg/cell the Foster-matching rate takes §2.2 to {benchmark:.4f} d; "
            f"D-218 measured {HANDOFF_DAYS_AT_FOSTER_RATE[pg]}"
        )
        verdicts[pg] = 5.0 <= benchmark <= 7.0

    assert verdicts == {18.0: False, 40.0: False, 50.0: False, 100.0: True}, (
        f"the window's survival per per-cell reading is now {verdicts}. D-218 measured exactly "
        "one survivor, at the 100 pg reading beer_generic.yaml calls ~2x textbook. A change here "
        "is a change to the ANSWER of D-216 §11, not to a number"
    )

    # The open end of Foster's sampling interval, on the one row that survived above.
    reached = _foster_days_to_target_gravity(
        _foster_pitch_gpl(100.0),
        22.0,
        days=12.0,
        q_sugar_max=Q_SUGAR_MAX_REACHING_FOSTER_AT_2D_100PG,
    )
    assert reached == pytest.approx(2.0, abs=0.05), (
        f"the 2 d arm reaches Foster's target in {reached:.4f} d, not 2.0 — re-bisect before "
        "reading the benchmark below"
    )
    open_end = _beer_days_to_target_gravity(q_sugar_max=Q_SUGAR_MAX_REACHING_FOSTER_AT_2D_100PG)
    assert not 5.0 <= open_end <= 7.0 and open_end == pytest.approx(3.71, abs=0.05), (
        f"reading Foster's endpoint at the open end of its sampling interval takes §2.2 to "
        f"{open_end:.4f} d; D-218 measured 3.71, outside 5-7. If this is now inside the window, "
        "the single surviving corner in the table above is no longer conditional on treating a "
        "72 h SAMPLE as an exact duration, and D-218 §3's central caveat is void"
    )


# ======================================================================================
# 13. D-219 — what the settled per-cell dry mass does to beer's two anchors
#
# The conversion itself, its corroboration and its band are guarded in tests/test_units.py.
# What is guarded HERE is the part that needs the engine: the biomass the beer scenario
# carries over Tyrell's counted pitch, and the verdict on §2.2's window once the fork
# §12 swept is down to one reading.
# ======================================================================================

#: Tyrell's counted pitch as a ``pitch_gpl``, at the settled conversion. Not a literal:
#: it comes through the same boundary function every other counted pitch does.
TYRELL_COUNTED_PITCH_GPL = cells_per_ml_to_pitch_gpl(9.96e6)

#: What the scenario ships, over what Tyrell counted. D-219's headline for the beer side.
TYRELL_SCENARIO_BIOMASS_EXCESS = 2.51

#: At Tyrell's counted pitch: the day-2 extract fraction, and the nitrogen drawn by 24 h.
#: Both measured at D-219 on the hourly grid §12's flux helper uses.
TYRELL_AT_COUNTED_PITCH = {"day2_fraction": 0.0903, "n_drawn_24h": 0.1446}

#: Tyrell's measured 24 h cell-count spread — the interval D-211 calls "what makes the
#: timing MEASURED rather than fitted", and which the counted pitch falls outside.
TYRELL_N_DRAWN_SPREAD = (0.234, 0.448)

#: The shipped model against Foster's endpoint at the SETTLED pitch, no knob touched.
#: D-218 §4's cleanest line was "3.33 d against a published <=3, an 11 % miss" — true, and
#: read on the 100 pg reading D-219 retires. These are the numbers that replace it.
FOSTER_AT_SETTLED_PITCH = {"d22": 4.8375, "d12": 10.7417}


def _tyrell_at_pitch(pitch_gpl: float, days: int = 7) -> tuple[dict[int, float], float]:
    """(daily fermented fraction, nitrogen fraction drawn by 24 h) at one pitch.

    Fixed hourly grid, matching :func:`_tyrell_flux_fraction`, so the day-2 number here is
    commensurable with the ones §12 and D-216 quote.
    """
    compiled = compile_scenario(
        Scenario(
            name="d219-tyrell-pitch",
            medium="beer",
            initial={**TYRELL_SCENARIO, "pitch_gpl": pitch_gpl},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=TYRELL_TRIAL_CELSIUS)],
            duration_days=float(days),
        )
    )
    grid = np.linspace(0.0, days * 24.0, days * 24 + 1)
    traj = simulate(
        compiled.process_set,
        dict(compiled.param_values),
        compiled.y0,
        compiled.t_span_h,
        t_eval=grid,
    )
    assert traj.success, traj.message
    y = np.asarray(traj.y, dtype=float)
    total = y[compiled.schema.slice("S"), :].sum(axis=0)
    s0 = float(total[0])
    frac = {d: float((s0 - total[d * 24]) / s0) for d in range(days + 1)}
    n = y[compiled.schema.slice("N"), :].sum(axis=0)
    return frac, float((n[0] - n[24]) / n[0])


def test_the_beer_scenario_carries_two_and_a_half_times_tyrells_counted_biomass():
    """``TYRELL_SCENARIO`` pitches 1.0 g/L; Tyrell counted 0.398 g/L (D-219).

    The ~100 pg/cell the archive carried since D-216 was never a conversion anyone chose.
    It is what you get by dividing this scenario's 1.0 g/L by Tyrell's counted
    9.96e6 cells/mL, so it is a RESIDUAL — it absorbs the true cell mass *and* every error
    in the model's per-gram uptake rate, which is why it came out three times any defensible
    cell mass. At the settled conversion the counted pitch is 0.398 g/L and the scenario
    carries **2.51x** the biomass Tyrell actually pitched.

    **The pitch is deliberately NOT corrected, and this test is where that refusal lives.**
    ``mu_max`` was fitted at 1.0 g/L against Tyrell's growth fraction, so moving the pitch
    inherits a ``mu_max`` refit AND a broken extract calibration in the same beat (D-218 §7).
    What ships is the measurement of what the excess is doing, which is large: the day-2
    extract shortfall goes 2.81x at the shipped pitch to **6.58x** at the counted one, and
    the nitrogen drawn by 24 h — D-211's own attribution for beer's pH timing — goes 0.360
    to **0.145**, outside Tyrell's measured 0.234-0.448 spread.

    Note this is the SAME direction D-216 §7 reported and a bigger number than it printed:
    that section says the 0.5 g/L shortfall is 3.51x where it is 5.31x (D-219 corrects it;
    prose only, never asserted). Its conclusion is strengthened, not weakened.

    A RED here is one of two things: the scenario pitch moved (in which case the ``mu_max``
    refit above is now owed), or the conversion moved.
    """
    excess = TYRELL_SCENARIO["pitch_gpl"] / TYRELL_COUNTED_PITCH_GPL
    assert excess == pytest.approx(TYRELL_SCENARIO_BIOMASS_EXCESS, abs=0.01), (
        f"the beer scenario carries {excess:.3f}x Tyrell's counted biomass; D-219 measured "
        f"{TYRELL_SCENARIO_BIOMASS_EXCESS}. Either the scenario pitch or the settled "
        "conversion has moved, and the two have very different consequences"
    )

    counted, n_drawn = _tyrell_at_pitch(TYRELL_COUNTED_PITCH_GPL)
    shipped, n_shipped = _tyrell_at_pitch(TYRELL_SCENARIO["pitch_gpl"])

    assert counted[2] == pytest.approx(TYRELL_AT_COUNTED_PITCH["day2_fraction"], abs=0.005), (
        f"at Tyrell's counted pitch the model ferments {counted[2]:.4f} of the wort by day 2; "
        f"D-219 measured {TYRELL_AT_COUNTED_PITCH['day2_fraction']}"
    )
    shortfall = TYRELL_FLUX_FRACTION[2] / counted[2]
    assert shortfall == pytest.approx(6.58, abs=0.15), (
        f"the day-2 shortfall at the counted pitch is {shortfall:.2f}x; D-219 measured 6.58. "
        "This is the cost of the honest pitch and it is the reason the correction is a "
        "measurement rather than a repair"
    )
    assert shortfall > TYRELL_FLUX_FRACTION[2] / shipped[2], (
        "the lag no longer gets WORSE at the counted pitch. That reverses D-216 §7's "
        "direction and would mean the biomass excess is masking a fast model, not a slow one"
    )

    assert n_drawn == pytest.approx(TYRELL_AT_COUNTED_PITCH["n_drawn_24h"], abs=0.01), (
        f"nitrogen drawn by 24 h at the counted pitch is {n_drawn:.4f}; D-219 measured "
        f"{TYRELL_AT_COUNTED_PITCH['n_drawn_24h']}"
    )
    assert not (TYRELL_N_DRAWN_SPREAD[0] <= n_drawn <= TYRELL_N_DRAWN_SPREAD[1]), (
        f"nitrogen drawn by 24 h at the counted pitch, {n_drawn:.4f}, is back inside Tyrell's "
        f"measured {TYRELL_N_DRAWN_SPREAD} spread. D-211's attribution would then hold at the "
        "honest pitch too, and D-216 §8's conditionality would be void"
    )
    # The designed contrast: the SHIPPED pitch is what lands inside that spread. Without
    # this the assert above cannot tell "outside" from "the predicate never fires".
    assert TYRELL_N_DRAWN_SPREAD[0] <= n_shipped <= TYRELL_N_DRAWN_SPREAD[1], (
        f"at the shipped 1.0 g/L the nitrogen drawn by 24 h is {n_shipped:.4f}, already "
        "outside Tyrell's spread. The contrast this test is built on has gone"
    )


def test_the_handoff_window_does_not_survive_the_settled_conversion():
    """§2.2's 5-7 d window fails, and D-219 is what makes that a verdict (D-218 §4).

    §12's fork table swept four readings of the per-cell mass and found the window surviving
    in exactly one cell — the 100 pg reading. D-219 settles the conversion at **40 pg/cell**,
    Coleman's own, which is a different row of that table: the Foster-matching rate takes the
    benchmark to **3.63 d**. So the surviving corner is the retired reading, and the window
    does not survive.

    Two separate things are asserted, and the order matters:

    1. the settled reading is IN the swept bracket and its verdict is False — the sweep is
       not being re-run here, it is being pointed at;
    2. the one reading whose verdict is True is a RETIRED one.

    **The window is still not retired in code, and that is not deference.** Retiring it means
    moving ``q_sugar_max``, and D-216 §4 priced that: a rate satisfying Foster leaves the
    model finishing a day early against both measured tails while still missing both measured
    day 2s. Fitting an acceptance criterion the model passes today to a rate that fits the
    data worse in SHAPE is not a fidelity gain. The verdict is recorded; the knob is not moved.

    A RED here is a change to the ANSWER, not to a digit.
    """
    from tests.test_units import RETIRED_READINGS_PG, SETTLED_BAND_PG

    settled = [
        pg for pg in PER_CELL_DRY_MASS_PG.values() if SETTLED_BAND_PG[0] <= pg <= SETTLED_BAND_PG[1]
    ]
    assert settled == [40.0], (
        f"the swept bracket now has {settled} inside the settled band, not exactly [40.0]. "
        "The verdict below reads one row of §12's table and needs that row to be unambiguous"
    )
    assert 40.0 in Q_SUGAR_MAX_REACHING_FOSTER, (
        "the settled reading is no longer one of the bisected rows, so §12's table no longer "
        "says anything about the reading that ships"
    )

    benchmark = _beer_days_to_target_gravity(q_sugar_max=Q_SUGAR_MAX_REACHING_FOSTER[40.0])
    assert benchmark == pytest.approx(HANDOFF_DAYS_AT_FOSTER_RATE[40.0], abs=0.05)
    assert not (5.0 <= benchmark <= 7.0), (
        f"at the settled conversion the Foster-matching rate takes §2.2 to {benchmark:.4f} d, "
        "which is now INSIDE the 5-7 d window. That reverses D-218 §4 and the brief's window "
        "would be corroborated by a third-party trial rather than refuted by one"
    )
    survivors = {pg for pg, days in HANDOFF_DAYS_AT_FOSTER_RATE.items() if 5.0 <= days <= 7.0}
    assert survivors <= set(RETIRED_READINGS_PG.values()), (
        f"the readings at which §2.2's window survives are {sorted(survivors)}, and at least "
        "one of them is NOT a retired reading. D-219's point is that the only corner the "
        "window survives in belongs to a conversion now known to be wrong; if a live reading "
        "has joined it, this record's §4 needs rewriting rather than citing"
    )
    assert survivors, (
        "no reading in the bracket survives 5-7 d at all. The verdict is unchanged but the "
        "CONTRAST this test is built on is gone -- 'the survivor is the retired one' becomes "
        "vacuous, and the guard would then pass on a broken sweep"
    )


def test_at_the_settled_conversion_the_shipped_model_is_slow_against_foster():
    """The number that replaces D-218 §4's "11 % miss" (D-219).

    D-218 §4 closed on its cleanest line: *"at the archive's own 100 pg reading, with no knob
    touched at all, the shipped model reaches S1's endpoint in 3.33 d against a published
    <=3. An 11 % miss."* That was measured, and it was measured on the reading D-219 retires.

    At the settled 40 pg/cell — Foster's counted 1.2e7 cells/mL is **0.48 g/L**, not 1.2 —
    the shipped model takes **4.84 d** against the same published <=3, a **1.61x** miss, and
    it also misses the 12 C ceiling (10.74 d against <=10). Both halves of D-218's conclusion
    therefore need restating rather than one: the brief's 5-7 d window is still refuted by
    every third-party endpoint found, AND the model is ~1.6x slow on a third-party trial
    rather than nearly right.

    Read against §12's other guards this is consistent, not new: D-215 measured the same
    engine fermenting Tyrell's wort ~2.8x too slowly *at a pitch already carrying 2.51x the
    counted biomass*. Those two numbers COMPOUND — they are not two routes agreeing on one
    figure — which is why the day-2 shortfall at the counted pitch is 6.58x and not 2.8x.

    Both durations run on the same 6-minute grid §12's temperature guards use, because a
    duration read off an hourly grid carries a quantum larger than the third digit here.
    """
    pitch = _foster_pitch_gpl(40.0)
    assert pitch == pytest.approx(cells_per_ml_to_pitch_gpl(FOSTER_PITCH_CELLS_PER_ML)), (
        "§12's own pitch helper at the settled reading disagrees with the boundary "
        "conversion; one of the two has drifted"
    )

    d22 = _foster_days_to_target_gravity(pitch, 22.0, days=12.0, per_hour=FOSTER_GRID_PER_HOUR)
    d12 = _foster_days_to_target_gravity(pitch, 12.0, days=25.0, per_hour=FOSTER_GRID_PER_HOUR)
    assert d22 == pytest.approx(FOSTER_AT_SETTLED_PITCH["d22"], abs=0.02)
    assert d12 == pytest.approx(FOSTER_AT_SETTLED_PITCH["d12"], abs=0.05)

    miss = d22 / FOSTER_DAYS_TO_FG_AT_22C
    assert miss == pytest.approx(1.61, abs=0.03), (
        f"the shipped model reaches Foster's endpoint in {d22:.4f} d against a published "
        f"<={FOSTER_DAYS_TO_FG_AT_22C}, a {miss:.2f}x miss; D-219 measured 1.61. D-218 §4's "
        "11 % was the same quantity at the retired 100 pg reading"
    )
    assert d12 > FOSTER_CEILING_DAYS_AT_12C, (
        f"the 12 C arm now finishes in {d12:.4f} d, inside Foster's {FOSTER_CEILING_DAYS_AT_12C} "
        "d incubation ceiling. The cold arm was the one place the model was NOT clearly slow, "
        "and if it now clears, the slowness is not uniform across temperature and D-217's "
        "refusal on E_a_uptake is worth re-opening"
    )
