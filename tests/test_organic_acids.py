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
from pathlib import Path

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.core.acidbase import charge_balance_is_populated
from fermentation.core.chemistry import (
    M_NITROGEN,
    M_TARTARIC,
    carbon_mass_fraction,
    sugar_species,
)
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
from fermentation.parameters.store import ParameterSet
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units import cells_per_ml_to_pitch_gpl
from fermentation.validation import BENCHMARKS
from tests.conftest import BEER_COUNTED_PITCH_CELLS_PER_ML

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

#: Tyrell's pitch, as they state it: a COUNT. Yeast was propagated, harvested at
#: "Hochkräusen" and pitched by cell concentration (MEBAK III 10.4.3); the paper prints no
#: pitch by weight, so this is the only pitch statement it makes. Fig. 4's cell-count panel
#: reads 9.98-9.94 x10^6/mL at day 0 across the four strains — one wort, one pitching rate —
#: and D-211 transcribed it as 9.96e6.
TYRELL_PITCH_CELLS_PER_ML = 9.96e6

#: What ``TYRELL_SCENARIO`` carried from D-178 to D-222: a flat 1.0 g/L, 2.51x Tyrell's counted
#: biomass. Named rather than left as a literal because three tests read it as the RETIRED arm,
#: and a retired value that survives only inside test bodies is one a later beat re-invents.
TYRELL_SCENARIO_RETIRED_PITCH_GPL = 1.0

#: Tyrell's wort as this engine states it: the derived fermentable extract split in the
#: customary all-malt proportions, pitched at their 15 °C, anchored at their measured wort pH.
#:
#: **The pitch is Tyrell's own count, converted (D-222).** It carried a flat ``1.0`` from
#: D-178 until D-222 — 2.51x the biomass Tyrell pitched, and a figure nothing sourced: D-219
#: showed it is a RESIDUAL, back-computed cell mass (~100 pg/cell) rather than a chosen one,
#: most likely the commercial dry-yeast dosing convention of a gram of PRODUCT per litre. It
#: now comes through the same boundary conversion every other counted pitch in this repo does,
#: so every Tyrell comparison in this file — pH course, acid courses, extract schedule,
#: nitrogen timing — is scored at the trial's own inoculum.
#:
#: The correction is not free and the cost is recorded rather than absorbed: it inherits the
#: ``mu_max`` refit D-219 §7 priced (0.034 -> 0.058 /h, D-211's method unchanged), the day-2
#: extract shortfall goes 2.81x -> 4.21x, and day 7 attenuates 0.782 against a measured 0.997.
#: What the refit buys back is D-211's nitrogen timing (0.298, inside Tyrell's spread) and the
#: pH course (7 of 8 days inside, as before). See §13.
TYRELL_SCENARIO = {
    "glucose_gpl": 0.15 * TYRELL_SUGAR_GPL,
    "maltose_gpl": 0.70 * TYRELL_SUGAR_GPL,
    "maltotriose_gpl": 0.15 * TYRELL_SUGAR_GPL,
    "yan_mgl": 200.0,
    "pitch_gpl": cells_per_ml_to_pitch_gpl(TYRELL_PITCH_CELLS_PER_ML),
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


@pytest.mark.xfail(
    strict=True,
    reason=(
        "D-215: the model ferments this wort too slowly - ~2.8x at day 2 when measured, ~4.2x "
        "once D-222 put the scenario on Tyrell's own counted pitch, ~3.2x since D-223 "
        "re-anchored the uptake rate to Foster's measured course"
    ),
)
def test_the_model_ferments_tyrells_wort_on_tyrells_schedule():
    """The extract panel of Fig. 4, scored against the model for the first time (D-215).

    **Why this is an xfail and not a pin.** It states the thing that is true of the source and
    false of the model, so a correct fix turns it GREEN and nothing has to be deleted — the D-208
    idiom. A plain assert on the model's *current* fraction would encode the defect instead, which
    is exactly why D-207 shipped the pH course as data no assert read.

    The gap is not subtle. Tyrell's wort is **59.4 %** fermented by day 2; this engine booked
    **~21 %** when D-215 measured it and books **~14 %** since D-222 corrected the scenario's
    pitch from 1.0 g/L to Tyrell's own counted 0.398, and it does not reach dryness until about
    day 10 where their extract curve is flat by day 5. The tolerance below is ±0.10 of the
    fermentable — roughly five times any plausible read error on the extract panel — so this
    cannot be argued down to a transcription quibble.

    **This is upstream of every acid course.** The acids are produced as ``Y · ΔS``, so an extract
    curve that is too gradual makes every flux-linked acid too gradual with it. Fix this and the
    acid courses move without any rate law changing; fix the acid rate laws while this stands and
    they will be fitted to compensate for it [[feedback-a-margin-can-be-borrowed-from-a-defect]].

    **SCOPE: this is measured on ONE scenario** — Tyrell's wort, 15 °C, YAN 200, and (since
    D-222) their own counted pitch of 0.398 g/L where it was 1.0. It
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
    assert shipped_rmse / flux_rmse == pytest.approx(0.632, abs=0.03), (
        f"the shape-error ratio is {shipped_rmse / flux_rmse:.3f}; D-183 measured 0.528, D-211 "
        "0.624 after slowing growth, D-222 0.581 after correcting the scenario pitch to Tyrell's "
        "own count, and D-223 0.632 after re-anchoring `q_sugar_max` to Foster's measured course. "
        "A move here means the growth path or the uptake rate changed again"
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
    assert measured_share_day1 / model_share_day1 == pytest.approx(2.60, abs=0.15), (
        f"the model books {model_share_day1:.3f} of acetic's rise by day 1 against Tyrell's "
        f"{measured_share_day1:.3f}, a {measured_share_day1 / model_share_day1:.2f}x shortfall; "
        "D-211 measured 2.15x once the growth rate was corrected and D-222 2.60x once the "
        "scenario carried Tyrell's own counted pitch. The gap D-211 FLAGGED is wider, not "
        "narrower: nothing in the file explains acetic's early rise, and the honest pitch makes "
        "that more visible. If it has closed, a producer faster than growth has been built and "
        "D-183's growth-linked choice should be re-read"
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
        # D-239: the three free amino-acid side chains are in the wort too, and they are a
        # READING of the nitrogen pool rather than a seed — so they enter at the run's own
        # starting N and LEAVE with it. Both ends are built here rather than only the start,
        # because the whole content of the term is that the wort buffers and the finished beer
        # does not [[feedback-gate-both-halves-of-a-pair]].
        amino_start = acidbase.amino_buffer_from_gpl(nitrogen_gpl_start, "beer", params)
        amino_end = acidbase.amino_buffer_from_gpl(nitrogen_gpl_end, "beer", params)
        seeded |= amino_start
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
        end |= amino_end  # the pool the yeast ate is not in the finished beer (D-239)
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
    assert min(degassed_nominal) == pytest.approx(0.941, abs=0.005), (
        f"the degassed prediction moved to {min(degassed_nominal):.1%} of Tyrell's measured "
        "0.81 pH drop at nominal yields; D-239 measures 94.1-129.2 % across the pKa band, "
        "D-209 measured 91.4-127.1 % for the same arm before the wort's three free "
        "amino-acid side chains were in the balance, and D-208 43.2-62.9 % before the "
        "nitrogen pool's charge was. Every one of those beats moved this arm UP by adding "
        "something the wort really has; none of them fitted it. This is the number that "
        "compares with a published beer pH"
    )
    assert max(degassed_nominal) == pytest.approx(1.292, abs=0.005), (
        f"the degassed prediction's high edge moved to {max(degassed_nominal):.1%}; D-209 "
        "measures 127.1 %, so at this arm's DAY-14 endpoint the high edge of the peptide pKa "
        "band now OVERSHOOTS the measured drop. Two things keep that honest and neither is a "
        "tuning knob: Tyrell's measurement stops at day 7 while this arm reads day 14, where "
        "the model is still producing acid (the day-7 comparison is the acceptance test below, "
        "101.7-107.1 %); and z-bar is DERIVED from published wort composition, never fitted"
    )
    assert min(vessel_nominal) == pytest.approx(1.093, abs=0.005), (
        f"the IN-VESSEL fraction moved to {min(vessel_nominal):.1%}; D-209 measures 106.8 % "
        "(D-183's 77.8 %, D-182's 77.6 %, and D-181's 42.7 % was this same model with no "
        "dissolved CO2 in its charge balance). Pinned as a model property: no published beer pH "
        "is measured in this frame, so a change here is a change to the vessel's chemistry, not "
        "to an agreement"
    )
    assert max(vessel_nominal) == pytest.approx(1.401, abs=0.005), (
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
    assert min(degassed_joint) == pytest.approx(0.734, abs=0.02), (
        f"the degassed joint low corner moved to {min(degassed_joint):.1%}; D-209 measures "
        "71.2 %, against D-208's 8.3 % and D-181's pre-carbonic 7.6 %. The low corner moved "
        "further than the high one for the same reason it did at D-182: a member predicting "
        "little acidification finishes at a higher pH, and the nitrogen term is a fixed charge "
        "removal, so it buys the most where the acids buy the least"
    )
    assert max(degassed_joint) == pytest.approx(1.425, abs=0.02), (
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
    assert min(vessel_joint) == pytest.approx(0.950, abs=0.02), (
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
    assert max(vessel_joint) == pytest.approx(1.517, abs=0.02), (
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
            acidbase._totals_molar(y, compiled.schema, params),
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


#: The eight wort acids Peyer's BC titration is computed on top of (decision D-181).
_WORT_BC_SEEDS = {
    "lactic": "lactic_typical_wort",
    "acetic": "acetic_typical_wort",
    "citrate": "citric_typical_wort",
    "malic": "malic_typical_wort",
    "succinic": "succinic_typical_wort",
    "pyruvic": "pyruvic_typical_wort",
    "formic": "formic_typical_wort",
    "oxalic": "oxalic_typical_wort",
}


def _peyer_wort_bc(params: ParameterSet, capacity_gpl: float, pka_peptide: float) -> float:
    """Peyer's fast buffering capacity for a wort carrying ``capacity_gpl`` at ``pka_peptide``.

    His own protocol (thesis sec 5.3.4): 375 µL of 1 M HCl into 25 mL,
    ``BC = log10[H+ added / H+ increase]``, reproduced on this engine's charge balance.
    Factored out at D-233 so the round-trip below and the drawn-member guard beside it are
    provably the SAME titration — two hand-copies would let one drift and still look agreed.

    **Since D-238 the titration itself lives in ``src``** (:func:`acidbase.peyer_fast_bc`),
    because ``y0_for_member`` roots on it per ensemble member. This wrapper stays because it
    builds the sample from the *parameter table* while the production path builds it from the
    *state*, and a test asserts those two agree. The round-trip below keeps its teeth either
    way: its two anchors are the shipped YAML literal and Peyer's published 1.18, and the
    solver produces neither.
    """
    values = params.resolve()
    pka = dict(acidbase.build_pka_map(values))
    pka["peptide_buffer"] = (pka_peptide,)
    totals = {
        slot: params[p].value / acidbase.ALL_ACIDS[slot].molar_mass
        for slot, p in _WORT_BC_SEEDS.items()
    }
    totals["peptide_buffer"] = capacity_gpl / acidbase.ALL_ACIDS["peptide_buffer"].molar_mass
    # D-239: Peyer titrated a real wort, so the sample carries the three free amino-acid side
    # chains too — at the CALIBRATION wort's own assimilable nitrogen, which is a coordinate this
    # sample did not have before (the eight organic-acid seeds are nitrogen-blind). Read from the
    # parameter store rather than typed, because the shipped capacity is a root taken ON it.
    totals |= acidbase.amino_buffer_from_gpl(
        values["peyer_control_wort_yan_gpl"], "beer", values
    )
    return acidbase.peyer_fast_bc(totals, pka)


def _beer_acid_params() -> ParameterSet:
    data = default_data_dir()
    return load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
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

    **It scores the NOMINAL pKa only, and D-233 measured what that leaves uncovered** — see
    the guard directly below.
    """
    params = _beer_acid_params()
    bc = _peyer_wort_bc(
        params,
        params["peptide_buffer_capacity_beer"].value,
        params["pKa_peptide_buffer"].value,
    )

    assert bc == pytest.approx(1.18, abs=1e-9), (
        f"the shipped peptide capacity reproduces BC = {bc:.6f}, not Peyer's published 1.18. "
        "The wort acid table changed without the back-solve being re-run."
    )


def test_holding_the_capacity_fixed_while_the_pka_moves_costs_peyers_1_18():
    """The MECHANISM behind D-238's repair — and it was a defect pin that outlived its defect.

    This was ``test_a_drawn_peptide_pka_carries_a_wort_that_is_not_peyers_1_18``, written at
    D-233 in the idiom *"a RED means it was FIXED — delete this guard and say so in the
    record"*. **D-238 fixed it and this test stayed GREEN**, which was predicted in advance:
    every number here is the SHIPPED CONSTANT scored through this file's own titration, and the
    repair does not touch that constant. It re-roots the *seed* per ensemble member, on a path
    this test never drives. That is D-236 §3's lesson arriving a second time — a defect pin is
    a statement about which path its arms reach, not about whether the defect survives — so it
    is re-scoped and renamed rather than deleted, exactly as copper's was.

    What it pins now is why the pair has to move TOGETHER: hold the capacity at the value
    rooted for pKa 4.25 and let the pKa go anywhere else in its sourced window, and the wort
    stops reproducing the BC = 1.18 the capacity exists to reproduce. That is the argument for
    rule 3 in ``CompiledScenario.y0_for_member``, and if these numbers move, that argument has
    changed and the record's 0.0100 pH day-14 cost needs re-measuring.

    The shape matters and is asserted, not just the span: BC is maximal AT the nominal by
    construction and falls off on BOTH sides, so the low edge is not the only offender and a
    one-sided claim about the band would be wrong.

    The repair itself is pinned by
    :func:`test_every_sampled_member_carries_peyers_wort_bc`, which reads the member's SEEDED
    SLOT — the parameter still holds the nominal after the repair, so a guard reading the
    parameter would be vacuous.
    """
    params = _beer_acid_params()
    shipped = params["peptide_buffer_capacity_beer"].value
    unc = params["pKa_peptide_buffer"].uncertainty
    nominal_pka = params["pKa_peptide_buffer"].value

    at_low = _peyer_wort_bc(params, shipped, unc.low)
    at_nominal = _peyer_wort_bc(params, shipped, nominal_pka)
    at_high = _peyer_wort_bc(params, shipped, unc.high)

    assert at_nominal == pytest.approx(1.18, abs=1e-9)
    assert at_low == pytest.approx(1.120250, abs=1e-5), (
        f"holding the shipped capacity at pKa {unc.low} gives BC = {at_low:.6f}; D-239 measures "
        "1.120250 and D-233 measured 1.116059 before the wort carried its three free amino-acid "
        "side chains. This is the shipped CONSTANT, which D-238's repair does not touch — so a "
        "move here is the acid table, the pKa window, or the D-239 split, never the repair."
    )
    assert at_high == pytest.approx(1.148657, abs=1e-5), (
        f"a member drawing pKa {unc.high} carries a wort at BC = {at_high:.6f}; D-239 measures "
        "1.148657, D-233 measured 1.145594."
    )
    # Both edges moved TOWARD the nominal, and that direction is the D-239 split showing up
    # where it should: 7.3 % of the wort's buffering left the lump, so the lump's pKa has less
    # of the wort to mis-place. The pair is less sensitive than it was, not more — a move in the
    # other direction would mean the split had added a second pKa-driven channel rather than
    # taking one away.
    assert (1.18 - at_low) < (1.18 - 1.116059), (
        "the low-edge miss grew; the split should have shrunk the peptide pKa's leverage"
    )
    assert at_low < at_high < at_nominal, (
        "BC is maximal at the NOMINAL pKa by construction and falls off on both sides. That "
        "ordering is the reason D-214's '1.1161-1.180' is a span and not a direction: the low "
        "edge is the worst point, but the high edge is wrong too."
    )


# -- the capacity/pKa pair, made coherent per member (decision D-238) ---------------------------
#
# D-233 §8 declined this repair because moving the BC back-solve into `src` "would make the
# round-trip test above compare the root-finder against itself". That objection is against
# DERIVING THE SHIPPED CONSTANT, and D-238 does not: the YAML literal is untouched and still
# scored against Peyer's published 1.18 by a titration neither side produces. What moves is the
# per-member SEED.

#: A beer that STATES its buffering protein. `peptide_buffer_gpl` makes the pool a scenario INPUT,
#: which D-24 excludes from sampling, so rule 3 must not fire — the branch is invisible from the
#: repaired side, so it owes its own arm (D-236's Arm F lesson).
TYRELL_NAMED_PEPTIDE_GPL = 1.4


def _beer_ensemble(scenario_extra: dict[str, float] | None = None, n_members: int = 8):
    """One beer ensemble over the peptide pair, returned with its compiled scenario.

    ``only=["pKa_peptide_buffer"]`` — with a single name drawn, every member's wort differs from
    the nominal's in exactly the quantity under test, so the per-member BC below is a statement
    about the pair and not about eighty-two other draws.
    """
    compiled = compile_scenario(
        Scenario(
            name="d238-peptide-pair",
            medium="beer",
            initial={**TYRELL_SCENARIO, **(scenario_extra or {})},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=2.0,
        )
    )
    ens = compiled.run_ensemble(n_members=n_members, seed=0, only=["pKa_peptide_buffer"])
    assert ens.sampled_names == ("pKa_peptide_buffer",), f"drew {ens.sampled_names}"
    assert ens.n_succeeded >= 2, "need at least two members"
    return compiled, ens


def test_the_runtime_solver_reproduces_the_shipped_capacity_at_the_nominal():
    """The control, and the reason the repair re-derives the shipped root rather than competing.

    [[feedback-the-setting-where-a-change-is-exact-is-the-control]]. Two claims, and the second
    is what makes the first mean anything:

    * :func:`acidbase.peptide_capacity_for_wort_bc`, run on a compiled beer ``y0`` **at the
      calibration wort's own assimilable nitrogen**, at the nominal pKa, returns the shipped
      ``peptide_buffer_capacity_beer`` **bit for bit**. Not to a tolerance — a tolerance would
      pass on a solver that merely landed nearby, and the whole argument for rule 3 is that it
      reproduces the offline back-solve rather than replacing it. (D-233 §1 reported these one
      ULP apart; that was a looser root-find, not a floor.)

      **The nitrogen qualifier is new at D-239 and it is not a loosening.** Three of the wort's
      buffering species now scale with ``N``, so the calibration wort acquired a coordinate it
      never had: this control used to run on Tyrell's scenario because ANY wort was Peyer's
      wort. Tyrell's 200 mg/L is 4.5 % more nitrogen than Peyer's control wort carries, and it
      roots 0.35 % lower — pinned below as a measurement rather than absorbed, because a beat
      that "fixed" it by re-rooting the shipped literal at Tyrell's YAN would be calibrating one
      man's wort to another's nitrogen [[feedback-right-number-wrong-condition]].
    * the state-built sample and the parameter-built sample are the SAME wort. The production
      path reads the acid slots off ``y0``; this file's titration builds them from the
      ``*_typical_wort`` parameters. Two ways of naming one composition is how they drift.

    Note what is NOT asserted here: that the solver is right. That is
    :func:`test_the_peptide_capacity_still_reproduces_peyers_published_wort_bc`'s job, and it
    scores the YAML literal against Peyer's published 1.18 — neither of which this produces.
    """
    calibration_yan = _beer_acid_params()["peyer_control_wort_yan_gpl"].value
    compiled = compile_scenario(
        Scenario(
            name="d238-control",
            medium="beer",
            initial={**TYRELL_SCENARIO, "yan_mgl": calibration_yan * 1000.0},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=1.0,
        )
    )
    resolved = compiled.parameters.resolve()
    shipped = compiled.parameters["peptide_buffer_capacity_beer"].value
    target = compiled.parameters["wort_buffering_capacity_peyer"].value

    solved = acidbase.peptide_capacity_for_wort_bc(compiled.y0, compiled.schema, resolved, target)
    assert solved == shipped, (
        f"the runtime back-solve returns {solved!r} against the shipped {shipped!r}. It is meant "
        "to RE-DERIVE the offline root, which is what makes the nominal member exact rather than "
        "close — do not paper over this with a tolerance"
    )

    # The YAN dependence itself, pinned — the price of the qualifier above. A wort that is not
    # the calibration wort roots elsewhere, and the size of "elsewhere" is what says whether the
    # coordinate matters. It is also the positive control on the line above: an implementation
    # that ignored `N` in the titration would return the shipped literal for BOTH worts and this
    # assert is the only thing that would notice
    # [[feedback-a-non-vacuity-check-can-itself-be-vacuous]].
    tyrell = compile_scenario(
        Scenario(
            name="d239-yan-dependence",
            medium="beer",
            initial=dict(TYRELL_SCENARIO),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=1.0,
        )
    )
    at_tyrell = acidbase.peptide_capacity_for_wort_bc(
        tyrell.y0, tyrell.schema, tyrell.parameters.resolve(), target
    )
    assert at_tyrell == pytest.approx(1.4300292930172551, rel=1e-9), (
        f"Tyrell's 200 mg/L wort roots at {at_tyrell!r}; D-239 measures 1.4300292930172551, "
        "0.35 % below the calibration wort's. If this has collapsed onto the shipped literal, "
        "the amino-buffer term has stopped reaching the titration"
    )
    assert at_tyrell < shipped, (
        "more wort nitrogen means more free amino-acid buffering, so LESS peptide is needed to "
        "reach Peyer's 1.18 — a root above the shipped literal has the sign backwards"
    )

    # ...and the two ways of building the sample agree, so the file's titration and the
    # production path are titrating one wort rather than two that happen to look alike.
    params = _beer_acid_params()
    from_state = acidbase._totals_molar(compiled.y0, compiled.schema, params.resolve())
    for slot, name in _WORT_BC_SEEDS.items():
        expected = params[name].value / acidbase.ALL_ACIDS[slot].molar_mass
        assert from_state[slot] == pytest.approx(expected, rel=1e-12), (
            f"the compiled wort's {slot} is not its `{name}` seed; the state-built and "
            "parameter-built samples have diverged and the control above is comparing two worts"
        )


def test_every_sampled_member_carries_peyers_wort_bc():
    """The repair (decision D-238) — every member's wort reproduces Peyer's 1.18, not just one.

    D-214 found the pair incoherent off-nominal and did not fix it; D-233 re-measured it (0.0100
    pH at day 14 on the low-pKa arm, 21 % of that arm's whole defect) and left it, because the
    fix looked like it would cost the round-trip guard its teeth. It does not — see this
    section's header.

    **Read the SEEDED SLOT, never the parameter.** ``peptide_buffer_capacity_beer`` still holds
    the nominal literal after the repair, by design, so a guard that read the parameter would
    pass identically before and after and pin nothing at all.

    The tolerance is the root-finder's, not a physical allowance: every member is rooted to
    ``xtol=1e-15`` on a target of exactly 1.18, so anything above 1e-9 is a wort that was never
    re-solved rather than one solved imprecisely.
    """
    compiled, ens = _beer_ensemble()
    peptide = ens.schema.slice("peptide_buffer")
    target = compiled.parameters["wort_buffering_capacity_peyer"].value
    molar_mass = acidbase.ALL_ACIDS["peptide_buffer"].molar_mass

    seeds = [float(ens.members[i][peptide, 0][0]) for i in range(ens.n_succeeded)]
    assert len(set(seeds)) == len(seeds), (
        "every member seeded the same capacity — either the draw is degenerate or rule 3 never "
        "fired, and the per-member BC below would then be vacuously the nominal's"
    )

    for i in range(ens.n_succeeded):
        values = ens.member_params[i]
        totals = acidbase._totals_molar(ens.members[i][:, 0], ens.schema, values)
        totals["peptide_buffer"] = seeds[i] / molar_mass
        bc = acidbase.peyer_fast_bc(totals, acidbase.build_pka_map(values))
        assert bc == pytest.approx(target, abs=1e-9), (
            f"member {i} drew pKa {values['pKa_peptide_buffer']:.4f} and carries a wort at "
            f"BC = {bc:.6f}, not Peyer's {target}. Its capacity was fitted to somebody else's "
            "pKa — the pair is incoherent again (D-214, D-233 §1)"
        )


def test_the_nominal_member_keeps_the_compiled_capacity_exactly():
    """The exact-nominal skip, which is what keeps D-24's byte-for-byte nominal claim structural.

    Rule 3 returns without touching the array when the member's ``pKa_peptide_buffer`` is the
    nominal one. The root-find *does* reproduce the shipped literal bit for bit (asserted
    directly above), so this branch changes no number today — it exists so that the nominal run's
    byte-for-byte reproducibility does not rest on a root-finder's tolerance surviving a SciPy
    upgrade. A guard clause, deliberately, and never a tolerance on the equality.
    """
    compiled, ens = _beer_ensemble()
    peptide = compiled.schema.slice("peptide_buffer")
    shipped = compiled.parameters["peptide_buffer_capacity_beer"].value

    builder = compiled.y0_for_member()
    assert builder is not None, "beer with an `initial_ph` must get a per-member builder"
    nominal = builder(compiled.parameters.resolve())

    assert float(nominal[peptide][0]) == shipped, "the nominal member's capacity moved"
    assert np.array_equal(nominal, compiled.y0), (
        "the whole rebuilt array must equal the compiled one at the nominal draw — this covers "
        "every present and future rule without being edited for each"
    )
    # ...and the ensemble's own nominal run is the compiled y0 too, not merely the builder's.
    assert float(ens.nominal[peptide, 0][0]) == shipped


def test_a_scenario_that_names_its_peptide_buffer_is_not_re_capacitated():
    """The other half of rule 3's branch — D-24's exclusion, kept intact.

    ``peptide_buffer_gpl`` states this wort's buffering protein. That is a scenario INPUT, and a
    parameter draw may never overwrite one; it is also the seam
    :func:`test_no_process_touches_the_peptide_buffer_pool` names as the RIGHT place to model
    less buffering protein, so a rule that quietly re-solved it would break that route.

    Guarded rather than trusted, because the branch is invisible from the repaired side: every
    assertion in :func:`test_every_sampled_member_carries_peyers_wort_bc` would still pass if the
    rule fired unconditionally (D-236's Arm F).

    The consequence is stated as well as the mechanism: a stated pool is NOT re-solved, so its
    members really do carry worts away from Peyer's 1.18 — and that is correct, because the
    scenario asked for a wort Peyer did not measure.
    """
    compiled, ens = _beer_ensemble({"peptide_buffer_gpl": TYRELL_NAMED_PEPTIDE_GPL})
    peptide = ens.schema.slice("peptide_buffer")
    target = compiled.parameters["wort_buffering_capacity_peyer"].value
    molar_mass = acidbase.ALL_ACIDS["peptide_buffer"].molar_mass

    seeds = [float(ens.members[i][peptide, 0][0]) for i in range(ens.n_succeeded)]
    assert seeds == [TYRELL_NAMED_PEPTIDE_GPL] * len(seeds), (
        "the per-member builder overwrote `peptide_buffer_gpl` with a back-solved capacity — "
        "D-24's exclusion (scenario inputs are never sampled) is breached"
    )
    drawn = [float(ens.member_params[i]["pKa_peptide_buffer"]) for i in range(ens.n_succeeded)]
    assert len(set(drawn)) == len(drawn), "the draw produced no spread; the arm is vacuous"

    # The positive control the "not re-solved" claim owes: at a stated pool the members' worts
    # genuinely do leave Peyer's target, which is what makes the assertion above a real branch
    # rather than a wort that happens to sit at 1.18 anyway.
    off = []
    for i in range(ens.n_succeeded):
        totals = acidbase._totals_molar(ens.members[i][:, 0], ens.schema, ens.member_params[i])
        totals["peptide_buffer"] = seeds[i] / molar_mass
        off.append(
            abs(
                acidbase.peyer_fast_bc(totals, acidbase.build_pka_map(ens.member_params[i]))
                - target
            )
        )
    assert max(off) > 1e-3, (
        "a stated peptide pool still lands on Peyer's BC for every member, so this arm cannot "
        "tell a skipped rule from a fired one"
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

    **Since D-238 something else writes this slot, and it is not a Process.**
    ``CompiledScenario.y0_for_member``'s rule 3 re-solves the capacity per ensemble member and
    seeds it BEFORE the cation anchor reads it — which is the pre-anchor seam this docstring
    already names as the right place, reached from the parameter side rather than the scenario
    key. It is not a counter-example to anything here: it moves the pool where the charge balance
    still sees it, so no member is left with a cation balancing a pool that is gone. The
    assertion below is unchanged and still means what it says — no *Process* drains this pool
    mid-ferment.
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

    **D-239 built a buffer-removal term and this test is what priced its shape in advance.** The
    three free amino-acid side chains are a real, sourced buffer that really does leave, and they
    land exactly where this test says such a term lands: 0.0038 pH at day 1 against 0.0202 at day
    7, a ratio of 5.3. So the term shipped on fidelity while agreeing WORSE with the one day beer
    misses [[feedback-closer-to-reality-decides]] — this test was never an argument that no such
    term exists, only that none of them answers D-211 sec 9's brief.
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
    assert early == pytest.approx(0.008688, abs=0.0005), (
        f"a 20 % pre-pitch protein loss moves day 1 by {early:.6f} pH; D-214 measured 0.017252, "
        "D-222 0.008954 at Tyrell's own counted pitch, D-223 0.009810 at the re-anchored uptake "
        "rate, and D-239 0.008688 once 7.3 % of the wort's buffering moved out of the peptide "
        "pool and into three amino-acid side chains that leave with the yeast — a 20 % loss of "
        "a SMALLER pool is a smaller loss"
    )
    assert late == pytest.approx(0.058608, abs=0.0015), (
        f"a 20 % pre-pitch protein loss moves day 7 by {late:.6f} pH; D-214 measured 0.058699, "
        "D-222 0.057136, D-223 0.059003 and D-239 0.058608"
    )
    assert late > 3.0 * early, (
        f"buffer removal is supposed to be LATE-weighted (D-214 measured day 7 at 3.4x day 1, "
        f"D-222 at 6.4x, D-223 at 6.0x, D-239 at 6.7x); here day 7 is {late / early:.2f}x day 1. "
        "If this "
        "ratio has fallen below 3, the shape "
        "argument that refuses trub settling as an answer to D-211 sec 9's brief no longer holds "
        "and the refusal needs re-measuring, not re-asserting."
    )


def test_the_trub_window_is_empty_at_the_edge_that_parks_it(tmp_path, beer_params):
    """The refusal as arithmetic, on the arm the parking question actually lives on.

    D-210 parked trub settling on the HIGH ``nitrogen_uptake_charge_beer`` edge and D-211 sec 9
    re-priced it there, so the nominal cannot answer it [[feedback-pin-the-band-not-the-nominal]].
    At that edge D-214 measured: day 1 needs a loss of **>= 27.6 %** to come inside its ceiling,
    while day 7 can afford **<= 3.1 %** before falling through its floor — empty by about
    ninefold, so no protein-loss fraction satisfies both ends.

    **D-222 RE-MEASURED both ends at Tyrell's own counted pitch and the window is emptier, not
    fuller.** Day 7 now affords **<= 12.6 %** (the course sits less acidic, so there is more
    headroom above the floor), but day 1 is now **unreachable at any loss whatsoever**: removing
    the ENTIRE peptide buffer leaves day 1 at 5.4448 against a 5.4010 ceiling. The refusal
    therefore stops being an arithmetic near-miss and becomes a saturation — which is a stronger
    form of the same verdict, and the reason this test asserts the 100 % arm rather than the old
    5 % one. Both halves are asserted so that a later beat cannot read "day 7 affords more now"
    as the window having opened.

    The old margins were tight and D-214 said so; D-222 loosened them and **D-223 gave them
    straight back**: the control's headroom above the floor went +0.0086 (D-214) to +0.0358
    (D-222) and +0.0033 at D-223. A faster engine ferments more completely and finishes more
    acidic, so the high nitrogen-charge edge went back onto the envelope's floor.

    **D-239 SPENT the last of it, and this test's premise moved with it.** The three free
    amino-acid side chains cost the day-7 course 0.020 pH, so the high edge's baseline is now
    **0.0176 BELOW** the floor rather than 0.0033 above it. The affordable further loss at that
    edge is therefore not small, it does not exist — which is why the quantitative end of this
    test moved to the NOMINAL edge, where a bracket still exists (1.72 %), and the high edge is
    kept as the pinned statement of what D-239 spent. Reading a bracket off an arm whose
    baseline is already outside would be scoring a margin against a deficit
    [[feedback-a-margin-is-a-claim-about-what-holds-it-open]].

    **The verdict is unchanged and is now stronger at both ends.** Day 1 stays a SATURATION at
    every edge — removing the entire peptide buffer still leaves it above its ceiling — and day
    7 affords 9.98 % / 1.72 % / nothing across the low, nominal and high edges. No protein-loss
    fraction satisfies both ends anywhere in the band.
    [[feedback-read-a-fast-curve-on-a-fixed-grid]] still applies
    to how this helper reads the curve (``np.interp`` onto the exact hour, never ``argmin`` over
    the solver's own output).
    """
    param = beer_params["nitrogen_uptake_charge_beer"]
    hi = param.uncertainty.high
    data_dir = _beer_data_dir_with_nitrogen_charge(tmp_path, hi)
    nominal_dir = _beer_data_dir_with_nitrogen_charge(tmp_path, param.value)
    floor_7 = 4.804 - 0.024
    ceiling_1 = TYRELL_PH_COURSE[1][1] + TYRELL_PH_READ_TOL

    high_baseline = _tyrell_ph_with_peptide_loss(data_dir, 0.0, 7.0)
    assert high_baseline == pytest.approx(4.7624, abs=0.005), (
        f"the high edge's day-7 baseline is {high_baseline:.4f}; D-239 measures 4.7624, which is "
        f"{floor_7 - high_baseline:.4f} pH BELOW the {floor_7:.3f} floor. Its history is "
        "+0.0086 headroom at D-214, +0.0358 at D-222, +0.0033 at D-223, and D-239's amino-acid "
        "split spent the rest. If this is back above the floor, something has returned "
        "buffering to the finished beer and the arithmetic below is scored against the wrong "
        "baseline"
    )
    assert high_baseline < floor_7, (
        "the high edge is expected OUTSIDE the day-7 envelope since D-239 — this is the cost "
        "that beat priced, not a regression to repair by weakening the term"
    )
    nominal_baseline = _tyrell_ph_with_peptide_loss(nominal_dir, 0.0, 7.0)
    assert nominal_baseline > floor_7, (
        f"the NOMINAL edge's day-7 baseline is {nominal_baseline:.4f}, at or below the "
        f"{floor_7:.3f} floor. D-239 measures 4.7846. The nominal is where this test's "
        "bracket now lives, so if it too has gone outside there is no arm left to score a "
        "window on and the refusal has to be re-derived rather than cited"
    )
    # The end that CLOSES the window: day 1 cannot be reached at all, so no loss fraction can
    # satisfy both ends however much room day 7 has gained.
    day1_at_total_loss = _tyrell_ph_with_peptide_loss(data_dir, 1.0, 1.0)
    assert day1_at_total_loss > ceiling_1, (
        f"removing the ENTIRE peptide buffer leaves day 1 at {day1_at_total_loss:.4f}, at or "
        f"below its {ceiling_1:.4f} ceiling. D-223 measures 5.4282 — unreachable. If a loss "
        "fraction can now reach day 1, D-214's refusal is no longer a saturation and the "
        "window has to be re-derived at both ends rather than cited"
    )
    # The end that is merely bounded: what day 7 can afford, re-bisected at D-223. The bracket
    # is deliberately wider than the bisection's own precision (8.935 %): a tight bracket here
    # would go red on solver noise rather than on the quantity moving.
    assert (
        _tyrell_ph_with_peptide_loss(nominal_dir, 0.015, 7.0)
        > floor_7
        > (_tyrell_ph_with_peptide_loss(nominal_dir, 0.020, 7.0))
    ), (
        "day 7's affordable loss at the NOMINAL edge is no longer bracketed by [1.5 %, 2.0 %]; "
        "D-239 bisects it at 1.72 %. The arm moved from the high edge to this one at D-239 "
        "because the high edge's own baseline left the envelope, so the series to compare "
        "against is this edge's: D-214 3.1 %, D-222 12.6 % and D-223 1.2 % were all read at the "
        "HIGH edge and are NOT continuous with this number "
        "[[feedback-right-number-wrong-condition]]"
    )
    # And the high edge stated as what it now is: nothing is affordable there, so the loss that
    # would be needed to satisfy day 1 is refused by day 7 before it starts. Asserted rather
    # than left implicit, because "the baseline is already out" is exactly the condition under
    # which a bisection would return a meaningless bracket instead of raising.
    assert _tyrell_ph_with_peptide_loss(data_dir, 0.001, 7.0) < floor_7, (
        "a 0.1 % protein loss at the high edge is inside the day-7 envelope, which would mean "
        "the baseline is back above the floor and the premise assert above is passing vacuously"
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


def test_the_nitrogen_charge_bands_high_edge_now_finishes_BELOW_tyrells_envelope(
    tmp_path, beer_params
):
    """BOTH band edges, one threshold each — not a shared floor the tight edge rides for free.

    ``test_the_model_reaches_tyrells_measured_beer_ph`` scores the NOMINAL value only, and a
    nominal that passes says nothing about a band whose edges the sampler actually draws
    [[feedback-pin-the-band-not-the-nominal]]. This runs the full model at the low, nominal and
    high edges of ``nitrogen_uptake_charge_beer`` and pins each day-7 pH separately.

    **D-222 SPENT the finding this test was built to hold visible, and that is recorded rather
    than deleted.** At the retired 1.0 g/L scenario pitch the high edge landed at 4.783 against
    an admissible floor of 4.780 — a margin of 0.003 pH, i.e. the band straddled the point of
    going past the measurement. At Tyrell's own counted pitch the engine attenuates less by day
    7, so the whole course sits less acidic and the high edge now clears the floor by **0.036**,
    twelvefold more room. The caution it was written for is therefore no longer live, and a beat
    adding further acidification should re-measure this edge rather than cite either number: the
    margin is a property of the ferment's extent, which is exactly what the pitch correction
    moved. The term itself is still derived from published wort composition, never fitted, and
    still a LOWER bound (the buffer-removal half of nitrogen uptake is inexpressible here and
    pushes the same way).

    **D-239 BUILT that buffer-removal half and this test's headline claim is what it cost.** The
    wort's three free amino-acid side chains — the ones Peyer names and this model had never
    carried — buffer at t=0 and leave with the yeast, worth 0.020 pH at day 7. The low and
    nominal edges stay inside Tyrell's envelope; **the HIGH edge does not**, finishing 0.018
    below the floor. The test is renamed for what it now forbids rather than what it used to
    assert [[feedback-name-guards-for-what-they-forbid]], and it keeps its teeth in the form
    that matters: each edge is still pinned separately, the ordering is still asserted, and the
    high edge's miss is bounded ABOVE so a further slide is caught.

    **This is not an xfail and must not be converted into one.** D-208's strict-xfail idiom is
    for something TRUE of the source and FALSE of the model. Here the model became more faithful
    and the AGREEMENT got worse, which is a statement that a different term is missing — most
    likely on the alkaline side, and D-232's open growth-extent residue is the standing
    candidate. Filing a correct beat as a defect is how a later reader reverts it.
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
    for label, value in (("low", lo_ph), ("nominal", nom_ph)):
        assert window[0] <= value <= window[1], (
            f"the {label} edge finishes day 7 at {value:.4f}, outside Tyrell's envelope "
            f"{lo_bound:.3f}-{hi_bound:.3f} widened by the {TYRELL_PH_READ_TOL} read tolerance. "
            "Since D-239 only the HIGH edge is expected outside; if a second edge has followed "
            "it, the term's cost has grown beyond what that beat priced"
        )
    assert hi_ph < window[0], (
        f"the high edge finishes day 7 at {hi_ph:.4f}, INSIDE Tyrell's envelope. D-239 measures "
        f"4.7625, {window[0] - hi_ph:.4f} below the {window[0]:.3f} floor. A return to the "
        "inside means the amino-acid buffering has stopped reaching the day-7 course — check "
        "that `_totals_molar` still carries the three side chains before re-pinning this"
    )
    # Each edge pinned on its own, so a shift that moved them together could not hide inside a
    # single containment check.
    assert lo_ph == pytest.approx(4.8069, abs=0.01)
    assert nom_ph == pytest.approx(4.7846, abs=0.01)
    assert hi_ph == pytest.approx(4.7625, abs=0.01)
    # ...and the high edge's margin, asserted as the number it is rather than described. It was
    # 0.003 pH at the retired scenario pitch (D-209), 0.036 at Tyrell's counted one (D-222)
    # because a less complete ferment finishes less acidic, and D-223 gave the whole of that back
    # -- 0.0033 -- because the re-anchored uptake rate finishes the ferment again. The edge is
    # still INSIDE Tyrell's envelope; it is inside by three thousandths of a pH unit, which is
    # what this assert now says out loud. Bounded ABOVE as well, so a further loss of
    # acidification is still caught.
    assert 0.010 < window[0] - hi_ph < 0.030, (
        f"the high edge now sits {window[0] - hi_ph:.4f} pH BELOW the bottom of the admissible "
        "window. Its whole history is a margin ABOVE that floor until this beat: 0.003 at the "
        "retired 1.0 g/L scenario pitch (D-209), 0.036 at Tyrell's counted 0.398 (D-222), 0.0033 "
        "once D-223 re-anchored the uptake rate, and D-239's amino-acid split spent it and 0.018 "
        "more. Bounded on BOTH sides on purpose: a smaller miss means the term is not reaching "
        "this edge, a larger one means something beyond D-239 is acidifying the day-7 course"
    )
    # The nominal's remaining headroom, pinned as the quantity a next beat has to work inside.
    # It is what stops "the high edge is out" being read as "the band is out".
    assert 0.001 < nom_ph - window[0] < 0.015, (
        f"the nominal edge's headroom above the floor is {nom_ph - window[0]:.4f} pH; D-239 "
        "measures 0.0046. This is the room any further same-sign term has at the NOMINAL, and "
        "it is now smaller than the read tolerance the envelope was widened by"
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
    assert day1 - hi == pytest.approx(0.172, abs=0.02), (
        f"the day-1 miss is {day1 - hi:.4f} pH above the envelope; D-208 measured 0.186 above, "
        "D-209 0.315 BELOW, D-211 0.070 above, and D-222 0.172 above once the scenario carried "
        "Tyrell's own counted pitch. All of them matter: the size and the SIDE are what say "
        "whether a beat moved the timing or the arithmetic. D-222 is a SIZE cost paid for a "
        "sourced inoculum — the side is unchanged, and the count is what sets the pitch, not "
        "this number [[feedback-fit-the-observable-not-the-consequence]]"
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
    assert fraction == pytest.approx(0.298, abs=0.03), (
        f"{fraction:.3f} of the nitrogen is drawn by 24 h; D-209 measured >0.99, D-211's "
        "re-derived rate put it at 0.363 and D-222's refit at Tyrell's own counted pitch at "
        "0.298. That agreement is what makes the timing MEASURED rather than fitted"
    )
    lo_spread, hi_spread = TYRELL_N_DRAWN_SPREAD
    assert lo_spread <= fraction <= hi_spread, (
        f"the drawn fraction {fraction:.3f} is outside Tyrell's measured {TYRELL_N_DRAWN_SPREAD} "
        "cell-count spread. This is the assert that makes the pin above a MEASUREMENT: D-219 "
        "priced the pitch correction as taking this to 0.145 and outside, and it is the "
        "`mu_max` refit at the corrected pitch that keeps it inside. If it has left the spread, "
        "the refit and the pitch have come apart"
    )


# ---------------------------------------------------------------------------------------
# D-216 — the two anchors on beer's fermentation SPEED, and what they forbid
# ---------------------------------------------------------------------------------------

#: The ``q_sugar_max`` that reproduces Tyrell's day-2 extract fraction exactly, found by
#: bisection on the shipped model on the hourly grid.
#:
#: **RE-BISECTED at D-222, and it LEFT THE BAND.** D-216 measured 1.397 and built the first
#: tier of its refusal on that value being INSIDE the parameter's printed 0.3-1.5 band —
#: *"the value is out of band is NOT the reason"* was the whole point of it. That was measured
#: against a scenario pitching 2.51x the biomass Tyrell counted. At Tyrell's own pitch the
#: day-2 shortfall widens 2.81x -> 4.21x and the knob that closes it is **2.3226**, which is
#: **1.55x the printed high edge**; the retired 1.397 now books 0.343 of the wort by day 2
#: against the measured 0.594.
#:
#: D-216 §3's framing is therefore spent, and the refusal it supported got SIMPLER rather
#: than weaker: the one knob that closes Tyrell's extract schedule is out of band on its own,
#: before the other anchor is consulted at all.
Q_SUGAR_MAX_MATCHING_TYRELL = 2.3226

#: The band ``q_sugar_max`` is drawn from, as printed in ``beer_generic.yaml``. Named here
#: because :data:`Q_SUGAR_MAX_MATCHING_TYRELL` having left it is a claim that must be RUN
#: against the shipped file rather than written into a comment that cannot go stale loudly.
#: RE-DERIVED at D-223: the retired (0.3, 1.5) spanned the translation uncertainty of a
#: rate-law transfer that no longer sources the value at all. What replaces it is the biological
#: spread across Foster's three Beer 1 ale controls crossed with two readings of the extrapolated
#: figure tail -- 4.0x wide down to 1.29x, and DRAWN, so the narrowing is a real change to every
#: beer ensemble rather than a documentation edit.
Q_SUGAR_MAX_PRINTED_BAND = (0.634, 0.818)

#: ``K_repression`` removed entirely — not a candidate value, the unbounded LIMIT of the
#: term that owns most of the lag (79 % at the retired pitch, 46 % at Tyrell's counted one —
#: D-216 §5, re-scored at D-222). Used to make the refusal two-tier.
K_REPRESSION_REMOVED = 1.0e6


def _beer_days_to_target_gravity(
    q_sugar_max: float | None = None, *, celsius: float | None = None, **overrides: float
) -> float:
    """Days for §2.2's 1.048 ale wort to reach 1.010 apparent, optionally re-rated.

    Deliberately reuses the benchmark's own wort, temperature and gravity construction rather
    than restating them, so this test cannot drift away from the anchor it is about.

    ``celsius`` defaults to the criterion's own temperature — **15 °C since D-221
    re-temperatured it**, not the 20 °C every claim in sections 10-13 was originally measured
    at. Pass 20.0 only to read a RETIRED frame deliberately, and say which record's frame it is.

    The 40 d span is D-221's: at 15 °C the shipped engine needs 9.00 d and a slow arm needs more,
    where the retired 20 °C frame fitted inside 14 d.
    """
    from tests.benchmarks.test_milestone1 import (
        _BEER_FERMENTABLE_S0,
        _BEER_OG_SG,
        TARGET_FG_SG,
        _apparent_gravity_series,
        _beer_scenario,
    )

    kwargs = {} if celsius is None else {"celsius": celsius}
    compiled = compile_scenario(_beer_scenario(duration_days=40.0, **kwargs))
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


def test_matching_tyrells_extract_schedule_overshoots_the_attenuation_benchmark():
    """Why D-215's extract xfail cannot be closed on the uptake rate (D-216 §4, D-221 §5).

    The two anchors on beer's fermentation speed are still not compatible under this rate law,
    but **the shape of the incompatibility inverted at D-221** and the old form of this test
    could not have seen it.

    * **Tyrell's measured extract course** — their wort is 59.4 % fermented by day 2, where the
      model books 14.1 % (``test_the_model_ferments_tyrells_wort_on_tyrells_schedule``). It was
      21.2 % until D-222 put the scenario on Tyrell's own counted pitch;
    * **§2.2's acceptance criterion** — a 1.048 ale wort reaching 1.010 apparent in 5-7 days, at
      **15 °C** since D-221 re-temperatured it (``test_beer_1048_og_attenuates_in_5_to_7_days``).

    Beer's uptake is ``q_sugar_max · X · Monod(S)``, so one constant scales both.

    **What D-216 measured, and the frame it was measuring in.** At the criterion's old 20 °C the
    shipped rate PASSED at 6.05 d and the admissible band was ``q`` in [0.425, 0.621]: the
    criterion sat just above the shipped 0.5 and **forbade** a faster engine, breaking at
    q ≈ 0.62 having closed under a fifth of Tyrell's gap. At the corrected 15 °C the shipped 0.5
    falls **below** the admissible band and the criterion now **demands** a faster engine —
    [0.667, 1.017] as D-221 bisected it, [0.612, 0.891] since D-222 refit ``mu_max`` at Tyrell's
    counted pitch. The band travels with the growth rate, so it is re-derived where it is used
    and never carried here as a constant.

    **The anchors agree on direction, and D-222 WIDENED their disagreement about magnitude.**
    Both want the engine faster than it ships. D-221 read the gap as ~1.37x — the criterion
    topping out at 1.017 against Tyrell's 1.397 — and called D-216 §11's open question a
    magnitude rather than a sign. At the corrected pitch the criterion tops out at 0.891 and
    Tyrell demands **2.3226**, a factor of **~2.61**. The narrowing D-221 recorded was itself
    borrowed from the pitch excess.

    **Tyrell's demand has left its own band, and that is a NEW reason rather than a louder one.**
    D-216 §3's first tier was explicitly *"the value is INSIDE 0.3-1.5, so out-of-band is NOT
    the reason"*. At Tyrell's counted pitch 2.3226 is **1.55x the printed high edge**, so the
    knob that closes his extract schedule is inadmissible before the other anchor is consulted
    at all. It is asserted below against the shipped parameter file rather than stated here, so
    that a band edit can turn it red.

    **The refusal survives the loss of its own premise.** D-216 refused 1.397 partly because the
    criterion forbade it. The criterion no longer forbids faster rates in general — but it
    still forbids the Tyrell-matching one, and by more than before: 2.3226 gives 2.38 d here and
    stays outside across the ENTIRE printed ``E_a_uptake`` band (30,000 -> 2.13 d, 63,000 ->
    2.50 d). Nothing in D-221 or D-222 licenses moving the rate.

    **The baseline is asserted first and deliberately**, and it is no longer "the criterion
    passes" — it cannot be, because D-221 established that it does not
    [[feedback-pair-the-red-with-an-ordering-preserving-baseline]]. What replaces it is the
    bracket: the shipped rate misses SLOW and the Tyrell rate misses FAST, so the override is
    shown to CROSS the window rather than merely to sit outside it on the same side as its own
    baseline, which would attribute nothing.
    """
    shipped = _beer_days_to_target_gravity()
    assert 5.0 < shipped < 7.0, (
        f"the benchmark wort attenuates in {shipped:.2f} d at the shipped q_sugar_max. D-223 "
        "re-anchored the rate to Foster's measured 15 °C course and measured 6.04 d, INSIDE "
        "§2.2's 5-7 d window; D-221 measured 9.00 d and D-222 8.50 d, both outside on the slow "
        "side. If this is outside again the engine has been re-rated and D-223's adjudication "
        "of beer's two speed anchors needs redoing"
    )
    # ...AT BOTH BAND EDGES, not only at the nominal. `q_sugar_max` is DRAWN, and this file's own
    # joint-corner test states the reason: a claim verified at a point where the sampler reads a
    # band is the archive's most repeated defect shape. D-223's record says the whole band clears
    # the window; that is a claim, and this is the guard for it
    # [[feedback-grep-finds-claims-not-guards]]. The LOW edge is the tight one -- 6.79 d against a
    # 7.0 ceiling, 0.21 d of room -- so a band widened downward turns this red before a draw does.
    edges = {q: _beer_days_to_target_gravity(q_sugar_max=q) for q in Q_SUGAR_MAX_PRINTED_BAND}
    assert edges[Q_SUGAR_MAX_PRINTED_BAND[0]] == pytest.approx(6.79, abs=0.05)
    assert edges[Q_SUGAR_MAX_PRINTED_BAND[1]] == pytest.approx(5.42, abs=0.05)
    for q, days in edges.items():
        assert 5.0 < days < 7.0, (
            f"q_sugar_max = {q} (a printed band EDGE) attenuates the benchmark wort in "
            f"{days:.2f} d, outside §2.2's 5-7 d window. D-223 measured the whole band inside "
            "(6.79 -> 5.42 d), which is what makes the criterion a property of the parameter "
            "rather than of the nominal the sampler happens to centre on"
        )

    matched = _beer_days_to_target_gravity(Q_SUGAR_MAX_MATCHING_TYRELL)
    assert matched < 5.0, (
        f"re-rated to the q_sugar_max that reproduces Tyrell's day-2 extract "
        f"({Q_SUGAR_MAX_MATCHING_TYRELL}), the benchmark wort attenuates in {matched:.2f} d, "
        "which is INSIDE §2.2's window. D-222 measured 2.38 d at 15 °C (D-221 measured 3.99 d "
        "at the retired pitch's 1.397, D-216 2.71 d at the retired 20 °C). If this is now "
        "inside, the two anchors no longer conflict and D-215's extract xfail is closable on "
        "this knob — which is a result, not a test failure: re-open D-216 §4"
    )

    # D-222's own tier, run against the SHIPPED file so that a band edit turns it red. D-216
    # §3 rested on this value being in band; at Tyrell's counted pitch it is not, and that is
    # a reason of a different KIND from the one above.
    # [[feedback-grep-finds-claims-not-guards]]
    printed = load_parameters(default_data_dir() / "beer_generic.yaml")["q_sugar_max"]
    assert (printed.uncertainty.low, printed.uncertainty.high) == Q_SUGAR_MAX_PRINTED_BAND, (
        f"q_sugar_max's printed band is now ({printed.uncertainty.low}, "
        f"{printed.uncertainty.high}); D-216 §3 and D-222 §4 both reason against "
        f"{Q_SUGAR_MAX_PRINTED_BAND}, so re-read them before re-pinning this"
    )
    assert printed.uncertainty.high < Q_SUGAR_MAX_MATCHING_TYRELL, (
        f"the Tyrell-matching q_sugar_max {Q_SUGAR_MAX_MATCHING_TYRELL} is back INSIDE its "
        f"printed band (high edge {printed.uncertainty.high}). D-216 §3's *'out of band is NOT "
        "the reason'* would then be live again and D-222 §4's simpler refusal would be gone"
    )
    ratio = Q_SUGAR_MAX_MATCHING_TYRELL / printed.uncertainty.high
    assert ratio == pytest.approx(2.839, abs=0.03), (
        f"the Tyrell-matching rate is {ratio:.3f}x the printed high edge; D-222 measured 1.548x "
        "against the retired 0.3-1.5 band and D-223 2.839x against the re-derived 0.634-0.818. "
        "The losing anchor did not move; the band it is judged against narrowed to what one "
        "measured trial actually spans"
    )

    # THE SEPARATION, which is what makes this an adjudication rather than a preference. Until
    # D-223 both rates missed the window and the bracket was "one misses slow, the other misses
    # fast". Since D-223 the shipped rate is INSIDE it and the Tyrell-matching one is outside on
    # the fast side, which is strictly stronger: the criterion now DISTINGUISHES the two anchors
    # instead of merely ordering them
    # [[feedback-pair-the-red-with-an-ordering-preserving-baseline]].
    assert matched < 5.0 < shipped < 7.0, (
        f"the criterion must SEPARATE the two anchors: shipped {shipped:.2f} d inside the "
        f"5-7 d window, re-rated {matched:.2f} d outside it on the fast side. If they now sit "
        "on the same side of the window, D-223's adjudication has lost the thing that decided it"
    )


def test_removing_catabolite_repression_entirely_still_misses_tyrells_schedule():
    """The refusal's second tier: not even the unbounded limit of the dominant term (D-216 §5).

    ``K_repression`` = 2.0 g/L is tier ``speculative``, source *"author estimate"* — the
    functional form is Gee & Ramirez's but ``beer_generic.yaml`` records that their numeric
    constants were not accessible in-source. With glucose at 12.3 g/L it holds maltose at 14 %
    of its rate on day 0 and maltotriose at 0.5 %, so the model's day 1 is essentially
    glucose-only. It is by far the largest single contributor to the early-limb lag, and it has
    the right SHAPE (a brake that vanishes as glucose clears, matching a lag that peaks at day 2
    and closes by day 7).

    **What removing it buys, in both frames.** At the retired 1.0 g/L pitch it took day 2 from
    0.212 to 0.514, **79 % of the gap** (D-216 §5). At Tyrell's own counted pitch it takes day 2
    from 0.141 to 0.350, **46 %** — the term did not shrink, the gap grew, because D-222 widened
    the day-2 shortfall 2.81x -> 4.21x.

    So the obvious next move is to re-source that constant. **This test says it would not be
    enough**, and says it more loudly than D-216 could. Removed ENTIRELY — not a candidate
    value, the limit — the model still falls short of Tyrell's day 2, while putting the
    benchmark at 4.54 d (3.42 d when D-216 measured it at the retired pitch), i.e. still on the
    fast side of §2.2's window. The refusal in D-216 §6 is therefore not "no in-band point
    works" but the stronger "not even the unbounded limit of the term that owns most of it".
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
    assert closed == pytest.approx(0.68, abs=0.05), (
        f"removing repression closes {closed:.0%} of the day-2 gap; D-216 measured 79 % at the "
        "retired scenario pitch, D-222 46 % at Tyrell's counted one, and D-223 68 % once the "
        "uptake rate was re-anchored to Foster's course. The share moves with the size of the "
        "gap it is a share OF, and that gap has now closed from 4.21x to 3.16x short: the term "
        "is still the largest single contributor and still cannot reach the measurement on its "
        "own, which is what the two-tier refusal in D-216 §5 rests on"
    )


def test_the_ph_course_endorses_NEITHER_pitch_and_the_count_is_what_decides():
    """The pH course cannot be used to pick beer's pitch, in either direction (D-222 §6).

    **What this test asserted before, and why it is spent.** D-216 §8 found beer's published pH
    agreement conditional on ``TYRELL_SCENARIO``'s 1.0 g/L — a figure nothing sourced — and read
    it generously: *"1.0 g/L is the value at which the model's biomass reproduces Tyrell's
    measured growth timing, and two independent observables endorse it against the per-cell
    arithmetic"*. D-219 then showed the 1.0 was a RESIDUAL rather than a choice, and D-222 moved
    the scenario onto Tyrell's own counted 9.96e6 cells/mL and refit ``mu_max`` there.

    **At the refit rate that endorsement is FALSE, and this test is the measurement.** The day-1
    miss does not merely shrink toward the retired pitch — it CROSSES ZERO between them:

    ===================  ===============  ================
    pitch                day-1 vs 5.377   days inside (/8)
    ===================  ===============  ================
    0.398 (counted)      +0.172 alkaline  7
    0.5                  +0.072 alkaline  7
    1.0 (retired)        -0.259 ACIDIC    7
    ===================  ===============  ================

    So the pH course has an optimum near ~0.55 g/L and neither the counted pitch nor the retired
    one sits on it; the retired pitch is not "the value two observables endorse", it is a value
    that overshoots on the other side. All three score 7 of 8 days, so the COUNT is what carries
    the decision and the score cannot referee it
    [[feedback-fit-the-observable-not-the-consequence]].

    **What this forbids** is the obvious future move: "restore the pitch, the pH course scores
    better". It does not, and a pitch chosen to zero a downstream residual would be booking every
    unbuilt beer-acid term into an inoculum. D-209 §8's buffer-removal half was the term the
    +0.172 was said to leave room for; **D-239 built it, and it took only 0.008 of the 0.172**,
    because a buffer that leaves with the yeast is late-weighted by construction (that shape is
    pinned by ``test_losing_wort_protein_acidifies_late_not_early``). So the day-1 miss is still
    open and it no longer has a named candidate.
    """

    def score(pitch: float) -> tuple[float, int]:
        compiled = compile_scenario(
            Scenario(
                name="d222-pitch",
                medium="beer",
                initial={**TYRELL_SCENARIO, "pitch_gpl": pitch},
                temperature_schedule=[TemperaturePoint(day=0.0, celsius=TYRELL_TRIAL_CELSIUS)],
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

        inside = sum(
            band_lo - TYRELL_PH_READ_TOL <= ph_at(day) <= band_hi + TYRELL_PH_READ_TOL
            for day, (band_lo, band_hi) in TYRELL_PH_COURSE.items()
        )
        return ph_at(1) - TYRELL_PH_COURSE[1][1], inside

    counted_miss, counted_inside = score(TYRELL_SCENARIO["pitch_gpl"])
    retired_miss, retired_inside = score(TYRELL_SCENARIO_RETIRED_PITCH_GPL)

    assert counted_miss == pytest.approx(0.164, abs=0.03) and counted_miss > 0.0, (
        f"at Tyrell's counted pitch day 1 misses by {counted_miss:+.4f} pH; D-222 measured "
        "+0.172 and D-239 +0.164, on the alkaline side. The SIDE is what D-222 said leaves room "
        "for D-209 §8's buffer-removal half — D-239 BUILT that half and it bought 0.008 pH of "
        "the 0.172, because the term is late-weighted by construction. A sign change here is a "
        "different finding, not a drift"
    )
    assert retired_miss == pytest.approx(-0.278, abs=0.03) and retired_miss < 0.0, (
        f"at the retired 1.0 g/L pitch day 1 misses by {retired_miss:+.4f} pH; D-222 measured "
        "-0.259, on the ACIDIC side. If the retired pitch is back to being the better-scoring "
        "one, D-216 §8's endorsement reading is live again and D-222 §6 needs re-measuring"
    )
    assert counted_miss > 0.0 > retired_miss, (
        "the day-1 miss no longer CROSSES ZERO between the two pitches. That crossing is the "
        "whole claim: it is what makes the pH course unable to referee the pitch"
    )
    # D-239 SEPARATED them, which this assert's own message asked for as a restatement rather
    # than a re-pin. The counted pitch keeps 7 of 8 days; the retired one drops to 5, losing days
    # 6 and 7 through the FLOOR — a heavier pitch ferments further and finishes more acidic, and
    # the amino-acid buffering that leaves with the yeast takes another 0.020 pH off the late
    # course on top of that.
    #
    # **The restatement is narrow and the claim it protects is unchanged.** The pH course now
    # agrees with the count instead of being neutral between the two pitches, so nothing about
    # "restore the pitch, the pH course scores better" has become true — it has become MORE
    # false. What is no longer available is the stronger form of D-222 §6's argument, that the
    # course cannot referee the pitch AT ALL: on the coarse measure it now can, and it referees
    # in favour of the counted inoculum. The reason that must not be read as an endorsement is
    # the one D-222 gave and this test still asserts above: the day-1 miss CROSSES ZERO between
    # them, so whichever pitch scores better on a day count is still choosing between two
    # different WRONG signs [[feedback-a-summary-statistic-is-not-the-curve]].
    assert counted_inside == 7, (
        f"the counted pitch keeps {counted_inside} of 8 days inside; D-222 and D-239 both "
        "measure 7, with day 1 the only miss"
    )
    assert retired_inside == 5, (
        f"the retired 1.0 g/L pitch keeps {retired_inside} of 8 days; D-222 measured 7 and D-239 "
        "measures 5, losing days 6 and 7 below the floor. If it is back to 7 the separation has "
        "closed and D-222 §6's stronger claim is live again — restate it, do not re-pin this"
    )
    assert counted_inside > retired_inside, (
        "the sourced inoculum must not score WORSE than the retired one on the coarse measure; "
        "if it does, the count-vs-score argument is being asked to defend a losing arm"
    )


# ======================================================================================
# 11. The temperature sensitivity of uptake — the lever D-216 named, and what holds it
#     open (decision D-217)
# ======================================================================================


def test_the_uptake_activation_energy_is_no_longer_inert_at_the_attenuation_benchmark():
    """D-216 §6's decoupling lever is GONE, not weakened (D-217 §1, corrected at D-221 §6).

    **What this test used to assert, and why it was true.** D-216 §6 named `E_a_uptake` as the
    only lever that decouples beer's two speed anchors — the one knob that could move Tyrell's
    15 °C course without moving §2.2 — on the grounds that the benchmark ran at exactly
    ``T_ref`` = 20 °C, so its Arrhenius factor was 1.0 whatever ``E_a`` was. D-217 §1 asserted
    that across the whole printed band and well outside it, on both signs, and measured it at
    **exactly 0.0000 d**. That measurement was correct and is preserved below as a statement
    about the retired frame.

    **Why it stops being true.** The freedom was never a fact about yeast. It was an artefact of
    the two anchors sitting at different temperatures, and §2.2's temperature was the one the
    literature contradicts. D-221 re-temperatured the criterion to 15 °C — which is
    :data:`TYRELL_TRIAL_CELSIUS` **exactly**. Both anchors now run at one temperature, so
    ``E_a_uptake`` enters both through an identical Arrhenius factor and can no longer move
    either one alone.

    **There is no magnitude argument to fall back on.** Across the printed 30,000-63,000 J/mol
    band the criterion moves **1.75 d** — against a window only 2.0 d wide. The parameter went
    from the criterion's most inert to one of its strongest levers in a single frame change.
    D-216 §6's argument does not survive in a weaker form; it is spent.

    The span is pinned on the criterion's own **1 h output grid**, whose quantum is 0.0417 d. On
    a 4x grid D-221 measured 1.7708 d; the 0.0208 d difference is half that quantum and not a
    model difference, so the tolerance below is set wider than one quantum deliberately
    [[feedback-pin-tolerance-vs-solver-tolerance]].

    **What is NOT affected.** D-217's refusal to re-source ``E_a_uptake`` rests on the corpus
    having nothing, not on this lever, and is untouched. And the low band edge does not rescue
    the criterion either: 30,000 J/mol gives 7.69 d, still outside 5-7.

    **A RED names which half moved.** The 20 °C arm is a claim about arithmetic (at ``T_ref`` the
    Arrhenius factor is 1.0 by construction) and can only fail if uptake stops being Arrhenius.
    The 15 °C arm is a claim about the shipped frame and fails if the criterion is re-temperatured
    again.
    """
    from tests.benchmarks.test_milestone1 import _BEER_BENCH_CELSIUS

    assert _BEER_BENCH_CELSIUS == TYRELL_TRIAL_CELSIUS, (
        f"the criterion runs at {_BEER_BENCH_CELSIUS} °C and Tyrell's trial at "
        f"{TYRELL_TRIAL_CELSIUS} °C. D-221's whole finding is that these coincide, which is what "
        "makes E_a_uptake degenerate between the two anchors. If they have separated again, the "
        "decoupling lever may be back and D-216 §6 needs re-reading rather than citing"
    )

    # The RETIRED 20 C frame: D-217's measurement, preserved as a claim about arithmetic.
    retired = _beer_days_to_target_gravity(celsius=20.0)
    for e_a in (-97000.0, 0.0, 30000.0, 55100.0, 63000.0, 80000.0):
        got = _beer_days_to_target_gravity(celsius=20.0, E_a_uptake=e_a)
        assert got == retired, (
            f"at the RETIRED 20 °C frame E_a_uptake = {e_a:.0f} J/mol moves the criterion to "
            f"{got:.4f} d against {retired:.4f} d. This must stay EXACTLY zero: 20 °C is T_ref, "
            "where the Arrhenius factor is 1.0 by construction. A RED means uptake has stopped "
            "being Arrhenius, not that the frame moved"
        )

    # The LIVE 15 C criterion: the lever D-216 s6 relied on being absent.
    lo = _beer_days_to_target_gravity(E_a_uptake=30000.0)
    hi = _beer_days_to_target_gravity(E_a_uptake=63000.0)
    span = hi - lo
    assert span == pytest.approx(1.2083, abs=0.05), (
        f"across the printed E_a_uptake band the live criterion moves {span:.4f} d; D-221 "
        f"measured 1.75 on this 1 h grid and D-223 1.2083 ({lo:.4f} -> {hi:.4f}) after "
        "re-anchoring the uptake rate -- a quicker engine spends less of the run in the "
        "temperature-sensitive limb, so the lever shrinks without the frame moving. TWO "
        "different changes land here and the endpoints above tell them apart. Collapsed toward "
        "ZERO with both endpoints near each other: the criterion has drifted back to T_ref and "
        "D-221 §6 needs redoing. Merely SMALLER with both endpoints faster: the uptake rate "
        "moved, and a quicker engine spends less of the run in the temperature-sensitive limb — "
        "check q_sugar_max before concluding anything about the frame"
    )
    window_width = BENCHMARKS["beer_attenuation"].high - BENCHMARKS["beer_attenuation"].low
    assert span > 0.5 * window_width, (
        f"the E_a_uptake band moves the criterion {span:.4f} d against a window {window_width:.1f} "
        "d wide. D-216 §6's decoupling argument required this to be ZERO; D-221 records that it "
        "is most of the window. A RED here would mean the lever is small enough to argue on "
        "magnitude again, which is a result — re-open D-216 §6"
    )
    # This assert used to read `lo > BENCHMARKS["beer_attenuation"].high` -- i.e. that even the
    # FASTEST edge of the E_a_uptake band left the criterion outside the window, so E_a could not
    # be a back door to closing D-221's miss without touching the rate. D-223 closed that miss on
    # the rate itself, and the back door closed with it: the whole printed band now sits inside
    # the window (5.1667 -> 6.3750 against 5-7). The claim to guard is therefore the opposite one,
    # and it is worth more -- the criterion's verdict on the shipped rate does not depend on which
    # E_a_uptake the sampler draws.
    lo_w, hi_w = BENCHMARKS["beer_attenuation"].low, BENCHMARKS["beer_attenuation"].high
    assert lo_w < lo <= hi <= hi_w, (
        f"the printed E_a_uptake band takes the criterion to [{lo:.4f}, {hi:.4f}] d against a "
        f"[{lo_w}, {hi_w}] d window. D-223 measured the whole band inside. If an edge has left "
        "it, the criterion has become E_a-conditional and §2.2's pass is no longer a property of "
        "the rate alone"
    )

    # E_a_growth's note carried the same 20 C inertness claim and loses it the same way, but
    # it is a much weaker lever and the numbers are asserted so the two notes cannot drift apart.
    g_retired = {
        e: _beer_days_to_target_gravity(celsius=20.0, E_a_growth=e) for e in (30000.0, 63000.0)
    }
    assert g_retired[30000.0] == g_retired[63000.0], (
        f"at the RETIRED 20 °C frame E_a_growth spans {g_retired}. Like E_a_uptake this must be "
        "exactly zero at T_ref; a RED means growth has stopped being Arrhenius"
    )
    g_span = _beer_days_to_target_gravity(E_a_growth=63000.0) - _beer_days_to_target_gravity(
        E_a_growth=30000.0
    )
    assert g_span == pytest.approx(0.1667, abs=0.05), (
        f"across its printed band E_a_growth moves the live criterion {g_span:.4f} d; D-221 "
        "measured 0.2917 and D-222 0.1667 at the refit growth rate. beer_generic.yaml's note "
        "for that parameter quotes this number"
    )
    assert g_span < 0.25 * span, (
        f"E_a_growth moves the criterion {g_span:.4f} d against E_a_uptake's {span:.4f} d. The "
        "file header's claim that growth is the weaker lever on attenuation is what makes the "
        "two corrected notes consistent; if they have converged, both notes need re-reading"
    )


def test_the_uptake_activation_energy_is_a_lever_only_because_the_trial_ran_cool():
    """The lever's whole size is Tyrell's distance from ``T_ref`` (D-217 §4, D-221 §6).

    **The sentence this docstring used to open with is now false and is kept here as the
    correction, not deleted.** It read: *"`E_a_uptake` is free at the benchmark, so it is the
    one parameter that can move Tyrell without moving §2.2."* That was true only while the two
    anchors sat at different temperatures, and the difference was an artefact — §2.2 was
    asserted at a 20 °C the literature contradicts. D-221 re-temperatured it to 15 °C, which is
    :data:`TYRELL_TRIAL_CELSIUS` exactly, so both anchors now take the SAME Arrhenius factor and
    the parameter moves both together. At the live criterion the printed band is worth 1.75 d
    against a 2.0 d window (see
    ``test_the_uptake_activation_energy_is_no_longer_inert_at_the_attenuation_benchmark``).

    What survives, and is what this test actually measures, is the arithmetic underneath: the
    lever exists at all only because a trial runs away from ``T_ref``. Had Tyrell's tube run AT
    ``T_ref``, sweeping the entire printed band would move nothing. That is still true, still
    asserted below on both arms, and is now a statement about Arrhenius rather than about a
    decoupling that no longer exists.

    Both arms are asserted because one alone would mislead. The 15 °C arm on its own reads
    as "there is a lever"; the 20 °C arm on its own reads as "there is no lever". Together
    they say what is true: the lever is worth exactly what the frame is worth, and the frame
    rests on §3.2's comparative sentence (see :data:`TYRELL_TRIAL_CELSIUS`), not on a printed
    fermentation temperature.

    The span is small either way. Across the whole 30,000-63,000 band the day-2 fraction moves
    0.027 at Tyrell's own counted pitch (0.045 at the retired 1.0 g/L one, which is what D-217
    measured) — under a tenth of the gap — which is why D-216 refused the low edge on
    provenance rather than on power.
    """

    def day2_span(celsius: float) -> tuple[float, float]:
        lo = _tyrell_flux_fraction(celsius=celsius, E_a_uptake=30000.0)[2]
        hi = _tyrell_flux_fraction(celsius=celsius, E_a_uptake=63000.0)[2]
        return lo, hi

    cool_lo, cool_hi = day2_span(TYRELL_TRIAL_CELSIUS)
    span_cool = cool_lo - cool_hi
    assert span_cool == pytest.approx(0.0382, abs=0.005), (
        f"at {TYRELL_TRIAL_CELSIUS:.0f} °C the printed E_a_uptake band moves Tyrell's day-2 "
        f"fraction by {span_cool:.4f}; D-217 measured 0.0449 at the retired scenario pitch, "
        "D-222 0.0267 at Tyrell's counted one, and D-223 0.0382 once the uptake rate was "
        "re-anchored — the lever scales with the flux it acts on, and it moves DOWN with the "
        "biomass and UP with the rate. If it has collapsed to zero the lever D-216 §6 named is "
        "gone and its refusal needs re-reading"
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
FOSTER_PITCH_CELLS_PER_ML = BEER_COUNTED_PITCH_CELLS_PER_ML  # one copy, in tests/conftest.py

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
    'retired: the unsourced upper end of D-216 §7\'s "textbook 40-60 pg"': 60.0,
    "retired: back-computed from the beer scenario pitch": 100.0,
}

#: ``q_sugar_max`` that lands the model exactly on Foster's 72 h endpoint at 22 °C, per
#: reading of the per-cell dry mass. Bisected on the shipped model on the hourly grid, and
#: verified live by the test below rather than trusted.
#:
#: **RE-DERIVED at D-222** at the refit ``mu_max``: faster growth means more catalytic biomass
#: earlier, so every row needs LESS uptake rate to reach the same endpoint. The retired values,
#: measured at ``mu_max`` = 0.034, were {18: unreachable, 40: 0.9242, 50: 0.8176, 100: 0.5602}.
#:
#: **18 pg/cell is now IN here, and that is a retired claim rather than a new finding.** At the
#: old growth rate the band's ceiling of 1.5 could not reach Foster's endpoint at that reading
#: at all (3.04 d against a published <=3, a saturation); at the refit rate it reaches it at
#: 1.0592 and the ceiling arrives in 2.38 d. The reading itself is still retired — D-219 settled
#: the per-cell mass at 40 pg — and its §2.2 row is still outside the window, so nothing about
#: the verdict moves. What moves is that every row is now a crossing.
Q_SUGAR_MAX_REACHING_FOSTER: dict[float, float] = {
    18.0: 1.0592,
    40.0: 0.7886,
    50.0: 0.7225,
    100.0: 0.5314,
}

#: What the band's CEILING reaches at 18 pg/cell, in days against a published <=3. At the
#: retired growth rate this was 3.04 — a saturation, the reason that reading had no crossing.
#: At the refit rate it is 2.38 and the row is an ordinary crossing above.
FOSTER_DAYS_AT_BAND_CEILING_18PG = 2.3833

#: The 2 d arm on the 100 pg row: the ``q_sugar_max`` reaching Foster's target at the OPEN
#: end of its sampling interval (:data:`FOSTER_SAMPLE_HOURS`) rather than at 72 h. It is the
#: whole difference between "one corner of the bracket survives" and "none does", so it is
#: named here rather than left as a literal in the one test that reads it.
Q_SUGAR_MAX_REACHING_FOSTER_AT_2D_100PG = 0.8228

#: What §2.2's benchmark then reads, per reading. Only one lands inside 5-7 d.
#: What §2.2's criterion reads when the engine is re-rated to Foster's 72 h endpoint, per
#: reading of the per-cell dry mass. **RE-MEASURED at D-221 against the criterion's corrected
#: 15 °C**; the retired 20 °C figures D-218 measured are kept below because the two together
#: are the finding — the surviving row moves from the RETIRED 100 pg reading to the SETTLED
#: 40 pg one, so the temperature repair and the cell-mass settlement agree where the old
#: pairing had them in conflict.
HANDOFF_DAYS_AT_FOSTER_RATE: dict[float, float] = {
    18.0: 4.3333,
    40.0: 5.5833,
    50.0: 6.0417,
    100.0: 8.0,
}

#: D-218's own figures at the criterion's retired 20 °C frame. Pinned, not deleted: a claim
#: measured in a frame that has moved is still a true claim about that frame, and the inversion
#: between these two dicts is what D-221 §7 asserts.
HANDOFF_DAYS_AT_FOSTER_RATE_RETIRED_20C: dict[float, float] = {
    18.0: 2.9167,
    40.0: 3.75,
    50.0: 4.0833,
    100.0: 5.375,
}

#: The lowest ``E_a_uptake`` at which the SETTLED 40 pg row sits inside §2.2's window (D-221 §7,
#: re-bisected at D-222). Below it the Foster-matching rate overshoots the criterion on the fast
#: side, so the compatibility is real but conditional. D-221 measured 40,165 J/mol — 69 % of the
#: printed 30,000-63,000 band — and at the refit ``mu_max`` it is 36,623, i.e. **80 %** of the
#: band, with the shipped 55,100 well inside. The conditionality LOOSENED; it did not go away,
#: and the claim is still *"compatible over most of the band including the shipped value"*.
E_A_UPTAKE_ADMITTING_FOSTER_AT_SETTLED_MASS = 36623.0


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

#: A finer grid for the ONE claim whose margin is smaller than the quantum above. D-222's refit
#: moved the Arrhenius crossing toward the band's high edge, leaving a margin of 0.0027 against
#: a 6-minute quantum of ~0.005 — unresolved. At 1 minute the quantum is ~0.0005 and the reading
#: reproduces to five decimals against a 15-second one, so the value is the model's rather than
#: the grid's. Used only by the crossing test; everything else stays commensurable on the coarse
#: grid [[feedback-read-a-fast-curve-on-a-fixed-grid]].
FOSTER_FINE_GRID_PER_HOUR = 60


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
    # The engine-state control. It used to read "§2.2 passes"; D-221 established that it does
    # not, at the criterion's corrected 15 °C, so the control now asserts the state D-221
    # measured. Its job is unchanged — to distinguish "the arms below say something" from "the
    # engine is in an unexpected state" — but a control asserting something false attributes
    # nothing at all.
    baseline = _beer_days_to_target_gravity()
    assert baseline == pytest.approx(6.04, abs=0.1), (
        f"the criterion wort attenuates in {baseline:.2f} d at the shipped parameters; D-221 "
        "measured 9.00 d at the corrected 15 °C, D-222 8.50 d after refitting `mu_max` at "
        "Tyrell's counted pitch, and D-223 6.04 d after re-anchoring `q_sugar_max` to Foster's "
        "measured course. The control failed, so nothing below is attributable"
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

    **D-222 moved the crossing and this test had to change its GRID, not its claim.** Refitting
    ``mu_max`` at Tyrell's counted pitch pushed the crossing from ~58 to ~61.6 kJ/mol, so the
    high-edge margin fell from 0.013 to **0.0027** — under the module fixture's 6-minute quantum
    of ~0.005, where the same reading on a 1-minute grid resolves it at a quantum of ~0.0005 and
    reproduces the value to five decimals (0.99733 against 0.99746). The edges are therefore
    integrated HERE at :data:`FOSTER_FINE_GRID_PER_HOUR` rather than read off the shared fixture.
    A margin shrinking toward the band edge is what a rate change does to a crossing; a margin
    read below its own quantum is what a stale grid does to a test.

    Asserting the two ENDS is what makes the middle mean something. Without them,
    ``|residual - 1| < 0.02`` reads as "the model is a bare Arrhenius response"; with them it
    reads "the model is not, and this is where the two curves cross".
    """
    gas_constant = 8.314462618
    t_cold, t_hot = 273.15 + 12.0, 273.15 + 22.0
    light_pitch = _foster_pitch_gpl(18.0)

    def residual(e_a: float, per_hour: int = FOSTER_FINE_GRID_PER_HOUR) -> float:
        arrhenius = float(np.exp((e_a / gas_constant) * (1.0 / t_cold - 1.0 / t_hot)))
        hot = _foster_days_to_target_gravity(
            light_pitch, 22.0, days=12.0, per_hour=per_hour, E_a_uptake=e_a
        )
        cold = _foster_days_to_target_gravity(
            light_pitch, 12.0, days=25.0, per_hour=per_hour, E_a_uptake=e_a
        )
        return (cold / hot) / arrhenius

    low, high = residual(30000.0), residual(63000.0)
    # The nominal keeps the shared fixture: its own claim carries a 0.02 tolerance, hundreds of
    # times the coarse quantum, so paying for two more fine integrations would buy nothing.
    nominal = float(foster_temperature_sweep[(18.0, 55100.0)]) / float(
        np.exp((55100.0 / gas_constant) * (1.0 / t_cold - 1.0 / t_hot))
    )

    assert low > 1.0 and high < 1.0, (
        f"the residual is {low:.4f} at the band's low edge and {high:.4f} at its high edge. "
        "D-218 §3 measured 1.108 and 0.987 and D-222 1.0686 and 0.9973 — a sweep that CROSSES "
        "1.0 inside the band. If both are now on one side there is no crossing, and the "
        "near-identity at the nominal below would be a structural property of the rate law "
        "rather than a coincidence"
    )
    assert min(low - 1.0, 1.0 - high) > 0.002, (
        f"the crossing's margins are {low - 1.0:.4f} and {1.0 - high:.4f}; D-218 measured 0.108 "
        f"and 0.013, D-222 0.0686 and 0.0027. The quantum on this ratio at "
        f"{FOSTER_FINE_GRID_PER_HOUR}/h is ~0.0005, so a margin under 0.002 is not resolved by "
        "a factor of four and the assert above would be reading the grid rather than the model. "
        "The high-edge margin is the one that shrank: at 0.0027 it is ~6x its own quantum, and "
        "a further rate increase would push the crossing past the band's high edge entirely"
    )
    assert abs(nominal - 1.0) < 0.02, (
        f"at the shipped E_a_uptake the residual is {nominal:.4f}; D-218 measured 1.0119 and "
        "D-222 1.0109. The nominal sitting within 1 % of the bare Arrhenius factor is what "
        "makes the model's temperature ratio LOOK like an identity, and §3's finding is that "
        "it is a crossing"
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


def test_fosters_endpoint_and_the_handoff_window_agree_at_the_settled_conversion():
    """The two literature anchors stop conflicting once the criterion is re-temperatured
    (D-218 §3, INVERTED at D-221 §7).

    D-216 §11 asked which anchor to calibrate beer's speed against. D-218 answered *"Foster,
    and the brief is wrong"* — because at the criterion's then-asserted 20 °C a rate reaching
    Foster's 72 h endpoint took §2.2 to 3.63 d at the settled cell mass, well outside 5-7 d.
    D-221 re-temperatured the criterion to 15 °C, and that reading inverts.

    Each row: the ``q_sugar_max`` that puts the model on Foster's 72 h endpoint at one reading
    of the per-cell dry mass, and what §2.2's criterion then reads in each frame.

    ======================  =========  ===============  ==============  ==============
    per-cell mass           pitch g/L  ``q_sugar_max``  §2.2 at 15 °C   §2.2 at 20 °C
    ======================  =========  ===============  ==============  ==============
    18 pg (retired)             0.216  1.5 = CEILING          3.83 d          2.58 d
    40 pg (SETTLED, D-219)      0.480  0.924            **5.42 d IN**         3.63 d
    50 pg                       0.600  0.818            **5.96 d IN**         4.00 d
    100 pg (retired)            1.200  0.560                  8.17 d    5.50 d IN
    ======================  =========  ===============  ==============  ==============

    **The survivor moves from the retired reading to the settled one.** At 20 °C exactly one
    row sat inside the window and it was the 100 pg conversion D-219 retired — the survivor was
    the reading known to be wrong. At 15 °C the surviving rows are the settled 40 pg one and its
    neighbour, and the retired 100 pg row is the one that now fails. Two independent repairs —
    D-219's cell-mass settlement and D-221's temperature correction — agree where the old
    pairing had them in conflict.

    **The agreement is real but CONDITIONAL, and the sweep is why that is known.** D-221 measured
    that the printed ``E_a_uptake`` band moves this criterion 1.75 d against a 2.0 d window, so
    a single in-band reading proves nothing on its own
    [[feedback-a-hit-can-be-two-errors-cancelling]]. Swept, the settled row sits inside the
    window from :data:`E_A_UPTAKE_ADMITTING_FOSTER_AT_SETTLED_MASS` upward — 69 % of the printed
    band, with the shipped 55,100 inside it — and overshoots on the fast side below that
    (30,000 -> 4.71 d). The claim is *"compatible over most of the band including the shipped
    value"*, never *"compatible"*.

    **D-218 §3's central caveat also inverts.** There the survivor held only because Foster's
    3 d was read as exact; at the open end of the sampling interval the 100 pg row went to
    3.71 d and the window broke in all eight bracket cells. At 15 °C that same open-end arm
    reads 5.54 d — inside. The window's survival is no longer conditional on treating a 72 h
    sample as an exact duration.

    **What this test does NOT say is that the model passes.** All four rows are RE-RATED engines.
    The shipped ``q_sugar_max`` = 0.5 takes the criterion to 9.00 d and that is D-221's headline
    xfail. What inverts here is which *literature* anchors can be satisfied together, not whether
    the engine satisfies them. Tyrell's extract schedule remains incompatible with both
    (``test_matching_tyrells_extract_schedule_overshoots_the_attenuation_benchmark``).

    A RED names one of the pins moving; check first whether the criterion's temperature moved.
    """
    from tests.benchmarks.test_milestone1 import TARGET_FG_SG

    assert TARGET_FG_SG == 1.010, (
        f"the benchmark's target gravity is {TARGET_FG_SG}; Foster's is 1.01 (they cite Parker "
        "2008). The comparison in this test is only sound while the two agree"
    )

    # D-218's saturated row is GONE and the retirement is asserted rather than deleted: at the
    # retired growth rate the band's ceiling of 1.5 could not reach Foster's endpoint at 18 pg
    # at all (3.0417 d against <=3), which is why that row had no crossing. At D-222's refit
    # rate the ceiling arrives in 2.38 d and the row bisects like every other. The reading is
    # still retired and its §2.2 verdict is still OUT, so no conclusion moves.
    at_ceiling = _foster_days_to_target_gravity(
        _foster_pitch_gpl(18.0), 22.0, days=12.0, q_sugar_max=1.5
    )
    assert at_ceiling == pytest.approx(FOSTER_DAYS_AT_BAND_CEILING_18PG, abs=0.05), (
        f"at 18 pg/cell the band CEILING q = 1.5 reaches Foster's endpoint in {at_ceiling:.4f} d; "
        f"D-222 measured {FOSTER_DAYS_AT_BAND_CEILING_18PG} where D-218 measured 3.0417 at the "
        "retired growth rate"
    )
    assert at_ceiling < FOSTER_DAYS_TO_FG_AT_22C, (
        f"the band ceiling at 18 pg is back to MISSING Foster's endpoint ({at_ceiling:.4f} d). "
        "That row would then be a saturation again and it must come out of "
        "Q_SUGAR_MAX_REACHING_FOSTER, which may hold crossings only"
    )
    verdicts = {}

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

    assert verdicts == {18.0: False, 40.0: True, 50.0: True, 100.0: False}, (
        f"the window's survival per per-cell reading is now {verdicts}. D-221 measured survivors "
        "at the SETTLED 40 pg reading and its 50 pg neighbour, with the retired 100 pg reading "
        "failing — the exact inverse of D-218's {18: F, 40: F, 50: F, 100: T} at the retired "
        "20 °C. A change here is a change to the ANSWER of D-216 §11, not to a number"
    )

    # The agreement above is CONDITIONAL on E_a_uptake, and unswept it would be a coincidence
    # wearing a result's clothes. [[feedback-a-hit-can-be-two-errors-cancelling]]
    settled_q = Q_SUGAR_MAX_REACHING_FOSTER[40.0]
    below = _beer_days_to_target_gravity(
        q_sugar_max=settled_q, E_a_uptake=E_A_UPTAKE_ADMITTING_FOSTER_AT_SETTLED_MASS - 5000.0
    )
    above = _beer_days_to_target_gravity(
        q_sugar_max=settled_q, E_a_uptake=E_A_UPTAKE_ADMITTING_FOSTER_AT_SETTLED_MASS + 5000.0
    )
    assert below < 5.0 <= above, (
        f"the settled row reads {below:.4f} d below the crossing and {above:.4f} d above it; "
        f"D-221 bisected the crossing at {E_A_UPTAKE_ADMITTING_FOSTER_AT_SETTLED_MASS:.0f} J/mol. "
        "If it no longer straddles, the compatibility is either unconditional or gone, and "
        "either way the docstring's '69 % of the band' is wrong"
    )
    assert _beer_days_to_target_gravity(q_sugar_max=settled_q, E_a_uptake=30000.0) < 5.0, (
        "the settled row is inside the window at the LOW printed E_a_uptake edge too, so the "
        "compatibility is unconditional and D-221 §7's careful hedge should be retired"
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
    assert 5.0 <= open_end <= 7.0 and open_end == pytest.approx(5.375, abs=0.05), (
        f"reading Foster's endpoint at the open end of its sampling interval takes §2.2 to "
        f"{open_end:.4f} d; D-221 measured 5.54 and D-222 5.375 at 15 °C, INSIDE 5-7, where "
        "D-218 measured 3.71 "
        "at the retired 20 °C and outside. This is the assert that says D-218 §3's central "
        "caveat inverted: the window's survival no longer depends on treating a 72 h SAMPLE as "
        "an exact duration"
    )
    retired_open_end = _beer_days_to_target_gravity(
        q_sugar_max=Q_SUGAR_MAX_REACHING_FOSTER_AT_2D_100PG, celsius=20.0
    )
    # This is a LIVE reading of the CURRENT model in D-218's retired frame, not a replay of
    # D-218's own number, and D-222 is where that distinction started to matter: the figure was
    # 3.71 d when the growth rate was 0.034 and is 3.625 d at the refit 0.058. What has to
    # survive is the SIDE — the open-end arm outside 5-7 at 20 C and inside at 15 C — because
    # that is what makes the inversion a frame change rather than an arithmetic one.
    assert retired_open_end == pytest.approx(3.625, abs=0.05), (
        f"the open-end arm re-reads as {retired_open_end:.4f} d in D-218's retired 20 °C frame; "
        "D-218 measured 3.71 at the growth rate of its day and D-222 3.625 at the refit one"
    )
    assert not (5.0 <= retired_open_end <= 7.0), (
        f"the open-end arm now lands INSIDE 5-7 d at the retired 20 °C too ({retired_open_end:.4f}"
        " d). The inversion above would then not be about the frame at all, and D-221 §7 would "
        "need re-deriving rather than citing"
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
TYRELL_COUNTED_PITCH_GPL = cells_per_ml_to_pitch_gpl(TYRELL_PITCH_CELLS_PER_ML)

#: What the scenario shipped from D-178 to D-222, over what Tyrell counted — D-219's headline
#: for the beer side, and the gap D-222 closed. Kept because the RETIRED pitch is still read as
#: an arm by three tests and the excess is what made it wrong.
TYRELL_SCENARIO_BIOMASS_EXCESS = 2.51

#: At Tyrell's counted pitch, which is now the SHIPPED pitch: the day-2 and day-7 extract
#: fractions and the nitrogen drawn by 24 h, on the hourly grid §12's flux helper uses.
#: **Measured at D-222, with ``mu_max`` refit at that pitch.** D-219 measured the same three
#: without the refit — 0.0903, 0.7010 and 0.1446 — which is what its price list quotes.
TYRELL_AT_COUNTED_PITCH = {
    "day2_fraction": 0.18773,
    "day7_fraction": 0.95867,
    "n_drawn_24h": 0.29743,
}

#: Tyrell's measured 24 h cell-count spread — the interval D-211 calls "what makes the
#: timing MEASURED rather than fitted", and which the counted pitch falls outside.
TYRELL_N_DRAWN_SPREAD = (0.234, 0.448)

#: The shipped model against Foster's endpoint at the SETTLED pitch, no knob touched.
#: D-218 §4's cleanest line was "3.33 d against a published <=3, an 11 % miss" — true, and
#: read on the 100 pg reading D-219 retires. These are the numbers that replace it.
FOSTER_AT_SETTLED_PITCH = {"d22": 3.2417, "d12": 7.1792}


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


def test_the_beer_scenario_now_carries_tyrells_own_counted_pitch():
    """``TYRELL_SCENARIO`` pitches what Tyrell counted — corrected at D-222, measured at D-219.

    **What this test asserted before.** D-219 measured the scenario carrying **2.51x** the
    biomass Tyrell pitched and declined to correct it, because ``mu_max`` had been fitted at
    1.0 g/L and moving the pitch inherited a refit. This test was where that refusal lived.
    D-222 took the refit and the correction together, so the excess is now 1.000 by construction
    and what needs guarding is the CONVERSION and the PRICE, not the gap.

    **Why the correction and not the measurement.** The 1.0 g/L was never chosen by anyone: it
    is what you get by dividing back from a pitch nobody sourced, so it absorbed the true cell
    mass *and* every error in the model's per-gram rate — a residual, which is why it implied
    ~100 pg/cell against a settled 40 (D-219 §1). Tyrell states a COUNT, 9.96e6 cells/mL, and
    the engine's biomass gram is defined by a count (Coleman's 4e-11 g). So the pitch is now a
    conversion of the paper's own number rather than a free parameter, and every Tyrell
    comparison in this file is scored at the trial's own inoculum.

    **The price, asserted here so it cannot be quietly forgotten.** D-219's price list was
    measured WITHOUT the refit and overstated the cost: it put the day-2 extract shortfall at
    6.58x and the nitrogen drawn by 24 h at 0.145, outside Tyrell's spread. With ``mu_max``
    refit at the corrected pitch by D-211's own method the shortfall is **4.21x** and the
    nitrogen fraction is back **inside** the spread. What does not come back is attenuation:
    day 7 reaches **0.782** of the fermentable extract against a measured 0.997, where the
    excess pitch reached 0.931 — so the engine is now visibly slow against Tyrell's ENDPOINT
    and not only its early limb, which is what D-220 measured independently against Foster.

    A RED here is the pitch or the conversion moving, and the two have very different causes.
    """
    assert TYRELL_SCENARIO["pitch_gpl"] == pytest.approx(TYRELL_COUNTED_PITCH_GPL, rel=1e-12), (
        f"the beer scenario pitches {TYRELL_SCENARIO['pitch_gpl']:.4f} g/L against Tyrell's "
        f"counted {TYRELL_COUNTED_PITCH_GPL:.4f}. D-222 corrected these to be the same number "
        "through the boundary conversion; a difference means one of them was edited to a literal"
    )
    assert (
        pytest.approx(TYRELL_SCENARIO_BIOMASS_EXCESS, abs=0.01)
        == TYRELL_SCENARIO_RETIRED_PITCH_GPL / TYRELL_COUNTED_PITCH_GPL
    ), (
        "the RETIRED 1.0 g/L pitch no longer works out to D-219's measured 2.51x excess over "
        "the counted pitch, so the settled conversion has moved and D-219 §1 needs re-reading"
    )

    counted, n_drawn = _tyrell_at_pitch(TYRELL_COUNTED_PITCH_GPL)
    retired, n_retired = _tyrell_at_pitch(TYRELL_SCENARIO_RETIRED_PITCH_GPL)

    assert counted[2] == pytest.approx(TYRELL_AT_COUNTED_PITCH["day2_fraction"], abs=0.005), (
        f"at Tyrell's counted pitch the model ferments {counted[2]:.4f} of the wort by day 2; "
        f"D-222 measured {TYRELL_AT_COUNTED_PITCH['day2_fraction']} at the refit growth rate "
        "(D-219 measured 0.0903 before it)"
    )
    shortfall = TYRELL_FLUX_FRACTION[2] / counted[2]
    assert shortfall == pytest.approx(3.16, abs=0.15), (
        f"the day-2 shortfall at the counted pitch is {shortfall:.2f}x; D-219 measured 6.58, "
        "D-222 4.21 with the growth refit, and D-223 3.16 after re-anchoring the uptake rate to "
        "Foster's measured course. It is still far short, which is why D-215's extract xfail "
        "stays xfail, and it is reported rather than tuned"
    )
    assert shortfall > TYRELL_FLUX_FRACTION[2] / retired[2], (
        "the lag no longer gets WORSE at the counted pitch than at the retired one. That "
        "reverses D-216 §7's direction and would mean the biomass excess was masking a fast "
        "model rather than a slow one"
    )
    assert counted[7] == pytest.approx(TYRELL_AT_COUNTED_PITCH["day7_fraction"], abs=0.01), (
        f"day 7 reaches {counted[7]:.4f} of the fermentable extract against Tyrell's measured "
        f"{TYRELL_FLUX_FRACTION[7]:.3f}; D-222 measured "
        f"{TYRELL_AT_COUNTED_PITCH['day7_fraction']}. The ENDPOINT half of the slowness is the "
        "part the excess pitch was hiding (it reached 0.931), and it is the half D-220 "
        "measured independently against Foster's course"
    )

    assert n_drawn == pytest.approx(TYRELL_AT_COUNTED_PITCH["n_drawn_24h"], abs=0.01), (
        f"nitrogen drawn by 24 h at the counted pitch is {n_drawn:.4f}; D-222 measured "
        f"{TYRELL_AT_COUNTED_PITCH['n_drawn_24h']}"
    )
    assert TYRELL_N_DRAWN_SPREAD[0] <= n_drawn <= TYRELL_N_DRAWN_SPREAD[1], (
        f"nitrogen drawn by 24 h is {n_drawn:.4f}, outside Tyrell's measured "
        f"{TYRELL_N_DRAWN_SPREAD} spread. D-219 priced the pitch correction as taking this to "
        "0.145 and outside; the refit is what keeps it inside, so a RED means the two have "
        "come apart and D-211's timing attribution no longer holds at the shipped pitch"
    )


def test_the_handoff_window_survives_the_settled_conversion_once_it_is_temperatured_right():
    """§2.2's window is corroborated, not refuted — the refutation was the temperature
    (D-218 §4, D-219 §5c, INVERTED at D-221 §7).

    **This test previously asserted the opposite, and its own RED message named this outcome.**
    It read: *"That reverses D-218 §4 and the brief's window would be corroborated by a
    third-party trial rather than refuted by one."* That is precisely what happened, and it
    happened without any rate moving.

    §12's fork table sweeps four readings of the per-cell mass. At the criterion's asserted
    20 °C the window survived in exactly one cell — the 100 pg reading D-219 retired — so the
    verdict was *the window does not survive*. D-220 then recovered Foster's own course and
    found the brief's 5-7 d duration real at **15 °C**, where the same three commercial ale
    controls take 5.06-6.26 d, and impossible at 20 °C, where they take 2.91-3.77 d. D-221
    re-temperatured the criterion accordingly, and at the settled **40 pg/cell** the
    Foster-matching rate takes it to **5.42 d — inside**.

    Two separate things are asserted, and the order matters:

    1. the settled reading is IN the swept bracket and its verdict is now True — the sweep is
       not being re-run here, it is being pointed at;
    2. the readings whose verdict is False now include the RETIRED 100 pg one, which is the
       only cell that survived in the old frame.

    **What died is the PAIRING, not the duration.** D-218 §4 and D-219 §5c both concluded that
    nothing supports 5-7 days. Corrected, the duration is supported and the 20 °C it was
    asserted at is not — [[feedback-right-number-wrong-condition]] for the second time on this
    same window.

    **The window is still not something the engine passes.** These rows are re-rated engines;
    the shipped rate takes the criterion to 9.00 d, which is D-221's xfail. Retiring or
    satisfying the window means moving ``q_sugar_max``, and D-216 §4's shape objection is
    untouched by any of this. The verdict is recorded; the knob is not moved.

    A RED here is a change to the ANSWER, not to a digit.
    """
    from tests.test_units import RETIRED_READINGS_PG, SETTLED_BAND_PG

    settled = [
        pg for pg in PER_CELL_DRY_MASS_PG.values() if SETTLED_BAND_PG[0] <= pg <= SETTLED_BAND_PG[1]
    ]
    assert settled == [40.0], (
        f"the swept bracket now has {settled} inside the settled band, not exactly [40.0]. The "
        "verdict below reads ONE row of §12's table and needs that row to be unambiguous. The "
        "likeliest cause is not the bracket but the BAND: widen biomass_N_fraction enough and the "
        "60 pg reading enters it, so check tests/test_units.py first"
    )
    assert 40.0 in Q_SUGAR_MAX_REACHING_FOSTER, (
        "the settled reading is no longer one of the bisected rows, so §12's table no longer "
        "says anything about the reading that ships"
    )

    benchmark = _beer_days_to_target_gravity(q_sugar_max=Q_SUGAR_MAX_REACHING_FOSTER[40.0])
    assert benchmark == pytest.approx(HANDOFF_DAYS_AT_FOSTER_RATE[40.0], abs=0.05)
    assert 5.0 <= benchmark <= 7.0, (
        f"at the settled conversion the Foster-matching rate takes §2.2 to {benchmark:.4f} d, "
        "which is OUTSIDE the 5-7 d window. D-221 measured 5.42 d at the corrected 15 °C. If "
        "this has fallen back outside, check the criterion's temperature before anything else"
    )
    retired_frame = _beer_days_to_target_gravity(
        q_sugar_max=Q_SUGAR_MAX_REACHING_FOSTER[40.0], celsius=20.0
    )
    assert not (5.0 <= retired_frame <= 7.0) and retired_frame == pytest.approx(
        HANDOFF_DAYS_AT_FOSTER_RATE_RETIRED_20C[40.0], abs=0.05
    ), (
        f"in the RETIRED 20 °C frame the same rate reads {retired_frame:.4f} d; D-218 measured "
        f"{HANDOFF_DAYS_AT_FOSTER_RATE_RETIRED_20C[40.0]}, outside 5-7. Both halves are asserted "
        "because the finding is the DIFFERENCE between the frames: if the retired frame now "
        "agrees too, the inversion is not the temperature and D-221 §7 is misattributed"
    )
    survivors = {pg for pg, days in HANDOFF_DAYS_AT_FOSTER_RATE.items() if 5.0 <= days <= 7.0}
    assert 40.0 in survivors and not survivors & set(RETIRED_READINGS_PG.values()), (
        f"the readings at which §2.2's window survives are {sorted(survivors)}. D-221 measured "
        "{40.0, 50.0}: the SETTLED conversion must be among them and no RETIRED one may be. "
        "At the criterion's old 20 °C this assert held the opposite — the sole survivor was the "
        "retired 100 pg reading — so a RED here says the frames have swapped back and D-221 §7 "
        "needs re-reading rather than citing"
    )
    assert survivors, (
        "no reading in the bracket survives 5-7 d at all. The verdict is unchanged but the "
        "CONTRAST this test is built on is gone -- 'the settled reading is the survivor' becomes "
        "vacuous, and the guard would then pass on a broken sweep"
    )


def test_at_the_settled_conversion_the_model_meets_fosters_endpoint_out_of_sample():
    """The number that replaced D-218 §4's "11 % miss", and where it ended up (D-219 -> D-223).

    D-218 §4 closed on its cleanest line: *"at the archive's own 100 pg reading, with no knob
    touched at all, the shipped model reaches S1's endpoint in 3.33 d against a published
    <=3. An 11 % miss."* That was measured, and it was measured on the reading D-219 retires.

    At the settled 40 pg/cell — Foster's counted 1.2e7 cells/mL is **0.48 g/L**, not 1.2 —
    the shipped model took **4.84 d** against the same published <=3 (a **1.61x** miss) and also
    missed the 12 C ceiling at 10.74 d. D-222 re-measured both at the refit growth rate: 4.48 d,
    a 1.49x miss, and 9.93 d — inside the 12 C ceiling. **D-223 re-anchored ``q_sugar_max`` to
    Foster's measured 15 C course and both land: 3.24 d, a 1.08x miss, and 7.18 d.**

    **This test changed its name because its claim changed sign, and the reading needs care.**
    22 C and 12 C are not the temperature the rate was fitted at, so these are out-of-sample
    checks and the agreement is not bought. But neither number is a duration measured against a
    duration: 3 d is a CEILING (Foster's sampling interval) and 10 d is an INCUBATION LENGTH.
    Clearing a ceiling is weaker than matching a course, which is why the load-bearing comparison
    lives in §14 against the recovered course itself, and why §13's own conclusion is now stated
    as "meets" rather than "beats". The brief's 5-7 d window is still refuted at 20 C by every
    third-party endpoint — that half is untouched by the re-anchoring.

    Read against §12's other guards this is consistent, not new: D-215 measured the same
    engine fermenting Tyrell's wort ~2.8x too slowly *at a pitch already carrying 2.51x the
    counted biomass*. Those two numbers COMPOUND — they are not two routes agreeing on one
    figure — which is why the day-2 shortfall at Tyrell's counted pitch was 6.58x before the
    refit and 4.21x after it.

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
    assert miss == pytest.approx(1.08, abs=0.03), (
        f"the shipped model reaches Foster's endpoint in {d22:.4f} d against a published "
        f"<={FOSTER_DAYS_TO_FG_AT_22C}, a {miss:.2f}x miss; D-219 measured 1.61, D-222 1.49 at "
        "the refit growth rate, and D-223 1.08 after re-anchoring the uptake rate. Note what "
        'this number is measured against: a CEILING ("<= 3 d"), not a duration, and 22 °C is '
        "not the temperature the rate was fitted at — so this is an out-of-sample check, and the "
        "residual 8 % is the model still being slightly slow rather than a fit residual"
    )
    # D-222's refit CLEARED Foster's 12 C ceiling — 9.93 d against <=10 — and that was a weaker
    # fact than it sounded, so the claim is anchored on the recovered course rather than on the
    # ceiling. The ceiling is the INCUBATION LENGTH of a run whose 12 C strains had not finished
    # inside it, i.e. a bound D-220's own course superseded: measured there the mean is 7.60 d,
    # and since D-223 the model reads 7.18 d against it — 0.94x, where D-222 was 1.31x slow.
    # Agreement across temperature is asserted in §14 on those ratios, not here on a ceiling.
    assert d12 < FOSTER_CEILING_DAYS_AT_12C, (
        f"the 12 C arm finishes in {d12:.4f} d, back OUTSIDE Foster's "
        f"{FOSTER_CEILING_DAYS_AT_12C} d incubation ceiling; D-222 measured 9.93 d inside it. "
        "D-219 measured 10.74 d at the retired growth rate, so a RED here means the refit has "
        "come undone"
    )
    # This used to assert `> 1.2`, i.e. that the cold arm was still slow against the measured
    # COURSE even after clearing the ceiling — the bound that survived. D-223 closed it: the
    # ratio is 0.94x. The concern the old assert encoded was that a NON-uniform closure would
    # mean the defect was a temperature response after all, and that is what is checked now,
    # against the same 12 C column: the cold arm has to land near the measured mean, not merely
    # somewhere. §14 asserts the same thing across all four temperatures; here it is pinned at
    # the one column that had a second, independent bound on it.
    assert 0.90 <= d12 / _foster_observed_mean(12) <= 1.10, (
        f"at 12 C the model takes {d12:.4f} d against D-220's measured mean of "
        f"{_foster_observed_mean(12):.2f} d, a ratio of {d12 / _foster_observed_mean(12):.2f}x. "
        "D-223 measured 0.94x, out of sample. If the cold arm has drifted away from the measured "
        "COURSE while 15 C still fits, the closure was not uniform across temperature after all "
        "and D-217's refusal on E_a_uptake is worth re-opening"
    )


# ======================================================================================
# 14. The SECOND measured beer fermentation course, and what it says about the brief's
#     window, the model's speed, and the 30 C agreement (decision D-220)
# ======================================================================================

#: Foster et al. 2022's **Supplementary Figure S1** -- the fermentation COURSE, recovered.
#:
#: §12 and §13 read this trial through main-text Figure 2, which is four timepoint panels
#: with temperature on the x-axis. The underlying course was never seen, so every duration
#: the archive took from Foster was a CEILING read off a sample (D-218's own caveat), and the
#: archive had exactly ONE measured beer extract course on disk -- Tyrell's (D-216 scope
#: note). Supp. Fig. S1 is *"the specific gravities of small-scale wort fermentations at
#: eight different temperatures"*: **10 strains x 8 temperatures x 9 timepoints**, 0-120 h.
#:
#: **These values are a transcription, not a plot reading.** The supplementary figure is
#: VECTOR (``pdfimages -list`` finds no raster), so every point below is a path vertex mapped
#: through the axis tick labels' own coordinates -- the fixed-grid lesson does not bite here.
#: The geometry also settles the layout rather than assuming it: the connecting polyline is
#: drawn UNDODGED at the true x, so its vertices are the means, while the error bars are
#: dodged per series. Receipts and the full 720-point extraction live in
#: ``M:\claud_projects\temp\ferment\d220-second-beer-course\``.
#:
#: Only the four temperatures inside this engine's beer range are carried here, for the three
#: **Beer 1** clade commercial ale controls. St. Lucifer is the Beer 2 control and *"never
#: reached an FG < 1.01, regardless of temperature"*; the six kveik isolates are a different,
#: admixed clade. Neither is what this engine models.
FOSTER_COURSE_HOURS = (0, 6, 12, 24, 36, 48, 72, 96, 120)
FOSTER_BEER1_STRAINS = ("Cali", "VT", "Kolsch")
FOSTER_COURSE_SG: dict[tuple[str, int], tuple[float, ...]] = {
    ("Cali", 12): (1.0430, 1.0407, 1.0411, 1.0405, 1.0387, 1.0367, 1.0311, 1.0252, 1.0192),
    ("Cali", 15): (1.0429, 1.0417, 1.0416, 1.0395, 1.0363, 1.0321, 1.0238, 1.0156, 1.0116),
    ("Cali", 22): (1.0420, 1.0397, 1.0368, 1.0293, 1.0220, 1.0173, 1.0111, 1.0077, 1.0065),
    ("Cali", 30): (1.0420, 1.0383, 1.0341, 1.0246, 1.0192, 1.0148, 1.0094, 1.0067, 1.0064),
    ("VT", 12): (1.0430, 1.0409, 1.0414, 1.0411, 1.0400, 1.0384, 1.0333, 1.0292, 1.0238),
    ("VT", 15): (1.0429, 1.0421, 1.0419, 1.0407, 1.0379, 1.0346, 1.0262, 1.0191, 1.0143),
    ("VT", 22): (1.0420, 1.0392, 1.0373, 1.0315, 1.0240, 1.0190, 1.0127, 1.0092, 1.0074),
    ("VT", 30): (1.0420, 1.0382, 1.0336, 1.0233, 1.0180, 1.0147, 1.0103, 1.0072, 1.0069),
    ("Kolsch", 12): (1.0430, 1.0415, 1.0418, 1.0414, 1.0403, 1.0379, 1.0344, 1.0297, 1.0256),
    ("Kolsch", 15): (1.0429, 1.0415, 1.0413, 1.0391, 1.0358, 1.0313, 1.0227, 1.0150, 1.0102),
    ("Kolsch", 22): (1.0420, 1.0398, 1.0366, 1.0293, 1.0208, 1.0151, 1.0095, 1.0070, 1.0061),
    ("Kolsch", 30): (1.0420, 1.0371, 1.0317, 1.0192, 1.0141, 1.0110, 1.0073, 1.0061, 1.0064),
    # 42 C is carried for ONE purpose: it is where the commercial ale controls stall, and
    # that stall is the paper's own headline result. Without it nothing in this table pins
    # the colour-to-temperature ordering the recovery had to assume -- 12/15/22/30 C alone
    # are monotone, and a monotone series cannot distinguish an ascending palette from a
    # rescaled one. No claim below reads it except the consistency check.
    ("Cali", 42): (1.0448, 1.0421, 1.0406, 1.0384, 1.0374, 1.0365, 1.0358, 1.0353, 1.0356),
    ("VT", 42): (1.0448, 1.0411, 1.0396, 1.0374, 1.0364, 1.0355, 1.0352, 1.0350, 1.0349),
    ("Kolsch", 42): (1.0448, 1.0406, 1.0400, 1.0379, 1.0374, 1.0369, 1.0366, 1.0365, 1.0365),
}

#: Days for each Beer 1 control to reach the 1.010 target, from the course above.
#:
#: **The two cold columns EXTRAPOLATE past the 120 h end of the experiment, and the bias has
#: a known SIGN.** The extrapolation runs on the 96-120 h slope, and that slope is already
#: decelerating (Cali at 15 C falls 3.43 SG-points/h over 72-96 h and 1.67 over 96-120 h), so
#: the true crossings are LATER than these numbers, never earlier. Every claim below is
#: written to be strengthened, not weakened, by that. At 15 C the extrapolation is also short
#: -- Kolsch is at 0.2378 extract fraction against a 0.2331 target at the last measured point.
FOSTER_DAYS_TO_TARGET = {
    ("Cali", 12): 6.54,
    ("VT", 12): 7.53,
    ("Kolsch", 12): 8.73,
    ("Cali", 15): 5.39,
    ("VT", 15): 5.91,
    ("Kolsch", 15): 5.04,
    ("Cali", 22): 3.33,
    ("VT", 22): 3.76,
    ("Kolsch", 22): 2.91,
    ("Cali", 30): 2.89,
    ("VT", 30): 3.10,
    ("Kolsch", 30): 2.26,
}

#: The model at Foster's own wort and counted pitch, through §12's own helper so the numbers
#: stay commensurable with everything §12 and §13 pin. Measured at D-220 on the 6-minute grid.
#:
#: **Recorded, but not what the tests below read.** The claims run off
#: :func:`foster_model_days`, which integrates live, because a pinned duration cannot see a
#: model change -- it fires on pin drift instead, which is precisely what D-216 conceded
#: about its own pitch test. This dict is the value at the time of writing and one test
#: checks the live numbers against it, so a drift is still reported; but the ~1.5x claim is
#: measured on every run.
FOSTER_MODEL_DAYS_AS_RECORDED = {12: 7.1792, 15: 5.6125, 22: 3.2417, 30: 1.8125}

#: How long each arm has to be integrated for the crossing to exist at all. The 12 C arm
#: needs 30 d because the model takes 10.7 there; a shorter span returns ``inf`` and the
#: ratio silently becomes meaningless rather than red.
FOSTER_ARM_DAYS = {12: 30.0, 15: 25.0, 22: 15.0, 30: 12.0}

#: Apparent activation energy of the ENDPOINT, kJ/mol, over each measured temperature step.
#: This is a lumped coefficient of the whole course, NOT a reading of ``E_a_uptake`` -- taking
#: it for one is D-183's lesson and the exact error D-217 refused de Andres-Toro's -97 for.
#: It is recorded because the model/measured DIVERGENCE across the two steps is the finding.
FOSTER_APPARENT_EA_KJ = {
    "model": {(15, 22): 55.4, (22, 30): 54.1},
    "measured": {(15, 22): 49.5, (22, 30): 17.9},
}


@pytest.fixture(scope="module")
def foster_model_days() -> dict[int, float]:
    """Days to 1.010 at Foster's wort and counted pitch, integrated LIVE, per temperature.

    Module-scoped because four tests read it and each entry is an integration (~0.8 s). The
    grid is §12's own :data:`FOSTER_GRID_PER_HOUR`, so every duration here is commensurable
    with the ones §12 and §13 pin rather than being a second, differently-quantised measure.
    """
    pitch = cells_per_ml_to_pitch_gpl(FOSTER_PITCH_CELLS_PER_ML)
    return {
        temp: _foster_days_to_target_gravity(
            pitch, float(temp), days=span, per_hour=FOSTER_GRID_PER_HOUR
        )
        for temp, span in FOSTER_ARM_DAYS.items()
    }


def _foster_observed_days(temp_c: int) -> list[float]:
    return [FOSTER_DAYS_TO_TARGET[(s, temp_c)] for s in FOSTER_BEER1_STRAINS]


def _foster_observed_mean(temp_c: int) -> float:
    days = _foster_observed_days(temp_c)
    return sum(days) / len(days)


def _apparent_ea_kj(t_cold: int, d_cold: float, t_hot: int, d_hot: float) -> float:
    """Arrhenius activation energy implied by two endpoint DURATIONS, kJ/mol."""
    inv = 1.0 / (t_cold + 273.15) - 1.0 / (t_hot + 273.15)
    return float(np.log(d_cold / d_hot) * 8.314 / inv / 1000.0)


def test_fosters_recovered_course_is_internally_consistent():
    """The transcription's own consistency check, and it is not a restatement.

    Foster ran ONE wort per temperature, so the t=0 gravity must be identical across every
    strain panel of a given temperature -- and the recovered values agree to four decimals
    without ever being told to. A calibration error, a mis-assigned panel or a mistyped row
    would not reproduce that, because each panel's axis was fitted independently.

    The rest is direction. Apparent gravity must not RISE (fermentation does not run
    backwards), and the 120 h endpoint must be non-monotone in TEMPERATURE, because these ale
    strains have an optimum near 22-30 C. That optimum is the paper's own headline result,
    and it is what pins the colour-to-temperature ordering the recovery had to assume.
    """
    for temp in (12, 15, 22, 30):
        ogs = {FOSTER_COURSE_SG[(s, temp)][0] for s in FOSTER_BEER1_STRAINS}
        assert len(ogs) == 1, (
            f"at {temp} C the three strain panels disagree about the wort they started from "
            f"({sorted(ogs)}). Foster pitched one wort per temperature, so this is a "
            "transcription or calibration fault, not a fact about the yeast"
        )
    for (strain, temp), course in FOSTER_COURSE_SG.items():
        assert len(course) == len(FOSTER_COURSE_HOURS)
        rises = [
            (FOSTER_COURSE_HOURS[i], course[i], course[i + 1])
            for i in range(len(course) - 1)
            if course[i + 1] > course[i] + 0.0006
        ]
        assert not rises, (
            f"{strain} at {temp} C has gravity RISING by more than the figure's own "
            f"resolution at {rises}. Small rises inside 0.0006 are the plot's line width; a "
            "larger one means two series were crossed during recovery"
        )
    endpoints = {
        t: sum(FOSTER_COURSE_SG[(s, t)][-1] for s in FOSTER_BEER1_STRAINS) / 3
        for t in (12, 15, 22, 30, 42)
    }
    assert endpoints[12] > endpoints[15] > endpoints[22], (
        f"the 120 h endpoints {endpoints} no longer improve monotonically from 12 to 22 C"
    )
    assert endpoints[42] > endpoints[12], (
        f"the 120 h endpoints {endpoints} do not TURN OVER at the hot end. 42 C must be the "
        "worst column of all -- worse even than 12 C -- because the commercial ale controls "
        "stall there while the kveik strains do not, which is the paper's own headline. That "
        "turnover is the only thing in this table that pins the colour-to-temperature map; "
        "12-30 C alone are monotone and cannot tell an ascending palette from a rescaled one"
    )
    assert abs(endpoints[22] - endpoints[30]) < 0.0005, (
        f"22 and 30 C are no longer within the figure's resolution of each other "
        f"({endpoints[22]:.5f} vs {endpoints[30]:.5f}). They are a TIE at 120 h, and no claim "
        "in this section may order them -- the fermentations are simply both finished"
    )


def test_one_rate_scale_fitted_at_fifteen_closes_the_level_error_at_twelve_and_twentytwo(
    foster_model_days: dict[int, float],
):
    """D-220's prediction, run out-of-sample — and the reason this test changed its name.

    Through D-222 this asserted that the engine was **uniformly ~1.4x too slow** at every
    temperature below 30 C: D-220 measured 1.41 / 1.54 / 1.45x at 12 / 15 / 22 C and D-222's
    growth refit narrowed it to 1.31 / 1.42 / 1.34x. The NEAR-CONSTANCY was always the claim
    rather than the size, and it carried a prediction: a level error is closed by a rate scale,
    so ONE rate constant should close all three columns at once.

    D-223 re-anchored ``q_sugar_max`` to Foster's measured 15 C course and the prediction holds:
    **0.94 / 1.03 / 0.97x**. That is the measurement this test now makes, and two of the three
    columns are OUT OF SAMPLE -- only 15 C was fitted, so 12 C and 22 C landing within 6 % of the
    measured mean is a check of the model's temperature RESPONSE that the fit could not buy. It
    is the strongest evidence the archive has that D-220 §4's reading (a level error, not a
    temperature-response one) was right, and it is why D-217's refusal to re-source
    ``E_a_uptake`` still stands: the response never needed fixing.

    **What it is NOT.** It is not corroboration of the 15 C column, which is fitted and is
    reported here only so the three read on one scale. And it is not an endorsement of the model
    above 22 C: the 30 C column is the designed contrast and it now runs the other way, 0.66x --
    the model is FASTER than measured there. Real ale yeast saturates above ~22 C (the apparent
    activation energy of Foster's own endpoint collapses from ~50 to ~18 kJ/mol across that
    step) and this model holds a near-constant ~54 kJ/mol, so closing the cold columns has
    necessarily opened the hot one. That trade is the honest description of what D-223 bought:
    agreement over 12-22 C, purchased against a saturation term the model does not have.
    """
    for temp, recorded in FOSTER_MODEL_DAYS_AS_RECORDED.items():
        assert foster_model_days[temp] == pytest.approx(recorded, abs=0.02), (
            f"the model's {temp} C duration has moved from the recorded {recorded} d to "
            f"{foster_model_days[temp]:.4f} d. That is not a failure by itself -- the ratio "
            "asserts below are the claim -- but the record's numbers are now stale and the "
            "D-220 entry should say so"
        )
    ratios = {t: foster_model_days[t] / _foster_observed_mean(t) for t in (12, 15, 22, 30)}
    for temp in (12, 15, 22):
        assert 0.90 <= ratios[temp] <= 1.10, (
            f"at {temp} C the model takes {foster_model_days[temp]:.2f} d against a measured "
            f"mean of {_foster_observed_mean(temp):.2f} d, a ratio of {ratios[temp]:.2f}x, "
            "outside the 0.90-1.10x D-223 measured across all three. D-220 and D-222 measured "
            "1.25-1.70x here, and the whole of that gap was closed by one rate constant fitted "
            "at 15 C alone. The cold columns extrapolate past 120 h on a DECELERATING tail, so "
            "the measured means are lower bounds and these ratios can only grow -- which makes "
            "a RED on the HIGH side the more likely reading artefact and a RED on the low side "
            "a real change in the model"
        )
    # OUT OF SAMPLE. 12 and 22 C were not fitted, and they are the half of this test that can
    # fail for a reason other than arithmetic [[feedback-fit-the-observable-not-the-consequence]].
    assert abs(ratios[12] - 1.0) < 0.10 and abs(ratios[22] - 1.0) < 0.10, (
        f"the two UNFITTED columns read {ratios[12]:.3f}x at 12 C and {ratios[22]:.3f}x at "
        "22 C. A rate fitted at 15 C alone is only evidence about the temperature response if "
        "these two land near 1.0 without being aimed at; if they have drifted, D-223's claim "
        "that the response never needed fixing has to be re-argued and D-217 re-opened"
    )
    assert ratios[30] < 0.80, (
        f"the 30 C contrast now reads {ratios[30]:.2f}x, i.e. the model is FASTER than measured "
        "there. It is what makes the three ratios above a measurement instead of a predicate "
        "that fires everywhere, and since D-223 it also carries the cost of the re-anchoring: "
        "real ale yeast saturates above ~22 C and this model does not. A RED here costs the "
        "section its control AND its statement of what the re-anchoring bought"
    )
    assert min(ratios[t] for t in (12, 15, 22)) > ratios[30], (
        "every sub-30 C ratio must exceed the 30 C one for the contrast to mean anything"
    )
    sub30 = [ratios[t] for t in (12, 15, 22)]
    assert max(sub30) - min(sub30) < 0.25, (
        f"the sub-30 C ratios span {max(sub30) - min(sub30):.3f} ({min(sub30):.2f}-"
        f"{max(sub30):.2f}); D-220 measured 0.13, D-222 0.11 and D-223 0.09. NEAR-CONSTANCY "
        "across the cold columns is what made this a LEVEL error rather than a "
        "temperature-response one, which is what licensed closing it with a single rate "
        "constant. It survives the closure, and it is still why none of this re-opens D-217's "
        "refusal to re-source E_a_uptake"
    )


def test_the_thirty_degree_overshoot_is_the_saturation_term_the_model_does_not_have(
    foster_model_days: dict[int, float],
):
    """The mechanism behind the hot column — and the reason this test was renamed at D-223.

    Real ale yeast SATURATES above ~22 C: Foster's endpoints barely improve from 22 to 30 C,
    and the apparent activation energy of the measured endpoint collapses from ~50 to ~18
    kJ/mol across that step. The model has no such term and holds a near-constant ~54-56
    kJ/mol throughout, so it runs down a straight line while the data flattens.

    **Through D-222 this test was called "the thirty degree AGREEMENT is a crossing of two
    errors"**, and its job was to stop anyone citing a 30 C match as the model getting beer's
    speed right: the model was ~1.4x slow everywhere below 30 C and landed inside the measured
    spread at 30 C only because two errors cancelled there. D-223 removed the first error --
    one rate fitted at 15 C closed 12 / 15 / 22 C to within 6 % — so the cancellation is gone
    and what is left is the second error, undisguised: at 30 C the model now takes 1.81 d
    against measured crossings of 2.26-3.10 d, i.e. it is **outside the spread on the FAST
    side**. There is no agreement here to caveat any more; there is an overshoot to explain,
    and the explanation is the same missing term it always was.

    The asserts below did not change, because the mechanism did not: they compare the two
    apparent activation energies, and it is their DIVERGENCE that is the claim. Only the name
    and the reading moved [[feedback-a-hit-can-be-two-errors-cancelling]].
    """
    # The overshoot itself, pinned — the fact that replaces the retired "agreement". The old
    # name survived D-223's first sweep because this test stayed GREEN: its asserts were about
    # the activation energies and never about where the model landed, so nothing went red when
    # the thing in the title stopped being true.
    at30 = _foster_observed_days(30)
    assert foster_model_days[30] < min(at30), (
        f"at 30 C the model takes {foster_model_days[30]:.2f} d against measured crossings "
        f"{sorted(at30)}. D-223 measured it OUTSIDE the spread on the fast side (1.81 vs a "
        "2.26 fastest strain). If it is back inside, the cancellation this test is named for "
        "has returned and the reading has to be redone rather than cited"
    )
    for who, expected in FOSTER_APPARENT_EA_KJ.items():
        for (cold, hot), pinned in expected.items():
            if who == "model":
                got = _apparent_ea_kj(cold, foster_model_days[cold], hot, foster_model_days[hot])
            else:
                got = _apparent_ea_kj(
                    cold, _foster_observed_mean(cold), hot, _foster_observed_mean(hot)
                )
            assert got == pytest.approx(pinned, abs=0.6), (
                f"the {who} apparent activation energy over {cold}-{hot} C is {got:.1f} "
                f"kJ/mol against a pinned {pinned}"
            )
    measured = FOSTER_APPARENT_EA_KJ["measured"]
    model = FOSTER_APPARENT_EA_KJ["model"]
    assert measured[(22, 30)] < measured[(15, 22)] / 2.0, (
        "the measured temperature response no longer COLLAPSES above 22 C, which is the "
        "whole reason the 30 C agreement is a crossing rather than a validation"
    )
    assert abs(model[(22, 30)] - model[(15, 22)]) < 5.0, (
        "the model's temperature response is no longer near-constant across the two steps. "
        "If the model has acquired a saturating term, the 30 C agreement may have stopped "
        "being a coincidence, and this section's reading of it needs redoing"
    )


def test_the_handoff_window_is_supported_at_fifteen_degrees_and_refuted_at_twenty(
    foster_model_days: dict[int, float],
):
    """The correction to D-218 §4 and D-219 §5(c), and it cuts both ways.

    Those records concluded, unconditionally, that *nothing found supports 5-7 days*. That
    was read off timepoint panels dominated by the moderate and high temperatures. The
    temperature-resolved course does not say it: at **15 C all three Beer 1 controls land
    INSIDE 5-7 d** (5.04, 5.39, 5.91), and at 12 C two of the three are slower than the
    window's slow edge. The brief's DURATION is a real brewing figure.

    What the evidence refutes is that duration TOGETHER WITH §2.2's 20 C. Foster brackets
    20 C, and the bracket alone settles it without needing an interpolation: at 22 C the same
    strains take 2.91-3.76 d, so at 20 C they cannot take 5-7. §2.2 asserts both halves at
    once, and only the pair is wrong.

    The bracket is used rather than an Arrhenius interpolation on purpose -- interpolating
    would put a hard-sounding point estimate on top of an input that is itself an
    extrapolation.
    """
    at15 = _foster_observed_days(15)
    assert all(5.0 <= d <= 7.0 for d in at15), (
        f"the 15 C crossings are {at15}; the correction to D-218 §4 rests on all three "
        "landing inside the brief's own 5-7 d window"
    )
    at22 = _foster_observed_days(22)
    assert max(at22) < 5.0, (
        f"the 22 C crossings are {at22}. The refutation of 5-7 d AT 20 C is a bracket "
        "argument, and it needs the hot side clear of the window's fast edge"
    )
    assert max(at22) < min(at15), (
        "the 15 C and 22 C crossing sets overlap, so 20 C is no longer bracketed and the "
        "claim about §2.2's temperature would need an interpolation this section refuses"
    )
    # This used to assert the model was SLOWER than every measured strain at both bracketing
    # temperatures -- the caveat that made §2.2's window "a criterion the model passes while
    # being ~1.5x slow, rather than an honest one". D-223 removed the caveat by removing the
    # 1.5x: the model now lands INSIDE the measured spread at both, which is what the assert
    # says now. The bracket argument itself is untouched -- it is a claim about Foster's
    # measured crossings, not about the model.
    for temp, crossings in ((15, at15), (22, at22)):
        assert min(crossings) <= foster_model_days[temp] <= max(crossings), (
            f"at {temp} C the model takes {foster_model_days[temp]:.2f} d against measured "
            f"crossings {sorted(crossings)}. D-223 put it inside the spread at both bracketing "
            "temperatures; outside it, §2.2's window goes back to being a criterion the model "
            "passes for the wrong reason and the caveat this assert retired has to come back"
        )


def test_the_measured_course_peaks_later_the_colder_it_runs():
    """D-218's *"every measured course peaks on day 1-2"* is temperature-conditional.

    That claim was assembled from Tyrell plus FITTED logistic curves and generalised into a
    structural verdict on the model's day-4 peak. The independent course says the peak time
    moves with temperature: the steepest interval ends at 24 h at 30 C, at 24-36 h at 22 C,
    at 48-72 h at 15 C, and at 12 C the ferment has not peaked by the 120 h end of the run.

    This does NOT rescue the model, and the section does not claim it does: run at Foster's
    OWN conditions the engine peaks later still at every temperature. What changes is the
    reason -- the defect is that the model is uniformly slow (test above), not that real
    ferments universally peak on day 1-2 while this one does not.

    Reported as a broad plateau rather than an argmax: at 15 C Cali's 36-48, 48-72 and
    72-96 h intervals run 3.53 / 3.45 / 3.43 SG-points/h, and picking a winner out of those
    is reading noise.
    """
    hours = FOSTER_COURSE_HOURS
    peak = {}
    for temp in (12, 15, 22, 30):
        spans = []
        for strain in FOSTER_BEER1_STRAINS:
            course = FOSTER_COURSE_SG[(strain, temp)]
            rates = [
                (course[i] - course[i + 1]) / (hours[i + 1] - hours[i])
                for i in range(len(course) - 1)
            ]
            best = max(rates[1:])  # the 0-6 h interval is start-up scatter, not fermentation
            within = [hours[i + 1] for i, r in enumerate(rates) if i >= 1 and r >= 0.90 * best]
            spans.append((min(within), max(within)))
        peak[temp] = spans
    assert all(lo >= 36 for lo, _ in peak[15]), (
        f"at 15 C the measured plateau now starts before 36 h ({peak[15]}), which is the "
        "reading that makes D-218's day-1-2 generalisation temperature-conditional"
    )
    assert all(hi <= 36 for _, hi in peak[30]), (
        f"at 30 C the measured plateau now runs past 36 h ({peak[30]}); the ordering "
        "cold-peaks-later is what this test asserts, and it needs the hot end to stay early"
    )
    assert all(hi >= 96 for _, hi in peak[12]), (
        f"at 12 C the ferment now peaks before the run ends ({peak[12]}). The 12 C column is "
        "the strongest form of the claim: not merely a later peak, but no peak inside five "
        "days"
    )


# ======================================================================================
# 12. The wort's free amino-acid side chains — the buffer that LEAVES with the yeast
#     (decision D-239)
# ======================================================================================
#
# `peptide_buffer_capacity_beer` is back-solved so the model's wort reproduces Peyer's published
# BC = 1.18. Until D-239 that back-solve ran on eight organic acids and the peptide lump and
# nothing else, so the free amino acids' share of a real wort's buffering had nowhere to go but
# into the lump — and the lump is permanent while the amino acids are eaten inside 24 h. The
# split holds the wort's TOTAL at 1.18 and gives the two halves their real fates.


def _d239_data_dir_without_the_term(tmp_path: Path) -> Path:
    """A parameter dir with D-239 UNWIRED — the pre-split model, reconstructed coherently.

    Both halves have to move together or the arm is not the old model, it is a double-count:
    zeroing the three ratios alone would leave `nitrogen_uptake_charge_beer` carrying the three
    side chains at their FULLY PROTONATED charge with nothing subtracting their dissociation, so
    the pool would shed 0.234 mol+ per mole N instead of 0.177 and the arm would be MORE acidic
    than either model [[feedback-gate-both-halves-of-a-pair]].
    """
    dest = tmp_path / "d239_off"
    shutil.copytree(default_data_dir(), dest)
    path = dest / "acidbase.yaml"
    text = path.read_text(encoding="utf-8")
    for name in acidbase.AMINO_BUFFER_RATIO_PARAMS.values():
        text, n = re.subn(
            rf"^({name}:\n  value: )[-0-9.eE]+$", r"\g<1>0.0", text, flags=re.MULTILINE
        )
        assert n == 1, f"{name}'s block shape moved; fix this helper"
        text, n = re.subn(
            rf"^({name}:\n(?:  (?!uncertainty)[^\n]*\n)*  uncertainty: )\{{[^\n]*\}}$",
            r'\g<1>{ low: 0.0, high: 0.0, note: "D-239 unwired" }',
            text,
            flags=re.MULTILINE,
        )
        assert n == 1, f"{name}'s band shape moved; fix this helper"
    text, n = re.subn(
        r"^(nitrogen_uptake_charge_beer:\n  value: )[-0-9.eE]+$",
        r"\g<1>0.1772",
        text,
        flags=re.MULTILINE,
    )
    assert n == 1
    text, n = re.subn(
        r"^(nitrogen_uptake_charge_beer:\n(?:  (?!uncertainty)[^\n]*\n)*  uncertainty: )"
        r"\{[^\n]*\}$",
        r'\g<1>{ low: 0.1665, high: 0.1880, note: "D-239 unwired" }',
        text,
        flags=re.MULTILINE,
    )
    assert n == 1
    path.write_text(text, encoding="utf-8")
    # The capacity moves back too: the lump re-absorbs what the split took out of it.
    cap = dest / "beer_acids.yaml"
    cap_text = cap.read_text(encoding="utf-8")
    cap_text = cap_text.replace("1.4350620340127729", "1.5480662315921656")
    cap.write_text(cap_text, encoding="utf-8")
    return dest


def test_the_wort_amino_buffer_is_a_reading_of_the_nitrogen_pool_and_leaves_with_it():
    """Proportional to ``N``, zero without it — the claim that makes this not a second pool.

    The three side chains hold no state slot; their concentration is a composition ratio times
    the ``N`` slot. Two consequences are asserted rather than assumed, because between them they
    are the whole reason no nitrogen is booked twice and no Process is needed:

    * doubling the wort's nitrogen doubles all three, EXACTLY (a ratio, not a fit);
    * a finished beer, whose nitrogen the yeast has taken, carries none of it.

    The second is the term's entire content. A guard that only checked the wort would pass on an
    implementation that seeded the three species and then never drained them — which is the
    permanent-lump defect this beat exists to fix, wearing the new term's name.
    """
    params = _beer_acid_params().resolve()
    one = acidbase.amino_buffer_from_gpl(0.1, "beer", params)
    two = acidbase.amino_buffer_from_gpl(0.2, "beer", params)
    assert set(one) == set(acidbase.AMINO_BUFFER_SPECS), "all three species or none"
    for name, conc in one.items():
        assert two[name] == pytest.approx(2.0 * conc, rel=1e-15), (
            f"{name} is not proportional to the nitrogen pool; it is a composition ratio times "
            "`N`, so twice the nitrogen must be exactly twice the species"
        )
    assert acidbase.amino_buffer_from_gpl(0.0, "beer", params) == dict.fromkeys(
        acidbase.AMINO_BUFFER_SPECS, 0.0
    ), "a wort with no assimilable nitrogen carries none of these"

    # ...and on a real ferment, which is where the claim actually bites.
    compiled = compile_scenario(
        Scenario(
            name="d239-drains",
            medium="beer",
            initial=dict(TYRELL_SCENARIO),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=14.0,
        )
    )
    res = compiled.run()
    resolved = compiled.parameters.resolve()
    states = np.asarray(res.y, dtype=float)
    t_h = np.asarray(res.t, dtype=float)
    at_start = acidbase.amino_buffer_molar(states[:, 0], compiled.schema, resolved)
    day7 = np.array([np.interp(7 * 24.0, t_h, states[i, :]) for i in range(states.shape[0])])
    at_day7 = acidbase.amino_buffer_molar(day7, compiled.schema, resolved)
    assert all(v > 1e-5 for v in at_start.values()), (
        f"the WORT carries no amino buffering ({at_start}); the term is inert from the start "
        "and every number this beat measured would be a null"
    )
    assert all(v < 1e-9 for v in at_day7.values()), (
        f"day 7 still carries {at_day7}; the pool must leave with the yeast, and a term that "
        "stays is the permanent-lump defect D-239 exists to remove"
    )


def test_the_amino_buffer_split_holds_the_pools_net_charge_where_D209_measured_it():
    """The decomposition is exact at the wort pH — the guard against a silent double-count.

    ``nitrogen_uptake_charge_beer`` now carries the three side chains at their fully-protonated
    charge and the balance carries their dissociation. At the pH the parameter is defined at,
    the two halves must cancel term for term and leave D-209's 0.1772 mol+ per mole N. Nothing
    downstream would raise if they did not: the anchor would simply absorb the difference and
    every beer would start at the right pH carrying the wrong charge, which is exactly the shape
    of defect D-234 found by census rather than by a red test.

    Asserted on a COMPILED state rather than on the parameters, because it is the balance's
    arithmetic that has to cancel, not the YAML's.
    """
    compiled = compile_scenario(
        Scenario(
            name="d239-identity",
            medium="beer",
            initial=dict(TYRELL_SCENARIO),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=1.0,
        )
    )
    resolved = compiled.parameters.resolve()
    schema = compiled.schema
    nitrogen_molar = float(compiled.y0[schema.slice("N")][0]) / M_NITROGEN

    cation_side = acidbase.nitrogen_charge_molar(compiled.y0, schema, resolved)
    pka = acidbase.build_pka_map(resolved)
    h_wort = 10.0 ** (-TYRELL_WORT_PH)
    anion_side = sum(
        conc * acidbase.mean_charge(h_wort, pka[name])
        for name, conc in acidbase.amino_buffer_molar(compiled.y0, schema, resolved).items()
    )
    net_per_mole_n = (cation_side - anion_side) / nitrogen_molar
    assert net_per_mole_n == pytest.approx(0.1772, abs=5e-4), (
        f"the pool's NET charge at the wort pH is {net_per_mole_n:.6f} per mole N; D-209 derived "
        "0.1772 and D-239 only re-partitioned it. A move here means the charge is being counted "
        "twice on one side of the balance"
    )
    # And it MOVES away from the wort pH, which is the physics the split buys. Frozen at 0.1772
    # before D-239; D-209 §8 measured the real pool at +0.188 by pH 4.86 and called the freeze an
    # understatement.
    h_beer = 10.0 ** (-4.86)
    anion_beer = sum(
        conc * acidbase.mean_charge(h_beer, pka[name])
        for name, conc in acidbase.amino_buffer_molar(compiled.y0, schema, resolved).items()
    )
    net_at_beer_ph = (cation_side - anion_beer) / nitrogen_molar
    assert net_at_beer_ph > net_per_mole_n, (
        "the pool's charge must RISE as pH falls (its side chains re-protonate); if it does not, "
        "the term has the wrong sign and is buffering the wrong way"
    )
    assert net_at_beer_ph == pytest.approx(0.1877, abs=5e-4), (
        f"at pH 4.86 the pool carries {net_at_beer_ph:.4f} per mole N; D-209 §8 measured the "
        "real 18-species pool at +0.188 there, so three species reproduce essentially all of a "
        "drift that used to be quoted as evidence the frozen value understated"
    )


def test_wine_carries_no_wort_amino_buffer_and_its_balance_is_untouched():
    """Beer-only, and the reason is measurement rather than scope convenience.

    The ratios are a malt wort's composition, and the same three species in Huang & Ough's must
    are worth 0.73 % of wine's own acid buffering against beer's 6.7 % — wine is an order of
    magnitude better buffered (D-209 §9's 48 vs 2-5 mEq/L/pH). Wine also speciates its amino
    acids as state slots already (D-100) while keeping them charge-inactive, so closing the same
    gap there is a different act on a different wiring.

    D-209 §9 priced wine BEFORE letting its term help beer; this asserts the structural half of
    that, which a number cannot: no wine state can reach these species at all.
    """
    params = _beer_acid_params().resolve()
    assert acidbase.amino_buffer_from_gpl(0.3, "wine", params) == {}, (
        "a must must never carry the wort's amino-buffer ratios — they are a malt composition"
    )
    compiled = compile_scenario(
        Scenario(
            name="d239-wine",
            medium="wine",
            initial={
                "brix": 22.0,
                "yan_mgl": 250.0,
                "initial_ph": 3.4,
                "pitch_gpl": 0.25,
                "tartaric_gpl": 6.0,
            },
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    resolved = compiled.parameters.resolve()
    totals = acidbase._totals_molar(compiled.y0, compiled.schema, resolved)
    assert not (set(totals) & set(acidbase.AMINO_BUFFER_SPECS)), (
        f"a wine charge balance carries {set(totals) & set(acidbase.AMINO_BUFFER_SPECS)}; "
        "D-239 is beer-only and a wine that reads these is a silently changed medium"
    )
    assert acidbase.ph_of_state(compiled.y0, compiled.schema, resolved) == pytest.approx(
        3.4, abs=1e-6
    ), "the wine anchor must be exactly what it was"


def test_an_unanchored_beer_gets_no_amino_buffer_either():
    """The D-179 gate, on the new term — an empty balance must stay empty.

    An un-anchored beer carries 200-odd mg/L in ``N`` with every acid slot at 0 (``_beer_acids``
    seeds them from ``initial_ph`` or not at all). ``nitrogen_charge_molar`` is gated on
    :func:`acidbase.charge_balance_is_populated` because ungated it would hand that empty balance
    +0.0025 mol/L of cation and answer pH ~11. The amino term is the mirror image — three ANIONS
    with no cation to meet them — and needs the same gate for the same reason.
    """
    compiled = compile_scenario(
        Scenario(
            name="d239-ungated",
            medium="beer",
            initial={k: v for k, v in TYRELL_SCENARIO.items() if k != "initial_ph"},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=1.0,
        )
    )
    resolved = compiled.parameters.resolve()
    assert acidbase.amino_buffer_molar(compiled.y0, compiled.schema, resolved) == {}, (
        "an un-anchored beer got amino buffering; nitrogen is not pH information on its own and "
        "an unpopulated balance must not be fabricated from it"
    )
    assert acidbase.ph_of_state(compiled.y0, compiled.schema, resolved) == pytest.approx(
        7.0, abs=1e-9
    ), "an empty balance still reads 7.0, byte-for-byte the pre-D-179 beer"


def test_the_amino_buffer_is_in_the_ph_solve_and_out_of_the_TA_equivalents_sum():
    """D-209 §8c's asymmetry, extended — and here the reason is arithmetic, not a scruple.

    :func:`acidbase.titratable_acidity` approximates each acid's contribution as
    ``protons - mean_charge``, i.e. every proton down to the fully dissociated species. On an
    amino acid that would count the α-amino proton at pKa 8.8-9.9, which does not come off before
    a titration's 8.2 endpoint — so including these would invent equivalents. They stay in the
    pH SOLVE, where they belong, and out of the SUM.

    Checked by reconstructing the sum from the same state and the same pH: a test that merely
    compared two TA numbers could not tell "excluded from the sum" from "worth very little".
    """
    compiled = compile_scenario(
        Scenario(
            name="d239-ta",
            medium="beer",
            initial=dict(TYRELL_SCENARIO),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
            duration_days=1.0,
        )
    )
    resolved = compiled.parameters.resolve()
    schema = compiled.schema
    shipped = acidbase.titratable_acidity(compiled.y0, schema, resolved)

    pka = acidbase.build_pka_map(resolved)
    totals = acidbase._totals_molar(compiled.y0, schema, resolved)
    byp = acidbase._byp_succinic_molar(compiled.y0, schema)
    ph = acidbase.solve_ph(
        totals, acidbase._cation(compiled.y0, schema, resolved), byp, 0.0, pka
    )
    h = 10.0 ** (-ph)
    rebuilt = byp * (acidbase.BYP_AS_SUCCINIC.protons - acidbase.mean_charge(h, pka["Byp"]))
    with_them = rebuilt
    for name, conc in totals.items():
        if name in acidbase.AMINO_BUFFER_SPECS:
            with_them += conc * (
                acidbase.AMINO_BUFFER_SPECS[name].protons - acidbase.mean_charge(h, pka[name])
            )
            continue
        term = conc * (acidbase.ALL_ACIDS[name].protons - acidbase.mean_charge(h, pka[name]))
        rebuilt += term
        with_them += term

    assert shipped == pytest.approx(rebuilt * (M_TARTARIC / 2.0), rel=1e-12), (
        "beer's TA is not the sum with the three side chains EXCLUDED; if they have joined the "
        "equivalents sum, the endpoint approximation is inventing an α-amino proton that no "
        "titration to pH 8.2 removes"
    )
    assert shipped != pytest.approx(with_them * (M_TARTARIC / 2.0), rel=1e-9), (
        "including them would give the SAME answer, so this test cannot tell the two apart and "
        "proves nothing [[feedback-a-non-vacuity-check-can-itself-be-vacuous]]"
    )


def test_the_amino_buffer_costs_the_day_7_course_and_almost_nothing_at_day_1(tmp_path):
    """The term's own price, against the model it replaced — both halves unwired together.

    The comparison arm restores `nitrogen_uptake_charge_beer` to D-209's 0.1772, zeroes the three
    ratios and puts the peptide capacity back to its pre-split root. That is the ONLY coherent
    way to switch this off: the split moved charge between two places, so unwiring one place
    alone leaves the pool shedding 0.234 per mole N instead of 0.177.

    The shape is the finding, and it is stronger than the shape
    ``test_losing_wort_protein_acidifies_late_not_early`` predicted for any buffer-removal term
    two records before this one was built. That test said such a term is worth little early and
    much later. Measured, this one is worth **+0.0023 pH at day 1 — the OTHER SIGN** — and the
    reason is a number this file already pins elsewhere: only **0.298** of the wort's nitrogen is
    drawn by 24 h (``test_the_day_1_pH_miss_survives_the_timing_fix_and_has_CHANGED_SIDES``, and
    it is inside Tyrell's own cell-count spread). So at day 1 the pool is still **70 % present**
    and its three side chains are still buffering, while the share of the permanent lump they
    replaced has already gone — the wort is briefly better buffered than before the split, not
    worse. Only from day 2, when uptake has run (9.2 % of the pool left), does the sign settle
    negative: -0.0140 at day 2, -0.0202 by day 7, flat after.

    **This depends on the uptake calendar and would invert if that moved.** D-209 §7 measured
    uptake as >99 % complete by 24 h, and on THAT calendar the pool would already be gone at day
    1 and the early sign would be negative. D-211 and D-222 re-derived it to 0.298. A beat that
    moves beer's uptake timing again must re-measure this sign rather than cite it
    [[feedback-a-calibrated-level-decays-when-anything-upstream-moves]].

    **So the term makes beer's one missed day slightly WORSE**, since day 1 misses on the
    alkaline side. That is recorded here rather than buried in the record, because it is the
    whole reason this beat shipped on fidelity rather than on agreement
    [[feedback-closer-to-reality-decides]] — and because a later reader looking for the day-1
    candidate D-222 reserved will find it spent, with the sign against them.
    """
    off = _d239_data_dir_without_the_term(tmp_path)
    for day, expected in ((0.0, 0.0), (1.0, +0.0023), (2.0, -0.0140), (7.0, -0.0202)):
        with_term = _tyrell_degassed_ph_at_day(default_data_dir(), day)
        without = _tyrell_degassed_ph_at_day(off, day)
        assert with_term - without == pytest.approx(expected, abs=0.0015), (
            f"day {day:.0f}: the amino-buffer split moves the course by "
            f"{with_term - without:+.4f} pH, D-239 measures {expected:+.4f}"
        )
    # Day 0 is EXACTLY zero, not merely small, and that is the structural claim: the t=0 cation
    # back-solve absorbs a species that is present when the anchor is taken (D-178's phosphate
    # result), so a re-partition of the wort cannot move the pH the scenario asked for. If this
    # ever became non-zero the anchor would be absorbing the term only approximately, which is
    # how a charge error hides [[feedback-the-setting-where-a-change-is-exact-is-the-control]].
    assert _tyrell_degassed_ph_at_day(default_data_dir(), 0.0) == pytest.approx(
        _tyrell_degassed_ph_at_day(off, 0.0), abs=1e-9
    ), "the anchored t=0 pH must be identical with and without the split"
