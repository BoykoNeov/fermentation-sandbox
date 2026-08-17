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
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

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
    **15 % attenuated** while acetic has already made **86 %** of its whole rise. So a yield on
    the sugar flux — which is what ``Y_acetic_sugar_beer`` is, a measured day7−day0 difference
    divided by a sugar divisor — must put the acid far too late. This recomputes exactly how
    late, from the flux the shipped run actually booked (recovered through *lactic*, which still
    rides that flux), rather than trusting a number in a comment.

    It is a claim about SHAPE, not size: the same yield reproduces the day-7 endpoint by
    construction, which is precisely why the error was invisible until the figure interiors were
    read.
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
    assert shipped_rmse < 0.6 * flux_rmse, (
        f"the growth-linked form scores {shipped_rmse:.1f} ppm RMSE against Tyrell's days 1-7 "
        f"and the retired flux-linked one {flux_rmse:.1f}; D-183 measured 32.5 vs 61.6"
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


def test_the_day_1_pH_is_still_missed_and_the_miss_is_uptake_timing():
    """THE REMAINING DEFECT, pinned with its attribution — and D-209 made it WORSE.

    Day 7 now lands inside Tyrell's envelope; day 1 does not, and it is further out than before
    this beat: **0.31 pH too acidic at nominal, against 0.19 too alkaline at D-208**. That is not
    the charge arithmetic being wrong, it is *when* the charge leaves. The cation change is the
    integral of nitrogen uptake, and this model empties its ``N`` slot inside about 20 hours, so
    the entire step is delivered before the day-1 reading is taken while the measured pH is only
    ~63 % of the way to its day-7 value.

    The attribution is measured here two ways rather than asserted: the uptake fraction by 24 h,
    and a counterfactual in which the SAME total charge is delivered on a slower linear ramp.
    **At a 60-hour ramp all eight measured days land inside the envelope**, and the window is
    narrow — 48 h leaves day 1 out and 72 h leaves day 2 out — so 60 h is a measurement on this
    trajectory and not a knob with slack in it. Against the model's own ~20 h that is roughly
    threefold too fast.

    **The ramp is NOT a fix and must not become one.** Nothing sources 60 hours; nitrogen-uptake
    timing in this engine is set by ``GrowthNitrogenLimited``'s constants, which were never
    calibrated against a wort FAN time course. Recorded as a located defect — the next beat on
    beer's pH is an uptake-timing question, not an acid-base one.
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
    assert day1 < lo - TYRELL_PH_READ_TOL, (
        f"day 1 reads {day1:.4f}, inside or above Tyrell's {lo:.3f}-{hi:.3f}. D-209 measures it "
        "0.31 pH BELOW — too acidic. If this passes the defect has been fixed and the docstring "
        "above is stale"
    )
    assert lo - day1 == pytest.approx(0.315, abs=0.03), (
        f"the day-1 miss is {lo - day1:.4f} pH; D-209 measures 0.315, and D-208 measured 0.186 "
        "in the other direction. Both numbers matter: the miss GREW, and a beat that shrinks it "
        "should say which"
    )

    # Attribution 1 — the nitrogen is gone before the reading, so the whole step has landed.
    nitrogen = states[schema.slice("N").start, :]
    consumed_by_24h = float(np.interp(24.0, t_h, nitrogen[0] - nitrogen))
    assert consumed_by_24h / (nitrogen[0] - min(nitrogen)) > 0.99, (
        "the attribution rests on the charge step being complete by the day-1 reading"
    )

    # Attribution 2 — the same total charge on a slower ramp. A COUNTERFACTUAL computed on the
    # shipped trajectory: it prices where the miss lives, and it is NOT a proposal. 60 hours is
    # unsourced; what makes it evidence rather than a fit is that the neighbouring ramps fail.
    total_charge = acidbase.nitrogen_charge_from_gpl(float(nitrogen[0]), "beer", params)

    def days_inside_on_ramp(ramp_hours: float) -> int:
        inside = 0
        for day, (band_lo, band_hi) in TYRELL_PH_COURSE.items():
            y = state_at(day * 24.0)
            delivered = total_charge * min(day * 24.0 / ramp_hours, 1.0)
            actual = acidbase.nitrogen_charge_from_gpl(
                float(np.interp(day * 24.0, t_h, nitrogen)), "beer", params
            )
            # The balance's cation side is ``slot + actual``; the counterfactual wants
            # ``slot + total_charge - delivered``, which at t=0 is the anchored total and once the
            # ramp completes is the bare slot. So the slot is nudged by the difference.
            y[schema.slice("cation_charge")] += (total_charge - delivered) - actual
            ramped = float(acidbase.degassed_ph_of_state(y, schema, params))
            inside += band_lo - TYRELL_PH_READ_TOL <= ramped <= band_hi + TYRELL_PH_READ_TOL
        return inside

    assert days_inside_on_ramp(60.0) == 8, (
        "on a 60 h counterfactual ramp all eight measured days must land inside the envelope — "
        "that is what locates the day-1 miss in uptake TIMING rather than in the charge "
        "arithmetic. The shipped model's own uptake finishes in ~20 h and reaches 7 of 8"
    )
    # The neighbours FAIL, which is what stops 60 h reading as a knob with slack in it: a
    # counterfactual that worked across a wide range would price nothing.
    assert days_inside_on_ramp(48.0) == 7, "48 h must still leave day 1 out"
    assert days_inside_on_ramp(72.0) == 7, "72 h must overshoot the other way and leave day 2 out"
