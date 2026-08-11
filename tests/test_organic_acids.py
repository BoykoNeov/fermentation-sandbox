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

**BOTH of D-180's omitted terms are now built, and that changes what a shortfall means.**
Against a measured drop of **0.81** pH — the mean of the extreme strains, which is what
``measured_drop`` below computes; D-180's prose quotes the four-strain mean 0.8125 and the two
must not be mixed — the model gives **77.6-97.0 %** at nominal across the sampled
``pKa_peptide_buffer`` band, and **63.8-109.4 % over the joint band** of all NINE drawn
quantities. The history of those two numbers is the whole story of this axis: 63-92 % and
41-105 % at D-180 (with a corner reaching the measurement), 42.7-62.2 % and 7.6-82.2 % at
D-181 (with nothing reaching), and now back above both — because D-181 removed an error that
was propping the agreement up and D-182 supplied the term that was genuinely missing.

**This still is not validation.** The nominal falls short by 3-22 %; the corner that reaches
is a corner of a 9-dimensional hypercube, not a draw anyone observed; and the two shape
failures section 9 records — acetic's mid-ferment transient and lactic's late rise — are
unmodelled and are not charge-balance terms. No test here is named or phrased as validating
the produced acids alone.
"""

import pytest

from fermentation.core import acidbase
from fermentation.core.acidbase import charge_balance_is_populated
from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.kinetics import (
    ORGANIC_ACID_SPECS,
    WORT_ACID_SINKS,
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
#: Fig. 4, same ferments: wort pH and the four strains' day-7 beer pH.
TYRELL_WORT_PH = 5.65
TYRELL_BEER_PH = (4.78, 4.90)
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

    yields = {spec.yield_param for spec in ORGANIC_ACID_SPECS}
    declared = {name for p in wine.active for name in p.reads}
    assert not (yields & declared), (
        f"a wine Process declares beer's organic-acid yields {sorted(yields & declared)}; "
        "that would put them in wine's sampled set and shift its draw sequence"
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
    for spec in ORGANIC_ACID_SPECS:
        assert float(res.series(spec.slot)[-1]) == 0.0, (
            f"{spec.slot} grew on a beer that supplied no initial_ph"
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
    # And it really does move all four acids plus sugar — an inert Process would pass the
    # subset check above vacuously.
    assert touched == {*(spec.slot for spec in ORGANIC_ACID_SPECS), "S"}


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
    Tyrell's Fig. 4 is a genuine external test. Two scopes, and conflating them is the whole
    trap this test exists to avoid:

    * **at nominal, across ``pKa_peptide_buffer``'s band: 77.6-97.0 %** of the measured
      0.81 pH drop (D-182; 42.7-62.2 % at D-181, 63-92 % at D-180). The model must still fall
      SHORT here.
    * **over the JOINT band — NINE drawn quantities, not one: 63.8-109.4 %.**

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

    **What each dimension is worth, measured** (all-nominal fraction 0.8559, whole band each):
    the peptide pKa moves it 0.776-0.970; the yields and floors comparably; the three seeds
    ~0.015; ``pKa_oxalic_2``/``pKa_pyruvic`` **0.0003**; and of D-182's three, the carbonic
    pKa is worth ~0.01 and the two solubility parameters ~0.01 between them — the CO2 term's
    SIZE is consequential but its band is not, because both edges sit within 10 % of a
    nominal that is itself a printed in-beer measurement.

    **A CORNER REACHES THE MEASUREMENT AGAIN, AND IT WAS PREDICTED IN ADVANCE.** D-180's did
    (104.5 %), and that reach belonged to the falling acids' absence rather than to the model;
    D-181 removed it and wrote that "a future change that makes one reach again is a signal to
    find out which omitted term arrived"
    [[feedback-a-margin-is-a-claim-about-what-holds-it-open]]. One has: dissolved CO2, the
    last of the two terms D-180 named, built at D-182 and pre-registered at 76-104 % before a
    line of it was written. **That does not make the model validated.** The nominal still
    falls short by 3-22 %, the reaching member is a corner of a 9-D hypercube nobody was seen
    to draw, and the acetic transient and lactic late rise (§9) remain unmodelled.

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
    # from what it produced, so this cannot drift from the run.
    flux = (float(res.series("acetic")[-1]) - float(res.series("acetic")[0])) / params[
        "Y_acetic_sugar_beer"
    ]

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
    ) -> float:
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
        cation = acidbase.solve_cation_charge(seeded, 0.0, 0.0, member, TYRELL_WORT_PH)
        start = acidbase.solve_ph(seeded, cation, 0.0, 0.0, member)
        assert start == pytest.approx(TYRELL_WORT_PH, abs=1e-6), (
            "every member must start at the supplied wort pH — that is what anchoring means"
        )
        end = dict(seeded)
        for spec in ORGANIC_ACID_SPECS:
            # Produced acids build on the RUN's seed, not the varied one: a yield is a measured
            # production and does not depend on how much of a DIFFERENT acid the wort carried.
            end[spec.slot] = (
                start_molar[spec.slot] + (edge(spec.yield_param, pick) * flux) / molar[spec.slot]
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
        return (start - acidbase.solve_ph(end, cation, 0.0, carbonic_molar, member)) / measured_drop

    pka_band = (
        beer_params["pKa_peptide_buffer"].uncertainty.low,
        beer_params["pKa_peptide_buffer"].value,
        beer_params["pKa_peptide_buffer"].uncertainty.high,
    )

    # Scope 1 — everything nominal but the peptide pKa. The "must fall short" claim lives HERE.
    at_nominal = [
        fraction(pka, "nom", "nom", "nom", "nom", "nom", "nom", "nom", "nom") for pka in pka_band
    ]
    assert min(at_nominal) > 0.70, (
        f"the predicted drop collapsed to {min(at_nominal):.0%} of Tyrell's measured one at "
        "nominal yields; D-182 measured 77.6-97.0 % across the pKa band (D-181's 42.7-62.2 % "
        "was the same model with no dissolved CO2 in its charge balance, and D-180's 63-92 % "
        "was that model also unable to lose the falling acids' charge)"
    )
    assert max(at_nominal) < 1.0, (
        f"at NOMINAL yields the predicted drop reached {max(at_nominal):.0%} of the measured "
        "one. It is still supposed to fall short: BOTH of the terms D-180 named as omitted "
        "have now been built (D-181's falling acids, D-182's dissolved CO2), so reaching "
        "100 % here no longer has a pending omission to explain it — it would mean the model "
        "gained acidification from somewhere unaccounted. What is still missing is named in "
        "this module's header: acetic's transient and lactic's late rise, neither of which is "
        "a charge-balance term."
    )

    # Scope 2 — the joint band the sampler can actually reach. NO upper bound is asserted
    # here, because a corner legitimately exceeds the measured drop; what is pinned is the
    # SPAN, so that a change which narrows or shifts it has to be looked at. D-181 adds the
    # three FLOORS as a third band dimension: their edges are named strains' own day-7 values,
    # so leaving them at nominal here would repeat the point-vs-band mistake one level out —
    # which is exactly what D-180's amendment had to correct in this very test.
    picks = ("lo", "nom", "hi")
    joint = [
        fraction(pka, pick, floor_pick, ox2, pyr, seed, carbonic, henry, vant_hoff)
        for pka in pka_band
        for pick in picks
        for floor_pick in picks
        for ox2 in picks
        for pyr in picks
        for seed in picks
        for carbonic in picks
        for henry in picks
        for vant_hoff in picks
    ]
    assert len(joint) == 3**9, "every drawn dimension must be varied, not a subset of them"
    assert min(joint) == pytest.approx(0.638, abs=0.02), (
        f"the joint low corner moved to {min(joint):.1%}; D-182 measured 63.8 % (D-181's was "
        "7.6 %) — yields at their low edge, peptide pKa HIGH, floors at their LOW edge (the "
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
    assert max(joint) == pytest.approx(1.094, abs=0.02), (
        f"the joint high corner moved to {max(joint):.1%}; D-182 measured 109.4 % (D-181's "
        "was 82.2 %) — yields high, peptide pKa low, floors high, seeds low, and the three "
        "CO2 parameters at the edges that dissolve the most and dissociate it hardest. NB "
        "this is a CORNER of a 9-D hypercube, not a member any ensemble was seen to draw. "
        "**A corner REACHES the measured drop again, and it was pre-registered.** D-180 had "
        "one at 104.5 %; D-181 removed it and wrote that a future change restoring it would "
        "signal that an omitted term had arrived. One has — dissolved CO2, D-180's own arm C, "
        "predicted at 76-104 % before D-182 was written. Reaching is therefore the EXPECTED "
        "outcome here and not evidence the model is right: the nominal still falls short, and "
        "no upper bound is asserted on this scope for the same reason it never was."
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

    **D-182 HALVED THE NUMBER, AND THAT IS A RESULT RATHER THAN A REGRESSION.** D-181 measured
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

    # The same two end states, re-solved with the carbonic term switched off. This must
    # reproduce D-181's +0.2094 pH: if it does, the shipped shrinkage is the CO2 buffer; if it
    # does not, something in the seeds, floors or charge arithmetic actually moved.
    def co2_free(compiled, res) -> float:
        params = compiled.parameters.resolve()
        y = res.y[:, -1]
        return acidbase.solve_ph(
            acidbase._totals_molar(y, compiled.schema),
            acidbase._cation(y, compiled.schema),
            acidbase._byp_succinic_molar(y, compiled.schema),
            0.0,
            acidbase.build_pka_map(params),
        )

    co2_free_gap = co2_free(c1, r1) - co2_free(c2, r2)
    assert 0.15 < co2_free_gap < 0.35, (
        f"with the carbonic term off the missing base is worth {co2_free_gap:.4f} pH, outside "
        "D-180's predicted +0.2-0.3 window that D-181 landed at 0.2094. The D-182 buffer "
        "explanation for the shipped 0.1128 only holds if this arm still reproduces D-181."
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
