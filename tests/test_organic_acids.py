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

**The headline is a compensation, not an agreement, and these tests are written to keep it
that way.** Against a measured drop of 0.75-0.87 pH the model gives 0.514-0.744 across the
sampled ``pKa_peptide_buffer`` band — but that is held open by two omitted terms of opposite
sign (three wort acids that fall and are not state slots, ~+0.2-0.3 pH; dissolved CO₂, ~−0.3
pH). No test here is named or phrased as validating the produced acids alone.
"""

import pytest

from fermentation.core import acidbase
from fermentation.core.acidbase import charge_balance_is_populated
from fermentation.core.chemistry import carbon_mass_fraction, sugar_species
from fermentation.core.kinetics import ORGANIC_ACID_SPECS, OrganicAcidExcretion
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
    cation = acidbase.solve_cation_charge(totals, 0.0, pka, 4.4)

    honest = acidbase.solve_ph(totals, cation, 0.0, pka)
    # The same succinic, half of it also booked through the Byp lump.
    doubled = acidbase.solve_ph(
        {"succinic": totals["succinic"] * 0.5}, cation, totals["succinic"] * 0.5, pka
    )
    assert honest == pytest.approx(doubled, abs=1e-9), (
        "sanity: splitting succinic between its slot and Byp must be pH-neutral"
    )
    # ...which is precisely why ADDING a Byp on top of a full succinic slot is not.
    inflated = acidbase.solve_ph(totals, cation, totals["succinic"], pka)
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


def test_the_predicted_ph_drop_is_short_of_the_measured_one_across_the_whole_pka_band(
    beer_params,
):
    """THE FREE PREDICTION — and the test is written to hold the shortfall, not the agreement.

    Nothing in ``beer_acids.yaml`` is fitted to pH, so comparing the modelled trajectory with
    Tyrell's Fig. 4 is a genuine external test. It lands at **63-92 %** of the measured drop
    depending on where ``pKa_peptide_buffer`` sits in its **sampled** band — a factor of ~1.45
    spread, which is why this is pinned over the band and not at the nominal
    [[feedback-pin-the-band-not-the-nominal]]. The peptide term, not any acid yield, is the
    parameter the prediction is most sensitive to.

    **The band arm RE-ANCHORS the cation, and getting that wrong is the trap.** Re-reading the
    shipped trajectory's pH at a different ``pKa_peptide_buffer`` while holding the
    ``cation_charge`` the compile back-solved at the NOMINAL pKa measures the wrong thing: it
    reports 72-80 %, and it also moves the START pH off the 5.65 the scenario supplied, which
    no ensemble member ever does. A sampled member draws its pKa and the compile then
    back-solves ITS cation to hit ``initial_ph``, so every member starts at 5.65 and they
    differ only in where they end. That is what this reproduces — the recurring
    verified-at-a-point / sampled-over-a-band shape, in the one place it would have quietly
    understated the spread.

    THE UPPER BOUND IS THE POINT. The model must NOT reproduce the measured drop, because two
    real terms of opposite sign are missing (three falling wort acids that are not state slots,
    ~+0.2-0.3 pH; dissolved CO₂, ~−0.3 pH). If a future change makes this agree exactly, that
    is a signal to go and find which of the two got added — not a success
    [[feedback-a-margin-is-a-claim-about-what-holds-it-open]].
    """
    compiled, res = _run(dict(TYRELL_SCENARIO))
    params = compiled.parameters.resolve()
    pka_map = acidbase.build_pka_map(params)
    slots = ("lactic", "acetic", "citrate", "malic", "succinic", "peptide_buffer")
    molar = {s: acidbase.ALL_ACIDS[s].molar_mass for s in slots}
    start_molar = {s: float(res.series(s)[0]) / molar[s] for s in slots}
    end_molar = {s: float(res.series(s)[-1]) / molar[s] for s in slots}
    band = (
        beer_params["pKa_peptide_buffer"].uncertainty.low,
        beer_params["pKa_peptide_buffer"].value,
        beer_params["pKa_peptide_buffer"].uncertainty.high,
    )
    measured_drop = TYRELL_WORT_PH - sum(TYRELL_BEER_PH) / 2.0

    fractions = []
    for pka in band:
        member = {**pka_map, "peptide_buffer": (pka,)}
        # What compile_scenario does for THIS member: anchor the cation to the wort pH it
        # was given, using the member's own pKa.
        cation = acidbase.solve_cation_charge(start_molar, 0.0, member, TYRELL_WORT_PH)
        start = acidbase.solve_ph(start_molar, cation, 0.0, member)
        end = acidbase.solve_ph(end_molar, cation, 0.0, member)
        assert start == pytest.approx(TYRELL_WORT_PH, abs=1e-6), (
            "every member must start at the supplied wort pH — that is what anchoring means"
        )
        fractions.append((start - end) / measured_drop)

    assert min(fractions) > 0.55, (
        f"the predicted pH drop collapsed to {min(fractions):.0%} of Tyrell's measured one; "
        "D-180 measured 63-92 % across this band"
    )
    assert max(fractions) < 1.0, (
        f"the predicted drop reached {max(fractions):.0%} of the measured one. It is supposed "
        "to fall SHORT: the falling wort acids (pyruvic/formic/oxalic) and dissolved CO2 are "
        "both absent, and they pull in opposite directions. Find out which one arrived."
    )
    # And the direction that makes it a beer at all: the pH must actually fall.
    assert fractions[1] > 0.0


def test_a_beer_seeded_at_finished_beer_levels_would_overshoot(beer_params):
    """Why the seeds HAD to move — kept as a test so the reason cannot be forgotten.

    D-179 dosed the ``*_typical_beer`` levels at pitch. Bolting a producer onto those seeds
    would add a whole ferment's production on top of an already-finished composition. This
    recomputes that counterfactual on the charge balance: it lands below any real beer, which
    is what makes "produce the acids" and "start from a wort" one decision rather than two.
    """
    params = beer_params.resolve()
    pka = acidbase.build_pka_map(params)
    molar = {s: acidbase.ALL_ACIDS[s].molar_mass for s in ("lactic", "acetic", "malic", "succinic")}

    for seeds, expect_below in (
        (
            {
                s: beer_params[f"{'citric' if s == 'citrate' else s}_typical_beer"].value
                for s in molar
            },
            4.4,
        ),
    ):
        totals = {s: g / molar[s] for s, g in seeds.items()}
        totals["peptide_buffer"] = (
            beer_params["peptide_buffer_capacity_beer"].value
            / acidbase.ALL_ACIDS["peptide_buffer"].molar_mass
        )
        cation = acidbase.solve_cation_charge(totals, 0.0, pka, expect_below)
        produced = dict(totals)
        for spec in ORGANIC_ACID_SPECS:
            produced[spec.slot] = (
                produced.get(spec.slot, 0.0)
                + (params[spec.yield_param] * TYRELL_SUGAR_GPL) / molar[spec.slot]
            )
        assert acidbase.solve_ph(produced, cation, 0.0, pka) < 4.3, (
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
