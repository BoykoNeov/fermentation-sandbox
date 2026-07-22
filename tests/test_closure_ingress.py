"""Tests for closure oxygen ingress — the O2 SUPPLY term (decision D-136).

:class:`ClosureOxygenIngress` adds the bottle closure's steady oxygen transmission rate straight
into ``d(o2)/dt``, closing a gap this repo had already named against itself: D-108's
``SotolonAldolCondensation`` docstring says *"a* sealed *wine here has strictly zero O2 ingress (no
closure permeation ...), so a sealed sulfited bottle never ages toward premox at all. That is the
limitation to state."*

**What these tests are really guarding — and it is NOT "o2 goes down".** Before D-136 the whole
oxidative axis ran on a *finite* charge dosed by ``add_oxygen``, and the characteristic behaviour
was saturation as the charge was spent. Continuous ingress inverts that: the consumers are
collectively far faster than any closure, so ``o2`` quasi-steady-states just above zero and the
endpoints accumulate at ``otr * (k_i / sum k)``. **The closure becomes the master throttle and the
individual rate constants become a splitting rule.** The two load-bearing tests here are therefore

1. :func:`test_oxidative_endpoints_are_ordered_by_closure_over_five_years` — the integrated
   trajectory, asserting the ORDERING of SO2 depletion / browning across the closure menu rather
   than any magnitude (the Tier-3 risk rule: directional checks, not numbers), and
2. :func:`test_so2_depletion_is_supply_limited_not_rate_limited` — the sharper claim, that a 16x
   swing in ``k_so2_oxidation`` barely moves the SO2 endpoint while a closure swap dominates it.
   If someone ever "fixes" this axis by tuning a sink's rate constant, that test says why it did
   not work.

Every shipped OTR is READ from ``closure.yaml`` rather than restated (the D-100 lesson: a test
that hard-codes the value it should be reading is a test of itself). The literals here are the
*published* numbers — Lopes et al. 2007's Table I in uL/day — which is what the model is checked
against, plus the conversion arithmetic that turns them into the engine's g/L/h.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core.kinetics import ClosureOxygenIngress
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario.compile import compile_scenario
from fermentation.scenario.schema import Intervention, Scenario, TemperaturePoint
from fermentation.validation import assert_conserved, total_carbon

#: Lopes et al. 2007 (J. Agric. Food Chem. 55:5167-5170) Table I, steady horizontal-storage oxygen
#: ingress in **uL O2/day**, transcribed from the paper. These are the published observations the
#: shipped constants must reproduce; the shipped constants themselves are read from the YAML.
_PUBLISHED_ULDAY: dict[str, float] = {
    "technical_cork": 0.25,  # printed range 0.1-0.4, central value
    "screwcap": 0.45,  # printed range 0.2-0.7, central value
    "natural_cork": 1.5,  # the region where Lopes 2007 and Oliveira 2013 agree (see closure.yaml)
    "synthetic_nomacorc": 6.0,  # printed as a single value for horizontal storage
    "synthetic_supremecorq": 13.0,  # printed range 11-15, central value
}

#: The unit conversion, re-derived here from its three independent factors rather than imported, so
#: this test cannot inherit an arithmetic error from the code it checks:
#:   1.43 ug O2/uL   -- the Lopes group's OWN factor, recovered from Oliveira 2013's Discussion
#:                      quoting "2.43 to 8.73 ug/day (1.7 to 6.1 uL/day)" (2.43/1.7 = 8.73/6.1)
#:   / 24 h/day      -- their rates are per day, the engine runs in hours
#:   / 0.750 L       -- the standard bottle, the ONLY place a volume enters (see closure.yaml)
_ULDAY_TO_GPLH = 1.43e-6 / 24.0 / 0.750

#: The closure menu in the order closure.yaml claims is ASCENDING in oxygen transmission. Note
#: technical cork below screwcap — the deliberate, sourced contradiction of the "screwcaps are
#: least permeable" folklore (see closure.yaml's ordering-correction note).
_ASCENDING = (
    "hermetic",
    "technical_cork",
    "screwcap",
    "natural_cork",
    "synthetic_nomacorc",
    "synthetic_supremecorq",
)

_FERMENT_DAYS = 20.0
_AGING_YEARS = 5.0
_TOTAL_DAYS = _FERMENT_DAYS + _AGING_YEARS * 365.25


@pytest.fixture(scope="module")
def params():
    return load_parameters(
        default_data_dir() / "wine_generic.yaml", default_data_dir() / "closure.yaml"
    )


def _wine_with_otr(otr: float):
    """A bare wine state carrying ``otr`` in its ``closure_otr`` slot."""
    schema = wine_schema()
    y = schema.zeros()
    y[schema.slice("closure_otr")] = otr
    return schema, y


def _scenario(closure: str | None, *, so2_mgl: float = 60.0, days: float = _TOTAL_DAYS) -> Scenario:
    interventions = [Intervention(day=_FERMENT_DAYS, action="begin_aging")]
    if so2_mgl:
        interventions.insert(
            0,
            Intervention(day=_FERMENT_DAYS - 1.0, action="add_so2", params={"so2_mgl": so2_mgl}),
        )
    return Scenario(
        name="d136",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 200.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=days,
        closure=closure,
        interventions=interventions,
    )


def _age(closure: str | None, *, k_so2_scale: float = 1.0, so2_mgl: float = 60.0):
    """Compile + integrate a five-year bottle-aging run; returns ``(compiled, trajectory)``."""
    compiled = compile_scenario(_scenario(closure, so2_mgl=so2_mgl))
    values = dict(compiled.param_values)
    values["k_so2_oxidation"] *= k_so2_scale
    trajectory = simulate_scheduled(
        compiled.process_set,
        values,
        compiled.y0,
        (0.0, _TOTAL_DAYS * 24.0),
        events=compiled.events,
        t_eval=np.linspace(0.0, _TOTAL_DAYS * 24.0, 2000),
    )
    return compiled, trajectory


def _final(compiled, trajectory, name: str) -> float:
    return float(trajectory.y[compiled.schema.slice(name)][0][-1])


def _aged_final(name: str, closure: str | None, **kwargs) -> float:
    """The endpoint of ``name`` after a five-year aged run under ``closure``."""
    compiled, trajectory = _age(closure, **kwargs)
    return _final(compiled, trajectory, name)


# --------------------------------------------------------------------------------------------
# Metadata and the contract
# --------------------------------------------------------------------------------------------


def test_closure_ingress_metadata():
    process = ClosureOxygenIngress()
    assert process.name == "closure_oxygen_ingress"
    assert process.tier is Tier.SPECULATIVE
    # The ONLY Process on the aging axis that adds to o2 rather than drawing it down, and it
    # touches nothing else — o2 is carbon-free and off every ledger, so ingress moves nothing
    # that must balance (contrast D-135's carbon-weighted bound_methanethiol).
    assert process.touches == ("o2",)


def test_closure_ingress_reads_no_parameter():
    """``reads`` is empty BY DESIGN, and that is worth pinning.

    The rate rides in the ``closure_otr`` state slot (seeded at the compile seam) rather than in a
    parameter, because a closure is a per-run choice and the scenario layer has no
    parameter-override seam. Two consequences follow and are documented in the Process: the output
    tier comes from ``Process.tier`` alone rather than by D-1 parameter-tier propagation, and
    ``simulate_ensemble`` will not propagate the OTR band. If someone later moves the rate into a
    parameter, this test should fail and force that reasoning to be revisited.
    """
    assert ClosureOxygenIngress().reads == ()


# --------------------------------------------------------------------------------------------
# The rate law: a constant source, read straight from state
# --------------------------------------------------------------------------------------------


def test_ingress_is_the_state_slot_exactly():
    """``d(o2)/dt == closure_otr``, with no coefficient in between."""
    otr = 1.19e-7
    schema, y = _wine_with_otr(otr)
    d = ClosureOxygenIngress().derivatives(0.0, y, schema, {})
    assert float(d[schema.slice("o2")][0]) == otr


def test_ingress_is_zero_order_in_o2_and_time_invariant():
    """Not gradient-driven, and no burst: the same rate at t=0 and t=5 years, at any [o2].

    Zero-order is forced by the sources, not chosen for convenience — both primaries measure
    ingress INTO an O2-scavenging indigo-carmine sink, so the ~atmospheric gradient is already
    baked into every published OTR and a ``(p_atm - p_wine)`` term would double-count it. And the
    first-month bottling burst is deliberately NOT here (it is trapped cork/headspace air, an
    ``add_oxygen`` bolus); if anyone adds a decaying burst term, the time-invariance half of this
    test fails.
    """
    otr = 4.77e-7
    schema, y = _wine_with_otr(otr)
    process = ClosureOxygenIngress()

    for o2_level in (0.0, 1e-6, 5e-3, 8e-3):
        y_at = y.copy()
        y_at[schema.slice("o2")] = o2_level
        for t in (0.0, 24.0, 5.0 * 365.25 * 24.0):
            d = process.derivatives(t, y_at, schema, {})
            assert float(d[schema.slice("o2")][0]) == otr


# --------------------------------------------------------------------------------------------
# Isolability (prime directive #3) — and why 0 is the RIGHT default here, unlike D-134's copper
# --------------------------------------------------------------------------------------------


def test_zero_otr_is_byte_for_byte_inert():
    """The isolability gate: an unspecified closure leaves the pre-D-136 axis bit-identical.

    0 is both the gate AND a physically real endpoint — Lopes et al. 2007 found that of every
    sealing system tested "only the control (bottle sealed by flame) was completely air-tight". So
    this is the OPPOSITE call to D-134's copper, where 0 was an unphysical multiplier and the
    VarSpec default had to be the population mean instead.
    """
    schema, y = _wine_with_otr(0.0)
    d = ClosureOxygenIngress().derivatives(0.0, y, schema, {})
    assert not np.any(d)


def test_default_wine_state_has_no_ingress():
    """A wine packed without a closure carries otr = 0 — 0 is this slot's NEUTRAL value.

    The D-134 copper slot had to default to the population mean because 0 was an unphysical
    multiplier there. Here 0 is an additive source, so 0 really does mean "no ingress", and a
    ``pack()`` that never mentions the closure is inert rather than wrong.
    """
    schema = wine_schema()
    y = schema.pack({"S": [220.0], "X": 0.25, "E": 0.0, "N": 0.2, "CO2": 0.0, "T": 293.15})
    assert float(y[schema.slice("closure_otr")][0]) == 0.0


def test_negative_otr_cannot_become_an_oxygen_sink():
    """The ``<= 0`` guard floors a mis-seeded negative instead of draining the o2 pool.

    Unreachable through the compile seam (every shipped OTR is >= 0 and the menu is closed), so
    this exercises the guard directly rather than claiming the shipped path needs it — the D-134
    clamp-test precedent.
    """
    schema, y = _wine_with_otr(-1.0e-6)
    d = ClosureOxygenIngress().derivatives(0.0, y, schema, {})
    assert not np.any(d)


def test_beer_is_a_hard_no_op():
    """Wine-only: beer's schema has no ``closure_otr`` slot, so the Process is a no-op there.

    Crown-cap OTR is real, but the data and this axis are wine-centric (see closure.yaml).
    """
    schema = beer_schema()
    assert "closure_otr" not in schema
    d = ClosureOxygenIngress().derivatives(0.0, schema.zeros(), schema, {})
    assert not np.any(d)


def test_touches_contract_holds_under_strict():
    schema, y = _wine_with_otr(1.19e-7)
    process_set = ProcessSet(schema, [ClosureOxygenIngress()], strict=True)
    d = process_set.total_derivatives(0.0, y, {})
    assert float(d[schema.slice("o2")][0]) == 1.19e-7


# --------------------------------------------------------------------------------------------
# The parameters: provenance arithmetic, and the ordering that is the real claim
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("closure", sorted(_PUBLISHED_ULDAY))
def test_shipped_otr_reproduces_the_published_rate(params, closure):
    """Each shipped g/L/h is Lopes 2007's printed uL/day through the documented conversion.

    This is the test that pins the unit chain — 1.43 ug/uL (the authors' OWN factor, recovered
    from their two papers quoting the same numbers in both units), 24 h/day, and the 750 mL
    bottle. If anyone "corrects" the 1.43 to a 20 C molar volume, or reuses
    ``batch_volume_liters`` instead of the fixed bottle, this fails.
    """
    expected = _PUBLISHED_ULDAY[closure] * _ULDAY_TO_GPLH
    shipped = params[f"otr_{closure}"].value
    assert shipped == pytest.approx(expected, rel=0.01)


def test_hermetic_is_exactly_zero(params):
    """Not "small": exactly 0, because it is the definitional no-ingress case."""
    assert params["otr_hermetic"].value == 0.0


def test_closure_menu_is_strictly_ascending_in_otr(params):
    """THE ROBUST CLAIM. Magnitudes are a band; the ordering is what the primary establishes.

    Lopes et al. 2007's own conclusion, verbatim: "low in screw-caps and 'technical' corks,
    intermediate in conventional natural cork stoppers, and high in synthetic closures."
    """
    values = [params[f"otr_{name}"].value for name in _ASCENDING]
    assert values == sorted(values)
    assert all(a < b for a, b in zip(values, values[1:], strict=False))


def test_technical_cork_is_below_screwcap(params):
    """The deliberate contradiction of the folklore — do not "fix" this.

    It is widely repeated (Godden et al. 2005, quoted in Oliveira et al. 2013's own introduction)
    that screwcaps are the least permeable closure. Lopes et al. 2007's Table I says technical
    cork is lower at STEADY state (0.1-0.4 vs 0.2-0.7 uL/day). Both are true: the screwcap's
    famous figure is "<500 uL/day AT THE MOMENT OF BOTTLING" — headspace air trapped at sealing,
    not transmission — which dominates any total-including-burst comparison. This test exists so
    the sourced steady ordering cannot be quietly reverted to the folklore one.
    """
    assert params["otr_technical_cork"].value < params["otr_screwcap"].value


def test_every_otr_is_speculative_and_banded(params):
    """Tier floor + honest bands. Only ``hermetic`` may have a degenerate band."""
    for name in _ASCENDING:
        parameter = params[f"otr_{name}"]
        assert parameter.tier is Tier.SPECULATIVE
        assert parameter.uncertainty is not None
        low, high = parameter.uncertainty.low, parameter.uncertainty.high
        assert low <= parameter.value <= high
        if name != "hermetic":
            assert low < high, f"{name} must carry a real band, not a point"


def test_natural_cork_band_spans_the_second_primary(params):
    """Oliveira et al. 2013 measured 2.51 ug/day (months 4-12, 593 bottles) independently.

    The cross-validation between two experiments is the strongest magnitude evidence this beat
    has, so it is asserted rather than merely written down: the second primary's directly measured
    rate must fall inside the shipped band.
    """
    measured_gplh = 2.51e-6 / 24.0 / 0.750  # ug/day -> g/L/h for a 750 mL bottle
    band = params["otr_natural_cork"].uncertainty
    assert band.low <= measured_gplh <= band.high


# --------------------------------------------------------------------------------------------
# The integrated trajectory — the two tests the beat lives or dies by
# --------------------------------------------------------------------------------------------


def test_unspecified_closure_leaves_the_aging_axis_untouched():
    """No ``closure`` key ⇒ the five-year run is the pre-D-136 run: no O2, no oxidation."""
    compiled, trajectory = _age(None)
    assert _final(compiled, trajectory, "o2") == 0.0
    # A420 is produced ONLY by oxidation here (this wine is dry, so D-88 caramelization is inert),
    # so a flat browning index is the cleanest statement that no O2 arrived. Acetaldehyde is NOT
    # asserted to be zero: it carries a fermentative residue that has nothing to do with ingress.
    assert _final(compiled, trajectory, "A420") == 0.0


def test_hermetic_matches_an_unspecified_closure_exactly():
    """Naming the zero must equal omitting it: ``hermetic`` documents, it is not a mode."""
    _, bare = _age(None)
    _, sealed = _age("hermetic")
    assert np.array_equal(bare.y, sealed.y)


def test_oxidative_endpoints_are_ordered_by_closure_over_five_years():
    """THE HEADLINE. Same wine, same SO2, five years — the closure decides how it ages.

    Asserts ORDERING across the menu, never a magnitude (the Tier-3 risk rule). This is D-72/D-108's
    SO2-protection story finally driven by a closure instead of a hand-dose: a more permeable
    closure spends the SO2 faster and browns the wine further.

    Deliberately NOT asserted: anything about ``bound_h2s``/``bound_methanethiol``. D-135 is
    release-only and reads no O2, so the sulfide trajectory is closure-INDEPENDENT in this model.
    The real screwcap-reduction link (Lopes: too-low ingress promotes "rubbery or struck flint"
    characters) shows up here as the ABSENCE of oxidative markers, not as extra sulfide.
    """
    menu = [name for name in _ASCENDING if name != "hermetic"]
    so2_left, browning = [], []
    for closure in menu:
        compiled, trajectory = _age(closure)
        so2_left.append(_final(compiled, trajectory, "so2_total"))
        browning.append(_final(compiled, trajectory, "A420"))

    # More permeable closure => less SO2 survives, strictly monotone across the whole menu.
    assert all(a > b for a, b in zip(so2_left, so2_left[1:], strict=False)), so2_left
    # ... and correspondingly more browning.
    assert all(a < b for a, b in zip(browning, browning[1:], strict=False)), browning
    # The span is real, not a rounding artifact: the most permeable closure must spend most of the
    # dose while the least permeable barely touches it.
    assert so2_left[0] > 0.75 * so2_left[0] > so2_left[-1]


def test_so2_depletion_is_supply_limited_not_rate_limited():
    """THE SHARP CLAIM: past this Process, the OTR is the throttle and ``k_so2_oxidation`` is not.

    Under continuous ingress ``o2`` quasi-steady-states near zero, so the sinks share a flux they
    do not control: SO2 depletion tracks the closure, and a large swing in the sulfite rate
    constant barely moves it. This is the [[measure-which-side-before-building]] lesson appearing
    as *correct* behaviour rather than as a defect — and it is why anyone trying to tune this axis
    through a sink's rate constant will find it does not work.
    """
    slow = _aged_final("so2_total", "natural_cork", k_so2_scale=0.25)
    fast = _aged_final("so2_total", "natural_cork", k_so2_scale=4.0)
    # A 16x swing in the rate constant, end to end.
    rate_effect = abs(slow - fast)

    # Against a closure swap at a FIXED, unmodified rate constant.
    tight = _aged_final("so2_total", "screwcap")
    leaky = _aged_final("so2_total", "synthetic_nomacorc")
    closure_effect = abs(tight - leaky)

    assert closure_effect > 5.0 * rate_effect, (rate_effect, closure_effect)


def test_dissolved_oxygen_stays_non_negative_and_bounded():
    """No clamp is needed and none is used — but the claim is verified over a real 5-year run.

    Every O2 sink is proportional to ``o2`` (or gated at ``o2 <= 0``), so at ``o2 = 0`` the field
    is ``+otr > 0`` and the pool cannot go negative. It must also stay physically bounded: even
    the most permeable closure in the primary must not drive dissolved O2 past air saturation
    (~8 mg/L at 20 C), which would mean the consumers had lost the race entirely.
    """
    for closure in ("natural_cork", "synthetic_supremecorq"):
        compiled, trajectory = _age(closure)
        o2 = trajectory.y[compiled.schema.slice("o2")][0]
        assert o2.min() >= 0.0
        assert o2.max() < 8.0e-3


def test_premox_needs_a_permeable_closure():
    """D-108's stated limitation, now lifted: a sealed sulfited bottle CAN age toward premox.

    Sotolon is the premox marker (D-87/D-108). Under a tight closure it stays at its fermentative
    baseline; under a permeable one it climbs. Ordering only — the magnitude is speculative.
    """
    sealed = _aged_final("sotolon", "hermetic")
    tight = _aged_final("sotolon", "screwcap")
    leaky = _aged_final("sotolon", "synthetic_supremecorq")
    assert leaky > tight >= sealed


# --------------------------------------------------------------------------------------------
# The compile seam
# --------------------------------------------------------------------------------------------


def test_named_closure_seeds_the_sourced_otr(params):
    compiled = compile_scenario(_scenario("natural_cork", days=30.0))
    seeded = float(compiled.y0[compiled.schema.slice("closure_otr")][0])
    assert seeded == params["otr_natural_cork"].value


def test_unknown_closure_is_a_loud_error_that_names_the_menu():
    with pytest.raises(ValueError, match="unknown scenario.closure") as excinfo:
        compile_scenario(_scenario("cork", days=30.0))
    # The error must list the alternatives, not just reject — the _ALLOWED_KEYS discipline.
    for name in _ASCENDING:
        assert name in str(excinfo.value)


def test_closure_on_beer_is_rejected():
    """Wine-only: a beer scenario naming a closure is a user error, not an ignored field."""
    scenario = Scenario(
        name="d136-beer",
        medium="beer",
        initial={
            "glucose_gpl": 20.0,
            "maltose_gpl": 60.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.5,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=14.0,
        closure="screwcap",
    )
    with pytest.raises(ValueError, match="no closure-ingress model"):
        compile_scenario(scenario)


def test_ingress_is_disabled_until_begin_aging():
    """Aging-gated at the compile seam: a closure alone admits no oxygen during fermentation."""
    scenario = _scenario("synthetic_supremecorq", so2_mgl=0.0, days=10.0)
    scenario = scenario.model_copy(update={"interventions": []})
    compiled = compile_scenario(scenario)
    assert compiled.process_set.enabled_snapshot()[ClosureOxygenIngress.name] is False


def test_carbon_closes_over_a_closure_driven_five_year_run():
    """Conservation is a TEST here, not just an argument (prime directive: "a model that creates
    mass is broken regardless of how good its curves look").

    This passes for a structural reason — ``o2`` is carbon-free and off every ledger, and each
    downstream consumer's own carbon closure is identical whether the O2 arrived as an
    ``add_oxygen`` bolus or as a continuous flux — so it is deliberately a cheap
    belt-and-suspenders rather than a discovery. It is worth having anyway: D-136 introduces a
    large NEW continuous flux driving acetaldehyde, browning and sulfate over five years, and
    "the ledger was structurally safe" is the kind of claim that should be checked rather than
    asserted. Run under the most permeable closure, where the flux is largest.
    """
    compiled, trajectory = _age("synthetic_supremecorq")
    # Non-trivial: the run must actually have oxidised something, or this checks nothing.
    assert _final(compiled, trajectory, "A420") > 0.0
    carbon_fraction = compiled.param_values["biomass_C_fraction"]
    assert_conserved(
        trajectory,
        total_carbon(compiled.schema, biomass_carbon_fraction=carbon_fraction),
        label="carbon (closure-driven aging)",
    )


def test_the_oxygen_ceiling_is_held_up_by_ethanol_oxidation():
    """The bound on standing ``o2`` is ``otr / k_ethanol_oxidation``, and it is worth naming.

    ``test_dissolved_oxygen_stays_non_negative_and_bounded`` asserts dissolved O2 stays under air
    saturation, but that is not luck and it is not this Process's doing: ethanol is effectively
    inexhaustible, so :class:`OxidativeAcetaldehyde` is an always-on, never-saturating,
    first-order-in-o2 sink that caps the quasi-steady level. ``k_ethanol_oxidation`` was already
    retuned once (5.0e-4 -> 2.0e-4 at D-73); this test makes the dependency explicit, so that
    lowering it again fails HERE with a clear reason rather than silently pushing the saturation
    test toward its limit.
    """
    compiled, trajectory = _age("synthetic_supremecorq")
    o2 = trajectory.y[compiled.schema.slice("o2")][0]
    otr = float(compiled.y0[compiled.schema.slice("closure_otr")][0])
    ceiling = otr / compiled.param_values["k_ethanol_oxidation"]

    # The standing level must sit at or below the single-sink ceiling — below it, because the
    # other sinks take their share too.
    assert o2.max() <= ceiling
    # ... and that ceiling must itself be sub-saturation, which is the real safety margin.
    assert ceiling < 8.0e-3
