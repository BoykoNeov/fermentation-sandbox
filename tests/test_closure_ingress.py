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
from fermentation.core.kinetics.o2_partition import o2_depletion_shares
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

#: The UNCERTAINTY BAND edges, in the same published uL/day units (decision D-162). The nominals
#: above have been pinned since D-136; these were not, and a mutation showed any edge could be
#: silently moved with the suite still green. Each entry records whether the edge is a number P1
#: PRINTED or a construction this file made — the two are not interchangeable and the difference
#: is exactly what the ordering-scoping note in closure.yaml turns on.
_PUBLISHED_BAND_ULDAY: dict[str, tuple[float, float]] = {
    # P1 prints horizontal 0.1-0.4 AND vertical 0.1-0.9 for technical corks. The high edge is
    # deliberately the VERTICAL ceiling so the band spans both orientations P1 measured.
    "technical_cork": (0.1, 0.9),
    # P1's full printed horizontal range. No vertical figure exists — Table I prints "--" — so
    # unlike technical cork this band covers ONE orientation. That asymmetry is not a defect
    # (it is what the source measured) but it is why the two bands are not commensurable.
    "screwcap": (0.2, 0.7),
    # The union of P1's two natural-cork windows: 1.7-6.1 (months 2-12) and 0.1-2.3 (12-36).
    "natural_cork": (0.1, 6.1),
    # P1 prints a single 6 horizontal and 8-9 vertical. The high edge is printed; the LOW edge
    # of 5 is this file's own modest extension below the single horizontal point, not a P1 number.
    "synthetic_nomacorc": (5.0, 9.0),
    # P1's printed horizontal range; its vertical 11-12 sits inside it.
    "synthetic_supremecorq": (11.0, 15.0),
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

#: Dissolved-O2 air saturation, ~8 mg/L at 20 C. Two DIFFERENT quantities are checked against it
#: here and they are not interchangeable (D-172 s6): the single-sink CEILING ``otr/k_ethanol`` is a
#: BOUND, and ``o2.max()`` is the STATE. A ceiling above this means the bound stops discriminating;
#: a state above it is unphysical.
_AIR_SATURATION_GPL = 8.0e-3


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


@pytest.mark.parametrize("closure", sorted(_PUBLISHED_BAND_ULDAY))
def test_shipped_band_edges_reproduce_the_published_range(params, closure):
    """The BAND edges carry provenance too — and until D-162 nothing pinned them.

    The nominals have been pinned to Table I since D-136 by the test above; the band edges were
    not pinned at all. ``test_every_otr_is_speculative_and_banded`` only checks the band is
    ordered and non-degenerate, so any edge could be moved to any other value and the suite
    stayed green — measured, not assumed: silently replacing ``otr_technical_cork``'s high edge
    (P1's *vertical* 0.9 uL/day) with the *horizontal* 0.4 passed all 1460 tests. That edge is
    load-bearing, because the closure bands are built to DIFFERENT ORIENTATION SCOPES and that
    asymmetry is what D-162's scoping note about the ordering rests on.

    Each edge below is annotated with whether it is a number P1 printed or a construction this
    file made from P1's numbers, so the distinction cannot be lost by a later edit.
    """
    low_ulday, high_ulday = _PUBLISHED_BAND_ULDAY[closure]
    uncertainty = params[f"otr_{closure}"].uncertainty
    assert uncertainty is not None
    assert uncertainty.low == pytest.approx(low_ulday * _ULDAY_TO_GPLH, rel=0.01)
    assert uncertainty.high == pytest.approx(high_ulday * _ULDAY_TO_GPLH, rel=0.01)


def test_hermetic_is_exactly_zero(params):
    """Not "small": exactly 0, because it is the definitional no-ingress case."""
    assert params["otr_hermetic"].value == 0.0


def test_closure_menu_is_strictly_ascending_in_otr(params):
    """The ordering at the NOMINALS — and that scope is the whole of what this pins.

    Lopes et al. 2007's own conclusion, verbatim: "low in screw-caps and 'technical' corks,
    intermediate in conventional natural cork stoppers, and high in synthetic closures." That
    sentence establishes three TIERS. This chain is strictly stronger: it adds an ordering
    *within* the low tier and one *within* the synthetic tier, neither of which the sentence
    makes, and the two are not equally well supported by the shipped bands. D-162 measured all
    three claims on joint draws — see ``closure.yaml``'s header for the scoping and the figures.
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

    SCOPE, per D-162: this is a claim about the two NOMINALS, and only about them. It is a
    refinement read off Table I, not part of the three-tier conclusion the primary states, and the
    two bands overlap heavily enough that independent draws invert it often. That is not a live
    defect — one run holds one closure and the sampler cannot reach these — but the assertion
    should not be read as carrying the bands with it.
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
        assert o2.max() < _AIR_SATURATION_GPL


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
    first-order-in-o2 sink that caps the quasi-steady level. The ethanol route's constant
    was already retuned once (5.0e-4 -> 2.0e-4 at D-73) and reparameterised again at D-172
    (it is now ``k_o2_depletion_total * f_ethanol_o2_share``); this test makes that explicit, so
    that lowering EITHER factor fails HERE with a clear reason rather than silently pushing the
    saturation test toward its limit.

    **That "fails HERE" is a claim about the NOMINALS, and the band breaches it without anyone
    editing anything** (D-172 amendment §3/§6; the recurring shape is a constraint verified at a
    POINT where the sampler reads a BAND). Two DIFFERENT quantities breach, with opposite
    attributions, and conflating them is the trap this paragraph exists to close:

    * **The CEILING below (a single-sink BOUND) — pre-existing, doubled by D-172.** With
      ``otr = 1.03e-6``: nominal ``k_ethanol = 2.00e-4`` gives 5.15e-3, PASS; the post-D-172 joint
      band edge (``k_o2_depletion_total`` low 1.0e-4 x ``f_ethanol_o2_share`` low 0.2 = 2.0e-5)
      gives 5.15e-2. But the retired ``k_ethanol_oxidation``'s own band low 4.0e-5 already gave
      2.58e-2, 3.2x over, so D-172 did not introduce this — it doubled it. A ceiling above 8.0e-3
      means the BOUND stops being informative, NOT that o2 gets there.
    * **The ACTUAL standing ``o2`` — introduced by D-172, REPAIRED at D-173.** ``o2.max()`` at the
      post-D-172 total-low was **9.68e-3, ABOVE air saturation**, because standing o2 is ``otr /
      SUM k``: ``f_ethanol_o2_share`` cancels and only the TOTAL reaches this slot, and D-172
      lowered the reachable minimum sum from ``4.0e-5 + 2.0e-4 = 2.4e-4`` to a flat ``1.0e-4``.
      "No declared edge moved" was true of every edge and still changed what the JOINT surface
      reaches. D-173 moved that low edge back to 2.4e-4 on the measurement, and
      :func:`test_the_low_edge_of_the_o2_total_cannot_draw_a_supersaturated_wine` is what now
      holds it — read that test, not this one, for the band-level claim.

    So this assert guards the NOMINAL against edits. It still does not establish the CEILING's
    margin across the band, and that breach is untouched: at the joint edge the bound reaches
    5.15e-2 and goes uninformative. Only the state claim was repairable by an edge move.
    """
    compiled, trajectory = _age("synthetic_supremecorq")
    o2 = trajectory.y[compiled.schema.slice("o2")][0]
    otr = float(compiled.y0[compiled.schema.slice("closure_otr")][0])
    # D-172: `k_ethanol_oxidation` is no longer an entry -- it is the ethanol half of the
    # always-on total, so the single-sink ceiling is formed from the same helper the Process uses.
    ceiling = otr / o2_depletion_shares(compiled.param_values)[0]

    # The standing level must sit at or below the single-sink ceiling — below it, because the
    # other sinks take their share too.
    assert o2.max() <= ceiling
    # ... and that ceiling must itself be sub-saturation, which is the real safety margin.
    assert ceiling < _AIR_SATURATION_GPL


def _o2_max_at_total(total: float, *, oxidative: str = "direct") -> float:
    """``o2.max()`` over a five-year run with the always-on O2 total forced to ``total``.

    Compiles fresh every call ON PURPOSE: ``begin_aging`` leaves the aging Processes ENABLED on
    the ``ProcessSet`` it fires against, so a reused compiled object starts the next integration
    from a different enabled set and every arm after the first is a different scenario.
    """
    compiled = compile_scenario(_scenario("synthetic_supremecorq"), oxidative=oxidative)
    values = dict(compiled.param_values)  # a property returning a fresh dict -- mutate the copy
    values["k_o2_depletion_total"] = total
    trajectory = simulate_scheduled(
        compiled.process_set,
        values,
        compiled.y0,
        (0.0, _TOTAL_DAYS * 24.0),
        events=compiled.events,
        t_eval=np.linspace(0.0, _TOTAL_DAYS * 24.0, 2000),
    )
    return float(trajectory.y[compiled.schema.slice("o2")][0].max())


def test_the_low_edge_of_the_o2_total_cannot_draw_a_supersaturated_wine():
    """D-173: the band's LOW edge must not put dissolved O2 above air saturation. The STATE, not
    the bound.

    This is the guard the D-172 amendment left open. ``k_o2_depletion_total``'s shipped low was
    1.0e-4, and at that draw ``o2.max()`` reached 9.68e-3 under the most permeable closure on the
    menu — a *sampled field able to draw a physically impossible state*, which no green suite
    catches because nothing was asserting off the nominal.

    **Three things make this a guard rather than a decoration.**

    1. The edge is READ from the compile seam, never restated. Re-lowering the YAML edge fails
       HERE, which is the whole point; a hard-coded 2.4e-4 would only pin this file to itself.
    2. MONOTONICITY is asserted, so "the low edge is the worst case" is *run* rather than assumed
       — and it doubles as the anti-self-confirmation control: if the override silently failed to
       apply, all three arms would return the same number and the ordering assert would fail.
       (`feedback-verify-the-restore-between-mutation-arms`.) It is asserted on the band's own
       GEOMETRIC INTERIOR, never on the nominal: with the nominal as the middle point, RAISING the
       low edge toward it collapses the ordering and this test would go RED on a change that makes
       the physics *safer* — a guard that forbids the safe direction is one that gets argued away
       (`feedback-name-guards-for-what-they-forbid`).
    3. Both shipped oxidative sets are covered. The cascade reads the same total as its activation
       floor and breached too (8.35e-3 at the retired low), so guarding ``direct`` alone would
       leave half the surface open.

    **The high edge moved under this test at D-174 (1.0e-3 -> 2.0e-3) and every arm here got
    safer**, which is the whole reason the guard is written off the seam: ``at_mid`` follows to
    6.93e-4 and ``at_high`` falls to 5.13e-4 g/L, 15.6x under saturation where it was 7.8x. D-174
    swept the total to 1.0e-1 — 200x the nominal — looking for a high-side channel to any
    impossible state on either oxidative set and on beer's dosed pool, and found none; the
    positive control on that same harness reproduced the 9.68e-3 breach below. So "raising this
    edge only makes the state safer" is now measured rather than conceded.

    **What this does NOT claim.** The margin at the low edge is 1.88x on this arm alone. The joint
    over the eleven OTHER sink parameters that reach ``d(o2)/dt`` closes it to **1.59x** (direct)
    and **1.41x** (cascade). Two further surfaces are not sink constants and are held out of that
    figure: ``T_ref``, which detunes every ``f(T)`` at once rather than perturbing one sink
    (including it gives 1.14x / 1.04x), and ``closure_otr``, the NUMERATOR — not sampled at all
    (``ClosureOxygenIngress.reads == ()`` by design), but at supremecorq's printed high of 15
    uL/day the corner reaches 8.16e-3, back above saturation. All measured at D-173, none guarded
    here (the corner costs 30+ integrations per arm), so what holds this margin open includes two
    decisions about the sampler surface and not only the edges below.
    """
    compiled = compile_scenario(_scenario("synthetic_supremecorq"))
    band = compiled.parameters["k_o2_depletion_total"].uncertainty
    assert band is not None and band.low is not None and band.high is not None

    at_low = _o2_max_at_total(band.low)
    assert at_low < _AIR_SATURATION_GPL, (
        f"k_o2_depletion_total's low edge {band.low:.2e} draws o2.max() = {at_low:.3e} g/L, at or "
        f"above air saturation {_AIR_SATURATION_GPL:.1e} -- an unphysical state in a sampled "
        "field. See D-173 before lowering this edge."
    )

    # Standing o2 is otr/SUM(k), so it DECREASES in the total: the low edge is the worst case, and
    # this ordering is what says so rather than assuming it. The middle point is the band's
    # GEOMETRIC centre, which is strictly interior for ANY band the schema admits -- so no edge
    # move in the safe direction can collapse the ordering (see the docstring). It also fails if
    # the override never took effect, since three equal numbers cannot be strictly ordered.
    at_mid = _o2_max_at_total(float(np.sqrt(band.low * band.high)))
    at_high = _o2_max_at_total(band.high)
    assert at_low > at_mid > at_high

    # The mildest arm is the paired GREEN control: it must clear saturation with room to spare, or
    # these arms are not separating the edge from the scenario.
    assert at_high < _AIR_SATURATION_GPL / 3.0

    # The cascade set reads the same total as its activation floor and breached identically.
    assert _o2_max_at_total(band.low, oxidative="cascade") < _AIR_SATURATION_GPL
