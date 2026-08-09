"""The scenario-override admissible range (decision D-164).

Two wine scenario knobs override a reference parameter at the compile seam:
``carrying_capacity_gpl`` -> ``biomass_carrying_capacity`` (D-30) and
``autolysis_rate_per_h`` -> ``k_autolysis`` (D-34). Both build the override as a new
``Parameter`` carrying the base's uncertainty band, so the band has always bounded how far
a "sweep" may go — but only as a SIDE EFFECT of ``Parameter._value_in_range``, whose error
names the parameter's epistemic band and never the knob the user set. D-163 §5 found this
channel; D-164 states the bound explicitly and pins it here.

What these tests forbid, in order:

1. A silent widening of the admissible range — an override outside the band must RAISE.
2. An error that does not name the knob. The pre-D-164 pydantic message named only
   ``biomass_carrying_capacity`` / ``k_autolysis``; a user who set ``carrying_capacity_gpl``
   was told about a parameter they had not heard of. Asserting only "it raises" would pass
   on that old message and pin nothing.
3. A moved band EDGE. The four edge values below are PRINTED — read verbatim from
   ``wine_generic.yaml`` — not constructed. D-163 measured that 652 of 678 archive band
   edges move in both directions with the full suite green; these four are guarded because
   they are load-bearing for a live gate, not merely descriptive.
4. A narrowed band that silently invalidates a shipped scenario. ``test_autolysis`` runs at
   EXACTLY the 1e-2 high edge, so the closed interval is itself the contract.
"""

from __future__ import annotations

import pytest

from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario

# PRINTED from src/fermentation/parameters/data/wine_generic.yaml (biomass_carrying_capacity
# :177, k_autolysis :438) — the shipped bands, transcribed, not derived.
_BANDS = {
    "carrying_capacity_gpl": ("biomass_carrying_capacity", 2.0, 5.0),
    "autolysis_rate_per_h": ("k_autolysis", 1.0e-4, 1.0e-2),
}


def _wine(**knobs: float) -> Scenario:
    """A minimal compilable wine scenario; ``knobs`` land in ``initial``."""
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25}
    initial.update(knobs)
    return Scenario(
        name="override-bounds",
        medium="wine",
        initial=initial,
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=14.0, celsius=25.0),
        ],
        interventions=[],
        duration_days=30.0,
    )


@pytest.mark.parametrize("knob", sorted(_BANDS))
def test_shipped_band_edges_are_what_the_gate_is_built_on(knob):
    """The PRINTED edges still match the YAML. Guards the band, not just the nominal."""
    pname, low, high = _BANDS[knob]
    base = compile_scenario(_wine()).parameters[pname]
    assert (base.uncertainty.low, base.uncertainty.high) == (low, high), (
        f"{pname}'s band moved. It is the admissible range for scenario knob {knob!r} "
        f"(D-164), so moving it silently changes what a scenario may request."
    )


@pytest.mark.parametrize("knob", sorted(_BANDS))
def test_override_at_either_exact_edge_is_admissible(knob):
    """The interval is CLOSED at both ends.

    Not pedantry: the shipped ``test_autolysis`` suite runs ``autolysis_rate_per_h=1e-2``,
    exactly the high edge. An open interval, or a band narrowed by a later re-sourcing,
    breaks that scenario for reasons that have nothing to do with autolysis physics.
    """
    _pname, low, high = _BANDS[knob]
    for edge in (low, high):
        compiled = compile_scenario(_wine(**{knob: edge}))
        assert compiled is not None


@pytest.mark.parametrize("knob", sorted(_BANDS))
def test_override_inside_the_band_compiles(knob):
    """Designed-GREEN control: without it, a harness broken for an unrelated reason
    (a missing temperature_schedule, say) would make every RAISED arm below confirm
    itself."""
    _pname, low, high = _BANDS[knob]
    compiled = compile_scenario(_wine(**{knob: 0.5 * (low + high)}))
    assert compiled is not None


@pytest.mark.parametrize("knob", sorted(_BANDS))
@pytest.mark.parametrize("side", ["below", "above"])
def test_override_outside_the_band_raises_and_names_the_knob(knob, side):
    pname, low, high = _BANDS[knob]
    bad = low * 0.5 if side == "below" else high * 2.0

    with pytest.raises(ValueError) as excinfo:
        compile_scenario(_wine(**{knob: bad}))
    msg = str(excinfo.value)

    # (2) above: the knob the USER set must appear. The pre-D-164 message named only the
    # parameter, so this assertion is the whole point of the test.
    assert knob in msg, f"error must name the scenario knob the user set; got: {msg}"
    assert pname in msg, f"error must name the overridden parameter; got: {msg}"
    assert "admissible" in msg, f"error must name what it forbids; got: {msg}"
    assert "D-164" in msg, f"error must cite the decision; got: {msg}"
    # The remedy must be actionable: widen the band in YAML, with provenance.
    assert "provenance" in msg


@pytest.mark.parametrize("knob", sorted(_BANDS))
def test_a_hair_outside_the_edge_is_rejected(knob):
    """The gate is on the edge itself, not somewhere loosely near it."""
    _pname, low, high = _BANDS[knob]
    with pytest.raises(ValueError):
        compile_scenario(_wine(**{knob: high * (1.0 + 1e-6)}))
    with pytest.raises(ValueError):
        compile_scenario(_wine(**{knob: low * (1.0 - 1e-6)}))


@pytest.mark.parametrize("knob", sorted(_BANDS))
def test_the_overridden_parameter_is_actually_ensemble_sampled(knob):
    """The gate's stated REASON must be true, or it gets argued away.

    ``_override_in_band``'s message justifies the bound by the ensemble sampler's
    ``triangular(low, value, high)`` precondition. That is only a live mechanism if the
    parameter is in fact drawn — ``sample_parameters`` skips any name whose band is
    zero-width. Both are (measured: 247 sampled names in a compiled wine run).
    """
    pname, low, high = _BANDS[knob]
    params = compile_scenario(_wine(**{knob: 0.5 * (low + high)})).parameters
    p = params[pname]
    assert p.uncertainty.high > p.uncertainty.low, (
        f"{pname} has a zero-width band, so the sampler would skip it and D-164's "
        f"stated justification for the gate would be false."
    )
    assert p.uncertainty.low <= p.value <= p.uncertainty.high
