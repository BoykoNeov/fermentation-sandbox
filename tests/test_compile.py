"""The scenario → core compile seam: unit conversion, validation, and assembly."""

from pathlib import Path

import pytest

from fermentation.core.media import MEDIA
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.units.convert import brix_to_sugar_gpl, celsius_to_kelvin
from fermentation.validation import assert_conserved, assert_nonnegative, total_carbon

DATA = Path(__file__).resolve().parents[1] / "src" / "fermentation" / "parameters" / "data"
WINE_PARAMS = DATA / "wine_generic.yaml"


def _wine_scenario(**overrides) -> Scenario:
    kwargs: dict[str, object] = {
        "name": "wine-benchmark",
        "medium": "wine",
        "initial": {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5},
        "temperature_schedule": [TemperaturePoint(day=0.0, celsius=20.0)],
        "duration_days": 14.0,
    }
    kwargs.update(overrides)
    return Scenario(**kwargs)


def _beer_scenario(**overrides) -> Scenario:
    kwargs: dict[str, object] = {
        "name": "beer-benchmark",
        "medium": "beer",
        "initial": {
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 1.0,
        },
        "temperature_schedule": [TemperaturePoint(day=0.0, celsius=20.0)],
        "duration_days": 7.0,
    }
    kwargs.update(overrides)
    return Scenario(**kwargs)


# -- wine: full end-to-end (parameters exist) ---------------------------------


def test_wine_compiles_to_canonical_initial_state():
    compiled = compile_scenario(_wine_scenario())
    state = compiled.schema.unpack(compiled.y0)

    # Wine has a single sugar slot, so the schema reads S back as a scalar
    # (size-1 vars collapse to a float on unpack); beer's length-3 S stays a vector.
    assert state["S"] == pytest.approx(brix_to_sugar_gpl(24.0))
    assert state["N"] == pytest.approx(250.0 / 1000.0)  # mg/L -> g/L
    assert state["X"] == pytest.approx(0.5)
    assert state["E"] == 0.0  # defaulted
    assert state["CO2"] == 0.0
    assert state["T"] == pytest.approx(celsius_to_kelvin(20.0))


def test_wine_time_span_is_hours():
    compiled = compile_scenario(_wine_scenario(duration_days=14.0))
    assert compiled.t_span_h == (0.0, 14.0 * 24.0)


def test_wine_loads_provenance_parameters_by_medium_and_strain():
    compiled = compile_scenario(_wine_scenario())
    assert "mu_max" in compiled.parameters
    # Sourced from Coleman 2007 (regression evaluated at the 20 C T_ref).
    assert compiled.param_values["mu_max"] == pytest.approx(0.095)
    # param_values is the plain hot-loop mapping; parameters keeps tiers/provenance.
    assert compiled.parameters["mu_max"].provenance.doi == "10.1128/aem.00670-07"


def test_optional_ethanol_can_be_supplied():
    compiled = compile_scenario(
        _wine_scenario(
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5, "ethanol_gpl": 3.0}
        )
    )
    assert compiled.schema.get(compiled.y0, "E") == pytest.approx(3.0)


# -- beer: structural seam (no sourced params yet) ----------------------------


def test_beer_sugar_vector_is_packed_in_uptake_order():
    # This test only asserts the state layout; the explicit wine-file override keeps
    # it independent of the beer parameter values (which beer_generic.yaml now holds).
    compiled = compile_scenario(_beer_scenario(), parameter_paths=[WINE_PARAMS])
    assert compiled.schema.get(compiled.y0, "S") == pytest.approx([15.0, 70.0, 20.0])
    assert compiled.schema.get(compiled.y0, "N") == pytest.approx(200.0 / 1000.0)


def test_beer_without_param_file_reports_the_missing_source():
    # beer_generic.yaml now exists; an unsourced strain still has no file, so the
    # missing-source path is exercised via a strain we have not added.
    with pytest.raises(FileNotFoundError, match="beer_saison.yaml"):
        compile_scenario(_beer_scenario(strain="saison"))


# -- validation at the boundary -----------------------------------------------


def test_unknown_initial_key_is_rejected():
    with pytest.raises(ValueError, match="unknown key.*brixx"):
        compile_scenario(
            _wine_scenario(initial={"brixx": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5})
        )


def test_missing_required_key_is_rejected():
    with pytest.raises(ValueError, match="missing required key 'brix'"):
        compile_scenario(_wine_scenario(initial={"yan_mgl": 250.0, "pitch_gpl": 0.5}))


def test_negative_value_is_rejected():
    with pytest.raises(ValueError, match="must be >= 0"):
        compile_scenario(_wine_scenario(initial={"brix": -1.0, "yan_mgl": 250.0, "pitch_gpl": 0.5}))


def test_missing_temperature_schedule_is_rejected():
    with pytest.raises(ValueError, match="temperature_schedule needs at least one point"):
        compile_scenario(_wine_scenario(temperature_schedule=[]))


def test_unknown_medium_is_rejected():
    with pytest.raises(KeyError, match="Unknown medium 'mead'"):
        compile_scenario(_wine_scenario(medium="mead"))


# -- the seam connects to the runtime -----------------------------------------


def test_compiled_wine_scenario_ferments_and_conserves_carbon():
    # The wired kinetics make a compiled wine scenario actually ferment — sugar
    # falls, ethanol and CO2 rise — and the run conserves carbon end-to-end,
    # proving compile() feeds the *full* Process+modifier set into simulate()
    # cleanly. (The dryness *timing* is the §2.2 benchmark's job and needs tuning;
    # here we only assert direction + conservation, so this test does not move with
    # the tuning.)
    compiled = compile_scenario(_wine_scenario(duration_days=21.0), strict=True)
    traj = simulate(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
    )
    assert traj.success
    sugar = traj.series("S")  # wine: 1-D (single slot)
    assert float(sugar[-1]) < float(sugar[0])  # sugar consumed
    assert traj.series("E")[-1] > 0.0  # ethanol produced
    assert traj.series("CO2")[-1] > 0.0  # CO2 evolved
    assert_conserved(
        traj,
        total_carbon(
            compiled.schema, biomass_carbon_fraction=compiled.parameters.value("biomass_C_fraction")
        ),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )
    assert_nonnegative(traj, ("X", "S", "N", "E", "CO2"), atol=1e-7)


def test_compiled_beer_scenario_ferments_and_conserves_carbon():
    # Same end-to-end check on beer's 3-sugar schema: the shared kinetic set
    # consumes the wort sugars (sequential uptake lives inside the uptake Process)
    # and closes carbon across all three slots.
    compiled = compile_scenario(_beer_scenario(duration_days=14.0), strict=True)
    traj = simulate(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
    )
    assert traj.success
    sugar = traj.series("S")  # beer: 2-D (three slots) x time
    assert float(sugar[:, -1].sum()) < float(sugar[:, 0].sum())  # total wort sugar consumed
    assert traj.series("E")[-1] > 0.0
    assert traj.series("CO2")[-1] > 0.0
    assert_conserved(
        traj,
        total_carbon(
            compiled.schema, biomass_carbon_fraction=compiled.parameters.value("biomass_C_fraction")
        ),
        rtol=1e-5,
        atol=1e-6,
        label="carbon",
    )


# -- the compile vocabulary stays in sync with the MEDIA registry -------------


@pytest.mark.parametrize("medium", sorted(MEDIA))
def test_every_registered_medium_compiles(medium):
    # Guards against a medium added to MEDIA without a matching allowed-keys +
    # initial-builder entry in compile.py: a new medium must be representable
    # here, and must then compile, or this fails at test time (not runtime).
    representative = {
        "wine": (_wine_scenario(), None),
        "beer": (_beer_scenario(), [WINE_PARAMS]),
    }
    assert medium in representative, (
        f"MEDIA registers {medium!r} but it has no representative scenario here; "
        "add one (and its compile vocabulary) so the seam stays covered."
    )
    scenario, paths = representative[medium]
    compiled = compile_scenario(scenario, parameter_paths=paths)
    assert compiled.y0.shape == (compiled.schema.size,)
