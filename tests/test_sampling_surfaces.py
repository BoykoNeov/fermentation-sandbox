"""The sampling-surface split, pinned (decision D-156, closing D-153 leg 1).

D-153 enumerated the archive's 18 parameter files and found they do **not** all reach a
sampler, and those that do are **not** all drawn the same way:

===================  ==========================================  =================
surface              files                                       distribution
===================  ==========================================  =================
compile seam         ``<medium>_generic`` + the 14 shared files   **triangular**
sensory              ``psychophysics.yaml``                       **uniform**
sensory              ``sensory.yaml``                             **never drawn**
===================  ==========================================  =================

This split is load-bearing for *records*, not only for runs. Every audit that quotes a
band's mass does arithmetic on one of these distributions and gets a different answer
from the other: at the copper band's own geometry the share of draws above 918 L/g reads
0.2895 triangular and 0.4477 uniform (D-153), a 16-point gap in the anti-conservative
direction. D-152's flag and D-154's band repair are both built on the triangular figure.

**Why these tests exist, measured rather than assumed.** D-156 ran a four-arm mutation
matrix against the shipped tree. Every arm was **GREEN** — the full suite stayed at 1444
passed while:

* ``sensory.yaml`` was added to ``_load_parameters``' shared list (arm A),
* ``_axis_draws`` was switched from ``rng.uniform`` to ``rng.triangular`` (arm B),
* ``sample_parameters``' default distribution was flipped to ``uniform`` (arm C1),
* ``simulate_ensemble``' default distribution was flipped to ``uniform`` (arm C2).

The four tests below are those four arms, closed. Each is written to fail on its own arm
and is checked against it, so none is a decoration — the D-155 standard.

**What is deliberately NOT guarded here: the removal direction.** Dropping any one of the
14 shared files from the seam is caught loudly by consumption — a compiled scenario reads
a name the file alone supplies and raises. D-156 measured this on all 14 rather than
sampling one: every drop is RED, from 8 failures (``hops``) to 276 (three of the
volatile-precursor files). Re-pinning it here would be the redundant half of the guard.
The addition direction has no such backstop, which is why it is the direction pinned.
"""

import numpy as np
import pytest

from fermentation.core.process import ProcessSet
from fermentation.parameters import default_data_dir
from fermentation.parameters.store import ParameterSet, load_parameters
from fermentation.runtime import sample_parameters, simulate_ensemble
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.sensory.compression import _axis_draws, _exponent_key, load_exponents
from fermentation.sensory.oav import load_thresholds

DATA = default_data_dir()

#: The 14 medium-agnostic files ``scenario.compile._load_parameters`` merges alongside the
#: medium file. Restated here ON PURPOSE: this list IS the pinned quantity, so deriving it
#: from the module under test would make the assertion vacuous (the D-108/D-109 shape). The
#: band CONTENTS are still read from the YAML, never restated — the D-100 lesson.
SHARED_FILES = (
    "acidbase",
    "vicinal_diketones",
    "acetaldehyde",
    "keto_acids",
    "hydrogen_sulfide",
    "additions",
    "hops",
    "aging",
    "oak",
    "polymerization",
    "thermal",
    "dms",
    "bound_sulfides",
    "closure",
)

#: The two files that load standalone and must stay off the seam (D-24, D-98).
SENSORY_FILES = ("psychophysics", "sensory")


def _wine_scenario() -> Scenario:
    return Scenario(
        name="surface-split-wine",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=14.0,
    )


def _beer_scenario() -> Scenario:
    return Scenario(
        name="surface-split-beer",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 70.0,
            "maltotriose_gpl": 20.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 1.0,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=7.0,
    )


def _names(*stems: str) -> set[str]:
    return {n for stem in stems for n in load_parameters(DATA / f"{stem}.yaml").names}


# -- surface 1: which files reach the compile seam (mutation arm A) -----------


@pytest.mark.parametrize(
    ("scenario", "medium_file"),
    [(_wine_scenario, "wine_generic"), (_beer_scenario, "beer_generic")],
)
def test_the_compile_seam_merges_exactly_the_declared_parameter_files(scenario, medium_file):
    """The seam's parameter set is the medium file plus the 14 shared files — no more.

    Asserted by NAME CONTENT rather than by reading ``_load_parameters``' source, so it
    fires in both directions: a new YAML joining the shared list adds names the union does
    not have, and a file leaving it removes names the union does. A run sees exactly ONE
    medium file, so the other medium's names must be absent too — that is what makes
    wine 98 + 149 shared / beer 32 + 149 the ceilings D-153 reports rather than 351.
    """
    compiled = compile_scenario(scenario())
    assert set(compiled.parameters.names) == _names(medium_file, *SHARED_FILES)


# -- surface 2 + 3: the sensory files are unreachable from the sampler (arm A) -


def test_the_sensory_files_cannot_be_reached_by_the_ensemble_sampler():
    """``psychophysics.yaml`` and ``sensory.yaml`` are invisible to ``simulate_ensemble``.

    The mechanism is name intersection, not intent: ``_resolve_sample_names`` computes
    ``chosen &= set(parameters.names)``, so a name absent from the compiled seam cannot be
    drawn even when a caller passes it explicitly via ``only=``. Disjointness is therefore
    the whole property, and it is what makes ``dominant_flip_sensitivity`` a *manual* Monte
    Carlo (D-24) and ``sensory.yaml``'s 36 bands genuinely never-drawn (D-153).
    """
    seam = set(compile_scenario(_wine_scenario()).parameters.names)
    seam |= set(compile_scenario(_beer_scenario()).parameters.names)
    assert seam.isdisjoint(load_exponents().names)
    assert seam.isdisjoint(load_thresholds().names)
    # ...and the same statement read off the files, so a rename cannot quietly satisfy it.
    assert seam.isdisjoint(_names(*SENSORY_FILES))


# -- the two distributions ----------------------------------------------------
#
# Tail mass above a point x is the statistic every band audit in this archive quotes, so
# it is the statistic these tests discriminate on. It is NOT usable at x == the mode:
# both laws put (mode-lo)/(hi-lo) of their mass below the mode, so the two agree there
# exactly. Away from the mode they diverge, and that is where these tests sample.


def _uniform_tail(lo: float, hi: float, x: float) -> float:
    return (hi - x) / (hi - lo)


def _triangular_tail(lo: float, mode: float, hi: float, x: float) -> float:
    if x >= mode:
        return (hi - x) ** 2 / ((hi - lo) * (hi - mode))
    return 1.0 - (x - lo) ** 2 / ((hi - lo) * (mode - lo))


def test_the_psychophysics_axis_is_drawn_uniform_not_triangular():
    """``_axis_draws`` draws UNIFORM over each exponent band (D-98, D-153 leg 1).

    Deliberate and defended in the docstring: the bands are honest ignorance, not panel
    spreads, so a uniform draw declines to claim the centre is likelier than the edge. The
    consequence for an audit is arithmetic — an excluded exponent's share of the draws is
    the WIDTH share, not the triangular mass — which is why it is pinned rather than left
    to the prose.

    Measured on a single-pool axis, which is the unconditioned marginal: with one pool
    ``preserve_order`` cannot bind (``width > 1`` is false), so this isolates the base draw
    from the rank conditioning D-98 layers on top of it. The conditioned law is that base
    conditioned on Cain's rank, so the base is the thing to pin.
    """
    exponents = load_exponents()
    band = exponents[_exponent_key("acetaldehyde")].uncertainty
    mode = exponents.value(_exponent_key("acetaldehyde"))
    lo, hi = band.low, band.high
    x = 0.7  # off the mode (0.65), where the two laws separate
    draws = _axis_draws(("acetaldehyde",), exponents, np.random.default_rng(0), 40_000, True)

    assert lo <= draws.min() and draws.max() <= hi
    tail = float(np.mean(draws[:, 0] > x))
    # Uniform 0.2439 vs triangular 0.1626 at this band — an 8-point separation, ~30 sigma
    # at 40k draws, so the tolerance below is noise-sized rather than gap-sized.
    assert tail == pytest.approx(_uniform_tail(lo, hi, x), abs=0.01)
    assert abs(tail - _triangular_tail(lo, mode, hi, x)) > 0.05


def test_the_compile_seam_is_drawn_triangular_by_default():
    """``sample_parameters`` draws ``triangular(low, value, high)`` unless told otherwise.

    This default is what D-152's 29 % and D-154's band repair are arithmetic ON, so a
    silent flip to uniform would not break a run — it would retroactively falsify records.
    Pinned behaviourally (the drawn mass), not by reading the signature's default, so the
    assertion is about what the sampler does.

    Measured on ``k_copper_multiplier``'s own shipped band, read from ``aging.yaml`` rather
    than restated, because that is the band the two records argue about.
    """
    parameters = load_parameters(DATA / "aging.yaml")
    band = parameters["k_copper_multiplier"].uncertainty
    mode = parameters.value("k_copper_multiplier")
    lo, hi = band.low, band.high
    x = 500.0  # below the mode, where the laws separate by 16 points
    rng = np.random.default_rng(0)
    drawn = np.array(
        [
            sample_parameters(parameters, rng, names=["k_copper_multiplier"])["k_copper_multiplier"]
            for _ in range(20_000)
        ]
    )

    tail = float(np.mean(drawn > x))
    assert tail == pytest.approx(_triangular_tail(lo, mode, hi, x), abs=0.02)
    assert abs(tail - _uniform_tail(lo, hi, x)) > 0.1


def test_the_ensemble_wrapper_inherits_that_default_rather_than_its_own(toy_schema, toy_process):
    """``simulate_ensemble`` carries a SECOND copy of the default, and it is pinned too.

    Arm C1 and arm C2 flip different lines: ``sample_parameters``' default is what a direct
    caller gets, ``simulate_ensemble``'s is what every compile-seam ensemble gets, and the
    test above catches only the first. ``simulate_ensemble`` threads one local into both
    ``_build_samples(distribution=...)`` and ``Ensemble(distribution=...)``, so the reported
    field IS the law the members were drawn from — which is what makes reading it here an
    assertion about the draw rather than about a signature.

    The mass check is not repeated per member: 20 000 members would be 20 000
    ``solve_ivp`` integrations for a statistic the shared draw path already pins above.
    """
    process_set = ProcessSet(toy_schema, [toy_process], strict=True)
    y0 = toy_schema.pack({"S": 100.0, "E": 0.0, "CO2": 0.0})
    ensemble = simulate_ensemble(process_set, ParameterSet({}), y0, (0.0, 1.0), n_members=2, seed=0)
    assert ensemble.distribution == "triangular"
