"""D-142: the SO2:O2 molar reaction ratio — Danilewicz's structural series, measured against data.

Sources, and which one backs which assertion:

* Danilewicz, J.C. (2016). "Reaction of Oxygen and Sulfite in Wine." Am. J. Enol. Vitic.
  67(1):13-17. doi:10.5344/ajev.2015.15069. **PAYWALLED — abstract only.** It backs the
  STRUCTURE and nothing else: "the O2:SO2 molar reaction ratio is close to 1:2 in ideal
  experimental conditions"; across eight wines (three red, five white) "the reaction ratio was
  found to be decreased down to 1:1.7 in most wines"; and with a white wine dosed with
  benzenesulfinic acid — "this substance reacts very efficiently with quinones and would
  therefore prevent their interaction with sulfite" — "the molar reaction ratio was then reduced
  to 1:1, as has been previously observed in model wine". No doses, temperatures or per-wine
  values are available, so none are asserted here.

* Miao, Y. & Waterhouse, A.L. (2025). "Rapid White Wine Shelf-Life Prediction by Forecasting
  Free SO2 Loss Post-Bottling." Am. J. Enol. Vitic. 76:0760008. doi:10.5344/ajev.2025.24057.
  **OPEN ACCESS (CC BY) — full text read.** This supplies the NUMBERS, the protocol this file
  reproduces, and a transcription of Danilewicz's own per-wine spread. Its Table 2 gives the
  SO2:O2 molar reaction ratio for eight commercial white wines; the seven with normal SO2 span
  **1.0972-1.6621**, and it reports Danilewicz 2016's commercial-white-wine range as
  **1.29-2.03**. Protocol: air-saturated to 7 mg/L O2, sealed BOD bottles with NO headspace,
  45 C for 5 days in the dark, free + total SO2 on days 0/1/3/5, across two rounds of air
  saturation; the ratio is the slope of total SO2 against total O2 consumed, times M(O2)/M(SO2).

**Why 1.7 is NOT the assertion.** D-138 and D-141 both say "1.7 must emerge, never be fitted",
and this file honours that by never pinning it: 1.7 is Danilewicz's mode, but it sits at the TOP
of Miao's independently measured range (whose mean over the seven normal wines is ~1.372). A
benchmark pinned at 1.7 would be pinned to one dataset. What is asserted instead is the
structure both datasets agree on — the two limits, the strict ordering between them, and the
real-wine band read as the union of the two reported ranges.

**The falsifier, and what it actually is.** D-141 claimed the direct set "structurally cannot
produce" this series. Measured, that claim needs a sharper statement: the direct set reaches a
real-wine *number* inside the envelope perfectly well (1.7707), because its fixed
``_SO2_PER_O2 = 2`` times bisulfite's share of a six-way O2 competition can land anywhere in
(0, 2). What it cannot do is *move between the limits by blocking a route*, because it has no
route to block. Blocking the quinone route changes the direct set's ratio by **0.0003** and the
cascade's by **0.859**. **That traverse — not the value 1.7 — is the rebuild's falsifiable
content**, and it is the one assertion here that no choice of operating point can soften.

**The exactness of the two limits is emergent, not fitted.** Both asymptotes fall out of
stoichiometry constants that are all 1 (``_H2O2_PER_O2``, ``_QUINONE_PER_O2``, ``_SO2_PER_H2O2``,
``_SO2_PER_QUINONE``): one O2 makes two oxidising equivalents, and bisulfite takes one of them or
both. Both limits are approached **from below** at finite SO2 for a physical reason, not a
numerical one — bisulfite shares the H2O2 node with ethanol, so its share of that node is <1
until SO2 is large. That is the same qualification Danilewicz makes when he calls the model-wine
quinone reduction "near quantitative" rather than quantitative.

**THE OPERATING POINT IS THE WHOLE BALLGAME, and this file got it wrong once.** Miao excludes
his own wine #1 because its free SO2 ran 8 -> 3 mg/L. The first version of this file stated that
bound in prose, never enforced it, and ran at a dose where free SO2 ended at **2.8 mg/L**
(direct) and 6.8 (cascade) — squarely inside the excluded regime. Every ratio it reported was
reading the curvature Miao removes, and the conclusion it drew was the opposite of the truth: it
shipped an ``xfail(strict=True)`` claiming the cascade fell BELOW the band by 4-7x on quinone
capture. Re-measured where the exclusion criterion actually holds, the cascade reads **1.1339 —
inside Miao's own 1.0972-1.6621** — and the calibrated direct set reads **1.7707, ABOVE it**.
``test_the_operating_point_satisfies_miaos_exclusion_criterion`` now enforces the bound in code.

**The direct half of that inversion is robust; the cascade half is dose-contingent.** The direct
set is above Miao at every valid dose and immovable. The cascade **straddles Miao's floor**:
1.0704 at 60 (out), 1.1035 at 70, 1.1339 at 80 (in) — and 80 is the shared operating point only
because the *direct* set exhausts its free SO2 at 60. So the cascade's in-band result is a
statement about a dose the other alternative forced, with ~3.4% headroom. Asserted as the
straddle it is (``test_cascade_brackets_miaos_lower_bound_across_its_valid_dose_range``), not as
a single in-band number — which is the mistake this file already made once, in the other
direction.

**But the agreement is at the wrong operating point, and that is the real finding.** The sim
under-binds SO2: its buffering capacity is below Miao's on both alternatives, so it has no bound
reservoir to replenish free SO2 under the O2 challenge. To stay above his floor it must be dosed
to ~59 mg/L free — roughly 2x his highest wine. **No dose matches his operating point at both
ends of the window.** So the in-band result is an envelope check at the wrong sulfite level, not
a validation, and the quinone-branching question D-141 left open is NOT settled by it. That
constraint is asserted, not merely conceded, by
``test_the_sim_cannot_actually_reach_miaos_operating_point``.
"""

from typing import NamedTuple

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario import Intervention, Scenario, TemperaturePoint
from fermentation.scenario.compile import compile_scenario

pytestmark = pytest.mark.benchmark

M_O2, M_SO2 = 32.00, 64.06

FERM, SETTLE, ROUND, N_ROUNDS = 20.0, 30.0, 5.0, 2
O2_DOSE, TEST_C = 7.0, 45.0

#: Miao Table 2, the seven wines with normal SO2 (wine #1 is excluded there too — it started at
#: 8 mg/L free SO2 and fell to 3, and "when SO2 approaches very low levels, oxidation reactions
#: must necessarily involve other components, and thus the ratio of reaction with O2 must
#: change", which is why it reads 0.8443).
MIAO_BAND = (1.0972, 1.6621)
#: Danilewicz 2016's commercial-white-wine spread, as transcribed by Miao & Waterhouse 2025.
DANILEWICZ_BAND = (1.29, 2.03)
#: The assertion band: the union of the two independently reported ranges. Deliberately NOT the
#: intersection and NOT a pin at 1.7 — the two datasets disagree, and the honest target is the
#: envelope they jointly support.
OBSERVED_BAND = (min(MIAO_BAND[0], DANILEWICZ_BAND[0]), max(MIAO_BAND[1], DANILEWICZ_BAND[1]))

#: Miao's sulfite buffering capacity (Δtotal SO2 / Δfree SO2) across the same eight wines.
MIAO_BUFFERING_BAND = (1.2526, 1.9882)
#: Miao excludes his wine #1 because its free SO2 ran 8 -> 3 mg/L: "when SO2 approaches very low
#: levels, oxidation reactions must necessarily involve other components, and thus the ratio of
#: reaction with O2 must change". Every ratio compared to his band must hold above this.
FREE_SO2_FLOOR_MGL = 10.0

#: The shared operating point: the lowest dose ON THE TESTED GRID {40, 60, 70, 80, 100, 120} at
#: which free SO2 stays above :data:`FREE_SO2_FLOOR_MGL` through the whole window on BOTH
#: alternatives — enforced by ``test_the_operating_point_satisfies_miaos_exclusion_criterion``
#: rather than asserted in prose. It was 40.0 when this file first shipped, which was a DEFECT:
#: free SO2 ended at 2.8 mg/L (direct) and 6.8 (cascade), i.e. inside the regime Miao's own band
#: excludes, and every ratio measured there was reading the curvature he removes.
#:
#: **It is forced by the DIRECT set, not the cascade.** The cascade alone is valid from 60
#: (free SO2 ends at 17.4 mg/L); the direct set is not (8.2 mg/L). So any statement about the
#: cascade *at this dose* is contingent on a dose the other alternative required — which is why
#: ``test_cascade_brackets_miaos_lower_bound`` measures its whole valid range instead.
WINE_REALISTIC_SO2 = 80.0
#: Doses at which the CASCADE alone satisfies the exclusion criterion, lowest first.
CASCADE_VALID_DOSES = (60.0, 70.0, 80.0)
#: Far above any legal wine. Present ONLY to show the two limits are asymptotes at exactly 1 and
#: 2 rather than approximate coincidences.
LIMIT_SO2 = 500.0

#: The quinone fates that COMPETE with bisulfite for the quinone pool. Zeroing them is the
#: "bisulfite wins both oxidising equivalents" limit that must read 1:2.
_COMPETING_QUINONE_CONSTANTS = (
    "k_quinone_polymerization",
    "k_strecker",
    "k_anthocyanin_fade",
    "k_ellagitannin_oxidation",
)
#: Danilewicz's benzenesulfinic acid, whose analogue here is an overwhelming competing quinone
#: sink. ``quinone_sulfonation`` cannot simply be switched off instead: it reads
#: ``k_so2_oxidation``, the SAME constant the H2O2 node uses as its bisulfite weight, so zeroing
#: that would block BOTH oxidising equivalents and drive the ratio to 0 rather than to 1.
_BLOCK_FACTOR = 1.0e6


def _scenario(so2_mgl: float) -> Scenario:
    """A finished white wine, then Miao's accelerated-oxidation protocol.

    ``closure="hermetic"`` is the sim's named zero-ingress closure (the D-45 "a named zero reads
    as a choice" idiom), which is the analogue of Miao's sealed BOD bottle — so "O2 consumed" is
    unambiguous rather than entangled with closure permeation.
    """
    t0 = FERM + SETTLE
    interventions = [
        Intervention(day=FERM - 1.0, action="add_so2", params={"so2_mgl": so2_mgl}),
        Intervention(day=FERM, action="begin_aging"),
    ]
    for r in range(N_ROUNDS):
        interventions.append(
            Intervention(day=t0 + r * ROUND, action="add_oxygen", params={"o2_mgl": O2_DOSE})
        )
    return Scenario(
        name="danilewicz-so2-o2",
        medium="wine",
        strain="generic",
        initial={
            "brix": 21.0,
            "yan_mgl": 200.0,
            "pitch_gpl": 0.3,
            "tannin_gpl": 0.3,  # white wine: low tannin
            "anthocyanin_gpl": 0.0,  # white wine: no anthocyanin
            "amino_acids_gpl": 0.5,
        },
        temperature_schedule=[
            TemperaturePoint(day=0.0, celsius=20.0),
            TemperaturePoint(day=t0 - 0.01, celsius=20.0),
            TemperaturePoint(day=t0, celsius=TEST_C),
        ],
        closure="hermetic",
        duration_days=t0 + N_ROUNDS * ROUND,
        interventions=interventions,
    )


class _Run(NamedTuple):
    """One accelerated-oxidation run, reduced to the three things Miao reports plus the guard."""

    #: Miao Table 2 — d(total SO2) / d(O2 consumed), molar.
    ratio: float
    #: Miao Table 3 — d(total SO2) / d(free SO2). The sulfite buffering capacity.
    buffering: float
    #: Free SO2 [mg/L] at each of the eight sample points, for the exclusion-criterion guard.
    free_mgl: tuple[float, ...]


def _run(so2_mgl: float, oxidative: str, mode: str) -> _Run:
    """Miao's SO2:O2 molar reaction ratio, by his own regression method.

    ``mode`` is ``"real"``, ``"blocked"`` (quinone route trapped, Danilewicz's benzenesulfinic
    acid) or ``"ideal"`` (bisulfite wins both oxidising equivalents).
    """
    compiled = compile_scenario(_scenario(so2_mgl), oxidative=oxidative)
    # NOTE ``CompiledScenario.param_values`` is a PROPERTY returning a fresh dict on each access,
    # so mutating it in place edits a throwaway and silently changes nothing. Bind it once.
    params = dict(compiled.param_values)
    if mode == "blocked" and "k_quinone_polymerization" in params:
        params["k_quinone_polymerization"] *= _BLOCK_FACTOR
    elif mode == "ideal":
        for name in _COMPETING_QUINONE_CONSTANTS:
            if name in params:
                params[name] = 0.0
    # The routes are blocked by their PARAMETERS rather than by ``ProcessSet.disable``, because
    # the ``begin_aging`` event re-enables every aging Process mid-run and silently undoes it.

    t0h = (FERM + SETTLE) * 24.0
    t_end = (t0h / 24.0 + N_ROUNDS * ROUND) * 24.0
    traj = simulate_scheduled(
        compiled.process_set,
        params,
        compiled.y0.copy(),
        (0.0, t_end),
        events=compiled.events,
        t_eval=np.linspace(0.0, t_end, 6000),
    )

    o2 = traj.y[compiled.schema.slice("o2")][0]
    total_so2 = traj.y[compiled.schema.slice("so2_total")][0]

    added = np.zeros_like(traj.t)
    for r in range(N_ROUNDS):
        added += np.where(traj.t >= t0h + r * ROUND * 24.0 - 1e-9, O2_DOSE / 1000.0, 0.0)
    o2_at_start = float(np.interp(t0h - 1e-6, traj.t, o2))
    consumed = o2_at_start + added - o2  # hermetic closure => no ingress term

    hours = [t0h + (r * ROUND + d) * 24.0 for r in range(N_ROUNDS) for d in (0.0, 1.0, 3.0, 5.0)]
    x = np.array([float(np.interp(h, traj.t, consumed)) for h in hours]) / M_O2
    y = np.array([float(np.interp(h, traj.t, total_so2)) for h in hours]) / M_SO2
    # Free SO2 is speciated at the eight SAMPLE points only, not mapped over all 6000 stored
    # columns: the coupled four-carbonyl binding solve is expensive and this file only ever reads
    # it at Miao's sampling schedule (mapping it wholesale tripled the file's runtime).
    # ``free_so2`` is g/L; keep that array for the mg/L guard and divide only for the molar
    # regression. (Scaling the mol/L array by 1000 gives mmol/L, not mg/L — that slip made the
    # exclusion guard read 0.3 where the truth was 17.8.)
    f_gpl = np.array(
        [
            acidbase.free_so2(
                traj.y[:, int(np.argmin(np.abs(traj.t - h)))], compiled.schema, params
            )
            for h in hours
        ]
    )
    f = f_gpl / M_SO2
    ratio = -float(np.polyfit(x, y, 1)[0])  # SO2 falls as O2 is consumed; report it positive
    free_factor = -float(np.polyfit(x, f, 1)[0])
    return _Run(
        ratio=ratio,
        buffering=ratio / free_factor if free_factor else float("nan"),
        free_mgl=tuple(float(v) * 1000.0 for v in f_gpl),
    )


@pytest.fixture(scope="module")
def runs() -> dict[tuple[str, str, float], _Run]:
    """Every run this file needs, integrated once each."""
    wanted = [
        ("cascade", "real", WINE_REALISTIC_SO2),
        ("cascade", "blocked", WINE_REALISTIC_SO2),
        ("cascade", "ideal", WINE_REALISTIC_SO2),
        ("cascade", "blocked", LIMIT_SO2),
        ("cascade", "ideal", LIMIT_SO2),
        ("direct", "real", WINE_REALISTIC_SO2),
        ("direct", "blocked", WINE_REALISTIC_SO2),
        ("direct", "ideal", WINE_REALISTIC_SO2),
        *[("cascade", "real", d) for d in CASCADE_VALID_DOSES if d != WINE_REALISTIC_SO2],
        ("direct", "real", CASCADE_VALID_DOSES[0]),
    ]
    return {key: _run(key[2], key[0], key[1]) for key in wanted}


# -- the operating point must be one the reference band actually covers ------------------------


def test_the_operating_point_satisfies_miaos_exclusion_criterion(runs):
    """Free SO2 stays above 10 mg/L at every sample point, on BOTH alternatives.

    This enforces in code the bound the first version of this file only stated in prose, and it
    is the guard that would have caught that defect: at the original 40 mg/L dose free SO2 ended
    at 2.8 mg/L (direct) and 6.8 (cascade), inside the very regime Miao excludes his wine #1 for.
    Every band comparison below is meaningless without it.
    """
    for alternative in ("direct", "cascade"):
        free = runs[(alternative, "real", WINE_REALISTIC_SO2)].free_mgl
        assert min(free) > FREE_SO2_FLOOR_MGL, f"{alternative} fell to {min(free):.1f} mg/L"


def test_the_sim_cannot_actually_reach_miaos_operating_point(runs):
    """Characterizes the gap that makes every band comparison here CONDITIONAL, not clean.

    Miao's wines carry a large bound-SO2 reservoir (total/free up to 153/23 = 6.6) that keeps
    free SO2 topped up under a 2 x 7 mg/L O2 challenge at free SO2 of only 9-28 mg/L. The sim
    under-binds SO2 (a limitation ``analysis.bound_so2_series`` already concedes in its own
    docstring, and which D-51's three extra carbonyls did not close), so its buffering capacity
    lands BELOW Miao's range on both alternatives. The consequence is structural: to keep free
    SO2 above his floor the sim must be dosed to ~59 mg/L free, roughly 2x his highest wine.
    **There is no dose that matches his operating point at both ends of the window.**

    Asserted as a regression guard on the CURRENT gap, in the idiom of
    ``test_validation_varela2004.py``: do not widen it to make CI green — close it by modelling
    the missing binders, and update D-142 when the direction changes.
    """
    for alternative in ("direct", "cascade"):
        buffering = runs[(alternative, "real", WINE_REALISTIC_SO2)].buffering
        assert buffering < MIAO_BUFFERING_BAND[0], (
            f"{alternative} buffering {buffering:.3f} has risen into Miao's "
            f"{MIAO_BUFFERING_BAND} — the under-binding gap has closed; update D-142"
        )


# -- the two limits: emergent stoichiometry, not fitted values ---------------------------------


def test_blocked_quinone_route_gives_one_so2_per_o2(runs):
    """Danilewicz's 1:1 — quinone trapped, so only H2O2 can oxidise bisulfite.

    Asserted at the asymptote, because at finite SO2 the limit is approached from BELOW:
    bisulfite shares the H2O2 node with ethanol, so its share of that node is <1.
    """
    assert runs[("cascade", "blocked", LIMIT_SO2)].ratio == pytest.approx(1.0, abs=0.01)


def test_uncontested_quinone_gives_two_so2_per_o2(runs):
    """Danilewicz's 1:2 — bisulfite takes BOTH oxidising equivalents one O2 makes.

    This is the number that must emerge rather than be asserted, and it does: every
    stoichiometry constant on the path is 1, so 2 is a consequence of "one O2 yields one H2O2
    and one quinone", never a fitted yield.
    """
    assert runs[("cascade", "ideal", LIMIT_SO2)].ratio == pytest.approx(2.0, abs=0.02)


def test_the_series_is_strictly_ordered_at_wine_realistic_sulfite(runs):
    """1:1 < real wine < 1:2, with real strictly interior — the whole point of partial capture."""
    blocked = runs[("cascade", "blocked", WINE_REALISTIC_SO2)].ratio
    real = runs[("cascade", "real", WINE_REALISTIC_SO2)].ratio
    ideal = runs[("cascade", "ideal", WINE_REALISTIC_SO2)].ratio
    assert blocked < real < ideal
    # and the interior point is not merely nudged off either limit
    assert real - blocked > 0.02
    assert ideal - real > 0.02


# -- the falsifier: only one of the two alternatives can traverse the series --------------------


def test_direct_set_cannot_move_between_the_limits(runs):
    """The rebuild's falsifiable content, stated as the thing the direct set CANNOT do.

    Blocking or freeing the quinone route leaves the direct set's ratio unchanged to <0.01,
    because it has no quinone to route. The residual is not exactly zero only because the
    competing constants zeroed for the "ideal" arm are themselves direct O2 sinks there, which
    perturbs the six-way O2 competition very slightly.
    """
    blocked = runs[("direct", "blocked", WINE_REALISTIC_SO2)].ratio
    ideal = runs[("direct", "ideal", WINE_REALISTIC_SO2)].ratio
    assert abs(ideal - blocked) < 0.01


def test_cascade_traverses_the_series(runs):
    """The same manipulation moves the cascade by three orders of magnitude more.

    Guarded as a ratio against the direct set's residual rather than as a bare threshold, so it
    keeps meaning if either alternative's absolute numbers move.
    """
    cascade_span = (
        runs[("cascade", "ideal", WINE_REALISTIC_SO2)].ratio
        - runs[("cascade", "blocked", WINE_REALISTIC_SO2)].ratio
    )
    direct_span = abs(
        runs[("direct", "ideal", WINE_REALISTIC_SO2)].ratio
        - runs[("direct", "blocked", WINE_REALISTIC_SO2)].ratio
    )
    assert cascade_span > 0.4
    assert cascade_span > 50.0 * direct_span


# -- agreement with the measured band ----------------------------------------------------------


def test_cascade_brackets_miaos_lower_bound_across_its_valid_dose_range(runs):
    """The cascade STRADDLES Miao's floor: out of band at its lowest valid dose, in above it.

    Measured: 1.0704 at 60 (below 1.0972), 1.1035 at 70, 1.1339 at 80 — so "the cascade lands in
    Miao's band" is true only from ~70 upward, and the shared operating point that makes it true
    was forced by the *direct* set's exhaustion, not by the cascade.

    Asserting the straddle rather than one in-band number is deliberate. The band's floor and the
    cascade's value differ by ~3.4% at the shared operating point, on a coarse dose grid; a bare
    ``lo <= ratio`` there would be a 3%-headroom assertion masquerading as agreement, and any
    unrelated model change would flip it into a confusing red. This states what was actually
    measured, and it is the second time in this file's short life that a single-dose verdict on
    this quantity turned out to be an artefact of the dose.
    """
    ratios = [runs[("cascade", "real", d)].ratio for d in CASCADE_VALID_DOSES]
    assert ratios == sorted(ratios), "ratio should rise monotonically with sulfite"
    assert ratios[0] < MIAO_BAND[0], "lowest valid dose should sit BELOW Miao's floor"
    assert MIAO_BAND[0] <= ratios[-1] <= MIAO_BAND[1], "shared operating point should be in band"


def test_the_shared_operating_point_is_forced_by_the_direct_set(runs):
    """Why the comparison is made at 80 and not at the cascade's own lowest valid dose.

    At 60 the cascade clears Miao's floor and the direct set does not, so 60 cannot be used for a
    like-for-like comparison. This records the asymmetry rather than leaving it as a bare
    constant, because it is exactly what makes the cascade's in-band result contingent.
    """
    low = CASCADE_VALID_DOSES[0]
    assert min(runs[("cascade", "real", low)].free_mgl) > FREE_SO2_FLOOR_MGL
    assert min(runs[("direct", "real", low)].free_mgl) < FREE_SO2_FLOOR_MGL


def test_direct_set_overshoots_miaos_range_but_stays_in_the_union_envelope(runs):
    """The SHIPPED default reads 1.7707 — above Miao's 1.6621, inside Danilewicz's 1.29-2.03.

    So at a valid operating point the two alternatives straddle the reference: the cascade sits
    inside Miao's own range and the calibrated direct set sits above it, in the part of the
    envelope only Danilewicz's dataset covers. That is a much weaker result for the direct set
    than this file first recorded, and it is not evidence for its mechanism either way — it
    reaches its number through a fixed ``_SO2_PER_O2 = 2`` scaled by bisulfite's share of a
    six-way competition for O2, the structure Gate 1 (D-137) found to be mechanistically false.
    """
    lo, hi = OBSERVED_BAND
    ratio = runs[("direct", "real", WINE_REALISTIC_SO2)].ratio
    assert lo <= ratio <= hi
    assert ratio > MIAO_BAND[1]
