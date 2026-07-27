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
produce" this series. Measured, that claim needs a sharper statement: the direct set reproduces
a real-wine *number* in band perfectly well (1.3875 at wine-realistic SO2 — see
``test_direct_set_reproduces_the_real_wine_number``), because its fixed ``_SO2_PER_O2 = 2`` times
bisulfite's share of a six-way O2 competition can land anywhere in (0, 2). What it cannot do is
*move between the limits by blocking a route*, because it has no route to block. Blocking the
quinone route changes the direct set's ratio by <0.001 and the cascade's by ~0.16. **That
traverse — not the value 1.7 — is the rebuild's falsifiable content.**

**The exactness of the two limits is emergent, not fitted.** Both asymptotes fall out of
stoichiometry constants that are all 1 (``_H2O2_PER_O2``, ``_QUINONE_PER_O2``, ``_SO2_PER_H2O2``,
``_SO2_PER_QUINONE``): one O2 makes two oxidising equivalents, and bisulfite takes one of them or
both. Both limits are approached **from below** at finite SO2 for a physical reason, not a
numerical one — bisulfite shares the H2O2 node with ethanol, so its share of that node is <1
until SO2 is large. That is the same qualification Danilewicz makes when he calls the model-wine
quinone reduction "near quantitative" rather than quantitative.

**What is red, and must stay red.** ``test_cascade_lands_in_the_observed_band`` is an
``xfail(strict=True)``: at wine-realistic SO2 the cascade returns 0.9888, BELOW Miao's
1.0972-1.6621. Decomposed (the identity ``r = H2O2 share to sulfite + quinone share to
sulfite``), the cascade's quinone share is 0.0584 where the band needs ~0.37-0.74. That is
D-141's "~100x fork" restated as a measurement: closing it needs roughly **4-7x** more quinone
reaching bisulfite, not 100x — a far smaller and better-posed gap than D-141 could state. It is
the same defect D-141 saw from the other side as A420 running 4.20x high in a sulfited wine.
Do NOT close this by fitting a constant; see D-142. When it is genuinely fixed the strict xfail
turns the suite red and forces this file to be updated.
"""

import numpy as np
import pytest

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

#: A dose whose free SO2 at test start (25.0 mg/L) sits inside Miao's own 9-28 mg/L spread.
WINE_REALISTIC_SO2 = 40.0
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


def _ratio(so2_mgl: float, oxidative: str, mode: str) -> float:
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
    slope, _ = np.polyfit(x, y, 1)
    return float(-slope)  # SO2 falls as O2 is consumed; report the positive ratio


@pytest.fixture(scope="module")
def ratios() -> dict[tuple[str, str, float], float]:
    """Every ratio this file needs, integrated once each."""
    wanted = [
        ("cascade", "real", WINE_REALISTIC_SO2),
        ("cascade", "blocked", WINE_REALISTIC_SO2),
        ("cascade", "ideal", WINE_REALISTIC_SO2),
        ("cascade", "blocked", LIMIT_SO2),
        ("cascade", "ideal", LIMIT_SO2),
        ("direct", "real", WINE_REALISTIC_SO2),
        ("direct", "blocked", WINE_REALISTIC_SO2),
        ("direct", "ideal", WINE_REALISTIC_SO2),
    ]
    return {key: _ratio(key[2], key[0], key[1]) for key in wanted}


# -- the two limits: emergent stoichiometry, not fitted values ---------------------------------


def test_blocked_quinone_route_gives_one_so2_per_o2(ratios):
    """Danilewicz's 1:1 — quinone trapped, so only H2O2 can oxidise bisulfite.

    Asserted at the asymptote, because at finite SO2 the limit is approached from BELOW:
    bisulfite shares the H2O2 node with ethanol, so its share of that node is <1.
    """
    assert ratios[("cascade", "blocked", LIMIT_SO2)] == pytest.approx(1.0, abs=0.01)


def test_uncontested_quinone_gives_two_so2_per_o2(ratios):
    """Danilewicz's 1:2 — bisulfite takes BOTH oxidising equivalents one O2 makes.

    This is the number that must emerge rather than be asserted, and it does: every
    stoichiometry constant on the path is 1, so 2 is a consequence of "one O2 yields one H2O2
    and one quinone", never a fitted yield.
    """
    assert ratios[("cascade", "ideal", LIMIT_SO2)] == pytest.approx(2.0, abs=0.02)


def test_the_series_is_strictly_ordered_at_wine_realistic_sulfite(ratios):
    """1:1 < real wine < 1:2, with real strictly interior — the whole point of "partial capture"."""
    blocked = ratios[("cascade", "blocked", WINE_REALISTIC_SO2)]
    real = ratios[("cascade", "real", WINE_REALISTIC_SO2)]
    ideal = ratios[("cascade", "ideal", WINE_REALISTIC_SO2)]
    assert blocked < real < ideal
    # and the interior point is not merely nudged off either limit
    assert real - blocked > 0.02
    assert ideal - real > 0.02


# -- the falsifier: only one of the two alternatives can traverse the series --------------------


def test_direct_set_cannot_move_between_the_limits(ratios):
    """The rebuild's falsifiable content, stated as the thing the direct set CANNOT do.

    Blocking or freeing the quinone route leaves the direct set's ratio unchanged to <0.01,
    because it has no quinone to route. The residual is not exactly zero only because the
    competing constants zeroed for the "ideal" arm are themselves direct O2 sinks there, which
    perturbs the six-way O2 competition very slightly.
    """
    blocked = ratios[("direct", "blocked", WINE_REALISTIC_SO2)]
    ideal = ratios[("direct", "ideal", WINE_REALISTIC_SO2)]
    assert abs(ideal - blocked) < 0.01


def test_cascade_traverses_the_series(ratios):
    """The same manipulation moves the cascade by more than an order of magnitude more.

    Guarded as a ratio against the direct set's residual rather than as a bare threshold, so it
    keeps meaning if either alternative's absolute numbers move.
    """
    cascade_span = (
        ratios[("cascade", "ideal", WINE_REALISTIC_SO2)]
        - ratios[("cascade", "blocked", WINE_REALISTIC_SO2)]
    )
    direct_span = abs(
        ratios[("direct", "ideal", WINE_REALISTIC_SO2)]
        - ratios[("direct", "blocked", WINE_REALISTIC_SO2)]
    )
    assert cascade_span > 0.4
    assert cascade_span > 50.0 * direct_span


# -- agreement with the measured band ----------------------------------------------------------


def test_direct_set_reproduces_the_real_wine_number(ratios):
    """The SHIPPED default lands inside the observed envelope, at wine-realistic SO2.

    This is a regression guard on the calibrated set, not a vindication of its mechanism: it
    reaches an in-band number through a fixed ``_SO2_PER_O2 = 2`` scaled by bisulfite's share of
    a six-way competition for O2, which is exactly the structure Gate 1 (D-137) found to be
    mechanistically false. Right number, wrong reason — and the test above is what says so.
    """
    lo, hi = OBSERVED_BAND
    assert lo <= ratios[("direct", "real", WINE_REALISTIC_SO2)] <= hi


@pytest.mark.xfail(
    strict=True,
    reason=(
        "D-142: the cascade returns 0.9888 at wine-realistic SO2, below Miao's 1.0972-1.6621. "
        "Its quinone-to-sulfite share is 0.058 where the band needs ~0.37-0.74, i.e. 4-7x too "
        "little quinone reaching bisulfite. This is D-141's open ~100x fork, now measured and "
        "much better posed. It must NOT be closed by fitting k_quinone_polymerization -- "
        "Nikolantonaki & Waterhouse 2012, the pull that would settle it, is paywalled and its "
        "abstract's rate constants are pseudo-first-order at unstated concentrations."
    ),
)
def test_cascade_lands_in_the_observed_band(ratios):
    """The gap that keeps the cascade non-default, asserted rather than described.

    Deliberately an ``xfail(strict=True)`` and not a widened band or a skip: when the quinone
    branching is genuinely fixed this turns the suite RED and forces the file to be updated,
    which a ``skip`` would not.
    """
    lo, hi = OBSERVED_BAND
    assert lo <= ratios[("cascade", "real", WINE_REALISTIC_SO2)] <= hi
