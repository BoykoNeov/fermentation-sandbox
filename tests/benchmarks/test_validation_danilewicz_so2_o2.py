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

**THE OPERATING POINT IS THE WHOLE BALLGAME, and this file got it wrong once.** The first version
ran at a dose where free SO2 ended at **2.8 mg/L** (direct) and 6.8 (cascade), deep in the
low-sulfite regime where — Miao's words — "oxidation reactions must necessarily involve other
components, and thus the ratio of reaction with O2 must change". Every ratio it reported was
reading that curvature, and the conclusion it drew was the opposite of the truth: it shipped an
``xfail(strict=True)`` claiming the cascade fell BELOW the band by 4-7x on quinone capture.
Re-measured clear of the curvature regime, the cascade reads **1.1339 — inside Miao's own
1.0972-1.6621** — and the calibrated direct set reads **1.7707, ABOVE it**.
``test_the_operating_point_clears_the_curvature_floor`` now enforces the bound in code.

**That floor is a sim-side conservatism and is NOT Miao's exclusion criterion**, which this file
also used to claim. His two printed criteria are a *starting* free SO2 level and an r2 cut, and
measured here neither does the job: the starting-level criterion would *select* the dose that
produced the inverted headline and reject the one that fixed it, and r2 does not discriminate at
all (0.98975 to 0.99987 across every dose). The floor is meanwhile stricter than his own data —
four of the seven wines defining his band would fail it. All three measurements, and why the
strictness is kept anyway, are recorded on :data:`SIM_CURVATURE_FLOOR_MGL`.

**The direct half of that inversion is robust; the cascade half is dose-contingent.** The direct
set is above Miao at every valid dose and immovable. The cascade **straddles Miao's floor**:
1.0704 at 60 (out), 1.1035 at 70, 1.1339 at 80 (in) — and 80 is the shared operating point only
because the *direct* set exhausts its free SO2 at 60. So the cascade's in-band result is a
statement about a dose the other alternative forced, with ~3.4% headroom. Asserted as the
straddle it is (``test_cascade_brackets_miaos_lower_bound_across_its_valid_dose_range``), not as
a single in-band number — which is the mistake this file already made once, in the other
direction.

**The agreement is at an operating point above Miao's wines — but NOT because the sim under-binds
SO2, which is what this file used to say (D-143 withdrew it).** Miao's buffering capacity is the
slope of an SO2 **addition** series on an unchanged wine (his Table 3 caption, n = 15, free SO2
swept ~4 -> ~52 mg/L), i.e. a **secant of the equilibrium locus**. The number this file compared
against it was a slope of the **oxidation** path, taken at the sim's own 37-59 mg/L free SO2 —
a different path at a different operating point. Measured like-for-like, on his span and the sim's
own four-carbonyl budget, the sim reads **1.3197, inside his 1.2526-1.9882**, so the buffering
gap did not exist. ``test_sulfite_buffering_matches_miao_on_his_own_statistic`` now measures his
statistic; ``test_the_oxidation_path_slope_is_a_different_quantity`` keeps the two apart so the
comparison cannot be silently re-made.

**Tier: that in-band result is ``speculative`` and must travel that way.** The secant here is set
almost entirely by pyruvate (0.3407 mM) and alpha-ketoglutarate (0.1369 mM); acetaldehyde is
~0 at this operating point and moves it ~0.01 anyway. Both keto-acid residuals are quasi-steady
ratios of ``tier: speculative`` author-estimate rate constants (D-49/D-50). At those pools' own
shipped uncertainty bounds the secant spans **1.1242 (low) - 1.3370 (nominal) - 2.1046 (high)**,
which BRACKETS Miao's band — the low corner falls below his floor. So "in band" is asserted at
nominal only, and nothing here may be written as though it survived the uncertainty.

**What IS still true is the operating point, stated without a cause.** To hold free SO2 clear of
the low-sulfite curvature regime through both rounds the sim must be dosed to ~59 mg/L free,
roughly 2x Miao's highest wine (Table 1: 9-28 mg/L), so the band comparisons here remain envelope
checks rather than a like-for-like validation, and the quinone-branching question D-141 left open
is NOT settled by them. ``test_the_operating_point_sits_above_miaos_wines`` asserts that, and
deliberately asserts **no explanation**: binding is ruled out above, and what replaces it — the
sim depleting free SO2 faster per mg of O2 than his Table 4 factors — is itself dose-dependent
here (mass-basis 1.47 at dose 40 cascade to 3.01 at dose 80 direct, against his 0.7456-1.8945),
so it is an open question (D-143), not a finding this file may encode.
"""

from typing import NamedTuple

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.core.chemistry import (
    M_5_OXOFRUCTOSE,
    M_ACETALDEHYDE,
    M_ALPHA_KETOGLUTARATE,
    M_O2,
    M_PYRUVATE,
    M_SO2,
)
from fermentation.runtime.schedule import simulate_scheduled
from fermentation.scenario import Intervention, Scenario, TemperaturePoint
from fermentation.scenario.compile import compile_scenario

pytestmark = pytest.mark.benchmark

# Molar masses come from ``core.chemistry`` (which builds them from atomic masses) rather than
# being retyped here. This file used to carry ``M_SO2 = 64.06`` against the model's 64.058 — a
# 3.1e-5 relative gap, far below any band asserted here, but every mg/L <-> mol/L conversion in
# this file is on the SAME locus the model solves, so agreeing to the last digit costs nothing and
# removes a units-edge discrepancy from the round-trip guard below.

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

#: Miao Table 3 — the sulfite buffering capacity (Δtotal SO2 / Δfree SO2) across his eight wines.
#: **Its method is load-bearing**: the caption reads "determined by linear regression analysis of
#: eight wine samples (n = 15)", and Figure 5 names it — "Total SO2 versus free SO2 **by the SO2
#: addition method**". So this is a SECANT of the equilibrium locus of an *unchanged* wine over
#: free SO2 ~4 -> ~52 mg/L, at fixed carbonyl pools. It is NOT a slope of the oxidation path, and
#: comparing one to the other is the category error D-143 withdrew (see the module docstring).
MIAO_BUFFERING_BAND = (1.2526, 1.9882)
#: The free-SO2 span Miao's addition series covers (Figure 5, axes free 0-60). The sim-side
#: secant must be evaluated over THIS span to be the same statistic.
MIAO_ADDITION_SPAN_MGL = (4.0, 52.0)
#: Miao Table 1 — the free SO2 his eight wines START the challenge at. Wine #1's 9 is the one he
#: excludes; his seven band-defining wines span 15-28.
MIAO_WINE_FREE_SO2_MGL = (15.0, 28.0)

#: The floor free SO2 must clear at EVERY sample point for a ratio measured here to be read as a
#: SO2:O2 ratio at all. **This is a sim-side conservatism, not Miao's exclusion criterion**, and
#: the distinction matters because the name it used to carry claimed otherwise.
#:
#: What it forbids: measuring the ratio in the low-sulfite regime where, in Miao's own words,
#: "oxidation reactions must necessarily involve other components, and thus the ratio of reaction
#: with O2 must change". That prose is his; the number 10.0 and the per-sample form are NOT.
#:
#: Miao's own printed criteria are two, and NEITHER is this. He excludes wine #1 on its *starting*
#: level ("the very low free SO2 of this wine, which started at 8 mg/L") and on its r2 of 0.8103
#: against >0.91 for the rest. Applied here, measured, both fail to do this guard's job:
#:
#: * **Starting level selects the opposite doses.** The sim starts at 25.0 mg/L free at dose 40
#:   (inside his 15-28) and 41.5-77.6 at doses 60-100 (all above his highest wine). His criterion
#:   would keep the dose this file measured its inverted headline at, and reject the one it fixed
#:   it with.
#: * **r2 does not discriminate at all.** The sim's ratio regression reads 0.98975 (dose 40) to
#:   0.99987 (dose 100) — every dose clears 0.91 by a wide margin, because a simulated series has
#:   none of the analytical scatter his criterion was screening for.
#:
#: And this floor is STRICTER than his data: read as free SO2 at O2 exhaustion, his Table 4
#: intercepts are 2.7369 / 6.3826 / 9.0102 / 7.2373 / 7.8348 / 12.258 / 13.544 / 17.576 for wines
#: 1-8, so **four of the seven wines that define the band would fail it** (2/3/4/5). That is
#: recorded, not fixed: the strictness is deliberate. The one defect this file has ever shipped
#: was measuring inside the curvature regime — free SO2 ending at 2.8 mg/L — and it inverted the
#: headline. A conservative sim-side floor is the guard against a repeat; it is simply not
#: entitled to Miao's name. [[feedback-name-guards-for-what-they-forbid]]
SIM_CURVATURE_FLOOR_MGL = 10.0

#: The shared operating point: the lowest dose ON THE TESTED GRID {40, 60, 70, 80, 100, 120} at
#: which free SO2 stays above :data:`SIM_CURVATURE_FLOOR_MGL` through the whole window on BOTH
#: alternatives — enforced by ``test_the_operating_point_clears_the_curvature_floor``
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


#: Every carbonyl the SHIPPED binder partitions bisulfite across (D-28 acetaldehyde, D-51's two
#: keto-acids, D-130's botrytis 5-oxofructose), as (state slot, dissociation-constant parameter,
#: molar mass). Read from ``acidbase`` rather than transcribed, so a fifth binder cannot be added
#: to the model and silently omitted from this file's statistic. 5-oxofructose is 0 in a
#: non-botrytis white, but dropping a binder because it happens to be zero HERE is how a scenario
#: with one of them non-zero would later read a wrong number without anything going red.
_SO2_BINDERS = (
    (acidbase.ACETALDEHYDE_KEY, acidbase.SO2_BINDING_PARAM, M_ACETALDEHYDE),
    (acidbase.PYRUVATE_KEY, acidbase.PYRUVATE_SO2_BINDING_PARAM, M_PYRUVATE),
    (acidbase.ALPHA_KG_KEY, acidbase.ALPHA_KG_SO2_BINDING_PARAM, M_ALPHA_KETOGLUTARATE),
    (acidbase.OXOFRUCTOSE_KEY, acidbase.OXOFRUCTOSE_SO2_BINDING_PARAM, M_5_OXOFRUCTOSE),
)


def _binder_pools(y, schema, params) -> tuple[tuple[float, float], ...]:
    """``(pool mol/L, K mol/L)`` for each of the four shipped SO2 binders, at one state."""
    pools = []
    for key, k_param, molar_mass in _SO2_BINDERS:
        amount = max(float(y[schema.slice(key)][0]), 0.0) / molar_mass if key in schema else 0.0
        pools.append((amount, float(params.get(k_param, 0.0))))
    return tuple(pools)


def _total_so2_at_free(free_mgl: float, beta: float, pools) -> float:
    """The equilibrium locus ``total(free)`` [mg/L as SO2] at fixed carbonyl pools.

    The same competitive-Langmuir system :func:`acidbase.bound_so2_molar` root-finds, evaluated
    from the other end: given free SO2, ``h = beta*u`` is known outright and every bound amount
    follows in closed form, so no solve is needed. That the two are one locus is asserted by
    ``test_the_analytic_locus_is_the_shipped_binding_solver`` rather than assumed.
    """
    u = free_mgl / 1000.0 / M_SO2
    h = beta * u
    bound = float(sum(a * h / (k + h) for a, k in pools if a > 0.0 and k > 0.0))
    return (u + bound) * M_SO2 * 1000.0


def _addition_method_secant(y, schema, params) -> float:
    """Miao's Table 3 statistic, computed on the sim: d(total)/d(free) as a SECANT over his span.

    An addition series holds the wine — hence the carbonyl pools and the pH — fixed and sweeps
    the SO2 that is added, so this reads the pools out of ONE state and moves only free SO2. That
    is what makes it comparable to :data:`MIAO_BUFFERING_BAND` and what the oxidation-path slope
    (where the pools, the pH and the free SO2 all move together, down a different path) is not.
    """
    ph = acidbase.ph_of_state(y, schema, params)
    pkas = tuple(float(params[n]) for n in acidbase.SO2_PKA_PARAM_NAMES)
    beta = acidbase.bisulfite_fraction(10.0**-ph, pkas)
    pools = _binder_pools(y, schema, params)
    lo, hi = MIAO_ADDITION_SPAN_MGL
    return (_total_so2_at_free(hi, beta, pools) - _total_so2_at_free(lo, beta, pools)) / (hi - lo)


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
    """One accelerated-oxidation run, reduced to the things Miao reports plus the guards."""

    #: Miao Table 2 — d(total SO2) / d(O2 consumed), molar.
    ratio: float
    #: Miao Table 3 — the sulfite buffering capacity, by HIS method: the addition-series secant
    #: of the equilibrium locus over :data:`MIAO_ADDITION_SPAN_MGL`, at the window-start state.
    #: This is the only quantity here comparable to :data:`MIAO_BUFFERING_BAND`.
    buffering_secant: float
    #: d(total SO2) / d(free SO2) along the OXIDATION path — the quantity D-142 mistook for the
    #: one above. Kept, and kept separate, because the whole finding is that they differ.
    oxidation_path_slope: float
    #: Free SO2 [mg/L] at each of the eight sample points, for the curvature-floor guard.
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
    # The addition-series secant is read at the state the challenge window OPENS at — the "wine
    # as bottled", which is what Miao adds SO2 to. It is a property of the carbonyl pools and the
    # pH, so it does not move with the O2 challenge and must not be regressed along it.
    start = traj.y[:, int(np.argmin(np.abs(traj.t - hours[0])))]
    return _Run(
        ratio=ratio,
        buffering_secant=_addition_method_secant(start, compiled.schema, params),
        oxidation_path_slope=ratio / free_factor if free_factor else float("nan"),
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


def test_the_operating_point_clears_the_curvature_floor(runs):
    """Free SO2 stays above :data:`SIM_CURVATURE_FLOOR_MGL` at every sample point, BOTH sets.

    This enforces in code the bound the first version of this file only stated in prose, and it
    is the guard that would have caught that defect: at the original 40 mg/L dose free SO2 ended
    at 2.8 mg/L (direct) and 6.8 (cascade), deep in the low-sulfite regime where the ratio's
    meaning changes. Every band comparison below is meaningless without it.

    Named for what it forbids rather than for Miao, because it is **not** his exclusion criterion
    — see :data:`SIM_CURVATURE_FLOOR_MGL`, which records both of his printed criteria, the
    measurement showing neither does this job here, and the four band-defining wines this floor
    is strict enough to exclude.
    """
    for alternative in ("direct", "cascade"):
        free = runs[(alternative, "real", WINE_REALISTIC_SO2)].free_mgl
        assert min(free) > SIM_CURVATURE_FLOOR_MGL, f"{alternative} fell to {min(free):.1f} mg/L"


def test_sulfite_buffering_matches_miao_on_his_own_statistic(runs):
    """The sim's addition-method secant lands INSIDE Miao's 1.2526-1.9882, on both alternatives.

    This replaces an assertion that claimed the opposite. D-142 asserted ``buffering <
    MIAO_BUFFERING_BAND[0]`` — "the sim under-binds SO2" — from a slope of the OXIDATION path
    taken at the sim's own 37-59 mg/L free SO2. Miao's number is the secant of an SO2 **addition**
    series over free 4-52 mg/L. Two mismatches compounded: wrong path, and free SO2 2-6x his.
    Correct both and the sim reads ~1.32, near his floor, matching his wines 8 (1.2526),
    3 (1.2734) and 4 (1.3955). There was no binding gap to close. (D-143)

    Identical on both alternatives, and at every dose, because the secant is a property of the
    carbonyl pools and the pH — neither of which the oxidative alternative or the SO2 dose moves.
    That invariance is asserted here too: it is the reason this statistic, unlike the ratio, is
    not an artefact of the operating point.

    **Tier: speculative** (Prime Directive 1). Pyruvate and alpha-ketoglutarate carry the whole
    secant and are quasi-steady ratios of ``tier: speculative`` author estimates; at their own
    shipped uncertainty bounds the secant spans 1.1242-2.1046, which brackets Miao's band. So the
    assertion is deliberately the nominal-pool one, and the tolerance is tight rather than
    generous — a wide band here would be pretending the corners had been checked. They have, in
    D-143, and the low one fails.
    """
    lo, hi = MIAO_BUFFERING_BAND
    secants = {
        alt: runs[(alt, "real", WINE_REALISTIC_SO2)].buffering_secant
        for alt in ("direct", "cascade")
    }
    for alternative, secant in secants.items():
        assert lo <= secant <= hi, f"{alternative} secant {secant:.4f} outside Miao's {lo}-{hi}"
    assert secants["direct"] == pytest.approx(secants["cascade"], abs=1e-9)


def test_the_oxidation_path_slope_is_a_different_quantity(runs):
    """Guards the category error itself: the two slopes are NOT interchangeable.

    Keeping the oxidation-path number in :class:`_Run` is only safe if something asserts it is a
    different quantity from Miao's, otherwise a later reader re-makes D-142's comparison with the
    field that is right there. The two disagree, they disagree in a direction (the oxidation slope
    runs LOW, which is exactly why the old assertion looked true), and the gap is far outside any
    tolerance either could be measured to.

    Measured at the shared operating point: secant 1.3197 against oxidation slopes of 1.1761
    (direct) and 1.1333 (cascade).
    """
    for alternative in ("direct", "cascade"):
        run = runs[(alternative, "real", WINE_REALISTIC_SO2)]
        assert run.oxidation_path_slope < run.buffering_secant - 0.05, (
            f"{alternative}: the oxidation-path slope {run.oxidation_path_slope:.4f} has "
            f"converged on the addition secant {run.buffering_secant:.4f} — if that is real, the "
            "two statistics may finally be comparable; re-read D-143 before assuming so"
        )


def test_the_analytic_locus_is_the_shipped_binding_solver():
    """The secant's closed form must BE ``acidbase.bound_so2_molar``, not a re-derivation of it.

    :func:`_total_so2_at_free` evaluates the competitive-Langmuir system from free SO2; the
    shipped code root-finds the same system from total. If the two ever parted, the secant would
    be measuring a private model of binding and reporting it as the sim's.

    Asserted as **pure algebra** — free -> total analytically, then total -> free through
    ``bound_so2_molar`` — with no state, no pH solve and no scenario, because that is the whole of
    the claim. (An earlier draft round-tripped through ``free_so2`` on a state vector instead and
    read errors up to 2.5e-4 mg/L, none of it from the binding solve: it was dragging in
    ``ph_of_state`` and the 64.06-vs-64.058 molar mass since removed. The lesson is the file's own
    recurring one — measure the statistic you mean, not one entangled with it.)

    Tolerance is the measured floor, not a round number: over a 120-point grid (5 pH x 3
    acetaldehyde pools x 2 keto-acid scalings x 4 free SO2 levels) the worst absolute error is
    **4.3e-10 mg/L**, so 1e-8 is ~20x the noise floor and still ~1e7 x tighter than anything this
    file concludes from. [[feedback-pin-tolerance-vs-solver-tolerance]]
    """
    params = dict(compile_scenario(_scenario(WINE_REALISTIC_SO2), oxidative="direct").param_values)
    k_values = [float(params[k_param]) for _, k_param, _ in _SO2_BINDERS]

    for ph in (2.90, 3.14, 3.35, 3.58, 3.80):
        beta = acidbase.bisulfite_fraction(
            10.0**-ph, tuple(float(params[n]) for n in acidbase.SO2_PKA_PARAM_NAMES)
        )
        # Every binder loaded, so the round trip exercises the coupled 4-carbonyl solve rather
        # than the degenerate keto-acids-only case this scenario happens to present.
        pools = tuple((0.5e-3, k) for k in k_values)
        for free_mgl in MIAO_ADDITION_SPAN_MGL:
            total_molar = _total_so2_at_free(free_mgl, beta, pools) / 1000.0 / M_SO2
            bound = acidbase.bound_so2_molar(total_molar, pools, beta)
            recovered = (total_molar - sum(bound)) * M_SO2 * 1000.0
            assert recovered == pytest.approx(free_mgl, abs=1e-8), f"pH {ph}, free {free_mgl}"


def test_the_operating_point_sits_above_miaos_wines(runs):
    """The band comparisons here are envelope checks, not a like-for-like validation.

    To hold free SO2 clear of :data:`SIM_CURVATURE_FLOOR_MGL` through both rounds the sim opens
    the window at ~59 mg/L free — roughly 2x wine #8's 28, the highest Miao reports. So no dose
    on this grid puts the sim inside his Table 1 range AND clear of the curvature regime, and the
    quinone-branching question D-141 left open is not settled by any agreement measured here.

    **This test asserts the fact and deliberately asserts no cause.** The cause this file used to
    give — a bound-SO2 reservoir too small to replenish free SO2 — is withdrawn: the buffering
    capacity is in band on Miao's own statistic (see
    ``test_sulfite_buffering_matches_miao_on_his_own_statistic``), which is asserted here as the
    negative half. What is left is that the sim spends free SO2 faster per mg of O2 than his
    Table 4 factors (0.7456-1.8945 mg/mg), but that gap is itself dose-dependent here — 1.47 at
    dose 40 on the cascade, 3.01 at dose 80 on the direct set — so it is an open question (D-143),
    and a test that named it as the reason would be the third dose-contingent verdict this file
    has recorded as though it were structural.
    """
    lo, hi = MIAO_BUFFERING_BAND
    for alternative in ("direct", "cascade"):
        run = runs[(alternative, "real", WINE_REALISTIC_SO2)]
        assert run.free_mgl[0] > MIAO_WINE_FREE_SO2_MGL[1], (
            f"{alternative} opens at {run.free_mgl[0]:.1f} mg/L free SO2, now within Miao's "
            f"{MIAO_WINE_FREE_SO2_MGL} — the operating-point gap has closed; update D-143"
        )
        assert lo <= run.buffering_secant <= hi, (
            "binding is NOT the reason: if the secant has left Miao's band the explanation "
            "above no longer holds and this docstring is stale"
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
    assert min(runs[("cascade", "real", low)].free_mgl) > SIM_CURVATURE_FLOOR_MGL
    assert min(runs[("direct", "real", low)].free_mgl) < SIM_CURVATURE_FLOOR_MGL


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
