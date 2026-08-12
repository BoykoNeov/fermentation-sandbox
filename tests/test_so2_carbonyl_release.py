"""The carbonyl-release side of the SO₂ equilibrium — what it reproduces, and where it stops
(decision D-190).

D-137 surfaced, and D-139/D-140 carried forward three times as "deserves its own D-record", the
claim from Bueno *et al.* (quoted in *Understanding Wine Chemistry* 2e Ch. 24) that in real wine
the rise in oxidised-aroma carbonyls is "well modeled by the expected release from bound SO₂
adducts" until molecular SO₂ < 0.05 mg/L, with de novo formation mattering only once the wine is
severely oxidised. The model builds two de-novo routes (D-71 Fenton, D-75 Strecker) and carries a
stateless competitive-Langmuir free/bound equilibrium (D-28, D-51, D-130). **Nothing tested the
ordering, and nothing tested the release.** D-190 measured it. These are the guards that measurement
licensed.

What is pinned here, and why each one is not covered elsewhere:

  * ``test_published_bound_fraction_*`` — the equilibrium reproduces a **published two-sided
    literal** (">99% bound when free SO₂ >30 mg/L; >95% when >2 mg/L", *Understanding Wine
    Chemistry* 2e Ch. 24 citing Ch. 17). Nothing in the suite referenced the bound *fraction*
    against any source before this file; the agreement was real but untested (D-187's lesson —
    test against published literals).
  * ``test_the_weak_edge_of_the_band_misses_the_95_percent_floor`` — and the band does **not**
    carry that agreement. The nominal clears the 95 % floor by ~0.2 pp; the band's weak edge
    misses it by ~1.5 pp. Pinned so the scope is a measured statement rather than a silent fact
    ([[feedback-pin-the-band-not-the-nominal]]).
  * ``test_the_asserted_pyruvate_alpha_kg_ordering_is_reachably_false`` — the two keto-acid Kd
    bands **overlap**, are sampled independently, and admit the reversal of the ordering the
    parameter file used to assert as fact. This guards the corrected prose.
  * ``test_no_aroma_carbonyl_binds_so2`` / ``test_the_aroma_readout_reads_total_not_free`` — the
    two scope statements D-190 refused to build past. Both are *documented limitations*, and both
    are pinned so that closing either is a deliberate act with a red test, not a silent drift.

**This file adds no physics and moves no value.** Every number in it is either published, derived
from the shipped parameters, or measured off the shipped equilibrium.
"""

import numpy as np
import pytest

from fermentation.core.acidbase import (
    _bound_molar_split,
    bisulfite_fraction,
    free_acetaldehyde,
    ph_of_state,
)
from fermentation.core.chemistry import M_ACETALDEHYDE, M_SO2
from fermentation.core.media import wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.ensemble import sample_parameters
from fermentation.runtime.integrate import Trajectory
from fermentation.sensory.oav import AROMA_COMPOUNDS, load_thresholds, oav_series
from fermentation.units.convert import ugl_to_gpl

#: The SO₂ acidity constants, by their real names in ``acidbase.yaml``.
SO2_PKA_NAMES = ("pKa_sulfurous_1", "pKa_sulfurous_2")

#: *Understanding Wine Chemistry* 2e, Ch. 24 (locally held), quoting Ch. 17 verbatim:
#: ">99% of acetaldehyde will be bound when free SO2 is >30 mg/L and >95% is bound when free SO2
#: is > 2 mg/L". Two published points, both one-sided floors. Transcribed from the text, not
#: from a summary of it.
PUBLISHED_FLOORS: tuple[tuple[float, float], ...] = ((30.0, 0.99), (2.0, 0.95))

#: The wine pH range these guards sweep. The published claim states no pH, so a pin that held at
#: only one would be reading a coincidence.
WINE_PH = (3.0, 3.3, 3.5, 3.8)


@pytest.fixture
def pset():
    return load_parameters(default_data_dir() / "acidbase.yaml")


@pytest.fixture
def params(pset):
    return pset.resolve()


def _constant_traj(schema: StateSchema, pools: dict[str, float], *, n: int = 4) -> Trajectory:
    """A synthetic constant-in-time trajectory with the named pools set (all else 0).

    Same idiom as ``tests/test_sensory_oav.py``: building the frozen Trajectory directly keeps
    this a test of the READOUT rather than of any solver run.
    """
    y: FloatArray = np.zeros((schema.size, n), dtype=np.float64)
    for pool, val in pools.items():
        y[schema.slice(pool), :] = val
    return Trajectory(
        schema=schema,
        t=np.linspace(0.0, 1.0, n),
        y=y,
        success=True,
        message="",
        tier_map=dict.fromkeys(schema.names, Tier.SPECULATIVE),
    )


def _bound_fraction(free_so2_mgl: float, k_acet: float, ph: float, pkas) -> float:
    """Bound share of acetaldehyde at a held FREE SO₂ — the basis the published claim uses.

    The claim is conditioned on free SO₂ rather than total, so this evaluates the Langmuir term
    at that free pool instead of partitioning a total. The acetaldehyde load cancels out of the
    ratio (``x/A = h/(K+h)``), which is why no load appears in the signature — and which is worth
    knowing, because the published claim is likewise quoted without one.
    """
    h = bisulfite_fraction(10.0 ** (-ph), pkas) * (free_so2_mgl / 1000.0) / M_SO2
    return h / (k_acet + h)


# ---------------------------------------------------------------------------------------------
# The published literal
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("ph", WINE_PH)
@pytest.mark.parametrize(("free_mgl", "floor"), PUBLISHED_FLOORS)
def test_published_bound_fraction_holds_at_the_nominal(pset, params, ph, free_mgl, floor):
    """At the shipped Kd, acetaldehyde is bound in the published proportion at both points.

    **This is an implementation check across a basis change, NOT independent corroboration of the
    constant** (D-190 amendment). Table 17.2 of the same book prints acetaldehyde Kd = 1.5e-6,
    bit-identical to the shipped value, and the Ch. 24 prose is Ch. 17's own arithmetic on it — so
    agreement on the *constant* is guaranteed and proves nothing. What is NOT guaranteed is the
    basis: this model references binding to bisulfite HSO₃⁻ with a pH-dependent β, while the
    textbook's apparent Kd is quoted per TOTAL free SO₂ at a stated pH. Passing at 30 and 2 mg/L
    across pH 3.0-3.8 is what pins that the two bases agree to within the ≤5 % the parameter
    file's own note claims for them.
    """
    pkas = tuple(params[n] for n in SO2_PKA_NAMES)
    frac = _bound_fraction(free_mgl, params["K_acetaldehyde_so2"], ph, pkas)
    assert frac > floor, (
        f"at pH {ph} and free SO2 {free_mgl} mg/L the model binds {frac:.4%} of acetaldehyde, "
        f"below the published floor of {floor:.0%}"
    )


@pytest.mark.parametrize("ph", WINE_PH)
def test_published_99_percent_floor_holds_across_the_WHOLE_band(pset, ph):
    """The 30 mg/L point survives both band edges — so *that* one is a band-wide property.

    Stated separately from the 2 mg/L point because the two behave differently, and collapsing
    them into one "the published claim holds" test would hide that.
    """
    p = pset["K_acetaldehyde_so2"]
    pkas = tuple(pset.resolve()[n] for n in SO2_PKA_NAMES)
    for label, k in (("low", p.uncertainty.low), ("value", p.value), ("high", p.uncertainty.high)):
        frac = _bound_fraction(30.0, k, ph, pkas)
        assert frac > 0.99, f"{label} edge binds only {frac:.4%} at 30 mg/L free SO2 (pH {ph})"


@pytest.mark.parametrize("ph", WINE_PH)
def test_the_weak_edge_of_the_band_misses_the_95_percent_floor(pset, ph):
    """The 2 mg/L point is a property of the NOMINAL, not of the band — measured, not asserted.

    ``K_acetaldehyde_so2``'s band high (2.1e-6, the pH-3.5 apparent Kd reported in the wine
    literature) binds 93.3-93.6 % at 2 mg/L free SO₂, against a published floor of 95 %. Both
    numbers are sourced and they disagree at the edge; the nominal clears the floor by only
    ~0.2 pp. This is pinned rather than fixed because moving either would mean choosing between
    two sourced quantities on no evidence — exactly the fit D-190 refused. The guard's job is to
    make a future band change surface the conflict instead of burying it.
    """
    p = pset["K_acetaldehyde_so2"]
    pkas = tuple(pset.resolve()[n] for n in SO2_PKA_NAMES)

    at_edge = _bound_fraction(2.0, p.uncertainty.high, ph, pkas)
    at_nominal = _bound_fraction(2.0, p.value, ph, pkas)

    assert at_edge < 0.95, (
        "the band's weak edge no longer misses the published 95% floor; if the band moved, "
        "re-read D-190 §5 before deleting this test"
    )
    assert 0.93 < at_edge < 0.94, f"edge bound fraction {at_edge:.4%} outside the measured range"
    # The nominal's headroom is thin, and saying so is the point: a pin at 95% is ~0.2 pp of margin.
    assert 0.95 < at_nominal < 0.96, f"nominal bound fraction {at_nominal:.4%} moved"


# ---------------------------------------------------------------------------------------------
# The keto-acid ordering — asserted in prose, contradicted by the bands
# ---------------------------------------------------------------------------------------------


def test_the_two_keto_acid_bands_overlap_at_all(pset):
    """The precondition for the next test, stated on its own so a band change names its cause."""
    pyr, akg = pset["K_pyruvate_so2"], pset["K_alpha_kg_so2"]
    lo = max(pyr.uncertainty.low, akg.uncertainty.low)
    hi = min(pyr.uncertainty.high, akg.uncertainty.high)
    assert lo < hi, (
        "the pyruvate and alpha-KG Kd bands no longer overlap; the ordering reversal D-190 "
        "measured is then unreachable and the parameter-file prose should be revisited"
    )


def test_the_asserted_pyruvate_alpha_kg_ordering_is_reachably_false(pset):
    """~2 % of independent draws make pyruvate the STRONGER binder — the reverse of the nominal.

    Until D-190 the ``K_alpha_kg_so2`` header comment stated the ordering as a fact ("slightly
    STRONGER than pyruvate ... the second-most-avid binder"). The nominal pair does say that; the
    *bands* do not, and they are sampled independently (``runtime/ensemble.py``). Measured on
    DRAWS rather than edges, which is the D-154 rule — an edge argument would have proved only
    that the corner exists, not that it is visited.

    The rate is pinned loosely (a band, not a literal) because it is a property of the two
    triangular distributions, not a measurement of anything physical.
    """
    rng = np.random.default_rng(190)
    n = 20_000
    reversals = sum(
        1
        for _ in range(n)
        if (s := sample_parameters(pset, rng, names=["K_pyruvate_so2", "K_alpha_kg_so2"]))[
            "K_pyruvate_so2"
        ]
        < s["K_alpha_kg_so2"]
    )
    frac = reversals / n
    assert 0.005 < frac < 0.05, (
        f"{frac:.4%} of draws reverse the pyruvate/alpha-KG ordering; D-190 measured ~2.0% over "
        f"200k draws. A large move here means a band moved — re-read D-190 §5."
    )


def test_a_pure_transposition_of_the_two_nominals_is_outside_the_shipped_bands(pset):
    """Swapping the two shipped nominals is *forbidden*, by pyruvate's own band low, by ~7 %.

    Recorded because it is the first hypothesis anyone will reach for on reading D-190 §5 (the
    shipped alpha-KG value is bit-identical to another book's pyruvate cell), and because it is
    why D-190's first probe of the swap errored at parameter load rather than producing physics.
    """
    pyr, akg = pset["K_pyruvate_so2"], pset["K_alpha_kg_so2"]
    assert akg.value < pyr.uncertainty.low, (
        "a straight transposition now lies inside pyruvate's band; D-190's argument that its "
        "~15% figures are unreachable no longer holds"
    )


# ---------------------------------------------------------------------------------------------
# The two scope statements D-190 refused to build past
# ---------------------------------------------------------------------------------------------


def test_no_aroma_carbonyl_binds_so2(params):
    """The carbonyls Bueno's claim is *about* bind nothing here — so release cannot exist for them.

    ``_bound_molar_split`` takes exactly four carbonyls (acetaldehyde, pyruvate, alpha-KG,
    5-oxofructose). Methional and phenylacetaldehyde — the oxidised-aroma aldehydes — are aroma
    pools with no adduct, so the model has no release term to compare against its de-novo routes
    and Bueno's ordering is not expressible for them. D-190 refused to build the adduct: no Kd for
    methional exists in any locally-held source, and *Understanding Wine Chemistry*'s Table 17.2,
    which does carry Kd values for other odorants, does not list it.

    This is a **limitation pinned as a limitation**. Closing it needs a sourced Kd, not a guess.
    """
    aroma = {c.pool for c in AROMA_COMPOUNDS["wine"]}
    binder_params = {
        "K_acetaldehyde_so2",
        "K_pyruvate_so2",
        "K_alpha_kg_so2",
        "K_5_oxofructose_so2",
    }
    assert binder_params <= set(params), "a binding constant vanished from the parameter set"
    for pool in ("methional", "phenylacetaldehyde"):
        assert pool in aroma, f"{pool} is no longer a wine aroma pool"
        assert f"K_{pool}_so2" not in params, (
            f"{pool} has acquired a binding constant — if that is deliberate, D-190's scope "
            f"section and this test both need rewriting, and the Kd needs a source"
        )


def test_the_split_takes_exactly_four_carbonyls(params):
    """A structural companion to the above: the arity is the scope.

    Pinned because the release story is entirely determined by which pools are in this call, and
    a fifth arriving without a D-record is the failure mode.
    """
    pkas = tuple(params[n] for n in SO2_PKA_NAMES)
    beta = bisulfite_fraction(10.0 ** (-3.4), pkas)
    out = _bound_molar_split(
        60e-3 / M_SO2, 5e-3 / M_ACETALDEHYDE, 0.0, 0.0, 0.0, beta, params
    )  # fmt: skip
    assert len(out) == 4, f"the carbonyl split returned {len(out)} shares, not 4"
    assert out[0] > 0.0, "acetaldehyde must stay at index 0 (free_acetaldehyde reads that share)"
    assert out[1] == out[2] == out[3] == 0.0, "absent carbonyls must contribute exactly zero"


def test_the_aroma_readout_reads_total_not_free(params):
    """The OAV divides the TOTAL pool by the threshold — including for acetaldehyde.

    Three chemistry consumers read ``free_acetaldehyde`` (the D-27 reduction, the D-80 bridging
    and one aging term); the sensory readout reads none of them. In a sulfited wine that makes the
    reported acetaldehyde OAV two-to-three orders of magnitude above the perceptible quantity.
    D-190 did **not** change this: with acetaldehyde at ~0.5 mg/L against a ~100 mg/L threshold
    the OAV is far below 1 either way, so switching the readout would move no verdict while
    quietly making every other pool's OAV inconsistent with it (only acetaldehyde has a free/bound
    split to read).

    Asserted on BEHAVIOUR, not on the source text: a `getsource` check would break on any
    harmless refactor and would pass if someone changed the behaviour while keeping the string
    ([[feedback-grep-finds-claims-not-guards]], one level up). Pinned so the asymmetry is a
    recorded choice. The pool where it would bite is methional, and that one is blocked above.
    """
    schema = wine_schema()
    acet_gpl, so2_gpl = 5.0e-3, 60.0e-3
    traj = _constant_traj(schema, {"acetaldehyde": acet_gpl, "so2_total": so2_gpl})
    thresholds = load_thresholds()
    threshold_gpl = ugl_to_gpl(thresholds.value("threshold_acetaldehyde_wine"))

    reported = float(oav_series(traj, thresholds, "acetaldehyde")[-1])
    assert reported == pytest.approx(acet_gpl / threshold_gpl, rel=1e-12), (
        "oav_series no longer reports total/threshold for acetaldehyde; if it now reads the free "
        "share, D-190's scope section is stale and this test should assert the new behaviour"
    )

    y = traj.y[:, -1]
    free = free_acetaldehyde(y, schema, params, ph_of_state(y, schema, params))
    assert free < acet_gpl, "no SO2 binding happened — the fixture is not exercising the split"
    # The gap is the finding: the readout is orders of magnitude above the perceptible share.
    assert reported / (free / threshold_gpl) > 100.0, (
        "the total-vs-free gap has collapsed; re-read D-190 §2 before trusting either number"
    )
