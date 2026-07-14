"""Odor-Activity-Value (OAV) sensory readout — the speculative aroma lens (decision D-67).

The handoff §4.2 sensory layer, built as a **pure readout** over aroma-active compounds the
chemistry already tracks (esters, fusels, diacetyl, acetaldehyde, H₂S, and — in wine —
4-ethylphenol, 4-ethylguaiacol, volatile thiols). It adds **no ODE state, no Process, no
ledger entry**: it consumes a finished :class:`~fermentation.runtime.integrate.Trajectory`
plus a threshold table and returns dimensionless OAVs. This mirrors how
:mod:`fermentation.analysis` maps the pH/SO₂/IBU series over a trajectory — a thin observable
sitting one layer up.

**The §4.2 cardinal rule / isolation firewall.** The sensory layer consumes the chemistry;
the chemistry *never* depends on the sensory layer. Nothing in ``core``/``runtime``/
``scenario`` imports this module. And the thresholds load **standalone** via
:func:`load_thresholds` — a new ``sensory.yaml`` that is **not** merged into any
``CompiledScenario`` at the compile seam (no RHS reads a perception threshold), so the
chemistry never even sees these numbers. Dosing/aroma chemistry is therefore byte-for-byte
unchanged by this readout.

**OAV.** ``OAV_i(t) = concentration_i(t) / threshold_i`` for aroma compound *i*. Above 1 the
compound is (individually) above its perception threshold. Both sides are compared in the
canonical g/L: every threshold in ``sensory.yaml`` is stored in µg/L (the human-readable
literature unit) and crossed to g/L at this boundary via
:func:`fermentation.units.convert.ugl_to_gpl`, so the ratio is directly dimensionless.

**The tier floor (§4.3 credibility firewall).** Every OAV output tier is
``combine(input_tier, threshold_tier, SPECULATIVE)`` → **speculative**, *even when the input
chemistry is validated*. The explicit ``SPECULATIVE`` is **not** redundant with the
threshold's own tier: the sensory *mapping itself* is the canonical speculative case (the
:class:`~fermentation.core.tiers.Tier` docstring names "sensory mapping"), so the floor must
hold even if a threshold were later mislabelled plausible. See :func:`oav_tier`.

**The lumped-pool call (D-66).** ``esters``/``fusels``/``mercaptans`` are single g/L pools
mixing several molecules; each is read against the threshold of one **named representative**
(the stand-in its ``VarSpec`` already names), with the "assumes fixed lump composition"
honesty cost flagged in that threshold's provenance ``notes``. The single-molecule pools
(diacetyl, acetaldehyde, H₂S, 4-EP, 4-EG) carry no such assumption.

**Medium-specificity.** The aroma set is medium-specific (mirroring the beer-only
``iso_alpha``): beer carries the 5 common pools + the 5 oak extractives (barrel-beer, D-86 + the
D-94 furaneol caramel); wine adds 4-EP, 4-EG, mercaptans, the two Strecker aldehydes (D-75) and the
four non-oxidative THERMAL Strecker aldehydes/sotolon (D-87) on top of the common + oak sets.
The five oak extractives (D-77 four + furaneol D-94) are shared by both media (D-86 — the oak axis
is a wood property).
Every threshold is matrix-specific too (``threshold_<pool>_beer`` vs ``_wine``) because ethanol and
the wine/beer matrix shift odor thresholds substantially. ``iso_alpha``/IBU is deliberately
excluded — it is a *taste* (bitterness), already read out by :func:`fermentation.analysis.
ibu_series` (D-64); an odor-threshold OAV does not apply to it.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass

import numpy as np

from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier, combine
from fermentation.parameters.store import ParameterSet, default_data_dir, load_parameters
from fermentation.runtime.integrate import Trajectory
from fermentation.units.convert import ugl_to_gpl


@dataclass(frozen=True)
class AromaCompound:
    """One aroma-active pool and how OAV maps onto it.

    ``pool`` is the state-variable name; ``representative`` is the named molecule whose
    threshold the pool is read against (identical to ``pool``'s molecule for the
    single-molecule pools, the stand-in for a ``lumped`` pool); ``descriptor`` is a **static**
    per-compound aroma label (a fixed vocabulary word, *not* a synthesised "this smells like
    …" — that projection is the deferred beat 1b). ``lumped`` flags that the OAV rides on a
    fixed-lump-composition assumption (D-66).
    """

    pool: str
    representative: str
    descriptor: str
    lumped: bool


#: The five aroma-active pools shared by both media (``_common_specs`` in ``core.media``).
_COMMON: tuple[AromaCompound, ...] = (
    AromaCompound("diacetyl", "2,3-butanedione", "buttery", lumped=False),
    AromaCompound("acetaldehyde", "acetaldehyde", "green apple / bruised", lumped=False),
    AromaCompound("h2s", "hydrogen sulfide", "rotten egg", lumped=False),
    AromaCompound("esters", "isoamyl acetate", "banana / fruity", lumped=True),
    AromaCompound("fusels", "isoamyl alcohol", "solventy / fusel", lumped=True),
)

#: The nine wine-only pools appended in ``wine_schema`` (Brett phenols + volatile thiols + the two
#: shared Strecker aldehydes + the four non-oxidative THERMAL Strecker aldehydes/sotolon of D-87).
#: ``methional`` and ``phenylacetaldehyde`` (D-75) are single-molecule pools with OPPOSITE valence —
#: methional the cooked-potato oxidative off-note, phenylacetaldehyde the honey note — and are
#: *shared* by the D-75 oxidative and D-87 thermal routes (same molecules, one pool + threshold
#: each). The four D-87-only pools are the sweet-wine/Madeira thermal suite: ``2_methylbutanal`` /
#: ``3_methylbutanal`` / ``2_methylpropanal`` (the malty branched-chain aldehydes) and ``sotolon``
#: (the curry/maple furanone). All read against their own matrix-specific thresholds. (The four oak
#: aroma extractives moved to the medium-agnostic ``_OAK`` tuple at D-86 — barrel-beer oak.)
_WINE_ONLY: tuple[AromaCompound, ...] = (
    AromaCompound("ethylphenols", "4-ethylphenol", "horse-sweat / barnyard", lumped=False),
    AromaCompound("ethylguaiacols", "4-ethylguaiacol", "clove / smoky", lumped=False),
    AromaCompound("mercaptans", "methanethiol", "reductive / drains", lumped=True),
    AromaCompound("methional", "methional", "cooked potato / oxidative", lumped=False),
    AromaCompound("phenylacetaldehyde", "phenylacetaldehyde", "honey / floral", lumped=False),
    AromaCompound("2_methylbutanal", "2-methylbutanal", "malty / almond", lumped=False),
    AromaCompound("3_methylbutanal", "3-methylbutanal", "malty / dark chocolate", lumped=False),
    AromaCompound("2_methylpropanal", "2-methylpropanal", "malty / grainy", lumped=False),
    AromaCompound("sotolon", "sotolon", "curry / maple / nutty", lumped=False),
)

#: The FIVE oak aroma extractives — the non-oxidative barrel/chip diffusion axis (decision D-77 for
#: the first four, D-94 for ``furaneol``), SHARED by wine and barrel-beer (decision D-86: both media
#: carry the oak slots via ``core.media._oak_specs``). Single-molecule pools; ``guaiacol`` here is
#: the OAK smoky note, DISTINCT from the Brett ``ethylguaiacols`` above (a different molecule);
#: ``furaneol`` (HDMF) is the caramel/toffee note of toasted/charred oak + ex-bourbon barrels, its
#: descriptor kept DISTINCT from ``sotolon``'s "maple". Each read against its
#: matrix-specific threshold (``threshold_<compound>_beer`` vs ``_wine``, ``sensory.yaml``). The oak
#: *ceiling* slots (set-and-hold, D-77) and the ``ellagitannin`` TASTE pool (D-78, read by
#: :func:`~fermentation.analysis.astringency_series`, not the odor lens) are NOT aroma pools — all
#: deliberately excluded here.
_OAK: tuple[AromaCompound, ...] = (
    AromaCompound("whiskey_lactone", "whiskey lactone", "coconut / oak", lumped=False),
    AromaCompound("vanillin", "vanillin", "vanilla", lumped=False),
    AromaCompound("guaiacol", "guaiacol", "smoky / toasty", lumped=False),
    AromaCompound("eugenol", "eugenol", "clove / spice", lumped=False),
    AromaCompound("furaneol", "furaneol", "caramel / toffee", lumped=False),
)

#: Medium -> its ordered aroma set. Beer = 5 common + 5 oak (10); wine = 5 common + 9 wine-only + 5
#: oak (19 — D-87 added the four thermal Strecker/sotolon pools to the wine-only set, D-94 added
#: furaneol to the oak set). Wine's order keeps the oak five last, so the D-87 pools slot in after
#: the two D-75 aldehydes, before oak.
AROMA_COMPOUNDS: Mapping[str, tuple[AromaCompound, ...]] = {
    "beer": _COMMON + _OAK,
    "wine": _COMMON + _WINE_ONLY + _OAK,
}


def load_thresholds() -> ParameterSet:
    """Load ``sensory.yaml`` standalone — the perception thresholds, provenance intact.

    Deliberately **not** routed through the compile seam's ``shared_files`` (no RHS reads a
    threshold), so these numbers never enter any ``CompiledScenario.param_values`` and cannot
    perturb the chemistry. Returns a :class:`ParameterSet` (each threshold in µg/L, all
    ``speculative``). NB (D-24): because they load standalone, thresholds sit **outside** the
    ensemble parameter sweep — ``simulate_ensemble`` samples only compiled-scenario params, so
    it does not propagate threshold uncertainty into an OAV band. Defensible for a readout
    already floored at speculative; stated here so it reads as a choice, not an oversight.
    """
    return load_parameters(default_data_dir() / "sensory.yaml")


def medium_of(schema: StateSchema) -> str:
    """Infer ``"beer"`` or ``"wine"`` from a state schema's signature slots.

    ``iso_alpha`` is beer-only, ``tartaric`` is wine-only (both present in their schema even
    at default 0), so they are unambiguous fingerprints. Raises if neither is present (some
    other/bare medium the sensory sets are not defined for).
    """
    names = set(schema.names)
    if "iso_alpha" in names:
        return "beer"
    if "tartaric" in names:
        return "wine"
    raise ValueError(
        "cannot infer medium for the sensory readout: schema has neither the beer-only "
        "'iso_alpha' nor the wine-only 'tartaric' signature slot"
    )


def _threshold_key(pool: str, medium: str) -> str:
    """The matrix-specific threshold parameter name, e.g. ``threshold_diacetyl_beer``."""
    return f"threshold_{pool}_{medium}"


def _compound(pool: str, medium: str) -> AromaCompound:
    """The :class:`AromaCompound` for ``pool`` in ``medium``, or a clear error.

    A wine-only pool requested on a beer trajectory raises here (rather than a downstream
    ``KeyError`` on the missing state slot) — the beer schema has no such slot at all.
    """
    for c in AROMA_COMPOUNDS[medium]:
        if c.pool == pool:
            return c
    available = ", ".join(c.pool for c in AROMA_COMPOUNDS[medium])
    raise ValueError(
        f"{pool!r} is not an aroma-active pool for medium {medium!r}; have: {available}"
    )


def oav_tier(input_tier: Tier, threshold_tier: Tier) -> Tier:
    """The tier floor for an OAV: ``combine(input_tier, threshold_tier, SPECULATIVE)``.

    Always :attr:`Tier.SPECULATIVE`, by construction — that is the point. The explicit
    ``SPECULATIVE`` term is **not** redundant with ``threshold_tier``: it encodes that the
    OAV *mapping itself* is speculative (a concentration ratio is only a proxy for perceived
    aroma), so the floor holds even for a validated input concentration and even if a
    threshold were later (mis)labelled plausible or validated. This is the §4.3 credibility
    firewall made executable — speculation may never borrow the validated core's confidence.
    """
    return combine([input_tier, threshold_tier, Tier.SPECULATIVE])


def oav_series(traj: Trajectory, thresholds: ParameterSet, pool: str) -> FloatArray:
    """OAV of one aroma pool at each stored time — dimensionless ``conc / threshold``.

    Selects the matrix-specific threshold for ``traj``'s medium, crosses it from µg/L to the
    canonical g/L, and divides the pool's g/L concentration by it. Monotone increasing in the
    pool; identically 0 when the pool is 0 (a clean run raises no false aroma). Raises for a
    pool not aroma-active in this medium (e.g. a wine-only pool on a beer trajectory).
    """
    medium = medium_of(traj.schema)
    _compound(pool, medium)  # validates membership; raises a clear error otherwise
    threshold_gpl = ugl_to_gpl(thresholds.value(_threshold_key(pool, medium)))
    conc_gpl = np.asarray(traj.series(pool), dtype=np.float64)
    return conc_gpl / threshold_gpl


@dataclass(frozen=True)
class OAVReading:
    """One aroma pool's OAV at a single time, with its descriptor, flag, and tier floor."""

    compound: str
    representative: str
    descriptor: str
    oav: float
    above_threshold: bool
    lumped: bool
    tier: Tier


@dataclass(frozen=True)
class SensoryProfile:
    """The aroma profile of a trajectory at one time — per-compound OAVs, not a summed scalar.

    Deliberately **not** a single aggregate number: summing OAVs assumes perceptual
    additivity, which is contested and would over-claim. Instead this reports each pool's
    :class:`OAVReading` and lets a caller ask which pools clear their threshold
    (:meth:`above_threshold`). Descriptor-space projection ("this smells like leather and
    banana") is the deferred, even-more-speculative beat 1b.
    """

    medium: str
    time_index: int
    readings: Mapping[str, OAVReading]

    def __iter__(self) -> Iterator[OAVReading]:
        return iter(self.readings.values())

    def above_threshold(self) -> list[str]:
        """Pools whose OAV exceeds 1 at this time (individually above perception)."""
        return [name for name, r in self.readings.items() if r.above_threshold]

    def tier(self) -> Tier:
        """The profile's tier — the floor across its readings, i.e. always speculative."""
        return combine(r.tier for r in self.readings.values())


def sensory_profile(
    traj: Trajectory, thresholds: ParameterSet, *, time_index: int = -1
) -> SensoryProfile:
    """Build the aroma :class:`SensoryProfile` at ``time_index`` (default: the finished state).

    Covers exactly the medium's aroma set (10 pools for beer, 19 for wine — D-94) — so the reported
    compounds *are* the medium's, never a wine-only pool on beer. Each reading's tier is the
    :func:`oav_tier` floor over the pool's chemistry tier (read from ``traj.tier_map``) and
    the threshold's tier; both fold under the mandatory speculative floor.
    """
    medium = medium_of(traj.schema)
    readings: dict[str, OAVReading] = {}
    for c in AROMA_COMPOUNDS[medium]:
        series = oav_series(traj, thresholds, c.pool)
        oav = float(series[time_index])
        input_tier = traj.tier_map.get(c.pool, Tier.VALIDATED)
        threshold_tier = thresholds.tier_of(_threshold_key(c.pool, medium))
        readings[c.pool] = OAVReading(
            compound=c.pool,
            representative=c.representative,
            descriptor=c.descriptor,
            oav=oav,
            above_threshold=oav > 1.0,
            lumped=c.lumped,
            tier=oav_tier(input_tier, threshold_tier),
        )
    return SensoryProfile(medium=medium, time_index=time_index, readings=readings)
