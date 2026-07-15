"""Descriptor-space projection — the OAV vector grouped into aroma vocabulary (decision D-95).

Beat **1b** of the handoff §4.2 sensory layer, deferred at D-66 and built here as slice 1 of
two. Where :mod:`fermentation.sensory.oav` (beat 1a, D-67) answers *"how many times over its
perception threshold is compound i?"*, this module answers *"which descriptor words does that
OAV vector light up, and which compound is driving each?"* — a **projection** of the 19 (wine)
/ 10 (beer) aroma pools onto 14 / 9 descriptor axes.

It consumes a :class:`~fermentation.sensory.oav.SensoryProfile` and returns a
:class:`DescriptorProfile`. Nothing in ``core``/``runtime``/``scenario`` imports it — the §4.2
cardinal rule (*the sensory layer consumes the chemistry; the chemistry never depends on the
sensory layer*) holds unchanged, one layer further up. It adds **no state, no Process, no
ledger entry, and no parameters** — see "no magic numbers" below.

**The additivity through-line — the load-bearing call.** Descriptor projection is inherently
many-to-many: ``malty`` collects three branched-chain aldehyde pools, ``smoky`` collects both
oak ``guaiacol`` and Brett ``ethylguaiacols``. But the layer directly below refused to
aggregate: :class:`~fermentation.sensory.oav.SensoryProfile` reports per-compound OAVs and
*never* a summed scalar, because summing assumes perceptual additivity — contested, and it
would over-claim (D-67). So the aggregation rule here is **not a free choice**. Summing OAVs
per descriptor would silently reintroduce the exact assumption the layer beneath rejected.
:class:`MaxRuleProjector` therefore takes the **max**, which asserts no additivity: it names
the loudest contributor and reports its OAV, nothing more. *We never assume additivity, at
any layer.* Weighting/compression/masking — the perceptual math that genuinely needs
parameters — is the explicitly deferred slice 2.

**Honest framing: what this layer does and does not add.** Under the max rule a descriptor
clears threshold **iff** one of its pools does, so :meth:`DescriptorProfile.above_threshold`
is a pure *regrouping* of beat 1a's pool-level flags — it carries no new above-threshold
information. Slice 1 delivers **vocabulary grouping + dominant attribution**, not a new
sensory claim. Making the per-descriptor number say something a regrouping cannot is exactly
what slice 2's weighting buys, and is why slice 2 is where the speculation lives.

**The swappable seam (the handoff's own requirement).** §4.2 asks for "a separate, swappable
sensory model ... with a clean seam so it can later be replaced by an ML model trained on real
sensory-panel data." :class:`DescriptorProjector` is that seam: any object with a
``project(SensoryProfile) -> DescriptorProfile`` method drops in. :class:`MaxRuleProjector` is
merely the one concrete v1 rule.

**No magic numbers — membership is structure, not parameters.** The pool→descriptor map is
*binary* (a pool either feeds an axis or does not), which pairs naturally with the max rule
and needs no weights. So slice 1 mints no constants and adds no YAML: it lives in code exactly
as :attr:`~fermentation.sensory.oav.AromaCompound.descriptor` already does (accepted at D-67),
and :class:`~fermentation.parameters.schema.Parameter` forbids extra fields anyway, so a
``descriptor:`` key could not join ``sensory.yaml``. Weights are precisely what will make
slice 2 need a provenance file (the ``thermal.yaml`` relative-weight precedent, D-87).

**Vocabulary vs. gloss (an accepted, documented redundancy).**
:attr:`~fermentation.sensory.oav.AromaCompound.descriptor` is a per-compound *prose gloss*
("banana / fruity"); :attr:`DescriptorAxis.name` is the machine-readable *vocabulary*
(``fruity``). They coexist and could drift. Deriving one from the other would mean parsing
slash-phrases (brittle) or churning beat 1a for no gain, so the redundancy is accepted and
recorded rather than tested.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from fermentation.core.tiers import Tier, combine
from fermentation.sensory.oav import AROMA_COMPOUNDS, SensoryProfile


@dataclass(frozen=True)
class DescriptorAxis:
    """One descriptor-vocabulary word and the aroma pools that feed it.

    ``pools`` is **binary many-to-many** membership: a pool either contributes to this axis or
    does not (no weights — see the module docstring's "no magic numbers"). A pool may feed
    several axes (``ethylguaiacols`` is both ``smoky`` and ``clove_spice`` — 4-EG genuinely
    smells of both), and an axis may collect several pools (``malty``). ``gloss`` is a human
    sentence for reports; it is documentation, never logic.

    Membership is declared **medium-agnostically**. A medium's axis set is *derived* by
    intersecting each axis with that medium's pool set (:func:`axes_for_medium`), so an axis
    with no pools in a medium simply does not exist there — beer can never report
    ``barnyard``, by construction rather than by a second hand-maintained list.
    """

    name: str
    pools: tuple[str, ...]
    gloss: str


#: The descriptor vocabulary — 14 axes over wine's 19 aroma pools, 9 over beer's 10 (D-95).
#:
#: Granularity is the owner's call (~12 many-to-many axes): coarse enough that the projection
#: does real work (``malty`` collapses three pools), fine enough to keep distinctions the
#: chemistry worked to earn. Two calls worth recording:
#:
#: * ``green_apple`` stays its own axis rather than joining ``cooked_potato`` under a shared
#:   "oxidative" — ``acetaldehyde`` is a fermentation intermediate (D-27) long before it is an
#:   oxidation product (D-71), so a young beer's green-apple note would be mislabelled.
#: * ``caramel`` (furaneol) and ``curry_maple`` (sotolon) stay split — toffee and curry are
#:   different smells, and D-94 kept those two compounds' descriptors deliberately distinct.
#:   NB this is a judgement, not a constraint: ``smoky`` merges ``guaiacol`` +
#:   ``ethylguaiacols``, whose distinctness the codebase flags just as loudly, so distinctness
#:   at the compound layer cannot by itself forbid a merge at this one. Collapsing distinct
#:   compounds is what a projection *is*; merging these two back into one ``nutty_caramel``
#:   axis would be equally defensible and costs exactly this one entry.
DESCRIPTOR_AXES: tuple[DescriptorAxis, ...] = (
    DescriptorAxis("fruity", ("esters",), "banana / pear-drop / general fruit esters"),
    DescriptorAxis("solventy", ("fusels",), "hot, solventy higher alcohols"),
    DescriptorAxis("buttery", ("diacetyl",), "butter / butterscotch"),
    DescriptorAxis("green_apple", ("acetaldehyde",), "green apple / bruised fruit"),
    DescriptorAxis("sulfidic", ("h2s", "mercaptans"), "rotten egg / drains / reductive"),
    DescriptorAxis("vanilla_oak", ("vanillin", "whiskey_lactone"), "vanilla / coconut / sweet oak"),
    DescriptorAxis("smoky", ("guaiacol", "ethylguaiacols"), "smoke / toast / campfire"),
    DescriptorAxis("clove_spice", ("eugenol", "ethylguaiacols"), "clove / warm spice"),
    DescriptorAxis("caramel", ("furaneol",), "caramel / toffee / burnt sugar"),
    DescriptorAxis("barnyard", ("ethylphenols",), "horse-sweat / barnyard / band-aid"),
    DescriptorAxis("floral_honey", ("phenylacetaldehyde",), "honey / rose / floral"),
    DescriptorAxis("cooked_potato", ("methional",), "cooked potato / oxidised savoury"),
    DescriptorAxis(
        "malty",
        ("2_methylbutanal", "3_methylbutanal", "2_methylpropanal"),
        "malt / almond / dark chocolate / grain",
    ),
    DescriptorAxis("curry_maple", ("sotolon",), "curry / maple / walnut"),
)


def axes_for_medium(medium: str) -> tuple[DescriptorAxis, ...]:
    """The axes that exist for ``medium``, each narrowed to that medium's pools.

    Derived from :data:`~fermentation.sensory.oav.AROMA_COMPOUNDS` rather than declared, so
    the vocabulary can never drift from the aroma set: an axis whose pools are all absent from
    the medium is dropped entirely (beer has no ``barnyard``), and an axis that keeps only
    some of its pools is narrowed (beer's ``sulfidic`` is ``h2s`` alone — ``mercaptans`` is
    wine-only; beer's ``smoky`` is oak ``guaiacol`` alone — Brett ``ethylguaiacols`` is
    wine-only). Wine → 14 axes, beer → 9 (a strict subset).
    """
    pools = {c.pool for c in AROMA_COMPOUNDS[medium]}
    narrowed = (
        DescriptorAxis(a.name, tuple(p for p in a.pools if p in pools), a.gloss)
        for a in DESCRIPTOR_AXES
    )
    return tuple(a for a in narrowed if a.pools)


def descriptor_tier(contributing: Iterable[Tier]) -> Tier:
    """The tier floor for a descriptor: ``combine(*contributing_oav_tiers, SPECULATIVE)``.

    Always :attr:`Tier.SPECULATIVE`, by construction — that is the point, and it repeats
    D-67's argument one layer up. The explicit ``SPECULATIVE`` is **not** redundant with the
    incoming OAV tiers (themselves already floored): it encodes that the *projection itself*
    is a further heuristic leap beyond the sourced ratio — grouping compounds under a word,
    and naming one of them dominant, is a claim about perception that no threshold measurement
    backs. The floor therefore holds even if an OAV were later (mis)labelled plausible.
    """
    return combine([*contributing, Tier.SPECULATIVE])


@dataclass(frozen=True)
class DescriptorReading:
    """One descriptor axis at a single time — its contributors, the loudest, and the floor.

    ``oav`` is the **maximum** over ``contributors``, never the sum (see the module docstring's
    additivity through-line): it is ``dominant``'s OAV, i.e. "the loudest thing making this
    smell malty, and how far over threshold it is". ``lumped`` propagates D-66's
    fixed-lump-composition honesty cost from the *dominant* contributor — a descriptor driven
    by a lumped pool inherits that pool's assumption, and the caveat must not vanish just
    because it crossed a layer.
    """

    descriptor: str
    contributors: Mapping[str, float]
    dominant: str
    oav: float
    above_threshold: bool
    lumped: bool
    tier: Tier


@dataclass(frozen=True)
class DescriptorProfile:
    """A trajectory's aroma profile in descriptor space at one time.

    The vocabulary-level counterpart of :class:`~fermentation.sensory.oav.SensoryProfile`, and
    like it deliberately **not** a single aggregate: no "total aroma intensity" scalar, no
    summed axis. Each axis reports its own :class:`DescriptorReading`.
    """

    medium: str
    time_index: int
    readings: Mapping[str, DescriptorReading]

    def __iter__(self) -> Iterator[DescriptorReading]:
        return iter(self.readings.values())

    def above_threshold(self) -> list[str]:
        """Descriptors with at least one pool over its perception threshold at this time.

        Under the max rule this is a pure **regrouping** of beat 1a's pool-level flags (a
        descriptor clears iff one of its pools does) — it adds vocabulary and attribution, not
        new above-threshold information. See the module docstring's honest framing.
        """
        return [name for name, r in self.readings.items() if r.above_threshold]

    def dominant_pools(self) -> dict[str, str]:
        """Descriptor → the pool driving it (the argmax contributor)."""
        return {name: r.dominant for name, r in self.readings.items()}

    def tier(self) -> Tier:
        """The profile's tier — the floor across its readings, i.e. always speculative."""
        return combine(r.tier for r in self.readings.values())


@runtime_checkable
class DescriptorProjector(Protocol):
    """The swappable seam: OAV vector → descriptor space (handoff §4.2).

    §4.2 asks for the projection to be "a separate, swappable sensory model ... with a clean
    seam so it can later be replaced by an ML model trained on real sensory-panel data should
    such data become available". This Protocol *is* that seam, and it is deliberately the
    narrowest thing that can be: consume a :class:`~fermentation.sensory.oav.SensoryProfile`,
    return a :class:`DescriptorProfile`. A future panel-trained model — which would emit
    genuine per-descriptor intensities rather than a max — satisfies it without touching one
    line of the chemistry, of beat 1a, or of any caller.
    """

    def project(self, profile: SensoryProfile) -> DescriptorProfile:
        """Project one OAV profile onto descriptor space."""
        ...


class MaxRuleProjector:
    """The v1 projector: each descriptor reads its **loudest** contributing pool's OAV.

    The deliberately-modest rule, and the only one consistent with the layer below (see the
    module docstring): ``max`` asserts no perceptual additivity, so no descriptor can claim an
    intensity that no single compound justifies — three pools at OAV 0.4 read 0.4, not 1.2,
    and the profile stays silent rather than inventing a below-threshold smell from a sum.
    Weighted/compressed/masked intensities are slice 2.

    Ties break toward the pool listed first in :attr:`DescriptorAxis.pools` (dict insertion
    order follows it), so ``dominant`` is deterministic.
    """

    def project(self, profile: SensoryProfile) -> DescriptorProfile:
        readings: dict[str, DescriptorReading] = {}
        for axis in axes_for_medium(profile.medium):
            contributors = {pool: profile.readings[pool].oav for pool in axis.pools}
            dominant = max(contributors, key=lambda p: contributors[p])
            oav = contributors[dominant]
            readings[axis.name] = DescriptorReading(
                descriptor=axis.name,
                contributors=contributors,
                dominant=dominant,
                oav=oav,
                above_threshold=oav > 1.0,
                lumped=profile.readings[dominant].lumped,
                tier=descriptor_tier(profile.readings[p].tier for p in axis.pools),
            )
        return DescriptorProfile(
            medium=profile.medium, time_index=profile.time_index, readings=readings
        )
