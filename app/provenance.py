"""From a line on a chart back to the papers underneath it.

Plenty of software draws a fermentation curve. What this engine can do that they cannot is
say, for any point on that line, which chemical steps produced it, which published numbers
those steps used, what conditions each number was measured under, how wide a range the paper
gave it, and therefore how much the line can be trusted. All of that is already stored — a
number cannot even be loaded into the engine without it — so this module is mostly assembly.

The trail runs in four steps, one function each:

1. :func:`mechanisms_for` — which chemical steps change this quantity, and which other steps
   speed them up or slow them down.
2. :func:`constants_for` — every published number those steps use.
3. :func:`card` — one number with its value, units, range, confidence and full citation.
4. :func:`limiting` — the numbers holding the confidence down, which is the answer to
   "why can I not trust this more?"
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from app.runner import MechanismInfo, RunResult
from fermentation.core.tiers import Tier


@dataclass(frozen=True)
class ConstantCard:
    """One parameter, flattened for display, with its whole provenance intact."""

    name: str
    value: float
    unit: str
    tier: Tier
    low: float
    high: float
    range_note: str
    source: str
    conditions: str
    doi: str | None
    notes: str
    #: Mechanisms in this run that read it.
    read_by: tuple[str, ...]

    @property
    def is_estimate(self) -> bool:
        """True when the source is an author's own guess rather than published work."""
        return "author estimate" in self.source.lower()

    @property
    def span_fraction(self) -> float | None:
        """Width of the stated range as a fraction of the value. ``None`` if undefined."""
        if self.value == 0.0:
            return None
        return abs(self.high - self.low) / abs(self.value)


def mechanisms_for(result: RunResult, variable: str) -> tuple[MechanismInfo, ...]:
    """The steps that change ``variable``, then the ones that scale those steps' rates."""
    return result.touching(variable)


def constants_for(
    result: RunResult, mechanisms: Sequence[MechanismInfo]
) -> tuple[ConstantCard, ...]:
    """Every number those steps use, weakest first.

    Sorted by confidence and then by name, so whatever is holding the answer down sits at the
    top of the list instead of being buried alphabetically.
    """
    read_by: dict[str, list[str]] = {}
    for m in mechanisms:
        for name in m.reads:
            read_by.setdefault(name, []).append(m.name)

    cards: list[ConstantCard] = []
    for name, readers in read_by.items():
        if name not in result.parameters:
            continue
        cards.append(card(result, name, tuple(sorted(readers))))
    cards.sort(key=lambda c: (int(c.tier), c.name))
    return tuple(cards)


def card(result: RunResult, name: str, read_by: tuple[str, ...] = ()) -> ConstantCard:
    """One number, with everything its data file was required to state about it."""
    p = result.parameters[name]
    return ConstantCard(
        name=p.name,
        value=p.value,
        unit=p.unit,
        tier=p.tier,
        low=p.uncertainty.low,
        high=p.uncertainty.high,
        range_note=p.uncertainty.note,
        source=p.provenance.source,
        conditions=p.provenance.conditions,
        doi=p.provenance.doi,
        notes=p.provenance.notes,
        read_by=read_by,
    )


def limiting(cards: Sequence[ConstantCard]) -> tuple[ConstantCard, ...]:
    """The numbers sitting at the weakest mark present — the ones holding the answer down.

    This is the half of the trail you can act on. A result marked speculative is not vaguely
    speculative: it is speculative because of *these* specific numbers, and finding a source
    for any one of them is a concrete job with a known payoff.
    """
    if not cards:
        return ()
    worst = min(int(c.tier) for c in cards)
    return tuple(c for c in cards if int(c.tier) == worst)


def why_this_tier(result: RunResult, variable: str) -> str:
    """One plain paragraph answering "how much can I trust this line, and why not more?"."""
    tier = result.traj.tier_map.get(variable)
    if variable not in result.touched_variables():
        return (
            f"Nothing in this run changes {variable}. It held whatever it was set to at the "
            "start, right through to the end, so the model has not actually made a claim "
            "about it. That is not the same as the value being confirmed correct."
        )
    mechs = mechanisms_for(result, variable)
    cards = constants_for(result, mechs)
    caps = limiting(cards)
    steps = [m for m in mechs if m.kind == "process"]
    scalers = [m for m in mechs if m.kind == "modifier"]

    parts = [
        f"{variable} is produced or consumed by {len(steps)} chemical step(s) in the model"
        + (f", with {len(scalers)} more that speed those up or slow them down" if scalers else "")
        + f". Between them they use {len(cards)} published numbers."
    ]
    if tier is not None:
        parts.append(WHY_TIER_SENTENCE[tier])
    if caps:
        names = ", ".join(c.name for c in caps[:5])
        more = f", and {len(caps) - 5} others" if len(caps) > 5 else ""
        parts.append(
            f"{len(caps)} of those numbers are the weakest link: {names}{more}. Finding a "
            "published measurement for any one of them is what would raise the mark."
        )
    estimates = [c for c in cards if c.is_estimate]
    if estimates:
        parts.append(
            f"{len(estimates)} of the {len(cards)} are the author's own estimates rather than "
            "measurements taken from a paper."
        )
    return " ".join(parts)


#: How the reported confidence is explained, per mark, in the paragraph above.
WHY_TIER_SENTENCE: dict[Tier, str] = {
    Tier.SPECULATIVE: (
        "The result is marked speculative, meaning at least one of those numbers is an "
        "estimate rather than a measurement."
    ),
    Tier.PLAUSIBLE: (
        "The result is marked plausible: every number behind it is published, but the result "
        "itself has never been compared against a real fermentation."
    ),
    Tier.VALIDATED: (
        "The result is marked validated, meaning it has been checked against real measured data."
    ),
}


def tier_census(result: RunResult) -> dict[Tier, int]:
    """How many of this run's numbers sit at each confidence mark."""
    counts: dict[Tier, int] = dict.fromkeys(Tier, 0)
    for name in result.parameters.names:
        counts[result.parameters.tier_of(name)] += 1
    return counts


def estimate_census(result: RunResult) -> tuple[int, int]:
    """``(numbers that are the author's estimates, total numbers)`` for this run."""
    total = 0
    estimates = 0
    for name in result.parameters.names:
        total += 1
        if "author estimate" in result.parameters[name].provenance.source.lower():
            estimates += 1
    return estimates, total
