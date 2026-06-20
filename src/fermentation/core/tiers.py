"""Confidence tiers — the project's prime directive made executable.

Every modelled quantity belongs to one of three tiers. The tier is metadata that
must travel with a value all the way to any output: the engine must never silently
blend a ``VALIDATED`` concentration with a ``SPECULATIVE`` one and present the
result as equally trustworthy.

Tier is a property of the *Processes* and *parameters* that produce a value, not
of the raw float that flows through the ODE solver. An output's tier is therefore
*derived* — it is the lowest (least trustworthy) tier of everything that fed it,
via :func:`combine`.
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import IntEnum


class Tier(IntEnum):
    """Confidence in a modelled quantity, ordered most → least trustworthy.

    Ordering matters: ``VALIDATED > PLAUSIBLE > SPECULATIVE``. A higher integer
    means *more* trustworthy, so the trustworthiness of a combination is the
    ``min`` of its parts (see :func:`combine`).
    """

    SPECULATIVE = 0
    """Real chemistry, but integrating it into a trustworthy prediction is not
    solved science (aging verdicts, sensory mapping). Isolate and label loudly."""

    PLAUSIBLE = 1
    """Sound mechanism with literature support, but harder to validate
    quantitatively. Validated directionally/qualitatively at best."""

    VALIDATED = 2
    """Established, published science that we check against benchmark curves."""

    @property
    def label(self) -> str:
        """Human-facing lowercase label, e.g. for plots and reports."""
        return self.name.lower()


def combine(tiers: Iterable[Tier]) -> Tier:
    """Return the lowest (least trustworthy) tier among ``tiers``.

    This is the rule for confidence propagation: an output is only as trustworthy
    as its weakest input. An empty input is treated as ``VALIDATED`` (the
    identity for ``min`` over this lattice) so that combining "nothing" with a
    value leaves that value's tier unchanged.

    >>> combine([Tier.VALIDATED, Tier.SPECULATIVE])
    <Tier.SPECULATIVE: 0>
    >>> combine([])
    <Tier.VALIDATED: 2>
    """
    return min(tiers, default=Tier.VALIDATED)
