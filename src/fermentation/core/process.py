"""The Process abstraction — the core's central idea.

A **Process** is anything that contributes to the time derivative of the state:
primary alcoholic fermentation, malolactic fermentation, oxidation, oak
extraction. The total derivative is the sum of the active Processes'
contributions. This compositionality is the key design bet — we model mechanisms
and let combinations emerge, rather than scripting outcomes — and it means a
speculative Process can be toggled off while the validated core still runs and
still passes its tests.

Each Process:
  * reads the current state and resolved parameters,
  * returns its contribution to ``d(state)/dt`` (same shape as the state array),
  * declares which state variables it touches and which :class:`Tier` it is,
  * can be individually enabled, disabled, or swapped (via :class:`ProcessSet`).

Parameters reach a Process as a plain ``Mapping[str, float]`` of *resolved*
numeric values. Provenance and tier metadata are handled at setup/analysis time
(see ``fermentation.parameters``), not inside this hot loop.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence

from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier, combine


class Process(ABC):
    """Base class for anything that contributes to ``d(state)/dt``.

    Subclasses set :attr:`name`, :attr:`tier`, and :attr:`touches`, and implement
    :meth:`derivatives`.
    """

    #: Stable identifier, unique within a :class:`ProcessSet`.
    name: str
    #: Confidence tier this Process belongs to.
    tier: Tier
    #: State variable names this Process is allowed to contribute to.
    touches: tuple[str, ...]

    @abstractmethod
    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        """Return this Process's contribution to ``d(state)/dt``.

        The returned array has the same shape as ``y``. Entries for variables not
        in :attr:`touches` must be zero — :class:`ProcessSet` can enforce this in
        ``strict`` mode.
        """
        raise NotImplementedError

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name!r}, tier={self.tier.label})"


class ProcessSet:
    """A collection of Processes that sums to the system's total derivative.

    Owns enable/disable state and derives per-variable output tiers from the
    active Processes that touch each variable. The ``total_derivatives`` callable
    is what the runtime hands to ``solve_ivp``.
    """

    def __init__(
        self,
        schema: StateSchema,
        processes: Sequence[Process],
        *,
        strict: bool = False,
    ) -> None:
        names = [p.name for p in processes]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate process names: {names}")
        for p in processes:
            unknown = set(p.touches) - set(schema.names)
            if unknown:
                raise ValueError(
                    f"Process {p.name!r} touches unknown variables {sorted(unknown)}; "
                    f"schema has {schema.names}"
                )
        self.schema = schema
        self._processes: dict[str, Process] = {p.name: p for p in processes}
        self._enabled: dict[str, bool] = {p.name: True for p in processes}
        #: When True, every contribution is checked against its ``touches`` set.
        self.strict = strict

    # -- membership / toggling ------------------------------------------------

    def enable(self, name: str) -> None:
        self._require(name)
        self._enabled[name] = True

    def disable(self, name: str) -> None:
        self._require(name)
        self._enabled[name] = False

    def is_enabled(self, name: str) -> bool:
        self._require(name)
        return self._enabled[name]

    @property
    def active(self) -> tuple[Process, ...]:
        return tuple(p for n, p in self._processes.items() if self._enabled[n])

    def _require(self, name: str) -> None:
        if name not in self._processes:
            raise KeyError(f"No process named {name!r}; have {list(self._processes)}")

    # -- the hot loop ---------------------------------------------------------

    def total_derivatives(self, t: float, y: FloatArray, params: Mapping[str, float]) -> FloatArray:
        """Sum the active Processes' contributions to ``d(state)/dt``."""
        total = self.schema.zeros()
        for p in self.active:
            contribution = p.derivatives(t, y, self.schema, params)
            if self.strict:
                self._check_touch_contract(p, contribution)
            total += contribution
        return total

    def _check_touch_contract(self, p: Process, contribution: FloatArray) -> None:
        if contribution.shape != (self.schema.size,):
            raise ValueError(
                f"Process {p.name!r} returned shape {contribution.shape}, "
                f"expected ({self.schema.size},)"
            )
        allowed = self.schema.zeros().astype(bool)
        for var in p.touches:
            allowed[self.schema.slice(var)] = True
        leaked = (contribution != 0.0) & ~allowed
        if leaked.any():
            bad = [n for n in self.schema.names if (contribution[self.schema.slice(n)] != 0).any()]
            raise ValueError(
                f"Process {p.name!r} contributed to undeclared variables; "
                f"nonzero on {bad}, declares touches={list(p.touches)}"
            )

    # -- tier propagation -----------------------------------------------------

    def tier_of(self, variable: str) -> Tier:
        """Derived output tier of ``variable``: the lowest tier among the active
        Processes that touch it. A variable nothing touches is ``VALIDATED`` (it
        is simply unchanging, which is trivially correct).

        NOTE (Milestone 1): this propagates *Process* tiers only. A Process must
        also be capped by the tiers of the parameters it consumes — otherwise a
        VALIDATED process running on speculative placeholder parameters would
        report a VALIDATED output, which is exactly the credibility-borrowing the
        tier system exists to prevent. Parameter-tier propagation is wired in when
        real Processes (which declare the parameters they read) arrive. See
        docs/DECISIONS.md D-1 and milestone-1-tasks.md.
        """
        if variable not in self.schema:
            raise KeyError(f"Unknown variable {variable!r}")
        contributing = [p.tier for p in self.active if variable in p.touches]
        return combine(contributing)

    def tier_map(self) -> dict[str, Tier]:
        """Derived tier for every state variable."""
        return {name: self.tier_of(name) for name in self.schema.names}

    def overall_tier(self) -> Tier:
        """Lowest tier across all active Processes — the tier of any output that
        mixes the whole state (e.g. a full time-series export)."""
        return combine([p.tier for p in self.active])
