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

**Rate modifiers (multiplicative, not additive).** Some mechanisms do not *add* a
flux — they *scale* an existing one: ethanol inhibition slows fermentative uptake,
temperature (Arrhenius) scales every rate constant. Because :class:`ProcessSet`
sums Processes, such a mechanism cannot be a summed Process (it would *add* to a
derivative, not multiply it). A :class:`RateModifier` instead returns a scalar
factor that :class:`ProcessSet` multiplies onto the *whole contribution vector* of
the Process(es) it targets, before summing. Scaling a conserving Process's entire
contribution by one scalar preserves its mass/atom balances, so a modifier never
breaks conservation. Modifiers are togglable and tier-tracked exactly like
Processes. (See decision D-10.)
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


class RateModifier(ABC):
    """A multiplicative scaling applied to one or more Processes' contributions.

    Where a :class:`Process` *adds* a flux, a ``RateModifier`` *scales* one: it
    returns a dimensionless :meth:`factor` that :class:`ProcessSet` multiplies onto
    the entire contribution vector of every Process named in :attr:`modifies`. This
    is how mechanisms that act *on* a rate rather than alongside it — ethanol
    inhibition (``factor`` in ``[0, 1]``), Arrhenius temperature dependence
    (``factor`` may exceed 1 above the reference temperature) — compose into the
    additive :class:`ProcessSet` without breaking it.

    Because the factor scales a conserving Process's *whole* contribution by one
    scalar, every balance that Process respects is preserved (a uniform slow-down
    of a carbon-neutral flux is still carbon-neutral). A modifier must therefore
    return a non-negative factor; a wall-type form (e.g. ``1 - E/E_max``) must
    clamp at zero so a state overshoot cannot flip the factor negative and *create*
    substrate. Modifiers are enabled/disabled and contribute to derived tiers
    exactly like Processes (see :class:`ProcessSet`).
    """

    #: Stable identifier, unique across the Processes *and* modifiers of a set.
    name: str
    #: Confidence tier this modifier belongs to.
    tier: Tier
    #: Names of the Processes whose contribution this modifier scales.
    modifies: tuple[str, ...]
    #: Parameters this modifier reads (for parameter-tier propagation, D-1).
    reads: tuple[str, ...] = ()

    @abstractmethod
    def factor(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> float:
        """Return the non-negative scalar multiplied onto each targeted Process."""
        raise NotImplementedError

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(name={self.name!r}, tier={self.tier.label}, "
            f"modifies={list(self.modifies)})"
        )


class ProcessSet:
    """A collection of Processes that sums to the system's total derivative.

    Owns enable/disable state and derives per-variable output tiers from the
    active Processes that touch each variable. The ``total_derivatives`` callable
    is what the runtime hands to ``solve_ivp``.

    Optionally holds :class:`RateModifier` objects. A modifier scales (multiplies)
    the contribution of every Process it names in :attr:`RateModifier.modifies`,
    applied per right-hand-side evaluation before the sum. Modifiers share the
    Processes' name space (names must be unique across both) and the same
    enable/disable machinery, so a speculative modifier toggles off cleanly. A
    disabled modifier contributes a factor of 1 and is excluded from tier
    derivation.
    """

    def __init__(
        self,
        schema: StateSchema,
        processes: Sequence[Process],
        *,
        modifiers: Sequence[RateModifier] = (),
        strict: bool = False,
    ) -> None:
        names = [p.name for p in processes]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate process names: {names}")
        mod_names = [m.name for m in modifiers]
        if len(mod_names) != len(set(mod_names)):
            raise ValueError(f"Duplicate modifier names: {mod_names}")
        clash = set(names) & set(mod_names)
        if clash:
            raise ValueError(f"Modifier name(s) clash with process names: {sorted(clash)}")
        for p in processes:
            unknown = set(p.touches) - set(schema.names)
            if unknown:
                raise ValueError(
                    f"Process {p.name!r} touches unknown variables {sorted(unknown)}; "
                    f"schema has {schema.names}"
                )
        for m in modifiers:
            missing = set(m.modifies) - set(names)
            if missing:
                raise ValueError(
                    f"Modifier {m.name!r} targets unknown process(es) {sorted(missing)}; "
                    f"known processes: {sorted(names)}"
                )
        self.schema = schema
        self._processes: dict[str, Process] = {p.name: p for p in processes}
        self._modifiers: dict[str, RateModifier] = {m.name: m for m in modifiers}
        self._enabled: dict[str, bool] = {
            **{p.name: True for p in processes},
            **{m.name: True for m in modifiers},
        }
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

    @property
    def active_modifiers(self) -> tuple[RateModifier, ...]:
        return tuple(m for n, m in self._modifiers.items() if self._enabled[n])

    def _require(self, name: str) -> None:
        if name not in self._enabled:
            raise KeyError(f"No process or modifier named {name!r}; have {list(self._enabled)}")

    # -- the hot loop ---------------------------------------------------------

    def total_derivatives(self, t: float, y: FloatArray, params: Mapping[str, float]) -> FloatArray:
        """Sum the active Processes' contributions to ``d(state)/dt``.

        Each active modifier's factor is evaluated once, then multiplied onto the
        contribution of every Process it targets before the sum. With no active
        modifiers this is exactly the plain additive path (no extra array ops).
        """
        total = self.schema.zeros()
        mods = self.active_modifiers
        factors = {m.name: m.factor(t, y, self.schema, params) for m in mods}
        for p in self.active:
            contribution = p.derivatives(t, y, self.schema, params)
            if self.strict:
                self._check_touch_contract(p, contribution)
            for m in mods:
                if p.name in m.modifies:
                    # Scale the whole vector by one scalar — preserves the
                    # Process's conservation balances (D-10). Strict checks the
                    # raw contribution above; scaling zeros leaves them zero.
                    contribution = contribution * factors[m.name]
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
        Processes that touch it *and* the active modifiers that scale any of those
        Processes. A variable nothing touches is ``VALIDATED`` (it is simply
        unchanging, which is trivially correct).

        A modifier shapes a variable's value just as the Process it scales does
        (an inhibition factor changes the realised flux), so a speculative
        modifier on a validated Process drags that Process's outputs down to
        speculative — the same least-trustworthy-input rule, extended to the
        multiplicative path.

        NOTE (Milestone 1): this propagates *Process and modifier* tiers only. A
        Process/modifier must also be capped by the tiers of the parameters it
        consumes — otherwise a VALIDATED process running on speculative placeholder
        parameters would report a VALIDATED output, which is exactly the
        credibility-borrowing the tier system exists to prevent. Parameter-tier
        propagation (each declares its ``reads``) is the next task. See
        docs/DECISIONS.md D-1 and milestone-1-tasks.md.
        """
        if variable not in self.schema:
            raise KeyError(f"Unknown variable {variable!r}")
        mods = self.active_modifiers
        tiers: list[Tier] = []
        for p in self.active:
            if variable not in p.touches:
                continue
            tiers.append(p.tier)
            tiers.extend(m.tier for m in mods if p.name in m.modifies)
        return combine(tiers)

    def tier_map(self) -> dict[str, Tier]:
        """Derived tier for every state variable."""
        return {name: self.tier_of(name) for name in self.schema.names}

    def overall_tier(self) -> Tier:
        """Lowest tier across all active Processes — the tier of any output that
        mixes the whole state (e.g. a full time-series export). Includes active
        modifiers that scale at least one active Process."""
        active_names = {p.name for p in self.active}
        tiers = [p.tier for p in self.active]
        tiers.extend(m.tier for m in self.active_modifiers if active_names & set(m.modifies))
        return combine(tiers)
