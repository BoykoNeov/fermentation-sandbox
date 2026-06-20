"""Pydantic models that enforce parameter provenance at load time."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from fermentation.core.tiers import Tier


class Uncertainty(BaseModel):
    """A plausible range for a parameter value.

    Required for every parameter — even a pure guess must state how unsure it is.
    Interpreted as a (rough) interval that should bracket the true value;
    ``value`` must lie within ``[low, high]``.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    low: float
    high: float
    note: str = ""

    @model_validator(mode="after")
    def _ordered(self) -> Uncertainty:
        if self.low > self.high:
            raise ValueError(f"uncertainty low ({self.low}) > high ({self.high})")
        return self


class Provenance(BaseModel):
    """Where a parameter came from and under what conditions it was measured."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    #: Citation, or the literal string "author estimate" for an honest guess.
    source: str
    #: Conditions the value was measured/estimated under (strain, T, medium, ...).
    conditions: str
    doi: str | None = None
    notes: str = ""

    @field_validator("source", "conditions")
    @classmethod
    def _nonempty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("must not be empty")
        return v


class Parameter(BaseModel):
    """A single kinetic/physical constant with full provenance.

    Constructing one with a missing required field raises ``ValidationError`` —
    this is the mechanism that makes "no magic numbers, ever" a hard guarantee
    rather than a guideline.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    value: float
    unit: str
    tier: Tier
    uncertainty: Uncertainty
    provenance: Provenance

    @field_validator("tier", mode="before")
    @classmethod
    def _coerce_tier(cls, v: Any) -> Tier:
        """Accept the human label ('validated'/'plausible'/'speculative') as well
        as the raw int, so YAML files can be written readably."""
        if isinstance(v, Tier):
            return v
        if isinstance(v, str):
            try:
                return Tier[v.strip().upper()]
            except KeyError:
                raise ValueError(
                    f"unknown tier {v!r}; expected one of {[t.label for t in Tier]}"
                ) from None
        return Tier(v)

    @field_validator("unit")
    @classmethod
    def _unit_nonempty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("unit must not be empty (use 'dimensionless' if so)")
        return v

    @model_validator(mode="after")
    def _value_in_range(self) -> Parameter:
        u = self.uncertainty
        if not (u.low <= self.value <= u.high):
            raise ValueError(
                f"parameter {self.name!r}: value {self.value} outside "
                f"uncertainty range [{u.low}, {u.high}]"
            )
        return self
