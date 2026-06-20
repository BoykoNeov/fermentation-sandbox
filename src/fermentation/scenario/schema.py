"""Declarative scenario schema (no physics lives here).

This is the I/O boundary, so fields are expressed in industry units (degrees
Brix/Plato, degrees C, days) and converted to canonical internal units when a
scenario is compiled into an initial state + Process set. That compilation step
arrives with Milestone 1; the schema is defined now so the seam is stable.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, field_validator


class TemperaturePoint(BaseModel):
    """One knot of a piecewise temperature schedule."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    day: float = Field(ge=0.0, description="Time since pitch, in days.")
    celsius: float


class Intervention(BaseModel):
    """A timed event that mutates state or changes the active Process set.

    ``action`` is a verb the runtime knows how to apply (e.g. ``add_dap``,
    ``pitch_mlf``, ``add_so2``, ``rack``); ``params`` carries its arguments in
    industry units. The set of recognised actions grows with the engine.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    day: float = Field(ge=0.0, description="Time since pitch, in days.")
    action: str
    params: dict[str, float | str] = Field(default_factory=dict)

    @field_validator("action")
    @classmethod
    def _nonempty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("action must not be empty")
        return v


class Scenario(BaseModel):
    """A fully declarative description of one fermentation run."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    #: Beverage family — selects default Process set and sugar layout downstream.
    medium: str = Field(description="e.g. 'wine', 'beer', 'cider', 'mead'")
    #: Strain identifier used to select a parameter overlay.
    strain: str = "generic"
    #: Initial composition in industry units; keys validated at compile time.
    initial: dict[str, float] = Field(default_factory=dict)
    temperature_schedule: list[TemperaturePoint] = Field(default_factory=list)
    interventions: list[Intervention] = Field(default_factory=list)
    duration_days: float = Field(default=14.0, gt=0.0)

    @field_validator("name", "medium")
    @classmethod
    def _nonempty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("must not be empty")
        return v
