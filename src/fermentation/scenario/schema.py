"""Declarative scenario schema (no physics lives here).

This is the I/O boundary, so fields are expressed in industry units (degrees
Brix/Plato, degrees C, days) and converted to canonical internal units when a
scenario is compiled into an initial state + Process set. That compilation step
arrives with Milestone 1; the schema is defined now so the seam is stable.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class TemperaturePoint(BaseModel):
    """One knot of a piecewise temperature schedule."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    day: float = Field(ge=0.0, description="Time since pitch, in days.")
    celsius: float


class HopAddition(BaseModel):
    """One hop addition to the boil (decision D-64).

    Bitterness (iso-alpha-acids) is created in the boil by thermal isomerization of the hop's
    alpha-acids, at a rate that depends on the contact time. Each addition names its
    alpha-acid content, mass, and boil contact time; the compile seam runs the Malowicki
    closed-form isomerization for each and sums the resulting iso-alpha into the initial
    ``iso_alpha`` state. Beer-only. Dry-hop / whirlpool (post-boil) additions are a documented
    v1 deferral — every addition here is treated as a kettle (boiling) addition.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    alpha_acid_percent: float = Field(
        gt=0.0, le=100.0, description="Alpha-acid content of the hop, % w/w (e.g. 5.0 for 5% AA)."
    )
    grams: float = Field(gt=0.0, description="Mass of hop added, grams.")
    boil_minutes: float = Field(
        ge=0.0, description="Contact time in the boil, minutes (0 ⇒ no isomerization)."
    )


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
    #: Hop schedule for bitterness (beer, decision D-64). Each addition isomerizes in the boil
    #: into the initial ``iso_alpha`` state at the compile seam. Empty ⇒ an unhopped beer (the
    #: bitterness beat is inert and the loss Process disabled). Requires ``batch_volume_liters``.
    hops: list[HopAddition] = Field(default_factory=list)
    #: Wort volume [L] — the batch size hop *mass* is diluted into (grams → g/L). Genuinely new
    #: to the engine, which is otherwise volume-agnostic (concentration-based); only hop dosing
    #: needs an absolute mass→concentration conversion. v1 uses ONE volume for boil and fermenter
    #: (kettle-loss/evaporation folded into ``hop_utilization_efficiency``). Required iff ``hops``.
    batch_volume_liters: float | None = Field(default=None, gt=0.0)
    #: Boil temperature [°C]. Default 100 (sea-level boil); lower for a whirlpool/high-altitude
    #: boil, which slows isomerization (the Malowicki Arrhenius temperature dependence).
    boil_celsius: float = Field(default=100.0, gt=0.0)
    #: Bottle closure, naming its steady oxygen transmission rate (wine, decision D-136). One of
    #: ``hermetic`` / ``technical_cork`` / ``screwcap`` / ``natural_cork`` / ``synthetic_nomacorc``
    #: / ``synthetic_supremecorq``; the compile seam looks the name up in ``closure.yaml`` and seeds
    #: the ``closure_otr`` state slot from the sourced Parameter, so the choice carries provenance
    #: instead of being a bare number. A NAMED MENU rather than a float precisely because prime
    #: directive #2 forbids inlining a constant — an arbitrary OTR would have no source.
    #:
    #: ``None`` (the default) is treated as ``hermetic``: no ingress, the whole aging axis
    #: byte-for-byte as it was before D-136. This is a genuine physical limiting case, not just a
    #: disable switch — Lopes et al. 2007 found a flame-sealed bottle to be the only truly air-tight
    #: seal they tested. Naming ``hermetic`` explicitly is the way to say "sealed, and I mean it".
    #:
    #: STEADY permeation only. The bottling burst (trapped cork/headspace air, 10–150× the steady
    #: rate over the first month) is a separate one-off dose — add an ``add_oxygen`` intervention at
    #: the ``begin_aging`` day for a freshly bottled wine. Scope: a standard 750 mL bottle.
    closure: str | None = Field(default=None)
    duration_days: float = Field(default=14.0, gt=0.0)

    @field_validator("name", "medium")
    @classmethod
    def _nonempty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("must not be empty")
        return v

    @model_validator(mode="after")
    def _hops_need_volume(self) -> Scenario:
        if self.hops and self.batch_volume_liters is None:
            raise ValueError(
                "scenario has 'hops' but no 'batch_volume_liters'; hop mass (grams) needs a "
                "wort volume to become a concentration (decision D-64)"
            )
        return self
