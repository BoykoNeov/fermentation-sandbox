"""Domain core: pure, deterministic, no I/O, no global state, no randomness.

Given a state and a parameter set, the core returns derivatives. That purity is
what makes it testable against benchmark curves and conservation laws.
"""

from fermentation.core.kinetics import (
    EthanolInhibition,
    GrowthNitrogenLimited,
    SugarUptakeToEthanolCO2,
)
from fermentation.core.media import MEDIA, Medium, beer_schema, get_medium, wine_schema
from fermentation.core.process import Process, ProcessSet, RateModifier
from fermentation.core.state import StateSchema, StateVector
from fermentation.core.tiers import Tier

__all__ = [
    "MEDIA",
    "EthanolInhibition",
    "GrowthNitrogenLimited",
    "Medium",
    "Process",
    "ProcessSet",
    "RateModifier",
    "StateSchema",
    "StateVector",
    "SugarUptakeToEthanolCO2",
    "Tier",
    "beer_schema",
    "get_medium",
    "wine_schema",
]
