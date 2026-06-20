"""Domain core: pure, deterministic, no I/O, no global state, no randomness.

Given a state and a parameter set, the core returns derivatives. That purity is
what makes it testable against benchmark curves and conservation laws.
"""

from fermentation.core.process import Process, ProcessSet
from fermentation.core.state import StateSchema, StateVector
from fermentation.core.tiers import Tier

__all__ = ["Process", "ProcessSet", "StateSchema", "StateVector", "Tier"]
