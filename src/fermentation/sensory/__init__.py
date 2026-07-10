"""The sensory readout layer — speculative aroma observables over a finished trajectory.

A top-layer sibling of :mod:`fermentation.analysis` in the dependency graph: it consumes a
``runtime.Trajectory`` and a standalone threshold table and returns Odor-Activity-Values.
Nothing in ``core``/``runtime``/``scenario`` imports it (the §4.2 cardinal rule). See
:mod:`fermentation.sensory.oav` for the OAV mapping, the tier floor, and the isolation
firewall. Everything here is ``speculative`` (decision D-67).
"""

from __future__ import annotations

from fermentation.sensory.oav import (
    AROMA_COMPOUNDS,
    AromaCompound,
    OAVReading,
    SensoryProfile,
    load_thresholds,
    medium_of,
    oav_series,
    oav_tier,
    sensory_profile,
)

__all__ = [
    "AROMA_COMPOUNDS",
    "AromaCompound",
    "OAVReading",
    "SensoryProfile",
    "load_thresholds",
    "medium_of",
    "oav_series",
    "oav_tier",
    "sensory_profile",
]
