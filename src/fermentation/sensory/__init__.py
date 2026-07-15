"""The sensory readout layer — speculative aroma observables over a finished trajectory.

A top-layer sibling of :mod:`fermentation.analysis` in the dependency graph: it consumes a
``runtime.Trajectory`` and a standalone threshold table and returns Odor-Activity-Values.
Nothing in ``core``/``runtime``/``scenario`` imports it (the §4.2 cardinal rule). See
:mod:`fermentation.sensory.oav` for the OAV mapping, the tier floor, and the isolation
firewall (beat 1a, decision D-67), and :mod:`fermentation.sensory.descriptors` for the
projection of that OAV vector onto a descriptor vocabulary (beat 1b slice 1, decision D-95).
Everything here is ``speculative``.
"""

from __future__ import annotations

from fermentation.sensory.descriptors import (
    DESCRIPTOR_AXES,
    DescriptorAxis,
    DescriptorProfile,
    DescriptorProjector,
    DescriptorReading,
    MaxRuleProjector,
    axes_for_medium,
    descriptor_tier,
)
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
    "DESCRIPTOR_AXES",
    "AromaCompound",
    "DescriptorAxis",
    "DescriptorProfile",
    "DescriptorProjector",
    "DescriptorReading",
    "MaxRuleProjector",
    "OAVReading",
    "SensoryProfile",
    "axes_for_medium",
    "descriptor_tier",
    "load_thresholds",
    "medium_of",
    "oav_series",
    "oav_tier",
    "sensory_profile",
]
