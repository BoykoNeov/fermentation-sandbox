"""Fermentation Sandbox — a tiered, provenance-backed fermentation simulation engine.

The package is layered with strictly one-directional dependencies (lower layers
know nothing of higher ones):

    scenario / validation   declarative recipes, benchmark comparison, analysis
    runtime                 time-stepping, events, phase switching
    core                    pure deterministic state + Process derivatives
    parameters / units      provenance-backed data and unit conversions

See ``docs/ARCHITECTURE.md`` for the full design.
"""

__version__ = "0.0.1"
