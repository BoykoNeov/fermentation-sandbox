"""Scenarios as data — a run described declaratively, with no physics.

A scenario captures initial composition, organism(s)/strain, temperature
schedule, vessel, and a timeline of interventions. Keeping physics out of this
layer makes parameter sweeps, Monte Carlo, and scenario sharing trivial, and
keeps the engine reusable across wine, beer, cider, and mead without code
changes. Scenarios are schema-validated YAML/JSON — deliberately not a custom DSL.
"""

from fermentation.scenario.compile import CompiledScenario, compile_scenario
from fermentation.scenario.schema import Intervention, Scenario, TemperaturePoint

__all__ = [
    "CompiledScenario",
    "Intervention",
    "Scenario",
    "TemperaturePoint",
    "compile_scenario",
]
