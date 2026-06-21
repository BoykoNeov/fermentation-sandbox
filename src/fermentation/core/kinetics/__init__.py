"""Validated-core kinetic Processes (Milestone 1).

Each primary-fermentation mechanism is a :class:`~fermentation.core.process.Process`
in its own module — growth, sugar uptake, ethanol inhibition, temperature
dependence. They are composed into a medium's ``ProcessSet`` once the full set
exists (until then they stay out of the ``MEDIA`` registry so the no-kinetics
baseline holds — see ``docs/plans/milestone-1-tasks.md``).
"""

from fermentation.core.kinetics.growth import GrowthNitrogenLimited

__all__ = ["GrowthNitrogenLimited"]
