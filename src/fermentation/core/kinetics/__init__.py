"""Validated-core kinetic mechanisms (Milestone 1).

Each primary-fermentation mechanism lives in its own module — growth, sugar
uptake, ethanol inhibition, temperature dependence. Most are additive
:class:`~fermentation.core.process.Process` objects; those that *scale* a rate
rather than add a flux (ethanol inhibition and Arrhenius temperature dependence)
are :class:`~fermentation.core.process.RateModifier` objects instead (see
decisions D-10, D-11). They are composed into a medium's ``ProcessSet`` once the
full set exists (until then they stay out of the ``MEDIA`` registry so the
no-kinetics baseline holds — see ``docs/plans/milestone-1-tasks.md``).
"""

from fermentation.core.kinetics.arrhenius import ArrheniusTemperature
from fermentation.core.kinetics.growth import GrowthNitrogenLimited
from fermentation.core.kinetics.inhibition import EthanolInhibition
from fermentation.core.kinetics.uptake import SugarUptakeToEthanolCO2

__all__ = [
    "ArrheniusTemperature",
    "EthanolInhibition",
    "GrowthNitrogenLimited",
    "SugarUptakeToEthanolCO2",
]
