"""Shared fixtures: a mass-conserving toy fermentation used to exercise the
runtime and the validation harness without committing to real kinetics.
"""

from collections.abc import Mapping

import pytest

from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema, VarSpec
from fermentation.core.tiers import Tier

# Gay-Lussac mass split: glucose (180) -> 2 ethanol (92) + 2 CO2 (88).
ETHANOL_FRACTION = 92.0 / 180.0  # 0.5111...
CO2_FRACTION = 88.0 / 180.0  # 0.4888...


class MassConservingFermentation(Process):
    """Saturating sugar uptake split into ethanol + CO2 by mass.

    No biomass growth, so total mass S + E + CO2 is conserved exactly — ideal for
    testing the conservation harness. Not real kinetics; just a clean invariant.
    """

    name = "toy_mass_conserving"
    tier = Tier.VALIDATED
    touches = ("S", "E", "CO2")

    def __init__(self, vmax: float = 5.0, ks: float = 5.0):
        self.vmax = vmax
        self.ks = ks

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        s = schema.get(y, "S")
        if s <= 0:
            return d
        consume = self.vmax * s / (self.ks + s)
        d[schema.slice("S")] = -consume
        d[schema.slice("E")] = consume * ETHANOL_FRACTION
        d[schema.slice("CO2")] = consume * CO2_FRACTION
        return d


@pytest.fixture
def toy_schema() -> StateSchema:
    return StateSchema(
        [
            VarSpec("S", "g/L", description="sugar"),
            VarSpec("E", "g/L", description="ethanol"),
            VarSpec("CO2", "g/L", description="evolved CO2"),
        ]
    )


@pytest.fixture
def toy_process() -> MassConservingFermentation:
    return MassConservingFermentation()
