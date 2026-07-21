"""Shared fixtures: a mass-conserving toy fermentation used to exercise the
runtime and the validation harness without committing to real kinetics.

Also home to :func:`seed_amino_acids`, the one place the D-100 must-spectrum seeding lives, so
that every amino-acid consumer's tests state the same thing by the same means.
"""

# Pin BLAS/OpenMP to a single thread PER PROCESS — this must run before numpy is
# first imported (below), so it lives at the very top of the root conftest. The
# suite is ~1250 independent solve_ivp integrations run process-parallel under
# ``pytest -n auto`` (pytest-xdist); with the default thread pools each of N
# workers would spawn N BLAS threads (N×N oversubscription) and the parallel run
# is *slower* than pinned — measured 382s unpinned vs 98s pinned on 16 cores.
# ``setdefault`` so an explicit outer override still wins.
import os

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

from collections.abc import Mapping

import pytest

from fermentation.core.chemistry import CO2_PER_HEXOSE, ETHANOL_PER_HEXOSE
from fermentation.core.kinetics.amino_acid_pools import AMINO_ACID_SPECS
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema, VarSpec
from fermentation.core.tiers import Tier


def seed_amino_acids(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], total: float
) -> FloatArray:
    """Load ``total`` g/L of assimilable amino acids at **must-spectrum composition** (D-100).

    The test-side twin of the compile seam's ``_wine_amino_acids`` split, and the state the
    D-100 gate algebra is designed around: at spectrum composition every per-species
    relative-depletion gate ``aa_i/(K·f_i + aa_i)`` collapses to the pre-split lumped gate
    ``aa/(K + aa)`` exactly. That is what lets the D-45/D-75/D-87/D-89 closed-form assertions
    keep asserting the *same numbers* across the split — they are not being loosened to
    accommodate speciation, they are being seeded in the state where speciation is provably a
    no-op on the rate. Anything a test still catches after this seeding is a real change.

    Mutates and returns ``y`` (the ``_wine``-builder idiom).
    """
    fractions = {spec.pool: params[spec.fraction_param] for spec in AMINO_ACID_SPECS}
    denominator = sum(fractions.values())
    for pool, fraction in fractions.items():
        y[schema.slice(pool)] = total * fraction / denominator
    return y


# Gay-Lussac mass split: glucose -> 2 ethanol + 2 CO2. Derived from the shared
# stoichiometry in fermentation.core.chemistry (single source of truth) so the
# toy's flows close to machine precision against total_carbon / total_mass.
ETHANOL_FRACTION = ETHANOL_PER_HEXOSE  # ~0.5114
CO2_FRACTION = CO2_PER_HEXOSE  # ~0.4886


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
