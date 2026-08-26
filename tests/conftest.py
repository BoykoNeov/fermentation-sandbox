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

from fermentation.core.chemistry import CO2_PER_HEXOSE, ETHANOL_PER_HEXOSE, M_NITROGEN
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


#: A COUNTED beer pitching rate, cells/mL — the one sourced inoculum two test modules both
#: need, kept here for the same reason :func:`seed_amino_acids` is (one statement, one means).
#: Foster 2022 (*Front. Microbiol.* 13:747546) pitched 1.2e7 cells/mL of ale yeast into a
#: 12.5 °P all-malt wort; ``test_organic_acids`` scores that trial's courses and
#: ``test_kinetics_byproducts`` uses it as the aroma calibration frame's inoculum (D-228).
#: Converted to the engine's g/L by :func:`fermentation.units.cells_per_ml_to_pitch_gpl`, the
#: same boundary D-222 put Tyrell's count through — never by a flat gram-per-litre, which
#: D-219 showed is a dry-yeast DOSING convention back-computed into a cell mass.
BEER_COUNTED_PITCH_CELLS_PER_ML = 1.2e7


#: Peyer 2017 Table 16 control column, mg/L, the 18 free amino acids of a malt wort (CW0.5+B,
#: diluted 50:50 — see :data:`PEYER_WORT_DILUTION`). Transcribed rather than trusted as a
#: literal, so both consumers re-derive from the cited source instead of from each other
#: [[feedback-transcribe-tables-not-prose]]. **Proline is absent from that table**, which is
#: exactly right and is the load-bearing property for the nitrogen half: proline is Jones &
#: Pierce Group D, brewing yeast does not assimilate it, and the ``N`` slot is assimilable
#: nitrogen by definition. So this composition IS the assimilable amino-acid set — there is no
#: FAN-to-YAN correction still owed on it, which is what D-230 checked before ruling on the
#: beer scenario's own ``yan_mgl``.
#:
#: Moved here from ``test_acidbase`` at D-230, when growth's extent audit became the SECOND
#: consumer. One transcription, one means — the reason :data:`BEER_COUNTED_PITCH_CELLS_PER_ML`
#: lives here too. The charge machinery that reads the pKa columns stays in ``test_acidbase``:
#: only the table is shared, not the use.
PEYER_WORT_AMINO_ACIDS_MGL = {
    "alanine": 36.9,
    "arginine": 47.6,
    "asparagine": 32.0,
    "aspartic": 27.5,
    "glutamic": 22.2,
    "glutamine": 41.5,
    "glycine": 11.3,
    "histidine": 22.0,
    "isoleucine": 23.4,
    "leucine": 50.7,
    "lysine": 30.2,
    "methionine": 10.2,
    "phenylalanine": 41.7,
    "serine": 23.6,
    "threonine": 20.1,
    "tryptophan": 14.0,
    "tyrosine": 30.8,
    "valine": 42.3,
}

#: ``(molar mass, nitrogen atoms, pKa_COOH, pKa_NH3, (side-chain pKa, sign) | None)``.
#: The nitrogen COUNT is the load-bearing column for both consumers: ``zbar``'s denominator is
#: ELEMENTAL nitrogen because that is what the ``N`` slot holds, so arginine's +1 spreads over
#: FOUR nitrogens and contributes +0.25 per mole N. Getting that convention backwards inflates
#: the cationic half roughly fourfold — and inflates D-230's wort-nitrogen total by the same
#: kind of factor, which is why both files assert against it rather than around it.
AMINO_ACID_CHEMISTRY = {
    "alanine": (89.09, 1, 2.34, 9.69, None),
    "arginine": (174.20, 4, 2.17, 9.04, (12.48, +1)),
    "asparagine": (132.12, 2, 2.02, 8.80, None),
    "aspartic": (133.10, 1, 1.99, 9.90, (3.90, -1)),
    "glutamic": (147.13, 1, 2.10, 9.47, (4.07, -1)),
    "glutamine": (146.15, 2, 2.17, 9.13, None),
    "glycine": (75.07, 1, 2.34, 9.60, None),
    "histidine": (155.16, 3, 1.82, 9.17, (6.04, +1)),
    "isoleucine": (131.17, 1, 2.36, 9.68, None),
    "leucine": (131.17, 1, 2.36, 9.60, None),
    "lysine": (146.19, 2, 2.18, 8.95, (10.53, +1)),
    "methionine": (149.21, 1, 2.28, 9.21, None),
    "phenylalanine": (165.19, 1, 1.83, 9.13, None),
    "serine": (105.09, 1, 2.21, 9.15, None),
    "threonine": (119.12, 1, 2.09, 9.10, None),
    "tryptophan": (204.23, 2, 2.83, 9.39, None),
    "tyrosine": (181.19, 1, 2.20, 9.11, (10.07, -1)),
    "valine": (117.15, 1, 2.32, 9.62, None),
}

#: Peyer Table 2 ("Malt wort nutrients (10-12 degP)"), row "Ammonia 25-30", read as mg N/L —
#: the reading the amino-acid units cross-check in ``test_acidbase`` establishes.
PEYER_WORT_AMMONIUM_MG_N_PER_L = (25.0, 30.0)

#: Table 16's wort is CW0.5, diluted 50:50 with water; that thesis's Table 15 confirms the 2x
#: (169 mg/L FAN undiluted against 78 diluted).
PEYER_WORT_DILUTION = 2.0


def peyer_wort_amino_acid_nitrogen_mmol_per_l() -> float:
    """Elemental nitrogen in FULL-STRENGTH Peyer wort's free amino acids, mmol/L.

    The amino-acid half of a malt wort's assimilable nitrogen. Proline is absent from the
    source table, so no non-assimilable fraction is included and none has to be subtracted.
    """
    return sum(
        mgl / AMINO_ACID_CHEMISTRY[name][0] * PEYER_WORT_DILUTION * AMINO_ACID_CHEMISTRY[name][1]
        for name, mgl in PEYER_WORT_AMINO_ACIDS_MGL.items()
    )


def peyer_wort_assimilable_nitrogen_mg_per_l() -> tuple[float, float]:
    """``(low, high)`` mg N/L of ASSIMILABLE nitrogen in a full-strength 10-12 degP malt wort.

    Amino acids (Table 16, proline-free by construction) plus ammonium (Table 2, whose printed
    25-30 range is the only band either half carries). This is the envelope D-230 scores the
    beer scenarios' assumed ``yan_mgl`` against — the assumption's first independent check in
    the ~50 records it has been carried.
    """
    amino_acid_mg = peyer_wort_amino_acid_nitrogen_mmol_per_l() * M_NITROGEN
    return tuple(amino_acid_mg + nh4 for nh4 in PEYER_WORT_AMMONIUM_MG_N_PER_L)  # type: ignore[return-value]


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
