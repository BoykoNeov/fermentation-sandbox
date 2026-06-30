"""Tests for the shared stoichiometric constants (fermentation.core.chemistry)."""

import pytest

from fermentation.core import chemistry as chem


def test_molar_masses_match_formulae():
    # Standard atomic masses C 12.011, H 1.008, O 15.999.
    masses = {
        chem.M_GLUCOSE: 180.156,
        chem.M_MALTOSE: 342.297,
        chem.M_MALTOTRIOSE: 504.438,
        chem.M_ETHANOL: 46.069,
        chem.M_CO2: 44.009,
        # Aroma-byproduct representative species (decision D-19).
        chem.M_ETHYL_ACETATE: 88.106,  # C4H8O2
        chem.M_ISOAMYL_OH: 88.150,  # C5H12O
    }
    for actual, expected in masses.items():
        assert actual == pytest.approx(expected, abs=1e-3)


def test_carbon_mass_fraction_known_species():
    # Glucose is 6 carbons in 180.156 g/mol -> 40% carbon by mass.
    assert chem.carbon_mass_fraction("glucose") == pytest.approx(0.4000, abs=1e-3)
    assert chem.carbon_mass_fraction("ethanol") == pytest.approx(0.5214, abs=1e-3)
    assert chem.carbon_mass_fraction("CO2") == pytest.approx(0.2729, abs=1e-3)
    # Maltose/maltotriose are carbon-richer per gram (less water of hydration).
    assert chem.carbon_mass_fraction("maltose") == pytest.approx(0.4211, abs=1e-3)
    assert chem.carbon_mass_fraction("maltotriose") == pytest.approx(0.4286, abs=1e-3)
    # Aroma-byproduct species the ester/fusel pools book against (decision D-19);
    # both are carbon-rich (reduced organics) — ethyl acetate 4 C / 88.11 g/mol,
    # isoamyl alcohol 5 C / 88.15 g/mol.
    assert chem.carbon_mass_fraction("ethyl_acetate") == pytest.approx(0.5453, abs=1e-3)
    assert chem.carbon_mass_fraction("isoamyl_alcohol") == pytest.approx(0.6813, abs=1e-3)


def test_carbon_mass_fraction_unknown_species_raises():
    with pytest.raises(KeyError, match="unknown species"):
        chem.carbon_mass_fraction("sucrose")


def test_hexose_fermentation_is_mass_balanced():
    # C6H12O6 -> 2 C2H5OH + 2 CO2: product mass equals substrate mass.
    product_mass = 2 * chem.M_ETHANOL + 2 * chem.M_CO2
    assert product_mass == pytest.approx(chem.M_GLUCOSE, abs=1e-9)
    # ...so the per-gram ethanol + CO2 yields sum to 1.
    yield_sum = chem.ETHANOL_PER_HEXOSE + chem.CO2_PER_HEXOSE
    assert yield_sum == pytest.approx(1.0, abs=1e-12)
    # The realised theoretical ethanol yield is the well-known ~0.511 g/g.
    ethanol_yield = chem.ETHANOL_PER_HEXOSE
    assert ethanol_yield == pytest.approx(0.511, abs=2e-3)


def test_hexose_fermentation_is_carbon_balanced():
    # Carbon leaving the sugar equals carbon entering ethanol + CO2 (per gram hexose).
    carbon_in = chem.carbon_mass_fraction("glucose")
    carbon_out = (
        chem.carbon_mass_fraction("ethanol") * chem.ETHANOL_PER_HEXOSE
        + chem.carbon_mass_fraction("CO2") * chem.CO2_PER_HEXOSE
    )
    assert carbon_out == pytest.approx(carbon_in, abs=1e-12)
