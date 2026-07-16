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

from fermentation.core.kinetics.acetaldehyde import (
    AcetaldehydeProduction,
    AcetaldehydeReduction,
)
from fermentation.core.kinetics.aging import (
    AcetaldehydeBridgedCondensation,
    AnthocyaninFading,
    Caramelization,
    EllagitanninOxidation,
    EsterHydrolysis,
    MaillardBrowning,
    MaillardStrecker,
    OakExtraction,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    SMMHydrolysis,
    StreckerDegradation,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
    TanninEthylTanninCondensation,
    TanninSelfPolymerization,
    ThermalAnthocyaninFade,
)
from fermentation.core.kinetics.amino_acids import AminoAcidAssimilation
from fermentation.core.kinetics.arrhenius import (
    ArrheniusTemperature,
    ColemanQuadraticDeathTemperature,
    arrhenius_factor,
)
from fermentation.core.kinetics.autolysis import YeastAutolysis, autolysis_flux
from fermentation.core.kinetics.brett import (
    BrettDeath,
    BrettDecarboxylation,
    BrettEthanolToxicity,
    BrettGrowth,
    BrettVinylphenolReduction,
    YeastPOFDecarboxylation,
    brett_environmental_gate,
    brett_ethanol_survival_factor,
)
from fermentation.core.kinetics.byproducts import (
    EsterSynthesis,
    EsterVolatilization,
    FuselAlcoholsEhrlich,
    FuselAminoAcidReroute,
    fusel_carbon_draw,
    fusel_production_rate,
    fusel_rate_shape,
)
from fermentation.core.kinetics.carrying_capacity import BiomassCarryingCapacity
from fermentation.core.kinetics.growth import GrowthNitrogenLimited, biomass_growth_rate
from fermentation.core.kinetics.hops import (
    IsoAlphaAcidLoss,
    boil_rate_constants,
    iso_alpha_fraction,
)
from fermentation.core.kinetics.hydrogen_sulfide import (
    AutolyticHydrogenSulfide,
    HydrogenSulfideProduction,
    HydrogenSulfideVolatilization,
)
from fermentation.core.kinetics.inactivation import EthanolInactivation
from fermentation.core.kinetics.inhibition import EthanolInhibition
from fermentation.core.kinetics.keto_acids import (
    AlphaKetoglutarateExcretion,
    AlphaKetoglutarateReassimilation,
    PyruvateExcretion,
    PyruvateReassimilation,
)
from fermentation.core.kinetics.malolactic import (
    MalolacticCitrateMetabolism,
    MalolacticConversion,
    MalolacticDeath,
    MalolacticGrowth,
    MalolacticSenescence,
    OenococcusDiacetylReduction,
    cardinal_temperature_factor,
    malolactic_environmental_gate,
    malolactic_toxicity_gate,
)
from fermentation.core.kinetics.mercaptans import AutolyticMercaptan
from fermentation.core.kinetics.precursor_fates import (
    NON_EHRLICH_FRACTION_PARAMS,
    PrecursorNonEhrlichFates,
    non_ehrlich_fraction_param,
)
from fermentation.core.kinetics.temperature import TemperatureRamp
from fermentation.core.kinetics.uptake import SugarUptakeToEthanolCO2
from fermentation.core.kinetics.vicinal_diketones import (
    AcetolactateDecarboxylation,
    AcetolactateExcretion,
    DiacetylReduction,
)

__all__ = [
    "AcetaldehydeBridgedCondensation",
    "AcetaldehydeProduction",
    "AcetaldehydeReduction",
    "AcetolactateDecarboxylation",
    "AcetolactateExcretion",
    "AlphaKetoglutarateExcretion",
    "AlphaKetoglutarateReassimilation",
    "AminoAcidAssimilation",
    "AnthocyaninFading",
    "AutolyticHydrogenSulfide",
    "AutolyticMercaptan",
    "BiomassCarryingCapacity",
    "BrettDeath",
    "BrettDecarboxylation",
    "BrettEthanolToxicity",
    "BrettGrowth",
    "BrettVinylphenolReduction",
    "DiacetylReduction",
    "ArrheniusTemperature",
    "ColemanQuadraticDeathTemperature",
    "Caramelization",
    "EllagitanninOxidation",
    "EsterHydrolysis",
    "MaillardBrowning",
    "MaillardStrecker",
    "OakExtraction",
    "OxidativeAcetaldehyde",
    "PhenolicBrowning",
    "SMMHydrolysis",
    "StreckerDegradation",
    "SulfiteOxidation",
    "TanninAnthocyaninCondensation",
    "TanninEthylTanninCondensation",
    "TanninSelfPolymerization",
    "ThermalAnthocyaninFade",
    "EsterSynthesis",
    "EsterVolatilization",
    "EthanolInactivation",
    "EthanolInhibition",
    "FuselAlcoholsEhrlich",
    "FuselAminoAcidReroute",
    "NON_EHRLICH_FRACTION_PARAMS",
    "PrecursorNonEhrlichFates",
    "non_ehrlich_fraction_param",
    "GrowthNitrogenLimited",
    "HydrogenSulfideProduction",
    "HydrogenSulfideVolatilization",
    "IsoAlphaAcidLoss",
    "MalolacticCitrateMetabolism",
    "MalolacticConversion",
    "MalolacticDeath",
    "MalolacticGrowth",
    "MalolacticSenescence",
    "OenococcusDiacetylReduction",
    "PyruvateExcretion",
    "PyruvateReassimilation",
    "SugarUptakeToEthanolCO2",
    "TemperatureRamp",
    "YeastAutolysis",
    "YeastPOFDecarboxylation",
    "arrhenius_factor",
    "autolysis_flux",
    "biomass_growth_rate",
    "boil_rate_constants",
    "brett_environmental_gate",
    "brett_ethanol_survival_factor",
    "cardinal_temperature_factor",
    "fusel_carbon_draw",
    "fusel_production_rate",
    "fusel_rate_shape",
    "iso_alpha_fraction",
    "malolactic_environmental_gate",
    "malolactic_toxicity_gate",
]
