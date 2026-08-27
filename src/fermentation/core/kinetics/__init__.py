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
    AntioxidantBurstOxidation,
    BoundHydrogenSulfideRelease,
    BoundMethanethiolRelease,
    Caramelization,
    ClosureOxygenIngress,
    EllagitanninOxidation,
    EsterHydrolysis,
    EthylAcetateEsterification,
    EthylHexanoateHydrolysis,
    MaillardBrowning,
    MaillardStrecker,
    OakExtraction,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    SMMHydrolysis,
    SotolonAldolCondensation,
    StreckerDegradation,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
    TanninEthylTanninCondensation,
    TanninSelfPolymerization,
    ThermalAnthocyaninFade,
)
from fermentation.core.kinetics.amino_acids import (
    AminoAcidAssimilation,
    AssimilableNitrogenUptake,
)
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
    EsterSynthesisGrowthCoupled,
    EsterVolatilization,
    FuselAlcoholsEhrlich,
    FuselAlcoholsEhrlichGrowthCoupled,
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
from fermentation.core.kinetics.inactivation import EthanolInactivation, EthanolToleranceDeath
from fermentation.core.kinetics.inhibition import EthanolInhibition
from fermentation.core.kinetics.keto_acids import (
    AlphaKetobutyrateExcretion,
    AlphaKetobutyrateReassimilation,
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
    malolactic_fatty_acid_gate,
    malolactic_toxicity_gate,
)
from fermentation.core.kinetics.mercaptans import AutolyticMercaptan
from fermentation.core.kinetics.organic_acids import (
    ACETIC_SLOT,
    ACETIC_SPECIES,
    ORGANIC_ACID_SPECS,
    WORT_ACID_SINKS,
    AceticAcidOverflow,
    OrganicAcidExcretion,
    OrganicAcidSpec,
    WortAcidRemoval,
    WortAcidSinkSpec,
    organic_acid_carbon_draw,
    organic_acid_rates,
)
from fermentation.core.kinetics.osmotic import OsmoticSubstrateInhibition
from fermentation.core.kinetics.oxidative_cascade import (
    OxygenActivation,
    PeroxideEthanolOxidation,
    PeroxideSulfiteOxidation,
    QuinoneAnthocyaninFading,
    QuinoneEllagitanninOxidation,
    QuinonePolymerization,
    QuinoneStreckerDegradation,
    QuinoneSulfonation,
    activation_rate,
    h2o2_branch_fraction,
)
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
from fermentation.core.kinetics.wort_oxygen import O2_SLOT, WortOxygenUptake

__all__ = [
    "AcetaldehydeBridgedCondensation",
    "AcetaldehydeProduction",
    "AcetaldehydeReduction",
    "AcetolactateDecarboxylation",
    "AcetolactateExcretion",
    "AlphaKetobutyrateExcretion",
    "AlphaKetobutyrateReassimilation",
    "AlphaKetoglutarateExcretion",
    "AlphaKetoglutarateReassimilation",
    "AminoAcidAssimilation",
    "AssimilableNitrogenUptake",
    "AnthocyaninFading",
    "AntioxidantBurstOxidation",
    "OxygenActivation",
    "PeroxideEthanolOxidation",
    "PeroxideSulfiteOxidation",
    "QuinoneAnthocyaninFading",
    "QuinoneEllagitanninOxidation",
    "QuinonePolymerization",
    "QuinoneStreckerDegradation",
    "QuinoneSulfonation",
    "activation_rate",
    "h2o2_branch_fraction",
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
    "EthylAcetateEsterification",
    "EthylHexanoateHydrolysis",
    "MaillardBrowning",
    "MaillardStrecker",
    "BoundHydrogenSulfideRelease",
    "BoundMethanethiolRelease",
    "ClosureOxygenIngress",
    "SotolonAldolCondensation",
    "ACETIC_SLOT",
    "ACETIC_SPECIES",
    "O2_SLOT",
    "ORGANIC_ACID_SPECS",
    "WORT_ACID_SINKS",
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
    "EsterSynthesisGrowthCoupled",
    "EsterVolatilization",
    "EthanolInactivation",
    "EthanolInhibition",
    "EthanolToleranceDeath",
    "FuselAlcoholsEhrlich",
    "FuselAlcoholsEhrlichGrowthCoupled",
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
    "AceticAcidOverflow",
    "OrganicAcidExcretion",
    "OrganicAcidSpec",
    "OsmoticSubstrateInhibition",
    "PyruvateExcretion",
    "PyruvateReassimilation",
    "SugarUptakeToEthanolCO2",
    "TemperatureRamp",
    "YeastAutolysis",
    "YeastPOFDecarboxylation",
    "WortAcidRemoval",
    "WortAcidSinkSpec",
    "WortOxygenUptake",
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
    "organic_acid_carbon_draw",
    "organic_acid_rates",
    "malolactic_fatty_acid_gate",
    "malolactic_toxicity_gate",
]
