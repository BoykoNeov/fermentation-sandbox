"""Boundary unit conversions.

Formulas are standard enology/brewing approximations; each is cited inline. They
are convenience conversions, not kinetic parameters, so they live in code rather
than the provenance store — but the citations make their assumptions auditable.
"""

from __future__ import annotations

# Density of pure ethanol at 20 C, g/mL. (CRC Handbook.)
_ETHANOL_DENSITY_GPML = 0.78924

# Density of pure water at 20 C, g/mL.
_WATER_DENSITY_GPML = 0.99820


def celsius_to_kelvin(celsius: float) -> float:
    return celsius + 273.15


def kelvin_to_celsius(kelvin: float) -> float:
    return kelvin - 273.15


def days_to_hours(days: float) -> float:
    return days * 24.0


def hours_to_days(hours: float) -> float:
    return hours / 24.0


def mgl_to_gpl(mg_per_l: float) -> float:
    """Milligrams per litre -> grams per litre (e.g. YAN reported in mg/L).

    A plain factor of 1000, but routed through ``fermentation.units`` so every
    industry -> canonical conversion crosses the same boundary (decision D-3)
    instead of being inlined as a bare literal at the call site.
    """
    return mg_per_l / 1000.0


def gpl_to_mgl(g_per_l: float) -> float:
    """Grams per litre -> milligrams per litre (the canonical -> industry inverse).

    Free and molecular SO₂ are conventionally reported in mg/L (the ~0.5-0.8 mg/L
    molecular-SO₂ stability target, the 30-60 mg/L free-SO₂ dose), so the canonical-unit
    ``acidbase.molecular_so2`` readout (g/L) crosses back to mg/L here rather than via a
    bare ``*1000`` at the call site (decision D-3).
    """
    return g_per_l * 1000.0


def gpl_to_ugl(g_per_l: float) -> float:
    """Grams per litre -> micrograms per litre (the canonical -> odor-threshold unit).

    Literature odor-detection (perception) thresholds are reported in µg/L–mg/L, so the
    sensory OAV readout (:mod:`fermentation.sensory.oav`) compares a state concentration
    (canonical g/L) against a µg/L threshold. Both sides must share a unit for the ratio to
    be dimensionless, so the state concentration crosses to µg/L here rather than via a bare
    ``*1e6`` at the call site (decision D-3). A plain factor of one million.
    """
    return g_per_l * 1_000_000.0


def ugl_to_gpl(ug_per_l: float) -> float:
    """Micrograms per litre -> grams per litre (the odor-threshold -> canonical inverse).

    The inverse of :func:`gpl_to_ugl`, for reading a µg/L literature figure back into the
    canonical g/L state unit (decision D-3).
    """
    return ug_per_l / 1_000_000.0


def brix_to_sg(brix: float) -> float:
    """Degrees Brix -> specific gravity (20/20 C).

    Inverse of the standard cubic Brix(SG) polynomial; expressed here as the
    widely used closed form (e.g. Wikipedia "Brix", attributed to the ASBC
    extract tables):

        SG = 1 + brix / (258.6 - (brix / 258.2) * 227.1)
    """
    return 1.0 + brix / (258.6 - (brix / 258.2) * 227.1)


def sg_to_brix(sg: float) -> float:
    """Specific gravity -> degrees Brix (ASBC cubic, accurate 0-40 Brix).

    Brix = -668.962 + 1262.45*SG - 776.43*SG^2 + 182.94*SG^3
    """
    return -668.962 + 1262.45 * sg - 776.43 * sg**2 + 182.94 * sg**3


def sg_to_plato(sg: float) -> float:
    """Specific gravity -> degrees Plato (ASBC cubic).

        P = -616.868 + 1111.14*SG - 630.272*SG^2 + 135.997*SG^3

    Degrees Plato (brewing) and degrees Brix (enology) both measure % sucrose by
    mass and are numerically near-identical; the small difference comes from the
    fitting polynomial used by each industry.
    """
    return -616.868 + 1111.14 * sg - 630.272 * sg**2 + 135.997 * sg**3


def plato_to_sg(plato: float) -> float:
    """Degrees Plato -> specific gravity.

    SG = 1 + plato / (258.6 - (plato / 258.2) * 227.1)
    """
    return 1.0 + plato / (258.6 - (plato / 258.2) * 227.1)


# Balling/Tabarie apparent-vs-real extract split. As wort ferments, a hydrometer
# reads an *apparent* extract below the true (real) dissolved-solids extract,
# because the ethanol present is less dense than water. Balling's classic relation
# (1843), also attributed to Tabarie, links the two through the original extract:
#
#     RE = 0.1808 * OE + 0.8192 * AE          (degrees Plato)
#
# 0.8192 (= 1 - 0.1808) is the share of the apparent reading that is true extract;
# the 0.1808 * OE term is the ethanol-density correction, which scales with how
# much extract has fermented (OE - RE). Standard brewing-science references:
# Balling (1843); de Clerck, "A Textbook of Brewing"; ASBC Methods of Analysis.
_TABARIE_OE_SHARE = 0.1808


def real_to_apparent_extract(real_extract_plato: float, original_extract_plato: float) -> float:
    """Real (true) extract -> apparent (hydrometer) extract, both in degrees Plato.

    Inverts Balling's ``RE = 0.1808*OE + 0.8192*AE`` for ``AE``. A fermenting
    beer's hydrometer reads low because the dissolved ethanol is lighter than
    water; this is the standard correction from the true dissolved-solids extract
    to the apparent reading. Before any fermentation (``RE == OE``) it returns
    ``OE`` unchanged.
    """
    return (real_extract_plato - _TABARIE_OE_SHARE * original_extract_plato) / (
        1.0 - _TABARIE_OE_SHARE
    )


def apparent_gravity(real_extract_plato: float, original_extract_plato: float) -> float:
    """Apparent specific gravity (the hydrometer reading) of a fermenting beer.

    Composes :func:`real_to_apparent_extract` with :func:`plato_to_sg`: the true
    dissolved-solids extract is depressed to its apparent value by the ethanol
    present, then expressed as specific gravity. This is the quantity brewers mean
    by a "final gravity ~1.010" — an apparent, ethanol-depressed reading, not the
    real extract (which for a 1.048 OG ale finishes nearer 1.016).
    """
    return plato_to_sg(real_to_apparent_extract(real_extract_plato, original_extract_plato))


def brix_to_sugar_gpl(brix: float, sg: float | None = None) -> float:
    """Degrees Brix -> dissolved sugar concentration in g/L.

    Brix is grams sucrose per 100 g solution, so the volumetric concentration is

        sugar [g/L] = brix [g/100 g] * density [g/mL] * 10

    If ``sg`` is not supplied it is derived from ``brix`` via :func:`brix_to_sg`
    (density of the solution ~= SG * water density, and SG is dimensionless
    relative to water at 20/20 C so density [g/mL] ~= SG * 0.9982).
    """
    if sg is None:
        sg = brix_to_sg(brix)
    density_gpml = sg * _WATER_DENSITY_GPML
    return brix * density_gpml * 10.0


def sugar_gpl_to_brix(sugar_gpl: float, sg: float = 1.0) -> float:
    """Dissolved sugar [g/L] -> degrees Brix at a known/assumed ``sg``.

    Inverse of :func:`brix_to_sugar_gpl` for a *given* solution density. Because
    Brix and density are mutually dependent, callers that need high accuracy
    should pass the measured ``sg``; the default ``sg=1.0`` is a dilute-solution
    approximation.
    """
    density_gpml = sg * _WATER_DENSITY_GPML
    return sugar_gpl / (density_gpml * 10.0)


def abv_from_ethanol(ethanol_gpl: float) -> float:
    """Ethanol concentration [g/L] -> alcohol by volume [% v/v].

        ABV [%] = (ethanol [g/L] / ethanol_density [g/mL]) / 10

    i.e. volume of ethanol per volume of solution, in percent. Uses pure-ethanol
    density at 20 C; this ignores volume contraction on mixing and so is a close
    approximation rather than an exact figure.
    """
    ethanol_ml_per_l = ethanol_gpl / _ETHANOL_DENSITY_GPML
    return ethanol_ml_per_l / 10.0


# -- Counted pitch -> pitch_gpl (decision D-219) -------------------------------
#
# Yeast dry mass per cell, in grams. This is NOT one estimate among several: it is
# the DEFINITION of the gram this engine's biomass state variable is expressed in.
#
# Coleman, Fish & Block 2007 (Appl. Environ. Microbiol. 73(18):5875-5884, doi:
# 10.1128/aem.00670-07), Materials and Methods, "Viable and total cell concentration
# procedures", VERBATIM:
#
#     "Each cell count was converted to grams per liter of cell mass, assuming that
#      each cell weighs 4 x 10^-11 g and that the 25 squares in the hemacytometer
#      grid were covered with 10 ul."
#
# Coleman measured CELLS with a hemacytometer and never weighed any. Every gram in
# that paper -- and therefore every wine parameter this engine fits to it (X, X_A,
# Y_X/N and its nitrogen regression, k'_d, mu_max) -- is a count multiplied by this
# constant. So a counted pitch converted at any other figure enters the model in
# units its own parameters do not use. That is the whole argument for this value;
# it is an identity, not a literature preference.
#
# WHAT IT IS NOT. Coleman writes "assuming", so nothing in the chain is a weighing,
# and the tier is plausible rather than validated. Do not let a later note turn this
# into a measurement.
#
# THE INDEPENDENT CHECK, which is what says the unit is also physically honest and
# fixes the frame Coleman left open (he says "cell mass", never "dry weight"):
# invert his yield back to the quantity he actually measured, cells per gram of
# nitrogen, and price it with an elemental composition he had no hand in.
#
#     Y_X/N(330 mg N/L) = exp(3.50 - 3.61e-3 * 330) = 10.06 g cell / g N   [his Fig. 4]
#     cells per g N     = 10.06 / 4e-11            = 2.515e11
#     dry mass per g N  = 1 / 0.114                = 8.77 g               [Roels
#                         CH1.8O0.5N0.2, biomass_N_fraction, DRY basis]
#     => dry mass per cell = 8.77 / 2.515e11       = 34.9 pg
#
# 34.9 against 40 is a 13 % agreement between an assumption made for a counting
# chamber and the elemental formula, and it settles wet-vs-dry: read as a WET mass,
# 40 pg would make yeast 33 % nitrogen on a dry basis, which is absurd against the
# sourced 6.4-8.3 % (decision D-270: wine yeast at 40-45 % crude protein on the
# source's own N x 6.25, and 50 % protein at one-sixth nitrogen -- the two statements
# D-267 sec 6 found; this line carried an unsourced wider range until then, and the
# sourced ceiling makes the wet reading MORE absurd, not less).
# NOTE the fraction below must stay the STATIC elemental one: wine's compiled
# biomass_N_fraction is 1/Y_X/N, Coleman's own inverse, and feeding it here returns
# 4e-11 g identically -- the check is only independent while its composition is
# (tests/test_biomass_nitrogen_frame.py). The engine's biomass_C_fraction / biomass_N_fraction are
# declared on dry cell weight in both medium files, so this is the frame they need.
# A geometric cross-check agrees on magnitude: a 100-150 fL wine/ale cell at
# rho ~ 1.11 g/mL and 34 % dry matter is 38-57 pg. BOTH constants are sourced at
# D-271 -- Klis, de Koster & Brul 2014 (Eukaryot. Cell 13(1):2-9) states the route
# and both values; this line previously carried the dry fraction as an unsourced
# "~30 %", which its own printed 30-57 pg edges show was doing duty as a 27-34 %
# range. The sourced 0.34 is the value the UPPER edge was already computed at, so
# the repair removes a spread rather than moving the check. It is still only a
# cross-check: the shipped band below comes from the elemental route.
# D-271 also hunted the one thing that would settle this properly -- a source
# pairing a cell COUNT with a gravimetric dry WEIGHT in one ferment at one
# timepoint, D-219's own open item -- and it does not exist as published. Four
# candidates, each reporting one currency and never both, are in
# tests/test_cell_mass_literature.py. Do not re-run that search casually.
#
# THE TWO READINGS THIS SUPERSEDES, and why each was what it was:
#   * 18 pg -- asserted, unsourced, in test_validation_varela2004.py and
#     test_validation_palma2012.py as "the standard ~18 pg/cell S. cerevisiae
#     dry-weight figure". It implies a ~50 fL cell: a small lab haploid, not the
#     diploid wine strain (EC1118) those very benchmarks run.
#   * ~100 pg -- never chosen by anyone. BACK-COMPUTED from the beer scenario's
#     pitch_gpl = 1.0 against Tyrell 2013's counted 9.96e6 cells/mL, so it is a
#     residual absorbing the true cell mass AND every error in the model's
#     per-gram uptake rate. It implies a ~300 fL cell, which no S. cerevisiae is.
#
# Band: 28-50 pg, the span the elemental route gives across biomass_N_fraction's own
# 0.08-0.14 uncertainty. Both superseded readings sit OUTSIDE it.
_YEAST_DRY_MASS_PER_CELL_G = 4.0e-11


def cells_per_ml_to_pitch_gpl(cells_per_ml: float) -> float:
    """A counted pitch [cells/mL] -> this engine's ``pitch_gpl`` [g/L dry biomass].

        pitch [g/L] = cells/mL * 1e3 mL/L * 4e-11 g/cell

    The literature states pitches as counts and this engine's state variable ``X``
    is a mass, so every published trial crosses here. See the block above for why
    the constant is Coleman's 4e-11 g and not a textbook cell weight: it is the
    definition of the gram the wine parameters were fitted in, so converting at
    anything else feeds the model a number in the wrong unit (decision D-219).
    """
    return cells_per_ml * 1e3 * _YEAST_DRY_MASS_PER_CELL_G
