"""Beer's organic-acid production — the acids yeast makes while it ferments (decision D-180).

Until D-180 beer's acid slots (D-179) were a **composition**: dosed at pitch, touched by no
Process, so a modelled beer's pH sat exactly where ``initial_ph`` put it for the whole run
while a real beer's falls by most of a pH unit. This module gives four of those five slots a
producer, which turns beer's pH from an input into a **prediction**.

The other half of that change lives at the compile seam, and the two are inseparable: the
seeds moved from finished-beer levels to **wort** levels (``*_typical_wort``,
``beer_acids.yaml``). A producer bolted onto finished-beer seeds would finish at roughly twice
the measured beer — verified at pH 4.26, below any real beer — so "produce the acids" and
"start from a wort" are one decision, not two.

**The source, and why it is unusually good for this engine.** Tyrell et al. (2013),
*BrewingScience* 66:75-84, fermented ONE 12 °Plato wort with FOUR strains in EBC tubes and
published, for the same ferments, both every organic acid's day-0-to-day-7 course (their
Figs 6-14) *and* the pH and extract curves (their Fig 4). So the yields are a measured
**difference on one wort**, not two studies subtracted — and the pH curve is then a **free
prediction**, because nothing in the parameter file is fitted to it.

**What the free prediction says (read the D-180 record before quoting it).** Against a
measured drop of 0.75-0.87 pH units the model produces 0.514-0.744 across the sampled
``pKa_peptide_buffer`` band. That is most of it, and it is **not** a validation, because the
agreement is held open by two omitted terms of opposite sign:

* three wort acids that FALL over a real ferment (pyruvic 22→~1, formic 26→~5, oxalic
  22→~5 ppm) are not beer state slots, so the model cannot lose that anion charge — worth
  about +0.2-0.3 pH of missing base;
* dissolved CO₂ at end-of-fermentation saturation is not in the charge balance — worth about
  −0.3 pH. (The engine's ``CO2`` pool is *cumulative evolved gas*, not dissolved
  concentration, so it must NOT be read as carbonic acid; that needs a saturation model.)

Each is a follow-up beat with a measured size. Neither is a reason to tune a yield.

**Mechanism — an extra sliver off ``S``, never uptake's yields.** There are two shipped ways
to fund a byproduct's carbon, and the choice matters:

* D-16's ``Y_glycerol_sugar``/``Y_byproduct_sugar`` rescale the ethanol/CO₂ split *inside*
  uptake. Beer holds both at 0 so its sugar→ethanol stays theoretical and its CO₂-ratio
  benchmark is byte-for-byte validated. Routing acids that way would move that benchmark.
* D-19's byproduct idiom pulls an *additional* sliver of ``S`` into the producing pool and
  leaves ``dX``/``dN``/``dE``/``dCO2`` untouched. All four acids together divert
  2.435e-3 g/g — about a fifth of wine's single ``Y_byproduct_sugar`` and ~0.24 % of the
  sugar — so the drift is of the same order as the fusel/ester draws already accepted.

This module takes the second. ``Byp`` therefore stays empty for beer and beer's
``Y_byproduct_sugar`` stays 0 — which is now **load-bearing rather than incidental**: beer
carries its own ``succinic`` slot *and* the charge balance reads ``Byp`` as
succinic-equivalent (:data:`~fermentation.core.acidbase.BYP_AS_SUCCINIC`), so a non-zero beer
``Y_byproduct_sugar`` would count succinic twice in the pH solve. Pinned by
``tests/test_organic_acids.py::test_beer_byproduct_yield_stays_zero_or_succinic_double_counts``.

**Flux-linked, and the shapes that costs.** Production tracks the fermentative sugar-uptake
rate, so the run integral is a clean yield (``Y · ΔS``) on any beer, and it stops at dryness.
Against Tyrell's own curves that is right for succinic and wrong in two named ways:

* **acetic** really peaks at ~200 ppm around day 2 (3.4× the wort level) and falls back to
  105-126 by day 7 as the yeast re-assimilate it. A monotone producer hits the endpoint and
  drops the peak — the same concession :mod:`~fermentation.core.kinetics.keto_acids` already
  makes for overflow pyruvate, and it needs the same excretion/re-assimilation *pair* to fix.
* **lactic** rises hardest *after* the extract curve goes flat (days 2-7). A flux-linked form
  makes it early instead. Tyrell name the candidate mechanism — end-of-ferment autolysis —
  and this engine has an autolysis Process, so a later beat has a route.

Reported, not tuned: one dataset cannot separate a late excretion from an autolytic release.

**Tier: plausible.** The yields are measured differences on a real wort with named strains at
both band edges, but the divisor passes through two derived steps (real-vs-apparent
attenuation, and a Plato→SG linearisation) and every level is a *figure read* rather than a
printed table. Parameter-tier propagation (D-1) does the rest.

**Isolable (prime directive 3).** Beer-only — no wine Process reads these yields, so a wine
ensemble never draws them — and disabled at the compile seam whenever ``initial_ph`` is
absent, exactly like the acid slots it fills. An un-anchored beer is therefore byte-for-byte
the pre-D-179 beer, and cannot acquire a charge balance mid-run by producing into empty slots.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from fermentation.core.chemistry import (
    M_ACETIC,
    M_LACTIC,
    M_MALIC,
    M_SUCCINIC,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.carbon_routing import (
    draw_carbon_from_sugar as _draw_carbon_from_sugar,
)
from fermentation.core.kinetics.carbon_routing import (
    fermentative_uptake_rates as _fermentative_uptake_rates,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


@dataclass(frozen=True)
class OrganicAcidSpec:
    """One produced organic acid: its state slot, its yield parameter, its carbon species.

    ``species`` is the key :func:`~fermentation.core.chemistry.carbon_mass_fraction` weights
    the sugar draw by, so each acid's carbon is booked at **its own** molecule — the D-99
    correction applied here from the start rather than after a lumped stand-in ships.
    """

    slot: str
    yield_param: str
    species: str
    molar_mass: float


#: The acids beer's yeast PRODUCES, and the single source of truth every layer derives from
#: (the ``ESTER_SPECS``/``FUSEL_SPECS`` registry discipline, D-96/D-99): the Process's
#: ``touches`` and ``reads``, the carbon draw, and the tests all read this tuple, so adding a
#: fifth produced acid is one entry rather than a new code path.
#:
#: **``citrate`` is deliberately absent and must stay absent.** Tyrell's section 1.5 ("citric
#: … basically depends on concentrations of wort"), their Table 2 ("final concentration mostly
#: determined by wort concentration") and their Fig. 12 (four strains scattering ±20 ppm around
#: the wort value with no trend) all say the same thing three ways: beer's citrate is
#: malt-derived, not a fermentation product. It keeps a wort seed and stays inert.
#:
#: The three acids that *fall* over a real ferment — pyruvic, formic, oxalic — are absent for a
#: different reason: none is a beer state slot (D-179 refused pyruvic on the name collision
#: with wine's dynamic ``pyruvate`` pool), so there is nothing to drain. See the module
#: docstring for what that omission is worth in pH.
ORGANIC_ACID_SPECS: tuple[OrganicAcidSpec, ...] = (
    OrganicAcidSpec("acetic", "Y_acetic_sugar_beer", "acetic_acid", M_ACETIC),
    OrganicAcidSpec("lactic", "Y_lactic_sugar_beer", "lactic_acid", M_LACTIC),
    OrganicAcidSpec("succinic", "Y_succinic_sugar_beer", "succinic_acid", M_SUCCINIC),
    OrganicAcidSpec("malic", "Y_malic_sugar_beer", "malic_acid", M_MALIC),
)


def organic_acid_rates(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float]
) -> list[tuple[OrganicAcidSpec, float]]:
    """Each produced acid's rate ``d(<slot>)/dt`` [g/L/h], paired with its spec.

    ``Y_<acid> · Σ_i r_i`` over the per-slot fermentative uptake rates — so the run integral
    is ``Y · ΔS``, a yield in the ordinary sense, on a wine-shaped single sugar slot or beer's
    three alike. Returns ``[]`` when nothing is fermenting, which is what stops production at
    dryness without a separate gate.
    """
    flux = sum(_fermentative_uptake_rates(y, schema, params))
    if flux <= 0.0:
        return []
    return [(spec, params[spec.yield_param] * flux) for spec in ORGANIC_ACID_SPECS]


def organic_acid_carbon_draw(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float]
) -> float:
    """Total carbon [g C/L/h] the producer books out of sugar, across all four acids.

    Each acid contributes ``rate · carbon_mass_fraction(species)`` at its **own** molecule's
    fraction. Factored out of the Process for the same reason the fusel draw is (D-33/D-99):
    a future consumer that refunds part of this carbon to another source must be able to
    compute the identical number rather than re-derive it.
    """
    return float(
        sum(
            rate * carbon_mass_fraction(spec.species)
            for spec, rate in organic_acid_rates(y, schema, params)
        )
    )


class OrganicAcidExcretion(Process):
    """Yeast-produced organic acids, flux-linked, carbon-routed out of sugar (decision D-180).

    ``d(<acid>)/dt = Y_<acid> · Σ_i r_i`` for each entry of :data:`ORGANIC_ACID_SPECS`, with
    the acids' carbon drawn out of ``S`` via
    :func:`~fermentation.core.kinetics.carbon_routing.draw_carbon_from_sugar` (option a1,
    D-19). So this touches the four acid slots and ``S`` — never ``E``/``CO2`` — and
    ``total_carbon`` closes exactly, because the carbon removed from sugar equals the carbon
    deposited in the pools by construction.

    The rates come from the shared
    :func:`~fermentation.core.kinetics.carbon_routing.fermentative_uptake_rates`, which the
    uptake Process also consumes. That sharing is what makes ``Y`` a yield rather than a
    coincidence, and it carries an obligation: **every RateModifier that scales uptake must
    also name this Process**, or the acid would be booked against a flux the solver never ran.
    Beer's medium definition does exactly that (``ArrheniusTemperature.for_uptake`` with this
    Process as an extra target) — the D-32 correctness coupling, in the form it takes here.

    **No CO₂ co-product**, matching the fusel sugar stand-in: routing sugar carbon into an
    acid pool asserts no named reaction, so no decarboxylation is claimed. Where this engine
    *does* assert a reaction (the D-106 Ehrlich branch) it charges the CO₂; here it does not.
    Real succinate formation via the reductive TCA branch and the GABA shunt draws on
    glutamate as much as on hexose, so the sugar source is a bookkeeping stand-in — the
    ledger closes, the metabolic claim is not made.
    """

    name = "organic_acid_excretion"
    tier = Tier.PLAUSIBLE
    #: The four produced acid slots plus ``S``, derived from the registry so a fifth acid
    #: cannot silently violate the ``touches`` contract.
    touches: tuple[str, ...] = (*(spec.slot for spec in ORGANIC_ACID_SPECS), "S")
    #: The four yields, plus the three uptake constants the shared rate helper reads —
    #: declared here too, because ``reads`` has two masters (tier propagation AND sampler
    #: scope, D-160): an undeclared read does not merely under-document a dependency, it
    #: narrows the reported ensemble spread below what the parameter's own provenance
    #: justifies.
    reads: tuple[str, ...] = (
        *(spec.yield_param for spec in ORGANIC_ACID_SPECS),
        "q_sugar_max",
        "K_sugar_uptake",
        "K_repression",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        rates = organic_acid_rates(y, schema, params)
        if not rates:
            return d
        for spec, rate in rates:
            d[schema.slice(spec.slot)] = rate
        # ONE draw for all four, each at its own molecule's carbon fraction.
        _draw_carbon_from_sugar(d, y, schema, organic_acid_carbon_draw(y, schema, params))
        return d
