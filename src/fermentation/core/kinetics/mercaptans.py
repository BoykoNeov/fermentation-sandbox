"""Mercaptans (volatile thiols) — the carbon-bearing autolytic reductive off-aroma (decision D-45).

The second reductive-fault sulfur beat, after the autolytic H₂S source (D-44). Beyond H₂S, the
other "reduction" off-aromas are the **mercaptans** (thiols): methanethiol ("cooked cabbage",
sensory threshold ~2–3 µg/L) and its onion/rubber sibling ethanethiol. Lumped into one produced
pool booked as **methanethiol** (:data:`~fermentation.core.chemistry.M_METHANETHIOL`) — the honest
single-species stand-in, the arginine-for-``amino_acids`` / p-coumaric-for-``hydroxycinnamics``
idiom (D-32/D-40).

**Formation — autolysis-linked, carbon from the amino-acid pool (Option A, owner-chosen).** Real
methanethiol comes chiefly from **methionine degradation** on the lees during self-digestion, so —
like the D-44 H₂S source — its rate is a **yield on the shared autolysis flux**
(:func:`~fermentation.core.kinetics.autolysis.autolysis_flux`, ``r = k_autolysis·f_T·X_dead``). But
unlike carbon-free H₂S, **methanethiol carries carbon**, so the Process cannot draw from nothing —
it must debit a tracked carbon pool. It draws that carbon from the ``amino_acids`` pool and
**deaminates** the nitrogen back to ammonium ``N`` — the exact D-33
:class:`~fermentation.core.kinetics.byproducts.FuselAminoAcidReroute` idiom (draw carbon from amino
acids, release their nitrogen to ``N``):

    r_merc      = y_mercaptan · autolysis_flux · [aa/(K_amino_acids + aa)]      [g MeSH/L/h]
    d[mercaptans] = +r_merc
    d[amino_acids] = −(r_merc · c_merc) / c_aa      (arginine mass carrying that carbon)
    d[N]           = +(that arginine mass) · y_N     (DEAMINATION → ammonium)

* **Availability gate** ``aa/(K_amino_acids + aa)`` (the swap/re-route gate): production ramps down
  smoothly to 0 as the pool empties, so the draw can never drive ``amino_acids`` negative (a
  solver-safe C¹ shadow, no hard clamp discontinuity). The pool is kept non-empty by the D-34
  autolysis refill, which the same flux runs — so mercaptan formation trails the refill.
* **Not flux-linked to fermentation** (first-order in ``X_dead``, via ``autolysis_flux``), so — like
  the D-44 H₂S source — it fires post-dryness and there is no CO₂ stream to strip it: mercaptans
  **accumulate as residual**, the reductive fault copper fining removes.

**PROVENANCE CAVEAT — the arginine lump, not literal methionine.** ``amino_acids`` is booked as
*arginine*, so the carbon and nitrogen this Process draws are **arginine's**, not methionine's: the
model releases ~0.66 mol N per mol MeSH (arginine's C:N through methanethiol's carbon) against real
methionine's ~1. Same order of magnitude, so no gross artifact — but this is the arginine-for-
``amino_acids`` stand-in (**exact on the carbon/nitrogen ledger, approximate on provenance**), *not*
faithful methionine chemistry. The "carbon from methionine" story motivates Option A; the ledger
sees arginine.

**Conservation — closes on both ledgers by construction (no new conservation code beyond weighting
``mercaptans`` in ``total_carbon``).** Carbon: the carbon into ``mercaptans`` (``r_merc·c_merc``)
equals the carbon out of ``amino_acids`` (``aa_mass·c_aa``) — the draw is sized to match. Nitrogen:
all the arginine nitrogen leaving ``amino_acids`` (``aa_mass·y_N``) lands in the ``N`` pool
(methanethiol is nitrogen-free), so ``total_nitrogen`` is unchanged. Both to machine precision,
since ``mercaptans`` is weighted in ``total_carbon`` (as methanethiol) and ``amino_acids``/``N`` are
already on both ledgers (D-32).

**TIER — a new structural drop on ``N`` (the D-27 ``E`` parallel, advisor-flagged).**
:class:`AutolyticMercaptan` is the **first autolysis-gated Process to write ``N``** (via the
deamination). So an autolysis-on run drops the *structural* ``tier_of("N")`` PLAUSIBLE→SPECULATIVE
— even on an autolysis-on / amino-acid-dose-*off* run, where the other N-writer
(:class:`~fermentation.core.kinetics.byproducts.FuselAminoAcidReroute`) stays disabled. The
param-aware tier users see was typically already speculative (growth reads speculative params), so
no headline change, and the ``N`` pool genuinely does carry a speculative deamination trace when
autolysis is opted in. Tier **speculative** overall (the yield magnitude is an author estimate).

**Isolability — opt-in and wine-only, with the D-34 autolysis Processes.** Touches only
``mercaptans``/``amino_acids``/``N`` (all wine slots), so it is wired into wine alone and disabled
**together with** :class:`~fermentation.core.kinetics.autolysis.YeastAutolysis` and
:class:`~fermentation.core.kinetics.hydrogen_sulfide.AutolyticHydrogenSulfide` at the compile seam
absent ``autolysis_rate_per_h`` — an undosed wine run is byte-for-byte the validated core and the
empty ``mercaptans`` slot keeps its tier. Guards: no dead biomass (``autolysis_flux ≤ 0``) or an
empty amino-acid pool ⇒ 0. Copper fining removes it (the ``add_copper`` verb, extended in D-45 to
bind mercaptans as copper mercaptide Cu(SR)₂).
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import carbon_mass_fraction
from fermentation.core.kinetics.amino_acid_pools import (
    SPEC_BY_SPECIES,
    depletion_gate,
    draw_precursor_carbon,
)
from fermentation.core.kinetics.autolysis import autolysis_flux
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The representative species for the lumped **mercaptans** thiol pool (decision D-45): methanethiol
#: (methyl mercaptan), the dominant reduction thiol. Named as a module constant so the carbon draw
#: here and the ``total_carbon`` weighting name one species (the D-19 idiom).
_MERCAPTAN_SPECIES = "methanethiol"

#: Methanethiol's actual precursor (decision D-100). D-45 had to draw the thiol's carbon from the
#: lumped ``amino_acids`` pool (arginine) and document the mismatch as a provenance caveat —
#: arginine contains no sulfur, so the molecule it books could not possibly make a mercaptan. With
#: the pool speciated, the draw finally *is* methionine: the sulfur-bearing amino acid whose
#: demethiolation genuinely releases methanethiol. **The D-45 caveat is retired, not restated.**
_PRECURSOR_SPECIES = "methionine"


class AutolyticMercaptan(Process):
    """Autolytic mercaptan (thiol) release — a carbon-bearing yield on the autolysis flux (D-45).

    ``r_merc = y_mercaptan · autolysis_flux(y) · [aa/(K_amino_acids+aa)]`` [g methanethiol/L/h];
    fills the lumped ``mercaptans`` pool, drawing the mercaptan carbon from the ``amino_acids`` pool
    and **deaminating** its nitrogen back to ammonium ``N`` (Option A, the D-33 fusel-reroute
    idiom). Carbon closes (C into ``mercaptans`` = C out of ``amino_acids``); nitrogen closes (arg N
    → ``N``, methanethiol is N-free). Not flux-linked, so it accumulates un-stripped post-dryness —
    the reductive fault. Opt-in and wine-only (rides the D-34 autolysis gate). The **first
    autolysis-gated ``N``-writer**, so it drops the structural ``tier_of("N")`` to speculative (the
    D-27 ``E`` parallel); tier speculative.
    """

    name = "autolytic_mercaptan"
    tier = Tier.SPECULATIVE
    #: Fills ``mercaptans``, debits ``methionine`` for its carbon (decision D-100 — the *actual*
    #: precursor, no longer the arginine lump), releases the nitrogen to ``N``.
    touches = ("mercaptans", _PRECURSOR_SPECIES, "N")
    #: ``y_mercaptan`` sets the g-MeSH-per-g-biomass-autolysed yield; ``k_autolysis``/
    #: ``E_a_autolysis``/``T_ref`` are the *same* autolysis constants
    #: :func:`~fermentation.core.kinetics.autolysis.autolysis_flux` reads (so all autolysis branches
    #: share one clock and one ``autolysis_rate_per_h`` override); ``K_amino_acids`` scaled by
    #: ``must_aa_fraction_methionine`` sets the relative-depletion gate (D-100). Their tiers cap the
    #: output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "y_mercaptan",
        "k_autolysis",
        "E_a_autolysis",
        "T_ref",
        "K_amino_acids",
        SPEC_BY_SPECIES[_PRECURSOR_SPECIES].fraction_param,
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        r_autolysis = autolysis_flux(y, schema, params)  # [g X_dead/L/h] — the shared D-34 flux
        if r_autolysis <= 0.0:
            return d  # no dead cells ⇒ nothing to autolyse (clamped, so no negative overshoot)
        # Methionine's own relative-depletion gate (decision D-100): → 0 as the pool empties, so
        # the draw can never drive it negative and an undosed wine is the byte-for-byte no-op.
        gate = depletion_gate(y, schema, params, (SPEC_BY_SPECIES[_PRECURSOR_SPECIES],))
        if gate <= 0.0:
            return d  # no methionine ⇒ no thiol source ⇒ no mercaptan (the D-33 no-op)

        r_merc = params["y_mercaptan"] * r_autolysis * gate  # [g methanethiol/L/h]
        # Draw the mercaptan carbon from METHIONINE and deaminate its nitrogen (Option A, D-33;
        # speciated at D-100): the carbon into mercaptans is sized to equal the carbon out of
        # methionine, so carbon closes; the released methionine nitrogen all goes to the N pool
        # since methanethiol is nitrogen-free, so nitrogen closes — both by construction.
        merc_carbon = r_merc * carbon_mass_fraction(_MERCAPTAN_SPECIES)  # [g C/L/h] in the thiol
        nitrogen = draw_precursor_carbon(d, schema, _PRECURSOR_SPECIES, merc_carbon)

        d[schema.slice("mercaptans")] = r_merc
        d[schema.slice("N")] = nitrogen  # DEAMINATION: methionine nitrogen → ammonium (D-33)
        return d
