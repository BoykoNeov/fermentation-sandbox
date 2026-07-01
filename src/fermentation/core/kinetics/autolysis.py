"""Yeast autolysis — the autolytic-peptide source that refills the amino-acid pool (D-34).

**What this closes (decisions D-23 → D-32 → D-34).** MLF with bacterial *growth* is blocked on
two prerequisites. The fusel Ehrlich re-route (D-33) closed the first. This closes the second:
*Oenococcus oeni* builds biomass from amino acids/peptides, but the ``amino_acids`` pool (D-32)
is **empty at the MLF pitch point** — the same yeast uptake that strips ``N`` to ~0 by day ~1.3
would strip any dosed amino acids too (the empirical finding that settles D-23). Something must
*refill* the pool post-AF. Real wine does this by **autolysis**: as the ferment ends, yeast die
(the ``X_dead`` pool the D-13 ethanol-inactivation Process fills) and slowly self-digest,
releasing intracellular amino acids and peptides into the wine — the basis of *sur lie* aging.
:class:`YeastAutolysis` is that flux: the **first consumer of ``X_dead``**, turning dead biomass
into assimilable ``amino_acids`` (feeding the D-32/D-33 consumers and, later, MLF-growth).

**The conservation problem, and why a debris pool (advisor-decided).** Dead biomass is
**carbon-rich** (mass C:N ``f_C/f_N`` ≈ 4–11 across Coleman's nitrogen range) while the
assimilable amino acids it releases are **nitrogen-rich** (arginine mass C:N ≈ 1.29). So per gram
of nitrogen liberated, biomass gives up 4–11 g of carbon but arginine can only hold ~1.3 g — **most
of the dead-cell carbon cannot leave as amino acids.** Physically that carbon is the yeast cell
wall (β-glucans / mannoproteins, ~30 % of dry mass), which is exactly the non-assimilable material
that stays as lees. Routing it to CO₂ would falsely claim autolysis *respires* the cell (it is
enzymatic self-digestion, not respiration) and would perturb a benchmarked pool. Instead the excess
carbon goes to a new carbon-only **``debris``** pool (booked as glucan; :data:`~fermentation.core.\
chemistry.M_GLUCAN`), the honest and physically-dominant fate — the ``esters_gas`` precedent (a
bookkeeping pool weighted in ``total_carbon`` only, carrying carbon that left the metabolite ledger
but not the atom balance).

**The flux — nitrogen-anchored, first-order in dead biomass.** With
``r = k_autolysis · f_T · X_dead`` [g X_dead/L/h] (``f_T = arrhenius_factor(T, E_a_autolysis,
T_ref)`` — autolysis is enzymatic, so warmer lees clear faster):

  * **liberate the dead-cell nitrogen as amino acids** — ``d[amino_acids] = +r·f_N/y_N`` (arginine
    carrying exactly the nitrogen ``r`` releases, ``r·f_N``),
  * **debit dead biomass** — ``d[X_dead] = −r``, and
  * **route the excess carbon to debris** — ``d[debris] = +(r·f_C − r·f_N·y_C/y_N) / c_debris``.

Carbon closes: the dead-cell carbon ``r·f_C`` splits into the amino acids' ``r·f_N·y_C/y_N`` and the
debris' remainder. Nitrogen closes: the dead-cell nitrogen ``r·f_N`` is exactly what the amino-acid
pool gains (debris is nitrogen-free). Both to machine precision, since the amino-acid and debris
pools are weighted in ``total_nitrogen``/``total_carbon``. The excess-carbon term is **structurally
non-negative** — biomass C:N always exceeds arginine's across the whole ``f_N`` range (0.039–0.114),
so ``f_C > f_N·y_C/y_N`` always and the split never flips (no clamp, no C⁰ kink for the BDF solver).

**Isolability — opt-in (the D-30 carrying-capacity pattern).** Unlike the always-on intrinsic aroma
pools, turning autolysis on *consumes* ``X_dead`` and fills ``amino_acids``/``debris`` — it
measurably perturbs the core, so it cannot be default-on without breaking the validated-core
byte-for-byte guarantee and the §2.2 benchmarks. It ships **wine-only and disabled by default**: the
compile seam enables it only when a scenario passes ``autolysis_rate_per_h`` (which also overrides
``k_autolysis``, letting a demonstration sweep the timescale). Disabled ⇒ excluded from the Process
set's derivatives and tier derivation, so an undosed wine run is byte-for-byte the validated core.
Its first guard is ``X_dead ≤ 0 ⇒ 0`` (no dead cells ⇒ nothing to autolyse; first-order in the
clamped pool cannot
overshoot negative). Wine-only, mirroring the wine-only ``amino_acids`` pool and nitrogen model
(D-30/D-32); beer deferred.

Tier: **speculative** — first-order autolysis of dead biomass is a standard lumped form, but
``k_autolysis`` and ``E_a_autolysis`` are author estimates and the single-amino-acid /
carbon-only-debris lumping is a simplification (real autolysate is a mix of amino acids, peptides,
nucleotides and mannoproteins, the last retaining some nitrogen). Parameter-tier propagation (D-1)
caps the ``X_dead``/``amino_acids``/``debris`` outputs accordingly.

SCOPE (v1): the refill flux only. The eventual *consumer* — an MLF-with-growth Process feeding a
growing ``X_mlf`` from this pool, plus the event loop that pitches bacteria post-AF — stays deferred
(decision D-23); this beat makes the pool it will draw on non-empty. A standalone excess-amino-acid
deamination flux (vs the D-33 fusel-coupled release) also remains future work.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The representative species for the non-assimilable cell-wall **debris** pool (decision D-34):
#: an anhydroglucose (glucan) repeat unit. Carbon-only — all the released nitrogen leaves as amino
#: acids, so the C-rich remainder carries none. Named as a module constant so the excess-carbon
#: routing here and the ``total_carbon`` weighting name one species (the D-19 idiom).
_DEBRIS_SPECIES = "glucan"


class YeastAutolysis(Process):
    """Autolysis of dead biomass into assimilable amino acids + cell-wall debris (decision D-34).

    ``r = k_autolysis·arrhenius(T, E_a_autolysis, T_ref)·X_dead`` [g/L/h]; liberates the dead-cell
    nitrogen as amino acids (``d[amino_acids] = +r·f_N/y_N``), debits dead biomass
    (``d[X_dead] = −r``), and routes the C-rich remainder to ``debris``. Carbon and nitrogen
    close by construction (the module docstring); the excess-carbon split is structurally
    non-negative (biomass is always more carbon-rich than arginine). Opt-in and wine-only — the
    first consumer of ``X_dead``, refilling the ``amino_acids`` pool post-AF for the deferred
    MLF-with-growth beat.
    """

    name = "yeast_autolysis"
    tier = Tier.SPECULATIVE
    #: Consumes dead biomass ``X_dead``; deposits amino acids and the carbon-only ``debris`` pool.
    touches = ("X_dead", "amino_acids", "debris")
    #: ``k_autolysis`` sets the rate, ``E_a_autolysis``/``T_ref`` its temperature shape;
    #: ``biomass_N_fraction``/``biomass_C_fraction`` (the same fractions the conservation checks
    #: read, D-8) partition the released mass between amino acids and debris. Their tiers cap the
    #: ``X_dead``/``amino_acids``/``debris`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_autolysis",
        "E_a_autolysis",
        "T_ref",
        "biomass_N_fraction",
        "biomass_C_fraction",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        x_dead = max(float(y[schema.slice("X_dead")][0]), 0.0)
        if x_dead <= 0.0:
            return d  # no dead cells ⇒ nothing to autolyse (clamped, so no negative overshoot)

        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_autolysis"], params["T_ref"])
        r = params["k_autolysis"] * f_t * x_dead  # [g X_dead/L/h] autolysed

        f_n = params["biomass_N_fraction"]
        f_c = params["biomass_C_fraction"]
        y_n = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
        y_c = carbon_mass_fraction(AMINO_ACID_SPECIES)
        c_debris = carbon_mass_fraction(_DEBRIS_SPECIES)

        # Nitrogen-anchored: the dead-cell nitrogen r·f_N leaves as amino acids (arginine carrying
        # exactly that nitrogen). The amino acids' carbon is r·f_N·(y_C/y_N); the rest of the
        # dead-cell carbon (r·f_C) is the non-assimilable cell-wall remainder → debris. Structurally
        # f_C > f_N·y_C/y_N (biomass C:N ≈ 4–11 ≫ arginine's ≈ 1.29), so the excess is always
        # positive — no clamp, no derivative kink (advisor; decision D-34).
        d[schema.slice("X_dead")] = -r
        d[schema.slice("amino_acids")] = r * f_n / y_n
        debris_carbon = r * (f_c - f_n * y_c / y_n)  # [g C/L/h] cell-wall carbon left behind
        d[schema.slice("debris")] = debris_carbon / c_debris
        return d
