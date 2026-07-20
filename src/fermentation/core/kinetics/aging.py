"""Aging chemistry — the slow, post-fermentation "years" axis (§4.1, decision D-69).

The first Tier-3 aging Process, opened in D-68 and built here. Unlike the Milestone-2
byproduct Processes (which *produce* aroma pools during active fermentation), an aging
Process acts on a *finished* wine/beer over months-to-years: fermentation is done, the
sugar is gone, the yeast racked or crashed, and the chemistry that remains is spontaneous
(hydrolysis, oxidation, condensation), not metabolic. This module holds the aging
Processes; :class:`EsterHydrolysis` is the first, :class:`OxidativeAcetaldehyde` the second,
:class:`SulfiteOxidation` the third, :class:`PhenolicBrowning` the fourth,
:class:`StreckerDegradation` the fifth, :class:`OakExtraction` (D-77) the sixth — the first
**non-oxidative** aging Process and a **separate axis** (barrel/chip aroma extraction, drawing no
O₂, orthogonal to the whole oxidative sub-axis below) — :class:`EllagitanninOxidation` (D-78) the
seventh (the oak-tannin O₂ sink that bridges the oak axis to the O₂ sub-axis), and
:class:`TanninAnthocyaninCondensation` (D-79) the eighth — the second **non-oxidative** Process and
a third **separate axis**: grape anthocyanin + grape tannin condense into stable polymeric pigment
(red-wine colour stabilization + astringency softening), drawing **neither O₂ nor oak** (see its
docstring for why it is oak- *and* O₂-independent) — and :class:`AcetaldehydeBridgedCondensation`
(D-80) the ninth: the **split-ledger** beat D-79 deferred, where dissolved-O₂ acetaldehyde (D-71)
bridges grape tannin to anthocyanin, the *first link from the oxidative sub-axis to red-wine
colour*.
It is the first aging colour Process to touch the **carbon ledger** — a new on-ledger
``ethyl_bridge``
slot captures the acetaldehyde carbon so it does not vanish into the off-ledger grape pigment (the
grape bulk off-ledger, the acetaldehyde-derived bridge on it — the "split ledger"; see its
docstring).

**The non-oxidative THERMAL axis (D-87/D-88/D-89).** :class:`MaillardStrecker` (D-87) is the
**O₂-independent thermal mirror** of :class:`StreckerDegradation` (D-75) — the beat D-75 deferred:
residual **sugar + heat** (α-dicarbonyls, no O₂) degrade amino acids to the sweet-wine / Madeira
Strecker suite (the two D-75 aldehydes + three branched-chain malty aldehydes + sotolon), sharing
the ``amino_acids`` limiting reagent with the oxidative route. :class:`Caramelization` (D-88) is the
matching **O₂-independent thermal mirror** of :class:`PhenolicBrowning` (D-74): **sugar-only**
browning to a new on-ledger ``melanoidin`` carbon-park pool, raising the shared ``A420`` with no O₂.
:class:`MaillardBrowning` (D-89) is the **amino-acid-incorporating** browning branch D-88 deferred:
**sugar + amino acid + heat** → a nitrogen-bearing ``maillard_melanoidin`` pool (the FIRST
non-biomass, non-arginine species on ``total_nitrogen``), retaining the amino-acid nitrogen in the
polymer and raising the same ``A420``. Together the three split the thermal amino-acid/browning fate
cleanly: **sugar-only** caramelization (D-88), the **N-retaining** melanoidin branch (D-89), and the
**N-releasing** volatile Strecker branch (D-87) — summed over the shared ``S``/``amino_acids``
reagents by ``ProcessSet``, no double-count (see :class:`MaillardBrowning`). All driven by residual
sugar (D-88/D-89 *consume* core ``S``), all the same thermal-mirror relationship
:class:`ThermalAnthocyaninFade` (D-83) has to :class:`AnthocyaninFading` (D-81). D-87/D-89 are
**wine-only** (they read the wine-only ``amino_acids`` pool, untracked in beer, D-32); D-88
:class:`Caramelization` is **medium-agnostic** (D-90) — sugar-only browning is physics, so beer's
residual dextrins caramelize too (the vectorized draw apportions across beer's 3-slot ``S``).

**The oxidative sub-axis (D-71).** :class:`OxidativeAcetaldehyde` opens the *oxidative* half of
the aging axis on a **dissolved-O₂ pool** (``o2``, a new carbon-free state slot, off every
conservation ledger like ``h2s``/``iso_alpha``). O₂ — not ethanol — is the rate-limiting
reactant (ethanol sits at ~100 g/L, effectively constant across aging), so the pool is *the*
substrate that bounds oxidation: acetaldehyde saturates as the O₂ charge is spent, the
bottle-aging reality a first-order-in-ethanol rate could never reproduce (it would grow
unbounded). O₂ enters via a dedicated ``add_oxygen`` dosing verb (one dose = a bottle's
ingress; repeated = micro-oxygenation / barrel), and a ``begin_aging`` run with **no** O₂ dosed
is purely *reductive* aging (screwcap/inert) — byte-for-byte the :class:`EsterHydrolysis`-only
aging, since the Process contributes exactly zero at ``o2 = 0``. **(D-88/D-90 supersede this for
residual-sugar beverages:** the "un-oxygenated aging is byte-for-byte the ester-only case" claim now
holds only for *dry* ones — a sealed **sweet** wine or **high-residual** beer browns thermally
(:class:`Caramelization`, medium-agnostic per D-90) with no O₂ at all, and a sweet wine additionally
develops the thermal Strecker suite (:class:`MaillardStrecker`, wine-only); the D-83-style mirror.)
Oxidative aging is
fundamentally a competition for a finite O₂ budget: the ``o2`` pool is the shared substrate the
whole oxidative sub-axis draws down, and **each O₂ consumer owns its own rate constant and draws
its own share** — ``ProcessSet`` sums them, so the pool depletes *once* and the O₂ splits among
the sinks by ``kᵢ / Σk`` (the additive pattern :class:`SulfiteOxidation` established at D-72,
extended to *always-on* sinks at **D-73**). So :class:`OxidativeAcetaldehyde` consumes only its
**ethanol-oxidation share** (``k_ethanol_oxidation``), not the whole flux, and
``y_acetaldehyde_per_o2`` is the *true* per-O₂ stoichiometric yield of that route alone — the
competition lives in the rate constants, not in a sub-unity yield. :class:`PhenolicBrowning`
(**D-74**) is the first always-on sink to land in that seam: the *dominant* O₂ consumer
(``k_browning``), it diverts most of the always-on flux to brown pigment and correspondingly
suppresses oxidative acetaldehyde — the reduction of ``k_ethanol_oxidation`` (5.0e-4 → 2.0e-4, so
``k_ethanol + k_browning`` holds the calibrated total O₂-depletion rate) that D-73 made possible and
D-71 could not express under "total rate". A further Strecker-degradation sink would slot in the
same way (D-73 reworked D-71's original "whole-flux / total-rate" framing so an always-on sink adds
cleanly, without double-counting).

**The first of those sinks: SO₂ scavenging (D-72).** :class:`SulfiteOxidation` is the first
sibling to claim its share of that ``o2`` budget. Dissolved O₂ oxidises free **bisulfite** (the
reactive antioxidant HSO₃⁻ — *not* molecular SO₂, which is the antimicrobial form) to sulfate,
so — competing for the same ``o2`` pool via ``ProcessSet`` summing — SO₂ diverts O₂ away from
ethanol oxidation: **while free SO₂ lasts, oxidative acetaldehyde is suppressed; once it is spent,
acetaldehyde climbs** (the classic wine threshold, emergent, nothing extra built). It decrements
the existing ``so2_total`` slot (no new pool) at the Danilewicz 2:1 mol SO₂:O₂ stoichiometry, and
self-throttles as D-47 acetaldehyde–SO₂ binding erodes the free pool.

**The first *always-on* sink: phenolic browning (D-74).** :class:`PhenolicBrowning` is the second
oxidative sibling and the first **always-on** claimant of the ``o2`` budget (SO₂ oxidation is
*substrate-gated* — zero without dosed SO₂; browning runs whenever O₂ is present). As a finished
wine/beer takes up O₂, dissolved O₂ oxidises **phenolics** (o-diphenols → o-quinones, which
polymerise to brown melanoidin/quinone pigment) — the gold→amber→brown of an aged/maderised white
wine, and oxidative darkening in beer. This is in fact the *dominant* O₂ consumer in wine oxidation
(phenol autoxidation is the primary O₂ sink; ethanol oxidation to acetaldehyde is a *secondary*
H₂O₂ fate), so ``k_browning`` is set the **larger** share and — competing for the same ``o2`` pool
via ``ProcessSet`` summing, exactly like SO₂ — it **diverts most of the always-on O₂ away from
ethanol oxidation, suppressing oxidative acetaldehyde**: the always-on analogue of SO₂'s protection
(SO₂ suppresses *until spent*; browning suppresses *permanently*, as a co-resident sink). Unlike the
other three aging Processes it has **no aroma product the D-67 OAV lens already reads**, so — the
D-68 selection criterion inverted — it needs **one new observable**: the ``A420`` browning index (an
optical **absorbance** at 420 nm, dimensionless AU, **not** a pigment mass). ``A420`` is a *state
slot* (the ``iso_alpha`` off-ledger-pool pattern, **not** the D-67 post-hoc OAV series): browning
pigment is **cumulative and irreversible** and its O₂ flux is *dynamic* (SO₂ competes, temperature
varies), so it must be **integrated** along the run — it cannot be reconstructed post-hoc from
(dosed − remaining) O₂. Because ``A420`` is an optical index rather than a mass, its carbon (which
would come from an *untracked* phenol pool) is sidestepped by construction: this Process touches
only ``{o2, A420}`` — **both off every ledger** — so it moves *nothing* conserved (cleaner even than
:class:`OxidativeAcetaldehyde`, which still borrows carbon E→acetaldehyde). ``d(A420)/dt ≥ 0``
always (monotonic; no clamp). Medium-agnostic (both media carry autoxidising polyphenols and brown
oxidatively — D-74 supersedes D-73's provisional "wine-only" parenthetical), so wired into both
media like :class:`OxidativeAcetaldehyde`.

**Off during the ferment, on during an aging segment (D-68/D-70).** These Processes ARE wired
into both media's ProcessSet (D-70) but **disabled at the compile seam** — a ``begin_aging``
scheduled event enables them for a long post-fermentation segment (the ``simulate_scheduled``
reconfigure mechanism, the ``pitch_mlf`` pattern minus the state mutation). So the validated
core and the Milestone-2 aroma beat stay byte-for-byte isolable (prime directive #3): a
compiled scenario with no ``begin_aging`` never activates this module (disabled ⇒ skipped by
``active``/``tier_of``/strict), and building a bare ProcessSet without this tuple *is* the
pre-aging model. During a post-dryness aging segment every OTHER producer of
``esters``/``fusels``/``Byp`` is fermentative-flux-gated and quiescent at ``S ≈ 0``, so the
aging ester/fusel signal is unconfounded — only :class:`EsterHydrolysis` moves those pools
(Stance A, D-70). The scenario-level ``begin_aging`` verb + span-via-``duration_days`` + the
§7 slow-phase integration (the segment restart lets the solver take large steps across the
quiescent aging segment) landed in D-70; in D-69 the Process was exercised directly via a
hand-built ``ProcessSet`` (the D-64 loss-Process test pattern).

----------------------------------------------------------------------------------------

**:class:`EsterHydrolysis` — young fruity esters fade with age.**

The chemistry: the acetate esters that dominate the young ``esters`` pool (isoamyl acetate,
the "banana" ester D-67 reads as the pool's sensory representative; ethyl acetate) form
*above* their hydrolysis equilibrium during fermentation and slowly **hydrolyse back toward
equilibrium** as the wine ages (Ramey & Ough 1980; Marais 1978). The isoamyl-acetate
hydrolysis this Process models is::

    isoamyl acetate + H₂O  →  isoamyl alcohol (a fusel)  +  acetic acid (a volatile acid)

so aging **fades the ester OAV**, **raises the fusel OAV**, and drifts **VA/pH** up — all
real, coherent aging phenomena, and the D-68 reason this was chosen as the first aging
Process (it moves OAVs the D-67 lens already reads, needing no new extraction driver and no
new state pool). The rate::

    d(esters)/dt = -k_ester_hydrolysis · f(T) · max(0, esters - isoamyl_acetate_eq)

is **net decay toward a lower equilibrium floor** ``isoamyl_acetate_eq``, *not* decay-to-zero
(D-68):
below ``isoamyl_acetate_eq`` the rate is zero. Framing this as "net decay toward a lower floor"
is the same fixed-composition honesty the D-67 sensory lump carries. ``f(T) = arrhenius_factor(T,
E_a_ester_hydrolysis, T_ref)`` gives the sourced **warmer-ages-faster** direction (cold
cellars preserve fruity esters). No fermentative-flux gate — aging runs when the flux is
zero — so unlike the M2 producers this Process is driven by temperature and the pool alone.

**The "deferred half" this docstring used to name DOES NOT EXIST TO BUILD (decision D-121).**
Until D-121 the text above claimed the bidirectional reality was that *"ethyl esters of fatty
acids sit below equilibrium young and slowly form on aging"*. **That is directionally wrong**,
and it conflated two families that move in *opposite* directions:

* **Straight-chain medium-chain fatty-acid ethyl esters** (yeast *lipid* metabolism — ethyl
  hexanoate/octanoate/decanoate) are secreted **above** equilibrium and **hydrolyse** on aging:
  the **same** direction as the acetates, not the opposite one.
* **Ethyl esters of branched and polyprotic/organic acids** (yeast *nitrogen* metabolism and
  chemical esterification — ethyl isobutanoate, ethyl 2-methylbutanoate, ethyl isovalerate,
  diethyl succinate, ethyl lactate) are the family that sits **below** equilibrium young and
  **forms**. The sim tracks **none** of them, so the deferred formation half as scoped had no
  pool to act on.

Sources for the DIRECTION SPLIT are three independent verbatim reads (Díaz-Maroto/Schneider/
Baumes 2005 SIDA on 1–5 y Muscadet; Makhotkina & Kilmartin 2012 on Sauvignon blanc at
5/10/18 °C; Marais & Pool 1980) — Ramey & Ough 1980 does not speciate the forming family, so
there it is cited for mechanism only. **The RATE is a different claim, and at D-123 Ramey & Ough
1980 supplied it:** the paper (paywalled at ACS but readable as an open scanned PDF on the
author's winery site) measures the isoamyl-acetate hydrolysis rate *in real wine* — Table IX/VIII
``k_obsd`` and Table X ``E_a`` — which now anchor ``k_ester_hydrolysis`` /
``E_a_ester_hydrolysis`` (previously author estimates). D-121's "canonical source is paywalled,
never the evidence" is superseded for the rate; it stands for the direction split.

``ethyl_acetate`` is the **third** case and the one exception among the acetates: it
**increases** on storage (ethanol at ~2 M plus accumulating acetic acid put it *below* its
esterification equilibrium). It is deliberately left inert here — see D-121/D-123 for why neither
it nor ``ethyl_hexanoate`` ships a term: for ``ethyl_acetate`` the direction is known but the
**rate is not** (it is absent from R&O's eight esters, and no canonical study covers a 51 mg/L
pool), and ``ethyl_hexanoate`` has a R&O *model-solution* k (Table I) but their own real-wine
data show it does not change appreciably — so a term would invent motion the source denies. Both
stay **blocked on sourcing**, not refused.

**Carbon — an on-ledger inter-pool transfer (conservation is back in force, D-68).** Unlike
the D-67 sensory readout (a pure diagnostic off the ledger), this is the first aging RHS
*on the carbon ledger*, so the carbon a decaying ester releases **must** be routed. The
carbon leaving ``esters`` per unit decayed is **ledger-fixed** at ``rate·c(ethyl_acetate)``
— the pool's D-19 mass weighting (ethyl acetate, C4) is immovable and untouched here. That
released-carbon budget is split between the two products and re-deposited via each product
pool's own carbon fraction, so ``total_carbon`` closes to machine precision for *any* split
summing to 1 (the split only re-partitions a fixed budget between two trace pools — it is
second-order on outputs; the ``esters → esters_gas`` transfer precedent, but C4 → C5-partial
+ C4-partial across two differently-weighted pools).

*Carbon is the invariant; mass carries a small documented gap.* Splitting a carbon-exact
budget across pools with **heterogeneous, fixed** mass weightings (ester debited as ethyl
acetate, products credited as isoamyl alcohol / succinic) is not mass-conserving: per gram of
ester decayed the products sum to ~0.955 g, a ~4.5 % gap. Real hydrolysis in fact *gains* mass
(it consumes water); this bookkeeping loses a little. That is the standard stand-in gap already
on record — beer's hydrolysis-water gap (D-8), the ``Gly``/``Byp`` redox-diversion gap that
scopes ``total_mass`` to ``{S, E, CO2}`` (D-16), and the VDK O₂/NAD(P)H gap (D-26): carbon has
no such term (water/redox H,O carry none) and is the rigorous invariant. ``total_mass`` does not
even see this transfer (it weights only ``{S, E, CO2}``, none of which this Process touches), so
the gap is scoped out by construction, not silently violated.

**The split ratio is 5:2 — fusels : Byp (decision D-69, the advisor-settled crux).** D-68
surfaced the choice and deferred it to this build's advisor pass. The pool's *mass* is
ledger-weighted as ethyl acetate (immovable, D-19), but the *split* of the released carbon
between the alcohol and acid products is the one free variable, and it is set by the
**isoamyl-acetate** stand-in reaction above: isoamyl alcohol carries **5** carbons, acetic
acid **2** — hence 5/7 of the carbon to ``fusels``, 2/7 to ``Byp``. This was chosen over the
ethyl-acetate-consistent 1:1 (ethyl 2C : acetyl 2C) because **this is a sensory Process**:
its entire reason to exist (D-68) is to fade the ester OAV and raise the fusel OAV, and D-67
already commits ``esters``' OAV to **isoamyl acetate** and ``fusels``' to **isoamyl
alcohol**. The coherent chemistry connecting those two committed representatives is exactly
isoamyl acetate → isoamyl alcohol + acetic acid = 5:2; the 1:1 alternative would route
*ethanol*-carbon (ethyl acetate's actual alcohol) into the isoamyl-alcohol-weighted
``fusels`` pool and read it through the isoamyl-alcohol OAV — fabricating the fusel-aroma
rise out of the wrong molecule, a bend on the exact quantity the Process exists to move. The
5:2 cost is narrative-only and invisible to every conservation test: the *debited* molecule
(ethyl acetate) and the *split* molecule (isoamyl acetate) differ — a documented stand-in
seam this Process **inherits** from the D-19/D-67 ethyl-acetate-mass / isoamyl-acetate-aroma
mismatch, not one it invents.

**§4.3 firewall tension — documented, owner-accepted (D-68 fork 2).** This speculative-tier
Process touches ``Byp`` (the acetic-acid product, booked as the succinic-acid stand-in),
which the *plausible*-tier pH/TA charge-balance readout reads — so a speculative aging
Process nudges a plausible output. The owner explicitly accepted the literal chemistry (the
VA/pH drift is a real aging phenomenon) in D-68 fork 2. Isolability is preserved regardless
(disable the Process and the drift vanishes). ``Byp`` is a succinic (C4 diprotic) stand-in
for acetic acid (C2 monoprotic) — the same D-16 bookkeeping stand-in the pool already is.

Tier: **speculative** — the aging axis is the Tier-3 speculative frontier; the hydrolysis
*form* is sourced (first-order approach to equilibrium, warmer-faster) but the magnitudes are
order-of-magnitude estimates. Parameter-tier propagation (D-1) caps the pool outputs at
speculative regardless of input tiers.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.acidbase import (
    SO2_STATE_KEY,
    bisulfite_so2_at_ph,
    free_acetaldehyde,
    ph_of_state,
)
from fermentation.core.chemistry import (
    CARBON_ATOMS,
    M_2_METHYLBUTANAL,
    M_2_METHYLPROPANAL,
    M_3_METHYLBUTANAL,
    M_ACETALDEHYDE,
    M_ALPHA_KETOBUTYRATE,
    M_CO2,
    M_ETHANOL,
    M_METHIONAL,
    M_O2,
    M_PHENYLACETALDEHYDE,
    M_SO2,
    M_SOTOLON,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
    sugar_species,
)
from fermentation.core.kinetics.amino_acid_pools import (
    ASSIMILABLE_SPECS,
    SPEC_BY_SPECIES,
    assimilable_carbon_per_nitrogen,
    depletion_gate,
    draw_assimilable_nitrogen,
    draw_precursor_carbon,
)
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import (
    ACETYL_CARBON_SHARE,
    ACETYLATION_ACETYL_CARBONS,
    ACETYLATION_ALCOHOL_CARBONS,
    ALCOHOL_CARBON_SHARE,
    HYDROLYSING_ESTER,
    ISOAMYL_ALCOHOL,
    VALINE_LABEL_TRACERS,
    labelled_fraction,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Representative species that carbon-accounts the acid *product* pool (``Byp``) the hydrolysis
#: fills — succinic acid (D-16), from the one chemistry source of truth. Neither of the other
#: two molecules is named here any more: the *debited* ester comes from ``HYDROLYSING_ESTER``
#: (since D-96 the real molecule, isoamyl acetate, rather than a stand-in), and the *alcohol
#: product* from ``ISOAMYL_ALCOHOL`` (since D-99 its own pool rather than a share of a lump).
#: Using these fractions both to release and to re-deposit the carbon is what makes the
#: transfer close in ``total_carbon`` exactly.
_BYP_SPECIES = "succinic_acid"

#: The 5:2 carbon split of the released ester carbon between ``fusels`` and ``Byp``
#: (decisions D-69/D-96), set by the hydrolysis reaction **isoamyl acetate → isoamyl alcohol +
#: acetic acid**: isoamyl alcohol carries 5 carbons, acetic acid 2. Stoichiometry of a named
#: reaction — a code-with-citation constant like the chemistry carbon counts, not an
#: empirical/uncertain YAML parameter.
#:
#: At D-69 this was a *stand-in* ratio: the pool it debited was ledger-weighted as ethyl acetate
#: (C4), so 5 + 2 = 7 could not equal the debited molecule's carbon and the mismatch had to be
#: documented instead of checked. Since D-96 the Process debits isoamyl acetate itself, so the
#: split **exactly partitions the debited molecule's seven carbons** — an invariant the module
#: now asserts at import rather than asserting in prose, and one that would fire if a future
#: edit re-pointed the hydrolysis at an ester whose stoichiometry these shares do not describe.
#: **Since D-115 the ratio lives in the registry, not here.** The synthesis side now runs the
#: same reaction in reverse (the 5:2-inverse re-route), so two modules read one stoichiometry —
#: the D-26/D-106 setup where private copies agree by luck until one of them changes. These are
#: aliases onto the shared constants, kept under their original names so the algebra below reads
#: as it always has.
_ISOAMYL_ALCOHOL_CARBONS = ACETYLATION_ALCOHOL_CARBONS  # the alcohol product → isoamyl_alcohol
_ACETIC_ACID_CARBONS = ACETYLATION_ACETYL_CARBONS  # the acid product → Byp
_FUSEL_CARBON_SHARE = ALCOHOL_CARBON_SHARE
_BYP_CARBON_SHARE = ACETYL_CARBON_SHARE

#: The two D-115 label tracers, resolved once. ``VALINE_LABEL_TRACERS`` is ordered
#: alcohol-then-ester (the direction the label flows through the acetylation); this Process runs
#: that flow backwards, crediting the alcohol from the ester.
_ALCOHOL_TRACER, _ESTER_TRACER = VALINE_LABEL_TRACERS

if CARBON_ATOMS[HYDROLYSING_ESTER.species] != _ISOAMYL_ALCOHOL_CARBONS + _ACETIC_ACID_CARBONS:
    raise AssertionError(  # pragma: no cover - structural invariant, D-96
        f"The {_ISOAMYL_ALCOHOL_CARBONS}:{_ACETIC_ACID_CARBONS} hydrolysis split must "
        f"partition every carbon of {HYDROLYSING_ESTER.species} "
        f"({CARBON_ATOMS[HYDROLYSING_ESTER.species]} C)."
    )

#: Ethanol and acetaldehyde are both C2, so the oxidative ``ethanol → acetaldehyde`` transfer is
#: mole-for-mole — exactly the (inverse of the) D-27 acetaldehyde reduction. Weighting the ethanol
#: debit by this molar-mass ratio makes the borrow carbon-exact: ``M_acet·cf_acet == M_eth·cf_eth·
#: (M_acet/M_eth) == 2·M_C`` per mole, so ``total_carbon`` closes to machine precision (D-71).
_ETHANOL_PER_ACETALDEHYDE = M_ETHANOL / M_ACETALDEHYDE

#: O₂ : SO₂ stoichiometry of sulfite oxidation — **2 mol SO₂ oxidised per mol O₂ consumed**
#: (decision D-72). The Danilewicz coupled-oxidation mechanism spends a bisulfite at *two* steps
#: per O₂ reduced: one HSO₃⁻ reduces the o-quinone back to the o-diphenol (regenerating the
#: catalyst) and one scavenges the resulting H₂O₂ (HSO₃⁻ + H₂O₂ → HSO₄⁻ + H₂O), so 2 SO₂ : 1 O₂.
#: This is also exactly the classic winemaking "~4 mg SO₂ consumed per mg O₂" mass rule of thumb
#: (2·M_SO2/M_O2 = 2·64/32 = 4). A code-with-citation constant like the chemistry carbon counts
#: (Danilewicz 2003/2007 oxygen-reduction mechanism; Boulton et al. 1996), NOT an uncertain YAML
#: parameter — it is reaction stoichiometry, not a rate. Distinct from the D-47 acetaldehyde–SO₂
#: *binding* (which reversibly sequesters ``so2_total`` without removing it): this route
#: *oxidises* SO₂ to sulfate and permanently removes it, so the two do not double-count.
_SO2_PER_O2 = 2.0  # mol SO₂ oxidised per mol O₂ consumed via the sulfite-scavenging route

#: The two Strecker-aldehyde product species (decision D-75), from the one chemistry source.
#: Each is a single-molecule pool (not a lump): methional is the methionine-derived "cooked-potato"
#: oxidative off-note (C4H8OS), phenylacetaldehyde the phenylalanine-derived "honey" note (C8H8O).
#: Naming them here keeps the carbon draw and the ``total_carbon`` weighting on one species (D-19).
_METHIONAL_SPECIES = "methional"
_PHENYLACETALDEHYDE_SPECIES = "phenylacetaldehyde"

#: The **oxidative** (quinone-driven) Strecker route's two products (decision D-75), speciated at
#: D-100: ``(pool, mol-fraction parameter or None for the remainder, precursor amino acid)``.
#: ``f_methional`` splits the aldehyde flux between them, and each half now draws **its own**
#: amino acid rather than the lumped arginine pool — methional really is methionine minus its
#: carboxyl, phenylacetaldehyde really is phenylalanine minus its. The thermal D-87 route makes
#: the same two molecules from the same two precursors (:data:`_MAILLARD_PRODUCTS`), so the two
#: routes correctly compete for methionine and phenylalanine while ``ProcessSet`` sums them.
_STRECKER_ROUTES: tuple[tuple[str, str | None, str], ...] = (
    (_METHIONAL_SPECIES, "f_methional", "methionine"),
    (_PHENYLACETALDEHYDE_SPECIES, None, "phenylalanine"),
)

#: The Strecker decarboxylation releases exactly **1 mol CO₂ per mol aldehyde** — the amino acid's
#: carboxyl carbon. On the carbon ledger (unlike ``o2``), so it is a genuine product term the carbon
#: bookkeeping must route, not an off-ledger emission (D-75).
_CO2_PER_STRECKER_ALDEHYDE = 1.0

#: The **five** products of the **non-oxidative THERMAL** Strecker route (decision D-87, the
#: sweet-wine / Madeira suite :class:`MaillardStrecker` produces) and their per-product config:
#: ``(pool name, molar mass, composition-weight parameter, precursor)``. The pool name is also the
#: :mod:`chemistry` species key (so ``carbon_mass_fraction`` weights the carbon draw +
#: ``total_carbon`` on one species, D-19). Each precursor's draw is sized to ITS OWN product carbon
#: (§ derivatives, D-100), and each releases **1 mol CO₂** (the amino acid's own carboxyl carbon).
#: ``total_carbon`` closes for **any** CO₂ attribution — a mis-keyed CO₂ term would pass every
#: conservation test silently — so the produced µg/L levels are anchored to literature (the D-75
#: follow-up fidelity lesson) and the 1:1 stoichiometry is pinned by a *driven* test (D-105/D-106).
#: Two are **shared** with the D-75 oxidative route (same molecules, same pools/thresholds); three
#: are D-87-only.
#:
#: **``precursor`` became load-bearing at D-100**: each product's carbon is drawn from the amino
#: acid it is actually made from, not from the lumped arginine pool. The Strecker map is the same
#: one the Ehrlich fusels read (leucine→3-methylbutanal as leucine→isoamyl alcohol), which is
#: exactly why the two routes compete for real.
#:
#: **THIS TABLE HAD SIX ROWS AND TWO FLAGS UNTIL D-107, AND SOTOLON WAS BOTH OF THEM.** ``sotolon``
#: carried ``decarboxylates=False`` (it is not a decarboxylation product) *and* ``de_novo=True``
#: (its
#: carbon is not its precursor's), and D-104 wrote that "the two flags travel together for the same
#: underlying reason: sotolon is not made the way the other five are". That reason was right and the
#: conclusion it licensed was too small: sotolon is not a Strecker degradation **at all** — it is an
#: aldol of α-ketobutyrate + acetaldehyde (Pham *et al.* 1995), so it did not belong in this
#: Process,
#: and the flags were the model reporting that fact in the only vocabulary it had. D-107 moved it to
#: :class:`SotolonAldolCondensation`; **both flags then had exactly one value across all five
#: remaining rows and were deleted**. What is left is five true Strecker aldehydes that all
#: decarboxylate and all take every gram from their own amino acid — a table with no exceptions in
#: it, which is the shape the physics always had.
#:
#: **Sotolon was the canary three times, always through threonine** — D-99's honest fusel rise
#: killed
#: it, D-100 revived it by speciating the pool, D-104's anabolic sink killed it again by (correctly)
#: finishing threonine off, and D-104's ``de_novo`` rescued it by ungating the rate. Every one of
#: those was a fight over a threonine dependence **sotolon never had**: threonine is its
#: *grand*parent (threonine → α-ketobutyrate → sotolon), and the model kept re-deriving the symptom
#: because the intermediate was missing. It is not missing now.
_MAILLARD_PRODUCTS: tuple[tuple[str, float, str, str], ...] = (
    ("methional", M_METHIONAL, "w_maillard_methional", "methionine"),
    (
        "phenylacetaldehyde",
        M_PHENYLACETALDEHYDE,
        "w_maillard_phenylacetaldehyde",
        "phenylalanine",
    ),
    ("2_methylbutanal", M_2_METHYLBUTANAL, "w_maillard_2_methylbutanal", "isoleucine"),
    ("3_methylbutanal", M_3_METHYLBUTANAL, "w_maillard_3_methylbutanal", "leucine"),
    ("2_methylpropanal", M_2_METHYLPROPANAL, "w_maillard_2_methylpropanal", "valine"),
)

#: The caramelization carbon-park species (decision D-88): the sugar carbon :class:`Caramelization`
#: consumes lands here (on ``total_carbon``, booked at melanoidin's caramelan-stand-in fraction), so
#: the sugar → melanoidin transfer closes exactly. Named here so the carbon draw + ``total_carbon``
#: weighting ride on one species (D-19), the ``_ESTER_SPECIES`` discipline.
_MELANOIDIN_SPECIES = "melanoidin"

#: The N-bearing Maillard melanoidin carbon+nitrogen-park species (decision D-89): the sugar carbon
#: AND the amino-acid carbon+nitrogen :class:`MaillardBrowning` consumes land here (on
#: ``total_carbon`` *and* ``total_nitrogen``, booked at the glucose–glycine stand-in fraction), so
#: the sugar + amino_acids → maillard_melanoidin transfer closes both ledgers exactly. Named here so
#: the carbon/nitrogen draws + the ``total_carbon``/``total_nitrogen`` weighting ride on one species
#: (D-19), the ``_MELANOIDIN_SPECIES`` discipline.
_MAILLARD_MELANOIDIN_SPECIES = "maillard_melanoidin"

#: The oak extractives and their set-and-hold ceiling slots (decisions D-77/D-78). Each extracted
#: pool (the first element) rises toward its own saturation ceiling (the second element); the
#: ceiling slots are **constant state** written *only* by the ``add_oak`` verb (``oak_gpl`` ×
#: toast-specific yield) and read — never written — here, the ``cation_charge`` set-and-hold idiom.
#: FIVE AROMA extractives — whiskey lactone (coconut, light-toast dominant), vanillin
#: (vanilla, medium-toast peak), guaiacol (smoky, heavy-toast), eugenol (clove, heavy-toast) (D-77)
#: and ``furaneol`` (caramel/toffee, rising with toast — the caramel furanone, D-94) — plus the
#: ``ellagitannin`` TASTE extractive (D-78, light-toast dominant / declining with toast —
#: thermolabile). Their *extraction* is identical diffusion-to-a-ceiling; ellagitannin additionally
#: feeds the O₂ sub-axis via the separate :class:`EllagitanninOxidation` sink (the aroma five draw
#: no O₂). ``furaneol`` is an ordinary aroma extractive here (off every ledger, wood/spirit-derived)
#: — it does NOT touch the on-ledger D-88 caramelization ``melanoidin``, so no collision (D-94).
_OAK_COMPOUND_CEILINGS: tuple[tuple[str, str], ...] = (
    ("whiskey_lactone", "whiskey_lactone_ceiling"),
    ("vanillin", "vanillin_ceiling"),
    ("guaiacol", "guaiacol_ceiling"),
    ("eugenol", "eugenol_ceiling"),
    ("furaneol", "furaneol_ceiling"),
    ("ellagitannin", "ellagitannin_ceiling"),
)

#: The two on-ledger species of the acetaldehyde-bridged condensation carbon transfer (decision
#: D-80). The acetaldehyde consumed is debited at its own carbon fraction and re-deposited into the
#: ``ethyl_bridge`` pool at the ethylidene fraction (the :class:`EsterHydrolysis` carbon-exact
#: split):
#: acetaldehyde (C2H4O) loses its carbonyl O as water on bridging, leaving the two-carbon ethylidene
#: (C2H4). Using these two fractions to release and re-deposit is what makes ``total_carbon`` close
#: to
#: machine precision through the transfer — the SPLIT-LEDGER capture that keeps the on-ledger
#: acetaldehyde carbon from vanishing into the off-ledger grape-phenolic pigment.
_BRIDGE_ACETALDEHYDE_SPECIES = "acetaldehyde"
_ETHYL_BRIDGE_SPECIES = "ethylidene"


class EsterHydrolysis(Process):
    """Aging hydrolysis of the fruity banana ester toward equilibrium (decisions D-69/D-96).

    ``d(isoamyl_acetate)/dt = -k_ester_hydrolysis · f(T) · max(0, isoamyl_acetate -
    isoamyl_acetate_eq)`` — first-order net decay of the ``isoamyl_acetate`` pool toward the lower
    equilibrium floor ``isoamyl_acetate_eq`` (not to zero), with ``f(T) = arrhenius_factor(T,
    E_a_ester_hydrolysis, T_ref)`` the sourced warmer-ages-faster factor. The released ester
    carbon (``rate·c(isoamyl_acetate)``) is split **5:2** into ``fusels`` (isoamyl alcohol,
    the alcohol product) and ``Byp`` (succinic-stand-in acetic acid, the acid product) — so
    aging fades the banana OAV, raises the fusel OAV, and drifts VA/pH up. See the module
    docstring for the full carbon algebra and the §4.3 firewall it inherits.

    **The 5:2 split is now EXACT — D-96 retired D-69's documented stand-in seam.** D-69 had
    to debit the *lumped* pool at its ledger-fixed **ethyl acetate** (C4) weighting while
    splitting the released carbon 5:2 as though the molecule were **isoamyl acetate** — the
    coherent stand-in reaction (isoamyl acetate → isoamyl alcohol C5 + acetic acid C2), chosen
    because this is a *sensory* Process and D-67 had committed the pool's OAV to isoamyl
    acetate. Debited molecule ≠ split molecule: a mismatch D-69 inherited from D-19/D-67,
    documented and invisible to every conservation test (closure holds for any split summing
    to 1). D-96 split the lump into single-molecule pools, so this Process now debits
    ``isoamyl_acetate`` **itself**: the molecule being hydrolysed, the molecule the 5:2 came
    from, and the molecule the OAV reads are finally all the same one. The ratio is unchanged
    — it was always isoamyl acetate's — but it is no longer a stand-in.

    **Scope — the banana ester only (D-96).** Of the three esters in ``ESTER_SPECS`` only
    ``isoamyl_acetate`` hydrolyses here, because the fruity→fusel fade *is* this Process's
    reason to exist (D-68) and that fade lives entirely in this pathway. Ethyl acetate
    hydrolysis is **deferred**: it is sensorially mute (its products, ethanol and acetic acid,
    are already bulk-present), so leaving that pool stable through aging is a documented v1
    simplification, not an oversight. Were it added, ethanol's honest destination is the core
    ``E`` slot — safe, because ``total_mass`` is scoped to ``{S, E, CO2}`` and asserted only on
    a byproduct-free configuration, where the ester pools are identically zero and this Process
    is inert. Ethyl hexanoate hydrolysis is likewise deferred.

    Off during the ferment (no fermentative-flux gate; it is temperature- and pool-driven);
    enabled only in a post-fermentation aging segment (D-68/D-70). Tier **speculative**.
    """

    name = "ester_hydrolysis"
    tier = Tier.SPECULATIVE
    #: Decays the ``isoamyl_acetate`` pool and routes the released carbon to the alcohol product
    #: (``isoamyl_alcohol``) and the acid product (``Byp``) — an on-ledger inter-pool transfer,
    #: so it touches those three and nothing else (no ``S``/``E``/``CO2``; aging draws no sugar).
    #: Since D-99 the alcohol lands in the isoamyl pool SPECIFICALLY: hydrolysing isoamyl
    #: acetate yields 3-methylbutan-1-ol and nothing else, so depositing into the old lump
    #: silently credited a share of it to four other molecules — and crediting the C5 ISOMER
    #: ``active_amyl_alcohol`` would be a different compound with a ~5.5× different potency.
    #: Since D-115 the two label tracers join them — the transfer carries its label as well as
    #: its carbon, or an aging segment would dilute the alcohol pool's enrichment with returned
    #: molecules booked as unlabelled.
    touches = (
        HYDROLYSING_ESTER.pool,
        ISOAMYL_ALCOHOL.pool,
        "Byp",
        *(tracer.tracer_pool for tracer in VALINE_LABEL_TRACERS),
    )
    #: ``k_ester_hydrolysis``/``E_a_ester_hydrolysis``/``isoamyl_acetate_eq`` are this Process's own
    #: (aging.yaml, D-69; ``k``/``E_a`` re-anchored to Ramey & Ough 1980's real-wine measurement at
    #: D-123, ``isoamyl_acetate_eq`` still author-estimate); ``T_ref`` is shared with every other
    #: Arrhenius rate. Their tiers cap the ``isoamyl_acetate``/``isoamyl_alcohol``/``Byp`` output
    #: tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_ester_hydrolysis",
        "E_a_ester_hydrolysis",
        "isoamyl_acetate_eq",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        ester = float(y[schema.slice(HYDROLYSING_ESTER.pool)][0])
        # Net decay toward the equilibrium floor: the excess above the floor, never below zero
        # (below the floor there is no net hydrolysis — the reverse formation is deferred, D-68).
        # max(0, ...) with isoamyl_acetate_eq > 0 also absorbs a solver undershoot (ester < 0 ⇒ 0).
        excess = max(0.0, ester - params["isoamyl_acetate_eq"])
        if excess <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_ester_hydrolysis"], params["T_ref"])
        rate = params["k_ester_hydrolysis"] * f_t * excess  # g isoamyl acetate/L/h decayed

        # The released carbon is now the REAL molecule's (isoamyl acetate, C7) — D-96 retired
        # D-69's debit-as-ethyl-acetate / split-as-isoamyl-acetate stand-in. Split it 5:2 and
        # re-deposit through each product pool's own carbon fraction, so total_carbon closes to
        # machine precision for any split summing to 1 — and 5:2 is now the true stoichiometry
        # of the molecule actually being debited (C7 → isoamyl alcohol C5 + acetic acid C2).
        carbon_released = rate * carbon_mass_fraction(HYDROLYSING_ESTER.species)  # g C/L/h
        d[schema.slice(HYDROLYSING_ESTER.pool)] = -rate
        d[schema.slice(ISOAMYL_ALCOHOL.pool)] = (
            _FUSEL_CARBON_SHARE * carbon_released / carbon_mass_fraction(ISOAMYL_ALCOHOL.species)
        )
        d[schema.slice("Byp")] = (
            _BYP_CARBON_SHARE * carbon_released / carbon_mass_fraction(_BYP_SPECIES)
        )

        # D-115: the label comes back with the C5. Hydrolysis is the exact reverse of the
        # acetylation, so an ester molecule that was valine-derived returns a valine-derived
        # alcohol molecule — mole for mole, hence the same 5:2 carbon algebra applied to the
        # tracer. Debited at the ESTER's own fraction (not the alcohol's): using the alcohol's
        # would assume the answer this route exists to compute, which is the vacuity trap the
        # second tracer slot was added to avoid.
        f_ester = labelled_fraction(y, schema, _ESTER_TRACER)
        if f_ester > 0.0:
            d[schema.slice(_ESTER_TRACER.tracer_pool)] = -rate * f_ester
            d[schema.slice(_ALCOHOL_TRACER.tracer_pool)] = (
                _FUSEL_CARBON_SHARE
                * carbon_released
                * f_ester
                / carbon_mass_fraction(ISOAMYL_ALCOHOL.species)
            )
        return d


class OxidativeAcetaldehyde(Process):
    """Oxidative aging: dissolved O₂ drives ethanol → acetaldehyde (decision D-71).

    The first **oxidative** aging Process and the head of the O₂ sub-axis. As a finished wine/beer
    takes up oxygen (bottle ingress, micro-oxygenation, barrel), the classic oxidative-aging
    reaction accumulates **acetaldehyde** — the "sherry"/bruised-apple/nutty oxidised note — which
    the D-67 OAV lens *already* reads (it is the same ``acetaldehyde`` pool the D-27 buffer fills,
    read as "green apple" fresh, "oxidised" when it climbs). So — like :class:`EsterHydrolysis`,
    and per the D-68 selection criterion — this Process moves an OAV the lens already reads and
    needs **no new aroma pool**; the one new slot is the ``o2`` *substrate*.

    ``d(o2)/dt = −r_O2`` with ``r_O2 = k_ethanol_oxidation · f(T) · [O2]`` — first-order in the
    dissolved-O₂ pool, where ``k_ethanol_oxidation`` is this route's **share** of the total
    O₂-depletion rate (D-73; *not* the whole rate — see the ledger paragraph below and
    :class:`SulfiteOxidation`), and ``f(T) = arrhenius_factor(T, E_a_ethanol_oxidation, T_ref)`` the
    sourced warmer-oxidises-faster factor. The oxidised carbon lands as acetaldehyde at a molar
    **yield**::

        d(acetaldehyde)/dt = +y_acetaldehyde_per_o2 · (r_O2 / M_O2) · M_acetaldehyde
        d(E)/dt            = −d(acetaldehyde)/dt · M_ethanol / M_acetaldehyde

    **O₂, not ethanol, is the rate-limiting reactant** (the D-71 design crux). Ethanol sits at
    ~100 g/L, essentially constant across aging, so a rate first-order in *ethanol* would be a
    constant rate in disguise — acetaldehyde rising linearly and **unbounded**, pinning the kinetic
    limit on the wrong species. Making the rate first-order in the finite ``o2`` pool instead gives
    the correct **saturating** behaviour: as the O₂ charge is consumed the pool decays toward zero
    and acetaldehyde plateaus, the bottle-aging reality. (Mechanistically the real path is *coupled*
    oxidation — O₂ oxidises o-diphenols → quinones + H₂O₂, then H₂O₂ oxidises ethanol → acetaldehyde
    (Wildenradt & Singleton 1974); the phenolic catalyst is folded into ``k_ethanol_oxidation`` in
    v1, a documented lump since no general phenol pool is tracked.)

    **Carbon — the clean reverse of the D-27 reduction.** Ethanol and acetaldehyde are both C2, so
    the ``E → acetaldehyde`` transfer is mole-for-mole; weighting the ethanol debit by
    ``M_ethanol/M_acetaldehyde`` (:data:`_ETHANOL_PER_ACETALDEHYDE`) makes the carbon that leaves
    ``E`` land exactly in ``acetaldehyde``, so ``total_carbon`` closes to machine precision. This is
    the mirror of :class:`~fermentation.core.kinetics.acetaldehyde.AcetaldehydeReduction`
    (acetaldehyde → ethanol), and during aging that reduction is **inert** — it is gated on *viable*
    ``X``, which is 0 in a racked/finished wine — so oxidation does not fight it: the acetaldehyde
    accumulates, correctly (a live-yeast ferment would instead reduce it straight back).

    **The O₂ pool is off every ledger; this Process draws its ethanol-oxidation *share* of it
    (D-73).** ``o2`` is carbon-free (``total_carbon``/``total_mass``/``total_nitrogen`` weight only
    their named pools, so ``o2`` contributes 0 to each, like ``h2s``/``iso_alpha``). The Process
    consumes ``r_O2 = k_ethanol_oxidation · f(T) · [O2]`` — its **own share** of the shared O₂
    budget, not the whole flux — and every mol it consumes yields ``y_acetaldehyde_per_o2`` mol
    acetaldehyde, the *true* per-O₂ stoichiometric yield of the ethanol route. Sibling sinks
    (:class:`SulfiteOxidation` and the dominant always-on :class:`PhenolicBrowning`, D-74; Strecker
    to
    come) each draw their **own** share via their own rate constant, and ``ProcessSet`` summing
    splits
    the finite O₂ among them by
    ``kᵢ / Σk`` — so the pool depletes *once*, and the competition that suppresses acetaldehyde
    lives in the rate constants, not in a shaded-down yield. Because O₂ carries no carbon,
    "spending" it is not a conservation violation; the carbon that *does* move (into acetaldehyde)
    is borrowed carbon-exactly from ``E``.

    *Ethanol oxidation is now the secondary always-on sink (D-74).* ``k_ethanol_oxidation`` is the
    ethanol-oxidation *share*, and as of **D-74** it is no longer the whole always-on flux:
    :class:`PhenolicBrowning` is a co-resident always-on O₂ sink, and the *dominant* one (phenol
    autoxidation is the primary O₂ consumer; ethanol oxidation is the secondary H₂O₂ fate). So
    ``k_ethanol_oxidation`` was **reduced 5.0e-4 → 2.0e-4** (browning takes the larger 3.0e-4 share)
    with ``k_ethanol + k_browning = 5.0e-4`` holding the empirical total O₂-depletion timescale —
    the
    anchor — unchanged. Aged acetaldehyde is therefore **lower** than the D-73 sole-sink estimate:
    with no SO₂ the ethanol route's share of a fully-consumed O₂ charge is ``k_ethanol / (k_ethanol
    +
    k_browning) = 0.4``, so the cumulative acetaldehyde is ~40 % of the sole-sink value (the
    "partitions
    down" D-73 promised, now realised). The D-72 substrate-gated :class:`SulfiteOxidation` needed no
    such re-baseline — it is simply zero without SO₂ — so the *always-on* re-baseline D-71 flagged
    and
    D-72 deferred is exactly what D-73 enabled and D-74 spent.

    *Mass carries the standing gap, scoped out by construction.* ``total_mass`` weights only
    ``{S, E, CO2}`` (the byproduct-free validated-core check, D-8/D-16): this Process debits ``E``
    into the unweighted ``acetaldehyde`` pool, so ``total_mass`` would drift *if asserted here* —
    but it is the same E↔acetaldehyde seam the D-27 buffer already carries, and ``total_mass`` is
    never asserted on an aging run (the aroma pools are all active). Carbon is the invariant.

    **Isolable + O₂-gated (prime directive #3).** Wired into both media's ``_AGING_PROCESSES`` tuple
    but **disabled at the compile seam** (aging is post-ferment); ``begin_aging`` enables it
    alongside :class:`EsterHydrolysis`. With no O₂ dosed the pool is 0 and the ``o2 <= 0`` guard is
    exact — the contribution is byte-for-byte zero, so a ``begin_aging`` run without ``add_oxygen``
    is purely *reductive* aging (the EsterHydrolysis-only case). Tier **speculative** (the aging
    axis is the Tier-3 frontier; the oxidation *form* is sourced, the magnitude an estimate).
    """

    name = "oxidative_acetaldehyde"
    tier = Tier.SPECULATIVE
    #: Consumes the dissolved-O₂ substrate and books the oxidised carbon as ``acetaldehyde``,
    #: borrowed carbon-exactly from ``E`` (the D-27 reduction reversed). Touches those three and
    #: nothing else — ``o2`` is off every ledger, so only the ``E → acetaldehyde`` transfer is on
    #: the carbon books, and it closes exactly.
    touches = ("o2", "acetaldehyde", "E")
    #: ``k_ethanol_oxidation``/``E_a_ethanol_oxidation``/``y_acetaldehyde_per_o2`` are this
    #: Process's own (aging.yaml, D-71); ``T_ref`` is shared with every Arrhenius rate. Tiers cap
    #: ``o2``/``acetaldehyde``/``E`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_ethanol_oxidation",
        "E_a_ethanol_oxidation",
        "y_acetaldehyde_per_o2",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        o2 = float(y[schema.slice("o2")][0])
        # No oxidant ⇒ no oxidative acetaldehyde: reductive aging (screwcap/inert) and the exact
        # isolability guard (an un-dosed begin_aging run is byte-for-byte the ester-only aging).
        # ``<= 0`` also absorbs a solver undershoot (o2 < 0 ⇒ no spurious production).
        if o2 <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_ethanol_oxidation"], params["T_ref"])
        # This route's SHARE of the O₂-depletion rate (D-73), not the whole flux: sibling sinks draw
        # their own shares and ProcessSet sums, so O₂ splits by kᵢ/Σk and the pool depletes once.
        r_o2 = params["k_ethanol_oxidation"] * f_t * o2  # g O2/L/h consumed by the ethanol route
        # Every mol O₂ this route consumes yields y_acetaldehyde_per_o2 mol acetaldehyde — the TRUE
        # per-O₂ yield of the route (D-73), not shaded for competitors. moles O2 = r_o2/M_O2.
        acet_rate = params["y_acetaldehyde_per_o2"] * (r_o2 / M_O2) * M_ACETALDEHYDE  # g/L/h
        d[schema.slice("o2")] = -r_o2  # this route's O₂ share is consumed (off every ledger)
        d[schema.slice("acetaldehyde")] = acet_rate
        # Carbon-exact C2 borrow from ethanol (the D-27 reduction reversed). No clamp needed: during
        # aging E ~ 100 g/L and acet_rate is trace, so this never drives E negative.
        d[schema.slice("E")] = -acet_rate * _ETHANOL_PER_ACETALDEHYDE
        return d


class SulfiteOxidation(Process):
    """Oxidative aging: dissolved O₂ oxidises free bisulfite → sulfate, spending SO₂ (D-72).

    The second **oxidative** aging Process and the first *sibling* on the O₂ sub-axis opened by
    :class:`OxidativeAcetaldehyde` (D-71). SO₂ is wine's antioxidant precisely because bisulfite is
    a *faster* O₂ scavenger than ethanol: as O₂ is taken up it is preferentially spent oxidising the
    free bisulfite pool (HSO₃⁻ → sulfate) rather than ethanol, so **while free SO₂ lasts, O₂ is
    diverted and little oxidative acetaldehyde forms; once SO₂ is exhausted, acetaldehyde climbs** —
    the celebrated wine-chemistry threshold, and the payoff of putting both sinks on one shared
    ``o2`` budget. There is nothing to build for the diversion itself: this Process and
    :class:`OxidativeAcetaldehyde` both draw down the same ``o2`` pool, so ``ProcessSet`` summing
    makes the O₂ split between them by their rates — the fraction reaching acetaldehyde is
    ``k_eth / (k_eth + k_so2·[HSO₃⁻])``, small while SO₂ is present, → 1 once it is gone. **No new
    pool** is needed (the D-68 selection criterion): it decrements the existing ``so2_total`` slot.

    ``d(o2)/dt = −r`` with ``r = k_so2_oxidation · f(T) · [O2] · [HSO₃⁻]`` (**bilinear** in the
    dissolved-O₂ pool and the free-bisulfite driver, ``f(T) = arrhenius_factor(T, E_a_so2_oxidation,
    T_ref)`` the sourced warmer-oxidises-faster factor), and::

        d(so2_total)/dt = −_SO2_PER_O2 · (r / M_O2) · M_SO2

    consumes **2 mol SO₂ per mol O₂** (:data:`_SO2_PER_O2` — the Danilewicz mechanism: one bisulfite
    reduces the o-quinone, one scavenges the H₂O₂; = the classic ~4 mg-SO₂-per-mg-O₂ mass rule).

    **The driver is free BISULFITE, not molecular SO₂ (the D-72 crux).** Molecular SO₂·H₂O is the
    reactive *antimicrobial* form (D-22); the HSO₃⁻ anion is the reactive *antioxidant* nucleophile
    (the reducer of quinones and scavenger of H₂O₂ — Danilewicz; this codebase's own
    :func:`~fermentation.core.acidbase.bisulfite_fraction` already names HSO₃⁻ "the reactive
    nucleophile"). So the rate reads :func:`~fermentation.core.acidbase.bisulfite_so2_at_ph` =
    ``free_SO₂ · bisulfite_fraction(pH)``, using only **free** SO₂ (bound bisulfite is already
    spent). Because bisulfite is ~0.94–0.99 of free across wine pH, the pH-coupling is mild — but a
    *stronger* coupling enters through **free** SO₂: as the sibling :class:`OxidativeAcetaldehyde`
    produces acetaldehyde, that acetaldehyde binds SO₂ (D-47), free SO₂ falls, and this scavenging
    rate **self-throttles**. Oxidation thus erodes SO₂'s protective capacity two ways — oxidative
    removal here + binding via D-47 — an emergent feedback the bilinear form buys over a plain
    SO₂-presence gate.

    **Off every ledger, no conservation term.** Both ``o2`` (D-71) and ``so2_total`` (a dosed,
    carbon-free input — D-22/D-28) are off ``total_carbon``/``total_mass``/``total_nitrogen`` (there
    is no sulfur ledger), so oxidising SO₂ to untracked sulfate moves no conserved quantity — this
    Process touches only those two slots and asserts nothing. Distinct from the D-47 binding, which
    *repartitions* ``so2_total`` (free ⇄ bound, reversible) without removing it: oxidation
    *removes* it (→ sulfate), so booking 2:1 here does not double-count the binding readout.

    **Wine-only + isolable + doubly substrate-gated (prime directive #3).** ``so2_total`` and the
    acid/cation pH slots are wine-only (beer's pH system is deferred, D-18), so — like the MLF and
    Brett Processes — this is wired into the *wine* medium only; the ``SO2_STATE_KEY not in schema``
    guard makes it a hard no-op on beer besides. Wired **disabled at the compile seam** (aging is
    post-ferment); ``begin_aging`` enables it alongside :class:`EsterHydrolysis` /
    :class:`OxidativeAcetaldehyde` (:data:`~fermentation.scenario.compile._AGING_GATED_PROCESSES`).
    With no O₂ *or* no SO₂ dosed the ``o2 ≤ 0`` / ``so2_total ≤ 0`` guards return byte-for-byte zero
    (and skip the pH solve), so a reductive aging (no ``add_oxygen``) or an unsulfited aging is
    exactly the case without this Process. Tier **speculative** (the aging axis is the Tier-3
    frontier; the *form* — O₂-limited, bisulfite-driven, warmer-faster, 2:1 stoichiometry — is
    sourced, the rate *magnitude* an order-of-magnitude estimate).
    """

    name = "sulfite_oxidation"
    tier = Tier.SPECULATIVE
    #: Consumes the dissolved-O₂ substrate and oxidises the free-bisulfite share of ``so2_total`` to
    #: (untracked) sulfate — both slots off every ledger, so nothing conserved moves; it touches
    #: those two and nothing else.
    touches = ("o2", "so2_total")
    #: ``k_so2_oxidation``/``E_a_so2_oxidation`` are its own (aging.yaml, D-72); ``T_ref`` is shared
    #: with every Arrhenius rate. The pKa/binding params read through ``acidbase`` (to
    #: derive free bisulfite at the solved pH) are omitted — all plausible, and the Process is
    #: already speculative, so they add no tier headline (the MalolacticConversion/brett rule).
    reads: tuple[str, ...] = ("k_so2_oxidation", "E_a_so2_oxidation", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (beer's pH/SO₂ system is deferred, D-18): a hard no-op on any schema
        # without them, belt-and-suspenders to the wine-only wiring.
        if SO2_STATE_KEY not in schema or "o2" not in schema:
            return d
        o2 = float(y[schema.slice("o2")][0])
        so2_total = float(y[schema.slice(SO2_STATE_KEY)][0])
        # No oxidant OR no SO₂ ⇒ no scavenging: reductive/unsulfited aging is byte-for-byte the
        # case without this Process, and neither guard pays a per-RHS pH solve for a zero result.
        # ``<= 0`` also absorbs solver undershoot (o2 < 0 / so2_total < 0 ⇒ no spurious use).
        if o2 <= 0.0 or so2_total <= 0.0:
            return d
        ph = ph_of_state(y, schema, params)
        # The reactive ANTIOXIDANT species is free bisulfite HSO₃⁻ (not molecular SO₂ — that is the
        # antimicrobial form); it self-throttles as acetaldehyde binds SO₂ (D-47) and free falls.
        bisulfite = bisulfite_so2_at_ph(y, schema, params, ph)  # g/L as SO₂
        if bisulfite <= 0.0:  # all free SO₂ sequestered by carbonyls ⇒ nothing reactive left
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_so2_oxidation"], params["T_ref"])
        r_o2 = params["k_so2_oxidation"] * f_t * o2 * bisulfite  # g O2/L/h via the SO₂ route
        d[schema.slice("o2")] = -r_o2
        # 2 mol SO₂ oxidised per mol O₂ (Danilewicz coupled oxidation): moles O₂ = r_o2 / M_O2.
        d[schema.slice(SO2_STATE_KEY)] = -_SO2_PER_O2 * (r_o2 / M_O2) * M_SO2
        return d


class PhenolicBrowning(Process):
    """Oxidative aging: dissolved O₂ oxidises phenolics → brown pigment (decision D-74).

    The second **oxidative** aging Process and the first **always-on** sink on the O₂ sub-axis (SO₂
    oxidation is *substrate-gated* — zero without dosed SO₂; browning runs whenever O₂ is present).
    As a finished wine/beer takes up oxygen, dissolved O₂ oxidises **phenolics** (o-diphenols →
    o-quinones, polymerising to brown melanoidin/quinone pigment) — the gold→amber→brown of an aged
    / maderised white wine, and the oxidative darkening of stale beer. This is the **dominant** O₂
    consumer in wine oxidation: phenol autoxidation is the *primary* O₂ sink, while ethanol
    oxidation
    to acetaldehyde (:class:`OxidativeAcetaldehyde`) is a *secondary* H₂O₂ fate — so ``k_browning``
    is
    the **larger** share of the shared O₂ budget, and browning **diverts most of the always-on O₂
    away from ethanol oxidation, suppressing oxidative acetaldehyde**. That suppression is the
    always-on analogue of SO₂'s protection (D-72): SO₂ suppresses *until it is spent*; browning, a
    co-resident always-on sink, suppresses *permanently* (the acetaldehyde partition ``k_ethanol /
    Σk`` emerges from ``ProcessSet`` summing, for free). Landing browning is what the D-73 rework
    enabled — ``k_ethanol_oxidation`` was reduced 5.0e-4 → 2.0e-4 so ``k_ethanol + k_browning``
    still
    holds the calibrated total O₂-depletion timescale (the anchor); under D-71's "total rate"
    framing
    this always-on sink could not have been added without double-counting.

    ``d(o2)/dt = −r_O2`` with ``r_O2 = k_browning · f(T) · [O2]`` — first-order in the dissolved-O₂
    pool (its **own share**, like :class:`OxidativeAcetaldehyde`), ``f(T) = arrhenius_factor(T,
    E_a_browning, T_ref)`` the sourced warmer-browns-faster factor. The O₂ it consumes accumulates
    the
    browning index::

        d(A420)/dt = +y_a420_per_o2 · (r_O2 / M_O2)

    **The new observable is an optical index, not an aroma or a mass (the D-74 crux).** Unlike the
    other three aging Processes, browning moves **no** pool the D-67 OAV lens already reads — its
    product is brown *pigment*, seen not smelled — so (the D-68 "reuse an existing pool" criterion
    inverted) it needs one new observable: ``A420``, the **absorbance at 420 nm** (dimensionless AU,
    1 cm path — the standard analytical browning measure), **not** a pigment concentration. Two
    consequences follow. (1) ``A420`` is a **state slot** (the ``iso_alpha`` off-ledger-pool
    pattern),
    **not** a post-hoc readout like the D-67 OAV series: browning pigment is *cumulative and
    irreversible* and its O₂ flux is *dynamic* (SO₂ competes for O₂, temperature varies the rate),
    so
    it must be **integrated** along the trajectory — it cannot be reconstructed after the run from
    (dosed − remaining) O₂. (2) Because ``A420`` is an **optical index rather than a mass**, the
    conservation question dissolves: the pigment's carbon would come from an *untracked* phenol
    pool,
    but an absorbance carries none, so ``A420`` is off every ledger by construction (like ``o2`` /
    ``iso_alpha``). This Process therefore touches only ``{o2, A420}`` — **both off every ledger** —
    and moves **nothing conserved at all**: it is the *cleanest* aging Process on the books (even
    :class:`OxidativeAcetaldehyde` still borrows carbon E→acetaldehyde). ``d(A420)/dt ≥ 0`` always
    (monotonic accumulation — no clamp; a solver undershoot ``o2 < 0`` is caught by the guard).

    **Medium-agnostic (D-74 supersedes D-73's provisional "wine-only").** D-73's worked drop-in
    tentatively marked browning "wine — o-diphenols are a wine pool", but there is no o-diphenol
    pool
    (the catalyst is lumped into ``k_browning``, as in :class:`OxidativeAcetaldehyde`), and both
    wine
    and beer carry autoxidising polyphenols that consume O₂ and brown oxidatively — so browning is a
    property of the molecules, not the biology (the shared-``aging.yaml`` discipline), and is wired
    into **both** media like :class:`OxidativeAcetaldehyde`. This is also *forced* to be
    consistent: the
    ``k_ethanol_oxidation`` reduction lives in the **shared** ``aging.yaml`` and applies to both
    media,
    so a wine-only browning sink would leave beer's total O₂-depletion rate silently halved below
    the
    anchor — the very in-tree inconsistency the D-73 rework existed to remove. Medium-agnostic
    browning
    keeps beer's O₂ budget whole (beer runs browning too, records its own ``A420``).

    **Isolable + O₂-gated (prime directive #3).** Wired into both media's ``_AGING_PROCESSES`` tuple
    but **disabled at the compile seam** (aging is post-ferment); ``begin_aging`` enables it with
    the other aging Processes. With no O₂ dosed the ``o2 ≤ 0`` guard is exact and the contribution
    is byte-for-byte zero (``A420`` stays 0), so a ``begin_aging`` run without ``add_oxygen`` is
    purely *reductive* aging — unchanged by this Process. Tier **speculative** (the aging axis is
    the Tier-3 frontier; the browning *form* — O₂-limited, warmer-faster — is sourced, the rate and
    per-O₂ absorbance yield are order-of-magnitude estimates).
    """

    name = "phenolic_browning"
    #: Consumes its share of the dissolved-O₂ substrate and books the oxidised phenol as the
    #: ``A420`` browning index — both slots off every ledger, so nothing conserved moves; it touches
    #: those two and nothing else (the cleanest aging Process — not even a carbon borrow).
    tier = Tier.SPECULATIVE
    touches = ("o2", "A420")
    #: ``k_browning``/``E_a_browning``/``y_a420_per_o2`` are this Process's own (aging.yaml, D-74);
    #: ``T_ref`` is shared with every Arrhenius rate. Tiers cap the ``o2``/``A420`` output tiers via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_browning", "E_a_browning", "y_a420_per_o2", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        o2 = float(y[schema.slice("o2")][0])
        # No oxidant ⇒ no browning: reductive/un-oxygenated aging is byte-for-byte the case without
        # this Process (A420 stays 0). ``<= 0`` also absorbs a solver undershoot (o2 < 0 ⇒ no
        # spurious browning), which keeps d(A420)/dt ≥ 0 (A420 monotonic, never reversed).
        if o2 <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_browning"], params["T_ref"])
        # This route's SHARE of the O₂-depletion rate (the larger, dominant share — D-74),
        # first-order
        # in o2 like the ethanol route; ProcessSet sums the sinks so the pool depletes once.
        r_o2 = params["k_browning"] * f_t * o2  # g O2/L/h consumed by the browning route
        d[schema.slice("o2")] = -r_o2  # this route's O₂ share is consumed (off every ledger)
        # Every mol O₂ this route consumes raises the A420 absorbance index by y_a420_per_o2 (AU per
        # mol O₂/L). moles O₂ = r_o2 / M_O2. A420 is an optical index (off every ledger), not a
        # mass —
        # so nothing conserved moves and no carbon is borrowed (unlike the ethanol route).
        d[schema.slice("A420")] = params["y_a420_per_o2"] * (r_o2 / M_O2)
        return d


class StreckerDegradation(Process):
    """Oxidative aging: O₂ (via quinones) degrades amino acids → Strecker aldehydes (D-75).

    The fifth aging Process and the **third** oxidative sibling on the O₂ sub-axis, after
    :class:`OxidativeAcetaldehyde` (D-71) and :class:`PhenolicBrowning` (D-74). As a finished wine
    takes up oxygen, the **o-quinones** that phenol autoxidation produces (the browning cascade)
    oxidatively deaminate and decarboxylate amino acids to **Strecker aldehydes** — **methional**
    (from methionine, the "cooked-potato" *oxidative off-note*, one of the sharpest markers of an
    oxidised/maderised white wine and of stale beer) and **phenylacetaldehyde** (from phenylalanine,
    the "honey/floral" note of aged white and dessert wines). Unlike the other four aging Processes,
    these products move **no** pool the D-67 OAV lens already reads, so this beat adds **two new
    aroma pools** — two, not one lumped, because the two aldehydes have **opposite sensory valence**
    (methional off, phenylacetaldehyde pleasant), the owner's D-75 fork.

    **Doubly substrate-gated — adds on top, NO re-baseline (the D-75 crux).** The rate is gated on
    the dissolved-O₂ pool **and** the amino-acid pool::

        gate = amino_acids / (K_amino_acids + amino_acids)              # smooth availability, [0,1)
        r_O2 = k_strecker · f(T) · [o2] · gate                         # O₂-limited AND aa-gated
        n_ald = y_strecker_per_o2 · (r_O2 / M_O2)                      # mol total aldehyde /L/h

    Because the O₂ draw itself carries the ``gate``, Strecker is **substrate-gated exactly like**
    :class:`SulfiteOxidation` (which is gated on ``o2`` AND SO₂). D-72 established the load-bearing
    rule: a substrate-gated sink **adds on top of the shared O₂ budget without any re-baseline** —
    zero without its substrate, so the default/beer trajectory is byte-for-byte preserved and
    ``k_ethanol_oxidation + k_browning = 5.0e-4`` is **untouched** (``k_strecker`` is a small extra
    wine-only draw that only fires when ``amino_acids`` is present — dosed nutrient, or future lees
    autolysis refill). This **supersedes** the D-71→D-74 forward-guess ("the next *always-on* sink —
    reduce ``k_ethanol_oxidation`` again to its share"), which wrongly assumed a significant,
    medium-agnostic sink. The ``amino_acids`` pool is the true **limiting reagent** (finite amino
    acid ⇒ finite Strecker aldehyde — the accumulation saturates as the pool is drawn down), so the
    aldehyde *level* is threshold-relevant (µg/L–mg/L vs ~0.5–1 µg/L thresholds) across the whole
    speculative parameter band while the O₂ draw stays a minor, in-band perturbation.

    **Carbon + nitrogen close by construction — the D-45 mercaptans idiom + a CO₂ term.** The
    aldehyde carbon is drawn from **its own precursor** — methionine for methional, phenylalanine
    for phenylacetaldehyde (decision D-100) — and that precursor's nitrogen is **deaminated** back
    to the ``N`` pool, exactly as
    :class:`~fermentation.core.kinetics.mercaptans.AutolyticMercaptan` does; the Strecker
    **decarboxylation** adds one product this idiom did not have — **1 mol CO₂ per mol aldehyde**
    (:data:`_CO2_PER_STRECKER_ALDEHYDE`, the acid's carboxyl carbon), on the carbon ledger. Each
    draw is *sized to that route's product carbon* (its aldehyde + its own carboxyl CO₂), so
    ``total_carbon`` closes to machine precision (the :class:`EsterHydrolysis` multi-product split
    idiom); all the precursor nitrogen lands in ``N`` and the products are nitrogen-free, so
    ``total_nitrogen`` closes. Before D-100 both routes drew the ``amino_acids`` arginine lump —
    **exact on the ledger, approximate on provenance**; that stand-in is **retired** here and in
    mercaptans alike. ``o2`` is off every ledger (D-71), so spending it moves
    nothing conserved. ``total_mass`` ({S,E,CO2}) sees the CO₂ term with no matching S/E debit, but
    it is never asserted on an aging run (the aroma pools are active) — the standing
    :class:`OxidativeAcetaldehyde` scope-out.

    **The inherited quinone double-count lump (documented, not fixed).** Mechanistically the O₂ is
    consumed at the phenol-oxidation step (:class:`PhenolicBrowning`'s draw), making the o-quinones
    that then do the Strecker deamination — so a separate ``k_strecker`` ``[o2]`` draw formally
    double-counts that shared quinone step. But :class:`PhenolicBrowning` and
    :class:`OxidativeAcetaldehyde` **already** double-count it against each other (both independent
    ``[o2]`` draws for one coupled cascade) — the additive-share v1 lump accepted at D-73. Strecker
    following suit is *consistent*; a two-stage (O₂ → quinone pool → {pigment, aldehyde,
    acetaldehyde}) rework is deliberately out of scope. **Scope:** this is the *oxidative*
    (quinone-driven) Strecker route only; the non-oxidative Maillard/sugar-dicarbonyl route (sweet
    wines, thermal) is deferred, keeping Strecker honestly on the ``o2`` sub-axis.

    **Wine-only + isolable + doubly O₂/aa-gated (prime directive #3).** ``amino_acids`` and the
    ``N``-deamination read wine-only slots (beer's amino-acid pool is not tracked, D-32), so — like
    :class:`SulfiteOxidation` — this is wired into the *wine* medium only; the
    ``"amino_acids" not in schema`` guard makes it a hard no-op besides. Wired **disabled at the
    compile seam** (aging is post-ferment); ``begin_aging`` enables it with the other aging
    Processes. With no O₂ or amino acids the ``o2 ≤ 0`` / ``aa ≤ 0`` guards return byte-for-byte
    zero, so a reductive (no ``add_oxygen``) or an amino-acid-free aging is exactly the case without
    this Process. **First aging Process to write ``N``** (via the deamination), so an enabled run
    drops structural ``tier_of("N")`` PLAUSIBLE→SPECULATIVE (the D-45 note). Tier **speculative**
    (the Strecker *form* — O₂-linked, amino-acid-driven, warmer-faster, aldehyde = acid − CO₂ —
    is sourced; every magnitude is an order-of-magnitude estimate).
    """

    name = "strecker_degradation"
    tier = Tier.SPECULATIVE
    #: Consumes its aa-gated share of the dissolved-O₂ substrate; books the two Strecker aldehydes
    #: (``methional``/``phenylacetaldehyde``) + the decarboxylation ``CO2``, drawing each one's
    #: carbon from **its own precursor** (methionine/phenylalanine, decision D-100 — the lumped
    #: arginine draw is retired) and deaminating that precursor's nitrogen to ``N``. ``o2`` is off
    #: every ledger; the C/N transfer closes exactly.
    touches = (
        "o2",
        "methional",
        "phenylacetaldehyde",
        "CO2",
        "N",
        *(precursor for _, _, precursor in _STRECKER_ROUTES),
    )
    #: ``k_strecker``/``E_a_strecker``/``y_strecker_per_o2``/``f_methional`` are this Process's own
    #: (aging.yaml, D-75); ``K_amino_acids`` is the *shared* availability half-saturation (the same
    #: constant the mercaptan/reroute gates read), scaled per-precursor by its
    #: ``must_aa_fraction_*`` share (D-100); ``T_ref`` is shared with every Arrhenius rate. Their
    #: tiers cap the output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_strecker",
        "E_a_strecker",
        "y_strecker_per_o2",
        "f_methional",
        "K_amino_acids",
        "T_ref",
        *(SPEC_BY_SPECIES[precursor].fraction_param for _, _, precursor in _STRECKER_ROUTES),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (beer tracks no amino acids, D-32): a hard no-op on any schema without
        # them, belt-and-suspenders to the wine-only wiring.
        if "methionine" not in schema or "o2" not in schema:
            return d
        o2 = float(y[schema.slice("o2")][0])
        # No O₂ ⇒ no Strecker: reductive aging is byte-for-byte the case without this Process.
        # ``<= 0`` also absorbs solver undershoot (o2 < 0 ⇒ no draw).
        if o2 <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_strecker"], params["T_ref"])
        # The un-gated O₂-depletion driver — a SMALL wine-only add-on (NOT in the 5.0e-4 always-on
        # anchor; substrate-gated, adds on top like SulfiteOxidation, D-72).
        driver = params["k_strecker"] * f_t * o2  # g O2/L/h before the availability gates
        f_meth = params["f_methional"]  # mol fraction methional; the rest is phenylacetaldehyde
        shares = {_METHIONAL_SPECIES: f_meth, _PHENYLACETALDEHYDE_SPECIES: 1.0 - f_meth}

        # Each aldehyde is throttled by ITS OWN precursor's relative-depletion gate (decision
        # D-100): methional stops when methionine runs out, phenylacetaldehyde when phenylalanine
        # does — independently, as the chemistry demands. At must-spectrum composition both gates
        # equal the old lumped gate, so this reduces byte-for-byte to the D-75 rate. The gates also
        # throttle the O₂ draw (each route takes only the share its substrate can support), so O₂,
        # carbon and nitrogen still vanish together and no precursor is driven negative.
        r_o2 = 0.0  # g O2/L/h actually consumed, summed over the two gated routes
        co2_mol = 0.0  # mol CO2/L/h — 1 per aldehyde (both routes decarboxylate)
        product_rates: list[tuple[str, float]] = []
        precursor_carbon: dict[str, float] = {}
        for pool, _fparam, precursor in _STRECKER_ROUTES:
            gate_i = depletion_gate(y, schema, params, (SPEC_BY_SPECIES[precursor],))
            if gate_i <= 0.0:
                continue  # this precursor is exhausted ⇒ its aldehyde (and its O₂ share) stop
            r_o2_i = shares[pool] * driver * gate_i  # this route's gated O₂ share
            n_i = params["y_strecker_per_o2"] * (r_o2_i / M_O2)  # mol aldehyde/L/h
            if n_i <= 0.0:
                continue
            rate_i = n_i * (M_METHIONAL if pool == _METHIONAL_SPECIES else M_PHENYLACETALDEHYDE)
            co2_i = _CO2_PER_STRECKER_ALDEHYDE * n_i * M_CO2  # the amino acid's carboxyl carbon
            r_o2 += r_o2_i
            co2_mol += _CO2_PER_STRECKER_ALDEHYDE * n_i
            product_rates.append((pool, rate_i))
            # Carbon THIS precursor must supply: its aldehyde + the CO₂ its own carboxyl released.
            precursor_carbon[precursor] = rate_i * carbon_mass_fraction(pool) + co2_i * (
                carbon_mass_fraction("CO2")
            )
        if not product_rates:
            return d

        # Draw each aldehyde's carbon from ITS OWN precursor, sized to match, and deaminate that
        # precursor's nitrogen to N (the D-45 idiom + CO₂ decarb, speciated at D-100): carbon out
        # of each amino acid == carbon into the product it made, and all their N lands in N
        # (products are N-free), so total_carbon and total_nitrogen close to machine precision.
        nitrogen = sum(
            draw_precursor_carbon(d, schema, precursor, carbon)
            for precursor, carbon in precursor_carbon.items()
        )

        d[schema.slice("o2")] = -r_o2  # this route's aa-gated O₂ share (off every ledger)
        for pool, rate_i in product_rates:
            d[schema.slice(pool)] = rate_i
        d[schema.slice("CO2")] = co2_mol * M_CO2
        d[schema.slice("N")] = nitrogen  # DEAMINATION: precursor N → ammonium (D-45)
        return d


class MaillardStrecker(Process):
    """Non-oxidative THERMAL Strecker aging: sugar dicarbonyls degrade amino acids, NO O₂ (D-87).

    The **O₂-independent thermal mirror** of :class:`StreckerDegradation` (D-75) — the beat that
    decision deferred ("the non-oxidative Maillard/sugar-dicarbonyl route (sweet wines, thermal) is
    deferred, keeping Strecker honestly on the ``o2`` sub-axis"). Where the D-75 oxidative route
    needs dissolved O₂ (its o-quinones are the amino-acid oxidant), THIS route is driven by
    **residual sugar + heat alone**: the sugar forms **α-dicarbonyls** (methylglyoxal, glyoxal,
    deoxyosones via Maillard/caramelization), and those dicarbonyls deaminate + decarboxylate amino
    acids to Strecker aldehydes with **no oxygen** — so a *sealed, sulfited, oxygen-free sweet wine*
    still ages, developing the Sauternes / Madeira / baked-wine aroma suite. This is the same
    O₂-independent-thermal-mirror relationship :class:`ThermalAnthocyaninFade` (D-83) has to
    :class:`AnthocyaninFading` (D-81).

    **The aroma suite — four new pools + two shared (the owner's full-scope D-87 choice).** The
    products are :data:`_MAILLARD_PRODUCTS`: **methional** (cooked-potato) and
    **phenylacetaldehyde** (honey) — the *same* two molecules the D-75 route makes, so the same two
    pools/thresholds are reused — plus four D-87-only pools: the three branched-chain **malty**
    Strecker aldehydes **2-methylbutanal** (isoleucine), **3-methylbutanal** (leucine),
    **2-methylpropanal** (valine), and **sotolon**, the curry/fenugreek/maple **furanone** that is
    *the* diagnostic marker of botrytized sweet wine / vin jaune / Madeira. (Scope, documented: the
    thermal-wine aldehyde signature is broader still — the branched-chain aldehydes + sotolon are
    the
    v1 set; other markers are lumped out. And sotolon is **not** a Strecker aldehyde — a
    threonine/acetaldehyde aldol furanone — so it carries no CO₂ term; its 2 acetaldehyde-derived
    carbons are lumped into its THREONINE draw, the acetaldehyde-coupled route deferred.)

    ``n_ald = k_maillard_strecker · f(T) · [S_total] · gate(aa)`` — first-order in the **residual
    sugar** (summed over the sugar vector, the dicarbonyl driver) and gated by the amino-acid
    availability ``gate = aa/(K_amino_acids + aa)`` (the D-33/D-75 smooth-Monod availability shape).
    ``f(T) = arrhenius_factor(T, E_a_maillard_strecker, T_ref)`` is the **strongly** warmer-faster
    factor: ``E_a_maillard_strecker ≈ 100 kJ/mol`` sits *above* the oxidative aging E_a's (~50), the
    sourced ordering that Maillard/caramelization out-accelerates oxidation with temperature — why a
    warm Madeira estufagem develops thermal character orders faster than cellar aging. ``n_ald`` is
    **mol total product/L/h directly** (the per-sugar dicarbonyl yield is folded into
    ``k_maillard_strecker``, avoiding a fake per-sugar molar conversion — the wine sugar is a
    glucose/fructose vector stand-in). It is split among the six products by their **normalized**
    composition weights ``w_maillard_*`` (relative amino-acid abundance × reactivity — a production
    flux, *not* potency, which the OAV thresholds already carry; the D-75 follow-up lesson), and the
    five decarboxylating aldehydes each contribute **1 mol CO₂** (sotolon none).

    **S is a read-only DRIVER, not consumed here (the forced carbon-closure choice).** The Strecker
    aldehyde's carbon skeleton **is** the amino acid (methional = methionine − COOH); the sugar
    dicarbonyl is only the electron-accepting oxidant, and its own carbon goes to melanoidin (booked
    by the separate :class:`Caramelization`, D-88). So the aldehyde carbon is drawn from
    ``amino_acids`` (the D-75 algebra), and ``S`` is **not** debited by this Process — which is not
    merely convenient but *forced*: booking a sugar draw here would break ``total_carbon``
    (melanoidin
    is off-ledger) *and* undercount real sugar loss (bulk thermal browning dominates depletion,
    D-88).
    The unbooked per-Strecker sugar consumption is µM-scale (µg/L aldehydes ⇒ trace dicarbonyl),
    negligible vs the g/L residual pool. So ``S`` is read but **not** in ``touches``.

    **Carbon + nitrogen close by construction (the D-75 idiom exactly).** Each precursor's draw
    is sized to **its own product's carbon** (its aldehyde, plus its own carboxyl CO₂ if it
    decarboxylates — sotolon does not), and all that precursor nitrogen is **deaminated** to ``N``
    (products N-free), so ``total_carbon`` and ``total_nitrogen`` close to machine precision.
    Because closure holds for *any* CO₂ attribution,
    the CO₂ is keyed to the :data:`_MAILLARD_PRODUCTS` ``decarboxylates`` flags explicitly and the
    produced µg/L levels are anchored to literature (no conservation test would catch a mis-key).
    Since D-100 each product draws its OWN precursor, so the pre-D-100 arginine-lump stand-in —
    exact on the ledger, approximate on provenance — is retired. ``total_mass`` ({S,E,CO2}) sees
    the CO₂ with no matching S/E debit, but is never asserted
    on an aging run (the standing :class:`OxidativeAcetaldehyde` scope-out).

    **Additive with the oxidative route over the shared ``amino_acids`` limiting reagent.** Both
    Strecker routes draw ``amino_acids`` and ``ProcessSet`` sums them, so the pool depletes *once*
    and
    splits by their rates — the o2-sharing pattern (D-73) applied to the amino-acid limiting
    reagent,
    no double-count. An O₂-**and**-sugar-rich aged sweet wine runs both; a sealed sweet wine runs
    only
    this one (the discriminating case, D-87); a dry oxidised wine runs only D-75.

    **Wine-only + isolable + doubly substrate-gated (prime directive #3).** ``amino_acids`` and the
    ``N``-deamination are wine-only (beer's amino-acid pool is not tracked, D-32), so — like
    :class:`StreckerDegradation` — this is wired into the *wine* medium only; the ``"amino_acids"
    not
    in schema`` guard is a hard no-op besides. Wired **disabled at the compile seam**;
    ``begin_aging``
    enables it with the other aging Processes. **Isolability rests on the ``amino_acids`` HARD
    gate**
    (undosed ⇒ exactly 0 ⇒ byte-for-byte the case without this Process — the default wine is
    unchanged); residual sugar is a **soft** driver (a "dry" wine still holds ~1–2 g/L, so the route
    is *negligible* there, not byte-for-byte zero — the physically-correct trace, not a claim of
    exact isolation on ``S``). **Writes ``N``** (deamination), so an enabled run drops structural
    ``tier_of("N")`` PLAUSIBLE→SPECULATIVE (the D-45/D-75 note). Tier **speculative** (the *form* —
    sugar-driven, amino-acid-gated, heat-accelerated, O₂-independent — is sourced; every magnitude
    is
    an order-of-magnitude estimate).
    """

    name = "maillard_strecker"
    tier = Tier.SPECULATIVE
    #: Writes the five thermal-route product pools + the decarboxylation ``CO2``, drawing the carbon
    #: from each product's OWN precursor (D-100) and deaminating its nitrogen to ``N``. There is NO
    #: ``o2`` term (the whole point). The C/N transfer closes exactly.
    #:
    #: **``S`` became a WRITE at D-104 and is READ-ONLY again at D-107.** D-104 gave sotolon a
    #: de-novo sugar draw because its carbon is not its precursor's; sotolon now lives in
    #: :class:`SotolonAldolCondensation` and takes its carbon from the tracked
    #: ``alpha_ketobutyrate``
    #: pool instead, so the sugar stand-in it needed is gone and every one of the five products left
    #: here takes every gram from its own amino acid. ``S`` is back to what it was through
    #: D-87/D-100: the **dicarbonyl driver**, read and never written. That is a fidelity statement,
    #: not a tidy-up — a Strecker degradation has no business drawing sugar carbon, and the only row
    #: that did was the row that was not a Strecker degradation.
    touches = (
        "methional",
        "phenylacetaldehyde",
        "2_methylbutanal",
        "3_methylbutanal",
        "2_methylpropanal",
        "CO2",
        "N",
        # Each product's OWN precursor (decision D-100) — no longer the lumped ``amino_acids``
        # pool. Five products over five distinct amino acids, one each.
        *dict.fromkeys(precursor for (_, _, _, precursor) in _MAILLARD_PRODUCTS),
    )
    #: ``k_maillard_strecker``/``E_a_maillard_strecker`` and the six ``w_maillard_*`` composition
    #: weights are this Process's own (thermal.yaml, D-87); ``K_amino_acids`` is the *shared*
    #: availability half-saturation (the same constant the mercaptan/reroute/D-75 gates read);
    #: ``T_ref`` is shared with every Arrhenius rate. Their tiers cap the output tiers via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_maillard_strecker",
        "E_a_maillard_strecker",
        "w_maillard_methional",
        "w_maillard_phenylacetaldehyde",
        "w_maillard_2_methylbutanal",
        "w_maillard_3_methylbutanal",
        "w_maillard_2_methylpropanal",
        "K_amino_acids",
        "T_ref",
        # Each precursor's must-spectrum share, which scales its relative-depletion gate (D-100).
        *dict.fromkeys(
            SPEC_BY_SPECIES[precursor].fraction_param for (_, _, _, precursor) in _MAILLARD_PRODUCTS
        ),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (beer tracks no amino acids, D-32; sotolon is a wine-only pool): a hard
        # no-op on any schema without them, belt-and-suspenders to the wine-only wiring.
        if "sotolon" not in schema or "leucine" not in schema:
            return d
        # Residual sugar is the dicarbonyl DRIVER (summed over the vector), a SOFT gate: a dry wine
        # (S ≈ 0) makes ~none, but the trace is physically real, so this is not an isolability
        # claim.
        s_total = max(float(y[schema.slice("S")].sum()), 0.0)
        if s_total <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_maillard_strecker"], params["T_ref"])
        driver = params["k_maillard_strecker"] * f_t * s_total  # mol product/L/h before gating
        if driver <= 0.0:
            return d

        # Split among the six products by their NORMALIZED composition weights (relative production
        # flux), each throttled by ITS OWN precursor's relative-depletion gate (decision D-100).
        # At must-spectrum composition every gate equals the old lumped gate, so this reduces
        # byte-for-byte to the D-87 rate; away from it, a product stops when the amino acid it is
        # made from runs out — which is the whole point (the D-100 pathology was sotolon dying
        # because a pool 38% ARGININE had been drained by fusels that eat leucine).
        weights = [params[wname] for (_, _, wname, _) in _MAILLARD_PRODUCTS]
        w_sum = sum(weights)  # > 0 (all speculative positive weights)
        co2_mol = 0.0  # mol CO2/L/h — every product here decarboxylates (D-107)
        product_rates: list[tuple[str, float]] = []
        # Carbon each precursor must supply — accumulated per species because the five products draw
        # from five different amino acids, and each product's CO2 must be attributed to the
        # precursor
        # that actually released it (its own carboxyl carbon), not to a lump.
        precursor_carbon: dict[str, float] = {}
        for (pool, m_i, _wname, precursor), w_i in zip(_MAILLARD_PRODUCTS, weights, strict=True):
            gate_i = depletion_gate(y, schema, params, (SPEC_BY_SPECIES[precursor],))
            if gate_i <= 0.0:
                continue  # this precursor is exhausted ⇒ its product stops (and cannot go negative)
            n_i = (w_i / w_sum) * driver * gate_i  # mol/L/h of product i
            if n_i <= 0.0:
                continue
            rate_i = n_i * m_i  # g/L/h
            product_rates.append((pool, rate_i))
            # 1 CO2 per Strecker decarboxylation — the amino acid's own carboxyl carbon, charged to
            # THIS product's precursor. Since D-107 this is unconditional: the one row that did not
            # decarboxylate was sotolon, which was never a Strecker aldehyde and now has its own
            # Process. Charging it is what makes the carbon-sized draw land on the true 1:1 molar
            # stoichiometry (D-105's signature) — and no conservation test would catch its absence,
            # so it is pinned by a driven test instead (D-87's trap, D-105's fix).
            co2_mol += _CO2_PER_STRECKER_ALDEHYDE * n_i
            carbon_i = rate_i * carbon_mass_fraction(pool) + (
                _CO2_PER_STRECKER_ALDEHYDE * n_i * M_CO2 * carbon_mass_fraction("CO2")
            )
            precursor_carbon[precursor] = precursor_carbon.get(precursor, 0.0) + carbon_i
        if not product_rates:
            return d

        # Draw each product's carbon from ITS OWN precursor, sized to match, and deaminate that
        # precursor's nitrogen to N (the D-45/D-75 idiom, speciated): carbon out of each amino acid
        # == carbon into the products + CO2 it supplied, and all their N lands in N (products are
        # N-free), so total_carbon and total_nitrogen both close to machine precision — now with
        # each molecule weighted at its own fractions rather than arginine's.
        nitrogen = sum(
            draw_precursor_carbon(d, schema, precursor, carbon)
            for precursor, carbon in precursor_carbon.items()
            if carbon > 0.0
        )
        for pool, rate_i in product_rates:
            d[schema.slice(pool)] = rate_i
        d[schema.slice("CO2")] = co2_mol * M_CO2
        d[schema.slice("N")] = nitrogen  # DEAMINATION: precursor N → ammonium (D-45)
        return d


class SotolonAldolCondensation(Process):
    """Sotolon by aldol condensation of α-ketobutyrate + acetaldehyde — the node's consumer (D-107).

    ``n_sot = k_sotolon_aldol · f(T) · [alpha_ketobutyrate] · [acetaldehyde]``  [mol sotolon/L/h],
    **bimolecular (mass-action) in its two real substrates**, with::

        d(sotolon)/dt             = +n_sot · M_SOTOLON
        d(alpha_ketobutyrate)/dt  = −n_sot · M_ALPHA_KETOBUTYRATE     (the C4 half)
        d(acetaldehyde)/dt        = −n_sot · M_ACETALDEHYDE           (the C2 half)

    **This is a purely chemical reaction and that is the whole argument.** Pham *et al.* 1995 formed
    sotolon from α-ketobutyric acid + acetaldehyde in a synthetic medium "under the conditions of
    temperature, pH, and alcoholic content similar to the ageing of wine in barrels … by a purely
    chemical mechanism". No enzyme, no yeast, no sugar. Sotolon accumulates in a **sealed bottle
    over
    years**, which is why its α-ketobutyrate must be the *excreted, extracellular* residual
    (:mod:`~fermentation.core.kinetics.keto_acids`, D-107) and not an intracellular flux
    intermediate — the same test D-49 applied to pyruvate, answered the other way because the
    reaction that consumes this pool happens where no intracellular pool can reach.

    **WHY THIS IS NOT IN :class:`MaillardStrecker` ANY MORE (the D-107 move).** Sotolon rode the
    thermal Strecker Process from D-87 to D-106 carrying two exception flags — ``decarboxylates=
    False`` and ``de_novo=True`` — and D-104 noted that "the two flags travel together for the same
    underlying reason: sotolon is not made the way the other five are". They were the model
    reporting,
    in the only vocabulary it had, that this compound was in the wrong Process. It is not a Strecker
    degradation: it is not a degradation *of an amino acid* at all, it releases no CO₂ because there
    is no carboxyl to lose, and its carbon was never its "precursor's". Moving it here deleted both
    flags and left :data:`_MAILLARD_PRODUCTS` five true Strecker aldehydes with no exceptions in it.

    **The sugar driver is gone, and Pons *et al.* 2010 is why.** Inside ``MaillardStrecker``
    sotolon's
    rate was ``k · f(T) · S`` — pseudo-first-order in residual sugar, the *dicarbonyl* driver. A
    dicarbonyl is what a Strecker degradation needs; an aldol needs its two substrates. The sugar
    driver made sotolon a **sweetness** marker, and that is empirically wrong: Pons *et al.* 2010
    identified this exact aldol as the source of sotolon in prematurely aged **DRY white wines**,
    where it is the *prémox* marker. Sugar does not appear in this rate law because it does not
    appear in the reaction.

    **THE EMERGENT PAYOFF — oxidation raises sotolon, with nothing scripted.** ``acetaldehyde`` is
    raised by :class:`OxidativeAcetaldehyde` (D-71) when O₂ is dosed, so an oxidised wine makes more
    sotolon *because* it made more acetaldehyde. That is Pons' prémox mechanism falling out of a
    rate law written from Pham's chemistry, and it is measured (D-107), not asserted. The same term
    covers vin jaune, where the flor's very high acetaldehyde is why 6 years under veil reaches
    120–268 µg/L (Pham).

    **SO₂ — READS FREE, NOT TOTAL (decision D-108, correcting D-107).** D-107 claimed here that "SO₂
    also enters correctly for free: it binds acetaldehyde (D-47), and this Process reads the pool
    the binding depletes". **That sentence was false and the code matched the sentence, not the
    chemistry.** The binding depletes *nothing*: the ``acetaldehyde`` slot holds TOTAL acetaldehyde
    and :func:`~fermentation.core.acidbase.free_acetaldehyde` derives the free share as
    ``free = total − bound`` — a read-only overlay on the D-51 equilibrium. So the rate read the
    total, *including the bound share*, and SO₂ — which strands acetaldehyde by protecting it from
    ADH — came out **raising** sotolon: measured, a dry wine dosed 60 mg/L SO₂ in the must went
    0.025 → 5.02 µg/L, a 200× rise to threshold. Reality runs the other way (Pons: low free SO₂ is
    the *prémox risk factor*). The adduct's carbonyl is blocked and an aldol **is** a nucleophilic
    attack on that carbonyl — the same argument :class:`AcetaldehydeBridging` (D-80) and the tannin
    polymerization already make for the ethylidene bridge, and :class:`~fermentation.core.kinetics.
    acetaldehyde.AcetaldehydeReduction` (D-47) makes for ADH. **This Process was the only reader of
    that pool in the file getting it wrong, and the file already carried the argument three times.**
    No literature was needed to convict — only the sibling Processes (the D-105 internal-
    contradiction shape).

    **THE SECOND EMERGENT PAYOFF — PONS' PRÉMOX MECHANISM ITSELF, AND IT ARRIVED BY CORRECTING A
    CLAIM THIS ENTRY FIRST GOT WRONG (D-108).** The suppression above is *not* permanent, and the
    first draft of this docstring said it was: ``so2_total`` looked like a conserved slot because a
    60 mg/L must dose measured **60.0000 mg/L at day 729** — but that probe dosed **no O₂**, and
    :class:`SulfiteOxidation` (which consumes ``so2_total``) needs O₂ as its substrate. It measured
    a sealed bottle and read that as a property of the slot: **a vacuous measurement that agreed
    with me** (D-106's lesson, one beat later). Dosed with O₂ the SO₂ *does* deplete, and sotolon
    recovers as it goes — measured, dry + 60 mg/L must SO₂, O₂ 0/5/20/60 mg/L::

        SO₂ end            60.000 → 44.597 → 30.369 → 24.508  mg/L
        sotolon             0.059 →  0.121 →  2.113 →  7.639  µg/L
        unsulfited control  0.025 →  0.677 →  2.632 →  7.837  µg/L
        PROTECTION RATIO        —  →  0.179 →  0.803 →  0.975   ← the evidence

    So the model now says what Pons says: **a sulfited wine is protected while its SO₂ lasts, and
    goes prémox as the free SO₂ fades** — because the O₂ that drives the sotolon has eaten the SO₂
    that was suppressing it. Nothing scripted; it is D-108's free-read composed with D-72's
    oxidation. **Read the RATIO row, not the sotolon row**: the raw ladder rises monotonically in
    the unsulfited control *too* (O₂ → acetaldehyde → sotolon is the dominant driver either way), so
    it would rise even if the protection never eroded. Only the sulfited/unsulfited ratio at matched
    O₂ divides that driver out, and it is what the test asserts (0.179 → 0.975 = protection spent).
    **The real bound is one layer out and is NOT this Process's**: a *sealed* wine
    here has strictly zero O₂ ingress (no closure permeation — the same gap D-102 named for DMS
    leaving through the closure), so a sealed sulfited bottle never ages toward prémox at all. That
    is the limitation to state, not "SO₂ is permanent", which was false. At O₂ = 0 the sulfited wine
    sits marginally *above* the unsulfited one (0.059 vs 0.025) — both far sub-threshold — because
    SO₂ stranded acetaldehyde that ADH would otherwise have cleared, and a small free share of it
    condenses. That is the D-47 protection showing through, not a residue of the D-107 bug.

    **Isolability is EXACT and it is why this is mass-action rather than gated.** The rate is the
    product of its two substrates, so a ProcessSet without ``_KETO_ACID_PROCESSES`` leaves
    ``alpha_ketobutyrate`` at 0 and this rate at **exactly** 0 — byte-for-byte, with no clamp, no
    availability constant, and no epsilon. Neither pool can be driven negative for the same reason
    (each draw is proportional to the pool it debits). A Michaelis-style availability gate would
    have
    needed a fabricated half-saturation constant per substrate; the honest rate law needs none — the
    **more faithful form is also the one with fewer invented numbers** (contrast the D-98 trap).

    **Carbon closes on ATOM COUNTS, not on a sized draw** — the strongest closure in the tree.
    Sotolon's 6 carbons are α-ketobutyrate's 4 + acetaldehyde's 2, so ``total_carbon`` closes
    because
    the *chemistry* balances, not because a mass was sized to make it. This route therefore has no
    D-105 signature to check: there is no carbon-sized draw here to be silently wrong about the mole
    count. It draws **1 mol of each substrate per mol of sotolon** because that is what is written.
    Nitrogen: untouched — both substrates are N-free, so this is the one aging aroma route that
    releases **no ammonium** (its nitrogen left earlier, at the threonine deamination that made the
    keto-acid). Mass carries a small gap: the aldol condenses out 1 mol H₂O
    (``C4H6O3 + C2H4O → C6H8O3 + H2O``), which is untracked — carbon is the invariant, as
    everywhere.

    **The temperature term is a CARRY-OVER, and is labelled one (the D-101/D-102 line).** It reads
    ``E_a_maillard_strecker`` — the constant sotolon already rode from D-87 through D-106 — so this
    is the status quo, not a new number. Pham sources the **direction** ("the formation of sotolon
    increases by increasing temperature") and nothing more; inventing a sotolon-specific activation
    energy to look precise is exactly the E_a D-101 fabricated and D-102 had to retract. Re-sourcing
    it is a parameter-file change, not a structural one.

    **D-107'S "LOAD-BEARING DEPENDENCY ON D-27" WAS A COMMENSURABILITY ERROR — RETIRED, MEASURED
    (decision D-108).** D-107 warned here that the dry arm's ``acetaldehyde ≈ 0.0`` mg/L was a D-27
    artifact because **"real dry whites hold ~30 mg/L"**, and that correcting D-27 would collapse
    the sweet-vs-dry separation. **The ~30 mg/L is a SULFITED figure and every run behind that
    warning was UNSULFITED.** Against the like-for-like target the model is nearly right: an
    unsulfited white really ends at **2.7 mg/L** (Herzan *et al.* 2020, Food Sci. Nutr.
    8:5850–5859, **PMC7684598**, Table 1, variant (0/0/0) — no SO₂ at any of must / maturation /
    bottling), and finished-wine acetaldehyde is driven by the **SO₂ regime**, not by the ADH gate
    — 2.7 → 6.5 → 17.2 → 25.9 → 51.6 mg/L as SO₂ is added, the paper concluding that to minimise
    acetaldehyde "it is recommended to exclude the use of SO₂". The model
    already carries that driver (D-47) and reproduces the ladder: dosed 60 mg/L in the must it
    strands **22.9 mg/L** against that study's 17.2. **D-27 is acquitted on the gate.** Its residual
    0.0-vs-2.7 floor is a real but small miss, named here and deliberately not patched — inventing a
    constant to hit 2.7 buys nothing observable (sotolon is sub-threshold either way), which is
    D-107's own lesson (v). The corrected sweet figure is **38.066**, not the 32.9 D-107 recorded.
    ⚠ **The separation is fragile for a different reason, and it is not latent — it is LIVE**: the
    collapse D-107 predicted for some future D-27 fix is reachable **today** by dosing SO₂, which
    every real winemaker does. That is what D-108 fixed above, and it was this Process's bug, not
    D-27's.

    **Scope — named, not hidden.** (1) Pham also measures sotolon rising with **decreasing pH** and
    **decreasing ethanol**; neither term is here, and the model has both quantities, so this is a
    real omission rather than an inexpressible one. (2) Botrytised sweet wines may carry an
    additional sugar/furanone contribution this Process does not model — what is modelled is the
    sourced aldol. (3) The **enantiomers are not split**: (S)-sotolon's threshold (~0.8 µg/L) is
    ~100× below (R)'s (~89 µg/L) (Pons *et al.*), so one racemic pool against one threshold is a
    lump. (4) Propanol is still made from threonine directly rather than from this pool, so the
    genuine propanol-vs-sotolon competition — which is over **α-ketobutyrate**, not threonine — is
    not yet expressible; that is the fusel side of the node (see D-107).

    **D-109 AUDITED (4) AND UPHELD IT — then rejected the fix it was read as prescribing.** The
    sentence above is *correct*: the competition really is absent, not merely mis-shaped, because
    :class:`~fermentation.core.kinetics.keto_acids.AlphaKetobutyrateExcretion`'s rate is flux-only
    (threonine's gate moves the carbon source, not the rate — measured bit-identical). But it was
    read as prescribing "draw propanol from ``alpha_ketobutyrate``", and **that pool is the wrong
    one, by this Process's own argument**. This aldol selects the *excreted, extracellular residual*
    because it runs in a sealed bottle where no intracellular pool can reach — D-49's test.
    Propanol is made *intracellularly, by living yeast, mid-ferment*: it IS D-49's flux
    intermediate, so the same test that chose this pool for sotolon **rejects it for propanol**.
    The arithmetic agrees and is blunter: propanol's molar demand is **~2.8× the total α-KB this
    pool ever excretes**, so the re-base is infeasible, not merely mis-attributed — it would starve
    propanol *and* collapse this Process's substrate, while looking like the promised competition
    at last. That same 2.8× is why the item is **relocated, not dissolved**: propanol's 2-KB demand
    is that much larger than *the excretion flux itself*, so an *intracellular* node — a
    quasi-steady **flux partition**, not a pool (D-49's own physics) — would couple propanol and
    sotolon materially. (Deliberately **not** "propanol is the node's dominant sink": 2-KB's
    committed route is isoleucine biosynthesis via KMV, which this model does not carry, so the
    share of *total* 2-KB synthesis is unmeasured. The ratio is what the measurement supports.)
    See D-109.

    Wine-only (both product and substrate slots are wine-only); wired **disabled at the compile
    seam** and enabled by ``begin_aging`` with the other aging Processes. Tier **speculative** (the
    *form* — bimolecular in α-ketobutyrate × acetaldehyde, purely chemical, warmer-faster — is
    sourced; the rate magnitude is an order-of-magnitude estimate).
    """

    name = "sotolon_aldol_condensation"
    tier = Tier.SPECULATIVE
    #: Fills ``sotolon`` from its two real substrates. No ``S`` (not sugar-driven), no ``CO2`` (no
    #: decarboxylation), no ``N`` (both substrates are nitrogen-free) — the three absences are the
    #: content of D-107.
    touches = ("sotolon", "alpha_ketobutyrate", "acetaldehyde")
    #: ``k_sotolon_aldol`` is this Process's own second-order rate constant (thermal.yaml, D-107);
    #: ``E_a_maillard_strecker`` is a **carry-over**, not a claim that this reaction shares the
    #: Strecker activation energy (see the class docstring); ``T_ref`` is shared with every
    #: Arrhenius rate. Their tiers cap the output tier via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_sotolon_aldol", "E_a_maillard_strecker", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots: a hard no-op on any schema without them (belt-and-suspenders to the
        # wine-only wiring, the MaillardStrecker idiom).
        if "sotolon" not in schema or "alpha_ketobutyrate" not in schema:
            return d
        keto = max(float(y[schema.slice("alpha_ketobutyrate")][0]), 0.0)
        acetaldehyde = max(float(y[schema.slice("acetaldehyde")][0]), 0.0)
        # Mass-action in BOTH substrates ⇒ either at zero means no reaction. This is the exact
        # isolability guarantee (keto-acid Processes off ⇒ pool 0 ⇒ byte-for-byte no-op), not a
        # defensive guard: the clamps only absorb solver undershoot.
        if keto <= 0.0 or acetaldehyde <= 0.0:
            return d
        # SO₂-BOUND ACETALDEHYDE CANNOT CONDENSE (decision D-108, the D-47/D-80 precedent). The
        # `acetaldehyde` slot holds TOTAL acetaldehyde — `free_acetaldehyde` derives the free share
        # from the D-51 equilibrium (`free = total − bound`); the binding does NOT deplete the slot.
        # The adduct's carbonyl is blocked, and an aldol condensation IS a nucleophilic attack on
        # that carbonyl — the identical argument AcetaldehydeBridging (D-80) and the tannin
        # polymerization already make for the ethylidene bridge, and AcetaldehydeReduction (D-47)
        # for ADH. The RATE reads free; the DRAW below still debits the total slot (consuming free
        # acetaldehyde removes it from the total and the equilibrium re-splits) — the D-47 idiom.
        # The guard is EXACT: an unsulfited run pays no per-RHS pH `brentq` and is byte-for-byte the
        # D-107 total-acetaldehyde case, so every output D-107 measured is unmoved.
        if SO2_STATE_KEY in schema and float(y[schema.slice(SO2_STATE_KEY)][0]) > 0.0:
            acetaldehyde = free_acetaldehyde(y, schema, params, ph_of_state(y, schema, params))
            if acetaldehyde <= 0.0:  # all acetaldehyde bound ⇒ no free carbonyl to condense
                return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_maillard_strecker"], params["T_ref"])
        # Second-order rate law (Pham et al. 1995): mol/L/h from the two molar concentrations.
        n_sot = (
            params["k_sotolon_aldol"]
            * f_t
            * (keto / M_ALPHA_KETOBUTYRATE)
            * (acetaldehyde / M_ACETALDEHYDE)
        )
        if n_sot <= 0.0:
            return d
        # 1 mol of each substrate per mol of sotolon. Carbon closes because 4 + 2 == 6 on the atom
        # counts in `chemistry`, so no draw is being *sized* to close it (D-105's blind spot cannot
        # exist here).
        d[schema.slice("sotolon")] = n_sot * M_SOTOLON
        d[schema.slice("alpha_ketobutyrate")] = -n_sot * M_ALPHA_KETOBUTYRATE
        d[schema.slice("acetaldehyde")] = -n_sot * M_ACETALDEHYDE
        return d


class Caramelization(Process):
    """Non-oxidative THERMAL browning: residual sugar caramelizes → melanoidin, NO O₂ (D-88).

    The **O₂-independent thermal mirror** of :class:`PhenolicBrowning` (D-74) — the browning half of
    the non-oxidative thermal axis :class:`MaillardStrecker` (D-87) opened. Where D-74 needs
    dissolved O₂ to oxidise phenolics brown, this route browns **residual sugar** by heat alone
    (thermal dehydration/caramelization to melanoidin), so a *sealed, oxygen-free sweet wine* still
    darkens with age — the amber-to-brown of an aged Sauternes / the deep colour of Madeira and
    baked/rancio wines. It raises the **same** ``A420`` browning index D-74 accumulates (both routes
    darken the wine — oxidative *and* thermal browning are one observable), so it needs **no** new
    observable, only a carbon-park pool for the sugar it consumes.

    **CARAMELIZATION, not Maillard (the D-88 scope honesty).** This is the **sugar-only** route: it
    touches ``{S, melanoidin, A420}`` and incorporates **no amino-acid nitrogen**. True Maillard
    melanoidins are nitrogen-bearing (sugar + amino acid); modelling that N-incorporating browning
    is
    deferred. So ``melanoidin`` here is a caramelization polymer (nitrogen-free), and this Process
    is
    honestly *caramelization*.

    ``d(S)/dt = −r`` with ``r = k_caramelization · f(T) · [S_total]`` — first-order in the residual
    sugar (summed over the vector), ``f(T) = arrhenius_factor(T, E_a_caramelization, T_ref)`` the
    **strongly** warmer-faster factor (``E_a_caramelization ≈ 100 kJ/mol``, above the oxidative
    aging
    E_a's — the same sourced Maillard/caramelization-≫-oxidation ordering :class:`MaillardStrecker`
    carries). The consumed sugar is booked as::

        d(melanoidin)/dt = +r · c(sugar) / c(melanoidin)          # carbon-exact transfer
        d(A420)/dt       = +y_a420_per_melanoidin · d(melanoidin)/dt   # the shared browning index

    **The FIRST aging Process to consume core ``S`` — an on-ledger carbon-park (the forced
    closure).**
    Every prior aging Process touches aroma pools / ``o2`` / ``amino_acids`` / ``N`` / ``E`` — none
    the core sugar. Because ``S`` is **on** ``total_carbon``, the sugar carbon this Process draws
    **must** land in a weighted pool or the transfer would read as carbon destroyed (unlike D-74's
    ``A420``, whose pigment carbon comes from an *untracked* phenol pool, so it is off-ledger). So
    ``melanoidin`` is an **on-ledger carbon-park** (the ``debris``/``glucan`` precedent, D-34): the
    carbon released from ``S`` (at the sugar's fraction) is redeposited into ``melanoidin`` (at its
    caramelan-stand-in fraction), so ``total_carbon`` closes to machine precision (the
    :class:`EsterHydrolysis` carbon-exact split). The water lost on dehydration is the standing
    aging-axis mass gap (``total_mass`` weights only ``{S, E, CO2}``, never asserted on an aging
    run);
    CO₂/volatile evolution of real caramelization is lumped into the polymer (a documented v1
    simplification). ``A420`` is the optical browning index (off every ledger, the D-74 slot), so it
    carries no carbon — only ``melanoidin`` parks it.

    **Medium-agnostic (D-90 supersedes D-88's provisional "wine-only v1").** Sugar-only browning is
    physics, not a grape property — beer/wort residual dextrins (unfermented maltose/maltotriose)
    caramelize and brown an aged/warm-stored beer exactly as sweet-wine sugar does — so D-90 wires
    this into **both** media (the ``melanoidin`` carbon-park is appended to both ``wine_schema`` and
    ``beer_schema``, the D-86 oak-to-beer pattern). The ``"melanoidin" not in schema`` guard is now
    only belt-and-suspenders. The **vectorized** draw (the D-90 rework) apportions the sugar debit
    across beer's 3-slot ``S`` by each component's share and releases carbon at each component's
    **own** fraction (glucose / maltose / maltotriose differ), so ``total_carbon`` closes for both
    media; wine's single slot reduces to the D-88 form byte-for-byte. The N-incorporating
    :class:`MaillardBrowning` (D-89) does **not** follow — beer's ``amino_acids`` pool is untracked
    (D-32), so beer thermal browning is caramelization only.

    **Isolable + a SOFT sugar gate (prime directive #3).** Wired **disabled at the compile seam**;
    ``begin_aging`` enables it with the other aging Processes. Residual sugar is a **SOFT** driver:
    a dry beverage (``S ≈ 0`` at the aging segment — a standard wine ferments to dryness, ``S ≤ 0``
    exactly, before ``begin_aging``; a standard beer to ``S ≈ 5e-11``, browning a negligible trace)
    is byte-for-byte (wine) or numerically inert (beer) via the per-component clamp + ``S ≤ 0``
    guard, so an ordinary *dry* aged beverage is unchanged. But a **sweet** wine / **high-residual**
    (under-attenuated, big-stout) beer now browns thermally even sealed and sulfited — so the
    D-71/D-74 "reductive aging is byte-for-byte the ester-only case" claim now holds only for *dry*
    beverages (the D-83-style supersession: a sealed sweet wine / residual-sugar beer is *not*
    inert). Tier **speculative** (the *form* — sugar-driven, heat-accelerated, O₂-independent
    browning — is sourced; the rate and per-melanoidin absorbance yield are order-of-magnitude
    estimates).
    """

    name = "caramelization"
    tier = Tier.SPECULATIVE
    #: Consumes core ``S`` and books the carbon into the on-ledger ``melanoidin`` carbon-park (so
    #: ``total_carbon`` closes), raising the shared off-ledger ``A420`` browning index. Touches
    #: those
    #: three and nothing else — NO ``o2`` (the whole point), no amino acids (sugar-only
    #: caramelization,
    #: not Maillard).
    touches = ("S", "melanoidin", "A420")
    #: ``k_caramelization``/``E_a_caramelization``/``y_a420_per_melanoidin`` are this Process's own
    #: (thermal.yaml, D-88); ``T_ref`` is shared with every Arrhenius rate. Their tiers cap the
    #: ``S``/``melanoidin``/``A420`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_caramelization",
        "E_a_caramelization",
        "y_a420_per_melanoidin",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Medium-agnostic (D-90): a hard no-op on any schema without the melanoidin carbon-park.
        if "melanoidin" not in schema:
            return d
        s_slice = schema.slice("S")
        # PER-COMPONENT clamp (the D-90 correctness pin), NOT ``max(sum, 0)``: on beer's 3-slot S a
        # solver undershoot can leave ONE component slightly negative while the sum stays positive —
        # then ``frac_i = y_i / s_total`` would be negative and the apportioned debit
        # ``-r_sugar·frac_i`` would flip POSITIVE, silently *creating* that sugar (carbon still
        # closes for either sign — the D-89-denominator trap family). Clamping each component first
        # makes a negative slot contribute a zero draw. For WINE's single slot this is identical to
        # ``max(sum, 0)``, so wine is byte-for-byte unchanged.
        s_clamped = y[s_slice].clip(min=0.0)  # per-component ≥ 0
        s_total = float(s_clamped.sum())  # residual sugar, summed over the (clamped) vector
        # Dry beverage (S ≈ 0 at the aging segment) ⇒ byte-for-byte inert (the SOFT sugar gate; the
        # clamp above already absorbed any solver undershoot). A sweet wine / high-residual beer
        # browns; a dry one is unchanged.
        if s_total <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_caramelization"], params["T_ref"])
        r_sugar = params["k_caramelization"] * f_t * s_total  # g sugar/L/h caramelized (total)
        # Carbon-exact transfer sugar → melanoidin: apportion the total draw across the S vector by
        # each component's share and release its carbon at that component's OWN fraction (glucose /
        # maltose / maltotriose differ), then redeposit into melanoidin at its caramelan-stand-in
        # fraction — so total_carbon closes to machine precision (the EsterHydrolysis split) in both
        # media. Wine's single hexose slot ⇒ the sum is one term at c(glucose), the D-88 form.
        species = sugar_species(schema)
        carbon_released = r_sugar * sum(
            float(s_clamped[i]) / s_total * carbon_mass_fraction(sp) for i, sp in enumerate(species)
        )  # g C/L/h
        mel_rate = carbon_released / carbon_mass_fraction(_MELANOIDIN_SPECIES)  # g melanoidin/L/h
        # Debit the residual sugar, APPORTIONED across the vector (each slot loses its share of the
        # total draw). ``s_clamped / s_total`` sums to 1, so the total debit is exactly ``r_sugar``;
        # for wine's single slot it is ``-r_sugar`` as before.
        d[s_slice] = -r_sugar * (s_clamped / s_total)
        d[schema.slice("melanoidin")] = mel_rate
        # The shared A420 browning index rises with the melanoidin formed (off every ledger, the
        # D-74
        # optical-index slot). Monotone (mel_rate ≥ 0), like the D-74 oxidative browning it joins.
        d[schema.slice("A420")] = params["y_a420_per_melanoidin"] * mel_rate
        return d


class MaillardBrowning(Process):
    """Non-oxidative THERMAL browning: sugar + amino acid → N-bearing melanoidin, NO O₂ (D-89).

    The **amino-acid-incorporating** thermal-browning route that the sugar-only
    :class:`Caramelization` (D-88) explicitly deferred ("modelling that N-incorporating browning is
    deferred … ``melanoidin`` here is a caramelization polymer, nitrogen-free"). Where D-88 browns
    **sugar alone** to a nitrogen-free caramelan, *true* Maillard browning condenses a **reducing
    sugar with an amino acid** (Amadori rearrangement → Strecker/Maillard cascade → brown polymer)
    and **retains the amino-acid nitrogen in the melanoidin** — that nitrogen is what makes a
    Maillard melanoidin *nitrogenous*. So this route consumes **both** core ``S`` and
    ``amino_acids``
    by heat with **no O₂**, and a sealed sweet wine on residual sugar + amino acids browns thermally
    this way too. It raises the **same** ``A420`` browning index D-74/D-88 accumulate (all browning
    is one observable), needing **no** new observable — only the new carbon+nitrogen-park pool.

    **The three thermal amino-acid/browning branches, and why they don't double-count.** With D-89
    the non-oxidative thermal axis splits the amino-acid fate cleanly into complementary branches
    that ``ProcessSet`` sums over the shared ``amino_acids`` (and ``S``) reagents — the O₂-sharing
    pattern (D-73) applied to two limiting reagents, so each pool depletes *once* and splits by the
    branches' rates:

    * :class:`Caramelization` (D-88) — **sugar-only** browning → nitrogen-free ``melanoidin``; runs
      even at zero amino acids.
    * :class:`MaillardBrowning` (**this**, D-89) — the **N-retaining** browning branch: sugar +
      amino acid → nitrogen-bearing ``maillard_melanoidin``, **all** drawn amino-acid nitrogen kept
      in the polymer.
    * :class:`MaillardStrecker` (D-87) — the **N-releasing / volatile** branch: the same sugar
      dicarbonyls deaminate + decarboxylate amino acids to volatile Strecker aldehydes + sotolon,
      refunding nitrogen to ``N`` and evolving CO₂.

    Real Maillard chemistry partitions amino-acid nitrogen between polymer-retention and
    Strecker-release; the **system** (D-87 + D-89) reproduces that partition while each branch stays
    internally pure — the N-retaining polymer branch keeps *its* nitrogen (D-89), the deaminating
    volatile branch releases *its* nitrogen (D-87). Putting a partial-deamination split inside D-89
    would double-count D-87's release (and add an un-pinnable free parameter — the D-75/D-87
    silent-mis-key hazard), so **all drawn nitrogen is retained here** (the owner's D-89
    closest-to-reality choice, given D-87 already owns the release branch).

    ``r_sugar = k_maillard_browning · f(T) · [S_total] · gate(aa)`` — first-order in the **residual
    sugar** (summed over the vector, the reducing-sugar substrate) and gated by amino-acid
    availability ``gate = aa/(K_amino_acids + aa)`` (the D-33/D-75 smooth-Monod shape), with ``f(T)
    =
    arrhenius_factor(T, E_a_maillard_browning, T_ref)`` the **strongly** warmer-faster factor
    (``E_a ≈ 100 kJ/mol``, above the oxidative aging E_a's — the same sourced
    Maillard/caramelization-≫-oxidation ordering :class:`MaillardStrecker`/:class:`Caramelization`
    carry). So it is amino-acid-limited: as amino acids exhaust the browning stops, while D-88
    caramelization keeps browning the residual sugar. The consumed sugar + amino acid are booked
    as::

        d(S)/dt                  = −r_sugar
        d(amino_acids)/dt        = −r_aa
        d(maillard_melanoidin)/dt = +r_m
        d(A420)/dt               = +y_a420_per_maillard_melanoidin · r_m

    **Carbon AND nitrogen close by construction — the draws are sized to the melanoidin formed (the
    first two-ledger aging transfer).** ``maillard_melanoidin`` is an on-ledger carbon+nitrogen-park
    (the ``melanoidin`` carbon-park extended to nitrogen — the FIRST non-biomass, non-arginine
    species on ``total_nitrogen``): its stand-in ``C8H12O5N`` fixes its carbon fraction ``c_m`` and
    nitrogen fraction ``n_m``. Requiring **all** the amino-acid nitrogen to land in the polymer and
    **all** the drawn carbon (sugar + amino acid) to land in the polymer gives two equations::

        nitrogen:  r_aa · n(arg)                 = r_m · n_m
        carbon:    r_sugar · c(sugar) + r_aa · c(arg) = r_m · c_m

    solved (given ``r_sugar`` from the rate law) as::

        r_m  = r_sugar · c(sugar) / (c_m − n_m · c(arg)/n(arg))
        r_aa = r_m · n_m / n(arg)

    so ``total_carbon`` **and** ``total_nitrogen`` close to machine precision for *any* stand-in
    formula. **The denominator is the one silent trap** (advisor's must-check): ``c_m − n_m·c(arg)/
    n(arg)`` must be comfortably positive or ``r_m`` flips sign and the Process would *create* sugar
    with no conservation test catching it (closure holds for either sign). The threshold is
    mass-ratio
    ``c_m/n_m > c(arg)/n(arg) = 72/56 ≈ 1.29`` (atomic C:N > ~1.5); the C-rich melanoidin (C:N ≈
    8:1,
    ``c_m/n_m ≈ 6.9``) clears it by ~5×, leaving the denominator ≈ 0.81·c_m (healthy, no blow-up).
    The amino-acid draw is the identity-agnostic ``{arginine, generic}`` blend (D-100) — melanoidin
    retains amino-acid nitrogen regardless of which molecule carried it, so this route never touches
    a precursor;
    the whole amino acid (carbon skeleton + nitrogen) is built into the polymer here — unlike D-87,
    where the carbon leaves as aldehyde+CO₂ and only the nitrogen returns. ``total_mass``
    ({S,E,CO2})
    sees the ``S`` debit with no matching credit (the melanoidin and the dehydration water are
    off it), the standing aging-axis mass gap, never asserted on an aging run.

    **Wine-only + isolable + doubly substrate-gated (prime directive #3).** ``amino_acids`` and the
    ``maillard_melanoidin`` park are wine-only (beer's amino-acid pool is not tracked, D-32; beer
    thermal browning is deferred, the D-86 oak-to-beer pattern), so — like :class:`MaillardStrecker`
    / :class:`Caramelization` — this is wired into the *wine* medium only; the guard on either slot
    is
    a hard no-op besides. Wired **disabled at the compile seam**; ``begin_aging`` enables it with
    the
    other aging Processes. **Isolability rests on the ``amino_acids`` HARD gate** (undosed ⇒ exactly
    0
    ⇒ byte-for-byte the case without this Process — the default wine is unchanged); residual sugar
    is
    a **soft** driver (a dry wine holds ~0 at the aging segment ⇒ inert via the ``S ≤ 0`` guard, but
    a
    "dry" wine's residual ~1–2 g/L would brown a negligible trace, not byte-for-byte zero — the
    physically-correct trace). Tier **speculative** (the *form* — sugar + amino acid + heat →
    N-bearing browning, O₂-independent, strongly warmer-faster — is sourced; every magnitude is an
    order-of-magnitude estimate).
    """

    name = "maillard_browning"
    tier = Tier.SPECULATIVE
    #: Consumes core ``S`` AND ``amino_acids`` and books both into the on-ledger
    #: ``maillard_melanoidin`` carbon+nitrogen-park (so ``total_carbon`` and ``total_nitrogen`` both
    #: close), raising the shared off-ledger ``A420`` browning index. Touches those four and nothing
    #: else — NO ``o2`` (the whole point), NO ``CO2``/``N`` (all carbon+nitrogen retained in the
    #: polymer; the deaminating/decarboxylating branch is :class:`MaillardStrecker`, D-87).
    touches = (
        "S",
        "maillard_melanoidin",
        "A420",
        # The identity-agnostic pools (decision D-100): melanoidin retains amino-acid NITROGEN
        # without regard to which molecule carried it, so this route draws {arginine, generic} —
        # the same substrate as the yeast/MLF/Brett swaps — and never touches a precursor. That
        # is what makes it immune to the Ehrlich re-route, resolving the D-100 starvation.
        *(spec.pool for spec in ASSIMILABLE_SPECS),
    )
    #: ``k_maillard_browning``/``E_a_maillard_browning``/``y_a420_per_maillard_melanoidin`` are this
    #: Process's own (thermal.yaml, D-89); ``K_amino_acids`` is the *shared* availability
    #: half-saturation (the same constant the mercaptan/reroute/D-75/D-87 gates read), scaled by the
    #: assimilable pools' combined must-spectrum share (D-100); ``T_ref`` is shared with every
    #: Arrhenius rate. Their tiers cap the output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_maillard_browning",
        "E_a_maillard_browning",
        "y_a420_per_maillard_melanoidin",
        "K_amino_acids",
        "T_ref",
        *(spec.fraction_param for spec in ASSIMILABLE_SPECS),
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (beer tracks no amino acids, D-32; the melanoidin park is a wine slot):
        # a hard no-op on any schema without them, belt-and-suspenders to the wine-only wiring.
        if "maillard_melanoidin" not in schema or "amino_acids" not in schema:
            return d
        # HARD amino-acid gate on the identity-agnostic substrate (decision D-100) — the
        # isolability guarantee: an undosed wine is byte-for-byte the case without this Process.
        # → 0 as the pools empty, so the draw can never drive either negative, and the ``<= 0``
        # short-circuit also absorbs solver undershoot.
        gate = depletion_gate(y, schema, params, ASSIMILABLE_SPECS)
        if gate <= 0.0:
            return d
        s_slice = schema.slice("S")
        # Residual sugar is the reducing-sugar substrate (summed over the vector), a SOFT gate: a
        # dry wine (S ≈ 0 at the aging segment) makes ~none, but the trace is physically real, so
        # this is not an isolability claim. ``<= 0`` also absorbs a solver undershoot (S < 0).
        s_total = max(float(y[s_slice].sum()), 0.0)
        if s_total <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_maillard_browning"], params["T_ref"])
        r_sugar = params["k_maillard_browning"] * f_t * s_total * gate  # g sugar/L/h into browning
        if r_sugar <= 0.0:
            return d

        # Size the melanoidin formed + the amino-acid draw so BOTH ledgers close (all nitrogen and
        # all carbon retained in the polymer — the N-retaining branch; the deaminating branch is
        # D-87). With ``R`` the drawn blend's mass C:N (decision D-100 — {arginine, generic}, split
        # by the nitrogen each holds, so ``R`` is what the draw will actually realise):
        #   nitrogen:  r_aa·n_aa                    = r_m·n_m
        #   carbon:    r_sugar·c(sugar) + r_m·n_m·R = r_m·c_m
        # =>  r_m = r_sugar·c(sugar) / (c_m − n_m·R).
        c_sugar = carbon_mass_fraction(sugar_species(schema)[0])  # wine's single hexose
        c_m = carbon_mass_fraction(_MAILLARD_MELANOIDIN_SPECIES)
        n_m = nitrogen_mass_fraction(_MAILLARD_MELANOIDIN_SPECIES)
        ratio = assimilable_carbon_per_nitrogen(y, schema)  # R [g C per g N drawn]
        # denom > 0 by construction: the C-rich melanoidin (C:N ≈ 8:1) far exceeds the assimilable
        # blend's C:N, which is bounded by its two members (arginine ≈ 1.29 .. glutamine ≈ 2.14)
        # whatever the pool composition — so ≈ 0.73..0.84·c_m, never near zero. (A negative denom
        # would silently create sugar, so this ordering is load-bearing; a metadata test pins it,
        # the class docstring's trap. D-100 widened the margin's worst case from 1.29 to 2.14 and
        # it still clears 8 by ~3.7×.)
        denom = c_m - n_m * ratio
        mel_rate = r_sugar * c_sugar / denom  # g maillard_melanoidin/L/h

        # WINE-ONLY (unlike sugar-only :class:`Caramelization`, which D-90 made medium-agnostic):
        # this route reads the amino-acid pools, untracked in beer (D-32), so it never runs on
        # 3-slot S. ``d[s_slice] = -r_sugar`` is therefore correct as written (wine's single slot).
        # Were beer amino acids ever tracked, this would need the D-90 vectorized apportionment
        # Caramelization now carries (apportion the draw across the S vector; per-component
        # carbon fractions).
        d[s_slice] = -r_sugar
        # Draw the retained nitrogen from {arginine, generic} (decision D-100). The carbon returned
        # is exactly ``mel_rate·n_m·R`` — the same ``R`` the melanoidin was sized with — so the
        # carbon ledger closes to machine precision by construction, not by coincidence.
        draw_assimilable_nitrogen(d, y, schema, mel_rate * n_m)
        d[schema.slice(_MAILLARD_MELANOIDIN_SPECIES)] = mel_rate
        # The shared A420 browning index rises with the melanoidin formed (off every ledger, the
        # D-74 optical-index slot). Maillard melanoidins brown more per mass than caramelan, so this
        # yield is its own (larger) parameter. Monotone (mel_rate ≥ 0), joining the D-74/D-88 route.
        d[schema.slice("A420")] = params["y_a420_per_maillard_melanoidin"] * mel_rate
        return d


class OakExtraction(Process):
    """Non-oxidative aging: oak extractives diffuse into the beverage toward a ceiling (D-77/D-86).

    The sixth aging Process (D-77), the **first non-oxidative** one. **OakExtraction itself draws no
    O₂** and takes no share of the shared ``o2`` budget — a pure diffusion process. As a finished
    wine **or beer** (barrel-beer oak, D-86) sits in oak (barrel or chips/staves), five **aroma**
    extractives diffuse in and rise toward
    a saturation ceiling: **whiskey lactone** (β-methyl-γ-octalactone, "coconut", the signature
    oak-lactone note, LIGHT-toast dominant), **vanillin** ("vanilla", MEDIUM-toast peak),
    **guaiacol** (a lignin-pyrolysis "smoky/toasty" phenol, HEAVY-toast dominant — the oak/toast
    note, *distinct* from the Brett 4-ethylguaiacol of D-55), **eugenol** ("clove/spice",
    HEAVY-toast) and **furaneol** (HDMF, "caramel/toffee", a thermal sugar-degradation furanone
    RISING with toast — the caramel note of toasted/charred oak and, via the D-93 spirit soak-back,
    of ex-bourbon barrels; decision D-94). These five move **no** pool the D-67 OAV lens already
    read, so D-77 added four new aroma pools and D-94 the fifth; those five are a **separate axis**,
    orthogonal to the browning/acetaldehyde/SO₂/Strecker competition. **``furaneol`` is off every
    ledger** (wood/spirit-derived, the ``iso_alpha`` treatment) — so it does *not* touch core ``S``
    or the on-ledger D-88 caramelization ``melanoidin``, and the two do not collide: ``melanoidin``
    is the caramelization *colour body* (on-ledger, ``A420``), ``furaneol`` the *volatile aroma* of
    the same browning chemistry (off-ledger, OAV). The genuinely deferred beat is caramel aroma from
    the *beverage's own* thermal caramelization (on-ledger — it would divert a sliver of sugar
    carbon out of the melanoidin park); this D-94 pool is oak/spirit-derived only (see D-94).

    **The ellagitannin bridge (D-78).** This Process *also* extracts a fifth pool, **ellagitannin**
    (oak's hydrolysable TASTE tannin, LIGHT-toast dominant / declining with toast — thermolabile),
    by the *identical* diffusion-to-a-ceiling form. But ellagitannin is **not** O₂-orthogonal: a
    *separate* Process, :class:`EllagitanninOxidation`, consumes it to scavenge dissolved O₂ (oak
    PROTECTION). So the ``ellagitannin`` pool has two Processes on it (extraction here, oxidation
    there — the ``o2`` two-Processes-one-pool precedent); *this* Process only pushes it up toward
    its ceiling. Ellagitannin is a taste (astringency), read out by ``analysis.astringency_series``,
    not the OAV odor lens.

    ``d(C_i)/dt = k_oak_extraction · f(T) · max(0, ceiling_i − C_i)`` per extractive ``i`` — a
    **first-order approach FROM BELOW** to a per-compound ceiling, the exact inverse of
    :class:`EsterHydrolysis`'s ``max(0, ester − isoamyl_acetate_eq)`` net decay toward a floor.
    ``f(T) =
    arrhenius_factor(T, E_a_oak_extraction, T_ref)`` is the *weak* warmer-extracts-faster factor
    (diffusion-limited, so ``E_a_oak_extraction`` is deliberately low — well below the reaction
    E_a's of the oxidative Processes). One **shared** ``k_oak_extraction`` across all four this beat
    (the ceilings carry the toast *profile*; per-compound rates are a documented refinement).

    **The ceiling is set at the dose, in a set-and-hold state slot (the ``cation_charge`` idiom).**
    Each ``ceiling_i`` lives in its own state slot that **no Process touches** — it is
    written *only* by the :func:`~fermentation.scenario.compile._verb_add_oak` verb, which computes
    ``ceiling_i = oak_gpl · oak_yield_<compound>_<toast>`` (the provenance-backed toast-specific
    yields in ``oak.yaml``) and holds it constant. So the *dose* (oak_gpl, toast) is a scenario
    choice — like every dosed input — while the *physics* (rate, activation energy, per-gram yields)
    is provenance-backed data. This Process reads the ceiling from state and rises ``C_i`` to it.

    **Off every ledger — the iso_alpha precedent, cleaner than the O₂ Processes.** The four
    extractives are **exogenous wood-derived** mass, tracked like the hop-derived ``iso_alpha``
    (D-64): their carbon comes from an *untracked* oak source, so booking them as a mass would need
    a wood carbon pool that does not exist. So — like ``iso_alpha``/``o2``/``A420`` — the extracted
    slots (and the ceiling slots) are off ``total_carbon``/``total_mass``/``total_nitrogen``, and
    this Process **moves nothing conserved**: it touches only the four extracted slots and, being a
    pure g/L transfer, needs **no** ``chemistry.py`` species registration (no molar-mass conversion
    in the RHS — unlike the O₂ Processes, which at least convert via ``M_O2``). ``d(C_i)/dt ≥ 0``
    always (monotone rise; ``C_i`` approaches but never exceeds its ceiling).

    **Isolable + gated on the ceiling (prime directive #3).** Medium-agnostic — wired into BOTH
    wine and beer (barrel-beer oak, D-86: the oak slots are carried by both schemas via
    ``core.media._oak_specs``, since oak extraction is a wood property, not a grape/pH one — unlike
    the wine-only :class:`SulfiteOxidation` / :class:`StreckerDegradation`); the ``"whiskey_lactone"
    not in schema`` guard makes it a hard no-op on any bare medium besides. Wired **disabled at the
    compile seam** (aging is post-ferment); ``begin_aging`` enables it with the other aging
    Processes. With **no** oak dosed every ``ceiling_i`` is 0, so — via the explicit ``ceiling_i ≤
    0`` guard — the contribution is byte-for-byte zero (the ``max(0, …)`` alone would not suffice:
    the floor here is **0**, so a solver undershoot ``C_i = −ε`` would give ``max(0, ε) > 0`` and
    fabricate extract; the guard blocks it, the o2≤0 idiom for a zero floor). So a ``begin_aging``
    run with no ``add_oak`` is byte-for-byte the case without oak — an aged beverage that never saw
    wood. Tier
    **speculative** (the extraction *form* — diffusion-limited approach to a ceiling, warmer-faster
    — is sourced; every magnitude, the yields especially, is an order-of-magnitude estimate).
    **Scope (v1):** ellagitannins are now modelled (D-78 — extracted here, O₂-scavenged by
    :class:`EllagitanninOxidation`); ``oak_gpl`` is the generalized oak-contact dose subsuming
    chips-g/L and barrel surface-to-volume ratio; barrel fill-number depletion is now modelled
    (D-91 — a reused barrel sets lower ceilings via ``oak_fill_retention``, at the ``add_oak`` dose;
    the *within-fill* finite-reservoir draw-down and per-compound retention remain deferred);
    whiskey lactone is a lumped cis+trans pool. A first-fill ex-**spirit** barrel's residual-spirit
    **ethanol** soak-back (the ~1% ABV a bourbon barrel donates) is also modelled — but as an
    ``add_oak`` ethanol *dose* (D-92), NOT by this Process: soak-back ethanol comes from the spirit,
    not the wood, and lands on the core ``E`` slot, so it is a discrete on-ledger dose, orthogonal
    to this off-ledger diffusion axis. The bourbon *aroma* soak-back (D-93) IS drawn by THIS
    Process: an ex-spirit ``add_oak`` dose BUMPS the vanillin/whiskey_lactone/guaiacol (D-93) and
    ``furaneol`` (caramel, D-94) ceilings by
    the residual spirit's congeners, and this Process leaches them in gradually on top of the wood
    diffusion — a
    ceiling bump being the only wood + spirit *additive* form (a bolus into the pool would be erased
    by the ``gap = ceiling − conc`` gate below). A gradual within-fill reservoir leach remains
    deferred (D-92).
    """

    name = "oak_extraction"
    tier = Tier.SPECULATIVE
    #: Writes only the five extracted-compound slots — the four aroma extractives (D-77) plus the
    #: ``ellagitannin`` taste extractive (D-78). The ceiling slots are read, never written (a
    #: set-and-hold constant the ``add_oak`` verb owns). Off every ledger (exogenous wood-derived
    #: mass, the iso_alpha precedent), so nothing conserved moves. (``ellagitannin`` is *also*
    #: consumed by the separate :class:`EllagitanninOxidation` O₂ sink — two Processes on one pool,
    #: the ``o2`` precedent — but this Process only *extracts* it.) ``furaneol`` (caramel/toffee,
    #: the caramel furanone, D-94) is the fifth aroma extractive, extracted identically (off every
    #: ledger, so no collision with the on-ledger D-88 caramelization ``melanoidin``).
    touches = ("whiskey_lactone", "vanillin", "guaiacol", "eugenol", "furaneol", "ellagitannin")
    #: ``k_oak_extraction``/``E_a_oak_extraction`` are this Process's own (oak.yaml, D-77); and
    #: ``T_ref`` is shared with every Arrhenius rate. The per-compound ceilings ride in *state* (by
    #: ``add_oak``), not params, so they are not in ``reads``. Tiers cap the four extracted pools'
    #: output tiers via parameter-tier propagation (D-1), flooring them at speculative.
    reads: tuple[str, ...] = ("k_oak_extraction", "E_a_oak_extraction", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # The oak extractives are appended to BOTH media (D-86, barrel-beer oak): a hard no-op on a
        # bare schema without them, belt-and-suspenders to the wiring.
        if "whiskey_lactone" not in schema:
            return d
        # Gate on STATE (the ceilings) BEFORE reading any oak param — so an un-oaked run (every
        # ceiling 0) is byte-for-byte inert even when oak.yaml is not loaded (the Strecker/Sulfite
        # substrate-gate-before-params discipline; an enabled-but-undosed Process mustn't KeyError).
        # The EXPLICIT ceiling ≤ 0 guard is load-bearing — the floor is 0 (unlike the ester
        # hydrolysis floor, which is > 0), so
        # ``max(0, ceiling − C)`` alone would let a solver undershoot C = −ε fabricate extract.
        active: list[tuple[str, float]] = []
        for compound, ceiling_name in _OAK_COMPOUND_CEILINGS:
            ceiling = float(y[schema.slice(ceiling_name)][0])
            if ceiling <= 0.0:  # no oak dosed for this compound ⇒ inert
                continue
            conc = float(y[schema.slice(compound)][0])
            gap = ceiling - conc  # remaining headroom below the saturation ceiling
            if gap <= 0.0:  # already at/above the ceiling ⇒ no further extraction (monotone rise)
                continue
            active.append((compound, gap))
        if not active:  # nothing to extract ⇒ return before touching oak params
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_oak_extraction"], params["T_ref"])
        k = params["k_oak_extraction"]
        for compound, gap in active:
            d[schema.slice(compound)] = k * f_t * gap  # first-order approach from below, off-ledger
        return d


class EllagitanninOxidation(Process):
    """Oxidative aging: dissolved O₂ oxidises oak ellagitannin → oak protects the beverage (D-78).

    The fourth **oxidative** sibling to claim a share of the shared ``o2`` budget (after
    :class:`OxidativeAcetaldehyde` D-71, :class:`PhenolicBrowning` D-74 and
    :class:`SulfiteOxidation` D-72), and the **bridge** from the D-77 oak extractive axis to the O₂
    sub-axis. Oak's hydrolysable tannin — the ``ellagitannin`` pool :class:`OakExtraction` fills by
    diffusion — is a potent **sacrificial antioxidant**: dissolved O₂ oxidises its galloyl/gallate
    groups, so the tannin **scavenges O₂**, diverting it from phenolic browning and ethanol
    oxidation. The emergent payoff — **the D-78 spine** — is oak **PROTECTION**: an oaked +
    oxygenated wine browns **less** (lower ``A420``) and accumulates **less** oxidative acetaldehyde
    than an un-oaked wine at the same O₂ dose. This is the :class:`SulfiteOxidation` "SO₂ protects
    until exhausted" threshold (D-72) with one difference: the buffer is **renewable**. While the
    wine is below the ceiling :class:`OakExtraction` re-supplies tannin as fast as this Process
    burns it, so oak buffers redox for **months-to-years** (an oaked+oxygenated wine's acetaldehyde
    may never climb) — contrast SO₂, whose finite pool is spent once. That renewability is
    physically correct (a barrel is a large tannin reservoir); the eventual **within-fill** wood
    exhaustion (the reservoir depleting as a single fill ages) is the deferred finite-reservoir
    refinement — distinct from D-91's *across-fill* ``fill_number``, which lowers a reused barrel's
    starting ceilings at dose time but does not draw them down dynamically here.

    ``d(o2)/dt = −r`` with ``r = k_ellagitannin_oxidation · f(T) · [o2] · [ellagitannin]`` —
    **bilinear** in the dissolved-O₂ pool and the ellagitannin driver (the :class:`SulfiteOxidation`
    form), ``f(T) = arrhenius_factor(T, E_a_ellagitannin_oxidation, T_ref)`` the sourced
    warmer-scavenges-faster factor at **reaction** scale (its own param, ~50 kJ/mol — distinct from
    the *weak* diffusion ``E_a_oak_extraction`` that governs the tannin's extraction: scavenging is
    a chemical oxidation, extraction is diffusion). The tannin is **consumed** as it scavenges::

        d(ellagitannin)/dt = −y_ellag_per_o2 · r

    at a **mass-based** yield (``g ellagitannin per g O₂``), *not* a molar stoichiometry:
    ellagitannin is a lumped hydrolysable-tannin macromolecule with no clean molar mass, so an
    ``M_ellagitannin`` would be fake precision (contrast :data:`_SO2_PER_O2`, a real-molecule molar
    ratio). So this consumption **softens the astringency** the ``ellagitannin`` pool carries — but
    one contributor among several now built: the grape-tannin softeners
    :class:`TanninAnthocyaninCondensation` (D-79), :class:`TanninSelfPolymerization` (D-84) and
    :class:`TanninEthylTanninCondensation` (D-85) draw the *grape* ``tannin`` pool while this D-78
    sink draws the *oak* ``ellagitannin`` pool — both feed
    :func:`~fermentation.analysis.astringency_series` (which sums the two harsh tannins).

    **Substrate-gated ⇒ adds on top, NO re-baseline (the D-72/D-75 rule).** The O₂ draw is bilinear
    in ``[ellagitannin]``, which is zero unless oak is dosed (``add_oak``), so — exactly like
    :class:`SulfiteOxidation` (gated on SO₂) and :class:`StreckerDegradation` (gated on amino acids)
    — this sink is **zero without its substrate** and therefore **adds on top** of the shared O₂
    budget with **no re-baseline**: ``k_ethanol_oxidation + k_browning = 5.0e-4`` (the always-on
    anchor) is **untouched**, and the no-oak trajectory (either medium) is byte-for-byte kept. A
    nice illustration that the substrate-gated / always-on distinction — not the magnitude — is
    what's load-bearing: ``k_ellagitannin_oxidation`` is banded so that, when oak *is* present, this
    is a **major** sink (it takes roughly a third-to-half of the O₂), yet it still needs no
    re-baseline (unlike the always-on :class:`PhenolicBrowning`, which forced the D-74 re-baseline).
    It is banded so the protection is **partial** — an oaked beverage still shows *some* oxidative
    character.

    **Off every ledger, no conservation term (the :class:`SulfiteOxidation` precedent).** Both
    ``o2`` (D-71) and ``ellagitannin`` (wood-derived, off
    ``total_carbon``/``total_mass``/``total_nitrogen`` like ``iso_alpha``/``A420``, D-77) are
    unweighted, so oxidising the tannin to untracked products moves **nothing conserved** — this
    Process touches only those two slots and asserts nothing. This is why the mass-based yield is
    legitimate: no ledger reads the ``ellagitannin`` mass, so the lump carries no fabricated carbon.

    **Medium-agnostic + isolable + doubly substrate-gated (prime directive #3).** The
    ``ellagitannin`` slot and the ``o2`` pool are both carried by BOTH media (barrel-beer oak, D-86:
    the oak axis is a wood property, so — like :class:`OakExtraction`, and unlike the wine-only
    :class:`SulfiteOxidation` / :class:`StreckerDegradation` — it is wired into wine *and* beer).
    The ``"ellagitannin" not in schema`` guard makes it a hard no-op on a bare medium. Wired
    **disabled at the compile seam** (aging is post-ferment); ``begin_aging`` enables it with the
    other aging Processes. With no O₂ *or* no oak dosed the ``o2 ≤ 0`` / ``ellagitannin ≤ 0`` guards
    return byte-for-byte zero, so a reductive (no ``add_oxygen``) or an un-oaked aging is exactly
    the case without this Process. Tier **speculative** (the aging axis is the Tier-3 frontier; the
    *form* — O₂-limited, tannin-driven, warmer-faster — is sourced, the rate/yield magnitudes
    order-of-magnitude estimates).
    """

    name = "ellagitannin_oxidation"
    tier = Tier.SPECULATIVE
    #: Consumes its share of the dissolved-O₂ substrate and the ``ellagitannin`` tannin it oxidises
    #: — both slots off every ledger, so nothing conserved moves; it touches those two and nothing
    #: else. (``ellagitannin`` is also *extracted* by :class:`OakExtraction` — two Processes on one
    #: pool, the ``o2`` precedent.)
    touches = ("o2", "ellagitannin")
    #: ``k_ellagitannin_oxidation``/``E_a_ellagitannin_oxidation``/``y_ellag_per_o2`` are this
    #: Process's own (oak.yaml, D-78); ``T_ref`` is shared with every Arrhenius rate. Their tiers
    #: cap the ``o2``/``ellagitannin`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_ellagitannin_oxidation",
        "E_a_ellagitannin_oxidation",
        "y_ellag_per_o2",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # ellagitannin is appended to BOTH media (D-86, barrel-beer oak): a hard no-op on any bare
        # schema without it, belt-and-suspenders to the wiring.
        if "ellagitannin" not in schema or "o2" not in schema:
            return d
        o2 = float(y[schema.slice("o2")][0])
        ellag = float(y[schema.slice("ellagitannin")][0])
        # No oxidant OR no tannin ⇒ no scavenging: reductive/un-oaked aging is byte-for-byte the
        # case without this Process. Gate on the ellagitannin STATE before reading any oak param
        # (the OakExtraction/Strecker substrate-gate-before-params discipline — an
        # enabled-but-undosed Process mustn't KeyError if oak.yaml is absent). ``<= 0`` also absorbs
        # solver undershoot.
        if o2 <= 0.0 or ellag <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_ellagitannin_oxidation"], params["T_ref"])
        # This route's SHARE of the O₂-depletion rate (bilinear in o2 and the tannin driver, the
        # SulfiteOxidation form) — substrate-gated on ellagitannin, so it adds on top of the anchor
        # with no re-baseline (D-72/D-75). ProcessSet sums the sinks, so o2 splits by kᵢ/Σk.
        r_o2 = (
            params["k_ellagitannin_oxidation"] * f_t * o2 * ellag
        )  # g O2/L/h via the tannin route
        d[schema.slice("o2")] = -r_o2
        # The sacrificial tannin is consumed at a MASS-based yield (g ellag / g O2), not a molar
        # stoichiometry — ellagitannin is a lumped macromolecule with no clean molar mass. Both
        # slots off every ledger, so this consumption moves nothing conserved (the SulfiteOxidation
        # idiom).
        d[schema.slice("ellagitannin")] = -params["y_ellag_per_o2"] * r_o2
        return d


class TanninAnthocyaninCondensation(Process):
    """Aging condensation: grape anthocyanin + grape tannin → stable polymeric pigment (D-79).

    The eighth aging Process, the **second non-oxidative** one (after :class:`OakExtraction`), and a
    **third separate axis**: the DOMINANT real red-wine astringency-softening and colour-stabilizing
    mechanism the D-77/D-78 oak beats deferred and the milestone-3 plan names. As a finished red
    wine
    ages, free monomeric **anthocyanin** (the bright, bleachable purple-red pigment) and condensed
    grape **tannin** (the harsh young flavan-3-ol astringency) combine into **polymeric pigment** —
    a
    softer-tasting, SO₂/pH-**stable** colour form (Somers 1971; Ribéreau-Gayon). Two emergent
    payoffs
    the D-79 readouts expose: **astringency softens**
    (:func:`~fermentation.analysis.astringency_series`
    now reads free tannin — ``tannin + ellagitannin``, both harsh; the polymeric pigment is soft and
    excluded — so it declines as tannin condenses) and **colour stabilizes**
    (:func:`~fermentation.analysis.color_series` counts free anthocyanin *and* the polymeric
    pigment,
    so total red colour is retained as its form shifts labile → stable — the young purple → aged
    brick-red evolution).

    ``r = k_polymerization · f(T) · [anthocyanin] · [tannin]`` — **bilinear** in the two grape pools
    (the :class:`SulfiteOxidation`/:class:`EllagitanninOxidation` form), ``f(T) =
    arrhenius_factor(T, E_a_polymerization, T_ref)`` the sourced warmer-condenses-faster factor at
    reaction scale. It consumes both reactants::

        d(anthocyanin)/dt = −r
        d(tannin)/dt      = −y_tannin_per_anthocyanin · r

    at a **mass-based** yield (``g tannin per g anthocyanin``), *not* a molar stoichiometry — both
    the anthocyanin pool and the condensed-tannin macromolecule are lumped pools with no clean molar
    mass, so an ``M_tannin`` would be fake precision (the ``y_ellag_per_o2`` / D-78 precedent).

    **The polymeric pigment is now an integrated ``polymeric_pigment`` slot (D-81 promotion).**
    Through D-79/D-80 condensation was anthocyanin's **sole** fate, so the stable pigment was
    exactly ``anthocyanin₀ − anthocyanin(t)`` and was reconstructed post-hoc (the ``iso_alpha``/IBU
    readout pattern). D-81's :class:`AnthocyaninFading` adds a **second** anthocyanin fate
    (oxidative → a colourless ``faded_anthocyanin`` slot), so that reconstruction no longer isolates
    the
    pigment — it would wrongly count the faded fraction — and the pigment had to become a real slot
    (the A420 discriminator, D-74: a driver with competing sinks is not reconstructible). So this
    Process now also writes ``d(polymeric_pigment)/dt = +r`` (anthocyanin-equivalents); the
    bridged route (D-80) writes the same. Both routes feed one shared pigment pool, so
    :func:`~fermentation.analysis.polymeric_pigment_series` reads the slot directly and
    :func:`~fermentation.analysis.color_series` (free anthocyanin + pigment) now **genuinely
    declines** as fading removes free anthocyanin without adding pigment (the colour-stability
    payoff D-79 predicted this promotion would unmask).

    **Oak-independent AND O₂-independent (the D-79 correctness crux).** Tannin–anthocyanin
    polymerization is a **grape**-tannin + **grape**-anthocyanin reaction: a steel-tank red with no
    oak and no oxygen still polymerizes, softens, and stabilizes its colour. So this Process gates
    on
    the grape-derived ``anthocyanin`` and ``tannin`` pools **only** — it draws **no** share of the
    shared ``o2`` budget (unlike every D-71..D-78 oxidative sink) and reads **no** oak pool. In
    particular ``tannin`` is the grape **condensed** (flavan-3-ol) tannin — a *different* molecule
    from oak's hydrolysable ``ellagitannin`` (D-78); reusing ellagitannin would wrongly make
    polymerization impossible without an ``add_oak`` dose. This is the grape ``tannin`` pool the
    D-78
    note deliberately left the namespace free for. Because it touches no ``o2`` it does not even
    interact with the ``k_ethanol_oxidation + k_browning`` always-on anchor — a wholly separate,
    non-oxidative axis (the :class:`OakExtraction` diffusion-axis precedent, on grape pools).

    **Off every ledger, no conservation term (the ``iso_alpha``/``ellagitannin`` precedent).** Both
    ``anthocyanin`` and ``tannin`` are grape-derived — their carbon comes from an *untracked*
    grape-solids source — so they are off ``total_carbon``/``total_mass``/``total_nitrogen`` (like
    ``iso_alpha``/``A420``/``ellagitannin``), and this Process moves **nothing conserved**: it
    touches only those two slots and asserts nothing (a pure g/L transfer, so — like
    :class:`OakExtraction` — no ``chemistry.py`` species registration). This is also why the yield
    is
    mass-based: no ledger reads either lump's mass, so it carries no fabricated carbon.

    **Doubly substrate-gated ⇒ adds on top, NO re-baseline (the D-72/D-75 rule).** The rate is
    bilinear in both pools, each zero unless dosed as a grape must input (``anthocyanin_gpl`` /
    ``tannin_gpl``), so a **white** wine (no anthocyanin) or a no-tannin run is byte-for-byte the
    case
    without this Process. Wine-only (both slots are appended to ``wine_schema``), like
    :class:`OakExtraction` / :class:`EllagitanninOxidation`; the ``"anthocyanin" not in schema``
    guard
    makes it a hard no-op besides. Wired **disabled at the compile seam** (aging is post-ferment);
    ``begin_aging`` enables it with the other aging Processes. Tier **speculative** (the aging axis
    is
    the Tier-3 frontier; the condensation *form* — bimolecular, saturating on the limiting
    anthocyanin, warmer-faster, softening + colour-stabilizing — is sourced, the magnitudes
    order-of-magnitude estimates). This is the **direct** condensation route; the sibling softeners
    are all now built: the acetaldehyde-bridged (ethylidene) anthocyanin route
    (:class:`AcetaldehydeBridgedCondensation`, D-80), SO₂/pH anthocyanin bleaching
    (:class:`AnthocyaninFading`, D-81 + :class:`ThermalAnthocyaninFade`, D-83), and grape-tannin
    self-polymerization (:class:`TanninSelfPolymerization`, D-84 direct +
    :class:`TanninEthylTanninCondensation`, D-85 bridged) — so astringency softening is **no longer
    one-directional**. See ``polymerization.yaml`` for the full scope + provenance.
    """

    name = "tannin_anthocyanin_condensation"
    tier = Tier.SPECULATIVE
    #: Consumes the two grape pools it condenses and DEPOSITS the stable ``polymeric_pigment`` (D-81
    #: promotion) — all three off every ledger (grape-derived colour-equivalents, the
    #: ``ellagitannin`` precedent), so nothing conserved moves. The pigment is now an integrated
    #: slot, not the old
    #: post-hoc readout: D-81's :class:`AnthocyaninFading` gives anthocyanin a second fate, so
    #: ``anthocyanin₀ − anthocyanin`` no longer isolates the pigment (the A420 discriminator, D-74).
    touches = ("anthocyanin", "tannin", "polymeric_pigment")
    #: ``k_polymerization``/``E_a_polymerization``/``y_tannin_per_anthocyanin`` are this Process's
    #: own (polymerization.yaml, D-79); ``T_ref`` is shared with every Arrhenius rate. Their tiers
    #: cap the ``anthocyanin``/``tannin`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_polymerization",
        "E_a_polymerization",
        "y_tannin_per_anthocyanin",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (anthocyanin/tannin are appended to wine_schema): a hard no-op on any
        # schema without them, belt-and-suspenders to the wine-only wiring.
        if "anthocyanin" not in schema or "tannin" not in schema:
            return d
        anthocyanin = float(y[schema.slice("anthocyanin")][0])
        tannin = float(y[schema.slice("tannin")][0])
        # No anthocyanin OR no tannin ⇒ no condensation: a white wine (no anthocyanin) or a
        # no-tannin
        # run is byte-for-byte the case without this Process. Gate on the grape STATE before reading
        # any polymerization param (the OakExtraction/Strecker substrate-gate-before-params
        # discipline — an enabled-but-undosed Process mustn't KeyError if polymerization.yaml is
        # absent). ``<= 0`` also absorbs solver undershoot.
        if anthocyanin <= 0.0 or tannin <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_polymerization"], params["T_ref"])
        # Bilinear direct condensation (the SulfiteOxidation/EllagitanninOxidation form): the rate
        # is
        # the anthocyanin consumption (g anthocyanin/L/h). No o2 term — a wholly separate, non-
        # oxidative grape axis (oak- AND O₂-independent), so it never touches the O₂ anchor.
        r = params["k_polymerization"] * f_t * anthocyanin * tannin  # g anthocyanin/L/h condensed
        d[schema.slice("anthocyanin")] = -r
        # Tannin consumed per anthocyanin condensed at a MASS-based yield (g tannin / g
        # anthocyanin),
        # NOT a molar stoichiometry — both are lumped pools with no clean molar mass (the
        # y_ellag_per_o2 idiom). Both slots off every ledger, so this moves nothing conserved.
        d[schema.slice("tannin")] = -params["y_tannin_per_anthocyanin"] * r
        # The condensed anthocyanin is DEPOSITED into the stable polymeric_pigment slot in
        # anthocyanin-equivalents (D-81 promotion): every unit anthocyanin loses to condensation
        # enters the pigment, so d(polymeric_pigment)/dt = +r exactly balances d(anthocyanin)/dt.
        # Off every ledger (grape-derived colour-equivalent), so this deposit moves nothing
        # conserved.
        d[schema.slice("polymeric_pigment")] = r
        return d


class AcetaldehydeBridgedCondensation(Process):
    """Aging condensation: acetaldehyde bridges grape tannin + anthocyanin → pigment (D-80).

    The ninth aging Process, the **third non-oxidative** one (after :class:`OakExtraction` and
    :class:`TanninAnthocyaninCondensation`), and the **split-ledger** beat D-79 explicitly deferred.
    It is the *second* formation pathway to polymeric pigment (after the D-79 *direct* condensation)
    and the **first link from the oxidative sub-axis to red-wine colour**: dissolved-O₂ acetaldehyde
    (:class:`OxidativeAcetaldehyde`, D-71) forms an **ethylidene bridge** ``—CH(CH₃)—`` linking a
    grape tannin unit to an anthocyanin unit (tannin–ethyl–anthocyanin), an
    *acetaldehyde-accelerated*
    condensation that stabilizes colour and softens astringency. So dosing O₂ (``add_oxygen``,
    micro-oxygenation) now **stabilizes red colour** — the emergent "controlled micro-ox" payoff
    D-79
    named — while simultaneously **drawing acetaldehyde down**, lowering the oxidised/bruised-apple
    acetaldehyde note (the real winemaking benefit of micro-ox: colour *and* reduced harshness).

    ``r = k_acetaldehyde_bridge · f(T) · [acetaldehyde_free] · [anthocyanin] · [tannin]`` — a
    **trilinear** lumped termolecular step (the D-79 bilinear form plus the acetaldehyde factor),
    with
    ``f(T) = arrhenius_factor(T, E_a_acetaldehyde_bridge, T_ref)`` the sourced warmer-bridges-faster
    factor at reaction scale (its **own** E_a, distinct from the direct-route ``E_a_polymerization``
    —
    two reactions, two E_a, prime directive #2, the :class:`EllagitanninOxidation` vs
    :class:`OakExtraction` precedent). Anchored on **anthocyanin** consumption (like the direct
    route),
    so ``r`` is the anthocyanin condensation rate and the pigment reconstruction stays valid
    (below).
    It consumes all three reactants::

        d(anthocyanin)/dt  = −r                                       # off-ledger (grape-derived)
        d(tannin)/dt       = −y_tannin_per_anthocyanin · r            # off-ledger (reuse D-79
        yield)
        d(acetaldehyde)/dt = −y_acetaldehyde_per_anthocyanin · r      # ON the carbon ledger
        d(ethyl_bridge)/dt = +(acetaldehyde carbon consumed) / c(ethylidene)   # ON-ledger, exact

    **The split ledger (the D-80 crux, and why D-79 deferred it).** One reaction straddles *two*
    ledgers. The grape-phenolic bulk (``anthocyanin`` + ``tannin``) is **off** every ledger
    (grape-derived, untracked — the ``iso_alpha``/``ellagitannin`` precedent), so consuming it moves
    nothing conserved, exactly as the direct route (D-79). But acetaldehyde's carbon is **on** the
    carbon ledger — it was borrowed carbon-exactly from ethanol ``E`` by
    :class:`OxidativeAcetaldehyde`
    (D-71) and by the D-27 fermentative buffer. Consuming on-ledger acetaldehyde into an
    *off*-ledger
    pigment would make carbon **vanish** and fail :func:`~fermentation.validation.conservation.\
    assert_conserved` — the trap D-79 named. The fix is a **new on-ledger ``ethyl_bridge`` state
    slot**
    (wine-only; weighted at ``carbon_mass_fraction("ethylidene")`` in
    :func:`~fermentation.validation.conservation.total_carbon`) that captures exactly the
    acetaldehyde
    carbon. The transfer uses the :class:`EsterHydrolysis` carbon-exact split — release at
    ``c(acetaldehyde)``, re-deposit at ``c(ethylidene)`` — so ``total_carbon`` closes to **machine
    precision** (the acetaldehyde C leaving equals the bridge C arriving), *non-trivially*: unlike
    the
    direct route (where carbon is flat because nothing on the ledger moves), here acetaldehyde↓ and
    ethyl_bridge↑ exactly cancel. Acetaldehyde (C2H4O) loses only its **carbonyl oxygen as water**
    on
    bridging, leaving the two-carbon ethylidene (C2H4); that lost O is the standing aging-axis
    **mass**
    gap (``total_mass`` weights only ``{S, E, CO2}`` and is never asserted on an aging run — the
    D-71
    ``E → acetaldehyde`` scope-out). Carbon is the invariant.

    **Why the bridge is an integrated slot, not a post-hoc readout (the A420 discriminator, D-74).**
    Acetaldehyde has **competing** fates — production (fermentative D-27, oxidative D-71), reduction
    back to ``E`` (D-27), SO₂ binding (D-47), and now bridging — so the bridged amount is **not**
    reconstructible from the acetaldehyde pool's drawdown. (As of D-81 the **pigment** is a slot for
    the same reason: :class:`AnthocyaninFading` gives anthocyanin a competing colourless fate,
    so ``anthocyanin₀ − anthocyanin`` no longer isolates the pigment either — both drivers
    acquired competing sinks, so both became slots.) And structurally ``total_carbon`` =
    ``weights @ y`` reads
    *state*, so the captured carbon must physically live in a slot for closure. Both reasons force
    the integrated ``ethyl_bridge`` slot.

    **Reads FREE acetaldehyde, not total (the D-47 precedent).** SO₂-bound acetaldehyde is the
    bisulfite adduct — its carbonyl is blocked, so it **cannot** form the ethylidene bridge (there
    is
    no free carbonyl for the flavanol to attack), exactly as :class:`AcetaldehydeReduction` reduces
    only the *free* share under SO₂ (D-47). So when SO₂ is dosed (``so2_total > 0``) the rate reads
    :func:`~fermentation.core.acidbase.free_acetaldehyde` rather than the total; the guard is exact
    —
    an unsulfited run pays no per-RHS pH ``brentq`` and is byte-for-byte the total-acetaldehyde
    case.
    **Emergent payoff:** SO₂ **delays** acetaldehyde-mediated colour stabilization (bound
    acetaldehyde
    is unavailable to bridge), a real winemaking fact that falls out of the shared binding
    equilibrium
    with nothing scripted — the flip side of D-72's "SO₂ protects against oxidation".

    **Both formation routes feed one shared ``polymeric_pigment`` slot (D-81).** By anchoring on
    anthocyanin (tannin–ethyl–anthocyanin, *not* tannin–ethyl–tannin) this route deposits pigment in
    anthocyanin-equivalents (``+r``) into the same slot the direct route fills, so
    :func:`~fermentation.analysis.polymeric_pigment_series` reads the total pigment straight from
    the slot and :func:`~fermentation.analysis.astringency_series` softens *more* (this route draws
    ``tannin`` down). :func:`~fermentation.analysis.color_series` (free anthocyanin + pigment) now
    **genuinely declines** once :class:`AnthocyaninFading` (D-81) removes free anthocyanin to the
    colourless ``faded_anthocyanin`` sink without adding pigment — the emergent micro-ox tension the
    O₂-coupled fade creates (some anthocyanin bridges to *stable* pigment here, some fades to
    *colourless*; SO₂ both protects against the fade, via the D-72 o2 draw, and delays this
    bridging). The tannin–ethyl–tannin branch — this route's *self*-reaction analogue, bridging two
    flavanols with no anthocyanin — is now built (:class:`TanninEthylTanninCondensation`, D-85,
    alongside the D-84 direct self-polymerization); it draws ``tannin`` without touching
    ``anthocyanin`` and deposits *no* pigment (colourless), so it softens astringency without
    perturbing this route's anchor-on-anthocyanin colour accounting.

    **Triply substrate-gated + wine-only + isolable (prime directive #3).** The rate is trilinear in
    ``acetaldehyde`` × ``anthocyanin`` × ``tannin``, so a white / no-tannin / no-acetaldehyde wine
    (and
    all of beer) is byte-for-byte the case without this Process. Wine-only (all four slots are
    wine-only: ``acetaldehyde`` is medium-agnostic but ``anthocyanin``/``tannin``/``ethyl_bridge``
    are
    appended to ``wine_schema``), gated on the grape STATE before reading any polymerization param
    (an enabled-but-undosed Process must not KeyError if ``polymerization.yaml`` is absent). Wired
    **disabled at the compile seam** (aging is post-ferment); ``begin_aging`` enables it with the
    other
    aging Processes. Tier **speculative** (the aging axis is the Tier-3 frontier; the *form* —
    acetaldehyde-bridged, trilinear, warmer-faster, SO₂-blocked, colour-stabilizing — is sourced,
    the
    magnitudes order-of-magnitude estimates). **Scope (v1):** the tannin–ethyl–tannin analogue is
    built separately (D-85, see above);
    the bridged pigment's colour is counted equal to free anthocyanin (the D-79 equal-absorptivity
    simplification). See ``polymerization.yaml`` for the full scope + provenance.
    """

    name = "acetaldehyde_bridged_condensation"
    tier = Tier.SPECULATIVE
    #: Consumes the two grape pools it bridges (``anthocyanin``/``tannin`` — off every ledger, the
    #: D-79 precedent) plus the on-ledger ``acetaldehyde``, whose carbon it captures in the
    #: on-ledger ``ethyl_bridge`` slot (the split-ledger transfer), and DEPOSITS the stable
    #: ``polymeric_pigment`` (D-81 promotion — the second formation route into the shared pigment
    #: pool). Touches those five and nothing else — the ``acetaldehyde → ethyl_bridge`` carbon
    #: closes exactly; the grape pair + pigment (all off-ledger) move nothing conserved.
    touches = ("acetaldehyde", "ethyl_bridge", "anthocyanin", "tannin", "polymeric_pigment")
    #: ``k_acetaldehyde_bridge``/``E_a_acetaldehyde_bridge``/``y_acetaldehyde_per_anthocyanin`` are
    #: this Process's own (polymerization.yaml, D-80); ``y_tannin_per_anthocyanin`` is shared with
    #: the
    #: direct route (D-79, same lumped adduct stoichiometry); ``T_ref`` is shared with every
    #: Arrhenius
    #: rate. The SO₂/pH params read inside :func:`free_acetaldehyde`/:func:`ph_of_state` are omitted
    #: —
    #: all plausible and the Process is already speculative, so they add no tier headline (the D-47
    #: :class:`AcetaldehydeReduction` / MLF-gate precedent). Tiers cap the outputs via D-1.
    reads: tuple[str, ...] = (
        "k_acetaldehyde_bridge",
        "E_a_acetaldehyde_bridge",
        "y_acetaldehyde_per_anthocyanin",
        "y_tannin_per_anthocyanin",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (anthocyanin/tannin/ethyl_bridge are appended to wine_schema): a hard
        # no-op
        # on any schema without them, belt-and-suspenders to the wine-only wiring.
        if "anthocyanin" not in schema or "tannin" not in schema:
            return d
        anthocyanin = float(y[schema.slice("anthocyanin")][0])
        tannin = float(y[schema.slice("tannin")][0])
        acetaldehyde = float(y[schema.slice("acetaldehyde")][0])
        # Triply substrate-gated: no acetaldehyde OR no anthocyanin OR no tannin ⇒ no bridging, so a
        # white / no-tannin / no-acetaldehyde wine is byte-for-byte the case without this Process.
        # Gate on the grape STATE before reading any polymerization param (the
        # OakExtraction/Strecker
        # substrate-gate-before-params discipline — an enabled-but-undosed Process mustn't KeyError
        # if
        # polymerization.yaml is absent). ``<= 0`` also absorbs solver undershoot.
        if anthocyanin <= 0.0 or tannin <= 0.0 or acetaldehyde <= 0.0:
            return d
        # SO₂-bound acetaldehyde is the bisulfite adduct — its carbonyl is blocked, so it CANNOT
        # form
        # the ethylidene bridge; read only the FREE share under SO₂ (the D-47 AcetaldehydeReduction
        # precedent). The ``so2_total > 0`` guard is EXACT — an unsulfited run pays no per-RHS pH
        # ``brentq`` and its rate is byte-for-byte the total-acetaldehyde case. Emergent: SO₂ delays
        # acetaldehyde-mediated colour stabilization (bound acetaldehyde is unavailable to bridge).
        bridging_acetaldehyde = acetaldehyde
        if SO2_STATE_KEY in schema and float(y[schema.slice(SO2_STATE_KEY)][0]) > 0.0:
            bridging_acetaldehyde = free_acetaldehyde(
                y, schema, params, ph_of_state(y, schema, params)
            )
            if bridging_acetaldehyde <= 0.0:  # all acetaldehyde bound ⇒ none available to bridge
                return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_acetaldehyde_bridge"], params["T_ref"])
        # Trilinear lumped termolecular rate (the D-79 bilinear form + the free-acetaldehyde
        # factor):
        # r is the anthocyanin condensation rate (g anthocyanin/L/h), anchored on anthocyanin so its
        # sole fate stays "→ pigment" and polymeric_pigment_series remains reconstructible.
        r = params["k_acetaldehyde_bridge"] * f_t * bridging_acetaldehyde * anthocyanin * tannin
        d[schema.slice("anthocyanin")] = -r
        # Tannin consumed per anthocyanin bridged at the SAME mass-based yield as the direct route
        # (same lumped adduct stoichiometry, D-79); off every ledger, moves nothing conserved.
        d[schema.slice("tannin")] = -params["y_tannin_per_anthocyanin"] * r
        # Deposit the bridged anthocyanin into the SHARED polymeric_pigment slot in
        # anthocyanin-equivalents (D-81 promotion), exactly as the direct route: both formation
        # pathways feed one pigment pool, so d(polymeric_pigment)/dt = +r. Off-ledger colour-
        # equivalent (distinct from the on-ledger ethyl_bridge carbon booked below — no
        # double-count: pigment tracks colour, ethyl_bridge the acetaldehyde carbon locked into it).
        d[schema.slice("polymeric_pigment")] = r
        # The ON-ledger half — the split-ledger carbon capture. Acetaldehyde consumed by the bridge
        # (g/L/h); its carbon must NOT vanish into the off-ledger pigment, so re-deposit it into the
        # on-ledger ethyl_bridge pool via the EsterHydrolysis carbon-exact split (release at
        # acetaldehyde's fraction, re-deposit at ethylidene's — acetaldehyde loses only its carbonyl
        # O
        # as water). total_carbon closes to machine precision: acetaldehyde C leaving == bridge C
        # arriving. The yield only sets the MAGNITUDE of the acetaldehyde drawdown; carbon balances
        # for any value.
        acet_consumed = params["y_acetaldehyde_per_anthocyanin"] * r  # g acetaldehyde/L/h
        d[schema.slice("acetaldehyde")] = -acet_consumed
        carbon_released = acet_consumed * carbon_mass_fraction(
            _BRIDGE_ACETALDEHYDE_SPECIES
        )  # g C/L/h
        d[schema.slice("ethyl_bridge")] = carbon_released / carbon_mass_fraction(
            _ETHYL_BRIDGE_SPECIES
        )
        return d


class AnthocyaninFading(Process):
    """Aging oxidative fade: dissolved O₂ degrades free anthocyanin → colourless (decision D-81).

    The tenth aging Process and the beat that finally makes :func:`~fermentation.analysis.\
    color_series` **genuinely decline**. Free monomeric ``anthocyanin`` has a **second**,
    irreversible fate besides condensation into stable pigment (D-79/D-80): dissolved O₂ drives it
    to **colourless** degradation products (quinone/peroxide-mediated oxidative fading; Somers &
    Evans; Ribéreau-Gayon). This is the *bleaching loss* the D-79/D-80 colour axis deferred — and
    the reason the pigment had to be promoted to an integrated slot (the ``anthocyanin₀ −
    anthocyanin`` reconstruction can no longer isolate the pigment once anthocyanin has two fates).

    ``r_o2 = k_anthocyanin_fade · f(T) · [o2] · [anthocyanin]`` — **bilinear** in the shared O₂ pool
    and the anthocyanin driver (the :class:`SulfiteOxidation`/:class:`EllagitanninOxidation` form),
    ``f(T) = arrhenius_factor(T, E_a_anthocyanin_fade, T_ref)`` the sourced warmer-fades-faster
    factor. It is this route's **share** of the O₂-depletion rate; anthocyanin is transferred to the
    colourless ``faded_anthocyanin`` slot at a mass-based yield::

        d(o2)/dt                = −r_o2
        d(anthocyanin)/dt       = −y_anthocyanin_per_o2 · r_o2
        d(faded_anthocyanin)/dt = +y_anthocyanin_per_o2 · r_o2

    a pure g/L **transfer** between two off-ledger slots (the faded pool gains exactly what
    anthocyanin loses), so the D-81 colour identity ``anthocyanin + polymeric_pigment +
    faded_anthocyanin ≡ anthocyanin₀`` closes by construction. ``faded_anthocyanin`` is colourless,
    so :func:`~fermentation.analysis.color_series` (free anthocyanin + pigment) falls by exactly the
    faded amount while the condensed **polymeric pigment survives** — the young-bleachable →
    aged-stable colour dynamic, the colour-stability payoff.

    **O₂-COUPLED, so SO₂ protection is EMERGENT (the D-81 correctness crux).** Fading draws the
    **shared** ``o2`` budget, exactly like the D-71..D-78 oxidative sinks — it is *not* a scripted
    ``g(SO₂, pH)`` decay. That matters: "SO₂ protects the colour" is true because SO₂ is an
    **antioxidant** — it scavenges O₂ (bisulfite oxidation, :class:`SulfiteOxidation`, D-72), so a
    sulfited wine simply has **less O₂ left** to fade the anthocyanin. ``ProcessSet`` sums the O₂
    sinks and splits ``o2`` by kᵢ/Σk, so this protection **falls out** of the shared pool with
    nothing scripted (the D-72/D-80 "SO₂ effect, emergent" signature). It also creates the real
    micro-ox tension: under O₂, some anthocyanin bridges to *stable* pigment
    (:class:`AcetaldehydeBridgedCondensation`, D-80) while some **fades to colourless** here — and
    SO₂ both *protects* against this fade (via the D-72 draw) and *delays* the bridging (bound
    acetaldehyde can't bridge, D-80), all from one shared equilibrium.

    **Substrate-gated ⇒ adds on top, NO re-baseline (the D-72/D-75/D-78 rule).** The O₂ draw is
    bilinear in ``[anthocyanin]``, zero unless a red must (``anthocyanin_gpl``) is dosed, so — like
    :class:`EllagitanninOxidation` — this sink is zero without its substrate and adds on top of the
    shared O₂ budget with no re-baseline: the ``k_ethanol_oxidation + k_browning = 5.0e-4`` anchor
    is untouched, and every white / no-anthocyanin (and all-beer) trajectory is byte-for-byte
    preserved. A *red* wine dosed with both anthocyanin and O₂ does now split its O₂ one more way
    (a new real sink competing with the D-71..D-78 siblings) — the physically-correct cost of a new
    oxidative pathway, documented, not silent.

    **Off every ledger, no conservation term (the :class:`EllagitanninOxidation` precedent).**
    ``o2`` (D-71) is carbon-free and both ``anthocyanin`` and ``faded_anthocyanin`` are
    grape-derived colour pools off ``total_carbon``/``total_mass``/``total_nitrogen`` (``iso_alpha``
    precedent), so fading moves **nothing conserved** — it touches only those three slots and
    asserts nothing. The mass-based yield is therefore legitimate (no ledger reads either lump's
    mass), and
    the anthocyanin→faded transfer is exact regardless of the yield's value.

    **Wine-only + isolable + doubly substrate-gated (prime directive #3).** ``anthocyanin`` /
    ``faded_anthocyanin`` are wine-only (appended to ``wine_schema``), so this is wired into the
    *wine* medium only; the ``"anthocyanin" not in schema`` guard makes it a hard no-op besides.
    Wired **disabled at the compile seam** (aging is post-ferment); ``begin_aging`` enables it with
    the other aging Processes. With no O₂ *or* no anthocyanin the ``o2 ≤ 0`` / ``anthocyanin ≤ 0``
    guards return byte-for-byte zero, so a reductive (no ``add_oxygen``) or a white aging is exactly
    the case without this Process. Tier **speculative** (the aging axis is the Tier-3 frontier; the
    *form* — O₂-limited, anthocyanin-driven, warmer-faster, SO₂-protected-emergently — is sourced,
    the rate/yield magnitudes order-of-magnitude estimates).

    **Scope:** the **oxidative** fade only — this route fades free anthocyanin *via O₂*. Two related
    phenomena, once deferred, are now built as their own beats (named, not smuggled in): (1) the
    **reversible SO₂/pH masking** of monomeric anthocyanin (the flavylium ⇌ colourless bisulfite
    adduct / carbinol equilibrium — the literal Somers "bleaching" assay) is a fast equilibrium
    *readout*, not a fate — :func:`~fermentation.analysis.observed_color_series` (D-82); (2) the
    **O₂-independent** (thermal/hydrolytic) bottle-aging fade — :class:`ThermalAnthocyaninFade`
    (D-83), a *separate* Process that fills the same ``faded_anthocyanin`` sink with **no** O₂ — so
    an anaerobic sealed red now fades thermally too (retiring the earlier "anaerobic red holds its
    colour" scope note). See ``polymerization.yaml`` for the full scope + provenance.
    """

    name = "anthocyanin_fading"
    tier = Tier.SPECULATIVE
    #: Consumes its share of the dissolved-O₂ substrate and TRANSFERS free anthocyanin into the
    #: colourless ``faded_anthocyanin`` slot — all three off every ledger, so nothing conserved
    #: moves; it touches those three and nothing else. (``anthocyanin`` is also drawn by the two
    #: condensation routes and the D-83 thermal fade — four Processes on one pool, the ``o2``
    #: precedent.)
    touches = ("o2", "anthocyanin", "faded_anthocyanin")
    #: ``k_anthocyanin_fade``/``E_a_anthocyanin_fade``/``y_anthocyanin_per_o2`` are this Process's
    #: own (polymerization.yaml, D-81); ``T_ref`` is shared with every Arrhenius rate. Their tiers
    #: cap the ``anthocyanin``/``faded_anthocyanin`` output tiers via parameter-tier propagation
    #: (D-1).
    reads: tuple[str, ...] = (
        "k_anthocyanin_fade",
        "E_a_anthocyanin_fade",
        "y_anthocyanin_per_o2",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (anthocyanin/faded_anthocyanin are appended to wine_schema): a hard no-op
        # on any schema without them, belt-and-suspenders to the wine-only wiring.
        if "anthocyanin" not in schema or "o2" not in schema:
            return d
        o2 = float(y[schema.slice("o2")][0])
        anthocyanin = float(y[schema.slice("anthocyanin")][0])
        # No oxidant OR no anthocyanin ⇒ no fading: reductive/white aging is byte-for-byte the case
        # without this Process. Gate on the anthocyanin STATE before reading any fade param (the
        # OakExtraction/Strecker substrate-gate-before-params discipline — an enabled-but-undosed
        # Process mustn't KeyError if polymerization.yaml is absent). ``<= 0`` also absorbs solver
        # undershoot.
        if o2 <= 0.0 or anthocyanin <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_anthocyanin_fade"], params["T_ref"])
        # This route's SHARE of the O₂-depletion rate (bilinear in o2 and the anthocyanin driver,
        # the SulfiteOxidation/EllagitanninOxidation form) — substrate-gated on anthocyanin, so it
        # adds on top of the anchor with no re-baseline (D-72/D-78). ProcessSet sums the sinks, so
        # o2 splits by kᵢ/Σk — this is why SO₂ protection is EMERGENT (SO₂'s D-72 draw leaves less).
        r_o2 = params["k_anthocyanin_fade"] * f_t * o2 * anthocyanin  # g O2/L/h via the fade route
        d[schema.slice("o2")] = -r_o2
        # Anthocyanin faded per g O2 at a MASS-based yield (g anthocyanin / g O2), not a molar
        # stoichiometry — anthocyanin is a lumped pool with no clean molar mass (the y_ellag_per_o2
        # idiom). It is a pure TRANSFER to the colourless faded slot: faded gains exactly what
        # anthocyanin loses, so the D-81 colour identity closes by construction. Both off every
        # ledger, so this moves nothing conserved.
        faded = params["y_anthocyanin_per_o2"] * r_o2  # g anthocyanin/L/h
        d[schema.slice("anthocyanin")] = -faded
        d[schema.slice("faded_anthocyanin")] = faded
        return d


class ThermalAnthocyaninFade(Process):
    """Aging thermal fade: heat degrades free anthocyanin → colourless, WITHOUT O₂ (decision D-83).

    The eleventh aging Process and the **second, O₂-independent** fate that fades free monomeric
    ``anthocyanin`` to colourless — the pathway :class:`AnthocyaninFading` (D-81) explicitly
    deferred. Beyond the O₂-driven oxidative bleaching, free anthocyanin also degrades by a
    **thermal/hydrolytic** route that needs **no oxygen at all**: the flavylium ring slowly opens
    and the pigment breaks down to colourless products purely as a function of temperature and time
    (Somers & Evans; Ribéreau-Gayon, *Handbook of Enology* — anthocyanin thermal degradation). This
    is why a **sealed, anaerobic** red still loses its bright monomeric colour on the shelf, and why
    **warm storage kills red colour** even in an inert, fully-reductive bottle.

    ``r = k_anthocyanin_thermal_fade · f(T) · [anthocyanin]`` — **first-order** in the anthocyanin
    driver alone (the :class:`EsterHydrolysis` first-order form, *not* the bilinear O₂-sink form of
    its D-81 sibling), ``f(T) = arrhenius_factor(T, E_a_anthocyanin_thermal_fade, T_ref)`` the
    sourced warmer-fades-faster factor. It transfers anthocyanin into the **same** colourless
    ``faded_anthocyanin`` slot the D-81 oxidative fade fills — one colourless sink, now with two
    contributing routes (oxidative + thermal)::

        d(anthocyanin)/dt       = −r
        d(faded_anthocyanin)/dt = +r

    a pure g/L **transfer** between two off-ledger slots (the faded pool gains exactly what
    anthocyanin loses), so the D-81 colour identity ``anthocyanin + polymeric_pigment +
    faded_anthocyanin ≡ anthocyanin₀`` still closes by construction. **No yield** (contrast D-81):
    the rate is already in ``g anthocyanin/L/h`` — there is no O₂ pool to convert *through*, so the
    anthocyanin→faded transfer is directly ``−r``/``+r`` (a cleaner sink than the bilinear O₂
    route, which needed ``y_anthocyanin_per_o2`` only because its ``r`` was in O₂ units). Two
    params, ``k`` and ``E_a``, and no more.

    **O₂-INDEPENDENT, so SO₂ does NOT protect (the D-83 correctness crux, the mirror of D-81).**
    This is the deliberate opposite of :class:`AnthocyaninFading`. That route draws the **shared**
    ``o2``
    budget, so SO₂ protects it *emergently* (SO₂ scavenges O₂ via :class:`SulfiteOxidation`, D-72,
    leaving less to fade). Thermal fade touches **no** ``o2`` at all — it is not an oxidation — so
    **SO₂ gives no protection** against it: a heavily-sulfited red *still* fades thermally, and only
    **cold storage** slows it (the ``E_a > 0`` temperature lever). That is the physically-honest
    split, and exactly why this is a *separate* Process rather than an SO₂-insensitive term bolted
    onto the fade: attributing thermal loss to the O₂/SO₂ pathway would mis-explain the mechanism.
    ``ProcessSet`` never routes ``o2`` here, so a reductive (no ``add_oxygen``) red — flat under
    D-81 alone — now genuinely declines, and the D-81 "anaerobic sealed red holds its colour"
    scope-note is **retired** (that was the acknowledged v1 gap this beat closes). Colour loss only,
    to the **colourless** ``faded_anthocyanin`` sink — *not* browning (oxidative browning is
    :class:`PhenolicBrowning`/``A420``, D-73; this adds no ``A420`` and no second browning pathway).

    **Off every ledger, no conservation term (the :class:`AnthocyaninFading` precedent).** Both
    ``anthocyanin`` and ``faded_anthocyanin`` are grape-derived colour pools off
    ``total_carbon``/``total_mass``/``total_nitrogen`` (``iso_alpha`` precedent), so thermal fading
    moves **nothing conserved** — it touches only those two slots and asserts nothing. The transfer
    is exact regardless of any parameter value.

    **Wine-only + isolable + substrate-gated (prime directive #3).** ``anthocyanin`` /
    ``faded_anthocyanin`` are wine-only (appended to ``wine_schema``), so this is wired into the
    *wine* medium only; the ``"anthocyanin" not in schema`` guard makes it a hard no-op besides.
    Wired **disabled at the compile seam** (aging is post-ferment); ``begin_aging`` enables it with
    the other aging Processes. With no anthocyanin the ``anthocyanin ≤ 0`` guard returns
    byte-for-byte zero, so a white aging is exactly the case without this Process. Tier
    **speculative** (the aging axis is the Tier-3 frontier; the *form* — O₂-independent,
    anthocyanin-driven, warmer-faster, SO₂-**un**protected — is sourced, the rate/E_a magnitudes
    order-of-magnitude estimates). **Scope (v1):** the aggregate thermal/hydrolytic fade to
    colourless. See ``polymerization.yaml`` for the full scope + provenance.
    """

    name = "thermal_anthocyanin_fade"
    tier = Tier.SPECULATIVE
    #: TRANSFERS free anthocyanin into the colourless ``faded_anthocyanin`` slot — both off every
    #: ledger, so nothing conserved moves; it touches those two and nothing else. Notably it does
    #: NOT touch ``o2`` (the D-83 crux — this is the O₂-INDEPENDENT fade, so SO₂ cannot protect it).
    #: (``anthocyanin`` is now drawn by FOUR Processes — the two condensation routes, the D-81
    #: oxidative fade, and this thermal fade — the shared-pool precedent.)
    touches = ("anthocyanin", "faded_anthocyanin")
    #: ``k_anthocyanin_thermal_fade``/``E_a_anthocyanin_thermal_fade`` are this Process's own
    #: (polymerization.yaml, D-83); ``T_ref`` is shared with every Arrhenius rate. No yield
    #: (contrast D-81): the rate is already g anthocyanin/L/h, a direct transfer. Their tiers cap
    #: the ``anthocyanin``/``faded_anthocyanin`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_anthocyanin_thermal_fade",
        "E_a_anthocyanin_thermal_fade",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (anthocyanin/faded_anthocyanin are appended to wine_schema): a hard no-op
        # on any schema without them, belt-and-suspenders to the wine-only wiring.
        if "anthocyanin" not in schema:
            return d
        anthocyanin = float(y[schema.slice("anthocyanin")][0])
        # No anthocyanin ⇒ no fading: a white aging is byte-for-byte the case without this Process.
        # Gate on the anthocyanin STATE before reading any fade param (the OakExtraction/Strecker
        # substrate-gate-before-params discipline — an enabled-but-undosed Process mustn't KeyError
        # if polymerization.yaml is absent). ``<= 0`` also absorbs solver undershoot. NOTE: unlike
        # the D-81 oxidative fade there is NO o2 gate — this route needs no oxygen, so it fires even
        # in a fully reductive (no add_oxygen) bottle (the whole point of the thermal pathway).
        if anthocyanin <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_anthocyanin_thermal_fade"], params["T_ref"])
        # First-order thermal fade (the EsterHydrolysis first-order form, NOT the D-81 bilinear
        # o2*anthocyanin form): r is the anthocyanin degradation rate (g anthocyanin/L/h). No o2
        # term — this is the O₂-INDEPENDENT pathway, so SO₂ (an antioxidant) gives no protection;
        # only cold storage (E_a > 0) slows it. It is a pure TRANSFER to the colourless faded slot:
        # faded gains exactly what anthocyanin loses, so the D-81 colour identity closes by
        # construction. Both off every ledger, so this moves nothing conserved.
        r = params["k_anthocyanin_thermal_fade"] * f_t * anthocyanin  # g anthocyanin/L/h
        d[schema.slice("anthocyanin")] = -r
        d[schema.slice("faded_anthocyanin")] = r
        return d


class TanninSelfPolymerization(Process):
    """Aging condensation: grape tannin self-polymerizes → soft polymer (decision D-84).

    The twelfth aging Process, the **fourth non-oxidative** one, and the first of the
    **tannin–tannin axis** that the D-79/D-80 condensation beats explicitly deferred (their
    "one-directional-per-pool" honesty note). Beyond condensing with anthocyanin (D-79/D-80),
    condensed grape **tannin** also reacts **with itself**: flavan-3-ol units link into larger
    polymers that taste **softer** — a direct (non-oxidative) self-condensation, the well-known
    "tannins polish/soften with age even with nothing else to react with" behaviour
    (Ribéreau-Gayon, *Handbook of Enology*; the
    proanthocyanidin-polymerization literature). This is why a **white** wine's tannin, or a red
    whose anthocyanin has been exhausted, **still softens** on aging — softening that the
    anthocyanin-dependent condensation routes alone could not produce.

    ``r = k_tannin_self_polymerization · f(T) · [tannin]²`` — **bimolecular** in the single tannin
    pool (a true *self*-reaction, so second-order in ``[tannin]`` — distinct from the *bilinear*
    two-pool D-79 form), ``f(T) = arrhenius_factor(T, E_a_tannin_self_polymerization, T_ref)`` the
    sourced warmer-polymerizes-faster factor. It is a **pure off-ledger tannin sink**::

        d(tannin)/dt = −r

    The polymerized tannin goes to **no destination slot** — deliberately, and consistently with the
    D-79/D-80 condensation routes, which already consume ``tannin`` as a pure sink (the shared
    ``polymeric_pigment`` slot they fill is deposited in **anthocyanin**-equivalents, *never* tannin
    mass; no ledger reads tannin mass). Adding a ``polymerized_tannin`` slot here but not for the
    condensation-consumed tannin would be asymmetric bookkeeping for a pool nothing conserved reads,
    so ``r`` folds the lumped self-condensation stoichiometry into one first-class sink rate and the
    tannin simply declines (the polymer is soft, so it drops out of astringency — see below). **No
    yield** (a self-reaction: there is no *second* pool to consume at a ratio).

    **The astringency payoff, and the honesty note it retires.** ``tannin`` is read as astringency
    (:func:`~fermentation.analysis.astringency_series`, mg/L free harsh tannin); the soft polymer is
    excluded, so drawing ``tannin`` down **softens** the wine — exactly as the D-79 route does, but
    **without needing anthocyanin**. Through D-80 the astringency readout carried a standing caveat
    ("grape tannin self-polymerization … a further-deferred beat, so anthocyanin is the limiting
    reagent and A–T condensation softens only modestly"); this Process builds that beat, so a
    no-anthocyanin (white / tannin-only) wine now genuinely softens and that note is **retired**.

    **Off every ledger, no conservation term (the :class:`TanninAnthocyaninCondensation`
    precedent).** ``tannin`` is grape-derived — off ``total_carbon``/``total_mass``/
    ``total_nitrogen`` (the ``iso_alpha``/``ellagitannin`` precedent) — so self-polymerization moves
    **nothing conserved**: it touches only that one slot and asserts nothing. Unlike the D-80
    bridged route it consumes **no**
    on-ledger acetaldehyde (this is the *direct*, acetaldehyde-free self-condensation — the
    acetaldehyde-bridged tannin–ethyl–tannin variant is the separate D-85 beat), so there is no
    split-ledger carbon capture: a wholly off-ledger sink.

    **Oak-independent AND O₂-independent.** Like :class:`TanninAnthocyaninCondensation`, this is a
    **grape**-tannin reaction: it draws **no** ``o2`` (not an oxidation) and reads **no** oak pool —
    ``tannin`` is grape **condensed** tannin, distinct from oak **hydrolysable** ``ellagitannin``
    (D-78), so a steel-tank red with no oak and no oxygen still self-polymerizes and softens.

    **Wine-only + isolable + substrate-gated (prime directive #3).** ``tannin`` is wine-only
    (appended to ``wine_schema``), so this is wired into the *wine* medium only; the ``"tannin" not
    in schema`` guard makes it a hard no-op besides. Wired **disabled at the compile seam** (aging
    is post-ferment); ``begin_aging`` enables it with the other aging Processes. With no tannin the
    ``tannin ≤ 0`` guard returns byte-for-byte zero, so a no-tannin run is exactly the case without
    this Process. Tier **speculative** (the aging axis is the Tier-3 frontier; the *form* —
    bimolecular self-condensation, warmer-faster, softening — is sourced, the rate/E_a magnitudes
    order-of-magnitude estimates). **Scope (v1):** the direct (acetaldehyde-free) self-condensation
    only; the acetaldehyde-bridged tannin–ethyl–tannin route is the separate D-85 beat. See
    ``polymerization.yaml`` for the full scope + provenance.
    """

    name = "tannin_self_polymerization"
    tier = Tier.SPECULATIVE
    #: Consumes the single grape ``tannin`` pool as a pure off-ledger sink (the soft polymer goes
    #: to no slot — the D-79/D-80 tannin-is-a-pure-sink precedent, since no ledger reads tannin
    #: mass), so nothing conserved moves; it touches that one slot and nothing else. No ``o2`` (not
    #: an oxidation), no ``acetaldehyde`` (the DIRECT route — bridged tannin–ethyl–tannin is D-85).
    #: (``tannin`` is now drawn by FOUR Processes — the two condensation routes, this, and the
    #: D-85 tannin–ethyl–tannin bridge.)
    touches = ("tannin",)
    #: ``k_tannin_self_polymerization``/``E_a_tannin_self_polymerization`` are this Process's own
    #: (polymerization.yaml, D-84); ``T_ref`` is shared with every Arrhenius rate. No yield (a
    #: self-reaction — one pool, no second reactant to consume at a ratio). Their tiers cap the
    #: ``tannin`` output tier via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_tannin_self_polymerization",
        "E_a_tannin_self_polymerization",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slot (tannin is appended to wine_schema): a hard no-op on any schema without it,
        # belt-and-suspenders to the wine-only wiring.
        if "tannin" not in schema:
            return d
        tannin = float(y[schema.slice("tannin")][0])
        # No tannin ⇒ no self-polymerization: a no-tannin run is byte-for-byte the case without this
        # Process. Gate on the tannin STATE before reading any param (the OakExtraction/Strecker
        # substrate-gate-before-params discipline — an enabled-but-undosed Process mustn't KeyError
        # if polymerization.yaml is absent). ``<= 0`` also absorbs solver undershoot.
        if tannin <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_tannin_self_polymerization"], params["T_ref"])
        # Bimolecular self-condensation ([tannin]², NOT the D-79 bilinear two-pool form): a true
        # self-reaction, so second-order in the single tannin pool. r is the tannin consumption rate
        # (g tannin/L/h); the lumped self-condensation stoichiometry is folded into k. No o2 term
        # (not an oxidation) and no acetaldehyde (the DIRECT route — bridged tannin–ethyl–tannin is
        # D-85). A pure OFF-LEDGER sink: the soft polymer goes to no slot (the D-79/D-80 precedent —
        # no ledger reads tannin mass), so this moves nothing conserved. Astringency (which reads
        # the free tannin) softens as the pool declines — WITHOUT anthocyanin (the D-80 honesty
        # note retired).
        r = params["k_tannin_self_polymerization"] * f_t * tannin * tannin  # g tannin/L/h
        d[schema.slice("tannin")] = -r
        return d


class TanninEthylTanninCondensation(Process):
    """Aging condensation: acetaldehyde bridges tannin + tannin → soft polymer (decision D-85).

    The thirteenth aging Process, the **fifth non-oxidative** one, and the second of the
    **tannin–tannin axis** — the acetaldehyde-bridged **tannin–ethyl–tannin** route that both the
    D-80 bridged beat and the D-84 direct self-polymerization explicitly deferred. It is the
    acetaldehyde-accelerated sibling of :class:`TanninSelfPolymerization` (D-84), exactly as
    :class:`AcetaldehydeBridgedCondensation` (D-80) is of :class:`TanninAnthocyaninCondensation`
    (D-79): dissolved-O₂ acetaldehyde (:class:`OxidativeAcetaldehyde`, D-71) forms an **ethylidene
    bridge** ``—CH(CH₃)—`` linking **two grape tannin** flavanols (tannin–ethyl–tannin), an
    acetaldehyde-accelerated self-condensation that **softens astringency**. So micro-oxygenation
    (``add_oxygen``) now softens even an **anthocyanin-free** tannin pool — a white / tannin-only
    wine's tannin polymerizes *faster* under O₂ — while simultaneously drawing acetaldehyde down
    (the real micro-ox benefit, mouthfeel *and* reduced harshness), with no colour involved.

    ``r = k_tannin_ethyl_tannin · f(T) · [acetaldehyde_free] · [tannin]²`` — the D-84
    **bimolecular** self-condensation form (second-order in ``[tannin]``, a true tannin–tannin
    reaction) **plus** the D-80 free-acetaldehyde factor,
    ``f(T) = arrhenius_factor(T, E_a_tannin_ethyl_tannin, T_ref)`` its **own** E_a (a distinct
    reaction — prime directive #2, the :class:`AcetaldehydeBridgedCondensation` vs
    :class:`OakExtraction` precedent). Anchored on **tannin** consumption (there is no anthocyanin
    here — both bridge ends are flavanols), so ``r`` is the tannin drawdown rate::

        d(tannin)/dt       = −r                                       # off-ledger (grape-derived)
        d(acetaldehyde)/dt = −y_acetaldehyde_per_tannin · r           # ON the carbon ledger
        d(ethyl_bridge)/dt = +(acetaldehyde carbon consumed) / c(ethylidene)   # ON-ledger, exact

    **The split ledger (reused verbatim from D-80).** The grape ``tannin`` bulk is **off** every
    ledger (grape-derived, untracked — the ``iso_alpha``/``ellagitannin`` precedent), so consuming
    it moves nothing conserved, exactly as the D-84 direct route. But acetaldehyde's carbon is
    **on** the carbon ledger (borrowed from ethanol ``E`` at D-71), so consuming it into an
    *off*-ledger polymer would make carbon vanish and fail
    :func:`~fermentation.validation.conservation.assert_conserved`. The fix is the **same on-ledger
    ``ethyl_bridge`` slot** D-80 introduced (weighted at ``carbon_mass_fraction("ethylidene")`` in
    :func:`~fermentation.validation.conservation.total_carbon`), filled by the **same**
    carbon-exact split — release at ``c(acetaldehyde)``, re-deposit at ``c(ethylidene)``
    (acetaldehyde loses only its carbonyl O as water). So ``total_carbon`` closes to **machine
    precision** non-trivially (acetaldehyde↓ exactly equals ``ethyl_bridge``↑ in carbon), and both
    bridged routes (D-80 anthocyanin, D-85 tannin) deposit into **one shared** ``ethyl_bridge``
    pool — its meaning is the ethylidene bridge carbon, whether the bridge terminates in pigment
    (D-80) or a tannin–tannin polymer (here).

    **A DIFFERENT acetaldehyde yield from D-80 (prime directive #2).** One acetaldehyde bridges
    **two flavanols** here — a different lumped stoichiometry from D-80's flavanol↔anthocyanin
    bridge — so this reads its **own** ``y_acetaldehyde_per_tannin`` (g acetaldehyde per g tannin
    consumed), *not* D-80's ``y_acetaldehyde_per_anthocyanin``. Reusing D-80's yield would fake a
    shared stoichiometry the two reactions do not have. The yield only sets the **magnitude** of
    the acetaldehyde drawdown; carbon balances for any value (the D-80 property).

    **Deposits NO pigment — the colour difference from D-80.** Both bridge ends are colourless
    flavanols (tannin–ethyl–tannin), so — unlike D-80 (tannin–ethyl–anthocyanin, which deposits
    ``polymeric_pigment`` in anthocyanin-equivalents) — this adds **nothing** to colour: it touches
    no ``anthocyanin`` and no ``polymeric_pigment``. It is a pure astringency softener (draws
    ``tannin``, like D-84), an O₂-driven one. The soft tannin–ethyl–tannin polymer goes to **no**
    destination slot (the D-84/D-79/D-80 tannin-is-a-pure-sink precedent — no ledger reads tannin
    mass), so besides the on-ledger acetaldehyde carbon it captures, it books nothing.

    **Reads FREE acetaldehyde, not total (the D-47/D-80 precedent).** SO₂-bound acetaldehyde is the
    bisulfite adduct — its carbonyl is blocked, so it **cannot** form the ethylidene bridge — so
    when SO₂ is dosed (``so2_total > 0``) the rate reads
    :func:`~fermentation.core.acidbase.free_acetaldehyde`; the guard is exact (an unsulfited run
    pays no per-RHS pH ``brentq`` and is byte-for-byte the total-acetaldehyde case). **Emergent:**
    SO₂ delays acetaldehyde-mediated tannin softening exactly as it delays the D-80 colour
    stabilization — the same shared binding equilibrium, nothing scripted.

    **Triply substrate-gated + wine-only + isolable (prime directive #3).** The rate is
    ``[acetaldehyde] · [tannin]²``, so a no-tannin / no-acetaldehyde wine is byte-for-byte the case
    without this Process (and all of beer — the slots are wine-only). Gated on the tannin STATE
    before reading any polymerization param (the substrate-gate-before-params discipline). Wired
    **disabled at the compile seam** (aging is post-ferment); ``begin_aging`` enables it with the
    other aging Processes. Tier **speculative** (the aging axis is the Tier-3 frontier; the *form* —
    acetaldehyde-bridged tannin–tannin, warmer-faster, SO₂-blocked, softening, colourless — is
    sourced, the magnitudes order-of-magnitude estimates). See ``polymerization.yaml`` for the full
    scope + provenance.
    """

    name = "tannin_ethyl_tannin_condensation"
    tier = Tier.SPECULATIVE
    #: Consumes the grape ``tannin`` pool (off every ledger — the D-84 pure-sink precedent, the
    #: soft tannin–ethyl–tannin polymer goes to no slot) plus the on-ledger ``acetaldehyde``, whose
    #: carbon it captures in the on-ledger ``ethyl_bridge`` slot (the D-80 split-ledger transfer,
    #: one shared bridge pool). Touches those three and nothing else — the ``acetaldehyde →
    #: ethyl_bridge`` carbon closes exactly; the tannin (off-ledger) moves nothing conserved.
    #: Notably it touches NO ``anthocyanin`` and NO ``polymeric_pigment`` (a colourless
    #: tannin–tannin polymer — the colour difference from D-80, which bridges to anthocyanin).
    touches = ("acetaldehyde", "ethyl_bridge", "tannin")
    #: ``k_tannin_ethyl_tannin``/``E_a_tannin_ethyl_tannin``/``y_acetaldehyde_per_tannin`` are this
    #: Process's own (polymerization.yaml, D-85 — its own acetaldehyde yield, distinct from D-80's
    #: ``y_acetaldehyde_per_anthocyanin`` since one acetaldehyde bridges TWO flavanols here);
    #: ``T_ref`` is shared with every Arrhenius rate. The SO₂/pH params read inside
    #: :func:`free_acetaldehyde`/:func:`ph_of_state` are omitted (all plausible, the Process is
    #: already speculative — the D-47/D-80 precedent). Tiers cap the outputs via D-1.
    reads: tuple[str, ...] = (
        "k_tannin_ethyl_tannin",
        "E_a_tannin_ethyl_tannin",
        "y_acetaldehyde_per_tannin",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (tannin/ethyl_bridge are appended to wine_schema): a hard no-op on any
        # schema without them, belt-and-suspenders to the wine-only wiring.
        if "tannin" not in schema:
            return d
        tannin = float(y[schema.slice("tannin")][0])
        acetaldehyde = float(y[schema.slice("acetaldehyde")][0])
        # Doubly substrate-gated: no tannin OR no acetaldehyde ⇒ no bridging, so a no-tannin /
        # no-acetaldehyde wine is byte-for-byte the case without this Process. Gate on the STATE
        # before reading any polymerization param (the substrate-gate-before-params discipline — an
        # enabled-but-undosed Process mustn't KeyError if polymerization.yaml is absent). ``<= 0``
        # also absorbs solver undershoot.
        if tannin <= 0.0 or acetaldehyde <= 0.0:
            return d
        # SO₂-bound acetaldehyde is the bisulfite adduct — its carbonyl is blocked, so it CANNOT
        # form the ethylidene bridge; read only the FREE share under SO₂ (the D-47/D-80 precedent).
        # The ``so2_total > 0`` guard is EXACT — an unsulfited run pays no per-RHS pH ``brentq`` and
        # its rate is byte-for-byte the total-acetaldehyde case. Emergent: SO₂ delays acetaldehyde-
        # mediated tannin softening (bound acetaldehyde is unavailable to bridge).
        bridging_acetaldehyde = acetaldehyde
        if SO2_STATE_KEY in schema and float(y[schema.slice(SO2_STATE_KEY)][0]) > 0.0:
            bridging_acetaldehyde = free_acetaldehyde(
                y, schema, params, ph_of_state(y, schema, params)
            )
            if bridging_acetaldehyde <= 0.0:  # all acetaldehyde bound ⇒ none available to bridge
                return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_tannin_ethyl_tannin"], params["T_ref"])
        # Bimolecular self-condensation ([tannin]², the D-84 form) ACCELERATED by free acetaldehyde
        # (the D-80 factor): r is the tannin drawdown rate (g tannin/L/h), anchored on tannin (no
        # anthocyanin here — both ends are flavanols). Consumes NO colour: touches neither
        # anthocyanin nor polymeric_pigment (a colourless tannin–tannin polymer — the D-80 colour
        # difference). The soft polymer goes to no slot (the D-84 pure-sink precedent).
        r = params["k_tannin_ethyl_tannin"] * f_t * bridging_acetaldehyde * tannin * tannin
        d[schema.slice("tannin")] = -r
        # The ON-ledger half — the D-80 split-ledger carbon capture, reused verbatim. Acetaldehyde
        # consumed by the bridge (g/L/h) at this route's OWN yield (one acetaldehyde per two
        # flavanols — distinct from D-80's per-anthocyanin yield); its carbon must NOT vanish into
        # the off-ledger tannin polymer, so re-deposit it into the on-ledger ethyl_bridge pool via
        # the carbon-exact split (release at acetaldehyde's C fraction, re-deposit at ethylidene's —
        # acetaldehyde loses only its carbonyl O as water). total_carbon closes to machine prec.:
        # acetaldehyde C leaving == bridge C arriving. The yield only sets the MAGNITUDE; carbon
        # balances for any value.
        acet_consumed = params["y_acetaldehyde_per_tannin"] * r  # g acetaldehyde/L/h
        d[schema.slice("acetaldehyde")] = -acet_consumed
        carbon_released = acet_consumed * carbon_mass_fraction(
            _BRIDGE_ACETALDEHYDE_SPECIES
        )  # g C/L/h
        d[schema.slice("ethyl_bridge")] = carbon_released / carbon_mass_fraction(
            _ETHYL_BRIDGE_SPECIES
        )
        return d


class SMMHydrolysis(Process):
    """DMS accumulation by hydrolysis of the grape-borne precursor pool (decision D-102).

    ``d(dms_potential)/dt = -k_smm_hydrolysis · f(T) · dms_potential`` and ``d(dms)/dt = +`` the
    same rate — a first-order Arrhenius decay of the precursor **straight into the product**, with
    ``f(T) = arrhenius_factor(T, E_a_smm_hydrolysis, T_ref)``. Fills the ``dms`` pool with the
    aged-wine "truffle / black olive / cooked corn" odorant, whose OAV the D-67 lens already reads.

    **A DISTINCT route — this is why the beat exists, and why it needed no ratio-split.** Every
    other sulfur pool in the model is autolysis-gated
    (:class:`~fermentation.core.kinetics.hydrogen_sulfide.AutolyticHydrogenSulfide` D-44,
    :class:`~fermentation.core.kinetics.mercaptans.AutolyticMercaptan` D-45). DMS is not: it
    accumulates by **spontaneous hydrolysis of S-methylmethionine (SMM) during bottle aging**, on
    the lees or off them, and so it carries its **own anchor** rather than dividing a shared
    autolytic yield into invented sub-yields. That is exactly the D-96 linchpin the ``mercaptans``
    speciation could not satisfy (D-101: you cannot split a pool that only ever had one molecule
    in it), and it is satisfiable here.

    **In DMS-EQUIVALENTS, so the conversion is 1:1 by construction.** ``dms_potential`` is g of
    *releasable DMS*, not g of SMM — the unit the wine literature reports DMSp in, because DMSp is
    an assay (cleave every precursor with hot alkali, measure the DMS that comes off). So there is
    **no stoichiometric yield parameter to guess** and SMM's molar-mass / iodide-salt-form
    ambiguity never enters the model. ``dms_potential + dms`` is invariant to machine precision —
    a real conservation law on two off-ledger pools, and a test.

    **No pH term — a SOURCED NULL, not an omission.** SMM hydrolysis in wine is a nucleophilic
    substitution of the sulfur by water, and below pH 5 that mechanism's rate is pH-insensitive:
    De La Burgade et al. 2025 measured SMM degradation at pH 2.8 vs 3.8 and found no significant
    difference. Wine's whole range sits inside the insensitive regime, so a pH term would be a
    knob with a measured value of zero. (This is also, exactly, why Scheuren's brewing activation
    energy is not used — see ``dms.yaml``'s header: his wort at pH 5.2 is in the *other*
    mechanism's regime, and E_a is a property of a mechanism. Transferred here it would predict a
    ~53-year cellar half-life, i.e. that aged wine never develops DMS at all.)

    **Not flux-linked, and cumulative — the D-42 contrast is the point.** Like the D-44/D-45
    reductive faults it fires post-dryness, so there is no CO₂ stream to strip it: unlike ``h2s``
    (whose D-42 stripping sink makes it a *residual*), ``dms`` only rises. ``d(dms)/dt >= 0``
    always (monotonic, no clamp) — a decay of a non-negative pool by a positive rate, so the
    precursor approaches 0 from above and needs no floor and no C¹ shadow gate.

    **Conservation — moves NOTHING conserved, by construction (the D-74 A420 argument).** Both
    slots are off every ledger, so the precursor's carbon comes from untracked SMM and lands in
    untracked DMS without ever touching ``total_carbon``/``total_nitrogen``. This is *cleaner*
    than :class:`~fermentation.core.kinetics.mercaptans.AutolyticMercaptan`, which must debit real
    methionine for its thiol's carbon. At µg/L the amounts are ~1e-6 of the carbon ledger anyway,
    so the off-ledger call costs nothing measurable and buys exactness.

    **SCOPE — three overstatements, ALL pushing the same way, and that way is conservative.**
    Stated here because together they bound every number this Process emits:

    * **DMSp ⊃ SMM.** The pool is a DMSp *assay* value (all precursors, incl. DMSPA and unknowns);
      the rate is SMM's. SMM is only 21–74% of DMS formed, so the substrate is over-stated by up
      to ~4×.
    * **Variety.** ``dms_potential_initial`` defaults generic wine to the only at-bottling dataset
      — six **Syrah** wines, a high-DMSp variety. A low-DMSp wine is over-predicted, measurably:
      see D-102 on the Amarone miss, which is **recorded, not tuned**.
    * **No closure permeation.** Real DMS escapes through cork (De La Burgade's headline result:
      up to 12% of initial). Modelling the bottle as sealed keeps DMS in. That is a closure-
      specific, orthogonal axis and a beat of its own — not a re-anchoring of ``k``.

    All three make the model cry fault **earlier** than reality. For an off-aroma that is the safe
    direction; a model that under-predicted a fault would be the dangerous one.

    **A FOURTH SCOPE, and it cuts the OTHER WAY — so "all conservative" is NOT the whole story.**
    This Process reports **DMS formed by aging** and seeds ``dms = 0``: it models no *at-bottling*
    DMS. Real wines have it, and substantially — the six Syrah wines carried **29.9–314.9 µg/L
    before any aging** (fermentative and yeast-mediated routes, D-102's scope). So the model
    **under-counts** an aged wine's *total* DMS, partly offsetting the three overstatements above.
    **The consequence that matters is a comparison rule, not a magnitude:** this pool is *not*
    commensurate with a measured aged-wine DMS total, and must not be validated against one (D-102
    made exactly that error and corrected it). Modelling the pre-aging DMS needs the fermentative
    and DMSO-reduction routes — a beat of its own.

    **What it cannot express, and it is DMS's most-studied property.** DMS *enhances* fruity-ester
    perception (it lowers the fruity pool's effective threshold) rather than adding to it. The OAV
    lens is structurally unable to carry a perceptual interaction (masking, blocked on ``cosα``
    since D-95/D-98), so the model tracks DMS's concentration honestly and stays **silent** on its
    best-documented effect. See ``threshold_dms_wine``'s notes.

    **Isolability — wine-only and aging-gated (prime directive #3).** Touches only
    ``{dms_potential, dms}`` — both new, both wine-only, both off every ledger — so no existing
    pool, ledger, tier or test can move whether it fires or not. Wired into wine's aging tuple and
    disabled at compile with the rest of the aging axis; an un-aged wine is byte-for-byte the
    pre-D-102 model. Tier **speculative** (the anchors are accelerated-aging extrapolations down
    to a cellar; see ``dms.yaml``).
    """

    name = "smm_hydrolysis"
    tier = Tier.SPECULATIVE
    #: Decays ``dms_potential`` into ``dms`` 1:1 in DMS-equivalents — an inter-pool transfer
    #: between two OFF-ledger slots, so it touches these two and nothing else. No carbon pool, no
    #: nitrogen pool, no ``o2`` (this is a hydrolysis, not an oxidation — a sealed, un-oxygenated,
    #: sulfited bottle makes DMS exactly as fast; the D-83 O₂-independent relationship).
    touches = ("dms_potential", "dms")
    #: ``k_smm_hydrolysis``/``E_a_smm_hydrolysis`` are this Process's own (dms.yaml, D-102);
    #: ``T_ref`` is shared with every other Arrhenius rate. Their tiers cap the ``dms`` output tier
    #: via parameter-tier propagation (D-1). Note ``dms_potential_initial`` is NOT read here — it
    #: seeds the slot at the compile seam (a must property); it does not drive the rate.
    reads: tuple[str, ...] = (
        "k_smm_hydrolysis",
        "E_a_smm_hydrolysis",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        potential = float(y[schema.slice("dms_potential")][0])
        if potential <= 0.0:
            # No precursor ⇒ no DMS (the un-seeded no-op). Also absorbs a solver undershoot: a
            # first-order decay cannot drive the pool negative, so this is a guard, not a clamp
            # on real physics.
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_smm_hydrolysis"], params["T_ref"])
        rate = params["k_smm_hydrolysis"] * f_t * potential  # g eq. DMS/L/h

        # 1:1 in DMS-equivalents (see the docstring): the precursor pool IS denominated in the DMS
        # it can release, so what leaves one slot enters the other with no yield and no molar-mass
        # conversion. dms_potential + dms is invariant to machine precision.
        d[schema.slice("dms_potential")] = -rate
        d[schema.slice("dms")] = rate
        return d
