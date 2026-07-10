"""Aging chemistry — the slow, post-fermentation "years" axis (§4.1, decision D-69).

The first Tier-3 aging Process, opened in D-68 and built here. Unlike the Milestone-2
byproduct Processes (which *produce* aroma pools during active fermentation), an aging
Process acts on a *finished* wine/beer over months-to-years: fermentation is done, the
sugar is gone, the yeast racked or crashed, and the chemistry that remains is spontaneous
(hydrolysis, oxidation, condensation), not metabolic. This module holds the aging
Processes; :class:`EsterHydrolysis` is the first, :class:`OxidativeAcetaldehyde` the second,
:class:`SulfiteOxidation` the third, and :class:`PhenolicBrowning` the fourth.

**The oxidative sub-axis (D-71).** :class:`OxidativeAcetaldehyde` opens the *oxidative* half of
the aging axis on a **dissolved-O₂ pool** (``o2``, a new carbon-free state slot, off every
conservation ledger like ``h2s``/``iso_alpha``). O₂ — not ethanol — is the rate-limiting
reactant (ethanol sits at ~100 g/L, effectively constant across aging), so the pool is *the*
substrate that bounds oxidation: acetaldehyde saturates as the O₂ charge is spent, the
bottle-aging reality a first-order-in-ethanol rate could never reproduce (it would grow
unbounded). O₂ enters via a dedicated ``add_oxygen`` dosing verb (one dose = a bottle's
ingress; repeated = micro-oxygenation / barrel), and a ``begin_aging`` run with **no** O₂ dosed
is purely *reductive* aging (screwcap/inert) — byte-for-byte the :class:`EsterHydrolysis`-only
aging, since the Process contributes exactly zero at ``o2 = 0``. Oxidative aging is
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

    d(esters)/dt = -k_ester_hydrolysis · f(T) · max(0, esters - esters_eq)

is **net decay toward a lower equilibrium floor** ``esters_eq``, *not* decay-to-zero (D-68):
below ``esters_eq`` the rate is zero. The bidirectional reality — ethyl esters of fatty
acids sit *below* equilibrium young and slowly *form* on aging — is the deferred half;
framing the acetate-dominated lump as "net decay toward a lower floor" is the same
fixed-composition honesty the D-67 sensory lump carries. ``f(T) = arrhenius_factor(T,
E_a_ester_hydrolysis, T_ref)`` gives the sourced **warmer-ages-faster** direction (cold
cellars preserve fruity esters). No fermentative-flux gate — aging runs when the flux is
zero — so unlike the M2 producers this Process is driven by temperature and the pool alone.

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

from fermentation.core.acidbase import SO2_STATE_KEY, bisulfite_so2_at_ph, ph_of_state
from fermentation.core.chemistry import (
    M_ACETALDEHYDE,
    M_CO2,
    M_ETHANOL,
    M_METHIONAL,
    M_O2,
    M_PHENYLACETALDEHYDE,
    M_SO2,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics.amino_acids import AMINO_ACID_SPECIES
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Representative species that carbon-account each pool the hydrolysis touches, from the one
#: chemistry source of truth (mirrors the ``_ESTER_SPECIES``/``_FUSEL_SPECIES`` discipline in
#: :mod:`~fermentation.core.kinetics.byproducts`). Esters book as ethyl acetate (D-19, the
#: immovable pool weighting), the fusel product as isoamyl alcohol, the acid product (``Byp``)
#: as succinic acid (D-16). Using these fractions both to release and to re-deposit the carbon
#: is what makes the transfer close in ``total_carbon`` exactly.
_ESTER_SPECIES = "ethyl_acetate"
_FUSEL_SPECIES = "isoamyl_alcohol"
_BYP_SPECIES = "succinic_acid"

#: The 5:2 carbon split of the released ester carbon between ``fusels`` and ``Byp`` (D-69),
#: set by the **isoamyl-acetate** stand-in reaction (see the class docstring for why isoamyl
#: acetate, not the ethyl-acetate mass species, sets the split): isoamyl alcohol carries 5
#: carbons, acetic acid 2. Stoichiometry of the named stand-in reaction — a code-with-citation
#: constant like the chemistry carbon counts, not an empirical/uncertain YAML parameter.
_ISOAMYL_ALCOHOL_CARBONS = 5  # the alcohol product → fusels
_ACETIC_ACID_CARBONS = 2  # the acid product → Byp
_FUSEL_CARBON_SHARE = _ISOAMYL_ALCOHOL_CARBONS / (_ISOAMYL_ALCOHOL_CARBONS + _ACETIC_ACID_CARBONS)
_BYP_CARBON_SHARE = _ACETIC_ACID_CARBONS / (_ISOAMYL_ALCOHOL_CARBONS + _ACETIC_ACID_CARBONS)

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

#: The Strecker decarboxylation releases exactly **1 mol CO₂ per mol aldehyde** — the amino acid's
#: carboxyl carbon. On the carbon ledger (unlike ``o2``), so it is a genuine product term the carbon
#: bookkeeping must route, not an off-ledger emission (D-75).
_CO2_PER_STRECKER_ALDEHYDE = 1.0


class EsterHydrolysis(Process):
    """Aging hydrolysis of fruity acetate esters toward equilibrium (decision D-69).

    ``d(esters)/dt = -k_ester_hydrolysis · f(T) · max(0, esters - esters_eq)`` — first-order
    net decay of the lumped ``esters`` pool toward the lower equilibrium floor ``esters_eq``
    (not to zero), with ``f(T) = arrhenius_factor(T, E_a_ester_hydrolysis, T_ref)`` the
    sourced warmer-ages-faster factor. The released ester carbon
    (``rate·c(ethyl_acetate)``, ledger-fixed by the D-19 esters weighting) is split **5:2**
    into ``fusels`` (isoamyl alcohol, the alcohol product) and ``Byp`` (succinic-stand-in
    acetic acid, the acid product) — so aging fades the ester OAV, raises the fusel OAV, and
    drifts VA/pH up. See the module docstring for the full carbon algebra, the 5:2 split
    rationale (the advisor-settled crux), and the §4.3 firewall / stand-in seams it inherits.

    Off during the ferment (no fermentative-flux gate; it is temperature- and pool-driven);
    enabled only in a post-fermentation aging segment (D-68/D-70). Tier **speculative**.
    """

    name = "ester_hydrolysis"
    tier = Tier.SPECULATIVE
    #: Decays its own ``esters`` pool and routes the released carbon to the alcohol product
    #: (``fusels``) and the acid product (``Byp``) — an on-ledger inter-pool transfer, so it
    #: touches those three and nothing else (no ``S``/``E``/``CO2``; aging draws no sugar).
    touches = ("esters", "fusels", "Byp")
    #: ``k_ester_hydrolysis``/``E_a_ester_hydrolysis``/``esters_eq`` are this Process's own
    #: (aging.yaml, D-69); ``T_ref`` is shared with every other Arrhenius rate. Their tiers cap
    #: the ``esters``/``fusels``/``Byp`` output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_ester_hydrolysis", "E_a_ester_hydrolysis", "esters_eq", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        esters = float(y[schema.slice("esters")][0])
        # Net decay toward the equilibrium floor: the excess above esters_eq, never below zero
        # (below the floor there is no net hydrolysis — the reverse formation is deferred, D-68).
        # max(0, ...) with esters_eq > 0 also absorbs a solver undershoot (esters < 0 ⇒ 0).
        excess = max(0.0, esters - params["esters_eq"])
        if excess <= 0.0:
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_ester_hydrolysis"], params["T_ref"])
        rate = params["k_ester_hydrolysis"] * f_t * excess  # g esters/L/h decayed

        # The released ester carbon is ledger-fixed by the pool's (immovable) ethyl-acetate
        # weighting; split it 5:2 and re-deposit through each product pool's own carbon
        # fraction, so total_carbon closes to machine precision for any split summing to 1.
        carbon_released = rate * carbon_mass_fraction(_ESTER_SPECIES)  # g C/L/h
        d[schema.slice("esters")] = -rate
        d[schema.slice("fusels")] = (
            _FUSEL_CARBON_SHARE * carbon_released / carbon_mass_fraction(_FUSEL_SPECIES)
        )
        d[schema.slice("Byp")] = (
            _BYP_CARBON_SHARE * carbon_released / carbon_mass_fraction(_BYP_SPECIES)
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
    aldehyde carbon is drawn from ``amino_acids`` (booked as arginine) and the amino-acid nitrogen
    is **deaminated** back to the ``N`` pool, exactly as
    :class:`~fermentation.core.kinetics.mercaptans.AutolyticMercaptan` does; the Strecker
    **decarboxylation** adds one product this idiom did not have — **1 mol CO₂ per mol aldehyde**
    (:data:`_CO2_PER_STRECKER_ALDEHYDE`, the acid's carboxyl carbon), on the carbon ledger. The
    arginine draw is *sized to the product carbon* (methional + phenylacetaldehyde + CO₂), so
    ``total_carbon`` closes to machine precision (the :class:`EsterHydrolysis` multi-product split
    idiom); all the arginine nitrogen lands in ``N`` and the products are nitrogen-free, so
    ``total_nitrogen`` closes. The arginine-for-``amino_acids`` stand-in is **exact on the ledger,
    approximate on provenance** (the drawn C/N is arginine's, not methionine/phenylalanine) — the
    same honest stand-in mercaptans carries. ``o2`` is off every ledger (D-71), so spending it moves
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
    #: (``methional``/``phenylacetaldehyde``) + the decarboxylation ``CO2``, drawing the carbon from
    #: ``amino_acids`` (arginine), deaminating its nitrogen to ``N``. Touches those six and nothing
    #: else — ``o2`` is off every ledger; the C/N transfer closes exactly.
    touches = ("o2", "methional", "phenylacetaldehyde", "CO2", "amino_acids", "N")
    #: ``k_strecker``/``E_a_strecker``/``y_strecker_per_o2``/``f_methional`` are this Process's own
    #: (aging.yaml, D-75); ``K_amino_acids`` is the *shared* availability half-saturation (the same
    #: constant the mercaptan/reroute gates read); ``T_ref`` is shared with every Arrhenius rate.
    #: Their tiers cap the output tiers via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = (
        "k_strecker",
        "E_a_strecker",
        "y_strecker_per_o2",
        "f_methional",
        "K_amino_acids",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Wine-only slots (beer's amino-acid pool is not tracked, D-32): a hard no-op on any schema
        # without them, belt-and-suspenders to the wine-only wiring.
        if "amino_acids" not in schema or "o2" not in schema:
            return d
        o2 = float(y[schema.slice("o2")][0])
        aa = max(float(y[schema.slice("amino_acids")][0]), 0.0)
        # No O₂ or amino acids ⇒ no Strecker: reductive/amino-acid-free aging is byte-for-byte
        # the case without this Process. ``<= 0`` also absorbs solver undershoot (o2 < 0 ⇒ no draw).
        if o2 <= 0.0 or aa <= 0.0:
            return d
        # Smooth availability gate (the D-33 swap/reroute idiom): throttles the O₂ draw down to
        # 0 as the pool empties, so O₂/carbon/N vanish together; amino_acids never goes negative.
        gate = aa / (params["K_amino_acids"] + aa)  # in [0, 1)
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_strecker"], params["T_ref"])
        # This route's aa-gated SHARE of the O₂-depletion rate — a SMALL wine-only add-on (NOT in
        # the 5.0e-4 always-on anchor; substrate-gated, adds on top like SulfiteOxidation, D-72).
        r_o2 = params["k_strecker"] * f_t * o2 * gate  # g O2/L/h consumed by the Strecker route
        n_ald = params["y_strecker_per_o2"] * (r_o2 / M_O2)  # mol total Strecker aldehyde/L/h
        f_meth = params["f_methional"]  # mol fraction methional; the rest is phenylacetaldehyde
        meth_rate = f_meth * n_ald * M_METHIONAL  # g methional/L/h
        phenyl_rate = (1.0 - f_meth) * n_ald * M_PHENYLACETALDEHYDE  # g phenylacetaldehyde/L/h
        co2_rate = _CO2_PER_STRECKER_ALDEHYDE * n_ald * M_CO2  # 1 CO₂ per aldehyde (decarb)

        # Draw the product carbon (aldehydes + CO₂, all on-ledger) from amino_acids sized to
        # match, and deaminate the arginine nitrogen to N (the D-45 idiom + CO₂ decarb): carbon
        # out of amino_acids == carbon into products, and all arginine N lands in N (products are
        # N-free), so total_carbon and total_nitrogen both close to machine precision.
        product_carbon = (
            meth_rate * carbon_mass_fraction(_METHIONAL_SPECIES)
            + phenyl_rate * carbon_mass_fraction(_PHENYLACETALDEHYDE_SPECIES)
            + co2_rate * carbon_mass_fraction("CO2")
        )  # g C/L/h into the products
        c_aa = carbon_mass_fraction(AMINO_ACID_SPECIES)
        y_n = nitrogen_mass_fraction(AMINO_ACID_SPECIES)
        aa_mass = product_carbon / c_aa  # arginine mass consumed to supply that carbon

        d[schema.slice("o2")] = -r_o2  # this route's aa-gated O₂ share (off every ledger)
        d[schema.slice("methional")] = meth_rate
        d[schema.slice("phenylacetaldehyde")] = phenyl_rate
        d[schema.slice("CO2")] = co2_rate
        d[schema.slice("amino_acids")] = -aa_mass
        d[schema.slice("N")] = aa_mass * y_n  # DEAMINATION: arginine N → ammonium (D-45)
        return d
