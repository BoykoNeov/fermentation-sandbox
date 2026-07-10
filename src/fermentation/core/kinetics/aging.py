"""Aging chemistry — the slow, post-fermentation "years" axis (§4.1, decision D-69).

The first Tier-3 aging Process, opened in D-68 and built here. Unlike the Milestone-2
byproduct Processes (which *produce* aroma pools during active fermentation), an aging
Process acts on a *finished* wine/beer over months-to-years: fermentation is done, the
sugar is gone, the yeast racked or crashed, and the chemistry that remains is spontaneous
(hydrolysis, oxidation, condensation), not metabolic. This module holds the aging
Processes; :class:`EsterHydrolysis` is the first and :class:`OxidativeAcetaldehyde` the
second.

**The oxidative sub-axis (D-71).** :class:`OxidativeAcetaldehyde` opens the *oxidative* half of
the aging axis on a **dissolved-O₂ pool** (``o2``, a new carbon-free state slot, off every
conservation ledger like ``h2s``/``iso_alpha``). O₂ — not ethanol — is the rate-limiting
reactant (ethanol sits at ~100 g/L, effectively constant across aging), so the pool is *the*
substrate that bounds oxidation: acetaldehyde saturates as the O₂ charge is spent, the
bottle-aging reality a first-order-in-ethanol rate could never reproduce (it would grow
unbounded). O₂ enters via a dedicated ``add_oxygen`` dosing verb (one dose = a bottle's
ingress; repeated = micro-oxygenation / barrel), and a ``begin_aging`` run with **no** O₂ dosed
is purely *reductive* aging (screwcap/inert) — byte-for-byte the :class:`EsterHydrolysis`-only
aging, since the Process contributes exactly zero at ``o2 = 0``. Because oxidative aging is
fundamentally a competition for a finite O₂ budget, the ``o2`` pool is the shared substrate the
whole future oxidative sub-axis (phenolic browning, Strecker degradation, SO₂ consumption) will
draw down — this Process claims only a speculative *yield* of it (``y_acetaldehyde_per_o2``,
below the mechanistic max), leaving the remainder as unmodeled oxidative sinks (D-71).

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

from fermentation.core.chemistry import (
    M_ACETALDEHYDE,
    M_ETHANOL,
    M_O2,
    carbon_mass_fraction,
)
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

    ``d(o2)/dt = −r_O2`` with ``r_O2 = k_ethanol_oxidation · f(T) · [O2]`` (first-order in the
    dissolved-O₂ pool, ``f(T) = arrhenius_factor(T, E_a_ethanol_oxidation, T_ref)`` the sourced
    warmer-oxidises-faster factor). The oxidised carbon lands as acetaldehyde at a molar **yield**::

        d(acetaldehyde)/dt = +y_acetaldehyde_per_o2 · (r_O2 / M_O2) · M_acetaldehyde
        d(E)/dt            = −d(acetaldehyde)/dt · M_ethanol / M_acetaldehyde

    **O₂, not ethanol, is the rate-limiting reactant** (the D-71 design crux). Ethanol sits at
    ~100 g/L, essentially constant across aging, so a rate first-order in *ethanol* would be a
    constant rate in disguise — acetaldehyde rising linearly and **unbounded**, pinning the kinetic
    limit on the wrong species. Making the rate first-order in the finite ``o2`` pool instead gives
    the correct **saturating** behaviour: as the O₂ charge is consumed the pool decays toward zero
    and acetaldehyde plateaus, the bottle-aging reality. (Mechanistically the real path is *coupled*
    oxidation — O₂ oxidises o-diphenols → quinones + H₂O₂, then H₂O₂ oxidises ethanol → acetaldehyde,
    Wildenradt & Singleton 1974; the phenolic catalyst is folded into ``k_ethanol_oxidation`` in
    v1, a documented lump since no general phenol pool is tracked.)

    **Carbon — the clean reverse of the D-27 reduction.** Ethanol and acetaldehyde are both C2, so
    the ``E → acetaldehyde`` transfer is mole-for-mole; weighting the ethanol debit by
    ``M_ethanol/M_acetaldehyde`` (:data:`_ETHANOL_PER_ACETALDEHYDE`) makes the carbon that leaves
    ``E`` land exactly in ``acetaldehyde``, so ``total_carbon`` closes to machine precision. This is
    the mirror of :class:`~fermentation.core.kinetics.acetaldehyde.AcetaldehydeReduction`
    (acetaldehyde → ethanol), and during aging that reduction is **inert** — it is gated on *viable*
    ``X``, which is 0 in a racked/finished wine — so oxidation does not fight it: the acetaldehyde
    accumulates, correctly (a live-yeast ferment would instead reduce it straight back).

    **The O₂ pool is off every ledger, and only a *yield* of it becomes acetaldehyde.** ``o2`` is
    carbon-free (``total_carbon``/``total_mass``/``total_nitrogen`` weight only their named pools,
    so ``o2`` contributes 0 to each, like ``h2s``/``iso_alpha``). The Process consumes the **whole**
    O₂ flux ``r_O2`` but routes only ``y_acetaldehyde_per_o2`` mol/mol into acetaldehyde — the
    remainder is the oxidative power spent on **unmodeled sinks** (phenolic browning, Strecker,
    SO₂ consumption) the future O₂ sub-axis will claim. Because O₂ carries no carbon, "spending" it
    without tracking every product is not a conservation violation — the carbon that *does* move
    (into acetaldehyde) is borrowed carbon-exactly from ``E``.

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
    #: ``k_ethanol_oxidation``/``E_a_ethanol_oxidation``/``y_acetaldehyde_per_o2`` are this Process's
    #: own (aging.yaml, D-71); ``T_ref`` is shared with every Arrhenius rate. Their tiers cap the
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
        r_o2 = params["k_ethanol_oxidation"] * f_t * o2  # g O2/L/h consumed (first-order in O2)
        # A molar yield of the consumed O₂ becomes acetaldehyde; the rest is unmodeled oxidative
        # sinks (the future O₂ sub-axis, D-71). moles O2 = r_o2/M_O2 ⇒ acetaldehyde mass rate below.
        acet_rate = params["y_acetaldehyde_per_o2"] * (r_o2 / M_O2) * M_ACETALDEHYDE  # g/L/h
        d[schema.slice("o2")] = -r_o2  # the whole O₂ flux is consumed (off every ledger)
        d[schema.slice("acetaldehyde")] = acet_rate
        # Carbon-exact C2 borrow from ethanol (the D-27 reduction reversed). No clamp needed: during
        # aging E ~ 100 g/L and acet_rate is trace, so this never drives E negative.
        d[schema.slice("E")] = -acet_rate * _ETHANOL_PER_ACETALDEHYDE
        return d
