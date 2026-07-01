"""Vicinal-diketone (VDK / diacetyl) pathway — the diacetyl rest, mechanistically.

The remaining §3.2 aroma beat after esters/fusels (decision D-26). Diacetyl
(2,3-butanedione, a buttery off-note) is the defining lager-quality parameter, and
unlike the monotone ester/fusel pools it is **produced then reabsorbed** — a
non-monotonic time course. Modelling that faithfully (the owner's call over the
simpler closure-only options) means the *real* three-step pathway, each step
carbon-closing on the existing ledger:

    sugar --excretion--> α-acetolactate --decarb--> diacetyl + CO2 --reduction--> 2,3-butanediol
             (D-26)        (C5 reservoir)  (C5→C4+C1)     (C4)         (C4→C4)      (flavourless)

**Why three pools, not two (the load-bearing modelling choice, decision D-26).**
The α-acetolactate *reservoir* is what makes the diacetyl rest a rest:

* :class:`AcetolactateExcretion` fills the reservoir **during active fermentation**
  (coupled to the fermentative flux; it stops at dryness), so the reservoir is full at
  the end of primary fermentation.
* :class:`AcetolactateDecarboxylation` (added in the decarb step) converts reservoir →
  diacetyl by a **spontaneous, non-enzymatic, strongly temperature-dependent** reaction
  that is **not gated on yeast** — so it keeps making diacetyl *after* fermentation, slowly,
  faster when warm. This is the rate-limiting, temperature-critical step (``E_a_decarb``
  held high) and the reason a rest takes days-to-weeks and a warm rest is faster.
* :class:`DiacetylReduction` (added last) is **fast, enzymatic, and gated on viable yeast**
  (``X``, not ``X_dead``) with **no flux term**, so it clears diacetyl as fast as it forms
  while live yeast is present but stops dead once the yeast is crashed / racked /
  ethanol-inactivated.

Together these make the defining behaviour *emerge* rather than being scripted: a warm
rest with live yeast empties the reservoir and cleans the beer up fast; crash or package
too early (reservoir still full, yeast gone) and diacetyl **rises** in the package with
nothing to reduce it. A two-pool model (diacetyl produced flux-linked, reduced by live
yeast) reproduces neither — its diacetyl generation dies with the sugar, so it cannot
strand a *rising* diacetyl and loses the temperature-criticality of the rest. Hence the
reservoir is load-bearing, not cosmetic.

**Carbon — routed through real species, closing on the existing ledger (decision D-26).**
The owner asked for something closer to reality than either a "reabsorbed-carbon
returns-to-sugar" stand-in or a carbon-unaccounted trace pool. Tracking the true
downstream product delivers that:

* Excretion draws α-acetolactate's carbon *out of ``S``* (via
  :func:`~fermentation.core.kinetics.carbon_routing.draw_carbon_from_sugar`, option a1,
  D-19), booked at the C5 α-acetolactate fraction. This stand-in is **better grounded**
  than the ester/fusel ones: α-acetolactate genuinely derives from pyruvate (sugar).
* Decarboxylation is a carbon-closing ``C5 → C4 + CO2`` step, exactly like malolactic
  ``malic → lactic + CO2`` (D-23) — carbon moves within the ledger, no draw.
* Reduction is a mole-for-mole ``C4 → C4`` transfer from ``diacetyl`` to ``butanediol``,
  like ``esters → esters_gas`` (D-20) — both pools weighted at their own carbon fraction.

So ``total_carbon`` (which weights all three pools, see
:mod:`fermentation.validation.conservation`) closes to **machine precision** through the
whole produce-then-reabsorb course. ``total_mass`` does *not* close: the oxidative decarb
consumes O2 and the reduction consumes NAD(P)H, both untracked — a small mass gap exactly
analogous to beer's hydrolysis water (D-8). Carbon is the invariant here.

**Isolability (prime directive #3).** The three Processes live in their own ``_VDK_
PROCESSES`` tuple (``fermentation.core.media``), so a ProcessSet built without them is the
prior core. Diacetyl is intrinsic yeast metabolism (not a dosed organism like MLF), so it
is wired into *both* media and runs on every default fermentation — like esters, turning
it on draws only a *trace* of sugar into the reservoir (α-acetolactate peaks ~mg/L, roughly
an order of magnitude below the ester draw — negligible on ``dS``), leaving
``dX``/``dE``/``dCO2``/``dN`` byte-for-byte until the decarb/reduction move that carbon on.

**Tiers.** ``E_a_decarb`` carries a **sourced ordering** (the α-acetolactate → diacetyl
conversion is non-enzymatic and accelerates with temperature — Haukeli & Lie 1978;
Krogerus 2013 review) with a speculative magnitude, mirroring the ester/fusel ``E_a``
orderings. Every rate constant is an order-of-magnitude estimate, so all three Processes
are **speculative**; parameter-tier propagation (D-1) caps the pool outputs at speculative
regardless. **Scope (v1):** yeast valine-pathway diacetyl only — MLF/citrate diacetyl
(*Oenococcus*) is deferred with the MLF-growth beat, so wine yeast-pathway diacetyl
understates real wine diacetyl. Acetoin is lumped into the terminal ``butanediol`` pool.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.chemistry import (
    M_ACETOLACTATE,
    M_BUTANEDIOL,
    M_CO2,
    M_DIACETYL,
    carbon_mass_fraction,
)
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import (
    draw_carbon_from_sugar,
    fermentative_flux_shape,
)
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The precursor species whose C5 formula carbon-accounts the ``acetolactate`` reservoir
#: (decision D-26). Its carbon mass fraction weights both the sugar draw here and the pool
#: in ``total_carbon`` — one chemistry source of truth, so draw and check cannot disagree.
_ACETOLACTATE_SPECIES = "alpha_acetolactate"


class AcetolactateExcretion(Process):
    """α-Acetolactate excretion — fills the VDK precursor reservoir during fermentation.

    ``d(acetolactate)/dt = k_acetolactate · X · S_total/(K_sugar_uptake + S_total)`` with
    the carbon drawn *out of ``S``* (booked as α-acetolactate, C5). Yeast overflows
    α-acetolactate from valine biosynthesis while it ferments, so production is tied to the
    biomass-catalysed sugar flux (sharing ``K_sugar_uptake`` with the uptake Process) and
    **stops at dryness** — leaving a full reservoir that the spontaneous decarboxylation
    then slowly converts to diacetyl (the rest; see :class:`AcetolactateDecarboxylation`).

    Held **temperature-flat** (no explicit Arrhenius factor) as a documented v1
    simplification: real α-acetolactate formation rises mildly with temperature, but the
    reservoir *size* is a weak lever on the rest — the load-bearing temperature dependence
    is the decarboxylation, not the excretion (decision D-26). Touches ``acetolactate`` and
    ``S`` only (never ``E``/``CO2``), so with the VDK pathway off the core is byte-for-byte
    and with it on only ``dS`` gains a *trace* negative term (α-acetolactate is ~mg/L).
    The carbon source is better grounded than the ester/fusel stand-ins — α-acetolactate
    genuinely derives from pyruvate (sugar). Tier **speculative** (rate magnitude estimate).
    """

    name = "acetolactate_excretion"
    tier = Tier.SPECULATIVE
    touches = ("acetolactate", "S")
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks;
    #: ``k_acetolactate`` sets the excretion magnitude. Their tiers cap ``acetolactate``'s
    #: output tier via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_acetolactate", "K_sugar_uptake")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        rate = params["k_acetolactate"] * flux
        d[schema.slice("acetolactate")] = rate
        draw_carbon_from_sugar(d, y, schema, rate * carbon_mass_fraction(_ACETOLACTATE_SPECIES))
        return d


class AcetolactateDecarboxylation(Process):
    """Spontaneous α-acetolactate → diacetyl + CO2 — the rate-limiting, T-critical step.

    ``d(acetolactate)/dt = −r·M_acetolactate``, ``d(diacetyl)/dt = +r·M_diacetyl``,
    ``d(CO2)/dt = +r·M_CO2`` with the molar turnover ``r = k_decarb · f(T) ·
    [acetolactate]/M_acetolactate`` and ``f(T) = arrhenius_factor(T, E_a_decarb, T_ref)``.
    The oxidative decarboxylation is C5 → C4 + CO2, so carbon closes mole-for-mole on the
    existing ledger exactly like malolactic ``malic → lactic + CO2`` (D-23); no sugar draw.

    **Non-enzymatic and NOT yeast-gated** — this is the whole point (decision D-26). The
    reaction proceeds outside the cell whether or not viable yeast is present, so the
    α-acetolactate reservoir keeps converting to diacetyl *after* fermentation ends. It is
    **first-order in α-acetolactate** and **strongly temperature-dependent** (``E_a_decarb``
    held high, above the reduction's ``E_a_reduction``), which makes it the rate-limiting
    step of VDK removal and the reason a diacetyl rest is temperature-critical: warm it up
    and the reservoir empties to diacetyl faster (Haukeli & Lie 1978; Krogerus 2013 review;
    the sourced ordering, magnitude speculative). ``acetolactate`` is clamped ≥ 0 so a
    solver undershoot cannot manufacture diacetyl. Mass carries a small gap (the oxidative
    decarb consumes untracked O2) — carbon is the invariant, as for beer's hydrolysis water
    (D-8). Tier **speculative** (rate magnitude estimate).
    """

    name = "acetolactate_decarboxylation"
    tier = Tier.SPECULATIVE
    touches = ("acetolactate", "diacetyl", "CO2")
    #: ``k_decarb`` sets the spontaneous conversion magnitude; ``E_a_decarb`` (the sourced,
    #: load-bearing temperature ordering) and ``T_ref`` set the temperature shape. Their
    #: tiers cap the diacetyl output tier via parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_decarb", "E_a_decarb", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        acetolactate = max(float(y[schema.slice("acetolactate")][0]), 0.0)
        if acetolactate <= 0.0:  # nothing in the reservoir to convert
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_decarb"], params["T_ref"])
        r = params["k_decarb"] * f_t * acetolactate / M_ACETOLACTATE  # molar turnover
        d[schema.slice("acetolactate")] = -r * M_ACETOLACTATE
        d[schema.slice("diacetyl")] = r * M_DIACETYL
        d[schema.slice("CO2")] = r * M_CO2
        return d


class DiacetylReduction(Process):
    """Enzymatic diacetyl → 2,3-butanediol by *viable* yeast — clears the buttery note.

    ``d(diacetyl)/dt = −L`` and ``d(butanediol)/dt = +L · M_butanediol/M_diacetyl`` with the
    mass loss ``L = k_reduction · X · f(T) · [diacetyl]`` and ``f(T) = arrhenius_factor(T,
    E_a_reduction, T_ref)``. The reduction is a mole-for-mole C4 → C4 transfer (diacetyl and
    butanediol both have four carbons), so weighting ``butanediol`` at its own carbon
    fraction keeps the transfer carbon-neutral (like ``esters → esters_gas``, D-20); no sugar
    draw. Mass carries a small gap (the reduction consumes untracked NAD(P)H) — carbon is
    the invariant. Acetoin, the intermediate, is lumped into the terminal ``butanediol`` pool.

    **Gated on VIABLE biomass ``X`` (not ``X_dead``), with NO fermentative-flux term
    (decision D-26).** These two choices are the whole game:

    * *Live-yeast gating* — reduction is enzymatic, so it stops the moment the yeast is
      crashed, racked or ethanol-inactivated (``X → X_dead``). That is what strands diacetyl
      when the beer is packaged too early: the reservoir keeps decarboxylating (that step is
      not yeast-gated) but nothing reduces the diacetyl it makes, so diacetyl **rises**.
    * *No flux term* — reduction must run during the *rest*, after sugar is gone (flux ≈ 0).
      Coupling it to the fermentative flux (as excretion is) would switch it off exactly when
      it is needed, and the diacetyl rest would never clear. So reduction reads only viable
      ``X`` and the diacetyl present.

    ``E_a_reduction`` is held **below** ``E_a_decarb`` so the spontaneous decarboxylation, not
    this reduction, is the temperature-critical rate-limiting step (the ordering that makes a
    warm rest work). ``diacetyl`` and ``X`` are clamped ≥ 0 against solver undershoot. Tier
    **speculative** (rate magnitude estimate).
    """

    name = "diacetyl_reduction"
    tier = Tier.SPECULATIVE
    touches = ("diacetyl", "butanediol")
    #: ``k_reduction`` sets the enzymatic reduction magnitude; ``E_a_reduction`` (held below
    #: ``E_a_decarb``) and ``T_ref`` set the temperature shape. Their tiers cap the butanediol
    #: output tier via parameter-tier propagation (D-1). Reads viable ``X`` from state.
    reads: tuple[str, ...] = ("k_reduction", "E_a_reduction", "T_ref")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        diacetyl = max(float(y[schema.slice("diacetyl")][0]), 0.0)
        if diacetyl <= 0.0:  # nothing to reduce
            return d
        x_viable = max(float(y[schema.slice("X")][0]), 0.0)
        if x_viable <= 0.0:  # no viable yeast ⇒ no reduction (diacetyl is stranded)
            return d
        temp = float(y[schema.slice("T")][0])
        f_t = arrhenius_factor(temp, params["E_a_reduction"], params["T_ref"])
        loss = params["k_reduction"] * x_viable * f_t * diacetyl  # mass loss of diacetyl
        d[schema.slice("diacetyl")] = -loss
        d[schema.slice("butanediol")] = loss * M_BUTANEDIOL / M_DIACETYL  # mole-for-mole C4→C4
        return d
