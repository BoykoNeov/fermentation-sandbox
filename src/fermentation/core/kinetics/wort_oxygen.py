"""Wort oxygen — the dissolved O₂ a pitched beer carries, and the yeast that strips it (D-213).

Beer's ``o2`` slot has existed since D-71, but **nothing ever seeded it**: a beer ferment ran at
exactly 0.000 mg/L dissolved oxygen from pitch to package. That is the one thing every real
brewery deliberately puts *into* the wort — it is aerated on the way to the fermenter precisely
so the yeast can rebuild the sterols it needs to divide. D-212 §3a found the gap while pricing an
early-acetic candidate, could not use it, and recorded it as an unseeded driver; this module
closes it.

**What ships is the STATE and its removal — not a coupling to anything.** The sourced passage
behind the seed (*Craft Beers*, see ``o2_wort_aeration_beer``) says in the same sentence that more
aeration grows more yeast, ferments harder, leaves less nitrogen, drops pH faster and makes more
acetaldehyde, VDK precursors, esters and higher alcohols. **None of those six is built here**, and
the omission is deliberate rather than partial: beer's ``mu_max`` was calibrated one beat earlier
(D-211), and a second growth limitation would re-open that calibration on an O₂→growth relation
the corpus does not quantify. See D-213 for the scope decision and the owner's call behind it.

**Consequence, stated plainly because it is the honest headline: this Process is INERT in beer's
default set.** The three Processes that consume ``o2`` — :class:`OxidativeAcetaldehyde`,
:class:`PhenolicBrowning` and :class:`EllagitanninOxidation` — are all wired into the *aging*
tuple and disabled at the compile seam until ``begin_aging``. So seeding O₂ and removing it moves
the ``o2`` column and **nothing else**.

**That aging-gating is NOT sufficient on its own, and saying so was this beat's one real defect**
(D-213 §9). As first shipped, this Process stayed enabled *past* the breakpoint and competed with
those sinks for a dosed ``add_oxygen``, eating ~45 % of it — an **aged** beer was measurably
changed. ``begin_aging`` therefore **disables** this Process, the first thing that verb has ever
switched off; before the breakpoint the sinks are off and this is on, after it the reverse.
**Inertness holds because that disable makes it hold**, not because nothing could ever read the
pool. Do not remove the disable on the grounds that the consumers are aging-gated — that is
precisely the reasoning the measurement refuted.

**Said precisely, because the exact claim differs between the RHS and the trajectory** — the
distinction D-42's ``h2s`` docstring draws, and the one
[[feedback-pin-tolerance-vs-solver-tolerance]] exists to keep honest. The **derivative** of every
other state column is byte-for-byte unchanged: this Process writes ``o2`` and nothing else, at
every hour, verified rather than asserted. The **integrated** trajectory then differs by ~2.5e-7
relative (worst absolute 8.2e-6 g/L on ~33 g/L of ethanol), which is a **pure adaptive-solver mesh
artifact** — adding the ``o2`` equation shifts error-controlled step selection, not any physical
pathway — and sits an order of magnitude inside the BDF ``rtol`` of 1e-6. **Do not write
"byte-for-byte" of the trajectory**; nothing downstream reads ``o2`` while fermentation runs, so
there is no physical coupling at all, but the integrator's mesh is not a physical pathway and not
a bitwise guarantee either. What the beat
buys is (a) a medium that no longer claims a brewery ferments anaerobically from t=0, and (b) a
pool that is *already at zero* by the time aging could begin — which is what stops a later
``begin_aging`` on beer from silently oxidising against a phantom 6.75 mg/L that a seed without a
sink would have left lying there. D-212 §7 named that hazard; this is the half that defuses it.

**The rate law, and why it is NOT growth-coupled.**

    d(o2)/dt = −k_o2_uptake_beer · X · o2

First-order in the dissolved oxygen, proportional to the biomass **present**. The *fate* the
sources describe is sterol and unsaturated-fatty-acid synthesis for new membranes during cell
division (four independent texts; see the parameter's provenance), which invites the growth-linked
form ``y · dX/dt`` that :class:`AceticAcidOverflow` uses. **That form was considered and rejected
on timing.** The same sources say the oxygen goes during the **lag phase** — *"pitched brewing
yeast will normally take several hours to adapt … before growth begins"*, and *"the dissolved
oxygen in beer rapidly disappears"* — i.e. **before** growth, not in step with it. Sizing a
growth-coupled yield from yeast sterol content (~1 % of dry weight, ~12 mol O₂ per mol ergosterol
⇒ ~0.010 g O₂/g X) empties the pool at ~26 h against this model's own growth curve, which
contradicts the sourced timescale. Yeast absorbs oxygen and lays down sterol *ahead* of dividing,
so biomass **present** is the right driver and biomass **formed** is not.

**That choice also removes a correctness coupling rather than creating one.** Because this Process
does not recompute the growth rate, it must **not** be a target of beer's growth Arrhenius modifier
— unlike :class:`AceticAcidOverflow`, which recomputes growth and therefore must be scaled by the
factor growth actually runs at (D-183, and D-32's rule before it). A Process that reads a state
slot it does not re-derive has no such obligation, which is why the ``reads`` tuple below is one
parameter and not four.

**Conservation is trivial and that is worth one line.** ``o2`` is carbon-free and sits off every
ledger (like ``h2s``/``iso_alpha``/``A420``), so removing it moves nothing that must balance —
no carbon-weighting subtlety of the kind D-135's ``bound_methanethiol`` needed, and no counterpart
slot to book the oxygen into. The oxygen genuinely leaves the modelled system: its destination
(membrane sterol) is not a tracked pool, exactly as D-29's sulfur has no sulfate state.
``touches = ("o2",)`` and nothing else.

Tier **speculative**: the seed level is a printed norm (``plausible``) but the removal constant is
an author estimate, and :meth:`ProcessSet.tier_of` takes the lower of the two.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: The only slot this Process moves. Named once so the class, its ``touches`` and its tests
#: cannot drift apart.
O2_SLOT = "o2"


class WortOxygenUptake(Process):
    """Yeast strips the wort's dissolved O₂ during the lag phase (decision D-213).

    ``d(o2)/dt = −k_o2_uptake_beer · X · o2`` — see the module docstring for why the driver is
    biomass *present* rather than biomass *formed*, and for the statement that this is inert in
    beer's default set.

    **Clamped at zero from below, and that clamp is load-bearing rather than defensive.** The
    form is exponential decay, which cannot reach zero in finite time and cannot cross it for a
    non-negative pool — but ``solve_ivp`` probes states the trajectory never visits, and a
    negative ``o2`` handed to the aging sinks would turn a *consumer* into a *source*. The guard
    is the same one :class:`ClosureOxygenIngress` states for the opposite sign.
    """

    name = "wort_oxygen_uptake"
    tier = Tier.SPECULATIVE
    #: Dissolved O₂ only. NOT ``X`` — reading the biomass is not contributing to it, and adding
    #: it would make ``ProcessSet(strict=True)`` permit a write that would double-count yeast.
    #: There is no destination slot either: the oxygen's fate is membrane sterol, which this
    #: model does not track (see the module docstring on conservation).
    touches: tuple[str, ...] = (O2_SLOT,)
    #: One parameter, and the shortness is the point — see the module docstring. This Process
    #: does not re-derive the growth rate, so unlike :class:`AceticAcidOverflow` it declares no
    #: growth constants and is not a growth-Arrhenius target.
    reads: tuple[str, ...] = ("k_o2_uptake_beer",)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        if O2_SLOT not in schema:
            # Hard no-op on a medium without the slot. Both current media carry it, so this is
            # symmetry with the other O2 Processes rather than a reachable branch.
            return d
        o2 = float(y[schema.slice(O2_SLOT)][0])
        if o2 <= 0.0:
            # Exhausted (the normal state after the first few hours) or a solver probe below
            # zero: contribute exactly nothing rather than a negative draw.
            return d
        biomass = float(y[schema.slice("X")][0])
        if biomass <= 0.0:
            # No yeast, no uptake — an unpitched wort holds its oxygen.
            return d
        d[schema.slice(O2_SLOT)] = -params["k_o2_uptake_beer"] * biomass * o2
        return d
