"""Osmotic/substrate inhibition of the fermentative flux at very high sugar (D-192).

A grape must does not ferment faster and faster the sweeter it gets. Past roughly
200 g/L the fermentation slows, and by 600-650 g/L a concentrated must is
*"practically unfermentable"* (Ribéreau-Gayon et al., **Handbook of Enology** vol. 1,
§3.4.1). The wired core had no term of that kind at all: growth's Monod
``S/(K_s + S)`` and uptake's ``S/(K_sugar_uptake + S)`` are both saturated far below
200 g/L, so before this module a must got *marginally faster* per litre as it got
sweeter (peak volumetric rate +14.5 % from 97 → 376 g/L, measured) and an 881 g/L
must — Tokaji-Eszencia concentration, which the source calls unfermentable —
fermented all the way to **19.8 % ABV**.

Like :mod:`~fermentation.core.kinetics.inhibition` this scales an existing flux rather
than adding one, so it is a :class:`~fermentation.core.process.RateModifier`, not a
summed :class:`~fermentation.core.process.Process`.

**Form — one-sided, thresholded, asymptotic.** With total sugar ``S``::

    f(S) = 1                                        for S <= S_thr
    f(S) = 1 / (1 + ((S - S_thr)/K_osm)**n_osm)      for S >  S_thr

Three properties are load-bearing and each was chosen against a measured alternative
(decision **D-192**):

* **Exactly 1 below the threshold — structurally, by early return.** Not "small", not
  "within tolerance": the factor is the literal constant ``1.0``, so every must whose
  sugar never reaches ``S_thr`` is *byte-for-byte* the pre-D-192 core. This is D-129's
  Gate 1 applied to the substrate axis, and it is what makes the term safe to default
  on despite being speculative.
* **C¹-smooth at the threshold** for ``n_osm >= 2`` (``f'(S_thr) = 0``), so BDF meets no
  derivative kink where the brake engages — the same reason
  :class:`~fermentation.core.kinetics.inactivation.EthanolToleranceDeath` uses a
  one-sided *quadratic* rather than Schenk's arctan sign-gate. ``n_osm < 2`` is rejected
  at read time.
* **Asymptotic, never zero.** A hard wall at the "unfermentable" concentration would
  make a very sweet must an *absorbing state*: uptake exactly 0, growth nitrogen-capped,
  and nothing left that can remove sugar — the must could never ferment at all, ever.
  Real must at that concentration ferments; it just takes years (Tokay Aszu "ferment for
  from two to five years or longer" — *Applied Wine Chemistry and Technology*). So the
  factor decays without reaching zero and the model produces a glacial ferment, not a
  frozen one.

**Why the threshold sits at 300 g/L and not at the source's 200 g/L.** The Handbook's
own onset is ~200 g/L, which is *below every normal wine must* (24 °Brix loads at
245 g/L). But the model's growth and death constants come from Coleman, Fish & Block
2007, whose fits are validated over **265-300 g/L initial sugar** and contain **no
substrate-inhibition term**. Whatever substrate inhibition operates below 300 g/L is
therefore already absorbed into those fitted constants, and an explicit term there
would double-count it. Placing ``S_thr`` at the top of that envelope is what makes the
inertness above *structural* rather than a tolerance argument.

**What that forfeits, stated rather than implied.** The Handbook also says alcohol
production *"can be lower in a must containing 300 g/l than in another containing only
200 g/l"*. That statement is **not reproduced** by this module, and D-192 measured why:
reaching it requires a brake that still bites at low sugar (a Haldane term with
``K_i ≈ 17 g/L``, i.e. a 92 % brake at 200 g/L), which is not "inhibition at high sugar"
at all but a global rate cut — and it puts the Coleman reconstruction's RMSE at 170.9
against a 2.0 g/L threshold. The 200-vs-300 reading is refused on the keystone's
evidence, not overlooked.

**Site of action.** The Handbook's mechanism sentence is growth-mediated: *"an elevated
amount of sugar hinders yeast growth and decreases the maximum population.
Consequently, fermentation slows."* So the modifier scales growth as well as uptake,
and in wine it must therefore also name
:class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation` — the swap funds
biomass from the amino-acid pool at growth's *realised* rate, and a modifier that scaled
growth without it would let the swap refund carbon and nitrogen against an unmodified
draw (the coupling decision D-32 records and ``media.py`` warns about). Including growth
is also the choice that costs the headline: at 70 °Brix it *raises* the arrested ABV from
5.36 % to 5.59 %, so the "unfermentable" result is not being bought with it.

**Conservation is automatic.** The factor scales each targeted Process's entire
contribution by one scalar, so every balance those Processes respect is preserved — a
uniformly slower carbon-neutral flux is still carbon-neutral, and a uniformly slower
growth draw still removes exactly ``f_N``/``f_C`` per gram of new biomass.

**Wine only.** The anchors are grape-must numbers and beer's own literature supplies
none; beer is therefore inert *by not being wired*, which is a stronger claim than
inertness by parameter value (the heaviest wort in the suite is ~250 g/L, below the
threshold, but nothing in the parameter file would guarantee that for a hypothetical
barleywine). Wiring it into beer would additionally have to name
``OrganicAcidExcretion``, which re-derives the fermentative flux itself.

Tier: **speculative**. WHERE the brake engages is sourced twice over (Coleman's envelope
for the onset, the Handbook's 600-650 g/L for the far anchor); HOW SHARP it is is not —
``n_osm`` is underdetermined by those two anchors, and across its band the arrested ABV
of an 881 g/L must spans 3.8-5.4 %. That split is the same one decision D-129 drew for
the ethanol ceiling.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.kinetics.amino_acids import AminoAcidAssimilation
from fermentation.core.kinetics.growth import GrowthNitrogenLimited
from fermentation.core.kinetics.uptake import SugarUptakeToEthanolCO2
from fermentation.core.process import RateModifier
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Smallest exponent that keeps ``f`` C¹-smooth at the threshold. Below 2 the factor has
#: a corner at ``S_thr`` (``f'(S_thr) = -1/K_osm`` at ``n = 1``) — precisely the kind of
#: derivative discontinuity a BDF step probes across and pays for
#: (``num_jac`` straddling a jump was measured at D-182).
_MIN_EXPONENT = 2.0


class OsmoticSubstrateInhibition(RateModifier):
    """One-sided high-sugar brake on the fermentative flux and on growth.

    ``factor = 1`` for ``S <= S_osmotic_threshold``; above it
    ``1/(1 + ((S - S_osmotic_threshold)/K_osmotic_inhibition)**n_osmotic_inhibition)``.
    ``S`` is the **total** over all sugar slots. Multiplied onto
    :class:`~fermentation.core.kinetics.uptake.SugarUptakeToEthanolCO2`,
    :class:`~fermentation.core.kinetics.growth.GrowthNitrogenLimited` and
    :class:`~fermentation.core.kinetics.amino_acids.AminoAcidAssimilation` by
    :class:`~fermentation.core.process.ProcessSet`. See the module docstring for the
    form, the threshold's provenance and what it deliberately does not reproduce.
    """

    name = "osmotic_substrate_inhibition"
    tier = Tier.SPECULATIVE
    #: Reference each target by its ``name`` (rename-safe) rather than a bare literal.
    #: ``AminoAcidAssimilation`` rides along because it funds biomass at growth's
    #: *realised* rate — see the module docstring and decision D-32.
    modifies = (
        SugarUptakeToEthanolCO2.name,
        GrowthNitrogenLimited.name,
        AminoAcidAssimilation.name,
    )
    reads: tuple[str, ...] = (
        "S_osmotic_threshold",
        "K_osmotic_inhibition",
        "n_osmotic_inhibition",
    )

    def factor(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> float:
        # Clamp the total >= 0 first: a negative solver excursion on one sugar slot must
        # not read as *less* inhibition than an exhausted must, and could not anyway
        # reach the branch below — but the clamp keeps the two readings identical.
        s_total = max(float(y[schema.slice("S")].sum()), 0.0)
        over = s_total - params["S_osmotic_threshold"]
        if over <= 0.0:
            # EXACTLY 1, by early return: below the threshold this modifier is not an
            # approximation of the pre-D-192 core, it *is* the pre-D-192 core.
            return 1.0
        exponent = params["n_osmotic_inhibition"]
        if exponent < _MIN_EXPONENT:
            raise ValueError(
                f"n_osmotic_inhibition = {exponent:g} is below {_MIN_EXPONENT:g}; the "
                "factor would have a derivative corner where the brake engages. See "
                "fermentation.core.kinetics.osmotic for why the exponent is bounded."
            )
        return float(1.0 / (1.0 + (over / params["K_osmotic_inhibition"]) ** exponent))
