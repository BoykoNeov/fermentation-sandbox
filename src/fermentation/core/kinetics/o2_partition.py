"""The always-on O2-depletion total, and the single fraction that splits it — decision **D-172**.

Wine and beer scavenge a dissolved-oxygen charge over weeks-to-months even with no dosed
sulfite, no grape phenolics and no oak: an *always-on* first-order sink on the ``o2`` pool.
D-71 shipped that as one constant. D-73/D-74 split it into an ethanol-oxidation route
(which makes acetaldehyde) and a phenolic-browning route (which makes A420), and D-141's
cascade re-summed it as an activation floor. Three entries, three bands — for one number.

**They were documented as one degree of freedom and drawn as three.** Every one of those
notes stated the identity ``k_ethanol_oxidation + k_browning_base = 5.0e-4``, and the
sampler drew each band independently, so the sum it asserted ranged over
``[2.4e-4, 1.2e-3]`` and equalled 5.0e-4 at exactly one point. D-171 §6 flagged it; this
module is the repair. The partition is now *parameterised* rather than *asserted*::

    k_ethanol_oxidation  ==  k_o2_depletion_total *      f_ethanol_o2_share
    k_browning_base      ==  k_o2_depletion_total * (1 - f_ethanol_o2_share)

so the identity holds by construction at every draw, in every oxidative set, and there is
nothing left for a test to assert about it — see :func:`o2_depletion_shares`.

**Why two products and not a subtraction.** ``total - ethanol`` is the more obvious way to
write "one degree of freedom", and it is wrong here for an arithmetic reason: at the shipped
nominals ``5.0e-4 - 2.0e-4`` is ``3.0000000000000003e-4``, not ``3.0e-4``. Both products are
exact in float64 (``5.0e-4 * 0.4 == 2.0e-4`` and ``5.0e-4 * 0.6 == 3.0e-4``), so the nominal
trajectory reproduces the pre-D-172 tree **bitwise** — which is the whole content of the
claim that this repair moved no number.

**And that requirement is invisible to the shipped suite — measured, not assumed.** Swapping
this function to the subtraction form leaves 361 tests GREEN (D-172 §9, arm A). The
``np.array_equal`` isolability pins do *not* catch it, because each compares two runs that
both use whichever formula is in the file, so a one-ULP shift moves them together and they
stay self-consistent. Exactness matters against the values the pre-repair tree produced, and
nothing in the suite compared against those. That is precisely why
``tests/test_o2_partition.py`` pins the arithmetic directly: without it, "simplify to
``total - total * f``" is a green-on-arrival change that silently un-does the bitwise claim.

**What this does NOT settle.** D-74's load-bearing claim is that browning is the *dominant*
share. That is now exactly ``f_ethanol_o2_share < 0.5`` — one parameter, one edge — and the
band's own construction rule licenses 0.6, so it still breaches, on 12.5 % of triangular
draws against the 55.20 % D-171 measured for the two independent bands. Closing it needs
that high edge narrowed to 0.5 on the strength of the author's own prose, which is the move
D-171 refused seven times. It is left open, and guarded at the nominal only.
"""

from __future__ import annotations

from collections.abc import Mapping

__all__ = ["o2_depletion_shares"]


def o2_depletion_shares(params: Mapping[str, float]) -> tuple[float, float]:
    """``(ethanol_share, browning_share)`` of the always-on O2-depletion rate, in 1/h.

    A pure function of the parameter map — the shared-gate idiom (D-31/D-100), so the
    direct set's two Processes and the cascade's H2O2 branch all recompute the same split
    independently instead of one of them owning it. The two returned constants are what
    ``k_ethanol_oxidation`` and ``k_browning_base`` were as YAML entries before D-172, at
    the same nominal values, bitwise.

    They sum to ``k_o2_depletion_total`` *mathematically* by construction; callers must not
    assume they do so *bitwise* for a drawn fraction, since ``t*f + t*(1-f)`` carries two
    roundings. The nominal is exact; an arbitrary draw is not guaranteed to be.
    """
    total = params["k_o2_depletion_total"]
    f = params["f_ethanol_o2_share"]
    return total * f, total * (1.0 - f)
