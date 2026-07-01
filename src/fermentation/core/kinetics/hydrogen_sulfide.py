"""Hydrogen sulfide (H₂S) — the low-nitrogen sulfidic off-aroma, as a produced-only pool.

The §3.2 aroma beat after the SO₂ free/bound split (decision D-29), and the
**accounting-easiest** metabolite yet: H₂S ("rotten egg", sensory threshold ~1–2 µg/L)
carries no carbon, so — like free SO₂ (D-22) — it sits on *no* conservation ledger and its
production draws from nothing the model tracks (there is no sulfate/sulfur state; the sulfur
is simply released).

**The mechanism (inverse-nitrogen gate).** Yeast reduces sulfate to sulfide via the
sulfate-reduction sequence during fermentation. Normally the sulfide is fixed onto
nitrogen skeletons (O-acetylserine / O-acetylhomoserine) to build cysteine and methionine;
when yeast-assimilable nitrogen (YAN) runs low, there is no nitrogen acceptor for the
sulfide, so it is excreted as H₂S. Production is therefore **de-repressed at low nitrogen** —
the opposite of the Ehrlich fusel gate (which *needs* amino-acid nitrogen, ``N/(K_n+N)``,
D-19). Modelled as the fermentative-flux shape times a decreasing inverse-N gate:

    d(h2s)/dt = k_h2s · X·S_total/(K_sugar_uptake + S_total) · K_h2s_n/(K_h2s_n + N)

* **Flux-linked** (shares ``K_sugar_uptake`` with the uptake Process), so H₂S is made only
  during active fermentation and stops at dryness — the sulfate-reduction machinery runs
  while the cell ferments. This is the :class:`~fermentation.core.kinetics.vicinal_diketones\
  .AcetolactateExcretion` / :class:`~fermentation.core.kinetics.acetaldehyde.\
  AcetaldehydeProduction` producer shape.
* **Inverse-N gate** ``K_h2s_n/(K_h2s_n + N)``: ~0 when N is replete, → 1 as N → 0. Its
  half-saturation ``K_h2s_n`` is a **separate parameter on the YAN scale** (~0.1 g/L), *not*
  the growth ``K_n`` (0.0088 g/L) — reusing the growth constant would make a razor-edge gate
  that only opens in a thin sliver at near-zero N. The YAN-scale constant instead makes the
  gate a smooth, physiologically-relevant repression across the must's nitrogen range
  (H₂S-management practice targets YAN ≳ 140–150 mg/L; Ugliano/Henschke).

**Temperature-flat (v1 simplification), like the excretion / acetaldehyde producers.** Real
H₂S has a strain- and temperature-dependence, but it is not clean or benchmarked, so v1
holds production temperature-independent and documents it. No Arrhenius factor.

**Carbon / conservation — nothing to close.** H₂S is carbon-free (registered with 0 carbon
in :mod:`fermentation.core.chemistry`), so it is absent from ``total_carbon`` and its
production perturbs no carbon/nitrogen/mass balance. The Process **touches only ``h2s``** and
*reads* ``X``/``S``/``N`` from state without writing them, so it is the most isolable beat in
the model: disabling it leaves the **derivative (RHS) of every other state column byte-for-byte
identical** (verified exactly — no other column's contribution changes). The *integrated*
trajectory then differs by only a ~1e-7 relative amount, which is a **pure adaptive-solver mesh
artifact** — adding the ``h2s`` equation shifts the error-controlled step selection, not any
physical pathway. This is cleaner than the acetaldehyde buffer (D-27), whose transient ``E``
dip feeds a *genuine* second-order ``E``→viability perturbation on top of the mesh effect;
here nothing downstream reads ``h2s``, so there is no physical coupling at all.

**Tier — no headline consequence.** The rate constants are order-of-magnitude estimates, so
the Process is **speculative** and parameter-tier propagation (D-1) caps the ``h2s`` output
tier at speculative. Unlike the diacetyl decarboxylation (which writes the shared ``CO2``
slot, D-26) and the acetaldehyde production (which writes ``E``, D-27), this Process writes a
**fresh pool nothing else reads**, so *no other column's* structural tier drops.

**Scope (v1) / honest caveats.**

* **Produced-only** — there is no CO₂-stripping volatilization sink yet, so the ``h2s`` pool
  is *cumulative H₂S produced*, which **overstates residual** H₂S: real fermentation sweeps
  most H₂S out with the CO₂ stream, leaving µg/L residuals. A stripping sink is the deferred
  follow-up, exactly the ester D-19 (produced-only) → D-20 (Henry's-law sink) precedent.
* **The cross-must YAN lever is muted by an upstream model gap.** The defining real behaviour
  — a low-YAN must makes far more H₂S than a high-YAN one — is only *partially* reproduced,
  because the nitrogen model strips ``N`` to ~0 early (by ~day 1.3) *regardless of dose* (no
  residual-N/proline floor). Once N = 0 the inverse gate is ~1 for the rest of the ferment
  for **every** must, so cumulative H₂S differs little between a low- and high-YAN must. What
  *does* emerge cleanly and is the acceptance anchor: **within a run, the H₂S production rate
  ramps up as N depletes** — the gate mechanism itself. The cross-must lever becomes real only
  once the deferred residual-N floor lands (see decision D-29 and the D-23 nitrogen-gap note).
* Yeast-pathway (sulfate-reduction) H₂S only; other sulfides / mercaptans and the
  copper-binding chemistry are out of scope.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.kinetics.carbon_routing import fermentative_flux_shape
from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier


class HydrogenSulfideProduction(Process):
    """H₂S release — flux-linked, de-repressed at low nitrogen; a carbon-free produced pool.

    ``d(h2s)/dt = k_h2s · X·S_total/(K_sugar_uptake + S_total) · K_h2s_n/(K_h2s_n + N)``.
    The fermentative-flux shape ties production to active fermentation (it stops at dryness);
    the inverse-N gate ``K_h2s_n/(K_h2s_n + N)`` de-represses production as YAN runs out —
    the low-nitrogen mechanism (the sulfate-reduction sequence outrunning the nitrogen-
    dependent assimilation of sulfide). Held **temperature-flat** (a documented v1
    simplification, like the α-acetolactate excretion, D-26).

    Touches **only ``h2s``** (carbon-free ⇒ no draw from ``S``, no ``CO2`` term, on no
    conservation ledger); reads ``X``/``S``/``N`` from state without writing them. So the
    Process is unconditionally isolable — disabling it leaves the RHS of every other column
    byte-for-byte identical (no ``h2s`` consumer exists to feed anything back); the integrated
    trajectory differs only by a ~1e-7 adaptive-solver mesh artifact, not a physical coupling.
    ``N`` is clamped ≥ 0 against solver undershoot; the flux shape clamps ``X``/``S``. Tier
    **speculative** (rate magnitude estimate).
    """

    name = "hydrogen_sulfide_production"
    tier = Tier.SPECULATIVE
    touches = ("h2s",)
    #: ``K_sugar_uptake`` is shared with the fermentative-uptake flux this tracks; ``k_h2s``
    #: sets the release magnitude and ``K_h2s_n`` the inverse-nitrogen half-saturation (on the
    #: YAN scale, distinct from the growth ``K_n``). Their tiers cap ``h2s``'s output tier via
    #: parameter-tier propagation (D-1).
    reads: tuple[str, ...] = ("k_h2s", "K_sugar_uptake", "K_h2s_n")

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        n = max(float(y[schema.slice("N")][0]), 0.0)
        k_n = params["K_h2s_n"]
        nitrogen_gate = k_n / (k_n + n)  # ~0 when N replete, -> 1 as N -> 0 (de-repression)
        d[schema.slice("h2s")] = params["k_h2s"] * flux * nitrogen_gate
        return d
