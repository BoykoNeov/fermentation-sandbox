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

**Scope / honest caveats.**

* **Residual vs cumulative (D-42).** The CO₂-stripping sink
  (:class:`HydrogenSulfideVolatilization`, decision D-42) now sweeps most H₂S out with the CO₂
  stream, so the ``h2s`` pool is the *residual* (dissolved, µg/L) H₂S reality shows, and
  ``h2s + h2s_gas`` is cumulative produced (the ester D-19 produced-only → D-20/D-21 Henry's-law
  sink precedent, but carbon-free). Dropping *just* the sink recovers this module's original
  D-29 produced-only ``h2s`` pool byte-for-byte (``h2s_gas`` stays 0).
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

from fermentation.core.kinetics.arrhenius import arrhenius_factor
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


class HydrogenSulfideVolatilization(Process):
    """CO₂-stripping loss of dissolved H₂S to the headspace — the Henry's-law sink (D-42).

    ``d(h2s)/dt = -k_h2s_volatil · X·S_total/(K_sugar_uptake+S_total) · f_gas(T) · f_part(T)
    · h2s`` and the equal-and-opposite ``+`` into ``h2s_gas``. The exact structural mirror of
    :class:`~fermentation.core.kinetics.byproducts.EsterVolatilization` (decisions D-20→D-21),
    but **carbon-free**, so it is *simpler*: neither pool is on any conservation ledger, so the
    transfer needs no weighting to stay neutral (contrast ``esters``→``esters_gas``, both weighted
    as ethyl acetate in ``total_carbon``). H₂S is far more volatile than ethyl acetate, so
    fermentation sweeps ~all of it out with the CO₂ stream, leaving the µg/L residual reality
    shows — which is exactly the D-29 produced-only overstatement this beat lifts.

    * ``f_gas(T) = arrhenius_factor(T, E_a_uptake, T_ref)`` — the **gas-flow** factor: the
      stripping rides the evolving-CO₂ stream, whose rate is the fermentative uptake flux scaled
      by the *same* ``E_a_uptake`` the uptake Process carries (the ester precedent).
    * ``f_part(T) = arrhenius_factor(T, dH_h2s_volatil, T_ref)`` — the **gas/liquid partition**
      (van't Hoff) factor for the H₂S Henry's-law constant, which *rises* with temperature
      (H₂S dissolution is exothermic ⇒ warmer ⇒ more volatile). ``dH_h2s_volatil`` ≈ 17.5 kJ/mol
      is the **sourced** partition enthalpy (Sander Henry's-law compilation, −d ln kH/d(1/T)
      ≈ 2100 K), a physical Q10 ≈ 1.3 — weaker than the ester's +45 kJ/mol, *not* a fitted lever.

    **First-order in the dissolved H₂S present, and stops when fermentation stops** (``flux → 0``
    at dryness), so all produced H₂S is co-temporal with a CO₂ stream that can strip it (our
    production is likewise flux-linked; the problematic post-fermentation / autolytic H₂S that
    persists *because* no CO₂ sweeps it is out of scope, like the ester sink's omission of slow
    passive evaporation after the cap goes on). Because both production and stripping share the
    flux, the **residual quasi-steady-state ``h2s_ss = k_h2s·gate / (k_h2s_volatil·f_gas·f_part)``
    has the flux cancel** — residual H₂S tracks the inverse-N gate and temperature, not the
    ferment speed: it **rises as ``N`` depletes** (the gate opens) then **freezes at dryness**
    (both terms gate off together). An emergent, honestly-flagged artifact: production is held
    temperature-flat (D-29) while stripping rises with T, so the model predicts *residual H₂S
    falls with a warmer ferment* — physically reasonable (warm ferments purge sulfide) but
    **unbenchmarked**, and reality is mixed (warmth also raises production / N-demand, held flat
    here); directional/speculative only.

    Touches ``h2s`` and ``h2s_gas`` only — never ``S``/``E``/``CO2``/``N`` — and since both are
    carbon-free and on no ledger, the sink is neutral on every conservation sum **by
    construction** (no closure test needed; the invariant is ``h2s + h2s_gas`` unchanged from the
    sink-off produced pool). ``h2s`` is clamped ≥ 0 so a solver undershoot cannot strip a negative
    pool into spurious gas. Tier **plausible** in form (CO₂-stripping by the evolving gas is
    well-understood Henry's-law physics), with speculative rate parameters that cap the pool
    outputs at speculative via parameter-tier propagation (D-1).
    """

    name = "hydrogen_sulfide_volatilization"
    tier = Tier.PLAUSIBLE
    touches = ("h2s", "h2s_gas")
    #: ``K_sugar_uptake``/``E_a_uptake`` are shared with the fermentative uptake whose CO₂ stream
    #: does the stripping (gas-flow factor); ``dH_h2s_volatil`` is the sourced H₂S Henry's-law
    #: partition enthalpy (gas/liquid factor); ``T_ref`` anchors both. ``k_h2s_volatil`` sets the
    #: absolute stripping magnitude. Their tiers cap the pool outputs via parameter-tier
    #: propagation (D-1). Medium-agnostic (one physical mechanism), so all live in the shared
    #: ``hydrogen_sulfide.yaml`` — no per-medium split (unlike the ester ``dH``, whose *synthesis*
    #: direction differs by beverage; H₂S has no such split).
    reads: tuple[str, ...] = (
        "k_h2s_volatil",
        "K_sugar_uptake",
        "E_a_uptake",
        "dH_h2s_volatil",
        "T_ref",
    )

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
        if flux <= 0.0:
            return d
        h2s_liquid = max(float(y[schema.slice("h2s")][0]), 0.0)
        if h2s_liquid <= 0.0:  # nothing dissolved to strip
            return d
        temp = float(y[schema.slice("T")][0])
        f_gas = arrhenius_factor(temp, params["E_a_uptake"], params["T_ref"])  # CO2 gas flow
        f_part = arrhenius_factor(temp, params["dH_h2s_volatil"], params["T_ref"])  # partition
        rate = params["k_h2s_volatil"] * flux * f_gas * f_part * h2s_liquid
        d[schema.slice("h2s")] = -rate
        d[schema.slice("h2s_gas")] = rate
        return d
