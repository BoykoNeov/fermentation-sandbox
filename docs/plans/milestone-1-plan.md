# Milestone 1 — Tier-1 validated core

> Status: **in progress**. Done: medium state schemas + the scenario→core
> `compile_scenario` seam (`fermentation.core.media`,
> `fermentation.scenario.compile`); the carbon/nitrogen/mass conservation quantity
> functions (`fermentation.validation.total_carbon`/`total_nitrogen`/`total_mass`
> over the shared `fermentation.core.chemistry` stoichiometry, decision D-8); the
> first kinetic Process, `GrowthNitrogenLimited` (`fermentation.core.kinetics`) —
> Monod growth co-limited by sugar and YAN, conserving carbon and nitrogen to
> machine precision. Next: sugar uptake → ethanol inhibition, then wire the set
> into the media and unskip the benchmarks.
> Goal: single-strain, isothermal, nitrogen-limited primary fermentation that
> passes the §2.2 wine **and** beer benchmarks (decision D-B).

## Definition of done

The skipped tests in `tests/benchmarks/test_milestone1.py` pass for real:

1. **Wine:** ~264 g/L sugar must (24 °Brix, density-corrected) at 20 °C ferments to
   dryness in **10–14 days**,
   with a visible lag → exponential → stationary biomass trajectory.
2. **Beer:** ~1.048 OG ale wort at 20 °C attenuates to ~1.010 in **5–7 days**.
3. **CO₂:** evolution rate rises to a peak then tails; its integral tracks sugar
   consumed (carbon balance) within ±5%.
4. Conservation (carbon, mass) holds to tolerance for every run; concentrations
   never go negative.

All of `pytest`, `ruff`, `mypy` stay green; new parameters carry real provenance.

## Model (coupled ODEs over {X, S, E, N, T})

Build as composable `Process` objects (tag tier honestly — most will start
`plausible` until validated against the benchmark curves, then promote):

- **Growth** — Monod/logistic on biomass, **nitrogen-limited** (cells grow while
  YAN lasts, then stop dividing but keep fermenting at declining rate). This is
  the textbook stuck/sluggish mechanism and is mandatory.
- **Sugar uptake → ethanol** — Gay-Lussac stoichiometry; realised yield
  ~0.46–0.48 g/g (already in the param store, `plausible`). Beer: sequential
  uptake over the `S` vector (glucose → maltose → maltotriose).
- **Ethanol inhibition** — growth/rate decline as `E` rises; viability collapses
  past a strain tolerance. Pick a literature functional form; tag the constant.
- **Temperature dependence** — Arrhenius on rate constants (isothermal for M1, but
  wire the hook); warmer = faster but earlier stress.
- **CO₂ evolution** — stoichiometric partner of sugar consumed; the primary
  measurable validation channel.

## Approach (test-driven)

1. ✅ Define the wine and beer `StateSchema`s and a `compile(Scenario) -> (y0,
   ProcessSet, params)` step (the scenario→core seam).
2. Write the carbon/mass conservation quantity functions for the real chemistry.
3. Implement Processes one at a time, each with its own unit test, behind the
   `strict=True` touch contract.
4. Source parameters (the time sink — see context doc) and replace placeholders.
5. Iterate functional forms + parameters until the benchmarks pass. Unskip them.

## Out of scope for M1

Byproducts, pH/acid system, SO₂, second organisms, aging, sensory — all Tier-2+.
Keep them out so the validated core stays small and protected.

## Risks

- **Parameter sourcing** dominates the effort (handoff §2.3): strain/condition
  specific, scattered, sometimes contradictory. Budget for reconciliation, not
  transcription.
- **Stiffness / step control** during the lag→exponential transition; lean on the
  implicit solver, don't fight it with explicit methods.
