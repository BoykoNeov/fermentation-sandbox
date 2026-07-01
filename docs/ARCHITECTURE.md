# Architecture

This describes the implemented skeleton (Milestone 0) and how the layers fit
together. The original brief is [`FERMENTATION_SIM_HANDOFF.md`](FERMENTATION_SIM_HANDOFF.md);
the design decisions and where we deviated from it are in [`DECISIONS.md`](DECISIONS.md).

## Layering

Four layers with strictly one-directional dependencies — a lower layer never
imports a higher one:

```
  scenario / validation   declarative recipes, benchmark comparison, analysis
        │  consumes time-series; owns no physics
  ──────┼──────────────────────────────────────────────────────────
  runtime                 time-stepping (solve_ivp), events, ensembles
        │  integrates the core; knows nothing of UI
  ──────┼──────────────────────────────────────────────────────────
  domain core             state vector + Process objects that contribute rates
        │  pure, deterministic, no I/O, no global state, no randomness
  ──────┼──────────────────────────────────────────────────────────
  parameters / units      versioned data (value + provenance + tier); conversions
```

Package map:

| Layer | Package | Key types |
|-------|---------|-----------|
| parameters | `fermentation.parameters` | `Parameter`, `Provenance`, `Uncertainty`, `ParameterSet`, `load_parameters`, `default_data_dir` |
| units | `fermentation.units` | `brix_to_sg`, `sg_to_plato`, `abv_from_ethanol`, … |
| core | `fermentation.core` | `Tier`, `StateSchema`, `VarSpec`, `StateVector`, `Process`, `ProcessSet`, `RateModifier`, `Medium`, `MEDIA`, `get_medium`, `wine_schema`, `beer_schema`; `chemistry` (molar masses, carbon fractions, Gay-Lussac split, `sugar_species`); `acidbase` (`solve_ph`, `ph_of_state`, `titratable_acidity`, `ph_tier`, charge balance); `kinetics` (`GrowthNitrogenLimited`, `SugarUptakeToEthanolCO2`, `EthanolInhibition`) |
| runtime | `fermentation.runtime` | `simulate`, `Trajectory` |
| scenario | `fermentation.scenario` | `Scenario`, `TemperaturePoint`, `Intervention`, `compile_scenario`, `CompiledScenario` |
| validation | `fermentation.validation` | `assert_conserved`, `assert_nonnegative`, `total_carbon`, `total_nitrogen`, `total_mass`, `BenchmarkSpec`, `ReferenceSeries`, `compare_series` |
| analysis | `fermentation.analysis` | `ph_series`, `titratable_acidity_series` (top-layer observables over a `Trajectory`) |

## The core

### State vector
A single contiguous `float64` numpy array (`StateSchema` maps names → index
slices). Keeping it a plain array is what lets `scipy.integrate.solve_ivp` drive
it efficiently. Variables can be vectors: `S` (sugar) is one slot for wine and
three (glucose/maltose/maltotriose) for beer, so beer is an *addition* not a
rewrite.

Initial v1 contents: viable biomass `X`, sugar(s) `S`, ethanol `E`, yeast
assimilable nitrogen `N`, temperature `T`, evolved CO₂ (the experimentally
measurable proxy — a primary validation channel). Plus *produced-only* pools that
are 0 at pitch and only accumulate during fermentation: inactivated biomass
`X_dead` (D-13), and the realised-yield byproduct sinks `Gly` (glycerol) and
`Byp` (lumped minor byproducts) (D-16). These declare a `VarSpec.default` so
`pack` fills them when omitted, while substrate/condition vars stay required.

### Process
A `Process` contributes to `d(state)/dt`. It declares `name`, `tier`, and the
state variables it `touches`, and implements
`derivatives(t, y, schema, params) -> contribution`. `ProcessSet` sums the active
Processes (that sum is what `solve_ivp` integrates) and, crucially, **derives the
output tier** of each variable as the lowest tier among the Processes that touch
it. Toggling a speculative Process off leaves the validated core intact.

In `strict=True` mode, `ProcessSet` verifies every Process only writes to the
variables it declared — a cheap guard used in tests.

### Rate modifiers
Some mechanisms *scale* an existing flux rather than *add* a new one — ethanol
inhibition slows fermentative uptake; Arrhenius temperature scales every rate
constant. Because `ProcessSet` **sums** Processes, these cannot be summed Processes
(they would add to a derivative, not multiply it). A `RateModifier` declares which
Processes it `modifies` (by name) and returns a scalar `factor`; `ProcessSet`
multiplies that factor onto the *entire contribution vector* of each targeted
Process before summing. Scaling a conserving Process's whole contribution by one
scalar preserves its mass/atom balances, so a modifier never breaks conservation,
and the `touches` contract still holds (scaling zeros stays zero). Modifiers are
enabled/disabled and feed `tier_of` exactly like Processes — a speculative modifier
drags the tier of the variables its target touches down to speculative.
`EthanolInhibition` is the first modifier; `ArrheniusTemperature` reuses the hook,
parameterised *per rate* (one instance per Process, each its own activation energy)
and reference-anchored so `f = 1` at `T_ref`. Stacked modifiers on one Process (e.g.
inhibition × Arrhenius on uptake) compose to a single scalar, so conservation still
holds. (See DECISIONS #10, #11.)

## Confidence tiers

`Tier` is an ordered enum (`VALIDATED > PLAUSIBLE > SPECULATIVE`). The trust of a
combination is the `min` (`Tier.combine`). Tiers are a property of *Processes and
parameters*, not of the raw floats flowing through the solver; an output's tier
is computed at the analysis boundary. This satisfies "the tier travels to every
output" without polluting the integration hot loop. (See DECISIONS #1.)

No Process or parameter is `VALIDATED` yet: that tier is reserved for checks against
independent *measured* time-series, which do not exist yet (DECISIONS #C, #17).
Passing the §2.2 benchmarks earns `PLAUSIBLE` — sound forms, sourced parameters,
reproduces the keystone model — not `VALIDATED`.

## Parameters with provenance

Every kinetic/physical constant is a `Parameter` requiring value, units, tier,
an `Uncertainty` range, and `Provenance` (source + measurement conditions). The
Pydantic models reject any entry missing these, so "no magic numbers" is a hard
load-time guarantee. Parameters live in YAML under
`src/fermentation/parameters/data/`; strain-specific overlays merge on top of
generic defaults (`ParameterSet.merge`).

## Units boundary

Canonical internal units: concentration **g/L** (≡ SI kg/m³), temperature **K**,
time **hours**. Industry units (°Brix, SG, °Plato, %ABV, °C, days) appear only on
the far side of `fermentation.units`. (See DECISIONS #3.)

## Runtime

`simulate(process_set, params, y0, t_span)` wraps `solve_ivp` with an implicit
adaptive method (BDF by default — fermentation is stiff) and returns a
`Trajectory` carrying the time grid, the state history, and the derived tier map.
The event-driven loop (interventions, phase switching) and the stochastic
ensemble wrapper layer on top of this without changing the core.

## Scenarios as data

A `Scenario` (schema-validated YAML/JSON, **not** a custom DSL) declares initial
composition, organism/strain, temperature schedule, vessel, and a timeline of
interventions. No physics lives here, which keeps sweeps, Monte Carlo, and
cross-beverage reuse trivial.

## Media and the compile seam

A **`Medium`** (`fermentation.core.media`) names a beverage family and fixes its
`StateSchema` plus the Processes that act on it. `wine_schema()` has one sugar
slot; `beer_schema()` has three (`glucose`/`maltose`/`maltotriose`, in uptake
order). The `MEDIA` registry maps name → `Medium`; processes are empty until the
M1 kinetics land, so a compiled medium currently integrates to a constant
baseline.

**`compile_scenario(scenario)`** (`fermentation.scenario.compile`) is the
scenario→core seam and the *only* place industry units cross into canonical ones
(°Brix → g/L, °C → K, days → hours). It validates the `scenario.initial`
vocabulary (per-medium allowed keys, non-negativity, required fields), seeds the
initial temperature from the schedule, loads `<medium>_<strain>.yaml`, assembles
the medium's `ProcessSet`, and returns a `CompiledScenario` record (`y0`,
`process_set`, `parameters` + resolved `param_values`, `schema`, `t_span_h`) that
drops straight into `simulate`. Beer's three sugars are supplied explicitly rather
than split from a single OG — that wort spectrum is a provenance-backed parameter,
not a magic constant in the seam. (See DECISIONS #7.)

## Validation

Two disciplines, both as code:
- **Conservation invariants** — `assert_conserved` / `assert_nonnegative` take a
  model-supplied conserved-quantity function and check it holds to tolerance along
  a trajectory. The chemistry-specific quantities are built by `total_carbon`,
  `total_nitrogen`, and `total_mass`, which weight each state variable using the
  shared stoichiometry in `fermentation.core.chemistry` — so a check can never
  disagree with the kinetics it audits. Carbon and nitrogen are the rigorous atom
  balances (the biomass C/N fraction is a passed-in Parameter); mass is scoped to
  the abiotic `S + E + CO₂` conversion. (See DECISIONS #8.)
- **Benchmark curves** — the §2.2 acceptance criteria are encoded as
  `BenchmarkSpec` data now; the `tests/benchmarks/` tests are skipped until the
  kinetics exist. `ReferenceSeries` + `compare_series` (RMSE/MAE) are the seam
  for scoring against *real* measured datasets when we obtain them.

## pH as a derived pure function (acid state + charge balance, D-18)

pH is **not** an integrated state — there is no `dpH/dt`. Like `total_carbon` and ABV it
is an instantaneous, pure algebraic function of state: `fermentation.core.acidbase` solves
electroneutrality `Σ charge = 0` for `[H⁺]` (a 1-D monotonic root-find in pH-space, via
`brentq`) given the charge-active acid concentrations and a pKa set, and reports
`pH = −log₁₀[H⁺]`. Building it as a full proton balance (not a tracked-pH approximation)
is what makes the Tier-2 couplings — MLF deacidification, SO₂ speciation — *emerge* rather
than be scripted (DECISIONS #18).

- **Acid state (wine only).** `wine_schema` appends four slots: `tartaric`, `malic`,
  `lactic` (diprotic/diprotic/monoprotic wine acids, carbon-weighted in `total_carbon` for
  a future MLF Process) and `cation_charge`, the net strong-cation charge density (mol⁺/L,
  K⁺-dominant). The cation is **mandatory** (weak acids alone give pH ≈ 2.3 vs a real ~3.3)
  and **back-solved from a measured `initial_ph`** at the compile seam (inverse anchoring),
  so the model predicts pH *changes*, not absolute initial pH. `beer_schema` is untouched —
  beer's acid system is deferred. The slots default to 0, so acid-free scenarios are inert
  and the validated core is unaffected (prime directive #3).
- **`Byp` include-by-reading.** The balance reads the existing `Byp` pool as a
  succinic-equivalent acid — zero new carbon, so `total_carbon` is unchanged and the D-16
  double-count is closed.
- **The observable layer.** Scalar `ph_of_state` / `titratable_acidity` are pure and live
  in core; the trajectory-series helpers (`ph_series`, `titratable_acidity_series`,
  `molecular_so2_series`) need `Trajectory`, so they sit one layer up in the new top-layer
  `fermentation.analysis` — mirroring how `units` provides scalar conversions and benchmarks
  map ABV over a series. Tier is reported via `acidbase.ph_tier` (computed explicitly as
  `plausible`, never the `VALIDATED` default of the inert acid slots).
- **SO₂ speciation = the first pH consumer, readout-only (D-22, D-28).** `wine_schema` appends
  a slot `so2_total` (total SO₂ as g/L SO₂-equivalent, dosed via `so2_total_mgl`, conserved/
  inert). `acidbase.speciate_so2` solves pH from the organic acids, then (D-28) splits the
  total into acetaldehyde-**bound** vs **free** via the bisulfite binding equilibrium
  (`bound_so2_molar` solves `(A−x)(C−x)β − Kx = 0`; new `bisulfite_fraction`), and returns the
  **molecular** (antimicrobial) fraction `1/(1+10^(pH−pKa₁))` (sulfurous pKa₁ 1.81) of *free* —
  the D-18 coupling *emerging*, not scripted. It is **readout-only**: SO₂ is kept out of the
  charge balance (the inverse anchoring makes in-balance vs readout identical at t=0) and out
  of titratable acidity (OIV excludes it), and is carbon-free — so dosing it leaves pH and
  `total_carbon` byte-for-byte (an isolability test pins this). At acetaldehyde=0 the split
  collapses to D-22 exactly (`free == total`). The lone RHS consumer is the MLF antimicrobial
  gate, which reads the *derived* free-molecular SO₂ (bound SO₂ is not antimicrobial), so the
  early acetaldehyde peak transiently sequesters SO₂ and relaxes suppression — an emergent
  competition. The bound-acetaldehyde-protected-from-ADH feedback stays deferred (readout-only).
- **Acetaldehyde** (`core/kinetics/acetaldehyde.py`, decision D-27) — the obligate main-
  pathway intermediate, modelled as a transient **ethanol-carbon buffer**: flux-linked
  production *borrows* a C2 slice of ethanol and viable-yeast-gated reduction *returns* it
  (both mole-for-mole C2→C2). It de-lumps the uptake Process's single sugar→ethanol step
  rather than adding a parallel pathway, so carbon closes touching neither `S` nor `CO2` and
  the `E` endpoint (hence the §2.2 benchmarks) is preserved to relative ~1e-8. The early
  produce-then-reabsorb peak emerges; a crash strands it (the D-26 live-yeast-gating shape).

## Testing & quality gates

`uv run pytest` (unit + integration + conservation; benchmarks skipped),
`uv run ruff check .`, `uv run mypy` (strict on `src`). CI runs all three on
Python 3.13 and 3.14.
