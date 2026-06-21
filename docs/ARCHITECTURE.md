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
| core | `fermentation.core` | `Tier`, `StateSchema`, `VarSpec`, `StateVector`, `Process`, `ProcessSet`, `Medium`, `MEDIA`, `get_medium`, `wine_schema`, `beer_schema`; `chemistry` (molar masses, carbon fractions, Gay-Lussac split, `sugar_species`); `kinetics` (`GrowthNitrogenLimited`) |
| runtime | `fermentation.runtime` | `simulate`, `Trajectory` |
| scenario | `fermentation.scenario` | `Scenario`, `TemperaturePoint`, `Intervention`, `compile_scenario`, `CompiledScenario` |
| validation | `fermentation.validation` | `assert_conserved`, `assert_nonnegative`, `total_carbon`, `total_nitrogen`, `total_mass`, `BenchmarkSpec`, `ReferenceSeries`, `compare_series` |

## The core

### State vector
A single contiguous `float64` numpy array (`StateSchema` maps names → index
slices). Keeping it a plain array is what lets `scipy.integrate.solve_ivp` drive
it efficiently. Variables can be vectors: `S` (sugar) is one slot for wine and
three (glucose/maltose/maltotriose) for beer, so beer is an *addition* not a
rewrite.

Initial v1 contents: viable biomass `X`, sugar(s) `S`, ethanol `E`, yeast
assimilable nitrogen `N`, temperature `T`, evolved CO₂ (the experimentally
measurable proxy — a primary validation channel).

### Process
A `Process` contributes to `d(state)/dt`. It declares `name`, `tier`, and the
state variables it `touches`, and implements
`derivatives(t, y, schema, params) -> contribution`. `ProcessSet` sums the active
Processes (that sum is what `solve_ivp` integrates) and, crucially, **derives the
output tier** of each variable as the lowest tier among the Processes that touch
it. Toggling a speculative Process off leaves the validated core intact.

In `strict=True` mode, `ProcessSet` verifies every Process only writes to the
variables it declared — a cheap guard used in tests.

## Confidence tiers

`Tier` is an ordered enum (`VALIDATED > PLAUSIBLE > SPECULATIVE`). The trust of a
combination is the `min` (`Tier.combine`). Tiers are a property of *Processes and
parameters*, not of the raw floats flowing through the solver; an output's tier
is computed at the analysis boundary. This satisfies "the tier travels to every
output" without polluting the integration hot loop. (See DECISIONS #1.)

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

## Testing & quality gates

`uv run pytest` (unit + integration + conservation; benchmarks skipped),
`uv run ruff check .`, `uv run mypy` (strict on `src`). CI runs all three on
Python 3.13 and 3.14.
