# Architecture

This describes how the layers fit together. Milestones 0 and 1 are complete
(the tested skeleton, then single-strain isothermal primary fermentation passing
both §2.2 benchmarks); Milestone 2 (pH, aroma byproducts, SO₂, MLF, temperature
scheduling, discrete interventions, stochastic ensembles) is in progress. The
original brief is [`FERMENTATION_SIM_HANDOFF.md`](FERMENTATION_SIM_HANDOFF.md);
the per-decision record and where we deviated from it are in
[`DECISIONS.md`](DECISIONS.md) (D-1 … D-37), which is the canonical archive.

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
| core | `fermentation.core` | `Tier`, `StateSchema`, `VarSpec`, `StateVector`, `Process`, `ProcessSet`, `RateModifier`, `Medium`, `MEDIA`, `get_medium`, `wine_schema`, `beer_schema`; `chemistry` (molar masses, carbon/nitrogen fractions, Gay-Lussac split, `sugar_species`); `acidbase` (`solve_ph`, `ph_of_state`, `titratable_acidity`, `ph_tier`, `speciate_so2`, charge balance); `kinetics` (growth, uptake, ethanol inhibition, Arrhenius, esters/fusels, acetaldehyde, vicinal diketones, H₂S, malolactic, amino-acid ledger, autolysis, temperature ramp, carrying-capacity cap) |
| runtime | `fermentation.runtime` | `simulate`, `Trajectory`; `simulate_scheduled`, `ScheduledEvent`, `ScheduledTrajectory` (event loop); `simulate_ensemble`, `Ensemble` (stochastic wrapper) |
| scenario | `fermentation.scenario` | `Scenario`, `TemperaturePoint`, `Intervention`, `compile_scenario`, `CompiledScenario` (`.run` / `.run_ensemble`); intervention verbs `add_dap` / `add_so2` / `rack` / `pitch_mlf` |
| validation | `fermentation.validation` | `assert_conserved`, `assert_nonnegative`, `total_carbon`, `total_nitrogen`, `total_mass`, `BenchmarkSpec`, `ReferenceSeries`, `compare_series` |
| analysis | `fermentation.analysis` | `ph_series`, `titratable_acidity_series`, `molecular_so2_series`, `ibu_series` (top-layer observables over a `Trajectory`) |
| sensory | `fermentation.sensory` | `oav_series`, `sensory_profile`, `oav_tier`, `load_thresholds`, `AROMA_COMPOUNDS` — the speculative Tier-3 OAV aroma readout over a `Trajectory` (D-67) |

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

Milestone 2 grows the vector further with additive, isolable pools introduced by
their respective Processes — aroma byproducts (`esters`/`fusels` and their gas
sinks, `acetaldehyde`, diacetyl's α-acetolactate/`diacetyl`/butanediol pools,
`h2s`), the wine acid/SO₂ system (`tartaric`/`malic`/`lactic`/`cation_charge`,
`so2_total`), and the MLF/nitrogen machinery (`X_mlf`, `citrate`, `amino_acids`,
`debris`). Each defaults to 0 (or is dosed at the compile seam) so a scenario that
doesn't use it stays byte-for-byte the validated core (prime directive #3); most
are detailed in the pH/aroma sections below and in `DECISIONS.md`.

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

Two wrappers layer on top of this without changing the pure core:
- **Event loop** (`simulate_scheduled`, D-35) — segments a run at `ScheduledEvent`
  breakpoints (mutate / reconfigure / param_update) and restarts `simulate` per
  segment (a dose is a real discontinuity; BDF order-restart is correct — not
  `solve_ivp(events=)`, which can't mutate-and-resume). It carries an external-flow
  ledger so conservation across a jump is `final == initial + Σ flows`, and
  min-combines the per-segment tier map. `events=()` is byte-for-byte plain
  `simulate`. Temperature scheduling (a driven `TemperatureRamp`) and the discrete
  intervention verbs (D-36) both ride this one mechanism; `CompiledScenario.run()`
  always dispatches through it.
- **Stochastic ensemble** (`simulate_ensemble`, D-24/25/37) — Monte-Carlo over the
  parameters' `Uncertainty` bands (triangular default; LHS/Sobol via `qmc`), scoped
  to the active Process set's reads, returning nominal + median + P5/P95 band and
  per-member conservation. Randomness lives *only* here (seeded), keeping the core
  pure and reproducible. `simulate_ensemble(events=…)` runs the ensemble over a
  scheduled run (D-37).

## Scenarios as data

A `Scenario` (schema-validated YAML/JSON, **not** a custom DSL) declares initial
composition, organism/strain, temperature schedule, vessel, and a timeline of
interventions. No physics lives here, which keeps sweeps, Monte Carlo, and
cross-beverage reuse trivial.

## Media and the compile seam

A **`Medium`** (`fermentation.core.media`) names a beverage family and fixes its
`StateSchema` plus the Processes that act on it. `wine_schema()` has one sugar
slot; `beer_schema()` has three (`glucose`/`maltose`/`maltotriose`, in uptake
order). The `MEDIA` registry maps name → `Medium`; each medium fixes the Processes
that act on it (the M1 kinetics core plus the M2 Tier-2 mechanisms, speculative
ones staying isolable/togglable per prime directive #3).

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
  `BenchmarkSpec` data; the wine and beer `tests/benchmarks/` now **pass** (5
  benchmark tests), gated behind the `benchmark` pytest marker so they run via
  `uv run pytest -m benchmark` rather than in the default suite. `ReferenceSeries`
  + `compare_series` (RMSE/MAE) are the seam for scoring against *real* measured
  datasets when we obtain them.

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
- **SO₂ speciation = the first pH consumer, readout-only (D-22, D-28, D-51).** `wine_schema`
  appends a slot `so2_total` (total SO₂ as g/L SO₂-equivalent, dosed via `so2_total_mgl`,
  conserved/inert). `acidbase.speciate_so2` solves pH from the organic acids, then splits the
  total into **bound** vs **free** via a competitive-Langmuir carbonyl-bisulfite equilibrium:
  `bound_so2_molar` takes a tuple of `(molar_concentration, Kd)` per carbonyl and solves one
  shared "reactive bisulfite" root `h` via `brentq`, from which each carbonyl's bound share is
  `Aᵢ·h/(Kᵢ+h)` — reducing exactly to the original D-28 single-carbonyl closed form
  `(A−x)(C−x)β − Kx = 0` when only one carbonyl is active. **D-51** (2026-07-07) generalised D-28
  from acetaldehyde alone to **acetaldehyde + pyruvate + α-ketoglutarate together**, all worked in
  **moles** (`M_ACETALDEHYDE`/`M_PYRUVATE`/`M_ALPHA_KETOGLUTARATE`), since bisulfite competition
  is molar and the three carbonyls have very different molar masses; `free_acetaldehyde` reads
  back only acetaldehyde's own bound share, so competing keto-acid pools measurably reduce
  acetaldehyde's SO₂ protection. Returns the **molecular** (antimicrobial) fraction
  `1/(1+10^(pH−pKa₁))` (sulfurous pKa₁ 1.81) of *free* —
  the D-18 coupling *emerging*, not scripted. It is **readout-only**: SO₂ is kept out of the
  charge balance (the inverse anchoring makes in-balance vs readout identical at t=0) and out
  of titratable acidity (OIV excludes it), and is carbon-free — so dosing it leaves pH and
  `total_carbon` byte-for-byte (an isolability test pins this). At all carbonyls=0 the split
  collapses to D-22 exactly (`free == total`). The lone RHS consumer is the MLF antimicrobial
  gate, which reads the *derived* free-molecular SO₂ (bound SO₂ is not antimicrobial), so the
  early acetaldehyde peak *and* the always-on keto-acid pools transiently/persistently sequester
  SO₂ and relax suppression — an emergent competition. The bound-acetaldehyde-protected-from-ADH
  feedback stays deferred (readout-only).
- **Acetaldehyde** (`core/kinetics/acetaldehyde.py`, decision D-27) — the obligate main-
  pathway intermediate, modelled as a transient **ethanol-carbon buffer**: flux-linked
  production *borrows* a C2 slice of ethanol and viable-yeast-gated reduction *returns* it
  (both mole-for-mole C2→C2). It de-lumps the uptake Process's single sugar→ethanol step
  rather than adding a parallel pathway, so carbon closes touching neither `S` nor `CO2` and
  the `E` endpoint (hence the §2.2 benchmarks) is preserved to relative ~1e-8. The early
  produce-then-reabsorb peak emerges; a crash strands it (the D-26 live-yeast-gating shape).

## The sensory / OAV readout — the speculative Tier-3 aroma lens (D-67)

`fermentation.sensory` is a **top-layer readout** (sibling of `analysis`) that maps
Odor-Activity-Values over a finished `Trajectory`: `OAV_i = concentration_i / threshold_i`
for each aroma-active pool the chemistry already tracks (`esters`, `fusels`, `diacetyl`,
`acetaldehyde`, `h2s`; wine adds `ethylphenols`/`ethylguaiacols`/`mercaptans`). It adds **no
state, no Process, no ledger entry** — the full suite stays byte-for-byte green (isolation by
construction). Opened as the first beat of Milestone 3 (`docs/plans/milestone-3-plan.md`).

- **The §4.2 cardinal rule / firewall.** The sensory layer consumes the chemistry; the
  chemistry never imports it back. Thresholds load **standalone** (`load_thresholds()` reads
  `parameters/data/sensory.yaml`) and are **never** merged into any `CompiledScenario` at the
  compile seam — no RHS reads a perception threshold, so the chemistry never even sees these
  numbers (a stronger isolation than any Tier-2 readout, which *is* merged because a Process
  reads it). A deliberate consequence (D-24): thresholds sit **outside** the ensemble sweep.
- **The tier floor (§4.3 credibility firewall).** `oav_tier(input, threshold)` returns
  `combine(input, threshold, SPECULATIVE)` → **always speculative, even for a validated
  input**. The explicit `SPECULATIVE` is not redundant with the threshold's tier: the sensory
  *mapping itself* is the canonical speculative case (`Tier` docstring). A pure-function test
  (`oav_tier(VALIDATED, VALIDATED) is SPECULATIVE`) proves the floor non-vacuously — a
  real-trajectory test would be a tautology since every aroma pool is speculative/plausible.
- **Matrix-specific thresholds, µg/L.** Keys are `threshold_<pool>_<beer|wine>` because
  ethanol/matrix shift odor thresholds; each `conditions` records the **measurement matrix**
  (not the same as the application medium — a water/model-solution measurement is flagged as a
  matrix gap in `notes`). Stored in µg/L (the literature unit), crossed to canonical g/L at
  the boundary via `units.convert.ugl_to_gpl`. `sensory_profile` reports **per-compound** OAVs
  + above-threshold flags (never a summed scalar — summing assumes contested additivity).
- **Lumped pools (D-66).** `esters`/`fusels`/`mercaptans` are read against one named
  representative's threshold (isoamyl acetate / isoamyl alcohol / methanethiol), with the
  "assumes fixed lump composition" honesty cost flagged in provenance. `iso_alpha`/IBU is
  excluded — a *taste*, already read out by `ibu_series` (D-64). Descriptor-space projection
  ("smells like leather and banana") is the deferred, even-more-speculative beat 1b.

## Testing & quality gates

`uv run pytest` (unit + integration + conservation; benchmarks skipped),
`uv run ruff check .`, `uv run mypy` (strict on `src`). CI runs all three on
Python 3.13 and 3.14.
