# Architecture

**What this document is:** a map of how the code is put together *right now* — layers, packages,
the state vector, the Process registry, and where each subsystem lives. It answers "where does X
live and what may it import."

**What it is not:** a history. *Why* anything is shaped this way, what was tried and rejected, and
what is still open all live in [`DECISIONS.md`](DECISIONS.md), which is the canonical archive.
Individual D-numbers are cited here as pointers into it. Deliberately **no `D-1 … D-n` range and no
milestone status line** — those pointers went stale twice; the archive's generated index carries
the live count, and [`plans/milestone-3-plan.md`](plans/milestone-3-plan.md) is a frozen log, not a
status. The original brief is [`FERMENTATION_SIM_HANDOFF.md`](FERMENTATION_SIM_HANDOFF.md)
(reference, not gospel).

Counts in this file are measured from the code, not remembered. Every one of them is reproducible
by the snippet in [Checking this document](#checking-this-document) at the end.

## Layering

Four layers with strictly one-directional dependencies — a lower layer never imports a higher one:

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

`analysis` and `sensory` are **top-layer readouts**, siblings of `validation`: they consume a
finished `Trajectory` and are imported by nothing lower. The chemistry never imports them back.

### Package map

| Layer | Package | Contents |
|-------|---------|----------|
| parameters | `fermentation.parameters` | `Parameter`, `Provenance`, `Uncertainty`, `ParameterSet`, `load_parameters`, `default_data_dir` — plus the 19 YAML data files |
| units | `fermentation.units` | `brix_to_sg`, `sg_to_plato`, `abv_from_ethanol`, `ugl_to_gpl`, … |
| core | `fermentation.core` | `Tier`; `StateSchema`, `VarSpec`, `StateVector`; `Process`, `ProcessSet`, `RateModifier`; `Medium`, `MEDIA`, `get_medium`, `wine_schema`, `beer_schema`; `chemistry`; `acidbase`; the `kinetics` subpackage |
| runtime | `fermentation.runtime` | `simulate`, `Trajectory`; `simulate_scheduled`, `ScheduledEvent`, `ScheduledTrajectory`; `simulate_ensemble`, `Ensemble` |
| scenario | `fermentation.scenario` | `Scenario`, `TemperaturePoint`, `Intervention`, `compile_scenario`, `CompiledScenario`; intervention verbs |
| validation | `fermentation.validation` | `assert_conserved`, `assert_nonnegative`, `total_carbon`, `total_nitrogen`, `total_mass`, `BenchmarkSpec`, `ReferenceSeries`, `compare_series` |
| analysis | `fermentation.analysis` | `ph_series`, `titratable_acidity_series`, `molecular_so2_series`, `free_so2_series`, `bound_so2_series`, `ibu_series`, `astringency_series`, `polymeric_pigment_series`, `color_series`, `observed_color_series`, `attribute_spread` |
| sensory | `fermentation.sensory` | `oav_series`, `sensory_profile`, `oav_tier`, `load_thresholds`, `AROMA_COMPOUNDS`; `MaxRuleProjector`/`DescriptorProjector`; `StevensProjector`, `dominant_flip_sensitivity` |

## The core

### State vector

A single contiguous `float64` numpy array; `StateSchema` maps names → index slices. Keeping it a
plain array is what lets `scipy.integrate.solve_ivp` drive it efficiently. Variables can be
vectors: `S` (sugar) is one slot for wine and three for beer
(glucose/maltose/maltotriose, in uptake order), so beer is an *addition*, not a rewrite.

Current size, from `get_medium(...).schema`:

| Medium | Named variables | Float slots | Sugar slots |
|--------|-----------------|-------------|-------------|
| wine   | 94 | 94 | 1 |
| beer   | 55 | 57 | 3 |

Tier and uncertainty do **not** ride inside these floats — they are properties of Processes and
parameters, derived at the analysis boundary (D-1).

The vector is grown only by *additive, isolable pools* introduced by their own Processes. Each
new pool declares a `VarSpec.default` (almost always 0, or is dosed at the compile seam) so a
scenario that does not use it is byte-for-byte the validated core — prime directive #3. Broadly
the wine-only slots cover the acid/SO₂ system, the speciated amino-acid and keto-acid pools, MLF
and Brett biomass, the oxidative-aging products (browning index, Strecker aldehydes, quinone), the
grape colour axis (anthocyanin/tannin/polymeric pigment), oak, closure and bound sulfides; the
beer-only slots are its own organic acids (`acetic`, `formic`, `oxalic`, `pyruvic`, `succinic`),
`peptide_buffer` and `iso_alpha`. `quinone` is present in **both** media regardless of which
oxidative set is wired, so slot indices do not move when the set changes.

### Process, ProcessSet, RateModifier

A `Process` contributes to `d(state)/dt`. It declares `name`, `tier`, and the state variables it
`touches`, and implements `derivatives(t, y, schema, params) -> contribution`. `ProcessSet` sums
the active Processes — that sum is what `solve_ivp` integrates — and **derives the output tier** of
each variable as the lowest tier among the Processes touching it. In `strict=True` mode it also
verifies every Process writes only to the variables it declared (used throughout the tests).

Some mechanisms *scale* an existing flux rather than *add* a new one. Because `ProcessSet` sums,
those cannot be Processes; a `RateModifier` names the Processes it `modifies` and returns a scalar
`factor`, which is multiplied onto the target's entire contribution vector before summing. Scaling
a conserving Process's whole contribution by one scalar preserves its atom balances, so a modifier
can never break conservation, and the `touches` contract still holds (scaling zeros stays zero).
Modifiers toggle and feed `tier_of` exactly like Processes. Stacked modifiers on one Process
compose to a single scalar (D-10, D-11).

`core/kinetics/` holds **80** concrete `Process`/`RateModifier` implementations across 23 modules
(the three base types live in `core/process.py`).

### Kinetics modules

| Module | Impls | What it covers |
|--------|------:|----------------|
| `growth.py` | 1 | biomass growth |
| `uptake.py` | 1 | fermentative sugar uptake |
| `inhibition.py` | 1 | ethanol inhibition (modifier) |
| `osmotic.py` | 1 | high-sugar osmotic/substrate brake (modifier, wine-only) |
| `arrhenius.py` | 2 | per-rate temperature scaling (modifier) |
| `temperature.py` | 1 | driven temperature ramp |
| `carrying_capacity.py` | 1 | biomass cap (modifier) |
| `inactivation.py` | 2 | death, ethanol-tolerance death |
| `autolysis.py` | 1 | autolysis → debris |
| `byproducts.py` | 4 | glycerol and the realised-yield byproduct sinks |
| `acetaldehyde.py` | 2 | the ethanol-carbon buffer intermediate (D-27) |
| `vicinal_diketones.py` | 3 | α-acetolactate → diacetyl → butanediol |
| `hydrogen_sulfide.py` | 3 | H₂S |
| `mercaptans.py` | 1 | methanethiol — the last lumped aroma pool |
| `amino_acids.py`, `amino_acid_pools.py` | 1 | the eight speciated amino-acid pools + ledger (D-100) |
| `keto_acids.py` | 6 | pyruvate / α-ketoglutarate / α-ketobutyrate, excretion + reassimilation |
| `carbon_routing.py` | — | shared carbon draw/refund helpers, ester and fusel route specs, label tracer |
| `precursor_fates.py` | 1 | precursor partitioning |
| `organic_acids.py` | 3 | beer's organic acids: excretion, acetic overflow, wort acid removal |
| `malolactic.py` | 6 | MLF conversion, citrate, diacetyl reduction, growth/death/senescence |
| `brett.py` | 6 | *Brettanomyces* growth/death/toxicity, decarboxylation, vinylphenol reduction, yeast POF |
| `hops.py` | 1 | iso-α-acid loss (IBU) |
| `aging.py` | 24 | the whole aging axis (below) |
| `oxidative_cascade.py` | 8 | the Fe(II)+O₂ activation alternative (below) |
| `o2_partition.py` | — | O₂ partitioning helper |

`aging.py` is the largest single module. Its 24 Processes cover ester hydrolysis and
esterification; the oxidative sinks (`OxidativeAcetaldehyde`, `SulfiteOxidation`,
`PhenolicBrowning`, `AntioxidantBurstOxidation`, `StreckerDegradation`); the non-oxidative thermal
routes (`MaillardStrecker`, `SotolonAldolCondensation`, `Caramelization`, `MaillardBrowning`); oak
(`OakExtraction`, `EllagitanninOxidation`); the colour axis
(`TanninAnthocyaninCondensation`, `AcetaldehydeBridgedCondensation`, `AnthocyaninFading`,
`ThermalAnthocyaninFade`, `TanninSelfPolymerization`, `TanninEthylTanninCondensation`); and
`SMMHydrolysis` (DMS), the two bound-sulfide releases, and `ClosureOxygenIngress`.

## Media, the process registry, and the compile seam

A **`Medium`** (`core/media.py`) is a plain record of four things: `name`, `schema`,
`process_factories`, `modifier_factories`. It is assembled from ~35 named tuples
(`_PRIMARY_FERMENTATION_PROCESSES`, `_AGING_PROCESSES`, `_MLF_PROCESSES`, `_BRETT_PROCESSES`, …),
each an independently toggleable group. Every one of those groups except the oxidative sets is
**additive**: off at the compile seam, switched on by an intervention, and an un-used group leaves
the run byte-for-byte unchanged.

`get_medium(name, *, oxidative="direct")` returns the wired medium. What each build contains:

| Medium | Oxidative set | Processes | Modifiers |
|--------|---------------|----------:|----------:|
| wine | `direct` (default) | 62 | 5 |
| wine | `cascade` | 64 | 5 |
| wine | `direct_burst` | 63 | 5 |
| beer | `direct` (default) | 27 | 3 |
| beer | `cascade` | 28 | 3 |
| beer | `direct_burst` | 27 | 3 |

### The three oxidative sets (D-141, extended D-147)

The oxidative axis is the one place where toggling is a *swap*, not an addition — both alternatives
draw on the same `o2` pool, so a build carrying both would silently double-count it.

- **`direct`** — the default. Six calibrated sinks each draw their own share straight from the
  shared `o2` pool, split by `k_i / Σk` through `ProcessSet`'s summing. It reads as competition and
  sums correctly, but asserts that ethanol, bisulfite, phenolics, amino acids, anthocyanin and oak
  tannin each react with dissolved O₂ — which they do not.
- **`cascade`** — routes all six behind a single Fe(II)+O₂ activation node that consumes the O₂ and
  produces two oxidants per mole; each former sink is re-homed onto whichever oxidant actually
  oxidises it. **Mutually exclusive with `direct`.**
- **`direct_burst`** — `direct` plus one further sink (`AntioxidantBurstOxidation`). A *superset of
  direct*, not a third mechanism, and opt-in rather than default. Wine-only in effect: beer's burst
  build is identical to its direct build, because the slot the sink needs is wine-only. There is
  deliberately no `cascade_burst`.

### The compile seam

`compile_scenario(scenario)` (`scenario/compile.py`) is the scenario→core seam and the **only**
place industry units cross into canonical ones (°Brix → g/L, °C → K, days → hours). It validates
the `scenario.initial` vocabulary per medium, seeds initial temperature from the schedule, loads
`<medium>_<strain>.yaml` over the shared parameter files, assembles the medium's `ProcessSet`, and
returns a `CompiledScenario` (`y0`, `process_set`, `parameters` + resolved `param_values`,
`schema`, `t_span_h`) that drops straight into `simulate`. Beer's three sugars are supplied
explicitly rather than split from a single OG — that wort spectrum is a provenance-backed
parameter, not a constant in the seam (D-7).

A `Scenario` is schema-validated YAML/JSON, **not** a custom DSL, and holds no physics.

**Intervention verbs dispatch through two tables (D-187).** Almost every verb is a function of its
own `Intervention` alone (`_INTERVENTION_VERBS`). `_SCENARIO_INTERVENTION_VERBS` holds the ones
whose *magnitude* comes from a scenario-level field, and `seal_bottle` is the first: it doses
`bottling_burst_<closure>` for whatever `scenario.closure` names, so its amount is sourced instead
of author-supplied. Both tables are searched before an action is called unknown. Cross-cutting
gates live in `_compile_interventions` rather than in the verbs, because a verb never sees the
scenario: `set_ph` needs `initial_ph` (D-186), and `seal_bottle` needs a `closure` and must not
precede `begin_aging` (the charge is net of that first month's steady ingress).

## Confidence tiers

`Tier` is an ordered enum (`VALIDATED > PLAUSIBLE > SPECULATIVE`); the trust of a combination is
the `min` (`Tier.combine`). Tiers belong to *Processes and parameters*, not to the floats flowing
through the solver, so the integration hot loop stays clean and an output's tier is computed at the
analysis boundary (D-1).

No Process or parameter is `VALIDATED`: that tier is reserved for checks against independent
*measured* time-series, which the project does not have. Passing the §2.2 benchmarks earns
`PLAUSIBLE` — sound forms, sourced parameters, reproduces the keystone model.

## Parameters with provenance

Every kinetic or physical constant is a `Parameter` requiring value, units, tier, an `Uncertainty`
range, and `Provenance` (source + measurement conditions). The Pydantic models reject any entry
missing these, so "no magic numbers" is a load-time guarantee, not a convention. Strain-specific
overlays merge on top of generic defaults (`ParameterSet.merge`).

The 19 files in `parameters/data/`: `wine_generic.yaml` and `beer_generic.yaml` (the per-medium
bases); `acidbase.yaml`, `beer_acids.yaml`, `closure.yaml`; `acetaldehyde.yaml`,
`keto_acids.yaml`, `vicinal_diketones.yaml`, `hydrogen_sulfide.yaml`, `bound_sulfides.yaml`,
`dms.yaml`; `aging.yaml`, `oak.yaml`, `polymerization.yaml`, `thermal.yaml`; `hops.yaml`,
`additions.yaml`; and the two that load **standalone**, outside any `CompiledScenario` —
`sensory.yaml` and `psychophysics.yaml`.

## Units boundary

Canonical internal units: concentration **g/L** (≡ SI kg/m³), temperature **K**, time **hours**.
Industry units (°Brix, SG, °Plato, %ABV, °C, days) appear only on the far side of
`fermentation.units` (D-3).

## Runtime

`simulate(process_set, params, y0, t_span)` wraps `solve_ivp` with an implicit adaptive method
(BDF by default — fermentation is stiff) and returns a `Trajectory` carrying the time grid, state
history, and derived tier map.

Two wrappers layer on top without changing the pure core:

- **Event loop** (`simulate_scheduled`, D-35) — segments a run at `ScheduledEvent` breakpoints
  (mutate / reconfigure / param_update) and restarts `simulate` per segment. A dose is a real
  discontinuity, so a BDF order-restart is correct — not `solve_ivp(events=)`, which cannot
  mutate-and-resume. It carries an external-flow ledger so conservation across a jump is
  `final == initial + Σ flows`, and min-combines the per-segment tier map. `events=()` is
  byte-for-byte plain `simulate`. Temperature scheduling and every discrete intervention verb ride
  this one mechanism; `CompiledScenario.run()` always dispatches through it.
- **Stochastic ensemble** (`simulate_ensemble`, D-24/25/37) — Monte-Carlo over the parameters'
  `Uncertainty` bands (triangular default; LHS/Sobol via `qmc`), scoped to the active Process set's
  reads, returning nominal + median + P5/P95 band and per-member conservation. Randomness lives
  **only** here and is seeded, keeping the core pure and reproducible.

Note that `only=`/`exclude=` shift the draw sequence, so an arm and its baseline are two different
random ensembles unless pinned to a fixed hypercube.

## pH as a derived pure function

pH is **not** integrated — there is no `dpH/dt`. Like `total_carbon` and ABV it is an
instantaneous, pure algebraic function of state: `core/acidbase.py` solves electroneutrality
`Σ charge = 0` for `[H⁺]` (a 1-D monotonic root-find in pH-space via `brentq`) given the
charge-active acids and a pKa set, and reports `pH = −log₁₀[H⁺]`. Building it as a full proton
balance rather than a tracked-pH approximation is what makes the couplings — MLF deacidification,
SO₂ speciation — *emerge* rather than be scripted (D-18).

**The acid registry is per-medium, not medium-agnostic (D-179).** `acidbase.acid_registry` selects
off `StateSchema.medium`:

- **`ACID_STATE` — wine's registry.** `tartaric`, `malic`, `lactic` plus `cation_charge`, the net
  strong-cation charge density (mol⁺/L, K⁺-dominant). The cation is **mandatory** (weak acids alone
  give pH ≈ 2.3 against a real ~3.3) and is **back-solved from a measured `initial_ph`** at the
  compile seam, so the model predicts pH *changes*, not absolute initial pH.
- **`BEER_ACIDS` — beer's registry**, beside it: `acetic`, `formic`, `oxalic`, `pyruvic`,
  `succinic`, plus `peptide_buffer`. **Beer's pH is a prediction**, not an inverse-anchored fit
  like wine's — which is the point, and also why it is the harder claim.

The `Byp` pool is read as a succinic-equivalent acid (zero new carbon, so `total_carbon` is
unchanged). Scalar `ph_of_state` / `titratable_acidity` are pure and live in core; the
trajectory-series helpers need a `Trajectory` and therefore sit one layer up in
`fermentation.analysis`.

**Anchoring runs in two places, and they are different functions.** `solve_cation_charge` anchors
at the *compile seam* from scenario inputs, with `Byp` and dissolved CO₂ structurally 0 because
nothing has fermented. `cation_charge_for_ph` anchors a *state*, reading those two off the vector —
the exact inverse of `ph_of_state`, term for term. The second exists because `initial_ph` fixes t=0
only and the ferment then drags pH somewhere the scenario never chose, which made an aging study at
a stated pH unwritable (D-150 measured the gap; D-186 closes it). The `set_ph` verb is its one
caller: no Process touches `cation_charge`, so it is constant across a plain run and moves only at a
scheduled adjustment. Both directions are real cellar operations on the same quantity — carbonate
deacidification raises the cation, cation-exchange resin lowers it — which is why the verb is
cation-moving rather than a pH dial; acidifying by *addition* is `add_acid` instead.

**SO₂ speciation** is the first pH consumer and is readout-only. `speciate_so2` solves pH from the
organic acids, then splits total SO₂ into bound vs free via a competitive-Langmuir
carbonyl-bisulfite equilibrium: `bound_so2_molar` takes `(molar_concentration, Kd)` per carbonyl
and solves one shared reactive-bisulfite root, so each carbonyl's bound share is `Aᵢ·h/(Kᵢ+h)`.
Competition is molar and the carbonyls differ greatly in molar mass, so this is worked in moles
(acetaldehyde + pyruvate + α-ketoglutarate together, D-51). It returns the **molecular**
(antimicrobial) fraction of *free* SO₂ — the coupling emerging, not scripted. SO₂ is kept out of
the charge balance and out of titratable acidity, and is carbon-free, so dosing it leaves pH and
`total_carbon` byte-for-byte. Its one RHS consumer is the MLF antimicrobial gate, which reads the
*derived* free-molecular value, so the early acetaldehyde peak and the always-on keto-acid pools
transiently sequester SO₂ and relax suppression.

## Readout layers

### `analysis` — chemistry observables

Trajectory-level series that need no perception model: pH, titratable acidity, the three SO₂
series, IBU, and the colour/mouthfeel axis (`astringency_series`, `polymeric_pigment_series`,
`color_series`, `observed_color_series`), plus `attribute_spread`.

### `sensory` — the speculative Tier-3 aroma lens

Maps Odor-Activity-Values over a finished `Trajectory` (`OAV = concentration / threshold`). It adds
**no state, no Process, no ledger entry**.

- **The firewall.** The sensory layer consumes the chemistry; the chemistry never imports it back.
  Thresholds load standalone (`load_thresholds()` reads `sensory.yaml`) and are **never** merged
  into a `CompiledScenario`, so no RHS ever sees a perception threshold — stronger isolation than
  any Tier-2 readout. A deliberate consequence: thresholds sit **outside** the ensemble sweep.
- **The tier floor.** `oav_tier` returns `combine(input, threshold, SPECULATIVE)` — always
  speculative, even for a validated input, because the sensory *mapping itself* is the canonical
  speculative case.
- **Matrix-specific thresholds, µg/L.** Keys are `threshold_<pool>_<beer|wine>` because
  ethanol/matrix shift odor thresholds; `conditions` records the measurement matrix, and a
  water/model-solution measurement is flagged as a matrix gap. `sensory_profile` reports
  **per-compound** OAVs and above-threshold flags, never a summed scalar.
- **Descriptor projection** (`sensory/descriptors.py`, D-95) projects the OAV vector onto 14 (wine)
  / 9 (beer) descriptor axes behind the `DescriptorProjector` Protocol, so a panel-trained model
  could swap in. It uses a **max rule, not a sum** — the layer beneath refuses to sum OAVs, so a
  summing projector would silently reintroduce contested additivity. Each descriptor reports its
  loudest contributor. Membership is structure (binary), so it lives in code and mints no
  constants; the axis set is derived per medium, so beer can never report a wine-only descriptor.
- **Stevens compression** (`sensory/compression.py`, D-98) compresses each contributor's OAV to a
  perceived intensity (`I = OAV ** n`) *before* the max rule. **Isolable and not the default** —
  delete `psychophysics.yaml` and the layer beneath is byte-for-byte unaffected. It can neither
  invent nor silence a detectable smell (`I > 1` iff `OAV > 1`), and its exponents are author
  estimates, so read it via `dominant_flip_sensitivity` rather than as a bare dominant.

`mercaptans` (methanethiol) is the last lumped pool in the project; the lump caveat derives from
the `AromaCompound.lumped` flag, not a hardcoded list, so it cannot linger on a pool that stopped
being lumped.

## Validation

Two disciplines, both as code:

- **Conservation invariants** — `assert_conserved` / `assert_nonnegative` take a model-supplied
  conserved-quantity function and check it along a trajectory. `total_carbon`, `total_nitrogen` and
  `total_mass` weight each state variable using the shared stoichiometry in `core/chemistry.py`, so
  a check can never disagree with the kinetics it audits. Carbon and nitrogen are rigorous atom
  balances; mass is scoped to the abiotic `S + E + CO₂` conversion (D-8).
- **Benchmark curves** — the §2.2 acceptance criteria are encoded as `BenchmarkSpec` data and
  **pass**, gated behind the `benchmark` pytest marker (`uv run pytest -m benchmark`).
  `ReferenceSeries` + `compare_series` are the seam for scoring against real measured datasets.

## Testing & quality gates

`uv run pytest -n auto` (76 test files; unit, integration, conservation, sampling-surface and
doc-consistency checks), `uv run ruff check .`, `uv run mypy` (strict on `src`). CI runs all three
on Python 3.13 and 3.14. Two of the test files guard documentation rather than physics:
`test_decisions_index.py` (the archive's generated index) and `test_memory_shape_hook.py`.

## Checking this document

Every count above is derived, not remembered. To re-derive them, from the repo root:

```
uv run python -c "
import pathlib, re
from fermentation.core.media import MEDIA, get_medium
for name in sorted(MEDIA):
    for ox in ('direct', 'cascade', 'direct_burst'):
        m = get_medium(name, oxidative=ox)
        print(name, ox, 'slots', m.schema.size, 'vars', len(m.schema.names),
              'procs', len(m.process_factories), 'mods', len(m.modifier_factories))
root = pathlib.Path('src/fermentation')
pat = re.compile(r'^class .*(?:Process|RateModifier)', re.M)
print('parameter files', len(list((root / 'parameters/data').glob('*.yaml'))))
print('test files', len(list(pathlib.Path('tests').rglob('test_*.py'))))
print('kinetics impls', sum(len(pat.findall(p.read_text(encoding='utf-8')))
                            for p in (root / 'core/kinetics').glob('*.py')))
"
```

It is one Python call on purpose. The obvious version of this check is a handful of shell
one-liners (`ls | wc -l`, `find`, `grep | awk`), and those run in Git Bash but **not** in
PowerShell, which is this project's primary shell — so for many readers the check would simply
error, while still *looking* like verification. A self-check that cannot execute is worse than none.

If a number here disagrees with that output, **the code is right and this document is stale** —
fix the document. Per `CLAUDE.md`, a beat that adds a Process, a state slot, a parameter file or a
package updates this file in the same commit.
