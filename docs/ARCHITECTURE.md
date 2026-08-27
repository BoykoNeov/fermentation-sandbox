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
| wine   | 97 | 97 | 1 |
| beer   | 57 | 59 | 3 |

Tier and uncertainty do **not** ride inside these floats — they are properties of Processes and
parameters, derived at the analysis boundary (D-1).

The vector is grown only by *additive, isolable pools* introduced by their own Processes. Each
new pool declares a `VarSpec.default` (almost always 0, or is dosed at the compile seam) so a
scenario that does not use it is byte-for-byte the validated core — prime directive #3. Broadly
the wine-only slots cover the acid/SO₂ system, the speciated amino-acid and keto-acid pools, MLF
and Brett biomass, the oxidative-aging products (browning index, Strecker aldehydes, quinone), the
grape colour axis (anthocyanin/tannin/polymeric pigment), oak, closure, bound sulfides and the
dosed antioxidant `ascorbate` (D-202); the
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

`core/kinetics/` holds **83** concrete `Process`/`RateModifier` implementations across 24 modules
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
| `byproducts.py` | 4 | glycerol and the realised-yield byproduct sinks; the two aroma producers each have a growth-coupled beer subclass (below) |
| `acetaldehyde.py` | 2 | the ethanol-carbon buffer intermediate (D-27) |
| `vicinal_diketones.py` | 3 | α-acetolactate → diacetyl → butanediol |
| `hydrogen_sulfide.py` | 3 | H₂S |
| `mercaptans.py` | 1 | methanethiol — the last lumped aroma pool |
| `amino_acids.py`, `amino_acid_pools.py` | 1 | the eight speciated amino-acid pools + ledger (D-100) |
| `keto_acids.py` | 6 | pyruvate / α-ketoglutarate / α-ketobutyrate, excretion + reassimilation |
| `carbon_routing.py` | — | shared carbon draw/refund helpers, ester and fusel route specs, label tracer, and the shared fermentative rate law (per-slot uptake rates D-180, evolved-CO₂ rate D-227) that every Process reading the flux must call rather than re-derive |
| `precursor_fates.py` | 1 | precursor partitioning |
| `organic_acids.py` | 3 | beer's organic acids: excretion, acetic overflow, wort acid removal |
| `wort_oxygen.py` | 1 | beer's wort aeration O₂, stripped by the yeast in the lag phase (D-213) |
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
| wine | `cascade` | 66 | 5 |
| wine | `direct_burst` | 63 | 5 |
| beer | `direct` (default) | 28 | 3 |
| beer | `cascade` | 29 | 3 |
| beer | `direct_burst` | 28 | 3 |

### The one place the media run different rate LAWS (D-226)

Everywhere else, a per-medium difference is a *parameter value* or a *group membership*; no
Process branches on a medium string. The aroma producers are the exception, and it is wired as
membership so that it stays visible in this table rather than hiding inside a `derivatives`:

| Group | wine | beer |
|-------|------|------|
| `_BYPRODUCT_PROCESSES` | `EsterVolatilization` | `EsterVolatilization` |
| `_WINE_AROMA_PRODUCERS` | `EsterSynthesis`, `FuselAlcoholsEhrlich` | — |
| `_BEER_AROMA_PRODUCERS` | — | `EsterSynthesisGrowthCoupled`, `FuselAlcoholsEhrlichGrowthCoupled` |

Wine's two producers ride the fermentative flux (biomass-**hours**); beer's ride `mu·X` (biomass
**formed**). The beer pair reads the *base* growth rate, so both are named as extra targets of
`ArrheniusTemperature.for_growth` in `_BEER_FERMENTATION_MODIFIERS` — the D-32 idiom, and the
condition for their integral to be the conserved ΔX. The totals in the table above are unchanged:
each medium still wires three byproduct Processes.

`DECISIONS.md` D-226 is the answer to *why* the two media differ; this table is only the answer to
*where*.

### The three oxidative sets (D-141, extended D-147)

The oxidative axis is the one place where toggling is a *swap*, not an addition — both alternatives
draw on the same `o2` pool, so a build carrying both would silently double-count it.

- **`direct`** — the default. Six calibrated sinks each draw their own share straight from the
  shared `o2` pool, split by `k_i / Σk` through `ProcessSet`'s summing. It reads as competition and
  sums correctly, but asserts that ethanol, bisulfite, phenolics, amino acids, anthocyanin and oak
  tannin each react with dissolved O₂ — which they do not.
- **`cascade`** — routes all six behind a single Fe(II)+O₂ activation node that consumes the O₂ and
  produces two oxidants per mole; each former sink is re-homed onto whichever oxidant actually
  oxidises it. **Mutually exclusive with `direct`.** It also carries one consumer that is *not* a
  re-home: `QuinoneHydrogenSulfideCapture` (D-201, wine-only), the quinone→H₂S sulfide sink. It is
  the case where a share of the quinone node is the wrong measure of worth — it takes 0.003 % of
  that node while removing ~10 % of the dissolved sulfide pool over one oxidation challenge, and
  it is the model's only passive post-fermentation sink on `h2s`. A second non-re-home joined it at
  D-202: `QuinoneAscorbateReduction` (wine-only), which completes the four-member top group of
  quinone nucleophiles. It is the only Process in any build that is inert by default *state* rather
  than by wiring — its `ascorbate` slot defaults to 0 because the source says new wine has
  negligible ascorbic acid, so it costs nothing until a scenario calls `add_ascorbate`.
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

**A dosing verb doses the whole compound, not the ion it was named for (D-210).** `add_dap` wrote
only `N` for six milestones, on the premise that the model tracks no phosphorus pool. The premise
was true and the inference was not: the charge balance needs a *total*, not a pool, and the dropped
phosphate was ~0.95 equivalents per mole of anion charge — so the verb was booking a salt as a
base. It now writes three slots (`N`, `phosphate`, `nitrogen_charge_excess`), all guarded on slot
presence so a schema predating D-210 behaves as it did. The dose-time pH *rise* is emergent from
the balance and must stay that way: +2 cation and ~−0.95 anion per mole leaves ~+1.05 net, which
*is* the protons the dosed HPO₄²⁻ takes up becoming H₂PO₄⁻ at wine pH.

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

A **counted pitch** (cells/mL, which is how the literature states one) crosses here too,
via `cells_per_ml_to_pitch_gpl`. That one is not a convenience factor: it carries the
per-cell dry mass defining the gram this engine's biomass is expressed in, so its
docstring — not this page — is where the argument lives (D-219).

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
  this one mechanism; `CompiledScenario.run()` always dispatches through it. A `mutate` is handed
  the **running** parameter map at its breakpoint (`StateMutation` is `(schema, state, params) ->
  state`, D-235), so a jump that back-solves state from parameters uses the member's own draw under
  an ensemble; 2 of the 13 verbs read it (`set_ph`, `add_dap`) and the rest ignore it.
- **Stochastic ensemble** (`simulate_ensemble`, D-24/25/37) — Monte-Carlo over the parameters'
  `Uncertainty` bands (triangular default; LHS/Sobol via `qmc`), scoped to the active Process set's
  reads, returning nominal + median + P5/P95 band and per-member conservation. Randomness lives
  **only** here and is seeded, keeping the core pure and reproducible.

Note that `only=`/`exclude=` shift the draw sequence, so an arm and its baseline are two different
random ensembles unless pinned to a fixed hypercube.

`y0` is shared across members with one exception: the parts of it a **parameter derives**
(D-233, D-236, D-238, D-241). `simulate_ensemble` takes an optional `y0_for_member` builder —
omitted ⇒ the fixed array, so a direct caller is byte-identical — and `CompiledScenario.run_ensemble`
supplies one via `CompiledScenario.y0_for_member()`. That method composes the rules of
`CompiledScenario._member_seed_rules()` and returns `None` when none applies. Four rule families
ship, and they run in this order because rule 3 writes an acid slot rule 1 reads:

1. **peptide buffer capacity** (beer, D-238) — re-rooted on the member's own map whenever the
   scenario did not name `peptide_buffer_gpl`;
2. **the `cation_charge` pH anchor** (D-233), whenever the scenario gave an `initial_ph` — the
   compile-seam back-solve reads the pKa map and every `pKa_*` is sampled, so without it a member
   starts at a pH the scenario never asked for;
3. **wine's `copper` seed** (D-236), whenever the scenario did *not* name `copper_gpl` — the slot
   is seeded from `copper_typical`, which is also `PhenolicBrowning`'s mean-centring reference, so
   without it `f(Cu) == 1` holds at the nominal draw alone;
4. **the D-45 fallback seeds** (D-241) — the `_SEED_FALLBACKS` table in `scenario/compile.py`
   (`dms_potential`, `bound_h2s`, `bound_methanethiol`, `burst_antioxidant`, beer's wort `o2`)
   plus `must_fermentable_fraction`'s `S[0]`. Each fires only while the compiled slot still holds
   the parameter's own value, which refuses a scenario-stated level, stops silently if the seam
   stops deriving the slot, and reproduces D-147's burst-wiring gate without knowing it exists.

It rebuilds those slots and nothing else — a full re-run of the initial builder would move seeds no
beat has measured, and the scenario-input axis D-24 excluded (Brix, YAN, a stated `copper_gpl`) is
unchanged.

**`CompiledScenario.seed_reads` is the sampling half of rule 4** (D-241). The sampler scopes itself
by `Process.reads` (`_schedule_reads`), and a parameter the compile seam consumed is read by no
Process at runtime — so it was invisible to every band the engine published. `seed_reads` is a
`tuple[str, ...]` **derived from the rules themselves**, threaded to `simulate_ensemble` by
`run_ensemble` and unioned in `_resolve_sample_names`: into the **default** branch only (`only=`
stays exact) and **before** `exclude` (a seed can still be pinned). Because the draw and the
re-seed come from one list, a name is drawable iff a rule re-seeds it. It carries **no tier claim**
— `reads` has two masters (D-160) and this channel deliberately takes only the sampling one; what
a seed's tier should do to a state slot is unmeasured.

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
- **`phosphate` is in BOTH registries (D-210)** — the one shared member, and shared because it is
  the same species in the same role: the counter-anion of an `add_dap` dose, which is a verb either
  medium can be given. Nothing seeds it, so it is 0 unless dosed. Diprotic on purpose: pKa₃ = 12.35
  is beyond both the beverage and `titratable_acidity`'s endpoint, and `AcidSpec.protons` is what
  TA subtracts from, so carrying it would invent an equivalent per mole rather than merely idle.
  **Not the malt phosphate D-178 refused** — that one is present at t=0 and the inverse anchor
  absorbs its near-constant charge; a dose lands after the anchor and is permanent.

The `Byp` pool is read as a succinic-equivalent acid (zero new carbon, so `total_carbon` is
unchanged). Scalar `ph_of_state` / `degassed_ph_of_state` / `titratable_acidity` are pure and live
in core; the trajectory-series helpers need a `Trajectory` and therefore sit one layer up in
`fermentation.analysis`.

**The `N` slot is on the CATION side (D-209).** The assimilable-nitrogen pool is not electrically
neutral — it is ammonium plus amino acids whose side chains are charged at fermentation pH — so
`acidbase.nitrogen_charge_molar` adds `z̄ · [N]` mol⁺/L to the cation term, `z̄` being the mean
charge per mole of **elemental** nitrogen (`nitrogen_uptake_charge_beer` / `_wine`, both in
`acidbase.yaml` and both in `PH_SYSTEM_READS`). Consequences worth knowing before touching the
module:

- **This is where beer's pH gets its fall.** Nitrogen uptake is what makes the balance move, so
  the term is the difference between reaching about half of beer's measured acidification and
  reaching it (D-208 → D-209).
- **It is a re-allocation at t=0, not an addition.** `cation_charge` was back-solved from
  `initial_ph` and so already contained this charge, lumped and frozen; the slot now holds the
  remainder and all three anchoring sites subtract the nitrogen term. Anchored pH and must TA at
  t=0 are unchanged.
- **`_cation` therefore takes `params`**, required and not defaulted, on the same reasoning that
  keeps `carbonic_molar` positional — an omitted nitrogen term would be invisible otherwise.
- **It rides D-179's opt-in gate.** `charge_balance_is_populated` must be true, or an un-anchored
  beer's empty balance would get cation charge with no acid to meet it and read ~11 instead of 7.
- **Three of that pool's species are now on the ANION side too (D-239).** Peyer §5.5 puts
  aspartate (pKa 3.86), glutamate (4.25) and histidine (6.04) at ~10 % of a wort's buffering
  capacity, and beer carried none of them — the back-solved `peptide_buffer` lump absorbed their
  share, and the lump is permanent while the amino acids are eaten. `acidbase.AMINO_BUFFER_SPECS`
  is the **third "include-by-reading" entry** after `Byp` (D-18) and `carbonic` (D-182): keys in
  the pKa map that belong to no acid registry and hold no state slot. Their concentration is
  `ratio · [N]` (`wort_*_per_n`, all in `PH_SYSTEM_READS`), so the pool drains with uptake and no
  nitrogen is booked twice. `z̄` above absorbed the split — it now carries the three at their
  fully-protonated charge, and the two halves cancel at wort pH to D-209's own 0.1772 exactly.
  **Beer-only**, measured not preferred: the same three in a must are worth 0.73 % of wine's acid
  buffering against beer's 6.7 %, and wine speciates its amino acids as slots already (D-100).
  `_totals_molar` therefore takes `params`, required, for the reason `_cation` does.
  Cost: 0.000 pH at t=0, **+0.0023 at day 1** (the pool is still 70 % present), −0.0202 at day 7,
  which takes the high edge of `nitrogen_uptake_charge_beer`'s band **outside** Tyrell's day-7
  envelope. That is a priced worsening of agreement in exchange for fidelity, not a regression.
- **It moves wine's SO₂ readout**, which is the non-obvious downstream reach: the −0.097 pH costs
  an anchored wine, molecular SO₂ rises **0.228 → 0.284 mg/L (+24.7 %)** while free SO₂ moves
  0.5 % — the pool is unchanged, its speciation is not. Anything scored against a molecular-SO₂
  threshold reads differently since D-209.
- It is the **charge half only**: no medium models the pool per species, so uptake removes its
  charge but not its buffering, and that omitted half pushes the same way. A lower bound.
- **`z̄` is a must/wort composition average, and a DOSE carries its own charge (D-210).**
  `add_dap` doses pure ammonium, whose `z̄` is exactly +1 — 3× wine's average — and the
  `nitrogen_charge_excess` slot records how far the pool sits above the average, re-mixed at each
  dose by `remix_nitrogen_charge_excess`. One dimensionless slot rather than a second nitrogen
  pool, because a charge-per-mole is invariant under proportional drawdown: only an addition of
  differently-charged nitrogen moves it, so no Process touches it. The stored quantity is the
  *excess* and not the mean charge so that 0.0 is both the default and the correct undosed value —
  a mean-charge slot would need a sentinel, and a sentinel compared against a state float is a gate
  `num_jac` straddles. The eight Processes that *add* to `N` keep the average, measured: ~88 % of
  that inflow is an un-draw (`AminoAcidAssimilation`, where composition is unchanged) and the rest
  are deaminations whose true charge is 0 or +1 depending on the keto acid's fate, bracketing the
  average at ≤0.033 pH transiently and ~2e-5 at the endpoint.
- The **negative-slot guard** is `cation_slot_after_nitrogen`. The subtraction runs after
  `solve_cation_charge`'s own negativity check and so escapes it; a high-YAN, low-acid,
  low-`initial_ph` must could ship a negative cation slot while its anchor round trip still
  passed. All three sites route through the guard. It raises
  `NitrogenExceedsCationDemandError`, a `ValueError` subclass, because `set_ph` after a dose
  reaches it too and `_verb_set_ph` used to rewrite every failure as "below the acid load's
  intrinsic pH" — blaming an acid load that was not the cause (D-210).

**Two pH frames, and which one a caller wants is not a detail (D-208).** `ph_of_state` is the pH
*inside the vessel*, dissolved CO₂ included — the only frame a rate may read, and what every
pH-reading Process does read. `degassed_ph_of_state` is the same balance with that term zeroed,
which is the frame *published* beer/wine pH values are measured in (Analytica-EBC 9.35: "pH at
20 °C of decarbonated beer"), and it exists for scoring against literature only. The gap is ~0.29
pH on a day-7 beer and ~0.0007 on a wine. `titratable_acidity` already excluded the term for the
same reason.

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
- **Benchmark curves** — the §2.2 acceptance criteria are encoded as `BenchmarkSpec` data, gated
  behind the `benchmark` pytest marker (`uv run pytest -m benchmark`). The wine and CO₂ criteria
  **pass**; the **beer attenuation criterion is a strict `xfail` since D-221**, which
  re-temperatured it from 20 to 15 °C on Foster 2022's measured course — the handoff's 5-7 d
  duration is real at 15 °C and impossible at the 20 °C it was asserted at, and at the corrected
  temperature the engine takes 9.00 d. Its 20 °C pass was a criterion calibrated to a cooler
  ferment certifying a model that is ~1.5× too slow (D-220).
  `ReferenceSeries` + `compare_series` are the seam for scoring against real measured datasets.

## Testing & quality gates

`uv run pytest -n auto` (83 test files; unit, integration, conservation, sampling-surface and
doc-consistency checks), `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy`
(strict on `src`). CI runs all **four** on Python 3.13 and 3.14 — the format check is a separate
gate from the lint check and fails independently of it, which is how four consecutive commits
shipped red after D-239 (each run died in ~17 s at the format step, not in the ~13 min test
step, which is what made the failures easy to misread). Two of the test files guard documentation rather than physics:
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
