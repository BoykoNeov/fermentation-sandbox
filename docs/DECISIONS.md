# Design decisions

Lightweight decision log. Each entry: the decision, the rationale, and (where
relevant) how it deviates from the handoff brief. The handoff explicitly states
"nothing is rock solid"; this file records where we reasoned past it.

## Process decisions (project setup)

These three were the handoff's §7 open questions, resolved with the project owner.

- **D-A — Repository:** public GitHub repo `fermentation-sandbox` under
  `BoykoNeov`, licensed under the Boyko Non-Commercial License v1.0 (BNCL-1.0):
  free for non-commercial use/modification with attribution, commercial use
  prohibited unless separately licensed. (Originally MIT; relicensed by the
  copyright holder.)
- **D-B — First validation target:** chase the wine (~24 °Brix → dry in 10–14 d)
  **and** beer (~1.048 OG → ~1.010 in 5–7 d) §2.2 benchmarks **in parallel** for
  Milestone 1. The architecture is a shared core regardless; this sets which
  benchmarks gate the milestone. Consequence: `S` is a sugar *vector* from day one
  (see D-4).
- **D-C — Real datasets:** none available yet. Validate against published
  benchmark curves + qualitative directional checks now; the validation harness
  (`ReferenceSeries`, `compare_series`) is built data-ready so real time-series
  drop in later without rework.

## Engineering decisions

### D-1 — Tier metadata is derived, not carried inside state floats
**Decision:** the integrated state is a plain `float64` array. Confidence tier is
a property of `Process` and `Parameter` objects; an output's tier is *computed* at
the analysis boundary (`ProcessSet.tier_of`, `Tier.combine`).
**Why:** `solve_ivp` needs a contiguous numeric array; wrapping each scalar in a
tier-carrying object (as a literal reading of handoff §1.2 suggests) would wreck
the integration hot loop and complicate the math. Deriving the tier from
contributors still guarantees "the tier travels to every output" — the actual
prime directive — without the cost.
**Deviation:** reinterprets handoff §1.2 ("each scalar should carry its tier").
**Status (M1):** closed. `ProcessSet.tier_of`/`tier_map`/`overall_tier` now take an
optional `param_tiers` map and fold in the tiers of the parameters each
Process/modifier declares it `reads` (Process gained a `reads` attribute matching
`RateModifier`). A VALIDATED process running on a speculative parameter therefore
reports speculative — the credibility-borrowing this entry warned about is gone.
The runtime path carries it end-to-end: `simulate(..., param_tiers=...)` forwards
into `Trajectory.tier_map` (build the map with `ParameterSet.tier_map()`). Two
honesty guards: a declared `read` absent from `param_tiers` raises `KeyError`
rather than defaulting to validated; and `param_tiers=None` yields the *structural*
(Process/modifier-only) tier — still useful, but narrower, so reporting paths pass
the map. See `tests/test_process.py` (parameter-tier propagation) and
`tests/test_integrate.py::test_trajectory_tier_map_caps_on_param_tiers`.

### D-2 — Provenance enforced by schema, not convention
**Decision:** parameters load through Pydantic models that *require*
value/units/tier/uncertainty/provenance; a missing field raises at load time.
**Why:** the handoff says "no magic numbers, no exceptions," but plain YAML can't
enforce that. Making it a load-time error turns the rule into a guarantee.

### D-3 — SI-ish canonical internal units; convert only at edges
**Decision:** concentration g/L (≡ kg/m³), temperature K, **time in hours**.
Conversions (Brix/SG/Plato/ABV/°C/days) live in `fermentation.units` and are
called only at I/O boundaries. No `pint` quantities in the hot loop.
**Why:** matches handoff §7's "single canonical internal representation."
Kelvin because Arrhenius needs absolute temperature. Hours (not SI seconds)
because kinetic constants are overwhelmingly reported per-hour and benchmarks are
quoted in days — human-scale numbers, fewer transcription errors. Documented so
the deviation from strict SI on the time axis is explicit.

### D-4 — Sugar `S` is always a vector
**Decision:** even wine uses a length-1 sugar vector; beer uses length-3
(glucose, maltose, maltotriose).
**Why:** honours the handoff's "expansion = addition, not rewrite." With D-B
(both benchmarks in parallel) this is required, not just nice-to-have.

### D-5 — Scenarios are schema-validated YAML, not a custom DSL
**Decision:** use Pydantic-validated YAML/JSON for scenarios.
**Why:** the handoff offered "YAML/JSON or a small DSL"; a DSL is premature
complexity. YAML gives us sweeps, sharing, and validation for free.

### D-6 — Tooling: uv + pytest/hypothesis + ruff + mypy(strict on src)
**Decision:** `uv` for env/deps; `pytest` (+ `hypothesis` for property tests like
unit round-trips and conservation); `ruff` lint+format; `mypy --strict` on `src`,
relaxed signature requirements for tests.
**Why:** fast, modern, reproducible. Strict types on the library catch real bugs;
forcing `-> None` on every pytest function is noise, so tests are exempt from that
one rule while still being type-checked.

### D-7 — Media live in the core; the compile seam owns the unit boundary
**Decision:** a `Medium` (`fermentation.core.media`) bundles a beverage family's
`StateSchema` with the Processes that act on it, held in a `MEDIA` registry.
`compile_scenario` (`fermentation.scenario.compile`) converts a `Scenario` into a
`CompiledScenario` record (`y0`, `process_set`, `parameters`/`param_values`,
`schema`, `t_span_h`).
**Why / shape:**
- *Schemas in core, not scenario.* Processes (core) reference variable names like
  `"S"`/`"N"`; putting the per-medium layout in the core gives them and the
  compile seam one source of truth, and keeps `Medium` pure (no I/O). The
  industry-unit *conversion* stays at the boundary (scenario), honouring "convert
  only at edges" (D-3) — so the `scenario.initial` key vocabulary and all
  Brix/°C/days conversions live in `compile`, keyed by medium name (kept in sync
  with `MEDIA` via explicit guards).
- *A named record, not a bare tuple.* The brief wrote the seam as
  `(y0, ProcessSet, params)`; a frozen `CompiledScenario` with named fields is
  less fragile and also carries the `schema`, `t_span_h`, and the full
  `ParameterSet` (tiers/provenance) alongside the resolved `param_values` the hot
  loop needs. The function is `compile_scenario`, not `compile`, to avoid
  shadowing the builtin.
- *Beer sugars are explicit.* `compile` does **not** split a single OG into
  glucose/maltose/maltotriose — that wort spectrum is a provenance-backed
  parameter (the M1 sourcing task), so baking a fixed split into the seam would be
  a magic number. Until kinetics land, `process_factories` is empty and a compiled
  medium integrates to a constant baseline (verified by test).
**Status (M1):** schemas + seam done; Processes register into each `Medium` as
they are implemented. Both `wine_generic.yaml` and `beer_generic.yaml` now exist
with sourced parameters (D-12), so `wine`/`beer` + the default `generic` strain
compile without an override; an unsourced strain still raises a clear
`FileNotFoundError`.

### D-8 — Conservation scope: carbon (+ nitrogen) are the rigorous invariants; mass is scoped to the abiotic conversion
**Decision (what each balance covers):**
- **Carbon** is the primary rigorous invariant. `total_carbon` sums grams of
  carbon over `{S (per sugar component), E, CO₂, X}`. To make it close *exactly*,
  M1's sugar→ethanol+CO₂ kinetics use the **theoretical** Gay-Lussac split
  (`C₆H₁₂O₆ → 2 C₂H₅OH + 2 CO₂`), which is carbon- and mass-balanced by atom
  count. The realised-yield gap (literature ~0.46–0.48 g ethanol/g sugar vs the
  0.511 theoretical) is real chemistry — carbon diverted to **glycerol and
  organic acids** — but those byproducts are Tier-2 and not tracked in M1. So
  that carbon sink is **deferred**, not lost.
  - *Visible consequence:* the M1 model's realised ethanol yield reads slightly
    high (~0.49–0.50, near theoretical) until the glycerol Process lands. This is
    acceptable because **none of the three M1 benchmarks gate on absolute ABV**
    (`wine_dryness` = days-to-dryness, `beer_attenuation` = days-to-gravity,
    `co2_peak_then_tail` = a CO₂/sugar *ratio*). The realised-yield parameter
    (`Y_ethanol_sugar = 0.47`) stays in the store for when glycerol arrives.
  - *Biomass carbon is routed from sugar, with no anabolic CO₂ (M1).*
    `GrowthNitrogenLimited` draws the new biomass's carbon skeleton straight from
    `S` (`carbon(S) removed = biomass_C_fraction · dX`), so `total_carbon` over
    `{S, E, CO₂, X}` closes to machine precision under growth alone. Respiratory/
    anabolic CO₂ is **not** modelled, so the biomass yield-on-sugar is carbon-cheap
    (~0.82 g/g in isolation). This is immaterial for M1: nitrogen caps biomass near
    `X₀ + N₀/f_N` (~2–3 g/L for wine), so only ~1–2 % of sugar is diverted to
    biomass. *Consequence to revisit:* that 1–2 % carbon never appears as CO₂, which
    eats into the `co2_peak_then_tail` ±5 % budget — a tuning note for when that
    benchmark is unskipped, not a problem now. Because biomass pulls H/O from the
    solvent (D-8's biomass-mass point), `total_mass` over `{S, E, CO₂}` does **not**
    close once growth is active — carbon, not mass, is the invariant to assert on a
    growth run.
- **Nitrogen** is the second rigorous invariant: `total_nitrogen` sums free YAN
  `N` plus nitrogen bound in biomass (`biomass_N_fraction · X`). Conserved once
  the nitrogen-limited growth Process exists.
- **Mass** closes only for a single **hexose** (wine): `C₆H₁₂O₆ → 2 C₂H₅OH + 2 CO₂`
  is mass-balanced (`180.156 = 92.138 + 88.018 g/mol`), so `total_mass` sums
  `{S, E, CO₂}` and is conserved to solver tolerance there. It does **not**
  generalise, by the same untracked-solvent-H/O mechanism in two places: (a)
  di-/trisaccharide uptake *hydrolyses*, pulling water into the product pool —
  maltose adds ~5.3% mass, maltotriose ~7.1% — so `{S,E,CO₂}` mass is **not** a
  beer invariant; and (b) dry biomass draws H/O from the solvent, so whole-system
  dry mass over `{X,S,E,N,CO₂}` does not close (~1–2%) either. **Carbon is the
  rigorous cross-medium invariant** (water carries no carbon — 12 C in maltose, 12
  C out), so `total_mass` *rejects a multi-component sugar* and beer relies on
  `total_carbon`. This narrows the CLAUDE.md "carbon/nitrogen/mass must balance"
  line: carbon and nitrogen are the enforced **atom** balances across media; mass
  is the wine/hexose abiotic-conversion check. Recorded here so the scoping is explicit, not silent.

**Why / where constants live:**
- Stoichiometric constants — molar masses and carbon-atom counts of glucose /
  maltose / maltotriose / ethanol / CO₂ — are exact consequences of the chemical
  formulae, so (like the conversion factors in `fermentation.units`, D-3) they
  live in code with citations: `fermentation.core.chemistry`. Putting them in the
  core makes them a **single source of truth** shared by the conservation checks
  *and* the sugar-uptake Process, so a check can never disagree with the kinetics
  it audits. The toy test fixture derives its split from the same module for the
  same reason. The S-slot→species map (`chemistry.sugar_species`) lives here too,
  for the same single-source-of-truth reason and because the core kinetics that
  draw carbon from sugar cannot import the validation layer (one-directional
  dependency) — `conservation.py` imports it back rather than duplicating it.
- **Biomass elemental composition** (C-fraction ≈ 0.48, N-fraction ≈ 0.11 from the
  canonical `CH₁.₈O₀.₅N₀.₂` formula) is *empirical and uncertain* and is consumed
  by both the conservation check and the growth Process — so it is a **Parameter**
  (provenance store), not a code constant. `total_carbon`/`total_nitrogen` take the
  biomass fraction as a **passed-in argument** (the caller resolves it from the
  store) rather than importing the loader into the core/validation math; if a
  schema has an `X` variable and no fraction is supplied, the builder raises rather
  than silently under-counting (which would report a *false* violation).

### D-9 — Sugar uptake is biomass-catalysed (decoupled from growth), with smooth catabolite repression for beer
**Decision:** `SugarUptakeToEthanolCO2` (`fermentation.core.kinetics.uptake`) makes
the fermentative flux a function of *standing biomass*, not of growth:
`r = q_sugar_max · X · S/(K_sugar_uptake + S)` per sugar slot. It is a separate
Process from `GrowthNitrogenLimited`, summed by `ProcessSet`. For beer's multi-sugar
`S`, slots are consumed in preference order via a **smooth** repression factor
`Π_{j<i} K_repression/(K_repression + S_j)` (each higher sugar suppressed while a
more-preferred one remains).
**Why:**
- *Decoupled from growth (not Pirt-style `q = μ/Y + m`).* Growth shuts off when YAN
  runs out (Monod on `N`), but most ethanol in a real primary ferment is made by
  *non-growing*, nitrogen-starved cells. A growth-coupled uptake would stall at high
  residual sugar the instant nitrogen ran out — it could never reach dryness. A
  maintenance term `m·X` "fixes" that only by reintroducing an independent
  biomass-catalysed flux under another name, with a poorly-constrained coefficient.
  So uptake is biomass-catalysed outright. Consequence: biomass yield-on-sugar is an
  *emergent* ratio of the two rates rather than a dialled coefficient — immaterial
  for M1 (no benchmark probes biomass yield; only ~1–2 % of sugar is diverted to
  biomass, D-8).
- *Smooth repression, not a hard switch.* A threshold gate ("don't touch maltose
  until glucose hits zero") puts a kink in the RHS that the BDF solver dislikes
  (tiny steps / chatter). A smooth repression factor is the actual mechanism
  (catabolite repression) *and* keeps the derivative continuous, for a couple of
  extra lines. Relies on the `S` slot order being the preference order, which
  `beer_schema` defines. `K_repression` is kept small (~2 g/L placeholder) so the
  switch is sharp; wine (one slot) never represses.
- *Theoretical Gay-Lussac yields.* Ethanol/CO₂ yields come from
  `chemistry.ethanol_yield`/`co2_yield` (theoretical 0.511/0.489 per hexose,
  generalised to di-/trisaccharides by `HEXOSE_UNITS`), **not** the realised
  `Y_ethanol_sugar = 0.47`, so carbon (wine+beer) and mass (wine) close exactly.
  This is the D-8 carbon-first scoping applied to the kinetics; `Y_ethanol_sugar`
  stays the Tier-2 glycerol-diversion hook, deliberately unread in M1.
- *Guards mirror `GrowthNitrogenLimited`.* Each `S_i` is clamped to ≥0 before it
  enters a Monod term or a repression denominator, and the Process returns zeros
  when `X ≤ 0` — without the clamp a negative solver excursion flips the uptake sign
  and *creates* sugar (and drives E/CO₂ negative), failing the carbon check.
**Consequence for the next task (`EthanolInhibition`):** `ProcessSet` is purely
*additive*, so ethanol inhibition cannot "multiply onto" uptake as a separate summed
Process. It must live either inside the uptake rate or in the modifier-hook
mechanism the `ArrheniusTemperature` task introduces. Uptake's rate computation is
kept isolated so a multiplicative modifier can wrap it. No inhibition is modelled
yet, so an M1 uptake-only run ferments to complete dryness.
**Resolved in D-10:** the modifier hook was built here (one task early), with
`EthanolInhibition` as its first consumer, wrapping uptake's *whole contribution* at
the `ProcessSet` level — so uptake needed no refactor after all.

### D-10 — Rate modifiers: multiplicative mechanisms scale a Process at the ProcessSet level
**Decision:** mechanisms that *scale* an existing flux rather than *add* one
(ethanol inhibition now; Arrhenius temperature next) are `RateModifier` objects, not
`Process` objects. A `RateModifier` declares `name`, `tier`, `modifies` (names of the
Processes it scales) and `reads`, and returns a scalar `factor(t, y, schema, params)`.
`ProcessSet` evaluates each active modifier's factor once per RHS call and multiplies
it onto the *entire contribution vector* of every Process it targets, before summing.
`EthanolInhibition` (`fermentation.core.kinetics.inhibition`) scales
`SugarUptakeToEthanolCO2`.
**Why this shape:**
- *Multiplicative, so it cannot be a summed Process.* `ProcessSet` is additive (D-9);
  an inhibition term that "multiplies onto" uptake cannot be a peer summed into the
  same total. The modifier hook is the mechanism D-9 anticipated.
- *Scale the whole vector at the `ProcessSet` level → conservation is free and uptake
  needs no refactor.* Multiplying a conserving Process's complete `(dS, dE, dCO2)` by
  one scalar preserves every balance it respects (a uniformly slower carbon-neutral
  flux is still carbon-neutral), so the carbon/mass checks pass on an inhibited run
  unchanged. Wrapping at the set level (not inside uptake) leaves uptake untouched and
  unaware it is being inhibited — cleaner than the in-rate wrap D-9 literally
  described. The `strict` touches contract still holds (scaling zeros stays zero).
- *Togglable and tier-tracked like a Process (prime directive #3).* Modifiers share
  the Processes' name space and enable/disable machinery; a disabled modifier
  contributes factor 1 and drops out of tier derivation. `tier_of` caps a variable by
  the tiers of the modifiers scaling any Process that touches it, so a speculative
  modifier on a validated Process reports speculative — the same weakest-input rule,
  extended to the multiplicative path. (Parameter-tier propagation — capping by the
  tiers of the `reads` params, including a modifier's own `reads` — is now wired into
  `tier_of`; see D-1's M1 status.)
**Deviation from D-9:** D-9 said inhibition would live "inside the uptake rate or in
the modifier-hook the `ArrheniusTemperature` task introduces" — i.e. it assumed the
hook would arrive *with* Arrhenius. We build it one task earlier, here, with
`EthanolInhibition` as its first consumer; Arrhenius will *reuse* it (targeting both
growth and uptake) rather than introduce it. Recorded so the reordering is explicit.
**Functional form — Levenspiel/Luong "toxic power".** `f = (1 - E/E_max)^n` for
`0 <= E < E_max`, else `0`, with `E_max = ethanol_tolerance` (existing param, read as
a *wall*: the flux reaches zero there, matching its "viability collapses past
tolerance" provenance) and `n = ethanol_inhibition_exponent` (new speculative param).
`n > 1` (placeholder 2.0) makes the touchdown C¹-smooth (`f'(E_max) = 0`), avoiding
the derivative kink a raw `n=1` linear form would put in the RHS for the BDF solver —
the same smoothness argument as D-9's catabolite repression. `E` is clamped `>= 0` and
`f` clamped at `0`, so a solver excursion cannot amplify the rate (factor > 1) or flip
it negative (which would *create* sugar).
**Known tension (tuning-task item; does not block this task):** the *placeholder*
`E_max = 110` g/L sits below a 24 °Brix must's ~124-135 g/L final ethanol, so an
inhibited wine run *stalls short of dryness* — opposite of benchmark #1. This is a
parameter-sourcing problem (a high-alcohol must implies a high-tolerance strain;
sourcing will likely push `E_max` to ~140-150, above `E_final`, so the ferment
slows-then-completes), not a flaw in the form: conservation is unaffected (uniform
scaling), the benchmark is skipped, and the unit tests assert the *mechanism* (smooth,
monotone, in `[0,1]`, conservation-preserving, togglable), never
dryness-under-inhibition. `EthanolInhibition` stays out of the `MEDIA` registry with
the other kinetics until the full set lands.

### D-11 — Arrhenius temperature dependence: a per-rate, reference-anchored RateModifier
**Decision:** temperature dependence is `ArrheniusTemperature`
(`fermentation.core.kinetics.arrhenius`), a `RateModifier` reusing the D-10 hook (no
new mechanism). It is **parameterised per rate**: each instance names the Process it
scales and the activation-energy Parameter it reads. The wine config uses two —
`ArrheniusTemperature.for_growth()` (reads `E_a_growth`, scales
`GrowthNitrogenLimited`) and `.for_uptake()` (reads `E_a_uptake`, scales
`SugarUptakeToEthanolCO2`) — sharing one `T_ref`. The factor is reference-anchored:

```
f(T) = exp( -(E_a / R) · (1/T - 1/T_ref) )
```

**Why this shape:**
- *Reference-anchored, no separate pre-exponential `A`.* Normalising to `T_ref` makes
  `f = 1` there, so the *measured* rate constant (`mu_max` / `q_sugar_max`) is used
  unscaled at its calibration temperature; above `T_ref` the factor exceeds 1
  (faster), below it is < 1 (slower). The measured constants *already* encode
  `A·exp(-E_a / R·T_ref)`, so carrying a standalone `A` would double-book the
  pre-exponential and could silently disagree with the rate constant it multiplies.
  Only `E_a` and `T_ref` are parameters; `A` is deliberately **not** one. (This is why
  `milestone-1-context.md`'s "Arrhenius A + E_a per rate" becomes *E_a + T_ref* per
  rate in practice — `T_ref` plays `A`'s role, anchored to the rate-constant
  provenance.)
- *Per-rate, not one shared `E_a`.* Growth and fermentation are distinct processes
  whose temperature sensitivities are not guaranteed equal, so collapsing them onto one
  `E_a` would bake in an unjustified assumption rather than let the data decide (prime
  directive #1). (The M1 placeholders are set equal pending sourcing — the *separate
  parameters* are the point, not a guessed ordering; and "fermentation continues at low
  T" is D-9's nitrogen decoupling, a separate effect, not a temperature one.) The
  codebase had
  already committed to this: `E_a_growth` is a per-process parameter name and the
  context doc says "per rate". The task line's "targets *both* growth and uptake"
  describes the *mechanism*, not an instance count — two instances of a parameterised
  modifier still target both. So this is the established design, not a deviation.
- *Conservation is free; no clamp.* `exp` is always positive, so the factor scales a
  targeted Process's whole contribution vector by a single positive scalar — every
  balance is preserved. Unlike the wall-type inhibition form there is no regime where
  the factor could go negative, so (unlike D-10) **no clamp is needed**; a defensive
  one would be inconsistent noise. Under **stacking** (uptake is scaled by ethanol
  inhibition *and* Arrhenius) the two factors compose to one combined scalar on a
  conserving vector, so carbon/nitrogen still close exactly (pinned by a 4-modifier
  full-run test).
- *Where the gas constant lives.* `R` is a *universal physical constant* (SI-exact
  since the 2019 redefinition), not a stoichiometric one — so it lives in code with a
  citation **local to the arrhenius module**, not in `core.chemistry` (whose docstring
  scopes it to molar masses / carbon counts) and not in the provenance store (which is
  for empirical, uncertain quantities). Same code-with-citation rule as D-3/D-8.
- *`name` is per-instance.* `ProcessSet` enforces unique names across Processes *and*
  modifiers, so `name`/`modifies`/`reads` are set in `__init__` (`"arrhenius_growth"`,
  `"arrhenius_uptake"`), not as class attributes — the one structural departure from
  the `EthanolInhibition` template.
- *Reads `T` from state, not params.* The factor reads `T` from the state vector
  (Kelvin, D-3), so it is already correct for the non-isothermal temperature dynamics
  of a later tier. In M1 no Process drives `T`, so a run is isothermal and the factor
  is constant within it; its job is to make *different-temperature* runs differ in rate
  (the directional "warmer ferments faster" check the unit tests assert).
**Tier:** the Arrhenius law is textbook → the *mechanism* is **plausible** (like
inhibition/growth/uptake). The placeholder `E_a`/`T_ref` are **speculative**;
parameter-tier propagation (D-1) caps the scaled outputs at speculative accordingly.
**New parameters:** `E_a_uptake` (60 kJ/mol placeholder) and `T_ref` (293.15 K, the
20 °C the rate-constant placeholders are anchored to); `E_a_growth` retained. All
speculative. **Held out of the `MEDIA` registry** with the other kinetics until the
full set lands. Tests in `tests/test_kinetics_arrhenius.py`.

### D-12 — Parameter sourcing: keystone literature, honest mapping, honest tiers
**Decision:** the placeholder kinetic constants are replaced with literature values
(`wine_generic.yaml` rewritten, `beer_generic.yaml` added), each carrying a real DOI
where the value traces to text actually read in-source. Keystone sources:
- **Wine — Coleman, Fish & Block 2007**, *Appl. Environ. Microbiol.* 73(18):5875-5884,
  `doi:10.1128/aem.00670-07` (PDF read directly). Strain Premier Cuvée (= EC-1118 /
  Prise de Mousse, *S. cerevisiae* var. *bayanus*), Chardonnay must, 11-35 °C. Its
  model is structurally close to ours (growth Monod-on-nitrogen, uptake
  Michaelis-Menten on sugar), so `mu_max`, `K_n`, `q_sugar_max`, `K_sugar_uptake`,
  and the temperature sensitivity map onto our parameters.
- **Beer — Zamudio Lara et al. 2022**, *Foods* 11(22):3602,
  `doi:10.3390/foods11223602` (open-access CC-BY, Tables 5/6 read directly). Real ale
  fermentation, Grainfather pilot plant, 17-26 °C. Supplies `mu_max` and
  `K_sugar_uptake`; corroborates the realised yield.

**Three reconciliations worth recording (the task was reconciliation, not transcription):**
- *Coleman's "Log" is the natural log*, not base-10 — confirmed by the paper's own
  statement that `mu_max ≈ 0.05/h` at 11 °C matching `exp(-3.92 + 0.0782·11) = 0.047`
  (base-10 gives 0.0009). All Table A2 regressions are evaluated at **T_ref = 20 °C**
  (the wine benchmark temperature): e.g. `mu_max = exp(-3.92 + 0.0782·20) = 0.095/h`.
- *Equivalent Arrhenius `E_a` from a log-linear regression.* Coleman models
  temperature as `ln(rate) = a0 + a1·T(°C)`, **not** Arrhenius. Matching the local
  sensitivity `d(ln rate)/dT` of our `f = exp(-(E_a/R)(1/T - 1/T_ref))` to Coleman's
  slope gives `E_a = a1·R·T_ref²` → growth 55.9 kJ/mol (a1=0.0782), uptake 55.1 kJ/mol
  (a1=0.0771). Transparent derivation, tier `plausible`. (These are **inert at the
  isothermal M1 benchmark** — `f = 1` at `T_ref` — so they are Tier-2 readiness only.)
- *`q_sugar_max` is `β_max / Y_E/S`, not `β_max`.* Coleman eq 5 gives
  `dS/dt = -(β_max/Y_E/S)·[S/(K_S+S)]·X_A`; our uptake's rate is *sugar* consumed, so
  `q_sugar_max = β_max/Y_E/S = 0.469/0.550 = 0.85 g/g/h` (β_max alone, eq 4, is the
  specific *ethanol* rate). Sanity: `0.511 × 0.85 ≈ 0.43 g/g/h` ethanol ≈ Coleman's
  observed β_max. The value was **not** selected to hit the benchmark timing (the #4
  trap); the eq-5 coefficient match settles it.

**`ethanol_tolerance` = 142 g/L (wine)** comes from the Premier Cuvée / EC-1118
technical data sheet (18% v/v × 0.789). This is the *exact strain Coleman used*, so the
value is sourced independently of the benchmark; it happening to exceed a 24 °Brix
must's ~135 g/L final ethanol (resolving the D-10 stall) is a consequence, not the
selection criterion. Tier `plausible` **with the caveat in-file** that the Luong-wall
*form* is our modelling choice (Coleman instead uses an ethanol-driven death term);
the value maps cleanly (max ABV achievable ≈ E_max where rate→0).

**Honest tiers (prime directive #1; do not inflate):** only parameters that a source
measures *in our functional form* are promoted to `plausible`. Staying `speculative`
even after the sweep: `K_s` (Coleman growth is Monod-on-N only — no sugar term, so no
analogue for our growth-stage co-limitation guard), `K_repression` (form matches
Gee-Ramirez catabolite repression but no numeric constant was accessible),
`ethanol_inhibition_exponent` (de Andrés-Toro use n=1; our n=2 is a C¹-smoothness
choice). `Y_ethanol_sugar` stays at the well-established realised 0.47 — Coleman's
fitted 0.55 g/g *exceeds* the 0.511 theoretical maximum (a fitting/measurement
artefact) so it is not adopted.

**Beer is honestly thinner.** Published beer models are structurally further from ours
(Zamudio growth is Droop-like; de Andrés-Toro is Monod-on-*sugar*; neither is
nitrogen-limited), so beer values transfer by magnitude, not identity, and more stay
`speculative`: `K_n` is transferred from the wine fit (no beer model fits a nitrogen
constant), `q_sugar_max` is derived from Zamudio's growth-coupled `k_S·mu_max`, and the
beer `E_a`'s carry the verifiable Coleman-derived value rather than de Andrés-Toro's
**~35 kJ/mol** — which is consistently *reported* in secondary sources but whose primary
table (`doi:10.1016/S0378-4754(98)00147-5`, paywalled) was **not read in-source**, so
its DOI is *not* minted onto an unread number (the uncertainty range admits it).

**Deviation from context doc:** `milestone-1-context.md` lists Coleman as
`10.1128/AEM.00845-07`; the correct DOI is **`10.1128/aem.00670-07`** (00845-07 is a
different paper). Corrected here and in the YAML.

### D-13 — Ethanol brake: cumulative cell inactivation (two-pool) replaces the Luong wall
**Decision:** the validated core's ethanol brake is **ethanol-driven cell inactivation**
(Coleman 2007 eqs 2/7: `dX_A/dt = μ·X_A − k_d·X_A`, `k_d = k'_d·E`), implemented as the
`EthanolInactivation` Process. It **replaces the Luong wall** (D-10) in the default
`wine`/`beer` media. `EthanolInhibition` is retained as an optional class (strain/study
use) but is no longer wired in — keeping both would double-count ethanol toxicity.

**Why the wall could not stay.** The Luong factor `(1 − E/E_max)ⁿ` is *instantaneous and
reversible*: it scales the present flux by the present ethanol, holds no memory, and
(for `E_max` below a 24 °Brix must's final ethanol) stalls the ferment short of dryness
forever. A wine's 10-14 day *timescale* is set by the **irreversible, cumulative** loss
of catalytic cells as ethanol kills them — a stateful integral of damage, not a function
of the instantaneous state. Only a cumulative mechanism both decelerates the tail and
still finishes.

**Two-pool representation (chosen over a φ viability-fraction).** `X` stays the *viable*
biomass it always was (growth and uptake are catalysed by `X`); an inactivated pool
`X_dead` is added. Inactivation moves mass `X → X_dead` at equal rate (`r = k'_d·E·X`).
Because both pools carry the *same* elemental composition, the transfer is **carbon- and
nitrogen-neutral by construction** — a gram leaving `X` arrives in `X_dead` with the same
`f_C`/`f_N`, so `total_carbon`/`total_nitrogen` (which weight both pools) are untouched by
death. A φ-fraction folded into `X` would have made the conservation checks read a
shrinking carbon pool as mass destruction. Tier `plausible` (sourced mechanism, not yet
validated against our own curves).

**`k'_d` sourcing — and a published-typo correction (Coleman Table A2).** `k'_d` is the
only *quadratic* Coleman parameter: `ln(k'_d) = a0 + a1·T + a2·T²` (T °C). Table A2 prints
the `a1` **mean** as `−1.08×10⁻³`, but its printed 95 %
credible region is `[−1.94×10⁻¹, −3.30×10⁻²]` (centre `−1.13×10⁻¹`, half-width `8.1×10⁻²`).
The corrected `−1.08×10⁻¹` sits **essentially at that centre** — where a near-symmetric
posterior's mean belongs — whereas the as-printed `−1.08×10⁻³` lands ~1.4 half-widths
**beyond the upper bound, on the opposite side** of the interval. The journal typesetting
dropped the `×10ⁿ` exponent from the `a1` mean column; the true value is `−1.08×10⁻¹`.
Three independent checks confirm it: (1) it reproduces the paper's stated
**~13× rise** in `k'_d` over 11→35 °C (the as-printed value gives 191×); (2) it keeps
`k'_d(35 °C) = 4.4×10⁻⁴` under Fig 3b's `6×10⁻⁴` axis (the as-printed value overshoots
to `1.8×10⁻²`, ~30× off-scale); (3) the **identical defect** appears in the `Log(Y_X/N)`
row (printed `a1` mean `−3.61` vs CR `[−4.35×10⁻³, −2.93×10⁻³]`; corrected `−3.61×10⁻³`
reproduces Fig 4) — a systematic fault, not a one-off. Corrected value at 20 °C:
`k'_d = exp(−9.81 − 0.108·20 + 0.00478·400) = 4.28×10⁻⁵ (g/L)⁻¹h⁻¹`. The as-printed
`3.64×10⁻⁴` stalls *Coleman's own* model at ~108 g/L residual; the corrected value
reproduces his Fig 6c completion at 20 °C. M1 is isothermal at 20 °C so no Arrhenius
modifier is attached to `k'_d` (the quadratic does not reduce to a single `E_a`).

**Two gaps left open (deliberately, both separate tasks).** With the corrected `k'_d` the
wired wine model **completes** (S → 0), but: (a) it dries in **~7.7 days, below the 10-14
day benchmark window** — at the time read as an *uptake-speed* gap (β_max/biomass).
**Superseded by D-14: that was a misdiagnosis.** At its conditions (250 mg/L YAN, ample
nitrogen) ~7.7 d is *correct* — Coleman calls anything over 7-10 d "sluggish," and our
engine matches his own model to ~10 % there. The 10-14 d window is the *nitrogen-limited*
regime; the real gaps were a benchmark fixture that wasn't N-limited and a missing
N-dependent biomass yield (see D-14). (b) ABV lands at **16.9 %** (E ≈ 133 g/L) from the
theoretical Gay-Lussac split — the realised-yield/glycerol-sink task. Neither was folded
into the `k'_d` decision.

## D-14 — Nitrogen-dependent biomass yield; the wine benchmark window re-anchored to Coleman

**Status: closed.** Task #7 ("calibrate the Fig 6c reconstruction") resolved — and it
overturned its own premise (the D-13 gap-(a) "uptake-speed gap").

**The reframe (evidence, not figure-reading).** A faithful re-implementation of Coleman's
comprehensive model (eqs 1-8, Table A2 @ 20 °C — the model the paper validates against the
measured Fig 6c curves) reproduces our engine **line-for-line** on biomass and sugar at
*both* 80 and 330 mg N/L (X, S within ~2 % across 12 days; tracked in
`tests/test_coleman_reconstruction.py`). Triangulated three ways: the reconstruction, our
engine, and Coleman's own text ("completion exceeding 7-10 days = sluggish/problem
fermentation"; midrange temperatures "reach dryness in the minimum amount of time"). So a
24-Brix/20 °C wine with ample nitrogen *should* finish in ~6-7 d; our engine is right, not
too fast. The 10-14 d figure was a **generic handoff heuristic**, never Coleman.

**The one real model gap — N-dependent yield.** Coleman Fig 4 / Table A2 show the
cell-mass-per-nitrogen yield `Y_X/N` is *not* constant: it rises sharply as initial YAN
falls (`ln Y_X/N = 3.50 − 3.61e-3·YAN_mgL`; nitrogen-starved cells are elementally
N-poorer, so a gram of N builds more dry mass). Our model used a **fixed** `Y_X/N = 1/f_N
= 8.77`, so at low nitrogen it built too little biomass and **stuck** (residual ~31 g/L)
exactly where Coleman finishes. Adopting Coleman's regression closes this. The `a1`
exponent carries the **identical published typo as `k'_d`** (D-13): printed `−3.61`, but
its credible region is `[−4.35e-3, −2.93e-3]`, so the true value is `−3.61e-3` (reproduces
Fig 4 at 80 → 24.8 g/g and 330 → 10.1 g/g).

**Where it lives — computed at the compile boundary (a deliberate new pattern).** In our
model all assimilated nitrogen enters biomass, so `Y_X/N = 1/f_N` identically; we therefore
**override `biomass_N_fraction`** (rather than add a separate yield the growth Process
inverts), preserving the single-source contract that keeps the nitrogen balance exact — the
`total_nitrogen` check reads the same per-run constant the growth Process does, so
`d/dt[N + f_N·X] = 0` regardless of `f_N`'s value. Unlike the temperature regressions
(pre-evaluated into the YAML at the fixed `T_ref`), this one's evaluation point is the
scenario's *initial nitrogen*, so it cannot be pre-baked: `compile_scenario` evaluates it
from the scenario's YAN and nowhere else. This puts a parameter *value* (not physics) at the
scenario boundary; `chemistry.py`'s charter explicitly excludes empirical/strain-dependent
quantities, so a documented compile-seam helper is the right home. `biomass_C_fraction`
stays fixed — biomass carbon is ~1 % of the sugar→ethanol flux, immaterial in M1 (the
growth Process docstring already scopes this). **Beer keeps its static `f_N`**: Coleman is a
wine model and there is no sourced beer `Y_X/N` regression, so the override is gated on the
regression coefficients being present (wine-only by construction, not by accident).

**Benchmark window re-anchored (the user's call on a guarded §2.2 spec).** Because the
validated core now reproduces the keystone source, the acceptance window should reflect that
source, not a generic heuristic. The wine fixture is anchored to Coleman's documented
conditions — **80 mg N/L (his low-N treatment), ~0.25 g/L pitch (25 g/hL, standard practice
and consistent with his Fig 2 inoculum of ~0.1-0.3 g/L)** — *not* tuned to the window; it
lands at **~9.2 d**. The `wine_dryness` window was lowered **10-14 → 8-14 d**: the floor
drops to the fast end of realistic pitching that the source supports (pitch is a real
~2.6-day lever at low N, so it is anchored to the source, never swept to fit), the sluggish
ceiling stays at 14. `tests/benchmarks/test_milestone1.py::test_wine_24brix_ferments_to_dryness_in_window`
is unskipped and passing.

**Beer is now unblocked** — it shared the same Coleman-framework parameters, so this had to
settle first.

## D-15 — Beer §2.2 benchmarks: apparent (ethanol-depressed) gravity, and q re-derived

**Status: closed.** The two beer acceptance criteria are live and passing:
`test_beer_1048_og_attenuates_in_5_to_7_days` and `test_co2_integral_tracks_sugar_consumed`.

**"1.010" is an *apparent* (hydrometer) gravity, not real extract — and that is load-bearing.**
A fermenting beer's hydrometer reads *below* the true dissolved-solids extract because the
ethanol present is lighter than water. A 1.048 OG ale that brewers call "FG 1.010" has a
*real* extract near 4.25 °P (~1.016); the 1.010 is the ethanol-depressed apparent reading. So
the model's `(sugar, ethanol)` state must be mapped to **apparent** gravity to be compared
against 1.010. We added the standard Balling/Tabarie relation `RE = 0.1808·OE + 0.8192·AE`
(degrees Plato) to `units/convert.py` (`real_to_apparent_extract`, `apparent_gravity`), cited
alongside the existing ASBC polynomials — it is a boundary unit conversion, not physics. This
is fidelity, not gold-plating: against a *real-extract* gravity the 1.010 target would demand
an unrealistic ~79 % real degree of fermentation; the apparent correction lets a realistic
**~66 % RDF all-malt wort** be consistent with 1.010.

**No new state or parameter — the unfermentable extract is implicit.** The model tracks only
fermentable sugars. Real extract at time *t* = `OG_extract − sugar_consumed(t)`, so the
unfermentable share is implicitly `OG_extract − S0` and never needs a state slot (it is
constant — inert to kinetics and conservation) nor a parameter. The wort spec lives in the
**test fixture** (sourced, like the wine benchmark hardcodes Brix/YAN/pitch): a 1.048 OG
all-malt ale, fermentable `S0 ≈ 88 g/L` of the ~125 g/L total extract (RDF ~70 %), sugar
spectrum glucose/maltose/maltotriose ≈ 15/62/23 % of fermentables (typical all-malt split),
YAN 200 mg/L and pitch 0.6 g/L (typical ale practice). `S0 ≈ 88 g/L` is the initial fermentable
sugar **measured in our beer source** (Zamudio Lara et al. 2022), not back-solved from 1.010 —
the discipline of D-14 applied to the wort. The wort finishes at apparent **~1.007**, well
*below* 1.010 (a ~3.5-point margin), so the crossing lands in the kinetic phase rather than at a
fragile asymptote where a small parameter nudge would flip the metric to a never-crossing `inf`.

**`q_sugar_max` re-derived 1.5 → 0.5 (still speculative).** At the old 1.5 a 1.048 wort
attenuated in ~2 d — far inside the 5-7 d window. The 1.5 came from Zamudio's
`k_S·mu_max = 15.3·0.098`, but that equates the **growth-coupled peak** flux with a sustained
**catalytic** rate. Zamudio's growth is Droop-like (`mu_X = mu_max(1 − S_min/S)`, sub-maximal
and declining as sugar falls), so `k_S·mu_max` is only a transient peak; our uptake is
*decoupled* (all biomass catalytic at `q`, no `mu` factor), whose realised-equivalent is a
factor ~3 lower, `q ≈ 0.5 g/g/h`. With that sourced `q` the run lands at **~5.5 d**, inside
the 5-7 d window. Stays **speculative**, uncertainty `[0.3, 1.5]` spanning the realised rate to
the growth-coupled peak. Beer's `q` is independent of wine's 0.85, so the green wine benchmark
is untouched.

**Honesty caveat — what this benchmark does and does not validate (recorded, not hidden).** The
two halves are not equally strong. The **endpoint** (apparent FG ~1.007, ABV ~5.8 %) genuinely
*falls out* of the sourced wort and the apparent-gravity mapping — that half is real validation.
The **timescale** is set by a *speculative* `q` chosen at the low end of its independently
derivable range, so the benchmark confirms `q ≈ 0.5` is *consistent with* 5-7 d, **not** that
the window emerges unforced: `q` is pinnable only to ~a factor of 2, so beer's timescale test is
a **weaker validation than wine's** — a plausibility check, consistent with D-12's "beer is
honestly thinner."

**CO2 benchmark — the measurable channel, with the biomass diversion made visible.** The
evolved-CO2 integral is compared to the Gay-Lussac CO2 predicted from sugar consumed, summed
**per species over all three slots** (so the maltose-2×/maltotriose-3× hexose factors are
exercised). The ratio is **0.977**, deliberately *not* 1.0: ~2-3 % of sugar carbon is routed
into biomass by growth (no anabolic CO2 in M1), so slightly less CO2 evolves than total sugar
consumed implies — the `[0.95, 1.05]` window accommodates exactly that diversion. This is the
§2.2 measurable-channel check, *not* the machine-precision carbon audit (that stays in the
conservation tests). The test also asserts the spec's qualitative shape with real kinetic
teeth: d(CO2)/dt rises to an interior peak then tails off.

## D-16 — Realised ethanol yield: an explicit glycerol/byproduct carbon sink, plus a must-fermentable-sugar correction

**The two gaps this closes.** Through D-15 a 24 Brix wine fermented to **ABV 16.9 %**
(E ≈ 134 g/L) — unrealistically high. Two distinct, independently-sourced effects were
missing (the open thread D-13 gap-(b) and D-14 flagged but did not fold in):
1. **Realised yield < theoretical.** Real ferments divert a few percent of sugar carbon to
   glycerol, organic acids and higher alcohols, so realised `Y_E ≈ 0.46–0.48` g/g, not the
   theoretical Gay-Lussac 0.511 the kinetics used.
2. **Brix overstates fermentable sugar.** `brix_to_sugar_gpl` treats *all* 24 Brix solids as
   fermentable hexose (263.8 g/L), but glucose+fructose are only ~90–95 % of ripe-must
   soluble solids (Ribéreau-Gayon 2006); the rest is acids/minerals/phenolics.

**Decision — source each effect, let ABV fall out (do NOT reverse-engineer `Y_E`).** The
realised-yield literature value (0.47) alone lands ABV at 15.7 %, *not* 14–15 %; forcing 14 %
by pushing `Y_E` to ~0.43 would sit below the literature **and** over-attribute carbon to
glycerol — the exact tuning D-14/D-15 refused. Instead both effects are sourced from *measured
quantities* and the ABV emerges:
- **Glycerol sink** — `Y_glycerol_sugar = 0.035` g/g (→ ~8.6 g/L, mid the 4–10 g/L dry-wine
  range; UC Davis Waterhouse Lab, Scanes 1998, Ribéreau-Gayon). **plausible** (magnitude well
  corroborated; the constant-fraction *form* is the simplification).
- **Minor-byproduct lump** — `Y_byproduct_sugar = 0.014` g/g (→ ~3.4 g/L succinic + acetic +
  2,3-butanediol + higher alcohols). **speculative** (a lump booked at one representative
  carbon fraction).
- **Must fermentable fraction** — `must_fermentable_fraction = 0.93` g/g, applied at the
  compile boundary so wine loads ~245 g/L not 264. **plausible** (Ribéreau-Gayon composition).

Result: realised `Y_E ≈ 0.482` (cross-checks the literature 0.46–0.48, **not** set to it),
**ABV ≈ 15.0 %**, glycerol ≈ 8.5 g/L, byproducts ≈ 3.4 g/L — all fallout, nothing fitted to a
target.

**Mechanism — fold the split into uptake's yields, not a competing flux.** A separate
glycerol Process would *add* sugar consumption and speed dryness toward the 8 d floor. Instead
`SugarUptakeToEthanolCO2` keeps the sugar flux `dS = −r` **unchanged** and scales the
theoretical ethanol/CO2 split by `(1 − f_C/c(species))`, depositing the diverted carbon into
two new state pools, `Gly` (carbon-accounted as glycerol C₃H₈O₃) and `Byp` (as succinic acid
C₄H₆O₄). The carbon placed in `Gly`/`Byp` **exactly equals** the carbon scaled out of
ethanol+CO2, so `total_carbon` (which now weights both pools) closes to machine precision for
*any* yields — algebra: `scale·c(species) + Y_gly·c(gly) + Y_byp·c(byp) = c(species)`,
identically, for hexose and di/trisaccharides alike.

**Togglable-off = validated core intact (prime directive 3).** Both yields **default to 0**,
and at 0 the Process *is* the theoretical Gay-Lussac core. So the byproduct diversion is a
parameter-gated speculative layer over a protected validated core: with it off, wine `{S,E,CO2}`
mass still closes exactly (`total_mass` is asserted only on a byproduct-off configuration);
with it on, glycerol/succinic are more reduced than the ethanol route and draw redox H/O from
the solvent (like biomass), so only **carbon** closes. **Beer carries both yields at 0** —
its sugar→ethanol stays theoretical and its CO2-ratio benchmark is byte-for-byte untouched.

**Where the fermentable fraction lives.** It is *must composition*, not yeast-strain kinetics,
but is resolved at the `compile_scenario` boundary like the D-14 nitrogen-dependent yield (its
evaluation is scenario-specific). It sits in `wine_generic.yaml` for now, flagged as a
must-constant that would need re-homing if a second wine strain file is added.

**Consequence to watch — the dryness window tightened (a finding, surfaced not tuned).** The
fermentable-fraction cut (264 → 245 g/L) plus slightly less ethanol (→ less inactivation → more
viable biomass) move days-to-dryness from **9.2 d to 8.33 d** — still inside the D-14 `[8, 14]`
window, but with thinner margin. This is reported, **not** tuned away: per D-14 the engine
matches Coleman's own model line-for-line (the reconstruction test now feeds Coleman the same
fermentable S₀, so it still tracks to RMSE ~1.3 g/L), so 8.33 d is the *correct* consequence of
sourced inputs. If a future change breaches 8 d, the question is whether the heuristic window or
the fraction needs re-examination against Coleman — not whether to nudge a yield.

**New state plumbing.** `Gly`/`Byp` (and, retroactively, `X_dead`) are *produced-only* pools —
always 0 at pitch — so `VarSpec` gained a `default` and `StateSchema.pack` fills defaulted
pools when omitted; substrate/condition vars (X, S, E, N, T, CO2) stay required, preserving the
typo guard. This let two state variables land without touching ~37 initial-condition call sites.

## D-17 — Tier-promotion sweep: VALIDATED is reserved for independent data; the §2.2 pass earns PLAUSIBLE

**Status: closed.** The final M1 task — a sweep of every Process, modifier and
parameter now that all three §2.2 benchmarks pass, to decide what moves up.
**Outcome: promote nothing.** Recorded here with the evidence, because "promote
nothing" is itself the honest decision (the user's call on the VALIDATED bar), not a
skipped task.

**The bar — why §2.2 does not clear it.** VALIDATED means "established published
science *checked against independent benchmark curves*." The §2.2 pass is necessary
but not sufficient:
- **No measured time-series exist yet** (D-C). We validate against a published
  *model* (Coleman 2007) and *benchmark windows*, not raw experimental curves.
- **The wine window is re-anchored to Coleman** (D-14) — the same source the wine
  constants come from. Clearing a window derived from your own source is a
  faithful-implementation cross-check (a strong one — the reconstruction tracks to
  RMSE ~1.3 g/L), not *independent* validation.
- **Beer is explicitly weaker** (D-15): the attenuation timescale is set by a
  speculative `q_sugar_max` chosen at the low end of its range — a plausibility
  check by the source's own admission.

So passing §2.2 *confirms the PLAUSIBLE tier is earned* (sound mechanism, sourced
parameters, reproduces the keystone model) but VALIDATED waits for real curves to
drop into the data-ready harness (`ReferenceSeries`/`compare_series`, D-C). The
pre-registered "promote once §2.2 passes" language in the growth/uptake/inactivation
docstrings is rewritten to say this.

**Why the call is also low-stakes — promotion is inert at the output level (verified
on the real compile path).** Flipping growth/uptake/inactivation to VALIDATED and
re-deriving tiers changes **nothing** on the param-aware path (the D-1 real guarantee
that reporting uses) for either medium, and on the structural (`param_tiers=None`)
path moves exactly one variable — `X_dead` (plausible→validated). Wine flux outputs
are param-capped: `X`/`S` by `K_s`, and `E`/`CO2`/`Gly`/`Byp` by `K_repression` +
`Y_byproduct_sugar` (all speculative, D-12); the structural path is held at plausible
for every flux variable by the two Arrhenius modifiers (D-11). So the tier system
already reports honestly *regardless* of the mechanism-axis label — promoting the
Processes would have been a semantic statement about the forms, capped away at the
outputs anyway. This is parameter-tier propagation (D-1) and modifier-tier capping
(D-10/D-11) working as designed.

**Clean calls that hold regardless of the bar (the sweep's actual content):**
- **Arrhenius modifiers stay PLAUSIBLE** — inert at the isothermal `T_ref` benchmark
  (`f = 1`), so §2.2 never exercises them; an untested mechanism cannot be promoted.
- **Beer `q_sugar_max` stays SPECULATIVE** (D-15, the weaker beer timescale check).
- **`K_s`, `K_repression`, `Y_byproduct_sugar`, `ethanol_inhibition_exponent` stay
  SPECULATIVE** (D-12: no source measures them in our form; `K_s`/`K_repression` are
  inert guards for wine yet still cap conservatively — the design, not a defect).
- Everything already PLAUSIBLE (the Coleman/Zamudio-sourced constants; the three
  mechanisms) stays PLAUSIBLE — earned, not inflated.

**Future promotion trigger.** A parameter/Process moves to VALIDATED when it is
checked against an *independent measured* dataset for our own functional form — the
first such time-series to land in `ReferenceSeries` is the cue to revisit this sweep.

## D-18 — Tier-2 scope: pH is a charge-balance solver (derived-algebraic), byproducts are built first

**Status: RESOLVED (solver built 2026-06-30; see "Resolution" at the end of this entry).**
This opens Milestone 2 (Tier-2). It records two
calls made by the project owner at the start of Tier-2 — the pH-richness one is the
handoff's explicit "open decision for the human" (§7), the build order deviates from
the handoff's suggested sequence (§6). Detail in `docs/plans/milestone-2-*.md`.

**Call 1 — pH/acid is a full proton/charge-balance solver, not a tracked-pH
approximation (resolves handoff §7 open decision #3).** Each weak acid in the system
(tartaric, malic, lactic, acetic, ± carbonic) is tracked as a state concentration;
at each RHS evaluation the charge-balance equation `Σ(charged species) = 0` is solved
for `[H⁺]` given those totals and a pKa set, and `pH = −log₁₀[H⁺]` is read out.

**Why charge-balance — the discriminator is prime-directive #-level compositionality,
not accuracy.** A tracked-pH-with-drift can only produce the two couplings Tier-2
actually needs — MLF deacidification (pH rises ~0.1–0.3 as malic→lactic) and SO₂
speciation (molecular fraction governed by pKa ≈ 1.81) — by *scripting* the pH
response to each event. That directly violates "compositionality over scripting; never
hardcode the outcome of a specific additive/organism combination" (handoff §5). The
charge-balance solver makes both *emerge*: MLF consumes malic → recompute `[H⁺]` →
pH rises as a *consequence*; dose SO₂ → speciation falls out of the current pH. The
handoff also flags pH as "core infrastructure, not a byproduct — many Tier-2 mechanisms
are wrong without it" (§3.4). Cost, stated honestly: a pKa set + **per-acid initial
concentrations** become sourced scenario inputs (like Brix/YAN), and the acids become
**carbon-accounted state variables**.

**Corollary — pH is a derived algebraic pure function, NOT an integrated state.**
The derived-vs-integrated question is *not* a separate fork; it falls out of richness.
Charge-balance ⟹ there is no `dpH/dt`: pH is an instantaneous algebraic function of
the acid state (a 1-D monotonic root-find, well-behaved for the BDF RHS), keeping the
core pure exactly as `total_carbon` etc. are pure functions of state. (A tracked
approximation would instead have made pH an integrated state with a drift Process —
recorded so this is not re-litigated when the solver is built.)

**Three couplings the pH beat must resolve (named now so they are not discovered late):**
1. **Evolved vs dissolved CO₂.** The existing `CO2` state is the *cumulative evolved*
   measurable proxy (D-15), **not** the dissolved pool that carbonic acid needs. The
   solver must either add/track dissolved CO₂ for carbonic, or justify omitting carbonic
   for wine (tartaric/malic dominate must buffering) and document the scope.
2. **Acid carbon vs the D-16 `Byp` sink.** Tracked organic acids carry carbon, and D-16
   already books `Byp` as succinic (C₄). When acids become explicit state, `total_carbon`
   weighting and the `Byp` lump must be reconciled so carbon is not double-counted.
3. **pKa(T).** pKa is temperature-dependent; once byproducts/Arrhenius push runs off
   `T_ref` the constant-pKa assumption needs either a T-correction or an explicit scoped
   caveat.

**Call 2 — build byproducts/temperature first, then pH; deviation from handoff §6.**
The handoff sequence is "pH first (it unblocks the rest), then SO₂, then byproducts."
We invert the first two: the **temperature-/metabolism-driven byproducts** (§3.2 —
esters & fusels) are built before the pH solver. Rationale:
- *It closes the one remaining skipped benchmark* (`test_lower_temperature_is_slower_but_cleaner`),
  keeping the project's test-driven discipline — every prior decision was anchored to a
  §2.2 test.
- *It finally exercises the dormant temperature axis.* The Arrhenius modifiers were built
  in M1 but are **inert at the isothermal `T_ref` benchmark** (D-11, D-17), so the
  "warmer ferments faster" machinery has never been exercised by an acceptance test. The
  benchmark's *"slower"* half works **today** (a constant non-`T_ref` run activates them);
  only the *"cleaner"* half needs new ester/fusel Processes.
- *It is the most self-contained Tier-2 physics — esters/fusels depend on T and N only,
  not on pH.* So building it first costs the pH chain (SO₂/MLF/Brett, which *do* need pH)
  nothing, and defers the heavy charge-balance commitment until its design is locked.

The **stochastic ensemble wrapper** (handoff §1.6/§6.3 "runtime maturation") is
physics-free and orthogonal to both; it can be built in parallel at any point. Its API
shape is an engineering choice, not a scoping gate, so it carries no DECISIONS entry —
just `docs/plans/milestone-2-*.md`.

### Resolution (built 2026-06-30) — the solver, and the choices Call 1 left open

The charge-balance solver is `fermentation.core.acidbase` (pure core, `brentq` in
pH-space) + the `fermentation.analysis` series layer (top-layer sibling of `validation`,
imports `Trajectory`). pH/TA are derived **pure functions of state**, exactly as Call 1's
corollary requires — no `dpH/dt`. Deliverable scope: **solver + post-hoc pH/TA readout,
no RHS consumer** (SO₂/MLF wire pH into rates in later beats). The owner-confirmed calls
that the open entry above did not yet fix:

1. **Wine-only acid state.** D-18 acids are all wine acids (`tartaric`/`malic`/`lactic`
   state slots, appended to `wine_schema` only). Beer pH is a phosphate-buffered
   different acid system with no sourced data — explicitly **deferred**; `beer_schema` is
   untouched, and `ACID_STATE` extends to it when the data lands.

2. **A strong-cation term is mandatory, not optional.** Weak acids alone give pH ≈ **2.3**
   at must tartaric levels (~33 mM, pKa₁ ≈ 3.04); real must is ≈ **3.3**. K⁺ as bitartrate
   supplies the counter-charge — without it the solver is *qualitatively* wrong. It is
   carried as a net strong-cation charge density (`cation_charge` state slot, mol⁺/L).

3. **Anchoring = inverse (now).** The scenario gives acid concentrations + a measured
   `initial_ph`; compile **back-solves the strong-cation charge** (closed form,
   `solve_cation_charge`) to reproduce it, then stores it as a constant state slot; pH
   evolves emergently as acids change. Honest claim: **D-18 predicts pH *changes*, not
   absolute initial pH** (initial pH is an input). This folds activity-coefficient and
   cation uncertainty into one fitted term (how Boulton's wine-pH model is anchored). The
   back-solved cation lands in a physical K⁺ range (~25–50 meq/L, i.e. 1–2 g/L ÷ 39.1 —
   pinned as the unit-conversion guard test, since the round-trip is tautological w.r.t.
   the g/L↔mol/L factor). *Forward-from-cation is a documented future option* — the core
   solver is anchoring-agnostic and the cation stays a state slot, so adding a forward
   `cation_meq_l` input later is additive.

4. **Coupling #2 (acid carbon vs `Byp`) = include-by-reading.** The charge balance reads
   the *existing* `Byp` pool as its succinic-equivalent (`BYP_AS_SUCCINIC`) — **zero new
   carbon**, so `total_carbon` is unchanged and the double-count is *closed, not deferred*.
   The new `tartaric`/`malic`/`lactic` slots are weighted in `total_carbon` (so a future
   MLF Process, malic C₄ → lactic C₃ + CO₂ C₁, stays carbon-closing) but are inert in
   D-18 (no Process touches them ⇒ derivatives 0 ⇒ constant), so carbon still closes to
   machine precision. Caveat: `Byp` lumps neutral 2,3-butanediol, slightly overstating
   acid charge (~1–1.5 mM vs a ~20 mM buffer — minor).

**The four scope caveats, with numbers (justified scope, not hand-waves):**
- **Coupling #1 — carbonic omitted.** At pH 3.3 bicarbonate charge ~0.03 mM vs a ~20 mM
  buffer (~0.1 %); correct to omit below pH ~4. `CO2` state stays the evolved proxy
  (D-15). Revisit threshold: deacidified/low-acid musts above pH ~4.
- **Coupling #3 — constant pKa.** Carboxylic ΔH_ionization ≈ 0; the pKa shift over
  10–30 °C is <0.05 units, inside the pKa uncertainty. (We omit carbonic — the one acid
  with real T-dependence.)
- **Ionic strength / activity.** Wine I ≈ 0.05–0.1 M; concentration-based *apparent* pKa
  is the standard plausible-tier simplification, and inverse anchoring folds the activity
  error into the fitted cation at t=0, leaving it to affect only the *slope* (buffer
  capacity), where we claim only directional fidelity.
- **Tier = `plausible`, computed explicitly.** CRC pKa values are measured (validated),
  but applying 25 °C / I=0 constants to wine is extrapolation. `acidbase.ph_tier` computes
  the derived pH/TA tier as `combine(pKa tiers, PLAUSIBLE)` — it must NOT inherit the
  `VALIDATED` default `tier_of` returns for the inert acid slots no Process touches.

**Known TA-series artifact (scoped, not a solver bug).** `titratable_acidity` is exact
given its inputs, and the *must* (t=0) TA lands in the textbook 6–9 g/L band. But the TA
*series* **rises** ~3–4 g/L over a ferment because the whole `Byp` pool is read as
fully-titratable diprotic succinic and `Byp` accumulates to ~3 g/L (D-16/D-19). Real wine
TA is flat-to-*declining* during ferment (tartrate precipitation, malic metabolism), so the
end-of-ferment TA is an **over-estimate, directional only** — trust the t=0 value. The
cause is upstream pool sizing/booking (`Byp` lumps neutral 2,3-butanediol yet is booked
diprotic; the pool itself exceeds real succinic 0.5–1.5 g/L), bounded as *minor for pH*
(~1–1.5 mM vs ~20 mM buffer) but *direct and larger for TA*. Fixing it belongs upstream
(speciate `Byp`, re-source the pool), not in the D-18 solver.

**Acceptance gate (proof-of-purpose, met):** on a malic-rich must (tartaric 4 / malic 4
g/L, anchored pH 3.4) the full malic→lactic substitution raises pH by **0.225**, inside
the required MLF band [0.1, 0.3] — MLF-enablement demonstrated *without* an MLF Process
built. Second, emergent demonstration: with acids constant, the core `Byp` realised-yield
diversion grows 0 → ~2.9 g/L over a wine ferment, and include-by-reading makes its
succinate charge count, so the pH *series* drifts mildly **down** (3.40 → 3.33, ~0.067)
with the cation frozen at pitch — the solver responds to acid dynamics with no scripting.
This keystone unblocks **SO₂ → MLF → mixed cultures**.

## D-19 — Aroma byproducts (esters/fusels): carbon routed from sugar (option a1)

**Status: settled (the carbon-accounting sub-decision of the byproducts beat).** The
ester (`EsterSynthesis`) and fusel (`FuselAlcoholsEhrlich`) Processes and their trace
produced-pool schema slots landed earlier in the beat under **interim option (b)** —
pools *outside* `total_carbon`, touching only their own slot, carbon closure
byte-for-byte. This entry records the agreed end state: **option (a), variant a1 —
route ester/fusel carbon *from sugar* and weight the pools in `total_carbon`**, so they
are real carbon-accounted state under one rule with `Gly`/`Byp` (D-16), not diagnostic
re-expressions. Project owner's call (2026-06-29), over the advisor/author lean toward
(b) and the closure-neutral a2 variant.

**What a1 does.** Each byproduct Process draws its species' carbon *out of `S`*
(`_draw_carbon_from_sugar`, splitting the draw across sugar slots in proportion to each
slot's carbon content, so wine's 1 slot and beer's 3 are handled by one routine), and
`total_carbon` weights `esters` as ethyl acetate (C₄H₈O₂) and `fusels` as isoamyl
alcohol (C₅H₁₂O). The per-RHS carbon removed from sugar exactly equals the carbon
deposited in the pool, so carbon closes to machine precision.

**The draw touches only `S` — never `E`/`CO2`.** This is the surgical part. The uptake
Process still ferments `S` to ethanol+CO2 unchanged; the byproducts pull an *additional*
sliver of `S`. So at the derivative level only `dS` gains a term — `dX`/`dN`/`dE`/`dCO2`
stay byte-for-byte identical with the byproducts off. The integrated core therefore
drifts only by the trace sugar they consume (~0.2 % of `S0`).

**The `Byp` double-count, resolved (the hard part).** `Byp` formerly lumped "organic
acids + higher alcohols" (booked as succinic acid). Fusels *are* higher alcohols, so
weighting a separate carbon-routed `fusels` pool on top would book that carbon twice.
Resolution: `Byp` is re-anchored to **organic acids / polyols only** —
`Y_byproduct_sugar` (wine) reduced 0.014 → 0.012, removing exactly the higher-alcohol
share (~0.0017 g/g); the higher alcohols now live solely in the `fusels` pool. Beer
needs no carve-out (its `Y_byproduct_sugar` is 0, so nothing was double-booked).

**Two bookkeeping caveats — the carbon source is accounting, not metabolism.**
(i) The Ehrlich pathway builds fusels from *amino-acid* skeletons, but `N` (YAN) carries
no carbon in `total_carbon`, so fusel carbon is sourced from sugar as a stand-in.
(ii) An ester's ethanol moiety is carbon *already counted in `E`*, so routing ester
carbon from sugar over-attributes fresh hexose. Both close the ledger exactly; neither
claims where the carbon physically came from. Fusels carry **no CO2 co-product** (the
Ehrlich decarboxylation is omitted) — a documented simplification keeping the draw a
clean 1:1 sugar→pool carbon transfer.

**Tier consequence (noted, not user-facing).** Because the byproduct Processes now
touch `S`, `ProcessSet.tier_of("S")` folds in their tiers; the *structural-only*
(`param_tiers=None`) tier of `S` drops PLAUSIBLE → SPECULATIVE when byproducts are on.
The **param-aware** tier users actually see is *already* SPECULATIVE today (growth reads
`K_s`, uptake reads `K_repression`/`Y_byproduct_sugar` — all speculative), so a1 changes
nothing on the headline path. This is the intrinsic price of "real carbon-accounted
state" and is **not** an a1-vs-a2 discriminator (a2 would drag `E`/`Byp` down the same
way by touching them). Isolability (prime directive #3) holds structurally: the
validated core is the ProcessSet built *without* the byproduct tuple.

**Why a1 over (b)/a2.** (b) keeps the pools as unaccounted diagnostics — fine for
closure but it never lets `total_carbon` *include* the aroma carbon, and it relies on
the fragile claim "their carbon is booked elsewhere" (which the `Byp` overlap shows was
only half-true). a2 (transfer carbon out of `E`/`Byp` with no sugar draw) is
closure-neutral but, by not drawing from sugar, sits functionally next to the rejected
(b); its only edge is a smaller blast radius. a1 is the most physically literal and
gives one consistent rule for every produced-only pool — the project's fidelity bar.

**Empirical results (verified, not assumed).** Carbon closes to **1.1×10⁻¹³** on a full
wine ferment with byproducts on. The §2.2 realism guards are unmoved: wine **ABV
14.99 %**, realised **Y_E 0.482**, **glycerol 8.49 g/L**, **Byp 2.91 g/L** (the
`Y_byproduct_sugar` carve and the trace fusel/ester sugar draw nearly cancel on ABV).
Beer **CO₂/sugar-consumed ratio 0.975** (was 0.977; still inside [0.95, 1.05]). Wine
aroma totals ~0.11 g/L esters + ~0.05 g/L fusels (trace, as expected). 213 tests green.

**Scope note.** This is the *carbon-accounting* half of the byproducts beat. The
ester/fusel rate + `E_a` placeholders are now sourced (see the sourcing-step record
below); unskipping the directional benchmark `test_lower_temperature_is_slower_but_\
cleaner` remains the final step of the beat.

### D-19 sourcing step — ester/fusel rate + E_a placeholders (2026-06-30)

Replaced the four placeholder constants (`k_ester`, `E_a_esters`, `k_fusel`,
`E_a_fusels`, both media) with literature-bounded values + honest provenance. The
load-bearing constraint (each `E_a` > `E_a_uptake` = 55,100 J/mol so the run-integrated
aroma total rises with temperature) is held. Headline: **the E_a ORDERING is now
sourced; the rate MAGNITUDES and exact E_a values stay speculative (directional only,
handoff §3.5).** Values: `E_a_esters` 75,000 → **80,000 J/mol**; `E_a_fusels` **70,000**
(unchanged); `k_ester` **4.0e-4 /h**, `k_fusel` **2.5e-3 /h** (unchanged, order-of-
magnitude targets). All four stay **speculative**.

**Sources read (all open / provided in-source — none recalled).**
- *de Andrés-Toro et al. 1998* (Math. Comput. Simul. 48(1):65-74), the canonical beer
  byproduct model, read IN-SOURCE via the open **CC-BY** reproduction *Pilarski &
  Gerogiorgis 2022* (Processes 10(11):2400, doi:10.3390/pr10112400) **Table 1**, which
  transcribes its parameters verbatim. Ethyl acetate (an ester; exactly our
  `ethyl_acetate` booking species) forms as `dC_EA/dt = Y_EA·μ_x·X_A` — tied to the
  **growth** rate, with `Y_EA = exp(89.92 − 26589/T)` and `μ_X0 = exp(108.31 −
  31934/T)` (form `μ = exp(A + B/T)`, T in K; apparent `E_a = −B·R`): apparent
  `E_a ≈ 221` and `≈ 265 kJ/mol`. **No fusel/higher-alcohol term exists** in this model.
- *Mouret et al. 2015* (Biochem. Eng. J. 103:211-218, doi:10.1016/j.bej.2015.07.017) and
  *Rollero/Mouret et al. 2014* (Appl. Microbiol. Biotechnol. 99:2291-2304,
  doi:10.1007/s00253-014-6210-9) — the wine aroma analog (the actual "Mouret 2014/2015"
  reading-list items; **provided by the project owner** mid-task). MODAPEC parameterises
  aroma as two-phase production *yields from sugar*, **linear in T and N₀** (not an
  Arrhenius per-flux rate), via gas–liquid balances that separate synthesis from
  evaporation.

**Ordering vs magnitude — and why de Andrés-Toro's magnitude does NOT transfer.** Its
ester rides on **growth**, while its own **sugar-uptake** term `μ_S0` (A=−41.92,
B=+11654) has a **NEGATIVE** apparent E_a (≈ −97 kJ/mol — sugar uptake *falls* with T in
that fit). So its internal ΔE_a (ester − flux) is ~480 kJ/mol *within a model whose flux
E_a is negative* — incommensurable with our +55,100 J/mol Coleman uptake E_a. Lifting its
ester E_a and differencing against Coleman would splice two incompatible models. The
**ordering survives** the mismatch (ester ≫ flux, robustly, in a real fitted model — the
citation); the **magnitude does not**. So E_a is held GENERIC, beer-grounded, ~80 kJ/mol,
banded wide (60,000–250,000, all > E_a_uptake). `k_ester`/`k_fusel` are order-of-magnitude
targets (de Andrés-Toro's `Y_EA·μ_x·X_A` and Mouret's yield form give no constant
transferable to our flux-coupled `k`). Verified totals at dryness: wine **14 °C → 137.5,
20 °C → 165.3, 25 °C → 191.8 mg/L** total aroma (esters 114 mg/L at 20 °C — in the
50–200 mg/L band; fusels ~51 mg/L). Cleaner when colder. 214 tests green, ruff + mypy
clean; §2.2 trio + carbon conservation unmoved (E_a is inert at the 20 °C benchmark,
f=1 at T_ref; no k changed).

**WINE ESTER finding — surfaced, not buried (the important correction).** The primary
wine data *contradicts* a naive "warmer ⇒ more wine esters": Rollero 2014 states
**"evaporation largely accounted for the effect of temperature on the accumulation of
esters in liquid,"** and the *total production* (synthesis) our non-volatile `esters`
pool represents is **weak and non-monotonic in T** (isoamyl acetate quadratic, lowest
~24 °C; ethyl hexanoate ~T-independent). So **no value of `E_a_esters` reproduces wine
ester behaviour — the missing physics is a volatilization / gas-stripping sink the model
does not yet simulate.** We therefore: (i) struck the earlier "+~75% esters per 15 °C"
brewing-folklore magnitude anchor from the *wine* ester provenance (it is a beer/general
number that does not transfer); (ii) kept one GENERIC, beer-grounded `E_a_esters` >
`E_a_uptake` (de Andrés-Toro's beer coupling is real); (iii) documented the wine truth in
the `E_a_esters` note and the `byproducts.py` tier docstring. **Citing Mouret/Rollero as
supporting a wine ester rise would be false provenance — they show the opposite for
liquid and ~flat for synthesis.** For *wine* the warmer⇒more-aroma benchmark direction is
carried by the **FUSELS**, whose total-production rise with T *is* supported (Mouret
2015). **Free fusel corroboration:** Mouret/Rollero confirm higher-alcohol synthesis is
optimal at ~200–300 mg N/L and **non-monotonic in nitrogen** — exactly the simplification
`FuselAlcoholsEhrlich` flags as the reason it is speculative; now cited in its provenance.

**M1 correction (flagged, not silently rewritten).** The beer file's M1 Arrhenius notes
cited a secondary "de Andrés-Toro ~35 kJ/mol for growth and ethanol." The in-source
Table 1 debunks it: growth apparent E_a ≈ **265**, ethanol ≈ **10.5**, sugar ≈ **−97**
kJ/mol — none is 35, and all are extreme lumped empirical-fit artifacts (which *is why*
we carry the clean Coleman-derived value, not de Andrés-Toro's). The beer `E_a_growth`/
`E_a_uptake` **values and bands are unchanged** (M1 not silently rewritten); only the
notes are corrected. **Open item for the owner:** the beer band low (30,000) was
justified by the now-debunked "~35 kJ/mol beer figure" — it is retained pending a
deliberate M1-band review.

**Two items beyond this checkbox (recorded for the owner; NOT built here).**
1. *Volatilization / gas-stripping sink.* The real mechanism behind "cleaner when colder"
   for wine esters is evaporative loss (warm, vigorous CO₂ evolution strips volatile
   acetate/ethyl esters), which this model omits. A gas–liquid balance term (cf. Mouret's
   MODAPEC, Morakul et al.) is the principled fix — **future work.**
2. *Benchmark premise.* `test_lower_temperature_is_slower_but_cleaner` (the next, final
   step of the beat) assumes warmer ⇒ more esters AND fusels. The *combined* esters+fusels
   total still rises with T in both media (beer esters + both media's fusels carry it), so
   the directional benchmark is passable as written. But the *wine-ester* half of its
   premise is confounded by evaporation; unskipping it honestly for wine may want the
   volatilization sink first. **Owner decision point** before that checkbox.

### D-20 — ester volatilization (gas-stripping) sink; benchmark unskipped (2026-06-30)

**Owner decision: option (B).** At the D-19 decision point the owner chose to **build the
volatilization / gas-stripping sink first**, then unskip
`test_lower_temperature_is_slower_but_cleaner` honestly — rather than pass the benchmark
on the combined esters+fusels total (option A), which would have hidden the wine-ester
inversion D-19 surfaced. This closes the byproducts beat.

**What was built.** A new produced-only bookkeeping pool **`esters_gas`** (volatilized
esters in the headspace) and a Process **`EsterVolatilization`** that strips liquid
`esters` into it:

```
d(esters)/dt   -= k_ester_volatil · X·S_total/(K_sugar_uptake+S_total) · f(T) · esters
d(esters_gas)/dt += (same)        with f(T) = arrhenius_factor(T, E_a_ester_volatil, T_ref)
```

It rides the **same fermentative-flux proxy** as the CO₂ evolution that does the stripping
(`_fermentative_flux_shape`), is **first-order in the liquid ester present**, and **stops
when fermentation stops** (`flux → 0` at dryness — a deliberate omission of slow passive
post-ferment evaporation, keeping the sink a clean function of the gas stream). Esters-only:
isoamyl alcohol (bp ~131 °C) is far less volatile than ethyl acetate (~77 °C), so fusels
stay the warmer⇒more-aroma carrier (Rollero 2014).

**Carbon — a neutral liquid→gas transfer (no sugar draw).** Unlike `EsterSynthesis`/
`FuselAlcoholsEhrlich` (which draw fresh sugar, a1/D-19), this Process moves carbon already
in `esters` into `esters_gas`, both booked as ethyl acetate. It touches `esters`/`esters_gas`
only — never `S`/`E`/`CO2`. `total_carbon` weights `esters_gas` at the same ethyl-acetate
fraction (the ester analogue of how evolved `CO2` stays counted: carbon leaves the liquid,
not the ledger), so closure stays at **machine precision** while wine's liquid esters
honestly fall with T. `esters` is clamped ≥ 0 so a solver undershoot can't strip a negative
pool into spurious gas.

**The per-medium E_a split (the load-bearing parameterisation, and the trap avoided).**
Near quasi-steady-state `[esters] ∝ f_synth(T)/f_volatil(T)` — the shared flux cancels, so
the *direction* is set purely by which activation energy is larger. With `E_a_ester_volatil`
sourced **per medium** (separate YAMLs), both directions are honest, captured by the E_a
balance not by two code paths:

| medium | `E_a_ester_volatil` vs `E_a_esters` | net liquid-ester direction | source |
|--------|-------------------------------------|----------------------------|--------|
| wine   | **above** (130k > 80k)              | **falls** with T (inversion) | Rollero 2014 — "evaporation largely accounted for the effect of T on liquid ester accumulation" |
| beer   | **below** (40k < 80k)               | **rises** with T            | de Andrés-Toro 1998 — ester rides the strongly-T-sensitive growth rate; warm ales are estery |

A *single global* stripping E_a above `E_a_esters` would have silently inverted **beer**
too (breaking the sourced warm-ale expectation) — the trap the per-medium split avoids.

**Honesty caveat on the wine magnitude (not buried).** A pure volatility/Henry's-law Q10 is
~2–3 (E_a ≈ 50–75 kJ/mol), which is *below* `E_a_esters` and would **not** invert on its
own. The model's ester *synthesis* (`E_a_esters` = 80k, generic-beer-grounded, monotone-
rising) is almost certainly **too T-sensitive for wine** (Rollero: wine ester synthesis is
weak/non-monotonic), so `E_a_ester_volatil` is set above it to reproduce the **net observed
liquid inversion given the rest of the model** — a lumped, *compensating* value (Q10 ~5.6),
not a first-principles Henry's constant. All four volatilization params stay **speculative**;
only the per-medium *ordering* relative to `E_a_esters` is sourced and load-bearing. This is
documented in the `E_a_ester_volatil` provenance note in both YAMLs.

**Empirical results (verified at 14/20/25 °C, carbon closing to machine precision each run).**
- *Wine* liquid esters **54 → 45 → 35 mg/L** (fall with T); volatilized `esters_gas`
  **39 → 69 → 101 mg/L** (rise — the stripped fraction); fusels **45 → 51 → 56 mg/L** (rise).
  Total *produced* (liquid+gas) still rises with T (synthesis), as claimed.
- *Beer* liquid esters **57 → 72 → 87 mg/L** (rise with T); fusels **37 → 41 → 46 mg/L**.

**Benchmark, rewritten honest per medium.** `test_lower_temperature_is_slower_but_cleaner`
is **unskipped** and asserts, reading the **liquid** pools only (the `esters_gas` headspace
is not aroma in the glass): both media slower-to-dryness + fewer **fusels** when colder
(the real "cleaner"); **beer** fewer liquid esters when colder; **wine** *more* liquid
esters when colder (the inversion). Asserting a combined total would hide the inversion the
sink was built to surface, so each pool's sourced direction is asserted explicitly. The unit
test `test_integrated_byproduct_total_falls_with_temperature` (which encoded the old
combined-total premise) is replaced by `test_integrated_wine_aroma_temperature_directions`
with the same per-pool checks as the E_a-ordering regression guard.

**Scope / impact.** Schema grows 11→12 (wine) and 13→14 (beer). §2.2 trio unmoved (the sink
is inert at the 20 °C benchmark relative to the bands; it moves carbon between two trace
pools, never touching `S`/`E`/`CO2`). Isolable (prime directive #3): `EsterVolatilization`
lives in the `_BYPRODUCT_PROCESSES` tuple, so the validated core is still the ProcessSet
built without it. **222 tests green** (was 214; +8 net incl. the now-live benchmark), ruff +
format + mypy clean.

**Still future work (recorded, not built here).** The flux-coupled stripping is a stand-in
for a full gas–liquid (Henry's-law) balance (cf. Mouret's MODAPEC, Morakul et al.); a
principled model would carry the partition coefficient explicitly and let passive
evaporation continue after the cap goes on. The `esters_gas` pool is the hook for that.

### D-21 — physical Henry's-law stripping + per-medium sourced synthesis E_a (2026-06-30)

**Owner decision: build the full Henry's-law balance (the rigorous option), then confirm
the unified build by prototyping.** This supersedes D-20's *parameterisation* (the
mechanism, gas pool, carbon bookkeeping, and benchmark structure from D-20 all stand);
what changed is *why* the wine/beer directions diverge and *which* parameters carry them.

**The reconcile that reframed it (advisor, then verified).** D-20 made the wine/beer
ester-direction split by *fudging the stripping* `E_a_ester_volatil` per medium (wine 130k
above `E_a_esters`, beer 40k below). But a **sourced Henry's-law stripping is a property of
the molecule, not the beverage** — the same partition K_H(T) in wine and beer (Morakul et
al. 2011 explicitly: the partition coefficient depends only on composition and temperature).
So a physical stripping *cannot* push opposite directions by itself; using one would have
**silently inverted beer too** (warm ales must stay estery). The direction therefore has to
live where it is genuinely sourced — in ester **synthesis**, which differs by medium in the
literature: beer strongly T-sensitive (de Andrés-Toro 1998, ester ride the growth rate,
apparent E_a ~221–265 kJ/mol), wine weak/non-monotonic (Mouret 2015; Rollero 2014). The two
options put to the owner (Henry's-law vs per-medium synthesis E_a) were thus **one build**.

**What changed in the model.**
- `EsterVolatilization` now reads `E_a_uptake` (gas-flow factor — the stripping rides the
  same Arrhenius-scaled fermentative flux as the CO₂ it travels on) and a new
  `dH_ester_volatil` (gas/liquid **partition** factor, van't Hoff), instead of the retired
  `E_a_ester_volatil`. Stripping T-sensitivity = `E_a_uptake + dH_ester_volatil` ≈ 100
  kJ/mol — **the same physical value in both media**.
- `dH_ester_volatil` = **45 000 J/mol**, *sourced*: ethyl-acetate Henry's-law solubility
  constant temperature dependence `d(ln kH)/d(1/T)` ≈ 5300–5700 K (NIST WebBook / Sander
  compilation, doi:10.5194/acp-15-4399-2015) ⇒ dissolution enthalpy ≈ −46 kJ/mol ⇒ the
  gas/liquid partition rises with T with effective enthalpy ≈ +45 kJ/mol, **Q10 ≈ 1.8** — a
  *physical* volatility value, not the fudged Q10 ≈ 5.6 D-20 needed. Identical in both YAMLs.
- `E_a_esters` is now **sourced per medium** (was a generic 80k both media under D-19):
  **beer 200 000 J/mol** (de Andrés-Toro steep ester-growth coupling, transferred as an
  ordering to our flux-coupled term) and **wine 55 100 J/mol** (= `E_a_uptake`). The wine
  value rests on a clean **mapping**: run-integrated synthesis scales as
  `arrh(E_a_esters)/arrh(E_a_uptake)` (the bare-flux integral to dryness is fixed by total
  sugar), so it is **T-independent exactly when `E_a_esters = E_a_uptake`** — the Arrhenius
  representation of Mouret's *flat/weak* wine ester production. Not a coincidence; the
  condition for flat integrated production.

**Why this is strictly more faithful (the point of choosing it).** Both directions now
emerge from **physical + sourced** parameters, with no compensating constant:
- **Wine:** synthesis flat (`E_a_esters = E_a_uptake`) + steeper physical stripping (~100k)
  ⇒ liquid esters **fall** with T (Rollero evaporation inversion), total production stays
  **flat**, and the stripped fraction (`esters_gas`) **rises** with T.
- **Beer:** synthesis steep (200k ≫ 100k stripping) ⇒ liquid esters **rise** with T
  (de Andrés-Toro warm-ale character).

D-20 additionally left wine *total* production rising with T (contra Rollero); D-21 fixes
that too — the `E_a_esters = E_a_uptake` mapping makes it exactly flat.

**Architecture (no new contract).** The only modifier on uptake is
`ArrheniusTemperature.for_uptake` (the Luong wall is unwired; `EthanolInactivation` is a
separate Process on `X`), so the gas flow is reproducible as `bare_flux ·
arrhenius(E_a_uptake)`; `EsterVolatilization` applies that factor itself and folds
`q_sugar_max·co2_yield·scale·(gas-volume/Henry-prefactor)` into `k_ester_volatil`. No
two-pass / derivative-passing contract was needed.

**Documented simplification.** The full Morakul (2011) partition is also *ethanol-dependent*
(`ln k_i = F1 + F2·E − (F3 + F4·E)·R·(1000/T − 1000/T_ref)`); we keep only the dominant
temperature (van't Hoff) lever via `dH_ester_volatil` and omit the ethanol terms (the `F`
coefficients are not openly available). All four volatilization/synthesis-E_a params stay
**speculative** in magnitude; the *orderings and the flat-production mapping* are sourced.

**Empirical results (verified, carbon closing to machine precision every run).** Wine liquid
esters **73 → 61 → 50 mg/L** (14/20/25 °C, fall), gas **41 → 53 → 64** (rise), total **flat
~114**; fusels **45 → 51 → 56** (rise). Beer liquid esters **22 → 72 → 181 mg/L** (rise);
fusels **37 → 41 → 46**. §2.2 trio unmoved (all at 20 °C where every Arrhenius factor = 1).
The directional benchmark `test_lower_temperature_is_slower_but_cleaner` passes per medium
on liquid pools; the unit guard `test_integrated_wine_aroma_temperature_directions` now also
asserts `esters_gas` **rises** with T. **222 tests green**, ruff + format + mypy clean.

## D-22 — SO₂ speciation: the pH-coupled molecular fraction, as a readout-only derived function

**Status: settled (built 2026-06-30).** The first consumer of the D-18 pH keystone, and
the payoff its "dose SO₂ → speciation falls out of the current pH" promise was written
against. Scope mirrors D-18's own deliverable boundary: **a derived pure-function readout,
no RHS consumer** (the antimicrobial suppression of MLF/spoilage growth wires in with those
organisms, exactly as pH had no consumer in D-18).

**What SO₂ does in wine, and what beat 1 covers.** Free SO₂ partitions by pH into
**molecular** SO₂·H₂O (the antimicrobial species), **bisulfite** HSO₃⁻ (dominant at wine
pH), and negligible **sulfite** SO₃²⁻; the molecular fraction is
`1/(1 + 10^(pH − pKa₁))` with pKa₁ ≈ 1.81, so it falls ~3× per 0.5 pH unit. Beat 1 builds
exactly this **free-SO₂ speciation readout**. The **free/bound split** (SO₂ reversibly
binds acetaldehyde and other carbonyls) is **deferred** — acetaldehyde is an unbuilt §3.2
byproduct — which is why the scenario input is **free SO₂ (mg/L)**, the variable winemakers
actually measure and target, not a total dose (treating a total addition as all-free would
overestimate molecular SO₂; framing the input as free makes the deferral honest, not a hole).

**The decision: readout-only — SO₂ is a state slot but NOT in the charge balance.** The fork
was whether sulfurous acid joins the proton/charge balance (so its bisulfite charge nudges
pH, and dosing SO₂ *acidifies* emergently) or pH is solved from the organic acids and free
SO₂ partitioned at that pH as a pure readout. **The D-18 inverse anchoring collapses the
fork at t=0:** `solve_cation_charge` back-solves the strong cation to reproduce `initial_ph`
*exactly*, so if SO₂ were in the balance at pitch, the fitted cation would simply absorb its
~0.6–0.8 meq/L of bisulfite charge and pH(t=0) would *still* be `initial_ph`. So the molecular
SO₂ number at t=0 — the only place fidelity is anchored — is **identical** in both designs;
the in-balance gain is ~zero where measured and second-order over the run (on top of an
already directional-only pH drift), while its cost (refactoring the freshly-landed D-18
signatures `charge_residual`/`solve_ph`/`solve_cation_charge` + the compile anchoring block)
is real. **Readout-only wins**, and it is still fully compositional — the forward coupling
D-18 promised is delivered, nothing is scripted. SO₂'s back-reaction on pH is a **scoped
caveat** (like carbonic in D-18, but smaller relative to its own effect): the reverse coupling
only becomes *visible* under a mid-ferment SO₂ *addition event* (unbuilt), and when wanted it
should be added by **generalizing `Byp`'s separate-arg into an `extra_acids: Mapping[str,float]`**
of non-carbon charge-active species (Byp + SO₂ both entries), not a 5th positional arg.

**What landed.**
- **`so2_free` state slot** on `wine_schema` only (g/L of SO₂-equivalent; `default=0.0`,
  inert — no Process touches it, so it is constant exactly like the D-18 acids). Beer is
  untouched (its acid/SO₂ system is deferred with its pH). Dosed via the optional scenario
  input `so2_free_mgl` (mg/L → g/L at compile); it does **not** enter the cation back-solve.
- **`acidbase.molecular_so2(y, schema, params)`** — the headline derived pure function:
  solves pH from the organic acids (`ph_of_state`), then returns `free_SO₂ × neutral_fraction(pH)`.
  Plus `molecular_so2_fraction(ph, pkas)` and a new `neutral_fraction(h, pkas)` (the
  undissociated-species share `h²/D`, the complement of `mean_charge`'s dissociation), and
  the `molecular_so2_series` analysis helper. Free SO₂ is expressed *as SO₂*, so the
  partition is mass-preserving and the readout needs no molar conversion; `units.gpl_to_mgl`
  reports the conventional mg/L.
- **`pKa_sulfurous_1` = 1.81, `pKa_sulfurous_2` = 7.20** in `acidbase.yaml`, sourced
  (Usseglio-Tomasset & Bosia 1984, carried in Boulton and Ribéreau-Gayon; CRC for pKa₂),
  tier **plausible**. **Deliberately kept out of `PKA_PARAM_NAMES`** (the pH-solver acid set):
  `build_pka_map`/`charge_residual` never see them — the structural guarantee that SO₂ is
  readout-only.
- **`M_SO2` = 64.06** chemistry constant (registered with **0 carbon atoms**, so
  `carbon_mass_fraction("sulfur_dioxide") = 0.0` and the slot is carbon-inert in every sum).

**Two caveats, both load-bearing, both scoped:**
- **Excluded from titratable acidity.** OIV TA explicitly excludes sulfurous (and carbonic)
  acid; readout-only gives this for free since SO₂ is not in `ACID_STATE`. This is *not*
  cosmetic — pKa₂ ≈ 7.2 means sulfite *is* partly formed at the pH-8.2 titration endpoint, so
  an SO₂-in-`ACID_STATE` design would have wrongly inflated TA.
- **Back-reaction on pH omitted** (the readout-only choice above); justified by the anchoring
  argument, additive to restore later.

**Tier = `plausible`, computed explicitly.** `acidbase.molecular_so2_tier` combines **both**
pKa sets — the pH-solver pKas (the readout solves pH) *and* the sulfurous pKas — floored at
`PLAUSIBLE`. SO₂ speciation is never `VALIDATED`: apparent constants applied to wine are
extrapolation, and the acceptance gate checks our implementation against Henderson-Hasselbalch
(the equation itself), a self-consistency check, not an independent dataset.

**Acceptance (met).** The molecular fraction lands on the textbook curve — **6.07 % / 2.00 %
/ 0.64 %** at pH 3.0 / 3.5 / 4.0 — and falls ~3× per 0.5 pH unit. The free SO₂ needed for the
**0.8 mg/L molecular** microbial-stability target reproduces the canonical winemaking table
(**~13 / 32 / 40 / 50 / 79 / 125 mg/L** at pH 3.0 / 3.4 / 3.5 / 3.6 / 3.8 / 4.0). Prime
directive #3 is pinned by an **isolability** test: on a shared time grid, dosing 60 mg/L SO₂
leaves every other state column byte-identical, the pH series identical, and carbon closing —
SO₂ is genuinely inert and outside both the charge balance and the carbon ledger. The series
also shows the molecular fraction **rising** late as the emergent `Byp` pH drift pulls pH
down — unscripted, the D-18 coupling working through SO₂. **249 tests green** (236 → +12 SO₂
+1 chemistry), ruff + format + mypy clean. This unblocks **MLF** (whose *O. oeni* growth is
SO₂-sensitive — the first RHS consumer of `molecular_so2`).

## D-23 — MLF v1 is conversion-only; the amino-acid ledger is a separate yeast/AF beat

**Status: scoped 2026-06-30; v1 IMPLEMENTED 2026-07-01 (see "Resolution" below).** Records
the design call for the beat — *Oenococcus oeni* malolactic fermentation — the empirical
evidence that settles it, and (Resolution) the open-knob choices made when v1 landed.

**The fork.** MLF converts L-malic acid (C4, diprotic) to L-lactic acid (C3, monoprotic) + CO₂,
mole-for-mole, deacidifying the wine (pH up ~0.1–0.3). The question was whether v1 should model
the *bacterium's growth* — and if so, where its biomass carbon comes from. Three paths surfaced:
(B2) **conversion-only** — run the malate→lactate flux with no bacterial biomass dynamics;
(B1-malate) growth funded from malate carbon; (B1-aa) growth funded from amino acids — the
biologically-right source, which requires making nitrogen carry carbon, a change to the protected
validated core.

**The amino-acid carbon problem, and the toggle that defuses it.** Path B1-aa is the honest one
— *O. oeni* builds biomass mostly from amino acids/peptides, not hexose — but `N` (YAN) is
deliberately carbon-free in `total_carbon` (D-19), so making amino acids a carbon source is a
*non-isolable* change to the core carbon ledger **and** the growth kinetic, violating prime
directive #3. The owner's proposal — a **toggleable amino-acid ledger** (a `default=0` pool that,
when populated, contributes to *both* the carbon and nitrogen ledgers) — restores isolability:
when the pool is empty the carbon term is additively zero and the core is byte-for-byte. The
advisor refined the *mechanism*: rather than a two-mode fork inside `GrowthNitrogenLimited` (a
permanent branch through the core's hottest kinetic, with a float-identical collapse you must
*prove*), implement it as a **separate isolable Process** — a pure *swap* that, for the
amino-acid-funded fraction of biomass, refunds sugar by the displaced biomass carbon, refunds the
ammonium `N` pool by the displaced biomass nitrogen, and debits the amino-acid pool by one
amino-acid mass carrying exactly that C and N. The swap is carbon-neutral **and** nitrogen-neutral
by construction, leaves growth (and the Coleman reconstruction) byte-for-byte untouched, and
contributes zero when the pool is empty — isolability is *structural*, not a tested coincidence.
Its one new input is the amino-acid pool's C:N ratio (a sourced, speculative `Parameter`).

**Why it is nonetheless a *separate* beat, not part of MLF — settled by running the model.** The
decisive question is whether the amino-acid pool has anything in it *at the MLF pitch point*. It
does not. A standard 24 Brix wine AF (the §2.2 Coleman anchor, 20 °C) was integrated to
completion and the lumped `N` trajectory inspected:

| Must | N first < 1 mg/L | N at dryness (pitch point) | Days to dryness |
|------|------------------|----------------------------|-----------------|
| 80 mg/L (Coleman low-N) | day 1.29 | ≈ 0 | 8.33 d |
| 300 mg/L (richly dosed) | day 1.33 | ≈ 0 | 5.17 d |

`N` is driven to the solver floor (~0) within ~1.3 days of pitch and sits there for the entire
post-AF period — *regardless of dose*. So at the MLF pitch (dryness, day 5–8) there is no
nitrogen, and the future amino-acid pool would be in exactly the same place (the same uptake that
drains `N` drains it). **MLF-growth is therefore structurally blocked until something replenishes
the pool post-AF** — an autolytic-peptide flux (yeast death → peptides → amino-acid pool),
unbuilt. The toggleable aa-ledger improves *primary-fermentation* (yeast) carbon honesty and is
the natural home to later re-route the D-19 fusel Ehrlich carbon off its sugar stand-in — but it
does not feed the bacteria. Hence: **MLF v1 = conversion-only; the amino-acid ledger is its own
yeast/AF beat; MLF-growth is a still-later composition of the two plus autolysis.**

**A model gap surfaced by the same run (flagged, not fixed).** The model drives even a 300 mg/L
must to *zero* nitrogen within ~1.3 days — it has no satiation cap, no luxury-uptake ceiling, no
residual-N floor. Real musts finish with 50–150 mg/L residual YAN plus an unusable **proline**
tail (yeast cannot assimilate proline anaerobically). So the model *overstates* nitrogen
exhaustion. This matters for the aa-ledger beat: doing it *honestly* means also modeling that
yeast stop assimilating when sated, otherwise the post-AF amino-acid residue is artificially
empty. More scope → more reason it is a careful separate beat, not a rider on MLF.

**MLF v1 scope (what the implementation session builds).**
- **Carbon closes on the existing ledger** — malic (C4) → lactic (C3) + CO₂ (C1) are already
  weighted in `total_carbon` (`chemistry.py`, anticipated since D-18); no new conservation code.
- **`X_mlf` as a dosed-but-inert catalyst slot** on `wine_schema` (`default=0.0`, isolable),
  dosed via a new scenario input `mlf_pitch_gpl`. In v1 *no Process grows or kills it* — it is a
  constant bacterial concentration scaling the conversion rate, so the later growth beat is a
  clean extension (add a growth Process touching `X_mlf`), not a refactor.
- **`MalolacticConversion` Process** — touches `malic`/`lactic`/`CO2`, reads `X_mlf`, pH
  (`ph_of_state`), molecular SO₂ (`molecular_so2`), ethanol `E`, and `T`. Flux is substrate-limited
  in malate, scaled by `X_mlf`, and gated by inhibition factors: low pH, high ethanol,
  **molecular SO₂** (the first RHS consumer of D-22), and a temperature optimum. Tier
  **speculative**.
- **Acceptance gate** — the existing hand-built `test_headline_malic_to_lactic_raises_ph`
  ΔpH ∈ [0.1, 0.3] (lands 0.225) becomes *emergent* from the Process on a malic-rich must.
- **Scope boundary** — runtime has no event mechanism, so v1 models **co-inoculation** MLF
  (bacteria present from t=0). **Sequential / post-AF MLF** (pitch at day N) needs the
  event-driven loop (deferred, see `runtime/integrate.py` docstring). Open knobs for the
  implementation session: the exact inhibition functional forms and their sourcing; whether
  `X_mlf` is explicit or folded into the rate constant.

**Resolution (v1 landed 2026-07-01).** `core/kinetics/malolactic.py`
(`MalolacticConversion`), `X_mlf` slot on `wine_schema`, `mlf_pitch_gpl` scenario input,
the *O. oeni* parameter block in `wine_generic.yaml`, and `tests/test_malolactic.py` (13
tests). 262 green, ruff + mypy clean, §2.2 trio unchanged. The molar turnover is

    r = k_mlf · X_mlf · [malate]/(K_mlf+[malate]) · g_pH · g_EtOH · g_SO₂ · γ(T)   [mol/L/h]

with `d(malic)=−r·M_malic`, `d(lactic)=+r·M_lactic`, `d(CO2)=+r·M_CO2`. Carbon *and* mass
close on the existing ledger (4 C = 3 C + 1 C; 134.087 = 90.078 + 44.009 g/mol, a clean
decarboxylation, no water term), so no new conservation code — verified at the RHS level
(weighted carbon rate ≈ 0) and over a full dosed run.

*The open knobs D-23 left open — chosen, all speculative-tier:*
- **`X_mlf` explicit** (scales the rate), not folded into `k_mlf` — keeps the later
  growth beat a clean add-a-Process extension.
- **Temperature = a cardinal-temperature optimum** (Rosso et al. 1993 CTMI,
  `cardinal_temperature_factor`; cardinals 8/23/37 °C), *not* a monotone Arrhenius — MLF
  genuinely declines in the warm, which Arrhenius cannot represent (the load-bearing reason
  D-23 named "a temperature optimum"). Peak 1 at `T_opt`, 0 outside `[T_min, T_max]`.
- **pH gate** = smooth logistic `1/(1+10^(pH_half−pH))` (midpoint pH 3.0): rises with pH, so
  malate→lactate deacidification is *self-reinforcing* (pH↑ ⇒ rate↑), bounded by 1 and
  self-limited by malate depletion — the emergent coupling the D-18 keystone exists for.
- **ethanol gate** = the Luong wall `max(0, 1−E/E_max)^n` reused from `EthanolInhibition`
  (`ethanol_tolerance_mlf` 110 g/L ≈ 14 % ABV, *below* the yeast's 142).
- **molecular-SO₂ gate** = `exp(−[SO₂]_molecular/s)`, partitioned at the *solved* pH — the
  first RHS consumer of the D-22 readout. Dosing ~80 mg/L free SO₂ arrests MLF (verified).

*Isolability (prime directive #3), two layers:* (a) **value** — the Process returns a zero
contribution *before* the per-RHS pH `brentq` whenever `X_mlf ≤ 0` or malate is gone, so an
undosed run is byte-for-byte the validated core and pays no solve; (b) **tier** — the
compile seam **disables** the Process when `mlf_pitch_gpl ≤ 0`, because `ProcessSet.tier_of`
counts *enabled* (not nonzero) processes, so an always-on-but-zero MLF would drag the inert
`malic`/`lactic` slots from VALIDATED to speculative on every undosed wine run. (`CO2` is
already speculative via the uptake Process, so it is unaffected either way.) When pitched,
`malic`/`lactic`/`CO2` correctly become speculative.

**Emergent finding — the ethanol "race-or-stall" (a genuine model behavior, flagged).** A
24-Brix must reaches ~135 g/L ethanol but `ethanol_tolerance_mlf` is 110, so the ethanol
gate **arrests MLF once AF ethanol crosses ~110 g/L (~day 4 at 20 °C)**. MLF must therefore
**complete in that early low-ethanol window or stall permanently** (ethanol never falls) —
which is *exactly why co-inoculation is used in practice*, and why in this model
co-inoculation is the only viable mode: post-AF (sequential) MLF is **doubly blocked** — no
event loop to pitch at day N *and* ethanol already past tolerance — reinforcing D-23's
co-inoculation scope. `k_mlf` (default 1.5e-2, speculative/order-of-magnitude) is tuned so a
realistic pitch (test uses 0.2 g/L) converts a malic-rich must to ~complete within that
window. Two honest caveats: (i) the 110 g/L wall is a speculative simplification — real
high-alcohol MLF strains tolerate ~15–16 % ABV; (ii) the **headline test is coupled to AF
timing** — a future change that speeds AF shrinks the MLF window, but the test (ΔpH ≥ 0.1)
would catch the regression, so the coupling is safe-but-explicit.

**Acceptance — added, not replaced (D-23 "becomes emergent").** The new headline
`test_headline_mlf_raises_ph_emergently` measures the **no-MLF control difference**
`pH_final(dosed) − pH_final(off)` = **0.1813** ∈ [0.1, 0.3]: robust because MLF touches only
`malic`/`lactic`/`CO2` and pH reads neither `CO2` (carbonic omitted, coupling #1) nor any AF
variable, so the two runs are byte-identical in X/S/E/N/Byp/cation and the gap is *purely*
the malic→lactic swap at the same final Byp. The original algebraic
`test_acidbase.test_headline_malic_to_lactic_raises_ph` (0.225) is **retained** — the two
prove different things (the solver responds to acid dynamics vs the Process *produces* those
dynamics).

**Minor (noted, not fixed).** The *O. oeni* parameters live in `wine_generic.yaml` (the
ester/fusel aroma set the precedent for non-yeast mechanisms there, and the wine compile
loads exactly that file so beer never sees them), but they are bacterium properties, not
yeast-strain ones — so a *second* wine-strain file would duplicate them, the same re-homing
caveat already flagged for `must_fermentable_fraction`.

## D-24 — Stochastic ensemble wrapper: Monte-Carlo over provenance bands, in the runtime

**Status: IMPLEMENTED 2026-07-01** (`runtime/ensemble.py`, `tests/test_ensemble.py`, 274 green).
The last big Milestone-2 item that carried no new physics — the parallel, physics-free beat
(`milestone-2-tasks.md`) the handoff §1.6 calls for: *"realism and replicate variation come
from a runtime layer that samples parameters within their provenance-declared uncertainty and
runs ensembles."* Every `Parameter` has always carried an `Uncertainty` band; until now nothing
at runtime read it.

**The seam.** `simulate_ensemble(process_set, parameters, y0, t_span, …)` takes the full
`ParameterSet` (it needs the bands) — the natural distinction from `simulate`, which takes
resolved floats. It draws `n_members` samples, integrates each with `simulate` on a shared
`t_eval` grid, and returns an `Ensemble`: the deterministic **nominal** run, the surviving
**members** `(n_succeeded, n_vars, n_times)`, each member's sampled param map, and the derived
`tier_map`. Randomness lives **only here**, behind an explicit `seed` — the core stays pure and
a single unsampled run stays byte-for-byte reproducible (the architecture rule + §1.6 split).

**Choices made (all revisited with the advisor):**

1. **Distribution = triangular `(low, mode=value, high)`**, `uniform` pluggable. "Bounds plus a
   most-likely value" is the textbook triangular case, and `value` *is* the sourced, benchmarked
   most-likely estimate — uniform would throw that away (extremes as likely as the best estimate).
   The reported band uses **outer percentiles (P5/P95 default)**, which keeps the full bracket
   visible and de-sensitises the result to the shape choice. Zero-width bands (`high ≤ low`) pin to
   `value` and consume no randomness.
2. **Plain Monte Carlo** by default, the method §1.6 names. Latin-hypercube / Sobol give better
   tail coverage per member; added as opt-in `sampler=` strategies in **D-25** (MC stays default).
3. **Sample only what the *active* Process set `reads`** (union of `Process.reads` +
   `RateModifier.reads`), intersected with the loaded params. Sampling anything else is a no-op on
   the trajectory and only dilutes the member count, so the spread means "sensitivity of *this*
   scenario". `only` overrides the set; `exclude` removes names from it (the pinning escape hatch).
   A neat consequence: on an undosed (MLF-off) wine run the pKa set is not read, so it is not
   sampled — the D-18 initial-pH anchor (back-solved at compile from nominal pKa) is untouched.
   When MLF *is* pitched the pKa set enters scope and the anchor holds only at nominal; that drift
   is *honest* (pKa uncertainty → uncertainty in the implied cation charge), and `exclude` pins it
   for a caller who wants the anchor preserved.
4. **Parameter uncertainty only** — scenario/initial-condition uncertainty (Brix, YAN) is a
   separate axis; `y0` is held fixed.
5. **Nominal ≠ median, and both are reported.** The median of nonlinear trajectories is not the
   trajectory of median parameters; the nominal is the deterministic reference, the median+band is
   the uncertainty summary.

**Independence caveat — checked against the actual bands, not hand-waved.** Parameters are sampled
independently, which ignores cross-parameter constraints. The two live groups were enumerated and
checked against their real `Uncertainty` bands (the advisor's decisive point: overlap decides
whether the caveat is vacuous, immaterial, or real):

- **Realised-yield partition — vacuous.** The uptake Process does *not* read `Y_ethanol_sugar`;
  ethanol/CO₂ use the theoretical Gay-Lussac split *scaled down*, and glycerol/byproduct carbon is
  **carved from** that same flux (`scale = 1 − diverted_c/c(species)`), with a hard `ValueError`
  guard if `scale < 0`. At band maxima `diverted_c ≈ 0.027` vs `c(glucose) ≈ 0.40` → `scale ≈ 0.93`;
  super-theoretical yield is structurally unreachable, the guard is a backstop (and a member that
  tripped it would be *counted as failed*, not silently dropped).
- **Load-bearing `E_a > E_a_uptake` byproduct ordering — immaterial.** Wine `E_a_esters` [40k,70k]
  fully overlaps `E_a_uptake` [47k,63k], but the wine ester T-direction is *intentionally null*
  (nominal `E_a_esters == E_a_uptake`, Mouret-flat, D-21) — scrambling it corrupts no demonstrated
  result. `E_a_fusels` [60k,250k] overlaps uptake only in [60k,63k], a tail-tail sliver where the
  triangular joint density ≈ 0. Beer `E_a_esters` [120k,265k] has *no* overlap → safe. Nominal
  orderings hold for the overwhelming majority; a stray inverted member is honest parameter
  uncertainty within a *speculative* band, and `exclude` pins the group for a strict ensemble.

**No silent truncation.** A sampled param set can make a member fail — `solve_ivp` returns
`success=False`, or the RHS *raises* (the uptake guard). Both are caught, recorded in `failures`,
and counted; the RNG advances one sample per member so reproducibility (including *which* members
fail) holds. Past `max_failure_fraction` (default 0.5) the driver **raises** rather than return a
survivorship-biased spread from the lucky survivors.

**Per-member conservation is the crown-jewel invariant.** `Ensemble.member_trajectory(i)`
reconstructs any member as a `Trajectory` so the deterministic harness (`assert_conserved`, …)
audits it. Carbon closes for *every* sampled member — but the check must use that member's **own**
accounting constants (e.g. its sampled `biomass_C_fraction`, which the growth Process draws sugar
carbon against), which is exactly why `member_params[i]` is stored; auditing with the nominal
constant reads genuine closure as drift.

## D-25 — Ensemble follow-ups: spread attribution, LHS/Sobol, per-member nitrogen

**Status: IMPLEMENTED 2026-07-01** (288 green). Three natural extensions of the D-24 ensemble —
*not gaps in it*, but the questions it makes askable. Built in the advisor-recommended order
(cheap probe first, refactor last), each committed separately.

1. **Per-member nitrogen conservation** (`tests/test_ensemble.py`). The D-24 crown-jewel
   (per-member carbon closure) extended to the nitrogen ledger. **Probed before trusting:** N
   closes to ~1e-12 across every member using that member's **own** sampled `biomass_N_fraction`
   (the growth Process draws N against it) — expected, since the aa-ledger is deferred (D-23) and
   fusels route *carbon*, not N, from sugar, so biomass is the only N sink. A failure here would
   have been a real N-leak finding, not a test to force green.

2. **Spread attribution by parameter and tier** (`analysis.attribute_spread`,
   `tests/test_attribution.py`). A first-order variance decomposition computed **post-hoc from one
   ensemble's stored `member_params`** — no extra integrations (OAT would need N extra ensembles and
   is a known-poor sensitivity method). Standardized-regression coefficients (SRC): because D-24
   samples parameters *independently*, the SRC² are near-orthogonal and ≈ sum to the regression R²,
   giving a genuine variance split; shares roll up by parameter `Tier`. **R² < 1 is expected** (the
   model is nonlinear — Monod/logistic/Arrhenius), so `1 − R²` is reported explicitly as the
   `unexplained` interaction/nonlinearity bucket — the budget never reads as "everything explained".
   `method="srrc"` rank-transforms first (robust fallback for monotone-but-curved responses). Needs
   n≳50–100 members for a stable fit (underdetermined fits raise). Lives one layer up in
   `analysis.py` (top-level observable over a runtime `Ensemble`), *not* core — attribution needs
   parameter tiers, passed in via `ParameterSet.tier_map()` (the Ensemble's `tier_map` is per state
   *variable*). On the wine ferment: ethanol spread is driven by `k_prime_d` (inactivation) and
   `q_sugar_max`; SRC R²≈0.6, SRRC≈0.72 surfacing the competing `Y_glycerol_sugar` sink.

3. **LHS / Sobol samplers** (`simulate_ensemble(sampler=…)`). `"mc"` stays the default and is
   **byte-identical** to before (same seeded PRNG sequence); `"lhs"` and `"sobol"` draw a stratified
   unit hypercube via `scipy.stats.qmc` then map it through each parameter's inverse CDF (triangular
   via `scipy.stats.triang`, `c=(value−low)/(high−low)`; or uniform). At a fixed member budget the
   estimator is ~8× more stable seed-to-seed than i.i.d. MC on the toy, with the **center unshifted**
   (the point: tighter tails, not a moved mean). Design constraints, all from the advisor:
   `only`/`exclude` scoping and the failed-member/survivorship accounting are **sampler-agnostic**;
   only *varying* parameters take a hypercube dimension (a pinned zero-width band stays at nominal —
   giving it a column wastes a dimension, unbalances Sobol, and divides `c` by zero); **Sobol requires
   a power-of-two `n_members`** and raises otherwise (no silent unbalanced sequence — the project's
   loud-failure ethos). Samples are drawn up front, so seed reproducibility holds for every sampler.

## D-26 — Diacetyl (vicinal diketones): the mechanistic 3-pool "diacetyl rest"

**Status: IMPLEMENTED 2026-07-01** (320 green). The flagship of the remaining §3.2
byproducts (diacetyl / acetaldehyde / H₂S). Diacetyl (2,3-butanedione, a buttery off-note)
is *the* defining lager-quality parameter, and unlike the monotone-accumulate ester/fusel
pools it is **produced then reabsorbed** — a non-monotonic time course (the "diacetyl rest").
Built as three commits (one Process each), one beat.

**The forks the owner decided (surfaced before building, per the "discuss disagreements"
rule).** Two were genuinely the owner's call:

1. **Sequencing:** diacetyl → acetaldehyde → H₂S, one Process per commit (owner chose the
   incremental order over one big beat). Diacetyl first: it is the flagship *and* the
   cleanest instance of the new produce-then-reabsorb shape, so it establishes the reusable
   kinetics before acetaldehyde (the thorniest — it sits on the main ethanol pathway).
2. **Carbon accounting — "something closer to reality"** than either offered default. The
   two easy options were (A) route production carbon from sugar and *return* reabsorbed
   carbon to sugar (a "returns-to-sugar" bookkeeping stand-in), or (B) a carbon-unaccounted
   trace pool outside `total_carbon`. The owner rejected both and asked for fidelity. The
   answer: **track the real downstream product.** The true VDK pathway is

   ```
   sugar → α-acetolactate → diacetyl + CO₂ → 2,3-butanediol
     (draw from S)   C5      (C5→C4+C1)  C4    (C4→C4)   flavourless
   ```

   Every step closes carbon on the *existing weighted ledger*: the α-acetolactate draw from
   sugar is the D-19 option-a1 routing; the decarboxylation `C5 → C4 + CO₂` is carbon-closing
   exactly like malolactic `malic → lactic + CO₂` (D-23); the reduction `C4 → C4` is a
   mole-for-mole transfer to a real tracked pool, like `esters → esters_gas` (D-20). No
   stand-in for the reabsorbed carbon, no vanished mass. `total_carbon` closes to machine
   precision through the whole produce-then-reabsorb course. (`total_mass` gains a small gap:
   the oxidative decarb consumes untracked O₂ and the reduction untracked NAD(P)H — carbon is
   the invariant, as for beer's hydrolysis water, D-8.) The α-acetolactate-from-sugar draw is
   *better* grounded than the ester/fusel stand-ins — α-acetolactate genuinely derives from
   pyruvate.

**The fidelity target (the second owner fork): C-full, not C-minimal.** The discriminator
put to the owner was: *must the model reproduce "crash/package too early ⇒ diacetyl rises"
and "a warm rest clears it faster"?* Yes ⇒ the **3-pool** model with the α-acetolactate
**reservoir**, not a 2-pool (diacetyl produced flux-linked, reduced by live yeast). The
reservoir is **load-bearing, not cosmetic**: in the 2-pool model diacetyl generation dies
with the sugar, so it can neither strand a *rising* diacetyl after a crash nor make the rest
temperature-critical. The advisor's earlier "defer the α-acetolactate lag for v1" was
explicitly reversed here for exactly this reason.

**Why the rest emerges (the three Processes, `core/kinetics/vicinal_diketones.py`):**

- **`AcetolactateExcretion`** fills the reservoir from the fermentative flux (shared
  `K_sugar_uptake`), so it stops at dryness — the reservoir is full at end of primary.
  **Temperature-flat** (a documented v1 simplification: the reservoir *size* is a weak lever;
  the temperature-criticality lives downstream). Draws its C5 carbon out of `S`.
- **`AcetolactateDecarboxylation`** converts reservoir → diacetyl + CO₂ by a **spontaneous,
  non-enzymatic, first-order, strongly temperature-dependent** reaction that is **NOT gated
  on yeast** — so it keeps making diacetyl *after* fermentation, faster when warm. This is
  the **rate-limiting, temperature-critical** step (`E_a_decarb` held high). Sourced ordering
  (Haukeli & Lie 1978; Krogerus 2013 review, doi:10.1002/jib.84 — "higher fermentation
  temperatures increase the conversion rate"); magnitude speculative.
- **`DiacetylReduction`** is **fast, enzymatic, gated on VIABLE `X` (not `X_dead`), with NO
  flux term** — so it clears diacetyl as fast as it forms while live yeast is present, but
  **stops dead** once the yeast is crashed / racked / ethanol-inactivated. The no-flux-term
  is essential: reduction must run during the rest (flux ≈ 0). `E_a_reduction` is held
  **below** `E_a_decarb` so decarb stays rate-limiting.

Together these make the defining behaviour *emerge*. **Verified empirically** (not asserted)
before the acceptance test was written:

| medium | 14/10 °C | 20/18 °C | 28/25 °C |
|---|---|---|---|
| **beer** final diacetyl | 0.195 (stranded, reservoir 4.7) | 0.040 | 0.001 mg/L |
| **wine** final diacetyl | 1.011 (stranded, reservoir 1.1) | 0.179 | 0.001 mg/L |

Warmer ⇒ monotonically cleaner (the headline "warm rest clears it faster"); a warm run shows
**peak-then-fall** (beer 25 °C peaks 0.076 @ day 4 → clears to 0.001); a cold run **strands**
diacetyl at its peak with a large **unconverted α-acetolactate reservoir** the warm run
consumes. The cold cases sit above the ~0.1 mg/L lager flavour threshold (a real off-note),
the warm cases well below.

**Isolability / wiring.** The three Processes live in their own `_VDK_PROCESSES` tuple. Unlike
MLF (a *dosed* organism, disabled at compile when unpitched), diacetyl is **intrinsic yeast
metabolism**, so it is wired into **both** media and runs on every default ferment — like the
ester/fusel byproducts. Turning it on draws only a *trace* of sugar (α-acetolactate peaks
~mg/L, roughly an order of magnitude below the ester draw), so `dX`/`dE`/`dCO₂`/`dN` stay
byte-for-byte until the decarb/reduction move that carbon on; the §2.2 trio is unmoved.
**Tiers:** all three Processes **speculative** (rate magnitudes are order-of-magnitude
estimates; only the `E_a_decarb > E_a_reduction` ordering is sourced), so parameter-tier
propagation (D-1) caps the pool outputs at speculative regardless.

**One honest tier consequence — the D-19 `S` parallel, made explicit (not silent).** The
decarboxylation is always-on, speculative, and the *first* such Process to write the shared
`CO2` slot (uptake aside; esters/fusels touch `S`, MLF is disabled unpitched). So on a default
run the *structural* `tier_of("CO2")` drops **PLAUSIBLE → SPECULATIVE** — exactly as
`tier_of("S")` did when the D-19 byproducts landed. But the **param-aware tier users actually
see was already SPECULATIVE** (the uptake Process reads speculative params — `E_a_uptake`,
realised-yield), so there is **no headline change**, and the drop is *honest*: the `CO2` pool
now genuinely contains a speculative decarb trace (real evolved CO₂ that belongs there —
sequestering it into a side pool to protect the tier would understate CO₂, a worse
dishonesty). Accepted as the correct behaviour, and **pinned by a test** (`test_vicinal_
diketones.py`) so it can never regress silently — the beer CO₂-ratio value stays in-band and
its user-facing tier is unchanged.

**Parameters** live in a new **shared, medium-agnostic** `vicinal_diketones.yaml` (merged at
the compile seam alongside `acidbase.yaml`), because the load-bearing decarboxylation is
*non-enzymatic* — a molecule property, not a beverage property (contrast the *per-medium*
ester `E_a`). Also promoted the shared `draw_carbon_from_sugar` / `fermentative_flux_shape`
helpers out of `byproducts.py` into `core/kinetics/carbon_routing.py` (one source of truth for
both the aroma and VDK Processes; behaviour unchanged).

**Scope (v1) / deferred.** Yeast valine-pathway diacetyl only — **MLF-derived diacetyl**
(*Oenococcus* from citrate, a real coupling now that MLF exists, D-23) is explicitly **out**,
so wine yeast-pathway diacetyl *understates* real wine diacetyl. The α-acetolactate
extracellular decarboxylation's ethanol/pH dependence (Kobayashi et al.) and its
excretion temperature dependence are omitted; acetoin is lumped into the terminal
`butanediol` pool. The acceptance gate demonstrates the rest via **isothermal** comparisons +
the natural end-of-ferment ethanol inactivation (a legitimate proxy — the mechanism, not the
temperature profile, produces the behaviour). A **temperature-ramp** test (cool ferment → warm
finish vs cool → cold hold, which `temperature_schedule` already supports) would demonstrate
the *literal* "warm rest" / "package early" scenarios and is a cheap deferred follow-up.
**Next in the beat (deferred):** acetaldehyde (produce-then-reabsorb on the *main* pathway —
reuses this shape; it is the carbonyl that binds SO₂, unlocking the D-22 free/bound split) —
**LANDED in D-27** — then H₂S (carbon-free, an inverse-low-N gate — the accounting-easiest,
following the SO₂ precedent).

## D-27 — Acetaldehyde: the main-pathway intermediate as a transient ethanol-carbon buffer

**Status: IMPLEMENTED 2026-07-01** (342 green). The second §3.2 aroma beat after diacetyl
(D-26). Acetaldehyde (ethanal, CH₃CHO) is the obligate intermediate on the *main* alcoholic-
fermentation pathway (sugar → … → pyruvate → acetaldehyde → ethanol) — the "green apple"
carbonyl that accumulates to an early peak during vigorous fermentation and is then reduced
to ethanol. Like diacetyl it is **produced then reabsorbed**, so it reuses the D-26 shape
(flux-linked production + viable-`X`-gated, no-flux-term reduction), but with **no middle
reservoir** (acetaldehyde is produced directly, not via a spontaneous-decarb precursor) — two
Processes, one commit.

**The load-bearing fork the owner decided (the advisor caught my error first).** I had
half-settled on the D-26 forward note's preview — *"acetaldehyde's carbon **draw** is an even
stronger **stand-in**"* — i.e. draw carbon from `S`, book it as acetaldehyde, reduce to
ethanol, mirroring the ester/fusel/acetolactate template. The advisor's decisive catch:
**that template does not apply here, because acetaldehyde's product is `E` itself, not a side
pool.** The uptake Process *already* performs the complete lumped sugar → ethanol + CO₂
conversion (which implicitly includes this intermediate). Drawing *fresh* sugar → acetaldehyde
→ *new* ethanol is therefore a **second, parallel** sugar→ethanol pathway — **net-new ethanol
that inflates ABV and realised yield by an amount scaling with pool *turnover*** (cumulative
acetaldehyde *produced*, not its peak). That is a genuine double-count, not the benchmark-
neutral trace the D-19 ester draw is (ester carbon lands in a side pool genuinely removed from
`E`; acetaldehyde carbon returns to `E`). The forward note had applied the side-pool template
before anyone noticed the product is `E`. Per the "specs aren't gospel / discuss disagreements"
rule this was surfaced to the owner as a fork, who chose the **buffer** model:

* Because acetaldehyde and ethanol are **both two-carbon**, the reduction acetaldehyde →
  ethanol is a mole-for-mole C2 → C2 transfer. So `AcetaldehydeProduction` **holds back** a
  transient slice of the ethanol the uptake just made — reclassifying it as the true
  intermediate: `d(acetaldehyde)/dt = +r`, `d(E)/dt = −r·M_eth/M_acet`, with `r =
  k_acetaldehyde · X · S/(K_sugar_uptake + S)`. No fresh sugar, no CO₂.
* `AcetaldehydeReduction` **returns** it: `d(acetaldehyde)/dt = −L`, `d(E)/dt =
  +L·M_eth/M_acet`, `L = k_acet_reduction · X_viable · f(T) · [acetaldehyde]`.

This **de-lumps** the existing pathway rather than duplicating it. It is *more* faithful, not
merely benchmark-safe: acetaldehyde genuinely **is** obligate in-transit ethanol carbon, so
borrowing from `E` asserts exactly the right provenance; a sugar draw would assert a parallel
pathway that does not exist.

**Carbon / benchmark consequences.** `total_carbon` (which now weights `acetaldehyde` at its
C2 fraction) closes to **machine precision** through the whole produce-then-reabsorb course,
touching **neither `S` nor `CO2`**. The `E` **endpoint** reconverges to the buffer-off core to
**relative ~1e-8** (the pool fully reduces back), so the §2.2 ABV / realised-yield / CO₂
benchmarks are preserved to far below any tolerance — verified, all 5 benchmarks unmoved.
Honest caveats made explicit and pinned by tests: (i) the isolability is **derivative-level**
(`dS`/`dCO2`/`dN` are byte-for-byte given the same state) — the *integrated* `S`/`CO2`/`N`
differ by a tiny ~1e-4 relative **second-order path perturbation**, because `E` feeds the
ethanol-inactivation viability brake, so the transient `E` dip nudges viability; (ii) `total_
mass` gains a small gap (the reduction moves untracked NAD(P)H) — carbon is the invariant, as
for the diacetyl reduction (D-26) and beer's hydrolysis water (D-8). One tier consequence (the
exact D-26 `CO2` parallel, pinned): `AcetaldehydeProduction` is the first always-on speculative
Process to *write* `E`, so the **structural** `tier_of("E")` drops PLAUSIBLE → SPECULATIVE, but
the **param-aware** tier users see was already SPECULATIVE (the uptake Process reads speculative
params), so there is no headline change.

**Emergent, verified empirically before the acceptance test (the D-26 checkpoint discipline).**

| medium | acetaldehyde peak | peak day (of run) | final |
|---|---|---|---|
| **wine** 20 °C | 37.5 mg/L | day 2.7 (of 21) | 0.00 mg/L |
| **beer** 18 °C | 38.2 mg/L | day 1.8 (of 14) | 0.00 mg/L |

The early peak *emerges* (production rides the flux and outruns the still-building reductive
capacity, then reduction — gated on viable yeast, no flux term — draws it back down as the
ferment slows), landing in the real range (wine ~30–80, beer peaks ~20–40 mg/L; threshold
~10–25 mg/L green apple). Warmer clears faster/lower (wine peak 55→37→23 mg/L at 14/20/28 °C,
via the Arrhenius on the enzymatic reduction). A crash before clearance **strands**
acetaldehyde (borrowed ethanol carbon un-returned) — the same live-yeast-gating structure as
the diacetyl rest; demonstrated at the unit level (`X = 0` ⇒ reduction 0).

**Isolability / wiring.** Both Processes live in their own `_ACETALDEHYDE_PROCESSES` tuple.
Like esters and the VDK pools (and unlike the *dosed* MLF organism), acetaldehyde is intrinsic
yeast metabolism, so it is wired into **both** media and runs on every default ferment.
Production is held **temperature-flat** (a documented v1 simplification, like the acetolactate
excretion, D-26); the enzymatic reduction carries the Arrhenius factor. Both Processes are
**speculative** (rate magnitudes are order-of-magnitude estimates; only the mechanism —
acetaldehyde is the obligate main-pathway intermediate reduced to ethanol by ADH — is
textbook, Boulton et al.; Ribéreau-Gayon et al.).

**Parameters** live in a new **shared, medium-agnostic** `acetaldehyde.yaml` (merged at the
compile seam alongside `acidbase.yaml`/`vicinal_diketones.yaml`), because acetaldehyde is
main-pathway yeast metabolism — a property of the pathway, not the beverage.

**Scope (v1) / deferred.** The acetaldehyde metabolite only. Acetaldehyde is the principal
SO₂-binder, so building it as real state **unlocks the deferred free/bound-SO₂ split** (D-22) —
but that is a separate **readout** commit (it only needs this state to exist, and carries its
own fork: does the dosed `so2_free` slot get reinterpreted as *total*, breaking D-22's
`molecular_so2`, or is a separate total/bound accounting added?), kept out of this beat per the
owner's one-Process-per-commit rhythm. **Next in the beat:** the SO₂ free/bound binding readout,
then H₂S (carbon-free, inverse-low-N gate — the accounting-easiest, following the SO₂ precedent).

## D-28 — SO₂ free/bound split: total conserved, free/bound/molecular derived at the solved pH

**Status: IMPLEMENTED 2026-07-01** (349 green). The readout the D-27 forward note anticipated,
unlocked now that acetaldehyde is real state. Acetaldehyde is the principal SO₂ binder in wine:
bisulfite HSO₃⁻ reacts with the carbonyl to a stable hydroxysulphonate adduct, so a share of
dosed SO₂ is **bound** (not antimicrobial, not analytically "free"). D-22 deferred this because
acetaldehyde was unbuilt and framed the dosed slot as *free* SO₂ to keep that deferral honest.

**The fork (D-27-flagged), decided by the owner: reinterpret the slot as TOTAL, derive free/bound.**
Two options were surfaced (per "discuss disagreements"): (1) rename `so2_free`→`so2_total`
(conserved, inert) and derive `bound = f(total, acetaldehyde, pH)`, `free = total − bound`,
`molecular = free × neutral_fraction(pH)`; or (2) keep `so2_free` pinned and add `bound`/`total`
additively. **Option 1 chosen — the decisive reason is conservation:** option 2 is non-conserving
(with free pinned and `bound = f(free, acetaldehyde)`, `total = free + bound` *grows as
acetaldehyde rises with no SO₂ added* — incoherent for a single dose, and it flattens molecular
instead of dipping, killing the payoff). Option 1 gives the real must chemistry — "added SO₂ gets
used up, then released": the early acetaldehyde peak sequesters SO₂ → free/molecular crash →
recover as acetaldehyde is reduced (D-27). At acetaldehyde = 0 the split collapses to D-22 exactly
(`free == total`), so the input-semantics change is invisible at the dosing moment (regression
anchor pinned; the D-22 6.07/2.00/0.64 % curve and the free-for-0.8-molecular table survive).

**The binding equilibrium (`acidbase.bound_so2_molar`, pure algebra).** Referenced to **bisulfite**
(the reactive nucleophile): `K = [free acetaldehyde]·[HSO₃⁻] / [adduct]` with `[HSO₃⁻] = free_SO₂ ·
bisulfite_fraction(pH)`, so pH enters mechanistically (new `acidbase.bisulfite_fraction`, the HA⁻
share `Ka₁·h/D`). With `A` = total acetaldehyde, `C` = total SO₂, `β` = bisulfite fraction, the 1:1
adduct `x` solves `(A−x)(C−x)·β − K·x = 0` — a quadratic whose *smaller* root is physical (clamped
to `[0, min(A,C)]`). pH is solved from the organic acids **first** (SO₂ still out of the charge
balance, D-22), so there is no circularity: `β` uses the organic-acid pH. Readouts:
`speciate_so2` (one pH solve → `So2Speciation(total, bound, free, molecular, …)`), thin scalar
wrappers `bound_so2`/`free_so2`/`molecular_so2`, and `molecular_so2_at_ph` for in-loop reuse.

**The one live consumer: the MLF antimicrobial gate.** MLF suppression is by *molecular* SO₂ —
the undissociated share of **free** SO₂ — so the gate (D-23) now reads the *derived* free-molecular
pool via `molecular_so2_at_ph` instead of the raw slot (bound SO₂ is not antimicrobial). This is a
correct consequence of the split, not new scope: it makes the emergent competition visible in a run
— dosing 80 mg/L SO₂ still *strongly* suppresses MLF (~0.13 g/L malic slips through, pH rise 0.12→
0.005), but the transient acetaldehyde peak (free crashes to ~0.9 mg/L near day 2) briefly relaxes
suppression, so it is not a perfect block. `test_so2_dose_suppresses_mlf_in_a_run` updated to this
faithful behaviour (threshold `>3.9`→`>3.8`, with the mechanism documented — *not* a weakening for CI).

**Readout-only, like D-22 (the deferred RHS coupling).** The split does **not** feed back into the
acetaldehyde reduction — bound acetaldehyde is notionally protected from ADH, but the D-27 reduction
still consumes it. That RHS coupling (and SO₂'s own bisulfite back-reaction on pH, still deferred
from D-22) is the scoped omission, caveated. Isolability holds **on a run with no live consumer**:
`so2_total` is inert (no Process touches it) and free/bound are pure readouts, so on an MLF-*off*
ferment dosing SO₂ leaves every other state column and pH byte-for-byte and carbon closing (the
D-22 isolability test survives verbatim under the rename). This is **conditional, not
unconditional**: once MLF is dosed, SO₂ *does* change the trajectory — that is the whole point of
the gate, pinned by `test_so2_dose_suppresses_mlf_in_a_run`. The two tests together are the honest
statement: SO₂ is inert until a consumer reads it, then it acts through that consumer alone.

**The parameter.** New `K_acetaldehyde_so2 = 1.5e-6 mol/L` in the shared `acidbase.yaml`, tier
**plausible**, band `[2e-7, 2.1e-6]` (order-of-magnitude literature scatter). Source: Burroughs &
Sparks (1973), the canonical carbonyl-bisulphite dissociation constants; apparent Kd 1.5e-6 (pH 3.3)
–2.06e-6 (pH 3.5) across the wine literature; Blouin (1966) "~0.04 % free acetaldehyde at 30 mg/L
free SO₂" as a shape anchor. **Basis pinned (advisor-flagged as load-bearing):** referenced to
bisulfite; the literature apparents are usually per *total free* SO₂ at a stated pH, but at wine pH
bisulfite is ~0.94–0.99 of free SO₂ (`bisulfite_fraction`), so the two bases differ ≤5 % — inside
the band. Honest overclaim caveated: acetaldehyde is the *principal* but not sole binder (pyruvate,
α-ketoglutarate, sugars also bind), so modelled `bound` under-estimates and the "total" slot is
really "free + acetaldehyde-bound" — free/molecular slightly over-estimate the protective pool.

**Tier.** `molecular_so2_tier` now folds in the binding `K` alongside both pKa sets, floored at
`PLAUSIBLE` (never `VALIDATED`; covers free/bound/molecular alike). **Emergent + verified before
the acceptance test:** dosing 50 mg/L total SO₂ into a 20 °C wine ferment, free SO₂ dips 50 → 0.9
mg/L at the acetaldehyde peak (day 1.7) and recovers to 50; `free + bound == total` to machine
precision at every column. 7 new tests (+1 MLF assertion tightened); **349 green**, ruff + mypy
clean. **Next in the beat:** H₂S (carbon-free, inverse-low-N gate — the accounting-easiest,
following the SO₂ precedent).

## Deferred (decide early in the relevant milestone)

- ~~**pH / acid model richness**~~ — **decided in D-18** (full charge-balance solver),
  built after the byproducts beat; **solver landed 2026-06-30** (`core.acidbase`,
  `fermentation.analysis`) — see D-18 "Resolution".
- ~~**Stochastic ensemble API**~~ — **decided in D-24 and IMPLEMENTED 2026-07-01**
  (`runtime/ensemble.py`): triangular Monte-Carlo over the `Uncertainty` bands, scoped to
  the active Process set's reads, nominal + median + P5/P95 band, per-member conservation.
- **Packaged parameter-data access:** tests read YAML via filesystem path. If we
  ship a wheel that must read its own data, switch to `importlib.resources`.
