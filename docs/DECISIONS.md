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

## D-29 — Hydrogen sulfide (H₂S): a carbon-free produced pool with an inverse-nitrogen gate

**Status: IMPLEMENTED 2026-07-01** (364 green + 5 benchmark). The §3.2 aroma beat after the SO₂
free/bound split (D-28), and — as the D-27/D-28 forward notes anticipated — the **accounting-
easiest** metabolite yet. H₂S ("rotten egg", sensory threshold ~1–2 µg/L) is released by the
yeast sulfate-reduction sequence: sulfate → sulfite → sulfide, where the sulfide is normally
fixed onto **nitrogen** skeletons (O-acetylserine/-homoserine) to build cysteine/methionine.
When yeast-assimilable nitrogen (YAN) runs low there is no acceptor, so sulfide is excreted as
H₂S — **de-repression at low nitrogen**, the exact inverse of the Ehrlich fusel gate (`N/(K_n+N)`,
D-19).

**The model — one Process, one produced-only pool, carbon-free.** New `h2s` state slot (g/L,
`default=0.0`, in `_common_specs` ⇒ **both** media). One additive Process
`HydrogenSulfideProduction` (`core.kinetics.hydrogen_sulfide`):

    d(h2s)/dt = k_h2s · X·S_total/(K_sugar_uptake + S_total) · K_h2s_n/(K_h2s_n + N)

flux-linked (shares `K_sugar_uptake`, so it stops at dryness — the sulfate-reduction machinery
runs while the cell ferments), inverse-N gated (`K_h2s_n/(K_h2s_n+N)`: ~0 when N replete, → 1 as
N → 0), and held **temperature-flat** (documented v1 simplification, like the α-acetolactate
excretion D-26 and the acetaldehyde production D-27). Intrinsic yeast metabolism, so wired into
both media (its own isolable `_H2S_PROCESSES` tuple; unlike the *dosed* MLF organism). Params in
a new shared, medium-agnostic `hydrogen_sulfide.yaml` (sulfate reduction is generic yeast
metabolism), merged at the compile seam alongside `acetaldehyde.yaml`/`vicinal_diketones.yaml`.

**Why a separate `K_h2s_n`, not the growth `K_n`.** The gate half-saturation is a **new
parameter on the YAN scale** (`0.1 g/L`, speculative, band `[0.05, 0.2]`), *deliberately distinct*
from the growth `K_n` (`0.0088 g/L`). Reusing the growth constant would make a razor-edge gate
that opens only in a thin sliver at near-zero N; the YAN-scale constant makes the repression a
smooth, physiologically-relevant function across a must's nitrogen range (H₂S-management practice
targets YAN ≳ 140–150 mg/L; Ugliano 2009; Jiranek/Henschke). `k_h2s = 2e-6 /h` (speculative, band
`[5e-7, 1e-5]`) sizes cumulative produced ~0.5 mg/L for a default low-YAN wine.

**The most isolable beat in the model — but precisely stated.** H₂S is **carbon-free** (registered
with 0 carbon in `chemistry`, like SO₂), so it sits on **no conservation ledger** (its sulfur is
untracked, exactly as free SO₂'s is — there is no sulfate/sulfur state) and needs **no new
conservation code**; carbon still closes to machine precision on a compiled run with H₂S wired in.
The Process **touches only `h2s`** and merely *reads* `X`/`S`/`N`, so disabling it leaves the
**RHS of every other column byte-for-byte identical — verified *exactly* (0.0)** across states
(`test_isolable_at_derivative_level`). The *integrated* trajectory then drifts by only ~1e-7
relative — a **pure adaptive-solver mesh artifact** (adding the `h2s` equation shifts the error-
controlled step selection), **not a physical coupling**; this is cleaner than the acetaldehyde
buffer (D-27), whose `E` write feeds a *genuine* second-order `E`→viability perturbation on top
of the mesh effect. The advisor predicted "byte-for-byte"; the empirics refined it to
byte-for-byte *at the RHS* + a ~1e-7 mesh artifact at the trajectory level (both pinned:
`test_isolable_at_derivative_level`, `test_trajectory_isolability_is_solver_mesh_only`). **No tier
headline either:** unlike the diacetyl decarb (writes shared `CO2`, D-26) and acetaldehyde
production (writes `E`, D-27), this writes a **fresh pool nothing reads**, so no other column's
structural tier drops. All-speculative.

**The load-bearing empirical finding (checked BEFORE writing the acceptance test — the D-26
checkpoint discipline).** The advisor flagged, and a run confirmed, that the defining real
behaviour — *low-YAN must ⇒ far more H₂S* — is only **partially** reproduced, because the
nitrogen model strips `N` to ~0 by **day ~1.3 regardless of dose** (the known no-residual-N-floor
gap, D-23). Once N = 0 the inverse gate is ~1 for the rest of the ferment for **every** must, so
the **cumulative endpoint lever is muted**: 80 / 150 / 300 mg/L YAN → 0.557 / 0.542 / 0.527 mg/L
(direction right, only ~5 %). So the acceptance test does **not** assert a hollow
`low_final ≫ high_final`. What *does* emerge cleanly and is the anchor: **the gate direction**,
tested two honest ways — (1) at the derivative level, rate(low N) > rate(high N) at fixed flux
(`test_inverse_nitrogen_gate_direction`); (2) integrated and cross-must, the low-YAN must produces
**~1.8× more H₂S by day 1** *even though it grows less biomass* (2.14 vs 1.70 g/L X, so the gate
wins over the higher flux — not a flux artifact). The muted endpoint is pinned as *small on
purpose* (`test_cross_must_endpoint_lever_is_muted`), documenting the gap rather than papering
over it.

**Scope (v1) / deferred.** **Produced-only** — no CO₂-stripping volatilization sink yet, so `h2s`
is *cumulative produced*, which **overstates residual** (real fermentation sweeps most H₂S out
with the CO₂ stream to µg/L residuals). The stripping sink is the deferred follow-up — the exact
ester **D-19 (produced-only) → D-20 (Henry's-law sink)** precedent. Yeast-pathway (sulfate-
reduction) H₂S only; other sulfides/mercaptans and copper-binding are out of scope. The full
cross-must YAN lever unlocks only when the **residual-N floor** lands (a separate nitrogen-model
beat; see Deferred). New `M_H2S` in `chemistry` (0 carbon). 15 new tests; **364 green** + 5
benchmark, ruff + mypy clean. **Next in §3.2:** the aroma beat is essentially complete (esters,
fusels, VDK/diacetyl, acetaldehyde, SO₂ speciation, H₂S); candidates are the H₂S CO₂-stripping
sink or the residual-N floor (which would make this beat's cross-must lever real).

## D-30 — Residual-nitrogen floor: an opt-in biomass carrying-capacity cap on growth

**Status: IMPLEMENTED 2026-07-01** (380 green + 5 benchmark, ruff + mypy clean). The D-29 forward
note's "residual-N floor" candidate, chosen to make the muted H₂S cross-must lever real. Closes
the nitrogen-model gap surfaced repeatedly since D-23: `GrowthNitrogenLimited` is the **sole**
nitrogen sink and its only shutoff is a tiny-`K_n` Monod term, so a wine ferment builds
`X ≈ X0 + N0/f_N` and strips yeast-assimilable nitrogen (YAN) to ~0 by **day ~1.3 regardless of
dose**. That erases every downstream low-N signal — most visibly the D-29 H₂S inverse-N gate,
which reads `N→0` for every must (lever muted to ~5 %).

**The mechanism — a logistic carrying-capacity RateModifier.** Real yeast populations saturate
*below* the nitrogen ceiling (oxygen/sterol limitation, density effects), leaving YAN unconsumed.
The textbook lumped form is a logistic cap: growth slows as biomass `X` nears a capacity `K` and
stops at it. Because this **scales** an existing flux rather than adding one, it is a
`RateModifier` (`core.kinetics.carrying_capacity.BiomassCarryingCapacity`), not a summed Process:

    factor(X) = clamp(1 − X/K,  0, 1)        K = biomass_carrying_capacity

multiplied onto `GrowthNitrogenLimited`'s **whole** contribution by `ProcessSet`. Linear `1−X/K`
(not the smoothed `(1−·)**n` ethanol wall) is deliberate — `X` self-limits (growth→0 as `X→K`),
so the state never gets driven past the wall and there is no derivative kink; the `[0,1]` clamp
still guards a solver overshoot `X>K` from flipping the factor negative (which would make growth a
biomass/nitrogen *source*). **Conservation is automatic:** scaling growth's whole contribution by
one scalar preserves `dN = −f_N·dX` and the proportional carbon-skeleton draw, so `total_nitrogen`
and `total_carbon` still close to solver tolerance with the cap on (`test_carbon_and_nitrogen_
close_with_the_cap_on`) — the nitrogen simply stays in the `N` pool once growth saturates. This is
the crux that makes a cap the right vehicle: less biomass, exact balances, residual N left behind.

**Why OPT-IN, not default — the fundamental Coleman conflict.** Coleman, Fish & Block (2007), the
keystone wine model, has **no** biomass cap: it consumes all YAN and builds full N-proportional
biomass at every dose, and `test_coleman_reconstruction` confirms our core reproduces that
line-for-line at 80 **and** 330 mg N/L. A pre-check (the D-26 checkpoint discipline: measure
before writing) established the tension is **not** a mechanism artifact but fundamental —
*restoring the H₂S lever requires residual **assimilable** N that differs by dose, which means not
consuming it, which means departing from Coleman's zero-residual biomass curve.* No mechanism
escapes this (a non-assimilable/proline split keeps Coleman intact but leaves assimilable N at
zero ⇒ lever still muted). The measured cost of turning the cap on in the default wine: Coleman
RMSE 1.35→up to 9.35 (80 mg/L) and 1.20→up to 27.84 (330 mg/L) vs the <2.0 gate. So per prime
directive #3 the cap ships **isolable and disabled by default**: wired into the wine medium but
the compile seam **disables** it unless a scenario opts in via `carrying_capacity_gpl`. Disabled ⇒
factor 1 **and** excluded from tier derivation (`ProcessSet` counts enabled, not nonzero,
modifiers — the wine-only MLF *tier* isolability argument, extended to the multiplicative path),
so an undosed wine run is **byte-for-byte the validated core** (verified *exactly* 0.0 across
states, `test_disabled_cap_equals_the_uncapped_rhs_exactly`) and growth stays PLAUSIBLE. Opt in
and growth's `X`/`S`/`N` **structural** tier drops PLAUSIBLE→SPECULATIVE, honestly flagging the
departure — no param-aware headline (growth already reads the speculative `K_s`, the D-26/D-27
pattern). Coleman reconstruction, §2.2 dryness/ABV, fusel/ester benchmarks all untouched.

**Provenance + seam.** New `biomass_carrying_capacity` in `wine_generic.yaml` (**speculative**,
`author estimate`, `2.5 g/L`, band `[2.0, 5.0]` — the cap must bite below the ~2.6–3.0 g/L
uncapped biomass to leave residual). The value 2.5 is the pre-check cap that restored the lever
while leaving ~0 residual at low YAN (the correct clinical picture). New optional wine scenario
key `carrying_capacity_gpl`: **presence enables** the modifier; its value **overrides** the YAML
reference (so a demonstration can sweep `K`), injected at the compile seam via
`_override_carrying_capacity` (mirrors the D-14 N-yield override).

**Emergent, verified.** With the cap on (K=2.5): the H₂S endpoint is monotone in dose
(80 > 150 > 300 mg/L YAN) and its **span widens materially versus the muted core**
(`test_cap_restores_the_h2s_cross_must_lever`, asserted as ordering + ratio, not brittle absolute
values); a **dose-dependent residual YAN** survives — low-YAN musts still (nearly) exhaust N while
high-YAN musts end well above (`test_cap_leaves_dose_dependent_residual_nitrogen`), the correct
clinical picture the core (~0 at every dose) cannot show. A capped wine still ferments to dryness
(`test_opt_in_wine_still_reaches_dryness` — less biomass slows the tail but per-cell uptake keeps
going).

**Scope (v1) / honesty.** **Wine-only** (the H₂S lever and the prospective MLF-with-growth model
are wine concerns), mirroring the wine-only MLF wiring; beer carrying capacity is deferred. The
**MLF unblock is PROSPECTIVE, not delivered**: MLF v1 is conversion-only with pH/ethanol/
molecular-SO₂/cardinal-T gates and **no nitrogen gate** (D-23), so residual N does *not* change
current MLF behaviour — it enables a *future* MLF-with-growth model. 16 new tests. Next §3.2
candidate remaining: the H₂S CO₂-stripping sink.

## D-31 — MLF-derived diacetyl: *Oenococcus oeni* citrate co-metabolism + bacterial reduction

**Status: IMPLEMENTED 2026-07-01** (395 green + 5 benchmark). The real coupling MLF (D-23) makes
available and the deferred half of the diacetyl story (D-26 built the *yeast* valine-pathway
diacetyl only). Alongside malate, *O. oeni* co-metabolises **citric acid**, overflowing
α-acetolactate that decarboxylates to **diacetyl** — the buttery note that defines many post-MLF
(esp. barrel-aged Chardonnay) wines, and a real winemaking control point (co-inoculation, lees
contact, and post-MLF SO₂ timing all move it). Two new *O. oeni* Processes in
`core/kinetics/malolactic.py`, wired into the wine-only `_MLF_PROCESSES` tuple and disabled with
the malate conversion at the compile seam when *O. oeni* is un-pitched.

**Owner decisions (the three forks put up front, discuss-before-build).** (a) **citrate is a must
input** (`citrate_gpl`, like `malic`/`tartaric`), so the level is a per-scenario lever; (b) the
carbon routes **via the shared α-acetolactate reservoir** (`citrate → α-acetolactate + CO₂`),
reusing the always-on D-26 decarboxylation + reduction so diacetyl *emerges* rather than being a
second pathway; (c) **add O. oeni's own diacetyl reduction now** (not deferred), so lees-contact
clearing is modelled.

**Why a citrate pool at all (the load-bearing scope decision).** MLF-diacetyl is a late-MLF,
often **post-dryness** phenomenon, so its carbon **cannot** come from sugar: the yeast VDK
stand-in draws α-acetolactate carbon out of `S` via `draw_carbon_from_sugar`, which correctly
**no-ops at `S=0`** — sourcing from an empty sugar pool would either strand carbon (breaking
`total_carbon` closure) or stop diacetyl production exactly when this beat needs it. Citrate is
present independent of sugar, so a dosed `citrate` slot (C6H8O7, added to `chemistry.py`,
`total_carbon`, the wine schema, and the compile vocabulary) is the **floor** for honest carbon
closure here, not scope creep (the advisor's decisive framing, confirming the finding).

**Stoichiometry is a lumped fiction — owned.** `MalolacticCitrateMetabolism`:
`d(citrate) = −r_c·M_citric`, `d(acetolactate) = +r_c·M_acetolactate`, `d(CO2) = +r_c·M_CO2` with
`r_c = k_citrate·X_mlf·[citrate]/(K_citrate+[citrate])·gate`. Citric acid (6 C) → α-acetolactate
(5 C) + CO₂ (1 C), so **carbon closes mole-for-mole (6 = 5 + 1)** on the existing ledger, exactly
like malic → lactic + CO₂ (D-23). *Mass* carries a small gap (192.124 ≠ 132.116 + 44.009), so
carbon is the invariant (as for the VDK decarb / beer hydrolysis water, D-8). CAVEAT: real citrate
metabolism is `citrate → acetate + oxaloacetate → pyruvate + CO₂`, ~2 citrate per α-acetolactate,
with **acetate** (a volatile-acidity contributor) the *dominant* co-product. The single-reaction
stand-in drops the acetate/lactate branches; `k_citrate` is held **low so citrate stays mostly
unconsumed** (~6 % at the reference dose) — the *trace diacetyl branch only*, which keeps the
fiction honest (we do not claim to resolve citrate's full fate).

**Rate — citrate's own Monod × the SHARED environmental gate (NOT malate's `r`).** A new helper
`malolactic_environmental_gate` factors out `g_pH·g_EtOH·g_SO₂·γ(T)`, now called by *both* the
malate conversion (a byte-equivalent refactor) and the citrate branch — so SO₂/ethanol/low-pH
arrest citrate metabolism just as they arrest MLF. Coupling to citrate (not the malate turnover)
is deliberate: malate's rate → 0 at malate depletion, which would kill exactly the post-malate
diacetyl peak this pool exists to capture. Each Process solves pH once (a second `brentq` only on
dosed runs — acceptably cheap, not optimised away).

**Bacterial reduction (owner's fork c).** `OenococcusDiacetylReduction`:
`L = k_mlf_diacetyl_reduction·X_mlf·f(T)·[diacetyl]` (shared `E_a_reduction`), a mole-for-mole
C4 → C4 transfer to `butanediol` like the yeast reducer (D-26). It **complements** the yeast
`DiacetylReduction`: in co-inoculation the yeast (higher rate × biomass) clears diacetyl fast
while viable; this bacterial reducer keeps clearing it after the yeast is ethanol-inactivated, as
long as *O. oeni* is present — the realistic lees-contact clean-up, and the reason removing the
bacteria (SO₂ / racking) locks diacetyl in. **Consequence flagged (advisor):** with *O. oeni*
dosed, MLF-diacetyl is **not permanently stranded** in v1 (`X_mlf` is a constant, never killed);
the "package/rack early ⇒ diacetyl locked in" case needs a racking event to remove `X_mlf`,
deferred to the event loop with the bacterial death/arrest gate (as for MLF conversion, D-23).

**Emergent, verified.** Dosing *O. oeni* + citrate lifts wine diacetyl clearly above the
yeast-only baseline (peak ~0.28 vs ~0.10 mg/L, ~2.8×, into the buttery range above the ~0.2 mg/L
threshold), with a **late peak** (~day 5–6, via the reservoir decarb lag, past the early
low-ethanol conversion window) that then **falls** as reduction clears it — the buttery-then-
cleaning-up MLF signature. A larger *O. oeni* dose leaves a lower final/peak ratio (bacterial
clearing). `total_carbon` closes to machine precision throughout.

**Isolability (prime directive #3).** Both Processes are in the dosed, disabled-when-unpitched
`_MLF_PROCESSES` tuple: an un-pitched (or citrate-free) wine run is **byte-for-byte** the prior
core, and citrate dosed *without* O. oeni sits inert (diacetyl matches the yeast-only baseline).
`citrate` keeps its **VALIDATED** tier when un-pitched (nothing active touches it) and drops to
**speculative** when dosed — the exact `malic`/`lactic` pattern (D-23). All new params
(`k_citrate`, `K_citrate`, `k_mlf_diacetyl_reduction`) are **speculative** order-of-magnitude
estimates in `wine_generic.yaml`; both Processes are speculative (`acetolactate`/`diacetyl` were
already speculative from the yeast VDK pathway, so no new tier headline). Citrate is **carbon-
active but not charge-active** — kept out of the D-18 pH balance in v1 (a scoped omission the
inverse anchoring absorbs at t=0, as for SO₂'s bisulfite charge, D-22). 14 new tests.

**Scope (v1) / deferred.** The dominant citrate → acetate/lactate branches and full citrate
depletion (the single-reaction stand-in); the bacterial arrest/death gate and a racking event
(so "SO₂ locks diacetyl in" and permanent stranding are not yet demonstrable); citrate in the pH
charge balance. These follow the MLF-with-growth beat (D-23) and the event loop.

## D-32 — Amino-acid ledger: a nitrogen-anchored, modifier-scaled biomass swap

**Status: IMPLEMENTED 2026-07-01** (406 green + 5 benchmark). Builds the toggleable
amino-acid ledger D-23 scoped and deferred (the separate yeast/AF beat). Yeast build biomass
mostly from amino acids, but the validated core sources *all* biomass carbon from sugar and
*all* biomass nitrogen from the lumped ammonium `N` pool, and `N` is deliberately carbon-free
in `total_carbon` (D-19). Making amino acids a carbon source is thus a change to the protected
carbon *and* nitrogen ledgers — restored to isolability by the owner's `default=0` `amino_acids`
pool, implemented (advisor's refinement, D-23) as a **separate isolable swap Process**, not a
branch in the core's hottest kinetic.

**The swap (`core/kinetics/amino_acids.py` `AminoAcidAssimilation`, wine-only).** For biomass
built at the shared `biomass_growth_rate` (extracted from `GrowthNitrogenLimited` so the swap
anchors to the *identical* rate), it consumes amino acids at `ρ` and **debits** the pool
(`d[amino_acids]=−ρ`), **refunds ammonium** (`d[N]=+ρ·y_N`), and **refunds sugar carbon**
(`d[S]+=+ρ·y_C`), leaving biomass `X` untouched. It is a pure transfer aa→S (carbon) and aa→N
(nitrogen), so **carbon- and nitrogen-neutral by construction** for any `ρ` — the pool is now
weighted in *both* `total_carbon` (arginine C-fraction) and `total_nitrogen` (arginine
N-fraction, the first per-species nitrogen accounting: new `NITROGEN_ATOMS` +
`nitrogen_mass_fraction` in `chemistry.py`). **Bookkeeping caveat (the D-19/D-31 stand-in
discipline):** mechanically the aa carbon is refunded to *sugar* (biomass carbon still comes
from growth's sugar draw, and the spared sugar ferments to ethanol) — arginine's carbon skeleton
is booked as spared hexose, not tracked through arginine catabolism. Carbon-closing and
defensible (aa-fed biomass really spares sugar for ethanol), but a stand-in; one consequence is
that dosing aa nudges ethanol up ~0.15–0.3 % of sugar. The §2.2 benchmarks run undosed, so they
are untouched.

**Nitrogen-anchored rate, N-rich representative (the load-bearing choices).** Amino acids *are*
part of YAN, so `ρ = ψ·gate(aa)·f_N·base_dx/y_N` with `gate(aa)=aa/(K_amino_acids+aa)` and
`ψ = amino_acid_assimilation_fraction ∈ [0,1]`. The advisor's decisive framing: **carbon
over-refund is non-physical** (creates hexose from amino acids = gluconeogenesis, which
fermenting yeast do not do) but **nitrogen over-refund is physical** (deamination of surplus aa
to ammonium). Anchoring on nitrogen makes the N refund `ρ·y_N = ψ·gate·f_N·base_dx ≤ f_N·base_dx`
(never over-refunds N, so no deamination branch in v1), and picking an **N-rich** representative
amino acid — **arginine** (C₆H₁₄N₄O₂, the dominant *assimilable* grape amino acid, mass C:N ≈
1.29 ≪ biomass's `f_C/f_N` ≈ 4.3) — makes the carbon refund `≈ 0.30·ψ·gate ≤ 0.30` of growth's
sugar-carbon draw for **all** ψ ≤ 1. So the carbon cap never binds: no clamp, no C⁰ kink for the
stiff BDF solver, no sugar creation. A carbon-rich amino acid (leucine ≈ 5.1) would sit at the
edge and force a clamp — the species choice *is* what keeps v1 clean.

**The correctness crux — modifier scaling (option 2, advisor-forced).** The safety above uses
growth's *pre-modifier* `base_dx`, but growth's realised biomass is `base_dx·M` where `M` is the
Arrhenius × (opt-in) carrying-capacity `RateModifier` product `ProcessSet` applies. A swap
refunding at `base_dx` while growth draws at `M·base_dx` would, at `M < 1` (cold ferment, or the
carrying cap near saturation with nitrogen still available — the D-30 residual-N regime), refund
more than the draw and **create sugar**. The fix: the wine growth Arrhenius (`for_growth` gains
an `*also_scales` target) *and* the carrying-capacity modifier (its `modifies` now names the
swap) scale the swap too, so refund and draw carry the same `M`:
`net dS = M·f_C·base_dx·(0.30·ψ·gate − 1) ≤ 0` and `net dN = M·f_N·base_dx·(ψ·gate − 1) ≤ 0`.
This was landed **fail-first** per the advisor: the guard tests (`net dS/dN ≤ 0` at a carrying-
saturation state; the swap refund scaling with the growth Arrhenius factor) were written to FAIL
with the unscaled swap and confirmed failing (`dS = +0.0279`, arrhenius ratio 1.0 vs 0.445),
then pass once the scaling landed. At `T_ref` `M = 1` and the mismatch never fires — the reason
a naive T_ref-only test would be vacuous.

**Isolability (undosed-only).** The compile seam disables the swap when `amino_acids_gpl ≤ 0`, so
an undosed wine run is byte-for-byte the validated core, the empty `amino_acids` slot keeps its
VALIDATED tier, and folding the swap into the two modifiers' `modifies` is transparent (a modifier
naming a zero-contribution/disabled Process is a no-op). **Dosed**, the swap *correctly* perturbs
the run: refunded N/S raise the pools growth reads on the next step, so dosing amino acids behaves
like **supplementary YAN** (nitrogen lasts longer ⇒ more biomass) — a second-order feedback, not a
first-order growth edit (growth's derivatives are untouched); the "byte-for-byte" claim is thus
undosed-only, and the swap's speculative tier drops growth's `S`/`N` outputs to speculative when
enabled (the D-26/D-30 structural-drop pattern). New speculative params
`amino_acid_assimilation_fraction=0.5`, `K_amino_acids=0.1 g/L` (author estimates) +
`amino_acids_gpl` scenario key. 11 new tests.

**Scope (v1) / deferred.** The swap only (primary-fermentation yeast carbon/nitrogen honesty).
The D-19 **fusel Ehrlich re-route** (drawing fusel carbon from this pool instead of its sugar
stand-in) is the natural later home but needs the deamination branch (Ehrlich releases aa-N), so
it is deferred with that carbon-anchored + explicit-deamination generalisation (D-23). MLF-growth
funded from this pool stays blocked on an **autolytic-peptide source** to refill it post-AF (the
pool is empty at the MLF pitch point, D-23). Wine-only; beer deferred with the wine-only nitrogen
model (D-30).

## D-33 — Fusel Ehrlich re-route: sourcing fusel carbon from amino acids, with deamination

**Status: IMPLEMENTED 2026-07-01** (417 green + 5 benchmark). Builds the first of the two
prerequisites the still-blocked MLF-with-growth beat was deferred on (D-23/D-32): the D-19 fusel
Ehrlich re-route, now that the amino-acid pool (D-32) gives the model a carbon- *and*
nitrogen-bearing amino-acid source and the **deamination branch** can therefore close.

**The gap.** :class:`FuselAlcoholsEhrlich` (D-19) books fusel carbon out of *sugar* — a documented
stand-in, because the real Ehrlich pathway builds higher alcohols from *amino-acid* skeletons
(transamination → decarboxylation → reduction) and releases the amino group as ammonium. Sugar was
used only because `N` (YAN) carries no carbon in `total_carbon` (D-19), so there was nowhere else to
draw it from. The `amino_acids` pool (arginine; D-32) removes that constraint.

**The mechanism — a separate wine-only swap (`core/kinetics/byproducts.py`
`FuselAminoAcidReroute`).** Mirroring the D-32 `AminoAcidAssimilation` swap, production stays entirely
in the producer; the re-route only moves the carbon *source*. For the amino-acid-sourced fraction
`g = aa/(K_amino_acids+aa)` (the same smooth availability gate the swap uses) of the fusel carbon
`F_c = rate·c_fusel`, it **refunds sugar** by `g·F_c` (undoing the producer's draw for that fraction),
**debits amino acids** by `g·F_c/c_aa`, and **releases ammonium** `N` by `(g·F_c/c_aa)·y_N` — the
deamination branch. Carbon closes (the fusel's `F_c` is now `(1−g)·F_c` from sugar + `g·F_c` from
amino acids); nitrogen closes (amino acids lose exactly the nitrogen `N` gains). Net sugar is
`−(1−g)·F_c ≤ 0` for all `g ≤ 1`, so it never creates sugar (spared-sugar→ethanol is the D-32
bookkeeping caveat). The producer and re-route share one `fusel_production_rate` helper (extracted
this beat) so the sugar refund matches the draw to machine precision — via a shared refund/draw pair
(`refund_carbon_to_sugar`, the inverse of `draw_carbon_from_sugar`, now the single source both the
swap and the re-route use).

**Why a separate Process was *forced*, not merely preferred (advisor).** Unlike the D-32 swap (whose
separation protected the *validated* growth kinetic), `FuselAlcoholsEhrlich` is already speculative
Tier-2. What forces the split is the **beer `touches` contract**: declaring `amino_acids`/`N` in the
both-media producer would raise at beer's `ProcessSet` construction (beer has no `amino_acids` slot).
So the re-route is wine-only and touches only `("S","amino_acids","N")` — **never `fusels`** (the
warm=more-fusel benchmark is untouched at the derivative level; verified).

**Not modifier-scaled (contrast D-32).** The swap is scaled by the growth Arrhenius/carrying modifiers
because it anchors to growth's *modified* rate. The re-route anchors to the *fusel* rate, which
carries its own `E_a_fusels` Arrhenius **inside** the Process and is scaled by no `RateModifier` — so
the re-route must also stay unmodified, and since both call the one shared helper and neither is a
modifier target, refund matches draw exactly with no D-32-style `M`-mismatch to guard.

**Documented lump — arginine over-releases nitrogen (advisor caveat).** Sourcing fusel carbon through
the N-rich representative amino acid deaminates `c_fusel/c_aa·y_N ≈ 0.78 g N per g fusel-carbon` —
roughly **4× the real leucine→isoamyl-alcohol N:C** (leucine carries one amino group over six
carbons). Conservation-exact, but a forced consequence of the single-species `amino_acids` lump
(arginine, chosen N-rich for the D-32 swap), the same class of stand-in as the sugar-carbon fiction it
replaces. The released N feeds back as supplementary YAN, but fusels are trace so the effect is
second-order and tiny.

**Isolability (undosed-only, paired with the producer).** The availability gate → 0 at `aa = 0`
(byte-for-byte the sugar-stand-in producer on an undosed run) and the compile seam disables the
re-route with the swap when `amino_acids_gpl ≤ 0` (tier isolability; the empty `amino_acids` slot
keeps VALIDATED). It is only valid while `FuselAlcoholsEhrlich` is active — it refunds sugar that
producer drew — so the two are kept paired (the same acceptable swap↔producer coupling as D-32's
swap↔growth; disabling the producer alone would let the re-route create sugar). No new parameters
(reuses `K_amino_acids`); 9 new tests (`tests/test_fusel_reroute.py`), 417 green + 5 benchmark, ruff
+ mypy clean, §2.2 undosed trio unchanged.

**STILL-DEFERRED for MLF-growth.** This closes the *fusel* half of the D-32-deferred pair. The other
prerequisite — an **autolytic-peptide source** to refill the amino-acid pool post-AF (it is empty at
the MLF pitch point, D-23) — is D-34. MLF-growth itself (a bacterial growth Process consuming the
pool + the event loop) stays deferred beyond both. The full deamination generalisation (a standalone
excess-aa deamination flux, vs this fusel-coupled release) also remains future work.

## D-34 — Yeast autolysis: the autolytic-peptide source that refills the amino-acid pool

**Status: IMPLEMENTED 2026-07-01** (428 green + 5 benchmark). Builds the *second* of the two
prerequisites the still-blocked MLF-with-growth beat was deferred on (D-23/D-32) — the first being
the D-33 fusel re-route. *O. oeni* builds biomass from amino acids/peptides, but the `amino_acids`
pool (D-32) is **empty at the MLF pitch point**: the same yeast uptake that strips `N` to ~0 by day
~1.3 would strip any dosed amino acids too (the empirical finding that settles D-23). Real wine
refills the pool by **autolysis** — dying yeast self-digest and release intracellular amino acids
(the basis of *sur lie* aging). :class:`YeastAutolysis` is that flux: the **first consumer of the
`X_dead` pool** (dead biomass from D-13 ethanol inactivation), turning it into assimilable
`amino_acids`.

**The conservation problem, and why a debris pool (advisor-decided; this was the one blocking fork).**
Dead biomass is **carbon-rich** (mass C:N `f_C/f_N` ≈ 4–11 across Coleman's nitrogen range) while the
assimilable amino acids it releases are **nitrogen-rich** (arginine mass C:N ≈ 1.29). So per gram of
nitrogen liberated, biomass gives up 4–11 g of carbon but arginine holds only ~1.3 g — **most of the
dead-cell carbon cannot leave as amino acids.** The advisor settled the excess-carbon sink decisively:
**not CO₂** (that would falsely claim autolysis *respires* the cell — it is enzymatic self-digestion,
not respiration — and would perturb a benchmarked pool; ~86 % of dead-cell carbon would be wrongly
mineralised), but a **carbon-only `debris` pool** (booked as glucan, C6H10O5). This is the physically
*dominant and correct* fate: yeast cell walls (β-glucans/mannoproteins) are ~30 % of dry mass and are
exactly the non-assimilable material that stays as lees. The `debris` pool is weighted in
`total_carbon` only (nitrogen-free — all released N goes to amino acids), the `esters_gas` idiom (a
bookkeeping pool carrying carbon that has left the metabolite pools but not the atom balance).

**The flux — nitrogen-anchored, first-order (`core/kinetics/autolysis.py`).** With `r = k_autolysis ·
arrhenius(T, E_a_autolysis, T_ref) · X_dead` [g X_dead/L/h] (autolysis is enzymatic, so warmer lees
clear faster): liberate the dead-cell nitrogen as amino acids (`d[amino_acids] = +r·f_N/y_N`, arginine
carrying exactly `r·f_N`), debit dead biomass (`d[X_dead] = −r`), and route the C-rich remainder to
debris (`d[debris] = (r·f_C − r·f_N·y_C/y_N)/c_debris`). Carbon closes (dead-cell carbon `r·f_C` splits
into the amino acids' carbon and the debris carbon); nitrogen closes (`r·f_N` is exactly what the
amino-acid pool gains; debris is N-free) — both to machine precision, verified at the RHS level and
over full runs. The excess-carbon split is **structurally non-negative** (biomass C:N always exceeds
arginine's over the whole `f_N` range 0.039–0.114), so `f_C > f_N·y_C/y_N` always and the split never
flips — **no clamp, no C⁰ kink** for the BDF solver (advisor-confirmed). The Process reads `f_N`/`f_C`
from params (so the compile-time Coleman override, D-14, flows through) and the conservation tests pull
them from `param_values`, not the raw YAML (advisor).

**Isolability — opt-in (the D-30 carrying-capacity pattern).** Unlike the always-on intrinsic aroma
pools, autolysis *consumes* core state (`X_dead`) and fills `amino_acids`/`debris`, so it measurably
perturbs the core and cannot be default-on without breaking the byte-for-byte guarantee and the §2.2
benchmarks. It ships **wine-only and disabled by default**: the compile seam enables it only when a
scenario passes `autolysis_rate_per_h` (which also overrides `k_autolysis`, letting a demonstration
sweep the *sur lie* timescale). Disabled ⇒ excluded from the derivatives *and* tier derivation (an
undosed wine run is byte-for-byte the validated core, verified). First guard `X_dead ≤ 0 ⇒ 0` (the
clamped first-order rate cannot overshoot negative). Wine-only, mirroring the wine-only `amino_acids`
pool / nitrogen model (D-30/D-32); beer deferred.

**Emergent (verified).** With autolysis on and amino acids un-dosed (so nothing consumes the pool —
the swap/re-route are compile-disabled), `X_dead` accumulates as the ferment ends and then feeds the
`amino_acids` pool, which **rises from empty** and keeps rising in the post-AF tail — the pool a later
MLF-with-growth model will draw on. `debris` outgrows `amino_acids` (most autolysed carbon is the
non-assimilable cell wall), the physically-right proportion.

**Tier: speculative** — first-order autolysis of dead biomass is a standard lumped form, but
`k_autolysis` (1e-3/h, ~29 d half-life, band [1e-4, 1e-2]) and `E_a_autolysis` (60 kJ/mol, band
[40k, 90k]) are author estimates and the single-amino-acid / carbon-only-debris lumping is a
simplification (real autolysate is a mix; mannoproteins retain some nitrogen). New species `glucan`
(C6H10O5) in `chemistry.py`; new wine-only `debris` slot (schema 25→26); new `autolysis_rate_per_h`
scenario key. 12 new tests (`tests/test_autolysis.py`) — including an advisor-ordered **three-way
composition** test (autolysis *feeds* the pool while the D-32 swap and D-33 re-route *drain* it, the
actual MLF-growth-prerequisite configuration every other test isolates apart: carbon + nitrogen close
over the full run) — 429 green + 5 benchmark, ruff + mypy clean, §2.2 undosed trio unchanged.

**STILL-DEFERRED — MLF-growth itself.** With both prerequisites now in hand (D-33 fusel re-route, D-34
autolysis refill), the remaining work is the *consumer*: an MLF-with-growth Process feeding a growing
`X_mlf` from the `amino_acids` pool, plus the **event loop** to pitch bacteria post-AF (runtime has no
event mechanism — the same block as sequential MLF, D-23). A standalone excess-amino-acid deamination
flux (vs the D-33 fusel-coupled release) also remains future work.

## D-35 — Event loop: segment-and-restart scheduling, and temperature as a driven ramp

**Status: IMPLEMENTED 2026-07-02** (449 green + 5 benchmark). The runtime gains its first
*time-driven* mechanism — the thing MLF-with-growth (D-23), a mid-ferment DAP/SO₂ dose, racking,
and a real temperature schedule all need but the model never had. Built in two parts on one
driver: the **verb-agnostic scheduling driver**, and the **temperature ramp** as its first client.
Discrete winemaking interventions (DAP/SO₂/racking/pitching) are the follow-up (D-36) on the same
driver.

**Scope fork (owner-decided, against the advisor's default).** Temperature scheduling was inert
too (only the earliest knot seeded the initial `T`; nothing drove it), so it was a live question
whether this beat includes it. The advisor recommended *deferring* temperature as a separate,
invasive "continuous forcing" beat. The owner chose to **do the ramp properly now** — and a
segmentation insight dissolved the invasiveness: a piecewise-*linear* schedule has a **constant
slope between knots**, so if the driver already restarts the integrator at breakpoints, temperature
is just a per-segment constant `dT/dt`. `T` stays an ordinary integrated state; every Process keeps
reading `y[T]` unchanged (the Arrhenius modifier was already written for a time-varying `T`, D-11);
nothing in core is refactored. The advisor endorsed this on reconsideration.

**The driver (`runtime/schedule.py`, `simulate_scheduled`).** Walks a run as segments separated by
`ScheduledEvent` breakpoints, calling the unchanged pure `simulate` on each segment and, at each
breakpoint, applying any of three opaque effects: a **state mutation** (`(schema, y) → y'`, a
dose/racking jump), an **in-place Process-set reconfiguration** (`enable`/`disable`, e.g. pitching
an organism mid-run), and/or a **parameter update** (a value in force from that time forward — how
the temperature slope changes per segment). It is **verb-agnostic**: it knows nothing about DAP or
temperature; the winemaking *vocabulary* + unit conversion live at the scenario compile boundary
(D-3), so runtime drives time, the boundary owns meaning, core stays pure physics.

Three properties are load-bearing:
- **Segment-and-restart, not `solve_ivp(events=…)`.** SciPy events detect zero-crossings and can
  terminate, but cannot mutate-and-resume. A dose is a genuine discontinuity, so the only correct
  approach is stop → jump → fresh `solve_ivp`. BDF re-initialising its order at each restart is
  *correct at a discontinuity*, not a perf bug. Because `dT/dt` is constant within a segment (a
  degree-1 polynomial), BDF integrates `T` **exactly** to round-off — verified `T(t)` matches the
  analytic line to `1e-10` — but *only because we segment at slope changes* (a segment spanning a
  slope discontinuity would not be exact).
- **External-flow ledger (conservation across a jump, a prime directive).** A dose injects mass
  from *outside* the system, so the single-run invariant becomes `final == initial + Σ inputs −
  Σ outputs`. Each mutation books its post-minus-pre **state delta** as an `ExternalFlow`; the
  continuous ODE still closes exactly *within* every segment, and the ledger is the correction term
  across the jumps. Booking the raw delta keeps the driver free of per-verb chemistry (the existing
  `total_carbon`/`total_nitrogen` weight it). The temperature ramp uses no mutations, so its ledger
  is empty — the machinery lands here, its winemaking payoff arrives in D-36.
- **Tier travels.** Per-segment `tier_map` snapshots are `combine`d (min) across segments, so a
  speculative Process enabled only for the back half of a run drags its variables to that tier for
  the *whole* trajectory (a run is only as trustworthy as its least-trustworthy segment).

Breakpoint times are emitted **once, post-mutation**, so a dose reads as a clean jump and the time
axis stays strictly monotone (downstream percentile/interp assume that). Same-instant events apply
in stable list order; events at `t0` seed the run before segment 0; events at/after `t_end` are
rejected (the boundary decides whether a late scenario intervention is an error). **Isolability:**
`events=()` is a single `simulate` call with identical arguments — byte-for-byte a plain run.

**The run chokepoint (`CompiledScenario.run`).** Storing `events` on the compiled scenario is not
enough — a hand-wired `simulate(cs.process_set, cs.param_values, cs.y0, cs.t_span_h)` *silently
ignores* them, and because the boundary injects `temperature_ramp_rate = slope_0` into
`param_values`, plain `simulate` would apply the *first* segment's slope for the whole run (correct
for a single-slope ramp, **wrong** for any multi-knot ramp or hold). So the compiled scenario grows
a single `run()` entry point that always dispatches through `simulate_scheduled(events=cs.events)`
— which, since `events=()` is byte-for-byte a plain `simulate`, is the right call for *every*
scenario (advisor-flagged gap; the same routing D-36 needs). **Caveat (deferred):** the stochastic
`simulate_ensemble` wraps the un-scheduled `simulate` and takes no `events`, so it shares the
multi-segment footgun; an ensemble-over-`simulate_scheduled` is a D-36 follow-up.

**The temperature ramp (`core/kinetics/temperature.py`, `TemperatureRamp`).** One Process,
`dT/dt = temperature_ramp_rate` (K/h), touching only `T`. Wired into **both** media (cellar
temperature is not a beverage property). The compile boundary (`_temperature_ramp_schedule`) turns
the `(day, °C)` knots into canonical hours/Kelvin, computes the piecewise-constant slope, and emits
a slope-change event only at interior knots where the slope **actually changes** — so **collinear
knots produce one segment** and a **flat/single-knot schedule produces none**. `T` is held (slope 0)
before the first knot and after the last. When (and only when) the schedule ramps, the boundary
mints a provenance-backed `temperature_ramp_rate` `Parameter` (the D-14/D-30/D-34 injection idiom)
for the first segment; later slopes ride the events. `CompiledScenario` gained an `events` field
carrying them.

**Reasoned deviation from the advisor on the disable-gate.** The advisor suggested *disabling*
`TemperatureRamp` when flat (mirroring the MLF/carrying/aa gates) for structural byte-for-byte. It
is instead **always enabled**, reading the rate with a `0.0` isothermal default and declaring **no
`reads`**. This gives the *same* two guarantees more simply: numerically, an un-ramped run adds
`0.0 + 0.0 == 0.0` to the `T` slot (byte-for-byte, verified against the untouched §2.2 benchmarks);
tier-wise, `tier_of("T")` is `combine([VALIDATED])` = VALIDATED — no drop, because both the Process
*and* the rate are VALIDATED (a set-point schedule is an exact input, not an empirical constant).
The advisor's disable-gate rationale ("an always-enabled ramp would drop `tier_of("T")`") only bites
when the Process/param is below VALIDATED, which is not the case here. Declaring no `reads` is
deliberate: the `reads` mechanism exists for D-1 *credibility* propagation, and a value exact by
construction borrows no credibility — declaring it would only force `temperature_ramp_rate` into
every `param_tiers` map (KeyError landmines across the bare-build test fixtures) and pointlessly
sweep it in the stochastic ensemble. The `0.0` `.get` default also shields hand-built param maps in
unit tests. Net: no gating logic, no injected-when-disabled parameter, and the isothermal path is
provably the pre-ramp core.

**Emergent (verified).** A run ramping 14 → 30 °C finishes with residual sugar **between** the
cold-held and hot-held isothermal bounds (`hot < ramp < cold`) — proof the Arrhenius kinetics read
the true time-varying `T`, not a constant, which is the whole point of activating the schedule.

**Tests.** `tests/test_schedule.py` (9) pins the verb-agnostic driver with toy Processes
(isolability, exact per-segment param integration, mutation + ledger, mid-run reconfiguration + tier
travel, day-0 seeding, same-instant ordering, out-of-window rejection). `tests/test_temperature_ramp.py`
(13) pins the temperature path (isothermal no-op, single-knot/flat → no events, collinear → one
segment, slope-change → one event, hold before/after, exact analytic line, scheduled==plain when
isothermal, the end-to-end `run()` multi-knot ramp→hold through the chokepoint, VALIDATED unsampled
rate, the emergent bound). `test_media` expects the always-on `temperature_ramp` in both media.
451 green + 5 benchmark, ruff + mypy clean.

**Deferred → D-36.** Discrete winemaking interventions (the verb registry at the compile boundary:
`add_dap`/`add_so2`/`rack`/`pitch_mlf`), the external-flow ledger's winemaking payoff (a DAP dose's
emergent H₂S-gate response, D-29), reconciling the compile-time MLF disable-gate with a *later*
pitch, and — separately — the stochastic ensemble wrapping `simulate_scheduled` (it wraps `simulate`
today).

## D-36 — Discrete winemaking interventions: the verb registry at the compile boundary

**Status: IMPLEMENTED 2026-07-02** (476 green + 5 benchmark). The winemaking payoff of the D-35
event loop. `Scenario.interventions` — a declarative timeline of verbs (`day`, `action`, `params`
in industry units) — was declared since Milestone 1 but *never consumed*: `compile_scenario` turned
only the temperature schedule into events. This activates it. Built one verb per commit on the
unchanged D-35 driver (`add_dap` → `add_so2` → `rack` → `pitch_mlf`); nothing in `runtime` or `core`
changed — all four verbs are pure vocabulary at the scenario→core compile seam (D-3).

**The registry (`scenario/compile.py`, `_INTERVENTION_VERBS`).** Each action name maps to a compiler
`(Intervention, StateSchema, ParameterSet) → ScheduledEvent`. A verb owns the *meaning*: which
canonical slot a dose lands on, which unit conversion applies, which Processes a pitch enables. The
driver stays verb-agnostic — it just segments-and-restarts and books each state jump as an
`ExternalFlow`. `_compile_interventions` dispatches the timeline, merges the resulting events with
the temperature-ramp events into the single `events` tuple `simulate_scheduled` sorts by time, and
enforces the `_ALLOWED_KEYS` discipline: an unknown verb, a day at/after the run duration, a missing
or unknown param, or a negative dose each raise at the boundary with a scenario-level message. New
verbs are added here and nowhere else. **Isolability:** no interventions ⇒ empty events ⇒ (absent a
ramp) byte-for-byte a plain `simulate`.

**`add_dap` — the headline, a *timing* effect a static dose cannot produce.** Doses diammonium
phosphate by mass (`dap_gpl`) and converts to the assimilable-N jump on the lumped `N` slot via a
new sourced `dap_nitrogen_fraction` (exact (NH₄)₂HPO₄ stoichiometry, 28.014/132.06 = 0.2121 g N/g,
VALIDATED with a zero-width band; new shared `additions.yaml`, the `must_fermentable_fraction`
precedent — a boundary conversion constant, not a magic number). **Phosphate is dropped** (no
phosphorus pool; P is non-limiting) — a scoped omission. The D-29 *static* N→H₂S lever was muted
(~5% span) because N strips to ~0 by day ~1.3 at every dose; a **mid-ferment** DAP dose is
categorically different — it restores N *while sugar (hence the flux the inverse gate
`K_h2s_n/(K_h2s_n+N)` multiplies) is still present*, so the H₂S production **rate drops immediately**
after the dose (verified: ~0.6× the undosed rate just after a day-2 dose) and recovers as the new N
is consumed. A competing effect is present and honest — the extra N feeds growth ⇒ more biomass ⇒
more flux later — but **net cumulative H₂S falls** (gate closure dominates), the realistic direction
(DAP is the standard H₂S-management lever). Emergent, not imposed: the model has no "DAP lowers H₂S"
term. Dose in the **active window** — a post-dryness dose lands where the flux is ~0 and shows
nothing.

**`add_so2` — rides neither elemental ledger.** Doses total SO₂ (`so2_mgl`) onto the conserved
`so2_total` slot (the same slot the initial `so2_total_mgl` addition uses, D-22/D-28); free/bound/
molecular SO₂ are re-derived at the solved pH, so a mid-ferment addition raises the antimicrobial
molecular fraction from that time forward (verified: molecular readout 0 pre-dose → positive
post-dose). SO₂ carries neither carbon nor nitrogen, so the flow perturbs **neither** balance —
both close with no correction term (contrast the DAP nitrogen jump). Raises on a medium without an
`so2_total` slot (beer).

**`rack` — the ledger's removal side.** Draws the wine off a fraction ∈ [0, 1] of its settled lees:
`X_dead` and (when autolysis is opted in, D-34) the cell-wall `debris` (`_LEES_SLOTS`, a single
source of truth). Viable biomass `X` and every dissolved species (sugar, ethanol, YAN, glycerol,
byproducts, acids, SO₂) are **left untouched** — a normal post-AF rack settles dead yeast, and a
concentration model has no volume change on racking, so touching the dissolved pools would be
physically wrong. Books the negative jump as an `ExternalFlow`. Both racked pools carry carbon (and
`X_dead` carries nitrogen), so the removal is a negative term in both ledgers.

**Crown-jewel ledger test (the payoff D-35's external-flow machinery was built for).** A run with
*both* an injection (DAP, +N) and a removal (rack, −C/−N) satisfies the run-wide identity
`final == initial + Σ external_flows` for carbon **and** nitrogen to machine precision — and the
ledger is non-trivial (rack removes carbon, DAP adds nitrogen). The continuous ODE closes exactly
within every segment; the ledger is the correction term across the jumps.

**`pitch_mlf` — the driver's third effect (in-place reconfiguration).** Inoculates *Oenococcus oeni*
mid-run: it both **mutates** `X_mlf` (the bacterial catalyst dose, `pitch_gpl`) and **reconfigures**
the Process set to enable `_MLF_GATED_PROCESSES` — malate→lactate conversion, the citrate
co-metabolism, and the bacterial diacetyl reduction (D-23/D-31). That tuple is now a **single source
of truth** shared with the compile-time disable-gate, so a sequential mid-run pitch is *symmetric*
with an initial co-inoculation and the two cannot drift. `X_mlf` is an inert carbon-/nitrogen-free
catalyst, so the pitch perturbs neither ledger. Because the Processes are enabled only from the
breakpoint, `simulate_scheduled` min-combines the per-segment tier maps (D-35): the malate/lactate/
citrate slots report **speculative for the whole run**, and revert to VALIDATED when unpitched
(disabled ⇒ inert). **Honest scope, verified:** a 22-Brix must (finishing ~107 g/L ethanol) converts
most of its malate under an **early** pitch (day 1: 3.0 → 0.69 g/L) but **stalls** under a *post-AF*
pitch (day 15: 3.0 → 2.97) — past the Luong ethanol wall (~110 g/L) the environmental gate keeps
conversion near zero. The verb makes pitch timing a *scenario* choice; it does not change the
kinetics (malolactic still completes only under co-inoculation / early pitch, D-23).

**Tests.** `tests/test_interventions.py` (25) pins all four verbs: the dose lands on the right slot
and books one flow; the H₂S rate-drop + net-suppression headline; SO₂ perturbs neither ledger and
raises the molecular readout; rack removes only the lees and leaves the wine; the combined DAP+rack
carbon/nitrogen crown-jewel; pitch_mlf enables exactly the gated set, catalyst is ledger-free, early
converts / late stalls, tier travels; the ramp+intervention **merge on one driver** (a multi-knot
temperature schedule *and* a DAP dose — the realistic scenario, and the only test with both sides of
`events` populated); and the vocabulary discipline (unknown verb, out-of-window day, bad params) +
isolability. 476 green + 5 benchmark, ruff + mypy clean.

**Deferred.** ~~The stochastic `simulate_ensemble` still wraps the un-scheduled `simulate`~~ —
**resolved in D-37** (ensemble over a multi-segment schedule). An **MLF-with-growth** consumer
Process composed with a `pitch_mlf` event is now unblocked (the loop exists, the tier travels,
autolysis + the amino-acid ledger refill the pool D-32/D-34) but stays future work. Other addition
verbs (acid/tannin/nutrient blends, chaptalization, cold-stabilisation racking with a viable-`X`
removal fraction) slot into `_INTERVENTION_VERBS` when needed.

## D-37 — Stochastic ensemble over a scheduled run: `simulate_ensemble(events=…)`

**Status: IMPLEMENTED 2026-07-02** (481 green + 5 benchmark, ruff + mypy clean). The last
D-35→D-36 follow-up. `simulate_ensemble` wrapped the *un*-scheduled `simulate`, so an ensemble
could not honour a temperature ramp or a dosing/pitching timeline — it shared the multi-segment
footgun D-35's `CompiledScenario.run()` was created to avoid. This routes the wrapper through
`simulate_scheduled` and adds an `events` parameter (default `()`), plus a
`CompiledScenario.run_ensemble(**kwargs)` that threads the compiled `events` — the stochastic
sibling of `run()`. With `events=()` it is byte-for-byte the previous ensemble (a no-event
`simulate_scheduled` is a single `simulate` segment), so every pre-existing ensemble test stays
green unchanged; the nominal run routes through the same path so its min-combined `tier_map` is
consistent with the members.

Three interactions the naive "just call `simulate_scheduled`" would get wrong, each handled:

* **Process-set isolation (the load-bearing one).** A `reconfigure` event (`pitch_mlf`) mutates the
  shared `ProcessSet._enabled` in place and is deliberately *not* self-restoring — a mid-run pitch
  persists for the rest of *that* run (D-35), a contract `test_interventions` pins on `cs.run()`.
  But an ensemble replays the schedule N times over the *same* set, so member i's `enable` would
  leak into member i+1's pre-pitch segments. Fix: a new public `ProcessSet.enabled_snapshot()` /
  `restore_enabled()` primitive; the ensemble captures the pristine state once and resets before
  every member (nominal included), and leaves the set pristine when done — a *batch* is side-effect-
  free on the set, distinct from a single *run* whose enable persists. Isolation lives in the
  wrapper, not in `simulate_scheduled` (moving it there would break the persist-within-a-run
  contract). *Subtlety, tested honestly:* for the **current** `pitch_mlf` verb the leak is
  numerically **inert** — the enabled MLF Processes are gated by the `X_mlf` catalyst, still 0 until
  the pitch *mutation*, so a leaked pre-pitch enable contributes zero flux and per-member conservation
  cannot catch it. The reset is therefore **defensive** for a future catalyst-free `reconfigure`, and
  its guard test uses a synthetic *ungated* enable so the leak is observable (byte-for-byte against
  independent fresh-set runs; the test fails if the per-member reset is removed).
* **Sampling scope must span the schedule.** `_resolve_sample_names` scoped to the *`t0`*-active
  reads, but a `pitch_mlf` enables the malolactic Processes only from the breakpoint — their
  kinetics (`k_mlf`, …) are disabled at `t0` and would be silently dropped from the sampled set,
  under-sampling exactly the parameters the pitched back half depends on. New `_schedule_reads`
  unions the active reads across every `reconfigure` in the schedule (replaying them onto a snapshot,
  then restoring). Over-covering is safe (sampling a param no active Process reads is a documented
  no-op, D-24); under-covering silently narrows the reported spread, so the union is strictly right.
* **The external-flow ledger is member-dependent.** DAP/SO₂/pitch inject fixed masses, but `rack`
  removes a *fraction of the settled lees*, whose mass at rack time depends on each member's sampled
  death/growth kinetics — so every member's removal `delta` differs. `Ensemble` gained `member_flows`
  (per member) + `segment_bounds` (scenario-fixed, stored once) + `nominal_flows`, and
  `member_trajectory(i)` / `nominal_trajectory()` now return a `ScheduledTrajectory` carrying them, so
  the across-jumps identity `final == initial + Σ flows` is auditable *per draw* — the crown-jewel
  D-36 conservation guard, extended to the whole ensemble rather than just the nominal.

**Type seam.** Because `member_trajectory` now returns a `ScheduledTrajectory`, the kinetics-agnostic
conservation helpers (`assert_conserved`, `max_drift`, `assert_nonnegative`) were re-typed against a
new structural `TrajectoryLike` Protocol (schema/`t`/`y`/`series`, read-only) — they never needed the
scheduling extras, only the state grid, so both trajectory types satisfy it and no call site changed.

**Caveat (documented, not built for).** A parameter that is *both* sampled and overwritten by an
event's `param_update` would use its sampled value pre-event and the fixed compile-time value after.
No current verb hits this: the only `param_update` payload is `temperature_ramp_rate`, which declares
no `reads` and is VALIDATED, so it is never sampled (D-35). A future forcing parameter that is both
uncertain and event-driven would need pinning via `exclude`.

**Tests.** `tests/test_ensemble.py` gained 5: un-scheduled isolability (empty ledger, single
`[t0, t_end]` segment, `member_/nominal_trajectory` still audit); the DAP+rack crown-jewel identity
*per member* with member-dependent rack removal (`np.std(rack_carbon) > 0`); schedule-union scope
(`k_mlf` sampled under a `pitch_mlf` despite being absent from the `t0` reads); process-set post-run
pristine + tier travel (`malic`/`lactic` speculative run-wide); and the discriminating **isolation**
guard (an *ungated* toy `reconfigure` so a leaked enable is numerically visible — each member equals
an independent fresh-set run byte-for-byte; verified to fail when the per-member reset is removed).
481 green + 5 benchmark, ruff + mypy clean.

## D-38 — MLF-growth: `X_mlf` becomes dynamic bacterial biomass (resolves the D-23 deferral)

**Status: IMPLEMENTED 2026-07-02** (496 green + 5 benchmark, ruff + mypy clean). The long-deferred
MLF-with-growth beat D-23 scoped out and D-32/D-33/D-34 built the prerequisites for. `X_mlf` was a
dosed-but-**inert** catalyst (constant, carbon-/nitrogen-free) that merely *scaled* the
`MalolacticConversion` rate; D-23 promised the growth beat would be "a clean add-a-Process
extension." This delivers exactly that: a new `MalolacticGrowth` Process makes `X_mlf` dynamic, so —
because conversion is **linear in `X_mlf`** — deacidification now *accelerates autocatalytically* as
the bacteria multiply. No refactor of the conversion kinetics.

**The growth law.** `dX_mlf/dt = μ_max_mlf · X_mlf · aa/(K_aa_mlf+aa) · S/(K_s+S) ·
g_pH·g_EtOH·g_SO₂·γ(T)` — Michaelis–Menten in the amino-acid fuel *and* in sugar (the fermentable
energy O. oeni co-metabolises), scaled by the **same** `malolactic_environmental_gate` conversion
uses. Reusing that gate is load-bearing: the Luong ethanol wall makes co-inoculation the *dominant*
MLF-growth mode **emergently** — a post-AF pitch into a high-ABV must lands past the O. oeni ethanol
tolerance so γ·g_EtOH ≈ 0 and bacteria cannot build up, while a normal-ABV sequential MLF (g_EtOH
small but nonzero) still grows. New speculative params `mu_max_mlf` (0.05/h) and `K_aa_mlf`
(0.05 g/L); it reuses `K_s` and the biomass fractions.

**Conservation — nitrogen-anchored, carbon shortfall from sugar (the anchoring fork).** New
bacterial biomass needs `f_N·dX_mlf` nitrogen and `f_C·dX_mlf` carbon. All the nitrogen comes from
the `amino_acids` pool (arginine), consuming `ρ = f_N·dX_mlf/y_N`; that arginine carries only `ρ·y_C`
carbon — *less* than the biomass needs, because arginine (mass C:N ≈ 1.29) is far more N-rich than
biomass (C:N ≈ 4–11) — so the **shortfall** `f_C·dX_mlf − ρ·y_C = dX_mlf·(f_C − f_N·y_C/y_N)` is
drawn from sugar. This is the mirror of yeast `GrowthNitrogenLimited` (N from a N-pool, C from
sugar) and the inverse of D-34 autolysis (which routes the *excess* carbon to debris). The shortfall
coefficient is **structurally positive** across Coleman's whole `f_N` range → no clamp, no C⁰ kink.
Touches `(X_mlf, amino_acids, S)` — notably **not** `N`.

- **Fork decided (owner away — advisor-recommended, higher-fidelity branch chosen; owner may
  revisit).** The alternative was **C-anchored**: consume arginine for the biomass *carbon* and
  deaminate the over-supplied nitrogen to ammonium `N` (the real O. oeni arginine-deiminase pathway).
  It conserves just as cleanly and works even after sugar is gone. It was **rejected** because the
  demo regime (co-inoculation) has abundant sugar (245→~45 g/L across the growth window, since the
  ethanol wall isn't crossed until ~day 4), so the "no sugar at MLF" premise that would have favoured
  it does not hold; and it carries two fictions the N-anchored branch avoids — booking *all* biomass
  carbon as arginine-derived (ADI actually *excretes* that carbon as ornithine/CO₂, it doesn't build
  biomass), and dumping a large artificial ammonium surplus. Gate-on-sugar then self-limits growth to
  the window it is physical in. This is a discuss-disagreements fork; the anchoring can be swapped
  later if the owner prefers the ADI reading.

**`X_mlf` promoted from inert catalyst to real biomass.** For the growth to conserve, `X_mlf` is now
**weighted in `total_carbon`/`total_nitrogen`** at the biomass fractions (bacterial ≈ yeast elemental
composition — a documented v1 simplification; the *same* fractions the growth stoichiometry draws
against, so closure is exact). Consequence, superseding the v1 "dosing X_mlf leaves total_carbon
byte-for-byte" claim: a co-inoculation dose / `pitch_mlf` flow now carries bacterial-biomass
carbon/nitrogen (booked on the D-36 external-flow ledger, which still closes; `test_interventions`
updated accordingly). On a conversion-only run `X_mlf` is constant, so the added ledger term is a
constant offset that still drifts to zero — every existing MLF/ensemble closure test stays green.

**Gating — its own tuple, keyed on amino acids alone, NOT the pitch (advisor-corrected).** The
compile seam disables `MalolacticGrowth` when `amino_acids_gpl ≤ 0` — the *same* gate as the D-32
swap / D-33 re-route, keyed on the feature (amino-acid-fed bacterial growth), not on the pitch. This
alone prevents the tier-isolability regression: every existing D-23/D-31 test pitches O. oeni but
doses *no* amino acids, so growth stays disabled and never drags the `amino_acids`/`S`/`X_mlf` tier
via `tier_of`. It is kept in its own tuple (not `_MLF_GATED_PROCESSES`) precisely because that gate
differs from conversion's. It is deliberately **not** additionally gated on the pitch: "bacteria
present" is runtime state the Process's own `X_mlf ≤ 0` guard handles, and whether post-pitch
bacteria then *grow* is left to the emergent environmental gate — mirroring how `MalolacticConversion`
trusts its ethanol gate rather than a compile rule. So co-inoculation dominance is **emergent** (a
high-ABV post-AF pitch is ethanol-arrested; a normal-ABV sequential MLF can still grow), not a
hard-coded co-inoc-only rule. *(An earlier draft gated on `pitch AND aa` and attributed co-inoc-only
to the ethanol wall as a compile fact — corrected here: that conflated the tier feature with runtime
state and mis-stated the physics.)* No existing test doses amino acids under an MLF pitch, so no
trajectory shifted.

**Scope / caveats (owned).** (1) Bacterial biomass composition = yeast's (incl. the Coleman `f_N`
override) — a lump; a dedicated bacterial `f_N`/`f_C` is deferred. (2) v1 is **growth-only**: no
bacterial death/decay (so the D-31 "SO₂/rack locks in diacetyl by killing the bacteria" case still
needs a death Process + the deferred rack-removes-`X_mlf`). (3) All-biomass-carbon-from-arginine-
or-sugar is a stand-in, not a claim about O. oeni's carbon metabolism.

**Tests.** `tests/test_mlf_growth.py` (15): the fail-first acceptance (same co-inoc+aa scenario,
growth on vs the growth Process disabled — the fixed-`X_mlf` control — day-3 malate halved,
`X_mlf` multiplied several-fold, gap vanishes without the Process); carbon+nitrogen closure over a
growing run; the `(X_mlf, amino_acids, S)` `touches` contract; derivative-level stoichiometry
closure; the aa-keyed gate matrix + the emergent mid-run-pitch case (early pitch grows, late
post-AF pitch is ethanol-arrested); the no-catalyst/no-fuel/no-sugar guards; never-creates-sugar
(sugar carbon drawn < biomass carbon built); the ethanol-wall arrest; and the speculative tier
capping `X_mlf` (discriminated on an aa-dosed *unpitched* run where growth is the only enabled
`X_mlf` toucher). `test_interventions` / `test_media` updated for the promotion + the new Process.

## D-39 — MLF death: `X_mlf` dies under SO₂ (`MalolacticDeath`) + rack removes it; the MLF arc closes

**Status: IMPLEMENTED 2026-07-02** (506 green, ruff + mypy clean) — two commits, one beat each:
commit 1 the `MalolacticDeath` Process, commit 2 the `rack`-removes-`X_mlf` extension. The counterpart to the
D-38 growth beat that completes the MLF arc (D-23 → D-31 → D-38 → D-39): a new `MalolacticDeath`
Process moves viable `X_mlf` into a new non-viable `X_mlf_dead` pool, so bacterial biomass now
*declines* and the *O. oeni* activities that scale with `X_mlf` — malate conversion, citrate →
diacetyl, and above all `OenococcusDiacetylReduction` — wind down as the bacteria die. This is the
mechanism the D-31 reducer flagged as deferred: **SO₂ (or a rack) removes the bacteria that clear
diacetyl on the lees, so it is locked in.** (Rack-removes-`X_mlf` is commit 2 of this decision.)

**The death law — SO₂-driven, Arrhenius temperature.**
`r_death = k_death_mlf · X_mlf · (1 − g_SO₂) · arrhenius(T, E_a_death_mlf, T_ref)` with
`g_SO₂ = exp(−[SO₂]_molecular / molecular_so2_inhib_mlf)` — **the same `g_SO₂` the conversion gate
uses** (D-22 antimicrobial readout, partitioned at the solved pH). Death is **exactly 0 without
SO₂** and rises toward its Arrhenius ceiling as molecular SO₂ accumulates. Temperature enters via
its **own Arrhenius factor** (warm accelerates the kill, cold slows it — the autolysis shape),
**not** the cardinal γ(T): γ(T) → 0 in the cold, which would spuriously make cold *kill*, whereas
cold in fact *preserves* bacteria. To supply that driver the shared MLF gate was split (no behaviour
change) into `malolactic_toxicity_gate` (pH·ethanol·SO₂) × `cardinal_temperature_factor` (γ(T)), with
the multiplication grouped exactly as before so the three growth/conversion consumers are
byte-for-byte unchanged (`test_environmental_gate_is_toxicity_times_gamma` pins the identity).

- **The crux — driver form, decided on empirical evidence (advisor-reconciled, owner's fidelity bar).**
  The first draft drove death by **`1 − toxicity`** (the full pH·ethanol·SO₂ gate), on the theory that
  accumulating ethanol would supply a natural post-AF die-off "for free." A probe **killed that form**:
  the Luong ethanol wall already drives `1 − toxicity` to **~0.92 at ordinary post-AF ethanol (~75
  g/L, no SO₂)**, so death was near-maximal *from ethanol alone* — *O. oeni* died in ~1 week, when in
  reality it persists for weeks-to-months in dry wine and is cleared deliberately by SO₂/racking. No
  power transform `(1 − tox)^p` rescues it: 0.92 cannot be mapped both ~0 (slow baseline) *and* kept
  clearly below the SO₂-elevated 0.996 (the SO₂:baseline ratio maxed at ~1.4× even at p=4). Ethanol's
  wall is a "can't grow" signal, not a "dying" one; coupling death to it was the bug. **Fix: a driver
  with no ethanol term — molecular SO₂ only.** This *decouples* `k_death_mlf` from the (unsulfited)
  early-pitch conversion test, so k was **re-tuned up** from the artifact 0.02/h to **0.05/h** — a
  full-SO₂-kill half-life ~14 h, so a stabilizing dose (~40 mg/L free ⇒ ~0.8 mg/L molecular) crashes
  the population ~90 % in ~2 d, verified directly. Co-inoc-vs-post-AF dominance now rests entirely on
  the *growth* gate's `g_EtOH`, where it belongs; a high-ABV post-AF pitch simply sits **inert**.

**v1 tradeoff (owned, not hidden; → RESOLVED in D-41).** Without SO₂ bacteria **never die** in v1 —
they persist and keep clearing diacetyl on the lees (the honest D-31 "leave on lees cleans up" case).
The slow ethanol/age decline of *O. oeni* over weeks-to-months was then **deferred to v2** — and is
now landed as the separate `MalolacticSenescence` baseline mortality (**D-41**), leaving this SO₂ kill
byte-for-byte unchanged. Less realistic than a slow decline, far more realistic than the 1-week ethanol
wipeout — and it makes the D-31 SO₂/rack lever **unconfounded**: only a deliberate winemaking action
removes viable bacteria. The kill-scale reuses `molecular_so2_inhib_mlf` (arrest-scale = kill-scale,
bacteriostatic ≈ bacteriocidal); a separate `molecular_so2_death_scale` is the v2 refinement.

**Conservation — a carbon/nitrogen-neutral transfer (the D-13 pattern).** Since D-38 both `X_mlf` and
`X_mlf_dead` are weighted in `total_carbon`/`total_nitrogen` at the *same* biomass fractions, so the
`X_mlf → X_mlf_dead` move (`d[X_mlf] = −r`, `d[X_mlf_dead] = +r`) is C- and N-neutral **by
construction** — the yeast `EthanolInactivation` `X → X_dead` precedent (D-13). No new conservation
code, no sugar draw; touches only `(X_mlf, X_mlf_dead)`. The new `X_mlf_dead` slot is the tenth wine
slot (schema size 26 → 27); `conservation.py` weights it at the biomass fractions guarded on presence.

**Gating.** `MalolacticDeath` is **pitch-gated** (enabled with the other `_MLF_PROCESSES` /
`_MLF_GATED_PROCESSES` when `mlf_pitch_gpl > 0`), NOT amino-acid-gated like growth — bacteria die
whether or not they were growing. An `so2_total ≤ 0` guard before the pH `brentq` makes an unsulfited
pitched run pay no solve and contribute byte-for-byte zero (mirrors the `total_so2 > 0` shortcut in
the toxicity gate). On a pitched run `X_mlf`/`X_mlf_dead` report **speculative** (honest: a population
that can be sulfited has a speculative trajectory). New speculative params `k_death_mlf` (0.05/h),
`E_a_death_mlf` (60 kJ/mol = `E_a_autolysis`); reuses `T_ref` + `molecular_so2_inhib_mlf`.

**Tests (commit 1, +8).** In `test_malolactic.py`: the SUPERSEDING integration test (no-SO₂ pitched
run is byte-for-byte inert — death exactly 0 — then a mid-run `add_so2` crashes `X_mlf` monotonically
to <10 %); RHS-level death-is-zero-without-SO₂, the neutral `d[X_mlf] = −d[X_mlf_dead]` transfer, the
`(X_mlf, X_mlf_dead)` `touches` contract, more-SO₂-kills-faster, the load-bearing
**cold-preserves-via-Arrhenius-not-γ(T)** case (dying below `T_min_mlf` where γ(T)=0, warm faster),
carbon+nitrogen closure over a death-active run, the speculative tier, and the gate-split identity.
`test_media` updated for the tenth slot + the new pitch-gated Process.

**Commit 2 — `rack` removes viable `X_mlf` + settled `X_mlf_dead` (the D-31 lever's physical half).**
Both *O. oeni* pools join `_LEES_SLOTS`, so a `rack` draws them off with the lees — the physical twin
of the SO₂ kill: racking early strands diacetyl (the deferred D-31 "rack ⇒ locked in" case). This is
a deliberate **asymmetry with yeast**, owned in the docstring: a rack leaves viable *yeast* `X`
untouched (it ferments in *suspension*, so racking gross lees leaves it working), but *O. oeni*
carries out MLF *on the lees* and goes with them. Both bacterial pools carry biomass C/N (weighted
since D-38), so — like `X_dead` — their removal books a negative C/N `ExternalFlow`; the run-wide
`final == initial + Σ flows` identity closes to machine precision for both elements (SO₂ carries
neither, so only the rack moves the ledger). No kinetics change — a `_LEES_SLOTS` + docstring edit.
The single-run "rack strands diacetyl" demo is confounded exactly as the death case is (removing
bacteria drops both the diacetyl sink *and* its citrate source), so it is validated on the `X_mlf`/
`X_mlf_dead` removal + conservation directly, not on a diacetyl curve. **+2 tests** (`test_
interventions`, 506 green): rack removes both pools while leaving viable yeast + dissolved species
untouched; and C/N closure across a dose-then-rack MLF run. **The MLF arc (D-23 → D-31 → D-38 → D-39)
is complete.**

## D-40 — Brettanomyces volatile phenols: the mixed-culture beat that closes Milestone 2

**Status: pt1 + pt2 + pt3 IMPLEMENTED 2026-07-02** (ruff + mypy clean; full suite green — see commits). The
last unchecked M2 physics beat ("Mixed cultures / Brett / sour consortium"). *Brettanomyces
bruxellensis* is the canonical wine spoilage yeast: it decarboxylates grape-must **hydroxycinnamic
acids** (p-coumaric, ferulic) to **vinylphenols**, then reduces those to the **ethylphenols**
(4-ethylphenol "barnyard", 4-ethylguaiacol "clove") that define Brett character. Built as a
multi-commit arc mirroring the MLF arc (conversion → growth → death): **pt1 = the phenol pathway with
a dosed catalyst**, **pt2 = `BrettGrowth`** (dynamic `X_brett`), **pt3 = `BrettDeath`** (the SO₂
kill); pt4 the POF+ yeast opt-in + emergent reservoir test — to follow.

**Two owner forks (decided by the user, pros/cons presented).** (1) *Pathway fidelity* → **3-pool +
POF+ yeast**: the `vinylphenols` intermediate earns its own state slot because it carries *emergent*
behaviour — a POF+ *S. cerevisiae* fills a shared reservoir it cannot clear (it has the decarboxylase
but not the reductase), and only Brett drains it, so "no Brett ⇒ vinylphenol strands" emerges (the
α-acetolactate-reservoir parallel, D-26/D-31). (2) *Phenol scope* → **lumped 4-EP + 4-EG**: one
`ethylphenols` pool from a lumped hydroxycinnamic precursor (booked as p-coumaric / 4-vinylphenol /
4-ethylphenol representative species). The two compose coherently — depth on the pathway (where
behaviour lives), lumping on the readout (the same mechanism twice with different sensory labels).

- **The advisor's blind-spot fix (fidelity, not preference).** The initial fork framing gated
  decarboxylation to *yeast only*, which would produce **nothing** for the canonical case the beat is
  named for: a **POF-negative wine spoiled by Brett alone** (yeast makes no vinylphenol → nothing to
  reduce). Reality: **Brett carries BOTH enzymes** — that is *why* it spoils normal wine unaided. So
  the Process set is the *union*: **Brett gets its own decarboxylase** (`BrettDecarboxylation`,
  `X_brett`-gated) *and* its reductase (`BrettVinylphenolReduction`), and the POF+ *yeast*
  decarboxylase becomes a separate **opt-in strain** Process (pt4, default OFF) — not gated on
  precursor presence (a POF- yeast in hydroxycinnamic-rich must must make no vinylphenol). The
  headline acceptance test is therefore the *canonical* case, not the POF+ reservoir.

**Carbon closes on the existing ledger — no new conservation code.** `BrettDecarboxylation`:
p-coumaric (C9) → vinylphenol (C8) + CO2 (C1), carbon-closing mole-for-mole (9 = 8 + 1, the malic →
lactic + CO2 idiom, D-23). `BrettVinylphenolReduction`: vinylphenol (C8) → ethylphenol (C8), a
mole-for-mole C8 → C8 transfer between two weighted pools (the diacetyl → butanediol idiom, D-26).
`total_carbon` weights all three phenol pools at their representative species, so the Processes touch
only `hydroxycinnamics`/`vinylphenols`/`ethylphenols`/`CO2` and add nothing to the harness (verified
closing to machine precision through the full precursor → intermediate → product chain).

**The Brett environmental gate — SO₂ and temperature only (the advisor's explicit warning).** Unlike
*O. oeni*, Brett is markedly **acid-tolerant** (spoils low-pH wine) and **ethanol-tolerant** (a
full-strength-wine barrel spoiler), so copying the MLF gate's pH logistic + Luong ethanol wall would
spuriously arrest Brett exactly where it thrives. So `gate = g_SO₂ · γ(T)` — **no pH, no ethanol
term**: molecular SO₂ (the D-22 antimicrobial readout) is the winemaker's lever, and a **cardinal
temperature optimum warmer than *O. oeni*'s** (`T_opt_brett` 32 °C vs MLF's 23 °C — Brett is a
warm-tolerant spoiler). The ethanol tolerance is asserted at the integration level:
`test_pitch_brett_post_af_at_high_ethanol` pitches Brett into a *finished* ~14 % ABV wine and confirms
4-EP still rises — the property that would silently die if anyone re-added an ethanol wall.

**Isolability + the compile seam (the MLF pattern).** `X_brett` is a constant, **carbon-free** dosed
catalyst in pt1 (weighted as real biomass only when `BrettGrowth` lands, pt2 — the exact `X_mlf`
D-23 → D-38 path). The Processes are wired into the wine medium but return zero before any pH work
when `X_brett ≤ 0` or the substrate is absent; the compile seam **disables** them unless Brett is
pitched (`brett_pitch_gpl` co-inoculation, or a mid-run `pitch_brett` intervention re-enabling the
same `_BRETT_GATED_PROCESSES` at its breakpoint), so an unpitched wine run is byte-for-byte the
validated core and the phenol slots keep their **VALIDATED** tier (`tier_of` counts enabled, not
nonzero, Processes). Both `X_brett`/`X_brett_dead` join `_LEES_SLOTS`, so racking draws Brett off the
lees (the spoilage twin of the SO₂ kill). Wine-only (beer has no phenol slots).

**Headline acceptance gate — a control-difference (parallels `test_headline_mlf_...`).** A POF-
wine + dosed hydroxycinnamics accumulates `ethylphenols` **only when Brett is pitched** (the no-Brett
control stays exactly 0); an SO₂ dose suppresses 4-EP >10× (metabolic arrest), and a rack removes
`X_brett` and halts production at the breakpoint. **+11 tests** (`test_brett.py`): headline, post-AF
`pitch_brett` verb + ethanol tolerance, SO₂/rack levers, carbon closure, per-Process
stoichiometry/`touches`, guards, unpitched tier isolability, the warm temperature optimum, and the
`speculative` tier. Two `test_media.py` composition assertions updated for the 5 new wine slots + 2
Brett Processes. All params `speculative` (author estimates; no per-catalyst kinetic model of this
flux form is sourced — Brett phenols are reported as bulk mg/L end-yields).

**pt2 — `BrettGrowth`: `X_brett` becomes dynamic (IMPLEMENTED 2026-07-02).** The Brett twin of
`MalolacticGrowth`, with one load-bearing difference: **Brett grows on ETHANOL, not sugar**, so it
builds up in a *dry, finished* wine — its real post-AF/barrel niche. Because the decarboxylase and
reductase are linear in `X_brett`, a growing population makes the volatile-phenol spoilage
**accelerate autocatalytically** over the months a barrel sits — the "it gets worse the longer you
leave it" dynamic a constant catalyst cannot produce. `dX_brett/dt = μ_max_brett · X_brett ·
aa/(K_aa_brett+aa) · E/(K_E_brett+E) · g_SO₂·γ(T) · (1 − X_brett/K)`.

- **Owner fork — carbon source → ETHANOL-drawn (decided by the user, pros/cons presented).** New
  biomass is nitrogen-anchored on the `amino_acids` pool (D-32, autolysis-refilled D-34), consuming
  `ρ = f_N·dX_brett/y_N` of arginine; the carbon **shortfall** `f_C·dX_brett − ρ·y_C` (arginine is
  N-rich, so it under-supplies carbon) is drawn from **ethanol `E`**, not sugar. That is the
  mechanistic reason Brett thrives where the wine is *dry*. Both ledgers close exactly; touches
  `(X_brett, amino_acids, E)` — **not** `S` (Brett skips sugar) and **not** `N` (no ammonium release,
  the D-38 anchoring choice). v1 models only the biomass-assimilation branch; the acetic-acid
  overflow (Brett's real ethanol-oxidation "volatile acidity" product) is a deferred pool, so the
  ethanol drawdown here is a lower bound on true consumption.

- **The carrying-capacity brake — required, because Brett has no self-arrest (the numeric crux).**
  `MalolacticGrowth` is self-limiting (its sugar Monod vanishes as sugar is consumed *and* its gate
  carries an ethanol wall). Brett deliberately has **neither** (dry-wine, ethanol-tolerant niche), so
  amino-acid Monod alone is *not* a ceiling: an autolysis-refilled aa pool would grow `X_brett`
  exponentially without bound. So `BrettGrowth` carries an intrinsic **logistic carrying capacity**
  `(1 − X_brett/K)` (`brett_carrying_capacity`), the same lumped form as the D-30 yeast
  `BiomassCarryingCapacity` (real Brett saturates at a finite cell density) — but *intrinsic and
  always-on*, not the opt-in isolable modifier D-30 is. Bounding `X_brett` small keeps the
  amino-acid draw rate small, so the pool depletes *smoothly* to a positive residual rather than
  overshooting negative.

- **The advisor-caught BDF blow-up + the fix (fidelity, not preference).** The first `BrettGrowth`
  drove `X_brett` → **23 g/L** and `amino_acids` → **−4.5 g/L** under the default **BDF** solver —
  yet **RK45 and LSODA both gave the correct bounded answer** (`X_brett` → ~0.1, aa ≥ 0). The RHS was
  *correct*; BDF was mis-integrating. Root cause (advisor's diagnosis): every hard guard must be
  shadowed by a *smooth* factor that reaches zero first — `aa` is shadowed by its Monod, the brake by
  `(1−X/K)`, but the **`E ≤ 0` guard had no shadow**, so `∂f_X/∂E` was a step at `E = 0`. BDF's
  finite-difference Jacobian straddled that step as ethanol rose through zero during primary AF,
  corrupting the Newton solve into an autocatalytic blow-up (the aa negativity is a *consequence* of
  the X blow-up, not an independent failure); RK45/LSODA build no Jacobian and never saw it. **Fix:
  the ethanol Monod `E/(K_E_brett+E)` is that missing smooth shadow** — and it is *also* physically
  right (Brett grows *on* ethanol, so growth scales with ethanol availability: ≈0 in an unfermented
  must, ≈1 in a finished wine, `K_E_brett` = 2 g/L kept small so it is ≈1 across the working range).
  So growth is now gated by ethanol availability, refining the pt1-era "amino-acid fuel + SO₂/temp"
  story. The regression is pinned **under BDF specifically** (`test_growth_bounded_..._under_bdf`
  asserts `assert_nonnegative` at `atol=1e-8` — the assertion that *caught* the bug — plus a
  BDF-vs-RK45-vs-LSODA agreement test that directly encodes "all three solvers agree").

- **Isolability (stricter gate than pt1).** `BrettGrowth` is wired into the wine medium but disabled
  at the compile seam unless a scenario **both** pitches Brett **and** doses amino acids (a stricter
  gate than the pt1 phenol Processes, so it is a separate tuple — avoids dragging the `amino_acids`/`E`
  tier onto pitched-but-not-aa-dosed runs, mirroring `MalolacticGrowth` vs `MalolacticConversion`).
  `X_brett` promotes from the pt1 carbon-free constant catalyst to weighted biomass with **no verb
  change** — the pitch/rack `ExternalFlow` auto-books its conservation flow (the exact `X_mlf`
  D-23 → D-38 path). **+10 tests** (`test_brett.py`): autocatalytic acceleration headline, the two
  BDF regressions, carbon+nitrogen closure, ethanol-drawn `touches`, the ethanol-availability Monod
  (+ its smoothness), the carrying-capacity brake, growth guards, the aa-gated compile-seam
  isolability, and the `speculative` tier. New params `mu_max_brett`, `K_aa_brett`, `K_E_brett`,
  `brett_carrying_capacity` — all `speculative` (Brett is a characteristically slow grower; no
  per-organism kinetic values sourced).

**pt3 — `BrettDeath`: the SO₂-driven kill (IMPLEMENTED 2026-07-02).** Completes the Brett arc
(pt1 pathway → pt2 growth → pt3 death), the twin of `MalolacticDeath` (D-39). It moves viable
`X_brett` into the non-viable `X_brett_dead` pool under molecular SO₂, so the spoilage population
*declines* when the wine is sulfited and the phenol activities that scale with `X_brett`
(decarboxylase + reductase) wind down. `r_death = k_death_brett · X_brett · (1 − g_SO₂) ·
arrhenius(T, E_a_death_brett, T_ref)`, `g_SO₂ = exp(−[SO₂]_molecular / molecular_so2_inhib_brett)`.

- **SO₂ alone is the *natural* driver for Brett (contrast the D-39 crux).** `MalolacticDeath` had to
  *drop* an ethanol/pH toxicity driver because *O. oeni*'s Luong ethanol wall spuriously made
  bacteria "die" from ordinary post-AF ethanol. Brett has **no such wall** — its gate
  (`brett_environmental_gate`) carries no ethanol or pH term at all, because Brett is ethanol- and
  acid-tolerant — so "molecular SO₂ alone kills Brett" is not a confounder-correction but the
  *directly correct* physics: the winemaker's ~0.5–0.8 mg/L molecular-SO₂ Brett-control target is the
  real-world expression of this term. Without SO₂ (or a rack) Brett persists indefinitely in v1 — an
  honest reflection of how tenacious a barrel Brett infection is; a slow benign-environment
  senescence is a deferred v2 refinement.

- **Arrhenius temperature, not the cardinal γ(T) (the load-bearing D-39 choice reused).** Death
  carries its own Arrhenius factor (warm accelerates the kill, cold slows it toward dormancy), **not**
  the metabolic gate's cardinal γ(T): γ(T) → 0 in the *cold*, which would make cold *kill* Brett,
  whereas cold in fact **preserves** it — part of why Brett is so hard to eradicate from a cool
  cellar, and why it is cleared by SO₂, not by chilling. `test_cold_preserves_brett_via_arrhenius_\
  not_gamma` pins this: below `T_min_brett` (where γ(T) = 0) death is still > 0 and rises with warmth.

- **Conservation — a carbon/nitrogen-neutral transfer (D-13), no new ledger code.** Since pt2 both
  `X_brett` and `X_brett_dead` are weighted in `total_carbon`/`total_nitrogen` at the same biomass
  fractions, so `d[X_brett] = −r`, `d[X_brett_dead] = +r` is neutral in both ledgers by construction
  (the yeast `X → X_dead` and bacterial `X_mlf → X_mlf_dead` precedent). Touches only
  `(X_brett, X_brett_dead)`.

- **Isolability + wiring.** Guards return zero *before* the pH `brentq` when `X_brett ≤ 0` or
  `so2_total ≤ 0` (the SO₂ guard is exact — death is identically 0 without SO₂), so a
  pitched-but-unsulfited run is byte-for-byte inert. `BrettDeath` is **pitch-gated** (in
  `_BRETT_PROCESSES`/`_BRETT_GATED_PROCESSES`), not amino-acid-gated like `BrettGrowth` — Brett dies
  whether or not it was growing, exactly as `MalolacticDeath` sits in `_MLF_PROCESSES` rather than the
  growth tuple. Racking already removes both pools (pt1 `_LEES_SLOTS`), so the physical twin of the
  SO₂ kill needed no new work. `reads` lists `molecular_so2_inhib_brett` explicitly (**not** the
  `_BRETT_GATE_READS` cardinals) — death uses Arrhenius, so it must not pull in `T_*_brett`.
  Consequence: on any *pitched* run `X_brett`/`X_brett_dead` report **speculative** (an enabled
  Process touches them) — honest, matching MLF; no test asserted them VALIDATED on a pitched run.

- **Headline + tests.** `test_so2_crashes_growing_brett_population` is the arc payoff and the
  advisor-sharpened discriminator: with amino acids dosed `X_brett` grows autocatalytically, then a
  mid-run SO₂ addition **kills** it — the unambiguous death signal (distinct from the growth gate's
  mere arrest) is that `X_brett_dead` *accumulates* **and** `X_brett` falls below its value at the
  dose, while the un-sulfited control keeps growing; ethylphenols end lower. **+8 tests** (`test_\
  brett.py`): the headline, the MLF-death-mirrored RHS suite (zero-without-SO₂, neutral transfer,
  `touches`, more-SO₂-kills-faster, cold-preserves-via-Arrhenius), integration-level carbon+nitrogen
  closure, and the `speculative` tier. `test_media.py` `BRETT_PROCESSES` gains `brett_death`. New
  params `k_death_brett` (0.03/h, below `k_death_mlf` — Brett is more SO₂-tolerant than *O. oeni*) and
  `E_a_death_brett` (60 kJ/mol, = `E_a_death_mlf`/`E_a_autolysis`), both `speculative` (no per-catalyst
  Brett mortality law is sourced; direction — SO₂ kills, cold preserves — is sound). **535 green** +
  5 benchmark, ruff + mypy clean.

**pt4 — `YeastPOFDecarboxylation`: the POF+ yeast opt-in + emergent reservoir (IMPLEMENTED 2026-07-06).**
Closes the Brett arc (and Milestone 2's last physics beat). A **POF+** (phenolic-off-flavour-positive)
primary *S. cerevisiae* strain carries the cinnamate decarboxylase — the *same* reaction as
`BrettDecarboxylation`, drawing must `hydroxycinnamics` into `vinylphenols` + CO2 (p-coumaric C9 →
vinylphenol C8 + CO2 C1, carbon-closing 9 = 8 + 1) — but **not** the reductase, so during AF it fills
the shared `vinylphenols` reservoir it cannot drain. With no Brett the vinylphenols **strand**
(`ethylphenols` stays 0); a Brett contamination arriving later gets a **head start** on the pre-filled
reservoir. This is the emergent yeast/Brett coupling the 3-pool design was chosen for (the
α-acetolactate-reservoir parallel, D-26/D-31), and the advisor's blind-spot fix realised: the union of
enzymes, with the POF+ yeast decarboxylase a separate opt-in.

- **Fork 1 (Process vs strain-flag) → separate opt-in Process** (settled at D-40; a strain flag baked
  into the always-on primary set would break byte-for-byte core isolability, prime directive #3). New
  `YeastPOFDecarboxylation` in `brett.py`, its own wine-only `_POF_PROCESSES` tuple, disabled by
  default at the compile seam.
- **Fork 2 (opt-in mechanism) → pure-enable key `pof_positive`** (owner-decided, pros/cons presented;
  advisor flagged it a genuine fork). POF+ is a *binary strain trait*, so the key only enables the
  Process (present/>0 ⇒ on); the rate stays the YAML `k_pof_decarb` — chosen over the D-34 autolysis
  rate-override idiom for fidelity (no physical "half-POF" strain) over pattern-uniformity. The gate is
  **wholly independent of `brett_pitch_gpl`** (a POF+ ferment need not have Brett; a POF-negative
  default wine must make no vinylphenol) — a distinct compile-seam branch, `test_pof_gate_is_\
  independent_of_the_brett_pitch` pins the orthogonality.
- **Fork 3 (carbon routing) → from `hydroxycinnamics`, identical to `BrettDecarboxylation`** (forced:
  it *is* the same chemical reaction, yeast-catalysed). Reuses `M_P_COUMARIC`/`M_VINYLPHENOL`/`M_CO2`,
  `touches=("hydroxycinnamics","vinylphenols","CO2")`, closes on the existing ledger with no new
  conservation code. When POF+ and Brett are both active they draw the *same* `hydroxycinnamics` pool
  (both close 9 = 8 + 1) — verified by `test_pof_carbon_closes` (POF+ alone, and POF+ with Brett).
- **Rate — flux-coupled, a graft of three tested precedents.** `r = k_pof_decarb · X · S/(K_sugar_\
  uptake+S) · [hc]/(K_hydroxycinnamic+[hc])`: rate structure ← `EsterSynthesis`/`AcetolactateExcretion`
  (`fermentative_flux_shape`, catalyst = viable yeast, NOT `X_brett`); carbon routing ← `BrettDecarbox\
  ylation`; gating ← the autolysis opt-in tuple. POF decarboxylation is a *primary-fermentation*
  phenomenon, so the flux term makes it track fermentative activity and **stop at dryness** (S→0 ⇒
  rate 0), leaving the reservoir for a later Brett — which is exactly what pre-fills it during AF. Reuses
  `K_hydroxycinnamic` (same whole-cell precursor affinity as Brett) and `K_sugar_uptake`. **No** Brett
  SO₂/temperature gate (this is yeast metabolism during AF).
- **Fork (temperature) → temperature-flat, no `E_a_pof`** (owner-decided). Cites `AcetolactateExcretion`
  (explicitly T-flat): temperature already enters through the AF-flux trajectory, and no pt4 behaviour
  needs POF's *intrinsic* temperature direction, so an unsourced `E_a` would buy nothing (prime
  directive #2). The ester beat carried its own Arrhenius only because temperature was *that* beat's
  subject.
- **The test-design crux (advisor-caught, load-bearing).** The **stranding** test is the PRIMARY
  headline (`test_pof_strands_vinylphenols_without_brett`, the pt1 control-difference parallel): POF+
  opted in, Brett never pitched ⇒ `vinylphenols` accumulate and strand, `ethylphenols` stays **exactly
  0 and VALIDATED** (no enabled Process touches it — the reductase is Brett's), while `vinylphenols`
  honestly reports speculative — timing-independent and unambiguous. The **head-start** comparison
  (`test_pof_gives_brett_a_head_start`) is the richer SECONDARY test, framed as an **early-time /
  time-to-threshold** claim, NOT an endpoint one: with the same total hydroxycinnamics in both arms,
  conservation forces the *asymptotic* ethylphenols **equal** (all hc → ep eventually), so asserting
  higher *final* ep would be wrong (the arms converge — measured: ~30× ahead at day-12, ~30 days sooner
  to threshold, but only 1.07× at day-120 as POF− catches up). Asserted at day pitch+3 (POF+ > 5× POF−)
  and via time-to-threshold.
- **Empirical tuning.** `k_pof_decarb` = 2.5e-6 mol/(g·h) (speculative) lands ~49 % conversion of a
  100 mg/L must hydroxycinnamic pool during AF (vp ≈ 36 mg/L stranded, hc_resid ≈ 50 mg/L) — a clean
  midpoint leaving a real stranded reservoir *and* residual precursor for the head-start arm. Sourced
  direction (POF+ yeast decarboxylates hydroxycinnamics to vinylphenols during AF — Chatonnet 1992/1997;
  Suárez 2007 review; PAD1/FDC1); magnitude an estimate.
- **Isolability + tests.** A POF-negative default run is byte-for-byte the validated core with all three
  phenol slots VALIDATED (`test_pof_negative_default_is_inert`). **+8 tests** (`test_brett.py`):
  stranding headline, head-start, decarboxylase stoichiometry/`touches`, flux-coupled guards (no
  precursor / dryness / dead yeast all → zero), carbon closure (alone + with Brett), the
  POF-independent-of-Brett gate + default isolability, and the `speculative` tier. `test_media.py`
  wine kinetic-set gains `yeast_pof_decarboxylation` (a new `POF_PROCESSES` set). New param
  `k_pof_decarb` (speculative). **543 green** + 5 benchmark, ruff + mypy clean. **D-40 (and the last M2
  physics beat) complete.** Deferred v2: POF conversion efficiency vs fermentation temperature (would
  add `E_a_pof`); vinylguaiacol/vinylphenol split (currently lumped, as for Brett).

## D-41 — MLF v2: benign senescence (`MalolacticSenescence`) — the slow baseline *O. oeni* decline

**Status: IMPLEMENTED 2026-07-06** (552 green + 5 benchmark, ruff + mypy clean). Lifts the owned
v1 tradeoff of D-39 (*"without SO₂, bacteria never die"*): a new **`MalolacticSenescence`** Process
gives *Oenococcus oeni* a small, always-on-when-pitched **baseline mortality** so a pitched, untreated
dry wine slowly loses its bacteria over **weeks-to-months** (age / ethanol / low-pH / nutrient stress)
even with no SO₂ and no rack — instead of holding a viable culture forever. It moves viable `X_mlf`
into the *same* non-viable `X_mlf_dead` pool the D-39 SO₂ kill uses, so the `X_mlf`-scaled activities
(conversion, citrate → diacetyl, lees-contact diacetyl reduction) fade as the population ages.

**The law — a constant baseline rate, Arrhenius temperature, and nothing else.**
`r_sen = k_senescence_mlf · X_mlf · arrhenius(T, E_a_death_mlf, T_ref)`. Total *O. oeni* mortality is
now `r_sen + r_death` (benign baseline + SO₂-induced), the two built as **separate isolable Processes**
so the D-39 SO₂ lever stays **byte-for-byte** as built and this baseline toggles off independently
(prime directive #3).

- **Environment-free — the load-bearing D-39 crux, reused (advisor-confirmed).** Senescence carries
  **no pH, ethanol, or SO₂ term**. "Benign" *means* environment-independent, and the reason is exactly
  the bug that deferred this to v2: coupling death to ethanol via the Luong wall drives `1 − g_EtOH` to
  ~0.92 at ordinary post-AF ethanol, wiping the culture out in ~1 week instead of the ~2 months reality
  shows. A **constant** baseline dodges it. The ethanol/starvation *modulation* of the baseline stays a
  documented deferral (a further v2 refinement), NOT reintroduced here.
- **Arrhenius, NOT the cardinal γ(T) (the D-39 temperature choice reused).** Warm accelerates the
  decline, cold slows it to dormancy — the physically correct direction. γ(T) peaks at `T_opt_mlf`
  (23 °C) and vanishes past `T_max_mlf`, which would make senescence *maximal at the growth optimum*
  and *switch off* in the warm — backwards for a decline. Reuses `E_a_death_mlf`/`T_ref` (no new
  temperature param); factor 1 at the 20 °C benchmark.
- **Magnitude.** New speculative `k_senescence_mlf` = **5e-4/h** ⇒ half-life ~**58 d** (~8 weeks) at
  `T_ref`, ~100× below the full-SO₂-kill `k_death_mlf` (0.05/h). Negligible over the ~4-day
  co-inoculation MLF window (~5 % `X_mlf` loss over 96 h), so the **§2.2 benchmarks and the D-23
  deacidification control-difference (asserted range [0.1, 0.3], nominal 0.1813) still pass**
  (verified: 5/5 benchmark green — the ~5 % loss leaves it essentially unmoved, well inside the
  range); a stabilizing SO₂ dose still crashes the population in ~1–3 d on top of it.

**Conservation — the carbon/nitrogen-neutral transfer, no new code (the D-13/D-39 pattern).** Both
`X_mlf` and `X_mlf_dead` are weighted at the *same* biomass fractions (since D-38/D-39), so
`d[X_mlf] = −r_sen`, `d[X_mlf_dead] = +r_sen` is C- and N-neutral by construction. **`X_mlf_dead` is a
terminal sink** (advisor blind-spot #3, verified): `YeastAutolysis` reads only the yeast `X_dead`
pool, so senescing bacteria do **not** refuel `amino_acids` — no self-cancelling recycling loop.
Touches only `(X_mlf, X_mlf_dead)`.

**Isolability + performance.** Reads no SO₂ and no pH, so it **never triggers a `brentq`** — strictly
cheaper than the SO₂ kill (`X_mlf ≤ 0` guard only). Pitch-gated: added to the *single source of truth*
`_MLF_GATED_PROCESSES` (compile) / `_MLF_PROCESSES` (media), so it is disabled unpitched and re-enabled
by a `pitch_mlf` intervention exactly like the D-39 death — bacteria age whether or not amino acids
were dosed. On a pitched run `X_mlf`/`X_mlf_dead` stay **speculative** (already so under D-39).

**Tests + the v1 assertions that flip (advisor blind-spot #1).** New `test_malolactic.py` section
(+9): the v2-headline neutral-transfer-without-SO₂, environment-free (SO₂- *and* ethanol-independent
rate), warm-accelerates/cold-preserves-via-Arrhenius, slow-relative-to-the-SO₂-kill (~100×), the
`(X_mlf, X_mlf_dead)` `touches`, no-pH-solve `reads` pin, C/N closure over a senescence-active
(no-SO₂, no-growth) run, and the speculative tier. **Flipped v1 assertions updated, not weakened:**
the ex-`..._no_so2_is_inert` integration test becomes *slow senescence decline + sharp SO₂ crash*; the
two growth-isolation tests (`test_mlf_growth`) disable senescence (or difference it out — the mid-run
`pitch_mlf` re-enables the whole gated set, so the growth signal is isolated as a growth-on−off
control difference in which senescence cancels); `test_media` wine kinetic-set gains
`malolactic_senescence`. **The MLF arc (D-23 → D-31 → D-38 → D-39 → D-41) closes its last deferral.**
Deferred further-v2: ~~ethanol/starvation modulation of the baseline~~ (**IMPLEMENTED in D-52**, see
below); a `BrettSenescence` twin (the same pattern) for the D-40 arc remains open — **deliberately
declined** in D-52's framing, see that entry.

## D-42 — H₂S CO₂-stripping sink (`HydrogenSulfideVolatilization`): residual vs cumulative produced

**Status: IMPLEMENTED 2026-07-06** (561 green + 5 benchmark, ruff + mypy clean). The D-29 forward
note's deferred follow-up, and the last open item on the aroma beat: H₂S production (D-29) was
**produced-only**, so the `h2s` pool was *cumulative produced* (~0.5–1 mg/L) and **overstated
residual** — real fermentation sweeps ~all H₂S out with the CO₂ stream, leaving the µg/L residuals
the sensory threshold (~1–2 µg/L) sits on. This beat adds the CO₂-stripping sink that lifts the
overstatement, the **exact ester D-20/D-21 precedent** (Henry's-law gas stripping) but **carbon-free**,
so *simpler*: neither pool is on any conservation ledger, so the liquid→gas transfer is neutral by
construction (no weighting, unlike `esters`→`esters_gas` in `total_carbon`).

**The mechanism — a flux-linked, first-order Henry's-law sink (the ester mirror).**
`HydrogenSulfideVolatilization` (in `core/kinetics/hydrogen_sulfide.py`) moves dissolved `h2s` into a
new carbon-free `h2s_gas` headspace pool:

    d(h2s)/dt = -k_h2s_volatil · X·S/(K_sugar_uptake+S) · f_gas(T) · f_part(T) · h2s   (into h2s_gas)
      f_gas(T)  = arrhenius(T, E_a_uptake)     — the CO₂ GAS-FLOW factor (stripping rides the CO₂ stream)
      f_part(T) = arrhenius(T, dH_h2s_volatil) — the gas/liquid PARTITION (van't Hoff Henry's-law)

* **First-order in dissolved H₂S, flux-linked, stops at dryness** (`flux → 0`), exactly the ester
  sink: all produced H₂S is co-temporal with a CO₂ stream that can strip it (production is likewise
  flux-linked). The problematic *post-fermentation / autolytic* H₂S that persists **because** no CO₂
  sweeps it is out of scope (the ester sink's omission of slow passive post-cap evaporation).
* **The flux cancels in the residual (the load-bearing structural point, advisor-confirmed).** Because
  production and stripping share the fermentative flux, the residual quasi-steady-state
  `h2s_ss = k_h2s·gate / (k_h2s_volatil·f_gas·f_part)` has the flux **cancel** — residual H₂S tracks
  the inverse-N gate and temperature, **not the ferment speed**. It **rises as `N` depletes** (the gate
  opens) then **freezes at dryness** (both terms gate off with the flux together). Verified empirically:
  residual rises monotonically to a plateau then holds (final == running max).

**Magnitude (prototyped to the physical anchor).** `k_h2s_volatil` = **1.0 L/(g·h)** (speculative)
sizes the stripping so residual sits at the µg/L sensory scale while cumulative produced stays at the
D-29 mg/L magnitude: at `T_ref` with the gate open, `h2s_ss = k_h2s/k_h2s_volatil = 2e-6/1.0 = 2 µg/L`
against ~0.5–1 mg/L produced ⇒ **~99.6–99.7 % stripped** (verified: residual 3.73 / 2.00 / 0.91 µg/L,
produced 0.89 / 0.56 / 0.31 mg/L at 14/20/28 °C). ~100× the ester coefficient (5e-3) — physically right,
H₂S is far more volatile than ethyl acetate.

**`dH_h2s_volatil` sourced, value AND sign (advisor sharpening #2 — the one figure not recalled).**
`dH_h2s_volatil` = **17 500 J/mol** (plausible-in-form/speculative-magnitude), from the Sander Henry's-law
compilation (doi:10.5194/acp-15-4399-2015): −d ln kH/d(1/T) ≈ **2000–2300 K** across sources (Wilhelm
1977, Carroll & Mather 1989, De Bruyn 1995), midpoint ~2100 K ⇒ dissolution enthalpy ~−17.5 kJ/mol
(**exothermic** ⇒ Henry volatility **rises** with T), so a **positive** dH in `arrhenius_factor` (same
sign as the ester's +45 kJ/mol, weaker lever, Q10 ≈ 1.3). **Honesty consequence flagged:** production is
held T-flat (D-29) while stripping rises with T, so the model emits an emergent *"residual H₂S falls with
a warmer ferment"* (3.73 → 0.91 µg/L, 14 → 28 °C). Physically reasonable (warm ferments purge sulfide)
but **unbenchmarked** and reality is mixed (warmth also raises production / N-demand, held flat here) —
tagged directional/speculative and named as an artifact of the T-flat production choice.

**Isolability + conservation (advisor sharpening #1 — the ledger trap avoided).** Both `h2s` and `h2s_gas`
are carbon-free and on **no** ledger (unlike `esters_gas`, which *is* weighted in `total_carbon`), so the
transfer is neutral on every conservation sum **by construction** — **no `conservation.py` change**. The
carbon-closure test is *not* ported; its replacement is the produced-total invariant:
`h2s + h2s_gas` (sink on) equals the sink-off `h2s` trajectory to ~1e-5 (`test_produced_total_is_invariant_
to_stripping`). Isolability holds two ways: dropping the whole `_H2S_PROCESSES` tuple leaves every other
column byte-for-byte (nothing reads `h2s`/`h2s_gas`); dropping **just** the sink recovers the D-29
produced-only `h2s` byte-for-byte (`h2s_gas` stays exactly 0). Both are **always-on in both media** (the
ester/VDK/acetaldehyde intrinsic-metabolism pattern) — the sink Process joins the producer in
`_H2S_PROCESSES`. **Params in the shared `hydrogen_sulfide.yaml`** (medium-agnostic — one physical
mechanism, no per-medium split, unlike the ester `dH` whose *synthesis* direction differs by beverage).

**Tier.** The sink Process is **plausible** in form (CO₂-stripping by the evolving gas is well-understood
Henry's-law physics, the standard explanation for the µg/L residual), with speculative rate params that
cap `h2s`/`h2s_gas` at speculative via parameter-tier propagation (D-1) — no headline change (`h2s` was
already speculative from production; `h2s_gas` is a fresh pool nothing reads, so no other column's tier
drops).

**Schema + the v1 assertions that flip (advisor sharpening #3).** New `h2s_gas` slot in `_common_specs`
(both media: wine schema 32→33, beer 19→20). New `test_hydrogen_sulfide.py` sink section (+8): the neutral
liquid→gas transfer, first-order-in-`h2s`, stop-at-dryness, the ≥0 guard, the physical T-partition lever,
the produced-total invariant, residual-rises-then-freezes-and-produced-plateaus, and residual-falls-with-a-
warmer-ferment. **Flipped run-level assertions updated, not weakened** — every place that read `h2s` meaning
"produced" now reads `h2s + h2s_gas`: the ex-`..._produced_only_and_plateaus` test becomes the residual/
produced split; the low-YAN-early and muted-cross-must levers (`test_hydrogen_sulfide`), the D-30
cap-restores-the-lever (`test_carrying_capacity`), and the two DAP-intervention H₂S tests
(`test_interventions` — production *rate* is now the gradient of the produced sum). **The §3.2 aroma beat
is complete** (esters, fusels, VDK/diacetyl, acetaldehyde, SO₂ free/bound speciation, H₂S production +
stripping); Milestone 2 physics closes. Deferred: the post-fermentation / autolytic H₂S source (persists
un-stripped); a copper-binding / mercaptan model.

## D-43 — Nitrogen redesign: a spike proves default-on residual *assimilable* N is Coleman-incompatible (the redesign as scoped is not worth building)

**Status: DECIDED (not built) 2026-07-06** — a decision-forcing pre-check (the D-26/D-30
"measure before writing" discipline) that resolves, and closes, the recurring "nitrogen-model
redesign" thread carried in the Deferred section since D-23/D-29/D-30. **Outcome: do not build
the large N-model refactor the backlog implied; record why, keep the D-30 opt-in cap, correct the
deferred-note framing.** No source change; the throwaway spike lives outside the repo
(`M:\claud_projects\temp\n_redesign_spike`).

**What was on the table.** The Deferred "residual-nitrogen / satiation floor" note called for "a
nitrogen-model redesign (explicit assimilable vs proline/non-assimilable pools + a satiation
floor)" to make a *default-on* residual-N model possible — motivated by the muted D-29 H₂S
cross-must lever (D-30) and a prospective MLF-with-growth N gate (D-23). The advisor's sharpening:
D-30's *opt-in biomass cap* conflated two separable levers — **total N consumed** (cutting it
reduces biomass → breaks Coleman's sugar curve) versus **N depletion *timing*** (untested). The
untested hypothesis: a preferential/two-pool or cell-quota model that consumes ~all N eventually
(preserving total biomass → preserving Coleman sugar) but drains it *later*, so N persists at
dose-dependent levels through days 1.5–4 and high-YAN musts suppress H₂S longer — a timing lever,
not a total lever.

**The spike.** Standalone scipy, built on the exact Coleman eqs 1–8 RHS the engine matches to
RMSE ~1.3 (`test_coleman_reconstruction`), sweeping two mechanisms: (A) preferential two-pool
(growth reads a fast pool, a slow pool refills it) and (B) a Droop cell-quota (uptake → internal
quota, growth from quota — the textbook decoupling). Measured Coleman sugar RMSE at 80 & 330 mg N/L
vs the H₂S cross-must span at 80/150/300. **Results:** (A) *stalls* — throttling N access throttles
biomass, ferment sticks (~135 g/L residual sugar). (B) finishes and *does* leave dose-dependent
residual N (N@day1.5 ≈ 44 vs 235 mg/L at 80 vs 300), but exposes a clean **anti-correlation** —
lower sugar RMSE ⟺ faster uptake ⟺ less residual ⟺ weaker lever; the Pareto frontier never enters
(RMSE<2, span≫muted). (Caveat, owned: the spike's Droop variant dropped Coleman's active-biomass/
ethanol-death submodel `k_d=k'_d·E`, so its absolute RMSE floor is *contaminated* and NOT citable
as proof — only the anti-correlation direction is clean. The real proof is the argument below.)

**The airtight refutation is mass balance, not the sweep.** Coleman builds biomass fast:
`μ ≈ μ_max ≈ 0.095/h` (`K_N=0.0088` is negligible whenever N>0), so growth is essentially done by
~day 1.3. To match Coleman's *sugar* curve you must match that biomass trajectory; mass balance then
pins `∫uptake = Δbiomass/Y_XN` on that same fast schedule ⟹ **external assimilable N ≈ 0 by ~day
1.3 for every dose**. A quota/luxury buffer can only make external N drain *faster*, never slower,
while biomass stays on Coleman's schedule. So the H₂S-flux window (days ~2–8, biomass high + sugar
present) sees N ≈ 0 regardless of dose — the dose-dependent residual N lives in days 1–3, the flux
weight lives in days 2–8. The only way to widen the lever is a *permanent* residual assimilable N
(the D-30 cap: N never reaches 0 at high dose), which necessarily means less biomass ⇒ breaks
Coleman (D-30 measured sugar RMSE 27.84 at 330). **Conclusion: you cannot hold Coleman AND widen
the H₂S lever via the N model, regardless of mechanism.**

**The reframe that dissolves the apparent D-30-vs-note contradiction.** The deferred note bundled
**two mechanisms with opposite Coleman-compatibility**: (1) an **assimilable-vs-proline split** is
Coleman-*safe* — Coleman's `n0` *is* YAN (assimilable); proline was never in it, so growth-on-YAN
leaves the sugar curve untouched and it is *default-on-able* — but proline does **not** feed the
H₂S gate (correct: not assimilable anaerobically) and nothing reads it today, so it is honest
bookkeeping / inert scaffolding until a consumer exists; (2) a **satiation floor leaving residual
*assimilable* N** breaks Coleman (the 27.84), is what the lever needs, and is inherently *opt-in*.
D-30 ("can't be default-on") is right for the **lever**; the note's "real fix" is right for
**fidelity**. Different axes, both correct.

**Strategic consequence.** The backlog premise "nitrogen redesign unblocks default residual-N /
MLF-with-growth" **largely does not hold on the assimilable axis** — the reason the big build was
declined. The genuine forks, for the record: (a) default-on proline/total-N accounting (honest but
inert until a consumer exists); (b) if the H₂S cross-must lever is the goal, re-point the *gate*
onto a dose-correlated proxy (initial-YAN / intracellular-N-status) — an **H₂S-model** change, not
an N-model one — the clean path to a default-on lever; (c) keep residual-assimilable-N / satiation
**opt-in** (the existing D-30 cap, possibly reformed) — **chosen**; (d) re-anchor away from Coleman
(the only route to default-on residual assimilable N; milestone-scale and data-gated — needs a
dataset with measured residual YAN + sugar). Owner picked (c): keep the D-30 cap as-is, no N-model
build. The negative result **is** the deliverable — it closes a question open since D-23.

## D-44 — Post-fermentation / autolytic H₂S source + copper fining (the two D-42 deferred items)

**Status: IMPLEMENTED 2026-07-06** (578 green + 5 benchmark, ruff + mypy clean). The two follow-ups
D-42 named at close (line "Deferred: the post-fermentation / autolytic H₂S source … a copper-binding /
mercaptan model"): the *reductive-fault* H₂S the flux-linked D-29/D-42 pair could not represent, and
its standard remediation. Ships in two parts; the **mercaptan pool is deliberately deferred** to a
scope decision (below).

**Part 1 — `AutolyticHydrogenSulfide` (`core/kinetics/hydrogen_sulfide.py`): a yield on the autolysis
flux.** As dead yeast self-digest they release intracellular sulfide (cysteine/methionine/GSH). This
is that release, coupled to the **same** first-order autolysis flux `YeastAutolysis` (D-34) runs:

    d(h2s)/dt = y_h2s_autolysis · (k_autolysis · arrhenius(T, E_a_autolysis, T_ref) · X_dead)

* **Yield-on-flux, not an independent rate (advisor steer).** It *recomputes* `YeastAutolysis`' own
  rate internally and scales by a yield `y_h2s_autolysis` [g H₂S / g biomass autolysed] — the D-33
  `FuselAminoAcidReroute` recompute-the-producer idiom. So the `autolysis_rate_per_h` opt-in (which
  *overrides* `k_autolysis` to sweep the sur-lie timescale, D-34) moves **peptide and sulfide release
  on one clock**; an independent constant would desynchronise the two halves of one self-digestion.
* **Why this is the reductive fault — the load-bearing contrast with D-29/D-42.** Unlike the D-29
  producer, this source is **not flux-linked** (first-order in `X_dead`, which persists post-dryness).
  The D-42 CO₂-stripping sink *is* flux-linked, so it **gates off at dryness** and cannot sweep this
  H₂S — post-fermentation autolytic sulfide **accumulates as residual**, the un-stripped "reduction"
  that develops on the lees and calls for racking / aeration / copper. Verified: with autolysis opted
  in the residual keeps rising deep post-dryness (day 15 → 40), > 5× the stripped-to-µg/L default;
  the default run still freezes at dryness (final == running max) — the rise is the new source, not
  the run length.
* **HONESTY CAVEAT #1 — emergent post-dryness, not an AF/post-AF switch.** `X_dead` accumulates
  *during* AF too (D-13 inactivation), so the Process fires whenever `X_dead > 0` and autolysis is on
  — including late AF, where the sink is still active and *does* strip the fresh release. The
  "persists un-stripped" character is **emergent** (flux→0 gates the sink off), not hard-coded.
* **HONESTY CAVEAT #2 — `h2s_gas` semantics broaden.** With autolysis on, `h2s + h2s_gas` ("cumulative
  produced", D-42) now sums sulfate-reduction *and* the stripped fraction of autolytic H₂S; the
  residual `h2s` additionally holds the un-stripped autolytic accumulation. Not a defect — the pools
  stay individually meaningful — but the D-42 "produced total" reading is no longer purely
  sulfate-reduction once the opt-in is set.
* **Magnitude anchored on biomass sulfur (advisor's provenance point).** `y_h2s_autolysis` = **2e-5**
  g H₂S / g biomass (speculative). Ceiling: all biomass S (~0.1–0.4 % dry wt) leaving as H₂S is
  y_max ≈ 1.06 · 2.5e-3 ≈ 2.6e-3; the model's 2e-5 is ~**1 %** of that — the trace released as free
  sulfide, the rest retained in released S-amino acids / peptides / GSH (untracked). Band 2e-6–2.6e-4.
* **Isolability + tier.** Carbon-free, touches **only** `h2s` (nothing reads it back), so the D-34
  isolability holds: opt-in and **wine-only**, disabled *together with* `YeastAutolysis` at the compile
  seam absent `autolysis_rate_per_h`. An undosed wine run is byte-for-byte the validated core — the 5
  §2.2 benchmarks pass unchanged (**run, not inferred**). Both autolysis Processes now sit in
  `_AUTOLYSIS_PROCESSES`. Speculative (already speculative on `h2s` from D-29 ⇒ no tier headline).

**Part 2 — `add_copper` intervention (`scenario/compile.py`): copper-fine H₂S out.** The remediation.
Copper (Cu²⁺, dosed as copper sulfate) precipitates dissolved sulfide as insoluble CuS (Cu²⁺ + H₂S →
CuS↓ + 2 H⁺, **1:1 mol**), settling out with the lees. The verb doses `copper_mgl` (mg/L Cu), converts
to the bindable H₂S mass via the sourced `copper_h2s_binding` = M_H2S/M_Cu = **0.536 g H₂S/g Cu**
(additions.yaml, plausible — stoichiometry exact, complete-binding an idealisation, banded to ~50 %
efficiency), and removes `min(h2s_present, capacity)` — copper in excess simply clears all dissolved
H₂S, the real outcome. **Ledger-neutral by construction** (H₂S carbon/nitrogen-free ⇒ removal books a
zero-weight external flow, the `add_so2` precedent — no `conservation.py` change). Guard: a negative
`h2s` undershoot is left untouched (copper never *adds* H₂S). Verified end-to-end: a day-38 fining of a
reductive (autolysing) wine collapses the accumulated residual to < 25 % of the un-fined run. **SCOPE
(v1):** the removal lever only — residual copper (excess Cu, a haze/toxicity concern) is untracked, and
copper binding of *mercaptans* is deferred with the mercaptan pool.

**Deferred — the mercaptan pool (taken to the owner, D-44).** Copper also binds mercaptans
(methanethiol / ethanethiol), the *other* reductive off-aromas. Two forks make this a scope decision,
not a detail: (a) mercaptans genuinely **carry carbon** (unlike H₂S), so a `mercaptans` pool must be a
real `total_carbon` species (+ copper removal then books a carbon flow like racking debris) — a
carbon-free lump would violate the exact-from-formula discipline and is **not on the table**; and (b)
formation is genuinely murky — methanethiol is mostly methionine-degradation / autolysis, *not* a clean
H₂S→thiol step, so a fabricated conversion rate is exactly what fidelity rejects. Owner to choose:
pool-with-carbon vs. copper-on-h2s-only for v1; and if pooled, autolysis-linked vs. H₂S-linked
formation.

## D-45 — Mercaptan (thiol) pool + copper mercaptide (the D-44 deferred fork; owner chose Option A)

**Status: IMPLEMENTED 2026-07-06** (594 passed incl. the 5 §2.2 benchmarks, ruff + mypy clean). The
D-44 deferred mercaptan-pool fork, resolved with the owner: **build a carbon-bearing pool** (not h2s-only
only), formed **autolysis-linked with carbon drawn from ``amino_acids`` and the nitrogen deaminated**
(Option A), and copper binding it **stoichiometrically** (Cu(SR)₂, 1 Cu : 2 thiol). The H₂S→thiol
formation route was rejected as chemically murky (owner + advisor).

**The pool — lumped ``mercaptans`` booked as methanethiol (CH₃SH, C1, N-free).** Methanethiol is the
dominant reduction thiol (cooked-cabbage, threshold ~2–3 µg/L; ethanethiol its sibling), the honest
single-species stand-in (the arginine-for-``amino_acids`` idiom). New ``M_METHANETHIOL`` in
``chemistry`` registered in all three dicts (1 carbon, 0 nitrogen). New ``mercaptans`` slot in
``wine_schema`` (wine 33→34; beer unchanged — wine-only, with ``amino_acids``), weighted in
``total_carbon`` (as methanethiol) and **absent from ``total_nitrogen``** (N-free).

**Part 1 — `AutolyticMercaptan` (`core/kinetics/mercaptans.py`), Option A.** A yield on the shared
autolysis flux, but — because methanethiol carries carbon (unlike H₂S) — it draws that carbon from
``amino_acids`` and **deaminates** the nitrogen to ``N`` (the exact D-33 ``FuselAminoAcidReroute``
idiom):

    r_merc         = y_mercaptan · autolysis_flux(y) · [aa/(K_amino_acids+aa)]     [g MeSH/L/h]
    d[mercaptans]  = +r_merc
    d[amino_acids] = −(r_merc·c_merc)/c_aa       (arginine mass carrying that carbon)
    d[N]           = +(that mass)·y_N            (DEAMINATION → ammonium)

* **Conservation closes on both ledgers by construction** (advisor-verified): carbon into
  ``mercaptans`` (``r_merc·c_merc``) equals carbon out of ``amino_acids`` (``aa_mass·c_aa``) — the
  draw is sized to match; all the arginine nitrogen leaving ``amino_acids`` lands in ``N`` (MeSH is
  N-free). Both to machine precision; **no new conservation code beyond weighting ``mercaptans`` in
  ``total_carbon``**. Pinned at the derivative level *and* on a full autolysis-on compiled run.
* **PROVENANCE CAVEAT — the arginine lump, not literal methionine (advisor #1).** Real methanethiol
  is from methionine degradation, and Option A's rationale was "methionine is a released amino acid"
  — but ``amino_acids`` is booked as *arginine*, so the carbon/N drawn are arginine's. The model
  deaminates ~0.66 mol N per mol MeSH vs. methionine's ~1 (same order, no gross artifact). Documented
  as the arginine-for-``amino_acids`` stand-in: **exact on the ledger, approximate on provenance** —
  *not* faithful methionine chemistry.
* **New TIER consequence — a structural drop on ``N`` (advisor #2, the D-27 ``E`` parallel).**
  ``AutolyticMercaptan`` is the **first autolysis-gated Process to write ``N``** (via deamination), so
  an autolysis-on run drops the *structural* ``tier_of("N")`` PLAUSIBLE→SPECULATIVE — **even on an
  autolysis-on / amino-dose-off run**, where the other N-writer (``FuselAminoAcidReroute``) stays
  disabled. Verified by run (default N tier = plausible; autolysis-on = speculative) and pinned
  (``test_is_the_first_autolysis_gated_n_writer``). The param-aware tier was typically already
  speculative, so no headline change.
* **Availability gate + not-flux-linked.** ``aa/(K_amino_acids+aa)`` ramps production to 0 as the
  pool empties (solver-safe, can't drive ``amino_acids`` negative); the D-34 refill keeps it
  non-empty. First-order in ``X_dead`` (not fermentation flux), so — like the D-44 H₂S source — it
  accumulates un-stripped post-dryness. New speculative ``y_mercaptan`` = **1e-5 g MeSH/g biomass**
  (~0.6 % of the biomass-methionine ceiling; set *below* the H₂S yield since reduction skews to H₂S).

**Part 2 — `add_copper` extended to bind mercaptans (Cu(SR)₂, 1 Cu : 2 thiol).** Copper binds **H₂S
first** (CuS Ksp ~10⁻³⁶ ≫ mercaptide, so sulfide is preferential), then binds mercaptans with the
**leftover** copper. New ``copper_mercaptan_binding`` = **1.514 g MeSH/g Cu** (= 2·M_MeSH/M_Cu;
plausible, **banded down hard** to ~20 % — copper is notably incomplete on mercaptans and useless on
the disulfides they oxidise to, advisor #3). **Ledger:** removing carbon-free H₂S is neutral (D-44),
but **removing mercaptans removes carbon** from the wine as the precipitated mercaptide — a *negative
external flow* the driver books (the racking-debris precedent), so ``final == initial + Σ flows``
still holds (verified; conservation-across-jump uses the flow identity, **not** ``assert_conserved``,
advisor #4). The pre-existing "copper is ledger-neutral" test still passes *because default wine has
``mercaptans ≡ 0``* — its comment now says so.

**Shared `autolysis_flux` helper (advisor, non-blocking).** Three Processes now recompute
``k_autolysis·f_T·X_dead``; extracted into one ``autolysis_flux(y, schema, params)`` in
``autolysis.py`` (the ``fusel_production_rate`` single-source idiom), so ``YeastAutolysis`` (D-34),
``AutolyticHydrogenSulfide`` (D-44) and ``AutolyticMercaptan`` (D-45) share one clock and one
``autolysis_rate_per_h`` override.

**Isolability + verification.** Opt-in and wine-only: disabled **together with** the other two
autolysis Processes at the compile seam absent ``autolysis_rate_per_h`` — an undosed wine run is
byte-for-byte the validated core (the **5 §2.2 benchmarks pass unchanged, run not inferred**; the
``mercaptans`` slot is a permanently-zero column there). +27 tests (new ``test_mercaptans.py`` +
copper-binds-both in ``test_interventions.py``); **594 passed (incl. the 5 §2.2 benchmarks)**.

## D-46 — Harden `solve_ph` to be total over ℝ: clamp when the electroneutral pH lies outside [0, 14]

**What broke.** After D-45 shipped, three `test_brett.py` integration tests
(`test_growth_accelerates_phenols`, `test_so2_crashes_growing_brett_population`,
`test_death_run_conserves_carbon_and_nitrogen`) went **red on `main`** — all three raise
`ValueError: f(a) and f(b) must have different signs` from `brentq` inside `acidbase.solve_ph`,
reached via `BrettDecarboxylation.derivatives → ph_of_state`. They were **green at the D-45
parent** (1241ba1). So D-45 regressed them.

**Root cause — a latent fragility D-45 exposed, not a new bug.** The failure is *not* in the D-44/D-45
derivatives: no-op'ing both `AutolyticHydrogenSulfide` and `AutolyticMercaptan` still reproduces it.
The only remaining change is that D-45 **appended the `mercaptans` state slot** (wine 33→34), and that
extra dimension shifts BDF's adaptive step sequence. The mechanism, isolated empirically:

- **RK45 and LSODA succeed**; the `cation_charge` slot is **constant at 0.0254 mol/L** (nothing writes
  it — it is a compile-time back-solve, D-18). BDF alone fails.
- The failing state has `cation_charge = 3.81 mol/L` — two orders of magnitude above physical. This is
  **BDF's `num_jac` Jacobian probe** perturbing the `cation_charge` slot far outside its physical range.
- At that unphysical cation, `charge_residual` is **positive across the whole [0, 14] bracket** (the cation
  swamps all acid buffering), so `brentq` finds no sign change and throws.

`solve_ph`'s fixed `[0, 14]` bracket implicitly assumed a physiological cation. `num_jac` probes *every*
state variable outside its physical range; `cation_charge` is the first that feeds a bracketed root-find,
so `solve_ph` was the first core helper to be **partial** (throwing on a valid-for-num_jac input) rather
than total. Any slot addition that reshuffles BDF stepping could trigger it — D-45 happened to be the one.

**Fix — make `solve_ph` total.** `charge_residual` is strictly monotone *decreasing* in pH, so
residual(0) is its max and residual(14) its min. Evaluate both ends first: both-positive ⇒ the
electroneutral pH is **above 14** (return 14.0); both-negative ⇒ **below 0** (return 0.0); otherwise the
single interior root exists and `brentq` finds it exactly as before. This is **exact, not a band-aid** —
returning the boundary *is* the correct "root lies outside the physical window," and a physiological cation
falls straight through to the identical `brentq` call ⇒ **bit-for-bit pH, byte-for-byte trajectories**
(RK45/LSODA prove the real trajectory never leaves the bracket). The clamp activates only on the num_jac
probe, which affects only the Jacobian (Newton convergence), never the solution — the RHS stays exact.

**Verification.** +3 direct unit tests in `test_acidbase.py` pin the totality at the *function* level
(huge cation → 14.0; strongly-negative cation → 0.0; physiological → unclamped interior root == anchored
target) — the Brett tests only catch it incidentally through a 120-day integration, so a future refactor
that stops triggering the probe must not silently un-total the solver. Full suite **600 passed (incl. the
5 §2.2 benchmarks, run not inferred)**; ruff + mypy clean. Validated core byte-for-byte preserved.

**Noted, not acted on.** (1) `_needs_ph_solve` fires on `so2_total > 0`, so num_jac probing `so2_total`
off exact zero triggers a pH solve in an *un*-sulfited run — a probe-only perf smell (the real trajectory
holds it at 0), not a correctness issue, and it would not fix the two genuinely SO₂-dosed tests anyway.
(2) num_jac probes every state var outside its physical range; a future Process with a `log`/`sqrt`/bracket
that assumes a physical domain could be the next `solve_ph` — harden reactively when exposed, not speculatively.

## D-47 — SO₂-bound acetaldehyde is protected from ADH: the D-28 free/bound split feeds back into the RHS

**The deferred coupling.** D-28 built the free/bound SO₂ split (acetaldehyde-bisulfite binding
equilibrium) but left it **readout-only**: the split did *not* feed back into the acetaldehyde
reduction. The deferred note said "bound acetaldehyde is notionally protected from ADH; that RHS
coupling is deferred." This decision lands it (owner-authorised beat; owner chose the wiring, below).

**Physics — reduce only the free (unbound) acetaldehyde.** Alcohol dehydrogenase reduces
acetaldehyde → ethanol, but the acetaldehyde-bisulfite adduct (1-hydroxyethanesulphonate) is not a
substrate. The literature is explicit: *"acetaldehyde bound to SO₂ could not be metabolized by yeast
during fermentation; only free acetaldehyde could impact metabolism"* — a stable 1:1 complex (Han et
al. 2020, *Food Chemistry*; S. Afr. J. Enol. Vitic. 2018). So `AcetaldehydeReduction` now reads the
**free** share, `free = total_acetaldehyde − bound`, with `bound` from the *same* `bound_so2_molar`
equilibrium the SO₂ readout uses (1:1 ⇒ bound SO₂ mol/L = bound acetaldehyde mol/L). Binding is fast
(98 % in ~90 min) relative to the enzymatic reduction, so the instantaneous-equilibrium (QSS) split is
justified. New pure helper `acidbase.free_acetaldehyde(y, schema, params, ph)`.

**The emergent consequence — SO₂ locks in acetaldehyde.** Because bound acetaldehyde is protected,
a sulfited wine no longer clears its acetaldehyde: it strands a residual. Measured (50 mg/L dose at
pitch): acetaldehyde peaks *higher* than the unsulfited run (72 vs 37 mg/L — reduction is throttled
from early on) and ends at **~27 mg/L stranded** (0.78 mol per mol SO₂), with free SO₂ pinned at ~22 %
of the dose. This is **near-stoichiometric at the stoichiometric edge** (50 mg/L SO₂ ≈ acetaldehyde
molar); at *sub*-stoichiometric field doses the model reduces to the observed ~0.76× degradation-rate
slowdown (Han 2020) and the ~366 µg-acetaldehyde-per-mg-SO₂ (~0.5:1 molar) field figure — so the
mechanism is grounded across the regime, not just qualitatively. The binding constant itself is unchanged
(`K_acetaldehyde_so2 = 1.5e-6`, the D-28 value, which is the literature K at pH 3.3 exactly).

**The retired invariant (the load-bearing change).** D-22/D-28 advertised "**SO₂ is readout-only** —
dosing it perturbs nothing else." **D-47 intentionally retires that for sulfited runs.** SO₂ now couples
into the acetaldehyde trajectory (and, through the SO₂ readout the MLF/Brett gates consume, into those
too). It is preserved *exactly* where it still holds:
- **Undosed runs are byte-for-byte the D-27 core.** The `so2_total > 0` guard is exact — no dose ⇒ no
  per-RHS pH `brentq`, no protection (the MLF/Brett SO₂-gate isolability idiom, D-39). **No §2.2
  benchmark doses SO₂**, so the acceptance suite is untouched.
- **Carbon still closes to machine precision.** The reduction only *throttles* the acetaldehyde→E
  transfer; it neither creates nor routes carbon.
- **pH is still not a charge actor** (~2e-6 drift): SO₂ couples only via acetaldehyde, which carries no
  charge — the D-22/D-28 "SO₂ not in the charge balance" claim is intact.
- **Footprint on the core ferment is second-order** (≤1e-3 of each column's scale): the only ripple is
  the borrowed-ethanol-carbon dip feeding the E→viability brake (the D-27 note). Acetaldehyde itself
  diverges order-unity (the intended stranding); everything else moves only at the E→viability level.

**Owner fork (surfaced before building): bake-in default-on vs opt-in toggle.** Owner chose **bake-in,
default-on** — protection lives in `AcetaldehydeReduction`, active whenever SO₂ is dosed, matching the
MLF SO₂-gate precedent and acetaldehyde's "intrinsic, always-on" framing (D-27). `touches` unchanged
(still `acetaldehyde`/`E`); `reads` unchanged too — the SO₂/pH params are read *inside*
`free_acetaldehyde`/`ph_of_state` and the acetaldehyde/E output is already speculative, so declaring them
would move no tier (the MLF-gate precedent). No new parameters.

**CAVEAT (speculative).** Bound acetaldehyde is treated **inert-to-ADH** — real adduct slowly
dissociates and degrades over months, so the stranding is an **upper bound on persistence** (the
literature's own "not metabolized *during fermentation*"). Dosing SO₂ at pitch is also the *maximal*-
stranding scenario; the common cellar case (SO₂ post-AF, into a wine the yeast already cleared) strands
almost nothing — pinned by `test_post_af_so2_dose_strands_far_less_than_a_pitch_dose`. **Second mechanism
(the induced-over-production half): now modelled in D-48.** D-47 elevates acetaldehyde by *protection
only* (throttled reduction); reality *superimposes* an **SO₂-induced over-production** — a redox pull
where the yeast excretes *more* acetaldehyde because SO₂ traps it (Han 2020). This note originally
scoped that out; **D-48 adds it as a total-SO₂-gated bump to `AcetaldehydeProduction`, scoped to the
transient PEAK** (the end state is already at/above the field slope from protection alone — see D-48 for
why an additive *end-state* term would overshoot).

**Downstream test consequences (all faithful, re-pinned to measured output).**
- `test_so2_dose_suppresses_mlf_in_a_run`: SO₂ dosed *during* AF is now only a **partial** MLF brake
  (retains malic 2.3/4.0, ΔpH +0.07) not a near-total one — stranded acetaldehyde sequesters most of the
  antimicrobial pool ("bound SO₂ is not antimicrobial", emergent and dynamic). The counterintuitive
  direction the change predicts.
- `test_molecular_so2_series_…`: molecular SO₂ now nets **down** over the run (free depressed by stranding
  dominates the pH-fraction rise) — the flip of the readout-only-era direction.
- The former byte-identical isolation test is refocused as
  `test_so2_coupling_strands_acetaldehyde_but_spares_the_core_ferment`.

**Verification.** New `test_acetaldehyde.py` D-47 section: unsulfited byte-for-byte closed form; SO₂
throttles the rate to the free share (comparable-molar and excess); post-AF strands ≪ pitch; carbon
closes on a stranding run; BDF vs RK45/LSODA agreement (the rate is now nonlinear in acetaldehyde/SO₂ via
the `bound_so2_molar` quadratic root on an always-on RHS). Full suite **606 passed (incl. the 5 §2.2
benchmarks)**; ruff + mypy clean. Validated core byte-for-byte preserved.

## D-48 — SO₂-induced over-production: a total-SO₂-gated bump on acetaldehyde production, scoped to the transient peak

**The task and the reversal it forced.** D-47's caveat scoped out the "induced over-production" half of
the SO₂ acetaldehyde elevation ("today it is out of scope"); this beat was authorised to add it — an
SO₂-gated bump to `AcetaldehydeProduction` capturing the redox pull where trapping the terminal electron
acceptor (acetaldehyde) makes the yeast intensify the glyceropyruvic pathway and excrete *more* of it
(Han 2020). Building it surfaced a finding that **reshaped the beat before any code shipped**, and it is
the crux of this decision.

**The empirical finding (why the naïve framing was wrong).** The premise behind the task — that D-47
captured only "protection" and a "production half" was missing from the finished-wine level — is
**contradicted by the model's own numbers**. With the induced bump OFF (D-47 protection only), end-state
acetaldehyde increments per SO₂ dose are **25.7 / 56.1 / 119.0 mg/L at 50 / 100 / 200 mg/L SO₂**, versus
the field correlation `W_acet = −4.4 + 0.39·W_tSO₂` (0.39 mg/mg ⇒ 19.5 / 39 / 78). D-47 protection
**alone already delivers 1.3–1.5× the full field slope** — there is no under-shoot for an additive
end-state term to fill; a bump there would only overshoot. The structural reason: **the finished-wine
level is capped by the SO₂-binding equilibrium (D-28), not by production.** An over-produced slice of
acetaldehyde is reduced back once flux → 0 (D-27 borrow-from-E); only the *bound* fraction survives, and
that pool is saturated — so *any* production driver (free or total SO₂) leaves the end state ~unchanged
(25.7 → 25.8). (The 1.3–1.5× overshoot is itself defensible: thermodynamic 1:1 adduct binding *should*
exceed a net-field regression whose wines bleed acetaldehyde and SO₂ to sinks the model omits; R = 0.837,
−4.4 intercept — a loose anchor. Whether to trim D-47's binding calibration toward the field slope is
left as a possible separate beat, not folded in here.)

**Owner fork (three options presented with pros/cons).** Given the finding, the owner chose **Option 3 —
scope D-48 to the transient PEAK, not the end state.** The mid-ferment peak *is* a real, distinct,
measurable phenomenon (active-ferment over-excretion, later cleared by ADH) with **no end-state literature
anchor**, so D-48 models it without double-counting the end state D-47 already delivers. (Option 1, don't
build, was the runner-up; Option 2, re-split D-47+D-48 to sum to 0.39, was rejected as *illusory* — you
cannot redistribute an end state that is set by the binding equilibrium, not production.)

**Driver = TOTAL SO₂ (reverses the owner's earlier free-SO₂ choice).** The owner had initially picked
**free** SO₂ (a stability premise: negative feedback, self-limiting). The data refuted that premise *for
this observable*: free SO₂ **collapses to ~0 at the peak** (nearly all sulfite bound to the rising
acetaldehyde), so a free-SO₂ driver is **empirically inert on the very peak it targets** (+0.1..+1.4 mg/L
across all doses/k — it self-quenches exactly when needed). Surfaced back to the owner, who switched to
**total** SO₂. Total SO₂ is open-loop but **stable in practice**: the term is flux-gated, so the
fermentative flux → 0 at dryness caps its time-integral (no runaway — verified empirically to 200 mg/L,
end-stranding flat). Feedback topology, for the record: free = negative feedback (stable *because* inert);
total = open-loop (stable via the flux cap); bound-acetaldehyde = positive feedback (most faithful but
reintroduces runaway — rejected).

**The term.** `d[acetaldehyde] += k_acet_so2_induced · flux · so2_total`, guarded on `so2_total > 0`, a
**carbon-exact borrow from E** exactly like the base production (the diverted glycerol carbon reality
shows parks here as acetaldehyde — no glycerol pool, a v1 simplification leaving the ethanol-yield
reduction slightly understated). It reads the **total SO₂ state slot directly** — no param, no per-RHS pH
`brentq` (cleaner than the free-SO₂ variant would have been). Net `dE` stays positive (base + induced
borrow is ≳100× below the uptake ethanol deposit). New parameter `k_acet_so2_induced` (shared
`acetaldehyde.yaml`), **value 4.0e-3 L/(g·h), tier speculative**. At 50 mg/L it lifts the peak **~+3.8
mg/L (71.6 → 75.4)**, dose-scaling to ~+15 at 200; end state and stranded residual **unchanged**. The
`so2_total > 0` guard is **exact** — an unsulfited run is byte-for-byte the D-27/D-47 core (and no §2.2
benchmark doses SO₂).

**Sizing the unanchored knob — a cross-process reality constraint (owner asked "what is closer to
reality").** There is no direct field anchor on the *peak* elevation (the 0.39 slope is end-state and
already met by D-47), so the initial value (1e-2) was unanchored — and building it exposed that D-48 is
**not contained to the peak**: the raised acetaldehyde sequesters SO₂ *during* the ferment, weakening the
molecular-SO₂ MLF brake (a real, textbook effect — "bound SO₂ is not antimicrobial"). At 1e-2 that brake
retained only 44 % of the malic (`test_so2_dose_suppresses_mlf_in_a_run` 2.3 → 1.76 g/L), crossing below
the literature "SO₂ remains a **partial** brake, *more than half* the malic retained" regime the D-47 work
established. This was surfaced to the owner (it is a bigger scope leak than the driver fork), who directed
sizing by reality. **Resolution: k is set to the largest value keeping the MLF brake in that
>half-retained regime.** Measured ceiling ≈ 5e-3 (malic 2.03, at the 2.0 floor); **nominal 4e-3** keeps a
safe margin (malic **2.09, ~48 % converted**) — the MLF test lands back **inside its original
`2.0 < malic < 3.0` band, so the D-48 change no longer perturbs it** (the earlier `1.3 < malic < 2.3`
re-pin is reverted). This represents *both* real phenomena — the induced peak over-production **and** the
partial MLF-weakening — at the largest self-consistent magnitude, rather than letting a free knob push a
validated observable out of its literature regime.

**Verification.** New `test_acetaldehyde.py` D-48 section: exact-when-undosed (base only, `==`); dosed
closed form (base + induced, carbon-exact, SO₂ read-only); peak lifts while end state is unchanged;
peak lift scales with dose (the total-SO₂ signature). `test_production_metadata` now pins
`k_acet_so2_induced` in `reads`; the MLF band is *unchanged* (D-48 sized to keep it). Full suite **610
passed (incl. the 5 §2.2 benchmarks)**; ruff + mypy clean. Validated core byte-for-byte preserved.

## D-49 — Excreted keto-acid overflow pool (pyruvate): the second SO₂-binding carbonyl, as a side pool not an on-pathway precursor

**The task and where it came from.** D-48 flagged that the model's finished-wine acetaldehyde stranding
overshoots the field regression `W_acet = −4.4 + 0.39·W_tSO₂` by **1.3–1.5×** and named it "a D-47/D-28
binding-calibration question." Investigating it established the overshoot is a **real missing mechanism,
not a mis-calibration**: the model routes *100 % of bound SO₂ onto acetaldehyde*, but real wine shares
dosed SO₂ with **competing carbonyls** — chiefly the excreted keto-acids **pyruvate** and
**α-ketoglutarate** (Jackowetz & Mira de Orduña 2013). A pre-check (multi-carbonyl partition at each
end-state free-SO₂ level, sourced competitor pools) confirmed that sharing SO₂ with persistent
finished-wine keto-acids pulls the slope from ~0.56 down toward the field ~0.39 at a typical 50 mg/L dose.
The owner chose to **build the competitor pools** (per-species, dynamic, D-19 side-pool idiom). D-49 is
the first: pyruvate. (α-KG = D-50; the coupled multi-carbonyl SO₂ equilibrium that reads them = D-51.)

**The load-bearing modelling choice: an EXCRETED SIDE POOL, not acetaldehyde's on-pathway precursor.**
The owner's initial instinct (and a first advisor recommendation) was the "maximum-fidelity" rework:
route acetaldehyde's carbon *through* pyruvate, its real metabolic precursor (pyruvate → acetaldehyde +
CO₂ → ethanol). This was **designed and then rejected** (the advisor retracted its own recommendation on
review). It conflates two physically distinct pools: acetaldehyde's precursor is the **intracellular flux
intermediate** (enormous flux, vanishing pool, never persists, never measured), whereas the SO₂-binding
pyruvate is the **extracellular excreted overflow residual** (small flux, persistent, measured). One pool
cannot be both, and the persistence mechanism the rework needed — "SO₂ shields pyruvate from pyruvate
decarboxylase" — is not real (PDC is intracellular; the excreted residual never meets it). Worse, that
shielding would make dosed SO₂ *sequester acetaldehyde's precursor* and **suppress** acetaldehyde — the
exact opposite of the SO₂-induced over-production D-48 just shipped. So the excreted-overflow side pool is
the **more** faithful structure for the quantity that matters here, and **acetaldehyde / D-27 / D-47 / D-48
stay entirely untouched.**

**The model.** `S --excretion--> pyruvate --reassimilation--> ethanol + CO₂`, the D-19/D-26 byproduct idiom.
`PyruvateExcretion` (flux-linked, temperature-flat) draws pyruvate's carbon *out of `S`* at the C3 fraction
while the yeast ferments — stops at dryness. `PyruvateReassimilation` returns it to `E`+`CO2`
(`C3 → C2 + C1`, one mole each — carbon-closing like malic → lactic + CO₂, D-23). Carbon returns to `E`/`CO2`
**not `S`** deliberately: post-dryness `S = 0`, so a refund-to-sugar would be a no-op that *destroys* carbon.
Both processes are **wine-only** (v1) — the SO₂ competition is a wine readout; no §2.2 beer benchmark asserts
a keto-acid level. Both **speculative** (rate magnitudes are estimates); the excreted-overflow *mechanism* is
textbook.

**The mid-build mechanism correction (the crux — flux-link the reassimilation, not viable-X-gate it).** The
first build borrowed the acetaldehyde-reduction template: reassimilation gated on **viable X with no flux
term**. It failed empirically — the finished-wine residual came out **0.0 mg/L**. Cause: a clean ferment
finishes with the yeast **still viable** (~0.4 g/L here), so a no-flux viable-X gate keeps clearing pyruvate
over the long post-dryness tail and drains the pool to ~0. The residual had been (wrongly) pegged to *yeast
death*; a normal ferment doesn't crash. The advisor confirmed this is a **mechanism bug, not tuning**, and
the diagnosis: overflow-pyruvate re-assimilation is **co-metabolic** (tracks active fermentation), the
*opposite* of ADH (which genuinely keeps reducing acetaldehyde through the post-ferment rest). Fix:
**flux-link the reassimilation** (share excretion's `X·S/(K+S)` shape). At dryness both terms die and the
pool **freezes** at its dryness value — a residual pegged to *end-of-fermentation*, hence **crash- and
duration-independent** (verified: 30.0 mg/L at both 21 and 40 days). Consequence, documented as a v1
simplification: because both terms ride the same flux shape the pool rises *monotonically* to the plateau
`k_pyruvate_excretion / k_pyruvate_reassimilation` rather than showing the real mid-ferment peak-then-decline
— but **nothing reads the peak** (D-51 reads only the residual), so the transient is dropped and the growth-
coupled excretion that would restore it (option B) is deferred. Sizing is by the **ratio** (3e-3 / 1e-1 =
0.03 g/L = 30 mg/L), in the real finished-wine range.

**Isolability (prime directive #3).** Own `_KETO_ACID_PROCESSES` tuple; a ProcessSet without it is the prior
core. Unlike the byte-for-byte acetaldehyde buffer, the pool routes a *trace* of sugar carbon on a detour to
ethanol (parking only the ~30 mg/L residual), so the ABV/CO₂ endpoints are **not** bit-identical to the
pool-off core — but the delta is **rel ~4.4e-5** (≪ 0.1 %), so the §2.2 benchmarks are preserved far below
tolerance. Carbon closes to **machine precision** (pool weighted at its C3 fraction in `total_carbon`).

**Verification.** New `test_keto_acids.py` (19 tests): closed forms + carbon-exact draw/release, the dryness-
freeze (reassimilation stops at `S=0` — the load-bearing difference from ADH), the persistent residual in the
finished-wine range *with the yeast still viable*, duration-independence (21 vs 40 days), machine-precision
carbon closure, the ≪0.1 % ABV/CO₂ isolability delta, wine-only wiring, and speculative tier propagation. Six
existing `full_params`/schema tests updated for the new shared YAML + state slot (mechanical). Full suite
**629 passed (incl. the 5 §2.2 benchmarks)**; ruff + mypy clean. **Next:** D-50 (α-KG, same structure), then
D-51 (the coupled multi-carbonyl SO₂ equilibrium that reads both pools — where the overshoot actually drops).

## D-50 — Excreted keto-acid overflow pool (alpha-ketoglutarate): the third SO₂-binding carbonyl, same structure as D-49 with one fix

**The task.** D-49 built pyruvate, the first excreted-overflow SO₂-binding keto-acid, and named α-KG
as "same structure" next. This beat builds it: `AlphaKetoglutarateExcretion` /
`AlphaKetoglutarateReassimilation` in `keto_acids.py`, a new wine-only `alpha_ketoglutarate` state slot
(schema 35→36), wired into `_KETO_ACID_PROCESSES` alongside pyruvate.

**The one design fork: the reassimilation carbon destination.** Two options were considered before
writing code (advisor-consulted): (a) mirror pyruvate exactly — return carbon to `E`/`CO2`; (b) the
"more faithful" route via the *real* α-ketoglutarate-dehydrogenase reaction, α-KG (C5) → succinate (C4)
+ CO2 (C1), landing in the existing `Byp` (succinic-acid-booked) pool. **Rejected (b).** The advisor's
key diagnosis: pyruvate's `C3 → C2(ethanol) + C1(CO2)` reassimilation is *nearly isolable* not because
"return to E" is inherently safe, but because that mole-for-mole split **happens to be exactly the
Gay-Lussac fermentation carbon ratio** (2 carbon to ethanol : 1 carbon to CO2) — so the detour is
stoichiometrically indistinguishable from the main pathway, and pool-on/off differs only by the frozen
residual (rel ~4e-5 endpoint delta, D-49). Routing to succinate/`Byp` instead would divert
reassimilation **throughput** — not just the residual, but ~10–20× more (the pool cycles many times
per ferment) — permanently away from ethanol, large enough to threaten both the §2.2 ABV/CO₂
benchmarks and any `Byp` assertion. Worse, the "fidelity" justification doesn't hold either way: α-KG
dehydrogenase is largely *repressed* under the anaerobic conditions that make α-KG overflow in the
first place, and the real dominant reassimilation fate is glutamate synthesis (α-KG + NH4+ →
glutamate, N-coupled, not modelled in v1) — so neither the ethanol/CO2 route nor the succinate route is
"more biochemically true"; both are lumped carbon-closing stand-ins (the fusel/ester idiom, D-19), and
fidelity is not the tiebreaker. Decision: **route to E+CO2, mirroring pyruvate**, but fix the ratio —
C5 does not divide 1:1 like pyruvate's C3, so mole-for-mole would give a CO2-heavy 1+4 split. Instead
the Process returns carbon at the **same 2:1 Gay-Lussac ratio**: `5/3` mol ethanol + `5/3` mol CO2 per
mole of α-KG (`C5 → C(10/3)` ethanol-carbon + `C(5/3)` CO2-carbon) — carbon-exact, and the general form
(`carbon_atoms/3` mol each) reduces to pyruvate's mole-for-mole case when `carbon_atoms == 3`.

**Everything else mirrors D-49 exactly.** Flux-linked excretion (temperature-flat, draws C5 from `S`)
+ flux-linked co-metabolic reassimilation (stops at dryness, freezing the residual — crash- and
duration-independent, verified at 21 vs 40 days). Both speculative; wine-only (v1, no §2.2 beer
benchmark asserts a keto-acid level). **Residual sized lower than pyruvate's ~30 mg/L**: nominal ratio
`k_alpha_kg_excretion / k_alpha_kg_reassimilation` = 2.0e-3 / 1.0e-1 = 0.02 g/L = 20 mg/L (α-KG is
typically somewhat less abundant in finished wine than pyruvate per the same D-49 sources, Jackowetz &
Mira de Orduña 2013). New `total_carbon` weighting term for `alpha_ketoglutarate` (own C5 fraction,
mirroring the pyruvate term) — caught before it could silently fail the carbon-conservation test.

**Isolability.** Same own `_KETO_ACID_PROCESSES` tuple as pyruvate; a ProcessSet without it is the
prior core. Not byte-for-byte (routes a trace of sugar carbon on a detour), but the ABV/CO2 endpoint
delta with both keto-acid pools on is **measured** (not just threshold-checked) at rel **~7.3e-5** —
roughly double pyruvate-alone's ~4e-5 (D-49), as expected from two detours, still ≪0.1 %. Carbon
closes to machine precision. The residual also lands exactly on the ratio's design target: α-KG
freezes at **20.0 mg/L** (pyruvate unchanged at 30.0 mg/L) on the standard 21-day acceptance run.

**CALIBRATION-PENDING flag for D-51 (advisor-raised).** Both keto-acid residuals (pyruvate 30 mg/L,
α-KG 20 mg/L) are honest order-of-magnitude author estimates, not fits — D-51 must not inherit them
as settled. Two things D-51 needs to re-derive, not assume: (1) the residual *ratio* between the two
pools may need to shift once the multi-carbonyl SO₂ equilibrium is actually fit to the field 0.39
mg/mg slope (D-48); (2) SO₂ binds molar concentration, not mass — α-KG's higher molar mass
(146.1 vs pyruvate's 88.06 g/mol) means 20 mg/L α-KG is only ≈0.137 mmol/L vs pyruvate's ≈0.341
mmol/L, i.e. α-KG's molar contribution to the binding competition is ~40% of pyruvate's despite
being ~67% of it by mass — D-51's equilibrium must work in moles, not the mg/L this beat reports.

**Verification.** Extended `test_keto_acids.py` with 17 new tests mirroring D-49's suite for α-KG
(metadata, closed forms incl. an explicit non-mole-for-mole regression guard, dryness freeze,
carbon-neutral draw/release, wine-only wiring, tier propagation), now 36 tests total; plus updated
acceptance tests to cover both pools together (persistent residual in range, duration-independence,
carbon closure, ABV/CO2 isolability). `test_media.py` schema-size/slot-tuple/`EXPECTED_PROCESSES`
updated (wine schema 35→36). A new `total_carbon` weighting term for `alpha_ketoglutarate` was needed
in `validation/conservation.py` (caught immediately by the carbon-conservation test, not silently).
**646 passed (incl. the 5 §2.2 benchmarks)**, ruff + mypy clean. **Next: D-51**, the coupled
multi-carbonyl SO₂ equilibrium that reads acetaldehyde + pyruvate + α-KG together — the beat both
keto-acid pools exist to feed.

## D-51 — Coupled multi-carbonyl SO₂ equilibrium, worked in moles: the actual D-48 overshoot fix — and an honest partial one

**The task and where it came from.** D-50 flagged both keto-acid residuals (pyruvate 30 mg/L,
α-KG 20 mg/L) as order-of-magnitude author estimates the multi-carbonyl equilibrium must
**re-derive against the field slope**, not inherit as settled — and flagged that the equilibrium
must work in **moles**, since SO₂ binds molar concentration and α-KG's higher molar mass
(146.1 vs pyruvate's 88.06 g/mol) means its 20 mg/L residual is only ~40% of pyruvate's molar
contribution despite being ~67% of it by mass. This beat does both: generalises D-28's
single-carbonyl equilibrium to the coupled three-carbonyl case, and empirically re-derives (rather
than assumes) whether the D-49/D-50 residual sizing actually closes the D-48 overshoot.

**The model: one shared root-find, not a multi-dimensional solve.** D-28's `bound_so2_molar` was a
closed-form quadratic for a single 1:1 carbonyl-bisulfite adduct. Generalising naively to N
competing carbonyls sharing one bisulfite pool would need an N-dimensional simultaneous solve —
but every carbonyl's bound fraction can be written as a Langmuir partition of one shared "reactive
bisulfite" variable `h`: `bound_i = A_i·h/(K_i+h)`, with `h` the unique root of the strictly
monotone-decreasing residual `β·total − β·Σᵢ(A_i·h/(K_i+h)) − h = 0` over `[0, β·total]` (guaranteed
sign change ⇒ `brentq` always finds it). This collapses the whole system to a single 1-D root-find
per RHS evaluation — consistent with the existing `solve_ph` precedent in `acidbase.py`. Verified
(20-trial random numeric check + a dedicated regression-anchor unit test) to reduce **exactly** to
the old D-28 quadratic when only one carbonyl is active — the isolability proof for prime directive
#3, done algebraically rather than via a toggle. `bound_so2_molar`'s signature changed from
mass-based scalars to a tuple of `(molar_concentration, Kd)` pairs, working natively in mol/L via
the existing `M_ACETALDEHYDE`/`M_PYRUVATE`/`M_ALPHA_KETOGLUTARATE` molar masses (chemistry.py, from
D-49/D-50) — resolving the D-50 calibration-pending mole-vs-mass flag. `free_acetaldehyde` reads
back only acetaldehyde's own bound share from the shared solve, so competing keto-acid pools
measurably reduce acetaldehyde's SO₂ protection — the mechanism by which D-51 addresses the
overshoot.

**New sourced parameters, cross-checked against the existing one.** `K_pyruvate_so2` (5.55e-4
mol/L) and `K_alpha_kg_so2` (1.4e-4 mol/L), both from Burroughs & Sparks (1973), *Sulphite-binding
power of wines and ciders I* — apparent dissociation constants at pH 3.3 for pyruvic acid and
2-ketoglutaric acid respectively. This is the same paper the pre-existing `K_acetaldehyde_so2`
traces to, and its acetaldehyde value (1.5e-6 mol/L) matches this codebase's independently-sourced
`K_acetaldehyde_so2` exactly — a direct sourcing cross-check, not a coincidence to lean on but a
confidence signal. Both new params tier `plausible`, uncertainty bands spanning a secondary review's
looser rounding plus the same pH-drift caveat `K_acetaldehyde_so2`'s band already carries.

**The empirical re-derivation — the honest finding.** Measured end-state total-acetaldehyde-vs-SO₂
increments at 50/100/200 mg/L SO₂ doses (the same dose ladder D-48 used) on the standard 21-day
acceptance run, at the **shipped nominal** D-49/D-50 residuals (30/20 mg/L):

| dose (mg/L) | D-48 (acetaldehyde-only) | D-51 (multi-carbonyl) | field target (0.39·dose) | D-48 overshoot | D-51 overshoot |
|---|---|---|---|---|---|
| 50  | 25.7  | 22.3  | 19.5 | 1.32× | 1.15× |
| 100 | 56.1  | 51.4  | 39.0 | 1.44× | 1.32× |
| 200 | 119.0 | 113.0 | 78.0 | 1.53× | 1.45× |

**D-51 is a real but PARTIAL fix.** Competition genuinely narrows the D-48 overshoot at every dose,
concentrated at the low end where the finite keto-acid capacity isn't yet saturated — but it does
not close it. A sensitivity check pushed both residuals to the **top of their already-sourced
literature uncertainty bands** (pyruvate → 100 mg/L via `k_pyruvate_excretion`=1e-2, α-KG → 70 mg/L
via `k_alpha_kg_excretion`=7e-3 — frozen state verified to land exactly there, not just threshold-
checked) and got 0.86×/1.10×/1.29×: the 50 mg/L point *undershoots* while 200 mg/L is still 1.29×
over. That crossover is the tell — no single scaling of finite-capacity Langmuir competitors can
match a response that stays linear (constant 0.39 mg/mg) across a dose range where binding sites
saturate. This is **structural, not a value not yet found**: more pool mass buys a bigger low-dose
win at the cost of a high-dose miss, it doesn't uniformly close the gap. Per the owner's explicit
guardrail — "do not force-fit beyond the literature-sourced pool ranges" — and advisor concurrence,
**the shipped D-49/D-50 residuals (30/20 mg/L) are unchanged.** The field's 0.39 mg/mg is an
*ensemble* regression over 237 wines with varying carbonyl levels and pH (Jackowetz & Mira de
Orduña 2013)
**[CITATION CORRECTED in D-61: this attribution conflates two different papers. The exact slope
equation `W_acetaldehyde = −4.4 + 0.39·W_tSO₂` (R = 0.837, p < 0.001) comes from Marrufo-Curtido,
Ferreira & Escudero 2022, *Foods* 11(3):476 — a 12-wine forced-oxidation study over a 20–124 mg/L
total-SO₂ range, NOT a within-wine titration. Jackowetz & Mira de Orduña 2013 is the separate
237-wine "Survey of SO₂ binding carbonyls" (Food Control 32(2):687–692), which reports average
binder *concentrations* (acetaldehyde 25/40, pyruvate 14/25, α-KG 74/31 mg/L red/white) — the
correct source for the finished-wine keto-acid *ranges* the D-49/D-50 residuals are anchored to,
but it does not report this slope. The "ensemble regression over 237 wines" phrasing is thus wrong
on both the wine count and the paper; the category-mismatch argument it supports is unaffected —
if anything strengthened, since the true anchor is a 12-wine cross-sectional survey to only 124
mg/L. See D-61.]**; tuning one ferment's pool size to chase it is a category mismatch, not a calibration,
and would trade a documented author-estimate for a fitted number whose only justification is
proximity to a plot — provenance this project ranks below a genuine fit. **This reshapes the task's
own premise** (named "the actual fix for the D-48 overshoot"): the data says D-51 is real, correct,
and load-bearing progress, not a closure — the same "task premise refined by data" shape as D-48
itself. Closing the remaining gap needs a different structure (e.g. a mechanism that scales with
dose rather than a fixed-capacity pool), deferred to a future milestone, not blocking M2.

**A genuine side effect, fixed honestly, not loosened blindly.** The always-on keto-acid pools now
also compete for bisulfite in `test_malolactic.py`'s SO₂-dosed MLF integration test (previously the
equilibrium only knew about acetaldehyde) — adding binding capacity lowers overall free/molecular
SO₂ (~21%→~15% of an 80 mg/L dose), weakening the MLF brake a little further and letting malic
conversion edge just past the halfway mark (~51%, was ~48%). Verified via a standalone run before
editing (not guessed); the test's band and docstring were updated to the measured value with an
explanation, following the same discipline as the D-47/D-48 caveat in `test_post_af_so2_dose_...`.

**Isolability (prime directive #3).** No new state, no new carbon flow — D-51 is a pure readout
generalisation over pools D-49/D-50 already built, exactly like D-28 itself. No `total_carbon`
change needed. The algebraic n=1 reduction (proven, not just tested at default params) is the
isolability guarantee: any ProcessSet lacking the keto-acid pools sees `bound_so2_molar` called
with zero-molar competitor entries, which fall out of the shared solve exactly as if they were
never passed.

**Verification.** `test_so2.py`: fixed the `bound_so2_molar` call sites for the new tuple API, added
a D-51 section (regression-anchor reduction-to-D-28 test, competition-conserves-and-binds-less
test, order-independence/clamping test — all pure algebra) plus a state-level integration test
(keto-acid pools present widen bound SO₂ and free more acetaldehyde). `test_acetaldehyde.py`'s
post-AF-dose test updated for the real consequence that residual keto acids now bind ~34% of a
late SO₂ dose that used to be assumed "nearly all free." `test_malolactic.py`'s SO₂-dosed MLF band
updated per the side effect above. **650 passed** (646 + 4 new D-51 tests, incl. the 5 §2.2
benchmarks), ruff + mypy clean. Validated core untouched; only the SO₂ speciation readout and its
three downstream consumers (`free_acetaldehyde`, the MLF gate, and the speciation dataclass) moved.

## D-52 — MLF v2 refinement: a bounded ethanol/starvation stress multiplier on `MalolacticSenescence`

**Status: IMPLEMENTED 2026-07-07** (654 passed incl. the 5 §2.2 benchmarks, ruff + mypy clean). With
M2 physics beats all complete through D-51, the owner asked for "whichever [MLF v2 deferred item] is
closer to reality," delegating a fidelity judgment among three open candidates: a `BrettSenescence`
twin (D-40 arc), a separate `molecular_so2_death_scale`, or lifting D-41's "environment-free"
ethanol/starvation-modulation deferral.

**The advisor caught a wrong initial pick, verified against the repo before building (the
"discuss disagreements" discipline).** The first-pass read favoured `BrettSenescence` — "Brett only
dies via SO₂ today, so it never declines on its own, which reads as a missing mechanism." An advisor
call reversed this: *Brett's defining real-world trait is persistence* (VBNC survival in barrel/bottle
for years, the textbook "low-and-slow" spoiler), and DECISIONS already says so explicitly (D-40 pt3:
"Without SO₂ (or a rack) Brett persists indefinitely in v1 — an **honest reflection** of how tenacious
a barrel Brett infection is"). A senescence twin would therefore be a fidelity *downgrade*, not a gain.
The advisor's counter-pick — ethanol/starvation modulation of `MalolacticSenescence` — was verified
against the actual functional form in `malolactic.py` before committing: `r_sen` is a *tiny* rate
(~100× below `k_death_mlf`) with no multiplier at all, so a *bounded* stress factor could scale it
without reproducing the D-39 wipeout (which came from multiplying a *large*, full-kill-calibrated rate
by an unbounded `1 − toxicity` ≈ 0.92). Both premises were confirmed by reading the source before any
code changed — the third-option `molecular_so2_death_scale` split was ruled out immediately as pure
parameter architecture with zero fidelity gain.

**The model — two smooth, capped Monod-type stress terms, not a re-run of the Luong wall.**

    r_sen = k_senescence_mlf · X_mlf · arrhenius(T, E_a_death_mlf, T_ref) · stress
    stress = 1 + k_senescence_ethanol_scale·[E/(E+ethanol_tolerance_mlf)]
               + k_senescence_starvation_scale·[K_aa_mlf/(K_aa_mlf+amino_acids)]

Each bracketed term is a Monod-type factor in **[0, 1)** by construction — no clamp needed, C¹ for the
BDF solver — unlike the Luong wall's near-binary "1 at zero stress, 0 at the tolerance wall" shape that
caused D-39's wipeout. `stress` is therefore hard-capped at
`1 + k_senescence_ethanol_scale + k_senescence_starvation_scale` regardless of how far ethanol or
nutrient depletion runs. **Reuses existing concentration scales** rather than adding new ones:
`ethanol_tolerance_mlf` and `K_aa_mlf` (already read by `MalolacticGrowth`/`MalolacticConversion`) are
the two terms' half-saturation points — the same "arrest-scale reused as a death-adjacent scale"
simplification `MalolacticDeath` already makes with `molecular_so2_inhib_mlf`. Only two new
*dimensionless ceiling* parameters are introduced. The starvation term reuses the growth fuel pool
(`amino_acids`) **inverted** (rises as the pool depletes): it is ≈1 (near-max) once amino acids are
exhausted — the D-23 finding places that at ~1.3 d post-pitch regardless of dose — and falls back
whenever autolysis (D-34) refills the pool, so it tracks the real nutrient-refill dynamic already in
the model rather than acting as a flat add-on.

**Magnitude sizing — empirically bounded, not fitted.** `k_senescence_ethanol_scale` = 1.0,
`k_senescence_starvation_scale` = 0.5 (both speculative, author estimates; direction — ethanol/
nutrient stress accelerates O. oeni decline — is sourced from the same Ribereau-Gayon/Bartowsky &
Henschke references `k_senescence_mlf` already cites). Combined ceiling 2.5× the baseline ⇒ a
worst-case half-life of ~23 d (~3.3 weeks) **at T_ref** even at simultaneously saturating ethanol and
full amino-acid exhaustion — verified directly at the RHS level (`test_senescence_no_wipeout_at_
worst_case_stress_at_benchmark_temperature`, E=1e4 g/L, amino_acids=0), never approaching the
~1-week D-39 wipeout regime. **Temperature is a separate stress axis** (advisor-caught: the first
version of this bound implicitly fixed T=20 °C and its name overclaimed "worst case" — a warm cellar
legitimately shortens the half-life further via the shared Arrhenius factor, e.g. to ~10 d at 30 °C;
that is correct physics, not a wipeout regression). Split into two tests: the T_ref half-life bound
above, plus `test_senescence_warm_worst_case_stress_stays_far_below_the_so2_kill`, which proves the
invariant that actually holds at *any* temperature — chronic senescence stays far below the acute
SO₂ kill because both share `arrhenius(T, E_a_death_mlf, T_ref)`, so their ratio is
temperature-invariant by construction (verified numerically at 30 °C). At typical post-AF dry-wine
conditions (E≈100–130 g/L, amino_acids≈0, T_ref) the measured stress factor is ≈2.0×, giving a ~29 d
effective half-life — comfortably "weeks," faster than D-41's flat ~58 d but nowhere near
catastrophic.

**Owner-flagged open question (advisor-raised, not resolved here — a fidelity call, not a
provenance footnote).** D-41 calibrated `k_senescence_mlf=5e-4` so a *typical* unsulfited wine loses
~half its O. oeni over ~2 months. But in an ordinary pitched run the starvation term saturates almost
immediately (amino acids ≈0 by ~day 1.3, the D-23 finding) and the ethanol term adds ~+0.5 post-AF,
so `stress≈2×` is close to the *typical* case, not the exception — meaning the *typical*-wine
half-life this decision actually produces is now ~29 d, not D-41's ~2 months. `k_senescence_mlf` was
left unchanged (5e-4), so it has silently shifted from "the typical-wine rate" to "a benign floor
that rarely applies in practice." Whether the literature's *typical* dry-unsulfited decline is closer
to ~2 months (in which case `k_senescence_mlf` should drop to ~2.5e-4 so the *typical*-stress case
re-anchors on D-41's original target) or closer to weeks-to-a-month (in which case the shipped value
is right or even conservative) is an empirical/owner judgment this decision does not make — flagged
for the owner, not decided here, following the same "re-derive, don't inherit" discipline D-48/D-50/
D-51 already established for calibration questions like this one.

**Isolability + performance preserved.** Still reads **no SO₂ and no pH** — `E`/`amino_acids` are read
directly off state, no equilibrium solve — so the Process remains strictly cheaper than the SO₂ kill,
exactly as D-41 built it. `touches` unchanged (`X_mlf`, `X_mlf_dead`); the carbon/nitrogen-neutral
transfer needs no new conservation code (D-13/D-39 pattern, unchanged). Pitch-gated at compile, not
amino-acid-gated (unchanged from D-41).

**A genuine, honestly-measured side effect on MLF-derived diacetyl clearing (the D-51 discipline
reused).** `OenococcusDiacetylReduction`'s lees-contact clearing scales with viable `X_mlf`; faster
senescence late in a long (30 d) run leaves less bacterial reductase around, so
`test_headline_citrate_lifts_and_then_clears_diacetyl`'s final/peak diacetyl ratio rose from
comfortably under its old 0.85 threshold to a measured **0.861** (X_mlf retains only ~0.49× its dose by
day 30, vs a slower D-41-only decline). This is real and expected — the more realistic senescence
means less bacteria on the lees to clear diacetyl late in the wine's life — not a bug; the test band
was widened to 0.90 with the measured value and explanation recorded, following the exact discipline
D-51 used for its own MLF SO₂-gate side effect (verify first, then band to what's actually measured,
never loosen blindly).

**Verification.** `tests/test_malolactic.py`: the single D-41-era `test_senescence_is_environment_
free`, which pinned the now-superseded "environment-free" invariant, was split into five tests rather
than deleted — `test_senescence_is_so2_independent` (SO₂ independence retained, unchanged),
`test_senescence_ethanol_stress_is_bounded` (monotone rise with ethanol, ratio against the
zero-stress floor strictly below the design ceiling), `test_senescence_starvation_stress_tracks_
amino_acid_depletion` (starved > replete), `test_senescence_no_wipeout_at_worst_case_stress_at_
benchmark_temperature` (the empirical wipeout guard at T_ref, half-life > 2 weeks), and — added after
an advisor pass caught that the wipeout guard's name overclaimed "worst case" while implicitly fixing
T=20 °C — `test_senescence_warm_worst_case_stress_stays_far_below_the_so2_kill` (proves the
temperature-invariant ratio against the acute SO₂ kill instead, verified at 30 °C) — a net +4 tests in
that file. Two more existing tests were updated in place, not weakened: `test_senescence_needs_no_ph_
solve`'s reads-tuple pin was extended for the four new params, and the integration-level `test_so2_
crashes_bacteria_over_the_slow_senescence_baseline` had its no-SO₂ decline bands re-measured and
tightened around the new values (day-21 ratio ~0.608 vs D-41's ~0.71; day-6 pre-dose ratio ~0.875 vs
~0.95) rather than loosened past a threshold. `tests/test_mlf_diacetyl.py`'s headline clearing test
band was widened per the side effect above (measured, not blind). **654 passed** (650 + 4 net new
tests, incl. the 5 §2.2 benchmarks), ruff + mypy clean.

## D-53 — Correction: `k_senescence_mlf` magnitude was wrong by ~50×, per real-wine literature

**Status: IMPLEMENTED 2026-07-07** (654 passed incl. the 5 §2.2 benchmarks, ruff + mypy clean). D-52
delegated a follow-up fidelity question to the owner ("re-anchor `k_senescence_mlf` to compensate for
typical stress, or leave it?") rather than deciding unilaterally. The owner asked for research before
deciding — this is that research's outcome, and it overturned the question's own premise.

**The deep-research finding (5 search angles, 22 sources fetched, 25 claims adversarially verified
3-vote).** Real, finished, unsulfited (SO₂-free) wine shows **no detectable spontaneous decline** in
O. oeni populations for 3–5 months post-MLF:

- **Windholtz, Miot-Sertier, Maupeu et al. 2025**, *OENO One* 59(3), doi:10.20870/
  oeno-one.2025.59.3.9346 — real Bordeaux red wine, 6 SO₂-management modalities tracked vatting →
  bottling (5 months). SO₂-free modalities: "high and stable population levels of around 10⁵
  CFU/mL" from end-of-MLF through 5 months.
- **Millet 2001** (Univ. Bordeaux 2 doctoral thesis, cited within Windholtz et al.) — 3 Bordeaux
  varieties in oak barrels, 0/30/50 mg/L SO₂ over 3 months. At 0 mg/L SO₂: population "maintained at
  around 10⁶ CFU/mL," even at pH 3.75–3.95; only 50 mg/L SO₂ was sufficient to inhibit it.
- **By contrast, Kioroglou, Mas & Portillo 2020**, *Frontiers in Microbiology*, doi:10.3389/
  fmicb.2020.562560 — the steep decline (10⁵–10⁶ → 10³–10⁴ CFU/mL by 3 months, undetectable by 12)
  is documented **only in wines that received SO₂**. That decline is `MalolacticDeath`'s (D-39)
  territory, not spontaneous senescence — a category error in the original D-41 framing.
- Acute ethanol/pH-shock mechanistic studies (da Silveira et al. 2002, doi:10.1128/
  aem.68.12.6087-6093.2002; Bastard et al. 2016, doi:10.3389/fmicb.2016.00613) confirm the *direction*
  (ethanol damages the membrane, worse at low pH) but operate on minutes-to-4-hours timescales under
  artificially harsh conditions (12–16% ethanol, pH 3.2) — several tempting extrapolations from this
  acute data to the weeks/months spontaneous-decline question were explicitly checked and **refuted**
  (0-3 adversarial votes) as overreach.

No survived source measures an actual first-order decay constant or CFU curve for spontaneous decline
beyond ~5 months — the evidence base is "no significant decline detected within the observed window,"
an upper bound, not a fitted point estimate.

**Diagnosis: D-41's original citations were misread, not wrong on their face.** Ribereau-Gayon
(Handbook of Enology) and Bartowsky & Henschke 2004 support the *general winemaking practice* that
"SO₂ is needed to reliably control spoilage LAB" — true, and the basis of `MalolacticDeath`. D-41
over-read that into a *specific* "O. oeni spontaneously declines over weeks-to-months without
intervention" claim, which the direct CFU evidence above does not support. The mistake propagated
into D-52's calibration target ("typical wine loses half its O. oeni over ~2 months") without being
independently re-checked — exactly the failure mode the project's "re-derive, don't inherit" discipline
(D-48/D-50/D-51) exists to catch, and this time it wasn't caught until the owner asked for it.

**The fix — magnitude only, mechanism untouched.** `k_senescence_mlf`: 5.0e-4 → **1.0e-5** (a round
number per advisor guidance — the data gives an upper bound, not a precision target; any value keeping
decline within CFU-measurement noise over ~5 months at typical D-52 stress is equally faithful to the
evidence). At `stress=1`, half-life is now ~2888 d (~7.9 y); at D-52's typical post-AF stress (~2×)
it's ~3.96 y; even D-52's worst-case combined-stress ceiling (2.5×, unchanged) gives ~3.16 y — all far
beyond the 3–5 month window the literature actually measured as "stable," a deliberately conservative
choice given no source pins the true value more precisely. D-52's stress-multiplier *mechanism*
(bounded Monod-type ethanol/starvation terms, the wipeout-avoidance structure) is **completely
unchanged** — only the baseline rate it scales was corrected. `k_senescence_ethanol_scale` (1.0) and
`k_senescence_starvation_scale` (0.5) are untouched as dimensionless ceilings; their provenance notes
were updated to point at the new baseline's resulting half-lives.

**Honest consequence, surfaced to the owner before proceeding (not buried in a re-band).** At this
magnitude, D-52's stress multiplier is now **empirically inert on every timescale this model
simulates** — even worst-case combined stress gives a multi-year half-life, invisible in any real run
(the model's longest integration test is 30 days). This is the *correct* closest-to-reality outcome
(spontaneous senescence genuinely is negligible at these timescales), but it changes what D-52 "does."
Two structural options were put to the owner — (a) keep the stress-multiplier structure as a
documented slow long-tail mechanism (decline beyond 5 months is genuinely unmeasured, so the mechanism
remains defensible even though it's unobservable at simulated timescales), or (b) simplify by
stripping the machinery back to D-41's flat-rate form, since two extra parameters now model something
no test can see. **Owner chose (a), the least-churn default** — keep + recalibrate.

**Test consequence — an assertion flip, not a re-band (advisor-caught: rerunning-and-rebanding would
have been the wrong instinct here).** The integration test asserting a *measurable* decline
(`test_so2_crashes_bacteria_over_the_slow_senescence_baseline`, D-52's day-21 ratio ~0.608 / day-6
~0.875) now directly contradicts the corrected evidence — those numbers described the wrong physics,
not just an imprecise band. Renamed to `test_so2_crashes_bacteria_over_the_near_stable_senescence_
baseline` and flipped to assert *near-stability* (measured day-21 ratio ~0.990, day-6 ~0.997),
still checking a nonzero (if tiny) monotone decline exists structurally, plus the SO₂ crash
mechanism is unaffected. `test_mlf_diacetyl.py`'s headline clearing test comment, which attributed
its ~0.861 final/peak ratio to D-52's faster senescence, was corrected: with X_mlf now ~98.6% viable
at day 30, the measured ratio reverts to ~0.742 (closer to D-41's original clean-clearing picture),
and the band tightened from 0.90 back to 0.80 to match. All other D-52 RHS-level tests (ethanol-bound,
starvation-tracks, no-wipeout, warm-vs-kill-ratio, SO₂-independence) are unaffected — they test
direction/ratios that are magnitude-independent of `k_senescence_mlf`. **654 passed** (unchanged
count — one test renamed and reassigned, none added/removed), ruff + mypy clean.

**Method beat worth remembering — the third `advisor()` call in this arc, and the value of asking
before assuming.** D-52 shipped a plausible-looking calibration (owner-delegated "closer to reality")
that turned out to rest on a misread citation. The owner declining to pick a number and asking for
research first — rather than accepting either of the two options originally offered — is what caught
it: neither "re-anchor to ~2 months" nor "leave it at ~29 days" was defensible once real CFU data was
checked. A third advisor call (post-research) then caught that the fix wasn't a simple re-band but an
assertion flip, and surfaced the "D-52 is now inert" honesty point before it could be silently
absorbed. Three advisor calls across one feature arc, each catching something the previous pass
missed — the discipline compounds.

## D-54 — POF v2 pt1: `E_a_pof` temperature dependence, direction-checked before calibrated

**Status: IMPLEMENTED 2026-07-07** (all Brett/POF tests green, ruff + mypy clean). With M2 physics
complete through D-53, the owner picked "POF v2" (temperature dependence for conversion efficiency,
plus splitting the lumped vinylphenol/vinylguaiacol pool) as the next work. The two pieces are
independent and sequenced separately (advisor guidance, matching the per-D-record discipline); this
entry is pt1 (`E_a_pof`) only — the pool split is a separate, larger decision (D-55+).

**The crux worth remembering — cloning a nearby `E_a_*` precedent would have picked the WRONG
ordering.** `YeastPOFDecarboxylation` was deliberately temperature-flat at D-40 pt4 ("no pt4 behaviour
needs POF's intrinsic direction"). The naive v2 move — add `arrhenius(T, E_a_pof)` to the rate,
magnitude cloned from a neighbouring decarboxylase `E_a` (e.g. `E_a_decarb` = 90 kJ/mol) — was
**caught by `advisor()` before any code was written**: `YeastPOFDecarboxylation`'s rate is
**flux-coupled** (`r ∝ fermentative_flux_shape`, which itself rides `E_a_uptake`), and this codebase
already has a named framework for exactly this interaction — the D-19 "KEY ORDERING CONSTRAINT"
governing `E_a_esters`/`E_a_fusels`: a flux-coupled byproduct's **net** (time-integrated-to-dryness)
total scales as `exp(-((E_a_byproduct − E_a_uptake)/R)(1/T−1/T_ref))`, because a warmer ferment also
finishes *faster*, shrinking the production window. So the net finished-wine direction is set by
`E_a_pof` **relative to** `E_a_uptake` (55,100 J/mol), not by `E_a_pof` in isolation — cloning a
positive value blind would have picked a direction by accident, the exact D-53 failure mode
(magnitude-by-analogy, direction unchecked) one decision later.

**Research resolved the direction, not an owner preference (the D-53 discipline applied
prospectively).** WebSearch found: (1) Edlin et al. 1998 (hydroxycinnamate decarboxylase purified
from *Brettanomyces anomalus* — the same enzyme family) puts the enzyme's own thermal optimum at
40 °C, well above any wine/beer ferment temperature, supporting a genuine positive intrinsic
`E_a_pof`; (2) brewing practice on this *exact* enzyme (Pad1/Fdc1, the same POF+ trait, well
corroborated across independent wheat-beer/Weizen fermentation-temperature sources — Brewing Science
Institute, Northern Brewer, Brülosophy trials) is unambiguous that **cooler fermentation retains more
clove/4-vinylguaiacol character; warmer fermentation favours esters over phenolics**. Net conversion
therefore **falls** with warmer temperature — the *opposite* ordering from esters/fusels (which need
`E_a > E_a_uptake` to *rise* with T), because the sourced real-world direction here is the reverse of
theirs. This literature is beer/Weizen-sourced, not wine-specific; extended to this model's wine POF+
yeast by the same enzyme-identity argument the module's own docstring already makes for `k_pof_decarb`.

**The fix.** `E_a_pof` = 25,000 J/mol (uncertainty 10,000–40,000, chosen so even the high end stays
below `E_a_uptake`'s own low uncertainty bound of 47,000 — the sourced direction must survive the
joint uncertainty band, not just the point estimate), embedded in `YeastPOFDecarboxylation.derivatives`
via `arrhenius_factor(T, E_a_pof, T_ref)` (the same `BrettDeath`/`AcetolactateDecarboxylation` embedded-
call idiom, not a `ProcessSet` `RateModifier`). **Honest continuity note:** v1's implicit `E_a_pof = 0`
already had `0 < E_a_uptake`, so the *emergent direction* was accidentally already correct before this
change — v2 replaces an implicit "enzyme rate is T-invariant" placeholder (a stronger, more obviously
false claim) with a genuine sourced-direction intrinsic term, sized to preserve and reinforce that same
direction rather than risk reversing it.

**Two new tests, split to isolate the two effects.** `test_pof_own_rate_rises_with_warmth` calls
`.derivatives()` directly at fixed flux/precursor, isolating the raw Arrhenius direction (pins
`E_a_pof > 0`). `test_pof_net_conversion_falls_with_warmer_fermentation` runs full POF+ (no Brett)
scenarios to dryness at 12 °C vs 28 °C over a shared 60-day window and compares frozen post-dryness
`vinylphenols` totals — empirically confirming the *net* direction the algebra predicts, not just
asserting it. Both pass. `test_pof_decarboxylation_stoichiometry_and_touches` (which calls
`.derivatives()` at the default `T = T_ref`, where `arrhenius_factor = 1` exactly) is numerically
unaffected — no re-pin needed. 39/39 `test_brett.py` green, ruff + mypy clean.

**Method beat worth remembering.** One `advisor()` call, before any code, caught a wrong-by-default
generalization (treat every new `E_a_*` like the nearest existing one) that this codebase's *own*
D-19 framework already had the tools to refute — the miss would have been not reading the codebase's
existing ordering-constraint machinery closely enough, the same category of miss D-52's pass 1 caught
(checking prior decisions before reasoning from first principles). Research then resolved the
direction empirically rather than picking a plausible number, continuing the D-53 discipline forward
rather than only applying it in hindsight.

## D-55 — POF v2 pt2: splitting the lumped vinylphenol/vinylguaiacol pool into a real ferulic branch

**Status: IMPLEMENTED 2026-07-07** (46 `test_brett.py` tests green, full suite + ruff + mypy clean,
across 3 commits). D-40's original design deliberately lumped both the precursor pair (p-coumaric +
ferulic hydroxycinnamic acids, booked as p-coumaric) and the product pair (4-vinylphenol +
4-vinylguaiacol, and downstream 4-ethylphenol + 4-ethylguaiacol) into three single pools
(`hydroxycinnamics`/`vinylphenols`/`ethylphenols`). The owner chose to split it, "your call" on how
to break the work into pieces — this closes that arc.

**The scope-collapsing fact, caught before any code was written.** The three candidate designs
initially considered were: (a) split only the product pools by a fixed ratio on the existing single
precursor, (b) keep the lumped ratio but relabel a fraction as "vinylguaiacol", or (c) build a
genuine second precursor pool. Options (a)/(b) collapse immediately on inspection of the actual
molar masses already in the codebase: `hydroxycinnamics` is *literally* booked as p-coumaric acid
(9 carbons; `M_P_COUMARIC` is used for every unit conversion), and ferulic acid is a **different,
10-carbon molecule** whose decarboxylation is `10 C → 9 C (vinylguaiacol) + 1 C (CO2)`, not
`9 C → 9 C + 0 C`. A 9-carbon precursor cannot yield a 9-carbon product plus a CO2 molecule without
manufacturing a carbon out of nothing — so any fixed-ratio split of the *existing* pool's output
breaks carbon closure by construction, the one invariant this codebase enforces as a test
(`assert_conserved`), not a suggestion. Only a genuine second precursor pool is carbon-exact and
species-faithful. This left one real design, not three co-equal options — surfaced to the owner as
a binary (full split vs. document-as-limit, the D-51 precedent) rather than offering fake choices.

**The split, mechanically (3 commits, each independently green).**

1. **Scaffolding (chemistry + state, no behaviour change).** New species `M_FERULIC` (C10H10O4,
   194.19 g/mol), `M_VINYLGUAIACOL` (C9H10O2, 150.18 g/mol), `M_ETHYLGUAIACOL` (C9H12O2, 152.19
   g/mol) in `chemistry.py`, verified against known real molar masses and carbon closure
   (`10 = 9 + 1`) by direct computation, not just formula arithmetic. Three new wine-only state
   slots (`ferulic_acid`/`vinylguaiacols`/`ethylguaiacols`, schema size 36 → 39), their
   `total_carbon` weighting, and the `test_media.py` schema-shape assertions updated in the same
   commit (tightly coupled to the slot count, unlike the later Process-behaviour tests).
2. **Decarboxylation branch.** `BrettDecarboxylation` and `YeastPOFDecarboxylation` both gained a
   second, independent branch (`ferulic_acid → vinylguaiacols + CO2`) via a shared
   `_decarboxylation_branch` helper (factored out to avoid 4×-duplicating the Monod/molar-mass
   arithmetic across 2 Processes × 2 branches). Both branches share the *same* catalyst/gate
   (`X_brett · gate` for Brett, `flux · arrhenius(T, E_a_pof)` for POF) — the enzyme and its
   environmental sensitivity don't depend on which substrate it happens to be processing.
3. **Reduction branch + scenario wiring.** `BrettVinylphenolReduction` gained the
   `vinylguaiacols → ethylguaiacols` branch via a shared `_reduction_branch` helper, and
   `ferulic_acid_gpl` was wired into the scenario compiler's `_ALLOWED_KEYS`/`_wine_initial` so
   scenarios can dose the new precursor exactly like `hydroxycinnamic_gpl`.

**Relative kinetics are sourced, not cloned — the same discipline D-54 established.** Edlin et al.
1998 (*Appl. Microbiol. Biotechnol.* 49:511-517) purified a hydroxycinnamate decarboxylase from
*Brettanomyces anomalus* (the same enzyme family as this model's decarboxylase) and report **paired**
Vmax/Km for both substrates in the *same* assay: Vmax 13,494 (ferulic) vs 22,256 (p-coumaric)
nmol/min/mg; Km 1.15 (ferulic) vs 1.55 (p-coumaric) mM. Those ratios (~0.606× rate, ~0.742×
half-saturation) are real, paired, sourced data — applied to this model's own already-speculative
absolute `k_brett_decarb`/`K_hydroxycinnamic` scale, so the *ratio* between the two branches carries
real evidentiary weight even though the absolute magnitude it scales remains an author estimate (the
same "ratio sourced, absolute speculative" pattern D-49/D-50's keto-acid pools used). New params:
`k_brett_decarb_ferulic`, `K_hydroxycinnamic_ferulic`, `k_pof_decarb_ferulic` (all ratio-derived,
uncertainty bands scaled by the same ratio as their point estimates). `E_a_pof` (D-54) and the Brett
environmental gate are **reused as-is** for the ferulic branch (same enzyme, same organism — no new
temperature/SO₂ parameters needed).

**One honest gap, surfaced rather than papered over.** Tchobanov et al. 2008 (*FEMS Microbiol.
Lett.* 284:213-217) directly confirm Brett's vinylphenol reductase acts on **both** 4-vinylguaiacol
and 4-vinylphenol — upgrading what D-40 had left as an unstated assumption to a sourced fact. But
that paper reports absolute kinetics for vinylguaiacol only (Km 0.14 mM, Vmax 1900 U/mg), with no
paired p-coumaric-branch number to derive a relative rate the way Edlin et al. 1998 allowed for the
decarboxylase. So `k_brett_reduction` is **reused unchanged** for both branches — a documented
simplification (enzyme identity sourced; relative rate not), distinct in kind from the
decarboxylase branches' sourced ratio. Similarly, no clean paired p-coumaric:ferulic *concentration*
ratio was found in must/wine literature, so `ferulic_acid` is dosed independently per-scenario
(default 0, like `hydroxycinnamics` itself) rather than forced to any fixed ratio of the p-coumaric
dose.

**New tests, mirroring the existing per-branch pattern at every level:** per-Process stoichiometry
and carbon-closure for the ferulic branch alone and composed with the p-coumaric branch (both
Brett's decarboxylase/reductase and POF+'s decarboxylase), a POF+-no-Brett stranding test for
vinylguaiacols (mirroring the pt4 headline), and a full-scenario end-to-end carbon-closure test
dosing both `hydroxycinnamic_gpl` and `ferulic_acid_gpl` together through the whole
compile→decarboxylate→reduce pipeline (a wiring-level check the per-Process unit tests can't catch
— e.g. a typo in the new scenario dosing key, or a slot the reduction step forgot to drain). All
existing `touches` assertions updated to the grown tuples (5 slots per decarboxylase Process, 4 per
the reductase). `test_wine_schema_has_single_sugar_slot`'s slot count and `WINE_BRETT_SLOTS` tuple
updated in the scaffolding commit. 44 `test_brett.py` tests total (8 net new: 2 from D-54's
`E_a_pof`, 6 from this split), 5 §2.2 benchmarks unaffected (undosed default runs stay
byte-for-byte the validated core — `ferulic_acid`/`vinylguaiacols`/`ethylguaiacols` default to 0
and no benchmark doses them).

**Closes the last D-40 pt4 deferral.** Both "POF v2" items from the deferred list — `E_a_pof`
temperature dependence (D-54) and the vinylguaiacol/vinylphenol split (D-55) — are now done.

## D-56 — First independent-data validation attempt (Varela et al. 2004): the model runs 2–4× too fast, diagnosed not fixed

**Status: DOCUMENTED 2026-07-07** (no core-code change; a real-data regression benchmark added). With
M2 physics complete through D-55, the owner picked "validation against real data" as the next
direction. This is the first time the project checked a core M1 output against a dataset genuinely
independent of the papers its own parameters were fit to.

**The independence discipline, fixed before any comparison ran.** Per `CLAUDE.md`'s tier definition
(`VALIDATED` = checked against independent measured data, `combine()`/`Tier.tiers.py`) and the
project's own prior honesty about it (D-C/D-46: reproducing Coleman, Fish & Block 2007 is a
consistency check, not validation, since that paper *is* where `mu_max`/`K_n`/`q_sugar_max`/
`K_sugar_uptake`/`biomass_N_yield_log_*` etc. come from) — any candidate dataset first had to pass the
test "were any of the model's parameters derived from this dataset or its source?" A deep-research
sweep (104 sub-agents, 21 sources, 25 adversarially-verified claims) surfaced two genuinely
independent wine candidates: Varela, Pizarro & Agosin 2004 (*Appl. Environ. Microbiol.* 70(6):3392-
3400, doi:10.1128/AEM.70.6.3392-3400.2004, Pontificia Universidad Católica de Chile — no author/lab
overlap with Coleman or any other cited source) and Palma et al. 2012 (Lisbon). Both turned out to be
figure-only for their raw time-series (direct WebFetch of both PMC articles found no numeric table for
the actual sugar/biomass/ethanol curves — Palma's "Table 1" is glucose-transport Km/Vmax, not a
fermentation time series), but **Varela's Table 1 gives exact endpoint values with real replicate-
based uncertainty** (3 independent experiments each): 300 mg N/L (well-fed) reaches dryness in
170 ± 12 h with 5.8 ± 0.1 g/L final biomass; 50 mg N/L (severely N-deficient) takes 700 ± 10 h with
1.5 ± 0.1 g/L biomass. Owner chose the endpoint check over digitizing figures — no digitization-error
uncertainty to carry, and a genuinely out-of-sample test (50 mg N/L sits below Coleman's fitted
70–350 mg N/L range).

**Setup: same strain, so this is a clean two-lab comparison.** Varela used *S. cerevisiae* EC1118
(Prise de Mousse) — the *same* strain `wine_generic.yaml`'s header already declares this model is
calibrated on (Premier Cuvee/EC-1118), removing the strain-difference hypothesis. Isothermal 28°C
(inside Coleman's 11–35°C fit range), synthetic must (120 g/L glucose + 120 g/L fructose = 240 g/L,
100% fermentable — `must_fermentable_fraction=0.93` is a real-grape-must correction that under-loads
a pure-sugar must by ~7%, noted but not fixed for this probe), 10⁶ cells/mL inoculum (converted to
`pitch_gpl≈0.018` via the standard ~18 pg/cell dry-weight figure — an order-of-magnitude conversion,
not exact, but the plausible range can't explain a 2–4× gap on its own).

**The result: the model runs 2–4× too fast, worse at low N — and this decomposes into THREE distinct,
separable findings, not one bug.**

1. **A uniform ~2× gap present even in-range** (N=300, both T and N inside Coleman's fit window): model
   83 h vs. Varela's measured 170 ± 12 h. This traces to an **already-documented** M1 simplification —
   `q_sugar_max`'s own provenance note (`wine_generic.yaml:73`) says the rate is "applied to TOTAL
   biomass with no active/inactive split, whereas Coleman's active X_A declines late — so M1 will
   over-catalyse the tail." Not a new discovery; an independent dataset confirming a known caveat has
   real, measurable cost.
2. **An additional ~2× gap specific to severe N-deficiency** (N=50, below Coleman's 70–350 mg N/L
   floor — genuine extrapolation): model 176 h vs. Varela's 700 ± 10 h, i.e. ~4× total vs. Varela.
   **Isolated cleanly via a biomass-hours integral:** because `K_sugar_uptake` (10.3 g/L) is tiny next
   to S (~100s g/L) for most of the run, sugar consumed ≈ `q_sugar_max_eff(T) · ∫X dt`, and since S₀ is
   identical (240 g/L) for both N conditions, the model's own structure forces `∫X dt` to dryness to be
   *nearly identical* between them (183.6 vs 183.7 g·h/L, confirmed numerically) — meaning **duration
   in this model is set entirely by how fast biomass X(t) builds**, i.e. by nitrogen-limited growth
   kinetics, not sugar-uptake rate. The model's N50/N300 duration ratio is 2.12×; Varela's real ratio is
   4.12×. Literature-consistent explanation (Bisson's stuck/sluggish-fermentation review, via search
   snippet — full-text PDF extraction failed, so this is *not* yet a citable primary source): hexose-
   transporter turnover/degradation accelerates under nitrogen deficiency, reducing per-cell
   fermentative capacity beyond what a lower biomass ceiling alone predicts — a mechanism absent from
   this model.
3. **A separate, genuine cross-study biomass-yield gap at N=300** (model 42% low: `Y_X/N` computed from
   the Coleman regression at YAN=300 gives 11.2 g cell/g N vs. Varela's implied 19.3 g/g; at N=50 the
   model is much closer, 27.7 vs. 30 g/g — only 8% low). Confirmed to be the model behaving exactly as
   designed (`X_max = Y_X/N · N0` reproduces the simulated peak biomass to 3 significant figures at both
   N levels), not a bug — a real difference between Coleman's Chardonnay-must lab strain-N-yield
   relationship and Varela's synthetic-must EC1118 fermentations. **Explicitly not to be "fixed" by
   raising the model's biomass** — more biomass would make the duration mismatch *worse*, not better,
   since the model already over-catalyses (finding 1).

**A single-term fix was prototyped and disproved — the firewall that stopped further tuning.** Per
Bisson's mechanism, a candidate fix is an ethanol-driven, nitrogen-gated decline in effective
`q_sugar_max` (`q_eff = q_sugar_max · exp(-k_decay · severity(N₀) · E)`, `severity = K_sev/(K_sev+N₀)`).
Monkeypatch-prototyped (no core files touched) and swept over a parameter grid: **no single-term fit
gets within 15% of both targets simultaneously** (best combined relative error ~60%), and there is a
structural reason, not just a sweep gap — narrowing `K_sev` to differentiate N=300 from N=50 leaves
N=300 under-corrected; widening `K_sev` to fix N=300's magnitude collapses the N-differentiation needed
to stretch the ratio to 4.12×. **At least two distinct effects are needed**, confirming findings 1 and
2 above are mechanistically separate, not one bug wearing two faces. The sweep was stopped there
deliberately: Varela is the project's only independent wine dataset, and it can only be a *validation*
set if it is never used as a *calibration* set — tuning ≥2 free parameters against 2 data points is a
guaranteed fit that proves nothing and burns the one check the project has. **If a two-mechanism build
is ever undertaken, the parameters must be sourced independently from Bisson's primary literature (the
review's cited 3.6×/10× specific-uptake fold-changes, transporter turnover rates) — not fit to
Varela — and then checked against a held-out condition or a third dataset**, preserving the
validation/calibration firewall. Not started; a candidate future task, not scheduled.

**What shipped: a real-data regression benchmark, not a physics fix.** `tests/benchmarks/
test_validation_varela2004.py` runs both conditions and asserts the model's *current* characterized
behavior (duration + biomass at each N level, and the gap ratio to Varela's measured values) stays
within the diagnosed bands — so a future change that silently widens *or* closes the gap gets caught
either way, and the honest "how far off are we" number stays live in the suite instead of decaying into
a stale doc comment. No `BENCHMARKS`/`ReferenceSeries` entry (Varela's data is two endpoints with
replicate uncertainty, not a fittable time series — the existing `compare_series` RMSE machinery
doesn't apply; a plain benchmark test in the `test_milestone1.py` "realism regression guard" style
fits better). No tier promotion: none of `growth.py`/`uptake.py`/`inhibition.py`/`arrhenius.py` moved
off `PLAUSIBLE` — matching an aggregate endpoint doesn't license per-parameter tier bumps (non-
identifiability: many different parameter combinations could reproduce the same duration/biomass pair),
and separately, `ProcessSet.tier_of`'s honest `param_tiers` path already floors wine `S`/`X` at
`SPECULATIVE` today via `K_s`/`K_repression`/`Y_byproduct_sugar` (all `speculative`, "author estimate"
placeholders) regardless of any Process-class tier — so an end-to-end `VALIDATED` output was never
reachable from this comparison alone, independent of the fit-quality question. §2.2 benchmarks
untouched; undosed default runs unaffected (this is a new, additional scenario, not a change to any
existing one).

**Method beat worth remembering: three advisor() passes, each catching a different failure mode in
real time.** Pass 1 (before running anything) caught that promoting a Process's tier is not what
"validation" mechanically does in this codebase — traced `ProcessSet.tier_of`'s actual `min()`-combine
behavior before writing a line of benchmark code. Pass 2 (after the first probe run) caught a
confounded 20°C-vs-28°C comparison that was about to misattribute the whole gap to Arrhenius
temperature extrapolation — the only valid same-temperature comparison was 28°C vs. 28°C, and D-14
already established the model reproduces Coleman's own 11–35°C shape line-for-line, which the wrong
framing would have silently contradicted. Pass 3 (mid-sweep) caught the validation/calibration firewall
before a "good enough" two-parameter fit could be mistaken for a validated mechanism. Each catch was a
premise correction the transcript shows in full, not a rubber stamp.

## D-57 — Correction: D-56 finding 1 was misdiagnosed (a stale note); the real bug was `k_prime_d`'s missing quadratic temperature scaling, fixed and sourced from Coleman's own regression

**Status: LANDED 2026-07-07.** Owner picked up D-56's "two-mechanism uptake-decline build" as the
next task: source Bisson transporter-turnover parameters independently to fix the extra ~2× gap at
severe nitrogen deficiency. Before building anything, the mechanism-1 premise ("`q_sugar_max`
applies to TOTAL biomass with no active/inactive split") was checked against the *current* code
rather than taken from the D-56 record — it did not hold up, and that check reshaped the entire task.

**Mechanism 1 does not exist as a fixable gap — it was already fixed in D-13, three commits before
the note that "diagnosed" it was even written.** `wine_generic.yaml`'s `q_sugar_max` caveat ("M1
applies this rate to TOTAL biomass with no active/inactive split") was added in commit `5da7725`
(D-12); `EthanolInactivation` — which splits `X` (viable) from `X_dead` and is what both
`GrowthNitrogenLimited` and `SugarUptakeToEthanolCO2` already read exclusively — landed in the very
next commit, `c244ae6` (D-13). Structural check: this model's `dX/dt = mu·X − k'_d·E·X` is
byte-for-byte Coleman's own eq. 2 for his active pool `X_A` (`test_coleman_reconstruction.py`
already proves line-for-line agreement); Coleman's separate "total biomass" `x` (his eq. 1, no
death term) is used by *nothing* in his own eqs. 1–8, so the model correctly never tracks it
either. The note describes a pre-D-13 model that no longer exists. Advisor-caught before any
Bisson literature search was spent chasing a mechanism that was never missing.

**The real, sourced bug: `k_prime_d` — Coleman's death-rate constant and the one parameter his fit
found QUADRATIC in temperature — shipped with no temperature modifier at all.** The D-12 provenance
note says so explicitly: "M1 is isothermal at 20 C so no Arrhenius modifier is attached (the
quadratic does not reduce to a single activation energy anyway)." Correct scoping *for M1* — but M2
added non-isothermal scenarios (temperature ramps, D-35/36) without anyone revisiting this, so every
non-20 C wine/beer run since has driven growth and uptake with Arrhenius scaling while leaving death
frozen at the 20 C rate. Decisive check (advisor-directed): integrate Coleman's own eqs. 1–8
(already sitting in `test_coleman_reconstruction.py`) at Varela's exact 28 C/S0/pitch/N0 inputs, with
`k_prime_d` frozen at its 20 C value exactly as the engine does — this reproduces the engine's
numbers almost exactly (N=300: 78.5 h vs engine 83 h; N=50: 164.5 h vs engine 176 h), while the
*correctly* temperature-scaled Coleman reference gives N=300: 84.5 h (barely different — short run,
death is a minor contributor by dryness) and N=50: 283 h (much longer — 40+ days of compounding
ethanol exposure at the wrong, too-gentle death rate). The asymmetry (fermentation-*driving*
processes correctly accelerate at 28 C, the fermentation-*braking* one doesn't) is exactly why D-56
read the gap as "worse at low N": it isn't a missing nitrogen-transporter mechanism, it's a death
rate quietly stuck at the wrong temperature on any long run.

**Fix: `ColemanQuadraticDeathTemperature`, a new `RateModifier` implementing the regression
directly, not an Arrhenius approximation.** `arrhenius.py` already had a per-rate `E_a` form
(D-11), but the D-12 note is explicit that a single activation energy cannot reproduce a quadratic's
curvature — so this modifier evaluates Coleman's `ln(k'_d) = a0 + a1·T_C + a2·T_C²` directly,
normalised to `T_ref` so the intercept `a0` cancels (`k_prime_d` itself already IS the T_ref-evaluated
value): `factor(T) = exp(a1·(T_C−T_ref_C) + a2·(T_C²−T_ref_C²))`, exactly 1 at `T = T_ref` (same
reference-anchored pattern as `ArrheniusTemperature`, D-11). Two new sourced parameters,
`k_prime_d_a1`/`k_prime_d_a2` (Coleman Table A2's linear/quadratic coefficients, tier PLAUSIBLE for
wine, transferred/SPECULATIVE for beer — same pattern as `k_prime_d` itself in each file), plus
`k_prime_d_t_floor` (11 C, Coleman's own studied-range floor): the quadratic's vertex sits at
~11.3 C, below which it unphysically predicts *more* death as it gets *colder* — an extrapolation
artifact outside Coleman's fitted range, not a real effect, so temperature is clamped to the floor
before the quadratic is evaluated (no ceiling clamp — the upward acceleration above 11.3 C is the
sourced, physically-correct "heat causes stuck fermentations" direction and Coleman's own fit runs
to 35 C). Wired into both wine and beer's shared `_PRIMARY_FERMENTATION_MODIFIERS` (`EthanolInactivation`
is a shared Process, D-13); at `T = T_ref = 20 C` the factor is exactly 1, so §2.2, the Coleman
reconstruction, and every other 20 C-anchored test are untouched by construction.

**Measured before/after against Varela (the D-56 comparison this was meant to improve):**

| condition | pre-D-57 | post-D-57 | Varela (real) |
|---|---|---|---|
| N=300 hours-to-dryness | 83.0 h | 89.0 h | 170 h |
| N=300 gap ratio | 2.05x | 1.91x | — |
| N=50 hours-to-dryness | 176.0 h | 314.0 h | 700 h |
| N=50 gap ratio | 3.98x | 2.23x | — |
| N50/N300 duration ratio | 2.12x | 3.53x | 4.12x |

The N=300 in-range comparison barely moves (short run, death immaterial by dryness — exactly the
"inert on short/high-N runs" prediction), confirming that residual ~1.9x gap is a genuine
Coleman-vs-Varela cross-study difference the engine faithfully reproduces (Coleman's own reference
model, run at 28 C with Varela's inputs, gives 84.5 h — matching the engine, not Varela). The N=50
gap narrows from ~4x to ~2.2x, and the central D-56 structural finding (model under-predicts how
much severe N-deficiency slows fermentation, relative to an in-range baseline) survives but shrinks:
the model's N50/N300 ratio was 1.94x too small relative to Varela's 4.12x pre-fix; it is now only
1.17x too small. **This residual is left as an open, honestly small gap — a Bisson-sourced
nitrogen-gated transporter-capacity mechanism (D-56's original proposal) is no longer clearly
warranted at this size, chasing a ~1.17x residual against a single out-of-range data point risks the
same overfitting the D-56 calibration/validation firewall was built to prevent. Owner's call whether
to pursue it further or accept this as a documented model limit** (see Deferred, below — updated
from D-56's framing).

**A second, independent correction surfaced while finishing this comparison properly (advisor-caught
before commit): the benchmark's biomass assertion was reading the wrong state variable, unrelated to
the `k_prime_d` fix itself.** `_run_varela_condition` compared Varela's biomass to viable `X` alone.
Checked directly against the paper (WebFetch of the primary source, not assumed): Varela measures
TOTAL dry cell weight by gravimetric filtration ("dried...to a constant weight at 85 C") — dead and
viable cells combined, not a viable count. Because `EthanolInactivation` only *transfers* mass
between `X`/`X_dead` (D-13), `X + X_dead` is exactly conserved once nitrogen-limited growth stops
(~40 h in, confirmed flat to 5 significant figures for the rest of both runs) — so it is both the
methodologically-correct comparison and a strictly more robust one than a viable-only reading, which
depends on exactly when death has progressed to at the dryness-crossing instant. Corrected: total
biomass comes out ~3.38 g/L at N=300 (42% below Varela's 5.8) and ~1.40 g/L at N=50 (7% below
Varela's 1.5) — reproducing D-56 finding 3's already-documented Y_X/N cross-study numbers almost
exactly, which the old viable-only reading had never actually been measuring. The biomass assertions
now cleanly guard that growth-yield finding, separate from the duration assertions' death/uptake
timing — a cleaner split than before, and independent of whether the `k_prime_d` fix above landed at
all (total biomass is unchanged by it, being mass-neutral under the death transfer).

**A related tension worth flagging, not fixing:** at N=50 the model's own viable/dead split implies
~94–98% of biomass is "dead" by the time dryness arrives, while Varela separately reports **>97%
viability throughout** (LIVE/DEAD membrane-integrity fluorescence staining). Read carefully before
treating this as a new crisis: Coleman's own reference model shows the *identical* near-total `X_A`
crash at N=50 (0.099 g/L of a ~1.4 g/L total, matching the engine) — so this is a Coleman-vs-Varela
divergence the engine faithfully reproduces, not a new model defect. More importantly, `X_dead` is
documented (`inactivation.py`) as loss of *catalytic* (fermentative) capacity — the classical yeast
**vitality** concept — which is a different quantity from LIVE/DEAD's **viability** (membrane
integrity); the two are not expected to agree, and `k_prime_d` was fit to Coleman's sugar curves
(D-13/D-14), never to a viability count. Changing `k_prime_d`'s magnitude to chase agreement with a
viability assay it was never fit against would break the Coleman line-for-line reconstruction and is
out of scope here — flagged for whoever next touches death-rate calibration or wants a user-facing
"% viable yeast" output, not actioned by D-57.

**Test consequences (measured, re-banded, not loosened blindly — the D-46/D-51/D-53 discipline):**
`test_validation_varela2004.py`'s three tests re-banded to the new measured values/ratios above (and
its docstring rewritten to state the corrected diagnosis, not the stale one). `test_media.py`'s
`EXPECTED_MODIFIERS` gained `coleman_death_temperature` for both media. One genuine downstream
consequence, not a bug: `test_vicinal_diketones.py::test_warmer_ferment_is_cleaner_the_diacetyl_rest`
asserted wine's 28 C/45-day run clears diacetyl below the ~0.1 mg/L lager-perceptibility threshold;
post-fix it measures 0.162 mg/L (was ~0.03) because the now-correctly-faster warm-ferment death
leaves less viable/reductase-capable biomass surviving to day 45 than the old, under-scaled death
rate did. The monotonic "warmer is cleaner" direction and a large (~3x) magnitude both still hold and
are what the re-banded assertion now checks; the sub-perceptibility claim was retired as no longer
true for isolated-yeast reductase at this exact duration (real wine also gets MLF bacterial diacetyl
reduction, unmodelled here). 664 passed (unchanged count — a fix, not new tests), ruff+mypy clean.

**Method beat worth remembering: two advisor() passes, each correcting a premise the transcript would
otherwise have carried forward uncritically.** Pass 1 caught that mechanism 1 (the task's whole
starting premise) was stale documentation, not a live bug — verified with a probe run showing `X`
already declines substantially via inactivation before the advisor call, then confirmed structurally
against Coleman's own eqs. via `test_coleman_reconstruction.py`. Pass 2, after the `k_prime_d`
discovery, directed the single decisive check (Coleman's own reference model at 28 C with Varela's
inputs, both with and without the temperature-scaling bug) that turned "this looks like a T-scaling
bug" into a demonstrated, quantified one — and flagged the blast-radius grep and the honest-residual
framing before declaring the fix complete. Both times the initial "two-mechanism build" framing was
half-wrong; the data reshaped it into "fix one sourced bug, measure, then let the owner decide if a
much smaller residual is worth a new mechanism" — the same D-48/D-49/D-51 pattern this project keeps
hitting when a delegated diagnosis is checked against current code rather than trusted at face value.

## D-58 — MLF v2 sub-items research: `BrettSenescence` twin re-confirmed declined; ethanol-toxicity death built

**Status: IMPLEMENTED (2026-07-08).** Picked up the two remaining D-52 "MLF v2 further refinements" sub-items
(`BrettSenescence` twin; a separate `molecular_so2_death_scale` for `MalolacticDeath`) as the next
task. Before building either, two independent literature-research agents (opposite angles — one
hunting for evidence a Brett senescence mechanism is needed, one hunting for evidence the model's
existing "persists indefinitely without SO₂" framing is correct) were run in parallel, mirroring the
D-53 method that overturned the analogous MLF senescence premise.

**Finding 1 — `BrettSenescence` twin: D-52's decision holds, converged from both angles.** Neither
agent found evidence for a generic, free-running, age-based decline mechanism. Every measured decline
in unsulfited Brett traces to a specific stressor: molecular SO₂ (Serpaggi et al. 2012 — VBNC loss of
culturability is SO₂-induced and *reversible* on pH shift, not an aging phenomenon), substrate
exhaustion (Vigentini et al. 2008 — decline only on fructose depletion, not with fructose present),
or ethanol toxicity (Barata et al. 2008, below). No source describes decline attributable to elapsed
time alone. **Do not build a generic `BrettSenescence` twin** — same conclusion as D-52, now
literature-checked rather than reasoned from folk wisdom alone.

**Finding 2 — the "persists indefinitely" wording is an overstatement; soften it.** Barata et al.
2008 (*Int. J. Food Microbiol.* 121(2):201–207, doi:10.1016/j.ijfoodmicro.2007.11.020 — full text
verified) directly contradicts *literal* indefinite persistence: in closed-system model wine (12%
v/v ethanol, pH 3.50, no residual sugar, 25 °C, no SO₂) Brett populations bloomed to ~10^8 CFU/mL
then declined to complete loss of culturability by ~1200 h (~50 days) — growth at 8% v/v ethanol,
death at 14%, upper growth ceiling ~14.5–15% (their Table 2). Plate-count "death" is complicated by a
VBNC state (counts ran >10× below methylene-blue-active cells, and resuscitation is strain-dependent)
and Cibrario et al. 2019's "decades" persistence (doi:10.1371/journal.pone.0222749) is cellar/genotype
re-isolation across vintages, not continuous single-population viability. **Two-layer verdict:**
literal "one population persists forever" is unsupported; operational/reservoir tenacity (VBNC +
resurrection + cellar/biofilm reservoir across vintages) is well-supported. D-40/D-52's "persists
indefinitely — the honest reflection of tenacity" should be read/quoted going forward as "no positive
evidence for spontaneous decline without SO₂; SO₂, ethanol toxicity, and substrate exhaustion account
for observed die-off" — tighter to the evidence, same practical conclusion at the model's ≤30-day,
cellar-temperature run horizon (D-53's "empirically inert at these timescales" logic applies here
too).

**Finding 3 — a genuine, sourced, currently-missing mechanism surfaced as a side effect: Brett has no
ethanol-toxicity upper gate.** Checked directly against `brett.py` (not asserted): `BrettGrowth`'s
only ceiling is the intrinsic logistic carrying-capacity brake `(1 − X_brett/K)`, which drives growth
to zero as `X_brett → K` — a **plateau**, never a decline. `BrettDeath` (D-40 pt3) is SO₂-driven only
(`total_so2 ≤ 0` returns identically zero). So today's model, run dry/unsulfited/high-ethanol, would
plateau at the carrying capacity — it structurally **cannot** reproduce Barata's bloom-then-death
dynamic, because nothing in the model currently gates on ethanol toxicity. This is a real gap, not
something the existing brake already covers (an advisor-flagged contradiction in the first-draft
research report, verified against the code before writing this record). `BrettGrowth` already treats
ethanol purely as a carbon source (Monod, D-40 pt2) with no upper wall — unlike MLF's
`ethanol_tolerance_mlf` Luong-wall gate, deliberately omitted from Brett per D-40's design warning
against copying the MLF gate (Brett is markedly more ethanol-tolerant, so an MLF-style wall would be
wrong — but Barata shows tolerance is bounded, not unlimited).

**Built (owner chose to build, not defer).** `BrettGrowth` already uses ethanol as a carbon
*source* (low-concentration regime); Barata's toxicity is a *high*-concentration effect on the
*same* state variable, so the death term is reconciled with the growth term rather than layered on
top. Implementation (`fermentation.core.kinetics.brett`):

- **`brett_ethanol_survival_factor(E, params)`** — a shared helper, ∈ [0, 1]. Deliberately NOT the
  standard whole-range Luong wall (`(1 − E/E_max)^n`, decaying continuously from `E = 0`) that MLF
  uses: a Luong wall centered near Barata's ~118 g/L ceiling would already suppress Brett
  substantially at ordinary wine strength (~90–105 g/L) — the exact mistake the Brett gate's
  no-ethanol-term design already avoids (D-40). Instead a **threshold** form: exactly 1 (no effect)
  for `E ≤ brett_ethanol_toxicity_onset`, easing smoothly (C1, `n = 2`, no BDF kink) to 0 by
  `brett_ethanol_toxicity_ceiling`. Sourced boundaries: onset 110 g/L (~14% v/v, Barata's death
  onset), ceiling 118 g/L (~15% v/v, Barata's growth ceiling) — both via the codebase's standard
  ethanol-density conversion, used as fixed values, not fit to Barata's curve.
- **`BrettGrowth`** multiplies this factor into its rate as an upper wall, alongside the existing
  low-concentration ethanol Monod — the combined shape is a *hump* (source at low E, flat across
  normal wine strength, arrested near the ceiling), the reconciliation the design fork called for.
  Verified this leaves ordinary wine strength byte-for-byte unaffected: a probe run of the standard
  22-Brix test scenario tops out at E ≈ 106.6 g/L, safely below the 110 g/L onset, so the existing
  `test_pitch_brett_post_af_at_high_ethanol` integration test (which explicitly asserts "no ethanol
  wall arrests Brett at full-strength wine ethanol") needed no change.
- **`BrettEthanolToxicity`** — a new sibling `Process` to `BrettDeath`, NOT an added term inside it
  (keeps `BrettDeath`'s existing "exactly 0 without SO₂" docstring/tests byte-for-byte true).
  `r_death = k_death_brett · X_brett · (1 − survival(E)) · arrhenius(T, E_a_death_brett, T_ref)` —
  the `BrettDeath` `1 − g_SO₂` idiom, reusing `k_death_brett`/`E_a_death_brett`/`T_ref` rather than
  sourcing new magnitude/temperature params (Barata measured at one fixed 25 °C, so no independent
  activation energy exists; reuse mirrors `BrettDeath`'s own documented arrest-scale = kill-scale
  simplification). Needs no SO₂ — the entire point. Pitch-gated alongside `BrettDeath` in
  `_BRETT_GATED_PROCESSES`/`_BRETT_PROCESSES`. Exact zero guard at/below onset (no pH solve ever).
- **Scope limitation, documented not silently dropped:** Barata's most-cited number (a 12% v/v,
  no-SO₂, 50-day crash) is explicitly confounded in the source — bloom-on-trace-carbon *then*
  starvation-plus-ethanol-stress — and 12% v/v (~95 g/L) sits below the onset, so this Process alone
  predicts no decline there. Only the distinct, unconfounded per-concentration boundary data (grow
  ~8%, death onset ~14%, ceiling ~14.5–15%) is modelled; a starvation-driven decline mechanism, if
  ever wanted, is separate and not scoped here.
- **New params (wine-only):** `brett_ethanol_toxicity_onset` (110 g/L), `brett_ethanol_toxicity_ceiling`
  (118 g/L), `brett_ethanol_toxicity_exponent` (2.0) — all speculative, sourced from Barata et al. 2008.
  No new state slots (reuses `X_brett`/`X_brett_dead`/`E`).
- **Tests (12 new, `tests/test_brett.py`):** exact-zero guard at/below onset, neutral transfer,
  touches, monotonicity onset→ceiling, reused-Arrhenius warm-accelerates direction, no-catalyst
  guard, speculative tier, the survival-factor helper's own boundary values (direct unit tests), the
  growth wall leaving normal wine strength unaffected, the growth wall arresting growth at the
  ceiling, and a headline integration test: an unsulfited 26-Brix (~13% ABV, above onset) scenario
  crashes a growing Brett population (`X_brett_dead` fills, `X_brett` declines from its peak) while
  a 22-Brix (~11% ABV, below onset) control keeps growing — plus a carbon/nitrogen conservation test
  with the mechanism active. The separate `molecular_so2_death_scale` split (the other original
  D-52 sub-item) remains available but still zero-fidelity-gain per D-52's own reasoning —
  deprioritized, not built here.

**Full-suite result:** 676 passed (664 + 12 new), ruff + mypy clean, no existing test's assertion
needed to change (the onset threshold sits above every existing scenario's finished-wine ethanol).

**Post-implementation advisor pass caught two more honest gaps, both fixed in docstrings (not
re-engineered):**

1. **A real C1 discontinuity at the onset, mis-described as fully smooth.** The shifted-threshold
   survival factor is C1 at the ceiling (derivative → 0 from both sides, verified numerically) but
   NOT at the onset — the flat pre-onset region (derivative 0) meets the power-law ramp (derivative
   `−n/span`, e.g. −0.25 at the shipped n=2/span=8) with a finite jump. The original docstring
   claimed blanket "C1, no BDF kink" — true only at the ceiling. Corrected to state both facts
   precisely rather than overclaim. Verified benign (not the D-40 pt2 C0-step pathology that
   actually caused a solver blow-up): a bounded Jacobian entry, and the full suite — including the
   headline test, which integrates straight through `E = onset` — passes without incident. Not
   re-engineered into a two-breakpoint smootherstep (would fix a cosmetic claim, not a real problem)
   — the docstring now says so explicitly instead.
2. **A magnitude tension, not re-tuned.** Reusing `k_death_brett` (0.03/h) gives a ~23 h full-kill
   half-life at the ceiling — much faster than real 14–15% ABV reds being famously Brett-prone
   (self-clearing in a day would contradict that folk wisdom) and than the multi-week timescale of
   Barata's own 12%-condition decline (though that number is the starvation-confounded result
   already excluded from this Process's scope, not a clean ethanol-only rate). Barata's Table 2
   reports boundary *concentrations* (grow/death-onset/ceiling), not a decline *rate* at any one
   level, so there is no sourced number to replace the reuse with — flagged in the docstring as the
   value to revisit if a future source supplies one, not silently tuned down now.

**Method beat:** two parallel Opus research agents, deliberately opposite-angle (one arguing
"decline exists", one arguing "persistence is real") to avoid one-sided confirmation, then a
same-session advisor() pass that caught a self-contradiction in the second agent's report (claiming
both "the existing brake already covers this" and "this is genuinely missing physics" — those can't
both be true) before it reached the owner — resolved by reading the actual `BrettGrowth`/`BrettDeath`
code rather than trusting the agent's synthesis. A second advisor-flagged risk (a standard Luong
wall would suppress Brett at ordinary wine strength, contradicting its established ethanol-tolerant
niche) was verified empirically (the 22-Brix E≈106.6 g/L probe) before committing to the threshold
functional form over the more obvious wholesale-reuse of the MLF Luong wall. A THIRD post-build
advisor() pass (above) caught the C1-claim overclaim and the death-rate magnitude tension — both
fixed by honest documentation, not by silently absorbing or re-tuning.

## D-59 — Validation-direction research sweep: strain collision found, SO₂ overshoot closed to a documented limit, N-gap and beer paths scoped

**Status: RESEARCH ONLY (2026-07-08), no code changes.** M2 physics + refinements complete through
D-58; owner was asked to pick the next direction among validation / UX / new physics scope (all three
had been sitting as open candidates since D-55). Owner picked **validation**, then asked to research
all three of its own open sub-threads before committing to any build: (1) the D-56/57 residual
N-specific rate gap (~1.17×), (2) the D-51 residual SO₂/acetaldehyde overshoot (1.15–1.45×), (3)
broadening validation coverage beyond the single Varela 2004 check. Run as **6 parallel Opus research
agents, 2 per thread, deliberately opposite-angle** (a "pursue"/mechanism-feasibility agent and a
"skeptic"/cost-and-firewall agent per thread) — the same method D-53 and D-58 used to avoid one-sided
confirmation. All agents were literature-research only (WebSearch/WebFetch), explicitly barred from
touching code.

**Finding 0 — the highest-leverage discovery, cuts across all three threads: Coleman 2007 (the
model's fit source) and Varela 2004 (the model's only independent validation) are the SAME strain
lineage (Prise de Mousse / EC1118 derivatives).** Surfaced by the coverage-research agent while
checking Palma 2012's strain against the other two. This means the project's "independent"
validation has, until now, effectively been on one strain twice — a materially weaker validation
posture than it looked. It reframes threads 1 and 2 (both checked only against Varela) as validated
against a narrower base than assumed, and makes strain-independent coverage (Palma, below) the
single highest-value forward move — ahead of either mechanism build.

**Finding 1 — N-specific gap (D-56/57): do NOT build a transporter mechanism yet; the real gap is
qualitative, not the 1.17× ratio.** The "pursue" agent found a genuinely firewall-safe candidate
mechanism (Salmon 1989 — sugar-transport catabolite inactivation triggered by N-exhaustion; Palma
2012 — per-cell glucose-uptake Vmax falls to ~20% of initial under N limitation, a threshold/
switching response, not proportional Monod decay), independent of Varela 2004. But both agents
converge on not building it now, for two separate reasons: (a) **Varela's own paper's central thesis,
and the model's Cramer 2002 → Coleman 2007 lineage, attribute the N-rate effect to viable-cell
biomass, not per-cell rate** — so per-cell inactivation may be the wrong lever entirely; (b) **there
is no third independent, in-regime, numeric severe-deficiency dataset** to check a new mechanism
against without re-using Varela, which would be circular (Varela is the project's only independent
wine-kinetics check). The skeptic agent additionally found the target itself is soft: Varela's 4.12×
ratio combines **two different endpoint definitions** — N=300's "170 h" is time-to-true-dryness,
N=50's "700 h" is time-to-93%-consumption of a fermentation that never went dry, arresting at **16
g/L residual sugar**. **The model does not reproduce this arrest at all — it always finishes dry.**
That qualitative miss, not the 1.17× ratio, is the real finding; the ratio is downstream of it. Cross-
strain variability in N-sensitivity is documented at >2.5× (Gutiérrez et al. 2012, 23 strains),
comfortably swamping the residual. **Recommended next step (not yet run): a zero-cost internal
diagnostic** — compare the model's N50 viable-biomass trajectory against Varela's measured cell
counts. If biomass matches and rate is still fast, a per-cell term is justified (source from Salmon
1989, NOT Palma — see Finding 3). If biomass runs too high, the fix is recalibrating existing death/
yield terms, not a new mechanism.

**Finding 2 — SO₂/acetaldehyde overshoot (D-51): real two-agent consensus, close it out as a
documented limit.** Two independent, mutually reinforcing arguments:
- **Affinity arithmetic is decisive on its own.** Acetaldehyde's SO₂-binding affinity is 100–370×
  tighter than its pyruvate/α-KG competitors (Burroughs & Sparks 1973 Kd values); at these bisulfite
  concentrations acetaldehyde is already ~99% bound. Freeing enough at the 200 mg/L dose to match the
  field slope would require pulling free bisulfite down ~200×, which would need a competing pool on
  the order of tens of mol/L — physically impossible. **No Langmuir-type binder pool, of any affinity
  or size, can close the high-dose end** — independent of whether the field anchor itself is right.
- **The field anchor is a category mismatch.** The "0.39 mg/mg, linear across 50–200 mg/L" reference
  traces to a **cross-sectional regression across ~12 heterogeneous commercial wines, measured only to
  ~124 mg/L** (Marrufo-Curtido, Ferreira & Escudero 2022, *Foods* 11(3):476 — the equation
  `−4.4 + 0.39·W_tSO2` appears there verbatim, quoted from prior survey work; the repo currently
  attributes it to Jackowetz & Mira de Orduña 2013, which should be reconciled). It is a population
  survey with a documented pH confound, not a within-wine SO₂ titration, and "linear to 200 mg/L" was
  never tested for curvature — it's an extrapolated straight-line fit. Separately, the broader
  controlled-dose-response literature (Cornell Research Focus 2011-3; Jackowetz et al. 2011; an OENO
  One industrial-strain study) reports a 0.2–0.5 mg/mg range (2.5× strain-driven spread) with study-to-
  study averages already disagreeing by ~1.2× — comparable to or larger than the model's 1.15–1.45×
  deviation. **The gap is at or below the reference data's own discriminating power.**
- **Decision: accept the D-51 residual (1.15–1.45×) as a documented, structurally-explained model
  limit. No fourth binder pool.** Two cheap loose ends worth doing, not building: reconcile the D-51
  citation (Jackowetz & Mira de Orduña 2013 → Marrufo-Curtido et al. 2022 for the exact equation), and
  optionally print acetaldehyde's bound fraction at the 200 mg/L dose (should read ~0.99) to confirm
  the "production, not binding, is the real remaining lever" claim — the D-48 `k_acet_so2_induced`
  coefficient, not the binding equilibrium, is where the residual slope actually lives if it's ever
  revisited.

**Finding 3 — broadening coverage: Palma 2012 is a genuine strain-independent validation candidate,
worth digitizing; beer-side independent validation is currently blocked by data access, not by
absence of a target.** Palma, Madeira, Mendes-Ferreira & Sá-Correia 2012, *Microbial Cell Factories*
11:99 (doi:10.1186/1475-2859-11-99) uses strain **PYCC 4072** — different from Coleman/Varela's Prise
de Mousse lineage (Finding 0) — at 320 vs 90 mg N/L, with **n=3 replicates and SD error bars on
Figure 1**, so digitization noise (~2–5% of axis range) stays well under the model's 1.17–2× gaps.
Glucose and ethanol curves (linear axes, ~8 points each) are directly usable; biomass is CFU/mL (log
scale) and not worth converting. A complementary no-digitization dataset (MDPI 2024, *Fermentation*
10(8):386, real density-time tables) was also found but has nitrogen as a *fitted* parameter, not
measured — useful only as a trajectory-shape sanity check, not a nitrogen-mechanism validator. **Beer
side:** no genuinely independent, in-regime (isothermal ale, ~1.048 OG), numeric time-series dataset
is publicly accessible — the two richest candidates (Zamudio Lara et al. 2022, de Andrés-Toro et al.
1998) are confirmed via `beer_generic.yaml` provenance and D-15/D-19 to be the model's own beer fit
sources (circular if reused). The right-regime data that exists (Reid et al. 2021's two ale datasets,
bracketing the benchmark almost exactly) is proprietary and never published numerically. The only
usable independent option found is a **lager** dataset (Speers et al. 2003, reconstructable via Reid
et al. 2021 Table 2's fitted logistic parameters) — off-regime, but usable as an explicit cross-regime
Arrhenius-scaling stress test rather than a same-regime validation. Beer independent validation is
recommended **deferred** pending an accessible in-regime dataset, mirroring how wine validation waited
for Varela.

**Advisor catch — Palma 2012 was assigned two mutually exclusive roles by two different agents, and
neither could see the conflict from inside its own scope.** The N-gap "pursue" agent (Finding 1)
proposed sourcing a per-cell transporter-inactivation mechanism's parameters from Palma's Vmax
measurements; the coverage agent (Finding 3) independently proposed the same paper's fermentation
curves as a validation dataset. Under the project's validation/calibration firewall these are
incompatible — sourcing a mechanism from Palma and then validating the resulting model against
Palma's own data would be self-confirming. **Resolved: reserve Palma for validation** (it's additive
— raises future validation power rather than spending it — and it's the only strain-independent
option available); if a per-cell N-transporter mechanism is ever built, it must source from Salmon
1989 instead (weaker parameterization, but firewall-clean). Recorded here so this constraint survives
across sessions: **Palma 2012's fermentation curves are earmarked for validation use; do not also
mine it for kinetic mechanism parameters.**

**Net outcome: no code changed this session.** Recommended next steps, in order of cost, none yet
started: (1) the N50 viable-biomass-vs-Varela diagnostic (internal, zero new code beyond a comparison
script); (2) Palma 2012 digitization as a second independent wine validation point (glucose + ethanol
only); (3) the D-51 citation reconciliation + optional bound-fraction print; (4) a decision on whether
to build the Speers/Reid lager cross-regime beer check or defer beer validation entirely. Owner has
not yet chosen which to start; UX and new-physics-scope remain the other two untouched top-level
directions from D-55's "next milestone" fork.

## D-60 — Palma 2012 digitization: second independent-data benchmark built, strain-independent N-gap corroborated, absolute timing gap flips direction (confounded, not a fidelity signal)

Owner picked up D-59's "Palma 2012 digitization" recommendation directly. Built
`tests/benchmarks/test_validation_palma2012.py`, the project's second independent-data
validation file (after Varela 2004, D-56/D-57) and its first against a genuinely
different strain (PYCC 4072, not Coleman/Varela's Prise de Mousse — D-59 Finding 0).

**Digitization:** Figure 1 (panels C glucose, D ethanol) was fetched as its original
CC-BY image via the PMC Open Access S3 mirror (`PMC3503800.1/1475-2859-11-99-1.jpg` —
the legacy FTP `oa_package` tarball route is now deprecated/404; PMC's own web viewer is
behind a proof-of-work JS gate that blocks `curl`; the S3 bucket, discovered via the
`oa.fcgi` API's stale-but-still-resolvable-through-the-new-layout link, was the working
path) and read off against a pixel grid calibrated to the panels' own axis ticks, at the
paper's confirmed real sampling times (0,6,24,48,72,80,96,144 h). Two of the paper's
three conditions digitized: CF (320 mg N/L) and LF (90 mg N/L); RF (LF refed with DAP at
72 h) deliberately deferred — a discrete mid-run intervention is a different validation
target (`add_dap` timing fidelity), out of scope for a first glucose+ethanol pass.

**Headline finding — the CF/LF absolute-timing gap not only persists on a second
dataset, it *flips direction*; the timing gap and the yield gap have DIFFERENT,
deliberately NOT conflated explanations, and neither is a fidelity signal:** at 20°C
(Palma's fermentation temp — exactly the engine's/Coleman's `T_ref`, so zero Arrhenius
extrapolation uncertainty, cleaner than Varela's 28°C), the engine reaches CF dryness at
~138 h against Palma's real ~72 h — ~1.9x *slower*, the opposite direction from Varela's
~1.9x too *fast* at 300 mg N/L/28°C. **Cross-checked against Coleman's own reference
model** (the same eqs-1-8 reconstruction `test_coleman_reconstruction.py` uses, re-run at
Palma's exact inputs: S0=200 g/L, N0=320 mg/L, pitch=0.018 g/L, 20°C): it dries at ~140 h,
~1.5% from the engine's ~138 h — **the engine faithfully reproduces Coleman at Palma's
own inputs**, so the gap to Palma is a genuine Coleman-vs-Palma difference, not an engine
defect (the exact D-57 argument, transplanted to a new dataset). **The timing gap's
best-supported explanation is strain, not protocol:** at 200 g/L glucose, S. cerevisiae is
strongly Crabtree-repressed and ferments even under full aeration, so respiratory carbon
diversion cannot explain a ~2x rate difference at this sugar level — PYCC 4072 (Palma)
and Prise de Mousse (Coleman/Varela) are simply different strains with different
fermentation rates, and this dataset's whole value is being the first strain-independent
check, so a gap here is expected, not a red flag. **Separately, the yield gap has its own,
narrower explanation:** Palma's real ethanol yield is only ~0.39-0.40 g/g glucose
consumed at both N levels (computed from the digitized endpoints: CF 78.9 g/L ethanol /
~199 g/L consumed; LF ~45.0 g/L / ~120 g/L) — well below the ~0.46-0.51 g/g anaerobic
range the engine itself uses (~0.48) — consistent with ethanol evaporating from a shaken,
cotton-stoppered 500 mL Erlenmeyer flask (120 rpm) over a multi-day shake; evaporation
affects the reported ethanol *level*, not the glucose-consumption *rate*, so it explains
the yield gap only, not the timing gap. A third, weaker data point — Varela's real CF
(28°C, warmer) took *longer* (170 h) than Palma's real CF (20°C, cooler, 72 h) — shows the
two "independent" datasets disagree with each other by ~2.4x, at least as much as either
disagrees with the engine (the same "gap is at or below the reference data's own
discriminating power" shape D-59 reached for the SO₂ overshoot); this is NOT read as a
clean temperature (anti-Arrhenius) comparison, since Varela and Palma are also different
strains — strain is confounded with temperature here, so no temperature-specific claim is
made. **Absolute CF/LF duration and ethanol level are therefore characterized as
regression guards (observed value + margin), never asserted as agreement targets against
Palma's raw numbers.**

**The regime-robust finding — corroborated on an independent strain:** comparing each
condition's own glucose-consumed *fraction* at 144 h (a ratio that cancels the yield/
evaporation confound, since both conditions share one flask protocol and only sugar
consumed — not ethanol produced — is compared) shows the engine still under-predicts how
much severe nitrogen limitation suppresses fermentation progress. Real Palma: CF ~99.5%
consumed, LF only ~60% (residual ~80 g/L, still visibly decelerating 122→80 g/L between
96-144 h — far from dry, deliberately NOT called "arrested", per D-58's overclaim lesson)
— ratio ~1.66. Engine: CF ~99.7%, LF ~79% (residual ~41 g/L) — ratio only ~1.26. Same
direction and shape as D-56/D-57's Varela finding and D-59's "model never reproduces
arrest" framing, now independent of strain — this is the load-bearing signal in this
dataset, not the absolute timing. Test
`test_palma2012_lf_vs_cf_progress_ratio_understates_palma` asserts the engine's ratio
stays below Palma's real ~1.66; the other two tests characterize CF's absolute duration
(band [125,150] h, gap-ratio band [1.7,2.15]x) and LF's absolute residual at 144 h (band
[35,48] g/L) as regression guards, matching the Varela file's established idiom
(observed + margin, not a loose pass; do not force-fit).

**Method note — two advisor() catches, not one:** the first (during the build, before
the test file was finalized) flagged that the draft summary was heading toward "model
runs too slow" as the headline, which reads as a fidelity gap rather than the
better-supported protocol/strain reading. The second (after commit, reviewing the
finished docstring) caught that the fix for the first catch had gone too far the other
way: it bundled the timing gap and the yield gap under one "shaken flask → respiration
and/or evaporation" story, but respiration cannot explain the timing gap at 200 g/L
glucose (Crabtree repression) — the two gaps needed separate explanations (strain for
timing, evaporation for yield only), which this entry and the test docstring now reflect.
This is the same overclaim discipline as D-58's "arrested" → "far from dry" softening,
applied a second time, to explanations rather than to the headline finding itself — a
reminder that the discipline has to be re-applied at each layer of the writeup, not just
the first one checked. The Coleman-reconstruction cross-check (run as a one-off probe
script, not added as a fourth permanent test — `test_coleman_reconstruction.py` already
carries the general Coleman-fidelity claim; re-deriving it a third time in-repo would be
redundant) was the decisive piece of evidence for the strain explanation, exactly
mirroring how the Varela file cites the same
reconstruction for its own 300 mg N/L point. 3 new tests, 679 passed (676+3), ruff+mypy
clean. No source code changed — this is a benchmark-only addition, no physics touched.

**Open / still not done:** RF (refeed) digitization, deferred as noted above; beer
independent validation (Speers/Reid lager cross-check, D-59) still undecided; the N50
viable-biomass-vs-Varela diagnostic (D-59's other recommended next step) not yet run.

## D-61 — N50 biomass diagnostic run (D-59 Finding 1 gate): biomass is not the culprit, so DON'T build the Salmon per-cell mechanism; + two D-59 loose ends closed

**Status: DIAGNOSTIC + DOC/PROVENANCE ONLY (2026-07-08), no physics changed.** Picked up
D-59's top-ranked cheap next step: the internal diagnostic D-59 Finding 1 set as the **gate**
before any stuck-fermentation / per-cell-rate mechanism build. Ran as a throwaway script
(`M:\claud_projects\temp\d61_varela_n50_biomass_diagnostic.py`, kept out of git per repo
etiquette — a one-shot diagnostic, not a permanent test; `test_validation_varela2004.py`
already guards the numbers it leans on). Advisor consulted before interpreting.

**The gate D-59 set:** compare the model's N50 *active/viable* biomass trajectory against
Varela's cells. If active biomass ≈ Varela's viable **and** the model still finishes ~2.2× fast
→ a per-cell N-gated rate term is justified (source Salmon 1989, firewall-clean; NOT Palma —
D-59 reserved Palma for validation). If active biomass runs **too high** → the fix is
recalibrating existing death/yield, not a new mechanism.

**Advisor's load-bearing correction to the gate (applied):** the comparison variable is the
whole diagnostic, not a footnote. Three *different physical quantities* are in play — model `X`
is **catalytic vitality** (fermentative capacity), Varela's >97% is **membrane viability**,
Varela's 1.5 g/L is **gravimetric total DCW**. Naively reading "model viable runs low → build
the per-cell term" off a vitality-vs-viability comparison would *manufacture* a false verdict
the Varela docstring itself warns against. So the diagnostic reports the **mean active `X`
across the 10%→90% sugar-consumption window** (the biomass actually fermenting), anchored on the
clean gravimetric-to-gravimetric **total** comparison, with the active figure as caveated
support only.

**Numbers (28 °C, research pitch 0.018 g/L, 240 g/L sugar — same setup as the Varela benchmark):**

| quantity | N=300 | N=50 | Varela |
|---|---|---|---|
| hours to dryness (model) | 89 h | 314 h (min S 0.85) | 170 h / 700 h* |
| model-vs-Varela speed gap | 1.91× fast | 2.23× fast | — |
| total biomass X+X_dead (model) | 3.38 | **1.40** | 5.8 / **1.5** DCW |
| mean active X in 10–90% window | 2.99 | **0.96** | ~5.63 / **~1.46** (DCW×97%) |
| active X across the window | 2.52→2.45 | **1.38→0.46** | (>97% viable throughout) |
| active fraction at 90% consumed | 72% | **33%** | >97% |

*Varela N50 never reaches dryness — it arrests at 16 g/L residual; "700 h" is time-to-93%-
consumption. The model always finishes dry (the qualitative miss D-59 Finding 1 flagged as the
*real* gap, downstream of the capacity-loss timing).

**Verdict — the gate *as D-59 posed it* is unevaluable, not un-fired; reframed onto the one clean
quantity it still hardens D-59's tentative "don't build yet" into "don't build":**
0. **D-59's gate rests on a category error.** It keyed both branches on *model-viable* vs
   *Varela-viable* cells (viable ≈ Varela → build per-cell; viable too high → recalibrate). But
   model `X` is **catalytic vitality** and Varela's >97% is **membrane viability** — the exact
   incomparable pair the advisor flagged. So the gate as literally posed **cannot be evaluated**;
   reading a verdict off model-`X`-vs-Varela-viability would manufacture an artifact. This is the
   most useful thing the diagnostic surfaced, and it is why the entry does not report a "biomass
   matches → build" firing even though total biomass does match (see next).
1. **Reframe onto the one clean, comparable quantity — total gravimetric biomass — and it MATCHES**
   at N50 (1.40 vs 1.5 g/L, ~7% low; reproduces D-56 finding 3's yield near-match). So the residual
   is **not a biomass-quantity error** in either direction: the model does not have *too much*
   biomass (ruling out "recalibrate death/yield to kill more," on evidence, not merely on the
   frozen-`k_prime_d` Coleman-reconstruction cost of D-57), and total is not low either.
2. **That localizes the residual to per-cell rate / capacity-loss timing** — nominally the "build"
   locus. Supporting (cross-strain-consistent) evidence: the model has *less* total biomass than
   Varela at BOTH N levels yet finishes faster (N300: total 3.38 vs 5.8, 1.91× fast; N50: 1.40 vs
   1.5, 2.23× fast), so the fast-finish is a per-cell *rate* feature present even in-range. That
   rate is **Coleman's own fitted rate** — the engine reproduces Coleman line-for-line
   (`test_coleman_reconstruction.py`; ~84.5 h at the N300 inputs, matching). Coleman and Varela are
   the **same strain** (D-59 Finding 0), so the gap is a genuine **Coleman-vs-Varela cross-lab
   difference**; a per-cell term calibrated to close it would calibrate the model *away from its own
   fit source* toward a different lab — a validation/calibration-firewall breach.
3. **The load-bearing reason not to build: that build locus is already occupied by `X_dead`, and
   Varela cannot adjudicate its timing.** Varela's arrest at 16 g/L is loss of fermentative capacity
   in membrane-viable cells — exactly what `X_dead` (catalytic-vitality loss, `inactivation.py`)
   represents. A Salmon 1989 catabolite-inactivation term would be **largely redundant** with
   existing machinery. The real open question is whether `X_dead`'s *timing/magnitude* under
   N-limitation is off (the model keeps enough active `X` — 1.38→0.46 g/L across the window — to
   grind to dryness where Varela arrests). But **model-vs-Varela cannot settle that**, because the
   only Varela quantity bearing on it (>97% membrane viability) is not comparable to `X_dead`
   (vitality) — the same category error as point 0 — and there is no third independent in-regime
   dataset to break the tie without re-using Varela.

**Decision: accept the D-56/D-57 residual N-gap as a documented model limit; do NOT build the
Salmon per-cell mechanism.** Reached from internal evidence (points 0–3), landing on the same
"documented limit, not a free action" outcome the advisor flagged as the honest possibility.

**Honest limitation of the diagnostic (per advisor's verify-first point):** it compares the model
trajectory against Varela's **endpoints + the single >97% viability figure**, not a digitized
viable-cell *time series* — the Varela benchmark only carries endpoints. The total-biomass anchor
is robust anyway (flat from ~40 h onward, so the endpoint is representative of the plateau), so
the gate is answerable now. A future strengthening could digitize Varela 2004's Figure-1 biomass-
over-time curve to check *when* biomass diverges; not blocking, and it would not change verdict
points 1–3 (they rest on the plateau + the in-range N300 result).

**Loose end A — D-59's acetaldehyde bound-fraction check (confirmed in-model).** Built a second
throwaway script (`M:\claud_projects\temp\d59_acetaldehyde_bound_fraction.py`) driving
`acidbase.speciate_so2` / `free_acetaldehyde` at a realistic finished-wine state (pH 3.4,
pyruvate ~30 / α-KG ~20 mg/L competitors) across the benchmark SO₂ doses. At the 200 mg/L dose
with acetaldehyde tracking the field increments (25.7/56.1/119 mg/L), the model reports
acetaldehyde **98.5% / 99.0% / 99.4% bound** — the ~0.99 D-59 Finding 2's affinity-arithmetic
argument predicted, confirming in-model *why* no fourth binder pool can free enough to move the
high-dose slope (acetaldehyde is already essentially fully sequestered; the residual slope lives
in D-48's `k_acet_so2_induced` production term, not the binding equilibrium). Only an unphysical
acetaldehyde overload (200 mg/L, carbonyl moles > SO₂ moles) drops it to 69%.

**Loose end B — D-59's D-51 citation reconciliation (done).** Verified against both primary
sources (WebFetch of the Marrufo-Curtido PMC full text; search-confirmed the Jackowetz survey
scope): the exact slope equation `W_acetaldehyde = −4.4 + 0.39·W_tSO₂` (R = 0.837, p < 0.001) is
from **Marrufo-Curtido, Ferreira & Escudero 2022, *Foods* 11(3):476** — a **12-wine** forced-
oxidation cross-sectional survey over a **20–124 mg/L** total-SO₂ range (NOT a within-wine
titration). The repo's D-51 entry attributed it to **Jackowetz & Mira de Orduña 2013**, which is
the *separate* 237-wine "Survey of SO₂ binding carbonyls" (Food Control 32(2):687–692) reporting
average binder *concentrations* (acetaldehyde 25/40, pyruvate 14/25, α-KG 74/31 mg/L red/white) —
the correct anchor for the finished-wine keto-acid *ranges* (D-49/D-50) but not for this slope.
Fixed: (a) a bracketed `[CITATION CORRECTED in D-61…]` note appended in place at the D-51 entry
(preserving the append-only log); (b) the correct Marrufo-Curtido attribution added to the two
live provenance strings in `acetaldehyde.yaml` where the 0.39 equation appears without a cite.
The category-mismatch argument D-51 built on the anchor is **unaffected — if anything strengthened**
(the true anchor is a 12-wine survey to only 124 mg/L, weaker discriminating power than an
imagined 237-wine regression). Jackowetz cites elsewhere (keto-acid ranges, pyruvate-as-second-
binder in `acidbase.yaml`) are legitimate and left as-is. `acetaldehyde.yaml` still loads; 31/31
acetaldehyde tests green.

**Net: no physics/source-code changed.** One DECISIONS correction note, two YAML provenance-string
fixes, two throwaway diagnostic scripts (out of git). The N-gap and SO₂-overshoot threads are now
both closed as documented, structurally-explained model limits. Remaining open validation threads
(unchanged from D-60): Palma RF digitization, beer independent validation (Speers/Reid lager cross-
regime check or defer), and the optional Varela Figure-1 biomass-time-series digitization noted above.

## D-62 — Palma 2012 RF (refeed) condition built: the DAP-refeed rescue is reproduced, but the engine INVERTS Palma's within-study RF-vs-CF ordering (same N-under-suppression gap, now via a dynamic intervention)

**Status: BENCHMARK-ONLY (2026-07-08), no physics/source-code changed.** Picked up the RF
digitization D-60 deferred — the third Palma 2012 condition, and the `add_dap` timing-fidelity
target D-60 explicitly flagged RF for. Two new tests in `tests/benchmarks/test_validation_palma2012.py`
(681 passed = 679+2, ruff+mypy clean). The D-60 digitization workspace survived in
`M:\claud_projects\temp\palma2012` (fig1.jpg, fulltext.xml, methods_dump.txt, calibrated panels),
so no re-fetch was needed. Advisor consulted before writing — it reframed the finding decisively.

**RF protocol (Palma Methods, verified against the fulltext, not assumed):** after 72 h the sluggish
LF broth (90 mg N/L) was split and one half refed with **230 mg N/L as 1.1 g/L (NH₄)₂HPO₄** (RF); the
other half stayed as the LF control. The engine reproduces this faithfully-to-the-additive: `add_dap`
1.1 g/L DAP, which the model's exact-stoichiometry `dap_nitrogen_fraction` (0.2121, VALIDATED) turns
into +233 mg N/L — Palma's stated 230 is the identical dose rounded.

**A real probe bug caught before it became a false finding.** The first probe ran RF through a bare
`simulate(process_set, …)` and found RF byte-identical to LF — which *looked* like "the engine
ignores the refeed", a headline-grade fidelity gap. It was wrong: `simulate` silently drops the
compiled `events`; the refeed must go through `compiled.run()` (→ `simulate_scheduled`). Verified by
inspecting the N slot: through `.run()`, N_RF jumps 0.00→0.233 g/L at 72 h; through bare `simulate`
it never moved. `_run_palma_condition` was unified onto `.run()` (byte-for-byte identical for the
CF/LF no-event case per the compile.py contract, and confirmed empirically — CF/LF bands unmoved).

**Advisor's load-bearing reframe (applied):** the engine RF dries at ~108 h vs Palma's real ~117 h —
tempting to read as "RF timing agrees" (unlike CF's ~1.9× gap). **The advisor flagged this closeness
as a CROSS-study comparison carrying the exact same ~1.9× strain confound that makes the engine's CF
timing untrustworthy — coincidental, not a signal. Do NOT build the test on it.** The confound-robust
axis is the **within-study RF-vs-CF ordering** (the same axis as the existing CF:LF ratio test):
- **Palma real:** RF finishes AFTER CF (RF/CF ~ 117/72 ~ **1.6**) — LF genuinely stalls, so the refed
  culture starts far behind and only catches up well after the never-stalled CF is done.
- **Engine:** RF finishes at-or-BEFORE CF (~108 vs ~138 h, RF/CF ~ **0.78**) — it under-penalizes the
  LF stall (engine LF is only mildly behind CF at 72 h, not stalled), so the refed culture is barely
  behind and beats CF.

**The engine INVERTS Palma's within-study ordering** (engine RF<CF; Palma RF>CF, both with wide margin
→ robust to digitization slop). This is the **same D-56/D-57/D-59/D-60 nitrogen-sensitivity shortfall
— the engine under-predicts how much severe N-limitation suppresses fermentation — now surfaced
through a DYNAMIC refeed intervention rather than a static contrast.** That the same gap reappears on
a fourth, mechanistically-different probe strengthens the D-61 verdict that it is a genuine, coherent
model limit (not an artifact of any one comparison).

**Mechanism VERIFIED, not inferred (the D-60 lesson applied):** re-ran the *corrected* (`.run()`) RF
and inspected viable biomass X directly — after the 72 h refeed X_RF rises **2.1 → ~7.6 g/L, peaking
~89 h**, while X_CF already declines from its **~61 h** peak (~3.3 g/L). So the engine's fast RF finish
is driven by a large *late* biomass burst on the refed nitrogen (an observation, not a story bolted on
after). Note: RF's ~7.6 g/L peak exceeds CF's ~3.3 g/L on nearly-equal total N (323 vs 320 mg/L)
because a late N dump lands when sugar is still abundant (~152 g/L) — a state-dependent biomass-yield
consequence of the timing, worth flagging but not itself a defect.

**Test design (matches the D-60/Varela idiom — observed + margin, never force-fit against Palma):**
`test_palma2012_rf_refeed_rescues_the_sluggish_lf_to_dryness` asserts the confound-robust rescue (RF
reaches dryness; the LF control is still ~41 g/L at 144 h — RF vs LF differ ONLY by the dose, so
strain/evaporation cancel) plus an absolute RF-dryness regression guard ([95,120] h, explicitly NOT
vs Palma's 117 h). `test_palma2012_rf_vs_cf_ordering_is_inverted_relative_to_palma` asserts the
discriminating within-study inversion (engine RF<CF; Palma's digitized RF>CF pinned as literals).

**Open / still not done (unchanged from D-60/D-61):** beer independent validation (Speers/Reid lager
cross-regime check or defer); optional Varela Figure-1 biomass-time-series digitization. The Palma 2012
dataset is now fully exercised (all three conditions built).

## D-63 — Beer-side independent check: the accessible lager data is single-temperature (confounded), so build an honest cross-regime Arrhenius stress test — NOT a lager validation — and defer the confound-cancelling ratio test pending Speers 2003's controlled series

**Status: BENCHMARK-ONLY (2026-07-09), no physics/source-code changed.** Owner picked up the D-59
beer-validation fork, choosing the Speers/Reid lager option over continued defer. New file
`tests/benchmarks/test_beer_temperature_response.py` (3 tests; 684 passed = 681+3, ruff+mypy clean).
Advisor consulted twice — once on framing before any web research, once on the concrete test design;
both reframes applied. This is the first benchmark ever to exercise the beer Arrhenius `E_a`'s, which
have been **inert in every prior benchmark** (all isothermal at `T_ref` = 20 °C, so f = 1).

**The data investigation decided the whole design.** The advisor's linchpin: the test forks entirely
on whether the source spans more than one temperature.
- **Accessible source is single-temperature.** The only freely reconstructable lager curve is Reid,
  Josey, MacIntosh, Maskell & Speers 2021 (*Fermentation* 7(1):13, doi:10.3390/fermentation7010013),
  Table 2 — Australian lager, OE **14.1 °P, single starting temperature 10 °C**, 3-parameter ADF
  logistic B = 0.06372 h⁻¹, midpoint M = 51.22 h (≈ 2.1 d). Fetched via the Heriot-Watt open-access
  mirror (`pure.hw.ac.uk`, `pdftotext`); MDPI/ResearchGate/academia.edu all 403 the fetcher.
- **Multi-temperature signal is paywalled AND likely confounded.** The temperature effect (rate ↑
  with starting temp, p<0.01) lives in Speers, Rogers & Smith 2003 (*J. Inst. Brew.* 109(3):229–235,
  doi:10.1002/j.2050-0416.2003.tb00163.x), which is Wiley-paywalled (no accessible free full text),
  and its effect is a regression across many **industrial** batches (brand/wort/pitch co-vary with
  temperature) — so even if obtained it may not be a clean controlled series.

**Why a single-temperature lager band would be dishonest (advisor's call, D-59's defer sharpened).**
Comparing "engine ale-yeast Arrhenius extrapolated to 10 °C" vs "real lager yeast (*S. pastorianus*)
in a 14.1 °P industrial wort at 10 °C" conflates the Arrhenius law with the organism + wort +
pitch-rate difference. Empirically the engine's low-pitch (0.6 g/L homebrew-like) 10 °C run hits its
attenuation midpoint at ~6.2 d — **~2.9× slower** than Speers' ~2.1 d industrial midpoint, a gap
dominated by pitch + organism. Guarding that gap as a regression band would guard the confound. The
file therefore **deliberately does not** assert the engine reproduces the 51 h midpoint.

**What was built instead — three claims from the engine's OWN 20 °C vs 10 °C runs** (midpoints 2.79 d
and 6.21 d; the 2.23× slowdown matches the E_a-predicted 2.22–2.25× almost exactly). The advisor's
key correction: split the recovered number (apparent E_a ≈ **55.3 kJ/mol**, round-tripping the input
E_a_uptake 55.1 / E_a_growth 55.9 kJ/mol) into two labeled claims, because one band conflates two
purposes:
1. **Wiring / regression guard** — apparent E_a ≈ input (band [50, 60] kJ/mol). Guards, on beer and
   over a full ferment composing BOTH growth and uptake, that the Arrhenius modifiers stay wired into
   fermentation timing — the **D-57 frozen-modifier bug class**. Not a strict duplicate of the
   existing Arrhenius tests, which are directional-only, uptake-only, and on wine.
2. **Reality check (the honest headline)** — the SAME E_a sits inside the range commonly reported for
   *S. cerevisiae* alcoholic fermentation (order ~40–90 kJ/mol, well under ~100; band [35, 100] a
   deliberately generous envelope, exact per-study figures NOT read in-source — the primary kinetic
   papers 403'd, only Reid 2021's B/M were read from an opened PDF). This is the **only reality-
   touching claim**, and it has teeth: the engine's ~55 kJ/mol lands inside while the ~265 kJ/mol de
   Andrés-Toro lumped-fit artifact the beer file rejects is excluded by an order of magnitude — a
   verdict robust to the exact edges, while staying humble about the organism gap.
3. **Cross-regime order-of-magnitude anchor (CONFOUNDED, loose)** — 10 °C reaches 90 % attenuation in
   a "cold lager ~1–2 weeks" window [5, 25] d (engine ~12 d). Deliberately loose; only catches an
   order-of-magnitude-wrong temperature model. The low-pitch assumption is why the engine sits at the
   slow end and misses Speers' fast industrial timing.

**Naming (advisor):** NOT `test_validation_speers2003` — Speers is not load-bearing in any assertion
(the [5,25] d band is our own cold-lager judgment, and the engine intentionally misses 51 h). Named
`test_beer_temperature_response` for what it does.

**The confound-cancelling ratio test — the version with genuine signal — is DEFERRED, now on
EVIDENCE from the primary source, not access failure.** A rate *ratio* across two temperatures cancels
the lager-vs-ale absolute-kinetics difference, isolating the temperature axis; it needs a *controlled*
temperature series (fitted rate/midpoint at ≥2 temperatures on one wort+yeast). **Owner obtained the
paywalled Speers 2003 PDF mid-session and it was read in-source (2026-07-09) — it is NOT such a
series, exactly as the advisor predicted.** Three disqualifying facts from the paper: (1) "the
starting temperature was **brand dependent**, the data for each brand was pooled" — temperature
co-varies with brand, i.e. with wort/gravity/yeast; (2) fermentations "were started at various
specified temperatures and **allowed to free rise to set temperatures**" — non-isothermal, so no
single temperature per curve for an Arrhenius fit; (3) Table I tabulates only P/B/M/P₀ with **no
per-brand temperature values at all** — temperature enters only as a regression factor (p<0.01 on
rate, p<0.001 on midpoint), so not even raw (T, rate) pairs are recoverable. A ratio built from it
would still carry the full brand confound. The deferral therefore stands on the data's structure, not
on access. The reusable helper `_apparent_activation_energy` remains the drop-in point should a
genuinely controlled isothermal series (same wort+yeast, ≥2 temperatures) ever surface.

**Firewall (prime directive 2): clean.** Engine E_a's derive from the Coleman 2007 wine fit; the
reference data is Speers/Reid lager — disjoint sources, so the comparison is not self-confirming.

## Deferred (decide early in the relevant milestone)

- ~~**pH / acid model richness**~~ — **decided in D-18** (full charge-balance solver),
  built after the byproducts beat; **solver landed 2026-06-30** (`core.acidbase`,
  `fermentation.analysis`) — see D-18 "Resolution".
- ~~**Stochastic ensemble API**~~ — **decided in D-24 and IMPLEMENTED 2026-07-01**
  (`runtime/ensemble.py`): triangular Monte-Carlo over the `Uncertainty` bands, scoped to
  the active Process set's reads, nominal + median + P5/P95 band, per-member conservation.
- ~~**H₂S CO₂-stripping volatilization sink**~~ (D-29 follow-up) — **decided + IMPLEMENTED in
  D-42 (2026-07-06)**: `HydrogenSulfideVolatilization` sweeps the volatile `h2s` into a new
  carbon-free `h2s_gas` headspace pool on the CO₂-evolution flux, so `h2s` is now the µg/L
  *residual* reality shows and `h2s + h2s_gas` is cumulative produced. The exact ester D-20/D-21
  precedent but simpler (carbon-free ⇒ no ledger weighting). See D-42.
- ~~**Post-fermentation / autolytic H₂S source + copper fining**~~ (the two D-42 deferred items) —
  **decided + IMPLEMENTED in D-44 (2026-07-06)**: `AutolyticHydrogenSulfide` is a yield on the D-34
  autolysis flux (opt-in, wine-only) whose non-flux-linked form makes it accumulate un-stripped as
  *residual* post-dryness — the reductive fault; `add_copper` precipitates it as CuS (stoichiometric,
  ledger-neutral). See D-44.
- ~~**Mercaptan (thiol) pool + copper mercaptide**~~ (the D-44 open fork) — **decided + IMPLEMENTED
  in D-45 (2026-07-06)**: owner chose a carbon-bearing `mercaptans` pool (methanethiol), formed
  autolysis-linked with carbon drawn from `amino_acids` + N deaminated (Option A, the D-33 idiom),
  and copper binding it stoichiometrically (Cu(SR)₂, 1 Cu:2 thiol, H₂S-first). See D-45. The
  reductive-sulfur beat (H₂S + mercaptans, autolytic sources + copper fining) is now **complete**.
- ~~**SO₂-bound acetaldehyde protected from ADH (the D-28 free/bound RHS coupling)**~~ — **decided +
  IMPLEMENTED in D-47 (2026-07-06)**: `AcetaldehydeReduction` reduces only the free (unbound) share
  (`acidbase.free_acetaldehyde`), so dosed SO₂ *locks in* acetaldehyde (near-stoichiometric stranding;
  ~0.76× degradation slowdown at field doses — literature-grounded). Owner chose bake-in default-on;
  the D-22/D-28 "SO₂ readout-only" invariant is **intentionally retired** for sulfited runs (undosed =
  byte-for-byte D-27, no benchmark doses SO₂, carbon still closes, pH still not a charge actor). See D-47.
- ~~**SO₂-induced acetaldehyde over-production (the D-47 caveat's deferred "production half")**~~ —
  **decided + IMPLEMENTED in D-48 (2026-07-06)**: a total-SO₂-gated bump to `AcetaldehydeProduction`
  (`k_acet_so2_induced`), scoped to the **transient peak** after the data showed D-47 protection *alone*
  already meets/exceeds the field 0.39 mg/mg end-state slope (end state is capped by the D-28 binding
  equilibrium, not production). Driver is **total** SO₂ — free SO₂ is empirically inert on the peak
  (collapses to ~0 there). Carbon-exact borrow from E; exact undosed guard; magnitude speculative and
  unanchored. See D-48.
- ~~**Residual-nitrogen / satiation floor**~~ — **addressed in D-30 (opt-in cap) and RESOLVED in
  D-43 (2026-07-06): the "default-on N redesign" is declined.** A spike + a mass-balance argument
  (D-43) proved that **default-on residual *assimilable* N is Coleman-incompatible regardless of
  mechanism** — Coleman builds biomass by ~day 1.3, which pins external assimilable N to ~0 by then
  for every dose, so no biomass-preserving N model (two-pool, cell-quota, satiation) can widen the
  H₂S lever or leave a late-window residual without cutting biomass and breaking the Coleman sugar
  curve. The deferred note's two mechanisms have *opposite* Coleman-compatibility: a proline/
  non-assimilable split is Coleman-safe but **inert** (nothing reads proline; it does not feed the
  assimilable H₂S gate), while a residual-*assimilable*-N floor is inherently **opt-in**. Decision:
  keep the D-30 opt-in `carrying_capacity_gpl` cap as-is; do not build the refactor. The residual-N
  lever stays opt-in. If the H₂S cross-must lever is ever wanted default-on, the clean route is
  re-pointing the *H₂S gate* onto a dose-correlated proxy (an H₂S-model change), not the N model —
  see D-43 forks (a)–(d).
- **Packaged parameter-data access:** tests read YAML via filesystem path. If we
  ship a wheel that must read its own data, switch to `importlib.resources`.
- ~~The residual D-51 overshoot (1.15–1.45× the field 0.39 mg/mg slope, worst at high SO₂ dose)~~
  — **RESOLVED (accepted, not closed) in D-59 (2026-07-08).** Two independent research angles
  converged: an affinity-arithmetic argument (acetaldehyde ~99% bound at these bisulfite levels,
  100–370× tighter than its competitors — no binder pool of any capacity/affinity can free enough
  at high dose) makes a fourth binder pool structurally unable to help regardless of the field
  anchor; separately the field "0.39, linear to 200 mg/L" anchor turned out to be a cross-sectional
  survey regression across ~12 wines measured only to ~124 mg/L (not a within-wine titration),
  sitting inside the broader literature's own 0.2–0.5 mg/mg / ~1.2×-study-disagreement envelope.
  **Decision: accept 1.15–1.45× as a documented, structurally-explained model limit — no new binder
  pool.** Both cheap follow-ups **DONE in D-61 (2026-07-08):** the D-51 citation is reconciled (the
  slope equation is Marrufo-Curtido et al. 2022, *Foods* 11(3):476, verified against the primary
  source; Jackowetz & Mira de Orduña 2013 is the separate 237-wine concentration survey — fixed in
  the D-51 note + `acetaldehyde.yaml` provenance), and the in-model bound fraction at 200 mg/L SO₂
  reads **98.5–99.4%** at realistic acetaldehyde levels, confirming the affinity-arithmetic argument.
  See D-59, D-61.
- ~~The D-56 Varela 2004 fermentation-rate gap~~ — **D-56's mechanism-1 diagnosis was WRONG (stale
  note; already fixed in D-13) and mechanism 2 was substantially CLOSED in D-57 (2026-07-07)** by
  fixing a real, sourced bug (`k_prime_d`'s missing quadratic temperature scaling) instead of
  building a novel Bisson mechanism. Gaps narrowed: N=300 ~2.05x→1.91x (barely moves, confirming a
  genuine Coleman-vs-Varela cross-study difference, not a model defect); N=50 ~3.98x→2.23x; the
  N50/N300 ratio shortfall against Varela's real 4.12x fell from ~1.94x-too-small to
  ~1.17x-too-small. **What remains is a small, honestly-documented residual** (that ~1.17x), which a
  Bisson-sourced nitrogen-gated transporter mechanism could still chase, but D-57 judged it no
  longer clearly worth the calibration/validation-firewall risk at this size — owner's call whether
  to pursue further or accept it as a documented model limit. **D-59 (2026-07-08) sharpened this:**
  the real gap is qualitative, not the ratio — Varela's N50 case is a *stuck* fermentation (arrested
  at 16 g/L residual sugar, never dry) and the model always finishes dry. A firewall-safe mechanism
  candidate exists (Salmon 1989 sugar-transport catabolite inactivation; NOT Palma 2012, which D-59
  earmarked for validation use instead — see below), but both D-59 research angles recommend running
  a cheap internal diagnostic first (model N50 viable-biomass vs. Varela's measured cells) before
  building anything. **DIAGNOSTIC RUN in D-61 (2026-07-08) — verdict: DON'T build the Salmon
  mechanism.** Total biomass matches Varela at N50 (1.40 vs 1.5 g/L) so biomass is not the culprit;
  the fast-finish is a per-cell *rate* feature present even in-range (model has *less* biomass than
  Varela at both N levels yet finishes faster), and that rate is Coleman's own fitted rate (same
  strain as Varela, D-59 Finding 0) — so the gap is a Coleman-vs-Varela cross-lab difference, and
  the model already encodes Varela's arrest phenomenon as `X_dead` (catalytic-vitality loss). A
  Salmon per-cell term would be redundant and would calibrate the model away from its own fit
  source. **Accepted as a documented model limit.** See D-59, D-61.
- ~~`BrettSenescence` twin~~ — **RE-CONFIRMED DECLINED in D-58 (2026-07-07)** via two independent
  literature-research agents: no source shows Brett declining from elapsed time alone (every observed
  decline traces to SO₂, ethanol toxicity, or substrate exhaustion). D-40/D-52's "persists
  indefinitely" wording should be read as "no positive evidence for spontaneous decline without SO₂,"
  not literal immortality. See D-58.
- ~~Brett ethanol-toxicity death gate~~ — **IMPLEMENTED in D-58 (2026-07-08).** `BrettEthanolToxicity`
  (a new sibling `Process` to `BrettDeath`) plus a `BrettGrowth` upper wall, both driven by a shared
  threshold survival factor sourced at Barata et al. 2008's boundaries (onset ~14% v/v/110 g/L,
  ceiling ~14.5–15%/118 g/L). No SO₂ needed — the point of the mechanism. See D-58.
- ~~A second independent wine validation dataset (Palma 2012)~~ — **BUILT in D-60
  (2026-07-08).** `tests/benchmarks/test_validation_palma2012.py` digitizes CF (320 mg
  N/L) and LF (90 mg N/L) glucose+ethanol curves (strain PYCC 4072, genuinely independent
  of Coleman/Varela's Prise de Mousse lineage). Corroborates the D-56/D-57 N-sensitivity
  shortfall on an independent strain; the absolute CF/LF timing gap flips direction from
  Varela and is protocol-confounded (shaken-flask yield ~0.39 g/g vs the engine's ~0.48),
  cross-checked engine-faithful-to-Coleman at Palma's inputs. **RF (refeed) BUILT in D-62
  (2026-07-08):** the DAP-refeed rescue is reproduced, but the engine inverts Palma's
  within-study RF-vs-CF ordering (engine RF<CF ~108<138 h; Palma RF>CF ~117>72 h) — the
  same N-under-suppression gap via a dynamic intervention; all three Palma conditions now
  built. See D-60, D-62. ~~**Beer-side independent check still open**~~ — **PARTIALLY ADDRESSED
  in D-63 (2026-07-09).** No publicly-accessible independent in-regime dataset exists (its two
  richest candidates are its own fit sources); the accessible off-regime lager reconstruction
  (Reid 2021 / Speers 2003) is single-temperature (10 °C), which is confounded by organism +
  pitch. Built `test_beer_temperature_response.py` — an honest cross-regime Arrhenius *stress
  test* (engine's own apparent E_a ~55 kJ/mol sits in the literature yeast range; excludes the
  265 kJ/mol artifact), NOT a lager validation. **The confound-cancelling ratio test stays
  deferred on EVIDENCE:** owner obtained the paywalled Speers 2003 PDF mid-session and it was
  read in-source (2026-07-09) — temperature is brand-dependent, the ferments free-rise
  (non-isothermal), and Table I tabulates no temperature values, so it is not the controlled
  series the ratio test needs. Deferral now rests on the data structure, not access. See D-59, D-63.

## D-64 — Hop bittering → IBU: the boil isomerization is a sourced wort-side compile-seam calc (Malowicki closed form), the fermentation loss is the only Process, iso-alpha is off the carbon ledger, and a utilization coefficient is ADDED (not fitted) to avoid a 2× IBU overprediction

**Status: BUILT (2026-07-10).** The §3.3 "additives with clear mechanisms" beat, owner-selected off
the post-D-63 menu (over the sensory/OAV Tier-3 capstone, aging chemistry, and the deferred-tail
options). New: `parameters/data/hops.yaml`, `core/kinetics/hops.py`, `analysis.ibu_series`,
`scenario.schema.HopAddition` + three `Scenario` fields, the compile-seam boil calc, a beer-only
`iso_alpha` state slot + `IsoAlphaAcidLoss` Process, and `tests/test_hops.py` (20 tests). Full
suite green (704 = 684 prior + 20 new, ruff+mypy clean); `test_media.py`'s three structural
assertions (beer schema names/size, canonical units, expected process set) were updated to include
the new slot/Process — a structural reflection, not a benchmark weakening. Advisor consulted once
before writing (the shape was endorsed; five sharpening points applied — see below).

**The physics has two regimes, handled in the two places they belong.** Bitterness is *iso*-alpha-
acids (isohumulones), which do not pre-exist in the hop:
1. **The boil** (~373 K, 60–90 min, PRE-fermentation, no yeast): a CONSECUTIVE first-order reaction
   `alpha --k1--> iso-alpha --k2--> degradation`, both constants measured by Malowicki & Shellhammer
   2005 (*J. Agric. Food Chem.* 53(11):4434-4439, doi:10.1021/jf0481296) over 90–130 °C:
   `k1 = 7.9e11·exp(-11858/T)` min⁻¹ (Ea 98.6 kJ/mol), `k2 = 4.1e12·exp(-12994/T)` min⁻¹ (Ea 108.0
   kJ/mol), T in K. Modelled by the CLOSED-FORM intermediate `[iso]/[a0] = k1/(k2-k1)·(e^{-k1 t} -
   e^{-k2 t})`, evaluated **once at the compile seam** per hop addition and summed — NOT a boil ODE
   phase (running the boil through the integrator would drive the yeast-free wort at 373 K). This is
   the same wort-side-input treatment `initial_ph` gets (D-18): only the *result* (iso-alpha
   delivered to the fermenter) enters the state. At 60 min/100 °C the closed form gives 47.6% of
   alpha as iso-alpha, still on the RISING limb (k2<k1, peak ~3 h) — matching brewing practice.
2. **Fermentation** (the engine's native regime): `IsoAlphaAcidLoss` removes iso-alpha by adsorption
   onto viable yeast (`d(iso_alpha)/dt = -k_iso_alpha_loss·X_viable·iso_alpha`, X-gated so a crashed/
   racked beer stops losing bitterness) — the ~5–20% wort-to-beer drop. This is the *dynamic content*
   of the beat and the reason hops touch the ODE at all; a crash mid-ferment strands the bitterness.

**Off the carbon ledger (the accounting choice).** Iso-alpha-acids are exogenous (they arrive via
hops, mg/L scale) and touch only `iso_alpha` — never S/E/CO2/N. Like dosed SO2 (D-22) they are
absent from `total_carbon`/`total_nitrogen` (an unreferenced slot gets weight 0 in `conservation`),
so the whole beat leaves the carbon invariant **byte-for-byte unchanged** — asserted directly:
`test_hopping_leaves_total_carbon_byte_for_byte` runs a hopped and an unhopped beer and checks the
carbon *series* are identical to 1e-9. The fermentation loss is adsorptive removal of hop-derived
mass, not a conversion within the fermentation carbon budget.

**The load-bearing modeling decision — a utilization coefficient is ADDED, not fitted (advisor point
2).** Malowicki's kinetics describe the isomerization of *dissolved* alpha faithfully (~48% at 60
min/100 °C in the kettle), but finished-beer utilization is only ~25–30% (typical brewing texts;
Tinseth ~23% at SG 1.050). The ~2× gap is a chain of physical losses NOT in Malowicki's pure-buffer
numbers: incomplete extraction from the hop material, break/trub adhesion, foam loss, and the
kettle→fermenter transfer. Reporting raw kettle iso-alpha would OVERPREDICT finished IBU ~2× — a
correspondence-with-reality failure (prime directive #1), not an acceptable simplification. So a
lumped `hop_utilization_efficiency` (0.55, banded [0.4, 0.75]) multiplies the end-of-boil iso-alpha
down to the fermenter-delivered value; `IsoAlphaAcidLoss` (~13% on a typical primary) then carries
it to the finished value. **Set from literature-typical utilization, NOT fitted to Tinseth** — the
Tinseth comparison is an independent cross-CHECK (fit-vs-fit, §3.5), which keeps the validation
firewall intact (the D-17/D-57 discipline). Composed effective utilization = 0.476·0.55·0.868 ≈ 22.7%
vs Tinseth 23.1%; the canonical recipe (1 oz 5% AA, 60 min, 5 gal, SG 1.050) finishes at ~17.0 IBU
vs Tinseth ~17.3. `test_finished_ibu_is_in_the_tinseth_ballpark` checks three recipes within ~30%.

**Volume is a genuinely new scenario quantity (advisor point 4).** Hop *mass* (grams) → g/L needs a
wort volume, which the otherwise concentration-based (volume-agnostic) engine did not track. Added
`Scenario.batch_volume_liters` (required iff `hops` is non-empty — a `model_validator` enforces it)
plus `Scenario.hops: list[HopAddition]` (alpha_acid_percent, grams, boil_minutes) and
`Scenario.boil_celsius` (default 100; lower for a whirlpool/altitude, which slows isomerization via
the Malowicki Arrhenius). v1 uses ONE volume for boil and fermenter (kettle-loss/evaporation folded
into the efficiency) — a documented simplification.

**Tiers derive, they are not asserted (advisor point 5; D-1).** The boil constants are sourced/
measured → **plausible** (not validated: the mapping to real wort — extraction, gravity, hop form —
is an honest-mapping step). But the finished `iso_alpha` also reads the speculative
`hop_utilization_efficiency` and the speculative `IsoAlphaAcidLoss`, so parameter-tier propagation
caps the finished-IBU readout at **speculative** — verified by `tier_of("iso_alpha")` = SPECULATIVE
on a hopped run, VALIDATED on an unhopped run (the loss Process is disabled at the compile seam when
no hops are scheduled, the MLF/Brett isolability pattern, so the empty slot keeps its tier and pays
no flux).

**Isolability (prime directive #3).** `_HOPS_PROCESSES` is wired into the BEER medium only (wine has
no `iso_alpha` slot); hops on a non-beer medium is a loud `ValueError`, not a silently-ignored field.
An unhopped beer is byte-for-byte the prior beer core (iso_alpha starts 0, loss disabled).

**Sourcing note (advisor point 1 — BLOCKED until resolved).** The Malowicki constants were taken
from the paper/corroborating sources, NOT recall: the ACS abstract + a secondary review + internal
consistency (`exponent·R = stated Ea` for BOTH k1 and k2: 11858·8.314 = 98.6 kJ/mol, 12994·8.314 =
108.0 kJ/mol) triangulate the values. The open-access thesis (Oregon State) is a scanned image and
unreadable via fetch, but the abstract-level values are unambiguous and independently reproduced.
The advisor's loose "Ea ~50 kJ/mol" sanity guess was wrong (the measured value is 98.6); the advisor
explicitly said to trust the paper if units check out, and they do — a case of primary-source
evidence correctly overriding an advisor heuristic.

**DEFERRED (v1 scope, documented in hops.yaml):** (a) the gravity-dependence of utilization (higher
wort gravity lowers hop utilization — Tinseth's bigness factor; Malowicki's pure buffer has no
gravity term and no mechanistic sourced form is available, so the efficiency is gravity-flat, anchored
at moderate gravity — and the "higher gravity → lower utilization" directional property is NOT
claimed); (b) dry-hop / whirlpool (post-boil, sub-100 °C) bitterness; (c) pH- and hop-form (pellet
vs whole vs extract) dependence; (d) polyphenol / oxidized-alpha (humulinone) bitterness. Only kettle
iso-alpha bitterness is modelled. See milestone-2-tasks.md.

## D-65 — §3.3 acid/sugar adjustments: the last §3.3 additive, two compile-seam verbs (`add_acid` general over the D-18 acids, `add_sugar` = sucrose inverted to hexose), both pure state mutations booked as external flows — no new Processes

**Status: BUILT (2026-07-10).** The closing beat of §3.3 "additives with clear mechanisms" —
owner-selected as the natural continuation of D-64 (hop bittering). The other three §3.3 additives
were already built: SO₂ (D-22/D-28, `add_so2` in D-36), nutrient DAP (`add_dap`, D-36), hop
bittering (D-64). This lands the fourth — **acid/sugar adjustments (tartaric acid additions,
chaptalization)** — which the handoff brief calls "simple state mutations via events," and that is
exactly what they are: two new intervention verbs at the compile→core seam, **no new Processes, no
new physics in the ODE**, riding the D-35/D-36 external-flow ledger the sibling verbs already use.
Full suite green (717 = 704 prior + 13 new tests in `test_interventions.py`), ruff + mypy clean.
One `advisor()` before writing (design endorsed, gotchas applied) and three owner forks decided by
`AskUserQuestion` before any code (below). No source file outside the verb registry + one new
`additions.yaml` param was touched; every prior benchmark unchanged.

**Two verbs, both the "add a species to its slot, book the external flow" idiom (the `add_dap` +N /
copper mercaptan −C precedent):**

1. **`add_acid {acid, gpl}`** — dose a charge-active organic acid. General over the D-18
   `acidbase.ACID_STATE` set (tartaric/malic/lactic): `params` names the `acid` and dose `gpl`, the
   whole mass lands on that acid's state slot. Those slots are wine-only (D-18), so the verb is
   **wine-only by slot presence** — a beer scenario raises ("needs a 'tartaric' slot"). The
   **load-bearing modelling choice**: the dose is the *pure acid* (it brings its own protons, no
   counter-cation), so it is added to the acid slot but **NOT to `cation_charge`**. The D-18 charge
   balance then re-solves the SAME back-anchored strong cation against MORE anion, so **pH drops and
   TA rises — emergently**, straight out of the keystone, not scripted. (Potassium bitartrate, which
   *does* add a counter-cation, would be a different verb — deferred.) Each acid carries carbon
   (tartaric/malic C4, lactic C3, all weighted in `total_carbon`), so the dose is a **positive**
   carbon external flow (opposite sign to the D-45 copper mercaptan −C removal) and nitrogen-free;
   the crown-jewel identity `final == initial + Σ flows` still closes to machine precision.

2. **`add_sugar {sugar_gpl}`** — chaptalize (and beer priming/adjunct). The dose is **sucrose** by
   mass, and the verb **inverts it AT THE DOSE** (a state mutation, NOT a kinetic pool — yeast
   invertase is fast vs the ferment) to hexose-equivalent via the exact new
   `sucrose_inversion_mass_ratio` (~1.0526). The +5.26% over the sucrose mass is hydrolysis water
   (C₁₂H₂₂O₁₁ + H₂O → 2 C₆H₁₂O₆; sucrose is an isomer of maltose, M = 342.30) — the SAME
   di-/tri-saccharide mass gain beer's wort sugars already carry (D-8, `chemistry.HEXOSE_UNITS` /
   `M_WATER`). The hexose lands on the fermentable sugar slot: wine's single lumped `S`, or beer's
   **glucose** component *specifically* (found by name via `chemistry.sugar_species(schema).index
   ("glucose")`, never broadcast across the maltose/maltotriose slots — the advisor-flagged 3-vector
   trap). Fructose from the inversion lumps as glucose-equivalent — exact on carbon and mass since
   they are isomers. More sugar ⇒ higher finished ethanol/ABV once it ferments out (emergent, no
   explicit ABV term). Carbon is conserved through inversion (water is carbon-free), so the flow
   books exactly the sucrose carbon (a positive flow); nitrogen-free; ledger closes.

**Three owner forks, decided up front by `AskUserQuestion` (owner chose the more capable option on
all three, over the recommended lighter defaults):**
- **Acid verb shape** → *general* `add_acid {acid, gpl}` over the ACID_STATE set (NOT a tartaric-only
  verb): any charge-active acid slot can be dosed; `test_add_acid_is_general_over_the_charge_active_
  acids` exercises malic.
- **Sugar dose basis** → *sucrose with explicit inversion* (NOT dose-as-hexose): the verb models the
  invertase mass gain via the new stoichiometric ratio, rather than treating `sugar_gpl` as
  already-hexose.
- **Sugar scope** → *wine + beer* (NOT wine-only): beer priming/adjunct is real, so the verb targets
  the glucose slot explicitly in both media; `test_add_sugar_on_beer_targets_glucose_only` pins that
  the inverted hexose lands on glucose alone.

**The one new parameter — `sucrose_inversion_mass_ratio` in `additions.yaml`.** Value 1.0526,
**VALIDATED with a zero-width uncertainty band** (exact stoichiometry, never swept by the ensemble) —
the `dap_nitrogen_fraction` precedent (prime directive #2 admits no magic numbers even for exact
stoichiometry; the value travels with its derivation in the provenance notes). It is a unit-conversion
constant read only by the verb at the compile boundary, not physics in the hot loop. No other new
params — `add_acid` needs none (the acids are already carbon-weighted and pKa-sourced from D-18).

**No tier movement (unlike `pitch_mlf`).** Neither verb enables a Process; both touch inert slots (the
acid slots and `S` have no derivative-touching Process gated on them here). So no tier drags —
`test_add_acid_moves_no_tier` asserts the acid/sugar/ethanol tiers are byte-identical dosed-vs-undosed.
pH's tier is already the PLAUSIBLE-floored pKa tier (D-18), unchanged by a dose.

**Isolability & conservation (the discipline every verb inherits).** A scenario with no interventions
is byte-for-byte a plain run (the pre-existing `test_no_interventions_is_byte_for_byte_plain_simulate`
still passes). Carbon closes across the jump for both verbs as a *positive* external flow
(`test_add_acid_books_a_positive_carbon_flow_and_no_nitrogen`,
`test_add_sugar_books_positive_carbon_and_no_nitrogen`), the mirror of the copper mercaptan −C removal.
The concentration-model no-volume-change caveat (shared by every verb) applies — stated once, not
re-litigated.

**The pH headline is asserted directionally, not to a magnitude (advisor point).** `acidbase.py` claims
directional/slope fidelity for its concentration-based *apparent* pKa simplification, so
`test_add_acid_lowers_ph_and_raises_ta` asserts pH↓ + TA↑ at a realistic ~2 g/L tartaric dose within a
sane band (< 1.0 pH unit), not a tight pH-delta that would over-claim. pH does not feed back into the
yeast kinetics in v1 (D-18), so the S/E/X trajectories are identical dosed-vs-undosed — only tartaric
(and the derived pH/TA readouts) move, a clean isolation. D-46 (`solve_ph` is total) guarantees even
an extreme acid dose cannot crash the solver.

**DEFERRED (v1 scope):** (a) potassium bitartrate / K-tartrate additions (deacidification via a
counter-cation — a different, cation-moving verb); (b) a kinetic sucrose pool with an explicit
invertase Process (instantaneous inversion at the dose is an excellent approximation — invertase is
fast relative to the ferment — and honours the brief's "simple state mutation" framing); (c) direct
glucose/fructose dosing (the sucrose form is the standard chaptalization sugar; a `form` param is a
trivial future extension); (d) volume change on addition (the engine is concentration-based;
volume-tracking is the D-64 `batch_volume_liters` frontier, not extended here). **§3.3 is now
COMPLETE** — all four "additives with clear mechanisms" (SO₂, DAP nutrient, hops, acid/sugar) built.
The next-direction frontier (Tier-3 sensory/OAV, aging, or UX) is the owner's call.

## D-66 — Milestone 3 (Tier-3) opened: two calls before any code — (1) build the sensory/OAV readout layer FIRST, aging chemistry second (inverting handoff §6-step-5 order); (2) lumped aroma pools get a representative-compound threshold per lump

**Status: SCOPING (2026-07-10).** Owner selected the next-direction frontier left open at D-65:
Tier-3 — the handoff §4 "frontier" where the chemistry is real but "integrating it into a
trustworthy prediction is *not solved science*." This entry records the milestone-opening design
calls (the plan is `docs/plans/milestone-3-plan.md`, mirroring how M2 opened with D-18); **no code
yet** — Tier-3 is `speculative`, isolated, and must never perturb the validated core (prime
directive #3). One `advisor()` pass shaped both calls before writing; the second call was put to the
owner by `AskUserQuestion` (the load-bearing fidelity fork). Preconditions confirmed clean: §3.3 /
all of Tier-2 settled, 717 green, `media.py` fully mature (Brett/MLF/keto-acids/autolysis all
landed) — the plan can state "Tier-2 settled" as fact, not assumption.

**Call 1 — sensory/OAV FIRST, aging chemistry SECOND (invert handoff §6 step 5's "aging then
sensory").** Rationale: the sensory layer is a **pure readout** over aroma-active compounds the
model *already* tracks (`esters`, `fusels`, `diacetyl`, `acetaldehyde`, `h2s`, `ethylphenols` =
4-EP, `ethylguaiacols` = 4-EG, `mercaptans`), so it adds **no new ODE physics and zero risk to the
validated core**, and once built it becomes the **acceptance lens for aging** — every aging Process's
effect on the aroma profile is then immediately visible. Aging chemistry is the heavier piece (new
speculative RHS Processes on a years-scale phase, phase-based integration per handoff §7, scattered
parameter sourcing), so it is second, one Process at a time behind its own tests. The handoff order
is reference-not-gospel (CLAUDE.md); the owner's own framing ("sensory/OAV, aging chemistry") put
sensory first too. Not burned on an `AskUserQuestion` — the architecture and the owner's phrasing
already agreed.

**Call 2 — lumped aroma pools use a representative-compound threshold per lump** (owner decision via
`AskUserQuestion`, over the single-compound-only alternative). OAV = concentration ÷ per-compound
perception threshold, but `esters`/`fusels`/`mercaptans` are single lumped g/L pools mixing molecules
whose thresholds span ~3 orders of magnitude (isoamyl acetate vs ethyl hexanoate vs ethyl acetate).
Options weighed: **(a)** assign each lump one *named representative* compound's threshold — the
stand-in its `VarSpec` already names (fusels → isoamyl alcohol; mercaptans → the "methanethiol
stand-in" it literally is; esters → isoamyl acetate) — compute OAV uniformly, and carry **"assumes
fixed lump composition"** loudly in provenance; **(b)** compute OAV only for the single-molecule pools
(diacetyl, acetaldehyde, h2s, 4-EP, 4-EG) and treat the lumps as descriptor-qualitative; **(c)** split
the lumps into constituent esters — **rejected**, that is a *chemistry*-layer change to serve sensory,
which inverts the §4.2 cardinal rule (chemistry never depends on the sensory layer). **Owner chose (a)**
— keeps the dominant young-product aromas (esters, fusels) in the numeric readout; the honesty cost is
the fixed-composition assumption, flagged at the source (the §4.3 "don't let speculation borrow the
core's credibility" concern is answered by the loud provenance note + the speculative tier floor).

**Architectural decisions baked into the plan (beat 1 = OAV ratio):**
- **Placement:** a new top-layer package `fermentation.sensory`, sibling of `fermentation.analysis`;
  consumes a `runtime.Trajectory` + a threshold table, imported by **nothing lower** (one-directional
  rule; §4.2 cardinal rule).
- **Thresholds load DIRECTLY into the sensory layer, NOT through the compile seam.** Unlike
  `acidbase.yaml` / `vicinal_diketones.yaml` (merged into every `CompiledScenario` at `compile.py`'s
  `shared_files` *because a Process reads them*), **no RHS reads a perception threshold** — so a new
  `sensory.yaml` is loaded by the sensory module standalone, never merged into `param_values`. A
  stronger isolation than any Tier-2 readout: the chemistry never even sees the sensory params.
- **Tier floor:** every OAV output tier = `Tier.combine(chemistry_input_tier, SPECULATIVE)` →
  **speculative even over a validated input** (the sensory mapping is itself the canonical speculative
  case named in the `Tier` docstring). Enforces the §4.3 firewall at the API.
- **Matrix-specificity is a provenance requirement:** ethanol/matrix shift most odor thresholds, so
  each threshold's `conditions` records the matrix (wine ≠ beer ≠ water/model), and any fallback is
  flagged as a matrix gap.
- **`iso_alpha`/IBU excluded** — it is a *taste* (bitterness), already a direct mg/L→IBU readout
  (D-64), not an odor threshold; not shoehorned into an OAV.
- **Descriptor-space projection is DEFERRED** to a separate, even-more-speculative sub-beat (1b): the
  OAV *ratio* is a defensible sourced number; "OAVs → smells like leather and banana" is a further
  heuristic leap, fenced behind a swappable seam so beat 1a stays honest.

**Reading list (to source at build, all `speculative`):** Guth 1997, Francis & Newton 2005, Meilgaard
1975, Ferreira et al. 2000 (odor thresholds); diacetyl ~0.1 mg/L (lager) and 4-EP/4-EG ~425/110 µg/L
(red wine) from the spoilage literature already cited in the VDK/Brett beats.

**Next:** beat 1a build — `fermentation.sensory.oav` + `sensory.yaml` + tests, recorded at D-67.

---

## D-67 — Beat 1a built: the OAV sensory readout (`fermentation.sensory`), first Tier-3 code

**Date:** 2026-07-10. **Milestone 3 / Tier-3, first build beat** (the scoping is D-66). Ships
the sensory/OAV layer the D-66 plan named — `fermentation.sensory.oav` + `sensory.yaml` +
`tests/test_sensory_oav.py` — as a **pure, isolated, speculative readout**. 729 tests green
(717 → +12), `ruff`/`mypy` clean. One `advisor()` pass before writing shaped the build; its one
blocking catch (below) is folded in.

**What landed.**
- **`fermentation.sensory`** (new top-layer package, sibling of `analysis`): `oav_series(traj,
  thresholds, pool)` = `conc / threshold` over a `Trajectory` (dimensionless); `sensory_profile`
  → a `SensoryProfile` of **per-compound** `OAVReading`s (OAV, static descriptor, above-threshold
  flag, lumped flag, tier) at a chosen time; `oav_tier`, `medium_of`, `load_thresholds`,
  `AROMA_COMPOUNDS`. Aroma set is medium-specific (beer = 5 common pools; wine = those + 4-EP/
  4-EG/mercaptans), medium inferred from the schema signature (`iso_alpha`→beer, `tartaric`→wine).
- **`parameters/data/sensory.yaml`** — 13 matrix-specific perception thresholds (5 beer + 8 wine),
  all µg/L, all `speculative`, sourced: diacetyl/acetaldehyde/h2s/isoamyl-acetate(esters)/isoamyl-
  alcohol(fusels) in beer (Meilgaard 1975); the same 5 + 4-EP/4-EG (Chatonnet 1992) + methanethiol
  (mercaptans) in wine (Guth 1997 model wine, Goniak & Noble 1987 H₂S, Martineau 1995 diacetyl).
- **`units.convert`** — added `gpl_to_ugl`/`ugl_to_gpl` (the g/L↔µg/L boundary, D-3); OAV crosses
  the *scalar* threshold µg/L→g/L so both sides compare in canonical g/L.

**The advisor's blocking catch — the tier-floor test was vacuous as first sketched.** Every aroma
pool is produced by a speculative/plausible Process, so `traj.tier_map[pool]` is never VALIDATED
and `combine([anything, SPECULATIVE])` is trivially SPECULATIVE — a floor test over a real
trajectory would prove nothing. **Fix:** factored a **pure** `oav_tier(input_tier, threshold_tier)
= combine([input, threshold, SPECULATIVE])` and assert `oav_tier(VALIDATED, VALIDATED) is
SPECULATIVE` **directly** — the only way to show the mapping caps a validated input. A second
end-to-end test pins an *untouched* pool's trajectory tier to VALIDATED and confirms the profile
reading still reads speculative. The explicit `SPECULATIVE` term is documented as **not** redundant
with the threshold's own tier: the sensory *mapping* is speculative, so the floor must hold even if
a threshold were later mislabelled plausible.

**Isolation, stated explicitly (so it reads as a choice, not luck).**
- **Byte-for-byte green by construction:** nothing lower imports `sensory`; the readout adds no
  state slot / RHS / ledger entry; `sensory.yaml` is **not** in `compile.py`'s `shared_files`, so
  it never enters any `CompiledScenario.param_values` and cannot perturb the chemistry. The full
  729-test suite passing is the end-to-end proof.
- **Thresholds sit outside the D-24 ensemble sweep** — a *deliberate* consequence of the standalone
  load: `simulate_ensemble` samples only compiled-scenario params, so it does not propagate
  threshold uncertainty into an OAV band. Defensible (the OAV is already floored speculative);
  recorded here so it never later reads as an oversight.

**Sourcing discipline (advisor (2), applied).** `conditions` records the **measurement matrix**,
which is *not* the application medium: the wine `esters`/`fusels` thresholds are Guth **model-wine**
(10% ethanol) values, flagged as a **matrix gap** in `notes`; beer thresholds are matrix-matched
(measured in beer). Widest uncertainty bands on the 3 lumped representatives (matrix *plus* the
fixed-composition assumption); firmest on the single-molecule wine phenols 4-EP/4-EG (Chatonnet).
The "matrix matches medium" test checks **set selection** (which pools the profile reports), not
measurement provenance. The golden test (diacetyl at 2× threshold → OAV ≈ 2) is named to make clear
it validates **plumbing** (arithmetic + the unit crossing), not the threshold magnitude.

**Deferred (unchanged from D-66):** sub-beat **1b** descriptor-space projection (kept out so the
sourced-ratio layer stays honest — `descriptor` is a *static* per-compound label, never a synthesised
"smells like X"); then the aging-chemistry beats (§4.1) on a years-scale phase, each validated by
this OAV lens.

**Next:** beat 1b (descriptor projection) *or* open the first aging Process — owner's call at the
next batch.

---

## D-68 — Aging axis OPENED: ester hydrolysis chosen as the first §4.1 Process (scoping + owner forks; no RHS yet)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, second beat opened** (after D-67 shipped the OAV
sensory readout). Owner picked "the first aging Process (§4.1)" as the direction; this entry records
the scoping — the Process chosen, two owner forks, the confirmed phase-attachment mechanism, and one
carbon-closure crux surfaced at design time and deferred to the build (D-69). **No RHS written yet**
(mirrors D-66 scoping → D-67 build): the crux below materially affects fidelity and wanted an advisor
pass that was rate-limited this turn.

**Process chosen — ester hydrolysis (advisor-affirmed).** The only §4.1 candidate that needs **no new
extraction driver and no new state pool to start**: it acts on the `esters` pool already tracked,
exercises the aging-phase pipeline on tractable chemistry, and — the payoff — moves an OAV the D-67
lens already reads (young fruity esters fade with age). Consistency win: D-67's sensory representative
for `esters` is **isoamyl acetate**, an acetate ester, exactly the class that hydrolyses and fades on
aging — so "net ester decay" is coherent with a choice already on record. The heavier candidates
(oxidation needs O₂-ingress modelling; oak extraction / tannin–anthocyanin need new pools) come later,
one Process at a time, validated by the OAV lens.

**Owner fork 1 (the direction):** first aging Process (§4.1) over beat 1b (descriptor projection) or
pausing — via `AskUserQuestion`.

**Owner fork 2 (carbon routing) — FAITHFUL SPLIT → `fusels` + `Byp`** (via `AskUserQuestion`, over "new
inert aging-products pool" and "→ Byp only"). Conservation is back in force (unlike the D-67 readout,
this is the first aging RHS *on the carbon ledger*): the carbon released by a decaying ester **must** be
routed. Owner chose the literal chemistry — isoamyl acetate + H₂O → isoamyl alcohol (→ `fusels`) +
acetic acid (→ `Byp`) — accepting that it (a) emergently **raises the fusel OAV** and drifts **pH/VA**
up with age (both real aging phenomena), a mild §4.3 firewall tension since a speculative aging Process
then touches the plausible-tier pH readout; and (b) uses `Byp` (succinic, C4 diprotic) as a stand-in
for acetic acid (C2 monoprotic). Isolability (togglability) is preserved regardless.

**Phase attachment — CONFIRMED, reuses the existing reconfigure mechanism (no new integration infra).**
`simulate_scheduled` already segments the timeline and a `ScheduledEvent.reconfigure` callback mutates
the `ProcessSet` in place; `ProcessSet.enable`/`disable` exist (the D-35/36 event precedent, e.g.
`pitch_mlf`). So an aging phase attaches as a **post-fermentation scheduled segment**: a `begin_aging`
event enables `EsterHydrolysis` (off during ferment) over a long span with the solver free to take
large steps (the §7 multi-scale concern — do not integrate years at ferment resolution — is answered by
the segment restart + large `max_step`, not new machinery). Open sub-questions for D-69: yeast state
during aging (racked/yeast-gone vs on-lees — decides whether ferment Processes are disabled or idle);
and the scenario-level expression of "then age N months" (extend `duration_days` + an `age`/`begin_aging`
verb, the D-36 intervention precedent).

**The carbon-closure crux surfaced at design time (the D-69 build must resolve).** The chemistry ledger
(`core.chemistry`) has **no `isoamyl_acetate` species** — the `esters` pool is carbon-weighted as
**ethyl acetate** (C4: 2C ethyl + 2C acetyl; `_ESTER_SPECIES` in `byproducts.py`), `fusels` as isoamyl
alcohol (C5), `Byp` as succinic (C4). So the owner's fork-2 framing (isoamyl acetate, C7 → C5 + C2) and
the pool's *ledger* stand-in (ethyl acetate, whose literal hydrolysis alcohol is **ethanol**, not a
fusel) disagree. Carbon leaving `esters` per gram decayed is **ledger-fixed** at `rate·c(ethyl_acetate)`;
the open question is the **split ratio** of that carbon between `fusels` and `Byp`:
- **5:2** (isoamyl-acetate molecular ratio, matches the owner's stated reaction + the D-67 OAV
  representative) — but mixes stand-ins (pool mass = ethyl acetate, split = isoamyl acetate);
- **1:1** (ethyl-acetate-consistent: acetyl 2C : alkyl 2C — matches the pool's own ledger structure,
  a single documented stand-in: ethanol-carbon routed to `fusels` rather than `E`).
Both close carbon by construction (the split ratio only re-partitions a fixed released-carbon budget
between two trace pools, so it is second-order on outputs); the choice is a fidelity/consistency call to
settle with the advisor at build. **Proposed RHS form** (advisor's framing, carried forward):
`d(esters)/dt = −k·f_T·max(0, esters − esters_eq)` — **net decay toward a lower equilibrium, decay-only**
(the bidirectional reality — ethyl esters of fatty acids slowly *form* on aging while acetates hydrolyse
— is the deferred half; framed as "net decay toward a lower equilibrium," the same fixed-composition
honesty the D-67 sensory lump carries — **not** decay-to-zero, which over-strips). Arrhenius `f_T`
(warm aging degrades faster). Tier **speculative**; the ethyl-acetate-pool / isoamyl-acetate-reaction
mismatch documented loudly (the D-19 "bookkeeping stand-in, not a metabolic claim" precedent).

**Decomposition:** **D-69 = the `EsterHydrolysis` physics** — the Process + a new `aging.yaml` params
file + direct unit/conservation tests (tested via `ProcessSet`, the D-64 loss-Process pattern), split
ratio resolved with the advisor. **D-70 = the aging-phase scenario wiring** — the `age N months` verb +
reconfigure enable + the §7 slow-phase integration end-to-end. **Next:** D-69 build (settle the 5:2-vs-1:1
split with the advisor first).

## D-69 — `EsterHydrolysis` built: the first aging RHS (§4.1), 5:2 carbon split (advisor-settled)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, aging beat built** (the D-68 scoping → this build,
mirroring D-66 → D-67). Ships the first §4.1 aging Process — `fermentation.core.kinetics.aging.
EsterHydrolysis` + a new shared `aging.yaml` + `tests/test_aging.py` — as a **speculative,
isolable, on-ledger** aging RHS. **742 tests green** (729 → +13), `ruff`/`mypy` clean. One
`advisor()` pass before writing settled the deferred split crux (below); its three build "musts"
(conservation test, don't touch the esters weighting, `esters_eq` a positive parameter) are all met.

**What landed.**
- **`EsterHydrolysis`** (new `core/kinetics/aging.py`, the home for the aging axis): `d(esters)/dt
  = −k_ester_hydrolysis·f(T)·max(0, esters − esters_eq)` — first-order **net decay toward a lower
  equilibrium floor** (not decay-to-zero; below the floor the rate is 0, the reverse ester-formation
  half deferred), `f(T) = arrhenius_factor(T, E_a_ester_hydrolysis, T_ref)` giving the sourced
  warmer-ages-faster direction. **No fermentative-flux gate** (aging runs when the flux is zero — it
  is temperature/pool-driven), unlike every M2 producer. `touches = ("esters","fusels","Byp")`,
  tier **speculative**.
- **`aging.yaml`** (shared, medium-agnostic like `vicinal_diketones.yaml` — ester hydrolysis is a
  molecule/pH property, not biology): `k_ester_hydrolysis` (1e-4/h, half-life ~3 mo–2.6 yr band),
  `E_a_ester_hydrolysis` (60 kJ/mol, Q10~2), `esters_eq` (5 mg/L floor). All **speculative**; the
  sourced parts are the *form* (first-order approach to equilibrium — Ramey & Ough 1980; Marais 1978)
  and the *direction* (E_a>0). **Not yet in the compile seam** — `EsterHydrolysis` is off the ferment
  ProcessSet and enabled only in a post-ferment segment (D-70), so the tests load `aging.yaml`
  directly; byte-for-byte isolation of the ferment is thereby preserved (prime directive #3).

**The split crux — resolved 5:2 (advisor flipped my initial 1:1 lean).** D-68 deferred the fusels:Byp
carbon split to this build's advisor pass. My going-in lean was **1:1** (ethyl-acetate-consistent),
on a "single documented stand-in" argument. The advisor showed that argument is **illusory**: the
esters pool's ethyl-acetate *mass* weighting is fixed by D-19 and immovable regardless of split, so
1:1 buys no reduction in stand-ins — the split is the one free variable, and there is no "clean"
choice, only *which representative it honors*. The discriminator is **what the Process is FOR**: it is
a **sensory** Process (its whole D-68 reason to exist is to fade the ester OAV and raise the fusel OAV),
and D-67 already commits `esters`'→isoamyl acetate, `fusels`'→isoamyl alcohol. The coherent chemistry
connecting those two committed representatives is **isoamyl acetate → isoamyl alcohol (5 C) + acetic
acid (2 C) = 5:2**. 1:1's hidden cost: its alcohol product is *ethanol* (ethyl acetate's real alcohol),
routed into the isoamyl-alcohol-weighted `fusels` pool and read through the isoamyl-alcohol OAV — it
would **fabricate the fusel-aroma rise out of the wrong molecule**, bending the exact quantity the
Process exists to move. 5:2's cost is **narrative-only** and invisible to every conservation test (the
*debited* molecule is ethyl acetate, the *split* molecule isoamyl acetate — a stand-in seam this
Process **inherits** from D-19/D-67, not one it invents). Bending narrative honesty to preserve sensory
honesty is the right trade for a sensory Process; 5:2 also gives the stronger fusel-OAV rise the owner
asked for (5/7 vs 1/2). **D-68 delegated this call to the advisor pass, so it commits without kicking
back to the owner** — documented loudly in the Process docstring.

**Carbon closure (the D-68 "conservation is back in force" requirement).** The carbon leaving `esters`
per unit decayed is ledger-fixed at `rate·c(ethyl_acetate)`; that budget is split 5:2 and re-deposited
via each product pool's *own* carbon fraction (`fusels`→isoamyl alcohol, `Byp`→succinic), so
`total_carbon` closes to **machine precision for any split summing to 1** (the `esters→esters_gas`
transfer precedent, but C4→C5-partial+C4-partial across two differently-weighted pools). Verified
per-RHS (`abs=1e-15`) and over an integrated ~1-year wine aging segment *and* a beer multi-slot run.
The 5:2 split constants are **code-with-citation** (stoichiometry of the named stand-in reaction, like
the chemistry carbon counts), not empirical YAML params.

**§4.3 firewall tension — documented, owner-accepted (D-68 fork 2), not relitigated.** The
speculative-tier `EsterHydrolysis` touches `Byp`, which the *plausible*-tier pH/TA readout reads (the
acetic-acid product drifts VA/pH up with age — a real aging phenomenon the owner chose the literal
chemistry for). Isolability preserved (disable ⇒ drift vanishes). `Byp` is the succinic (C4 diprotic)
stand-in for acetic acid (C2 monoprotic), the same D-16 pool stand-in.

**Next: D-70** — the aging-phase scenario wiring: an `age N months` verb, the `begin_aging` reconfigure
that enables `EsterHydrolysis` over a long segment, `aging.yaml` into `compile.py`'s `shared_files`, and
the §7 slow-phase (large-`max_step`) integration end-to-end.

## D-70 — Aging-phase scenario wiring: the `begin_aging` verb + `EsterHydrolysis` into the compile seam (§4.1)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, aging beat wired end-to-end** (D-69 built the RHS; this
wires it into the *scenario* pipeline). Ships the `begin_aging` intervention verb, `EsterHydrolysis`
into both media (disabled at compile), `aging.yaml` into `compile.py`'s `shared_files`, and a
scenario-level test file `tests/test_aging_scenario.py`. **753 tests green** (742 → +11 aging-scenario
tests, incl. a beer-path smoke test; the five bare-wine-RHS isolability fixtures gained `aging.yaml`),
`ruff`/`mypy` clean. One
`advisor()` pass before writing reframed the work around the deferred D-68 sub-question below.

**The advisor reframe — the load-bearing fork was NOT verb naming, but *what runs during aging*.**
D-68 deferred "yeast state during aging"; `begin_aging` *enables* `EsterHydrolysis` but *disables*
nothing, so the aging segment runs the **full** wine/beer set, not the clean `X=0, S=0, only-
EsterHydrolysis` envelope D-69 was tested against. The advisor named the one Process that could
confound the aging ester signal: **`EsterVolatilization`** (also moves `esters`, Arrhenius-driven). I
verified its RHS: it is **fermentative-flux-gated** (`_fermentative_flux_shape(y, …); if flux <= 0: return`),
and `fermentative_flux_shape` returns 0 when sugar OR biomass is 0. So at dryness (`S ≈ 0`) it — and
`EsterSynthesis`, `FuselAlcoholsEhrlich`, and the `Byp` uptake routing, **every** producer of the three
aging pools — is quiescent. **The aging ester/fusel/Byp signal is therefore unconfounded: only
`EsterHydrolysis` moves those pools during a post-dryness aging segment.** This settles the deferred
call as **Stance A** (leave the ferment set on; the aging effect emerges) over Stance B (reconfigure
disables the ferment set) — no need to disable anything, because the flux gate already does it, and the
one non-flux-gated draw across aging (`EthanolInactivation`, X→X_dead, carbon-neutral) drives the state
*toward* the D-69 `X≈0` envelope. Recorded as a first-class invariant in the code + a test, not left implicit.

**What landed.**
- **`EsterHydrolysis` wired into both media** (`_AGING_PROCESSES = (EsterHydrolysis,)` in `media.py`),
  medium-agnostic like the shared VDK/H₂S kinetics (ester hydrolysis is a molecule/pH property, and
  `esters`/`fusels`/`Byp` exist in both schemas) — but **disabled unconditionally at the compile seam**.
  Unlike the pitch-gated MLF/Brett tuples (which can co-inoculate at t0), aging is **inherently
  post-ferment** — there is no aging at t0 — so there is no t0-enable path; the *only* way to turn it on
  is a `begin_aging` event. Disabled ⇒ skipped by `active`/`tier_of`/strict, so an un-aged scenario is
  byte-for-byte the pre-aging core and `esters`/`fusels`/`Byp` keep their pre-aging tier (prime directive #3).
- **`begin_aging` verb** (`_verb_begin_aging`): the `pitch_mlf` reconfigure pattern **minus the state
  mutation** — a pure phase switch that `ps.enable("ester_hydrolysis")` at its `day` and injects/removes
  no mass (aging inoculates nothing). Takes **no params**; guards that the aging params are loaded (the
  `add_dap`/`additions.yaml` discipline) so a caller-supplied `parameter_paths` without `aging.yaml`
  fails loudly at compile, not as a bare `KeyError` mid-integration.
- **`aging.yaml` into `shared_files`**: every compiled scenario now carries `k_ester_hydrolysis`/
  `E_a_ester_hydrolysis`/`esters_eq` — inert (read by nothing) until `begin_aging` fires.

**Verb design — bare `begin_aging`, aging span via `duration_days` (advisor-endorsed over `age {months}`).**
The aging span is expressed by `duration_days` (put `begin_aging` at the ferment/aging boundary day and
extend the duration to cover the tail) — the `duration_days`-is-the-single-span-source invariant stays
clean, zero schema change. The `age {months}`-as-*intervention* alternative was **rejected**: an
intervention at `day == duration_days` is rejected by `_compile_interventions`' at/beyond-duration check,
so that framing fights an existing invariant. If "N months" ergonomics are wanted later, the clean shape
is a *top-level* `age_months` field the compiler uses to extend `t_span` and auto-insert the `begin_aging`
event — deferred as reversible, low-stakes sugar.

**§7 slow-phase integration — no new machinery, exactly as D-68 predicted.** `simulate_scheduled`
already segments the timeline; the BDF solver re-initialises its order at the `begin_aging` breakpoint,
and with the fermentative flux gone at dryness it takes large steps across the quiescent aging segment
(default `max_step=∞`). The end-to-end wine ferment→age run (30 d ferment + 150 d warm aging, the **full**
wine set active) integrates to `success=True` and `total_carbon` closes end-to-end (`begin_aging`
mutates no state ⇒ **no external flow** ⇒ the plain `final == initial` invariant, verified). The one real
wrinkle is **output resolution, not accuracy**: a whole-span default `t_eval=linspace(0, span, 200)`
under-samples the ferment; integration stays accurate (dense output), and callers wanting a fine ferment
curve pass `t_eval` — flagged so a coarse ferment plot never reads as a bug.

**Tier travels across the reconfigure (D-35 min-combine).** `EsterHydrolysis` is enabled only for the
aging back half, but `simulate_scheduled` min-combines the per-segment tier maps, so the speculative
aging Process drags `esters`/`fusels`/`Byp` to **speculative for the whole run** — a run is only as
trustworthy as its least-trustworthy segment. No `KeyError` risk: the aging params ride in every
`tier_map` now, and a disabled Process's `reads` are never consulted.

**Regression surface (both changes perturb defaults, neither perturbs trajectories).** Adding
`EsterHydrolysis` to the medium factories bumps the process-set *membership* (`test_media.py`'s
`EXPECTED_PROCESSES` gained `AGING_PROCESSES`); adding `aging.yaml` bumps every `ParameterSet`'s keys.
Both are count/contents changes only — disabled ⇒ skipped, so every un-aged *trajectory* stays
byte-for-byte. **Next:** beat 1b (descriptor projection) or the next §4.1 aging Process (oxidation / oak
extraction), each on the same `begin_aging` segment, validated by the D-67 OAV lens.

## D-71 — `OxidativeAcetaldehyde` built: the oxidative aging axis opens on a dissolved-O₂ pool (§4.1)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, the second §4.1 aging Process and the head of the OXIDATIVE
sub-axis** (D-69/D-70 built the *hydrolytic* half). Ships `OxidativeAcetaldehyde` (`core/kinetics/aging.py`),
a new dissolved-oxygen state slot `o2` (both media), an `add_oxygen` dosing verb, three `aging.yaml`
oxidation params, and oxidation tests in `test_aging.py` + `test_aging_scenario.py`. **768 tests green**
(752 → +16), `ruff`/`mypy` clean. One `advisor()` pass before writing settled the design crux, and the
**O₂-pool-vs-unbounded fork was put to the owner** (per "surface design decisions before building") — the
owner chose the dissolved-O₂ pool, so this beat *opens an axis*, not just a leaf Process.

**The chemistry (owner-endorsed as the right first oxidation Process).** As a finished wine/beer takes up
oxygen (bottle ingress, micro-oxygenation, barrel), O₂ oxidises ethanol → **acetaldehyde** — the
'sherry'/bruised-apple/nutty **oxidised** note. Like `EsterHydrolysis` (and per the D-68 selection
criterion), it moves an OAV the **D-67 lens already reads** (the same `acetaldehyde` pool the D-27 buffer
fills, 'green apple' fresh vs 'oxidised' when it climbs) and needs **no new aroma pool** — the one new slot
is the `o2` *substrate*.

**The advisor crux — O₂, not ethanol, is the rate-limiting reactant (this set the whole design).** My
first instinct was a rate first-order in ethanol (mirroring `EsterHydrolysis`). The advisor caught this as
a **fidelity defect**: ethanol sits at ~100 g/L, essentially constant across aging, so a rate first-order
in ethanol is a *constant rate in disguise* — acetaldehyde would rise **linearly and unbounded**, pinning
the kinetic limit on the wrong species. Mechanistically it is **coupled oxidation** (Wildenradt & Singleton
1974): O₂ oxidises o-diphenols → quinones + H₂O₂, then H₂O₂ oxidises ethanol → acetaldehyde — so O₂ is both
the driver and the natural bound. Making the rate **first-order in a finite `o2` pool** gives the correct
**saturating** behaviour (acetaldehyde plateaus as the O₂ charge is spent — bottle-aging reality). The
phenolic catalyst is folded into `k_ethanol_oxidation` in v1 (a documented lump; no general phenol pool
tracked).

**The fork put to the owner — dissolved-O₂ pool (Approach B) vs unbounded ethanol-first (Approach A).**
Because "oxidation done right" is *bigger* than `EsterHydrolysis` (a new state slot + a dosing verb + a
yield param = the **foundation of an O₂ sub-axis**, not a leaf Process), and per the owner's
"discuss-disagreements / surface design decisions before building" norm, the choice was surfaced. **Owner
chose B.** Rationale: correct driver + saturating bound, *and* the `o2` pool is the **shared substrate** the
whole future oxidative sub-axis (phenolic browning, Strecker, SO₂ consumption) will draw down — build it
now and those slot in as extra O₂ sinks; build ethanol-first and the foundation gets redone.

**What landed.**
- **`o2` state slot** (`_common_specs`, both media, default 0, g/L): the dissolved-oxygen aging substrate.
  **Carbon-free and off EVERY ledger** — `total_carbon`/`total_mass`/`total_nitrogen` weight only their
  explicitly-named pools, so `o2` contributes 0 to each with no registration (the `h2s`/`iso_alpha`
  precedent). `M_O2` added to `chemistry.py` as a plain constant (like `M_WATER`), used only for the
  g/L-O₂ → moles conversion that sets the yield.
- **`OxidativeAcetaldehyde`** (`_AGING_PROCESSES`, both media): `d(o2)/dt = −k_ethanol_oxidation·f(T)·[o2]`
  (first-order in O₂, Arrhenius warmer-faster), `d(acetaldehyde)/dt = +y_acetaldehyde_per_o2·(r_O2/M_O2)·
  M_acetaldehyde`, `d(E)/dt = −that·M_ethanol/M_acetaldehyde`. The `E → acetaldehyde` transfer is the
  **clean reverse of the D-27 reduction** (both C2, mole-for-mole), so `total_carbon` closes to **machine
  precision**; the standing E↔acetaldehyde mass gap is scoped out (`total_mass` = `{S,E,CO2}` is never
  asserted on an aging run). During aging `X=0`, so `AcetaldehydeReduction` (viable-X-gated) is inert —
  oxidation does not fight it, acetaldehyde accumulates.
- **The whole O₂ flux is consumed, only a *yield* becomes acetaldehyde.** `y_acetaldehyde_per_o2` (~1 mol/mol,
  banded 0.5–2, below the mechanistic max) — the remainder is the oxidative power spent on **unmodeled
  sinks** the future sub-axis will claim. Because O₂ is carbon-free, "spending" it without tracking every
  product is *not* a conservation violation; the carbon that does move (into acetaldehyde) is borrowed
  carbon-exactly from `E`. Sanity: ~40 mg/L cumulative O₂ × yield ≈ 55 mg/L acetaldehyde (fresh ~10–40,
  oxidised ~100–300) — in range, verified against literature, not hardcoded.
- **The seam the NEXT oxidative Process inherits (flagged for the next author).** `k_ethanol_oxidation`
  is presently the **total** O₂-depletion rate — `OxidativeAcetaldehyde` alone drains the *entire* `o2`
  flux (`d(o2)/dt = −k·f(T)·o2`), with the sub-unity yield absorbing "unmodeled fate." When browning /
  Strecker / direct-SO₂-consumption Processes are added, they must **not** each independently drain the
  full pool (that would over-consume O₂). The clean refactor at that point: make `k_ethanol_oxidation` the
  *ethanol-oxidation share* of a common O₂-depletion rate and have each O₂ consumer draw its own share, so
  the pool depletes once across all sinks. Recorded here so the seam is explicit, not a surprise.
- **`add_oxygen` verb** (the `add_so2` pattern): doses `o2_mgl` → g/L onto the `o2` slot, carbon-free. One
  dose = a bottle's ingress; repeated = micro-ox/barrel. The runtime books an external flow for the mutate
  delta, but it is **carbon/nitrogen-free** (o2 off every ledger), so the single-run carbon ledger still
  closes with no correction term.
- **`begin_aging` now gates BOTH aging Processes** (`_AGING_GATED_PROCESSES = (EsterHydrolysis,
  OxidativeAcetaldehyde)`): one tuple drives the enable (verb reconfigure) and the compile-seam disable, so
  they stay symmetric as the axis grows; the param guard covers both Processes' params.

**Isolability — a second gate makes reductive aging free.** Oxidation is inert at `o2 = 0` (exact guard),
so a `begin_aging` run **without** `add_oxygen` is purely **reductive** aging (screwcap/inert — a real
case) — byte-for-byte the `EsterHydrolysis`-only aging. A test pins that acetaldehyde ends exactly where
the un-aged run leaves it (the oxidation Process cannot move acetaldehyde without O₂). An un-aged run stays
byte-for-byte the pre-aging core (both aging Processes disabled at compile).

**§4.3 firewall / tier.** Speculative in FORM (Tier-3 frontier); the oxidation *form* (O₂-limited,
warmer-faster) is sourced, magnitudes are estimates. `o2`/`acetaldehyde`/`E` floor at speculative when the
Process is enabled (non-vacuous — proven for all three). `acetaldehyde`'s tier was **already** speculative
(the D-27 buffer), so the only *new* tier consequence is `o2` going speculative in an aged-with-oxygen run.

**Regression surface.** The new `o2` slot bumped the schema golden tests (`test_media.py`: SHARED +`o2`,
wine size 39→40, beer 21→22) and `EXPECTED_PROCESSES` gained `oxidative_acetaldehyde`; both count/contents
changes only — every default (un-aged) *trajectory* stays byte-for-byte (the slot is 0, the Process
disabled). **Next:** the next oxidative sub-axis Process drawing the same `o2` budget (phenolic browning /
Strecker / direct SO₂ consumption), oak extraction, or the deferred beat 1b (descriptor projection) — each
on the `begin_aging` segment, validated by the D-67 OAV lens.

## D-72 — `SulfiteOxidation` built: SO₂ scavenging is the first sink on the shared O₂ budget (§4.1)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, the third §4.1 aging Process and the first *sibling* on the
D-71 oxidative sub-axis** — the first O₂ sink to claim its share of the shared `o2` budget. Ships
`SulfiteOxidation` (`core/kinetics/aging.py`, **wine-only**), a new `bisulfite_so2_at_ph` helper in
`acidbase.py`, two `aging.yaml` params + one code-with-citation stoichiometry constant, and SO₂-oxidation
tests in `test_aging.py` + `test_aging_scenario.py`. **782 tests green** (+13 SO₂-oxidation, full suite incl.
benchmarks confirmed), `ruff`/`mypy` clean. **Two
`advisor()` passes before writing** — one on the design, one *reconciling a chemistry-species correction I
raised against the advisor's own earlier framing* — and the axis + rate-form fork was **put to the owner**
(who chose SO₂ consumption + the bilinear form). The reason this was the right first pick: it reuses the
existing `so2_total` pool (no new aroma pool), is non-regressive, and delivers a **celebrated wine-chemistry
threshold for free**.

**The chemistry — "SO₂ protects until exhausted, then acetaldehyde climbs."** SO₂ is wine's antioxidant
because **bisulfite (HSO₃⁻) is a faster O₂ scavenger than ethanol**. Both `SulfiteOxidation` and the D-71
`OxidativeAcetaldehyde` draw down the *same* `o2` pool, so `ProcessSet` summing splits the O₂ between them by
their rates: the fraction reaching acetaldehyde is `k_eth / (k_eth + k_so2·[HSO₃⁻])` — small while free SO₂
lasts (O₂ diverted to SO₂, oxidative acetaldehyde suppressed), → 1 once SO₂ is spent (acetaldehyde climbs).
End-to-end verified: at ~40 mg/L O₂, 0/30/100/300 mg/L SO₂ give ~55/46/22/0.7 mg/L acetaldehyde, with SO₂
consumed at the classic **~4 mg SO₂ per mg O₂** mass rule — all emergent, nothing extra built for the
diversion itself.

**The D-71 "refactor needed" prediction turned out UNNECESSARY (the key design finding).** D-71 flagged a
seam: `k_ethanol_oxidation` is presently the *total* O₂-depletion rate (`OxidativeAcetaldehyde` drains the
whole flux), and it warned that adding a second O₂ sink would over-consume O₂ unless `k_ethanol_oxidation`
were refactored into an "ethanol share of a common rate." **The advisor showed this refactor is a phantom
for a *substrate-gated* sink.** SO₂ oxidation is gated on its own substrate (`so2_total > 0`), so: (a) with
no SO₂ dosed it contributes **byte-for-byte zero** → D-71's reductive *and* oxidative curves are unchanged
(nothing to regress against); (b) with SO₂ present, the two first-order-in-O₂ rates simply **sum**, which is
*physically correct* — competing reactions split a finite pool by `kᵢ/Σk`, and O₂ (off every ledger) is
consumed exactly once, so summing is **not** double-counting. So I did **not** refactor D-71; the sub-axis
grows by *adding gated sinks*, not by re-partitioning a shared rate. (The refactor only becomes real for an
*always-on* sink like phenolic browning — which is exactly why those are worse first picks.)

**The advisor reconcile — bisulfite, not molecular SO₂ (a correction I raised against the advisor).** The
owner-facing fork I wrote (following the advisor's first pass) framed the bilinear driver as "molecular
SO₂." Digging into `acidbase.py` I found its own `bisulfite_fraction` docstring already names **HSO₃⁻ "the
reactive nucleophile"**, and the primary literature (Danilewicz) is explicit that molecular SO₂·H₂O is the
reactive *antimicrobial* form while **bisulfite is the reactive *antioxidant*** (the reducer of o-quinones
and scavenger of H₂O₂). I surfaced the conflict in a second advisor pass; the advisor **took the correction**
(its "molecular is the reactive form" parenthetical conflated the two). This corrected a *label in the
option I presented*, not the owner's actual decision (bilinear-over-gated = faithful-over-simple still
holds), so per the "adapt on primary-source contradiction" norm I proceeded without a blocking re-ask and
flagged it visibly. Net fidelity *gain*: `bisulfite_fraction` is ~0.94–0.99 across wine pH (mild pH
coupling), but a **stronger** coupling enters through *free* SO₂ — as `OxidativeAcetaldehyde` makes
acetaldehyde that binds SO₂ (D-47), free SO₂ falls and this scavenging **self-throttles**. Oxidation erodes
SO₂'s protective capacity two ways (oxidative removal here + D-47 binding) — the emergent feedback the
bilinear form buys over a plain SO₂-presence gate.

**What landed.**
- **`SulfiteOxidation`** (wine-only, `_OXIDATIVE_SO2_PROCESSES`, wired into the *wine* medium only like
  `_MLF_PROCESSES`/`_BRETT_PROCESSES`): `d(o2)/dt = −k_so2_oxidation·f(T)·[o2]·[HSO₃⁻]` (bilinear, Arrhenius
  warmer-faster), `d(so2_total)/dt = −2·(r_O2/M_O2)·M_SO2`. Touches only `o2`/`so2_total` — **both off every
  ledger** (no sulfur ledger; oxidising SO₂ to untracked sulfate moves nothing conserved), so no conservation
  term and nothing asserted. Wine-only because `so2_total` + the acid/cation pH slots are wine-only (D-18);
  an `SO2_STATE_KEY not-in-schema` guard makes it a hard no-op on beer besides.
- **Wine-only, not the shared `_AGING_PROCESSES`.** Unlike `EsterHydrolysis`/`OxidativeAcetaldehyde` (both
  media), this reads wine-only state, so it follows the MLF/Brett wine-only wiring. It still rides the aging
  gate: added to `_AGING_GATED_PROCESSES`, disabled at compile, enabled by `begin_aging` (both loops guard
  `name in process_set`, so listing a wine-only Process there is beer-safe).
- **`bisulfite_so2_at_ph`** (`acidbase.py`): `free_SO₂ · bisulfite_fraction(pH)` at an already-solved pH —
  the reactive antioxidant driver, mirroring `molecular_so2_at_ph` (the antimicrobial one), using *free* SO₂
  (bound bisulfite is already spent). Solves pH once via the hot-loop `_at_ph` discipline.
- **Stoichiometry as a code constant, not a parameter.** `_SO2_PER_O2 = 2.0` (mol SO₂ per mol O₂) — the
  Danilewicz coupled-oxidation mechanism spends one bisulfite reducing the o-quinone and one scavenging the
  H₂O₂ per O₂, = the classic ~4 mg-SO₂-per-mg-O₂ mass rule (`2·M_SO2/M_O2 = 4`). Reaction stoichiometry, so a
  code-with-citation constant like the chemistry carbon counts, not an uncertain YAML magic number. Distinct
  from D-47 **binding** (reversible free↔bound repartition of `so2_total`): this **oxidises** SO₂ to sulfate
  and removes it, so the two do not double-count.
- **`aging.yaml`** gains `k_so2_oxidation` (0.2 L/(g·h), banded — the load-bearing claim is the ORDERING,
  bisulfite out-competing ethanol for O₂, not the magnitude) and its own `E_a_so2_oxidation` (a separate
  param from the ethanol one per prime-directive #2 — distinct reaction, distinct provenance). Both
  speculative.

**Isolability + tier.** Doubly substrate-gated: inert at `o2 ≤ 0` **or** `so2_total ≤ 0` (both return
byte-for-byte zero and skip the pH solve), so a reductive (no-O₂) *or* an unsulfited aging is exactly the
case without this Process, and an un-aged run stays byte-for-byte the pre-aging core. Speculative in FORM
(Tier-3 frontier; the oxidation *form* — O₂-limited, bisulfite-driven, warmer-faster, 2:1 — is sourced, the
rate *magnitude* an estimate); `o2`/`so2_total` floor at speculative when enabled (non-vacuous).

**§4.3 firewall.** Decrementing `so2_total` nudges the *plausible*-tier molecular-SO₂ / antimicrobial and
pH/free-SO₂ readouts — the same accepted precedent as `EsterHydrolysis → Byp → pH` (D-68 fork 2). Isolable
(disable the Process and the drift vanishes); documented, owner-precedented.

**Regression surface.** `test_media.py`'s `EXPECTED_PROCESSES[wine]` gains `sulfite_oxidation` (a new
wine-only entry); beer is untouched. Every default (un-aged / no-SO₂ / reductive) trajectory stays
byte-for-byte (the Process is disabled at compile, and gated to zero without both substrates). **Next:** the
remaining O₂ sinks (phenolic browning / Strecker — each an *always-on* sink that WOULD need the D-71 rate
refactor, so build with that in view), oak extraction (a separate axis, no O₂), or the deferred beat 1b
(descriptor projection) — each on the `begin_aging` segment, validated by the D-67 OAV lens.

## D-73 — O₂ sub-axis reworked for *always-on* sinks: `k_ethanol_oxidation` is a **share**, not the total (§4.1)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, an enabling rework of the D-71 oxidative sub-axis — no new
Process, no new state, no new verb.** Reworks how `OxidativeAcetaldehyde` (D-71) accounts for the shared
`o2` budget so the next **always-on** O₂ sink (phenolic browning, Strecker degradation) can be added without
double-counting. Redefines `k_ethanol_oxidation` from the *total* O₂-depletion rate → the **ethanol-oxidation
share**, reframes `y_acetaldehyde_per_o2` as the route's *true* per-O₂ yield (re-baselined 1.0 → 1.5), and
rewrites the "whole-flux / unmodeled-sinks" docstrings + `aging.yaml` provenance. **782 tests green**
(byte-for-byte on structure; only the `y` value moved, and every aging test reads `y` from params or asserts
*relative* magnitudes, so all stay green), `ruff`/`mypy` clean. **One `advisor()` pass before writing** (it
caught the validation hole and the yield-coherence trap below), and the **scope fork was put to the owner**
(per "surface design decisions before building"): *seam-prep only* vs *re-baseline the numbers now* — the
owner chose **re-baseline now, with a note, until browning is built**.

**Why this is a defect fix, not a rename (the load-bearing point).** The two framings were already
**inconsistent in-tree**: D-71's docstring called `k_ethanol_oxidation` "the *total* O₂-depletion rate"
(`OxidativeAcetaldehyde` drains the whole flux), but D-72's `k_so2_oxidation` provenance *already* compared
against it as "the ethanol *route* rate that SO₂ out-competes 8–16×" — i.e. a **share**. "Ethanol share"
reconciles them. So this is not preparing-for-the-future cosmetics; it removes a live contradiction.

**The architecture — keep it additive; do NOT build a shared-total-rate.** Each O₂ consumer owns its own
rate constant (first-order or bilinear in `[o2]`); `ProcessSet` sums them, so the pool depletes **once** and
the O₂ splits among the sinks by `kᵢ / Σk`. This is exactly the pattern `SulfiteOxidation` established at
D-72 — D-73 only extends it to *always-on* sinks and fixes the naming/magnitude so they compose. The
tempting alternative — a single `k_o2_total` with per-sink fractions — was **rejected**: it would couple the
Processes (each needing the others' rates), break the independent-derivative contract, and can't represent a
substrate-gated sink whose share varies with its substrate. No new machinery; the summing already exists.

**The yield-coherence fix (the advisor's must-fix).** The old `y_acetaldehyde_per_o2 = 1.0` carried a
*muddled* rationale — "kept below the mechanistic max to leave headroom for unmodeled sinks." But cutting the
yield never freed O₂ for anything: the O₂ **partition lives in the rate constants, not in `y`**. Under the
clean reframe, `y` is the ethanol route's *own* per-O₂ conversion (competing fates are now explicit sibling
Processes), so it should sit near the coupled-oxidation yield, not at the floor. Re-baselined **1.0 → 1.5**
(mid-to-upper of the sourced ~1–2 band: the H₂O₂ arm reliably gives ~1 acetaldehyde per O₂, the catalytic
o-quinone arm adds a partial second equivalent). Sanity re-anchored: ~40 mg/L cumulative O₂ × 1.5 ≈ **82
mg/L** acetaldehyde (fresh ~10–40, moderately oxidised ~100–300) — in range.

**The owner's "with a note, until browning is built" — the interim caveat, recorded loudly.** Re-baselining
`y` *upward* while ethanol oxidation is still the **sole always-on** O₂ sink means it transiently receives the
whole always-on O₂ flux → aged-with-O₂ acetaldehyde is an **upper estimate**, higher than before this rework.
This is accepted, on the owner's call, *with the caveat flagged in three live spots* (the `OxidativeAcetaldehyde`
docstring, the `y_acetaldehyde_per_o2` provenance, and the module docstring): when the always-on browning /
Strecker sinks land, `k_ethanol_oxidation` is **reduced to its true share** (their sum holds the empirical
total O₂-depletion timescale — the anchor) and the acetaldehyde partitions down. `k_ethanol_oxidation`'s
**value stays 5.0e-4** now (its "share" currently *equals* the "total", the sole-sink identity), so O₂
*depletion* is unchanged; only the acetaldehyde *yield* moved.

**The acceptance test is a worked drop-in, not pytest (the advisor's hole).** "782 byte-for-byte green"
cannot verify this rework — nothing functional changed and **no always-on sink exists yet** to exercise the
new property. The seam is closed iff Browning slots in cleanly *on paper*. It does:

```python
class PhenolicBrowning(Process):          # always-on O₂ sink, wine (o-diphenols are a wine pool)
    name = "phenolic_browning"; tier = Tier.SPECULATIVE
    touches = ("o2", ...)                  # see "where the O₂ goes" below
    def derivatives(self, t, y, schema, params):
        d = schema.zeros(); o2 = float(y[schema.slice("o2")][0])
        if o2 <= 0.0: return d
        f_t = arrhenius_factor(..., params["E_a_browning"], params["T_ref"])
        d[schema.slice("o2")] = -params["k_browning"] * f_t * o2   # its OWN share, first-order in o2
        ...                                                        # + its product term
        return d
```

- **No double-count.** `ProcessSet` sums `OxidativeAcetaldehyde` + `SulfiteOxidation` + `PhenolicBrowning`,
  each drawing its own `−kᵢ·f(T)·[o2]` (SO₂'s bilinear), so total depletion is `(k_ethanol + k_browning +
  k_so2·[HSO₃⁻])·f(T)·[o2]` and O₂ is consumed **once**. To hold the calibrated total O₂-scavenging
  timescale, `k_ethanol_oxidation` is **reduced** at that build so `k_ethanol + k_browning ≈ 5.0e-4` (the
  present sole-sink value; browning then diverts O₂ from acetaldehyde exactly as SO₂ does — the *always-on*
  analogue of SO₂'s protection, and the acetaldehyde partition `k_ethanol/Σk` emerges for free). **This is
  the reduction D-73 made possible and D-71 could not express under "total rate".**
- **Where the O₂ goes (scoped honestly).** This rework closes the **O₂-accounting** seam — an always-on sink
  can now consume its `o2` share without perturbing the others. Whether browning's *product* is representable
  is a **separate** question: it produces brown melanoidin/quinone polymers, for which **no state pool
  exists**. So Browning lands cleanly today as a **pure O₂ diverter** (no product pool — it suppresses
  acetaldehyde, the visible/sensory payoff being the acetaldehyde it *prevents*), and gains a browning-index
  (`A420`) pool only if/when that readout is wanted (a D-67-style diagnostic, off the ledger). Strecker
  likewise needs new aldehyde aroma pools for its products. **D-73 permits the O₂ accounting, not the product
  side** — the product pools are their own future beats.

**Supersession discipline (the advisor's must-not).** The D-71/D-72 entries are **left as written** (true when
written); this D-73 entry *supersedes* D-71's "total rate" framing. The **live** docstrings + `aging.yaml`
provenance are updated to the share framing (they must describe the code as it is now); the decision log is
**appended, never rewritten**. The stale "consumes the WHOLE flux / remainder = unmodeled sinks" language was
removed from the `aging.py` module + class docstrings, the `derivatives` comments, and `aging.yaml` (three
regions), and the `add_oxygen` verb docstring now names both O₂ consumers.

**§4.3 firewall / tier.** Unchanged — still speculative in FORM (Tier-3 frontier). The `y` bump is a
speculative-magnitude re-estimate within the sourced band; the sourced load-bearing claims (O₂ is the
rate-limiter; warmer oxidises faster) are untouched.

**Regression surface.** **Zero structural churn** — no schema change, no new/removed Process, no `touches`
change, so `test_media.py` goldens are untouched. The only numerical move is `y_acetaldehyde_per_o2` 1.0 →
1.5; every aging test reads it from params (closed-form, saturation-ceiling) or asserts *relative* magnitudes
(oxidative > reductive; more-SO₂ ⇒ less acetaldehyde; OAV climbs), so all 782 stay green — the change is
*invisible to the suite by construction*, which is exactly why the worked drop-in above is the real
acceptance artifact. **Next:** `PhenolicBrowning` / Strecker as *always-on* O₂ sinks (now unblocked — reduce
`k_ethanol_oxidation` to its share as each lands, per the drop-in), oak extraction (a separate axis, no O₂),
or the deferred beat 1b (descriptor projection).

## D-74 — `PhenolicBrowning` built: oxidative browning is the first *always-on* O₂ sink, and it makes the browning **visible** (§4.1)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, the fourth aging Process — the first *always-on* sink the D-73
rework enabled, and the first aging Process to add a NEW observable.** `PhenolicBrowning` (medium-agnostic)
draws its share of the shared `o2` budget (`d(o2)/dt = −k_browning·f(T)·[o2]`, first-order like ethanol
oxidation) and accumulates a new state slot, `A420` — the **oxidative-browning index** (absorbance at 420 nm,
dimensionless AU) — at `d(A420)/dt = +y_a420_per_o2·(r_o2/M_O2)`. It is the **dominant** always-on O₂ consumer
(phenol autoxidation is the primary O₂ sink; ethanol oxidation is a secondary H₂O₂ fate), so it **diverts most
of the always-on O₂ away from ethanol oxidation and suppresses oxidative acetaldehyde** — the always-on
analogue of D-72's SO₂ protection, but *permanent* (a co-resident sink, not a spent one). This is exactly the
worked drop-in D-73 published as its acceptance artifact, now realised. New tests pass (browning unit +
scenario), the D-73 worked drop-in's O₂-accounting seam is exercised for the first time, `ruff`/`mypy` clean.
**One `advisor()` pass before writing** (it caught the medium-scope blocker and reframed the product fork),
and **the product fork was put to the owner** (per "surface design decisions before building"): pure O₂
diverter vs an observable browning index — the owner chose *"do what is closer to reality"*, i.e. the
observable (a Process named "browning" that produces no brown fails the correspondence bar).

**The owner's fork — build the OBSERVABLE, not a pure diverter (the load-bearing decision).** The advisor's
sharp reframe: a pure O₂ diverter with `E_a_browning == E_a_ethanol_oxidation` and no product is *algebraically
identical to just lowering `y_acetaldehyde_per_o2`* — both sinks are first-order in `[o2]` with the same
`f(T)`, so the split is a constant and cumulative acetaldehyde = `y·(k_eth/Σk)·O₂ = y_eff·O₂`. So
`PhenolicBrowning` earns *independent existence* only via a distinct `E_a` (temperature-dependent partition) OR
an observable. The owner's "closer to reality" settles it toward the observable: aged white wine's single most
visible signature is the gold→amber→brown (the A420 index), so **the observable is what makes the Process
faithful** — and it makes browning non-degenerate regardless of `E_a` (a `y`-cut produces no A420).

**`A420` is the `iso_alpha` STATE-SLOT pattern, NOT the D-67 post-hoc OAV (the advisor's precision).** Browning
pigment is **cumulative and irreversible**, and its O₂ flux is **dynamic** (SO₂ competes for O₂, temperature
varies), so `A420` must be **integrated along the run** — it cannot be reconstructed after the fact from
(dosed − remaining) O₂ the way a D-67 OAV series is computed from a finished trajectory. Two corollaries,
both load-bearing: (1) `A420` is documented as an **optical absorbance index (AU, dimensionless), NOT a pigment
mass** — which is *why* it is legitimately off every ledger (the pigment's carbon would come from an *untracked*
phenol pool; an optical index sidesteps conservation entirely). So `PhenolicBrowning` touches only `{o2, A420}`
— **both off every ledger** — and moves **nothing conserved at all**, the cleanest aging Process on the books
(cleaner even than `OxidativeAcetaldehyde`, which still borrows carbon E→acetaldehyde). (2) `d(A420)/dt ≥ 0`
always — monotonic accumulation, no clamp needed (the `o2 ≤ 0` guard also absorbs a solver undershoot).

**Medium-agnostic — D-74 SUPERSEDES D-73's provisional "wine-only" parenthetical (supersession discipline).**
D-73's worked drop-in tentatively wrote `class PhenolicBrowning(Process): # ... wine (o-diphenols are a wine
pool)`. That is superseded here, for two reasons the advisor surfaced. **Physics first:** there is *no*
o-diphenol pool (the catalyst is lumped into `k_browning`, as in `k_ethanol_oxidation`), and **both** wine and
beer carry autoxidising polyphenols that consume O₂ and brown oxidatively — so browning is a property of the
molecules, not the biology (the shared-`aging.yaml` discipline), and belongs in **both** media like ethanol
oxidation. **And it is forced to be consistent:** the `k_ethanol_oxidation` reduction (below) lives in the
**shared** `aging.yaml` and applies to both media, so a *wine-only* browning sink would leave **beer's** total
O₂-depletion rate silently **halved** below the 5.0e-4 anchor — the exact in-tree inconsistency the D-73 rework
existed to remove. The `ProcessSet` touches-contract *surfaced* this: a both-media Process cannot own a
one-media slot (`touches ⊆ schema.names` is validated at construction), so a medium-agnostic browning **forces**
`A420` into both schemas — which is also its correct architectural home (a general oxidation product, with
`esters`/`acetaldehyde`/`o2`, not the wine-only pH/SO₂ cluster).

**The `k_ethanol_oxidation` reduction — spending the D-73 seam (5.0e-4 → 2.0e-4).** D-73 redefined
`k_ethanol_oxidation` as the ethanol *share* (not the total) precisely so an always-on sibling could be added
without double-counting; D-74 spends that. Each sink owns its first-order-in-`[o2]` rate, `ProcessSet` sums
them, so the pool depletes **once** and O₂ splits by `kᵢ/Σk`. `k_browning = 3.0e-4` (the **dominant** ~60%
share) + `k_ethanol_oxidation = 2.0e-4` (the secondary ~40%) = **5.0e-4**, the calibrated total O₂-depletion
rate (the anchor), unchanged — so the O₂ *timescale* is unchanged; only the *partition* moved. `E_a_browning`
is set **equal** to `E_a_ethanol_oxidation` (50 kJ/mol): browning and ethanol oxidation are the **same
coupled-oxidation cascade**, so equal `E_a` is the honest default AND keeps the partition temperature-
**independent** (the sum is exactly `5.0e-4·f(T)` at every T, not just T_ref). A distinct/higher browning `E_a`
(maderization is arguably more T-sensitive) was considered and **rejected**: it would encode a partition-shift
direction that is not clearly sourced, and — because the A420 observable already makes browning non-degenerate
— it is not needed (its own param regardless, per prime directive #2, the `E_a_so2_oxidation` precedent).

**Acetaldehyde partitions DOWN — resolving the three D-73 "interim, until browning is built" caveats.** D-73
planted, in three live spots (the `OxidativeAcetaldehyde` docstring, the `y_acetaldehyde_per_o2` provenance,
the module docstring) plus `k_ethanol_oxidation`'s own "currently == the total (sole always-on sink)" note, the
flag that aged-with-O₂ acetaldehyde was an **upper estimate** pending browning. All four are flipped here from
"interim/pending" → "browning now takes its share." With browning as the dominant sink, the ethanol route's
share of a fully-consumed O₂ charge (no SO₂) is `k_ethanol/(k_ethanol+k_browning) = 0.4`, so cumulative aged
acetaldehyde is **~40 %** of the D-73 sole-sink value (~82 → ~33 mg/L at a 40 mg/L O₂ dose) — the "partitions
down" D-73 promised, realised — with the balance of the O₂ now going to **visible browning** (A420). Also
reconciled: `k_so2_oxidation`'s provenance note hardcoded "~8–16× the ethanol-oxidation `k_ethanol_oxidation`
(5e-4/h)" — a **live inconsistency** the moment `k_ethanol` dropped to 2e-4; corrected to compare against the
**combined** always-on total (`k_ethanol + k_browning = 5e-4`), which is what SO₂ actually out-competes.
`y_acetaldehyde_per_o2` **value is unchanged** (1.5, the route's own yield) — only its provenance updated; the
partition moved via the `k`'s, exactly as the shares design intends.

**§4.3 firewall / tier.** Unchanged — speculative in FORM (Tier-3 frontier). The browning *form* (O₂-limited,
warmer-faster) is sourced; `k_browning`, `E_a_browning` and (especially) `y_a420_per_o2` are order-of-magnitude
estimates. `y_a420_per_o2` is flagged the **most speculative parameter in `aging.yaml`** — A420-per-O₂ is not a
tabulated quantity, so it is an author estimate anchored to observed white-wine A420 ranges (fresh ~0.05–0.1;
browned ~0.3–0.6+ AU) over a plausible O₂ exposure, banded an order of magnitude; only the *monotonic visible
browning* is load-bearing, not the exact absorbance.

**Regression surface.** A new state slot (`A420`) in **both** schemas moved the `test_media.py` size goldens
(wine 40→41, beer 22→23), the `SHARED` name tuple, the units dict (A420 = "AU"), and the `AGING_PROCESSES`
set. The only kinetic-value move is `k_ethanol_oxidation` 5.0e-4 → 2.0e-4 (+ the three new params); every
existing aging test reads params or asserts *relative* magnitudes (oxidative > reductive; more-SO₂ ⇒ less
acetaldehyde), so all stay green. New tests: `PhenolicBrowning` unit (closed form, first-order-in-O₂, dominant
share, monotonic A420 + saturation, medium-agnostic on beer, reductive isolability, BOTH ledgers flat, the
headline O₂-diversion suppressing acetaldehyde to the ~40 % share, tier floor) + scenario (compile-seam gate,
A420 climbs oxidative / 0 reductive, browns the beer path, carbon closes end to end). **Next:** Strecker
degradation (the next always-on O₂ sink — reduce `k_ethanol_oxidation` again to its share, per the drop-in;
needs new aldehyde aroma pools for its products), oak extraction (a separate axis, no O₂), or beat 1b
(descriptor projection). A `PhenolicBrowning` product-pool beyond the A420 index (real melanoidin speciation)
would need a tracked phenol pool — deliberately out of scope (the optical index is the faithful v1 observable).

## D-75 — `StreckerDegradation` built: the O₂/amino-acid Strecker aldehydes (cooked-potato + honey), a WINE-ONLY *substrate-gated* sink (§4.1)

**Date:** 2026-07-10. **Milestone 3 / Tier-3, the fifth aging Process and the third oxidative sibling on
the shared `o2` budget** (after `OxidativeAcetaldehyde` D-71 and `PhenolicBrowning` D-74). `StreckerDegradation`
(WINE-ONLY) models the oxidative Strecker route: dissolved O₂ — via the o-quinones of phenol autoxidation (the
browning cascade) — oxidatively deaminates and decarboxylates amino acids to **Strecker aldehydes**,
**methional** (from methionine, the "cooked-potato" *oxidative off-note*, the marker of an oxidised/maderised
white wine and of stale beer) and **phenylacetaldehyde** (from phenylalanine, the "honey/floral" note of aged
white/dessert wines). It is the first aging Process to add **aroma pools the D-67 OAV lens did not previously
read** (two new single-molecule state slots + two new `sensory.yaml` thresholds). All new tests pass (13 unit +
6 scenario), `ruff`/`mypy` clean, the full suite green. **One `advisor()` pass before writing** (it surfaced the
`amino_acids`-availability blocker, reframed the substrate-gate justification, and prescribed the CO₂ decarb
term), and **two forks were put to the owner** (per "surface design decisions before building"): the pool
granularity and the O₂-accounting divergence.

**The verified blocker — `amino_acids` is 0 post-ferment unless dosed (must-verify, not assume).** The advisor
flagged that the whole design produces *nothing* if the `amino_acids` pool is spent when `begin_aging` fires.
Checked on a representative trajectory: with **no** amino-acid dose the pool is **0** at the aging-segment start
(`AminoAcidAssimilation` draws it down during ferment); with a modest must dose (0.5 g/L) **~11.6 mg/L** survives
into the aging segment and holds constant across it. So Strecker is **substrate-gated exactly like `mercaptans`**
(needs autolysis dosed) **and `SulfiteOxidation`** (needs SO₂ dosed): silent by default (physically correct — a
wine with no residual amino acids / no lees makes no Strecker), exercisable by dosing `amino_acids_gpl` (the
nutrient-rich / aged-on-lees case where Strecker aldehydes actually form). A future lees-autolysis refill would
make it fire on an un-dosed sur-lie aging; deferred.

**Owner fork 1 — TWO pools, not one lumped (opposite sensory valence).** Methional is an off-note (cooked
potato), phenylacetaldehyde is *pleasant* (honey) — lumping them under one threshold would be sensorially
incoherent. The owner chose **two single-molecule pools** (`methional`, `phenylacetaldehyde`) over the
esters/fusels-style single lump. Cost: two new wine slots + two thresholds + a composition split parameter
(`f_methional`, the methional mol share). Booked at each aldehyde's own carbon fraction (methional C4H8OS,
phenylacetaldehyde C8H8O), read by the D-67 lens with descriptors "cooked potato / oxidative" and "honey /
floral".

**Owner fork 2 — "closer to reality" ⇒ SUBSTRATE-GATED, add on top, NO re-baseline (supersedes the D-71→D-74
forward-guess).** D-71→D-74 repeatedly forecast Strecker as "the next *always-on* O₂ sink — reduce
`k_ethanol_oxidation` again to its share." **That guess is wrong for Strecker.** Gating the O₂ draw on
`amino_acids` (`r_o2 = k_strecker·f(T)·[o2]·[aa/(K+aa)]`) makes Strecker **doubly substrate-gated** (on `o2` AND
`amino_acids`), exactly like `SulfiteOxidation` (on `o2` AND SO₂). D-72's load-bearing rule: **a substrate-gated
sink adds on top of the shared O₂ budget without any re-baseline** — zero without its substrate ⇒ the
default/beer trajectory is byte-for-byte preserved. So `k_ethanol_oxidation + k_browning = 5.0e-4` is
**untouched**; `k_strecker` is a small wine-only draw that fires only when `amino_acids` is present. The
alternative offered (re-baseline the shared `k_ethanol` for Strecker) was flagged as re-introducing the exact
in-tree inconsistency the D-73/D-74 rework removed (it would wrongly cut *beer's* O₂ budget for a wine-only sink
that is zero by default). The owner chose "closer to reality" — i.e. the substrate-gated add-on. The
D-71→D-74 forward-notes are **retired** here (supersession discipline): the `k_ethanol_oxidation` provenance note
that said "a further always-on Strecker sink would reduce this share again" is corrected in-tree to point at this
decision. Gating the **O₂ draw itself** (not just the aldehyde) on `aa` is load-bearing (advisor): O₂/carbon/N
all vanish together as `aa` empties, so the sink reverts cleanly and never "assigns" O₂ to a product that
cannot form.

**Carbon + nitrogen close by construction — the D-45 mercaptan idiom + a decarboxylation CO₂ term.** The aldehyde
carbon is drawn from `amino_acids` (booked as arginine) and the amino-acid nitrogen **deaminated** to the `N`
pool, exactly as `AutolyticMercaptan` (D-45) does; the Strecker **decarboxylation** adds a product that idiom did
not have — **1 mol CO₂ per mol aldehyde** (the carboxyl carbon, on the carbon ledger — do **not** skip it, the
advisor's must-fix). The arginine draw is *sized to the total product carbon* (methional + phenylacetaldehyde +
CO₂), so `total_carbon` closes to machine precision (the `EsterHydrolysis` multi-product split idiom); all the
arginine N lands in `N` and the products are N-free, so `total_nitrogen` closes. Verified per-RHS (residual <
1e-18) and end-to-end over the full ferment+aging trajectory (both ledgers flat; the O₂ dose flow is carbon- AND
nitrogen-free). The **arginine-for-`amino_acids` stand-in** is exact on the ledger, approximate on provenance
(the drawn C/N is arginine's, not methionine's/phenylalanine's) — the same honest stand-in `mercaptans` carries.
Tier consequence (the D-45 note): Strecker **writes `N`** (deamination), so an enabled run drops structural
`tier_of("N")` PLAUSIBLE→SPECULATIVE. `total_mass` ({S,E,CO2}) sees the CO₂ term with no matching S/E debit, but
is never asserted on an aging run (the standing `OxidativeAcetaldehyde` scope-out).

**The inherited quinone double-count lump (documented, not fixed).** Mechanistically the O₂ is consumed at the
phenol-oxidation step (browning's draw), making the o-quinones that then do the Strecker deamination — so a
separate `k_strecker` `[o2]` draw formally double-counts that shared quinone step. But browning and
ethanol-oxidation **already** double-count it against each other (both independent `[o2]` draws for one coupled
cascade) — the additive-share v1 lump accepted at D-73. Strecker following suit is *consistent*; a two-stage
(O₂ → quinone pool → {pigment, aldehyde, acetaldehyde}) rework is deliberately **out of scope** (a larger
structural beat). **Scope:** this is the *oxidative* (quinone-driven) Strecker route only; the non-oxidative
Maillard/sugar-dicarbonyl route (sweet wines, thermal) is deferred, keeping Strecker honestly on the `o2`
sub-axis.

**Magnitudes (all speculative, Tier-3 frontier).** `k_strecker = 5.0e-5 /h` (a small add-on, ~10 % of the 5.0e-4
always-on total at full aa gate, but aa-throttled to a minor in-band perturbation at cellar residual-aa);
`E_a_strecker = 50 kJ/mol` (its own param per prime directive #2; warmer-faster, the canonical beer-staling
direction); `y_strecker_per_o2 = 0.5 mol/mol` (the quinone-mediated per-O₂ aldehyde yield, discounted for
competing quinone fates); `f_methional = 0.6` (methional-dominant mol split, an empirical composition estimate —
hence a YAML param, unlike the stoichiometric 5:2 ester split). Thresholds: `threshold_methional_wine = 0.5
µg/L` (very potent), `threshold_phenylacetaldehyde_wine = 1.0 µg/L`. Verified end-to-end: a 40–60 mg/L O₂ +
0.5 g/L amino-acid aged wine reaches **~350 µg/L methional / ~270 µg/L phenylacetaldehyde** (OAV ~700 / ~270) —
strongly threshold-positive, and the `amino_acids` pool is the hard cap on total aldehyde.

**§4.3 firewall.** Speculative in FORM (the Strecker *form* — O₂-linked, amino-acid-driven, warmer-faster,
aldehyde = amino acid − CO₂ — is sourced; the magnitudes are estimates). Writing `N`/`CO2` (plausible-tier pools)
from a speculative Process is the accepted D-27/D-45 precedent; isolable (disable the Process and the drift
vanishes).

**Regression surface.** Two new wine-only state slots (wine 41→43, beer untouched), two new chemistry species,
two new `sensory.yaml` thresholds, two new `AromaCompound`s (wine aroma set 8→10), four new `aging.yaml` params.
`test_media.py` goldens updated (wine size + `WINE_STRECKER_SLOTS` + `WINE_AGING_PROCESSES` gains
`strecker_degradation`); `k_ethanol_oxidation` **unchanged** (the whole point — no re-baseline), so every default
/ un-aged / reductive / no-amino-acid trajectory stays byte-for-byte. New tests: `StreckerDegradation` unit
(closed form, carbon + nitrogen closure per-RHS, double-substrate isolability, first-order-in-O₂ + saturating aa
gate, methional-dominant split, wine-only no-op on beer, warmer-faster, integrated saturation alongside the
co-resident always-on sinks, the speculative tier floor incl. the structural N-write) + scenario (compile-seam
gate wine-only, aldehydes climb with O₂+amino-acid dose / 0 without either substrate, carbon **and** nitrogen
close end to end, both Strecker OAVs climb through the D-67 lens). **Next:** oak extraction (a separate aging
axis, no O₂ — diffusion-limited vanillin / whiskey lactones / gallotannins), the deferred beat 1b (descriptor
projection), a lees-autolysis `amino_acids` refill (would make Strecker fire on un-dosed sur-lie aging), or the
deferred non-oxidative Maillard Strecker route.
