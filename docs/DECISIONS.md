# Design decisions

Lightweight decision log. Each entry: the decision, the rationale, and (where
relevant) how it deviates from the handoff brief. The handoff explicitly states
"nothing is rock solid"; this file records where we reasoned past it.

## Process decisions (project setup)

These three were the handoff's §7 open questions, resolved with the project owner.

- **D-A — Repository:** public GitHub repo `fermentation-sandbox` under
  `BoykoNeov`, MIT licensed.
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

**Status: open (scoping decision).** This opens Milestone 2 (Tier-2). It records two
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

## Deferred (decide early in the relevant milestone)

- ~~**pH / acid model richness**~~ — **decided in D-18** (full charge-balance solver),
  built after the byproducts beat.
- ~~**Stochastic ensemble API**~~ — **scoped into Milestone 2** (`milestone-2-plan.md`);
  physics-free runtime wrapper over the existing `Uncertainty` ranges, buildable in
  parallel with the byproducts beat.
- **Packaged parameter-data access:** tests read YAML via filesystem path. If we
  ship a wheel that must read its own data, switch to `importlib.resources`.
