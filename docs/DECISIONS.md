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
the `a1` **mean** as `−1.08×10⁻³`, but that value lies **outside its own printed 95 %
credible region** `[−1.94×10⁻¹, −3.30×10⁻²]` — a posterior mean cannot. The journal
typesetting dropped the `×10ⁿ` exponent from the `a1` mean column; the true value is
`−1.08×10⁻¹`. Three independent checks confirm it: (1) it reproduces the paper's stated
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
day benchmark window** — an *uptake-speed* gap, since the Table A2 reconstruction runs
~2× faster than Coleman's measured 20 °C data even with *zero* death; this is the next
tuning question (β_max/biomass), not a `k'_d` matter. (b) ABV lands at **16.9 %** (E ≈ 133
g/L) from the theoretical Gay-Lussac split — the realised-yield/glycerol-sink task.
Neither was folded into the `k'_d` decision.

## Deferred (decide early in the relevant milestone)

- **pH / acid model richness** (handoff §3.4, §7): full proton/charge balance vs.
  a tracked-pH approximation upgraded later. Decide at the **start of Tier-2** —
  pH feeds SO₂ speciation and microbial growth, so it unblocks much of Tier-2.
  Not needed for Milestone 1.
- **Stochastic ensemble API** (handoff §1.6): parameter sampling within
  provenance bounds as a runtime wrapper. The `Uncertainty` ranges already exist
  to feed this; design the wrapper during "runtime maturation" (handoff §6.3).
- **Packaged parameter-data access:** tests read YAML via filesystem path. If we
  ship a wheel that must read its own data, switch to `importlib.resources`.
