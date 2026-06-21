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
**Status (M0):** `ProcessSet.tier_of` currently propagates **Process** tiers
only. A Process must also be capped by the tiers of the *parameters* it consumes
(a VALIDATED process on speculative params must report speculative, not
validated). That propagation lands in Milestone 1, when real Processes declare
the parameters they read — see `milestone-1-tasks.md`. Until then the code's
guarantee is narrower than this entry's "Processes and parameters" intent.

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
they are implemented. `beer_generic.yaml` does not exist yet, so beer compiles
only with an explicit `parameter_paths=` override (a clear `FileNotFoundError`
otherwise).

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
  same reason.
- **Biomass elemental composition** (C-fraction ≈ 0.48, N-fraction ≈ 0.11 from the
  canonical `CH₁.₈O₀.₅N₀.₂` formula) is *empirical and uncertain* and is consumed
  by both the conservation check and the growth Process — so it is a **Parameter**
  (provenance store), not a code constant. `total_carbon`/`total_nitrogen` take the
  biomass fraction as a **passed-in argument** (the caller resolves it from the
  store) rather than importing the loader into the core/validation math; if a
  schema has an `X` variable and no fraction is supplied, the builder raises rather
  than silently under-counting (which would report a *false* violation).

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
