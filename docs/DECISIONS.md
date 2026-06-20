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
