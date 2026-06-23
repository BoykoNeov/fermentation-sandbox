# Milestone 2 — context (key files & decisions)

## Opening decisions (DECISIONS D-18)

- **pH = full charge-balance solver**, derived-algebraic (no `dpH/dt`). Resolves the
  handoff §7 open decision; chosen over tracked-pH because tracked-pH can only do MLF/
  SO₂ by *scripting* (violates compositionality, handoff §5).
- **Byproducts/temperature beat is built first**, ahead of pH (inverts handoff §6).

## Where Tier-2 extends the M1 core (what already exists)

- **State / schema:** `src/fermentation/core/state.py` — `StateSchema`, `VarSpec`.
  New byproduct pools reuse the `VarSpec.default = 0.0` **produced-only pool** pattern
  (D-16) that `Gly`/`Byp`/`X_dead` use, so they land without touching pack call sites.
- **Media / schemas:** `src/fermentation/core/media.py` — `wine_schema`/`beer_schema`,
  `Medium`, `MEDIA`. Add `esters`/`fusels` to both schemas; wire the new Processes via
  each `Medium`'s `process_factories` (additive — no `modifier_factories` needed).
- **Processes:** `src/fermentation/core/process.py` — subclass `Process`
  (`name`/`tier`/`touches`/`reads` + `derivatives`). Byproduct Processes are **additive**
  (they produce a compound), so unlike M1's inhibition/Arrhenius they need **no**
  `RateModifier`.
- **Temperature axis:** `src/fermentation/core/kinetics/arrhenius.py` — the
  `ArrheniusTemperature` modifiers already scale growth & uptake; **inert at `T_ref`**
  (D-11/D-17), so a sub-`T_ref` run is what finally exercises them. Byproduct rates
  reuse the same `f(T)` shape or their own sourced T-sensitivity.
- **Stoichiometry / carbon:** `src/fermentation/core/chemistry.py` +
  `src/fermentation/validation/conservation.py` — if byproduct pools are carbon-routed
  (recommended, D-16 style), extend `total_carbon` to weight them, drawing carbon from
  the same shared chemistry constants the kinetics use (single source of truth).
- **Runtime:** `src/fermentation/runtime/integrate.py` — `simulate`/`Trajectory`. The
  **stochastic ensemble wrapper** layers here (over `simulate`), outside the pure core.
- **Scenario / compile seam:** `src/fermentation/scenario/compile.py` —
  `compile_scenario`. Scenario-specific *values* evaluated at the boundary already have
  a precedent (D-14 N-dependent yield, D-16 fermentable fraction); the pH beat's
  per-acid initial concentrations will key in here too.
- **Benchmarks:** `tests/benchmarks/test_milestone1.py` — the still-skipped
  `test_lower_temperature_is_slower_but_cleaner` is this beat's gate.
- **Param store:** `src/fermentation/parameters/data/*.yaml` — `wine_generic.yaml`,
  `beer_generic.yaml`. `Uncertainty` ranges on every entry already feed the stochastic
  wrapper.

## Parameters to source (each needs full provenance, handoff §2.3)

**Byproducts beat:** ester-synthesis rate + T-sensitivity; fusel (Ehrlich) rate + its
non-monotonic N dependence. Reconcile against de Andrés-Toro et al. 1998 (dynamic beer
model with byproduct terms), Malherbe et al. 2004 (wine), and ester/higher-alcohol
kinetic literature. Tag `plausible` only where a source measures *our* form, else
`speculative` (directional validation only — handoff §3.5).

**pH beat (later):** pKa set for tartaric/malic/lactic/acetic (± carbonic) — textbook
but T-dependent; per-acid initial concentrations become scenario inputs (must/wort
composition). SO₂ pKa ≈ 1.81. MLF: *O. oeni* growth/decay constants, malic→lactic yield.

## Open questions to resolve at beat start

- **Byproduct carbon source** — route from sugar (D-16 style, exact carbon closure,
  recommended) vs trace-unaccounted. Fusels are Ehrlich (amino-acid-derived), so their
  carbon skeleton is arguably N-pool-linked, not sugar-linked — pick one cleanly.
- **Fusel vs N non-monotonicity** — first pass may approximate; the handoff insists on
  modelling the *pathway*, not fitting a slope. Decide how much pathway detail v1 carries.
- **pH beat (deferred to its own start):** the three D-18 couplings — evolved-vs-dissolved
  CO₂ for carbonic, acid carbon vs the `Byp`=succinic lump, and pKa(T).
- **Stochastic wrapper API** — sampling strategy (LHS vs plain MC), how tier/uncertainty
  bands are reported on ensemble outputs. Engineering choice, not a scoping gate.
