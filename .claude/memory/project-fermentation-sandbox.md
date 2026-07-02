---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation
engine in Python (uv, scipy/numpy/pydantic). Repo:
https://github.com/BoykoNeov/fermentation-sandbox (default branch `main`).

**This file is session-boot context, not a changelog.** The full per-decision
narrative lives in the repo and is the single source of truth — do NOT copy it
here (see [[feedback-batch-end-ritual]]):

- `docs/DECISIONS.md` — every engineering decision D-1 … D-37, with the fork,
  the choice, and the reasoning. **The canonical archive.**
- `docs/ARCHITECTURE.md` — layering, package map, the core/runtime/scenario seams.
- `docs/plans/milestone-*.md` — task checklists per milestone.
- `CLAUDE.md` — prime directives (tiers, provenance, one-directional deps).

**Status (2026-07-02):**
- **Milestone 0 & 1 COMPLETE.** M1 = single-strain isothermal nitrogen-limited
  primary fermentation passing both §2.2 benchmarks (wine ~24°Brix→dry 8–14d,
  beer ~1.048→~1.010 in ~5.5d). Tier sweep (D-17): VALIDATED reserved for
  independent measured data; the benchmark pass earns PLAUSIBLE.
- **Milestone 2 in progress** — currently at **D-37** (481 tests green + 5
  benchmark; ruff + mypy strict clean). Built so far: pH charge-balance solver
  (D-18), aroma byproducts esters/fusels (D-19–21), SO₂ speciation + free/bound
  split (D-22/28), MLF conversion + O.oeni diacetyl (D-23/31), stochastic
  ensemble + spread attribution + LHS/Sobol (D-24/25/37), diacetyl/acetaldehyde/
  H₂S (D-26/27/29), residual-N carrying cap (D-30), amino-acid ledger + fusel
  re-route + autolysis (D-32/33/34), event loop + temperature ramp + discrete
  intervention verbs `add_dap`/`add_so2`/`rack`/`pitch_mlf` (D-35/36).

**Open / next candidates** (details in DECISIONS "Deferred" + latest D-record):
- MLF-with-growth consumer Process (bacterial growth on `X_mlf` funded from the
  aa-ledger + autolysis) composed with a `pitch_mlf` event — now unblocked.
- H₂S CO₂-stripping volatilization sink (D-29 follow-up).
- Default-on residual-N model (needs a nitrogen-model redesign, not a cap; D-30).
- Mixed cultures / Brett / sour (resource competition, extended Monod).
