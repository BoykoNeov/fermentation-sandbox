# Fermentation Sandbox — project guide

A research-grade fermentation simulation engine. The bar is **correspondence
with reality**, not fun or convenience. Read `docs/ARCHITECTURE.md` for the
design and `docs/DECISIONS.md` for *why* it is shaped this way before making
structural changes.

## Prime directives (enforced in code + tests, not just honoured)

1. **Fidelity is tiered.** Every quantity is `validated` / `plausible` /
   `speculative` (`fermentation.core.tiers.Tier`). The tier must travel with the
   value to every output. Never blend tiers silently — an output's tier is the
   *lowest* of its inputs (`Tier.combine`, `ProcessSet.tier_of`).
2. **Parameters are data with provenance, never magic numbers in code.** Every
   constant lives in a YAML file under `src/fermentation/parameters/data/` and
   loads through the `Parameter` schema (value, units, source, conditions,
   uncertainty, tier — all required). A guess is tagged `speculative` with
   `source: "author estimate"`. If you need a number, add a provenance entry;
   do not inline it.
3. **The validated core is built first and protected.** Speculative Processes
   must stay isolable — togglable off without breaking the core or its tests.

## Architecture rule (one-directional dependencies)

```
scenario / validation  →  runtime  →  core  →  parameters / units
```

Lower layers must not import higher ones. The **core is pure**: no I/O, no global
state, no randomness; given state + params it returns derivatives. Randomness and
ensembles live in `runtime` as a wrapper. Physics never lives in `scenario`.

## Conventions

- **State is a plain `float64` numpy array** driven by `solve_ivp`; the
  name→index map is `StateSchema`. Tier/uncertainty do NOT ride inside these
  floats (see DECISIONS). Sugar `S` is always a vector (1 slot for wine, 3 for
  beer).
- **Canonical internal units:** concentration g/L, temperature K, time hours.
  Convert only at I/O edges via `fermentation.units` (Brix/SG/Plato/ABV/°C/days).
- **A `Process`** declares `name`, `tier`, `touches`, and returns its
  contribution to `d(state)/dt`. `ProcessSet` sums active processes and derives
  output tiers. Run `ProcessSet(..., strict=True)` in tests to enforce the
  `touches` contract.
- **Conservation laws are tests.** Carbon/nitrogen/mass must balance to
  tolerance (`fermentation.validation.assert_conserved`). A model that creates
  mass is broken regardless of how good its curves look.

## Commands

```bash
uv sync                 # install deps + dev tools
uv run pytest -q        # tests (must stay green)
uv run pytest -m benchmark   # §2.2 acceptance benchmarks (skipped until kinetics)
uv run ruff check .     # lint    (also: ruff format .)
uv run mypy             # types (strict on src; tests exempt from signature reqs)
```

## Repo etiquette

- Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `chore:`).
- Each commit should pass `ruff`, `mypy`, and `pytest`. CI runs all three.
- **Do not weaken or delete the skipped benchmark tests** in
  `tests/benchmarks/` to make CI green — implement the model until they pass.
- Keep `docs/plans/milestone-*.md` updated as work progresses.
- The original brief is `docs/FERMENTATION_SIM_HANDOFF.md` (reference, not gospel).
