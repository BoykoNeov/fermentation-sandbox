# Fermentation Sandbox

A research-grade simulation engine for **wine and beer fermentation**, grounded in
published science where possible and clearly labelled as speculative where not.

This is a **research sandbox, not a game or a homebrew calculator.** The bar is
correspondence with reality, not fun or convenience. Three principles are enforced
in code and tests, not just honoured in spirit:

1. **Fidelity is tiered and explicit.** Every modelled quantity is `validated`,
   `plausible`, or `speculative`. The tier travels with the value all the way to
   any output. The engine never silently blends a validated concentration with a
   speculative one and presents them as equally trustworthy.
2. **Parameters are data with provenance, never magic numbers in code.** Every
   kinetic constant carries its source, the conditions it was measured under,
   units, an uncertainty range, and a tier — enforced at load time.
3. **The validated core is built first and protected.** Speculative layers are
   isolated so they cannot contaminate the core's numerics or its tests.

The growth path is **validated → plausible → speculative**, and the architecture
makes each expansion an *addition* rather than a rewrite.

## Architecture

Four layers, strictly one-directional dependencies (lower layers know nothing of
higher ones):

```
  scenario / validation   declarative recipes, benchmark comparison, analysis
  runtime                 time-stepping, events, phase switching, ensembles
  domain core             pure deterministic state + Process derivatives
  parameters / units      versioned data (value + provenance + tier), conversions
```

The **domain core has no UI, no file I/O, no global state, and no randomness.**
Given a state and a parameter set it returns derivatives — which is exactly what
makes it testable against benchmark curves and conservation laws.

A **Process** is anything that contributes to `d(state)/dt` (primary fermentation,
malolactic fermentation, oxidation, oak extraction). The total derivative is the
sum of the active Processes, so a speculative Process can be toggled off and the
validated core still runs and still passes its tests.

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the full design and
[`docs/DECISIONS.md`](docs/DECISIONS.md) for the design decisions and their
rationale.

## Status

**Milestone 0 — honest skeleton (current).** Layering, the Process interface, the
provenance-backed parameter store, the unit-conversion boundary, a runnable
deterministic integrator, and a validation harness (conservation checks +
benchmark specs) — all tested, *before* any real kinetics. The §2.2 acceptance
benchmarks are encoded as skipped tests waiting for the kinetics.

**Milestone 1 (next).** Single-strain, isothermal, nitrogen-limited primary
fermentation; pass the wine (~24 °Brix → dry in 10–14 days) and beer (~1.048 OG
→ ~1.010 in 5–7 days) benchmarks. See [`docs/plans/`](docs/plans/).

## Getting started

Requires Python ≥ 3.13 and [uv](https://docs.astral.sh/uv/).

```bash
uv sync                  # create the venv and install deps (incl. dev tools)
uv run pytest            # run the test suite
uv run ruff check .      # lint
uv run ruff format .     # format
uv run mypy              # type-check
```

Benchmarks (skipped until kinetics land) are marked; run only them with:

```bash
uv run pytest -m benchmark
```

## License

Boyko Non-Commercial License v1.0 (BNCL-1.0) — see [`LICENSE`](LICENSE) and
[`NOTICE`](NOTICE). Free to use, modify, and share for **non-commercial
purposes** with attribution; **commercial use is prohibited** unless separately
licensed by the copyright holder.
