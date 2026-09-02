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

Five layers, strictly one-directional dependencies (lower layers know nothing of
higher ones):

```
  app (interface)         local console + written report — NOT part of the package
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

Milestones 0 (skeleton), 1 (single-strain primary fermentation) and 2 (pH as a charge
balance, byproducts, MLF/Brett, the amino-acid ledger) are complete; Milestone 3 (the
oxidative aging axis and the sensory readout) is where the work is now. The engine carries
around 80 Process implementations across two media, a stochastic ensemble runtime, a local
console, and a test suite of roughly 2100 tests. Nothing is tiered `validated`: that tier is
reserved for checks against independent measured time-series, and the §2.2 acceptance
benchmarks earn `plausible`.

Three documents carry the state, and none of them is written by hand from memory:

- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) — the structure map: where every subsystem
  lives, the state vector, the Process registry, with every count derived from the code.
- [`docs/DECISIONS.md`](docs/DECISIONS.md) — the append-only archive of every engineering
  decision, with a generated index and correction map.
- [`docs/OPEN.md`](docs/OPEN.md) — what is scientifically open right now, generated from the
  suite's `xfail` markers and the archive's declared-but-unshipped reversals.

`docs/plans/` holds the original milestone plans as frozen logs; they are history, not status.

## Running it

The console is the way in: a local page where you describe a batch, run it, and follow any
line on the chart back to the papers the numbers came from. It runs entirely on your own
machine and uploads nothing.

**Windows** — double-click `start-console.cmd`.
**macOS / Linux** — run `./start-console.sh`.

Either one checks what is needed, installs the interface dependencies on first use, picks a
port that is free, and opens the page in your browser. The first start takes a minute or two
while the pieces are downloaded; after that it is a few seconds. Leave the terminal window
open while you use it — closing it shuts the console down.

The only thing to install first is [uv](https://docs.astral.sh/uv/getting-started/installation/),
which fetches the right Python (≥ 3.13) and everything else by itself. The launcher says so
if it is missing; it will not install anything for you.

The equivalent typed by hand, if you would rather:

```bash
uv sync --group ui
uv run streamlit run app/main.py
```

The console is optional and separate: a plain `uv sync` still installs a four-dependency
research library with no interface code in it. See [`app/README.md`](app/README.md) for what
the page shows and why it is arranged that way.

## Working on the code

```bash
uv sync                 # create the venv and install deps (incl. dev tools)
uv run pytest -n auto   # the full test suite, in parallel
uv run ruff check .     # lint
uv run ruff format .    # format
uv run mypy             # type-check
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
