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
uv run pytest -n auto   # FULL suite in parallel — ~1.6 min vs ~11.5 min serial (must stay green)
uv run pytest tests/test_<module>.py   # iteration loop: just the file you're editing (~2-5 s)
uv run pytest -n auto --lf             # after a red run: re-run only what failed
uv run pytest -m benchmark   # §2.2 acceptance benchmarks (skipped until kinetics)
uv run ruff check .     # lint
uv run ruff format --check .   # SEPARATE CI gate -- `ruff check` passing does not imply it
uv run mypy             # types (strict on src; tests exempt from signature reqs)

uv run python tools/nicepytest.py -n auto   # same, at BelowNormal — when sharing the box
```

The suite is ~1250 independent `solve_ivp` integrations, so it is embarrassingly
parallel: `pytest-xdist`'s `-n auto` gives a ~7× wall-clock win on a many-core box.
`tests/conftest.py` pins BLAS/OpenMP to one thread per worker (before numpy imports) —
without that pin, N workers each spawn N BLAS threads and the parallel run is *slower*
than pinned. Plain `uv run pytest -q` still works (serial, ~11.5 min); prefer `-n auto`.

**Sharing the machine.** `-n auto` takes every logical CPU at normal priority, which
makes the desktop — and any other suite running concurrently — crawl. Measured on the
8-core/16-thread dev box: **119 s under light background load vs 363 s** against a
competing 26-worker suite from another project. Neither figure is a quiet-box
measurement — other agent sessions run their own suites on this machine, so enumerate
competing processes before trusting any wall-clock number here.

`tools/nicepytest.py` drops to `BELOW_NORMAL_PRIORITY_CLASS` before importing pytest and
forwards every argument verbatim; the class is inherited, so all N workers run at base
priority 6 (verified by walking the process tree). Priority and worker count are
complementary, not alternatives — priority makes the suite *yield* a core it already
holds, a smaller `-n` never takes it. Under a CPU-bound competitor only the latter
preserves that competitor's throughput, so a reduced `-n` is the lever to reach for when
sharing the box matters more than finishing fast. **Which `-n` is untested**: the sweep
needs a quiet box, so `-n auto` remains the only measured configuration.

**Do not switch to `--dist worksteal`.** The default `--dist load` hands each worker a
contiguous chunk (~46 tests at 16 workers), which keeps a whole test file together and
lets expensive module-scoped fixtures be paid once.
`tests/benchmarks/test_validation_danilewicz_so2_o2.py` is the case that matters: 12 of
its 14 tests share one `runs` fixture costing ~50 s, and the durations profile confirms
it is instantiated exactly once. Worksteal re-queues tests across workers and would
re-pay that fixture per stealing worker — it buys balance and sells it back at a loss.
`--dist loadfile`/`loadscope` have the *opposite* failure mode: they pin a module to one
worker, which guarantees the fixture reuse above but serialises the largest files
(`test_aging.py` alone accounts for 231 of the recorded setups). They are a balance
regression, not a fixture one — don't cite the worksteal reason against them.

## Repo etiquette

- Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `chore:`).
- Each commit should pass `ruff`, `mypy`, and `pytest`. CI runs all three.
- **Do not weaken or delete the skipped benchmark tests** in
  `tests/benchmarks/` to make CI green — implement the model until they pass.
- **`docs/ARCHITECTURE.md` is the structure map — update it in the SAME commit** as any
  change to the layer/package layout, the Process registry, the state vector, or the
  parameter-file set. It is the answer to "where does X live"; `DECISIONS.md` remains the
  answer to "why". Its counts are *derived* — re-run the snippet at the foot of that file
  instead of hand-editing a number, and if the two disagree the code is right.
- **`docs/plans/milestone-*.md` are FROZEN LOGS, not live status.** Never read them for what
  is built, open, or next — that is `ARCHITECTURE.md` and `DECISIONS.md`. Do not resume
  maintaining them. The old "keep the plans updated" rule is retired: it was followed for a
  while, then silently stopped, and the half-updated result caused two re-proposals of
  shipped work. A doc nobody updates is worse than one that says it is history.
- The original brief is `docs/FERMENTATION_SIM_HANDOFF.md` (reference, not gospel).

## Navigating and appending to DECISIONS.md

The archive is very large (tens of thousands of lines), so **never read it start to finish**. Its top block
is generated by `tools/gen_decisions_toc.py` and gives three ways in: a
**subsystem cut** (grep a bucket name — records appear under every bucket they
match, so it is a search aid, not a partition), the **ordered record list** with
anchors, and a **correction map**. Then `Grep '^## D-<n> — '` to jump to the record.

**Check the index row before trusting a record you jumped to.** A record carries
markers for what *it* corrects, but nothing for corrections made *against* it —
those live only as ⚠ in the index. So a `Grep` straight to D-56 lands on it with
no sign that D-57 overturned it, which is precisely the trap the map exists for.

When appending a record:

- **Keep the `## D-<n> — …` heading on ONE physical line.** A wrapped ATX heading
  renders as a heading plus a stray paragraph and its anchor becomes a dead link;
  D-133…D-136 all shipped that way. `--check` now refuses this.
- **If the record corrects an earlier one, say so in a marker line** right under
  the heading — this is the only thing that grows the correction map:

  ```
  **Corrects:** D-56 — what was wrong, in a clause.
  **Flags:** D-71, D-74 — a reversal that is agreed but NOT yet in the code.
  ```

  `Corrects` = the fix shipped; `Flags` = identified, code still carries the old
  behaviour. The generator derives the back-edge, so the corrected record grows a
  ⚠ without being edited (the archive stays append-only). **The map only knows
  what a marker declares** — absence of ⚠ is not a guarantee a record is current.
- Then run `uv run python tools/gen_decisions_toc.py`. CI and
  `tests/test_decisions_index.py` both fail on a stale index, a wrapped heading,
  and a marker pointing at a nonexistent or later record.
