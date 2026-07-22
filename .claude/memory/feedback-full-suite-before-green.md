---
name: feedback-full-suite-before-green
description: "A new Process in a shared registry needs the FULL pytest suite green before claiming \"all green\" — the domain suite misses exact-set and end-to-end tests"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d345979b-b042-43e6-b68f-31e19f6eb0ef
  modified: 2026-07-21T09:40:46.108Z
---

When a change adds a `Process` to a shared registry (`_AGING_PROCESSES` /
`_AGING_GATED_PROCESSES` in `media.py`/`compile.py`), running only the domain
suite (`test_aging.py`) is NOT enough to declare it green. Two classes of test
live elsewhere and WILL break:

- **Exact-set assertions** — `test_media.py::test_registered_media_wire_the_full_kinetic_set`
  compares `{p.name for p in pset.active}` against a hardcoded `EXPECTED_PROCESSES`
  set (+ the `AGING_PROCESSES` set literal). A new process is an "Extra item".
  Update the expected set.
- **End-to-end interaction tests** — `test_aging_scenario.py` runs the FULL aging
  set in a compiled scenario, so a new aging Process that touches a shared pool
  (e.g. `Byp`) can flip an A/B assertion another Process owned. At D-127
  `EthylAcetateEsterification` consuming `Byp` (forming EtOAc) flipped
  `byp_aged > byp_plain`.

**Why:** At D-127 I ran the aging suite (229 green) + mypy + ruff, reported "all
green / regression risk low", and pushed — then the full suite caught 3 failures.
The commit itself was fine (see [[feedback-always-commit-push]]), but the CLAIM
was premature.

**How to apply:** For a shared-registry/Process change, don't say "all green"
until `uv run pytest -q` (the WHOLE suite, ~15 min) passes. The suite buffers
output to the end, so a 0-byte task-output file means still running, not passing.
Commit per the always-commit rule, but phrase the status honestly ("domain suite
+ mypy + ruff green; full suite still confirming") until the full run returns.
