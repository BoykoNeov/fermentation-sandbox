---
name: feedback-run-tests-at-below-normal-priority
description: "Always run the test suite through tools/nicepytest.py (BelowNormal priority), never bare `uv run pytest`"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c5d2411f-695e-478f-bbae-e0b90d259b55
  modified: 2026-08-28T09:21:11.787Z
---

Run the Fermentation suite as `uv run python tools/nicepytest.py -n auto` — **never** bare
`uv run pytest -n auto`. Applies to every run, not just when the box looks busy.

**Why:** `-n auto` takes every logical CPU at normal priority, which makes the desktop and any
other suite on the machine crawl. Other agent sessions run their own suites here. `nicepytest.py`
drops to `BELOW_NORMAL_PRIORITY_CLASS` before importing pytest and forwards every argument
verbatim; the class is inherited, so all N workers run at base priority 6. Priority makes the
suite *yield* a core it already holds, which a smaller `-n` never does.

**How to apply:** substitute `python tools/nicepytest.py` wherever `pytest` would go — the flags
are identical (`-n auto`, `--lf`, a single file path, `-m benchmark`). Asked for on 2026-08-28
after a bare `uv run pytest -n auto` full run. See `CLAUDE.md`'s Commands section, which already
documents the tool. Related: [[feedback-full-suite-before-green]],
[[feedback-enumerate-competitors-before-timing]].
