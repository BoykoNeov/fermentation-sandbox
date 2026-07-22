---
name: feedback-never-pipe-checks-to-tail
description: "Never verify ruff/mypy/pytest by piping to tail in an && chain — the pipe returns tail's exit code and hides failures"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 4c5542b6-994d-42ff-9b5f-a6dbc7d14d50
  modified: 2026-07-20T09:06:50.065Z
---

When running the project's checks, **never** write
`uv run ruff check . 2>&1 | tail -5 && uv run mypy ... && uv run pytest -q 2>&1 | tail -12`.
The pipeline's exit status is **`tail`'s**, which is always 0, so the `&&` chain
reports success while ruff and pytest are failing. Capture to files and echo `$?`
per command instead:

```bash
uv run ruff check . > ruff.txt 2>&1; echo "ruff exit=$?"
uv run mypy       > mypy.txt 2>&1; echo "mypy exit=$?"
uv run pytest -q  > pytest.txt 2>&1; echo "pytest exit=$?"
```

**Why:** this exact pattern reported "exit code 0" at D-117 while 4 ruff errors
and 2 failing conservation tests sat in the output — a red tree that would have
been committed as green. The Fermentation suite runs ~11 minutes, which is
precisely what makes truncating the output tempting.

**How to apply:** the full suite is slow, so run it backgrounded and read the
captured files; judge pass/fail on the echoed exit codes, never on the tail of
the output. Related: [[feedback-always-commit-push]] — each commit is supposed to
pass all three, so a masked failure defeats the rule.
