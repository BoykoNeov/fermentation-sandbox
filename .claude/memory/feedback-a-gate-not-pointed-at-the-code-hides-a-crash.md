---
name: feedback-a-gate-not-pointed-at-the-code-hides-a-crash
description: "New code outside the type checker's `files` is unchecked in a way lint does not cover — pointing mypy at it found an AttributeError on a branch that had never executed. Extend the gate, don't document the gap"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 86025c03-b1ae-45b8-b84a-18332d969846
  modified: 2026-09-01T12:15:47.255Z
---

**A quality gate covers a path list, not "the repo", and new code lands outside it silently.**
`ruff check .` walks the tree, so a new top-level package is linted the day it appears; `mypy`
runs on `files = [...]` in `pyproject.toml`, so the same package is *invisible* to it and nothing
says so. The gap reads as "we chose not to type this" long after nobody chose anything.

D-261. `app/` was outside `files`, and the obvious move — reviewed and recommended as such — was
to state the exemption in the commit and move on. Pointing mypy at it instead cost one empty
`src/fermentation/py.typed` (so the engine's own types are visible from outside the package,
which a library this strict should ship regardless) and one `ignore_missing_imports` override for
the three UI dependencies, exactly as `scipy` already had. It found **two real defects in
minutes**, one of them `ens.ensemble.sampled_names` on a bare `Ensemble` — an `AttributeError`
on the branch that warns when the uncertainty ranking is underdetermined, i.e. **a line that had
never executed and that no amount of clicking would have reached before a user did**.

**Why:** tests only reach branches someone thought to enter; a type checker reaches every branch
whether or not it runs. That is precisely its value on code with rarely-taken warning paths,
error paths and optional-dependency paths — which is most of an interface layer. And the cost of
adding it is usually one marker file and one override, not the stub-writing project it looks like
from outside.

**How to apply:** when adding a package, check every gate's *scope* — `pyproject.toml`'s mypy
`files`, ruff's `src`, pytest's `testpaths`, CI's own invocation — and extend each before writing
much code. `ruff check app/` passing is not `ruff check .` passing: run the repo-wide form, since
per-directory runs miss the tests you just wrote. If a dependency ships no types, add an
`ignore_missing_imports` override beside the existing ones rather than exempting your own code.
Prefer extending the gate over documenting the gap: a documented gap is a permanent one.
Related: [[feedback-full-suite-before-green]], [[feedback-never-pipe-checks-to-tail]].
