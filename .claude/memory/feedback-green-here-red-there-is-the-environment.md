---
name: feedback-green-here-red-there-is-the-environment
description: The same checker green locally and red in CI is a statement about the environment, and it reports consequences rather than the cause
metadata:
  type: feedback
---

A strict type checker run against a **partial** environment does not report "I am missing a
package". It reports the *consequences* of missing it, in the vocabulary of the code under test, at
every affected site. Here `streamlit` was absent on CI, so `st.stop()` lost its `NoReturn`
annotation — the only thing telling the checker that code below a guard is unreachable — and 44
narrowings collapsed at once into "Item None of X has no attribute", all of which read as real
defects in correct code. Nothing in that output named the cause.

**Why:** the count is what misleads. One missing package produced 44 errors across one file, which
looks like a broken module rather than a broken environment, and invites 44 local repairs to code
that is fine.

**How to apply:** when *the same command* is green on the dev box and red in CI, that difference is
about the environment first and the code second — diff the installed packages before reading the
errors. Reproduce by pointing `UV_PROJECT_ENVIRONMENT` at a throwaway directory rather than
mutating the working `.venv`. And check what the checker is pointed at: mypy runs on a `files` list,
so it can be checking a package whose dependencies the install step never installs.
See [[feedback-a-gate-not-pointed-at-the-code-hides-a-crash]] and [[feedback-a-green-matrix-may-vary-nothing]] (D-265).
