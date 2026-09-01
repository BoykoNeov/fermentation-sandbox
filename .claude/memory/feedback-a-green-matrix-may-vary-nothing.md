---
name: feedback-a-green-matrix-may-vary-nothing
description: A CI matrix claims it varies something; nothing in a green run verifies that it did
metadata:
  type: feedback
---

A CI matrix **asserts** that it varies something. A green run does not verify the assertion — a
harness that silently does less than it claims still comes back green. The two-version matrix here
ran one version for the entire life of the matrix: `uv sync --python 3.13` honoured the pin, then
every bare `uv run` after it took the newest interpreter it could find, **deleted the environment
and rebuilt it** as 3.14 with the default dependency groups. No test had ever executed on 3.13.

**Why:** the faults a harness reports are faults in the code. A fault in the harness reports
nothing, so it is invisible in exactly the state everyone treats as proof. Eleven red runs were
needed to expose this; a hundred green ones would have concealed it. It surfaced only because
something unrelated forced the two jobs to disagree.

**How to apply:** read the harness's own log lines, not just its verdict — the evidence here was
"Removed virtual environment" inside a step whose output was "All checks passed!". When a matrix
axis matters, make each step name the axis rather than trusting inheritance, and check once that
the axis is really varying. Treat flags that pin an axis as load-bearing, never boilerplate:
deleting them restores the bug silently, and what comes back does not look like a harness fault.
See [[feedback-a-fast-gate-masks-the-slow-one-behind-it]] and [[feedback-measure-the-surface-in-the-unit-that-fails-it]] (D-265).
