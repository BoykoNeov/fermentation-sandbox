---
name: feedback-a-fast-gate-masks-the-slow-one-behind-it
description: "A pipeline that runs a cheap check before an expensive one turns the expensive one's result into an unknown while the cheap one is red — and 'still failing' reads as one problem across many commits"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ea3f3770-ac9d-4299-8476-1248aca8b862
  modified: 2026-08-27T09:26:06.371Z
---

When CI runs a fast gate (lint, format, types) ahead of a slow one (the suite), a red fast gate
does not mean the slow gate would pass — it means **nobody knows**. Fix the fast one and the slow
one can surface a second, unrelated failure that has been sitting there the whole time. Several
commits of "CI is still red" read as one continuing problem and can be two.

**Why:** the Fermentation repo's CI was red from D-239 through D-241. The cause was
`ruff format --check`, a gate the project guide never listed separately from `ruff check` — four
unformatted files, three commits inheriting them. Repairing that revealed a **second** cause the
format gate had been hiding for the same span: a D-238 guard asserting that a root-finder returns
a shipped literal bit for bit, true on the author's Windows box and 1–3 ULP off on both CI
Pythons. It had never passed on CI. Nobody noticed because pytest had not run since D-238.

**How to apply:** **read the duration, not just the colour.** In this repo a ~17 s failure is a
gate and a ~15 min failure is the suite, so the two are trivially distinguishable in
`gh run list` — and their *durations changing* between commits is the signal that the failure
changed identity even though the status did not. Concretely: after fixing a CI break, **watch the
next run to completion and read the conclusion** rather than assuming green; treat "red again" as
a fresh diagnosis, not a continuation. And when you fix a gate, ship the fix so the *intermediate*
commit is green too — splitting a format repair so one file lands in a later commit leaves a red
commit in the history for no benefit, which is a mistake I made in this very sequence.

Reordering the gates is usually the wrong fix — putting the suite first costs its full wall-clock
on every trivial lint slip. The mitigation is the reading habit, not the config.
