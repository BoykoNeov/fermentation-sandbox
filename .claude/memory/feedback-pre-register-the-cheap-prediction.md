---
name: feedback-pre-register-the-cheap-prediction
description: "Before an expensive campaign, write the free static prediction to disk and score it afterwards — it is what makes a surprise claimable"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 92f44b73-555b-4c7e-a475-7e54c360aef1
  modified: 2026-07-29T06:21:18.550Z
---

Before a long measurement campaign, spend the **cheap static pass first**, write its
prediction to disk, and score it against the results afterwards. Without that, "the
campaign found X" is indistinguishable from "I read the output and built a story around
whatever pattern it happened to show".

**Why:** D-163 ran 73 mutation arms (~2.5 h) over every band edge in the archive. A free
`grep '\.uncertainty' tests/` — 39 hits across 10 files — was written to `PREDICTION.md`
*before any arm ran*. It scored its negative half **exactly**: 13 files / 141 bands
predicted GREEN in both directions, 26 consecutive GREENs measured. It also called 4 of 5
RED files, including `sensory`, where getting it right required predicting the **direction**
(green on shrink, red on widen), not just the file.

That pre-registration bought two things nothing else could:

- **The one miss became legible.** `psychophysics` was predicted RED and came back GREEN.
  Because the prediction named *why* it should go red (a disjointness constraint), the miss
  could be diagnosed as an **operator limitation** — scaling each half-width about its own
  nominal preserves overlap by construction, so that constraint *cannot* fire — rather than
  being silently absorbed as "measured GREEN, therefore unguarded". Post-hoc, that cell would
  have read as coverage.
- **The surprise became claimable.** The campaign's real finding (an uncertainty band
  silently gating a scenario override) was an *unpredicted* RED. Claiming "static analysis
  could not have seen this" is only honest if you wrote down beforehand what static analysis
  did see.

**How to apply:** Do the free pass first and commit it to disk before the expensive one
starts. State expected outcomes **per arm and per direction** whenever the operator is
asymmetric — a prediction that only names files is half a prediction. Then score every row,
misses included, in the record itself. The interesting cell is a **predicted-GREEN that comes
back RED**: that is precisely what the cheap method cannot see, and it is where the finding
usually is. Name the operator's blind spot in advance too, so its scope is a stated bound
rather than a later concession [[feedback-conceded-caveats-are-not-coverage]]. Related:
[[feedback-mutate-the-premise-before-building-the-guard]] (run the mutation before building
the guard), [[feedback-count-and-print-your-skips]] (the denominator the prediction is scored
against), [[feedback-a-majority-is-not-a-direction]] (fix the stopping rule before the run
lands).
