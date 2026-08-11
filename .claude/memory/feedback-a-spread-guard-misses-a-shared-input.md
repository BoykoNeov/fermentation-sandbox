---
name: feedback-a-spread-guard-misses-a-shared-input
description: "A wrong shared input moves the MEAN, not the spread, so a spread-ratio guard is blind to the defect it was built for"
metadata:
  node_type: memory
  type: feedback
---

Before pinning a **spread**, check which moment of the distribution the defect you fear actually
moves. At D-186 the hazard was an anchor solving from a *different* parameter map than the run
integrates against, and the obvious guard was "this anchor's off-nominal spread may not exceed the
pre-existing anchor's by more than 30 %". It was sized across three member counts (0.962 / 0.988 /
1.015) and then falsified by planting exactly that defect. A pKa shift that threw the **nominal**
anchor **0.29 pH** off target moved the ratio only **0.988 → 1.119** — green at any threshold that
would not also fire on noise.

**Why:** the failure is structural, not a badly chosen constant. A wrong input that is **shared by
every member** displaces them all by nearly the same amount, so the error lands in the **mean** and
the **spread** barely notices. A spread statistic is therefore blind to precisely the class of bug
that motivates it, and blind in the reassuring direction — the guard sits green in the suite and
reads as coverage. What did detect the planted defect was the plain nominal-exactness assertion,
which already existed. See [[feedback-grep-finds-claims-not-guards]] and
[[feedback-mutate-the-premise-before-building-the-guard]] for the same shape at other stages.

**How to apply:** falsify every guard against the specific defect it names, *before* shipping it,
and read the verdict rather than the redness of something nearby. **If it stays green, the honest
output is the measurement without the guard** — record the number and the rejected design, do not
weaken the threshold until it fires, and do not ship it anyway "for coverage". Per-member or
per-arm statistics need the same question asked in reverse: an error that varies per member will
show in the spread and can hide in the mean.
