---
name: feedback-size-a-set-with-the-instrument-that-can-see-it
description: Estimating how big an unenumerated set is with the instrument that already failed to find its known members re-imports the blind spot into the estimate
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1a6ea7a3-7bda-4049-b556-1015cdbe0204
  modified: 2026-08-26T12:24:18.948Z
---

Before enumerating a set, it is natural to size it cheaply so the work can be scoped. **The cheap
sizing instrument must not be the one whose blind spot created the set.** If it were reliable, the
set would already be enumerated.

**Why:** at D-234 the task was "every parameter read at compile that is also sampled". I sized it
with a grep for `parameters["…"]`, got **26** compile reads and predicted **12-20** members. The
answer was **32**. The 21 names the grep could not see — the 19 `pKa_*` and both
`nitrogen_uptake_charge_*` — are read off `parameters.resolve()` inside a generator expression,
and that block contains the *entire* defect the previous beat (D-233) had just repaired. So the
estimate was wrong by more than the census was large, and it was wrong in the one direction that
would have justified skipping the beat as small. A sizing that reports "one beat, not five" is a
scoping decision, and this one was made by the blind instrument.

**How to apply:** when a set exists *because* something went unnoticed, ask what noticed it late
and why. Then size with an instrument of a different kind — instrumentation over source-matching,
runtime recording over static scan — or state the estimate as a **floor** and say which mechanism
it cannot see. A floor is honest and still scopes the work; a point estimate from the failed
instrument is the failure re-imported. Pair the sizing with the **known** members: if your
instrument cannot re-find the two members you already know about, it cannot bound the rest
[[feedback-grep-finds-claims-not-guards]] [[feedback-count-and-print-your-skips]].
