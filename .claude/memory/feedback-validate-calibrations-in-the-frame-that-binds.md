---
name: feedback-validate-calibrations-in-the-frame-that-binds
description: A constraint that only binds under competition cannot be validated on an isolated Process — solve-in-isolation plus test-in-isolation is self-sealing
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6b5be9a9-1775-4a76-9e7e-75bdeeec42f6
  modified: 2026-07-28T03:31:30.811Z
---

**A calibration solved in isolation and then tested in isolation proves nothing about the model
that ships.** If a constraint only binds *because other consumers compete for the same pool*,
removing the competition removes the constraint's entire content — and the test still passes.

**Why:** D-133 jointly solved two constants from two Ferreira constraints, doing the arithmetic on
a single Process that got all the O₂ it needed, then wrote unit tests the same way. The exhaustion
test set `o2 = 0.5 g/L` — 500 mg/L, ~62× air saturation — explicitly so the pool's decay was "not
confounded by O₂ itself running out". Running out **was** the physics: end-to-end the sink wins
~35% of a real 8 mg/L charge and plateaus at 39.2% left, so the constraint fails structurally. It
stayed green for four decisions. The failure mode is self-sealing — the operating point that makes
the test pass is the one that deletes what it was testing.

**How to apply:** Before trusting a fitted constant, run its own source's constraints **through the
shipped ProcessSet at a reachable operating point**, and at the *source's* operating point too
(judging a calibration outside its stated condition is a separate error). Watch for a test comment
that justifies an extreme input as removing a "confound" — check whether the thing being removed is
the binding constraint. When constraints jointly pin N constants and one fails, say precisely what
survives: usually a **product** is pinned and the **split** is not, and re-solving the split is a
second unmeasured number propping up the first. Relabel such a test rather than deleting it — its
rate-law content is usually real, only its headline claim is false.
Related: [[feedback-rejected-values-must-be-unreachable]], [[feedback-name-guards-for-what-they-forbid]],
[[feedback-conceded-caveats-are-not-coverage]].
