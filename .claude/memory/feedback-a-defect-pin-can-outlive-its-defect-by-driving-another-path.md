---
name: feedback-a-defect-pin-can-outlive-its-defect-by-driving-another-path
description: "A guard written as 'a RED here means it was fixed' stayed GREEN through its own repair, because its arms drove a path the fix does not sit on"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 217f4393-f3f4-48e3-9bef-4bb21f5faae0
  modified: 2026-08-26T13:36:37.391Z
---

This repo pins known defects on purpose, in D-233's idiom: *a RED here means the defect was
repaired — delete the guard and say so in the record, never revert the repair.* That idiom is
only sound for a pin whose arms actually **reach** the code the repair changes. Check that
before you read either colour.

D-234 pinned two LIVE defects. Repairing both at D-235/D-236, one pin went red exactly as
written and was deleted for a positive-form replacement. The other **stayed green** — its arms
drove `simulate_scheduled` by hand on the compiled `y0` with a patched parameter map, and the
repair lives in `y0_for_member`, which that path never calls. Nothing was wrong with the guard
and nothing was wrong with the repair; they simply did not intersect.

**Why:** two opposite errors sit on that green, and both look reasonable in the moment. Reading
it as *the repair was inert* is wrong — the defect was measured in the ensemble frame (16.65 %
→ bit-identical) and the hand-wired path is just a different caller. Deleting it anyway on the
strength of the instruction is also wrong — its mechanism claim is still true, and it carried
the only guard on a null (`f(Cu)` moves nothing when both roles move coherently) that would
otherwise have gone unwatched.

**How to apply:** when you write a defect pin, name in its docstring **which entry point** it
drives, so a later beat can tell at a glance whether a repair should move it. When a repair
lands and the pin stays green, do not treat the colour as a verdict on the repair: trace the
call path, then either re-scope the guard to the mechanism it really pins (renaming it, so the
name stops promising a defect that is gone) and add a **new** guard in the frame the repair
lives in, or delete it if the claim died with the defect. Evidence for a repair has to come
from a test that exercises the repaired path — never from an old guard's colour
([[feedback-a-control-needs-mechanical-reach]], [[feedback-run-the-mutation-the-claim-names]]).
