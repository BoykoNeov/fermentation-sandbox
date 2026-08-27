---
name: feedback-a-guard-must-be-scored-where-its-subject-lives
description: "A guard scored in the one configuration where its subject cannot exist forbids nothing, and passes hardest exactly when the thing it watches for happens"
metadata:
  node_type: memory
  type: feedback
---

A guard that watches for a change to some Process, wiring or medium must be **scored in a
configuration where that thing is actually present**. Score it in the default configuration out of
habit and it can be structurally blind: the subject is absent, the property is vacuously true, and
the test is greenest precisely when the change it exists to catch has landed.

**Why:** at D-240 I shipped `test_the_priced_names_are_still_undrawn`, whose whole job is to go RED
if a later beat wires `burst_antioxidant_initial` into a Process's `reads`. The mutation arm did
exactly that — and **all 46 tests passed**. The guard compiled the scenario with the *default*
`oxidative="direct"` set, where `AntioxidantBurstOxidation` is not wired into the medium at all
(D-147 ships the burst non-default), so its `reads` never reach `_schedule_reads` and the name
stayed undrawn no matter what the Process declared. The pin existed, was well argued, was
committed to a record, and forbade nothing. Only the arm found it; a green suite never would.

**How to apply:** for every guard, name the configuration its subject lives in and score it there —
here, parametrise the compile on the wiring the row belongs to rather than taking the default. Two
habits follow. **Run the arm that performs the change the guard's own error message describes**, in
the words of that message; a guard whose stated trigger has never been exercised is a claim, not a
test. And when an arm comes back **all-green where you predicted RED, suspect the guard before the
prediction** — a passing mutation arm is usually the harness failing to reach, not the code being
robust. Cf. [[feedback-a-defect-pin-can-outlive-its-defect-by-driving-another-path]], the same
blindness arriving through the code path rather than the configuration.
