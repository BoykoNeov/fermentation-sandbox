---
name: feedback-a-rename-can-restore-green-without-the-claim
description: "Repairing a test that broke on a rename can restore GREEN while leaving it unable to see the change the rename was part of; check the fixed test still discriminates"
metadata:
  node_type: memory
  type: feedback
---

When a test breaks because a symbol it imports was renamed, **fixing the name is not finishing
the repair.** Ask what the rename was PART of, and whether the repaired test can still tell that
change from its opposite. Often it cannot, and the red was masking a guard that forbade nothing.

**Why:** D-229 moved an index-row cap from characters to BYTES — the whole point of the record,
because the harness loader counts bytes and was silently truncating — and renamed the constant.
Its test module still imported the old name and raised `AttributeError` (1 failure + 9 fixture
errors) for two records. The obvious repair is the rename. But **every row those tests build is
ASCII**, so `len(row) == len(row.encode())` and all of them pass identically under either unit: a
revert to counting characters would have been invisible to the module written to guard it.
D-231 shipped the rename **with** a row of em-dashes at 0.75x the cap in characters and 2.25x in
bytes. Falsification: reverting the cap to characters turned exactly that one test red and left
28 green — and the 28 ARE the finding.

**How to apply:** after fixing a rename-induced break, construct the input on which the old and
new semantics DISAGREE and assert the new one. If you cannot construct one, the test does not
test the rename's subject and you should say so rather than bank the green. Assert the probe's
preconditions (under one measure, over the other) before the behaviour, so a probe that stops
straddling fails loudly instead of passing vacuously. Related:
[[feedback-name-guards-for-what-they-forbid]], [[feedback-full-suite-before-green]],
[[feedback-a-red-set-firing-wider-is-a-finding]].
