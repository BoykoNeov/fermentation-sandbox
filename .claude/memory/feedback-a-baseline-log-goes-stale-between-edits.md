---
name: feedback-a-baseline-log-goes-stale-between-edits
description: "A saved suite log stops being a baseline the moment any other edit lands; diffing against it credited an unrelated change to my fix, which had in fact repaired nothing"
metadata:
  type: feedback
---

A saved test-suite log is a baseline **only for the exact tree that produced it**. At D-223 I
predicted a fix would take the suite 16 failures → 15 and scored it a hit off a diff against a
`suite_after_q.txt` written earlier in the same session. It was written *before* an unrelated
edit (a benchmark's strict `xfail` decorator being removed), and that edit is what closed the
one test. **The fix repaired zero tests** — the guard it was aimed at correctly still failed —
and for a while it was invisible to the entire suite.

**Why:** a count is not an attribution. "16 → 15" is consistent with "my change fixed one",
"my change fixed two and broke one", and "something else fixed one and my change did nothing",
and the diff of failure *names* only separates those if the baseline is from the same tree. A
prediction that lands on the right number for the wrong reason is worse than a miss, because it
retires the question. Compare [[feedback-a-hit-can-be-two-errors-cancelling]].

**How to apply:** before diffing against a saved log, check `git status`/`git diff --stat`
against what the log's tree contained; if anything else moved, re-run the baseline. Prefer
diffing failure NAMES over counts, and when a name flips, confirm the mechanism (which assert,
which value) rather than accepting the flip. If a fix turns out to repair nothing the suite can
see, that is not a null result — it is the fix **owing a guard**
([[feedback-prefer-the-variant-your-guards-can-see]]), and the guard should be falsified with
mutation arms before you believe it.
