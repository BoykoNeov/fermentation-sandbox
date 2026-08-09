---
name: feedback-verify-the-restore-between-mutation-arms
description: "In a mutation matrix, prove the baseline came back between arms — a silent restore failure makes the arms you expect to be RED confirm themselves"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1289a7da-873a-4fc1-882a-f8c7f961f6e7
  modified: 2026-08-09T11:30:48.069Z
---

**Between mutation arms, verify the file actually got restored — and include at least one arm whose
expected outcome is the opposite of the others.** In D-158 the backup used
`cp $F /tmp/f.bak 2>/dev/null || cp $F <durable>`; git-bash has its own `/tmp`, so the *primary* `cp`
succeeded, the fallback never ran, the later restore `cp` failed on a nonexistent path, and arms B and
C both ran against a file still carrying arm A's mutation. All three reported RED — and RED was the
expected answer for two of them.

**Why:** a mutation matrix is graded against expectations, so any arm whose expectation is "red"
cannot distinguish a working guard from a broken harness. The corruption was caught only because one
arm's expected answer was GREEN (a consistent re-sourcing must pass) and it came back red. Without
that arm the matrix would have shipped looking complete. Same family as
[[feedback-count-and-print-your-skips]] — a harness that fails quietly returns the answer you already
believed.

**How to apply:** back up to an explicit, verified path under `M:\claud_projects\temp` (never a bare
`/tmp`, and never a `||` fallback whose primary can succeed somewhere unintended); `ls` the backup
before trusting it. Re-assert the baseline between arms — `git status --porcelain` clean, or grep the
value back — and run the untouched baseline as its own arm at both ends. Design at least one arm to be
green: for a guard that pins a derivation, the green arm is a *consistent* re-sourcing (move the source
number and the shipped number together), which is also the arm that proves the guard pins the
derivation rather than the literals. See [[feedback-mutate-the-premise-before-building-the-guard]].

**Generalises past mutation matrices to any scored harness (D-164).** The first run of a
*pre-registered repro* reported 9 of 9 arms raising, with two predictions marked AS PREDICTED — and
every arm had raised for one unrelated reason (the scenario carried no `temperature_schedule`, so
nothing ever reached the code under test). Two predictions confirmed themselves against a harness that
never ran the experiment. The three designed-GREEN control arms are the only thing that showed it.
So: **any arm whose expected verdict matches the harness's failure mode is unscored until a control
proves the harness works** — for exception-counting harnesses the failure mode is "raises", so the
control must be an arm that must NOT raise. Write the controls into the harness from the first run,
not after a suspicious result. See [[feedback-compute-the-clean-fix-before-adopting-it]].
