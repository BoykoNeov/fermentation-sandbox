---
name: feedback-always-commit-push
description: "User wants commit + push to main on every completed change, not just at batch/session end"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ddd3bb59-ccfa-42d4-beb5-35acc4eeb6b4
  modified: 2026-07-20T15:07:20.991Z
---

When a change is complete, **always commit and push to main** — don't wait for
an explicit "commit" request or for batch/session end.

**Why:** User said "always commit and push" after a standalone license-relicense
task, generalizing it beyond the [[feedback-batch-end-ritual]] (which fires on
batch/session end). They want durable, pushed state per completed task.
**Repeated with visible frustration on 2026-07-20** ("i say this again and
again") after a session opened by *reporting* an uncommitted finished D-118 tree
as a next-step option instead of just committing it. **Repeated again the same
day** ("i have repeated this countless times ... dont wait for test suites")
after a session opened by launching `ruff && mypy && pytest` on an uncommitted
finished D-121 tree — the suite ran past its 600s timeout and the commit was
still unmade.

**How to apply:** The trigger is "work is done / decision made / memory
updated" — commit and push to `main` **at that moment**. Do **not** gate the
commit on `pytest`, `mypy`, `ruff`, or any build: run checks *after* pushing,
and push a follow-up fix if they fail. A long suite must never sit between
finished work and a durable commit. Still run the full
[[feedback-batch-end-ritual]] (memory + docs update) at batch/session end.
Branch first only if explicitly asked. See [[feedback-never-pipe-checks-to-tail]]
for how to run those checks once they are unblocking.

**Session-start corollary:** a dirty tree holding *completed* work is not a
decision to surface — verify and ship it. Only genuinely unfinished or
ambiguous work belongs in a "what's next" list.
