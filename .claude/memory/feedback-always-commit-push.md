---
name: feedback-always-commit-push
description: "User wants commit + push to main on every completed change, not just at batch/session end"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ddd3bb59-ccfa-42d4-beb5-35acc4eeb6b4
---

When a change is complete, **always commit and push to main** — don't wait for
an explicit "commit" request or for batch/session end.

**Why:** User said "always commit and push" after a standalone license-relicense
task, generalizing it beyond the [[feedback-batch-end-ritual]] (which fires on
batch/session end). They want durable, pushed state per completed task.

**How to apply:** After finishing a discrete unit of work and verifying it
(tests/build green), stage the relevant files, make a Conventional Commit, and
push to `main`. Still run the full [[feedback-batch-end-ritual]] (memory + docs
update) at batch/session end. Branch first only if explicitly asked.
