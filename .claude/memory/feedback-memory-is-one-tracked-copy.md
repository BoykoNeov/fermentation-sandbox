---
name: feedback-memory-is-one-tracked-copy
description: "The memory store is ONE copy, tracked in the repo; the C: harness path is a junction, not a second store"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 28f4dd80-13f3-4d05-a1a4-7323c3a2bb8a
  modified: 2026-07-28T11:27:28.117Z
---

Every memory file must exist as a **single tracked copy inside the project
folder**. Since 2026-07-28 the harness path
`C:\Users\boiko\.claude\projects\M--claud-projects-Fermentation\memory` is a
**directory junction** pointing at `M:\claud_projects\Fermentation\.claude\memory`.
The two paths are the same bytes on disk. Writing through either one lands in
the repo and dirties `git status` immediately.

**Why:** they were two independent copies that silently diverged twice. Commit
`bc77e4f` hand-synced them once; by 2026-07-28 the tracked mirror was 6 days
stale — 6 memories missing outright, 4 more stale, and
`project-fermentation-sandbox.md` read 178 lines against the live 250. Opening
the in-repo copy therefore showed week-old prohibitions as current, which is
exactly the failure [[feedback-verify-latest-state-not-breadcrumbs]] exists to
catch. A mirror kept in step by hand is not a mechanism.

**How to apply:** never "sync", copy, or mirror the memory folder — there is
nothing to sync. Do not add a copy hook. If `git status` ever shows
`.claude/memory` and the C: path disagreeing, the junction has been replaced by
a real directory (the harness recreating the path would do this silently) —
restore it with `New-Item -ItemType Junction` rather than resuming the copying.
Memory writes now land in the working tree, so they get committed with the batch
under [[feedback-batch-end-ritual]] and [[feedback-always-commit-push]].

Scope: Fermentation only. The other ~60 project stores under
`C:\Users\boiko\.claude\projects\` are untouched real directories.
