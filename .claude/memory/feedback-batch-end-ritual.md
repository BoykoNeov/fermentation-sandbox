---
name: feedback-batch-end-ritual
description: "What to do when a work batch or planning ends, or the user says \"session end\""
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
---

When a work batch **or a planning session** ends — and explicitly when the user
says **"session end"** — run the full close-out ritual: **update memory, update
docs, commit, and push to `main`.** The user reaffirmed this directly on
2026-06-21, adding "planning" as a trigger alongside batch end / "session end".

**Why:** the user stated this directly as a standing instruction for the project.
It keeps the repo, the docs, and cross-session memory in sync so the next session
resumes cleanly. In the first session I completed docs+commit+push but initially
forgot the memory step — don't repeat that.

**How to apply:** at batch end, in order: (1) write/refresh memory files here +
update `MEMORY.md`. **The full per-decision narrative goes in `docs/DECISIONS.md`
ONLY — memory gets a *distilled status update / pointer*, never a copy.** Do not
append the batch write-up to [[project-fermentation-sandbox]] or to a `MEMORY.md`
index line; instead update the project file's short status block and bump the
current D-number. `MEMORY.md` index entries stay one line (~200 chars). This is
the guardrail that keeps memory from re-bloating — the decision log is the archive,
memory is the boot context. **This rule as prose is NOT sufficient — it has already
failed once.** It was added by `acd3ce1` (2026-07-02) claiming to fix "the cause, not
just the symptom"; the project file still regrew 2.4KB → 277KB (**114×**, a full
D-38→D-111 changelog in two formats) by 2026-07-17. The mechanism that now backs it:
`.claude/hooks/check_memory_size.py`, a PostToolUse hook (project `.claude/settings.json`)
that warns when `project-fermentation-sandbox.md` exceeds **150 lines**. It *detects*;
it cannot enforce distillation — that is still a judgement call at ritual time. If the
warning fires, distil rather than raise the cap. (1b) **sync the repo's tracked copy** — the memory files are
version-controlled under `.claude/memory/` (since commit `1c095ab`), so refresh
them from the live dir before committing:
`cp ~/.claude/projects/M--claud-projects-Fermentation/memory/*.md .claude/memory/`
then `git add .claude/memory/`. **This step drifts silently — check it every time.**
On 2026-07-17 the tracked copy was found **52 decisions stale** (last committed at
D-59, `33cb5ab`, 2026-07-08) while the live file had grown 6× past it. Nothing catches
this: the live memory lives *outside* the repo, so a clean `git status` is not evidence
the tracked copy matches — it only means the stale copy is committed. Diff the two, don't
trust the tree. This is the durable mechanism (the user asked on
2026-06-23 that memory be tracked *always, with the rest*, not a one-off snapshot)
— so memory is committed alongside docs/code every checkpoint, behaving like any
other tracked file; (2) update affected docs (`docs/ARCHITECTURE.md`,
`docs/DECISIONS.md`, the `docs/plans/milestone-*.md` trio); (3) `git commit` with
Conventional Commits; (4) `git push`. Run the ritual even when there is no code to
push. See [[project-fermentation-sandbox]].

**Always push directly to `main`, and do NOT ask first** — this is a solo public
repo with CI on `main` and no PR flow; do **not** branch or open a PR for routine
work. The user stated "always push to main" (2026-06-20) and reaffirmed **"always
commit and push"** (2026-06-23) after I paused to confirm a push at the end of the
D-16 work. So: completing a self-contained piece of work *is* a batch end — commit
and push without a confirmation prompt. Don't treat "the user is present" or "this
was my inference of batch end" as reasons to hold the push; the standing
authorization covers it. NOTE: the auto-mode guardrail blocked a direct `main` push
once because it read "commit and pause" as not authorizing it; a `Bash(git push:*)`/
`Bash(git push origin main:*)` allow-rule in settings is the durable fix so it
doesn't re-prompt.
