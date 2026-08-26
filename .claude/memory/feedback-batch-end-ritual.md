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

**How to apply:** at batch end, in order: (1) write/refresh memory files here, then add the
row. **A new epistemics lesson gets its row in `.claude/memory/lessons/<group>.md`, NOT in
`MEMORY.md`** (2026-08-26 split). `MEMORY.md` keeps ten rows only: the user, the project
ledger, and the workflow rules that fire every session — the hook's `BOOT_ROWS` is the list,
and it flags a `feedback-*` row that lands in the index instead. **Note the harness's own
memory instructions say "add a one-line pointer in MEMORY.md" — that is the pre-split rule and
following it verbatim regrows the file; the hook exists because this sentence is not enough.** **The full per-decision narrative goes in `docs/DECISIONS.md`
ONLY — memory gets a *distilled status update / pointer*, never a copy.** Do not
append the batch write-up to [[project-fermentation-sandbox]] or to a `MEMORY.md`
index line; instead update the project file's short status block and bump the
current D-number. Index rows stay one line, 320 bytes hard. This is
the guardrail that keeps memory from re-bloating — the decision log is the archive,
memory is the boot context. **This rule as prose is NOT sufficient — it has already
failed once.** It was added by `acd3ce1` (2026-07-02) claiming to fix "the cause, not
just the symptom"; the project file still regrew 2.4KB → 277KB (**114×**, a full
D-38→D-111 changelog in two formats) by 2026-07-17. The mechanism that now backs it:
`.claude/hooks/check_memory_size.py`, a PostToolUse hook (project `.claude/settings.json`).
It caps **per block** (8 lines in the project memory, 14 in `CLAUDE.md`, 320 BYTES per index
row in `MEMORY.md` *and* `lessons/`) and **REPORTS every surface's total without capping it**
— including `MEMORY.md`'s headroom against the harness load limit, which it silently truncates
at; that went unmeasured for months because the hook counted lines and the loader counts bytes — the whole-file
cap was removed at **D-177** after being raised four times and re-pinned at each value; never
put one back [[feedback-a-cap-being-written-to-cannot-be-raised]]. It *detects*; it cannot
enforce distillation — that is still a judgement call at ritual time. If the warning fires,
**distil the NEW block**, never evict an old prohibition. (1b) **the tracked copy needs no
sync step** — the harness path `~/.claude/projects/M--claud-projects-Fermentation/memory/` is a
**junction into `.claude/memory/`**, so writing memory writes the repo file directly; there is
exactly one copy and nothing to `cp` [[feedback-memory-is-one-tracked-copy]]. (The 2026-07-17
"52 decisions stale" incident predates the junction; the `cp ... .claude/memory/` step this file
used to prescribe is now a no-op that can only mask a real problem.) `git add .claude/memory/`
still applies: the user asked on 2026-06-23 that memory be tracked *always, with the rest*,
not a one-off snapshot, so it is committed alongside docs/code every checkpoint and behaves
like any other tracked file; (2) update affected docs (`docs/ARCHITECTURE.md`,
`docs/DECISIONS.md` — **never `docs/plans/milestone-*.md`, frozen logs since D-184**); (3) `git commit` with
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
