---
name: lessons-records
description: "The archive and the session-boot surfaces: stale breadcrumbs, refusals that read as gaps, doc rot, caps that become targets, and shell/commit mechanics"
metadata:
  node_type: memory
  type: feedback
---

**Lessons — Records, docs, memory & tooling.** Split out of `.claude/memory/MEMORY.md` on 2026-08-26; that file's
index points here by path. These rows carry **no `MEMORY.md` row of their own**, so they cost
nothing until this file is read — the same arrangement `prohibitions/` has had since D-185.

**Read this file before you write the code, not after the review.** Each row is *the trap* +
*what to do instead*; the measurement that earned it is in the linked file. If a row looks
like pedantry, open its file — every one of them cost a beat.

- [CRLF joins inflate line count](../feedback-crlf-join-inflates-line-count.md) — splitting a CRLF file on `
- [Commit messages corrupted by shells](../feedback-no-powershell-heredoc-in-bash.md) — PS here-strings in Bash, and `Out-File -Encoding utf8`'s BOM; exit 0 proves nothing — and verify with a binary read, since `git cat-file | od -c` invents CRLF in the pipe and reported 118 CRs in a message with zero
- [Enumerate competitors before timing](../feedback-enumerate-competitors-before-timing.md) — other agent sessions run 26-worker pytest on this box; the same suite measured 119 s and 363 s, so prefer counts over durations
- [Memory is one tracked copy](../feedback-memory-is-one-tracked-copy.md) — the C: harness path is a junction into `.claude/memory`; never sync or mirror it, there is nothing to sync
- [A cap being written to can't be raised](../feedback-a-cap-being-written-to-cannot-be-raised.md) — pinned at exactly 250 for 13 commits, then at 300 for 8 more with every block under cap: a per-item cap can't rescue a total and retirement couldn't pay for it, so delete the threshold and report the number (D-177)
- [A doc rots where it duplicates](../feedback-a-doc-rots-where-it-duplicates.md) — ARCHITECTURE.md decayed only where it restated DECISIONS.md; ask what surface OWNS each claim. Suspect the maintenance rule before the discipline (D-184)
- [Measure which END is growing](../feedback-measure-which-end-is-growing.md) — the live frontier held at 6-17 for two months while the settled tail went 0→22; four mechanisms failed, all aimed at the SUM. Age ≠ finishedness, so when the growing end is incompressible the lever is granularity, not eviction (D-185)
- [A refusal reads like a gap](../feedback-a-refusal-reads-like-a-gap.md) — in a copied-forward "still open" list, "nothing sources it" and "REFUSED because the fix would be a guess" look identical; open the ORIGINATING record. Nearly led with a decision D-149 had already taken (D-186)
- [Check the blocker is still blocking](../feedback-check-the-blocker-is-still-blocking.md) — unrelated work demolishes a premise and leaves NO ⚠; D-108's item rode the open lists 20 times after its blocker shipped. Probe the code first: it turned a build into a measurement
- [A scope note can carry a mechanism claim](../feedback-a-scope-note-can-carry-a-mechanism-claim.md) — "not built because <mechanism>" is TWO claims, and the mechanism inherits the scope decision's immunity from review
