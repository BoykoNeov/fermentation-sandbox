# Memory index — Fermentation Sandbox

**The 99 epistemics lessons are NOT rows in this file.** They live in `.claude/memory/lessons/`
— five files reached BY PATH, carrying no index row of their own, so they cost nothing until
read. Same arrangement `prohibitions/` has had since D-185, applied here on 2026-08-26.
**Open the matching lessons file BEFORE you build, verify, measure or source — not after the
review.** A new lesson gets a row in a lessons file, **never here**: row COUNT was the channel
that overflowed, and this file only stays small while it stops growing with the record count.

Why not just trim the rows: that was tried on 2026-08-19 (17 rows, 25840 → 24301 bytes) and was
back over the limit in seven days. 35 % of the old file was link-and-title text no edit can
shrink, and every new row adds ~90 bytes of it whatever the prose discipline.

## Lessons — by path, not loaded at boot
- [Guards, tests & mutation arms](lessons/guards.md) — What a test actually catches: designing a mutation, reading what a RED names, pinning bands and tolerances, and the guards that quietly forbid nothing (29 lessons)
- [Measuring a run & attributing the number](lessons/measurement.md) — Reading a number off the model and deciding what caused it: baselines, controls, censuses, sampling grids, summary statistics and null results (23 lessons)
- [Parameters, bands, fits & identifiability](lessons/parameters.md) — Where a constant's freedom lives: per-parameter bands vs joint claims, calibrated levels that decay, rate-law coordinates, and fitting vs scoring (22 lessons)
- [Sources, literature & transcription](lessons/sources.md) — Getting a number out of a paper intact: access, tables vs prose, units, measurement frames, the scope of a negative, and notes that drift from the source (14 lessons)
- [Records, docs, memory & tooling](lessons/records.md) — The archive and the session-boot surfaces: stale breadcrumbs, refusals that read as gaps, doc rot, caps that become targets, and shell/commit mechanics (11 lessons)

## Loaded every session
- [User: BoykoNeov](user-boykoneov.md) — research engineer, owns the fermentation sim; fidelity over convenience
- [Batch-end ritual](feedback-batch-end-ritual.md) — on batch/planning/session end: distilled status to memory, full entry to DECISIONS.md, commit + push
- [Always commit + push](feedback-always-commit-push.md) — commit and push to main the moment work is done; never gate the commit on pytest/mypy/ruff
- [Discuss disagreements](feedback-discuss-disagreements.md) — surface design disagreements before building; specs aren't gospel
- [Closer to reality decides](feedback-closer-to-reality-decides.md) — never offer "hold it unchanged" as a peer to a more faithful option; ship the faithful one and record the size. A surprising magnitude means re-check the reasoning, not retreat
- [Project: Fermentation Sandbox](project-fermentation-sandbox.md) — status + repo + a **LEDGER of settled subjects**, at D-229. **Read it before proposing any Milestone-3 work** — closed work reads as unbuilt from outside; detail sits in `prohibitions/`, by path, not indexed here. Archive: docs/DECISIONS.md
- [Never pipe checks to tail](feedback-never-pipe-checks-to-tail.md) — `cmd | tail && ...` returns tail's exit 0 and hides ruff/pytest failures
- [Full suite before "green"](feedback-full-suite-before-green.md) — a new Process in a shared registry breaks exact-set + end-to-end tests outside the domain suite; don't claim green until full pytest passes
- [Verify latest state, not breadcrumbs](feedback-verify-latest-state-not-breadcrumbs.md) — old "Next:" lists get burned down by later decisions; check the latest D-record + code before proposing work
- [Best-practices reference](reference-claude-best-practices.md) — rosmur claudecode-best-practices URL to apply

