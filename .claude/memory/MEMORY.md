# Memory index — Fermentation Sandbox

- [User: BoykoNeov](user-boykoneov.md) — research engineer, owns the fermentation sim; fidelity over convenience
- [Batch-end ritual](feedback-batch-end-ritual.md) — on batch/planning/session end: distilled status to memory, full entry to DECISIONS.md, commit + push
- [Always commit + push](feedback-always-commit-push.md) — commit and push to main the moment work is done; never gate the commit on pytest/mypy/ruff
- [Discuss disagreements](feedback-discuss-disagreements.md) — surface design disagreements before building; specs aren't gospel
- [Project: Fermentation Sandbox](project-fermentation-sandbox.md) — status + repo; M0/M1/M2 complete, Milestone-3/Tier-3 aging in progress at D-127. Per-decision archive: docs/DECISIONS.md
- [Never pipe checks to tail](feedback-never-pipe-checks-to-tail.md) — `cmd | tail && ...` returns tail's exit 0 and hides ruff/pytest failures
- [Full suite before "green"](feedback-full-suite-before-green.md) — a new Process in a shared registry breaks exact-set + end-to-end tests outside the domain suite; don't claim green until full pytest passes
- [Verify latest state, not breadcrumbs](feedback-verify-latest-state-not-breadcrumbs.md) — old "Next:" lists get burned down by later decisions; check the latest D-record + code before proposing work
- [Rejected values must be unreachable](feedback-rejected-values-must-be-unreachable.md) — an unphysical value in a sampled field is a live defect a green suite won't catch
- [Name guards for what they forbid](feedback-name-guards-for-what-they-forbid.md) — a guard mislabelled with a real mechanism invites a fair objection and gets argued away
- [Measure which side before building](feedback-measure-which-side-before-building.md) — a one-directional corrective only helps in one direction; check the sign, the reach, and whether it's a rate knob on a supply-limited quantity
- [No PowerShell here-strings in Bash](feedback-no-powershell-heredoc-in-bash.md) — `@'...'@` in the Bash tool silently corrupts commit messages; exit 0 proves nothing
- [Best-practices reference](reference-claude-best-practices.md) — rosmur claudecode-best-practices URL to apply
