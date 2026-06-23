---
name: reference-claude-best-practices
description: Claude Code best-practices guide the user wants applied to this project
metadata: 
  node_type: memory
  type: reference
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
---

The user asked to apply reasonable best practices from
https://rosmur.github.io/claudecode-best-practices/ to the Fermentation project.

Applied in M0: concise root `CLAUDE.md` (<2k tokens, conventions + commands),
three-file planning docs (`docs/plans/milestone-1-{plan,context,tasks}.md`),
test-driven benchmarks, Conventional Commits, CI quality gates (ruff/mypy/pytest).

Deferred / not done: a `.claude/settings.json` permission allowlist — the harness
auto-mode classifier blocks me from writing one (treated as self-modification).
The user can add it themselves or run `/fewer-permission-prompts`. See
[[project-fermentation-sandbox]].
