---
name: feedback-no-powershell-heredoc-in-bash
description: "PowerShell here-strings (@'...'@) silently corrupt commit messages when used in the Bash tool — use a real bash heredoc instead"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 4c5542b6-994d-42ff-9b5f-a6dbc7d14d50
  modified: 2026-07-20T09:29:50.930Z
---

This environment exposes **two shells with different syntax**: the PowerShell
tool (Windows PowerShell 5.1) and the Bash tool (Git Bash / POSIX sh). A
PowerShell here-string `@'...'@` passed to the **Bash** tool does not error — bash
treats the `@'` and `'@` as literal text, so the command *succeeds* and the
corruption lands in the artifact.

**Why:** at D-117 a `git commit -m @'...'@` in the Bash tool produced the subject
`@ fix: D-117 follow-up -- ...` with a stray trailing `@` in the body. Exit code
0, push succeeded, nothing failed. The leading `@` broke the repo's Conventional
Commits rule, and fixing a *pushed* commit message needs `--amend` plus a
force-push to `main` — a history rewrite for a formatting slip. Same failure
family as [[feedback-never-pipe-checks-to-tail]]: **a zero exit code is not
evidence the command did what was intended.**

**How to apply:** in the Bash tool use a real bash heredoc —
`git commit -F - <<'EOF' ... EOF` (quoted `EOF` keeps `$` and backticks literal).
In the PowerShell tool use `@'...'@` with the closing `'@` at column 0. Never
mix. For multi-line commit messages, prefer writing the message to a file under
`M:\claud_projects\temp` and using `git commit -F <file>` — it sidesteps shell
quoting entirely and the message is reviewable before it is committed.
Verify with `git log -1 --format='%s'` **before** pushing, not after.
