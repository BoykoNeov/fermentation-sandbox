---
name: feedback-no-powershell-heredoc-in-bash
description: "Commit messages get corrupted by shell mechanics: PS here-strings in the Bash tool, Out-File's BOM, and embedded double quotes splitting the arg to native git — write the message to a file and use -F"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 4c5542b6-994d-42ff-9b5f-a6dbc7d14d50
  modified: 2026-07-28T11:40:23.359Z
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

**A second corruption in the same family, found at D-135 (2026-07-22): the BOM.**
`... | Out-File -FilePath x.txt -Encoding utf8` in **Windows PowerShell 5.1
writes a UTF-8 BOM**, so `git commit -F x.txt` produced a subject line starting
with an invisible `U+FEFF` (`efbbbf` under `xxd`). `git log --oneline` looked
almost right; the Conventional Commits prefix was no longer at the start of the
line. Exit 0 again, and it was **already pushed** before I noticed — so fixing it
needed the force-push this memory was written to avoid. **Use
`-Encoding utf8NoBOM` (PS 7+) or `[IO.File]::WriteAllText(path, text)` (any
version), or avoid the file entirely and pass `-m` a here-string directly.**

**Also at D-135: a large bash heredoc in the Bash tool failed to parse** —
`cat >> file <<'EOF'` with a ~160-line markdown payload died with
`unexpected EOF while looking for matching '`, despite the quoted delimiter.
Nothing was written (verified), so it failed loudly rather than silently. For
big payloads, skip the shell: use the **Write tool** to create the file, then
append with a short `python -c` — no quoting layer to get wrong.

**RECURRED at D-138 (2026-07-27), identically — this memory existed and did not
stop it.** Same `git commit -m @'...'@` in the Bash tool, same stray `@` subject,
same force-push to fix. The reason it did not stop it: the failure happens at
*compose* time, but the memory's guard fires at *verify* time, and I ran the
verify **after `git push`** — in one `&&` chain that committed and pushed
together. **The chain is the bug.** `git commit && git push` in a single command
removes the only window in which the check is worth anything.

**A THIRD mechanic, 2026-07-28: embedded double quotes split the argument.**
`git commit -m @'...'@` in the **PowerShell** tool — closing `'@` at column 0,
exactly as this memory prescribed — still failed. PS 5.1 does **not** escape a
`"` inside the string when handing it to a native exe, so `git`'s own
command-line parser re-splits the argument at the quote. A message containing
`D-72's "right-looking 1:2"` reached git as several args and died with
`fatal: Invalid path '1:/2 trap was folded back in...'`. **Nothing was
committed** — and I only noticed because the fragment happened to resemble a
path. Had the tail parsed cleanly, git would have committed a **truncated
message at exit 0**, the silent failure this whole memory is about. So the
`@'...'@` advice above is **necessary but not sufficient**: it protects `$` and
backticks, not `"`.

**How to apply:** for any commit message containing a double quote — or just by
default for multi-paragraph messages — **write it with the Write tool and use
`git commit -F <file>`**. That has now survived where `-m` failed twice. Never
put `git commit` and `git push` in the same command —
commit, verify `%s`, then push as a separate call. In the Bash tool use a real
bash heredoc —
`git commit -F - <<'EOF' ... EOF` (quoted `EOF` keeps `$` and backticks literal)
— for *short* messages only. In the PowerShell tool use `@'...'@` with the
closing `'@` at column 0. Never mix. If routing through a file, write it with
the Write tool, not `Out-File`. **Verify before pushing** — a visual check does not catch a
BOM, and the whole point is that every failure in this family exits 0.

**But verify with a BINARY read, never through a Git Bash pipe.** `git cat-file commit HEAD |
od -c` (and `git log --format=%B | xxd`) reported **118 CR bytes** in a message that has
**zero** — MSYS translates LF → CRLF *inside the pipe*, so the instrument invents exactly the
corruption it is being used to rule out, and it does so for **every** commit in the repo, which
makes it look like a confirmed repo-wide defect rather than a measurement artefact. Read the
object in Python instead: `subprocess.run(["git","cat-file","commit","HEAD"],
capture_output=True).stdout` then `.count(b"\r")` and check `[:3] != b"\xef\xbb\xbf"`. A check
whose false-positive rate is 100 % is worse than no check — it would have sent me to "fix"
clean history [[feedback-name-the-field-your-predicate-read]].
