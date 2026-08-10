---
name: feedback-crlf-join-inflates-line-count
description: "Splitting a CRLF file on \\n and rejoining leaves a lone \\r mid-line, which splitlines() counts as a break — so a line-count reduction goes UP"
metadata:
  type: feedback
---

When compressing a **line-count-capped** file (the session-boot memory, `.claude/hooks/check_memory_size.py`),
never join lines from `text.split("\n")` without stripping `\r` first. On a CRLF file every element
still ends in `\r`, so `a + " " + b` embeds a **lone `\r` mid-line** — and `str.splitlines()`, which the
hook uses, treats a bare `\r` as a line break. The join that was supposed to *remove* a line **adds** one.

**Measured (D-175 batch-end):** joining two bullets took the file from 302 → 302 by the hook while
`count("\n")+1` said 301. Two different counters disagreeing by exactly the number of joins is the
signature. Fix: `t.replace("\r\n", "\n").replace("\r", " ")` before measuring, and measure with
`len(t.splitlines())` — the same function the checker uses — never `count("\n")+1`.

**Why:** the repo is `text=auto eol=lf` (`.gitattributes`) so the *committed* copy is LF, but the
working copy under Windows can be CRLF. A mixed file therefore produces **no git diff** for the
line endings, so the corruption is invisible in review — you only see it as a count that will not
go down. Distinct from [[feedback-no-powershell-heredoc-in-bash]], which is about encoding/BOM, not
line-count inflation.

**How to apply:** at batch end, when the memory hook says "file is N lines (backstop 300)", normalize
line endings first, then compress, then re-measure with `splitlines()`. If your count and the hook's
count differ at all, stop and find the stray separator before editing further —
[[feedback-count-and-print-your-skips]] is the same discipline applied to a denominator.
