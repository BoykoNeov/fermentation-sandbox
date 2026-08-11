#!/usr/bin/env python
"""Keep the three session-boot surfaces visible, and warn when one drifts into a changelog.

The boot surfaces are `CLAUDE.md`, `.claude/memory/MEMORY.md` and
`.claude/memory/project-fermentation-sandbox.md`. The per-decision narrative belongs in
docs/DECISIONS.md alone (see .claude/memory/feedback-batch-end-ritual.md).

That rule already existed as prose when commit acd3ce1 (2026-07-02) cut the project memory
29KB -> 2.4KB and added the guardrail "to fix the cause, not just the symptom". It regrew to
277KB / 2699 lines in 15 days -- 114x. A prose rule in a memory file is not a mechanism, so
this is the mechanism: a PostToolUse hook that makes the regression visible at the moment it
is written.

It warns; it cannot enforce. Distilling a status block is a judgement call, and no line count
can make it. The point is that the drift stops being silent.

WHY THERE IS NO WHOLE-FILE CAP ANY MORE (2026-08-11)
----------------------------------------------------
LINE_CAP was raised 150 -> 200 -> 250 -> 300, each time at the owner's request and each time
for the same reported reason: the cap was evicting live prohibitions to make room for newer
ones. D-169 measured why raising never worked -- the file was being written *to* the cap, not
bounded by it (pinned at exactly 250 for 13 commits, then +47 lines in one day when the cap
moved). D-169 kept 300 as a backstop behind the new per-block check and pre-committed to a
test: "If it bites anyway, that is evidence worth having BEFORE it is raised a fifth time."

It bit. Measured 2026-08-11, the same two ways:

    SHAPE  77 blocks, median 3 lines, max 8 -- ZERO over BLOCK_LINE_CAP. This is verbatim
           the healthy file this docstring described ("300 lines made of 75 short
           prohibitions is healthy").
    TOTAL  exactly 300 on 8 of the last 9 commits (one at 299) -- the identical
           round-number signature that diagnosed 250.

So the per-block check works and the total is a target again, at 300. The arithmetic is why,
and it is not a discipline failure: BLOCK_LINE_CAP bounds what each new record ADDS (<= 8
lines); nothing bounds how many records there are. At ~1 record per session and ~3.5 lines
each, any fixed total is reached every ~10 sessions, at any value.

D-169's licensed relief was retirement, not a raise: a record superseded per the correction
map collapses to a pointer. That was sized before acting. Of the 84 D-records the project
memory cites, 39 carry a correction/flag marker; the 18 blocks citing ONLY corrected records
occupy 59 of 300 lines -- but reading them, nearly all state prohibitions their corrector
SHARPENED rather than replaced ("Mechanism B CLOSED", "a pH term is REFUSED", "600 L/g, never
re-open"). Those blocks are already collapsed. Realistic recovery is 5-15 lines: one to four
sessions. The escape hatch cannot pay for the check.

Two replacements were designed and rejected before landing on the third:

  * A DERIVED total (e.g. 8 x blocks) has no round number to write toward -- but the median
    block is 3 lines, so the file runs at ~3.8x blocks and the check can never fire. A guard
    that forbids nothing is worse than none: it reads as coverage. Cf. D-172's sum assert
    becoming an identity.
  * A FIFTH RAISE is what the whole record forbids.

So the total is now REPORTED, never capped. A number with no threshold cannot be a target
(there is nothing to hit) and cannot be vacuous (it always shows). Accumulation stays visible
on every write, and calling a retirement pass goes back to being the owner's judgement --
which is what the standing rule wanted, minus the round number to fill.

THE THIRD SURFACE IS NOW MEASURED -- AND WHAT THAT CHECK DOES NOT CATCH
-----------------------------------------------------------------------
The previous version of this docstring named its own gap: "NOT COVERED: CLAUDE.md is the
third session-boot surface and is unmeasured." Measured 2026-08-11: it has grown 66 -> 138
lines (3.3KB -> 8.0KB) since 2026-06-20, and took +30 lines on 2026-08-09 -- the same day the
memory cap moved to 300. What landed is prohibition-shaped, not convention-shaped: a
measurement narrative on suite parallelism plus "Do not switch to --dist worksteal", with its
reasoning attached. One topic is now ~25% of the file.

That is NOT proof of displacement -- the perf session that wrote it had its own reason, and
the content is legitimately about this project's test suite. What is proven is that the third
surface holds the same content type, grows at a comparable rate, and nothing looked at it.

GUIDE_BLOCK_LINE_CAP is sized from CLAUDE.md's OWN distribution (34 blocks, median 3, max 12
-- the prime-directives list item at 12 and the worksteal note at 11), not borrowed from the
memory file's 8, which would fire on legitimate documentation prose and get argued away.
14 clears both live blocks by 2 and still trips the 20-30 line changelog shape that is the
disease.

IT WOULD NOT HAVE CAUGHT 2026-08-09. That day arrived as two ~10-line blocks, both under any
honest cap for this file. Catching it needs a GROWTH check (+28% in a day), which needs git
history a PostToolUse hook does not have. Not shipped; named here so the gap is not conceded
in prose and then treated as covered.

A digit-density check was designed and REJECTED: digits-per-line is ~10 for the 21-line
narrative bullet and ~20 for a 3-line one, so an absolute digit budget reads bullet SIZE, not
evidence density -- and the highest-value guardrails are the ones that are nothing but
corrected values ("never re-narrow 0.084 to 0.08", "247 not 279"). It would have penalised
precisely what must be retained.

STILL NOT COVERED: the global ~/.claude/CLAUDE.md is a fourth boot surface, lives outside the
repo, and is deliberately out of scope here. MEMORY.md is capped per ROW, so row COUNT
remains an open channel -- 5 rows at 2026-06-23, 40 at 2026-08-11, +1 per record.

Reads the PostToolUse payload on stdin; emits hook JSON on stdout.
"""

from __future__ import annotations

import json
import pathlib
import re
import sys
from typing import NamedTuple

PROJECT_NAME = "project-fermentation-sandbox.md"
INDEX_NAME = "MEMORY.md"
GUIDE_NAME = "CLAUDE.md"

MEMORY_DIR = (".claude", "memory")

# One top-level block == one distilled record. Measured 2026-08-09 across the 49 bullets then
# live: median 4 lines, 44 of 49 at or under 8. The five over were the five most recent beats
# (21, 14, 10, 10, 9) -- the recency skew that makes each new record eat the budget an older
# guardrail was holding.
#
# It caps BLOCKS, not just "- " bullets: a cap that bound only bullets would let narrative
# displace into a column-0 paragraph and escape, which is the same arbitrage as
# INDEX_ROW_CHAR_CAP one level in. The 20 non-bullet paragraphs then live measured 6 lines at
# most, so this costs prose nothing.
#
# This is now the ONLY threshold on the project memory. See the docstring for why the total
# went away rather than up.
BLOCK_LINE_CAP = 8

# CLAUDE.md is documentation, not distilled records, so it gets its own number rather than
# the memory file's. Sized from its own shape (measured 2026-08-11: 34 blocks, median 3,
# max 12) to be inert today and trip the changelog shape. See the docstring for what it
# does NOT catch.
GUIDE_BLOCK_LINE_CAP = 14

# One "- [Title](file.md) -- hook" row in the index. Measured 2026-08-09: median 211 chars,
# next-longest 308, and the project-memory row at 950 -- a status paragraph that had been
# displaced out of the capped file into the uncapped one.
INDEX_ROW_CHAR_CAP = 320

MAX_REPORTED = 10


class Finding(NamedTuple):
    file: str
    line: int
    detail: str


def block_spans(lines: list[str]) -> list[tuple[int, list[str]]]:
    """Group the body into top-level blocks: bullets and paragraphs alike.

    A block runs until a blank line, and a "- " always starts a new one, so a bullet with
    wrapped continuations and a heading with its lead-in prose each count once. YAML
    frontmatter is skipped -- it is schema, not content, and is the only block in the live
    memory file that exceeds the cap.
    """
    start_at = 0
    if lines and lines[0].strip() == "---":
        closing = next((i for i, line in enumerate(lines[1:], 1) if line.strip() == "---"), None)
        if closing is not None:
            start_at = closing + 1

    spans: list[tuple[int, list[str]]] = []
    current: tuple[int, list[str]] | None = None
    for number, line in enumerate(lines[start_at:], start_at + 1):
        if not line.strip():
            if current is not None:
                spans.append(current)
                current = None
        elif current is None or line.startswith("- "):
            if current is not None:
                spans.append(current)
            current = (number, [line])
        else:
            current[1].append(line)
    if current is not None:
        spans.append(current)
    return spans


def _cited(body: list[str]) -> str:
    """The records a block cites, so the warning names what to go and distil."""
    records = sorted(set(re.findall(r"D-\d+", " ".join(body))), key=lambda r: int(r[2:]))
    return f" -- {', '.join(records)}" if records else ""


def block_findings(name: str, text: str, cap: int) -> list[Finding]:
    return [
        Finding(name, start, f"block is {len(body)} lines (cap {cap}){_cited(body)}")
        for start, body in block_spans(text.splitlines())
        if len(body) > cap
    ]


def index_findings(text: str) -> list[Finding]:
    return [
        Finding(INDEX_NAME, number, f"index row is {len(line)} chars (cap {INDEX_ROW_CHAR_CAP})")
        for number, line in enumerate(text.splitlines(), 1)
        if line.startswith("- [") and len(line) > INDEX_ROW_CHAR_CAP
    ]


def _shape(text: str) -> str:
    """The reported totals: size with no threshold attached, so it cannot become a target."""
    lines = text.splitlines()
    sizes = sorted(len(body) for _, body in block_spans(lines))
    if not sizes:
        return f"{len(lines)} lines"
    median = sizes[len(sizes) // 2]
    return f"{len(lines)} lines, {len(sizes)} blocks, median {median}, max {sizes[-1]}"


def project_root(path: pathlib.Path) -> pathlib.Path | None:
    """Walk up to the directory holding both CLAUDE.md and the memory dir.

    This is what keeps the global ~/.claude/CLAUDE.md out of scope: nothing on its parent
    chain carries a memory dir, so it resolves to None and the hook stays silent.
    """
    for candidate in (path, *path.parents):
        if (candidate / GUIDE_NAME).is_file() and (candidate.joinpath(*MEMORY_DIR)).is_dir():
            return candidate
    return None


def collect(root: pathlib.Path) -> tuple[list[Finding], list[str]]:
    """Findings that need action, and the reported shape of every surface."""
    memory = root.joinpath(*MEMORY_DIR)
    surfaces = (
        (GUIDE_NAME, root / GUIDE_NAME, GUIDE_BLOCK_LINE_CAP),
        (PROJECT_NAME, memory / PROJECT_NAME, BLOCK_LINE_CAP),
        (INDEX_NAME, memory / INDEX_NAME, None),
    )

    findings: list[Finding] = []
    reported: list[str] = []
    for name, path, cap in surfaces:
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        findings.extend(block_findings(name, text, cap) if cap else index_findings(text))
        reported.append(f"  {name}: {_shape(text)}")
    return findings, reported


def render(findings: list[Finding], reported: list[str]) -> str:
    body = "\n".join(reported)
    if not findings:
        return f"MEMORY SHAPE: clean. Session-boot surfaces:\n{body}"

    shown = findings[:MAX_REPORTED]
    detail = "\n".join(f"  {f.file}:{f.line}  {f.detail}" for f in shown)
    elided = len(findings) - len(shown)
    if elided:
        detail += f"\n  ... and {elided} more"
    return (
        f"MEMORY SHAPE: {len(findings)} finding(s) in session-boot context.\n"
        f"{detail}\n"
        f"Surfaces:\n{body}\n"
        "Session-boot memory is PROHIBITIONS + POINTERS, not a changelog. Distil each block "
        "to trigger + verdict + anchor; the measurements it cites are already held losslessly "
        "in docs/DECISIONS.md and are reachable by grep. Numbers that ARE the prohibition "
        "stay. Do NOT satisfy this by deleting an older guardrail -- the per-block cap exists "
        "so the newest record distils instead. There is no whole-file cap: the totals above "
        "are reported so accumulation stays visible, and a retirement pass is a judgement "
        "call, never a number to trim to."
    )


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    tool_input = payload.get("tool_input") or {}
    tool_response = payload.get("tool_response") or {}
    raw = tool_response.get("filePath") or tool_input.get("file_path") or ""
    if not raw:
        return 0

    # Windows paths arrive with backslashes; the memory dir is reached through a junction
    # from the harness path, so resolve to the one real location before matching.
    path = pathlib.Path(str(raw).replace("\\", "/"))
    if path.name not in (PROJECT_NAME, INDEX_NAME, GUIDE_NAME):
        return 0
    try:
        path = path.resolve()
    except OSError:
        return 0

    root = project_root(path.parent)
    if root is None:
        return 0

    # All three load at boot and overflow moves between them, so an edit to any one reports
    # on all three -- checking only the edited file is what let the displacement into the
    # index row, and then into CLAUDE.md, go unseen.
    findings, reported = collect(root)
    if not reported:
        return 0

    message = render(findings, reported)
    json.dump(
        {
            "systemMessage": message,
            "hookSpecificOutput": {
                "hookEventName": "PostToolUse",
                "additionalContext": message,
            },
        },
        sys.stdout,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
