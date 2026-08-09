#!/usr/bin/env python
"""Warn when the session-boot memory drifts back into a changelog.

The project memory is session-boot context; the per-decision narrative belongs
in docs/DECISIONS.md alone (see .claude/memory/feedback-batch-end-ritual.md).

That rule already existed as prose when commit acd3ce1 (2026-07-02) cut the file
29KB -> 2.4KB and added the guardrail "to fix the cause, not just the symptom".
It regrew to 277KB / 2699 lines in 15 days -- 114x. A prose rule in a memory file
is not a mechanism, so this is the mechanism: a PostToolUse hook that makes the
regression visible at the moment it is written.

It warns; it cannot enforce. Distilling a status block is a judgement call, and
no line count can make it. The point is that the drift stops being silent.

WHY THE TOTAL LINE COUNT IS NO LONGER THE PRIMARY CHECK
-------------------------------------------------------
LINE_CAP was raised 150 -> 200 -> 250 -> 300, each time at the owner's request and
each time for the same reported reason: the cap was evicting live prohibitions to
make room for newer ones. The git history of the target file shows why raising it
never worked -- it was being written to, not merely bounded:

    2026-07-22  178          2026-07-28  250  |
    2026-07-28  211 .. 237   2026-07-29  250  | 13 commits, 12 days,
                             ...        250  | pinned EXACTLY at 250
                             2026-08-09  250  |
    ------------------- cap raised 250 -> 300 -------------------
    2026-08-09  252 -> 267 -> 285 -> 297   (+47 lines in ONE day)

Content does not land on a round number thirteen times. The cap was functioning as
a *target*: the fill rate was set by the cap, not by the work. And because a total
was the ONLY check, the sole way to satisfy it was eviction -- which is exactly the
harm the raises were meant to undo.

So the pressure moves to SHAPE. Guardrails accumulate slowly and are 1-4 lines each;
a changelog entry is 20-30. A file of 300 lines made of 75 short prohibitions is
healthy; the same 300 lines made of 12 narratives is the disease, and a total is
blind to the difference. BULLET_LINE_CAP sees it exactly, and cannot be satisfied by
evicting an old record -- only by distilling the new one.

LINE_CAP stays at 300 deliberately. It is now a backstop behind the per-bullet check,
not the binding constraint, so it should stop being the thing that bites. If it bites
anyway, that is evidence worth having BEFORE it is raised a fifth time.

A digit-density check was designed and REJECTED: digits-per-line is ~10 for the
21-line narrative bullet and ~20 for a 3-line one, so an absolute digit budget reads
bullet SIZE, not evidence density -- and the highest-value guardrails are the ones
that are nothing but corrected values ("never re-narrow 0.084 to 0.08", "247 not
279"). It would have penalised precisely what must be retained.

NOT COVERED: CLAUDE.md is the third session-boot surface and is unmeasured. Capping
project memory alone already displaced a status paragraph into the (then unchecked)
MEMORY.md index row -- 950 chars against a 211-char median. The same arbitrage exists
one level out.

Reads the PostToolUse payload on stdin; emits hook JSON on stdout when over cap.
"""

from __future__ import annotations

import json
import pathlib
import re
import sys
from typing import NamedTuple

PROJECT_NAME = "project-fermentation-sandbox.md"
INDEX_NAME = "MEMORY.md"

# Whole-file backstop. See the module docstring: no longer the primary check, and
# not to be raised without first showing that the per-bullet check did not bite.
LINE_CAP = 300

# One top-level block == one distilled record. Measured 2026-08-09 across the 49
# bullets then live: median 4 lines, 44 of 49 at or under 8. The five over were the
# five most recent beats (21, 14, 10, 10, 9) -- the recency skew that makes each new
# record eat the budget an older guardrail was holding.
#
# It caps BLOCKS, not just "- " bullets: a cap that bound only bullets would let
# narrative displace into a column-0 paragraph and escape, which is the same
# arbitrage as INDEX_ROW_CHAR_CAP one level in. The 20 non-bullet paragraphs then
# live measured 6 lines at most, so this costs prose nothing.
BLOCK_LINE_CAP = 8

# One "- [Title](file.md) -- hook" row in the index. Measured 2026-08-09: median 211
# chars, next-longest 308, and the project-memory row at 950 -- a status paragraph
# that had been displaced out of the capped file into the uncapped one.
INDEX_ROW_CHAR_CAP = 320

MAX_REPORTED = 10


class Finding(NamedTuple):
    file: str
    line: int
    detail: str


def block_spans(lines: list[str]) -> list[tuple[int, list[str]]]:
    """Group the body into top-level blocks: bullets and paragraphs alike.

    A block runs until a blank line, and a "- " always starts a new one, so a bullet
    with wrapped continuations and a heading with its lead-in prose each count once.
    The YAML frontmatter is skipped -- it is schema, not content, and is the only
    block in the live file that exceeds the cap.
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


def project_findings(text: str) -> list[Finding]:
    lines = text.splitlines()
    findings = [
        Finding(
            PROJECT_NAME,
            start,
            f"block is {len(body)} lines (cap {BLOCK_LINE_CAP}){_cited(body)}",
        )
        for start, body in block_spans(lines)
        if len(body) > BLOCK_LINE_CAP
    ]
    if len(lines) > LINE_CAP:
        findings.append(
            Finding(PROJECT_NAME, len(lines), f"file is {len(lines)} lines (backstop {LINE_CAP})")
        )
    return findings


def index_findings(text: str) -> list[Finding]:
    return [
        Finding(INDEX_NAME, number, f"index row is {len(line)} chars (cap {INDEX_ROW_CHAR_CAP})")
        for number, line in enumerate(text.splitlines(), 1)
        if line.startswith("- [") and len(line) > INDEX_ROW_CHAR_CAP
    ]


def collect(directory: pathlib.Path) -> list[Finding]:
    findings: list[Finding] = []
    for name, check in ((PROJECT_NAME, project_findings), (INDEX_NAME, index_findings)):
        try:
            text = (directory / name).read_text(encoding="utf-8")
        except OSError:
            continue
        findings.extend(check(text))
    return findings


def render(findings: list[Finding]) -> str:
    shown = findings[:MAX_REPORTED]
    body = "\n".join(f"  {f.file}:{f.line}  {f.detail}" for f in shown)
    elided = len(findings) - len(shown)
    if elided:
        body += f"\n  ... and {elided} more"
    return (
        f"MEMORY SHAPE: {len(findings)} finding(s) in session-boot context.\n"
        f"{body}\n"
        "Session-boot memory is PROHIBITIONS + POINTERS, not a changelog. Distil each "
        "block to trigger + verdict + anchor; the measurements it cites are already "
        "held losslessly in docs/DECISIONS.md and are reachable by grep. Numbers that "
        "ARE the prohibition stay. Do NOT satisfy this by deleting an older guardrail "
        "-- the per-block cap exists so the newest record distils instead."
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

    # Windows paths arrive with backslashes; the hook must match either form.
    path = pathlib.Path(str(raw).replace("\\", "/"))
    if path.name not in (PROJECT_NAME, INDEX_NAME):
        return 0

    # Both files load at boot and overflow moves between them, so an edit to either
    # reports on both -- checking only the edited one is what let the displacement
    # into the index row go unseen.
    findings = collect(path.parent)
    if not findings:
        return 0

    message = render(findings)
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
