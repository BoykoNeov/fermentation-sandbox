#!/usr/bin/env python
"""Warn when the project memory file grows back into a changelog.

The project memory is session-boot context; the per-decision narrative belongs
in docs/DECISIONS.md alone (see .claude/memory/feedback-batch-end-ritual.md).

That rule already existed as prose when commit acd3ce1 (2026-07-02) cut the file
29KB -> 2.4KB and added the guardrail "to fix the cause, not just the symptom".
It regrew to 277KB / 2699 lines in 15 days -- 114x. A prose rule in a memory file
is not a mechanism, so this is the mechanism: a PostToolUse hook that makes the
regression visible at the moment it is written.

It warns; it cannot enforce. Distilling a status block is a judgement call, and
no line count can make it. The point is that the drift stops being silent.

Reads the PostToolUse payload on stdin; emits hook JSON on stdout when over cap.
"""

from __future__ import annotations

import json
import pathlib
import sys

# The distilled file is ~60 lines. 200 leaves generous headroom for a real status
# update while still catching the changelog shape (the bloated copies ran 550+).
# Raised 150 -> 200 on 2026-07-21 at the owner's request (the live-threads block grew
# with the D-130/D-131 aging builds; still well below the changelog regime).
LINE_CAP = 200
TARGET_NAME = "project-fermentation-sandbox.md"


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
    if path.name != TARGET_NAME:
        return 0

    try:
        lines = len(path.read_text(encoding="utf-8").splitlines())
    except OSError:
        return 0

    if lines <= LINE_CAP:
        return 0

    message = (
        f"MEMORY BLOAT: {TARGET_NAME} is now {lines} lines (cap {LINE_CAP}). "
        "It is session-boot context, NOT a changelog -- the full per-decision "
        "narrative belongs in docs/DECISIONS.md only. This file regrew 114x after "
        "the 2026-07-02 re-cut because the guardrail was prose. Distil it back to a "
        "status block + pointers before committing (see feedback-batch-end-ritual)."
    )
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
