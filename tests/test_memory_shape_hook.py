"""The session-boot memory shape hook (``.claude/hooks/check_memory_size.py``).

The hook warns when project memory drifts back into a changelog. It was tuned four
times (150 -> 200 -> 250 -> 300) before the git history showed the total line count
was being *written to* rather than bounded -- 250 exactly for thirteen commits across
twelve days, then +47 lines the day the cap moved. The pressure therefore moved to
per-block shape, which cannot be satisfied by evicting an older guardrail.

**These tests pin the hook's LOGIC on synthetic input, never the live files'
compliance.** A test asserting the real memory is under cap would convert a warning
into enforcement, and the hook's whole stated premise is that distilling a status
block is a judgement call no line count can make. The live files are allowed to be
over cap; the hook is required to *say so*.

The hook lives in ``.claude/hooks/`` (not an installed package), so it is loaded by
path, the same way ``test_decisions_index.py`` loads the TOC generator.
"""

from __future__ import annotations

import importlib.util
import io
import json
import pathlib
import sys
from types import ModuleType

import pytest

HOOK = pathlib.Path(__file__).resolve().parent.parent / ".claude" / "hooks" / "check_memory_size.py"


def _load_hook() -> ModuleType:
    spec = importlib.util.spec_from_file_location("check_memory_size", HOOK)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def hook() -> ModuleType:
    return _load_hook()


def _bullet(lines: int, text: str = "prohibition") -> str:
    """A top-level bullet occupying ``lines`` physical lines (continuations indent 2)."""
    return "\n".join([f"- **{text}**"] + ["  continuation" for _ in range(lines - 1)])


# --------------------------------------------------------------------------- shape


def test_bullet_at_cap_is_silent_and_one_over_is_not(hook: ModuleType) -> None:
    """The cap is inclusive: it is a ceiling on the block, not on the block minus one."""
    assert hook.project_findings(_bullet(hook.BLOCK_LINE_CAP)) == []
    over = hook.project_findings(_bullet(hook.BLOCK_LINE_CAP + 1))
    assert len(over) == 1
    assert f"{hook.BLOCK_LINE_CAP + 1} lines" in over[0].detail


def test_many_short_bullets_are_healthy_shape(hook: ModuleType) -> None:
    """The point of the metric change: accumulating guardrails must not trip it.

    Forty 4-line prohibitions is what the file is *for*. Only the backstop may
    speak here, and only about the total -- never about a block.
    """
    text = "\n\n".join(_bullet(4) for _ in range(40))
    assert [f for f in hook.project_findings(text) if "block" in f.detail] == []


def test_one_narrative_bullet_trips_it_at_a_fraction_of_the_total(hook: ModuleType) -> None:
    """A 21-line record (D-168's measured footprint) is caught in a 60-line file.

    This is the case the total line count could not see: well under LINE_CAP, and
    still the changelog shape.
    """
    text = "\n\n".join([_bullet(21, "D-168 narrative"), _bullet(3), _bullet(3)])
    assert len(text.splitlines()) < hook.LINE_CAP
    findings = hook.project_findings(text)
    assert len(findings) == 1
    assert findings[0].line == 1


def test_finding_names_the_records_to_distil(hook: ModuleType) -> None:
    """The usability win over the previous generic prose: which record, and where."""
    body = "- **D-168 and D-24**\n" + "\n".join("  cites D-163" for _ in range(20))
    (finding,) = hook.project_findings(body)
    assert finding.detail.endswith("-- D-24, D-163, D-168")  # numeric order, deduped
    assert finding.file == hook.PROJECT_NAME


def test_bullets_are_delimited_by_blank_lines_and_lead_ins(hook: ModuleType) -> None:
    """Two 5-line bullets under an axis heading are three blocks, not one."""
    text = (
        "**Axis (D-1)**\n" + _bullet(5) + "\n\n" + _bullet(5) + "\n\n**Next axis**\n" + _bullet(5)
    )
    assert [len(body) for _, body in hook.block_spans(text.splitlines())] == [1, 5, 5, 1, 5]
    assert hook.project_findings(text) == []


def test_a_column_zero_paragraph_cannot_escape_the_cap(hook: ModuleType) -> None:
    """The loophole a bullet-only cap would leave open.

    Narrative evicted from a bullet can be re-typed as an axis lead-in at column 0.
    That is the same arbitrage that moved a status paragraph into the uncapped
    ``MEMORY.md`` index row, one level in, so blocks are capped, not bullets.
    """
    paragraph = "\n".join(f"**Axis** narrative line {i} citing D-168" for i in range(12))
    (finding,) = hook.project_findings(paragraph)
    assert "block is 12 lines" in finding.detail


def test_yaml_frontmatter_is_exempt(hook: ModuleType) -> None:
    """It is schema, not content -- and the only live block over cap."""
    front = "---\nname: p\ndescription: d\nmetadata:\n" + "\n".join(
        f"  key{i}: v" for i in range(9)
    )
    assert hook.project_findings(front + "\n---\n\n" + _bullet(3)) == []


# ----------------------------------------------------------------------- backstop


def test_total_backstop_still_fires_on_a_file_of_short_bullets(hook: ModuleType) -> None:
    """Shape is the primary check, but the 2699-line regression must stay catchable."""
    text = "\n\n".join(_bullet(2) for _ in range(hook.LINE_CAP))
    totals = [f for f in hook.project_findings(text) if "backstop" in f.detail]
    assert len(totals) == 1


def test_backstop_is_not_the_binding_constraint_for_a_distilled_file(hook: ModuleType) -> None:
    """LINE_CAP was left at 300 on purpose -- if it bites, that is evidence to keep."""
    assert hook.project_findings("\n\n".join(_bullet(4) for _ in range(60))) == []


# -------------------------------------------------------------------------- index


def test_index_row_cap_catches_a_displaced_status_paragraph(hook: ModuleType) -> None:
    """The 950-char row: overflow squeezed out of the capped file into the uncapped one."""
    short = "- [Title](f.md) — hook"
    long = "- [Project](p.md) — " + "x" * hook.INDEX_ROW_CHAR_CAP
    findings = hook.index_findings(f"# Index\n\n{short}\n{long}\n")
    assert len(findings) == 1
    assert findings[0].line == 4
    assert findings[0].file == hook.INDEX_NAME


def test_index_cap_ignores_prose_that_is_not_an_index_row(hook: ModuleType) -> None:
    assert hook.index_findings("a very long paragraph " * 40) == []


# ------------------------------------------------------------------ payload wiring


def _run(hook: ModuleType, monkeypatch: pytest.MonkeyPatch, payload: dict[str, object]) -> str:
    monkeypatch.setattr(sys, "stdin", io.StringIO(json.dumps(payload)))
    out = io.StringIO()
    monkeypatch.setattr(sys, "stdout", out)
    assert hook.main() == 0
    return out.getvalue()


@pytest.fixture
def memory_dir(tmp_path: pathlib.Path, hook: ModuleType) -> pathlib.Path:
    (tmp_path / hook.PROJECT_NAME).write_text(_bullet(21), encoding="utf-8")
    (tmp_path / hook.INDEX_NAME).write_text(
        "- [P](p.md) — " + "x" * hook.INDEX_ROW_CHAR_CAP, encoding="utf-8"
    )
    return tmp_path


@pytest.mark.parametrize("edited", ["PROJECT_NAME", "INDEX_NAME"])
def test_editing_either_boot_file_reports_on_both(
    hook: ModuleType, memory_dir: pathlib.Path, monkeypatch: pytest.MonkeyPatch, edited: str
) -> None:
    """Checking only the edited file is what let the displacement go unseen."""
    target = memory_dir / getattr(hook, edited)
    emitted = json.loads(_run(hook, monkeypatch, {"tool_input": {"file_path": str(target)}}))
    message = emitted["systemMessage"]
    assert "2 finding(s)" in message
    assert hook.PROJECT_NAME in message and hook.INDEX_NAME in message
    assert emitted["hookSpecificOutput"]["hookEventName"] == "PostToolUse"


def test_windows_backslash_paths_are_matched(
    hook: ModuleType, memory_dir: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    raw = str(memory_dir / hook.PROJECT_NAME).replace("/", "\\")
    assert _run(hook, monkeypatch, {"tool_response": {"filePath": raw}}) != ""


@pytest.mark.parametrize(
    "payload",
    [
        {},
        {"tool_input": {"file_path": "src/fermentation/core/kinetics.py"}},
        {"tool_input": {"file_path": "other-memory.md"}},
    ],
    ids=["no-path", "unrelated-source-file", "another-memory-file"],
)
def test_silent_on_everything_else(
    hook: ModuleType, monkeypatch: pytest.MonkeyPatch, payload: dict[str, object]
) -> None:
    assert _run(hook, monkeypatch, payload) == ""


def test_malformed_payload_never_breaks_the_edit(
    hook: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A PostToolUse hook that raises would make memory edits fail; it must not."""
    monkeypatch.setattr(sys, "stdin", io.StringIO("not json"))
    assert hook.main() == 0


def test_report_is_capped_and_says_how_many_it_elided(
    hook: ModuleType, tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    count = hook.MAX_REPORTED + 5
    (tmp_path / hook.PROJECT_NAME).write_text(
        "\n\n".join(_bullet(9) for _ in range(count)), encoding="utf-8"
    )
    target = tmp_path / hook.PROJECT_NAME
    message = json.loads(_run(hook, monkeypatch, {"tool_input": {"file_path": str(target)}}))[
        "systemMessage"
    ]
    assert f"{count} finding(s)" in message
    assert message.count("block is") == hook.MAX_REPORTED
    assert f"and {count - hook.MAX_REPORTED} more" in message
