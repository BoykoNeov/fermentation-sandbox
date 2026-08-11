"""The session-boot memory shape hook (``.claude/hooks/check_memory_size.py``).

The hook warns when a session-boot surface drifts back into a changelog. Its total line
count was tuned four times (150 -> 200 -> 250 -> 300) before the git history showed the
number was being *written to* rather than bounded -- 250 exactly for thirteen commits
across twelve days, then +47 lines the day the cap moved (D-169). The pressure moved to
per-block shape, which cannot be satisfied by evicting an older guardrail.

D-177 removed the whole-file total entirely rather than raising it a fifth time: at 300
the same round-number signature returned (exactly 300 on eight of nine commits) while
every block sat under cap, because a per-block cap bounds what each record ADDS and
nothing bounds how many records there are. Totals are now REPORTED, never capped -- a
number with no threshold cannot be a target -- and CLAUDE.md joined as the third measured
surface with a cap sized from its own shape.

**These tests pin the hook's LOGIC on synthetic input, never the live files'
compliance.** A test asserting the real memory is under cap would convert a warning into
enforcement, and the hook's whole stated premise is that distilling a status block is a
judgement call no line count can make. The live files are allowed to be over cap; the
hook is required to *say so*.

The hook lives in ``.claude/hooks/`` (not an installed package), so it is loaded by path,
the same way ``test_decisions_index.py`` loads the TOC generator.
"""

from __future__ import annotations

import importlib.util
import io
import json
import pathlib
import sys
from types import ModuleType
from typing import Any, cast

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


def _memory_findings(hook: ModuleType, text: str) -> list[Any]:
    """Block findings as the project memory file would produce them."""
    return cast("list[Any]", hook.block_findings(hook.PROJECT_NAME, text, hook.BLOCK_LINE_CAP))


# --------------------------------------------------------------------------- shape


def test_bullet_at_cap_is_silent_and_one_over_is_not(hook: ModuleType) -> None:
    """The cap is inclusive: it is a ceiling on the block, not on the block minus one."""
    assert _memory_findings(hook, _bullet(hook.BLOCK_LINE_CAP)) == []
    over = _memory_findings(hook, _bullet(hook.BLOCK_LINE_CAP + 1))
    assert len(over) == 1
    assert f"{hook.BLOCK_LINE_CAP + 1} lines" in over[0].detail


def test_many_short_bullets_are_healthy_shape(hook: ModuleType) -> None:
    """The point of the metric change: accumulating guardrails must not trip it.

    Forty 4-line prohibitions is what the file is *for*.
    """
    assert _memory_findings(hook, "\n\n".join(_bullet(4) for _ in range(40))) == []


def test_one_narrative_bullet_trips_it_at_a_fraction_of_any_total(hook: ModuleType) -> None:
    """A 21-line record (D-168's measured footprint) is caught in a 27-line file.

    This is the case a total line count could not see: far under every value LINE_CAP
    ever held, and still the changelog shape.
    """
    text = "\n\n".join([_bullet(21, "D-168 narrative"), _bullet(3), _bullet(3)])
    assert len(text.splitlines()) < 150  # the lowest LINE_CAP ever shipped
    findings = _memory_findings(hook, text)
    assert len(findings) == 1
    assert findings[0].line == 1


def test_finding_names_the_records_to_distil(hook: ModuleType) -> None:
    """The usability win over generic prose: which record, and where."""
    body = "- **D-168 and D-24**\n" + "\n".join("  cites D-163" for _ in range(20))
    (finding,) = _memory_findings(hook, body)
    assert finding.detail.endswith("-- D-24, D-163, D-168")  # numeric order, deduped
    assert finding.file == hook.PROJECT_NAME


def test_bullets_are_delimited_by_blank_lines_and_lead_ins(hook: ModuleType) -> None:
    """Two 5-line bullets under an axis heading are three blocks, not one."""
    text = (
        "**Axis (D-1)**\n" + _bullet(5) + "\n\n" + _bullet(5) + "\n\n**Next axis**\n" + _bullet(5)
    )
    assert [len(body) for _, body in hook.block_spans(text.splitlines())] == [1, 5, 5, 1, 5]
    assert _memory_findings(hook, text) == []


def test_a_column_zero_paragraph_cannot_escape_the_cap(hook: ModuleType) -> None:
    """The loophole a bullet-only cap would leave open.

    Narrative evicted from a bullet can be re-typed as an axis lead-in at column 0. That
    is the same arbitrage that moved a status paragraph into the uncapped ``MEMORY.md``
    index row, one level in, so blocks are capped, not bullets.
    """
    paragraph = "\n".join(f"**Axis** narrative line {i} citing D-168" for i in range(12))
    (finding,) = _memory_findings(hook, paragraph)
    assert "block is 12 lines" in finding.detail


def test_yaml_frontmatter_is_exempt(hook: ModuleType) -> None:
    """It is schema, not content -- and the only live block over cap."""
    front = "---\nname: p\ndescription: d\nmetadata:\n" + "\n".join(
        f"  key{i}: v" for i in range(9)
    )
    assert _memory_findings(hook, front + "\n---\n\n" + _bullet(3)) == []


# ------------------------------------------------------- the total is reported, not capped


def test_no_length_of_short_bullets_is_ever_a_finding(hook: ModuleType) -> None:
    """D-177: the whole-file total is gone, not raised.

    A file of 400 distilled prohibitions is healthy at any length -- that was already
    this hook's stated position ("300 lines made of 75 short prohibitions is healthy"),
    and keeping a total meant the file was trimmed to a round number instead.
    """
    text = "\n\n".join(_bullet(3) for _ in range(400))
    assert len(text.splitlines()) > 1000
    assert _memory_findings(hook, text) == []


def test_no_threshold_constant_survives_for_a_whole_file(hook: ModuleType) -> None:
    """The regression guard: re-adding a total is what D-177 forbids.

    Reintroducing ``LINE_CAP`` would restore the number the file was written to for
    thirteen commits at 250 and eight at 300.
    """
    assert not hasattr(hook, "LINE_CAP")


def test_shape_report_states_the_total_without_judging_it(hook: ModuleType) -> None:
    """Accumulation stays visible; the visibility is the replacement for the cap."""
    report = hook._shape("\n\n".join(_bullet(3) for _ in range(100)))  # + 99 separators
    assert "399 lines" in report and "100 blocks" in report and "median 3" in report
    assert "cap" not in report and "backstop" not in report


# ------------------------------------------------------------------- the third surface


def test_guide_carries_its_own_cap_not_the_memory_file_s(hook: ModuleType) -> None:
    """CLAUDE.md is documentation, so its cap comes from its own measured shape.

    Its live prime-directives item is 12 lines. Borrowing the memory file's 8 would fire
    on legitimate documentation prose and get argued away.
    """
    assert hook.GUIDE_BLOCK_LINE_CAP > hook.BLOCK_LINE_CAP
    twelve = _bullet(12, "Fidelity is tiered")
    assert hook.block_findings(hook.GUIDE_NAME, twelve, hook.GUIDE_BLOCK_LINE_CAP) == []
    assert _memory_findings(hook, twelve) != []  # the same block, in the memory file


def test_guide_cap_still_trips_the_changelog_shape(hook: ModuleType) -> None:
    """What it is for: the 20-30 line narrative entry, wherever it lands."""
    (finding,) = hook.block_findings(
        hook.GUIDE_NAME, _bullet(25, "D-168 narrative"), hook.GUIDE_BLOCK_LINE_CAP
    )
    assert finding.file == hook.GUIDE_NAME


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
def repo(tmp_path: pathlib.Path, hook: ModuleType) -> pathlib.Path:
    """A project root shaped like the real one: guide at the top, memory two levels in."""
    memory = tmp_path.joinpath(*hook.MEMORY_DIR)
    memory.mkdir(parents=True)
    (tmp_path / hook.GUIDE_NAME).write_text(_bullet(3), encoding="utf-8")
    (memory / hook.PROJECT_NAME).write_text(_bullet(21), encoding="utf-8")
    (memory / hook.INDEX_NAME).write_text(
        "- [P](p.md) — " + "x" * hook.INDEX_ROW_CHAR_CAP, encoding="utf-8"
    )
    return tmp_path


def _emit(
    hook: ModuleType, monkeypatch: pytest.MonkeyPatch, target: pathlib.Path
) -> dict[str, Any]:
    raw = _run(hook, monkeypatch, {"tool_input": {"file_path": str(target)}})
    return cast("dict[str, Any]", json.loads(raw))


@pytest.mark.parametrize("edited", ["PROJECT_NAME", "INDEX_NAME", "GUIDE_NAME"])
def test_editing_any_boot_file_reports_on_all_three(
    hook: ModuleType, repo: pathlib.Path, monkeypatch: pytest.MonkeyPatch, edited: str
) -> None:
    """Checking only the edited file is what let the displacement go unseen -- twice.

    First into the uncapped ``MEMORY.md`` row, then into the unmeasured ``CLAUDE.md``.
    """
    name = getattr(hook, edited)
    target = repo / name if name == hook.GUIDE_NAME else repo.joinpath(*hook.MEMORY_DIR, name)
    message = _emit(hook, monkeypatch, target)["systemMessage"]
    assert "2 finding(s)" in message
    assert all(n in message for n in (hook.GUIDE_NAME, hook.PROJECT_NAME, hook.INDEX_NAME))


def test_a_clean_surface_still_reports_its_shape(
    hook: ModuleType, repo: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Silence on clean is what let accumulation go unnoticed between the raises."""
    memory = repo.joinpath(*hook.MEMORY_DIR)
    (memory / hook.PROJECT_NAME).write_text(_bullet(4), encoding="utf-8")
    (memory / hook.INDEX_NAME).write_text("- [P](p.md) — hook", encoding="utf-8")
    message = _emit(hook, monkeypatch, memory / hook.PROJECT_NAME)["systemMessage"]
    assert "clean" in message
    assert "4 lines, 1 blocks" in message


def test_a_guide_outside_a_project_is_out_of_scope(
    hook: ModuleType, tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The global ``~/.claude/CLAUDE.md`` shares the name and is deliberately not measured."""
    stray = tmp_path / hook.GUIDE_NAME
    stray.write_text(_bullet(30), encoding="utf-8")
    assert _run(hook, monkeypatch, {"tool_input": {"file_path": str(stray)}}) == ""


def test_windows_backslash_paths_are_matched(
    hook: ModuleType, repo: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    raw = str(repo.joinpath(*hook.MEMORY_DIR, hook.PROJECT_NAME)).replace("/", "\\")
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
    hook: ModuleType, repo: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    count = hook.MAX_REPORTED + 5
    target = repo.joinpath(*hook.MEMORY_DIR, hook.PROJECT_NAME)
    target.write_text("\n\n".join(_bullet(9) for _ in range(count)), encoding="utf-8")
    message = _emit(hook, monkeypatch, target)["systemMessage"]
    assert f"{count + 1} finding(s)" in message  # + the fixture's long index row
    assert message.count("block is") == hook.MAX_REPORTED
    assert f"and {count + 1 - hook.MAX_REPORTED} more" in message
