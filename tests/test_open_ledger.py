"""The generated ``docs/OPEN.md`` ledger, and the marker grammar it depends on.

``docs/OPEN.md`` is the derived answer to "what is scientifically open right now",
built by ``tools/gen_open_ledger.py`` from two signals only: xfail/skip markers under
``tests/`` and ``**Flags:**`` / ``**Unflags:**`` markers in ``docs/DECISIONS.md``. The
failure modes pinned here are the index's, transplanted:

* **The ledger goes stale.** ``--check`` is only a guard if something runs it, so the
  same assertion is made here for a local ``pytest`` run.
* **A flag pointer dangles or runs backwards.** A ``Flags:`` on a nonexistent record
  would put a confident open item on nothing; one on a *later* record is the transposed-
  digit typo that an existence check cannot catch.
* **An ``Unflags:`` retires a flag nobody declared.** The marker is new, so its grammar
  (must postdate D-F; D-F must actually flag the named target) is exercised on a
  synthetic archive rather than trusted.
* **A reason goes uncaptured.** The two strict xfails in ``test_organic_acids.py`` and the
  three in the Herzan benchmark are the ledger's reason for existing; if the AST walk
  loses them, the ledger is quietly empty.

The generator lives in ``tools/`` (not an installed package), so it is loaded by path.
This file deliberately imports nothing from ``fermentation``.
"""

from __future__ import annotations

import importlib.util
import pathlib
import re
import sys
from types import ModuleType

import pytest

TOOL = pathlib.Path(__file__).resolve().parent.parent / "tools" / "gen_open_ledger.py"


def _load_generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gen_open_ledger", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# A path-loaded module is untyped as far as mypy is concerned, so every attribute
# read off `gen` is `Any` -- hence the explicit `str(...)` coercions below.
gen = _load_generator()


@pytest.fixture(scope="module")
def text() -> str:
    return str(gen.DECISIONS.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def body(text: str) -> str:
    """The archive with the generated index removed -- what the generator parses."""
    return str(gen.toc.body_without_index(text))


def test_ledger_is_current(text: str) -> None:
    """The committed ledger matches what the generator would write.

    Same assertion CI makes via ``--check``; kept here so a local ``pytest`` run catches
    a new xfail or ``Flags:`` marker whose ledger was not regenerated.
    """
    assert gen.OPEN.exists(), (
        "docs/OPEN.md is missing -- run: uv run python tools/gen_open_ledger.py"
    )
    rebuilt = str(gen.generate(text))
    assert rebuilt == gen.OPEN.read_text(encoding="utf-8"), (
        "docs/OPEN.md is stale -- run: uv run python tools/gen_open_ledger.py"
    )


def test_generating_twice_changes_nothing(text: str) -> None:
    assert gen.generate(text) == gen.generate(text)


def test_every_marker_target_exists_and_is_earlier(body: str) -> None:
    """Checked two ways: the generator must accept the real archive, and the raw
    ``Flags:``/``Unflags:`` pointers must each name an existing, earlier record."""
    gen.parse_archive(body)  # raises LedgerError on any refused marker

    existing = {int(m.group(2)) for m in gen.HEADING.finditer(body)}
    bad: list[str] = []
    # Only the TARGET groups: a clause may legitimately cite its own record (D-159's does),
    # and that is prose, not a pointer.
    pointers = [
        ("flags", m.start(), m.group(2)) for m in gen.MARKER.finditer(body) if m.group(1) == "Flags"
    ] + [
        ("unflags", m.start(), f"{m.group(1)} {m.group(2) or ''}")
        for m in gen.UNFLAGS.finditer(body)
    ]
    for verb, offset, targets in pointers:
        preceding = list(gen.HEADING.finditer(body, 0, offset))
        assert preceding, "a marker appeared before any D-heading"
        source = int(preceding[-1].group(2))
        for target in re.findall(r"D-(\d+)", targets):
            if int(target) not in existing:
                bad.append(f"D-{source} {verb} nonexistent D-{target}")
            elif int(target) > source or (int(target) == source and verb != "flags"):
                bad.append(f"D-{source} {verb} D-{target}, which is not earlier")
    assert not bad, "; ".join(bad)


def test_a_record_may_flag_its_own_section_but_not_unflag_itself() -> None:
    """D-164..D-167 each carry `**Flags:** D-16n §k` -- a reversal identified in that very
    record and deliberately not shipped. That is an open item, so a self-flag is accepted;
    a self-Unflags would retire a flag the same record declares and is refused."""
    archive = gen.parse_archive(_archive("", "**Flags:** D-2 §6 — own section, not shipped.\n"))
    assert [(f.source, f.target) for f in archive.open_flags] == [(2, 2)]

    with pytest.raises(gen.LedgerError, match="not an earlier record"):
        gen.parse_archive(
            _archive("", "**Flags:** D-2 — own.\n**Unflags:** D-2 — and retired by itself.\n")
        )


def _archive(*records: str) -> str:
    """A synthetic body: each argument is one record's text under `## D-n — n`."""
    return "".join(f"## D-{i} — record {i}\n\n{extra}\n" for i, extra in enumerate(records, 1))


def test_unflags_names_a_declared_flag() -> None:
    """The Unflags grammar, on a four-record archive where D-2 flags D-1."""
    flagged = "**Flags:** D-1 — the clause.\n"

    # Valid: a later record retires D-2's flag, and exactly that pair moves.
    archive = gen.parse_archive(_archive("", flagged, "", "**Unflags:** D-2 — done.\n"))
    assert archive.open_flags == []
    assert [(r.flag.source, r.flag.target, r.retired_by) for r in archive.retired] == [(2, 1, 4)]

    # Valid, narrowed: `on D-1` names a target D-2 really flags.
    archive = gen.parse_archive(_archive("", flagged, "", "**Unflags:** D-2 on D-1 — done.\n"))
    assert archive.open_flags == []
    assert len(archive.retired) == 1

    # Nothing retired when no Unflags exists: the pair stays open.
    archive = gen.parse_archive(_archive("", flagged, "", ""))
    assert [(f.source, f.target) for f in archive.open_flags] == [(2, 1)]
    assert archive.retired == []

    # Refused: Unflags on a LATER record (D-2 cannot retire what D-4 will declare).
    with pytest.raises(gen.LedgerError, match="not an earlier record"):
        gen.parse_archive(_archive("", "**Unflags:** D-4 — early.\n", "", flagged))

    # Refused: Unflags on a record that declares no flag at all.
    with pytest.raises(gen.LedgerError, match="declares no"):
        gen.parse_archive(_archive("", flagged, "", "**Unflags:** D-1 — nothing.\n"))

    # Refused: narrowing to a target the flagging record never flagged.
    with pytest.raises(gen.LedgerError, match="never flagged D-3"):
        gen.parse_archive(_archive("", flagged, "", "**Unflags:** D-2 on D-3 — wrong.\n"))

    # Refused: a Flags marker on a later or nonexistent record.
    with pytest.raises(gen.LedgerError, match="not an earlier record"):
        gen.parse_archive(_archive("**Flags:** D-2 — forward.\n", ""))
    with pytest.raises(gen.LedgerError, match="does not exist"):
        gen.parse_archive(_archive("", "**Flags:** D-9 — dangling.\n"))


def test_unflags_regex_parses_both_forms() -> None:
    whole = gen.UNFLAGS.match("**Unflags:** D-137 — the cascade shipped.")
    assert whole is not None
    assert whole.group(1) == "D-137" and whole.group(2) is None

    narrowed = gen.UNFLAGS.match("**Unflags:** D-137 on D-71, D-74 — two of five.")
    assert narrowed is not None
    assert narrowed.group(1) == "D-137"
    assert re.findall(r"D-\d+", narrowed.group(2)) == ["D-71", "D-74"]


def test_xfail_reasons_are_captured() -> None:
    """The five strict xfails the ledger exists for, with the records their reasons name."""
    markers, _ = gen.collect_test_markers()
    by_name = {str(m.name): m for m in markers}

    organic = {
        "test_the_model_ferments_tyrells_wort_on_tyrells_schedule",
        "test_the_three_flux_linked_acid_courses_are_mistimed",
    }
    herzan = {
        "test_maturation_sulfiting_raises_total_acetaldehyde",
        "test_bottling_sulfiting_raises_total_acetaldehyde",
        "test_the_later_stages_separate_the_sulfited_must_variants",
    }
    missing = (organic | herzan) - set(by_name)
    assert not missing, f"xfails the AST walk lost: {sorted(missing)}"

    for name in organic:
        m = by_name[name]
        assert m.file == "tests/test_organic_acids.py"
        assert m.kind_label == "strict xfail" and not m.benchmark
        assert "D-215" in m.records
        assert m.reason.startswith("D-215:")
    for name in herzan:
        m = by_name[name]
        assert m.file == "tests/benchmarks/test_validation_herzan_acetaldehyde.py"
        assert m.kind_label == "strict xfail" and m.benchmark
        assert "D-188" in m.records
        assert m.reason.startswith("D-188:")

    # Benchmarks sort first, and the section is file-then-line within each half.
    kinds = [bool(m.benchmark) for m in markers]
    assert kinds == sorted(kinds, reverse=True)


def test_body_skips_are_listed_not_hidden() -> None:
    """A `pytest.skip()` inside a test body is conditional: counted, never an open item."""
    markers, body_skips = gen.collect_test_markers()
    census = [b for b in body_skips if b.file == "tests/test_banded_undrawn_census.py"]
    assert len(census) == 2
    assert all(b.kind == "skip" and "parameter set" in b.reason for b in census)
    assert not any(m.file == "tests/test_banded_undrawn_census.py" for m in markers)
