"""The generated ``docs/DECISIONS.md`` index, and the invariants it depends on.

The archive is ~14.5k lines and gains a record most sessions, so the index is
generated (``tools/gen_decisions_toc.py``). Three failure modes have actually
happened and are pinned here rather than trusted:

* **The index went stale.** ``--check`` existed and was wired to nothing, so the
  index silently drifted behind the records it indexes.
* **Four headings (D-133..D-136) were hard-wrapped across physical lines.** A
  Markdown ATX heading is single-line: GitHub renders only line 1 as the heading
  and the rest as a stray paragraph, so every computed anchor for those four was a
  dead link.
* **A correction pointer could dangle.** The ``**Corrects:**`` / ``**Flags:**``
  markers are hand-written, so a typo'd D-number would put a confident
  "corrected by D-9999" warning on a record and mislead every later reader.

The generator lives in ``tools/`` (not an installed package), so it is loaded by
path.
"""

from __future__ import annotations

import importlib.util
import pathlib
import re
import sys
from types import ModuleType

import pytest

TOOL = pathlib.Path(__file__).resolve().parent.parent / "tools" / "gen_decisions_toc.py"


def _load_generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gen_decisions_toc", TOOL)
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
    """The archive with the generated block removed -- what the generator parses."""
    return str(gen.body_without_index(text))


def test_index_is_current(text: str) -> None:
    """The committed index matches what the generator would write.

    Same assertion CI makes via ``--check``; kept here so a local ``pytest`` run
    catches an appended record whose index was not regenerated.
    """
    rebuilt, _ = gen.build_index(text)
    expected = re.sub(
        re.escape(gen.BEGIN) + r".*?" + re.escape(gen.END),
        lambda _: rebuilt,
        text,
        flags=re.DOTALL,
    )
    assert expected == text, (
        "DECISIONS.md index is stale -- run: uv run python tools/gen_decisions_toc.py"
    )


def test_no_heading_is_hard_wrapped(body: str) -> None:
    assert gen.find_wrapped_headings(body) == []


def test_wrapped_heading_guard_actually_fires() -> None:
    """The guard is only worth having if it catches the thing it was written for."""
    wrapped = "## D-9 — a title that keeps\ngoing onto the next line\n\nbody\n"
    assert gen.find_wrapped_headings(wrapped) == ["D-9 (line 1)"]

    single = "## D-9 — a title on one line\n\nbody\n"
    assert gen.find_wrapped_headings(single) == []


def test_every_record_number_is_unique(body: str) -> None:
    numbers = [int(m.group(2)) for m in gen.HEADING.finditer(body)]
    duplicates = {n for n in numbers if numbers.count(n) > 1}
    assert not duplicates, f"duplicate D-numbers: {sorted(duplicates)}"


def test_correction_markers_point_at_records_that_exist(body: str) -> None:
    """A dangling pointer would render a confident warning about nothing."""
    existing = {m.group(1) for m in gen.HEADING.finditer(body)}
    dangling = [
        target
        for match in gen.MARKER.finditer(body)
        for target in re.findall(r"D-\d+", match.group(2))
        if target not in existing
    ]
    assert not dangling, f"markers point at nonexistent records: {sorted(set(dangling))}"


def test_a_correction_always_postdates_what_it_corrects(body: str) -> None:
    """D-53 may correct D-52; nothing may correct a record written after it.

    Catches the likely typo class (transposed digits) that
    :func:`test_correction_markers_point_at_records_that_exist` cannot, because a
    transposition usually lands on a record that does exist.
    """
    backwards: list[str] = []
    for match in gen.MARKER.finditer(body):
        preceding = list(gen.HEADING.finditer(body, 0, match.start()))
        assert preceding, "a marker appeared before any D-heading"
        source = int(preceding[-1].group(2))
        for target in re.findall(r"D-(\d+)", match.group(2)):
            if int(target) >= source:
                backwards.append(f"D-{source} claims to correct D-{target}")
    assert not backwards, "; ".join(backwards)


def test_marker_regex_parses_both_single_and_comma_list_targets() -> None:
    single = "**Corrects:** D-52 — the clause."
    listed = "**Flags:** D-71, D-74, D-132 — the shared clause."

    got_single = gen.MARKER.match(single)
    assert got_single is not None
    assert got_single.group(1) == "Corrects"
    assert re.findall(r"D-\d+", got_single.group(2)) == ["D-52"]

    got_listed = gen.MARKER.match(listed)
    assert got_listed is not None
    assert got_listed.group(1) == "Flags"
    assert re.findall(r"D-\d+", got_listed.group(2)) == ["D-71", "D-74", "D-132"]


def test_so2_titles_do_not_leak_into_the_oxygen_bucket() -> None:
    """`o[2₂]` must not match inside "so₂"/"co₂".

    Without the lookbehind this swept every SO2 record into "Oxidation, O₂ &
    aging" -- 43 of 137 records in one bucket, including MLF and ester ones.
    """
    sulfur_only = "SO₂ free/bound split: total conserved, free/bound/molecular derived"
    assert "Oxidation, O₂ & aging" not in gen.topics_of(sulfur_only)
    assert "Sulfur, SO₂ & sulfides" in gen.topics_of(sulfur_only)

    genuinely_oxygen = "closure oxygen ingress: a zero-order per-closure OTR"
    assert "Oxidation, O₂ & aging" in gen.topics_of(genuinely_oxygen)


def test_buckets_are_multi_membership_not_a_partition() -> None:
    """D-137 is an O₂ record *and* an SO₂ one; both buckets must claim it."""
    d137 = (
        "the O2 sink partition audit: SO2 takes ~90% of the whole 5-year O2 budget "
        "and never exhausts"
    )
    topics = gen.topics_of(d137)
    assert "Oxidation, O₂ & aging" in topics
    assert "Sulfur, SO₂ & sulfides" in topics


def test_every_record_lands_in_at_least_one_bucket(body: str) -> None:
    """An unbucketed record is a TOPIC_RULES gap, not a fact about the record."""
    unbucketed = [m.group(1) for m in gen.HEADING.finditer(body) if not gen.topics_of(m.group(3))]
    assert not unbucketed, f"no TOPIC_RULES match: {unbucketed}"


def test_generating_twice_changes_nothing(text: str) -> None:
    once, count = gen.build_index(text)
    applied = re.sub(
        re.escape(gen.BEGIN) + r".*?" + re.escape(gen.END),
        lambda _: once,
        text,
        flags=re.DOTALL,
    )
    twice, count_again = gen.build_index(applied)
    assert once == twice
    assert count == count_again
