---
name: feedback-a-note-can-state-its-span-twice
description: One note stated its band twice — kJ/mol and a Q10 — and the two disagreed; cross-check a note against ITSELF before sourcing it against the world
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6fb0d3b7-0f29-4c9d-9dab-fca8a841b25f
  modified: 2026-08-09T18:46:12.821Z
---

**A provenance note can state the same span in two units or two derived forms. Convert them to
a common basis and check them against each other BEFORE going to the literature** — the internal
check is free, needs no source, and can invalidate the whole external question.

D-168: four of D-167 §5's five "contradicted" bands stated their span twice — once in kJ/mol,
once as a Q10 at 20 °C. Q10 is not independent: the shipped `arrhenius_factor` determines it from
the band. Converted:

| | account A | account B (Q10) → kJ/mol | overlap | shipped |
|---|---|---|---|---|
| thermal ×3 | ~100–150 | ~3–4.5 → [81.2, 111.1] | [100, 111] | [80, 130] — **neither** |
| oak | ~10–25 | ~1.3–1.6 → [19.4, 34.7] | [19.4, 25] | [10, 35] — **neither** |

**Why:** D-167 offered "move the band to its cited range **or** correct the note". The first option
is **ill-posed** when there are two cited ranges and they are near-disjoint — no band satisfies
both. That is not discoverable by sourcing harder; only the internal cross-check finds it. The
thermal three also carried a *third* claim ("banded around a Q10 ~3.5" = 92.6 kJ/mol) against a
nominal of Q10 3.87.

**How to apply:** when a note gives a number in any derived form (Q10, half-life, a ratio, a
percentage of something else), recompute it **through the shipped law** and compare. Disagreement
means the span is AUTHOR-CONSTRUCTED, whatever the `source:` field says.

**Report the denominator and a positive control, or it is not a result.** Archive-wide the check
ran **29** notes mentioning Q10 → **21 CHECKABLE** (`unit: J/mol` + varying band + a Q10 *range*;
8 were point claims constraining only the nominal) → **6 agree, 15 disagree**, and stratified by
the worse endpoint gap that is **6 GROSS / 2 moderate / 7 rounding** — a flat "15 disagree" would
have sold a 0.14 notational miss as a 1.6 defect [[feedback-count-and-print-your-skips]]. The
control matters as much: six notes DO reproduce their band exactly — but they are **one band
`[30k,70k]` reused six times, not six independent successes**, so cite it at that weight
[[feedback-a-null-result-needs-a-positive-control]]. Three GROSS cases were invisible to the
citation-shaped screen that started the beat, because they fail on the **derived axis** while
their citation is fine.

**Run the sweep on the PRE-edit tree.** Rewriting the notes corrupted the harness's own input two
ways: reformatted claims fell out of the regex, and a rewritten note that *quotes its retired
claim* scored as carrying it live — 21 checkable collapsed to 18. Extract with
`git show <pre-commit>:<path>`.

**Do NOT then assign edges to accounts.** It is tempting to read "low edge from account B, high
from A" off proximity — thermal's 80,000 beside B's 81,180, oak's 35,000 beside B's 34,730. That is
inference from a round number adjacent to a bound, and writing it into a shipped note *invents*
provenance: D-167 §4's SELF-RESTATEMENT anti-pattern committed from the authoring side. See
[[feedback-nominal-on-a-band-edge-is-not-inertness]]. Only claim it when the note **names** the
account for that edge, as beer's did. Sibling to [[feedback-name-the-field-your-predicate-read]]
(that one: the predicate read the wrong field; this one: the right field held two answers) and to
[[feedback-a-text-screen-has-units-and-self-reference]] (units again, opposite direction).
