---
name: feedback-a-null-needs-the-predicate-the-claim-makes
description: "A search null is only as broad as the pattern SHAPE it ran; D-217's 75-hit census never queried the second half of its own claim"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 44f9c73d-c9c3-43fa-889b-feafa52223be
  modified: 2026-08-18T07:40:54.703Z
---

**A null result is a claim about a predicate, and the predicate must match every clause of the
sentence you write.** D-217 searched 26 files with 10 patterns — `activation energ`, `arrhenius`,
`Q10`, `kJ/mol`, `van't Hoff`, … — got **75 hits, read all 75**, and wrote: *"zero give an
activation energy, a Q10, **or a convertible two-temperature rate pair**."* The first two clauses
were earned. The third was not: a convertible pair reads *"primary fermentation, 18-22 °C, 3-5
days"* and contains **none of those ten patterns**. The census could not have found one, so the
claim was asserted over a category the predicate never touched. A second census built for that
shape — temperature-near-duration in both orders, attenuation schedules, "ferments in N days" —
returned **28 more distinct hits** the first had never seen.

The null survived (every candidate paired a lager below 12 °C with an ale above 20 °C — two
species, duration attached to neither), so nothing shipped was wrong. **That is the trap.** A big
hit count and a careful read make a census *feel* exhaustive, and the missing clause is invisible
precisely because it produced no hits to notice.

**Why:** the "all 75 read" figure measures diligence on the hits you *found*, never coverage of
the claim you *make*. And the claim outlives the probe — this one went into a shipped
`provenance.notes` field, where [[feedback-a-notes-field-is-unchecked-storage]] applies: unchecked
storage, quoted by later beats as though measured.

**How to apply:** before writing a null, take the sentence apart clause by clause and ask, for each
one, *which pattern would have hit if this clause were false?* If no pattern answers, either add
one and re-run, or narrow the sentence to what was searched. Prefer patterns built from the SHAPE
the evidence would take (a number near a unit near a duration) over the vocabulary you expect an
author to use — an author who has the fact rarely uses your word for it. Then report both
denominators separately, so the next reader can see which half of the claim each one paid for.
Sibling of [[feedback-name-the-field-your-predicate-read]] (that one is the wrong FIELD, this one
is the missing CLAUSE) and [[feedback-count-and-print-your-skips]].
