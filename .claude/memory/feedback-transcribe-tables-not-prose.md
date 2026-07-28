---
name: feedback-transcribe-tables-not-prose
description: "Source numbers from a paper's TABLES, transcribed — the prose disagreed with its own table twice, once on a value and once on the direction"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d9476cef-cba5-404b-a92e-ba999b450f5b
  modified: 2026-07-27T18:27:55.788Z
---

When a paper is going to back a `source:` field, **transcribe its tables** rather
than sourcing from its prose, abstract, or a search summary. Two measured
failures, one beat apart:

- **D-135** — transcribing the tables **changed 4 numbers** that had been taken
  from the same paper's prose.
- **D-136** — a paper's prose contradicted **its own table's direction**. The
  table wins. (Consequence that survives: technical cork ships *below* screwcap.)
- **D-102** — the extreme case: `E_a = 128 ± 37 kJ/mol` was banked from a **search
  summary** and shipped, flagged in its own entry as "must be read before it
  ships". The paper's equation (6) says **186 ± 12**, and forbids the transfer
  anyway.
- **D-142 → D-143** — the worst kind, because it looked cited: a band was pinned
  in a test cited to "Miao Table 3" while Tables 1/2 were transcribed and **3 was
  not**. Transcribing it showed the caption named the METHOD ("by the SO₂ addition
  method"), i.e. a secant of the equilibrium locus — and the sim-side number it was
  compared against was an oxidation-path slope at 2–6× the reference's free SO₂.
  Correcting the statistic moved the sim from "low against every one of eight
  wines" to **inside the band**, and withdrew a blocking prerequisite.
  **A number cited to an untranscribed table is not sourced.** Transcribe captions
  and axis labels too — the method lives there, and comparing two different
  statistics is a bigger error than getting a digit wrong.

**Why:** prose rounds, summarizes, and restates from memory; a table is the
measurement. A number sourced from prose carries a provenance claim the paper does
not actually support, which is a Prime-Directive-2 violation wearing a citation.
And "the paper says X" is unfalsifiable in review unless the table is written down
somewhere a reader can check.

**How to apply:** before a number enters a YAML `value:`/`source:`, extract the
table it came from into the working notes (`M:\claud_projects\temp\...`) and cite
table + row. When prose and table disagree, **the table wins** — and say in the
D-record that they disagreed, because it usually means the prose is a
generalization that does not transfer. Watch for mixed semantics *within* one
figure: D-119's Fig. 6A had 2-PE as a TOTAL while the other bars were UNLABELLED
segments; asserting a uniform rule and then sanding a "2.57% ≈ 2.5%" mismatch to
protect it is the tell — derive each bar from its own printed constraints.
Getting the paper at all is [[feedback-paywalled-is-one-host]].
