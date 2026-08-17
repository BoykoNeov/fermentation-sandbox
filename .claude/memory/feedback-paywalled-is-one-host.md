---
name: feedback-paywalled-is-one-host
description: "\"Paywalled\" is a property of ONE host, not of a paper — check author-hosted, institutional, thesis and trade-reprint copies before recording a source as blocked"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d9476cef-cba5-404b-a92e-ba999b450f5b
  modified: 2026-07-28T10:09:13.860Z
---

Never record a source as "unreadable" / "paywalled" / "source-blocked" on the
evidence of one host returning 402/403. In this project that call has been wrong
**six times**, each time on a source that was openly available elsewhere:

- **D-123** — Ramey & Ough 1980: author-hosted scanned PDF at `rameywine.com`.
- **D-135** — Franco-Luesma: institutional repository (`zaguan.unizar.es/record/56225`).
  Reading it also revealed the citation itself was wrong — *Food Chemistry*
  199:42-50, **not** JAFC, as D-101 and `mercaptans.py` both had it.
- **D-136** — Lopes 2007: ACS paywalled, but the authors' own trade reprint prints
  the same table.
- **D-137** — Nguyen & Waterhouse 2021: **the first author's thesis on
  eScholarship IS the paper**, chapter 4 verbatim.
- **D-151** — Carrasco-Quiroz 2022 (*Foods*, doi 10.3390/foods11131961): D-150 recorded
  its Table 2 "not recoverable in this beat's fetches". **MDPI 403s automated fetches**,
  but it is gold OA — the **PMC deposit** (PMC9266014) and the **EuropePMC
  `fullTextXML`** endpoint both serve the full tables, and transcribing from both
  independently is its own cross-check. This one was worse than a false block: the
  record had already *named* the table as the one unlock the beat needed.

- **D-208** — a **STANDARDS COLLECTION**, not a paper: MEBAK Band II 2.14, priced at ~60 %
  unreachable in the pre-registration. The publisher's own site serves the method page, and
  its decarbonation sentence is what closed the beat; the equivalent Analytica-EBC 9.35 scope
  line is public too. **A standard is published by a body with a website** — check it before
  pricing a purchased collection as a paywall.

**Why:** a false block is not a neutral gap — it silently redirects the work into
guessing. D-101 recorded bottle reduction as source-blocked and *guessed* a
thioacetate/disulfide mechanism; the paper that would have corrected it was open
the whole time, and D-135 had to overturn the mechanism. A deferral can be right
while its stated reason is wrong, so the reason has to be checked on its own.

**How to apply:** before writing "blocked" into a D-record, a `source:` field, or
a deferral, try in order: the author's own site, the institutional/university
repository, **the first author's thesis** (chapters are often the paper verbatim),
and trade-journal reprints. `WebFetch` returns binary for several of these — fall
back to `pdftotext -layout`. Only after those fail is "blocked" a finding, and
say which hosts were tried. Pairs with
[[feedback-transcribe-tables-not-prose]]: once you do get the paper, the tables
are what you source from.
