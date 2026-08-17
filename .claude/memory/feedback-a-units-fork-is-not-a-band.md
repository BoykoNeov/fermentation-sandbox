---
name: feedback-a-units-fork-is-not-a-band
description: A parameter band built from two readings of one printed number is an unresolved unit disguised as measurement error - resolve it and name the loser
metadata: 
  node_type: memory
  type: feedback
  originSessionId: fb2737c8-b92e-460e-9204-f352f202c353
  modified: 2026-08-17T10:49:54.472Z
---

At D-209 the load-bearing input was one printed table row — *"Ammonia 25-30"* under a column
headed only `mg/L`. I could not tell whether it meant mg **N**/L or mg **NH₃**/L, so I carried
**both** readings forward and crossed them with the printed 25-30 and with three pH values: a
tidy 12-cell grid that I reported as the parameter's band, **0.1465-0.1988**. It was not a band.
Nearly the whole spread *was* the units fork, and the derived-value nominal would have been the
midpoint of two mutually exclusive readings of the same number.

**Why:** an uncertainty band is a claim that the true value is somewhere inside because the
measurement is imprecise. A units fork is a claim that I don't know what the source said. Folding
the second into the first makes ignorance look like precision, hides the one thing a later reader
would need to check, and — worst — lets the nominal be a value **no reading of the source
supports**. It also inflates the band, which makes the parameter look conservatively handled when
it is the opposite. D-209's fork was only ~13 % wide; the shape is what matters, not the size.

**How to apply:** when a band's width traces to *which reading you took* rather than to
*measurement scatter*, stop and resolve the reading. Usually a neighbouring row in the same table
settles it: at D-209, summing a per-species table as elemental nitrogen landed **inside** the same
column's printed free-amino-acid figure while the compound-mass reading fell far outside, so the
column was mg N/L and the ammonia row two lines down read the same way. Then ship **one** reading
as nominal with its own real band, and record the loser as a **named scope alternative** in the
`uncertainty.note` — not as extra width. Related: [[feedback-a-text-screen-has-units-and-self-reference]]
(units are what a text screen misses), [[feedback-a-note-can-state-its-span-twice]] (a note whose
two statements of its own span are near-disjoint), [[feedback-transcribe-tables-not-prose]].

A band and a fork also fail differently, which is the tell: widen a band and every member stays
physically possible; widen a fork and half the members describe a source that does not exist.
