---
name: feedback-a-generic-partition-is-not-a-defined-medium
description: Substituting a typical-composition figure for a defined synthetic medium is a category error — and check what a quoted share INCLUDES before deriving a range from it
metadata: 
  node_type: memory
  type: feedback
  originSessionId: b530ceda-f935-4502-a0e8-cac9593dd384
  modified: 2026-08-27T12:46:24.463Z
---

A published *typical* composition (a handbook's "ammonium is 3-10 % of must nitrogen") describes a
class of natural media. A study's **defined synthetic medium** has a composition the paper states.
Filling the second from the first is a category error: it produces a plausible-looking must that is
not the one the measurement was made on, and any comparison built on it scores the model against a
medium nobody used.

**Why:** D-244 needed a composition for two fixtures anchored on published trials, and briefly
composed one from Ribéreau-Gayon's grape-must partition. Two independent faults: the substitution
itself, and a transcription error inside it — the Handbook's "amino acids 25-30 %" **includes
proline**, which is ~48 % of must amino acids and is excluded from assimilable nitrogen by
definition. The derived amino share was therefore wrong in the direction that decided the
conclusion, and the parameter this repo already derives from the same passage carries the same
un-subtracted tension in its uncertainty note.

**How to apply:** when a fixture needs a source's medium, the composition comes from the source or
the fixture keeps the state it was validated on and the gap is recorded. Before deriving a range
from any quoted share, ask what the share's denominator and numerator **include** — a fraction of
"total nitrogen" and a fraction of "assimilable nitrogen" are different numbers, and the
non-assimilable part is often the larger half. Related: [[feedback-a-notes-field-is-unchecked-storage]],
[[feedback-migrate-a-fixture-by-its-state-not-its-intent]].
