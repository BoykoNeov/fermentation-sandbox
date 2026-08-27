---
name: feedback-migrate-a-fixture-by-its-state-not-its-intent
description: "When a semantics change forces fixtures to be re-declared, find the re-declaration that preserves their state bit-for-bit; re-authoring what you think they meant invents data"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: b530ceda-f935-4502-a0e8-cac9593dd384
  modified: 2026-08-27T12:46:15.074Z
---

A change to what a scenario field *means* forces every existing fixture to be re-declared. There is
usually one re-declaration that leaves the fixture's actual state **bit-identical** — find it and
apply it everywhere, mechanically. Re-authoring a fixture's intent instead ("what composition did
this study really have?") invents inputs, and every number downstream then rests on them.

**Why:** D-244 made `yan_mgl` the total assimilable nitrogen rather than the ammonium share. Adding
each dose's nitrogen to its declared YAN preserves the pitch state exactly, so the only thing that
moves is the quantity the beat set out to correct — which is what makes the resulting reds
readable. For two source-tied fixtures a different treatment was tried first: hold the declared
total at the paper's own number and compose the must from a generic partition. It produced a
headline ("the model fails its sourced floors on any realistic must") that rested on a dose nobody
published, and the sourced range it leaned on was itself wrong. The mechanical migration reproduced
the same qualitative result at a third of the magnitude.

**How to apply:** ship a public helper that computes the migration (`amino_acid_dose_nitrogen_mgl`
is one) so fixtures call it instead of carrying hand-computed constants, and so real users can
migrate too. Where the mechanical migration is knowingly imperfect — here it preserves a
commensurability violation the fixture's own comment forbids — **record the residue in the fixture
and the decision** rather than repairing it from data you do not have. Related:
[[feedback-a-generic-partition-is-not-a-defined-medium]], [[feedback-closer-to-reality-decides]].
