---
name: feedback-check-the-schema-not-the-caller
description: "\"The model can't represent X\" is a claim about the whole schema — never infer it from the one Process that happens not to read X"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2e1a5b62-b816-4f57-a993-5f119bd75a95
  modified: 2026-07-22T14:35:26.985Z
---

When about to record that the model **cannot represent** something — no state
variable, no pool, no carrier — verify it against the **schema** (`StateSchema`
in `src/fermentation/core/state.py`, plus the `chemistry.py` molar-mass table),
not against the one Process I happened to be reading.

**Why:** at D-137 I concluded "the counter-sign term on acetaldehyde has no
state variable to live in," reasoning from a verified fact — the browning
Process reads `tannin + anthocyanin`, and Gislason's radical quench needs a
*cinnamate side-chain double bond*, which neither has. The inference from "this
Process doesn't read it" to "the 93-slot schema can't carry it" was never
checked, and it was wrong: `hydroxycinnamics` (D-40) and `ferulic_acid` (D-55)
are live slots carrying exactly that structure, at ~10–200 mg/L. One grep would
have settled it.

The failure is asymmetric, which is what makes it worth a rule. A false
*positive* ("the model has X") dies at the first test. A false *negative* ("the
model can't do X") is self-sealing — it silently converts a cheap change into a
schema extension, gets written into a decision record as a constraint, and
nobody re-tests it because it reads as settled. It also hides couplings: here
the real finding was not a missing slot but that `brett.py` already consumes
both pools, so an antioxidant quench would **compete with Brett** for the same
substrate. That is a live design constraint I would have missed entirely.

**How to apply:** before writing "cannot represent" / "no carrier" / "would need
a new slot," grep the schema and `chemistry.py` for the species *and its
chemical class*, then check who already reads and writes it — an existing
consumer is a coupling that must be declared, not discovered. Reading the schema
is not building; it does not breach a "no code before the gates" hold, so there
is never a reason to skip it. Related: [[feedback-rejected-values-must-be-unreachable]]
(a green suite won't catch it either) and
[[feedback-verify-latest-state-not-breadcrumbs]] (same root: trusting a stale
secondhand summary over the current artefact).
