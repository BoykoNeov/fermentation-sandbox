---
name: feedback-name-the-field-your-predicate-read
description: A census counts whatever FIELD its predicate read — say which one in the headline, or the count silently becomes a claim about a neighbouring property
metadata:
  type: feedback
---

D-163 selected its "buildable" set with one line —
`if "author estimate" in p.provenance.source.lower(): continue` — and reported the
result as **"55 bands / 110 EDGES have a real citation behind them."**

`provenance.source` is the **value's** citation field. The predicate never touched
`uncertainty.low`, `uncertainty.high` or `uncertainty.note`. So the property measured
was *"the value is cited"*; the property reported was *"the edges are cited"*. D-167
reclassified all 110 edges and found **4** with a parameter-specific published span,
and **73** with no external account of the span at all — including all 8 remaining
`acidbase` pKa bands, the slice D-163's own "Next" named as the obvious first target.

The clinching detail: **4 bands inside the "externally sourced" set carry the note
`AUTHOR-ESTIMATED band (x0.3/x3) around a MEASURED MEAN`.** The predicate was
searching for exactly that string, in the wrong field. Archive-wide, **44 of 339**
live bands cite a source while their note declares the span is the author's.

**Why:** a one-line predicate is easy to write and easy to *narrate past*. The gap
between "what the field holds" and "what I want to claim" closes silently in prose,
and the number then travels — into a Next list, into a memory bullet, into the next
beat's plan — carrying a scope it never had. Nobody re-reads the predicate; they read
the sentence. This is [[feedback-check-the-schema-not-the-caller]] from the other end:
there the claim outran the schema, here it outran the field.

**How to apply:** when a census produces a headline number, write the **field name**
into the sentence — "55 bands whose *`source` field* is not `author estimate`", not
"55 bands with a citation behind them". If the property you actually want lives in a
different field, say whether that field exists at all: D-167's real finding was that
**band-provenance has nowhere to live but free text**, which is why the classifier
*could not* have been written correctly, and which D-164 §6 had already flagged from
the override side. Two records reaching the same missing field from opposite
directions is the signal to stop patching predicates. Related:
[[feedback-count-and-print-your-skips]], [[feedback-pin-the-band-not-the-nominal]].
