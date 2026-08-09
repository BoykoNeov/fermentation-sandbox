---
name: feedback-a-screen-is-not-idempotent-under-its-own-repair
description: Re-running a text screen after fixing what it found can score the FIXED entries worse — the count is not a progress metric
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 9db1bf02-daf7-493f-be8a-c4e84aa60976
  modified: 2026-08-09T21:15:25.286Z
---

A text screen that flags "this note disagrees with its own data" will often flag
the **repaired** note too — because the honest repair is to write *both* accounts
down, and a worst-gap classifier then sees the wider disagreement. At D-170 the
sweep returned **15 of 21 on both trees, the same 15 names**, but four of them had
been *resolved* at D-168: repairing them made their measured gap **grow** (thermal
1.309 → 1.801; `E_a_decarb` 1.609 → 2.809). A naive re-run reads "no progress".

**Why:** the screen measures text against data; a repair that adds context adds
text. Nothing in the string distinguishes a live external claim from a derived
restatement from a quoted, explicitly-retired one — that is the missing
band-provenance schema field (D-164 §6 / D-167 §10 / D-168 §7 / D-170 §2), and its
absence makes the *measurement* read backwards, not just the note read ambiguous.

**How to apply:** never quote a screen's count as a progress metric across a repair
boundary. Run it against the **pre-repair tree** for the denominator (that also
validates the harness — D-170 reproduced D-168's 21/15 exactly, which is what
licensed everything after it), and diff the **membership**, not the number. When
the count is unchanged, check whether the *reasons* are, before reporting either
way. Same family as [[feedback-a-text-screen-has-units-and-self-reference]] and
[[feedback-name-the-field-your-predicate-read]]; pairs with
[[feedback-pre-register-the-cheap-prediction]] — D-170's pre-registered prediction
(that the field-complete re-run would move the count) **failed**, and saying so was
the result.
