---
name: feedback-conceded-caveats-are-not-coverage
description: "A record that names its own gap and ships the conclusion anyway: the concession reads as rigour and functions as cover — run the test on every branch the criterion claims to cover"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6c4983d6-15b1-49d6-af70-a34fd86dc75e
  modified: 2026-07-27T10:26:58.680Z
---

When a record states a criterion and then applies it to **some** of the things the
criterion covers, the unmeasured branches must be labelled **asserted**, not written
in the same voice as the measured ones. Worse: naming the weakness in prose does
**not** discharge it.

**Why:** D-138 (Gate 2) built a representability criterion — `[X]_ss = flux / k` vs
`solve_ivp`'s `atol` — and ran it on **one of three** intermediates. `quinone` was
asserted ("lifetimes hours-to-days"), and `Fe(III)` was excluded on a *sourcing*
argument while the record's own sentence conceded: "total-iron surplus and Fe(III)
quasi-steady-state are **different claims**." That concession looked like rigour and
worked as cover — it let an untested branch ship inside a record whose headline was a
measurement. The conclusion happened to survive, but for a **different reason** than
the one recorded: Fe(III) *passes* representability by 400–1500x and is QSS on a 63x
timescale separation, not on iron surplus. A conclusion that is right by luck is not
sourced.

**The selection variant, caught one pass later in the same record.** The fix for the
above quoted Nguyen Table 3.1 as a "63x" timescale separation, taking the *slow* limb
explicitly for conservatism — while sitting on the **no-copper row of a six-cell grid**.
The real range was **3.3x–63x**. A conservative-sounding choice on one axis (slow limb)
concealed a non-conservative one on another (favourable row), and the conservatism
narrated in the prose is exactly what stopped it being checked. **When a number comes
from a table, transcribe the whole table** — the defect and its fix are the same shape,
so "I already corrected this record" is not evidence the next number in it is sound.

Two sub-lessons that generalise past this record:

- **Necessary is not sufficient.** Failing `atol` *forces* QSS (no choice exists);
  passing it does not *force* a state slot. Two intermediates can both end up QSS for
  structurally different reasons, and collapsing them under one label lends the weaker
  verdict the stronger one's certainty. Say which argument carries which.
- **A derived number sitting next to citations that measured it reads as corroborated.**
  D-138's Fe(II):Fe(III) ~2000:1 was derived, printed beside Danilewicz 2016b/2018 which
  *measured* that ratio, neither pulled. Label it unchecked or pull the paper — same
  failure shape as [[feedback-paywalled-is-one-host]] and
  [[feedback-transcribe-tables-not-prose]].

**The detection signal, which is narrower than the fix.** All three defects were caught
by review, none by a check that could have been run. What they shared: the prose
**narrated its own rigour** at exactly the point the evidence was thinnest — "taking the
*slow* limb, the one that governs", "that is a different claim", "**both** fail". So the
cheap trigger is not re-checking everything; it is: **wherever a draft argues for its own
carefulness, check the number under it.** That fired three times in one record.

**How to apply:** before closing a record built on a criterion, list every branch the
criterion claims to cover and check each has a run or an explicit *asserted* label with
a named pull. If a sentence in the draft concedes a gap ("that is a different claim",
"strictly this doesn't follow"), treat it as a **blocking to-do**, not as having handled
it. Relatedly: the fix for an unmeasured branch is often one table away —
[[feedback-transcribe-tables-not-prose]] supplied the Fe(III) constant here, third time
that has been the answer. See also [[feedback-check-the-schema-not-the-caller]] for the
other self-sealing shape.
