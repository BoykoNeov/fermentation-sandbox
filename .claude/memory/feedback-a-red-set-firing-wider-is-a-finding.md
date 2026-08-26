---
name: feedback-a-red-set-firing-wider-is-a-finding
description: "When a mutation arm turns more tests red than pre-registered, the extra REDs are evidence about coupling — record them, never widen the prediction to match"
metadata:
  node_type: memory
  type: feedback
---

A pre-registered RED set that fires **wider** than predicted is a result, not a bookkeeping
error. The extra REDs name a coupling you did not know the knob had. **Record the miss and what
it revealed; never quietly widen the prediction to make the arm "match".**

**Why:** D-230's arm A set beer's assumed wort nitrogen to the value that would land Tyrell's
counts, pre-registering four REDs. It produced **eight**. The four extra were every rate and
timing test in the file — because `mu_max` is fitted on a growth fraction *normalised on the
peak*, and the peak IS the nitrogen-limited ceiling. So wort nitrogen is **not separable from
the rate fit**: moving it re-opens a calibration two records back. That was the strongest reason
not to move the constant, and it existed only in the gap between the prediction and the run.
Widening the prediction to "8 tests" would have deleted it.

**How to apply:** when observed REDs ⊃ predicted, do not edit the prediction. Name each
unpredicted test, say what it scores, and ask what shared quantity links it to the mutation —
that answer is usually a structural claim worth more than the arm. Then put it in the guard's own
docstring, where the next beat reads it. Observed ⊂ predicted is the opposite problem (a guard
that cannot see) and is not this. Related:
[[feedback-pre-register-the-cheap-prediction]], [[feedback-full-suite-before-green]],
[[feedback-verify-the-restore-between-mutation-arms]].
