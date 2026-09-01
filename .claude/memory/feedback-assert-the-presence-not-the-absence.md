---
name: feedback-assert-the-presence-not-the-absence
description: A guard that asserts a wrong shape is ABSENT passes on the fault that removes the whole shape; assert the right value is present
metadata:
  type: feedback
---

A guard written as "the axis type is not `log`" was **true while the bug was live**: the
refusal path returned an empty spec, so the axis had no type at all *and* no range — Plotly's
own −1…1 came back, which is the exact fault the feature existed to fix. Asserting the
absence of the wrong shape cannot distinguish "the right thing happened" from "nothing
happened".

**Why:** the two failure modes a guard must separate are usually *wrong value* and *no value*,
and a negative assertion collapses them. It also passes vacuously whenever the code path is
skipped, which is the commonest way a feature silently stops running.

**How to apply:** assert the positive fact — `range[0] == 0.0`, not `type != "log"`. And when
one setting can turn another feature off, loop the guard over both settings rather than
testing the default: this one only reddened once it ran `for log_y in (False, True)`. See
[[feedback-a-guard-can-be-blind-to-the-mutation-it-names]].
