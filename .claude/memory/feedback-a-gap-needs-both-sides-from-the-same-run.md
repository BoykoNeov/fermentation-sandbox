---
name: feedback-a-gap-needs-both-sides-from-the-same-run
description: "A published number is a baseline ANCHOR, never the other side of a live difference — pairing a new arm against an archived value measured the change I had just made, not the effect I was attributing"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 08a1d6fe-e428-435f-8525-f1307096e639
  modified: 2026-09-04T14:48:31.077Z
---

**When you report a gap between two arms, both sides must come out of the same run.** A number the
archive published is the right thing to *reproduce*; it is the wrong thing to *subtract from*.

At D-268 the joint fusel arm's gate-closure fraction had to stay close to the blend arm's — the
claim being that swapping the precursor sink does not move production. I read the joint arm's new
36.59 / 48.50 / 57.03 % and subtracted D-257's **published** blend numbers, 36.3 / 48.1 / 55.9,
got 0.29 / 0.40 / **1.13** points, decided the sink had started moving production at the highest
nitrogen, **widened that guard's bound from 1.0 to 1.5**, and wrote a receipts paragraph
explaining the widening. Every part of that was wrong. The blend arm measured in the same run
reads 36.48 / 48.39 / 56.92, so the real gaps are 0.107 / 0.110 / 0.114 — flat, and unchanged by
a composition change that moved the Crépin splits by 6-11 points. What I had "measured" was
almost entirely my own edit leaking through the published number's older composition.

**Why:** it is the exact inverse of
[[feedback-reproduce-a-published-number-before-trusting-the-new-column]], and the two are easy to
confuse because both involve an archived value. There, the published number is the *anchor* you
check the harness against before trusting anything. Here it was silently promoted into one
*operand* of a live difference — and an archived number carries the whole configuration it was
measured under, so any difference against it contains every change made since. The failure is
invisible: the arithmetic runs, the magnitude is plausible, and a loosened bound is exactly the
kind of edit nobody re-derives.

**How to apply:** before writing `a - b`, ask where `b` came from *in this process*. If it came
from a document, it is not `b` — go and read the corresponding arm from the same fixture. Keep
published values for two jobs only: reproducing them on the *superseded* configuration (D-268
keeps `growth_superseded` / `joint_superseded` arms live for precisely this), and pinning what a
record claims. And treat "I am about to widen a bound to admit a number I just measured" as the
tripwire it is — the widening is the moment to re-measure, not to explain. What caught it here
was the *extra* assert I added to pin the shape of the gap I thought I had found: it failed and
printed 0.11 / 0.11 / 0.11. Related: [[feedback-pair-the-arm-with-its-baseline]],
[[feedback-a-baseline-log-goes-stale-between-edits]].
