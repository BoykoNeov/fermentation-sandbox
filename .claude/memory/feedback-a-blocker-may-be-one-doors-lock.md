---
name: feedback-a-blocker-may-be-one-doors-lock
description: "A recorded \"we cannot do X because Y\" is usually an argument against ONE implementation of X; before inheriting the decline, ask whether Y is intrinsic to the goal or belongs to the route"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ea3f3770-ac9d-4299-8476-1248aca8b862
  modified: 2026-08-27T08:45:45.032Z
---

When the archive says a repair was declined, read the reason as a claim about a **route**, not
about the goal, until you have checked which it is. A well-argued blocker names a real cost of the
one implementation its author had in mind. A second door may have none of that cost — and because
the decline is written down, articulate, and repeated, it is easy to inherit as a fact about the
problem.

**Why:** D-237 and then D-240 both declined to make a compile-time seed drawable, and both gave
the same reason: the sampler is scoped by `Process.reads`, `reads` also propagates tiers (D-160),
so declaring a seed on a Process asserts something about tiers that nothing had measured. Every
word true. But it is an argument against *declaring the name on a Process* — the sampling job and
the tier job are one field doing two things, and nothing required buying the second to get the
first. D-241 published the names on the compiled scenario instead and unioned them into the
sampled set directly: same repair, no tier claim, and the objection simply does not arise. Two
records had carried the decline forward, the second of them stating it as the reason a priced list
was being *handed over* rather than acted on.

**How to apply:** when you meet a recorded decline, write the blocker as a sentence and ask what
it is quantified over — "we cannot draw this seed" or "we cannot draw this seed *by declaring it on
a Process*". If the second, enumerate the other doors before accepting the first. The tell is a
blocker phrased in terms of a specific mechanism, field or call site rather than in terms of the
physics or the measurement. Note this is the *converse* failure to
[[feedback-check-the-blocker-is-still-blocking]], which is about a blocker that has since been
lifted; this one was never load-bearing in the first place. And when you do take the other door,
say in the record which half of the original objection you avoided and which you still owe —
D-241 took the sampling half and left the tier half explicitly unmeasured, because a repair that
quietly took both would be D-160's own error with the axes swapped.
