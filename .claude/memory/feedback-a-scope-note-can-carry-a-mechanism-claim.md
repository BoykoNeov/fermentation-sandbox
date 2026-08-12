---
name: feedback-a-scope-note-can-carry-a-mechanism-claim
description: "A deliberate omission justified by a mechanism welds a world-claim onto a scope decision, and the claim inherits the decision's immunity from review"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 42d75f5f-6caf-4f0a-904b-2781dd61f15c
  modified: 2026-08-12T07:18:19.609Z
---

When a comment explains why something was **left unbuilt**, separate the two sentences it is
really making: **what we chose not to build** (revisit on cost/benefit) and **what we believe
about the world** (revisit on evidence). Welded together, the second inherits the first's
immunity — scope notes get re-read asking *"is this still worth deferring?"*, never *"is this
still true?"*

**Why:** D-44 chose not to track residual copper — a fair author's call — but justified it with
*"the CuS drops out with the lees"*. That is not a scope decision; it is a claim about wine, and
the literature had **already retracted it** ("incorrectly assumed in older textbooks", UWC 2nd
ed. Ch. 24: the products stay dispersed). It survived ~147 records purely because of *where* it
was written. Worse, `add_copper`'s docstring had promoted it to a blocker — *"closing it means
choosing a residual-copper fraction, and **nothing sources one**"* — so a false mechanism was
doing active work keeping the fix out, while the source refuting it sat in the project's own
library, already mined for other beats. The fix took one afternoon once the sentence was read as
a claim rather than as bookkeeping.

**How to apply:** Treat "we didn't build X because <mechanism>" as two reviewable items. If the
mechanism is load-bearing enough to justify the omission, it is load-bearing enough to need a
citation — and a citation is checkable, which a bare mechanism sentence in a scope note is not.
When you meet one, **open the source before accepting the deferral**, especially when the claim
is phrased as unsourceable: "nothing sources one" is a claim about the world's literature, and
it is the single easiest kind to be wrong about. Related but distinct from
[[feedback-check-the-blocker-is-still-blocking]] (there the blocker was silently demolished by
unrelated work and left no ⚠; here the blocker was **never true**, and nobody looked because it
read as a scope note). See also [[feedback-re-read-the-source-you-already-mined]] and
[[feedback-paywalled-is-one-host]]. Landed at D-191.
