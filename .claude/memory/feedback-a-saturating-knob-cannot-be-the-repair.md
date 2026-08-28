---
name: feedback-a-saturating-knob-cannot-be-the-repair
description: A rate proportional to a state variable cannot outrun that variable however large its coefficient — sweep the knob to absurdity and show it saturates before pricing or deferring a fix
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6025169e-c573-4f1a-af77-db80fef3cfd0
  modified: 2026-08-28T07:44:59.938Z
---

A shipped rate law read `rate ∝ r · X` (X = biomass). A record had pinned the *extent* as
insensitive to `r` across a 200× sweep and carefully said the **timing** was a separate matter
that missed by ~1.5×, which read as "the knob is untouched, so the timing repair is still
available". It was not: swept to **1000×** the shipped value, the exhaustion time saturated at
~16.4 h and the gap topped out at 1.71×, never overtaking the enclosing run's 1.94×. The reason
is structural, not numerical — uptake proportional to the biomass performing it can never consume
faster than cells accumulate, so the knob cannot reach the quantity it would be fitted to.

**Why:** "insensitive across a sweep" is normally read as *good* — it is what makes a shipped
value a bound rather than a fit. But the same insensitivity is also a statement that the
parameter cannot deliver a repair, and the two readings live in different sections of the same
record. Absent an explicit sweep on the *missing* observable, a knob stays on the open list as a
deferred repair for as long as nobody proves it inert there too.

**How to apply:** when a record says "extent is insensitive, timing misses", sweep the same knob
against the **timing** and take it to absurdity — three orders of magnitude, not the tidy 200×
that made the extent point. If it saturates, the finding is stronger than a refusal on cost: the
repair is *unreachable*, which is a one-way answer and should be recorded as such rather than
deferred [[feedback-price-the-refusal-not-just-the-error]]. Then read the saturation as
information about the **functional form**: a coefficient that cannot move an outcome means the
shape is wrong, and that is the thing to name for the next beat. Related:
[[feedback-measure-which-side-before-building]],
[[feedback-read-a-channels-timing-against-its-own-run]].
