---
name: feedback-a-scope-note-can-size-the-wrong-half
description: "A prior record's scope note named two omissions and sized only one; the sized half moved the endpoint by exactly 0.0 and the unsized one owned all of it"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 41ccdfd6-1486-4600-8cb4-43288814b4cf
  modified: 2026-08-17T12:43:33.806Z
---

When an earlier record's scope note names **more than one** omission and attaches a magnitude to
one of them, do not inherit that magnitude as the beat's headline. Ask, before building, **which
omission the observable actually depends on** — and check the ones it did *not* size.

**Why:** D-209 §8c named both of `add_dap`'s gaps and sized the nitrogen one: *"a DAP-dosed wine's
acidification is understated by roughly 3× on the dosed fraction"*. Measured at D-210, that half
moves a dry wine's **final pH by exactly 0.0** — the dosed nitrogen is fully consumed either way,
so the pool's charge returns to ~0 whatever its charge-per-mole was. What it moves is the
**excursion** (+0.235 pH at the dose instant). The half the note filed second and did not size,
the dropped phosphate counter-anion, owns **100 %** of the permanent change (−0.162 pH). The two
also oppose at the dose instant, so shipping the sized half alone would have reported +0.235 pH as
the fix when the other cancels 71 % of it. The note was not wrong about anything it asserted; it
attached its number to the half that cannot reach the endpoint.

**How to apply:** for each named omission, ask what it does to the quantity you will report at the
*horizon* you will report it at, not at the instant it acts. A term carried by a pool that empties
is an excursion term at a dry endpoint and an endpoint term only where the pool is left standing
(stuck ferment, sweet wine, short horizon — check that case explicitly, it is a real scenario, not
an edge case). Then isolate each half in **its own** slot/flag rather than measuring
with-and-without the whole beat: a single number credits each half with the other's cancellation
[[feedback-build-the-term-that-makes-agreement-worse-first]]. Related: a scope note's mechanism
claim inherits the scope decision's immunity from review
[[feedback-a-scope-note-can-carry-a-mechanism-claim]].
