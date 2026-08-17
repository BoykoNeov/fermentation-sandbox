---
name: feedback-gate-both-halves-of-a-pair
description: "Two halves of one physical thing rode different gates, so a dose booked the anion alone AND opened the gate that suppressed the cation"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 41ccdfd6-1486-4600-8cb4-43288814b4cf
  modified: 2026-08-17T13:02:17.498Z
---

When you add the second half of a physically inseparable pair, check **what gates each half** —
not just that each half is individually right. If the existing half rides an opt-in condition and
the new one rides nothing, you have built a mechanism that can fire at half strength, and the half
that fires is whichever the guard does not cover.

**Why:** at D-210 a dose's ammonium charge went through a function gated on "did a scenario supply
pH information?", while its phosphate counter-anion went in as an ordinary acid slot that the
charge balance reads unconditionally. Two failures, and the second was not on anyone's list:

1. on a beverage with an empty charge balance the dose booked the **anion alone**, which is more
   acidic than the chemistry allows — the exact mirror of the artefact that gate exists to prevent;
2. writing an acid slot **opened the gate**, because acid slots are what the gate tests. So a
   *nutrient addition* switched on a whole unrelated term for a beverage whose pH no scenario had
   supplied: 1.43 pH at the dose instant, 0.65 pH at the end.

The fix was atomicity — both writes ride the gate together, or neither happens. And the scenario it
broke was a **benchmark's own shape**, whose assertions stayed green because they score a quantity
that barely reads the affected one.

**How to apply:** when a change adds a second write to a coupled system, enumerate the gates on
*every* member and make the group atomic. Prefer "all or nothing" over "each is separately
correct": if a later beat opens the gate, the all-or-nothing version is missing a term (recoverable)
while the split version has a half-booked one (an artefact). And check the *predicate* the gate
tests against the *slots* you are now writing — a guard whose condition is "is any of this
populated" turns any new write into a way to satisfy it. A green benchmark is not evidence here
[[feedback-grep-finds-claims-not-guards]]: check whether its assertions read the quantity you moved
[[feedback-a-summary-statistic-is-not-the-curve]].
