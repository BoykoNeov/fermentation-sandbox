---
name: feedback-a-derived-yield-encodes-its-rate-law
description: A constant derived as Δmeasured/Δdivisor is not a measurement of a yield — it silently encodes the rate law you assumed; check that before defending anything it implies
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2d0a7825-2cd3-4bf8-ac09-9bfa58d5add5
  modified: 2026-08-11T18:04:41.144Z
---

A parameter derived as **`Δ(measured quantity) / Δ(divisor)`** is a measured *difference* over
a chosen denominator. It is a **yield** only if the production really tracks that denominator.
That assumption rides inside the number with no provenance field of its own — so the constant
looks measured, its band looks measured (both edges can be real named-strain deltas), and
**every behaviour it implies inherits an authority it never earned.**

**The case (D-183).** `Y_acetic_sugar_beer` = (day-7 − day-0 acetic) / (sugar fermented), both
ends real figure reads on one wort. It implies "finished acetic scales with gravity". When a
new rate law would have changed that scaling, I treated the inversion as a **regression against
a measurement** and spent three design iterations trying to preserve it — a growth-linked pair,
a Luedeking–Piret three-term form, a second state slot — all contorted around protecting it.

The premise was false. Mapping the source's *extract* figure onto its *acid* figure showed
**86 % of the rise happens inside the first 15 % of the sugar flux**, so production never
tracked that divisor. What was measured at day 7 is a **level** (105-126 ppm), not a yield. The
gravity scaling was a consequence of the refuted rate law, not an observation, and there was
nothing to protect.

**Why:** the division is done once, in a beat that cares about the endpoint, and the divisor is
picked because it is the flux the code already has. Afterwards nothing distinguishes
"measured yield" from "measured difference + assumed rate law" — least of all a provenance
`source:` field, which correctly cites the figures both ends came from. So the assumption is
laundered into a citation, and later beats defend it as data. This is
[[feedback-name-the-field-your-predicate-read]] in the parameter layer: the thing that is
sourced is not the thing being claimed.

**How to apply:** when a change would alter what a derived constant implies, first ask **what
the denominator asserts** and whether the source can test it. Write the ratio out longhand —
"grams of X per gram of Y" — and check the source for the *time course* of both X and Y, not
just their endpoints; if one moves mostly while the other is still flat, the ratio is not a
yield. Do this **before** designing around the implication, because the cost is not a wrong
number (the endpoint is reproduced by construction either way) but the blind alleys spent
defending it. State the finding in the parameter's own note so the retired constant carries its
own retirement, and keep it read by a test that **recomputes** the counterfactual
[[feedback-grep-finds-claims-not-guards]]. Related to
[[feedback-re-read-the-source-you-already-mined]] — the check usually needs figure interiors
that are already on disk.
