---
name: a-guard-that-hardcodes-an-input-cannot-price-it
description: "A guard that writes a fixture's own input in as a literal reads a false cost the moment that input moves - it is not measuring the model, it is measuring the literal"
metadata:
  node_type: memory
  type: feedback
---

The Coleman-yield guard formed its prediction as `0.25 + initial/f_N`, with the fixture's
inoculum written in. At the only other candidate inoculum it subtracted a starting biomass the run
did not have and read 0.925 against its own band of 0.984 +/- 0.02. A whole record (D-249 section 3)
then published "moving the pitch sells the anchor" as a priced cost of a trade offered to the owner.
Against the run's own pitch the same guard reads 0.9848 - inside the band. The cost was the literal.

**Why:** a guard is correct at the configuration it runs at, so a hardcoded input is invisible in a
green suite. It only fires when someone changes that input - which is exactly when someone is
pricing a trade, and exactly when a wrong number does the most damage.

**How to apply:** any quantity a guard's *prediction* needs must be read off the run, never written
in beside it. Before quoting a guard's verdict as the price of moving X, grep the guard for a
literal value of X. And when the repo has an idiom for this (here: `x0 + ...`, twice, in a sibling
file), a diverging line is the bug, not the idiom.
Related: [[feedback-verify-latest-state-not-breadcrumbs]], [[feedback-a-flat-sweep-can-be-a-blind-observable]].
