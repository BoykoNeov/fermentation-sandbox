---
name: feedback-a-shared-fixture-has-two-consumers
description: "Price a fixture change on every study that fixture feeds, not just the one the record is about — the unpriced consumer is where the regression hides"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ab06c4ba-e5e4-4ff7-a3de-6478c4835471
  modified: 2026-08-31T13:02:41.996Z
---

When a record prices a change to a shared test fixture, it prices it on the paper it happens to be
arguing about. Enumerate what else that fixture builds before calling the price complete.

**Why:** D-249 and D-252 both priced moving `commensurate_scenario`'s pitch entirely on Crépin —
her exhaustion timing, her peak biomass, the Coleman anchor. That function builds **both** the
Crépin and the Minebois musts, and D-248's strongest external result (Minebois's two in-study fusel
shares going 1.73×/1.68× to ~1.0×) rides on the same pitch. Two records' worth of pricing said
nothing about it. Measured at D-253 it was fine — 1.009×/0.977× → 1.014×/0.982× — but that was
luck, not diligence: had it regressed, the trade the owner approved would not have been the trade
they got.

The pleasing half of the same fact: the value being moved *to* was sourced from **Minebois**
(1×10⁶ cells/mL), so the shared fixture made the second consumer the one with the citation. A
shared fixture cuts both ways and neither direction is visible from one caller.

**How to apply:** before pricing a fixture edit, grep for the fixture's callers and list the
external comparisons each one feeds. Put the unpriced consumers in the record explicitly, with
their measured before/after, even when the answer is "no change" — a stated null is what stops the
next beat re-opening the question. If a consumer cannot be measured, say the price is partial
rather than letting the record read as complete.

Related: [[feedback-compare-at-the-moment-the-source-sampled]], [[feedback-a-guard-that-hardcodes-an-input-cannot-price-it]].
