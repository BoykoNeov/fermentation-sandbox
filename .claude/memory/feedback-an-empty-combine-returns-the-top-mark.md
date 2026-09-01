---
name: feedback-an-empty-combine-returns-the-top-mark
description: "A min-combine over an empty list returns the BEST tier, so any quantity nothing touches reads `validated` — the one mark nothing has earned. Every output surface must ask what is DRIVEN before it reads a tier"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 86025c03-b1ae-45b8-b84a-18332d969846
  modified: 2026-09-01T12:15:12.464Z
---

**`min` over an empty list is the top of the scale, so "nothing made a claim about this" and
"this is confirmed" come out of the tier map as the same word.** `ProcessSet.tier_of` combines
with `min`, and the identity of a min is the maximum element — so every state slot no Process
touches reports `VALIDATED`. That is the one mark nothing in this engine has earned, and it is
handed out automatically to exactly the quantities the model says the least about.

D-261. Building the console, the summary badges came back reading *validated* for an unhopped
beer's bitterness, for browning in a wine that never ages, and for pH in a run with no acid
mechanism live. All three are flat lines at whatever the compile seam seeded, and all three were
about to be shown to a user as *checked against real measured data*. The fix is not a
special-case list: `RunResult.touched_variables()` collects what any active mechanism actually
writes, and a readout whose sources are all outside that set reports **inert — no tier at all**,
with copy saying it held its starting value throughout, which is not the same as being confirmed.

Same family, one layer down: **`Tier.SPECULATIVE` is the enum's zero and is therefore falsy**, so
`if tier:` silently reports the least confident mark as *no mark*. Every check in `app/` is
`is None`, and a test pins the falsiness so the discipline has a stated reason rather than a
habit.

**Why:** the tier system exists to stop confidence being borrowed, and this is the one place the
engine's own arithmetic hands it out for free. It cannot be caught by reading a chart — a flat
line badged `validated` looks like a well-behaved constant, not like a hole. Any surface that
reads `tier_map` without asking what is driven will make the same claim, and the more inert the
run the louder the claim gets.

**How to apply:** before reading a tier for display, compute the set of variables the active
mechanisms write and intersect. Report the empty case as its own state with its own word, never
as a tier. Make the guard **two-sided** — assert both that the surface says "inert" *and* that
the raw engine value it is protecting against really is `VALIDATED`, so the test fails loudly if
that convention ever changes rather than passing for a reason that no longer holds. Related:
[[feedback-ask-the-engine-for-a-scope-never-re-derive-it]],
[[feedback-a-guard-can-be-blind-to-the-mutation-it-names]].
