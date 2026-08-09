---
name: feedback-count-and-print-your-skips
description: A measurement harness that silently drops inputs it cannot parse reports a clean result it never measured; count and print every skip
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 390b406d-8aa0-41a7-bf5d-21877a7769f1
  modified: 2026-08-09T21:25:53.742Z
---

When a harness filters, parses, or evaluates its way to a population, **count and print
what it dropped**. A `continue` with no counter turns "I could not decide this one" into
"this one is fine", and the summary line reports a denominator that was never measured.

**Why:** D-157's sweep reported its decidable class **5 of 5 clean**. It wasn't. `literal()`
evaluated assertion bounds with `eval`, so any bound that was not a bare numeric constant
raised, returned `None`, and the whole assertion was `continue`d — uncounted. Adding
module-constant resolution and a skip counter took the population 5 → 7 and exposed a live
defect (`ethyl_acetate_eq`, 0.99% of draws escaping the interval its own test asserts). The
same harness also produced the opposite failure — an ignored comparison-**operator
direction** turned `x > 0` into an upper bound of 0 and confidently reported a 100% breach.
Both were found by making the harness *print*, not by reasoning about it.

**How to apply:** Every `continue`/`except: pass`/filter in a measurement script appends to
a `SKIPPED` list that gets printed with the results. State the denominator as
"N decided, M skipped", never as "N, all clean". A skipped item is an **undecided** item,
and undecided is not a synonym for safe — false negatives are the self-sealing direction
[[feedback-check-the-schema-not-the-caller]], and naming a gap in prose does not discharge
it [[feedback-conceded-caveats-are-not-coverage]]. Related:
[[feedback-mutate-the-premise-before-building-the-guard]],
[[feedback-rejected-values-must-be-unreachable]].

**Amendment (D-163): printing the ledger is necessary, not sufficient — reconcile the
headline against it.** The band-edge campaign *did* print its skip ledger, and the prose two
paragraphs above it still read "all 678 edges moved, in both directions". Both came from the
same run. The trap is that `failed_to_mutate == 0` and `edges_moved == 678` are **different
claims**: the first says every edge the operator *targeted* landed, the second says every
edge *exists* in the targeted set. Six half-edges were never targeted at all (value sitting
on its own band edge, so a rescale-about-the-nominal moves nothing), making the real
denominator 672. So: after a run, take each headline number and re-derive it **from the
ledger**, not from the intent of the code that produced it. A ledger nobody subtracts from
the headline is decoration. The blind spot is also worth naming rather than just counting —
here those six turned out to be exactly the six the archive had already enumerated at D-153
for an unrelated reason, which made the gap a *characterised* bound instead of a loose end.

**Amendment (D-170): a RESOLUTION rule drops things too, and it leaves no `continue` to
count.** The ordering screen matched a claim's sibling by walking candidate names
longest-first and `break`ing on the first hit — so a sentence naming two parameters silently
resolved to the **longer string**, and the other claim was never emitted at all. There was no
skip to print, because from the harness's point of view nothing was skipped. It dropped
**D-168 §3's own ordering 2** that way, and the mis-resolved row was then hand-killed as a
false positive, which hid the loss twice. Removing the `break` took 49 → 69 candidates and
24 → 36 evaluable, and exposed a note that was **false at the shipped nominals**.
**So: measure RECALL against a known-true set, not just skips.** Any screen re-covering
ground a previous record worked by hand should re-find that record's findings *with the
correct binding*, and print the ambiguous-input count (here: 21 windows naming >1 sibling —
the denominator nobody had). Until that control exists, every count off the screen is a
**lower bound**, and must be written as one — calling it a census is the D-163 §7 error
arriving from the recall side. Pairs with
[[feedback-a-null-result-needs-a-positive-control]] (same discipline, other direction) and
[[feedback-a-screen-is-not-idempotent-under-its-own-repair]].
