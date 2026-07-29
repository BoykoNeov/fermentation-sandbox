---
name: feedback-pin-the-band-not-the-nominal
description: "A parameter's nominal being pinned to its source says nothing about its band edges; test the edges separately, and record which are PRINTED vs CONSTRUCTED"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: f80af2b2-1bcc-4313-a84f-86ba7a57c1e8
  modified: 2026-07-29T02:01:23.736Z
---

A guard that pins a parameter's **nominal** against its source leaves the **uncertainty band
edges** completely unguarded — and a generic "the band is ordered and non-degenerate" test
does not close that gap. Pin the edges to their published numbers too, through the same
conversion, and annotate each edge as **PRINTED** (a number the source states) or
**CONSTRUCTED** (something the file built from the source's numbers).

**Why:** D-162 measured this on `closure.yaml`. Every `otr_*` nominal had been pinned to
Lopes 2007's Table I since D-136, and `test_every_otr_is_speculative_and_banded` checked
`low <= value <= high` and `low < high`. So the bands *looked* covered. A mutation silently
replacing `otr_technical_cork`'s high edge — P1's **vertical** 0.9 µL/day — with the
**horizontal** 0.4 passed all **1460** tests. The edge was load-bearing for a whole archive
argument (the closure-ordering scoping) and nothing asserted it. This is the
point-vs-band shape [[feedback-validate-calibrations-in-the-frame-that-binds]] one level
further out: not "a constraint checked at a point where the sampler reads a band," but
"provenance checked on the nominal where the band carries its own, different provenance."

The PRINTED/CONSTRUCTED distinction matters because it is invisible once the numbers are in
YAML: `otr_technical_cork`'s band spans **two storage orientations** (its high edge is
vertical) while `otr_screwcap`'s spans **one** (Table I prints `--` for screwcap vertical —
never measured), and `otr_synthetic_nomacorc`'s low edge is the file's own extension, not a
P1 number at all. Bands built to different scopes are **not commensurable**, and comparing
them on joint draws manufactures apparent disagreement: the 40.7 % inversion fell to 4.7 %
measured like-for-like. That is a reason to *state the scope*, never to narrow a band that
honestly spans what the source measured.

**How to apply:** When adding or auditing a provenance parameter, treat `value`, `low` and
`high` as three separately-sourced quantities. Before trusting that a band is guarded, run
the mutation — move an edge and see whether the suite goes red
[[feedback-mutate-the-premise-before-building-the-guard]] — with one arm designed to stay
green [[feedback-verify-the-restore-between-mutation-arms]]. D-162's "Next" names the open
sweep: band edges elsewhere under `parameters/data/` are plausibly as unpinned as these were.
