---
name: feedback-a-gap-can-be-held-open-by-a-second-anchor
description: "D-216 - the knob that closes a measured gap was in band; what forbade it was a SECOND anchor riding the same knob, so check what else the knob holds before calling the gap a defect"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: f6587b9b-5f34-4f36-af85-d37c405530dc
  modified: 2026-08-17T19:55:05.598Z
---

**Before treating a measured gap as a defect, ask what ELSE the knob that would close it is
holding open.** D-215 measured this engine fermenting a published wort ~2.8× too slowly and
left it unattributed. D-216 found the knob immediately — the specific sugar-uptake rate — and
found the required value **inside its own printed band** (1.397 against a shipped 0.5, band
0.3-1.5). Nothing about the parameter forbade the fix. What forbade it was the *other* anchor
on the same constant: an acceptance criterion for a different wort at a different temperature,
which the same knob also sets, and which breaks at **q ≈ 0.6** — after closing under a fifth of
the gap. Searched jointly with the second-largest term, **no in-band pair satisfies both**; at
the most generous corner a **1.79×** shortfall survives.

**Why:** a gap is measured against one anchor, but a shared parameter is a claim about every
anchor at once. Score only the anchor that motivated the beat and the fix looks free — it is in
band, it is one line, the test turns green — and you discover the cost when an unrelated
benchmark goes red and gets "adjusted". The conflict is the finding, and it is a finding about
which anchor the model should be calibrated to, not about the parameter. That question belongs
to the owner: here one anchor is a published trial and the other a rule-of-thumb window from the
project brief, which the repo's own guide calls "reference, not gospel".

**How to apply:** when a knob would close a gap, enumerate every test and benchmark that reads
that knob and score them all in the same probe — not afterwards. Then make the refusal **two-tier**
if you can: mine held at "no in-band point works" *and* "not even the unbounded limit of the
dominant term reaches it" (removing catabolite repression entirely closed 79 % of the gap and
still fell short). The second tier is what kills the obvious re-proposal — *re-source that
constant and the gap goes away*. And watch for a parameter that is free for one anchor and not
the other: the activation energy here is **exactly** inert at the benchmark's temperature because
that temperature IS the reference, so it moved one anchor and not the other — a real decoupling
lever, refused only because its low band edge is a value the archive had already debunked.
Related: [[feedback-fit-the-observable-not-the-consequence]],
[[feedback-a-band-is-per-parameter-a-claim-is-joint]], [[feedback-closer-to-reality-decides]].
