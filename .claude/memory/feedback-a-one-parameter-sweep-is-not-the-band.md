---
name: feedback-a-one-parameter-sweep-is-not-the-band
description: "A ratio measured by sweeping ONE parameter with the rest nominal is that parameter's contribution, never a property of the reported band"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e703dc5b-abaf-4a2b-a0de-748460abea2f
  modified: 2026-08-26T11:17:18.984Z
---

Sweeping one parameter across its band with every other parameter held at nominal measures **that
parameter's own contribution**. It does **not** measure what the reported ensemble band does, and
the two can differ by two orders of magnitude — because the other parameters dominate the spread
and a shared error largely cancels across members.

**Why:** at D-233 I swept `pKa_peptide_buffer` alone and found the shipped day-14 pH span was
**1.287×** the coherent one. I wrote that into a findings file as *"the reported beer pH band is
~29 % wider"* — a claim about the band produced by all **83** sampled parameters, from a
measurement of one. Run properly over the full sampled set the band moves **1.008× — 0.8 %**. The
headline was wrong by ~36×, and it was wrong in the direction that flattered the beat I was about
to propose. The real finding was a different shape entirely: a **per-member trajectory error**
(worst member 0.0346 pH at day 14) with an essentially unchanged band, so the case for the fix had
to rest on the t=0 contract instead of on the spread.

**How to apply:** before quoting a ratio as a property of "the band", ask **how many parameters
varied when you measured it**. If the answer is one, the number is scoped to that parameter and
must be phrased that way — "this parameter's own contribution is inflated 1.287×", never "the band
is 29 % wide". To make a claim about the band, vary everything the sampler varies. And prefer the
observable where the artefact has **no legitimate component to net out**: at t=0 the state IS `y0`
and the anchor guarantees the answer, so any spread there is 100 % defect and needs no inflation
arithmetic at all — which is both the airtight claim and the one that survives every scoping
objection. This is [[feedback-a-magnitude-defence-can-be-scoped-to-the-smallest-member]] run in
reverse: that row is about defending with the smallest member, this one about *accusing* with the
loudest. Record the retraction rather than quietly deleting the number
[[feedback-a-band-is-per-parameter-a-claim-is-joint]].
