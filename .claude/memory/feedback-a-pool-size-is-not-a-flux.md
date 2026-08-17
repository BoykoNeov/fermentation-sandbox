---
name: feedback-a-pool-size-is-not-a-flux
description: "Killing a candidate because an intermediate's POOL is too small is unsound - measure the throughput; the right cause matters because a wrong one gets overturned"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 319bf67b-4e83-4c11-9ffb-068bbcb68911
  modified: 2026-08-17T15:40:09.839Z
---

I refused a candidate mechanism (acetate from acetaldehyde via ALD) by noting the modelled
**acetaldehyde pool was 32.2 mg/L** at day 1 while the gap needed **~47.7 mg/L** routed away —
"it dies on mass". **That argument is unsound.** Acetaldehyde is an intermediate with turnover;
the pool is a *level*, and what the mechanism draws on is a *flux*. Measured properly: 3.178 g/L
of ethanol formed over day 0-1, so throughput through that node is **~3038 mg/L** and the gap
needs **1.57 %** of it. The flux is ample. The candidate still dies — on **identifiability**
(nothing sources what fraction the enzyme diverts) and on **shape** (the pool peaks day 3, not
day 1) — but for entirely different reasons.

**Why:** a refusal record's value is the *enumeration of causes* — each candidate dying of a
named, different thing is what stops a later beat re-proposing it. A cause that is wrong is worse
than no cause: it is the thing a later beat overturns, and overturning it makes the whole refusal
look reopenable when only one bullet was bad. A standing pool says nothing about how much passes
through a node per unit time; for any intermediate on a main pathway the throughput is typically
orders of magnitude larger than the pool.

**How to apply:** never price a candidate against a **standing concentration** when the mechanism
consumes a **rate**. Integrate the flux — often available for free from a downstream product
(here, ethanol formed × the stoichiometric ratio). Then state the *fraction of throughput*
required: a mechanism needing 1.5 % of a flux is cheap and dies on sourcing, one needing 150 %
dies on mass, and those are different records. Same discipline as
[[feedback-check-the-schema-not-the-caller]] — verify the quantity you are actually claiming
about, and cf. [[feedback-a-margin-is-a-claim-about-what-holds-it-open]].
