---
name: feedback-a-doc-rots-where-it-duplicates
description: "A doc decays exactly where it restates another doc's content; and when one goes stale, suspect the maintenance rule before the discipline"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 246a44fa-2457-4d66-9173-cbb8f0ece5c9
  modified: 2026-08-11T22:10:45.144Z
---

`ARCHITECTURE.md` went **127 `DECISIONS.md` commits** stale (D-184). But not uniformly: its
*structural* sections — layering, `Process`/`RateModifier`, tiers, units — were still accurate
after four months, while every section that restated the archive's reasoning had rotted. The
`docs/plans/milestone-*.md` files, which duplicate the archive **by construction**, had rotted
completely — one still offered shipped work as "deferred", which had already caused two
re-proposals ([[feedback-verify-latest-state-not-breadcrumbs]]).

**Why:** a document that is the *sole owner* of its content stays true for free — nothing else can
contradict it. One that restates another surface decays at the rate that other surface advances,
and the decay is invisible because both copies look authoritative. So the diagnostic question is
not "is this doc old?" but **"does anything else own what it says?"**

The second half is about the guardrail. The rule meant to prevent this — "keep the plans updated as
work progresses" — was *followed for a while and then silently stopped*. Re-affirming it would have
bought one more cycle of the same decay. **When a doc goes stale, check the rule before blaming the
discipline**: a rule that depends on remembering, with nothing measuring it, has already failed once
by the time you notice. Cf. [[feedback-a-cap-being-written-to-cannot-be-raised]] — a prose rule is
not a mechanism.

**How to apply:** when writing or repairing a doc, first ask what surface *owns* each claim, and cut
anything a lower/other surface owns down to a pointer. Give a doc derived numbers plus a runnable
snippet that re-derives them, and state which side wins on disagreement (the code). Prefer deleting
a failed maintenance rule over restating it, and prefer bannering a superseded doc as history over
half-updating it — a half-maintained page is worse than a frozen one, because its unstruck items
read as open work. Verify the snippet actually runs before shipping it: D-184's first draft ended in
`bc`, which is not installed on this box, and a self-check that cannot execute reads as verification
([[feedback-grep-finds-claims-not-guards]]).
