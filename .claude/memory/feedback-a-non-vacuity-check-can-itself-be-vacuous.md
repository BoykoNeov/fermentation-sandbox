---
name: feedback-a-non-vacuity-check-can-itself-be-vacuous
description: The assertion added to prove two things really differ can compare a projection on which they never differ — assert on membership, not on the filtered view
metadata:
  type: feedback
---

D-237 asserts the compile-time census is the same under all three oxidative wirings. An equality
is satisfiable by a harness that changed nothing, so it was given a companion assertion: the three
wirings must really produce different Process sets.

Written the obvious way — comparing `process_set.active` — that companion **passed on all four
scenarios while comparing nothing**. Every oxidative Process is disabled until `begin_aging`, so at
compile time the three wirings have identical *active* sets by construction. The guard added to
prevent vacuity was itself vacuous, and only fired once it compared *membership*
(`enabled_snapshot()`) instead.

**Why:** a non-vacuity check feels like the end of the audit, so nobody audits it. It is ordinary
code and fails the ordinary way — it read a filtered view of the object, and the filter removed
exactly the difference it was looking for.

**How to apply:** Make the non-vacuity assertion FAIL once, deliberately, before trusting it — the
same discipline as a mutation arm, applied to the control rather than the claim. Prefer the widest
representation of the thing (all members, not the enabled ones; the raw set, not a projection), and
when a legitimate identity makes one arm genuinely equal, exempt it **by name** rather than
weakening the check ([[feedback-count-and-print-your-skips]]).

Related: [[feedback-a-null-result-needs-a-positive-control]], [[feedback-grep-finds-claims-not-guards]].
