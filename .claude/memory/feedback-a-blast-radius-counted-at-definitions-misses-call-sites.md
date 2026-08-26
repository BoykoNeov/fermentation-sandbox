---
name: feedback-a-blast-radius-counted-at-definitions-misses-call-sites
description: "Pricing a signature change by grepping definitions missed 28 call sites; count where the break lands, not where the change is written"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 217f4393-f3f4-48e3-9bef-4bb21f5faae0
  modified: 2026-08-26T13:36:15.748Z
---

When you price the cost of changing a **signature**, count the sites that will *break*, not the
sites you will *edit*. Those are different sets and the second is the smaller one.

D-234 priced widening `StateMutation` by one argument as: 13 `def mutate` in `src`, 13 `mutate=`
sites, 1 in `tests`, 3 type references. All four numbers were about where the callable is
*defined or handed over*. Shipping it at D-235, the count that actually cost the work was
**28 direct `event.mutate(schema, y)` call sites across six test files** — every one a test
applying a verb's jump in isolation, every one a hard failure under the new arity, and none
of them matched by either grep. `mypy` found them. The type-reference count was also wrong (2,
not 3), which is the cheap tell that nobody re-ran the search.

**Why:** the same instrument failure D-234 itself had just diagnosed one level up — it found
21 census members invisible to a `parameters["…"]` grep because the *reads* happened inside a
helper, not at the pattern. Here a `def mutate` grep missed the *calls* for the same reason.
Both times the search was written for the thing being changed rather than for the thing that
would break, and both times the missed set was larger than the found one.

**How to apply:** before quoting a blast radius, ask what a *consumer* of this thing looks
like, and search for that shape too — `\.name\(` for a call, not just `def name`. Better, let
the type checker or the compiler enumerate it: `mypy` / a full build gives the exact set for
free and cannot have a blind spot the way a regex does. If you must quote a grep-derived
number, quote it as a **floor** and say which instrument produced it
([[feedback-size-a-set-with-the-instrument-that-can-see-it]],
[[feedback-grep-finds-claims-not-guards]]). A priced radius that a later beat has to correct
is worse than an unpriced one, because it reads as measured.
