---
name: feedback-check-the-blocker-is-still-blocking
description: "A record's premise can be demolished by unrelated work, which leaves NO correction marker - probe the code before treating an old item as blocked"
metadata:
  node_type: memory
  type: feedback
---

When picking a beat off an old open item, the **first probe is not "how do I build this" but
"is the reason this was closed still true?"** — run it against the **code**, not the record.

**Why:** D-108 named its own blocker in plain words (closure O₂ ingress) and D-136 built it 28
records later **without naming D-108**, because D-136 was about oxygen and not about
acetaldehyde. No `Corrects:` marker was owed, so **D-108 grew no ⚠** and the item rode the
copy-forward open lists **20 more times after its blocker shipped** — still described as out of
scope for a reason that had stopped being true. The correction map has an edge for a later beat
that **disagrees** with an earlier one, and **no edge at all** for a later beat that **changes
the world the earlier one was describing**. That case is the commoner one, because it needs no
disagreement and nobody involved is wrong at the time. Sibling of
[[feedback-a-refusal-reads-like-a-gap]]: there a REFUSED item read as open, here a genuinely
blocked one stayed listed after its blocker was removed. Same hole — **an open list is a claim
about the present, stored as prose written in the past.**

**A third shape, and the worst one: the stale blocker living in CODE (D-254).** The reason two
fusel fixtures could not be scored against their own papers was written as a comment inside the
fixtures — "the composition is in the paper and not in this repo … closing the gap needs the
paper". D-246 put both musts in the repo eight records later and neither comment moved. A comment
is worse than an open list because it is read *at the moment someone is deciding*, it carries no
date, and it reads as a fact about the world rather than as a note from a past beat. **Two more
things that shape found.** The comment was duplicated **verbatim in two files** and only one copy
had expired, so "fix the comment" was not a find-and-replace but a per-site re-check — and the
live one (Rollero) narrowed from "not in this repo" to a single answerable question once the dead
one was understood. And a *different* comment in the same beat, on a clamp, said it "never binds
in practice" while the record that shipped it had measured it binding in the same document. So
run the probe against code comments too, and when one expires, grep the repo for its twin.

**How to apply:** One probe, before any design. At D-188 it was a sealed sulfited bottle read at
the end: non-zero ⇒ the beat was a **measurement**, not a build, and the whole pitch to the owner
changed before a line was written. Cheap, and it is what makes "still blocked" a finding rather
than an inheritance. Pairs with [[feedback-verify-latest-state-not-breadcrumbs]] (which is about
*later* records burning down a "Next:" list) — this is the inverse direction, an *earlier*
record's premise quietly expiring. See `.claude/memory/prohibitions/acetaldehyde-ladder.md`.
