---
name: feedback-reproduce-a-published-number-before-trusting-the-new-column
description: "An A/B harness over two source trees must reproduce a number the archive already published before its new column is trusted — mine silently ran one tree's code twice"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ae98d8d8-bbb4-4e60-be15-48e20c8776e7
  modified: 2026-08-25T19:40:17.931Z
---

An A/B measurement across two versions of the code must first reproduce a number the **old**
version already published. Only then is its new column evidence.

**Why:** at D-227 I built a harness to compare the shipped tree against HEAD and it printed a
clean, plausible table — every arm ~1.24, the "before" and "after" columns clearly different. It
was wrong. `sys.path.insert` cannot *un*-import a package: the first `import fermentation` won,
both columns ran that same code, and only the parameter files differed. The output looked exactly
like a real result. What caught it was not inspection but arithmetic — the "before" column did not
reproduce D-226's own published 0.9549 and 13.9 %, which that record had measured on that exact
tree. Re-run with one **subprocess per tree**, both numbers landed on the nose, and the new column
became trustworthy in the same instant.

**How to apply:** before reading an A/B table, find a value the archive already states for the
baseline arm and check the harness returns it. If the baseline column has no published anchor,
make one — run a number the old record reports, even if it is not the quantity you care about.
Then, for anything comparing two working trees: separate **processes**, never separate `sys.path`
entries, and have each process print the tree it loaded. The same applies to a tree's parameter
data — pair each tree with its *own* data dir, since crossing them is the mistake that produced
the fake table above. Related: [[feedback-verify-the-mutation-applied-not-just-the-restore]] and
[[feedback-a-baseline-log-goes-stale-between-edits]] — this is the third way a baseline can lie,
and the only one where the harness itself is the liar.
