---
name: feedback-a-clamp-sits-between-designed-and-realised
description: "An estimator that credits a parameter's DESIGNED share is wrong wherever a clamp truncates the branch - and the error can be invisible to every guard that reads a concentration"
metadata:
  node_type: memory
  type: feedback
---

When a helper attributes an outcome by multiplying a **declared share parameter** by an amount
consumed, it is asserting that nothing sits between the design and what the run actually did.
Check that. A `min()`, a headroom clamp, a `max(0, …)` or any saturation between the two makes the
helper a **biased estimator**, and if the clamped quantity is part of a conserved partition the
bias is not a loss — it is a *transfer*, so one output reads high by exactly what another reads
low. Correcting it can therefore move a verdict from one alcohol/species/pool to another rather
than making it go away.

**Why:** at D-254 `_amino_acid_share` credited each valine branch at its designed share (0.23
secondary / 0.15 primary). `ehrlich_draws` had been truncating the secondary branch at `headroom`
since D-111 — and **D-111's own record measured it doing so** ("realised 0.2233 against a sourced
0.23, isobutanol 0.1567 against 0.15"). The helper, written ~130 records later, re-assumed what
that record had already falsified. Because the non-Ehrlich lump scales against the Ehrlich draw,
the total was pinned, so the 4.5 % taken off isoamyl reappeared on isobutanol: the corrected
measurement moved an over-attribution from one alcohol to the other instead of removing it.

**The second half is why nobody caught it.** The clamp was **invisible in the state** — the
precursor was supply-limited and each product's amount was set by its own rate law, so deleting
the whole route moved every concentration by ≤3.4e-6 while moving provenance by 4.5 %. Every guard
in the repo read concentrations. A quantity no test can see will not be found by adding more of
the same tests, and its code comment said it "never binds in practice" for eleven years of records.

**How to apply:** when writing or reading an attribution helper, grep the path from the parameter
to the applied value for anything that can truncate, and if one exists either integrate the
applied draw or state the bias with its direction. Ask whether the truncated quantity is part of a
partition whose total is pinned — if it is, expect an equal-and-opposite error elsewhere and
measure both legs, never just the one that helps. And when a clamp is documented as never firing,
that sentence is a claim needing a receipt, not a reason to skip the check.

Related: [[feedback-a-toggle-measures-nothing-on-an-exhausted-pool]],
[[feedback-check-the-blocker-is-still-blocking]], [[feedback-a-ratio-guard-can-pass-on-an-overproduction]].
