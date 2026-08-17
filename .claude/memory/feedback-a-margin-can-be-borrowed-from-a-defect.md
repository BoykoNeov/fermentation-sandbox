---
name: feedback-a-margin-can-be-borrowed-from-a-defect
description: "Fixing a defect can degrade a neighbouring 'win' whose margin was being held open by that same defect — re-read what the win was scored against, and pin the shortfall it was hiding"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: dcd1ccd2-410e-4bfa-9d67-bd5e911c9dec
  modified: 2026-08-17T14:46:38.029Z
---

When a beat corrects a defect, expect a **neighbouring settled result to get worse**, and go
looking for it rather than waiting for a test to find it. A margin measured while a defect was
live may have been **borrowed from that defect**.

**The case (D-211).** D-183 retired a flux-linked acetic-acid yield in favour of a growth-linked
producer, pinning the margin at 32.5 ppm RMSE against 61.6 — ratio 0.528, guarded as `< 0.6`.
D-211 corrected beer's growth rate (2.88× too fast) and that test went red at **40.7 vs 65.3**,
ratio 0.624. The retirement itself survived, but its *evidence* was partly artificial: when
D-183 measured it, growth finished inside ~20 h, so a growth-linked producer delivered ~100 % of
its acid by day 1 and **appeared** to explain Tyrell's early rise. The real crop is ~36 % grown
by day 1. Acetic therefore rises faster than growth **and** faster than flux — the model books
0.360 of its rise by day 1 against a measured 0.773 — so neither driver in the file explains it.
The old defect had been concealing a 2.15× shortfall in the very quantity the choice was made on.

**Why:** a comparative claim ("A fits better than B") is scored *through* the model, so every
defect upstream of both is in both numbers. A defect that happens to flatter A makes A's win
look like evidence for A's mechanism when it is partly evidence for the defect. The failing test
is the good outcome here — but the temptation at that moment is to read it as "my change broke
something" and relax the threshold, which deletes the finding.

**How to apply:** when a beat changes a shared upstream quantity, grep the archive for results
**scored through** it, not just code that reads it. When such a test goes red, ask first whether
the old margin was borrowed. Then: keep the claim the test actually protects (here, the
ordering), re-pin the moved number **with both the old and new values in the message**, and add
a **new assert for the shortfall the correction exposed** — otherwise the beat quietly converts
a discovered defect into a loosened threshold ([[feedback-a-cap-being-written-to-cannot-be-raised]]).
A `Flags:` marker on the older record, not a `Corrects:`, is usually right: the verdict stands,
the reasoning behind it does not.
