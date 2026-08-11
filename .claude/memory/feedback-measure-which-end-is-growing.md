---
name: feedback-measure-which-end-is-growing
description: "Before bounding a growing surface, measure which END of it grows — and if that end is incompressible, the lever is granularity, not eviction"
metadata:
  node_type: memory
  type: feedback
---

The project status memory grew 300 → 406 lines with **zero** cap findings. Splitting the growth by
age settled it (D-185): entries about the newest five decisions have held at **6-17 for two months**,
while entries citing records 40+ decisions back went **0 → 22, monotone**. Two ends behaving
oppositely — and *every* mechanism tried across D-169, D-177 and D-185 had been aimed at their
**sum**, which is why four of them failed in four different ways (a total gets written up to; a
derived total can never fire; a fifth raise is what the record forbids; a one-off distillation
regrew 114×).

**Why:** a whole-surface number averages over parts with different dynamics, so it cannot tell you
where to push. A guard aimed at a sum lands on whatever is being edited — which is the live work,
the one end that was already fine.

Then the second half, which killed the obvious fix. I designed a budget on the growing end and it
had every property the earlier candidates lacked: it could fire, it had no round number, it never
touched live work. It died on **reading the 22 entries it would act on** — all 22 are live
prohibitions (the O₂ gate is inverted; direction is the owner's call; do not rebuild bottle
reduction). **Age does not measure finishedness**, and by survivorship the oldest guardrails are
the ones that keep earning their place, so retiring by age deletes by exactly the wrong key. Cf.
[[feedback-a-cap-being-written-to-cannot-be-raised]].

**How to apply:** partition the surface by whatever dimension plausibly separates its dynamics (age
of the thing referenced, subject, author) and print the series per partition across git history
before designing any bound — the table is cheap and it is what makes the design decidable
([[feedback-pre-register-the-cheap-prediction]]). Predict the shape first and be willing to lose:
mine ("old entries collapse as they settle") was false — mean size was flat at 4.4/4.0/4.0/3.0 lines
and I had generalised from one collapsed bullet I happened to read. When the growing end turns out
to be **incompressible**, stop looking for something to evict: the remaining lever is
**granularity** — stop loading all of it to use one part — and the cost then becomes conditional
rather than resident, which is a different win from making it smaller. Say plainly that the rate is
unchanged ([[feedback-conceded-caveats-are-not-coverage]]).
