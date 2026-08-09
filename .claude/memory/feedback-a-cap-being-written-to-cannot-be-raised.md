---
name: feedback-a-cap-being-written-to-cannot-be-raised
description: Before raising any threshold, check whether the measured value sits exactly AT it across commits — that means the cap sets the fill rate, so raising it is futile and the overflow just relocates
metadata:
  type: feedback
---

**Before raising a cap, plot the measured value's history. If it sits *exactly* at the cap
across many commits, the cap is a TARGET, not a limit — it sets the fill rate, and raising it
is provably futile.** And when a cap covers one of several equivalent surfaces, the overflow
does not disappear; it **relocates to the uncapped sibling**.

D-169: the project-memory line cap was raised four times (150 → 200 → 250 → 300), each time
because it was evicting live prohibitions. The git history settled it — the file sat at
**exactly 250 for 13 consecutive commits across 12 days**, then took **47 of the next 50 lines
the day the cap moved to 300**. Content does not land on a round number 13 times. Meanwhile the
*uncapped* `MEMORY.md` index row for that same file had grown to **950 chars** against a
211-char median: content squeezed out of the capped file had reappeared next door.

**Why:** a total that is the *only* check can be satisfied in exactly one way — eviction — so
each raise reproduces the harm it was meant to undo. The number was never the defect; being a
single poolable budget was. Raising it buys one beat and re-arms the same failure.

**How to apply:**
- **Measure before tuning.** `git log` the file, print the metric per commit. A flat line at
  the cap value is the diagnosis; a rising line indifferent to the cap is a different disease
  (real accumulation) needing a retirement policy, not a bigger number.
- **Move the pressure to shape, not size.** Cap the *per-item* footprint so the only way to
  comply is to distil the new item — eviction of an old one then cannot satisfy the check.
  Cap **blocks/paragraphs, not just list items**, or the content re-types itself one level in.
- **Enumerate every surface the budget covers** before capping one, and don't overstate
  what a cap covers: `CLAUDE.md` is a boot surface and is unmeasured, and `MEMORY.md` is
  capped per **row**, so row **count** is still an open channel.
- **Run the criterion that licenses the cut.** "It's held losslessly elsewhere" justified
  every measurement dropped and was asserted three times before being checked
  [[feedback-conceded-caveats-are-not-coverage]]. 20 of 20 sampled figures were present —
  but one only matched in the archive's notation (`[35000, 75000]`, not the memory's
  `[35k,75k]`), so grep the target's units, not your own
  [[feedback-a-text-screen-has-units-and-self-reference]].
- **Distrust a metric that a stronger-looking signal explains.** A digit-density check looked
  clean at 217-vs-65 and was rejected: normalised per line the "evidence-dense" block ran ~10
  digits/line against a short prohibition's ~20, so it was reading item SIZE
  [[feedback-name-the-field-your-predicate-read]].
- **Test the checker's logic, never the live file's compliance** — the latter converts a
  warning into enforcement [[feedback-name-guards-for-what-they-forbid]].
