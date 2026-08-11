---
name: feedback-a-cap-being-written-to-cannot-be-raised
description: Before raising any threshold, check whether the measured value sits exactly AT it across commits — that means the cap sets the fill rate, so raising it is futile, the overflow relocates, and the honest end state may be to delete the threshold and report the number instead
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

**The sequel (D-177) is the part to remember: a per-item cap does NOT rescue a total.** Keeping
300 as a "backstop" behind the new per-block check lasted 10 sessions. It re-pinned at **exactly
300 on 8 of the last 9 commits** while **every block sat under cap** — the file was simultaneously
the healthy shape the design described and glued to the ceiling, because a per-item cap bounds
what each new item **ADDS** and nothing bounds **how many items there are**. At ~1 record/session
× ~3.5 lines, any fixed total is reached every ~10 sessions, at any value. **Size the escape
hatch before trusting it**: the licensed retirement move covered 18 blocks / 59 of 300 lines, but
reading them, nearly all were prohibitions their corrector *sharpened* rather than replaced, so
real recovery was 5–15 lines — the hatch could not pay for the check. **A derived cap is not the
fix either** (8×blocks against a file running 3.8× can never fire: a guard that forbids nothing
reads as coverage). **The end state is to delete the threshold and REPORT the number** — with
nothing to hit, it cannot be a target; being always printed, it cannot be vacuous; and the
judgement call goes back to a person instead of being trimmed to a round number.

**How to apply:**
- **Measure before tuning.** `git log` the file, print the metric per commit. A flat line at
  the cap value is the diagnosis; a rising line indifferent to the cap is a different disease
  (real accumulation) needing a retirement policy, not a bigger number.
- **Move the pressure to shape, not size.** Cap the *per-item* footprint so the only way to
  comply is to distil the new item — eviction of an old one then cannot satisfy the check.
  Cap **blocks/paragraphs, not just list items**, or the content re-types itself one level in.
- **Enumerate every surface the budget covers** before capping one, and don't overstate what a
  cap covers. Both named gaps proved real: `CLAUDE.md` was the unmeasured third boot surface and
  had grown 66 → 138 lines (+30 the day the memory cap moved), now measured at D-177 with a cap
  sized from **its own** distribution — importing the memory file's 8 would have fired on
  legitimate documentation prose [[feedback-name-guards-for-what-they-forbid]]. `MEMORY.md` is
  capped per **row**, so row **count** is still open (5 → 40).
- **Say what the guard cannot see.** D-177's block cap would NOT have caught the +30-line day it
  was written for — that arrived as two ~10-line blocks, and catching it needs a growth rate a
  PostToolUse hook has no history for [[feedback-conceded-caveats-are-not-coverage]].
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
