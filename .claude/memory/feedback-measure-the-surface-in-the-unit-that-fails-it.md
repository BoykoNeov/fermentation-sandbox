---
name: feedback-measure-the-surface-in-the-unit-that-fails-it
description: A monitor that reports a different unit from the one the consumer enforces is blind by construction — it can run clean for months while the thing it watches is already failing
metadata:
  node_type: memory
  type: feedback
---

`check_memory_size.py` exists to keep the session-boot surfaces visible, runs on every write to
them, and reported `MEMORY.md` as "105 lines, 21 blocks, median 1, max 4" — clean, for months.
Meanwhile the harness had been **silently truncating that file at load** because it loads on
**bytes** (25910 against ~24.4KB) and every number the hook printed was a **line** count. The
one mechanism watching the surface could not see the only failure that surface has. Its per-row
cap had the same fault: 320 **chars** on a file full of `—`, `→`, `×`, `₂`.

**Why:** a monitor is usually written against the shape the *author* worries about (bloat,
sprawl, changelog drift) rather than the constraint the *consumer* actually enforces. Those
coincide often enough that the mismatch never shows — until the surface approaches the real
limit, at which point the monitor's silence reads as a pass. This is worse than no monitor:
the hook's clean report is *why* nobody looked, and the failure is silent on both sides.
Note the failure was NOT a discipline failure — the rows were being written carefully the
whole time; the wrong unit is a design defect in the check.

**How to apply:** for any guard on a resource, name the consumer and the unit **it** rejects
in — bytes for a loader, tokens for a context window, wall-clock for a timeout, rows for a
quota — and report that number with its headroom, not a proxy that correlates with it. If the
limit is external, say in the comment that the value is *derived* rather than measured (this
one is `24.4 × 1024` off a rounded display; two of its digits are not real). Report headroom
always; make it a finding only when actually exceeded, so it cannot become the round-number
target [[feedback-a-cap-being-written-to-cannot-be-raised]] diagnoses.

**Corollary, and it is the half that nearly shipped broken.** The structural fix here moved 91
rows out of `MEMORY.md` into `lessons/`, and left the *write* rule as prose in three places —
the index header, the batch-end ritual, and the harness's own memory instructions. **The third
one is outside the repo and says the opposite** ("add a one-line pointer in `MEMORY.md`"), so
a session following it verbatim reverts the split. Before trusting a rule that lives in more
than one surface, enumerate the surfaces and check they *agree*; where one is not yours to
edit, the rule needs a mechanism, not a fourth sentence. Same disease as
[[feedback-a-doc-rots-where-it-duplicates]], with the twist that the rotted copy is upstream.
Related: [[feedback-name-the-field-your-predicate-read]] (which field a predicate read),
[[feedback-conceded-caveats-are-not-coverage]] — the docstring had already conceded that this
channel was uncovered, and a concession is not a check.
