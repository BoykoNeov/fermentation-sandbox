---
name: feedback-a-majority-is-not-a-direction
description: "A 6-of-8 majority in the predicted direction flipped to 10-of-24 when the seeds tripled; fix the stopping rule and both outcomes before the run lands, and label rows with the run that produced them"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 95d03d1d-9761-4bf2-a538-1fd5388b2bbc
  modified: 2026-07-29T01:12:10.227Z
---

**A majority across seeds is not a direction until it is compared against the
statistic's own seed-to-seed range.** Measured, D-161:

| seeds | widened | median POST/PRE | reading it invites |
|---|---|---|---|
| 8 | **6 / 8** | **1.036** | "widens modestly, as predicted" |
| 24 | **10 / 24** | **0.970** | no detectable direction |

Same arms, same members, same design — only more seeds. The effect the archive
*expected* (D-160: "the expectation is widening where the output is pH-sensitive")
went from a 3:1 majority to a coin flip, and the median crossed to the other side
of 1. The tell was available at 8 seeds and I nearly missed it: the per-seed ratio
already ranged **0.957–1.348**, so an effect of 1.036 was a fifth of its own noise.

**Why:** an 8-point sign test has almost no power (6/8 is p = 0.29 — nowhere near
evidence), but a table of eight rows *looks* like a measurement, and a majority in
the direction you predicted is the easiest thing in the world to write up as
confirmation. That is fitting the story to the number.

**How to apply:**

- **Write the stopping rule and BOTH outcomes down before the extra run lands.**
  "If it stays ~16/24 the answer is a small effect below significance; if it moves
  toward 12/24 the answer is no detectable direction" — committed in advance, either
  result is reportable and neither is a save.
- **Put the effect next to the noise, in the same units.** If the per-seed spread
  of the statistic straddles the effect, say so in the record instead of quoting a
  median. Same discipline as [[feedback-pin-tolerance-vs-solver-tolerance]]: measure
  the floor before you claim to be above it.
- **Report the p-value even when it is embarrassing**, and flag multiplicity: 24
  uncorrected per-slot sign tests produced two nominal p < 0.05 rows, both with a
  median ratio of exactly 1.000 — sign flips on float noise, which is what chance
  produces at that count.
- **Do not extend again to chase p < 0.05.** Once the effect is small relative to
  its own noise, more seeds buy precision on a question already answered.

**And label rows with the run that produced them.** D-159 printed
"13.99 / 16.69 / 20.09 % across seeds 0/1/2"; those are the values **sorted
ascending**. The real order is 16.69 / 20.09 / 13.99. Nothing in the archive could
catch it — **re-running D-159's own harness verbatim was the only thing that did**,
and it is the cheapest check available when a record hands you a numbered row.
Related: [[feedback-transcribe-tables-not-prose]] (the same failure on the reading
side), [[feedback-count-and-print-your-skips]], and
[[feedback-a-null-result-needs-a-positive-control]] — this null carries its
positive control (`only=['mu_max']` moves the target 64 %) and its negative
(`ethyl_acetate` flat at 4.8e-09).
