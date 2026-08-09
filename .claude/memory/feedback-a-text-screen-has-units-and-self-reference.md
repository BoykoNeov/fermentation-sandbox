---
name: feedback-a-text-screen-has-units-and-self-reference
description: "Is my number in the text?" fails two ways — unit-shifted restatement hides real hits, and a note quoting its own constructed edge is anti-evidence, not evidence
metadata:
  type: feedback
---

D-167 screened 110 band edges for *"does this edge's literal value appear in the
parameter's own text?"*. The screen was wrong twice, in opposite directions, and
**hand-verification caught both — the screen never noticed either.**

**Unit shift (false NEGATIVE, and it dominated).** `E_a` bands ship in **J/mol**;
every note discusses them in **kJ/mol**. `E_a_polymerization` is `[35000, 75000]`
and its note reads *"activation energies ~35-75 kJ/mol"* — the edges are right
there, scaled by 10³, and the screen called both ABSENT. Fixing it took the count
**14 → 42**. A range dash cost one more: `literature 1.77-1.86` tokenised as
`[1.77, -1.86]`, so a genuine hit scored as a near-miss.

**Self-restatement (false POSITIVE, and it inverts the meaning).** The
`uncertainty.note` often quotes its own edge *while declaring the author invented it*:
*"banded down to ~50 % efficiency **(0.27)**"*, *"spans the realised rate **(0.3)** up
to the growth-coupled peak **(1.5)**"*. The screen scores these as sourced; the
sentence says the opposite. Six edges. **A "the note contains the edge" pin is
anti-correlated with external sourcing here.** Plus dumber positives: `Table **6**`,
the formula subscript `CH1.8O**0.5**N0.2`, `pH 5.2 model wort, **90**-130 C`.

Of 42 machine hits, **14 were rejected on reading**.

**Why:** a text screen looks like a measurement and is really a proposal. Its recall
is set by formatting conventions the data was never written to satisfy, and its
precision by the fact that prose repeats numbers for many reasons. Reporting the raw
count as the finding is [[feedback-count-and-print-your-skips]] with the denominator
inverted — here it was the *numerator* that lied, both ways.

**How to apply:** name the failure modes **in the pre-registration**, before running —
radix/format, unit shift, range punctuation, rounding — then hand-verify every hit
*and* every near-miss, and print the walk from raw hits to adjudicated ones so the
rejects are visible. Make the matcher report the *scale factor* (10³ = J→kJ is
credible, 10⁴ is coincidence) so the reading is cheap. And treat "the repo's text says
so" as [[feedback-transcribe-tables-not-prose]] demands: it licenses *"the repo asserts
X"*, never *"the source prints X"* — an upper bound, the same shape as D-166's
`reads`-reach. Related: [[feedback-pre-register-the-cheap-prediction]].
