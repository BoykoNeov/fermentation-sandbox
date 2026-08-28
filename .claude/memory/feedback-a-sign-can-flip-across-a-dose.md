---
name: feedback-a-sign-can-flip-across-a-dose
description: "Two errors in one observable can REVERSE its sign across a knob that scales only one of them; sweep that knob to separate them"
metadata:
  node_type: memory
  type: feedback
---

D-250 measured whether yeast nitrogen uptake starves a co-inoculated malolactic bacterium. On the
2 g/L amino-acid must the starved arm converted **more** malate (0.5526 vs 0.6478 left) — the
opposite of starvation. On a 0.5 g/L must the same comparison left **4.0x more** malate. Same
model, same Processes, opposite sign.

**Because two errors shared the observable and only one scaled with the dose.** Losing the
bacterium's nitrogen costs conversion at both doses. The second error — a pH excursion from
mis-booked ammonium — was +0.045 at 0.5 g/L and +0.215 at 2.0, and it feeds MLF's own pH logistic
gate. At the high dose it over-compensated and reversed the reading.

**The decomposition is free when one knob scales one error.** No new machinery: run the existing
probe at two doses. The confound is a *function* of the dose, the effect under test is not, so the
sweep separates them where a single point cannot.

**And check which way a confound points before dismissing an attribution.** Brett's gate partitions
molecular SO2 at the solved pH, so the same pH rise *helped* Brett. Its 96 % growth loss is
attributable **despite** a tailwind, which is stronger than a null confound.

**How to apply:** when a reading contradicts the mechanism you expect, do not report either the
number or the reversal until you have named every path from your change to that observable. If one
of them scales with a knob and the other does not, sweep it. Reporting the single point would have
published the wrong sign. Related: [[feedback-a-hit-can-be-two-errors-cancelling]],
[[feedback-one-knob-decomposition-beats-an-unbuildable-control]].
