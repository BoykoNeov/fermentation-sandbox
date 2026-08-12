---
name: feedback-a-hit-can-be-two-errors-cancelling
description: "A pre-registered prediction that lands is only evidence if its INPUTS were right - D-200's P1 hit from an 8.2x wrong anchor times a 10x wrong concentration"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ddd93c4d-92cd-49f0-aa3a-10aaf21a1123
  modified: 2026-08-12T18:37:01.775Z
---

**Score a pre-registered prediction on its INPUTS, not just its output.** D-200 predicted
glutathione would take ~0.47 % of the quinone node and measured 0.32-0.53 %. The hit was
worthless: it was built from `sulfonation = 1.86 %`, which is the **low edge** of a band whose
**nominal is 15.21 %** (8.2× too low), multiplied by a concentration **10× too high** (a
loose bound from the citing chapter instead of the real level in the cited one). Two wrong
inputs of opposite sign, cancelling to within 1.5×. **Either error alone would have missed by
an order of magnitude**, and a scorecard that only compares predicted-vs-measured would have
recorded a clean hit and learned nothing.

**Why:** pre-registration earns its keep by making a miss *diagnosable*
([[feedback-pre-register-the-cheap-prediction]]). A hit is the case where nobody looks again —
so a hit is exactly where a wrong input survives, gets cited as validated, and propagates. The
same beat's *miss* (P2, 6.2× off) is what exposed the shared anchor error; had P1 been the only
prediction, the bad anchor would have shipped as confirmed.

**How to apply:** when a prediction lands, re-derive it from the measured quantities before
scoring it. If the inputs disagree with what you now know, the prediction **missed** — record it
as a miss with the cancellation named, however close the number came. Two specific traps this
beat hit, both worth checking by reflex:

- **The anchor.** A sensitivity written as a span ("1.86 % → 15.21 % → 59.69 %") has its nominal
  in the MIDDLE, not the front. Reading position 1 as the nominal is
  [[feedback-pin-the-band-not-the-nominal]] in its cheapest form.
- **The citing sentence vs the cited source.** A book can be ~10× loose against its own
  reference: UWC Ch. 24 argued from "GSH <30 mg/L, <0.1 mM; Chapter 5" where Ch. 5 says "no more
  than a few mg/L". Both supported the author's qualitative point, so the looseness was invisible
  until the number was used quantitatively. **Open the cited chapter, not the citing one** —
  the sibling of [[feedback-re-read-the-source-you-already-mined]] and
  [[feedback-a-notes-field-is-unchecked-storage]].

Related: [[feedback-a-named-pull-may-not-answer-the-question]],
[[feedback-conceded-caveats-are-not-coverage]].
