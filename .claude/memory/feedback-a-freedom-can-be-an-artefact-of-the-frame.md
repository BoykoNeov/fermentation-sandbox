---
name: feedback-a-freedom-can-be-an-artefact-of-the-frame
description: "An 'inert / free / decoupled' property is a claim about a frame, not the model; correct the frame and it can vanish entirely - and it takes every claim measured in that frame with it"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 47cfa586-43a7-4fdb-9cf6-b1f266627a26
  modified: 2026-08-18T15:15:22.627Z
---

When a parameter is recorded as **inert**, **free**, or **decoupled from an anchor**, that is a
claim about the **frame the anchor was asserted in**, not a property of the model. Correct the
frame and the freedom can vanish outright — not shrink, vanish — and there may be no weaker
version of the claim to fall back on.

At D-221 the beer acceptance criterion was re-temperatured 20 → 15 °C, because the literature
supports its duration at 15 and contradicts it at 20. `E_a_uptake` had been named "the ONE knob
that moves one anchor without moving the other" and that was **measured at exactly 0.0000 d** —
correct, and correct only because the two anchors sat at different temperatures. The criterion's
new temperature is the other anchor's temperature **exactly**, so one Arrhenius factor now enters
both: the printed band moves the criterion **1.75 d against a 2.0 d window**. Most inert → one of
its strongest levers, in one frame change. I had assumed a magnitude argument would survive
("restate the exactness as a small number"); there was none.

**Why:** exactness is seductive — a measured 0.0000 reads as structural, so two records built a
refusal architecture on it. But it was arithmetic about `T_ref`, and the frame it was arithmetic
*in* was the thing under review. A freedom that exists "by construction" is only as durable as
the construction.

**How to apply:**
- Before citing something as inert/free/decoupled, ask **what frame makes it so**, and whether
  that frame is itself a live claim. Write the frame into the note, not just the number.
- **Measure the blast radius before designing the edit.** Changing a shared anchor moves every
  claim ever measured against it. I predicted one test would move; **six** did, and one carried
  the decoupling argument. Patch the value temporarily, run the affected files, read the list —
  it costs one run and it changes what you build.
- **Preserve the retired frame's measurement, don't delete it.** A claim measured in a frame that
  has moved is still true *about that frame*; assert both, because the DIFFERENCE between them is
  usually the finding.
- **A guard staying green is not evidence its prose survived.** One test here passed throughout
  while its docstring asserted the dead claim — it never read the anchor that moved. Grep the
  prose, not just the asserts. [[feedback-grep-finds-claims-not-guards]]
- When a frame correction makes two sources suddenly **agree**, sweep before claiming it: the
  agreement here held over 69 % of a band and failed on the rest.
  [[feedback-a-hit-can-be-two-errors-cancelling]]

See [[prohibition-beer-criterion-temperature]], and
[[feedback-agreement-can-be-a-frame-difference]] for the same lesson from the other direction.
