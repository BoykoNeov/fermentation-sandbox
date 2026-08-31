---
name: price-a-transient-against-a-ceiling-of-the-same-kind
description: "Before a band is used as a ceiling, check it describes the same QUANTITY and the same MOMENT - an end-state composition caps nothing transient"
metadata:
  node_type: memory
  type: feedback
---

D-251 declined a calibration because it loaded cells to 0.139 g N/g dry weight against
`biomass_N_fraction`'s [0.08, 0.14]. That band is a canonical N-*replete* **end-state elemental
composition**. The quantity being priced was a **transient store**, at the model's own maximum, at
a moment the source never sampled. Crepin's own whole-run value sits at or below **0.0786** - under
the band's LOW edge - so the band was never a ceiling on anything in that ferment, and at the one
moment she did sample, the calibration was carrying LESS nitrogen per gram than her data imply.

**Why:** a declared band looks like a ceiling because it has a high edge. It is only a ceiling on
the quantity it was measured as, at the state it was measured in. Pricing a peak against a
steady-state band is a category error that reads as diligence.

**How to apply:** name the quantity, the basis (per gram of NEW cell vs per gram of ALL biomass -
they differ by exactly the inoculum) and the moment, on BOTH sides, before comparing. Then check the
source's own value against the band: if the measurement falls outside the band entirely, the band is
the wrong object and no comparison against it will mean anything.
Related: [[feedback-a-reachable-target-can-be-unaffordable]], [[feedback-compare-at-the-moment-the-source-sampled]].
