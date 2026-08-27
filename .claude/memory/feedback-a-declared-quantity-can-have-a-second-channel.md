---
name: feedback-a-declared-quantity-can-have-a-second-channel
description: "A constant derived at compile from a declared input stays fitted at that input while a SECOND input adds to the same physical quantity — so the run leaves the point its own constant was fitted at, and the two errors compound"
metadata:
  node_type: memory
  type: feedback
---

When a model derives a constant by evaluating a fit at one scenario input, that derivation silently
assumes the input **is** the physical quantity. Add a second knob that contributes to the same
quantity and the assumption breaks without any code changing: the run now sits somewhere the
constant was not fitted for, and nothing in the ledger complains because conservation is about
totals, not about evaluation points.

**Why:** wine's `yan_mgl` seeds the ammonium slot **and** is where Coleman's yield regression is
evaluated (D-14). `amino_acids_gpl` seeds eight pools that are also on the nitrogen ledger (D-100).
They **add**; they do not partition. At the suite's commonest dose a scenario declaring 250 mg N/L
carries **362.7**, and the yield stays fitted for a 250 must — a poorer must than the run is — so
more nitrogen and a higher yield-per-nitrogen compound in the same direction. D-32's own text says
*"amino acids are part of YAN"*, which is the premise the seam does not implement; D-14 predates it
and its "all assimilated nitrogen enters biomass" was true when written. Nitrogen conservation
closes to 3e-14 throughout. **This is a declaration defect that no conservation test can see.**

**How to apply:** for every constant computed at a scenario boundary, write down the physical
quantity its fit is a function of, then grep for **every** input that moves that quantity — not just
the one the derivation reads. Where you find a second, measure the gap between declared and actual
at a realistic setting before proposing anything. Beware the obvious repair: summing the channels
into the fit took wine's commonest dose outside the regression's own fitted range and reached
absurd values at larger ones, so it trades a wrong evaluation point for an extrapolated one. Two
bad conventions with prices is the honest deliverable; pick neither silently
[[feedback-a-scope-note-can-carry-a-mechanism-claim]].
