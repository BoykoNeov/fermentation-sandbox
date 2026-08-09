---
name: feedback-nominal-on-a-band-edge-is-not-inertness
description: "If a parameter's nominal EQUALS a band endpoint, the run at that endpoint is the nominal run — bitwise identical by construction, not by inertness — so an edge screen manufactures a straddle"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 10a4c2c5-2d8d-4887-b416-58ee04aef9d0
  modified: 2026-08-09T13:34:31.167Z
---

When classifying a parameter by pinning it at its band `low` and `high` and comparing
against the nominal run, any parameter whose `value` **equals** one of its own endpoints
returns bitwise identical at that edge **by construction**. If the other edge moves, a
two-point screen reports "inert at one edge, active at the other" — a STRADDLE that is
purely an artefact of the band's geometry.

D-166 found 3 such live bands in wine (`f_non_ehrlich_phenylalanine`, `copper_h2s_binding`,
`copper_mercaptan_binding` — all on the *high* edge, all unread, so harmless) and **2 in
beer that were not harmless**: `Y_glycerol_sugar` and `Y_byproduct_sugar` sit on the *low*
edge at value 0.0, are read by an active Process, and are genuinely ACTIVE across the band
interior. Uncorrected, the census would have reported exactly two straddles in beer that do
not exist. These are the same two bands D-165 §2's ratio filter had to drop for a zero
nominal — the same degeneracy surfacing from the other side.

**Why:** the screen's premise is "same output ⇒ the parameter did nothing". That premise
fails when the perturbed map is not actually perturbed. The failure is silent and it points
the wrong way — it *invents* a finding rather than hiding one, which is worse, because an
invented straddle is exactly the kind of result that gets written up.

**How to apply:** before any band-edge sweep, compute the set
`{n : value(n) in (low(n), high(n))}` and print it. Classify those over the band
**interior** instead (e.g. 2/25/50/75/98 % of the span), where no sample point coincides
with the nominal. Pin the set in a test so a new one cannot appear unnoticed — and check it
**per medium**, because the dangerous ones were in beer while wine's were all benign.

Related: [[feedback-pin-the-band-not-the-nominal]],
[[feedback-a-null-result-needs-a-positive-control]], [[feedback-count-and-print-your-skips]].
