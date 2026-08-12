---
name: feedback-a-magnitude-defence-can-be-scoped-to-the-smallest-member
description: "A conceded lump defended with a number: check the number covers the WHOLE lump, not one member — D-75's was true and 1000x too small"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 720cec99-de09-483a-bafb-7296d361e840
  modified: 2026-08-12T12:12:27.646Z
---

When a record concedes a known approximation and defends it with a magnitude ("it is only ~2 %,
a minor perturbation"), **verify the number spans the whole thing being defended, not one member
of it.** A true number about the smallest member reads exactly like a bound on the lump.

D-75 conceded that the direct oxidative set charges oxygen twice for one shared chemical step, and
defended the size with `k_strecker` — "~2 % of the always-on total, well under 1 % of the O₂".
Measured at D-197, both figures are **right**: Strecker takes 0.012–0.033 % of consumed oxygen, and
the 2 % is an exact ratio of two constants. But the lump is **five** draws, and the one the
concession never names takes **~1000× more**; with sulfite present the five take **62–76 %** of the
oxygen. Nothing was wrong — and the paragraph still implied a bound two orders of magnitude too
small, for 122 records.

**Why:** a magnitude defence is a *scope* claim wearing a number. The number gets audited (it was
correct); the scope never does, because the sentence names one term and the reader supplies the
"…and the rest is like it". The same shape as [[feedback-a-scope-note-can-carry-a-mechanism-claim]]
and [[feedback-a-notes-field-is-unchecked-storage]], but harder to catch: there is no wrong value to
find.

**How to apply:** before accepting a conceded approximation as small, enumerate every term the
concession covers — from the registry, not the prose — and measure **each** share. Ask "is the term
this number describes the biggest one?" If it is not, the defence is unmade even though the number
survives. Then check what the deviation is worth against data before treating removal as an
improvement: D-197's naive de-duplication moved the model from 15 % below the measured band to 2.3×
below it, so the lump was load-bearing and the honest output was a measurement plus
[[feedback-closer-to-reality-decides]] in reverse — record the size, keep the behaviour.
