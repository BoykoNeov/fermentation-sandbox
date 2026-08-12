---
name: feedback-a-band-is-per-parameter-a-claim-is-joint
description: "An uncertainty band is stated per parameter but inhabited jointly, so a prose claim about TWO parameters' relationship is covered by neither band"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: f829484e-0ed2-4cb0-97ec-aaac6b46ce0c
  modified: 2026-08-12T06:26:30.160Z
---

**A claim about the *relationship* between two parameters is not covered by either parameter's
uncertainty band.** At D-190 two SO2 binding constants each carried a diligent note naming the
competing secondary value and a band widened to span it. Read one at a time, nothing was missing.
But the two alternatives were **the same alternative** — one source's pair — and adopting it flips
which carbonyl binds more strongly. The parameter file asserted that ordering as a fact
("slightly STRONGER than pyruvate ... the second-most-avid binder") while its own overlapping
bands, sampled independently, reversed it in **1.996 % of 200 000 draws**.

**Why:** per-parameter review is the natural unit and it is structurally blind here. Each band was
*correct*; the sentence spanning them was wrong. No amount of checking either entry on its own
surfaces it, because the defect lives in the region of the joint space neither note describes.

**How to apply:** when a comment or docstring states a **relationship** — stronger than, faster
than, always exceeds, dominates — treat it as a claim about a **region of the joint parameter
space**, not as a gloss on two values. Then sample that space through the shipped sampler and
**count how often the sentence is false** (draws, not band edges — an edge argument proves the
corner exists, not that it is visited). If the count is not zero, the sentence is wrong even when
every value and every band is right; fix the sentence, and do not move a value to rescue it.

Corollary from the same beat: **check that your mutation is inside the sampled space before
reading its red.** Swapping the two nominals produced 461 failed / 401 errors that looked like
physics and was a parameter *load* failure — one value sat 6.7 % below its own band low. An
unreachable configuration's magnitude is not a consequence of anything.

Related: [[feedback-pin-the-band-not-the-nominal]],
[[feedback-nominal-on-a-band-edge-is-not-inertness]],
[[feedback-pair-the-red-with-an-ordering-preserving-baseline]],
[[feedback-a-notes-field-is-unchecked-storage]],
[[carbonyl-release-and-binding]].
