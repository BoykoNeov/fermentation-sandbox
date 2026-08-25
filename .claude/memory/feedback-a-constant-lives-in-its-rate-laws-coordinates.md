---
name: feedback-a-constant-lives-in-its-rate-laws-coordinates
description: "A parameter's printed value and its prose are in the coordinates of the rate law that reads it; change the law and the same number means something else — both of beer's aroma activation energies described the rate constant while their notes were read as describing the output"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 9d970785-c90f-4a56-ae65-f4ba3014180b
  modified: 2026-08-25T18:50:35.141Z
---

Before changing a rate law, work out **what each constant it reads will then MEAN**, and check the
existing prose against the OBSERVABLE rather than against the constant.

D-226 moved beer's aroma pools from flux coupling to growth-extent coupling. Under the flux form,
run-integrated synthesis scaled as `arrh(E_a)/arrh(E_a_uptake)`, because the bare-flux integral to
dryness is fixed by total sugar. So:

* `E_a_esters` printed **200,000 J/mol** while the apparent activation energy of the thing you can
  measure was **144,900** (= 200,000 − 55,100).
* `E_a_fusels` printed **70,000** with a note claiming **"Q10 ~ 2.6"**, while the higher alcohols
  actually rose **1.2183×** over 15-25 °C — **Q10 1.22**. The note was true *of the rate constant*
  and had been read for 200-odd records as though it described the beer.

Under the extent form the growth Arrhenius cancels against the nitrogen limit, so each parameter
becomes the apparent activation energy **directly**. The proof that this was a coordinate change and
not a discovery: reverting `E_a_fusels` to 70,000 under the NEW coupling makes the pools span
**2.6644×** — exactly the Q10 the note always claimed. The note finally became true when the law
changed underneath it.

**Why:** a note that describes a constant reads identically to one that describes an output, and
nothing in a parameter file distinguishes them. The discrepancy is invisible while the law is
stable, and it surfaces as a "surprise" the moment the law moves — at which point it is easy to
mistake a bookkeeping identity for a physical finding.

**How to apply:**
- Derive the map from parameter to observable **algebraically first, then confirm it on the engine**.
  For D-226 the map was exact: flux ⇒ `E_a − E_a_uptake`, extent ⇒ `E_a`. `arrh(200,000)` = 16.44×
  and `arrh(144,900)` = 7.59× reproduced the two measured columns to four figures.
- When re-anchoring to preserve behaviour, **solve on the engine, not on the algebra**. The
  synthesis-only algebra said 144,900; the answer was 152,000, and the 4.9 % gap was a downstream
  sink re-timing. [[feedback-fit-the-observable-not-the-consequence]]
- **Re-read every note on every constant the law touches**, and ask of each sentence: is this about
  `k`, or about what comes out of the vessel? [[feedback-a-notes-field-is-unchecked-storage]]
- An ordering constraint between constants is usually the *algebraic form, in one coupling*, of a
  physical claim. D-19's "every `E_a` must exceed `E_a_uptake`" is the flux form of "output rises
  with temperature"; the extent form is "`E_a` > 0". Restate it per coupling — and assert **which
  coupling the medium is wired to**, or the test quietly checks the wrong inequality forever.
- Guard the observable. Every level guard here scored the calibration frame, which is `T_ref`, where
  every Arrhenius factor is exactly 1.0 — so **1,857 tests were blind to a 2× change in the
  temperature response**. [[feedback-a-threshold-cannot-separate-same-sign-regimes]]
