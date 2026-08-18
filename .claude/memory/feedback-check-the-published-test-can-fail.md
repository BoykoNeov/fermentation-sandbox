---
name: feedback-check-the-published-test-can-fail
description: "Before reporting that the model PASSES a published test, sweep the band and check the test can fail; mine was cleared by every in-band value and rested on a sampling-grid point read as exact"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: bd744d4f-5794-466f-81ff-deb7e49f9bff
  modified: 2026-08-18T09:13:32.284Z
---

Before writing "the model passes <published test>", **sweep the parameter band and confirm some
in-band value FAILS it.** A bound no admissible configuration can cross is not evidence about the
configuration that ships.

D-218 found the archive's first same-yeast two-temperature pair — the exact shape D-217 had
searched 26 files for and not found. Confound-free by construction: same strains, same wort, same
pitch, so the duration **ratio** is immune to every conversion and species confound in the
problem. The model came in at 2.21 against a published ≤3.33 and I reported a pass, and told the
owner the temperature response was vindicated.

It was cleared by the **entire** printed band — 1.59 at the low edge, 2.44 at the high — and the
first value that fires is 1.43× **out of band**. The archive already had this rule twice, stated
about *mutation arms* (D-216 §10, D-217 §6: an out-of-band arm tests nothing). I did not recognise
it in the mirror-image form, where the non-discriminating thing is the **literature**.

Two things it cost, past the wrong headline. It nearly licensed re-opening a refusal D-217 had
closed on a counted census. And a second bug hid behind the first: **the bound itself was looser
than I read it.** The source's "3 days" cites a figure panel, and the panels are *timepoints* —
12/48/72/120 h. So "3 days" is the 72 h **sample**, a ceiling, not a duration; an upper bound on a
ratio needs a *lower* bound on the denominator, which the paper never gives. Read at the open end
the bound is 5.0. I only went looking for that after the sweep showed the pass was empty.

**How to apply.** When a source hands you a number to score against: (1) sweep your parameter's
printed band and record the value at which the test first fires — if that value is out of band,
the result is "this literature does not constrain us", never "we pass"; (2) ship the out-of-band
arm as the guard's **positive control**, or "everything clears" is indistinguishable from a broken
predicate ([[feedback-a-null-result-needs-a-positive-control]]); (3) before trusting the bound's
tightness, find where in the paper the number was *read* — a figure panel is a sampling grid and
its numbers are ceilings ([[feedback-transcribe-tables-not-prose]]); and (4) check the direction:
a bound on a ratio needs bounds on both halves, and a ceiling on the denominator is worthless.

Sibling failure the same beat: the ratio agreed with the bare Arrhenius factor to 0.7 % at the
nominal, which looks structural. Swept, the residual crosses 1.0 *inside* the band — a coincidence,
not a law. Same discipline, same fix: **sweep before you attribute.**
