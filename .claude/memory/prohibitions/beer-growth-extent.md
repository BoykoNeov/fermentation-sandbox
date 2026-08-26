---
name: beer-growth-extent
description: "D-230/D-232 — beer's 1.55x growth-extent overshoot is NOT a nitrogen error; both D-222 candidates closed, BOTH ways of deciding the residue refused, and the residue is now THREE-way"
metadata:
  node_type: memory
  type: project
---

# D-230 — THE EXTENT GAP IS NOT IN THE NITROGEN BUDGET. Read this before proposing anything about beer's growth extent, wort nitrogen, or the count-to-gram conversion.

D-222 entered a **1.546x** extent overshoot (model **5.404x** against Tyrell's counted
**2.918-3.483x**) with two candidates and scored neither; D-228 restated it unchanged. **Both
candidates are now closed and the residue is a different kind of thing.**

## What must never be re-argued

* **`yan_mgl = 200` is CORROBORATED, not assumed-and-unchecked. Do not "fix" it.** This repo
  already held a sourced 10-12 °P malt wort's full free-amino-acid composition — Peyer 2017
  Table 16, transcribed at D-209 to derive `nitrogen_uptake_charge_beer`, whose **ratios** were
  used and whose **total** nobody had summed. Summed: **164.0 mg N/L** amino acids (reproducing
  the YAML's own printed cross-check) **+ 25-30 ammonium = 189-194 mg N/L assimilable**. The
  assumed 200 is **1.031-1.058x** that. The overshoot is 1.546x. **The nitrogen assumption
  cannot carry it.**
* **No FAN-to-YAN correction is owed — checked, not assumed.** Peyer's Table 16 contains **no
  proline** (18 acids, re-read from the full thesis), so the transcription IS the assimilable
  set. A FAN figure in an assimilable-N slot would have been a real 20-30 % scope error; it
  isn't one. **Never re-open this as "200 is a FAN number".**
* **The partition candidate ("not all assimilated N reaches suspended biomass") is REFUSED BY
  ARITHMETIC.** Hold the engine's 40 pg gram and ask what cell nitrogen would reproduce Tyrell's
  crop: **0.202-0.262 g N/g dry cell, i.e. 20-26 %**, against real yeast at 7-12 % and this
  engine's own `biomass_N_fraction` band topping out at 0.14 — **1.77x outside** it. No
  admissible partition works. Do not park it, do not re-propose it as unbuilt.
* **`yan_mgl` IS NOT AN ISOLABLE KNOB — measured, and it was NOT predicted.** Falsification arm A
  (set it to the 113 mg/L that lands the counts) turned **8** tests red where **4** were
  pre-registered. The four extra are every rate/timing test in the file: `mu_max` is fitted on a
  growth fraction **normalised on the peak, and the peak IS the nitrogen-limited ceiling**.
  Moving wort nitrogen re-opens D-222's `mu_max` refit. The full cascade is **five-way** (eight
  aroma constants, `Y_acetic_biomass_beer`, the pH cation seed, the second calibration frame's
  own 200, and the rate fit).
* **The residue is a TWO-WAY FRAME AMBIGUITY and neither branch may be picked silently.**
  Branch 1: Tyrell's cells are **70.9-91.9 pg** against the settled 40 (**1.77-2.30x**).
  Branch 2: **43.6-56.5 %** of the crop had already left suspension at the day-3 peak, needing
  no cell-mass change — live because Tyrell's own counts fall **22.4-45.1 %** in the single
  following day (an ENVELOPE-EDGE fall, not per-strain). **A later beat that closes the gap must
  say which branch it closed**; a repair that lands the fold without naming one has assumed the
  other away, and a guard enforces exactly that.
* **NO `Flags:` MARKER IS OWED ON D-219 AND NONE WAS WRITTEN.** Branch 1 is an **independent
  third** estimate of per-cell dry mass (Peyer's wort N × Roels elemental fraction ÷ Tyrell's
  counts — **no Coleman input**) and it lands near the ~100 pg D-219 retired as a back-computed
  residual. That is **not** a refutation: branch 2 explains the same gap with the gram untouched,
  and D-219's own "still open" (a count paired with a gravimetric dry weight in one ferment at
  one timepoint) is still unsupplied. What changed is that the engine's gram now has **two
  independent estimates ~2x apart** where it had one assumption. Prose only.
* **NEVER harmonise beer onto wine's nitrogen-dependent yield.** `_apply_nitrogen_dependent_yield`
  is gated OFF for beer deliberately. Coleman's `Y_X/N` regression (a0 3.50, a1 −3.61e-3) at
  200 mg/L gives `f_N ≈ 0.062` — **more** biomass per nitrogen than the elemental 0.114 — and
  would take the overshoot **5.40x → ~9x**. Wine earns it (Monod-on-N fitted directly on must);
  beer has no such fit. Same shape as wine's `mu_max`: "NOT an inconsistency to harmonise".
* **Settling stays rate-blocked, and that was GREPPED not guessed.** The five beer texts in the
  corpus describe flocculation only qualitatively (trigger, strain-dependence, cropping use) —
  no constant, no time course. What D-230 adds to that absence is **how much it would have to be
  worth** (branch 2's 44-56 %), which it did not have before.
* **One premise of D-213's O₂ decline HAS CHANGED — owner's call, not a beat's.** D-213 declined
  the oxygen→growth coupling partly because "none of the six directional predictions is reachable
  in the default set". Growth **extent** now is — three tests score it. The decline's other half
  (a second growth limitation re-opens a freshly calibrated `mu_max`, on a relation the corpus
  does not quantify) is **untouched and strengthened** by the isolability finding above.

## Shared-state note

Peyer's composition + the amino-acid chemistry table now live in `tests/conftest.py` as the
single sourced copy (joining `BEER_COUNTED_PITCH_CELLS_PER_ML`, D-228). `test_acidbase` keeps its
charge machinery and reads the shared table; its zbar re-derivation passing unchanged is the proof
the move preserved the table. **Do not re-transcribe either table anywhere else.**

# D-232 — BOTH WAYS OF DECIDING THE RESIDUE ARE REFUSED, AND THERE IS NOW A THIRD BRANCH.

The owner authorised **measuring which branch the evidence favours**, explicitly NOT moving
`cells_per_ml_to_pitch_gpl`. Nothing moved. What must never be re-invented:

* **The day-1-vs-day-3 settling profile is REFUSED — and NOT for being inert.** It **reverses its
  verdict** across `mu_max`'s own band: at 0.053 day 1 needs LESS settled than day 3 (flocculation
  possible), at 0.075 MORE (impossible). Span at day 1 is **0.348**; day 3 is inert at **0.0028**
  because it is the peak the fit normalises on. **Inert would mean "no information"; flipping means
  "information about the FIT".** Do not re-propose this test unless both band edges agree.
* **The pH panel as an independent clock is REFUSED, by a defect already on record.** The count panel
  measures biomass IN SUSPENSION, the pH panel measures nitrogen uptake = biomass MADE, and inverting
  a STATE map keeps `mu_max` out. But the whole signal sits at **day 1** (days 2-3 are at saturation),
  and day 1 carries the model's known **+0.162** alkaline miss — D-209 §8's unbuilt buffer-removal
  half. In the clock's own units that defect is worth **0.148** in nitrogen fraction against a signal
  of **+0.13 to +0.21**: it accounts for ALL of it, and is **70 %** of the width of Tyrell's own day-1
  count envelope. **Re-runnable ONLY once buffer removal is built** — a guard says so.
* **THERE IS A THIRD BRANCH and this beat ADDED it rather than closing one.** Wine yeast in must at
  330 mg N/L and ale yeast in wort at ~190 are not the same organism in the same medium. **Its
  nitrogen-level half is REFUSED and goes the WRONG WAY** — Coleman's own regression at 200 mg N/L
  deepens the disagreement from **2.03-2.63× to 3.25-4.21×**, which is the harmonisation fence above
  arriving from the count side on a quantity it was not derived from. The **organism/medium half is
  untouched and unsourced**, and is the residue.
* **The disagreement is best stated as a COUNT, and NOT as "two trials".** `biomass_N_fraction` and
  the 4e-11 gram both cancel (Coleman's `Y_X/N` IS a count × 4e-11), leaving **cells per gram of
  assimilated N**: engine-wine **2.515e11** vs Tyrell **0.955-1.236e11**. Central **2.03-2.63×**;
  headline **1.32-4.03×** across the regression's own credible regions; **1.45-1.73×** on the wrong
  (total-cell) convention. **The sign survives every band — nothing reaches 1.0.** But the Coleman
  side is **THIS REPO's fit**, so never call it two independent trials.
* **Branch 1 cannot be SIZED, only signed.** Across `biomass_N_fraction`'s own band the demand is
  **58-131 pg**; wort nitrogen moves it 3-6 %. The lowest edge anywhere is **57.8 pg = 1.44×** the
  engine's 40, so it never collapses onto the gram either.
* **The gram is NOT isolable from the rate fit.** Falsification arm A doubled it: **4 REDs predicted,
  10 got** — same structural reason `yan_mgl` isn't isolable. Restore SHA-256 verified, arm B green.
* **Sourcing is EMPTY and the negative HAS A SCOPE.** Patterns run over all five beer texts AND
  Peyer's 243-page thesis in full: `dry (cell )?weight`, `dry matter`, `g dry`, `per cell`,
  `cells per g`, `cells/g`, `pg`, `10^1[0-2] cells`, `pitching rate`, `million cells`, `cells/mL`,
  `g/hL`. Nothing. What turns up is a **third instance of the PRODUCT-dosing frame D-219 rejected**
  (*Concepts in Wine Chemistry*: 20 g/hL ↔ 1e6 cells/cc ⇒ **200 pg**). It is not a rival reading.

## Still open after D-232

1. **The extent gap itself** — attributed and bounded, NOT closed, and now **THREE-way**. Both
   published ways of deciding it are refused above.
2. **No settling / flocculation Process**, in either medium (rate-blocked).
3. **D-209 §8's buffer-removal half** — unbuilt, and now also the gate on the pH clock.
4. `mu_max`'s timing residue, `E_a_fusels`, wine's two calibration frames, `k_ester_volatil`'s
   sourcing — all where D-227/D-228 left them. Wine's nitrogen budget was never audited this way.
