---
name: prohibition-beer-growth-cutoff
description: "D-228 — the growth cutoff blamed for beer's aroma taper was the calibration frame's inoculum; what is corrected, what is guarded, and what stays open"
metadata: 
  node_type: memory
  type: project
  originSessionId: 870e1cf7-7f86-43eb-b8df-19e2af8180be
  modified: 2026-08-26T07:38:04.844Z
---

# D-228 — THE CUTOFF WAS AN INOCULUM. Read this before proposing anything about beer's growth cutoff or the aroma taper.

D-226 §8, D-227 §10 and this repo's own ledger all named "growth stops dead at nitrogen
exhaustion (day 0.92) with 81 % of the sugar still to ferment" as the inherited limitation
behind beer's aroma taper, and as the next thing to build. **Half of it is a property of the
frame it was measured in, and that half is fixed.**

## What must never be re-argued

* **The aroma calibration frame pitched a flat 1.0 g/L and now pitches a COUNT.**
  `_beer_calibration_scenario` (21 d / 20 C / YAN 200, where all EIGHT beer aroma constants are
  defined) carried the same back-computed residual D-219 retired and D-222 removed from
  `TYRELL_SCENARIO` — a gram of dry-yeast PRODUCT per litre, ~100 pg/cell against a settled 40.
  At 40 pg it is **2.5e7 cells/mL, 2.5× a counted ale pitch**. Now Foster 2022's counted
  **1.2e7 cells/mL** through `cells_per_ml_to_pitch_gpl`. **Never put a flat g/L back.**
* **Nitrogen-limited growth reaches a ceiling fixed in ABSOLUTE terms**, so a heavier inoculum
  reaches it sooner. That is the whole mechanism of the "cutoff": 90 % of the gain takes
  **0.785 d at 1.0 g/L and 1.199 d counted (1.53×)**; D-226's "day 0.92" is the retired frame's
  99 % point. The five higher alcohols sat at **0.9999 of finished by day 1** and now sit at
  **0.672**; isoamyl acetate 1.312 → 0.613.
* **The correction is FREE and that is measured on the SHIPPED rate law.** All eight levels,
  enumerated from `ESTER_SPECS`/`FUSEL_SPECS`, move **6e-7 relative (6e-5 %)**. D-226's
  "pitch is inert" was measured when the sink read the flux SHAPE; **D-227 moved it onto evolved
  CO₂**, so it was re-measured rather than inherited. **Every level guard is BLIND to this**, so
  the beat's only visible guard is the growth WINDOW.
* **On the trial that HAS counts the growth curve was ALREADY right.** Tyrell Fig. 4, scored with
  the repo's midpoint-denominator helper: day 1 **0.299** in 0.235-0.448, day 2 **0.910** in
  0.781-1.012, day 3 **0.993** in 0.872-1.128. **And the peak day agrees when the model is
  sampled DAILY as the trial was** — continuous `argmax` says 2.32 d, daily sampling says **day 3**,
  measured day 3. Never quote the 2.32-vs-3 gap as a discrepancy.
* **The taper stays REFUSED and the half that died is the DIAGNOSIS, not the refusal.**
  Luedeking-Piret is still two free parameters against one observable with **no beer ester
  time-course on disk** (re-grepped at D-228). What is dead is "the cutoff is why".
* **The 1.55× extent overshoot stays D-222's open item.** Entered, not closed. **The counting
  method does not settle it**: MEBAK III 10.4.3 counts cell CONCENTRATION with no viability
  discrimination, so the comparable quantity is `X + X_dead` — **5.404 vs 5.378 fold, 1.55×
  either way**. D-222's two candidates (the assumed 200 mg/L YAN, Tyrell prints no FAN; and the
  booking of one lumped N pool wholly into suspended biomass) are unchanged.
* **The day-4 decline is a SEPARATE absence.** Tyrell's counts fall 30-45 % from day 3 to day 4;
  ethanol death is far too slow at ~4 % ABV and `X + X_dead` is monotone. That is yeast leaving
  suspension and **there is no settling Process in either medium**. Never fold it into the taper.
* **Scope fences, deliberate:** §2.2's benchmark pitch (0.6 g/L, "homebrew-like") is UNTOUCHED and
  has its own reasoning; wine's pitches are out of scope — this was not a repo-wide pitch audit;
  `mu_max`'s timing residue (D-227 §7) is unentered.

## Still open after D-228

1. **The 1.55× extent overshoot** (D-222's, above) — needs Tyrell's FAN or a nitrogen-partition
   change, not a rate.
2. **No settling / flocculation Process**, in either medium.
3. **`mu_max`'s timing residue**, `E_a_fusels`'s magnitude, wine's two calibration frames,
   `k_ester_volatil`'s sourcing — all where D-227 left them.
