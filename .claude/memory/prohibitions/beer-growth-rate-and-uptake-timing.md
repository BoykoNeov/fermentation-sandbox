---
name: beer-growth-rate-and-uptake-timing
description: "D-211 - beer's mu_max is re-derived from Tyrell Fig. 4's cell-count panel; the uptake-timing defect is CLOSED and the two parked terms are still parked"
metadata: 
  node_type: memory
  type: project
  originSessionId: dcd1ccd2-410e-4bfa-9d67-bd5e911c9dec
  modified: 2026-08-17T14:58:10.650Z
---

**Live prohibitions — beer's growth RATE and uptake TIMING (D-211).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read
it before proposing anything about beer's growth rate, nitrogen-uptake timing, or the two
terms D-210 parked. Every bullet is *what it forbids* + the record to read for *why*. **If a
prohibition looks unconvincing, go read D-211 — do not argue past it from this file.**

**The uptake-timing defect is CLOSED and its parameter is CALIBRATED. Do not re-propose it.**
- **`mu_max` (beer) is 0.034 /h, band 0.031-0.040, MEASURED** against Tyrell 2013 Fig. 4's
  **cell-count panel** — the third panel of the figure D-180 and D-207 had already cropped for
  extract and pH. The retired 0.098 was **2.88× too fast**. **Never restore it, and never
  re-cite Zamudio's 0.098 as this parameter's value** — Zamudio still supplies the FORM and the
  temperature scaling, not the magnitude: 0.098 is a **Droop** growth rate transferred by
  magnitude into **Monod-on-N**, the same unfinished translation D-15 fixed on the sugar side.
- **The EXTENT was already right; only the RATE was wrong.** Model multiplication 2.75× against
  a measured 2.92-3.48×, a **6.3 % shortfall pinned as its own recorded deviation** — do not
  fold it into the timing claim or "fix" it by moving the rate. The comparison that carries the
  result is **NORMALISED on each curve's own peak**, which is what cancels the cell-mass
  conversion — and that matters, because the 1.0 g/L pitch implies **~100 pg/cell against a
  textbook 40-60**. Day 4 is the **settling** limb (counts fall, the model cannot); never
  normalise on it. Both curves peak and fall — normalising on the FINAL X is wrong, `X` is
  drained by two inactivation Processes.
- **The pH course was NOT fitted and must not become the fit.** Counts admit 0.031-0.040;
  the pH course is 8/8 only at 0.038-0.048. **0.040 scores better and was REFUSED**: D-209 §8
  says its charge term is a LOWER bound (buffer removal unbuilt, pushes pH DOWN), so a day-1
  residual that is **too alkaline is what that term would close**. Day 1 went **0.315 too acidic
  → 0.070 too alkaline**; the 24 h drawdown fraction **>0.99 → 0.363**, inside the measured
  0.234-0.448. The band **spans a verdict change** (7/8 at the low edge, 8/8 at the high).
- **The `E_a_growth` arm is measured and REFUSED.** One temperature cannot separate `mu_max`
  from `E_a` — at 15 °C only their product is observable. **E_a = 204.6 kJ/mol reproduces the
  same counts from the OLD `mu_max`** and would leave 20 °C bitwise untouched, but it is
  **3.2× outside E_a's own band** and takes 10 °C attenuation 12.50 → **17.71 d**. The shipped
  arm also worsens 10 °C (→ 13.75 d); **neither number is Speers' midpoint statistic**, so
  neither scores against D-63's 2.1 d.
- **The band narrowed 7.00× → 1.29×, FROM THE FAST END ONLY, and `mu_max` is DRAWN** (it is the
  exemplar in `test_drawability_surface.py`). Attenuation spread across the band **0.79 → 0.25 d**.
  Edges are **CONSTRUCTED** — the mu keeping days 1-2 inside the four-strain spread **at the
  nominal E_a**. **Never fold E_a's band into this one**: both are drawn, it double-counts. The
  retired band was a 17-26 °C temperature range, which double-counted the Arrhenius axis.
- **Wine's `mu_max` 0.095 is now 2.8× higher and that is NOT to be harmonised** — Coleman fitted
  Monod-on-N **directly** (`test_coleman_reconstruction.py` rebuilds those curves with this same
  form), so wine never went through a Droop transfer. **The old note's "matches the wine value"
  corroboration is RETIRED.** No flag is owed on the wine parameter.
- **The Chemistry of Beer's "Class I amino acids within the first 20 hours" is NOT a refutation** —
  it is the diagnosis. This model has **one lumped N pool** and was applying the **fastest
  class's** timescale to all of it. Do not cite it to restore the old rate.

**SUPERSEDED AT D-214: both parked terms are now CLOSED, not parked** — antiport has zero beer-text
sourcing, trub is a pre-pitch event already inside the calibration and a charge violation after the
anchor. §9's own two numbers are corrected to **0.0274** and **0.0086** (an adaptive-grid read; see
`.claude/memory/prohibitions/trub-settling-and-the-peptide-pair.md`). The paragraph below is kept as
written because the *brief* it states — beer wants acidification EARLY and none late — still stands
and is still unanswered; only its two candidate answers are spent.

**D-210's two parked terms are STILL PARKED — measured, not inferred.** K⁺/H⁺ antiport and trub
protein settling were parked on the **high `z̄` edge**, not the nominal. Under the new rate the
day-7 headroom above the floor goes **0.0028 → 0.0082 pH**: tripled, still closed to a constant
same-sign term. **Never answer this from the nominal.** What DID change is the shape: day 1 now
sits **0.034 pH ABOVE** its ceiling at that same edge, where it sat 0.48 below. So beer's pH
wants **acidification early and none late** — a term with a **time profile**, not a constant.
That is the brief for the next beer-pH beat.

**Three prose numbers were conditioned on the old rate and are now fixed — do not "re-correct".**
**`q_sugar_max`'s high edge 1.5 is STILL RIGHT**: its `15.3 * 0.098` uses **Zamudio's** growth
rate for **their** peak, not this file's. `15.3 * 0.034 = 0.52` lands on the nominal 0.5 and is a
**trap**, not a discovery — and `q_sugar_max` is **banded and drawn**, so it is a live sampling
surface. Its "1.5 attenuates a 1.048 wort in ~2 d" is **re-measured at 2.58 d** (D-15 unchanged).
`test_beer_temperature_response.py`'s 10 °C midpoint is **7.46 d / 3.55×**, was 6.2 d / 2.9×.

**D-183 is FLAGGED, not corrected.** Its growth-linked acetic producer still beats the retired
flux-linked one (RMSE 40.7 vs 65.3), but the margin moved 0.528 → **0.624** because it had been
scored against the wrong growth rate. Acetic rises faster than growth **and** than flux: the
model books **0.360** of its rise by day 1 against a measured **0.773**, a **2.15×** shortfall
now pinned. Its *"86 % of acetic's entire rise"* is a **units slip corrected in prose only** —
86.00 is the day-1 rise in **mg/L**; the share is **77.3 %**. No shipped value moved.
