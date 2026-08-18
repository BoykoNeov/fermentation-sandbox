---
name: beer-second-measured-course
description: "D-220 - the archive's second measured beer fermentation course, recovered from Foster's vector supplementary figure: sec 2.2's window is mis-temperatured rather than refuted, the engine is ~1.5x slow below 30 C, and the 30 C agreement is a crossing"
metadata: 
  node_type: memory
  type: project
  originSessionId: d9b6cdf3-1baf-48f5-bebb-63068e2c2118
  modified: 2026-08-18T11:12:37.113Z
---

**SUPERSEDED IN PART BY D-223 — read this block first.** Beer's two speed anchors were
ADJUDICATED and the rate MOVED. `q_sugar_max` 0.5 -> **0.72 g/g/h**, band 0.3-1.5 -> **0.634-0.818**
(1.29x, and it is DRAWN, so every beer ensemble moved).
1. **Foster's measured 15 C course + sec 2.2's window WIN; Tyrell's day-2 extract point LOSES.**
   They are ONE anchor, not two: every rate reproducing Foster's course lands inside the
   criterion's admissible interval and vice versa. The Tyrell-matching 2.3226 reads 2.38 d
   against the criterion (outside, FAST) and is **2.839x** the re-derived printed high edge.
2. **sec 2.2's criterion is 6.04 d, INSIDE 5-7 for the first time since D-221.** Its strict xfail
   is gone. The pass is PARTLY SELF-REFERENTIAL (D-221 set the criterion's temperature from the
   same paper the rate is fitted to) and must never be cited as third-party corroboration.
3. **"The engine is uniformly ~1.4x too slow below 30 C" is DEAD.** One rate fitted at 15 C alone
   lands **0.945 / 1.030 / 0.973x** at 12 / 15 / 22 C — and 12 and 22 were NOT fitted, so that is
   an out-of-sample check of the temperature response. D-220 sec 4's reading was right and
   D-217's refusal to re-source `E_a_uptake` still stands.
4. **The 30 C column INVERTED**, 0.659x: the model is now 1.52x too FAST there. Real ale yeast
   saturates above ~22 C and this model does not. That is the price the re-anchoring paid and it
   is NOT repaired — do not "fix" a single column by inventing a saturation term.
5. **Prices, measured:** day-7 pH margin 0.036 -> **0.0033**; affordable trub peptide loss
   12.6 % -> **1.2 %**; beer's whole VDK ladder and both ester pools fall ~26-31 % (they are
   biomass-hour-linked, not flux-linked). Tyrell's own day-7 attenuation IMPROVES 0.782 ->
   **0.959** and his day-2 shortfall 4.21x -> **3.16x**; D-215's extract xfail stays xfail.
6. **A defect the faster engine EXPOSED (not created):** `EthylAcetateEsterification`'s formation
   half was funding itself out of beer's `Byp`, which is exactly 0 by construction (D-16). It now
   scales that half by `clip(Byp / acetic_acid_typical, 0, 1)` — one-sided, wine bitwise
   unchanged (4.77x clear of the clip). **Beer can re-form ester only from acid it released
   itself.** D-176's "beer hydrolyses so the pool is safe" was a MARGIN and it is spent.
**OPEN, and first in line:** beer's two ester calibrations now disagree about whether the faster
engine is better — ethyl acetate moves AWAY from its published mean (21.3 -> 15.9 vs 23.7) while
isoamyl acetate moves TOWARD its target (3.83 -> 2.80 vs ~2.2). Both `k`s were fitted at the
retired rate; one is wrong and D-223 could not say which.

**Live prohibitions — the SECOND measured beer course (D-220).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path and
does **not** index it in `MEMORY.md`. Read it before proposing anything about beer's ferment
speed, §2.2's attenuation window, Foster 2022, or the model's peak timing. **If a prohibition
looks unconvincing, go read D-220 — do not argue past it from this file.**

**MEASURED, three corrections, NOTHING BUILT.** No parameter moved, no knob touched.

- **The course is ON DISK and was inside a paper the archive had already mined TWICE.** Foster
  2022's **Supplementary Figure S1** — 10 strains x 8 temperatures (12/15/22/30/35/37/40/42 C)
  x 9 timepoints, 0-120 h — arrived in the Europe PMC bundle (`PMC8966892/supplementaryFiles`)
  nobody had requested. D-218/D-219 read only main-text Fig. 2, so every duration the archive
  quoted from Foster was a **CEILING off a sample**. **D-216's "no second measured beer extract
  course on disk" is now FALSE.** Full 720-point recovery in
  `M:\claud_projects\temp\ferment\d220-second-beer-course\`.
- **The values are a TRANSCRIPTION, not a plot reading — never re-derive them by eye.** The
  figure is VECTOR; every point is a path vertex mapped through the axis ticks. The layout was
  established from geometry: the polyline is UNDODGED (its vertices are the means), the error
  bars are dodged. Two checks it was not told about: the t=0 gravity agrees across all ten
  panels per temperature to four decimals, and the **40-42 C stall** pins the colour ordering
  (12-30 C alone are monotone and cannot).
- **§2.2's 5-7 d window is NOT refuted — it is MIS-TEMPERATURED, and only the PAIR is wrong.**
  At **15 C all three Beer 1 controls land INSIDE 5-7 d** (5.04/5.39/5.91); at 12 C two are
  slower than the slow edge. The brief's **DURATION is corroborated**. What cannot be true is
  that duration **at §2.2's 20 C**: the same strains take 2.91-3.76 d at 22 C. **Corrects
  D-218 §4 and D-219 §5(c)'s "nothing found supports 5-7 days"** — read off timepoint panels
  dominated by the warm end. **Always argue this as a BRACKET, never the interpolation** (3.40-
  4.27 d at 20 C is recorded and asserted NOWHERE — its cold input is itself an extrapolation).
- **The cold columns EXTRAPOLATE past 120 h and the bias has a KNOWN SIGN — never treat it as
  softness.** The tail DECELERATES (Cali 15 C: 3.43 then 1.67 SG-points/h), so true crossings
  are **LATER**, and Foster's 1.045 vs §2.2's 1.048 pushes the same way. Every claim is written
  to be strengthened by it.
- **The engine is ~1.5x SLOW and it is a LEVEL error, not a temperature-response one:**
  **1.41x / 1.54x / 1.45x at 12 / 15 / 22 C** (model 10.742 / 8.379 / 4.838 d vs measured means
  7.60 / 5.45 / 3.34), at Foster's own wort and counted pitch through §12's helper. Near-constant
  across the range ⇒ **D-217's refusal to re-source `E_a_uptake` is VINDICATED on new evidence**;
  the lever is a rate SCALE (D-216's `q_sugar_max`), whose refusal rests on §2.2 — **the premise
  this beat moves, not the arithmetic.** **§2.2's window is what HIDES the slowness.**
- **NEVER cite the 30 C agreement as the model getting beer's speed right — it is a CROSSING.**
  Measured apparent activation energy **COLLAPSES 49.5 → 17.9 kJ/mol** across 22→30 C (real yeast
  saturates); the model holds **55.5 → 53.7**. Straight line vs flattening curve. That endpoint
  E_a is a **LUMPED coefficient of the whole course, NEVER a reading of `E_a_uptake`** (D-183's
  error, D-217's refusal) — and at 45.5-55.5 it is **in band**, so it discriminates nothing.
  **D-218's "never cite Foster's temperatures as a temperature test" STANDS, for a NEW reason.**
- **D-218's "every measured course peaks day 1-2, and that is STRUCTURAL" is TEMPERATURE-
  CONDITIONAL.** Plateau of steepest fall ends 24 h at 30 C, 24-36 h at 22 C, **48-72 h at 15 C**,
  and **no peak inside 120 h at 12 C**. Report as a **PLATEAU, never an argmax** (15 C Cali runs
  3.53/3.45/3.43 across three intervals). **This does NOT rescue the model** — at Foster's OWN
  conditions the engine peaks later still at every temperature. The first draft claimed it did,
  on a cross-frame comparison against Tyrell's wort/pitch/temperature; the defect survives, its
  stated reason does not.
- **Foster's OG statements FORK** (12.5 °P ≈ 1.0505 vs the printed 1.045). **Neither is used** —
  the figure's MEASURED per-run t=0 gravity is what attenuation is scored against. Named, not split.
- **Cramer 2002 is CLOSED on five enumerated hosts** (Wiley abstract-only, Semantic Scholar
  CLOSED, Unpaywall not-OA, ResearchGate request-only, no institutional/thesis copy). **D-219's
  tier stays `plausible`.** **de Andrés-Toro 1998 is still blocked, now CHARACTERISED**: the
  browser route — the one D-218 never tried — returns a **bot-detection challenge**, a REFUSED
  route, not a failed one; and PMC9689312 does not reprint the data (its own batches, 19-28 C).
- **Guards: `tests/test_organic_acids.py` §14, five tests, ~4 s.** Model durations are
  **INTEGRATED LIVE, not pinned** (a pinned duration fires on drift and cannot see a model
  change — D-216 conceded exactly that about its own pitch test). Falsified on four in-band arms
  with restore verified between each: `q_sugar_max`→1.397 reddens 3, `E_a_uptake`→63000 reddens
  2, `mu_max`→0.040 reddens 1, `E_a_esters`→210000 is the **designed GREEN**. The two
  literature-only tests stay green under every model arm **by design** — they are claims about
  Foster's data.
- **Pre-registration scored 3 MISSES of 4** (P2 the shape confirmation, P5 the 15 C endpoint,
  P7 "a second course will not rescue the window") — and the misses are the entire beat.
- **Scope NOT entered:** `q_sugar_max` not moved; §2.2's window **not retired** (and §3 changes
  what retiring would MEAN — the honest repair may be to the benchmark's **TEMPERATURE**, a
  different edit from the one D-218/D-219 priced, and the owner's); 35-42 C recovered but the
  engine has **no high-temperature stress term** (that is WHY §5 crosses) and building one is
  **not licensed**; kveik and the Beer 2 control recovered and unused.
