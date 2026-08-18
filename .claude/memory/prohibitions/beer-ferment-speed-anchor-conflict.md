---
name: beer-ferment-speed-anchor-conflict
description: "D-216 - beer's early-limb ferment lag is REFUSED as a parameter fix: the knob is in band but a second anchor forbids it, not even removing repression reaches the measurement, and the pH agreement turns out conditional on the scenario pitch"
metadata: 
  node_type: memory
  type: project
  originSessionId: f6587b9b-5f34-4f36-af85-d37c405530dc
  modified: 2026-08-17T19:55:35.490Z
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

**Live prohibitions — beer's fermentation SPEED and the pitch its pH rides on (D-216).** Detail
split out of `.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by
path. Read it before proposing anything about beer's uptake rate, catabolite repression, the
§2.2 attenuation benchmark, or the pitch in `TYRELL_SCENARIO`. **If a prohibition looks
unconvincing, go read D-216 — do not argue past it from this file.**

**SUPERSEDED IN PART BY D-222 — read this block first.** Three things below are now FALSE.
1. **The pitch prohibition is SPENT and the pitch is CORRECTED.** `TYRELL_SCENARIO` now carries
   Tyrell's own counted 9.96e6 cells/mL = **0.3984 g/L** through `cells_per_ml_to_pitch_gpl`.
   The "two observables endorse 1.0 g/L" argument was measured **without refitting `mu_max`**,
   which D-216 itself named as inherited: refit at the corrected pitch (0.034 → **0.058 /h**,
   band 0.031-0.040 → **0.053-0.075**), the nitrogen attribution returns to 0.298 (inside
   Tyrell's 0.234-0.448) and the pH course returns to **7 of 8** days inside. D-216's 6-of-8 was
   measured at 0.5 g/L, a reading D-219 retired; at the counted pitch without a refit it is 5/8.
2. **Growth EXTENT does NOT "own 0.8 %" and is not near.** That was a pitch-1.0 statement. The
   gain is nitrogen-limited and therefore fixed in absolute terms (~1.75 g/L either way), so a
   lighter inoculum multiplies further: **5.385× against a measured 2.92-3.48×, 1.55× above the
   high edge**. No growth RATE can repair it (the fold moves 0.010 across the whole band) — the
   candidates are the assumed 200 mg/L YAN and the single-lumped-pool assumption, neither
   sourced. Do not re-run this as if D-211 had it right.
3. **"1.397 is INSIDE the band so out-of-band is NOT the reason" is DEAD.** Re-bisected at the
   counted pitch the Tyrell-matching `q_sugar_max` is **2.3226 = 1.55× the printed high edge**.
   The refusal got SIMPLER: the knob is inadmissible on its own band before the second anchor is
   consulted. The second tier moved the same way — removing `K_repression` entirely closes 46 %
   of the day-2 gap (was 79 %) and puts the benchmark at 4.54 d.
The day-2 shortfall is now **4.21×** and §2.2's criterion reads **8.50 d** at its corrected 15 °C.

**THE LAG IS NOT D-211'S DOING, AND NOT THE GROWTH EXTENT.**
- Reverting `mu_max` to the pre-D-211 0.098 gives day 2 = **0.289** against a measured 0.594.
  The lag **pre-dates D-211** — 2.05× before, 2.81× after. D-211 worsened an existing defect and
  is **not impugned**; its fit was to measured cell counts. Never re-open `mu_max` for this.
- Growth EXTENT owns **0.8 %**. Measured, not argued: the fold sits 6-7 % below Tyrell's
  envelope at days 2-3 (D-211 called that "near"), and forcing it mid-envelope via YAN 255 moves
  day-2 flux **0.212 → 0.215**. Do not re-run this as if it were open.

**THE FIX IS REFUSED, AND THE REFUSAL IS TWO-TIER. Do not re-propose either knob.**
- `q_sugar_max` = **1.397** reproduces Tyrell's day 2 and is **INSIDE** the printed 0.3-1.5 band
  — so "the value is out of band" is NOT the reason. The reason is §2.2's benchmark, which the
  same constant sets: it breaks at **q ≈ 0.6**, having closed under a fifth of the gap, and lands
  **2.71 d** at 1.397. **No in-band (q, K_repression) pair satisfies both anchors**; best
  benchmark-passing Tyrell day 2 is 0.269 against 0.594.
- `K_repression` = 2.0 is `speculative`/"author estimate" and owns **79 %** of the lag (removing
  it takes day 2 to 0.514) with the right SHAPE. **Re-sourcing it is still not enough** — the
  unbounded LIMIT falls short of 0.594 and puts the benchmark at 3.42 d. This tier exists
  precisely to kill *"re-source the constant and the lag goes away"*.
- The generous corner (every rate parameter at its fastest in-band edge, benchmark held at
  exactly 5.00 d) leaves **1.79×**, and that is a **lower bound**: it uses `E_a_uptake` = 30 kJ/mol,
  the one lever that **decouples** the anchors (the benchmark runs at exactly `T_ref`, so its
  Arrhenius factor is 1.0 by construction — free — while Tyrell at 15 °C is scaled, worth
  2.81× → 2.41×). **That edge is REFUSED**: it is "retained from the now-debunked ~35 kJ/mol beer
  figure" (D-19). Named, measured, not adopted.

**~~THE PITCH IS LOAD-BEARING — NEVER "CORRECT" IT~~ — RETIRED AT D-222, see the block at the
top of this file. Kept because the reasoning is instructive, not because it is live.**
- `TYRELL_SCENARIO` pitches **1.0 g/L** against Tyrell's 9.96e6 cells/mL (~100 pg/cell, ~2× the
  textbook 40-60). At an "honest" 0.5 g/L the extract lag gets **WORSE (2.81× → 3.51×)** and the
  pH course drops **7/8 → 6/8** days inside, with D-211's pinned day-1 miss going **0.070 → 0.354**
  and its nitrogen attribution **0.363 → 0.181, outside** the 0.234-0.448 spread it cites as what
  makes the timing "MEASURED rather than fitted". **D-211 is FLAGGED, not corrected** — two
  observables endorse 1.0 g/L against the per-cell arithmetic. What is forbidden is reading
  D-211's 0.070 as unconditional.

**SUPERSEDED IN PART BY D-220 — read `beer-second-measured-course.md` before using this block.**
Two claims below are now FALSE: there **IS** a second measured beer extract course on disk (Foster's
Supp. Fig. S1, 8 temperatures), and §2.2's window is **not** simply "the anchor that is not a
measurement" — the brief's DURATION is corroborated by a published trial at 15 °C, and it is the
pairing with 20 °C that fails. **The refusal to move `q_sugar_max` STANDS** — it rests on the
benchmark and on the shape objection, neither of which D-220 spends.

**SCOPE, and what was NOT measured.** There is **no second measured beer extract curve on disk** —
Zamudio Lara 2022 is cited all through `beer_generic.yaml` but is **not in the local corpus**. So
nothing is claimed about worts other than Tyrell's and D-215's SCOPE paragraph stands. The frame
was checked and is clean (the extract panel is *apparent* extract, already converted at D-183 §2).
Days 3-4 of the measured course live in D-215's §4 table, **not** in the shipped constant.

**What shipped: three tests, no `src/` or parameter change.** The pitch test is honestly a
**conditional PIN, not a mechanism guard** — doubling `q_sugar_max` moves its number only
0.354 → 0.323, so it fires on pin drift. Also learned: **a mutation arm must be IN BAND**, because
the parameter store validates at load and an out-of-band arm dies in pydantic
([[feedback-verify-an-xfail-fails-for-its-stated-reason]]).

**OPEN, and the owner's call — not a default.** The two anchors conflict and the model sits at
the one that is not a measurement: §2.2's 5-7 d window is an acceptance criterion from the
handoff brief (which `CLAUDE.md` calls *"reference, not gospel"*), Tyrell's course is a published
trial. Which one beer's kinetics should be calibrated against is **recorded, not taken**.
[[feedback-a-gap-can-be-held-open-by-a-second-anchor]]
