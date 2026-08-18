---
name: prohibition-beer-aroma-calibration-levels
description: "D-224 — beer's seven calibrated aroma levels: what is settled, what must never be re-argued, and the one question left open"
metadata: 
  node_type: memory
  type: project
  originSessionId: f44fbd81-53a8-4e00-ac37-9275b8fd1ab0
  modified: 2026-08-18T22:44:16.620Z
---

# Beer's aroma calibration levels (D-224) — SETTLED

Reached by path from the ledger in [[project-fermentation-sandbox]]. Full record: `docs/DECISIONS.md` D-224.

## The verdict — do not re-open "which ester calibration is wrong"

**Neither.** D-223 §8's "one of the two is wrong and nothing here adjudicates it" is **ANSWERED**.
Isoamyl acetate is first-order in the `isoamyl_alcohol` pool (D-97 ATF1 coupling) and inherited a
**1.61×** error in its precursor; ethyl acetate's precursor is ethanol, which saturates the same
enzyme, so it did not. Put the precursor back on its own published mean and **both esters read 0.79×**
their targets — one common factor, the ferment-speed change D-223 shipped.

**Quote the 3 %, never the 0.4 %.** The 0.793 / 0.790 agreement is two disagreements cancelling: the
D-99 slacks disagree −2.35 %, the moves-since disagree +2.86 %. Defensible claim = *both moved by one
common factor 0.72-0.74, to within 3 %*. [[feedback-a-hit-can-be-two-errors-cancelling]]

## What is BUILT (7 values + 5 bands + 2 band rescalings) — never re-anchor without reading D-224

`k_propanol` 3.55e-4, `k_isobutanol` 3.41e-4, `k_active_amyl_alcohol` 3.65e-4, `k_isoamyl_alcohol`
1.13e-3, `k_2_phenylethanol` 9.12e-4, `k_ethyl_acetate` 1.39e-4, `k_isoamyl_acetate` 5.36e-4. All seven
land the level their own `conditions:` field states, **within 0.4 %**, at the calibration frame
**21 d / 20 °C / YAN 200 / pitch 1.0** — and all seven were inside their pre-existing bands, so no band
decision is hiding in a value.

**The targets are the STATED ones (20.0, 2.2, and Wang's four means), not what D-99 happened to land**
(21.32, 2.402). Absorbing that 6.6 % / 9.2 % slack is deliberate: it is what makes the two esters agree
**by construction** instead of by cancellation.

**Two figures already printed in `beer_generic.yaml` are RESTORED, not computed, by this:** the pool at
"~0.34 mM, ~87× below ATF1's Km" (shipped state: 0.547 mM / 54.4×) and "OAV ~0.6" against Meilgaard
(shipped: 0.965). That is the corroboration that 30.0 mg/L is the level the file was written for.

## The mechanism — the ONE sentence this reduces to

`FuselAlcoholsEhrlich` is **GATED on nitrogen (`N/(K_n+N)`) but draws its carbon from `S`**, so its
output is not **BOUNDED** by the nitrogen it is nominally made from: a slower-growing yeast holds the
gate open longer and makes more higher alcohol from the same YAN (every run ends at `N ≈ 0`). Therefore:

* **`mu_max` is the higher-alcohol knob** (band edges ×1.094 / ×0.774) and leaves ethyl acetate alone.
* **`q_sugar_max` is the ester knob** (×1.111 / ×0.897) and leaves the higher alcohols alone.
* **`isoamyl_acetate` reads BOTH.** That is why it looked like the disagreement.

**The separation is a property of the two KNOBS, not of the pools** — they share `S`. A 1.68× change in
the Ehrlich `k` moves packaged ethyl acetate **0.077 %** (0.28 % vs 0.19 % of sugar drawn), enough to
break D-223's corner pin. Never say the pools are decoupled.

## The drift, and that it was never scored

D-99 reproduced **all four** Wang/Frank/Steinhaus 2024 Table 1 beer means to **<0.5 %**. **D-211's
`mu_max` 0.098 → 0.034 multiplied every one by 2.87; D-222's → 0.058 halved it to 1.68. NEITHER RECORD
MENTIONS HIGHER ALCOHOLS** (grepped). Wine is the control and never moved (`mu_max` is per medium):
isobutanol 32.997 against 33.0, 2-phenylethanol 28.713 against 28.7.

## The sensory statement — this is why it is a DEFECT, not a preference

`isoamyl_alcohol` is the **only** beer pool of the five with a sourced in-matrix threshold (Meilgaard
1975, ~50 mg/L). The shipped model ran **48.261 mg/L, OAV 0.965** — within 3.6 % of claiming a solventy
note the file says a sound ale must not have — and **over the threshold (52.809) at `mu_max`'s low
edge**; 29 % of triangular draws on that `k` alone crossed it. **The repair moves the NOMINAL 0.965 →
0.602. It does NOT move the tail** — the corrected band still crosses at its top (33 %), because ×3 the
mean IS 90 mg/L and such ales exist. [[feedback-a-margin-is-a-claim-about-what-holds-it-open]]

## The band fork — the loser is named, do not re-argue it

D-99's five Ehrlich bands are ×0.3/×3 (propanol ×0.2/×5) of a centre **2.05× below** the nominal shipped
in the same commit, so the multipliers in force were **×0.145/×1.45** and the notes were false from
birth (both halves landed at 955ebbc; checked with `git log -S`). [[feedback-pin-the-band-not-the-nominal]]

* **SHIPPED:** edges = the stated multiple of the corrected nominal (isoamyl alcohol spans 9-90 mg/L).
* **REJECTED:** rescale the shipped edges, preserving ×0.145/×1.45 — it puts Meilgaard's 50 mg/L
  **outside the band entirely** (29 % → 0), i.e. asserts no ale can be fusel-y. Nothing sources that.
* **The two ESTER bands get the OPPOSITE treatment on purpose** (rescaled with the nominal, the
  D-97/D-99 convention) because their stated concentration span and their actual one agree. **Both rules
  are now tests** — do not apply one where the other belongs.

## The cost, and why the third ester option was refused

D-223's `Byp` funding constraint is exercised by a thinner slice: joint-corner flip fraction **5.37 % →
0.0767 %** (~1 draw in 1300; the corner still FORMS, +0.518 mg/L). **Anchoring ethyl acetate at Wang's
23.7 mg/L would take it to exactly 0.00 %** and leave that constraint invisible again — plus D-176's
independent reason (that survey mean carries a sour-beer tail). **Never move ethyl acetate to 23.7.**

**The 20 °C frame is load-bearing and now written down:** `E_a_esters` = 200 kJ/mol, so the same run at
15 °C lands ethyl acetate at **6.10 mg/L, below its own 10 mg/L floor**. 20 °C is a typical ale ferment;
§2.2's 15 °C is Foster's cool trial (D-221), a **different** frame. Duration is irrelevant past dryness.

## Guards (5 new) — the durable half

`test_the_finished_beer_lands_the_aroma_levels_its_rate_constants_are_defined_by` (equality against the
TARGET, never a snapshot), `..._each_drawn_speed_knob_moves_only_the_half_of_the_aroma_set_it_is_coupled_to`
(at the drawn band EDGES), `..._isoamyl_alcohol_stays_below_its_only_sourced_threshold_across_the_growth_band`,
and the two band-arithmetic tests. Falsified in 3 arms, designed GREEN in each, restore verified by
SHA-256. **Arm A missed its prediction and the miss corrected an over-claim in a test written that
hour** ("leaves ethyl acetate EXACTLY alone" → 0.03 %).

## OPEN — named, not licensed

**The coupling itself.** Beer's aroma pools are coupled to **biomass-hours**; both cited mechanisms are
**extent**-coupled — de Andrés-Toro 1998 forms ethyl acetate as `Y_EA·mu_x·X_A` (biomass FORMED,
nitrogen-limited ⇒ invariant to `q_sugar_max`), and the Ehrlich pathway's substrate is amino acids
(⇒ invariant to `mu_max`). **Under an extent-coupled rate law none of the seven levels would have moved
at D-211/D-222/D-223 and D-224 would not exist.** NOT licence to build it: the flux coupling is
D-19/D-21/D-96/D-97/D-99's, carries the sourced temperature ORDERING, and reaches wine too. Its own beat.

Also open: **Wang Table 1's beer means for the two ACETATE esters were never sought** (paper not on
disk); `q_sugar_max`/`mu_max` are NOT re-visited — they are fitted to measured courses and the aroma
levels are the consequence.
