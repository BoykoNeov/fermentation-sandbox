---
name: prohibition-beer-criterion-temperature
description: "D-221 - §2.2's beer attenuation criterion is re-temperatured 20 -> 15 C and is now a strict xfail; E_a_uptake's decoupling freedom is dead, and the two literature anchors start agreeing"
metadata: 
  node_type: memory
  type: project
  originSessionId: 47cfa586-43a7-4fdb-9cf6-b1f266627a26
  modified: 2026-08-18T15:15:03.651Z
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

# §2.2's beer criterion: 15 °C, xfail, and what the temperature change killed (D-221)

**SETTLED. Do not re-propose any of this as unbuilt, and do not "fix" the xfail by moving a
temperature back.** Read this file before touching the beer attenuation criterion, `E_a_uptake`,
`E_a_growth`, `T_ref`, or any claim that something is "inert at the benchmark".

## The edit

`BENCHMARKS["beer_attenuation"]`'s `conditions` moved **20 → 15 °C**. **`low`/`high` are the
handoff's own 5-7 d and did NOT move** — D-220 found that duration real at 15 °C (Foster's three
commercial ale controls: 5.06-6.26 d, conservatively re-derived at D-221 from the 720-point
recovery) and impossible at 20 °C (2.91-3.77 d at 22 °C). **Foster's measured band sits INSIDE
5-7**, so keeping the wider window is the reading generous to the model.

`test_beer_1048_og_attenuates_in_5_to_7_days` is now **`xfail(strict=True)` at 9.000 d** — a
2.74 d miss against Foster's own conservative slow strain. **Never "fix" this by narrowing the
window to Foster's band** (a tightening nothing mandates) **or by re-temperaturing back**. The
20 °C pass was a criterion calibrated to a cooler ferment certifying a ~1.5× too slow model.

Reading is grid-stable and is the **quieter** of the two frames: 0.0052 d spread over 1-8
points/h, against 0.0365 d at the retired 20 °C.

## What DIED, and must never be cited again

**D-216 §6's decoupling lever is GONE, not weakened.** `E_a_uptake` was "the ONE parameter that
can move Tyrell's 15 °C course without moving §2.2" **only because the two anchors sat at
different temperatures**. `TYRELL_TRIAL_CELSIUS` is **15.0 exactly** and the criterion is now
15.0, so one Arrhenius factor enters both and the parameter moves them together.

| frame | printed 30,000-63,000 J/mol band moves the criterion by |
|---|---|
| 20 °C (retired) | **0.0000 d** — exact, by construction (`T_ref`) |
| 15 °C (live) | **1.75 d** — against a window only **2.0 d** wide |

**There is no magnitude argument to retreat to.** Most inert → one of its strongest levers.
The low band edge does not rescue the criterion either (30,000 → 7.69 d, still outside).
`E_a_growth` lost the same claim more weakly: 0.0000 → **0.2917 d**, a sixth of uptake's.
**All three yaml notes (`E_a_uptake`, `E_a_growth`, `T_ref`) are corrected and both spans are
now ASSERTED**, so the notes cannot drift from the code.

**D-217's refusal to re-source `E_a_uptake` is UNTOUCHED** — it rests on the corpus having
nothing, not on this lever. `T_ref` itself is unmoved: it is where the rate constants are READ,
a property of the fits, not of any acceptance criterion.

## The conflict INVERTED in direction

| frame | admissible `q_sugar_max` for 5-7 d | shipped 0.5 | Tyrell-matching 1.397 |
|---|---|---|---|
| 20 °C (retired) | 0.4250 – 0.6214 | 6.05 d PASS | 2.69 d |
| 15 °C (live) | 0.6670 – 1.0165 | 9.00 d FAIL | 3.99 d |

At 20 °C the criterion **forbade** a faster engine; at 15 °C it **demands** one. Beer's two speed
anchors no longer disagree about the **sign** — both want the engine faster — only about
magnitude, by ~1.37×. **D-216 §11's open question is now a magnitude question, never a direction
one.**

**D-216's refusal of 1.397 SURVIVES its own premise being spent** (blocking check, not an
assumption): 1.397 stays outside across the **ENTIRE** printed `E_a_uptake` band (30,000 → 3.58,
63,000 → 4.17), reaching 5-7 only at ~120,000 J/mol, 1.9× the high edge. **Nothing in D-221
licenses a rate change**; D-216 §4's shape objection is untouched.

## D-218 §3 and D-219 §5c INVERTED

| per-cell mass | `q` | §2.2 at 15 °C | §2.2 at 20 °C |
|---|---|---|---|
| 18 pg (retired) | 1.500 = ceiling | 3.83 | 2.58 |
| **40 pg (SETTLED)** | 0.9242 | **5.42 INSIDE** | 3.63 |
| 50 pg | 0.8176 | **5.96 INSIDE** | 4.00 |
| 100 pg (retired) | 0.5602 | 8.17 | 5.50 INSIDE |

**The surviving corner moved from the RETIRED reading to the SETTLED one.** D-219's cell-mass
settlement and D-221's temperature correction agree where the old pairing had them in conflict.
D-218 §3's "survives only if 72 h is read as exact" also inverts: the open-end arm goes
3.71 d (outside) → **5.54 d (inside)**.

**The agreement is CONDITIONAL and must always be quoted that way.** Swept, not spot-checked:
inside from **40,165 J/mol** up — **69 %** of the printed band, shipped 55,100 inside it —
and overshooting below (30,000 → 4.71 d). Say *"compatible over most of the band including the
shipped value"*, **never "compatible"**.

**None of it says the model passes** — all four rows are RE-RATED engines. Tyrell's extract
schedule stays incompatible with both. **Never cite §5's agreement as a temperature-response
result**: the `E_a_uptake` band moves that same reading 1.75 d, so it discriminates nothing.
D-218's and D-220's rule against reading Foster's temperatures as a response test STANDS.

## Falsification and process notes

Four in-band arms, hash-verified restore between each, designed GREEN (`E_a_esters` → 210,000).
**`T_ref` → 288.15 (its printed low edge = 15 °C, which would restore the inertness) kills seven
of eight** — the strongest evidence §4 is about the frame, not a number.

Pre-registration: **2 hits, 4 misses**. P3 predicted one test would move; **six** did. P5
predicted 1.397 would land inside; it did not, and the pre-registered consequence is what the
record asserts. P2's 0.246 d frame gap is a **different wort**, not a different helper —
cross-frame beer durations are not comparable below ~0.25 d.

`test_the_uptake_activation_energy_is_a_lever_only_because_the_trial_ran_cool` stayed **GREEN
with a false docstring** — it never reads the criterion. Corrected in place.
See [[feedback-a-freedom-can-be-an-artefact-of-the-frame]].
