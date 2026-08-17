---
name: trub-settling-and-the-peptide-pair
description: "D-214 - trub settling is REFUSED as a fermentation-phase term (pre-pitch, already in the calibration, charge-violating after the anchor), and the peptide capacity/pKa pair is measured incoherent but deliberately not fixed"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — trub protein settling and the peptide pair (D-214).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about trub, beer's protein buffer, `peptide_buffer_capacity_beer`, or
K⁺/H⁺ antiport. Every bullet is *what it forbids* + the record to read for *why*. **If a
prohibition looks unconvincing, go read D-214 — do not argue past it from this file.**

**REFUSED and MEASURED. D-209 §8's two parked terms are now BOTH closed — do not re-propose
either as unbuilt, and do not re-park them.**

- **K⁺/H⁺ antiport is DEAD ON SOURCING and was dropped before probing.** `antiport` returns
  **ZERO hits in every beer text in the corpus**; its only 8 hits are lactic-acid-bacteria
  malate/citrate antiporters in WINE books — different organism, different reaction. D-209 §8's
  wort K⁺ = 550 mg/L is a **pool size, not a flux**. Building it invents the rate constant that
  governs the timing of the curve it would be fitted against — D-213 §7's refusal, one beat later.
- **Trub is a PRE-PITCH event and is ALREADY INSIDE THE CALIBRATION.** Peyer's 1.34 → 1.21 drop is
  **at the boil**; cold break is **at chilling**; *Chemistry of Beer* §2.9: the coagulate is
  *"removed before the wort is fermented"*. `peptide_buffer_capacity_beer` is back-solved from
  Peyer's **1.18 CONTROL WORT — already post-boil**. **Never call trub an omission.** No text in
  the corpus measures buffering capacity after pitching; the 4 beer texts do not discuss BC at all.
- **A fermentation-phase Process draining `peptide_buffer` is a CHARGE-BALANCE VIOLATION, not a
  small acidification.** That pool rides `_BEER_ACID_SEEDS`, so the t=0 cation back-solve is fitted
  **with it present as the counter-anion**, and it is the biggest one (organic acids alone are ~1/6
  of wort's buffering). Cutting 20 % at 6 h takes day 1 to **pH 7.08**; cutting all of it takes day
  7 to **11.66**. Day 7 is IDENTICAL for a 6 h and a 24 h cut — the effect is on the **permanent
  charge**, so it has no time profile at all. `test_no_process_touches_the_peptide_buffer_pool`
  guards this and was **falsified before shipping** (RED naming `wort_acid_removal`).
- **D-209 §8's "same-sign acidification" classification is CORRECTED.** Pre-anchor the sign is
  right; post-anchor it is **opposite**. Both halves were named in one clause and they disagree.
- **Pre-anchor it is refused on SHAPE, and that is the cleanest kill** — it needs no magnitude the
  corpus withholds. A 20 % loss moves day 1 by **−0.0172** and day 7 by **−0.0587**: the effect
  **GROWS 3.40×** with time, because removing buffer *amplifies cumulative acid production* rather
  than adding acid. D-211 §9's brief wants **early and none late**; this is the exact inverse.
- **The window is EMPTY at the arm that parks it, by ~9×.** At the HIGH `z̄` edge day 1 needs a
  loss **≥ 27.59 %** while day 7 affords **≤ 3.09 %**. Never answer this from the nominal.
- **The peptide PAIR is incoherent off-nominal — MEASURED, RECORDED, NOT FIXED, and that scope
  line is deliberate.** The capacity is pinned "because it is derived from `pKa_peptide_buffer`",
  but the pKa is read at RUNTIME and drawn (12 distinct values / 12 members) while the capacity is
  a COMPILE-time seed solved at the nominal pKa only. So a member carries BC **1.1161-1.180**, not
  1.18. **It cannot move the parking verdict**: **+0.0099 pH** on the low-pKa arm — already
  **0.145 BELOW** the day-7 floor, out of band whatever the capacity — and **exactly 0.0000** at
  the nominal, so pricing it against the 0.0086 headroom **crosses two scopes**. **It needs its OWN
  beat — never let it ride out on a refusal's immunity from review.**
- **`acidbase.yaml`'s capacity triple was STALE and is corrected**: 10.75/10.28/10.45 mM was
  computed on the FIVE-acid table, before D-181 moved the nominal 11.36 → 10.52. Re-measured
  **10.99 / 10.52 / 10.70**; the "stays stable" claim survives (4.5 % vs the quoted 5 %).
- **D-211 §9's two numbers are CORRECTED to 0.0274 (day 1 above ceiling) and 0.0086 (day-7
  headroom).** They were an `argmin` read over the solver's ADAPTIVE output, which landed on
  **t = 23.6382 h** — 22 min early on the steepest part of the curve, reproducing 0.034/0.0082 to
  the digit. **Verdicts unchanged.** **No shipped test was affected** — `_tyrell_degassed_ph_at_day`
  uses `np.interp` onto the exact hour; the artefact lived in the probe.
