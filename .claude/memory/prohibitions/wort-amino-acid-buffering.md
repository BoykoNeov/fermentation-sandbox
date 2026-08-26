---
name: wort-amino-acid-buffering
description: "D-239 - beer's three free amino-acid side chains BUILT as a reading of the N pool, the peptide lump re-partitioned at constant total, and the day-7 high band edge deliberately outside Tyrell's envelope"
metadata: 
  node_type: memory
  type: project
  originSessionId: 80611b07-f4ed-492d-8b4b-03c3866a0e21
  modified: 2026-08-26T18:41:56.174Z
---

**Live prohibitions — the wort's free amino-acid buffering (D-239).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about beer's buffering capacity, the peptide lump, D-209 §8's
"buffer-removal half", or beer's day-7 pH agreement. Every bullet is *what it forbids* + the
record to read for *why*. **If a prohibition looks unconvincing, go read D-239 — do not argue
past it from this file.**

**D-209 §8's BUFFER-REMOVAL HALF IS BUILT. Never re-propose it as unbuilt or as a candidate for
beer's day-1 miss.**
- **The three species are Asp/Glu/His SIDE CHAINS ONLY, monoprotic, at Peyer's own printed pKas
  (3.86 / 4.25 / 6.04).** The α groups are **DECLINED, not deferred** — the source prints them
  only as RANGES over "most amino acids" (1.7-2.2, 8.8-10.6) and a per-species value invented
  from a range is the transcription error. The other 15 amino acids and these three's α groups
  stay inside the lump, one rule applied uniformly. **Never root the split on all 18**: 57 % of
  that effect is buffering below pH 3 (Peyer's titration ends near 3.0) and ~6 % in the vessel —
  it is D-178's phosphate refusal applied backwards, and it costs 0.056 pH instead of 0.020.
- **NOT a state slot and NOT a Process.** Third "include-by-reading" entry after `Byp` (D-18) and
  `carbonic` (D-182); concentration is `ratio · [N]`, so the pool drains for free and **no
  nitrogen is booked twice**. `_totals_molar` now REQUIRES `params` (D-182's forcing function) —
  an omitted term does not raise, it returns a wort that buffers like a finished beer.
- **The charge is a RE-PARTITION, never new charge.** `nitrogen_uptake_charge_beer` 0.1772 →
  **0.23418212124814006** (both edges by the same constant; **the band WIDTH is unchanged** — it
  is the ammonium range and this is a composition). At wort pH the two halves cancel term for
  term to D-209's **0.1772 exactly**. **Never compare 0.234 with D-209's 0.177 and call it
  drift** — the comparable quantity is the NET, and a test holds it. The lump re-roots
  **1.5480662315921656 → 1.4350620340127729** with the wort still at Peyer's 1.18 exactly.

**THE COST IS PRICED AND THE DAY-7 HIGH EDGE IS OUTSIDE ON PURPOSE. Do not "repair" it.**
- Δ pH vs a coherently unwired arm: **0.000000 at t=0** (the anchor absorbs it — asserted as
  EXACT), **+0.0023 at day 1**, −0.0140 day 2, **−0.0202 day 7**. **The day-1 sign is POSITIVE
  and that CORRECTS D-209 §8's "same-sign, can only make the acidification larger"**: only 0.298
  of the nitrogen is drawn by 24 h, so the pool is still **70 % present** and buffering while the
  lump's share of it has gone. **It depends on the uptake calendar and would invert on D-209's
  own >99 %-by-24 h** — re-measure the sign if beer's uptake timing ever moves again.
- **`nitrogen_uptake_charge_beer`'s HIGH edge now finishes day 7 at 4.7625, 0.0203 BELOW Tyrell's
  floor** (5/8 days inside; low and nominal stay inside at 7/8, nominal headroom **0.0046**).
  That edge's margin history: 0.003 (D-209) → 0.036 (D-222) → 0.0033 (D-223) → spent. The test is
  RENAMED for what it forbids and keeps every tooth. **NEVER convert this to an xfail**: D-208's
  idiom is for something true of the source and false of the model; here the model got MORE
  faithful and the agreement got worse, so a different term is missing (D-232's growth-extent
  residue is the standing candidate, on the alkaline side).
- **Beer's day-1 miss is now OPEN WITH NO NAMED CANDIDATE.** D-222 reserved its +0.172 for this
  term; it took **0.008**, with the sign against it on day 1.

**What else moved — re-measure, never cite the old number.**
- Predicted pH-drop band **94.1-129.2 %** nominal / **73.4-142.5 %** joint (was 91.4-127.1 /
  71.2-140.1). Peptide-pKa leverage FELL (BC 1.1203/1.1487 vs D-233's 1.1161/1.1456) — less of
  the wort is in the lump. **D-214's window arm moved off the high edge** (its baseline is now
  outside; a bracket there would score a margin against a deficit): affordable day-7 loss at the
  NOMINAL is **1.72 %**, and D-214's 3.1 % / D-222's 12.6 % / D-223's 1.2 % are HIGH-edge numbers,
  **not continuous with it**. Day 1 stays a SATURATION at every edge, so trub stays refused.
- **The two pitches SEPARATED** — 7 days inside at Tyrell's counted pitch vs 5 at the retired one
  (D-222 measured 7 at both). D-222 §6's stronger claim (the course cannot referee the pitch at
  all) is RETIRED; the claim it protects survives because the day-1 miss still CROSSES ZERO.
- **The calibration wort gained a NITROGEN coordinate** (`peyer_control_wort_yan_gpl` =
  0.19147259991812246 g N/L). Tyrell's 200 mg/L legitimately roots **0.35 % lower** (1.4300292930172551).
  **Never re-root the shipped literal at Tyrell's YAN** — rule 3 re-roots per member, which is how
  a wort that is not Peyer's still reaches 1.18.

**Falsification, and the one thing it repaired ahead of the band.**
- Arm A (ratios zeroed, `z̄` left high) turns the identity guard RED at **0.234182** — that arm is
  what makes it non-vacuous. Arm B (trio out of the compile anchor) gives **pH 5.3448 vs 5.65** —
  and my predicted DIRECTION was wrong, the anchor under-supplies cation. Arm C (gate removed)
  gives an un-anchored beer **pH 3.644**, the mirror of D-179's pH-11 artefact.
- **Rule 3's exact-nominal skip WAS keyed on `pKa_peptide_buffer` alone** and now compares every
  name the back-solve reads. Measured before closing: banding `wort_aspartate_per_n` ±5 % put 8
  members at BC **1.1831-1.1883** with every seed at the nominal — D-233's defect re-entering by
  another door. After: **1.180000000** for all eight. All six D-239 names ship **PINNED**, so it
  is bit-identical today. **The ratios are a LOCKED PAIR with `nitrogen_uptake_charge_beer`'s band
  (2.61 %, shared ammonium denominator) — band one and you must band the others.**
- **WINE IS UNTOUCHED and that is MEASURED, not preferred**: the same three species in a must are
  **0.73 %** of its acid buffering vs beer's 6.7 %, and wine speciates amino acids as slots
  already (D-100) while keeping them charge-inactive. Closing it there is a different act.
  The ionic-strength caveat on all three new pKas is **UNPRICED**.
