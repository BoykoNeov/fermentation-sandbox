---
name: beer-early-acetic-rise
description: "D-212 - beer's early acetic shortfall is REFUSED; the day-1 pH admits an acetic WINDOW that Tyrell's own measured 145 sits outside, and four candidates die of four different causes"
metadata: 
  node_type: memory
  type: project
  originSessionId: 319bf67b-4e83-4c11-9ffb-068bbcb68911
  modified: 2026-08-17T15:40:51.610Z
---

**Live prohibitions — beer's EARLY acetic rise and the day-1 pH window (D-212).** Detail split
out of `.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by
path. Read it before proposing anything about beer's acetic timing, beer's wort oxygen, or
beer's glycerol. Every bullet is *what it forbids* + the record to read for *why*. **If a
prohibition looks unconvincing, go read D-212 — do not argue past it from this file.**

**The beat is a REFUSAL and it was PRE-REGISTERED as 55 % likely to be one. Do not re-propose
it as unbuilt work.**

- **NEVER aim at Tyrell's measured day-1 acetic of 145.0 mg/L.** The day-1 pH admits an acetic
  **WINDOW** — **94.24-141.22 / 90.20-136.82 / 86.11-132.38** mg/L at the lo/nom/hi arms of
  `nitrogen_uptake_charge_beer` — and **145.0 is ABOVE it at every arm** (by 3.78/8.18/12.62).
  The model's 80.13 is short by only **+6 to +14**, a **fifth** of the +64.87 the acetic curve
  implies. **A build that hit the acetic target would score as a success while moving the pH
  OUT of band.** Both numbers come from the same figure of the same paper.
  [[feedback-a-hit-can-be-two-errors-cancelling]]
- **It is NOT a general day-1 acid-profile problem** — checked, not assumed: lactic 111.9
  (87-142), succinic 56.6 (32-76), malic 99.8 (81-116) and all three falling acids are inside
  their day-7 bands, `citrate` static at 205.
- **An ADDITIVE early source is forbidden — it breaks day 7.** D-183 refused the sink, so every
  mg/L made early survives: **117.75 + 64.87 = 182.6 against a measured 105-126.** Only a
  **redistribution** of the existing producer preserves the total by construction. Its cap is
  **f = 2.781** (the whole 58.75 mg/L rise by day 1, producer silent after); entering the band
  needs **f = 1.668/1.477/1.283**, so **any f in 1.67-2.78 works at all three arms at once.**
- **The refusal is on MECHANISM, not identifiability — do not cite D-203/205/206's reason.**
  The verdict is *insensitive* to the split over most of its range (contrast D-205, where
  3.4 vs 3.0 gave "two oxidised wines or none"). The pH half is **buildable and robust**; what
  is missing is a **named driver**. Inferring the shape from Tyrell's own curve is
  [[feedback-fit-the-observable-not-the-consequence]].
- **Four candidates, FOUR DIFFERENT deaths — never collapse them to "unsourced".**
  **(a) O₂/PDH-bypass** dies on an **unseeded driver**: beer's `o2` is **identically 0.000** and
  its three Processes are **all non-default** — verified in the REGISTRY, not one run.
  **(b) ALD/acetaldehyde** dies on **identifiability, NOT mass** — throughput **~3038 mg/L** over
  day 0-1, gap needs **1.57 %**; the pool argument (32.2 vs 47.7) is **UNSOUND**, and the pool
  peaks **day 3** anyway. **(d) a bare lag term** is curve-fitting.
- **(c) glycerol/redox is the coupling the corpus DOES support** (Handbook Vol 1: GPD1 with
  ALD2/ALD3, acetate regenerates NADH) and it is blocked by **D-16's deliberate
  `Y_glycerol_sugar` = 0.0**, protecting the attenuation and CO₂-ratio benchmarks — a **SCOPE
  decision, NOT a sourcing block.** Its sources are **high-sugar musts**, not a 12 °P wort.
- **The beer corpus has NO yeast-side acetic timing** — all four beer texts give spoilage
  bacteria (*Acetobacter*, *Pectinatus*, lambic). **A LOCAL-CORPUS NULL, recorded as one**
  ([[feedback-paywalled-is-one-host]] — wrong 6× so far).
- **The NEXT beat is named and its anchor is SOURCED, but it is NOT started.** *Craft Beers*
  gives wort O₂ at **5.5-8.0 mg/L** *and* a direction on this beat's observable ("increased
  aeration causes … **a faster pH drop**, more acetaldehyde…"), predicting the sign of **five**
  quantities already modelled — a pair, so it constrains a **response**
  [[feedback-a-pair-constrains-a-response]]. **What it does NOT give is an O₂→acetate
  stoichiometry**, and that is what would make it a build.
- **Three SCOPE questions belong to the owner, not to a beat** (D-212 §7): seeding beer's O₂
  (a new default-set surface — the `o2` Processes are non-default, so seeding changes what the
  oxidative cascade does in beer); un-zeroing beer's glycerol (**benchmark-touching**, D-16);
  and the **wine** side of the same coupling, which **IS expressible today** (`Y_glycerol_sugar`
  0.035, its own note already conceding *"glycerol is front-loaded with growth in reality"*) —
  **recorded, not begun.**
- **D-183 is unchanged**: its producer choice still wins (RMSE **40.7 vs 65.3**, re-measured
  here) and its refused re-assimilation sink **STANDS** — this beat leaned on that refusal
  rather than reopening it. **D-211's flag on D-183 stays OPEN**: the 2.15× day-1 shortfall is
  now measured, priced and explained as far as the corpus allows, still not built.
- **D-211 §13's sweep claim is CORRECTED.** *"A sweep … returns nothing else"* missed
  `Y_acetic_biomass_beer`'s own note — the parameter for the very acid its `Flags:` names — which
  held **two** numbers conditioned on the retired growth rate: the plateau day (**1 → 3**) and
  the RMSE pair (**61.6 → 32.5** becoming **65.3 → 40.7**, a ~38 % cut, not a halving). Both
  repaired; no shipped value moved, which is why a green suite could not catch them.
