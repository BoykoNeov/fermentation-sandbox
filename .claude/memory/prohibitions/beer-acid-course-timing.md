---
name: beer-acid-course-timing
description: "D-215 - the corpus holds no missing beer acidifier; the three flux-linked acid timing errors OPPOSE and cancel; and the engine ferments Tyrell's wort ~2.8x too slowly, unscored until now"
metadata: 
  node_type: memory
  type: project
  originSessionId: c639a01c-94f5-44c3-a91b-491833d1c5c9
  modified: 2026-08-17T18:41:36.380Z
---

**Live prohibitions — beer's acid COURSES and the extract schedule (D-215).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about beer's flux-linked acid producers, beer's fermentation speed, or
"what the literature says is missing" from beer's pH. Every bullet is *what it forbids* + the
record to read for *why*. **If a prohibition looks unconvincing, go read D-215 — do not argue past
it from this file.**

**THE CORPUS HAS NO MISSING BEER ACIDIFIER. Do not re-open the search.** All 24 texts swept.
- **Calcium phosphate precipitation is a BOIL effect** (*The Chemistry of Beer*: 5.6-5.8 → 5.2-5.4,
  *"largely because of the precipitation of calcium phosphate"*). **Pre-pitch, inside Peyer's
  control-wort anchor — exactly like trub at D-214.** It is the mechanism a brewer names first and
  it is **not an omission**. Never propose it as one.
- **Yeast H⁺-ATPase / proton extrusion is a LOCAL-CORPUS NULL.** Every hit is lactic-acid bacteria
  or wine; the only two *S. cerevisiae* hits are **amino-acid/H⁺ symport**, already shipped as
  `nitrogen_uptake_charge_beer`. **I reached for this from general recall and the corpus refused
  it** — do not restore it from memory either.
- **No beer text treats buffering capacity during fermentation** (one hit, an interview about water
  hardness). Peyer's thesis is still mash/boil/wort only, as D-214 recorded.
- **Tyrell Fig. 4 is still the ONLY beer pH curve on disk.** Other texts give endpoint levels
  (3.83-4.49) and one mash statement (5.0 → 4.3). Never claim a second curve without producing it.
- **Autolysis is the one unbuilt mechanism named — and it is spoiled, not pending.** Tyrell say it
  raises organic acids late; *Craft Beers* says it **raises pH** — **opposite signs** — and that
  passage is weeks-scale cone yeast against a 7-day window. `YeastAutolysis` is wine-only/disabled
  ("beer deferred"). `Y_lactic_sugar_beer`'s note already settled it: *"one dataset cannot separate
  a late excretion from an autolytic release."* **If ever built, BOTH halves ship together** — the
  acid release AND the amino-acid release that pushes pH the other way ([[feedback-gate-both-halves-of-a-pair]]).

**THE THREE `Y·ΔS` ACID TIMING ERRORS OPPOSE AND NEARLY CANCEL. Never "fix the shape" with one knob.**
- Fraction of each acid's day-0→day-7 rise done **by day 2**: **succinic 45.9 % measured vs 20.5 %
  modelled (LATE by 25 pts)**, **malic −4.1 % vs 20.5 % (EARLY by 25 pts)**, lactic 14.6 % vs 20.5 %
  (early by 6). All three share one rate law so the modelled fraction is **identical by
  construction**. Forcing all three onto measured nets **−0.0083 pH** at day 1 because malic
  (+0.0182) and succinic (−0.0205) annihilate. **There is no single correction that helps all
  three**; one that made them agree would be fitting a compromise, not a mechanism.
- **Malic is NOT non-monotone at the four-strain mean** — I claimed that and it was over-read. The
  mean dips 2.25/1.0/1.25 ppm against a **±2 ppm read tolerance**; only strain 15 falls clearly, and
  it is the strain recorded as producing essentially none. **No malic sink is warranted.** (Beer has
  none; the only `malic` consumer anywhere is wine's pitched `MalolacticFermentation`.)
- **The pH consequence is UNRESOLVED and must not be quoted.** ±2 ppm of the acids is worth
  **0.0316 pH at day 1 / 0.0266 at day 7** against a day-7 headroom of **0.0082** — the noise floor
  is **3.2×** the headroom and **3.8×** the effect. Reported inconclusive by pre-registered rule.

**THE ENGINE FERMENTS TYRELL'S WORT ~2.8× TOO SLOWLY, AND IT IS NOT THE DAY-1 pH CAUSE.**
- Fraction of fermentable consumed — model **8.6 / 21.2 / 37.5 / 56.0 / 72.8 / 93.1 %** on days
  1/2/3/4/5/7 against Tyrell's **15.0 / 59.4 / 80 / 90 / 100.3 / 99.7 %**. **Nothing had ever scored
  it**: Fig. 4's three panels were each read by a different beat (D-180 extract *endpoints*, D-207
  pH, D-211 cell counts) and no test asked whether the engine ferments on the source's schedule.
  **Total attenuation is INSIDE §2.2's 5-7 d window (6.08 d, D-211), so the endpoint check passed
  while the SHAPE went unlooked-at** — [[feedback-a-summary-statistic-is-not-the-curve]], one panel over.
- **DO NOT record this as the cause of D-211 §9's day-1 miss.** Re-scoring pH matched on *extent*
  instead of the calendar **OVERSHOOTS**: at Tyrell's day-1 flux the model gives **5.195** against a
  measured **5.258-5.377**, flipping day 1 from **0.070 too alkaline to 0.063 too acidic**. Two
  causes are live — a slow ferment and a too-generous acid yield per gram — and they partly cancel
  on a time-matched read. Flux-matching is not clean either (it shifts growth-linked N uptake).
  **D-211 §6 reserved this same residual for an unbuilt buffer term; that reservation now has an
  unscored competitor. FLAGGED, not resolved — it needs its own beat.**

**What shipped: tests only, no `src/` or parameter change.** `TYRELL_ACID_COURSE_PPM` (Figs 9/10/14
interiors, four-strain means, ±2 ppm, anchored at both ends to numbers recorded beats earlier),
`TYRELL_FLUX_FRACTION` (Fig. 4's extract panel), two **strict xfails** in D-208's idiom (state what
is true of the source and false of the model, so a fix turns them GREEN — never pin the current
wrong value, which is what D-207 refused), and one anchoring test that passes today and must.
**Both xfails were re-verified under `--runxfail`** and that caught one failing on an
`AttributeError` in my own test — [[feedback-verify-an-xfail-fails-for-its-stated-reason]].
