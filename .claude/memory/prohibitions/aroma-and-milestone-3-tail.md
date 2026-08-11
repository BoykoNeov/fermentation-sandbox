---
name: aroma-and-milestone-3-tail
description: "Esters, fusels/2-PE, sulfides, closures and the Milestone-3 tail (D-117 to D-136, D-146, D-176)"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — aroma compounds and the Milestone-3 tail.** Detail split out of
`.claude/memory/project-fermentation-sandbox.md` at D-185; that file's ledger points here by
path. Read it when working on this subject. Every bullet is *what it forbids* + the D-record
to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue
past it from this file.** **Never evict an old prohibition to buy a line.**

**Milestone-3 tail (D-146) — two REFUSED, one BLOCKED; do not re-open as unbuilt**
- **Methionine sink + `methionol` is BLOCKED, not deferred.** Don't build it, don't pick a value:
  `f_non_ehrlich_methionine` is **ill-posed** (leaves its own `[0,1)` domain across one panel). The D-118 gate
  **passes** ⇒ the block is on the parameter; Crépin/Rollero **structurally incapable**, not silent. Unlock:
  **¹³C-methionine tracer**. **Sotolon's OAV is NOT enantiomer-split** — `threshold_sotolon_wine` is the
  **racemate's**, pool racemic; D-107's "soft by ~100×" is **withdrawn**.
- **Do not "finish" the `oav`→`magnitude` rename.** Scoped to `DescriptorReading.oav` alone —
  `OAVReading.oav`/`oav_series`/`oav_tier` are genuine OAVs and a test fails if swept. **`TOPIC_RULES`
  (`tools/gen_decisions_toc.py`) fixes a new record's bucket — edit the rules, not the heading.** Has `methion`/`sotolon`/`assertion`; a `magnitude` rule **misfiles D-53** and `surface` **misfiles D-128**.

**Sulfides (D-135) + closures (D-136)**
- Do NOT rebuild bottle reduction as thioacetate/disulfide precursors — D-101's mechanism guess was WRONG.
  **Release only** (white MeSH ~4× under-predicted: a limitation, not a bug); **no copper coupling** (a future one must be asymmetric); **temperature-FLAT** — never ship my ~158 kJ/mol.
- **Technical cork ships BELOW screwcap** (nominals only — D-162); µL→µg is the authors' own **1.43** (STP); seven
  "do not fix" traps in D-136. `ClosureOxygenIngress.reads` is **`()` by design**, pinned. **Its 10 band edges are
  the archive's ONLY two-sided provenance pins besides `ethyl_acetate_eq`** (D-162/D-163)
  [[feedback-pin-the-band-not-the-nominal]].

**Fusels / 2-PE (D-117 → D-120)**
- `f_non_ehrlich_phenylalanine` = 0.975, `f_de_novo_2_phenylethanol` = 0.9827. **Never put 0.963 in a sampled
  field** — different denominator [[feedback-rejected-values-must-be-unreachable]]. **A de-novo fraction is not scale-invariant.**
- **The blocker moved, it did not lift** — now "does the Phe flux scale with total 2-PE?"; the **"only a T4
  snapshot" defence is REFUTED** (flat). The shipped 2-PE cap is **INERT at the realistic dose** — a carbon-refund
  guard, **not** the 18.9% fix; never conflate. The live lever on isoamyl is `f_non_ehrlich_leucine`, not a route.
  Crépin's 0.815 vs Minebois's implied ~0.29 = open **D-103** conflict — **never averaged**. Rayne 2016 CALC ester k is 6–18× off R&O MEASURED.

**Esters (D-123 → D-127, D-176)**
- Three ester Processes ship, all SPECULATIVE. `k_H2T` NEGATIVE, shipped faithfully (wine-only). **No pH factor**
  on hexanoate/EtOAc — R&O T-VII ratios are isoamyl's. MCFA esters hydrolyse **WITH** the acetates; the *forming*
  family is branched/polyprotic ethyl esters, **none tracked**; **no uniform C4–C10 constant**. **D-131:**
  book-sweep note WRONG, primary **T6** calibrates. **D-130:** DHA omitted. **D-129:** WHERE/HOW-SHARP.
- **EtOAc eq is BERTHELOT-COUPLED (D-176, corrects D-127) — `ethyl_acetate_eq` is NOT the equilibrium**, it is
  its value at an anchor, ×`acetic_acid_typical`×ethanol. **`ethanol_ref_ester_eq` 94.68 = 12 % v/v E-rate matrix,
  NEVER the 14 % RATE matrix**; refs PINNED ⇒ moving one is refused **at load**. **Beer HYDROLYSES** ⇒ kills the
  negative `Byp` (pool **0, no producer**, D-16 **Flagged**) — but that is a **MARGIN, +0.741 mg/L at the JOINT
  corner of TWO bands** (flips at **+3.6 % ethanol**; 0 of 180 members reach it), **never an impossibility**;
  **wine +24 %**. Acetic band IARC ale **[12,155]
  PRINTED**, nominal **CONSTRUCTED**; **never Wang's 311 as nominal** (sour-inflated; its 10 % backs **FORM** only). Beer lands **below** published — **accepted**. **Beer pH 4.30 is a COINCIDENCE** (empty balance; 7.0 at t0 *and* packaging).
