---
name: aroma-and-milestone-3-tail
description: "Esters, fusels/2-PE, sulfides, closures, the fusel NODE (D-245) and the Milestone-3 tail (D-117 to D-136, D-146, D-176)"
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

**D-259 — the growth-anchored sink, RE-MEASURED. Read before citing D-104's refusal or D-257 §7.**
- **D-104's "leucine 20.9 %" is ONE EDGE of an input recorded NOWHERE.** No yeast protein
  composition exists in `src/`, in the D-104 record, or in any receipts folder (`must_aa_fraction_*`
  is the MUST spectrum, a different thing). Re-measured across a stated bracket leucine spans
  **13.1-22.0 %** and D-104's number sits at the **top**. The refusal survives at every edge — the
  point number does not. **Never quote 20.9 % as "the model's" growth-anchored split.**
- **"Exactly reversed" is CORRECTED.** Only the two ENDS are inverted now; **isoleucine > valine
  matches Crépin**. Valine moved 45.8 → ~25 % and it is the only precursor that gained a second
  Ehrlich branch after D-104 (D-111).
- **THE POSITIVE RESULT: growth-anchoring LANDS the sourced fate where a de-novo route exists.**
  Phenylalanine reads **97.2-98.5 %** against its sourced 0.975. **Verified by mutation** — zero the
  de-novo share and it falls to **37.1 %**, so it is the ROUTE, not supply limitation (leucine's
  ~17 % is the control). **The inversion is a property of the four precursors WITHOUT a de-novo
  route, not of anchoring to growth** — so do not cite D-104 as refusing the FORM.
- **THE FENCE, AND IT IS SCOPED — do not quote it wider.** Un-inverting leucine **by cutting the
  Ehrlich draw** needs 11.9-40.7×, which takes the model's leucine share of isoamyl from **1.51 %**
  (already BELOW Rollero's measured 3.4-8.2 %) to **0.037-0.127 %**. Refused. It does **NOT** refuse
  un-inversion in general and does **NOT** fence D-116: the split is a RATIO, and the
  **numerator-side repair — raising growth's own draw — is UNTESTED**. The real defect is the
  denominator: the model eats **173 µM** leucine against a protein demand of **580-1088 µM**.
- **D-257 §7's blocker DOES clear** (phenylalanine 45.45 → 0.09 % under the `k` rescale) — but that
  is **within-fixture**; never quote it against D-257's own 20.3/65.8 %. **The build stays owner-gated.**

**D-260 — the NUMERATOR side. Read before proposing ANY growth-driven precursor work.**
- **D-259 §5's "that side is untested" is ANSWERED and the answer is a REFUSAL — do not re-open
  it as open.** The split and Rollero's leucine tracer are **ONE knob**: the pool ends 0.00 %
  in every arm and isoamyl moves < 0.01 %, so `tracer = (1−split)·consumed/isoamyl` binds them.
  Shipped = Crépin end (81.50 % / 1.507 %); growth-anchored = Rollero end (27.58 % / 5.900 %).
  **Neither is the other's repair.** The lever WORKS (λ=5 → 73.29 %) and pays 2.176 % on the
  tracer. **Never say growth-anchoring "misses the tracer"** — it LANDS it; the shipped form misses it.
- **The tracer gain is NOT growth-anchoring's** — the shipped form at a matched `f`=0.174 reads
  6.730 % against the counterfactual's 6.726 %. It belongs to the SPLIT, whatever produces it.
- **D-259's bracket was the PRE-MODIFIER frame — corrected to 21.3-33.7 % (leucine).** A
  growth-anchored draw carries **none** of growth's rate modifiers (they attach by Process NAME;
  D-32 is the rule). **D-104's 20.9 % now sits BELOW the corrected low edge, not at the top.**
  Refusal + order correction SURVIVE; every number moves, and "monotonically inverted" does not
  describe the corrected form (valine brackets Crépin's 41, threonine sits above her 38).
- **The live question is the isoamyl DENOMINATOR, and it is commensurability, NOT calibration.**
  Joint satisfaction needs ≤ 1170 µM; the fixture makes 2123; Rollero's own print 793-1365.
  **`k_isoamyl_alcohol` is RIGHT (D-112) — this is not a licence to re-fit it.**

**Fusel node (D-245) — read before touching the de-novo helpers or citing D-120**
- **`_de_novo_share` must use `ehrlich_primary_share`, never `(1-f)`.** Valine's Ehrlich carbon
  splits AGAIN (0.15 isobutanol / 0.23 isoamyl via KIC), so `(1-f)=0.38` charges isobutanol with
  isoamyl's carbon — **2.533x, live D-111→D-245**, and it read as a model defect the moment D-244
  pushed it under the floor. **ONE helper, four callers**; the sibling had the rule right in a
  comment the whole time. Its result is a **BOUND**: the `ehrlich_draws` headroom cap binds for
  the first ~1.35 h (≤12.8 % of the valine pool), so isobutanol is ≥89.6 %, point estimate 90.5 %.
  **Never "tighten" it by integrating the branch** — that is D-103's quadrature.
- **D-120's no-cap refusal has LOST BOTH measured legs — but do NOT build the cap here.** Direction
  flipped (isoamyl 5.42 vs 5.34 %, isobutanol 9.47 vs 8.78 %) and the instrument now bites (Phe
  stops exhausting, 12.8 % left; the cap moves 2-PE's realised share 12.7 %). **The isoamyl trip is
  inside the harness's own cap-window systematic (~5.30 %, i.e. UNDER) — never build on it**;
  isobutanol's survives. A sourced `f_de_novo_isoamyl` does not exist here and the 2-PE closure
  algebra re-run on the model's own abundances is **refused at D-206**. Sourcing beat.
- **α-KB: the claim survives, the margin does not.** Propanol demand / α-KB throughput 2.60 → 1.358,
  re-pinned two-sided [1.30, 1.42] with `> 1.0` asserted separately as the claim. **Do not restore
  a one-sided floor** — it cannot catch the excretion flux growing toward the demand. The supply
  argument now rests on 36 % headroom over an AUTHOR-ESTIMATED rate.
- **Never mark a multi-assert test xfail whole, and never mark a LOOP over a registry whole** —
  D-244 did it twice and D-245's own first commit did it a third time: the Minebois guard looped
  2-PE → isoamyl → isobutanol, died on isoamyl (the leg the record calls NOISE) and never evaluated
  **isobutanol, the leg it calls the real flip**. Split per member; the **12.7 % cap bite is a GREEN
  pin** now ([10 %, 16 %], the receipt for `Flags: D-120`). Five xfails, not four
  [[feedback-an-xfail-buries-the-asserts-after-it]]. **D-112 finding 4's mechanism is
  CORRECTED** (the evaluation point, not "the dose's deamination-N sustaining the gate"); its
  conclusion — `k_isoamyl_alcohol` is right — **stands**, and the undosed anchor test is the control
  that proves it (bit-for-bit through D-244).

**D-266 — the JOINT repair (D-257's blend + D-259/D-260's growth sink), MEASURED. Read before citing either refusal.**
- **Neither refusal was a measurement of the pair, and the pair was never run until D-266.** It clears
  D-257's blocker (Phe 72.9 % → 0.000 % left) and keeps the blend's timing (36.6 / 48.5 / 57.0 % at NT vs
  42-54; response 1.30×). **The pre-registered "pays on the tracer" was WRONG**: leucine split 27.6 → 54.5 %
  (46-61 across the composition bracket), tracer 3.62 % INSIDE Rollero's 3.4-8.2. Isoamyl on Crépin's must
  2123 → 2177 µM, so **D-260's line does not move — the joint arm is a different point on it**, ceiling 1170 unchanged.
- **What the pair costs, and no half could show it:** ile/val/thr over-shoot Crépin (66/65/81 vs 51/41/38);
  the growth sink ALONE breaks the propanol 0.80 floor (0.714) and the blend restores it (0.875); Rollero-must
  leucine enrichment brackets (1.78/5.44/8.33 vs 3.4/4.2-4.7/6.8-8.2) while the labelled AMOUNT is 1.6-2.5× over;
  valine's tracer moves the WRONG way. **Per-axis closer, not overall — the fork stays the owner's.** Receipts
  are IN the repo: `docs/receipts/d266-joint-fusel-repair/`; guards in `tests/test_fusel_joint_repair.py`.
