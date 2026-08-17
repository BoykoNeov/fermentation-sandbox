---
name: beer-acid-base
description: "All six beer acid-base beats (D-178 to D-183) - the beat is complete and beer's pH is a prediction"
metadata: 
  node_type: memory
  type: project
  originSessionId: 09a91935-982e-429a-ba2b-c094a06612d5
  modified: 2026-08-17T08:32:47.404Z
---

**Live prohibitions — beer acid-base.** Detail split out of
`.claude/memory/project-fermentation-sandbox.md` at D-185; that file's ledger points here by
path. Read it when working on this subject. Every bullet is *what it forbids* + the D-record
to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue
past it from this file.** **Never evict an old prohibition to buy a line.**

**The FRAME is CLOSED and beer's pH validation HALVES (D-208). The term stays; the comparison moved.**
- **Tyrell's pH is a DECARBONATED reading** — MEBAK II 2.14 is *"pH (EBC)"*, EBC 9.35's scope is
  *"pH at 20 °C of DECARBONATED beer"*. **Never score `ph_of_state` against a published pH**: use
  `degassed_ph_of_state` (comparison ONLY — every pH-reading Process needs the in-vessel one).
  Headline **43.2-62.9 %** nominal, **8.3-82.7 %** joint, **0 corners reach**. D-182's rise was
  ~35 pp of FRAME, and the shipped `>0.70` floor needed the sample to keep **65 %** of saturation.
- **The frame is a BOUND, so walk `s` ∈ [0,1]**: days 4-7 unreachable at EVERY `s`, and the `s`
  fitting day 1 (0.150) leaves day 7 at 5.17. **Never re-anchor the floor to 43.2 %** — that pins
  `s`=0; the claim is an `xfail(strict)` on the day-7 LEVEL. **~0.4 pH of acidification is MISSING.**

**Beer pH was validated on ONE NUMBER and the CURVE read (D-207). Test-only; nothing moved.**
- **Never call the pH-drop test a TRAJECTORY test** — it scores ONE endpoint fraction (87.1 %, inside
  its own 77.6-97.0) while day 1 is **0.195 pH too acidic, 8.1× the read tolerance**. Model is BELOW
  the band early, ABOVE it days 4-7 (~2×, these do NOT carry it): it **overshoots then stalls**.
  **85.1 % of the day-1 drop is D-182's CO₂ term**, saturating by day 1 then flat 13 d — its docstring
  said "within hours"; unexamined ≠ hidden. Uptake confound RULED OUT (14.1 % vs measured 15 %).
- **Its day-1 SIGN is CORRECTED by D-208** — the fork it called unclosable closed CO₂-FREE, so day 1
  is 7.8× out the OTHER way (not acidic enough) and §5's CO₂ attribution is off-frame. **MEBAK = ZERO
  corpus hits was the wrong search, not a block** — the standard publishes its own scope line. Day 7
  is ABOVE the envelope either way. Days 1-6 stay **data no assert reads** (a shape pin encodes the
  defect); **nothing in the suite pins any pH SHAPE** (88 refs, 11 asserts, 2 walks).

**Beer acid-base, part 6 — acetic's RATE LAW (D-183). The SPIKE is still NOT modelled.**
- **NEVER re-propose the keto-acid excretion/re-assimilation pair for acetic** (D-180 §9's own
  proposal): the figure INTERIORS refute BOTH halves — 86 % of the rise in the first **15 %** of the
  flux, and **half the fall at ZERO flux**. Removal REFUSED on measurement: fit gaps all under the
  ±3 ppm read tol, floor **unidentifiable** (mean fits the bound at 0), endpoint would be
  **horizon-dependent** (108/65/20/**0** ppm at d7/14/30/400). Unlock = a 2nd wort's time course.
- Acetic LEFT `ORGANIC_ACID_SPECS` (now **3**) for growth-linked `AceticAcidOverflow`; modifier moved
  `for_uptake`→**`for_growth`** (mis-move is SILENT at 15 °C). Shipped curve **MONOTONE, asserted**;
  RMSE 61.6→32.5 only fixes *when*. **`Y_acetic_sugar_beer` was NEVER a measured yield** — retired,
  kept as the counterfactual [[feedback-a-derived-yield-encodes-its-rate-law]]. Joint band **10 dims,
  59049 corners**, in the SHIPPING commit. Headline **77.8-97.3 %**: +0.2 pp, PRE-REGISTERED and
  **headline-neutral by construction** (0.1233 pp/ppm) — never credit it. Endpoint scaling moved
  **gravity → YAN**; neither direction is measured.

**Beer acid-base, part 5 — dissolved CO2 (D-182). Both omitted terms BUILT; do not re-open.**
- **NEVER read the `CO2` slot as dissolved** — it is cumulative EVOLVED GAS (~40 g/L beer, ~100 wine, vs
  ~2 saturation); the balance reads **`min(evolved, C_sat(T))`**, unsmoothed. **No state slot, no registry
  entry** (`_totals_molar` would skip it SILENTLY) — own positional arg beside `byp`, **REQUIRED not
  defaulted**. Anchor is CO2-free by construction ⇒ D-178's phosphate absorption does NOT apply.
- **MEDIUM-AGNOSTIC and measured, not preferred**: wine **0.0007 pH** (sits 3 units BELOW the apparent
  pKa 6.43 — D-178's own geometry inverted), 400-d aged wine **≤1e-5**. **TA EXCLUDES it** (degassed
  sample; bitwise-pinned). **NEVER gate it** [[feedback-a-gate-is-a-discontinuity-the-solver-probes]].
- Headline **77.6-97.0 %** nominal / **63.8-109.4 %** joint (**9** dims, 19683 corners) — a corner reaches
  again and that was PRE-REGISTERED. **The two terms are NOT additive**: D-181's +0.2094 is **+0.1128**
  beside CO2, because carbonic dissociates MORE as pH rises. `pKa_carbonic_1` needed a HAND entry in
  `PKA_PARAM_NAMES` (not a registry member ⇒ derivation skipped it). **1 atm only, no vessel pressure.**

**Beer acid-base, part 4 — the SINK (D-181). The headline now agrees WORSE, on purpose.**
- **NEVER flux-link `WortAcidRemoval`** — Tyrell Table 2 scores all three `--` at LOWER Krausen and `0` after,
  and the house idiom peaks MID-ferment. First-order to a **measured floor**, temperature-flat. **NO mechanism
  asserted** ⇒ the 3 slots are **OFF every ledger** (`iso_alpha` idiom) and **absent from `MOLAR_MASS`** so a
  future producer RAISES; never weight them "for completeness". **Corrects D-179**: buffering index ranks
  these last, **charge REMOVED** ranks them with the produced four — a fully dissociated acid loses the MOST.
- Nominal **42.7-62.2 %** (was 63-92), joint **7.6-82.2** (was 41-105): **nothing in the band reaches** now.
  Denominator is **0.81** (extreme-strain mean) in code, **0.8125** in D-180 prose — **never mix them**.
  **All 4 new pKas are consequence-free** — `pKa_oxalic_2` was shipped as "the sensitive one" and is 0.027 pp:
  the claim was about a REAL beer (4.78-4.90), the model stops at **5.24**. Re-measure when CO2 lands.
  Missing base **+0.2094 pH** — **but only CO2-FREE; it is +0.1128 since D-182**, never mix them.
  Capacity re-anchored **1.6708→1.5481** (5-acid back-solve reproduces D-180
  **bitwise**). ONE shared `k` (3x/⅓x moves pH **0.00011**). `pyruvic` ≠ wine's `pyruvate` — **reconcile before
  beer ever wires keto-acids**. **Dissolved CO2 BUILT at D-182** — nothing of opposite sign left [[feedback-build-the-term-that-makes-agreement-worse-first]].

**Beer acid-base, part 3 — the PRODUCER (D-180). Beat COMPLETE; beer's pH is now a PREDICTION.**
- **Tyrell 2013 is TWO datasets — D-179 read only Table 1.** Its own EBC trials (Figs 4, 6-14) give ONE wort's
  acid course d0-d7 for 4 strains **plus the pH+extract of the SAME ferments** ⇒ wort seeds, yields and the
  divisor. **FIGURE READS, never table-grade.** Seeds are `*_typical_wort` now: a producer on D-179's
  finished-beer seeds lands **pH 4.26** ⇒ seeds+producer are **ONE decision**. Citrate gets **NO yield** (sourced 3×).
- **The pH-drop agreement is a MARGIN, not validation — and its RANGE IS SCOPED.** **63-92 % only at NOMINAL
  yields**; over the **JOINT** band (4 sampled yields × pKa) it is **41-105 %**, so "must fall short" is TRUE at
  nominal and **FALSE band-wide** — the archive's point-vs-band shape, **5th instance** (amendment). Upper-bound
  guard is scoped to nominal; joint corners pinned with **NO** upper bound. Pre-reg's 3 arms: adding the 3 acids
  that FALL (pyruvic/formic/oxalic, **no state slots**) makes it **WORSE** (32-70 %); +dissolved CO2 gives
  76-104 % ⇒ two omitted terms of OPPOSITE sign. **Never read `CO2` as carbonic acid: cumulative EVOLVED GAS.**
  Band arm must **RE-ANCHOR the cation per member** (fixed nominal cation reports 72-80 % *and* moves the start pH).
- **Beer `Y_byproduct_sugar`=0 is now LOAD-BEARING** — succinic would double-count (own slot + `Byp`). Yield rides
  the **shared** `fermentative_uptake_rates` (bitwise-pinned) and **must** be named in `for_uptake`'s targets (D-32
  coupling) else the yield breaks with temperature. Booked on the FERMENTATIVE flux = **97.12 %** of ΔS (growth
  takes the rest) ⇒ ~2.9 % under the anchor, **stated not tuned**. Peptide capacity re-anchored on WORT acids
  (1.5125→**1.6708**), closing mismatch #3 and **worsening #1** (now a traverse); costs 0.020 pH.

**Beer acid-base, part 2 — the STATE (D-179). Beat COMPLETE; do not re-propose any of it as unbuilt.**
- **NEVER merge `WINE_ACIDS`/`BEER_ACIDS` into one registry.** Both media carry a **`citrate` slot** and only
  beer's is charge-active (wine's is carbon-only, **D-31 Flagged not fixed**) ⇒ a union silently changes wine.
  Keyed by **`StateSchema.medium`, an explicit LABEL — never sniff slots** (that is D-178's bug). Wine resolves
  to the identical dict: pKa values **bitwise** (`float.hex()`), 94 slots, `ph_tier` PLAUSIBLE.
- **`ph_tier`/`molecular_so2_tier` take a schema and MUST be passed one**; unscoped = **conservative union**
  (SPECULATIVE), never wine — beer's peptide buffer is speculative and would drag wine down. **`PH_SYSTEM_READS`
  is the union ON PURPOSE**: it shifts a wine ensemble's draw SEQUENCE (nominals/bands unchanged), and the
  alternative is D-160's silent band-narrowing. **Beer's acids are INERT** — composition, not fate (D-16 open).
- **The EtOAc gate is a POPULATED BALANCE (`charge_balance_is_populated`), never slot presence and never
  `cation>0`.** Slot presence gave beers an **empty balance ⇒ pH 7.0 ⇒ 5000×**; but `cation>0` is wrong the
  OTHER way — an acid-dosed unanchored WINE is pH **2.23**, a real `h≈11.8`, so that gate changes wine MORE
  than the beer bug it fixes. **0 suite scenarios hit it ⇒ either error ships green.**
  `initial_ph` is beer's **opt-in gate** (absent ⇒ all acids 0, byte-for-byte pre-D-179). EtOAc **lingers 2.03×**
  (6.73→13.69 mg/L/400 d) — **rate 12.6× ≠ outcome 2.03×**, the pool relaxes to a floor. pH drifts **−0.0041 via Byp**.
- **`peptide_buffer_capacity_beer` is PINNED zero-width ON PURPOSE** — it is a **function of the pKa** via the
  back-solve, so banding both samples pairs reproducing **no measurement**. Band rides **`pKa_peptide_buffer`
  [3.86,4.50]**, which IS sampled (11.96-14.61 mg/L; pinned control 0.000000). BC **0.309-0.544** with the peptide
  term vs **0.103** without. **A compile-time dose is in no `reads` ⇒ never sampled** — check before calling it a band.

**Beer acid-base (D-178) — the SOLVER half; the PREMISE was wrong**
- **NEVER re-propose malt phosphate as beer's buffer.** pKas **2.15/7.20** vs beer ~4.3 ⇒ charge **0.9867→0.9990
  FLAT** over 4.0-4.6; 0.15 mEq/L/pH even at 700 mg/L vs organic acids' ~2.3 — and an **inverse anchor absorbs a
  constant charge**, so it is a **near no-op, not a weak buffer**. Source questions it too. **Citrate** (triprotic,
  2 pKas in range, source-NAMED) is why the n-protic branch exists. Pinned `test_acidbase_polyprotic.py` (17).
- **Brewing's BC is a DIMENSIONLESS LOG RATIO, not mEq/L/pH** — that UNIT, not the paywall, is why **9 hosts
  failed**. Open host = **Peyer 2017 UCC thesis, CORA `10468/4694`** (= the 402'd jib.447). Beer today **0.001**,
  organic acids alone **~0.20**, published wort **1.18** ⇒ **peptides are the majority**; free amino acids ~10 %.
- **n-protic dispatches at `len>=3` ONLY and is NOT bitwise BY DESIGN** — never route 1/2-pKa through it, never
  delete the fast paths as duplication. **`EsterHydrolysis`'s gate is `"tartaric" in schema`** — the old
  `cation_charge` "wine marker" would **CRASH** beer once beer has a cation. **EtOAc `pH_ref`=3.3, beer is ABOVE
  it ⇒ 5-20× SLOWER, not faster** (my pre-reg had the sign inverted). **Succinic 36-166 vs ~900-3500 — never
  average.** D-176's acetic [12,155]: IARC calls it ALE, Coote/Kirsop calls the same range LAGER — **transposed**.
