---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
  modified: 2026-07-28T11:00:41.661Z
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation engine in Python (uv, scipy/numpy/pydantic). Repo: https://github.com/BoykoNeov/fermentation-sandbox (branch `main`).

**Session-boot context: PROHIBITIONS and POINTERS only** — not a changelog; reasoning lives in the
D-record. A bullet explaining *what happened* is too long: cut it to what it forbids. **Cap 250**
(`.claude/hooks/check_memory_size.py`; [[feedback-batch-end-ritual]]).

## Where the records are

- `docs/DECISIONS.md` — canonical archive, **~16.7k lines: never read it linearly.** Its generated top block
  gives a subsystem cut, the ordered list, and a **correction map** (⚠). **The ⚠ lives only in the index** — a
  record carries no marker for corrections *against* it, so check its index row first. Append per `CLAUDE.md`'s
  `Corrects:`/`Flags:` convention, then run `tools/gen_decisions_toc.py`. **File is LF.**
- `docs/ARCHITECTURE.md` (seams); `docs/plans/milestone-*.md` (its "Active beat: sensory" header is **STALE** —
  that beat closed at D-95/D-98); `CLAUDE.md` (prime directives + archive conventions).

## Status (2026-07-28)

M0/M1/M2 **complete**. **Milestone 3** (sensory/OAV + Tier-3 aging, owner's pick at D-66) in progress, at
**D-152**. Aging build order **built** — `aging.py` carries 24 Processes; sensory 1a/1b closed. **D-139's
leftovers ALL closed** (§2.4 D-148, §2.5 D-149). Suite **1443 passed**. Wine schema **94 slots** / beer
**47** — `quinone` in both regardless of set. **Three** oxidative sets (`direct` default / `cascade` /
`direct_burst`). Most remaining Milestone-3 work is **blocked on external sourcing**, not on building.

## Do NOT re-propose — I did, twice, from stale "Next:" breadcrumbs
[[feedback-verify-latest-state-not-breadcrumbs]].
- **All lumps are speciated**: esters→3 (D-96), fusels→5 (D-99), amino_acids→8 (D-100), mercaptans = a
  methanethiol false-lump (D-110); `lumped` kept **dormant**.
- **Beat 1b (descriptor projection) is COMPLETE** (D-95 + D-98); only masking remains, `cosα`-blocked.
- **Shipped and spent, do not re-open as unbuilt:** D-128, D-129 (`EthanolToleranceDeath`), D-130, D-131,
  D-132/133/134 (O2 rate), D-135, D-136.
- **Isoamyl de-novo entry — REFUSED at D-120, measured not built.** A one-directional ceiling can't fix an
  under-attribution; every alcohol is on the wrong side; it can't reach isoamyl's valine branch; it's a rate
  knob on a supply-limited quantity. [[feedback-measure-which-side-before-building]]
- **Closed:** leucine shortfall (D-112); shared-BAT parsimony (D-116); Rollero (D-115); ester-aging "deferred
  half" (D-121). **Beer 3-sugar kinetics are NOT in the 5 beer books** — don't re-sweep.

## Live prohibitions, by axis
**Oxidation (D-132 → D-137)**
- D-132's phenolic boost is **additive, never proportional**, browning-side only.
- D-133's `burst_antioxidant` is **EXCESS** over D-132's steady rate and must **read none of**
  `tannin`/`anthocyanin`/`so2_total` — the anti-double-count, and it binds. (*Medium* wiring is a separate
  question, settled at D-147 below.)
- `k_copper_multiplier` = **600** L/g. **§2.5 CLOSED (D-149)** — never re-open as "the printed table says copper
  is stronger": Nguyen T3.1's 3.93× converts to **2092 L/g**, within 5% of the rejected 2000 and failing the
  same budget (2.32×). 600 is held by the **source** (free CuSO₄ model wine; initial- vs steady-rate), **NOT
  Ferreira's ceiling** — at 2000 direct reads 1.72×, inside; it binds only via the cascade (k=1927), which
  post-dates D-134. Cu share of O₂ uptake **61%/21% direct · 28%/15% burst · 100% cascade** ⇒ isolated `f_Cu`
  **is** the cascade's number, an upper bound elsewhere; whole band clears. Re-open needs **real wine, ≥3 Cu
  levels**. Guard: §Guard 4. **BOUNDED now, not null — D-152 converted it; never re-open as "only a null".**
  The copper-**orthogonal** L16 bounds **k ≤ 918 L/g** (1003→663 over `copper_typical`'s band) — pH-**3.3
  simple effect**, **TWO-SIDED** 95%. Excludes **2092, 2000 AND band-high 1500 under ALL SIX arms**. **Never
  rebuild it on Table 2's SDs** — replicate noise is **1.07%** of D-151's residual, SE understates **9.67×**,
  reads `k ≤ 0`; under that term Cu's F = 269 vs crit 3.89 though the paper calls Cu *not* significant.
  **Keep condition 12** — dropping it TIGHTENS 1.754→1.417, **opposite** to D-151's call. Transfers because
  the frame **matches the source** (Danilewicz free-CuSO₄ model wine ↔ dosed-CuCl₂ reconstitution). Guard 7;
  resolution **~1148**, does NOT separate 1500 from 2000. **NEVER cite it as evidence FOR 600** — only **1 of
  6** arms admits 600 (the shipped, loosest); both *marginal* arms read **58/133, under band-low 200**
  (amendment 1). One-directional: against higher k, never for 600. **LIVE DEFECT:** 1500 excluded yet
  `ensemble.py` draws `triangular(200,600,1500)` ⇒ **~29% of draws excluded**
  [[feedback-rejected-values-must-be-unreachable]]. **Next beat = the BAND, not the value**; sweep every
  sampled band. **No Fe(III)
  state** — NOT D-134's "iron in surplus": QSS rests on a **~18×** separation (Nguyen T3.1 3.3–63×).
- **A pH term on activation/`k_browning_eff` is REFUSED (D-150) — never re-open as "Nguyen's table shows strong
  pH dependence"; it was converted.** (a) **No mechanism** — phenolate (pKa 9–10) predicts **10×/unit** vs
  3.867× measured (needs 68% pH-flat in a catechol-only wine). (b) **NOT separable from copper**: Cu ratio
  5.200/3.929/2.586 across one pH unit (**2.011× swing**; D-10 needs 1.000×) ⇒ per row `k_copper_multiplier` =
  2425/**2092**/1511, so `f_pH` **re-opens D-149**. **(b) is the STRONGEST leg now (D-151), never "an artefact
  of Nguyen's non-orthogonal dosing"** — Carrasco-Quiroz 2022's copper-**ORTHOGONAL** L16 (pH×Cu unaliased)
  reproduces it: swing **1.826×**, F=5.95, same direction, **stronger** without the outlier. Its `ΔO₉₀₋₁₀/Δt`
  null is **NOT a partial unblock** — r=0.988 with `1/Δt`, a **duration statistic in rate units**. (a) is
  **worse**: four exponents, no mechanism. (c) **Wrong statistic** — Carrascón 2018 (8 real reds) splits
  **initial 4.585×/unit vs steady 1.619×**; Nguyen's 3.867× is *initial* on a *steady* node ⇒ **the pH term's
  home is the burst (D-133)**, D-142/D-149 both pointed wrong. (d) 0.2093 dec/pH is affordable but
  **between-wine**, pH⟂variety collinear (r=0.608, n=8, ns); a sim term is **within-wine**; model
  **1.0000×/1.0060×** over pH 3.26–3.61 vs real 1.158×. Guard 5's bound is an **OBSERVED spread (1.420×), NOT a
  ceiling**. Re-open needs a **within-wine** real-wine steady-OCR series, ≥3 levels, **now also measuring Cu**.
  **`initial_ph` anchors t=0 only** (→3.2084/3.5135 at the dose): **no way to SET an aging pH.**
- First-order in `[o2]` RE-CONFIRMED; **no MM/`Km`** (Ferreira's "linear" headline is a cross-cycle
  re-saturation artifact). **Do NOT re-attempt the Ferreira/Carrascón PLS extraction** — all three blocked
  (`_findings/D-134-copper-ocr-sourcing.md`).
- **The O2 gate is Fe(II)+O2** — polyphenols and sulfite react with O2 negligibly. SO2's ~90% budget share is
  right in *size*, **INVERTED in role** (enabler, not competitor) — why D-72's wrong mechanism yields a
  right-looking 1:2. Copper re-homes to the O2 step.
- **Never re-fix the acetaldehyde/phenolic inversion inside the parallel frame — it provably cancels.** Approve
  the cascade on MECHANISM only; it will **not** flip the sign (total O2 is phenol-INDEPENDENT under D-136).
  **Never tune to phenol ⇒ more MeCHO.**
- Danilewicz 2011 is **PARAPHRASE** — pull the abstract before it backs a `source:`. Its one in-wine branch
  measurement is **mixing-limited**: bounds the Fenton share, calibrates nothing (its no-SO2 constants are
  fine). The atmospheric "H2O2 + S(IV) is pH-independent above pH 2" does **NOT** transfer (in wine free SO2 is
  a fixed pool, so rate ∝ [H⁺], ~10×/pH unit).

**Burst — WIRED and NON-DEFAULT (D-147). Not dead code, not unbuilt, not refused.**
- `oxidative="direct_burst"` = direct **+** the burst; a **superset of direct**, not a mechanism. **No
  `cascade_burst`** (D-138's node stays UNDETERMINED; its "transient modifier" is FLAGGED).
- **Never make it default**: moves all 31 D-140 pins (`o2` −37.4%, `A420` −36.4% @2y, floor ~1e-8) and re-opens
  the Danilewicz direct arm — **moot under this placement, do not re-derive**.
- **The split is UNPINNED and stays so.** C1 (day-1 excess) HOLDS (0.93–0.95 vs 1.0) ⇒ pins only the **product**
  `k_burst_oxidation · burst_antioxidant_initial`. C2 (95% spent in 10 d) FAILS structurally: needs **3.3 mg/L
  O₂ through that route alone**, wins ~35% of an 8 mg/L charge ⇒ plateaus at 39.2% left. **Never re-solve the
  split** — a second unmeasured number (D-146's trap). Unlock: **Ferreira 2015 per-cycle O₂ curves**, not the
  averages.
- **Self-exhaustion is Ferreira's PROTOCOL, not the sink** (~8 mg/L/10-d cycle vs a cork's 2.09 mg/L/2 y,
  ~1400×) ⇒ under a cork a **permanent ~37% tax that GROWS to a plateau**. Two guards forbid calling it a
  transient — **never relax them to a "pool depletes" check**.
- **Seed follows the consumer**: sourced under `direct_burst`, **0.0** elsewhere; dosing `burst_antioxidant_gpl`
  with no consumer **raises** (D-45's "absent ≠ 0" **inverts** without one).
- Ferreira's **2.7× does not discriminate** — flips with an unfixed window AND with SO₂; **with SO₂ the base
  model already meets it unwired**. Burst narrows SO₂-dependence 1.56×→1.33× (its one point in favour). 2.7× is
  D-133's **paraphrase** — pull it before it backs a `source:`.

**Cascade — BUILT and NON-DEFAULT (D-141). Do not re-build it; do not flip it silently.**
- `core/kinetics/oxidative_cascade.py`, 8 Processes, `quinone` the one new slot, OFF the ledger (`h2o2` QSS:
  mean lifetime 28.5 s = `1/k` NOT `ln2/k`). **REPLACEMENT** via `_OXIDATIVE_SETS`;
  `get_medium`/`compile_scenario(..., oxidative=)`. **`"direct"` is default and stays default.**
- **Do NOT close the fate gap by tuning.** Fates move up to **25×** because re-homing `k·[o2]·[S]` →
  `k·[quinone]·[S]` rescales by `[quinone]_ss/[o2]_ss` (~0.06). Both settling sources are **paywalled** (nine
  hosts checked) — the abstract K values are pseudo-first-order at unstated concentrations, **never a `source:`
  field** ([[feedback-conceded-caveats-are-not-coverage]]).
- **Never re-derive `k_activation_floor` as a fit** — it **is** `k_ethanol_oxidation + k_browning_base`, so
  activation re-expresses **that PAIR only**, at `copper_typical`; it does **NOT** reproduce direct's TOTAL. The
  O2 budget agreeing is **supply limitation (D-136), NOT a rate-law check.** Activation **reads the reductant
  pools** — never lump.
- The 31 D-140 guards stand — **never re-derive their pins**; edit only the two seams. Beer's O2 is a **floor
  `>=5 mg/L`**, never 5.71. Pin rtol **1e-4** (BDF runs 1e-6); `quinone == 0.0` under direct is **exact**.
- **Benchmark EXISTS** (`tests/benchmarks/test_validation_danilewicz_so2_o2.py`, D-142+amend). **Never pin 1.7**
  — one dataset's mode, ABOVE the other's whole range; assert the limits 1/2 (exact asymptotes, emergent) +
  bands. The falsifier is the **traverse**: direct 0.0003 vs cascade 0.859 (D-141's "structurally cannot
  produce" was wrong — direct lands in-envelope).
- **Operating point is load-bearing — enforce the >10 free-SO₂ floor, never state it in prose.** At dose 40 free
  SO₂ ended 2.8/6.8 (invalid) and the headline INVERTED; the "4-7× quinone shortfall" + `xfail` are
  **WITHDRAWN**. The floor is `SIM_CURVATURE_FLOOR_MGL`, **NOT Miao's criterion** — **never re-encode his** (his
  *starting*-level cut SELECTS 40 and rejects 80; his r² can't discriminate, 0.98975–0.99987). **A borrowed
  criterion can be inapplicable, not just strict.** Keep the value + excluded-wines table; never restore the
  Miao name. **Never report a single-dose verdict on this ratio** — only direct's is unconditional (**1.7707**);
  cascade **straddles** (1.0704 @60 out, 1.1338 @80).
- **"The sim under-binds SO₂" is WITHDRAWN (D-143) — never re-assert.** Miao's buffering is an **ADDITION-method
  secant**, D-142 compared an **oxidation-path slope** — **never compare those two**, and never
  `oxidation_path_slope` vs `MIAO_BUFFERING_BAND`. Like-for-like the 4-carbonyl budget reads **1.32, inside
  1.2526–1.9882**; the positive claim is speculative (**never re-run a pool sweep**). **Not** a schema gap,
  **no** K≫h binder, **no species named**. Acetaldehyde is **0.0000 mM at D-142's operating point** (SO₂ dosed
  once, day 19) — a scenario artefact; a must dose gives 0.29–1.00 mM, and **D-145 refused it**.
- **The four D-142 artefacts are FIXED (D-144)** — Miao's **Table 4 intercepts ARE free SO₂ at O₂ exhaustion**:
  **read the intercept, never recompute from an assumed dose** (that flipped wine 5); D-144's replacement test
  **NAMES NO CAUSE — keep it that way**. Pins: pH does **not** move with total SO₂ (frozen-β 0.0005%); the locus
  guard is **pure algebra** (1e-8) — a state route re-entangles `ph_of_state`. Masses from `core.chemistry`
  (**M_SO2 = 64.058**); `_SO2_BINDERS` reads all **four**. Secant 1.3197 invariance is the **schedule's**, not
  the statistic's (D-145).
- **A green suite proves nothing about a quoted decimal — every assertion is a band.** Re-measure quoted digits
  after any units change (the mass fix moved the straddle). Quinone branching **still NOT settled** — the ratio
  comparison was never sensitive enough.

**Dosing schedule — REFUSED (D-145). Never re-propose must-dosing; benchmark unchanged.**
- **D-143's "five for five" is three for five.** **Import statistics from the shipped benchmark module** —
  re-derived helpers carried both errors.
- **The reservoir CANNOT move the Table 3 secant**: acetaldehyde K=1.5e-6 vs h=6.1e-5–7.9e-4 ⇒ **h/K 41–528,
  saturated both ends** ⇒ moves the **INTERCEPT**, not the slope (**weak-binder** statistic). Grid interior
  17/58, **none satisfies all six**; `free0*` 45.5 direct / 30.5 cascade is **floor-CONTINGENT — never state
  bare**; floor-robust depletion is **not** ratio-independent.
- **Depletion gap RETIRED, not reframed** — D-143 omitted `M_SO2/M_O2`; Table 2's slope is **mass-basis**, Eq. 1
  closes to 0.6492–0.9696, only the secant is invariant and it is in band.
- Four traps that each produce a *plausible* green: `ProcessSet.disable()` is silently undone by `begin_aging`;
  `param_values` is a **property returning a fresh dict**; mol/L ×1000 is mmol/L not mg/L; **`trajectory.y` is
  `(n_states, n_times)` — `y[-1]` is the LAST SLOT's series** (`quinone`, all-zero), not the final state, so
  `array_equal(a.y[-1], b.y[-1])` compares zeros (D-147).
- D-139's L9 over-draw **never fired**; **§2.4 CLOSED (D-148) — never re-open as "the Brett/quench over-draw":
  no quench draw exists** (D-141 D2). Live pair was Brett/POF; **summing is NOT the mechanism** — one drawer
  over-draws alone at K→0, in-band two ≈ one (−5e-11 vs the 1e-9 floor). **Never build the shared depletion
  gate**: `depletion_gate`'s 8 sites buy per-draw first-orderness, which `_decarboxylation_branch` has.
  Out-of-band numbers are **BDF step artefacts** (LSODA/Radau 5 orders apart, non-monotone in K) — never quote
  as calibrated. Live but unreachable: an undershot pool **freezes negative**, `conservation.py` weights both
  unclamped. Cross-domain trap: aging Process **names** are enumerated by `tests/test_fusel_*.py`.

**Milestone-3 tail (D-146) — two REFUSED, one BLOCKED; do not re-open as unbuilt**
- **Methionine sink + `methionol` is BLOCKED, not deferred.** Don't build it, don't pick a value:
  `f_non_ehrlich_methionine` is **ill-posed** (**+0.421 → −1.355** across ONE 22-strain panel — it leaves its
  own `[0,1)` domain). The D-118 gate **passes** (20–39%), so the block is on the parameter. Crépin never
  **detected** methionine; Rollero fed U-¹³C leu/val — both **structurally incapable**, not silent. The de-novo
  rescue needs a second guess. Unlock: ¹³C-methionine tracer.
- **Sotolon's OAV is NOT enantiomer-split.** `threshold_sotolon_wine` is the **racemate's** and the pool is
  racemic (non-enzymatic aldol) ⇒ commensurate. D-107's "soft by ~100×" is **withdrawn**; the 0.8/89 pair is
  model-solution, excluded by the file's own matrix rule. Cost is a scope limit and **asymmetric**: ≤2× under,
  ~56× over. [[feedback-rejected-values-must-be-unreachable]].
- **Do not "finish" the `oav`→`magnitude` rename.** Scoped to `DescriptorReading.oav` alone;
  `OAVReading.oav`/`oav_series`/`oav_tier` are genuine OAVs and a test fails if swept. `TOPIC_RULES` gained
  `methion`/`sotolon`; a `magnitude` rule misfiles D-53.

**Sulfides (D-135) + closures (D-136)**
- Do NOT rebuild bottle reduction as thioacetate/disulfide precursors — D-101's mechanism guess was WRONG.
  **Release only** (white MeSH ~4× under-predicted: a limitation, not a bug); **no copper coupling** (a future
  one must be asymmetric); **temperature-FLAT** — never ship my ~158 kJ/mol.
- **Technical cork ships BELOW screwcap.** µL→µg factor is the authors' own **1.43** (STP). Seven "do not fix"
  traps are listed in D-136.

**Fusels / 2-PE (D-117 → D-120)**
- `f_non_ehrlich_phenylalanine` = 0.975, `f_de_novo_2_phenylethanol` = 0.9827. **Never put 0.963 in a sampled
  field** — different denominator (109 µM total 2-PE vs the model's ~235 µM).
  [[feedback-rejected-values-must-be-unreachable]]. **A de-novo fraction is not scale-invariant.**
- **The blocker moved, it did not lift** — now "does the Phe flux scale with total 2-PE?" The **"only a T4
  snapshot, it climbs later" defence is REFUTED** (flat). The shipped 2-PE cap is **INERT at the realistic
  dose** — a carbon-refund guard, **not** the 18.9% fix; never conflate.
- The live lever on isoamyl is `f_non_ehrlich_leucine`, not a route. Crépin's 0.815 (shipped) vs Minebois's
  implied ~0.29 = open **D-103** conflict: 29.3% vs 77–86% to protein, **never averaged**; her `f` and de-novo
  share are **one study**. Rayne 2016 CALC ester k is 6–18× off R&O MEASURED.

**Esters (D-123 → D-127)**
- Three ester Processes ship, all SPECULATIVE. `k_H2T` is NEGATIVE and shipped faithfully (wine-only). **No pH
  factor** on hexanoate/EtOAc — R&O Table VII ratios are isoamyl's.
- **Direction facts (still true):** MCFA esters hydrolyse **WITH** the acetates; the *forming* family is
  branched/polyprotic ethyl esters (diethyl succinate, ethyl lactate) — none tracked. Ethyl decanoate trap: fast
  decrease, **no uniform C4–C10 constant**.
- **Sourcing corrections. D-131:** book-sweep note **WRONG** (C10 23 µM = "marked delay", not total; C12 white =
  total); primary **Table 6** is the calibration, yeast MCFA deferred. **D-130:** DHA omitted, gluconolactone
  deferred. **D-129:** WHERE = sourced 142, HOW-SHARP = spec. `k_d2`.

## Accepted deviations — recorded, NOT tuned (do not re-litigate as bugs)
Realised Phe share under-shoots (guard-safe); static share ignores feedback inhibition; de-novo decarb
CO₂ uncharged; ester/alcohol ratio marginally >1; `acidbase.py`'s docstring concession costs nothing. **"Bound SO₂ under-modelled" is NOT one — D-143.**

## Open asks / external

- **Ask Querol** (`aquerol@iata.csic.es`) for raw SI: Phe dose vs total 2-PE. ¹³C Ile rides along.
- **Single-host obligation OPEN** — Minebois rests on one PMC deposit, two live parameters on one figure
  [[feedback-paywalled-is-one-host]]. **PMC + EuropePMC are ONE deposit, not two sources** (D-152 amd 1).
- **D-104's un-inversion** — scoped, UNSOURCED, not started, owner's call. D-116 moved its gate onto **in-situ
  [E] + de-novo-KIC + decarboxylase fluxes**; also prices D-103's leucine conflict.
- Durable findings under `M:\claud_projects\temp\ferment\`: `_findings\`, `d13{5..9}-*\`, `d14{1..9}-*\` — incl.
  `d142-pulls\`+`d143-so2-binding\` (Miao **T2/3/4**), `d149-copper-refit\` (Nguyen **T3.1**), `d151-l16-ph\`
  (Carrasco-Quiroz **T1+2**), `d152-copper-bound\` (**SDs**+bound); `_txt\carrascon-red-kinetics-2018.txt`
  = **Carrascón 2018 reds** (D-150).

## Not started (deferred tail; D-110's narrowing still unconfirmed by owner)

Pham's pH + ethanol terms; growth-linked excretion (D-49 opt B); peptide pool; variety-specific DMSp;
yeast-autolysate spectrum; re-anchor `f_methional`; masking (cosα-blocked); D-55's stale Brett prose;
acetaldehyde in maturation + the 0-vs-2.7 floor; ester `_eq` floors; pH factor for hexanoate/EtOAc; osmotic
inhibition >~200 g/L; `k_d2`; adduct release; closure OTR(T) + bottling burst; no post-Fenton O₂ draw
(D-142). **`add_copper` never writes the `copper` slot** — needs a residual-Cu fraction. D-143/4 ← D-145.

## Standing rule
**Direction is the owner's call, every time** — ask before picking the next milestone or beat (D-66, [[feedback-discuss-disagreements]]).
