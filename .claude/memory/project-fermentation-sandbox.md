---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
  modified: 2026-07-28T19:31:52.305Z
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation engine in Python (uv, scipy/numpy/pydantic). Repo: https://github.com/BoykoNeov/fermentation-sandbox (branch `main`).

**Session-boot context: PROHIBITIONS and POINTERS only** — not a changelog. Every bullet is *what it forbids* +
the D-record to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue past
it from this file.** **Cap 250** (`.claude/hooks/check_memory_size.py`; [[feedback-batch-end-ritual]]).

## Where the records are

- `docs/DECISIONS.md` — canonical archive, **~18.5k lines: never read it linearly.** Generated top block gives a
  subsystem cut, ordered list, and **correction map (⚠)**. **The ⚠ lives only in the index** — check a record's
  index row before trusting it. Append per `CLAUDE.md`'s `Corrects:`/`Flags:`, then `tools/gen_decisions_toc.py`
  (**edit `TOPIC_RULES` if a new record buckets nowhere**). **File is LF.**
- `docs/ARCHITECTURE.md` (seams); `docs/plans/milestone-*.md` ("Active beat: sensory" header is **STALE**, closed
  at D-95/D-98); `CLAUDE.md` (prime directives + archive conventions).

## Status (2026-07-28)

M0/M1/M2 **complete**. **Milestone 3** (sensory/OAV + Tier-3 aging, owner's pick at D-66) in progress, at
**D-160**. Aging build order **built** — `aging.py` carries 24 Processes; sensory 1a/1b closed. **D-139's
leftovers ALL closed** (§2.4 D-148, §2.5 D-149). Suite **1460 passed**. Wine schema **94 slots** / beer **47**
— `quinone` in both regardless of set. **Three** oxidative sets (`direct` default / `cascade` / `direct_burst`).
Most remaining Milestone-3 work is **blocked on external sourcing**, not on building.

## Do NOT re-propose — I did, twice, from stale "Next:" breadcrumbs
[[feedback-verify-latest-state-not-breadcrumbs]]. **A D-record's own "Next:" is itself a breadcrumb list** —
D-156's still named the withdrawn "under-bound SO₂ pool" (D-143) as open.
- **All lumps are speciated**: esters→3 (D-96), fusels→5 (D-99), amino_acids→8 (D-100), mercaptans a
  methanethiol false-lump (D-110); `lumped` stays **dormant**. **Beat 1b (descriptor projection) COMPLETE**
  (D-95 + D-98); only masking remains, `cosα`-blocked.
- **Shipped and spent, not unbuilt:** D-128, D-129 (`EthanolToleranceDeath`), D-130 … D-136. **Isoamyl de-novo
  entry — REFUSED at D-120, measured not built**: a rate knob on a supply-limited quantity
  [[feedback-measure-which-side-before-building]]. **Closed:** leucine shortfall (D-112); shared-BAT parsimony
  (D-116); Rollero (D-115); ester-aging (D-121). **Beer 3-sugar kinetics are NOT in the 5 beer books.**

## Live prohibitions, by axis

**Sampled bands (D-153 → D-157) — BOTH archive-wide sweeps are DONE. Do not re-run either.** The **sampler**
surface (which bands exist) = D-153/D-156; the **assertion** surface (constraints checked at a point) = D-157.
- **THE RECURRING SHAPE, 4 instances (D-118, D-154, D-155, D-157): a constraint verified at a POINT
  where the sampler reads a BAND.** Whenever a guard or bound uses a nominal, check whether that
  quantity is itself sampled — and take the **joint** worst case over every band involved.
- **FOUR surfaces (D-157), TWO distributions.** compile-seam **247 DISTINCT** varying, **triangular** — **NOT
  279**: that is a per-*file* sum double-counting `beer_generic` names `wine_generic` also defines (149+98+32);
  wine 247 / beer 181 per medium are right. **Structural = 61, NOT D-157's 66 (D-159)** — 66 was a union over 13
  scenarios; 5 (`E_a_pof`, `k_pof_decarb{,_ferulic}`, `biomass_carrying_capacity`, `k_iso_alpha_loss`) are merely
  **scenario-inert**. Predicate = **declared `reads`**, not a census. **FOUR** modes: Process-read DRAWN;
  **compile-consumed** (`_closure_otr`, `_iso_alpha_at_pitch`, `*_initial`), **verb-consumed**
  (`_verb_add_copper/_oak`), **declared-by-an-UNREGISTERED-class** (`ethanol_inhibition_exponent`) **NEVER**.
  A 5th — runtime-read via a helper, declared nowhere — existed and is **closed at D-160**.
  **19 of the 33 names in BOTH medium files carry DIFFERENT bands, some disjoint** (`E_a_esters` wine 40–70k vs
  beer 120–265k). `psychophysics.yaml` **24, UNIFORM** (`compression.py::_axis_draws` — deliberate, never "fix"
  to triangular); `sensory.yaml` **36, NEVER sampled** (`load_thresholds`) — **do not re-audit those 36**.
  **Never apply the triangular mass statistic to psychophysics** — 16 points off, anti-conservative.
  **PINNED (D-156, `tests/test_sampling_surfaces.py`) — do not re-audit, do not "simplify".** Its
  `SHARED_FILES` is restated **on purpose** (deriving it = D-108/D-109 vacuity — the same trap that made
  D-159's census version vacuous). **TWO triangular defaults** — `ensemble.py:108` and `:444`, separate
  paths/tests. **A distribution test at `x == mode` is vacuous** — sample **off-mode**. **Removal** is not
  *separately* pinned (all 14 drops already RED by consumption); **addition** is what D-156 closed.
- **D-157's live contradiction CLOSED (D-158) — the band WON; never re-narrow 0.084 to 0.08.** Resolved
  **INTERNALLY, no fetch**: ask which number is *sourced*, not which is right — 84 = Shinohara's 16.4% E-rate;
  **30–80 occurred ONCE in the repo, in the test comment asserting it**. Corrects **D-127**. Test **recomputes**
  all three points (`abs=5e-4` = 3-dp half-width; never `rel`/`round(x,3)` — pins *formatting*). Band = E-rate
  spread at **FIXED acetic 0.35** ⇒ documented narrowing, **NOT** the point-vs-band shape.
- **`reads` has TWO masters — tier propagation AND sampler scope (D-160, fixes D-159's defect).**
  `PH_SYSTEM_READS`/`SO2_BINDING_READS` (`acidbase.py`) → 19 members/5 modules. **Keep DISJOINT** (pinned:
  `EsterHydrolysis`/`EthylAcetateEsterification` solve pH, never partition SO₂) and **derived**, never re-listed.
  **`temperature_ramp_rate` stays undeclared BY DESIGN** (`temperature.py:32-37`) — "declare everything the
  sweep records" is WRONG. No tier moved. **Every pH-dependent aging band predates this and is NOT COMPARABLE**
  (sorted draw order moved, not just width) — re-measure; never reason about the delta. Next beat.
- **Closure ordering: ~42% inverted on joint draws, but declaration-level ONLY — do not "fix" it.** Two `otr_*`
  never coexist in a run (one `scenario.closure` → the **`closure_otr` STATE slot**, `compile.py:2254`; the
  sampler draws parameters, not state). **12** assertions still undecided
  [[feedback-count-and-print-your-skips]]; D-159 pins consumption (`test_drawability_surface.py`).
- **The prose flag is REFUSED, measured not asserted** — 44 of 51 hits matched `ceiling` alone (oak's
  *extraction asymptote*). **Never re-run or "tighten" the regex**: the archive uses `RECALIBRATED`/`ceiling`/
  `upper bound` for HEALTHY provenance too. Only a **measured per-parameter bound** finds it.
- **Exactly TWO parameters have a real bound**: `k_copper_multiplier` (**band clamped to it, D-154**) and
  `f_de_novo_2_phenylethanol` (**D-118, guarded** — breach point *recomputed*, the template both follow).
- **`f_non_ehrlich_phenylalanine`'s HIGH edge is load-bearing** for D-118's floor — joint margin
  **3.07e-5**. **D-153's "unguarded" was WRONG; D-155 REFUSED the guard** (4-arm matrix caught every
  arm ⇒ decoration); hardened instead — the breach test reads the band **edge**, not the nominal.
  **Run mutants BEFORE building a guard** [[feedback-mutate-the-premise-before-building-the-guard]] —
  now decided three ways: D-155 REFUSED, D-156 LICENSED, D-157 licensed nothing. Its top-mode is
  deliberate (0.531 = hard *measured* protein floor). **0.963 inside its band is a COINCIDENCE, not a
  second defect** — D-119 forbids *assigning* it (different denominator), not drawing near it
  [[feedback-check-the-schema-not-the-caller]].

**Oxidation (D-132 → D-137, D-149 → D-152)**
- D-132's phenolic boost is **additive, never proportional**, browning-side only.
- D-133's `burst_antioxidant` is **EXCESS** over D-132's steady rate and must **read none of**
  `tannin`/`anthocyanin`/`so2_total` — the anti-double-count, and it binds.
- `k_copper_multiplier` = **600** L/g, **§2.5 CLOSED (D-149)** — never re-open on "the printed table says copper
  is stronger" (Nguyen T3.1 → 2092, failing the rejected 2000's budget). Held by the **source**, not Ferreira's
  ceiling; re-open needs **real wine, ≥3 Cu levels**. Guard 4.
- **BOUNDED, not null (D-152):** copper-orthogonal L16 gives **k ≤ 918 L/g**, excluding 2092, 2000 and band-high
  1500 under all six arms (Guard 7). **Never rebuild on Table 2's SDs; keep condition 12; NEVER cite it as
  evidence FOR 600** (only 1 of 6 arms admits 600) — one-directional, against higher k only.
- **BAND FIXED (D-154) — was the live defect; do not re-propose it.** High **1500 → 662.8** L/g: the bound at
  `copper_typical`'s **MAXIMUM**, not the shipped-centring 918, because `copper_typical` is **itself sampled**
  and `k_bound` *decreases* in it — joint draws bind at its max. **Never take 918**: over 200k joint draws it
  still violates **5.01%** of pairs (1500 violates **37.6%**, so D-152's 29% — at a *fixed* centring —
  **understates** it). Verify sampling claims **on draws, not edges**. **Never adopt D-152's printed "663"** —
  exact is 662.802522, so 663 ships red. Value/low edge/`copper_typical` **untouched**; Ferreira's 2.2× budget
  is **no longer** the high-side constraint. Guard **recomputes** the bound + asserts monotonicity, never reads
  the note [[feedback-rejected-values-must-be-unreachable]].
- **No Fe(III) state** — not D-134's "iron in surplus": QSS rests on a **~18× separation**.
- **A pH term on activation/`k_browning_eff` is REFUSED (D-150)** — never re-open on "Nguyen's table shows strong
  pH dependence". Strongest leg is **inseparability from copper**, corroborated by Carrasco-Quiroz's
  copper-**orthogonal** L16 (D-151) — never dismiss that as Nguyen's dosing. His is an *initial* statistic on a
  *steady* node ⇒ **the pH term's home is the burst**. Guard 5's bound is an **OBSERVED spread, NOT a ceiling**.
  Re-open needs a **within-wine** steady-OCR series, ≥3 levels, measuring Cu.
- **`initial_ph` anchors t=0 only** — **no way to SET an aging pH**. First-order in `[o2]`; **no MM/`Km`**. **Do
  NOT re-attempt the Ferreira/Carrascón PLS extraction** — all three blocked (`_findings/D-134-*.md`).
- **The O2 gate is Fe(II)+O2**; SO2 is right in *size*, **INVERTED in role** (enabler, not competitor) — which is
  why **D-72's wrong mechanism yields a right-looking 1:2**; never read that ratio as confirmation.
- **Never re-fix the acetaldehyde/phenolic inversion inside the parallel frame — it provably cancels.** Approve
  the cascade on MECHANISM only; **never tune to phenol ⇒ more MeCHO.**
- Danilewicz 2011 is **PARAPHRASE** — pull the abstract before it backs a `source:`; its in-wine branch is
  **mixing-limited**, and its atmospheric pH-independence does **NOT** transfer to wine.

**Burst — WIRED and NON-DEFAULT (D-147). Not dead code, not unbuilt, not refused.**
- `oxidative="direct_burst"` = direct **+** burst, a **superset**, not a mechanism. **No `cascade_burst`**
  (D-138 stays UNDETERMINED; its "transient modifier" is FLAGGED).
- **Never make it default** — moves all 31 D-140 pins and re-opens the Danilewicz direct arm, **moot: do not
  re-derive.** **The split is UNPINNED and stays so** — C1 pins only the **product**, C2 fails **structurally**.
  **Never re-solve it** (D-146's second-unmeasured-number trap). Unlock: **Ferreira 2015 per-cycle O₂ curves**.
- **Self-exhaustion is Ferreira's PROTOCOL, not the sink** (~1400× the cork flux) ⇒ a permanent **~37% tax that
  GROWS to a plateau**; two guards forbid calling it transient — **never relax to a "pool depletes" check**.
  **Seed follows the consumer**: **0.0** outside `direct_burst`; dosing `burst_antioxidant_gpl` with no consumer
  **raises** (D-45's "absent ≠ 0" **inverts**). Ferreira's **2.7× does not discriminate** and is D-133
  **paraphrase** — pull it before it backs a `source:`.

**Cascade — BUILT and NON-DEFAULT (D-141). Do not re-build; do not flip it silently.**
- `core/kinetics/oxidative_cascade.py`, 8 Processes, `quinone` the one new slot, OFF the ledger (`h2o2` QSS is
  **`1/k` NOT `ln2/k`**). **REPLACEMENT** via `_OXIDATIVE_SETS` (`get_medium` /
  `compile_scenario(..., oxidative=)`). **`"direct"` is default and stays default.**
- **Do NOT close the fate gap by tuning** — fates move up to **25×** purely from re-homing the rate law. Both
  settling sources are **paywalled (nine hosts)**; abstract K values are **never a `source:` field**
  [[feedback-conceded-caveats-are-not-coverage]].
- **Never re-derive `k_activation_floor` as a fit** — it **is** `k_ethanol_oxidation + k_browning_base` and does
  **not** reproduce direct's TOTAL; the O2 budget agreeing is **supply limitation (D-136), NOT a rate-law
  check**. Activation **reads the reductant pools** — never lump.
- The 31 D-140 guards stand — **never re-derive their pins**; edit only the two seams. Beer's O2 is a **floor
  `>=5 mg/L`**, never 5.71. Pin rtol **1e-4** [[feedback-pin-tolerance-vs-solver-tolerance]]; `quinone == 0.0`
  under direct is **exact**.
- **Benchmark EXISTS and is ACTIVE** (`tests/benchmarks/test_validation_danilewicz_so2_o2.py`, 15 tests, no
  xfail/skip — **not open work**). **Never pin 1.7** — one dataset's mode, above the other's whole range; assert
  the limits 1/2 + bands. The falsifier is the **traverse** (D-141's "structurally cannot produce" was wrong).
- **Operating point is load-bearing — enforce the >10 free-SO₂ floor, never state it in prose.** It is
  `SIM_CURVATURE_FLOOR_MGL`, **NOT Miao's criterion — never re-encode his**; keep the value + excluded-wines
  table, never restore the name. **Never report a single-dose verdict on this ratio** (direct unconditional,
  cascade **straddles**). The **"4–7× quinone shortfall" + its `xfail` are WITHDRAWN**.
- **"The sim under-binds SO₂" is WITHDRAWN (D-143) — never re-assert**; never compare an addition-method secant
  with an oxidation-path slope (concretely `oxidation_path_slope` vs `MIAO_BUFFERING_BAND`); **never re-run a
  pool sweep**. Acetaldehyde 0.0000 mM at D-142's operating point is a **scenario artefact**, not a schema gap.
- **The four D-142 artefacts are FIXED (D-144)** — Miao's **Table 4 intercepts ARE free SO₂ at O₂ exhaustion**:
  **read the intercept, never recompute from an assumed dose**. D-144's test **NAMES NO CAUSE — keep it that
  way**. Locus guard is **pure algebra** (a state route re-entangles `ph_of_state`); masses from `core.chemistry`
  (**M_SO2 = 64.058**); `_SO2_BINDERS` reads all **four**.
- **A green suite proves nothing about a quoted decimal — every assertion is a band.** Quinone branching **NOT
  settled**.

**Dosing schedule — REFUSED (D-145). Never re-propose must-dosing; benchmark unchanged.**
- **D-143's "five for five" is three for five.** **Import statistics from the shipped benchmark module** —
  re-derived helpers carried both errors. **The reservoir CANNOT move the Table 3 secant** — saturated both
  ends ⇒ moves the **INTERCEPT**, not the slope; `free0*` is **floor-CONTINGENT, never state it bare**.
  **Depletion gap RETIRED, not reframed** — D-143 omitted `M_SO2/M_O2`; only the secant is invariant, in band.
- **Four traps that each produce a plausible green:** `ProcessSet.disable()` is silently undone by `begin_aging`;
  `param_values` is a **property returning a fresh dict**; mol/L ×1000 is mmol/L not mg/L; **`trajectory.y` is
  `(n_states, n_times)` — `y[-1]` is the LAST SLOT's series**, not the final state (D-147).
- **§2.4 CLOSED (D-148) — never re-open as "the Brett/quench over-draw": no quench draw exists.** Live pair was
  Brett/POF; **summing is NOT the mechanism**. **Never build the shared depletion gate** — `depletion_gate`'s 8
  sites buy per-draw first-orderness `_decarboxylation_branch` already has. Out-of-band numbers are **BDF step
  artefacts**. Live but unreachable: an undershot pool **freezes negative**, `conservation.py` weights both
  unclamped. Cross-domain trap: aging Process **names** are enumerated by `tests/test_fusel_*.py`.

**Milestone-3 tail (D-146) — two REFUSED, one BLOCKED; do not re-open as unbuilt**
- **Methionine sink + `methionol` is BLOCKED, not deferred.** Don't build it, don't pick a value:
  `f_non_ehrlich_methionine` is **ill-posed** (leaves its own `[0,1)` domain across one panel). The D-118 gate
  **passes** ⇒ the block is on the parameter; Crépin/Rollero **structurally incapable**, not silent. Unlock:
  **¹³C-methionine tracer**.
- **Sotolon's OAV is NOT enantiomer-split** — `threshold_sotolon_wine` is the **racemate's**, pool racemic; D-107's "soft by ~100×" is **withdrawn**.
- **Do not "finish" the `oav`→`magnitude` rename.** Scoped to `DescriptorReading.oav` alone —
  `OAVReading.oav`/`oav_series`/`oav_tier` are genuine OAVs and a test fails if swept. **`TOPIC_RULES`
  (`tools/gen_decisions_toc.py`) is where a new record's bucket is fixed — edit the rules, not the heading.**
  Has `methion`/`sotolon`/`assertion`; a `magnitude` rule **misfiles D-53** and `surface` **misfiles D-128**.

**Sulfides (D-135) + closures (D-136)**
- Do NOT rebuild bottle reduction as thioacetate/disulfide precursors — D-101's mechanism guess was WRONG.
  **Release only** (white MeSH ~4× under-predicted: a limitation, not a bug); **no copper coupling** (a future
  one must be asymmetric); **temperature-FLAT** — never ship my ~158 kJ/mol.
- **Technical cork ships BELOW screwcap** (nominals only — D-157); µL→µg is the authors' own **1.43** (STP); seven
  "do not fix" traps in D-136. `ClosureOxygenIngress.reads` is **`()` by design**, already pinned.

**Fusels / 2-PE (D-117 → D-120)**
- `f_non_ehrlich_phenylalanine` = 0.975, `f_de_novo_2_phenylethanol` = 0.9827. **Never put 0.963 in a sampled
  field** — different denominator [[feedback-rejected-values-must-be-unreachable]]. **A de-novo fraction is not
  scale-invariant.**
- **The blocker moved, it did not lift** — now "does the Phe flux scale with total 2-PE?"; the **"only a T4
  snapshot" defence is REFUTED** (flat). The shipped 2-PE cap is **INERT at the realistic dose** — a carbon-refund
  guard, **not** the 18.9% fix; never conflate. The live lever on isoamyl is `f_non_ehrlich_leucine`, not a route.
  Crépin's 0.815 vs Minebois's implied ~0.29 = open **D-103** conflict — **never averaged**. Rayne 2016 CALC
  ester k is 6–18× off R&O MEASURED.

**Esters (D-123 → D-127)**
- Three ester Processes ship, all SPECULATIVE. `k_H2T` is NEGATIVE and shipped faithfully (wine-only). **No pH
  factor** on hexanoate/EtOAc — R&O Table VII ratios are isoamyl's. **Direction facts:** MCFA esters hydrolyse
  **WITH** the acetates; the *forming* family is branched/polyprotic ethyl esters — none tracked. Ethyl decanoate
  trap: **no uniform C4–C10 constant**. **Sourcing corrections. D-131:** book-sweep note **WRONG**, primary
  **Table 6** is the calibration. **D-130:** DHA omitted, gluconolactone deferred. **D-129:** WHERE/HOW-SHARP.

## Accepted deviations — recorded, NOT tuned (do not re-litigate as bugs)
Realised Phe share under-shoots (guard-safe); static share ignores feedback inhibition; de-novo decarb CO₂
uncharged; ester/alcohol ratio marginally >1; `acidbase.py` docstring concession. **"Bound SO₂
under-modelled" is NOT one — D-143.** (The pKa sampling gap WAS a live defect; fixed at D-160.)

## Open asks / external

- **Ask Querol** (`aquerol@iata.csic.es`) for raw SI: Phe dose vs total 2-PE. ¹³C Ile rides along.
- **Single-host obligation OPEN** — Minebois rests on one PMC deposit, two live parameters on one figure
  [[feedback-paywalled-is-one-host]]. **PMC + EuropePMC are ONE deposit, not two sources** (D-152 amd 1).
- **D-104's un-inversion** — scoped, UNSOURCED, not started, owner's call. D-116 moved its gate onto **in-situ
  [E] + de-novo-KIC + decarboxylase fluxes**; also prices D-103's leucine conflict.
- Durable findings under `M:\claud_projects\temp\ferment\`: `_findings\`, `d13{5..9}-*\`, `d14{1..9}-*\`,
  `d15{7,8,9}-*\`, `d160-undeclared-reads\` — incl. `d142-pulls\`+`d143-so2-binding\` (Miao **T2/3/4**),
  `d149-copper-refit\` (Nguyen **T3.1**), `d151-l16-ph\` (Carrasco-Quiroz **T1+2**), `d152-copper-bound\`
  (**SDs**+bound); `_txt\carrascon-red-kinetics-2018.txt` = **Carrascón 2018 reds** (D-150).

## Not started (deferred tail; D-110's narrowing still unconfirmed by owner)

Pham's pH + ethanol terms; growth-linked excretion (D-49 opt B); peptide pool; variety-specific DMSp;
yeast-autolysate spectrum; re-anchor `f_methional`; masking (cosα-blocked); D-55's stale Brett prose;
acetaldehyde in maturation + the 0-vs-2.7 floor; ester `_eq` floors; pH factor for hexanoate/EtOAc; osmotic
inhibition >~200 g/L; `k_d2`; adduct release; closure OTR(T) + bottling burst; no post-Fenton O₂ draw (D-142).
**`add_copper` never writes the `copper` slot** — needs a residual-Cu fraction (binding constants are
verb-consumed ⇒ never sampled). D-143/4 ← D-145.

## Standing rule
**Direction is the owner's call, every time** — ask before picking the next milestone/beat, offering only
UNBLOCKED options (D-66, [[feedback-discuss-disagreements]]).
