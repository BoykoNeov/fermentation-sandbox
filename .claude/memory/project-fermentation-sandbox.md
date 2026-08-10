---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
  modified: 2026-08-10T13:57:34.582Z
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation engine in Python (uv, scipy/numpy/pydantic). Repo: https://github.com/BoykoNeov/fermentation-sandbox (branch `main`).

**Session-boot context: PROHIBITIONS and POINTERS only** — not a changelog. Every bullet is *what it forbids* +
the D-record to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue past
it from this file.** **Caps: 8 lines per BLOCK, 320 chars per `MEMORY.md` index row, 300 lines total as a
BACKSTOP** (`.claude/hooks/check_memory_size.py`, D-169; [[feedback-batch-end-ritual]]). Distil the NEW block;
**never evict an old prohibition to buy a line** — the per-block cap exists so eviction cannot satisfy it.

## Where the records are
- `docs/DECISIONS.md` — canonical archive, **~20.5k lines: never read it linearly.** Generated top block gives a
  subsystem cut, ordered list, and **correction map (⚠)**. **The ⚠ lives only in the index** — check a record's
  index row before trusting it. Append per `CLAUDE.md`'s `Corrects:`/`Flags:`, then `tools/gen_decisions_toc.py`
  (**edit `TOPIC_RULES` if a new record buckets nowhere**). **File is LF.**
- `docs/ARCHITECTURE.md` (seams); `docs/plans/milestone-*.md` ("Active beat: sensory" header is **STALE**, closed
  at D-95/D-98); `CLAUDE.md` (prime directives + archive conventions).

## Status (2026-08-10)
M0/M1/M2 **complete**. **Milestone 3** (sensory/OAV + Tier-3 aging, owner's pick at D-66) in progress, at
**D-173**. Aging build order **built** — `aging.py` carries 24 Processes; sensory 1a/1b closed. **D-139's
leftovers ALL closed** (§2.4 D-148, §2.5 D-149). Suite **1532 passed**. Wine **94 slots** / beer **47**, `quinone`
in both regardless of set; **three** oxidative sets (`direct` default/`cascade`/`direct_burst`). Most remaining
Milestone-3 work is **blocked on external sourcing**, not on building.

## Do NOT re-propose — I did, twice, from stale "Next:" breadcrumbs
[[feedback-verify-latest-state-not-breadcrumbs]]. **A D-record's own "Next:" is a breadcrumb list too** — D-156's
still named the withdrawn "under-bound SO₂ pool" (D-143) as open.
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
- **FOUR surfaces (D-157), TWO distributions — PINNED (D-156, `tests/test_sampling_surfaces.py`): do not re-audit,
  do not "simplify".** Use these counts: compile-seam **246 DISTINCT** varying (**247 pre-D-172**) — **NOT 279**
  (a per-*file* sum double-counting shared names); **structural 61, NOT D-157's 66** (D-159, 5 merely
  scenario-inert). Predicate = **declared `reads`**; **declared-by-UNREGISTERED-class is NEVER drawn**. **A shared
  name can carry DIFFERENT, even disjoint, bands per medium** (19 of 33). `psychophysics.yaml` is **UNIFORM** —
  never "fix" to triangular, and **never apply the triangular mass statistic to it** (anti-conservative);
  `sensory.yaml`'s **36 are NEVER sampled — do not re-audit them**. `SHARED_FILES` restated **on purpose**
  (deriving it = D-108/D-109 vacuity). **A distribution test at `x == mode` is vacuous** — sample **off-mode**.
- **D-157's live contradiction CLOSED (D-158) — the band WON; never re-narrow 0.084 to 0.08.** Resolved
  **INTERNALLY, no fetch**: which number is *sourced*? 84 = Shinohara's 16.4% E-rate; **30–80 occurred ONCE, in the
  test comment asserting it**. Corrects **D-127**. Test **recomputes** all three (`abs=5e-4`; never `rel`/`round(x,3)` — pins *formatting*). Band = E-rate spread at **FIXED acetic 0.35**, a documented narrowing.
- **`reads` has TWO masters — tier propagation AND sampler scope (D-160, fixes D-159's defect).**
  `PH_SYSTEM_READS`/`SO2_BINDING_READS` (`acidbase.py`) → 19 members/5 modules. **Keep DISJOINT** (pinned) and
  **derived**, never re-listed. **`temperature_ramp_rate` stays undeclared BY DESIGN**; no tier moved.
  **RESTATEMENT DONE (D-161) — never re-run it.** Affected class = **ONE** row, a 9-name `only=` **ISOLATION** not a
  band; its "13.99/16.69/20.09 across seeds 0/1/2" is **SORTED** [[feedback-a-majority-is-not-a-direction]]. Through the shipped sampler the fix
  is **undetectable**; predicted **widening NOT observed** (10/24, p=0.54). D-159 pins consumption (`test_drawability_surface.py`).
- **Closure ordering SCOPED (D-162) — do not "fix", narrow or re-measure it.** **Three** claims: P1's
  **three-tier** grouping is its conclusion verbatim (breach **2.4%**); `technical<screwcap` (**40.7%**) and
  `nomacorc<supremecorq` (**0.0%**) are Table-I **nominals**, not that sentence. 94% of the 41.9% chain breach is
  that one pair, **~9/10 of THAT band SCOPE** (like-for-like **4.7%**). Declaration-level ONLY: one
  `scenario.closure` → the **`closure_otr` STATE slot**; forcing all 5 moves **0 of 94** vs **54** (`mu_max`).
  **12** assertions undecided [[feedback-count-and-print-your-skips]] — **the list was never persisted**.
- **Band EDGES measured ARCHIVE-WIDE (D-163), 73 arms — do NOT re-run the sweep.** Of **678** live edges (339
  bands) only **19 guarded** (→**23**, D-164); **652** move both ways suite-green; **14 of 18 files wholly
  unguarded**. **Scale-only** provably **cannot** fire `psychophysics`'s disjointness (needs **translation**); **6
  half-edges immovable** = D-153 Leg 4's. **Flags D-162.** Its "**55 externally-sourced**" is **CORRECTED by D-167**.
  **Prose flag REFUSED, measured** — 44 of 51 hits were `ceiling` alone; **never re-run or "tighten" the regex**.
- **CLASS (d) CLOSED at the seam (D-164) — do NOT recentre the band.** The override mints a `Parameter` carrying
  the base's band ⇒ `_value_in_range` gated `carrying_capacity_gpl`/`autolysis_rate_per_h` by an *epistemic* band,
  naming the **parameter** not the **knob**; now stated (`_override_in_band`, **no new constants**) + pinned
  (`test_scenario_override_bounds.py`, 14). **Multiplicative recentre BUILT, MEASURED, REJECTED** — `k_autolysis`
  log-symmetric ⇒ **3.71×** the request vs shipped **0.67×**: structurally clean, arithmetically wrong. **The
  ensemble/nominal gap is NOT an override defect** — cap non-biting **76.8%** with **NO** override vs 68.6% with
  ⇒ **wide-band** artefact. Range **still** = the band (**Flags**). `test_autolysis`
  sits on the **exact** 1e-2 high edge ⇒ narrowing it breaks a scenario test for non-autolysis reasons.
- **WIDE-BAND artefact SPLIT, Mechanism A CLOSED (D-165) — do NOT re-census bands, do NOT ship a shape field.**
  **TWO** mechanisms: **A** log-width (`(lo+m+hi)/3m`), **B** threshold proximity — D-164's 76.8% is **B, not
  width**; censusing together = the D-162 muddle. **A:** 351→**339 live**, `r≥10` **123 (36.5%)**, `lo≤0` **6**, **no
  VALIDATED band live** (**338/121/83** log-sym MEASURED at D-173; a whole-dir dedupe gives **306**, a DIFFERENT
  denominator); worst `k_d2_ethanol_tolerance_death` **r=300**; **27** say "order of
  magnitude" **in prose no code reads**. **Corrects D-24**: outer percentiles de-sensitise across the **two linear
  shapes offered** (1.30×), **not** log-scale (**3.71×**) — pinned `tests/test_band_shape_sampling.py` (8), **6 arms
  killing 6 of 8**. **Attribution NOT monotone in `r`.** Log-tri **NOT shipped** [[feedback-pair-the-arm-with-its-baseline]].
- **Mechanism B CLOSED — the switch-site census (D-166); B's surface is PARAMETERS, not code sites** (corrects
  D-165's "Next": a `params[` grep misses `ethanol_tolerance`). Classify **bitwise** (fixed `t_eval`), **never by
  epsilon**; **populate the must** first. Only **2 of 66** read params are switch-gated ⇒ **`reads` reach is an
  UPPER BOUND**, **and inertness is a MARGIN** — held open by two byproduct-yield bands and closed by a
  **documented toggle** (Gay-Lussac diversion OFF) [[feedback-a-margin-is-a-claim-about-what-holds-it-open]]. **A
  nominal ON a band edge is inert by construction** [[feedback-nominal-on-a-band-edge-is-not-inertness]].
  **`ethanol_inhibition` is NOT in wine's active set**; **beer is UNCENSUSABLE** (`_ALLOWED_KEYS["beer"]`: no acid
  input, pH 7.0). Pinned `tests/test_switch_site_census.py`.
- **The "buildable 55" is RESCOPED, first slice CLOSED (D-167) — corrects D-163's §7 headline**, whose predicate
  tested the **VALUE's** `provenance.source` and was reported about **EDGES**
  [[feedback-name-the-field-your-predicate-read]]. Of 110 edges only **4** are a parameter-specific published span,
  **73 have NO account**; the **strong** class is **DERIVED** ⇒ **no source-field predicate can find it**.
  **`acidbase`'s 13 are the EMPTIEST slice — TERMINAL** (CRC is a reference book ⇒ [[feedback-paywalled-is-one-host]]
  does **NOT** transfer). Screen wrong twice [[feedback-a-text-screen-has-units-and-self-reference]]. **NO GUARD SHIPS.**
- **Those 5 are CLOSED (D-168) — do NOT re-open them as "bands that contradict their sources".** Four carry a
  **SECOND** account (a **Q10 at 20 °C**) D-167 never read, near-disjoint from the first ⇒ **"move it to its cited
  range" is ILL-POSED** [[feedback-a-note-can-state-its-span-twice]]; thermal's matches neither ⇒
  **AUTHOR-CONSTRUCTED**. Decider = **the ORDERINGS at JOINT band edges** (D-118's shape) ⇒ **`E_a_oak_extraction`
  high 35k → 25k SHIPPED** (justified by the **citation**, ordering repair a **consequence**
  [[feedback-name-guards-for-what-they-forbid]]); wine `E_a_esters` ships **55,100 == `E_a_uptake`** (D-21). Beer's
  **DISSOLVES**, leaving the axis's **tightest joint margin, +2,000 J/mol**. **NEVER assign edges to accounts from
  proximity**, and **do not re-draft the WITHDRAWN `Corrects: D-77`** (nominal-scoped, never states the band).
- **D-168 §1b's Q10 sweep is CLOSED (D-170) — never re-run it expecting a defect count.** It returns **15 again,
  for the OPPOSITE reason**: repairing a note makes a worst-gap screen score it **WORSE** ⇒ **NOT idempotent under
  its own repair**; a naive re-run reads "no progress". **Field-complete (all FOUR fields) reproduced 21/15
  EXACTLY — the predicted delta did NOT happen.** **7 rounding-class = notational, MEASURED**, left unedited on
  purpose. `E_a_acet_reduction`/`E_a_smm_hydrolysis`: **wrong when WRITTEN, not stale** (D-172 inverts this).
- **The recurring shape: 7 invariants BREACH and NOT ONE can be guarded where it breaches (D-171; corrects
  D-170's "5").** Recall holes found the extra two: **a `>` is not an ordering WORD**, the **AMBIGUOUS-PARSE
  quarantine was never revisited**, *"steeper"* **inverts on a retention fraction**. `E_a_decarb`'s shape does NOT
  transfer (margin **exactly 0**, intent UNRECORDED) while all 7 are **NEGATIVE** ⇒ `low >= high` needs an edge
  moved and **an `author estimate` cannot license a narrowing ⇒ ZERO edges moved**. Guards **nominal-scoped ONLY**,
  each naming the breach it does NOT forbid (`y_acetaldehyde_per_tannin` 34.2 %, `k_death_brett` 23.1, `E_a_fusels`
  0.005); `y_methanethiol` 19.4 stays **UNGUARDED** — a guard would pin my own prose.
- **A guard can NAME the ordering and still not fire on it (D-171 §4)** — `test_integrated_wine_aroma_
  temperature_directions` survives `E_a_fusels`/`E_a_uptake` **fully inverted**: its sensitivity sits
  **above** the −3,000 J/mol its own bands permit.
- **O2 partition REPARAMETERISED (D-172) — never re-add `k_ethanol_oxidation`/`k_browning_base`/`k_activation_floor`.**
  ⇒ `k_o2_depletion_total` × `f_ethanol_o2_share`, **two PRODUCTS, never `total-total*f`** (1 ULP; 361 GREEN).
  Breach **55.20 → 12.54 %**, **ZERO edges moved**. **Never 3 co-drawn — max 2.** Sum assert now an **IDENTITY ⇒
  vacuous**. **The total's high 1.0e-3 has NO account** (printed alt **2.0e-3**, D-71) — **OWNER'S CALL.**
- **D-172 amdt — SAME defect one level up, FLAGGED not repaired:** `E_a_ethanol_oxidation`/`E_a_browning` are **TWO
  entries**, identical value AND band ⇒ split T-independent **at `T_ref` ONLY**; D-74's *"exact at every T"* is nominal-only.
- **D-172 §6 SHIPPED (D-173) — the ONE measurement-licensed edge move; NEVER re-lower to 1.0e-4** (it drew
  9.68e-3 > 8.0e-3 air sat). `k_o2_depletion_total` = **[2.4e-4, 1.0e-3]**: **licence is the MEASUREMENT, not the
  [2.4e-4, 1.2e-3] account** ⇒ **one edge of two ON PURPOSE, HIGH still unaccounted/owner's**. **The joint is NOT
  (total×share) — share CANCELS**: 38 params swept, **12 raise `o2` ⇒ margin 1.14×/cascade 1.04×, NOT 1.9× — do
  NOT re-run.** Guard = compile seam + **monotonicity**. **Residual `closure_otr` UNSAMPLED — Flags D-136.**
- **D-89's sotolon caution is FLAGGED, not resolved (D-170 §6) — do NOT treat it as a live bound.**
  `k_maillard_browning`'s **whole band** moves sotolon OAV by **<1e-4** in the scenario its own note names
  (positive control moving; pool ends **0.63 of 0.8 g/L** ⇒ competition not binding). **Not a value
  change**; deciding it needs D-89's calibration scenario, which it **never pins**.
- **Exactly TWO parameters have a real bound**: `k_copper_multiplier` (**band clamped to it, D-154**) and
  `f_de_novo_2_phenylethanol` (**D-118, guarded** — breach point *recomputed*, the template both follow).
- **`f_non_ehrlich_phenylalanine`'s HIGH edge is load-bearing** for D-118's floor (joint margin **3.07e-5**);
  top-mode deliberate (0.531 = hard *measured* protein floor). **D-153's "unguarded" was WRONG; D-155 REFUSED
  the guard** (4-arm matrix caught every arm ⇒ decoration); hardened instead — the breach test reads the band
  **edge**, not the nominal. **Its high edge is D-163's one *asserted-but-never-moved* pin** (value==high ⇒ no
  operator can shift it). **Run mutants BEFORE building a guard**
  [[feedback-mutate-the-premise-before-building-the-guard]] — D-155 REFUSED, D-156 LICENSED, D-157 nothing,
  D-161 REFUSED. **0.963 inside its band is a COINCIDENCE** — D-119 forbids *assigning* it, not drawing near it [[feedback-check-the-schema-not-the-caller]].

**Oxidation (D-132 → D-137, D-149 → D-152)**
- D-132's phenolic boost is **additive, never proportional**, browning-side only. D-133's `burst_antioxidant` is
  **EXCESS** over it and must **read none of** `tannin`/`anthocyanin`/`so2_total` — anti-double-count, and binds.
- `k_copper_multiplier` = **600** L/g, **§2.5 CLOSED (D-149)** — never re-open on "the printed table says copper
  is stronger" (Nguyen T3.1 → 2092, failing the rejected 2000's budget). Held by the **source**, not Ferreira's
  ceiling; re-open needs **real wine, ≥3 Cu levels**. Guard 4.
- **BOUNDED, not null (D-152):** copper-orthogonal L16 gives **k ≤ 918 L/g**, excluding 2092, 2000 and band-high
  1500 under all six arms (Guard 7). **Never rebuild on Table 2's SDs; keep condition 12; NEVER cite it as
  evidence FOR 600** (only 1 of 6 arms admits 600) — one-directional, against higher k only.
- **BAND FIXED (D-154) — was the live defect; do not re-propose it.** High **1500 → 662.8** L/g: the bound at
  `copper_typical`'s **MAXIMUM**, not the shipped-centring 918, because `copper_typical` is **itself sampled**
  and `k_bound` *decreases* in it. **Never take 918**: over 200k joint draws it still violates **5.01%** of pairs
  (1500 violates **37.6%**, so D-152's 29% **understates** it). Verify sampling claims **on draws, not edges**.
  **Never adopt D-152's printed "663"** — exact is 662.802522, so 663 ships red. Value/low
  edge/`copper_typical` **untouched**. Guard **recomputes** the bound + asserts monotonicity, never reads the
  note [[feedback-rejected-values-must-be-unreachable]]. **No Fe(III) state** — not D-134's "iron in surplus": QSS rests on a **~18× separation**.
- **A pH term on activation/`k_browning_eff` is REFUSED (D-150)** — never re-open on "Nguyen's table shows strong
  pH dependence". Strongest leg is **inseparability from copper**, corroborated by Carrasco-Quiroz's
  copper-**orthogonal** L16 (D-151) — never dismiss that as Nguyen's dosing. His is an *initial* statistic on a
  *steady* node ⇒ **the pH term's home is the burst**. Guard 5's bound is an **OBSERVED spread, NOT a ceiling**;
  re-open needs a **within-wine** steady-OCR series, ≥3 levels, measuring Cu.
- **`initial_ph` anchors t=0 only** — **no way to SET an aging pH**. First-order in `[o2]`; **no MM/`Km`**. **Do
  NOT re-attempt the Ferreira/Carrascón PLS extraction** — all three blocked (`_findings/D-134-*.md`).
- **The O2 gate is Fe(II)+O2**; SO2 is right in *size*, **INVERTED in role** (enabler, not competitor) — which is
  why **D-72's wrong mechanism yields a right-looking 1:2**; never read that ratio as confirmation. **Never
  re-fix the acetaldehyde/phenolic inversion inside the parallel frame — it provably cancels**; approve the
  cascade on MECHANISM only, **never tune to phenol ⇒ more MeCHO**. Danilewicz 2011 is **PARAPHRASE** — pull the
  abstract before it backs a `source:`; in-wine **mixing-limited**, atmospheric pH-independence does **NOT**
  transfer to wine.

**Burst — WIRED and NON-DEFAULT (D-147). Not dead code, not unbuilt, not refused.**
- `oxidative="direct_burst"` = direct **+** burst, a **superset**, not a mechanism. **No `cascade_burst`**
  (D-138 stays UNDETERMINED; its "transient modifier" is FLAGGED).
- **Never make it default** — moves all 31 D-140 pins and re-opens the Danilewicz direct arm, **moot: do not
  re-derive.** **The split is UNPINNED and stays so** — C1 pins only the **product**, C2 fails **structurally**.
  **Never re-solve it** (D-146's second-unmeasured-number trap). Unlock: **Ferreira 2015 per-cycle O₂ curves**.
- **Self-exhaustion is Ferreira's PROTOCOL, not the sink** (~1400× the cork flux) ⇒ a permanent **~37% tax that
  GROWS to a plateau**; two guards forbid calling it transient — **never relax to a "pool depletes" check**.
  **Seed follows the consumer**: **0.0** outside `direct_burst`; dosing `burst_antioxidant_gpl` with no consumer
  **raises** (D-45's "absent ≠ 0" **inverts**). Ferreira's **2.7× does not discriminate**, and is paraphrase.

**Cascade — BUILT and NON-DEFAULT (D-141). Do not re-build; do not flip it silently.**
- `core/kinetics/oxidative_cascade.py`, 8 Processes, `quinone` the one new slot, OFF the ledger (`h2o2` QSS is
  **`1/k` NOT `ln2/k`**). **REPLACEMENT** via `_OXIDATIVE_SETS` (`get_medium` /
  `compile_scenario(..., oxidative=)`). **`"direct"` is default and stays default.**
- **Do NOT close the fate gap by tuning** — fates move up to **25×** purely from re-homing the rate law. Both
  settling sources are **paywalled (nine hosts)**; abstract K values are **never a `source:` field**
  [[feedback-conceded-caveats-are-not-coverage]].
- **Never re-derive the activation floor as a fit** (name retired at D-172, above) — the O2 budget
  agreeing is **supply limitation (D-136), NOT a rate-law check**. Activation **reads the reductant
  pools** — never lump. The 31 D-140 guards stand — **never re-derive their pins**; edit only the two
  seams. Beer's O2 is a **floor `>=5 mg/L`**, never 5.71. Pin rtol **1e-4**
  [[feedback-pin-tolerance-vs-solver-tolerance]]; `quinone == 0.0` under direct is **exact**.
- **Benchmark EXISTS and is ACTIVE** (`tests/benchmarks/test_validation_danilewicz_so2_o2.py`, 15 tests, no
  xfail/skip — **not open work**). **Never pin 1.7** — one dataset's mode, above the other's whole range; assert
  the limits 1/2 + bands. The falsifier is the **traverse** (D-141's "structurally cannot produce" was wrong). Quinone branching **NOT settled**.
- **Operating point is load-bearing — enforce the >10 free-SO₂ floor, never state it in prose.** It is
  `SIM_CURVATURE_FLOOR_MGL`, **NOT Miao's criterion — never re-encode his**; keep the value + excluded-wines
  table. **Never report a single-dose verdict on this ratio** (direct unconditional, cascade **straddles**). The
  **"4–7× quinone shortfall" + its `xfail` are WITHDRAWN**.
- **"The sim under-binds SO₂" is WITHDRAWN (D-143) — never re-assert**; never compare an addition-method secant
  with an oxidation-path slope (concretely `oxidation_path_slope` vs `MIAO_BUFFERING_BAND`); **never re-run a
  pool sweep**. Acetaldehyde 0.0000 mM at D-142's operating point is a **scenario artefact**, not a schema gap. **A green suite proves nothing about a quoted decimal — every assertion is a band.**
- **The four D-142 artefacts are FIXED (D-144)** — Miao's **Table 4 intercepts ARE free SO₂ at O₂ exhaustion**:
  **read the intercept, never recompute from an assumed dose**. D-144's test **NAMES NO CAUSE — keep it that
  way**. Locus guard is **pure algebra** (a state route re-entangles `ph_of_state`); masses from `core.chemistry`
  (**M_SO2 = 64.058**); `_SO2_BINDERS` reads all **four**.

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
  artefacts**; an undershot pool **freezes negative**. Aging Process **names** are enumerated by `test_fusel_*`.

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

**Esters (D-123 → D-127)**
- Three ester Processes ship, all SPECULATIVE. `k_H2T` is NEGATIVE and shipped faithfully (wine-only). **No pH
  factor** on hexanoate/EtOAc — R&O Table VII ratios are isoamyl's. **Direction facts:** MCFA esters hydrolyse
  **WITH** the acetates; the *forming* family is branched/polyprotic ethyl esters — none tracked. Ethyl decanoate
  trap: **no uniform C4–C10 constant**. **Sourcing corrections. D-131:** book-sweep note **WRONG**, primary **Table 6** is the calibration. **D-130:** DHA omitted, gluconolactone deferred. **D-129:** WHERE/HOW-SHARP.

## Accepted deviations — recorded, NOT tuned (do not re-litigate as bugs)
Realised Phe share under-shoots (guard-safe); static share ignores feedback inhibition; de-novo decarb CO₂
uncharged; ester/alcohol ratio marginally >1; `acidbase.py` docstring concession. **"Bound SO₂ under-modelled" is NOT one — D-143.** (pKa sampling gap: was a live defect, fixed D-160, restated D-161.)

## Open asks / external
- **Ask Querol** (`aquerol@iata.csic.es`) for raw SI: Phe dose vs total 2-PE. ¹³C Ile rides along.
- **Single-host obligation OPEN** — Minebois rests on one PMC deposit, two live parameters on one figure
  [[feedback-paywalled-is-one-host]]. **PMC + EuropePMC are ONE deposit, not two sources** (D-152 amd 1).
- **D-104's un-inversion** — scoped, UNSOURCED, not started, owner's call. D-116 moved its gate onto **in-situ [E]
  + de-novo-KIC + decarboxylase fluxes**; also prices D-103's leucine conflict.
- Durable findings under `M:\claud_projects\temp\ferment\`: `_findings\`, `d13{5..9}-*\`, `d14{1..9}-*\`,
  `d15{7,8,9}-*\`, `d16{0..8}-*\`, `d170-q10-generalise\`/`d171-ordering-guards\` — incl. `d142-pulls\`+`d143-so2-binding\` (Miao **T2/3/4**), `d149-copper-refit\`
  (Nguyen **T3.1**), `d151-l16-ph\` (Carrasco-Quiroz **T1+2**), `d163-band-edges\`/`d165-wide-band\`/`d166-switch-census\`/`d167-edge-provenance\` (reusable harnesses); `_txt\carrascon-red-kinetics-2018.txt` = **Carrascón 2018 reds**.

## Not started (deferred tail; D-110's narrowing still unconfirmed by owner)
Pham's pH + ethanol terms; growth-linked excretion (D-49 opt B); peptide pool; variety-specific DMSp;
yeast-autolysate spectrum; re-anchor `f_methional`; masking (cosα-blocked); D-55's stale Brett prose; acetaldehyde
in maturation + the 0-vs-2.7 floor; ester `_eq` floors; pH factor for hexanoate/EtOAc; osmotic inhibition >~200 g/L; `k_d2`; adduct release; closure OTR(T) + bottling burst; no post-Fenton O₂ draw (D-142); **`add_copper` never writes the `copper` slot** (needs a residual-Cu fraction); D-143/4 ← D-145.

## Standing rule
- **NEVER raise the total cap again without first showing the per-block cap did not bite (D-169).** It was a
  **TARGET, not a limit**: the file sat at **exactly 250 for 13 consecutive commits / 12 days**, then took **47 of
  the next 50 lines the day it moved to 300** — the fill rate followed the cap, so all four raises were futile.
  Total is now a **BACKSTOP**; **shape** binds. **When it binds, the licensed move is ⚠-COLLAPSE RETIREMENT** —
  a record superseded per the correction map drops to a pointer — **never a raise, never eviction**. A
  **digit-density check was REJECTED on measurement** (confounded by block size; it penalises hardest the
  guardrails that are nothing but corrected values). **`CLAUDE.md` is a boot surface and is UNMEASURED**, and
  `MEMORY.md` is capped **per row, not in total** — row COUNT is the open channel; expect displacement there.

**Direction is the owner's call, every time** — ask before picking the next milestone/beat, offering only UNBLOCKED options (D-66, [[feedback-discuss-disagreements]]).
